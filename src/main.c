#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

// Domain
#define NX 600
#define NY 520
#define DX 50.0
#define DY 50.0
#define TOTAL_TIME 600.0  

// Physics
#define G 9.81
#define MANNING_N 0.035

// Urrá Conditions
#define H_UPSTREAM 130.5
#define H_DOWNSTREAM 68.0
#define DAM_X_POS 314
#define DAM_HEIGHT 73.0
#define SLOPE 0.00015

// Breaking of Dam
#define BREACH_WIDTH 200.0
#define BREACH_START_TIME 5.0
#define BREACH_FORMATION_TIME 60.0

// Control
#define USE_REAL_TOPOGRAPHY 1
#define CFL_NUMBER 0.08  

double calculate_dt_adaptive(double *h, double *u, double *v, int nx, int ny, double dx, double cfl);

int main() {
    printf("\n========================================\n");
    printf("SIMULACIÓN ROTURA DE PRESA - URRÁ\n");
    printf("========================================\n\n");
    
    printf("Configuración del dominio:\n");
    printf("  Distancia total: %.1f km\n", (NX * DX) / 1000.0);
    printf("  Ancho del valle: %.1f km\n", (NY * DY) / 1000.0);
    printf("  Posición de la presa: celda %d (%.0f m)\n", DAM_X_POS, DAM_X_POS * DX);
    printf("  Resolución: %.0f x %.0f metros por celda\n", DX, DY);
    printf("  Total de celdas: %d (%.2f millones)\n\n", NX * NY, (NX * NY) / 1e6);

    int total_cells = NX * NY;
    double *h = (double*)malloc(total_cells * sizeof(double));
    double *u = (double*)malloc(total_cells * sizeof(double));
    double *v = (double*)malloc(total_cells * sizeof(double));
    double *h_new = (double*)malloc(total_cells * sizeof(double));
    double *u_new = (double*)malloc(total_cells * sizeof(double));
    double *v_new = (double*)malloc(total_cells * sizeof(double));
    double *zb = (double*)malloc(total_cells * sizeof(double));

    if (!h || !u || !v || !h_new || !u_new || !v_new || !zb) {
        printf("Error reservando memoria.\n");
        return -1;
    }

    if (USE_REAL_TOPOGRAPHY) {
        printf("Cargando topografía real...\n");
        charge_real_topography(zb, NX, NY, "/home/nicolas/UrraDamBreachSimulation/data/urra_specific_topography.bin");
    } else {
        printf("Usando topografía sintética\n");
    }

    printf("Inicializando condiciones de la presa.\n");
    initialize_simulation_parameters(h, u, v, zb, NX, NY, DX, DY, 
                                    H_UPSTREAM, H_DOWNSTREAM, 
                                    DAM_X_POS, DAM_HEIGHT, 
                                    SLOPE, USE_REAL_TOPOGRAPHY);
    
    double initial_volume = calculate_volume(h, total_cells, DX, DY);
    printf("\nVolumen inicial: %.3f millones de m³\n", initial_volume / 1e6);

    printf("\nIniciando Simulación con CFL adaptativo...\n");
    double start_time = omp_get_wtime();

    double current_time = 0.0;
    int step = 0;
    
    while (current_time < TOTAL_TIME) {
        double dt = calculate_dt_adaptive(h, u, v, NX, NY, DX, CFL_NUMBER);
        if (dt < 0.001) dt = 0.001;
        if (current_time + dt > TOTAL_TIME) dt = TOTAL_TIME - current_time;
        
        apply_breach(h, u, v, zb, NX, NY, DX, DY, 
                    current_time, DAM_X_POS, BREACH_WIDTH, 
                    BREACH_START_TIME, BREACH_FORMATION_TIME);
        
        solver_pass(h, u, v, h_new, u_new, v_new, zb,
                    NX, NY, DX, DY, dt, MANNING_N);
        
        #pragma omp parallel for
        for (int i = 0; i < total_cells; i++) {
            h[i] = h_new[i];
            u[i] = u_new[i];
            v[i] = v_new[i];
        }
        
        current_time += dt;
        step++;
        
        static double last_save_time = 0;
        if (current_time - last_save_time >= 30.0 || step == 1) {
            save_snapshot(h, step, NX, NY, "/home/nicolas/UrraDamBreachSimulation/snapshots");
            last_save_time = current_time;
            
            double max_h = 0.0, min_h = 1e9, max_u = 0.0, max_v = 0.0;
            for (int i = 0; i < total_cells; i++) {
                if (h[i] > max_h) max_h = h[i];
                if (h[i] < min_h) min_h = h[i];
                if (fabs(u[i]) > max_u) max_u = fabs(u[i]);
                if (fabs(v[i]) > max_v) max_v = fabs(v[i]);
            }
            printf("T=%.1fs | paso=%d | dt=%.4f | h=[%.4f,%.2f] | u_max=%.2f | v_max=%.2f\n", 
                   current_time, step, dt, min_h, max_h, max_u, max_v);
        }
        
        if (step % 500 == 0) {
            int error = 0;
            for (int i = 0; i < total_cells && !error; i++) {
                if (isnan(h[i]) || isnan(u[i]) || isnan(v[i])) {
                    printf("ERROR: NaN en paso %d, celda %d\n", step, i);
                    error = 1;
                }
            }
            if (error) break;
        }
    }

    double end_time = omp_get_wtime();
    printf("\n========================================\n");
    printf("RESUMEN FINAL\n");
    printf("========================================\n");
    printf("Simulación completada en %.2f segundos\n", end_time - start_time);
    printf("Pasos totales: %d\n", step);
    printf("Tiempo simulado: %.1f segundos (%.1f minutos)\n", current_time, current_time/60.0);
    
    double final_volume = calculate_volume(h, total_cells, DX, DY);
    printf("\nVolumen inicial: %.3f millones m³\n", initial_volume / 1e6);
    printf("Volumen final:   %.3f millones m³\n", final_volume / 1e6);
    printf("Diferencia: %.3f millones m³ (%.1f%%)\n", 
           (initial_volume - final_volume) / 1e6,
           (initial_volume - final_volume) / initial_volume * 100);
    
    save_snapshot(h, step, NX, NY, "snapshots");

    free(h); free(u); free(v); free(h_new); 
    free(u_new); free(v_new); free(zb);

    printf("\nSimulación finalizada.\n");
    return 0;
}

double calculate_dt_adaptive(double *h, double *u, double *v, int nx, int ny, double dx, double cfl) {
    double max_speed = 0.0;
    
    #pragma omp parallel for reduction(max:max_speed)
    for (int i = 0; i < nx * ny; i++) {
        double c = sqrt(G * fmax(h[i], 0.05));
        double vel = sqrt(u[i]*u[i] + v[i]*v[i]);
        double wave_speed = c + vel;
        if (wave_speed > max_speed) max_speed = wave_speed;
    }
    
    double dt = cfl * dx / (max_speed + 0.001);
    return fmin(dt, 0.5);
}