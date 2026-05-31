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
#define DT 0.05
#define TOTAL_STEPS 10000
#define SAVE_EACH 300

// Physics
#define G 9.81
#define MANNING_N 0.035

// Urrá Conditions
#define H_UPSTREAM 130.5
#define H_DOWNSTREAM 68.0
#define DAM_X_POS 0
#define DAM_HEIGHT 73.0
#define SLOPE 0.00015

// Breaking of Dam
#define BREACH_WIDTH 200.0
#define BREACH_START_TIME 5.0
#define BREACH_FORMATION_TIME 60.0

// Control
#define USE_REAL_TOPOGRAPHY 1

int main() {
    printf("\n========================================\n");
    printf("SIMULACIÓN ROTURA DE PRESA - URRÁ\n");
    printf("========================================\n\n");
    
    printf("Configuración del dominio:\n");
    printf("  Distancia total: %.1f km\n", (NX * DX) / 1000.0);
    printf("  Ancho del valle: %.1f km\n", (NY * DY) / 1000.0);
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
        charge_real_topography(zb, NX, NY, "data/urra_specific_topography.bin");
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

    printf("\nIniciando Simulación...\n");
    double start_time = omp_get_wtime();

    for (int step = 0; step < TOTAL_STEPS; step++) {
        double actual_time = step * DT;
        apply_breach(h, u, v, zb, NX, NY, DX, DY, 
                    actual_time, DAM_X_POS, BREACH_WIDTH, 
                    BREACH_START_TIME, BREACH_FORMATION_TIME);
        solver_pass(h, u, v, h_new, u_new, v_new, zb,
                    NX, NY, DX, DY, DT, MANNING_N);

        #pragma omp parallel for
        for (int i = 0; i < total_cells; i++) {
            h[i] = h_new[i];
            u[i] = u_new[i];
            v[i] = v_new[i];
        }

        if (step % SAVE_EACH == 0) {
            save_snapshot(h, step, NX, NY, "snapshots");
            if (step % (SAVE_EACH * 10) == 0) {
                double progress = (double)step / TOTAL_STEPS * 100;
                double sim_time = step * DT;
                printf("Progreso: %.1f%% | Tiempo: %.0f s (%.1f min)\n", 
                    progress, sim_time, sim_time / 60.0);
            }
        }

        if (step % 1000 == 0) {
            int error = 0;
            for (int i = 0; i < total_cells && !error; i++) {
                if (isnan(h[i]) || h[i] > 1000.0) {
                    printf("ERROR: Inestabilidad en paso %d, celda %d\n", step, i);
                    error = 1;
                }
            }
            if (error) break;
        }
    }

    double end_time = omp_get_wtime();
    printf("\nSimulación completada en %.2f segundos\n", end_time - start_time);
    
    double final_volume = calculate_volume(h, total_cells, DX, DY);
    printf("Volumen final: %.3f millones m³\n", final_volume / 1e6);
    
    save_snapshot(h, TOTAL_STEPS, NX, NY, "snapshots");

    free(h); free(u); free(v); free(h_new); 
    free(u_new); free(v_new); free(zb);

    printf("\nSimulación finalizada.\n");
    return 0;
}