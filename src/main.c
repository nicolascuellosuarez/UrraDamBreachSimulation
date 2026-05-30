#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "init.h"

// Domain
#define NX 600
#define NY 520
#define DX 50.0
#define DY 50.0
#define DT 0.05
#define TOTAL_STEPS 144000
#define SAVE_EACH 300

// Physics
#define G 9.81
#define MANNING_N 0.035

// Urrá Conditions.
#define H_UPSTREAM 130.5
#define H_DOWNSTREAM 68.0
#define DAM_X_POS 0
#define DAM_HEIGHT 73.0
#define SLOPE 0.00015

// Breaking of Dam.
#define BREACH_WIDTH 200.0
#define BREACH_START_TIME 5.0
#define BREACH_FORMATION_TIME 60.0

// Control of the use of real topography data.
#define USE_REAL_TOPOGRAPHY 0

double calculate_volume(double *h, int n);
void solver_pass(double *h, double *u, double *v, double *h_new, double *u_new, double *v_new,
                double *zb, int nx, int ny, double dx, double dy, double dt, double n_manning);
void apply_breach(double *h, double *u, double *v, double *zb, int nx, int ny, double dx, 
                double dy, double time, int dam_x, double breach_width, double brach_start, 
                double breach_time);
void save_snapshot(double *h, int pass, int nx, int ny);
void charge_real_topography(double *zb, int nx, int ny, const char *filename);


int main() {
    printf("\n========================================\n");
    printf("SIMULACIÓN ROTURA DE PRESA - URRÁ\n");
    printf("========================================\n\n");
    
    printf("Configuración del dominio:\n");
    printf("  Distancia total: %.1f km (Urrá → Montería)\n", (NX * DX) / 1000.0);
    printf("  Ancho del valle: %.1f km\n", (NY * DY) / 1000.0);
    printf("  Resolución: %.0f x %.0f metros por celda\n", DX, DY);
    printf("  Total de celdas de simulación: %d (%.2f millones)\n\n", NX * NY, (NX * NY) / 1e6);

    // Reserving memory for arrays.

    int total_of_shells = NX * NY;
    double *h = (double*)malloc(total_of_shells * sizeof(double));
    double *u = (double*)malloc(total_of_shells * sizeof(double));
    double *v = (double*)malloc(total_of_shells * sizeof(double));
    double *h_new = (double*)malloc(total_of_shells * sizeof(double));
    double *u_new = (double*)malloc(total_of_shells * sizeof(double));
    double *v_new = (double*)malloc(total_of_shells * sizeof(double));
    double *zb = (double*)maloc(total_of_shells * sizeof(double));

    if (!h || !u || !v || !h_new || !u_new || !v_new || !zb) {
        printf("Error reservando memoria para las matrices.\n");
        return -1;
    }

    if (USE_REAL_TOPOGRAPHY) {
        printf("Cargando topografía real desde archivo .tif... \n");
        charge_real_topography(zb, NX, NY, "data/topography.tif");
    } else {
        printf("Usando topografía sintética (pendiente %.6f)\n", SLOPE);
    }

    printf("Inicializando condiciones de la presa. \n");
    initialize_simulation_parameters(h, u, v, zb, NX, NY, DX, DY, 
                                    H_UPSTREAM, H_DOWNSTREAM, 
                                    DAM_X_POS, DAM_HEIGHT, 
                                    SLOPE, USE_REAL_TOPOGRAPHY);
    
    double initial_volume = calculate_volume(h, total_of_shells);
    printf("\nVolumen inicial del embalsado: %.3f millones de m³\n", initial_volume / 1e6);

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
        for (int i = 0; i < total_of_shells; i++) {
            *(h + i) = h_new[i];
            *(u + i) = u_new[i];
            *(v + i) = v_new[i];
        }

        if (step % SAVE_EACH == 0) {
            save_snapshot(h, step, NX, NY);
            if (step % (SAVE_EACH * 10) == 0) {
                double progress = (double)step / TOTAL_STEPS * 100;
                double simulation_time = step * DT;
                printf("Progreso: %.1f%% | Tiempo Simulado: %.0f segundos (%.1f minutos)\n", 
                    progress, simulation_time, simulation_time / 60.0);
            }
        }

        if (step % 1000 == 0) {
            int is_there_an_error = 0;
            for (int i = 0; i < total_of_shells && !is_there_an_error; i++) {
                if (isnan(h[i]) || h[i] > 1000.0) {
                    printf("ERROR: Inestabilidad numérica en el paso %d, celda %d\n", step, i);
                    is_there_an_error = 1;
                }
            }
            if (is_there_an_error) {
                printf("Finalizando la simulación por inestabilidad numérica.\n");
                break;
            }
        }
    }

    double final_time = omp_get_wtime();
    printf("\n Simulación completada en %.2f segundos \n", final_time - start_time);
    double final_volume = calculate_volume(h, total_of_shells);
    printf("El volúmen final: %.3f millones m³\n", final_volume / 1e6);
    printf("DIferencia (pérdida por borde): %.3f millones de m³ (%.1f%%)\n", (initial_volume - final_volume) / 1e6, (initial_volume - final_volume) / initial_volume * 100);
    save_snapshot(h, TOTAL_STEPS, NX, NY);

    free(h); free(u); free(v); free(h_new); 
    free(u_new); free(v_new), free(zb);

    printf("\nSimulación finalizada. Snapshots guardados en carpeta 'snapshots/'\n");
    printf("Total de snapshots: ~%d\n", TOTAL_STEPS / SAVE_EACH);
    printf("Ejecuta el script de Python para visualizar los resultados.\n");
    return 0;
}
