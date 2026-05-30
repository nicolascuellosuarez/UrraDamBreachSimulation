#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init.h"

void initialize_simulation_parameters(double *h, double *u, double *v, double *zb, 
                                    int NX, int NY, double dx, double dy, double h_upstream,
                                    double h_downstream, double dam_x_pos, double dam_height,
                                    double slope, int real_topography) 
{
    int i, j;
    int idx;
    double x, y;

    for (i = 0; i < NX * NY; i++) {
        // Setting Velocities to zero.
        *(u + i) = 0.0;
        *(v + i) = 0.0;
    }

    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            idx = i + j * NX;
            x = i * dx;
            y = j * dy;

            if (real_topography && zb != NULL) {
                // Nothing to do, because "zb" already has values
            } else {
                zb[idx] = 50.0 - slope * x + 2.0 * sin(0.001 * y);
            }

            if (i < dam_x_pos) {
                h[idx] = h_upstream - zb[idx];
                if (h[idx] < 0) h[idx] = 0;
            } else {
                h[idx] = h_downstream - zb[idx];
                if (h[idx] < 0) h[idx] = 0;
            }

            if (h[idx] < 0.01) h[idx] = 0.01; // Minimum Value of water level
        }
    }

    printf("La inicialización de la Simulación fue completada:\n");
    printf("Profundidad máxima aguas arriba: %.2d m\n", h_upstream);
    printf("Profundidad aguas abajo: %.2d m\n", h_downstream);
    printf("Posición de la presa: celda %d (%.0f m)\n", dam_x_pos, dam_x_pos * dx);
    printf("Pendiente del Terreno: %.6f\n", slope);
}