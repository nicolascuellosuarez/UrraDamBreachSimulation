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


        }
    }
}