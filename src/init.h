#ifndef INIT_H
#define INIT_H

void initialize_simulation_parameters(double *h, double *u, double *v, double *zb, 
                                    int NX, int NY, double dx, double dy, double h_upstream,
                                    double h_downstream, double dam_x_pos, double dam_height,
                                    double slope, int real_topography);

#endif