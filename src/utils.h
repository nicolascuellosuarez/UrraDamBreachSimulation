#ifndef UTILS_H
#define UTILS_H

double calculate_volume(double *h, int n, double dx, double dy);
void save_snapshot(double *h, int step, int nx, int ny, const char *folder);
void charge_real_topography(double *zb, int nx, int ny, const char *filename);
void apply_breach(double *h, double *u, double *v, double *zb, int nx, int ny, 
                double dx, double dy, double time, int dam_x, double breach_width, 
                double breach_start, double breach_time);
void solver_pass(double *h, double *u, double *v, double *h_new, double *u_new, double *v_new,
                double *zb, int nx, int ny, double dx, double dy, double dt, double n_manning);

#endif