#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "utils.h"
#define G 9.81

double calculate_volume(double *h, int n, double dx, double dy) {
    double volume = 0.0;
    #pragma omp parallel for reduction(+:volume)
    for (int i = 0; i < n; i++) {
        volume += h[i];
    }
    return volume * dx * dy;
}

void save_snapshot(double *h, int step, int nx, int ny, const char *folder) {
    char filename[256];
    snprintf(filename, sizeof(filename), "%s/h_%06d.bin", folder, step);
    FILE *f = fopen(filename, "wb");
    if (f) {
        fwrite(h, sizeof(double), nx * ny, f);
        fclose(f);
    }
}

void charge_real_topography(double *zb, int nx, int ny, const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) {
        printf("ERROR: No se pudo abrir %s\n", filename);
        return;
    }
    size_t elements = fread(zb, sizeof(double), nx * ny, f);
    fclose(f);
    if (elements == nx * ny) {
        printf("Topografía real cargada: %d x %d\n", nx, ny);
    } else {
        printf("ERROR: Solo se leyeron %zu de %d elementos\n", elements, nx * ny);
    }
}

void apply_breach(double *h, double *u, double *v, double *zb, int nx, int ny, double dx, 
                  double dy, double time, int dam_x, double breach_width, double breach_start, 
                  double breach_time) {
    if (time < breach_start) return;
    
    double factor = 1.0;
    if (time < breach_start + breach_time) {
        factor = (time - breach_start) / breach_time;
    }
    
    int cells_breach = (int)(breach_width / dx);
    if (cells_breach < 1) cells_breach = 1;
    
    int start_breach = (ny / 2) - (cells_breach / 2);
    int end_breach = start_breach + cells_breach;
    
    int j, idx;
    #pragma omp parallel for private(j, idx)
    for (j = start_breach; j < end_breach; j++) {
        if (j >= ny) continue;
        idx = dam_x + j * nx;
        if (zb[idx] > 50.0) {
            zb[idx] = 50.0 + (zb[idx] - 50.0) * (1.0 - factor);
        }
    }
}

void solver_pass(double *h, double *u, double *v, double *h_new, double *u_new, double *v_new,
                 double *zb, int nx, int ny, double dx, double dy, double dt, double n_manning) {
    int total = nx * ny;
    
    #pragma omp parallel for
    for (int i = 0; i < total; i++) {
        h_new[i] = h[i];
        u_new[i] = u[i];
        v_new[i] = v[i];
    }
    
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < ny - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            int idx = i + j * nx;
            
            double h_left = h[(i-1) + j*nx];
            double h_right = h[(i+1) + j*nx];
            double u_avg = (u[(i+1)+j*nx] + u[idx]) / 2.0;
            
            double h_bottom = h[i + (j-1)*nx];
            double h_top = h[i + (j+1)*nx];
            double v_avg = (v[i+(j+1)*nx] + v[idx]) / 2.0;
            
            double dhx = (h_right - h_left) / (2.0 * dx);
            double dhy = (h_top - h_bottom) / (2.0 * dy);
            
            double new_h = h[idx] - dt * (u_avg * dhx + v_avg * dhy);
            
            h_new[idx] = (new_h > 0.01) ? new_h : 0.01;
            
            u_new[idx] = u[idx] - dt * G * (h_new[idx] - h_new[(i-1)+j*nx]) / dx;
            v_new[idx] = v[idx] - dt * G * (h_new[idx] - h_new[i+(j-1)*nx]) / dy;
        }
    }
}