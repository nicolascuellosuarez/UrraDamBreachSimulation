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
            
            double h_lim = fmax(h[idx], 0.05);
            double h_left = fmax(h[(i-1) + j*nx], 0.05);
            double h_right = fmax(h[(i+1) + j*nx], 0.05);
            double h_bottom = fmax(h[i + (j-1)*nx], 0.05);
            double h_top = fmax(h[i + (j+1)*nx], 0.05);
            
            double hu_left = h_left * u[(i-1)+j*nx];
            double hu_right = h_lim * u[idx];
            double dhu_dx = (hu_right - hu_left) / dx;
            
            double hv_bottom = h_bottom * v[i+(j-1)*nx];
            double hv_top = h_lim * v[idx];
            double dhv_dy = (hv_top - hv_bottom) / dy;
            
            double new_h = h_lim - dt * (dhu_dx + dhv_dy);
            
            if (new_h < 0.05) new_h = 0.05;
            if (new_h > 200.0) new_h = 200.0;
            h_new[idx] = new_h;
            
            double eta_center = h_lim + zb[idx];
            double eta_left = h_left + zb[(i-1)+j*nx];
            double eta_right = h_right + zb[(i+1)+j*nx];
            double eta_bottom = h_bottom + zb[i+(j-1)*nx];
            double eta_top = h_top + zb[i+(j+1)*nx];
            
            double deta_dx = (eta_right - eta_left) / (2.0 * dx);
            double deta_dy = (eta_top - eta_bottom) / (2.0 * dy);
            
            double max_grad = 0.15;
            if (fabs(deta_dx) > max_grad) deta_dx = (deta_dx > 0) ? max_grad : -max_grad;
            if (fabs(deta_dy) > max_grad) deta_dy = (deta_dy > 0) ? max_grad : -max_grad;
            
            double u_new_val = u[idx] - dt * G * deta_dx;
            double v_new_val = v[idx] - dt * G * deta_dy;
            
            double max_vel = 30.0;
            if (u_new_val > max_vel) u_new_val = max_vel;
            if (u_new_val < -max_vel) u_new_val = -max_vel;
            if (v_new_val > max_vel) v_new_val = max_vel;
            if (v_new_val < -max_vel) v_new_val = -max_vel;
            
            double speed = sqrt(u_new_val * u_new_val + v_new_val * v_new_val);
            if (speed > 0.1 && new_h > 0.1) {
                double Sf = (n_manning * n_manning * speed * speed) / pow(new_h, 4.0/3.0);
                double fric_factor = dt * G * Sf / (speed + 0.1);
                if (fric_factor > 0.5) fric_factor = 0.5;
                u_new_val = u_new_val * (1.0 - fric_factor);
                v_new_val = v_new_val * (1.0 - fric_factor);
            }
            
            u_new[idx] = u_new_val;
            v_new[idx] = v_new_val;
        }
    }
    
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        int idx = 0 + j * nx;
        h_new[idx] = h_new[1 + j * nx];
        u_new[idx] = 0.0;
        v_new[idx] = v_new[1 + j * nx];
    }
    
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        int idx = (nx-1) + j * nx;
        h_new[idx] = h_new[(nx-2) + j * nx];
        u_new[idx] = 0.0;
        v_new[idx] = v_new[(nx-2) + j * nx];
    }
    
    #pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        int idx = i + 0 * nx;
        h_new[idx] = h_new[i + 1 * nx];
        u_new[idx] = u_new[i + 1 * nx];
        v_new[idx] = 0.0;
    }
    
    #pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        int idx = i + (ny-1) * nx;
        h_new[idx] = h_new[i + (ny-2) * nx];
        u_new[idx] = u_new[i + (ny-2) * nx];
        v_new[idx] = 0.0;
    }
}