#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "compute.h"
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

int row, col, _iter;
double direct_n_c = sqrt(2) / (sqrt(2) + 1) / 4;
double diagnol_n_c = 1 / (sqrt(2) + 1) / 4;

static inline double calculate_t(int r, int c, const double* tinit, const double* conductivity) {
    int index = r*row + c;
    if (r==0 || r==row-1) return tinit[index]; // boundry points
    if ((r==1 || r==row-2) && _iter > 0) return tinit[index]; // halo grid points
    
    // calculate neighbors location
    int u_index = (r - 1)*row + c;
    int d_index = (r + 1)*row + c;
    int l_index = r*row + (c - 1 + col) % col;
    int r_index = r*row + (c + 1 + col) % col;
    int lu_index = (r - 1)*row + (c - 1 + col) % col;
    int ru_index = (r - 1)*row + (c + 1 + col) % col;
    int ld_index = (r + 1)*row + (c - 1 + col) % col;
    int rd_index = (r + 1)*row + (c + 1 + col) % col;

    // calculate current grid point temperature
    double remainting_c = 1 - conductivity[index];
    double res = conductivity[index] * tinit[index] + 
        remainting_c * direct_n_c * (tinit[u_index] + tinit[d_index] + tinit[l_index] + tinit[r_index]) +
        remainting_c * diagnol_n_c * (tinit[lu_index] + tinit[ru_index] + tinit[ld_index] + tinit[rd_index]);
    
    return res;
}

static void draw_picture(const struct parameters* p, const double* world, size_t key){
    int N = p->N, M = p->M;
    begin_picture(key, M, N, p->io_tmin, p->io_tmax);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=M;j++){
            draw_point(j-1,i-1,world[j*(N+2)+i]);
        }
    }
    end_picture();

}

void do_compute(const struct parameters* p, struct results *r)
{
    double tmin = p->io_tmax + 1, tmax = p->io_tmin - 1, tavg, time;
    double maxdiff = abs(tmin - tmax);
    size_t iter = p->maxiter;

    const double* tinit = p->tinit;
    const double* conductivity = p->conductivity;
    struct timespec start, end;

    // calculate temperature
    double* tmp = 0;
    double sum = 0.0;
    bool terminate = false;
    row = p->M, col = p->N;
    double* _tinit = (double*)malloc(row * col * sizeof(double));
    
    clock_gettime(CLOCK_MONOTONIC, &start); // time begin
    // interate maxiter times
    for(_iter = 0; _iter < iter; _iter++) {
        if(terminate) break; // check convergence needed...
        
        sum = 0.0;
        for(size_t r = 0; r < row; r++) {
            for(size_t c = 0; c < col; c++) {
                int index = r*row + c;
                _tinit[index] = calculate_t(r, c, tinit, conductivity);
                tmin = min(tmin, _tinit[index]);
                tmax = max(tmax, _tinit[index]);
                maxdiff = max(maxdiff, fabsf(_tinit[index] - tinit[index]));
                sum += _tinit[index];

                // if (maxdiff > p->threshold ) {
                //     terminate = false;
                // }
            }
        }
        tavg = sum/(row * col);
        r->niter = _iter+1;
        r->tmin = tmin;
        r->tmax = tmax;
        r->maxdiff = maxdiff;
        r->tavg = tavg;
        
        tmp = tinit;
        tinit = _tinit;
        _tinit = tmp;
    }

    // calculate characteristic values
    // int total = row * col;
    // double sum = 0.0;
    // for(size_t i = 0; i < total; i++) {
    //     tmin = min(tmin, tinit[i]);
    //     tmax = max(tmax, tinit[i]);
    //     sum += tinit[i];
    // }
    // maxdiff = abs(tmin - tmax);
    // tavg = sum/total;
    clock_gettime(CLOCK_MONOTONIC, &end); // time end
    time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
    // end of calculation, save results
    // r->niter = _iter;
    // r->tmin = tmin;
    // r->tmax = tmax;
    // r->maxdiff = maxdiff;
    // r->tavg = tavg;
    r->time = time;

    draw_picture(p, tinit, 1);
    free(_tinit);
}
