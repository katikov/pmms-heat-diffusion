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

static inline double calculate_diff(int r, int c, const double* tinit, double* new_tinit, const double* conductivity) {
    int old_index = r*row + c;
    int new_index = (r+1) * (col+2) + (c+1);
    // if (r==0 || r==row-1) return tinit[index]; // boundry points
    // if ((r==1 || r==row-2) && _iter > 0) return tinit[index]; // halo grid points
    
    // calculate neighbors location
    // int u_index = (r - 1)*row + c;
    // int d_index = (r + 1)*row + c;
    // int l_index = r*row + (c - 1 + col) % col;
    // int r_index = r*row + (c + 1 + col) % col;
    // int lu_index = (r - 1)*row + (c - 1 + col) % col;
    // int ru_index = (r - 1)*row + (c + 1 + col) % col;
    // int ld_index = (r + 1)*row + (c - 1 + col) % col;
    // int rd_index = (r + 1)*row + (c + 1 + col) % col;
    int u_index = (r-1) * (col+2) + (c+1);
    int d_index = (r+2) * (col+2) + (c+1);
    int l_index = (r+1) * (col+2) + c;
    int r_index = (r+1) * (col+2) + (c+2);
    int lu_index = (r-1) * (col+2) + c;
    int ru_index = (r-1) * (col+2) + (c+2);
    int ld_index = (r+2) * (col+2) + c;
    int rd_index = (r+2) * (col+2) + (c+2);

    // calculate current grid point temperature
    double remainting_c = 1 - conductivity[old_index];
    double old_temperature = tinit[new_index];
    double new_temperature = conductivity[old_index] * tinit[new_index] + 
        remainting_c * direct_n_c * (tinit[u_index] + tinit[d_index] + tinit[l_index] + tinit[r_index]) +
        remainting_c * diagnol_n_c * (tinit[lu_index] + tinit[ru_index] + tinit[ld_index] + tinit[rd_index]);

    new_tinit[new_index] = new_temperature;
    double* tmp = new_tinit;
    new_tinit = tinit;
    tinit = tmp;

    return fabs(old_temperature - new_temperature);
}

static void draw_picture(const struct parameters* p, const double* world, size_t key) {
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
    double tmin = p->io_tmax + 1, tmax = p->io_tmin - 1, tavg, time, maxdiff;
    size_t iter = p->maxiter;

    const double* tinit = p->tinit;
    const double* conductivity = p->conductivity;
    struct timespec start, end;

    // calculate temperature
    row = p->M, col = p->N;
    double* tmp = 0;
    bool terminate = false;
    double* _tinit = (double*)malloc((row+2) * (col+2) * sizeof(double));
    double* new_tinit = (double*)malloc((row+2) * (col+2) * sizeof(double));

    /* Initiate new tinit array
        +-----+-----+-----+-----+-----+-----+-----+
        |(0,0)|     |     |     |     |     |     |
        +-----+-----+-----+-----+-----+-----+-----+
        |     |(1,1)|     |     |     |     |     |
        +-----+-----+-----+-----+-----+-----+-----+
        |     |     |     |     |     |     |     |
        +-----+-----+-----+-----+-----+-----+-----+
        |Right|     |     |     |     |     | Left|
        +-----+-----+-----+-----+-----+-----+-----+
        |     |     |     |     |     |(M,N)|     |
        +-----+-----+-----+-----+-----+-----+-----+
        |     |     |     |     |     |     |(M+1,|
        |     |     |     |     |     |     | N+1)|
        +-----+-----+-----+-----+-----+-----+-----+
    */
    // initial rows: 0 and row + 1
    _tinit[0] = tinit[col-1];
    _tinit[col+1] = tinit[0];
    _tinit[(row+2) * (col+2) - col] = tinit[row*col - 1];
    _tinit[(row+2) * (col+2) - 1] = tinit[row*col - col];
    for(int i = 0; i < col; i++) {
        _tinit[i+1] = tinit[i];
        _tinit[col*(row+1) + i+1] = tinit[col*(row-1) + i];
    }
    // initial cols: 0 and col + 1
    for(int i = 0; i < row; i++) {
        _tinit[(i+1)*(col+2)] = tinit[i*col + col - 1];
        _tinit[(i+1)*(col+2) + col + 1] = tinit[i*col];
    }
    // inital tinit
    for(int i = 0; i < row; i++)
        for(int j = 0; j < col; j++)
            _tinit[(i+1)*col + (j+1)] = tinit[i*col+j];

    for(int i = 0; i < (row+2)*(col+2); i++) {
        new_tinit[i] = _tinit[i];
    }

    clock_gettime(CLOCK_MONOTONIC, &start); /* time begin */
    
    for(_iter = 0; _iter < iter; _iter++) {
        if(terminate) break; // check convergence
        
        for(size_t r = 0; r < row; r++) {
            for(size_t c = 0; c < col; c++) {
                double diff = calculate_diff(r, c, _tinit, new_tinit, conductivity);
                maxdiff = max(maxdiff, diff);
                
                if (maxdiff > p->threshold ) {
                    terminate = false;
                }
            }
        }
    }

    // calculate characteristic values
    int total = row * col;
    double sum = 0.0;
    int idx;
    for(int i=0; i<row; i++) {
        for(int j=0; j<col; j++) {
            idx = (i+1) * (col+2) + (j+1);
            tmin = min(tmin, tinit[idx]);
            tmax = max(tmax, tinit[idx]);
            sum += tinit[idx];
        }
    }
    tavg = sum/total;

    clock_gettime(CLOCK_MONOTONIC, &end); /* time end */
    time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
    
    // end of calculation, save results
    r->niter = _iter;
    r->tmin = tmin;
    r->tmax = tmax;
    r->maxdiff = maxdiff;
    r->tavg = tavg;
    r->time = time;

    // draw_picture(p, tinit, 1);
    free(_tinit);
}
