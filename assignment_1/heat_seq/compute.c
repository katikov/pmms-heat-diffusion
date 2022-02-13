#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

int row, col, _iter;

inline double calculate_t(int r, int c, const double *tinit, const double *conductivity)
{
    int index = r * row + c;
    if (r == 0 || r == row - 1)
        return tinit[index]; // boundry points
    if ((r == 1 || r == row - 2) && _iter > 0)
        return tinit[index]; // halo grid points

    // calculate conductivity coefficients
    double remainting_c = 1 - conductivity[index];
    double direct_n_c = remainting_c * sqrt(2) / (sqrt(2) + 1) / 4;
    double diagnol_n_c = remainting_c / (sqrt(2) + 1) / 4;

    // calculate neighbors location
    int u_index = (r - 1) * row + c;
    int d_index = (r + 1) * row + c;
    int l_index = r * row + c - 1;
    int r_index = r * row + c + 1;
    int lu_index = (r - 1) * row + c - 1;
    int ru_index = (r - 1) * row + c + 1;
    int ld_index = (r + 1) * row + c - 1;
    int rd_index = (r + 1) * row + c + 1;

    // calculate current grid point temperature
    double res = conductivity[index] * tinit[index] +
                 direct_n_c * (tinit[u_index] + tinit[d_index] + tinit[l_index] + tinit[r_index]) +
                 diagnol_n_c * (tinit[lu_index] + tinit[ru_index] + tinit[ld_index] + tinit[rd_index]);

    return res;
}

void do_compute(const struct parameters *p, struct results *r)
{
    double tmin = p->io_tmax + 1, tmax = p->io_tmin - 1, tavg, time;
    double maxdiff = fabs(tmin - tmax);
    size_t iter = p->maxiter;

    const double *tinit = p->tinit;
    const double *conductivity = p->conductivity;

    // calculate temperature
    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    double *tmp;
    row = p->M, col = p->N;
    double *_tinit = (double *)malloc(row * col * sizeof(double));
    int flag = 0;
    for (_iter = 0; _iter < iter && flag == 1; _iter++)
    {
        // interate maxiter times
        // check convergence needed...
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < col; c++)
            {
                int index = r * row + c;
                double res = _tinit[index] = calculate_t(r, c, tinit, conductivity);
                if (fabs(res - tinit[index]) > p->threshold)
                    flag = 1;
            }
        }

        tmp = tinit;
        tinit = _tinit;
        _tinit = tmp;
    }

    // calculate characteristic values
    int total = row * col;
    double sum = 0.0;
    for (size_t i = 0; i < total; i++)
    {
        tmin = min(tmin, tinit[i]);
        tmax = max(tmax, tinit[i]);
        sum += tinit[i];
    }
    tavg = sum / total;
    struct timespec now_time;
    clock_gettime(CLOCK_MONOTONIC, &now_time);
    

    // end of calculation, save results
    r->niter = _iter;
    r->tmin = tmin;
    r->tmax = tmax;
    r->maxdiff = maxdiff;
    r->tavg = tavg;
    r->time = now_time.tv_sec - start_time.tv_sec 
        + 1e-9*(now_time.tv_nsec - start_time.tv_nsec);;

    free(_tinit);
}