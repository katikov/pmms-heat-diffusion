#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "compute.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

int row, col;
double direct_n_c;
double diagnol_n_c;
static inline double calculate_diff(int r, int c, const double *tinit, double *new_tinit, const double *conductivity)
{
    int old_index = r * row + c;
    int new_index = r * (col + 2) + (c + 1);

    int u_index = new_index - (col + 2);
    int d_index = new_index + (col + 2);
    int l_index = new_index - 1;
    int r_index = new_index + 1;
    int lu_index = u_index - 1;
    int ru_index = ru_index + 1;
    int ld_index = d_index - 1;
    int rd_index = d_index + 1;

    // calculate current grid point temperature
    double remainting_c = 1 - conductivity[old_index];
    double old_temperature = tinit[new_index];
    double new_temperature = conductivity[old_index] * tinit[new_index] +
                             remainting_c * direct_n_c * (tinit[u_index] + tinit[d_index] + tinit[l_index] + tinit[r_index]) +
                             remainting_c * diagnol_n_c * (tinit[lu_index] + tinit[ru_index] + tinit[ld_index] + tinit[rd_index]);

    new_tinit[new_index] = new_temperature;
    double *tmp = new_tinit;
    new_tinit = tinit;
    tinit = tmp;

    return fabs(old_temperature - new_temperature);
}

void do_compute(const struct parameters *p, struct results *r)
{
    // initialization
    row = p->N;
    col = p->M;
    const double *tinit = p->tinit;
    double tmin = p->io_tmax + 1, tmax = p->io_tmin - 1, tavg, time;
    size_t iter = p->maxiter;
    size_t threshold = p->threshold;
    double *conductivity = p->conductivity;
    double *_tinit = (double *)malloc((row + 2) * (col + 2) * sizeof(double));
    double *new_tinit = (double *)malloc((row + 2) * (col + 2) * sizeof(double));
    double *tmp_ptr;
    struct timespec start, end;
    direct_n_c = sqrt(2) / (sqrt(2) + 1) / 4;
    diagnol_n_c = 1 / (sqrt(2) + 1) / 4;
    // first initailize the actual matrix
    _tinit = _tinit + col + 2;
    new_tinit = new_tinit + col + 2;
    for (int i = 0; i < row; i++)
    {
        int current_row = i * col;
        for (int j = 1; j < col + 2; j++)
        {
            int index = current_row + j;
            _tinit[index] = p->tinit[index];
        }
        _tinit[current_row + 0] = p->tinit[current_row + col - 1];
        _tinit[current_row + col + 1] = p->tinit[current_row];
    }
    // second1 initialize halo grid
    for (int j = 0; j < col + 2; j++)
    {
        _tinit[-1 * col + j] = _tinit[j];
        _tinit[row * col + j] = _tinit[(row - 1) * col + j];

        new_tinit[-1 * col + j] = _tinit[j];
        new_tinit[row * col + j] = _tinit[(row - 1) * col + j];
    }

    // iterator parameters
    bool terminate = false;
    double maxdiff = 0;
    clock_gettime(CLOCK_MONOTONIC, &start); /* time begin */

    for (int _iter = 0; _iter < iter; _iter++)
    {
        if (terminate)
            break; // check convergence

        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < col; c++)
            {
                double diff = calculate_diff(r, c, _tinit, new_tinit, conductivity);
                maxdiff = max(maxdiff, diff);

                if (maxdiff > p->threshold)
                {
                    terminate = false;
                }
            }
        }
    }

    // calculate characteristic values
    int total = row * col;
    double sum = 0.0;
    int idx;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            idx = (i + 1) * (col + 2) + (j + 1);
            tmin = min(tmin, tinit[idx]);
            tmax = max(tmax, tinit[idx]);
            sum += tinit[idx];
        }
    }
    tavg = sum / total;

    clock_gettime(CLOCK_MONOTONIC, &end); /* time end */
    time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);

    // end of calculation, save results
    r->niter = iter;
    r->tmin = tmin;
    r->tmax = tmax;
    r->maxdiff = maxdiff;
    r->tavg = tavg;
    r->time = time;

    // draw_picture(p, tinit, 1);
    // free(tinit);
}