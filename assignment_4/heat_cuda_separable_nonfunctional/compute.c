#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "cuda_compute.h"
#include "compute.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define M_SQRT2    1.41421356237309504880

#ifdef FAST
typedef float real;
#else
typedef double real;
#endif

#ifdef GEN_PICTURES
static void do_draw(const struct parameters *p,
             size_t key, size_t h, size_t w,
             real (*restrict g)[h][w])
{
    begin_picture(key, w-2, h-2, p->io_tmin, p->io_tmax);
    size_t i, j;
    for (i = 1; i < h-1; ++i)
        for (j = 1; j < w-1; ++j)
            draw_point(j-1, i-1, (*g)[i][j]);
    end_picture();
}
#endif

static void inline do_copy(size_t h, size_t w,
             real (*restrict g)[h][w])
{
    size_t i;

    // copy left and right column to opposite border 
    for (i = 0; i < h; ++i) {
        (*g)[i][w-1] = (*g)[i][1];
        (*g)[i][0] = (*g)[i][w-2];
    }
}

// Does the reduction step and return if the convergence has setteled 
static int fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, 
                        real (*restrict a)[h][w],
                        real (*restrict b)[h][w],
                        real iter,
                        struct timespec *before)
{
    // compute min/max/avg 
    real tmin = INFINITY, tmax = -INFINITY;
    real sum = 0.0;
    real maxdiff = 0.0;
    struct timespec after;

    //We have said that the final reduction does not need to be included.
    clock_gettime(CLOCK_MONOTONIC, &after);
 
    for (size_t i = 1; i < h - 1; ++i)
        for (size_t j = 1; j < w - 1; ++j) 
        {
            real v = (*a)[i][j];
            real v_old = (*b)[i][j];
            sum += v;
            if (tmin > v) tmin = v;
            if (tmax < v) tmax = v;
        }

    r->niter = iter;
    r->tmin = tmin;
    r->tmax = tmax;
    r->tavg = sum / (p->N * p->M);

    r->time = (real)(after.tv_sec - before->tv_sec) + 
        (real)(after.tv_nsec - before->tv_nsec) / 1e9;
 return 0;
}


void do_compute(const struct parameters* p, struct results *r)
{
    // alias input parameters 
    const double (*restrict tinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->tinit;
    const double (*restrict cinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->conductivity;

    // transpose matrix for optimal coalesced copy parallelism
    // allocate grid data 
    const size_t h = p->M + 2;
    const size_t w = p->N + 2;
    real (*restrict h_current)[h][w] = malloc(h * w * sizeof(real));
    real (*restrict h_next)[h][w] = malloc(h * w * sizeof(real));

    // allocate halo for conductivities 
    real (*restrict cond)[h][w] = malloc(h * w * sizeof(real));
    real(*restrict cond_p)[h][w] = malloc(h * w * sizeof(real));

    static const real c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
    static const real c_cdiag = 0.25 / (M_SQRT2 + 1.0);

    // set initial temperatures and conductivities 
    for (int i = 1; i < h - 1; ++i)
        for (int j = 1; j < w - 1; ++j) 
        {
            (*h_current)[i][j] = (*tinit)[j-1][i-1];
            real c = (*cinit)[j - 1][i - 1];
            (*cond)[i][j] = c - (1-c)* c_cdiag * c_cdiag /c_cdir;
            (*cond_p)[i][j] = 1-c;
        }

    // smear outermost column (row pre-transposition) to border 
    for (int i = 1; i < h-1; ++i) {
        (*h_current)[i][0] = (*h_next)[i][0] = (*h_current)[i][1];
        (*h_current)[i][-1] = (*h_next)[i][w-1] = (*h_current)[i][w-2];
    }
    // loop outermost row to innermost
    memcpy(&(*h_current)[0][0], &(*h_current)[h-2][0], w*sizeof(real));
    memcpy(&(*h_current)[h-1][0], &(*h_current)[1][0], w * sizeof(real));
    warmup();

    // compute 
    size_t iter;
    cuda_init(h, w);
    struct timespec before, after;
    clock_gettime(CLOCK_MONOTONIC, &before);

    real* d_current  = createDeviceMemory(h_current, sizeof(real)*h*w);
    real* d_next = createDeviceMemory(h_next, sizeof(real)*h*w);
    real* d_temp = createDeviceMemory(NULL, sizeof(real) * h * w);
    real* d_cond = createDeviceMemory(cond, sizeof(real)*h*w);
    real* d_cond_p = createDeviceMemory(cond_p, sizeof(real) * h * w);

    // pre-emptive swap that gets undone on first loop iteration
    { void* tmp = d_current; d_current = d_next; d_next = tmp; }

    clock_gettime(CLOCK_MONOTONIC, &after);
    // printf("%f\n", (real)(after.tv_sec - before.tv_sec) + 
    //     (real)(after.tv_nsec - before.tv_nsec) / 1e9);

    for (iter = 1; iter <= p->maxiter; ++iter)
    {
        /* swap source and destination */
        { void *tmp = d_current; d_current = d_next; d_next = tmp; }
        real maxdiff = cuda_do_compute_step(h, w, d_current, d_next, d_temp, d_cond, d_cond_p);
        r->maxdiff = maxdiff;
        if(maxdiff<p->threshold){iter++;break;} 

        /* conditional reporting */
        if (iter % p->period == 0) {
            r->niter = iter;
            //fill_report(p, r, h, w, dst, src, iter, &before);
            cuda_fill_report(p, r, h, w, d_next, &before);
            if(p->printreports) report_results(p, r);
        }
    }
    //getDeviceMemory(src, src_device, sizeof(real)*h*w);
    //getDeviceMemory(dst, dst_device, sizeof(real)*h*w);

    /* report at end in all cases */
    iter--;
    //fill_report(p, r, h, w, dst, src, iter, &before);
    r->niter = iter;
    cuda_fill_report(p, r, h, w, d_next, &before);

    free(h_current);
    free(h_next);
    free(cond);
    free(cond_p);

    freeDeviceMemory(d_current);
    freeDeviceMemory(d_next);
    freeDeviceMemory(d_cond);
    freeDeviceMemory(d_cond_p);
    cuda_finalize();
}



