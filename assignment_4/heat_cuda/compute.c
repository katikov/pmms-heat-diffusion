#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cuda_compute.h"
#include "compute.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifdef GEN_PICTURES
static void do_draw(const struct parameters *p,
             size_t key, size_t h, size_t w,
             double (*restrict g)[h][w])
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
             double (*restrict g)[h][w])
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
                        double (*restrict a)[h][w],
                        double (*restrict b)[h][w],
                        double iter,
                        struct timespec *before)
{
    // compute min/max/avg 
    double tmin = INFINITY, tmax = -INFINITY;
    double sum = 0.0;
    double maxdiff = 0.0;
    struct timespec after;

    //We have said that the final reduction does not need to be included.
    clock_gettime(CLOCK_MONOTONIC, &after);
 
    for (size_t i = 1; i < h - 1; ++i)
        for (size_t j = 1; j < w - 1; ++j) 
        {
            double v = (*a)[i][j];
            double v_old = (*b)[i][j];
            sum += v;
            if (tmin > v) tmin = v;
            if (tmax < v) tmax = v;
        }

    r->niter = iter;
    r->tmin = tmin;
    r->tmax = tmax;
    r->tavg = sum / (p->N * p->M);

    r->time = (double)(after.tv_sec - before->tv_sec) + 
        (double)(after.tv_nsec - before->tv_nsec) / 1e9;
 return 0;
}


void do_compute(const struct parameters* p, struct results *r)
{
        // alias input parameters 
    const double (*restrict tinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->tinit;
    const double (*restrict cinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->conductivity;

    // allocate grid data 
    const size_t h = p->N + 2;
    const size_t w = p->M + 2;
    double (*restrict g1)[h][w] = malloc(h * w * sizeof(double));
    double (*restrict g2)[h][w] = malloc(h * w * sizeof(double));

    // allocate halo for conductivities 
    double (*restrict c)[h][w] = malloc(h * w * sizeof(double));

    

    static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
    static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

    // set initial temperatures and conductivities 
    for (int i = 1; i < h - 1; ++i)
        for (int j = 1; j < w - 1; ++j) 
        {
            (*g1)[i][j] = (*tinit)[i-1][j-1];
            (*c)[i][j] = (*cinit)[i-1][j-1];
        }

    // smear outermost row to border 
    for (int j = 1; j < w-1; ++j) {
        (*g1)[0][j] = (*g2)[0][j] = (*g1)[1][j];
        (*g1)[h-1][j] = (*g2)[h-1][j] = (*g1)[h-2][j];
    }

    // compute 
    size_t iter;
    double (*restrict src)[h][w] = g2;
    double (*restrict dst)[h][w] = g1;

    do_copy(h, w, src);
    do_copy(h, w, dst);
    cuda_init(h, w);

    struct timespec before, after;
    clock_gettime(CLOCK_MONOTONIC, &before);

    double *src_device  = createDeviceMemory(src, sizeof(double)*h*w);
    double *dst_device = createDeviceMemory(dst, sizeof(double)*h*w);
    double *c_device = createDeviceMemory(c, sizeof(double)*h*w);
    clock_gettime(CLOCK_MONOTONIC, &after);
    // printf("%f\n", (double)(after.tv_sec - before.tv_sec) + 
    //     (double)(after.tv_nsec - before.tv_nsec) / 1e9);

    for (iter = 1; iter <= p->maxiter; ++iter)
    {
        /* swap source and destination */
        { void *tmp = src_device; src_device = dst_device; dst_device = tmp; }
        double maxdiff = cuda_do_compute_step(h, w, src_device, dst_device, c_device);
        r->maxdiff = maxdiff;
        if(maxdiff<p->threshold){iter++;break;} 

        /* conditional reporting */
        if (iter % p->period == 0) {
            r->niter = iter;
            //fill_report(p, r, h, w, dst, src, iter, &before);
            cuda_fill_report(p, r, h, w, dst_device, &before);
            if(p->printreports) report_results(p, r);
        }
    }
    //getDeviceMemory(src, src_device, sizeof(double)*h*w);
    //getDeviceMemory(dst, dst_device, sizeof(double)*h*w);

    /* report at end in all cases */
    iter--;
    //fill_report(p, r, h, w, dst, src, iter, &before);
    r->niter = iter;
    cuda_fill_report(p, r, h, w, dst_device, &before);

    free(c);
    free(g2);
    free(g1);

    freeDeviceMemory(src_device);
    freeDeviceMemory(dst_device);
    freeDeviceMemory(c_device);
    cuda_finalize();
}



