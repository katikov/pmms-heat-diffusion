#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include <pthread.h>
#include <stdio.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// extern int omp_heat_chunk_size;
// extern int omp_heat_parallel_type;
// extern int omp_heat_percentual;
// extern int omp_heat_columns;
// extern int omp_heat_pure_seq;

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

static void inline do_copy_parallel(size_t h, size_t w,
             double (*restrict g)[h][w], size_t s, size_t e)
{
    size_t i;

    // copy left and right column to opposite border 
    for (i = s; i <= e; ++i) {
        (*g)[i][w-1] = (*g)[i][1];
        (*g)[i][0] = (*g)[i][w-2];
    }
}


// Does the reduction step and return if the convergence has setteled 
// TODO: parallel fill
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


static int proc_count;
// h, w, c_cdir, c_cdiag, r,src,dst,p,c
static struct timespec before;
pthread_mutex_t mutex_diff = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barrier;
static double _maxdiff;
static size_t _h, _w;
static double _c_cdir,_c_cdiag;
static void* _src, *_dst, *_p, *_c, *_r;
static long _iter = 1;

pthread_cond_t barrier_cond = PTHREAD_COND_INITIALIZER;
static int barrier_count=0;
pthread_spinlock_t spin;

static int inline fill_report_parallel(const struct parameters *p, struct results *r,
                        size_t h, size_t w, 
                        double (*restrict a)[h][w],
                        double (*restrict b)[h][w],
                        double iter,
                        struct timespec *before,
                        size_t s, size_t e)
{
    // compute min/max/avg 
    double tmin = INFINITY, tmax = -INFINITY;
    double sum = 0.0;
    double maxdiff = 0.0;
    struct timespec after;

    for (size_t i = s; i <= e; ++i)
        for (size_t j = 1; j < w - 1; ++j) 
        {
            double v = (*a)[i][j];
            double v_old = (*b)[i][j];
            sum += v;
            if (tmin > v) tmin = v;
            if (tmax < v) tmax = v;
        }

    pthread_mutex_lock(&mutex_diff);
    r->tmin = MIN(r->tmin,tmin);
    r->tmax = MAX(r->tmax,tmax);
    r->tavg += sum / (p->N * p->M);
    pthread_mutex_unlock(&mutex_diff);

 return 0;
}


void* do_compute_parallel(void* params){
    const struct parameters* p = _p;
    const size_t h = _h, w = _w;
    double (*restrict src)[h][w] = _src;
    double (*restrict dst)[h][w] = _dst;

    double (*restrict c)[h][w] = _c;

    const double c_cdir = _c_cdir, c_cdiag = _c_cdiag;
    struct results *r = _r;

    int tid = *(int *)params;
    int start = tid*(p->N)/proc_count+1;
    int end = (tid+1)*(p->N)/proc_count;
    
    size_t iter;
//printf("%d %d\n",start, end);
    for (iter = 1; iter <= p->maxiter; ++iter)
    {
        { void *tmp = src; src = dst; dst = tmp; }

        // initialize halo on source
        //do_copy(h, w, src);
        do_copy_parallel(h, w, src, start, end);

        double maxdiff = 0.0;
        pthread_barrier_wait(&barrier);
        if(tid==0) _maxdiff = 0.0;
        pthread_barrier_wait(&barrier);
        // compute 
        
        for (int i = start; i <= end; ++i) {
            for (int j = 1; j < w - 1; ++j)
            {
                double w = (*c)[i][j];
                double restw = 1.0 - w;

                (*dst)[i][j] = w * (*src)[i][j] +

                    ((*src)[i + 1][j] + (*src)[i - 1][j] +
                        (*src)[i][j + 1] + (*src)[i][j - 1]) * (restw * c_cdir) +

                    ((*src)[i - 1][j - 1] + (*src)[i - 1][j + 1] +
                        (*src)[i + 1][j - 1] + (*src)[i + 1][j + 1]) * (restw * c_cdiag);

                double diff = fabs((*dst)[i][j] - (*src)[i][j]);
                if (diff > maxdiff) maxdiff = diff;
            }
        }
        pthread_mutex_lock(&mutex_diff);
        _maxdiff = MAX(_maxdiff, maxdiff);
        pthread_mutex_unlock(&mutex_diff);
        if (tid==0 && iter % p->period == 0) {
            r->niter = iter;
            r->tmin = INFINITY;
            r->tmax = -INFINITY;
            r->tavg = 0;
            struct timespec after;
            clock_gettime(CLOCK_MONOTONIC, &after);
            r->time = (double)(after.tv_sec - before.tv_sec) + 
                 (double)(after.tv_nsec - before.tv_nsec) / 1e9;
        }
        pthread_barrier_wait(&barrier);
        if(tid==0)
            r->maxdiff = _maxdiff;
        if(_maxdiff<p->threshold){iter++;break;} 

        // conditional reporting 
        if (iter % p->period == 0) {
           fill_report_parallel(p, r, h, w, dst, src, iter, &before, start, end);
           pthread_barrier_wait(&barrier);
           //fill_report(p, r, h, w, dst, src, iter, &before);
            if(tid == 0 && p->printreports) report_results(p, r);
        }
        
        #ifdef GEN_PICTURES
        do_draw(p, iter, h, w, src);
        #endif

        // swap source and destination

    }

    if(tid==0)_iter = iter;
    return NULL;
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

    (*src)[0][w-1] = (*src)[0][1];
    (*src)[0][0] = (*src)[0][w-2];
    (*src)[h-1][w-1] = (*src)[h-1][1];
    (*src)[h-1][0] = (*src)[h-1][w-2];

    (*dst)[0][w-1] = (*dst)[0][1];
    (*dst)[0][0] = (*dst)[0][w-2];
    (*dst)[h-1][w-1] = (*dst)[h-1][1];
    (*dst)[h-1][0] = (*dst)[h-1][w-2];

    proc_count = MAX(p->nthreads, 1);

    _r = r; _src = src; _dst = dst;
    _p = (void *)p; _c = c; _h = h; _w = w;
    _c_cdiag = c_cdiag; _c_cdir = c_cdir;

    int params[proc_count];
    for(int i=0;i<proc_count;i++) params[i]=i;
    pthread_t thread_ids[proc_count];

    //{ void *tmp = src; src = dst; dst = tmp; }

    // initialize halo on source
    //do_copy(h, w, src);
    pthread_barrier_init(&barrier, NULL, proc_count);
    // pthread_attr_t attr;
    // pthread_attr_init ( & attr );
    // pthread_attr_setscope ( & attr , PTHREAD_SCOPE_SYSTEM );


    clock_gettime(CLOCK_MONOTONIC, &before);
    for(int i=0;i<proc_count;i++){
        pthread_create(&thread_ids[i], NULL, &do_compute_parallel, &params[i]);
    }
    for(int i=0;i<proc_count;i++){
        void* result;
        pthread_join(thread_ids[i], &result);
    }

    // report at end in all cases 
    _iter--;
    if(_iter%2){
        void *tmp = src; src = dst; dst = tmp; 
    }
    fill_report(p, r, h, w, dst, src, _iter, &before);

    free(c);
    free(g2);
    free(g1);
}
