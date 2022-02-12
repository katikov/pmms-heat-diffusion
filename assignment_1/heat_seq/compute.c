#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include "fail.h"


static inline void get_result(const struct parameters* p, const double* world, struct results *r, struct timespec start_time){
    //tavg, tmax, tmin, time
    //double maxdiff = 0;
    double tsum = 0;
    double tmax = p->io_tmin;
    double tmin = p->io_tmax;

    int N = p->M, M = p->N;
    for(int i=1;i<=N;i++){
        for(int j=1;j<=M;j++){
            double temp = world[i*(M+2)+j];
            //double temp_old = old_world[i*(M+2)+j];
            tsum += temp;
            if(temp>tmax)tmax = temp;
            if(temp<tmin)tmin = temp;
            //double diff = fabs(temp - temp_old);
            //if(diff>maxdiff)maxdiff = diff;
        }
    }
    struct timespec now_time;
    clock_gettime(CLOCK_MONOTONIC, &now_time);
    r->time = now_time.tv_sec - start_time.tv_sec 
        + 1e-9*(now_time.tv_nsec - start_time.tv_nsec);
    //r->maxdiff = maxdiff; 
    r->tmax = tmax;
    r->tmin = tmin;
    r->tavg = tsum / (N*M);
    
}

static inline double get_maxdiff(const struct parameters* p, const double* world, const double* old_world){
    double maxdiff = 0;
    int N = p->M, M = p->N;

    for(int i=1;i<=N;i++){
        for(int j=1;j<=M;j++){
            double temp = world[i*(M+2)+j];
            double temp_old = old_world[i*(M+2)+j];
            double diff = fabs(temp - temp_old);
            maxdiff = (maxdiff>diff)?maxdiff:diff;
        }
    }
    return maxdiff;
    
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
    int N = p->N;
    int M = p->M;
    const double weight_diag = 1/(sqrt(2)+1)/4;
    const double weight_direct = sqrt(2)/(sqrt(2)+1)/4;
    double* cur_world;
    double* next_world;
    double* conductivity;

    if (!(conductivity = calloc((N) * (M), sizeof(double)))) die("calloc");
    for(int i=0;i<N;i++)for(int j=0;j<M;j++){
        conductivity[j*N+i] = p->conductivity[i*M+j];
    }

    
    if (!(cur_world = calloc((N+2) * (M+2), sizeof(double)))) die("calloc");
    if (!(next_world = calloc((N+2) * (M+2), sizeof(double)))) die("calloc");
    
    // TODO: transpose the matrix
    for(int j=0;j<p->M;j++){
        cur_world[(j+1)*(p->N + 2)] = p->tinit[j];
        next_world[(j+1)*(p->N + 2)] = p->tinit[j];

        cur_world[(j+2)*(p->N + 2)-1] = p->tinit[(p->N-1)*(p->M)+j];
        next_world[(j+2)*(p->N + 2)-1] = p->tinit[(p->N-1)*(p->M)+j];
    }
    cur_world[0] = next_world[0] = cur_world[(p->M)*(p->N + 2)];
    cur_world[(p->M+1)*(p->N + 2)] = next_world[(p->M+1)*(p->N + 2)] = cur_world[(p->N + 2)];
    cur_world[p->N + 1] = next_world[p->N + 1] = cur_world[(p->M + 1)*(p->N + 2)-1];
    cur_world[(p->M + 2)*(p->N + 2)-1] = next_world[(p->M + 2)*(p->N + 2)-1] = cur_world[2*(p->N + 2)-1];
    for(int i=0;i<p->N;i++){
        for(int j=0;j<p->M;j++){
            cur_world[(j+1)*(p->N+2) + i + 1] = p->tinit[i*p->M + j];              
        }
    }

    // for(int j=0;j<M;j++){
    //     cur_world[j+1] = p->tinit[j];
    //     next_world[j+1] = p->tinit[j];

    //     cur_world[(N+1)*(M + 2) + j+1] = p->tinit[(N-1)*M+j];
    //     next_world[(N+1)*(M + 2)+ j+1] = p->tinit[(N-1)*M+j];
    // }
    // cur_world[0] = next_world[0] = cur_world[M];
    // cur_world[M + 1] = next_world[M + 1] = cur_world[1];
    
    // cur_world[(N+1)*(M + 2)] = next_world[(N+1)*(M + 2)] = cur_world[(N+2)*(M+2)-2];
    // cur_world[(N+2)*(M+2)-1] = next_world[(N+2)*(M+2)-1] = cur_world[(N+1)*(M+2)+1];

    // for(int i=0;i<N;i++){
    //     for(int j=0;j<M;j++){
    //         cur_world[(i+1)*(M+2) + j + 1] = p->tinit[i*M + j];              
    //     }
    // }
    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    

    for(int iter=0; iter<p->maxiter; iter++){
        double maxdiff = get_maxdiff(p,cur_world, next_world);
        r->maxdiff = maxdiff;
        if(maxdiff < p->threshold){
            break;
        }
        if(p->printreports && iter%(p->period)==0 && iter){
            get_result(p, cur_world, r, start_time);
            report_results(p,r);
        }
        for(int i=1;i<=N;i++){
            // cur_world[i*(M+2)] = cur_world[(i+1)*(M+2)-2];
            // cur_world[(i+1)*(M+2)-1] = cur_world[i*(M+2)+1];
            cur_world[i] = cur_world[M*(N+2)+i];
            cur_world[(M+1)*(N+2)+i] = cur_world[(N+2)+i];
        }

        for(int i=1;i<=M;i++){
            for(int j=1;j<=N;j++){
                next_world[i*(N+2)+j] = cur_world[i*(N+2)+j] * conductivity[(i-1)*N+j-1]
                    + (     (cur_world[(i-1)*(N+2)+j-1] + cur_world[(i-1)*(N+2)+j+1] 
                            + cur_world[(i+1)*(N+2)+j-1] + cur_world[(i+1)*(N+2)+j+1])*weight_diag
                        +   (cur_world[(i-1)*(N+2)+j] + cur_world[i*(N+2)+j-1]
                            + cur_world[i*(N+2)+j+1] + cur_world[(i+1)*(N+2)+j])*weight_direct
                    ) 
                    * (1-conductivity[(i-1)*N+j-1]);
            }
        }

        double* temp_world = cur_world;
        cur_world = next_world;
        next_world = temp_world;

        r->niter = iter + 1;
    }
    
    r->maxdiff = get_maxdiff(p,cur_world, next_world);
    get_result(p, cur_world, r, start_time);
    
    // export pgm file
    //draw_picture(p, cur_world, 1);
}
