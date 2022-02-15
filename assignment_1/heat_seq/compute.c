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

    int N = p->N, M = p->M;
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


static void draw_picture(const struct parameters* p, const double* world, size_t key){
    int N = p->N, M = p->M;
    begin_picture(key, M, N, p->io_tmin, p->io_tmax);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=M;j++){
            draw_point(j-1,i-1,world[i*(M+2)+j]);
        }
    }
    end_picture();

}

void do_compute(const struct parameters* p, struct results *r)
{
    int N = p->N;
    int M = p->M;
    const double weight_diag_base = 1/(sqrt(2)+1)/4;
    const double weight_direct_base = sqrt(2)/(sqrt(2)+1)/4;
    double* cur_world;
    double* next_world;
    const double* conductivity = p->conductivity;

    double* weight_diag;
    double* weight_direct;
    if (!(weight_diag = calloc((N) * (M), sizeof(double)))) die("calloc");
    if (!(weight_direct = calloc((N) * (M), sizeof(double)))) die("calloc");
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            weight_diag[i*M+j] = weight_diag_base * (1-conductivity[i*M+j]);
            weight_direct[i*M+j] = weight_direct_base * (1-conductivity[i*M+j]);
        }
    }

    // double* sum_ud;
    // if (!(sum_ud = calloc((N+2) * (M+2), sizeof(double)))) die("calloc");

    if (!(cur_world = calloc((N+2) * (M+2), sizeof(double)))) die("calloc");
    if (!(next_world = calloc((N+2) * (M+2), sizeof(double)))) die("calloc");
    

    for(int j=0;j<M;j++){
        cur_world[j+1] = p->tinit[j];
        next_world[j+1] = p->tinit[j];

        cur_world[(N+1)*(M + 2) + j+1] = p->tinit[(N-1)*M+j];
        next_world[(N+1)*(M + 2)+ j+1] = p->tinit[(N-1)*M+j];
    }
    cur_world[0] = next_world[0] = cur_world[M];
    cur_world[M + 1] = next_world[M + 1] = cur_world[1];
    
    cur_world[(N+1)*(M + 2)] = next_world[(N+1)*(M + 2)] = cur_world[(N+2)*(M+2)-2];
    cur_world[(N+2)*(M+2)-1] = next_world[(N+2)*(M+2)-1] = cur_world[(N+1)*(M+2)+1];

    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            cur_world[(i+1)*(M+2) + j + 1] = p->tinit[i*M + j];              
        }
    }
    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    double maxdiff = p->threshold+1;
    for(int iter=0; iter<p->maxiter; iter++){
        //double maxdiff = get_maxdiff(p,cur_world, next_world);
        r->maxdiff = maxdiff;
        if(maxdiff < p->threshold){
            break;
        }
        maxdiff=0;
        if(p->printreports && iter%(p->period)==0 && iter){
            get_result(p, cur_world, r, start_time);
            report_results(p,r);
        }
        for(int i=1;i<=N;i++){
            cur_world[i*(M+2)] = cur_world[(i+1)*(M+2)-2];
            cur_world[(i+1)*(M+2)-1] = cur_world[i*(M+2)+1];
        }

        for(int i=1;i<=N;i++){
            double d1 = cur_world[(i-1)*(M+2)] + cur_world[(i+1)*(M+2)];
            double d2 = cur_world[(i-1)*(M+2)+1] + cur_world[(i+1)*(M+2)+1];
            for(int j=1;j<=M;j++){
                double d3 = cur_world[(i-1)*(M+2)+j+1] + cur_world[(i+1)*(M+2)+j+1];
                double temp = cur_world[i*(M+2)+j] * conductivity[(i-1)*M+j-1]
                    +    (d1+d3)*weight_diag[(i-1)*M+j-1]
                    +   (d2 + cur_world[i*(M+2)+j-1] + cur_world[i*(M+2)+j+1])*weight_direct[(i-1)*M+j-1]
                    ;
                d1 = d2;
                d2 = d3;
                double diff = fabs(temp - cur_world[i*(M+2)+j]);
                maxdiff = (maxdiff>diff)?maxdiff:diff;
                next_world[i*(M+2)+j] = temp;
            }
        }

        double* temp_world = cur_world;
        cur_world = next_world;
        next_world = temp_world;

        r->niter = iter + 1;
    }
    
    //r->maxdiff = get_maxdiff(p,cur_world, next_world);
    r->maxdiff = maxdiff;
    get_result(p, cur_world, r, start_time);
    free(weight_diag);
    free(weight_direct);
    free(next_world);
    free(cur_world);
    
    // export pgm file
    //draw_picture(p, cur_world, 1);
}