#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
// #include <avxintrin.h>
// #include <avx2intrin.h>
#include "compute.h"
#include "fail.h"


static inline void get_result(const struct parameters* p, const double* world, struct results *r, struct timespec start_time){
    //tavg, tmax, tmin, time
    //double maxdiff = 0;
    double tsum = 0;
    double tmax = p->io_tmin;
    double tmin = p->io_tmax;

    int N = p->N, M = p->M;
    //int num_cols = (M+1)/2*2+2;
    int num_cols = ((M+3)/4)*4+2;
    for(int i=1;i<=N;i++){
        for(int j=1;j<=M;j++){
            double temp = world[i*num_cols+j];
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
    int num_cols = ((M+3)/4)*4+2;
    //int num_cols = (M+1)/2*2+2;
    begin_picture(key, M, N, p->io_tmin, p->io_tmax);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=M;j++){
            draw_point(j-1,i-1,world[i*num_cols+j]);
        }
    }
    end_picture();

}

void do_compute(const struct parameters* p, struct results *r)
{
    int N = p->N;
    int M = p->M;
    int num_cols = ((M+3)/4)*4+2;  
    //int num_cols = (M+1)/2*2+2; // M/2 upper_bound

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


    if (!(cur_world = calloc((N+2) * num_cols, sizeof(double)))) die("calloc");
    if (!(next_world = calloc((N+2) * num_cols, sizeof(double)))) die("calloc");
    


    for(int j=0;j<M;j++){
        cur_world[j+1] = p->tinit[j];
        next_world[j+1] = p->tinit[j];

        cur_world[(N+1)*num_cols + j+1] = p->tinit[(N-1)*M+j];
        next_world[(N+1)*num_cols+ j+1] = p->tinit[(N-1)*M+j];
    }
    cur_world[0] = next_world[0] = cur_world[M];
    cur_world[M + 1] = next_world[M + 1] = cur_world[1];
    
    cur_world[(N+1)*num_cols] = next_world[(N+1)*num_cols] = cur_world[(N+1)*num_cols+M];
    cur_world[(N+1)*num_cols+M+1] = next_world[(N+1)*num_cols+M+1] = cur_world[(N+1)*num_cols+1];

    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            cur_world[(i+1)*num_cols + j + 1] = p->tinit[i*M + j];              
        }
    }
    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    double maxdiff = p->threshold+1;
    __m256d maxdiff_vec = _mm256_set1_pd(maxdiff);
    for(int iter=0; iter<p->maxiter; iter++){
        r->maxdiff = maxdiff;
        if(maxdiff < p->threshold){
            break;
        }
        maxdiff=0;
        maxdiff_vec = _mm256_set1_pd(0.0);
        if(p->printreports && iter%(p->period)==0 && iter){
            get_result(p, cur_world, r, start_time);
            report_results(p,r);
        }
        for(int i=1;i<=N;i++){
            cur_world[i*num_cols] = cur_world[i*num_cols+M];
            cur_world[i*num_cols+M+1] = cur_world[i*num_cols+1];
        }



        for(int i=1;i<=N;i++){
            double d0 = cur_world[(i-1)*num_cols] + cur_world[(i+1)*num_cols];
            double d1 = cur_world[(i-1)*num_cols+1] + cur_world[(i+1)*num_cols+1];
            int j;
            for(j=1;j<=M-3;j+=4){

                // vectorize with flags -ftree-vectorize
/*                
                double d2 = cur_world[(i-1)*num_cols+j+1] + cur_world[(i+1)*num_cols+j+1];
                double d3 = cur_world[(i-1)*num_cols+j+2] + cur_world[(i+1)*num_cols+j+2];
                double d4 = cur_world[(i-1)*num_cols+j+3] + cur_world[(i+1)*num_cols+j+3];
                double d5 = cur_world[(i-1)*num_cols+j+4] + cur_world[(i+1)*num_cols+j+4];


                double temp0 = cur_world[i*num_cols+j] * conductivity[(i-1)*M+j-1]
                    +    (d0+d2)*weight_diag[(i-1)*M+j-1]
                    +   (d1 + cur_world[i*num_cols+j-1] + cur_world[i*num_cols+j+1])*weight_direct[(i-1)*M+j-1]
                    ;
                double temp1 = cur_world[i*num_cols+j+1] * conductivity[(i-1)*M+j]
                    +    (d1+d3)*weight_diag[(i-1)*M+j]
                    +   (d2 + cur_world[i*num_cols+j] + cur_world[i*num_cols+j+2])*weight_direct[(i-1)*M+j]
                    ;
                double temp2 = cur_world[i*num_cols+j+2] * conductivity[(i-1)*M+j+1]
                    +    (d2+d4)*weight_diag[(i-1)*M+j+1]
                    +   (d3 + cur_world[i*num_cols+j+1] + cur_world[i*num_cols+j+3])*weight_direct[(i-1)*M+j+1]
                    ;
                double temp3 = cur_world[i*num_cols+j+3] * conductivity[(i-1)*M+j+2]
                    +    (d3+d5)*weight_diag[(i-1)*M+j+2]
                    +   (d4 + cur_world[i*num_cols+j+2] + cur_world[i*num_cols+j+4])*weight_direct[(i-1)*M+j+2]
                    ;

                d0 = d4;
                d1 = d5;

                next_world[i*num_cols+j] = temp0;
                next_world[i*num_cols+j+1] = temp1;
                next_world[i*num_cols+j+2] = temp2;
                next_world[i*num_cols+j+3] = temp3;

                maxdiff = (maxdiff >= temp0)? maxdiff:temp0;
                maxdiff = (maxdiff >= temp1)? maxdiff:temp1;
                maxdiff = (maxdiff >= temp2)? maxdiff:temp2;
                maxdiff = (maxdiff >= temp3)? maxdiff:temp3;
            
*/
                __m256d cur = _mm256_loadu_pd(&cur_world[i*num_cols+j]);
                __m256d cond = _mm256_loadu_pd(&conductivity[(i-1)*M+j-1]);
                __m256d sum = _mm256_mul_pd(cur,cond);

                __m256d up = _mm256_loadu_pd(&cur_world[(i-1)*num_cols+j+1]);
                __m256d down = _mm256_loadu_pd(&cur_world[(i+1)*num_cols+j+1]);
                __m256d sum_ud = _mm256_add_pd(up,down);
                __m256d left = _mm256_set_pd(sum_ud[1],sum_ud[0],d1,d0);
                __m256d sum_diag = _mm256_add_pd(left, sum_ud);
                sum = _mm256_fmadd_pd(sum_diag, _mm256_loadu_pd(&weight_diag[(i-1)*M+j-1]),sum);

                //sum = _mm256_add_pd(sum, sum_diag);

                __m256d sum_direct_l = _mm256_loadu_pd(&cur_world[i*num_cols+j-1]);
                __m256d sum_direct_r = _mm256_loadu_pd(&cur_world[i*num_cols+j+1]);
                __m256d sum_direct_lr = _mm256_add_pd(sum_direct_l, sum_direct_r);
                __m256d sum_direct_ud = _mm256_set_pd(sum_ud[2],sum_ud[1],sum_ud[0],d1);
                __m256d sum_direct = _mm256_add_pd(sum_direct_lr, sum_direct_ud);
                sum = _mm256_fmadd_pd(sum_direct, _mm256_loadu_pd(&weight_direct[(i-1)*M+j-1]),sum);

                _mm256_storeu_pd(&next_world[i*num_cols+j],sum);
                d0 = sum_ud[2];
                d1 = sum_ud[3];

                __m256d diff_vec1 = _mm256_sub_pd(cur,sum);
                __m256d diff_vec2 = _mm256_sub_pd(sum,cur);
                __m256d diff_vec = _mm256_max_pd(diff_vec1,diff_vec2);
                maxdiff_vec = _mm256_max_pd(maxdiff_vec, diff_vec);
                
            }
            
            for(;j<=M;j++){
                double d2 = cur_world[(i-1)*num_cols+j+1] + cur_world[(i+1)*num_cols+j+1];
                double temp0 = cur_world[i*num_cols+j] * conductivity[(i-1)*M+j-1]
                    +    (d0+d2)*weight_diag[(i-1)*M+j-1]
                    +   (d1 + cur_world[i*num_cols+j-1] + cur_world[i*num_cols+j+1])*weight_direct[(i-1)*M+j-1]
                    ;
                next_world[i*num_cols+j] = temp0;
                double diff = fabs(cur_world[(i-1)*num_cols+j]-next_world[(i-1)*num_cols+j]);
                maxdiff = (maxdiff>=diff)?maxdiff:diff;
                d0=d1;d1=d2;
            }
            for(int k=0;k<4;k++) maxdiff = (maxdiff>=maxdiff_vec[k])?maxdiff:maxdiff_vec[k];
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
