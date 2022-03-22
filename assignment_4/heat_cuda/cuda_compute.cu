#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <time.h>

#include "input.h"
#include "output.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

__constant__ __device__ double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
__constant__ __device__ double c_cdiag = 0.25 / (M_SQRT2 + 1.0);
static double* maxdiff_block;
static double* min_block;
static double* max_block;
static double* sum_block;
const int threadBlockSize = 128;

static void checkCudaCall(cudaError_t result) {
    if (result != cudaSuccess) {
        printf("cuda error \n");
        exit(1);
    }
}

extern "C" 
void* createDeviceMemory(void* p, size_t length){
    double* res;
    checkCudaCall(cudaMalloc((void **) &res, length));
    if(res == NULL){
        printf("cudaMalloc error\n");
        exit(1);
    }
    checkCudaCall(cudaMemcpy(res, p, length, cudaMemcpyHostToDevice));
    return res;
}

extern "C" 
void getDeviceMemory(void* host, void* device, size_t length){
    checkCudaCall(cudaMemcpy(host, device, length, cudaMemcpyDeviceToHost));
}

extern "C" 
void freeDeviceMemory(void* p){
    checkCudaCall(cudaFree(p));
}



__global__ void vectorAddKernel(float* deviceA, float* deviceB, float* deviceResult) {
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
// insert operation here
    deviceResult[i] = deviceA[i]+deviceB[i];
}

__global__ void cudaCopyKernel(unsigned h, unsigned w,
             double* g)
{
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < h){
        g[i*w+w-1] = g[i*w+1];
        g[i*w] = g[i*w+w-2];
    }
}

__global__ void cudaReduceKernel(unsigned h, unsigned w, double* data, 
                            double* global_min, double* global_max, double* global_sum) {
    unsigned tid = threadIdx.x;
    unsigned i = blockIdx.x+1;
    unsigned j = blockIdx.y*blockDim.x + threadIdx.x + 1;
    unsigned idx = i*w+j;
    __shared__ double min[threadBlockSize];
    __shared__ double max[threadBlockSize];
    __shared__ double sum[threadBlockSize];
    min[tid] = INFINITY; max[tid] = -INFINITY; sum[tid] = 0;
    
    if(j<w-1){
        min[tid] = data[idx];
        max[tid] = data[idx];
        sum[tid] = data[idx];
    }
    __syncthreads();
    for(unsigned s = blockDim.x/2; s>0; s/=2){
        if(tid<s){
            max[tid] = MAX(max[tid], max[tid+s]);
            min[tid] = MIN(min[tid], min[tid+s]);
            sum[tid] += sum[tid+s];
        }
        __syncthreads();
    }

    // for(unsigned s=1; s<blockDim.x; s*=2){
    //     if(tid%(2*s)==0){
    //         max[tid] = MAX(max[tid], max[tid+s]);
    //         min[tid] = MIN(min[tid], min[tid+s]);
    //         sum[tid] += sum[tid+s];
    //     }
    //     __syncthreads();
    // }
    
    if(tid==0){
        unsigned bid = ((w-2+threadBlockSize-1)/threadBlockSize)*blockIdx.x + blockIdx.y;
        global_min[bid] = min[0];
        global_max[bid] = max[0];
        global_sum[bid] = sum[0];
    }

}

__global__ void cudaMaxdiffKernel(int step, int max_id, double* data){
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i<max_id)
        data[i] = MAX(data[i], data[i+step]);
}

// __global__ void cudaPrecomputeKernel(unsigned w, 
//                         double* src, double* temp_result)
// {
//     //size_t tid = threadIdx.x;
//     unsigned i = blockIdx.x+1;
//     unsigned j = blockIdx.y*blockDim.x + threadIdx.x;
//     unsigned idx = i*w+j;
//     if(j<w){
//         temp_result[idx] = src[idx-w] + src[idx+w];
//     }
// }


__global__ void cudaComputeKernel(unsigned w, 
                        double* src, double* dst, double* c, double* global_maxdiff)
{
    unsigned tid = threadIdx.x;
    unsigned i = blockIdx.x+1;
    unsigned j = blockIdx.y*blockDim.x + threadIdx.x + 1;
    unsigned idx = i*w+j;
    __shared__ double diff[threadBlockSize];
    diff[tid] = 0.0;

    if(j<w-1){

        double weight = c[idx];
        double restw = 1.0 - weight;

        // dst[idx] = weight * src[idx] +
        //      (temp_result[idx] + src[idx+1] + src[idx-1]) * (restw * c_cdir) +
        //      (temp_result[idx-1] + temp_result[idx+1]) * (restw * c_cdiag);

        dst[idx] = weight * src[idx] +
            (src[idx+w] + src[idx-w] + src[idx+1] + src[idx-1]) * (restw * c_cdir) +

            (src[idx-w-1] + src[idx-w+1] + src[idx+w-1] + src[idx+w+1]) * (restw * c_cdiag);
        
        diff[tid] = fabs(dst[idx] - src[idx]);
    }
    
    __syncthreads();

    for(unsigned s = blockDim.x/2; s>0; s/=2){
        if(tid<s){
            diff[tid] = MAX(diff[tid], diff[tid+s]);
        }
        __syncthreads();
    }

    // for(unsigned s=1; s<blockDim.x; s*=2){
    //     if(tid%(2*s)==0)
    //         diff[tid] = MAX(diff[tid], diff[tid+s]);
    //     __syncthreads();
    // }

    if(tid==0){
        unsigned blocksPerLine = (w-2+threadBlockSize-1)/threadBlockSize;
        global_maxdiff[blockIdx.x * blocksPerLine + blockIdx.y] = diff[0];
    }
    
}

extern "C"
void cuda_init(size_t h, size_t w){
    size_t blocksPerLine = (w-2+threadBlockSize-1)/threadBlockSize;
    size_t num_blocks = (h-2)*blocksPerLine;
    //checkCudaCall(cudaMalloc((void **) &maxdiff_block, h*w * sizeof(double)));
    checkCudaCall(cudaMalloc((void **) &maxdiff_block, num_blocks * sizeof(double)));
    checkCudaCall(cudaMalloc((void **) &min_block, num_blocks * sizeof(double)));
    checkCudaCall(cudaMalloc((void **) &max_block, num_blocks * sizeof(double)));
    checkCudaCall(cudaMalloc((void **) &sum_block, num_blocks * sizeof(double)));
}

extern "C"
void cuda_finalize(){
    checkCudaCall(cudaFree(maxdiff_block));
    checkCudaCall(cudaFree(min_block));
    checkCudaCall(cudaFree(max_block));
    checkCudaCall(cudaFree(sum_block));
}

extern "C"
double cuda_do_compute_step(size_t h, size_t w, 
                        double* src, double* dst, double* c){

        
    cudaCopyKernel<<<(h+threadBlockSize-1)/threadBlockSize, threadBlockSize>>>(h, w, src);
    size_t blocksPerLine = (w-2+threadBlockSize-1)/threadBlockSize;
    dim3 numBlocks(h-2,blocksPerLine);
    dim3 numBlocks2(h-2, (w+threadBlockSize-1)/threadBlockSize);

    ///cudaPrecomputeKernel<<<numBlocks2 ,threadBlockSize>>>(w, src, temp_sum);

    
    cudaComputeKernel<<<numBlocks, threadBlockSize>>>(w, src, dst, c, maxdiff_block);  
    

    int len;
    for(len=blocksPerLine*(h-2); len>=threadBlockSize*8; len=(len+1)/2){
        cudaMaxdiffKernel<<<(len+threadBlockSize-1)/threadBlockSize, threadBlockSize>>>((len+1)/2, len/2, maxdiff_block);
    }
    cudaDeviceSynchronize();

    
    // double r[blocksPerLine*(h-2)];
    // getDeviceMemory(r, maxdiff_block, blocksPerLine*(h-2)*sizeof(double));
    // double res=0;
    // for(int i=0;i<blocksPerLine*(h-2);i++){
    //     res = MAX(res, r[i]);
    // }

    double r[len];
    getDeviceMemory(r, maxdiff_block, len*sizeof(double));
    double res=0;
    for(int i=0;i<len;i++){
        res = MAX(res, r[i]);
    }
    return res;
    
}

extern "C"
void cuda_fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, double* data,
                        struct timespec *before)
{
    // compute min/max/avg 
    double tmin = INFINITY, tmax = -INFINITY;
    double sum = 0.0;
    struct timespec after;

    size_t blocksPerLine = (w-2+threadBlockSize-1)/threadBlockSize;
    dim3 numBlocks(h-2,blocksPerLine);

    //We have said that the final reduction does not need to be included.
    clock_gettime(CLOCK_MONOTONIC, &after);
    
    cudaReduceKernel<<<numBlocks, threadBlockSize>>>(h, w, data, min_block, max_block, sum_block);
    cudaDeviceSynchronize();

    double recv[blocksPerLine*(h-2)];
    getDeviceMemory(recv, min_block, blocksPerLine*(h-2)*sizeof(double));
    for(int i=0;i<blocksPerLine*(h-2);i++){
        tmin = MIN(tmin, recv[i]);
    }

    getDeviceMemory(recv, max_block, blocksPerLine*(h-2)*sizeof(double));
    for(int i=0;i<blocksPerLine*(h-2);i++){
        tmax = MAX(tmax, recv[i]);
    }

    getDeviceMemory(recv, sum_block, blocksPerLine*(h-2)*sizeof(double));
    for(int i=0;i<blocksPerLine*(h-2);i++){
        sum += recv[i];
    }
    r->tmin = tmin;
    r->tmax = tmax;
    r->tavg = sum / (p->N * p->M);
    
    r->time = (double)(after.tv_sec - before->tv_sec) + 
        (double)(after.tv_nsec - before->tv_nsec) / 1e9;
    
}

