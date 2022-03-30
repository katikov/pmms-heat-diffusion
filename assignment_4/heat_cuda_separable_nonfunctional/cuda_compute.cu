#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <time.h>

#include "input.h"
#include "output.h"

#ifdef FAST
typedef float real;
#else
typedef double real;
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define M_SQRT2    1.41421356237309504880

__constant__ __device__ real c_cdir = M_SQRT2 / (M_SQRT2 + 1.0);
__constant__ __device__ real c_cdiag = (M_SQRT2 + 1.0);

static real* maxdiff_block;
static real* min_block;
static real* max_block;
static real* sum_block;
const int threadBlockSize = 128;
const int blockLimit = 65535;

static void checkCudaCall(cudaError_t result) {
    if (result != cudaSuccess) {
        printf("cuda error \n");
        exit(1);
    }
}

extern "C" 
void* createDeviceMemory(void* p, size_t length){
    real* res;
    checkCudaCall(cudaMalloc((void **) &res, length));
    if(res == NULL){
        printf("cudaMalloc error\n");
        exit(1);
    }
    if (p) checkCudaCall(cudaMemcpy(res, p, length, cudaMemcpyHostToDevice));
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

__global__ void cudaCopyKernel(unsigned h, unsigned w,
             real* g)
{
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < h){
        g[i*w+w-1] = g[i*w+1];
        g[i*w] = g[i*w+w-2];
    }
}

__global__ void cudaReduceKernel(unsigned h, unsigned w, real* data, 
                            real* global_min, real* global_max, real* global_sum) {
    unsigned tid = threadIdx.x;
    unsigned i = blockIdx.x+1;
    unsigned j = blockIdx.y*blockDim.x + threadIdx.x + 1;
    unsigned idx = i*w+j;
    __shared__ real min[threadBlockSize];
    __shared__ real max[threadBlockSize];
    __shared__ real sum[threadBlockSize];
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

__global__ void cudaMaxdiffKernel(int step, int max_id, real* data){
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i<max_id)
        data[i] = MAX(data[i], data[i+step]);
}

// __global__ void cudaPrecomputeKernel(unsigned w, 
//                         real* src, real* temp_result)
// {
//     //size_t tid = threadIdx.x;
//     unsigned i = blockIdx.x+1;
//     unsigned j = blockIdx.y*blockDim.x + threadIdx.x;
//     unsigned idx = i*w+j;
//     if(j<w){
//         temp_result[idx] = src[idx-w] + src[idx+w];
//     }
// }

__global__ void warm_up_gpu() {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float ia, ib;
    ia = ib = 0.0f;
    ib += ia + tid;
}

extern "C"
void warmup() {
    warm_up_gpu <<< blockLimit, threadBlockSize >>> ();
}

__global__ void cudaComputeKernel(unsigned w, unsigned h,
                        real* src, real* dst, real* c, real* global_maxdiff)
{
    unsigned tid = threadIdx.x;
    unsigned i = blockIdx.x+1;
    unsigned j = blockIdx.y*blockDim.x + threadIdx.x + 1;
    unsigned idx = i*w+j;
    __shared__ real diff[threadBlockSize];
    diff[tid] = 0.0;

    if(j<w-1){

        real weight = c[idx];
        real restw = 1.0 - weight;

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

// =====================================================================================
// =====================================================================================
// =====================================================================================
// Taken from NVIDIA White Paper on separable convolutions and adapted
// https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_64_website/projects/convolutionSeparable/doc/convolutionSeparable.pdf
/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

////////////////////////////////////////////////////////////////////////////////
// Common host and device functions
////////////////////////////////////////////////////////////////////////////////

//Round a / b to nearest higher integer value
int iDivUp(int a, int b) {
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Round a / b to nearest lower integer value
int iDivDown(int a, int b) {
    return a / b;
}

//Align a to nearest higher multiple of b
int iAlignUp(int a, int b) {
    return (a % b != 0) ? (a - a % b + b) : a;
}

//Align a to nearest lower multiple of b
int iAlignDown(int a, int b) {
    return a - a % b;
}

////////////////////////////////////////////////////////////////////////////////
// Loop unrolling templates, needed for best performance
////////////////////////////////////////////////////////////////////////////////
#ifdef FAST
#define REAL_SIZE 1
#else
#define REAL_SIZE 2
#endif
#define KERNEL_RADIUS 1
#define KERNEL_W (2 * KERNEL_RADIUS + 1)
#define IMUL(a, b) __mul24(a, b)
#define ROW_TILE_W 128/REAL_SIZE
#define KERNEL_RADIUS_ALIGNED 16 // Half warps are 16 elemenets big
#define COLUMN_TILE_W 16/REAL_SIZE
#define COLUMN_TILE_H 48/REAL_SIZE

__global__ void cudaRowComputeKernel(
    unsigned dataW,
    unsigned dataH,
    real* d_Result,
    real* d_Data
){
        //Data cache
        __shared__ real data[KERNEL_RADIUS + ROW_TILE_W + KERNEL_RADIUS];

        //Current tile and apron limits, relative to row start
        const int         tileStart = IMUL(blockIdx.x, ROW_TILE_W);
        const int           tileEnd = tileStart + ROW_TILE_W - 1;
        const int        apronStart = tileStart - KERNEL_RADIUS;
        const int          apronEnd = tileEnd + KERNEL_RADIUS;

        //Clamp tile and apron limits by image borders
        const int    tileEndClamped = min(tileEnd, dataW - 1);
        const int apronStartClamped = max(apronStart, 0);
        const int   apronEndClamped = min(apronEnd, dataW - 1);

        //Row start index in d_Data[]
        const int          rowStart = IMUL(blockIdx.y, dataW);

        //Aligned apron start. Assuming dataW and ROW_TILE_W are multiples 
        //of half-warp size, rowStart + apronStartAligned is also a 
        //multiple of half-warp size, thus having proper alignment 
        //for coalesced d_Data[] read.
        const int apronStartAligned = tileStart - KERNEL_RADIUS_ALIGNED;

        const int loadPos = apronStartAligned + threadIdx.x;
        //Set the entire data cache contents
        //Load global memory values, if indices are within the image borders,
        //or initialize with zeroes otherwise
        if (loadPos >= apronStart) {
            const int smemPos = loadPos - apronStart;

            data[smemPos] =
                ((loadPos >= apronStartClamped) && (loadPos <= apronEndClamped)) ?
                d_Data[rowStart + loadPos] : 0;
        }


        //Ensure the completness of the loading stage
        //because results, emitted by each thread depend on the data,
        //loaded by another threads
        __syncthreads();
        const int writePos = tileStart + threadIdx.x;
        //Assuming dataW and ROW_TILE_W are multiples of half-warp size,
        //rowStart + tileStart is also a multiple of half-warp size,
        //thus having proper alignment for coalesced d_Result[] write.
        if (writePos <= tileEndClamped) {
            const int smemPos = writePos - apronStart;
            real sum = 0;
            sum = data[smemPos - 1] * c_cdiag + data[smemPos] * c_cdir + data[smemPos + 1] * c_cdiag;
            d_Result[rowStart + writePos] = sum;
        }
    }

__global__ void cudaColumnComputeKernel(
    int dataW,
    int dataH,
    real* d_Src,
    real* d_Result,
    real* d_Data,
    int smemStride,
    int gmemStride,
    real* cond_self,
    real* cond_neighbors,
    real* maxdiff
) {
    //Data cache
    __shared__ real data[COLUMN_TILE_W * (KERNEL_RADIUS + COLUMN_TILE_H + KERNEL_RADIUS)];

    //Current tile and apron limits, in rows
    const int         tileStart = IMUL(blockIdx.y, COLUMN_TILE_H);
    const int           tileEnd = tileStart + COLUMN_TILE_H - 1;
    const int        apronStart = tileStart - KERNEL_RADIUS;
    const int          apronEnd = tileEnd + KERNEL_RADIUS;

    //Clamp tile and apron limits by image borders
    const int    tileEndClamped = min(tileEnd, dataH - 1);
    const int apronStartClamped = max(apronStart, 0);
    const int   apronEndClamped = min(apronEnd, dataH - 1);

    //Current column index
    const int       columnStart = IMUL(blockIdx.x, COLUMN_TILE_W) + threadIdx.x;

    //Shared and global memory indices for current column
    int smemPos = IMUL(threadIdx.y, COLUMN_TILE_W) + threadIdx.x;
    int gmemPos = IMUL(apronStart + threadIdx.y, dataW) + columnStart;
    //Cycle through the entire data cache
    //Load global memory values, if indices are within the image borders,
    //or initialize with zero otherwise
    for (int y = apronStart + threadIdx.y; y <= apronEnd; y += blockDim.y) {
        data[smemPos] =
            ((y >= apronStartClamped) && (y <= apronEndClamped)) ?
            d_Data[gmemPos] : 0;
        smemPos += smemStride;
        gmemPos += gmemStride;
    }

    //Ensure the completness of the loading stage
    //because results, emitted by each thread depend on the data, 
    //loaded by another threads
    __syncthreads();
    //Shared and global memory indices for current column
    smemPos = IMUL(threadIdx.y + KERNEL_RADIUS, COLUMN_TILE_W) + threadIdx.x;
    gmemPos = IMUL(tileStart + threadIdx.y, dataW) + columnStart;
    //Cycle through the tile body, clamped by image borders
    //Calculate and output the results
    int start = max(tileStart + threadIdx.y, 1);
    int end = min(tileEndClamped, dataH - 2);
    real res;
    real source;
    for (int y = start; y <= end; y += blockDim.y) {
        real sum = 0;
        sum = data[smemPos - 1] * 0.25 + data[smemPos] * (0.25 * c_cdir / c_cdiag) + data[smemPos + 1] * 0.25;
        source = d_Src[gmemPos];
        res = sum * cond_neighbors[gmemPos] + source * cond_self[gmemPos];
        d_Result[gmemPos] = res;
        maxdiff[gmemPos] = abs(source - res);
        smemPos += smemStride;
        gmemPos += gmemStride;
    }
}

// End of Code for separable convolution taken and adapte from NVIDIA White Paper
// =====================================================================================
// =====================================================================================
// =====================================================================================

extern "C"
void cuda_init(size_t h, size_t w){
    size_t blocksPerLine = (w-2+threadBlockSize-1)/threadBlockSize;
    size_t num_blocks = (h-2)*blocksPerLine;
    //checkCudaCall(cudaMalloc((void **) &maxdiff_block, h*w * sizeof(double)));
    checkCudaCall(cudaMalloc((void **) &maxdiff_block, num_blocks * sizeof(real)));
    checkCudaCall(cudaMalloc((void **) &min_block, num_blocks * sizeof(real)));
    checkCudaCall(cudaMalloc((void **) &max_block, num_blocks * sizeof(real)));
    checkCudaCall(cudaMalloc((void **) &sum_block, num_blocks * sizeof(real)));
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
                        real* src, real* dst, real* tmp,
                        real* cond_self, real * cond_neighbors){

    dim3 blockGridRows(iDivUp(w, ROW_TILE_W), h);
    dim3 blockGridColumns(iDivUp(w, COLUMN_TILE_W), iDivUp(h, COLUMN_TILE_H));
    dim3 threadBlockRows(KERNEL_RADIUS_ALIGNED + ROW_TILE_W + KERNEL_RADIUS);
    dim3 threadBlockColumns(COLUMN_TILE_W, 8);
    size_t blocksPerLine = (w - 2 + threadBlockSize - 1) / threadBlockSize;
    dim3 numBlocks(h - 2, blocksPerLine);
    dim3 numBlocks2(h - 2, (w + threadBlockSize - 1) / threadBlockSize);

    cudaRowComputeKernel <<<blockGridRows, threadBlockRows >>> (w, h, src, tmp);
    cudaColumnComputeKernel <<<blockGridColumns, threadBlockColumns >>>
        (w, h, 
         src, tmp, dst, 
         COLUMN_TILE_W * threadBlockColumns.y, 
         w * threadBlockColumns.y,
         cond_self, cond_neighbors, maxdiff_block);
    cudaMemcpy(&dst[0], &dst[(h - 2) * w], w * sizeof(real), cudaMemcpyDeviceToDevice);
    cudaMemcpy(&dst[(h - 1)], &dst[1 * w], w * sizeof(real), cudaMemcpyDeviceToDevice);

    int len;
    for(len=blocksPerLine*(h-2); len>=threadBlockSize*8; len=(len+1)/2){
        cudaMaxdiffKernel<<<(len+threadBlockSize-1)/threadBlockSize, threadBlockSize>>>((len+1)/2, len/2, maxdiff_block);
    }
    cudaDeviceSynchronize();

    
    // real r[blocksPerLine*(h-2)];
    // getDeviceMemory(r, maxdiff_block, blocksPerLine*(h-2)*sizeof(real));
    // real res=0;
    // for(int i=0;i<blocksPerLine*(h-2);i++){
    //     res = MAX(res, r[i]);
    // }

    real* r = (real * )malloc(sizeof(real)*len);
    getDeviceMemory(r, maxdiff_block, len*sizeof(real));
    real res=0;
    for(int i=0;i<len;i++){
        res = MAX(res, r[i]);
    }
    return res;
    
}

extern "C"
void cuda_fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, real* data,
                        struct timespec *before)
{
    // compute min/max/avg 
    real tmin = INFINITY, tmax = -INFINITY;
    real sum = 0.0;
    struct timespec after;

    size_t blocksPerLine = (w-2+threadBlockSize-1)/threadBlockSize;
    dim3 numBlocks(h-2,blocksPerLine);

    //We have said that the final reduction does not need to be included.
    clock_gettime(CLOCK_MONOTONIC, &after);
    
    cudaReduceKernel<<<numBlocks, threadBlockSize>>>(h, w, data, min_block, max_block, sum_block);
    cudaDeviceSynchronize();

    real * recv = (real * )malloc(sizeof(real)*blocksPerLine*(h-2));
    getDeviceMemory(recv, min_block, blocksPerLine*(h-2)*sizeof(real));
    for(int i=0;i<blocksPerLine*(h-2);i++){
        tmin = MIN(tmin, recv[i]);
    }

    getDeviceMemory(recv, max_block, blocksPerLine*(h-2)*sizeof(real));
    for(int i=0;i<blocksPerLine*(h-2);i++){
        tmax = MAX(tmax, recv[i]);
    }

    getDeviceMemory(recv, sum_block, blocksPerLine*(h-2)*sizeof(real));
    for(int i=0;i<blocksPerLine*(h-2);i++){
        sum += recv[i];
    }
    r->tmin = tmin;
    r->tmax = tmax;
    r->tavg = sum / (p->N * p->M);
    
    r->time = (real)(after.tv_sec - before->tv_sec) + 
        (real)(after.tv_nsec - before->tv_nsec) / 1e9;
    
}

