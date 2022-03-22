#ifndef CUDA_COMPUTE_H
#define CUDA_COMPUTE_H
#include <time.h>
#include "cuda_compute.h"
#include "compute.h"

#ifdef __cplusplus
extern "C" 
#endif 

void cuda_do_compute();
void* createDeviceMemory(void* p, size_t length);
void getDeviceMemory(void* host, void* device, size_t length);
void freeDeviceMemory(void* p);
void cuda_init(size_t h, size_t w);
void cuda_finalize();

double cuda_do_compute_step(size_t h, size_t w, double* src, double* dst, double* c);
void cuda_fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, double* data,
                        struct timespec *before);
#endif
