#include <stdio.h>
#include <string.h>
#include "timer.h"

#define image_height 1024
#define image_width 1024
#define filter_height 5
#define filter_width 5

#define border_height ((filter_height/2)*2)
#define border_width ((filter_width/2)*2)
#define input_height (image_height + border_height)
#define input_width (image_width + border_width)

#define block_size_x 32
#define block_size_y 16

#define SEED 1234

using namespace std;
__constant__ float _filter[filter_height*filter_width];
void convolutionSeq(float *output, float *input, float *filter) {
    //for each pixel in the output image

  timer sequentialTime = timer("Sequential");
  
  sequentialTime.start();

    for (int y=0; y < image_height; y++) {
        for (int x=0; x < image_width; x++) { 
	    output[y*image_width+x]=0;
            //for each filter weight
            for (int i=0; i < filter_height; i++) {
                for (int j=0; j < filter_width; j++) {
                    output[y*image_width+x] += input[(y+i)*input_width+x+j] * filter[i*filter_width+j];
                }
            }
	    output[y*image_width+x] /= 35;
        }
    }
  
  sequentialTime.stop(); 
  cout << "convolution (sequential): \t\t" << sequentialTime << endl;

}


__global__ void convolution_kernel_naive(float *output, float *input, float *filter) {
    unsigned y = blockIdx.y*blockDim.y + threadIdx.y;
    unsigned x = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned idx = y*input_width+x;
    if(x<image_width && y<image_height){
        float res = 0.0;
        for(int i=0;i<filter_height;i++){
            for(int j=0;j<filter_width;j++){
                res += input[idx+i*input_width+j] * filter[i*filter_width+j];
            }
        }
        output[y*image_width+x] = res/35.0;
    }
}

__global__ void convolution_kernel_constant(float *output, float *input) {
    unsigned y = blockIdx.y*blockDim.y + threadIdx.y;
    unsigned x = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned idx = y*input_width+x;
    if(x<image_width && y<image_height){
        float res = 0.0;
        for(int i=0;i<filter_height;i++){
            for(int j=0;j<filter_width;j++){
                res += input[idx+i*input_width+j] * _filter[i*filter_width+j];
            }
        }
        output[y*image_width+x] = res/35.0;
    }
}

__global__ void convolution_kernel(float *output, float *input, int start_y) {
    unsigned y = blockIdx.y*blockDim.y + threadIdx.y + start_y;
    unsigned x = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned idx = y*input_width+x;
    if(x<image_width && y<image_height){
        float res = 0.0;
        for(int i=0;i<filter_height;i++){
            for(int j=0;j<filter_width;j++){
                res += input[idx+i*input_width+j] * _filter[i*filter_width+j];
            }
        }
        output[y*image_width+x] = res/35.0;
    }
}

void convolutionCUDA(float *output, float *input, float *filter) {
    float *d_input; float *d_output; // float *d_filter;
    cudaError_t err;
    timer kernelTime = timer("kernelTime");
    timer memoryTime = timer("memoryTime");

    // memory allocation
    err = cudaMalloc((void **)&d_input, input_height*input_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaMalloc((void **)&d_output, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString( err )); }

kernelTime.start();
    memoryTime.start();

    err = cudaMemset(d_output, 0, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemset output: %s\n", cudaGetErrorString( err ));  }
    cudaMemcpyToSymbol(_filter, filter, sizeof(float)*filter_height*filter_width);
    memoryTime.stop();

    dim3 threads(block_size_x, block_size_y);
    const int nstreams = 2;
    int grid_x = (image_width+threads.x-1)/threads.x;
    int grid_y = ((image_height+threads.y-1)/threads.y+nstreams-1)/nstreams;
    dim3 grid(grid_x, grid_y);
    cudaStream_t streams[nstreams];

    for(int i=0;i<nstreams;i++){
        cudaStreamCreate(&streams[i]);
        int start_pos = i*threads.y*grid_y;
        int length = threads.y*grid_y + border_height;
        if(start_pos+length > input_height) length = input_height - start_pos;
        cudaMemcpyAsync(&d_input[start_pos*input_width], &input[start_pos*input_width], length * input_width * sizeof(float), cudaMemcpyHostToDevice, streams[i]);
        convolution_kernel<<<grid, threads, 0, streams[i]>>>(d_output, d_input, start_pos);
    }
    cudaDeviceSynchronize();
    memoryTime.start();
    cudaMemcpy(output, d_output, image_height*image_width*sizeof(float), cudaMemcpyDeviceToHost);
    memoryTime.stop();
kernelTime.stop();


    err = cudaFree(d_input);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaFree(d_output);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString( err )); }

    cout << "convolution (cuda): \t\t" << kernelTime << endl;
    //cout << "convolution (memory): \t\t" << memoryTime << endl;
}


void convolutionCUDAConstant(float *output, float *input, float *filter) {
    float *d_input; float *d_output; // float *d_filter;
    cudaError_t err;
    timer kernelTime = timer("kernelTime");
    timer memoryTime = timer("memoryTime");
    timer totalTime = timer("totalTime");

    // memory allocation
    err = cudaMalloc((void **)&d_input, input_height*input_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaMalloc((void **)&d_output, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString( err )); }
    //err = cudaMalloc((void **)&d_filter, filter_height*filter_width*sizeof(float));
    //if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_filter: %s\n", cudaGetErrorString( err )); }

    totalTime.start();
    memoryTime.start();
    cudaMemcpyToSymbol(_filter, filter, sizeof(float)*filter_height*filter_width);
    // host to device 
    err = cudaMemcpy(d_input, input, input_height*input_width*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device input: %s\n", cudaGetErrorString( err ));  }
    //err = cudaMemcpy(d_filter, filter, filter_height*filter_width*sizeof(float), cudaMemcpyHostToDevice);
    //if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device filter: %s\n", cudaGetErrorString( err ));  }
    
    // zero the result array 
    err = cudaMemset(d_output, 0, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemset output: %s\n", cudaGetErrorString( err ));  }
    memoryTime.stop();
    //setup the grid and thread blocks
    //thread block size
    dim3 threads(block_size_x, block_size_y);
    //problem size divided by thread block size rounded up
    dim3 grid(int(ceilf(image_width/(float)threads.x)), int(ceilf(image_height/(float)threads.y)) );

    //measure the GPU function
    kernelTime.start();
    convolution_kernel_constant<<<grid, threads>>>(d_output, d_input);
    //dim3 test((image_width+255)/256, image_height);
    //convolution_kernel_naive2<<<test, 256>>>(d_output, d_input, d_filter);
    cudaDeviceSynchronize();
    kernelTime.stop();
 
    //check to see if all went well
    err = cudaGetLastError();
    if (err != cudaSuccess) { fprintf(stderr, "Error during kernel launch convolution_kernel: %s\n", cudaGetErrorString( err )); }

    //copy the result back to host memory
    memoryTime.start();
    err = cudaMemcpy(output, d_output, image_height*image_width*sizeof(float), cudaMemcpyDeviceToHost);
    memoryTime.stop();
    totalTime.stop();
    
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy device to host output: %s\n", cudaGetErrorString( err )); }
 
    err = cudaFree(d_input);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaFree(d_output);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString( err )); }
    //err = cudaFree(d_filter);
    //if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_filter: %s\n", cudaGetErrorString( err )); }

    cout << "convolution (kernel): \t\t" << kernelTime << endl;
    cout << "convolution (memory): \t\t" << memoryTime << endl;
    cout << "convolution (total): \t\t" << totalTime << endl;

}

void convolutionCUDANaive(float *output, float *input, float *filter) {
    float *d_input; float *d_output; float *d_filter;
    cudaError_t err;
    timer kernelTime = timer("kernelTime");
    timer memoryTime = timer("memoryTime");
    timer totalTime = timer("totalTime");

    // memory allocation
    err = cudaMalloc((void **)&d_input, input_height*input_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaMalloc((void **)&d_output, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString( err )); }
    err = cudaMalloc((void **)&d_filter, filter_height*filter_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_filter: %s\n", cudaGetErrorString( err )); }

    totalTime.start();
    memoryTime.start();
    //cudaMemcpyToSymbol(_filter, filter, sizeof(float)*filter_height*filter_width);
    // host to device 
    err = cudaMemcpy(d_input, input, input_height*input_width*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device input: %s\n", cudaGetErrorString( err ));  }
    err = cudaMemcpy(d_filter, filter, filter_height*filter_width*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device filter: %s\n", cudaGetErrorString( err ));  }
    
    // zero the result array 
    err = cudaMemset(d_output, 0, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemset output: %s\n", cudaGetErrorString( err ));  }
    memoryTime.stop();
    //setup the grid and thread blocks
    //thread block size
    dim3 threads(block_size_x, block_size_y);
    //problem size divided by thread block size rounded up
    dim3 grid(int(ceilf(image_width/(float)threads.x)), int(ceilf(image_height/(float)threads.y)) );

    //measure the GPU function
    kernelTime.start();
    convolution_kernel_naive<<<grid, threads>>>(d_output, d_input, d_filter);
    cudaDeviceSynchronize();
    kernelTime.stop();
 
    //check to see if all went well
    err = cudaGetLastError();
    if (err != cudaSuccess) { fprintf(stderr, "Error during kernel launch convolution_kernel: %s\n", cudaGetErrorString( err )); }

    //copy the result back to host memory
    memoryTime.start();
    err = cudaMemcpy(output, d_output, image_height*image_width*sizeof(float), cudaMemcpyDeviceToHost);
    memoryTime.stop();
    totalTime.stop();

    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy device to host output: %s\n", cudaGetErrorString( err )); }
 
    err = cudaFree(d_input);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaFree(d_output);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString( err )); }
    err = cudaFree(d_filter);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_filter: %s\n", cudaGetErrorString( err )); }

    cout << "convolution (kernel): \t\t" << kernelTime << endl;
    cout << "convolution (memory): \t\t" << memoryTime << endl;
    cout << "convolution (total): \t\t" << totalTime << endl;

}

int compare_arrays(float *a1, float *a2, int n) {
    int errors = 0;
    int print = 0;

    for (int i=0; i<n; i++) {

        if (isnan(a1[i]) || isnan(a2[i])) {
            errors++;
            if (print < 10) {
                print++;
                fprintf(stderr, "Error NaN detected at i=%d,\t a1= %10.7e \t a2= \t %10.7e\n",i,a1[i],a2[i]);
            }
        }

        float diff = (a1[i]-a2[i])/a1[i];
        if (diff > 1e-6f) {
            errors++;
            if (print < 10) {
                print++;
                fprintf(stderr, "Error detected at i=%d, \t a1= \t %10.7e \t a2= \t %10.7e \t rel_error=\t %10.7e\n",i,a1[i],a2[i],diff);
            }
        }

    }

    return errors;
}
        

int main() {
    int i; 
    int errors=0;

    //allocate arrays and fill them
    float *input = (float *) malloc(input_height * input_width * sizeof(float));
    float *output1 = (float *) calloc(image_height * image_width, sizeof(float));
    float *output2 = (float *) calloc(image_height * image_width, sizeof(float));
    float *filter = (float *) malloc(filter_height * filter_width * sizeof(float));

    for (i=0; i< input_height * input_width; i++) {
        input[i] = (float) (i % SEED);
    }

//THis is specific for a W==H smoothening filteri, where W and H are odd.
    for (i=0; i<filter_height * filter_width; i++) { 
      filter[i] = 1.0;
    }

    for (i=filter_width+1; i<(filter_height - 1) * filter_width; i++) {
	if (i % filter_width > 0 && i % filter_width < filter_width-1) filter[i]+=1.0; 
    }

    filter[filter_width*filter_height/2]=3.0;
//end initialization
   
    //measure the CPU function
    convolutionSeq(output1, input, filter);
    //measure the GPU function
    convolutionCUDA(output2, input, filter);
    //convolutionCUDAConstant(output2, input, filter);
    //convolutionCUDANaive(output2, input, filter);

    //check the result
    errors += compare_arrays(output1, output2, image_height*image_width);
    if (errors > 0) {
        printf("TEST FAILED! %d errors!\n", errors);
    } else {
        printf("TEST PASSED!\n");
    }


    free(filter);
    free(input);
    free(output1);
    free(output2);

    return 0;
}


