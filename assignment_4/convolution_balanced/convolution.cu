#include <stdio.h>
#include <string.h>
#include "timer.h"

//#include "config.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>


#define block_size_x 32
#define block_size_y 16

#define SEED 1234

using namespace std;

struct image_params {
    int input_h;
    int input_w;
    int image_h;
    int image_w;
    int cpu_h;
    int cpu_w;
    int gpu_h;
    int gpu_w;
    int filter_h;
    int filter_w;
    int ratio;
    float* filter;
    float* input;
    float* output;
};

void convolutionSeq(image_params* params) {
    //for each pixel in the output image

    timer sequentialTime = timer("Sequential");
    int input_height = params->input_h;
    int input_width = params->input_w;
    int image_height = params->cpu_h;
    int image_width = params->cpu_w;
    int filter_height = params->filter_h;
    int filter_width = params->filter_w;
    float* input = params->input;
    float* output = params->output;
    float* filter = params->filter;

    sequentialTime.start();

    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {
            output[y * image_width + x] = 0;
            //for each filter weight
            for (int i = 0; i < filter_height; i++) {
                for (int j = 0; j < filter_width; j++) {
                    params -> output[y * image_width + x] += input[(y + i) * input_width + x + j] * filter[i * filter_width + j];
                }
            }
            params -> output[y * image_width + x] /= 35;
        }
    }
    sequentialTime.stop();
    cout << sequentialTime << ",";
}

__constant__ int d_input_width;
__constant__ int d_input_height;
__constant__ int d_filter_width;
__constant__ int d_filter_height;
__constant__ int d_image_width;
__constant__ int d_image_height;

__global__ void convolution_kernel(float* output, float* input, float *d_filter, int start_y) {
    extern __shared__ int s_filter[];
    int linid = threadIdx.x + threadIdx.y * blockDim.x;
    if (linid < d_filter_height * d_filter_width) {
        s_filter[linid] = d_filter[linid];
    }

    unsigned y = blockIdx.y * blockDim.y + threadIdx.y + start_y;
    unsigned x = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned idx = y * d_input_width + x;
    if (x < d_image_width && y < d_image_height) {
        float res = 0.0;
        for (int i = 0; i < d_filter_height; i++) {
            for (int j = 0; j < d_filter_width; j++) {
                res += input[idx + i * d_input_width + j] * s_filter[i * d_filter_width + j];
            }
        }
        output[y * d_image_width + x] = res / 35.0;
    }
}

void convolutionCUDA(image_params *params) {
    int image_height = params->gpu_h;
    int image_width = params->gpu_w;
    int filter_height = params->filter_h;
    int filter_width = params->filter_w;
    int input_height = image_height + filter_height;
    int input_width = params->input_w;
    float* input = params->input + params->cpu_h*params->input_w;
    float * output = params->output + params->cpu_h * params->image_w;
    float * filter = params->filter;

    float* d_input; float* d_output; float *d_filter;
    cudaError_t err;
    timer kernelTime = timer("kernelTime");
    timer memoryTime = timer("memoryTime");

    // memory allocation
    err = cudaMalloc((void**)&d_input, input_height * input_width * sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString(err)); }
    err = cudaMalloc((void**)&d_output, image_height * image_width * sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString(err)); }
    err = cudaMalloc((void**)&d_filter, filter_height * filter_width * sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString(err)); }

    kernelTime.start();
    memoryTime.start();

    err = cudaMemset(d_output, 0, image_height * image_width * sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemset output: %s\n", cudaGetErrorString(err)); }
    cudaMemcpyAsync(d_filter, filter, sizeof(float) * filter_height * filter_width, cudaMemcpyHostToDevice);
    memoryTime.stop();

    dim3 threads(block_size_x, block_size_y);
    const int nstreams = 2;
    int grid_x = (image_width + threads.x - 1) / threads.x;
    int grid_y = ((image_height + threads.y - 1) / threads.y + nstreams - 1) / nstreams;
    dim3 grid(grid_x, grid_y);
    cudaStream_t streams[nstreams];
    
    cudaMemcpyToSymbol("d_input_width", &input_width, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("d_input_height", &input_height, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("d_filter_width", &filter_width, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("d_filter_height", &filter_height, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("d_image_width", &image_width, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("d_image_height", &image_height, sizeof(int), 0, cudaMemcpyHostToDevice);

    int border_height = (filter_height / 2) * 2;
    for (int i = 0; i < nstreams; i++) {
        cudaStreamCreate(&streams[i]);
        int start_pos = i * threads.y * grid_y;
        int length = threads.y * grid_y + border_height;
        if (start_pos + length > input_height) length = input_height - start_pos;
        cudaMemcpyAsync(&d_input[start_pos * input_width], &input[start_pos * input_width], length * input_width * sizeof(float), cudaMemcpyHostToDevice, streams[i]);
        convolution_kernel <<<grid, threads, filter_height * filter_width, streams[i] >>> (d_output, d_input, d_filter, start_pos);
    }

    cudaDeviceSynchronize();
    memoryTime.start();
    cudaMemcpy(output, d_output, image_height * image_width * sizeof(float), cudaMemcpyDeviceToHost);
    memoryTime.stop();
    kernelTime.stop();


    if (cudaFree(d_input) != cudaSuccess) { fprintf(stderr, "Error in freeing d_input: %s\n", cudaGetErrorString(err)); }
    if (cudaFree(d_output) != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString(err)); }

    cout << memoryTime << ",";
    cout << kernelTime << ",";
}

static void usage(const char* pname)
{
    printf("Usage: %s [OPTION]...\n"
        "\n"
        "  -n NUM     Height.\n"
        "  -m NUM     Width.\n"
        "  -p NUM     GPU Workload Ratio %.\n"
        "  -k NUM     Filter Size.\n"
        , pname);
    exit(0);
}

void read_parameters_f(struct image_params* p, int argc, char** argv)
{
    /* set defaults */
    p->image_h = 1024;
    p->image_w = 1024;
    p->filter_h = 5;
    p->ratio = 0;   
    int ch;

    while ((ch = getopt(argc, argv, "hH:k:m:M:n:N:p:")) != -1)
    {
        switch (ch) {
        case 'm': case 'M': p->image_w = strtol(optarg, 0, 10); break;
        case 'n': case 'N': p->image_h = strtol(optarg, 0, 10); break;
        case 'k': p->filter_h = strtol(optarg, 0, 10); break;
        case 'p': p->ratio = strtol(optarg, 0, 10); break;
        case 'h': default: usage(argv[0]);
        }
    }
    p->ratio = p->ratio < 0 ? 0 : p->ratio;
    p->ratio = p->ratio > 100 ? 100 : p->ratio;

    /*printf("Parameters:\n"
        "  -n %d # number of rows\n"
        "  -m %d # number of columns\n"
        "  -p %d # GPU work share\n"
        "  -k %d # kernel size\n",
        p->image_h, p->image_w, p->ratio, p->filter_h);
    */
}

int main(int argc, char** argv)
{
    // Input Handling
    struct image_params ip;

    read_parameters_f(&ip,argc,argv);

    size_t image_height = ip.image_h;
    size_t image_width = ip.image_w;
    size_t filter_width = ip.filter_h;
    size_t filter_height = ip.filter_h;
    uint32_t ratio = ip.ratio;


    size_t border_height = ((filter_height / 2) * 2);
    size_t border_width = ((filter_width / 2) * 2);
    size_t input_height = (image_height + border_height);
    size_t input_width = (image_width + border_width);

    int i;

    //allocate arrays and fill them
    float* input = (float*)malloc(input_height * input_width * sizeof(float));
    float* output = (float*)calloc(image_height * image_width, sizeof(float));
    float* filter = (float*)malloc(filter_height * filter_width * sizeof(float));

    for (i = 0; i < input_height * input_width; i++) {
        input[i] = (float)(i % SEED);
    }

    // Gaussian Filter Generation

    //This is specific for a W==H smoothening filteri, where W and H are odd.
    for (i = 0; i < filter_height * filter_width; i++) {
        filter[i] = 1.0;
    }
    for (i = filter_width + 1; i < (filter_height - 1) * filter_width; i++) {
        if (i % filter_width > 0 && i % filter_width < filter_width - 1) filter[i] += 1.0;
    }

    filter[filter_width * filter_height / 2] = 3.0;
    //end initialization

    int lower = (image_height * ratio) / 100;

    ip.ratio = ratio;
    ip.input_h = input_height;
    ip.input_w = input_width;
    ip.image_h = image_height;
    ip.image_w = image_width;
    ip.cpu_h = image_height - lower;
    ip.cpu_w = image_width;
    ip.gpu_h = lower;
    ip.gpu_w = image_width;
    ip.filter_h = filter_height;
    ip.filter_w = filter_width;
    ip.filter = filter;
    ip.input = input;
    ip.output = output;

    timer totalTime = timer("Sequential");
    // Run the mixed model
    totalTime.start();
    convolutionSeq(&ip);
    convolutionCUDA(&ip);
    totalTime.stop();
    cout << totalTime << endl;

    free(filter);
    free(input);
    free(output);

    return 0;
}





