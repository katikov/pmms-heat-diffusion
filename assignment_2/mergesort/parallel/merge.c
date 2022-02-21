#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <omp.h>
#define max(a,b) (a)>=(b)?(a):(b)
#define min(a,b) (a)<=(b)?(a):(b)

/* Ordering of the vector */
typedef enum Ordering {ASCENDING, DESCENDING, RANDOM} Order;

int debug = 0;
int num_threads = 1;

void inline merge(int *a, int *b, int *mid_b, int num_lhs, int num_rhs, int l){
    int i=0,j=0;
    for(int k = 0; k < l; k++) {
        if (i < num_lhs && (j >= num_rhs || b[i] <= mid_b[j])) {
            a[k] = b[i];
            i++;
        } else {
            a[k] = mid_b[j];
            j++;
        }
    }

}
void top_down_mergesort(int *b, long l, int *a){
    if(l<=1)
        return;
    int num_rhs = l/2;
    int num_lhs = l-num_rhs;
    int *mid_a = a+num_lhs;
    int *mid_b = b+num_lhs;
    top_down_mergesort(a, num_lhs, b);
    top_down_mergesort(mid_a, num_rhs, mid_b);
    merge(a, b, mid_b, num_lhs, num_rhs, l);
}

const int static threshold = 1024;
void top_down_mergesort_parallel(int *b, long l, int *a){
    if(l<=threshold) {
        // if(flag) memcpy(a,b,l*sizeof(int));
        // else memcpy(b,a,l*sizeof(int));
        top_down_mergesort(b,l,a);
        return;    
    }else{
        int num_rhs = l/2;
        int num_lhs = l-num_rhs;
        int *mid_a = a+num_lhs;
        int *mid_b = b+num_lhs;
        #pragma omp task
        {
            top_down_mergesort_parallel(a, num_lhs, b);
        }
        #pragma omp task
        {
            top_down_mergesort_parallel(mid_a, num_rhs, mid_b);
        }
        #pragma omp taskwait
        merge(a, b, mid_b, num_lhs, num_rhs, l);

    }

}

/* Sort vector v of l elements using mergesort */
void msort(int *v, long l, int *b){

    //memcpy(b,v,l*sizeof(int));
    int copy_threads = min(4,num_threads);
    #pragma omp parallel for num_threads(copy_threads) if(num_threads>1)
    for(int i=0;i<l;i++){
        b[i]=v[i];
    }
    #pragma omp parallel
    {
        #pragma omp single
        {
            //if(debug)
            top_down_mergesort_parallel(b,l,v);
        }
        
    }
    
}


void print_v(int *v, long l) {
    printf("\n");
    for(long i = 0; i < l; i++) {
        if(i != 0 && (i % 10 == 0)) {
            printf("\n");
        }
        printf("%d ", v[i]);
    }
    printf("\n");
}

int main(int argc, char **argv) {

    int c;
    int seed = 42;
    long length = 1e4;
    
    Order order = ASCENDING;
    int *vector;

    struct timespec before, after;

    /* Read command-line options. */
    while((c = getopt(argc, argv, "adrgp:l:s:")) != -1) {
        switch(c) {
            case 'a':
                order = ASCENDING;
                break;
            case 'd':
                order = DESCENDING;
                break;
            case 'r':
                order = RANDOM;
                break;
            case 'l':
                length = atol(optarg);
                break;
            case 'g':
                debug = 1;
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'p':
                num_threads = atoi(optarg);
                break;
            case '?':
                if(optopt == 'l' || optopt == 's') {
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                }
                else if(isprint(optopt)) {
                    fprintf(stderr, "Unknown option '-%c'.\n", optopt);
                }
                else {
                    fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                }
                return -1;
            default:
                return -1;
        }
    }

    /* Seed such that we can always reproduce the same random vector */
    srand(seed);

    /* Allocate vector. */
    vector = (int*)malloc(length*sizeof(int));
    if(vector == NULL) {
        fprintf(stderr, "Malloc failed...\n");
        return -1;
    }

    /* Fill vector. */
    switch(order){
        case ASCENDING:
            for(long i = 0; i < length; i++) {
                vector[i] = (int)i;
            }
            break;
        case DESCENDING:
            for(long i = 0; i < length; i++) {
                vector[i] = (int)(length - i);
            }
            break;
        case RANDOM:
            for(long i = 0; i < length; i++) {
                vector[i] = rand();
            }
            break;
    }

    if(debug) {
        print_v(vector, length);
    }

    int *b = (int*)malloc(length*sizeof(int));
    if(b == NULL) {
        fprintf(stderr, "Malloc failed...\n");
        exit(-1);
    }

    clock_gettime(CLOCK_MONOTONIC, &before);

    /* Sort */
    msort(vector, length,b);

    clock_gettime(CLOCK_MONOTONIC, &after);
    free(b);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Mergesort took: % .6e seconds \n", time);

    if(debug) {
        print_v(vector, length);
    }

    return 0;
}

/*
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=active
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
*/