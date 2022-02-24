#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#define max(a,b) (a)>=(b)?(a):(b)
#define min(a,b) (a)<=(b)?(a):(b)

/* Ordering of the vector */
typedef enum Ordering {ASCENDING, DESCENDING, RANDOM} Order;
int debug = 0;
int num_threads = 0;
int inner_threads = 0;
int outer_threads = 0;

int inline binary_search(int *a, int len, int r){
    int L=0,R=len-1;
    int ans=0;
    while(L<=R){
        int mid = (L+R)/2;
        if(a[mid]<=r){
            ans=mid+1;
            L=mid+1;
        }else{
            R=mid-1;
        }
    }
    return ans;
}

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

int threshold = 16384;
void merge_parallel(int *a, int *b, int *c, int num_lhs, int num_rhs){
    int l = num_lhs + num_rhs;
    if(l<=threshold){
        merge(a,b,c,num_lhs,num_rhs,l);
        return;
    }
    if(num_lhs<num_rhs){
        int* temp=b;
        b=c;
        c=temp;
        int t=num_lhs;
        num_lhs=num_rhs;
        num_rhs=t;
    }
    int rb = num_lhs/2;
    int rc = binary_search(c,num_rhs,b[rb-1]);
    //printf("%d %d\n",rb,rc);
    a[rb+rc-1] = b[rb-1];
    #pragma omp task
    {
        merge_parallel(a,b,c,rb-1,rc);
    }
    #pragma omp task
    {
        merge_parallel(a+rb+rc,b+rb,c+rc,num_lhs-rb,num_rhs-rc);
    }
    #pragma omp taskwait


}

void top_down_mergesort(int *b, int l, int *a){
    if(l<=1)
        return;
    int num_rhs = l/2;
    int num_lhs = l-num_rhs;
    int *mid_a = a+num_lhs;
    int *mid_b = b+num_lhs;
    top_down_mergesort(a, num_lhs, b);
    top_down_mergesort(mid_a, num_rhs, mid_b);
    merge(a,b,mid_b,num_lhs,num_rhs,l);
}

void top_down_mergesort_parallel(int *b, int l, int *a, bool flag){
    if(l<=threshold) {
        if(flag) memcpy(a,b,l*sizeof(int));
        else memcpy(b,a,l*sizeof(int));
        top_down_mergesort(b,l,a);
        return;    
    }else{
        int num_rhs = l/2;
        int num_lhs = l-num_rhs;
        int *mid_a = a+num_lhs;
        int *mid_b = b+num_lhs;
        #pragma omp task
        {
            top_down_mergesort_parallel(a, num_lhs, b, !flag);
        }
        #pragma omp task
        {
            top_down_mergesort_parallel(mid_a, num_rhs, mid_b, !flag);
        }
        #pragma omp taskwait
        merge_parallel(a, b, mid_b, num_lhs, num_rhs);
        //merge(a, b, mid_b, num_lhs, num_rhs, l);

    }

}

/* Sort vector v of l elements using mergesort */
void vecsort(int **vector_vectors, int *vector_lengths, long length_outer, int length_inner_max){

    int *b;

    omp_set_nested(1);

#pragma omp parallel num_threads(num_threads) private(b) 
    {
        b = (int*)malloc(sizeof(int)*length_inner_max);
#pragma omp for schedule(guided)
        for(long i = 0; i < length_outer; i++) {
            
            memcpy(b,vector_vectors[i],vector_lengths[i]*sizeof(int));
// #pragma omp parallel num_threads(inner_threads) 
//             {
// #pragma omp single 
//                 {
//                     top_down_mergesort_parallel(b, vector_lengths[i], vector_vectors[i],0);
            top_down_mergesort(b, vector_lengths[i], vector_vectors[i]);
//                 };
//             }
//         }
        free(b);
    }
}

void print_v(int **vector_vectors, int *vector_lengths, long length_outer) {
    printf("\n");
    for(long i = 0; i < length_outer; i++) {
        for (int j = 0; j < vector_lengths[i]; j++){
            if(j != 0 && (j % 10 == 0)) {
                printf("\n");
            }
            printf("%d ", vector_vectors[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {

    int c;
    int seed = 42;
    long length_outer = 1e4;
    Order order = ASCENDING;
    int length_inner_min = 100;
    int length_inner_max = 1000;

    int **vector_vectors;
    int *vector_lengths;

    struct timespec before, after;


    /* Read command-line options. */
    while ((c = getopt(argc, argv, "adrgn:x:l:p:s:t:o:i:")) != -1) {
        switch (c) {
            case 't': //threshold
                threshold = atol(optarg);
                break;
            case 'o': //outer threads
                outer_threads = atoi(optarg);
                break;
            case 'i': // inner threads
                inner_threads = atoi(optarg);
                break;
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
                length_outer = atol(optarg);
                break;
            case 'n':
                length_inner_min = atoi(optarg);
                break;
            case 'x':
                length_inner_max = atoi(optarg);
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
                if (optopt == 'l' || optopt == 's') {
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                } else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option '-%c'.\n", optopt);
                } else {
                    fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                }
                return -1;
            default:
                return -1;
        }
    }

    // set number of threads
    if(num_threads==0){
        if(outer_threads && inner_threads){
            num_threads = outer_threads*inner_threads; 
        }
        else if(outer_threads){
            inner_threads = 32/outer_threads;
            num_threads = outer_threads*inner_threads; 
        }else{
            inner_threads = 2;
            outer_threads = 16;
            num_threads = 32;
        }
    }else{
        if(outer_threads){
            inner_threads = num_threads / outer_threads;
        }else if(inner_threads){
            outer_threads = num_threads / outer_threads;
        }else{
            if(num_threads%2==0) inner_threads = 2;
            else inner_threads = 1;
            outer_threads = num_threads/inner_threads;
        }
    }
    omp_set_num_threads(num_threads);

    /* Seed such that we can always reproduce the same random vector */
    srand(seed);

    /* Allocate vector. */
    vector_vectors = (int **) malloc(length_outer * sizeof(int*));
    vector_lengths = (int *) malloc(length_outer * sizeof(int));
    if (vector_vectors == NULL || vector_lengths == NULL) {
        fprintf(stderr, "Malloc failed...\n");
        return -1;
    }

    assert(length_inner_min < length_inner_max);

    /* Determine length of inner vectors and fill them. */
    for (long i = 0; i < length_outer; i++) {
        int length_inner = (rand() % (length_inner_max + 1 - length_inner_min)) + length_inner_min ; //random number inclusive between min and max
        vector_vectors[i] = (int *) malloc(length_inner * sizeof(int));
        vector_lengths[i] = length_inner;

        /* Allocate and fill inner vector. */
        switch (order) {
            case ASCENDING:
                for (long j = 0; j < length_inner; j++) {
                    vector_vectors[i][j] = (int) j;
                }
                break;
            case DESCENDING:
                for (long j = 0; j < length_inner; j++) {
                    vector_vectors[i][j] = (int) (length_inner - j);
                }
                break;
            case RANDOM:
                for (long j = 0; j < length_inner; j++) {
                    vector_vectors[i][j] = rand();
                }
                break;
        }
    }

    if(debug) {
        print_v(vector_vectors, vector_lengths, length_outer);
    }

    clock_gettime(CLOCK_MONOTONIC, &before);

    /* Sort */
    vecsort(vector_vectors, vector_lengths, length_outer, length_inner_max);

    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
              (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Vecsort took: % .6e seconds \n", time);

    if(debug) {
        print_v(vector_vectors, vector_lengths, length_outer);
    }

    return 0;
}

/*
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=active
export OMP_NUM_THREADS=32
export OMP_PLACES=cores
*/