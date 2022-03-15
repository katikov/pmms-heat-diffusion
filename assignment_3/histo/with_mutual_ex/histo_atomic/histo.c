#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>


void die(const char *msg){
    if (errno != 0) 
        perror(msg);
    else
        fprintf(stderr, "error: %s\n", msg);
    exit(1);
}   

void generate_image(int num_rows, int num_cols, int * image){
    for (int i = 0; i < num_cols * num_rows; ++i)
    {
        image[i] = rand() % 256; //255 + 1 for num bins
    }
}

void read_image(const char * image_path, int num_rows, int num_cols, int * image){
	char format[3];
    FILE *f;
    unsigned imgw, imgh, maxv, v;
    size_t i;

	printf("Reading PGM data from %s...\n", image_path);

	if (!(f = fopen(image_path, "r"))) die("fopen");

	fscanf(f, "%2s", format);
    if (format[0] != 'P' || format[1] != '2') die("only ASCII PGM input is supported");
    
    if (fscanf(f, "%u", &imgw) != 1 ||
        fscanf(f, "%u", &imgh) != 1 ||
        fscanf(f, "%u", &maxv) != 1) die("invalid input");

    if (imgw != num_cols || imgh != num_rows) {
        fprintf(stderr, "input data size (%ux%u) does not match cylinder size (%zux%zu)\n",
                imgw, imgh, num_cols, num_rows);
        die("invalid input");
    }

    for (i = 0; i < num_cols * num_rows; ++i)
    {
        if (fscanf(f, "%u", &v) != 1) die("invalid data");
        image[i] = ((int)v * 255) / maxv; //255 for num bins
    }
    fclose(f);
}

void print_histo(int * histo){
	for (int i = 0; i < 256; ++i)
	{	
		if(i != 0 && (i % 10 == 0)) {
            printf("\n");
        }
		printf("%d ", histo[i]);
	}
    printf("\n");
}

void print_image(int num_rows, int num_cols, int * image){
	int index = 0;
	for (int i = 0; i < num_rows; ++i){	
		for (int j = 0; j < num_cols; ++j){
	        index = i * num_cols + j;
			printf("%d ", image[index]);
		}
	}
	printf("\n");
}

typedef struct histo_thread_params
{
    int* data_start; // start of the sequence chunk
    int* data_end; // end of the sequence chunk
    int* histo; // min and max range element values
} histo_thread_params;

void* histogram(void* tparams){
    histo_thread_params* params = (histo_thread_params*)tparams;
    int* image_start = params->data_start;
    int* image_end = params->data_end;
    int* histo = params->histo;

    int* iter = image_start;
    while(iter < image_end) {
        //histo[*iter] += 1;
        __atomic_fetch_add(&(histo[*iter]), 1, __ATOMIC_SEQ_CST);
        ++iter;
    }
    pthread_exit(NULL);
}

int main(int argc, char *argv[]){
    int c;
    int seed = 42;
    const char *image_path = 0;
    image_path ="../../../../images/pat1_100x150.pgm";
    int gen_image = 0;
    int debug = 0;

    int num_rows = 150;
    int num_cols = 100;
    int num_threads = 1;

    struct timespec before, after;

    int * histo = (int *) calloc(256, sizeof(int));

    /* Read command-line options. */
    while((c = getopt(argc, argv, "s:i:rp:n:m:g")) != -1) {
        switch(c) {
            case 's':
                seed = atoi(optarg);
                break;
            case 'i':
            	image_path = optarg;
            	break;
            case 'r':
            	gen_image = 1;
            	break;
            case 'p':
                num_threads = atoi(optarg);
                break;
            case 'n':
            	num_rows = strtol(optarg, 0, 10);
            	break;
            case 'm':
				num_cols = strtol(optarg, 0, 10);
				break;
			case 'g':
				debug = 1;
				break;
            case '?':
                fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                return -1;
            default:
                return -1;
        }
    }

    int * image = (int *) malloc(sizeof(int) * num_cols * num_rows);
    int* counts = (int*)malloc(sizeof(int) * num_threads);
    int* image_start = image;

    /* Seed such that we can always reproduce the same random vector */
    if (gen_image){
    	srand(seed);
    	generate_image(num_rows, num_cols, image);
    }else{
    	read_image(image_path,num_rows, num_cols, image);
    }

    // Threading Boilerplate
    pthread_attr_t attr;
    pthread_attr_init ( &attr ); 
    pthread_attr_setscope ( &attr , PTHREAD_SCOPE_SYSTEM );
    histo_thread_params* args = malloc(num_threads * sizeof(*args));
    pthread_t threads[num_threads];
    int excess = (num_cols * num_rows) % num_threads;
    int base = (num_cols * num_rows) / num_threads;
    for (int i = 0; i < num_threads; ++i) {
        counts[i] = base + (i < excess?1:0);
    }

    clock_gettime(CLOCK_MONOTONIC, &before);
    /* Do your thing here */
    for (int i = 0; i < num_threads; ++i){
        args[i] = (histo_thread_params){
            .data_start = image_start,
            .data_end = image_start + counts[i],
            .histo = histo
            };
        image_start += counts[i];
        pthread_create(&threads[i], &attr, histogram, (void*)&args[i]);
        }
    for (int i = 0; i < num_threads; ++i)
        pthread_join(threads[i], NULL);

    /* Do your thing here */

    if (debug){
    	print_histo(histo);
    }

    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Histo took: % .6e seconds \n", time);
    //print_histo(histo);
}
