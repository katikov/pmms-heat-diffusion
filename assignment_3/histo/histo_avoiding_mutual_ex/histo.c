#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>
int num_rows;
int num_cols;
int num_threads;
int workload_per_thread;
struct thread_data
{
   int	thread_id;
   int  *histogram;
   int  *image;
   int  start;
   int  end;
};

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
		// printf("[%d]: %d\n",i+1, histo[i]);
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

void *histogram(void* threadarg){
    //TODO: For Students
    int taskid, start, end;
    int *histogram, *image;
    struct thread_data *my_data;
    my_data = (struct thread_data *) threadarg;

    taskid  = my_data->thread_id;
    start   = my_data->start;
    end     = my_data->end;
    image   = my_data->image;
    histogram = my_data->histogram;

    for (int i = start; i <= end; ++i){	
		histogram[image[i]]++;
	}
    pthread_exit((void*) histogram);
}

int main(int argc, char *argv[]){
    int c;
    int seed = 42;
    const char *image_path = 0;
    image_path ="../../../images/pat1_100x150.pgm";
    int gen_image = 0;
    int debug = 0;

    num_rows = 150;
    num_cols = 100;
    num_threads = 1;

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
    /* Create pthread and attributes */
    pthread_t p_threads[ num_threads ];
    pthread_attr_t attr;
    pthread_attr_init ( &attr ); 
    pthread_attr_setscope ( &attr , PTHREAD_SCOPE_SYSTEM );
    
    /* Create a list of pointers to store thread return values */
    void *results[num_threads];
    struct thread_data thread_data_array[num_threads];
    workload_per_thread = num_cols * num_rows / num_threads;
    int * image = (int *) malloc(sizeof(int) * num_cols * num_rows);

    /* Seed such that we can always reproduce the same random vector */
    if (gen_image){
    	srand(seed);
    	generate_image(num_rows, num_cols, image);
    }else{
    	read_image(image_path,num_rows, num_cols, image);
    }

    /* Intial thread data structure */
    int count = 0;
    int i;
    for(i=0; i<num_threads; i++){
        thread_data_array[i].thread_id = i;
        thread_data_array[i].histogram = (int*) calloc(256, sizeof(int));
        thread_data_array[i].image = image;
        thread_data_array[i].start = count;
        if(i==num_threads-1)
            thread_data_array[i].end = num_cols*num_rows-1;
        else
            thread_data_array[i].end = count+workload_per_thread-1;
        
        // printf("%d %d\n", thread_data_array[i].start, thread_data_array[i].end);
        count+=workload_per_thread;
    }

    clock_gettime(CLOCK_MONOTONIC, &before);
    /* Do your thing here */

    for(i=0;i<num_threads;i++){
        pthread_create(&p_threads[i], &attr, &histogram, (void *) 
       &thread_data_array[i]); 
    }
    for(i=0;i<num_threads;i++){
        pthread_join(p_threads[i],
                     &results[i]);
    }

    /* Do your thing here */

    /* Combine thread return histogram */
    for(i=0; i<num_threads; i++){
        for(int j=0; j<256; j++){
            histo[j] = histo[j] + ((int*)results[i])[j];
        }
    }

    if (debug){
    	print_histo(histo);
    }

    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Histo took: % .6e seconds \n", time);
}