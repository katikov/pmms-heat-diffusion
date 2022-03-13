#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <sys/time.h>
#include <errno.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define END_FLAG -1
#define MAX_THREADS 4080



void die(const char *msg){
    if (errno != 0) 
        perror(msg);
    else
        fprintf(stderr, "error: %s\n", msg);
    exit(1);
}   

typedef struct ComparatorPackage
{
    pthread_mutex_t m;
    int *buffer;
    pthread_cond_t c_pro;
    pthread_cond_t c_cons;
    pthread_t thread;
    int num;   // number of elements in buffer
    struct ComparatorPackage* prev; // linked list, prev->thread is the predecessor
} ComparatorPackage;

pthread_attr_t attr;


int buffer_size = 50;
int length = 1e4;
int debug = 0;

// join the output thread in main
ComparatorPackage* p_final = NULL;
pthread_mutex_t final_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t final_cond = PTHREAD_COND_INITIALIZER;


// total number of existing threads <= 4080
int total_threads=0;
pthread_mutex_t create_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t create_cond = PTHREAD_COND_INITIALIZER;

ComparatorPackage* newComparatorPackage(){
    pthread_mutex_lock(&create_mutex);
    while(total_threads>=MAX_THREADS)
        pthread_cond_wait(&create_cond, &create_mutex);
    total_threads++;
    pthread_mutex_unlock(&create_mutex);

    ComparatorPackage* res = malloc(sizeof(ComparatorPackage));
    if(res==NULL) die("malloc");
    res->buffer = malloc(sizeof(int) * buffer_size);
    if(res->buffer==NULL) die("malloc");
    pthread_mutex_init(&(res->m),NULL);
    pthread_cond_init(&(res->c_pro), NULL);
    pthread_cond_init(&(res->c_cons), NULL);
    res->num = 0;
    res->prev = NULL;
    return res;
}

void deleteComparatorPackage(ComparatorPackage* p){
    if(p->buffer) free(p->buffer);
    void* res;
    pthread_join(p->thread, &res);
    free(p);
    pthread_mutex_lock(&create_mutex);
    total_threads--;
    pthread_cond_signal(&create_cond);
    pthread_mutex_unlock(&create_mutex);
}

int push_buffer(ComparatorPackage* p, int pos, int value)
{

    pthread_mutex_lock(&(p->m));
    if (p->num > buffer_size)
        die("buffer error");
    
    while (p->num >= buffer_size)
        pthread_cond_wait(&(p->c_pro), &(p->m));
    p->buffer[pos] = value;
    pos = (pos + 1) % buffer_size;
    p->num++;
    pthread_cond_signal(&(p->c_cons));
    pthread_mutex_unlock(&(p->m));
    return pos;
}

int pop_buffer(ComparatorPackage* p, int pos, int *value)
{

    pthread_mutex_lock(&(p->m));
    if (p->num > buffer_size)
        die("buffer error");
    while (p->num<= 0)
        pthread_cond_wait(&(p->c_cons), &(p->m));
    *value = p->buffer[pos];
    pos = (pos + 1) % buffer_size;
    p->num--;
    pthread_cond_signal(&(p->c_pro));
    pthread_mutex_unlock(&(p->m));
    return pos;
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


// output thread, only receive 1 EXIT
void *output(void *p)
{
    // no producer in this thread
    ComparatorPackage *p_receive = (ComparatorPackage *)p;
    
    int current_value;
    int next_in = 0;
    int print_cnt = 0;
    for(;;)
    {
        next_in = pop_buffer(p_receive, next_in, &current_value);
        if (current_value == END_FLAG)
            break;
        if(debug){
            if(print_cnt!=0 && (print_cnt%10==0))
                printf("\n");
            printf("%d ", current_value);
        }
        print_cnt++;
    }
    if(debug)
        printf("\n");

    deleteComparatorPackage(p_receive->prev);

    pthread_mutex_lock(&final_mutex);
    p_final = p;
    pthread_cond_signal(&final_cond);
    pthread_mutex_unlock(&final_mutex);
    return NULL;
}


void *comparator(void *p)
{
    ComparatorPackage *p_receive = (ComparatorPackage *)p;
    ComparatorPackage *p_next = newComparatorPackage();
    p_next->prev = p_receive;


    int next_in = 0;
    int next_out = 0;

    int store; // store 1 integer local data
    int current_value; 
    int state = 0; // 0: compare; 1: forward
    int thread_created = 0; // 0: tail thread; 1: not tail thread

    next_in = pop_buffer(p_receive, next_in, &store);
    

    for(;;){
        next_in = pop_buffer(p_receive, next_in, &current_value);
        if(state == 0){ 
            if(current_value == END_FLAG){ //first end
                state = 1;
                if(thread_created==0){ // create output
                    thread_created = 1;
                    pthread_create(&(p_next->thread), &attr, output, p_next);
                }else{
                    next_out = push_buffer(p_next, next_out, current_value);
                }
                next_out = push_buffer(p_next, next_out, store);
            }else{ // compare and forward
                if(thread_created==0){
                    thread_created = 1;
                    pthread_create(&(p_next->thread), &attr, comparator, p_next);
                }
                if(current_value>store){
                    int temp = store; store = current_value; current_value = temp;
                }
                next_out = push_buffer(p_next, next_out, current_value);
            }
        }else{ // forward mode
            next_out = push_buffer(p_next, next_out, current_value);
            if(current_value == END_FLAG)// second end
                break;
        }

    }
    void* result;
    deleteComparatorPackage(p_receive->prev);
    return NULL;
}


void *generator(void *p)
{
    ComparatorPackage *p_prev = p;
    ComparatorPackage *p_next = newComparatorPackage();
    p_next->prev = p_prev;

    int next_out = 0;

    pthread_create(&(p_next->thread), &attr, comparator, p_next);
    //pthread_create(&(p_next->thread), &attr, output, p_next);

    for (int i = 0; i < length; i++){
        int value = rand();
        next_out = push_buffer(p_next, next_out, value);
    }
    for (int i = 0; i < 2; i++)
        next_out = push_buffer(p_next, next_out, END_FLAG);

    return NULL;
}

int main(int argc, char *argv[]){
    int c;
    int seed = 42;
    //long length = 1e4;
    long stack_size = 16384;

    struct timespec before;
    struct timespec  after;



    /* Read command-line options. */
    while((c = getopt(argc, argv, "l:s:b:t:g:")) != -1) {
        switch(c) {
            case 'l':
                length = atol(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'b':
                buffer_size = atoi(optarg);
                break;
            case 't':
                stack_size = atol(optarg);
                break;
            case 'g':
                debug = atoi(optarg);
                break;
            case '?':
                fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                return -1;
            default:
                return -1;
        }
    }

    buffer_size = MAX(buffer_size, length/(MAX_THREADS-2));


    /* Seed such that we can always reproduce the same random vector */
    srand(seed);

    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    pthread_attr_setstacksize(&attr, stack_size);
    //pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

    clock_gettime(CLOCK_MONOTONIC, &before);

    ComparatorPackage *p = newComparatorPackage();

    pthread_create(&(p->thread), &attr, generator, p);
    pthread_mutex_lock(&final_mutex);
    while(p_final == NULL)
        pthread_cond_wait(&final_cond, &final_mutex);
    deleteComparatorPackage(p_final);
    pthread_mutex_unlock(&final_mutex);

    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Pipesort took: % .6e seconds \n", time);

}
