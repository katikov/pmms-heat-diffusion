#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <sys/time.h>
#include <errno.h>

#define END_FLAG -1

typedef enum Order {ASCENDING, DESCENDING, RANDOM} Order;
typedef enum State{
    COMPARE,
    FORWARD
} State;

void die(const char *msg){
    if (errno != 0) 
        perror(msg);
    else
        fprintf(stderr, "error: %s\n", msg);
    exit(1);
}   

typedef struct ComparatorPackage
{
    pthread_mutex_t *m;
    int *buffer;
    pthread_cond_t *c_pro;
    pthread_cond_t *c_cons;
    int num;
    
} ComparatorPackage;

Order order = RANDOM;
pthread_attr_t attr;


int buffer_size = 50;


// ComparatorPackage* newComparatorPackage(){
//     ComparatorPackage* res = malloc(sizeof(ComparatorPackage));
//     if(res==NULL) die("malloc");
//     res->buffer = malloc(sizeof(int) * buffer_size);
//     if(res->buffer==NULL) die("malloc");
//     pthread_mutex_init(&(res->m),NULL);
//     pthread_cond_init(&(res->c_pro), NULL);
//     pthread_cond_init(&(res->c_cons), NULL);
//     res->num = 0;
    
//     return res;
// }

int push_buffer(ComparatorPackage* p, int pos, int value)
{

    pthread_mutex_lock(p->m);
    if (p->num > buffer_size)
        die("buffer error");
    
    while (p->num >= buffer_size)
        pthread_cond_wait(p->c_pro, p->m);
    p->buffer[pos] = value;
    pos = (pos + 1) % buffer_size;
    p->num++;
    pthread_cond_signal(p->c_cons);
    pthread_mutex_unlock(p->m);
    return pos;
}

int pop_buffer(ComparatorPackage* p, int pos, int *value)
{

    pthread_mutex_lock(p->m);
    if (p->num > buffer_size)
        die("buffer error");
    while (p->num<= 0)
        pthread_cond_wait(p->c_cons, p->m);
    *value = p->buffer[pos];
    pos = (pos + 1) % buffer_size;
    p->num--;
    pthread_cond_signal(p->c_pro);
    pthread_mutex_unlock(p->m);
    return pos;
}


// output thread, only receive 1 EXIT
void *output(void *p)
{
    // no producer in this thread
    ComparatorPackage *p_receive = (ComparatorPackage *)p;
    int current_value;
    int next_in = 0;
    for(;;)
    {
        next_in = pop_buffer(p_receive, next_in, &current_value);
        if (current_value == END_FLAG)
            break;
        printf("%d ", current_value);
    }
    printf("\n");
    return NULL;
}

void *comparator(void *p)
{
    ComparatorPackage *p_receive = (ComparatorPackage *)p;


    int buffer[buffer_size];
    int next_in = 0;
    int next_out = 0;

    pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t c_cons = PTHREAD_COND_INITIALIZER;
    pthread_cond_t c_pro = PTHREAD_COND_INITIALIZER;
    ComparatorPackage p_next;
    p_next.buffer = buffer;
    p_next.c_cons = &c_cons;
    p_next.c_pro = &c_pro;
    p_next.num = 0;
    p_next.m = &m;
    pthread_t thread;


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
                    pthread_create(&thread, &attr, output, &p_next);
                }else{
                    next_out = push_buffer(&p_next, next_out, current_value);
                }
                next_out = push_buffer(&p_next, next_out, store);
            }else{ // compare and forward
                if(thread_created==0){
                    thread_created = 1;
                    pthread_create(&thread, &attr, comparator, &p_next);
                }
                if(current_value>store){
                    int temp = store; store = current_value; current_value = temp;
                }
                next_out = push_buffer(&p_next, next_out, current_value);
            }
        }else{ // forward mode
            next_out = push_buffer(&p_next, next_out, current_value);
            if(current_value == END_FLAG)// second end
                break;
        }

    }
    void* result;
    if(thread_created)
        pthread_join(thread, &result);
    return NULL;
}


void *generator(void *p)
{

    int length = *(int *)p;
    int buffer[buffer_size];
    int next_out = 0;

    pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t c_cons = PTHREAD_COND_INITIALIZER;
    pthread_cond_t c_pro = PTHREAD_COND_INITIALIZER;
    ComparatorPackage p_next;
    p_next.buffer = buffer;
    p_next.c_cons = &c_cons;
    p_next.c_pro = &c_pro;
    p_next.num = 0;
    p_next.m = &m;
    pthread_t thread;

    pthread_create(&thread, &attr, comparator, &p_next);
    //pthread_create(&thread, &attr, output, &p_next);

    for (int i = 0; i < length; i++){
        int value = rand();
        next_out = push_buffer(&p_next, next_out, value);
    }
    for (int i = 0; i < 2; i++)
        next_out = push_buffer(&p_next, next_out, END_FLAG);

    void* result;
    pthread_join(thread, &result);
    return NULL;
}

int main(int argc, char *argv[]){
    int c;
    int seed = 42;
    long length = 1e4;
    long stack_size = 16384;

    struct timespec before;
    struct timespec  after;



    /* Read command-line options. */
    while((c = getopt(argc, argv, "l:s:b:t:")) != -1) {
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
            case '?':
                fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                return -1;
            default:
                return -1;
        }
    }

    /* Seed such that we can always reproduce the same random vector */
    srand(seed);

    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    pthread_attr_setstacksize(&attr, stack_size);
    //pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    int l = length;

    clock_gettime(CLOCK_MONOTONIC, &before);

    pthread_t thread;
    pthread_create(&thread, &attr, generator, &l);
    void* result;
    pthread_join(thread, &result);


    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Pipesort took: % .6e seconds \n", time);

}
