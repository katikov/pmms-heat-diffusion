#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <sys/time.h>
// #include "pipesort.h"
enum
{
    INITIALIZE,
    COMPARE,
    END
} State;
// global thread attributes
pthread_attr_t attr;
pthread_t thread;
// size of bounded buffer between threads
int buffer_size = 1;
int add = 0; /* place to add next element */
int rem = 0; /* place to remove next element */
int num = 0;
pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;    /* mutex lock for buffer */
pthread_cond_t c_cons = PTHREAD_COND_INITIALIZER; /* consumer waits on this cond var */
pthread_cond_t c_prod = PTHREAD_COND_INITIALIZER; /* producer waits on this cond var */
sem_t sort;
void *successor(void *p)
{
}
void *comparator(void *p)
{
    int *receive_buffer = (int *)p;
    int receive_pos = 0;
    int receive_buffer_num = 0;
    int out_buffer_num = 0;
    int *out_buffer = malloc(sizeof(int) * buffer_size);
    int out_pos = 0;
    int store;
    int current_value;
    int state = INITIALIZE;
    while (1)
    {
        // get the value from buffer
        pthread_mutex_lock(&m);
        if (receive_buffer_num > buffer_size)
            exit(1);                              /* overflow */
        while (receive_buffer_num == buffer_size) /* block if buffer is full */
            pthread_cond_wait(&c_cons, &m);
        current_value = receive_buffer[receive_pos];
        receive_pos = (receive_pos + 1) % buffer_size;
        receive_buffer_num++;
        pthread_mutex_unlock(&m);
        pthread_cond_signal(&c_cons);

        switch (state)
        {
        case INITIALIZE:
            store = current_value;
            state = COMPARE;
            break;
        case COMPARE:
            pthread_create(&thread, &attr, comparator, (void*)&out_buffer);
            if (store > current_value)
            {
                pthread_mutex_lock(&m);
                if (out_buffer_num > buffer_size)
                    exit(1);
                while (out_buffer_num == buffer_size)
                    pthread_cond_wait(&c_prod, &m);
                out_buffer[out_pos] = current_value;
                out_pos = (out_pos + 1) % buffer_size;
                out_buffer_num++;
                pthread_cond_signal(&c_cons);
                pthread_mutex_unlock(&m);
            }
            else
            {
                pthread_mutex_lock(&m);
                if (out_buffer_num > buffer_size)
                    exit(1);
                while (out_buffer_num == buffer_size)
                    pthread_cond_wait(&c_prod, &m);
                out_buffer[out_pos] = store;
                out_pos = (out_pos + 1) % buffer_size;
                out_buffer_num++;
                pthread_cond_signal(&c_cons);
                pthread_mutex_unlock(&m);
                store = current_value;
            }
            break;
        case END:
            break;
        default:
            break;
        }
    }
}

void *generator(void *p)
{
    // create buffer
    int *buffer = malloc(sizeof(int) * buffer_size);
    // create counter to write to buffer
    int next_in = 0;
    long length = (long)p;
    pthread_create(&thread, &attr, comparator, (void *)&buffer);

    for (int i = 0; i < length; i++)
    {
        pthread_mutex_lock(&m);
        if (num > buffer_size)
            exit(1);               /* overflow */
        while (num == buffer_size) /* block if buffer is full */
            pthread_cond_wait(&c_prod, &m);
        /* if executing here, buffer not full so add element */
        buffer[next_in] = rand();
        next_in = (next_in + 1) % buffer_size;
        num++;
        pthread_cond_signal(&c_cons);
        pthread_mutex_unlock(&m);
    }
    num = 0;
    for (int i = 0; i < 2; i++)
    {
        pthread_mutex_lock(&m);
        if (num > buffer_size)
            exit(1);               /* overflow */
        while (num == buffer_size) /* block if buffer is full */
            pthread_cond_wait(&c_prod, &m);
        /* if executing here, buffer not full so add element */
        buffer[next_in] = -1;
        next_in = (next_in + 1) % buffer_size;
        num++;
        pthread_cond_signal(&c_cons);
        pthread_mutex_unlock(&m);
    }
}
int main(int argc, char *argv[])
{
    int c;
    int seed = 42;
    long length = 1e4;

    struct timespec before;
    struct timespec after;

    /* Read command-line options. */
    while ((c = getopt(argc, argv, "l:s:")) != -1)
    {
        switch (c)
        {
        case 'l':
            length = atol(optarg);
            break;
        case 's':
            seed = atoi(optarg);
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

    clock_gettime(CLOCK_MONOTONIC, &before);
    /* Do your thing here */

    // initialize pthread attributes
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    long l = length;
    sem_init(&sort, 0, 0);
    pthread_create(&thread, &attr, generator, (void *)l);
    sem_wait(&sort);
    /* Do your thing here */
    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Pipesort took: % .6e seconds \n", time);
}
