#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <sys/time.h>
typedef struct
{
    pthread_mutex_t *m;
    int *num;
    int *buffer;
    pthread_cond_t *c_pro;
    pthread_cond_t *c_cons;
} ComparatorPackage;

enum
{
    INITIALIZE,
    RECEIVE_BEFORE_SUCCESSOR,
    NORMAL_COMPARE,
    END
} State;

// global thread attributes
pthread_attr_t attr;
pthread_t thread;
int buffer_size = 50;

int push_buffer(int *buffer, pthread_mutex_t *m, pthread_cond_t *c_cons, pthread_cond_t *c_pro, int value, int *num, int pos)
{

    pthread_mutex_lock(m);
    if (*num > buffer_size)
        exit(1);
    while (*num == buffer_size)
        pthread_cond_wait(c_pro, m);
    buffer[pos] = value;
    pos = (pos + 1) % buffer_size;
    (*num)++;
    pthread_cond_signal(c_cons);
    pthread_mutex_unlock(m);
    return pos;
}
void *output(void *p)
{
    // no producer in this thread
    ComparatorPackage *p_receive = (ComparatorPackage *)p;
    int current_value;
    int receive_pos;
    int *num = p_receive->num;
    int state = INITIALIZE;
    while (1)
    {
        pthread_mutex_lock(p_receive->m);
        if (*num > buffer_size)
            exit(1);                /* overflow */
        while (*num == 0) /* block if buffer is full */
            pthread_cond_wait(p_receive->c_cons, p_receive->m);
        current_value = p_receive->buffer[receive_pos];
        receive_pos = (receive_pos + 1) % buffer_size;
        (*num)--;
        pthread_mutex_unlock(p_receive->m);
        pthread_cond_signal(p_receive->c_pro);

        if (current_value == -1)
        {
            if (state == INITIALIZE)
            {
                state = END;
                continue;
            }
            else
            {
                break;
            }
        }
        printf("%i\n", current_value);
    }
}
void *comparator(void *p)
{
    fprintf(stderr, "c");
    ComparatorPackage *p_receive = (ComparatorPackage *)p;
    int receive_pos = 0;
    int *num = p_receive->num;
    int *out_buffer = malloc(sizeof(int) * buffer_size);
    int *out_num = malloc(sizeof(int));
    int out_pos = 0;
    int store; // local store data
    int current_value;
    int state = INITIALIZE; // first state
    pthread_cond_t *out_c_pro = malloc(sizeof(pthread_cond_t));
    pthread_cond_t *out_c_cons = malloc(sizeof(pthread_cond_t));
    pthread_mutex_t *out_m = malloc(sizeof(pthread_mutex_t));
    pthread_cond_init(out_c_pro, NULL);
    pthread_cond_init(out_c_cons, NULL);
    pthread_mutex_init(out_m, NULL);
    ComparatorPackage p_next;
    p_next.buffer = out_buffer;
    p_next.c_cons = out_c_cons;
    p_next.c_pro = out_c_pro;
    p_next.num = out_num;
    p_next.m = out_m;
    while (1)
    {
        fprintf(stderr, "c2");
        pthread_mutex_lock(p_receive->m);
        if (*num > buffer_size)
            exit(1);
        while (*num == 0)
            pthread_cond_wait(p_receive->c_cons, p_receive->m);
        current_value = p_receive->buffer[receive_pos];
        receive_pos = (receive_pos + 1) % buffer_size;
        (*num)--;
        pthread_mutex_unlock(p_receive->m);
        pthread_cond_signal(p_receive->c_pro);
        switch (state)
        {
        case INITIALIZE:
            store = current_value;
            state = NORMAL_COMPARE;
            break;
        case RECEIVE_BEFORE_SUCCESSOR:
        {
            if (current_value == -1)
            {
                pthread_create(&thread, &attr, output, &p_next);
                state = END;
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, current_value, p_next.num, out_pos);
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, store, p_next.num, out_pos);
                continue;
            }
            pthread_create(&thread, &attr, comparator, &p_next);
            state = NORMAL_COMPARE;
            if (store > current_value)
            {
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, current_value, p_next.num, out_pos);
            }
            else
            {
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, store, p_next.num, out_pos);
                store = current_value;
            }
        }
        case NORMAL_COMPARE:
        {
            pthread_create(&thread, &attr, comparator, &p_next);
            if (current_value == -1)
            {
                state = END;
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, current_value, p_next.num, out_pos);
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, store, p_next.num, out_pos);
                store = current_value;
                continue;
            }
            if (store > current_value)
            {
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, current_value, p_next.num, out_pos);
            }
            else
            {
                out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, store, p_next.num, out_pos);
                store = current_value;
            }
            break;
        }
        case END:
        {
            out_pos = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, current_value, p_next.num, out_pos);
            if (current_value == -1)
            {
                break;
            }
            break;
        }
        default:
            break;
        }
    }
}

void *generator(void *p)
{

    int length = *(int *)p;
    int *buffer = malloc(sizeof(int) * buffer_size);
    int next_out = 0;
    int num = 0;//represent current number in buffer

    pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t c_cons = PTHREAD_COND_INITIALIZER;
    pthread_cond_t c_pro = PTHREAD_COND_INITIALIZER;
    ComparatorPackage p_next;
    p_next.buffer = buffer;
    p_next.c_cons = &c_cons;
    p_next.c_pro = &c_pro;
    p_next.num = &num;
    p_next.m = &m;
    pthread_create(&thread, &attr, comparator, &p_next);

    for (int i = 0; i < length; i++)
    {
        next_out = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, rand(), p_next.num, next_out);
    }
    for (int i = 0; i < 2; i++)
    {

        next_out = push_buffer(p_next.buffer, p_next.m, p_next.c_cons, p_next.c_pro, -1, p_next.num, next_out);
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
    // initialize pthread attributes
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
<<<<<<< HEAD

=======
>>>>>>> 65a0fb904c61457880ab7cf7117e6dcdec192144
    int l = length;

    clock_gettime(CLOCK_MONOTONIC, &before);

    pthread_create(&thread, &attr, generator, &l);
<<<<<<< HEAD

=======
>>>>>>> 65a0fb904c61457880ab7cf7117e6dcdec192144
    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Pipesort took: % .6e seconds \n", time);
}
