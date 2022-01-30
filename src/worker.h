#ifndef WORKER_H
#define WORKER_H

#include <pthread.h>

#include "bnb.h"
#include "consts.h"

pthread_t threads[MAX_N_THREADS];
int n_working_threads;

struct thread_work_arg_t
{
    int thread_idx;
    int n;
    int m;
    double *x;
} thread_work_args[MAX_N_THREADS];

double *thread_problem_copies_a[MAX_N_THREADS][MAX_DIMENSION];
double thread_problem_copies_a_data[MAX_N_THREADS][MAX_DIMENSION * MAX_DIMENSION];
double thread_problem_copies_b[MAX_N_THREADS][MAX_DIMENSION];
double thread_problem_copies_c[MAX_N_THREADS][MAX_DIMENSION];

bool waiting_for_work[MAX_N_THREADS];
int sent_n_jobs[MAX_N_THREADS];
pthread_cond_t work_cond_vars[MAX_N_THREADS];
pthread_mutex_t work_mutex[MAX_N_THREADS];
int done;
pthread_mutex_t done_lock;
node_t *next_node[MAX_N_THREADS];
struct thread_output
{
    node_t *node1;
    node_t *node2;
} thread_outputs[MAX_N_THREADS];

void thread_work_loop(void *arg);

#endif /* WORKER_H */