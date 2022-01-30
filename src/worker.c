#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "bnb.h"
#include "simplex.h"
#include "worker.h"

void thread_work_loop(void *arg)
{
    struct thread_work_arg_t t = *(struct thread_work_arg_t *)arg;
    int thread_idx = t.thread_idx;
    int n = t.n;
    int m = t.m;
    double **a = &thread_problem_copies_a[thread_idx][0];
    double *b = &thread_problem_copies_b[thread_idx][0];
    double *c = &thread_problem_copies_c[thread_idx][0];
    double *x = t.x;

    simplex_t *thread_simplex = &s[thread_idx];
    pthread_mutex_lock(&work_mutex[thread_idx]);
    while (true) {
        waiting_for_work[thread_idx] = true;
        while (waiting_for_work[thread_idx]) {
            pthread_cond_wait(&work_cond_vars[thread_idx], &work_mutex[thread_idx]);
        }

        if (done == 0) {
            pthread_mutex_unlock(&work_mutex[thread_idx]);
            break;
        }

        node_t *p = next_node[thread_idx];
        next_node[thread_idx] = NULL;
        pthread_mutex_unlock(&work_mutex[thread_idx]);
        node_t *q1 = succ(thread_simplex, p, m, n, a, b, c, p->h, 1, floor(p->xh), x);
        node_t *q2 = reuse_succ(thread_simplex, p, m, n, a, b, c, p->h, -1, -ceil(p->xh), x);
        pthread_mutex_lock(&work_mutex[thread_idx]);
        thread_outputs[thread_idx] = (struct thread_output){q1, q2};
    }
}