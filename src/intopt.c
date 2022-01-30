#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bnb.h"
#include "consts.h"
#include "heap.h"
#include "simplex.h"
#include "worker.h"

static inline void spawn_thread(int n, int m, double **a, double *b, double *c, double *x)
{
    int thread_idx = n_working_threads;

    for (int j = 0; j < m; j++) {
        memcpy(thread_problem_copies_a[thread_idx][j], a[j], (n + 1) * sizeof(double));
    }
    memcpy(thread_problem_copies_b[thread_idx], b, m * sizeof(double));
    memcpy(thread_problem_copies_c[thread_idx], c, n * sizeof(double));

    thread_work_args[thread_idx] = (struct thread_work_arg_t){thread_idx, n, m, x};

    pthread_create(&threads[thread_idx], NULL, (void *(*)(void *))thread_work_loop, &thread_work_args[thread_idx]);
    n_working_threads++;
}

static inline bool solve_done()
{
    for (int i = 0; i < n_working_threads; i++) {
        pthread_mutex_lock(&work_mutex[i]);
        if (!waiting_for_work[i] ||
            thread_outputs[i].node1 != NULL ||
            thread_outputs[i].node2 != NULL) {
            pthread_mutex_unlock(&work_mutex[i]);
            return false;
        }
        pthread_mutex_unlock(&work_mutex[i]);
    }

    done = 0;
    for (int i = 0; i < n_working_threads; i++) {
        pthread_mutex_lock(&work_mutex[i]);
        waiting_for_work[i] = false;
        pthread_cond_signal(&work_cond_vars[i]);
        pthread_mutex_unlock(&work_mutex[i]);
    }

    return true;
}

static inline void reduce_set()
{
    pthread_mutex_lock(&bound_set_mutex);
    if (bound_z > -INFINITY) {
        bound_set(bound_z);
        bound_z = -INFINITY;
    }
    pthread_mutex_unlock(&bound_set_mutex);
}

static inline void harvest_nodes_and_feed_threads()
{
    for (int i = 0; i < n_working_threads; i++) {
        pthread_mutex_lock(&work_mutex[i]);
        bool is_waiting = waiting_for_work[i];
        node_t *output_node1 = thread_outputs[i].node1;
        node_t *output_node2 = thread_outputs[i].node2;
        thread_outputs[i].node1 = NULL;
        thread_outputs[i].node2 = NULL;
        pthread_mutex_unlock(&work_mutex[i]);

        if (is_waiting) {
            if (output_node1 != NULL)
                add_to_set(output_node1);

            if (output_node2 != NULL)
                add_to_set(output_node2);

            if (set_nonempty()) {
                pthread_mutex_lock(&work_mutex[i]);
                next_node[i] = next_from_set();
                waiting_for_work[i] = false;
                pthread_cond_signal(&work_cond_vars[i]);
                pthread_mutex_unlock(&work_mutex[i]);
            }
        }
    }
}

bool first_solve = true;
double intopt(int m, int n, double **a, double *b, double *c, double *x)
{
    node_t *p = initial_node(m, n, a, b, c);
    best_yet_z = -INFINITY;
    bound_z = -INFINITY;
    p->z = simplex(&s[0], p->m, p->n, p->a, p->b, p->c, p->x, 0);
    if (integer(p) || !isfinite(p->z)) {
        best_yet_z = p->z;
        if (integer(p))
            memcpy(x, p->x, p->n * sizeof(double));
        remove_node(p);
        return best_yet_z;
    }

    if (first_solve) {
        for (int i = 0; i < MAX_N_THREADS; i++)
            pthread_cond_init(&work_cond_vars[i], NULL);

        for (int i = 0; i < MAX_N_THREADS; i++)
            pthread_mutex_init(&work_mutex[i], NULL);
        for (int i = 0; i < MAX_N_THREADS; i++)
            for (int j = 0; j < MAX_DIMENSION; j++)
                thread_problem_copies_a[i][j] = &thread_problem_copies_a_data[i][MAX_DIMENSION * j];
        first_solve = false;
    }

    if (h == NULL) {
        /*
         * NOTE: the heap is not freed properly when the code is compiled without
         * the custom main function since we don't have control over the code in
         * the main function, and it would be inefficient to reallocate the heap
         * every time we call intopt.
         */
        h = calloc(INIT_HEAP_SIZE, sizeof(struct heap_elem));
        heap_size = INIT_HEAP_SIZE;
    }
    heap_idx = 0;

    int n_original_threads = 1;
    switch (n) {
    case 1:
        n_original_threads = 1;
        break;
    case 2:
        n_original_threads = 1;
        break;
    case 3:
        n_original_threads = 1;
        break;
    case 4:
        n_original_threads = 1;
        break;
    case 5:
        n_original_threads = 1;
        break;
    case 6:
        n_original_threads = 1;
        break;
    case 7:
        n_original_threads = 1;
        break;
    case 8:
        n_original_threads = 1;
        break;
    case 9:
        n_original_threads = 6;
        break;
    case 10:
        n_original_threads = 6;
        break;
    case 11:
        n_original_threads = 6;
        break;
    case 12:
        n_original_threads = 8;
        break;
    case 13:
        n_original_threads = 10;
        break;
    case 14:
        n_original_threads = 14;
        break;
    case 15:
        n_original_threads = 15;
        break;
    case 16:
        n_original_threads = 16;
        break;
    case 17:
        n_original_threads = 34;
        break;
    case 18:
        n_original_threads = 45;
        break;
    case 19:
        n_original_threads = MAX_N_THREADS - 29;
        break;
    default:
        n_original_threads = MAX_N_THREADS - 23;
        break;
    }

    done = 1;

    for (int i = 0; i < n_original_threads; i++)
        spawn_thread(n, m, a, b, c, x);

    add_to_set(p);
    branch(p);
    while (true) {
        if (n_working_threads == MAX_N_THREADS)
            goto max_thread_loop;
        if (heap_idx > 8 * n_working_threads)
            spawn_thread(n, m, a, b, c, x);

        harvest_nodes_and_feed_threads();
        if (!set_nonempty())
            if (solve_done())
                break;
        reduce_set();
    }

    if (false) {
    max_thread_loop:
        while (true) {
            harvest_nodes_and_feed_threads();
            if (!set_nonempty())
                if (solve_done())
                    break;
            reduce_set();
        }
    }

    for (int i = 0; i < n_working_threads; i++)
        pthread_join(threads[i], NULL);
    n_working_threads = 0;

    if (best_yet_z == -INFINITY)
        return NAN;
    else
        return best_yet_z;
}