#ifndef BNB_H
#define BNB_H

#include <pthread.h>
#include <stdbool.h>

#include "simplex.h"

typedef struct node_t node_t;
struct node_t
{
    int m;
    int n;
    int k;
    int h;
    double xh;
    double ak;
    double bk;
    double *min;
    double *max;
    double **a;
    double *b;
    double *x;
    double *c;
    double z;
};

/*
 * For storing the current best solution and bounding the set
 */
extern double best_yet_z;
extern pthread_mutex_t best_yet_z_mutex;
extern double bound_z;
extern pthread_mutex_t bound_set_mutex;

node_t *initial_node(int m, int n, double **a, double *b, double *c);
node_t *succ(simplex_t *s, node_t *p, int m, int n, double **a, double *b, double *c, int k,
             double ak, double bk, double *x);
node_t *reuse_succ(simplex_t *s, node_t *p, int m, int n, double **a, double *b, double *c, int k,
                   double ak, double bk, double *x);
bool branch(node_t *q);
bool integer(node_t *p);

static inline void remove_node_contents(node_t *p)
{
    int i, m;
    m = p->m;

    /*
     * Remove all array variables
     */
    for (i = 0; i < m + 1; i++)
        free(p->a[i]);
    free(p->a);
    p->a = NULL;
    free(p->b);
    free(p->x);
    free(p->c);
}

static inline void remove_node(node_t *p)
{
    if (p->a != NULL)
        remove_node_contents(p);
    free(p->min);
    free(p->max);
    free(p);
}

#endif /* BNB_H */