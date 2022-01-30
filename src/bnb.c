#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bnb.h"
#include "consts.h"
#include "simplex.h"

double best_yet_z;
pthread_mutex_t best_yet_z_mutex = PTHREAD_MUTEX_INITIALIZER;
double bound_z = -INFINITY;
pthread_mutex_t bound_set_mutex = PTHREAD_MUTEX_INITIALIZER;

/* 
 * Branch and bound algorithm
 */

node_t *initial_node(int m, int n, double **a, double *b, double *c)
{
    node_t *p = calloc(1, sizeof(node_t));

    p->a = calloc(m + 1, sizeof(double *));
    NULL_CHECK(p->a);
    for (int i = 0; i < m + 1; i++) {
        p->a[i] = calloc(n + 1, sizeof(double));
        NULL_CHECK(p->a[i]);
    }

    p->b = calloc(m + 1, sizeof(double));
    NULL_CHECK(p->b);

    p->c = calloc(n + 1, sizeof(double));
    NULL_CHECK(p->c);

    p->x = calloc(n + 1, sizeof(double));
    NULL_CHECK(p->x);

    p->min = calloc(n, sizeof(double));
    NULL_CHECK(p->min);

    p->max = calloc(n, sizeof(double));
    NULL_CHECK(p->max);

    p->m = m;
    p->n = n;

    for (int i = 0; i < m; i++)
        memcpy(p->a[i], a[i], (n + 1) * sizeof(double));

    memcpy(p->b, b, n * sizeof(double));
    memcpy(p->c, c, n * sizeof(double));

    for (int i = 0; i < n; i++) {
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY;
    }

    return p;
}

node_t *extend(node_t *p, int m, int n, double **a, double *b, double *c, int k,
               double ak, double bk)
{
    node_t *q = calloc(1, sizeof(node_t));
    int i, j;

    q->k = k;
    q->ak = ak;
    q->bk = bk;
    if (ak > 0 && p->max[k] < INFINITY)
        q->m = p->m;
    else if (ak < 0 && p->min[k] > 0)
        q->m = p->m;
    else
        q->m = p->m + 1;
    q->n = p->n;
    q->h = -1;
    q->z = p->z;
    q->a = calloc(q->m + 1, sizeof(double *));
    NULL_CHECK(q->a);
    for (int i = 0; i < q->m + 1; i++) {
        q->a[i] = calloc(q->n + 1, sizeof(double));
        NULL_CHECK(q->a[i]);
    }

    q->b = calloc(q->m + 1, sizeof(double));
    NULL_CHECK(q->b);

    q->c = calloc(q->n + 1, sizeof(double));
    NULL_CHECK(q->c);

    q->x = calloc(q->n + 1, sizeof(double));
    NULL_CHECK(q->x);

    q->min = calloc(n, sizeof(double));
    NULL_CHECK(q->min);

    q->max = calloc(n, sizeof(double));
    NULL_CHECK(q->max);

    memcpy(q->min, p->min, n * sizeof(double));
    memcpy(q->max, p->max, n * sizeof(double));

    for (int i = 0; i < m; i++)
        memcpy(q->a[i], a[i], (n + 1) * sizeof(double));

    memcpy(q->b, b, m * sizeof(double));
    memcpy(q->c, c, n * sizeof(double));

    if (ak > 0) {
        if (q->max[k] == INFINITY || bk < q->max[k])
            q->max[k] = bk;
    } else if (q->min[k] == -INFINITY || -bk > q->min[k]) {
        q->min[k] = -bk;
    }

    for (i = m, j = 0; j < n; j++) {
        if (q->min[j] > -INFINITY) {
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i++;
        }

        if (q->max[j] < INFINITY) {
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i++;
        }
    }

    return q;
}

node_t *reuse_expand(node_t *p, int m, int n, double **a, double *b, double *c, int k,
                     double ak, double bk)
{
    int i, j;
    int p_m = p->m;
    p->k = k;
    p->ak = ak;
    p->bk = bk;
    if ((ak > 0 && p->max[k] < INFINITY) || (ak < 0 && p->min[k] > 0)) {
        p->a = realloc(p->a, (p->m + 1) * sizeof(double *));
        NULL_CHECK(p->a);

        for (int i = m; i < p_m + 1; i++) {
            memset(p->a[i], 0, (n + 1) * sizeof(double));
        }
    } else {
        p->m++;
        p->a = realloc(p->a, (p->m + 1) * sizeof(double *));
        NULL_CHECK(p->a);

        for (int i = m; i < p_m + 1; i++) {
            memset(p->a[i], 0, (n + 1) * sizeof(double));
        }
        p->a[p_m + 1] = calloc(n + 1, sizeof(double));
    }

    p->h = -1;

    free(p->b);
    p->b = calloc(p->m + 1, sizeof(double));
    NULL_CHECK(p->b);

    free(p->c);
    p->c = calloc(p->n + 1, sizeof(double));
    NULL_CHECK(p->c);

    free(p->x);
    p->x = calloc(p->n + 1, sizeof(double));
    NULL_CHECK(p->x);

    memcpy(p->b, b, m * sizeof(double));
    memcpy(p->c, c, n * sizeof(double));
    for (int i = 0; i < m; i++)
        memcpy(p->a[i], a[i], (n + 1) * sizeof(double));

    if (ak > 0) {
        if (p->max[k] == INFINITY || bk < p->max[k])
            p->max[k] = bk;
    } else if (p->min[k] == -INFINITY || -bk > p->min[k]) {
        p->min[k] = -bk;
    }

    for (i = m, j = 0; j < n; j++) {
        if (p->min[j] > -INFINITY) {
            p->a[i][j] = -1;
            p->b[i] = -p->min[j];
            i++;
        }

        if (p->max[j] < INFINITY) {
            p->a[i][j] = 1;
            p->b[i] = p->max[j];
            i++;
        }
    }

    return p;
}

bool isinteger(double *xp)
{
    double x = *xp;
    double r = lround(x);
    if (fabs(r - x) < EPSILON) {
        *xp = r;
        return true;
    } else {
        return false;
    }
}

bool integer(node_t *p)
{
    for (int i = 0; i < p->n; i++)
        if (!isinteger(&p->x[i]))
            return 0;
    return 1;
}

void bound(node_t *p, double *x)
{

    pthread_mutex_lock(&best_yet_z_mutex);
    /*
     * The if statement below might seem redundant, but it is not.
     * Especially on smaller simplices, sometimes two threads find
     * a better solution at the same time, and if the first thread
     * has the better value, it would otherwise be overwritten by
     * the second thread. 
     */
    if (p->z <= best_yet_z) {
        pthread_mutex_unlock(&best_yet_z_mutex);
        return;
    }
    best_yet_z = p->z;
    pthread_mutex_unlock(&best_yet_z_mutex);
    memcpy(x, p->x, p->n * sizeof(double));
    pthread_mutex_lock(&bound_set_mutex);
    bound_z = p->z;
    pthread_mutex_unlock(&bound_set_mutex);
}

bool branch(node_t *q)
{
    double min, max;

    for (int h = 0; h < q->n; h++) {
        if (!isinteger(&q->x[h])) {
            if (q->min[h] == -INFINITY)
                min = 0;
            else
                min = q->min[h];
            max = q->max[h];
            if (floor(q->x[h]) < min || ceil(q->x[h]) > max)
                continue;
            q->h = h;
            q->xh = q->x[h];
            return 1;
        }
    }
    return 0;
}

node_t *succ(simplex_t *s, node_t *p, int m, int n, double **a, double *b, double *c, int k,
             double ak, double bk, double *x)
{
    node_t *q = extend(p, m, n, a, b, c, k, ak, bk);
    if (q == NULL) {
        return NULL;
    }
    q->z = simplex(s, q->m, q->n, q->a, q->b, q->c, q->x, 0);

    pthread_mutex_lock(&best_yet_z_mutex);
    bool might_be_best = q->z > best_yet_z;
    pthread_mutex_unlock(&best_yet_z_mutex);
    if (might_be_best) {
        if (integer(q)) {
            bound(q, x);
        } else if (branch(q)) {
            return q;
        }
    }
    remove_node(q);
    return NULL;
}

node_t *reuse_succ(simplex_t *s, node_t *p, int m, int n, double **a, double *b, double *c, int k,
                   double ak, double bk, double *x)
{
    node_t *q = reuse_expand(p, m, n, a, b, c, k, ak, bk);
    if (q == NULL) {
        return NULL;
    }
    q->z = simplex(s, q->m, q->n, q->a, q->b, q->c, q->x, 0);
    pthread_mutex_lock(&best_yet_z_mutex);
    bool might_be_best = q->z > best_yet_z;
    pthread_mutex_unlock(&best_yet_z_mutex);

    if (might_be_best) {
        if (integer(q)) {
            bound(q, x);
        } else if (branch(q)) {
            return q;
        }
    }
    remove_node(q);
    return NULL;
}