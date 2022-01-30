#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "consts.h"

typedef struct simplex_t simplex_t;
struct simplex_t
{
    int m;
    int n;
    int *var;
    double **a;
    double *b;
    double *x;
    double *c;
    double y;
} s[MAX_N_THREADS];

double simplex(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x,
               double y);

#endif /* SIMPLEX_H */