
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "consts.h"
#include "heap.h"
#include "intopt.h"

#ifdef MAIN

#define SAFE_SCAN(fstr, ...)                           \
    {                                                  \
        if (scanf(fstr, __VA_ARGS__) == EOF) {         \
            fprintf(stderr, "Failed to read input\n"); \
            exit(0);                                   \
        }                                              \
    }

void print_sol(int n, double y, double *x)
{
    printf("z = %.1lf\n", y);

    if (y < INFINITY) {
        for (int i = 0; i < n; i++) {
            printf("x%d = %.1lf\n", i, x[i]);
        }
    }
}

int main(int argc, const char **argv)
{
    int m, n;

    SAFE_SCAN("%d %d", &m, &n);

    double *c = calloc(n, sizeof(double));
    NULL_CHECK(c);

    for (int i = 0; i < n; i++) {
        SAFE_SCAN("%lf", &c[i]);
    }

    double **a = calloc(m, sizeof(double *));
    NULL_CHECK(a);

    for (int i = 0; i < m; i++) {
        a[i] = calloc(n + 1, sizeof(double));
        NULL_CHECK(a[i]);
        for (int j = 0; j < n; j++) {
            SAFE_SCAN("%lf", &a[i][j]);
        }
    }

    double *b = calloc(m, sizeof(double));
    NULL_CHECK(b);

    for (int i = 0; i < m; i++)
        SAFE_SCAN("%lf", &b[i]);

    double *x = calloc(n + 1, sizeof(double));
    NULL_CHECK(x);

    double z = intopt(m, n, &a[0], &b[0], &c[0], &x[0]);
    print_sol(n, z, x);

    /*
     * Free the heap variables
     */
    free(c);
    for (int i = 0; i < m; i++)
        free(a[i]);
    free(a);
    free(b);
    free(x);
    free(h);
    return 0;
}
#endif