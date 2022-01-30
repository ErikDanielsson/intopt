#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "consts.h"
#include "simplex.h"

int init(simplex_t *s, int m, int n, double **a, double *b, double *c,
         double *x, double y, int *var)
{
    s->m = m;
    s->n = n;
    s->a = a;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;

    if (var == NULL) {
        var = calloc(m + n + 1, sizeof(int));
        for (int i = 0; i < m + n + 1; i++)
            var[i] = i;
    }

    s->var = var;

    int k = 0;
    for (int i = 0; i < m; i++) {
        if (b[i] < b[k])
            k = i;
    }

    return k;
}

static inline int select_nonbasic(simplex_t *s)
{
    int n = s->n;
    for (int i = 0; i < n; i++) {
        if (s->c[i] > EPSILON) {
            return i;
        }
    }
    return -1;
}

void pivot(simplex_t *s, int row, int col);
void prepare(simplex_t *s, int k)
{
    int m = s->m;
    int n = s->n;
    for (int i = m + n; i > n; i--) {
        s->var[i] = s->var[i - 1];
    }
    s->var[n] = n + m;

    n++;
    for (int i = 0; i < m; i++) {
        s->a[i][n - 1] = -1;
    }

    s->x = calloc(n + m, sizeof(double));
    s->c = calloc(n, sizeof(double));

    s->c[n - 1] = -1;

    s->n = n;
    pivot(s, k, n - 1);
}

double xsimplex(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h);

bool initial(simplex_t *s, int m, int n, double **a, double *b, double *c,
             double *x, double y, int *var)
{
    int i, j, k;
    double w;
    k = init(s, m, n, a, b, c, x, y, var);
    if (b[k] >= 0)
        return 1;
    prepare(s, k);
    n = s->n;

    s->y = xsimplex(s, m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

    for (i = 0; i < m + n; i++) {
        if (s->var[i] == m + n - 1) {
            if (fabs(s->x[i]) > EPSILON) {
                free(s->x);
                free(s->c);
                return 0;
            } else {
                break;
            }
        }
    }

    if (i >= n) {
        for (j = k = 0; k < n; k++)
            if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j]))
                j = k;
        pivot(s, i - n, j);
        i = j;
    }

    if (i < n - 1) {
        k = s->var[i];
        s->var[i] = s->var[n - 1];
        s->var[n - 1] = k;
        for (k = 0; k < m; k++) {
            w = s->a[k][n - 1];
            s->a[k][n - 1] = s->a[k][i];
            s->a[k][i] = w;
        }
    }
    free(s->c);
    s->c = c;
    s->y = y;
    for (k = n - 1; k < n + m - 1; k++)
        s->var[k] = s->var[k + 1];
    n = s->n = s->n - 1;

    double *t = calloc(n, sizeof(double));
    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            if (k == s->var[j]) {
                t[j] += s->c[k];
                goto outer;
            }
        }
        for (j = 0; j < m; j++) {
            if (k == s->var[n + j]) {
                break;
            }
        }
        s->y += s->c[k] * s->b[j];
        for (i = 0; i < n; i++)
            t[i] -= s->c[k] * s->a[j][i];
    outer:;
    }
    for (i = 0; i < n; i++)
        s->c[i] = t[i];
    free(t);
    free(s->x);
    return 1;
}

void pivot(simplex_t *s, int row, int col)
{
    double **a = s->a;
    double *b = s->b;
    double *c = s->c;
    int n = s->n;
    int m = s->m;
    double piv_val = 1 / a[row][col];
    int tmp;
    tmp = s->var[col];
    s->var[col] = s->var[n + row];
    s->var[n + row] = tmp;

    double b_row = b[row], c_col = c[col];
    b[row] *= piv_val;
    c[col] *= -piv_val;
    s->y += c_col * b_row * piv_val;

    /* This makes the loops (somewhat) branchless */
    a[row][col] = 0;

    for (int i = 0; i < m; i++) {
        if (i != row) {
            double s = -a[i][col] * piv_val;
            for (int j = 0; j < n; j++) {
                a[i][j] += s * a[row][j];
            }
            a[i][col] = s;
            b[i] += s * b_row;
        }
    }
    for (int i = 0; i < n; i++) {
        double s = a[row][i] * piv_val;
        a[row][i] = s;
        c[i] -= c_col * s;
    }

    a[row][col] = piv_val;
}

double xsimplex(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h)
{
    int i, row, col;
    if (!initial(s, m, n, a, b, c, x, y, var)) {
        free(s->var);
        return NAN;
    }

    while ((col = select_nonbasic(s)) >= 0) {
        row = -1;
        double a_p = 1;
        double b_p = 1;

        for (i = 0; i < m; i++) {
            if (a[i][col] > EPSILON) {
                row = i;
                a_p = a[i][col];
                b_p = b[i];
                break;
            }
        }
        /*
         * We are allowed to replace the divisions with multiplications
         * in the inequality since all variables are positive
         */
        for (; i < m; i++) {
            if (a[i][col] > EPSILON && b[i] * a_p < b_p * a[i][col]) {
                row = i;
                a_p = a[i][col];
                b_p = b[i];
            }
        }

        if (row < 0) {
            free(s->var);
            return INFINITY;
        }
        pivot(s, row, col);
    }

    if (h == 0) {
        for (int i = 0; i < n; i++) {
            if (s->var[i] < n) {
                x[s->var[i]] = 0;
            }
        }

        for (int i = 0; i < m; i++) {
            if (s->var[n + i] < n) {
                x[s->var[n + i]] = s->b[i];
            }
        }
        free(s->var);
    } else {
        for (int i = 0; i < n; i++) {
            x[i] = 0;
        }

        for (int i = n; i < n + m; i++) {
            x[i] = s->b[i - n];
        }
    }
    return s->y;
}

double simplex(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x,
               double y)
{
    return xsimplex(s, m, n, a, b, c, x, y, NULL, 0);
}