#ifndef CONSTS_H
#define CONSTS_H

#define MAIN
#define INIT_HEAP_SIZE 1024
#define EPSILON (1.0e-6)
#define MAX_N_THREADS 85
#define MAX_DIMENSION 22

#define NULL_CHECK(var)                                                               \
    {                                                                                 \
        if (var == NULL) {                                                            \
            fprintf(stderr, "Error: Heap allocation returned null-pointer, exiting"); \
            exit(EXIT_FAILURE);                                                       \
        }                                                                             \
    }

#endif /* CONSTS_H */