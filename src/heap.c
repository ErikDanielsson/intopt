#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bnb.h"
#include "consts.h"
#include "heap.h"

struct heap_elem *h = NULL;
unsigned long heap_size = 0;
int heap_idx = 0;

static inline void sift_up(double z, node_t *p)
{
    int i = heap_idx;
    double current_z = h[(i - 1) / 2].z;
    while (i != 0 && z > current_z) {
        h[i] = h[(i - 1) / 2];
        i = (i - 1) / 2;
        current_z = h[(i - 1) / 2].z;
    }
    h[i] = (struct heap_elem){z, p};
}

void add_to_set(node_t *p)
{
    pthread_mutex_lock(&best_yet_z_mutex);
    if (p->z < best_yet_z) {
        pthread_mutex_unlock(&best_yet_z_mutex);
        return;
    }
    pthread_mutex_unlock(&best_yet_z_mutex);
    double z = p->z;
    if (heap_size == heap_idx) {
        heap_size *= 2;
        h = realloc(h, heap_size * sizeof(struct heap_elem));
        NULL_CHECK(h);
    }

    sift_up(z, p);
    heap_idx++;
}

bool set_nonempty()
{
    return heap_idx > 0;
}

static inline void heapify(int i)
{
    int s = i, largest = i;
    while (s < heap_idx) {
        int l = 2 * s + 1;
        int r = l + 1;
        largest = s;
        if (l < heap_idx && h[l].z > h[largest].z)
            largest = l;
        if (r < heap_idx && h[r].z > h[largest].z)
            largest = r;
        if (largest != s) {
            struct heap_elem tmp = h[s];
            h[s] = h[largest];
            h[largest] = tmp;
            s = largest;
        } else {
            break;
        }
    }
}

node_t *next_from_set()
{
    node_t *next_n = h[0].node;
    heap_idx--;
    h[0] = h[heap_idx];
    h[heap_idx] = (struct heap_elem){0, NULL};

    heapify(0);
    return next_n;
}

void bound_set(double z)
{
    int i, j, k;

    int new_idx = heap_idx;
    k = heap_idx;
    j = heap_idx;

    bool last_is_empty = false;
    for (i = heap_idx - 1; i >= 0; i--) {
        if (h[i].z < z) {
            remove_node(h[i].node);
            h[i] = (struct heap_elem){0, NULL};
            new_idx--;
            last_is_empty = true;
        } else {
            if (!last_is_empty) {
                j--;
            } else {
                int mv_dist = j - (i + 1);
                int mem_dist = k - j;
                if (mv_dist >= mem_dist) {
                    memcpy(&h[i + 1], &h[j], mem_dist * sizeof(struct heap_elem));
                } else {
                    memcpy(&h[i + 1], &h[j + mem_dist - mv_dist], mv_dist * sizeof(struct heap_elem));
                }
                k = i + mem_dist + 1;
                j = i;
            }
            last_is_empty = false;
        }
    }
    heap_idx = new_idx;

    for (i = heap_idx / 2; i >= 0; i--) {
        heapify(i);
    }
}