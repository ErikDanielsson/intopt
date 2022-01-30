#ifndef HEAP_H
#define HEAP_H

#include <stdbool.h>

#include "bnb.h"

struct heap_elem
{
    double z;
    node_t *node;
};
extern struct heap_elem *h;
extern unsigned long heap_size;
extern int heap_idx;

void add_to_set(node_t *p);
bool set_nonempty();
node_t *next_from_set();
void bound_set(double z);

#endif /* HEAP_H */