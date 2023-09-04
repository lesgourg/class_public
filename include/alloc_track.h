// Header guards
#ifndef ALLOC_TRACK_H
#define ALLOC_TRACK_H

#include <stdlib.h>

#ifdef _OPENMP
#include "omp.h"
#endif

// Enable if you want to standalone test
// #define DEBUG_ALLOC_TRACK

void * tracked_malloc(size_t size);
void * tracked_calloc(size_t nmemb, size_t size);
void tracked_free(void *ptr);
void * tracked_realloc(void *ptr, size_t size);
void tracked_free_all(void);
void walk_chunks(void);
int init_memory_tracking(void);
void shutdown_memory_tracking(void);

#endif
