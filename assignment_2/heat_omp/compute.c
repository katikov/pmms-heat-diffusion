#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include <omp.h>
#include <stdio.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

extern int omp_heat_chunk_size;
extern int omp_heat_parallel_type;
extern int omp_heat_percentual;
extern int omp_heat_columns;
extern int omp_heat_pure_seq;

/// Replace/erase the following line:
#include "ref1.c"

void do_compute(const struct parameters* p, struct results *r)
{
/// Replace/erase the following line:
#include "ref2.c"
}