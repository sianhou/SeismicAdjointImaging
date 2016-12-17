//
// Created by hsa on 07/12/16.
//

#ifndef SJI_SJINC_H
#define SJI_SJINC_H

#include <sys/malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef GFDOPENMP_
#include <omp.h>
#endif

#include "sjabc.h"
#include "sjfile.h"
#include "sjmalloc.h"
#include "sjmath.h"
#include "sjfd.h"

//! complex
typedef struct fcomplex_ {
    float r, i;
} fcomplex;

typedef struct dcomplex_ {
    double r, i;
} dcomplex;

#define pi 3.141592653589793

#endif //SJI_SJINC_H
