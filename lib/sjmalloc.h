//
// Created by hsa on 07/12/16.
//

#ifndef SJI_SJMALLOC_H
#define SJI_SJMALLOC_H

#include "sjinc.h"

//! function for allocate the 1 demension array and initialize with zero
void *sjalloc1d(int n1, int size);

//! function for allocate the 2 demension array and initialize with zero
void **sjalloc2d(int n2, int n1, int size);

//! function for allocate the 3 demension array and initialize with zero
void ***sjalloc3d(int n3, int n2, int n1, int size);

//! function for free the 1 demension array
void sjfree1d(void *p);

void sjcheckfree1d(void *p);

//! function for free the 2 demension array
void sjfree2d(void **p);

void sjcheckfree2d(void **p);

//! function for free the 3 demension array
void sjfree3d(void ***p);

void sjcheckfree3d(void ***p);

//! Fast allocate the 1 demension using macro
#define sjmilloc1d(p, n1) (int *)sjalloc1d(n1, sizeof(int))
#define sjmflloc1d(p, n1) (float *)sjalloc1d(n1, sizeof(float))
#define sjmdlloc1d(p, n1) (double *)sjalloc1d(n1, sizeof(double))
#define sjmclloc1d(p, n1) (fcomplex *)sjalloc1d(n1, sizeof(fcomplex))
#define sjmzlloc1d(p, n1) (dcomplex *)sjalloc1d(n1, sizeof(dcomplex))

//! Fast allocate the 2 demension using macro
#define sjmilloc2d(p, n2, n1) (int **)sjalloc2d(n2, n1,  sizeof(int))
#define sjmflloc2d(p, n2, n1) (float **)sjalloc2d(n2, n1, sizeof(float))
#define sjmdlloc2d(p, n2, n1) (double **)sjalloc2d(n2, n1, sizeof(double))
#define sjmclloc2d(p, n2, n1) (fcomplex **)sjalloc2d(n2, n1, sizeof(fcomplex))
#define sjmzlloc2d(p, n2, n1) (dcomplex **)sjalloc2d(n2, n1, sizeof(dcomplex))

//! Fast allocate the 3 demension using macro
#define sjmilloc3d(p, n3, n2, n1) (int ***)sjalloc3d(n3, n2, n1, sizeof(int))
#define sjmflloc3d(p, n3, n2, n1) (float ***)sjalloc3d(n3, n2, n1, sizeof(float))
#define sjmdlloc3d(p, n3, n2, n1) (double ***)sjalloc3d(n3, n2, n1, sizeof(double))
#define sjmclloc3d(p, n3, n2, n1) (fcomplex ***)sjalloc3d(n3, n2, n1, sizeof(fcomplex))
#define sjmzlloc3d(p, n3, n2, n1) (dcomplex ***)sjalloc3d(n3, n2, n1, sizeof(dcomplex))

//! Fast free the 1 demension using macro
#define sjmfree1d(p) sjfree1d((void *) p)
#define sjmfree2d(p) sjfree2d((void **) p)
#define sjmfree3d(p) sjfree3d((void ***) p)

#define sjmcheckfree1d(p) sjcheckfree1d((void *) p); p=NULL;
#define sjmcheckfree2d(p) sjcheckfree2d((void **) p); p=NULL;
#define sjmcheckfree3d(p) sjcheckfree3d((void ***) p); p=NULL;

#endif //SJI_SJMALLOC_H
