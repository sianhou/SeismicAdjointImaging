//
// Created by hsa on 07/12/16.
//

#ifndef SJI_SJINC_H
#define SJI_SJINC_H

#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef GFDOPENMP_
#include <omp.h>
#endif

//! complex
typedef struct fcomplex_ {
    float r, i;
} fcomplex;

typedef struct dcomplex_ {
    double r, i;
} dcomplex;

//! survey
typedef struct {
    int sy, sx, sz;     //! source position in y(n3)- x(n2)- and z(n1)-direction
    int y0, x0, z0;     //! local model start in y(n3)- x(n2)- and z(n1)-direction
    int ny, nx, nz;     //! local model size in y(n3)- x(n2)- and z(n1)-direction
    int gnx, gny, gnz;  //! global model size in y(n3)- x(n2)- and z(n1)-direction
    int nr, tr;

    int ns, maxnr;
    int nparas;
    int *ry, *rx, *rz;
    char *surveyfile;
} sjssurvey;

//! source
typedef struct {
    int k1, srcrange, srctrunc;

    float dt, fp, amp, srcdecay;
    float *wavelet;
} sjssource;

//! model
typedef struct {
    int nb;
    float ds;
    float **gvp2d, **gvs2d;
    float **vp2d, **vs2d;
    float **gipp2d;
    float **ipp2d, **nipp2d;
    char *vpfile, *vsfile, *ippfile;
} sjsgeo;

//! wave
typedef struct {
    int nt, jsnap, nsnap, ycutdirect;
    float **recy, **recx, **recz;
    float ***snapy2d, ***snapx2d, ***snapz2d;
    char *recyfile, *recxfile, *reczfile;
    int yadjointbc; //! == 1: process adjoint wavefield as boundary condition
                    //! != 1: process adjoint wavefield as source
} sjswave;

#define pi 3.141592653589793

#include "sjsimulation.h"

#include "sjabc.h"
#include "sjfile.h"
#include "sjmalloc.h"
#include "sjmath.h"


#endif //SJI_SJINC_H
