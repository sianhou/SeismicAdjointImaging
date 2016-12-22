//
// Created by hsa on 07/12/16.
//

#ifndef GJI_SJTOOL_H
#define GJI_SJTOOL_H

#include "sjinc.h"

int sjricker1d(float *ricker, int nt, int t0, float dt, float fp, float amp);

void sjextend2d(float **input, int nx, int nz, int ex0, int ex1, int ez0, int ez1, float **output);

void sjextract2d(float **input, int x0, int z0, int nx, int nz, float **output);

void sjprojaddeq2d(float **input0, float **input1, int x0, int z0, int nx, int nz);

void sjprojdiveq2d(float **input0, float **input1, int x0, int z0, int nx, int nz);

void sjimage2d(float ***input1, float ***input2, int nt, int nx, int nz, int mode, float **output);

void sjimagefilter2d(float **input, int n2, int n1, int mode);

#endif //GJI_SJTOOL_H
