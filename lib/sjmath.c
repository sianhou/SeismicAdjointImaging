//
// Created by hsa on 07/12/16.
//

#include "sjmath.h"

int sjricker1d(float *ricker, int nt, int t0, float dt, float fp, float amp) {

    int it = 0;

    float tmpf = (float) (pi * pi * fp * fp), tmpt = dt * dt;

    double tmpr;

    for (it = 0; it < nt; ++it) {
        tmpr = (1.0 - 2.0 * tmpf * (it - t0) * (it - t0) * tmpt) * exp(-tmpf * (it - t0) * (it - t0) * tmpt);
        ricker[it] = ((float) tmpr) * amp;
    }

    return 1;
}

void sjextend2d(float **input, int nx, int nz,
                int ex0, int ex1, int ez0, int ez1, float **output) {
    int ix, iz;

    //! Centering
    for (ix = 0; ix < nx; ++ix)
        for (iz = 0; iz < nz; ++iz)
            output[ix + ex0][iz + ez0] = input[ix][iz];

    //! Left
    for (ix = 0; ix < ex0; ++ix)
        for (iz = 0; iz < nz; ++iz)
            output[ix][iz + ez0] = input[0][iz];

    //! Right
    for (ix = 0; ix < ex1; ++ix)
        for (iz = 0; iz < nz; ++iz)
            output[ex0 + nx + ix][iz + ez0] = input[nx - 1][iz];

    //! Upper
    for (ix = 0; ix < ex0 + nx + ex1; ++ix)
        for (iz = 0; iz < ez0; ++iz)
            output[ix][iz] = output[ix][ez0];

    //! Below
    for (ix = 0; ix < ex0 + nx + ex1; ++ix)
        for (iz = 0; iz < ez1; ++iz)
            output[ix][ez0 + nz + iz] = output[ix][ez0 + nz - 1];
}

void sjextract2d(float **input, int x0, int z0, int nx, int nz, float **output) {
    int ix, iz;
    for (ix = 0; ix < nx; ++ix)
        for (iz = 0; iz < nz; ++iz)
            output[ix][iz] = input[x0 + ix][z0 + iz];
}

void sjprojaddeq2d(float **input0, float **input1, int x0, int z0, int nx, int nz) {
    int ix, iz;
    for (ix = 0; ix < nx; ++ix)
        for (iz = 0; iz < nz; ++iz)
            input0[x0 + ix][z0 + iz] += input1[ix][iz];
}

void sjprojdiveq2d(float **input0, float **input1, int x0, int z0, int nx, int nz) {
    int ix, iz;
    for (ix = 0; ix < nx; ++ix)
        for (iz = 0; iz < nz; ++iz)
            input0[x0 + ix][z0 + iz] /= (input1[ix][iz] + 1.0e-4);
}

void sjimage2d(float ***input1, float ***input2, int nt, int nx, int nz, int mode, float **output) {
    int ix, iz, it;
    memset(output[0], 0, nx * nz * sizeof(float));
    for (it = 0; it < nt; ++it)
        for (ix = 0; ix < nx; ++ix)
            for (iz = 0; iz < nz; ++iz)
                output[ix][iz] += input1[it][ix][iz] * input2[it][ix][iz];
}

void sjlaplcefilter2d(float **input, int n2, int n1) {
    int ix, iz;

    float **ptr = (float **) sjalloc2d(n2, n1, sizeof(float));

    memcpy(ptr[0], input[0], n2 * n1 * sizeof(float));

    for (ix = 2; ix < n2 - 2; ++ix)
        for (iz = 2; iz < n1 - 2; ++iz)
            input[ix][iz] = -4.0f * ptr[ix][iz] + ptr[ix - 1][iz] + ptr[ix + 1][iz] + ptr[ix][iz - 1] + ptr[ix][iz + 1];
}

//! Process standard input
