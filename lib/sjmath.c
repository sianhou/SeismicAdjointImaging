//
// Created by hsa on 07/12/16.
//

#include "sjmath.h"

#define B60 (-2.982778e+0f)
#define B61 ( 1.714286e+0f)
#define B62 (-2.678571e-1f)
#define B63 ( 5.291005e-2f)
#define B64 (-8.928571e-3f)
#define B65 ( 1.038961e-3f)
#define B66 (-6.012506e-5f)

#define SATOFDD2N1(a, ix, iz)(B60* a[ix][iz]+ \
                            B61*(a[ix][iz+1]+a[ix][iz-1]) + \
                            B62*(a[ix][iz+2]+a[ix][iz-2]) + \
                            B63*(a[ix][iz+3]+a[ix][iz-3]) + \
                            B64*(a[ix][iz+4]+a[ix][iz-4]) + \
                            B65*(a[ix][iz+5]+a[ix][iz-5]) + \
                            B66*(a[ix][iz+6]+a[ix][iz-6]) )

#define SATOFDD2N2(a, ix, iz)(B60* a[ix][iz]+ \
                            B61*(a[ix+1][iz]+a[ix-1][iz]) + \
                            B62*(a[ix+2][iz]+a[ix-2][iz]) + \
                            B63*(a[ix+3][iz]+a[ix-3][iz]) + \
                            B64*(a[ix+4][iz]+a[ix-4][iz]) + \
                            B65*(a[ix+5][iz]+a[ix-5][iz]) + \
                            B66*(a[ix+6][iz]+a[ix-6][iz]) )


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
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
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
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz, it)
#endif
    for (it = 0; it < nt; ++it)
        for (ix = 0; ix < nx; ++ix)
            for (iz = 0; iz < nz; ++iz)
                output[ix][iz] += input1[it][ix][iz] * input2[it][ix][iz];
}

void sjimagefilter2d(float **input, int n2, int n1, int mode) {
    int ix, iz;

    float **ptr = (float **) sjalloc2d(n2, n1, sizeof(float));

    memcpy(ptr[0], input[0], n2 * n1 * sizeof(float));

#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
    for (ix = 2; ix < n2 - 2; ++ix)
        for (iz = 2; iz < n1 - 2; ++iz)
            input[ix][iz] = -4.0f * ptr[ix][iz] + ptr[ix - 1][iz] + ptr[ix + 1][iz] + ptr[ix][iz - 1] + ptr[ix][iz + 1];
}

//! Process standard input
