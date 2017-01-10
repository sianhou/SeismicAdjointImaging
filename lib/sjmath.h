//
// Created by hsa on 07/12/16.
//

#ifndef GJI_SJMATH_H
#define GJI_SJMATH_H

#include "sjinc.h"

void sjguasssmoothf2d(float **z, int n2, int n1, float alpha, int length, float **x);

float sjvecdotf(int n, float a, float *x, float *y);

void sjvecaddf(float *z, int n, float a, float *x, float b, float *y); //! c = a + b

void sjvecsubf(float *z, int n, float a, float *x, float b, float *y); //! c = a - b

//! z[] = a * x[] * y[]
void sjvecmulf(float *z, int n, float a, float *y, float *x);

void sjvecdivf(float *z, int n, float a, float *x, float *y, float ep); //! c = a / (b + ep)

//! Z[] = 0.0f
void sjveczerof(float *z, int n);

void sjcgdirection(float *d, int n, float *g1, float *g0, int iter);

float sjcgstepsize(int n, float *s, float *x, float err, int iter);

#endif //GJI_SJMATH_H
