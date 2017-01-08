//
// Created by hsa on 07/12/16.
//

#include "sjmath.h"

void sjvecaddf(float *z, int n, float a, float *x, float b, float *y) {
    //! z[] = a*x[] + b*y[]
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] + b * y[ii];
}

void sjvecsubf(float *z, int n, float a, float *x, float b, float *y) {
    //! z[] = a*x[] - b*y[]
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] - b * y[ii];
}

void sjvecmulf(float *z, int n, float a, float *y, float *x) {
    //! z[] = a * x[] * y[]
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] * y[ii];
}

void sjvecdivf(float *z, int n, float a, float *x, float *y, float ep) {
    //! z[] = a * x[] / (b[] + ep)
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] / (y[ii] + ep);
}

void sjveczerof(float *z, int n) {
    memset(z, 0, n * sizeof(float));
}

//! Process standard input
