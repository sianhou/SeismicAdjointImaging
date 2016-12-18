//
// Created by hsa on 07/12/16.
//

#ifndef SJI_SJFD_H
#define SJI_SJFD_H

#include "sjinc.h"

//! Two dimension acoustic simulation based on constant velocity-stress equation
void sjawsgfd2d(int nt, int sx, int sz, int srcrange, int srctrunc, //! Source
                float dt, float srcdecay, float *wav,
                int nx, int nz, //! Model
                float ds, float **vp,
                int nb, //! Boundary condition
                int nr, int *rx, int *rz, //! Survey
                int ysnap, int jsnap,
                float **record, float ***snap, //! Wavefield
                int ycutdirect);

//! Two dimension acoustic reverse time simulation based on constant velocity-stress equation
void sjawrtsgfd2d(int nt, float dt, //! Source
                  int nx, int nz, //! Model
                  float ds, float **vp,
                  int nb, //! Boundary condition
                  int nr, //! Survey
                  int *rx, int *rz,
                  int ysnap, int jsnap, //! Wavefield
                  float **rec, float ***snap);

#endif //SJI_SJFD_H
