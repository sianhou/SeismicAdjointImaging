//
// Created by hsa on 07/12/16.
//

#include "sjinc.h"

#define C50 ( 1.239407e+0f)
#define C51 (-1.105315e-1f)
#define C52 ( 2.496329e-2f)
#define C53 (-5.804879e-3f)
#define C54 ( 9.358680e-4f)

#define sjmsgfd2dn1(p, ix, iz) (C50*(p[ix][iz+1]-p[ix][iz-0]) + \
                                C51*(p[ix][iz+2]-p[ix][iz-1]) + \
                                C52*(p[ix][iz+3]-p[ix][iz-2]) + \
                                C53*(p[ix][iz+4]-p[ix][iz-3]) + \
                                C54*(p[ix][iz+5]-p[ix][iz-4]) )

#define sjmsgfd2dn2(p, ix, iz) (C50*(p[ix+1][iz]-p[ix-0][iz]) + \
                                C51*(p[ix+2][iz]-p[ix-1][iz]) + \
                                C52*(p[ix+3][iz]-p[ix-2][iz]) + \
                                C53*(p[ix+4][iz]-p[ix-3][iz]) + \
                                C54*(p[ix+5][iz]-p[ix-4][iz]) )

//! Two dimension acoustic simulation based on constant velocity-stress equation
void sjawsgfd2d(int nt, int sx, int sz, int srcrange, int srctrunc, //! Source
                float dt, float srcdecay, float *wav,
                int nx, int nz, //! Model
                float ds, float **vp,
                int nb, //! Boundary condition
                int nr, //! Survey
                int *rx, int *rz,
                int ysnap, int jsnap, //! Wavefield
                float **record, float ***snap) {

    //------------------------ Main loop ------------------------//
    //! Define parameters
    int it, ir, ix, iz;

    //------------------------ Finite difference ------------------------//
    //! Define parameters
    const int marg = 6;

    //------------------------ Source ------------------------//
    //! Calculate parameters
    sx += nb + marg;
    sz += nb + marg;

    //------------------------ Model ------------------------//
    //! Define parameters
    int nxb, nzb;
    float ids;
    float **cp;
    //! Calculate parameters
    ids = -dt / ds;
    nxb = nx + 2 * marg + 2 * nb;
    nzb = nz + 2 * marg + 2 * nb;
    //! Allocate memory
    cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    //! Extend the model
    sjextend2d(vp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //------------------------ Boundary condition ------------------------//
    //! Define parameters
    float **gxl, **gxr, **gzu, **gzb;
    //! Allocate memory
    gxl = (float **) sjalloc2d(nzb, 3, sizeof(float));
    gxr = (float **) sjalloc2d(nzb, 3, sizeof(float));
    gzu = (float **) sjalloc2d(nxb, 3, sizeof(float));
    gzb = (float **) sjalloc2d(nxb, 3, sizeof(float));

    //------------------------ Wavefield ------------------------//
    //! Define parameters
    float **vx0, **vx1, **vz0, **vz1, **p0, **p1;
    //! Allocate memory
    vx0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    vx1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    vz0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    vz1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    p0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    p1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitohabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    for (ix = 0; ix < nxb; ix++)
        for (iz = 0; iz < nzb; iz++)
            cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids;
    //! Wavefield
    memset(vx0[0], 0, nxb * nzb * sizeof(float));
    memset(vx1[0], 0, nxb * nzb * sizeof(float));
    memset(vz0[0], 0, nxb * nzb * sizeof(float));
    memset(vz1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));

    /**********************************************************************************************/
    /* ! Wavefield                                                                                */
    /**********************************************************************************************/

    //! Wavefield exploration
    for (it = 0; it < nt; it++) {
        //! Source
        if (it < srctrunc)
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p0[ix + sx][iz + sz] += wav[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate veloctiy
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
        for (ix = marg; ix < nxb - marg; ix++) {
            for (iz = marg; iz < nzb - marg; iz++) {
                vx1[ix][iz] = sjmsgfd2dn2(p0, ix, iz) * ids + vx0[ix][iz];
                vz1[ix][iz] = sjmsgfd2dn1(p0, ix, iz) * ids + vz0[ix][iz];
            }
        }

        //! Calculate stress
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p1[ix][iz] = cp[ix][iz] * (sjmsgfd2dn2(vx1, ix - 1, iz) + sjmsgfd2dn1(vz1, ix, iz - 1)) + p0[ix][iz];

        //! Boundary condition
        sjapplyohabc2d(vx1, vx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(vz1, vz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
#ifdef GFDOPENMP_
#pragma omp parallel for private(ir)
#endif
        for (ir = 0; ir < nr; ir++)
            record[ir][it] = p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if (ysnap == 1 && ((it % jsnap) == 0))
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    snap[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(vx0[0], vx1[0], nxb * nzb * sizeof(float));
        memcpy(vz0[0], vz1[0], nxb * nzb * sizeof(float));
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(vx0);
    sjmfree2d(vx1);
    sjmfree2d(vz0);
    sjmfree2d(vz1);
    sjmfree2d(p0);
    sjmfree2d(p1);
}

//! Two dimension acoustic reverse time simulation based on constant velocity-stress equation
void sjawrtsgfd2d(int nt, float dt, //! Source
                  int nx, int nz, //! Model
                  float ds, float **vp,
                  int nb, //! Boundary condition
                  int nr, //! Survey
                  int *rx, int *rz,
                  int ysnap, int jsnap, //! Wavefield
                  float **rec, float ***snap) {

    //------------------------ Runtime ------------------------//
    //! Define parameters
    int it, ir, ix, iz;

    //------------------------ Finite difference ------------------------//
    //! Define parameters
    const int marg = 6;

    //------------------------ Model ------------------------//
    //! Define parameters
    int nxb, nzb;
    float ids;
    float **cp;
    //! Calculate parameters
    ids = -dt / ds;
    nxb = nx + 2 * marg + 2 * nb;
    nzb = nz + 2 * marg + 2 * nb;
    //! Allocate memory
    cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    //! Extend the model
    sjextend2d(vp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //------------------------ Boundary condition ------------------------//
    //! Define parameters
    float **gxl, **gxr, **gzu, **gzb;
    //! Allocate memory
    gxl = (float **) sjalloc2d(nzb, 3, sizeof(float));
    gxr = (float **) sjalloc2d(nzb, 3, sizeof(float));
    gzu = (float **) sjalloc2d(nxb, 3, sizeof(float));
    gzb = (float **) sjalloc2d(nxb, 3, sizeof(float));

    //------------------------ Wavefield ------------------------//
    //! Define parameters
    float **vx0, **vx1, **vz0, **vz1, **p0, **p1;
    //! Allocate memory
    vx0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    vx1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    vz0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    vz1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    p0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    p1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitohabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    for (ix = 0; ix < nxb; ix++)
        for (iz = 0; iz < nzb; iz++)
            cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids;
    //! Wavefield
    memset(vx0[0], 0, nxb * nzb * sizeof(float));
    memset(vx1[0], 0, nxb * nzb * sizeof(float));
    memset(vz0[0], 0, nxb * nzb * sizeof(float));
    memset(vz1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));

    /**********************************************************************************************/
    /**/ //! Wavefield                                                                          /**/
    /**********************************************************************************************/

    //! Wavefield revese time exploration
    for (it = nt - 1; it >= 0; --it) {
        //! Source
#ifdef GFDOPENMP_
#pragma omp parallel for private(ir)
#endif
        for (ir = 0; ir < nr; ir++)
            p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = rec[ir][it];

        //! Calculate velocity
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
        for (ix = marg; ix < nxb - marg; ix++) {
            for (iz = marg; iz < nzb - marg; iz++) {
                vx0[ix][iz] = sjmsgfd2dn2(p1, ix, iz) * ids + vx1[ix][iz];
                vz0[ix][iz] = sjmsgfd2dn1(p1, ix, iz) * ids + vz1[ix][iz];
            }
        }

        //! Calculate stress
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p0[ix][iz] = cp[ix][iz] * (sjmsgfd2dn2(vx0, ix - 1, iz) + sjmsgfd2dn1(vz0, ix, iz - 1)) + p0[ix][iz];

        //! Boundary condition
        sjapplyohabc2d(vx0, vx1, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(vz0, vz1, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(p0, p1, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if (ysnap == 1 && ((it % jsnap) == 0))
#ifdef GFDOPENMP_
#pragma omp parallel for private(ix, iz)
#endif
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    snap[it / jsnap][ix - nb - marg][iz - nb - marg] = p0[ix][iz];

        //! Update
        memcpy(&p1[0][0], &p0[0][0], nxb * nzb * sizeof(float));
        memcpy(&vx1[0][0], &vx0[0][0], nxb * nzb * sizeof(float));
        memcpy(&vz1[0][0], &vz0[0][0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(vx0);
    sjmfree2d(vx1);
    sjmfree2d(vz0);
    sjmfree2d(vz1);
    sjmfree2d(p0);
    sjmfree2d(p1);
}