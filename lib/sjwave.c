//
// Created by hsa on 07/01/17.
//

#include "sjwave.h"
#include "sjinc.h"

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

void sjfilter2d(float **a, int n2, int n1, char *mode) {
    int ix, iz;
    if (strcmp(mode, "laplace") == 0) {
        float **p = (float **) sjalloc2d(n2, n1, sizeof(float));
        memcpy(p[0], a[0], n2 * n1 * sizeof(float));
        for (ix = 2; ix < n2 - 2; ++ix)
            for (iz = 2; iz < n1 - 2; ++iz)
                a[ix][iz] = -4.0f * p[ix][iz] + p[ix - 1][iz] + p[ix + 1][iz] + p[ix][iz - 1] + p[ix][iz + 1];

        sjmfree2d(p);
    }
}

void sjsetsurface(float **a, int n2, int n1, float val) {
    int ix, iz;
    for (ix = 0; ix < n2; ++ix)
        for (iz = 0; iz < n1; ++iz)
            a[ix][iz] = val;
}

/**********************************************************************************************/
/* ! Finite Difference                                                                        */
/**********************************************************************************************/

/**********************************************************************************************/
/* ! Finite Difference                                                                        */
/**********************************************************************************************/

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

#define B60 (-2.982778e+0f)
#define B61 ( 1.714286e+0f)
#define B62 (-2.678571e-1f)
#define B63 ( 5.291005e-2f)
#define B64 (-8.928571e-3f)
#define B65 ( 1.038961e-3f)
#define B66 (-6.012506e-5f)

#define B611  0.562500000000
#define B612 -0.112500000000
#define B613  0.012500000000
#define B622  0.022500000000
#define B623 -0.002500000000
#define B633  0.000277777778

#define sjmfd2dn1(a, ix, iz)( B60* a[ix][iz]+ \
                            B61*(a[ix][iz+1]+a[ix][iz-1]) + \
                            B62*(a[ix][iz+2]+a[ix][iz-2]) + \
                            B63*(a[ix][iz+3]+a[ix][iz-3]) + \
                            B64*(a[ix][iz+4]+a[ix][iz-4]) + \
                            B65*(a[ix][iz+5]+a[ix][iz-5]) + \
                            B66*(a[ix][iz+6]+a[ix][iz-6]) )

#define sjmfd2dn2(a, ix, iz)( B60* a[ix][iz]+ \
                            B61*(a[ix+1][iz]+a[ix-1][iz]) + \
                            B62*(a[ix+2][iz]+a[ix-2][iz]) + \
                            B63*(a[ix+3][iz]+a[ix-3][iz]) + \
                            B64*(a[ix+4][iz]+a[ix-4][iz]) + \
                            B65*(a[ix+5][iz]+a[ix-5][iz]) + \
                            B66*(a[ix+6][iz]+a[ix-6][iz]) )

#define sjmfd2dnc(a, ix, iz)( B611*(a[ix+1][iz+1]-a[ix-1][iz+1]-a[ix+1][iz-1]+a[ix-1][iz-1]) + \
                            B612*(a[ix+1][iz+2]-a[ix-1][iz+2]-a[ix+1][iz-2]+a[ix-1][iz-2]+a[ix+2][iz+1]-a[ix-2][iz+1]-a[ix+2][iz-1]+a[ix-2][iz-1]) + \
                            B613*(a[ix+1][iz+3]-a[ix-1][iz+3]-a[ix+1][iz-3]+a[ix-1][iz-3]+a[ix+3][iz+1]-a[ix-3][iz+1]-a[ix+3][iz-1]+a[ix-3][iz-1]) + \
                            B622*(a[ix+2][iz+2]-a[ix-2][iz+2]-a[ix+2][iz-2]+a[ix-2][iz-2]) + \
                            B623*(a[ix+2][iz+3]-a[ix-2][iz+3]-a[ix+2][iz-3]+a[ix-2][iz-3]+a[ix+3][iz+2]-a[ix-3][iz+2]-a[ix+3][iz-2]+a[ix-3][iz-2]) + \
                            B633*(a[ix+3][iz+3]-a[ix-3][iz+3]-a[ix+3][iz-3]+a[ix-3][iz-3]))

//! Two dimension constant density acoustic forward simulation
void sjawfd2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    int srctrunc = opt->srctrunc;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    int ycutdirect = opt->ycutdirect;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = 1.0f / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {
        //! Source
        if (it < srctrunc)
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate veloctiy
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] =
                        cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) * dt2 + 2.0f * p1[ix][iz] -
                        p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wav->recz[ir][it] = p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wav->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
    }

    //------------------------ Cut direct wav ------------------------//
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];

        //------------------------ Initialization ------------------------//
        //! Boundary condition
        sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);

        //! Model
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids2;

        //! Wavefield
        memset(p2[0], 0, nxb * nzb * sizeof(float));
        memset(p1[0], 0, nxb * nzb * sizeof(float));
        memset(p0[0], 0, nxb * nzb * sizeof(float));

        //! Wavefield exploration
        for (it = 0; it < nt; it++) {

            //! Source
            if (it < srctrunc)
                for (ix = -srcrange; ix <= srcrange; ix++)
                    for (iz = -srcrange; iz <= srcrange; iz++)
                        p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

            //! Calculate veloctiy
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++)
                    p2[ix][iz] =
                            cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) * dt2 + 2.0f * p1[ix][iz] -
                            p0[ix][iz];

            //! Boundary condition
            sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++)
                wav->recz[ir][it] -= p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

            //! Update
            memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
            memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
}

//! Two dimension constant density acoustic scatter forward simulation
void sjaswfd2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    int srctrunc = opt->srctrunc;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = 1.0f / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **ipp = sjmflloc2d(nxb, nzb);
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);
    sjextend2d(geo->ipp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, ipp);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);
    float **s2 = sjmflloc2d(nxb, nzb);
    float **s1 = sjmflloc2d(nxb, nzb);
    float **s0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    memset(p2[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(s2[0], 0, nxb * nzb * sizeof(float));
    memset(s1[0], 0, nxb * nzb * sizeof(float));
    memset(s0[0], 0, nxb * nzb * sizeof(float));

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {
        //! Source
        if (it < srctrunc)
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Laplace operator
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz));

        //! Scatter source
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                s1[ix][iz] += p2[ix][iz] * ipp[ix][iz];

        //! Scatter wavefield
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                s2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(s1, ix, iz) + sjmfd2dn2(s1, ix, iz)) * dt2
                             + 2.0f * s1[ix][iz] - s0[ix][iz];

        //! Primary wavefield
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = p2[ix][iz] * dt2 + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(s2, s1, s0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wav->recz[ir][it] = s1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wav->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));

        memcpy(s0[0], s1[0], nxb * nzb * sizeof(float));
        memcpy(s1[0], s2[0], nxb * nzb * sizeof(float));
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);
    sjmfree2d(ipp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
    sjmfree2d(s2);
    sjmfree2d(s1);
    sjmfree2d(s0);
}

//! Two dimension acoustic reverse time simulation based on constant density equation
void sjawrtmfd2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //! Runtime
    int it, ir, ix, iz;
    const int marg = 6;

    //! Option
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //! Survey
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = 1.0f / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {
        //! Source#
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->recz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->recz[ir][it];
        }

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p0[ix][iz] =
                        cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) * dt2 + 2.0f * p1[ix][iz] -
                        p2[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p0, p1, p2, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0) {
            for (ix = nb + marg; ix < nxb - nb - marg; ix++) {
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    geo->ipp2d[ix - nb - marg][iz - nb - marg] +=
                            wav->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p1[ix][iz];
                    geo->nipp2d[ix - nb - marg][iz - nb - marg] +=
                            wav->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] *
                            wav->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg];
                }
            }
        }

        //! Update
        memcpy(p2[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p0[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p0);
    sjmfree2d(p1);
    sjmfree2d(p2);
}

