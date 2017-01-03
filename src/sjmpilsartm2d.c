//
// Created by 侯思安 on 2016/12/28.
//

#include "../lib/sjinc.h"
#include <mpi.h>

void sjascfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave);

void sjawrtfd2dx(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave);

int main(int argc, char *argv[]) {

    //! Runtime
    int is = 0, iter = 0, flag = 1, ix, iz;
    double tstart, tend, Tstart, Tend;

    //! MPI
    int rankid, nrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //! Source
    sjssource source;
    flag &= sjssource_init(&source);
    flag &= sjssource_getparas(&source, argc, argv);

    //! Survey
    sjssurvey survey;
    flag &= sjssurvey_init(&survey);
    flag &= sjssurvey_getparas(&survey, argc, argv);

    //! Model
    sjsgeo geo;
    flag &= sjsgeo_init(&geo);
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "vp");
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "ipp");
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "lsipp");

    //! Wave
    sjswave wave;
    flag &= sjswave_init(&wave);
    flag &= sjswave_getparas(&wave, argc, argv, "recz");

    if (flag) {

        //! Source
        source.wavelet = (float *) sjalloc1d(source.srctrunc, sizeof(float));
        sjricker1d(source.wavelet, source.srctrunc, source.k1, source.dt, source.fp, source.amp);

        //! Model
        geo.gvp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        geo.gipp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        sjreadsuall(geo.gvp2d[0], survey.gnx, survey.gnz, geo.vpfile);
        sjreadsuall(geo.gipp2d[0], survey.gnx, survey.gnz, geo.ippfile);
        float **gipp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));

        //! ------------------------ Lsrtm ------------------------
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("---------------- 2D Acoustic LSRTM start ----------------\n");
        }

        for (iter = 0; iter < 10; ++iter) {

            //! Bcast model
            MPI_Bcast(geo.gvp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Bcast(geo.gipp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);


            for (is = rankid; is < survey.ns; is += nrank) {
                //! Time
                tstart = (double) clock();

                //! Survey
                sjssurvey_readis(&survey, is);

                //! Model
                geo.vp2d = (float **) sjalloc2d(survey.nx, survey.nz, sizeof(float));
                geo.ipp2d = (float **) sjalloc2d(survey.nx, survey.nz, sizeof(float));
                sjextract2d(geo.gvp2d, survey.x0, survey.z0, survey.nx, survey.nz, geo.vp2d);
                sjextract2d(geo.gipp2d, survey.x0, survey.z0, survey.nx, survey.nz, geo.ipp2d);

                //! Forward simulation
                wave.recz = (float **) sjalloc2d(survey.nr, wave.nt, sizeof(float));
                wave.snapz2d = (float ***) sjalloc3d(wave.nsnap, survey.nx, survey.nz, sizeof(float));
                sjascfd2d(&source, &survey, &geo, &wave);

                //! Difference wavefield
                float **recz = (float **) sjalloc2d(survey.nr, wave.nt, sizeof(float));
                sjreadsu(recz[0], survey.nr, wave.nt, sizeof(float), survey.tr, 0, wave.reczfile);
                for (ix = 0; ix < survey.nr; ++ix)
                    for (iz = 0; iz < wave.nt; ++iz)
                        wave.recz[ix][iz] -= recz[ix][iz];

                //! Backward simulation
                memset(geo.ipp2d[0], 0, survey.nx * survey.nz * sizeof(float));
                sjawrtfd2dx(&source, &survey, &geo, &wave);

                //! Stacking
                sjprojaddeq2d(gipp2d, geo.ipp2d, survey.x0, survey.z0, survey.nx, survey.nz);

                //! Free
                sjmcheckfree2d(geo.vp2d);
                sjmcheckfree2d(geo.ipp2d);
                sjmcheckfree2d(wave.recz);
                sjmcheckfree2d(recz);
                sjmcheckfree3d(wave.snapz2d);

                //! Time
                tend = (double) clock();
                printf("Single shot LSARTM complete - %d/%d - time=%fs.\n", is + 1, survey.ns,
                       (tend - tstart) / CLOCKS_PER_SEC);
                printf("Rankid=%d, iter=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n", rankid, iter + 1,
                       survey.sx + survey.x0, survey.sz + survey.z0,
                       survey.rx[0] + survey.x0, survey.rx[survey.nr - 1] + survey.x0,
                       survey.rz[0] + survey.z0, survey.rz[survey.nr - 1] + survey.z0);
            }

            //! Communication
            if (rankid == 0) {
                MPI_Reduce(MPI_IN_PLACE, gipp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(gipp2d[0], gipp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            //! Add
            if (rankid == 0)
                for (ix = 0; ix < survey.gnx; ++ix)
                    for (iz = 0; iz < survey.gnz; ++iz)
                        geo.gipp2d[ix][iz] += gipp2d[ix][iz];
        }
        
        if (rankid == 0) {
            //! Output
            sjwritesuall(geo.gipp2d[0], survey.gnx, survey.gnz, geo.ds, geo.lsippfile);

            //! Time
            Tend = (double) clock();
            printf("Acoustic LSRTM complete - time=%fs.\n", (Tend - Tstart) / CLOCKS_PER_SEC);
        }

    } else {
        printf("\nExamples:   sjmpilsartm2d survey=survey.su vp=vp.su recz=recz.su ipp=mig.su lsipp=lsmig.su\n");
        sjbasicinformation();
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}

//! Two dimension acoustic simulation based on constant density equation
void sjascfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

    //! Runtime
    int it, ir, ix, iz;

    //! Finite difference
    const int marg = 6;

    //! Source
    int srctrunc = source->srctrunc;
    int srcrange = source->srcrange;
    float dt = source->dt;
    float srcdecay = source->srcdecay;
    float *wav = source->wavelet;

    //! Survey
    int nx = survey->nx;
    int nz = survey->nz;
    int sx = survey->sx + geo->nb + marg;
    int sz = survey->sz + geo->nb + marg;
    int nr = survey->nr;
    int *rx = survey->rx;
    int *rz = survey->rz;

    //! Model
    int nb = geo->nb;
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ds = geo->ds;
    float ids = source->dt * source->dt / ds / ds;
    float **cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **ipp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);
    sjextend2d(geo->ipp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, ipp);

    //! Wave
    int nt = wave->nt;
    int jsnap = wave->jsnap;

    //! Boundary condition
    float **gxl = (float **) sjalloc2d(nzb, 8, sizeof(float));
    float **gxr = (float **) sjalloc2d(nzb, 8, sizeof(float));
    float **gzu = (float **) sjalloc2d(nxb, 8, sizeof(float));
    float **gzb = (float **) sjalloc2d(nxb, 8, sizeof(float));

    //! Wavefield
    float **p2 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));

    float **s2 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **s1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **s0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    for (ix = 0; ix < nxb; ix++)
        for (iz = 0; iz < nzb; iz++)
            cp[ix][iz] = cp[ix][iz] * cp[ix][iz];
    //! Wavefield
    memset(p2[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));

    memset(s2[0], 0, nxb * nzb * sizeof(float));
    memset(s1[0], 0, nxb * nzb * sizeof(float));
    memset(s0[0], 0, nxb * nzb * sizeof(float));

    /**********************************************************************************************/
    /* ! Wavefield                                                                                */
    /**********************************************************************************************/

    //! Wavefield exploration
    for (it = 0; it < nt; it++) {

        //! Source
        if (it < srctrunc)
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p1[ix + sx][iz + sz] += wav[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate laplace operator
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) * ids;

        //! Calculate scatter wavefield
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                s2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(s1, ix, iz) + sjmfd2dn2(s1, ix, iz)) * ids
                             + 2.0f * s1[ix][iz] - s0[ix][iz] + p2[ix][iz] * cp[ix][iz] * ipp[ix][iz];

        //! Calculate
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = p2[ix][iz] + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(s2, s1, s0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wave->recz[ir][it] = s1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));

        memcpy(s0[0], s1[0], nxb * nzb * sizeof(float));
        memcpy(s1[0], s2[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

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
void sjawrtfd2dx(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

    //! Runtime
    int it, ir, ix, iz;

    //! Finite difference
    const int marg = 6;

    //! Source
    float dt = source->dt;

    //! Survey
    int nx = survey->nx;
    int nz = survey->nz;
    int nr = survey->nr;
    int *rx = survey->rx;
    int *rz = survey->rz;

    //! Model
    int nb = geo->nb;
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ds = geo->ds;
    float ids = dt * dt / ds / ds;
    float **cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //! Wave
    int nt = wave->nt;
    int jsnap = wave->jsnap;

    //! Boundary condition
    float **gxl = (float **) sjalloc2d(nzb, 8, sizeof(float));
    float **gxr = (float **) sjalloc2d(nzb, 8, sizeof(float));
    float **gzu = (float **) sjalloc2d(nxb, 8, sizeof(float));
    float **gzb = (float **) sjalloc2d(nxb, 8, sizeof(float));

    //! Wavefield
    float **p2 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    for (ix = 0; ix < nxb; ix++)
        for (iz = 0; iz < nzb; iz++)
            cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids;
    //! Wavefield
    memset(p2[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));

    /**********************************************************************************************/
    /* ! Wavefield                                                                                */
    /**********************************************************************************************/

    //! Wavefield reverse time exploration
    for (it = nt - 1; it >= 0; --it) {

        for (ir = 0; ir < nr; ir++)
            p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wave->recz[ir][it];

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p0[ix][iz] =
                        cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) + 2.0f * p1[ix][iz] - p2[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p0, p1, p2, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    geo->ipp2d[ix - nb - marg][iz - nb - marg] +=
                            wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] *
                            (p0[ix][iz] - 2.0f * p1[ix][iz] + p2[ix][iz]) * dt * dt;

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