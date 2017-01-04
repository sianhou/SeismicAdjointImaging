//
// Created by hsa on 12/12/16.
//
#include <mpi.h>
#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    //! Runtime
    int is = 0, ix, iz, flag = 1;
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

    //! Wave
    sjswave wave;
    flag &= sjswave_init(&wave);
    flag &= sjswave_getparas(&wave, argc, argv, "recz");

    //! ------------------------ RTM2D ------------------------
    if (flag) {

        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic RTM start ------------------------\n");
        }

        //! Source
        source.wavelet = (float *) sjalloc1d(source.srctrunc, sizeof(float));
        sjricker1d(source.wavelet, source.srctrunc, source.k1, source.dt, source.fp, source.amp);

        //! Model
        float **nmig = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        geo.gipp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        geo.gvp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        sjreadsuall(geo.gvp2d[0], survey.gnx, survey.gnz, geo.vpfile);
        MPI_Bcast(geo.gvp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

        //! Rtm
        for (is = rankid; is < survey.ns; is += nrank) {
            //! Time
            tstart = (double) clock();

            //! Survey
            sjssurvey_readis(&survey, is);

            //! Model
            geo.vp2d = (float **) sjalloc2d(survey.nx, survey.nz, sizeof(float));
            geo.ipp2d = (float **) sjalloc2d(survey.nx, survey.nz, sizeof(float));
            geo.nipp2d = (float **) sjalloc2d(survey.nx, survey.nz, sizeof(float));
            sjextract2d(geo.gvp2d, survey.x0, survey.z0, survey.nx, survey.nz, geo.vp2d);

            //! Wavefield
            wave.recz = (float **) sjalloc2d(survey.nr, wave.nt, sizeof(float));
            wave.snapz2d = (float ***) sjalloc3d(wave.nsnap, survey.nx, survey.nz, sizeof(float));

            //! Forward simulaion
            sjawfd2d(&source, &survey, &geo, &wave);

            //! Read record
            sjreadsu(wave.recz[0], survey.nr, wave.nt, sizeof(float), survey.tr, 0, wave.reczfile);

            //! Adjoint image
            sjawrtmfd2d(&source, &survey, &geo, &wave);

            //! Laplace filter
            sjlaplcefilter2d(geo.ipp2d, survey.nx, survey.nz);

            //! Stacking
            sjprojaddeq2d(geo.gipp2d, geo.ipp2d, survey.x0, survey.z0, survey.nx, survey.nz);
            sjprojaddeq2d(nmig, geo.nipp2d, survey.x0, survey.z0, survey.nx, survey.nz);

            //! Free
            sjmfree2d(geo.vp2d);
            sjmfree2d(geo.ipp2d);
            sjmfree2d(geo.nipp2d);
            sjmfree2d(wave.recz);
            sjmfree2d(wave.snapz2d);

            //! Time
            tend = (double) clock();
            printf("Single shot RTM complete - %d/%d - time=%fs.\n", is + 1, survey.ns,
                   (tend - tstart) / CLOCKS_PER_SEC);
            printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n", rankid,
                   survey.sx + survey.x0, survey.sz + survey.z0,
                   survey.rx[0] + survey.x0, survey.rx[survey.nr - 1] + survey.x0,
                   survey.rz[0] + survey.z0, survey.rz[survey.nr - 1] + survey.z0);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (rankid == 0) {
            //! Reduce
            MPI_Reduce(MPI_IN_PLACE, geo.gipp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, nmig[0], survey.gnx * survey.gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

            //! Source
            sjprojdiveq2d(geo.gipp2d, nmig, 0, 0, survey.gnx, survey.gnz);

            //! Cut surface
            for (ix = 0; ix < survey.gnx; ++ix)
                for (iz = 0; iz < 50; ++iz)
                    geo.gipp2d[ix][iz] = 0.0f;

            //! Output
            sjwritesuall(geo.gipp2d[0], survey.gnx, survey.gnz, geo.ds, geo.ippfile);

            //! Time
            Tend = (double) clock();
            printf("Acoustic RTM complete - time=%fs.\n\n\n", (Tend - Tstart) / CLOCKS_PER_SEC);
        } else {
            //! Reduce
            MPI_Reduce(geo.gipp2d[0], geo.gipp2d[0], survey.gnx * survey.gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(nmig[0], nmig[0], survey.gnx * survey.gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        //! Free
        sjmfree1d(source.wavelet);
        sjmfree2d(geo.gipp2d);
        sjmfree2d(geo.gvp2d);
        sjmfree2d(nmig);

    } else {
        printf("\nExamples:   sjmpiartm2d survey=survey.su vp=vp.su recz=recz.su ipp=mig.su\n");
        sjbasicinformation();
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}
