//
// Created by hsa on 12/12/16.
//
#include <mpi.h>
#include "../lib/sjinc.h"

int sjartm2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

    //! Runtime
    int is = 0;
    double tstart, tend, Tstart, Tend;
    //! MPI
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //! Time
    if (rankid == 0) {
        Tstart = (double) clock();
        printf("---------------- 2D Acoustic RTM start  ----------------\n");
    }

    float **mig = (float **) sjalloc2d(survey->gnx, survey->gnz, sizeof(float));
    float **nmig = (float **) sjalloc2d(survey->gnx, survey->gnz, sizeof(float));
    float **gnmig = (float **) sjalloc2d(survey->gnx, survey->gnz, sizeof(float));

    MPI_Bcast(geo->gvp2d[0], survey->gnx * survey->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //! Rtm
    for (is = rankid; is < survey->ns; is += nrank) {
        //! Time
        tstart = (double) clock();

        //! Survey
        sjssurvey_readis(survey, is);

        //! Model
        geo->vp2d = (float **) sjalloc2d(survey->nx, survey->nz, sizeof(float));
        geo->ipp2d = (float **) sjalloc2d(survey->nx, survey->nz, sizeof(float));
        geo->nipp2d = (float **) sjalloc2d(survey->nx, survey->nz, sizeof(float));
        sjextract2d(geo->gvp2d, survey->x0, survey->z0, survey->nx, survey->nz, geo->vp2d);

        //! Wavefield
        wave->recz = (float **) sjalloc2d(survey->nr, wave->nt, sizeof(float));
        wave->snapz2d = (float ***) sjalloc3d(wave->nsnap, survey->nx, survey->nz, sizeof(float));

        //! Forward simulaion
        sjawsgfd2d(source, survey, geo, wave);

        //! Read record
        sjreadsu(wave->recz[0], survey->nr, wave->nt, sizeof(float), survey->tr, 0, wave->reczfile);

        //! Adjoint image
        wave->yadjointbc = 1;
        wave->ycutdirect = 0;
        sjawrtfd2d(source, survey, geo, wave);

        //! Laplace filter
        sjlaplcefilter2d(geo->ipp2d, survey->nx, survey->nz);

        //! Stacking
        sjprojaddeq2d(mig, geo->ipp2d, survey->x0, survey->z0, survey->nx, survey->nz);
        sjprojaddeq2d(nmig, geo->nipp2d, survey->x0, survey->z0, survey->nx, survey->nz);

        //! Free
        sjmcheckfree2d(geo->vp2d);
        sjmcheckfree2d(geo->ipp2d);
        sjmcheckfree2d(geo->nipp2d);
        sjmcheckfree2d(wave->recz);
        sjmcheckfree3d(wave->snapz2d);

        //! Time
        tend = (double) clock();
        printf("Single shot RTM complete - %d/%d - time=%fs.\n", is + 1, survey->ns,
               (tend - tstart) / CLOCKS_PER_SEC);
        printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n", rankid,
               survey->sx + survey->x0, survey->sz + survey->z0,
               survey->rx[0] + survey->x0, survey->rx[survey->nr - 1] + survey->x0,
               survey->rz[0] + survey->z0, survey->rz[survey->nr - 1] + survey->z0);
    }

    //! Communication
    MPI_Reduce(mig[0], geo->gipp2d[0], survey->gnx * survey->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(nmig[0], gnmig[0], survey->gnx * survey->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //! Source compenate
    sjprojdiveq2d(geo->gipp2d, gnmig, 0, 0, survey->gnx, survey->gnz);

    //------------------------ Information ------------------------//
    if (rankid == 0) {
        Tend = (double) clock();
        printf("Acoustic RTM complete - time=%fs.\n", (Tend - Tstart) / CLOCKS_PER_SEC);
    }

    //! Free
    sjmcheckfree2d(mig);
    sjmcheckfree2d(nmig);
    sjmcheckfree2d(gnmig);

    return 1;
}

int main(int argc, char *argv[]) {

    //! Runtime
    int flag = 1;

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
    sjsgeo geo2d;
    flag &= sjsgeo_init(&geo2d);
    flag &= sjsgeo_getparas2d(&geo2d, argc, argv, "vp");
    flag &= sjsgeo_getparas2d(&geo2d, argc, argv, "ipp");

    //! Wave
    sjswave wave2d;
    flag &= sjswave_init(&wave2d);
    flag &= sjswave_getparas(&wave2d, argc, argv, "recz");

    if (flag) {

        //! Source
        source.wavelet = (float *) sjalloc1d(source.srctrunc, sizeof(float));
        sjricker1d(source.wavelet, source.srctrunc, source.k1, source.dt, source.fp, source.amp);

        //! Model
        geo2d.gvp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        sjreadsuall(geo2d.gvp2d[0], survey.gnx, survey.gnz, geo2d.vpfile);

        //! Migration
        geo2d.gipp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));
        sjartm2d(&source, &survey, &geo2d, &wave2d);

        //! Output
        sjwritesuall(geo2d.gipp2d[0], survey.gnx, survey.gnz, geo2d.ds, geo2d.ippfile);

    } else {
        printf("\nExamples:   sjmpiartm2d survey=survey.su vp=vp.su recz=recz.su ipp=mig.su\n");
        sjbasicinformation();
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}
