//
// Created by 侯思安 on 2016/12/28.
//

#include "../lib/sjinc.h"
#include <mpi.h>

//! Two dimension constant density acoustic LSRTM gradient
void sjlsartmgrad2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

int main(int argc, char *argv[]) {

    //! Runtime
    int is = 0, iter = 0, flag = 1, ix, iz;
    double tstart, tend, Tstart, Tend;

    //! MPI
    int rankid, nrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //! Survey
    sjssurvey sur;
    flag &= sjssurvey_init(&sur);
    flag &= sjssurvey_getparas(&sur, argc, argv);

    //! Model
    sjsgeology geo;
    flag &= sjsgeo_init(&geo);
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "vp");
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "ipp");
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "lsipp");

    //! Wave
    sjswave wav;
    flag &= sjswave_init(&wav);
    flag &= sjswave_getparas(&wav, argc, argv, "recz");

    //! Option
    sjsoption opt;
    flag &= sjsoption_init(&opt);
    flag &= sjsoption_getparas(&opt, argc, argv);

    //------------------------ LSARTM2D ------------------------//
    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic LSRTM start ------------------------\n");
        }

        //! Allocate memory
        geo.gvp2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gipp2d = sjmflloc2d(sur.gnx, sur.gnz);

        //! Set model
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);
        sjreadsuall(geo.gipp2d[0], sur.gnx, sur.gnz, geo.ippfile);

        //! Process migration
        for (iter = 0; iter < 2; ++iter) {

            //! Calculate gradient
            sjlsartmgrad2d(&sur, &geo, &wav, &opt);
        }

        if (rankid == 0) {
            //! Output
            sjwritesuall(geo.gipp2d[0], sur.gnx, sur.gnz, opt.ds, geo.lsippfile);

            //! Time
            Tend = (double) clock();
            printf("Acoustic LSRTM complete - time=%fs.\n", (Tend - Tstart) / CLOCKS_PER_SEC);
        }

    } else {
        if (rankid == 0) {
            printf("\nExamples:   sjmpilsartm2d sur=sur.su vp=vp.su recz=recz.su ipp=mig.su lsipp=lsmig.su\n");
            sjbasicinformation();
        }
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}

//! Two dimension constant density acoustic LSRTM gradient
void sjlsartmgrad2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int is = 0;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Bcast model ------------------------//
    float **nmig = sjmflloc2d(sur->gnx, sur->gnz);
    MPI_Bcast(geo->gipp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(geo->gvp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //------------------------ Calculating gradient ------------------------//
    //! Informaiton
    if (rankid == 0) {
        tstart = (double) clock();
        printf("Calculating gradient -                     ");
    }

    //! Optimization
    for (is = rankid; is < sur->ns; is += nrank) {
        //! Memory
        geo->vp2d = sjmflloc2d(sur->nx, sur->nz);
        geo->ipp2d = sjmflloc2d(sur->nx, sur->nz);
        geo->nipp2d = sjmflloc2d(sur->nx, sur->nz);
        float **recz = sjmflloc2d(sur->nr, opt->nt);
        wav->recz = sjmflloc2d(sur->nr, opt->nt);
        wav->snapz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Set model
        sjextract2d(geo->gvp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->vp2d);
        sjextract2d(geo->gipp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->ipp2d);

        //! Simulation
        sjreadsu(recz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->reczfile);
        sjawfd2d(sur, geo, wav, opt);

        //! Difference wavefield
        sjvecsubf(wav->recz[0], sur->nr * opt->nt, 1.0f, wav->recz[0], 1.0, recz[0]);

        //! Adjoint image
        opt->ystacksrc = 1;
        memset(geo->ipp2d[0], 0, sur->nx * sur->nz * sizeof(float));
        memset(geo->nipp2d[0], 0, sur->nx * sur->nz * sizeof(float));
        sjawrtmfd2d(sur, geo, wav, opt);

        //! Stacking
        sjvecaddf(geo->gipp2d[sur->x0], sur->nx * sur->nz, 1.0f, geo->gipp2d[sur->x0], 1.0f, geo->ipp2d[0]);
        sjvecaddf(nmig[sur->x0], sur->nx * sur->nz, 1.0f, nmig[sur->x0], 1.0f, geo->nipp2d[0]);

        //! Free
        sjmcheckfree2d(geo->vp2d);
        sjmcheckfree2d(geo->ipp2d);
        sjmcheckfree2d(geo->nipp2d);
        sjmcheckfree2d(wav->recz);
        sjmcheckfree2d(recz);
        sjmcheckfree3d(wav->snapz2d);

        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //! Communication
    if (rankid == 0) {

        MPI_Reduce(MPI_IN_PLACE, geo->gipp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, nmig[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        //! Source
        sjvecdivf(geo->gipp2d[0], sur->gnx * sur->gnz, 1.0, geo->gipp2d[0], nmig[0], 0.00001f);

        //! Cut surface
        sjsetsurface(geo->gipp2d, sur->gnx, 30, 0.0f);

        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        MPI_Reduce(geo->gipp2d[0], geo->gipp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(nmig[0], nmig[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    sjmcheckfree2d(nmig);
}