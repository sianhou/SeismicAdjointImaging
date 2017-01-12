//
// Created by hsa on 12/01/17.
//

#include "../lib/sjinc.h"
#include <mpi.h>

//! Two dimension constant density acoustic WTI gradient
void sjawtigrad2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt, float **grad, float *err);

int main(int argc, char *argv[]) {

    //! Runtime
    int iter = 0, flag = 1;
    double Tstart, Tend;

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

    //! Wave
    sjswave wav;
    flag &= sjswave_init(&wav);
    flag &= sjswave_getparas(&wav, argc, argv, "profz");

    //! Option
    sjsoption opt;
    flag &= sjsoption_init(&opt);
    flag &= sjsoption_getparas(&opt, argc, argv);

    //------------------------ LSAWTI2D ------------------------//
    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic WTI start ------------------------\n");
        }
        //! Set model
        geo.gvp2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gipp2d = sjmflloc2d(sur.gnx, sur.gnz);
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);
        memset(geo.gipp2d[0], 0, sur.gnx * sur.gnz * sizeof(float));

        opt.niter = 0;

        //! Process migration
        float **g0 = sjmflloc2d(sur.gnx, sur.gnz);
        float **g1 = sjmflloc2d(sur.gnx, sur.gnz);
        float **cg = sjmflloc2d(sur.gnx, sur.gnz);
        float *err = sjmflloc1d(opt.niter);

        //! Inversion
        do {
            //! Calculate gradient
            sjawtigrad2d(&sur, &geo, &wav, &opt, g1, err + iter);

            //! Optimization
            if (rankid == 0) {
                //! CG
                sjcgsolver(geo.gipp2d[0], sur.gnx * sur.gnz, cg[0], g1[0], g0[0], iter);

                /*

                char *file = (char *) malloc(1024 * sizeof(char));
                sprintf(file, "g1-%d.su", iter);
                sjwritesuall(g1[0], sur.gnx, sur.gnz, opt.ds, file);
                sprintf(file, "cg-%d.su", iter);
                sjwritesuall(cg[0], sur.gnx, sur.gnz, opt.ds, file);
                sprintf(file, "li-%d.su", iter);
                sjwritesuall(geo.gipp2d[0], sur.gnx, sur.gnz, opt.ds, file);

                //! Information
                Tend = (double) clock();
                printf("Acoustic WTI complete - %2d/%2d - residual=%e - time=%6.2fs.\n",
                       iter + 1, opt.niter, err[iter], (Tend - Tstart) / CLOCKS_PER_SEC);
                       */
            }

            iter += 1;

        } while (iter < opt.niter);

        if (rankid == 0) {
            printf("Acoustic WTI completed.\n\n");
        }

        sjmfree2d(geo.gvp2d);
        sjmfree2d(geo.gipp2d);
        sjmfree2d(cg);
        sjmfree2d(g1);
        sjmfree2d(g0);
    } else {
        if (rankid == 0) {
            printf("\nExamples:   sjmpilsartm2d sur=sur.su vp=vp.su profz=profz.su ipp=lsipp.su\n");
            sjbasicinformation();
        }
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}

//! Two dimension constant density acoustic LSRTM gradient
void sjawtigrad2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt, float **grad, float *err) {

    //------------------------ Runtime ------------------------//
    int is = 0;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Bcast model ------------------------//
    *err = 0.0f;
    memset(grad[0], 0, sur->gnx * sur->gnz * sizeof(float));
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
        geo->spp2d = sjmflloc2d(sur->nx, sur->nz);
        wav->profz = sjmflloc2d(sur->nr, opt->nt);
        wav->fwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);
        wav->bwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);

        //! Forward exploration
        sjaswfd2d(sur, geo, wav, opt);

        //! Backward exploration
        opt->rtmopt = 2;
        sjawrtmfd2d(sur, geo, wav, opt);

        //! Time shift imaging
        int tshift = 20;
        float ***tsimage = sjmflloc3d(2 * tshift + 1, sur->nx, sur->nz);
        float ***tsillum = sjmflloc3d(2 * tshift + 1, sur->nx, sur->nz);
        int it, ishift, ix, iz;
        for (ishift = -tshift; ishift <= tshift; ++ishift) {
            for (it = tshift; it < opt->nt - tshift; ++it) {
                for (ix = 0; ix < sur->nx; ++ix) {
                    for (iz = 0; iz < sur->nz; ++iz) {
                        tsimage[ishift + tshift][ix][iz] +=
                                wav->fwz2d[it - ishift][ix][iz] * wav->bwz2d[it + ishift][ix][iz];
                        tsillum[ishift + tshift][ix][iz] +=
                                wav->fwz2d[it - ishift][ix][iz] * wav->fwz2d[it - ishift][ix][iz];
                    }
                }
            }
        }

        sjwritesuall(tsimage[0][0], (2 * tshift + 1) * sur->nx, sur->nz, 10.0f, "tsimage.su");
        sjwritesuall(tsimage[0][0], (2 * tshift + 1) * sur->nx, sur->nz, 10.0f, "tsillum.su");
    }
}