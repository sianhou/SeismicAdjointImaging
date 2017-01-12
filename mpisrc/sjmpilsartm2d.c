//
// Created by 侯思安 on 2016/12/28.
//

#include "../lib/sjinc.h"
#include <mpi.h>

//! Two dimension constant density acoustic LSRTM gradient
void sjlsartmgrad2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt, float **grad, float *err);

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

    //------------------------ LSARTM2D ------------------------//
    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic LSRTM start ------------------------\n");
        }
        //! Set model
        geo.gvp2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gipp2d = sjmflloc2d(sur.gnx, sur.gnz);
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);
        memset(geo.gipp2d[0], 0, sur.gnx * sur.gnz * sizeof(float));

        opt.niter = 60;

        //! Process migration
        float **g0 = sjmflloc2d(sur.gnx, sur.gnz);
        float **g1 = sjmflloc2d(sur.gnx, sur.gnz);
        float **cg = sjmflloc2d(sur.gnx, sur.gnz);
        float *err = sjmflloc1d(opt.niter);

        //! Inversion
        do {
            //! Calculate gradient
            if (iter == 0) {
                opt.ystacksrc = 0;
                sjlsartmgrad2d(&sur, &geo, &wav, &opt, g1, err + iter);
            } else {
                opt.ystacksrc = 1;
                sjlsartmgrad2d(&sur, &geo, &wav, &opt, g1, err + iter);
            }

            //! Optimization
            if (rankid == 0) {
                //! CG
                sjcgsolver(geo.gipp2d[0], sur.gnx * sur.gnz, cg[0], g1[0], g0[0], iter);

                char *file = (char *) malloc(1024 * sizeof(char));
                sprintf(file, "g1-%d.su", iter);
                sjwritesuall(g1[0], sur.gnx, sur.gnz, opt.ds, file);
                sprintf(file, "cg-%d.su", iter);
                sjwritesuall(cg[0], sur.gnx, sur.gnz, opt.ds, file);
                sprintf(file, "li-%d.su", iter);
                sjwritesuall(geo.gipp2d[0], sur.gnx, sur.gnz, opt.ds, file);

                //! Information
                Tend = (double) clock();
                printf("Acoustic LSRTM complete - %2d/%2d - residual=%e - time=%6.2fs.\n",
                       iter + 1, opt.niter, err[iter], (Tend - Tstart) / CLOCKS_PER_SEC);
            }

            iter += 1;

        } while (iter < opt.niter);

        if (rankid == 0) {
            printf("Acoustic LSRTM completed.\n\n");
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
void sjlsartmgrad2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt, float **grad, float *err) {

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
        float **recz = sjmflloc2d(sur->nr, opt->nt);

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);
        sjextract2d(geo->ipp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gipp2d);

        //! Simulation
        sjreadsu(recz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->profzfile);
        sjaswfd2d(sur, geo, wav, opt);

        //! Difference wavefield
        sjvecsubf(wav->profz[0], sur->nr * opt->nt, 1.0f, wav->profz[0], 1.0f, recz[0]);
        *err += sjvecdotf(sur->nr * opt->nt, 1.0f, wav->profz[0], wav->profz[0]);

        //! Adjoint image
        memset(geo->ipp2d[0], 0, sur->nx * sur->nz * sizeof(float));
        memset(geo->spp2d[0], 0, sur->nx * sur->nz * sizeof(float));
        sjawrtmfd2d(sur, geo, wav, opt);

        //! Stacking
        sjvecaddf(grad[sur->x0], sur->nx * sur->nz, 1.0f, grad[sur->x0], 1.0f, geo->ipp2d[0]);
        sjvecaddf(nmig[sur->x0], sur->nx * sur->nz, 1.0f, nmig[sur->x0], 1.0f, geo->spp2d[0]);

        //! Free
        sjmcheckfree2d(geo->vp2d);
        sjmcheckfree2d(geo->ipp2d);
        sjmcheckfree2d(geo->spp2d);
        sjmcheckfree2d(wav->profz);
        sjmcheckfree2d(recz);
        sjmcheckfree3d(wav->fwz2d);

        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    if (rankid == 0) {
        //! Communication
        MPI_Reduce(MPI_IN_PLACE, err, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, grad[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, nmig[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        //! Source
        sjvecdivf(grad[0], sur->gnx * sur->gnz, 1.0, grad[0], nmig[0], 0.00001f);

        //! Cut source
        sjsetsurface(grad, sur->gnx, 30, 0.0);

        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        //! Communication
        MPI_Reduce(err, err, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(grad[0], grad[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(nmig[0], nmig[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    sjmcheckfree2d(nmig);
}