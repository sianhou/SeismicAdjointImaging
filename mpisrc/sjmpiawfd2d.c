//
// Created by hsa on 16/12/16.
//

#include <mpi.h>
#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    //! Runtime
    int flag = 1, is = 0;
    double tstart, tend, Tstart, Tend;

    //! MPI
    int mpiid, rankid, nrank;
    MPI_Status stauts;
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

    //! Wave
    sjswave wave2d;
    flag &= sjswave_init(&wave2d);
    flag &= sjswave_getparas(&wave2d, argc, argv, "recz");

    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf(" ---------------- 2D Acoustic simulation start  ---------------- \n");
        }

        //! Source
        source.wavelet = (float *) sjalloc1d(source.srctrunc, sizeof(float));
        sjricker1d(source.wavelet, source.srctrunc, source.k1, source.dt, source.fp, source.amp);

        //! Model
        geo2d.gvp2d = (float **) sjalloc2d(survey.gnx, survey.gnz, sizeof(float));

        sjreadsuall(geo2d.gvp2d[0], survey.gnx, survey.gnz, geo2d.vpfile);

        //! Simulation
        for (is = rankid; is < survey.ns; is += nrank) {
            //! Time
            tstart = (double) clock();

            //! Survey
            sjssurvey_readis(&survey, is);

            //! Model
            geo2d.vp2d = (float **) sjalloc2d(survey.nx, survey.nz, sizeof(float));
            sjextract2d(geo2d.gvp2d, survey.x0, survey.z0, survey.nx, survey.nz, geo2d.vp2d);

            //! Wavefield
            wave2d.recz = (float **) sjalloc2d(survey.nr, wave2d.nt, sizeof(float));
            wave2d.snapz2d = (float ***) sjalloc3d(wave2d.nsnap, survey.nx, survey.nz, sizeof(float));

            //! Simulation
            sjawfd2d(&source, &survey, &geo2d, &wave2d);

            //! Output
            if (rankid == 0) {
                sjwritesu(wave2d.recz[0], survey.nr, wave2d.nt, sizeof(float), source.dt, is, wave2d.reczfile);
                tend = (double) clock();
                printf("Single shot simulation complete - %d/%d - time=%fs.\n", is + 1, survey.ns,
                       (tend - tstart) / CLOCKS_PER_SEC);
                printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n", rankid,
                       survey.sx + survey.x0, survey.sz + survey.z0,
                       survey.rx[0] + survey.x0, survey.rx[survey.nr - 1] + survey.x0,
                       survey.rz[0] + survey.z0, survey.rz[survey.nr - 1] + survey.z0);

                for (mpiid = 1; mpiid < nrank; mpiid++) {
                    if ((is + mpiid) < survey.ns) {
                        MPI_Recv(wave2d.recz[0], survey.nr * wave2d.nt, MPI_FLOAT, mpiid, 99, MPI_COMM_WORLD, &stauts);
                        //! Output in rank != 0
                        sjwritesu(wave2d.recz[0], survey.nr, wave2d.nt, sizeof(float), source.dt, is + mpiid,
                                  wave2d.reczfile);
                    }
                }
            } else {
                MPI_Send(wave2d.recz[0], survey.nr * wave2d.nt, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
                tend = (double) clock();
                printf("Single shot simulation complete - %d/%d - time=%fs.\n", is + 1, survey.ns,
                       (tend - tstart) / CLOCKS_PER_SEC);
                printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n", rankid,
                       survey.sx + survey.x0, survey.sz + survey.z0,
                       survey.rx[0] + survey.x0, survey.rx[survey.nr - 1] + survey.x0,
                       survey.rz[0] + survey.z0, survey.rz[survey.nr - 1] + survey.z0);
            }

            //! Free
            sjcheckfree2d((void **) geo2d.vp2d);
            sjcheckfree2d((void **) wave2d.recz);
            sjcheckfree3d((void ***) wave2d.snapz2d);
        }

        //------------------------ Information ------------------------//
        if (rankid == 0) {
            Tend = (double) clock();
            printf("2D acoustic simulation complete - time=%fs.\n", (Tend - Tstart) / CLOCKS_PER_SEC);
        }

    } else {
        printf("\nExamples:   sjmpiawfd2d survey=survey.su vp=vp.su recz=recz.su nt=3001\n");
        sjbasicinformation();
    }


    MPI_Finalize();

    return 0;
}
