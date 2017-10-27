//
// Authors: Hou, Sian - sianhou1987@outlook.com
//          Wang, Guangchao - wgcupc@163.com

#include "../lib/sjinc.h"
#include <mpi.h>

int main(int argc, char *argv[]) {

    //! Runtime
    int is = 0, flag = 1;
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
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "izz");

    //! Wave
    sjswave wav;
    flag &= sjswave_init(&wav);
    flag &= sjswave_getparas(&wav, argc, argv, "profz");

    //! Option
    sjsoption opt;
    flag &= sjsoption_init(&opt);
    flag &= sjsoption_getparas(&opt, argc, argv);

    //! ------------------------ RTM2D ------------------------
    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic Illumation start ------------------------\n");
        }

        //! Model
        geo.gizz2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gnzz2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gvp2d = sjmflloc2d(sur.gnx, sur.gnz);
        
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);
        MPI_Bcast(geo.gvp2d[0], sur.gnx * sur.gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

        //! Rtm
        for (is = rankid; is < sur.ns; is += nrank) {
            //! Time
            tstart = (double) clock();

            //! Survey
            sjssurvey_readis(&sur, is);

            //! Model
            geo.vp2d = sjmflloc2d(sur.nx, sur.nz);
            geo.izz2d = sjmflloc2d(sur.nx, sur.nz);
            geo.nzz2d = sjmflloc2d(sur.nx, sur.nz);
            sjextract2d(geo.vp2d, sur.x0, sur.z0, sur.nx, sur.nz, geo.gvp2d);

            //! Forward simulaion
            wav.profz = sjmflloc2d(sur.nr, opt.nt);
            wav.fwz2d = sjmflloc3d(opt.nsnap, sur.nx, sur.nz);
            sjafor2d(&sur, &geo, &wav, &opt);

            //! Read record
            sjreadsu(wav.profz[0], sur.nr, opt.nt, sizeof(float), sur.tr, 0, wav.profzfile);

            //! Adjoint image
            opt.ystacksrc = 0;
            sjaillumation2d(&sur, &geo, &wav, &opt);
            
            //! Stacking
	    memset(geo.gizz2d[0], 0, sur.gnx * sur.gnz * sizeof(float));
            memset(geo.gnzz2d[0], 0, sur.gnx * sur.gnz * sizeof(float));
            sjvecaddf(&geo.gizz2d[sur.x0][sur.z0], sur.nx * sur.nz, 1.0f, &geo.gizz2d[sur.x0][sur.z0], 1.0f, &geo.izz2d[0][0]);
            sjvecaddf(&geo.gnzz2d[sur.x0][sur.z0], sur.nx * sur.nz, 1.0f, &geo.gnzz2d[sur.x0][sur.z0], 1.0f, &geo.nzz2d[0][0]);

	    char *file = (char *) malloc(1024 * sizeof(char));
            sprintf(file, "%s-source-%d", geo.izzfile, is);
	    sjwritesu(geo.gizz2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, 0, file);
            sprintf(file, "%s-receiver-%d", geo.izzfile, is);
            sjwritesu(geo.gnzz2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, 0, file);

            //! Free
            sjmfree2d(geo.vp2d);
            sjmfree2d(geo.izz2d);
            sjmfree2d(geo.nzz2d);
            sjmfree2d(wav.profz);
            sjmfree3d(wav.fwz2d);

            //! Time
            tend = (double) clock();
            printf("Single shot Illumation complete - %d/%d - time=%fs.\n", is + 1, sur.ns,
                   (tend - tstart) / CLOCKS_PER_SEC);
            printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n", rankid,
                   sur.sx + sur.x0, sur.sz + sur.z0,
                   sur.rx[0] + sur.x0, sur.rx[sur.nr - 1] + sur.x0,
                   sur.rz[0] + sur.z0, sur.rz[sur.nr - 1] + sur.z0);
        }

        MPI_Barrier(MPI_COMM_WORLD);
	
	if (rankid == 0) {
            printf("Acoustic Illumation complete - time=%fs.\n\n", (Tend - Tstart) / CLOCKS_PER_SEC);
        }

        //! Free
        sjmfree2d(geo.gizz2d);
        sjmfree2d(geo.gnzz2d);
        sjmfree2d(geo.gvp2d);
    } else {
        printf("\nExamples:   sjmpiartm2d survey=sur.su vp=vp.su profz=profz.su ipp=mig.su\n");
        sjbasicinformation();
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}
