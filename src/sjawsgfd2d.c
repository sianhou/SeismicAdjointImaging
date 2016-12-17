//
// Created by hsa on 12/12/16.
//

#include "../lib/sjinc.h"
#include "../lib/sjfile.h"

int main(int argc, char *argv[]) {

    if (argc == 1) {
        printf("\nSimulate 2D acoustic wavefield using constant density velocity-stress equation.\n\n");
        printf("Parameters:\n");
        printf("svy:          Input filename of seismic survey.\n");
        printf("vp:           Input filename of seismic P-model.\n");
        printf("rec:          Output filename of seismic profile.\n");
        printf("nt:           Total time of acoustic simulation.\n");
        printf("dt:           Interval of time(s), default = 0.001.\n");
        printf("ds:           Interval of space(m), default = 10.0.\n");
        printf("Other parameters:\n");
        printf("k1:           Peak position of wavelet, default = 30.\n");
        printf("srcrange:     Total range of source, default = 10.\n");
        printf("srctrunc:     Total time of source, default = 301.\n");
        printf("fp:           Peak frequency of wavelet, default = 20.\n");
        printf("amp:          Peak amplitude of wavelet, default = 1.0.\n");
        printf("srcdecay:     Decay of source, default = 0.4.\n");
        printf("nb:           Range of ABC, default = 15.\n");
        printf("ompnum:       Number of OpenMP threads, default = 4.\n");
        printf("\nExamples:   sjawsgfd2d svy=survey.su vp=vp.su rec=profile.su nt=3001\n");
        sjbasicinformation();
    } else {
        //------------------------ Initlization ------------------------//
        //! Define parameters
        int is = 0, ompnum;
        double tstart, tend, Tstart, Tend, runtime, Runtime;
        //! Read parameters
        if (!sjmgeti("ompnum", ompnum)) ompnum = 4;
        //! Set parameters
#ifdef GFDOPENMP_
        omp_set_num_threads(ompnum);
#endif

        //------------------------ Source ------------------------//
        //! Define parameters
        int nt, k1, srcrange, srctrunc;
        float dt, fp, amp, srcdecay;
        float *wavelet = NULL;
        //! Read parameters
        if (!sjmgeti("nt", nt)) {
            printf("ERROR: Should input nt in program sjawsgfd2d!\n");
            exit(0);
        }
        if (!sjmgeti("k1", k1)) k1 = 30;
        if (!sjmgeti("srcrange", srcrange)) srcrange = 10;
        if (!sjmgeti("srctrunc", srctrunc)) srctrunc = 301;
        if (!sjmgetf("dt", dt)) dt = 0.001;
        if (!sjmgetf("fp", fp)) fp = 20.0;
        if (!sjmgetf("amp", amp)) amp = 10.0;
        if (!sjmgetf("srcdecay", srcdecay)) srcdecay = 0.4;
        //! Allocate memory
        wavelet = (float *) sjalloc1d(nt, sizeof(float));
        //! Calculate parameters
        sjricker1d(wavelet, nt, k1, dt, fp, amp);

        //------------------------ Boundary condition ------------------------//
        //! Define parameters
        int nb;
        if (!sjmgeti("nb", nb)) nb = 15;

        //------------------------ Survey ------------------------//
        //! Define parameters
        sury svy;
        char *svyfile;
        int ns, nr;
        int *ry = NULL, *rx = NULL, *rz = NULL;
        //! Read parameters
        if (!sjmgets("svy", svyfile)) {
            printf("ERROR: Should input survey in program sjawsgfd2d!\n");
            exit(0);
        }
        ns = sjgetsvyns(svyfile);
        nr = sjgetsvynr(svyfile);
        //! Allocate memory
        ry = (int *) sjalloc1d(nr, sizeof(int));
        rx = (int *) sjalloc1d(nr, sizeof(int));
        rz = (int *) sjalloc1d(nr, sizeof(int));

        //------------------------ Model ------------------------//
        //! Define parameters
        char *vpfile;
        float ds;
        float **gvp = NULL;
        //! Initialize parameters
        if (!sjmgetf("ds", ds)) ds = 10.0;
        if (!sjmgets("vp", vpfile)) {
            printf("ERROR: Should input vp in program sjawsgfd2d!\n");
            exit(0);
        }
        //! Read parameters
        sjreadsurvey(0, &svy, ry, rx, rz, svyfile); //! Read the first shot's survey to get global model parameters
        //! Allocate memory
        gvp = (float **) sjalloc2d(svy.gxl, svy.gzl, sizeof(float));
        //! Read model
        sjreadsuall(gvp[0], svy.gxl, svy.gzl, vpfile);

        //------------------------ Wavefield ------------------------//
        //! Define parameters
        char *recfile;
        //! Initialize parameters
        if (!sjmgets("rec", recfile)) {
            printf("ERROR: Should output rec in program sjawsgfd2d!\n");
            exit(0);
        }
        //------------------------ Start ------------------------//
        printf("\nConstant density acoustic simulation start.\n\n");
#ifdef GFDOPENMP_
        Tstart = omp_get_wtime();
#else
        Tstart = (double)clock();
#endif

        for (is = 0; is < ns; ++is) {
            //------------------------ Information ------------------------//
            //! Time
#ifdef GFDOPENMP_
            tstart = omp_get_wtime();
#else
            tstart = (double)clock();
#endif

            //------------------------ Survey ------------------------//
            //! Read parameters
            sjreadsurvey(is, &svy, ry, rx, rz, svyfile);

            //------------------------ Model ------------------------//
            //! Define parameters
            float **lvp = NULL;
            //! Allocate memory
            lvp = (float **) sjalloc2d(svy.lxl, svy.lzl, sizeof(float));
            //! Extract model
            sjextract2d(gvp, svy.lx0, svy.lz0, svy.lxl, svy.lzl, lvp);

            //------------------------ Wavefield ------------------------//
            //! Define parameters
            float **profile = NULL, ***snap = NULL;
            //! Allocate memory
            profile = (float **) sjalloc2d(svy.nr, nt, sizeof(float));

            //------------------------ Simulation ------------------------//
            sjawsgfd2d(nt, svy.sx, svy.sz, srcrange, srctrunc, //! Source
                       dt, srcdecay, wavelet,
                       svy.lxl, svy.lzl, ds, lvp, //! Model
                       nb, //! Boundary condition
                       svy.nr, rx, rz, //! Survey
                       0, 1, profile, snap); //! Wavefield

            //------------------------ Output ------------------------//
            sjwritesu(profile[0], svy.nr, nt, sizeof(float), dt, is, recfile);

            //------------------------ Information ------------------------//
#ifdef GFDOPENMP_
            tend = omp_get_wtime();
            runtime = tend - tstart;
#else
            tend = (double)clock();
            runtime = (tend - tstart) / CLOCKS_PER_SEC;
#endif
            printf("Single shot simulation complete - %d/%d - time=%fs.\n", is + 1, ns, runtime);
            printf("sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n",
                   svy.sx + svy.lx0, svy.sz + svy.lz0, rx[0] + svy.lx0, rx[svy.nr - 1] + svy.lx0, rz[0] + svy.lz0,
                   rz[svy.nr - 1] + svy.lz0);

            //------------------------ Free memory ------------------------//
            sjmcheckfree2d(lvp);
            sjmcheckfree2d(profile);
            sjmcheckfree3d(snap);
        }

        //------------------------ Information ------------------------//
#ifdef GFDOPENMP_
        Tend = omp_get_wtime();
        Runtime = Tend - Tstart;
#else
        Tend = (double)clock();
        Runtime = (Tend - Tstart) / CLOCKS_PER_SEC;
#endif
        printf("Constant density acoustic simulation complete - time=%fs.\n", Runtime);

        //------------------------ Free memory ------------------------//
        sjmcheckfree1d(wavelet);
        sjmcheckfree1d(ry);
        sjmcheckfree1d(rx);
        sjmcheckfree1d(rz);
        sjmcheckfree2d(gvp);
    }
    return 0;
}
