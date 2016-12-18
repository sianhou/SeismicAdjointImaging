//
// Created by hsa on 12/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    if (argc == 1) {
        printf("\n2D acoustic RTM using constant density velocity-stress equation.\n\n");
        printf("Parameters:\n");
        printf("svy:          Input filename of seismic survey.\n");
        printf("vp:           Input filename of seismic P-model.\n");
        printf("rec:          Input filename of seismic profile.\n");
        printf("mig:          Output filename of RTM image (SU).\n");
        printf("Other parameters:\n");
        printf("dt:           Interval of time(s), default = 0.001.\n");
        printf("ds:           Interval of space(m), default = 10.0.\n");
        printf("k1:           Peak position of wavelet, default = 30.\n");
        printf("srcrange:     Total range of source, default = 10.\n");
        printf("srctrunc:     Total time of source, default = 301.\n");
        printf("fp:           Peak frequency of wavelet, default = 20.\n");
        printf("amp:          Peak amplitude of wavelet, default = 1.0.\n");
        printf("srcdecay:     Decay of source, default = 0.4.\n");
        printf("nb:           Range of ABC, default = 15.\n");
        printf("jsnap:        Interval of snap, default = 2.\n");
        printf("ycutdirect:   Cut direct wave , 0: didn't cut ;\n");
        printf("                                1: cut direct (default).\n");
        printf("ompnum:       Number of OpenMP threads, default = 4.\n");
        printf("\nExamples:   sjartm2d svy=survey.su vp=vp.su rec=profile.su mig=image.su\n");
        sjbasicinformation();
    } else {
        //------------------------ Initlization ------------------------//
        //! Define parameters
        int is = 0, ir = 0, ompnum;
        double tstart, tend, Tstart, Tend, runtime, Runtime;
        //! Read parameters
        if (!sjmgeti("ompnum", ompnum)) ompnum = 4;
        //! Set parameters
#ifdef GFDOPENMP_
        omp_set_num_threads(ompnum);
#endif

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
            printf("ERROR: Should input vp in program sgartm2d!\n");
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
        int ycutdirect;
        char *recfile;
        //! Read parameters
        if (!sjmgeti("ycutdirect", ycutdirect)) ycutdirect = 1;
        if (!sjmgets("rec", recfile)) {
            printf("ERROR: Should input rec in program sgartm2d!\n");
            exit(0);
        }

        //------------------------ Source ------------------------//
        //! Define parameters
        int nt, k1, srcrange, srctrunc, jsnap, nsnap;
        float dt, fp, amp, srcdecay;
        float *wavelet = NULL;
        //! Read parameters
        nt = sjgetsun1(sizeof(float), recfile);
        if (!sjmgeti("k1", k1)) k1 = 30;
        if (!sjmgeti("srcrange", srcrange)) srcrange = 10;
        if (!sjmgeti("srctrunc", srctrunc)) srctrunc = 301;
        if (!sjmgeti("jsnap", jsnap)) jsnap = 2;
        if (!sjmgetf("dt", dt)) dt = 0.001;
        if (!sjmgetf("fp", fp)) fp = 20.0;
        if (!sjmgetf("amp", amp)) amp = 10.0;
        if (!sjmgetf("srcdecay", srcdecay)) srcdecay = 0.4;
        //! Allocate memory
        wavelet = (float *) sjalloc1d(nt, sizeof(float));
        //! Calculate parameters
        nsnap = (nt - 1) / jsnap + 1;
        sjricker1d(wavelet, nt, k1, dt, fp, amp);

        //------------------------ RTM ------------------------//
        //! Define parameters
        char *rtmfile;
        float **rtm = NULL;
        rtm = (float **) sjalloc2d(svy.gxl, svy.gzl, sizeof(float));
        //! Calculate parameters
        memset(rtm[0], 0, svy.gxl * svy.gzl * sizeof(float));
        //! Initialize parameters
        if (!sjmgets("mig", rtmfile)) {
            printf("ERROR: Should output mig in program sgartm2d!\n");
            exit(0);
        }

        //------------------------ Start ------------------------//
        printf("\nAcoustic RTM start.\n\n");
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
            float **profilef = NULL, **profileb = NULL, ***snapf = NULL, ***snapb = NULL;
            //! Allocate memory
            profilef = (float **) sjalloc2d(svy.nr, nt, sizeof(float));
            profileb = (float **) sjalloc2d(svy.nr, nt, sizeof(float));
            snapf = (float ***) sjalloc3d(nsnap, svy.lxl, svy.lzl, sizeof(float));
            snapb = (float ***) sjalloc3d(nsnap, svy.lxl, svy.lzl, sizeof(float));
            //! Read parameters
            sjreadsu(profileb[0], svy.nr, nt, sizeof(float), svy.tr, 0, recfile);

            //------------------------ RTM ------------------------//
            //! Define parameters
            float **image = NULL;
            image = (float **) sjalloc2d(svy.lxl, svy.lzl, sizeof(float));
            //! Calculate parameters
            memset(image[0], 0, svy.lxl * svy.lzl * sizeof(float));

            //------------------------ Forward ------------------------//
            sjawsgfd2d(nt, svy.sx, svy.sz, srcrange, srctrunc, //! Source
                       dt, srcdecay, wavelet,
                       svy.lxl, svy.lzl, ds, lvp, nb, //! Model
                       svy.nr, rx, rz, //! BC & survey
                       1, jsnap, profilef, snapf, //! Wavefield
                       ycutdirect);

            //------------------------ Backward ------------------------//
            sjawrtsgfd2d(nt, dt, //! Source
                         svy.lxl, svy.lzl, ds, lvp, //! Model
                         nb, //! Boundary condition
                         svy.nr, rx, rz, //! Survey
                         1, jsnap, profileb, snapb); //! Wavefield

            //------------------------ Image ------------------------//
            sjimage2d(snapf, snapb, nsnap, svy.lxl, svy.lzl, 1, image);

            sjimagefilter2d(image, svy.lxl, svy.lzl, 1);

            sjprojaddeq2d(rtm, image, svy.lx0, svy.lz0, svy.lxl, svy.lzl);

            //------------------------ Information ------------------------//
#ifdef GFDOPENMP_
            tend = omp_get_wtime();
            runtime = tend - tstart;
#else
            tend = (double)clock();
            runtime = (tend - tstart) / CLOCKS_PER_SEC;
#endif

            printf("Single shot RTM complete - %d/%d - time=%fs.\n", is + 1, ns, runtime);
            printf("sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n\n",
                   svy.sx + svy.lx0, svy.sz + svy.lz0, rx[0] + svy.lx0, rx[svy.nr - 1] + svy.lx0, rz[0] + svy.lz0,
                   rz[svy.nr - 1] + svy.lz0);

            //------------------------ Free memory ------------------------//
            sjmcheckfree2d(lvp);
            sjmcheckfree2d(profilef);
            sjmcheckfree2d(profileb);
            sjmcheckfree3d(snapf);
            sjmcheckfree3d(snapb);
            sjmcheckfree2d(image);
        }

        //------------------------ Output ------------------------//
        sjwritesuall(rtm[0], svy.gxl, svy.gzl, ds, rtmfile);

        //------------------------ Information ------------------------//
#ifdef GFDOPENMP_
        Tend = omp_get_wtime();
        Runtime = Tend - Tstart;
#else
        Tend = (double)clock();
        Runtime = (Tend - Tstart) / CLOCKS_PER_SEC;
#endif
        printf("Acoustic RTM complete - time=%fs.\n", Runtime);

        //------------------------ Free memory ------------------------//
        sjmcheckfree1d(wavelet);
        sjmcheckfree1d(rx);
        sjmcheckfree1d(rz);
        sjmcheckfree2d(gvp);
        sjmcheckfree2d(rtm);
    }
}