//
// Created by hsa on 20/12/16.
//

#include "sjsimulation.h"

//! Survey
int sjssurvey_init(sjssurvey *ptr) {
    ptr->sy = -1;
    ptr->sx = -1;
    ptr->sz = -1;

    ptr->y0 = -1;
    ptr->x0 = -1;
    ptr->z0 = -1;

    ptr->ny = -1;
    ptr->nx = -1;
    ptr->nz = -1;

    ptr->gny = -1;
    ptr->gnx = -1;
    ptr->gnz = -1;

    ptr->nr = -1;
    ptr->tr = -1;

    ptr->ns = -1;
    ptr->maxnr = -1;

    ptr->nparas = 14;

    ptr->ry = NULL;
    ptr->rx = NULL;
    ptr->rz = NULL;
    ptr->surveyfile = NULL;

    return 1;
}

int sjssurvey_display(sjssurvey *ptr) {
    printf("Display survey information:\n");
    printf("sy=%d\n", ptr->sy);
    printf("sx=%d\n", ptr->sx);
    printf("sz=%d\n", ptr->sz);
    printf("y0=%d\n", ptr->y0);
    printf("x0=%d\n", ptr->x0);
    printf("z0=%d\n", ptr->z0);
    printf("ny=%d\n", ptr->ny);
    printf("nx=%d\n", ptr->nx);
    printf("nz=%d\n", ptr->nz);
    printf("gny=%d\n", ptr->gny);
    printf("gnx=%d\n", ptr->gnx);
    printf("gnz=%d\n", ptr->gnz);
    printf("nr=%d\n", ptr->nr);
    printf("tr=%d\n", ptr->tr);
    printf("ns=%d\n", ptr->ns);
    printf("maxnr=%d\n", ptr->maxnr);
    printf("nparas=%d\n", ptr->nparas);
    printf("ry=%p\n", ptr->ry);
    printf("rx=%p\n", ptr->rx);
    printf("rz=%p\n", ptr->rz);
    printf("surveyfile=%s\n", ptr->surveyfile);
    return 1;
}

int sjssurvey_readis(sjssurvey *ptr, int is) {

    int nr = ptr->maxnr;
    int np = ptr->nparas;

    if (is < ptr->ns) {
        //! Allocate memory
        sjcheckfree1d((void *) ptr->ry);
        sjcheckfree1d((void *) ptr->rx);
        sjcheckfree1d((void *) ptr->rz);
        ptr->ry = (int *) sjalloc1d(nr, sizeof(int));
        ptr->rx = (int *) sjalloc1d(nr, sizeof(int));
        ptr->rz = (int *) sjalloc1d(nr, sizeof(int));

        //! Get parameters
        sjreadsu(ptr, 1, np, sizeof(int), is, 0, ptr->surveyfile);

        //! get ry rx rz
        sjreadsu(ptr->ry, 1, nr, sizeof(int), is, np + 0 * nr, ptr->surveyfile);
        sjreadsu(ptr->rx, 1, nr, sizeof(int), is, np + 1 * nr, ptr->surveyfile);
        sjreadsu(ptr->rz, 1, nr, sizeof(int), is, np + 2 * nr, ptr->surveyfile);
        return 1;
    } else {
        printf("ERROR: Is exceed ns in function sjssurvey_readis()!");
        return 0;
    }
}

int sjssurvey_write(sjssurvey *ptr, int ifappend) {
    int nr = ptr->maxnr;
    int np = ptr->nparas;
    int *tmp = (int *) sjalloc1d(np + 3 * nr, sizeof(int));

    memcpy(tmp + 0 * np + 0 * nr, ptr, np * sizeof(int));
    memcpy(tmp + 1 * np + 0 * nr, ptr->ry, nr * sizeof(int));
    memcpy(tmp + 1 * np + 1 * nr, ptr->rx, nr * sizeof(int));
    memcpy(tmp + 1 * np + 2 * nr, ptr->rz, nr * sizeof(int));

    sjwritesu((void *) tmp, 1, np + 3 * nr, sizeof(int), 1.0f, ifappend, ptr->surveyfile);

    return 1;
}

int sjssurvey_getparas(sjssurvey *ptr, int argc, char **argv) {

    if (argc == 1) {
        printf("* survey:     Input filename of seismic survey.\n");
        return 0;
    } else {
        if (!sjmgets("survey", ptr->surveyfile)) {
            printf("ERROR: Should set survey file!\n");
            return 0;
        };

        //! Get ns
        ptr->ns = sjgetsun2(sizeof(int), ptr->surveyfile);

        //! Get maxnr
        ptr->maxnr = (sjgetsun1(sizeof(int), ptr->surveyfile) - ptr->nparas) / 3;

        //! Get parameters
        sjreadsu(ptr, 1, ptr->nparas, sizeof(int), 0, 0, ptr->surveyfile);

        return 1;
    }
}

//! Source
int sjssource_init(sjssource *ptr) {

    ptr->k1 = -1;
    ptr->srcrange = -1;
    ptr->srctrunc = -1;

    ptr->dt = -1.0f;
    ptr->fp = -1.0f;
    ptr->amp = -1.0f;
    ptr->srcdecay = -1.0f;
    ptr->wavelet = NULL;

    return 1;
}

int sjssource_display(sjssource *ptr) {
    printf("Display source information:\n");
    printf("k1=%d\n", ptr->k1);
    printf("srcrange=%d\n", ptr->srcrange);
    printf("srctrunc=%d\n", ptr->srctrunc);
    printf("dt=%f\n", ptr->dt);
    printf("fp=%f\n", ptr->fp);
    printf("amp=%f\n", ptr->amp);
    printf("srcdecay=%f\n", ptr->srcdecay);
    printf("wavelet=%p\n", ptr->wavelet);
    return 1;
}

int sjssource_getparas(sjssource *ptr, int argc, char **argv) {
    if (argc == 1) {
        printf("  dt:         Interval of time(s), default = 0.001.\n");
        printf("  k1:         Peak position of wavelet, default = 50.\n");
        printf("  srcrange:   Total range of source, default = 10.\n");
        printf("  srctrunc:   Total time of source, default = 301.\n");
        printf("  fp:         Peak frequency of wavelet, default = 20.\n");
        printf("  amp:        Peak amplitude of wavelet, default = 10.0.\n");
        printf("  srcdecay:   Decay of source, default = 0.4.\n");
        return 0;
    } else {
        if (!sjmgetf("dt", ptr->dt)) ptr->dt = 0.001;
        if (!sjmgeti("k1", ptr->k1)) ptr->k1 = 50;
        if (!sjmgeti("srcrange", ptr->srcrange)) ptr->srcrange = 10;
        if (!sjmgeti("srctrunc", ptr->srctrunc)) ptr->srctrunc = 301;
        if (!sjmgetf("fp", ptr->fp)) ptr->fp = 20.0;
        if (!sjmgetf("amp", ptr->amp)) ptr->amp = 10.0;
        if (!sjmgetf("srcdecay", ptr->srcdecay)) ptr->srcdecay = 0.4;
        return 1;
    }
}

//! geo-model
int sjsgeo_init(sjsgeo *ptr) {
    ptr->nb = -1;
    ptr->ds = -1.0f;

    ptr->gvp2d = NULL;
    ptr->gvs2d = NULL;

    ptr->vp2d = NULL;
    ptr->vp2d = NULL;

    ptr->vpfile = NULL;
    ptr->vsfile = NULL;

    return 1;
}

int sjsgeo_display(sjsgeo *ptr) {
    printf("Display geo information:\n");
    printf("nb=%d\n", ptr->nb);
    printf("ds=%f\n", ptr->ds);
    printf("gvp2d=%p\n", ptr->gvp2d);
    printf("gvs2d=%p\n", ptr->gvs2d);
    printf("vp2d=%p\n", ptr->vp2d);
    printf("vpfile=%s\n", ptr->vpfile);
    printf("vsfile=%s\n", ptr->vsfile);
    return 1;
}

int sjsgeo_getparas2d(sjsgeo *ptr, int argc, char **argv, char *info) {
    if (strcmp(info, "vp") == 0) {
        if (argc == 1) {
            printf("* vp:         2D seismic P-velocity.\n");
            printf("  ds:         Interval of space(m), default = 10.0.\n");
            printf("  nb:         Range of ABC, default = 15.\n");
            return 0;
        } else {
            if (!sjmgets("vp", ptr->vpfile)) {
                printf("ERROR: Should input 2D vp file!\n");
                exit(0);
            }
            if (!sjmgetf("ds", ptr->ds)) ptr->ds = 10.0;
            if (!sjmgeti("nb", ptr->nb)) ptr->nb = 15;

            return 1;
        }
    }
}

//! wave
int sjswave_init(sjswave *ptr) {
    ptr->nt = -1;
    ptr->jsnap = -1;
    ptr->nsnap = -1;
    ptr->ycutdirect = -1;
    ptr->recy = NULL;
    ptr->recx = NULL;
    ptr->recz = NULL;

    ptr->snapy2d = NULL;
    ptr->snapx2d = NULL;
    ptr->snapz2d = NULL;

    ptr->recyfile = NULL;
    ptr->recxfile = NULL;
    ptr->reczfile = NULL;

    return 1;
}

int sjswave_display(sjswave *ptr) {
    printf("Display wave information:\n");
    printf("nt=%d\n", ptr->nt);
    printf("jsnap=%d\n", ptr->jsnap);
    printf("nsnap=%d\n", ptr->nsnap);
    printf("ycutdirect=%d\n", ptr->ycutdirect);

    printf("recy=%p\n", ptr->recy);
    printf("recx=%p\n", ptr->recx);
    printf("recz=%p\n", ptr->recz);
    printf("snapy2d=%p\n", ptr->snapy2d);
    printf("snapx2d=%p\n", ptr->snapx2d);
    printf("snapz2d=%p\n", ptr->snapz2d);

    printf("recyfile=%s\n", ptr->recyfile);
    printf("recxfile=%s\n", ptr->recxfile);
    printf("reczfile=%s\n", ptr->reczfile);

    return 1;
}

int sjswave_getparas(sjswave *ptr, int argc, char **argv, char *info) {
    if (strcmp(info, "recy") == 0) {
        if (argc == 1) {
            printf("* recy:       Seicmic record in y (cxline) - direction.\n");
            printf("* nt:         Total time.\n");
            printf("  ycutdirect: Cut direct wave, 0: didn't cut;\n");
            printf("                               1: cut direct (default).\n");
            printf("  jsnap:      Interval of snap, default = 2.\n");
            return 0;
        } else {
            if (!sjmgets("recy", ptr->recyfile)) {
                printf("ERROR: Should set recy!\n");
                exit(0);
            }
            if (!sjmgeti("nt", ptr->nt)) {
                printf("ERROR: Should set nt!\n");
                exit(0);
            }
            if (!sjmgeti("ycutdirect", ptr->ycutdirect)) ptr->ycutdirect = 1;
            if (!sjmgeti("jsnap", ptr->jsnap)) ptr->jsnap = 0;
            if (!sjmgeti("jsnap", ptr->jsnap)) ptr->jsnap = 2;
            ptr->nsnap = (ptr->nt - 1) / ptr->jsnap + 1;
            return 1;
        }
    }
    if (strcmp(info, "recx") == 0) {
        if (argc == 1) {
            printf("* recx:       Seicmic record in x (inline) - direction.\n");
            printf("* nt:         Total time.\n");
            printf("  ycutdirect: Cut direct wave, 0: didn't cut;\n");
            printf("                               1: cut direct (default).\n");
            printf("  jsnap:      Interval of snap, default = 2.\n");
            return 0;
        } else {
            if (!sjmgets("recx", ptr->recxfile)) {
                printf("ERROR: Should set recx!\n");
                exit(0);
            }
            if (!sjmgeti("nt", ptr->nt)) {
                printf("ERROR: Should set nt!\n");
                exit(0);
            }
            if (!sjmgeti("ycutdirect", ptr->ycutdirect)) ptr->ycutdirect = 1;
            if (!sjmgeti("jsnap", ptr->jsnap)) ptr->jsnap = 2;
            ptr->nsnap = (ptr->nt - 1) / ptr->jsnap + 1;

            return 1;
        }
    }
    if (strcmp(info, "recz") == 0) {
        if (argc == 1) {
            printf("* recz:       Seicmic record in z (depth) - direction.\n");
            printf("* nt:         Total time.\n");
            printf("  ycutdirect: Cut direct wave, 0: didn't cut;\n");
            printf("                               1: cut direct (default).\n");
            printf("  jsnap:      Interval of snap, default = 2.\n");
            return 0;
        } else {
            if (!sjmgets("recz", ptr->reczfile)) {
                printf("ERROR: Should set recz!\n");
                exit(0);
            }
            if (!sjmgeti("nt", ptr->nt)) {
                printf("ERROR: Should set nt!\n");
                exit(0);
            }
            if (!sjmgeti("ycutdirect", ptr->ycutdirect)) ptr->ycutdirect = 1;
            if (!sjmgeti("jsnap", ptr->jsnap)) ptr->jsnap = 2;
            ptr->nsnap = (ptr->nt - 1) / ptr->jsnap + 1;
            return 1;
        }
    }

    return 0;
}

/**********************************************************************************************/
/* ! Finite Difference                                                                        */
/**********************************************************************************************/

//! Two dimension acoustic simulation based on constant velocity-stress equation
void sjawsgfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

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
    float ids = -dt / ds;
    float **cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //! Wave
    int nt = wave->nt;
    int jsnap = wave->jsnap;
    int ycutdirect = wave->ycutdirect;

    //! Boundary condition
    float **gxl = (float **) sjalloc2d(nzb, 3, sizeof(float));
    float **gxr = (float **) sjalloc2d(nzb, 3, sizeof(float));
    float **gzu = (float **) sjalloc2d(nxb, 3, sizeof(float));
    float **gzb = (float **) sjalloc2d(nxb, 3, sizeof(float));

    //! Wavefield
    float **vx0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **vx1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **vz0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **vz1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));


    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitohabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    for (ix = 0; ix < nxb; ix++)
        for (iz = 0; iz < nzb; iz++)
            cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids;
    //! Wavefield
    memset(vx0[0], 0, nxb * nzb * sizeof(float));
    memset(vx1[0], 0, nxb * nzb * sizeof(float));
    memset(vz0[0], 0, nxb * nzb * sizeof(float));
    memset(vz1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));

    /**********************************************************************************************/
    /* ! Wavefield                                                                                */
    /**********************************************************************************************/

    //! Wavefield exploration
    for (it = 0; it < nt; it++) {
        //! Source
        if (it < srctrunc)
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p0[ix + sx][iz + sz] += wav[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate veloctiy
        for (ix = marg; ix < nxb - marg; ix++) {
            for (iz = marg; iz < nzb - marg; iz++) {
                vx1[ix][iz] = sjmsgfd2dn2(p0, ix, iz) * ids + vx0[ix][iz];
                vz1[ix][iz] = sjmsgfd2dn1(p0, ix, iz) * ids + vz0[ix][iz];
            }
        }

        //! Calculate stress
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p1[ix][iz] = cp[ix][iz] * (sjmsgfd2dn2(vx1, ix - 1, iz) + sjmsgfd2dn1(vz1, ix, iz - 1)) + p0[ix][iz];

        //! Boundary condition
        sjapplyohabc2d(vx1, vx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(vz1, vz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wave->recz[ir][it] = p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(vx0[0], vx1[0], nxb * nzb * sizeof(float));
        memcpy(vz0[0], vz1[0], nxb * nzb * sizeof(float));
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
    }

    /**********************************************************************************************/
    /* ! Cut direct wave                                                                          */
    /**********************************************************************************************/
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];

        //------------------------ Initialization ------------------------//
        //! Boundary condition
        sjinitohabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
        //! Model
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids;
        //! Wavefield
        memset(vx0[0], 0, nxb * nzb * sizeof(float));
        memset(vx1[0], 0, nxb * nzb * sizeof(float));
        memset(vz0[0], 0, nxb * nzb * sizeof(float));
        memset(vz1[0], 0, nxb * nzb * sizeof(float));
        memset(p0[0], 0, nxb * nzb * sizeof(float));
        memset(p1[0], 0, nxb * nzb * sizeof(float));

        //! Wavefield exploration
        for (it = 0; it < nt; it++) {
            //! Source
            if (it < srctrunc)
                for (ix = -srcrange; ix <= srcrange; ix++)
                    for (iz = -srcrange; iz <= srcrange; iz++)
                        p0[ix + sx][iz + sz] += wav[it] * expf(-srcdecay * (ix * ix + iz * iz));

            //! Calculate veloctiy
            for (ix = marg; ix < nxb - marg; ix++) {
                for (iz = marg; iz < nzb - marg; iz++) {
                    vx1[ix][iz] = sjmsgfd2dn2(p0, ix, iz) * ids + vx0[ix][iz];
                    vz1[ix][iz] = sjmsgfd2dn1(p0, ix, iz) * ids + vz0[ix][iz];
                }
            }

            //! Calculate stress
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++)
                    p1[ix][iz] =
                            cp[ix][iz] * (sjmsgfd2dn2(vx1, ix - 1, iz) + sjmsgfd2dn1(vz1, ix, iz - 1)) + p0[ix][iz];

            //! Boundary condition
            sjapplyohabc2d(vx1, vx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            sjapplyohabc2d(vz1, vz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            sjapplyohabc2d(p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++)
                wave->recz[ir][it] -= p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

            //! Update
            memcpy(vx0[0], vx1[0], nxb * nzb * sizeof(float));
            memcpy(vz0[0], vz1[0], nxb * nzb * sizeof(float));
            memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(vx0);
    sjmfree2d(vx1);
    sjmfree2d(vz0);
    sjmfree2d(vz1);
    sjmfree2d(p0);
    sjmfree2d(p1);
}

//! Two dimension acoustic reverse time simulation based on constant velocity-stress equation
void sjawrtsgfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

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
    float ids = -dt / ds;
    float **cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //! Wave
    int nt = wave->nt;
    int jsnap = wave->jsnap;

    //! Boundary condition
    float **gxl = (float **) sjalloc2d(nzb, 3, sizeof(float));
    float **gxr = (float **) sjalloc2d(nzb, 3, sizeof(float));
    float **gzu = (float **) sjalloc2d(nxb, 3, sizeof(float));
    float **gzb = (float **) sjalloc2d(nxb, 3, sizeof(float));

    //! Wavefield
    float **vx0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **vx1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **vz0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **vz1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p0 = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    float **p1 = (float **) sjalloc2d(nxb, nzb, sizeof(float));

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitohabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    for (ix = 0; ix < nxb; ix++)
        for (iz = 0; iz < nzb; iz++)
            cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids;
    //! Wavefield
    memset(vx0[0], 0, nxb * nzb * sizeof(float));
    memset(vx1[0], 0, nxb * nzb * sizeof(float));
    memset(vz0[0], 0, nxb * nzb * sizeof(float));
    memset(vz1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));

    /**********************************************************************************************/
    /**/ //! Wavefield                                                                          /**/
    /**********************************************************************************************/

    //! Wavefield revese time exploration
    for (it = nt - 1; it >= 0; --it) {
        //! Source
        for (ir = 0; ir < nr; ir++)
            p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wave->recz[ir][it];

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++) {
            for (iz = marg; iz < nzb - marg; iz++) {
                vx0[ix][iz] = sjmsgfd2dn2(p1, ix, iz) * ids + vx1[ix][iz];
                vz0[ix][iz] = sjmsgfd2dn1(p1, ix, iz) * ids + vz1[ix][iz];
            }
        }

        //! Calculate stress
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p0[ix][iz] = cp[ix][iz] * (sjmsgfd2dn2(vx0, ix - 1, iz) + sjmsgfd2dn1(vz0, ix, iz - 1)) + p0[ix][iz];

        //! Boundary condition
        sjapplyohabc2d(vx0, vx1, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(vz0, vz1, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplyohabc2d(p0, p1, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    geo->ipp2d[ix - nb - marg][iz - nb - marg] +=
                            wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p0[ix][iz];


        //! Update
        memcpy(&p1[0][0], &p0[0][0], nxb * nzb * sizeof(float));
        memcpy(&vx1[0][0], &vx0[0][0], nxb * nzb * sizeof(float));
        memcpy(&vz1[0][0], &vz0[0][0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(vx0);
    sjmfree2d(vx1);
    sjmfree2d(vz0);
    sjmfree2d(vz1);
    sjmfree2d(p0);
    sjmfree2d(p1);
}

//! Two dimension acoustic simulation based on constant density equation
void sjawfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

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
    float ids = dt * dt / ds / ds;
    float **cp = (float **) sjalloc2d(nxb, nzb, sizeof(float));
    sjextend2d(geo->vp2d, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, cp);

    //! Wave
    int nt = wave->nt;
    int jsnap = wave->jsnap;
    int ycutdirect = wave->ycutdirect;

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

    //! Wavefield exploration
    for (it = 0; it < nt; it++) {

        //! Source
        if (it < srctrunc)
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p1[ix + sx][iz + sz] += wav[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate veloctiy
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] =
                        cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) + 2.0f * p1[ix][iz] - p0[ix][iz];


        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wave->recz[ir][it] = p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
    }

    /**********************************************************************************************/
    /* ! Cut direct wave                                                                          */
    /**********************************************************************************************/
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];

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

        //! Wavefield exploration
        for (it = 0; it < nt; it++) {

            //! Source
            if (it < srctrunc)
                for (ix = -srcrange; ix <= srcrange; ix++)
                    for (iz = -srcrange; iz <= srcrange; iz++)
                        p1[ix + sx][iz + sz] += wav[it] * expf(-srcdecay * (ix * ix + iz * iz));

            //! Calculate veloctiy
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++)
                    p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) + 2.0f * p1[ix][iz] -
                                 p0[ix][iz];

            //! Boundary condition
            sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++)
                wave->recz[ir][it] -= p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

            //! Update
            memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
            memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
}

//! Two dimension acoustic reverse time simulation based on constant density equation
void sjawrtmfd2d(sjssource *source, sjssurvey *survey, sjsgeo *geo, sjswave *wave) {

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

        //! Source
        for (ir = 0; ir < nr; ir++)
            p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wave->recz[ir][it];

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p0[ix][iz] =
                        cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) + 2.0f * p1[ix][iz] - p2[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p0, p1, p2, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0) {
            for (ix = nb + marg; ix < nxb - nb - marg; ix++) {
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    geo->ipp2d[ix - nb - marg][iz - nb - marg] +=
                            wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p1[ix][iz];
                    geo->nipp2d[ix - nb - marg][iz - nb - marg] +=
                            wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg] *
                            wave->snapz2d[it / jsnap][ix - nb - marg][iz - nb - marg];
                    }
                }
        }

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