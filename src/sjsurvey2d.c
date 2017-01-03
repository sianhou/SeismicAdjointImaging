//
// Created by hsa on 09/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("\nCreative a sample 2d survey for adjoint imaging.\n\n");
        printf("vel:     Input a velocity file.\n");
        printf("sx0:     Begin of shot in inline(x or n2) direction.\n");
        printf("sz0:     Begin of shot in vertical(z or n1) direction.\n");
        printf("ns:      Number of shots in survey.\n");
        printf("dsx:     Interval of each shot in inline(x or n2) direction.\n");
        printf("rx0:     Begin of receiver in inline(x or n2) direction.\n");
        printf("rz0:     Begin of receiver in vertical(z or n1) direction.\n");
        printf("nr:      Number of receivers in one shot in inline(x or n2) direction.\n");
        printf("drx:     Interval of receiver in inline(x or n2) direction in each shot.\n");
        printf("drz:     Interval of receiver in vertical(z or n1) direction in each shot.\n");
        printf("x0:     Begin of local model in inline(x or n2) direction.\n");
        printf("nx:     Local model length in inline(x or n2) direction.\n");
        printf("survey:  Output name of survey file.\n\n");
        printf("Explain: for (is = 0; is < ns; ++is)\n");
        printf("             sx = sx0 + is*dsx\n");
        printf("             sz = sz0\n");
        printf("             x0 = x0 + is*dsx\n");
        printf("             z0 = z0\n");
        printf("             for (ir = 0; ir < nr; ++ir)\n");
        printf("                 rx = rx0 + ir*drx + is*dsx\n");
        printf("                 rz = rz0 + ir*drz\n\n");
        printf("Example: sjsurvey2d vel=vel.su sx0=11 sz0=5 ns=2 dsx=10 dsz=0\n");
        printf("                    rx0=51 rz0=5 nr=5 drx=2 drz=0 x0=0 lxl=201\n");
        printf("                    survey=svy.su\n");
        sjbasicinformation();
    }
    //! Read parameters
    char *inputname;
    if (!sjmgets("vel", inputname)) {
        printf("ERROR: Should input model in program sjsurvey2d!\n");
        exit(0);
    };
    char *outputname;
    if (!sjmgets("survey", outputname)) {
        printf("ERROR: Should output survey in program sjsurvey2d!\n");
        exit(0);
    };

    int sx0, sz0, ns, dsx;
    if (!sjmgeti("sx0", sx0)) {
        printf("ERROR: Should intput sx0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("sz0", sz0)) {
        printf("ERROR: Should intput sxz in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("ns", ns)) {
        printf("ERROR: Should intput ns in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("dsx", dsx)) {
        printf("ERROR: Should intput dsx in program sjsurvey2d!\n");
        exit(0);
    };
    int rx0, rz0, nr, drx, drz;
    if (!sjmgeti("rx0", rx0)) {
        printf("ERROR: Should intput rx0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("rz0", rz0)) {
        printf("ERROR: Should intput rz0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("nr", nr)) {
        printf("ERROR: Should intput nr in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("drx", drx)) {
        printf("ERROR: Should intput drx in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("drz", drz)) {
        printf("ERROR: Should intput drz in program sjsurvey2d!\n");
        exit(0);
    };
    int x0, nx;
    if (!sjmgeti("x0", x0)) {
        printf("ERROR: Should intput x0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("nx", nx)) {
        printf("ERROR: Should intput nx in program sjsurvey2d!\n");
        exit(0);
    };

    int is, ir;

    //! Get Size of model
    int vpnx, vpnz;
    vpnx = sjgetsun2(sizeof(float), inputname);
    vpnz = sjgetsun1(sizeof(float), inputname);

    //! Check parameters
    for (is = 0; is < ns; ++is) {
        if ((x0 + is * dsx) < 0 || (x0 + is * dsx) > (vpnx - 1)) {
            printf("ERROR: The begin of local model exceed the velocity model!");
            exit(0);
        }
        if ((x0 + is * dsx + nx - 1) < 0 || (x0 + is * dsx + nx - 1) > (vpnx - 1)) {
            printf("ERROR: The end of local model exceed the velocity model!");
            exit(0);
        }
        if (sx0 < x0 || sx0 > (x0 + nx - 1)) {
            printf("ERROR: The sx exceed the local model!");
            exit(0);
        }
        if (sz0 < 0 || sz0 > (vpnz - 1)) {
            printf("ERROR: The sz exceed the local model!");
            exit(0);
        }
        if (rx0 < x0 || rx0 > (x0 + nx - 1)) {
            printf("ERROR: The rx exceed the local model!");
            exit(0);
        }
        if ((rx0 + (nr - 1) * drx) < x0 || (rx0 + (nr - 1) * drx) > (x0 + nx - 1)) {
            printf("ERROR: The rx exceed the local model!");
            exit(0);
        }
        if (rz0 < 0 || rz0 > (vpnz - 1)) {
            printf("ERROR: The rz exceed the local model!");
            exit(0);
        }
        if ((rz0 + (nr - 1) * drz) < 0 || (rz0 + (nr - 1) * drz) > (vpnz - 1)) {
            printf("ERROR: The rz exceed the local model!");
            exit(0);
        }
    }

    //! Output survey
    sjssurvey survey;
    sjssurvey_init(&survey);
    survey.surveyfile = outputname;
    survey.maxnr = nr;
    survey.ry = (int *) sjalloc1d(nr, sizeof(int));
    survey.rx = (int *) sjalloc1d(nr, sizeof(int));
    survey.rz = (int *) sjalloc1d(nr, sizeof(int));
    for (is = 0; is < ns; ++is) {
        //! sy,sx,sz
        survey.sy = 0;
        survey.sx = sx0 - x0;
        survey.sz = sz0;
        //! y0,x0,z0
        survey.y0 = 0;
        survey.x0 = x0 + is * dsx;
        survey.z0 = 0;
        //! yl,xl,zl
        survey.ny = 1;
        survey.nx = nx;
        survey.nz = vpnz;
        //! gyl,gxl,gzl
        survey.gny = 1;
        survey.gnx = vpnx;
        survey.gnz = vpnz;
        //! nr
        survey.nr = nr;
        //! tr
        if (is != 0) {
            survey.tr += nr;
        } else {
            survey.tr = 0;
        }
        //! receiver
        for (ir = 0; ir < nr; ++ir) {
            survey.ry[ir] = 0;
            survey.rx[ir] = rx0 - x0 + ir * drx;
            survey.rz[ir] = rz0 + ir * drz;
        }
        //! Output
        sjssurvey_write((void *) &survey, is);
    }

    return 0;
}