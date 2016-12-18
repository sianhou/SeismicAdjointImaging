//
// Created by hsa on 09/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("\nCreative a sample 2d survey for RTM (writed in int).\n\n");
        printf("vel:     Input a velocity file.\n");
        printf("sx0:     Begin of shot in inline(x or n2) direction.\n");
        printf("sz0:     Begin of shot in vertical(z or n1) direction.\n");
        printf("ns:      Number of shots in survey.\n");
        printf("dsx:     Interval of each shot in inline(x or n2) direction.\n");
        printf("rx0:     Begin of receiver in inline(x or n2) direction.\n");
        printf("rz0:     Begin of receiver in vertical(z or n1) direction.\n");
        printf("nrx:     Number of receivers in one shot in inline(x or n2) direction.\n");
        printf("drx:     Interval of receiver in inline(x or n2) direction in each shot.\n");
        printf("drz:     Interval of receiver in vertical(z or n1) direction in each shot.\n");
        printf("lx0:     Begin of local model in inline(x or n2) direction.\n");
        printf("lxl:     Local model length in inline(x or n2) direction.\n");
        printf("sjot:    Output name of survey file.\n\n");
        printf("Explain: for (is = 0; is < ns; ++is)\n");
        printf("             sx = sx0 + is*dsx\n");
        printf("             sz = sz0\n");
        printf("             lx0 = lx0 + is*dsx\n");
        printf("             lz0 = lz0\n");
        printf("             for (ir = 0; ir < nr; ++ir)\n");
        printf("                 rx = rx0 + ir*drx + is*dsx\n");
        printf("                 rz = rz0 + ir*drz\n\n");
        printf("Example: sjsurvey2d vel=vel.su sx0=11 sz0=5 ns=2 dsx=10 dsz=0\n");
        printf("                    rx0=51 rz0=5 nr=5 drx=2 drz=0 lx0=0 lxl=201\n");
        printf("                    sjot=svy.su\n");
        sjbasicinformation();
    } else {
        char *inputname;
        if (!sjmgets("vel", inputname)) {
            printf("ERROR: Should input model in program sjsurvey2d!\n");
            exit(0);
        };
        char *outputname;
        if (!sjmgets("sjot", outputname)) {
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
        int lx0, lxl;
        if (!sjmgeti("lx0", lx0)) {
            printf("ERROR: Should intput lx0 in program sjsurvey2d!\n");
            exit(0);
        };
        if (!sjmgeti("lxl", lxl)) {
            printf("ERROR: Should intput mxl in program sjsurvey2d!\n");
            exit(0);
        };

        int is, ir;

        //! Get Size of model
        int vpnx, vpnz;
        vpnx = sjgetsun2(sizeof(float), inputname);
        vpnz = sjgetsun1(sizeof(float), inputname);

        //! Calculate survey
        int **survey = (int **) sjalloc2d(ns, sjgetsvynum() + 3 * nr, sizeof(int));

        for (is = 0; is < ns; ++is) {
            //! sy sx sz
            survey[is][0] = 0;
            survey[is][1] = sx0 - lx0;
            survey[is][2] = sz0;
            //! vyl vxl vzl
            survey[is][3] = 0;
            survey[is][4] = vpnx;
            survey[is][5] = vpnz;
            //! ly0 lx0 lz0
            survey[is][6] = 0;
            survey[is][7] = lx0 + is * dsx;
            survey[is][8] = 0;
            //! lyl lxl lzl
            survey[is][9] = 1;
            survey[is][10] = lxl;
            survey[is][11] = vpnz;
            //! nr
            survey[is][12] = nr;
            //! tr
            if (is != 0) {
                survey[is][13] = survey[is - 1][13] + nr;
            } else {
                survey[is][13] = 0;
            }
            //! rec
            for (ir = 0; ir < nr; ++ir) {
                //! recy recx recz
                survey[is][14 + ir] = 0;
                survey[is][14 + 1 * nr + ir] = rx0 - lx0 + ir * drx;
                survey[is][14 + 2 * nr + ir] = rz0 + ir * drz;
            }

            //! Check data
            if ((lx0 + is * dsx) < 0 || (lx0 + is * dsx) > (vpnx - 1)) {
                printf("ERROR: The begin of local model exceed the velocity model!");
                exit(0);
            }
            if ((lx0 + is * dsx + lxl) < 0 || (lx0 + is * dsx + lxl) > (vpnx - 1)) {
                printf("ERROR: The end of local model exceed the velocity model!");
                exit(0);
            }
            if (sx0 < lx0 || sx0 > (lx0 + lxl - 1)) {
                printf("ERROR: The sx exceed the local model!");
                exit(0);
            }
            if (sz0 < 0 || sz0 > (vpnz - 1)) {
                printf("ERROR: The sz exceed the local model!");
                exit(0);
            }

            if (rx0 < lx0 || rx0 > (lx0 + lxl - 1)) {
                printf("ERROR: The rx exceed the local model!");
                exit(0);
            }
            if ((rx0 + (nr - 1) * drx) < lx0 || (rx0 + (nr - 1) * drx) > (lx0 + lxl - 1)) {
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

        sjwritesu(survey[0], ns, sjgetsvynum() + 3 * nr, sizeof(int), 1.0, 0, outputname);

        return 0;
    }
}