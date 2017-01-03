
//
// Created by hsa on 16/12/16.
//

#include "../lib/sjinc.h"
#include "../lib/sjfile.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("\nDisplay 2D survey.\n\n");
        printf("sjin:   Input a survey file.\n");
        printf("is:      The sequence number of shot to display, default = 1.\n");
        printf("         1<= is <= ns.\n\n");
        printf("flag:    1 - Display local information, default = 1.\n");
        printf("Example: sjdisplaysurvey2d sgin=svy.su\n");
        sjbasicinformation();
    } else {
        char *inputname;
        if (!sjmgets("sjin", inputname)) {
            printf("ERROR: Should input survey in program sgdisplaysurvey2d!\n");
            exit(0);
        };

        int is, ir, ns, nr, flag;
        sury svy;
        int *ry = NULL, *rx = NULL, *rz = NULL;
        //! Get parameters
        if (!sjmgeti("is", is)) is = 1;
        if (!sjmgeti("flag", flag)) flag = 1;
        ns = sjgetsvyns(inputname);
        nr = sjgetsvynr(inputname);
        //! Allocate parameters
        ry = (int *) sjalloc1d(nr, sizeof(int));
        rx = (int *) sjalloc1d(nr, sizeof(int));
        rz = (int *) sjalloc1d(nr, sizeof(int));

        if (is > ns) {
            printf("ERROR: Shot sequence number exceed ns!\n");
            exit(0);
        } else {

            if (flag) {
                sjreadsurvey(is - 1, &svy, ry, rx, rz, inputname);
                printf("\n");
                printf("Display local survey information(in local coordinate): %d\n", ns);
                printf("\n");
                printf("Total shot number: %d\n", ns);
                printf("Single shot sequence: %d\n", is);
                printf("Single shot position: sx=%d sz=%d\n", svy.sx, svy.sz);
                printf("\n");
                printf("Total model size: gxl=%d gzl=%d\n", svy.gxl, svy.gzl);
                printf("Local model size: lxl=%d lzl=%d\n", svy.lxl, svy.lzl);
                printf("Local model begin: lx0=%d lz0=%d\n", svy.lx0, svy.lz0);
                printf("\n");
                printf("Local receiver number: nr=%d\n", svy.nr);
                printf("Local receiver position: rx=%-4d ", rx[0]);

                for (ir = 1; ir < 3; ++ir) {
                    if (ir < nr)
                        if (ir != nr - 1)
                            printf("%-4d ", rx[ir]);
                }
                printf("... %-4d\n", rx[nr - 1]);
                printf("Local receiver position: rz=%-4d ", rz[0]);
                for (ir = 1; ir < 3; ++ir) {
                    if (ir < nr)
                        if (ir != nr - 1)
                            printf("%-4d ", rz[ir]);
                }
                printf("... %-4d\n", rz[nr - 1]);
            } else {
                sjreadsurvey(is - 1, &svy, ry, rx, rz, inputname);
                printf("\n");
                printf("Display global survey information(in global coordinate): %d\n", ns);
                printf("\n");
                printf("Total shot number: %d\n", ns);
                printf("Single shot sequence: %d\n", is);
                printf("Single shot position: sx=%d sz=%d\n", svy.sx + svy.lx0, svy.sz + svy.lz0);
                printf("\n");
                printf("Total model size: gxl=%d gzl=%d\n", svy.gxl, svy.gzl);
                printf("Local model size: lxl=%d lzl=%d\n", svy.lxl, svy.lzl);
                printf("Local model begin: lx0=%d lz0=%d\n", svy.lx0, svy.lz0);
                printf("\n");
                printf("Local receiver number: nr=%d\n", svy.nr);
                printf("Local receiver position: rx=%-4d ", rx[0] + svy.lx0);

                for (ir = 1; ir < 3; ++ir) {
                    if (ir < nr)
                        if (ir != nr - 1)
                            printf("%-4d ", rx[ir] + svy.lx0);
                }
                printf("... %-4d\n", rx[nr - 1] + svy.lx0);
                printf("Local receiver position: rz=%-4d ", rz[0] + svy.lz0);
                for (ir = 1; ir < 3; ++ir) {
                    if (ir < nr)
                        if (ir != nr - 1)
                            printf("%-4d ", rz[ir] + svy.lz0);
                }
                printf("... %-4d\n", rz[nr - 1] + svy.lz0);
            }
        }
    }
    return 0;
}