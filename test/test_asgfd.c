//
// Created by hsa on 07/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    //!  data

    int **data1 = (int **) sjalloc2d(301, 201, sizeof(int));
    int **data2 = (int **) sjalloc2d(301, 201, sizeof(int));

    int ii, jj;

    for (ii = 0; ii < 301; ++ii) {
        for (jj = 0; jj < 201; ++jj) {
            data1[ii][jj] = ii + jj;
            data2[ii][jj] = ii - jj;
        }
    }

    sjwritesu(data1[0], 301, 201, sizeof(int), 0.001, 0, "data1");
    sjwritesu(data1[0], 301, 201, sizeof(int), 0.001, 0, "data12");
    sjwritesu(data2[0], 301, 201, sizeof(int), 0.001, 1, "data12");

    //sjreadsu(data1[0], 301, 201, sizeof(int), 301 * (201 * sizeof(int) + 240L), "data12");

    for (ii = 0; ii < 301; ++ii) {
        for (jj = 0; jj < 201; ++jj) {
            if (data1[ii][jj] != data2[ii][jj]) printf("error: %d, %d!\n",ii,jj);
        }
    }
    printf("ok!\n");

    printf("%d,%d!\n",data1[0][0],data2[0][0]);

    printf("n2=%d\n", sjgetsun2(sizeof(float), "data12"));
    printf("n1=%d\n", sjgetsun1(sizeof(float), "data12"));

    return 0;
}
