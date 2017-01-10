//
// Created by hsa on 07/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    //!  data

    float **data1 = (float **) sjalloc2d(100, 150, sizeof(float));
    float **data2 = (float **) sjalloc2d(100, 150, sizeof(float));

    int ix,iz;
    for(ix=0; ix<50;ix++)
        for(iz=0;iz<75;iz++)
            data1[ix][iz] = 1000.0f;
    for(ix=50; ix<100;ix++)
        for(iz=0;iz<75;iz++)
            data1[ix][iz] = 1500.0f;
    for(ix=0; ix<50;ix++)
        for(iz=75;iz<150;iz++)
            data1[ix][iz] = 2000.0f;
    for(ix=50; ix<100;ix++)
        for(iz=75;iz<150;iz++)
            data1[ix][iz] = 2500.0f;
    sjguasssmoothf2d(data1, 100, 150, 0.4, 10, data1);
    sjwritesuall(data1[0],100,150,10.0,"data1.su");
    sjwritesuall(data2[0],100,150,10.0,"data2.su");


    return 0;
}
