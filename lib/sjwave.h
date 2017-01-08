//
// Created by hsa on 07/01/17.
//

#ifndef SJI_SJWAVE_H
#define SJI_SJWAVE_H

#include "sjinc.h"

int sjricker1d(float *ricker, int nt, int t0, float dt, float fp, float amp);

void sjextend2d(float **input, int nx, int nz,
                int ex0, int ex1, int ez0, int ez1, float **output);

void sjextract2d(float **input, int x0, int z0, int nx, int nz, float **output);

void sjfilter2d(float **a, int n2, int n1, char *mode);

void sjsetsurface(float **a, int n2, int n1, float val);

/**********************************************************************************************/
/* ! Finite Difference                                                                        */
/**********************************************************************************************/

//! Two dimension constant density acoustic forward simulation
void sjawfd2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic scatter forward simulation
void sjaswfd2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension acoustic reverse time simulation based on constant density equation
void sjawrtmfd2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

#endif //SJI_SJWAVE_H
