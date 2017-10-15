//
// Creat by wgc
// 
#ifndef SJI_SJENVELOPEMISFIT_H
#define SJI_SJENVELOPEMISFIT_H

#include "sjinc.h" 

//�� Convolution For one Trace
void convolve_cwp (int lx, int ifx, float *x,int ly, int ify, float *y,int lz, int ifz, float *z);

//! Hilbert Transform 
void hilbert (int n, float x[], float y[]);

// Get Envelope For one Trace
void envelope(int n,float *x,float *ex,float *eh);

//! Misfit for envelope misfit function
//! N=number of time interval P is envelope power always set 2 
void sjmisfit_evelope(float **syn, float **obs, float **misfit, int n, int trace_num, int p);

#endif
