#! /bin/bash

binpath=../../bin
outpath=.
nthd=1

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=evp.bin n2=2721 n1=236 su=$outpath/evp.su
$binpath/sjbin2su binary=evps.bin n2=2721 n1=236 su=$outpath/evps.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=1 nr=351 vel=$outpath/evps.su x0=310 nx=741 dx0=10 sx0=371 sz0=5 dsx=0 rx0=20 rz0=5 drx=2 drz=0  survey=$outpath/survey.d

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.d vp=$outpath/evps.su profz=$outpath/recz.su nt=3001 dt=0.002 fp=20 ds=20

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiartm2d survey=survey.d vp=evps.su profz=recz.su nt=3001 dt=0.002 izz=mig.su fp=20 ds=20.0
