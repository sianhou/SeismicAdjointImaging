#! /bin/bash

binpath=../../bin
outpath=.
nthd=1

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=model.bin n2=301 n1=301 su=$outpath/vp.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=137 nr=291 vel=$outpath/vp.su x0=0 nx=301 dx0=0 sx0=151 sz0=20 dsx=10 rx0=5 rz0=20 drx=1 drz=0 survey=$outpath/survey.d


#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.d vp=$outpath/vp.su profz=$outpath/recz.su0 nt=2001 dt=0.002 fp=20 ds=10 ycutdirect=0 yfreebc=0
mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.d vp=$outpath/vp.su profz=$outpath/recz.su1 nt=2001 dt=0.002 fp=20 ds=10 ycutdirect=0 yfreebc=1

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

# mpirun -np $nthd $binpath/sjmpiartm2d survey=survey.d vp=evps.su profz=recz.su nt=3001 dt=0.002 izz=mig.su fp=20 ds=20.0
