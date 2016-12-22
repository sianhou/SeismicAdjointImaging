#! /bin/bash

path=../../

#-----------------------------------------------------------------------------
# Creative a survey data
#-----------------------------------------------------------------------------

$path/sjsurvey2d vel=marvel.su sx0=361 sz0=5 ns=2 dsx=10 dsz=0 rx0=211 rz0=5 nr=301 drx=1 drz=0 x0=201 nx=321 survey=marsvy.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np 2 $path/sjmpiawsgfd2d survey=marsvy.su vp=marvel.su recz=marrec.su nt=3001 dt=0.001

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

# Openmp
#../../bin/sjartm2d svy=marsvy.su vp=marvel.su rec=marrec.su mig=marmig.su ompnum=8

# MPI
mpirun -np 2 $path/sjmpiartm2d survey=marsvy.su vp=marvel.su recz=marrec.su mig=marmig.su nt=3001 dt=0.001
