#! /bin/bash

#-----------------------------------------------------------------------------
# Creative a survey data
#-----------------------------------------------------------------------------

../../bin/sgsurvey2d sgin=marvel.su sx0=361 sz0=5 ns=2 dsx=10 dsz=0 rx0=211 rz0=5 nr=301 drx=1 drz=0 lx0=201 lxl=321 sgot=marsvy.su

#-----------------------------------------------------------------------------
# Simulation with GFDXYZ
#-----------------------------------------------------------------------------

# Openmp
#../../bin/sjawsgfd2d svy=marsvy.su vp=marvel.su rec=marrec.su nt=3001 dt=0.001 ompnum=8

# MPI
mpirun -np 2 ../../bin/sgmpiawsgfd2d svy=marsvy.su vp=marvel.su rec=marrec.su nt=3001 dt=0.001

#-----------------------------------------------------------------------------
# RTM with GFDXYZ
#-----------------------------------------------------------------------------

# Openmp
#../../bin/sgartm2d svy=marsvy.su vp=marvel.su rec=marrec.su mig=marmig.su ompnum=8

# MPI
mpirun -np 2 ../../bin/sgmpiartm svy=marsvy.su vp=marvel.su rec=marrec.su mig=marmig_mpi.su
