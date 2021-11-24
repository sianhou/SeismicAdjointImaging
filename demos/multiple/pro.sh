#! /bin/bash

binpath=../../bin
outpath=.
nthd=2

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=model.bin n2=301 n1=301 su=$outpath/vp.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=4 nr=291 vel=$outpath/vp.su x0=0 nx=301 dx0=0 sx0=151 sz0=20 dsx=10 rx0=5 rz0=20 drx=1 drz=0 survey=$outpath/survey.d


#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.d vp=$outpath/vp.su profz=$outpath/recz.su0 nt=2001 dt=0.002 fp=20 ds=10 ycutdirect=0 yfreebc=0 ysnap=1
mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.d vp=$outpath/vp.su profz=$outpath/recz.su1 nt=2001 dt=0.002 fp=20 ds=10 ycutdirect=0 yfreebc=1 ysnap=1

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

# mpirun -np $nthd $binpath/sjmpiartm2d survey=survey.d vp=evps.su profz=recz.su nt=3001 dt=0.002 izz=mig.su fp=20 ds=20.0

# sfsuread < vp.su endian=0 | sfput d1=10 d2=10 label1="Depth" label2="Distance" unit1="m" unit2="m" title="model" | sfgrey color=j > model.vpl
# vpconvert bgcolor=l format=png *.png

# sfsuread < recz.su0 endian=0 | sfput d1=0.002 d2=10 label1="Time" label2="Distance" unit1="ms" unit2="m" title="Absorbing BC"| sfwindow n2=291 | sfgrey > recz0.vpl
# sfsuread < recz.su1 endian=0 | sfput d1=0.002 d2=10 label1="Time" label2="Distance" unit1="ms" unit2="m" title="free surface" | sfwindow n2=291 | sfgrey > recz1.vpl
# vpconvert bgcolor=white format=jpg *.vpl

# sfsuread < recz.su0-snap-0 endian=0 | sfput n1=301 n2=301 n3=2001 d1=10 d2=10 label1="Depth" label2="Distance" unit1="m" unit2="m" title="Absorbing" | sfwindow j3=40 | sfgrey gainpanel=a > recz0-snap.vpl
# sfsuread < recz.su1-snap-0 endian=0 | sfput n1=301 n2=301 n3=2001 d1=10 d2=10 label1="Depth" label2="Distance" unit1="m" unit2="m" title="Free" | sfwindow j3=40 | sfgrey gainpanel=a > recz1-snap.vpl
# vpconvert bgcolor=white format=gif *.vpl