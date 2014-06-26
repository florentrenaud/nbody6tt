
Based on Nbody6 downloaded on 29 January 2013.

List of changes made to Sverre's Nbody6. Discard this document if you don't want to dig in the code.
All the changes are commented in the .f and .h files. They start with the tag '*** FlorentR' and end with '*** FRenaud'.
The number(s) in curly brackets indicates the line number(s) where change(s) have(s) been made.

For future updates: make sure that KZ(14)=9 is considered in xtrnld.f (for the moment, it is with KZ(14).GE.3 on line 41)

In Ncode
========

### Makefile
* add the ttcal.f, ttinit.f, ttforce.f, ttgalaxy.f files to the source list
* add a target all to make cpu version (not related to tides)

### adjust.f	{241}
* set the tidal radius to 10x the half-mass radius

### common6.h	{7,76}
* add the TT block that contains tidal tensor data
* add the boolean TTMODE

### define.f	{127}
* (mode B) add the reading of the initial position and velocity of the cluster


### escape.f	{157}
* include the case of tidal tensor when setting RTIDE

### gcint.f		{27,55}
* (mode B) call ttforce to get the force from the galaxy on the cluster

### intgrt.f		{132}
* (mode B) add the case of the tidal tensor in mode B

### jacobi.f	{8,13,18,20,39,42,44}
* (only used if KZ(21)>1)
* use of the bound mass as the cluster mass (also for KZ(14).NE.9)
* use of potential+kinetic energy as membership criterion

### mydump.f	{14,45,92,95,98,146}
* add the TT block to the saving and restoring procedures.
* use the restart.dat file. (not related to tides)
* add the TTMODE boolean

### modify.f	{75}
* call ttinit if the tensor file needs to be re-read

### nbody6.f	{31}
* print a little header to the output, in order to identify the version used to produce the simulation.

### output.f	{205}
* (mode B) output the orbit

### param.h {11}
* add the nbttmax parameter (this implies to modify the previous line too)

### regint.f	{136,150,350}
* include the case of the tidal tensor

### start.f {42}
* call ttinit if KZ(14) = 9

### ttcal.f (new file)
* Computes the tidal tensor at current time by means of interpolation.
* (mode B) update the tensor if needed

### ttforce.f	(new file)
* (mode B) computes the force and the tidal tensor at the position of the cluster, using the user-defintion of the galactic potential

### ttgalaxy.f	(new file)
* (mode B) user-definition of the galactic potenial

### ttinit.f	(new file)
* determine the TT mode (A or B)
* force the center of mass correction (KZ(31)=1)
* (mode A) read the tensor data from file 'tt.dat', scale and initialize
* (mode B) read the initial coordinates of the cluster and initialize

### xtrnlf.f	{28}
* include the computation of the external force from the tidal tensor

### xtrnlv.f	{22,42,84}
* include the potential energy from the tidal tensor (only at setup)


In GPU2
=====

### Makefile_gpu
* add the ttcal.f, ttinit.f, ttforce.f, ttgalaxy.f files to the source list

### adjust.f	{262}
* set the tidal radius to 10x the half-mass radius as in Ncode/adjust.f

### gpucor.f	{118,130,346}
* include the case of the tidal tensor, as in Ncode/regint.f

### intgr.omp.f	{257}
* (mode B) add the case of the tidal tensor in mode B, as in Ncode/intgr.f

### start_omp.F	{42}
* call ttinit if KZ(14) = 9 as in Ncode/start.f
