
Based on Nbody6 downloaded on 29 January 2015.

List of changes made to Sverre's Nbody6. Discard this document if you don't want to dig in the code.
All the changes are commented in the .f and .h files. They start with the tag '*** FlorentR' and end with '*** FRenaud'.

In Ncode
========

### Makefile
* add the ttcal.f, ttinit.f, ttforce.f, ttgalaxy.f files to the source list
* add a target all to make cpu version (not related to tides)

### adjust.f
* set the tidal radius to 10x the half-mass radius

### common6.h
* add the TT block that contains tidal tensor data
* add the boolean TTMODE

### define.f
* (mode B) add the reading of the initial position and velocity of the cluster

### escape.f
* include the case of tidal tensor when setting RTIDE

### gcint.f
* (mode B) call ttforce to get the force from the galaxy on the cluster

### intgrt.f
* (mode B) add the case of the tidal tensor in mode B

### jacobi.f
* (only used if KZ(21)>1)
* use of the bound mass as the cluster mass (also for KZ(14).NE.9)
* use of potential+kinetic energy as membership criterion

### mydump.f
* add the TT block to the saving and restoring procedures.
* use the restart.dat file. (not related to tides)
* add the TTMODE boolean

### modify.f
* call ttinit if the tensor file needs to be re-read

### nbody6.f
* print a little header to the output, in order to identify the version used to produce the simulation.

### output.f
* (mode B) output the orbit

### param.h
* add the nbttmax parameter (this implies to modify the previous line too)

### regint.f
* include the case of the tidal tensor

### start.f
* call ttinit if KZ(14) = 9

### ttcal.f (new file)
* (mode A) computes the tidal tensor at current time by means of interpolation.

### ttforce.f	(new file)
* (mode B) computes the force and its time derivative (using tidal tensor) at the position of the stars or the cluster, using the user-defintion of the galactic potential

### ttgalaxy.f	(new file)
* (mode B) user-definition of the galactic potential

### ttinit.f	(new file)
* determine the TT mode (A or B)
* force the center of mass correction (KZ(31)=1)
* (mode A) read the tensor data from file 'tt.dat', scale and initialize
* (mode B) read the initial coordinates of the cluster and initialize

### xtrnlf.f
* include the computation of the external force from the tidal tensor

### xtrnlv.f
* include the potential energy from the tidal tensor (only at setup)

### xtrnlt.f
* (mode B) Compute the galactic force on tail members using the user-definition of the potential



In GPU2
=====

### Makefile.build
* add the ttcal.f, ttinit.f, ttforce.f, ttgalaxy.f files to the source list

### adjust.f
* set the tidal radius to 10x the half-mass radius as in Ncode/adjust.f

### gpucor.f
* include the case of the tidal tensor, as in Ncode/regint.f

### intgr.omp.f
* (mode B) add the case of the tidal tensor in mode B, as in Ncode/intgr.f

### start_omp.F
* call ttinit if KZ(14) = 9 as in Ncode/start.f
