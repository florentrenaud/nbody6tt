Version 3.0	- 6 June 2012

Based on Nbody6 downloaded on the 5 June 2012 from http://www.ast.cam.ac.uk/~sverre/web/pages/nbody.htm (files nbody6.tar.gz and gpu2.tar.gz). Bug fixes and patches from Sverre since this date are NOT taken into account.


List of changes made to Sverre's Nbody6. Discard this document if you don't want to dig in the code.
All the changes are commented in the .f and .h files. They start with the tag '*** FlorentR' and end with '*** FRenaud'.
The number(s) in curly brackets indicates the line number(s) where change(s) have(s) been made.

For future updates: make sure that KZ(14)=9 is considered in xtrnld.f (for the moment, it is with KZ(14).GE.3 on line 41)


In Ncode
========

### Makefile
* add the ttcal.f and ttinit.f files to the source list
* add a target all to make cpu version (not related to tides)

### adjust.f	{241}
* set the tidal radius to 10x the half-mass radius

### common6.h	{76}
* add the TT block that contains tidal tensor data

### escape.f	{157}
* include the case of tidal tensor when setting RTIDE

### jacobi.f	{8,13,18,20,39,42,44}
* (only used if KZ(21)>1)
* use of the bound mass as the cluster mass (also for KZ(14).NE.9)
* use of potential+kinetic energy as membership criterion

### mydump.f	{14,45,92,95,98,146}
* add the TT block to the saving and restoring procedures.
* use the restart.dat file. (not related to tides)

### modify.f	{75}
* call ttinit if the tensor file needs to be re-read

### nbody6.f	{31}
* print a little header to the output, in order to identify the version used to produce the simulation.

### param.h {12}
* add the nbttmax parameter (this implies to modify the previous line too)

### regint.f	{136}
* include the case of the tidal tensor

### start.f {42}
* call ttinit if KZ(14) = 9

### ttcal.f (new file)
* Computes the tidal tensor at current time by means of interpolation.

### ttinit.f	(new file)
* force the center of mass correction (KZ(31)=1)
* read the tensor data from file 'tt.dat', scale and initialize

### xtrnlf.f	{28}
* include the computation of the external force from the tidal tensor

### xtrnlv.f	{84}
* include the potential energy from the tidal tensor



In GPU2
=====

### Makefile_gpu
* add the ttcal.f and ttinit.f files to the source list

### adjust.f	{262}
* set the tidal radius to 10x the half-mass radius as in Ncode/adjust.f

### start_omp.F	{42}
* call ttinit if KZ(14) = 9 as in Ncode/start.f

