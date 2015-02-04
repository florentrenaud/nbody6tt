nbody6tt
========

Version 5.4.2

nbody6tt is based on Sverre Aarseth's nbody6, including the GPU2 package. The original version has been downloaded on the 29 January 2015 from http://www.ast.cam.ac.uk/~sverre/web/pages/nbody.htm (files nbody6.tar.gz and gpu2.tar.gz). Bug fixes and patches from Sverre since this date are NOT taken into account.


Description
===========

nbody6tt includes the treatment of external tidal fields in nbody6. Theory and method are detailed in Renaud, Gieles & Boily (2011) and Renaud & Gieles (2015).

Since v3.0, nbody6tt offers two modes:
* (A) a table a tidal tensors computed separatly is provided to the code. Mode A allows arbitrarily complex configurations (e.g. galaxy mergers). Tidal tensors are usually extracted from galaxy simulations and thus are known at a frequency larger than the smaller timestep. Therefore their values are interpolated in time (quadratic interpolation) at the times required by the simulation. Before the timestamp of the second tensor is reached, the interpolation is linear.
* (B) the user defines a potential and initial position/velocity for the cluster. The orbit of the cluster (its guiding center, to be exact) and the corresponding tidal forces on stars are computed 'live' by nbody6tt. Only the analytical form of the galactic potential is required. The tidal force is computed as the difference between the force at the position of a star, and that at the position of the guiding center of the cluster. Tensors are not used to compute the force. However, tensors are used to get the time derivative of the tidal force (see Ncode/ttforce.f).

An IMPORTANT note on the tidal radius: in complex tidal fields, it is very involved to compute the centrifugal and Euler effects that must be taken into account when estimating the tidal radius. When using the tidal tensor, the tidal radius of nbody6tt is only a scalelength used to set the stripping radius. To be clear: RTIDES has no physical meaning in nbody6tt.



User guide
==========

Before using nbody6tt, it is STRONGLY recommended to get familiar with nbody6 and make sure it installs and run on your machine. If nbody6tt does not compile or crashes, it is very likely that you would have had the same problem with nbody6. See the help of nbody6 for solving the problem.

The tidal tensor version is only supposed to cause you science trouble!


How to install the CPU version
------------------------------
How to install the CPU version:
(for the GPU version, skip this part, see below)

Get the nbody6tt directory, including the subdirectories
	Chain[, GPU], Nchain, Ncode
Link common6.h and params.h in the Nchain directory to common6.h and params.h in the Ncode directory:
	cd Nchain
	ln -s ../Ncode/common6.h common6.h 
	ln -s ../Ncode/params.h params.h    
	cd ..
Check that Ncode/params.h contains values that matches your requirements (e.g. NMAX).
Go to the Ncode directory:
	cd Ncode
Edit Makefile by setting the compiler options and paths.
Run the makefile:
	make
This will create the executable nbody6, that can be used the same way than nbody6, i.e.:
	nbody6 < input > output


How to install the GPU version
------------------------------
(for the CPU version, skip this part, see above)

Get the nbody6tt directory, including the subdirectories
	Chain, GPU2, Nchain, Ncode
Link common6.h and params.h in the GPU2 and Nchain directories to common6.h and params.h in the Ncode directory:
	cd GPU2
	ln -s ../Ncode/common6.h common6.h 
	ln -s ../Ncode/params.h params.h    
	cd ../Nchain
	ln -s ../Ncode/common6.h common6.h 
	ln -s ../Ncode/params.h params.h
	cd ..
Check that Ncode/params.h contains values that matches your requirement (e.g. NMAX).
Go to the GPU2 directory:
	cd GPU2
Edit Makefile, Makefile_gpu and ../Ncode/Makefile by setting the compiler options and paths.
Run the makefile with the gpu target:
	make gpu
This will create the executable nbody6.gpu, that can be used the same way than nbody6, i.e.:
	nbody6.gpu < input > output



How to run nbody6tt in mode A (CPU or GPU)
--------------------------------

(A1) Create a tensor file named 'tt.dat' in the running directory: its syntax is as follow:
	NBTT TTUNIT TTOFFSET
	TTTIME1 TTENS1(9)
	TTTIME2 TTENS2(9)
	TTTIME3 TTENS3(9)
	...
where NBTT is the number of tensors given in the file (i.e. the total number of rows in 'tt.dat' minus one), TTUNIT is the timescale for the galactic run in Myr (= how many Myr corresponds to t=1), TTOFFSET is an offset (in Myr) to be added to the timestamps of the tensors (it can be positive, 0.0, or negative). TTTIMEx is a timestamp in galactic run units, TTENSx(9) are the 9 components of the tensor (although only 6 are used), in the unit system of the galactic run (i.e. time^-2). The lines of the tensors must be ordered chronologically, with no duplicates.

Therefore, to convert the times of the tensors into Nbody6 units, Nbody6tt computes (in Ncode/ttinit.f)
    TTTIME(K)*(TTUNIT/TSTAR)+TTOFFSET/TSTAR
(with TSTAR being the scaling factor of nbody units into Myr).

Important: make sure that NBTT is smaller than NBTTMAX in Ncode/params.h.

(A2) Set the option KZ(14) to 9 in the input to use the tidal tensor treatment.

(A3) Run the code:
	nbody6[.gpu] < input > output

Important: the code will stop when the tidal tensor is not defined (before the first TTTIME*TTUNIT and after the last TTTIME*TTUNIT).

Important: when setting tensors manually, be sure to have a sufficient number of tensors in your tt.dat file. (See the interpolation method in Ncode/ttcal.f)

Note that setting KZ(14)=0 is not strictly equivalent to having KZ(14)=9 with a null tensor, because of a different behavior of the code (mainly for stripping escapers). To compare the evolution of a cluster in a given tidal field with that of an isolated cluster, it is safer to make both runs with KZ(14)=9.

The no-tides 'tt.dat' file is:
	2 100000.0 0.0
	0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0



How to run nbody6tt in mode B (CPU or GPU)
--------------------------------

(B1) Edit the 'Ncode/ttgalaxy.f' file: create a fortran function that returns the value of the galactic potential at the position x,y,z of the cluster center, and at the time t. Recompile nbody6tt:

CPU version:
	cd Ncode ; make

GPU version:
	cd GPU2 ; make GPU

("make clean" should not be required.)

This step (B1) is to be done each time that the galactic potential is changed.
   
(B2) In the input file, set option KZ(14) to 9, and the RG(1:3) (initial position of the cluster in kpc) and VG(1:3) (initial velocity of the cluster in km/s) values, as you would do in the usual KZ(14)=3 mode. Avoid setting a initial position strictly equal to 0,0,0. (A small offset would be sufficient.)

(B3) Make sure that your running directory does NOT contains a file called 'tt.dat'. If a 'tt.dat' file is found, the code will run in Mode A. Run the code:
	nbody6[.gpu] < input > output

Note that setting KZ(14)=0 is not strictly equivalent to having KZ(14)=9 with a null potential, nor to KZ(14)=3 with the same potential and same orbit. This is because of a different behavior of the code (mainly for stripping escapers). To compare the evolution of a cluster in a given tidal field with that of an isolated cluster, it is safer to make both runs with KZ(14)=9. The no-tides 'ttgalaxy' function is:
	TTPHIG = 0.0


How to restart nbody6tt (A or B mode)
-------------------------------------

Whether the code crashed or had been killed (runtime limit, killed by user...), you should have a file named 'restart.tmp' in the running directory. In the very unlikely case where 'restart.tmp' does not exist, you should have a 'fort.2' file. In the extremely unlikely case where none exist, restarting is not possible.

(1) You probably want to restart in a new directory, otherwise your previous data will be lost.

(2) Rename 'restart.tmp' (if you don't have it, use 'fort.2') as 'restart.dat' in the new directory.

(3a) To restart by keeping the parameters and the tensors as they were initially, create the input file with the usual syntax, i.e.
	2 TCOMP
where TCOMP is the maximum CPU time in minutes.

(3b) To restart with a modification of some parameter (e.g. by reading a new tensor file (A) or a new position/velocity for the cluster (B)), create the input file with the usual syntax, i.e.:
	KSTART TCOMP                                (KSTART = 3, 4, 5)
	DTADJ DELTAT TADJ TNEXT TCRIT QE J KZ(J)    (if > 0 & KSTART = 3 or 5).
	ETAI ETAR ETAU DTMIN RMIN NCRIT             (if > 0 & KSTART = 4 or 5).   

Setting KSTART= 3 or 5, J=14 and KZ(J)=9 will make the code read the tensor file again and the interpolation scheme will be reset (i.e. the first step is linear again).
Example: to read a new table of tensors:
	3 1000000000
	0.0 0.0 0.0 0.0 0.0 0.0 14 9

(3c) In mode B, you can restart with a new galactic potential: define it in 'ttgalaxy.f', compile it, and restart normally (see point 3a, above, to use the last position and velocity stored in restart.tmp, or see point 3b, also above, with a new setting).

Restarting from a file created with another version of the code is not possible (in principle).
