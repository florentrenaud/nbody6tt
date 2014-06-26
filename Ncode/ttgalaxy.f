      SUBROUTINE TTGALAXY(X,Y,Z,T,RSCALE,MSCALE,VSCALE,TSCALE,TTPHIG)
*
*       User-defined galactic potential
*       -------------------------------
*** FlorentR - new subroutine

* Compute the galactic potential as a function of the position X, Y, Z and
* time T.
*
* All numerical value should be given in the nbody6 units:
* To convert a value in physical units into nbody6 units, *divide* the
* physical quantity by: RSCALE,MSCALE,VSCALE,TSCALE for the length,
* mass, velocity and time respectively.
* Scale your definition of the potential to G = 1.

      REAL*8 X,Y,Z,T, RSCALE, MSCALE, VSCALE, TSCALE, TTPHIG
      LOGICAL FIRST
      
***vvvvvvvvvvvvvvvvvvvvvv Modify below vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv***

* Declare your variables and constants
      REAL*8 R,R2,Rcirc,Rcirc2
      REAL*8 MASS

* Declare which variables are constants (i.e. to be computed only once)
      SAVE MASS

***^^^^^^^^^^^^^^^^^^^^^^ Modify above ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^***

* Do not change this part.
      SAVE FIRST
      DATA FIRST /.TRUE./
      IF(FIRST) THEN

***vvvvvvvvvvvvvvvvvvvvvv Modify below vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv***

* Compute your constants
        MASS = (-1e10 / MSCALE)
	
***^^^^^^^^^^^^^^^^^^^^^^ Modify above ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^***

* Do not change this part.
        FIRST = .FALSE.
      ENDIF

***vvvvvvvvvvvvvvvvvvvvvv Modify below vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv***

* Compute your variables
!      Rcirc2 = X**2+Y**2
!      Rcirc = SQRT(Rcirc2)
!      R2 = Rcirc2+Z**2
!      R = SQRT(R2)
      R = SQRT(X**2+Y**2+Z**2)

* Compute the value of the galactic potential TTPHIG
*      TTPHIG = 0.0D0 ! no galaxy
      TTPHIG = MASS / R ! point-mass      

! put some examples here

***^^^^^^^^^^^^^^^^^^^^^^ Modify above ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^***

* Do not change this part.
      RETURN
      END
