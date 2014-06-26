      SUBROUTINE TTTAIL(XI,FM,FMD)
*
*
*       Compute the gravitational force on tidal tail (mode B)
*       ------------------------------------------------------
*
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      
      REAL*8 XI(3),FM(3),FMD(3)
      REAL*8 TT12DX, TTTMP
      REAL*8 TTX(36), TTPHI(12)
      LOGICAL FIRST
      
      SAVE TTX, TT12DX, FIRST
      DATA FIRST /.TRUE./
* init

      IF(FIRST) THEN 
        TT12DX = 12.0 * TTDX

* see ttforce.f for details
* (only points 1 4 7 9 11 12 14 15 17 19 22 25 are used, but renamed here)
       TTX = (/ 0,0,0,0,-2,-1,1,2,0,0,0,0, 
     &          0,0,-2,-1,0,0,0,0,1,2,0,0,
     &          -2,-1,0,0,0,0,0,0,0,0,1,2 /)

        DO I=1,36
          TTX(I) = TTX(I)*TTDX
        END DO

C* Reset the array of old force on tail members
C      DO I=1,NMAX
C        OLDFM(1,K) = 0.0
C        OLDFM(2,K) = 0.0
C        OLDFM(3,K) = 0.0
C      END DO

      END IF

* compute potential using the same time as for the tidal tensor
      DO I=1,12
        CALL TTGALAXY( XI(1)+TTX(I), XI(2)+TTX(I+12), XI(3)+TTX(I+24),
     &        TTTIME(2), RBAR, ZMBAR, VSTAR, TSTAR, TTTMP)
        TTPHI(I) = TTTMP
      END DO

* force on the cluster:
*  second-order difference
* (minus sign because F = -grad phi)
      FM(1) = (TTPHI(8)-8.*TTPHI(7)+8.*TTPHI(6)-TTPHI(5))/TT12DX
      FM(2) = (TTPHI(10)-8.*TTPHI(9)+8.*TTPHI(4) -TTPHI(3) )/TT12DX
      FM(3) = (TTPHI(12)-8.*TTPHI(11)+8.*TTPHI(2) -TTPHI(1) )/TT12DX

* for now, we neglect the time-derivative of the force
      FMD(1) = 0.0
      FM(2) = 0.0
      FM(3) = 0.0

      RETURN
      END
