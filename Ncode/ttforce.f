      SUBROUTINE TTFORCE(XI,FM,FD,DT)
*
*
*       Compute the force and the tidal tensor (mode B)
*       -----------------------------------------------
*
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      
      REAL*8 XI(3),FM(3),FD(3),DT
      REAL*8 TT2DX, TT12DX, TT4DX2, TTTMP
      REAL*8 TTX(75), TTPHI(25)
      REAL*8 OLDFM(3)
      LOGICAL FIRST
      
      SAVE OLDFM, TTX, TT2DX, TT12DX, TT4DX2, FIRST
      DATA FIRST /.TRUE./
* init

      IF(FIRST) THEN 
        TT2DX = 2.0 * TTDX
        TT12DX = 6.0 * TT2DX
        TT4DX2 = TT2DX**2

* Evaluate the potential at the position of the cluster, plus
* at 6 points around the cluster and 6 points around these 6 points
* = 25 points (some are counted twice).

* TTX contains the position of the 25 points, relative to the cluster.
* Points are ordered by increasing z, then y, then x. (See below)
* The cluster is number 13.
* TTX(1:25)  = x positions
* TTX(26:50) = y positions
* TTX(51:75) = z positions

        TTX = (/ 0, 
     &             0, -1, 0, 1, 0,
     &             0, -1, 0, 1, -2, -1, 0, 1, 2, -1, 0, 1, 0,
     &             0, -1, 0, 1, 0,
     &             0,
     &           0,
     &             -1, 0, 0, 0, 1,
     &             -2, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 2,
     &             -1, 0, 0, 0, 1,
     &             0,
     &           -2,
     &             -1, -1, -1, -1, -1,
     &             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &             1, 1, 1, 1, 1,
     &             2 /)

        DO I=1,75
          TTX(I) = TTX(I)*TTDX
        END DO

        OLDFM(1) = 0.0
        OLDFM(2) = 0.0
        OLDFM(3) = 0.0

        TTTIME(1) = 0.0
        TTTIME(2) = 0.0

      ELSE

        TTTIME(1) = TTTIME(2)
        TTTIME(2) = TTTIME(2) + DT
        
      END IF
     
      DO I=1,25
        CALL TTGALAXY( XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TTTIME(2), RBAR, ZMBAR, VSTAR, TSTAR, TTTMP)
        TTPHI(I) = TTTMP
      END DO

* force on the cluster:
*  second-order difference using the points 4,9,12,14,17,22
* (minus sign because F = -grad phi)
      FM(1) = (TTPHI(15)-8.*TTPHI(14)+8.*TTPHI(12)-TTPHI(11))/TT12DX
      FM(2) = (TTPHI(19)-8.*TTPHI(17)+8.*TTPHI(9) -TTPHI(7) )/TT12DX
      FM(3) = (TTPHI(25)-8.*TTPHI(22)+8.*TTPHI(4) -TTPHI(1) )/TT12DX

*  first-order difference using the points 4,9,12,14,17,22
*      FM(1) = ( TTPHI(12) - TTPHI(14) ) / TT2DX
*      FM(2) = ( TTPHI(9)  - TTPHI(17) ) / TT2DX
*      FM(3) = ( TTPHI(4)  - TTPHI(22) ) / TT2DX

* time derivative of the force:
      IF(FIRST) THEN
        FD(1) = 0.0
        FD(2) = 0.0
        FD(3) = 0.0
        FIRST=.FALSE.      
      ELSE
        FD(1) = ( FM(1) - OLDFM(1) ) / DT
        FD(2) = ( FM(2) - OLDFM(2) ) / DT
        FD(3) = ( FM(3) - OLDFM(3) ) / DT
      ENDIF

      OLDFM(1) = FM(1)
      OLDFM(2) = FM(2)
      OLDFM(3) = FM(3)

* tidal tensor on the cluster:
*  first-order difference on the force at the points 4,9,12,14,17,22

      DO I=1,3
        DO J=1,3
          TTENS(I,J,1) = TTENS(I,J,2)
        END DO
      END DO
	
      TTENS(1,1,2) = (TTPHI(13)-TTPHI(15)-TTPHI(11)+TTPHI(13)) / TT4DX2
      TTENS(1,2,2) = (TTPHI(16)-TTPHI(18)-TTPHI(8)+TTPHI(10)) / TT4DX2
      TTENS(1,3,2) = (TTPHI(21)-TTPHI(23)-TTPHI(3)+TTPHI(5)) / TT4DX2
      TTENS(2,1,2) = TTENS(1,2,2)
      TTENS(2,2,2) = (TTPHI(13)-TTPHI(19)-TTPHI(7)+TTPHI(13)) / TT4DX2
      TTENS(2,3,2) = (TTPHI(20)-TTPHI(24)-TTPHI(2)+TTPHI(6)) / TT4DX2
      TTENS(3,1,2) = TTENS(1,3,2)
      TTENS(3,2,2) = TTENS(2,3,2)
      TTENS(3,3,2) = (TTPHI(13)-TTPHI(25)-TTPHI(1)+TTPHI(13)) / TT4DX2

* if output is required, do it here
*      WRITE(666,*) TTTIME(2), (XI(I),I=1,3)
*      CALL FLUSH(666)

*      WRITE(667,*) TTTIME(2), DT 
*      CALL FLUSH(667)

*      WRITE(668,*) TTTIME(2), ((TTENS(I,J,2),I=1,3), J=1,3)
*      WRITE(668,*) TTTIME(2), TTENS(1,2,2), TTENS(2,1,2)
*      CALL FLUSH(668)

      RETURN
      END
*          
*  Positions of the TTX points, wrt the cluster (number 13)
*
*                          TTDX
*                       <-------->
*
*                      +--------+--------+--------+--------+       ^
*                     /        /        /        /        /       / TTDX
*                    /        /        /        /        /       /
*                   +--------+--------+--------+--------+       V
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +--------+------- 25 ------+--------+     z = +2 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+--------+--------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*          
*          
*                      +--------+--------+--------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +--------+------- 24 ------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +------- 21 ----- 22 ----- 23 ------+     z = +1 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+------- 20 ------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*          
*          
*                      +--------+------- 19 ------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +------- 16 ----- 17 ----- 18 ------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                11 ----- 12 ----- 13 ----- 14 ---- 15     z = 0
*               /        /        /        /        /
*              /        /        /        /        /
*             +------- 8 ------ 9 ------ 10 ------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+------- 7 -------+--------+
*          
*          
*                      +--------+--------+--------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +--------+------- 6 -------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +------- 3 ------ 4 ------ 5 -------+     z = -1 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+------- 2 -------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*          
*          
*                      +--------+--------+--------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +--------+--------+--------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +--------+------- 1 -------+--------+     z = -2 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+--------+--------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*
*     Force at 4
*        Fx = ( TTPHI(3)  - TTPHI(5)  ) / TT2DX
*        Fy = ( TTPHI(2)  - TTPHI(6)  ) / TT2DX
*        Fz = ( TTPHI(1)  - TTPHI(13) ) / TT2DX
*     Force at 9
*        Fx = ( TTPHI(8) - TTPHI(10)  ) / TT2DX
*        Fy = ( TTPHI(7) - TTPHI(13)  ) / TT2DX
*        Fz = ( TTPHI(2) - TTPHI(20)  ) / TT2DX
*     Force at 12
*        Fx = ( TTPHI(11) - TTPHI(13) ) / TT2DX
*        Fy = ( TTPHI(8)  - TTPHI(16) ) / TT2DX
*        Fz = ( TTPHI(3)  - TTPHI(21) ) / TT2DX
*     Force at 14
*        Fx = ( TTPHI(13) - TTPHI(15) ) / TT2DX
*        Fy = ( TTPHI(10) - TTPHI(18) ) / TT2DX
*        Fz = ( TTPHI(5)  - TTPHI(23) ) / TT2DX
*     Force at 17
*        Fx = ( TTPHI(16) - TTPHI(18) ) / TT2DX
*        Fy = ( TTPHI(13) - TTPHI(19) ) / TT2DX
*        Fz = ( TTPHI(6)  - TTPHI(24) ) / TT2DX
*     Force at 22
*        Fx = ( TTPHI(21) - TTPHI(23) ) / TT2DX
*        Fy = ( TTPHI(20) - TTPHI(24) ) / TT2DX
*        Fz = ( TTPHI(13) - TTPHI(25) ) / TT2DX
*
