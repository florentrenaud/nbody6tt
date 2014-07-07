      SUBROUTINE TTFORCE(XI,XIDOT,FM,FD)
*
*
*       Compute the galactic force and its time derivative (mode B)
*       -----------------------------------------------------------
*
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      
      REAL*8 XI(3),XIDOT(3),FM(3),FD(3)
      REAL*8 TT2DX, TT12DX, TT4DX2, TTTMP
      REAL*8 TTX(75), TTPHI(25)

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


      TTDX = (2.2D-16)**0.25D0 * 0.33D0
     &  * SQRT((XI(1)**2 + XI(2)**2 + XI(3)**2)

      TT2DX = 2.0 * TTDX
      TT12DX = 6.0 * TT2DX
      TT4DX2 = TT2DX**2
      DO I=1,75
        TTX(I) = TTX(I)*TTDX
      END DO

* Get the galactic potential around the position XI using
* user-defined formula in ttgalaxy.f

      DO I=1,25
         CALL TTGALAXY( XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TG, RBAR, ZMBAR, VSTAR, TSTAR, TTTMP)
         TTPHI(I) = TTTMP
      END DO

* Get the galactic force at the position XI:   F = -grad phi
* by second-order difference using the points 1,4,7,9,11,12,14,15,17,19,22,25
      FM(1) = (TTPHI(15)-8.*TTPHI(14)+8.*TTPHI(12)-TTPHI(11))/TT12DX
      FM(2) = (TTPHI(19)-8.*TTPHI(17)+8.*TTPHI(9) -TTPHI(7) )/TT12DX
      FM(3) = (TTPHI(25)-8.*TTPHI(22)+8.*TTPHI(4) -TTPHI(1) )/TT12DX

*  first-order difference using the points 4,9,12,14,17,22
*      FM(1) = ( TTPHI(12) - TTPHI(14) ) / TT2DX
*      FM(2) = ( TTPHI(9)  - TTPHI(17) ) / TT2DX
*      FM(3) = ( TTPHI(4)  - TTPHI(22) ) / TT2DX


* Get the tidal tensor at XI
* Note that, contrary to mode A,
* the tidal tensor is not used to compute the tidal force in mode B.

*  first-order difference on the force at the points 4,9,12,14,17,22
      TTENS(1,1,1) = (TTPHI(13)-TTPHI(15)-TTPHI(11)+TTPHI(13)) / TT4DX2
      TTENS(1,2,1) = (TTPHI(16)-TTPHI(18)-TTPHI(8)+TTPHI(10)) / TT4DX2
      TTENS(1,3,1) = (TTPHI(21)-TTPHI(23)-TTPHI(3)+TTPHI(5)) / TT4DX2
      TTENS(2,1,1) = TTENS(1,2,1)
      TTENS(2,2,1) = (TTPHI(13)-TTPHI(19)-TTPHI(7)+TTPHI(13)) / TT4DX2
      TTENS(2,3,1) = (TTPHI(20)-TTPHI(24)-TTPHI(2)+TTPHI(6)) / TT4DX2
      TTENS(3,1,1) = TTENS(1,3,1)
      TTENS(3,2,1) = TTENS(2,3,1)
      TTENS(3,3,1) = (TTPHI(13)-TTPHI(25)-TTPHI(1)+TTPHI(13)) / TT4DX2

* Get the time derivative of the tidal force:
* dF/dt = dF/dx * dx/dt = T * v
* where T is the tidal tensor
* In 3D: dF_x/dt = dF_x/dx v_x + dF_x/dy v_y + dF_x/dz v_z
*                = T_{x,x} v_x + T_{x,y} v_y + T_{x,z} v_z

      FD(1) = TTENS(1,1,1) * XIDOT(1) + TTENS(1,2,1) * XIDOT(2) 
     &   + TTENS(1,3,1) * XIDOT(3)
      FD(2) = TTENS(2,1,1) * XIDOT(1) + TTENS(2,2,1) * XIDOT(2)
     &   + TTENS(2,3,1) * XIDOT(3)
      FD(3) = TTENS(3,1,1) * XIDOT(1) + TTENS(3,2,1) * XIDOT(2)
     &   + TTENS(3,3,1) * XIDOT(3)

      RETURN
      END
          
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
