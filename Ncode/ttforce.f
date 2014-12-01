      SUBROUTINE TTFORCE(XI,XIDOT,FM,FD,TTGET)
*
*
*       Compute the galactic force and its time derivative (mode B)
*       -----------------------------------------------------------
*
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)

      INTEGER TTTDEP, TTGET, TTBUILD
      REAL*8 XI(3),XIDOT(3),FM(3),FD(3), TTTMP

      REAL*8 TTCX(75), TTCDX(4), TTC2DX(4), TTC12DX(4)
      REAL*8 TTX(75), TTDX(4), TT2DX(4), TT12DX(4), TTDXS

      REAL*8 TTPHI(25), TTPHIT(6)


      SAVE TTCX, TTCDX, TTC2DX, TTC12DX, TTDXS, TTCOUNT


* evaluate potential at current position
      CALL TTGALAXY(XI(1), XI(2), XI(3), TG,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
      TTPHI(1) = TTTMP

* Restore the stencil used for the guiding centre (previously computed)
      IF(TTGET .EQ. 0) THEN
        DO K=1,4
          TTDX(K) = TTCDX(K)
          TT2DX(K) = TTC2DX(K)
          TT12DX(K) = TTC12DX(K)
        ENDDO
        DO K=1,75
         TTX(K) = TTCX(K)
        ENDDO
      ENDIF


* treat time-dependence, if necessary (only for guiding centre)
      IF(TTTDEP .GT. 0 .AND. TTGET .GT. 0) THEN

* define initial guess for the step size to be big.
        IF(TG .EQ. 0) THEN
          TTDX(4) = 1.0D0
        ELSE
          TTDX(4) = 1.1D-2 * ABS(TG)
        ENDIF

        CALL TTSTEP(XI, TTDX, 4, TTPHI(1), 3D0)
        TT2DX(4) = 2.0 * TTDX(4)
      ENDIF


      TTBUILD = 0
      DO K=1,3
* compute step size 
* for guiding centre and stars far from guiding centre
        IF(TTGET .GT. 0 .OR. XI(K) .GT. 10.0*TTDXS(K)) THEN
          TTCOUNT = TTCOUNT + 1
* ste flog to rebuild the stencil
          TTBUILD = 1

* define initial guess for the step size to be big.
          IF(XI(K) .EQ. 0) THEN
            TTDX(K) = 1.0D0
          ELSE
            TTDX(K) = 1.1D2 * ABS(XI(K))
          ENDIF

          CALL TTSTEP(XI, TTDX, K, TTPHI(1), 5D0)
          TT2DX(K) = 2.0 * TTDX(K)
          TT12DX(K) = 6.0 * TT2DX(K)

        ENDIF
        
      ENDDO

      
* Build the stencil using the value of TTDX
      IF(TTBUILD .GT. 0) THEN

* TTX contains the position of the 25 points, relative to XI
* The position XI is number 1
* TTX(1:25)  = x positions
* TTX(26:50) = y positions
* TTX(51:75) = z positions

        TTX = (/ 0, 
     &        0, 0,-1, 1, 0, 0,
     &        0, 0,-1, 1, 0, 0,-1, 1,-2, 2,-1, 1, 0, 0,-1, 1, 0, 0,
     &        0, 
     &        0,-1, 0, 0, 1, 0,
     &        0,-1, 0, 0, 1,-2,-1,-1, 0, 0, 1, 1, 2,-1, 0, 0, 1, 0,
     &        0,
     &        -1, 0, 0, 0, 0, 1,
     &        -2,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2 /)
        
        DO I=1,25
          TTX(I) = TTX(I)*TTDX(1)
          TTX(I+25) = TTX(I+25)*TTDX(2)
          TTX(I+50) = TTX(I+50)*TTDX(3)
        END DO

* save the step size and stencil used for the guiding centre
        IF(TTGET .EQ. 1) THEN
          write(666,*) TTCOUNT
          TTCOUNT = 0
          DO K=1,4
            TTCDX(K) = TTDX(K)
            TTC2DX(K) = TT2DX(K)
            TTC12DX(K) = TT12DX(K)
            TTDXS(K) = TTDX(K)**5 / 2.220446D-16
          ENDDO
          DO K=1,75
            TTCX(K) = TTX(K)
          ENDDO           
        ENDIF
      ENDIF
* at this point, the step sizes and the stencil have the values that
* will be used in the computation.



* Get the galactic potential around the position XI
      DO I=2,25
         CALL TTGALAXY( XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TG, RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
         TTPHI(I) = TTTMP
      END DO

* Get the galactic force at XI:
      FM(1) = (TTPHI(17)-8.*TTPHI(5)+8.*TTPHI(4)-TTPHI(16))/TT12DX(1)
      FM(2) = (TTPHI(20)-8.*TTPHI(6)+8.*TTPHI(3)-TTPHI(13))/TT12DX(2)
      FM(3) = (TTPHI(25)-8.*TTPHI(7)+8.*TTPHI(2)-TTPHI(8))/TT12DX(3)


* Get the time derivative of the tidal force:
* dF/dt = partial F / partial x * partial x / partial t
*             + partial F / partial t
*       = T * v + partial F / partial t
* where T is the tidal tensor. 
* The last term (partial F / partial t) is 0 for constant potentials 
* and can be neglected for adiabatic changes of the potential.
* This saves a lot of CPU time.
*
* In 3D: dF_x/dt = partial F_x / partial x v_x
*                  + partial F_x / partial y v_y
*                  + partial F_x / partial z v_z
*                  + partial F_x / partial t
*                = T_{x,x} v_x + T_{x,y} v_y + T_{x,z} v_z
*                  + partial F_x / partial t


* Get the tidal tensor at XI:
      TTENS(1,1,1) = (TTPHI(1)-TTPHI(17)-TTPHI(16)+TTPHI(1)) 
     &        / (TT2DX(1)*TT2DX(1))
      TTENS(1,2,1) = (TTPHI(18)-TTPHI(19)-TTPHI(14)+TTPHI(15))
     &        / (TT2DX(1)*TT2DX(2))
      TTENS(1,3,1) = (TTPHI(22)-TTPHI(23)-TTPHI(10)+TTPHI(11))
     &        / (TT2DX(1)*TT2DX(3))
      TTENS(2,1,1) = TTENS(1,2,1)
      TTENS(2,2,1) = (TTPHI(1)-TTPHI(20)-TTPHI(13)+TTPHI(1))
     &        / (TT2DX(2)*TT2DX(2))
      TTENS(2,3,1) = (TTPHI(21)-TTPHI(24)-TTPHI(9)+TTPHI(12))
     &        / (TT2DX(2)*TT2DX(3))
      TTENS(3,1,1) = TTENS(1,3,1)
      TTENS(3,2,1) = TTENS(2,3,1)
      TTENS(3,3,1) = (TTPHI(1)-TTPHI(25)-TTPHI(8)+TTPHI(1))
     &        / (TT2DX(3)*TT2DX(3))


      FD(1) = TTENS(1,1,1) * XIDOT(1) + TTENS(1,2,1) * XIDOT(2) 
     &   + TTENS(1,3,1) * XIDOT(3)
      FD(2) = TTENS(2,1,1) * XIDOT(1) + TTENS(2,2,1) * XIDOT(2)
     &   + TTENS(2,3,1) * XIDOT(3)
      FD(3) = TTENS(3,1,1) * XIDOT(1) + TTENS(3,2,1) * XIDOT(2)
     &   + TTENS(3,3,1) * XIDOT(3)



* if the potential is time dependent
      IF(TTTDEP .GT. 0) THEN

* compute partial F / partial t at t+dt and t-dt.
        DO I=2,7
           CALL TTGALAXY(XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TG-TTDX(4), RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
           TTPHI(I) = TTTMP

           CALL TTGALAXY(XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TG+TTDX(4), RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
           TTPHI(I+6) = TTTMP
        END DO

* sum T*v and partial F / partial t
        FD(1) = FD(1)+(-TTPHI(4)+TTPHI(5)+TTPHI(10)-TTPHI(11)) 
     &        / (TT2DX(4)*TT2DX(1))
        FD(2) = FD(2)+(-TTPHI(3)+TTPHI(6)+TTPHI(9)-TTPHI(12))
     &        / (TT2DX(4)*TT2DX(2))
        FD(3) = FD(3)+(-TTPHI(2)+TTPHI(7)+TTPHI(8)-TTPHI(13))
     &        / (TT2DX(4)*TT2DX(3))
      ENDIF

      RETURN
      END
