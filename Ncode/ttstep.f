      SUBROUTINE TTSTEP(XI, TTH, TTAXIS, TTP, TTORDER)
*
*
*       Compute the step size used in the evaluation of the galactic 
*       acceleration and jerk (mode B)
*       -----------------------------------------------------------
*
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)

      INTEGER NTAB, J,K,L,TTAXIS,TTTDEP, BESTJ, BESTK, TTDIR
      REAL*8 TTORDER, TTH(4), TTHH, XI(3), TTCON, TTCON2, TTSAFE
      PARAMETER (TTCON=1.4,TTCON2=TTCON*TTCON,NTAB=10,TTSAFE=2.)
      REAL*8 TTFAC, TTHA(4), TTDF(NTAB,NTAB), TTERR(NTAB), TTERRT
      REAL*8 TTP, TTTMP, TTPHI(6), TTPOS(4)


      TTPOS(1) = XI(1)
      TTPOS(2) = XI(2)
      TTPOS(3) = XI(3)
      TTPOS(4) = TG
      TTHA = (/ 0, 0, 0, 0 /)
      TTHA(TTAXIS) = 1
      TTDIR = 1

! loop on L until the error on the derivative is satisfactory      
      DO L=1,NTAB

! define the direction of update of the initial guess for the step size
! (< 0 : decrease step size)
! (> 0 : increase step size)
        IF(TTDIR .LT. 0) THEN
          TTDIR = -1
        ELSE
          TTDIR = 1
        ENDIF
      
        BESTJ = 1
        BESTK = 1
! set the initial error to a big value
        TTERR(L) = 1.0D30 

! force x and h+x to differ by an exactly representable number
! temp = x + h
! h = temp − x
        TTTMP = TTPOS(TTAXIS) + TTH(TTAXIS)
        TTHH = TTTMP - TTPOS(TTAXIS)

        CALL TTGALAXY(XI(1)+TTHA(1)*TTHH, XI(2)+TTHA(2)*TTHH, 
     &        XI(3)+TTHA(3)*TTHH, TG+TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
        TTPHI(1) = TTTMP
        CALL TTGALAXY(XI(1)-TTHA(1)*TTHH, XI(2)-TTHA(2)*TTHH, 
     &        XI(3)-TTHA(3)*TTHH, TG-TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
        TTPHI(2) = TTTMP
        CALL TTGALAXY(XI(1)+2.*TTHA(1)*TTHH, XI(2)+2.*TTHA(2)*TTHH, 
     &        XI(3)+2.*TTHA(3)*TTHH, TG+2.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
        TTPHI(3) = TTTMP
        CALL TTGALAXY(XI(1)-2.*TTHA(1)*TTHH, XI(2)-2.*TTHA(2)*TTHH, 
     &        XI(3)-2.*TTHA(3)*TTHH, TG-2.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
        TTPHI(4) = TTTMP

! form the first estimate of the 3rd order derivative
        IF(TTORDER .EQ. 3.0) THEN
          TTDF(1,1) = (-TTPHI(4)+2.*TTPHI(2)-2.*TTPHI(1)+TTPHI(3))
     &              /(2.*TTHH**3)
        ENDIF

! form the first estimate of the 5th order derivative
        IF(TTORDER .EQ. 5.0) THEN
          CALL TTGALAXY(XI(1)+3.*TTHA(1)*TTHH, XI(2)+3.*TTHA(2)*TTHH, 
     &        XI(3)+3.*TTHA(3)*TTHH, TG+3.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
          TTPHI(5) = TTTMP
          CALL TTGALAXY(XI(1)-3.*TTHA(1)*TTHH, XI(2)-3.*TTHA(2)*TTHH, 
     &        XI(3)-3.*TTHA(3)*TTHH, TG-3.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
          TTPHI(6) = TTTMP

          TTDF(1,1) = (-TTPHI(6) + 4.*TTPHI(4) - 5.*TTPHI(2)
     &      + 5.*TTPHI(1) - 4.*TTPHI(3) + TTPHI(5)) / (120.*TTHH**5)
        ENDIF


        DO 12 K=2, NTAB
! adjust step size
          TTHH = TTHH / TTCON

! re-evaluate derivative at the relevant order
          CALL TTGALAXY(XI(1)+TTHA(1)*TTHH, XI(2)+TTHA(2)*TTHH, 
     &        XI(3)+TTHA(3)*TTHH,TG+TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
          TTPHI(1) = TTTMP
          CALL TTGALAXY(XI(1)-TTHA(1)*TTHH, XI(2)-TTHA(2)*TTHH, 
     &        XI(3)-TTHA(3)*TTHH, TG-TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
          TTPHI(2) = TTTMP
          CALL TTGALAXY(XI(1)+2.*TTHA(1)*TTHH, XI(2)+2.*TTHA(2)*TTHH, 
     &        XI(3)+2.*TTHA(3)*TTHH, TG+2.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
          TTPHI(3) = TTTMP
          CALL TTGALAXY(XI(1)-2.*TTHA(1)*TTHH, XI(2)-2.*TTHA(2)*TTHH, 
     &        XI(3)-2.*TTHA(3)*TTHH, TG-2.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
          TTPHI(4) = TTTMP

          IF(TTORDER .EQ. 3.0) THEN
            TTDF(1,K) = (-TTPHI(4)+2.*TTPHI(2)
     &               -2.*TTPHI(1)+TTPHI(3))/(2.*TTHH**3)
          ENDIF

          IF(TTORDER .EQ. 5.0) THEN
            CALL TTGALAXY(XI(1)+3.*TTHA(1)*TTHH, XI(2)+3.*TTHA(2)*TTHH, 
     &        XI(3)+3.*TTHA(3)*TTHH, TG+3.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
            TTPHI(5) = TTTMP
            
            CALL TTGALAXY(XI(1)-3.*TTHA(1)*TTHH, XI(2)-3.*TTHA(2)*TTHH, 
     &        XI(3)-3.*TTHA(3)*TTHH, TG-3.*TTHA(4)*TTHH,
     &        RBAR, ZMBAR, VSTAR, TSTAR, TTTMP, TTTDEP)
            TTPHI(6) = TTTMP

            TTDF(1,K) = (-TTPHI(6) + 4.*TTPHI(4) - 5.*TTPHI(2)
     &      + 5.*TTPHI(1) - 4.*TTPHI(3) + TTPHI(5)) / (120.*TTHH**5)
          ENDIF

          TTFAC = TTCON2

! apply Neville's method
          DO 11 J=2, K 
            TTDF(J,K)=(TTDF(J-1,K)*TTFAC-TTDF(J-1,K-1))/(TTFAC-1.)
            TTFAC = TTCON2*TTFAC
          
            TTERRT = MAX(ABS(TTDF(J,K)-TTDF(J-1,K)),
     &          ABS(TTDF(J,K)-TTDF(J-1,K-1)))

            IF (TTERRT .LE. TTERR(L)) THEN
              TTERR(L) = TTERRT
              BESTJ = J
              BESTK = K
            ENDIF
   11     CONTINUE

          IF (ABS(TTDF(K,K)-TTDF(K-1,K-1)) .GE. TTSAFE*TTERR(L)) EXIT
   12   CONTINUE
    

! reset error to large value if derivative is zero
        IF(TTDF(BESTJ,BESTK) .EQ. 0.0) THEN
          TTERR = 1.0D30
        ELSE
! temp = x + h
! h = temp − x
          TTTMP = TTPOS(TTAXIS) + 
     &       ABS(2.220446D-16*TTP/TTDF(BESTJ,BESTK))**(1./TTORDER)
          TTH(TTAXIS) = TTTMP - TTPOS(TTAXIS)

! compute relative error
          TTERR(L) = ABS(TTERR(L) / TTDF(BESTJ,BESTK))
        ENDIF

        IF(TTERR(L) .LT. 1.0D-4) THEN
          EXIT
        ELSE
! if relative error is to big, update step size
! by either increasing it or decreasing it
! (this actually depends on the potential at current position)
          IF(L .GT. 1) THEN
            IF(TTERR(L) .GT. TTERR(L-1)) TTDIR = -2*TTDIR
          ENDIF
          TTH(TTAXIS) = TTH(TTAXIS) * (10.0)**TTDIR
        ENDIF

      ENDDO
! end loop over L

! if derivative is zero, return a dummy step size
      IF(TTDF(BESTJ,BESTK) .EQ. 0.0) THEN
          TTTMP = TTPOS(TTAXIS) + 
     &       ABS(1d-5 * TTPOS(TTAXIS))
          TTH(TTAXIS) = TTTMP - TTPOS(TTAXIS)
      ENDIF
      
C      write(1,*) XI(1), TTH(TTAXIS), TTERR(L), L

      RETURN
      END

