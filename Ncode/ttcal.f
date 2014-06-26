      SUBROUTINE TTCAL
*
*
*       Computation of tidal elements by quadratic interpolation
*       --------------------------------------------------------
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      REAL*8 DT21,DT32,DT31,A,B,T1,T2,T3, TNOW
      REAL*8 SAVETNOW
      INTEGER INDTT
      LOGICAL FIRST

      SAVE FIRST, INDTT, SAVETNOW
      DATA FIRST /.TRUE./
* this seems to be compiler-dependent but it is used elsewhere in NBODY6,
*   so it should work the same here

* INDTT saves the id of the tensor whose timestamp is immediatly after
*   the current time, set in TNOW
*   We assume that the update of the tidal tensor is only to be done forward
*   in time. If it is not the case (i.e. a tensor before the current one is
*   needed) the interpolation scheme will correct for this by making a sligtly
*   bigger error.

* SAVETNOW saves the value of TNOW for the last time the tensor has been updated.

      TNOW = TIME+TOFF

      IF(FIRST) THEN
       INDTT = 2
       FIRST=.FALSE.
      ELSE
* If the time has not been updated, the tensor is up-to-date: do nothing and exit
        IF(TNOW.EQ.SAVETNOW) RETURN 
      END IF

      SAVETNOW = TNOW
      
      IF(TNOW.GT.TTTIME(NBTT)) THEN
* When TNOW is after the timestamp of the last tensor:
*  Terminate after optional COMMON save.
        WRITE(6,*) 'ERROR: the tidal tensor is not defined anymore'
        IF (KZ(1).GT.0) CALL MYDUMP(1,1)
        STOP
      ELSE
        
* Zero or one iteraton (usually):
        DO WHILE(TNOW.GT.TTTIME(INDTT))
          INDTT = INDTT + 1
        END DO

* Treat first step as simple linear interpolation
        IF(INDTT.EQ.2) THEN
          DT32 = TTTIME(INDTT) - TTTIME(INDTT-1)
          DO I=1,3
            DO J=1,3
              T2 = TTENS(I,J,INDTT-1)
              T3 = TTENS(I,J,INDTT)        
              A = (T3 - T2) / DT32
* Tensor (at t = TNOW)
              TTEFF(I,J)  = T2 + A*(TNOW - TTTIME(INDTT-1))
* First derivative
              DTTEFF(I,J) = A
            END DO
          END DO
        ELSE
      
* Use quadratic interpolation for other steps
          DT21 = TTTIME(INDTT) - TTTIME(INDTT-1)
          DT32 = TTTIME(INDTT+1) - TTTIME(INDTT)
          DT31 = TTTIME(INDTT+1) - TTTIME(INDTT-1)
          DO I=1,3
            DO J=1,3
              T1= TTENS(I,J,INDTT-1)
              T2= TTENS(I,J,INDTT)
              T3= TTENS(I,J,INDTT+1)
              A = (T2-T1)/DT21
              B = ((T3-T1)*DT21-(T2-T1)*DT31) / (DT31*DT32*DT21)

              TTEFF(I,J) = T1 + A * (TNOW-TTTIME(INDTT-1)) 
     &          + B * (TNOW-TTTIME(INDTT-1))*(TNOW-TTTIME(INDTT))
              DTTEFF(I,J) = A + B *(2.0*TNOW-TTTIME(INDTT-1)
     &          -TTTIME(INDTT))
            END DO
          END DO
        END IF
      END IF

      RETURN
      END
