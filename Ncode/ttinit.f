      SUBROUTINE TTINIT
*
*       Read, scale and initialize the tidal tensors
*       --------------------------------------------------------
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
              
* Force the center of mass correction when the tensor is used
      IF (KZ(31).NE.1) THEN
        WRITE (6,*)'Warning! The center of mass will be corrected',
     &        ' after energy check. (KZ(31) forced to 1).'
        KZ(31) = 1
      ENDIF
      
* Reset
      NBTT = 0
      TTUNIT = 0.0
      TTOFFSET = 0.0
      DO I=1,NBTTMAX
        TTTIME(I) = 0.0
        DO J=1,3
          TTENS(J,1,I) = 0.0
          TTENS(J,2,I) = 0.0
          TTENS(J,3,I) = 0.0
        END DO
      END DO


* Read the tt.dat file                
      OPEN (UNIT=20,STATUS='UNKNOWN',FORM='FORMATTED',FILE='tt.dat')
      READ(20,*), NBTT, TTUNIT, TTOFFSET

      IF (NBTT.GT.NBTTMAX) THEN
        WRITE (6,*) '*** ERROR: increase NBTTMAX in param.h'
        STOP
      END IF
      
      DO I=1,NBTT
        READ(20,*), TTTIME(I),
     &        (TTENS(J,1,I),TTENS(J,2,I),TTENS(J,3,I),J=1,3)
      ENDDO
      CLOSE(20)

      WRITE (6,10) NBTT, TTOFFSET, TTUNIT
   10 FORMAT (/, 12X,'TIDAL TENSORS READ: ',I10,
     &      '   TIME OFFSET: ',F9.3,
     &     '   GALACTIC TIME: ',F9.3)
      CALL FLUSH(6)

* Scale the tidal tensor, using TTUNIT i.e. the time unit of the
** galactic run, in Myr  
      DO K=1,NBTT
* scale time
        TTTIME(K) = TTTIME(K)*(TTUNIT/TSTAR)+TTOFFSET/TSTAR
        DO I=1,3
          DO J=1,3
* scale tensor components
            TTENS(I,J,K)= TTENS(I,J,K)*(TSTAR/TTUNIT)**2
          END DO
        END DO
      END DO

* Check the tensors are defined for the entire run      
      IF(TIME+TOFF.LT.TTTIME(1)) THEN
        WRITE(6,*) '*** ERROR: the tidal tensor is not defined at t0'
        STOP
      END IF

      IF(TCRIT.GT.TTTIME(NBTT)) THEN
          WRITE (6,*) '*** WARNING: the tidal tensor is undefined ',
     &         'after t=', TTTIME(NBTT)
      END IF


* no coriolis force with the tidal tensor (is this needed?)
      TIDAL(4) = 0.0D0

* set the tidal radius arbitrarily
      RTIDE = 10.0*RSCALE
      RTIDE0 = RTIDE
      
* Initialisation of tidal tensor and derivatives
      CALL TTCAL

*** FRenaud

      RETURN
      END
