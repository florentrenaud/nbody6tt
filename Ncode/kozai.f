      SUBROUTINE KOZAI(IPAIR,IM,ECC1,SEMI1,ITERM)
*
*
*       Kozai relations.
*       ----------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  EI(3),HI(3),HO(3),A1(3),A2(3)
      SAVE IT
      DATA IT /0/
*
*
*       Define KS & c.m. indices for outer orbit with ECC1 & SEMI1.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      ITERM = 0
*
*       Evaluate inner semi-major axis and eccentricity.
      RI = 0.0
      TD2 = 0.0
      DO 1 K = 1,3
          RI = RI + UM(K,IM)**2
          TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
    1 CONTINUE
      ZMB = CM(1,IM) + CM(2,IM)
      SEMI = -0.5*ZMB/HM(IM)
      ECC2 = (1.0 - RI/SEMI)**2 + TD2**2/(ZMB*SEMI)
      ECC = SQRT(ECC2)
*
*       Form useful scalars (XREL & VREL consistent with UM & UMDOT).
      A12 = 0.0
      A22 = 0.0
      A1A2 = 0.0
      RI2 = 0.0
      VI2 = 0.0
      RV = 0.0
      DO 5 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A1(K) = XREL(K1,IM)*VREL(K2,IM) - XREL(K2,IM)*VREL(K1,IM)
          A2(K) = (X(K1,I2) - X(K1,I1))*(XDOT(K2,I2) - XDOT(K2,I1))
     &          - (X(K2,I2) - X(K2,I1))*(XDOT(K1,I2) - XDOT(K1,I1))
          A12 = A12 + A1(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A1(K)*A2(K)
          RI2 = RI2 + XREL(K,IM)**2
          VI2 = VI2 + VREL(K,IM)**2
          RV = RV + XREL(K,IM)*VREL(K,IM)
    5 CONTINUE
*
*       Construct the Runge-Lenz vector (Heggie & Rasio 1995, Eq.(5)).
      EI2 = 0.0
      DO 10 K = 1,3
          EI(K) = (VI2*XREL(K,IM) - RV*VREL(K,IM))/BODY(I1) -
     &                              XREL(K,IM)/SQRT(RI2)
          EI2 = EI2 + EI(K)**2
   10 CONTINUE
      EI2 = MIN(EI2,0.999999d0)
*
*       Define unit vectors for inner eccentricity and angular momenta.
      COSJ = 0.0
      SJSG = 0.0
      DO 15 K = 1,3
          EI(K) = EI(K)/SQRT(EI2)
          HI(K) = A1(K)/SQRT(A12)
          HO(K) = A2(K)/SQRT(A22)
          COSJ = COSJ + HI(K)*HO(K)
          SJSG = SJSG + EI(K)*HO(K)
   15 CONTINUE
*
*       Evaluate the expressions A & Z.
      A = COSJ*SQRT(1.0 - EI2)
      Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
      Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
      EMAX = ONE6*(Z + 1.0 - 4.0*A**2 + SQRT(Z2))
      EMAX = MAX(EMAX,0.0001d0)
      EMAX = SQRT(EMAX)
*       Skip maximum eccentricity below 0.90.
      IF (EMAX.LT.0.0) GO TO 50
*
*       Form minimum eccentricity (Douglas Heggie, Sept. 1996).
      AZ = A**2 + Z - 2.0
      IF (AZ.GE.0.0) THEN
          AZ1 = 1.0 + Z - 4.0*A**2
          EMIN2 = ONE6*(AZ1 - SQRT(AZ1**2 - 12.0*AZ))
      ELSE
          EMIN2 = 1.0 - 0.5*(A**2 + Z)
      END IF
      EMIN2 = MAX(EMIN2,0.0001D0)
      EMIN = SQRT(EMIN2)
*
*       Estimate Kozai time-scale (N-body units).
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I1))
      TK1 = TWOPI*ABS(SEMI1)*SQRT(ABS(SEMI1)/BODY(I))
      TKOZ = TK1**2*BODY(I)*(1.0 - ECC1**2)**1.5/(BODY(I2)*TK)
*
*       Evaluate numerical precession factor in Myr (elliptic integral).
      CONST = PFAC(A,Z)
*       Note published factor 2/(3*pi) (Kiseleva, Eggleton & Mikkola MN/98).
      CONST = CONST*4.0/(1.5*TWOPI*SQRT(6.0))
      TKOZ = CONST*TKOZ
*
*       Obtain inclination (see routine INCLIN).
      FAC = MIN(COSJ,1.0D0)
      ANGLE = ACOS(FAC)
*
*       Check optional output (every 100 orbits or each time for #42 > 1).
      IT = IT + 1
      IF ((MOD(IT,100).EQ.0.AND.EMAX.GT.0.99).OR.KZ(42).GT.1) THEN
          ZI = ANGLE*360.0/TWOPI
          WRITE (42,30)  TIME+TOFF, NAME(I1), ZI, EMIN, EMAX, ECC,
     &                   SEMI1/SEMI, TKOZ, TMDIS(IM)
   30     FORMAT (' KOZAI    T NAM1 INC EM EX E A1/A TKOZ TMD ',
     &                       F9.3,I7,F7.1,F6.2,2F8.4,F6.1,1P,2E9.1)
          CALL FLUSH(42)
          IF (IT.GT.2000000000) IT = 0
      END IF
*
   50 RETURN
*
      END
