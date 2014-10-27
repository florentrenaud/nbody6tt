      SUBROUTINE TTGALAXY(X,Y,Z,T,RSCALE,MSCALE,VSCALE,TSCALE,TTPHIG)
*
*       User-defined galactic potential
*       -------------------------------

*** FlorentR - new subroutine

************ IMPORTANT NOTES:    (please read carefully)
*
* Compute the galactic potential TTPHIG as a function of the position
* X, Y, Z and time T.
*
* common6 is not included to avoid the user to name its variables at
* those of the common6. Instead, the parameters are passed as arguments.
*
* To convert a value in physical units into nbody6 units, *divide* the
* physical quantity in Msun, pc, Myr, km/s
* by: MSCALE,RSCALE,TSCALE,VSCALE for the length, mass, velocity
* and time respectively.
*
* Keep in mind that this routine is executed 25 times for every
* star and tail member, at every timestep. That's a lot! So make sure
* you optimize your code!
* Avoid doing the same computation twice. Put as much as possible
* in the "FIRST" block (see below).
*
* Friendly advice: do not convert X, Y, Z, T
* into physical units (kpc, Myr). Convert the parameters instead, in the
* "FIRST" block if possible.
*
* Second friendly advice (because I am very friendly!): in subroutines,
* do not overwrite the potential. Do this instead: ttphig = ttphig + ...
* This way, you can easily sum up multiple components!
*

      REAL*8 X,Y,Z,T, MSCALE, RSCALE, TSCALE, VSCALE, TTPHIG

      TTPHIG = 0.D0 ! init (= no galaxy)

!      CALL pointmass(ttphig, x, y, z, mscale,rscale,tscale,vscale)
!      CALL nfw(ttphig, x, y, z, mscale,rscale,tscale,vscale)
!      CALL nfw2(ttphig, x, y, z, mscale,rscale,tscale,vscale)
!      CALL isothermal(ttphig, x, y, z, mscale,rscale,tscale,vscale)
!      CALL bulgediskhalo(ttphig, x, y, z, mscale,rscale,tscale,vscale)
!      CALL spiralarms(ttphig, x, y, z, mscale,rscale,tscale,vscale)


      RETURN
      END

************************************************************************
************************************************************************
************************************************************************
      SUBROUTINE pointmass(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale

      real*8 mass
      save mass

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        mass = 1e10 / mscale ! Msun -> Nbody Units

        first = .false.
      ENDIF

* compute variables here
      ttphig = ttphig - mass / sqrt(x**2+y**2+z**2)

      RETURN
      END

************************************************************************
      SUBROUTINE nfw(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale

      real*8 m, r0, r
      save m, r0
 
      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        m = 1e13 / mscale  ! Msun -> Nbody Units
        r0 = 30e3 / rscale ! pc -> Nbody Units

        first = .false.
      ENDIF

* compute variables here
      r = sqrt(x**2+y**2+z**2)
      ttphig = ttphig - m / r * log(1D0 + r/r0)

      RETURN
      END

************************************************************************
      SUBROUTINE nfw2(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale

      real*8 mvir, rs, H0, c, r200, pi, r
      save mvir, rs, c

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        mvir = 1.23805D12 / mscale  ! Msun -> Nbody units
        H0 = 71.0 * 1D-6 / vscale   ! km/s/Mpc -> km/s/pc ->  Nbody
        c = 15.8163D0 ! dimensionless
        pi = atan(1.0)*4D0

        r200 = (mvir / ( (4D0/3D0)*pi*200D0 * 
     &     3D0*H0**2/(8D0*pi) ))**(1D0/3D0)
        rs = r200/c ! already in Nbody units

        first = .false.
      ENDIF

* compute variables here
      r = sqrt(x**2+y**2+z**2)
      ttphig = ttphig - mvir / r * log(1D0 + r/rs) 
     &  / (log(1D0 + c) - c/ (1D0 + c))

      
      RETURN
      END

************************************************************************
      SUBROUTINE isothermal(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale

      real*8 v02H
      save v02H

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        v02H = 0.5 * (220D0 / vscale)**2   ! Km/s -> Nbody units

        first = .false.
      ENDIF

* compute variables here
      ttphig = ttphig + v02H * log(x**2+y**2+z**2)
      
      RETURN
      END

************************************************************************
      SUBROUTINE bulgediskhalo(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale

      integer i
      real*8 mb1, rb1, mb2, rb2, phib
      real*8 md(3), a(3), b2, phid
      real*8 vh2h, r02, phih
      real*8 r2
      save mb1, rb1, mb2, rb2, md, a, b2, vh2h, r02

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

*     Bulge parameters
        mb1 = 3.0e9 / mscale   ! Msun -> Nbody units
        rb1 = 2700.0 / rscale  ! pc -> Nbody units
      
        mb2 = 1.6e10 / mscale  ! Msun -> Nbody units
        rb2 = 420.0 / rscale   ! pc -> Nbody units

*     Disk parameters
        md(1) = 8.9e10 / mscale   ! Msun -> Nbody units
        md(2) = -6.9e10 / mscale  ! Msun -> Nbody units
        md(3) = 2.8e10 / mscale   ! Msun -> Nbody units
        a(1) = 5d3 / rscale       ! pc -> Nbody units
        a(2) = 15.8d3 / rscale    ! pc -> Nbody units
        a(3) = 33d3 / rscale      ! pc -> Nbody units
        b2 = (300.0 / rscale)**2  ! pc -> Nbody units

*     Halo parameters
        vh2h = 0.5 * (225.0 / vscale)**2   ! km/s -> Nbody units
        r02 = (8.4d3 / rscale)**2          ! pc -> Nbody units

        first = .false.
      ENDIF

* compute variables here
      r2 = x**2 + y**2 + z**2
      
      phib = -mb1/sqrt(r2 + rb1**2) - mb2/sqrt(r2+rb2**2)
      
      phid = 0D0
      do i=1,3
        phid = phid - md(i)/sqrt(x**2+y**2+(a(i)+sqrt(z**2+b2))**2)
      enddo

      phih = vh2h * log(r2 + r02)

      ttphig = ttphig + phib + phid + phih
      
      RETURN
      END

************************************************************************
      SUBROUTINE spiralarms(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale


      integer N
      real*8 pi, alpha, H, Rs, r0, A, omega
      real*8 theta, xp, yp, r, KH, beta, D

      save N, pi, alpha, H, Rs, r0, A, omega


      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

        pi = 4.0D0*atan(1.0)
*     Spiral arm parameters for Cox & Gomez (2002), A&A
        N = 2
        alpha = 15.5*pi/180.0
        H = 300.0 / rscale   ! pc -> Nbody units
        Rs = 2600.0 / rscale ! pc -> Nbody units
        r0 = 5600.0 /rscale  ! pc -> Nbody units
        A = 0.0336 / mscale * rscale**3  ! Msun / pc^3 -> Nbody units
        omega = 20.0 * 0.001023 * tscale ! km/s/kpc -> Myr^-1 -> nbody

        first = .false.
      ENDIF

* compute variables here
      theta = -omega*t
      xp = cos(theta)*x - sin(theta)*y
      yp = sin(theta)*x + cos(theta)*y
      r = sqrt(xp**2 + yp**2 + z**2)

      KH = H*real(N)/(r*sin(alpha))
      beta= KH*(1.0 + 0.4*KH)
      D = (1.0 + KH + 0.3*(KH)**2)/(1+0.3*KH)
      
      ttphig = ttphig + ( -4.0*pi*H*A*exp(-(r-r0)/Rs) *
     &    cos(2.0*(atan2(yp,xp) - log(r/r0)/tan(alpha))) ) /
     &    (cosh(KH*z/beta)**beta*(KH*D))
      
      RETURN
      END



************************************************************************
************************************************************************
************************************************************************
**** Subroutine template
      SUBROUTINE mynewpot(ttphig,x,y,z,mscale,rscale,tscale,vscale)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, mscale, rscale, tscale, vscale


      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here


        first = .false.
      ENDIF

* compute variables here

      ttphig = ttphig - 0.0
      
      RETURN
      END

************************************************************************
************************************************************************
************************************************************************
