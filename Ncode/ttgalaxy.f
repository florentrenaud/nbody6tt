      SUBROUTINE TTGALAXY(X,Y,Z,T,RSCALE,MSCALE,VSCALE,TSCALE,TTPHIG)
*
*       User-defined galactic potential
*       -------------------------------

*** FlorentR - new subroutine

* Compute the galactic potential TTPHIG as a function of the position
* X, Y, Z and time T.
*
* common6 is not included to avoid the user to name its variables at those
* of the common6. Instead, the parameters are passed as arguments.
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
!      SAVE MASS


*** Initialization (done only once)
*** Put as many thing here as possible to speed up the computation
      SAVE FIRST
      DATA FIRST /.TRUE./
      IF(FIRST) THEN

* Compute your constants
!        MASS = (-1e10 / MSCALE)


         FIRST = .FALSE.
      ENDIF
*** End of initialization




*** The rest is done as often as required:

* Compute your variables
!      Rcirc2 = X**2+Y**2
!      Rcirc = SQRT(Rcirc2)
!      R2 = Rcirc2+Z**2
!      R = SQRT(R2)
!      R = SQRT(X**2+Y**2+Z**2)



!     Simple NFW Vmax = 210 km/s, rmax = 30 kpc
      CALL nfw_simple(x*rscale, y*rscale, z*rscale, t, ttphig)

!     NFW Denis
!      CALL nfw(x*rscale, y*rscale, z*rscale,t,ttphig)

!     Singular isotherm halo
!      r = sqrt(x**2 + y**2 + z**2)
!      ttphig = 220.**2*log(r*rscale)/vscale**2

!     Combination




* Compute the value of the galactic potential TTPHIG
*      TTPHIG = 0.0D0 ! no galaxy
*      TTPHIG = MASS / R ! point-mass      

      CALL GET_POT(x*rscale, y*rscale, z*rscale, t, TTPHIG)


      TTPHIG = TTPHIG / vscale**2

*** At this point, TTPHIG contains the potential of the galaxy at the position
* x, y, z and time t, in the units of Nbody6

      RETURN
      END

****************************
*** Subroutines
****************************

      subroutine get_pot(x,y,z,t,pot)
      implicit none
      real*8 x,y,z,t,pot
      real*8 spiralpot, barpot, nfwpot

      call bulge_disk_halo(x,y,z,t,pot)
      call nfw_simple(x,y,z,t,nfwpot)
      
!     call bar(x,y,z,t,barpot)
!     call spiralarm(x,y,z,t,spiralpot)
      
      pot = pot + nfwpot

      return
      end


      subroutine spiralarm(x,y,z,t,pot)
      implicit none
      integer N
      real*8 alpha, omega, A, r0, H, Rs ! parameters
      real*8 K, beta, D, KH, G, pi, gamma
      real*8 x,y,z,t,pot,r,phi, theta, xp, yp

Cf2py intent(out) pot
      G = 0.0043
      pi = 4.0*atan(1.0)

c     Spiral arm parameters for Cox & Gomez (2002), A&A
      N = 2
      alpha = 15.5*pi/180.0 ! rad
      H = 300.0             ! pc
      Rs = 2600.0           ! pc
      r0 = 5600.0           ! pc
      A = 0.0336            ! Msun/pc^3
      omega = 20.0          ! km/s/kpc

c     Rotate coordinate with respect to bar
      omega = omega*0.001023 ! Myr^-1
      theta = -omega*t

      xp = cos(theta)*x - sin(theta)*y
      yp = sin(theta)*x + cos(theta)*y

      r = sqrt(xp**2 + yp**2 + z**2)

c     Potential
      K = real(N)/(r*sin(alpha))
      KH = K*H
      beta= KH*(1.0 + 0.4*KH)
      D = (1.0 + KH + 0.3*(KH)**2)/(1+0.3*KH)

      phi = atan2(yp,xp)
      gamma = 2.0*(phi - log(r/r0)/tan(alpha))
      pot = -4.0*pi*G*H*A*exp(-(r-r0)/Rs)*cos(gamma)
      pot = pot/(cosh(K*z/beta)**beta*(K*D))

      return
      end


*****
      subroutine bar(x,y,z,t,barpot)
      implicit none
      logical first
      real*8 x,y,z,t,barpot
      real*8 a,b,c,mbar,omega
      real*8 xp, yp, theta, pi
      save first
      data first /.true./

      pi = 4.0*atan(1.0)

Cf2py intent(out) barpot
      if (first) then
         a = 3140.0             ! pc
         b = a*0.375
         c = a*0.256
         mbar = 0.98e10         ! Msun
         
!         call inpo3(a,b,c,0.0043*mbar)
         first = .false.
      endif

      omega  = 60.0             ! km/s/kpc

c     Rotate coordinate with respect to bar
      omega = omega*0.001023 ! Myr^-1
      theta = -omega*t

      xp = cos(theta)*x - sin(theta)*y
      yp = sin(theta)*x + cos(theta)*y

! 1    call PBAR3(xp,yp,z,barpot)

      return
      end


*****

      subroutine bulge_disk_halo(x,y,z,t,pot)
      implicit none
      integer i

      real*8 x,y,z,t, pot
      real*8 G, r2
      real*8 mc1, rc1, mc2, rc2, potc   ! bulge
      real*8 md(3),a(3), b2, potd       ! disk
      real*8 vh2, r02, poth             ! halo
      real*8 f
      
      G = 0.0043   ! pc Msun^-1 (km/s)^2

c     Bulge parameters
      mc1 = 3.0e9               ! Msun
      rc1 = 2700.0              ! pc
      
      mc2 = 1.6e10              ! Msun
      rc2 = 420.0               ! pc
      
c     Disk parameters
      md(1) = 8.9e10            ! Msun
      md(2) = -6.9e10 
      md(3) = 2.8e10 
      
      a(1) = 5000.0             ! pc
      a(2) = 15800.0
      a(3) = 33000.0
      b2 = 300.0**2
      
c     Halo parameters
      vh2 = 225.0**2            ! km/s
      r02 = 8400.0**2           ! pc

c     Compute potentials
      r2 = x**2 + y**2 + z**2
      
      potc = -G*mc1/sqrt(r2 + rc1**2) - G*mc2/sqrt(r2+rc2**2)
      
      potd = 0.0
      do i=1,3
         potd = potd-G*md(i)/sqrt(x**2+y**2+(a(i)+sqrt(z**2+b2))**2)
      enddo

      poth = 0.5*vh2*log(r2 + r02)

      pot = potd
      return
      end

*****

      subroutine nfw_simple(x,y,z,t,pot)
!     Note that input units are assumed to be physical
!     Simple NFW based on rmax = 30 kpc and vmax = 210 km/s
!     Values for rs and GMs need to be solved numerically outside of this code

      implicit none

      real*8 x, y, z, t, pot, r
      real*8 GMs, r0
      
      r0  =   13.8723089938D0          ! kpc
      GMs = 2829425.85977              ! kpc (km/s)**2

      r = sqrt(x**2 + y**2 + z**2)/1D3 ! kpc

      pot = -GMs/r * log(1+r/r0)

      return
      end


*****

      subroutine nfw(x,y,z,t,pot)
!     Note that input units are assumed to be physical

      implicit none
      real*8 x,y,z,t,pot,r
      real*8 Mvir, H0, c
      real*8 G, pi
      real*8 rho_crit, r200, rs, xx
      
      r = sqrt(x**2 + y**2 + z**2)

      Mvir = 1.23805D12
      H0 = 71D-6
      c = 15.8163D0
      
      G = 0.0043D0
      pi = atan(1.0)*4D0

      rho_crit = 3D0*H0**2/(8D0*pi*G)

      r200 = (Mvir/((4D0/3D0)*pi*200D0*rho_crit))**(1D0/3D0)
      
      rs = r200/c
      xx = r/rs

      pot =  - G*Mvir/r * log(1+xx)/(log(1+c) - c/(1+c))

      return
      end


*****

      subroutine nfwt(x,y,z,t,pot)
!     Compute time dependent NFW potential according to Buist & Helmi
!     Note that input units are assumed to be physical

      implicit none
      real*8 x,y,z,t,pot, cluster_age
      real*8 H0, OM, t0, zz, tBB, dt
      real*8 ag, gamma, Ms0, Rs0, Ms, Rs, G, xx

!     Parameters
      ag = 0.1
      gamma = 2.0
      Ms0 = 1.5e11                   ! Msun
      Rs0 = 16.0e3                   ! pc
      cluster_age = 12.              ! Gyr
      tBB = 13.575                   ! Gyr

      G = 0.0043                     ! pc/Msun (km/s)^2

!     Compute redshift
      dt = 1.0e3*(tBB - cluster_age) ! Time offset
      H0 = 73.0*1.023*1.0e-6         ! 1/Myr
      OM = 0.25                      ! Omega_matter
      t0 = 2.0/(3.0*H0*sqrt(1.0-OM)) ! Typical time-scale [Myr]
      zz = ((1.0-OM)/OM)**(1.0/3.0)*(sinh((t+dt)/t0))**(-2.0/3.0) -1.0
      zz = 0.0

      Ms = Ms0*exp(-2.0*ag*zz)
      Rs = Rs0*exp(-2.0*ag*zz/gamma)

      xx = sqrt(x**2 + y**2 + z**2)/Rs
      pot = -G*Ms/Rs * log(1.0+xx)/xx/(log(2.)-0.5)
   
      return
      end
