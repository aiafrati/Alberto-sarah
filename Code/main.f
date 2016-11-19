      program main

C     +++++++++++++++++++++++++++++++++++++++++++++++++++++
c     ++                                                 ++
c     +  Fully-nonlinear BEM solver for the water entry   +
c     +    of a 2D/Axisymmetric, QUASI-arbitrary shaped   +
c     +                      body                         +
c     +   It is suitable for both, pure water entry or    +
c     +                 2D+t simulations                  +
c     ++                                                 ++
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

c     definition of variables and arrays

      implicit none

c     npamx = max number of panels
c     nfil  = max order of the Dold filter

      integer npamx,nfil
      parameter (npamx=5000, nfil=7)

c     scon, sve, spre, spot: first two letters of the output files for
c     configurations, velocity, pressure and potential

      character*2 scon,svel,spre,spot

c     yv,zv; yn,zn: coordinates of the panels vertices at the two 
c                   Runge-Kutta levels
c     yce,zce; ycn,zcn: centroids of the panels at the two 
c                   Runge-Kutta levels
c     tmy,tmz: components of the unit tangent vector
c     rny,rnz: components of the unit normal vector (ORIENTED INWARDS) 
c     amp:     panel length
c     vym1,vzm1;vym2,vzm2: velocity components of the panel centroids 
c                    at the two Runge-Kutta levels
c     yg,zg;ygn,zgn: coordinate of the body geometry (2D contour in the
c                    2D+t case) at the original and actual positions
c     tg,ygs2,zgs2: variables for spline reconstruction of the body
c                   contour

      real*8 yv(npamx),zv(npamx),yn(npamx),zn(npamx)
      real*8 yce(npamx),zce(npamx),ycn(npamx),zcn(npamx)
      real*8 tmy(npamx),tmz(npamx),rny(npamx),rnz(npamx),amp(npamx)
      real*8 phi(npamx),phin(npamx),dphi(npamx),dpht(npamx)
      real*8 vym1(npamx),vym2(npamx),vzm1(npamx),vzm2(npamx)
      real*8 yg(npamx),zg(npamx),ygn(npamx),zgn(npamx)
      real*8 tg(npamx),ygs2(npamx),zgs2(npamx)


c     vfall: vertical velocity (it can be either constant of variable
c            depending on kvfall)
c     pro0: initial submergence of the apex of the body
c     ampp: initial FS panel size near the body contour
c     pfraz: ratio ampp/smallest panel size
c     ancut: angle of cut of the jet
c     escr: growth factor of panel size (bulk region)
c     estr: domain size (maximum x-coordinate on the right)
c     tend: simulation time
c     frdt: CFL limit (generally 0.25) 
c     amplim: minimum panel amplitude for discretization

      real*8 vfall,pro0,ampp,pfraz,ancut,escr,estr,tend,frdt,amplim,
     &     pi,t

c     npc: number of panel on the body contour (not including the
c          shallow water part)
c     npf: number of panel on the free surface (not including the
c          shallow water part)
c     ngo: number of points of the body geometry (2D contour in case of
c          2D+t)
c     krest: restart option (0: start from undisturbed solution, 1:
c          start from a previous solution given as a "restart" file)
c     k2dt: 2D+t option (0: vertical water entry, 1: 2D+t solution)
c     kvfall: type of entry velocity (0: constant entry velocity, 1:
c          variable entry velocity). In the 2D+t case, 0 can be used for
c          steady planing and constant inclination of the keel 
c          (e.g. wedge shaped), otherwise 1 is needed
c     ksta: number of time steps between two successive outputs
c     ift: time steps for the action of the Dold filter (typically 4)
c     iford: order of the Dold filter (typically 3)
c     llf: index of the output files

      integer npc,npf,ngo,npt
      integer krest,k2dt,kvfall,ksta,ift,iford,llf

c     kphi: index 0 for Neumann BC, 1 for Diriclet, other values for
c          specific use (e.g. jet modelling)

      integer kphi(npamx)

c     Dold filter variables

      real*8   cdold(1:nfil,0:nfil),ddold(1:nfil)

c     Local Variables

      integer i

c     First executable statement

c     Input of data for time integration and body shape

      call input (krest,k2dt,vfall,kvfall,pro0,ampp,pfraz,ancut,
     &    escr,estr,tend,frdt,ksta,scon,svel,spot,spre,ift,iford,
     &    npamx,ngo,yg,zg)

      call initial (krest,pro0,ampp,pfraz,escr,estr,npamx,ngo,
     &   nfil,yg,zg,npc,npf,npt,yv,zv,yce,zce,amp,amplim,tmy,
     &   tmz,rny,rnz,phi,dphi,kphi,cdold,ddold,t,llf,tg,ygs2,zgs2,
     &   ygn,zgn,pi)
   

      do i = 1,npc+1
        write(58,*) yv(i),zv(i)
      enddo
      do i = npc+1,npc+npf+1
        write(59,*) yv(i),zv(i)
      enddo




      stop
      end
