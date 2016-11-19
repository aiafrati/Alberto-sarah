
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine input (krest,k2dt,vfall,kvfall,pro0,ampp,pfraz,ancut,
     &    escr,estr,tend,frdt,ksta,scon,svel,spot,spre,ift,iford,
     &    npamx,ngo,yg,zg)

      implicit none

c     npamx = max number of panels on entry

      integer npamx

c     yg,zg; coordinate of the body geometry (2D contour in the
c                    2D+t case) at the original and actual positions

      real*8 yg(npamx),zg(npamx)

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

      real*8 vfall,pro0,ampp,pfraz,ancut,escr,estr,tend,frdt

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

      integer ngo
      integer krest,k2dt,kvfall,ksta,ift,iford

c     scon, sve, spre, spot: first two letters of the output files for
c     configurations, velocity, pressure and potential

      character*2 scon,svel,spre,spot

c     local variables

      integer ig

c     First executive statement

      write(*,*) '--> read input data for discretization and time 
     & integration '

      open(8,file='slam.inp',status='UNKNOWN')
        read(8,*) krest
        read(8,*) k2dt
        read(8,*) vfall
        read(8,*) kvfall
        read(8,*) pro0 
        read(8,*) ampp 
        read(8,*) pfraz
        read(8,*) ancut
        read(8,*) escr
        read(8,*) estr
        read(8,*) tend 
        read(8,*) frdt 
        read(8,*) ksta
        read(8,'(a)') scon
        read(8,'(a)') svel
        read(8,'(a)') spot
        read(8,'(a)') spre
        read(8,*) ift 
        read(8,*) iford
      close(8)

      write(*,*) '--> read body geometry '

      open(unit=7,file='geo.in',status='unknown')
        read(7,*) ngo
        do ig=1,ngo
          read(7,*) yg(ig), zg(ig)
        enddo
      close(7)

      return
      end
