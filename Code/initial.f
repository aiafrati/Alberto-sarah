
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine initial(

c      entry variables
     & krest,pro0,ampp,pfraz,escr,estr,npamx,ngo,
     & nfil,yg,zg,
c      variables on exit
     & npc,npf,npt,yv,zv,yce,zce,amp,amplim,tmy,tmz,rny,rnz,phi,dphi,
     & kphi,cdold,ddold,t,llf,tg,ygs2,zgs2,ygn,zgn,pi)

      implicit none

      integer krest,npamx,npc,npf,npt,ngo,llf

      real*8 pro0,ampp,pfraz,escr,estr,amplim,t,pi,vfall

      real*8 yv(npamx),zv(npamx),yce(npamx),zce(npamx)
      real*8 tmy(npamx),tmz(npamx),rny(npamx),rnz(npamx),amp(npamx)
      real*8 phi(npamx),dphi(npamx),kphi(npamx)
      real*8 yg(npamx),zg(npamx),ygn(npamx),zgn(npamx)
      real*8 tg(npamx),ygs2(npamx),zgs2(npamx)

c     Dold Filter variables

      integer  nfil
      real*8   cdold(1:nfil,0:nfil),ddold(1:nfil)

c     Local variables

      integer i,kint
      logical iint
      real*8 yp1,ypn,dy,dz,dl,tai,am,amy,amz,dt,dyi,dzi,rm1,rq1,
     &       rm2,rq2,tlo,yin,zin,ampl

c     First executable statement

c     initialization of the Dold filter variables 

      call filiniz(nfil,cdold,ddold)

c     initialization of some generally used variables

      pi    = acos(-1.d0)
      amplim = ampp/pfraz  ! va inizializzata anche in restart

      if (krest.eq.0) then

c       initialization of time integration variables at the beginning

        t      = 0.d0 ! time
        llf    = 0    ! output configuration number

c       initialization of the spline representation for the body contour
c       starting from the geometry in input
c       NOTE: be sure that the dataset of points in geo.in is large
c       enough to assure an accurate reconstruction

        tg(1) = 0.d0
        do i=2,ngo
           dy = yg(i) - yg(i-1)
           dz = zg(i) - zg(i-1)
           dl = sqrt(dy**2 + dz**2)
           tg(i) = tg(i-1) + dl
        enddo

        yp1 = 1.e31
        ypn = 1.e31
        call spline(yg,zg,ygs2,zgs2,tg,yp1,ypn,ngo,npamx)

c       move the body to the actual position (zgn(1)=pro0). Note that
c       the second derivatives are independent of the vertical translation

        do i =1,ngo
          zgn(i) = zg(i)+pro0
          ygn(i) = yg(i)
        end do

c       define the wetted portion of the body by intersecting the still
c       liquid surface with the body contour. By using the notation
c       z=m*y+q, the equation of the still liquid surface is m=0, q=pro0

        rm1 = 0.d0
        rq1 = 0.d0 
        i  = 0
        iint = .false.    ! becomes 1 when intersection is found
        do while (.not.iint.and.i.lt.ngo)
          i = i+1
          if (abs(ygn(i)-ygn(i+1)).gt.1.d-8) then
            rm2 = (zgn(i)-zgn(i+1))/(ygn(i)-ygn(i+1))
            rq2 = zgn(i)-rm2*ygn(i)
            if (abs(rm1-rm2).gt.1.d-8) then
              yin = - (rq2-rq1)/(rm1-rm2)
              zin = rm1*yin+rq1
              iint = .true.
            else
              stop 'Intersection not found: Parallel lines '
            end if
          else   ! vertical line
            yin = ygn(i)
            zin = rm1*yin+rq1
            iint = .true.
          end if
        enddo

        if (.not.iint) stop 'Intersection not found within ngo '

        kint = i   ! is the last body node before the intersection
        dyi  = yin-ygn(kint)
        dzi  = zin-zgn(kint)

c       curvilinear abscissa at the intersection point (needed for
c       spline discretization)

        tai  = tg(kint) + sqrt( dyi**2 + dzi**2 )

c       discretization of the wetted portion of the body contour
c       at the initial time step a minimum of 6 panels is used

        npc = 6
        dt  = tai/float(npc)
        do i = 1,npc+1
           tlo = float(i-1)/float(npc)*tai
           call splint(tlo,yin,zin,ygn,zgn,ygs2,zgs2,tg,ngo,npamx)
           yv(i)  = yin
           zv(i)  = zin 
        enddo        

c       discretization of the free surface

        yin = yv(npc+1)
        npf  = log(1.d0+(estr-yin)/ampp*(escr-1.d0))/log(escr) + 1
        ampl = (estr-yin)*(1.d0-escr)/(1.d0-escr**npf)

        do i = npc+2,npc+npf+1
          yv(i) = yv(i-1) + ampl
          zv(i) = 0.d0
          ampl  = ampl*escr
        enddo
        npt  = npc + npf  

        if (npt.gt.npamx) stop 'Maximum panel number exceeded '

c       initialization of panel data and boundary conditions

        do i = 1,npc
          yce(i)  = (yv(i+1)+yv(i))/2.d0
          zce(i)  = (zv(i+1)+zv(i))/2.d0
          amy     =  yv(i+1) - yv(i)
          amz     =  zv(i+1) - zv(i)
          am      = sqrt(amy*amy+amz*amz)
          amp(i)  = am
          tmy(i)  = amy/am     
          tmz(i)  = amz/am     
          rny(i)  =  tmz(i)
          rnz(i)  = -tmy(i)
          dphi(i) = -vfall*rnz(i)
          kphi(i) = 0
        enddo
        do i = npc+1,npt   
          yce(i)  = (yv(i+1)+yv(i))/2.d0
          zce(i)  = (zv(i+1)+zv(i))/2.d0
          amy     =  yv(i+1) - yv(i)
          amz     =  zv(i+1) - zv(i)
          am      = sqrt(amy*amy+amz*amz)
          amp(i)  = am
          tmy(i)  = amy/am     
          tmz(i)  = amz/am     
          rny(i)  =  tmz(i)
          rnz(i)  = -tmy(i)
          phi(i)  = 0.d0
          kphi(i) = 1
        enddo

      else    ! Restart option
        stop ' Restart option not included yet '
      end if


c   ---- here below some quantities (mainly related to the jet and its
c   modelling) which were in the old version and have not been initialized
c   yet.
c   To be kept until the full funcionalities are restored

CKK* -- variabili del getto separato
CKK
CKK        kmed = 0 
CKK        ksup = 0 
CKK        ngo1=ngo
CKK        nsep=0
CKK        nsepo=0
CKK        do i=1,npamx
CKK          ksep(i)=0
CKKc          kord(i)=1
CKK          kord(i)=2
CKK        enddo
CKK
CKK* - getto
CKK        ng   = 0
CKK        kget = 0
CKK        ampli = 0.d0
CKK        epsg = epsgg*ampli 

      return
      end
