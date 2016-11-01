      program main

C     +++++++++++++++++++++++++++++++++++++++++++++++++++++
c     ++                                                 ++
c     +  Fully-nonlinear BEM solver for the water entry   +
c     +    of a 2D/Axisymmetric, QUASI-arbitrary shaped   +
c     +                      body                         +
c     +       It also has some 2D+t capabilities          +
c     ++                                                 ++
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

c     - Include files:

c       -- "slam_p.h" contains the definition of reals and the parameters 
c                     used throughout the code
c       -- "slam_f.h" contains the definition of the main variables and
c                     the parameters used throughout the code

      include"slam_p.h"
      include"slam_f.h"

c     yv,zv; yn,zn: coordinates of the panels vertices along the solid
c                   body contour  at the two Runge-Kutta levels

c     ysl,zsl,ynsl,znsl: coordinates of vertices of the free-surface
c                   panels at the two Runge-Kutta levels

c     yce,zce; ycn,zcn: centroids of the panels lying on the body
c                   contour at the two Runge-Kutta levels

c     ycesl,zcesl,ycnsl,zcnsl: centroids of the free-surface
c                   panels at the two Runge-Kutta levels

c     phi,dphi: velocity potential and normal derivative on the body
c               panels

c     phisl,dphisl: velocity potential and normal derivative on the
c               free-surface panels

c     phinsl: free-surface velocity potential at the second R-K level

c     phin: velocity potential in the jet region at the second R-K level

c     ygb,zgb: coordinates of shallow water vertices lying on the body

c     dpht,dphtsl: tangential derivative of the velocity potential on
c                  body and free-surface, respectively

c     tmy,tmz: unit tangent vector on the body surface
c     tmysl,tmzsl: unit tangent vector on free-surface panels

c     rny,rnz: unit normal vector on the body surface
c     rnysl,rnzsl: unit normal vector on free-surface panels

c     amp,ampsl: panel length on body (amp) and free surface (ampsl)
 
c     vym1,vzm1,vym2,vzm2: velocity components on the body at the two
c                          RK levels
c     vym1sl,vzm1sl,vym2sl,vzm2sl: velocity components at free surface
c                          at the two RK levels 

c     depn1,depn2: time derivative of the velocity potential on the free
c                  surface at the two RK levels

c     depn1s,depn2s: time derivative of the velocity potential in the
c                  jet region at the two RK levels

c     dpt,dpnt: tangential derivative of the tangential (dpt) and normal 
c               (dpnt) velocity on the body and jet region

c     dptsl,dpntsl: tangential derivative of the tangential (dpt) and normal 
c               (dpnt) velocity on the free surface

c     yco,zco: centroids coordinates at the previous step used in prefin
c            to compute \partial phi/\partial t by finite differences
c     vyo,vzo,phio: velocity and potential at the previous step used in 
c            prefin to compute \partial phi/\partial t by finite differences

c     rl: matching factor in the matching region
c     a1,b1,c1,d1,e1: coefficients of the FEM model in the modelled
c                region of the jet

c     a2,b2,c2,d2,e2: as before but for the time derivative of the
c                velocity potential (pressure solution)

c     xig, zeb,zef: coordinates (in local frame of reference) of the 
c               vertices of the modelled part of the jet. The xig (or
c               \xi) coordinate is the same for the two vertices located 
c               along the edge orthogonal to the free surface, whereaz 
c               zeb is the \zeta coordinate of the vertex on the body 
c               countour (and it is zero) and zef is the coordinate at the 
c               free surface

c     xigs, zegs: as above but they are the coordinate of the midpoint
c               of the jet element

c     xigs,zegs,xigb,zegb,xigf,zegf: coordinates of the vertices of the
c               jet element in the Glabal frame of reference

c     ty,tz,ry,rz: component of the tangent and normal unit vectors in
c               the modelled jet region

c     ygn,zgn: body geometry (yg,zg) translated vertically

c     tn: is the curvilinear abscissa associated with ygn and zgn

c     ygs2,zgs2: second derivatives needed for the spline reconstruction
c              of the body contour

c     tg: curvilinear abscissa along the body (defined on the basis of
c          yg,zg) which is used in the spline interpolation

c     tgb: curvilinear abscissa used on the jet surfaces (both on the
c          body contour and on the free surface side)

c     kse, ksep: indices related which identifies the elements in the
c          jet region which are detached (separated) from the body contour

c     tcb,ycb,zcb: curvilinear abscissa and coordinates of the centroids
c          of the elements in the modelled part of the jet lying along the 
c          body contour

c     tc: same as tcb but at a different RK level

c ---------------------- to be better understood

c     dpt2, dpht2: terms used to compute the pressure in different ways
c     pre,pre2,pres: pressure computed using different methods but
c              the differences have to be better understood

c     dphtbsl,dphn,dphnb,phb,ph,phib: these arrays are used in an overlapping
c            region between the bulk of the fluid and the modelled part of 
c            the jet in order to smooth the transition between the two regions
c            They have to be further understood but it seems they are
c            only used inside solv22 and solv22p

c     hp,xj,ze,xis,zes: these variables seems still related to the jet
c               modelling but not all of them are used at present. Better 
c               keep them in the debugging phase

c     dpttsl,dpb,dptb,dpntbsl: these variables are used in solv22p and
c            seems to be related to the solution of the Laplace equation 
c            for the time derivative of the velocity potential but their 
c            use isn't clear

c     vxi,dpntt: vxi is the tangential velocity in the jet and dpntt is
c            the second derivative of the velocity potential. They seem
c            to be used for the pressure

c     kord, kor: are indices used in the modelled part of the jet and/or
c     in the transition region, but it is not clear how they work

c ----------------------

      character*2 scon,svel,spre,spot
       
      dimension yv(npamx),zv(npamx),yn(npamx),zn(npamx)
      dimension ysl(npamx),zsl(npamx),ynsl(npamx),znsl(npamx)
      dimension yce(npamx),zce(npamx),ycn(npamx),zcn(npamx)
      dimension ycsl(npamx),zcsl(npamx),ycnsl(npamx),zcnsl(npamx)
      dimension phi(npamx),phisl(npamx),phinsl(npamx),phin(npamx)
      dimension dphi(npamx),dphisl(npamx)
      dimension ygb(npamx),zgb(npamx)
      dimension dpht(npamx),dphtsl(npamx)
      dimension tmy(npamx),tmz(npamx),rny(npamx),rnz(npamx)
      dimension tmysl(npamx),tmzsl(npamx),rnysl(npamx),rnzsl(npamx)
      dimension amp(npamx),ampsl(npamx)
      dimension vym1(npamx),vzm1(npamx),vym2(npamx),vzm2(npamx)
      dimension vymsl1(npamx),vzmsl1(npamx),vymsl2(npamx),vzmsl2(npamx)
      dimension depn1(npamx),depn2(npamx)
      dimension depn1s(npamx),depn2s(npamx)
      dimension dpt(npamx),dpnt(npamx),dptsl(npamx),dpntsl(npamx)

      dimension yco(npamx),zco(npamx)
      dimension vyo(npamx),vzo(npamx),phio(npamx) 
      dimension rl(npamx) 
      dimension a1(npamx),b1(npamx),c1(npamx),d1(npamx),e1(npamx)
      dimension xig(-npamx:npamx),xigs(-npamx:npamx)
      dimension zeb(-npamx:npamx),zegs(-npamx:npamx),zef(-npamx:npamx)
      dimension xigb(-npamx:npamx),zegb(-npamx:npamx)
      dimension xigf(-npamx:npamx),zegf(-npamx:npamx)
      dimension a2(npamx),b2(npamx),c2(npamx),d2(npamx),e2(npamx)
      dimension ry(npamx),rz(npamx),ty(npamx),tz(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx),tn(npamx)
      dimension tg(npamx),tgb(npamx),ksep(npamx),kse(npamx) 
      dimension tcb(npamx),ycb(npamx),zcb(npamx),tc(npamx)

c----------------- to be better understood
      dimension dpt2(npamx),dpht2(npamx)
      dimension pre(npamx),pre2(npamx),pres(npamx)
      dimension dphtbsl(npamx),dphn(npamx),dphnb(npamx)
      dimension phb(npamx),ph(npamx),phib(npamx)
      dimension hp(npamx),xj(npamx),ze(npamx),xis(npamx),zes(npamx)
      dimension dpttsl(npamx),dpb(npamx),dptb(npamx),dpntbsl(npamx)
      dimension vxi(npamx),dpntt(npamx)
      dimension kord(npamx),kor(npamx)

c-----------------

*
*
      call input(vfall0,ancut,iiget,jjget,frint,rmg,epsgg,eskg,
     #       gfrac,pro0 ,frdt,ampp,pfraz,escr,estr,tend,ksta,scon,svel,
     #       spot,spre,idis,ift,iford,frtend,ramii,ramiii,eskkk)
      call initial(t,amplim,vfall,llf,spo0,pro0,npc,npt,npsl,
     #           kget,ng,yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,amp,ampsl,
     #           rny,rnz,rnysl,rnzsl,tmy,tmz,tmysl,tmzsl,ampp,vfall0,
     #           dphi,phisl,pfraz,estr,escr,epsg,epsgg,ampli,
     #           ygn,zgn,ygs2,zgs2,tg,ngo,ngo1,nsep,nsepo,ksep,kord,
     #         nnold,nn1old,frin,frfi,jind,tin,jfid,tfi,jend,kmed,ksup)
*
* - soluzione, calcolo velocita
      call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,
     #                  phi,dphi,phisl,dphisl,jt)
*
      call calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,
     #            phi,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phisl,
     #            rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,
     #            amp,vym1,vzm1,rny,rnz,tmy,tmz,
     #            ampsl,vymsl1,vzmsl1,rnysl,rnzsl,tmysl,tmzsl,vxi,
     #            ry,rz,ty,tz,kse)
*
* --------- CALCOLO LUNGHEZZA STRISCIA DI CONTROLLO
          rrl=0.d0
          do i=jind+1,jfid-1
            dl= sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
            rrl= rrl + dl
          enddo
          dlin= (1.d0-tin)*sqrt( (ycsl(jind+1)-ycsl(jind))**2 +
     #                    (zcsl(jind+1)-zcsl(jind))**2 )
          dlfi= tfi*sqrt( (ycsl(jfid+1)-ycsl(jfid))**2 +
     #                    (zcsl(jfid+1)-zcsl(jfid))**2 )
          rrl= rrl + dlin + dlfi
          yin = ycsl(jind)+tin*(ycsl(jind+1)-ycsl(jind))
          zin = zcsl(jind)+tin*(zcsl(jind+1)-zcsl(jind))
          yfi = ycsl(jfid)+tfi*(ycsl(jfid+1)-ycsl(jfid))
          zfi = zcsl(jfid)+tfi*(zcsl(jfid+1)-zcsl(jfid))
          write(77,'(i4,2d15.6,2i4,4d15.6)') 0 ,0.d0,rrl,jind,jfid+1,
     #               yin,zin,yfi,zfi                             
           
*-----------
*
      call stampa(vfall,ng,kget,npc,npsl,jt,t,dt,frdt,llf,
     #                  scon,svel,spot,spre,phi,dphi,phisl,dphisl, 
     #                 yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,
     #                 vym1,vzm1,vymsl1,vzmsl1,
     #                 dpht,dphtsl,dpt,dpnt,pre,pre2,pres,vxi,
     #                 jend,jind,jfid,yin,zin,yfi,zfi)
*
*----------------------------------------------------------------------
*- inizio integrazione temporale
      proat = zv(1)
      jt    = 0
      do while (t.lt.tend)
        jt=jt+1
        vfall = vfall0
        if(jt.lt.30) vfall=vfall0*float(jt)/float(30) !E LA PRIMA ITERAZIONE ?
*- calcolo passo di integrazione, e volume
        write(*,*) 'prima check, nsep',nsep
        call check(yv,zv,npc,npsl,dt,t,jt,ampsl,frdt,vymsl1,vzmsl1,
     #             ycsl,zcsl,kget,frtend,tend,zgn,ngo,vfall,ngo1,nsep)
        write(*,*) 'dt jt ',dt,jt,vfall
*- spostamento (predictor) dei centroidi 
        write(*,*) 'predictor',npsl
ccc        write(15,*) '#jt ',jt,ng
ccc        write(16,*) '#jt ',jt,ng
c      if(jt.eq.254)then
c        ksep(npc+ng)   =1
c        ksep(npc+ng-1) =1
c        ksep(npc+ng-2) =1
c        ksep(npc+ng-3) =1
c        ksep(npc+ng-4) =1
c      endif
        do ip = 1,npc+ng
ccc          write(15,'(2d15.7,i4)') vym1(ip),vzm1(ip),ksep(ip)
          if(nsep.eq.1)then
          if(ksep(ip).ne.0)then
          deph2      = vym1(ip)**2+ vzm1(ip)**2
          depn1s(ip) = deph2/2.d0
          ycn(ip)    = yce(ip) + vym1(ip)*dt/2.d0
          zcn(ip)    = zce(ip) + vzm1(ip)*dt/2.d0
          phin(ip)   = phi(ip) + depn1s(ip)*dt/2.d0
          else 
          phin(ip)   = phi(ip) 
          endif 
ccc          write(16,'(2d15.7,i4)') phi(ip),phin(ip),ksep(ip)
          if(ksep(ip).eq.2)then
            ay  = ycn(ip)  - ygn(ngo)
            ay1 = ycn(ip-1)- ygn(ngo)
            if(ay1.gt.0.d0)then
              ksep(ip) = 1
            else
              ycn(ip)  = yce(ip)
              zcn(ip)  = zce(ip)
              phin(ip) = phi(ip) 
            endif 
          endif 
          else
              ycn(ip)  = yce(ip)
              zcn(ip)  = zce(ip)
              phin(ip) = phi(ip)
          endif 
        enddo
ccc        write(15,*) 
ccc        write(16,*) 
        do ip = 1,npsl
ccc          write(15,'(2d15.7,i4)') vymsl1(ip),vzmsl1(ip),ksep(ip)
          deph2        = vymsl1(ip)**2 + vzmsl1(ip)**2
          depn1(ip)    = deph2/2.d0
          ycnsl(ip)    = ycsl(ip) + vymsl1(ip)*dt
          zcnsl(ip)    = zcsl(ip) + vzmsl1(ip)*dt
          phinsl(ip)   = phisl(ip)+ depn1(ip)*dt
ccc          write(16,'(2d15.7)') phisl(ip),phinsl(ip)
        enddo
ccc        write(15,*) 
ccc        write(15,*) 
ccc        write(16,*) 
ccc        write(16,*) 
*- sposto vertice corpo
        proat = zv(1) - vfall*dt
*
ccckkk        write(30,*) '# jt ',jt,ng
ccckkk        do i=1,npc+ng
ccckkk          write(30,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(30,*)
ccckkk        write(30,*)
ccckkk        write(31,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(31,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(31,*)
ccckkk        write(31,*)
ccckkk        write(32,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(32,'(4d15.7)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(32,*)
ccckkk        write(32,*)
ccckkk        write(33,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(33,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(33,*)
ccckkk        write(33,*)
ccckkk        write(34,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(34,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(34,*)
ccckkk        write(34,*)
*- ricostruisco la configurazione di tentativo dei vertici
        call splver2(ycnsl,zcnsl,ynsl,znsl,npsl,proat,
     #               kget,estr,ng,ygn,zgn,ygs2,zgs2,tg,ngo,iint,jt,
     #               ampli,ngo1,nsep,nsepo,ycn,zcn,npc,
     #               phin,ksep,di,ang,tc)
*
ccckkk        write(40,*) '# jt ',jt,ng
ccckkk        do i=1,npc+ng
ccckkk          write(40,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(40,*)
ccckkk        write(40,*)
ccckkk        write(41,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(41,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(41,*)
ccckkk        write(41,*)
ccckkk        write(42,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(42,'(4d15.7)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(42,*)
ccckkk        write(42,*)
ccckkk        write(43,*) '# jt ',jt,ng
ccckkk        do i=1,npsl+1   
ccckkk          write(43,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(43,*)
ccckkk        write(43,*)
ccckkk        write(44,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(44,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(44,*)
ccckkk        write(44,*)
*- rigriglio corpo
        call ridis5(0,ng,proat,kget,ynsl,znsl,ygb,zgb,
     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
     #                  nsep,nsepo,ksep,phin,
     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)
*- normali tangenti  ampiezze e boundary conditions
        call nortan(vfall,yn,zn,amp,tmy,tmz,rny,rnz,npc,
     #            ynsl,znsl,ampsl,tmysl,tmzsl,rnysl,rnzsl,npsl,
     #            dphi,kget,ng,ycn,zcn,yce,zce,yv,zv,
     #            ycnsl,zcnsl,ycsl,zcsl,ysl,zsl,phinsl,phisl,
     #            phin,phi,ygb,zgb,0,ksep,jt)
*- trattamneto getto
        if(kget.eq.1)then
            write(*,*) ' AAA Chiamo get 1 '
            call  get(ng,npc,ycnsl,zcnsl,ycn,zcn,yn,zn,  
     #               xigs,zegs,xigb,zegb,xigf,zegf) 
*
        endif
*
ccckkk        write(50,*) '# jt ',jt
ccckkk        do i=1,npc+ng
ccckkk          write(50,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(50,*)
ccckkk        write(50,*)
ccckkk        write(51,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(51,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(51,*)
ccckkk        write(51,*)
ccckkk        write(52,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(52,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(52,*)
ccckkk        write(52,*)
ccckkk        write(53,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(53,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(53,*)
ccckkk        write(53,*)
ccckkk        write(54,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(54,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(54,*)
ccckkk        write(54,*)
        
*- risolvo il problema nella configurazione intermedia
        print*,'t_old ',t
        t   = t + dt
        nng = kget*ng
        print*,'t_new ',t
        if(kget.eq.0)then
          call solv1(npc,nng,npsl,yn,zn,ygb,zgb,ynsl,znsl,amp,ampsl,
     #                  phin,dphi,phinsl,dphisl,jt)
        elseif(kget.eq.1)then
c          call solv2(frint,ng,npc,npsl,ycn,zcn,yn,zn,ynsl,znsl,
c     #            ycnsl,zcnsl,dphi,phinsl,dphtsl,dphtbsl,phb,
c     #           xigs,zegs,xigb,zegb,xigf,zegf,
c     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
c     #          ph,dphn,dphnb,a1,b1,c1,d1,e1,phi,phib,dphisl,dphibsl,rl,
c     #          mb,mf,mt,m,n,nt,ntt,xj,ze,xis,zes,jt)
            write(*,*) 'AAA chiamo solv22 - 1'
          call solv22(frint,ng,npc,npsl,ycn,zcn,yn,zn,ynsl,znsl,
     #            ycnsl,zcnsl,dphi,phinsl,dphtsl,dphtbsl,phb,
     #           xigs,zegs,xigb,zegb,xigf,zegf,
     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
     #        ph,dphn,dphnb,a1,b1,c1,d1,e1,phin,rl,
     #          mb,mf,mt,m,n,nt,ntt,xj,ze,xis,zes,jt,ksep,kse,kord,kor)
          call calsol(mb,mf,mt,m,n,nt,ntt,phin,phib,phb,dphisl,
     #           dphnb,rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,ry,rz,kse,dphi)
        endif
*- calcolo velocita
          call calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,
     #            phin,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phinsl,
     #              rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,
     #              amp,vym2,vzm2,rny,rnz,tmy,tmz,
     #              ampsl,vymsl2,vzmsl2,rnysl,rnzsl,tmysl,tmzsl,vxi,
     #              ry,rz,ty,tz,kse)
*- spostamento definitivo dei centroidi
        write(*,*) 'corrector '
        tysl = 1.d31
        do ip = 1,npc+ng
          if(ksep(ip).eq.1)then
          deph2      = vym2(ip)**2 + vzm2(ip)**2
          depn2s(ip) = deph2/2.d0
          ycn(ip)    = yce(ip)+(vym1(ip)+vym2(ip))*dt/2.d0
          zcn(ip)    = zce(ip)+(vzm1(ip)+vzm2(ip))*dt/2.d0
          phin(ip)   = phi(ip)+(depn1s(ip)+depn2s(ip))*dt/2.d0
          tysl       = min(tysl,ycn(ip))
          endif 
        enddo
        do ip = 1,npsl
          deph2       = vymsl2(ip)**2 + vzmsl2(ip)**2
          depn2(ip)   = deph2/2.d0
          ycnsl(ip)   = ycsl(ip)+(vymsl1(ip)+vymsl2(ip))*dt/2.d0
          zcnsl(ip)   = zcsl(ip)+(vzmsl1(ip)+vzmsl2(ip))*dt/2.d0
          phinsl(ip)  = phisl(ip)+(depn1(ip)+depn2(ip))*dt/2.d0
        enddo
*- ricostruisco la configurazione della superficie libera
        call splver2(ycnsl,zcnsl,ynsl,znsl,npsl,proat,
     #               kget,estr,ng,ygn,zgn,ygs2,zgs2,tg,ngo,iint,jt,
     #               ampli,ngo1,nsep,nsepo,ycn,zcn,npc,
     #               phin,ksep,di,ang,tc)
*- calcolo angolo intersezion, event. separo (e/o rigriglio) il getto
ccckkk        write(60,*) '# jt ',jt,ng,ngo1
ccckkk        do i=1,npc+ng
ccckkk          write(60,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(60,*)
ccckkk        write(60,*)
ccckkk        write(61,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(61,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(61,*)
ccckkk        write(61,*)
ccckkk        write(62,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(62,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(62,*)
ccckkk        write(62,*)
ccckkk        write(63,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(63,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(63,*)
ccckkk        write(63,*)
ccckkk        write(64,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(64,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(64,*)
ccckkk        write(64,*)

        write(*,*) 'jind jifid 0',jind,jfid
        call shallo(ynsl,znsl,ygb,zgb,kget,jt,npsl,
     #               ng,ampli,rmg,iiget,jjget,
     #               ycnsl,zcnsl,phinsl,dphisl,epsg,epsgg,eskg,gfrac,
     #               ygn,zgn,ygs2,zgs2,tg,tgb,ngo,ngo1,nsep,
     #               ycb,zcb,tcb,nngo,jind,jfid,jend)
        write(*,*) 'jind jifid 1',jind,jfid
ccckkk        write(70,*) '# jt ',jt,ng
ccckkk        do i=1,npc+ng
ccckkk          write(70,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(70,*)
ccckkk        write(70,*)
ccckkk        write(71,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(71,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(71,*)
ccckkk        write(71,*)
ccckkk        write(72,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(72,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(72,*)
ccckkk        write(72,*)
ccckkk        write(73,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(73,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(73,*)
ccckkk        write(73,*)
ccckkk        write(74,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(74,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(74,*)
ccckkk        write(74,*)
*- redistribuzione centroidi e vertici SL
        if (mmm.eq.0.or.jt.le.jjget) then
        call disun2(jt,yn,zn,ng,ynsl,znsl,ycnsl,zcnsl,phinsl,npsl,
     #    ygb,zgb,escr,kget,estr,amplim,ampli,npt,npc,iiget,
     #    ycb,zcb,jind,jfid,tin,tfi,jend)
*- ridefinizione (e ridiscretizzazione) dei pannelli sul corpo
        write(*,*) 'jind jifid 2',jind,jfid
        nng=ng*kget
ccckkk        write(80,*) '# jt ',jt,ng
ccckkk        do i=1,npc+nngo
ccckkk          write(80,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(80,*)
ccckkk        write(80,*)
ccckkk        write(81,*) '# jt ',jt
ccckkk        do i=1,npc+nngo+1
ccckkk          write(81,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(81,*)
ccckkk        write(81,*)
ccckkk        write(82,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(82,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(82,*)
ccckkk        write(82,*)
ccckkk        write(83,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(83,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(83,*)
ccckkk        write(83,*)
ccckkk        write(84,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(84,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(84,*)
ccckkk        write(84,*)
c        if(jt.eq.iiget+1) stop
        nng=ng*kget
        mmm=mod(jt,idis)
        call ridis5(1,ng,proat,kget,ynsl,znsl,ygb,zgb,
     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
     #                  nsep,nsepo,ksep,phin,
     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)
        else
*- ridefinizione (e ridiscretizzazione) dei pannelli sul corpo
        call ridis5(0,ng,proat,kget,ynsl,znsl,ygb,zgb,
     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
     #                  nsep,nsepo,ksep,phin,
     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)
        end if
        nng=ng*kget
*- filtro la superficie libera
        mm = mod(jt,ift)
        if (mm.eq.0) then
          if(nsep.eq.1)then
            call doldfil1(yn,npc+2,npc+nng-2,npc+nng+1)
            call doldfil1(zn,npc+2,npc+nng-2,npc+nng+1)
          endif
          call doldfil1(ynsl,1+iford,npsl-10,npsl+1)
          call doldfil1(znsl,1+iford,npsl-10,npsl+1)
          call doldfil1(phinsl,1+iford,npsl-10,npsl)
        end if
          call nortan(vfall,yn,zn,amp,tmy,tmz,rny,rnz,npc,
     #            ynsl,znsl,ampsl,tmysl,tmzsl,rnysl,rnzsl,npsl,
     #            dphi,kget,ng,ycn,zcn,yce,zce,yv,zv,
     #            ycnsl,zcnsl,ycsl,zcsl,ysl,zsl,phinsl,phisl,
     #            phin,phi,ygb,zgb,1,ksep,jt)
ccckkk        write(90,*) '# jt ',jt,ng
ccckkk        do i=1,npc+ng
ccckkk          write(90,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(90,*)
ccckkk        write(90,*)
ccckkk        write(91,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(91,*) yn(i),zn(i),ksep(i)
ccckkk        enddo
ccckkk        write(91,*)
ccckkk        write(91,*)
ccckkk        write(92,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(92,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(92,*)
ccckkk        write(92,*)
ccckkk        write(93,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(93,*) ynsl(i),znsl(i)
ccckkk        enddo
ccckkk        write(93,*)
ccckkk        write(93,*)
ccckkk        write(94,*) '# jt ',jt
ccckkk        do i=1,ngo1     
ccckkk          write(94,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
ccckkk        enddo
ccckkk        write(94,*)
ccckkk        write(94,*)
*- trattamento getto 
        if(kget.eq.1)then
            write(*,*) ' AAA Chiamo get 2 '
            call  get(ng,npc,ycsl,zcsl,ycn,zcn,yn,zn,  
     #               xigs,zegs,xigb,zegb,xigf,zegf) 
        endif
        nng = kget*ng
* - chiamo il solutore:
* - prima ricalcolo dphtsl (non so se serve, dipende dalle BCs
        do i=1,npsl
          if(i.eq.1)then
            tf = sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
            dff= (phisl(i+1)-phisl(i))/tf
            dphtsl(i) = dff
          elseif(i.eq.npsl)then
            tb = sqrt( (ycsl(i-1)-ycsl(i))**2 + (zcsl(i-1)-zcsl(i))**2 )
            dfb= (phisl(i)-phisl(i-1))/tb
            dphtsl(i) = dfb
          else
            tf = sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
            tb = sqrt( (ycsl(i-1)-ycsl(i))**2 + (zcsl(i-1)-zcsl(i))**2 )
            dff= (phisl(i+1)-phisl(i))/tf
            dfb= (phisl(i)-phisl(i-1))/tb
            dphtsl(i) = 0.5d0*(dff+dfb)
          endif
        enddo
        if(jt.le.jjget)then
          call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,
     #                  phi,dphi,phisl,dphisl,jt)
        else
            write(*,*) 'AAA chiamo solv22 - 2'
          call solv22(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,
     #           ycsl,zcsl,dphi,phisl,dphtsl,dphtbsl,phb,
     #           xigs,zegs,xigb,zegb,xigf,zegf,
     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
     #           ph,dphn,dphnb,a1,b1,c1,d1,e1,phi,rl,
     #           mb,mf,mt,m,n,nt,ntt,
     #           xj,ze,xis,zes,jt,ksep,kse,kord,kor)
          call calsol(mb,mf,mt,m,n,nt,ntt,phi,phib,phb,dphisl,
     #                  dphnb,rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,
     #                 ry,rz,kse,dphi)
        endif
*- velocita   
        call calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,
     #                 phi,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phisl,
     #                 rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,
     #                 amp,vym1,vzm1,rny,rnz,tmy,tmz,
     #                 ampsl,vymsl1,vzmsl1,rnysl,rnzsl,tmysl,tmzsl,vxi,
     #                 ry,rz,ty,tz,kse)
*
ccckkk        write(10,*) '# jt ',jt,ng
ccckkk        do i=1,npc+ng
ccckkk          write(10,'(4d16.8,i4)') yce(i),zce(i),phi(i),dphi(i),ksep(i)
ccckkk        enddo
ccckkk        write(10,*)
ccckkk        write(10,*)
ccckkk        write(11,*) '# jt ',jt
ccckkk        do i=1,npc+ng+1
ccckkk          write(11,*) yv(i),zv(i),ksep(i)
ccckkk        enddo
ccckkk        write(11,*)
ccckkk        write(11,*)
ccckkk        write(12,*) '# jt ',jt,ng
ccckkk        do i=1,npsl  
ccckkk          write(12,'(4d16.8)') ycsl(i),zcsl(i),phisl(i),dphisl(i)
ccckkk        enddo
ccckkk        write(12,*)
ccckkk        write(12,*)
ccckkk        write(13,*) '# jt ',jt
ccckkk        do i=1,npsl+1   
ccckkk          write(13,*) ysl(i),zsl(i)
ccckkk        enddo
ccckkk        write(13,*)
ccckkk        write(13,*)
*
* -------- LUNGHEZZA STRISCIA DI CONTROLLO -----
*
        if(jend.eq.0)then
          rrl=0.d0
          do i=jind+1,jfid-1
            dl= sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
            rrl= rrl + dl
          enddo
          dlin= (1.d0-tin)*sqrt( (ycsl(jind+1)-ycsl(jind))**2 +
     #                    (zcsl(jind+1)-zcsl(jind))**2 )
          dlfi= tfi*sqrt( (ycsl(jfid+1)-ycsl(jfid))**2 +
     #                    (zcsl(jfid+1)-zcsl(jfid))**2 )
          rrl= rrl + dlin + dlfi
          yin = ycsl(jind)+tin*(ycsl(jind+1)-ycsl(jind))
          zin = zcsl(jind)+tin*(zcsl(jind+1)-zcsl(jind))
          yfi = ycsl(jfid)+tfi*(ycsl(jfid+1)-ycsl(jfid))
          zfi = zcsl(jfid)+tfi*(zcsl(jfid+1)-zcsl(jfid))
          write(77,'(i4,d15.6,d25.15,2i4,4d15.6)') jt,t,rrl,jind,jfid+1,
     #               yin,zin,yfi,zfi                             
        endif
* ---------
*
*
*
* - PRESSIONE ------------------------------------------------------ 
* -- calcolo derivata tangenziale della velocita tangenziale
      nng   = ng*kget

* !! dphtu is the tangential velocity at the apex. 
* !! ATTENZIONE: la variabile 'san" non è inizializzata. Deve essere
* l'angolo limite all'apice del corpo. Definirlo nella routine initial
* (o in ridis se si vuole variarlo nel tempo) in base alla geometria

      write(*,*) ' ATTENZIONE - SISTEMARE variabile san !!!'
      dphtu    = -vfall*san
      dphtd    = (phi(2)-phi(1))/(0.5d0*(amp(1)+amp(2)))
      dpht2(1) = ( dphtd - dphtu )/amp(1)
      dpnt(1)  = dpht2(1)*dphi(1)
      dpntt(1)  = dpht2(1)*dphi(1)
      i = 1 
      do i = 2,npc+nng-1
        amu = 0.5d0*(amp(i-1)+amp(i))
        amd = 0.5d0*(amp(i+1)+amp(i))
        dpht2u    = (dpht(i)-dpht(i-1))/amu
        dpht2d    = (dpht(i+1)-dpht(i))/amd
        dpht2(i)  = 0.5d0*(dpht2u+dpht2d)
        vxiu      = (vxi(i)-vxi(i-1))/amu
        vxid      = (vxi(i+1)-vxi(i))/amd
        dpht22    = 0.5d0*(vxiu+vxid)
        dpnt(i)   = dpht2(i)*dphi(i)
        if(i.ge.npc+1)then 
         dpntt(i)   = dpht22*dphi(i)
         else
         dpntt(i)   = dpnt(i)
        endif  
      enddo     
      dpht2(npc+nng)=(dpht(npc+nng)-dpht(npc+nng-1))/
     #               (0.5d0*(amp(npc+nng)+amp(npc+nng-1)))
      dpnt(npc+nng) = dpht2(npc+nng)*dphi(npc+nng) 
      dpntt(npc+nng)= dpnt(npc+nng) 
* -- calcolo dphi/dt sulla SL
      do i = 1,npsl
        dptsl(i) = -(dphisl(i)**2 + dphtsl(i)**2)/2.d0
      enddo
      do i = 1,npc+nng 
        dpt(i)  = -(dphi(i)**2 + dpht(i)**2)/2.d0
        dpt2(i) = dpt(i)
      enddo
*
      call bcpre(ng,kget,vfall,npc,npsl,tmz,amp,phi,dphi,dpht,
     #                 vxi,dphisl,dphtsl,yce,zce,dpnt,dpntt,dptsl,jt,tn,
     #                 ygn,zgn,ygs2,zgs2,tg,ngo1,tc)
*
      if(kget.eq.0)then
        call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,
     #                  dpt2,dpnt,dptsl,dpntsl,jt)
      else
        write(*,*) ' ATTENZIONE: controllare dpttsl,dpttsl nella
     #              chiamata alla SOLV22p   !!!!! '
        call solv22p(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,
     #              ycsl,zcsl,dpnt,dptsl,dpttsl,dpttsl,dpb,
     #              xigs,zegs,xigb,zegb,xigf,zegf,
     #              ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
     #              ph,dphn,dphnb,a2,b2,c2,d2,e2,dpt,
     #              rl,mb,mf,mt,m,n,nt,ntt,xj,ze,xis,
     #              zes,jt,ksep,kse,kord,kor)
*
        call calsolp(mb,mf,mt,m,n,nt,ntt,dpt,dptb,dpb,dpntsl,dphnb,
     #             rl,a2,b2,c2,d2,e2,xj,ze,xis,zes,ry,rz,kse,dpnt)
        call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,
     #                  dpt2,dpnt,dptsl,dpntsl,jt)
      endif
*
* - calcolo con Dphi/Dt
      call prefin(npc,nng,npco,nngo,vfall,dt,dpht,dphi,tmy,tmz,
     #                  rny,rnz,vyo,vzo,yce,zce,yco,zco,phi,phio,
     #                  vym1,vzm1,pres)
* -- calcolo pressione e forza
      ff  = 0.d0
      ff2 = 0.d0
      do i=1,npc+nng
        p1 = -dpt(i)
        p11 = -dpt2(i)
        p2 = -(dphi(i)**2+dpht(i)**2)/2.d0
        pre(i) = p1+p2
        pre2(i) = p11+p2
        ff = ff + pre(i)*amp(i)*rnz(i)
        ff2 = ff2 + pres(i)*amp(i)*rnz(i)
c        write(41,'(i6,6d15.6)') i,yce(i),zce(i),p1,-p2,pre(i),pre2(i)
c        if(i.eq.npc) write(41,*)
      enddo
c      write(41,*)
c      write(41,*)
      ff  = -2.d0*ff
      ff2 = -2.d0*ff2
      write(99,'(i10,3d15.7)') jt,t,ff,ff2
* FINE PRESSIONE -------------------------------------------------------
* - stampo
      if(mod(jt,ksta).eq.0.or.jt.le.200)then
      call stampa(vfall,ng,kget,npc,npsl,jt,t,dt,frdt,llf,
     #                  scon,svel,spot,spre,phi,dphi,phisl,dphisl, 
     #                 yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,
     #                 vym1,vzm1,vymsl1,vzmsl1,
     #                 dpht,dphtsl,dpt,dpnt,pre,pre2,pres,vxi,
     #                 jend,jind,jfid,yin,zin,yfi,zfi)

      endif
*
      nng = ng*kget
c      if(jt.eq.jjget+1)stop
* ------------------------------------  CONTROLLI
cc - B+
c      write(31,*) '# ',jt,zv(1)
c      do i=1,npc+1
c        write(31,'(2d15.7)') yv(i),zv(i)
c      enddo
c      write(31,*)  
c      write(31,*)  
c      write(32,*) '# ',jt
c      do i=1,npc
c        write(32,'(8d15.7)') yce(i),
c     #       zce(i),phi(i),dphi(i),
c     #            dpht(i),vym2(i),vzm2(i) ,amp(i)  
c      enddo
c      write(32,*)
c      do i =1,nng
c        write(32,'(8d15.7)') 0.5d0*(ygb(i)+ygb(i+1)),
c     #        0.5d0*(zgb(i)+zgb(i+1)),phi(npc+i),dphi(npc+i),
c     #            dpht(npc+i),vym2(npc+i),vzm2(npc+i) ,amp(npc+i)  
c      enddo
c      write(32,*)  
c      write(32,*)  
c      write(33,*) '# ',jt
c      do i=1,npsl+1
c        write(33,'(4d15.7)') ysl(i),zsl(i)
c      enddo
c      write(33,*)  
c      write(33,*)  
c      write(34,*) '# ',jt
c      do i=1,npsl
c        write(34,'(8d15.7)') 0.5d0*(ysl(i)+ysl(i+1)),
c     #        0.5d0*(zsl(i)+zsl(i+1)),phisl(i),dphisl(i),
c     #            dphtsl(i),vymsl2(i),vzmsl2(i) ,ampsl(i)  
c      enddo
c      write(34,*)  
c      write(34,*)
*-----------------------------------------------  FINE CONTROLLI
      enddo
*
      stop
      end
