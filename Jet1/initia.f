
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine initia
*
      include"slam_p.h"
      include"slam_v.h"
      include"slam_f.h"
*
c     - inizializzazione variabili
*
      pi    = acos(-1.d0)
      ande  = alfa*pi/180.d0
      call filiniz
*
      if (lrest.ne.1) then
*
c       - inizializzazione variabili per inteagrazione temporale
*
        t      = 0.d0 
        tust   = 0.d0
        amplim = ampp/pfraz
        vfall  = vfall0
        vsopr  = vfall0
        dzet   = 0.d0
        llf    = 0
        kcut   = 0
        nppr   = 0
*
c       - parametri per la dinamica
*
        weir = rmu*(vfall0**2/9.81)**2
        rkk  = vfall0**2/vea
*
*        spo0  =  ampp*cos(ande)/rappi
*
*--------- per ora lascio cosi pro0 anche se ande non ha + senso
*
        pro0  = -ampp*sin(ande)/rappi
c        pro0  = -0.001
*
*  ---- interpolazione spline della geometria iniziale (-> ygs2,zgs2)
*
* ---   calcolo ascissa curvilinea ta
*
      write(*,*) 'geometria input'
      ta(1) = 0.d0
      do 20, ip=2,ngo
         dy = yg(ip) - yg(ip-1)
         dz = zg(ip) - zg(ip-1)
         dl = sqrt(dy**2 + dz**2)
         ta(ip) = ta(ip-1) + dl
  20  continue
*
* ---   interpolazione spline su ng punti
*
      yp1 = 1.d31
      ypn = 1.d31
      call splone(ygs2,ta,yg,yp1,ypn,ngo,npamx)
      call splone(zgs2,ta,zg,yp1,ypn,ngo,npamx)
*
* --- infittisco gli offset iniziali (ngo), mettendone ng
c      ng=201
      r=0.d0
      dr=ta(ngo)/float(ng-1)
      tag(1)=0.d0       
      do ig=1,ng
           r = r+dr
           if(ig.eq.1)  r=0.d0
           if(ig.eq.ng) r=ta(ngo)
           call splont(r,yy,ta,yg,ygs2,ngo,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           call splont(r,zz,ta,zg,zgs2,ngo,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           ygn(ig) = yy
           zgn(ig) = zz
	   write(*,*) yy,zz
        if(ig.gt.1)then 
        tag(ig)=tag(ig-1)+sqrt((yy-ygn(ig-1))**2+(zz-zgn(ig-1))**2)     
        endif
      end do
      do i=1,ng
        yg(i)=ygn(i)
        zg(i)=zgn(i)
        ta(i)=tag(i)
      end do
      yp1 = 1.d31
      ypn = 1.d31
      call splone(ygs2,ta,yg,yp1,ypn,ng,npamx)
      call splone(zgs2,ta,zg,yp1,ypn,ng,npamx)
*
* -------------------------------------------------------
*     sposto ygn e zgn a pro0
      do i =1,ng
c         ygn(i) = ygn(i)
         zgn(i) = zgn(i) + pro0
      end do
*
* - leggo  config SL e potenziale
*
      open(unit=40,file='conf.inp',status='unknown')
*      read(40,*) ampk	
      read(40,*) nsl
      do iv=1,nsl+1
         read(40,*) ysl(iv),zsl(iv)
         write(*,*) ysl(iv),zsl(iv)
      end do
      write(*,*)
      write(*,*)
      close(40)
      open(unit=40,file='phi.inp',status='unknown')
      read(40,*) nsl
      do iv=1,nsl
         read(40,*) ycsl(iv),zcsl(iv),phisl(iv),xnull 
      end do
      close(40)
* - dimensionalizzo
      dim1  = abs(pro0)
      dim2  = abs(pro0*vfall0)
      do iv =1,nsl
         ysl(iv)   = ysl(iv)*dim1
         zsl(iv)   = zsl(iv)*dim1 
         ycsl(iv)   = ycsl(iv)*dim1
         zcsl(iv)   = zcsl(iv)*dim1 
         phisl(iv) = phisl(iv)*dim2
         write(*,*) ysl(iv),zsl(iv)
      end do
      write(*,*)
      write(*,*)
      ysl(nsl+1) = ysl(nsl+1)*dim1
      zsl(nsl+1) = zsl(nsl+1)*dim1
* -taglio o aggiungo pezzo
      if (ysl(nsl+1).gt.estr) then
        do i=1,nsl
          if(ysl(i).gt.estr)then 
            icut = i-1
            goto 222
          endif
        enddo
 222    continue
        nsl = icut
      else
*        yl  = ysl(nsl)
*        yr  = estr
*        am = ysl(nsl)-ysl(nsl-1)
*        npr  = log(1.d0+(yr-yl)/am*(escr-1.d0))/log(escr) + 1
*        ampl = (yr-yl)*(1.d0-escr)/(1.d0-escr**npr)
*        do i = 1,npr
*          ysl(nsl+i)   = ysl(nsl+i-1) + ampl
*          zsl(nsl+i)   = 0.d0
*          phisl(nsl+i) = 0.d0
*          ampl         = ampl*escr
*        end do
*        nsl = nsl+npr-1
         nsl = nsl
      endif
* 
* - punto di attacco della SL sul corpo (zv(npc+1))
      cm = 0.d0
      cq = zsl(1) 
      call izero(kint,yi,zi,ygn,zgn,1,ng,cm,cq)
      write(*,*) 'yi zi ',yi,zi
      tas = tag(kint) + sqrt( (zi-zgn(kint))**2 + (yi-ygn(kint))**2 )
c      call splont(tas,yi,ta,yg,ygs2,ng,npamx,
c     #                 xnull,0,xnull,0,xnull,0)
c      call splont(tas,zi,ta,zg,zgs2,ng,npamx,
c     #                 xnull,0,xnull,0,xnull,0)
c      zi = zi + pro0
*
* - traslo tutta la SL a partire dal primo pto sul corpo
*
      write(*,*) 'yi zi ',yi,zi
*
      dy = yi - ysl(1)
      do i =1,nsl+1
         ysl(i) = ysl(i)+dy
         write(*,*) ysl(i),zsl(i)
      end do
      write(*,*)
      write(*,*)
      do i =1,nsl
         ycsl(i) = ycsl(i)+dy
      end do

* - ricalcolo bene il tappo normale al corpo
c      yk = ysl(2)
c      zk = zsl(2)
c      call milan(kint,ygn,zgn,1,ng,yk,zk)
c      w1y = ygn(kint+1) - ygn(kint)
c      w1z = zgn(kint+1) - zgn(kint)
c      call intret(yi,zi,xnull,xnull2,ygn(kint),zgn(kint),w1y,w1z,
c     #            yk,zk,-w1z,w1y)
c      dyi  = yi-ygn(kint)
c      dzi  = zi-zgn(kint)
c      tai  = ta(kint)+sqrt(dyi**2+dzi**2)
c      tai = tas
      spo0 = yi
*
      ysl(1)  = yi
      zsl(1)  = zi
      ycsl(1) = 0.5*(ysl(1)+ysl(2)) 
      zcsl(1) = 0.5*(zsl(1)+zsl(2)) 
*
*
* - dimensione primo pannello corpo, calcolo npc e grigliatura corpo
      amii  = sqrt((ysl(2)-ysl(3))**2 + (zsl(2)-zsl(3))**2)
*
      tai=tas
*      
      amdi  = tai
      eskk  = escr
      if (eskk.gt.1.d0)then
        npci = int( log( 1.d0+(eskk-1.d0)*amdi/amii)/log(eskk) )
      else
        npci = int(amdi/amii)
      endif 
      if (npci.gt.1)then
        amii = amdi*(1.d0-eskk)/(1.d0-eskk**npci)
      endif  
*        
      amii1=amii
      r  = tai
      do iv =npci+1,1,-1
         r = r - amii
         if(iv.eq.1)      r = 0.d0
         if(iv.eq.npci+1) r = tai
         call splont(r,yyy,ta,yg,ygs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
         call splont(r,zzz,ta,zg,zgs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
         yv(iv)  = yyy
         zv(iv)  = zzz + pro0
         if(iv.le.npci) then
            ampy       = yv(iv+1) - yv(iv)
            ampz       = zv(iv+1) - zv(iv)
            amp(iv)    = sqrt(ampy**2 + ampz**2)
            rnz        = -ampy/amp(iv)
            dphi(iv)   = -vfall0*rnz
            kphi(iv)   = 0
         endif
         if (iv.le.npci) amii=amii*eskk
      end do 
      npc = npci
* - centroidi corpo
      do i = 1,npc
        yce(i) = 0.5*(yv(i+1)+yv(i))
        zce(i) = 0.5*(zv(i+1)+zv(i))
      end do
*
* - SL
*
        write(*,*)
        write(*,*)
        write(*,*) '# SL'
        do isl = 2,nsl
           yv(npc+isl)  = ysl(isl)
           zv(npc+isl)  = zsl(isl)
           yce(npc+isl-1)  = ycsl(isl-1)
           zce(npc+isl-1)  = zcsl(isl-1)
           phi(npc+isl) = phisl(isl)
           kphi(npc+isl)= 1
           amy  = (yv(npc+isl)-yv(npc+isl-1))
           amz  = (zv(npc+isl)-zv(npc+isl-1))
           amp(npc+isl-1) = sqrt(amy**2 + amz**2)
           write(*,*) yv(npc+isl),zv(npc+isl) 
        end do
        yv(npc+nsl+1) = ysl(nsl+1)
        zv(npc+nsl+1) = zsl(nsl+1)
        yce(npc+nsl)  = ycsl(nsl)
        zce(npc+nsl)  = zcsl(nsl)
        amy  = (yv(npc+nsl+1)-yv(npc+nsl))
        amz  = (zv(npc+nsl+1)-zv(npc+nsl))
        amp(npc+nsl) = sqrt(amy**2 + amz**2)
*        
        npt = npc + nsl  
        if (npt.gt.npamx) stop 'aumenta npamx'
*
        phi(npc+1)  = phisl(1)
        kphi(npc+1) = 2
        kcut        = 1


*---------------------------------------------------------
*       ampk
        ampk = sqrt((yv(npc+2)-yv(npc+1))**2+(zv(npc+2)-zv(npc+1))**2) 

*---------------------------------



        volb0 = 0.d0 
        volj  = 0.d0
        do i = 1,npt
          volb0 = volb0 - (yv(i+1)-yv(i))*(zv(i+1)+zv(i))/2.d0
        end do
        cmm0 = -volb0
      else

        stop ' OPZIONE RESTART NON DISPONIBILE '

      end if

      return
      end
*
