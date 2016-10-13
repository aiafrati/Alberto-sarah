c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine initial(frin,frfi,jind,tin,jfid,tfi)

      include"slam_p.h"
      include"slam_v.h"
      include"slam_f.h"
*
      write(*,*) '.......... :  initial'
c     - inizializzazione variabili

      pi    = acos(-1.d0)
      ande  = alfa*pi/180.d0
      call filiniz
*
      frin  = 0.015d0
      frfi  = 0.025d0
      if (lrest.ne.1) then
*
c       - inizializzazione variabili per integrazione temporale
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
*        pro0  = -ampp*sin(ande)/rappi
*
        if(iwig.eq.2)then
* ----- GENERIC SHIP -------------------------
*        ux=1.d0
        alpha = alpha/180.d0*pi
        write(*,*) 'alpha ', alpha/pi*180.d0
        call geoini(alpha,ztr,wl,x00,dztr,al,bl,tl)
        write(*,'(a9,3d15.5)') 'al bl tl ' ,al,bl,tl
        write(*,'(a22,4d15.5)') 'alpha, ztr, xs00, dztr,wl',
     #                       alpha,ztr,x00,dztr,wl
*        xs0 = -pro0
*        dx0   = -pro0/tan(alpha)
        xs0 =  x00 + dx0
        tend = (wl - dx0)/ux
        write(*,*) 'iwig= 2, xs0 , tend',xs0,tend
        call geopla(xs0,dztr,alpha,ygn,zgn,rnxg,ng,al,bl,tl)
*
        tag(1)=0.d0
        do ip=2,ng
           dy = ygn(ip)-ygn(ip-1)
           dz = zgn(ip)-zgn(ip-1)
           tag(ip)=tag(ip-1)+sqrt(dy**2+dz**2)
        end do
*  ---- taglio la geometria a z = 0. 
        cm0 = 0.d0
        cq0 = 0.d0 
        call izero(kint,yi,zi,ygn,zgn,1,ng-1,cm0,cq0)
        dyi  = yi-ygn(kint)
        dzi  = zi-zgn(kint)
        tai  = tag(kint) + sqrt( dyi**2 + dzi**2 )
        spo0 = yi
        yp1 = 1.d31
        ypn = 1.d31
        call splone(ygs2,tag,ygn,yp1,ypn,ng,npamx)
        call splone(zgs2,tag,zgn,yp1,ypn,ng,npamx)
        call splone(rnxgs2,tag,rnxg,yp1,ypn,ng,npamx)
* ----- GENERIC SHIP END ---------------------        
        elseif(iwig.eq.1)then
* ----- WIGLEY SHIP --------------------------
*    - ATTENZIONE : vale anche per savander,
*      scommentando le righe zmaxs =...,,  ztr0 = ..., if(xi0...) ,
*                            wl = ... , *      call savan ....
*
*        ux   =  1.d0
*  - misure wigley
*        al   =  .8d0
*        bl   =  .16d0
*        tl   =  .04d0
* - misure savander
         al = 7.5d0
         bl = 2.d0
         tl = 0.d0
        alpha = alpha/180.d0*pi
*        alpha= 6.d0/180.d0*pi
        ztr0 = -al/2.d0*sin(alpha) - tl*cos(alpha)
        zmaxs = 0.5d0
        ztr0 = -al*sin(alpha) 
        dztr = -(ztr - abs(ztr0))
        write(*,*) 'p ', ztr, ztr0 ,dztr,alpha
*        x000 = -tl*cos(alpha)**2/sin(alpha) - tl*sin(alpha)
        xi0  = -tl/tan(alpha) + dztr/sin(alpha)
        x00  = xi0*cos(alpha) - tl*sin(alpha)
       if(xi0.lt.-0.5d0*al)STOP '...initial: STOP, scafo poco assettato'
       if(xi0.lt.-0.d0)STOP '...initial: STOP, scafo poco assettato'
        wl   = 0.5d0*al*cos(alpha) +  tl*cos(alpha)**2/sin(alpha)
     #         - dztr/sin(alpha)*cos(alpha)
        wl   = al*cos(alpha) - dztr/sin(alpha)*cos(alpha)
*        aa   = -pro0/tan(alpha)
*        xs0  = x00 + aa  
*        dx0   = -pro0/tan(alpha)
        xs0  = x00 + dx0  
        tend = (wl - dx0)/ux
        write(*,*) 'iwig =1 , xs0, tend',xs0,tend,alpha/pi*180.,dztr
*        call wigl(alpha,al,bl,tl,xs0,ygn,zgn,rnxg,ng,dztr)
       write(*,*) 'ppppp ' , al , bl , alpha 
       write(*,*) 'ppppp ' , dztr , zmaxs , xs0 
       call savan(al,bl,alpha,dztr,zmaxs,xs0,ygn,zgn,rnxg,ng)
        tag(1)=0.d0
        do ip=2,ng
           dy = ygn(ip)-ygn(ip-1)
           dz = zgn(ip)-zgn(ip-1)
           tag(ip)=tag(ip-1)+sqrt(dy**2+dz**2)
        end do
*  ---  taglio la geometria a z = 0. 
        cm0 = 0.d0
        cq0 = 0.d0 
        call izero(kint,yi,zi,ygn,zgn,1,ng-1,cm0,cq0)
        dyi  = yi-ygn(kint)
        dzi  = zi-zgn(kint)
        tai  = tag(kint) + sqrt( dyi**2 + dzi**2 )
        spo0 = yi
        yp1 = 1.d31
        ypn = 1.d31
        write(*,*) 'geommmm' 
        do j = 1,ng+1
          write(*,'(i4,3d15.7)') j, ygn(j),zgn(j),rnxg(j)
        end do
        call splone(ygs2,tag,ygn,yp1,ypn,ng,npamx)
        call splone(zgs2,tag,zgn,yp1,ypn,ng,npamx)
        call splone(rnxgs2,tag,rnxg,yp1,ypn,ng,npamx)
*------ WIGLEY SHIP end--------------
*
        elseif(iwig.eq.0) then
*
*
* --- lettura geometria corpo. il vertice e' a y=0, z=0
*
        open(unit=7,file='geo.in',status='unknown')
        read(7,*) ngo
        do ig=1,ngo
          read(7,*) yg(ig), zg(ig)
        enddo
        close(7)
*
*   --- tend : considero un wetting factor  limite > 1.3 o 2.d0 
        ztot = zg(ngo) - zg(1)
        zm0  = -pro0
        if(i2d.eq.0) then
           wf = 1.3d0
        else
           wf = 2.d0 
        endif
        zt   = ztot/wf
        zt2  = zt - zm0
        tend = zt2/vfall0
        write(*,*) '   tend  <= ' ,tend  
*
**  --- interpolazione spline della geometria iniziale (-> ygs2,zgs2)
**
** ---   calcolo ascissa curvilinea ta
**
        ta(1) = 0.d0
        do 20, ip=2,ngo
           dy = yg(ip) - yg(ip-1)
           dz = zg(ip) - zg(ip-1)
           dl = sqrt(dy**2 + dz**2)
           ta(ip) = ta(ip-1) + dl
  20    continue
**
** ---   interpolazione spline su ng punti
**
        yp1 = 1.d31
        ypn = 1.d31
        call splone(ygs2,ta,yg,yp1,ypn,ngo,npamx)
        call splone(zgs2,ta,zg,yp1,ypn,ngo,npamx)
**
*c       - discretizzazione iniziale della superficie del corpo immerso
**
**  ---- infittisco gli offset iniziali (ngo), mettendone ng.
*c      ng=201
**
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
           if(ig.gt.1)then 
            tag(ig)=tag(ig-1)+sqrt((yy-ygn(ig-1))**2+(zz-zgn(ig-1))**2)
           endif
        end do
        do i=1,ng
          yg(i) = ygn(i)
          zg(i) = zgn(i)
          ta(i) = tag(i)
        end do
        yp1 = 1.d31
        ypn = 1.d31
        call splone(ygs2,ta,yg,yp1,ypn,ng,npamx)
        call splone(zgs2,ta,zg,yp1,ypn,ng,npamx)
**
**  ---  taglio la geometria a z = -pro0 
        cm0 = 0.d0
        cq0 = -pro0 
        call izero(kint,yi,zi,ygn,zgn,1,ng-1,cm0,cq0)
        dyi  = yi-ygn(kint)
        dzi  = zi-zgn(kint)
        tai  = tag(kint) + sqrt( dyi**2 + dzi**2 )
        spo0 = yi
*
** -- sposto ygn a pro0
*
        do i=1,ng
          zgn(i)=zgn(i)+pro0
        end do
*
*
        endif
*
**
* -- genero la griglia sulla parte tagliata della geometria e la traslo
*      di pro0
*      la posizione  r = float(iv)/float(npci+1)*tai e' come 
*      in Code (il primo
*      pannello viene lungo doppio cosi)
*
        npci = 6 
        r    = 0.d0
        dr   = tai/float(npci)
        do 30, iv = 1,npci+1
*           r = r + dr
           r = float(iv)/float(npci+1)*tai
           if (iv.eq.1)      r = 0.d0
           if (iv.eq.npci+1) r = tai
           call splont(r,yy,tag,ygn,ygs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           call splont(r,zz,tag,zgn,zgs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           if(iwig.gt.0)then
           call splont(r,rnx,tag,rnxg,rnxgs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           endif
           yv(iv) = yy
           zv(iv)   = zz 
           rnxx(iv) = rnx
  30    continue        
        npc = npci
*
c       - discretizzazione iniziale della superficie libera
*
        npr  = log(1.d0+(estr-spo0)/ampp*(escr-1.d0))/log(escr) + 1
        ampl = (estr-spo0)*(1.d0-escr)/(1.d0-escr**npr)
        do i = npc+1,npc+npr
          yv(i+1) = yv(i) + ampl
          zv(i+1) = 0.d0
          ampl    = ampl*escr
        enddo
*
        npt = npc + npr 
        if (npt.gt.npamx) stop 'aumenta npamx'
*  
c       - inizializzo potenziale di perturbazione e ampiezze dei 
c         pannelli sulla sup. libera

        do i = 1,npt
          yce(i) = (yv(i+1)+yv(i))/2.d0
          zce(i) = (zv(i+1)+zv(i))/2.d0
          amx    =  yv(i+1) - yv(i)
          amz    =  zv(i+1) - zv(i)
          amp(i) = sqrt(amx*amx+amz*amz)
          rnx    =  amz/amp(i)
          rnz    = -amx/amp(i)
          if (i.le.npc) then
            kphi(i)= 0
            if(iwig.gt.0)then
            rnxm   = 0.5d0*(rnxx(i)+rnxx(i+1))
            dphi(i) = -ux*rnxm/sqrt(1.d0 - rnxm**2)
            else
            dphi(i)= -vfall0*rnz
            dphi2  = vfall*cos(ande)
            endif
          else
            kphi(i)= 1
            phi(i) = 0.d0
          end if 
        end do
        volb0 = 0.d0 
        volj  = 0.d0
        do i = 1,npt
          volb0 = volb0 - (yv(i+1)-yv(i))*(zv(i+1)+zv(i))/2.d0
        end do
        cmm0 = -volb0
*
* ---  striscia di SL da monitorare
*
        yin = frin*estr
        yfi = frfi*estr
        do j=npc+1,npt
          if(yin.lt.yce(j))goto 98
        enddo
 98     jind = j-1
        do j=npc+1,npt
          if(yfi.lt.yce(j))goto 99
        enddo
 99     jfid = j-1
        tin  = (yin-yce(jind))/(yce(jind+1)-yce(jind))
        tfi  = (yfi-yce(jfid))/(yce(jfid+1)-yce(jfid))
*
*
      else

        stop ' OPZIONE RESTART NON DISPONIBILE '

      end if

      return
      end

