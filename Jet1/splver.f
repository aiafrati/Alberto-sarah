
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine splver (iint,xg,zg,ng,xc,zc,xp,zp,npc,npt,proat,
     #                   ande,kcut,ampk,tag,xgs2,zgs2)
*
c     - questa routine determina le posizioni dei vertici tramite
c       interpolazione spline dei centroidi.
*
      include"slam_p.h"
      dimension xp(npamx+1),zp(npamx+1)
      dimension xc(npamx),zc(npamx),xx(npamx),zz(npamx)
      dimension tp(npamx+1),xs2(npamx),zs2(npamx+1)
      dimension xg(npamx+1),zg(npamx+1)
      dimension xgs2(npamx+1),zgs2(npamx+1),tag(npamx+1)
*
      write(*,*) '.......... :  splver'
*
      epsp = 1.d-10 
      epst = 0.1d0
*
c    - definizione della ascissa curvilinea per le spline
*
      write(*,*) '   kcut = ' ,kcut
      tp(1) = 0.d0
      xx(1) = xc(npc+kcut+1)
      zz(1) = zc(npc+kcut+1)
      do i = 1,npt-npc-kcut-1
        ax      = (xc(i+1+npc+kcut)-xc(i+npc+kcut))
        az      = (zc(i+1+npc+kcut)-zc(i+npc+kcut))
        dis     = sqrt(ax*ax+az*az)
        tp(i+1) = tp(i) + dis
        xx(i+1) = xc(i+1+npc+kcut)
        zz(i+1) = zc(i+1+npc+kcut)
      end do
* - interpolazione spline
      yp1 = 1.d+31
      ypn = 1.d+31
      call splone(xs2,tp,xx,yp1,ypn,npt-npc-kcut,npamx)
      call splone(zs2,tp,zz,yp1,ypn,npt-npc-kcut,npamx)
      do i = npc+kcut+2,npt
        tts = (tp(i-npc-kcut)+tp(i-npc-kcut-1))/2.d0
        call splont(tts,xp(i),tp,xx,xs2,npt-npc-kcut,npamx,
     #             xnull,0,xnull,0,xnull,0)
        call splont(tts,zp(i),tp,zz,zs2,npt-npc-kcut,npamx,
     #             xnull,0,xnull,0,xnull,0)
      end do

c     - trovo il nodo npc+1 con l'intersezione retta-retta
*
* -- trovo pannello e punto di intersezione corpo  con retta
*   ( xs,zs,wsx,wsz) : trovo cioe' il nodo npc+1
*
      cms = (zp(npc+kcut+2)-zc(npc+kcut+1))/
     &      (xp(npc+kcut+2)-xc(npc+kcut+1))
      cqs = zc(npc+kcut+1) - cms*xc(npc+kcut+1)
      wsx = xp(npc+kcut+2) - xc(npc+kcut+1)
      wsz = zp(npc+kcut+2) - zc(npc+kcut+1)
      xs  = xc(npc+kcut+1)
      zs  = zc(npc+kcut+1)   
      x1 = xc(npc+1)
      z1 = zc(npc+1)
*
        do i=1,ng-1
           w2x = xg(i+1) - xg(i)
           w2z = zg(i+1) - zg(i)
           x2  = xg(i)
           z2  = zg(i)
           if (kcut.eq.0) then
             w1x = wsx
             w1z = wsz
           else
             w1x = zg(i+1) - zg(i)
             w1z =-( xg(i+1) - xg(i) )
           endif
           call intret(xi,zi,r1,r2,x1,z1,w1x,w1z,x2,z2,w2x,w2z)
           if(r2.ge.0.d0.and.r2.le.1.d0)then
             iint = i
             write(*,*) 'iint r2 ',iint,r2
             if(kcut.gt.0)then
*     ora trovo il nodo npc+2
               call intret(xin,zin,r1,r2,x1,z1,w1x,w1z,xs,zs,wsx,wsz)
               xp(npc+2) = xin
               zp(npc+2) = zin
             endif
*
* -- intersezione spline corpo con retta w1x,w1z:
           dt   = tag(i+1)-tag(i)
           xt1  = xg(i)
           zt1  = zg(i)
           xt2  = xg(i+1)
           zt2  = zg(i+1)
           xd1  = xgs2(i)
           zd1  = zgs2(i)
           xd2  = xgs2(i+1)
           zd2  = zgs2(i+1)
           f1   = w1x*zt1 - w1z*xt1
           f2   = w1x*zt2 - w1z*xt2
           fs1  = w1x*zd1 - w1z*xd1
           fs2  = w1x*zd2 - w1z*xd2
           fc   = w1x*z1  - w1z*x1
           den  = fs1 - fs2 
           dd   = abs(1.d0/6.d0*(dt**2)*den) 
*   vedo se il coeff. di x**3 e' < eps
           if(dd.le.epsp)then
           write(*,*) 'dd eps piccolo'
             a = (dt**2)*fs2/2.d0
             b = f1 - f2 - (dt**2)*fs1/6.d0 - (dt**2)*fs2/3.d0 
             c = f2 - fc
             if(abs(a).le.epsp)then
           write(*,*) 'aaa! eps piccolo'
               if(abs(b).le.epsp)then
                 write(*,*) 'ATTENZIONE!!! x=-c/b, b molto piccolo !!!'
               endif 
                 tt1 = -c/b
                 tt2 = 1.d31
                 tt3 = 1.d31
             else
               delta = b**2 - 4.d0*a*c
           write(*,'(a9,4d15.6)') 'abcdelta ', a,b,c,delta
               if(delta.lt.0.d0)then
                 tt1 = 1.d31
                 tt2 = 1.d31
                 tt3 = 1.d31
               else
                 tt1 = (-b + sqrt(delta))/(2.d0*a)
                 tt2 = (-b - sqrt(delta))/(2.d0*a)
                 tt3 = 1.d31
               endif
             endif
           else
           write(*,*) 'dd eps grande'
             a = 3.d0*fs2/den
             b = ( 6.d0/(dt**2)*(f1-f2) - (fs1+2.d0*fs2) )/(fs1-fs2) 
             c = 6.d0/(dt**2)*(f2-fc)/(fs1-fs2)
             call poli3(a,b,c,tt1,tt2,tt3)
           endif
           write(*,*) 'tt1 .. ',tt1,tt2,tt3
           im = 0
           if(tt1.ge.-epst.and.tt1.le.(1.d0+epst)) then
             ttt = tt1
             im  = im + 1
           endif
           if(tt2.ge.-epst.and.tt2.le.(1.d0+epst)) then
             im  = im + 1
             ttt = tt2
           endif
           if(tt3.ge.-epst.and.tt3.le.(1.d0+epst)) then
             im  = im + 1
             ttt = tt3
           endif
           if(im.eq.0) write(*,*) 'splver, intersezione non trovata'
           if(im.gt.1) write(*,*) 'splver, intersezione doppia '
*
           c  = (ttt**3 - ttt)*(dt**2)/6.d0  
           d  = (-(ttt**3) + 3.d0*(ttt**2) - 2.d0*ttt)*(dt**2)/6.d0
           xip = ttt*xt1 + (1.d0-ttt)*xt2 + c*xd1 + d*xd2
           zip = ttt*zt1 + (1.d0-ttt)*zt2 + c*zd1 + d*zd2
*
*  -------------------------------------
             goto 999
           endif
        end do
  999   continue     
        write(*,*) 'iint  xx ' ,iint
        write(*,'(2d25.10)')  xi,zi
        write(*,'(2d25.10)')  xip,zip
        xi = xip
        zi = zip
        xp(npc+1) = xi
        zp(npc+1) = zi
*
c       - verifica
*
      if(kcut.eq.0)then
*
        if (zp(npc+1).lt.proat) stop ' intersezione non fisica .. '
*
*        diff2 = zc(npc+1) - (cm1*xc(npc+1)+cq1)
        ck    = xc(npc+1)
        call izero2(kint,xi,zi,xg,zg,1,ng,ck)
        if (kint.eq.-999) then
           write(*,*) '   OKKIO!, izero2 in splver non trova intersez.'
           write(*,*) '   corpo troppo verticale? '
           write(*,*) '   SL troppo lontano?'
        endif
        diff  = zc(npc+1) - zi 
        if (diff.ge.0.d0) write(*,*) '   !!! CENTROIDE DENTRO CORPO !!!'
*
        else
* ----------------------------------------------------------
*
*
        amtt = sqrt((xin-xp(npc+1))**2+(zin-zp(npc+1))**2)
*
*        if (amtt.gt.ampk) then
          amcc = ampk
*        else
          amcc = amtt
          ampk = amcc
*        end if
*
*        ampy       = xp(npc+1) - xp(npc)
*        ampz       = zp(npc+1) - zp(npc)
*
*        am         = sqrt(ampy**2 + ampz**2)
*        rny        = ampz/am
*        rnz        = -ampy/am
*        xp(npc+2)  = xp(npc+1) + amcc*rny
*        zp(npc+2)  = zp(npc+1) + amcc*rnz
*       
      end if
*
      xp(npt+1)  = xc(npt) + (xc(npt)-xp(npt))
      zp(npt+1)  = zc(npt) + (zc(npt)-zp(npt))
*
      return
      end
*
