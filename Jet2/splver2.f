
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine splver2(xc,zc,xp,zp,npsl,proat,kget,
     #             estr,ng,xg,zg,xgs2,zgs2,tg,ngo,iint,jt,
     #                 ampli,ngo1,nsep,nsepo,xcn,zcn,npc,phg,phgs2,
     #                phin,ksep,di,ang,tc)
c     - questa routine determina le posizioni dei vertici tramite
c       interpolazione spline dei centroidi.

      include"slam_p.h"
      dimension xp(npamx),zp(npamx)
      dimension xc(npamx),zc(npamx),xx(npamx),zz(npamx)
      dimension tp(npamx),xs2(npamx),zs2(npamx)
      dimension xg(npamx),zg(npamx),xgs2(npamx),zgs2(npamx),tg(npamx) 
      dimension xcn(npamx),zcn(npamx),phg(npamx),phgs2(npamx) 
      dimension phin(npamx),ksep(npamx),tc(npamx)
*
      write(*,*) '--> splver2........'
*
      epsp = 1.d-10
      ngo1 = ngo
*--------------------------------------------
c    - definizione della ascissa curvilinea per le spline
      tp(1) = 0.d0
      xx(1) = xc(1)
      zz(1) = zc(1)
      do i = 1,npsl-1
        ax      = (xc(i+1)-xc(i))
        az      = (zc(i+1)-zc(i))
        dis     = sqrt(ax*ax+az*az)
        tp(i+1) = tp(i) + dis
        xx(i+1) = xc(i+1)
        zz(i+1) = zc(i+1)
      enddo

      yp1 = 1.d+31
      ypn = 1.d+31
      call spline(xx,zz,xs2,zs2,tp,yp1,ypn,npsl,npamx)
      do i = 2,npsl
        tts = (tp(i)+tp(i-1))/2.d0
        call splint(tts,xp(i),zp(i),xx,zz,xs2,zs2,tp,npsl,npamx,0)
      enddo
*--------------------------------------------
c   -  estrapolo i vertici a monte 
c      zz2 = zp(2)
c      xx2 = xp(2)
c      cms = (zz2-zc(1))/(xx2-xc(1))
c      cqs = zc(1) - cms*xc(1)
c      cm1 = tan(ande)
c      cq1 = proat
c      xp(1) = (cqs-cq1)/(cm1-cms)
c      zp(1) = cm1*xp(1) + cq1
*******************************-------------- xg zg ngo tg xgs2 zgs2
*
      x1  = xc(2)
      z1  = zc(2)
      x1  = xc(1)
      z1  = zc(1)
      w1x = xc(1)-xc(2)
      w1z = zc(1)-zc(2)
c      w1x = xc(1)-xp(2)
c      w1z = zc(1)-zp(2)
*
c -- test se cercare vertice tramite intersezione corpo-SL
      noint = 0
      if(nsep.eq.1)then 
         noint=1 
      endif
      write(*,*) 'splver, noint=',noint
        if(noint.eq.0)then
*
        ng0 = iint-10
        if(ng0.le.0) ng0 = 1
        iint=0
        do i = ng0,ngo-1
         w2x = xg(i+1) - xg(i)
         w2z = zg(i+1) - zg(i)
         x2  = xg(i)
         z2  = zg(i)
         call intret(xi,zi,r1,r2,x1,z1,w1x,w1z,x2,z2,w2x,w2z)
c         write(*,*) xg(i),zg(i),r2
         if(r2.ge.0.d0.and.r2.le.1.d0)then
           iint = i
           write(*,*) 'iint r2 ',iint,r2
* -- intersezione spline corpo con retta w1x,w1z:
           dt   = tg(i+1)-tg(i)
           eps0 = 0.01d0
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
           im = 0
           t0 = -eps0
           t1 = 1.d0+eps0
           write(*,*) 'tt1 .. ',tt1,tt2,tt3
           write(*,*) 't1 t0.. ',t1,t0
           if(tt1.ge.t0.and.tt1.le.t1) then
             ttt = tt1
             im  = im + 1
           endif
           if(tt2.ge.t0.and.tt2.le.t1) then
             im  = im + 1
             ttt = tt2
           endif
           if(tt3.ge.t0.and.tt3.le.t1) then
             im  = im + 1
             ttt = tt3
           endif
           if(im.gt.1) write(*,*) 'splver, intersezione doppia '
           if(im.eq.0) write(*,*) 'splver, intersezione non trovata'
           write(*,*) ' tt123 ',tt1,tt2,tt3
*
           c  = (ttt**3 - ttt)*(dt**2)/6.d0  
           d  = (-(ttt**3) + 3.d0*(ttt**2) - 2.d0*ttt)*(dt**2)/6.d0
           xip = ttt*xt1 + (1.d0-ttt)*xt2 + c*xd1 + d*xd2
           zip = ttt*zt1 + (1.d0-ttt)*zt2 + c*zd1 + d*zd2
*
*  -------------------------------------
             goto 999
         endif
        enddo
  999   continue     
        xi = xip
        zi = zip
        xp(1) = xi
        zp(1) = zi
        if(iint.eq.0)then
          write(*,*) 'SPLVER, intersezione non trovata ! ',jt
          write(*,*) 'SPLVER, prolungo corpo ! ',jt,ammax
          xp(1)  = x2 + w2x*r2 
          zp(1)  = z2 + w2z*r2
          iint  = ngo 
          ngo1  = ngo+1      
          xg(ngo1) = xp(1)
          zg(ngo1) = zp(1)
          dx = xg(ngo1)-xg(ngo1-1)
          dz = zg(ngo1)-zg(ngo1-1)
          dd = sqrt(dx**2+dz**2)
          tg(ngo1) = tg(ngo1-1)+dd
          tc(ngo1) = tg(ngo1)
        endif
*
* -- se noint = 1 estrapolo vertice e basta
        else
*
c          xp(1)    = xc(1) + 0.5d0*(xc(1)-xc(2))
c          zp(1)    = zc(1) + 0.5d0*(zc(1)-zc(2))
          ddd      = sqrt((xc(1)-xc(2))**2+(zc(1)-zc(2))**2)
          xp(1)    = xc(1) + di*(xc(1)-xc(2))/ddd
          zp(1)    = zc(1) + di*(zc(1)-zc(2))/ddd
          xp(1)    = xc(1) + 0.5d0*(xc(1)-xc(2))
          zp(1)    = zc(1) + 0.5d0*(zc(1)-zc(2))
c          dx = xc(1) - xcn(npc+ng)
c          dz = zc(1) - zcn(npc+ng)
c          dd = sqrt(dx**2 + dz**2)
c          ang1 = asin(dx/dd)
c          ang2 = ang-ang1
c          xp(1) = xcn(npc+ng)+di*cos(ang2)
c          zp(1) = zcn(npc+ng)+di*sin(ang2)
          write(*,*) 'PIZZO',ang,ang1,ang2,di
          write(*,*) xp(1),zp(1)
          write(*,*) xcn(npc+ng),zcn(npc+ng)
*
        endif     

* -- aggiungo punti separati alla spline corpo
c
          write(*,*) 'xp zp r2, ', xp(1),zp(1),r2
          do i=1,npc+ng
ccc            if(ksep(i).eq.2)then
ccc              dx  = xcn(i) - xg(ngo)
ccc              w1x = xg(ngo)- xg(ngo-1)
ccc              w1z = zg(ngo)- zg(ngo-1)
ccc              ngo1= ngo1+1
ccc              xg(ngo1) = xcn(i)
ccc             zg(ngo1) = zg(ngo)+w1z/w1x*dx
ccc              dx = xg(ngo1)-xg(ngo1-1)
ccc              dz = zg(ngo1)-zg(ngo1-1)
ccc              dd = sqrt(dx**2+dz**2)
ccc              tg(ngo1) = tg(ngo1-1)+dd
ccc            elseif(ksep(i).eq.1)then
            if(ksep(i).eq.1)then
              ngo1=ngo1+1
              xg(ngo1) = xcn(i)
              zg(ngo1) = zcn(i)
              dx = xg(ngo1)-xg(ngo1-1)
              dz = zg(ngo1)-zg(ngo1-1)
              dd = sqrt(dx**2+dz**2)
              tg(ngo1) = tg(ngo1-1)+dd
              tc(i)    = tg(ngo1)
            endif
          enddo 
          if(nsep.eq.1)then
c            ddd = sqrt((xc(1)-xc(2))**2+(zc(1)-zc(2))**2)
            ddd =sqrt((xg(ngo1)-xg(ngo1-1))**2+(zg(ngo1)-zg(ngo1-1))**2)
            xg(ngo1+1) = xg(ngo1) + di*(xg(ngo1)-xg(ngo1-1))/ddd 
            zg(ngo1+1) = zg(ngo1) + di*(zg(ngo1)-zg(ngo1-1))/ddd 
            xg(ngo1+1) = xg(ngo1) + 0.5d0*(xg(ngo1)-xg(ngo1-1)) 
            zg(ngo1+1) = zg(ngo1) + 0.5d0*(zg(ngo1)-zg(ngo1-1)) 
            ngo1 = ngo1+1      
            dx = xg(ngo1)-xg(ngo1-1)
            dz = zg(ngo1)-zg(ngo1-1)
            dd = sqrt(dx**2+dz**2)
            tg(ngo1) = tg(ngo1-1)+dd
            tc(ngo1) = tg(ngo1)
c          else
c            ngo1 = ngo1+1      
c            xg(ngo1) = xp(1)
c            zg(ngo1) = zp(1)
c            dx = xg(ngo1)-xg(ngo1-1)
c            dz = zg(ngo1)-zg(ngo1-1)
c            dd = sqrt(dx**2+dz**2)
c            tg(ngo1) = tg(ngo1-1)+dd
c            tc(ngo1) = tg(ngo1)
          endif
* -- ricalcolo coefficienti spline corpo
            yp1 = 1.d31
            ypn = 1.d31
            call spline(xg,zg,xgs2,zgs2,tg,yp1,ypn,ngo1,npamx)
        write(*,*) 'iint  xx ' ,iint,ngo1,ngo,noint
        write(*,'(2d25.10)')  xi,zi
        write(*,'(2d25.10)')  xip,zip
*

*************************************
c    - definizione della ascissa curvilinea per le spline
cckk*-------------------------------------------
cckk      tp(1) = 0.d0
cckk      xx(1) = xp(1)
cckk      zz(1) = zp(1)
cckk      npa = npsl+1      
cckk      do i = 1,npa-1
cckk        if (i.eq.1) then
cckk          ax      = (xc(i)-xp(i))
cckk          az      = (zc(i)-zp(i))
cckk        else
cckk          ax      = (xc(i)-xc(i-1))
cckk          az      = (zc(i)-zc(i-1))
cckk        end if
cckk        dis     = sqrt(ax*ax+az*az)
cckk        tp(i+1) = tp(i) + dis
cckk        xx(i+1) = xc(i)
cckk        zz(i+1) = zc(i)
cckk      enddo
cckk      yp1 = 1.d+31
cckk      ypn = 1.d+31
cckk      call spline(xx,zz,xs2,zs2,tp,yp1,ypn,npa,npamx)
cckk      do i = 2,npsl
cckk        tts = (tp(i+1)+tp(i))/2.d0
cckk        call splint(tts,xp(i),zp(i),xx,zz,xs2,zs2,tp,npa,npamx,0)
cckk      x1  = xc(2)
cckk      z1  = zc(2)
cckk      enddo
c--------------------------------------------------------------
c  - trovo nodo xp() xp() nel getto come media centroidi
c      do i=2,ng+1
      nn = 5 
c      nn = ng
c      nn = npsl
      do i=2,nn
c        ss  = float(i)/float(nn)
        ss  = 0.d0 
        xp2 = 0.5d0*(xc(i-1)+xc(i))
        zp2 = 0.5d0*(zc(i-1)+zc(i))
        xp(i) = xp2*(1.d0-ss)+ xp(i)*ss
      enddo

*
c       - verifica

c        if (zp(1).lt.proat) stop ' intersezione non fisica .. '

c        diff = zc(1) - (cm1*xc(1)+cq1)
c        if (diff.ge.0.d0) write(*,*) ' centroide dentro corpo ... '

*-------------------------------------------
c   -  estrapolo i vertici a valle 

      xp(npsl+1)  = xc(npsl) + (xc(npsl)-xp(npsl))
      zp(npsl+1)  = zc(npsl) + (zc(npsl)-zp(npsl))
c      xp(npsl+1)  = estr
c      ff = (estr-xc(npsl))/(xc(npsl)-xp(npsl)) 
c      zp(npsl+1)  = zc(npsl) + ff*(zc(npsl)-zp(npsl))
c      zp(npsl+1)  = 0.d0 
c      if(zp(npsl+1).lt.0.d0) zp(npsl+1) = 0.d0     
*
      return
      end
