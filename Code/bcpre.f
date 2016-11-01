      subroutine bcpre(ng,kget,vfall,npc,npsl,tmz,amp,phi,dphi,dpht,
     #                 vxi,dphisl,dphtsl,yce,zce,dpnt,dpntt,dptsl,jt,tn,
     #                 ygn,zgn,ygs2,zgs2,tg,ngo1,tc)
      include "slam_p.h"
      dimension tmz(npamx),phi(npamx),dphi(npamx),dpht(npamx)
      dimension dphisl(npamx),dphtsl(npamx)
      dimension yce(npamx),zce(npamx),amp(npamx),vxi(npamx)  
      dimension dpnt(npamx),dptsl(npamx),dpht2(npamx),dpht2b(npamx) 
      dimension dpntt(npamx),tn(npamx),tg(npamx),tc(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx) 
*
      write(*,*) '--> bcpre......'
* -- calcolo derivata tangenziale della velocita tangenziale
* - primo punto , interpolaz y = a + b t + c t**2 , t asciss curv
c      t0 = 0.d0
c      y0 = -vfall*tmz(1)
c      t1 = amp(1)
c      y1 = (phi(2)-phi(1))/(0.5d0*(amp(1)+amp(2)))
c      t2 = amp(1)+amp(2)
c      y2 = (phi(3)-phi(2))/(0.5d0*(amp(2)+amp(3)))
c      call der2gr(y0,y1,y2,t0,t1,t2,a,b,c)
* -------  calcolo derivata analiticamente in t = t1/2
c      dpht2(1) = b + 2.d0*(t1/2.d0)*c
      nng      = ng*kget
      dphtu    = -vfall*tmz(1)
      dphtd    = (phi(2)-phi(1))/(0.5d0*(amp(1)+amp(2)))
      dpht2(1)  = ( dphtd - dphtu )/amp(1)
c      write(*,*) 'bcpre2 ' , dpht2(1),dpht222
      i = 1 
      do i = 2,npc+nng-1
        amu = 0.5d0*(amp(i-1)+amp(i))
        amd = 0.5d0*(amp(i+1)+amp(i))
        dpht2u    = (dpht(i)-dpht(i-1))/amu
        dpht2d    = (dpht(i+1)-dpht(i))/amd
        dpht2(i)  = 0.5d0*(dpht2u+dpht2d)
        vxiu      = (vxi(i)-vxi(i-1))/amu
        vxid      = (vxi(i+1)-vxi(i))/amd
        dpht2b(i) = 0.5d0*(vxiu+vxid)
      enddo     
      dpht2(npc+nng) = (dpht(npc+nng)-dpht(npc+nng-1))/
     #                 (0.5d0*(amp(npc+nng)+amp(npc+nng-1)))
      dpht2b(npc+nng)= (vxi(npc+nng)-vxi(npc+nng-1))/
     #                 (0.5d0*(amp(npc+nng)+amp(npc+nng-1)))
* - calcolo ws d(un)/ds = ws w.dn/ds e gli altri termini
      wy  = 0.d0
      wz  = -vfall
ccc - cilindro 
c      rho = 0.2d0
c      cur = 1.d0/rho
c      cur = 0.d0
ccc
 
ccc      write(88,*) '#' , jt
ccc     write(47,*) '#' , jt
ccc    write(46,*) '#' , jt
      do i=1,npc+nng
        ws     = wz*tmz(i)
        wn     = dphi(i)
c         tt = 0.5d0*(tn(i)+tn(i+1))
         tt = tc(i) 
        yy     = yce(i)
        zz     = zce(i)
ccc        write(46,'(i4,3d15.7)') i,yy,zz,tt
        call splont(tt,yy,tg,ygn,ygs2,ngo1,npamx,yt,1,ytt,1,yttt,1)
        call splont(tt,zz,tg,zgn,zgs2,ngo1,npamx,zt,1,ztt,1,zttt,1)
        call curv(yy,yt,ytt,yttt,zz,zt,ztt,zttt,cur,dsrny,dsrnz,conc)
ccc        write(47,'(i4,6d15.7)') i,yy,zz,cur,dsrny,dsrnz,conc
ccc - cilindro 
c        thet   = asin(yy/rho)
c         thet  = 0.d0
c        dsrny  = cur*cos(thet)
c        dsrnz  = cur*sin(thet)
ccc 
        wsdun  = ws*(wy*dsrny + wz*dsrnz)
        wsus   = ws*dpht(i)*cur
        wnwn   = wn**2*cur
* wsus e wnwn cambiano  segno con la concavita (+ per cilindro)
        if(conc.le.0.d0)then
          dpnt(i) = wsus + wnwn - wsdun + wn*dpht2(i)
          dpntt(i)= wsus + wnwn - wsdun + wn*dpht2(i)
        else
          dpnt(i) = -wsus - wnwn - wsdun + wn*dpht2(i)
          dpntt(i)= -wsus - wnwn - wsdun + wn*dpht2(i)
        endif
ccc        write(88,'(6d15.7)') yce(i),wsus,wnwn,-wsdun,wn*dpht2(i)
      enddo
ccc      write(46,*)
ccc      write(46,*)
ccc      write(47,*)
ccc      write(47,*)
ccc      write(88,*)
ccc      write(88,*)
* -- calcolo dphi/dt sulla SL
      do i = 1,npsl
        dptsl(i) = -(dphisl(i)**2 + dphtsl(i)**2)/2.d0
      enddo
*
      return
      end
*


*---------------------------------------------------------------
*
      subroutine der2gr(y0,y1,y2,t0,t1,t2,a,b,c)
      include "slam_p.h"
c      implicit double precision (a-h,o-z)        
*
      eps  = 1.d-12
*
      dt10 = t1 - t0
      dt20 = t2 - t0
      dt21 = t2 - t1
      if(dt10.le.eps) stop 't1 - t0 troppo piccolo'
      if(dt20.le.eps) stop 't2 - t0 troppo piccolo'
      if(dt21.le.eps) stop 't2 - t1 troppo piccolo'
*
      a = y0*t1*t2/(dt20*dt10) - y1*t0*t2/(dt21*dt10)
     #   +y2*t0*t1/(dt20*dt21)
*     
      b = -y0*(t2+t1)/(dt20*dt10) + y1*(t2+t0)/(dt21*dt10)
     #    -y2*(t1+t0)/(dt20*dt21)
*
      c = y0/(dt20*dt10) - y1/(dt21*dt10) + y2/(dt20*dt21)
*
      return
      end 
* -------------------------------------------------------
*
       subroutine curv(x,xt,xtt,xttt,y,yt,ytt,yttt,cur,rn2sx,rn2sy,
     #                 conc)
       include "slam_p.h"
*
*  curv calcola le derivate rispetto all'ascissa curvilinea, la
*  curvatura e vettori tangente e normale, e la derivata tangente
*  del vettore normale. tutto cio' a pertire dai vori di
*  x,xt,xtt,xttt, y,yt,ytt,yttt
*  cioe' il pto e le sue derivate fino alla terza, fatte
*  rispetto al generico paramtro t che descrive la curva x(t), y(t)
*
*
       st  = sqrt(xt**2 + yt**2) 
       stt = (xt*xtt + yt*ytt)/st
*
       xs  = xt/st
       ys  = yt/st
* -cur = curvatura
       xss = (xtt - xs*stt)/st**2
       yss = (ytt - ys*stt)/st**2
       cur = sqrt(xss**2 + yss**2)
*
       sttt = -stt**2/st + (xtt**2 + xt*xttt + ytt**2 + yt*yttt)/st
       xsss = (xttt - 3.*st*stt*xss - sttt*xs)/st**3 
       ysss = (yttt - 3.*st*stt*yss - sttt*ys)/st**3 
* -tau = dx/ds, versore tangente 
       taux = xs
       tauy = ys
       tau  = sqrt(xs**2 + ys**2) 
*
* -rn1 = normale definita come dtau/ds; cambia segno con la concavita',
*        e non e' definita nei pti di flesso (dove xss=yss=0, curv=0!!).
*        inoltre non e' definita nel primo e ultimo pto se si e'
*        interpolato fissando i coefficienti y2(1)=y2(n)=0 in spline,
*        ovvero ponendo y1n=ypn=1.e+31. infatti ytt = y2 in ogni pto ti!
* -rn2 = normale definita come R.tau, R rotazione oraria di pi/2. non
*        cambia segno con la concavita', ed e' definita ovunque.
*
       conc = 0.d0
       rn1x = 0.d0
       rn1y = 0.d0
       if(cur.gt.0)then
       rn1x = xss/cur
       rn1y = yss/cur
*
* -rn1s = drn1/ds   , rn2s=drn2/ds 
*       
       rn1sx = xsss/cur - (xss*xsss + yss*ysss)*xss/cur**3
       rn1sy = ysss/cur - (xss*xsss + yss*ysss)*yss/cur**3
       endif
       rn2x  = tauy
       rn2y  = -taux
       rn2sx = yss
       rn2sy = -xss
       conc  = rn1x*rn2x + rn1y*rn2y
*
*       write(*,'(2e14.6)')  rn1sx , rn1sy 
*       write(*,*)
*
*
*
*
       return
       end
*



* ----------------------------------------------------------------------
*
      subroutine splont(tpo,ypo,t,y,y2,n,nmax,
     #                  yt,kyt,ytt,kytt,yttt,kyttt)
      include "slam_p.h"
*
*  splint calcola il valore ypo del punto interpolato spline, che
*  e' individuato sulla curva yi(ti), i =1,n dal valore tpo del
*  parametro t che descrive la curva. 
*
      dimension y(nmax),y2(nmax),t(nmax)
*
      klo = 1
      khi = n
 1    if (khi-klo.gt.1) then
        k = (khi+klo)/2
        if(t(k).gt.tpo)then
          khi = k
        else
          klo = k
        endif
      goto 1
      endif
      h  = t(khi)-t(klo)
      if (h.eq.0.d0) stop 'bad xa input.'
*      uh=1.d0/h
      a  = (t(khi)-tpo)/h
      b  = (tpo-t(klo))/h
      a1 = (a**3-a)
      b1 = (b**3-b)
      c1 = (h**2)/6.d0
*
*     - coordinate punto interpolato
*
      ypo = a*y(klo)+b*y(khi) + (a1*y2(klo)+b1*y2(khi))*c1
*      
      if(kyt+kytt+kyttt.eq.0)then
        return
      endif
      if(kyt.eq.1)then
*   -   derivata prima   (rispetto a t)
        coy = (y(khi) - y(klo))/h
        yt  = coy - (3.d0*a**2-1.d0)/6.d0*h*y2(klo) 
     #            + (3.d0*b**2-1.d0)/6.d0*h*y2(khi)
      endif
      if(kytt.eq.1)then
*   -   derivata seconda (rispetto a t)
        ytt = a*y2(klo) + b*y2(khi)
      endif
      if(kyttt.eq.1)then
*   -   derivata terza   (rispetto a t)
        yttt = (y2(khi)-y2(klo))/h
      endif
*
      return
      end
*
*
*---------------------------------------------------
*  kmon va inserito tra gli argomenti della subroutine....
*
c     - verifica monotonicita' se richiesto
*
*      if (kmon.ne.1.or.klo.eq.1.or.khi.eq.n) return
*
c     - verifica monotonicita' (solo per la componente zs)
c             .. normale al pannello precedente
*   qui zs sta al posto di y, tp al posto di t, tts al posto di tpo
*   zpo al posto di ypo
*
*      amx =  tp(klo) - tp(klo-1)
*      amz =  zs(klo) - zs(klo-1)
*      amt = sqrt(amx*amx+amz*amz)
*      rnzp = -amx/amt
*      rnxp =  amz/amt
*
c           .. normale al pannello i
*
*      amx =  tp(khi) - tp(klo)
*      amz =  zs(khi) - zs(klo)
*      amt = sqrt(amx*amx+amz*amz)
*      rnzl = -amx/amt
*      rnxl =  amz/amt
*
c          .. normale al pannello i+1
*
*      amx =  tp(khi+1) - tp(khi)
*      amz =  zs(khi+1) - zs(khi)
*      amt = sqrt(amx*amx+amz*amz)
*      rnzs = -amx/amt
*      rnxs =  amz/amt
*
c           .. calcolo distanze normali
*
*      dins = (tp(khi)-tts)*rnxs+(zs(khi)-xpo)*rnzs
*      dinl = (tp(khi)-tts)*rnxl+(zs(khi)-xpo)*rnzl
*      dinp = (tp(klo)-tts)*rnxp+(zs(klo)-xpo)*rnzp
*  
*      if (dins.gt.0.d0.and.dinp.gt.0.d0.and.dinl.lt.0.d0) then 
*        return
*      else
*        zpo = zs(klo)+
*     &        (zs(khi)-zs(klo))/(tp(khi)-tp(klo))*(tts-tp(klo))
*      end if
*
c     -- fine verifica
*
*
*
*      return
*      end


