
*---------------------------------------------------------------
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
       rn2x = tauy
       rn2y = -taux
       rn2sx = yss
       rn2sy = -xss
       conc = rn1x*rn2x + rn1y*rn2y
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
