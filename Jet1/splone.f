
*
      subroutine splone(y2,t,y,yp1,ypn,n,nmax)
      include "slam_p.h"
*
*    questa subroutine calcola i coefficienti y2 per l'interpolazi one
*    spline della curva yi(ti),  i=1,n
*
*
      dimension t(nmax),y(nmax),y2(nmax),u(npamx)
*
*     - interpolazione componente y 
*
      if (yp1.gt.1.d+30) then
        y2(1) = 0.d0
        u(1)  = 0.d0
      else
        y2(1) = -0.5d0
        u(1)  = (3.d0/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1)) - yp1)
      end if
      do 11 i = 2,n-1
        sig     = (t(i)-t(i-1))/(t(i+1)-t(i-1))
        p       = sig*y2(i-1)+2.d0
        y2(i)   = (sig-1.d0)/p
        u(i)    = (6.d0*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1))
     &            /(t(i)-t(i-1)))/(t(i+1)-t(i-1)) - sig*u(i-1))/p
   11 continue
      if (ypn.gt.1.d+30) then
        qn = 0.d0
        un = 0.d0
      else
        qn = 0.5d0
        un = (3.d0/(t(n) - t(n-1)))*(ypn - (y(n)-y(n-1))/(t(n)-t(n-1)))
      end if
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k = n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
   12 continue
*
      return
      end
*
