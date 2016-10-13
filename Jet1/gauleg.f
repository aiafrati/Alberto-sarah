
*---------------------------------------------------------*
*
c
      subroutine gauleg(x1,x2,x,w,n)
*
* calcolo dei pesi e punti di integrazione di gauss
* x1,x2: intervallo di integrazione
*     n : numero di punti
*    x(n): vettore dei punti
*    w(n): vettore dei pesi
*
*      implicit double precision (a-h,o-z)
      include "slam_p.h"
      parameter (eps=3.d-15)
      dimension x(n),w(n)
*
      m =(n+1)/2
      pi=acos(-1.d0)
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(pi*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j - 1.d0)*z*p2-(j - 1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z - 1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end
*
