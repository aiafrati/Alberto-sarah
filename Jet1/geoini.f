      subroutine geoini(alpha,ztr,wl,xs00,dztr,al,bl,tl)
*      implicit double precision (a-h,o-z)
*      parameter(nmax=1000)
      include "slam_p.h"
      dimension xk(4*npamx),yk(4*npamx),zk(4*npamx)
      write(*,*) '.......... :  geoini'
*
      ca = cos(alpha)
      sa = sin(alpha)
      write(*,*) 'al    ',alpha,pi
      write(*,*) 'alpha ', alpha/pi*180.d0,ca,sa
      ximi =  1.0d+30
      xima = -1.0d+30
      etma = -1.0d+30
      zemi =  1.0d+30
      open(unit=7, file='geo.in')
      read(7,*) nsec
      do 10, i=1,nsec
        read(7,*) np
        do 20, j=1,np
           read(7,*) xi,et,ze
           ximi = min(xi,ximi)
           xima = max(xi,xima)
           etma = max(et,etma)
           zemi = min(ze,zemi) 
           if(j.eq.1)then
             xk(i) =  xi*ca + ze*sa
             yk(i) =  et
             zk(i) = -xi*sa + ze*ca
*             write(*,'(a3,4d15.5,i4)') 'ze ', ze,sa,ca,zk(i),i
           endif
  20    continue
  10  continue
      rewind(7)
      write(*,*) 'zk(nsec) ',nsec, zk(nsec)

*
      dztr = -(ztr - abs(zk(nsec)))
      al   = xima - ximi
      bl   = 2.d0*etma
      tl   = -zemi
*
      do 30, i = 1,nsec
         zk(i) = zk(i) + dztr
         if(zk(i).le.0.d0)then
           if(i.eq.1)then
             write(*,*) 'OKKIO, geoini : la scafo affonda!!!'
             stop
           endif
           w1x = 1.d0
           w1z = 0.d0
           x1  = 0.d0
           z1  = 0.d0
           w2x = xk(i) - xk(i-1)
           w2z = zk(i) - zk(i-1)
           x2  = xk(i-1)
           z2  = zk(i-1)
           call intret(xin,zin,r,s,x1,z1,w1x,w1z,x2,z2,w2x,w2z)
           xs00 = x2 + s*w2x
           wl   = xk(nsec) - xs00
           goto 99
         endif
  30  continue
  99  continue
*
      close(7)
      return
      end


