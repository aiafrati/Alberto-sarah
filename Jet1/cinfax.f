

*---------------------------------------------------------*
*
      subroutine cinfax(eps2,xc,zc,x1,z1,x2,z2,t,w,n,cinfg,cinfdg)
*      implicit double precision (a-h,o-z)
      include "slam_p.h"
*
* calcolo dei coeff. di influenza per problema a simm assiale
*
*      parameter(nmax=20)
      dimension t(ngmax),w(ngmax)
*      pi=4.d0*atan(1.d0)
      epsh2=0.d0
* - tx e tz non sono normalizzati
      tx = x2-x1
      tz = z2-z1
      amp = sqrt(tx**2 + tz**2)
      rnx = tz/amp
      rnz = -tx/amp
      cinfg  = 0.d0
      cinfdg = 0.d0
c      write(*,*) 'cionf '
      do i = 1,n
        ww = w(i)*amp
        xq = x1 + t(i)*tx
        zq = z1 + t(i)*tz
        d2 = (xc-xq)**2 + (zc-zq)**2 + eps2
        f  = 4.d0*xc*xq
        h2 = f/(d2+f)
        rk = sqrt(1.d0 - h2)
        rk2 = 1.d0 - h2
        if(h2.lt.epsh2)then
*          write(*,*) 'h2 < 10-2 !!! '
          g = - pi*2.d0*xq/(d2+f)**0.5d0
         dg = - pi*2.d0*xq*(rnz*(zc-zq)-rnx*xq)/(d2+f)**1.5d0
        else
          eh2 = cel(rk,1.d0,1.d0,rk2)/rk2
          fh  = cel(rk,1.d0,1.d0,1.d0)
          gh  = eh2*(2.d0/h2 - 1.d0) - 2.d0/h2*fh
          g =-4.d0*xq/(d2+f)**0.5d0*fh
        dg=4.d0*xq/(d2+f)**1.5d0*(rnx*(xc*gh - xq*eh2)+rnz*(zc-zq)*eh2)
        endif
        cinfg  = cinfg  + ww*g
        cinfdg = cinfdg + ww*dg
*        write(*,'(5d15.6)') xq,x2,xc,gh,h2
      end do
* 
      return
      end      
*
