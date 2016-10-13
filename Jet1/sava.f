* ------------------------------------------------------
      subroutine savan(al,bl,alpha,dztr,zmax,xx,yg,zg,rnxg,ng)
*      implicit double precision (a-h,o-z)
      include "slam_p.h"
***********************************************************************
* - dati  al (lunghezza) e yb (larghezza), e dati gli angoli          *
*   di deadrise agli estremi, calcola la carena, in base              *
*   all'equazione  z(x,y) = a(x)*y + b(x)*y**2, e la normale, prodotto*
*   vettore (normalizzato a 1) dei due vettori tangenti               *
*   (0,1, a + 2by)  e (1,0, (dxa)*y + (dxb)*y**2)  , che si           *
*   ottengono differenziando rispett. in y 2 x la                     *
*   rappresentazione parametrica della superficie di carena           *
*   (x,y,z(x,y))                                                      *
*                                                                     *
*   ATTENZIONE : xx yy zz coordinate nel sistema 'assoluto',          *
*                x   y  z coordinate nel sistema della nave           * 
***********************************************************************
      dimension yg(npamx),zg(npamx),rnxg(npamx)
* 
      pi = acos(-1.d0)
      epsp = 1.d-8
*
      yb = 0.5d0*bl
* - angoli agli estremi : 0 = prua, 1 = poppa, a = centro, b = murata
      theta0 =  14.d0/180.d0*pi 
      thetb0 =  14.d0/180.d0*pi 
      theta1 =  14.d0/180.d0*pi
      thetb1 =  14.d0/180.d0*pi 
*
      ca = cos(alpha)
      sa = sin(alpha)
*
      z  = 0.d0
      dz = zmax/float(ng) 
*
      write(*,*) 'ng ',ng
      do 10, i=1,ng+1 
*
*      z = 0.d0
      x  = (xx - z*sa)/ca
      zz = -x*sa + z*ca + dztr
*
      a = zpa(al,theta0,theta1,x)
      b =(zpb(al,thetb0,thetb1,x) - zpa(al,theta0,theta1,x))/(2.d0*yb)
      delta = a**2 + 4.d0*b*z
      if (delta.lt.0.d0) STOP 'Savander delta < 0!'
*
      if(abs(b).le.epsp)then
        y = z/a
      else
        y  = (-a + sqrt(a**2 + 4.d0*b*z))/(2.d0*b)
      endif
      yy = y
      zt = a*y + b*y**2
*      zz = z
*      write(*,'(6d13.6)') a,b,a**2 + 4.d0*b*z,yy,z,zt
      dxa = dxzpa(al,theta0,theta1,x)
      dxb = ( dxzpb(al,thetb0,thetb1,x) - dxa ) / (2.d0*yb)
      rnx = dxa*y + dxb*y**2
*      rnx = dxzpa(al,theta0,theta1,x)*y+dxzpb(al,thetb0,thetb1,x)*y**2
      rny = a + 2.d0*b*y
      rnz = -1.d0
*
      rn  = sqrt(rnx**2 + rny**2 + rnz**2)
      rnx = (rnx*ca + rnz*sa)/rn
      rny = rny/rn
      rnz = (-rnx*sa + rnz*ca)/rn
*
      yg(i)   = yy
      zg(i)   = zz
      rnxg(i) = rnx
*      write(*,*) yg(i),zg(i),rnxg(i)
*
      z = z + dz
  10  continue
*
      return
      end
*--------------------------------------------------------
      function zpa(al,theta0,theta1,x)
      implicit double precision (a-h,o-z)
*
      pi = acos(-1.d0)
*
      zpa0 = tan( theta0)
      zpa1 = tan( theta1)
* -- lineare in zpa da prua a poppa
      zpa  = zpa0*(1.d0 - x/al) + zpa1*x/al
* -- parabolica in zpa (con derivata nulla a poppa)
*      zpa  = zpa1 + (zpa0 - zpa1)*(x/al - 1.d0)**2
* -- lineare in thet
      thet = theta0*(1.d0 - x/al) + theta1*x/al
*      zpa  = tan(thet)
* -- parabolica in thet
      thet = theta1 + (theta0 - theta1)*(x/al - 1.d0)**2
      zpa  = tan(thet)
*
      return
      end
* --------------------------------------------------------
      function zpb(al,thetb0,thetb1,x)
      implicit double precision (a-h,o-z)
*
      pi = acos(-1.d0)
*
      zpb0 = tan(thetb0)
      zpb1 = tan(thetb1)
* -- lineare da prua a poppa
      zpb  = zpb0*(1.d0 - x/al) + zpb1*x/al
* -- parabolica (con derivata nulla a poppa)
*      zpb  = zpb1 + (zpb0 - zpb1)*(x/al - 1.d0)**2
* -- lineare in thet
      thet = thetb0*(1.d0 - x/al) + thetb1*x/al
*      zpb  = tan(thet)
* -- parabolica in thet
      thet = thetb1 + (thetb0 - thetb1)*(x/al - 1.d0)**2
      zpb  = tan(thet)
*
      return
      end
* ---------------------------------------------------------
      function dxzpa(al,theta0,theta1,x)
      implicit double precision (a-h,o-z)
*
      pi = acos(-1.d0)
*
      zpa0 = tan(theta0)
      zpa1 = tan(theta1)
* -- derivata del caso lineare
      dxzpa = (zpa1 - zpa0)/al
* -- derivata del caso parabolico
*      dxzpa = (zpa0 - zpa1)*2.d0/al*(x/al - 1.d0) 
* -- derivata caso lineare in thet
      thet = theta0*(1.d0 - x/al) + theta1*x/al
*      dxzpa = 1.d0/cos(thet)**2 * (theta1 - theta0)/al
* -- derivata caso parabolico in thet
      thet = theta1 + (theta0 - theta1)*(x/al - 1.d0)**2
      dxzpa = 1.d0/cos(thet)**2*(theta0-theta1)*2.d0/al*(x/al-1.d0)
      return
      end
* ---------------------------------------------------------
      function dxzpb(al,thetb0,thetb1,x)
      implicit double precision (a-h,o-z)
*
      pi = acos(-1.d0)
*
      zpb0 = tan(thetb0)
      zpb1 = tan(thetb1)
* -- derivata del caso lineare
      dxzpb = (zpb1 - zpb0)/al
* -- derivata del caso parabolico
*      dxzpb = (zpb0 - zpb1)*2.d0/al*(x/al - 1.d0) 
* -- derivata caso lineare in thet
      thet = thetb0*(1.d0 - x/al) + thetb1*x/al
*      dxzpb = 1.d0/cos(thet)**2 * (thetb1 - thetb0)/al
* -- derivata caso parabolico in thet
      thet = thetb1 + (thetb0 - thetb1)*(x/al - 1.d0)**2
      dxzpb = 1.d0/cos(thet)**2*(thetb0-thetb1)*2.d0/al*(x/al-1.d0)
*
      return
      end
* 

