
*
      subroutine intret(xi,yi,s,t,x1,y1,w1x,w1y,x2,y2,w2x,w2y)
*----------------------------------------------------------------------
*  intersezione (xi,yi) tra le due rette (x1,y1) + s(w1x,w1y)
*                                      e (x2,y2) + t(w2x,w2y)
*  w1 e w2 non sono necessariamente di modulo unitario
*----------------------------------------------------------------------
      include "slam_p.h"
      parameter(eps=1.d-6)
*
* - scan e' il prod scal tra il versore di una e il versore normale
*   all'altra, ed e' nullo in caso di rette parallele
*
      a1 = sqrt(w1x**2+w1y**2) 
      a2 = sqrt(w2x**2+w2y**2)
      a = min(a1,a2)
      if (a.lt.eps**2) then
       write(*,*) 'intret, vettore nullo! '
       return
      endif  
      sca = w1y*w2x - w1x*w2y
      scan= sca/(sqrt(w1x**2+w1y**2)*sqrt(w2x**2+w2y**2))
      if (abs(scan).le.eps) then
        write(*,*) 'intret: okkio rette parallele!!!! '
        return
      endif
*
      s = ((w2y*x1 - w2x*y1) - (w2y*x2 - w2x*y2))/sca
*
      t = ((w1y*x1 - w1x*y1) - (w1y*x2 - w1x*y2))/sca
*
      xi = x1 + s*w1x
      yi = y1 + s*w1y
*
*      xi = x2 + t*w2x
*      yi = y2 + t*w2y
*
      return
      end

