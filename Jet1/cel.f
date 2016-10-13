
*---------------------------------------------------------*
*
c
      function cel(qqc,pp,cc,bb)
c
c     questa function calcola gli integrali ellittici 'completi'
c     del primo e secondo tipo; in particolare e(k) =
c     cel(kc,1,1,kc^2) e f(k) = cel(kc,1,1,1) con kc =sqrt(
c     1-k^2). la precisione e' il quadrato di ca.
c     (numerical recipes!!)
c
*      implicit double precision (a-h,o-z)
      include "slam_p.h"
      parameter (ca = 3.d-15)
c
      if(qqc.eq.0.d0) pause 'failure in cel'
      pig = acos(- 1.d0)
      pio = pig/2.d0
      qc = abs(qqc)
      a  = cc
      b  = bb
      p  = pp
      e  = qc
      em = 1.d0
      if(p.gt.0.d0) then
      p  = sqrt(p)
      b  = b/p
      else
      f  = qc*qc
      q  = 1.d0 - f
      g  = 1.d0 - p
      f  = f  - p
      q  = q*(b-a*p)
      p  = sqrt(f/g)
      a  = (a-b)/g
      b  = - q/(g*g*p) + a*p
      end if
   10 f  = a
      a  = a+b/p
      g  = e/p
      b  = b+f*g
      b  = b+b
      p  = g+p
      g  = em
      em = qc+em
      if(abs(g-qc).gt.g*ca) then
      qc = sqrt(e)
      qc = qc+qc
      e  = qc*em
      go to 10
      end if
c
      cel = pio*(b+a*em)/(em*(em+p))
      return
      end

 
