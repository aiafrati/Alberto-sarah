
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

