
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine splint(tts,xpo,zpo,xs,zs,xs2,zs2,tp,n,nmax,kmon)

      include"slam_p.h"

      dimension xs(nmax),zs(nmax),xs2(nmax),zs2(nmax),tp(nmax)

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(tp(k).gt.tts)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=tp(khi)-tp(klo)
      if (h.eq.0.d0) stop 'bad xa input.'
      uh=1.d0/h
      a=(tp(khi)-tts)*uh
      b=(tts-tp(klo))*uh
      a1=(a**3-a)
      b1=(b**3-b)
      c1=(h**2)/6.d0

c     - coordinate punto interpolato

      xpo=a*xs(klo)+b*xs(khi)+
     &        (a1*xs2(klo)+b1*xs2(khi))*c1
      zpo=a*zs(klo)+b*zs(khi)+
     &        (a1*zs2(klo)+b1*zs2(khi))*c1

c     - verifica monotonicita' se richiesto

      if (kmon.ne.1.or.klo.eq.1.or.khi.eq.n) return

c     - verifica monotonicita' (solo per la componente zs)
c             .. normale al pannello precedente

      amx =  tp(klo) - tp(klo-1)
      amz =  zs(klo) - zs(klo-1)
      amt = sqrt(amx*amx+amz*amz)
      rnzp = -amx/amt
      rnxp =  amz/amt

c           .. normale al pannello i

      amx =  tp(khi) - tp(klo)
      amz =  zs(khi) - zs(klo)
      amt = sqrt(amx*amx+amz*amz)
      rnzl = -amx/amt
      rnxl =  amz/amt

c          .. normale al pannello i+1

      amx =  tp(khi+1) - tp(khi)
      amz =  zs(khi+1) - zs(khi)
      amt = sqrt(amx*amx+amz*amz)
      rnzs = -amx/amt
      rnxs =  amz/amt

c           .. calcolo distanze normali

      dins = (tp(khi)-tts)*rnxs+(zs(khi)-xpo)*rnzs
      dinl = (tp(khi)-tts)*rnxl+(zs(khi)-xpo)*rnzl
      dinp = (tp(klo)-tts)*rnxp+(zs(klo)-xpo)*rnzp
  
      if (dins.gt.0.d0.and.dinp.gt.0.d0.and.dinl.lt.0.d0) then 
        return
      else
        zpo = zs(klo)+
     &        (zs(khi)-zs(klo))/(tp(khi)-tp(klo))*(tts-tp(klo))
      end if

c     -- fine verifica



      return
      end
