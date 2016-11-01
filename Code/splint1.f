
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine splint1(tts,xpo,xs,xs2,tp,n,nmax)

      include"slam_p.h"

      dimension xs(nmax),xs2(nmax),tp(nmax)

      klo=1
      khi=N
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


      return
      end
