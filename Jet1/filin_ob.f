
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine filin( xc,zc,xp,zp,xs,zs,amp,phip,phis,fint3)

      include"slam_p.h"

      amx = xs - xp
      amz = zs - zp
      ds2 = amp/2.d0

      xme = (xs + xp)/2.d0
      zme = (zs + zp)/2.d0
      xic = ( (xc-xme)*amx + (zc-zme)*amz )/amp
      etc = ( (xc-xme)*amz - (zc-zme)*amx )/amp
      et2 = etc*etc
      tmm = (xic - ds2)
      tpp = (xic + ds2)
      tm2 = tmm*tmm
      tp2 = tpp*tpp
      elm = log(tm2+et2)
      elp = log(tp2+et2)

      if (et2.gt.0.d0) then
        fint2 = - (atan(tmm/etc)-atan(tpp/etc))
      else
        fint2 = 0.d0
      end if

      phime = (phip+phis)/2.d0
      dephi = (phis-phip)/amp
      fint3 = (phime+xic*dephi)*fint2 + dephi*etc*(elm-elp)/2.d0

      return
      end

