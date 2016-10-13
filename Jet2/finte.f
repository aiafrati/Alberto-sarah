
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine finte( xc,zc,xp,zp,xs,zs,amp,fint1,fint2,eppa)

      include"slam_p.h"
      eps = 1.d-20

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

      if (et2.gt.eps) then
        fint1 = -(tmm*elm-tpp*elp)/2.d0 +
     &           (tmm-tpp)-etc*(atan(tmm/etc)-atan(tpp/etc))
        fint2 = -(atan(tmm/etc)-atan(tpp/etc))
      else 
        fint1 = -(tmm*elm-tpp*elp)/2.d0 + (tmm-tpp)
        fint2 = 0.d0
      end if

      return
      end
