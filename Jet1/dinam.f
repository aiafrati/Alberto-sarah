
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dinam(proat,fot,devi,dems)

      include"slam_p.h"
      include"slam_v.h"

c     - calcolo termine di galleggiamento

cc      vbou = proat**2/tan(ande)

c     - calcolo accelerazione

cc      devi = - fot/weir + 9.81*(1.d0 - vbou/weir)
      devi  = - fot/weir - dzet*rkk/weir
      dems = dzet*rkk/(rga*weir)

      return
      end
*
