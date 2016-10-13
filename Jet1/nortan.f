

      subroutine nortan(nstart,nend,yv,zv,aamy,aamz,amp,
     #                  ttmy,ttmz,rrny,rrnz,yc,zc,icen)
*
      include "slam_p.h"
      dimension yv(npamx+1),zv(npamx+1),aamy(npamx),aamz(npamx)
      dimension amp(npamx),ttmy(npamx),ttmz(npamx)
      dimension rrny(npamx),rrnz(npamx),yc(npamx),zc(npamx)
*
      write(*,*) '.......... :  nortan'
*
      do ip = nstart,nend
        amy =  yv(ip+1) - yv(ip) 
        amz =  zv(ip+1) - zv(ip) 
        am  =  sqrt(amy**2+amz**2)
        tmy =  amy/am
        tmz =  amz/am
*  
        amp(ip)  =  am
        aamy(ip) =  amy 
        aamz(ip) =  amz 
        ttmy(ip) =  amy/am
        ttmz(ip) =  amz/am
        rrny(ip) =  tmz
        rrnz(ip) = -tmy
      end do
*
      if (icen.gt.0) then
      do ip = nstart,nend
        yc(ip)   =  0.5d0*(yv(ip+1)+yv(ip))
        zc(ip)   =  0.5d0*(zv(ip+1)+zv(ip))
      end do
      endif
*
      return
      end
*
