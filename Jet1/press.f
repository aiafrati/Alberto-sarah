*
* -----------------------------------------------------------
*
      subroutine press(iwig,dt,ux,vfall,phi,dphi,dpht,phio,
     #                     npc,npco,yc,zc,ty,tz,rny,rnz,
     #                     yco,zco,vyo,vzo,tto,amp,
     #                     pres )
      include "slam_p.h"
      dimension phi(npamx),dphi(npamx),dpht(npamx)
      dimension phio(npamx)
      dimension yc(npamx),zc(npamx),tt(npamx)
      dimension yco(npamx),zco(npamx),tto(npamx)
      dimension ty(npamx),tz(npamx),rny(npamx),rnz(npamx)
      dimension vyo(npamx),vzo(npamx)
      dimension amp(npamx),pres(npamx)
*
      write(*,*) '.......... :  press'
*
* - ascissa curvilinea
*
      if(iwig.eq.0)then
        vf2 = 0.5d0*vfall**2
      else
        vf2 = 0.5d0*ux**2
      endif 
      tt(1)=0.5d0*amp(1)
      do i=2,npc
        tt(i)=tt(i-1)+0.5d0*(amp(i-1)+amp(i))
      enddo
*
      do i=1,npc
*
* -- calcolo velocita' fluido media tra vecchia e nuova
*
        vyn  = dpht(i)*ty(i) + dphi(i)*rny(i)
	vzn  = dpht(i)*tz(i) + dphi(i)*rnz(i)
	vy   = 0.5d0*(vyn+vyo(i))
	vz   = 0.5d0*(vzn+vzo(i))
	v2   = 0.5d0*(vy**2 + vz**2)/vf2
*
* -- velocita' pti corpo 
*
        wwy = (yc(i)-yco(i))/dt
        wwz = (zc(i)-zco(i))/dt
*
* -- calcolo pressione 
*
	xphi    = (phi(i)-phio(i))/dt
        wvv      = wwy*vy + wwz*vz   
        dpdtt    = (xphi - wvv)/vf2
        pres(i)  = (-dpdtt - v2)
*
* --aggiorno
*
*        write(36,'(17d15.7)') yc(i),zc(i),dpht(i),dphi(i),v2,wwy,wwz,
*     #        xphi,wvv,phi(i),phio(i),dt,pres(i),vyo(i),vzo(i),vyn,vzn
*
* 
*        phio(i) = phi(i)
        tto(i)  = tt(i)
*        vyo(i)  = vyn
*        vzo(i)  = vzn
* 
* 
      end do
*      write(36,*)
*      write(36,*)
*      npco=npc
*
      return
      end
* ----------------------------------------------------------------------
