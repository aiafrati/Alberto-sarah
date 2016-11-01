*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine prefin(npc,nng,npco,nngo,vfall,dt,dpht,dphi,tmy,tmz,
     #                  rny,rnz,vyo,vzo,yce,zce,yco,zco,phi,phio,
     #                  vym1,vzm1,pres)
*
      include"slam_p.h"
      dimension dpht(npamx),dphi(npamx),tmy(npamx),tmz(npamx)
      dimension rny(npamx),rnz(npamx),vyo(npamx),vzo(npamx)
      dimension yce(npamx),zce(npamx),yco(npamx),zco(npamx)  
      dimension phi(npamx),phio(npamx),pres(npamx) 
      dimension vym1(npamx),vzm1(npamx) 
      write(*,*)'--> prefin.... '
*
* calcolo  pressione con Dphi/Dt
*
      vf2 = 0.5d0*vfall**2
      vf2 = 1.d0
      do i=1,npc+nng
*
* -- calcolo velocita' fluido media tra vecchia e nuova
*
        vyn  = dpht(i)*tmy(i) + dphi(i)*rny(i)
        vzn  = dpht(i)*tmz(i) + dphi(i)*rnz(i)
c        vy   = 0.5d0*(vyn+vyo(i))
c        vz   = 0.5d0*(vzn+vzo(i))
        vy   = vyn
        vz   = vzn
        v2   = 0.5d0*(vy**2 + vz**2)/vf2
*
* -- velocita' pti corpo 
*
        
        wwy = (yce(i)-yco(i))/dt
        wwz = (zce(i)-zco(i))/dt
*
* -- calcolo pressione 
*
        xphi    = (phi(i)-phio(i))/dt
        wvv      = wwy*vy + wwz*vz   
        dpdtt    = (xphi - wvv)/vf2
        pres(i)  = (-dpdtt - v2)
*
*
        if (.not.(nng.eq.nngo.and.npc.eq.npco)) then
          pres(i)=0.d0
        endif
* 
      end do
*
      write(*,*) 'npco npc nngo nng',npco,npc,nngo,nng
      npco = npc
      nngo = nng
*
*   aggiorno
*
      do i=1,npc+nng
        yco(i) = yce(i)
        zco(i) = zce(i)
        phio(i)= phi(i)
        vyo(i) = vym1(i) 
        vzo(i) = vzm1(i) 
      enddo
*
      return
      end
