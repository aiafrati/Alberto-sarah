
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine calsolp(mb,mf,mt,m,n,nt,ntt,phi,phib,phb,dphisl,
     #                  dphnb,rl,a1,b1,c1,d1,e1,xi,ze,xis,zes,
     #                  ry,rz,kse,dphi)
*
      include"slam_p.h"
      dimension phi(npamx),phib(npamx),phb(npamx)
      dimension dphisl(npamx)
      dimension dphnb(npamx)
      dimension rl(npamx),ry(npamx),rz(npamx)
      dimension a1(npamx),b1(npamx),c1(npamx),d1(npamx),e1(npamx)
      dimension xi(npamx),ze(npamx),xis(npamx),zes(npamx)
      dimension kse(npamx),dphi(npamx)
      write(*,*) '--> calsol...........'
*
*      nngb = nget*kget
* - superficie corpo
      do i=1,mb
        if(kse(i).eq.0)then
        phi(i)  = phb(i)
        phib(i) = phb(i)
        else
        dphi(i) = dphnb(i)
        endif 
      enddo
      do i=mt+1,mt+m
        xi1   = xi(i)
        ze1   = ze(i)
        rli   = rl(i)
        ii    = mb+(i-mt)
        if(kse(i).eq.0)then
        phgi  = phg(a1,b1,c1,d1,e1,xi1,ze1,xis,zes,i)
        phi(ii)   = (1.d0-rli)*phb(i) + rli*phgi
        phib(ii)  = phb(i)
        else
        ryi   = ry(i)
        rzi   = rz(i)
        dphxii  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphzei  = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
        dphngi  = dphxii*ryi + dphzei*rzi
        dphi(ii)  = (1.d0-rli)*dphnb(i) + rli*dphngi
        endif 
      enddo
      do i = nt+1,nt+n
        xi1   = xi(i)
        ze1   = ze(i)
        rli   = rl(i)
        ii    = mb+m+(i-nt)
        if(kse(i).eq.0)then
        phgi  = phg(a1,b1,c1,d1,e1,xi1,ze1,xis,zes,i)
        phi(ii)   = phgi
        else
        ryi   = ry(i)
        rzi   = rz(i)
        dphxii  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphzei  = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
        dphngi  = dphxii*ryi + dphzei*rzi
        dphi(ii)  = dphngi
        endif 
      enddo 
*
* - sup. lib.
      do i = mb+1,mb+mf
        ii = m+n+(i-mb)
        dphisl(ii) = dphnb(i)
      enddo
      do i = mt+m+1,nt
        rli   = rl(i)
c        dhi   = hp(i)
c        sqh   = sqrt(1.d0+dhi**2)
        xi1   = xi(i)
        ze1   = ze(i)
        ryi   = ry(i)
        rzi   = rz(i)
        dphxii  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphzei  = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
c        dphngi  = (dphxii*dhi-dphzei)/sqh
c        dphtgi  = (-dphxii*-dhi*dphzei)/sqh
        dphngi  = dphxii*ryi + dphzei*rzi
        ii = n + (i-mt-m)
        dphisl(ii)  = (1.d0-rli)*dphnb(i) + rli*dphngi
      enddo
      do i = nt+n+1,ntt
        rli   = rl(i)
c        dhi   = hp(i)
c        sqh   = sqrt(1.d0+dhi**2)
        xi1   = xi(i)
        ze1   = ze(i)
        ryi   = ry(i)
        rzi   = rz(i)
        dphxii = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphzei = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
c        dphngi  = (dphxii*dhi-dphzei)/sqh
c        dphtgi  = (-dphxii*-dhi*dphzei)/sqh
        dphngi  = dphxii*ryi + dphzei*rzi
        ii =  (i-nt-n)
        dphisl(ii)  = dphngi
      enddo
*
      return 
      end
*
