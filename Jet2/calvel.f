
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,
     #                 phi,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phisl,
     #                 rl,a1,b1,c1,d1,e1,xi,ze,xis,zes,
     #                 amm,vy,vz,rny,rnz,tmy,tmz,
     #                 ammsl,vysl,vzsl,rnysl,rnzsl,tmysl,tmzsl,
     #                 vxi,ry,rz,ty,tz,kse)
*
      include"slam_p.h"
      dimension phi(npamx),phib(npamx),phb(npamx)
      dimension dpht(npamx),dphi(npamx),dphisl(npamx)
      dimension dphnb(npamx),dphtsl(npamx),phisl(npamx)
      dimension rl(npamx)
      dimension a1(npamx),b1(npamx),c1(npamx),d1(npamx),e1(npamx)
      dimension xi(npamx),ze(npamx),xis(npamx),zes(npamx)
      dimension amm(npamx),vy(npamx),vz(npamx),rny(npamx),rnz(npamx)
      dimension tmy(npamx),tmz(npamx)
      dimension ammsl(npamx),vysl(npamx),vzsl(npamx)
      dimension tmysl(npamx),tmzsl(npamx),rnysl(npamx),rnzsl(npamx) 
      dimension vxi(npamx)
      dimension ry(npamx),rz(npamx),ty(npamx),tz(npamx) 
      dimension kse(npamx)
      write(*,*) '--> calvel...........'
*
*      nngb = nget*kget
* - superficie corpo
c      do i=1,mb
c        phi(i)  = phb(i)
c        phib(i) = phb(i)
c      enddo
      do i=mt+1,mt+m
        xi1   = xi(i)
        ze1   = ze(i)
        rli   = rl(i)
        tyi   = ty(i)
        tzi   = tz(i)
        ii    = mb+(i-mt)
        if(i.eq.mt+m)then
          i2 = nt+1
        else
          i2 = i+1
        endif
        if(kse(i).ne.1)then
        vxixm   = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        vxixp   = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i2)
        vxix    = 0.5d0*(vxixm+vxixp)
        vxizm  = dphze(b1,d1,e1,xi1,ze1,xis,zes,i)
        vxizp  = dphze(b1,d1,e1,xi1,ze1,xis,zes,i2)
        vxiz    = 0.5d0*(vxizm+vxizp)
        vxi(ii)= tyi*vxix  + tzi*vxiz
        phib(ii)  = phb(i)
        else
        ryi   = ry(i)
        rzi   = rz(i)
        dphxiim  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphxiip  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i2)
        dphxii   = 0.5d0*(dphxiim+dphxiip)
        dphxiim = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphzeip = dphze(c1,d1,e1,xi1,ze1,xis,zes,i2)
        dphzei   = 0.5d0*(dphzeim+dphzeip)
        dphngi  = dphxii*ryi + dphzei*rzi
        dphtgi  = dphxii*tyi + dphzei*tzi
        dphi(ii)  = (1.d0-rli)*dphnb(i) + rli*dphngi
        endif
      enddo
      do i = nt+1,nt+n
        xi1   = xi(i)
        ze1   = ze(i)
        rli   = rl(i)
        tyi   = ty(i)
        tzi   = tz(i)
        ii    = mb+m+(i-nt)
        if(i.lt.nt+n)then
          i2 = i+1
        else
          i2 = i
        endif
        if(kse(i).ne.1)then
        vxixm = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        vxixp = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i2)
        vxix    = 0.5d0*(vxixm+vxixp)
        vxizm = dphze(b1,d1,e1,xi1,ze1,xis,zes,i)
        vxizp = dphze(b1,d1,e1,xi1,ze1,xis,zes,i2)
        vxiz    = 0.5d0*(vxizm+vxizp)
        vxi(ii)= tyi*vxix  + tzi*vxiz
        else
        ryi   = ry(i)
        rzi   = rz(i)
        dphxiim = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphxiip = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i2)
        dphxii   = 0.5d0*(dphxiim+dphxiip)
        dphzeim = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
        dphzeip = dphze(c1,d1,e1,xi1,ze1,xis,zes,i2)
        dphzei   = 0.5d0*(dphzeim+dphzeip)
        dphngi  = dphxii*ryi + dphzei*rzi
        dphtgi  = dphxii*tyi + dphzei*tzi
        dphi(ii)  = (1.d0-rli)*dphnb(i) + rli*dphngi
        endif
      enddo 
c      write(*,*)
c      write(*,*)
c
*
      if(kget.eq.1)then
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
        tyi   = ty(i)
        tzi   = tz(i)
        if(i.eq.mt+m+1)then
          i2 = i 
        else
          i2 = i-1
        endif
        dphxiip  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphxiim  = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i2)
        dphxii   = 0.5d0*(dphxiip + dphxiim)
        dphzeip  = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
        dphzeim  = dphze(c1,d1,e1,xi1,ze1,xis,zes,i2)
        dphzei   = 0.5d0*(dphzeip + dphzeim)
c        dphngi  = (dphxii*dhi-dphzei)/sqh
c        dphtgi  = (-dphxii*-dhi*dphzei)/sqh
        dphngi  = dphxii*ryi + dphzei*rzi
        dphtgi  = dphxii*tyi + dphzei*tzi
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
        tyi   = ty(i)
        tzi   = tz(i)
        if(i.eq.nt+n+1)then
          i2 = i
        else
          i2 = i-1
        endif
        dphxiip = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i)
        dphxiim = dphxi(b1,d1,e1,xi1,ze1,xis,zes,i2)
        dphxii   = 0.5d0*(dphxiip + dphxiim)
        dphzeip = dphze(c1,d1,e1,xi1,ze1,xis,zes,i)
        dphzeim = dphze(c1,d1,e1,xi1,ze1,xis,zes,i2)
        dphzei   = 0.5d0*(dphzeip + dphzeim)
c        dphngi  = (dphxii*dhi-dphzei)/sqh
c        dphtgi  = (-dphxii*-dhi*dphzei)/sqh
        dphngi  = dphxii*ryi + dphzei*rzi
        dphtgi  = dphxii*tyi + dphzei*tzi
        ii =  (i-nt-n)
        dphisl(ii)  = dphngi
      enddo
*
      endif
*
      nng=kget*ng
* -  dpht corpo
*
      do i=1,npc+nng
        if(i.eq.1)then
          amf      = 0.5d0*(amm(i)+amm(i+1))
          dphtff   = (phi(i+1)-phi(i))/amf 
          dpht(i)  = dphtff
        elseif(i.eq.npc+nng)then
          amb      = 0.5d0*(amm(i)+amm(i-1))
          dphtbb   = (phi(i)-phi(i-1))/amb 
          dpht(i)  = dphtbb 
        else
          amf      = 0.5d0*(amm(i)+amm(i+1))
          amb      = 0.5d0*(amm(i)+amm(i-1))
          dphtff   = (phi(i+1)-phi(i))/amf 
          dphtbb   = (phi(i)-phi(i-1))/amb 
          dpht(i)  = 0.5d0*(dphtff+dphtbb)
        endif
        vy(i) = dpht(i)*tmy(i)+ dphi(i)*rny(i)
        vz(i) = dpht(i)*tmz(i)+ dphi(i)*rnz(i)
      enddo
* - dphtsl sup lib
      do i=1,npsl
        if(i.eq.1)then
          amf       = 0.5d0*(ammsl(i)+ammsl(i+1))
          dphtff    = (phisl(i+1)-phisl(i))/amf 
          dphtsl(i) = dphtff
        elseif(i.eq.npsl)then
          amb        = 0.5d0*(ammsl(i)+ammsl(i-1))
          dphtbb     = (phisl(i)-phisl(i-1))/amb 
          dphtsl(i)  = dphtbb
        else
          amb       = 0.5d0*(ammsl(i)+ammsl(i-1))
          amf       = 0.5d0*(ammsl(i)+ammsl(i+1))
          dphtbb    = (phisl(i)-phisl(i-1))/amb 
          dphtff    = (phisl(i+1)-phisl(i))/amf 
          dphtsl(i) = 0.5d0*(dphtff+dphtbb)
        endif
        vysl(i) = dphtsl(i)*tmysl(i) + dphisl(i)*rnysl(i)
        vzsl(i) = dphtsl(i)*tmzsl(i) + dphisl(i)*rnzsl(i)
      enddo
*
      return
      end
*
