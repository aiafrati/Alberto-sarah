
*
      function rintli(ss,s,y,n)
      include "slam_p.h"
      dimension y(npamx),s(npamx)
*
      klo = 1
      khi = n
 1    if (khi-klo.gt.1) then
        k = (khi+klo)/2
        if(s(k).gt.ss)then
          khi = k
        else
          klo = k
        endif
      goto 1
      endif
c      if (s(khi)-s(klo).eq.0.d0) STOP 'bad  input'
      if (s(khi).lt.ss) then
c         yy = y(khi)
	 c  = (y(khi)-y(khi-1))/(s(khi)-s(khi-1))
	 ds = ss-s(khi)
         yy = y(khi)+c*ds
      elseif(s(klo).gt.ss)then
c         yy = y(klo)
	 c  = (y(klo+1)-y(klo))/(s(klo+1)-s(klo))
	 ds = s(klo)-ss
         yy = y(klo)+c*ds
      else
         c  = (y(khi)-y(klo))/(s(khi)-s(klo))
	 ds = ss-s(klo)
         yy = y(klo)+c*ds
      endif
*      
      rintli=yy     
*
      return
      end
*---------------------------------------------------
      function infind(ss,s,n)
      include "slam_p.h"
      dimension s(npamx)
*
*
      klo = 1
      khi = n
      if (s(khi).lt.ss) then
        infind=khi
	return
      elseif(s(klo).gt.ss)then
        infind=klo
        return
      else
*      
 1    if (khi-klo.gt.1) then
        k = (khi+klo)/2
        if(s(k).gt.ss)then
          khi = k
        else
          klo = k
        endif
      goto 1
      endif
*
      infind = klo
      endif
*
      return
      end
      subroutine prein(yv,zv,npc,yc,zc,tto,yvo,zvo,yvo2,zvo2,so,npco,
     #           pho,dpho,dpto,pho2,dpho2,dpto2)
      include "slam_p.h"
*
* -- mi servono: tto,yvo,zvo,yvo2,zvo2,so,npco
*                pho,dpho,dpto,pho2,dpho2,dpto2
*
      dimension phi(npamx),dphi(npamx),dpht(npamx)
      dimension pho(npamx),dpho(npamx),dpto(npamx)
      dimension pho2(npamx),dpho2(npamx),dpto2(npamx)
      dimension yc(npamx),zc(npamx),yv(npamx),zv(npamx)
      dimension tto(npamx),yvo(npamx),zvo(npamx),yvo2(npamx),zvo2(npamx)
      dimension so(npamx)
      
*     


* - ascissa curvilinea
*
      so(1) = 0.d0
      so(2) = sqrt( (yc(1)-yv(1))**2 + (zc(1)-zv(1))**2 )
      npco=npc
      do i=2,npco
        dy  = yc(i)-yc(i-1)
        dz  = zc(i)-zc(i-1)
        so(i+1)= so(i) + sqrt(dy**2 + dz**2)
      end do
      dy = yv(npco+1)-yc(npco)
      dz = zv(npco+1)-zc(npco)
      so(npco+2)=so(npco+1) + sqrt( dy**2 + dz**2 )
*
*
* - estrapolazione lineare a monte e a valle
*
      dy = yv(npc+1)-yc(npc)
      dz = zv(npc+1)-zc(npc)
      ds = sqrt(dy**2 + dz**2)
      c     = (phi(2)-phi(1))/(so(3)-so(2))
      phi0  = phi(1) - c*so(2)
      c     = (phi(npc)-phi(npc-1))/(so(npc+1)-so(npc))
      phi1  = phi(npc) + c*ds
      c     = (dphi(2)-dphi(1))/(so(3)-so(2))
      dphi0 = dphi(1) - c*so(2)
      c     = (dphi(npc)-dphi(npc-1))/(so(npc+1)-so(npc))
      dphi1 = phi(npc) + c*ds
      c     = (dpht(2)-dpht(1))/(so(3)-so(2))
      dpht0 = dpht(1) - c*so(2)
      c     = (dpht(npc)-dpht(npc-1))/(so(npc+1)-so(npc))
      dpht1 = dpht(npc) + c*ds
*
*
* - sistemo nel nuovo array
* 
      pho(1)     = phi0
      pho(npc+2) = phi1
      dpho(1)    = dphi0
      dpho(npc+2)= dphi1
      dpto(1)    = dpht0
      dpto(npc+2)= dpht1
      tto(1)     = 0.d0
      yvo(1)     = yv(1)
      zvo(1)     = zv(1)
      do i = 1,npc
        dy = yv(i+1)-yv(i)
	dz = zv(i+1)-zv(i)
	dtt= sqrt(dy**2+dz**2)
	tto(i+1)  = tto(i)+dtt
	yvo(i+1)  = yv(i+1)
        zvo(i+1)  = zv(i+1)
        pho(i+1)  = phi(i)
        dpho(i+1) = dphi(i)
        dpto(i+1) = dpht(i)
      end do
*
* - interpolo il potenziale e le vel tangenti 
* 
      yp1=1.d+31
      ypn=1.d+31
      call splone( yvo2,tto, yv,yp1,ypn,npco+1,npamx )  
      call splone( zvo2,tto, zv,yp1,ypn,npco+1,npamx )  
      call splone( pho2,so, pho,yp1,ypn,npco+2,npamx )  
      call splone( dpho2,so, dpho,yp1,ypn,npco+2,npamx )  
      call splone( dpto2,so, dpto,yp1,ypn,npco+2,npamx )  
*
      do i=1,npco+1
        write(*,*) i,tto(i),yvo(i),yvo2(i)
      end do
      return
      end
*-----------------------------------------------------------
