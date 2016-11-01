*
      subroutine get(nget,npc,ycnsl,zcnsl,ycn,zcn,yn,zn,  
     #               xigs,zegs,xigb,zegb,xigf,zegf) 
*
      include"slam_p.h"
      dimension  ycnsl(npamx),zcnsl(npamx)
      dimension  ycn(npamx),zcn(npamx),yn(npamx),zn(npamx)
      dimension  xigs(-npamx:npamx), zegs(-npamx:npamx)
      dimension  xigb(-npamx:npamx), zegb(-npamx:npamx)
      dimension  xigf(-npamx:npamx), zegf(-npamx:npamx)
      write(*,*) '--> get..........'
* - aggiorno ze, h , dh, phit, dphtt
******************************************************************
c
c - rotazione sistema di riferimento (okkio ze ha segno invertito!)
      nb = 2
c      do i=1,nget+nb
c        za(i) = znsl(i)-proat
c        zac(i) = zcnsl(i)-proat 
c      enddo
c      za(nget+nb+1) = znsl(nget+nb+1)-proat 
c      call rot(ande,ynsl,za,yr,zr,nget+nb+1)
c      call rot(ande,ycnsl,zac,ycr,zcr,nget+nb)
* - calcolo h e dh e xi (nodi e centroidi)
c      do i=nget+nb,1,-1
c        h(nget-i+1)   = -zcr(i)
c        hh(nget-i+1)  = -zr(i+1)
c        dh(nget-i+1)  = -(zr(i)-zr(i+1))/(yr(i)-yr(i+1))
c - sup. corpo
c        xi(nget-i+1) = ycr(i)-yr(nget+nb+1)
c        xin(nget-i+1)= yr(i+1)-yr(nget+nb+1)
c      enddo
c      hh(nget+1) = -zr(1)
c      xin(nget+1)= yr(1)-yr(nget+nb+1)
* - calcolo xig  xigs  zegs zeb zef (per solv2)
      do i= -nb,nget
c        xig(i) = xi(i)
c        zeb(i) = 0.d0
c        zef(i) = h(i)
        xigb(i) = ycn(npc+i)
        zegb(i) = zcn(npc+i)
c        xigb(i) = 0.5d0*(yn(npc+i)+yn(npc+i+1))
c        zegb(i) = 0.5d0*(zn(npc+i)+zn(npc+i+1))
        xigf(i) = ycnsl(nget-i+1)
        zegf(i) = zcnsl(nget-i+1)
      enddo
      write(95,*) '# ',jt
      do i = -nb+1,nget
c        xigs(i) = xi(i)
c        zegs(i) = 0.d0 
c        xigs(i) = 0.5d0*(xi(i)+xi(i-1))
c        zegs(i) = 0.25d0*(h(i)+h(i-1))
        xigs(i) = 0.25d0*(xigb(i)+xigb(i-1)+xigf(i)+xigf(i-1))
        zegs(i) = 0.25d0*(zegb(i)+zegb(i-1)+zegf(i)+zegf(i-1)) 
        write(95,'(i6,10d15.6)') 
     #      i,xigs(i),zegs(i),
     #      xigb(i),zegb(i),xigf(i),zegf(i),
     #      xigb(i-1),zegb(i-1),xigf(i-1),zegf(i-1)
      enddo  
      write(95,*)
      write(95,*)
*
      return
      end
