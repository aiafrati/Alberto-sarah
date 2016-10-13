
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine caveta(phi,amm,dfit,npc,npt,kcut,dphin)

      include"slam_p.h"
      
      dimension phi(npamx),dfit(npamx),amm(npamx),dphin(npamx)
*
      write(*,*) '.......... :  caveta'
*
c     superficie libera
*
      do i = npc+kcut+1,npt
        if (i.eq.npc+kcut+1) then
          dfit(i) = 2.d0*(phi(i+1)-phi(i))/(amm(i+1)+amm(i))
        else if (i.eq.npt) then
          dfit(i) = 2.d0*(phi(i)-phi(i-1))/(amm(i-1)+amm(i))
        else
          dfif = 2.d0*(phi(i+1)-phi(i))/(amm(i+1)+amm(i))
          dfib = 2.d0*(phi(i)-phi(i-1))/(amm(i-1)+amm(i))
          dfit(i) = (dfif+dfib)/2.d0
        end if
      enddo
*
c     superficie corpo 
*
      do i = 1,npc
        if (i.eq.1) then
          dfit(i) = 2.d0*(phi(i+1)-phi(i))/(amm(i+1)+amm(i))
        else if (i.eq.npc) then
          dfit(i) = 2.d0*(phi(i)-phi(i-1))/(amm(i-1)+amm(i))
        else
          dfif = 2.d0*(phi(i+1)-phi(i))/(amm(i+1)+amm(i))
          dfib = 2.d0*(phi(i)-phi(i-1))/(amm(i-1)+amm(i))
          dfit(i) = (dfif+dfib)/2.d0
        end if
      enddo
*
c       - velocita' tang. tappo 
*
        if (kcut.eq.1) then
          dfit(npc+1) = dphin(npc)
        end if
*
      return
      end
*
