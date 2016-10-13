      subroutine nortan(vfall,yn,zn,amp,tmy,tmz,rny,rnz,npc,
     #            ynsl,znsl,ampsl,tmysl,tmzsl,rnysl,rnzsl,npsl,
     #            dphi,kget,ng,ycn,zcn,yce,zce,yv,zv,
     #            ycnsl,zcnsl,ycsl,zcsl,ysl,zsl,phinsl,phisl,
     #            phin,phi,ygb,zgb,k2,ksep,jt)
*
      include "slam_p.h"
      dimension yn(npamx),zn(npamx),amp(npamx),tmy(npamx),tmz(npamx)
      dimension rny(npamx),rnz(npamx),ynsl(npamx),znsl(npamx)
      dimension ampsl(npamx),tmysl(npamx),tmzsl(npamx),rnysl(npamx)
      dimension rnzsl(npamx),kphi(npamx),ycn(npamx),zcn(npamx)
      dimension yce(npamx),zce(npamx),yv(npamx),zv(npamx),dphi(npamx)
      dimension ycnsl(npamx),zcnsl(npamx),ycsl(npamx),zcsl(npamx)
      dimension ysl(npamx),zsl(npamx),phinsl(npamx),phisl(npamx)
      dimension ygb(npamx),zgb(npamx),phi(npamx),phin(npamx)
      dimension ksep(npamx)
*
      write(*,*) '--> nortan......'
*
*- aggiorno ampiezze, normali tangenti e  boundary conditions
*
        nng=kget*ng
        do ip = 1,npc
          amy = yn(ip+1) - yn(ip)
          amz = zn(ip+1) - zn(ip)
          am  = sqrt(amy*amy + amz*amz)
          amp(ip) =  am
          tmy(ip) =  amy/am
          tmz(ip) =  amz/am
          rny(ip) =  tmz(ip)
          rnz(ip) = -tmy(ip)
          dphi(ip)= -vfall*rnz(ip)
c          endif
        enddo
* -getto , solo ampiezze
        do ip = 1,nng
          amy = ygb(ip+1)-ygb(ip)
          amz = zgb(ip+1)-zgb(ip)
          am  = sqrt(amy*amy + amz*amz)
          amp(npc+ip) = am
          tmy(npc+ip) =  amy/am
          tmz(npc+ip) =  amz/am
          rny(npc+ip) =  tmz(npc+ip)
          rnz(npc+ip) = -tmy(npc+ip)
          dphi(npc+ip)=  -vfall*rnz(npc+ip)
c          endif
        enddo
* - SL
        do ip = 1,npsl
          amy = ynsl(ip+1) - ynsl(ip)
          amz = znsl(ip+1) - znsl(ip)
          am  = sqrt(amy*amy + amz*amz)
          ampsl(ip) = am
          tmysl(ip) =  amy/am
          tmzsl(ip) =  amz/am
          rnysl(ip) =  tmzsl(ip)
          rnzsl(ip) = -tmysl(ip)
        enddo
*- eventualmetnte riposiziono array nuovi in quelli vecchi
        if(k2.eq.1)then
          do ip=1,npc+nng
            yv(ip) = yn(ip)
            zv(ip) = zn(ip)
            yce(ip) = ycn(ip)
            zce(ip) = zcn(ip)
            phi(ip) = phin(ip)
          enddo
          yv(npc+nng+1) = yn(npc+nng+1)
          zv(npc+nng+1) = zn(npc+nng+1)
          yv(npc+1) = yn(npc+1)
          zv(npc+1) = zn(npc+1)
          do ip=1,npsl
            ysl(ip) = ynsl(ip)
            zsl(ip) = znsl(ip)
            ycsl(ip) = ycnsl(ip)
            zcsl(ip) = zcnsl(ip)
            phisl(ip) = phinsl(ip)
          enddo
          ysl(npsl+1) = ynsl(npsl+1)
          zsl(npsl+1) = znsl(npsl+1)
        endif
*
        return
        end 
