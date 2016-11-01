
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine solv1(npc,nng,npsl,xv,zv,xgb,zgb,xsl,zsl,amp,ampsl,
     #                  phi,dphi,phisl,dphisl,jt)
*
      include"slam_p.h"
      dimension xv(npamx),zv(npamx),xsl(npamx),zsl(npamx)
      dimension xgb(npamx),zgb(npamx), amp(npamx),ampsl(npamx)
      dimension aa(npamx,npamx),bb(npamx),appo(npamx)
      dimension indx(npamx)
      dimension phi(npamx),dphi(npamx),phisl(npamx),dphisl(npamx)
      parameter (naux=3*npamx)
      character*1 trans
      write(*,*) '--> solv1...........'
c     - inizializzo matrice coefficienti e vettore termini noti
*
        npe = npsl+npc+nng
        write(*,*) npsl,npc,nng
*
        do  ip = 1,npe
          if(ip.le.npc)then
            xx = (xv(ip+1)+xv(ip))/2.d0
            zz = (zv(ip+1)+zv(ip))/2.d0
          elseif(ip.gt.npc.and.ip.le.npc+nng)then
            xx = (xgb(ip+1-npc)+xgb(ip-npc))/2.d0
            zz = (zgb(ip+1-npc)+zgb(ip-npc))/2.d0
          else
            xx = (xsl(ip+1-npc-nng)+xsl(ip-npc-nng))/2.d0
            zz = (zsl(ip+1-npc-nng)+zsl(ip-npc-nng))/2.d0
          endif
          bb(ip) = 0.d0
          do j = 1,npe
            aa(ip,j) = 0.d0
          enddo
          do jp = 1,npe
            if(jp.le.npc)then
            xp = xv(jp)
            xs = xv(jp+1)
            zp = zv(jp)
            zs = zv(jp+1)
            am = amp(jp)
            elseif(jp.gt.npc.and.jp.le.npc+nng)then
            xp = xgb(jp-npc)
            xs = xgb(jp+1-npc)
            zp = zgb(jp-npc)
            zs = zgb(jp+1-npc)
            am = amp(jp)
            else
            xp = xsl(jp-npc-nng)
            xs = xsl(jp+1-npc-nng)
            zp = zsl(jp-npc-nng)
            zs = zsl(jp+1-npc-nng)
            am = ampsl(jp-npc-nng)
            endif
            call finte( xx,zz,xp,zp,xs,zs,am,fint1,fint2,0.d0)
c           -- effetto pannello simmetrico
            if(jp.le.npc)then
            xs = -xv(jp)
            xp = -xv(jp+1)
            zs = zv(jp)
            zp = zv(jp+1)
            am = amp(jp)
            elseif(jp.gt.npc.and.jp.le.npc+nng)then
            xs = -xgb(jp-npc)
            xp = -xgb(jp+1-npc)
            zs = zgb(jp-npc)
            zp = zgb(jp+1-npc)
            am = amp(jp)
            else
            xs = -xsl(jp-npc-nng)
            xp = -xsl(jp+1-npc-nng)
            zs = zsl(jp-npc-nng)
            zp = zsl(jp+1-npc-nng)
            am = ampsl(jp-npc-nng)
            endif
            call finte( xx,zz,xp,zp,xs,zs,am,fins1,fins2,0.d0)
*
c           -- definizione coefficienti matrice e termini noti
            if(jp.le.npc+nng)then
              aa(ip,jp) = - fint2 - fins2
              bb(ip)    = bb(ip) + dphi(jp)*(fint1+fins1)
            else 
              aa(ip,jp) = - fint1 - fins1
              bb(ip)    = bb(ip) + phisl(jp-npc-nng)*(fint2+fins2)
            endif
            
*
          enddo
c       aggiungo il contributo autoindotto al termine noto o sulla
c       matrice dei coefficienti a seconda dei casi
          if(ip.le.npc+nng)then
            aa(ip,ip) = aa(ip,ip) + pi
          else
            bb(ip) = bb(ip)-pi*phisl(ip-npc-nng)
          endif
c         if(jt.eq.40)then
c          do i=1,npe
c          do j=1,npe
c          write(70,*)i,j,aa(i,j)
c          enddo 
c          write(70,*)
c          write(70,*)
c          enddo
c          stop
c        endif 
*
        enddo

*
* ----------------- SOLUZIONE CON ROUTINE GAUSS ------------------------
c      do i=1,npe
c      aa(i,npe+1)=bb(i)
c      end do
c      call gauss(aa,npamx,npe,npe+1,1.d-12)
c      do i=1,npe
c      bb(i)=aa(i,npe+1)
c      end do 
c


c
c     -------------   SOLUZIONE SISTEMA LINEARE CON ESSL   -------------

c      call dgefcd(aa,npamx,npe,indx,1,rcond,det,appo,npamx)
c      call dges(aa,npamx,npe,indx,bb,0) 

c     -------------   SOLUZIONE SISTEMA LINEARE CON LAPACK -------------

      call DGETRF(npe,npe,aa,npamx,indx,idum)
      trans='N'
      call DGETRS(trans,npe,1,aa,npamx,indx,bb,npamx,idum)

c     -------------   FINE SOLUZIONE SISTEMA LINEARE       -------------

      do ip = 1,npe
        if (ip.le.npc+nng) then 
          phi(ip)  = bb(ip)
        else
          dphisl(ip-npc-nng) = bb(ip)
        endif
      enddo
*
      return
      end
