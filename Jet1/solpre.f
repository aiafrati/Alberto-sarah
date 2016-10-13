
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine solpre(xv,zv,amp,npc,npt,dpnt,dpt,coef)

      include"slam_p.h"
      dimension xv(npamx+1),zv(npamx+1)
      dimension amp(npamx),dpt(npamx),dpnt(npamx)
      dimension aa(npamx,npamx),bb(npamx),appo(npamx)
      dimension indx(npamx)
      parameter (naux=3*npamx)
      character*1 trans

c     - inizializzo matrice coefficienti e vettore termini noti

      npe = npt+npc

      do ip = 1,npt
        xx = (xv(ip+1)+xv(ip))/2.d0
        zz = (zv(ip+1)+zv(ip))/2.d0
        bb(ip) = 0.d0
 
        do j = 1,npt+npc
          aa(ip,j) = 0.d0
        enddo

        do jp = 1,npc
          xp = xv(jp)
          xs = xv(jp+1)
          zp = zv(jp)
          zs = zv(jp+1)
  
          call finte( xx,zz,xp,zp,xs,zs,amp(jp),fint1,fint2,0.d0)

c         -- effetto pannello simmetrico
  
          xp = -xv(jp+1)
          xs = -xv(jp)
          zp = zv(jp+1)
          zs = zv(jp)

          call finte( xx,zz,xp,zp,xs,zs,amp(jp),fins1,fins2,0.d0)

c         -- definizione coefficienti matrice e termini noti

          aa(ip,jp) = - fint2 - fins2
          aa(ip,npt+jp) = - fint1 - fins1
        enddo

        do jp = npc+1,npt
          xp = xv(jp)
          xs = xv(jp+1)
          zp = zv(jp)
          zs = zv(jp+1)
  
          call finte( xx,zz,xp,zp,xs,zs,amp(jp),fint1,fint2,0.d0)

c         -- effetto pannello simmetrico
  
          xp = -xv(jp+1)
          xs = -xv(jp)
          zp = zv(jp+1)
          zs = zv(jp)

          call finte( xx,zz,xp,zp,xs,zs,amp(jp),fins1,fins2,0.d0)

c         -- definizione coefficienti matrice e termini noti

          aa(ip,jp) = - fint1 - fins1
          bb(ip)    = bb(ip) + dpt(jp)*(fint2+fins2)
        enddo

        if (ip.le.npc) then
          aa(ip,ip) = aa(ip,ip) + pi
        else
          bb(ip) = bb(ip)-pi*dpt(ip)
        end if
      enddo

c     - parte "implicita"

      do ip = npt+1,npt+npc
        bb(ip) = 0.d0
 
        do j = 1,npt+npc
          aa(ip,j) = 0.d0
        enddo

        do jp = 1,npc
          aa(ip,jp) = -coef*amp(jp)
        enddo
        aa(ip,ip) = 1.d0

        bb(ip) = dpnt(ip-npt)
      enddo
*
* ----------------- SOLUZIONE CON ROUTINE GAUSS ------------------------
      do i=1,npe
      aa(i,npe+1)=bb(i)
      end do
      call gauss(aa,npamx,npe,npe+1,1.d-12)
      do i=1,npe
      bb(i)=aa(i,npe+1)
      end do 
c
c     -------------   SOLUZIONE SISTEMA LINEARE CON ESSL   -------------

c      call dgefcd(aa,npamx,npe,indx,1,rcond,det,appo,npamx)
c      call dges(aa,npamx,npe,indx,bb,0) 

c     -------------   SOLUZIONE SISTEMA LINEARE CON LAPACK -------------

c      call DGETRF(npe,npe,aa,npamx,indx,idum)
c      trans='N'
c      call DGETRS(trans,npe,1,aa,npamx,indx,bb,npamx,idum)

c     -------------   FINE SOLUZIONE SISTEMA LINEARE       -------------

      do ip = 1,npc
        dpt(ip) = bb(ip)
      enddo
      do ip = npc+1,npt
        dpnt(ip) = bb(ip)
      enddo
      do ip = npt+1,npt+npc
        dpnt(ip-npt) = bb(ip)
      enddo

      return
      end

