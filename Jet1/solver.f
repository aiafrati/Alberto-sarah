
* ---------------------------------------------------------------- ------


c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
    
      subroutine solver(
     #   jt,xv,zv,amp,phi,dphi,kphi,npt,npc,kcut,t,w,n,i2d)
*
      include"slam_p.h"
      dimension t(ngmax),w(ngmax)
      dimension xv(npamx+1),zv(npamx+1)
      dimension amp(npamx),phi(npamx),dphi(npamx),kphi(npamx+1)
      dimension aa(npamx,npamx),bb(npamx),appo(npamx)
      dimension indx(npamx),det(2)
      parameter (naux=3*npamx)
      character*1 trans
*
      write(*,*) '.......... :  solver'
*
c     - inizializzo matrice coefficienti e vettore termini noti
*     
      if (kcut.eq.0) then
        npe = npt
      if(i2d.eq.1)then
        cost = 1.d0
        do ip = 1,npt
          xx = (xv(ip+1)+xv(ip))/2.d0
          zz = (zv(ip+1)+zv(ip))/2.d0
          bb(ip) = 0.d0
          do j = 1,npt
            aa(ip,j) = 0.d0
          enddo
*          
           do jp = 1,npt
            xp = xv(jp)
            xs = xv(jp+1)
            zp = zv(jp)
            zs = zv(jp+1)
            call finte( xx,zz,xp,zp,xs,zs,amp(jp),fint1,fint2,0.d0)
c   -- effetto pannello simmetrico
            xp = -xv(jp+1)
            xs = -xv(jp)
            zp = zv(jp+1)
            zs = zv(jp)
            call finte( xx,zz,xp,zp,xs,zs,amp(jp),fins1,fins2,0.d0)
c   -- definizione coefficienti matrice e termini noti
            if (kphi(jp).gt.0) then
              aa(ip,jp) = - fint1 - fins1
              bb(ip)    = bb(ip) + phi(jp)*(fint2+fins2)
            else
              aa(ip,jp) = - fint2 - fins2
              bb(ip)    = bb(ip) + dphi(jp)*(fint1+fins1)
            end if
           enddo
c   -- aggiungo il contributo autoindotto al termine noto o sulla
c   -- matrice dei coefficienti a seconda dei casi
          if (kphi(ip).eq.0) then
            aa(ip,ip) = aa(ip,ip) + cost*pi
          else 
            bb(ip) = bb(ip)-cost*pi*phi(ip)
          end if
        end do
       elseif(i2d.eq.0)then
         cost=2.d0 
         ammin=amp(1)
         do i = 1,npt
           ammin=min(amp(i),ammin)
         end do
         eps = 0.d0*ammin
         eps2=eps**2
         do ip = 1,npt
          xx = (xv(ip+1)+xv(ip))/2.d0
          zz = (zv(ip+1)+zv(ip))/2.d0
          bb(ip) = 0.d0
          do j = 1,npt
            aa(ip,j) = 0.d0
          enddo
          do jp = 1,npt
            xp = xv(jp)
            xs = xv(jp+1)
            zp = zv(jp)
            zs = zv(jp+1)
            call cinfax(eps2,xx,zz,xp,zp,xs,zs,t,w,n,cinfg,cinfdg)
            if (kphi(jp).eq.0) then
              aa(ip,jp) = - cinfdg 
              bb(ip)    = bb(ip) + dphi(jp)*cinfg
            else
              aa(ip,jp) = - cinfg 
              bb(ip)    = bb(ip) + phi(jp)*cinfdg
            end if
          end do
c --- aggiungo il contributo autoindotto al termine noto o sulla
c ---  matrice dei coefficienti a seconda dei casi
          if (kphi(ip).eq.0) then
            aa(ip,ip) = aa(ip,ip) + cost*pi
          else 
            bb(ip) = bb(ip)-cost*pi*phi(ip)
          end if
        enddo
      endif
*
      else if (kcut.eq.1) then
        npe = npt+1
*
      if(i2d.eq.1)then
        cost=1.d0
        do ip = 1,npt
          xx = (xv(ip+1)+xv(ip))/2.d0
          zz = (zv(ip+1)+zv(ip))/2.d0
          bb(ip) = 0.d0
          do j = 1,npe
            aa(ip,j) = 0.d0
          enddo
          do jp = 1,npt
            xp = xv(jp)
            xs = xv(jp+1)
            zp = zv(jp)
            zs = zv(jp+1)
            call finte( xx,zz,xp,zp,xs,zs,amp(jp),fint1,fint2,0.d0)
c    -- effetto pannello simmetrico
            xp = -xv(jp+1)
            xs = -xv(jp)
            zp = zv(jp+1)
            zs = zv(jp)
            call finte( xx,zz,xp,zp,xs,zs,amp(jp),fins1,fins2,0.d0)
c   -- definizione coefficienti matrice e termini noti
            if (kphi(jp).eq.2) then
              aa(ip,jp)    = - fint2 - fins2
              aa(ip,npt+1) = - fint1 - fins1
            else if (kphi(jp).eq.1) then
              aa(ip,jp) = - fint1 - fins1
              bb(ip)    = bb(ip) + phi(jp)*(fint2+fins2)
            else
              aa(ip,jp) = - fint2 - fins2
              bb(ip)    = bb(ip) + dphi(jp)*(fint1+fins1)
            end if
          end do
c -- aggiungo il contributo autoindotto al termine noto o sulla
c -- matrice dei coefficienti a seconda dei casi
          if (kphi(ip).eq.1) then
            bb(ip) = bb(ip)-cost*pi*phi(ip)
          else 
            aa(ip,ip) = aa(ip,ip) + cost*pi
          end if
        end do
      elseif(i2d.eq.0)then
        cost=2.d0 
        ammin=amp(1)
        do i = 1,npt
         ammin=min(amp(i),ammin)
        end do
        eps = 0.d0*ammin
        eps2=eps**2
        do ip = 1,npt
          xx = (xv(ip+1)+xv(ip))/2.d0
          zz = (zv(ip+1)+zv(ip))/2.d0
          bb(ip) = 0.d0
          do j = 1,npe
            aa(ip,j) = 0.d0
          enddo
          do jp = 1,npt
            xp = xv(jp)
            xs = xv(jp+1)
            zp = zv(jp)
            zs = zv(jp+1)
            call cinfax(eps2,xx,zz,xp,zp,xs,zs,t,w,n,cinfg,cinfdg)
            if (kphi(jp).eq.0) then
              aa(ip,jp) = - cinfdg 
              bb(ip)    = bb(ip) + dphi(jp)*cinfg
            elseif(kphi(jp).eq.1) then
              aa(ip,jp) = - cinfg 
              bb(ip)    = bb(ip) + phi(jp)*cinfdg
	    else
	      aa(ip,jp)    = - cinfdg
              aa(ip,npt+1) = - cinfg
            end if
          end do
c -- aggiungo il contributo autoindotto al termine noto o sulla
c -- matrice dei coefficienti a seconda dei casi
          if (kphi(ip).eq.1) then
            bb(ip) = bb(ip)-cost*pi*phi(ip)
          else 
            aa(ip,ip) = aa(ip,ip) + cost*pi
          end if
        end do
      endif

        do j = 1,npt+1
          aa(npe,j) = 0.d0
        end do

c       impongo derivata normale tappo come media dei valori 
c       estrapolati dal lato superficie libera e corpo

        ag1 = amp(npc)/2.d0
        ag2 = amp(npc) + amp(npc-1)/2.d0
        ag3 = amp(npc) + amp(npc-1) + amp(npc-2)/2.d0
        bg1 = amp(npc+2)/2.d0
        bg2 = amp(npc+2) + amp(npc+3)/2.d0
        bg3 = amp(npc+2) + amp(npc+3) + amp(npc+4)/2.d0
        ca1 = (ag2**2-ag3**2)/( (ag1-ag3)*(ag1-ag2)*(ag2-ag3) )
        ca2 = (ag1**2-ag2**2)/( (ag1-ag3)*(ag1-ag2)*(ag2-ag3) )
        cb1 = (bg2**2-bg3**2)/( (bg1-bg3)*(bg1-bg2)*(bg2-bg3) )
        cb2 = (bg1**2-bg2**2)/( (bg1-bg3)*(bg1-bg2)*(bg2-bg3) )

        amy =  xv(npc+3) - xv(npc+2)
        amz =  zv(npc+3) - zv(npc+2)
        ty2 =  amy/amp(npc+2)
        tz2 =  amz/amp(npc+2)
        ry2 =  tz2
        rz2 = -ty2
        amy =  xv(npc+2) - xv(npc+1)
        amz =  zv(npc+2) - zv(npc+1)
        ty1 =  amy/amp(npc+1)
        tz1 =  amz/amp(npc+1)
        ry1 =  tz1
        rz1 = -ty1
        aaa = ry1*ty2 + rz1*tz2
        bbb = ry1*ry2 + rz1*rz2

        aa(npe,npe)   = 1.d0

c       - velocita' tappo da continuita'

*        do i = npc+2,npt
*          aa(npe,i) = amp(i)/amp(npc+1)
*        enddo

*        bba = 0.d0
*        do i = 1,npc
*          bba = bba-amp(i)*dphi(i)
*        enddo
*        bb(npe) = bba/amp(npc+1)

c       - velocita' tappo come valore estrapolato lato corpo

*        aa(npe,npc)   = ca1 
*        aa(npe,npc-1) = -(ca1+ca2)
*        aa(npe,npc-2) = ca2
*        bb(npe)       = 0.d0 

c       - velocita' tappo come valore lato corpo

*        aa(npe,npc)   =  2.d0/(amp(npc)+amp(npc-1))
*        aa(npe,npc-1) = -2.d0/(amp(npc)+amp(npc-1))
*        bb(npe)       = 0.d0 

c       - velocita' tappo come valore estrapolato lato sup. lib.

        aa(npe,npc+2) = -bbb

        bb(npe)       = -aaa*( cb1*(phi(npc+2)-phi(npc+3)) -
     &                         cb2*(phi(npc+3)-phi(npc+4)) )

c       - velocita' tappo come media valore lato corpo e sup. lib

*        aa(npe,npc)   = ca1/2.d0
*        aa(npe,npc-1) = -(ca1+ca2)/2.d0
*        aa(npe,npc-2) = ca2/2.d0
*        aa(npe,npc+2) = -bbb/2.d0
*        bb(npe)       = -aaa*( cb1*(phi(npc+2)-phi(npc+3)) -
*     &                         cb2*(phi(npc+3)-phi(npc+4)) )/2.d0

      end if

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
c     -------------   SOLUZIONE SISTEMA LINEARE CON ESSL   -------------

c      call dgefcd(aa,npamx,npe,indx,1,rcond,det,appo,npamx)
c      call dges(aa,npamx,npe,indx,bb,0)
c      write(*,*) 'rcond, det ',sngl(rcond),sngl(det(1)),sngl(det(2)) 

c     -------------   SOLUZIONE SISTEMA LINEARE CON LAPACK -------------

      call DGETRF(npe,npe,aa,npamx,indx,idum)
      trans='N'
      call DGETRS(trans,npe,1,aa,npamx,indx,bb,npamx,idum)

c     -------------   FINE SOLUZIONE SISTEMA LINEARE       -------------

      do ip = 1,npt
        if (kphi(ip).eq.1) then 
          dphi(ip)  = bb(ip)
        else
          phi(ip)  = bb(ip)
        endif
      enddo

      if (kcut.eq.1) dphi(npc+1) = bb(npe)

      return
      end

