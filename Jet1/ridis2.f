*
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine ridis2(k,proat,escr,yn,zn,ygn,ygs2,zgn,zgs2,ng,
     #                  amp,iint,tag, kcut,npt,npc,phin,dphi,
     #                   ycn,zcn,vfall)
*
      include"slam_p.h"
*      include"slam_v.h"
      dimension yn(npamx),zn(npamx),ygn(npamx+1),zgn(npamx+1)
      dimension ygs2(npamx+1),zgs2(npamx+1),tag(npamx)
      dimension phin(npamx),dphi(npamx),amp(npamx)
      dimension ycn(npamx),zcn(npamx)
      write(*,*)  '.......... :  ridis2'
*
      kk   = k
      eskk = escr
*
      yy   = yn(npc+1)
      zz   = zn(npc+1)
      yf   = 0.d0
      zf   = proat
      dy   = yy-ygn(iint)
      dz   = zz-zgn(iint)
      tai  = tag(iint)+sqrt(dy**2+dz**2)
*
      amdi = tai
*
      amii = sqrt((yn(npc+2)-yn(npc+1))**2+
     &            (zn(npc+2)-zn(npc+1))**2 )
*
      if (kcut.eq.1) then
        amsu = sqrt((yn(npc+3)-yn(npc+2))**2+
     &              (zn(npc+3)-zn(npc+2))**2 )
        amii = min(amsu,amii)
      end if
*
      if (amdi.lt.(1.d0+eskk)*amii) kk = 0
*       
      if (eskk.gt.1.d0) then
        npcn = int( log( 1.d0+(eskk-1.d0)*amdi/amii)/log(eskk) )
      else 
        npcn = int( amdi/amii )
      end if
*
      if (npcn.gt.npc.and.npt.ge.npamx) then
        write(*,*) '  !!! RAGGIUNTE DIMENSIONI MASSIME, REGRID BLOCCATO'
        kk = 0
      end if
*
      if (kk.eq.0) then
        npcn = npc
        amii = amdi
        if (npcn.gt.1.and.eskk.gt.1.d0) then
          amii = amdi*(1.d0-eskk)/(1.d0-eskk**npcn)
        end if
*
        yn(npcn+1) = yy
        zn(npcn+1) = zz
        write(*,*) 'pip ', yy,zz
*
        r  = tai
        amii1 = amii
*        do iv =npcn+1,1,-1
        do iv =npcn,1,-1
          r = r - amii
          if(iv.eq.1)      r = 0.d0
*          if(iv.eq.npcn+1) r = tai
          call splont(r,yyy,tag,ygn,ygs2,ng,npamx,
     #                   xnull,0,xnull,0,xnull,0)
          call splont(r,zzz,tag,zgn,zgs2,ng,npamx,
     #                    xnull,0,xnull,0,xnull,0)
          yn(iv) = yyy 
          zn(iv) = zzz 
          if(iv.le.npcn) then
              ampy       = yn(iv+1) - yn(iv)
              ampz       = zn(iv+1) - zn(iv)
              amp(iv)    = sqrt(ampy**2 + ampz**2)
              rnz        = -ampy/amp(iv)
              dphi(iv)   = -vfall*rnz
          endif
          if (iv.le.npcn) amii=amii*eskk
        end do
*
*---------------------------------------------------------------------
*
      else
*
        ninc = npcn - npc
*
        if (kcut.eq.0.and.ninc.lt.0) then
          ninc = 0 ! non li riduco prima del taglio
          npcn = npc
        end if
*
        amii = amdi 
        if (npcn.gt.1) amii = amdi*(1.d0-eskk)/(1.d0-eskk**npcn)
*
        if (ninc.gt.0) then
c         nrid=0
          write(*,*) '   aggiungo pannelli sul corpo '
          yn(npt+ninc+1) = yn(npt+1)
          zn(npt+ninc+1) = zn(npt+1)
          do ip = npt,npc+1,-1
            ic        = ip + ninc
            phin(ic)  = phin(ip)
            dphi(ic)  = dphi(ip)
            ycn(ic)   = ycn(ip)
            zcn(ic)   = zcn(ip)
            yn(ic)    = yn(ip) 
            zn(ic)    = zn(ip) 
          enddo
        else if (ninc.lt.0) then
c         nrid=0
          write(*,*) '   taglio pannelli dal corpo '
          do ip = npc+1,npt
            ic        = ip + ninc
            phin(ic)  = phin(ip)
            dphi(ic)  = dphi(ip)
            ycn(ic)   = ycn(ip)
            zcn(ic)   = zcn(ip)
            yn(ic)    = yn(ip) 
            zn(ic)    = zn(ip) 
          enddo
          yn(npt+ninc+1) = yn(npt+1)
          zn(npt+ninc+1) = zn(npt+1)
        end if
*        
c       - ridefinizione dei pannelli sul corpo
*
        yn(npcn+1) = yy
        zn(npcn+1) = zz
*        
*        write(*,*) 'pi2p ', yy,zz
        amii1=amii
        r  = tai
*        do iv =npcn+1,1,-1
        do iv =npcn,1,-1
           r = r - amii
           if(iv.eq.1)      r = 0.d0
*           if(iv.eq.npcn+1) r = tai
           call splont(r,yyy,tag,ygn,ygs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           call splont(r,zzz,tag,zgn,zgs2,ng,npamx,
     #                 xnull,0,xnull,0,xnull,0)
           yn(iv) = yyy
           zn(iv) = zzz
           if(iv.le.npcn) then
              ampy       = yn(iv+1) - yn(iv)
              ampz       = zn(iv+1) - zn(iv)
              amp(iv)    = sqrt(ampy**2 + ampz**2)
              rnz        = -ampy/amp(iv)
              dphi(iv)   = -vfall*rnz
           endif
           if (iv.le.npcn) amii=amii*eskk
        end do 
*
        npc = npcn
        npt = npt + ninc
      end if 
*
      return
      end
*
