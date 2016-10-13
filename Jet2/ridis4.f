
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ridis4(k,ng,proat,kget,ynsl,znsl,ygb,zgb,
     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
     #                  nsep,nsepo,ksep,phg,phgs2,phin)
*
      include"slam_p.h"
      dimension ygb(npamx),zgb(npamx),ynsl(npamx),znsl(npamx)
      dimension yn(npamx),zn(npamx),ycn(npamx),zcn(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx)
      dimension tg(npamx),tgb(npamx),ksep(npamx) 
      dimension phg(npamx),phgs2(npamx),phin(npamx)
      dimension tn(npamx),tpo(npamx),phpo(npamx),phpo2(npamx)
*
      write(*,*) '--> ridis4...........'
*
      kk   = k
      eskk = escr
      ninc = 0
      nng  = ng*kget
      epss = 1.d-10
*
      yf = 0.d0
      zf = proat
      if(ng.eq.0)then
        yy   = ynsl(1)
        zz   = znsl(1) 
        amii = sqrt((ynsl(2)-ynsl(1))**2+ (znsl(2)-znsl(1))**2 )
        ay   = yy - ygn(iint) 
        az   = zz - zgn(iint)
        aa   = sqrt(ay**2+az**2)
        amdi = tg(iint) + aa
        tt   = amdi 
      write(*,*) '--> ammi amdi 1 .',amii,amdi
      else
        if(kk.eq.0)then
* - con kk=0 ricalcolo vertici getto, con kk=1 no
*
        do ii=2,ng+1
          y1  = ynsl(ii)
          z1  = znsl(ii)
          jj  = ng-ii+2
          jjj  = npc+ng+2-ii
          if(ksep(jjj).ne.1)then
c          else
          do i=1,ngo1-1
             w1y = zgn(i+1) - zgn(i)
             w1z = -(ygn(i+1) - ygn(i))
             w2y = ygn(i+1) - ygn(i)
             w2z = zgn(i+1) - zgn(i)
             y2  = ygn(i)
             z2  = zgn(i)
             call intret(yi,zi,r1,r2,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
             if(r2.ge.(0.d0-epss).and.r2.le.(1.d0+epss))then
               kint = i
               ay = yi - ygn(kint)
               az = zi - zgn(kint)
               aa = sqrt(ay**2+az**2)
               ta = tg(kint) + aa
               call splint(ta,yb,zb,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
               ygb(jj) = yb 
               zgb(jj) = zb 
               tgb(jj) = ta
               goto 999
             endif
          enddo
C - DA INSERIRE LA DEFINIZIONE DEI PUNTI GETTO NON RIGRIGLIATI (ksep=1)
C          else
          endif
  999     continue
        enddo
*
        endif
* - se kk=1 ygb e zgb  (e ampli) e tgb gia' calcolati in shallo
        ygb(ng+1) = ynsl(1)
        zgb(ng+1) = znsl(1)
        yy = ygb(1)
        zz = zgb(1)
        tt = tgb(1)
        amii = ampli
        amdi = tgb(1) 
      endif
* - !! eskk =! 1.d0  !!!
      nn = int( log( 1.d0+(eskk-1.d0)*amdi/amii)/log(eskk) )
      no = npc - ng*(1-kget)
* - calcolo differenza 
      ninc = nn+ng*(1-kget)-npc
c - non cambio se kk=0
      if(kk.eq.0.and.ng.gt.0) ninc = 0
* - non li riduco se ng=0 (cioe prima di iiget) 
      if(ninc.lt.0.and.ng.eq.0)then
        nn   = no
        ninc = 0
      elseif(ninc.lt.0)then
        write(*,*) 'taglio   ',ninc,' pannelli dal corpo'
      elseif(ninc.gt.0)then
        write(*,*) 'aggiungo ',ninc,' pannelli sul corpo'
      else
        write(*,*) 'npc costante, ninc= ',ninc,ng,kget
        write(*,*) yy,zz
      endif
* - ricalcolo amii
      amii = amdi*(1.d0-eskk)/(1.d0-eskk**nn)
      write(*,*) 'amii amdi2 ',amii,amdi,nn
* - griglio bulk
      yn(nn+1) = yy
      zn(nn+1) = zz
      r  = tt
      tn(nn+1) = r
      dr = amii
      write(*,'(i4,3d15.7)') nn+1,r, yy,zz
      write(*,*) 'npc ng' , npc,ng,nn
      do ip = nn,1,-1
        r        = r-dr
        if(nn.eq.1) r = 0.d0
        tn(ip) = r 
        call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        yn(ip)   = y 
        zn(ip)   = z
        dr       = dr*eskk
        write(*,'(i4,3d15.7)') ip,r, y,z
      enddo
      write(*,*) 'PIPPO',npc
      npco = npc
      npc  = nn  + ng*(1-kget)
      npt  = npt + ninc  
      write(*,*) 'PIPPO2',npc
* - aggiungo getto a yn zn, se ng>0 e kget = 0
      if(ng.gt.0)then
        do i=2,ng+1
          yn(nn+i) = ygb(i)  
          zn(nn+i) = zgb(i)
          tn(nn+i) = tgb(i)
        write(*,'(i4,3d15.7)') npc+i,tgb(i), ygb(i),zgb(i)
        enddo
c        if(nsepo.eq.1) ksep(nn+ng)=0
      endif
*
      if(kk.eq.0)then
*
        write(*,*) 'RIDE 1 kk=0',ngo1
c - aggiorno solo centroidi del corpo, se kk=0
      tt = 0.d0
      do i = 1,npc+nng
       if(ksep(i).ne.1)then
c        if(i.le.npc)then
c        ycn(i) = (yn(i+1)+yn(i))/2.d0
c        zcn(i) = (zn(i+1)+zn(i))/2.d0
c        else
c        ii = i-npc
c        ycn(i) = (ygb(ii+1)+ygb(ii))/2.d0
c        zcn(i) = (zgb(ii+1)+zgb(ii))/2.d0
c        endif
c       else
        tt = 0.5d0*(tn(i)+tn(i+1)) 
        call splint(tt,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        ycn(i)  = y 
        zcn(i)  = z 
        write(*,*) i,y,z
       endif
      enddo
*
      else
*
c --  ascissa curvilinea
      tpo(npco+nng+1) = tn(npc+nng+1)
      npo  = npco+nng+1
      tt   = 0.d0  
      do i=1,npco+nng
        if(i.eq.1)then
          dy=ycn(i)-yn(i)
          dz=zcn(i)-zn(i)
        else
          dy=ycn(i)-ycn(i-1)
          dz=zcn(i)-zcn(i-1)
        endif
        dd     = sqrt(dy**2+dz**2)
        tt     = tt+dd
        tpo(i) = tt
        phpo(i)= phin(i)
        write(*,*) i,tt,phin(i)
      enddo
* -- estrapolo potenziale
      dtt  = tpo(npco+nng)-tpo(npco+nng-1)
      phia = phin(npco+nng)
      phib = phin(npco+nng-1)
      delphi = (phia-phib)/dtt 
      phie = phia+delphi*(tpo(npco+nng+1)-tpo(npco+nng))
      phpo(npco+nng+1) = phie
      write(*,*) npco+nng+1,tpo(npco+nng+1),phie
        write(*,*) 'RIDE 1',i
cc INTERPOLAZ SPLINE POTENZIALKE>>
      yp1 = 1.d31
      ypn = 1.d31
      call spline1(phpo,phpo2,tpo,yp1,ypn,npo,npamx)
c -- aggiorno divisione corpo SL, e tutti i centroidi ,se kk=1
      tt = 0.d0
      write(*,*) 'ppp '
      do i = 1,npc+nng
        tt = 0.5d0*(tn(i)+tn(i+1))
        write(*,*) ,i,tt 
        call splint(tt,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        if(tt.gt.tg(ngo))then
          if(ksep(i-1).eq.0)then
            ksep(i) = 2
            ichin   = i
          else
            ksep(i) = 1
          endif 
        else
          ksep(i) = 0
        endif
        ycn(i)  = y 
        zcn(i)  = z
        if(i.gt.1)then
        call splint1(tt,ph,phpo,phpo2,tpo,npo,npamx)
        phin(i) = ph 
        endif
c        write(*,*) i,y,z,ph
      enddo
*
      endif     

        
      return
      end
