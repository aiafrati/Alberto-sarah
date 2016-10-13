
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ridis5(k,ng,proat,kget,ynsl,znsl,ygb,zgb,
     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
     #                  nsep,nsepo,ksep,phg,phgs2,phin,
     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)
*
      include"slam_p.h"
      dimension ygb(npamx),zgb(npamx),ynsl(npamx),znsl(npamx)
      dimension yn(npamx),zn(npamx),ycn(npamx),zcn(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx)
      dimension tg(npamx),tgb(npamx),ksep(npamx) 
      dimension phg(npamx),phgs2(npamx),phin(npamx)
      dimension tn(npamx),tpo(npamx),phpo(npamx),phpo2(npamx)
      dimension ycnsl(npamx),zcnsl(npamx)
      dimension ycb(npamx),zcb(npamx),tcb(npamx)
      dimension phinsl(npamx),tc(npamx),kord(npamx)
      dimension rft(npamx)
*
      write(*,*) '--> ridis5...........'
*
      ne   = 2 
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
* -  ricalcolo centroidi getto
ccc        write(25,*) ' # jt ng',jt,ng
        do ii=1,ng
          y1  = ycnsl(ii)
          z1  = zcnsl(ii)
          jj  = ng-ii+1
          jjj  = npc+jj
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
               ycb(jj) = yb 
               zcb(jj) = zb 
               tcb(jj) = ta
ccc               write(25,'(i4,3d15.7,i4)') jj,yb,zb,ta,0
               goto 999
             endif
          enddo
C - DA INSERIRE LA DEFINIZIONE DEI PUNTI GETTO NON RIGRIGLIATI (ksep=1)
          else
             ycb(jj) = ycn(jjj) 
             zcb(jj) = zcn(jjj) 
             tcb(jj) = tc(jjj) 
ccc               write(25,'(i4,3d15.7,i4)') jj,ycb(jj),zcb(jj),tcb(jj),1
          endif
  999     continue
        enddo
ccc         write(25,*)
ccc         write(25,*)
* - calcolo vertici getto
ccc        write(26,*) ' # jt ng',jt,ng
c        write(26,'(i4,3d15.7)') 1,ygb(1),zgb(1),tgb(1)
        do i=1,ng-1
          tb = 0.5d0*(tcb(i)+tcb(i+1))
          call splint(tb,yb,zb,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          ygb(i+1) = yb 
          zgb(i+1) = zb 
          tgb(i+1) = tb 
ccc          write(26,'(i4,3d15.7)') i+1,yb,zb,tb
        enddo
        if(nsep.eq.0)then
          ygb(ng+1) = ynsl(1) 
          zgb(ng+1) = znsl(1)
        else
          ygb(ng+1) = ygn(ngo1) 
          zgb(ng+1) = zgn(ngo1)
        endif
        dy = ygb(ng+1)-ygb(ng) 
        dz = zgb(ng+1)-zgb(ng)
        dd = sqrt(dy**2+dz**2)
        tgb(ng+1) = tgb(ng)+dd 
ccc        write(26,'(i4,3d15.7)') ng+1,ygb(ng+1),zgb(ng+1),tgb(ng+1)
ccc        write(26,*)
ccc        write(26,*)
*
        endif
* - se kk=1 ygb e zgb  (e ampli) e tgb gia' calcolati in shallo
        yy = ycb(1)
        zz = zcb(1)
        tt = tcb(1)
        amii = ampli
        amdi = tcb(1)-0.5d0*amii 
      endif
********************
      nn1 = 0
      if(amdi-20.d0*ramii*amii.gt.tg(ngo).or.ksup.eq.1)then
* - non voglio cohe cominci con kk=0
        if(ksup.eq.0.and.kk.eq.0) goto 98
        ksup = 1
        write(*,*) 'POPOPO' 
        t1 = tg(ngo)
        t0 = 0.5d0*ramiii*amii
        tl = t1-t0
        a1 = ramiii*amii/tl
        a2 =  ramii*amii/tl
        et1 = eskkk
        et2 = eskkk
        n1 = nn/2
        n2 = nn-n1
        if(kk.eq.0) then
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,kkk)
          if(nn.ne.nnold)then
            nn  =  nnold
            n2  =  nn - n1
            if(n2.lt.0)then 
              write(*,*) 'ATTENTO A DISTRI, ridis'
            endif 
            call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,0)
          endif
        else
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,kkk)
        endif 
        tn(1) = 0.d0
        yn(1) = ygn(1) 
        zn(1) = zgn(1) 
        write(55,*) '# jt ' ,jt,nn,kkk 
          write(55,'(i4,3d15.7)') 1,yn(1),zn(1), 0.
        do i = 1,nn
          ip     = i+1
          dtn    = 0.5d0*(rft(i)+rft(i+1))*tl  
          tn(ip) = t0 + dtn
          call splint(tn(ip),y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip) = y 
          zn(ip) = z 
          write(55,'(i4,3d15.7)') ip,y,z,rft(i)
        enddo
        write(55,*)
        r1 = tcb(1)
        r0 = tg(ngo)
        rl = r1-r0
        a1 = ramii*amii/rl
        a2 = 1.d0*amii/rl
        et1 = eskkk
        et2 = eskkk
        n1 = nn/2
        n2 = nn-n1
        if(kk.eq.0) then
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn1,kkk)
          if(nn1.ne.nn1old)then
            nn1  =  nn1old
            n2   =  nn1 - n1
            if(n2.lt.0)then 
              write(*,*) 'ATTENTO A DISTRI 2, ridis'
            endif 
            call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn1,0)
          endif
        else
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn1,kkk)
        endif 
        write(55,*) ' # ciao' ,a1,a2,nn1,kkk
        do i = 1,nn1
          ip     = nn+1+i
          dtn    = 0.5d0*(rft(i)+rft(i+1))*rl  
          tn(ip) = r0 + dtn
          call splint(tn(ip),y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip) = y 
          zn(ip) = z 
          write(55,'(i4,3d15.7)') ip,y,z,rft(i)
        enddo
        write(55,*) 
        write(55,*) 
        ninc   = nn + nn1 + 1  + ng*(1-kget) - npc
        nn1old = nn1 
        nnold  = nn 
*
      elseif(amdi-4.d0*amii.gt.tg(ngo).or.kmed.eq.1)then
*  - non voglio che cominci con kk=0
  98    continue
        if(kmed.eq.0.and.kk.eq.0) goto 99
        write(*,*) 'PIPIPI'
        kmed  = 1
        amdi  = tcb(1) - amii
        amdi1 = amdi-tg(ngo)
        nn1   = int( log( 1.d0+(eskk-1.d0)*amdi1/amii)/log(eskk) )
        if(kk.eq.0.and.kmed.eq.1) nn1 = nn1old
        if(nn1.ne.0) then
           amii1 = amdi1*(1.d0-eskk)/(1.d0-eskk**nn1)
          write(*,*) 'PIPIPI nn1 = ',nn1
        else
          write(*,*) 'PIPIPI nn1 = 0'
        endif
        nn1old = nn1
        if(nn1.ne.0)  nn1  = nn1+1
*
        amii0  = amii1*eskk**(nn1-1) 
        amdi0  = tg(ngo)-0.5d0*amii0
        nn     = int( log( 1.d0+(eskk-1.d0)*amdi0/amii0)/log(eskk) )
        if(kk.eq.0) nn = nnold
        amii0  = amdi0*(1.d0-eskk)/(1.d0-eskk**nn)
        ninc   = nn + nn1 + ng*(1-kget) - npc
        nnold  = nn
*
        dr           = amii
        r            = amdi + 0.5d0*dr 
        tn(nn+nn1+1) = r
        call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        yn(nn+nn1+1) = y 
        zn(nn+nn1+1) = z
        ygb(1)       = y
        zgb(1)       = z
        tgb(1)       = r
        r = r - dr
        dr = amii1 
        write(*,*) 'RID ', amdi1,amii1,nn1
        write(*,*) 'RID ', amdi0,amii0,nn  
        write(28,'(2i4,3d15.7)') jt,nn+nn1+1,r, y,z
        do ip = nn+nn1,nn+1,-1
          tn(ip)   = r 
          call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip)   = y 
          zn(ip)   = z
          ksep(ip) = 1
          if(ip.eq.nn+1) ksep(ip)=2
          r        = r - dr
          dr       = dr*eskk
          write(28,'(2i4,3d15.7)') jt,ip,r,y,z
        enddo
          write(28,'(2i4,3d15.7)') 
*
        r             = amdi0 
        dr            = amii0 
        do ip = nn,1,-1
          r        = r-dr
          if(ip.eq.1) r = 0.d0
          tn(ip)   = r 
          call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip)   = y 
          zn(ip)   = z
          ksep(ip) = 0
          dr       = dr*eskk
          write(28,'(2i4,3d15.7)') jt,ip,r,y,z
        enddo
          write(28,*) 
          write(28,*) 
*
      else
 
 99     write(*,*) 'PUPUPU'
*
********************     
* - !! eskk =! 1.d0  !!!
      nn = int( log( 1.d0+(eskk-1.d0)*amdi/amii)/log(eskk) )
      no = npc - ng*(1-kget)
* - calcolo differenza 
      ninc = nn+ng*(1-kget)-npc
c - non cambio se kk=0
      if(kk.eq.0.and.ng.gt.0) then
        ninc = 0
        nn   = no
      endif
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
      nnold = nn
* - ricalcolo amii
      amii = amdi*(1.d0-eskk)/(1.d0-eskk**nn)
      write(*,*) 'amii amdi2 ',amii,amdi,nn
* - griglio bulk
      r        = amdi 
      tn(nn+1) = r
      call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
      yn(nn+1) = y 
      zn(nn+1) = z
      ygb(1)   = y
      zgb(1)   = z
      tgb(1)   = r
      dr       = amii
ccc      write(27,'(2i4,3d15.7)') jt,nn+1,r, y,z
      do ip = nn,1,-1
        r        = r-dr
        if(ip.eq.1) r = 0.d0
        tn(ip)   = r 
        call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        yn(ip)   = y 
        zn(ip)   = z
        dr       = dr*eskk
ccc        write(27,'(2i4,3d15.7)') jt,ip,r,y,z
      enddo
*
      endif
**************************
ccc        write(27,*)
      npco = npc
      npc  = nn  + nn1 +  ng*(1-kget)
      npt  = npt + ninc  
* - aggiungo getto a yn zn, se ng>0 e kget = 0
      if(ng.gt.0)then
        do i=2,ng+1
          yn(nn+nn1+i) = ygb(i)  
          zn(nn+nn1+i) = zgb(i)
          tn(nn+nn1+i) = tgb(i)
ccc        write(27,'(2i4,3d15.7)') jt,npc+i,tgb(i), ygb(i),zgb(i)
        enddo
      endif
ccc      write(27,*) 
ccc      write(27,*)
*
      if(kk.eq.0)then
*
c - aggiorno solo centroidi del corpo, se kk=0
      tt = 0.d0
      do i = 1,npc+nng
       if(ksep(i).ne.1)then
        tt = 0.5d0*(tn(i)+tn(i+1)) 
        call splint(tt,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        ycn(i)  = y 
        zcn(i)  = z
ccc        write(25,'(2i4,2d15.7)') jt,i,y,z 
       endif
      enddo
ccc      write(25,*) 
ccc      write(25,*)
*
      else
*
c --  ascissa curvilinea
      tpo(npco+nngo+1) = tn(npc+nng+1)
c      npo    = npco+nngo+1
      npo    = npco+nngo
      dy     = ycn(1)-yn(1)
      dz     = zcn(1)-zn(1)
      tt     = sqrt(dy**2+dz**2) 
      tpo(1) = tt
      phpo(1)= phin(1) 
      write(*,*) 'ppppp0 ',npco,nngo  
ccc      write(24,*) '# jt',jt,nsep
      do i=2,npco+nngo
        dy = ycn(i)-ycn(i-1)
        dz = zcn(i)-zcn(i-1)
        dd     = sqrt(dy**2+dz**2)
        tt     = tt+dd
        tpo(i) = tt
        phpo(i)= phin(i)
ccc        write(24,'(i4,4d15.7)') i,ycn(i),zcn(i),tpo(i),phpo(i)
      enddo
ccc      write(24,*)
ccc      write(24,*)
* -- estrapolo potenziale
      write(*,*) 'ppppp0 '  
      dtt  = tpo(npco+nngo)-tpo(npco+nngo-1)
      write(*,*) dtt,tpo(npco+nngo+1),tpo(npco+nngo)
      phia = phin(npco+nngo)
      phib = phin(npco+nngo-1)
      write(*,*) 'ppppp1 '  
      write(*,*) phia,phib
      delphi = (phia-phib)/dtt 
      write(*,*) 'ppppp2 ' ,ng 
      write(*,*) phia,delphi,tpo(npco+nngo+1),tpo(npco+nngo)
      phie = phia+delphi*(tpo(npco+nngo+1)-tpo(npco+nngo))
      phpo(npco+nngo+1) = phie
      write(*,*) npco+nngo+1,tpo(npco+nngo+1),phie
cc INTERPOLAZ SPLINE POTENZIALKE>>
      yp1 = 1.d31
      ypn = 1.d31
      call spline1(phpo,phpo2,tpo,yp1,ypn,npo,npamx)
c -- aggiorno divisione corpo SL, e tutti i centroidi ,se kk=1
      tt  = 0.d0
c      tysl = 1d31
c      if(ng.eq.0)then
cc      ycn(npc+nng)  = ycn(npco+nngo)
cc      zcn(npc+nng)  = zcn(npco+nngo)
cc      phin(npc+nng) = phin(npco+nngo)
c      else
c      ycn(npc+nng)  = ycb(nng)
c      zcn(npc+nng)  = zcb(nng)
c      endif
c      tt = tcb(nng)
c      call splint1(tt,ph,phpo,phpo2,tpo,npo,npamx)
c      phin(npc+nng) = ph
ccc      write(23,*) '# jt',jt,nsep
      m  = int(frint*ng)
      n  = ng - m
      if(nsep.eq.1) ksep(npc+nng) = 1
      tto = 0.d0
      do i = 1,npc+nng
        if(i.le.npc)then
        tt = 0.5d0*(tn(i)+tn(i+1))
        else
        tt = tcb(i-npc)
        endif 
        call splint(tt,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        if(nsep.eq.1)then
        epst = (tt-tto)/3.d0
        if(tt.ge.tg(ngo)-epst)then
          if(ksep(i-1).eq.0)then
            ksep(i) = 2
            kord(i) = 2
            ichin   = i
c          elseif(y.lt.tysl)then
c            ksep(i) = 2
          elseif(i-npc.le.m+0.25d0*n)then
            ksep(i) = 1
            kord(i) = 2
          else
            ksep(i) = 1
            kord(i) = 1
          endif 
        else
          ksep(i) = 0
          kord(i) = 2
        endif
        endif
c        endif
        ycn(i)  = y 
        zcn(i)  = z
        tc(i)   = tt
        tto     = tt
        if(i.gt.1)then
        call splint1(tt,ph,phpo,phpo2,tpo,npo,npamx)
        phin(i) = ph 
        endif
ccc        write(23,'(i4,5d15.7,i4)') i,tt,tg(ngo),y,z,ph,ksep(i)
      enddo
ccc      write(23,'(i4,5d15.7,i4)') npc+nng,tpo(npo),tg(ngo),ycn(npc+nng),
ccc     #      zcn(npc+nng),phin(npc+nng),ksep(npc+nng)
ccc      write(23,*)
ccc      write(23,*)
* - controllo sulla separazione
      if(nsep.eq.0)then
        if(ycn(npc+nng-ne-1).gt.ygn(ngo))then
           write(*,*) 'RIDIS SEPARA .. ',nsep,npc+nng
           write(*,*) ne,ygn(ngo)
           write(*,*) ycn(npc+nng-ne-1),zcn(npc+nng-ne-1) 
           write(*,*) ycn(npc+nng-ne),zcn(npc+nng-ne) 
           write(*,*) ycn(npc+nng-ne+1),zcn(npc+nng-ne+1) 
           write(*,*) ycn(npc+nng),zcn(npc+nng) 
           nsep = 1
c           ng   = ng-ne
c           nng  = nng-ne
c           yn(npc+nng+1) = yn(npc+nng+1+ne)
c           zn(npc+nng+1) = zn(npc+nng+1+ne)
c           ksep(npc+ng)   = 2
c           ksep(npc+ng-1) = 2
c           npsl = npsl - ne
c           do i = 1,npsl
c             ycnsl(i)  = ycnsl(i+ne)
c             zcnsl(i)  = zcnsl(i+ne)
c             ynsl(i+1) = ynsl(i+1+ne)
c             znsl(i+1) = znsl(i+1+ne)
c             phinsl(i) = phinsl(i+ne)
c           enddo
            do i=npc+ng-ne,npc+ng
              ksep(i) = 1
              kord(i) = 1
            enddo
            ksep(npc+ng-ne-1) = 2
            kord(npc+ng-ne-1) = 2
c            kord(npc+ng-ne-1) = 1
* - lunghezza e angolo del vertice 
           ay = ynsl(1)-ycn(npc+ng)      
           az = znsl(1)-zcn(npc+ng)
           aa = sqrt(ay**2+az**2)      
           by = ycnsl(1)-ycn(npc+ng)      
           bz = zcnsl(1)-zcn(npc+ng) 
           bb = sqrt(by**2+bz**2)      
           sc = (ay*by + az*bz)/(aa*bb)
           th = acos(sc)
           di = aa
           ang= th
           write(*,*) 'e OORA ' ,di,ang*180./pi
        endif
      endif
*
      endif     
*
      return
      end
