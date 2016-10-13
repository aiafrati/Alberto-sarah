      subroutine shallo(ynsl,znsl,ygb,zgb,kget,jt,npsl,
     #               ng,ampli,rmg,iiget,jjget,
     #               ycnsl,zcnsl,phinsl,dphisl,epsg,epsgg,eskg,gfrac,
     #               ygn,zgn,ygs2,zgs2,tg,tgb,ngo,ngo1,nsep,
     #               ycb,zcb,tcb,nngo,jind,jfid,jend)
*
      include "slam_p.h"
      dimension ynsl(npamx),znsl(npamx)
      dimension ygb(npamx),zgb(npamx),ycb(npamx),zcb(npamx)
      dimension ycnsl(npamx),zcnsl(npamx),phinsl(npamx),dphisl(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx) 
      dimension tg(npamx),tgb(npamx),tcb(npamx) 
      write(*,*) '--> shallo......'
*
* --  vedo se buttare pannelli troppo vicini
      if(nsep.eq.0)then
      epsg = epsgg*ampli
      else
c      epsg = max(epsgg*ampli,0.0d0*ygn(ngo))
      epsg = 3.d0*epsgg*ampli
      endif
      epss = 1.d-10 
c      eskg = 1.08d0
*
      nngo = ng*kget
c      nsep = 1
c      write(37,*) '#  jt ',jt
c      do j=1,ngo1
c        write(37,*) ygn(j),zgn(j)
c      enddo
c      write(37,*)
c      write(37,*)
      if(nsep.lt.1)then
*
      write(*,*) 'epsg nsep 0',epsg,epsgg,ampli,nsep
      kk   = 0
      kbi  = 0
        do i=1,ng
          y1  = ycnsl(i)
          z1  = zcnsl(i)
          do j=1,ngo1-1
            y2   = ygn(j)
            z2   = zgn(j)
            w2y  = ygn(j+1)-ygn(j)
            w2z  = zgn(j+1)-zgn(j)
            ww   = sqrt(w2y**2+w2z**2)
            w1y  = -w2z/ww
            w1z  =  w2y/ww
            call intret(yi,zi,si,ti,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
            if(ti.ge.0.d0.and.ti.le.1.d0)then
              if(abs(si).lt.epsg) kk=i
            endif
          enddo
        enddo
      if(kk.gt.0)then
        write(*,*) 'shallo 0, butto ',kk,' pannelli'
        do i=kk+1,npsl
          ycnsl(i-kk) = ycnsl(i)
          zcnsl(i-kk) = zcnsl(i)
          phinsl(i-kk)= phinsl(i)
          dphisl(i-kk)= dphisl(i) 
          ynsl(i-kk)  = ynsl(i) 
          znsl(i-kk)  = znsl(i)
        enddo
        ynsl(npsl+1-kk)  = ynsl(npsl+1) 
        znsl(npsl+1-kk)  = znsl(npsl+1)
* --  shift indici striscia di controllo
        if(jind.lt.kk+1.and.jend.eq.0)then
           write(*,*)'attenzione, la striscia se ne va!!!'
           write(*,*) 'jend shallo = 1'
           jend=1
        elseif(jend.eq.0)then 
          jind=jind-kk
          jfid=jfid-kk
        endif
      endif
      npsl = npsl - kk
*
c      if(nsep.lt.1)then
* - punto di intersezione corpo SL
      y1 = ynsl(1)
      z1 = znsl(1)
      kkin = 9999
      do j=1,ngo1-1
        y2   = ygn(j)
        z2   = zgn(j)
        w2y  = ygn(j+1)-ygn(j)
        w2z  = zgn(j+1)-zgn(j)
        ww   = sqrt(w2y**2+w2z**2)
        w1y  = -w2z/ww
        w1z  =  w2y/ww
        call intret(yi,zi,si,ti,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
        if(ti.ge.(0.d0-epss).and.ti.le.(1.d0+epss))then
          kb = j
          yb = yi 
          zb = zi
          ay = yi-ygn(kb)
          az = zi-zgn(kb)
          tb = tg(kb)+sqrt(ay**2+az**2)
          kkin = 0
          goto 999
        endif
      enddo
 999  continue
c      if(kkin.gt.0)then
c        write(*,*) 'shallo, no intersezione punta' 
c        yb = ygn(ngo1)
c        zb = zgn(ngo1)
c        tb = tg(ngo1)
c        kb = ngo1-1
c      endif
      
*
      else
*
      write(*,*) 'epsg nsep 1',epsg,epsgg,ampli,nsep
      kk   = 0
      kbi  = 0
        do i=1,ng
          y1  = ycnsl(i)
          z1  = zcnsl(i)
          do j=1,ngo1-1
            y2   = ygn(j)
            z2   = zgn(j)
            w2y  = ygn(j+1)-ygn(j)
            w2z  = zgn(j+1)-zgn(j)
            ww   = sqrt(w2y**2+w2z**2)
            w1y  = -w2z/ww
            w1z  =  w2y/ww
            call intret(yi,zi,si,ti,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
            if(ti.ge.0.d0.and.ti.le.1.d0)then
              if(abs(si).lt.epsg) kk=i
            endif
          enddo
        enddo
      if(kk.gt.0)then
        write(*,*) 'shallo 1, butto ',kk,' pannelli'
        do i=kk+1,npsl
          ycnsl(i-kk) = ycnsl(i)
          zcnsl(i-kk) = zcnsl(i)
          phinsl(i-kk)= phinsl(i)
          dphisl(i-kk)= dphisl(i) 
          ynsl(i-kk)  = ynsl(i) 
          znsl(i-kk)  = znsl(i)
        enddo
        ynsl(npsl+1-kk)  = ynsl(npsl+1) 
        znsl(npsl+1-kk)  = znsl(npsl+1)
      endif
      npsl = npsl - kk
      ngo1 = ngo1 - kk
* - punto di intersezione corpo SL
      y1 = ynsl(1)
      z1 = znsl(1)
      kkin = 9999
      j = ngo1-1
      y2   = ygn(j)
      z2   = zgn(j)
      w2y  = ygn(j+1)-ygn(j)
      w2z  = zgn(j+1)-zgn(j)
      ww   = sqrt(w2y**2+w2z**2)
      w1y  = -w2z/ww
      w1z  =  w2y/ww
      call intret(yi,zi,si,ti,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
      kb = ngo1-1 
      yb = yi 
      zb = zi
      ay = yi-ygn(kb)
      az = zi-zgn(kb)
      tb = tg(kb)+sqrt(ay**2+az**2)
      tg(ngo1)  = tb
      ygn(ngo1) = yb
      zgn(ngo1) = zb
      kkin = 0
      call spline(ygn,zgn,ygs2,zgs2,tg,yp1,ypn,ngo1,npamx)
      write(*,*) 'punti usat i',j,ngo,ti
      write(*,*) ygn(j),zgn(j)
      write(*,*) ygn(j+1),zgn(j+1)
      write(*,*) ynsl(1),znsl(1)
      write(*,*) yb,zb 
******************
c        if(kk.gt.0)then
c        ng   = ng - kk
c        ngo1 = ngo1 - kk
c        tb = 0.5d0*(tg(ngo1)+tg(ngo1-1))
c        call splint(tb,yb,zb,ygn,zgn,ygs2,zgs2,tg,ngo1+kk,npamx,0)
c        ygn(ngo1) = yb 
c        zgn(ngo1) = zb
c        tg(ngo1)  = tb
c        yp1 = 1.d31
c        ypn = 1.d31
c        call spline(ygn,zgn,ygs2,zgs2,tg,yp1,ypn,ngo1,npamx)
c        endif
***********
ccc        tb = 0.5d0*(tg(ngo1-1) + tg(ngo1-2))
ccc        yb = 0.5d0*(ygn(ngo1-1) + ygn(ngo1-2))
ccc        zb = 0.5d0*(zgn(ngo1-1) + zgn(ngo1-2))
c        call splint(tb,yb,zb,ygn,zgn,ygs2,zgs2,tg,ngo1+kk,npamx,0)
***********
      endif
ccc      write(38,*) '#  jt ',jt
ccc      do j=1,ngo1
ccc        write(38,'(3d15.7,i4)') ygn(j),zgn(j),tg(j),j
ccc      enddo
ccc      write(38,*)
ccc      write(38,*)

* - calcolo della zona getto
        if(jt.gt.iiget)then
          ii = 2 
          dyo =1.d0
          do i=2,npsl
            ii=ii+1
            dy = ynsl(i+1)-ynsl(i)
            if (dy.gt.0.d0.and.dyo.le.0.d0) ng=ii
            dyo = dy
          enddo 
          write(*,*) 'ng! ',ng,ngo1,ngo
          if(ng.gt.0.and.jt.gt.jjget) kget = 1 
c - separo getto (shallow water) dal resto
c - ampiezza (lungo il corpo) della zona getto
          if(ng.gt.0)then
            y1  = ynsl(ng+1)
            z1  = znsl(ng+1)
            do j=1,ngo1-1
              y2   = ygn(j)
              z2   = zgn(j)
              w2y  = ygn(j+1)-ygn(j)
              w2z  = zgn(j+1)-zgn(j)
              w1y  = -w2z
              w1z  =  w2y
              call intret(yi,zi,si,ti,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
              if(ti.ge.0.d0.and.ti.le.1.d0) then 
                kbb = j 
                ybb = yi 
                zbb = zi
                ay  = yi-ygn(kbb)
                az  = zi-zgn(kbb)
                tbb = tg(kbb)+sqrt(ay**2+az**2)
                write(*,*) 'BIGO.. ',kbb,tbb,tg(kbb),tb,kb,tg(kb)
                goto 9
              endif
            enddo
  9         continue
            ampg = gfrac*(tb-tbb)
            tbi  = tb-ampg
            call splint(tbi,ybi,zbi,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
            ygb(1) = ybi
            zgb(1) = zbi
            tgb(1) = tbi
* - intersezione con la SL (per determinare ampli e quindi ng)
            do j=1,ngo1
              if(tg(j).gt.tbi)then
                kbi = j-1
                goto 88
              endif
            enddo
  88        continue
            if(kbi.eq.0) then
              write(*,*) 'intersezione non trovata SJALLO'
              write(*,*) 'kbi,',kbi,tbi,ampg,tb
              stop
            endif
            wby = -(zgn(kbi+1)-zgn(kbi))           
            wbz =   ygn(kbi+1)-ygn(kbi)            
            do i=1,ng
              y1  = ynsl(i)
              z1  = znsl(i)
              w1y = ynsl(i+1)-ynsl(i) 
              w1z = znsl(i+1)-znsl(i) 
              call intret(yi,zi,si,ti,y1,z1,w1y,w1z,ybi,zbi,wby,wbz)
              if(si.ge.0.d0.and.si.le.1.d0) goto 888  
            enddo
  888       continue       
            ki = i
            ay = ybi-yi
            az = zbi-zi
            write(*,*) 'ng2-1! ',ng,ampli
            amplio = ampli
            ampli = sqrt(ay**2+az**2)/rmg
            write(*,*) 'ng2-0! ',ng,ampli
            ampli = 0.5*ampli + 0.5*amplio
*- calcolo ng
            if(eskg.ne.1.d0)then
              ng = int( log( 1.d0+(eskg-1.d0)*ampg/ampli)/log(eskg) )
              if(ng.le.1) ng = 2
              ampli = ampg*(1.d0-eskg)/(1.d0-eskg**ng)
            else
              ng   = int(ampg/ampli)
              if(ng.le.1) ng = 2
              ampg = float(ng)*ampli 
            endif
            write(*,*) 'ng2! ',ng,ampli
c - array pti corpo nel getto, vertici e centroidi
            dr = ampli
            r  = tbi
            write(29,*) '# jt ',jt,nsep
            write(28,*) '# jt ',jt,nsep
            write(29,'(i4,3d15.7)') 1,ygb(1),zgb(1),tgb(1)
            do i = 1,ng
              r  = r+dr
              rc = r-0.5d0*dr
c              if(i.eq.ng) r = tb
              call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
              call splint(rc,yy,zz,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
              ygb(i+1) = y
              zgb(i+1) = z
              tgb(i+1) = r
              ycb(i)   = yy
              zcb(i)   = zz
              tcb(i)   = rc
              dr  = dr*eskg
            write(29,'(i4,3d15.7)') i+1,ygb(i+1),zgb(i+1),tgb(i+1)
            write(28,'(i4,3d15.7)') i,ycb(i),zcb(i),tcb(i)
            enddo
            if(nsep.eq.1)then
             ygb(ng+1) = ygn(ngo1)
             zgb(ng+1) = zgn(ngo1)
             tgb(ng+1) = tg(ngo1) 
             ycb(ng) = ygn(ngo1-1)
             zcb(ng) = zgn(ngo1-1)
             tcb(ng) = tg(ngo1-1)
            endif 
            write(28,*)
            write(29,*)
            write(29,*)
            write(28,*)
          endif
        endif
*
        return
        end 
