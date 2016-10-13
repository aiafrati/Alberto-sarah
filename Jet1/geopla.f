* -------------------------------------------------------
      subroutine geopla(x0,dzm,alpha,yg,zg,rnxg,ng,al,bl,tl)
* -------------------------------------------------------
*      implicit double precision (a-h,o-z)
      include "slam_p.h"
      parameter (nmax=npamx,nsecm=300,kmax=nsecm)
      dimension yg(nmax),zg(nmax),rnxg(nmax)
      dimension rnyg(nmax),rnzg(nmax)
      dimension x(4*nmax),z(4*nmax),xi(4*nmax),ze(4*nmax)
      dimension xb(4*nmax),zb(4*nmax),xib(4*nmax),zeb(4*nmax)
      dimension xi1(nsecm,nmax),et1(nsecm,nmax),ze1(nsecm,nmax)
      dimension t1(nsecm,nmax)
      dimension t1l(nsecm),np1(nsecm)
      write(*,*) '.......... :  geopla' 
* --- 
      eps = 5.d-6
      pi    = acos(-1.d0)
      ca    = cos(alpha)
      sa    = sin(alpha)
      tga   = tan(alpha)
* --- linea di chiglia e insellatura fino a sezione  dopo x0 
      open(unit=7, file='geo.in')
      rewind(7)
      read(7,*) nsec
      do 10, i=1,nsec
        read(7,*) np
        read(7,*) xi(i),ynull,ze(i)
        x(i) =  xi(i)*ca + ze(i)*sa
        z(i) = -xi(i)*sa + ze(i)*ca + dzm
        if(x(i).gt.x0) goto 11
        do 20, j=2,np
          read(7,*)
  20    continue
  10  continue
  11  rewind(7) 
      isec = i
*
      read(7,*)
      do 15, i=1,nsec
        read(7,*) np
        do 25, j=1,np-1
          read(7,*)
  25    continue
        read(7,*) xib(i),ynull,zeb(i)
        xb(i) =  xib(i)*ca + zeb(i)*sa
        zb(i) = -xib(i)*sa + zeb(i)*ca + dzm
        if(xb(i).gt.x0) goto 16
  15  continue
  16  rewind(7)
      isecb = i
* -- DIAGNOSTICA --
      write(*,*) 'isec, isecb, nsec ',isec,isecb,nsec
      if(isecb.gt.isec)then
        isecb = isec
        write(*,*)  'geopla   : OKKIO, isecb > isec '
        write(*,*)  '           provo con isecb = isec'
      endif
      if(isec.eq.1)then
        dx = abs(x(1)-x0)
        if(dx.le.eps)then
          x0   = x(1) + 0.5d0*eps
          isec = 2
          isecb= 2
          write(*,*) 'geopla   : OKKIO, x0 < x(1)!'
          write(*,*) '           di poco, lo sposto dentro'
        else
          write(*,*) 'geopla   : ALT, x0 < x(1) !!!!'
          STOP
        endif
      endif
      if(isec.gt.nsec)then
        dx = abs(x(nsec)-x0)
        if(dx.le.eps)then
          x0   = x(nsec) - 0.5d0*eps
          isec = nsec
          if(isecb.gt.nsec) isecb = nsec
          write(*,*) 'geopla   : OKKIO, x0 > x(nsec)!'
          write(*,*) '           di poco, lo sposto dentro'
        else
          write(*,*) 'geopla   : ALT, x0 > x(nsec) !!!!'
          STOP
        endif
      endif
*       
* --- vado a prendere le sezioni tra isecb-1 e isec
      read(7,*)
      do 30, i=1,isecb-2
        read(7,*) np
        do 40,j=1,np
           read(7,*)
  40    continue
  30  continue
      do 35, i=1,isec-isecb+1+1
      read(7,*) np1(i)
      read(7,*) xi1(i,1),et1(i,1),ze1(i,1)
*      write(*,*) xi1(i,1),et1(i,1),ze1(i,1)
      t1(i,1) = 0.d0   
        do 50, j=2,np1(i)
           read(7,*) xi1(i,j),et1(i,j),ze1(i,j)   
*           write(*,*) xi1(i,j),et1(i,j),ze1(i,j)   
           det = et1(i,j) - et1(i,j-1)
           dze = ze1(i,j) - ze1(i,j-1)
           t1(i,j) = t1(i,j-1) + sqrt(det**2 + dze**2)
  50    continue
      t1l(i) = t1(i,np1(i))
  35  continue
* --- calcolo z0 come intersezione tra linea di chiglia e x=x0
      w1x  = 0.d0
      w1z  = 1.d0
      x2   = x(isec-1)
      z2   = z(isec-1)
      w2x  = x(isec)-x(isec-1)
      w2z  = z(isec)-z(isec-1)
      call intret(xin,zin,r,s,x0,0.d0,w1x,w1z,x2,z2,w2x,w2z)
      z0   = z2 + s*w2z
      xi0  = x0*ca - (z0-dzm)*sa
      zi0  = x0*sa + (z0-dzm)*ca
* --- calcolo coordinate sezione 
      wxso  = w2x
      wyso  = 0.d0 
      wzso  = w2z
      yg(1) = 0.d0
      zg(1) = z0     
*      ng = 21
* --------- ESATTA -----------------------------------
*          zs     = -xis*sa + (zi+z0)*ca
        fwx       = fwigx(xi0,zi0,al,bl,tl)
        fwz       = fwigz(xi0,zi0,al,bl,tl)
        rn        = sqrt(1.d0 + fwx**2 + fwz**2)
*        rnxg(1)   = -(fwx*ca + fwz*sa)/rn
* --------- ESATTA END --------------------------------       
      tsa = 0.d0 
      dtsa= 1.d0/float(ng-1)
      ii  = isec-isecb+1+1 
      do 70, j=2,ng
        tsa = tsa + dtsa
* --- interpolo su sezione 1 e trovo pto (xi1s..) di ascissa tsa
        do 75, k=1,kmax 
*        write(*,*) 'ciao 3 ',k,ii
        ts1 = tsa*t1l(ii)
        do 80, j1=1,np1(ii)
          if(t1(ii,j1).gt.ts1) goto 88
  80    continue
  88    if(j1.gt.np1(ii)) j1 = np1(ii) 
*        write(*,*) j1,t1(ii,j1),ts1
        w1   = (ts1-t1(ii,j1-1))/(t1(ii,j1)-t1(ii,j1-1))
        xi1s =  xi1(ii,j1-1)*(1.d0 - w1) + xi1(ii,j1)*w1
        et1s =  et1(ii,j1-1)*(1.d0 - w1) + et1(ii,j1)*w1
        ze1s =  ze1(ii,j1-1)*(1.d0 - w1) + ze1(ii,j1)*w1
        t1xi =  xi1(ii,j1) - xi1(ii,j1-1)
        t1et =  et1(ii,j1) - et1(ii,j1-1)
        t1ze =  ze1(ii,j1) - ze1(ii,j1-1)
*        write(*,*) w1,xi1s,et1s,ze1s
* --- interpolo su sezione 2 e trovo pto (xi2s..) di ascissa tsa 
        ts2 = tsa*t1l(ii-1)
        do 90, j2=1,np1(ii-1)
          if(t1(ii-1,j2).gt.ts2) goto 99
  90    continue
  99    if(j2.gt.np1(ii-1)) j2 = np1(ii-1) 
        w2   = (ts2-t1(ii-1,j2-1))/(t1(ii-1,j2)-t1(ii-1,j2-1))
*        write(*,'(a3,i3,3d15.5)') 'j2 .',j2,t1(ii-1,j2),ts2,w2
        xi2s =  xi1(ii-1,j2-1)*(1.d0 - w2) + xi1(ii-1,j2)*w2
        et2s =  et1(ii-1,j2-1)*(1.d0 - w2) + et1(ii-1,j2)*w2
        ze2s =  ze1(ii-1,j2-1)*(1.d0 - w2) + ze1(ii-1,j2)*w2
        t2xi =  xi1(ii-1,j2) - xi1(ii-1,j2-1)
        t2et =  et1(ii-1,j2) - et1(ii-1,j2-1)
        t2ze =  ze1(ii-1,j2) - ze1(ii-1,j2-1)
* --- retta da (xi1s,..) a (xi2s,...)  (=tau1)
        wxs = (xi1s-xi2s)*ca + (ze1s-ze2s)*sa        
        wys = (et1s-et2s)
        wzs =-(xi1s-xi2s)*sa + (ze1s-ze2s)*ca        
        x2s = xi2s*ca + ze2s*sa
        y2s = et2s             
        z2s =-xi2s*sa + ze2s*ca + dzm
* --- interseco con x=x0
        call intret(xin,zin,r,s,x0,z0,w1x,w1z,x2s,z2s,wxs,wzs)
        if(s.lt.0.d0.and.ii.gt.2)then
          write(*,*) 'prendo sezione prima, s=',s
* qui va aggiunta la parte di ricalcolo sezioni 1 e 2
          ii = ii-1
        else
*          write(*,*) 'sezione OK, s=',s
          goto 77
        endif
  75    continue
  77    continue
        xs = x2s + wxs*s
        ys = y2s + wys*s
        zs = z2s + wzs*s
        yg(j) = ys
        zg(j) = zs
* --------- ESATTA -----------------------------------
        xis = xs*ca - (zs-dzm)*sa
        ets = ys
        zes = xs*sa + (zs-dzm)*ca
*          xis    = x0/ca - (zi+z0)*tga
        yge  = fwig(xis,zes,al,bl,tl)
        zge  = -xs*tga + zes*sa*tga + zes*ca + dzm
*          zs     = -xis*sa + (zi+z0)*ca
        fwx       = fwigx(xis,zes,al,bl,tl)
        fwz       = fwigz(xis,zes,al,bl,tl)
        rn        = sqrt(1.d0 + fwx**2 + fwz**2)
        rnxge     = -(fwx*ca + fwz*sa)/rn
*        rnxg(j)   = rnxge 
* --------- ESATTA END --------------------------------       
        tauxn = 0.d0
        if(j.gt.2)then
          tauyn = yg(j)-yg(j-1)
          tauzn = zg(j)-zg(j-1)
          tauy = 0.5d0*(tauyn+tauyo)
          tauz = 0.5d0*(tauzn+tauzo)
        else
          tauyn = yg(j)-yg(j-1)
          tauzn = zg(j)-zg(j-1)
          tauy = tauyn
          tauz = tauzn
        endif  
        taux = 0.d0 
*        write(*,'(i4,5d15.5)') j-1,tauy,tauz,wxso,wyso,wzso
* -- A: normale come media tra calcolo con tauo e taun.
        rnx = tauy*wzso - tauz*wyso
        rny = tauz*wxso - taux*wzso 
        rnz = taux*wyso - tauy*wxso
        rn  = sqrt(rnx**2 + rny**2 + rnz**2) 
        rnxg(j-1) = rnx/rn
        rnxga     = rnx/rn
        rnyg(j-1) = rny/rn
        rnzg(j-1) = rnz/rn
* -- B: altro modo: media tra wxs e wxso, e assegno rn al centroide!! 
        wxsm = 0.5d0*(wxso+wxs)
        wysm = 0.5d0*(wyso+wys)
        wzsm = 0.5d0*(wzso+wzs)
        rnx = tauyn*wzsm - tauzn*wysm
        rny = tauzn*wxsm - tauxn*wzsm 
        rnz = tauxn*wysm - tauyn*wxsm
        rn  = sqrt(rnx**2 + rny**2 + rnz**2)
        rnmx = rnx/rn
        rnmy = rny/rn
        rnmz = rnz/rn
*        rnxg(j-1) = rnx/rn
*        rnyg(j-1) = rny/rn
*        rnzg(j-1) = rnz/rn
* -- C: altro modo : media tra normali calcolate con t1 e t2.
*                    manca da fare la normale per j=1
        t1x = t1xi*ca + t1ze*sa 
        t1y = t1et
        t1z = -t1xi*sa + t1ze*ca
        t2x = t2xi*ca + t2ze*sa 
        t2y = t2et
        t2z = -t2xi*sa + t2ze*ca
        rn1x = t1y*wzs - t1z*wys
        rn1y = t1z*wxs - t1x*wzs 
        rn1z = t1x*wys - t1y*wxs
        rn2x = t2y*wzs - t2z*wys
        rn2y = t2z*wxs - t2x*wzs 
        rn2z = t2x*wys - t2y*wxs
        rn1 = sqrt(rn1x**2 + rn1y**2 + rn1z**2) 
        rn2 = sqrt(rn2x**2 + rn2y**2 + rn2z**2)
        rn1x = rn1x/rn1 
        rn1y = rn1y/rn1 
        rn1z = rn1z/rn1 
        rn2x = rn2x/rn2 
        rn2y = rn2y/rn2 
        rn2z = rn2z/rn2
*        rnxg(j) = 0.5d0*(rn1x+rn2x)
*        rnxg(j) = s*rn1x + (1.d0-s)*rn2x
*        rnyg(j) = s*rn1y + (1.d0-s)*rn2y
*        rnzg(j) = s*rn1z + (1.d0-s)*rn2z
        wxso = wxs 
        wyso = wys 
        wzso = wzs 
        tauxo = 0.d0 
        tauyo = tauyn
        tauzo = tauzn
*
*        write(41,'(8d15.7)') x0,xs,yg(j),zg(j),yge,zge,rnxga,rnxge
  70  continue     
*      write(41,*) 
*      write(41,*) 
      rnx = tauyn*wzs - tauzn*wys
      rny = tauzn*wxs - tauxn*wzs 
      rnz = tauxn*wys - tauyn*wxs 
      rn = sqrt(rnx**2 + rny**2 + rnz**2) 
      rnxg(ng) = rnx/rn 
      rnyg(ng) = rny/rn 
      rnzg(ng) = rnz/rn 
*
      close(7)
      return
      end
