
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine disuni(phisu,anint,jind,jfid,tin,tfi)

      include"slam_p.h"
      include"slam_v.h"
      dimension xx(npamx+1),zz(npamx+1),pp(npamx+1)
      dimension xs2(npamx+1),zs2(npamx+1),ps2(npamx+1)
      dimension tp(npamx+1),tv(npamx+1)
*
      write(*,*) '.......... :  disuni'
*
* - calcolo max(abs(dphi))
c      dphim=0.d0
c      do i = npc+kcut+1,npt-1
c         dphim=max( dphim,abs(dphi(i)) )
c     end do
*- vorrei redistribuire i vertici tenendo conto di un limite tra
*  la differenza (stimate) sul valore di dphi su due centrodi 
*  successivi  
*
c     Ridistribuzione dei punti di controllo.
*
c     - definizione della ascissa curvilinea
*
      kfi = 0
      kin = 0
*   
      write(*,*) '  aniinit ', anint 
      if (anint.lt.90.d0.and.kcut.eq.0) then
        amy =  yn(npc+1) - yn(npc) 
        amz =  zn(npc+1) - zn(npc)
        am  = sqrt(amy**2 + amz**2)
        rny =  amz/am
        rnz =  -amy/am 
        am1 = (yn(npc+2)-yn(npc+1))*rny + 
     &        (zn(npc+2)-zn(npc+1))*rnz
        write(*,*) 'am1 4 ',am1
        am1b = (yn(npc+2)-yn(npc+1))*sin(ande) -
     &        (zn(npc+2)-zn(npc+1))*cos(ande)
        am1 = max(am1,amplim)
        write(*,*) 'am1 5 ',am1
      else
        ams = sqrt((yn(npc+2+kcut)-yn(npc+1+kcut))**2 +
     &             (zn(npc+2+kcut)-zn(npc+1+kcut))**2  )
        amt = sqrt((yn(npc+2)-yn(npc+1))**2 +
     &             (zn(npc+2)-zn(npc+1))**2  )
        am1 = min(amt,ams)
*        write(*,*) 'am1 1 ', am1
        am1 = max(am1,amplim)
*        write(*,*) 'am1 2 ', am1
      end if
*
c     controllo l'angolo massimo tra due pannelli:
c     se e' maggiore di un valore limite riduco l'ampiezza dei pannelli
c     di conseguenza
*
      amma = 0.d0
      do i = npc+kcut+1,(npt+npc)/2
        ss = (yn(i+1)-yn(i))*(yn(i+2)-yn(i+1)) +
     &       (zn(i+1)-zn(i))*(zn(i+2)-zn(i+1))
        ammp = sqrt( (yn(i+1)-yn(i))**2 + (zn(i+1)-zn(i))**2 )
        amms = sqrt( (yn(i+2)-yn(i+1))**2 + (zn(i+2)-zn(i+1))**2 )
        sc = ss/(ammp*amms)
        if (abs(sc).gt.1.d0) sc = sc/abs(sc)
        an = abs(acos(sc)*180.d0/pi)
        amma = max(amma,an)
      end do
      if (amma.gt.10.d0) am1 = am1*10.d0/amma
*        write(*,*) 'am1 3 ', am1
*
c     - inizializzazione dati per interpolazione spline
*
      npa   = npt-npc-kcut+1
      tp(1) = 0.d0
      xx(1) = yn(npc+kcut+1)
      zz(1) = zn(npc+kcut+1)
      ihat  = npc+kcut
*
      ama   = sqrt( (yn(ihat+2)-yn(ihat+1))**2 +
     &              (zn(ihat+2)-zn(ihat+1))**2 )
      amb   = sqrt( (yn(ihat+3)-yn(ihat+2))**2 +
     &              (zn(ihat+3)-zn(ihat+2))**2 )
      pp(1) = phin(ihat+1)+(phin(ihat+1)-phin(ihat+2))/
     &        (ama+amb)*ama
*
      do i = 1,npa-1
        if (i.eq.1) then
          ax      = (ycn(i+npc+kcut)-yn(i+npc+kcut))
          az      = (zcn(i+npc+kcut)-zn(i+npc+kcut))
        else
          ax      = (ycn(i+npc+kcut)-ycn(i+npc+kcut-1))
          az      = (zcn(i+npc+kcut)-zcn(i+npc+kcut-1))
        end if
        dis     = sqrt(ax*ax+az*az)
        tp(i+1) = tp(i) + dis
        xx(i+1) = ycn(i+npc+kcut)
        zz(i+1) = zcn(i+npc+kcut)
        pp(i+1) = phin(i+npc+kcut)
      end do
      ax      = (yn(npt+1)-ycn(npt))
      az      = (zn(npt+1)-zcn(npt))
      dis     = sqrt(ax*ax+az*az)
      tp(npa+1) = tp(npa) + dis
      xx(npa+1) = yn(npt+1)
      zz(npa+1) = zn(npt+1)
      pp(npa+1) = phin(npt) ! sarebbe meglio estrapolazione lineare
*
* --  assegno i punti della striscia di controllo
*
      iind = jind - (npc+kcut-1)
      ifid = jfid - (npc+kcut-1)
      ttin = tp(iind) + tin*(tp(iind+1)-tp(iind))
      ttfi = tp(ifid) + tfi*(tp(ifid+1)-tp(ifid))
*
c     interpolazione spline di y,z e potenziale
*
      yp1 = 1.d+31
      ypn = 1.d+31
      call splone(xs2,tp,xx,yp1,ypn,npa+1,npamx)
      call splone(zs2,tp,zz,yp1,ypn,npa+1,npamx)
      call splone(ps2,tp,pp,yp1,ypn,npa+1,npamx)
*
c     ricostruzione posizione vertici
*
      amme = (tp(2)+tp(3))/2.d0
      if (am1.gt.amme) then
        if (amme.lt.amplim) then
          ttl = min(am1,amplim)
        else
          ttl = amme
        end if
      else
        ttl = am1
      end if
      write(*,*) 'det ',ttl
*
c     - assegno il primo vertice sulla sup. lib.
*
      call splont(ttl,xpo,tp,xx,xs2,npa+1,npamx,xnull,0,xnull,0,xnull,0)
      call splont(ttl,zpo,tp,zz,zs2,npa+1,npamx,xnull,0,xnull,0,xnull,0)
*
      yn(npc+kcut+2) = xpo
      zn(npc+kcut+2) = zpo
      ivp   = 2
      det   = ttl
*
c     - assegno il centroide e potenziale
*
      pen = ( phin(npc+kcut+1)-phin(npc+kcut+2) )/
     &      ( tp(2)-tp(3) )
      trn = phin(npc+kcut+1)-pen*tp(2)
      tcc = ttl/2.d0
      call splont(tcc,xo,tp,xx,xs2,npa+1,npamx,xnull,0,xnull,0,xnull,0)
      call splont(tcc,zo,tp,zz,zs2,npa+1,npamx,xnull,0,xnull,0,xnull,0)
      call splont(tcc,po,tp,pp,ps2,npa+1,npamx,xnull,0,xnull,0,xnull,0)
*
      ycn(npc+kcut+1) = xo
      zcn(npc+kcut+1) = zo
      phi(npc+kcut+1) = pen*tcc + trn ! lineare nel primo tratto
c      phi(npc+kcut+1) = po
*
      do while (ttl.lt.tp(npa+1))
*
c       - vedo quanto manca ed in base a questo fisso il det
*
        dtma = tp(npa+1) - ttl 
        if (dtma.gt.2.d0*escr*det) then
          det  = det*escr
        else if (dtma.gt.escr*det) then
          det  = dtma/2.d0
        else
          det  = dtma
        end if 
*
c       - assegno il vertice
*
        ttt = ttl + det
        call splont(ttt,xpo,tp,xx,xs2,npa+1,npamx,
     #              xnull,0,xnull,0,xnull,0)
        call splont(ttt,zpo,tp,zz,zs2,npa+1,npamx,
     #              xnull,0,xnull,0,xnull,0)
*
        ivp = ivp + 1
        yn(npc+kcut+ivp) = xpo
        zn(npc+kcut+ivp) = zpo
*
c       - assegno il centroide e potenziale
*
        tcc = ttl + det/2.d0
        call splont(tcc,xo,tp,xx,xs2,npa+1,npamx,
     #              xnull,0,xnull,0,xnull,0)
        call splont(tcc,zo,tp,zz,zs2,npa+1,npamx,
     #              xnull,0,xnull,0,xnull,0)
        call splont(tcc,po,tp,pp,ps2,npa+1,npamx,
     #              xnull,0,xnull,0,xnull,0)
*
        ycn(npc+kcut+ivp-1) = xo
        zcn(npc+kcut+ivp-1) = zo
        phi(npc+kcut+ivp-1) = po
*
* --  riposiziono la striscia di controllo
*
        if(ttin.lt.tcc.and.kin.eq.0)then
          jind = npc+kcut+ivp-2
          tin  = (ttin-tcco)/(tcc-tcco)
          kin  = 1
        endif
        if(ttfi.lt.tcc.and.kfi.eq.0)then
          jfid = npc+kcut+ivp-2
          tfi  = (ttfi-tcco)/(tcc-tcco)
          kfi  = 1
        endif
*
*
*
c       - aggiorno ttl e riprendo
*
        ttl = ttt
        tcco= tcc
      end do
*
      npt = npc+kcut+ivp-1
*       
      do i = npc+kcut+1,npt
        phin(i) = phi(i)
      end do
*
    
      return
      end
