
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine disun2(jt,yn,zn,ng,ysl,zsl,ycsl,zcsl,phisl,npsl,
     #             ygb,zgb,escr,kget,estr,amplim,ampli,npt,npc,iiget,
     #             ycb,zcb,jind,jfid,tin,tfi,jend)
*
*
      include"slam_p.h"
      dimension xx(npamx),zz(npamx),pp(npamx)
      dimension xs2(npamx),zs2(npamx),ps2(npamx)
      dimension tp(npamx),tv(npamx)
c
      dimension ysl(npamx),zsl(npamx),ycsl(npamx),zcsl(npamx) 
      dimension phisl(npamx),ygb(npamx),zgb(npamx)
      dimension ygfc(npamx),zgfc(npamx),ygf(npamx),zgf(npamx)
      dimension yn(npamx),zn(npamx) ,ycb(npamx),zcb(npamx) 
*
      write(*,*) '--> disun2........',iiget,jt
        write(*,*) ' ng '  ,ng
*       
      kin = 0
      kfi = 0
      tcc  = 0.d0
      tcco = 0.d0
*
        amy =  yn(npc+1) - yn(npc) 
        amz =  zn(npc+1) - zn(npc)
        am  = sqrt(amy**2 + amz**2)
        rny =  amz/am
        rnz =  -amy/am 
        a2y = ysl(2)-ysl(1)
        a2z = zsl(2)-zsl(1)
        amc = (ysl(2)-ysl(1))*rny +  (zsl(2)-zsl(1))*rnz 
        amt = sqrt((ysl(2)-ysl(1))**2 + (zsl(2)-zsl(1))**2 )
        sca = (-amy*a2y-amz*a2z)/(amt*am)
        anint = acos(sca)*180.d0/pi
        write(*,*) 'am1 4 ',am1,anint
        anint=90.d0
        if(anint.lt.90.d0) then
          am1=amc
        else
          am1 = amt
        endif
        am1 = max(am1,amplim)
        write(*,*) 'am1 5 ',am1
        if(jt.gt.iiget)      am1 = ampli 
        write(*,*) 'am1 55 ',am1
       
c     controllo l'angolo massimo tra due pannelli:
c     se e' maggiore di un valore limite riduco l'ampiezza dei pannelli
c     di conseguenza
        if(jt.le.iiget)then
        amma = 0.d0
        do i = 1,npsl/2
          ss = (ysl(i+1)-ysl(i))*(ysl(i+2)-ysl(i+1)) +
     &         (zsl(i+1)-zsl(i))*(zsl(i+2)-zsl(i+1))
          ammp = sqrt( (ysl(i+1)-ysl(i))**2 + (zsl(i+1)-zsl(i))**2 )
          amms =sqrt((ysl(i+2)-ysl(i+1))**2+(zsl(i+2)-zsl(i+1))**2 )
          sc = ss/(ammp*amms)
          if (abs(sc).gt.1.d0) sc = sc/abs(sc)
          an = abs(acos(sc)*180.d0/pi)
          amma = max(amma,an)
        enddo
        if (amma.gt.10.d0) am1 = am1*10.d0/amma
        endif
          write(*,*) 'am1 3 ',am1
c      endif
*
c     - inizializzazione dati per interpolazione spline
*
      npa   = npsl+1
      tp(1) = 0.d0
      xx(1) = ysl(1)
      zz(1) = zsl(1)
      ama   = sqrt( (ysl(2)-ysl(1))**2 + (zsl(2)-zsl(1))**2 )
      amb   = sqrt( (ysl(3)-ysl(2))**2 + (zsl(3)-zsl(2))**2 )
c - estrapolo phi nel pto di intersezione, lato SL
      pp(1) = phisl(1)+(phisl(1)-phisl(2))/(ama+amb)*ama
*
      do i = 1,npa-1
        if (i.eq.1) then
          ax      = (ycsl(i)-ysl(i))
          az      = (zcsl(i)-zsl(i))
        else
          ax      = (ycsl(i)-ycsl(i-1))
          az      = (zcsl(i)-zcsl(i-1))
        end if
        dis     = sqrt(ax*ax+az*az)
        tp(i+1) = tp(i) + dis
        xx(i+1) = ycsl(i)
        zz(i+1) = zcsl(i)
        pp(i+1) = phisl(i)
      enddo
c      ysl(npa+1) = estr
      ax      = (ysl(npsl+1)-ycsl(npsl))
      az      = (zsl(npsl+1)-zcsl(npsl))
      dis     = sqrt(ax*ax+az*az)
      tp(npa+1) = tp(npa) + dis
      xx(npa+1) = ysl(npsl+1)
      zz(npa+1) = zsl(npsl+1)
      pp(npa+1) = phisl(npsl) ! sarebbe meglio estrapolazione lineare
*
* --  assegno i punti della striscia di controllo
*
      iind = jind + 1 
      ifid = jfid + 1
      ttin = tp(iind) + tin*(tp(iind+1)-tp(iind))
      ttfi = tp(ifid) + tfi*(tp(ifid+1)-tp(ifid))
      write(*,*) 'jjjind ' ,ttin, ttfi, iind, ifid
*
c - interpolazione spline di y,z e potenziale
      yp1 = 1.d+31
      ypn = 1.d+31
      call spline3(xx,zz,pp,xs2,zs2,ps2,tp,yp1,ypn,npa+1,npamx)
      if(ng.gt.0)then
c - cerco intersezione normale corpo  con SL
c  ricostruzione posizione vertici  getto
      kk  = 1
      nn  = ng
      ttlo= 0.d0
      write(36,*) '# jt ' ,jt 
      write(37,*) '# jt ' ,jt 
      write(37,'(i4,4d15.7)')
     #              j,ygf(1),zgfc(1),ygb(ng+1),zgb(ng+1)
      do   j = 1,ng
        yk  = ygb(ng+1-j)
        zk  = zgb(ng+1-j)
        wky = -(zgb(ng+2-j)-zgb(ng+1-j))
        wkz =   ygb(ng+2-j)-ygb(ng+1-j)
        do ii = 1,npsl-1
          y1 = xx(ii)
          z1 = zz(ii)
          w1y = xx(ii+1)-xx(ii)
          w1z = zz(ii+1)-zz(ii)
          call intret(yin,zin,si,ti,yk,zk,wky,wkz,y1,z1,w1y,w1z)
          if(ti.gt.0.d0.and.ti.le.1.d0) goto 999
        enddo
  999   kk  = ii
        ay  = yin - xx(kk) 
        az  = zin - zz(kk)
        ttl = tp(kk)+sqrt(ay**2+az**2) 
* ---- interpolo linearmente nel getto
        if(j.le.4)then
        ygf2 = yin
        zgf2 = zin
c        yk = 0.5d0*(ygb(ng+1-j)+ygb(ng+2-j))
c        zk = 0.5d0*(zgb(ng+1-j)+zgb(ng+2-j))
        yk = ycb(ng+1-j)
        zk = zcb(ng+1-j)
        do ii = 1,npsl-1
          y1 = xx(ii)
          z1 = zz(ii)
          w1y = xx(ii+1)-xx(ii)
          w1z = zz(ii+1)-zz(ii)
          call intret(yiin,ziin,si,ti,yk,zk,wky,wkz,y1,z1,w1y,w1z)
          if(ti.gt.0.d0.and.ti.le.1.d0) goto 99
        enddo
  99    kkc  = ii
        if(j.le.4)then
        write(*,*) 'DISDIS ', j,kkc,kk,ti
        endif
        ygfc2   = yiin
        zgfc2   = ziin
        phisl2  = pp(ii)*(1.d0-ti) + ti*pp(ii+1)
        ygf(j+1)= ygf2
        zgf(j+1)= zgf2
        ygfc(j) = yiin
        zgfc(j) = ziin
        phisl(j)= phisl2
        write(37,'(i4,4d15.7)')
     #              j,ygf(j+1),zgfc(j+1),ygb(ng+1-j),zgb(ng+1-j)
        write(36,'(i4,4d15.7)')
     #               j,ygfc(j),zgfc(j),ycb(ng+1-j),zcb(ng+1-j)
        else
*-----   interpolaz spline nel getto     
        call splint(ttl,xpo,zpo,xx,zz,xs2,zs2,tp,npa+1,npamx,0)
        if(j.le.nn)then
c          ss=float(j)/float(nn)
          ss=1.d0
        else
          ss=1.d0
        endif
        ygf(j+1) = ygf2*(1.d0-ss)   + ss*xpo
        zgf(j+1) = zgf2*(1.d0-ss)   + ss*zpo
        write(37,'(i4,4d15.7)')
     #              j,ygf(j+1),zgfc(j+1),ygb(ng+1-j),zgb(ng+1-j)
        tcc = ttlo+(ttl-ttlo)*0.5d0
        call splint3(tcc,xo,zo,po,xx,zz,pp,xs2,zs2,ps2,tp,npa+1,npamx,0)
        ygfc(j)  = ygfc2*(1.d0-ss)  + ss*xo
        zgfc(j)  = zgfc2*(1.d0-ss)  + ss*zo
        phisl(j) = phisl2*(1.d0-ss) + ss*po
        write(36,'(i4,4d15.7)')
     #               j,ygfc(j),zgfc(j),ycb(ng+1-j),zcb(ng+1-j)
        endif
*------------------*        
        ttlo = ttl
      enddo
       write(37,*) 
       write(37,*) 
       write(36,*) 
       write(36,*) 
        ay  = yin - xx(kk) 
        az  = zin - zz(kk)
        ttl = tp(kk)+sqrt(ay**2+az**2) 
c -riposiziono nell'array ysl e ycsl
      tcc =0.d0
      do j = 1,ng
        ysl(j+1) = ygf(j+1)
        zsl(j+1) = zgf(j+1)
        ycsl(j) = ygfc(j)
        zcsl(j) = zgfc(j)
        if(j.eq.1)then
          ay = ycsl(j)-ysl(1)
          az = zcsl(j)-zsl(1)
        else
          ay = ycsl(j)-ycsl(j-1)
          az = zcsl(j)-zcsl(j-1)
        endif
        tcc = tcc+sqrt(ay**2+az**2)
* --  riposiziono la striscia di controllo
*
        if(ttin.lt.tcc.and.kin.eq.0)then
          jind = j-1 
          tin  = (ttin-tcco)/(tcc-tcco)
          kin  = 1
          if(jind.eq.0) then 
            jend=1
            write(*,*) 'disun2, striscia sparisce!'
            write(*,*) 'jend dis 1 =1'
          endif
        endif
        if(ttfi.lt.tcc.and.kfi.eq.0)then
          jfid = j-1 
          tfi  = (ttfi-tcco)/(tcc-tcco)
          kfi  = 1
        endif
*
        write(78,'(3i5,3d15.7)') jt,j,jind,ttin,tcc,tcco
        tcco = tcc
      enddo
       write(78,*) 
*
      endif
c  ricostruzione posizione vertici bulk
c      if(kget.eq.0)then
        amme = (tp(2)+tp(3))/2.d0
        write(*,*) 'amme   ',amme
        write(*,*) 'amplim ',amplim
        write(*,*) 'am1    ',am1
 
        if (am1.gt.amme) then
          if (amme.lt.amplim) then
            det = min(am1,amplim)
          else
            det = amme
          end if
        else
          det = am1
        end if
        write(*,*) 'det ',det
        if(ng.eq.0) then
           det=det/escr
c           det=det
           ttl = 0.d0
        else
           det = det
        endif
        if(jt.gt.iiget)then
          escr2 = 1.d0
          det = am1
        endif
        ivp = 1 
c      elseif(kget.eq.1)then
c        det = ampli/escr
c        ivp = 1 
c      endif
      write(*,*) 'ttl '
      write(*,*) ttl,det,npa,tp(npa+1),escr
      kt = 1
      do while(ttl.lt.tp(npa+1)) 
c       - vedo quanto manca ed in base a questo fisso il det

        dtma = tp(npa+1) - ttl 
        if (dtma.gt.2.d0*escr*det) then
          if(jt.gt.iiget.and.ttl.lt.tp(ng+1))then
           det =  det*escr2
          else
           det  = det*escr
          endif
        else if (dtma.gt.escr*det) then
          det  = dtma/2.d0
        else
          det  = dtma
        end if 
*
        ivp = ivp+1
        tcc = ttl + det/2.d0
        ttt = ttl + det
* - iterpol. spline centroidi e potenziale 
        call splint3(tcc,xo,zo,po,xx,zz,pp,xs2,zs2,ps2,tp,npa+1,npamx,0)
        ycsl(ng+ivp-1) = xo
        zcsl(ng+ivp-1) = zo
        phisl(ng+ivp-1) = po
* - interpol lineare centroidi
c        do j = kt,npa
c          if(tp(j).gt.tcc) then
c            kc = j-1
c            goto 88
c          endif
c        enddo
c  88    continue
c        am = sqrt( (xx(kc+1)-xx(kc))**2 + (zz(kc+1)-zz(kc))**2)
c        wx = xx(kc+1)-xx(kc)
c        wz = zz(kc+1)-zz(kc)
c        wp = pp(kc+1)-pp(kc)
c        s  = (tcc-tp(kc))/am
c        ycsl(ng+ivp-1)  = xx(kc) + s*wx 
c        zcsl(ng+ivp-1)  = zz(kc) + s*wz 
c        phisl(ng+ivp-1) = pp(kc) + s*wp
c        kt = kc 
* - iterpol. spline vertici
c        write(*,*) 'disun splint', ng,ivp,ttt 
        call splint(ttt,xpo,zpo,xx,zz,xs2,zs2,tp,npa+1,npamx,0)
c        write(*,*) 'disun splint' ,ng,ivp,xpo
        ysl(ng+ivp) = xpo
        zsl(ng+ivp) = zpo
* - interpol lineare vertici
c        if(ttt.eq.tp(npa+1))then
c          ysl(ng+ivp) = xx(npa+1)  
c          zsl(ng+ivp) = zz(npa+1) 
c        else 
c          do j = kt,npa
c            if(tp(j).gt.ttt) then
c              kv = j-1
c              goto 888
c            endif
c          enddo
c  888     continue
c          am = sqrt( (xx(kv+1)-xx(kv))**2+(zz(kv+1)-zz(kv))**2)
c          wx = xx(kv+1)-xx(kv)
c          wz = zz(kv+1)-zz(kv)
c          s  = (ttt-tp(kv))/am
c          ysl(ng+ivp) = xx(kv) + s*wx 
c          zsl(ng+ivp) = zz(kv) + s*wz 
c          kt = kv
c        endif
c       - aggiorno ttl e riprendo
* --  riposiziono la striscia di controllo
*
        if(ttin.lt.tcc.and.kin.eq.0)then
          jind = ng+ivp-2
          tin  = (ttin-tcco)/(tcc-tcco)
          kin  = 1
          if(jind.eq.0) then
            jend = 1
            write(*,*) 'disun, la striscia va nel pizzo!'
            write(*,*) 'jend dis 2 = 1'
          endif
        endif
        if(ttfi.lt.tcc.and.kfi.eq.0)then
          jfid = ng+ivp-2
          tfi  = (ttfi-tcco)/(tcc-tcco)
          kfi  = 1
        endif
        write(78,'(3i5,3d15.7)') jt,ng+ivp-1,jind,ttin,tcc,tcco
*
        tcco = tcc

        ttl = ttt
      enddo
      nng = ng*kget
      npsl= ng + ivp-1
      npt = npsl - ng*kget + npc 
*
      return
      end      
