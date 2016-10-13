
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine initial(t,amplim,vfall,llf,spo0,pro0,npc,npt,npsl,
     #             kget,ng,yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,amp,ampsl,
     #             rny,rnz,rnysl,rnzsl,tmy,tmz,tmysl,tmzsl,ampp,vfall0,
     #             dphi,phisl,pfraz,estr,escr,epsg,epsgg,ampli,
     #             ygn,zgn,ygs2,zgs2,tg,ngo,ngo1,nsep,nsepo,ksep,kord,
     #         nnold,nn1old,frin,frfi,jind,tin,jfid,tfi,jend,kmed,ksup)
      include"slam_p.h"
      include"slam_f.h"
      dimension yv(npamx),zv(npamx),yce(npamx),zce(npamx)
      dimension ysl(npamx),zsl(npamx),amp(npamx),ampsl(npamx)
      dimension rny(npamx),rnz(npamx),tmy(npamx),tmz(npamx)
      dimension rnysl(npamx),rnzsl(npamx),tmysl(npamx),tmzsl(npamx)
      dimension ycsl(npamx),zcsl(npamx)
      dimension dphi(npamx),phisl(npamx)
      dimension yg(npamx),zg(npamx),ygn(npamx),zgn(npamx)
      dimension tg(npamx),ygs2(npamx),zgs2(npamx),ksep(npamx)
      dimension kord(npamx)
      write(*,*) '--> initial...........'
*
c     - inizializzazione variabili

      pi    = acos(-1.d0)
      call filiniz
      write(*,*) '--> PIPPOOOOO.........'
*
c       - inizializzazione variabili per inteagrazione temporale
*
        t      = 0.d0 
        amplim = ampp/pfraz
        vfall  = vfall0
        llf    = 0
*
* - lettura file offset
* --- lettura geometria corpo. il vertice e' a y=0, z=0
*
        write(*,*) 'pippo' 
        open(unit=7,file='geo.in',status='unknown')
        read(7,*) ngo
        do ig=1,ngo
          read(7,*) yg(ig), zg(ig)
        enddo
        close(7)
        write(*,*) 'pippo2' 
*
* -- variabili del getto separato
*
        kmed = 0 
        ksup = 0 
        ngo1=ngo
        nsep=0
        nsepo=0
        do i=1,npamx
          ksep(i)=0
c          kord(i)=1
          kord(i)=2
        enddo
*
*  --- interpolazione spline della geometria iniziale (-> ygs2,zgs2)
*
* ---   calcolo ascissa curvilinea ta
*
        tg(1) = 0.d0
        do 20, ip=2,ngo
           dy = yg(ip) - yg(ip-1)
           dz = zg(ip) - zg(ip-1)
           dl = sqrt(dy**2 + dz**2)
           tg(ip) = tg(ip-1) + dl
  20    continue
*
        write(*,*) 'pippo3'
        yp1 = 1.d31
        ypn = 1.d31
        call spline(yg,zg,ygs2,zgs2,tg,yp1,ypn,ngo,npamx)
        write(*,*) 'pippo4'
*
        ngon1 = 20
        ngon2 = 200
        dtg1 = 0.75d0*tg(ngo)/float(ngon1)
        dtg2 = 0.25d0*tg(ngo)/float(ngon2)
        write(*,*) 'pippo5',dtg1,dtg2
        ygn(1)=yg(1)
        zgn(1)=zg(1)
        tai1 = 0.d0
        do i=2,ngon1+1
          tai1 = tai1+dtg1
          if(i.eq.ngon1+1) tai1=0.75d0*tg(ngo)
          call splint(tai1,ygn(i),zgn(i),yg,zg,ygs2,zgs2,tg,ngo,npamx,0)
ccc          write(55,*) ygn(i),zgn(i),tai1
        enddo 
        tai2 = 0.75d0*tg(ngo)
        write(*,*) 'pippo6'
        do i=ngon1+2,ngon1+ngon2+1
          tai2 = tai2+dtg2
          if(i.eq.ngon1+ngon2+1) tai2=tg(ngo)
          call splint(tai2,ygn(i),zgn(i),yg,zg,ygs2,zgs2,tg,ngo,npamx,0)
ccc          write(55,*) ygn(i),zgn(i),tai2
        enddo 
ccc          write(55,*) 
ccc          write(55,*) 
        write(*,*) 'pippo7'
        ngo  =ngon1+ngon2+1
        ngo1 =ngo
ccc          write(55,*) 1,yg(1),zg(1)
        do ip=2,ngo
          yg(ip) = ygn(ip)
          zg(ip) = zgn(ip)
          dy = yg(ip) - yg(ip-1)
          dz = zg(ip) - zg(ip-1)
          dl = sqrt(dy**2 + dz**2)
ccc          write(55,*),ip,yg(ip),zg(ip)
          tg(ip) = tg(ip-1) + dl
        enddo
        call spline(yg,zg,ygs2,zgs2,tg,yp1,ypn,ngo,npamx)
        write(*,*) 'pippo8'
*              
* - taglio la geometria a z = -pro0
c        pro0 = -6.25d-3
        y1  = 0.d0
        z1  = -pro0
        w1y = 1.d0
        w1z = 0.d0
        write(*,*) '-pro0 ' ,-pro0 
        do i=2,ngo
          y2  = yg(i-1)
          z2  = zg(i-1)
          w2y = yg(i)-yg(i-1)
          w2z = zg(i)-zg(i-1)
          call intret(yi,zi,si,ti,y1,z1,w1y,w1z,y2,z2,w2y,w2z)
          if(ti.ge.0.d0.and.ti.le.1.d0) goto 999
          write(*,'(4d15.6)') si,ti,yg(i-1),zg(i-1)
        enddo
 999    kint = i-1
        dyi  = yi-yg(kint)
        dzi  = zi-zg(kint)
        tai  = tg(kint) + sqrt( dyi**2 + dzi**2 )
        write(*,*) 'kint ',kint,yg(kint),zg(kint)
*
* -- sposto ygn a pro0
*
        do i=1,ngo
          zgn(i) = zg(i)+pro0
          ygn(i) = yg(i)
        end do
*
**
* -- genero la griglia sulla parte tagliata della geometria e la traslo
*      di pro0
*      la posizione  r = float(iv)/float(npci+1)*tai e' come 
*      in Code (il primo
*      pannello viene lungo doppio cosi)
*
        npci = 6 
        r    = 0.d0
        dr   = tai/float(npci)
        do 30, iv = 1,npci+1
*           r = r + dr
           r = float(iv)/float(npci+1)*tai
           if (iv.eq.1)      r = 0.d0
           if (iv.eq.npci+1) r = tai
           call splint(r,yy,zz,ygn,zgn,ygs2,zgs2,tg,ngo,npamx,0)
           yv(iv)   = yy
           zv(iv)   = zz 
c           write(40,*) iv,yv(iv),zv(iv)
  30    continue        
        npc = npci
        nnold  = npc
        nn1old = 0
*
c       - discretizzazione iniziale della superficie libera
*
        spo0 = yv(npc+1)
        npr  = log(1.d0+(estr-spo0)/ampp*(escr-1.d0))/log(escr) + 1
        ampl = (estr-spo0)*(1.d0-escr)/(1.d0-escr**npr)
        ysl(1) = yv(npc+1)
        zsl(1) = zv(npc+1)
c        write(41,*) 1,ysl(1),zsl(1)
        do i = 1,npr
          ysl(i+1) = ysl(i) + ampl
          zsl(i+1) = 0.d0
          ampl    = ampl*escr
c        write(41,*) i+1,ysl(i+1),zsl(i+1)
        enddo
        npsl = npr
        npt  = npc + npsl 
* - getto
        ng   = 0
        kget = 0
        ampli = 0.d0
        epsg = epsgg*ampli 

*
        if (npt.gt.npamx) stop 'aumenta npamx'

c       - inizializzo potenziale di perturbazione e ampiezze dei 
c         pannelli sulla sup. libera

        do i = 1,npc
          yce(i) = (yv(i+1)+yv(i))/2.d0
          zce(i) = (zv(i+1)+zv(i))/2.d0
          amy    =  yv(i+1) - yv(i)
          amz    =  zv(i+1) - zv(i)
          am     = sqrt(amy*amy+amz*amz)
          amp(i) = am
          tmy(i) = amy/am     
          tmz(i) = amz/am     
          rny(i) =  tmz(i)
          rnz(i) = -tmy(i)
          dphi(i)= -vfall*rnz(i)
        enddo
        do i = 1,npsl  
          ycsl(i) = (ysl(i+1)+ysl(i))/2.d0
          zcsl(i) = (zsl(i+1)+zsl(i))/2.d0
          amy     =  ysl(i+1) - ysl(i)
          amz     =  zsl(i+1) - zsl(i)
          am      = sqrt(amy*amy+amz*amz)
          ampsl(i)= am
          tmysl(i) = amy/am     
          tmzsl(i) = amz/am     
          rnysl(i) =  tmzsl(i)
          rnzsl(i) = -tmysl(i)
          phisl(i) = 0.d0
        enddo
        cmm0 = 0.d0 
        do i = 1,npc
          cmm0 = cmm0 + (yv(i+1)-yv(i))*(zv(i+1)+zv(i))/2.d0
        enddo
        do i = 1,npsl
          cmm0 = cmm0 + (ysl(i+1)-ysl(i))*(zsl(i+1)+zsl(i))/2.d0
        enddo
        write(*,*) cmm0
*
*
* ---  striscia di SL da monitorare
*
        frin  = 0.002d0
        frfi  = 0.005d0
        jend  = 0
*
        yin = frin*estr
        yfi = frfi*estr
        do j=1,npsl
          if(yin.lt.ycsl(j))goto 98
        enddo
 98     jind = j-1
        do j=1,npsl
          if(yfi.lt.ycsl(j))goto 99
        enddo
 99     jfid = j-1
        tin  = (yin-ycsl(jind))/(ycsl(jind+1)-ycsl(jind))
        tfi  = (yfi-ycsl(jfid))/(ycsl(jfid+1)-ycsl(jfid))
*
*

*
      return
      end
