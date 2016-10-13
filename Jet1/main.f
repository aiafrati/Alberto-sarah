*
      program main
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
c     ++                                                 ++
c     +  programma di calcolo del flusso generato dalla   +
c     +    caduta di un                                   +
c     + corpo simmetrico 2d a assialsimmetrico            +
c     +                                                   +
c     ++                                                 ++
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

c     - File inclusi: 

c       -- "slam_p.h" in cui sono riportati la definizione 
c                     delle variabili e i parametri per il 
c                     dimensionamento,
c       -- "slam_v.h" nel quale sono riportate le aree common 
c                     e sono dimensionati i vettori.

      include"slam_p.h"
      include"slam_v.h"
      include"slam_f.h"

      dimension taa(npamx),yv2(npamx),zv2(npamx)
      dimension rrny(npamx),rrnz(npamx),ttmy(npamx),ttmz(npamx)
      dimension accnn(npamx)
      dimension aamy(npamx),aamz(npamx)
      dimension dpht(npamx),dph2(npamx) 
      dimension tg(ngmax),wg(ngmax) 
*
      dimension pres(npamx)
      dimension fosl2(ntmx),yco(npamx),zco(npamx),tto(npamx)
      dimension ycoo(npamx),zcoo(npamx),rrnyo(npamx),rrnzo(npamx)
      dimension rrnyoo(npamx),rrnzoo(npamx)
      dimension ttmzo(npamx),ttmyo(npamx),ampo(npamx)
      dimension phio(npamx),dphio(npamx),dphto(npamx)
      dimension vyo(npamx),vzo(npamx),yvo(npamx),zvo(npamx)
*
      dimension yyn(npamx),zzn(npamx),pphi(npamx)
      dimension ddphi(npamx),ddpht(npamx),yycn(npamx),zzcn(npamx)
      ggr  = 0.d0
      ipre = 0
      call input(i2d,kg) ! lettura dati necessari per integrazione 
*
      open(unit=30,file='mass')
      open(unit=35,file=sfor,status='NEW')
      open(unit=36,file='ppp')
      open(unit=41,file='pp1')
      open(unit=42,file='pp2')
      open(unit=43,file='pp3')
      open(unit=44,file='pp4')
      open(unit=45,file='pp5')
      open(unit=46,file='pp6')
      open(unit=47,file='pp7')
      open(unit=77,file='striscia')
*
      if(i2d.eq.0)then
        call gauleg(0.d0,1.d0,tg,wg,kg)
      endif
      if(kini.eq.1)then
        call initia   ! inizializzazione delle variabili
      else
        call initial(frin,frfi,jind,tin,jfid,tfi)
      endif
*
*
* --- inizializzo array per subroutine press
      do i = 1,npc
         vyo(i) = 0.d0 
         vzo(i) = 0.d0 
         yco(i) = 0.d0 
         zco(i) = 0.d0 
         tto(i) = 0.d0
        phio(i) = 0.d0 
      end do
      t0 = abs(pro0/vfall0)
*
c - chiamo il solutore e calcolo le velocita tangenti
*
      adim = abs(zv(1))
      adim2= abs(vfall*zv(1))
*
      call solver(jt,yv,zv,amp,phi,dphi,kphi,npt,npc,kcut,tg,wg,kg,i2d)
      call caveta(phi,amp,dpht,npc,npt,kcut,dphi)
      call nortan(1,npt,yv,zv,aamy,aamz,amp,ttmy,ttmz,
     #            rrny,rrnz,yce,zce,0)
*
      do ip = 1,npt
        vym(ip,1) = dphi(ip)*rrny(ip) + dpht(ip)*ttmy(ip)
        vzm(ip,1) = dphi(ip)*rrnz(ip) + dpht(ip)*ttmz(ip)
      end do
* --------- CALCOLO LUNGHEZZA STRISCIA DI CONTROLLO
          rl=0.d0
          do i=jind+1,jfid-1
            dl= sqrt( (yce(i+1)-yce(i))**2 + (zce(i+1)-zce(i))**2 )
            rl= rl + dl
          enddo
          dlin= (1.d0-tin)*sqrt( (yce(jind+1)-yce(jind))**2 +
     #                    (zce(jind+1)-zce(jind))**2 )
          dlfi= tfi*sqrt( (yce(jfid+1)-yce(jfid))**2 +
     #                    (zce(jfid+1)-zce(jfid))**2 )
          rl= rl + dlin + dlfi
          yin = yce(jind)+tin*(yce(jind+1)-yce(jind))
          zin = zce(jind)+tin*(zce(jind+1)-zce(jind))
          yfi = yce(jfid)+tfi*(yce(jfid+1)-yce(jfid))
          zfi = zce(jfid)+tfi*(zce(jfid+1)-zce(jfid))
          write(77,'(i4,2d15.6,2i4,4d15.6)') 0 ,0.d0,rl,jind,jfid+1,
     #               yin,zin,yfi,zfi                             
           
*-----------
*
      if(lrest.ne.1)call stampa(proat,pres,ntagl,nrid,npco,yco,zco,ampo
     #                ,jend,jind,jfid,yin,zin,yfi,zfi)         
*
      proat = zv(1)
      kfor  = 1
      sca   = (aamy(npc)*aamy(npc+1)+aamz(npc)*aamz(npc+1))/
     #                  (amp(npc+1)*amp(npc))
      anint = (pi-acos(sca))*180.d0/pi
* -- centroidi (x calcolo pressione ?? )
        do i = 1,npc
          yce(i) = (yv(i+1)+yv(i))/2.d0
          zce(i) = (zv(i+1)+zv(i))/2.d0
        end do
*
c - inizio integrazione temporale
*
      frdtt = frdt
      jt    = 0
      jjj   = 0
      acce  = 0.d0
      ekj   = 0.d0
      jend= 0 
*
      do while (t.lt.tend)
*
        write(*,*) '                        ------------------'
        write(*,*)
*
c - aggiorno il passo di integrazione temporale e sposto il corpo di
*           vfall*dt
*
        jt = jt + 1
        ntagl=1
        nrid =1
        if (jt.lt.10) then
          frdt=frdtt*jt/10
        else
          frdt = frdtt
        endif
        call check(proat,jt)
        write(*,*) '   anint = ',sngl(anint)
*
*        write(36,*) 
*        write(36,*)
*
c - spostamento (predictor) dei centroidi 
*
        write(*,*) '   ---> PREDICTOR <---'
        do ip = npc+1,npt
          deph2      = vym(ip,1)**2 + vzm(ip,1)**2
          ycn(ip)    = yce(ip) + vym(ip,1)*dt
          zcn(ip)    = zce(ip) + vzm(ip,1)*dt
          depn(ip,1) = deph2/2.d0 - ggr*zcn(ip)
          phin(ip)   = phi(ip) + depn(ip,1)*dt
        end do
*
c - ricostruisco la config. predictor dei vertici SL:splver...
c - ridefinizione dei pannelli sul corpo: ridis...
c - ampiezze , normali e tangenti...   : nortan...
*
        call splver(iint,ygn,zgn,ng,ycn,zcn,yn,zn,npc,npt,proat,
     #              ande,kcut,ampk,tag,ygs2,zgs2) 

ccc*----------------------------------------------------- CONTROLLI
      write(31,*) '# ',jt,zv(1)
      do i=1,npc+1
        write(31,'(2d15.7)') yn(i),zn(i)
      enddo
      write(31,*)  
      write(31,*)  
      write(32,*) '# ',jt
      do i=1,npc
        write(32,'(8d15.7)') ycn(i),zcn(i),phi(i),dphi(i),
     #            dpht(i),vym(i,1),vzm(i,1) ,amp(i)  
      enddo
      write(32,*)
      write(32,*)  
      write(33,*) '# ',jt
      do i=npc+1,npt+1
        write(33,'(4d15.7)') yn(i),zn(i)
      enddo
      write(33,*)  
      write(33,*)  
      write(34,*) '# ',jt
      do i=npc+1,npt
        write(34,'(8d15.7)')ycn(i),
     #        zcn(i),phin(i),dphi(i),
     #            dpht(i),vym(i,1),vzm(i,1) ,amp(i)  
      enddo
      write(34,*)  
      write(34,*)
ccc*-------------------------------------------------- FINE CONTROLLI


        call ridis(0,proat,nrid,jind,jfid)
        call nortan(1,npt,yn,zn,aamy,aamz,amp,ttmy,ttmz,
     #              rrny,rrnz,ycn,zcn,1)
*       
c - risolvo il problema nella configurazione intermedia
*
        do ip = 1,npc
            kphi(ip) = 0
            if(iwig.eq.0)then
            dphi(ip) = -vfall*rrnz(ip)
            else
            rnxm = 0.5d0*(rnxx(ip)+rnxx(ip+1))
            dphi(ip) = -ux*rnxm/sqrt(1.d0-rnxm**2)
            endif
        end do
        do ip = npc+1,npt 
            kphi(ip) = 1
        enddo
        if (kcut.eq.1) then
            kphi(npc+1)= 2
        end if
*
c        if(jt.ge.150)then
c        kcut=0
c        kphi(npc+1)=1 
c        endif
        call solver(jt,
     #        yn,zn,amp,phin,dphi,kphi,npt,npc,kcut,tg,wg,kg,i2d)
        call caveta(phin,amp,dpht,npc,npt,kcut,dphi)
c        if(jt.ge.150) kcut=1

*
c - spostamento definitivo dei centroidi
*
        write(*,*) '   ---> CORRECTOR <---'
        do ip = npc+1,npt
          vym(ip,2) = dphi(ip)*rrny(ip) + dpht(ip)*ttmy(ip)
          vzm(ip,2) = dphi(ip)*rrnz(ip) + dpht(ip)*ttmz(ip)
          deph2     = vym(ip,2)**2 + vzm(ip,2)**2
          ycn(ip)   = yce(ip)+(vym(ip,1)+vym(ip,2))*dt/2.d0
          zcn(ip)   = zce(ip)+(vzm(ip,1)+vzm(ip,2))*dt/2.d0
          depn(ip,2)= deph2/2.d0  - ggr*zcn(ip)
          phin(ip)  = phi(ip)+(depn(ip,1)+depn(ip,2))*dt/2.d0
        end do
*
c - ricostruisco la configurazione della SL             : splver... 
* - normali, tangenti....                               : nortan...
c - calcolo dell'angolo all'interse ev. taglio pannelli : taglio...
* - redistribuzione vertici e centroidi sulla SL        : disuni...
c - redistribuzione vertici sul corpo                   : ridis...
* - normali, tangenti... (non tocco ycn zcn sulla SL!)  : nortan...
*
c        write(*,*) 'ampk 0                  ',ampk
*
        call splver(iint,ygn,zgn,ng,ycn,zcn,yn,zn,npc,npt,proat,
     #              ande,kcut,ampk,tag,ygs2,zgs2)
        call nortan(npc,npt,yn,zn,aamy,aamz,amp,
     #                  ttmy,ttmz,rrny,rrnz,ycn,zcn,0)
*----------------------------------------------------------------
*               CALCOLO ENERGIA
*----------------------------------------------------------------
c        nnpc = npc
c        nnpt = npt
c        write(68,*) '#',jt
c        do i=1,nnpt
c           yyn(i)  = yn(i)
c           zzn(i)  = zn(i)
c           yycn(i) = ycn(i)
c           zzcn(i) = zcn(i)
c           pphi(i) = phin(i)
c           ddphi(i)= dphi(i)
c           write(68,*) i,ddphi(i),pphi(i)
c        enddo
c        write(68,*)
c        write(68,*)
c        call ridis2(1,proat,escr,yyn,zzn,ygn,ygs2,zgn,zgs2,ng,
c     #               amp,iint,tag,kcut,nnpt,nnpc,pphi,ddphi,
c     #               yycn,zzcn,vfall)
c        write(69,*) '#',jt
c        do i =1,nnpt
c           write(69,*) i,ddphi(i),pphi(i)
c        enddo
c        write(69,*)
c        write(69,*)
c  
c        call nortan(1,nnpt,yyn,zzn,aamy,aamz,amp,
c     #                  ttmy,ttmz,rrny,rrnz,yce,zce,1)
c        call solver(jt,
c     #        yyn,zzn,amp,pphi,ddphi,kphi,nnpt,nnpc,kcut,tg,wg,kg,i2d)
c        call caveta(pphi,amp,ddpht,nnpc,nnpt,kcut,ddphi)
c        write(79,*) '#',jt
c        do i =1,nnpt
c           write(79,*) i,ddphi(i),pphi(i)
c        enddo
c        write(79,*)
c        write(79,*)
c        eko = ekn
c        ekb = 0.d0
c        write(67,*) '#',jt
c        do jj=1,nnpt
c          write(67,*) jj,pphi(jj),ddphi(jj)
c          ekb = ekb + amp(jj)*pphi(jj)*ddphi(jj)
c        enddo
c        write(67,*)
c        write(67,*)
c        call tagli2(kcut,ttmy,ttmz,rrny,rrnz,aamy,aamz,anint,
c     #                  nnpt,nnpc,yyn,zzn,ygn,zgn,ng,amp,yycn,zzcn,
c     #                  pphi,ddphi,ddpht,ek,amplim,ancut)
c        ekn = ekb + ek
cc        write(66,'(i4,5d15.6)') jt,t,ek,ekb,ekn,eko
c        call nortan(1,npt,yn,zn,aamy,aamz,amp,
c     #                  ttmy,ttmz,rrny,rrnz,yce,zce,1)
*----------------------------------------------------------------
*               FINE         CALCOLO ENERGIA
*----------------------------------------------------------------
*
c        if(jt.lt.150) then
        call taglio(k2cut,ttmy,ttmz,rrny,rrnz,aamy,aamz,anint,ntagl,
     #              jind,jfid,jend)
c        endif
        call disuni(phisu,anint,jind,jfid,tin,tfi)
        call ridis(1,proat,nrid,jind,jfid)
        call nortan(1,npt,yn,zn,aamy,aamz,amp,
     #                  ttmy,ttmz,rrny,rrnz,yce,zce,1)
*
*
c - rimetto le coordinate in xv,zv e aggiorno potenziale e 
c   ampiezze e centroidi_SL
*       
       do ip=1,npc
         yv(ip) = yn(ip)
         zv(ip) = zn(ip)
         kphi(ip)= 0
         if(iwig.eq.0)then
          dphi(ip)= -vfall*rrnz(ip)
         else
          rnxm = 0.5d0*(rnxx(ip)+rnxx(ip+1))
          dphi(ip) = -ux*rnxm/sqrt(1.d0-rnxm**2)
         endif
       end do
       yv(npc+1) = yn(npc)
       zv(npc+1) = zn(npc)
       do ip=npc+1,npt   
         yv(ip)  = yn(ip)
         zv(ip)  = zn(ip)
         yce(ip) = ycn(ip)
         zce(ip) = zcn(ip) 
         kphi(ip)= 1
         phi(ip) = phin(ip)
       end do
       yv(npt+1)=yn(npt+1)
       zv(npt+1)=zn(npt+1)
*
        if (kcut.eq.1) then
          kphi(npc+1) = 2
        end if
*
c - filtro la superficie libera
*
        mm = mod(jt,ift)
        if (mm.eq.0) then
          call doldfil1(yv,npc+1+kcut+iford,npt-10,npt+1)
          call doldfil1(zv,npc+1+kcut+iford,npt-10,npt+1)
          call doldfil1(phi,npc+1+kcut+iford,npt-10,npt)
        end if
*
* - normali, tangenti , ampiezze (no centroidi)....
*
        call nortan(1,npt,yv,zv,aamy,aamz,amp,
     #                  ttmy,ttmz,rrny,rrnz,yce,zce,0)
*
c - ridefinizione centroidi tappo e primo pannello SL dopo tappo.
*   NB: disuni ha ridefinito tutti i nodi e centroidi, dal nodo
*   ncp+kcut+1 e dal centroide npc+kcut+1.
*   Percio' qui a rigore servirebbe definire solo il centroide tappo.
*   in effetti cosi facendo il codice va bene, solo fa piu' passi che
*   ridefinendo linearmente anche il centroide npc+kcut+1
*    ................SERVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!...
*   cambia poco tra il ridefinire solo il tappo o anche il primo dopo il
*   tappo
*
        do i = npc+1,npc+2*kcut
          yce(i) = (yv(i+1)+yv(i))/2.d0
          zce(i) = (zv(i+1)+zv(i))/2.d0
        end do
*
c - chiamo il solutore
*
c        if(jt.ge.150)then
c        kcut=0
c        kphi(npc+1)=1
c        endif
        call solver(jt,
     #        yv,zv,amp,phi,dphi,kphi,npt,npc,kcut,tg,wg,kg,i2d)
        call caveta(phi,amp,dpht,npc,npt,kcut,dphi)
c        if(jt.ge.150) kcut=1
        
*
c - definisco le variabili per la stampa
*
        do ip = 1,npt
          vym(ip,1) = dphi(ip)*rrny(ip) + dpht(ip)*ttmy(ip) 
          vzm(ip,1) = dphi(ip)*rrnz(ip) + dpht(ip)*ttmz(ip)
        end do
*
* ---  ENERGIA --------------  
*
*
        ekb  = 0.d0
        if(i2d.eq.1)then
          do jj=1,npt-3
            ekb = ekb + amp(jj)*phi(jj)*dphi(jj)
          enddo
        else
          do jj=1,npt-3
            ekb = ekb + amp(jj)*phi(jj)*dphi(jj)*yce(jj)
          enddo
          ekb = ekb*pi
        endif
        ekb = -ekb
        if(jt.gt.2)then
        if(ntagl.eq.1.and.ntaglo.ne.0)then
          dek = (ekb - ekbo)/dt
          write(66,'(i4,6d15.6)') jt,t,ekb,ekj,dek
        elseif(ntagl.eq.0)then
c  ---   estrapolazione da told a t 
c        con parabola ekbe = ekbo + b*dt + c*dt**2. cosi ho l'energia
c        che dovrebbe avere ekb  al tempo  t senza taglio
          dtco= dtv + dtvv
          a  = ekbo
          b1 = dtco**2*(ekboo-ekbo) - dto**2*(ekbooo-ekbo)
          b2 = -dtco**2*dtv + dtv**2*dtco
          b  = b1/b2
          c1 = dtco*(ekboo-ekbo) - dto*(ekbooo-ekbo)
          c2 = dtv**2*dtco - dtco**2*dtv
          c3  = c1/c2
          ekbe = a + b*dt + c3*dt**2
          ekjj = ekbe - ekb
          ekj  = ekj + ekjj  
          write(66,'(i4,6d15.6)') jt,t,ekb,ekj,dek
        endif
        endif
        ntaglo = ntagl
        dtvv  = dtv
        dtv   = dt
        ekbooo = ekboo
        ekboo  = ekbo
        ekbo   = ekb
*
* -------   FINE ENERGIA ---------
*
*
* -------- LUNGHEZZA STRISCIA DI CONTROLLO -----
*
        if(jend.eq.0)then
          rl=0.d0
          do i=jind+1,jfid-1
            dl= sqrt( (yce(i+1)-yce(i))**2 + (zce(i+1)-zce(i))**2 )
            rl= rl + dl
          enddo
          dlin= (1.d0-tin)*sqrt( (yce(jind+1)-yce(jind))**2 +
     #                    (zce(jind+1)-zce(jind))**2 )
          dlfi= tfi*sqrt( (yce(jfid+1)-yce(jfid))**2 +
     #                    (zce(jfid+1)-zce(jfid))**2 )
          rl= rl + dlin + dlfi
          yin = yce(jind)+tin*(yce(jind+1)-yce(jind))
          zin = zce(jind)+tin*(zce(jind+1)-zce(jind))
          yfi = yce(jfid)+tfi*(yce(jfid+1)-yce(jfid))
          zfi = zce(jfid)+tfi*(zce(jfid+1)-zce(jfid))
          write(77,'(i4,d15.6,d25.15,2i4,4d15.6)') jt,t,rl,jind,jfid+1,
     #               yin,zin,yfi,zfi                             

        endif
* ---------
        write(*,*)  '   ---> PRESSIONE <---'
*
* --------------------------------------------------------------
* la prebem calcola la pressione all'instante di tempo
* precedente (xche' uso schema centered per l'accelerazione),
* invece press all'istante corrente,
* sarebbe + logico calcolarle allo stesso istante
* --------------------------------------------------------------
*
*        if(jt.le.2)then
*          write(36,*) 10.d0
*          write(36,*)
*          write(36,*)
*        endif

        if(jt.gt.2)then
*
*  - CALCOLO PRESSIONE alla 'zhao-faltinsen'
*
        call press(iwig,dt,ux,vfall,phi,dphi,dpht,phio,
     #                     npc,npco,yce,zce,ttmy,ttmz,rrny,rrnz,
     #                     yco,zco,vyo,vzo,tto,amp,
     #                     pres )
*
c - CALCOLO PRESSIONE  con BEM (updated solo esplicito!!!)
*
        if((npcoo.eq.npco.and.npco.eq.npc).and.ntagl.eq.1)then
        ipre=ipre+1
        call prebem(jt,proat,acce,accl,ttmz,rrnz,dpht,cax,i2d,
     #                  tg,wg,kg,pres,yco,zco,ycoo,zcoo,rrny,
     #                  rrnyo,rrnyoo,rrnzo,rrnzoo,phio,dphio,
     #                  dphto,vyo,vzo,ttmzo,ttmyo,ampo,yvo,zvo,zvoo0,
     #                  npco,npto,dto,ipre)
        endif
*
c       - Aggiorno la velocita' di caduta
*
          write(*,*) 'VFALL ',vfall
          vfall = vfall + acce*dt
          vsopr = vsopr + accl*dt
          dzet  = dzet + (vfall-vsopr)*dt
*
          voat(jt) = vfall
          voso(jt) = vsopr
          foel(jt) = rkk*dzet
          write(*,*) 'VFALL ',vfall
*
        endif
*
c       - stampa configurazioni
*
        mm = mod(jt,ksta)
*
*        if((npcoo.eq.npco.and.npco.eq.npc).and.ntagl.eq.1)then
        if ( (tsta.gt.0.d0.and.isp.eq.1).or.
     &       (tsta.lt.0.d0.and.mm.eq.0)     ) then
          call stampa(proat,pres,ntagl,nrid,npco,yco,zco,ampo,
     #                jend,jind,jfid,yin,zin,yfi,zfi)         
        end if
*        endif
*
        isp = ist
*  -  Aggiorno yco ycoo rrnyo rrnyoo rrnzo rrnzoo
        npcoo= npco
        npco = npc
        npto = npt
        dto  = dt
        zvoo0= zvo0
        zvo0 = zv(1)
        do i = 1,npc
          yvo(i)  = yv(i)
          zvo(i)  = zv(i)
          ycoo(i) = yco(i)
          zcoo(i) = zco(i)
          yco(i)  = yce(i)
          zco(i)  = zce(i)
          rrnyoo(i) = rrnyo(i)
          rrnzoo(i) = rrnzo(i)
          rrnyo(i)  = rrny(i)
          rrnzo(i)  = rrnz(i)
          ttmyo(i)  = ttmy(i)
          ttmzo(i)  = ttmz(i) 
          ampo(i)   = amp(i)
          vyo(i)   = vym(i,1)
          vzo(i)   = vzm(i,1)
          phio(i)  = phi(i)
          dphio(i) = dphi(i)
          dphto(i) = dpht(i)
        enddo
        yvo(npc+1)  = yv(npc+1)
        zvo(npc+1)  = zv(npc+1)
        do i = npc+1,npt
          yvo(i)  = yv(i)
          zvo(i)  = zv(i)
          ampo(i) = amp(i)
          zco(i)  = zce(i)
          vyo(i)   = vym(i,1)
          vzo(i)   = vzm(i,1)
        end do         
*
*
*
      end do
*
      close(30)
      close(35)
      close(36)
      close(41)
      close(42)
      close(43)
      close(44)
      close(45)
      close(46)
      close(47)
      close(77)
      stop
      end
