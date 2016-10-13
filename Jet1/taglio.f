*
*------------------------------------------------------------------
      subroutine taglio(k2cut,ttmy,ttmz,rrny,rrnz,aamy,aamz,anint,ntagl,
     #                 jind,jfid,jend)
*
      include "slam_p.h"
      include "slam_v.h"
      dimension ttmy(npamx),ttmz(npamx),rrny(npamx),rrnz(npamx)
      dimension aamy(npamx),aamz(npamx)
*
          write(*,*) '.......... :  taglio',kcut
*
          t1y   = ttmy(npc) 
          t1z   = ttmz(npc) 
          rny   = rrny(npc)
          rnz   = rrnz(npc)
*
      if (kcut.eq.0) then
*
          t2y   = ttmy(npc+1) 
          t2z   = ttmz(npc+1) 
          sca   = (t1y*t2y + t1z*t2z) 
          anint = (pi-acos(sca))*180.d0/pi
          segno = rny*t2y + rnz*t2z
          segno = segno/abs(segno)
          anint = anint*segno
          dinom = rrny(npc)*aamy(npc+1)+rrnz(npc)*aamz(npc+1)
          kpac  = npc+1
*
      else
*
c         - individuo il pannello che forma l'angolo minore col corpo
*           (ed ev. la minima distanza normale )
*
          iin   = npc+2
          anint = 90.d0
          anloc = 0.d0
          dinom = 100.d0*amplim
          tttt = sqrt((yn(npc+3)-yn(npc+2))**2+(zn(npc+3)-zn(npc+2))**2)
          write(*,*) 'kcut anloc anint', kcut, sngl(anloc),sngl(anint)
          do while (anloc.lt.90.d0)
            t2y   = ttmy(iin) 
            t2z   = ttmz(iin)
            sca   = (t1y*t2y + t1z*t2z) 
            anloc1 = (pi-acos(sca))*180.d0/pi
            segno = rny*t2y + rnz*t2z
            segno = segno/abs(segno)
            anloc = anloc1*segno
            write(*,*) 'anloc ', t, iin, sngl(anloc),sngl(anint)
            yk    = yn(iin)
            zk    = zn(iin)
            dinloc= 1.d30
            do j=1,ng-1
              dy = yk-ygn(j)
              dz = zk-zgn(j)
              dis= sqrt(dy**2+dz**2)
              dinloc= min(dis,dinloc)
            enddo  
            if (anloc.lt.anint)then
*            if (anloc.lt.ancut)then
*            if (anloc.lt.anint.or.dinloc.lt.amplim) then
              anint = anloc
              dinom = dinloc
              kpac  = iin
            end if
            iin = iin + 1
          end do
          write(*,*)
          write(*,*)
      end if
*
      dino  = aamy(npc+1)*rrny(npc) +  
     &          aamz(npc+1)*rrnz(npc)
      dino  = dino + float(kcut)
*
      anint=90.d0
      if (anint.lt.ancut.or.dino.lt.amplim.and.anint.lt.60.d0) then
*      if (anint.lt.ancut.or.dinom.lt.amplim.and.anint.lt.90.d0) then
*      
        ntagl=0
        k2cut=k2cut+1
        write(*,*) '   taglio ',kpac-npc-kcut,' pannellij'
        yor   = yn(npc+1)
        zor   = zn(npc+1)
        aor   = 0.d0
        if (kcut.eq.1) then
          aor   = amp(npc+1)
        endif
        yk  = yn(kpac+1)
        zk  = zn(kpac+1) 
        call milan(iint,ygn,zgn,1,ng-1,yk,zk)
* - iint usato in ridis...(per interp. spline su yg zg...)
* -  calcolo volume tagliato. meglio usare sempre ynapp etc...!
        call milan(kint,yn,zn,1,npc,yk,zk)
        cm1 = (zn(kint+1)-zn(kint))/(yn(kint+1)-yn(kint))
        cq1 = zn(kint)-cm1*yn(kint)
        cm2 = -(yn(kint+1)-yn(kint))/(zn(kint+1)-zn(kint))
        cq2 = zk - cm2*yk
        ynap  = (cq2-cq1)/(cm1-cm2)
        znap  = cm1*ynap + cq1
        voljj = (yn(kint+1)-ynap)*(zn(kint+1)+znap)/2.d0
        do ii = kint+2, kpac+1
           dvoljj = (yn(ii)-yn(ii-1))*(zn(ii)+zn(ii-1))/2.d0
           voljj  = voljj + dvoljj 
        end do
        voljj = voljj + (ynap-yn(kpac+1))*(znap+zn(kpac+1))/2.d0
        volj  = volj + voljj
*          
c         - shift delle variabili in caso di ritaglio
*
        if (kcut.eq.1) then
          ish = kpac-npc-1
          do ip = kpac+1,npt
            ic        = ip - ish
            phin(ic)  = phin(ip)
            ycn(ic)   = ycn(ip)
            zcn(ic)   = zcn(ip)
            yn(ic)    = yn(ip)
            zn(ic)    = zn(ip)
          end do
          yn(npt-ish+1) = yn(npt+1)
          zn(npt-ish+1) = zn(npt+1)
          npt = npt - ish
* --  shift indici striscia di controllo
          if(jind.lt.kpac+3.and.jend.eq.0)then
            write(*,*)'attenzione, la striscia se ne va!!!'
            jend=1
          elseif(jend.eq.0)then 
            jind=jind-ish
            jfid=jfid-ish
          endif
* --  
        end if
*
        yn(npc+1)   = ynap
        zn(npc+1)   = znap
*
c         - tolgo dal volume originale la parte tagliata
*
        ampk = sqrt((yn(npc+2)-yn(npc+1))**2+
     &                (zn(npc+2)-zn(npc+1))**2)
*
        dkta  = sqrt( (yn(npc+1)-yor)**2 + (zn(npc+1)-zor)**2 )
        aknu  = ampk
        voljjj= dkta*(aknu+aor)/2.d0
        cmm0  = cmm0 - dkta*(aknu+aor)/2.d0
        kcut  = 1
*
*       write(*,*) 'voljj,voljjj     ', voljj,voljjj 
       write(*,*) '   anint = ', anint     
  
*
      end if
*
*
      return
      end
*
