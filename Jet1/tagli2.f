*
*------------------------------------------------------------------
      subroutine tagli2(kcut,ttmy,ttmz,rrny,rrnz,aamy,aamz,anint,
     #                  npt,npc,yn,zn,ygn,zgn,ng,amp,ycn,zcn,
     #                  phin,dphi,dpht,ek,amplim,ancut)
*
      include "slam_p.h"
c      include "slam_v.h"
      dimension ttmy(npamx),ttmz(npamx),rrny(npamx),rrnz(npamx)
      dimension aamy(npamx),aamz(npamx)
      dimension yn(npamx+1),zn(npamx+1),ygn(npamx+1),zgn(npamx+1)
      dimension amp(npamx),ycn(npamx),zcn(npamx)
      dimension phin(npamx),dphi(npamx),dpht(npamx)
*
          write(*,*) '.......... :  tagli2',kcut
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
      if (anint.lt.ancut.or.dino.lt.amplim.and.anint.lt.60.d0) then
*      if (anint.lt.ancut.or.dinom.lt.amplim.and.anint.lt.90.d0) then
*      
        write(*,*) '   taglio ',kpac-npc-kcut,' pannelliii '
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
* ---  calcolo posizione ynap rispetto al centroide kint
        if(znap.gt.zcn(kint))then
          kmid = kint
        else
          kmid = kint-1
        endif
* ---  dpht e dphi come interpol. tra i vicini        
        d1  = sqrt((znap - zcn(kmid))**2 + (ynap-ycn(kmid))**2)
        d2  = sqrt((znap - zcn(kmid+1))**2 + (ynap-ycn(kmid+1))**2)
        sb  = d2
        write(*,*) 'sb = ' ,sb
        d1  = d1/(d1+d2)
        d2  = 1.d0 - d1
        vyuu = dpht(kmid+1)*ttmy(kmid+1)+dphi(kmid+1)*rrny(kmid+1) 
        vyud = dpht(kmid)*ttmy(kmid)+dphi(kmid)*rrny(kmid) 
        vyu  = d2*vyud + d1*vyuu
        vzuu = dpht(kmid+1)*ttmz(kmid+1)+dphi(kmid+1)*rrnz(kmid+1) 
        vzud = dpht(kmid)*ttmz(kmid)+dphi(kmid)*rrnz(kmid) 
        vzu  = d2*vzud + d1*vzuu
        phiu = d2*phin(kmid) + d1*phin(kmid+1)
        kmi2 = kpac
        d1  = sqrt((zk - zcn(kmi2))**2 + (yk-ycn(kmi2))**2)
        d2  = sqrt((zk - zcn(kmi2+1))**2 + (yk-ycn(kmi2+1))**2)
        d1  = d1/(d1+d2)
        d2  = 1.d0 - d1
        write(*,*) 'd1,d2 ',d1,d2
        vydu = dpht(kmi2+1)*ttmy(kmi2+1)+dphi(kmi2+1)*rrny(kmi2+1) 
        vydd = dpht(kmi2)*ttmy(kmi2)+dphi(kmi2)*rrny(kmi2) 
        vyd  = d2*vydu + d1*vydd
        vzdu = dpht(kmi2+1)*ttmz(kmi2+1)+dphi(kmi2+1)*rrnz(kmi2+1) 
        vzdd = dpht(kmi2)*ttmz(kmi2)+dphi(kmi2)*rrnz(kmi2) 
        vzd  = d2*vzdu + d1*vzdd
        phid = d2*phin(kmi2) + d1*phin(kmi2+1)
        vy  = 0.5d0*(vyu + vyd)
        vz  = 0.5d0*(vzu + vzd)
        phic = 0.5d0*(phiu+phid)
* ---- calcolo energia getto
c        ek = 0.d0
        do kk=1,kpac-kint
          jj  = kint + kk
          vvn = dphi(jj)*phin(jj)
          ek  = ek + amp(jj)*vvn
          write(*,*) 'jj = ',jj, dphi(jj),phin(jj),amp(jj),ek
        enddo
         write(*,*) 'ek 1 ',ek
* ---  aggiungo pezzetto sul corpo
        vvn  = dphi(kmid)*phin(kmid)
        write(*,*) 'vvn,sb ',vvn,sb,kmid,dphi(kmid),phin(kmid) 
        ek   = ek + vvn*sb
         write(*,*) 'ek 2 ',ek
* --- aggiungo tappo nuovo
        st   = sqrt((ynap - yk)**2 + (znap - zk)**2)
        rny  = -(zk-znap)/st
        rnz  = (yk-ynap)/st
        ty   = (yk-ynap)/st
        tz   = (zk-znap)/st
        vn   = vy*rny + vz*rnz
        vn1   = vyuu*rny + vzuu*rnz
        vn2   = vyud*rny + vzud*rnz
        vn3   = vydu*rny + vzdu*rnz
        vn4   = vydd*rny + vzdd*rnz
        write(*,*) 'tap nuov ',vn,phic
        write(*,*) 'tap nuov ',vn1,vn2
        write(*,*) 'tap nuov ',vn3,vn4
        vvn  = phic*vn
        vvn  = -dphi(npc+1)*phin(npc+1)
        write(*,*) dphi(npc+1),phin(npc+1),st
        ek   = ek + vvn*st         
         write(*,*) 'ek 3 ',ek
*           
*          
c         - shift delle variabili in caso di ritaglio
*
        write(*,*) 'kpac,npc',kpac,npc,npt
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
        end if
*
        write(*,*) 'npt',npt
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
