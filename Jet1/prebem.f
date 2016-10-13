      subroutine prebem(jt,proat,acce,accl,ttmz,rrnz,dpht,cax,i2d,
     #                  tg,wg,kg,pres,yco,zco,ycoo,zcoo,rrny,
     #                  rrnyo,rrnyoo,rrnzo,rrnzoo,phio,dphio,
     #                  dphto,vyo,vzo,ttmzo,ttmyo,ampo,yvo,zvo,zvoo0,
     #                  npco,npto,dto,ipre)
*
      include "slam_p.h"
      include "slam_v.h"
*      parameter
      dimension taa(npamx),yv2(npamx),zv2(npamx)
      dimension ttmz(npamx),rrnz(npamx),accnn(npamx)
      dimension dpht(npamx),fosl2(ntmx)
      dimension tg(ngmax),wg(ngmax)
      dimension pres(npamx)
      dimension yco(npamx),zco(npamx),ycoo(npamx),zcoo(npamx)
      dimension rrny(npamx),rrnyo(npamx),rrnyoo(npamx)
      dimension rrnzo(npamx),rrnzoo(npamx)
      dimension ttmzo(npamx),ttmyo(npamx),ampo(npamx)
      dimension phio(npamx),dphio(npamx),dphto(npamx)
      dimension vyo(npamx),vzo(npamx),yvo(npamx),zvo(npamx)
      write(*,*) '.......... :  prebem ' 
*   
*       
*        if (kitm.ge.0) then  ! Processo Iterativo
*
          kit  = 0
          if (jt.lt.6) then
            acce = acce
          else  ! assegno accelerazione di tentativo come media
            acce = +ggr - (0.25*( fosl(jt-1)+fosl(jt-2)+fosl(jt-3)+
     &              fosl(jt-4) ) +dzet*rkk) /weir
          end if
  100     continue
          aten = acce
*
c         -- Calcolo d/dt di dphi/dn sul corpo:
*
*       --- interpolazione spline
            taa(1) = 0.d0
            do i = 2,npco+1
              dy = yvo(i) - yvo(i-1)
              dz = zvo(i) - zvo(i-1)
              taa(i) = taa(i-1) + sqrt(dy**2 + dz**2)
            end do
*
            yp1=1.d+31
            ypn=1.d+31
            call splone(yv2,taa,yvo,yp1,ypn,npco+1,npamx)  
            call splone(zv2,taa,zvo,yp1,ypn,npco+1,npamx) 
            cax = 1.d0 - float(i2d) 
*
            write(*,*) 'npco ',npco
*            write(36,*) '# ipre, jt,t ',ipre,jt,t
* - vel wz0 del primo vertice:
            wz0u = (zv(1) - zvo(1))/dt
            wz0d = (zvo(1) - zvoo0)/dto
            wz0  = 0.5d0*(wz0u + wz0d)
            do i = 1,npco
*
              dy  = 0.5d0*(yvo(i+1) - yvo(i))
              dz  = 0.5d0*(zvo(i+1) - zvo(i))
              tpo = taa(i) + sqrt(dy**2 + dz**2)
              call splont(tpo,ypo,taa,yvo,yv2,npco+1,npamx,
     #                  yt,1,ytt,1,yttt,1)
              call splont(tpo,zpo,taa,zvo,zv2,npco+1,npamx,
     #                  zt,1,ztt,1,zttt,1)
              call curv(ypo,yt,ytt,yttt,zpo,zt,ztt,zttt,cur,rn2sy,rn2sz,
     #                  conc)
              thet = atan(abs(dz/dy))
              curp = sin(thet)/(yco(i))
              if(iwig.gt.-1)then
                cwig = 1.d0
* - w: vel corpo; variabile da punto a punto 
                wuy = (yce(i) - yco(i))/dt
                wuz = (zce(i) - zco(i))/dt
                wdy = (yco(i) - ycoo(i))/dto
                wdz = (zco(i) - zcoo(i))/dto
                wy  = 0.5d0*(wuy + wdy)
                wz  = 0.5d0*(wuz + wdz)
                wt = wy*ttmyo(i) + wz*ttmzo(i)
                wn = wy*rrnyo(i) + wz*rrnzo(i)
* - accelerazione corpo
                acy = (wuy - wdy)/(0.5d0*(dt+dto))
                acz = (wuz - wdz)/(0.5d0*(dt+dto))
                accnn(i) = acy*rrnyo(i) + acz*rrnzo(i)
* - derivata della normale, Dn / Dt
                dnuy = (rrny(i) -rrnyo(i))/dt 
                dnuz = (rrnz(i) -rrnzo(i))/dt 
                dndy = (rrnyo(i) -rrnyoo(i))/dto 
                dndz = (rrnzo(i) -rrnzoo(i))/dto
                dny  = 0.5d0*(dnuy + dndy) 
                dnz  = 0.5d0*(dnuz + dndz)
              else
                cwig = 0.d0
* - w: vel corpo; okkio vfall varia con la dinamica 
                wy = 0.d0
                wz = -vfall
                wt = wz*ttmzo(i)
                wn = wz*rrnzo(i)
* - accelerazione corpo
                accnn(i) = aten*rrnzo(i) ! calcolo esplicito
              endif
* - uso potenziale vecchio
              phii = phio(i)
              dphii= dphio(i)
              dphti= dphto(i)
              ampi = ampo(i)
              ttmzi= ttmzo(i)
              yci  = yco(i)
              if(i.lt.npco)then
* - w: vel corpo; variabile da punto a punto 
                wupy = (yce(i+1) - yco(i+1))/dt
                wupz = (zce(i+1) - zco(i+1))/dt
                wdpy = (yco(i+1) - ycoo(i+1))/dto
                wdpz = (zco(i+1) - zcoo(i+1))/dto
                wpy  = 0.5d0*(wupy + wdpy)
                wpz  = 0.5d0*(wupz + wdpz)
                wpt = wpy*ttmyo(i+1) + wpz*ttmzo(i+1)
                wpn = wpy*rrnyo(i+1) + wpz*rrnzo(i+1)
                phiu = phio(i+1)
                dphiu= dphio(i+1)
                dphtu= dphto(i+1)
                ampu = ampo(i+1)
              endif
              if(i.gt.1)then
* - w: vel corpo; variabile da punto a punto 
                wumy = (yce(i-1) - yco(i-1))/dt
                wumz = (zce(i-1) - zco(i-1))/dt
                wdmy = (yco(i-1) - ycoo(i-1))/dto
                wdmz = (zco(i-1) - zcoo(i-1))/dto
                wmy  = 0.5d0*(wumy + wdmy)
                wmz  = 0.5d0*(wumz + wdmz)
                wmt = wmy*ttmyo(i-1) + wmz*ttmzo(i-1)
                wmn = wmy*rrnyo(i-1) + wmz*rrnzo(i-1)
                phid  = phio(i-1)
                dphid = dphio(i-1)
                dphtd = dphto(i-1)
                ampd  = ampo(i-1)
              endif
              if(i.eq.npco)then
                ampdd = ampo(i-2)
                phidd = phio(i-2)
              endif
* - termine wt*dun/dtau  , s ascissa curvilinea
* - se w e' globale (iwig=0), wt*dun/dtau = wt*w.dn/dtau 
              if(iwig.gt.-1)then
                if(i.eq.1)then
                  dwns = der(wpn,wn,0.5d0*(ampi+ampu))
                elseif (i.eq.npco)then
                  dwns = der(wn,wmn,0.5d0*(ampi+ampd))
                else
                  dwnsu = der(wpn,wn,0.5d0*(ampi+ampu))
                  dwnsd = der(wn,wmn,0.5d0*(ampi+ampd))
                  dwns = 0.5d0*(dwnsu + dwnsd)
                endif
              endif
* - termine wn*dut/dtau = dphi*dut/dtau; dut/dtau:= dfft (o dfft2)
              if (i.eq.1) then
c                 write(41,*) 'i=1, wz,wz0 ',
c     #               sngl(t),sngl(wz),sngl(wz0u),sngl(wz0d),sngl(wz0)
                 dff2 = der(phiu,phii,0.5d0*(ampi+ampu) )
                 dff1 = wz0*ttmzi
                 dfft = (dff2 - dff1)/ampi
                 dfft2= (dff2 - dff1)/ampi
                 dphit=0.5d0*(wz*ttmzi+der(phiu,phii,0.5d0*(ampu+ampi)))
              else if (i.eq.npco) then
                 dfff = 0.d0
                 dfft = der(dphti,dphtd,0.5d0*(ampi+ampd) )
                 dfft2 = der2( phii,phid,phidd,
     #           0.5d0*(ampi+ampd),0.5d0*(ampd+ampdd) )
                 dphit = dphti
              else
                 dfff = der(dphtu,dphti,0.5d0*(ampi+ampu) )
                 dffb = der(dphti,dphtd,0.5d0*(ampi+ampd) )
                 dfft = 0.5d0*(dfff+dffb)
                 dfft2 = der2( phiu,phii,phid,
     #           0.5d0*(ampu+ampi),0.5d0*(ampi+ampd) )
                 dphit = dphti
              endif
*
*              accn = aten*cos(ande) ! calcolo esplicito
* ------  TANIZAWA ----- OBSOLETO, solo per iwig=0-
              ter1 = (-cur)*(dphit-wt)**2
              ter2 =  (cur)*dphii**2
              ter3 = cur*dphit**2
              ter4 = wn*dfft
              ter4b= wn*dfft2
              ter5 = -dphit*(wy*rn2sy+wz*rn2sz)
              ter6 =-(dphii*sin(thet)+dphit*cos(thet))/yci*dphii
* ------- home made --------------------------
              if(iwig.gt.-1)then
                terr6b = -wt*dwns
              else
                terr6b = -wt*(wy*rn2sy+wz*rn2sz)
              endif   
              terr1= cur*wn**2
              terr2= curp*wn**2*cax
              terr3= cur*dphit*wt
              terr4= wn*dfft2
              terr5=-dphit*cos(thet)/yce(i)*wn*cax
c     #               *sin(thet)
c              terr5= 0.d0
              terr6 = -dphit*(wy*rn2sy+wz*rn2sz)
              terr7 = dny*(wy-vyo(i)) + dnz*(wz-vzo(i))
* 
              if(conc.le.0.d0)then
c                dpnt(i) = accnn(i)+ter1+ter2+ter3+ter4b+ter5
c     #                    +ter6 + cax*curp*dphi(i)**2+ter6
                dpnt(i) = accnn(i)+terr1+terr2+terr3+terr4-terr5
     #                    +terr6b + cwig*terr7
              else
                dpnt(i) = accnn(i)-terr1+terr2-terr3+terr4-terr5
     #                    +terr6b + cwig*terr7
              endif
              kphi(i) = 0
*
*              write(36,'(9d15.7)')yco(i),zco(i),accnn(i),
*     # dpnt(i)-accnn(i)-cwig*terr7,terr7,dny,dnz,wy-vyo(i),wz-vzo(i)
            end do
*            write(36,*)
*            write(36,*)
*
c           -- Assegno d phi/dt sulla superficie libera
*
            do ip = npco+1,npto
              dpt(ip) = -(vyo(ip)**2+vzo(ip)**2)/2.d0 - ggr*zco(ip)
*              dpt(ip) = -(vym(ip,1)**2+vzm(ip,1)**2)/2.d0 - ggr*zce(ip)
              kphi(ip) = 1
            end do 
*
c           -- chiamo il solutore
*
           
           call solver(jt,
     #        yvo,zvo,ampo,dpt,dpnt,kphi,npto,npco,0,tg,wg,kg,i2d)
c           call solver(jt,
c     #        yv,zv,amp,dpt2,dpnt2,kphi,npt,npco,0,tg,wg,kg,i2d)
*
   
c           -- calcolo forze 
*
            tt(jt)    = t
            fosl(jt)  = 0.d0
            fosl2(jt) = 0.d0
            fosll     = 0.d0
            fosll2    = 0.d0
            if(iwig.eq.0)then
              vf2       = 0.5d0*vfall**2
            else
              vf2       = 0.5d0*ux**2
            endif 
              area1     = 0.d0
              area2     = 0.d0
              zw1       = 0.d0
              zw2       = 0.d0
*	    
              if(i2d.eq.1)then
c	      do i=1,npco
c	        area=area+amp(i)
c	      enddo
c	      area=area*2.d0
              do i = 1,npco
*              uu2      = (vym(i,1)**2 + vzm(i,1)**2)/2.d0
              uu2      = (vyo(i)**2 + vzo(i)**2)/2.d0
              pre1(i)  = - dpt(i)/vf2
c              pre12(i) = - dpt2(i)/vf2
              pre2(i)  = - uu2/vf2
              ptot     = pre1(i)+pre2(i) - ggr*zco(i)/vf2
*	      if(ptot.gt.0.d0)then
              fosl(jt) = fosl(jt) + ptot*ampo(i)*(-rrnzo(i))
              fosll    = fosll    + ptot*ampo(i)*(-rrnzo(i))
              area1    = area1 + ampo(i)
              zw1      = zco(i)
*	      endif
*	      if(pres(i).gt.0.d0)then
              fosl2(jt)= fosl2(jt) + pres(i)*ampo(i)*(-rrnzo(i))
              fosll2   = fosll2    + pres(i)*ampo(i)*(-rrnzo(i))
              area2    = area2 + ampo(i)
              zw2      = zco(i)
*	      endif
              end do
c     - effetto parte simmetrica
              fosl(jt)  = 2.d0*fosl(jt)
              fosl2(jt) = 2.d0*fosl2(jt) 
              fosll     = fosll/area1
              fosll2    = fosll2/area2
              zw1       = (zvo(1)-zw1)/zvo(1)
              zw2       = (zvo(1)-zw2)/zvo(1) 
              else
c	      do i=1,npco
c	        area=area+amp(i)*yce(i)
c	      enddo 
c	      area=area*2.d0*pi
              do i = 1,npco
*              uu2      = (vym(i,1)**2 + vzm(i,1)**2)/2.d0
              uu2      = (vyo(i)**2 + vzo(i)**2)/2.d0
              pre1(i)  = - dpt(i)/vf2
c              pre12(i) = - dpt2(i)/vf2
              pre2(i)  = - uu2/vf2
              ptot     = pre1(i)+pre2(i) - ggr*zco(ip)/vf2
*	      if(ptot.gt.0.d0)then
              fosl(jt) = fosl(jt) + ptot*ampo(i)*(-rrnzo(i))*yco(i)
              fosll    = fosll  + ptot*ampo(i)*(-rrnzo(i))
              area1    = area1+ampo(i)*yco(i)
              zw1      = zco(i)
*              endif
*	      if(pres(i).gt.0.d0)then
              fosl2(jt)= fosl2(jt) + pres(i)*ampo(i)*(-rrnzo(i))*yco(i)
              fosll2   = fosll2 + pres(i)*ampo(i)*(-rrnzo(i))
              area2    = area2+ampo(i)*yco(i)
              zw2      = zco(i)
*              endif
              end do
              fosl(jt)  = 2.d0*pi*fosl(jt)
              fosl2(jt) = 2.d0*pi*fosl2(jt)
              fosll     = fosll/area1
              fosll2    = fosll2/area2
              zw1       = (zvo(1)-zw1)/zvo(1)
              zw2       = (zvo(1)-zw2)/zvo(1) 
              endif
           write(35,*) sngl(tt(jt)),sngl(fosl(jt)),sngl(fosl2(jt)),
     #                   sngl(fosll),sngl(fosll2)
*	    
c   - - - Dinamica corpo rigido
*
            call dinam(proat,fosl(jt),acce,accl)
            
*
c         -- fine processo iterativo
*
            kit = kit + 1
            dva = abs(aten-acce)
*
          if (dva.gt.abs(0.01d0*acce).and.kit.lt.kitm) go to 100
*
*        else       ! Calcolo Implicito

*c         - Costruzione Vettore Termini Noti Lato Corpo
*c         - Parte Elastica (costante)
*
*          felt = rkk*dzet*cos(ande)/weir
*
*c         - Parte Idrodinamica da u**2 (costante)
*
*          fdit = 0.d0
*          do i = 1,npco
*            uu2     = (vym(i,1)**2 + vzm(i,1)**2)/2.d0
*            fdit    = fdit + uu2*amp(i)
*          end do
*          fdit = 2.d0*fdit*cos(ande)**2/weir
*        
*c         - Parte Idrodinamica Un d(Utau)/dtau
*
*          do i = 1,npco
*            if (i.eq.1) then
*              dff1 = - vfall*sin(ande)
*              dff2 = 2.d0*(phi(i+1)-phi(i))/(amp(i+1)+amp(i))
*              dfft = (dff2-dff1)/amp(i)
*            else if (i.eq.npco) then
*              dfft = 2.d0*(dpht(i)-dpht(i-1))/(amp(i-1)+amp(i))
*            else
*              dfff = 2.d0*(dpht(i+1)-dpht(i))/(amp(i+1)+amp(i))
*              dffb = 2.d0*(dpht(i)-dpht(i-1))/(amp(i-1)+amp(i))
*              dfft = (dfff+dffb)/2.d0 
*            end if
*
*            dpnt(i)  = fdit-felt+dphi(i)*dfft ! attento:solo per wedge
*          enddo
*
*c         - Costruzione Vettore Termini Noti Lato Superficie Libera
*
*          do i = npco+1,npt
*            dpt(i)  = -(vym(i,1)**2 + vzm(i,1)**2)/2.d0
*          enddo 
*
*c         - Chiamo il solutore
*
*          coef = 2.d0*cos(ande)**2/weir
*          call solpre(yv,zv,amp,npco,npt,dpnt,dpt,coef)
*
*c         - calcolo forze 
*
*          tt(jt)   = t
*          fosl(jt) = 0.d0
*          do i = 1,npco
*            uu2     = (vym(i,1)**2 + vzm(i,1)**2)/2.d0
*            pre1(i) = - dpt(i)
*            pre2(i) = - uu2
*            ptot    = pre1(i)+pre2(i)
*            fosl(jt)= fosl(jt) + ptot*amp(i)*cos(ande)
*          enddo
*
*c         - effetto parte simmetrica
*
*          fosl(jt) = 2.d0*fosl(jt)
*
*c         - Dinamica corpo rigido
*
*          call dinam(proat,fosl(jt),acce,accl)
*        end if
**
        return
        end
