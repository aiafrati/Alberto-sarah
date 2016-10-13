       kint=npc.'
        kint= npc
        yi  = yv(npc+1)
        zi  = zv(npc+1)
      endif
      volb   = 0.d0 
      do i = 1,kint-1
        volb = volb - (yv(i+1)-yv(i))*(zv(i+1)+zv(i))/2.d0
      end do
      volb   = volb - (yi -yv(kint))*(zi+zv(kint))/2.d0
      volb1  = volb - volb0
      if(jt.gt.1)then
      dvolab = dvol
      dvolabt= dvolab-dvolabo
      dvolre = dvol/volb1
      dvolret= dvolabt/(volb1-volb1o)
*      
*      write(*,*)
*      write(*,*) '   variazione volume abs e rel: ', dvolab, dvolre
*      write(*,*) '   variazione vol  da told a t: ', dvolabt, dvolret
*      write(*,*) '   volume spostato            : ', volb1
*      write(*,*) '   volume spostato da told a t: ', volb1-volb1o
*      write(*,*) '   perc. volume tagliata      : ', volj/volb1
*      write(*,*)
*
      endif
*
c     - calcolo passo integrazione temporale
*
      dtk = dt
      dta = 1.d+10
*
      do i = npc+kcut+1,npt
        vvy = vym(i,1)
        vvz = vzm(i,1)
        vvm = sqrt(vvy*vvy + vvz*vvz)
        dtc = frdt*amp(i)/vvm
        dta = min(dtc,dta)
      end do
*
c     - verifica 
*     - solo prima del taglio!
*
      if (kcut.eq.0) then    
        cm1   = (zv(npc+1)-zv(npc))/(yv(npc+1)-yv(npc))
        cq1   = zv(npc)-cm1*yv(npc) 
        dista = cq1+cm1*yce(npc+1)-zce(npc+1)
        velco = vzm(npc+1,1)-cm1*vym(npc+1,1)
        if (velco.gt.0.d0) then
          dtsu = 0.2d0*dista/velco
          if (dtsu.lt.dta) then
            dta = 0.8d0*dtsu
          endif
        endif
      endif
*
      dts = tust + tsta - t
      if (tsta.lt.0.d0) then
        dt = dta
        isp = 1
        ist = 1
      else if ((dts/dta).le.1.0) then
        dt = dts
        ist = 1
        tust = t + dts
      else if ((dts/dta).lt.2.0) then
        dt = dts/2.d0
        ist = 0
      else
        dt = dta
        ist = 0
      end if
*
* ------------------
*
        write(*,*) '   jt    = ',jt
        write(*,*) '   t old = ',t
        if (dt.lt.1.d-15) stop ' time step nullo o negativo '
        t = t + dt
        write(*,*) '   t new = ',t
        write(*,*) '   dt    = ',dt
*
*       - traslo il corpo a  proat  
*
        if(iwig.eq.0)then
*
        proat = zv(1) - vfall*dt 
        write(*,*) '   vfall = ',vfall
        write(*,*) '   proat = ',proat
        do ig=1,ng
           zgn(ig) = zgn(ig) - vfall*dt
        end do
        do iv=1,npc+1
           yn(iv)  = yv(iv) 
           zn(iv)  = zv(iv)  - vfall*dt
        end do
*
        else
*
        xxs0 = ux*t
        xsec = xxs0 + xs0
        wlfr = (xsec - x00)/wl
        write(*,*) '    xsec = ' , xsec 
        write(*,*) '    wlfr = ' , wlfr
        if(iwig.eq.1)then
          call wigl(alpha,al,bl,tl,xsec,ygn,zgn,rnxg,ng,dztr)
        elseif(iwig.eq.2)then
          call wigl(alpha,al,bl,tl,xsec,ygn,zgn,rnxg,ng,dztr)
          call geopla(xsec,dztr,alpha,ygn,zgn,rnxgf,ng,al,bl,tl)
        endif
* 
*        write(36,*) '# ',jt 
*        do i=1,ng
*          write(36,*) ygn(i),zgn(i),i
*        enddo
*        write(36,*)
*        write(36,*)
*
        tag(1) = 0.d0
        do ip=2,ng
           dy = ygn(ip)-ygn(ip-1)
           dz = zgn(ip)-zgn(ip-1)
           tag(ip)=tag(ip-1)+sqrt(dy**2+dz**2)
        end do
*
        yp1 = 1.d31
        ypn = 1.d31
        call splone(ygs2,tag,ygn,yp1,ypn,ng,npamx)
        call splone(zgs2,tag,zgn,yp1,ypn,ng,npamx)
        call splone(rnxgs2,tag,rnxg,yp1,ypn,ng,npamx)
*
        endif
*
      return
      end
