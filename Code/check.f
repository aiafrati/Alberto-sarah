
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine check(yv,zv,npc,npsl,dt,t,jt,ampsl,frdt,
     #    vymsl1,vzmsl1,ycsl,zcsl,kget,frtend,tend,zgn,ngo,vfall,ngo1,
     #     nsep)

      include"slam_p.h"
      dimension yv(npamx),zv(npamx),ampsl(npamx)
      dimension ycsl(npamx),zcsl(npamx),vymsl1(npamx),vzmsl1(npamx)
      dimension zgn(npamx)
      write(*,*) '--> check......... '
c     - calcolo volume

c      cmm = -cmm0
c      do i = 1,npc
c        cmm = cmm + (yv(i+1)-yv(i))*(zv(i+1)+zv(i))/2.d0
c      enddo
c      do i = 1,npsl
c        cmm = cmm + (ysl(i+1)-ysl(i))*(zsl(i+1)+zsl(i))/2.d0
c      enddo
c      write(*,*) ' variazione volume = ',cmm

c     - calcolo passo integrazione temporale

      dtk = dt
      dta = 1.d+10
c      if(jt.lt.iiget)then
c        dyb = yv(npc)-yv(npc+1)
c        dzb = zv(npc)-zv(npc+1)
c        amb = sqrt(dyb**2+dzb**2)
c        dyb = dyb/amb
c        dzb = dzb/amb
c        dyf = ysl(2)-ysl(1)
c        dzf = zsl(2)-zsl(1)
c        amf = sqrt(dyf**2+dzf**2)
c        dyf = dyf/amf
c        dzf = dzf/amf
c        arc = acos(dyb*dyf+dzb*dzf)
c        if(arc.gt.0.5d0*pi)then
c          frdtt= 0.02d0 + (frdt-0.02d0)*(pi-arc)/(0.5d0*pi)
c          write(*,*) 'frdt frdtt ' ,frdt,frdtt,arc/pi*180.d0
c        endif  
c      endif
      if(jt.le.10)then
        frdtt=frdt*jt/10
      else
        frdtt = frdt
      endif
      do i = 1,npsl
        vvy = vymsl1(i)
        vvz = vzmsl1(i)
        vvm = sqrt(vvy*vvy + vvz*vvz)
        dtc = frdtt*ampsl(i)/vvm
c        write(24,'(i5,6d15.6)') i,frdtt,ampsl(i),vvy,vvz,dtc,dta
        dta = min(dtc,dta)
      enddo
      dtc = frtend*tend
      write(*,*) 'dtc dta ',dtc,dta
      dta = min(dta,dtc)
c     - verifica

      cm1   = (zv(npc+1)-zv(npc))/(yv(npc+1)-yv(npc))
      cq1   = zv(npc)-cm1*yv(npc) 
      dista = cq1+cm1*ycsl(1)-zcsl(1)
      velco = vzmsl1(1)-cm1*vymsl1(1)
      write(*,*) 'cm1 cq1..',cm1,cq1,dista,velco

      if (velco.gt.0.d0.and.kget.eq.0) then
        dtsu = 0.2*dista/velco
        if (dtsu.lt.dta) then
          dta = 0.8*dtsu
        end if
      end if
      dt = dta
      write(*,*) 'dta2 ',dta
      if(nsep.eq.1) dt=dt/2.d0
*
      if (dt.lt.1.d-15) stop ' time step nullo o negativo '
      write(*,*) ' passo,dt,angolo,velocita '
      write(*,*) jt,sngl(dt),sngl(vfall)
* - aggiorno posizione corpo 
      do i=1,ngo
        zgn(i)=zgn(i)-vfall*dt
      enddo
*

      return
      end
