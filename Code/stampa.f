
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      subroutine stampa(vfall,nget,kget,npc,npsl,jt,t,dt,frdt,llf,
     #                  scon,svel,spot,spre,phi,dphi,phisl,dphisl, 
     #                 yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,
     #                 vym1,vzm1,vymsl1,vzmsl1,
     #                 dpht,dphtsl,dpt,dpnt,pre,pre2,pres,vxi,
     #                 jend,jind,jfid,yin,zin,yfi,zfi)
*                      
      include "slam_p.h"   

      character*2 scon,svel,spre,spot

      dimension yv(npamx),zv(npamx),yce(npamx),zce(npamx)
      dimension vym1(npamx),vzm1(npamx),vymsl1(npamx),vzmsl1(npamx)
      dimension ycsl(npamx),zcsl(npamx),ysl(npamx),zsl(npamx)
      dimension phisl(npamx),dphisl(npamx)
      dimension dpht(npamx),dphtsl(npamx),phi(npamx),dphi(npamx)
      dimension dpt(npamx),dpnt(npamx)
      dimension pres(npamx),pre(npamx),pre2(npamx),vxi(npamx) 
      character*6 vro
      write(*,*) '--> stampa...........'
*
      sca   = 20.d0*dt/frdt
      nnget = nget*kget
      prof  = -zv(1)
    1 format('#  tempo = ',e16.8,i7)
* -- vertici 
      call ropn (llf,scon,vro)
      open(21,file=vro,status='UNKNOWN')
      write(21,1) t,jt
      do iv = 1,npc+nnget+1
        write(21,'(4e15.7)') sngl(yv(iv)),sngl(zv(iv)),
     #               sngl(yv(iv)/prof),sngl(zv(iv)/prof)
        if(iv.eq.npc+1) write(21,*)
      enddo
      if(nnget.gt.0)write(21,*)
      do iv = 1,npsl+1
        write(21,'(4e15.7)') sngl(ysl(iv)),sngl(zsl(iv)),
     #                sngl(ysl(iv)/prof),sngl(zsl(iv)/prof)
        if(iv.eq.nget) write(21,*)
      enddo
      write(21,*)
      write(21,*)
      if(jend.eq.0)then
        write(21,*) sngl(yin),sngl(zin)
        do iv = jind+1,jfid
          write(21,*) sngl(ycsl(iv)),sngl(zcsl(iv))
        enddo
          write(21,*) sngl(yfi),sngl(zfi)
      endif
      close(21)
* -- velocita
      call ropn (llf,svel,vro)
      open(22,file=vro,status='UNKNOWN')
      write(22,1) t,jt
      do ip = 1,npc+nnget
        write(22,*) sngl(yce(ip)),sngl(zce(ip))
        write(22,*) sngl(yce(ip)+sca*vym1(ip)),
     #              sngl(zce(ip)+sca*vzm1(ip))
        write(22,*)
      enddo
      do ip = 1,npsl
        write(22,*) sngl(ycsl(ip)),sngl(zcsl(ip))
        write(22,*) sngl(ycsl(ip)+sca*vymsl1(ip)),
     #              sngl(zcsl(ip)+sca*vzmsl1(ip))
        write(22,*)
      enddo
      close(22)
* -- centroidi e potenziale
      call ropn (llf,spot,vro)
      open(22,file=vro,status='UNKNOWN')
      write(22,1) t,jt
      do ip = 1,npc+nnget
        if(ip.eq.npc+1) write(22,*) 
        v2 = 0.5d0*(dphi(ip)**2+dpht(ip)**2)
        write(22,'(8e15.7)') sngl(yce(ip)),sngl(zce(ip)),
     #         sngl(yce(ip)/prof),sngl(zce(ip)/prof),
     #         sngl(phi(ip)),sngl(dphi(ip)),sngl(dpht(ip)),sngl(vxi(ip))
      enddo
      write(22,*)
      do ip = 1,npsl
        v2 = 0.5d0*(dphisl(ip)**2+dphtsl(ip)**2)
        write(22,'(8e15.7)') sngl(ycsl(ip)),sngl(zcsl(ip)),
     #         sngl(ycsl(ip)/prof),sngl(zcsl(ip)/prof),
     #         sngl(phisl(ip)),sngl(dphisl(ip)),
     #         sngl(dphtsl(ip)),sngl(v2)
        if(ip.eq.nget) write(22,*)
      enddo
      write(22,*)
* -- centroidi e pressione
      call ropn (llf,spre,vro)
      open(22,file=vro,status='UNKNOWN')
      write(22,1) t,jt
      do ip = 1,npc+nnget
        if(ip.eq.npc+1) write(22,*) 
        write(22,'(9e15.7)') sngl(yce(ip)),sngl(zce(ip)),
     #         sngl(yce(ip)/prof),sngl(zce(ip)/prof),
     #        sngl(dpt(ip)),sngl(dpnt(ip)),
     #        sngl(pre(ip)),sngl(pre2(ip)),sngl(pres(ip))
      enddo
      close(22)
* 
* 
      llf = llf+1
*
      return
      end
