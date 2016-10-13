* ----------------------------------------------------------
*
      subroutine stampa(proat,pres,ntagl,nrid,npco,yco,zco,ampo,
     #                  jend,jind,jfid,yin,zin,yfi,zfi)

      include"slam_p.h"
      include"slam_v.h"
      character*6 vro
      dimension pres(npamx),yco(npamx),zco(npamx),ampo(npamx)
*
      write(*,*) '.......... :  stampa'
*
    1 format('#  tempo = ',e16.8)
*
      call ropn (llf,scon,vro)
      open(21,file=vro,status='UNKNOWN')
        write(21,1) t
        do iv = 1,npt+1
          write(21,*) sngl(yv(iv)),sngl(zv(iv)),sngl(rnxx(iv))
        enddo
        write(21,*)
        write(21,*)
        if(jend.eq.0)then
          write(21,*) sngl(yin),sngl(zin)
          do iv = jind+1,jfid
            write(21,*) sngl(yce(iv)),sngl(zce(iv))
          enddo
          write(21,*) sngl(yfi),sngl(zfi)
        endif

      close(21)
*
      sca = 20.*dt/frdt
      call ropn (llf,svel,vro)
      open(22,file=vro,status='UNKNOWN')
        write(22,1) t
        do ip = 1,npt
          write(22,*) sngl(yce(ip)),sngl(zce(ip)),
     &                sngl(yce(ip)+sca*vym(ip,1)),
     &               sngl(zce(ip)+sca*vzm(ip,1))
        enddo
      close(22)
*
c      if(ntagl*nrid.gt.0)then
      call ropn (llf,sprt,vro)
      open(22,file=vro,status='UNKNOWN')
        write(22,1) t-dt
        tts=0.d0
        dtts=0.5d0*ampo(1)
        do ip = 1,npco
          tts = tts + dtts
c          write(22,*) sngl(zad),sngl(tts),sngl(pre1(ip)),
c     &                sngl(-pre2(ip)),sngl(pre1(ip)+pre2(ip))
          write(22,'(8e15.6)') sngl(tts),sngl(yco(ip)),sngl(zco(ip)),
     #        sngl(pre1(ip)),sngl(-pre2(ip)),sngl(pre12(ip)+pre2(ip)),
     &        sngl(pres(ip)),sngl(pre1(ip)+pre2(ip))
           dtts=0.5d0*(ampo(ip)+ampo(ip+1))
        enddo
      close(22)
c      endif
*
      call ropn (llf,sadi,vro)
      open(22,file=vro,status='UNKNOWN')
      adim  = abs(vfall0*zv(1))  
      adim2 = abs(vfall0*(t+t0))
        do iv = 1,npt
          yad    = yce(iv)/abs(zv(1))
          zad    = zce(iv)/abs(zv(1))
          phiad  = phi(iv)/adim
          dphiad = dphi(iv)/adim
          write(22,'(4e16.6)') sngl(yad),sngl(zad),
     #                  sngl(phiad) , sngl(dphiad)        
        enddo
      close(22)
*
      call ropn (llf,sad2,vro)
      open(22,file=vro,status='UNKNOWN')
        do iv = 1,npt+1
          yad    = yv(iv)/abs(zv(1))
          zad    = zv(iv)/abs(zv(1))
          write(22,'(2e16.6)') sngl(yad),sngl(zad)
        enddo
      close(22)
*
      call ropn (llf,spot,vro)
      open(22,file=vro,status='UNKNOWN')
        do ip = 1,npt
          ycc = (yv(ip)+yv(ip+1))/2.d0
          zcc = (zv(ip)+zv(ip+1))/2.d0
          write(22,*) sngl(ycc),sngl(zcc),sngl(phi(ip)),sngl(dphi(ip))
          if (ip.eq.npc) write(22,*)
        enddo
      close(22)
*
      llf = llf + 1
*
      return
      end

