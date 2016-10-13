
      

c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine input(i2d,kg)

      include"slam_p.h"
      include"slam_v.h"
      include"slam_f.h"
*
      write(*,*) '.......... :  input'
*
c     - dati letti dal file "slam.inp":
*
c       --  iwig, = 1 per 2d+t wigley, 2 per 2d+t, 0 per 2d (o axi)
c       -- i2d,kg    = 2d , pti gauss (se i2d=0) 
*       -- kini    = 1= inizia da sol esistente, 0= inizia da 0
c       -- alpha,ztr = (iwig>0) angolo trim, immersione transom in assetto 
c       -- vfall0, ux  = velocita' di caduta o di avanzamento (iwig=0, o > 0)
c       -- pro0,dx0 = immersione iniziale o dx per sezione iniziale (iwig=0, >0)
c       -- ampp    = ampiezza pannelli in prossimita' del corpo
c       -- pfraz   = rapporto ampp e ampiezza limite inferiore
c       -- ancut   = angolo minimo di taglio jet in gradi
*       -- ng    = n. pti ricostruiti sul corpo 
c       -- frdt    = "cfl" massimo, o dt
c       -- escr    = coefficiente dilatazione pannelli
c       -- estr    = estensione a monte del dominio a destra
c       -- tend    = tempo finale simulazione
c       -- tsta    = salto temporale tra due stampe (<0 ogni passo)
c       -- ksta    = salto passi tra due stampe
c       -- scon    = estensione files salvati configurazioni (2 lettere)
c       -- svel    = estensione files salvati velocita' (2 lettere)
c       -- sprt    = estensione files salvati pressione (2 lettere)
c       -- sadi    = estensione files salvati conf. adim. (2 lettere)
c       -- sad2    = estensione files salvati conf. adim. (2 lettere)
c       -- spot    = estensione files salvati potenziale (2 lettere)
c       -- sfor    = estensione file forze (5 lettere)
c       -- lrest   = opzione restart (Si=1)
c       -- rmu     = coefficiente di massa adimensionale
c       -- kitm    = numero it. pressione (0=Esplicito,-1=Implicito)
c       -- ift     = salto passi nel tempo per filtro di Dold
c       -- iford   = ordine filtro di Dold
c       -- rga     = rapporto (m sopra/m sotto)
c       -- vea     = velocita elastica adimensionale (rho*v0^2/rkk)
*
      open(8,file='slam.inp',status='UNKNOWN')
        read(8,*) iwig
        read(8,*) i2d,kg 
        read(8,*) kini 
        read(8,*) alpha,ztr 
        read(8,*) vfall0,ux
        read(8,*) pro0,dx0
        read(8,*) ampp 
        read(8,*) pfraz
*        read(8,*) alfa 
        read(8,*) ancut 
        read(8,*) ng 
*        read(8,*) rappi
        read(8,*) frdt 
        read(8,*) escr 
        read(8,*) estr 
*        read(8,*) tend 
        read(8,*) tsta
        read(8,*) ksta
        read(8,'(a)') scon
        read(8,'(a)') svel
        read(8,'(a)') sprt
        read(8,'(a)') sadi
        read(8,'(a)') sad2
        read(8,'(a)') spot
        read(8,'(a)') sfor
        read(8,*) lrest
        read(8,*) rmu
        read(8,*) kitm
        read(8,*) ift 
        read(8,*) iford
        read(8,*) rga
        read(8,*) vea
      close(8)
*
        pro0 = -pro0
*
        write(*,*) '   iwig   = ',iwig
        write(*,*) '   i2d    = ',i2d 
        write(*,*) '   alfa   = ',alfa 
        write(*,*) '   vfall0 = ',vfall0 
        write(*,*) '   ancut  = ',ancut
*
*
*
      return
      end

