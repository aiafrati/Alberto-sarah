
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine input(vfall0,ancut,iiget,jjget,frint,rmg,
     #       epsgg,eskg,gfrac,pro0 ,frdt,ampp,pfraz,escr,estr,
     #       tend,ksta,scon,svel,spot,spre,idis,ift,iford,frtend,ramii,
     #       ramiii,eskkk)

      include"slam_p.h"

      character*2 scon,svel,spre,spot

      write(*,*) '--> input.........'
c      include"slam_f.h"
c     - dati letti dal file "slam.inp":
*
c       -- vfall0  = velocita' di caduta
c       -- ancut   = angolo minimo di taglio jet in gradi
c       -- iiget   = passo attivazione griglia x modello getto
c       -- jjget   = passo attivazione modello getto
c       -- frint   = matching region (0=no matching, frint<0.9))
c       -- rmg    = num. pannelli base getto
c       -- ramii  = multiplo di amii, griglia bordi punto separazione 
c       -- ramiii  = multiplo di amii, griglia vertice corpo 
c       -- eskkk  = stretching pannelli dal vertice                   
c       -- espgg  =limite inferiore pannelli SL vicino corpo (fraz. di ampli
c       -- eskg    = stretching pannelli getto 
c       -- gfrac   = frazione getto con angolo<ancut (gfrac<=1.)
c       -- pro0    = immersione iniziale                  
c       -- frdt    = "cfl" massimo, o dt
c       -- ampp    = ampiezza pannelli in prossimita' del corpo
c       -- pfraz   = rapporto ampp e ampiezza limite inferiore
c       -- escr    = coefficiente dilatazione pannelli bulk
c       -- estr    = estensione a monte del dominio a destra
c       -- tend    = tempo finale simulazione
c       -- frtend  = dtmax (in frazione tempo finale simulazione)
c       -- ksta    = salto passi tra due stampe
c       -- scon    = estensione files salvati configurazioni (2 lettere)
c       -- svel    = estensione files salvati velocita' (2 lettere)
c       -- spot    = estensione files salvati potenziale (2 lettere)
c       -- spre    = estensione files salvati pressione (2 lettere)
c       -- idis    = intervallo azione disun2 
c       -- ift     = salto passi nel tempo per filtro di Dold
c       -- iford   = ordine filtro di Dold
*
      open(8,file='slam.inp',status='UNKNOWN')
        read(8,*) vfall0
        read(8,*) ancut
        read(8,*) iiget
        read(8,*) jjget
        read(8,*) frint
        read(8,*) rmg   
        read(8,*) ramii 
        read(8,*) ramiii 
        read(8,*) eskkk 
        read(8,*) epsgg 
        read(8,*) eskg 
        read(8,*) gfrac
        read(8,*) pro0 
        read(8,*) frdt 
        read(8,*) ampp 
        read(8,*) pfraz
        read(8,*) escr 
        read(8,*) estr 
        read(8,*) tend 
        read(8,*) frtend 
        read(8,*) ksta
        read(8,'(a)') scon
        read(8,'(a)') svel
        read(8,'(a)') spot
        read(8,'(a)') spre
        read(8,*) idis
        read(8,*) ift 
        read(8,*) iford
      close(8)

      return
      end
