c     --  file da includere nel programma slam.f  --

c     - contiene la definizione delle variabili e i parametri per il
c       dimensionamento dei vettori

      implicit real*8 (a-h,o-z)
      character*2 scon,svel,sprt,sadi,spot,sad2
      character*5 sfor

      parameter (npamx=1100, ntmx = 400000 , ngmax=65 )

      common/costanti/pi
