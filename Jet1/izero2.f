
*
**********************************************************
*
      subroutine izero2(iint,yi,zi,yn,zn,nstart,nend,ck)
*****************************************************************************
*
* questa subroutine calcola il pto di intersezione tra una retta
* y=ck e un corpo
* definito da nc pannelli (enc+1 nodi).
*
* OUT: iint : n. pannello su cui giace il pto di intersezione
*      yi   : coord. y del pto di intersezione 
*      zi   : coord. z del pto di intersezione
*
*****************************************************************************
      include "slam_p.h"
*      implicit real*8 (a-h,o-z)
*      parameter(npamx=1000)
      parameter (eps=1.d-5)
      dimension yn(npamx),zn(npamx)
*
      ipar= -1000
      iint= -999
      do 10,i=nstart,nend
         y1 = yn(i)
*         write(*,*) 'y1 ',y1
         y2 = yn(i+1)
         z1 = zn(i)
         z2 = zn(i+1)
         cm2 = (z1-z2)/(y1-y2)
*         cc = abs((cm2-cm1)/cm1)
*         if(cc.le.eps)then
*           write(*,*) ' izero: rette parallele !!! ', i
*           ipar = i
*         endif
         if(abs((y2-y1)/max(abs(y2),abs(y1))).le.eps)then
           write(*,*) 'izero2, rette parallele!!! '
           ipar=i
         endif
         v  = ( ck - y1 )/( y2 - y1 ) 
         if(v.ge.0.d0.and.v.le.1.d0)then
           yi = ck 
           zi = (1.d0 - v)*z1 + v*z2
           iint = i
           goto 99
         endif
 10   continue
 99   continue
      if(iint.eq.-999)then
         write(*,*) 'izero2, intersezione non trovata!' 
      endif
*
      return
      end

