
*----------------------------------------------------------------------------
*
 
      subroutine izero(iint,yi,zi,yn,zn,nstart,nend,cm,cq)
*****************************************************************************
*
* questa subroutine calcola il pto di intersezione tra una retta e un corpo
* definito da nc pannelli (e nc+1 nodi).
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
         y2 = yn(i+1)
         z1 = zn(i)
         z2 = zn(i+1)
         cm2 = (z1-z2)/(y1-y2)
*         write(*,*) 'cm2 ' ,cm2
         if(cm.gt.eps)then
           cc = abs((cm2-cm)/cm)
           if(cc.le.eps)then
             write(*,*) ' izero: rette parallele !!! ', i
             ipar = i
           endif
         endif
         v  = ( z1 - cm*y1 - cq )/( z1 - z2 - cm*(y1-y2) ) 
         if(v.ge.0.d0.and.v.le.1.d0)then
           yi = (1. - v)*y1 + v*y2
           zi = (1. - v)*z1 + v*z2
           iint = i
           goto 99
         endif
 10   continue
 99   continue 
      if(ipar.eq.iint) then
        write(*,*) ' intersezione con rette quasi parallele ', iint
      endif
      if(iint.eq.-999)then
        write(*,*) ' izero: intersezione non trovata!  '
      endif
*
      return
      end
