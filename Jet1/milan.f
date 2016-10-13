
*-----------------------------------------------------------------
*
      subroutine milan(icross,xn,yn,nstart,nend,xk,yk)
*
* questa subroutine calcola l'intersezione di due rette:
* Q = P(i+1) + tau*sip e Q=Pk + xi*si
*
* si e sip sono due vettori (unitari). sip e' perpendicolare a si, si
* e' il vettore definito come si=P(i+1)-P(i).
* si = (dlix,dliy)/dli qui sotto.
*
* xi e tau sono i parametri che definiscono il punto della retta.
*
* dopo l'intersezione c'e' un test per sapere se il punto Pk
* resta sempre dalla stessa parte della retta Pi + tau*sip,
* o se e' passato dall'altra parte rispetto al passo precedente
*
      include"slam_p.h"
      dimension xn(npamx),yn(npamx)
*
      xiold = -1000.
*      write(*,*) 'mailan ,,' ,nstart,nend
      do 10,i=nstart,nend
         dlix  = xn(i+1)- xn(i)
         dliy  = yn(i+1)- yn(i)
         dljkx = xn(i+1)  - xk
         dljky = yn(i+1)  - yk
         dli   = sqrt(dlix**2 + dliy**2)
         xi    = (dlix*dljkx + dliy*dljky)/dli
         tau   = (dlix*dljky - dliy*dljkx)/dli
*         write(*,*) 'xi i',xi,i
         if(xiold*xi.lt.0.)then
*         write(*,*) 'icrosss!!!! '
         icross= i
         goto 99
         endif
         xiold = xi
 10   continue
 99   continue
*
      return
      end
    
