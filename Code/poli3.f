*
      subroutine poli3(a,b,c,x1,x2,x3)

      include "slam_p.h"
* ---------------------------------------------------------------     
*
* -- calcolo degli zeri (reali!) di un` eq. di 3o grado
*
*     x**3 + a*(x**2) + b*x + c = 0
*     
*     y = x + a/3 
*
*     Referenza : G.Gorni, L'equazione di terzo grado,
*                 dispensa di 7 pp. trovata in rete.
*                 G.Gorni e' ordinario a Matematica a Udine.
*-----------------------------------------------------------------
      pi = acos(-1.d0)
* --- soluzione 
      p = -(a**2)/3.d0 + b
      q =  2.d0*(a**3)/27.d0 - a*b/3.d0 + c
      write(*,*) ' p q ',p,q
*
      delta = (q**2)/4.d0 + (p**3)/27.d0
*      
*    
      if(delta.ge.0.d0)then
        u1  = (-q/2.d0 + sqrt(delta))
        if(u1.lt.0d0)then
          u1 = -(-u1)**(1.d0/3.d0)
        else
          u1 = u1**(1.d0/3.d0)
        endif
        v1  = (-q/2.d0 - sqrt(delta)) 
        if(v1.lt.0d0)then
          v1 = -(-v1)**(1.d0/3.d0)
        else
          v1 = v1**(1.d0/3.d0)
        endif
        x1  = -a/3.d0 + u1 + v1
        x2  = 1.d31
        x3  = 1.d31
      else
        if(q.le.0.d0)then
          theta = atan(-2.d0*sqrt(-delta)/q)
        else
          theta = pi + atan(-2.d0*sqrt(-delta)/q)
        endif
        x1 = -a/3.d0 + 2.d0*sqrt(-p/3.d0)*cos(theta/3.d0)       
        x2 = -a/3.d0 + 2.d0*sqrt(-p/3.d0)*cos((theta+2.d0*pi)/3.d0)
        x3 = -a/3.d0 + 2.d0*sqrt(-p/3.d0)*cos((theta+4.d0*pi)/3.d0)
      endif
*
*      write(*,*) 'delta   : ',delta
*      write(*,*) 'u1  v1  : ',u1,v1
*      write(*,*) 'x1      : ',x1
*      write(*,*) 'x2      : ',x2
*      write(*,*) 'x3      : ',x3
*
      return
      end


      
