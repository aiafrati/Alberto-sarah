      subroutine gristr(x0,x1,a1,an,c,n) 
      implicit double precision(a-h,o-z)
* IN = x0,x1, a1,an
* OUT= a1 (modificato!), c,n
*
*      x0 = 0.d0
*      x1 = 1.d0
*      at = x1-x0
*      a1 = at/float(100)
*      an = 5.516d-2 
*      c  = 1.05d0
*
* -- calcolo a partire da a1 e c --> ricavo n 
*      n = ngri1(at,a1,c)
*      x = x0
*      dx = a1
*      do 10, i=1,n
*        x=x+dx
*        write(*,*) i,x,dx
*        dx=dx*c 
* 10   continue
*
* --- calcolo a partire da a1 e an ---.  ricavo n e c, poi risistemo a1
*    per non avere avanzi di pannello alla fine
*
      write(*,*) 
      write(*,*) 
      n = ngri2(at,a1,an)
      c = cgri(at,a1,an)
      a1 = at/(c**n - 1.d0)*(c-1.d0)
*      x = x0
*      dx = a1
*      do 20, i=1,n-1
*        x=x+dx
*        write(*,*) i,x,dx
*        dx=dx*c 
* 20   continue
*      x=x1
*      write(*,*) i,x,dx
      return 
      end
*----------------------------------------------------------
      function ngri1(at,a1,c)
      implicit double precision(a-h,o-z)
*     at      amp.totale segmento 
*     a1      amp primo pannello
*     c       coeff. stretching
*     ngrist  numero pannelli 
*
      ngri1=log(1.d0 + at*(c-1.d0)/a1)/log(c)      
*
      return
      end
*----------------------------------------------------------
      function ngri2(at,a1,an)
      implicit double precision(a-h,o-z)
*     at amp.totale segmento 
*     a1  amp primo pannello
*     an  amp primo pannello
*     c   coeff. stretching
*     ngri2   numero pannelli 
*
      c = (at-a1)/(at-an)
      ngri2 = 1+log(an/a1)/log(c)
*
      return
      end
*----------------------------------------------------------
      function cgri(at,a1,an)
      implicit double precision(a-h,o-z)
      cgri = (at-a1)/(at-an)
      return
      end
*----------------------------------------------------------
 
