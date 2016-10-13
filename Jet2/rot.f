*---------------------------------------------------------*
      subroutine rot(beta,x,y,rx,ry,n)
*
*  rx,ry sono le coord. del punto in un sistema di assi ruotati di
*  angolo beta in senso antiorario 
*
      implicit double precision(a-h,o-z)
      dimension x(*),y(*),rx(*),ry(*)
*
      cb=cos(beta)
      sb=sin(beta)
      do 10,i=1,n
        rx(i) = x(i)*cb + y(i)*sb
        ry(i) =-x(i)*sb + y(i)*cb
  10  continue
*
      return
      end
*---------------------------------------------------------*

