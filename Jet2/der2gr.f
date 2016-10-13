      program main
      implicit double precision (a-h,o-z)         
c      parameter
c      dimension 
*
      aa = 1.d0
      bb = 2.d0
      cc = 3.d0
*
      t0  = 0.d0
      yy0 = aa + bb*t0  + cc*t0**2
      t1  = 4.d0 
      yy1 = aa + bb*t1  + cc*t1**2
      t2  = 6.d0 
      yy2 = aa + bb*t2  + cc*t2**2
*
      call der2gr(yy0,yy1,yy2,t0,t1,t2,a,b,c)
*
      write(*,*) 'abc0 ' ,aa,bb,cc
      write(*,*) 'abc1 ' ,a,b,c


      stop
      end
        
* ----------------------------------------------------
      subroutine der2gr(y0,y1,y2,t0,t1,t2,a,b,c)
      implicit double precision (a-h,o-z)        
*
      eps  = 1.d-12
*
      dt10 = t1 - t0
      dt20 = t2 - t0
      dt21 = t2 - t1
      if(dt10.le.eps) stop 't1 - t0 troppo piccolo'
      if(dt20.le.eps) stop 't2 - t0 troppo piccolo'
      if(dt21.le.eps) stop 't2 - t1 troppo piccolo'
*
      a = y0*t1*t2/(dt20*dt10) - y1*t0*t2/(dt21*dt10)
     #   +y2*t0*t1/(dt20*dt21)
*     
      b = -y0*(t2+t1)/(dt20*dt10) + y1*(t2+t0)/(dt21*dt10)
     #    -y2*(t1+t0)/(dt20*dt21)
*
      c = y0/(dt20*dt10) - y1/(dt21*dt10) + y2/(dt20*dt21)
*
      return
      end 
