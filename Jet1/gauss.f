
      subroutine gauss(a,nmax,n,m,eps)
*************************************************
*soluzione sistema lineare con met. gauss pivot
*
************************************************ 
      include "slam_p.h"
      dimension a(nmax,*)
      do 10 j=1,n-1
      k=j
      do 20 i=j+1,n
      if(abs(a(k,j))-abs(a(i,j)))30,20,20
   30 k=i
   20 continue
      if(k-j)50,40,50
   50 do 60 l=j,m
      c=a(j,l)
      a(j,l)=a(k,l)
   60 a(k,l)=c
   40 if(abs(a(j,j))-eps)120,120,70
   70 do 80 k=j+1,m
      p=a(j,k)/a(j,j)
      do 80 i=j+1,n
   80 a(i,k)=a(i,k)-a(i,j)*p
   10 continue
      if(abs(a(n,n))-eps)120,120,90
   90 do 100 il=n+1,m
      do 100 j=n,1,-1
      a(j,il)=a(j,il)/a(j,j)
      do 100 i=1,j-1
  100 a(i,il)=a(i,il)-a(i,j)*a(j,il)
      return
  120 write(*,500)eps
  500 format(5x,'pivot inferieur a ',1p,e16.6)
      stop
      end
***************************+++
**************************
