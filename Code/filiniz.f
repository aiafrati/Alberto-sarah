
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine filiniz ()

c     - Inizializza i coefficienti del filtro di Dold

      include 'slam_p.h'
c     include 'slam_v.h'
      include 'slam_f.h'

      real*8 a(1:7,0:7)
      data a(1,0),a(1,1),a(1,2),a(1,3),a(1,4),a(1,5),a(1,6),a(1,7)
     &    /2.d0  ,-1.d0 ,0.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0/
      data a(2,0),a(2,1),a(2,2),a(2,3),a(2,4),a(2,5),a(2,6),a(2,7)
     &    /6.d0  ,-4.d0 ,1.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0/
      data a(3,0),a(3,1),a(3,2),a(3,3),a(3,4),a(3,5),a(3,6),a(3,7)
     &    /20.d0 ,-15.d0,6.d0  ,-1.d0 ,0.d0  ,0.d0  ,0.d0  ,0.d0/
      data a(4,0),a(4,1),a(4,2),a(4,3),a(4,4),a(4,5),a(4,6),a(4,7)
     &    /70.d0 ,-56.d0,28.d0 ,-8.d0 ,1.d0  ,0.d0  ,0.d0  ,0.d0/
      data a(5,0),a(5,1) ,a(5,2),a(5,3),a(5,4),a(5,5),a(5,6),a(5,7)
     &    /252.d0,-210.d0,120.d0,-45.d0,10.d0 ,-1.d0 ,0.d0  ,0.d0/
      data a(6,0),a(6,1) ,a(6,2),a(6,3) ,a(6,4),a(6,5),a(6,6),a(6,7)
     &    /924.d0,-792.d0,495.d0,-220.d0,66.d0 ,-12.d0,1.d0  ,0.d0/
      data a(7,0) ,a(7,1)  ,a(7,2) ,a(7,3)  ,a(7,4),a(7,5),a(7,6),a(7,7)
     &    /3432.d0,-3003.d0,2002.d0,-1001.d0,364.d0,-91.d0,14.d0 ,-1.d0/

      integer*4 e(1:7)
      data e(1),e(2),e(3),e(4),e(5),e(6),e(7)
     &    /2   ,4   ,6   ,8   ,10  ,12  ,14/

      integer*4 i,j

      do i=1,7
         d(i) = 2.d0**e(i)
         do j=0,7
            c(i,j) = a(i,j)
         end do
      end do

      return
      end