
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine doldfil1 (v,ni,nf,ntt)

c     - filtra i valori di v

      include 'slam_p.h'
      include 'slam_f.h'

      real*8    v(npamx+1),vaux(npamx+1)
      real*8    delta

      integer*4 i,j

      do i=1,ntt
         vaux(i) = v(i)
      enddo

      do 110 i=ni,nf
         v(i) = v(i) - delta (vaux(i-iford))
110   continue

      return
      end
