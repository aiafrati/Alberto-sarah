
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      function delta (f)

c     - Filtro di Dold

      include 'slam_f.h'

      real*8    f(-iford:iford),delta

      integer*4 i

      delta = f(0)*c(iford,0)
      do i = 1,iford
         delta = delta + (f(i)+f(-i))*c(iford,i)
      end do
      delta = delta/d(iford)

      return
      end
