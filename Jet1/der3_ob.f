
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++

      function der3(d11,d12,d13,d21,d22,d23,d31,d32,d33,t1,t2,t3)

      include"slam_p.h"

      dete = d11*d22*d33 + d12*d23*d31 + d21*d32*d13 -
     &       d13*d22*d31 - d23*d32*d11 - d12*d21*d33

      rnum = t1*d22*d33 + d12*d23*t3 + t2*d32*d13 -
     &       d13*d22*t3 - d23*d32*t1 - d12*t2*d33

      der3 = rnum/dete

      return
      end

