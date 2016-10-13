
*--------------------------------------------------------------
*
      double precision function distan(q1,q2,w1,w2,p1,p2)
*
      include"slam_p.h"
*
* distanza  di un punto (p1,p2) da una retta M=Q+xi*w, xi parametro,
* w vettore che definisce la direzione della retta
*
      w   = sqrt(w1**2 + w2**2)
      wn1 = w1/w
      wn2 = w2/w
      wp1 = wn2
      wp2 = -wn1
      qp1 = q1 - p1
      qp2 = q2 - p2
*
      distan=abs(wp1*qp1 + wp2*qp2)
*
      return
      end

