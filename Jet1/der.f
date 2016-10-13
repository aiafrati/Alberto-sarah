
* ----------------------------------------------
*
      double precision  function der(ff,fb,dx)
      implicit double precision (a-h,o-z)
*
* -  derivata alle diff. finite : der = (f(xb)-f(xa))/(xb-xa) 
*
      der = (ff - fb)/dx
*
      return
      end
*
