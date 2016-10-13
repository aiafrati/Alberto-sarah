
* ----------------------------------------------
*
      double precision function der2(ff,fc,fb,dxf,dxb)
      implicit double precision (a-h,o-z)
*
* - derivata seconda alle diff finite. dxf = xf-xc,  dxb = xc-xb,
*   der2 = 2.*( (f(xf)-f(xb))/dxf - (f(xc)-f(xb))/dxb ) / (dxf+dxb)
*   in pratica der2 = der( der(ff,fc,dxf), der(fc,fb,dxb), 0.5*(dxf+dxb) )
*
      der2 = 2.*( der(ff,fc,dxf) - der(fc,fb,dxb) )/( dxf + dxb )
*
      return
      end

