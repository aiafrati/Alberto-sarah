*
      subroutine wigl(alpha,al,bl,tl,x0,yg,zg,rnxg,ng,dzm) 
*      implicit double precision (a-h,o-z)
      include "slam_p.h"
      dimension yg(npamx),zg(npamx),rnxg(npamx)
      pi = acos(-1.d0)
*
*      alpha = 8.d0/180.d0*pi
      z0    = 0.d0*tl
      ca    = cos(alpha)
      sa    = sin(alpha)
      tga   = tan(alpha)
*
*      tl   = .05d0
*      al   = .8d0
*      bl   = .08d0
*
      zi   = -tl
      dzi  = 3.d0*tl/float(ng)
* 
*      write(*,*) 'wigl.f '
        do 20, ig = 1,ng+1
*
          xis     = x0/ca - (zi+z0)*tga
          yg(ig)  = fwig(xis,zi,al,bl,tl)
          zg(ig)  = -x0*tga + (zi+z0)*sa*tga + (zi+z0)*ca +dzm
*          zs     = -xis*sa + (zi+z0)*ca
          fwx       = fwigx(xis,zi,al,bl,tl)
          fwz       = fwigz(xis,zi,al,bl,tl)
          rn        = sqrt(1.d0 + fwx**2 + fwz**2)
          rnxg(ig)  = -(fwx*ca + fwz*sa)/rn
*
*          write(*,*) yg(ig), zg(ig)
          zi= zi + dzi
 20     continue
*
      return
      end
* -----------------------------------------------------------*
      function fwig(x,z,al,bl,tl)
      implicit double precision (a-h,o-z)
       z0=z
       if (z.gt.0.d0) z0=0.d0
       fwig = 0.5d0*bl*( 1.d0 - (2.d0*x/al)**2 )*( 1.d0 - (z0/tl)**2 )
      return
      end      
      function fwigx(x,z,al,bl,tl)
      implicit double precision (a-h,o-z)
       z0=z
       if (z.gt.0.d0) z0=0.d0
       fwigx = -bl*x/(0.5d0*al)**2 *( 1.d0 - (z0/tl)**2)
      return
      end      
      function fwigz(x,z,al,bl,tl)
      implicit double precision (a-h,o-z)
      z0=z
      if (z.gt.0.d0) then
        fwigz  = 0.d0
      else
        fwigz = -bl*z0/(tl**2) *( 1.d0 - (2.d0*x/al)**2 )
      endif
      return
      end      
*
