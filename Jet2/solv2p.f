
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine solv2p(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,
     #           ycsl,zcsl,dphi,phisl,dphtsl,dphtbsl,phb,
     #           xigs,zegs,xigb,zegb,xigf,zegf,
     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
     #           ph,dphn,dphnb,a,b,c,d,e,phi,phib,dphisl,dphibsl,rl,    
     #           mb,mf,mt,m,n,nt,ntt,xi,ze,xis,zes,jt)
*
      include"slam_p.h"
      dimension dphn(npamx),dphi(npamx),ph(npamx),dphtb(npamx)
      dimension dphtsl(npamx)
      dimension phisl(npamx),phb(npamx),dphnb(npamx),dphtbsl(npamx)  
      dimension dphisl(npamx),dphibsl(npamx),phi(npamx),phib(npamx)
      dimension yc(npamx),zc(npamx),yce(npamx),zce(npamx)
      dimension yb(npamx),zb(npamx),yf(npamx),zf(npamx)
      dimension yv(npamx),zv(npamx),ysl(npamx),zsl(npamx)
      dimension ycsl(npamx),zcsl(npamx),am(npamx)
      dimension xig(-npamx:npamx),xigs(-npamx:npamx)
      dimension zeb(-npamx:npamx),zegs(-npamx:npamx)
      dimension zef(-npamx:npamx)
      dimension xi(npamx),xis(npamx),ze(npamx),zes(npamx)
      dimension aa(npamx,npamx),bb(npamx),appo(npamx)
      dimension ff(npamx),rl(npamx)
      dimension a(npamx),b(npamx),c(npamx),d(npamx),e(npamx)
      dimension indx(npamx)
      dimension xigb(-npamx:npamx),zegb(-npamx:npamx)
      dimension xigf(-npamx:npamx),zegf(-npamx:npamx)
      dimension ty(npamx) ,tz(npamx),ry(npamx),rz(npamx) 
      dimension tmy(npamx),tmz(npamx),rny(npamx),rnz(npamx) 
      dimension tmysl(npamx),tmzsl(npamx),rnysl(npamx),rnzsl(npamx) 
      parameter (naux=3*npamx)
      character*1 trans
c
c      dphn dphi yc zc yce zce yb zb yf zf yv zv ysl zsl ycsl zcsl
c      ph dphtb phisl
c      dh xig xigs zeb zebs zef zefs hp xi xis ze zes 
c      am aa ff bb rl
c      phb dphnb a b c d e 
      write(*,*) '--> solv2.............'
*
* - riorganizzo array punti yf yb zy zb yc zc am e phi dphn dpht
*                           dphnb = deriv normale bulk
*                           dphtb = deriv tangenz bulk
*                           phb   = phi nel bulk
      m  = int(frint*ng)
      n  = ng - m
      mb = npc+ng-(m+n)
      mf = npsl-(m+n)
      mt = mb + mf
      nt = mt + 2*m
      ntt= nt + 2*n
      write(*,*) npc,ng,npsl
      write(*,*) mb,m,n
      write(*,*) mf,m,n
      write(*,*) mt,nt,ntt
*
      do i=1,mb+m+n
        if(i.le.mb)then
          dphn(i) = dphi(i)
          dphnb(i)= dphi(i)
c          yc(i)   = yce(i)
c          zc(i)   = zce(i)
          yb(i)   = yv(i) 
          zb(i)   = zv(i) 
          yf(i)   = yv(i+1) 
          zf(i)   = zv(i+1)
          yc(i)   = 0.5d0*(yf(i)+yb(i))
          zc(i)   = 0.5d0*(zf(i)+zb(i))
          ry(i)   = rny(i)
          rz(i)   = rnz(i)
          ty(i)   = tmy(i)
          tz(i)   = tmz(i)
c          write(71,'(i6,6d15.7)') i,yc(i),zc(i),
c     #                yb(i),zb(i),yf(i),zf(i) 
c     #                ry(i),rz(i),ty(i),tz(i) 
        elseif(i.ge.mb+1.and.i.le.mb+m)then
          dphn(mt+i-mb)= dphi(i)
          dphnb(mt+i-mb)= dphi(i)
c          yc(mt+i-mb ) = yce(i)
c          zc(mt+i-mb)  = zce(i)
          yb(mt+i-mb)  = yv(i) 
          zb(mt+i-mb)  = zv(i) 
          yf(mt+i-mb)  = yv(i+1) 
          zf(mt+i-mb)  = zv(i+1) 
          yc(mt+i-mb)  = 0.5d0*(yf(mt+i-mb)+yb(mt+i-mb))
          zc(mt+i-mb)  = 0.5d0*(zf(mt+i-mb)+zb(mt+i-mb))
          ry(mt+i-mb)= rny(i)
          rz(mt+i-mb)= rnz(i)
          ty(mt+i-mb)= tmy(i)
          tz(mt+i-mb)= tmz(i)
c          write(72,'(i6,6d15.7)') mt+i-mb,yc(mt+i-mb),zc(mt+i-mb),
c    #                yb(mt+i-mb),zb(mt+i-mb),yf(mt+i-mb),zf(mt+i-mb)   
c     #                dphn(mt+i-mb),dphnb(mt+i-mb)  
c     #                ry(mt+i-mb),rz(mt+i-mb),ty(mt+i-mb),tz(mt+i-mb) 
                       
        else
          dphn(nt+i-(mb+m))= dphi(i) 
          dphnb(nt+i-(mb+m))= dphi(i) 
c          yc(nt+i-(mb+m))  = yce(i)
c          zc(nt+i-(mb+m))  = zce(i)
          yb(nt+i-(mb+m))  = yv(i) 
          zb(nt+i-(mb+m))  = zv(i) 
          yf(nt+i-(mb+m))  = yv(i+1) 
          zf(nt+i-(mb+m))  = zv(i+1) 
          yc(nt+i-(mb+m))  = 0.5d0*(yf(nt+i-(mb+m))+yb(nt+i-(mb+m)))
          zc(nt+i-(mb+m))  = 0.5d0*(zf(nt+i-(mb+m))+zb(nt+i-(mb+m)))
          ry(nt+i-(mb+m))  = rny(i) 
          rz(nt+i-(mb+m))  = rnz(i) 
          ty(nt+i-(mb+m))  = tmy(i) 
          tz(nt+i-(mb+m))  = tmz(i) 
c          write(73,'(i6,6d15.7)')
c     #    nt+i-(mb+m),yc(nt+i-(mb+m)),zc(nt+i-(mb+m)), 
c     #  yb(mt+i-(mb+m)),zb(mt+i-(mb+m)),yf(mt+i-(mb+m)),zf(mt+i-(mb+m)),
c     #           dphn(nt+i-(mb+m)),dphnb(nt+i-(mb+m))   
c     #  ry(nt+i-(mb+m)),rz(nt+i-(mb+m)),ty(nt+i-(mb+m)),tz(nt+i-(mb+m)) 
        endif  
      enddo    
c      write(71,*)
c      write(71,*)
c      write(72,*)
c      write(72,*)
c      write(73,*)
c      write(73,*)
      do i=1,mf+m+n
        if(i.le.n)then
          dphtbsl(i)=dphtsl(i)
          dphtb(nt+n+i)=dphtbsl(i)
c          dpht(nt+n+i)= dphtsl(i)
          ph(nt+n+i)  = phisl(i)
          phb(nt+n+i)  = phisl(i)
c          yc(nt+n+i)  = ycsl(i)
c          zc(nt+n+i)  = zcsl(i)
          yb(nt+n+i)  = ysl(i) 
          zb(nt+n+i)  = zsl(i) 
          yf(nt+n+i)  = ysl(i+1) 
          zf(nt+n+i)  = zsl(i+1)
          yc(nt+n+i)  = 0.5d0*(yf(nt+n+i)+yb(nt+n+i))
          zc(nt+n+i)  = 0.5d0*(zf(nt+n+i)+zb(nt+n+i))
          ry(nt+n+i)  = rnysl(i)
          rz(nt+n+i)  = rnzsl(i)
          ty(nt+n+i)  = tmysl(i)
          tz(nt+n+i)  = tmzsl(i)
c          write(83,'(i6,6d15.7)') nt+n+i,yc(nt+n+i),zc(nt+n+i),
c     #                yb(nt+n+i),zb(nt+n+i),yf(nt+n+i),zf(nt+n+i)    
c     #                dphtb(nt+n+i),ph(nt+n+i),phb(nt+n+i) 
c     #             ry(nt+n+i),rz(nt+n+i),ty(nt+n+i),tz(nt+n+i)    
        elseif(i.ge.n+1.and.i.le.n+m)then
          dphtbsl(i)=dphtsl(i)
          dphtb(mt+m+(i-n))= dphtbsl(i)
c          dpht(mt+m+(i-n)) = dphtsl(i)
          ph(mt+m+(i-n))   = phisl(i)
          phb(mt+m+(i-n))   = phisl(i)
c          yc(mt+m+(i-n))   = ycsl(i)
c          zc(mt+m+(i-n))   = zcsl(i)
          yb(mt+m+(i-n))   = ysl(i) 
          zb(mt+m+(i-n))   = zsl(i) 
          yf(mt+m+(i-n))   = ysl(i+1) 
          zf(mt+m+(i-n))   = zsl(i+1) 
          yc(mt+m+(i-n))   = 0.5d0*(yf(mt+m+(i-n))+yb(mt+m+(i-n)))
          zc(mt+m+(i-n))   = 0.5d0*(zf(mt+m+(i-n))+zb(mt+m+(i-n)))
          ry(mt+m+(i-n))   = rnysl(i)
          rz(mt+m+(i-n))   = rnzsl(i)
          ty(mt+m+(i-n))   = tmysl(i)
          tz(mt+m+(i-n))   = tmzsl(i)
c         write(82,'(i6,6d15.7)') mt+m+i-n,yc(mt+m+i-n),zc(mt+m+i-n),  
c     #             yb(mt+m+i-n),zb(mt+m+i-n),yf(mt+m+i-n),zf(mt+m+i-n) 
c     #            dphtb(mt+m+i-n), ph(mt+m+i-n), phb(mt+m+i-n)     
c     #          ry(mt+m+i-n),rz(mt+m+i-n),ty(mt+m+i-n),tz(mt+m+i-n)    
        else
          dphtbsl(i)=dphtsl(i)
          dphtb(mb+i-(n+m))= dphtbsl(i)
c          dpht(mb+i-(n+m)) = dphtsl(i)
          ph(mb+i-(n+m))   = phisl(i)
          phb(mb+i-(n+m))   = phisl(i)
c          yc(mb+i-(n+m))   = ycsl(i)
c          zc(mb+i-(n+m))   = zcsl(i)
          yb(mb+i-(n+m))   = ysl(i) 
          zb(mb+i-(n+m))   = zsl(i) 
          yf(mb+i-(n+m))   = ysl(i+1) 
          zf(mb+i-(n+m))   = zsl(i+1) 
          yc(mb+i-(n+m))   = 0.5d0*(yf(mb+i-(n+m))+yb(mb+i-(n+m)))
          zc(mb+i-(n+m))   = 0.5d0*(zf(mb+i-(n+m))+zb(mb+i-(n+m)))
          ry(mb+i-(n+m))   = rnysl(i)
          rz(mb+i-(n+m))   = rnzsl(i)
          ty(mb+i-(n+m))   = tmysl(i)
          tz(mb+i-(n+m))   = tmzsl(i)
c          write(81,'(i6,6d15.7)')
c     #        mb+i-(n+m),yc(mb+i-(n+m)),zc(mb+i-(m+n)),  
c     #      yb(mb+i-(n+m)),zb(mb+i-(m+n)),yf(mb+i-(n+m)),zf(mb+i-(m+n))
c     #       dphtb(mb+i-(m+n)), ph(mb+i-(m+n)), phb(mb+i-(m+n))     
c     #      ry(mb+i-(n+m)),rz(mb+i-(n+m)),ty(mb+i-(n+m)),tz(mb+i-(n+m))
        endif
      enddo
c      write(81,*)
c      write(81,*)
c      write(82,*)
c      write(82,*)
c      write(83,*)
c      write(83,*)
* - ampiezze
      do i=1,ntt
        am(i)=sqrt( (yf(i)-yb(i))**2 + (zf(i)-zb(i))**2 ) 
      enddo
* - getto
c      write(75,*) '# jt ',jt
c      write(76,*) '# jt ',jt
      i0 = ng-(m+n)
      do i=i0-1,ng
        if(i.le.i0)then
          ii = mb+i-i0
          jj = mb-(i-i0)+1 
        elseif(i.ge.i0+1.and.i.le.i0+m)then
          ii = mt+i-i0
          jj = nt-(i-i0)+1
        else
          ii = nt+(i-i0-m)
          jj = ntt-(i-i0-m)+1
        endif
c        hp(ii) = dh(i)
        xi(ii)  = xigb(i)
        xis(ii) = xigs(i)
        ze(ii)  = zegb(i)
        zes(ii) = zegs(i)
c        hp(jj)  = dh(i)
        xi(jj)  = xigf(i)
        xis(jj) = xigs(i)
        ze(jj)  = zegf(i)
        zes(jj) = zegs(i)
c        write(75,'(2i6,4d15.5)') ii,i,xi(ii),xis(ii),
c     #              ze(ii),zes(ii)
c        write(76,'(2i6,4d15.5)') jj,i,xi(jj),xis(jj),
c     #              ze(jj),zes(jj)
      enddo
c      write(75,*)
c      write(75,*)
c      write(75,*) mb-1
c      write(76,*)
c      write(76,*)
* - fattore di matching rl non va con n=0 m>1
      do i=1,nt
        rl(i)=0.d0
      enddo
      do i=nt+1,ntt
        rl(i)=1.d0
      enddo
      do i=1,m
        ii  = mt+i
        if(i.eq.1)then
          dam = 0.5d0*(am(ii)+am(mb))
        else
          dam = 0.5d0*(am(ii)+am(ii-1))
        endif
        rl(ii)   = rl(ii-1) + dam 
      enddo
      if(n.eq.0)then 
        dam = am(mt+m)
      else
        dam = 0.5d0*( am(mt+m) + am(nt+1) )
      endif
      amt = rl(ii) + dam
      do i=1,m
        rl(mt+i)   = rl(mt+i)/amt 
        rl(nt-i+1) = rl(mt+i)
c        write(67,'(3i6,d15.6)') i,mt+i,nt-i+1,rl(mt+i)
      enddo
c      write(67,*)
c      write(67,*)
*      
*
*
* - costituisco la matrice:
*
*1 - inizializzo matrice
      npe = nt+5*(n+m)
      do i=1,npe
        ff(i) = 0.d0
        bb(i) = 0.d0
      do j=1,npe
        aa(i,j) = 0.d0
      enddo
      enddo
*2 - parte BEM, righe =1,nt
* - qui i e indice di riga  NB: ntt = NT negli appunti
      do i = 1,nt
        yy = yc(i)
        zz = zc(i)
        do k = 1,ntt
          yp = yb(k)
          ys = yf(k)
          zp = zb(k)
          zs = zf(k)
          yps= -yf(k)
          yss= -yb(k)
          zps= zf(k)
          zss= zb(k)
* - coeff influenza, pannello + simmetrico rispetto a y=0
          call finte(yy,zz,yp,zp,ys,zs,am(k),fint1,fint2,0.d0)
          call finte(yy,zz,yps,zps,yss,zss,am(k),fins1,fins2,0.d0)
          gik  =  (fint1+fins1)/(2.d0*pi)
          dgik = -(fint2+fins2)/(2.d0*pi)
          if(k.le.mb)then
            aa(i,k) = aa(i,k) + dgik    !/1000.d0
            ff(i)   = ff(i) + dphn(k)*gik
          elseif(k.ge.mb+1.and.k.le.mt)then
            aa(i,k) = aa(i,k) - gik
            ff(i)   = ff(i) - ph(k)*dgik 
          elseif(k.ge.mt+1.and.k.le.mt+m)then
            rlk  = rl(k)
c            hpk  = hp(k)
c            sqh  = sqrt(1.d0+hpk**2)
            xik  = xi(k)
            xisk = xis(k)
            zek  = ze(k)
            zesk = zes(k)
            aik = dgik
            bik = dgik*(xik-xisk) 
            cik = dgik*(zek-zesk) 
            dik = 0.5d0*dgik*( (xik-xisk)**2-(zek-zesk)**2 )  
            eik = dgik*(xik-xisk)*(zek-zesk)  
            aa(i,k)                = aa(i,k)+ (1.d0-rlk)*dgik !/1000.d0 
            aa(i,nt+k-mt)          = aa(i,nt+k-mt)         + rlk*aik
            aa(i,nt+(n+m)+k-mt)    = aa(i,nt+(n+m)+k-mt)   + rlk*bik
            aa(i,nt+2*(n+m)+k-mt)  = aa(i,nt+2*(n+m)+k-mt) + rlk*cik
            aa(i,nt+3*(n+m)+k-mt)  = aa(i,nt+3*(n+m)+k-mt) + rlk*dik
            aa(i,nt+4*(n+m)+k-mt)  = aa(i,nt+4*(n+m)+k-mt) + rlk*eik
            ff(i)  = ff(i) + dphn(k)*gik
          elseif(k.ge.mt+m+1.and.k.le.nt)then
            kk   = mt+(nt-k)+1
            rlk  = rl(k)
c            hpk  = hp(k)
c            sqh  = sqrt(1.d0+hpk**2)
            xik  = xi(k)
            xiskk= xis(kk)
            zek  = ze(k)
            zeskk= zes(kk)
            ryk  = ry(k)
            rzk  = rz(k)
            bik  = ryk*gik
            cik  = rzk*gik
            dik  = gik*( ryk*(xik-xiskk) - rzk*(zek-zeskk) )  
            eik  = gik*( ryk*(zek-zeskk) + rzk*(xik-xiskk) )   
c           bik  = + hpk/sqh*gik
c           cik  = - 1.d0/sqh*gik 
c            dik  = + gik/sqh*( hpk*(xik-xiskk) + (zek-zeskk) )
c            eik  = + gik/sqh*( hpk*(zek-zeskk) - (xik-xiskk) )
            aa(i,nt+(n+m)+kk-mt)    = aa(i,nt+(n+m)+kk-mt)   - rlk*bik
            aa(i,nt+2*(n+m)+kk-mt)  = aa(i,nt+2*(n+m)+kk-mt) - rlk*cik
            aa(i,nt+3*(n+m)+kk-mt)  = aa(i,nt+3*(n+m)+kk-mt) - rlk*dik
            aa(i,nt+4*(n+m)+kk-mt)  = aa(i,nt+4*(n+m)+kk-mt) - rlk*eik
            aa(i,k) = aa(i,k)-(1.d0-rlk)*gik
            ff(i)   = ff(i) - ph(k)*dgik 
          elseif(k.ge.nt+1.and.k.le.nt+n)then
c            hpk  = hp(k)
c            sqh  = sqrt(1.d0+hpk**2)
            xik  = xi(k)
            xisk = xis(k)
            zek  = ze(k)
            zesk = zes(k)
            aik = dgik
            bik = dgik*(xik-xisk)
            cik = dgik*(zek-zesk) 
            dik = 0.5d0*dgik*( (xik-xisk)**2-(zek-zesk)**2 )  
            eik = dgik*(xik-xisk)*(zek-zesk)  
            aa(i,nt+m+k-nt)          = aa(i,nt+m+k-nt)         + aik
            aa(i,nt+(n+m)+m+k-nt)    = aa(i,nt+(n+m)+m+k-nt)   + bik
            aa(i,nt+2*(n+m)+m+k-nt)  = aa(i,nt+2*(n+m)+m+k-nt) + cik
            aa(i,nt+3*(n+m)+m+k-nt)  = aa(i,nt+3*(n+m)+m+k-nt) + dik
            aa(i,nt+4*(n+m)+m+k-nt)  = aa(i,nt+4*(n+m)+m+k-nt) + eik
            ff(i)  = ff(i) + dphn(k)*gik
          elseif(k.ge.nt+n+1)then
            kk   = nt+(ntt-k)+1
c            hpk  = hp(k)
c            sqh  = sqrt(1.d0+hpk**2)
            xik  = xi(k)
            xiskk= xis(kk)
            zek  = ze(k)
            zeskk= zes(kk)
            ryk  = ry(k)
            rzk  = rz(k)
            bik  = ryk*gik
            cik  = rzk*gik
            dik  = gik*( ryk*(xik-xiskk) - rzk*(zek-zeskk) )  
            eik  = gik*( ryk*(zek-zeskk) + rzk*(xik-xiskk) )   
c            bik  = + hpk/sqh*gik
c            cik  = - 1.d0/sqh*gik 
c            dik  = + gik/sqh*( hpk*(xik-xiskk) + (zek-zeskk) )
c            eik  = + gik/sqh*( hpk*(zek-zeskk) - (xik-xiskk) )
            aa(i,nt+(n+m)+m+kk-nt)  = aa(i,nt+(n+m)+m+kk-nt)   - bik
            aa(i,nt+2*(n+m)+m+kk-nt)= aa(i,nt+2*(n+m)+m+kk-nt) - cik
            aa(i,nt+3*(n+m)+m+kk-nt)= aa(i,nt+3*(n+m)+m+kk-nt) - dik
            aa(i,nt+4*(n+m)+m+kk-nt)= aa(i,nt+4*(n+m)+m+kk-nt) - eik
            ff(i)   = ff(i) - ph(k)*dgik 
          endif
        enddo
* - termini autoindotti
        if (i.le.mb) then 
          aa(i,i) = aa(i,i) + 0.5d0  !/1000.d0 
        elseif(i.ge.mb+1.and.i.le.mt)then
          ff(i)   = ff(i) - 0.5d0*ph(i)
        elseif(i.ge.mt+1.and.i.le.mt+m)then
          rli = rl(i)
          xii = xi(i)
          xisi= xis(i)
          zei = ze(i)
          zesi= zes(i)
          aa(i,i)             = aa(i,i) + (1.d0-rli)*0.5d0  !/1000.d0
          aa(i,nt+i-mt)       = aa(i,nt+i-mt)  +  rli*0.5d0 
          aa(i,nt+(m+n)+i-mt) = aa(i,nt+(m+n)+i-mt) + 
     #                                             rli*0.5d0*(xii-xisi) 
          aa(i,nt+2*(m+n)+i-mt)= aa(i,nt+2*(m+n)+i-mt) + 
     #                                            rli*0.5d0*(zei-zesi) 
          aa(i,nt+3*(m+n)+i-mt)= aa(i,nt+3*(m+n)+i-mt) + 
     #                          rli*0.25d0*((xii-xisi)**2-(zei-zesi)**2)
          aa(i,nt+4*(m+n)+i-mt)=aa(i,nt+4*(m+n)+i-mt) + 
     #                          rli*0.5d0*(xii-xisi)*(zei-zesi)
        elseif(i.ge.mt+m+1.and.i.le.nt)then
          ff(i)   = ff(i) - 0.5d0*ph(i)
        endif 
      enddo
*
*3 - parte getto ,  righe =nt+1,nt+5(m+n)
* - NB : qui i non e indice di riga, bensi indice di pannello (cella) del getto
*
      do 100, i = mt+1, nt+n 
*
        if(i.ge.mt+1.and.i.le.mt+m)then 
          if(i.eq.mt+1)then
            i1 = mb
            i2 = mb-1
            j  = mb+1
            j1 = nt
          elseif(i.eq.mt+2)then
            i1 = i-1
            i2 = mb
            j  = nt-(i-mt)+2
            j1 = j-1
          else
            i1 = i-1
            i2 = i-2
            j  = nt-(i-mt)+2
            j1 = j-1
          endif
          ii  = nt+i-mt
          iii = ii
          iij = ii
          mn  = m+n
          k1  = m
*
        elseif(i.ge.nt+1.and.i.le.nt+n)then
*
          if(i.eq.nt+1)then
            if(m.eq.0)then
              i1 = mb
              i2 = mb-1
              j  = mb+1
              j1 = ntt
            else
              i1 = mt+m
              i2 = mt+m-1
              j  = mt+m+1
              j1 = ntt
            endif
          elseif(i.eq.nt+2)then
            if(m.eq.0)then
              i1 = i-1
              i2 = mb
              j  = ntt-(i-nt)+2
              j1 = j-1
            else
              i1 = i-1
              i2 = mt+m
              j  = ntt-(i-nt)+2
              j1 = j-1
            endif
          else
            i1 = i-1
            i2 = i-2 
            j  = ntt-(i-nt)+2
            j1 = j-1
          endif
          ii  = nt+i-nt
          iii = ii+5*m
          iij = ii+m
          mn  = m+n
          k1  = n
*
        else 
*
          goto 100
        endif
*
        rli  = rl(i)
        rli1 = rl(i1)
        drl  = rli1-rli
        dxif = (xi(i)-xi(i1))     !*1000.d0
        dxib = (xi(i1)-xi(i2))    !*1000.d0
        xii  = xi(i)
        xii1 = xi(i1)
        xisi = xis(i)
        xisi1= xis(i1) 
        zei  = ze(i)
        zei1 = ze(i1)
        zesi = zes(i)
        zesi1= zes(i1) 
        xij  = xi(j)
        xij1 = xi(j1)
        zej  = ze(j)
        zej1 = ze(j1)
c        sqh  = sqrt(1.d0+hp(i1)**2)
c        hpi1 = hp(i1)
        ryj  = ry(j)
        rzj  = rz(j)
        ryi1 = ry(i1)
        rzi1 = rz(i1)
        ryi  = ry(i)
        rzi  = rz(i)
c        write(99,'(''drl '',3i6,7d15.6)')i,i1,i2,rli,
c     #                      ryj,rzj,ryi1,rzi1,ryi,rzi
c        write(98,'(''drl '',3i6,5d15.6)')i,i1,i2,drl,zesi,zesi1,zei,zei1
c        write(97,'(''drl '',3i6,5d15.6)')i,i1,i2,drl,xisi,xisi1,xij,xij1
c        write(96,'(''drl '',3i6,5d15.6)')i,i1,i2,drl,zesi,zesi1,zej,zej1
* - continuita phi sul corpo (i1)
        aa(iii,i1)         = aa(iii,i1)         + drl
        aa(iii,iij)        = aa(iii,iij)        + rli  
        aa(iii,iij-1)      = aa(iii,iij-1)      - rli1 
        aa(iii,iij+mn)     = aa(iii,iij+mn)     + rli*(xii1-xisi) 
        aa(iii,iij-1+mn)   = aa(iii,iij-1+mn)   - rli1*(xii1-xisi1) 
        aa(iii,iij+2*mn)   = aa(iii,iij+2*mn)   + rli*(zei1-zesi) 
        aa(iii,iij-1+2*mn) = aa(iii,iij-1+2*mn) - rli1*(zei1-zesi1) 
        aa(iii,iij+3*mn)   = aa(iii,iij+3*mn)   + 
     #                       0.5d0*rli*((xii1-xisi)**2-(zei1-zesi)**2)
        aa(iii,iij-1+3*mn) = aa(iii,iij-1+3*mn) - 
     #                      0.5d0*rli1*((xii1-xisi1)**2-(zei1-zesi1)**2)
        aa(iii,iij+4*mn)   = aa(iii,iij+4*mn)   +
     #                      rli*(xii1-xisi)*(zei1-zesi)
        aa(iii,iij-1+4*mn) = aa(iii,iij-1+4*mn) - 
     #                      rli1*(xii1-xisi1)*(zei1-zesi1)
* - continuita vn sulla SL (j)
c        aa(iii,j)         = aa(iii,j)         + drl           
c        aa(iii,iij+mn)    = aa(iii,iij+mn)    + rli*ryj 
c        aa(iii,iij-1+mn)  = aa(iii,iij-1+mn)  - rli1*ryj 
c        aa(iii,iij+2*mn)  = aa(iii,iij+2*mn)  + rli*rzj 
c        aa(iii,iij-1+2*mn)= aa(iii,iij-1+2*mn)- rli1*rzj 
c        aa(iii,iij+3*mn)  = aa(iii,iij+3*mn)  +
c     #                         rli*( ryj*(xij-xisi) - rzj*(zej-zesi) ) 
c        aa(iii,iij-1+3*mn)= aa(iii,iij-1+3*mn)- 
c     #                        rli1*( ryj*(xij-xisi1) - rzj*(zej-zesi1) )
c        aa(iii,iij+4*mn)  = aa(iii,iij+4*mn)  +
c     #                         rli*( ryj*(zej-zesi) + rzj*(xij-xisi) ) 
c        aa(iii,iij-1+4*mn)= aa(iii,iij-1+4*mn)- 
c     #                        rli1*( ryj*(zej-zesi1)+ rzj*(xij-xisi1) )
*- continuita vtau sul corpo (i1)  NON FUNZIONA
c        if(i.eq.nt+1)then
c          aa(iii,i2)         = aa(iii,i2)         + 0.5d0/dxib
c          aa(iii,i1)         = aa(iii,i1)         +
c     #                                  0.5d0/dxif-0.5d0/dxib
c          aa(iii,iij)       = aa(iii,iij)       - 0.5d0/dxif 
c          aa(iii,iij+mn)    = aa(iii,iij+mn)    - 0.5d0/dxif*(xii1-xisi)
c          aa(iii,iij+2*mn)  = aa(iii,iij+2*mn)  - 0.5d0/dxif*(zei1-zesi)
c          aa(iii,iij+3*mn)  = aa(iii,iij+3*mn)  -
c     #                  0.5d0/dxif*0.5d0*((xii1-xisi)**2-(zei1-zesi)**2)
c          aa(iii,iij+4*mn)  = aa(iii,iij+4*mn)  -
c     #                      0.5d0/dxif*(xii1-xisi)*(zei1-zesi)
c          aa(iii,iij+mn)      = aa(iii,iij+mn)      + rli 
c          aa(iii,iij+3*mn)    = aa(iii,iij+3*mn)    + rli*(xii1-xisi) 
c          aa(iii,iij+4*mn)    = aa(iii,iij+4*mn)    + rli*(zei1-zesi) 
c          aa(iii,i2)         = aa(iii,i2)         + 0.5d0/dxib
c          aa(iii,i1)         = aa(iii,i1)         - 0.5d0/dxib 
c        elseif(i.ge.mt+1.and.i.le.mt+m)then
c          aa(iii,i)          = aa(iii,i)          - 0.5d0/dxif
c        else 
c          aa(iii,iij+mn)      = aa(iii,iij+mn)      + rli 
c          aa(iii,iij+3*mn)    = aa(iii,iij+3*mn)    + rli*(xii1-xisi) 
c          aa(iii,iij+4*mn)    = aa(iii,iij+4*mn)    + rli*(zei1-zesi) 
c          aa(iii,iij-1+mn)    = aa(iii,iij-1+mn)    - rli1 
c          aa(iii,iij-1+3*mn)  = aa(iii,iij-1+3*mn)  - rli1*(xii1-xisi1) 
c          aa(iii,iij-1+4*mn)  = aa(iii,iij-1+4*mn)  - rli1*(zei1-zesi1) 
c        endif
* - impongo vtau sulla SL (j)
c        aa(iii+k1,iij+mn)    = aa(iii+k1,iij+mn)    - 1.d0/sqh 
c        aa(iii+k1,iij+2*mn)  = aa(iii+k1,iij+2*mn)  - hpi1/sqh 
c        aa(iii+k1,iij+3*mn)  = aa(iii+k1,iij+3*mn)  -
c     #                          1.d0/sqh*( (xij-xisi)-hpi1*(zej-zesi) ) 
c        aa(iii+k1,iij+4*mn)  = aa(iii+k1,iij+4*mn)  -
c     #                          1.d0/sqh*( (zej-zesi)+hpi1*(xij-xisi) )
* - continuita phi sulla SL (j)
        aa(iii+k1,iij)       = aa(iii+k1,iij)        + rli  
        aa(iii+k1,iij-1)     = aa(iii+k1,iij-1)      - rli1 
        aa(iii+k1,iij+mn)    =aa(iii+k1,iij+mn)    +rli*(xij-xisi) 
        aa(iii+k1,iij-1+mn)  =aa(iii+k1,iij-1+mn)  -rli1*(xij-xisi1)
        aa(iii+k1,iij+2*mn)  =aa(iii+k1,iij+2*mn)  +rli*(zej-zesi) 
        aa(iii+k1,iij-1+2*mn)=aa(iii+k1,iij-1+2*mn)-rli1*(zej-zesi1)
        aa(iii+k1,iij+3*mn)  =aa(iii+k1,iij+3*mn)  + 
     #                       0.5d0*rli*((xij-xisi)**2-(zej-zesi)**2)
        aa(iii+k1,iij-1+3*mn) = aa(iii+k1,iij-1+3*mn) - 
     #                      0.5d0*rli1*((xij-xisi1)**2-(zej-zesi1)**2)
        aa(iii+k1,iij+4*mn)   = aa(iii+k1,iij+4*mn)   +
     #                           rli*(xij-xisi)*(zej-zesi)
        aa(iii+k1,iij-1+4*mn) = aa(iii+k1,iij-1+4*mn) - 
     #                           rli1*(xij-xisi1)*(zej-zesi1)
* 
* - phi sulla SL  (j)
c        aa(iii+k1,iij)       = aa(iii+k1,iij)       + 1.d0 
c        aa(iii+k1,iij+mn)    = aa(iii+k1,iij+mn)    + (xij-xisi)
c        aa(iii+k1,iij+2*mn)  = aa(iii+k1,iij+2*mn)  + (zej-zesi)
c        aa(iii+k1,iij+3*mn)  = aa(iii+k1,iij+3*mn)  +
c     #                             0.5d0*((xij-xisi)**2-(zej-zesi)**2)
c        aa(iii+k1,iij+4*mn)  = aa(iii+k1,iij+4*mn)  +
c     #                             (xij-xisi)*(zej-zesi)
* - vn sul corpo (i1)
        aa(iii+2*k1,iij+mn)   = aa(iii+2*k1,iij+mn)    + ryi1 
        aa(iii+2*k1,iij+2*mn) = aa(iii+2*k1,iij+2*mn)  + rzi1 
        aa(iii+2*k1,iij+3*mn) = aa(iii+2*k1,iij+3*mn)  +
     #                            ryi1*(xii1-xisi) - rzi1*(zei1-zesi)
        aa(iii+2*k1,iij+4*mn) = aa(iii+2*k1,iij+4*mn)  +
     #                            ryi1*(zei1-zesi) + rzi1*(xii1-xisi)
* - vn sul corpo (i)
        aa(iii+3*k1,iij+mn)   = aa(iii+3*k1,iij+mn)    + ryi 
        aa(iii+3*k1,iij+2*mn) = aa(iii+3*k1,iij+2*mn)  + rzi 
        aa(iii+3*k1,iij+3*mn) = aa(iii+3*k1,iij+3*mn)  +
     #                            ryi*(xii-xisi) - rzi*(zei-zesi)
        aa(iii+3*k1,iij+4*mn) = aa(iii+3*k1,iij+4*mn)  +
     #                            ryi*(zei-zesi) + rzi*(xii-xisi)
* - phi sulla SL (j1)
        aa(iii+4*k1,iij)       = aa(iii+4*k1,iij)       + 1.d0 
        aa(iii+4*k1,iij+mn)    = aa(iii+4*k1,iij+mn)    + (xij1-xisi)
        aa(iii+4*k1,iij+2*mn)  = aa(iii+4*k1,iij+2*mn)  + (zej1-zesi)
        aa(iii+4*k1,iij+3*mn)  = aa(iii+4*k1,iij+3*mn)  +
     #                             0.5d0*((xij1-xisi)**2-(zej1-zesi)**2)
        aa(iii+4*k1,iij+4*mn)  = aa(iii+4*k1,iij+4*mn)  +
     #                             (xij1-xisi)*(zej1-zesi)
*
c         if(i.eq.nt+1) ff(iii) = ff(iii) + 0.2429340d-02
c         if(i.eq.nt+1) ff(iii)=ff(iii)+1.94288d0
         ff(iii)       = ff(iii)      + 0.0d0
c        ff(iii+ k1 )  = ff(iii+k1)   + 0.d0 
c        ff(iii+ k1 )  = ff(iii+k1)   + dphtb(j)
        ff(iii+ k1 )  = ff(iii+k1)   - drl*ph(j)
c        write(99,'(3i6,3d15.6)') i,j,iii+k1,ff(iii+k1),ph(j),drl
c        ff(iii+k1)  = ff(iii+k1) + ph(j)       
c        ff(iii+2*k1)  = ff(iii+2*k1) + 0.d0 
        ff(iii+2*k1)  = ff(iii+2*k1) + dphn(i1)
        ff(iii+3*k1)  = ff(iii+3*k1) + dphn(i)
        ff(iii+4*k1)  = ff(iii+4*k1) + ph(j1)       
 100  continue
*    
      do i=1,npe
        bb(i) = ff(i)
      enddo
*
c      do i=1,npe
c      sss = 0.d0
c      do j=1,npe
c        sss = max(sss,abs(aa(i,j)))
c        write(90,'(2i6,d15.7)') i,j,aa(i,j)
c      enddo
c      write(90,*)
c      write(90,*)
c      write(91,'(i6,3d15.6)') i,yc(i),bb(i),sss
c      enddo
c      write(91,*)
c      write(91,*)
*
c ----------------- SOLUZIONE CON ROUTINE GAUSS ------------------------
c      do i=1,npe
c      aa(i,npe+1)=bb(i)
c      end do
c      call gauss(aa,npamx,npe,npe+1,1.d-12)
c      do i=1,npe
c      bb(i)=aa(i,npe+1)
c      end do 



c
c     -------------   SOLUZIONE SISTEMA LINEARE CON ESSL   -------------

c      call dgefcd(aa,npamx,npe,indx,1,rcond,det,appo,npamx)
c      call dges(aa,npamx,npe,indx,bb,0) 

c     -------------   SOLUZIONE SISTEMA LINEARE CON LAPACK -------------

      call DGETRF(npe,npe,aa,npamx,indx,idum)
      trans='N'
c      rmax = 0.d0
c      rmin = 1.d31
c      do i=1,npe
c        rmax = max(rmax,abs(aa(i,i)))
c        rmin = min(rmin,abs(aa(i,i)))
c      enddo
c      write(*,*) 'rr ',rmax,rmin
      call DGETRS(trans,npe,1,aa,npamx,indx,bb,npamx,idum)

c     -------------   FINE SOLUZIONE SISTEMA LINEARE       -------------

*    risistemo arrays dopo soluzione sistema lineare
*
      do i=1,mb
        ph(i) = bb(i)
        phb(i) = bb(i)
c        write(76,'(i6,4d15.7)') i,yc(i),zc(i),ph(i),dphn(i)
      enddo
c        write(76,*) 
      do i=mb+1,mb+mf
        dphn(i) = bb(i)
        dphnb(i) = bb(i)
c        write(77,'(i6,4d15.7)') i,yc(i),zc(i),dphn(i),ph(i)
      enddo
      do i=mt+1,mt+m
        phb(i) = bb(i)
c        write(76,'(i6,4d15.7)') i,yc(i),zc(i),phb(i),dphn(i)
      enddo 
c        write(77,*) 
      do i=mt+m+1,nt
        dphnb(i) = bb(i)
c        write(77,'(i6,4d15.7)') i,yc(i),zc(i),dphnb(i),ph(i)
      enddo
c      write(76,*) 
c      write(76,*) 
c      write(77,*) 
c      write(77,*) 
c      do i =nt+1,nt+(m+n)
c        ii = i-nt
c        a(ii) = bb(i)
c        b(ii) = bb(i+m+n)
c        c(ii) = bb(i+2*(m+n))
c        d(ii) = bb(i+3*(m+n))
c        e(ii) = bb(i+4*(m+n))
c        write(80,'(i6,5d25.17)') ii,a(ii),b(ii),c(ii),d(ii),e(ii)
c      enddo 
c      write(80,*) '# jt ',jt
      do i = mt+1,mt+m
        j  = nt-(i-mt)+1
        ii = i-mt 
        a(i) = bb(nt+ii)
        b(i) = bb(nt+(m+n)+ii)
        c(i) = bb(nt+2*(m+n)+ii)
        d(i) = bb(nt+3*(m+n)+ii)
        e(i) = bb(nt+4*(m+n)+ii)
        a(j) = a(i)
        b(j) = b(i)
        c(j) = c(i)
        d(j) = d(i)
        e(j) = e(i)
c        write(80,'(i6,5d15.7)') i,a(i),b(i),c(i),d(i),e(i)
      enddo
c      write(80,*)
      do i = nt+1,nt+n
        j  = ntt-(i-nt)+1
        ii = i-nt+m 
        a(i) = bb(nt+ii)
        b(i) = bb(nt+(m+n)+ii)
        c(i) = bb(nt+2*(m+n)+ii)
        d(i) = bb(nt+3*(m+n)+ii)
        e(i) = bb(nt+4*(m+n)+ii)
        a(j) = a(i)
        b(j) = b(i)
        c(j) = c(i)
        d(j) = d(i)
        e(j) = e(i)
c        write(80,'(i6,5d15.7)') i,a(i),b(i),c(i),d(i),e(i)
      enddo
c      write(80,*)
c      write(80,*) 
      
* riposiziono negli array d uscita 
c      do i = 1,mb+m
c        if(i.le.mb)then
c          phi(i) = ph(i)
c          phib(i) = ph(i)
c        elseif(i.ge.mb+1.and.i.le.mb+m)then
c          phib(i) = phb(mt+i-mb)  
c        endif  
c      enddo    
c      do i = n+1,mf+m+n
c        if(i.ge.n+1.and.i.le.n+m)then
c          dphibsl(i) = dphnb(mt+m+(i-n))    
c        else
c          dphibsl(i)  = dphn(mb+i-(n+m))    
c          dphisl(i)  = dphn(mb+i-(n+m))    
c        endif
c      enddo

ccc* -  scrivo un po di controlli
ccc      do i=1,m
ccc        if(i.eq.1)then
ccc         ii1 = mb
ccc         jj  = mt+1
ccc        else
ccc         ii1 = mt+i-1
ccc         jj  = nt-i+2
ccc        endif
ccc        ii = mt+i
ccc        dxi  = xi(ii1)-xis(ii)
ccc        dze  = ze(ii1)-zes(ii)
ccc        phig = a(ii)+b(ii)*dxi+c(ii)*dze+0.5d0*d(ii)*(dxi**2-dze**2)+
ccc     #         e(ii)*dxi*dze 
cccc        write(92,'(4d15.7)') yc(ii1),zc(ii1),phb(ii1),phig 
ccc        dxi  = xi(jj)-xis(ii)
ccc        dze  = ze(jj)-zes(ii)
ccc        phig = a(ii)+b(ii)*dxi+c(ii)*dze+0.5d0*d(ii)*(dxi**2-dze**2)+
ccc     #         e(ii)*dxi*dze 
cccc        write(93,'(4d15.7)') yc(jj ),zc(jj ),ph(jj ),phig 
ccc       enddo         
ccc       do i=1,m+n
ccc        if(i.eq.1.and.m.eq.0)then
ccc         ii1 = mb
ccc         jj  = mb+1
ccc         jj1 = ntt
ccc         ii  = nt+i
ccc         i1  = i-1
ccc        elseif(i.eq.1.and.m.ne.0)then
ccc         ii1 = mb
ccc         jj  = mb+1
ccc         jj1 = nt 
ccc         ii  = mt+i
ccc         i1  = i-1
ccc        elseif(i.ge.2.and.i.le.m)then
ccc         ii  = mt+i
ccc         ii1 = ii-1 
ccc         jj  = nt-i+2
ccc         jj1 = jj-1 
ccc         i1  = i-1
ccc        else
ccc         ii  = nt+i-m
ccc         ii1 = ii-1 
ccc         jj  = ntt-(i-m)+2
ccc         jj1 = jj-1 
ccc         i1  = i-1
ccc        endif
ccc        sqh=sqrt(1.d0+hp(ii)**2)
ccc        rlii = rl(ii)
ccc        a1 =  a(ii)
ccc        b1 =  b(ii)
ccc        c1 =  c(ii)
ccc        d1 =  d(ii)
ccc        e1 =  e(ii)
ccc        dxi = xi(ii)-xis(ii)
ccc        dze = ze(ii)-zes(ii)
ccc        phig = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
ccc        vxig  = b1+d1*dxi+e1*dze
ccc        vzeg  = c1-d1*dze+e1*dxi
ccc        dxi = xi(jj1)-xis(ii)
ccc        dze = ze(jj1)-zes(ii)
ccc        phif = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
ccc        vxif  = b1+d1*dxi+e1*dze
ccc        vzef  = c1-d1*dze+e1*dxi
ccc        vnf   = (vxif*hp(ii)-vzef)/sqh
ccc        vtf   = (-vxif-hp(ii)*vzef)/sqh
ccc        a1 =  a(ii+1)
ccc        b1 =  b(ii+1)
ccc        c1 =  c(ii+1)
ccc        d1 =  d(ii+1)
ccc        e1 =  e(ii+1)
ccc        dxi = xi(ii)-xis(ii+1)
ccc        dze = ze(ii)-zes(ii+1)
ccc        phigg = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
ccc        vxigg = b1+d1*dxi+e1*dze
ccc        vzegg = c1-d1*dze+e1*dxi
ccc        dxi = xi(jj1)-xis(ii+1)
ccc        dze = ze(jj1)-zes(ii+1)
ccc        phiff = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
ccc        vxiff = b1+d1*dxi+e1*dze
ccc        vzeff = c1-d1*dze+e1*dxi
ccc        vnff  = (vxiff*hp(ii)-vzeff)/sqh
ccc        vtff  = (-vxiff-hp(ii)*vzeff)/sqh
ccc        phii = (1.d0-rlii)*phb(ii)+rlii*phig
ccc        dphii= (1.d0-rlii)*dphnb(jj1)+rlii*vnff
ccc        write(62,'(6d15.6)')   yc(ii),zc(ii), phb(ii),phig,phigg,phii
ccc        write(64,'(8d15.7)') 
ccc     #      yc(ii),zc(ii), dphtb(ii),vxig,vxigg,dphn(ii),vzeg,vzegg
ccc        write(63,'(5d15.6)') yc(jj1),zc(jj1),ph(jj1),phif,phiff
ccc        write(65,'(9d15.7)') 
ccc     #      yc(jj1),zc(jj1), dphtb(jj1),vtf,vtff,dphnb(jj1),vnf,vnff,
ccc     #      dphii
cccc     #      yc(jj1),zc(jj1), dphtb(jj1),vxif,vxiff,dphn(jj1),vzef,vzeff
ccc      enddo
ccc      write(62,*)
ccc      write(62,*)
ccc      write(63,*)
ccc      write(63,*)
ccc      write(64,*)
ccc      write(64,*)
ccc      write(65,*)
ccc      write(65,*)
ccc            
ccc*
      return
      end
