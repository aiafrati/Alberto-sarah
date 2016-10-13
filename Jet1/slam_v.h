c     --  file da includere nel programma slam.f  --

c     - contiene le aree common da utilizzare nel programma 

      common/char/scon,svel,sprt,sadi,spot,sfor,sad2
*
      common/varint/npt,npc,lrest,llf,ist,isp,kcut,npci,ksta,kitm,ng,
     &              iint,kini,ngo,iwig
*
      common/varrea/dt,ampp,pfraz,estr,alfa,ande,vfall0,t,tend,tsta,
     &              tust,frdt,escr,pro0,dtpa,ampk,
     &              angma,cmm0,ancut,spo0,rapma,amplim,rappi,
     &              weir,vfall,rmu,tfor,rga,vea,rkk,dzet,vsopr,
     &              volb0,volj,dvolab,dvolabo,volb1,volb1o,t0,ggr,
     &              xs0,ux,al,bl,tl,alpha, x00,wl,dztr,ztr,dx0,zmaxs
*
      common/vecint/kphi(npamx+1)
*
      common/vecrea/yv(npamx+1)  , zv(npamx+1) ,
     &              yn(npamx+1)  , zn(npamx+1) ,
     &              amp(npamx)   , dphi(npamx) ,
     &              phi(npamx)   , phin(npamx) ,
     &              depn(npamx,2),
     &              yce(npamx)   , zce(npamx)  ,
     &              ycn(npamx)   , zcn(npamx)  ,
     &              vym(npamx,2) , vzm(npamx,2),
     &              dpnt(npamx)  , dpt(npamx)  ,
     &              pre1(npamx)  , pre2(npamx) ,
     &              pred1(npamx) , pred2(npamx),
     &              tt(ntmx)     , fosl(ntmx)  ,
     &              voat(ntmx)   , voso(ntmx)  ,
     &              foel(ntmx)   ,
     &              yg(npamx+1)  , zg(npamx+1) , ta(npamx),
     &              ygn(npamx+1) , zgn(npamx+1),
     &              ygs2(npamx+1),zgs2(npamx+1),
     &              ysl(npamx)   , zsl(npamx)  , phisl(npamx),
     &              ycsl(npamx)  , zcsl(npamx) ,
     &              pre12(npamx),dpnt2(npamx),dpt2(npamx), 
     &              tag(npamx)   , rnxg(npamx),rnxgf(npamx),
     &              rnxgs2(npamx), rnxx(npamx)
c-----------------------------------------------------------------
