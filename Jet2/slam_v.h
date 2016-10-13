c     --  file da includere nel programma slam.f  --

c     - contiene le aree common da utilizzare nel programma 

      common/char/scon,svel,sprt,sadi,spot,sfor

      common/varint/npt,npc,lrest,llf,ist,isp,kcut,npci,ksta,kitm

      common/varrea/dt,ampp,pfraz,estr,alfa,ande,vfall0,t,tend,tsta,
     &              tust,frdt,escr,pro0,dtpa,ampk,
     &              angma,cmm0,ancut,spo0,rapma,amplim,rappi,
     &              weir,vfall,rmu,tfor,rga,vea,rkk,dzet,vsopr

      common/vecint/kphi(npamx+1)

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
     &              foel(ntmx)
c-----------------------------------------------------------------
