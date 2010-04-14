      implicit none
      double precision x,d,gradx,fx,xnorm,dnorm,sig,x0,d0,gradx0,
     1       fx0,x0norm, d0norm,sig0,difx,ftest,dirder,fxbest
      logical logit
      common/domin_xdat/x(NX),d(NX),gradx(NX),fx,xnorm,dnorm,sig,
     1          x0(NX),d0(NX),gradx0(NX),fx0,x0norm,d0norm,sig0,
     2          difx(NX),ftest,dirder,fxbest,logit


      double precision stepterm,sigsm,sigla,alpha,beta,delta,theta,dskal
      common/domin_stepparam/stepterm,sigsm,sigla,alpha,beta,delta,
     1               theta,dskal


      double precision epsx,epsg,epsf
      logical test1,test2,test3,cold,test4
      common/domin_param/epsx,epsg,epsf,test1,test2,test3,test4,cold


      double precision epsmach,tolmach,deldif
      common/domin_mach/epsmach,tolmach,deldif


      integer n
      common/domin_dim/n


      integer icf,icgf,iup,inoup,ires,imod
      common/domin_count/icf,icgf,iup,inoup,ires,imod


      double precision a,hessian,a0,diaga0,diaga,z,y,tracea,
     1                 traceb,accinf
      common/domin_bfgsdat/a(NX,NX),hessian(NX,NX),
     1               diaga(NX),a0(NSMALL,NSMALL),
     2               diaga0(NSMALL),accinf(MAXIT,14),z(NX),y(NX),
     3               tracea,traceb


      logical analyt
      common/domin_gradinf/analyt


      double precision eta,thet2,bfgsterm
      integer niter,itermax
      common/domin_bfgsparam/niter,itermax,eta,thet2,bfgsterm


      double precision xst
      common/domin_xst/xst(NX)


      character*32 ident
      common/domin_ident/ident


      integer strlen
      external strlen

