c
c     made some changes to support openmp
c     May 2007, H.Rue
c      
#include "dominparam.inc"
      integer function domin_fileid()
c     return the fileid for each thread
      integer omp_get_thread_num
      external omp_get_thread_num
      domin_fileid = 10+omp_get_thread_num()
      end
c
      subroutine domin_seteps(u_epsx,u_epsf,u_epsg)
#include "domincommon.inc"
      double precision u_epsx,u_epsf,u_epsg
      if (u_epsx .ge. 0.0) then
         epsx = u_epsx
      endif
      if (u_epsf .ge. 0.0) then
         epsf = u_epsf
      endif
      if (u_epsg .ge. 0.0) then
         epsg = u_epsg
      endif
      end
c
      subroutine domin
c************************************************************
c     initialization program for bfgs-minimization
c************************************************************
#include "domincommon.inc"
      character*2 fileidtxt
      integer domin_fileid
      external domin_fileid
      double precision term,tol1
      integer i,j,k
      tol1=99999999.0
c     we compute the machine precision first
      epsmach = 2.d0**(-20)
      term=1.d0+epsmach
      do while ( term .ne. 1.d0 )
         epsmach=epsmach/2.d0
         term=1.d0+epsmach
      enddo
      epsmach=epsmach+epsmach
      tolmach=epsmach
      do while ( tolmach .ne. 0.d0 )
         tol1=tolmach
         tolmach=tolmach/16.d0
      enddo
      tolmach=tol1
c     accinf contains information to be printed in the short protocol
      do i=1,MAXIT
         do j=1,14
            accinf(i,j)=0.d0
         enddo
      enddo
      test1=.false.
      test2=.false.
      test3=.false.
      test4=.false.
      epsx=epsmach**(2.d0/3.d0)
      epsg=epsx
      epsf=epsmach*1.d3
      epsf=dsqrt(epsmach)
      eta=epsmach**.5d0
      sigsm=epsx
      sigla=2048.d0
      beta=2.d0
      alpha=.1d0
      delta=.01d0
      theta=.99d0
      deldif=dexp(dlog(epsmach)/7.d0)*.25d0
      itermax=MAXIT
      call domin_setup
      if (strlen(ident) .eq. 0) then
         logit = .false.
      else
         write(fileidtxt,'(I2.2)') domin_fileid() -10
         open(domin_fileid(),
     !        file='domin'//ident//'-'//fileidtxt//'.res',
     !        status='unknown')
         logit = .true.
      endif
      icf=0
      icgf=0
      iup=0
      inoup=0
      ires=0
      imod=0
c***********************
      call domin_bfgs(cold)
c     start may change control parameters
c     param cold for bfgs means fresh start
      if ( .not. analyt )      icf=icf-n*icgf*6
      if (logit) then
         write(domin_fileid(),fmt='(/1x,a40,
     1        /1x,''termination reason '',i4,/1x,''fxopt = '',d21.14)')
     2        ident,int(bfgsterm),fx
         do i=1,n
            write(domin_fileid(),fmt='(1x,''x('',i3,'')='',d21.14,
     1           '' grad('',i3,'')='',d21.14)') i,x(i),i,gradx(i)
         enddo
      endif
c     anzahl der berechneten funktionswerte ohne
c     die fuer die numerische differentiation benoetigten
      if (logit) then
         write(domin_fileid(),3000) icf,icgf
         write(domin_fileid(),
     f        fmt='(1x,''last estimate of condition number of'',
     f        '' updated hessian is '',d15.7)') accinf(niter,8)
         write(domin_fileid(), fmt='(1x,''steps='',i5,'' updates='',i5,
     f        '' noupdates='',i5,
     f        '' restarts='',i5,'' modsteps='',i5)') 
     f        niter,iup,inoup,ires,imod
      endif
      if ( bfgsterm .lt. 0.d0 ) test1=.true.
      if ( test1 ) then
         if (logit) then
            write(domin_fileid(),fmt='('' bfgs short protocol '')')
            do j=1,niter
               write(domin_fileid(),2000) (accinf(j,k),k=1,14)
            enddo
         endif
      endif
c     ausgabe der kurzinformation
c     it = schrittnummer , f = funktionswert, delx = laenge der aende-
c     rung in x , ngrad = norm des gradienten, (nicht monoton) ,
c     n(a)= norm der quasi-newton-matrix, entspricht norm der hesse-
c     matrix, thet2 klein gegen 1 bedeutet, dass die funktion sich
c     stark nichtkonvex verhaelt , resultiert in verringeter effizienz
c     sig = schrittweite x(neu)=x(alt)-sig*richtung (d),
c     sig klein gegen eins bedeutet, dass man sich noch in grosser
c     entfernung von einem lokalen minimum befindet, dskal ungl 1
c     wenn die berechnete abstiegsrichtung d zu gross ist. 
c     verhindert unnoetrige probleme mit grossen argumenten x bei step
c     n(d)=norm korrektur d, icf+ = anzahl f-auswertungen fuer diesen
c     schritt. in gutartigen faellen 1,2. icf+ gross bedeutet, dass
c     d schlecht, d.h. a schlecht (eventuell restart versuchen)
c     oder x noch sehr ungenau. restart = 1 bedeutet, dass a=einheitsmatrix
c     gesetzt wurde. update = -1 bedeutet, dass die quasi-newton-korrektur
c     unterdrueckt wurde, d.h. die funktion verhaelt sich start nicht-
c     konvex
c     abbruch : bfgsterm = 0 , 1 erwuenschte genauigkeit erreicht
c     =2 richtungsableitung wird extrem klein
c     in der regel optimal erreichbare genauigkeit
c     erzielt, aber epsx , epsg zu klein vorgegeben
c     (tritt ein bei schlecht konditionierten problemen,
c     wo die parameter stark von kleinen aenderungen der
c     daten abhaengen)
c     -2 itermax schritte erreicht, ohne dass genauigkeit 
c     erreicht,-1 bedeutet, dass trotz restart mit der 
c     einheitsmatrix schrittweiten zu klein. tritt nur
c     ein bei extrem schlecht konditionierten problemen,
c     wo die "breite" eines "tales" schmaler als sigsm ist.
c     -3 nichtbehebbare fehler bei aufruf f oder gradf
      if (logit) then
         close(domin_fileid())
      endif
      return
 2000 format(
     1     (/1x,'it=',d10.4,'     f=',d10.4,'  delx=',d10.4,
     2     ' ngrad=',d10.4,'  n(a)=',d10.4,' thet2=',d10.4
     3     /14x,'dirder=',d10.4,'  cond=',d10.4
     4     /14x,'   sig=',d10.4,' dskal=',d10.4,'  n(d)=',
     5     d10.4,'  icf+=',d10.4,
     6     /13x,'restart=',d8.1,'  update=',d8.1))
 3000 format(/,' function calls',i11,' gradient calls  ',
     1     i10)
      end
c****************************************************************
c     determination of stepsize by an armijo-like test for 
c     descent
c****************************************************************
      subroutine domin_step
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
c********************************************************************
c     n = tatsaechliche dimension
c     x = laufender punkt
c     d = richtung
c     gradx = gradient von f bei x
c     x,d,gradx = eingabe
c     x0,d0,gradx0 etc. information aus schritt zuvor (wird automatisch
c     umgespeichert)
c     xnorm,dnorm = euklidische laenge von x bzw. d
c     stepterm = 1 bei erfolg , =-1 sonst 
c     sig = errechnete schrittweit
c     es wird davon ausgegangen, dass 1 asymptotisch optimale schrittweite ist
c     andernfalls muss sig anders initialisiert werden 
c     sigsm = minimale zulaessige schrittweite
c     sigla = maximale zulaessige schrittweite
c     alpha = maximal zulaessige schrittweitenreduktion pro teilung
c     delta = parameter im abstiegstest
c     beta = max. zugel. wachstumsfaktor fuer xnorm bei schrittweite 1
c     theta = schranke fuer cos(winkel(laufende richtung,vorhergehende
c     richtung). falls ueberschritten, wird schrittweite groesser
c     1 versucht
c*********************************************************************
      integer i
      double precision x1(NX),domin_vecnor,t,sigmin,fmin,
     f     cosphi,f1,domin_skal1
      integer ierr
      sigmin = 9999999.0
      fmin = 9999999.0
      xnorm=domin_vecnor(n,x)
      dnorm=domin_vecnor(n,d)
      dskal=1.d0
      cosphi=domin_skal1(1,n,d,d0)
      t=1.d0
      do i=1,n
         if ( dabs(d(i)) .gt. (dabs(x(i))+1.d0)*(beta-1.d0) )
     f        t=dmin1(t,(beta-1.d0)*(dabs(x(i))+1.d0)/dabs(d(i)))
      enddo
      if ( t .lt. 1.d0  ) then
c*****reskalierung von d        
         dnorm=t*dnorm
         dirder=t*dirder
         cosphi=cosphi*t
         dskal=t
         do i= 1,n
            d(i)=d(i)*t
         enddo
      endif
      cosphi=dabs(cosphi)/(dnorm*d0norm)
c**** arithmetische fehler werden in f gemeldet und hier als versagen 
c     der schrittweitenbestimmung interpretiert, d.h. es erfolgt
c     schrittweitenverkleinerung
      sig=1.d0
 100  continue
c*****naechstes x berechnen und abstieg testen
      do i =1,n
         x1(i)=x(i)-sig*d(i)
      enddo
      ierr=0
      call domin_f(x1,f1,ierr)
c***  abfrage, ob arithmetischer fehler auftrat
      if ( ierr.ne.0 ) then
         if ( sig .gt. 1.d0 ) then
c***  fehler trat bei vergroesserung auf, vergroesserung abbrechen
            ierr=0
            sig=sigmin
            f1=fmin
            do i=1,n
               x1(i)=x(i)-sig*d(i)
            enddo
            goto 200
         endif
         sig=sig*.5d0
         if ( sig .ge. sigsm ) then
            goto 100
         else
            stepterm = -1.d0
            return
         endif
      endif
      ftest=f1
      call domin_info(6)
c*****neuer funktionswert
      if ( sig .gt. 1.d0 ) then
         if ( f1 .gt. fmin ) then
c*****f nimmt  nicht weiter ab
            sig=sigmin
            f1=fmin
            do i=1,n
               x1(i)=x(i)-sig*d(i)
            enddo
            goto 200
         else
            if ( sig .lt. sigla ) then
               fmin = f1
               sigmin = sig
               sig=dmin1(sigla,sig+sig)
               goto 100
            else
               goto 200
            endif
         endif
      endif
      if ( fx - f1 .ge. sig*delta*dirder ) then
c*****goldstein-armijo abstiegstest erfuellt
         if ( sig .eq. 1.d0 .and. cosphi .ge. theta .and. sig0 
     1        .ge. 1.d0  ) then
c*****versuche vergroesserung der schrittweite
            fmin = f1
            sigmin = sig
            if ( fx-f1 .lt. sig*dirder ) then
               sig=dmin1(sigla,
     1              dmax1(sig+sig,dirder*sig**2/
     1              (2.d0*(f1-fx+sig*dirder))))
            else
               sig=dmin1(sigla,sig+sig)
            endif
            goto 100
         endif
         goto 200
      else
         sig=dmax1(alpha*sig,
     1        dmin1(dirder*sig**2/(2.d0*(f1-fx+sig*dirder)),
     2        .5d0*sig))
         if ( sig .lt. sigsm ) then
            stepterm = -1.d0
            return
         else
            goto 100
         endif
      endif
 200  continue
c******akzeptiere neues x 
      fx0 = fx
      fx = f1
      stepterm = 1.d0 
      x0norm = xnorm
      d0norm = dnorm
      sig 0 = sig
      do i = 1,n
         gradx0(i)=gradx(i)
         x0(i)=x(i)
         d0(i)=d(i)
         x(i)=x1(i)
         difx(i)=x(i)-x0(i)
      enddo
      if ( analyt ) then
         call domin_gradf(x,gradx ,ierr)
      else
         call domin_numgradf(x,gradx,ierr)
      endif
      xnorm=domin_vecnor(n,x)
      if ( ierr.ne.0 ) then
c**** irregulaer: f auswertbar, gradf nicht
         stepterm=-2.d0
         return
      endif
      return
      end
c*******************************************************************
c     scalar product of two vectors or parts of vectors
c*******************************************************************
      double precision function domin_skal1(i,j,a,b)
      implicit none
      integer i,j,k
      double precision a(*),b(*)
      double precision s
      if ( i .gt. j ) then
         domin_skal1=0.d0
         return
      else
         s=0.d0
         do   k=i,j
            s=s+a(k)*b(k)
         enddo
         domin_skal1=s
         return
      endif
      end
c*******************************************************************
c     bfgs-method using modified update and cholesky-decomposition
c     inaccurate line search possible
c*******************************************************************
      subroutine domin_bfgs(cold1)
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
      logical cold1
      integer icfold,icgfold,icsf,i,j,ifail
      double precision v(NX),s1,s2,s3,s4,s5,difxin,wi,domin_vecnor,
     f     domin_skal1,den1
      double precision dif1,dif2,h
      logical restart
      integer ierr
      external domin_vecnor,domin_skal1
c*******initialization *********************
      restart = .false.
      icsf=0
      if ( cold1 ) then
         do   i=1,n
            x(i)=xst(i)
            x0(i)=x(i)
            difx(i)=1.d0+dabs(x(i))
c************to avoid termination in testx in step one
         enddo
         xnorm=domin_vecnor(n,x)
         call domin_initia
c********initia initializes the quasi-newton-matrix a
      endif
c********allows for modification of problem parameters
      sig0=0.d0
      d0norm=1.d0
c********no stepsize larger than 1 in the first step
      niter=0
      thet2=1.d0
      bfgsterm=0.d0
      do   i=1,n
         d0(i)=0.d0
      enddo
      ierr=0
      call domin_f(x,fx,ierr)
      if ( ierr.ne.0 ) then
         bfgsterm=-3.d0
         return
      endif
      ierr=0
      if ( analyt ) then
         call domin_gradf(x,gradx,ierr)
      else
         call domin_numgradf(x,gradx,ierr)
      endif
      if ( ierr.ne.0 ) then
         bfgsterm=-3.d0
         return
      endif
c******************************************************
c     main iteration loop starts here
c******************************************************
 100  continue
      niter=niter+1
      if ( niter .gt. itermax ) then
         niter=niter-1
         bfgsterm=-2.d0
         return
      endif
      do i=1,14
         accinf(niter,i)=0.d0
      enddo
      accinf(niter,1)=niter
      accinf(niter,2)=fx
      accinf(niter,3)=domin_vecnor(n,difx)
      accinf(niter,4)=domin_vecnor(n,gradx)
      accinf(niter,13)=-1.d0
      call domin_info(1)
      s1=diaga(1)
      s2=a(1,1)
      do i=2,n
         s1=s1+diaga(i)
         s2=dmin1(s2,a(i,i))
      enddo
c*****s1 = spur (a), s2=minimales element auf der diagonale des cholesky
c*****faktors von a
      if ( niter .gt. 1 ) then
         if ( (tracea*epsmach .gt. 1.d0 .and. traceb*epsmach .gt. 1.d0)
     f        .or. s1*epsmach .gt. s2**2 .or. 
     f        accinf(niter-1,8) .gt. 1.d0 /eta  ) then
            call domin_refresh(niter)
            if ( bfgsterm .eq. -3.d0 ) return
            restart=.true.
         else
            restart=.false.
         endif
      endif
      call domin_testx
      icfold=icf
      icgfold=icgf
c********check accuracy
      if ( bfgsterm .eq. 1.d0 )  return
 200  continue
      call domin_cholsol(a,n,gradx,d)
c********d = quasi-newton direction (uphill)
      dnorm=domin_vecnor(n,d)
      if ( dnorm .eq. 0.d0 ) then
         bfgsterm=0.d0
         return
      endif
      dirder=domin_skal1(1,n,d,gradx)
      accinf(niter,7)=dirder
      if ( delta*dirder .le. epsmach*(dabs(fx)+1.d0)*100.d0 ) then
         bfgsterm=2.d0
         return
      endif
      call domin_info(2)
      call domin_step
      accinf(niter,9)=sig
      accinf(niter,10)=dskal
      accinf(niter,11)=dnorm
      i=icf-icfold
      if ( .not. analyt ) i=i-6*n*(icgf-icgfold)
      accinf(niter,12)=dble(float(i))
      if ( test4 ) then
         write(*,fmt='(i5,'' f='',d15.7,'' n(g)='',d15.7,'' dird='',
     1        d15.7,'' sig='',d15.7)')
     2        niter,fx,accinf(niter,4),accinf(niter,7),accinf(niter,9)
      endif
      call domin_info(3)
      if ( stepterm .ne. 1.d0 ) then
         accinf(niter,13)=1.d0
         if ( restart ) then
            bfgsterm=-1.d0
            return
         else
            call domin_refresh(niter)
            if ( bfgsterm .eq. -3.d0 ) return
            restart=.true.
            goto 200
         endif
      endif
      if ( fx0-fx .gt. epsf*(dabs(fx)+epsf)) then
         icsf=0
      else
         icsf=icsf+1
      endif
c     make a change from 4 to 3 there; it's enough (hrue)
      if ( icsf .eq. 3 ) then
         bfgsterm=3.d0
         return
      endif
      if ( domin_vecnor(n,difx) .eq. 0 ) then
         bfgsterm=4.d0
         return
      endif
c*********************************************
c     compute update vectors z and y
c*********************************************
      thet2=1.d0
      s1=0.d0
      s2=0.d0
      s3=0.d0
      den1=domin_skal1(1,n,d,gradx0)/dskal
      h=1.d0/dsqrt(den1)
      do   i=1,n
         y(i)=gradx0(i)*h
         dif1=difx(i)
         dif2=gradx(i)-gradx0(i)
         s1=dif1**2 + s1
         s2=dif2**2 + s2
         s3=dif1*dif2 + s3
         z(i)=dif2
      enddo
      if ( s3 .lt. eta*s2 ) then
c**********************************************
c     bfgs-updating impossible from this information
c     f not uniformly convex of course
c**********************************************
         s4=((sig*dskal)**2)*den1
c********s(k)(transpose)*a(k)*s(k) in the usual notation
         if ( s3 .lt. .2d0*s4 ) thet2=.8d0*s4/(s4-s3)
         if ( thet2 .ne. 1.d0 ) imod=imod+1
         s3=0.2d0*s4
         h=(1.d0-thet2)*sig*dskal
c********powells modification
c********h*gradx0 = a(k)*s(k)*(1-thet2)
         do   i=1,n
            z(i)=thet2*z(i)+h*gradx0(i)
         enddo
      endif
      accinf(niter,6)=thet2
      s3=1.d0/dsqrt(s3)
      do   i=1,n
         z(i)=z(i)*s3
      enddo
      call domin_info(4)
      accinf(niter,14)=1.d0
      call domin_cholsol(a,n,z,v)
c**** v = b(k)*y(k)/sqrt(y's) in the usual notation
      s4=0.d0
      s5=0.d0
      do i=1,n
         difxin=difx(i)*s3
         wi=difxin+v(i)
c**** w=(sig(i)+(b*y)(i))/sqrt(y's)
         s5=s5+z(i)*wi
         s4=s4+v(i)*difxin
      enddo
c********update a now and compute its cholesky-decomposition
      if ( n .le. NSMALL ) then
c********direct updating, save old information
         do    i=1,n
            diaga0(i)=diaga(i)
            a0(i,i)=a(i,i)
            diaga(i)=diaga(i)-y(i)**2+z(i)**2
            a(i,i)=diaga(i)
c********the strict lower part of a(a0) holds the quasi-newton
c     matrix, its diagonal being diaga (diaga0)
c********the upper triangular part holds the cholesky-factor  
            do    j=1,i-1
               a0(i,j)=a(i,j)
               a0(j,i)=a(j,i)
               a(i,j)=a(i,j)-y(i)*y(j)+z(i)*z(j)
            enddo
         enddo
c------hrue: store the hessian explicite
         do i=1,n
            do j=1,i
               hessian(i,j) = a(i,j)
               hessian(j,i) = a(i,j)
            enddo
         enddo
c         print *,' domin.F: this is the current estimate of the Hessian'
c         do i=1,n
c            write(*,*) (hessian(i,j),j=1,n)
c         enddo
c      print *,' '
c------
         call domin_chol(a,n,a,ifail)
         if ( ifail .ne. 0 ) then
c********restore old matrices 
            accinf(niter,14)=-1.d0
            inoup=inoup+1
            do    i=1,n
               diaga(i)=diaga0(i)
               do    j=1,n
                  a(i,j)=a0(i,j)
               enddo
            enddo
         else
            iup=iup+1 
         endif
      else
c*********large n: quasi-newton matrix given as its 
c     cholesky-factor only. store/modify/restore
c*********this information only
         do    i=1,n
            diaga(i)=a(i,i)
            do    j=1,i-1
               a(i,j)=a(j,i)
            enddo
         enddo
         call domin_aufdat(a,z,y,n,ifail)
         if ( ifail .ne. 0 ) then
            inoup=inoup+1
            accinf(niter,14)=-1.d0
            do    i=1,n
               a(i,i)=diaga(i)
               do    j=1,i-1
                  a(j,i)=a(i,j)
               enddo
            enddo
         else
            iup=iup+1
         endif
      endif
c****************************************************
c     next step
c****************************************************
      if ( accinf (niter,14) .eq. 1.d0 ) then
c**** recursive computation of trace(a) and trace(b)
c**** it is stopped, if one trace becomes negative, which indicates
c**** unreliability
         if ( tracea .gt. 0.d0 )
     f        tracea=tracea+domin_vecnor(n,z)**2-domin_vecnor(n,y)**2
         if ( traceb .gt. 0.d0 )
     f        traceb=traceb-2.d0*s4+s1*s3**2*s5 
         call domin_info(7)
      endif
      call domin_info(5)
      goto 100
      end
c*******************************************************************
c     subprograms for computation of an updated cholesky-
c     decomposition, method of stewart
c*******************************************************************
      subroutine domin_leftelim(a,b,y,yl,n)
      implicit none
      integer domin_fileid
      external domin_fileid
      integer i,j,n
      double precision a(NX,NX),b(NX),y(NX),yl,h
c     leftelim assumes that the cholesky-factor of a
c     a=r(transpose)*r is stored in the upper half of a.
c     b is a right hand side. leftelim solves 
c     y(transpose)*r = b(transpose)
c     yl=norm(y)**2
      yl=0.d0
      do   i=1,n
         h=b(i)
         do    j=1,i-1
            h = h - a(j,i)*y(j)
         enddo
         h=h/a(i,i)
         y(i)=h
         yl = h**2 + yl
      enddo
      return
      end
c**********************************************************
      double precision function domin_dsq1(a,b)
      implicit none
      double precision a,b,a1,b1
      a1=dabs(a)
      b1=dabs(b)
      if ( a1 .gt. b1 ) then 
         domin_dsq1=a1*dsqrt(1.d0+(b1/a1)**2)
      else
         if ( b1 .gt. a1 ) then
            domin_dsq1=b1*dsqrt(1.d0+(a1/b1)**2)
         else
            domin_dsq1=a1*dsqrt(2.d0)
         endif
      endif
      return
      end
c*******************************************************************
      subroutine domin_aufdat(r,z,y,n,ifail)
      implicit none
      integer domin_fileid
      external domin_fileid
      integer ifail,i,j,n,i1
      double precision r(NX,NX),z(NX),y(NX)
      double precision sdiag(NX),rn1(NX),w(NX)
      double precision yl,zl,wl,wn1,ai,bi,h,domin_dsq1
c     aufdat computes  the upper triangular cholesky-factor
c     r1 of
c     r(transpose)*r+z*z(transpose)-y*y(transpose)
c     and restores it in r. the strict lower triangle of r re-
c     mains unchanged.
c     ifail=1 if the decomposition does'nt exist, ifail=2 on dimension
c     error, ifail=0 on success.
      if ( n .gt. NX ) then
         ifail=2
         return
      endif
      ifail=0
c     save subdiagonal
      do    i=1,n-1
         sdiag(i)=r(i+1,i)
         r(i+1,i)=0.d0
      enddo
c     step one: include z*z(transpose)
      zl=0.d0
      do   i=1,n
         zl = zl + z(i)**2
      enddo
      if ( zl .ne. 0.d0 ) then
c     solve w(transpose)*r=z(transpose)
         call domin_leftelim(r,z,w,wl,n)
         wl=dsqrt(wl+1.d0)
c     u(2)*u(3)*...*u(n)*w = ( norm(w),0,..,0)(transpose)
c     u(i) rotations
         do   i=n,2,-1
            if ( w(i) .ne. 0.d0 ) then
               i1=i-1
               ai=w(i1)
               bi=w(i)
               w(i1)=domin_dsq1(ai,bi)
               ai=ai/w(i1)
               bi=-bi/w(i1)
               r(i,i1)=bi*r(i1,i1)
               r(i1,i1)=ai*r(i1,i1)
               do   j=i,n
                  h = ai*r(i1,j) - bi*r(i,j)
                  r(i,j) = bi*r(i1,j) + ai*r(i,j)
                  r(i1,j) = h
               enddo
            endif
         enddo
c     r=d*r, d=diag(wl,1,...,1), r now hessenberg
         do    i=1,n
            r(1,i)=r(1,i)*wl
         enddo
c     r=u(n-1)*...*u(1)*r now upper triangular, 
c     u(i)  givens-rotations
         do    i=1,n-1
            i1=i+1
            ai=r(i,i)
            bi=-r(i1,i)
            h=domin_dsq1(ai,bi)
            ai=ai/h
            bi=bi/h
            r(i,i)=h
            r(i1,i)=0.d0
            do   j=i+1,n
               h = ai*r(i,j) - bi*r(i1,j)
               r(i1,j) = bi*r(i,j) + ai*r(i1,j)
               r(i,j) = h
            enddo
         enddo
      endif
c     step two :  include -y*y(transpose)
      yl=0.d0
      do   i=1,n
         yl = yl + y(i)**2
      enddo
      if ( yl .ne. 0.d0 ) then
         call domin_leftelim(r,y,w,wl,n)
         if ( wl .ge. 1.d0 ) then
            ifail=1
         else
            wl=dsqrt(1.d0-wl)
            wn1=wl
c******************************************************
c     ( r(new) ,0 )                (    r  , w )
c     (-----------) = u(1)*...u(n)*(-----------)
c     (y(transp),1)                ((0,..,0),wl)             
c******************************************************
            do    i=n,1,-1
               ai=wn1
               bi=w(i)
               wn1=domin_dsq1(ai,bi)
               ai=ai/wn1
               bi=bi/wn1
               rn1(i)=bi*r(i,i)
               r(i,i)=ai*r(i,i)
               do    j=i+1,n
                  h = ai*r(i,j) - bi*rn1(j)
                  rn1(j) = bi*r(i,j) + ai*rn1(j)
                  r(i,j) = h
               enddo
            enddo
         endif
      endif
c     restore subdiagonal
      do   i=1,n-1
         r(i+1,i)=sdiag(i)
      enddo
      return
      end
c*******************************************************************
c     cholesky decomposition of a symmetric positive matrix a
c     the strict lower triangle remains unaffected
c     the upper triangle including the diagonal holds the cholesky-
c     factor. initially the strict upper triangle may be undefined
c*******************************************************************
      subroutine domin_chol(a,n,r,ifail)
      implicit none
      integer domin_fileid
      external domin_fileid
c     copyright = spellucci 
      integer i,j,k,ifail,n
      double precision a(NX,NX),r(NX,NX),h,s
c     computes the cholesky-factor r of a=r(transp)*r
c     and stores it in the upper triangle of r
c     a and r may be identical. the strict lower triangle
c     of a remains unchanged anyway.
c     ifail .ne. 0 if the decomposition does'nt exist,
c     =0 otherwise
      ifail=0
      do    i=1,n
         h=a(i,i)
         do   j=1,i-1
            h=h-r(j,i)**2
         enddo
         if ( h .le. 0.d0 ) then
            ifail=1
            return
         endif
         h=dsqrt(h)
         r(i,i)=h
         h=1.d0/h
         do    k=i+1,n
            s=0.d0
            do    j=1,i-1
               s=s+r(j,i)*r(j,k)
            enddo
            s=(a(k,i)-s)*h
            r(i,k)=s
         enddo
      enddo
      return
      end
      subroutine domin_cholsol(a,n,b,x)
      implicit none
      integer domin_fileid
      external domin_fileid
c     copyright = spellucci 
      integer i,j,n
      double precision a(NX,NX),b(NX),c(NX),x(NX)
      double precision s
c     solves the linear system a*x=b, a symmetric positive definite
c     it is assumed that the cholesky-factor r of
c     a = r(transpose)*r
c     is stored in the upper triangle of a by chol.
      do   i=1,n
         s=0.d0
         do    j=1,i-1
            s=s+a(j,i)*c(j)
         enddo
         c(i)=(b(i)-s)/a(i,i)
      enddo
      do    i=n,1,-1
         s=0.d0
         do    j=i+1,n
            s=s+a(i,j)*x(j)
         enddo
         x(i)=(c(i)-s)/a(i,i)
      enddo
      return
      end
c*******************************************************************
c     initia is used for initializing the quasi newton-update
c*******************************************************************
      subroutine domin_initia
c     copyright = spellucci 
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
      integer i,j
      double precision xd(NX),xst1(NX),gradxs(NX),
     f     gradx1(NX),gn,domin_vecnor
c*******initialize with scaled identity
      integer ierr
      external domin_vecnor
      do i=1,n
c***  xst1 should be in the range of definiton of gradf
         xst1(i)=x(i)*(1.d0+epsx)+dsign(epsx**2,x(i))
         xd(i)=x(i)-xst1(i)
      enddo 
      ierr=0
      if ( analyt ) then
         call domin_gradf(x,gradxs,ierr)
      else
         call domin_numgradf(x,gradxs,ierr)
      endif
      if ( ierr.ne.0 ) then
         bfgsterm=-3.d0
         return
      endif
      gn=domin_vecnor(n,gradxs)
      if ( accinf(1,1) .eq. 0.d0 ) epsg=epsg*(gn+epsg)
      if ( analyt ) then
         call domin_gradf(xst1,gradx1,ierr)
      else
         call domin_numgradf(xst1,gradx1,ierr)
      endif
      if ( ierr.ne.0 ) then
         bfgsterm=-3.d0
         return
      endif
      do i=1,n
         gradxs(i)=gradxs(i)-gradx1(i)
      enddo
      gn=domin_vecnor(n,gradxs)
      gn=gn/domin_vecnor(n,xd)/dble(float(n))
      if ( gn .le. dsqrt(epsmach) ) gn = 1.d0
c**** gn should be an estimate of the order of magnitude of the hessian
c**** in case of a very small gn the computation of gn itself is expected
c**** to be unreliable
      if ( n .le. NSMALL ) then
         do    i=1,n
            do    j=1,n
               a(i,j)=0.d0
               a0(i,j)=0.d0
            enddo
            a(i,i)=dsqrt(gn)
            diaga(i)=gn
            diaga0(i)=gn
         enddo
      else
         do    i=1,n
            do    j=1,n
               a(j,i)=0.d0
            enddo
            a(i,i)=dsqrt(gn)
            diaga(i)=gn
         enddo
      endif
      tracea=dble(float(n))*gn
      traceb=dble(float(n))/gn
      return
      end
c*******************************************************************
c     refresh is used on restart in bfgs
c     restart is done only on occasion of a too small stepsize
c     sig .lt. sigsm
c*******************************************************************
      subroutine domin_refresh(niter1)
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
      integer niter1
c     copyright = spellucci 
      ires=ires+1
      accinf(niter1,13)=1.d0
      call domin_initia
c********restart with scaled unit matrix
      return
      end
c*******************************************************************
c     testx checks for sufficient accuracy
c*******************************************************************
      subroutine domin_testx
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
      logical xsatisfy,gsatisfy
      double precision an,sx,sg,damin,domin_vecnor
      integer i
      an=0.d0
      sx=domin_vecnor(n,difx)
      sg=domin_vecnor(n,gradx)
      damin=a(1,1)
      do   i=1,n
         an=an+diaga(i)
         damin=dmin1(damin,a(i,i))
      enddo
      accinf(niter,5)=an
      if ( tracea .gt. 0.d0 .and. traceb .gt. 0.d0 ) then
         accinf(niter,8)=tracea*traceb/dble(float(n*n))
      else
         accinf(niter,8)=an/damin**2
      endif
      xsatisfy = sx .le. epsx*(xnorm+epsx)
      gsatisfy = sg .le. dmax1(epsg,epsmach*accinf(niter,8)) 
      if ( xsatisfy .and. gsatisfy ) bfgsterm=1.d0
      return
      end
c*******************************************************************
c     info provides for intermediate output for 
c     testing purposes
c     info may have an empty body of course
c*******************************************************************
      subroutine domin_info(itestloc)
c     copyright = spellucci 
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
      integer itestloc,i,j
      if ( .not. test2 ) return
      goto ( 1000,2000,3000,4000,5000,6000,7000 ) ,itestloc
 1000 continue
      if (logit) then
         write(domin_fileid(),1100) niter,fx,accinf(niter,4)
 1100    format(1h0,'******************************************'
     1        /1x,i3,'-ter iteration step',
     2        /1x,'f = ',d21.14,' norm(gradx)=',d21.14)
         write(domin_fileid(),*) ' x = '
         write(domin_fileid(),fmt='(1x,5(d21.14,2x))') (x(j),j=1,n)
      endif
      if ( test3 .and. n .le. NSMALL) then 
         call domin_matdru(a,NX,NX,n,n,'a/r-mat ',10,.false.)
         if (logit) then
            write(domin_fileid(),*)
     f           'diagonal of updated quasi newton matrix'
            write(domin_fileid(),1200) (diaga(j),j=1,n)
 1200       format(1x,5(d21.14,2x))
         endif
      endif
      return
 2000 continue
      if (logit) then
         write(domin_fileid(),2100)    accinf(niter,5),accinf(niter,8)
 2100    format(/'norm(a)=',d10.4,'  cond(a)=',d10.4)
         write(domin_fileid(),*) ' d = '
         write(domin_fileid(),fmt='(1x,5(d21.14,2x))') (d(j),j=1,n)
      endif
      return
 3000 continue
      if (logit) then
         write(domin_fileid(),3100) (accinf(niter,j),j=9,12)
 3100    format(/,'sigma=',d10.4,'  dskal=',d10.4,' norm(d)=',d10.4,
     1        ' number of function calls=',d10.4)
      endif
      return
 4000 continue
      if (logit) then
         write(domin_fileid(),4100) thet2
         write(domin_fileid(),*) ' update vektor z '
         write(domin_fileid(),4200) (z(i),i=1,n)
         write(domin_fileid(),*) ' update vektor y '
         write(domin_fileid(),4200) (y(i),i=1,n)
      endif
      return
 4100 format(/,'parameter theta2 =',d10.4)
 4200 format(1x,5(d21.14,2x))
 5000 continue
      if (logit) then
         write(domin_fileid(),5100) accinf(niter,13),accinf(niter,14)
 5100    format(/'restart =',d8.1,' update = ',d8.1)
      endif
      return
 6000 continue
      if (logit) then
         write(domin_fileid(),6100) sig,ftest
 6100    format(1x,'step***  sig = ',d10.4,'  f(x-sig*d)=',d21.14)
      endif
      return
 7000 continue
      if (logit) then
         write(domin_fileid(),
     f        fmt='(1x,''trace a='',d15.7,'' trace b='',d15.7)')
     f        tracea,traceb
      endif
      end
c*****************************************************************
c     subprogram for structured output of a matrix 
c*****************************************************************
      subroutine domin_matdru(a,me,ne,ma,na,heading,lognum,fix)
      implicit none
      integer domin_fileid
      external domin_fileid
      integer lognum,me,ma,ne,na,ju,jo,i,j
      double precision a(me,ne)
      logical fix
      character*(*) heading
      integer anz
      if ( ma .gt. 999 .or. na .gt. 999 ) then
         write(lognum
     f        ,fmt='(''call to matdru suppressed, dim too large'')')
         return
      endif
      write(lognum,10)  heading
 10   format(/1x,a40)
      anz=8
      jo=0
 20   continue
      ju=jo+1
      jo=min(ju+anz-1,na)
      write(lognum,30)  (j,j=ju,jo)
 30   format(/1x,'row/column',8(6x,i3,6x))
      do i=1,ma
         if ( fix ) then
            write(lognum,40) i,(a(i,j),j=ju,jo)
 40         format(1x,3x,i4,3x,8(g14.7,1x))
         else
            write(lognum,50) i,(a(i,j),j=ju,jo)
 50         format(1x,3x,i4,3x,8(d14.7,1x))
         endif
      enddo
      if ( jo .lt. na ) goto 20
      return
      end
c*****************************************************************
      double precision function domin_vecnor(n,x)
      implicit none
c***  euclidean norm of x , avoid overflow
      integer n,i
      double precision x(*)
      double precision xm,h
      if ( n .le. 0 ) then
         domin_vecnor=0.d0
         return
      endif
      xm=dabs(x(1))
      do i=2,n
         xm=dmax1(xm,dabs(x(i)))
      enddo
      if ( xm .eq. 0.d0 ) then
         domin_vecnor=0.d0
         return
      else
         h=0.d0
         do i=1,n
            h=h+(x(i)/xm)**2
         enddo
         domin_vecnor=xm*dsqrt(h)
         return
      endif
      end    
c***************************************************************
c****************************************************************
c     computes gradient by high precision numerical differentiation
c*******************************************************************
c     copyright = spellucci
      subroutine domin_numgradf(xl,gradxl,ierr)
#include "domincommon.inc"
      integer domin_fileid
      external domin_fileid
      integer ierr
      double precision xl(*),gradxl(*),fxl(6)
      double precision sd1,sd2,sd3,xincr,zz,d1,d2,d3
      integer i
c****************************************************
c     high precision numerical differentiation
c     by sixth order extrapolation
c****************************************************
      ierr=0
      do   i=1,n
         zz=xl(i)
         xincr=deldif*(dabs(xl(i))+1.d0)
         xl(i)=zz-xincr
         call domin_f(xl,fxl(1),ierr)
         if ( ierr.ne.0 ) return
         xl(i)=zz+xincr
         call domin_f(xl,fxl(2),ierr)
         if (ierr.ne.0) return
         xincr=xincr+xincr
         d1=xincr
         xl(i)=zz-xincr
         call domin_f(xl,fxl(3),ierr)
         if (ierr.ne.0) return
         xl(i)=zz+xincr
         call domin_f(xl,fxl(4),ierr)
         if (ierr.ne.0) return
         xincr=xincr+xincr
         d2=xincr
         xl(i)=zz-xincr
         call domin_f(xl,fxl(5),ierr)
         if (ierr.ne.0) return
         xl(i)=zz+xincr
         call domin_f(xl,fxl(6),ierr)
         if (ierr.ne.0) return
         xl(i)=zz
         d3=xincr+xincr
         sd1=(fxl(2)-fxl(1))/d1
         sd2=(fxl(4)-fxl(3))/d2
         sd3=(fxl(6)-fxl(5))/d3
         sd3=sd2-sd3
         sd2=sd1-sd2
         sd3=sd2-sd3
         gradxl(i)=sd1+.4d0*sd2+sd3/45.d0
      enddo
      icgf=icgf+1
      return
      end
