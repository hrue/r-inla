c/*!
*     \file GMRFLib-fortran.F
*     \brief Fortran routines used in GMRFLib
*
*     This file contains some fortran routines used internally
*     in GMRFLib.
*
*     See the source-code for details.
*/      
*
*      
      subroutine cmsd (cmean, csd, indx, n, k, a, lda, x)
      integer            k, lda, n,indx
      double precision   a( lda, * ), x( * ), cmean,csd
*
*  purpose
*  =======
*
*  return the conditional mean and stddev for x(indx)|x(indx+1:n)
*  for current values of x, in  'cmean' and 'csd'
*
*  note that 'indx' is in fortran-index's.      
*
*  the remaining parameters 'n,k,a,lda,x' is as is dtbsv.f when
*  called as
*  dtbsv_("L", "T", "N", &n, &nband, bchol, &ldim, rhs, &stride,1,1,1)
*     
      integer i,l
      double precision temp
      intrinsic min

      if( n.eq.0 ) then
         cmean = 0.0
         csd   = 0.0
         return
      endif

      temp = 0.0
      l    = 1 - indx
      do i = min(n, indx+k), indx+1, -1
         temp = temp - a(l+i, indx)*x(i)
      enddo
      cmean = temp/a(1,indx)
      csd = 1/a(1,indx)
      return
      end


      subroutine dtbsvSpecial(uplo,trans,diag,n,k,a,lda,x,incx,
     $     findx,toindx)
      integer            incx, k, lda, n,findx,toindx
      character*1        diag, trans, uplo
      double precision   a( lda, * ), x( * )
*  
*       
*  dtbsv  solves [a is lower band matrix]
*     a'*x = b,
*  but only from findx to toindx, both included.
*  note that it solves backward, so findx >= toindx, and that
*  findx and toindx are fortran-index's going from 1..n.      
*      
      parameter        ( zero = 0.0d+0 )
      double precision   temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            nounit
      logical            lsame
      external           lsame
      external           xerbla
      intrinsic          max, min
      kx = 0
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( k.lt.0 )then
         info = 5
      else if( lda.lt.( k + 1 ) )then
         info = 7
      else if( incx.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'dtbsv ', info )
         return
      end if
*
*     quick return if possible.
*
      if( n.eq.0 )
     $   return
*
      nounit = lsame( diag, 'n' )
*
*     set up the start point in x if the increment is not unity. this
*     will be  ( n - 1 )*incx  too small for descending loops.
*
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
*
*     start the operations. in this version the elements of a are
*     accessed by sequentially with one pass through a.
*
      if( lsame( trans, 'n' ) )then
*
*        form  x := inv( a )*x.
*
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     l = kplus1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( kplus1, j )
                     temp = x( j )
                     do 10, i = j - 1, max( 1, j - k ), -1
                        x( i ) = x( i ) - temp*a( l + i, j )
   10                continue
                  end if
   20          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 40, j = n, 1, -1
                  kx = kx - incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = kplus1 - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( kplus1, j )
                     temp = x( jx )
                     do 30, i = j - 1, max( 1, j - k ), -1
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      - incx
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     l = 1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( 1, j )
                     temp = x( j )
                     do 50, i = j + 1, min( n, j + k )
                        x( i ) = x( i ) - temp*a( l + i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  kx = kx + incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = 1  - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( 1, j )
                     temp = x( jx )
                     do 70, i = j + 1, min( n, j + k )
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      + incx
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
*
*        form  x := inv( a')*x.
*
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  l    = kplus1 - j
                  do 90, i = max( 1, j - k ), j - 1
                     temp = temp - a( l + i, j )*x( i )
   90             continue
                  if( nounit )
     $               temp = temp/a( kplus1, j )
                  x( j ) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  l    = kplus1  - j
                  do 110, i = max( 1, j - k ), j - 1
                     temp = temp - a( l + i, j )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit )
     $               temp = temp/a( kplus1, j )
                  x( jx ) = temp
                  jx      = jx   + incx
                  if( j.gt.k )
     $               kx = kx + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               if (findx.lt.1) return
************************************************************************
* her er endringene: loop gaar fra findx til toindx istedet for n til 1
c               print *,findx,toindx,' HER'
               do 140, j = findx, toindx, -1
************************************************************************
                  temp = x( j )
                  l    = 1      - j
                  do 130, i = min( n, j + k ), j + 1, -1
                     temp = temp - a( l + i, j )*x( i )
c                     print *, a(l+i,j),x(i)
  130             continue
                  if( nounit )
     $               temp = temp/a( 1, j )
                  x( j ) = temp
c                  print *,'set x(j)',j,x(j),a(1,j)
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  l    = 1       - j
                  do 150, i = min( n, j + k ), j + 1, -1
                     temp = temp - a( l + i, j )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit )
     $               temp = temp/a( 1, j )
                  x( jx ) = temp
                  jx      = jx   - incx
                  if( ( n - j ).ge.k )
     $               kx = kx - incx
  160          continue
            end if
         end if
      end if
*
      return
*
*     end of dtbsv .
*
      end

      double precision function erfi(x)
c     this function computes
c     erfi(x) = 1/sqrt(Pi)*int(exp(t^2),t=0..x);
      double precision x,y,u,v
      logical flag
      y=0.0
      call wofz(x,y,u,v,flag)
      erfi = exp(x**2)*v
      end

      subroutine wofz (xi, yi, u, v, flag)
c     algorithm 680, collected algorithms from acm.  this work published
c     in transactions on mathematical software, vol. 16, no. 1, pp. 47.
c
c     given a complex number z = (xi,yi), this subroutine computes the
c     value of the faddeeva-function w(z) = exp(-z**2)*erfc(-i*z), where
c     erfc is the complex complementary error-function and i means
c     sqrt(-1).  the accuracy of the algorithm for z in the 1st and 2nd
c     quadrant is 14 significant digits; in the 3rd and 4th it is 13
c     significant digits outside a circular region with radius 0.126
c     around a zero of the function.  all real variables in the program
c     are double precision.
c
c
c  the code contains a few compiler-dependent parameters :
c     rmaxreal = the maximum value of rmaxreal equals the root of
c                rmax = the largest number which can still be
c                implemented on the computer in double precision
c                floating-point arithmetic
c     rmaxexp  = ln(rmax) - ln(2)
c     rmaxgoni = the largest possible argument of a double precision
c                goniometric function (dcos, dsin, ...)
c  the reason why these parameters are needed as they are defined will
c  be explained in the code by means of comments
c
c
c  parameter list
c     xi     = real      part of z
c     yi     = imaginary part of z
c     u      = real      part of w(z)
c     v      = imaginary part of w(z)
c     flag   = an error flag indicating whether overflow will
c              occur or not; type logical;
c              the values of this variable have the following
c              meaning :
c              flag=.false. : no error condition
c              flag=.true.  : overflow will occur, the routine
c                             becomes inactive
c  xi, yi      are the input-parameters
c  u, v, flag  are the output-parameters
c
c  furthermore the parameter factor equals 2/sqrt(pi)
c
c  the routine is not underflow-protected but any variable can be
c  put to 0 upon underflow;
c
c  reference - gpm poppe, cmj wijers; more efficient computation of
c  the complex error-function, acm trans. math. software.
c
      implicit double precision (a-h, o-z)
      logical a, b, flag
      parameter (factor   = 1.12837916709551257388d0,
     *           rmaxreal = 0.5d+154,
     *           rmaxexp  = 708.503061461606d0,
     *           rmaxgoni = 3.53711887601422d+15)

c     get rid of warnings
      h2   = 0.0
      qlambda = 0.0
      v2 = 0.0
      u2 = 0.0

      flag = .false.
      xabs = dabs(xi)
      yabs = dabs(yi)
      x    = xabs/6.3
      y    = yabs/4.4
c
c     the following if-statement protects
c     qrho = (x**2 + y**2) against overflow
c
      if ((xabs.gt.rmaxreal).or.(yabs.gt.rmaxreal)) goto 100
      qrho = x**2 + y**2
      xabsq = xabs**2
      xquad = xabsq - yabs**2
      yquad = 2*xabs*yabs
      a     = qrho.lt.0.085264d0
      if (a) then
c
c  if (qrho.lt.0.085264d0) then the faddeeva-function is evaluated
c  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
c  n is the minimum number of terms needed to obtain the required
c  accuracy
c
        qrho  = (1-0.85*y)*dsqrt(qrho)
        n     = idnint(6 + 72*qrho)
        j     = 2*n+1
        xsum  = 1.0/j
        ysum  = 0.0d0
        do 10 i=n, 1, -1
          j    = j - 2
          xaux = (xsum*xquad - ysum*yquad)/i
          ysum = (xsum*yquad + ysum*xquad)/i
          xsum = xaux + 1.0/j
 10     continue
        u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0
        v1   =  factor*(xsum*xabs - ysum*yabs)
        daux =  dexp(-xquad)
        u2   =  daux*dcos(yquad)
        v2   = -daux*dsin(yquad)
        u    = u1*u2 - v1*v2
        v    = u1*v2 + v1*u2
      else
c
c  if (qrho.gt.1.o) then w(z) is evaluated using the laplace
c  continued fraction
c  nu is the minimum number of terms needed to obtain the required
c  accuracy
c
c  if ((qrho.gt.0.085264d0).and.(qrho.lt.1.0)) then w(z) is evaluated
c  by a truncated taylor expansion, where the laplace continued fraction
c  is used to calculate the derivatives of w(z)
c  kapn is the minimum number of terms in the taylor expansion needed
c  to obtain the required accuracy
c  nu is the minimum number of terms of the continued fraction needed
c  to calculate the derivatives with the required accuracy
c
        if (qrho.gt.1.0) then
          h    = 0.0d0
          kapn = 0
          qrho = dsqrt(qrho)
          nu   = idint(3 + (1442/(26*qrho+77)))
        else
          qrho = (1-y)*dsqrt(1-qrho)
          h    = 1.88*qrho
          h2   = 2*h
          kapn = idnint(7  + 34*qrho)
          nu   = idnint(16 + 26*qrho)
        endif
        b = (h.gt.0.0)
        if (b) qlambda = h2**kapn
        rx = 0.0
        ry = 0.0
        sx = 0.0
        sy = 0.0
        do 11 n=nu, 0, -1
          np1 = n + 1
          tx  = yabs + h + np1*rx
          ty  = xabs - np1*ry
          c   = 0.5/(tx**2 + ty**2)
          rx  = c*tx
          ry  = c*ty
          if ((b).and.(n.le.kapn)) then
            tx = qlambda + sx
            sx = rx*tx - ry*sy
            sy = ry*tx + rx*sy
            qlambda = qlambda/h2
          endif
 11     continue
        if (h.eq.0.0) then
          u = factor*rx
          v = factor*ry
        else
          u = factor*sx
          v = factor*sy
        end if
        if (yabs.eq.0.0) u = dexp(-xabs**2)
      end if
c
c  evaluation of w(z) in the other quadrants
c
      if (yi.lt.0.0) then
        if (a) then
          u2    = 2*u2
          v2    = 2*v2
        else
          xquad =  -xquad
c
c         the following if-statement protects 2*exp(-z**2)
c         against overflow
c
          if ((yquad.gt.rmaxgoni).or.
     *        (xquad.gt.rmaxexp)) goto 100
          w1 =  2*dexp(xquad)
          u2  =  w1*dcos(yquad)
          v2  = -w1*dsin(yquad)
        end if
        u = u2 - u
        v = v2 - v
        if (xi.gt.0.0) v = -v
      else
        if (xi.lt.0.0) v = -v
      end if
      return
  100 flag = .true.
      return
      end

      subroutine dchdc(a,lda,p,work,jpvt,job,info,eps)
      integer lda,p,jpvt(*),job,info
      double precision a(lda,*),work(*),eps
c
c     dchdc computes the cholesky decomposition of a positive definite
c     matrix.  a pivoting option allows the user to estimate the
c     condition of a positive definite matrix or determine the rank
c     of a positive semidefinite matrix.
c
c     on entry
c
c         a      double precision(lda,p).
c                a contains the matrix whose decomposition is to
c                be computed.  onlt the upper half of a need be stored.
c                the lower part of the array a is not referenced.
c
c         lda    integer.
c                lda is the leading dimension of the array a.
c
c         p      integer.
c                p is the order of the matrix.
c
c         work   double precision.
c                work is a work array.
c
c         jpvt   integer(p).
c                jpvt contains integers that control the selection
c                of the pivot elements, if pivoting has been requested.
c                each diagonal element a(k,k)
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      element.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free element.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final element.
c
c                before the decomposition is computed, initial elements
c                are moved by symmetric row and column interchanges to
c                the beginning of the array a and final
c                elements to the end.  both initial and final elements
c                are frozen in place during the computation and only
c                free elements are moved.  at the k-th stage of the
c                reduction, if a(k,k) is occupied by a free element
c                it is interchanged with the largest free element
c                a(l,l) with l .ge. k.  jpvt is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c         a      a contains in its upper half the cholesky factor
c                of the matrix a as it has been permuted by pivoting.
c
c         jpvt   jpvt(j) contains the index of the diagonal element
c                of a that was moved into the j-th position,
c                provided pivoting was requested.
c
c         info   contains the index of the last positive diagonal
c                element of the cholesky factor.
c
c     for positive definite matrices info = p is the normal return.
c     for pivoting with positive semidefinite matrices info will
c     in general be less than p.  however, info may be greater than
c     the rank of a, since rounding error can cause an otherwise zero
c     element to be positive. indefinite systems will always cause
c     info to be less than p.
c
c     linpack. this version dated 08/14/78 .
c     j.j. dongarra and g.w. stewart, argonne national laboratory and
c     university of maryland.
c
c
c     blas daxpy,dswap
c     fortran dsqrt
c
c     internal variables
c
      integer pu,pl,plp1,j,jp,jt,k,kb,km1,kp1,l,maxl
      double precision temp
      double precision maxdia
      logical swapk,negk
c
      pl = 1
      pu = 0
      info = p
      if (job .eq. 0) go to 160
c
c        pivoting has been requested. rearrange the
c        the elements according to jpvt.
c
         do 70 k = 1, p
            swapk = jpvt(k) .gt. 0
            negk = jpvt(k) .lt. 0
            jpvt(k) = k
            if (negk) jpvt(k) = -jpvt(k)
            if (.not.swapk) go to 60
               if (k .eq. pl) go to 50
                  call dswap(pl-1,a(1,k),1,a(1,pl),1)
                  temp = a(k,k)
                  a(k,k) = a(pl,pl)
                  a(pl,pl) = temp
                  plp1 = pl + 1
                  if (p .lt. plp1) go to 40
                  do 30 j = plp1, p
                     if (j .ge. k) go to 10
                        temp = a(pl,j)
                        a(pl,j) = a(j,k)
                        a(j,k) = temp
                     go to 20
   10                continue
                     if (j .eq. k) go to 20
                        temp = a(k,j)
                        a(k,j) = a(pl,j)
                        a(pl,j) = temp
   20                continue
   30             continue
   40             continue
                  jpvt(k) = jpvt(pl)
                  jpvt(pl) = k
   50          continue
               pl = pl + 1
   60       continue
   70    continue
         pu = p
         if (p .lt. pl) go to 150
         do 140 kb = pl, p
            k = p - kb + pl
            if (jpvt(k) .ge. 0) go to 130
               jpvt(k) = -jpvt(k)
               if (pu .eq. k) go to 120
                  call dswap(k-1,a(1,k),1,a(1,pu),1)
                  temp = a(k,k)
                  a(k,k) = a(pu,pu)
                  a(pu,pu) = temp
                  kp1 = k + 1
                  if (p .lt. kp1) go to 110
                  do 100 j = kp1, p
                     if (j .ge. pu) go to 80
                        temp = a(k,j)
                        a(k,j) = a(j,pu)
                        a(j,pu) = temp
                     go to 90
   80                continue
                     if (j .eq. pu) go to 90
                        temp = a(k,j)
                        a(k,j) = a(pu,j)
                        a(pu,j) = temp
   90                continue
  100             continue
  110             continue
                  jt = jpvt(k)
                  jpvt(k) = jpvt(pu)
                  jpvt(pu) = jt
  120          continue
               pu = pu - 1
  130       continue
  140    continue
  150    continue
  160 continue
      do 270 k = 1, p
c
c        reduction loop.
c
         maxdia = a(k,k)
         kp1 = k + 1
         maxl = k
c
c        determine the pivot element.
c
         if (k .lt. pl .or. k .ge. pu) go to 190
            do 180 l = kp1, pu
               if (a(l,l) .le. maxdia) go to 170
                  maxdia = a(l,l)
                  maxl = l
  170          continue
  180       continue
  190    continue
c
c        quit if the pivot element is not positive.
c
         if (maxdia .gt. eps) go to 200
c           print *, 'SINGULAR, threshold ',eps
            info = k - 1
c     ......exit
            go to 280
  200    continue
         if (k .eq. maxl) go to 210
c
c           start the pivoting and update jpvt.
c
            km1 = k - 1
            call dswap(km1,a(1,k),1,a(1,maxl),1)
            a(maxl,maxl) = a(k,k)
            a(k,k) = maxdia
            jp = jpvt(maxl)
            jpvt(maxl) = jpvt(k)
            jpvt(k) = jp
  210    continue
c
c        reduction step. pivoting is contained across the rows.
c
         work(k) = dsqrt(a(k,k))
         a(k,k) = work(k)
         if (p .lt. kp1) go to 260
         do 250 j = kp1, p
            if (k .eq. maxl) go to 240
               if (j .ge. maxl) go to 220
                  temp = a(k,j)
                  a(k,j) = a(j,maxl)
                  a(j,maxl) = temp
               go to 230
  220          continue
               if (j .eq. maxl) go to 230
                  temp = a(k,j)
                  a(k,j) = a(maxl,j)
                  a(maxl,j) = temp
  230          continue
  240       continue
            a(k,j) = a(k,j)/work(k)
            work(j) = a(k,j)
            temp = -a(k,j)
            call daxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
  250    continue
  260    continue
  270 continue
  280 continue
      return
      end

      integer function strlen(st)
      integer i
      character st*(*)
      i = len(st)
      do while (st(i:i) .eq. ' ')
         i = i - 1
      enddo
      strlen = i
      return
      end

      SUBROUTINE DGEMV_FAILSAFE(TRANS,M,N,ALPHA,A,LDA,
     +                          X,INCX,BETA,Y,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  END IF
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
