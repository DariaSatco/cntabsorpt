!*******************************************************************************
!*******************************************************************************
! Project      : libMath.f90
!===============================================================================
! Purpose      :
! Library for some mathematical utilites (functions and subroutines)
!-------------------------------------------------------------------------------
! Author       : ART Nugraha (nugraha@flex.phys.tohoku.ac.jp)
! Latest Vers. : 2012.10.18
!-------------------------------------------------------------------------------
! Reference(s) :
! Numerical Recipes in Fortran
!-------------------------------------------------------------------------------
! Contents     :
! - FUNCTION igcd(n1,n2)
! - FUNCTION vecLength(n,V)
! - FUNCTION pythagoras(a,b)
! - FUNCTION chebev(a,b,c,m,x)
! - FUNCTION erf(x)
! - FUNCTION fexp(x)
! - FUNCTION diracDelta(x1,y1,x2,y2,fwhm)
! - FUNCTION rtbis(func,x1,x2,xacc) 
! - SUBROUTINE outProd(n1,V1,n2,V2,ldu,U)
! - SUBROUTINE linArray(nx,xmin,xmax,xarray)
! - SUBROUTINE dos1Dgauss(nb,nk,rka,ldenk,Enk,E,fwhm,nn,DSn)
! - SUBROUTINE solveHam(n,ldh,ham,ldo,ovlp,matz,il,iu,nout,w,ldz,z)
! - SUBROUTINE zbrac(func,x1,x2,succes)
! - SUBROUTINE sort2(N,Ra,Rb)
! - SUBROUTINE indexx(n,arr,indx)
! - SUBROUTINE ode2(vector,method,nvar,t1,t2,eps,h1,derivs,maxstp,nstp,hnex)
! - SUBROUTINE rkqs(y,dvdt,n,x,htry,eps,yscal,hdid,hnext,derivs)
! - SUBROUTINE rkck(y,dvdt,n,x,h,yout,yerr,derivs)

!*******************************************************************************
!*******************************************************************************
INTEGER FUNCTION igcd(n1,n2)
!===============================================================================
! Calculate greatest common divisor of two integers n1 and n2
!===============================================================================
  IMPLICIT NONE

! input variables:
  INTEGER, INTENT(in)    :: n1, n2
! working variables:      
  INTEGER                :: i, j, ir
      
  i = MAX(ABS(n1),ABS(n2))
  j = MIN(ABS(n1),ABS(n2))
  
  IF (j == 0) THEN
     igcd = i
     RETURN
  END IF
  
  DO
     ir = MOD(i,j)
     IF (ir == 0) THEN
        igcd = j
        RETURN
     ELSE
        i = j
        j = ir
     END IF
  END DO
  
END FUNCTION igcd
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION vecLength(n,V)
!===============================================================================
! Calculate length of n-dimensional vector V in cartesian coordinates
!-------------------------------------------------------------------------------
! Input        :
!  n             dimension of vector V
!  V(n)          cartesian coordinates of vector V
! Output:
!  vecLength     length of vector V = sqrt(dot_product(V,V))
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n
  REAL(8), DIMENSION(n), INTENT(in) :: V !(n)
  
! working variable
  INTEGER :: i

! function to be used (location: math.f90)
  REAL(8) :: pythagoras
      
  IF (n < 1) THEN
     WRITE (*,*) 'vLength err: n < 1 :', n
     STOP
  ELSE IF (n == 1) THEN
     vecLength = DABS(V(1))
  ELSE IF (n == 2) THEN
     vecLength = pythagoras(V(1),V(2))
  ELSE     
     vecLength = pythagoras(V(1),V(2))
     DO i = 3,n
        vecLength = pythagoras(vecLength,V(i))
     END DO
  END IF
  
END FUNCTION vecLength
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION pythagoras(a,b)
!===============================================================================
! Calculate hypotenuse without destructive overflow or underflow
!-------------------------------------------------------------------------------
! Input        :
!  a,b           two sides of a triangle
! Output:
!  pythagoras    sqrt(a^2 + b^2)
!===============================================================================
  IMPLICIT NONE

! input variables  
  REAL(8), INTENT(in)    :: a,b

! working variables
  REAL(8)                :: absa, absb
      
  absa = DABS(a)
  absb = DABS(b)
  
  IF (absa > absb) THEN
     pythagoras = absa * DSQRT(1.D0 + (absb/absa)**2)
  ELSE
     IF(absb <= EPSILON(1.D0))THEN
        pythagoras = 0.D0
     ELSE
        pythagoras = absb * DSQRT(1.D0 + (absa/absb)**2)
     END IF
  END IF

END FUNCTION pythagoras
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION chebev(a,b,c,m,x)
!===============================================================================
! Calculate Chebyshev function
!-------------------------------------------------------------------------------
! Input        :
!  a, b          starting and end points
!  c(m)          chebyshev polynomials
!  m             order of the polynomials
!  x             increment
! Output       :
!  chebev        chebyshev function up to mth order
!===============================================================================
  IMPLICIT NONE

! input variables  
  INTEGER, INTENT(IN)    :: m
  REAL(8), INTENT(IN)    :: a, b, x, c(m)

! working variables
  INTEGER :: j
  REAL(8) :: d, dd, sv, y, y2
  IF ((x-a)*(x-b) > 0.D0) THEN
     WRITE (*,*) 'chebev err: x not in range:', x, a, b
     STOP
  END IF

  d  = 0.D0
  dd = 0.D0
  y  = (2.D0*x - a - b) / (b - a)
  y2 = 2.D0*y

  DO j = m, 2, -1
     sv = d
     d  = y2*d - dd + c(j)
     dd = sv
  END DO

! evaluate chebyshev function
  chebev = y*d - dd + 0.5*c(1)

END FUNCTION chebev
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION erf(x)
!===============================================================================
! Calculate error function
!-------------------------------------------------------------------------------
! Input        :
!  x             any real number
! Output       :
!  erf           erf(x)
! Note         :
!  This function calls other functions and subroutines obtained from
!  Numerical Recipes, i.e.
!  - gammp
!  - gcf
!  - gser
!  - gammln
!===============================================================================
  IMPLICIT NONE

! input variable  
  REAL(8) :: x

! use function gammp
  REAL(8) :: gammp

  IF (x < 0.D0) THEN
     erf = -gammp(.5,x**2)
  ELSE
     erf = gammp(.5,x**2)
  ENDIF

END FUNCTION erf

!-------------------------------------------------------------------------------

REAL(8) FUNCTION gammp(a,x)
  
  IMPLICIT NONE

  REAL(8) :: a, x

! use functions gcf, gser
  REAL(8) :: gammcf, gamser, gln

  IF (x < 0.D0 .OR. a <= 0.D0) PAUSE 'bad arguments in gammp'

  IF (x < a + 1.D0) THEN
     CALL gser(gamser,a,x,gln)
     gammp = gamser
  ELSE
     CALL gcf(gammcf,a,x,gln)
     gammp = 1.D0 - gammcf
  ENDIF

END FUNCTION gammp

!-------------------------------------------------------------------------------

SUBROUTINE gcf(gammcf,a,x,gln)

  IMPLICIT NONE

  INTEGER :: ITMAX
  REAL(8) :: a,gammcf,gln,x,EPS,FPMIN
  PARAMETER (ITMAX=100,EPS=3.D-7,FPMIN=1.D-30)

! use function gammln
  INTEGER :: i
  REAL(8) :: an, b, c, d, del, h, gammln

  gln = gammln(a)
  b   = x + 1.D0 - a
  c   = 1.D0/FPMIN
  d   = 1.D0/b
  h   = d

  DO i = 1, ITMAX

     an = -i * (i-a)
     b  = b + 2.D0
     d  = an*d + b

     IF (DABS(d) < FPMIN) d = FPMIN
     c  = b + an/c

     IF (DABS(c) < FPMIN) c = FPMIN
     d   = 1.D0/d
     del = d*c
     h   = h*del

     IF (DABS(del-1.D0) < EPS) THEN
        gammcf = DEXP(-x + a*DLOG(x) - gln)*h
        RETURN
     END IF

  END DO

  PAUSE 'a too large, ITMAX too small in gcf'

END SUBROUTINE gcf

!-------------------------------------------------------------------------------

SUBROUTINE gser(gamser,a,x,gln)

  IMPLICIT NONE

  INTEGER :: ITMAX
  REAL(8) :: a, gamser, gln, x, EPS
  PARAMETER (ITMAX=100,EPS=3.D-7)

! use function gammln
  INTEGER :: n
  REAL(8) :: ap, del, sum, gammln
  
  gln = gammln(a)
  IF (x <= 0.D0) THEN
     IF (x < 0.D0) PAUSE 'x < 0 in gser'
     gamser = 0.D0
     RETURN
  END IF

  ap  = a
  sum = 1.D0/a
  del = sum

  DO n = 1, ITMAX

     ap  = ap + 1.D0
     del = del*x / ap
     sum = sum+del

     IF (DABS(del) < DABS(sum)*EPS) THEN
        gamser = sum*DEXP(-x+a*DLOG(x)-gln)
        RETURN
     END IF

  END DO

  PAUSE 'a too large, ITMAX too small in gser'

END SUBROUTINE gser

!-------------------------------------------------------------------------------

REAL(8) FUNCTION gammln(xx)

  IMPLICIT NONE

  REAL(8) :: xx
  INTEGER :: j
  REAL(8) :: ser, stp, tmp, x, y, cof(6)
  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
       24.01409824083091d0,-1.231739572450155d0,&
       .1208650973866179d-2,-.5395239384953d-5,&
       2.5066282746310005d0/
  x   = xx
  y   = x
  tmp = x + 5.5D0
  tmp = (x + 0.5D0)*DLOG(tmp) - tmp
  ser = 1.000000000190015d0

  DO j = 1,6
     y   = y + 1.D0
     ser = ser + cof(j)/y
  END DO
  gammln = tmp + DLOG(stp*ser/x)

END FUNCTION gammln
!*******************************************************************************
!*******************************************************************************
FUNCTION fexp(x)

  REAL(8), INTENT(in)    :: x
  REAL(8), PARAMETER     :: cutoff=-15.D0
      
  IF (x <= cutoff) THEN
     fexp = 0.D0
  ELSE
     fexp = EXP(x)
  END IF

END FUNCTION fexp
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION diracDelta(x1,y1,x2,y2,fwhm)
!===============================================================================
! Approximate Dirac delta function
!-------------------------------------------------------------------------------
!..diracDelta=Int( Lorentzian(alpha+beta*x), from y1 to y2 )/(y2-y1)
!..where y1=alpha+beta*x1 and y2=alpha+beta*x2
!===============================================================================
  IMPLICIT NONE

  REAL(8), INTENT(in)    :: x1, x2, y1, y2, fwhm      
  REAL(8), PARAMETER     :: pi  = 3.141592654D0
  REAL(8), PARAMETER     :: tol = 1.D-7
      
  REAL(8)                :: dk, bb, a1, a2
      
  dk = x2-x1
  bb = y2 - y1
                
  a1 = 2.D0*y1 / fwhm
  a2 = 2.D0*y2 / fwhm
  
  diracDelta = (DATAN(DBLE(a2)) - DATAN(DBLE(a1))) / (pi*bb)
  
  RETURN

END FUNCTION diracDelta
!*******************************************************************************
!*******************************************************************************
SUBROUTINE outProd(n1,V1,n2,V2,ldu,U)
!===============================================================================
! Calculate outer product matrix of two real vectors V1 and V2
!-------------------------------------------------------------------------------
! Input        :
!  n1            dimension of vector V1
!  V1(n1)        vector V1
!  n2            dimension of vector V2
!  V2(n2)        vector V2
!  ldu           leading dimension of U
! Output       :
!  U(n1,n2)      outer product matrix U(i,j) = V1(i) x V2(j)
!=======================================================================
  IMPLICIT NONE
      
! input variables
  INTEGER, INTENT(in)    :: n1,n2,ldu
  REAL(8), INTENT(in)    :: V1(n1)
  REAL(8), INTENT(in)    :: V2(n2)

! output variables
  REAL(8),  INTENT(out)   :: U(ldu,n2) !(n1,n2)
  
! working variables
  INTEGER                :: i,j
  
! safety check
  IF (n1 > ldu) THEN
     WRITE (*,*) 'Outer Product err: n1 > ldu: ', n1, ldu
     STOP
  END IF
  
  DO i = 1, n1
     DO j = 1, n2
        U(i,j) = V1(i)*V2(j)
     END DO
  END DO
  
END SUBROUTINE outProd
!*******************************************************************************
!*******************************************************************************
SUBROUTINE linArray(nx,xmin,xmax,xarray)
!===============================================================================
! Make a linearly spaced array from xarray(1) = xmin to xarray(nx) = xmax
!-------------------------------------------------------------------------------
! Input        :
!  nx   (int)    number of array points
!  xmin (real)   minimum value of the array
!  xmax (real)   maximum value of the array
! Output       :
!  xarray        linearly spaced array
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: nx
  REAL(8), INTENT(in)    :: xmin,xmax

! output variable
  REAL(8), INTENT(out)   :: xarray(nx)
      
! working variables
  REAL(8)                :: dx
  INTEGER                :: i
      
  dx = (xmax - xmin) / (DBLE(nx) - 1.D0)
  DO i = 1, nx
     xarray(i) = xmin + (DBLE(i) - 1.d0)*dx
  END DO

END SUBROUTINE linArray
!*******************************************************************************
!*******************************************************************************
SUBROUTINE dos1Dgauss(nb,nk,rka,ldenk,Enk,E,fwhm,nn,DSn)
!===============================================================================
! one dimensional density of states per unit length (eV*A)^(-1)
! NOTE: A spin degeneracy factor of two is NOT included
!-------------------------------------------------------------------------------
! Input        :
!  nb            number of bands
!  nk            number of k-points
!  rka(nk)       k points (1/A)
!  ldenk         leading dimension of Enk
!  Enk(nk,nb)    Energy bands (eV)
!  E             energy at which DOS is desired (eV)
!  fwhm          fwhm width for broadened delta function
!  nn            band index for desired density of states
! Output       :
!  DSn           density of states due to n-th band (1/eV)
!===============================================================================
  IMPLICIT NONE
      
! input variables
  INTEGER, INTENT(in)    :: nb, nk, ldenk, nn
  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(nk,nb)
  REAL(8), INTENT(in)    :: E,fwhm

! output variable
  REAL(8), INTENT(out)   :: DSn

! working variables and parameter
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0     
  REAL(8)                :: gamma, rgamma, ss, A, B, erf, erf1, erf2
  INTEGER                :: ier, k

! check input for errors
  ier = 0
  IF (nk > ldenk) THEN
     ier = 1
     WRITE (*,*) 'dos1Dgauss err: nk.gt.ldenk:', nk, ldenk 
  END IF
  IF (nn > nb) THEN
     ier=1
     WRITE (*,*) 'dos1Dgauss err: nn.gt.nb:   ', nn, nb  
  END IF
  IF (ier.NE.0) STOP
      
! gaussian broadening parameter (1/eV**2)
  gamma  = 2.77259D0 / fwhm**2
  rgamma = DSQRT(gamma)
      
! density of states
  ss = 0.D0
  DO k = 1, nk-1
      
     A = (rka(k+1)*Enk(k,nn)-rka(k)*Enk(k+1,nn))/(rka(k+1)-rka(k))-E
     B = (Enk(k+1,nn)-Enk(k,nn))/(rka(k+1)-rka(k))
        
     IF (B /= 0.D0) THEN
        erf2 = erf(rgamma*(Enk(k+1,nn)-E))
        erf1 = erf(rgamma*(Enk(k  ,nn)-E))
        ss   = ss + (erf2-erf1)/B
     ELSE
        ss   = ss+DSQRT(gamma/pi)*EXP(-gamma*A**2)*(rka(k+1)-rka(k)) 
     END IF
      
  END DO
  
  DSn = ss/(2.D0*pi)
      
END SUBROUTINE dos1Dgauss
!*******************************************************************************
!*******************************************************************************
SUBROUTINE solveHam(n,ldh,ham,ldo,ovlp,matz,il,iu,nout,w,ldz,z)
!===============================================================================
! Solves generalized hermitian eigenvalue problem Ham*z = w*Ovlp*z
! The eigenvectors are normalized so that conjg(z)*Ovlp*z = 1 and
! the sign of z is fixed so the largest real component is positive.
! LAPACK routine ZHEGV is used.
!-------------------------------------------------------------------------------
! Input        :
!  n             order of hermitian matrix (Ham and Ovlp are n x n)
!  ldh           leading dimension of ham
!  ham(ldh,n)    complex hermitian matrix to be diagonalized
!  ldo           leading dimension of ovlp
!  ovlp(ldo,n)   complex hermitian overlap matrix
!  matz          option flag (0=evalues, 1=evalues+evectors)
!  il,iu         upper and lower indices of desired eigenvalues
!                if il or iu <= 0, all eigenvalues are returned
! Output       :
!  nout          number of eigenvalues returned
!  w(n)          eigenvalues in ascending order
!  ldz           leading dimension of z
!  z(ldz,n)      complex eigenvectors if matz = 1
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER                :: n, ldh, ldo, matz, il, iu
  COMPLEX(8)             :: ham(ldh,n), ovlp(ldo,n)
  
! output variables
  INTEGER                :: nout, ldz
  REAL(8)                :: w(n)
  COMPLEX(8)             :: z(ldz,n)
  
! working variables      
  INTEGER                :: lwork, info, i, j, ier, itype
  REAL(8)                :: rez, rezmax, rwork(3*n-2), ww(n)
  COMPLEX(8)             :: a(n,n), b(n,n), work(2*n-1)
  CHARACTER              :: jobz, uplo
  
! check input for errors
  ier = 0
  IF (matz /= 0 .AND. matz /= 1) THEN
     ier = 1
     WRITE (*,*) 'solveHam err: invalid matz: ', matz
  END IF
  
  IF (n > ldh) THEN
     ier = 1
     WRITE (*,*) 'solveHam err: n.gt.ldh: ', n, ldh
  END IF
  
  IF (n > ldo) THEN
     ier = 1
     WRITE (*,*) 'solveHam err: n.gt.ldo: ', n, ldo
  END IF
  
  IF (matz == 1 .AND. n > ldz) THEN
     ier = 1
     WRITE (*,*) 'solveHam err: n.gt.ldz: ', n, ldz
  END IF
  
  IF (il > 0 .AND. iu > 0 .AND. il > iu) THEN
     ier = 1
     WRITE (6,*) 'solveHam err: il.gt.iu: ', il, iu
  END IF
  
  IF (iu > n) THEN
     ier = 1
     WRITE (6,*) 'solveHam err: iu.gt.n: ', iu, n
  END  IF
  
  IF  (ier /= 0) STOP
  
! number of eigenvalues returned      
  IF (il <= 0 .OR. iu <= 0) THEN
     nout = n
  ELSE
     nout = iu - il + 1
  END IF
  
! solve generalized hermitian eigenvalue problem using lapack routine zhegv
  itype = 1
  jobz  = 'N'
  IF (matz == 1) jobz='V'
  uplo  = 'U'
  
  DO i = 1, n
     DO j = 1, n
        a(i,j) = ham(i,j)
        b(i,j) = ovlp(i,j)
     END DO
  END DO
  
  lwork = MAX(1,2*n-1) 
  CALL zhegv(itype,jobz,uplo,n,a,n,b,n,ww,work,lwork,rwork,info)
  
  IF(info /= 0) THEN
     WRITE (*,*) 'solveHam err: zhegv returns info =',info
     STOP
  END IF
  
! eigenvalues and eigenvectors in the desired range
  IF (nout == n) THEN
     
     DO i = 1, n
        w(i) = ww(i)
     END DO
     
     IF (matz.EQ.1) THEN
        DO j = 1 ,n
           DO i = 1, n
              z(i,j) = a(i,j)
           END DO
        END DO
     END IF
     
  ELSE
     
     DO i = 1, nout
        w(i) = ww(i+il-1)
     END DO
     
     IF (matz == 1) THEN
        DO j = 1, nout
           DO i = 1, n
              z(i,j) = a(i,j+il-1)
           END DO
        END DO
     END IF
     
  END IF
  
! phase convention for eigenvectors
! biggest real part is positive
  IF (matz == 1) THEN
     DO i = 1, nout
        
        rezmax = 0.D0
        DO j = 1, n
           IF (DABS(DBLE(z(j,i))) > rezmax) THEN
              rez = DBLE(z(j,i))
              rezmax = DABS(DBLE(z(j,i)))
           END IF
        END DO
        
        IF (rez < 0.D0) THEN
           DO j = 1, n
              z(j,i) = -z(j,i)
           END DO
        END IF
        
     END DO
  END IF
  
  RETURN
  
END SUBROUTINE solveHam
!*******************************************************************************
!*******************************************************************************
SUBROUTINE zbrac(func,x1,x2,succes)
!===============================================================================
! Bracketing routine from numerical recipe
!===============================================================================
  IMPLICIT NONE

  REAL(8)                :: x1, x2
  REAL(8), EXTERNAL      :: func 
  REAL(8), PARAMETER     :: FACTOR = 1.6D0
  INTEGER, PARAMETER     :: NTRY = 50
  INTEGER                :: j 
  REAL(8)                :: f1, f2 
  LOGICAL                :: succes 
  
  IF (x1 == x2) PAUSE 'you have to guess an initial range in zbrac' 
  f1 = func(x1) 
  f2 = func(x2) 
  succes = .TRUE. 
  DO j = 1, NTRY 
     IF (f1*f2 < 0.D0) RETURN 
     IF (ABS(f1) < ABS(f2)) THEN 
        x1 = x1 + FACTOR*(x1-x2) 
        f1 = func(x1) 
     ELSE 
        x2 = x2 + FACTOR*(x2-x1) 
        f2 = func(x2) 
     END IF
  END DO
  succes = .FALSE. 
END SUBROUTINE zbrac
!*******************************************************************************
!*******************************************************************************
FUNCTION rtbis(func,x1,x2,xacc) 
!===============================================================================
! Bisection routine from numerical recipe
!===============================================================================
  IMPLICIT NONE
  
  REAL(8)                :: rtbis, x1, x2, xacc
  REAL(8), EXTERNAL      :: func
  INTEGER, PARAMETER     :: JMAX = 40 
  INTEGER                :: j 
  REAL(8)                :: dx, f, fmid, xmid 
  
  fmid = func(x2) 
  f    = func(x1) 
  IF (f*fmid >= 0.D0) PAUSE 'root must be bracketed in rtbis' 
  IF (f < 0.D0) THEN 
     rtbis = x1 
     dx    = x2-x1 
  ELSE 
     rtbis = x2 
     dx    = x1-x2 
  END IF
  
  DO j = 1, JMAX 
     dx = dx * .5D0
     xmid = rtbis + dx 
     fmid = func(xmid) 
     IF (fmid <= 0.D0) rtbis = xmid 
     IF (ABS(dx) < xacc .OR. fmid == 0.D0) RETURN 
  END DO

  PAUSE 'too many bisections in rtbis'

END FUNCTION rtbis
!*******************************************************************************
!*******************************************************************************
SUBROUTINE sort2(N,Ra,Rb)
!===============================================================================
! Sorting routine from numerical recipe
!===============================================================================
  IMPLICIT NONE

! Dummy arguments
  INTEGER                :: N
  REAL(8), DIMENSION(N)  :: Ra, Rb
  INTENT (IN) N
  INTENT (INOUT) Ra , Rb

! Local variables
  INTEGER                :: i, ir, j, l
  REAL(8)                :: rra, rrb
  
  l  = N/2 + 1
  ir = N
  DO WHILE ( .TRUE. )
     IF ( l>1 ) THEN
        l = l - 1
        rra = Ra(l)
        rrb = Rb(l)
     ELSE
        rra = Ra(ir)
        rrb = Rb(ir)
        Ra(ir) = Ra(1)
        Rb(ir) = Rb(1)
        ir = ir - 1
        IF ( ir==1 ) THEN
           Ra(1) = rra
           Rb(1) = rrb
           EXIT
        END IF
     END IF
     i = l
     j = l + l
     DO WHILE ( .TRUE. )
        IF ( .NOT..TRUE. ) THEN
           RETURN
        ELSE IF ( j<=ir ) THEN
           IF ( j<ir ) THEN
              IF ( Ra(j)<Ra(j+1) ) j = j + 1
           END IF
           IF ( rra<Ra(j) ) THEN
              Ra(i) = Ra(j)
              Rb(i) = Rb(j)
              i = j
              j = j + j
           ELSE
              j = ir + 1
           END IF
           CYCLE
        END IF
        Ra(i) = rra
        Rb(i) = rrb
        GOTO 100
     END DO
     EXIT
100 END DO

END SUBROUTINE sort2
!*******************************************************************************
!*******************************************************************************
SUBROUTINE indexx(n,arr,indx)
!===============================================================================
! Indexing routine from numerical recipe
!===============================================================================
  IMPLICIT NONE

  INTEGER, PARAMETER     :: M = 7, NSTACK = 50
  INTEGER                :: n, indx(n)
  REAL(8)                :: arr(n), a
  
  INTEGER                :: i, indxt, ir, itemp, j, jstack, k, l
  INTEGER                :: istack(NSTACK)

  DO j = 1, n
     indx(j)=j
  END DO
  
  jstack=0
  l  = 1
  ir = n

1 IF (ir-l < M) THEN
     
     DO j = l+1, ir
        indxt = indx(j)
        a = arr(indxt)
        DO i = j-1, 1, -1
           IF (arr(indx(i)) <= a) GOTO 2
           indx(i+1)=indx(i)
        END DO
        i=0
2       indx(i+1) = indxt
     END DO
     
     IF (jstack == 0) RETURN
     ir = istack(jstack)
     l  = istack(jstack-1)
     jstack = jstack-2
  
  ELSE
     k = (l+ir)/2
     itemp     = indx(k)
     indx(k)   = indx(l+1)
     indx(l+1) = itemp
     
     IF (arr(indx(l+1)) > arr(indx(ir))) THEN
        itemp     = indx(l+1)
        indx(l+1) = indx(ir)
        indx(ir)  = itemp
     END IF
     
     IF (arr(indx(l)) > arr(indx(ir))) THEN
        itemp     = indx(l)
        indx(l)   = indx(ir)
        indx(ir)  = itemp
     END IF
     
     IF (arr(indx(l+1)) > arr(indx(l))) THEN
        itemp     = indx(l+1)
        indx(l+1) = indx(l)
        indx(l)   = itemp
     ENDIF
     i = l+1
     j = ir
     indxt = indx(l)
     a = arr(indxt)
3    CONTINUE
     i = i + 1
     IF (arr(indx(i)) < a) GOTO 3
4    CONTINUE
     j = j-1
     IF (arr(indx(j)) > a) GOTO 4
     IF (j < i) GOTO 5
     itemp   = indx(i)
     indx(i) = indx(j)
     indx(j) = itemp
     GOTO 3
5    indx(l) = indx(j)
     indx(j) = indxt
     jstack  = jstack+2
     
     IF (jstack > NSTACK) THEN
        WRITE (*,*) 'NSTACK too small in indexx'
        WRITE (*,*) 'NSTACK:', NSTACK
        STOP
     END IF
        
     IF (ir-i+1.GE.j-l) THEN
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir = j-1
     ELSE
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l = i
     END IF
  END IF

 GOTO 1

END SUBROUTINE indexx
!*******************************************************************************
!*******************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Added by Daria Satco !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION diracDelta_eps1(x1,y1,x2,y2,fwhm)
!===============================================================================
! Approximate Dirac delta function
!-------------------------------------------------------------------------------
!..diracDelta=Int( f(alpha+beta*x), from y1 to y2 )/(y2-y1)
!..where y1=alpha+beta*x1 and y2=alpha+beta*x2
!.. f = x/( x**2 + fwhm**2 )
!===============================================================================
  IMPLICIT NONE

  REAL(8), INTENT(in)    :: x1, x2, y1, y2, fwhm
  REAL(8), PARAMETER     :: pi  = 3.141592654D0
  REAL(8), PARAMETER     :: tol = 1.D-1

  REAL(8)                :: dk, bb, a1, a2

  dk = x2-x1
  bb = (y2 - y1)

  a1 = y1**2 + fwhm**2
  a2 = y2**2 + fwhm**2
  diracDelta_eps1 = DLOG(ABS(a2)/ABS(a1)) / (2*bb)

  RETURN

END FUNCTION diracDelta_eps1

!*******************************************************************
