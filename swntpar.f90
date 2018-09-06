!*******************************************************************************
!*******************************************************************************
! Project      : tubepar.f90
!===============================================================================
! Purpose      :
! Library to generate important nanotube parameters
!-------------------------------------------------------------------------------
! Author       : ART Nugraha (nugraha@flex.phys.tohoku.ac.jp)
! Latest Vers. : 2012.11.05
!-------------------------------------------------------------------------------
! Reference    :
! Physical Properties of Carbon Nanotubes
! R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
!-------------------------------------------------------------------------------
! Contents     :
! - FUNCTION nHexagon(n,m)
! - FUNCTION chTheta(n,m)
! - FUNCTION tubeDiam(n,m)
! - FUNCTION chLength(n,m)
! - FUNCTION trLength(n,m)
! - SUBROUTINE chVecXY(n,m,Ch)
! - SUBROUTINE trVecXY(n,m,it1,it2,T)
! - SUBROUTINE unitVecXY(a1,a2)
! - SUBROUTINE printTubeClass(n,m,iunit)
!*******************************************************************************
!*******************************************************************************
INTEGER FUNCTION nHexagon(n,m)
!===============================================================================
! Calculate number of hexagons in (n,m) carbon nanotube unit cell
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  nHexagon      number of hexagons
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n,m
   
! working variables
  INTEGER, SAVE          :: nch   = 0
  INTEGER, SAVE          :: mch   = 0
  INTEGER, SAVE          :: nhex  = 0 
  INTEGER :: idr, igcd
  
  IF (n <= 0 .OR. m > n) THEN
     WRITE (*,*) 'nHexagon err: invalid n,m: ', n, m
     STOP
  END IF
      
  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m

! calculate d_R = gcd(2n+m,2m+n)
     idr   = igcd(2*n+m,2*m+n)
! calculate number of hexagons in a tube unit cell
     nhex  = 2*(n**2 + n*m + m**2) / idr     
  END IF

! return nHexagon
  nHexagon = nhex

END FUNCTION nHexagon
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION chTheta(n,m)
!===============================================================================
! Calculate chiral angle of (n,m) nanotube (radians)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  chTheta       tube chiral angle (radians)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  
! working variables
  INTEGER, SAVE          :: nch   = 0
  INTEGER, SAVE          :: mch   = 0
  REAL(8), SAVE          :: theta = 0.D0
  REAL(8), PARAMETER     :: sqrt3 = 1.732050807568877D0
  REAL(8)                :: rm, rn, sint

  IF (n <= 0 .OR. m > n) THEN
     WRITE (*,*) 'chTheta err: invalid n,m: ', n, m
     STOP
  END IF
  
  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m
     
! real type of input variables
     rm = DBLE(m)
     rn = DBLE(n)

! calculate theta
     sint  = sqrt3 * rm / ( 2.D0 * DSQRT(rn**2 + rn*rm + rm**2) )
     theta = DASIN(sint)
  END IF
  
! return chiral angle
  chTheta = theta

END FUNCTION chTheta
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION tubeDiam(n,m)
!===============================================================================
! Calculate diameter of (n,m) nanotube (angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  tubeDiameter  nanotube diameter (A)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER,  INTENT(in)   :: n, m

! working variables and parameters
  REAL(8), PARAMETER     :: pi  = 3.141592653589793D0
  INTEGER, SAVE          :: nch = 0
  INTEGER, SAVE          :: mch = 0
  REAL(8), SAVE          :: diameter = 0.D0
      
! variable for calling function chLength
  REAL(8)                :: chLength

  IF (n <= 0 .OR. m > n) THEN
     WRITE (*,*) 'tubeDiameter err: invalid n,m: ', n, m
     STOP
  END IF
      
  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m
! dt = |Ch| / pi (circumference formula)
     diameter = chLength(n,m) / pi
  END IF

! return tube diameter      
  tubeDiam = diameter            
      
END FUNCTION tubeDiam
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION chLength(n,m)
!===============================================================================
! Calculate length of (n,m) chiral vector (angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  chLength      length of chiral vector (A)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m

! working variables
  REAL(8), DIMENSION(3)  :: Ch     
  INTEGER, SAVE          :: nch     = 0
  INTEGER, SAVE          :: mch     = 0
  REAL(8), SAVE          :: vlength = 0
      
! function called from math.f90
  REAL(8) :: vecLength

  IF (n <= 0 .OR. m > n) THEN
     WRITE (*,*) 'chLength err: invalid n,m: ', n, m
     STOP
  END IF

  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m

! define chiral vector (input: n, m and output: Ch)
     CALL chVecXY(n,m,Ch)

! calculate vector length
     vlength = vecLength(2,Ch)
  END IF

! return chiral vector length
  chLength = vlength
      
END FUNCTION chLength
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION trLength(n,m)
!===============================================================================
! Calculate Length of Translation vector in (n,m) nanotube (angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  trLength      length of nanotube translation vector (A)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m

! working variables
  REAL(8), DIMENSION(2)  :: T !(2)
  INTEGER, SAVE          :: nch = 0
  INTEGER, SAVE          :: mch = 0
  REAL(8), SAVE          :: tlength = 0.D0
  INTEGER                :: it1, it2

! function called from math.f90
  REAL(8) ::  vecLength

  IF (n <= 0 .OR. m > n) THEN
     WRITE (*,*) 'trLength err: invalid n,m: ', n, m
     STOP
  END IF
     
  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m

! define translation vector
     CALL trVecXY(n,m,it1,it2,T)

! calculate translation vector length
     tlength = vecLength(2,T)      
  END IF

! return length of translation vector      
  trLength = tlength      
      
END FUNCTION trLength
!*******************************************************************************
!*******************************************************************************
SUBROUTINE chVecXY(n,m,Ch)
!===============================================================================
! Define chiral vector (n,m) in xy coordinates (angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  Ch(2)         xy coordinates of chiral vector (A)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)  :: n,m

! output variables
  REAL(8), INTENT(out) :: Ch(2)
      
! working variables
  REAL(8) :: a1(2), a2(2)     
  INTEGER, SAVE        ::  nch = 0
  INTEGER, SAVE        ::  mch = 0
  REAL(8), SAVE        ::  ChSave(2)
      
  IF (n <= 0 .OR. m > n) THEN
     WRITE (*,*) 'chVecXY err: invalid n,m: ', n, m
     STOP
  END IF

  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m
     
     CALL unitVecXY(a1,a2)
     ChSave(1) = DBLE(n)*a1(1) + DBLE(m)*a2(1)
     ChSave(2) = DBLE(n)*a1(2) + DBLE(m)*a2(2)
  END IF

! return chiral vector in xy coordinates  
  Ch(1) = ChSave(1)
  Ch(2) = ChSave(2)
  
END SUBROUTINE chVecXY
!*******************************************************************************
!*******************************************************************************
SUBROUTINE trVecXY(n,m,it1,it2,T)
!===============================================================================
! Define translation vector T parallel to nanotube axis
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
! Output       :
!  it1,it2       translation vector in (a1,a2)
!  T(2)          translation vector in xy coordinates (A)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m

! output variables
  INTEGER, INTENT(out)   :: it1, it2
  REAL(8), INTENT(out)   :: T(2)         ! (T(1),T(2)) -> (x,y)

! working variables      
  REAL(8)                :: a1(2), a2(2) ! graphene unit vectors
  INTEGER                :: idr,igcd
      
  CALL unitVecXY(a1,a2)
            
  idr = igcd(2*n+m, 2*m+n)
      
! define t1 and t2
  it1 = (2*m+n)/idr
  it2 =-(2*n+m)/idr
      
! calculate the coordinates
  T(1) = DBLE(it1)*a1(1) + DBLE(it2)*a2(1)
  T(2) = DBLE(it1)*a1(2) + DBLE(it2)*a2(2)      
      
END SUBROUTINE trVecXY
!*******************************************************************************
!*******************************************************************************
SUBROUTINE unitVecXY(a1,a2)
!===============================================================================
! Calculate graphene unit vectors in xy coordinates (angstroms)
!-------------------------------------------------------------------------------
! Input        : none
! Output       :
!  Ch(2)         xy coordinates of chiral vector (A)
!===============================================================================
  IMPLICIT NONE

  REAL(8), PARAMETER     :: a = 2.49 ! angstrom, graphene lattice constant
  REAL(8), PARAMETER     :: sqrt3 = 1.732050807568877D0
            
  REAL(8) :: a1(2), a2(2)
      
  a1(1) =  a*sqrt3 / 2.D0
  a1(2) =  a/2.D0
        
  a2(1) =  a*sqrt3 / 2.D0
  a2(2) = -a / 2.D0
  
END SUBROUTINE unitVecXY
!*******************************************************************************
!*******************************************************************************
SUBROUTINE printTubeClass(n,m,iunit)
  IMPLICIT NONE
  INTEGER, INTENT(in)    :: n, m, iunit
  
  REAL(8), PARAMETER     :: hbar = 6.582D-4 !(eV-ps)
  REAL(8), PARAMETER     :: pi   = 3.14159265358979D0
  REAL(8)                :: Em(6)
  
  INTEGER                :: id, idr, igcd, metal, nhex, nHexagon, mu, i
  REAL(8)                :: diam, tubeDiam
  REAL(8)                :: Tlength, trLength
  REAL(8)                :: chTheta
  
  id      = igcd(n,m)
  idr     = igcd(2*n+m,2*m+n)
  metal   = MOD(n-m,3)
  
  nhex    = nHexagon(n,m)
  diam    = tubeDiam(n,m)
  Tlength = trLength(n,m)
   
! tube classification
  i = iunit
  
  WRITE (i,*) '              Tube classification'
  WRITE (i,*) ' Chiral indices (n,m) :', n, m
  
  IF (n == m) THEN
     WRITE (i,*) ' Class                : Armchair'
  ELSE IF (m == 0) THEN
     WRITE (i,*) ' Class                : Zigzag'
  ELSE
     WRITE (i,*) ' Class                : Chiral'      
  ENDIF
  
! metallicity
  SELECT CASE (metal)
     
  CASE(0)
     IF (id == idr) THEN
     WRITE (i,*) ' Metallicity          : Metal-1 (M0)'      
     ELSE
     WRITE (i,*) ' Metallicity          : Metal-2 (M0)'      
     END IF
     
  CASE(1)
     WRITE (i,*) ' Metallicity          : Mod 1 Semiconductor (S2)'
     
  CASE(2)
     WRITE (i,*) ' Metallicity          : Mod 2 Semiconductor (S1)'      
     
  CASE default
     WRITE (i,*) ' invalid metallicity  : ', metal
     
  END SELECT
  
  WRITE (i,*) ' Family (2 n + m)     : ', 2*n+m
  WRITE (i,*)
  WRITE (i,*) ' Chiral angle (Radian): ', chTheta(n,m)
  WRITE (i,*)
  WRITE (i,*) ' Tube diameter (A)    : ', diam
  WRITE (i,*) ' Unit cell length (A) : ', Tlength
  WRITE (i,*) ' Unit cell hexagons   : ', nhex                  
  WRITE (i,*) ' Unit cell atoms      : ', 2*nhex
  WRITE (i,*)      
            
END SUBROUTINE printTubeClass
!*******************************************************************************