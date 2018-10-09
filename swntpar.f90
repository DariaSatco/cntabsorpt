!*******************************************************************************
!*******************************************************************************
! Project      : tubepar.f90
!===============================================================================
! Purpose      :
! Library to generate important nanotube parameters
!-------------------------------------------------------------------------------
! Authors      :
! - ART Nugraha (nugraha@flex.phys.tohoku.ac.jp)
! - Daria Satco  (dasha.shatco@gmail.com)
! Latest Vers. : 2018.10.02
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Daria Satco added (autumn 2018) !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! - SUBROUTINE CutLineK(n,m,muk11,muk22,mukp11,mukp22)
! - SUBROUTINE solvePQeq(n,m,p,q)
! - SUBROUTINE CutLineii(n,m,nhex,muii)
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
!===============================================================================
! sends to the terminal the main information about CNT
! ------------------------------------------------------------------------------
! Input      :
!  n,m           chiral vector coordinates in (a1,a2)
! iunit          unit connected to the output file
! Output     : none
!===============================================================================
  IMPLICIT NONE
  INTEGER, INTENT(in)    :: n, m, iunit
  
  REAL(8), PARAMETER     :: hbar = 6.582D-4 !(eV-ps)
  REAL(8), PARAMETER     :: pi   = 3.14159265358979D0
  
  INTEGER                :: id, idr, igcd, metal, nhex, nHexagon, i
  REAL(8)                :: diam, tubeDiam
  REAL(8)                :: Tlength, trLength
  REAL(8)                :: chTheta
  INTEGER                :: muk11,muk22,mukp11,mukp22
  
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
  
  CALL CutLineK(n,m,muk11,muk22,mukp11,mukp22)
  ! CutLineK gives the result assuming 0 .. N-1 range for cutting lines
  ! in the code we use mu = 1 .. N
  ! the point N in our calculations is equal to 0 point in another notation

! metallicity
  SELECT CASE (metal)
     
  CASE(0)
     IF (id == idr) THEN
        WRITE (i,*) ' Metallicity          : Metal-1 (M0)'
        WRITE (i,*) ' Cutting lines corresponding to K points: ', muk11, mukp11
     ELSE
        WRITE (i,*) ' Metallicity          : Metal-2 (M0)'
        WRITE (i,*) ' Cutting lines corresponding to K point: ', muk11
     END IF
     

  CASE(1)
     WRITE (i,*) ' Metallicity          : Mod 1 Semiconductor (S2)'
     WRITE (i,*) ' Cutting lines corresponding to E11, E22 near K point: ', muk11, muk22
     WRITE (i,*) ' Cutting lines corresponding to E11, E22 near Kp point: ', mukp11, mukp22
     
  CASE(2)
     WRITE (i,*) ' Metallicity          : Mod 2 Semiconductor (S1)'
     WRITE (i,*) ' Cutting lines corresponding to E11, E22 near K point: ', muk11, muk22
     WRITE (i,*) ' Cutting lines corresponding to E11, E22 near Kp point: ', mukp11, mukp22
     
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
!*******************************************************************************
SUBROUTINE CutLineK(n,m,muk11,muk22,mukp11,mukp22)
!===============================================================================
! calculates the indeces of the closest to the K point cutting lines
! the analytical expressions are taken from:
! Saito, R., et al.
! "Cutting lines near the Fermi energy of single-wall carbon nanotubes."
! Physical Review B 72.15 (2005): 153413.
! ************
! the cutting lines are counted from 0 to N-1
! ************
! ------------------------------------------------------------------------------
! Input      :
!  n,m             chiral vector coordinates in (a1,a2)
! Output     :
! muk11, muk22
! Semiconducting tubes    : cutting line indeces for E11 and E22 transitions near the K point
! Metallic tubes          : muk11 - cutting line index for Fermi energy level near K point, muk22 is put to zero

! mukp11, mukp22
! Semiconducting tubes    : cutting line indeces for E11 and E22 transitions near the K' point
! Metallic tubes          : muk11 - cutting line index for Fermi energy level near K' point, mukp22 is put to zero
!===============================================================================
  IMPLICIT NONE
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(out)   :: muk11, muk22, mukp11, mukp22

  REAL(8), PARAMETER     :: hbar = 6.582D-4 !(eV-ps)
  REAL(8), PARAMETER     :: pi   = 3.14159265358979D0

  INTEGER                :: id, idr, igcd, metal, nhex, nHexagon
  INTEGER                :: it1, it2, p, q
  REAL(8)                :: T(2)
  INTEGER                :: num, denum

  id      = igcd(n,m)
  idr     = igcd(2*n+m,2*m+n)
  metal   = MOD(n-m,3)

  nhex    = nHexagon(n,m)

  CALL trVecXY(n,m,it1,it2,T)

! metallic
  IF ( metal == 0 ) THEN
     IF (id == idr) THEN
     ! Metal-1 CNT
         muk11 = nhex/3
         mukp11 = 2*nhex/3
     ELSE
     ! Metal-2 CNT
        CALL solvePQeq(n,m,p,q)

        num = (3*p + 2)*nhex
        denum = -3*it2

        IF (MOD(3*m/idr, 3) == 1) THEN
            muk11 = (num - m)/denum
        ELSEIF (MOD(3*m/idr, 3) == 2)  THEN
            muk11 = (num + m)/denum
        ELSE
            WRITE(*,*) 'The eq.13 is wrong'
        END IF

        mukp11 = nhex - muk11
     END IF

     muk22  = 0
     mukp22 = 0

  END IF


  IF ( metal == 1 .or. metal == 2 ) THEN
      IF ( MOD(nhex,3) == 1 ) THEN

            muk11 = (nhex - 1)/3
            muk22 = (nhex + 2)/3

            mukp11 = nhex - (nhex - 1)/3
            mukp22 = nhex - (nhex + 2)/3

      ELSEIF ( MOD(nhex,3) == 2 ) THEN

            muk11 = (nhex + 1)/3
            muk22 = (nhex - 2)/3

            mukp11 = nhex - (nhex + 1)/3
            mukp22 = nhex - (nhex - 2)/3
      END IF
  ENDIF

END SUBROUTINE CutLineK
!*******************************************************************************
!*******************************************************************************
SUBROUTINE solvePQeq(n,m,p,q)
!===============================================================================
! solves the eq.(13) from:
! Saito, R., et al.
! "Cutting lines near the Fermi energy of single-wall carbon nanotubes."
! Physical Review B 72.15 (2005): 153413.
! ------------------------------------------------------------------------------
! Input      :
!  n,m             chiral vector coordinates in (a1,a2)
! Output     :
! p,q              solution
!===============================================================================
  IMPLICIT NONE
  INTEGER, INTENT(in)    :: n,m
  INTEGER, INTENT(out)   :: p,q

  INTEGER                :: it1, it2, idr, igcd
  REAL(8)                :: T(2)         ! (T(1),T(2)) -> (x,y)
  INTEGER                :: iq, ip
  INTEGER                :: cond

  CALL trVecXY(n,m,it1,it2,T)
  idr = igcd(2*n+m,2*m+n)

  p = 0
  q = 0

  IF (MOD(3*m/idr, 3) == 1) THEN
    cond = (-3*m + idr)/(3*idr)
  ELSEIF (MOD(3*m/idr, 3) == 2)  THEN
    cond = (-3*m - idr)/(3*idr)
  ELSE
    WRITE(*,*) 'The eq.13 is wrong'
  END IF

  !PRINT*, 'cond =', cond
  DO ip = 0, -it2-1
    DO iq = 0, it1-1

        IF ( (it1*ip + it2*iq) == cond ) THEN
            p = ip
            q = iq
        END IF

    END DO
  END DO

END SUBROUTINE solvePQeq
!*******************************************************************************
!*******************************************************************************
SUBROUTINE CutLineii(n,m,nhex,muii)
!===============================================================================
! calculates the indeces of cutting lines, which correspond to particular transitions
! ------------------------------------------------------------------------------
! Input      :
!  n,m             chiral vector coordinates in (a1,a2)
!  nhex            number of hexagons
! Output     :
! muii             array of cutting line ideces which correspond to certain Eii (4,nhex/2 + 1)
! comments on dimensionality: for semiconducting-1 and -2 CNTs as well as for metal-1 the degeneracy of
! energy bands is equal to 2, then muii is non zero only in range (1:2, :)
! for metal-2 CNT energy bands are four fold degenerate, then non zero values are (1:4, 1: ..)
!===============================================================================
IMPLICIT NONE

!input variables
  INTEGER, INTENT(in)    :: n,m, nhex

!output variables
  INTEGER, INTENT(out)   :: muii(4,nhex/2+1)

!working variables
  INTEGER                :: muk11, muk22, mukp11, mukp22
  INTEGER                :: even, odd
  INTEGER                :: mu11, mu22, musize
  INTEGER                :: i, dk, dk1, dk2, max_position(1), mpos
  INTEGER                :: metal, cond, denom
  INTEGER                :: id, idr, igcd

  id      = igcd(n,m)
  idr     = igcd(2*n+m,2*m+n)

  CALL CutLineK(n,m,muk11,muk22,mukp11,mukp22)
  muii = 0

  metal   = MOD(n-m,3)

  even = 2
  odd = 1

  IF ( metal == 0 ) THEN
! consider metallic  CNT
      dk1 = 1
      dk2 = -1

      IF ( muk11 < nhex/2 ) THEN
          mu11 = muk11
      ELSE
          mu11 = nhex - muk11
      END IF

      muii(1,1) = mu11
      muii(2,1) = nhex - muii(1,1)

      muii(3,1) = mu11
      muii(4,1) = nhex - muii(1,1)

      ! check if we first reach the center or the edge
      IF ( mu11 > INT(nhex/4) .and. mu11 < nhex/2 ) THEN
        cond = 1
        denom = 2
      ELSEIF ( mu11 < INT(nhex/4) ) THEN
        cond = 4
        denom = 1
      ELSEIF ( mu11 == nhex/2 ) THEN
        cond = 1
        denom = 1
      ELSE
        PRINT*, "No condition was found"
      END IF

!      PRINT*, cond
!      PRINT*, muii(1,1), muii(2,1), muii(3,1), muii(4,1)

      i = 1
      DO WHILE ( MAXVAL(muii(cond,:)) .lt. nhex/denom )
          i = i + 1
          muii(1,i) = muii(1,i-1) + dk1
          muii(2,i) = nhex - muii(1,i)

          muii(3,i) = muii(3,i-1) + dk2
          muii(4,i) = nhex - muii(3,i)
          !PRINT*, muii(1,i), muii(2,i), muii(3,i), muii(4,i), 'max', MAXVAL(muii(cond,:))
      END DO

      PRINT*, i
      IF ( i .ge. nhex/2+1 ) THEN
        GOTO 1
      ELSE
          IF ( cond == 1 ) THEN
            dk = -1
            muii(1,i+1) = MINVAL(muii(1,1:i)) + dk
            muii(2,i+1) = nhex - muii(1,i+1)
          ELSE
            dk = 1
            muii(1,i+1) = MAXVAL(muii(1,1:i)) + dk
            muii(2,i+1) = nhex - muii(1,i+1)
          END IF

          i = i + 1
          DO WHILE ( i < nhex/2+1 )
              IF ( muii(1,i) == muii(2,i) ) EXIT

              muii(1,i+1) = muii(1, i) + dk
              muii(2,i+1) = nhex - muii(1,i+1)
              i = i + 1
              !PRINT*, i, muii(1,i), muii(2,i)
          END DO
      END IF

  ELSE
!consider semiconducting CNT

! fix ordering of cutting lines numbers
! such as muii(1,:) <-> 1, nhex/2
! muii(2,:) <-> nhex/2, nhex
      IF (muk11 < mukp11) THEN
          mu11 = muk11
          mu22 = muk22
      ELSE
          mu11 = mukp11
          mu22 = mukp22
      END IF

! muii(:,k) <-> Ekk
      IF (mu11 > mu22) THEN
          dk1 = 1
          dk2 = -1
      ELSE
          dk1 = -1
          dk2 = 1
      END IF

! check if we first reach the center or the edge
      IF (mu11 > INT(nhex/4)) THEN
        cond = 1
      ELSE
        cond = 2
      END IF

! assign starting values
      muii(1,1) = mu11
      muii(2,1) = nhex - mu11

      muii(1,2) = mu22
      muii(2,2) = nhex - mu22

      DO WHILE ( MAXVAL(muii(cond,:)) .lt. nhex/(2-(cond-1)) )
          even = even + 2
          odd = odd + 2

          muii(1,odd) = muii(1,odd-2) + dk1
          muii(2,odd) = nhex - muii(1,odd)
          !PRINT*, muii(1,odd)
          muii(1,even) = muii(1,even-2) + dk2
          muii(2,even) = nhex - muii(1,even)
          !PRINT*, muii(1,even), 'max', MAXVAL(muii(1,1:nhex/2+1))
      END DO

      max_position = MAXLOC(muii(1,:))
      !PRINT*, max_position(1)

      mpos = max_position(1)
      IF ( cond == 1 ) THEN
        dk = -1
        muii(1,mpos+1) = MINVAL(muii(1,1:mpos)) + dk
        muii(2,mpos+1) = nhex - muii(1,mpos+1)
      ELSE
        dk = 1
        muii(1,mpos+1) = MAXVAL(muii(1,1:mpos)) + dk
        muii(2,mpos+1) = nhex - muii(1,mpos+1)
      END IF

      DO i = mpos + 1, nhex/2
          muii(1,i+1) = muii(1, i) + dk
          muii(2,i+1) = nhex - muii(1,i+1)
      END DO

   END IF

1   IF (MINVAL(muii) < 0) STOP "WARNING: Negative value in cutting line"

    PRINT*, "Cutting lines numbers matched with ii transitions"
    IF ( metal == 0 ) THEN
        DO i = 1, nhex/2+1
            PRINT*, i-1,i-1, muii(1,i), muii(2,i), muii(3,i), muii(4,i)
        END DO
    ELSE
        DO i = 1, nhex/2+1
            PRINT*, i,i, muii(1,i), muii(2,i)
        END DO
    END IF

 END SUBROUTINE Cutlineii
 !*******************************************************************************


