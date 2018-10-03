!*******************************************************************************
!*******************************************************************************
! Project      : tubestruct.f90
!===============================================================================
! Purpose      :
! Library to generate nanotube structures and related vector operations
!-------------------------------------------------------------------------------
! Authors      :
! - ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
! - Daria Satco  (dasha.shatco@gmail.com)
! Latest Vers. : 2018.09.30
!-------------------------------------------------------------------------------
! Reference(s) :
! Physical Properties of Carbon Nanotubes
! R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
!-------------------------------------------------------------------------------
! Contents     :
! - FUNCTION NNatom(iatom,nn)
! - SUBROUTINE phij(n,m,iatom,ivec,nn,q,mu,phi)
! - SUBROUTINE phaseFactor(n,m,j1,j2,q,mu,phi)
! - SUBROUTINE rxyzVec(n,m,iatom,ivec,nn,Rxyz)
! - SUBROUTINE rxyzJ1J2Vec(n,m,iatom,j1,j2,Rxyz)
! - SUBROUTINE NNj1j2(iatom,ivec,nn,j1,j2)
! - SUBROUTINE tubetoXYZ(n,m,Vtube,Vxyz)
! - SUBROUTINE XYtoTube(n,m,Vin,Vout)
! - SUBROUTINE tauAjVecXY(ivec,nn,tau)
! - SUBROUTINE tauBjVecXY(ivec,nn,tau)
! - SUBROUTINE rhoVec(n,m,iatom,rho)
! - SUBROUTINE rho2Vec(n,m,iatom,ivec,rho)
! - SUBROUTINE tubeRadUnitVec(Vxyz,Ur)
! - SUBROUTINE r2xyzVec(n,m,iatom,ivec1,ivec2,Rxyz)
! - SUBROUTINE drVec(n,m,iatom,Dr)
! - SUBROUTINE dr2Vec(n,m,iatom,ivec,Dr)
! - SUBROUTINE dN2xyzVec(n,m,iatom,ivec1,ivec2,dN2xyz)
! - SUBROUTINE dNxyzVec(n,m,iatom,ivec,nn,dNxyz)
! - SUBROUTINE srmatrix(n,m,iatom,Sr)
! - SUBROUTINE sr2matrix(n,m,iatom,ivec,Sr)
! - SUBROUTINE sqmatrix(n,m,iatom,ivec,nn,q,mu,Sq)
! - SUBROUTINE sjmatrix(n,m,iatom,ivec,nn,Sj)
! - SUBROUTINE reducedCutLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Daria Satco added (autumn 2018) !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! - SUBROUTINE getHexagonPosition(n,m,nhex,j1j2)
!*******************************************************************************
!*******************************************************************************
INTEGER FUNCTION NNatom(iatom,nn)
!===============================================================================
! label neighboring atoms
!===============================================================================
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    :: iatom, nn
      
  IF(nn.EQ.0 .OR. nn.EQ.2) THEN
     NNatom = iatom
  ELSE
     IF(iatom.EQ.1) NNatom = 2
     IF(iatom.EQ.2) NNatom = 1
  END IF
      
  RETURN

END FUNCTION NNatom
!*******************************************************************************
!*******************************************************************************
SUBROUTINE phij(n,m,iatom,ivec,nn,q,mu,phi)
!===============================================================================
! Compute the phase factor phi = Q \cdot Rj (dimensionless), where
! Rj is the center of the two atom unit cell to which the neighbor atom
! at R(n)_(iatom,ivec,nn) belongs and Q is the 2D Q vector that satisfies
! translational and rotational boundary conditions
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            neighbor index nn = 0,1,2,3,4
!  q             wavevector along tube axis (1/Angstroms)
!  mu            labels manifolds (0...N_hex-1)
! Output       :
!  phi           dot product of Q and Rj (dimensionless)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, mu, iatom, ivec, nn
  REAL(8), INTENT(in)    :: q

! output variable
  REAL(8), INTENT(out)   :: phi
  
! working variables
  INTEGER                :: j1,j2
  
! near neighbor two atom unit cell indices j1 and j2
  CALL NNj1j2(iatom,ivec,nn,j1,j2)  
    
! dot product of Q and Rj (dimensionless)      
  CALL phaseFactor(n,m,j1,j2,q,mu,phi)
  
END SUBROUTINE phij
!*******************************************************************************
!*******************************************************************************
SUBROUTINE phaseFactor(n,m,j1,j2,q,mu,phi)
!===============================================================================
! Compute the phase factor phi = Q.Rj (dimensionless) where Rj is the
! center of the two atom unit cell in unrolled tube coordinates and Q
! is the 2D Q vector that satisfies translational and rotational
! boundary conditions
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  j1,j2         center of two atom unit cell in (a1,a2)
!  q             wavevector along tube axis (1/Angstroms)
!  mu            labels manifolds (0...N_hex-1)
! Output       :
!  phi           dot product of Q and Rj (dimensionless)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n,m,j1,j2,mu
  REAL(8), INTENT(in)    :: q

! output variable
  REAL(8), INTENT(out)   :: phi
  
! working variables and parameters
  REAL(8), PARAMETER     :: a     = 2.49D0
  REAL(8), PARAMETER     :: pi    = 3.14159265358979D0
  REAL(8), PARAMETER     :: sqrt3 = 1.73205080756888D0
  
  REAL(8)                :: rn, rm, qq, rmu, rj1, rj2
  REAL(8)                :: theta1, theta2, tau1, tau2, theta, tau
  
  rn  = DBLE(n)
  rm  = DBLE(m)
  qq  = DBLE(q)      
  rmu = DBLE(mu)
  
  rj1 = DBLE(j1)
  rj2 = DBLE(j2)
      
! screw operator rotations for basis vectors a1 and a2 (radians)
  theta1 = (2.D0*rn + rm)*pi / (rn**2 + rn*rm + rm**2)
  theta2 = (2.D0*rm + rn)*pi / (rn**2 + rn*rm + rm**2)
      
! screw operator translations for basis vectors a1 and a2 (Angstroms)      
  tau1   =  a*sqrt3*rm / (2.D0 * DSQRT(rn**2 + rn*rm + rm**2))
  tau2   = -a*sqrt3*rn / (2.D0 * DSQRT(rn**2 + rn*rm + rm**2))
      
! screw operator taking two atom unit cell at origin to near neighbor
  theta  = rj1*theta1 + rj2*theta2
  tau    = rj1*tau1   + rj2*tau2
      
! phase factor (dimensionless)
  phi    = rmu*theta + qq*tau
          
END SUBROUTINE phaseFactor
!*******************************************************************************
!*******************************************************************************
SUBROUTINE rxyzVec(n,m,iatom,ivec,nn,Rxyz)
!===============================================================================
! Compute the xyz coordinates for A and B atoms and their neighbors
! R(nn)_(iatom,ivec). If nn=0, return R_A or R_B.
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            nearest neighbor index nn = 0,1,2,3,4
! Output       :
!  Rxyz(3)       R(nn)_(iatom,ivec) in xyz coordinates (Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables     
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variable
  REAL(8), INTENT(out)   :: Rxyz(3)
      
! working variables
  REAL(8)                :: tau0(2), tau(2), Rtube(2)
  INTEGER                :: ier, nnn, iiatom, iivec,ii
  
  INTEGER, SAVE, DIMENSION(0:4)         :: nvecs = (/ 1, 3, 6, 3, 6 /)
  REAL(8), SAVE, DIMENSION(2,6, 0:4, 3) :: Rsave !(iatom,ivec,nn, 3)
  INTEGER, SAVE                         :: nch  = 0
  INTEGER, SAVE                         :: mch  = 0          

! check input for errors
  ier = 0
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ier = 1
     WRITE (*,*) 'rxyzVec err: invalid iatom: ', iatom
  END IF
  IF (ier /= 0) STOP
      
! update Rxyz lookup table if (n,m) has changed
  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m
        
     Rsave = 0.D0
     DO iiatom = 1, 2
        DO nnn = 0, 4
           DO iivec = 1, nvecs(nnn)
      
! tube coordinates of A or B atom in two atom unit cell (Angstroms)
              CALL tauAjVecXY(1,1,tau0)
              IF (iiatom == 1) THEN
                 tau0 =-tau0/2.D0 !tauA = -tau(1)_(A1)/2 (graphene xy coord.)
              ELSE
                 tau0 = tau0/2.D0 !tauB =  tau(1)_(A1)/2
              END IF
              CALL XYtoTube(n,m,tau0,tau0)      

! if nnn = 0, return xyz coordinates of A or B in two atom unit cell
              IF (nnn == 0) THEN
                 Rtube = tau0
                 CALL tubetoXYZ(n,m,Rtube,Rxyz)

              ELSE             
! tube coordinates of neighbor atom (iatom,ivec,nn) (Angstroms)
                 IF (iiatom == 1) THEN
                    CALL tauAjVecXY(iivec,nnn,tau)
                 ELSE
                    CALL tauBjVecXY(iivec,nnn,tau)      
                 END IF
                 CALL XYtoTube(n,m,tau,tau)
                 Rtube = tau0+tau
      
! go from tube to xyz coordinates (Angstroms)            
                 CALL tubetoXYZ(n,m,Rtube,Rxyz)
              END IF
        
! save Rxyz (Angstroms)
              DO ii = 1,3
                 Rsave(iiatom,iivec,nnn,ii) = Rxyz(ii)
              END DO
      
           END DO
        END DO
     END DO
  END IF
      
! return Rxyz using lookup table (Angstroms)
  iivec = ivec
  IF(nn == 0) iivec = 1
  DO ii = 1,3
     Rxyz(ii) = Rsave(iatom,iivec,nn,ii)
  END DO
      
END SUBROUTINE rxyzVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE rxyzJ1J2Vec(n,m,iatom,j1,j2,Rxyz)
!===============================================================================
! compute xyz atomic position vector for an atom in the J = (j1,j2)
! hexagon in (n,m) nanotubes. The A and B atoms are iatom = 1 or 2.
! If iatom is otherwise, the center of the hexagon is returned.
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         atom in the unit cell (1=A,2=B)
!  j1,j2         center of two atom unit cell in a1 a2 basis
! Output       :
!  Rxyz(3)       atomic position in tube xyz coordinates (Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, j1, j2

! output variable
  REAL(8), INTENT(out)   :: Rxyz(3)

! working variables
  REAL(8)                :: a1(2), a2(2), tau(2), Rj(2)     
  REAL(8)                :: rj1, rj2
      
!..center of hexagonal unit cell in graphene coordinates (Angstroms)
  CALL unitVecXY(a1,a2)
  rj1 = float(j1)
  rj2 = float(j2)
  Rj  = rj1*a1 + rj2*a2
      
! offset in graphene coordinates to get atomic position (Angstroms)      
  CALL tauAjVecXY(1,1,tau)
  IF (iatom == 1) THEN
     Rj = Rj-tau/2.D0
  END IF
  IF(iatom == 2) THEN
     Rj = Rj + tau/2.D0
  END IF
      
! convert to xyz coordinates (Angstroms)      
  CALL XYtoTube(n,m,Rj,Rj)
  CALL tubetoXYZ(n,m,Rj,Rxyz)            
      
  RETURN

END SUBROUTINE rxyzJ1J2Vec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE NNj1j2(iatom,ivec,nn,j1,j2)
!===============================================================================
! for an atom in the two atom unit cell, find the two atom unit cell
! indices j1 and j2 that its neighbors belong to. The indices j1 and j2
! specify the center of the unit cell as R = j1*a1 + j2*a2.
!===============================================================================
! Input        :
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            neighbor index nn = 0,1,2,3,4
! Output       :
!  j1,j2         center of two atom unit cell in (a1,a2)
!===============================================================================
  IMPLICIT NONE
  
! input variables    
  INTEGER, INTENT(in)    :: iatom, ivec, nn

! output variables
  INTEGER, INTENT(out)   :: j1, j2
  
! working variable
  INTEGER                :: ierr
      
! check input for errors
  ierr = 0
      
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ierr = 1
     WRITE (*,*) 'NNj1j2 err: invalid iatom: ', iatom 
  END IF
      
  IF (nn < 0 .OR. nn > 4) THEN
     ierr = 1
     WRITE (*,*) 'NNj1j2 err: invalid nn:    ', nn
  END IF
      
  IF (nn == 1 .AND. ivec > 3) THEN
     ierr = 1
     WRITE (*,*) 'NNj1j2 err: nn = 1 and ivec > 3'
  END IF
      
  IF (nn == 2 .AND. ivec > 6) THEN
     ierr = 1
     WRITE (*,*) 'NNj1j2 err: nn = 2 and ivec > 6'
  END IF
      
  IF (nn == 3 .AND. ivec > 3) THEN
     ierr = 1
     WRITE (*,*) 'NNj1j2 err: nn = 3 and ivec > 3'
  END IF
      
  IF (nn == 4 .AND. ivec > 6) THEN
     ierr = 1
     WRITE (*,*) 'NNj1j2 err: nn = 4 and ivec > 6'
  END IF
      
  IF (ierr.NE.0) STOP
      
! find j1 and j2 for an A atom (iatom = 1)

! onsite case
  IF (nn == 0) THEN
     j1 = 0
     j2 = 0
  END IF

! nearest neighghbor J1 J2 indices
  IF (nn == 1) THEN
      
     IF (ivec == 1) THEN
        j1 = 0
        j2 = 0
     END IF

     IF (ivec == 2) THEN
        j1 = 0
        j2 =-1
     END IF
      
     IF (ivec == 3) THEN
        j1 =-1
        j2 = 0
     END IF

  END IF
      
! second nearest neighbor J1 J2 indices
  IF (nn == 2) THEN
     
     IF (ivec == 1) THEN
        j1 = 1
        j2 = 0
     END IF

     IF (ivec == 2) THEN        
        j1 = 0
        j2 = 1
     END IF
      
     IF (ivec == 3) THEN
        j1 = 1
        j2 =-1
     END IF
      
     IF (ivec == 4) THEN
        j1 =-1
        j2 = 1
     END IF

     IF (ivec == 5) THEN
        j1 = 0
        j2 =-1      
     END IF
     
     IF (ivec == 6) THEN
        j1 =-1
        j2 = 0       
     END IF

  END IF

! third nearest neighghbor J1 J2 indices
  IF (nn == 3) THEN
      
     IF (ivec == 1) THEN
        j1 = -1
        j2 = -1
     END IF

     IF (ivec == 2) THEN
        j1 = 1
        j2 =-1
     END IF
      
     IF (ivec == 3) THEN
        j1 = -1
        j2 = 1
     END IF
      
  END IF
      
! fourth nearest neighghbor J1 J2 indices
  IF (nn == 4) THEN
      
     IF (ivec == 1) THEN
        j1 = 1
        j2 = 0
     END IF

     IF (ivec == 2) THEN
        j1 = 0
        j2 = 1
     END IF
      
     IF (ivec == 3) THEN
        j1 = 1
        j2 =-2
     END IF
      
     IF (ivec == 4) THEN
        j1 =-2
        j2 = 1
     END IF

     IF (ivec == 5) THEN
        j1 = 0
        j2 =-2
     END IF
      
     IF (ivec == 6) THEN
        j1 =-2
        j2 = 0      
     END IF
      
  END IF

! change sign of j1 and j2 if this is a B atom (iatom=2)
  IF (iatom == 2) THEN
     j1 = -j1
     j2 = -j2
  END IF
      
END SUBROUTINE NNj1j2
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tubetoXYZ(n,m,Vtube,Vxyz)
!===============================================================================
! Convert 2D vector in tube coordinates to 3D xyz coordinates
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in a1,a2 basis
!  Vtube(2)      2D vector in tube coordinates (Angstroms)
! Output       :
!  Vxyz(3)       3D vector in xyz coordinates (Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  REAL(8), INTENT(in)    :: Vtube(2)

! output variable
  REAL(8), INTENT(out)   :: Vxyz(3)

! working variables      
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0
  REAL(8)                :: theta, tubeRadius, tubeDiam, chLength
      
  tubeRadius = tubeDiam(n,m)/2.D0
  theta      = 2.D0 * pi * Vtube(1)/chLength(n,m)
      
  Vxyz(1) = tubeRadius*DCOS(theta)
  Vxyz(2) = tubeRadius*DSIN(theta)
  Vxyz(3) = Vtube(2)
      
END SUBROUTINE tubetoXYZ
!*******************************************************************************
!*******************************************************************************
SUBROUTINE XYtoTube(n,m,Vin,Vout)
!===============================================================================
! Convert an xy vector to unrolled tube coordinates.
! In tube coordinates, x is taken parallel to the chiral vector (Ch)
! and y is parallel to the translation vector (T)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in a1,a2 basis
!  Vin(2)        vector in xy coordinates (A)
! Output       :
!  Vout(2)       vector in tube coordinates (A)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m
  REAL(8), INTENT(in)    :: Vin(2)

! output variable
  REAL(8), INTENT(out)   :: Vout(2)
      
! working variables
  REAL(8)                :: Ch(2), T(2)
  REAL(8)                :: Vtemp(2)
  
  INTEGER                :: it1,it2
  REAL(8)                :: rL,chLength
  REAL(8)                :: rT,trLength
      
  CALL chVecXY(n,m,Ch)      
  CALL trVecXY(n,m,it1,it2,T)
      
  rL = chLength(n,m)
  rT = trLength(n,m)            
      
  Vtemp(1) = DOT_PRODUCT(Vin,Ch)
  Vtemp(2) = DOT_PRODUCT(Vin, T)
  
  Vtemp(1) = Vtemp(1)/rL
  Vtemp(2) = Vtemp(2)/rT
      
  Vout(1)  = Vtemp(1)
  Vout(2)  = Vtemp(2)      

END SUBROUTINE XYtoTube
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tauAjVecXY(ivec,nn,tau)
!===============================================================================
! vectors from a graphene A atom to nearest neighbors up to the fourth
! nearest neighbor shell in graphene xy coordinates
!-------------------------------------------------------------------------------
! Input        :
!  ivec          index for nearest neighbor vector in the shell
!  nn            nearest neighbor index nn = 1,2,3,4
! Output       :
!  tau(2)        vector from A to nearest neighbor in xy coordinates (A)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: ivec, nn

! output variable
  REAL(8), INTENT(out)   :: tau(2)

! parameters
  REAL(8), PARAMETER     :: a     = 2.49D0
  REAL(8), PARAMETER     :: sqrt3 = 1.73205080756888D0

! check input for errors
  IF (nn <= 0 .OR. nn > 4) THEN
     WRITE (*,*) 'tauAjVecXY err: invalid nn: ', nn
     STOP
  END IF
      
  IF (nn == 1 .AND. ivec > 3) THEN
     WRITE (*,*) 'nn = 1 and ivec > 3'
     STOP
  END IF
      
  IF (nn == 2 .AND. ivec > 6) THEN
     WRITE (*,*) 'nn = 2 and ivec > 6'
     STOP
  END IF
      
  IF (nn == 3 .AND. ivec > 3) THEN
     WRITE (*,*) 'nn = 3 and ivec > 3'
     STOP
  END IF
      
  IF (nn == 4 .AND. ivec > 6) THEN
     WRITE (*,*) 'nn = 4 and ivec > 6'
     STOP
  END IF
      
! nearest neighghbor vectors (A)
  IF (nn == 1) THEN
      
     IF (ivec == 1) THEN
        tau(1) = a/sqrt3
        tau(2) = 0.D0
        RETURN
     END IF

     IF (ivec == 2) THEN
        tau(1) =-a/(2.D0*sqrt3)
        tau(2) = a/2.D0
        RETURN
     END IF
      
     IF (ivec == 3) THEN
        tau(1) =-a/(2.D0*sqrt3)
        tau(2) =-a/2.D0
        RETURN
     END IF

     WRITE (*,*) 'tauAjVecXY err: ivec,nn:', ivec, nn
     STOP
  END IF
      
! second nearest neighghbor vectors (A)
  IF (nn == 2) THEN
      
     IF (ivec == 1) THEN
        tau(1) = a*sqrt3/2.
        tau(2) = a/2.
        RETURN
     END IF

     IF (ivec == 2) THEN
        tau(1)= a*sqrt3/2.
        tau(2)=-a/2.        
        RETURN
     END IF
      
     IF (ivec == 3) THEN
        tau(1) = 0.D0
        tau(2) = a
        RETURN
     END IF
      
     IF (ivec == 4) THEN
        tau(1) = 0.D0
        tau(2) =-a
        RETURN
     END IF

     IF (ivec == 5) THEN
        tau(1) =-a*sqrt3/2.D0
        tau(2) = a/2.D0
        RETURN
     END IF
      
     IF (ivec == 6) THEN
        tau(1) = -a*sqrt3/2.D0
        tau(2) =-a/2.D0 
        RETURN
     END IF

     WRITE(*,*) 'tauAjVecXY err: ivec,nn:',ivec,nn      
     STOP
  END IF

! third nearest neighghbor vectors (A)
  IF (nn == 3) THEN
      
     IF(ivec == 1) THEN
        tau(1) = -2.D0*a/sqrt3
        tau(2) = 0.D0
        RETURN
     END IF

     IF (ivec == 2) THEN
        tau(1) = a/sqrt3
        tau(2) = a
        RETURN
     END IF
     
     IF (ivec == 3) THEN
        tau(1) = a/sqrt3
        tau(2) =-a
        RETURN
     END IF
     
     WRITE (*,*) 'tauAjVecXY err: ivec,nn:',ivec,nn
     STOP      
  END IF
      
! fourth nearest neighghbor vectors (A)
  IF (nn == 4) THEN
      
     IF (ivec == 1) THEN
        tau(1) = 5.*a/(2.D0*sqrt3)
        tau(2) = a/2.D0
        RETURN
     END IF

     IF (ivec == 2) THEN
        tau(1) = 5.*a/(2.*sqrt3)
        tau(2) =-a/2.        
        RETURN
     END IF
      
     IF (ivec == 3) THEN
        tau(1) =-a/(2.D0*sqrt3)
        tau(2) = 3.D0*a/2.D0
        RETURN
     END IF
      
     IF (ivec == 4) THEN
        tau(1) =-a/(2.D0*sqrt3)
        tau(2) =-3.D0*a/2.D0
        RETURN
     END IF

     IF (ivec == 5) THEN
        tau(1) =-2.D0*a/sqrt3
        tau(2) = a
        RETURN
     END IF
     
     IF (ivec == 6) THEN
        tau(1) = -2.D0*a/sqrt3
        tau(2)=-a        
        RETURN
     END IF

     WRITE (*,*) 'tauAjVecXY err: ivec,nn:', ivec, nn
     STOP
  END IF
          
  WRITE (*,*) 'tauAjVecXY err: never get here'
  WRITE (*,*) 'ivec,nn:', ivec, nn
  STOP
      
END SUBROUTINE tauAjVecXY
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tauBjVecXY(ivec,nn,tau)
!===============================================================================
! vectors from a graphene B atom to nearest neighbors up to the fourth
! nearest neighbor shell in graphene xy coordinates.
!-------------------------------------------------------------------------------
! Input        :
!  ivec          index for nearest neighbor vector in the shell
!  nn            nearest neighbor index nn = 1,2,3,4
! Output       :
!  tau(2)        vector from B to nearest neighbor in xy coordinates (A)
!===============================================================================
  IMPLICIT NONE
      
! input variables
  INTEGER, INTENT(in)    :: ivec, nn

! output
  REAL(8), INTENT(out)   :: tau(2)
      
  CALL tauAjVecXY(ivec,nn,tau)
  tau(1) = -tau(1)
  tau(2) = -tau(2)
      
END SUBROUTINE tauBjVecXY
!*******************************************************************************
!*******************************************************************************
SUBROUTINE rhoVec(n,m,iatom,rho)
!===============================================================================
! out-of-surface unit xyz vector defined in Jiang et al for the
! A or B atom in an (n,m) carbon nanotube.
! Reference: J.-W, Jiang et al, PRB, 73, 235434 (2006). 
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
! Output       :
!  rho(3)        xyz coordinates of out-of-surface unit vector (none) 
!=======================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)  :: n, m, iatom

! output variable
  REAL(8), INTENT(out) :: rho(3)
      
! working variables
  REAL(8)              :: r0(3), r1(3)     
  INTEGER              :: ier, nn, ivec
  REAL(8)              :: rnorm, vecLength
      
! check input for errors
  ier = 0
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ier = 1
     WRITE (*,*) 'rhoVec err: invalid iatom:', iatom
  END IF
  
  IF(ier /= 0) STOP
      
! z component of A or B atom position (Angstroms)
  nn = 0
  CALL rxyzVec(n,m,iatom,1,nn,r0)
  r0(1) = 0.D0
  r0(2) = 0.D0
      
! sum over ivec = 1..3
  nn  = 1
  rho = 0.D0
  DO ivec = 1, 3
     CALL rxyzVec(n,m,iatom,ivec,nn,r1)
     rho = rho + r1-r0      
  END DO
      
! unit rho vector (dimensionless)                 
  rnorm = vecLength(3,rho)
  rho = rho/rnorm 

END SUBROUTINE rhoVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE rho2Vec(n,m,iatom,ivec,rho)
!===============================================================================
! out-of-surface unit xyz vector defined in Jiang et al for the three
! nearest neighbors of the A or B atom in an (n,m) carbon nanotube.
! Reference: J.-W, Jiang et al, PRB, 73, 235434 (2006).
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in a1,a2 basis
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for three nearest neighbors (ivec=1..3)
! Output       :
!  rho(3)        xyz coordinates of out-of-surface unit vector (none)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec

! output variable
  REAL(8), INTENT(out)   :: rho(3)
      
! working variables
  REAL(8)                :: r0(3), r1(3)
  REAL(8)                :: rnorm, vecLength ! function vecLength
  INTEGER                :: ier, nn, ivec1

! check input for errors
  ier = 0
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ier = 1
     WRITE (*,*) 'rho2Vec err: invalid iatom:', iatom
  END IF
  IF (ivec < 1 .OR. ivec > 3) THEN
     ier = 1
     WRITE (*,*) 'rho2Vec err: invalid ivec: ', ivec
  END IF
  IF (ier /= 0) STOP
      
! z component of R(1)_(iatom,ivec) position vector (Angstroms)
  nn = 1
  CALL rxyzVec(n,m,iatom,ivec,nn,r0)
  r0(1) = 0.D0
  r0(2) = 0.D0
      
! sum over ivec = 1..3
  rho = 0.D0
  DO ivec1 = 1, 3
     CALL r2xyzVec(n,m,iatom,ivec,ivec1,r1)        
     rho = rho + (r1 - r0)
  END DO
      
! unit rho vector (dimensionless)                 
  rnorm = vecLength(3,rho)
  rho   = rho/rnorm                         
      
END SUBROUTINE rho2Vec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tubeRadUnitVec(Vxyz,Ur)
!===============================================================================
! unit outward radial vector in a carbon nanotube passing through Vxyz
!-------------------------------------------------------------------------------
! Input        :
!  Vxyz(3)       position vector in xyz coordinates (Angstroms)
! Output       :
!  Ur(3)         radial unit vector passing through Vxyz (dimensionless)
!===============================================================================
  IMPLICIT NONE

! input variable
  REAL(8), INTENT(in)    :: Vxyz(3)

! output variable
  REAL(8), INTENT(out)   :: Ur(3)

! working variable
  REAL(8)                :: Vlength, pythagoras ! from libMath.f90
      
  Vlength = pythagoras(Vxyz(1),Vxyz(2))
  
  Ur(1) = Vxyz(1)/Vlength
  Ur(2) = Vxyz(2)/Vlength
  Ur(3) = 0.D0
      
END SUBROUTINE tubeRadUnitVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE r2xyzVec(n,m,iatom,ivec1,ivec2,Rxyz)
!===============================================================================
! Compute the xyz coordinates for R(1)_(iatom,ivec1,ivec2) which is the
! ivec2 nearest neighbor of the ivec1 nearest neighbor of atom A or B
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec1         index for nearest neighbor of R(1)_(iatom)
!  ivec2         index for nearest neighbor of R(1)_(iatom,ivec1)
! Output       :
!  Rxyz(3)       R(1)_(iatom,ivec1,ivec2) in xyz coordinates (Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec1, ivec2

! output variable
  REAL(8), INTENT(out)   :: Rxyz(3)
      
! working variables
  REAL(8)                :: tau0(2), tau1(2), tau2(2), Rtube(2)      

  INTEGER                :: ier,nn

! check input for errors
  ier = 0
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ier = 1
     WRITE (*,*) 'r2xyzVec err: invalid iatom: ', iatom
  END IF
  
  IF (ivec1 < 1 .OR. ivec1 > 3) THEN
     ier = 1
     WRITE (*,*) 'r2xyzVec err: invalid ivec1: ', ivec1
  END IF
  
  IF (ivec2 < 1 .OR. ivec2 > 3) THEN
     ier = 1
     WRITE (*,*) 'r2xyzVec err: invalid ivec2: ', ivec2
  END IF
  
  IF (ier /= 0) STOP 'r2xyzVec unknown err'
      
! tube coordinates of A or B atom in two atom unit cell (Angstroms)
  CALL tauAjVecXY(1,1,tau0)
  IF (iatom == 1) THEN
     tau0 =-tau0/2.D0 !tau_A = -tau(1)_(A1)/2 (graphene xy coordinates)
  ELSE
     tau0 = tau0/2.D0 !tau_B =  tau(1)_(A1)/2
  END IF
  CALL XYtoTube(n,m,tau0,tau0)
      
! tube coordinates of pointer to neighbor of R_(iatom) (A)
  nn = 1
  IF(iatom == 1) THEN
     CALL tauAjVecXY(ivec1,nn,tau1)
  ELSE
     CALL tauBjVecXY(ivec1,nn,tau1)      
  END IF
  CALL XYtoTube(n,m,tau1,tau1)
      
! tube coordinates of pointer to neighbor of R_(iatomiivec1) (A)
  nn = 1
  IF (iatom == 2) THEN
     CALL tauAjVecXY(ivec2,nn,tau2)
  ELSE
     CALL tauBjVecXY(ivec2,nn,tau2)      
  END IF
  CALL XYtoTube(n,m,tau2,tau2)
      
! go from tube to xyz coordinates (Angstroms)
  Rtube = tau0 + tau1 + tau2        
  CALL tubetoXYZ(n,m,Rtube,Rxyz)

END SUBROUTINE r2xyzVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE drVec(n,m,iatom,Dr)
!===============================================================================
! compute the xyz coordinates of the vector
! Dr = sum(ivec=1..3, Nhat(1)_(iatom,ivec))
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
! Output:
!  Dr(3)         sum(ivec=1..3, Nhat(1)_(iatom,ivec))
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom
  REAL(8), INTENT(out)   :: Dr(3)
  
! working variables
  REAL(8)                :: Dr1(3), Dr2(3), Dr3(3)
  INTEGER                :: nn,ivec

  nn   = 1
  ivec = 1
  CALL dNxyzVec(n,m,iatom,ivec,nn,Dr1)
  
  ivec = 2
  CALL dNxyzVec(n,m,iatom,ivec,nn,Dr2)
  
  ivec = 3
  CALL dNxyzVec(n,m,iatom,ivec,nn,Dr3)
  
  Dr = Dr1 + Dr2 + Dr3
         
END SUBROUTINE drVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE dr2Vec(n,m,iatom,ivec,Dr)
!===============================================================================
! compute the xyz coordinates of the vector
! Dr = sum( ivec2=1..3, Nhat(1)_((iatom,ivec),ivec2) )
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in a1,a2 basis
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          near neighbor index (1,2,3)
! Output       :
!  Dr(3)         sum( ivec2=1..3, Nhat(1)_((iatom,ivec),ivec2) )
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec

! output variable
  REAL(8), INTENT(out)   :: Dr(3)

! working variables      
  INTEGER                :: ivec2
  REAL(8)                :: Dr1(3), Dr2(3), Dr3(3)
  
  ivec2 = 1
  CALL dN2xyzVec(n,m,iatom,ivec,ivec2,Dr1)      
  
  ivec2 = 2
  CALL dN2xyzVec(n,m,iatom,ivec,ivec2,Dr2)      
  
  ivec2 = 3
  CALL dN2xyzVec(n,m,iatom,ivec,ivec2,Dr3)                  
      
  Dr = Dr1 + Dr2 + Dr3
  
END SUBROUTINE dr2Vec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE dN2xyzVec(n,m,iatom,ivec1,ivec2,dN2xyz)
!===============================================================================
! compute the xyz coordinates for the dimensionless unit vector
! Nhat(1)_(iatom,ivec1,ivec2) which points from R(1)_(iatom,ivec1) to
! R(1)_((iatom,ivec1),ivec2)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec1,ivec2   near neighbor index labels (1,2,3)
! Output       :
!  dNxyz(3)      xyz unit vector from R(1)_(iatom,ivec1) to
!                R(1)_((iatom,ivec1),ivec2)
!=======================================================================      
  IMPLICIT NONE

! input variables
  INTEGER,  INTENT(in)   :: n,m,iatom,ivec1,ivec2

! output variable
  REAL(8),  INTENT(out)  :: dN2xyz(3)

! working variables      
  REAL(8)                :: R0(3), Rxyz(3)     
  INTEGER                :: ier, m1, m2
  REAL(8)                :: rnorm, vecLength

! check for errors      
  ier = 0
  IF (ivec1 < 1 .OR. ivec1 > 3) THEN
     WRITE (*,*) 'dN2xyzVec err: invalid ivec1:', ivec1
     ier = 1
  END IF
  
  IF (ivec2 < 1 .OR. ivec2 > 3) THEN
     WRITE (*,*) 'dN2xyzVec err: invalid ivec2:',ivec2
     ier = 1
  END IF
  
  IF(ier /= 0) STOP      
      
! start definition
  IF (ivec1 == 1 .AND. ivec2 == 1) THEN
     CALL dNxyzVec(n,m,iatom,ivec1,1,dN2xyz)
     dN2xyz = -dN2xyz
  END IF
      
  IF (ivec1 == 1 .AND. ivec2 == 2) THEN
     m2 = 2
     m1 = 1
     CALL rxyzVec(n,m,iatom,m2,2,Rxyz)      
     CALL rxyzVec(n,m,iatom,m1,1,R0)
     dN2xyz = Rxyz - R0
     rnorm  = vecLength(3,dN2xyz)
     dN2xyz = dN2xyz/rnorm 
  END IF
  
  IF (ivec1 == 1 .AND. ivec2 == 3) THEN
     m2 = 1
     m1 = 1
     CALL rxyzVec(n,m,iatom,m2,2,Rxyz)      
     CALL rxyzVec(n,m,iatom,m1,1,R0)
     dN2xyz = Rxyz-R0
     rnorm  = vecLength(3,dN2xyz)
     dN2xyz = dN2xyz/rnorm      
  END IF
           
  IF (ivec1 == 2 .AND. ivec2 == 1) THEN
     m2 = 5
     m1 = 2
     CALL rxyzVec(n,m,iatom,m2,2,Rxyz)      
     CALL rxyzVec(n,m,iatom,m1,1,R0)
     dN2xyz = Rxyz - R0
     rnorm  = vecLength(3,dN2xyz)
     dN2xyz = dN2xyz/rnorm      
  END IF
      
  IF (ivec1 == 2 .AND. ivec2 == 2) THEN
     CALL dNxyzVec(n,m,iatom,ivec1,1,dN2xyz)
     dN2xyz = -dN2xyz
  END IF
      
  IF (ivec1 == 2 .AND. ivec2 == 3) THEN
     m2 = 3
     m1 = 2
     CALL rxyzVec(n,m,iatom,m2,2,Rxyz)      
     CALL rxyzVec(n,m,iatom,m1,1,R0)
     dN2xyz = Rxyz - R0
     rnorm  = vecLength(3,dN2xyz)
     dN2xyz = dN2xyz/rnorm      
  END IF
  
  IF (ivec1 == 3 .AND. ivec2 == 1) THEN
     m2 = 6
     m1 = 3
     CALL rxyzVec(n,m,iatom,m2,2,Rxyz)      
     CALL rxyzVec(n,m,iatom,m1,1,R0)
     dN2xyz = Rxyz-R0
     rnorm  = vecLength(3,dN2xyz)
     dN2xyz = dN2xyz/rnorm      
  END IF
      
  IF (ivec1 == 3 .AND. ivec2 == 2) THEN
     m2 = 4
     m1 = 3
     CALL rxyzVec(n,m,iatom,m2,2,Rxyz)      
     CALL rxyzVec(n,m,iatom,m1,1,R0)
     dN2xyz = Rxyz - R0
     rnorm  = vecLength(3,dN2xyz)
     dN2xyz = dN2xyz/rnorm      
  END IF
      
  IF (ivec1 ==  3 .AND. ivec2 == 3) THEN
     CALL dNxyzVec(n,m,iatom,ivec1,1,dN2xyz)
     dN2xyz = -dN2xyz
  END IF

  RETURN
END SUBROUTINE dN2xyzVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE dNxyzVec(n,m,iatom,ivec,nn,dNxyz)
!===============================================================================
! compute the xyz coordinates for the dimensionless unit vector
! Nhat(nn)_(iatom,ivec) which points from R_iatom to R(nn)_(iatom,ivec).
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            nearest neighbor index nn = 1,2,3,4
! Output:
!  dNxyz(3)      xyz unit vector from R_iatom to R(nn)_(iatom,ivec)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variable
  REAL(8), INTENT(out)   :: dNxyz(3)
  
! working variables
  REAL(8)                :: R0(3), Rxyz(3)     
  INTEGER                :: i
  REAL(8)                :: rnorm, vecLength
      
  CALL rxyzVec(n,m,iatom,ivec,nn,Rxyz)      
  CALL rxyzVec(n,m,iatom,ivec,0,R0)
  
  dNxyz = Rxyz - R0
  rnorm = vecLength(3,dNxyz)
  DO i = 1, 3
     dNxyz(i) = dNxyz(i)/rnorm
  END DO
                        
END SUBROUTINE dNxyzVec
!*******************************************************************************
!*******************************************************************************
SUBROUTINE srmatrix(n,m,iatom,Sr)
!===============================================================================
! compute the sum of outer product matrices given by
! Sr = sum(ivec=1,3  Nhat(1)_(iatom,ivec) x Nhat(1)_(iatom,ivec) )
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
! Output       :
!  Sr(3,3)       sum(ivec=1,3  Hhat(1)_(iatom,ivec) x Hhat(1)_(iatom,ivec) ) 
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom

! output variable
  REAL(8), INTENT(out)   :: Sr(3,3)
  
! working variables
  REAL(8)                :: dNxyz(3), U(3,3)
  INTEGER                :: nn, ivec
      
  Sr = 0.D0     
  nn = 1      
  DO ivec = 1,3
     CALL dNxyzVec(n,m,iatom,ivec,nn,dNxyz)
     CALL outProd(3,dNxyz,3,dNxyz,3,U)    
     Sr = Sr + U ! sum up
  END DO
  
END SUBROUTINE srmatrix
!*******************************************************************************
!*******************************************************************************
SUBROUTINE sr2matrix(n,m,iatom,ivec,Sr)
!===============================================================================
! compute the sum of outer product matrices given by Sr = sum(ivec2=1,3
! Nhat(1)_((iatom,ivec),ivec2)) x Nhat(1)_((iatom,ivec),ivec2) )
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
! Output:
!  Sr(3,3)       sum_i Nhat(1)_(iatom,ivec),i)) x Nhat(1)_(iatom,ivec),i) ) 
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec

! output variable
  REAL(8), INTENT(out)   :: Sr(3,3)
  
! working variables
  REAL(8)                :: dNxyz(3), U(3,3)
  INTEGER                :: ivec2
      
  Sr = 0.D0
            
  DO ivec2 = 1, 3
     CALL dN2xyzVec(n,m, iatom,ivec,ivec2, dNxyz)
     CALL outProd(3,dNxyz,3,dNxyz,3,U)
! sum up
     Sr = Sr+U
  END DO
      
END SUBROUTINE sr2matrix
!*******************************************************************************
!*******************************************************************************
SUBROUTINE sqmatrix(n,m,iatom,ivec,nn,q,mu,Sq)
!===============================================================================
! compute the matrix S(j)*exp(i*phi(j)) (dimensionless)
!------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            nearest neighbor index nn = 1,2,3,4
!  q             wavevector along tube axis (1/Angstroms)
!  mu            labels phonon manifolds (0...N_hex-1)
! Output       :
!  Sq(3,3)       complex S(j)*exp(i*phi(j))   
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn, mu
  REAL(8), INTENT(in)    :: q

! output variable
  COMPLEX(8),INTENT(out) :: Sq(3,3)

! working variables              
  REAL(8)                :: Sj(3,3)
  INTEGER                :: ier, i, j
  REAL(8)                :: phi     
  COMPLEX(8)             :: ci, cphase

  ci = (0.D0,1.D0)
      
! check input for errors
  ier = 0
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ier = 1
     WRITE (*,*) 'Sqmatrix err: invalid iatom: ', iatom
  END IF
  IF (ier.NE.0) STOP
      
! rotation matrix Sj (dimensionless)
  CALL sjmatrix(n,m,iatom,ivec,nn,Sj)

! phase angle phi
  CALL phij(n,m,iatom,ivec,nn,q,mu,phi)      
      
! multiply Sj by phase factor exp(I*phi)
  cphase = CDEXP(ci*phi)
  DO i = 1, 3
     DO j = 1, 3
        Sq(i,j) = Sj(i,j) * cphase
     END DO
  END DO
      
  RETURN
END SUBROUTINE sqmatrix
!*******************************************************************************
!*******************************************************************************
SUBROUTINE sjmatrix(n,m,iatom,ivec,nn,Sj)
!===============================================================================
! the rotation matrix about the tube axis after the application of
! two screw operations a1(j1) and a2(j2).
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            neighbor index nn = 1,2,3,4
! Output       :
!  Sj(3,3)       rotation matrix (dimensionless)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn
  REAL(8), INTENT(out)   :: Sj(3,3)
  
! working variables
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0
  INTEGER                :: j1, j2
  REAL(8)                :: rn, rm, denom, theta1, theta2, theta, cost, sint
      
! rotation angles about tube axis for a1 and a2 screw operators
  rn = float(n)
  rm = float(m)
  denom  = rn**2 + rn*rm + rm**2
  theta1 = (2.D0*rn + rm)*pi / denom
  theta2 = (2.D0*rm + rn)*pi / denom      
      
! total rotation angle about tube axis
  CALL NNj1j2(iatom,ivec,nn,j1,j2)
  theta = DBLE(j1)*theta1 + DBLE(j2)*theta2
      
! rotation matrix      
  cost = DCOS(theta)
  sint = DSIN(theta)
      
  Sj(1,1) = cost
  Sj(1,2) =-sint
  Sj(1,3) = 0.D0
      
  Sj(2,1) = sint
  Sj(2,2) = cost
  Sj(2,3) = 0.D0
      
  Sj(3,1) = 0.D0
  Sj(3,2) = 0.D0
  Sj(3,3) = 1.D0        
      
END SUBROUTINE sjmatrix
!*******************************************************************************
!*******************************************************************************
SUBROUTINE reducedCutLine(n,m,muin,muout)
!===============================================================================
! for (n,m) tube put cutting line index mu in range mu = 1 ... nhex
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in a1,a2 basis 
!  mu_in         cutting line index (unreduced)
! Output       :
!  mu_out        cutting line index (reduced)
!===============================================================================
  IMPLICIT NONE      

! input variables
  INTEGER, INTENT(in)    :: n, m, muin

! output variable
  INTEGER, INTENT(out)   :: muout

! working variables
  INTEGER                :: nhex, nHexagon
  
  nhex  = nHexagon(n,m)
  muout =MOD(ABS(muin),nhex)
      
  IF (muin <= 0 ) muout = nhex - muout
  IF (muout == 0) muout =nhex
      
END SUBROUTINE reducedCutLine
!*******************************************************************************
!*******************************************************************************
!-------------------------------------------------------------------------------
!--------------- added by Daria Satco ------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE getHexagonPosition(n,m,nhex,j1j2)
!==============================================================================
! calculates J=(j1,j2) which correspond to all hexagons within CNT cell
! Input
! n,m           chiral vector coordinates in a1,a2 basis
! Output
! j1j2          array of (j1,j2)   (nhex)
!==============================================================================
IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, nhex

! output variable
  INTEGER, INTENT(out)   :: j1j2(nhex,2)

! working variables
  INTEGER                :: it1, it2
  REAL(8)                :: T(2)
  REAL(8)                :: Ch(2)
  REAL(8)                :: a1(2), a2(2)
  INTEGER                :: j1,j2, nn
  REAL(8)                :: cond1, cond2, cond3, cond4

  CALL  trVecXY(n,m,it1,it2,T)
  CALL  chVecXY(n,m,Ch)
  CALL  unitVecXY(a1,a2)

  nn = 0

  DO j1 = 0, it1+n
    DO j2 = it2, m

        cond1 = real(it2*j1)/it1
        cond2 = real(m*j1)/n
        cond3 = real(it2*(j1-n))/it1
        cond4 = real(m*(j1-it1))/n

        IF ( (cond1 .le. REAL(j2)) .and. (cond2 .ge. REAL(j2)) .and. (cond3 > REAL(j2-m)) .and. (cond4 < REAL(j2-it2)) ) THEN
            nn = nn + 1
            j1j2(nn,1) = j1
            j1j2(nn,2) = j2

!            PRINT*, j1j2(nn,1), j1j2(nn,2)

        END IF

    END DO
  END DO

IF ( nhex .ne. nn ) STOP 'wrong j1j2 array'

END SUBROUTINE getHexagonPosition
!*****************************************************************************
!*****************************************************************************
