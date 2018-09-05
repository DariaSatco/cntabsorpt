!*******************************************************************************
!*******************************************************************************
! Project      : libswntPhon.f90
!===============================================================================
! Purpose      :
! Calculate vibrational states for a SWNT using force-constant model
!-------------------------------------------------------------------------------
! Authors      : ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
!                Gary Sanders (sanders@phys.ufl.edu) 
! Latest Vers. : 2012.11.09
!-------------------------------------------------------------------------------
! Reference(s) :
! [1] Physical Properties of Carbon Nanotubes
!     R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
! [2] J.-W, Jiang et al, PRB, 73, 235434 (2006)
!-------------------------------------------------------------------------------
! Contents     :
! - SUBROUTINE fcTubeEm(n,m,mu,q,homega)
! - SUBROUTINE tubePhDOS(n,m,ne,Earray, DOS)
! - SUBROUTINE fcTubeDm(n,m,mu,q,Ds)
! - SUBROUTINE grphpar(alpha1,alpha2,alpha3,alpha4,beta,gamma1,gamma2)
! - SUBROUTINE springF(n,m,iatom,ivec,nn,F)
! - SUBROUTINE bondBendingF(n,m,iatom,ivec,nn,F)
! - SUBROUTINE danglingBondF(n,m,iatom,ivec,nn,F)
! - SUBROUTINE bondTwistF(n,m,iatom,ivec,nn,F)
! - SUBROUTINE bondOutNormal(n,m,iatom,ibond,Ur)
!*******************************************************************************
!*******************************************************************************
SUBROUTINE fcTubeEm(n,m,mu,q,homega)
!===============================================================================
! compute the phonon energies hw(q,mu) in an (n,m) carbon nanotube
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  mu            labels phonon manifolds (1...N_hex)
!  q             wave vector along tube axis (1/A) (|q| < pi/T)
! Output       :
!  homega(6)     phonon energies (eV)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, mu
  REAL(8), INTENT(in)    :: q

! output variable
  REAL(8), INTENT(out)   :: homega(6) !(6)
  
! working variables
  COMPLEX(8)             :: Dq(6,6), S(6,6) 
  COMPLEX(8)             :: dDq(6,6), dZq(6,6) ! eigenvectors
  REAL(8)                :: w(6)               ! eigenvalues 
  
! lapack driver variables
  INTEGER                :: matz, nmat, il, iu, nout, i
      
! compute 6x6 dynamical matrix (Joules/meter**2)
  CALL fcTubeDm(n,m,mu,q,Dq)
  dDq=Dq

  nmat = 6
  S = 0.D0
  DO i = 1, nmat
     S(i,i) = 1.D0
  END DO

! option flag (0=evalues, 1=evalues+evectors)
  matz = 0 ! calculate evalues only

! if il or iu <= 0, all eigenvalues are returned
  il=0
  iu=0

! solveHam(n,ldh,ham,ldo,ovlp,matz,il,iu,nout,w,ldz,z)
! nout     ->  number of eigenvalues returned
! w(n)     ->  eigenvalues in ascending order
! ldz      ->  leading dimension of z
! z(ldz,n) ->  complex eigenvectors if matz = 1
  CALL solveHam(nmat,nmat,dDq,nmat,S,matz,il,iu,nout,w,nmat,dZq)
       
! return nanotube phonon energies in eV
  DO i = 1, 6
     homega(i) = 4.64327D-3 * DSQRT(DABS(w(i)))
     IF(w(i) < 0.D0) homega(i) = -homega(i)
  END DO

! uncomment below if want to return nanotube phonon dispersion in cm^{^1}
!!$  DO i = 1, 6
!!$     homega(i) = 8065.5447732 * homega(i)
!!$     IF(w(i) < 0.D0) homega(i) = -homega(i)
!!$  END DO
      
END SUBROUTINE fcTubeEm
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tubePhDOS(n,m,ne,Earray, DOS)
!===============================================================================
! Phonon density of states per carbon atom for an (n,m) carbon nanotube
! (states/carbon atom/eV)
!
! Note: Since there are three phonon degrees of freedom per carbon atom
! site, the total number modes per atom = 3. We thus have the sum rule:
!
!            Int( DOS(E), E=-infinity..infinity) = 3
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in a1,a2 basis
!  ne            number of energies
!  Earray(ne)    array of energies (eV)
! Output       :
!  DOS(ne)       zone folded phonon density of states
!===============================================================================
  IMPLICIT NONE

! parameters
  INTEGER, PARAMETER     :: nk = 301
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0     

! input variables
  INTEGER, INTENT(in)    :: n, m, ne
  REAL(8), INTENT(in)    :: Earray(ne)

! output variable
  REAL(8), INTENT(out)   :: DOS(ne)

! working variables
  REAL(8)                :: rka(nk), homega(6)
  REAL(8), ALLOCATABLE   :: Enk(:,:) !(nk,6*nhex)

  INTEGER                :: nout, nhex, nHexagon
  INTEGER                :: mu, k, indexx, ii, ie, nn

  REAL(8)                :: T, trLength
  REAL(8)                :: rk, rkmin, rkmax, dk, fwhm, fwhm1, E, DSn
      
! allocate storage
  nhex = nHexagon(n,m)
  nout = 6*nhex      
  ALLOCATE(Enk(nk,nout))      
      
! evenly spaced k points from 0 to pi/T (1/Angstroms)
  T     = trLength(n,m)
  rkmin = 0.D0
  rkmax = pi/T
  dk    = (rkmax-rkmin) / (DBLE(nk)-1.D0)
  DO k = 1, nk
     rka(k) = rkmin + (k-1)*dk
  END DO
      
! phonon energy bands at evenly spaced k points (eV)
  DO k = 1, nk
     rk = rka(k)
     indexx=0
     DO mu = 1, nhex
        CALL fcTubeEm(n,m,mu,rk,homega)
        DO ii = 1, 6
           indexx = indexx + 1
           Enk(k,indexx) = homega(ii)
        END DO
     END DO
  END DO
      
! find FWHM linewidth based on energy array (eV)
  fwhm = ABS(Earray(2) - Earray(1))
  DO ie = 1, ne-1
     fwhm1 = ABS(Earray(ie) - Earray(ie+1))
     IF(fwhm1 < fwhm) fwhm = fwhm1
  END DO
  fwhm = 5.D0*fwhm
      
! accumulate phonon density of states per unit length
  DO ie = 1, ne
     DOS(ie)=0.D0
     DO nn = 1, nout      
        E = Earray(ie)
        CALL dos1Dgauss(nout,nk,rka,nk,Enk,E,fwhm,nn,DSn)
        DOS(ie) = DOS(ie) + DSn
     END DO
  END DO
      
! convert to density of modes/Carbon_atom/eV
  DO ie = 1, ne
     DOS(ie) = 3.D0*(T/DBLE(nout)) * DOS(ie)
  END DO
      
  DEALLOCATE(Enk)      

END SUBROUTINE tubePhDOS
!*******************************************************************************
!*******************************************************************************
SUBROUTINE fcTubeDm(n,m,mu,q,Ds)
!===============================================================================
! Calculate nanotube symmetry adapted dynamical matrix (Joules/meter**2)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  mu            labels phonon manifolds (0...N_hex-1)
!  q             wave vector along tube axis (1/A) (- pi/T < q < pi/T)
! Output       :
!  Ds(6,6)       Hermitian dynamical matrix D(q,mu) (Joules/meter**2)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, mu
  REAL(8), INTENT(in)    :: q

! output variables
  COMPLEX(8),INTENT(out) :: Ds(6,6)

! working variables
  REAL(8)                :: F1(3,3), F2(3,3), F3(3,3), F4(3,3), F(3,3)
  REAL(8)                :: Diag(3,3)
  COMPLEX(8)             :: Daa(3,3),Dbb(3,3),Dab(3,3),Dba(3,3)
  COMPLEX(8)             :: Sq(3,3), Dtemp(6,6)    
  INTEGER                :: iatom, nn, ivec, i, j
  INTEGER, SAVE          :: nvecs(4) 

  nvecs = (/ 3, 6, 3, 6 /) 
  Ds = 0.D0

! Daa contribution
  Daa   = 0.D0
  iatom = 1

  Diag  = 0.D0
  DO nn = 1, 4
     DO ivec = 1, nvecs(nn)
        CALL springF(n,m,iatom,ivec,nn,F1)
        CALL bondBendingF(n,m,iatom,ivec,nn,F2)
        CALL danglingBondF(n,m,iatom,ivec,nn,F3)
        CALL bondTwistF(n,m,iatom,ivec,nn,F4)
        F = F1 + F2 + F3 + F4
        Diag = Diag + F
     END DO
  END DO
      
  nn = 2
  DO ivec = 1, nvecs(nn)
     CALL springF(n,m,iatom,ivec,nn,F1)
     CALL bondBendingF(n,m,iatom,ivec,nn,F2)
     CALL danglingBondF(n,m,iatom,ivec,nn,F3)
     CALL bondTwistF(n,m, iatom,ivec,nn,F4)
     F = F1 + F2 + F3 + F4
     CALL sqmatrix(n,m,iatom,ivec,nn,q,mu,Sq)
     Daa = Daa - MATMUL(F,Sq)
  END DO
  
  DO i = 1, 3
     DO j = 1, 3
        Ds(i,j) = Diag(i,j) + Daa(i,j)
     END DO
  END DO
      
! Dbb contribution
  Dbb   = 0.D0
  iatom = 2
      
  Diag  = 0.D0
  DO nn = 1, 4
     DO ivec = 1, nvecs(nn)
        CALL springF(n,m,iatom,ivec,nn,F1)
        CALL bondBendingF(n,m,iatom,ivec,nn,F2)
        CALL danglingBondF(n,m,iatom,ivec,nn,F3)
        CALL bondTwistF(n,m,iatom,ivec,nn,F4)
        F = F1 + F2 + F3 + F4
        Diag = Diag + F
     END DO
  END DO
      
  nn = 2
  DO ivec = 1, nvecs(nn)
     CALL springF(n,m,iatom,ivec,nn,F1)
     CALL bondBendingF(n,m,iatom,ivec,nn,F2)
     CALL danglingBondF(n,m,iatom,ivec,nn,F3)
     CALL bondTwistF(n,m,iatom,ivec,nn,F4)
     F = F1 + F2 + F3 + F4
     CALL sqmatrix(n,m,iatom,ivec,nn,q,mu,Sq)
     Dbb = Dbb - MATMUL(F,Sq)
  END DO
      
  DO i = 1, 3
     DO j = 1, 3
        Ds(i+3,j+3) = Diag(i,j) + Dbb(i,j)
     END DO
  END DO
      
! Dab contribution
  Dab   = 0.D0
  iatom = 1
  
  DO nn = 1, 4
     IF (nn /= 2) THEN
        DO ivec = 1, nvecs(nn)
           CALL springF(n,m,iatom,ivec,nn,F1)
           CALL bondBendingF(n,m,iatom,ivec,nn,F2)
           CALL danglingBondF(n,m,iatom,ivec,nn,F3)
           CALL bondTwistF(n,m,iatom,ivec,nn,F4)
           F = F1 + F2 + F3 + F4
           CALL sqmatrix(n,m,iatom,ivec,nn,q,mu,Sq)
           Dab = Dab-MATMUL(F,Sq)
        END DO
     END IF
  END DO
  
  DO i = 1, 3
     DO j = 1, 3
        Ds(i,j+3) = Dab(i,j)
     END DO
  END DO
      
! Dba contribution
  Dba   = 0.D0
  iatom = 2
  
  DO nn = 1, 4
     IF (nn /= 2) THEN
        DO ivec = 1, nvecs(nn)
           CALL springF(n,m,iatom,ivec,nn,F1)
           CALL bondBendingF(n,m,iatom,ivec,nn,F2)
           CALL danglingBondF(n,m,iatom,ivec,nn,F3)
           CALL bondTwistF(n,m,iatom,ivec,nn,F4)
           F = F1 + F2 + F3 + F4
           CALL sqmatrix(n,m, iatom,ivec,nn,q,mu,Sq)
           Dba = Dba - MATMUL(F,Sq)
        END DO
     END IF
  END DO
      
  DO i = 1, 3
     DO j = 1, 3
        Ds(i+3,j) = Dba(i,j)
     END DO
  END DO
      
! correct roundoff errors
  Dtemp = Ds
  DO i = 1, 6
     DO j = 1, 6
        Ds(i,j)= ( Dtemp(i,j) + CONJG(Dtemp(j,i)) ) / 2.D0
     END DO
  END DO
            
END SUBROUTINE fcTubeDm
!*******************************************************************************
!*******************************************************************************
SUBROUTINE grphpar(alpha1,alpha2,alpha3,alpha4,beta,gamma1,gamma2)
!===============================================================================
! Give graphene phonon parameters
!-------------------------------------------------------------------------------
! Output only  :
!  alphai        FC bond stretching parameters for shells 1..4 (J/m^2)
!  beta          in-plane bond bending
!  gamma1        dangling bond interaction
!  gamma2        bond-twisting interaction
!-------------------------------------------------------------------------------
  IMPLICIT NONE

! output variables      
  REAL(8), INTENT(out) :: alpha1, alpha2, alpha3, alpha4
  REAL(8), INTENT(out) :: beta, gamma1, gamma2
            
! bond stretching parameters
  alpha1 = 237.625D0
  alpha2 = 17.4688D0
  alpha3 = 12.8952D0
  alpha4 = 7.16627D0
  
! in-plane bond bending
  beta   = 18.3652D0

! dangling bond interaction
  gamma1 = 16.6685D0

! bond-twisting interaction
  gamma2 = 6.60915D0

  RETURN
  
END SUBROUTINE grphpar
!*******************************************************************************
!*******************************************************************************
SUBROUTINE springF(n,m,iatom,ivec,nn,F)
!===============================================================================
! Force constant matrices for spring interactions 
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          neighbor index
!  nn            near neighbor shell
! Output       :
!  F(3,3)        force constant matrix (Joules/meter**2)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variable
  REAL(8), INTENT(out)   :: F(3,3)
      
! working variables (force constant)
  REAL(8)                :: alpha(4), dNxyz(3), U(3,3)     
  REAL(8)                :: alpha1, alpha2, alpha3, alpha4
  REAL(8)                :: beta, gamma1, gamma2      
      
  CALL grphpar(alpha1,alpha2,alpha3,alpha4,beta,gamma1,gamma2)
     
! bond stretching parameters
  alpha(1) = alpha1
  alpha(2) = alpha2
  alpha(3) = alpha3
  alpha(4) = alpha4
  
  CALL dNxyzVec(n,m,iatom,ivec,nn,dNxyz)
  CALL outProd(3,dNxyz, 3,dNxyz, 3,U)
      
! fill in the matrix
  F = 2.D0 * alpha(nn) * U
      
END SUBROUTINE springF
!*******************************************************************************
!*******************************************************************************
SUBROUTINE bondBendingF(n,m,iatom,ivec,nn,F)
!===============================================================================
! Force constant matrices for bond bending interactions 
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          neighbor index
!  nn            near neighbor shell
! Output       :
!  F(3,3)        force constant matrix (Joules/meter**2)
!===============================================================================
  IMPLICIT NONE
     
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variable
  REAL(8), INTENT(out)   :: F(3,3)
      
! working variables
  INTEGER, DIMENSION(6)  :: nvec1 =  (/ 1, 1, 2, 3, 2, 3 /)
  INTEGER, DIMENSION(6)  :: nvec2 =  (/ 3, 2, 3, 2, 1, 1 /)
  REAL(8)                :: Dr(3)   , Dr2(3)
  REAL(8)                :: Sr(3,3) , Sr2(3,3)
  REAL(8)                :: dNxyz(3), dN2xyz(3)
  REAL(8)                :: U(3,3), U2(3,3), U3(3,3)    
  INTEGER                :: nn1, ivec1, ivec2
  REAL(8)                :: alpha1, alpha2, alpha3, alpha4
  REAL(8)                :: beta, gamma1, gamma2

  IF (nn > 2) THEN
     F = 0.D0
     RETURN
  END IF

  CALL grphpar(alpha1,alpha2,alpha3,alpha4,beta,gamma1,gamma2)
  
  IF(nn == 1) THEN      
     CALL drVec(n,m,iatom,Dr)
     CALL srmatrix(n,m,iatom,Sr)
     CALL dr2Vec(n,m,iatom,ivec,Dr2)
     CALL sr2matrix(n,m,iatom,ivec,Sr2)
     CALL dNxyzVec(n,m,iatom,ivec,nn,dNxyz)
     
     CALL outProd(3,dNxyz, 3,dNxyz, 3,U)
     CALL outProd(3,dNxyz, 3,Dr,    3,U2)
     CALL outProd(3,Dr2,   3,dNxyz, 3,U3)
     
     F = 2.D0*beta*(Sr+Sr2) - 8.D0*beta*U + 2.D0*beta*(U2-U3)

  END IF
      
  IF (nn == 2) THEN
     nn1 = 1
     ivec1 = nvec1(ivec)
     ivec2 = nvec2(ivec)
     CALL dNxyzVec(n,m,iatom,ivec1,nn1,dNxyz)
     CALL dN2xyzVec(n,m,iatom,ivec1,ivec2,dN2xyz)
     CALL outProd(3,dN2xyz,3,dNxyz,3,U)
    
     F = 2.D0*beta*U      

  END IF
      
  RETURN
END SUBROUTINE bondBendingF
!*******************************************************************************
!*******************************************************************************
SUBROUTINE danglingBondF(n,m,iatom,ivec,nn,F)
!===============================================================================
! Force constant matrices for dangling bond interactions 
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          neighbor index
!  nn            near neighbor shell
! Output       :
!  F(3,3)        force constant matrix (Joules/meter**2)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variable
  REAL(8), INTENT(out)   :: F(3,3)

! working variables
  INTEGER, DIMENSION(6)  :: nvec1 =  (/ 1, 1, 2, 3, 2, 3 /)
  REAL(8)                :: rho(3), U1(3,3), U2(3,3)   
  INTEGER                :: ivec1
  REAL(8)                :: alpha1, alpha2, alpha3, alpha4
  REAL(8)                :: beta, gamma1, gamma2      
  
  IF (nn > 2) THEN
     F = 0.D0
     RETURN
  END IF

  CALL grphpar(alpha1,alpha2,alpha3,alpha4,beta,gamma1,gamma2)
     
  IF (nn == 1) THEN
     CALL rhoVec(n,m,iatom,rho)
     CALL outProd(3,rho,3,rho,3,U1)       
     CALL rho2Vec(n,m,iatom,ivec, rho)
     CALL outProd(3,rho, 3,rho, 3,U2)
     
     F = 3.D0 * gamma1 * (U1+U2)

  END IF
      
  IF (nn == 2) THEN
     ivec1 = nvec1(ivec)
     CALL rho2Vec(n,m,iatom,ivec1, rho)
     CALL outProd(3,rho, 3,rho, 3,U1)
        
     F = -gamma1 * U1

  END IF
            
  RETURN
END SUBROUTINE danglingBondF
!*******************************************************************************
!*******************************************************************************
SUBROUTINE bondTwistF(n,m,iatom,ivec,nn,F)
!===============================================================================
! Force constant matrices for bond twisting interactions 
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          neighbor index
!  nn            near neighbor shell
! Output       :
!  F(3,3)        force constant matrix (Joules/meter**2)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variable
  REAL(8), INTENT(out)   :: F(3,3)
  
! working variables
  INTEGER, DIMENSION(6)  :: ibonda =  (/ 2, 1, 5, 6, 3, 4 /)      
  INTEGER, DIMENSION(3)  :: ibond1a = (/ 5, 1, 2 /)
  INTEGER, DIMENSION(3)  :: ibond2a = (/ 6, 3, 4 /)                      
  REAL(8)                :: Ur1(3), Ur2(3), F1(3,3), F2(3,3)     
  INTEGER                :: ibond, ibond1, ibond2     
  REAL(8)                :: alpha1, alpha2, alpha3, alpha4
  REAL(8)                :: beta, gamma1, gamma2
      
  IF (nn == 1) THEN
     F = 0.D0
     RETURN
  END IF
  
  CALL grphpar(alpha1,alpha2,alpha3,alpha4,beta,gamma1,gamma2)
      
  IF (nn == 2) THEN
     ibond = ibonda(ivec)
     CALL bondOutNormal(n,m,iatom,ibond,Ur1)
     CALL outProd(3,Ur1,3,Ur1,3,F1)

     F = gamma2*F1

  END IF
  
  IF (nn == 3) THEN
     ibond1 = ibond1a(ivec)
     CALL bondOutNormal(n,m,iatom,ibond1,Ur1)
     CALL outProd(3,Ur1,3,Ur1,3,F1)
        
     ibond2=ibond2a(ivec)
     CALL bondOutNormal(n,m,iatom,ibond2,Ur2)
     CALL outProd(3,Ur2,3,Ur2,3,F2) 

     F = gamma2*(F1+F2)                                  

  END IF
      
  IF (nn.EQ.4) THEN
     ibond=ivec
     CALL bondOutNormal(n,m, iatom,ibond, Ur1)        
     CALL outProd(3,Ur1, 3,Ur1, 3,F1)

     F = -gamma2*F1              

  END IF
  
  RETURN
END SUBROUTINE bondTwistF
!*******************************************************************************
!*******************************************************************************
SUBROUTINE bondOutNormal(n,m,iatom,ibond,Ur)
!===============================================================================
! outward normal vector at bond midpoints (dimensionless)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ibond         bond index (ibond=1..6)
! Output       :
!  Ur(3)         ourward normal at bond midpoint (dimensionless)
!=======================================================================
  IMPLICIT NONE
     
! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ibond

! output variable
  REAL(8), INTENT(out)   :: Ur(3)

! working variables      
  INTEGER, DIMENSION(6)  :: nvec1 = (/ 1, 1, 2, 3, 2, 3 /)
  REAL(8)                :: R1(3), R2(3), Rbond(3)     
  INTEGER                :: ier, ivec1
      
  ier = 0

! check for errors
  IF (iatom /= 1 .AND. iatom /= 2) THEN
     ier = 1
     WRITE (*,*) 'bondOutNormal err: invalid iatom:', iatom
  END IF
      
  IF (ibond < 1 .OR. ibond > 6) THEN
     ier=1
     WRITE(*,*) 'bondOutNormal err: invalid ibond:', ibond
  ENDIF
  
  IF (ier /= 0) STOP
      
! define the vector
  ivec1 = nvec1(ibond)
  CALL rxyzVec(n,m,iatom,ivec1,1,R1)
  CALL rxyzVec(n,m,iatom,ibond,2,R2)
  Rbond = (R1+R2)/2.D0
  CALL tubeRadUnitVec(Rbond,Ur)             
                                            
  RETURN
END SUBROUTINE bondOutNormal
!***********************************************************************
