!*******************************************************************************
!*******************************************************************************
! Project      : libswntElec.f90
!===============================================================================
! Purpose      :
! Calculate electronic states for a SWNT using extended tight-binding model (up to 3rd NN)
!-------------------------------------------------------------------------------
! Authors      : ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
!                Gary Sanders (sanders@phys.ufl.edu)
!                Daria Satco  (dasha.shatco@gmail.com)
! Latest Vers. : 2018.11.30
!-------------------------------------------------------------------------------
! Reference(s) :
! [1] Physical Properties of Carbon Nanotubes
!     R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
! [2] Carbon Nanotube Photophysics, G. G. Samsonidze MIT Ph.D. Thesis (2006)
!-------------------------------------------------------------------------------
! Contents     :
! - SUBROUTINE etbTubeEn(n,m,mu,rk,Ek)
! - SUBROUTINE etbTubeBand(n,m,mu,rk,Ek,Zk)
! - SUBROUTINE etbTubeEii(n,m,ii,rkii,eii,ierr)
! - SUBROUTINE etbTubeEgap(n,m,mu,rk,egap)
! - SUBROUTINE tubeElDOS(n,m,ne,Earray,DOS)
! - SUBROUTINE piHamOvlp(n,m,rk,mu,H,S)
! - SUBROUTINE tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)
! - SUBROUTINE etbPiTB3(n,m,rk,nout,Ek)
! - SUBROUTINE ChargeDensity(ne, Tempr, Efermi, Earray, DOS, charge)
! - SUBROUTINE QuantumCapacitance(ne, Tempr, Efermi, Earray, DOS, capacitance)
! - FUNCTION HppPi(r)
! - FUNCTION OppPi(r)
! - FUNCTION fermiLevel(n,m,Tempr,density)
! - FUNCTION fermiLevelFunc(Ef)
! - FUNCTION elecDensity(n,m,Tempr,Ef)
! - FUNCTION fermi(E,Ef,rkT)
!*******************************************************************************
!*******************************************************************************
SUBROUTINE etbTubeEn(n,m,mu,rk,Ek)
!===============================================================================
! Energy bands for (n,m) carbon nanotube in extended tight binding model
! WITHOUT eigenvector (wavefunctions) calculations
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  mu            cutting lines (0...N_hex-1)
!  rk            electron wavevector (1/A) (0 < k < pi/T)
! Output       :
!  Ek(2)         electronic energies in ascending order (eV)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, mu
  REAL(8), INTENT(in)    :: rk    !(1/A)

! output variable
  REAL(8), INTENT(out)   :: Ek(2) !(eV)      
      
! working variables
  COMPLEX(8)             :: H(2,2), S(2,2), Zk(2,2)

! lapack driver variables
  INTEGER                :: matz, il, iu, nout1

! option flag (0=evalues, 1=evalues+evectors)
  matz = 0 ! calculate evalues only

! if il or iu <= 0, all eigenvalues are returned
  il = 0   ! lower indices of desired eigenvalues
  iu = 0   ! upper indices of desired eigenvalues
        
  CALL piHamOvlp(n,m,rk,mu,H,S)

! solveHam(n,ldh,ham,ldo,ovlp,matz,il,iu,nout,w,ldz,z)
! nout     ->  number of eigenvalues returned
! w(n)     ->  eigenvalues in ascending order
! ldz      ->  leading dimension of z
! z(ldz,n) ->  complex eigenvectors if matz = 1
  CALL solveHam(2,2,H,2,S,matz,il,iu,nout1,Ek,2,Zk)

END SUBROUTINE etbTubeEn
!*******************************************************************************
!*******************************************************************************
SUBROUTINE etbTubeBand(n,m,mu,rk,Ek,Zk)
!===============================================================================
! Energy bands for (n,m) carbon nanotube in extended tight binding model
! WITH eigenvector (wavefunction) calculations
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  mu            cutting lines (0...N_hex-1)
!  rk            electron wavevector (1/A) (0 < k < pi/T)
! Output       :
!  Ek(2)         electronic energies in ascending order (eV)
!  Zk(2,2)       electronic wavefunctions (dimensionless)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, mu
  REAL(8), INTENT(in)    :: rk    !(1/A)

! output variable
  REAL(8), INTENT(out)   :: Ek(2) !(eV)      
      
! working variables
  COMPLEX(8)            :: H(2,2), S(2,2), Zk(2,2)

! lapack driver variables
  INTEGER                :: matz, il, iu, nout1

! option flag (0=evalues, 1=evalues+evectors)
  matz = 1 ! calculate evalues+evectors

! if il or iu <= 0, all eigenvalues are returned
  il = 0   ! lower indices of desired eigenvalues
  iu = 0   ! upper indices of desired eigenvalues
    
  CALL piHamOvlp(n,m,rk,mu,H,S)
! solveHam(n,ldh,ham,ldo,ovlp,matz,il,iu,nout,w,ldz,z)
! nout     ->  number of eigenvalues returned
! w(n)     ->  eigenvalues in ascending order
! ldz      ->  leading dimension of z
! z(ldz,n) ->  complex eigenvectors if matz = 1
  
  CALL solveHam(2,2,H,2,S,matz,il,iu,nout1,Ek,2,Zk)
            
END SUBROUTINE etbTubeBand
!*******************************************************************************
!*******************************************************************************
SUBROUTINE etbTubeEii(n,m,ii,rkii,eii,ierr)
!===============================================================================
! find the E_{ii} energy gaps and the k value at which they occur
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  ii            value of i = 1,2,3...
! Output       :
!  rkii          k value for Eii transition (1/Angstrom) (0 ... pi/T)
!  eii           Eii transition energy (eV)
!  ierr          0=normal completion, 1=invalid ii            
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)  :: n, m, ii

! output variables
  REAL(8), INTENT(out) :: rkii, eii
  INTEGER, INTENT(out) :: ierr

! working variables 
  REAL(8), SAVE, DIMENSION(:), ALLOCATABLE :: rks
  REAL(8), SAVE, DIMENSION(:), ALLOCATABLE :: egaps      

  INTEGER, SAVE        :: nch = 0
  INTEGER, SAVE        :: mch = 0
  INTEGER, SAVE        :: metal
  INTEGER, SAVE        :: nmu
  INTEGER              :: nhex, nHexagon, mu, indx
  REAL(8)              :: rk, egap
  
  IF (n /= nch .OR. m /= mch) THEN
     nch   = n
     mch   = m
     metal = MOD(n-m,3)      
     nhex  = nHexagon(n,m)
     
     IF (ALLOCATED(rks)) DEALLOCATE(rks)
     ALLOCATE(rks(nhex))        
     
     IF (ALLOCATED(egaps)) DEALLOCATE(egaps)
     ALLOCATE(egaps(nhex))
     
     nmu = nhex/2
     DO mu = 1, nmu
        CALL etbTubeEgap(n,m,mu,rk,egap)
        rks(mu)   = rk
        egaps(mu) = egap
     END DO
     
     CALL sort2(nmu,egaps,rks)
  ENDIF
            
  IF (metal == 0) THEN
     indx = ii+1
  ELSE
     indx = ii
  END IF
  
  IF (indx < 1 .OR. indx > nmu) THEN
     ierr = 1
     rkii =-999.
     eii  =-999.
  ELSE
     ierr = 0
     rkii = rks(indx)
     eii  = egaps(indx)
  END IF
  
  RETURN

END SUBROUTINE etbTubeEii
!*******************************************************************************
!*******************************************************************************
SUBROUTINE etbTubeEgap(n,m,mu,rk,egap)
!===============================================================================
! find the energy gap in the mu'th manifold and k value at which it occurs
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  mu            labels electronic manifolds (0...N_hex-1)
! Output       :
!  rk            k at energy gap (0 ... pi/T) (1/Angstroms)
!  egap          pi band energy gap (eV)          
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, mu

! output variables
  REAL(8), INTENT(out)   :: rk, egap

! working variables and parameters
  INTEGER, PARAMETER     :: nk = 41
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0
  REAL(8), DIMENSION(nk) :: rka, Ega  !(nk)
  REAL(8)                :: rkmin, rkmax, dk
  REAL(8)                :: trLength, fgap
  INTEGER                :: n1, m1, mu1, k
  COMMON /fgapcom/ n1,m1,mu1      

! pass variables to fgap using common block
  n1  = n
  m1  = m
  mu1 = mu

! bracket a minimum
  rkmax = pi/trLength(n,m)
  rkmin = -rkmax
      
  dk = (rkmax - rkmin) / (nk-1.D0)
  DO k = 1, nk
     rk = rkmin + (k-1)*dk
     rka(k) = rk
     Ega(k) = fgap(rk)
  END DO
      
  egap = MINVAL(Ega)
  rk   = ABS(rka(MINLOC(Ega,dim=1)))

END SUBROUTINE etbTubeEgap
!-------------------------------------------------------------------------------
REAL(8) FUNCTION fgap(rk)

  IMPLICIT NONE

  REAL(8)                :: Ek(2), rk
  INTEGER                :: n1, m1, mu1
  COMMON /fgapcom/ n1,m1,mu1
      
  CALL etbTubeEn(n1,m1,mu1,rk,Ek)      
      
  fgap = Ek(2) - Ek(1)
  RETURN
END FUNCTION fgap
!-------------------------------------------------------------------------------
FUNCTION golden(ax,bx,cx,tol,xmin)
  REAL(8)                :: golden, ax, bx, cx, tol, xmin, C
  REAL(8)                :: f1, f2, x0, x1, x2, x3
  REAL(8)                :: fgap
  REAL(8), PARAMETER     :: R = .61803399D0

  C  = 1.D0 - R
  x0 = ax
  x3 = cx
  IF (ABS(cx-bx) > ABS(bx-ax))THEN
     x1 = bx
     x2 = bx + C*(cx-bx)
  ELSE
     x2 = bx
     x1 = bx - C*(bx-ax)
  END IF
  f1 = fgap(x1)
  f2 = fgap(x2)

1 IF ( ABS(x3-x0) > tol*(ABS(x1)+ABS(x2)) )THEN
     IF (f2 < f1) THEN
        x0 = x1
        x1 = x2
        x2 = R*x1 + C*x3
        f1 = f2
        f2 = fgap(x2)
     ELSE
        x3 = x2
        x2 = x1
        x1 = R*x2 + C*x0
        f2 = f1
        f1 = fgap(x1)
     END IF
     GOTO 1
  END IF
  
  IF (f1 < f2)THEN
     golden = f1
     xmin   = x1
  ELSE
     golden = f2
     xmin   = x2
  END IF
  RETURN

END FUNCTION golden
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tubeElDOS(n,m,ne,Earray,DOS)
!===============================================================================
! Electron density of states per carbon atom for an (n,m) carbon
! nanotube (states/carbon atom/eV)
!
! Note: Since there is one p_z orbital per carbon atom site, the total
! number of pi electrons per atom is 2. We thus have the sum rule:
!
!            Int( DOS(E), E=-infinity..infinity) = 2
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  ne            number of energies
!  Earray(ne)    array of energies (eV)
! Output       :
!  DOS(ne)       electron density of states (states/atom/eV)
!=======================================================================
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
  REAL(8)                :: rka(nk), Ek(2)
  REAL(8), ALLOCATABLE   :: Enk(:,:) !(nk,2*nhex)
     
  INTEGER                :: nout, nhex, nHexagon
  INTEGER                :: mu, k, indexx, i, ie, nn
      
  REAL(8)                :: T, trLength
  REAL(8)                :: rk, rkmin, rkmax, dk, fwhm, fwhm1, E, DSn
            
! allocate storage
  nhex = nHexagon(n,m)
  nout = 2*nhex      
  ALLOCATE(Enk(nk,nout))
      
! evenly spaced k points from 0 to pi/T (1/Angstroms)
  T     = trLength(n,m)
  rkmin = 0.D0
  rkmax = pi/T
  dk    = (rkmax-rkmin) / (DBLE(nk) - 1.D0)
  DO k = 1, nk
     rka(k) = rkmin + (k-1)*dk
  END DO
      
! electron energy bands at evenly space k points (eV)
  DO k = 1, nk
     rk = rka(k)
     indexx = 0
     DO mu = 1, nhex
        CALL etbTubeEn(n,m,mu,rk,Ek)
        DO i = 1, 2
           indexx = indexx + 1
           Enk(k,indexx) = Ek(i)
        END DO
     END DO
  END DO
      
! find FWHM linewidth based on energy array (eV)
  fwhm = ABS(Earray(2) - Earray(1))
  DO ie = 1, ne-1
     fwhm1 = ABS(Earray(ie) - Earray(ie+1))
     IF(fwhm1 < fwhm) fwhm = fwhm1
  END DO
      
! accumulate density of states per unit length
  DO ie = 1, ne
     DOS(ie) = 0.D0
     DO nn = 1, nout      
        E = Earray(ie)
        CALL dos1Dgauss(nout,nk,rka,nk,Enk,E,fwhm,nn,DSn)
        DOS(ie) = DOS(ie) + DSn
     END DO
  END DO
      
! convert to density of states/Carbon_atom/eV
  DO ie = 1, ne
     DOS(ie) = 2.D0*(T/DBLE(nout)) * DOS(ie)
  END DO
  
  DEALLOCATE(Enk)      
 
END SUBROUTINE tubeElDOS
!*******************************************************************************
!*******************************************************************************
SUBROUTINE piHamOvlp(n,m,rk,mu,H,S)
!===============================================================================
! Hamiltonian and Overlap Matrices for mu'th cutting line for Pi bands
! of carbon nanotubes considering long-range interactions 
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  rk            nanotube wavevector (1/A)
!  mu            labels electronic cutting lines (0...N_hex-1)
! Output       :
!  H(2,2)        complex hamiltonian matix (eV)
!  S(2,2)        complex overlap matrix (dimensionless)
!===============================================================================
  IMPLICIT NONE
            
! input variables
  INTEGER, INTENT(in)    :: n,m,mu
  REAL(8), INTENT(in)    :: rk

! output variables
  COMPLEX(8),INTENT(out) :: H (2,2), S(2,2)

! working variables
  COMPLEX(8), SAVE       :: ci = (0.D0, 1.D0)   
  COMPLEX(8)             :: css, csh, expphi              
            
  COMPLEX(8)             :: Htemp(2,2), Stemp(2,2)
      
  INTEGER                :: i, j, iatom, jatom, nn, ivec
  REAL(8)                :: ham, ovlp, phi
            
! on-site contributions to hamiltonian and overlap
  H = 0.D0
  S = 0.D0

  DO iatom = 1, 2
     CALL tbAtomHamOvlp(n,m,iatom,1,0,ham,ovlp) 
     H(iatom,iatom) = ham
     S(iatom,iatom) = ovlp
  END DO
      
! nearest neighbor contributions to hamiltonian and overlap
  nn = 1      
      
! H(AB) and S(AB)      
  iatom = 1
  jatom = 2
  csh   = 0.D0
  css   = 0.D0

  DO ivec = 1, 3
     CALL tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)      
     CALL phij(n,m,iatom,ivec,nn,rk,mu,phi)        
     expphi = CDEXP(ci*phi)
     csh    = csh + expphi*ham
     css    = css + expphi*ovlp
  END DO
  H(iatom,jatom) = H(iatom,jatom) + csh
  S(iatom,jatom) = S(iatom,jatom) + css
      
! H(BA) and S(BA)
  iatom = 2
  jatom = 1
  csh   = 0.D0
  css   = 0.D0
  DO ivec = 1, 3
     CALL tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)      
     CALL phij(n,m,iatom,ivec,nn,rk,mu,phi)
     expphi = CDEXP(ci*phi)
     csh    = csh + expphi*ham
     css    = css + expphi*ovlp      
  END DO
  H(iatom,jatom) = H(iatom,jatom) + csh
  S(iatom,jatom) = S(iatom,jatom) + css
      
! second neighbor contributions to hamiltonian and overlap
  nn = 2
      
! H(AA) and S(AA)      
  iatom = 1
  jatom = 1
  csh   = 0.D0
  css   = 0.D0
  DO ivec = 1, 6
     CALL tbAtomHamOvlp(n,m, iatom,ivec,nn, ham,ovlp)      
     CALL phij(n,m,iatom,ivec,nn,rk,mu,phi)
     expphi = CDEXP(ci*phi)
     csh    = csh+expphi*ham
     css    = css+expphi*ovlp      
  END DO
  H(iatom,jatom) = H(iatom,jatom) + csh
  S(iatom,jatom) = S(iatom,jatom) + css
      
! H(BB) and S(BB)      
  iatom = 2
  jatom = 2
  csh   = 0.D0      
  css   = 0.D0
  DO ivec = 1, 6
     CALL tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)      
     CALL phij(n,m,iatom,ivec,nn,rk,mu,phi)
     expphi = CDEXP(ci*phi)
     csh    = csh + expphi*ham
     css    = css + expphi*ovlp      
  END DO
  H(iatom,jatom) = H(iatom,jatom) + csh
  S(iatom,jatom) = S(iatom,jatom) + css
      
! third neighbor contributions to hamiltonian and overlap
  nn = 3      
      
! H(AB) and S(AB)      
  iatom = 1
  jatom = 2
  csh   = 0.D0
  css   = 0.D0
  DO ivec = 1, 3
     CALL tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)      
     CALL phij(n,m,iatom,ivec,nn,rk,mu,phi)
     expphi = CDEXP(ci*phi)
     csh    = csh + expphi*ham
     css    = css + expphi*ovlp      
  END DO
  H(iatom,jatom) = H(iatom,jatom) + csh
  S(iatom,jatom) = S(iatom,jatom) + css
      
! H(BA) and S(BA)
  iatom = 2
  jatom = 1
  csh   = 0.D0
  css   = 0.D0
  DO ivec = 1, 3
     CALL tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)      
     CALL phij(n,m,iatom,ivec,nn,rk,mu,phi)
     expphi = CDEXP(ci*phi)
     csh    = csh + expphi*ham
     css    = css + expphi*ovlp      
  END DO
  H(iatom,jatom) = H(iatom,jatom) + csh
  S(iatom,jatom) = S(iatom,jatom) + css
      
! correct roundoff errors
  Htemp = H
  Stemp = S
  DO i = 1, 2
     DO j = 1, 2
        H(i,j) = (Htemp(i,j) + CONJG(Htemp(j,i))) / 2.D0
        S(i,j) = (Stemp(i,j) + CONJG(Stemp(j,i))) / 2.D0
     END DO
  END DO
                                          
END SUBROUTINE piHamOvlp
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbAtomHamOvlp(n,m,iatom,ivec,nn,ham,ovlp)
!===============================================================================
! Calculate the atomic hamiltonian and overlap matrix elements between an A
! or B atom in the two-atom unit cell and an atom in a near neighbor shell
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iatom         specifies atom in two atom unit cell (1=A,2=B)
!  ivec          index for nearest neighbor vector in the shell
!  nn            neighbor index nn = 0,1,2,3,4
! Output       :
!  ham           hamiltonian matrix element (eV)
!  ovlp          overlap matrix element (dimensionless)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iatom, ivec, nn

! output variables
  REAL(8), INTENT(out)   :: ham, ovlp
      
! working variables
  REAL(8)                :: R0(3),R1(3),dR(3)

  REAL(8), SAVE, DIMENSION(2,6, 0:4) :: Htable !(iatom,ivec,nn)
  REAL(8), SAVE, DIMENSION(2,6, 0:4) :: Stable !(iatom,ivec,nn)
  INTEGER, SAVE, DIMENSION(0:4)      :: nvecs = (/ 1, 3, 6, 3, 6 /)     

  INTEGER, SAVE          :: nch = 0
  INTEGER, SAVE          :: mch = 0
      
  INTEGER                :: iiatom, nnn, iivec

! function declaration
  REAL(8)                :: r, vecLength     
  REAL(8)                :: HppPi, OppPi
      
! update hamiltonian and overlap lookup tables if (n,m) has changed
  IF(n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m
      
     DO iiatom = 1, 2
        DO nnn = 0, 4
           DO iivec = 1, nvecs(nnn)
      
              IF(nnn == 0) THEN
                 ham  = -.7078D0
                 ovlp = 1.D0
              ELSE
                 CALL rxyzVec(n,m,iiatom,1,0,R0)
                 CALL rxyzVec(n,m,iiatom,iivec,nnn,R1)
                 dR   = R1 - R0
                 r    = vecLength(3,dR)
                 ham  = HppPi(r)
                 ovlp = OppPi(r)
              END IF
        
              Htable(iiatom,iivec,nnn) = ham
              Stable(iiatom,iivec,nnn) = ovlp        
              
           END DO
        END DO
     END DO
      
  END IF

! return hamiltonian and overlap matrix elements from lookup tables
  iivec = ivec
  IF(nn == 0) iivec = 1

  ham  = Htable(iatom,iivec,nn)
  ovlp = Stable(iatom,iivec,nn)      
      
END SUBROUTINE tbAtomHamOvlp
!*******************************************************************************
!*******************************************************************************
SUBROUTINE etbPiTB3(n,m,rk,nout,Ek)
!===============================================================================
! Carbon nanotube energy bands in extended tight binding model
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2) basis
!  rk            nanotube wavevector (1/A)
! Output       :
!  nout          number of electronic energies returned
!  Ek(nout)      electronic energies in ascending order (eV)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  REAL(8), INTENT(in)    :: rk       !(1/A)

! output variables
  INTEGER, INTENT(out)              :: nout
  REAL(8), DIMENSION(*), INTENT(out):: Ek !(nout) !(eV) 
  
! working variables
  INTEGER, ALLOCATABLE   :: indx(:)   !(nout)
  REAL(8), ALLOCATABLE   :: Earray(:) !(nout)
      
  REAL(8)                :: w(2)
  COMPLEX(8)             :: H(2,2),S(2,2),z(2,2)
      
  INTEGER                :: nhex, nHexagon, i
  INTEGER                :: matz, il, iu, mu, nout1, id1, id2

! number of hexagons in nanotube unit cell
  nhex = nHexagon(n,m)
  nout = 2*nhex
      
! loop over electronic manifolds
  matz = 0
  il = 0
  iu = 0
  DO mu = 1, nhex       
     CALL piHamOvlp(n,m,rk,mu,H,S)
     CALL solveHam(2,2,H,2,S,matz,il,iu,nout1,w,2,z)
     
     id1 = 2*mu-1
     id2 = 2*mu
        
     Ek(id1) = w(1)
     Ek(id2) = w(2)        
                      
  END DO
      
! sort energies and manifolds
  ALLOCATE(indx(nout))
  ALLOCATE(Earray(nout))

  DO i = 1, nout
     Earray(i) = Ek(i)
  END DO
  CALL indexx(nout,Earray,indx)
  DO i = 1, nout
     Ek(i) = Earray(indx(i))
  END DO

  DEALLOCATE(indx)
  DEALLOCATE(Earray)
      
  RETURN
END SUBROUTINE etbPiTB3
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION HppPi(r)
!===============================================================================
! pp_pi hopping matrix element vs interatomic separation
!-------------------------------------------------------------------------------
! Input        :
!  r             interatomic separation (Angstroms)
! Output       :
!  HppPi         pp-pi hopping matrix (eV)
!===============================================================================
  IMPLICIT NONE

! input variable
  REAL(8), INTENT(in)    :: r
      
! working variable and parameter
  REAL(8), PARAMETER     :: a0 = .52917721D0    !(angstroms)
  REAL(8), PARAMETER     :: Eh = 27.21138D0    
  REAL(8)                :: a, b, x, ss, chebev !chebyshev function     
  REAL(8)                :: c(10)

  c(1) =-.3793837D0
  c(2) = .3204470D0
  c(3) =-.1956799D0
  c(4) = .0883986D0
  c(5) =-.0300733D0
  c(6) = .0074465D0
  c(7) =-.0008563D0
  c(8) =-.0004453D0
  c(9) = .0003842D0
  c(10)=-.0001855D0
            
  a = 1.D0
  b = 7.D0
      
  x = r / a0
  IF (x > b) THEN
     ss = 0.D0
  ELSE
     ss = Eh*chebev(a,b,c,10,x)
  END IF

! return HppPi
  HppPi = ss
    
END FUNCTION HppPi
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION OppPi(r)
!===============================================================================
! pp_pi overlap matrix element vs interatomic separation
!-------------------------------------------------------------------------------
! Input        :
!  r             interatomic separation (Angstroms)
! Output       :
!  OppPi         pp-pi overlap matrix (dimensionless)
!===============================================================================
  IMPLICIT NONE

! input variable
  REAL(8), INTENT(in)    :: r

! working variable and parameter
  REAL(8), PARAMETER     :: a0 = .52917721D0
  REAL(8)                :: a, b, x, ss, chebev
  REAL(8)                :: c(10)

  c(1) = .3715732D0
  c(2) =-.3070867D0
  c(3) = .1707304D0
  c(4) =-.0581555D0
  c(5) = .0061645D0
  c(6) = .0051460D0
  c(7) =-.0032776D0
  c(8) = .0009119D0
  c(9) =-.0001265D0
  c(10)=-.0000227D0
            
  a = 1.D0
  b = 7.D0
  
  x = r / a0
  IF (x > b) THEN
     ss = 0.D0
  ELSE
     ss = chebev(a,b,c,10,x)
  END IF

! return OppPi
  OppPi = ss
      
END FUNCTION OppPi
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION fermiLevel(n,m,Tempr,density)
!===============================================================================
! Fermi energy (eV) for (n,m) carbon nanotube
! vs temperature (deg K) and net electron density (electrons/Angststrom)
!===============================================================================
! Input        :
!  n,m           chiral vector coordinates in (a1,a2) basis
!  Tempr         lattice temperature (deg K)
!  density       net electron density per unit length (1/Angstroms)
! Output       :
!  fermiLevel    tube Fermi energy, Ef(eV)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: n,m
  REAL(8), INTENT(in)    :: Tempr,density
  
! working variables
  LOGICAL                :: success
  
  INTEGER                :: n1,m1
  REAL(8)                :: Tempr1, density1

  COMMON /fermilevelcom/ n1,m1,Tempr1,density1 ! common variables
      
  REAL(8)                :: x1, x2, xacc
           
  REAL(8), EXTERNAL      :: fermiLevelFunc
  REAL(8), EXTERNAL      :: rtbis      
      
  n1 = n
  m1 = m
  Tempr1   = Tempr
  density1 = density            

! bracket a root
  x1 =-5.D0
  x2 = 5.D0
  CALL zbrac(fermiLevelFunc,x1,x2,success)
  IF(success .EQV. .FALSE.) THEN
     WRITE (*,*) 'fermiLevel err:'
     WRITE (*,*) 'zbrac fails to bracket a root:'
     WRITE (*,*) 'x1,f(x1):', x1, fermiLevelFunc(x1)
     WRITE (*,*) 'x2,f(x2):', x2, fermiLevelFunc(x2)        
     STOP
  END IF
      
! find fermi energy by bisection
  xacc = 1.D-6           
  fermiLevel = rtbis(fermiLevelFunc,x1,x2,xacc)      

END FUNCTION fermiLevel
!-------------------------------------------------------------------------------
REAL(8) FUNCTION fermiLevelFunc(Ef)
  IMPLICIT NONE
  REAL(8), INTENT(in)    :: Ef
      
  INTEGER                :: n1, m1
  REAL(8)                :: Tempr1, density1
  COMMON /fermilevelcom/ n1,m1,Tempr1,density1 ! common variables
      
  REAL(8)   elecDensity 
  
  fermiLevelFunc = elecDensity(n1,m1,Tempr1,Ef) - density1

END FUNCTION fermiLevelFunc
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION elecDensity(n,m, Tempr,Ef)
!===============================================================================
! net electron density (electrons/Angststrom) for (n,m) carbon nanotube
! vs temperature (deg K) and Fermi level (eV)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2) basis
!  Tempr         lattice temperature (deg K)
!  Ef            Fermi level (eV)
! Output       :
!  elecDensity   net electron density per unit length (1/Angstroms)
!=======================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  REAL(8), INTENT(in)    :: Tempr, Ef

! working variables and parameters
  INTEGER, PARAMETER     :: nk = 81
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0

  REAL(8), ALLOCATABLE   :: En(:)        !(2*Nhex)
  REAL(8), SAVE, ALLOCATABLE :: Enk(:,:) !(2*Nhex,nk)
  REAL(8), SAVE, ALLOCATABLE :: rka(:)   !(nk)
      
  INTEGER, SAVE          :: ifirst = 1      
  INTEGER, SAVE          :: mch = 0
  INTEGER, SAVE          :: nch = 0
  INTEGER, SAVE          :: nhex
  INTEGER, SAVE          :: nout
      
  INTEGER                :: nHexagon,k,iband,i
  
  REAL(8)                :: rkmin, rkmax, dk, rk, rkT, dkk, ss, Ekk, fermi
  REAL(8)                :: trLength
     
! nanotube energy bands (eV)
! check for errors then allocate
  IF (n /= nch .OR. m /= mch) THEN
     nch = n
     mch = m        
        
     nhex = nHexagon(n,m)
     nout = 2*nhex

     IF (ifirst == 0) DEALLOCATE(rka)
     ALLOCATE(rka(nk))
        
     IF (ifirst == 0) DEALLOCATE(Enk)
     ALLOCATE(Enk(nout,nk))     
      
! define k point array (1/A)
     rkmax = pi/trLength(n,m)
     rkmin = 0.D0
     dk = (rkmax - rkmin) / (nk - 1.D0)
     
     DO k = 1, nk
        rka(k) = rkmin + (k-1)*dk
     END DO

! compute energy bands En(k) (eV)
     ALLOCATE(En(nout))
     DO k = 1, nk
        rk = rka(k)
        CALL etbPiTB3(n,m,rk,nout,En)
        DO iband = 1, nout
           Enk(iband,k) = En(iband)
        END DO
     END DO
     DEALLOCATE(En)      
  END IF
      
! integrate over k to obtain electron density per unit length
  rkT = .025853D0*(Tempr/300.D0)       
  dk = ABS(rka(2)-rka(1))
      
  ss = 0.D0
  DO k = 1, nk
     dkk = dk
     IF(k == 1 .OR. k == nk) dkk = dk/2.D0
     DO i = 1, nhex
        Ekk = Enk(i,k)
        ss = ss + dkk*(fermi(Ekk,Ef,rkT)-1.D0)
     END DO
      
     DO i = nhex + 1, nout
        Ekk = Enk(i,k)
        ss = ss + dkk*fermi(Ekk,Ef,rkT)
     END DO
  END DO
  elecDensity = 2.D0*ss/pi      
             
  ifirst = 0

END FUNCTION elecDensity
!*******************************************************************************
!*******************************************************************************
REAL(8) FUNCTION fermi(E,Ef,rkT)
!===============================================================================
! fermi-dirac distribution function for electrons
!             f(E,Ef,kT) = 1 / [ 1+exp( (E-Ef)/kT) ]
!-------------------------------------------------------------------------------
! Input        :
!  E             electron energy (eV)
!  Ef            fermi level (eV)
!  rkT           thermal energy (eV)
! Output       : 
!  fermi         electron occupation probability (dimensionless)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  REAL(8), INTENT(in)    :: E, Ef, rkT
  
! working variables
  REAL(8),  PARAMETER    :: etol = 80.D0
      
  REAL(8)                :: arg
  
  IF (rkT == 0.D0) THEN
     IF (E > Ef) THEN
        fermi = 0.D0
     ELSE IF (E == Ef) THEN
        fermi = .5D0
     ELSE
        fermi = 1.D0
     END IF
  ELSE
     arg = (E - Ef) / rkT
     IF (arg > etol) THEN
        fermi = 0.D0
     ELSE IF (arg.LT.-etol) THEN
        fermi = 1.D0
     ELSE
        fermi=1.D0/(1.D0+EXP(arg))
     END IF
  END IF
  
  RETURN
  
END FUNCTION fermi
!*******************************************************************************
!*******************************************************************************
SUBROUTINE ChargeDensity(ne, Tempr, Efermi, Earray, DOS, charge)
!===============================================================================
! calculates charge (eFermi > 0) and hole (eFermi < 0) density in CNT
!-------------------------------------------------------------------------------
! Input        :
!  ne            number of energies in Earray
!  Tempr         lattice temperature (deg K)
!  Efermi        fermi level (eV)
!  Earray       vector of energies - x axis (eV)
!  DOS           vector of DOS (states/atom/eV)
! Output       :
!  charge        number of electrons (1/atom)
!===============================================================================
IMPLICIT NONE
! input variables
  INTEGER, INTENT(in)    :: ne

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in)    :: Earray(ne)
  REAL(8), INTENT(in)    :: DOS(ne)

! output variable
  REAL(8), INTENT(out)   :: charge

! working variables
  REAL(8)                :: fermi
  REAL(8)                :: fermiDist, rkT, dE
  INTEGER                :: i

  rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
  dE = abs(Earray(2) - Earray(1))
  charge = 0.D0

  DO i=1,ne
    IF (Efermi > 0.D0) THEN
        IF (Earray(i) .ge. 0.D0) THEN
            fermiDist = fermi(Earray(i),Efermi,rkT)
            charge = charge + dE * DOS(i) * fermiDist
        END IF
    ELSE
        IF (Earray(i) .le. 0.D0) THEN
            fermiDist = fermi(Earray(i),Efermi,rkT)
            charge = charge + dE * DOS(i) * (1.D0 - fermiDist)
        END IF
    END IF
  END DO

END SUBROUTINE ChargeDensity
!*******************************************************************************
!*******************************************************************************
SUBROUTINE QuantumCapacitance(ne, Tempr, Efermi, Earray, DOS, capacitance)
!===============================================================================
! calculates quantum capacitance in CNT
!-------------------------------------------------------------------------------
! Input        :
!  ne            number of energies in Earray
!  Tempr         lattice temperature (deg K)
!  Efermi        fermi level (eV)
!  Earray        vector of energies - x axis (eV)
!  DOS           vector of DOS (states/atom/eV)
! Output       :
!  capacitance   quantum capacitance
!===============================================================================
IMPLICIT NONE
! input variables
  INTEGER, INTENT(in)    :: ne

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in)    :: Earray(ne)
  REAL(8), INTENT(in)    :: DOS(ne)

! output variable
  REAL(8), INTENT(out)   :: capacitance

! working variables
  REAL(8)                :: fermiDistDeriv, ch, rkT, dE
  INTEGER                :: i

  rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
  dE = abs(Earray(2) - Earray(1))
  capacitance = 0.D0

  DO i=1,ne
     ch = exp( (Earray(i) - Efermi)/(2*rkT) ) + exp( - (Earray(i) - Efermi)/(2*rkT) )
     fermiDistDeriv = 1.D0/rkT * 1/ch**2
     capacitance = capacitance + dE * DOS(i) * fermiDistDeriv
  END DO

END SUBROUTINE QuantumCapacitance
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tubeJDOS(n,m,nhex,nk,Enk,rka,ne,Earray,JDOSinter,JDOSintraV, JDOSintraC )
!===============================================================================
! Joint density of states per carbon atom for an (n,m) carbon
! nanotube (states/carbon atom/eV)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  ne            number of energies
!  Earray(ne)    array of energies (eV)
! Output       :
!  JDOSinter(ne)       joint density of states (states/atom/eV)
!  JDOSintra(ne)       joint density of states (states/atom/eV)
!=======================================================================
  IMPLICIT NONE

! parameters
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  REAL(8), INTENT(in)    :: rka(nk)
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: Earray(ne)

! output variable
  REAL(8), INTENT(out)   :: JDOSinter(ne), JDOSintraC(ne), JDOSintraV(ne)

! working variables
  REAL(8), ALLOCATABLE   :: EnkInter(:,:)  !(nk,nhex)
  REAL(8), ALLOCATABLE   :: EnkIntraC(:,:) !(nk,2*(nhex-1))
  REAL(8), ALLOCATABLE   :: EnkIntraV(:,:) !(nk,2*(nhex-1))

  INTEGER                :: noutInter, noutIntra
  INTEGER                :: mu1, mu2, k, indexx1, indexx2, i, ie, nn, ii

  REAL(8)                :: T, trLength
  REAL(8)                :: rk, rkmin, rkmax, dk, fwhm, fwhm1, E, DSn


! allocate storage
  noutInter = nhex
  noutIntra = 2*(nhex-1)
  ALLOCATE(EnkInter(nk,noutInter))
  ALLOCATE(EnkIntraC(nk,noutIntra))
  ALLOCATE(EnkIntraV(nk,noutIntra))

  EnkInter = 0.D0
  EnkIntraC = 0.D0
  EnkIntraV = 0.D0

! electron energy bands at evenly space k points (eV)
  DO k = 1, nk
    indexx1 = 1
    indexx2 = 1
    DO mu1 = 1, nhex
        DO i = -1, 1
            mu2 = mu1 + i

            IF ( (mu2 < 1) .OR. (mu2 > nhex) ) CYCLE

            ! interband
            IF ( mu1 == mu2 ) THEN
                EnkInter(k,indexx1) = Enk(2,mu2,k) - Enk(1,mu1,k)
                indexx1 = indexx1 + 1
            ELSE
            ! intraband
                EnkIntraV(k,indexx2) = Enk(2,mu2,k) - Enk(2,mu1,k) ! valence
                EnkIntraC(k,indexx2) = Enk(1,mu2,k) - Enk(1,mu1,k) ! conduction
                indexx2 = indexx2 + 1
            END IF

        END DO
    END DO
  END DO

! find FWHM linewidth based on energy array (eV)
  fwhm = ABS(Earray(2) - Earray(1))
  DO ie = 1, ne-1
     fwhm1 = ABS(Earray(ie) - Earray(ie+1))
     IF(fwhm1 < fwhm) fwhm = fwhm1
  END DO

! accumulate density of states per unit length
  DO ie = 1, ne
     JDOSinter(ie) = 0.D0
     JDOSintraC(ie) = 0.D0
     JDOSintraV(ie) = 0.D0
     DO nn = 1, noutInter
        E = Earray(ie)
        CALL dos1Dgauss(noutInter,nk,rka,nk,EnkInter,E,fwhm,nn,DSn)
        JDOSinter(ie) = JDOSinter(ie) + DSn
     END DO
     DO nn = 1, noutIntra
        E = Earray(ie)
        CALL dos1Dgauss(noutIntra,nk,rka,nk,EnkIntraC,E,fwhm,nn,DSn)
        JDOSintraC(ie) = JDOSintraC(ie) + DSn
     END DO
     DO nn = 1, noutIntra
        E = Earray(ie)
        CALL dos1Dgauss(noutIntra,nk,rka,nk,EnkIntraV,E,fwhm,nn,DSn)
        JDOSintraV(ie) = JDOSintraV(ie) + DSn
     END DO
  END DO

! convert to density of states/Carbon_atom/eV
  T = trLength(n,m)
  DO ie = 1, ne
     JDOSinter(ie) = (T/DBLE(noutInter)) * JDOSinter(ie)
     JDOSintraC(ie) = (T/DBLE(noutIntra)) * JDOSintraC(ie)
     JDOSintraV(ie) = (T/DBLE(noutIntra)) * JDOSintraV(ie)
  END DO

  DEALLOCATE(EnkInter, EnkIntraC, EnkIntraV)

END SUBROUTINE tubeJDOS
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tubeJDOSij(mu1,mu2,n,m,nhex,nk,Enk,rka,ne,Earray,JDOSinter,JDOSintraV, JDOSintraC )
!===============================================================================
! Joint density of states per carbon atom for an (n,m) carbon
! nanotube (states/carbon atom/eV)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  ne            number of energies
!  Earray(ne)    array of energies (eV)
! Output       :
!  JDOSinter(ne)       joint density of states (states/atom/eV)
!  JDOSintra(ne)       joint density of states (states/atom/eV)
!=======================================================================
  IMPLICIT NONE

! parameters
  REAL(8), PARAMETER     :: pi = 3.14159265358979D0

! input variables
  INTEGER, INTENT(in)    :: mu1, mu2
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  REAL(8), INTENT(in)    :: rka(nk)
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: Earray(ne)

! output variable
  REAL(8), INTENT(out)   :: JDOSinter(ne), JDOSintraC(ne), JDOSintraV(ne)

! working variables
  REAL(8), ALLOCATABLE   :: EnkInter(:,:)  !(nk,nhex)
  REAL(8), ALLOCATABLE   :: EnkIntraC(:,:) !(nk,2*(nhex-1))
  REAL(8), ALLOCATABLE   :: EnkIntraV(:,:) !(nk,2*(nhex-1))

  INTEGER                :: noutInter, noutIntra
  INTEGER                :: k, indexx1, indexx2, i, ie, nn, ii

  REAL(8)                :: T, trLength
  REAL(8)                :: rk, rkmin, rkmax, dk, fwhm, fwhm1, E, DSn


! allocate storage
  noutInter = 1
  noutIntra = 2
  ALLOCATE(EnkInter(nk,noutInter))
  ALLOCATE(EnkIntraC(nk,noutIntra))
  ALLOCATE(EnkIntraV(nk,noutIntra))

  EnkInter = 0.D0
  EnkIntraC = 0.D0
  EnkIntraV = 0.D0

! find FWHM linewidth based on energy array (eV)
  fwhm = ABS(Earray(2) - Earray(1))
  DO ie = 1, ne-1
     fwhm1 = ABS(Earray(ie) - Earray(ie+1))
     IF(fwhm1 < fwhm) fwhm = fwhm1
  END DO

! energy bands at evenly space k points (eV)
  DO k = 1, nk
    indexx1 = 1
    indexx2 = 1
    ! interband
    IF ( mu1 == mu2 ) THEN
        EnkInter(k,indexx1) = Enk(2,mu2,k) - Enk(1,mu1,k)
        indexx1 = indexx1 + 1
    ELSE
    ! intraband
        EnkIntraV(k,indexx2) = Enk(2,mu2,k) - Enk(2,mu1,k) ! valence
        EnkIntraC(k,indexx2) = Enk(1,mu2,k) - Enk(1,mu1,k) ! conduction
        indexx2 = indexx2 + 1
    END IF
  END DO

! accumulate density of states per unit length
  DO ie = 1, ne
     JDOSinter(ie) = 0.D0
     JDOSintraC(ie) = 0.D0
     JDOSintraV(ie) = 0.D0
     IF ( mu1 == mu2 ) THEN
         DO nn = 1, noutInter
            E = Earray(ie)
            CALL dos1Dgauss(noutInter,nk,rka,nk,EnkInter,E,fwhm,nn,DSn)
            JDOSinter(ie) = JDOSinter(ie) + DSn
         END DO
     END IF
     DO nn = 1, noutIntra
        E = Earray(ie)
        CALL dos1Dgauss(noutIntra,nk,rka,nk,EnkIntraC,E,fwhm,nn,DSn)
        JDOSintraC(ie) = JDOSintraC(ie) + DSn
     END DO
     DO nn = 1, noutIntra
        E = Earray(ie)
        CALL dos1Dgauss(noutIntra,nk,rka,nk,EnkIntraV,E,fwhm,nn,DSn)
        JDOSintraV(ie) = JDOSintraV(ie) + DSn
     END DO
  END DO

! convert to density of states/Carbon_atom/eV
  T = trLength(n,m)
  DO ie = 1, ne
     JDOSinter(ie) = (T/DBLE(noutInter)) * JDOSinter(ie)
     JDOSintraC(ie) = (T/DBLE(noutIntra)) * JDOSintraC(ie)
     JDOSintraV(ie) = (T/DBLE(noutIntra)) * JDOSintraV(ie)
  END DO

  DEALLOCATE(EnkInter, EnkIntraC, EnkIntraV)

END SUBROUTINE tubeJDOSij
!*******************************************************************************
!*******************************************************************************

