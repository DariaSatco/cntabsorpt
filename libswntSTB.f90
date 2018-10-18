!*******************************************************************************
!*******************************************************************************
! Project      : libswntSTB.f90
!===============================================================================
! Purpose      :
! Calculate electronic states for a SWNT using simple tight-binding model
!-------------------------------------------------------------------------------
! Authors      : ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
!                Daria Satco (dasha.shatco@gmail.com)
! Latest Vers. : 2018.10.18
!-------------------------------------------------------------------------------
! Reference(s) :
! [1] Physical Properties of Carbon Nanotubes
!     R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
! [2] Carbon Nanotube Photophysics, G. G. Samsonidze MIT Ph.D. Thesis (2006)
!-------------------------------------------------------------------------------
! Contents     :
! - SUBROUTINE tbTubeBand(n,m,mu,rk,Ek,Zk)
! - SUBROUTINE piHamTB(n,m,rk,mu,H,S)
! - SUBROUTINE stbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole)
! - SUBROUTINE stbDipolZ(n,m,n1,mu1,n2,mu2,rk,zDipole)
! - SUBROUTINE stbDipolZ2(n,m,n1,mu1,n2,mu2,rk,zDipole)
! - SUBROUTINE stbDipolXY(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
!*******************************************************************************
!*******************************************************************************
SUBROUTINE stbTubeBand(n,m,mu,rk,Ek,Zk)
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

END SUBROUTINE stbTubeBand
!*******************************************************************************
!*******************************************************************************
SUBROUTINE piHamTB(n,m,rk,mu,H,S)
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

! correct roundoff errors
  Htemp = H
  Stemp = S
  DO i = 1, 2
     DO j = 1, 2
        H(i,j) = (Htemp(i,j) + CONJG(Htemp(j,i))) / 2.D0
        S(i,j) = (Stemp(i,j) + CONJG(Stemp(j,i))) / 2.D0
     END DO
  END DO

END SUBROUTINE piHamTB
!*******************************************************************************
!*******************************************************************************
SUBROUTINE stbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole)
!===============================================================================
! optical dipole matrix element (1/Angstrom) for (n,m) carbon nanotube
!       D = < n1,mu1,rk | gradient | n2,mu2,rk > (1/Angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  n1            bra vector electronic state (n=1,2)
!  mu1           bra vector electronic state manifold (0...NHexagon-1)
!
!  n2            ket vector electronic state (n=1,2)
!  mu2           ket vector electronic state manifold (0...NHexagon-1)
!
!  rk            electronic state k (1/A) (-pi/T < k < pi/T)
! Output       :
!  cDipole(3)    complex dipole matrix element (1/Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, n1, mu1, n2, mu2
  REAL(8), INTENT(in)    :: rk

! output variable
  COMPLEX(8),INTENT(out) :: cDipole(3)

! working variables
  INTEGER                :: mmu1, mmu2
  REAL(8)                :: rkk
  COMPLEX(8)             :: xDipole, yDipole, zDipole

  cDipole = 0.D0

  CALL reducedCutLine(n,m,mu1,mmu1)
  CALL reducedCutLine(n,m,mu2,mmu2)

  rkk = rk
  CALL stbDipolXY(n,m,n1,mmu1,n2,mmu2,rkk,xDipole,yDipole)
  CALL stbDipolZ (n,m,n1,mmu1,n2,mmu2,rkk,zDipole)

  cDipole(1) = xDipole
  cDipole(2) = yDipole
  cDipole(3) = zDipole

  RETURN

END SUBROUTINE stbDipoleMX
!*******************************************************************************
!*******************************************************************************
SUBROUTINE stbDipolZ(n,m,n1,mu1,n2,mu2,rk,zDipole)
!===============================================================================
! this subroutine calls tbDipolZ2
! z component of optical dipole matrix element for (n,m) carbon nanotube
!        Dz = < n1,mu1,rk | d/dz | n2,mu2,rk > (1/Angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  n1            bra vector electronic state (n=1,2)
!  mu1           bra vector electronic state manifold (0...NHexagon-1)
!
!  n2            ket vector electronic state (n=1,2)
!  mu2           ket vector electronic state manifold (0...NHexagon-1)
!
!  rk            electronic state k (1/A) (-pi/T < k < pi/T)
! Output       :
!  zDipole       complex dipole matrix element (1/Angstroms)
!===============================================================================
  IMPLICIT NONE
  REAL(8), PARAMETER     :: rktol = .001D0

! input variables
  INTEGER, INTENT(in)    :: n, m, n1, mu1, n2, mu2
  REAL(8), INTENT(in)    :: rk

! output variable
  COMPLEX(8),INTENT(out) :: zDipole

! working variable
  REAL(8)                :: rd, rd1, rd2
  COMPLEX(8)             :: zD1, zD2

! selection rule
  IF (mu1 /= mu2) THEN
     zDipole = 0.D0
     RETURN
  END IF

! calculate z component of dipole vector (1/Angstroms)
  CALL stbDipolZ2(n,m,n1,mu1,n2,mu2,rk,zDipole)

! fix sign convention relative to neighboring points in k space
  CALL stbDipolZ2(n,m,n1,mu1,n2,mu2,rk-rktol,zD1)
  CALL stbDipolZ2(n,m,n1,mu1,n2,mu2,rk+rktol,zD2)

  rd  = REAL(zDipole)
  rd1 = REAL(zD1)
  rd2 = REAL(zD2)

  IF (rd < 0.D0 .AND. rd1 > 0.D0 .AND. rd2 > 0.D0) THEN
     zDipole = -zDipole
  END IF

  IF (rd > 0.D0 .AND. rd1 < 0.D0 .AND. rd2 < 0.D0) THEN
     zDipole = -zDipole
  END IF

! return zdipole, imaginary part is zero
  zDipole = CMPLX(REAL(zDipole), 0.D0)

END SUBROUTINE stbDipolZ
!*******************************************************************************
!*******************************************************************************
SUBROUTINE stbDipolZ2(n,m,n1,mu1,n2,mu2,rk,zDipole)
!===============================================================================
! this is the kernel of tbDipolZ subroutine
! z component of optical dipole matrix element for (n,m) carbon nanotube
!        Dz = < n1,mu1,rk | d/dz | n2,mu2,rk > (1/Angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  n1            bra vector electronic state (n=1,2)
!  mu1           bra vector electronic state manifold (0...NHexagon-1)
!
!  n2            ket vector electronic state (n=1,2)
!  mu2           ket vector electronic state manifold (0...NHexagon-1)
!
!  rk            electronic state k (1/A) (-pi/T < k < pi/T)
! Output       :
!  zDipole       complex dipole matrix element (1/Angstroms)
!===============================================================================
  IMPLICIT NONE
  REAL(8), PARAMETER     :: rktol = .0005D0

! input variables
  INTEGER, INTENT(in)    :: n, m, n1, mu1, n2, mu2
  REAL(8), INTENT(in)    :: rk

! output variable
  COMPLEX(8),INTENT(out) :: zDipole

! working variables
  INTEGER, SAVE, DIMENSION(0:4) :: nvecs = (/ 1, 3, 6, 3, 6 /)
  COMPLEX(8), SAVE       :: ci = (0.D0,1.D0)

  REAL(8)                :: Ek(2)   ! energy band
  COMPLEX(8)             :: Zk(2,2) ! wf. coeff.
  INTEGER                :: ier, iatom, jatom, nn, ivec, j1, j2, NNatom
  REAL(8)                :: rkk, phi, dipole
  COMPLEX(8)             :: c1, c2

! check input for errors
  ier = 0
  IF (n1 /= 1 .AND. n1 /= 2) THEN
     ier = 1
     WRITE (*,*) 'tbDipolZ2 err: invalid n1:', n1
  END IF

  IF (n2 /= 1 .AND. n2 /= 2) THEN
     ier = 1
     WRITE (*,*) 'tbDipolZ2 err: invalid n2:', n2
  END IF

  IF (ier /= 0) STOP

! selection rule
  IF (mu1 /= mu2) THEN
     zDipole = 0.D0
     RETURN
  END IF

! electronic pi orbitals (mu1,k)
  rkk = rk
  IF (ABS(rkk) < rktol) rkk = rktol
  CALL stbTubeBand(n,m,mu1,rkk,Ek,Zk)

! compute z component of dipole vector (1/Angstroms)
  zDipole = 0.D0
  DO iatom = 1, 2
        nn = 1
        DO ivec = 1, nvecs(nn)

           CALL NNj1j2(iatom,ivec,nn, j1,j2)
           CALL phaseFactor(n,m,j1,j2,rk,mu1,phi)

           jatom = NNatom(iatom,nn)
           CALL atomDipolZ(n,m,iatom,0,0,jatom,j1,j2,dipole)

           c1 = CONJG( Zk(iatom,n1) )
           c2 = Zk(jatom,n2) * CDEXP(ci*phi)

           zDipole = zDipole + c1*c2*dipole

        END DO

  END DO

  RETURN

END SUBROUTINE stbDipolZ2
!*******************************************************************************
!*******************************************************************************
SUBROUTINE stbDipolXY(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
!===============================================================================
! xy component of optical dipole matrix element for (n,m) carbon tube
!        Dx = < n1,mu1,rk | d/dx | n2,mu2,rk > (1/Angstroms)
!        Dy = < n1,mu1,rk | d/dy | n2,mu2,rk > (1/Angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  n1            bra vector electronic state (n=1,2)
!  mu1           bra vector electronic state manifold (0...NHexagon-1)
!
!  n2            ket vector electronic state (n=1,2)
!  mu2           ket vector electronic state manifold (0...NHexagon-1)
!
!  rk            electronic state k (1/A) (-pi/T < k < pi/T)
! Output       :
!  xDipole       x component of dipole matrix element (1/Angstroms)
!  yDipole       y component of dipole matrix element (1/Angstroms)
!===============================================================================
  IMPLICIT NONE
  REAL(8), PARAMETER     :: rktol = .0005D0

! input variables
  INTEGER, INTENT(in)    :: n, m, n1, mu1, n2, mu2
  REAL(8), INTENT(in)    :: rk

! output variables
  COMPLEX(8),INTENT(out) :: xDipole, yDipole

! working variables
  INTEGER, SAVE, DIMENSION(0:4) :: nvecs = (/ 1, 3, 6, 3, 6 /)
  COMPLEX(8), SAVE              :: ci = (0.D0,1.D0)

  REAL(8),    DIMENSION(2)      :: Ek1, Ek2  !(2)
  COMPLEX(8), DIMENSION(2,2)    :: Zk1, Zk2  !(2,2)

  COMPLEX(8)                    :: c1, c2
  REAL(8)                       :: dipole(3)

  INTEGER                :: nhex, nHexagon, ier, isel, iflip
  INTEGER                :: mmu1,mmu2,nn1,nn2
  INTEGER                :: iatom, jatom, nn, ivec, j1, j2, NNatom

  REAL(8)                :: rkk,phi

! check input for errors
  nhex=nHexagon(n,m)

  ier=0
  IF (n1 /= 1 .AND. n1 /= 2) THEN
     ier = 1
     WRITE (*,*) 'tbDipolXY err: invalid n1:', n1
  END IF

  IF (mu1 < 1 .OR. mu1 > nhex) THEN
     ier = 1
     WRITE (*,*) 'tbDipolXY err: invalid mu1:', mu1
  END IF

  IF (n2 /= 1 .AND. n2 /= 2) THEN
     ier = 1
     WRITE (*,*) 'tbDipolXY err: invalid n2:', n2
  END IF

  IF(mu2 < 1 .OR. mu2 > nhex) THEN
     ier = 1
     WRITE(*,*) 'tbDipolXY err: invalid mu2:', mu2
  END IF

  IF (ier /= 0) STOP

! initialize xDipole and yDipole
  xDipole = 0.D0
  yDipole = 0.D0

! selection rule
! Define states (nn1,mmu1) and (nn2,mmu2) so that mmu2 = mmu1+1
  isel  = 0
  iflip = 0

  IF (mu2 == mu1+1) THEN
     isel = 1
     mmu1 = mu1
     mmu2 = mu2
     nn1  = n1
     nn2  = n2
  END IF

  IF (mu2 == mu1-1) THEN
     isel  = 1
     iflip = 1 ! exchange bra and ket vectors
     mmu1  = mu2
     mmu2  = mu1
     nn1   = n2
     nn2   = n1
  END IF

  IF (mu1 == nhex .AND. mu2 == 1) THEN
     isel = 1
     mmu1 = mu1
     mmu2 = mu2
     nn1  = n1
     nn2  = n2
  END IF

  IF (mu1 == 1 .AND. mu2 == nhex) THEN
     isel  = 1
     iflip = 1 ! exchange bra and ket vectors
     mmu1  = mu2
     mmu2  = mu1
     nn1   = n2
     nn2   = n1
  END IF

  IF (isel /= 1) RETURN

! mmu2 = mmu1+1 case
! electronic pi orbitals (mu1,k) and (mu2,k)
  rkk = rk
  IF (ABS(rkk) < rktol) rkk = rktol
  CALL stbTubeBand(n,m,mmu1,rkk,Ek1,Zk1)
  CALL stbTubeBand(n,m,mmu2,rkk,Ek2,Zk2)

! compute x and y components of dipole vector (1/Angstroms)
  DO iatom = 1, 2
        nn = 1
        DO ivec = 1, nvecs(nn)

           CALL NNj1j2(iatom,ivec,nn,j1,j2)
           CALL phaseFactor(n,m,j1,j2,-rkk,-mmu1,phi)

           jatom = NNatom(iatom,nn)
           CALL atomDipoleMX(n,m,jatom,j1,j2,iatom,0,0,dipole)

           c1 = CONJG( Zk1(jatom,nn1) )
           c2 = Zk2(iatom,nn2)*CDEXP(ci*phi)

           xDipole = xDipole+c1*c2*(dipole(1)-ci*dipole(2))

        END DO

  END DO

  xDipole = xDipole/2.D0
  yDipole = ci*xDipole

! use symmetry relation to reverse the exchange bra and ket vectors
  IF (iflip == 1) THEN
     xDipole = -CONJG(xDipole)
     yDipole = -CONJG(yDipole)
  END IF

  RETURN
END SUBROUTINE stbDipolXY
!*******************************************************************************
!*******************************************************************************
