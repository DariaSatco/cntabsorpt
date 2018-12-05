!*******************************************************************************
!*******************************************************************************
! Project      : libDrudeOpt.f90
!===============================================================================
! Purpose      :
! Optical properties of a SWNT, i.e. Opt. Matrix Elements and Opt. Absorption
!-------------------------------------------------------------------------------
! Authors      :
! - ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
! - Daria Satco  (dasha.shatco@gmail.com)
! Latest Vers. : 2018.11.30
!-------------------------------------------------------------------------------
! Reference(s) :
! [1] Physical Properties of Carbon Nanotubes
!     R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
! [2] J.-W, Jiang et al, PRB, 73, 235434 (2006)
! [3] J. Jiang et al, PRB, 71, 205420 (2005).
! [4] A. Gruneis, PhD thesis, Tohoku University (2004).
! [5] Lin, M. F., and Kenneth W-K. Shung., PRB, 50, 17744 (1994).
! [6] Sasaki, Ken-ichi, and Yasuhiro Tokura, Physical Review Applied 9.3 034018 (2018).
!-------------------------------------------------------------------------------
! Contents     :
! - SUBROUTINE tbDipoleMXDr(n,m,n1,mu1,n2,mu2,rk,cDipole)
! - SUBROUTINE tbDipolZDr(n,m,n1,mu1,n2,mu2,rk,zDipole)
! - SUBROUTINE tbDipolZ2Dr(n,m,n1,mu1,n2,mu2,rk,zDipole)
! - SUBROUTINE tbDipolXYDr(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
! - SUBROUTINE DielPermittivityDr(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,ebg,fwhm,ne,hw,eps1,eps2)
! - SUBROUTINE DynConductivityDr(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm1,sigm2)
!*******************************************************************************
!*******************************************************************************
SUBROUTINE DielPermittivityDr(n,m,nhex,nk,rka,Enk,Tempr,Efermi,epol,ebg,fwhm,ne,hw,eps1,eps2)
!===============================================================================
! Compute the real and imaginary parts of the DRUDE dielectric function as a function
! of probe photon energy
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  nhex          number of hexagons
!  nk            number of k points
!  rka           array of k points (1/A)
!  Enk           array of energies (eV)
!  cDipole       array of complex matrix elements (1/A)
!  Tempr         lattice temperature (deg K)
!  Efermi        Fermi level
!  epol(3)       complex unit electric polarization vector (none)
!  ebg           background permittivity (dimensionless)
!  fwhm          fwhm probe linewidth (eV)
!  ne            number of probe photon energies
!  hw(ne)        array of probe photon energies (eV)
! Output       :
!  eps1(ne)      real part of dielectric function (none)
!  eps2(ne)      imaginary part of dielectric function (none)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in)    :: epol(3)
  REAL(8), INTENT(in)    :: ebg

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)


! output variable
  REAL(8),  INTENT(out)  :: eps1(ne), eps2(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre
  REAL(8)                :: eps0(ne)
  REAL(8)                :: reint(ne), imagint(ne)
  REAL(8)                :: fnk(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ress(:), imss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  INTEGER                :: ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: diameter, tubeDiam, area, eps1Part, eps2Part

!background dielectric constant (dimensionless)
eps0(:) = ebg

! dielectric function prefactor (eV**3 Angstroms**3)
! --------------------------------------------------
! see for the prefactor expression paper:
! Sanders, G. D., et al.
! "Resonant coherent phonon spectroscopy of single-walled carbon nanotubes."
! Physical Review B 79.20 (2009): 205434.
  diameter = tubeDiam(n,m)        !(Angstroms)
  area = pi*(diameter/2.D0)**2    !(Angstroms**2)
  pre  = 8.D0*pi*e2*hbarm**2/area !(eV**3 Angstroms**3)

! calculate Fremi distributions
  CALL FermiDistributionArray(nhex,nk,Enk,Tempr,Efermi,fnk)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ress) .EQV. .TRUE.) DEALLOCATE(ress)
  IF (ALLOCATED(imss) .EQV. .TRUE.) DEALLOCATE(imss)
  ALLOCATE(ress(ne))
  ALLOCATE(imss(ne))
  ress = 0.D0
  imss = 0.D0

  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction
     DO mu1 = 1, nhex
        n2 = n1
        mu2 = mu1

        CALL RealImagPartIntegralDr(n,m,n1,mu1,nhex,nk,rka,Enk,fnk,epol,fwhm,ne,hw,reint,imagint)
! accumulate dielectric function vs photon energy
        DO ie = 1, ne
            IF (hw(ie) .le. ptol) THEN
                eps1Part = imagint(ie)/1.D-3
                eps2Part = reint(ie)/1.D-3
            ELSE
                eps1Part = imagint(ie)/hw(ie)  ! (1/Angstroms**3 1/eV**3)
                eps2Part = reint(ie)/hw(ie)    ! (1/Angstroms**3 1/eV**3)
            END IF

            ress(ie) = ress(ie) + eps1Part
            imss(ie) = imss(ie) + eps2Part
        END DO

     END DO
  END DO

! real and imaginary part of dielectric function (dimensionless)
  eps1 = eps0 + pre*ress
  eps2 = pre*imss

END SUBROUTINE DielPermittivityDr
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipolZDr(n,m,n1,mu1,rk1,rk2,zDipole)
!===============================================================================
! this subroutine calls tbDipolZ2
! z component of optical dipole matrix element for (n,m) carbon nanotube
!        Dz = < n,mu,rk1 | d/dz | n,mu,rk2 > (1/Angstroms)
! initial and final state differ only by rk, n and mu are conserved
! these matrix elements are designed to consider Drude conductivity
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  n1            1 or 2 correspond to valence or conduction band
!  mu1           cutting line number
!
!  rk1           bra electronic state k (1/A) (-pi/T < k < pi/T)
!  rk2           ket electronic state k (1/A) (-pi/T < k < pi/T)
! Output       :
!  zDipole       complex dipole matrix element (1/Angstroms)
!===============================================================================
  IMPLICIT NONE
  REAL(8), PARAMETER     :: rktol = .001D0

! input variables
  INTEGER, INTENT(in)    :: n, m, n1, mu1
  REAL(8), INTENT(in)    :: rk1,rk2

! output variable
  COMPLEX(8),INTENT(out) :: zDipole

! working variable
  REAL(8)                :: rd, rd1, rd2
  COMPLEX(8)             :: zD1, zD2

  IF (rk1 == rk2) THEN
     zDipole = 0.D0
     RETURN
  END IF

! calculate z component of dipole vector (1/Angstroms)
  CALL tbDipolZ2Dr(n,m,n1,mu1,n1,mu1,rk1,rk2,zDipole)

! fix sign convention relative to neighboring points in k space
  CALL tbDipolZ2Dr(n,m,n1,mu1,n1,mu1,rk1-rktol,rk2-rktol,zD1)
  CALL tbDipolZ2Dr(n,m,n1,mu1,n1,mu1,rk1+rktol,rk2+rktol,zD2)

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

END SUBROUTINE tbDipolZDr
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipolZ2Dr(n,m,n1,mu1,n2,mu2,rk1,rk2,zDipole)
!===============================================================================
! this is the kernel of tbDipolZDr subroutine
! z component of optical dipole matrix element for (n,m) carbon nanotube
!        Dz = < n1,mu1,rk1 | d/dz | n2,mu2,rk2 > (1/Angstroms)
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  n1            bra vector electronic state (n=1,2)
!  mu1           bra vector electronic state manifold (0...NHexagon-1)
!
!  n2            ket vector electronic state (n=1,2)
!  mu2           ket vector electronic state manifold (0...NHexagon-1)
!
!  rk1           bra electronic state k (1/A) (-pi/T < k < pi/T)
!  rk2           ket electronic state k (1/A) (-pi/T < k < pi/T)
! Output       :
!  zDipole       complex dipole matrix element (1/Angstroms)
!===============================================================================
  IMPLICIT NONE
  REAL(8), PARAMETER     :: rktol = .0005D0

! input variables
  INTEGER, INTENT(in)    :: n, m, n1, mu1, n2, mu2
  REAL(8), INTENT(in)    :: rk1, rk2

! output variable
  COMPLEX(8),INTENT(out) :: zDipole

! working variables
  INTEGER, SAVE, DIMENSION(0:4) :: nvecs = (/ 1, 3, 6, 3, 6 /)
  COMPLEX(8), SAVE       :: ci = (0.D0,1.D0)

  REAL(8)                :: Ek1(2),Ek2(2)   ! energy band
  COMPLEX(8)             :: Zk1(2,2),Zk2(2,2) ! wf. coeff.
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

  IF (n1 /= n2 .and. rk1 /= rk2) THEN
     zDipole = 0.D0
     RETURN
  END IF

! electronic pi orbitals (mu1,k1)
  rkk = rk1
  IF (ABS(rkk) < rktol) rkk = rktol
  CALL etbTubeBand(n,m,mu1,rkk,Ek1,Zk1)

! electronic pi orbitals (mu2,k2)
  rkk = rk2
  IF (ABS(rkk) < rktol) rkk = rktol
  CALL etbTubeBand(n,m,mu2,rkk,Ek2,Zk2)

! compute z component of dipole vector (1/Angstroms)
  zDipole = 0.D0
  DO iatom = 1, 2
     DO nn = 1, 4
        DO ivec = 1, nvecs(nn)

           CALL NNj1j2(iatom,ivec,nn,j1,j2)
           CALL phaseFactor(n,m,j1,j2,-rk1,-mu1,phi)

           jatom = NNatom(iatom,nn)
           CALL atomDipolZ(n,m,jatom,j1,j2,iatom,0,0,dipole)

           c1 = CONJG( Zk1(jatom,n1) )
           c2 = Zk2(iatom,n2) * CDEXP(ci*phi)

           zDipole = zDipole + c1*c2*dipole

        END DO
     END DO
  END DO

  RETURN

END SUBROUTINE tbDipolZ2Dr
!*******************************************************************************
!*******************************************************************************
SUBROUTINE DynConductivityDr(n,m,nhex,nk,rka,Enk,Tempr,Efermi,epol,fwhm,ne,hw,sigm1,sigm2)
!===============================================================================
! Compute the real and imaginary parts of the DRUDE dynamical conductivity as a function
! of probe photon energy
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  nhex          number of hexagons
!  nk            number of k points
!  rka           array of k points (1/A)
!  Enk           array of energies (eV)
!  cDipole       array of complex matrix elements (1/A)
!  Tempr         lattice temperature (deg K)
!  Efermi        Fermi level
!  epol(3)       complex unit electric polarization vector (none)
!  ebg           background permittivity (dimensionless)
!  fwhm          fwhm probe linewidth (eV)
!  ne            number of probe photon energies
!  hw(ne)        array of probe photon energies (eV)
! Output       :
!  sigm1(ne)      real part of dynamical conductivity (e^2/h)
!  sigm2(ne)      imaginary part of dynamical conductivity (e^2/h)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in)    :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: sigm1(ne), sigm2(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), ALLOCATABLE   :: ress(:), imss(:) !(ne)
  REAL(8)                :: reint(ne), imagint(ne)
  REAL(8)                :: fnk(2,nhex,nk)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: h = 4.135D-15         !(eV-s)

  INTEGER                :: ie
  INTEGER                :: n1, mu1, n2, mu2

  REAL(8)                :: diameter, tubeDiam


! conductivity function prefactor (eV**2 Angstroms**3) * [e^2/h] = (A/s)
  diameter = tubeDiam(n,m)              !(Angstroms)
  pre  = 32.D0*hbarm**2/diameter * e2/h !(eV**2 Angstroms**3) * [e^2/h] = (A/s)

  IF (ALLOCATED(ress) .EQV. .TRUE.) DEALLOCATE(ress)
  IF (ALLOCATED(imss) .EQV. .TRUE.) DEALLOCATE(imss)
  ALLOCATE(ress(ne))
  ALLOCATE(imss(ne))
  ress = 0.D0
  imss = 0.D0

! calculate Fremi distributions
  CALL FermiDistributionArray(nhex,nk,Enk,Tempr,Efermi,fnk)

! sum over n1, mu1
  DO n1 = 1, 2         ! 1 <-> valence, 2 <-> conduction
     DO mu1 = 1, nhex
        n2 = n1
        mu2 = mu1

        CALL RealImagPartIntegralDr(n,m,n1,mu1,nhex,nk,rka,Enk,fnk,epol,fwhm,ne,hw,reint,imagint)

        DO ie = 1, ne
           ress(ie) = ress(ie) + reint(ie)
           imss(ie) = imss(ie) + imagint(ie)
        END DO


     END DO
  END DO

! real and imaginary parts of conductivity [e^2/h] = (A/s)
  sigm1 = pre*ress
  sigm2 = -pre*imss

END SUBROUTINE DynConductivityDr
!*******************************************************************************
!*******************************************************************************
SUBROUTINE RealImagPartIntegralDr(n,m,n1,mu1,nhex,nk,rka,Enk,fnk,epol,fwhm,ne,hw,reint,imagint)
!===============================================================================
! Compute the integral over dk1,dk2 for particular n, mu which corresponds to
! the Drude dynamical conductivity as a function of probe photon energy
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral indeces
!  n1            1 or 2 corresponds to valence or conduction band
!  mu1           cutting line number
!  nhex          number of hexagons
!  nk            number of k points
!  rka           array of k points (1/A)
!  Enk           array of energies (eV)
!  fnk           array of fermi distributions (dimensionless)
!  epol(3)       complex unit electric polarization vector (none)
!  ebg           background permittivity (dimensionless)
!  fwhm          fwhm probe linewidth (eV)
!  ne            number of probe photon energies
!  hw(ne)        array of probe photon energies (eV)
! Output       :
!  reint(ne)      real part of integral over dk (1/Angstroms**3 1/eV**2)
!  imagint(ne)    imaginary part of integral over dk (1/Angstroms**3 1/eV**2)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: n1, mu1
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  REAL(8), INTENT(in)    :: fnk(2,nhex,nk)

  REAL(8), INTENT(in)    :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: reint(ne), imagint(ne)

! working variables and parameter
  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: ptol  =  1.D-15
  REAL(8), PARAMETER     :: hbarvfermi = 6.582119 !(eV-A) !hbar*vfermi, vfermi = 10**6 m/s

  INTEGER                :: k1, k2, ie

! for calling some functions
  REAL(8)                :: dk, gammaDr, drudeBroadening, dkDr, rk1, rk2
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvgRe, diracAvgIm, Eab, diracDelta, diracDelta_im

  COMPLEX(8)             :: zDipole
  COMPLEX(8)             :: css

 reint = 0.D0
 imagint = 0.D0

 gammaDr = 0.05 ! eV
 dkDr = gammaDr/hbarvfermi

! sum over k for particular n1, mu1, n2
DO k1 = 1,nk
    DO k2 = 1,nk
! define Drude broadening distribution
    rk1 = rka(k1)
    rk2 = rka(k2)
    IF (ABS(rk1-rk2) < 5*dkDr) THEN
        !drudeBroadening = 1.D0/dkDr
        drudeBroadening = dkDr/((rk1 - rk2)**2 + dkDr**2) * 1.D0/pi
    ELSE
        CYCLE
    END IF
!energy difference
    Eab = Enk(n1,mu1,k2) - Enk(n1,mu1,k1)
! squared optical dipole matrix element (1/Angstroms**2)
    CALL tbDipolZDr(n,m,n1,mu1,rk1,rk2,zDipole)
    css = epol(3)*zDipole
! square of matrix element
    p2   = CDABS(css)**2
! multiply by distribuion function
    p2df = p2*(fnk(n1,mu1,k1) - fnk(n1,mu1,k2)) !(1/Angstroms**2)

! if small then skip
    IF ( ABS(Eab) .lt. ptol ) THEN
        IF (ABS(p2df) > ptol) STOP 'WARNING: danger of division by zero'
    ENDIF

    IF (ABS(p2df) <= ptol) CYCLE

! k-cell boundaries (1/A)
    IF (k1 == 1) THEN
        x1 = rka(1)
    ELSE
        x1 = (rka(k1-1) + rka(k1))/2.D0
    END IF

    IF (k2 == nk) THEN
        x2 = rka(nk)
    ELSE
        x2 = (rka(k2+1) + rka(k2))/2.D0
    END IF

        dk = x2-x1

! band energies at k-cell boundaries (eV)
    IF (k1 == 1) THEN
        enk1n = Enk(n1,mu1,1)
        enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
    ELSE IF (k1 == nk) THEN
        enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
        enk1p = Enk(n1,mu1,nk)
    ELSE
        enk1n = (Enk(n1,mu1,k1-1)+Enk(n1,mu1,k1  ))/2.D0
        enk1p = (Enk(n1,mu1,k1  )+Enk(n1,mu1,k1+1))/2.D0
    END IF

     IF (k2 == 1) THEN
        enk2n = Enk(n1,mu1,1)
        enk2p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
     ELSE IF (k2 == nk) THEN
        enk2n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
        enk2p = Enk(n1,mu1,nk)
     ELSE
        enk2n = (Enk(n1,mu1,k2-1)+Enk(n1,mu1,k2  ))/2.D0
        enk2p = (Enk(n1,mu1,k2 )+Enk(n1,mu1,k2+1))/2.D0
     END IF

! accumulate function vs photon energy
    DO ie = 1, ne
        y1 = enk2n - enk1n - hw(ie)
        y2 = enk2p - enk1p - hw(ie)

        ! calculate real part of integral
        diracAvgRe = diracDelta(x1,y1,x2,y2,fwhm)    !(1/eV)
        IF(diracAvgRe == 0.) THEN
            reint(ie) = reint(ie) + 0.D0
        ELSE
            reint(ie) = reint(ie) + dk/2*p2df*diracAvgRe/Eab * dk *drudeBroadening   ! (1/Angstroms**3 1/eV**2)
        END IF

        ! calculate imaginary part of integral
        diracAvgIm = diracDelta_im(x1,y1,x2,y2,fwhm)    !(1/eV)
        IF(diracAvgIm == 0.) THEN
            imagint(ie) = imagint(ie) + 0.D0
        ELSE
            imagint(ie) = imagint(ie) + dk/(2*pi)*p2df*diracAvgIm/Eab * dk *drudeBroadening    ! (1/Angstroms**3 1/eV**2)
        END IF

    END DO
  END DO
END DO

END SUBROUTINE RealImagPartIntegralDr
!*******************************************************************************
!*******************************************************************************
