!*******************************************************************************
!*******************************************************************************
! Project      : libswntOpt.f90
!===============================================================================
! Purpose      :
! Optical properties of a SWNT, i.e. Opt. Matrix Elements and Opt. Absorption
!-------------------------------------------------------------------------------
! Authors      :
! - ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
! - Daria Satco  (dasha.shatco@gmail.com)
! Latest Vers. : 2018.09.30
!-------------------------------------------------------------------------------
! Reference(s) :
! [1] Physical Properties of Carbon Nanotubes
!     R. Saito, G. Dresselhaus, M. S. Dresselhaus (ICP, 1998)
! [2] J.-W, Jiang et al, PRB, 73, 235434 (2006)
! [3] J. Jiang et al, PRB, 71, 205420 (2005).
! [4] A. Gruneis, PhD thesis, Tohoku University (2004).
! [5] Lin, M. F., and Kenneth W-K. Shung., PRB, 50, 17744 (1994).
!-------------------------------------------------------------------------------
! Contents     :
! - SUBROUTINE polVector(theta,epol)
! - SUBROUTINE imagDielAlpha(ne,hw,eps2,refrac,alpha)
! - SUBROUTINE imagDielEn(n,m,Tempr,doping,epol,fwhm,ne,hw,eps2)
! - SUBROUTINE tbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole)
! - SUBROUTINE tbDipolZ(n,m,n1,mu1,n2,mu2,rk,zDipole)
! - SUBROUTINE tbDipolZ2(n,m,n1,mu1,n2,mu2,rk,zDipole)
! - SUBROUTINE tbDipolXY(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
! - SUBROUTINE atomDipolZ(n,m,iiatom,jj1,jj2,iatom,j1,j2,zdipole)
! - SUBROUTINE atomDipoleMX(n,m,iiatom,jj1,jj2,iatom,j1,j2,dipole)
! - FUNCTION gx (cg, cg1, cg3, cg5, cg7, cg9, cg11, cg13)
! - FUNCTION gy (cg, cg1, cg3, cg5, cg7, cg9, cg11, cg13)
! - FUNCTION gz (cg, cg1, cg3, cg5, cg7, cg9, cg11, cg13)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Daria Satco added (autumn 2018) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! - SUBROUTINE realDielEn(n,m,Tempr,doping,epol,fwhm,ne,hw,eps1)
! - SUBROUTINE realDielEn_met2(ne,hw,eps2,eps1)
! - SUBROUTINE imagDielEn_met2(ne,hw,eps1,eps2)
! - SUBROUTINE EELS(ne,eps1,eps2,eelspec)
! - SUBROUTINE ImagDynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm2)
! - SUBROUTINE RealDynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm1)
! - SUBROUTINE ImagDynConductivityInter(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm2)
! - SUBROUTINE ImagDynConductivityIntra(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm2)
! - SUBROUTINE Func_test(nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,difFermiDist,matrElementSq,diracAvgFunc)
! - SUBROUTINE tbDipolXYOrt(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
!*******************************************************************************
!*******************************************************************************
SUBROUTINE polVector(theta,epol)
!===============================================================================
! unit polarization vector (dimensionless)
!-------------------------------------------------------------------------------
! Input        :
!  theta         angle between tube axis and electric field (degree)
! Output       :
!  epol(3)       complex unit electric polarization vector (dimensionless) 
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  
! input variable
  REAL(8), INTENT(in)    :: theta

! output variable
  REAL(8),INTENT(out) :: epol(3)
      
! working variable
  REAL(8)                :: angle
      
! linear polarization vector
  IF (theta == 0.) THEN
     epol(1) = 0.D0
     epol(2) = 0.D0
     epol(3) = 1.D0
  ELSE IF (theta == 90.) THEN
     epol(1) = 1.D0
     epol(2) = 0.D0
     epol(3) = 0.D0    
  ELSE
     angle = DBLE(theta) / 57.29577951308233D0
     epol(1) = DSIN(angle)
     epol(2) = 0.D0
     epol(3) = DCOS(angle)
  END IF

  RETURN

END SUBROUTINE polVector
!*******************************************************************************
!*******************************************************************************
SUBROUTINE imagDielAlpha(ne,hw,eps2,refrac,alpha)
!===============================================================================
! Calculate absorption coefficient (1/cm) from imaginary dielectric function
!-------------------------------------------------------------------------------
! Input        :
!  ne            number of probe photon energies
!  hw(ne)        array of probe photon energies (eV)
!  eps2(ne)      imaginary part of dielectric function (none)  
!  refrac        refractive index (dimensionless)
! Output       :
!  alpha(ne)     absorption coefficient (1/cm)
!===============================================================================
  IMPLICIT NONE
  
! input variables
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne), eps2(ne)
  REAL(8), INTENT(in)    :: refrac

! output variable
  REAL(8), INTENT(out)   :: alpha(ne)

! working variable and parameter
  REAL(8), PARAMETER     :: hbarc = 1.97D-5 !(eV cm)
  INTEGER                :: ie
      
  DO ie = 1, ne
     alpha(ie) = eps2(ie)*hw(ie) / (refrac*hbarc)
  END DO

END SUBROUTINE imagDielAlpha
!*******************************************************************************
!*******************************************************************************
SUBROUTINE imagDielEn(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,eps2)
!===============================================================================
! Compute the imaginary part of the dielectric function as a function
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
!  eps2(ne)      imaginary part of dielectric function (none)           
!===============================================================================
  IMPLICIT NONE
      
! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in)    :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

  
! output variable
  REAL(8),  INTENT(out)  :: eps2(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: ptol  =  1.D-15
      
  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: diameter, tubeDiam, area
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta
  
  COMPLEX(8)             :: css

! first pass
! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! dielectric function prefactor (eV**3 Angstroms**3)
! --------------------------------------------------
! see for the prefactor expression paper:
! Sanders, G. D., et al.
! "Resonant coherent phonon spectroscopy of single-walled carbon nanotubes."
! Physical Review B 79.20 (2009): 205434.

     diameter = tubeDiam(n,m)        !(Angstroms)
     area = pi*(diameter/2.D0)**2    !(Angstroms**2)
     pre  = 8.D0*pi*e2*hbarm**2/area !(eV**3 Angstroms**3)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss=0.D0
      
  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction
     DO mu1 = 1, nhex              
        DO n2 = 1, 2
           DO mu2 = 1, nhex
              
              IF (n1 == n2 .AND. mu1 == mu2) CYCLE
              
              DO k=1,nk

!energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k))  !(1/Angstroms**2)

                 ! if small then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

                 IF (ABS(p2df) <= ptol) CYCLE        
                 
! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1
                 
! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk) 
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! accumulate dielectric function vs photon energy
                 DO ie = 1, ne        
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)
                    diracAvg = diracDelta(x1,y1,x2,y2,fwhm)  ! (1/eV)
                    IF (diracAvg == 0.) CYCLE

                    IF (hw(ie) .le. ptol) THEN
                        ss(ie) = ss(ie) + dk/2 * p2df * diracAvg / (1.D-3 * Eab)
                    ELSE
                        ss(ie) = ss(ie) + dk/2 * p2df * diracAvg / (hw(ie) * Eab)  ! (1/Angstroms**3 1/eV**3)
                    END IF
                 END DO
                 
              END DO
           END DO
        END DO
     END DO
  END DO
  
! imaginary part of dielectric function (dimensionless)      
  eps2 = pre*ss                
  
END SUBROUTINE imagDielEn
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole)
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
  CALL tbDipolXY(n,m,n1,mmu1,n2,mmu2,rkk,xDipole,yDipole)
  CALL tbDipolZ (n,m,n1,mmu1,n2,mmu2,rkk,zDipole)      
      
  cDipole(1) = xDipole
  cDipole(2) = yDipole
  cDipole(3) = zDipole
  
  RETURN

END SUBROUTINE tbDipoleMX
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipolZ(n,m,n1,mu1,n2,mu2,rk,zDipole)
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
  CALL tbDipolZ2(n,m,n1,mu1,n2,mu2,rk,zDipole)

! fix sign convention relative to neighboring points in k space
  CALL tbDipolZ2(n,m,n1,mu1,n2,mu2,rk-rktol,zD1)
  CALL tbDipolZ2(n,m,n1,mu1,n2,mu2,rk+rktol,zD2)
      
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
      
END SUBROUTINE tbDipolZ
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipolZ2(n,m,n1,mu1,n2,mu2,rk,zDipole)
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
  CALL etbTubeBand(n,m,mu1,rkk,Ek,Zk) 
      
! compute z component of dipole vector (1/Angstroms)      
  zDipole = 0.D0
  DO iatom = 1, 2
     DO nn = 1, 4
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
  END DO
      
  RETURN

END SUBROUTINE tbDipolZ2
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipolXY(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
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
  CALL etbTubeBand(n,m,mmu1,rkk,Ek1,Zk1)
  CALL etbTubeBand(n,m,mmu2,rkk,Ek2,Zk2)
        
! compute x and y components of dipole vector (1/Angstroms)
  DO iatom = 1, 2
     DO nn = 1, 4
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
  END DO
      
  xDipole = xDipole/2.D0
  yDipole = ci*xDipole
      
! use symmetry relation to reverse the exchange bra and ket vectors
  IF (iflip == 1) THEN
     xDipole = -CONJG(xDipole)
     yDipole = -CONJG(yDipole)
  END IF
                              
  RETURN
END SUBROUTINE tbDipolXY
!*******************************************************************************
!*******************************************************************************
SUBROUTINE atomDipolZ(n,m,iiatom,jj1,jj2,iatom,j1,j2,zdipole)
!===============================================================================
! z component of dipole vector <r'J'| d/dz |r J> (1/Angstroms)
! between two pi orbitals centered at R_{r'J'} and R_{rJ} on the
! surface of an (n,m) carbon nanotube.
!
! References:
! J. Jiang et al, PRB, 71, 205420 (2005).
! A. Gruneis, PhD thesis, Tohoku University (2004).
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iiatom        first atom in the unit cell (1=A,2=B)
!  jj1,jj2       center of first two atom unit cell in a1 a2 basis
!  iatom         second atom in the unit cell (1=A,2=B)
!  j1,j2         center of second two atom unit cell in a1 a2 basis
! output
!  zdipole       z component of atomic dipole vector (1/Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m, iiatom, jj1, jj2, iatom, j1, j2

! output variable
  REAL(8), INTENT(out)   :: zdipole

! gaussian expansion parameters for carbon 2p_{pi} orbitals
! Reference: J. Jiang, PRB, 71, 205420 (2005)
  REAL(8), DIMENSION(4)  :: Im = (/  .050, .413 ,  1.061, 1.046  /)
  REAL(8), DIMENSION(4)  :: a  = (/ .1067, .6078, 29.586, 3.338  /)

  REAL(8)                :: R1(3), R2(3)
      
  INTEGER                :: ier, m1, m2
  REAL(8)                :: R1x, R1y, R1z, R2x, R2y, R2z
  REAL(8)                :: ss, pre, a1, a2, gz
      
! check input for errors
  ier = 0
  IF (iiatom < 1 .OR. iiatom > 2) THEN
     ier = 1
     WRITE (*,*) 'atomDipolZ err: invalid iiatom:', iiatom
  END IF
  
  IF (iatom  < 1 .OR. iatom  > 2) THEN
     ier = 1
     WRITE (*,*) 'atomDipolZ err: invalid iatom:', iatom
  END IF
  
  IF (ier.NE.0) STOP
      
! xyz coordinates of first and second atoms (Angstroms)
  CALL rxyzJ1J2Vec(n,m, iiatom,jj1,jj2, R1)  
  CALL rxyzJ1J2Vec(n,m,  iatom, j1, j2, R2)
            
! convert xyz coordinates from Angstroms to Bohr radii
  R1  = R1/.5292D0
  R2  = R2/.5292D0

  R1x = R1(1)
  R1y = R1(2)
  R1z = R1(3)
  
  R2x = R2(1)
  R2y = R2(2)
  R2z = R2(3)
  
! z component of atomic diple vector (1/Bohr radii)
  ss = 0.D0
  DO m1 = 1,4
     DO m2 = 1,4
        pre = Im(m1)*Im(m2)
        a1 = a(m1)
        a2 = a(m2)
        
        ss = ss + pre*gz(R1x,R1y,R1z,R2x,R2y,R2z,a1,a2)
                        
     END DO
  END DO
  
  zdipole = ss
      
! convert dipole vector from (1/Bohr radii) to (1/Angstroms)
  zdipole = zdipole/.5292D0
      
END SUBROUTINE atomDipolZ
!*******************************************************************************
!*******************************************************************************
SUBROUTINE atomDipoleMX(n,m,iiatom,jj1,jj2,iatom,j1,j2,dipole)
!===============================================================================
! compute the dipole vector <r'J'| gradient |r J> (1/Angstroms)
! between two pi orbitals centered at R_{r'J'} and R_{rJ} on the
! surface of an (n,m) carbon nanotube.
!
! References:
! J. Jiang et al, PRB, 71, 205420 (2005).
! A. Gruneis, PhD thesis, Tohoku University (2004).
!-------------------------------------------------------------------------------
! Input        :
!  n,m           chiral vector coordinates in (a1,a2)
!  iiatom        first atom in the unit cell (1=A,2=B)
!  jj1,jj2       center of first two atom unit cell in (a1,a2)
!  iatom         second atom in the unit cell (1=A,2=B)
!  j1,j2         center of second two atom unit cell in (a1,a2)
! Output
!  dipole(3)     atomic dipole vector (1/Angstroms)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: iiatom, jj1, jj2, iatom, j1, j2

! output variable
  REAL(8), INTENT(out)   :: dipole(3)

! gaussian expansion parameters for carbon 2p_{pi} orbitals
! Reference: J. Jiang, PRB, 71, 205420 (2005). 
  REAL(8), DIMENSION(4)  :: Im = (/  .050, .413 ,  1.061, 1.046  /)
  REAL(8), DIMENSION(4)  :: a  = (/ .1067, .6078, 29.586, 3.338  /)
  
  REAL(8)                :: R1(3), R2(3), ss(3)
  
  INTEGER                :: ier, m1, m2

  REAL(8)                :: gx,gy,gz      
  REAL(8)                :: R1x, R1y, R1z, R2x, R2y, R2z
  REAL(8)                :: pre,a1,a2
      
! check input for errors
  ier = 0
  IF (iiatom.LT.1 .OR. iiatom.GT.2) THEN
     ier = 1
     WRITE (*,*) 'atomic dipole err: invalid iiatom:', iiatom
  END IF
  
  IF (iatom < 1 .OR. iatom > 2) THEN
     ier = 1
     WRITE (*,*) 'atomic dipole err: invalid iiatom:', iatom
  END IF
  
  IF(ier.NE.0) STOP
      
! xyz coordinates of first and second atoms (Angstroms)
  CALL rxyzJ1J2Vec(n,m,iiatom,jj1,jj2, R1)  
  CALL rxyzJ1J2Vec(n,m, iatom, j1, j2, R2)
            
! convert xyz coordinates from Angstroms to Bohr radii
  R1  = R1/.5292D0
  R2  = R2/.5292D0

  R1x = R1(1)
  R1y = R1(2)
  R1z = R1(3)
      
  R2x = R2(1)
  R2y = R2(2)
  R2z = R2(3)
      
! atomic diple vector (1/Bohr radii)
  ss = 0.D0
  DO m1 = 1, 4
     DO m2 = 1, 4
        pre = Im(m1) * Im(m2)
        a1 = a(m1)
        a2 = a(m2)
        
        ss(1) = ss(1) + pre*gx(R1x,R1y,R1z,R2x,R2y,R2z,a1,a2)
        ss(2) = ss(2) + pre*gy(R1x,R1y,R1z,R2x,R2y,R2z,a1,a2)
        ss(3) = ss(3) + pre*gz(R1x,R1y,R1z,R2x,R2y,R2z,a1,a2)
                        
     END DO
  END DO
  
  dipole = ss
  
! convert dipole vector from (1/Bohr radii) to (1/Angstroms)
  dipole = dipole/.5292D0
      
END SUBROUTINE atomDipoleMX
!*******************************************************************************
!*******************************************************************************
! Gaussian orbital expansion
!-------------------------------------------------------------------------------
REAL(8) FUNCTION gx (cg, cg1, cg3, cg5, cg7, cg9, cg11, cg13)
!-------------------------------------------------------------------------------
  IMPLICIT REAL(8) (c,t)
  ! maple 11 file: Gvector.mw  
  t1 = ((cg13) * (cg7))
  t5 = ((cg5) * (cg1))
  t8 = (cg1 ** 2)
  t9 = ((cg11) * (t8))
  t10 = (cg7 ** 2)
  t11 = ((cg13) * (t10))
  t12 = ((t11) * (cg))
  t15 = (cg5 ** 2)
  t16 = ((cg13) * (t15))
  t17 = ((t16) * (cg))
  t20 = (cg ** 2)
  t22 = ((cg11) * (t20) * (cg))
  t25 = (t20 ** 2)
  t27 = ((cg13) * (cg5))
  t30 = ((cg11) * (cg))
  t31 = (t15 ** 2)
  t38 = ((t15) * (cg5))
  t39 = ((cg13) * (t38))
  t42 = ((cg11) * (t20))
  t45 = ((cg1) * (cg7))
  t53 = ((cg11) * (cg1))
  t55 = ((t10) * (cg7) * (cg13))
  t60 = ((cg11) * (t8) * (cg1))
  t67 = -2 * t1 * cg * cg1 + 2 * t1 * t5 - 4 * t9 * t12 - 4 * t9 *&
       & t17 - 6 * t22 * t16 + 2 * cg11 * t25 * t27 - 2 * t30 * t31 * cg13&
       & + 2 * cg11 * cg7 * t5 + 2 * t9 * t39 + 6 * t42 * t39 - 2 * t30 * &
       &t45 - 2 * t22 * t11 - 2 * t30 * t16 * t10 + 2 * t53 * t55 * cg - 2&
       & * t60 * t1 * cg5 + 2 * t60 * t1 * cg
  t68 = cg5 * cg11
  t69 = t68 * t20
  t74 = t27 * t20
  t77 = cg * t15 * cg11
  t79 = t11 * cg5
  t85 = t1 * cg1
  t102 = -3 * t69 + cg * t10 * cg11 + 3 * t17 + t12 - 3 * t74 + 3 &
       &* t77 + 4 * t9 * t79 - 2 * t53 * t1 * t38 - 6 * t69 * t85 + 2 * t9&
       & * t74 - t68 * t8 - cg5 * t8 * cg13 - 2 * t53 * t55 * cg5 + 4 * t4&
       &2 * t79 + 6 * t77 * t85 + 2 * t22 * t85
  t105 = 1.772453851
  t109 = cg11 + cg13
  t110 = t109 ** 2
  t111 = t110 ** 2
  t112 = SQRT((t109))
  t116 = SQRT((t20 + t8))
  t120 = SQRT((t15 + t10))
  t124 = (cg3 ** 2)
  t125 = (cg9 ** 2)
  t134 = EXP(((2 * cg3 * cg9 - t124 - t20 - t15 - t10 - t125 +&
       & 2 * t45 - t8 + 2 * cg * cg5) * cg13 * cg11 / t109))
  gx = (t67 + t102) * (cg11) * (cg13) * t105 * 0.31415&
       &92654E1 / t112 / (t111) / t116 / t120 * t134
  RETURN
END FUNCTION gx

!-------------------------------------------------------------------------------
REAL(8) FUNCTION gy (cg, cg1, cg3, cg5, cg7, cg9, cg11, cg13)
!-------------------------------------------------------------------------------
  IMPLICIT REAL(8) (c,t)
  ! maple 11 file: Gvector.mw         
  t1 = (cg1 ** 2)
  t3 = ((cg11) * (t1) * (cg1))
  t4 = (cg5 ** 2)
  t5 = ((cg13) * (t4))
  t8 = (t1 ** 2)
  t10 = ((cg13) * (cg7))
  t13 = (cg7 ** 2)
  t14 = ((t13) * (cg13))
  t17 = ((cg13) * (cg5))
  t18 = ((cg7) * cg)
  t21 = ((cg13) * cg)
  t25 = ((cg11) * (cg1))
  t26 = (t13 ** 2)
  t30 = (cg * (cg5))
  t39 = (cg ** 2)
  t40 = ((cg11) * (t39))
  t41 = ((t13) * (cg7))
  t42 = ((cg13) * (t41))
  t46 = ((cg11) * (t39) * cg)
  t50 = ((t14) * (cg1))
  t53 = ((cg11) * cg)
  t55 = ((t4) * (cg5) * (cg13))
  t59 = ((t53) * (cg5))
  t60 = ((t10) * (t1))
  t65 = ((t5) * (cg1))
  t68 = -2 * t3 * t5 + 2 * cg11 * t8 * t10 - 6 * t3 * t14 + 2 * t17&
       & * t18 - 2 * t21 * cg5 * cg1 - 2 * t25 * t26 * cg13 - 2 * t25 * t30&
       & - 2 * t25 * t14 * t4 + 2 * cg11 * cg5 * t18 + 2 * t40 * t42 + 2   &
       & * t46 * t17 * cg1 - 4 * t40 * t50 - 2 * t53 * t55 * cg7 - 6 * t59 &
       & * t60 + 6 * t59 * t50 - 4 * t40 * t65
  t81 = cg11 * t1
  t86 = t5 * cg7
  t101 = 2 * t3 * t21 * cg5 + 2 * t53 * t55 * cg1 - 3 * t60 - cg11&
       & * cg7 * t39 + t65 + 3 * t25 * t13 + 3 * t50 - 3 * t81 * cg7 - t10&
       & * t39 + t25 * t4 + 4 * t40 * t86 + 4 * t81 * t86 - 2 * t53 * t17 &
       & * t41 + 2 * t40 * t60 + 6 * t81 * t42 - 2 * t46 * t17 * cg7
  t104 = 1.772453851
  t108 = cg11 + cg13
  t109 = t108 ** 2
  t110 = t109 ** 2
  t111 = SQRT((t108))
  t115 = SQRT((t4 + t13))
  t119 = SQRT((t39 + t1))
  t123 = (cg3 ** 2)
  t124 = (cg9 ** 2)
  t133 = EXP(((2 * cg3 * cg9 - t123 - t39 - t4 - t13 - t124 + &
       &2 * cg1 * cg7 - t1 + 2 * t30) * cg13 * cg11 / t108))
  gy = (t68 + t101) * (cg11) * (cg13) * t104 * 0.3141592654E1 &
       & / t111 / (t110) / t115 / t119 * t133
  RETURN
END FUNCTION gy

!-------------------------------------------------------------------------------
REAL(8) FUNCTION gz (cg, cg1, cg3, cg5, cg7, cg9, cg11, cg13)
!-------------------------------------------------------------------------------
  IMPLICIT REAL(8) (c,t)      
  ! maple 11 file: Gvector.mw       
  t1 = ((cg13) * (cg11))
  t2 = (cg5 ** 2)
  t9 = (cg7 ** 2)
  t14 = (cg ** 2)
  t25 = (cg1 ** 2)
  t37 = ((cg) * (cg5))
  t43 = ((cg1) * (cg7))
  t63 = 2 * t1 * cg * t2 * cg5 - cg5 * cg11 * cg + 2 * t1 * cg1 * &
       &t9 * cg7 - 2 * t1 * t14 * t9 - 4 * t1 * t14 * t2 + 2 * t1 * t14 * &
       &cg1 * cg7 - 2 * t1 * t25 * t2 - 4 * t1 * cg * cg5 * cg1 * cg7 - 4 &
       &* t1 * t25 * t9 + 2 * t1 * t37 * t9 - cg11 * cg1 * cg7 + 2 * t1 * &
       &t43 * t2 + 2 * t1 * t25 * cg1 * cg7 + 2 * t1 * t25 * cg * cg5 + 2 &
       &* t1 * t14 * cg * cg5 - cg1 * cg13 * cg7 - cg * cg13 * cg5
  t68 = 1.772453851
  t70 = cg11 + cg13
  t71 = t70 ** 2
  t72 = t71 ** 2
  t73 = SQRT((t70))
  t78 = SQRT((t14 + t25))
  t81 = SQRT((t2 + t9))
  t86 = (cg3 ** 2)
  t87 = (cg9 ** 2)
  t95 = EXP(((2 * cg3 * cg9 - t86 - t14 - t2 - t9 - t87 + 2 * &
       &t43 - t25 + 2 * t37) * cg13 * cg11 / t70))
  gz = (t63) * (cg13) * (cg11) * (cg3 - cg9) * t68&
       & * 0.3141592654E1 / t73 / (t72) / t78 / t81 * t95
  RETURN
END FUNCTION gz
!*******************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional subroutines written be Daria Satco (2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*******************************************************************************
SUBROUTINE realDielEn(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,ebg,fwhm,ne,hw,eps1)
!===============================================================================
! Compute the real part of the dielectric function as a function
! of probe photon energy
! expression from Lin, M. F., and Kenneth W-K. Shung. "Plasmons and optical properties of carbon nanotubes."
! Physical Review B 50.23 (1994): 17744.
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
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in) :: epol(3)

  REAL(8), INTENT(in)    :: fwhm
  REAL(8), INTENT(in)    :: ebg

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: eps1(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  REAL(8)                :: eps0(ne)

  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: diameter, tubeDiam, area
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta_eps1

  COMPLEX(8)             :: css

! first pass

!background dielectric constant (dimensionless)
eps0(:) = ebg

! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! dielectric function prefactor (eV**3 Angstroms**3)
! --------------------------------------------------
! see for the prefactor expression paper:
! Sanders, G. D., et al.
! "Resonant coherent phonon spectroscopy of single-walled carbon nanotubes."
! Physical Review B 79.20 (2009): 205434.

     diameter = tubeDiam(n,m)        !(Angstroms)
     area = pi*(diameter/2.D0)**2    !(Angstroms**2)
     pre  = 8.D0*pi*e2*hbarm**2/area !(eV**3 Angstroms**3)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss = 0.D0

  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction
     DO mu1 = 1, nhex
        DO n2 = 1, 2
           DO mu2 = 1, nhex

              IF (n1 == n2 .AND. mu1 == mu2) CYCLE

              DO k=1,nk

!energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k))  !(1/Angstroms**2)

                 ! if small p2df then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

                 IF (ABS(p2df) <= ptol) CYCLE

!**********************************************************
! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1

! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk)
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! accumulate dielectric function vs photon energy
                 DO ie = 1, ne
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)

                    diracAvg = diracDelta_eps1(x1,y1,x2,y2,fwhm)  ! (1/eV)
                    IF(diracAvg == 0.) CYCLE

                    IF (hw(ie) .le. ptol) THEN
                        ss(ie) = ss(ie) + dk/(2*pi) * p2df * diracAvg / ( 1.D-3 * Eab )
                    ELSE
                        ss(ie) = ss(ie) + dk/(2*pi) * p2df * diracAvg / ( hw(ie) * Eab )   ! eV**(-3) * Angstroms**(-3)
                    END IF

                 END DO

              END DO
           END DO
        END DO
     END DO
  END DO

! real part of dielectric function (dimensionless)
  eps1 = eps0 + pre*ss

END SUBROUTINE realDielEn
!*******************************************************************************
!*******************************************************************************

SUBROUTINE realDielEn_met2(ne,ebg,hw,eps2,eps1)
!===============================================================================
! Compute the real part of the dielectric function as a function
! of probe photon energy
! according to Kramers-Kronig relations for dielectric function
!-------------------------------------------------------------------------------
! Input        :
!  ne            number of probe photon energies
!  ebg           background dielectric permittivity
!  hw(ne)        array of probe photon energies (eV)
!  eps2(ne)      imaginary part of dielectric function (none)
! Output       :
!  eps1(ne)      real part of dielectric function (none)
!===============================================================================
  IMPLICIT NONE

  ! input variables
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: ebg
  REAL(8), INTENT(in)    :: hw(ne)
  REAL(8), INTENT(in)    :: eps2(ne)

  ! output variable
  REAL(8),  INTENT(out)  :: eps1(ne)

  ! working variables and parameters
  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), ALLOCATABLE   :: ss(:)
  REAL(8)                :: eps0(ne)
  INTEGER :: ie, ii

! begin subroutine

! background dielectric constant (dimensionless)
eps0(:) = ebg

IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss = 0.D0

DO ie = 1, ne
    DO ii = 1,ne
        IF (ii == ie) CYCLE
        ss(ie) = ss(ie) + ABS( hw(1)-hw(2) ) * 2.D0 / pi * eps2(ii) * hw(ii)/( hw(ii)**2 - hw(ie)**2 )
    END DO
END DO

eps1 = eps0 + ss

END SUBROUTINE realDielEn_met2
!*******************************************************************************
!*******************************************************************************

SUBROUTINE imagDielEn_met2(ne,hw,eps1,eps2)
!===============================================================================
! Compute the imaginary part of the dielectric function as a function
! of probe photon energy
! according to Kramers-Kronig relations for dielectric function
!-------------------------------------------------------------------------------
! Input        :
!  ne            number of probe photon energies
!  hw(ne)        array of probe photon energies (eV)
!  eps1(ne)      real part of dielectric function (none)
! Output       :
!  eps2(ne)      imaginary part of dielectric function (none)
!===============================================================================
  IMPLICIT NONE

  ! input variables
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)
  REAL(8), INTENT(in)    :: eps1(ne)

  ! output variable
  REAL(8),  INTENT(out)  :: eps2(ne)

  ! working variables and parameters
  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), ALLOCATABLE   :: ss(:)
  INTEGER :: ie, ii

! begin subroutine

IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss = 0.D0

DO ie = 1, ne
    DO ii = 1,ne
        IF (ii == ie) CYCLE
        ss(ie) = ss(ie) - ABS( hw(1)-hw(2) ) * 2.D0 / pi * (eps1(ii) - eps1(ne))/( hw(ii) - hw(ie) )
    END DO
END DO

eps2 = ss

END SUBROUTINE imagDielEn_met2
!*******************************************************************************
!*******************************************************************************

SUBROUTINE EELS(ne,eps1,eps2,eelspec)
!===============================================================================
! Compute the electron energy loss spectra
!  eels ~ Im(-1/eps) = eps2 / eps^2
!-------------------------------------------------------------------------------
! Input        :
!  ne            number of probe photon energies
!  eps1(ne)      real part of dielectric function (none)
!  eps2(ne)      imaginary part of dielectric function (none)
! Output       :
!  eels(ne)      eels spectrum (none)
!===============================================================================
  IMPLICIT NONE

  ! input variables
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: eps1(ne)
  REAL(8),  INTENT(in)   :: eps2(ne)

  ! output variable
  REAL(8),  INTENT(out)  :: eelspec(ne)

  ! working variables and parameters
  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  INTEGER :: ie

! begin subroutine

  DO ie = 1, ne
        eelspec(ie) = eps2(ie)/( eps1(ie)**2 + eps2(ie)**2 )
  END DO

END SUBROUTINE EELS
!*******************************************************************************
!*******************************************************************************
SUBROUTINE RealDynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm1)
!===============================================================================
! Compute the real part of the dynamical conductivity as a function
! of probe photon energy
! expression from Sasaki, Ken-ichi, and Yasuhiro Tokura. "Theory of a Carbon-Nanotube Polarization Switch."
! Physical Review Applied 9.3 (2018): 034018.
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
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in) :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: sigm1(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: h = 4.135D-15         !(eV-s)
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: diameter, tubeDiam
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta

  COMPLEX(8)             :: css

! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! conductivity function prefactor (eV**2 Angstroms**3) * [e^2/h] = (A/s)
     diameter = tubeDiam(n,m)        !(Angstroms)
     pre  = 32.D0*hbarm**2/diameter * e2/h !(eV**2 Angstroms**3) * [e^2/h] = (A/s)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss=0.D0

  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction
     DO mu1 = 1, nhex
        DO n2 = 1, 2
           DO mu2 = 1, nhex

              IF (n1 == n2 .AND. mu1 == mu2) CYCLE

              DO k=1,nk

!energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k)) !(1/Angstroms**2)

                 ! if small then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

                 IF (ABS(p2df) <= ptol) CYCLE

! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1

! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk)
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! accumulate dielectric function vs photon energy
                 DO ie = 1, ne
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)
                    diracAvg = diracDelta(x1,y1,x2,y2,fwhm)    !(1/eV)
                    IF(diracAvg == 0.) CYCLE

                    ss(ie) = ss(ie) + dk/2*p2df*diracAvg/Eab     ! (1/Angstroms**3 1/eV**2)
                 END DO

              END DO
           END DO
        END DO
     END DO
  END DO

! imaginary part of dielectric function [e^2/h] = (A/s)
  sigm1 = pre*ss

END SUBROUTINE RealDynConductivity
!*******************************************************************************
!*******************************************************************************
SUBROUTINE ImagDynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm2)
!===============================================================================
! Compute the imaginary part of the dynamical conductivity as a function
! of probe photon energy
! expression from Sasaki, Ken-ichi, and Yasuhiro Tokura. "Theory of a Carbon-Nanotube Polarization Switch."
! Physical Review Applied 9.3 (2018): 034018.
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
!  sigm2(ne)     imaginary part of dynamical conductivity (e^2/h)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in) :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: sigm2(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: h     = 4.135D-15     !(eV-s)
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: diameter, tubeDiam
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta_eps1

  COMPLEX(8)             :: css


! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! conductivity function prefactor (eV**2 Angstroms**3) * [e^2/h] = (A/s)
     diameter = tubeDiam(n,m)        !(Angstroms)
     pre  = 32.D0*hbarm**2/diameter * e2/h !(eV**2 Angstroms**3) * [e^2/h] = (A/s)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss = 0.D0

  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction
     DO mu1 = 1, nhex
        DO n2 = 1, 2
           DO mu2 = 1, nhex

              IF (n1 == n2 .AND. mu1 == mu2) CYCLE

              DO k=1,nk

!energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k))

                 ! if small p2df then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

                 IF (ABS(p2df) <= ptol) CYCLE

!**********************************************************
! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1

! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk)
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! accumulate dielectric function vs photon energy
                 DO ie = 1, ne
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)

                    diracAvg = diracDelta_eps1(x1,y1,x2,y2,fwhm) ! (1/eV)
                    IF(diracAvg == 0.) CYCLE

                    ss(ie) = ss(ie) + dk/(2*pi) * p2df * diracAvg/Eab  ! eV**(-3) * Angstroms**(-3)
                 END DO

              END DO
           END DO
        END DO
     END DO
  END DO

! real part of dielectric function (dimensionless)
  sigm2 = -pre*ss

END SUBROUTINE ImagDynConductivity
!*******************************************************************************
!*******************************************************************************
SUBROUTINE Absorption(ne,eps1,eps2,sigm1,sigm2,absorpt)
!===============================================================================
! Calculates absorption as Re(sigma/epsilon)
!-------------------------------------------------------------------------------
! Input        :
! eps1(ne)      real part of dielectric function (dimensionless)
! eps2(ne)      imaginary part of dielectric function (dimensionless)
! sigm1(ne)     real part of dynamical conductivity (e^2/h)
! sigm2(ne)     imaginary part of dynamical conductivity (e^2/h)
! Output       :
! absorpt(ne)   light absorption (e^2/h)
!===============================================================================
  IMPLICIT NONE

! Input variables
  INTEGER, INTENT(in)   :: ne

  REAL(8),  INTENT(in)  :: eps1(ne)
  REAL(8),  INTENT(in)  :: eps2(ne)
  REAL(8),  INTENT(in)  :: sigm1(ne)
  REAL(8),  INTENT(in)  :: sigm2(ne)

! Output variables
  REAL(8),  INTENT(out) :: absorpt(ne)

! working variables
  INTEGER               :: ie

  DO ie = 1, ne
        absorpt(ie) = ( sigm1(ie) * eps1(ie) + sigm2(ie) * eps2(ie) )/( eps1(ie)**2 + eps2(ie)**2 )
  END DO

END SUBROUTINE Absorption
!*******************************************************************************
!*******************************************************************************
SUBROUTINE ImagDynConductivityInter(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm2)
!===============================================================================
! Compute the imaginary part of the dynamical conductivity as a function
! of probe photon energy
! only interband transition are considered in this subroutine
! expression from Sasaki, Ken-ichi, and Yasuhiro Tokura. "Theory of a Carbon-Nanotube Polarization Switch."
! Physical Review Applied 9.3 (2018): 034018.
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
!  sigm2(ne)     imaginary part of interband dynamical conductivity (e^2/h)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in) :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: sigm2(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: h     = 4.135D-15     !(eV-s)
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: diameter, tubeDiam
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta_eps1

  COMPLEX(8)             :: css

! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! conductivity function prefactor (eV**2 Angstroms**3) * [e^2/h] = (A/s)
     diameter = tubeDiam(n,m)        !(Angstroms)
     pre  = 32.D0*hbarm**2/diameter * e2/h !(eV**2 Angstroms**3) * [e^2/h] = (A/s)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss = 0.D0

  DO n1 = 1,2   ! 1 <-> valence, 2 <-> conduction
    !interband transitions
    IF (n1 == 1) THEN
      n2 = 2
    ELSE
      n2 = 1
    END IF

     DO mu1 = 1, nhex
           DO mu2 = 1, nhex

              DO k=1,nk

! energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k))

                 ! if small p2df then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

                 IF (ABS(p2df) <= ptol) CYCLE

!**********************************************************
! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1

! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk)
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! accumulate dielectric function vs photon energy
                 DO ie = 1, ne
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)

                    diracAvg = diracDelta_eps1(x1,y1,x2,y2,fwhm) ! (1/eV)
                    IF(diracAvg == 0.) CYCLE

                    ss(ie) = ss(ie) + dk/(2*pi) * p2df * diracAvg/Eab  ! eV**(-3) * Angstroms**(-3)
                 END DO

              END DO
           END DO
        END DO
      END DO

! real part of dielectric function (dimensionless)
  sigm2 = -pre*ss

END SUBROUTINE ImagDynConductivityInter
!*******************************************************************************
!*******************************************************************************
SUBROUTINE ImagDynConductivityIntra(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,sigm2)
!===============================================================================
! Compute the imaginary part of the dynamical conductivity as a function
! of probe photon energy
! only intraband transition are considered in this subroutine
! expression from Sasaki, Ken-ichi, and Yasuhiro Tokura. "Theory of a Carbon-Nanotube Polarization Switch."
! Physical Review Applied 9.3 (2018): 034018.
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
!  sigm2(ne)     imaginary part of intraband dynamical conductivity (e^2/h)
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: n, m
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in) :: epol(3)

  REAL(8), INTENT(in)    :: fwhm

  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)

! output variable
  REAL(8),  INTENT(out)  :: sigm2(ne)

! working variables and parameter

  REAL(8), SAVE          :: pre

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), ALLOCATABLE   :: ss(:) !(ne)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)     e2 = e^2
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: h = 4.135D-15         !(eV-s)
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: diameter, tubeDiam
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta_eps1


  COMPLEX(8)             :: css

! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! conductivity function prefactor (eV**2 Angstroms**3) * [e^2/h] = (A/s)
     diameter = tubeDiam(n,m)        !(Angstroms)
     pre  = 32.D0*hbarm**2/diameter * e2/h !(eV**2 Angstroms**3) * [e^2/h] = (A/s)

! sum over n1, mu1, n2, mu2 and k
  IF (ALLOCATED(ss) .EQV. .TRUE.) DEALLOCATE(ss)
  ALLOCATE(ss(ne))
  ss = 0.D0

  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction

     n2 = n1     ! intraband transitions

     DO mu1 = 1, nhex
           DO mu2 = 1, nhex

              IF (mu1 == mu2) CYCLE

              DO k=1,nk

! energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k))

                 ! if small p2df then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

                 IF (ABS(p2df) <= ptol) CYCLE

!**********************************************************
! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1

! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk)
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! accumulate dielectric function vs photon energy
                 DO ie = 1, ne
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)

                    diracAvg = diracDelta_eps1(x1,y1,x2,y2,fwhm) ! (1/eV)
                    IF(diracAvg == 0.) CYCLE

                    ss(ie) = ss(ie) + dk/(2*pi) * p2df * diracAvg/Eab  ! eV**(-3) * Angstroms**(-3)
                 END DO

              END DO
           END DO
        END DO
      END DO

! real part of dielectric function (dimensionless)
  sigm2 = -pre*ss

END SUBROUTINE ImagDynConductivityIntra
!*******************************************************************************
!*******************************************************************************
SUBROUTINE Func_test(nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,fwhm,ne,hw,&
difFermiDist,matrElementSq,diracAvgFunc)
!===============================================================================
! Subroutine to calculate functions of main importance
! Just to control if everything is reasonable
!-------------------------------------------------------------------------------
! Input        :
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
!  difFermiDist(2,2,nhex,nhex,nk)   Fermi distribution as a function of k
!  matrElementSq(2,2,nhex,nhex,nk)  matrix element square as a function of k
!  diracAvgFunc(2,2,nhex,nhex,nk,ne)   Lorentzian-like function as a function of k
!===============================================================================
  IMPLICIT NONE

! input variables
  INTEGER, INTENT(in)    :: nhex
  INTEGER, INTENT(in)    :: nk

  REAL(8), INTENT(in)    :: rka(nk)
  REAL(8), INTENT(in)    :: Enk(2,nhex,nk)
  COMPLEX(8), INTENT(in) :: cDipole(3,nk,2,nhex,2,nhex)

  REAL(8), INTENT(in)    :: Tempr
  REAL(8), INTENT(in)    :: Efermi

  REAL(8), INTENT(in)    :: epol(3)

  REAL(8), INTENT(in)    :: fwhm
  INTEGER, INTENT(in)    :: ne
  REAL(8), INTENT(in)    :: hw(ne)


! output variables
  REAL(8)                :: difFermiDist(2,2,nhex,nhex,nk)
  REAL(8)                :: matrElementSq(2,2,nhex,nhex,nk)
  REAL(8)                :: diracAvgFunc(2,2,nhex,nhex,nk,ne)

! working variables and parameter

  REAL(8), SAVE, ALLOCATABLE :: fnk(:,:,:)  !(2,nhex,nk)

  REAL(8), PARAMETER     :: pi    =  3.14159265358979D0
  REAL(8), PARAMETER     :: e2    = 14.4          !(eV-A)
  REAL(8), PARAMETER     :: hbarm =  7.62         !(eV-A**2)
  REAL(8), PARAMETER     :: h     = 4.13D-15      !(eV-s)
  REAL(8), PARAMETER     :: ptol  =  1.D-15

  INTEGER                :: k, mu, ii, ie
  INTEGER                :: n1, mu1, n2, mu2

! for calling some functions
  REAL(8)                :: dk, rkT
  REAL(8)                :: fermi, Ekk
  REAL(8)                :: p2, p2df, x1, x2, enk1n, enk1p, enk2n, enk2p
  REAL(8)                :: y1, y2, diracAvg, Eab, diracDelta_eps1

  COMPLEX(8)             :: css

! initial distribution functions (dimensionless)
     IF (ALLOCATED(fnk) .EQV. .TRUE.) DEALLOCATE(fnk)
     ALLOCATE(fnk(2,nhex,nk))
     rkT = .025853D0 * (Tempr/300.D0) ! thermal energy
     DO ii = 1, 2
        DO mu = 1, nhex
           DO k = 1, nk
              Ekk = Enk(ii,mu,k)
              fnk(ii,mu,k) = fermi(Ekk,Efermi,rkT)
           END DO
        END DO
     END DO

! LOOP over n1, mu1, n2, mu2 and k
  difFermiDist = 0.D0
  matrElementSq = 0.D0
  diracAvgFunc = 0.D0

  DO n1 = 1, 2   ! 1 <-> valence, 2 <-> conduction
    DO n2 = 1,2

     DO mu1 = 1, nhex
           DO mu2 = 1, nhex

              IF (n1 == n2 .and. mu1 == mu2) CYCLE

              DO k=1,nk

! k-cell boundaries (1/A)
                 IF (k == 1) THEN
                    x1 = rka(1)
                    x2 = (rka(k+1) + rka(k))/2.D0
                 ELSE IF (k == nk) THEN
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = rka(nk)
                 ELSE
                    x1 = (rka(k-1) + rka(k))/2.D0
                    x2 = (rka(k+1) + rka(k))/2.D0
                 END IF
                 dk = x2-x1

! band energies at k-cell boundaries (eV)
                 IF (k == 1) THEN
                    enk1n = Enk(n1,mu1,1)
                    enk1p = (Enk(n1,mu1,1)+Enk(n1,mu1,2))/2.D0
                    enk2n = Enk(n2,mu2,1)
                    enk2p = (Enk(n2,mu2,1)+Enk(n2,mu2,2))/2.D0
                 ELSE IF (k == nk) THEN
                    enk1n = (Enk(n1,mu1,nk-1)+Enk(n1,mu1,nk))/2.D0
                    enk1p = Enk(n1,mu1,nk)
                    enk2n = (Enk(n2,mu2,nk-1)+Enk(n2,mu2,nk))/2.D0
                    enk2p = Enk(n2,mu2,nk)
                 ELSE
                    enk1n = (Enk(n1,mu1,k-1)+Enk(n1,mu1,k  ))/2.D0
                    enk1p = (Enk(n1,mu1,k  )+Enk(n1,mu1,k+1))/2.D0
                    enk2n = (Enk(n2,mu2,k-1)+Enk(n2,mu2,k  ))/2.D0
                    enk2p = (Enk(n2,mu2,k  )+Enk(n2,mu2,k+1))/2.D0
                 END IF

! energy difference
                 Eab = Enk(n2,mu2,k) - Enk(n1,mu1,k)

! squared optical dipole matrix element (1/Angstroms**2)
                 css = 0.D0
                 !cycle for scalar product of polarization vector and dipole matrix element
                 DO ii = 1, 3
                    css = css + epol(ii)*cDipole(ii,k,n1,mu1,n2,mu2)
                 END DO
                 ! square of matrix element
                 p2   = CDABS(css)**2
                 matrElementSq(n1,n2,mu1,mu2,k) = sqrt(p2)
                 ! multiply by distribuion function
                 p2df = p2*(fnk(n1,mu1,k) - fnk(n2,mu2,k))
                 difFermiDist(n1,n2,mu1,mu2,k) = fnk(n1,mu1,k) - fnk(n2,mu2,k)

                 ! if small p2df then skip
                 IF ( ABS(Eab) .lt. ptol ) THEN
                    IF (ABS(p2df) > ptol) STOP 'something strange'
                 ENDIF

!**********************************************************
                DO ie = 1, ne
                    y1 = enk2n - enk1n - hw(ie)
                    y2 = enk2p - enk1p - hw(ie)

                    diracAvg = diracDelta_eps1(x1,y1,x2,y2,fwhm) ! (1/eV)
                    diracAvgFunc(n1,n2,mu1,mu2,k,ie) = diracAvg/Eab
                END DO

              END DO
           END DO
        END DO
      END DO
    END DO

END SUBROUTINE Func_test
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tbDipolXYOrt(n,m,n1,mu1,n2,mu2,rk,xDipole,yDipole)
!===============================================================================
! another version of tbDipolXY -> different algorithm but yields same results
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
  REAL(8), PARAMETER     :: tol = 1.D-15

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
  INTEGER,ALLOCATABLE           :: j1j2(:,:)

  INTEGER                :: nhex, nHexagon, ier, isel, iflip
  INTEGER                :: mmu1,mmu2,nn1,nn2
  INTEGER                :: iatom, jatom, nn, ivec, j1, j2, NNatom, jj

  REAL(8)    rkk, phi1, phi2

! check input for errors
  nhex=nHexagon(n,m)
  IF (ALLOCATED(j1j2) .EQV. .TRUE.) DEALLOCATE(j1j2)
  ALLOCATE(j1j2(nhex,2))

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
  CALL etbTubeBand(n,m,mmu1,rkk,Ek1,Zk1)
  CALL etbTubeBand(n,m,mmu2,rkk,Ek2,Zk2)
  CALL getHexagonPosition(n,m,nhex,j1j2)

! compute x and y components of dipole vector (1/Angstroms)
DO jj = 1,nhex
  DO iatom = 1, 2
     DO nn = 1, 4
        DO ivec = 1, nvecs(nn)

           CALL NNj1j2(iatom,ivec,nn,j1,j2)
           CALL phaseFactor(n,m,j1j2(jj,1),j1j2(jj,2),-rk,-mmu1,phi1)
           CALL phaseFactor(n,m,j1j2(jj,1)+j1,j1j2(jj,2)+j2,rk,mmu2,phi2)

           jatom = NNatom(iatom,nn)
           CALL atomDipoleMX(n,m,iatom,j1j2(jj,1),j1j2(jj,2),jatom,j1j2(jj,1)+j1,j1j2(jj,2)+j2,dipole)

           c1 = CONJG( Zk1(iatom,nn1) )
           c2 = Zk2(jatom,nn2)*CDEXP(ci*(phi1+phi2))

           xDipole = xDipole+c1*c2*(dipole(1)-ci*dipole(2))

        END DO
     END DO
  END DO
END DO

  xDipole = xDipole/2.D0*1.D0/nhex
  yDipole = ci*xDipole

! use symmetry relation to reverse the exchange bra and ket vectors
  IF (iflip == 1) THEN
     xDipole = -CONJG(xDipole)
     yDipole = -CONJG(yDipole)
  END IF

  RETURN
END SUBROUTINE tbDipolXYOrt
!*******************************************************************************
!*******************************************************************************
