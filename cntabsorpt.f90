!*******************************************************************************
!*******************************************************************************
! Project      : cntabsorpt.f90
!===============================================================================
! Purpose      :
! Electronic dispersion, dielectric properties and absorption spectra of SWNTs
!-------------------------------------------------------------------------------
! Method       :
! [ Electron ]   Third nearest-neighbor ETB model (but w/o sigma band effect)
! [ Optic    ]   Dielectric functions and absorption acording to Kubo formula
!-------------------------------------------------------------------------------
! Authors      :
! - ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
! - Daria Satco  (dasha.shatco@gmail.com)
! Latest Vers. : 2018.09.30
!-------------------------------------------------------------------------------
! Required files :
! - globvar         -- global variables
! - tubepar.f90     -- nanotube parameters
! - tubestruct.f90  -- nanotube structure libraries
! - libMath.f90     -- mathematical libraries
! - libswntElec.f90 -- electronic/excitonic states
! - libswntOpt.f90  -- optical matrix element, dielectric function, conductivity
!*******************************************************************************
!*******************************************************************************
PROGRAM cntabsorpt
!===============================================================================
  USE globvar
  IMPLICIT NONE
  
! parameters
  REAL(8), PARAMETER     :: pi    = 3.14159265358979D0
  REAL(8), PARAMETER     :: hbar  = 6.582D-4 !(eV-ps)
  REAL(8), PARAMETER     :: h     = 4.13D-15 !(eV-s)
  REAL(8), PARAMETER     :: e2    = 14.4     !(eV-A)
  REAL(8), PARAMETER     :: hbarm =  7.62    !(eV-A**2)  hbarm = h^2 / m
  REAL(8), PARAMETER     :: hbarc = 1.97D-5  !(eV cm)
  REAL(8), PARAMETER     :: hbarvfermi = 6.582119 !(eV-A) !hbar*vfermi, vfermi = 10**6 m/s
  REAL(8), PARAMETER     :: ptol  =  1.D-15
!-------------------------------------------------------------------------------
! variables for calling tube structure function
  INTEGER                :: n, m, nHexagon
  REAL(8)                :: trLength
  REAL(8)                :: fermiLevel
!-------------------------------------------------------------------------------
! calculation parameter input
  INTEGER                :: k, mu, ii, ie, nee, nq, nep
  INTEGER                :: nhw_laser
  REAL(8)                :: Tempr, rkT, doping, emine, emaxe, eminp, emaxp
  REAL(8)                :: epmin, epmax
  REAL(8)                :: laser_theta, laser_fwhm, dummy
  REAL(8)                :: ebg
  REAL(8)                :: Efermi
!-------------------------------------------------------------------------------
! variables for electronic states
  REAL(8)                :: rkmin, rkmax, dk, rk, de, eii
  INTEGER                :: nk1
  REAL(8), ALLOCATABLE   :: rka(:)           !(nk)
  REAL(8), ALLOCATABLE   :: Enk(:,:,:)       !(2,nhex,nk)
  COMPLEX(8), ALLOCATABLE:: Znk(:,:,:,:)
  REAL(8), DIMENSION(2)  :: En               !(unit eV)
  REAL(8), ALLOCATABLE   :: Enkt(:,:,:)      !(2,nhex,nk)
  COMPLEX(8)             :: Zk(2,2)
  
  REAL(8), ALLOCATABLE   :: eEarray(:)       !(nee)
  REAL(8), ALLOCATABLE   :: eDOS(:)          !(nee)
!-------------------------------------------------------------------------------
! variables for optical properties
  INTEGER                 :: mmu
  INTEGER                 :: n1, mu1, n2, mu2
  REAL(8)                 :: dep
  REAL(8), ALLOCATABLE    :: Px2k(:,:)       !(nk,nhex)
  REAL(8), ALLOCATABLE    :: Pz2k(:,:)       !(nk,nhex)
 
  REAL(8)                 :: epol(3)

  REAL(8), ALLOCATABLE    :: hw_laser(:)     !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps2(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps2kk(:)       !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps1(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps1kk(:)       !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eelspec(:)      !(nhw_laser)
  REAL(8), ALLOCATABLE    :: alpha(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm1(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm1_intra(:)  !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm1_inter(:)  !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2_intra(:)  !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2_inter(:)  !(nhw_laser)
  REAL(8), ALLOCATABLE    :: absorpt(:)      !(nhw_laser)
  COMPLEX(8), ALLOCATABLE :: cDipole(:,:,:,:,:,:) !(3,nk,2,nhex,2,nhex)
  REAL(8), ALLOCATABLE    :: absorptPart(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: eps1Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: eps2Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: sigm1Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: sigm2Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: eps1Dr(:), eps2Dr(:)
  REAL(8), ALLOCATABLE    :: sigm1Dr(:), sigm2Dr(:)

  REAL(8)                 :: charge
!-------------------------------------------------------------------------------
! variables for input and output files 
  CHARACTER(40)           :: infile, outfile, path
  CHARACTER(5)            :: fermistr, thetastr

! working variales
  INTEGER                 :: zeroNumber
  REAL(8), ALLOCATABLE    :: plasmonFreq(:), eps1Zeros(:)
  INTEGER, ALLOCATABLE    :: plasmonPosition(:)
  REAL(8)                 :: diameter, area, prediel, precond, tubeDiam, rkii
  REAL(8)                 :: ratio
  INTEGER, ALLOCATABLE    :: muii(:,:)
  INTEGER                 :: place_max, zeroExists, metal, loc_mu1(2), loc_mu2(2), plasmonExists, truePlasmon, ierr
  INTEGER                 :: i, mn, j, mnj, l, s, k1, k2, b1, b2

!------------------------------------------------------------------------------
!temp variables
  REAL(8), ALLOCATABLE   :: difFermiDist(:,:,:,:,:)
  REAL(8), ALLOCATABLE   :: matrElementSq(:,:,:,:,:)
  REAL(8), ALLOCATABLE   :: diracAvgFunc(:,:,:,:,:,:)
  INTEGER                :: max_position(5), min_position(5)

  REAL                   :: divergence(9), kCoef(10), maxAbsDif(9)
  REAL(8), ALLOCATABLE   :: eps2aii(:,:)        !(nhw_laser)
  REAL(8)                :: rk1, rk2
  COMPLEX(8)             :: test(3)

!*************************************************************
!!!!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!
!*************************************************************

     WRITE (*,*) 'Input file name MUST be: tube.param.**** '
     WRITE (*,*) '**** can be up to 20 symbols length string'
     WRITE (*,*) 'Enter suffix of input data set (1a20):'
     READ  (*,2002) infile
     outfile = infile

     CALL EXECUTE_COMMAND_LINE( 'mkdir -p tube'//TRIM(outfile) )
     path = './tube'//TRIM(outfile)//'/'
!----------------------------------------------------------------------
!                       read data from file
!----------------------------------------------------------------------
     OPEN(unit=22,file='tube.param.'//infile, status = 'old')

     READ (22,*) Tempr
     READ (22,*) n,m
     READ (22,*) nk
     READ (22,*) doping
     READ (22,*) Efermi
     READ (22,*) refrac
     READ (22,*) ebg
     READ (22,*) nhw_laser
     READ (22,*) epmin,epmax
     READ (22,*) laser_fwhm
     READ (22,*) nee,emine,emaxe
     READ (22,*) laser_theta

     CLOSE(unit=22)
      
  WRITE (*,*) '====================================================='
  CALL printTubeClass(n,m,6)
  WRITE (*,*) '====================================================='
  WRITE (*,*) ' Here is the input data:'
  WRITE (*,*) '====================================================='
  WRITE (*,*) 'Temperature (deg K)            :',Tempr
  WRITE (*,*) 'Chiral indices (n,m)           :',n,m
  WRITE (*,*) 'Doping (n-p) per length (1/A)  :',doping
  WRITE (*,*) 'Fermi level                    :',Efermi
  WRITE (*,*) 'Refractive index               :',refrac
  WRITE (*,*) 'Background dielectric permittivity:', ebg
  WRITE (*,*) '-----------------------------------------------------'
  WRITE (*,*) 'Number of laser photon energies:',nhw_laser
  WRITE (*,*) 'Laser photon energy range (eV) :',epmin,epmax
  WRITE (*,*) 'Laser linewidth (eV)           :',laser_fwhm
  WRITE (*,*) 'Laser polarization angle (deg) :',laser_theta
  WRITE (*,*) '-----------------------------------------------------'      
  WRITE (*,*) 'Electron k points, nk          :',nk
  WRITE (*,*) 'Electron DOS energies, nee     :',nee
  WRITE (*,*) 'Electron DOS energy range (eV) :',emine,emaxe
  WRITE (*,*) '-----------------------------------------------------' 

!----------------------------------------------------------------------
!                           write log file
!----------------------------------------------------------------------
    OPEN(unit=22,file=TRIM(path)//'tube.log.'//outfile)
      
    rkT = .025853D0 * (Tempr/300.D0)         

    WRITE (22,*) '===================================================='
    WRITE (22,*) '           tube.log.'//outfile
    WRITE (22,*) '===================================================='
    CALL printTubeClass(n,m,22)
    WRITE (22,*) '===================================================='
    WRITE (22,*) 'Temperature (deg K)            :',Tempr
    WRITE (22,*) 'Thermal energy (eV)            :',rkT
    WRITE (22,*) 'Chiral indices (n,m)           :',n,m
    WRITE (22,*) 'Doping (n-p) per length (1/A)  :',doping
    WRITE (22,*) 'Fermi level                    :',Efermi
    WRITE (22,*) 'Refractive index               :',refrac
    WRITE (22,*) 'Background dielectric permittivity:', ebg
    WRITE (22,*) '----------------------------------------------------'
    WRITE (22,*) 'Number of laser photon energies:',nhw_laser
    WRITE (22,*) 'Laser photon energy range (eV) :',epmin,epmax
    WRITE (22,*) 'Laser linewidth (eV)           :',laser_fwhm
    WRITE (22,*) 'Laser polarization angle (deg) :',laser_theta              
    WRITE (22,*) '----------------------------------------------------'  
    WRITE (22,*) 'Electron k points, nk          :',nk
    WRITE (22,*) 'Electron DOS energies, nee     :',nee
    WRITE (22,*) 'Electron DOS energy range (eV) :',emine,emaxe
    WRITE (22,*) '----------------------------------------------------'

    CLOSE(unit=22)

!----------------------------------------------------------------------
!                          begin execution
!----------------------------------------------------------------------
  WRITE (*,*) '====================================================='
  WRITE (*,*) '                begin execution'
  WRITE (*,*) '====================================================='      

!----------------------------------------------------------------------
!                   equilibrium nanotube energy bands (eV)
!----------------------------------------------------------------------
  WRITE(*,*) '..Electronic bands in ETB model'
    
! allocate storage for energy bands
  nhex=nHexagon(n,m)

! choose the correct number of k-points
! it depends on broadening factor and unit cell translational vector length
! the dk is based on broadening, but also multiplied by 1/5
! the prefactor 1/5 was chosen after several checks of calculational convergence
! the resulting accuracy should be at least 2 signs
  dk = laser_fwhm/(5*hbarvfermi)
  nk1 = INT(2*pi/(trLength(n,m)*dk))

  IF (nk1 > nk) THEN
    nk = nk1
    PRINT*, 'Number of k points was changed, nk = ', nk
  END IF

  ALLOCATE(rka(nk))
  ALLOCATE(Enk(2,nhex,nk))
  ALLOCATE(Znk(2,2,nhex,nk))
  
! define k point array (1/A): -pi/T .. pi/T
  rkmax=pi/trLength(n,m)
  rkmin=-rkmax
  CALL linArray(nk,rkmin,rkmax,rka)
  dk=rka(2)-rka(1)

! metallicity
  IF (MOD(n-m,3) == 0) THEN
    metal = 1
  ELSE
    metal = 0
  END IF

! compute the array of cutting line indeces corresponding to certain ii transitions
  ALLOCATE(muii(4,nhex/2+1))
  CALL CutLineii(n,m,nhex,muii)
  OPEN(unit=22,file=TRIM(path)//'tube.cutline.'//outfile)
    IF ( metal == 1 ) THEN
        DO i = 1, nhex/2+1
            WRITE(22,*) i-1,i-1, muii(1,i), muii(2,i), muii(3,i), muii(4,i)
        END DO
    ELSE
        DO i = 1, nhex/2+1
            WRITE(22,*) i,i, muii(1,i), muii(2,i)
        END DO
    END IF
  CLOSE(unit=22)
  WRITE(*,*) 'cutting line info in tube.cutline.'//outfile

! compute energy bands En(k) (eV)
! NOTE: cutting lines are always taken from 1 to N
! the notations mu = 0..N-1 and mu = 1..N are consistent
! when 0 <-> N and other points are the same
  DO mu=1,nhex 
     DO k=1,nk
        rk=rka(k)
        CALL etbTubeBand(n,m,mu,rk,En,Zk)
        !CALL stbTubeBand(n,m,mu,rk,En,Zk)
        DO ii=1,2
           Enk(ii,mu,k)=En(ii)
           Znk(ii,1:2,mu,k)=Zk(ii,1:2)
        END DO
     END DO
  END DO

  OPEN(unit=22,file=TRIM(path)//'tube.Enk.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk),((Enk(ii,mu,k),ii=1,2),mu=1,nhex)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'electronic En(k) in tube.Enk.'//outfile

!Write wave function
  OPEN(unit=22,file=TRIM(path)//'tube.Znk.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk),((aimag(Znk(ii,1,mu,k)),ii=1,2),mu=1,nhex)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'electronic wavefunctions Zn(k) in tube.Znk.'//outfile

  DEALLOCATE(Znk)

!----------------------------------------------------------------------
!            electron density of states (states/atom/eV)
!----------------------------------------------------------------------

  WRITE (*,*) '====================================================='
  WRITE (*,*) '..electron density of states'
  WRITE (*,*) '..number of energies:',nee
  
  ALLOCATE(eEarray(nee))
  ALLOCATE(eDOS(nee))
  
  CALL linArray(nee,emaxe,emine,eEarray)
  de = eEarray(2) - eEarray(1)
  
  CALL tubeElDOS(n,m,nee,eEarray,eDOS)
      
  OPEN(unit=22,file=TRIM(path)//'tube.elecDOS.'//outfile)
  DO ie = 1, nee
     WRITE(22,1001) eDOS(ie),eEarray(ie)
  END DO
  CLOSE(unit=22)
  WRITE(*,*) 'DOS in tube.elecDOS.'//outfile

!  DEALLOCATE(eDOS)
!  DEALLOCATE(eEarray)

!----------------------------------------------------------------------
!             squared optical matrix elements (eV)
!----------------------------------------------------------------------
  WRITE (*,*) '====================================================='
  WRITE (*,*) '..Squared optical matrix elements (eV)'

  ALLOCATE(cDipole(3,nk,2,nhex,2,-1:1))
  DO n1 = 1,2
    DO n2 = 1,2
        DO mu1 = 1, nhex
            DO i=-1,1
                mu2 = mu1 + i ! mu1-1, mu1, mu1+1

                DO k = 1, nk
                rk = rka(k)

                    CALL tbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole(1:3,k,n1,mu1,n2,i)) ! (1/A)
                    !CALL stbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole(1:3,k,n1,mu1,n2,mu2)) ! (1/A)

                END DO

            END DO
         END DO
     END DO
  END DO

! Write output for matrix elements square: x and z components
! ------------------------------------------------------------
! x polarization      
! IF(INT(laser_theta) .NE. 0) THEN
!    WRITE(*,*) '..x polarization'
!
!    ALLOCATE(Px2k(nk,nhex))
!
! downward cutting line transitions from mu --> mu-1      
!    DO mu = 1, nhex
!       DO k = 1, nk
!          mmu = mu-1
!          IF (mmu.LT.1) mmu = nhex
!          Px2k(k,mu) = 3.81D0 * CDABS(cDipole(1,k,1,mu,2,mmu))**2
!          ! 3.81 = (hbar^2 / (2 m_e)) (eV-A**2)
!       END DO
!    END DO
!
!    OPEN(unit=22,file=TRIM(path)//'tube.Px2k_dn.'//outfile)
!    DO k = 1, nk
!       WRITE(22,1001) rka(k)/rka(nk),(Px2k(k,mu),mu=1,nhex)
!    END DO
!    CLOSE(unit=22)
!    WRITE(*,*) 'Px2_dn(k) in tube.Px2k_dn.'//outfile
!
! upward cutting line transitions from mu --> mu+1
!    DO mu = 1, nhex
!       DO k = 1, nk
!          mmu = mu+1
!          IF (mmu > nhex) mmu=1
!          Px2k(k,mu) = 3.81D0 * CDABS(cDipole(1,k,1,mu,2,mmu))**2
!       END DO
!    END DO
!
!    OPEN(unit=22,file=TRIM(path)//'tube.Px2k_up.'//outfile)
!    DO k = 1, nk
!       WRITE(22,1001) rka(k)/rka(nk),(Px2k(k,mu),mu=1,nhex)
!    END DO
!    CLOSE(unit=22)
!    WRITE(*,*) 'Px2_up(k) in tube.Px2k_up.'//outfile
!
!    DEALLOCATE(Px2k)
!
! END IF
!
! z polarization
! WRITE(*,*) '..z polarization'
! ALLOCATE(Pz2k(nk,nhex))
!
! cutting line transitions from mu --> mu 
! DO mu = 1, nhex
!    DO k = 1, nk
!       Pz2k(k,mu) = 3.81D0 * CDABS(cDipole(3,k,1,mu,2,mu))**2
!    END DO
! END DO
!
! OPEN(unit=22,file=TRIM(path)//'tube.Pz2k.'//outfile)
! DO k = 1, nk
!    WRITE(22,1001) rka(k)/rka(nk),(Pz2k(k,mu),mu=1,nhex)
! END DO
! CLOSE(unit=22)
! WRITE(*,*) 'Pz2(k) in tube.Pz2k.'//outfile
!
! DEALLOCATE(Pz2k)
            
!----------------------------------------------------------------------
!          real and imaginary part of dielectric function (dimensionless)
!              absorption coefficient (1/cm), conductivity (e^2/h)
!----------------------------------------------------------------------
  WRITE (*,*) '====================================================='
  WRITE (*,*) '..Real and Imaginary part of dielectric function'
  WRITE (*,*) '  absorption coefficient (1/cm), conductivity (e^2/h)'
  WRITE (*,*) '..Number of laser energies :',nhw_laser
  WRITE (*,*) '..Laser fwhm linewidth (eV):',laser_fwhm
  WRITE (*,*) '..Laser polarization angle (deg) :',laser_theta
  
  IF ( doping .ne. 0.) THEN
    Efermi = fermiLevel(n,m,Tempr,doping)
    WRITE (*,*) '..Fermi level calculated from doping', Efermi
  ENDIF

  ALLOCATE(hw_laser(nhw_laser))
  CALL linArray(nhw_laser,epmin,epmax,hw_laser)
  dep = hw_laser(2)-hw_laser(1)

! unit polarization vectors in complex notation (dimensionless)
  CALL polVector(laser_theta,epol)

  WRITE (thetastr, 360) INT(laser_theta/10.)

! PLASMON output ***********************************************************
  OPEN(unit=24,file=TRIM(path)//'tube.plasmon_intra.'//'theta'//TRIM(thetastr)//'.'//TRIM(outfile)//'.csv')
  OPEN(unit=25,file=TRIM(path)//'tube.epsZero.'//'theta'//TRIM(thetastr)//'.'//outfile)
  OPEN(unit=26,file=TRIM(path)//'tube.plasmon_inter.'//'theta'//TRIM(thetastr)//'.'//TRIM(outfile)//'.csv')
  OPEN(unit=27,file=TRIM(path)//'tube.chargeDens.'//outfile)

  WRITE (24,*) "Efermi ", "mu1 ", "mu2 ", "Ei ", "Ej ", "frequency ", "ratio ", "Absorpt(w) "
  WRITE (26,*) "Efermi ", "mu1 ", "mu2 ", "Ei ", "Ej ", "frequency ", "ratio ", "Absorpt(w) "
  WRITE(25,*) "Efermi ", "zero"

! ***********************************************************************
! FERMI LEVEL LOOP
! ***********************************************************************
! ---------------------------------------------------------------------
! cycle over fermi level position
! ---------------------------------------------------------------------
!  WRITE (*,*) '--------------------------------------------------------'
!  WRITE (*,*) '..begin DO loop over Fermi level in range -2.5..2.5 eV'
  WRITE (*,*) '..begin DO loop over Fermi level in range 1..2 eV'
  DO i = 1, 5
  WRITE (*,*) '--------------------------------------------------------'

!  Efermi = -2.5 + (i-1) * 0.1D0
  Efermi = 1. + (i-1) * 0.25D0
  WRITE (*,*) '..Fermi level: ', Efermi
  WRITE (fermistr, 350) Efermi
  CALL EXECUTE_COMMAND_LINE( 'mkdir -p tube'//TRIM(outfile)//'/pol_'//TRIM(thetastr)//'_fl_'//TRIM(ADJUSTL(fermistr)) )
  path = './tube'//TRIM(outfile)//'/pol_'//TRIM(thetastr)//'_fl_'//TRIM(ADJUSTL(fermistr))//'/'
  WRITE (*,*) '--------------------------------------------------------'

  CALL ChargeDensity(nhex, nk, nee, Enk, Tempr, Efermi, eEarray, eDOS, charge)
  WRITE(27,1001) Efermi, charge
  WRITE (*,*) 'charge density in tube.chargeDens.'//outfile
  WRITE (*,*) '--------------------------------------------------------'

! ======================================================================
! ========================== permittivity ==============================
! real and imaginary parts of dielectric permittivity
! ======================================================================
  WRITE (*,*) '--------------------------------------------------------'
! interband part ==============
  ALLOCATE(eps1(nhw_laser))
  ALLOCATE(eps2(nhw_laser))
  ALLOCATE(eps1Part(2,nhex,2,nhex,nhw_laser))
  ALLOCATE(eps2Part(2,nhex,2,nhex,nhw_laser))
  CALL DielPermittivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,ebg,laser_fwhm,nhw_laser,hw_laser,&
  eps1Part,eps2Part,eps1,eps2)
! DRUDE part ==================
!  ALLOCATE(eps1Dr(nhw_laser))
!  ALLOCATE(eps2Dr(nhw_laser))
!  CALL DielPermittivityDr(n,m,nhex,nk,rka,Enk,Tempr,Efermi,epol,ebg,laser_fwhm,nhw_laser,hw_laser,eps1Dr,eps2Dr)

! plot eps1(hw) *********************************
  OPEN(unit=22,file=TRIM(path)//'tube.eps1.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eps1(ie) !, eps1Dr(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'real part of dielectric function in tube.eps1.'//outfile

! plot eps2(hw) ********************************
  OPEN(unit=22,file=TRIM(path)//'tube.eps2.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eps2(ie) !, eps2Dr(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of dielectric function in tube.eps2.'//outfile
!  DEALLOCATE(eps1Dr)
!  DEALLOCATE(eps2Dr)

! =========== KRAMERS-KRONIG =====================
!  ALLOCATE(eps1kk(nhw_laser))
!  ALLOCATE(eps2kk(nhw_laser))
!  CALL DielPermittivityKrKr(nhw_laser,ebg,hw_laser,eps1,eps2,eps1kk,eps2kk)   !Kramers-Kronig
! plot eps1(hw) (Kramers-Kronig) ****************
! OPEN(unit=22,file=TRIM(path)//'tube.eps1kk.'//outfile)
! DO ie = 1, nhw_laser
!    WRITE(22,1001) hw_laser(ie),eps1kk(ie)
! ENDDO
! CLOSE(unit=22)
! WRITE(*,*) 'Kramers-Kronig real part of dielectric function in tube.eps1kk.'//outfile
!
! plot eps2(hw) (Kramers-Kronig) ***************
! OPEN(unit=22,file=TRIM(path)//'tube.eps2kk.'//outfile)
! DO ie = 1, nhw_laser
!    WRITE(22,1001) hw_laser(ie),eps2kk(ie)
! ENDDO
! CLOSE(unit=22)
! WRITE(*,*) 'Kramers-Kronig imaginary part of dielectric function in tube.eps2kk.'//outfile
!  DEALLOCATE(eps1kk)
!  DEALLOCATE(eps2kk)
!! =======================================================================
!! ======================= absorption coefficient ========================
!! alpha
!! =======================================================================
!  WRITE (*,*) '--------------------------------------------------------'
!   ALLOCATE(alpha(nhw_laser))
!  CALL imagDielAlpha(nhw_laser,hw_laser,eps2,refrac,alpha)
!
!! plot absorption alpha(hw) ********************
!  OPEN(unit=22,file=TRIM(path)//'tube.alpha.'//outfile)
!  DO ie = 1,nhw_laser
!     WRITE(22,1001) hw_laser(ie),alpha(ie)
!  ENDDO
!  CLOSE(unit=22)
!  WRITE(*,*) 'alpha(hw) in tube.alpha.'//outfile
!   DEALLOCATE(alpha)

! =======================================================================
! ======================= Electron energy loss spectra  =================
! eels spectra
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  ALLOCATE(eelspec(nhw_laser))
  CALL EELS(nhw_laser,eps1,eps2,eelspec)

! plot eels(hw) *******************************
  OPEN(unit=22,file=TRIM(path)//'tube.eels.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eelspec(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'eels in tube.eels.'//outfile
  DEALLOCATE(eelspec)
! =======================================================================
! ========================== conductivity ===============================
! real and imaginary parts part
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
! interband part =============
  ALLOCATE(sigm1(nhw_laser))
  ALLOCATE(sigm2(nhw_laser))
  ALLOCATE(sigm1Part(2,nhex,2,nhex,nhw_laser))
  ALLOCATE(sigm2Part(2,nhex,2,nhex,nhw_laser))
  CALL DynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,&
  sigm1Part,sigm2Part,sigm1,sigm2)
! DRUDE part =================
!  ALLOCATE(sigm1Dr(nhw_laser))
!  ALLOCATE(sigm2Dr(nhw_laser))
!  CALL DynConductivityDr(n,m,nhex,nk,rka,Enk,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm1Dr,sigm2Dr)

! plot sigm1(hw) ******************************
  OPEN(unit=22,file=TRIM(path)//'tube.sigm1.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), sigm1(ie)/(e2/h) !, sigm1Dr(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'real part of conductivity in tube.sigm1.'//outfile

! plot sigm2(hw) *******************************
  OPEN(unit=22,file=TRIM(path)//'tube.sigm2.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), sigm2(ie)/(e2/h) !, sigm2Dr(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of conductivity in tube.sigm2.'//outfile

! ============= Intra- and interband contributions to conductivity ===============
! imaginary part of intraband conductivity
! =======================================================================
! WRITE (*,*) '--------------------------------------------------------'
!  ALLOCATE(sigm1_intra(nhw_laser))
!  ALLOCATE(sigm2_intra(nhw_laser))
!  ALLOCATE(sigm1_inter(nhw_laser))
!  ALLOCATE(sigm2_inter(nhw_laser))
! CALL DynConductivityIntra(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm1_intra,sigm2_intra)
!
! plot sigm2_intra(hw) *******************************
! OPEN(unit=22,file=TRIM(path)//'tube.sigm2_intra.'//outfile)
! DO ie = 1, nhw_laser
!    WRITE(22,1001) hw_laser(ie), sigm2_intra(ie)/(e2/h)
! ENDDO
! CLOSE(unit=22)
! WRITE(*,*) 'imaginary part of intraband conductivity in tube.sigm2_intra.'//outfile
!
! imaginary part of interband conductivity
! =======================================================================
! CALL DynConductivityInter(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm1_inter,sigm2_inter)
!
! plot sigm2_inter(hw) *******************************
! OPEN(unit=22,file=TRIM(path)//'tube.sigm2_inter.'//outfile)
! DO ie = 1, nhw_laser
!    WRITE(22,1001) hw_laser(ie), sigm2_inter(ie)/(e2/h)
! ENDDO
! CLOSE(unit=22)
! WRITE(*,*) 'imaginary part of interband conductivity in tube.sigm2_inter.'//outfile
!  DEALLOCATE(sigm1_intra)
!  DEALLOCATE(sigm2_intra)
!  DEALLOCATE(sigm1_inter)
!  DEALLOCATE(sigm2_inter)

! =======================================================================
! ============================ absorption ===============================
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  ALLOCATE(absorpt(nhw_laser))
  CALL Absorption(nhw_laser,eps1,eps2,sigm1,sigm2,absorpt)

! plot absorpt(hw) *******************************
  OPEN(unit=22,file=TRIM(path)//'tube.absorpt.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), absorpt(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'absorption in tube.absorpt.'//outfile
  DEALLOCATE(sigm1, sigm2)
!! =======================================================================
!! ======================= looking for plasmon ===========================
  WRITE(*,*) '-------------------------------------------------------'
  ALLOCATE(absorptPart(2,nhex/2,2,nhex/2,nhw_laser))
  WRITE(*,*) '..LOOP over all possible transitions'

  DO n1 = 1,2
    DO n2 = 1,2
        DO mu1 = 1, nhex/2
            DO mu2 = 1, nhex/2

            CALL Absorption(nhw_laser,eps1,eps2,sigm1Part(n1,mu1,n2,mu2,:),sigm2Part(n1,mu1,n2,mu2,:),absorptPart(n1,mu1,n2,mu2,:))

            END DO
         END DO
      END DO
  END DO

! CONTRIBUTIONS output *****************************************************
! OPEN(unit=122,file=TRIM(path)//'tube.eps1Part.'//outfile)
! OPEN(unit=123,file=TRIM(path)//'tube.eps2Part.'//outfile)
! OPEN(unit=124,file=TRIM(path)//'tube.sigm1Part.'//outfile)
! OPEN(unit=125,file=TRIM(path)//'tube.sigm2Part.'//outfile)
! OPEN(unit=126,file=TRIM(path)//'tube.absorptPart.'//outfile)
! **************************************************************************

  plasmonExists = 0
  zeroNumber = 0

! look for zeros in real part of dielectric function
! plasmon is present if there is zero in dielectric function
  DO ie = 1, nhw_laser-1
      IF ( eps1(ie) * eps1(ie+1) < 0.D0 .and. eps1(ie) < 0.D0 ) THEN
        zeroNumber = zeroNumber + 1
      END IF
  END DO

  IF (zeroNumber /= 0) THEN
    ALLOCATE(plasmonPosition(zeroNumber))
    ALLOCATE(eps1Zeros(zeroNumber))

    plasmonExists = 1
    plasmonPosition = 0
    eps1Zeros = 0.0

      l = 1
      place_max = 0
      DO ie = 1, nhw_laser-1
          IF ( eps1(ie) * eps1(ie+1) < 0.D0 .and. eps1(ie) < 0.D0 ) THEN
             ! save the information about zeros of dielectric function
             eps1Zeros(l) = (hw_laser(ie)+hw_laser(ie+1))/2
             WRITE(25,1001) Efermi, eps1Zeros(l)

             ! look for plasmon position in proximity of eps=0
             b1 = ie
             b2 = ie
             DO WHILE ( b1 < nhw_laser-1 )
                b1 = b1-1
                b2 = b2+1
                place_max = b1 + MAXLOC(absorpt(b1:b2),1) - 1
                IF (place_max /= b1 .and. place_max /= b2) EXIT
             END DO
             plasmonPosition(l) = place_max

             l = l+1
          END IF
      END DO

  END IF


! consider only one K point
  DO n1 = 1,2
    DO n2 = 1,2
        DO mu1 = 1, nhex/2
            DO mu2 = 1, nhex/2

! look for zeros in real part of ij contribution to dielectric function to choose important contributions
            zeroExists = 0
            DO ie = 1, nhw_laser-1
                IF ( eps1Part(n1,mu1,n2,mu2,ie) * eps1Part(n1,mu1,n2,mu2,ie+1) < 0.D0 ) THEN
                    zeroExists = 1
                END IF
            END DO

            loc_mu1 = MINLOC(ABS(muii - mu1))
            loc_mu2 = MINLOC(ABS(muii - mu2))

! check if zeros are present in both n1,mu1,n2,mu2 contribution to dielectric function
! and also in the total dielectric function
            IF ( zeroExists == 1 .and. plasmonExists == 1 ) THEN
            DO l = 1,zeroNumber
            ratio = absorptPart(n1,mu1,n2,mu2,plasmonPosition(l))/MAXVAL(absorptPart(:,:,:,:,plasmonPosition(l)))
            ! check if contribution is essential for plasmon
               IF ( ratio > 1.D-1 ) THEN
                    ! the output differs for metallic and semiconducting nanotubes:
                    ! for SC we start counting of Ei from 1
                    ! for M we start counting of Ei from 0
                        IF ( n1 == n2 ) THEN
                            WRITE(24,3003) Efermi, mu1, mu2, loc_mu1(2) - metal, &
                            loc_mu2(2) - metal, hw_laser(plasmonPosition(l)), &
                            ratio, absorpt(plasmonPosition(l))/(e2/h)
                        ELSE
                            WRITE(26,3003) Efermi, mu1, mu2, loc_mu1(2) - metal, &
                            loc_mu2(2) - metal, hw_laser(plasmonPosition(l)), &
                            ratio, absorpt(plasmonPosition(l))/(e2/h)
                        END IF
                END IF

            END DO
            END IF
! ------------------------------------------------------------------------
! HERE the output STARTS for contributions *******************************
! ------------------------------------------------------------------------
!              IF ( zeroExists == 1 ) THEN
!                  DO ie = 1, nhw_laser
!                      WRITE(122,1001) n1,n2, mu1,mu2, loc_mu1(2) - metal, loc_mu2(2) - metal,&
!                      hw_laser(ie), eps1Part(n1,mu1,n2,mu2,ie)
!                      WRITE(123,1001) n1,n2, mu1,mu2, loc_mu1(2) - metal, loc_mu2(2) - metal,&
!                      hw_laser(ie), eps2Part(n1,mu1,n2,mu2,ie)
!                      WRITE(124,1001) n1,n2, mu1,mu2, loc_mu1(2) - metal, loc_mu2(2) - metal,&
!                      hw_laser(ie), sigm1Part(n1,mu1,n2,mu2,ie)/(e2/h)
!                      WRITE(125,1001) n1,n2, mu1,mu2, loc_mu1(2) - metal, loc_mu2(2) - metal,&
!                      hw_laser(ie), sigm2Part(n1,mu1,n2,mu2,ie)/(e2/h)
!                      WRITE(126,1001) n1,n2, mu1,mu2, loc_mu1(2) - metal, loc_mu2(2) - metal,&
!                      hw_laser(ie), absorptPart(n1,mu1,n2,mu2,ie)/(e2/h)
!                  ENDDO
!                  WRITE(122,1001) ' '
!                  WRITE(123,1001) ' '
!                  WRITE(124,1001) ' '
!                  WRITE(125,1001) ' '
!                  WRITE(126,1001) ' '
!              END IF

! ------------------------------------------------------------------------
! HERE the output ENDS for contributions *********************************
! ------------------------------------------------------------------------

            END DO
        END DO
    END DO
  END DO

  IF (ALLOCATED(plasmonPosition) .EQV. .TRUE.) DEALLOCATE(plasmonPosition)
  IF (ALLOCATED(eps1Zeros) .EQV. .TRUE.) DEALLOCATE(eps1Zeros)

  DEALLOCATE(absorpt,absorptPart)
  DEALLOCATE(eps1,eps2)
  DEALLOCATE(eps1Part,eps2Part)
  DEALLOCATE(sigm1Part,sigm2Part)

!  CLOSE(unit=122)
!  CLOSE(unit=123)
!  CLOSE(unit=124)
!  CLOSE(unit=125)
!  CLOSE(unit=126)

  WRITE(*,*) '..End of LOOP '

  WRITE(*,*) '-------------------------------------------------------'

!  ! plot eps1Part(hw) *******************************
!  WRITE(*,*) 'different contributions in eps1 in tube.eps1Part.'//outfile
!  ! plot eps2Part(hw) *******************************
!  WRITE(*,*) 'different contributions in eps2 in tube.eps2Part.'//outfile
!  ! plot sigm1Part(hw) *******************************
!  WRITE(*,*) 'different contributions in sigm1 in tube.sigm1Part.'//outfile
!  ! plot sigm2Part(hw) *******************************
!  WRITE(*,*) 'different contributions in sigm2 in tube.sigm2Part.'//outfile
!  ! plot absorptPart(hw) *******************************
!  WRITE(*,*) 'different contributions in absorption in tube.absorptPart.'//outfile

  END DO

  WRITE (*,*) '--------------------------------------------------------'
  WRITE (*,*) '..end of DO loop over Fermi level'
  WRITE (*,*) '--------------------------------------------------------'

! ***********************************************************************
! END OF FERMI LEVEL LOOP
! ***********************************************************************

  CLOSE(unit=24)
  CLOSE(unit=25)
  CLOSE(unit=26)
  CLOSE(unit=27)

  WRITE(*,*) 'plasmon information in tube.plasmon_intra.'//TRIM(outfile)//'.csv ', 'and ', &
  'tube.plasmon_inter.'//TRIM(outfile)//'.csv '
  WRITE (*,*) '--------------------------------------------------------'


! ============= part of code to check the calculations ================================
! *************** please, remove it later *********************************************
! ---------------------------------------------------------
! check convergence of the result according to dk chosen
! ---------------------------------------------------------
!  DEALLOCATE(rka)
!  DEALLOCATE(Enk)
!  DEALLOCATE(Znk)
!  DEALLOCATE(cDipole)
!
!  ALLOCATE(eps2aii(10,nhw_laser))
!
!  kCoef = (/1., 2., 3., 4., 5., 6., 7., 8., 9., 10./)
!
!  DO i=1,10
!
!      dk = laser_fwhm/hbarvfermi*1.D0/kCoef(i)
!      nk = INT(pi/(trLength(n,m)*dk))
!
!      PRINT*, 'Coef = ', kCoef(i), 'Number of k points, nk = ', nk
!
!  ALLOCATE(rka(nk))
!  ALLOCATE(Enk(2,nhex,nk))
!  ALLOCATE(Znk(2,2,nhex,nk))
!  ALLOCATE(cDipole(3,nk,2,nhex,2,nhex))
!
!! define k point array (1/A): -pi/T .. pi/T
!  rkmax=pi/trLength(n,m)
!  rkmin=-rkmax
!  CALL linArray(nk,rkmin,rkmax,rka)
!  dk=rka(2)-rka(1)
!
!! compute energy bands En(k) (eV)
!  DO mu=1,nhex
!     DO k=1,nk
!        rk=rka(k)
!        CALL etbTubeBand(n,m,mu,rk,En,Zk)
!        DO ii=1,2
!           Enk(ii,mu,k)=En(ii)
!           Znk(ii,1:2,mu,k)=Zk(ii,1:2)
!        END DO
!     END DO
!  END DO
!
!  DO n1 = 1,2
!    DO n2 = 1,2
!        DO mu1 = 1, nhex
!            DO mu2 = 1, nhex
!
!                DO k = 1, nk
!                rk = rka(k)
!
!                    CALL tbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole(1:3,k,n1,mu1,n2,mu2)) ! (1/A)
!
!                END DO
!
!            END DO
!         END DO
!     END DO
!  END DO
!
!  CALL imagDielEn(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,eps2)
!
!  eps2aii(i,1:nhw_laser) = eps2(1:nhw_laser)
!
!  DEALLOCATE(rka)
!  DEALLOCATE(Enk)
!  DEALLOCATE(Znk)
!  DEALLOCATE(cDipole)
!
!  END DO
!
!  DO i=1,9
!    divergence(i) = 0.0
!    DO ie = 1, nhw_laser
!        divergence(i) = divergence(i) + (eps2aii(i,ie) - eps2aii(10,ie))**2
!    END DO
!    divergence(i) = SQRT(divergence(i))/nhw_laser
!    maxAbsDif(i) = MAXVAL(ABS(eps2aii(i,1:nhw_laser) - eps2aii(10,1:nhw_laser)))
!  END DO
!
!  OPEN(unit=22,file='tube.divergence.'//outfile)
!  DO i=1,9
!     WRITE(22,1001) kCoef(i), divergence(i), maxAbsDif(i)
!  ENDDO
!  CLOSE(unit=22)
!  WRITE(*,*) 'divergence in file tube.divergence.'//outfile
!
!  DEALLOCATE(eps2aii)

! ============= part of code to check the calculations ================================
! ************************************* END *******************************************

! WRITE(*,*) '====================================================='
! WRITE(*,*) 'Fermi level: ', 'from input', Efermi, 'from doping', fermiLevel(n,m,Tempr,doping)

  DEALLOCATE(cDipole)
  DEALLOCATE(Enk,rka)
  DEALLOCATE(hw_laser)

  WRITE(*,*) 'Program successfully completed'

! TEST cycle to find chyralities **********************************
! OPEN(unit=22,file='tube.chyralities')
!DO m = 0,30
!   DO n = m,30
!
!    IF ( n == 0 ) CYCLE
!    diameter = tubeDiam(n,m)/10 !nm
!    IF ( diameter .ge. 0.5D0 .and. diameter .le. 2.D0 ) THEN
!       WRITE(22,*) n,m
!    END IF
!
!   END DO
!END DO
!CLOSE(unit=22)
!WRITE(*,*) 'tube chyralities in tube.chyralities '


! ****************************************************************

!----------------------------------------------------------------------
!                          format statements
!----------------------------------------------------------------------
      WRITE(*,*) '====================================================='
      STOP

1001  FORMAT(901(1X,G13.5))
1002  FORMAT(4F8.4)
2002  FORMAT(1A20)
350   FORMAT(F4.1)
360   FORMAT(I1)
3003  FORMAT(F8.2,4I8,6F8.3)

END PROGRAM cntabsorpt
!*******************************************************************************
