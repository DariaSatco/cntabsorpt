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
! Latest Vers. : 2020.06.18
!-------------------------------------------------------------------------------
! Required files :
! - globvar         -- global variables
! - tubepar.f90     -- nanotube parameters
! - tubestruct.f90  -- nanotube structure libraries
! - libMath.f90     -- mathematical libraries
! - libswntElec.f90 -- electronic/excitonic states
! - libswntOpt.f90  -- optical matrix element, dielectric function, conductivity
! - libswntSTB.f90  -- simple tight-binding for energy bands and dipole matrix
! - libDrudeOpt.f90 -- Drude conductivity
!*******************************************************************************
!*******************************************************************************
PROGRAM cntabsorpt
!===============================================================================
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
  INTEGER                :: k, mu, ii, ie, nee
  INTEGER                :: nhex
  INTEGER                :: nk
  INTEGER                :: nhw_laser
  REAL(8)                :: Tempr, rkT, doping, emine, emaxe
  REAL(8)                :: refrac
  REAL(8)                :: epmin, epmax
  REAL(8)                :: laser_theta, laser_fwhm
  REAL(8)                :: ebg
  REAL(8)                :: Efermi
!-------------------------------------------------------------------------------
! variables for electronic states
  REAL(8)                :: rkmin, rkmax, dk, rk, de
  INTEGER                :: nk1
  REAL(8), ALLOCATABLE   :: rka(:)           !(nk)
  REAL(8), ALLOCATABLE   :: Enk(:,:,:)       !(2,nhex,nk)
  COMPLEX(8), ALLOCATABLE:: Znk(:,:,:,:)
  REAL(8), DIMENSION(2)  :: En               !(unit eV)
  COMPLEX(8)             :: Zk(2,2)
  
  REAL(8), ALLOCATABLE   :: eEarray(:)       !(nee)
  REAL(8), ALLOCATABLE   :: eDOS(:)          !(nee)
  REAL(8), ALLOCATABLE   :: JDOSinter(:), JDOSintraC(:), JDOSintraV(:)      !(nee)
!-------------------------------------------------------------------------------
! variables for optical properties
  INTEGER                 :: n1, mu1, n2, mu2
  REAL(8)                 :: dep
  REAL(8)                 :: epol(3)

  REAL(8), ALLOCATABLE    :: hw_laser(:)     !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps2(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps1(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm1(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eelspec(:)      !(nhw_laser)

  REAL(8), ALLOCATABLE    :: absorpt(:)      !(nhw_laser)
  COMPLEX(8), ALLOCATABLE :: cDipole(:,:,:,:,:,:) !(3,nk,2,nhex,2,nhex)
  REAL(8), ALLOCATABLE    :: absorptPart(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: eps1Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: eps2Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: sigm1Part(:,:,:,:,:)
  REAL(8), ALLOCATABLE    :: sigm2Part(:,:,:,:,:)

  REAL(8)                 :: charge
  REAL(8)                 :: capacitance
!-------------------------------------------------------------------------------
! variables for input and output files 
  CHARACTER(40)           :: infile, outfile, path
  CHARACTER(5)            :: fermistr, thetastr

! working variales
  INTEGER                 :: zeroNumber
  REAL(8), ALLOCATABLE    :: eps1Zeros(:)
  INTEGER, ALLOCATABLE    :: plasmonPosition(:)
  REAL(8)                 :: ratio
  INTEGER, ALLOCATABLE    :: muii(:,:)
  INTEGER                 :: place_max, zeroExists, metal, loc_mu1(2), loc_mu2(2), plasmonExists
  INTEGER                 :: i, l, b1, b2


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

!----------------------------------------------------------------------
!            joint density of states (states/atom/eV)
!----------------------------------------------------------------------
  ALLOCATE(JDOSinter(nee))
  ALLOCATE(JDOSintraC(nee))
  ALLOCATE(JDOSintraV(nee))

  CALL tubeJDOS(n,m,nhex,nk,Enk,rka,nee,eEarray,JDOSinter,JDOSintraV, JDOSintraC)

  OPEN(unit=22,file=TRIM(path)//'tube.JDOS.'//outfile)
  DO ie = 1, nee
     WRITE(22,1001) eEarray(ie), JDOSinter(ie), JDOSintraV(ie), JDOSintraC(ie)
  END DO
  CLOSE(unit=22)
  WRITE(*,*) 'JDOS in tube.JDOS.'//outfile

!  DEALLOCATE(eDOS)
!  DEALLOCATE(eEarray)

!----------------------------------------------------------------------
!             squared optical matrix elements (eV)
!----------------------------------------------------------------------
  WRITE (*,*) '====================================================='
  WRITE (*,*) '..Squared optical matrix elements (eV)'

  ALLOCATE(cDipole(3,nk,2,nhex,2,-1:1))
  cDipole = 0.D0
  DO n1 = 1,2
    DO n2 = 1,2
        DO mu1 = 1, nhex
            DO i=-1,1
                mu2 = mu1 + i ! mu1-1, mu1, mu1+1
                IF(mu2 < 1 .OR. mu2 > nhex) CYCLE

                DO k = 1, nk
                rk = rka(k)

                    CALL tbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole(1:3,k,n1,mu1,n2,i)) ! (1/A)

                END DO

            END DO
         END DO
     END DO
  END DO
            
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
  OPEN(unit=28,file=TRIM(path)//'tube.quantCapac.'//outfile)

  WRITE (24,*) "Efermi ", "mu1 ", "mu2 ", "Ei ", "Ej ", "frequency ", "ratio ", "Absorpt(w) "
  WRITE (26,*) "Efermi ", "mu1 ", "mu2 ", "Ei ", "Ej ", "frequency ", "ratio ", "Absorpt(w) "
  WRITE(25,*) "Efermi ", "zero"

! ***********************************************************************
! FERMI LEVEL LOOP
! ***********************************************************************
! ---------------------------------------------------------------------
! cycle over fermi level position
! ---------------------------------------------------------------------
!  Efstart = 0.
!  Efend = 2.
!  dEf = 0.1
!  NEf = int((Efend - Efstart)/dEf) + 1
!  WRITE (*,*) '--------------------------------------------------------'
!  WRITE (*,*) '..begin DO loop over Fermi level in range ', Efstart, '..', Efend, ' eV'
!  DO i = 1, Nef
!  WRITE (*,*) '--------------------------------------------------------'
!
!  Efermi = Efstart + (i-1) * dEf
  WRITE (*,*) '..Fermi level: ', Efermi
  WRITE (fermistr, 350) Efermi
  CALL EXECUTE_COMMAND_LINE( 'mkdir -p tube'//TRIM(outfile)//'/pol_'//TRIM(thetastr)//'_fl_'//TRIM(ADJUSTL(fermistr)) )
  path = './tube'//TRIM(outfile)//'/pol_'//TRIM(thetastr)//'_fl_'//TRIM(ADJUSTL(fermistr))//'/'
  WRITE (*,*) '--------------------------------------------------------'

! calculate charge density =======================
  CALL ChargeDensity(nee, Tempr, Efermi, eEarray, eDOS, charge)
  WRITE(27,1001) Efermi, charge
  WRITE (*,*) 'charge density in tube.chargeDens.'//outfile
  WRITE (*,*) '--------------------------------------------------------'

! quantum capacitance ============================
  CALL QuantumCapacitance(nee, Tempr, Efermi, eEarray, eDOS, capacitance)
  WRITE(28,1001) Efermi, capacitance
  WRITE (*,*) 'quantum capacitance in tube.quantCapac.'//outfile
  WRITE (*,*) '--------------------------------------------------------'

!  END DO
!
!  WRITE (*,*) '--------------------------------------------------------'
!  WRITE (*,*) '..end of DO loop over Fermi level'
!  WRITE (*,*) '--------------------------------------------------------'

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


! =======================================================================
! ======================= looking for plasmon ===========================
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

! CONTRIBUTIONS output - uncomment if you want to save different ij plasmonic contributions ******
! OPEN(unit=122,file=TRIM(path)//'tube.eps1Part.'//outfile)
! OPEN(unit=123,file=TRIM(path)//'tube.eps2Part.'//outfile)
! OPEN(unit=124,file=TRIM(path)//'tube.sigm1Part.'//outfile)
! OPEN(unit=125,file=TRIM(path)//'tube.sigm2Part.'//outfile)
! OPEN(unit=126,file=TRIM(path)//'tube.absorptPart.'//outfile)
! *************************************************************************************************

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

!  END DO
!
!  WRITE (*,*) '--------------------------------------------------------'
!  WRITE (*,*) '..end of DO loop over Fermi level'
!  WRITE (*,*) '--------------------------------------------------------'

! ***********************************************************************
! END OF FERMI LEVEL LOOP
! ***********************************************************************

  CLOSE(unit=24)
  CLOSE(unit=25)
  CLOSE(unit=26)
  CLOSE(unit=27)
  CLOSE(unit=28)

  WRITE(*,*) 'plasmon information in tube.plasmon_intra.'//TRIM(outfile)//'.csv ', 'and ', &
  'tube.plasmon_inter.'//TRIM(outfile)//'.csv '
  WRITE (*,*) '--------------------------------------------------------'

  DEALLOCATE(cDipole)
  DEALLOCATE(Enk,rka)
  DEALLOCATE(hw_laser)

  WRITE(*,*) 'Program successfully completed'



! ****************************************************************
!----------------------------------------------------------------------
!                          format statements
!----------------------------------------------------------------------
      WRITE(*,*) '====================================================='
      STOP

1001  FORMAT(901(1X,G13.5))
!1002  FORMAT(4F8.4)
2002  FORMAT(1A20)
350   FORMAT(F4.1)
360   FORMAT(I1)
3003  FORMAT(F8.2,4I8,6F8.3)

END PROGRAM cntabsorpt
!*******************************************************************************
