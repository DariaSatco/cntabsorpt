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
  REAL(8), PARAMETER     :: hbarc = 1.97D-5  !(eV cm)
  REAL(8), PARAMETER     :: hbarvfermi = 6.582119 !(eV-A) !hbar*vfermi, vfermi = 10**6 m/s
!-------------------------------------------------------------------------------
! variables for calling tube structure function
  INTEGER                :: n, m, nHexagon
  REAL(8)                :: trLength
  REAL(8)                :: fermiLevel
!-------------------------------------------------------------------------------
! calculation parameter input
  INTEGER                :: iprofile, imenu, ichange, iedit, ierr
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
  REAL(8), ALLOCATABLE   :: EiiOpt(:), rkiiOpt(:)
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
  REAL(8), ALLOCATABLE    :: eps2a(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps2a1(:)       !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps1a(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eps1a1(:)       !(nhw_laser)
  REAL(8), ALLOCATABLE    :: eelspec(:)      !(nhw_laser)
  REAL(8), ALLOCATABLE    :: alpha(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm1(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2(:)        !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2_intra(:)  !(nhw_laser)
  REAL(8), ALLOCATABLE    :: sigm2_inter(:)  !(nhw_laser)
  REAL(8), ALLOCATABLE    :: absorpt(:)      !(nhw_laser)
  COMPLEX(8), ALLOCATABLE :: cDipole(:,:,:,:,:,:) !(3,nk,2,nhex,2,nhex)
!-------------------------------------------------------------------------------
! variables for input and output files 
  CHARACTER(40)           :: infile, outfile

!temp variables
  REAL(8), ALLOCATABLE   :: plotfuncmaxtest(:), plotfuncmintest(:), kaxis(:)
  REAL(8), ALLOCATABLE   :: ss0(:,:,:,:)
  REAL(8), ALLOCATABLE   :: difFermiDist(:,:,:,:)
  REAL(8), ALLOCATABLE   :: matrElementSq(:,:,:,:)
  REAL(8), ALLOCATABLE   :: diracAvgFunc(:,:,:,:)
  INTEGER                :: prefer_position(4), max_position(4), min_position(4)

  REAL                   :: divergence(9), kCoef(10)
  INTEGER                :: i
  REAL(8), ALLOCATABLE   :: eps2aii(:,:)        !(nhw_laser)
!----------------------------------------------------------------------
!                              Main Menu
!----------------------------------------------------------------------
  iprofile = 0          
  WRITE (*,*) '====================================================='
  WRITE (*,*) '          Main Menu for coherent.f90'
1 WRITE (*,*) '====================================================='
  WRITE (*,*) 'Select an option:'
  WRITE (*,*) '(1)  Start with default data set'
  WRITE (*,*) '(2)  Read and edit an old data set'
  WRITE (*,*) '(3)  Delete old data set but save log & param files'
  WRITE (*,*) '(4)  Totally delete old data set'
  WRITE (*,*) '(5)  Stop'
  WRITE (*,*) '====================================================='
  imenu = 1

  IF (iprofile /= 1) READ (*,*) imenu
  IF (imenu < 1 .OR. imenu > 5 ) GOTO 1
  IF (imenu == 5) STOP
  
!----------------------------------------------------------------------
!                 menu option: delete an old data set
!----------------------------------------------------------------------
  IF (imenu == 3 .OR. imenu == 4) THEN
     WRITE (*,*) 'Enter suffix of data set to be deleted:'
     READ  (*,2002) infile
     outfile = infile
      
     IF (imenu == 4) THEN      
        OPEN(unit=22,file='tube.param.'//outfile)
        CLOSE(unit=22,status='delete')
      
        OPEN(unit=22,file='tube.log.'//outfile)
        CLOSE(unit=22,status='delete')      
     END IF
     
     OPEN(unit=22,file='tube.Enk.xyy.'//outfile)
     CLOSE(unit=22,status='delete')                                         
     
     OPEN(unit=22,file='tube.elecDOS.xyy.'//outfile)
     CLOSE(unit=22,status='delete')      
     
     OPEN(unit=22,file='tube.Px2k_up.xyy.'//outfile)
     CLOSE(unit=22,status='delete') 
     
     OPEN(unit=22,file='tube.Px2k_dn.xyy.'//outfile)
     CLOSE(unit=22,status='delete') 
     
     OPEN(unit=22,file='tube.Pz2k.xyy.'//outfile)
     CLOSE(unit=22,status='delete')   
     
     OPEN(unit=22,file='tube.eps2.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.eps2kk.xyy.'//outfile)
     CLOSE(unit=22,status='delete')     
     
     OPEN(unit=22,file='tube.eps1.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.eps1kk.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.alpha.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.absorpt.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.sigm1.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.sigm2.xyy.'//outfile)
     CLOSE(unit=22,status='delete')

     OPEN(unit=22,file='tube.sigm2_intra.xyy.'//outfile)
     CLOSE(unit=22,status='delete')
     
     OPEN(unit=22,file='tube.sigm2_inter.xyy.'//outfile)
     CLOSE(unit=22,status='delete')
     
     GOTO 1
  END IF
      
!----------------------------------------------------------------------
!                    menu option: default data set
!----------------------------------------------------------------------
  IF (imenu == 1) THEN
     
     infile  = 'default'
     outfile = 'default'        
     
     Tempr   = 300.D0
     
     n = 11
     m = 0
     
     nk      = 100

     doping  = 0.D0
     Efermi  = 0.0D0
     
     refrac      = 1.3D0
     ebg         = 2.4D0

     nhw_laser   = 501
     epmin       = .5D0
     epmax       = 5.D0
     laser_fwhm  =.15D0

     nee     = 501                                 
     emine   = -10.D0
     emaxe   = 15.D0
     
     laser_theta = 0.D0
          
  END IF
      
!----------------------------------------------------------------------
!              menu option: read and edit old data set
!----------------------------------------------------------------------
  IF (imenu == 2) THEN
     WRITE (*,*) 'Enter suffix of old data set (1a20):'
     READ  (*,2002) infile
     
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
      
  END IF
      
  ichange = 0      
2 WRITE (*,*) '====================================================='
  CALL printTubeClass(n,m,6)
  WRITE (*,*) '====================================================='
  WRITE (*,*) ' Here is the input data:'
  WRITE (*,*) '====================================================='
  WRITE (*,*) '(0)  OK. Start simulation'
  WRITE (*,*) '-----------------------------------------------------'
  WRITE (*,*) '(1)  Temperature (deg K)            :',Tempr
  WRITE (*,*) '(2)  Chiral indices (n,m)           :',n,m
  WRITE (*,*) '(3)  Doping (n-p) per length (1/A)  :',doping
  WRITE (*,*) '(4)  Fermi level                    :',Efermi
  WRITE (*,*) '(5)  Refractive index               :',refrac
  WRITE (*,*) '(6)  Background dielectric permittivity:', ebg
  WRITE (*,*) '-----------------------------------------------------'
  WRITE (*,*) '(14) Number of laser photon energies:',nhw_laser
  WRITE (*,*) '(15) Laser photon energy range (eV) :',epmin,epmax
  WRITE (*,*) '(16) Laser linewidth (eV)           :',laser_fwhm
  WRITE (*,*) '(17) Laser polarization angle (deg) :',laser_theta        
  WRITE (*,*) '-----------------------------------------------------'      
  WRITE (*,*) '(18) Electron k points, nk          :',nk
  WRITE (*,*) '(19) Electron DOS energies, nee     :',nee
  WRITE (*,*) '(20) Electron DOS energy range (eV) :',emine,emaxe
  WRITE (*,*) '-----------------------------------------------------' 
  WRITE (*,*) 'Select an item to edit or enter 0 to continue:'
  iedit = 0
  IF (iprofile /= 1) READ (*,*) iedit
  IF (iedit == 0) GOTO 3
  
  ichange = 1
  WRITE(*,*) 'Enter new value(s):'
  
  IF (iedit == 1)  READ (*,*) Tempr
  IF (iedit == 2)  READ (*,*) n,m
  IF (iedit == 3)  READ (*,*) doping
  IF (iedit == 4)  READ (*,*) Efermi
  IF (iedit == 5)  READ (*,*) refrac
  IF (iedit == 6)  READ (*,*) ebg

  IF (iedit == 14) READ (*,*) nhw_laser
  IF (iedit == 15) READ (*,*) epmin,epmax
  IF (iedit == 16) READ (*,*) laser_fwhm
  IF (iedit == 17) READ (*,*) laser_theta  

  IF (iedit == 18) READ (*,*) nk
  IF (iedit == 19) READ (*,*) nee
  IF (iedit == 20) READ (*,*) emine,emaxe        

  GOTO 2
3 CONTINUE
  
!----------------------------------------------------------------------
!                        write parameter file
!----------------------------------------------------------------------
    IF (ichange /= 0) THEN
       IF (imenu == 2) WRITE(*,*) 'Suffix for input data set:'//infile
       WRITE (*,*) 'Enter suffix for output data set (20 char max):'
       READ (*,2002) outfile
    ELSE
       outfile = infile
    END IF

    OPEN(unit=22,file='tube.param.'//outfile)
    WRITE (22,*) Tempr
    WRITE (22,*) n,m
    WRITE (22,*) nk
    WRITE (22,*) doping
    WRITE (22,*) Efermi
    WRITE (22,*) refrac
    WRITE (22,*) ebg
    WRITE (22,*) nhw_laser
    WRITE (22,*) epmin,epmax
    WRITE (22,*) laser_fwhm
    WRITE (22,*) nee,emine,emaxe
    WRITE (22,*) laser_theta

    CLOSE(unit=22)

!----------------------------------------------------------------------
!                           write log file
!----------------------------------------------------------------------
    OPEN(unit=22,file='tube.log.'//outfile)
      
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
  dk = laser_fwhm/(4*hbarvfermi)
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

! compute energy bands En(k) (eV)      
  DO mu=1,nhex 
     DO k=1,nk
        rk=rka(k)
        CALL etbTubeBand(n,m,mu,rk,En,Zk)
        DO ii=1,2
           Enk(ii,mu,k)=En(ii)
           Znk(ii,1:2,mu,k)=Zk(ii,1:2)
        END DO
     END DO
  END DO

  OPEN(unit=22,file='tube.Enk.xyy.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk),((Enk(ii,mu,k),ii=1,2),mu=1,nhex)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'electronic En(k) in tube.Enk.xyy.'//outfile

!Write wave function
  OPEN(unit=22,file='tube.Znk.xyy.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk),((aimag(Znk(ii,1,mu,k)),ii=1,2),mu=1,nhex)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'electronic En(k) in tube.Znk.xyy.'//outfile

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
      
  OPEN(unit=22,file='tube.elecDOS.xyy.'//outfile)
  DO ie = 1, nee
     WRITE(22,1001) eDOS(ie),eEarray(ie)
  END DO
  CLOSE(unit=22)
  WRITE(*,*) 'DOS in tube.elecDOS.xyy.'//outfile

!----------------------------------------------------------------------
!             squared optical matrix elements (eV)
!----------------------------------------------------------------------
  WRITE (*,*) '====================================================='
  WRITE (*,*) '..Squared optical matrix elements (eV)'

  ALLOCATE(cDipole(3,nk,2,nhex,2,nhex))

  DO n1 = 1,2
    DO n2 = 1,2
        DO mu1 = 1, nhex
            DO mu2 = 1, nhex

                DO k = 1, nk
                rk = rka(k)

                    CALL tbDipoleMX(n,m,n1,mu1,n2,mu2,rk,cDipole(1:3,k,n1,mu1,n2,mu2)) ! (1/A)

                END DO

            END DO
         END DO
     END DO
  END DO

          
! x polarization      
  IF(INT(laser_theta) .NE. 0) THEN
     WRITE(*,*) '..x polarization'
     
     ALLOCATE(Px2k(nk,nhex))

! downward cutting line transitions from mu --> mu-1      
     DO mu = 1, nhex 
        DO k = 1, nk
           mmu = mu-1
           IF (mmu.LT.1) mmu = nhex
           Px2k(k,mu) = 3.81D0 * CDABS(cDipole(1,k,1,mu,2,mmu))**2
           ! 3.81 = (hbar^2 / (2 m_e)) (eV-A**2)
        END DO
     END DO

     OPEN(unit=22,file='tube.Px2k_dn.xyy.'//outfile)
     DO k = 1, nk
        WRITE(22,1001) rka(k)/rka(nk),(Px2k(k,mu),mu=1,nhex)
     END DO
     CLOSE(unit=22)
     WRITE(*,*) 'Px2_dn(k) in tube.Px2k_dn.xyy.'//outfile      

! upward cutting line transitions from mu --> mu+1
     DO mu = 1, nhex 
        DO k = 1, nk
           mmu = mu+1
           IF (mmu > nhex) mmu=1
           Px2k(k,mu) = 3.81D0 * CDABS(cDipole(1,k,1,mu,2,mmu))**2
        END DO
     END DO

     OPEN(unit=22,file='tube.Px2k_up.xyy.'//outfile)
     DO k = 1, nk
        WRITE(22,1001) rka(k)/rka(nk),(Px2k(k,mu),mu=1,nhex)
     END DO
     CLOSE(unit=22)
     WRITE(*,*) 'Px2_up(k) in tube.Px2k_up.xyy.'//outfile
      
     DEALLOCATE(Px2k)

  END IF

! z polarization
  WRITE(*,*) '..z polarization'
  ALLOCATE(Pz2k(nk,nhex))

! cutting line transitions from mu --> mu 
  DO mu = 1, nhex 
     DO k = 1, nk
        Pz2k(k,mu) = 3.81D0 * CDABS(cDipole(3,k,1,mu,2,mu))**2
     END DO
  END DO

  OPEN(unit=22,file='tube.Pz2k.xyy.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk),(Pz2k(k,mu),mu=1,nhex)
  END DO
  CLOSE(unit=22)
  WRITE(*,*) 'Pz2(k) in tube.Pz2k.xyy.'//outfile
      
  DEALLOCATE(Pz2k)
            
!----------------------------------------------------------------------
!          real and imaginary part of dielectric function (dimensionless)
!              absorption coefficient (1/cm), conductivity (e^2/h)
!----------------------------------------------------------------------
! allocate storage
  ALLOCATE(hw_laser(nhw_laser))
  ALLOCATE(eps2a(nhw_laser))
  ALLOCATE(eps2a1(nhw_laser))
  ALLOCATE(eps1a(nhw_laser))
  ALLOCATE(eps1a1(nhw_laser))
  ALLOCATE(alpha(nhw_laser))
  ALLOCATE(eelspec(nhw_laser))
  ALLOCATE(sigm1(nhw_laser))
  ALLOCATE(sigm2(nhw_laser))
  ALLOCATE(sigm2_intra(nhw_laser))
  ALLOCATE(sigm2_inter(nhw_laser))
  ALLOCATE(absorpt(nhw_laser))
  
  WRITE (*,*) '====================================================='
  WRITE (*,*) '..Real and Imaginary part of dielectric function'
  WRITE (*,*) '  absorption coefficient (1/cm), conductivity (e^2/h)'
  WRITE (*,*) '..Number of laser energies :',nhw_laser
  WRITE (*,*) '..Laser fwhm linewidth (eV):',laser_fwhm
  
  IF ( doping .ne. 0.) THEN
    Efermi = fermiLevel(n,m,Tempr,doping)
    WRITE (*,*) '..Fermi level calculated from doping', Efermi
  ENDIF

  CALL linArray(nhw_laser,epmin,epmax,hw_laser)
  dep = hw_laser(2)-hw_laser(1)

! unit polarization vectors in complex notation (dimensionless)
  CALL polVector(laser_theta,epol)

! ========================== permittivity ==============================
! real part of dielectric permittivity
! ======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  CALL realDielEn(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,ebg,laser_fwhm,nhw_laser,hw_laser,eps1a) !From Lin's paper
  ! n,m,nhex,nk,rka,Enk,Tempr,Efermi,epol,laser_fwhm

! plot eps1(hw) (Lin's) ************************
  OPEN(unit=22,file='tube.eps1.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eps1a(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'real part of dielectric function in tube.eps1.xyy.'//outfile


  CALL realDielEn_met2(nhw_laser,ebg, hw_laser,eps2a,eps1a1)   !Kramers-Kronig

! plot eps1(hw) (Kramers-Kronig) ****************
  OPEN(unit=22,file='tube.eps1kk.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eps1a1(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'real part of dielectric function in tube.eps1kk.xyy.'//outfile

! imaginary part of dielectric permittivity
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  CALL imagDielEn(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,eps2a)  !From Lin's paper

! plot eps2(hw) ********************************
  OPEN(unit=22,file='tube.eps2.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eps2a(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of dielectric function in tube.eps2.xyy.'//outfile


  CALL imagDielEn_met2(nhw_laser,hw_laser,eps1a,eps2a1)  !Kramers-Kronig

! plot eps2(hw) (Kramers-Kronig) ***************
  OPEN(unit=22,file='tube.eps2kk.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eps2a1(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of dielectric function in tube.eps2kk.xyy.'//outfile
! =======================================================================

! ======================= absorption coefficient ========================
! alpha
  WRITE (*,*) '--------------------------------------------------------'
  CALL imagDielAlpha(nhw_laser,hw_laser,eps2a,refrac,alpha)

! plot absorption alpha(hw) ********************
  OPEN(unit=22,file='tube.alpha.xyy.'//outfile)
  DO ie = 1,nhw_laser
     WRITE(22,1001) hw_laser(ie),alpha(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'alpha(hw) in tube.alpha.xyy.'//outfile
! =======================================================================

! ======================= Electron energy loss spectra  =================
! eels spectra
  WRITE (*,*) '--------------------------------------------------------'
  CALL EELS(nhw_laser,eps1a,eps2a,eelspec)

! plot eels(hw) *******************************
  OPEN(unit=22,file='tube.eels.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie),eelspec(ie)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'eels in tube.eels.xyy.'//outfile
! =======================================================================

! ========================== conductivity ===============================
! real part
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  CALL RealDynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm1)

! plot sigm1(hw) ******************************
  OPEN(unit=22,file='tube.sigm1.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), sigm1(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'real part of conductivity in tube.sigm1.xyy.'//outfile

! imaginary part
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  CALL ImagDynConductivity(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm2)

! plot sigm2(hw) *******************************
  OPEN(unit=22,file='tube.sigm2.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), sigm2(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of conductivity in tube.sigm2.xyy.'//outfile

! absorption
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  CALL Absorption(nhw_laser,eps1a,eps2a,sigm1,sigm2,absorpt)

! plot absorpt(hw) *******************************
  OPEN(unit=22,file='tube.absorpt.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), absorpt(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'absorption in tube.absorpt.xyy.'//outfile

! imaginary part of intraband conductivity
! =======================================================================
  WRITE (*,*) '--------------------------------------------------------'
  CALL ImagDynConductivityIntra(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm2_intra)

! plot sigm2_intra(hw) *******************************
  OPEN(unit=22,file='tube.sigm2_intra.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), sigm2_intra(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of intraband conductivity in tube.sigm2_intra.xyy.'//outfile

! imaginary part of interband conductivity
! =======================================================================
  CALL ImagDynConductivityInter(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,sigm2_inter)

! plot sigm2_inter(hw) *******************************
  OPEN(unit=22,file='tube.sigm2_inter.xyy.'//outfile)
  DO ie = 1, nhw_laser
     WRITE(22,1001) hw_laser(ie), sigm2_inter(ie)/(e2/h)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'imaginary part of interband conductivity in tube.sigm2_inter.xyy.'//outfile


! ============= part of code to check the calculations ================================
! *************** please, remove it later *********************************************

  WRITE (*,*) '====================================================='
  WRITE (*,*) '..test under integral function'

  ALLOCATE(ss0(2,nhex,nhex,nk))
  ALLOCATE(difFermiDisT(2,nhex,nhex,nk))
  ALLOCATE(matrElementSq(2,nhex,nhex,nk))
  ALLOCATE(diracAvgFunc(2,nhex,nhex,nk))

  CALL ImagDynConductivityIntra_test(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,ss0,difFermiDist, &
  matrElementSq,diracAvgFunc)

    max_position = maxloc(ss0)
    min_position = minloc(ss0)

    PRINT*, 'maximum at', max_position
    PRINT*, 'minimum at', min_position

  OPEN(unit=22,file='tube.test_max.xyy.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk), ss0(max_position(1), max_position(2), max_position(3), k), &
     difFermiDist(max_position(1), max_position(2), max_position(3), k), &
     matrElementSq(max_position(1), max_position(2), max_position(3), k), &
     diracAvgFunc(max_position(1), max_position(2), max_position(3), k)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'test max in file tube.test_max.xyy.'//outfile

  OPEN(unit=22,file='tube.test_min.xyy.'//outfile)
  DO k = 1, nk
     WRITE(22,1001) rka(k)/rka(nk), ss0(min_position(1), min_position(2), min_position(3), k), &
     difFermiDist(min_position(1), min_position(2), min_position(3), k), &
     matrElementSq(min_position(1), min_position(2), min_position(3), k), &
     diracAvgFunc(min_position(1), min_position(2), min_position(3), k)
  ENDDO
  CLOSE(unit=22)
  WRITE(*,*) 'test min in file tube.test_min.xyy.'//outfile

  DEALLOCATE(ss0,difFermiDist,matrElementSq,diracAvgFunc)

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
!  kCoef = (/10., 8., 6., 4., 2., 1., 0.8, 0.5, 0.2, 0.1/)
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
!  CALL imagDielEn(n,m,nhex,nk,rka,Enk,cDipole,Tempr,Efermi,epol,laser_fwhm,nhw_laser,hw_laser,eps2a)
!
!  eps2aii(i,1:nhw_laser) = eps2a(1:nhw_laser)
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
!        divergence(i) = divergence(i) + (eps2aii(i+1,ie) - eps2aii(i,ie))**2
!    END DO
!    divergence(i) = SQRT(divergence(i))
!  END DO
!
!  OPEN(unit=22,file='tube.divergence.xyy.'//outfile)
!  WRITE(22,1001) kCoef(1), 0.D0, 1.D0/kCoef(1)
!  DO i=1,9
!     WRITE(22,1001) kCoef(i+1), divergence(i), 1.D0/kCoef(i+1)
!  ENDDO
!  CLOSE(unit=22)
!  WRITE(*,*) 'divergence in file tube.divergence.xyy.'//outfile
!
!  DEALLOCATE(eps2aii)

! ============= part of code to check the calculations ================================
! ************************************* END *******************************************

 WRITE(*,*) '====================================================='
 WRITE(*,*) 'Fermi level: ', 'from input', Efermi, 'from doping', fermiLevel(n,m,Tempr,doping)


!----------------------------------------------------------------------
!                          format statements
!----------------------------------------------------------------------
      WRITE(*,*) '====================================================='
      STOP

1001  FORMAT(901(1X,G13.5))
1002  FORMAT(2F8.4)
2002  FORMAT(1A20)

END PROGRAM cntabsorpt
!*******************************************************************************
