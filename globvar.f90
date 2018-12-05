!*******************************************************************************
MODULE globvar
! dummy module
!===============================================================================
! Global variables
!-------------------------------------------------------------------------------
! variables passed to external program
!  INTEGER                :: nhex
!  INTEGER                :: nvar
!  INTEGER                :: nk
!  REAL(8)                :: refrac
!
!  INTEGER, ALLOCATABLE   :: nmk2index(:,:,:) !(2,nhex,nk)
!  REAL(8), ALLOCATABLE   :: Rnk(:,:,:,:,:)   !(2,nhex,2,nhex,nk)
!!-------------------------------------------------------------------------------
!! variables passed to function sm(t)
!  INTEGER                :: nt
!  REAL(8), ALLOCATABLE   :: times(:)         !(nt)
!  REAL(8), ALLOCATABLE   :: Smt              !(4,nt)
!!-------------------------------------------------------------------------------
!! variables passed to fk(t)
!  REAL(8), ALLOCATABLE   :: fkt(:,:,:,:)     !(2,nhex,nk,nt)
!
END MODULE globvar
!*******************************************************************************
