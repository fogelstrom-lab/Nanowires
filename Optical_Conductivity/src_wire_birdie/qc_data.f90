!-------------------------------------------------------------------
! Modernized module replacing COMMON blocks from qc.dat
!-------------------------------------------------------------------
module qc_data
  use iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: nop = 3*2 ! Number of (real-valued) selfenergies
  complex(dp) :: czero, cone, w, en, pauli0(2,2), pauli1(2,2), pauli2(2,2), pauli3(2,2)
  real(dp) :: pi, Rx, temp, delta0, etaR, phase, error_tol, Ngap, ndeltas, Ec, emax, de, ex_old
  integer :: sx, iemax

  complex(dp), allocatable :: gr_1(:,:), gr_2(:,:), delx(:)
  complex(dp), allocatable :: gimp0(:), fimp0(:), timp0(:), uimp0(:)
  complex(dp), allocatable :: gimpP(:,:), fimpP(:,:), timpP(:,:), uimpP(:,:) 
  complex(dp), allocatable :: gimpM(:,:), fimpM(:,:), timpM(:,:)

  real(dp),    allocatable :: grid(:), er(:)
  real(dp)                 :: tau(2),sigma(2), srate(0:2)
  complex(dp) :: delL, delR

  complex(dp), allocatable :: avj(:,:), avg(:,:), avf(:,:), avt(:,:)
  complex(dp), allocatable :: gr_1p(:,:), gr_2p(:,:), gr_1m(:,:), gr_2m(:,:)
  complex(dp), allocatable ::  gimp(:,:),  fimp(:,:),  timp(:,:), uimp(:,:)

  complex(dp), allocatable :: alpha1p(:,:), beta1p(:,:)
  complex(dp), allocatable :: alpha1m(:,:), beta1m(:,:)
  complex(dp), allocatable :: alpha2p(:,:), beta2p(:,:)
  complex(dp), allocatable :: alpha2m(:,:), beta2m(:,:)


! MPI variables

  integer :: myrank, nprocs, ierr

contains

subroutine initialize_constants()
    ! Initialize parameters and (re)allocate arrays based on sx and epts.
    implicit none

    czero = (0.0_dp, 0.0_dp)
    cone  = (1.0_dp, 0.0_dp)
    w     = (0.0_dp, 1.0_dp)
    en    = (0.0_dp, 0.0_dp)
    pi    = acos(-1.0_dp)
    Rx    = 0.0_dp
    temp  = 0.0_dp
    delta0 = 0.0_dp
    phase  = 0.0_dp
    tau    = 0.0_dp
    sigma  = 0.0_dp
    error_tol    = 0.0_dp
    emax   = 0
    ndeltas = 0

    pauli0=reshape([(1.0_dp,0.0_dp),czero,czero,(1.0_dp,0.0_dp)],[2,2])
    pauli1=reshape([czero,(1.0_dp,0.0_dp),(1.0_dp,0.0_dp),czero],[2,2])
    pauli2=reshape([czero,w,-w,czero],[2,2])
    pauli3=reshape([(1.0_dp,0.0_dp),czero,czero,(-1.0_dp,0.0_dp)],[2,2])

  end subroutine initialize_constants
end module qc_data
