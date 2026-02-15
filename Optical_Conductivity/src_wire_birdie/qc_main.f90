!-------------------------------------------------------------------
program qc_main
  use mpi_f08
  use qc_data
  use qc_mod
  use qc_cond
  use iso_fortran_env, only: dp => real64

  implicit none
  integer :: itype

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

! write(*,*) ' Make initialisations',myrank
  call initialize_constants()
! write(*,*) ' Getting input',myrank
  call get_input
! write(*,*) ' Got input',myrank
!
! > Below for a simple calculation based on the input in the file input.dat
!  
  itype = 0
  call get_equilibrium(itype)
  call qc_conductance
!
! > Below for maping the spectral gap as a function of potential scattering rate
!  
! do while (tau(1) <= 0.05*delta0)
!    itype = 1
!    call get_equilibrium(itype)
!    tau(1) = tau(1) + 0.0025*delta0
!    call MPI_Bcast(tau   , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! enddo    
! 
! do while (tau(1) <= 2.0*delta0) 
!    itype = 1
!    call get_equilibrium(itype)
!    tau(1) = tau(1) + 0.025*delta0
!    call MPI_Bcast(tau   , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! enddo    
!
! > Below for maping the spectral gap as a function of magnetic scattering rate
!  
! do while (tau(2) <= 0.05*delta0)
!    itype = 1
!    call get_equilibrium(itype)
!    tau(2) = tau(2) + 0.0025*delta0
!    call MPI_Bcast(tau   , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! enddo    
! 
! do while (tau(2) <= 2.0*delta0) 
!    itype = 1
!    call get_equilibrium(itype)
!    tau(2) = tau(2) + 0.025*delta0
!    call MPI_Bcast(tau   , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! enddo    
!
! > Below for maping Andreev-level spectra vs junction length
!
!  Ngap=0.0
!  call MPI_Bcast(Ngap   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!  do while (Ngap <= 0.80*Rx) 
!     itype = 3
!     call get_equilibrium(itype)
!     Ngap = Ngap + 0.1
!     call MPI_Bcast(Ngap   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!  enddo    
!
! > Below for maping Andreev-level spectra vs phase difference over the junction
!
! phase=0.0
! call MPI_Bcast(phase   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! do while (phase <= 2.002*pi) 
!    itype = 4
!    call get_equilibrium
!    phase = phase + 0.010*pi
!    call MPI_Bcast(phase   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! enddo    
 
  call MPI_Finalize(ierr)
end program qc_main
