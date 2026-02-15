module mpi_utils
  use mpi_f08
  use math_utils, only: dp
  implicit none
  private
  public :: sync_cols_gather_indexed_rpairs

contains

  subroutine sync_cols_gather_indexed_rpairs(a, comm)
    !! Assemble complex(dp) a(:,:) distributed over the SECOND index (columns).
    !! Each rank computes whole columns (fixed j, all i).
    !! Safe for any nprocs/bounds; safe to call multiple times.
    complex(dp), intent(inout) :: a(:,:)          ! a(iL:iU, jL:jU)
    type(MPI_Comm), intent(in) :: comm

    integer :: myrank, nprocs
    integer :: iL, iU, jL, jU, n1, n2
    integer :: i, j, c, mycols, total_cols
    integer, allocatable :: j_local(:), j_all(:)
    integer, allocatable :: cnt_c(:), disp_c(:)
    integer, allocatable :: cnt_r(:), disp_r(:)
    real(dp), allocatable :: send_buf(:), recv_buf(:)
    integer :: rr

    call MPI_Comm_rank(comm, myrank)
    call MPI_Comm_size(comm, nprocs)

    iL = lbound(a,1); iU = ubound(a,1); n1 = iU - iL + 1
    jL = lbound(a,2); jU = ubound(a,2); n2 = jU - jL + 1

    ! 1) Build local list of column indices (energy indices) this rank owns
    mycols = 0
    do j = jL, jU
      if (modulo(j - jL, nprocs) == myrank) mycols = mycols + 1
    end do

    allocate(j_local(max(1,mycols)))
    if (mycols > 0) then
      mycols = 0
      do j = jL, jU
        if (modulo(j - jL, nprocs) == myrank) then
          mycols = mycols + 1
          j_local(mycols) = j
        end if
      end do
    end if

    ! 2) Allgather the column index lists
    allocate(cnt_c(nprocs), disp_c(nprocs))
    call MPI_Allgather(mycols, 1, MPI_INTEGER, cnt_c, 1, MPI_INTEGER, comm)

    disp_c(1) = 0
    do c = 2, nprocs
      disp_c(c) = disp_c(c-1) + cnt_c(c-1)
    end do
    total_cols = sum(cnt_c)

    allocate(j_all(max(1,total_cols)))
    call MPI_Allgatherv(j_local, max(0,mycols), MPI_INTEGER, &
                        j_all, cnt_c, disp_c, MPI_INTEGER, comm)

    ! 3) Pack each owned COLUMN (length n1) into REAL(dp) [Re,Im]
    allocate(send_buf(max(1, 2*n1*mycols)))
    rr = 0
    if (mycols > 0) then
      do c = 1, mycols
        j = j_local(c)
        do i = iL, iU
          rr = rr + 1
          send_buf(2*rr-1) = real(a(i,j), dp)
          send_buf(2*rr  ) = aimag(a(i,j))
        end do
      end do
    end if

    ! 4) Exchange all columns
    allocate(cnt_r(nprocs), disp_r(nprocs))
    do c = 1, nprocs
      cnt_r(c) = 2 * n1 * cnt_c(c)
    end do
    disp_r(1) = 0
    do c = 2, nprocs
      disp_r(c) = disp_r(c-1) + cnt_r(c-1)
    end do

    allocate(recv_buf(max(1, 2*n1*total_cols)))
    call MPI_Allgatherv(send_buf, 2*n1*max(0,mycols), MPI_DOUBLE_PRECISION, &
                        recv_buf, cnt_r, disp_r, MPI_DOUBLE_PRECISION, comm)

    ! 5) Unpack back into a(:, j)
    if (total_cols > 0) then
      rr = 0
      do c = 1, total_cols
        j = j_all(c)
        do i = iL, iU
          rr = rr + 1
          a(i, j) = cmplx(recv_buf(2*rr-1), recv_buf(2*rr), kind=dp)
        end do
      end do
    end if

    deallocate(j_local, j_all, cnt_c, disp_c, cnt_r, disp_r, send_buf, recv_buf)
  end subroutine sync_cols_gather_indexed_rpairs

end module mpi_utils
