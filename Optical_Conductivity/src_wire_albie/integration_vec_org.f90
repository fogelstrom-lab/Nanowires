module integration_vec
  use iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: integrate_vec_adaptive

  ! Interface: vector-valued complex integrand
  abstract interface
    subroutine f_vec_interface(x, f)
      import dp
      real(dp),    intent(in)  :: x
      complex(dp), intent(out) :: f(:)
    end subroutine f_vec_interface
  end interface

contains

  !--------------------------------------------------------------------
  ! Adaptive Simpson integration of a complex vector-valued function.
  !
  ! f_vec:   user-supplied routine: given x, returns f(x)(:)
  ! a, b:    integration limits
  ! f_int:   result vector; size determines number of components
  ! tol:     requested absolute tolerance (on max-norm of error)
  ! max_depth (optional): maximum recursion depth (default: 20)
  !--------------------------------------------------------------------
  subroutine integrate_vec_adaptive(f_vec, a, b, f_int, tol, max_depth)
    procedure(f_vec_interface)       :: f_vec
    real(dp),          intent(in)    :: a, b
    complex(dp),       intent(out)   :: f_int(:)   ! caller allocates, sets size
    real(dp),          intent(in)    :: tol
    integer, optional, intent(in)    :: max_depth

    integer      :: depth_max, ncomp
    real(dp)     :: h, m
    complex(dp)  :: fa(size(f_int)), fm(size(f_int)), fb(size(f_int)), S0(size(f_int))

    if (abs(b - a) < tiny(1.0_dp)) then
       f_int = (0.0_dp, 0.0_dp)
       return
    end if

    depth_max = 20
    if (present(max_depth)) depth_max = max_depth

    ncomp = size(f_int)

    ! Initial function evaluations
    call f_vec(a, fa)
    call f_vec(b, fb)
    m = 0.5_dp*(a + b)
    call f_vec(m, fm)

    h  = b - a
    S0 = (h / 6.0_dp) * (fa + 4.0_dp*fm + fb)

    call asimpson_vec(a, b, fa, fm, fb, S0, tol, 0, depth_max, f_int, ncomp)

  contains

    !--------------------------------------------------------------
    ! Internal recursive adaptive Simpson.
    ! Note: has host access to f_vec, so we don't pass it as arg.
    !--------------------------------------------------------------
    recursive subroutine asimpson_vec(a, b, fa, fm, fb, S, tol, depth, depth_max, result, ncomp)
      real(dp),    intent(in) :: a, b, tol
      integer,     intent(in) :: depth, depth_max, ncomp
      complex(dp), intent(in) :: fa(ncomp), fm(ncomp), fb(ncomp)
      complex(dp), intent(in) :: S(ncomp)
      complex(dp), intent(out):: result(ncomp)

      real(dp)    :: h, m, m1, m2
      complex(dp) :: fl(ncomp), fr(ncomp)
      complex(dp) :: S_left(ncomp), S_right(ncomp), S_lr(ncomp)
      complex(dp) :: resL(ncomp), resR(ncomp)
      real(dp)    :: err

      h  = b - a
      m  = 0.5_dp*(a + b)
      m1 = 0.5_dp*(a + m)
      m2 = 0.5_dp*(m + b)

      call f_vec(m1, fl)
      call f_vec(m2, fr)

      S_left  = (h / 12.0_dp) * (fa + 4.0_dp*fl + fm)
      S_right = (h / 12.0_dp) * (fm + 4.0_dp*fr + fb)
      S_lr    = S_left + S_right

      err = maxval( abs(S_lr - S) )

      if (err <= 15.0_dp*tol .or. depth >= depth_max) then
         result = S_lr
      else
         call asimpson_vec(a, m, fa, fl, fm, S_left,  0.5_dp*tol, depth+1, depth_max, resL, ncomp)
         call asimpson_vec(m, b, fm, fr, fb, S_right, 0.5_dp*tol, depth+1, depth_max, resR, ncomp)
         result = resL + resR
      end if
    end subroutine asimpson_vec

  end subroutine integrate_vec_adaptive

end module integration_vec
