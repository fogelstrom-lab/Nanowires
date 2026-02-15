module integration_scalar
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: integrate_complex_adaptive

  abstract interface
    subroutine f_complex_sub(x, fx)
      import dp
      real(dp),    intent(in)  :: x
      complex(dp), intent(out) :: fx
    end subroutine f_complex_sub
  end interface

contains

  subroutine integrate_complex_adaptive(f, a, b, I, atol, rtol, max_depth)
    procedure(f_complex_sub)         :: f
    real(dp),    intent(in)          :: a, b
    complex(dp), intent(out)         :: I
    real(dp),    intent(in)          :: atol, rtol
    integer,     intent(in), optional :: max_depth

    integer     :: depth_max
    real(dp)    :: m, h
    complex(dp) :: fa, fm, fb, S0

    if (abs(b - a) < tiny(1.0_dp)) then
      I = (0.0_dp, 0.0_dp)
      return
    end if

    depth_max = 20
    if (present(max_depth)) depth_max = max_depth

    m = 0.5_dp*(a + b)
    h = b - a

    call f(a, fa)
    call f(m, fm)
    call f(b, fb)

    S0 = (h / 6.0_dp) * (fa + 4.0_dp*fm + fb)

    call asimpson(a, b, fa, fm, fb, S0, atol, rtol, 0, depth_max, I)

  contains

    recursive subroutine asimpson(a, b, fa, fm, fb, S, atol, rtol, depth, depth_max, result)
      real(dp),    intent(in) :: a, b, atol, rtol
      integer,     intent(in) :: depth, depth_max
      complex(dp), intent(in) :: fa, fm, fb, S
      complex(dp), intent(out):: result

      real(dp)    :: m, m1, m2, h, err, scale, tol_eff
      complex(dp) :: fl, fr, S_left, S_right, S_lr, resL, resR

      h  = b - a
      m  = 0.5_dp*(a + b)
      m1 = 0.5_dp*(a + m)
      m2 = 0.5_dp*(m + b)

      call f(m1, fl)
      call f(m2, fr)

      S_left  = (h / 12.0_dp) * (fa + 4.0_dp*fl + fm)
      S_right = (h / 12.0_dp) * (fm + 4.0_dp*fr + fb)
      S_lr    = S_left + S_right

      err     = abs(S_lr - S)
      scale   = max(abs(S_lr), abs(S))
      tol_eff = max(atol, rtol * scale)

      if (err <= 15.0_dp*tol_eff .or. depth >= depth_max) then
        result = S_lr
      else
        call asimpson(a, m, fa, fl, fm, S_left,  0.5_dp*atol, rtol, depth+1, depth_max, resL)
        call asimpson(m, b, fm, fr, fb, S_right, 0.5_dp*atol, rtol, depth+1, depth_max, resR)
        result = resL + resR
      end if
    end subroutine asimpson

  end subroutine integrate_complex_adaptive

end module integration_scalar
