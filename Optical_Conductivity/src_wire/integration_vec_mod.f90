module integration_vec
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: integrate_vec_adaptive

  abstract interface
    subroutine f_vec_interface(x, f)
      import dp
      real(dp),    intent(in)  :: x
      complex(dp), intent(out) :: f(:)
    end subroutine f_vec_interface
  end interface

contains

  subroutine integrate_vec_adaptive(f_vec, a, b, f_int, atol, rtol, max_depth)
    procedure(f_vec_interface)       :: f_vec
    real(dp),          intent(in)    :: a, b
    complex(dp),       intent(out)   :: f_int(:)
    real(dp),          intent(in)    :: atol, rtol
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

    call f_vec(a, fa)
    call f_vec(b, fb)
    m = 0.5_dp*(a + b)
    call f_vec(m, fm)

    h  = b - a
    S0 = (h / 6.0_dp) * (fa + 4.0_dp*fm + fb)

    call asimpson_vec(a, b, fa, fm, fb, S0, atol, rtol, 0, depth_max, f_int, ncomp)

  contains

    recursive subroutine asimpson_vec(a_loc, b_loc, fa_loc, fm_loc, fb_loc, S, atol_loc, &
                                      rtol_loc, depth, depth_max_loc, result, ncomp_loc)
      real(dp),    intent(in) :: a_loc, b_loc, atol_loc, rtol_loc
      integer,     intent(in) :: depth, depth_max_loc, ncomp_loc
      complex(dp), intent(in) :: fa_loc(ncomp_loc), fm_loc(ncomp_loc), fb_loc(ncomp_loc)
      complex(dp), intent(in) :: S(ncomp_loc)
      complex(dp), intent(out):: result(ncomp_loc)

      real(dp)    :: h_loc, m_loc, m1, m2
      complex(dp) :: fl(ncomp_loc), fr(ncomp_loc)
      complex(dp) :: S_left(ncomp_loc), S_right(ncomp_loc), S_lr(ncomp_loc)
      complex(dp) :: resL(ncomp_loc), resR(ncomp_loc)
      real(dp)    :: err, scale, tol_loc_eff

      h_loc  =         b_loc - a_loc
      m_loc  = 0.5_dp*(a_loc + b_loc)
      m1     = 0.5_dp*(a_loc + m_loc)
      m2     = 0.5_dp*(m_loc + b_loc)

      call f_vec(m1, fl)
      call f_vec(m2, fr)

      S_left  = (h / 12.0_dp) * (fa_loc + 4.0_dp*fl + fm_loc)
      S_right = (h / 12.0_dp) * (fm_loc + 4.0_dp*fr + fb_loc)
      S_lr    = S_left + S_right

      err   = maxval( abs(S_lr - S) )
      scale = max( maxval(abs(S_lr)), maxval(abs(S)) )
      tol_loc_eff = max(atol_loc, rtol_loc * scale)

      if (err <= 15.0_dp*tol_loc_eff .or. depth >= depth_max_loc) then
         result = S_lr
      else
         ! Dela tol_loceransen mellan delintervallen (konservativt)
         call asimpson_vec(a_loc, m, fa_loc, fl, fm_loc, S_left,  0.5_dp*atol_loc, rtol_loc, depth+1, depth_max, resL, ncomp_loc)
         call asimpson_vec(m, b_loc, fm_loc, fr, fb_loc, S_right, 0.5_dp*atol_loc, rtol_loc, depth+1, depth_max, resR, ncomp_loc)
         result = resL + resR
      end if
    end subroutine asimpson_vec

  end subroutine integrate_vec_adaptive

end module integration_vec
