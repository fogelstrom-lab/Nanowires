module math_utils
  use iso_fortran_env, only: dp => real64
  implicit none
contains
!
! --- Safe complex division: res = d / n  ----------------------
!
  function cdiv(d, n) result(res)
    
    complex(dp), intent(in) :: d, n
    complex(dp)             :: res
    complex(dp)             :: cn
    real(dp)                :: denom

    cn    = conjg(n)
    denom = real(n * cn, kind=dp)
    if (denom < tiny(1.0_dp)) error stop "cdiv: division by zero"
    res = (d * cn) / denom
  end function cdiv
!
! --- Safe complex 2x2-matrix division (note non-commutative: res = d * n^-1  ----------------------
!
  pure function mdivR(d, n) result(res)
    
    complex(dp), intent(in) :: d(2,2), n(2,2)
    complex(dp)             :: res(2,2)
    complex(dp)             :: cn(2,2)      ! Inverse of n

    cn  = inv2x2(n)
    res = matmul(d,cn)
  end function mdivR
!
! --- Safe complex 2x2-matrix division (note non-commutative: res = n^-1 * d  ----------------------
!
  pure function mdivL(d, n) result(res)
    
    complex(dp), intent(in) :: d(2,2), n(2,2)
    complex(dp)             :: res(2,2)
    complex(dp)             :: cn(2,2)      ! Inverse of n

    cn  = inv2x2(n)
    res = matmul(cn,d)
  end function mdivL
!
! --- Inverse of a 2x2 real or complex matrix  ----------------------
!
  pure function inv2x2(A) result(Ainv)
    
    implicit none
    complex(dp), intent(in) :: A(2,2)
    complex(dp) :: Ainv(2,2)
    complex(dp) :: det

    det = det2x2(A)

    Ainv(1,1) =  A(2,2) / det
    Ainv(1,2) = -A(1,2) / det
    Ainv(2,1) = -A(2,1) / det
    Ainv(2,2) =  A(1,1) / det

  end function inv2x2
!
! --- Determinant of a 2x2 real or complex matrix  ----------------------
!
  pure function det2x2(A) result(det)
    
    implicit none
    complex(dp), intent(in) :: A(2,2)
    complex(dp) :: det

    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

  end function det2x2
!
! --- Computes the principal square root of a general 2x2 complex matrix A. --
! --- DReference: Higham (2008), Functions of Matrices, §6.2                --
!
  function sqrtm2x2(M) result(SQ)

    implicit none
    complex(dp), intent(in) :: M(2,2)
    complex(dp) :: SQ(2,2)
    complex(dp) :: a,b,c,d, t, delta, s, r

    a = M(1,1);  b = M(1,2)
    c = M(2,1);  d = M(2,2)

    t      = a + d
    delta  = a*d - b*c
    s      = sqrt(delta)
    r      = sqrt(t + 2.0_dp*s)

    if (abs(r) == 0.0_dp) then
       error stop "sqrtm2x2: singular or nilpotent matrix (r=0)"
    end if

    SQ(1,1) = (a + s) / r
    SQ(1,2) =  b      / r
    SQ(2,1) =  c      / r
    SQ(2,2) = (d + s) / r
  end function sqrtm2x2
!
!------------------------------------------------------------------------
!
end module math_utils
!
!------------------------------------------------------------------------
