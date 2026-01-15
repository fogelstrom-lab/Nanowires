Module qc_cond
  use mpi_f08
  use qc_data
  use qc_farm
  use iso_fortran_env, only: dp => real64
  use mpi_utils, only: sync_cols_gather_indexed_rpairs
  use math_utils, only: cdiv   ! <--- import the function
  use integration_scalar,  only: integrate_complex_adaptive

  implicit none
  integer :: nfreq, counter, icall_lin
  real(dp) :: omega, omh
  complex(dp) :: cond
  real(dp), allocatable :: frequencies(:)

  public :: qc_conductance, integral, linear_resp
contains
!-----------------------------------------------------------------------
!
   subroutine qc_conductance
!
!-----------------------------------------------------------------------
!
   integer :: i, ir
   complex(dp), allocatable :: condv(:), lcondv(:)
   real(dp) :: sumsrate, oga, d2
   integer :: ulog, ios, pstep
   character(len=128) :: fname

      sumsrate = tau(1)+2.0*tau(2)
      call setfrec(sumsrate)
      if(.not.allocated( condv)) allocate( condv(1:nfreq))
      if(.not.allocated(lcondv)) allocate(lcondv(1:nfreq))

      condv  =  (0.0_dp, 0.0_dp)
      lcondv =  (0.0_dp, 0.0_dp)
      d2 = delta0*delta0

      ulog = 100 + myrank
      write(fname,'("qc_cond_log_rank",i0,".dat")') myrank

       open(ulog, file=fname, status="replace", action="write", iostat=ios)
       if (ios /= 0) then
           write(*,*) "rank", myrank, "failed to open", trim(fname), " iostat=", ios
           call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       end if

! later inside loop:
      if(myrank == 0) write(*,*) ' going with ',nfreq,' frequencies'
      do i=1+myrank, nfreq, nprocs
         ex_old = 3000.0
         counter = 0
         icall_lin = 0
         omega = frequencies(i)
         omh = omega/2.0
         oga = omega/delta0
         call integral
         lcondv(i) = real(cond)/sumsrate + w*omega*aimag(cond)/2.0
         if(myrank == 0) write(*,1100) myrank, oga, lcondv(i), 1.0/(omega**2+d2*sumsrate**2)
         write(ulog,*) oga, counter, icall_lin ! logg how many self-energy iterations per frecquency
      enddo
      close(ulog)

      call MPI_Reduce(lcondv, condv, nfreq, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
          write(*,*) "rank", myrank, "MPI_Reduce failed ierr=", ierr
          call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
      end if
      
      if(myrank == 0) then

         open(10,file="avcond.dat", status="replace", action="write", iostat=ios)
         if (ios/=0) call MPI_Abort(MPI_COMM_WORLD, 10, ierr)

         do i=1,nfreq
            omega = frequencies(i)
            oga = omega/delta0
            write(10,1000) oga, condv(i) ,1.0/(omega**2+sumsrate**2)
         enddo
         close(10)
      endif
!
 1000 format(2000(1x,e14.8))
 1100 format(1x,i2,200(1x,e14.8))
!
!--------------------------------------------------------------------------
!
   end subroutine qc_conductance
!
!--------------------------------------------------------------------------
!
   subroutine integral
!
!--------------------------------------------------------------------------
!
      integer  :: iestep, ie, maxdepth
      real(dp) :: epsl, epsu, dom, atol, rtol
      complex(dp) :: icond

      atol = 1.0e-7_dp
      rtol = 1.0e-5_dp
      maxdepth = 30
      icond = (0.0_dp, 0.0_dp)
      cond  = (0.0_dp, 0.0_dp)
      dom = max(omega,5.0*etaR)

      iestep = iemax/100
!     if(myrank == 0) write(*,*) ' integrating :', omega

      do ie = -iemax, iemax-iestep, iestep

         epsl = er(ie)
         epsu = er(ie+iestep)
         call integrate_complex_adaptive(linear_resp, epsl, epsu, icond, atol, rtol, maxdepth)
         cond = cond + icond
      end do
      cond = cond/dom

!1000 format(30(1x,e14.8))
!
!--------------------------------------------------------------------------
!
   end subroutine integral
!
!--------------------------------------------------------------------------
!
   subroutine linear_resp(energy, icond)
!
!-----------------------------------------------------------------------
!
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp),    intent(in)  :: energy
      complex(dp), intent(out) :: icond

      integer :: i, iep, iem, icont, icmp, iprint, imploop, maxloop
      real(dp) :: d, thp, thm, dth, err, epsp, epsm
      complex(dp) :: gold, told, xold, aold, gexp, wd, kshift, vp, vt

      complex(dp)  NRp, NRm, NKp, NKm
      complex(dp)  gRp, tRp, Xxp, Xap
      complex(dp)  gRm, tRm, Xxm, Xam
      complex(dp)  R1p, R2p, R1m, R2m
      complex(dp)  X1p, X2p, X1m, X2m

      complex(dp)  dgRp, dgRm, dgXp, dgXm
      complex(dp)  dhRp, dhRm, dhXp, dhXm
      complex(dp)  dfRp, dfRm, dfXp, dfXm
      complex(dp)  dtRp, dtRm, dtXp, dtXm
      complex(dp)   agR, ahR, afR, atR
      complex(dp)  agX, ahX, afX, atX

      complex(dp)  lcurrR, lcurrK
!
! -- Get the right leading-order self energies
!
      iprint = 1
      icmp = 0
      icall_lin = icall_lin+1
!
! --  en = er+omega/2
!
      iep = 1
      epsp = energy+omh
      thp=tanh(0.5*epsp/Temp)
      call get_imp_at_e(epsp, iep, icmp, iprint)
      counter = counter + icmp
!
! -- en = er-omega/2
!
      iem = 2
      epsm = energy-omh
      thm=tanh(0.5*epsm/Temp)
      call get_imp_at_e(epsm, iem, icmp, iprint)
      counter = counter + icmp
!
! -- Generate the linear response Green's function
!
      kshift = 2.0*w*etaR
      dth=(thp-thm)
      err =100.0
      maxloop = 1 
      imploop = 0  
!
!  Note that advanced coherence functions are computed by symmetry: a =-conjg(\tilde r)
!  and similarly for \tilde a. 
!
!-- p.v > 0 --
!
      NRp=-w/((1.0+gr_1p(iep)*gr_2p(iep))*(1.0+gr_1p(iem)*gr_2p(iem)))
      NKp=-w/((1.0+gr_1p(iep)*gr_2p(iep))*conjg(1.0+gr_1p(iem)*gr_2p(iem)))
 
      R1p = alpha1p(iep)+beta1p(iem)
      R2p = alpha2p(iep)+beta2p(iem)

      X1p = alpha1p(iep)-conjg(alpha1p(iem))-kshift
      X2p = alpha2p(iep)-conjg(alpha2p(iem))-kshift
!
!-- p.v < 0 --
!
      NRm=-w/((1.0+gr_1m(iep)*gr_2m(iep))*(1.0+gr_1m(iem)*gr_2m(iem)))
      NKm=-w/((1.0+gr_1m(iep)*gr_2m(iep))*conjg(1.0+gr_1m(iem)*gr_2m(iem)))      

      R1m = alpha1m(iep)+beta1m(iem)
      R2m = alpha2m(iep)+beta2m(iem)

      X1m = alpha1m(iep)-conjg(alpha1m(iem))-kshift
      X2m = alpha2m(iep)-conjg(alpha2m(iem))-kshift
!
!-- p.v > 0 --
!
         vp = cone
         vt =-cone
         gRp = (vp*gr_1p(iem)-gr_1p(iep)*vt)/R1p
         tRp =-(vt*gr_2p(iem)-gr_2p(iep)*vp)/R2p

         vp = cone*dth
         vt =-cone*dth
         Xxp =-(vp-gr_1p(iep)*vt*conjg(gr_1p(iem)))/X1p
         Xap =-(vt-gr_2p(iep)*vp*conjg(gr_2p(iem)))/X2p

!
!-- p.v < 0 --  The vertex corrections ~ p and thus the sign change in vp and vt
!
         vp =-cone
         vt = cone
         gRm = (vp*gr_1m(iem)-gr_1m(iep)*vt)/R1m
         tRm =-(vt*gr_2m(iem)-gr_2m(iep)*vp)/R2m

         vp =-cone*dth
         vt = cone*dth
         Xxm =-(vp-gr_1m(iep)*vt*conjg(gr_1m(iem)))/X1m
         Xam =-(vt-gr_2m(iep)*vp*conjg(gr_2m(iem)))/X2m

!
!-- Make the Green's function p.v > 0, a factor two less due to averaging over direction
!
         dgRp =-NRp*(gRp*gr_2p(iem)+gr_1p(iep)*tRp)
         dfRp = NRp*(gRp-gr_1p(iep)*tRp*gr_1p(iem))
         dtRp = NRp*(tRp-gr_2p(iep)*gRp*gr_2p(iem))
         dhRp = NRp*(tRp*gr_1p(iem)+gr_2p(iep)*gRp)

         dgXp = NKp*(Xxp-gr_1p(iep)*Xap*conjg(gr_1p(iem)))
         dfXp =-NKp*(gr_1p(iep)*Xap-Xxp*conjg(gr_1p(iem)))
         dtXp = NKp*(gr_2p(iep)*Xxp-Xap*conjg(gr_2p(iem)))
         dhXp = NKp*(Xap-gr_2p(iep)*Xxp*conjg(gr_2p(iem)))
!
!-- Make the Green's funct:on p.v < 0
!
         dgRm =-NRm*(gRm*gr_2m(iem)+gr_1m(iep)*tRm)
         dfRm = NRm*(gRm-gr_1m(iep)*tRm*gr_1m(iem))
         dtRm = NRm*(tRm-gr_2m(iep)*gRm*gr_2m(iem))
         dhRm = NRm*(tRm*gr_1m(iem)+gr_2m(iep)*gRm)

         dgXm = NKm*(Xxm-gr_1m(iep)*Xam*conjg(gr_1m(iem)))
         dfXm =-NKm*(gr_1m(iep)*Xam-Xxm*conjg(gr_1m(iem)))
         dtXm = NKm*(gr_2m(iep)*Xxm-Xam*conjg(gr_2m(iem)))
         dhXm = NKm*(Xam-gr_2m(iep)*Xxm*conjg(gr_2m(iem)))
!
!-- Current contribution
!
      lcurrR =  (dgRp-dhRp-(dgRm-dhRm))*thm-thp*conjg(dgRp-dhRp-(dgRm-dhRm))
      lcurrK =  (dgXp-dhXp-(dgXm-dhXm))

      icond = 0.25*(lcurrR+lcurrK)

!1000 format(30(1x,e14.8))
!
!-----------------------------------------------------------------------
!
   end subroutine linear_resp
!
!-----------------------------------------------------------------------
!
   subroutine setfrec(sumsrate)
!
!--------------------------------------------------------------------------
!
      real(dp), intent(in) :: sumsrate 
      real(dp) :: dom0, dom, oga, intlim
      integer :: ifreq
   
      intlim = 30.0_dp * max(delta0,delta0*sumsrate)
      dom0 = 0.02_dp * delta0

      oga = dom0/2.0_dp
      ifreq = 0

      do while (oga < intlim)
         ifreq = ifreq+1

         if(oga < 5.0_dp*dom0) then
            dom = dom0/10.0_dp
         elseif(oga < 2.0_dp*delta0 - 10.0_dp*dom0) then
            dom = dom0/2.0_dp
         elseif(oga < 2.0_dp*delta0 + 10.0_dp*dom0) then
            dom = dom0/5.0_dp
         elseif(oga < 3.0_dp*delta0) then
            dom = dom0
         else
            dom = 0.1_dp*delta0
         endif

         oga = oga + dom
      end do

      do while(mod(ifreq,nprocs) /=0 ) ! See that number of frequencies divides evenly on processors 
         ifreq = ifreq+1
      enddo

      nfreq = ifreq
!      write(*,*) '# frequencies :', nfreq,ifreq
      allocate(frequencies(1:nfreq))

      oga = dom0/2.0_dp
      ifreq = 0

      do ifreq = 1,nfreq, 1

         if(oga < 5.0_dp*dom0) then
            dom = dom0/10.0_dp
         elseif(oga < 2.0_dp*delta0 - 10._dp*dom0) then
            dom = dom0/2.0_dp
         elseif(oga < 2.0_dp*delta0 + 10._dp*dom0) then
            dom = dom0/5.0_dp
         elseif(oga < 3.0_dp*delta0) then
            dom = dom0
         else
            dom = 0.1_dp*delta0
         endif

         frequencies(ifreq) = oga
!        write(2,*) ifreq,oga
         oga = oga + dom
      end do
!
!--------------------------------------------------------------------------
!
   end subroutine setfrec
!
!--------------------------------------------------------------------------
!
end module qc_cond
!
!-----------------------------------------------------------------------
