Module qc_cond
  use mpi_f08
  use qc_data
  use qc_farm
  use iso_fortran_env, only: dp => real64
  use mpi_utils, only: sync_cols_gather_indexed_rpairs
  use math_utils, only: cdiv   ! <--- import the function
  use integration_vec,  only: integrate_vec_adaptive

  implicit none
  integer :: nfreq, counter, icall_lin
  real(dp) :: omega, omh
  real(dp), allocatable :: frequencies(:)
  complex(dp), allocatable :: cond(:)
  public :: qc_conductance, integral, linear_resp
contains
!-----------------------------------------------------------------------
!
   subroutine qc_conductance
!
!-----------------------------------------------------------------------
!
   integer :: i, ir, ic
   complex(dp), allocatable :: condv(:,:), lcondv(:,:), avcond(:) 
   real(dp) :: sumsrate, oga
   integer :: ulog, ios, pstep
   character(len=128) :: fname

      sumsrate = tau(1)+2.0*tau(2)
      call setfrec(sumsrate)
      allocate(condv(1:sx,1:nfreq),lcondv(1:sx,1:nfreq))
      allocate(cond(1:sx))
      allocate(avcond(1:nfreq))

      condv  =  (0.0_dp, 0.0_dp)
      lcondv =  (0.0_dp, 0.0_dp)
      ic =  sx/2


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
         lcondv(:,i) = real(cond(:))/sumsrate + w*omega*aimag(cond(:))/2.0
         if(myrank == 0) write(*,1100) myrank, oga, lcondv(ic,i), lcondv(1,i)
         write(ulog,*) oga, counter, icall_lin ! logg how many self-energy iterations per frecquency
      enddo
      close(ulog)

      call MPI_Reduce(lcondv, condv, sx*nfreq, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
          write(*,*) "rank", myrank, "MPI_Reduce failed ierr=", ierr
          call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
      end if
      
      if(myrank == 0) then

         open(10,file="avcond.dat", status="replace", action="write", iostat=ios)
         if (ios/=0) call MPI_Abort(MPI_COMM_WORLD, 10, ierr)
         open(11,file="cond_real.dat", status="replace", action="write", iostat=ios)
         if (ios/=0) call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
         open(12,file="cond_imag.dat", status="replace", action="write", iostat=ios)
         if (ios/=0) call MPI_Abort(MPI_COMM_WORLD, 12, ierr)

         avcond(:)=sum(condv, dim=1)/sx
         pstep = max(1,sx/10)
         do i=1,nfreq
            omega = frequencies(i)
            oga = omega/delta0
            write(10,1000) oga, avcond(i)
            write(11,1000) oga, (real(condv(ir,i)),ir=1,sx,pstep)
            write(12,1000) oga,(aimag(condv(ir,i)),ir=1,sx,pstep)
         enddo
         close(10)
         close(11)
         close(12)
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
      real(dp) :: epsl, epsu, dom, tol
      complex(dp), allocatable :: icond(:) 

      if (.not. allocated(icond)) allocate(icond(1:sx)) 

      tol = 1.0e-7_dp
      maxdepth = 30
      icond = (0.0_dp, 0.0_dp)
      cond  = (0.0_dp, 0.0_dp)
      dom = max(omega,5.0*etaR)

      iestep = iemax/100
!     if(myrank == 0) write(*,*) ' integrating :', omega

      do ie = -iemax, iemax-iestep, iestep

         epsl = er(ie)
         epsu = er(ie+iestep)
         call integrate_vec_adaptive(linear_resp, epsl, epsu, icond, tol, maxdepth)
         cond = cond + icond
      end do

      cond = cond/dom
!     if(myrank == 0) write(*,1000) omega, cond(sx/2)

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
      complex(dp), intent(out) :: icond(:)

      integer :: i, iep, iem, icont, icmp, iprint
      real(dp) :: d, thp, thm, err, epsp, epsm
      complex(dp) :: NRpp, NRpm,  NRmp, NRmm
      complex(dp) :: dgRp, dgRm, dgXp, dgXm, kshift
      complex(dp) :: gold, told, xold, aold, gexp, wd

      complex(dp), allocatable :: gRpL(:), tRpL(:), XxpL(:), XapL(:)
      complex(dp), allocatable :: gRmL(:), tRmL(:), XxmL(:), XamL(:)
      complex(dp), allocatable :: gRp(:), tRp(:), Xxp(:), Xap(:)
      complex(dp), allocatable :: gRm(:), tRm(:), Xxm(:), Xam(:)
      complex(dp), allocatable :: R1p(:), R2p(:), R1m(:), R2m(:)
      complex(dp), allocatable :: X1p(:), X2p(:), X1m(:), X2m(:)

      complex(dp), allocatable :: lcurrR(:), lcurrK(:)
!
! --
!
      if(.not.allocated(gRp)) allocate(gRp(1:sx))
      if(.not.allocated(gRm)) allocate(gRm(1:sx))
      if(.not.allocated(tRp)) allocate(tRp(1:sx))
      if(.not.allocated(tRm)) allocate(tRm(1:sx))
      if(.not.allocated(Xxp)) allocate(Xxp(1:sx))
      if(.not.allocated(Xxm)) allocate(Xxm(1:sx))
      if(.not.allocated(Xap)) allocate(Xap(1:sx))
      if(.not.allocated(Xam)) allocate(Xam(1:sx))

      if(.not.allocated(gRpL)) allocate(gRpL(1:sx))
      if(.not.allocated(gRmL)) allocate(gRmL(1:sx))
      if(.not.allocated(tRpL)) allocate(tRpL(1:sx))
      if(.not.allocated(tRmL)) allocate(tRmL(1:sx))
      if(.not.allocated(XxpL)) allocate(XxpL(1:sx))
      if(.not.allocated(XxmL)) allocate(XxmL(1:sx))
      if(.not.allocated(XapL)) allocate(XapL(1:sx))
      if(.not.allocated(XamL)) allocate(XamL(1:sx))

      if(.not.allocated(R1p)) allocate(R1p(1:sx))
      if(.not.allocated(R1m)) allocate(R1m(1:sx))
      if(.not.allocated(R2p)) allocate(R2p(1:sx))
      if(.not.allocated(R2m)) allocate(R2m(1:sx))
      if(.not.allocated(X1p)) allocate(X1p(1:sx))
      if(.not.allocated(X1m)) allocate(X1m(1:sx))
      if(.not.allocated(X2p)) allocate(X2p(1:sx))
      if(.not.allocated(X2m)) allocate(X2m(1:sx))

      if(.not.allocated(lcurrR)) allocate(lcurrR(1:sx))
      if(.not.allocated(lcurrK)) allocate(lcurrK(1:sx))
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
      wd=w*abs(grid(1)-grid(2))/2.0
      kshift = 2.0*w*etaR
!
!-- p.v > 0 --
!
      R1p(:) = alpha1p(:,iep)+beta1p(:,iem)
      R2p(:) = alpha2p(:,iep)+beta2p(:,iem)

      gRp(:) =-(gr_1p(:,iep)+gr_1p(:,iem))/R1p(:)
      tRp(:) =-(gr_2p(:,iep)+gr_2p(:,iem))/R2p(:)

      X1p(:) = alpha1p(:,iep)-conjg(alpha1p(:,iem))-kshift
      X2p(:) = alpha2p(:,iep)-conjg(alpha2p(:,iem))-kshift

      Xxp(:) = (1.0+gr_1p(:,iep)*conjg(gr_1p(:,iem)))/X1p(:)
      Xap(:) =-(1.0+gr_2p(:,iep)*conjg(gr_2p(:,iem)))/X2p(:)

      gRpL(:) = gRp(:)
      tRpL(:) = tRp(:)
         
      XxpL(:) = Xxp(:)
      XapL(:) = Xap(:)
!
!-- p.v < 0 --
!
      R1m(:) = alpha1m(:,iep)+beta1m(:,iem)
      R2m(:) = alpha2m(:,iep)+beta2m(:,iem)

      gRm(:) = (gr_1m(:,iep)+gr_1m(:,iem))/R1m(:)
      tRm(:) = (gr_2m(:,iep)+gr_2m(:,iem))/R2m(:)

      X1m(:) = alpha1m(:,iep)-conjg(alpha1m(:,iem))-kshift
      X2m(:) = alpha2m(:,iep)-conjg(alpha2m(:,iem))-kshift

      Xxm(:) =-(1.0+gr_1m(:,iep)*conjg(gr_1m(:,iem)))/X1m(:)
      Xam(:) = (1.0+gr_2m(:,iep)*conjg(gr_2m(:,iem)))/X2m(:)

      gRmL(:) = gRm(:)
      tRmL(:) = tRm(:)
      
      XxmL(:) = Xxm(:)
      XamL(:) = Xam(:)
!
!-- Solve for the spatial dependence
!
      gold = gRpL(1)
      told = tRmL(1)
      xold = XxpL(1)
      aold = XamL(1)
 
      do i=2,sx,1
         gexp   = exp(wd*R1p(i-1))
         gold   = gold*gexp+(1.0-gexp)*gRpL(i-1)
         gexp   = exp(wd*R1p(i))
         gold   = gold*gexp+(1.0-gexp)*gRpL(i)
         gRp(i) = gold
 
         gexp   = exp(wd*X1p(i-1))
         xold   = xold*gexp+(1.0-gexp)*XxpL(i-1)
         gexp   = exp(wd*X1p(i))
         xold   = xold*gexp+(1.0-gexp)*XxpL(i)
         Xxp(i) = xold
 
         gexp   = exp(wd*R2m(i-1))
         told   = told*gexp+(1.0-gexp)*tRmL(i-1)
         gexp   = exp(wd*R2m(i))
         told   = told*gexp+(1.0-gexp)*tRmL(i)
         tRm(i) = told
 
         gexp   = exp(wd*X2m(i-1))
         aold   = aold*gexp+(1.0-gexp)*XamL(i-1)
         gexp   = exp(wd*X2m(i))
         aold   = aold*gexp+(1.0-gexp)*XamL(i)
         Xam(i) = aold
      enddo
 
      gold = gRmL(sx)
      told = tRpL(sx)
      xold = XxmL(sx)
      aold = XapL(sx)
 
      do i=sx-1,1,-1
         gexp   = exp(wd*R1m(i+1))
         gold   = gold*gexp+(1.0-gexp)*gRmL(i+1)
         gexp   = exp(wd*R1m(i))
         gold   = gold*gexp+(1.0-gexp)*gRmL(i)
         gRm(i) = gold
 
         gexp   = exp(wd*X1m(i+1))
         xold   = xold*gexp+(1.0-gexp)*XxmL(i+1)
         gexp   = exp(wd*X1m(i))
         xold   = xold*gexp+(1.0-gexp)*XxmL(i)
         Xxm(i) = xold
 
         gexp   = exp(wd*R2p(i+1))
         told   = told*gexp+(1.0-gexp)*tRpL(i+1)
         gexp   = exp(wd*R2p(i))
         told   = told*gexp+(1.0-gexp)*tRpL(i)
         tRp(i) = told
 
         gexp   = exp(wd*X2p(i+1))
         aold   = aold*gexp+(1.0-gexp)*XapL(i+1)
         gexp   = exp(wd*X2p(i))
         aold   = aold*gexp+(1.0-gexp)*XapL(i)
         Xap(i) = aold
      enddo
!
!-- Make the current
!
      do i = 1,sx ,1
!
!-- Make the diagonal part of the Green's function p.v > 0, a factor two less due to averaging over direction
!
         NRpp=0.25/(1.0+gr_1p(i,iep)*gr_2p(i,iep))
         NRpm=0.25/(1.0+gr_1p(i,iem)*gr_2p(i,iem))

         dgRp = w*NRpp*((gr_2p(i,iep)+gr_2p(i,iem))*gRp(i) &
                       +(gr_1p(i,iep)+gr_1p(i,iem))*tRp(i))*NRpm
         dgXp =-w*NRpp*((Xxp(i)-gr_1p(i,iep)*Xap(i)*conjg(gr_1p(i,iem))) &
                       -(Xap(i)-gr_2p(i,iep)*Xxp(i)*conjg(gr_2p(i,iem))))*conjg(NRpm)
!
!-- Make the diagonal part of the Green's function p.v > 0
!
         NRmp=0.25/(1.0+gr_1m(i,iep)*gr_2m(i,iep))
         NRmm=0.25/(1.0+gr_1m(i,iem)*gr_2m(i,iem))

         dgRm = w*NRmp*((gr_2m(i,iep)+gr_2m(i,iem))*gRm(i) &
                       +(gr_1m(i,iep)+gr_1m(i,iem))*tRm(i))*NRmm
         dgXm =-w*NRmp*((Xxm(i)-gr_1m(i,iep)*Xam(i)*conjg(gr_1m(i,iem))) &
                       -(Xam(i)-gr_2m(i,iep)*Xxm(i)*conjg(gr_2m(i,iem))))*conjg(NRmm)

         lcurrR(i) = -(dgRp-dgRm)*thm+thp*conjg(dgRp-dgRm)
         lcurrK(i) = -(thp-thm)*(dgXp-dgXm) 
      enddo         

      icond = lcurrR+lcurrK

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
