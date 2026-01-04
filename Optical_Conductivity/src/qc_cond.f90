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
   real(dp) :: sumsrate, oga, d2
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
         lcondv(:,i) = real(cond(:))/sumsrate + w*omega*aimag(cond(:))/2.0
         if(myrank == 0) write(*,1100) myrank, oga, lcondv(ic,i), lcondv(1,i), 1.0/(omega**2+d2*sumsrate**2)
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
            write(10,1000) oga, avcond(i) ,1.0/(omega**2+sumsrate**2)
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
      real(dp) :: epsl, epsu, dom, atol, rtol
      complex(dp), allocatable :: icond(:) 

      if (.not. allocated(icond)) allocate(icond(1:sx)) 

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
         call integrate_vec_adaptive(linear_resp, epsl, epsu, icond, atol, rtol, maxdepth)
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

      integer :: i, iep, iem, icont, icmp, iprint, imploop
      real(dp) :: d, thp, thm, err, epsp, epsm
      complex(dp) :: gold, told, xold, aold, gexp, wd, kshift
      complex(dp) :: aR(2,2), bR(2,2), tmp(2,2), cR(2,2), cX(2,2), gR(2,2), gX(2,2)

      complex(dp), allocatable :: NRp(:), NRm(:), NKp(:), NKm(:)
      complex(dp), allocatable :: gRpL(:), tRpL(:), XxpL(:), XapL(:)
      complex(dp), allocatable :: gRmL(:), tRmL(:), XxmL(:), XamL(:)
      complex(dp), allocatable :: gRp(:), tRp(:), Xxp(:), Xap(:)
      complex(dp), allocatable :: gRm(:), tRm(:), Xxm(:), Xam(:)
      complex(dp), allocatable :: R1p(:), R2p(:), R1m(:), R2m(:)
      complex(dp), allocatable :: X1p(:), X2p(:), X1m(:), X2m(:)

      complex(dp), allocatable :: dgRp(:), dgRm(:), dgXp(:), dgXm(:)
      complex(dp), allocatable :: dhRp(:), dhRm(:), dhXp(:), dhXm(:)
      complex(dp), allocatable :: dfRp(:), dfRm(:), dfXp(:), dfXm(:)
      complex(dp), allocatable :: dtRp(:), dtRm(:), dtXp(:), dtXm(:)
      complex(dp), allocatable :: agR(:), ahR(:), afR(:), atR(:)
      complex(dp), allocatable :: agX(:), ahX(:), afX(:), atX(:)

      complex(dp), allocatable :: n_gR0(:), n_gX0(:)
      complex(dp), allocatable :: n_gR1(:), n_gX1(:)
      complex(dp), allocatable :: n_gR2(:), n_gX2(:)
      complex(dp), allocatable :: n_gR3(:), n_gX3(:)

      complex(dp), allocatable :: o_gR0(:), o_gX0(:)
      complex(dp), allocatable :: o_gR1(:), o_gX1(:)
      complex(dp), allocatable :: o_gR2(:), o_gX2(:)
      complex(dp), allocatable :: o_gR3(:), o_gX3(:)

      complex(dp), allocatable :: lcurrR(:), lcurrK(:)
!
! --
!
      if(.not.allocated(NRp)) allocate(NRp(1:sx))
      if(.not.allocated(NRm)) allocate(NRm(1:sx))
      if(.not.allocated(NKp)) allocate(NKp(1:sx))
      if(.not.allocated(NKm)) allocate(NKm(1:sx))

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

      if(.not.allocated(dgRp)) allocate(dgRp(1:sx))
      if(.not.allocated(dhRp)) allocate(dhRp(1:sx))
      if(.not.allocated(dfRp)) allocate(dfRp(1:sx))
      if(.not.allocated(dtRp)) allocate(dtRp(1:sx))
      if(.not.allocated(dgRm)) allocate(dgRm(1:sx))
      if(.not.allocated(dhRm)) allocate(dhRm(1:sx))
      if(.not.allocated(dfRm)) allocate(dfRm(1:sx))
      if(.not.allocated(dtRm)) allocate(dtRm(1:sx))

      if(.not.allocated(dgXp)) allocate(dgXp(1:sx))
      if(.not.allocated(dhXp)) allocate(dhXp(1:sx))
      if(.not.allocated(dfXp)) allocate(dfXp(1:sx))
      if(.not.allocated(dtXp)) allocate(dtXp(1:sx))
      if(.not.allocated(dgXm)) allocate(dgXm(1:sx))
      if(.not.allocated(dhXm)) allocate(dhXm(1:sx))
      if(.not.allocated(dfXm)) allocate(dfXm(1:sx))
      if(.not.allocated(dtXm)) allocate(dtXm(1:sx))

      if(.not.allocated(agR)) allocate(agR(1:sx))
      if(.not.allocated(ahR)) allocate(ahR(1:sx))
      if(.not.allocated(afR)) allocate(afR(1:sx))
      if(.not.allocated(atR)) allocate(atR(1:sx))

      if(.not.allocated(agX)) allocate(agX(1:sx))
      if(.not.allocated(ahX)) allocate(ahX(1:sx))
      if(.not.allocated(afX)) allocate(afX(1:sx))
      if(.not.allocated(atX)) allocate(atX(1:sx))

      if(.not.allocated(n_gR0)) allocate(n_gR0(1:sx))
      if(.not.allocated(n_gR1)) allocate(n_gR1(1:sx))
      if(.not.allocated(n_gR2)) allocate(n_gR2(1:sx))
      if(.not.allocated(n_gR3)) allocate(n_gR3(1:sx))

      if(.not.allocated(n_gX0)) allocate(n_gX0(1:sx))
      if(.not.allocated(n_gX1)) allocate(n_gX1(1:sx))
      if(.not.allocated(n_gX2)) allocate(n_gX2(1:sx))
      if(.not.allocated(n_gX3)) allocate(n_gX3(1:sx))

      if(.not.allocated(o_gR0)) allocate(o_gR0(1:sx))
      if(.not.allocated(o_gR1)) allocate(o_gR1(1:sx))
      if(.not.allocated(o_gR2)) allocate(o_gR2(1:sx))
      if(.not.allocated(o_gR3)) allocate(o_gR3(1:sx))

      if(.not.allocated(o_gX0)) allocate(o_gX0(1:sx))
      if(.not.allocated(o_gX1)) allocate(o_gX1(1:sx))
      if(.not.allocated(o_gX2)) allocate(o_gX2(1:sx))
      if(.not.allocated(o_gX3)) allocate(o_gX3(1:sx))

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
      err =100.0
      imploop = 0  

      n_gR0 = czero
      n_gR1 = czero
      n_gR2 = czero
      n_gR3 = czero

      n_gX0 = czero
      n_gX1 = czero
      n_gX2 = czero
      n_gX3 = czero
!
!-- p.v > 0 --
!
      NRp(:)=w/((1.0+gr_1p(:,iep)*gr_2p(:,iep))*(1.0+gr_1p(:,iem)*gr_2p(:,iem)))
      NKp(:)=w/((1.0+gr_1p(:,iep)*gr_2p(:,iep))*conjg(1.0+gr_1p(:,iem)*gr_2p(:,iem)))
 
      R1p(:) = alpha1p(:,iep)+beta1p(:,iem)
      R2p(:) = alpha2p(:,iep)+beta2p(:,iem)

      X1p(:) = alpha1p(:,iep)-conjg(alpha1p(:,iem))-kshift
      X2p(:) = alpha2p(:,iep)-conjg(alpha2p(:,iem))-kshift
         
      gRpL(:) =-(gr_1p(:,iep)+gr_1p(:,iem))/R1p(:)
      tRpL(:) =-(gr_2p(:,iep)+gr_2p(:,iem))/R2p(:)

      XxpL(:) =-(1.0+gr_1p(:,iep)*conjg(gr_1p(:,iem)))/X1p(:)
      XapL(:) = (1.0+gr_2p(:,iep)*conjg(gr_2p(:,iem)))/X2p(:)
!
!-- p.v < 0 --
!
      NRm(:)=w/((1.0+gr_1m(:,iep)*gr_2m(:,iep))*(1.0+gr_1m(:,iem)*gr_2m(:,iem)))
      NKm(:)=w/((1.0+gr_1m(:,iep)*gr_2m(:,iep))*conjg(1.0+gr_1m(:,iem)*gr_2m(:,iem)))      

      R1m(:) = alpha1m(:,iep)+beta1m(:,iem)
      R2m(:) = alpha2m(:,iep)+beta2m(:,iem)

      X1m(:) = alpha1m(:,iep)-conjg(alpha1m(:,iem))-kshift
      X2m(:) = alpha2m(:,iep)-conjg(alpha2m(:,iem))-kshift

      gRmL(:) =-(gr_1m(:,iep)+gr_1m(:,iem))/R1m(:)
      tRmL(:) =-(gr_2m(:,iep)+gr_2m(:,iem))/R2m(:)

      XxmL(:) =-(1.0+gr_1m(:,iep)*conjg(gr_1m(:,iem)))/X1m(:)
      XamL(:) = (1.0+gr_2m(:,iep)*conjg(gr_2m(:,iem)))/X2m(:)
!
!-- Solve for the spatial dependence
!
      gold = gRpL(1)
      told = tRmL(1)
      xold = XxpL(1)
      aold = XamL(1)

      gexp   = exp(wd*R1p(1))
      gRp(1) = gold*gexp+(1.0-gexp)*gRpL(1)
      gexp   = exp(wd*X1p(1))
      Xxp(1) = xold*gexp+(1.0-gexp)*XxpL(1)
      gexp   = exp(wd*R2m(1))
      tRm(1) = told*gexp+(1.0-gexp)*tRmL(1)
      gexp   = exp(wd*X2m(1))
      Xam(1) = aold*gexp+(1.0-gexp)*XamL(1)

      gold = gRp(1)
      told = tRm(1)
      xold = Xxp(1)
      aold = Xam(1)
 
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

      gexp    = exp(wd*R1m(sx))
      gRm(sx) = gold*gexp+(1.0-gexp)*gRmL(sx)
      gexp    = exp(wd*X1m(sx))
      Xxm(sx) = xold*gexp+(1.0-gexp)*XxmL(sx)
      gexp    = exp(wd*R2p(sx))
      tRp(sx) = told*gexp+(1.0-gexp)*tRpL(sx)
      gexp    = exp(wd*X2p(sx))
      Xap(sx) = aold*gexp+(1.0-gexp)*XapL(sx)

      gold = gRm(sx)
      told = tRp(sx)
      xold = Xxm(sx)
      aold = Xap(sx)
 
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
!     do i = 1,sx ,1
!
!-- Make the diagonal part of the Green's function p.v > 0, a factor two less due to averaging over direction
!
!        dgRp = NRp(i)*((gr_2p(i,iep)+gr_2p(i,iem))*gRp(i) &
!                      +(gr_1p(i,iep)+gr_1p(i,iem))*tRp(i))
!        dgXp =-NKp(i)*((Xxp(i)-gr_1p(i,iep)*Xap(i)*conjg(gr_1p(i,iem))) &
!                      -(Xap(i)-gr_2p(i,iep)*Xxp(i)*conjg(gr_2p(i,iem))))
!        
!-- Make the diagonal part of the Green's function p.v > 0
!        
!        dgRm = NRm(i)*((gr_2m(i,iep)+gr_2m(i,iem))*gRm(i) &
!                      +(gr_1m(i,iep)+gr_1m(i,iem))*tRm(i))
!        dgXm =-NKm(i)*((Xxm(i)-gr_1m(i,iep)*Xam(i)*conjg(gr_1m(i,iem))) &
!                      -(Xam(i)-gr_2m(i,iep)*Xxm(i)*conjg(gr_2m(i,iem))))

!        lcurrR(i) = -(dgRp-dgRm)*thm+thp*conjg(dgRp-dgRm)
!        lcurrK(i) = -(thp-thm)*(dgXp-dgXm)
!     enddo      
!
!-- Make the Green's function p.v > 0, a factor two less due to averaging over direction
!
         dgRp(:) =-NRp(:)*(gRp(:)*gr_2p(:,iem)+gr_1p(:,iep)*tRp(:))
         dfRp(:) =-NRp(:)*(gRp(:)+gr_1p(:,iep)*tRp(:)*gr_1p(:,iem))
         dtRp(:) = NRp(:)*(tRp(:)+gr_2p(:,iep)*gRp(:)*gr_2p(:,iem))
         dhRp(:) = NRp(:)*(tRp(:)*gr_1p(:,iem)+gr_2p(:,iep)*gRp(:))

         dgXp(:) =-NKp(:)*(Xxp(:)-gr_1p(:,iep)*Xap(:)*conjg(gr_1p(:,iem)))
         dfXp(:) = NKp(:)*(gr_1p(:,iep)*Xap(:)-Xxp(:)*conjg(gr_1p(:,iem)))
         dtXp(:) = NKp(:)*(gr_2p(:,iep)*Xxp(:)-Xap(:)*conjg(gr_2p(:,iem)))
         dhXp(:) =-NKp(:)*(Xap(:)-gr_2p(:,iep)*Xxp(:)*conjg(gr_2p(:,iem)))
!
!-- Make the Green's funct:on p.v < 0
!
         dgRm(:) =-NRm(:)*(gRm(:)*gr_2m(:,iem)+gr_1m(:,iep)*tRm(:))
         dfRm(:) =-NRm(:)*(gRm(:)+gr_1m(:,iep)*tRm(:)*gr_1m(:,iem))
         dtRm(:) = NRm(:)*(tRm(:)+gr_2m(:,iep)*gRm(:)*gr_2m(:,iem))
         dhRm(:) = NRm(:)*(tRm(:)*gr_1m(:,iem)+gr_2m(:,iep)*gRm(:))

         dgXm(:) =-NKm(:)*(Xxm(:)-gr_1m(:,iep)*Xam(:)*conjg(gr_1m(:,iem)))
         dfXm(:) = NKm(:)*(gr_1m(:,iep)*Xam(:)-Xxm(:)*conjg(gr_1m(:,iem)))
         dtXm(:) = NKm(:)*(gr_2m(:,iep)*Xxm(:)-Xam(:)*conjg(gr_2m(:,iem)))
         dhXm(:) =-NKm(:)*(Xam(:)-gr_2m(:,iep)*Xxm(:)*conjg(gr_2m(:,iem)))
!
!-- Averages over directions, the minus sign due to the perturbation ~ p
!
         agR = (dgRp - dgRm)
         afR = (dfRp - dfRm)
         atR = (dtRp - dtRm)
         ahR = (dhRp - dhRm)

         agX = (dgXp - dgXm)
         afX = (dfXp - dfXm)
         atX = (dtXp - dtXm)
         ahX = (dhXp - dhXm)

         o_gR0 = n_gR0
         o_gR1 = n_gR1
         o_gR2 = n_gR2
         o_gR3 = n_gR3

         o_gX0 = n_gX0
         o_gX1 = n_gX1
         o_gX2 = n_gX2
         o_gX3 = n_gX3

         err = 0.0

         do i = 1, sx, 1

            gR(1,1) = agR(i)
            gR(1,2) = afR(i)
            gR(2,1) = atR(i)
            gR(2,2) = ahR(i)

            gX(1,1) = agX(i)
            gX(1,2) = afX(i)
            gX(2,1) = atX(i)
            gX(2,2) = ahX(i)

!-- potential scatterers

            aR(1,1) = gimpP(i,1)
            aR(1,2) = fimpP(i,1)
            aR(2,1) = timpP(i,1)
            aR(2,2) =-gimpP(i,1)

            bR(1,1) = gimpP(i,2)
            bR(1,2) = fimpP(i,2)
            bR(2,1) = timpP(i,2)
            bR(2,2) =-gimpP(i,2)

            tmp = czero
            tmp = matmul(gR,bR)
            cR = srate(1)*matmul(aR,tmp)

            bR = conjg(transpose(bR))
            tmp = czero
            tmp = matmul(gX,bR)
            cX = srate(1)*matmul(aR,tmp)

            n_gR0(i) = cR(1,1)
            n_gR1(i) = cR(1,2)
            n_gR2(i) = cR(2,1)
            n_gR3(i) = cR(2,2)

            n_gX0(i) = cX(1,1)
            n_gX1(i) = cX(1,2)
            n_gX2(i) = cX(2,1)
            n_gX3(i) = cX(2,2)

!-- magnetic scatterers

            aR(1,1) = gimpM(i,1)
            aR(1,2) = fimpM(i,1)
            aR(2,1) = timpM(i,1)
            aR(2,2) =-gimpM(i,1)

            bR(1,1) = gimpM(i,2)
            bR(1,2) = fimpM(i,2)
            bR(2,1) = timpM(i,2)
            bR(2,2) =-gimpM(i,2)

            tmp = czero
            tmp = matmul(gR,bR)
            cR = srate(2)*matmul(aR,tmp)

            bR = conjg(transpose(bR))
            tmp = czero
            tmp = matmul(gX,bR)
            cX = srate(2)*matmul(aR,tmp)

            n_gR0(i) = n_gR0(i)+cR(1,1)
            n_gR1(i) = n_gR1(i)+cR(1,2)
            n_gR2(i) = n_gR2(i)+cR(2,1)
            n_gR3(i) = n_gR3(i)+cR(2,2)

            n_gX0(i) = n_gX0(i)+cX(1,1)
            n_gX1(i) = n_gX1(i)+cX(1,2)
            n_gX2(i) = n_gX2(i)+cX(2,1)
            n_gX3(i) = n_gX3(i)+cX(2,2)

            err = err+abs(n_gR0(i)-o_gR0(i))**2+abs(n_gR1(i)-o_gR1(i))**2+ &
                      abs(n_gR2(i)-o_gR2(i))**2+abs(n_gR3(i)-o_gR3(i))**2

            err = err+abs(n_gX0(i)-o_gX0(i))**2+abs(n_gX1(i)-o_gX1(i))**2+ &
                      abs(n_gX2(i)-o_gX2(i))**2+abs(n_gX3(i)-o_gX3(i))**2
         enddo
         err = sqrt(err)

!-- Current contribution
!
      lcurrR =  (dgRp-dhRp+(dgRm-dhRm))*thm-thp*conjg(dgRp-dhRp+(dgRm-dhRm))
      lcurrK =  (thp-thm)*(dgXp-dhXp+(dgXm-dhXm))

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
