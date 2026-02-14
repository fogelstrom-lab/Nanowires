!-------------------------------------------------------------------
! Modernized module with procedures from qc.f and riccati.f
!-------------------------------------------------------------------
module qc_mod
  use mpi_f08
  use qc_data
  use write_stuff
  use qc_farm
! use qc_cond
  use iso_fortran_env, only: dp => real64
  use math_utils, only: cdiv   ! <--- import the function
  implicit none
  private
  public :: get_input, get_equilibrium, check_tilde
contains
!
!---------------------------------------------------------------------
!
   subroutine get_input
!
!---------------------------------------------------------------------
!
      complex(dp) :: op
      real(dp) :: x, dx, halfway, uu, testi, x1, x2, p1, p2, ex, em, f, taub
      real(dp), allocatable :: testf(:)
      integer :: ii, i1, i2

      delta0=0.2808

      if(myrank ==0) then
         read(*,*) temp
         write(*,1510)' The temperature                         : ',temp
         temp = temp*delta0
         read(*,*) Rx
         write(*,1510)' The total wire length (in xi:s)         : ',Rx
         read(*,*) Ngap
         write(*,1510)' Length of the uncovered part (in xi:s)  : ',Ngap
         if(Ngap > Rx) stop ' The uncovered part is longer than the nanowire'
         read(*,*) dx
         write(*,1510)' The discretisation step (in xi:s)       : ',dx
         read(*,*) phase
         write(*,1510)' The phase difference (in pi)            : ',phase
         phase=phase*pi
         read(*,*) tau(1) 
         write(*,1510)' The pot. impurity scatt rate (in deltas): ',tau(1)
         read(*,*) sigma(1)
         write(*,1510)' The pot. impurity scatt. phase shift    : ',sigma(1)
         read(*,*) tau(2) 
         write(*,1510)' The mag. impurity scatt rate (in deltas): ',tau(2)
         read(*,*) sigma(2)
         write(*,1510)' The mag. impurity scatt. phase shift    : ',sigma(2)
         read(*,*) emax
         write(*,1510)' Energy window (in deltas)               : ',emax
         read(*,*) iemax
         write(*,1500)' Energy steps in this window             : ',iemax
         read(*,*) error_tol
         write(*,1520)' The error tolerance                     : ',error_tol
      endif

      call MPI_Bcast(temp      , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(Rx        , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dx        , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(Ngap      , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(phase     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(tau       , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(sigma     , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(Emax      , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iemax     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(error_tol , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      
!
! - set up the spatial grid and the order parameter profile
!
      Ngap = 0.5_dp*Ngap

      sx      = int(Rx/dx)        
      if(mod(sx,2).ne. 0) sx = sx+1   ! even number of space steps
      if(myrank==0) write(*,*) max(1,sx/10)

      x       = Rx/float(sx)
      halfway = Rx/2.0
      if(myrank ==0) write(*,*) ' # steps', sx,' per length ',Rx,x

      if(.not.allocated(grid)) allocate(grid(1:sx))
      do ii=1,sx,1
         grid(ii)=(float(ii)-0.5)*x-halfway
      enddo

      op   = delta0*exp(w*phase)
      DelR =       op 
      DelL =  delta0

      if(.not.allocated(delx)) allocate(delx(1:sx))
      delx = czero
      delx(1) = delL 
      delX(sx) = delR

!       do ii=1,sx
!          if(grid(ii) < -Ngap) then
!             delx(ii)= DelL
!          endif
!          if(grid(ii) > Ngap) then
!             delx(ii)= DelR 
!          endif
!          if(myrank == 0) write(2,1000) grid(ii),delx(ii)
!       enddo
!       if(myrank == 0) close(2)
!
! - set up the energy grid
!
      if(.not.allocated(er)) allocate(er(-iemax:iemax))
      etaR = delta0/100.0
      de   = 0.495*pi/iemax
      Ec   = delta0*emax/tan(0.495*pi)

      do ii=-iemax,iemax,1
         er(ii)=Ec*tan(de*ii)
      enddo
!
! - set up the impurity properties
!
      uu = 0.0_dp
      if(abs(sigma(1)) .gt. 1.0e-10) uu = sqrt((1.0-sigma(1))/sigma(1)) 
      srate(0) = uu 
      srate(1) = tau(1)*delta0
      srate(2) = tau(2)*delta0

      if(.not.allocated(uimp)) allocate(uimp(1:sx, -iemax:iemax))
      if(.not.allocated(gimp)) allocate(gimp(1:sx, -iemax:iemax))
      if(.not.allocated(fimp)) allocate(fimp(1:sx, -iemax:iemax))
      if(.not.allocated(timp)) allocate(timp(1:sx, -iemax:iemax))
      uimp = czero
      gimp = czero
      fimp = czero
      timp = czero

      if(.not.allocated(uimp0)) allocate(uimp0(1:sx))
      if(.not.allocated(gimp0)) allocate(gimp0(1:sx))
      if(.not.allocated(fimp0)) allocate(fimp0(1:sx))
      if(.not.allocated(timp0)) allocate(timp0(1:sx))
      gimp0 = czero
      fimp0 = czero
      timp0 = czero

      if(.not.allocated(gimpP)) allocate(gimpP(1:sx,2))
      if(.not.allocated(uimpP)) allocate(uimpP(1:sx,2))
      if(.not.allocated(fimpP)) allocate(fimpP(1:sx,2))
      if(.not.allocated(timpP)) allocate(timpP(1:sx,2))
      gimpP = czero
      uimpP = czero
      fimpP = czero
      timpP = czero

      if(.not.allocated(gimpM)) allocate(gimpM(1:sx,2))
      if(.not.allocated(fimpM)) allocate(fimpM(1:sx,2))
      if(.not.allocated(timpM)) allocate(timpM(1:sx,2))
      gimpM = czero
      fimpM = czero
      timpM = czero

!
!-- Leading-order quantities that are saved for later
!
      if(.not. allocated(gr_1)) allocate(gr_1(1:sx,2))
      if(.not. allocated(gr_2)) allocate(gr_2(1:sx,2))

      if(.not.allocated(avj)) allocate(avj(1:sx, -iemax:iemax))
      if(.not.allocated(avg)) allocate(avg(1:sx, -iemax:iemax))
      if(.not.allocated(avf)) allocate(avf(1:sx, -iemax:iemax))
      if(.not.allocated(avt)) allocate(avt(1:sx, -iemax:iemax))

      avj = czero
      avg = czero
      avf = czero
      avt = czero
!
!-- For the linear response part (1 refers to non-tilded, 2 to tilded quanities, p/m is p.v >/< 0 
!
      if(.not.allocated(gr_1p)) allocate(gr_1p(1:sx,2))
      if(.not.allocated(gr_2p)) allocate(gr_2p(1:sx,2))

      if(.not.allocated(alpha1p))  allocate(alpha1p(1:sx,2))
      if(.not.allocated(alpha2p))  allocate(alpha2p(1:sx,2))
      if(.not.allocated( beta1p))  allocate( beta1p(1:sx,2))
      if(.not.allocated( beta2p))  allocate( beta2p(1:sx,2))

      if(.not.allocated(gr_1m)) allocate(gr_1m(1:sx,2))
      if(.not.allocated(gr_2m)) allocate(gr_2m(1:sx,2))

      if(.not.allocated(alpha1m))  allocate(alpha1m(1:sx,2))
      if(.not.allocated(alpha2m))  allocate(alpha2m(1:sx,2))
      if(.not.allocated( beta1m))  allocate( beta1m(1:sx,2))
      if(.not.allocated( beta2m))  allocate( beta2m(1:sx,2))
!
! small test of the interpolation skeme
!
      if(myrank == 0) then
         allocate(testf(-iemax:iemax))
         testf = 0.0

         do ii = -iemax, iemax
            p1 = exp(-0.3*(er(ii)-3.0*delta0)**2)
            p2 = exp(-2.1*(er(ii)+2.0*delta0)**2)
            testf(ii) = p1 + p2
            write(1,*) ii,er(ii)
            write(3,1000) er(ii),testf(ii),p1,p2
         enddo
       
         em = delta0*emax
         ex = -1.5*em
         do while(ex < 1.5*em) 
            if(abs(ex) < em) then
               x= atan(ex/Ec)
               if(x .ge. 0) then
                  i1 = int(x/de)
                  i2 = i1 + 1
                  x1 = er(i1)
                  x2 = er(i2)
                  p1 = (x2-ex)/(x2-x1)
                  p2 = (ex-x1)/(x2-x1)
               else
                  i2 = -int(abs(x)/de)
                  i1 = i2 - 1
                  x1 = er(i1)
                  x2 = er(i2)
                  p1 = (x2-ex)/(x2-x1)
                  p2 = (ex-x1)/(x2-x1)
               endif
               testi = p1*testf(i1)+p2*testf(i2)
            else
               if(ex > 0) i1 = iemax
               if(ex < 0) i1 =-iemax
               testi = testf(i1)
            endif
            write(4,1000) ex,testi
            ex = ex + 0.01*delta0
         enddo
         close(1)
         close(3)
         close(4)
      endif

 1000 format(9(1x,e12.6))
 1500 format(a,i5)
 1510 format(a,f8.4)
 1520 format(a,e12.4)
!
!---------------------------------------------------------------------
!
   end subroutine get_input
!
!---------------------------------------------------------------------
!
   subroutine check_tilde
!
!---------------------------------------------------------------------
!
      integer :: i,ie
      real(dp) :: errt

      errt=0.0
      i=max(sx/3,1)
      do ie=-iemax,iemax,1
         if(ie.ne.0) then
            errt=errt+abs(gimp(i,ie)+conjg(gimp(i,-ie))) &
                     +abs(fimp(i,ie)-conjg(timp(i,-ie))) &
                     +abs(timp(i,ie)-conjg(fimp(i,-ie)))
         endif
      enddo
      write(*,*) ' Tilde symmetry to ', errt
!
!---------------------------------------------------------------------
!
   end subroutine check_tilde
!
!---------------------------------------------------------------------
!
   subroutine get_equilibrium(itype)
!
!---------------------------------------------------------------------
!
      integer :: counter, itype, ien, i, idum, iprint
      real(dp) :: err, gp, energy, etaR_org

      avj = (0.0_dp, 0.0_dp)
      avg = (0.0_dp, 0.0_dp)
      avf = (0.0_dp, 0.0_dp)
      avt = (0.0_dp, 0.0_dp)

      etaR_org = etaR

      do while (etaR >= etaR_org/(2048.0))
         call get_ses(counter,err)
         gp = etaR / delta0
         if(myrank == 0) write(*,1000) ' Did round with etaR = ',gp,counter,' tau =', tau(1),tau(2),' Ngap =',Ngap*2.0_dp
         etaR = etaR/2.0
      enddo
      etaR = etaR*2.0
 
      if(myrank == 0) then
         call check_tilde
         call write_cuts(itype)
      endif

      deallocate(avj,avg,avf,avt)

      idum = 0
      iprint = 1 
      ex_old = 3000
      do ien = -iemax + myrank, iemax, 200*nprocs
         if(ien .le. 0) energy = 0.5*(er(ien+1)+er(ien))
         if(ien .ge. 0) energy = 0.5*(er(ien)+er(ien-1))
         i = 1
         call get_imp_at_e(energy,i,idum,iprint)     
         if(myrank == 0) write(*,*) ien, energy/delta0, idum
      end do

      if(myrank == 0 .and. iprint== 0) then
         close(500)
         close(501)
         close(502)
         close(503)
         close(600)
         close(601)
         close(602)
         close(603)
         close(700)
         close(701)
         close(702)
         close(703)
         close(800)
         close(801)
         close(802)
         close(803)
      endif

 1000 format(a,e12.6,' ( iterations=',i8,')',a,2(1x,f8.4),a,f8.4)
!1100 format(a,2(1x,f7.4),a,10(1x,e14.8))
!1200 format(10(1x,e14.8))
! 
!---------------------------------------------------------------------
!
   end subroutine get_equilibrium
!
!---------------------------------------------------------------------
!
end module qc_mod
!
!--------------------------------------------------------------------
