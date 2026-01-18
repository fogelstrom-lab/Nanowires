!-------------------------------------------------------------------
! Modernized module with procedures from qc.f and riccati.f
!-------------------------------------------------------------------
module qc_farm
  use mpi_f08
  use qc_data
  use mpi_utils, only: sync_cols_gather_indexed_rpairs
  use write_stuff
  use iso_fortran_env, only: dp => real64
  use math_utils, only: cdiv   ! <--- import the function
  implicit none
  private
  public :: get_ses, get_rimps, get_gammas, get_averages, get_imp_at_e
  
contains
!
!---------------------------------------------------------------------
!
   subroutine get_ses(counter,err)
!
!---------------------------------------------------------------------
!
      integer, intent(inout) :: counter
      real(dp), intent(inout) :: err
      integer :: ii, icmptyp, icmp, ien

      counter = 0
      icmptyp = 1

      do ii = -iemax + myrank, iemax, nprocs
         ien = ii
         en   = er(ien)+w*etaR
         call get_rimps(ien, err, icmp, icmptyp)
         counter = counter + icmp
         if (icmp .ge. 50000) write(*,*) ' Fail (',myrank,')', ien, err, aimag(lgimp(ii)), icmp
      end do

      call MPI_Allreduce(MPI_IN_PLACE, counter, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

!-- MPI synchronize arrays so all ranks have a full copy for subsequent rounds

!-- Averaged Green's function components and equilibrium current contribution (e)

      call MPI_Reduce(lavj, avj, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(lavg, avg, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(lavf, avf, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(lavt, avt, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

!-- Impurity self energies (e)

      call MPI_Reduce(lgimp, gimp, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(lfimp, fimp, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(ltimp, timp, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(luimp, uimp, 2*iemax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!-- 
 1000 format(1x,e12.6,3(1x,e12.6))
!---------------------------------------------------------------------
!
   end subroutine get_ses
!
!---------------------------------------------------------------------
!
   subroutine get_rimps(ien, err, icmp, icmptyp)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: ien, icmptyp
      integer, intent(out) :: icmp
      real(dp), intent(inout) :: err 
 
      complex(dp) :: gimp1,fimp1,timp1,dh
      integer :: i1, i2, ic, ioshift    
      real(dp) :: x, x1, x2, p1, p2, ex, em, ah 
 
      save

      icmp = 0
      ic = sx/2
      ioshift = ien

      if(icmptyp == 1) then

         gimp0 = gimp(ien)
         fimp0 = fimp(ien)
         timp0 = timp(ien)

      elseif(icmptyp==2) then

         ex = real(en)                     ! current frequency     ex = Ec*tan(x)
         em = delta0*emax                  ! max frequency stored  Ec = emax*delta0/tan(de*iemax)

         if(abs(ex) .le. em) then
            x = atan(abs(ex)/Ec)                ! running variable x_i1 < x < x_i2
            if(ex .ge. 0) then
               i1 = int(x/de)              ! left bound index  x_i1 = Ec*tan(de * i1)
               i2 = i1 + 1                 ! right bound index x_i2 = Ec*tan(de * i2)
               x1 = er(i1)
               x2 = er(i2)
               p1 = (x2-ex)/(x2-x1)
               p2 = (ex-x1)/(x2-x1)
            else
               i2 = -int(x/de)             ! right bound
               i1 = i2 - 1                 ! left bound
               x1 = er(i1)
               x2 = er(i2)
               p1 = (x2-ex)/(x2-x1)
               p2 = (ex-x1)/(x2-x1)
            endif
            gimp0 = p1*gimp(i1)+p2*gimp(i2)
            fimp0 = p1*fimp(i1)+p2*fimp(i2)
            timp0 = p1*timp(i1)+p2*timp(i2)
            
!           if(myrank == 0) write(29,1000) ex, gimp0(ic), fimp0(ic), timp0(ic), icmp
         else
            if(ex > 0) i1 = iemax
            if(ex < 0) i1 =-iemax
            gimp0 = gimp(i1)
            fimp0 = fimp(i1)
            timp0 = timp(i1)
         endif
      else
         stop ' Fucking up again...'
      endif

      err = 1000.0

      do while(err > error_tol .and. icmp < 10000) 

         gimp1 = gimp0
         fimp1 = fimp0
         timp1 = timp0

         call get_gammas
         call get_averages(ioshift,icmptyp)
      
         err=0.0
         dh=gimp0-gimp1
         ah=     conjg(gimp1)*gimp1
         err=    conjg(dh)*dh/(ah+error_tol)

         dh=fimp0-fimp1
         ah=     conjg(fimp1)*fimp1
         err=err+conjg(dh)*dh/(ah+error_tol)

         dh=timp0-timp1
         ah=     conjg(timp1)*timp1
         err=err+conjg(dh)*dh/(ah+error_tol)

         err=sqrt(err)

         icmp=icmp+1
      enddo
!
! -- Save the leading-order impurity self-energies
!
      if(icmptyp == 1) then
         lgimp(ien) = gimp0
         lfimp(ien) = fimp0
         ltimp(ien) = timp0
         luimp(ien) = uimp0
      endif
!
! -- Save the t-matrices for linear response
!
      if(icmptyp == 2) then
         if(abs(srate(1)).gt. 1e-14) then
             gimpP(ioshift) = gimpP(ioshift)/srate(1)
             fimpP(ioshift) = fimpP(ioshift)/srate(1)
             timpP(ioshift) = timpP(ioshift)/srate(1)
          endif

          if(abs(srate(2)).gt. 1e-14) then
             gimpM(ioshift) = gimpM(ioshift)/srate(2)
             fimpM(ioshift) = fimpM(ioshift)/srate(2)
             timpM(ioshift) = timpM(ioshift)/srate(2)
          endif
      endif

 1000 format(7(1x,e10.4),1x,i7)
!
!---------------------------------------------------------------------
!
   end subroutine get_rimps
!
!---------------------------------------------------------------------
!
   subroutine get_averages(ien,icmptyp)
!
!----------------------------------------------------------------------
!
      integer, intent(in) :: icmptyp, ien
      integer :: it, ii, ioshift
      real(dp) :: magn, magd, half
      complex(dp) :: g1, g2, g3, g32, norm, dig, cj, Denom, D1, pg
!
!----------------------------------------------------------------------
!
!-- Perform the required integrations
!
      if(icmptyp == 1) ioshift = 1
      if(icmptyp == 2) ioshift = ien

      half = 0.5_dp
      magd = 4.0_dp*sigma(2)*(1.0_dp-sigma(2))
      magn = 2.0_dp*(1.0_dp-sigma(2))*srate(2)

      it=1
      pg = gr_1(it)*gr_2(it)             

      norm= -w/(cone+pg)
      g1  = 2.0_dp*gr_1(it)*norm           
      g2  = 2.0_dp*gr_2(it)*norm
      dig =       (cone-pg)*norm
 
      g3 = dig
      cj = dig
 
      it=2
      pg = gr_1(it)*gr_2(it)

      norm=-w/(cone+pg)
      g1  = g1+2.0_dp*gr_1(it)*norm
      g2  = g2+2.0_dp*gr_2(it)*norm
      dig =          (cone-pg)*norm

      g1 = g1*half
      g2 = g2*half
      g3 = (g3+dig)*half

      cj= cj-dig
!
! -- If iterating for equilibrium self enegies at given energy-grid point save the averages
!
      if(icmptyp == 1) then
         lavj(ien) = cj
         lavg(ien) = g3
         lavf(ien) = g1
         lavt(ien) = g2
      endif

      g32 = g3*g3
      norm = g32+g1*g2
!
! -- Potential scattering
!
      Denom = cone/(cone-sigma(1)*(cone+norm))
      cj = srate(1)*Denom
      fimpP(ioshift) = cj*g1
      timpP(ioshift) =-cj*g2
      gimpP(ioshift) = cj*g3

      uimp0 = srate(0)*cj
!
! -- Magnetic scattering
!
      D1 = (cone-sigma(2)*(cone-norm))
      Denom = cone/(D1*D1-magd*g32)
      cj = srate(2)*D1*Denom
      pg = magn*Denom

      fimpM(ioshift) = -cj*g1
      timpM(ioshift) =  cj*g2
      gimpM(ioshift) = -(cj-pg)*g3

      gimp0 = gimpP(ioshift) + gimpM(ioshift)
      fimp0 = fimpP(ioshift) + fimpM(ioshift)
      timp0 = timpP(ioshift) + timpM(ioshift)
!
!---------------------------------------------------------------------
!
   end subroutine get_averages
!
!---------------------------------------------------------------------
!
   subroutine get_gammas
!
!---------------------------------------------------------------------
!
      integer :: i, it
      real(dp) :: dx
      complex(dp) :: zen, del1, del2, lam
      complex(dp) :: ga0, gb0
      complex(dp) :: gold, w2, gC, g0, gaL, gaR, gbL, gbR
      save
!
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!
      w2 = 2.0_dp*w

        zen =        en  -gimp0
       del1 =       DelL +fimp0 
       del2 = conjg(DelL)+timp0 
        lam = sqrt(del1*del2-zen*zen)

        ga0 = -del1/(zen+w*lam)
        gb0 =  del2/(zen+w*lam)
!
! -- reservoir coherence functions
!      
      gaL = ga0
      gbL = gb0
      gaR = ga0
      gbR = gb0
!
!-----------------  gamma1 (p.v > 0) ----------------------------------
!
      it=1

      gold  = gaL
      gr_1(it)=gold
!
!-----------------  gamma2 (p.v > 0) ------------------------------------
!
      gold  = gbR
      gr_2(it)=gold
!
!-----------------  gamma1 (p.v < 0) ------------------------------------
!
      it=2

      gold = gaR
      gr_1(it)=gold
!
!-----------------  gamma2 (p.v <0 )------------------------------------
!
      gold = gbL
      gr_2(it)=gold
!
!--------------------------------------------------------------------
!
   end subroutine get_gammas
!
!---------------------------------------------------------------------
!
   subroutine get_imp_at_e(energy, iep, icmp, iprint)
!
!--------------------------------------------------------------------
!
      integer, intent(in) :: iep, iprint
      real(dp), intent(in) :: energy
      integer, intent(inout) :: icmp

      integer :: i, icmptyp,ic
      real(dp) :: err, rep, aip
      complex(dp) ::s11, s12, s21, s22  
!
! -- Get the right leading-order self energies
!
      icmptyp = 2
      icmp = 0
!
      en = energy+w*etaR
      call get_rimps(iep, err, icmp, icmptyp)

      s11 =    en-gimp0-uimp0
      s22 =    en-gimp0+uimp0
      s12 =       DelL +fimp0
      s21 = conjg(DelL)+timp0

      gr_1p(iep) = gr_1(1)
      gr_2p(iep) = gr_2(1)

      gr_1m(iep) = gr_1(2)
      gr_2m(iep) = gr_2(2)
!
!-- v.p > 0
!
      alpha1p(iep) = s11 + s21*gr_1p(iep)
       beta1p(iep) = s22 + s21*gr_1p(iep)

      alpha2p(iep) = s22 - s12*gr_2p(iep)
       beta2p(iep) = s11 - s12*gr_2p(iep)
!
!-- v.p > 0
!
      alpha1m(iep) = s11 + s21*gr_1m(iep)
       beta1m(iep) = s22 + s21*gr_1m(iep)

      alpha2m(iep) = s22 - s12*gr_2m(iep)
       beta2m(iep) = s11 - s12*gr_2m(iep)
!
!--------------------------------------------------------------------
!
   end subroutine get_imp_at_e
!
!--------------------------------------------------------------------
!
end module qc_farm
!
!--------------------------------------------------------------------
