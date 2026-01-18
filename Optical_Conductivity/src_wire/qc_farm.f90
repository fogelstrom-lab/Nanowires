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
         if (icmp .ge. 50000) write(*,*) ' Fail (',myrank,')', ien, err, aimag(gimp0(0)), icmp
      end do

      call MPI_Allreduce(MPI_IN_PLACE, counter, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

!     if(myrank==0) write(*,*) '========'
!     if(myrank==0) write(*,*) ' summed up',counter,2*iemax+1,float(counter)/float(2*iemax+1)
!     if(myrank==0) write(*,*) '========'

!-- MPI synchronize arrays so all ranks have a full copy for subsequent rounds

!-- Averaged Green's function components and equilibrium current contribution (r,e)
      call sync_cols_gather_indexed_rpairs(avj, MPI_COMM_WORLD)
      call sync_cols_gather_indexed_rpairs(avg, MPI_COMM_WORLD)
      call sync_cols_gather_indexed_rpairs(avf, MPI_COMM_WORLD)
      call sync_cols_gather_indexed_rpairs(avt, MPI_COMM_WORLD)

!-- Impurity self energies (r,e)
      call sync_cols_gather_indexed_rpairs(gimp, MPI_COMM_WORLD)
      call sync_cols_gather_indexed_rpairs(fimp, MPI_COMM_WORLD)
      call sync_cols_gather_indexed_rpairs(timp, MPI_COMM_WORLD)
      call sync_cols_gather_indexed_rpairs(uimp, MPI_COMM_WORLD)
!-- 
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
 
      complex(dp), allocatable :: gimp1(:),fimp1(:),timp1(:),dh(:)
      integer :: i1, i2, ic, ioshift    
      real(dp) :: x, x1, x2, p1, p2, ex, em, ah 
 
      save

      icmp = 0
      ic = sx/2
      ioshift = ien

      if(.not.allocated(dh)) allocate(dh(1:sx))
      if(.not.allocated(gimp1)) allocate(gimp1(1:sx))
      if(.not.allocated(fimp1)) allocate(fimp1(1:sx))
      if(.not.allocated(timp1)) allocate(timp1(1:sx))

      if(icmptyp == 1) then

         gimp0 = gimp(:,ien)
         fimp0 = fimp(:,ien)
         timp0 = timp(:,ien)

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
            gimp0 = p1*gimp(:,i1)+p2*gimp(:,i2)
            fimp0 = p1*fimp(:,i1)+p2*fimp(:,i2)
            timp0 = p1*timp(:,i1)+p2*timp(:,i2)
            
!           if(myrank == 0) write(29,1000) ex, gimp0(ic), fimp0(ic), timp0(ic), icmp
         else
            if(ex > 0) i1 = iemax
            if(ex < 0) i1 =-iemax
            gimp0 = gimp(:,i1)
            fimp0 = fimp(:,i1)
            timp0 = timp(:,i1)
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
         ah=     abs(dot_product(conjg(gimp1),gimp1))
         err=    abs(dot_product(conjg(dh),dh))/(ah+error_tol)

         dh=fimp0-fimp1
         ah=     abs(dot_product(conjg(fimp1),fimp1))
         err=err+abs(dot_product(conjg(dh),dh))/(ah+error_tol)

         dh=timp0-timp1
         ah=     abs(dot_product(conjg(timp1),timp1))
         err=err+abs(dot_product(conjg(dh),dh))/(ah+error_tol)

         err=sqrt(err)

         icmp=icmp+1
      enddo
!
! -- Save the leading-order impurity self-energies
!
      if(icmptyp == 1) then
         gimp(:,ien) = gimp0
         fimp(:,ien) = fimp0
         timp(:,ien) = timp0
         uimp(:,ien) = uimp0
      endif
!
! -- Save the t-matrices for linear response
!
      if(icmptyp == 2) then
         if(abs(srate(1)).gt. 1e-14) then
             gimpP(:,ioshift) = gimpP(:,ioshift)/srate(1)
             fimpP(:,ioshift) = fimpP(:,ioshift)/srate(1)
             timpP(:,ioshift) = timpP(:,ioshift)/srate(1)
          endif

          if(abs(srate(2)).gt. 1e-14) then
             gimpM(:,ioshift) = gimpM(:,ioshift)/srate(2)
             fimpM(:,ioshift) = fimpM(:,ioshift)/srate(2)
             timpM(:,ioshift) = timpM(:,ioshift)/srate(2)
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

      do ii=1,sx,1
!        if(abs(grid(ii))-Ngap > 0) srateb=srate(3)

         it=1
         pg = gr_1(ii,it)*gr_2(ii,it)

         norm= -w/(cone+pg)
         g1  = 2.0_dp*gr_1(ii,it)*norm
         g2  = 2.0_dp*gr_2(ii,it)*norm
         dig =        (cone-pg)*norm
 
         g3 = dig
         cj = dig
 
         it=2
         pg = gr_1(ii,it)*gr_2(ii,it)

         norm=-w/(cone+pg)
         g1  = g1+2.0_dp*gr_1(ii,it)*norm
         g2  = g2+2.0_dp*gr_2(ii,it)*norm
         dig =           (cone-pg)*norm

         g1 = g1*half
         g2 = g2*half
         g3 = (g3+dig)*half

         cj= cj-dig
!
! -- If iterating for equilibrium self enegies at given energy-grid point save the averages
!
         if(icmptyp == 1) then
            avj(ii,ien) = cj
            avg(ii,ien) = g3
            avf(ii,ien) = g1
            avt(ii,ien) = g2
         endif

         g32 = g3*g3
         norm = g32+g1*g2
!
! -- Potential scattering
!
         Denom = cone/(cone-sigma(1)*(cone+norm))
         cj = srate(1)*Denom
         fimpP(ii,ioshift) = cj*g1
         timpP(ii,ioshift) =-cj*g2
         gimpP(ii,ioshift) = cj*g3

         uimp0(ii) = srate(0)*cj
!
! -- Magnetic scattering
!
         D1 = (cone-sigma(2)*(cone-norm))
         Denom = cone/(D1*D1-magd*g32)
         cj = srate(2)*D1*Denom
         pg = magn*Denom

         fimpM(ii,ioshift) = -cj*g1
         timpM(ii,ioshift) =  cj*g2
         gimpM(ii,ioshift) = -(cj-pg)*g3

      enddo

      gimp0(1:sx) = gimpP(1:sx,ioshift) + gimpM(1:sx,ioshift)
      fimp0(1:sx) = fimpP(1:sx,ioshift) + fimpM(1:sx,ioshift)
      timp0(1:sx) = timpP(1:sx,ioshift) + timpM(1:sx,ioshift)
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
      complex(dp), allocatable :: zen(:), del1(:), del2(:), lam(:)
      complex(dp), allocatable :: ga0(:), gb0(:), lam2w(:), eX(:)
      complex(dp) :: gold, w2, gC, g0, gaL, gaR, gbL, gbR
      save
!
!--------------------------------------------------------------------
!
      if(.not.allocated(  zen)) allocate(  zen(1:sx))
      if(.not.allocated( del1)) allocate( del1(1:sx))
      if(.not.allocated( del2)) allocate( del2(1:sx))
      if(.not.allocated(  lam)) allocate(  lam(1:sx))
      if(.not.allocated(lam2w)) allocate(lam2w(1:sx))
      if(.not.allocated(  ga0)) allocate(  ga0(1:sx))
      if(.not.allocated(  gb0)) allocate(  gb0(1:sx))
      if(.not.allocated(   eX)) allocate(   eX(1:sx))
!
!--------------------------------------------------------------------
!
      w2 = 2.0_dp*w
      dx  = abs(grid(1)-grid(2)) ! *2.0_dp/2.0_dp

        zen(:) =           en  -gimp0(:)
       del1(:) =       delx(:) +fimp0(:) 
       del2(:) = conjg(delx(:))+timp0(:) 
        lam(:) = sqrt(del1(:)*del2(:)-zen(:)*zen(:))
      lam2w(:) = w2*lam(:)
         ex(:) = exp(-lam(:)*dx)

        ga0(:) = -del1(:)/(zen(:)+w*lam(:))
        gb0(:) =  del2(:)/(zen(:)+w*lam(:))
!
! -- reservoir coherence functions
!      
      gaL = ga0(1)
      gbL = gb0(1)
      gaR = ga0(sx)
      gbR = gb0(sx)
!
!-----------------  gamma1 (p.v > 0) ----------------------------------
!
      it=1

      gold  = gaL
      g0    = ga0(1)
      gC    = (gold-g0)/(lam2w(1)+del2(1)*(gold-g0))
      gold  = g0+lam2w(1)*gC*eX(1)/(cone-del2(1)*gC*eX(1))
      gr_1(1,it)=gold

      do i=2,sx,1
         g0    = ga0(i-1)
         gC    =  (gold-g0)/(lam2w(i-1)+del2(i-1)*(gold-g0))
         gold  =  g0+lam2w(i-1)*gC*eX(i-1)/(cone-del2(i-1)*gC*eX(i-1))

         g0    = ga0(i)
         gC    =  (gold-g0)/(lam2w(i)+del2(i)*(gold-g0))
         gold  =  g0+lam2w(i)*gC*eX(i)/(cone-del2(i)*gC*eX(i))

         gr_1(i,it)=gold
      enddo
!
!-----------------  gamma2 (p.v > 0) ------------------------------------
!
      gold  = gbR
      g0    = gb0(sx)
      gC    = (gold-g0)/(lam2w(sx)-del1(sx)*(gold-g0))
      gold  = g0+lam2w(sx)*gC*eX(sx)/(cone+del1(sx)*gC*eX(sx))
      gr_2(sx,it)=gold

      do i=sx-1,1,-1
         g0    = gb0(i+1)
         gC    = (gold-g0)/(lam2w(i+1)-del1(i+1)*(gold-g0))
         gold  = g0+lam2w(i+1)*gC*eX(i+1)/(cone+del1(i+1)*gC*eX(i+1))

         g0    = gb0(i)
         gC    = (gold-g0)/(lam2w(i)-del1(i)*(gold-g0))
         gold  = g0+lam2w(i)*gC*eX(i)/(cone+del1(i)*gC*eX(i))
         gr_2(i,it)=gold

      enddo
!
!-----------------  gamma1 (p.v < 0) ------------------------------------
!
      it=2

      gold = gaR
      g0   = ga0(sx)
      gC   = (gold-g0)/(lam2w(sx)+del2(sx)*(gold-g0))
      gold = g0+lam2w(sx)*gC*eX(sx)/(cone-del2(sx)*gC*eX(sx))
      gr_1(sx,it)=gold

      do i=sx-1,1,-1
         g0   = ga0(i+1)
         gC   = (gold-g0)/(lam2w(i+1)+del2(i+1)*(gold-g0))
         gold = g0+lam2w(i+1)*gC*eX(i+1)/(cone-del2(i+1)*gC*eX(i+1))

         g0   = ga0(i)
         gC   = (gold-g0)/(lam2w(i)+del2(i)*(gold-g0))
         gold = g0+lam2w(i)*gC*eX(i)/(cone-del2(i)*gC*eX(i))
         gr_1(i,it)=gold
      enddo
!
!-----------------  gamma2 (p.v <0 )------------------------------------
!
      gold = gbL
      g0   = gb0(1)
      gC   = (gold-g0)/(lam2w(1)-del1(1)*(gold-g0))
      gold = g0+lam2w(1)*gC*eX(1)/(cone+del1(1)*gC*eX(1))
      gr_2(1,it)=gold

      do i=2,sx,1
         g0   = gb0(i-1)
         gC   = (gold-g0)/(lam2w(i-1)-del1(i-1)*(gold-g0))
         gold = g0+lam2w(i-1)*gC*eX(i-1)/(cone+del1(i-1)*gC*eX(i-1))

         g0   = gb0(i)
         gC   = (gold-g0)/(lam2w(i)-del1(i)*(gold-g0))
         gold = g0+lam2w(i)*gC*eX(i)/(cone+del1(i)*gC*eX(i))
         gr_2(i,it)=gold
      enddo
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
      complex(dp), allocatable ::s11(:), s12(:), s21(:), s22(:)  
!
! -- 
!
      if(.not.allocated(s11)) allocate(s11(1:sx))
      if(.not.allocated(s12)) allocate(s12(1:sx))
      if(.not.allocated(s21)) allocate(s21(1:sx))
      if(.not.allocated(s22)) allocate(s22(1:sx))
!
! -- Get the right leading-order self energies
!
      icmptyp = 2
      icmp = 0
!
      en = energy+w*etaR
      call get_rimps(iep, err, icmp, icmptyp)

      s11(:) =    en-gimp0(:)-uimp0(:)
      s22(:) =    en-gimp0(:)+uimp0(:)
      s12(:) =       delx(:) +fimp0(:)
      s21(:) = conjg(delx(:))+timp0(:)

      gr_1p(:,iep) = gr_1(:,1)
      gr_2p(:,iep) = gr_2(:,1)

      gr_1m(:,iep) = gr_1(:,2)
      gr_2m(:,iep) = gr_2(:,2)
!
!-- v.p > 0
!
      alpha1p(:,iep) = s11(:) + s21(:)*gr_1p(:,iep)
       beta1p(:,iep) = s22(:) + s21(:)*gr_1p(:,iep)

      alpha2p(:,iep) = s22(:) - s12(:)*gr_2p(:,iep)
       beta2p(:,iep) = s11(:) - s12(:)*gr_2p(:,iep)
!
!-- v.p > 0
!
      alpha1m(:,iep) = s11(:) + s21(:)*gr_1m(:,iep)
       beta1m(:,iep) = s22(:) + s21(:)*gr_1m(:,iep)

      alpha2m(:,iep) = s22(:) - s12(:)*gr_2m(:,iep)
       beta2m(:,iep) = s11(:) - s12(:)*gr_2m(:,iep)
!
!-- Check that Im-part > 0
!
      if(myrank == 0 .and. iprint == 0) then
         ic = sx/2
         do i=1,sx,1
            rep =  real(alpha1p(i,iep))
            aip = aimag(alpha1p(i,iep))
            write(500,1000) grid(i), real(s21(i))
            write(600,1000) grid(i),aip
            write(700,1000) grid(i), real(gr_1p(i,iep))
            write(800,1000) grid(i),aimag(gr_1p(i,iep))

            rep =  real(alpha2p(i,iep))
            aip = aimag(alpha2p(i,iep))
            write(501,1000) grid(i),aimag(s21(i))
            write(601,1000) grid(i),aip
            write(701,1000) grid(i), real(gr_2p(i,iep))
            write(801,1000) grid(i),aimag(gr_2p(i,iep))

            rep =  real(alpha1m(i,iep))
            aip = aimag(alpha1m(i,iep))
            write(502,1000) grid(i), real(s12(i))
            write(602,1000) grid(i),aip
            write(702,1000) grid(i), real(gr_1m(i,iep))
            write(802,1000) grid(i),aimag(gr_1m(i,iep))

            rep =  real(alpha2m(i,iep))
            aip = aimag(alpha2m(i,iep))
            write(503,1000) grid(i),aimag(s12(i))
            write(603,1000) grid(i),aip
            write(703,1000) grid(i), real(gr_2m(i,iep))
            write(803,1000) grid(i),aimag(gr_2m(i,iep))
         enddo

         do i = 1,sx,1
            write(500,*)
            write(501,*)
            write(502,*)
            write(503,*)
            write(600,*)
            write(601,*)
            write(602,*)
            write(603,*)
            write(700,*)
            write(701,*)
            write(702,*)
            write(703,*)
            write(800,*)
            write(801,*)
            write(802,*)
            write(803,*)
         enddo
      endif

  1000 format(30(1x,e12.6))
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
