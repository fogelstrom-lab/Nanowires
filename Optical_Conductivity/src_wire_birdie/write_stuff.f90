!-------------------------------------------------------------------
! Modernized module with procedures from qc.f and riccati.f
!-------------------------------------------------------------------
module write_stuff
  use qc_data
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: write_results, write_cuts
contains
!
!---------------------------------------------------------------------
!
   subroutine write_results
!
!--------------------------------------------------------------------
!
     integer :: i, ii, ixstep
     real(dp) :: E0,R

     ixstep= 1 !sx/40

     open(80,file='av_ReJ',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(avj(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImJ',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,aimag(avj(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ReS0',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(uimp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImS0',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,aimag(uimp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ReS3',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(gimp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImS3',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,aimag(gimp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ReSu',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(fimp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImSu',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,aimag(fimp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ReSd',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(timp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImSd',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,aimag(timp(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ReG',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(avg(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImG',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,-aimag(avg(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ReF',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,real(avf(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

     open(80,file='av_ImF',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        do i=1,sx,ixstep
           R=float(i)
           write(80,1000) R,E0,aimag(avf(i,ii))
        enddo
        write(80,*)
     enddo
     close(80)

1000 format(50(1x,e15.7))
!
!--------------------------------------------------------------------
!
   end subroutine write_results
!
!--------------------------------------------------------------------
!
   subroutine write_cuts(itype)
!
!--------------------------------------------------------------------
!
      integer :: ii, ix, itype
      real(dp) :: E0, x

      ix=sx/2

      if(itype ==0) then
         open(80,file='DoS_L.dat',status='unknown')
         open(81,file='DoS_C.dat',status='unknown')
         open(82,file='DoS_R.dat',status='unknown')

         open(90,file='Gimp_x.dat',status='unknown')
         open(91,file='Dimp_x.dat',status='unknown')
         open(92,file='Timp_x.dat',status='unknown')
         open(93,file='Himp_x.dat',status='unknown')

         do ii=-iemax,iemax,1
            E0=er(ii)/delta0
            write(80,1000) E0,-aimag(avg( 1,ii))
            write(81,1000) E0,-aimag(avg(ix,ii))
            write(82,1000) E0,-aimag(avg(sx,ii))
         
            write(90,1000) E0,gimp( 1,ii)+uimp( 1,ii)
            write(91,1000) E0,fimp( 1,ii)
            write(92,1000) E0,timp( 1,ii)
            write(93,1000) E0,gimp( 1,ii)-uimp( 1,ii)
         enddo
         write(90,*)
         write(91,*)
         write(92,*)
         write(93,*)
         do ii=-iemax,iemax,1
            E0=er(ii)/delta0
         
            write(90,1000) E0,gimp(ix,ii)+uimp(ix,ii)
            write(91,1000) E0,fimp(ix,ii)
            write(92,1000) E0,timp(ix,ii)
            write(93,1000) E0,gimp(ix,ii)-uimp(ix,ii)
         enddo
         write(90,*)
         write(91,*)
         write(92,*)
         write(93,*)
         do ii=-iemax,iemax,1
            E0=er(ii)/delta0
         
            write(90,1000) E0,gimp(sx,ii)+uimp(sx,ii)
            write(91,1000) E0,fimp(sx,ii)
            write(92,1000) E0,timp(sx,ii)
            write(93,1000) E0,gimp(sx,ii)-uimp(sx,ii)
         enddo

         close(80)
         close(81)
         close(82)

         close(90)
         close(91)
         close(92)
         close(93)
      else
         if(itype == 1) x = tau(1)
         if(itype == 2) x = tau(2)
         if(itype == 3) x = Ngap*2.0_dp
         if(itype == 4) x = phase/pi

         do ii=-iemax,iemax,1
            E0=er(ii)/delta0

            write(80,1000) x,E0,-aimag(avg( 1,ii))
            write(81,1000) x,E0,-aimag(avg(ix,ii))
            write(82,1000) x,E0,-aimag(avg(sx,ii))

         enddo

         write(80,*)
         write(81,*)
         write(82,*)
      endif

1000 format(50(1x,e15.7))

!
!--------------------------------------------------------------------
!
   end subroutine write_cuts
!
!--------------------------------------------------------------------
!
end module write_stuff 
!
!--------------------------------------------------------------------
