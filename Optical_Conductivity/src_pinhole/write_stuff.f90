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
     integer :: ii
     real(dp) :: E0

     open(80,file='av_ReJ',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(avj(ii))
     enddo
     close(80)

     open(80,file='av_ImJ',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,aimag(avj(ii))
     enddo
     close(80)

     open(80,file='av_ReS0',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(uimp(ii))
     enddo
     close(80)

     open(80,file='av_ImS0',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,aimag(uimp(ii))
     enddo
     close(80)

     open(80,file='av_ReS3',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(gimp(ii))
     enddo
     close(80)

     open(80,file='av_ImS3',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,aimag(gimp(ii))
     enddo
     close(80)

     open(80,file='av_ReSu',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(fimp(ii))
     enddo
     close(80)

     open(80,file='av_ImSu',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,aimag(fimp(ii))
     enddo
     close(80)

     open(80,file='av_ReSd',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(timp(ii))
     enddo
     close(80)

     open(80,file='av_ImSd',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,aimag(timp(ii))
     enddo
     close(80)

     open(80,file='av_ReG',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(avg(ii))
     enddo
     close(80)

     open(80,file='av_ImG',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,-aimag(avg(ii))
     enddo
     close(80)

     open(80,file='av_ReF',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,real(avf(ii))
     enddo
     close(80)

     open(80,file='av_ImF',status='unknown')
     do ii=-iemax,iemax,1
        E0=er(ii)/delta0
        write(80,1000) E0,aimag(avf(ii))
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
         open(80,file='DoS.dat',status='unknown')

         open(90,file='Gimp.dat',status='unknown')
         open(91,file='Dimp.dat',status='unknown')
         open(92,file='Timp.dat',status='unknown')
         open(93,file='Himp.dat',status='unknown')

         do ii=-iemax,iemax,1
            E0=er(ii)/delta0
            write(80,1000) E0,-aimag(avg(ii))
         
            write(90,1000) E0,gimp(ii)+uimp(ii)
            write(91,1000) E0,fimp(ii)
            write(92,1000) E0,timp(ii)
            write(93,1000) E0,gimp(ii)-uimp(ii)
         enddo
         close(80)

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

            write(80,1000) x,E0,-aimag(avg(ii))

         enddo

         write(80,*)
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