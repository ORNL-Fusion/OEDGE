      subroutine prep_rtcs (ir, iflg, ji, je, p1, cf)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, iflg, ji, je
      real(dp), intent(in) :: p1
      real(dp), intent(out) :: cf(9)
      real(dp) :: dum
      type(poly_data), pointer :: rp

      select case (iflg)
      case (3)
        if (.not.reacdat(ir)%lrtc) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR RATE COEFFICIENT',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%rtc%poly

      case (4)
        if (.not.reacdat(ir)%lrtcmw) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' MOMENTUM-WEIGHTED RATE COEFFICIENT',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%rtcmw%poly

      case (5)
        if (.not.reacdat(ir)%lrtcew) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' ENERGY-WEIGHTED RATE COEFFICIENT',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%rtcew%poly

      case (6)
        if (.not.reacdat(ir)%loth) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' OTHER REACTION',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%oth%poly

      case default
        write (iunout,*) ' call to prep_rtcs with wrong flag, iflg = ',
     .                     iflg
        write (iunout,*) ' 3 <= iflg <= 6 assumed '
        call exit_own(1)
      end select

      call dbl_poly (rp%dblpol,p1,0._dp,dum,cf,ji,je,
     .               rp%rcmn, rp%rcmx, rp%fparm, rp%ifexmn, rp%ifexmx)

      return
      end subroutine prep_rtcs
