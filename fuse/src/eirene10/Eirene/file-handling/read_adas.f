!pb  21.11.06: index error corrected in defintion of ap%dte

      subroutine read_adas (ir,reac,isw,iz1)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, isw, iz1
      character(len=*), intent(in) :: reac
      integer :: nz, nde, nte, iza, ize, io, lc, ind, ian, ien,
     .           ide, ite, iz
      character(132) :: zeile
      type(adas_data), pointer :: ap

      read (29,*,iostat=io) nz, nde, nte, iza, ize

      if (io .ne. 0) then
        write (iunout,*) ' ERROR READING FILE FROM ADAS DATABASE '
        write (iunout,*) ' DIRECTORY IS ',reac
        call exit_own(1)
      end if

      if ((iz1 < iza) .or. (iz1 > ize)) then
        write (iunout,*) ' ERROR READING FILE FROM ADAS DATABASE '
        write (iunout,*) ' REQUESTED Z1 IS NOT AVAILABLE '
        write (iunout,*) ' Z1, ZA, ZE ',IZ1, IZA, IZE
        call exit_own(1)
      end if
      
      allocate (ap)
      allocate (ap%dens(nde))
      allocate (ap%temp(nte))
      allocate (ap%dde(nde))
      allocate (ap%dte(nte))
      allocate (ap%fit(nte,nde))

      ap%ndens = nde
      ap%ntemp = nte

      read (29,*)

      lc = len_trim(reac)
      if (reac(lc:lc) == 'r') then
        read (29,*)
        read (29,*)
      end if

! read densities
      read (29,*) (ap%dens(ide),ide=1,nde)

! read temperatures
      read (29,*) (ap%temp(ite),ite=1,nte)
      
! find appropriate Z1-block

      do

! read line between data blocks
        read (29,'(A132)') zeile
        if (zeile(1:5) == '-----') then
          ind = index(zeile,'Z1')
          if (ind == 0) cycle
          ian = ind + scan(zeile(ind+1:),' ')
          ien = ian + scan(zeile(ian+1:),'/') - 1
          read (zeile(ian:ien),*) iz
          if (iz == iz1) exit
        end if  
 
      end do

      do ite = 1, nte
        read (29,*) (ap%fit(ite,ide), ide = 1,nde)
      end do

      close (29)

! set up differenz arrays for densities and temperatures
      
      do ide=1,nde-1
        ap%dde(ide) = 1._dp / (ap%dens(ide+1) - ap%dens(ide))
      end do
      
      do ite=1,nte-1
        ap%dte(ite) = 1._dp / (ap%temp(ite+1) - ap%temp(ite))
      end do


      select case (isw)

      case (0)
        IF (REACDAT(IR)%LPOT) THEN
          WRITE (IUNOUT,*) ' POTENTIAL ALREADY SPECIFIED FOR REACTION',
     .                       IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lpot = .true.
        allocate (reacdat(ir)%pot)
        nullify (reacdat(ir)%pot%line)
        nullify (reacdat(ir)%pot%poly)
        nullify (reacdat(ir)%pot%hyd)
        reacdat(ir)%pot%adas => ap
        reacdat(ir)%pot%ifit = 3

      case (1)
        IF (REACDAT(IR)%LCRS) THEN
          WRITE (IUNOUT,*) ' CROSS SECTION ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lcrs = .true.
        allocate (reacdat(ir)%crs)
        nullify (reacdat(ir)%crs%line)
        nullify (reacdat(ir)%crs%poly)
        nullify (reacdat(ir)%crs%hyd)
        reacdat(ir)%crs%adas => ap
        reacdat(ir)%crs%ifit = 3

      case (2:4)
        IF (REACDAT(IR)%LRTC) THEN
          WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtc = .true.
        allocate (reacdat(ir)%rtc)
        nullify (reacdat(ir)%rtc%line)
        nullify (reacdat(ir)%rtc%poly)
        nullify (reacdat(ir)%rtc%hyd)
        reacdat(ir)%rtc%adas => ap
        reacdat(ir)%rtc%ifit = 3

      case (5:7)
        IF (REACDAT(IR)%LRTCMW) THEN
          WRITE (IUNOUT,*) ' MOMEMTUM WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtcmw = .true.
        allocate (reacdat(ir)%rtcmw)
        nullify (reacdat(ir)%rtcmw%line)
        nullify (reacdat(ir)%rtcmw%poly)
        nullify (reacdat(ir)%rtcmw%hyd)
        reacdat(ir)%rtcmw%adas => ap
        reacdat(ir)%rtcmw%ifit = 3

      case (8:10)
        IF (REACDAT(IR)%LRTCEW) THEN
          WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtcew = .true.
        allocate (reacdat(ir)%rtcew)
        nullify (reacdat(ir)%rtcew%line)
        nullify (reacdat(ir)%rtcew%poly)
        nullify (reacdat(ir)%rtcew%hyd)
        reacdat(ir)%rtcew%adas => ap
        reacdat(ir)%rtcew%ifit = 3

      case (11:12)
        IF (REACDAT(IR)%LOTH) THEN
          WRITE (IUNOUT,*) ' OTHER POLYNOMIAL FIT COEFFICIENTS',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%loth = .true.
        allocate (reacdat(ir)%oth)
        nullify (reacdat(ir)%oth%line)
        nullify (reacdat(ir)%oth%poly)
        nullify (reacdat(ir)%oth%hyd)
        reacdat(ir)%oth%adas => ap
        reacdat(ir)%oth%ifit = 3

      case default
        WRITE (IUNOUT,*) ' WRONG REACTION TYPE SPCIFIED '
        WRITE (IUNOUT,*) ' REACTION NO. ', IR
        WRITE (IUNOUT,*) ' REACTION TYPE H.', ISW
        CALL EXIT_OWN(1)
      end select
      
      return
      end subroutine read_adas
      
