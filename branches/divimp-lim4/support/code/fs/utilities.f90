module utilities



contains



  subroutine find_free_unit_number(unit)
    implicit none
    integer :: unit
    !
    !     FIND_FREE_UNIT_NUMBER:
    !
    !     This routine scans through unit numbers looking for one that
    !     is not currently in use. This number is returned. This code
    !     is based on the assumption that any unit numbers returned will
    !     be used before this routine is called again asking for another 
    !     number - otherwise it will likely return the previous value.
    !
    integer :: test_unit
    logical :: unit_open

    test_unit = 10
    unit_open = .false.

    ! Check for unit number assignment.  
    Do While (Unit_open.or.test_unit.gt.98)
       test_unit=test_unit + 1
       Inquire (Unit = test_unit, Opened = Unit_open)
    End Do

    if (unit_open) then 
       write(0,'(a)') 'ERROR: No unit numbers available:',test_unit
       stop
    else

       unit = test_unit
    endif

    return
  end subroutine find_free_unit_number


  function upcase(string) result(upper)
    character(len=*), intent(in) :: string
    character(len=len(string)) :: upper
    integer :: j

    do j = 1,len(string)
       if(string(j:j) >= "a" .and. string(j:j) <= "z") then
          upper(j:j) = achar(iachar(string(j:j)) - 32)
       else
          upper(j:j) = string(j:j)
       end if
    end do

  end function upcase




end module utilities
