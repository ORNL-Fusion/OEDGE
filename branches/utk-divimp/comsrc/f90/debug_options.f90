module debug_options

  implicit none

  !
  ! The purpose of the trace options is to allow debugging of crashes that do not report
  ! the exact location of the crash. Whenever pr_trace is called and debugging is on a 
  ! message indicating the current position in the code will be printed. Hopefully this will
  ! enable faster debugging of these errors in future. 
  !
  ! trace_execution - logical controlling on/off of function
  ! trace_unit - unit number for debugging output
  !


  integer :: trace_unit = 0
  logical :: trace_execution =.false.

  interface pr_trace

     module procedure s1write_trace,s2write_trace,iwrite_trace,rwrite_trace

  end interface


  public init_trace,pr_trace,toggle_trace


contains

  subroutine toggle_trace 
     implicit none
     trace_execution = .not.trace_execution

     if (trace_execution) then 
        write(0,'(a,i8)') 'TRACE TURNED ON: OUTPUT UNIT =',trace_unit
     else
        write(0,'(a)') 'TRACE TURNED OFF:'
     endif

  end subroutine toggle_trace


  subroutine init_trace(unitid,trace_state)
    implicit none
    integer :: unitid
    logical :: trace_state

    trace_execution = trace_state
    trace_unit = unitid

    if (trace_state) then 
       write(trace_unit,*) 'TRACE INITIALIZED:',trace_unit,trace_execution
    endif

  end subroutine init_trace


  subroutine s1write_trace(string1)
    implicit none

    character*(*) string1

    if (trace_execution) then 

       write(trace_unit,'(a,a)') 'TRACE:',trim(string1)

    endif

  end subroutine s1write_trace

  subroutine s2write_trace(string1,string2)
    implicit none

    character*(*) string1,string2

    if (trace_execution) then 

       write(trace_unit,'(a,a,a,a)') 'TRACE:',trim(string1),' DATA:',trim(string2)

    endif

  end subroutine s2write_trace


  subroutine iwrite_trace(string1,data1)
    implicit none

    character*(*) :: string1
    integer :: data1

    if (trace_execution) then 

       write(trace_unit,'(a,a,a,i8)') 'TRACE:',trim(string1),' DATA:',data1

    endif

  end subroutine iwrite_trace

  subroutine rwrite_trace(string1,data1)
    implicit none

    character*(*) :: string1
    real :: data1

    if (trace_execution) then 

       write(trace_unit,'(a,a,a,g12.5)') 'TRACE:',trim(string1),' DATA:',data1

    endif

  end subroutine rwrite_trace



end module debug_options
