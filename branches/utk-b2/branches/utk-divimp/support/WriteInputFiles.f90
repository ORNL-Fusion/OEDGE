
module common_data

  integer,parameter :: inunit = 10
  integer,parameter :: outunit = 11

  integer,parameter :: maxlinelength = 1000
  integer,parameter :: maxlines = 2000

  integer,parameter :: max_parameter_blocks = 10

  integer,parameter :: stdin = 5
  integer,parameter :: stdout = 6
  integer,parameter :: errunit = 6


end module common_data



module parameter_block_data

  use common_data

  type parameter_block
     character*1000 :: doc
     logical :: option_tag_spec
     character*1000 :: option_tag
     character*1000 :: option_tag_off_data
     character*1000 :: option_tag_on_data
     character*1000 :: option_tag_switch_off_value
     character*1000 :: parameter_tag
     integer ::        num_parameter_values 
     character*1000 :: parameter_value_list
     character*1 :: parameter_sep 
  end type parameter_block


  type (parameter_block) :: parameter_data(max_parameter_blocks)

  integer :: parameter_cycles(max_parameter_blocks)
  integer :: count_array(max_parameter_blocks) 
  integer :: pb_count 

  !  logical :: cases_done =.false.

contains


  subroutine init_count_array
    implicit none

    count_array = 1

  end subroutine init_count_array


  subroutine read_parameter_blocks(ierr)
    implicit none
    integer :: ierr

    integer :: in

    do in = 1,pb_count
       call read_parameter_block(parameter_data(in),ierr)
       parameter_cycles(in) = parameter_data(in)%num_parameter_values
    end do

  end subroutine read_parameter_blocks



!  subroutine read_number_parameter_blocks(ierr)
!    implicit none
!    integer :: ierr
!
!    read(stdin,*,iostat=ierr) pb_count
!
!  end subroutine read_number_parameter_blocks


  subroutine read_parameter_block(pb,ierr)
    implicit none
    integer ierr

    type (parameter_block) :: pb

    logical :: readingpb,finishedpb
    character*1000 :: inputline
    integer :: len,ios,start

    !
    ! Assign default values
    !
    !
    pb%option_tag_spec = .false.
    pb%parameter_sep = ' '
    pb%doc='CHANGED = '

    readingpb = .false.
    finishedpb = .false.
    ierr = 0

    do while (.not.finishedpb)

       read(stdin,'(a1000)',iostat=ios) inputline

       if (ios.ne.0) then 
          ! covers eof exit at end of last pb - also possible errors
          finishedpb=.true.
          ierr = ios
          exit
       endif

       len = len_trim(inputline)

       write(0,*) 'READ_PARAMETER_BLOCK:',trim(inputline)

       select case (inputline(1:3))

       case ('PB:')
          if (readingpb) then
             finishedpb = .true.
             backspace (stdin)
          else
             readingpb = .true.
          endif

       case ('PD:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%doc=trim(inputline(3+start:))
!          read(inputline(4:len),'(a)') pb%doc

       case ('OT:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%option_tag=trim(inputline(3+start:))
          pb%option_tag_spec=.true.
!          read(inputline(4:len),'(a)') pb%option_tag

       case ('OF:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%option_tag_off_data=trim(inputline(3+start:))
!          read(inputline(4:len),'(a)') pb%option_tag_off_data

       case ('ON:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%option_tag_on_data=trim(inputline(3+start:))
!          read(inputline(4:len),'(a)') pb%option_tag_on_data

       case ('OS:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%option_tag_switch_off_value=trim(inputline(3+start:))

!          read(inputline(4:len),'(a)') pb%option_tag_switch_off_value

       case ('PT:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%parameter_tag=trim(inputline(3+start:))

!          read(inputline(4:len),'(a)') pb%parameter_tag

       case ('PN:')
          ! Numeric read doesn't need adjustment
          read(inputline(4:len),*) pb%num_parameter_values

       case ('PV:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%parameter_value_list=trim(inputline(3+start:))
          
!          read(inputline(4:len),'(a)') pb%parameter_value_list

       case ('PS:')
          start = next_non_space(trim(inputline(4:)),1)
          pb%parameter_sep=trim(inputline(3+start:))
!          read(inputline(4:len),'(a)') pb%parameter_sep

       case default

       end select

    end do

  end subroutine read_parameter_block


  subroutine get_parameter_value(parameter_index,parameter_value_data,parameter_sep,parameter_value)
    implicit none
    integer :: parameter_index
    character*(*) :: parameter_value_data,parameter_sep,parameter_value

    !
    ! This routine extracts the parameter value at the parater_index position of the parameter_value_data where the 
    ! parameter fields are separated by the parameter_sep character
    !
    !

    integer :: parameter_count,parameter_start,next_sep
    !integer,external :: next_non_space
    integer :: maxlen
    logical :: found

    found = .false.

    maxlen = len_trim(parameter_value_data)

    parameter_count = 0

    parameter_start = next_non_space(parameter_value_data,1)

    do while (.not.found)

       parameter_count = parameter_count + 1

       next_sep = (parameter_start-1) + index(parameter_value_data(parameter_start:),parameter_sep)

       if (next_sep.eq.(parameter_start-1)) then 
          ! No more separators found - assign pointer to end of string
          next_sep = maxlen+1
       endif

       if (parameter_count.eq.parameter_index) then 
          ! Parameter found
          parameter_value = parameter_value_data(parameter_start:next_sep-1)
          found = .true.
       elseif (next_sep.eq.maxlen+1) then 
          ! No more parameters to be read - assign value of last parameter read
          found = .true.
          parameter_value = parameter_value_data(parameter_start:next_sep-1)
       else
          parameter_start = next_non_space(parameter_value_data,next_sep+1)
       endif
    end do

!    write(0,'(a,i6,1x,,a,a1,a)') 'PARAM:',parameter_index,trim(parameter_value_data),':',trim(parameter_value)


  end subroutine get_parameter_value



  integer function next_non_space(string,start)
    implicit none
    character*(*) :: string
    integer,intent(in) :: start

    logical :: found
    integer :: maxlen

    maxlen = len_trim(string)

    if (maxlen.le.0) then 
       next_non_space = 0
       return
    endif

    next_non_space = start
       
    found = .false.

    do while (.not.found)

       if (next_non_space.gt.maxlen) then 
          found = .true.
          next_non_space = 0 
       elseif (string(next_non_space:next_non_space).ne.' ') then 
          found = .true.
       else
          next_non_space = next_non_space + 1
       endif

    end do

  end function next_non_space

  recursive subroutine update_count(in,cases_done)
    implicit none
    integer :: in
    logical :: cases_done

    if (count_array(in).eq.parameter_cycles(in)) then 
       if (in+1.gt.pb_count) then 
          cases_done = .true.
       else
          call update_count(in+1,cases_done)
          count_array(in) = 1
       endif
    else
       count_array(in) = count_array(in) + 1
    endif
    return

  end subroutine update_count



end module parameter_block_data


module case_data

  use common_data
  use error_handling
  use parameter_block_data

  private

  type case_type
     character*1000 :: data(maxlines)
     character*1000 :: name
     integer :: line_count
  end type case_type

  type (case_type) :: basecase,outcase

  public :: read_case_name,read_case,write_case,read_string,create_case,create_case_name,read_input_data

contains

  subroutine read_string(casename,ierr)
    character*(*) casename

    read(stdin,'(a)',iostat=ierr) casename

    if (ierr.ne.0) then 

       call errmsg('READ_STRING:','ERROR READING CASE NAME FROM INPUT FILE:'//casename(1:len_trim(casename)))

       return
    end if

  end subroutine read_string


  subroutine read_case(casename,ierr)
    use common_data

    character*(*) :: casename
    integer :: ierr

    integer :: len,line_count,ios

    basecase%name = trim(casename)

    ierr = 0

    open(inunit,file=trim(casename),form='formatted',iostat=ierr) 

    if (ierr.ne.0) then 
       call errmsg('READ_CASE','ERROR READING INPUT CASE:'//trim(casename))
       stop
    endif

    ios = 0
    line_count = 0 

    do while (ios.eq.0)

       line_count = line_count+1

       if (line_count.gt.maxlines) then 
          call errmsg('READ_CASE','ERROR: TOO MANY LINES IN INPUT CASE')
          stop
       endif

       read(inunit,'(a1000)',iostat=ios) basecase%data(line_count)
!       write(0,'(a,i6,a)') 'BASECASE:',line_count,trim(basecase%data(line_count))
       
    end do

    basecase%line_count = line_count-1


    close(inunit)

    return

  end subroutine read_case



  subroutine write_case(casename,ierr)

    character*(*) :: casename
    integer :: ierr


    integer :: len,line_count

    ierr = 0

    open(outunit,file=trim(casename),form='formatted',iostat=ierr) 

    if (ierr.ne.0) then 
       call errmsg('WRITE_CASE','ERROR OPENING CASE NAME FOR OUTPUT:'//casename(1:len))
       stop
    endif

    do line_count = 1,outcase%line_count
       len = len_trim(outcase%data(line_count))
       write(outunit,'(a)') outcase%data(line_count)(1:len)
    end do

    close(outunit)

  end subroutine write_case


  subroutine create_case(count_array,pb_count,outcasename,outcasebasename,ierr)
    implicit none
    character*(*) :: outcasename,outcasebasename
    integer :: pb_count,ierr
    integer :: count_array(pb_count),in
    !type (case_type) :: tmpcase
    
    character*1000 :: desc1,desc2
    character*1000 :: pb_desc
    character*1000 :: case_desc


    !
    ! Apply the parameter blocks to the base input file using the parameter indices specified
    ! by the count_array
    !
    case_desc = ' '

    call get_base_case(outcase)

    !
    ! Add the base case name to the title line (first line of the input file)
    !

!    write(0,*) 'LINE1:',trim(outcase%data(1))

    read(outcase%data(1),*) desc1,desc2

!    write(0,*) 'DESC1:',trim(desc1)
!    write(0,*) 'DESC2:',trim(desc2)


    write(outcase%data(1),'(a1,a,a1,1x,a1,a,a1)') "'",trim(desc1),"'","'",trim(desc2)//' '//trim(outcasebasename),"'"

   ! write(0,*) 'TITLE:',trim(outcase%data(1))

    !
    ! Apply the parameter blocks to the case
    !

    do in = 1,pb_count
       call apply_pb(outcase,parameter_data(in),count_array(in),pb_desc)
       case_desc = trim(case_desc)//' : '//trim(pb_desc)
    end do

    write(0,'(4a)') 'WRITING CASE: ',trim(outcasebasename),' ',trim(case_desc)



    call write_case(outcasename,ierr)

  end subroutine create_case


  subroutine get_base_case(tmpcase)
    type (case_type) :: tmpcase

    tmpcase = basecase

  end subroutine get_base_case

  subroutine create_case_name(outcasename,outcaseext,outcasebasename,final_outcasename,case_index)
    implicit none
    integer :: case_index
    character*(*) :: outcasename,outcaseext,final_outcasename,outcasebasename

    character*1 :: d1
    character*2 :: d2
    character*3 :: d3
    character*4 :: d4

!    write(0,*) 'CASE_INDEX:',case_index,int(log10(real(case_index)))

    select case (int(log10(real(case_index))))

    case (0)
       write(d1,'(i1)') case_index
       final_outcasename = trim(outcasename) // d1

    case(1)
       write(d2,'(i2)') case_index
       final_outcasename = trim(outcasename) // d2

    case(2)

       write(d3,'(i3)') case_index
       final_outcasename = trim(outcasename) // d3

    case(3)
       write(d4,'(i4)') case_index
       final_outcasename = trim(outcasename) // d4

    case default

       final_outcasename = trim(outcasename) // '-error'

    end select
    
    outcasebasename = trim(final_outcasename)

    final_outcasename = trim(final_outcasename) // trim(outcaseext)


  end subroutine create_case_name




  subroutine apply_pb(tmpcase,pb,parameter_index,pb_desc)
    implicit none

    type (case_type) :: tmpcase
    type (parameter_block) :: pb
    integer :: parameter_index
    integer :: in
    character*(*) :: pb_desc

    !
    ! Now - need to apply the parameter variation to the input base file 
    ! INCLUDING any option switch values as well as parameter values. 
    !

    character*1000 :: doc,desc1,desc2
    character*1000 :: parameter_value

    call get_parameter_value(parameter_index,pb%parameter_value_list,pb%parameter_sep,parameter_value)

    ! This code uses a string comparison to evaluate whether the parameter value is equal to the option off value

    !
    ! Loop through the input file -
    ! 1) Adjust the comment line at the beginning 
    ! 2) Adjust any specified option tag
    ! 3) Adjust any specified value input
    !

    do in = 1,tmpcase%line_count
       !
       ! Description
       !
       if (tmpcase%data(in)(2:5).eq.'+A02') then 
            read(tmpcase%data(in),*) desc1,desc2
            write(tmpcase%data(in),'(a1,a,a1,1x,a1,a,a1)') "'",trim(desc1),"'","'",trim(desc2)//' :'// trim(pb%doc) //'='//trim(parameter_value),"'"
            pb_desc=trim(pb%doc)//'='//trim(parameter_value)
       endif

       !
       ! Option tag
       !

       if (tmpcase%data(in)(2:5).eq.trim(pb%option_tag)) then 
          if (pb%option_tag_spec.and.(trim(parameter_value).ne.trim(pb%option_tag_switch_off_value))) then 
            read(tmpcase%data(in),*) desc1
            write(tmpcase%data(in),'(a1,a,a1,4x,a)') "'",trim(desc1),"'",trim(pb%option_tag_on_data)
         else
            read(tmpcase%data(in),*) desc1
            write(tmpcase%data(in),'(a1,a,a1,4x,a)') "'",trim(desc1),"'",trim(pb%option_tag_off_data)
         endif
      endif
      !
      ! Parameter value tag
      !

       if (tmpcase%data(in)(2:5).eq.trim(pb%parameter_tag)) then 
            read(tmpcase%data(in),*) desc1
            write(tmpcase%data(in),'(a1,a,a1,4x,a)') "'",trim(desc1),"'",trim(parameter_value)
         endif

     end do


  end subroutine apply_pb

  subroutine read_input_data(casename,outcasename,outcaseext,ierr)
    implicit none

    character*(*) :: casename, outcasename,outcaseext
    integer :: ierr


    logical :: reached_parameter_blocks
    character*1000 :: inputline
    integer :: ios,start

    reached_parameter_blocks = .false.

    !
    ! Set default case name extension 
    !
    outcaseext = '.d6i'


    do while (.not.reached_parameter_blocks)

       read(stdin,'(a1000)',iostat=ios) inputline

       write(0,*) 'READ_INPUT_FILE:',trim(inputline)

       if (ios.ne.0) then 
          ierr = ios
          call errmsg('READ_INPUT_DATA: READ ERROR',ios)
          stop
       endif

       select case (inputline(1:3))


       case ('CN:')
          start = next_non_space(inputline(4:),1)
          casename = trim(inputline(3+start:))
          

       case ('ON:')
          start = next_non_space(inputline(4:),1)
          outcasename = trim(inputline(3+start:))


       case ('OE:')
          start = next_non_space(inputline(4:),1)
          outcaseext = trim(inputline(3+start:))

       case ('PN:')

          read(inputline(4:),*) pb_count

       case ('PB:')
          reached_parameter_blocks = .true.
          backspace (stdin)

       case default

          call errmsg('READ_INPUT_DATA:','UNEXPECTED TAG '//trim(inputline))


       end select


    end do

    write(0,*) trim(casename)
    write(0,*) trim(outcasename)
    write(0,*) trim(outcaseext)

    !
    ! Load parameter block data
    !
    call read_parameter_blocks(ierr)

    !
    ! load base case data
    !
    call read_case(casename,ierr)


  end subroutine read_input_data

end module case_data





program writeinputfiles
use error_handling
use parameter_block_data
use case_data

implicit none
!
!     This program reads in a data file with a list of specifications
!     for creating a set of simulation runs. 
!
!     It reads in a base case name then stores that base case - all 
!     other input files start as this base case. 
!
!     It then modifies the contents of the base case varying each 
!     parameter in order from the list of those to be changed.
!
!     Each change is documented in the comment block at the beginning of the file. 
!     The form of the documentation is read from the input file. 
!
!     The code can optionally create a new series of input files for each base
!     parameter change or just number them sequentially. 
!
!  Input file structure:
!
!  'CN:'  <Case name to read>
!  'ON:'  <Base output name for case> 
!  'OE:'  <Base output name extension - defaults to .d6i>
!
!  'NPC:' <number of parameters to vary>
!
!  'PB:' <parameter variation data block>
!       'PD:' <comment for documentation - default 'CHANGE = <parameter value>' > 
!       'OT:' <option parameter tag> 
!       'OF:' <option off data>
!       'ON:' <option on data>
!       'OS:' <value for option switch off> 
!       'PT:' <parameter tag>
!       'PN:' <number of values for this parameter>
!       'PV:' <parameter value list> 
!       'PS:' <parameter separator - defaults to space>
!        
! ... repeat parameter blocks ...
!
!
      character*1000 :: casename,outcasename,outcaseext,final_outcasename,outcasebasename

      integer :: ierr  ! error parameter

      integer :: in
      integer :: case_index_start = 1
      integer :: case_count, case_index
      logical :: cases_done
      


      call read_input_data(casename,outcasename,outcaseext,ierr)

      cases_done = .false.
      
      call init_count_array

      case_count = 0

      do while (.not.cases_done)

         case_index = case_count + case_index_start

         call create_case_name(outcasename,outcaseext,outcasebasename,final_outcasename,case_index)

         call create_case(count_array,pb_count,final_outcasename,outcasebasename,ierr)

         in =1 
         call update_count(in,cases_done)
        
         case_count= case_count+1
    
      end do

    end program writeinputfiles


