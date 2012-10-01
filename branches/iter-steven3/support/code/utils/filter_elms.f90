module filter_elms
  use error_handling
  use utilities
  implicit none

  ! turn get_elm_time into an interface so I can add functionality without breaking existing code
  interface get_elm_time

     module procedure get_elm_time_ext,get_elm_time_base

  end interface

  save
  private

  real,allocatable :: elm_data(:,:)
  real :: elm_start, elm_end, elm_effect_start,elm_effect_end
  integer :: elm_cnt


  public :: load_elm_times, set_elm_criteria,get_elm_time


contains


  subroutine load_elm_times(filename,ierr)
    implicit none
    character*(*) filename


    integer unit,ierr,in
    character*256 :: line


    ! open elm time file

    call find_free_unit_number(unit)

    open(unit,file=trim(filename),status='old',form='formatted',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('ELM time file not found',filename)
       return
    endif

    ! preread file to determine number of entries
    elm_cnt = 0

    do while (ierr.eq.0)

       elm_cnt=elm_cnt + 1

       read(unit,*,iostat=ierr) line

    end do


    ierr = 0 
    elm_cnt = elm_cnt -1

    write(0,*) 'Allocating space for ELM data: N_elms =',elm_cnt


    ! allocate storage for ELM time and magnitude data

    allocate(elm_data(elm_cnt,3),stat=ierr)

    if (ierr.ne.0) then 
       call errmsg('Problem allocating ELM data array',ierr)
       return
    endif

    ! initialize to 0.0
    elm_data = 0.0


    ! read in ELM time data
    ! column 1 - index , column 2 - elm time, column 3 - elm magnitude

    rewind(unit)

    do in = 1,elm_cnt
       read(unit,*,iostat=ierr) elm_data(in,1),elm_data(in,2),elm_data(in,3)
    end do


    !do in = 1,elm_cnt
    !   write(6,'(a,i8,10(1x,g18.8))')  'Elm data:',in, elm_data(in,1),elm_data(in,2),elm_data(in,3)
    !end do



    close(unit)


  end subroutine load_elm_times


  subroutine set_elm_criteria(elm_start_input,elm_end_input,elm_effect_start_input,elm_effect_end_input)
    implicit none
    real :: elm_start_input, elm_end_input, elm_effect_start_input,elm_effect_end_input


    elm_start = elm_start_input
    elm_end   = elm_end_input

    elm_effect_start = elm_effect_start_input
    elm_effect_end   = elm_effect_end_input

    ! make start inputs negative if they aren't already
    if (elm_start.gt.0.0) elm_start = -elm_start
    if (elm_effect_start.gt.0.0) elm_effect_start = -elm_effect_start


  end subroutine set_elm_criteria

  subroutine get_elm_time_ext(time,elm_time_offset,elmref,inelm,elmdata,iflag)
    implicit none

    real :: time,elm_time_offset
    integer :: elmref,iflag
    logical :: inelm
    real :: elmdata

    call get_elm_time_data(time,elm_time_offset,elmref,inelm,elmdata,iflag)
    

  end subroutine get_elm_time_ext


  subroutine get_elm_time_base(time,elm_time_offset,elmref,inelm)
    implicit none
    real :: time,elm_time_offset
    integer :: elmref
    logical :: inelm

    ! local

    real :: elmdata
    integer :: iflag
    iflag = 0

    call get_elm_time_data(time,elm_time_offset,elmref,inelm,elmdata,iflag)

  end subroutine get_elm_time_base


  subroutine get_elm_time_data(time,elm_time_offset,elmref,inelm,elmdata,iflag)

    implicit none

    real :: time,elm_time_offset,elmdata
    integer :: elmref,iflag
    logical :: inelm

    integer :: tpos
    real :: deltat1,deltat2

    ! two time ranges around an ELM are defined
    !
    ! 1) elm_start : elm_end ... in this time frame the point is associated with the ELM
    ! 2) elm_effect_start : elm_effect:end ... if not in band (1) then the point is NOT in the ELM however, if it 
    !    is in band (2) it may be affected by the ELM and so data may be excluded. elmref and time offset are still set
    !    but inelm is set to .false. 
    ! 3) Not associated with an ELM ... outside of the given time windows around an ELM ... inelm is false. 
    ! 
    ! the input value iflag controls the nature of the information that will be loaded into elmdata
    ! iflag = 0 ... stores fraction of time between consecutive ELMs
    !

    elm_time_offset = 0.0
    elmref = 0
    elmdata = 0.0

    ! special cases for boundary regions
    if (time.lt.(elm_data(1,2)+elm_effect_start)) then 
       inelm = .false.
       elm_time_offset = time - elm_data(1,2)
       elmref = 0
       if (iflag.eq.0) then 
          elmdata = -1.0
       endif
       return
    elseif (time.lt.(elm_data(1,2)+elm_start)) then 
       inelm = .false.
       elm_time_offset = time - elm_data(1,2)
       elmref = 1
       if (iflag.eq.0) then 
          elmdata = -1.0
       endif
       return
    elseif (time.le.elm_data(1,2)) then 
       inelm = .true.
       elm_time_offset = time - elm_data(1,2)
       elmref = 1
       if (iflag.eq.0) then 
          elmdata = -1.0
       endif
       return
    elseif (time.gt.elm_data(elm_cnt,2)+elm_effect_end) then 
       inelm = .false.
       elm_time_offset = time - elm_data(elm_cnt,2)
       elmref = 0
       if (iflag.eq.0) then 
          elmdata = -1.0
       endif
       return
    elseif (time.gt.elm_data(elm_cnt,2)+elm_end) then 
       inelm = .false.
       elm_time_offset = time - elm_data(elm_cnt,2)
       elmref = elm_cnt
       if (iflag.eq.0) then 
          elmdata = -1.0
       endif
       return
    elseif (time.ge.elm_data(elm_cnt,2)) then 
       inelm = .true.
       elm_time_offset = time - elm_data(elm_cnt,2)
       elmref = elm_cnt
       if (iflag.eq.0) then 
          elmdata = -1.0
       endif
       return
    endif


    ! finds nearest higher value
    tpos = ipos(time,elm_data(:,2),elm_cnt)

    ! time is between the ELMs at tpos-1 and tpos

    deltat1 = time - elm_data(tpos-1,2)
    deltat2 = time - elm_data(tpos,2)

    if (iflag.eq.0) then 
       ! fraction of time since last ELM
       elmdata = deltat1/(deltat1+abs(deltat2))
    endif

    !write(6,'(a,i10,10(1x,g18.8))') 'Data:',tpos,time,elm_data(tpos-1,1),elm_data(tpos-1,2),elm_data(tpos,2),deltat1,deltat2

    ! need to consider cases of closely spaced ELMs ... in this case associate the data point with the later ELM.
    ! decision tree - check for association with later ELM first

    if (deltat2.gt.elm_start) then 
       ! associate with later ELM
       inelm = .true.
       elmref = tpos
       elm_time_offset = deltat2
       return
    elseif (deltat1.lt.elm_end) then 
       ! associate with earlier ELM
       inelm = .true.
       elmref = tpos-1
       elm_time_offset = deltat1
       return
    elseif (deltat2.gt.elm_effect_start) then
       inelm = .false.
       elmref = tpos
       elm_time_offset = deltat2
       return
    elseif (deltat1.lt.elm_effect_end) then 
       ! associate with earlier ELM
       inelm = .false.
       elmref = tpos-1
       elm_time_offset = deltat1
       return
    else
       ! not associated with an ELM
       ! assign time to closest preceding ELM 
       inelm = .false.
       elmref = 0
       elm_time_offset = deltat1
       return
    endif



  end subroutine get_elm_time_data





end module filter_elms
