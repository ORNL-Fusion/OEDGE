module filter_elms
  use error_handling
  use utilities
  implicit none


  save
  private

  real,allocatable :: elm_data(:,:)
  real :: elm_start_offset, elm_end_offset, elm_effect_start_offset,elm_effect_end_offset
  integer :: elm_cnt


  public :: load_elm_times, set_elm_criteria


contains


  subroutine load_elm_times(filename,ierr)
    implicit none
    character*(*) filename


    integer unit,ierr
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

    close(unit)


  end subroutine load_elm_times


  subroutine set_elm_criteria(elm_start_offset,elm_end_offset,elm_effect_start_offset,elm_effect_end_offset)
    implicit none
    real :: elm_start_offset, elm_end_offset, elm_effect_start_offset,elm_effect_end_offset


    elm_start = elm_start_offset
    elm_end   = elm_end_offset

    elm_effect_start = elm_effect_start_offset
    elm_effect_end   = elm_effect_end_offset

    ! make start offsets negative
    if (elm_start.gt.0.0) elm_start = -elm_start
    if (elm_effect_start.gt.0.0) elm_effect_start = -elm_effect_start


  end subroutine set_elm_criteria


  subroutine get_elm_time(time,elm_time_offset,elmref,inelm)
    implicit none
    real :: time,elm_time_offset
    integer :: elmref
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



    elm_time_offset = 0.0
    elmref = 0

    ! special cases for boundary regions
    if (time.lt.(elm_data(1,2)+elm_effect_start)) then 
       inelm = .false.
       elm_time_offset = time - elm_data(1,2)
       elmref = 0
       return
    elseif (time.lt.(elm_data(1,2)+elm_start)) then 
       inelm = .false.
       elm_time_offset = time - elm_data(1,2)
       elmref = 1
       return
    elseif (time.le.elm_data(1,2)) then 
       inelm = .true.
       elm_time_offset = time - elm_data(1,2)
       elmref = 1
       return
    elseif (time.gt.elm_data(elm_cnt,2)+elm_effect_end) then 
       inelm = .false.
       elm_time_offset = time - elm_data(elm_cnt,2)
       elmref = 0
       return
    elseif (time.gt.elm_data(elm_cnt,2)+elm_end) then 
       inelm = .false.
       elm_time_offset = time - elm_data(elm_cnt,2)
       elmref = elm_cnt
       return
    elseif (time.ge.elm_data(elm_cnt,2)) then 
       inelm = .true.
       elm_time_offset = time - elm_data(elm_cnt,2)
       elmref = elm_cnt
       return
    endif


    ! finds nearest higher value
    tpos = ipos(time,elm_data(:,2),elm_cnt)

    ! time is between the ELMs at tpos-1 and tpos

    deltat1 = time - elm_data(tpos-1,2)
    deltat2 = time - elm_data(tpos,2)

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
    elseif (deltat2.gt.elm_effect_start)
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

  end subroutine get_elm_time



end module filter_elms
