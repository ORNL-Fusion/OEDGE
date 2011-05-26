module castem_field_line_data

  use error_handling
  !use utilities
  use common_utilities

  ! ribbon grid input options
  use ribbon_grid_options

  implicit none

  private

  save


  public :: read_identifier_data, read_castem_intersection_data,read_ray_intersection_data,print_field_line_summary,&
       calculate_castem_limiter_surface,calculate_ray_limiter_surface,generate_grid,write_grid,&
       assign_grid_to_divimp,deallocate_castem_storage


  !  type node_struc
  !     integer :: fl
  !     integer :: in
  !     real*8  :: s,r
  !  end type node_struc

  type intersection_struc
     integer :: line_id
     integer :: intsect_id
     real*8  :: xi,yi,zi,lc,metric,angle
     integer :: int_type
     integer :: bump
     logical :: int_used 
     !     type (node_struc) :: next
     !     type (node_struc) :: last
  end type intersection_struc

  type field_line_struc
     integer :: line_id
     real*8 :: xs,ys,zs,dist
     integer :: int_up, int_down, int_tot,int_tan
     integer :: int_stored
     type (intersection_struc), allocatable :: int_data(:)
  end type field_line_struc


  type (field_line_struc),allocatable :: field_line(:)

  integer :: n_tangency, n_enter, n_leave, n_field_lines, tot_n_intsects
  real*8 :: max_field_line_len,min_field_line_len
  real*8 :: max_lc,min_lc
  real*8 :: min_dist, max_dist
  !real*8 :: xstart,ystart,zstart,xend,yend,zend


  character*256 :: date_castem,date_process

  logical :: header_has_been_loaded

  !  type (node_struc) :: init_node

  !  type (node_struc),allocatable :: node_list

  integer :: n_nodes,node_cnt
  real*8,allocatable :: surf_r(:),surf_s(:),surf_fl(:),surf_int(:),surf_sep(:)
  real*8,allocatable :: av_s(:),av_r(:),av_type(:),av_min_r(:),av_max_r(:),av_angle(:), av_metric(:)
  real*8,allocatable :: av_tan_r(:),av_tan_s(:),av_wall_r(:),av_wall_s(:),tan_ord_r(:),av_tan_ind(:)


  real*8 :: r_limiter_max,r_limiter_min,s_limiter_max,s_limiter_min,min_tan_sep

  integer :: av_group_cnt
  integer :: av_wall_cnt, av_tan_cnt

  integer :: intsec_cnt


  ! arrays to help with defining the limiter surface from RAY data
  integer :: nbumps
  integer, allocatable :: bump_inf(:,:)

  real*8,allocatable :: av_tan_bump(:), av_tan_fl(:)



  ! grid generation
  integer, allocatable :: nvp(:)
  real *8, allocatable :: rvp(:,:), zvp(:,:),rcen(:,:),zcen(:,:)
  integer, allocatable :: nknots(:)
  integer :: npoly, nrings, max_nrings, max_nknots,max_npoly
  integer,allocatable :: poly_ref(:,:)

  real*8,parameter :: NEW_VERTEX = 5.0
  real*8,parameter :: FIXED_VERTEX = 6.0
  real*8,parameter :: SURFACE_START = 1.0
  real*8,parameter :: SURFACE_END = -1.0   ! 2.0
  real*8,parameter :: TANGENCY = 4.0
  real*8,parameter :: WALL = 3.0
  real*8,parameter :: DELETE_POINT = 10.0

  integer,parameter :: outunit = 6

  ! RAY intersection point type codes
  integer,parameter :: RAY_TAN = 0
  integer,parameter :: RAY_ENTER=1
  integer,parameter :: RAY_EXIT= 2
  integer,parameter :: RAY_END = 5


  ! Wall definition
  integer :: n_wall_segments,max_wall_segments
  real*8,allocatable :: wall_segments(:,:)
  real*8 :: wall_start_r,wall_start_s

  integer :: nwall
  real*8,allocatable :: wall_r(:),wall_s(:)

  ! Intersection subset selection
  logical :: filter_intersections



  ! options

  integer :: opt_block_av
  integer :: grid_option


  logical,parameter :: debug=.true.


contains


  subroutine read_identifier_data(file,ierr)
    implicit none
    character*(*) :: file
    integer :: ierr

    integer :: iunit

    ! local variables
    integer :: id,inter_up,inter_down,inter_tot,inter_tan
    real*8 :: xinit,yinit,zinit,dist
    integer :: in

    !    init_node%fl = -1
    !    init_node%in = -1

    ! initialize min and max r coordinate for grid
    min_dist =  1e24
    max_dist = -1e24

    ! initialize min and max s coordinate for grid
    min_lc =  1e24
    max_lc = -1e24


    ! find free unit number and open the file
    call find_free_unit_number(iunit)

    open(iunit,file=trim(file),status='old',iostat=ierr)

    if (ierr.ne.0) return

    header_has_been_loaded = .false.

    call read_file_header(iunit,ierr)

    if (ierr.ne.0) return

    header_has_been_loaded = .true.

    if (n_field_lines.gt.0) then 
       if (allocated(field_line)) deallocate(field_line)
       allocate(field_line(n_field_lines),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:FIELD_LINE:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    tot_n_intsects = n_enter + n_leave + n_tangency
    intsec_cnt = 0


    write(outunit,'(a,5(1x,i8))') 'Intsects:',tot_n_intsects,n_enter,n_leave,n_tangency,n_field_lines
    write(0,'(a,5(1x,i8))')       'Intsects:',tot_n_intsects,n_enter,n_leave,n_tangency,n_field_lines


    ! field line array has been allocated ... now read in and allocate the intersection
    ! arrays

    do in = 1,n_field_lines

       read(iunit,*) id,xinit,yinit,zinit,inter_up,inter_down,inter_tot,dist

       !write(0,*) 'READ:=',in,id,dist
       ! extract starting location from first field line 
       !if (id.eq.1) then 
       !   xstart = xinit
       !   ystart = yinit
       !   zstart = zinit
       !elseif (id.eq.n_field_lines) then 
       !   xend = xinit
       !   yend = yinit
       !   zend = zinit
       !endif

       ! Assign 0 for number of tangency points - included in RAY but not CASTEM listings
       inter_tan = 0

       call assign_field_line_data(field_line(in),id,xinit,yinit,zinit,inter_up,inter_down,inter_tan,inter_tot,dist,ierr)

    end do

    !write(outunit,'(a,10(1x,g18.8))') 'END POINTS:',xstart,ystart,zstart,xend,yend,zend


  end subroutine read_identifier_data


  subroutine assign_field_line_data(fl,id,xinit,yinit,zinit,inter_up,inter_down,inter_tan,inter_tot,dist,ierr)
    implicit none
    type(field_line_struc) :: fl
    integer :: id, inter_up,inter_down,inter_tot,inter_tan,ierr
    real*8 :: xinit,yinit,zinit,dist


    min_dist = min(min_dist,dist)
    max_dist = max(max_dist,dist)

    fl%line_id = id
    fl%xs = xinit
    fl%ys = yinit
    fl%zs = zinit

    ! the following line assumes that the first line in the sheaf is the first read in ... if not this 
    ! will need to be a post processed calculation
    ! This arbitrarily uses the R/X coordinate as the starting reference 
    !fl%dist = xstart+ sqrt((xinit-xstart)**2 + (yinit-ystart)**2 + (zinit-zstart)**2)
    !
    ! dist has been added to the IDENTIFIER file so calculation is not required
    fl%dist = dist

    fl%int_up = inter_up
    fl%int_down = inter_down 
    fl%int_tot = inter_tot

    fl%int_stored = 0

    ! number of tangency points on a line is included in RAY but not CASTEM output
    fl%int_tan  = inter_tan


    if (inter_tot.ne.(inter_up+inter_down+inter_tan)) then 
       call errmsg('Error: Intersections do not add up:',id)
    end if

    if (allocated(fl%int_data)) then 
       deallocate (fl%int_data)
    endif

    ! Allocate space to hold the total number of recorded intersections for this line
    ! only allocate if the line HAS any intersections
    if (inter_tot.gt.0) then 
       allocate(fl%int_data(inter_tot),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:FL%INT_DATA:IERR =',ierr)
          stop
       endif
    endif

    write(outunit,'(a,3(1x,i8))') 'Assign FL:',id,inter_tot,ierr

  end subroutine assign_field_line_data



  subroutine read_file_header(iunit,ierr)
    implicit none
    integer :: iunit,ierr

    integer :: n_key_strings,in
    character*256 :: key_strings(10)
    character*512 :: line
    logical :: header_read


    ! set up strings to scan for in input block

    n_key_strings = 9

    key_strings(1) = 'DATA'
    key_strings(2) = 'number of field lines followed :'
    key_strings(3) = 'number of tangency points :'
    key_strings(4) = 'number of points entering a surface :'
    key_strings(5) = 'number of points leaving a surface :'
    key_strings(6) = 'Date CASTEM calculations :'
    key_strings(7) = 'Date post-processing data :'
    key_strings(8) = 'maximum connection length toward top :'
    key_strings(9) = 'maximum connection length toward bottom :'

    ! if header data has been previously loaded then just scan for DATA and position 
    ! in file to read the data
    ! Note: the DATA tag is followed by a line with column headers that needs to be removed

    ierr = 0
    header_read = .false.


    do while (ierr.eq.0.and.(.not.header_read))


       read(iunit,'(a512)',iostat=ierr) line

       !write(outunit,'(a)') trim(line)

       if (header_has_been_loaded) then 

          ! Scan for the data tag
          if (trim(key_strings(1)).eq.line(1:len_trim(key_strings(1)))) then
             ! read the extra line
             read(iunit,'(a512)',iostat=ierr) line
             header_read = .true.
          endif

       else

          ! scan for various header data and perform appropriate actions
          do in = 1,n_key_strings

             if (trim(key_strings(in)).eq.line(1:len_trim(key_strings(in)))) then

                select case (in)
                case (1) 
                   ! key_strings(1) = 'DATA'
                   ! DATA tag
                   ! read the extra line 
                   read(iunit,'(a512)',iostat=ierr) line
                   header_read = .true.

                case (2)
                   ! key_strings(2) = 'number of field lines followed :'
                   read(line(len_trim(key_strings(in))+1:),*) n_field_lines
                case (3)
                   ! key_strings(3) = 'number of tangency points :'
                   read(line(len_trim(key_strings(in))+1:),*) n_tangency

                case (4)
                   ! key_strings(4) = 'number of points entering a surface :'
                   read(line(len_trim(key_strings(in))+1:),*) n_enter

                case (5)
                   ! key_strings(5) = 'number of points leaving a surface :'
                   read(line(len_trim(key_strings(in))+1:),*) n_leave

                case (6)
                   ! key_strings(6) = 'Date CASTEM calculations :'
                   read(line(len_trim(key_strings(in))+1:),'(a20)') date_castem

                case (7)
                   ! key_strings(7) = 'Date post-processing data :'
                   read(line(len_trim(key_strings(in))+1:),'(a20)') date_process

                case (8)
                   ! key_strings(8) = 'maximum connection length toward top :'
                   read(line(len_trim(key_strings(in))+1:),'(a20)') max_field_line_len

                case (9)
                   ! key_strings(9) = 'maximum connection length toward bottom :'
                   read(line(len_trim(key_strings(in))+1:),'(a20)') min_field_line_len

                case default
                   call errmsg('Error reading file header: Unknown tag number',in)
                   header_read = .true.
                   return
                end select

             endif

          end do

       endif

    end do


  end subroutine read_file_header




  subroutine read_castem_intersection_data(file,ierr)
    implicit none
    character*(*) :: file
    integer :: ierr


    integer :: iunit
    character*512 :: line

    integer :: line_id,inter_id,itype,iend,bump_id
    real*8 :: xint,yint,zint,lc,metric,angle

    integer :: sect_cnt


    ! find free unit number and open the file
    call find_free_unit_number(iunit)

    open(iunit,file=trim(file),status='old',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('Error opening Intersection file:'//trim(file),ierr)
       return
    endif

    ! this call positions the file at the beginning of the data listing
    call read_file_header(iunit,ierr)

    if (ierr.ne.0) then 
       call errmsg('Error reading Intersection file header:',ierr)
       return
    endif

    ! load the intersection data into the allocated arrays

    sect_cnt = 0

    do while (iend.eq.0) 

       read(iunit,'(a512)',iostat=iend) line
       !write(outunit,'(a)') trim(line)

       if (iend.eq.0) then 

          read(line,*) line_id,inter_id,xint,yint,zint,lc,itype

          !
          ! Assign default values for data not available in CASTEM files
          !
          bump_id = 0
          metric  = 0.0
          angle = 0.0

          call assign_intsect_data(line_id,field_line(line_id),inter_id,xint,yint,zint,lc,itype,bump_id,metric,angle,ierr)
          sect_cnt = sect_cnt + 1

       endif

    end do

    ! set error condition if no valid intersection data found
    if (sect_cnt.eq.0) then
       ierr = -2
    endif

    if (sect_cnt.ne.(n_tangency+n_enter+n_leave)) then 
       call errmsg('Incorrect number of intersection points read:',sect_cnt)
    endif

    !write(outunit,'(a,i10)') 'Total intersections read:', sect_cnt


  end subroutine read_castem_intersection_data


  subroutine read_ray_intersection_data(file,ierr)
    implicit none
    character*(*) :: file
    integer :: ierr


    integer :: iunit
    character*512 :: line



    integer :: line_id,inter_id,itype,iend
    real*8 :: xint,yint,zint,lc

    integer :: sect_cnt

    logical :: finished = .false. 

    real :: ver


    ! initialize min and max r coordinate for grid
    min_dist =  1e24
    max_dist = -1e24

    ! initialize min and max s coordinate for grid
    min_lc =  1e24
    max_lc = -1e24


    ! find free unit number and open the file
    call find_free_unit_number(iunit)

    open(iunit,file=trim(file),status='old',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('Error opening Intersection file:'//trim(file),ierr)
       return
    endif

    ! The RAY field line data file is tagged and divided into sections. 
    ! Comments in the file are denoted by a * in the first character position
    ! Each section has a fixed format. 

    !write(0,*) '1:'

    !    {VERSION}
    !  2.0
    !
    !    {TRACE SUMMARY}
    !
    !
    !* number of traces
    !       300
    !*
    !*  INDEX      - field line trace number
    !*  RHO_REL    - cross-field coordinate for the ribbon grid, with the origin on the inner most surface
    !*  RHO_ABS    - distance from the separatrix mapped to the outer midplane, from the true magnetic equilibrium
    !*  ORIGIN     - index of origin point on each field line (only useful for the _full data file)
    !*  X,Y,Z      - location of the origin point in Cartesian coordinates
    !*  # TANGENT  - number of tangency points
    !*  # INTER_LO - number of intersection points for S < 0
    !*  # INTER_HI - number of intersection points for S > 0
    !
    !*  CODE:
    !*
    !*    0 - tangency point
    !*    1 - wall intersection point, where the field line goes from inside the vacuum vessel volume to outside
    !*    2 - wall intersection point, where the field line goes from outside the vacuum vessel volume to inside
    !*    3 - point is inside the vacuum vessel volume
    !*    4 - point is outside
    !
    !*
    !*   index     rho_rel     rho_abs  origin    code           x           y           z     # tangent  # inter_lo  # inter_hi
    !*                 (m)         (m)                                   (m)         (m)         (m)
    !        1   0.0000000  -1.0000000     340       3   8.3000000   0.6000000   0.0000000             0           0           0
    !
    !
    !    {TRACE DATA}
    !
    !*  S      - distance along the field the trace, with S=0 at the origin point
    !*  IMPACT - impact angle at the surface for each intersection point
    !*  WIDTH  - for radial grid morphing
    !*  METRIC - for radial grid morphing
    !*
    !*   trace    npts
    !        1       2                                                                      (       1     930)
    !*   index    code             s             x           y           z        impact         width      metric          dphi
    !*                           (m)           (m)         (m)         (m)     (degrees)                               (degrees)
    !        1       5   -21.6504992    -5.9383600  -3.3095700  -1.0995400    -1.0000000    -1.0000000  -1.0000000  -169.510  -0.500       930  0
    !      930       5    33.4137013     2.1597300   4.6354500  -4.7369100    -1.0000000    -1.0000000  -1.0000000   -65.490   0.000         1  0
    !*   trace    npts
    !        2       2                                                                      (     931    1860)
    !*   index    code             s             x           y           z        impact         width      metric          dphi
    !*                           (m)           (m)         (m)         (m)     (degrees)                               (degrees)
    !        1       5   -21.6518355    -5.9388400  -3.3108800  -1.0996300    -1.0000000    -1.0000000  -1.0000000  -169.510  -0.500      1860  0
    !      930       5    33.4186278     2.1604400   4.6389600  -4.7384700    -1.0000000    -1.0000000  -1.0000000   -65.490   0.000       931  0
    !*   trace    npts
    !        3       3                                                                      (    1861    2790)
    !*   index    code             s             x           y           z        impact         width      metric          dphi
    !*                           (m)           (m)         (m)         (m)     (degrees)                               (degrees)
    !        1       5   -21.6004287    -5.9330700  -3.3081800  -1.1522000    -1.0000000    -1.0000000  -1.0000000  -169.010  -0.500      2790  0
    !      910       0    32.5711013     1.3734468   4.6236794  -5.0580945    -1.0000000    -1.0000000  -1.0000000   -74.799  -0.309      1881  3
    !      930       5    33.4216754     2.1608600   4.6413200  -4.7393800    -1.0000000    -1.0000000  -1.0000000   -65.490   0.000      1861  0
    !
    !
    ierr = 0


    do while (.not.finished) 

       read(iunit,'(a)',iostat=iend) line

       if (iend.eq.0) then 

          if (line(1:9).eq.'{VERSION}') then 

             !write(0,*) '2:'
             call read_ray_version(iunit,ver,ierr)
             if (ierr.ne.0) then 
                call errmsg('READ_RAY_VERSION: ERROR READING VERSION NUMBER:',ierr)
                stop 'ERROR:read_ray_version'
             endif

          elseif  (line(1:15).eq.'{TRACE SUMMARY}') then 

             !write(0,*) '3:'
             call read_ray_trace_summary(iunit,ierr)

          elseif (line(1:12).eq.'{TRACE DATA}') then 

             !write(0,*) '4:'
             call read_ray_trace_data(iunit,ierr)

          endif

       else
          finished = .true.
       endif


    end do


  end subroutine read_ray_intersection_data

  subroutine read_ray_version(iunit,ver,ierr)
    implicit none
    integer :: iunit,ierr
    real :: ver
    character*512 :: line

    logical :: done

    done = .false. 

    !write(0,*) '1a:'

    do while (.not.done) 
       read (iunit,'(a)',iostat=ierr) line

       if (line(1:1) .ne.'*') then 
          read(line,*) ver
          done = .true.
       endif

    end do

    write(outunit,'(a,g12.5)') 'RAY FILE VERSION #:',ver


  end subroutine read_ray_version


  subroutine read_ray_trace_summary(iunit,ierr)
    implicit none
    integer :: iunit,ierr

    character*512 :: line

    ! The trace summary header line has been read - now read in the number of traces followed by 
    ! the trace data lines. 

    logical :: first,finished

    ! variables to read trace summary data
    integer :: id, origin, code, inter_up,inter_down,inter_tan,inter_tot
    real*8 :: rho_rel, rho_abs,xinit,yinit,zinit,metric

    integer :: field_line_cnt

    first = .true.
    finished = .false.

    field_line_cnt = 0

    do while (.not.finished)

       read (iunit,'(a)') line

       !write(0,'(a)') ':'//trim(line)//':'

       ! All data lines start with a blank
       if (line(1:1).eq.' ') then 
          ! the first blank data line is expected to contain the number of field line traces
          if (first) then 
             read(line,*) n_field_lines

             !write(0,*) 'N_FIELD_LINES:',n_field_lines
             if (n_field_lines.gt.0) then 
                if (allocated(field_line)) deallocate(field_line)
                allocate(field_line(n_field_lines),stat=ierr)
                if (ierr.ne.0) then 
                   call errmsg('ALLOCATION ERROR:FIELD_LINE:IERR =',ierr)
                   stop
                endif

                first = .false. 
             else
                ierr=1
                call errmsg('READ_RAY_TRACE_SUMMARY','N_FIELD_LINES=0')
             endif

          else
             ! read in field line summary lines:
             !*
             !*   index     rho_rel     rho_abs  origin    code           x           y           z     # tangent  # inter_lo  # inter_hi
             !*                 (m)         (m)                                   (m)         (m)         (m)


             ! Notes: rho_rel contains coordinate that will be used for grid generation in this case. 
             !        May want to add some of the other data to the field_line_struc later

             read(line,*)  id,rho_rel,rho_abs,origin,code,xinit,yinit,zinit,inter_tan,inter_down,inter_up 
             !write(0,'(a)') 'FL:'//trim(line)//':'

             ! The field line listing in RAY files contains the end points but these are not totalled in the
             ! intersection summary numbers ... so 1 has to be added to both the n_int_up and n_int_down numbers to get the correct
             ! number for reading the intersection data (note: the intersection data contains a comment with the number of intersections
             ! and this data will be checked
             inter_up = inter_up+1
             inter_down = inter_down+1
             inter_tot = inter_down+inter_up+inter_tan

             call assign_field_line_data(field_line(id),id,xinit,yinit,zinit,inter_up,inter_down,inter_tan,inter_tot,rho_rel,ierr)

             field_line_cnt = field_line_cnt + 1

          endif

       endif

       if (field_line_cnt.ge.n_field_lines.and.(.not.first)) then 
          ! Set exit condition since all field lines should have been read in
          finished = .true.
       endif

    end do

    return 

  end subroutine read_ray_trace_summary



  subroutine read_ray_trace_data(iunit,ierr)
    implicit none
    integer :: iunit, ierr
    ! This routine reads in and assigns the individual intersection data

    logical :: finished
    integer :: field_line_cnt
    character*512 :: line

    integer :: node_id,itype,id2,code2
    real*8 :: s,xint,yint,zint,angle_to_surf,width,metric,toroidal_angle,dphi
    integer :: bump_id

    integer :: line_id,int_cnt
    integer :: sect_cnt

    finished = .false.
    field_line_cnt = 0 
    ierr = 0

    n_tangency = 0
    tot_n_intsects = 0
    n_enter = 0
    n_leave = 0


    !{TRACE DATA}
    !*
    !*  S      - distance along the field the trace, with S=0 at the origin point
    !*  IMPACT - impact angle at the surface for each intersection point
    !*  WIDTH  - for radial grid morphing
    !*  METRIC - for radial grid morphing
    !*
    !*   trace    npts
    !        1       2                                                                      (       1     932)
    !*   index    code    bump             s             x           y           z        impact         width      metric          dphi
    !*                                                       (m)           (m)         (m)         (m)     (degrees)                               (degrees)
    !        1       5       0   -21.6504992    -5.9383600  -3.3095700  -1.0995400    -1.0000000    -1.0000000  -1.0000000  -169.510  -0.500       932  0
    !      932       5       0    33.4137013     2.1597300   4.6354500  -4.7369100    -1.0000000    -1.0000000  -1.0000000   -65.490   0.000         1  0
    !*   trace    npts

    ! After the {TRACE DATA} tag we expect a compilation of intersections and tangency points for all of the field line traces 
    ! These will be listed in the following format - a trace identifier with the trace index and the number of intersections followed by the 
    ! listing of the intersections. There will be one block for each field line. 

    do while (.not.finished)

       read(iunit,'(a)') line
       !write(0,'(a)') 'LINE:',trim(line),':'

       if (line(1:1).eq.' ') then 
          ! read trace index and intersection count
          read(line,*) line_id, int_cnt
          !write (0,'(a,2i6,a)') 'READING:',line_id,int_cnt,':'//trim(line)//':'

          sect_cnt = 0

          do while (sect_cnt.lt.int_cnt)


             read(iunit,'(a)') line
             !write(0,'(a)') 'LINE:',trim(line),':'

             if (line(1:1).eq.' ') then 

                sect_cnt = sect_cnt + 1

                ! S and lc are the same quantity
                read(line,*) node_id,itype,bump_id,s,xint,yint,zint,angle_to_surf,width,metric,toroidal_angle,dphi,id2,code2
                ! note: the numerical values assigned to itype are different for RAY and CASTEM output
                !       this should be consolidated for use here to an internal typing system
                !       Additional intersection data from RAY will be assigned later. 

                call assign_intsect_data(line_id,field_line(line_id),sect_cnt,xint,yint,zint,s,itype,bump_id,metric,angle_to_surf,ierr)

                select case (itype)
                case(RAY_TAN) ! tangency
                   n_tangency = n_tangency + 1
                case(RAY_ENTER) ! enter surface
                   n_enter = n_enter + 1
                case(RAY_EXIT) ! leave surface
                   n_leave = n_leave + 1
                case(RAY_END) ! end point (counts as enter surface)
                   n_enter = n_enter + 1
                case default
                   call errmsg('Unexpected intersection type code =',itype)
                end select

             endif

          end do

          field_line_cnt = field_line_cnt + 1

       endif

       if (field_line_cnt.ge.n_field_lines) then 
          ! when intersection data has been read for all of the field lines - set exit condition
          finished = .true.
       endif

       tot_n_intsects = n_tangency + n_enter + n_leave

    end do

    return
  end subroutine read_ray_trace_data





  subroutine assign_intsect_data(line_id,fl,inter_id,xint,yint,zint,lc,itype,bump_id,metric,angle_to_surf,ierr)
    implicit none
    integer :: line_id
    type(field_line_struc) :: fl
    integer :: inter_id,itype,ierr,bump_id
    real*8 :: xint,yint,zint,lc,metric,angle_to_surf
    integer :: in

    ! Find min and max for lc

    min_lc = min(min_lc,lc)
    max_lc = max(max_lc,lc)

    ! Check for space in array
    if (fl%int_stored.lt.fl%int_tot) then 
       fl%int_stored = fl%int_stored + 1

       ! Assign index variable just to keep the code neater
       in = fl%int_stored

       fl%int_data(in)%line_id = fl%line_id
       fl%int_data(in)%intsect_id = inter_id
       fl%int_data(in)%xi = xint
       fl%int_data(in)%yi = yint
       fl%int_data(in)%zi = zint
       fl%int_data(in)%lc = lc
       fl%int_data(in)%metric = metric
       fl%int_data(in)%angle = angle_to_surf
       fl%int_data(in)%int_type = itype
       fl%int_data(in)%int_used = .false. 
       fl%int_data(in)%bump = bump_id


       !       fl%int_data(in)%next = node_init
       !       fl%int_data(in)%last = node_init

       intsec_cnt = intsec_cnt + 1


       if (in.ne.inter_id) then 
          call errmsg('Intersection count not equal to intersection ID',fl%line_id)
       endif

    else
       ! error - too many intersections
       call errmsg(' Too many intersections found for line',fl%line_id)
    endif


  end subroutine assign_intsect_data


  subroutine print_field_line_summary
    implicit none
    integer :: outunit, in,ierr
    character*100 :: fname

    outunit = 6

    write(outunit,'(a,i10)') ' Summary of Field line data read:',n_field_lines
    write(outunit,'(a,4(1x,i10))') ' Total Intersections expected: ',n_tangency, n_enter, n_leave, n_tangency+n_enter+n_leave

    do in = 1, n_field_lines
       write(outunit,'(a,3(1x,i12))') 'FL:',field_line(in)%line_id,field_line(in)%int_tot,field_line(in)%int_stored
    end do

    !fname = 'lim.txt'
    ! find free unit number and open the file
    !call find_free_unit_number(outunit)
    !open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(outunit,*) 'Writing out limiter surface: ' ,node_cnt
    write(0,*) 'Writing out limiter surface: ' ,node_cnt
    do in = 1,node_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,surf_r(in),surf_s(in),surf_fl(in),surf_int(in)
    end do

    !close(outunit)

    !fname = 'lim-av.txt'
    ! find free unit number and open the file
    !call find_free_unit_number(outunit)
    !open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(outunit,*) 'Writing out averaged limiter surface: ' ,av_group_cnt
    write(0,*) 'Writing out averaged limiter surface: ' ,av_group_cnt

    do in = 1,av_group_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,av_r(in),av_s(in)
       !write(0,'(i8,5(1x,g18.8))') in,av_r(in),av_s(in)
    end do

    !close(outunit)


    !fname = 'tan-av.txt'

    ! find free unit number and open the file
    !call find_free_unit_number(outunit)

    !open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(0,*) 'Writing out tangency points: ' ,av_tan_cnt
    write(outunit,*) 'Writing out tangency points: ' ,av_tan_cnt

    do in = 1,av_tan_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,av_tan_r(in),av_tan_s(in)
       !write(0,'(i8,5(1x,g18.8))') in,av_tan_r(in),av_tan_s(in)
    end do

    !close(outunit)

    !fname = 'wall-av.txt'

    ! find free unit number and open the file
    !call find_free_unit_number(outunit)

    !open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(0,*) 'Writing out wall points: ' ,av_wall_cnt
    write(outunit,*) 'Writing out wall points: ' ,av_wall_cnt

    do in = 1,av_wall_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,av_wall_r(in),av_wall_s(in)
       !write(0,'(i8,5(1x,g18.8))') in,av_wall_r(in),av_wall_s(in)
    end do

    !close(outunit)

  end subroutine print_field_line_summary



  subroutine calculate_castem_limiter_surface
    !use common_utilities
    implicit none

    integer :: in,if,ierr
    !integer :: node_cnt
    real,allocatable :: surf_sep(:),av_group(:)
    !real,allocatable :: surf_ang(:)

    real :: min_sep
    real,parameter :: max_fact = 5.0
    integer,parameter :: max_block_size = 100
    integer :: it, ix
    real :: tmp_r,tmp_s,tmp_fl,tmp_int
    real :: test_sep,last_sep
    logical :: leave
    integer :: group_id
    integer :: av_cnt

    integer :: tmp_wall_cnt, tmp_tan_cnt,last_face
    integer :: double_tan,double_wall
    real :: dir, last_dir

    integer :: index_offset

    !real :: ,last_ang,test_ang
    !real,parameter :: PI=3.141592654
    !real,parameter :: ang_limit = PI/18.0
    !real,parameter :: raddeg=57.29577952 

    !
    ! This routine goes through the intersection data and orders them so that when connected they form a piece-wise continuous description of the 
    ! intersection surface. The basic assumption is that the intersection data will be sortable by S coordinate ... that for a given specified S there is 
    ! only one surface intersection - thus no "islands" in the mesh formed by the intersections. Due to lack of data or less common situations this assumption
    ! may not always be true so the code will adjust the intersections as required. 
    !
    ! In addition, this code will also determine how to connect intersections on different limiters ... these would usually be identified by a "tangency" point 
    ! derivative change in the surface definition forming the bottom side of the tile as opposed to the peaks at limiter tips. Three methods can be used. 
    ! 1) connect the surface intersections
    ! 2) add a surface intersection point at the farthest out location
    ! 3) remove all intersections from the surface until the farther in intersection. 
    !
    !
    ! A second function of this code is to collect, order and create a reference list for all tangency points. If the rho (R) separation of tangency points is less 
    ! than a minimum value these will be moved onto the same tangency line. This data is used in grid generation to guide the choice of rho (R) values for the sides
    ! of the grid polygons
    !
    ! D.Elder 
    !
    ! run through all intersections - sort node_list based on S values
    !
    ! Alternatively ... go through the intersection points and join each to its nearest neighbour ... this should work given the grid resolution ... then double check 
    ! by verifying the Lc coordinate monotonicity 
    !
    ! Start with an initial point at min_dist, min_lc and end at min_dist, max_lc
    ! This means the number of nodes in the list should be the total number of intersections + 2 ... the first and last nodes will have no cell reference. 
    ! Using just the closest will not work - especially for the big gaps ... need to go back to sorting by S value then check distances to see if there are any fix ups
    ! needed. 

    ! Initialize averaging from global variable

    if (rg_int_win_mins.ne.rg_int_win_maxs) then 
       filter_intersections = .true.
    else
       filter_intersections = .false.
    endif

    opt_block_av = rg_block_av

    n_nodes = tot_n_intsects

    write(outunit,*) 'n_nodes:',n_nodes,tot_n_intsects

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_s)) deallocate(surf_s)
       allocate(surf_s(n_nodes),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:SURF_S:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_r)) deallocate(surf_r)
       allocate(surf_r(n_nodes),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:SURF_R:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_fl)) deallocate(surf_fl)
       allocate(surf_fl(n_nodes),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:SURF_FL:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_int)) deallocate(surf_int)
       allocate(surf_int(n_nodes),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:SURF_INT:IERR =',ierr)
          stop
       endif
    else
       return
    endif


    ! Populate the node list with intersection data 

    node_cnt = 0

    do if = 1,n_field_lines
       do in = 1,field_line(if)%int_tot

          if (.not.filter_intersections.or.&
               (filter_intersections.and. &
               (field_line(if)%int_data(in)%lc.ge.rg_int_win_mins).and.&
               (field_line(if)%int_data(in)%lc.le.rg_int_win_maxs))) then 
             node_cnt = node_cnt + 1
             surf_fl(node_cnt)=if
             surf_int(node_cnt)=in
             surf_s(node_cnt)=field_line(if)%int_data(in)%lc
             !
             ! jdemod - change this to dist - it will be the same in the case
             !          of a radial grid but will be better for a grid with a diagonal
             !          starting line. 
             !        - what is really needed here is PSIn or equivalent - need to 
             !          come up with something for starting lines not perpendicular to 
             !          the field lines in the poloidal plane. 
             !
             !surf_r(node_cnt)= field_line(if)%xs
             !
             surf_r(node_cnt)= field_line(if)%dist

             !          node_used(node_cnt) =0

          endif

       end do
    end do


    write(outunit,*) 'node_cnt:',node_cnt

    call sort_arrays(0,node_cnt,surf_s,surf_r,surf_fl,surf_int)


    ! if block averaging is on
    ! calculate intersection data separations

    if (opt_block_av.eq.1) then 

       if (tot_n_intsects.gt.0) then 
          if (allocated(surf_sep)) deallocate(surf_sep)
          allocate(surf_sep(node_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:SURF_SEP:IERR =',ierr)
             stop
          endif
       else
          return
       endif

       if (tot_n_intsects.gt.0) then 
          if (allocated(av_group)) deallocate(av_group)
          allocate(av_group(node_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_GROUP:IERR =',ierr)
             stop
          endif
       else
          return
       endif

       av_group_cnt = 1

       !if (tot_n_intsects.gt.0) then 
       !   if (allocated(surf_ang)) deallocate(surf_ang)
       !   allocate(surf_ang(node_cnt),stat=ierr)
       !else
       !   return
       !endif


       do in = 1,node_cnt-1

          surf_sep(in) = sqrt((surf_s(in+1)-surf_s(in))**2 + (surf_r(in+1)-surf_r(in))**2)
          !surf_ang(in) = atan2c(surf_s(in+1)-surf_s(in),surf_r(in+1)-surf_r(in))

          !             surf_sep(in1,in2) = sqrt((surf_s(in1)-surf_s(in2))**2 + (surf_r(in1)-surf_r(in2))**2)

       end do

       surf_sep(node_cnt) = surf_sep(node_cnt-1)


       ! run through grouping and reordering

       min_sep = minval(surf_sep)

       write(outunit,*) 'Minimum intersection separation = ', min_sep

       do in = 1,node_cnt
          ! look through the list for points that are separated by larger distances

          av_group(in) = av_group_cnt


          if (in.eq.1) then 
             last_sep = 10.0 * min_sep
             !last_ang = surf_ang(in)
          else
             last_sep = surf_sep(in-1)

             !if (last_sep.gt.max_fact * min_sep) then 
             !   last_ang = surf_ang(in)
             !else
             !   last_ang = surf_ang(in-1)
             !endif

          endif

          !write (6,'(a,i9,2l8,10(1x,g18.6))') 'SEP:',in,surf_sep(in).gt.max_fact * last_sep,abs(surf_ang(in)-last_ang).gt.ang_limit,surf_sep(in),max_fact*last_sep,surf_s(in),surf_r(in),surf_ang(in)*raddeg,last_ang*raddeg,ang_limit*raddeg

          !write (outunit,'(a,i9,l8,10(1x,g18.6))') 'SEP:',in,surf_sep(in).gt.max_fact * last_sep,surf_sep(in),max_fact*last_sep,surf_s(in),surf_r(in),av_group(in)

          ! Test to see if next point could be out of series
          !if (surf_sep(in).gt.max_fact * last_sep.or.abs(surf_ang(in)-last_ang).gt.ang_limit) then 
          if (surf_sep(in).gt.max_fact * last_sep) then 
             it = 1
             leave = .false.

             ! Look for additional points that are in series
             do while (.not.leave)
                it = it+1
                if (it.ge.max_block_size.or.(it+in+1).ge.node_cnt) then 
                   leave = .true.
                   av_group_cnt = av_group_cnt + 1
                else

                   test_sep = sqrt((surf_s(in+it)-surf_s(in))**2 + (surf_r(in+it)-surf_r(in))**2)
                   ! found a closer point - swap it with the next one up

                   !test_ang = atan2c(surf_s(in+it)-surf_s(in),surf_r(in+it)-surf_r(in))

                   !                   write(6,'(a,3i8,l8,10(1x,g18.6))') 'test:',in,it,in+it,test_sep.lt.max_fact*last_sep,test_sep,max_fact*last_sep,surf_s(in+it),surf_r(in+it),test_ang*raddeg,last_ang*raddeg, abs(test_ang-last_ang)

                   !write(outunit,'(a,3i8,l8,10(1x,g18.6))') 'test:',in,it,in+it,test_sep.lt.max_fact*last_sep,test_sep,max_fact*last_sep,surf_s(in+it),surf_r(in+it)

                   ! test for additional in series points 
                   !if (test_sep.lt.max_fact*last_sep.and.abs(test_ang-last_ang).lt.ang_limit) then 
                   if (test_sep.lt.max_fact*last_sep) then 
                      tmp_r = surf_r(in+it)
                      tmp_s = surf_s(in+it)
                      tmp_fl = surf_fl(in+it)
                      tmp_int = surf_int(in+it)

                      do ix = it-1,1,-1
                         surf_r(in+ix+1)   = surf_r(in+ix)
                         surf_s(in+ix+1)   = surf_s(in+ix)
                         surf_fl(in+ix+1)  = surf_fl(in+ix)
                         surf_int(in+ix+1) = surf_int(in+ix)
                         write(outunit,'(a,2i8,5(1x,g18.6))') 'Moving:',in,ix,in+ix,surf_r(in+ix+1), surf_s(in+ix+1)
                      end do

                      surf_r(in+1)   = tmp_r
                      surf_s(in+1)   = tmp_s
                      surf_fl(in+1)  = tmp_fl
                      surf_int(in+1) = tmp_int

                      do ix = 0,it
                         surf_sep(in+ix) = sqrt((surf_s(in+ix+1)-surf_s(in+ix))**2 + (surf_r(in+ix+1)-surf_r(in+ix))**2)
                      end do
                      leave = .true.
                   end if
                end if


             end do

          end if

       end do





       !
       ! Print out grouped data
       ! 


       ! slmod begin
       do in = 2,node_cnt-1
          !
          !do in = 2,node_cnt
          ! slmod end
          test_sep = sqrt((surf_s(in+1)-surf_s(in))**2 + (surf_r(in+1)-surf_r(in))**2)
          !write(outunit,'(a,i8,l8,10(1x,g18.6))') 'Nodes:', in,surf_sep(in).gt.max_fact*surf_sep(in-1),surf_r(in),surf_s(in),surf_sep(in),test_sep,av_group(in)


       end do


       !
       ! Allocate temp arrays for averaging
       !

       !write(outunit,*) 'Avgroup:',av_group_cnt, av_group(node_cnt)

       if (filter_intersections) then 
          ! need to add 2 points at each end when using a subset
          av_group_cnt = av_group(node_cnt) +4
          index_offset = 2
       else
          av_group_cnt = av_group(node_cnt) +2
          index_offset = 1
       endif


       if (av_group_cnt.gt.0) then 
          if (allocated(av_s)) deallocate(av_s)
          allocate(av_s(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_S:IERR =',ierr)
             stop
          endif
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_r)) deallocate(av_r)
          allocate(av_r(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_R:IERR =',ierr)
             stop
          endif
       else
          return
       endif

       ! record max and min r values in each grouping
       if (av_group_cnt.gt.0) then 
          if (allocated(av_min_r)) deallocate(av_min_r)
          allocate(av_min_r(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_MIN_R:IERR =',ierr)
             stop
          endif
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_max_r)) deallocate(av_max_r)
          allocate(av_max_r(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_MAX_R:IERR =',ierr)
             stop
          endif
       else
          return
       endif

       av_min_r = max_dist + 1.0
       av_max_r = min_dist - 1.0


       group_id = 0

       do in = 1,node_cnt


          if (av_group(in).ne.group_id) then 
             ! finish up last group
             if (group_id.ne.0) then
                av_s(group_id+index_offset) = av_s(group_id+index_offset)/av_cnt
                av_r(group_id+index_offset) = av_r(group_id+index_offset)/av_cnt
             endif

             ! start next group
             group_id = av_group(in)
             av_cnt = 1
             av_s(group_id+index_offset) = av_s(group_id+index_offset) + surf_s(in)
             av_r(group_id+index_offset) = av_r(group_id+index_offset) + surf_r(in)

             av_min_r(group_id+index_offset) = min(av_min_r(group_id+index_offset),surf_r(in))
             av_max_r(group_id+index_offset) = max(av_max_r(group_id+index_offset),surf_r(in))

          else
             ! continue counting current group
             av_cnt = av_cnt + 1
             av_s(group_id+index_offset) = av_s(group_id+index_offset) + surf_s(in)
             av_r(group_id+index_offset) = av_r(group_id+index_offset) + surf_r(in)

             av_min_r(group_id+index_offset) = min(av_min_r(group_id+index_offset),surf_r(in))
             av_max_r(group_id+index_offset) = max(av_max_r(group_id+index_offset),surf_r(in))

          endif

          !          write(6,'(a,2i10,10(1x,g18.8))') 'Nodes2:',in,group_id,surf_s(in),surf_r(in),av_min_r(group_id+1),av_max_r(group_id+1)

       end do

       ! finish off last group

       av_s(group_id+index_offset) = av_s(group_id+index_offset)/av_cnt
       av_r(group_id+index_offset) = av_r(group_id+index_offset)/av_cnt


       ! add first and last points
       ! 
       !av_r(1) = min_dist
       !av_s(1) = min_lc
       !av_min_r(1) = min_dist
       !av_max_r(1) = min_dist

       !av_r(av_group_cnt) = min_dist
       !av_s(av_group_cnt) = max_lc
       !av_min_r(av_group_cnt) = min_dist
       !av_max_r(av_group_cnt) = min_dist

    else

       ! copy over all nodes without block averaging - in case they fix the bugs in the data

       if (filter_intersections) then 
          ! filtered data sets will have two extra points at each end
          av_group_cnt = node_cnt +4
          index_offset = 2
       else
          av_group_cnt = node_cnt +2
          index_offset = 1
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_s)) deallocate(av_s)
          allocate(av_s(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_S:IERR =',ierr)
             stop
          endif
          av_s = 0.0
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_r)) deallocate(av_r)
          allocate(av_r(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_R:IERR =',ierr)
             stop
          endif
          av_r = 0.0
       else
          return
       endif
       ! record max and min r values in each grouping
       if (av_group_cnt.gt.0) then 
          if (allocated(av_min_r)) deallocate(av_min_r)
          allocate(av_min_r(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_MIN_R:IERR =',ierr)
             stop
          endif
          av_min_r = 0.0
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_max_r)) deallocate(av_max_r)
          allocate(av_max_r(av_group_cnt),stat=ierr)
          if (ierr.ne.0) then 
             call errmsg('ALLOCATION ERROR:AV_MAX_R:IERR =',ierr)
             stop
          endif
          av_max_r = 0.0
       else
          return
       endif


       do in = 1,node_cnt
          av_s(in+index_offset) = surf_s(in)
          av_r(in+index_offset) = surf_r(in)
          av_min_r(in+index_offset) = surf_r(in)
          av_max_r(in+index_offset) = surf_r(in)
       end do

    endif



    if (filter_intersections) then 

       ! add end faces - defined in input file

       av_r(1) = rg_minr
       av_s(1) = rg_mins
       av_min_r(1) = min_dist
       av_max_r(1) = min_dist

       av_r(2) = rg_maxr
       av_s(2) = rg_mins
       av_min_r(2) = min_dist
       av_max_r(2) = min_dist


       av_r(av_group_cnt-1) = rg_maxr
       av_s(av_group_cnt-1) = rg_maxs
       av_min_r(av_group_cnt-1) = min_dist
       av_max_r(av_group_cnt-1) = min_dist

       av_r(av_group_cnt) = rg_minr
       av_s(av_group_cnt) = rg_maxs
       av_min_r(av_group_cnt) = min_dist
       av_max_r(av_group_cnt) = min_dist


    else

       ! add end points

       av_r(1) = min_dist
       av_s(1) = min_lc
       av_min_r(1) = min_dist
       av_max_r(1) = min_dist

       av_r(av_group_cnt) = min_dist
       av_s(av_group_cnt) = max_lc
       av_min_r(av_group_cnt) = min_dist
       av_max_r(av_group_cnt) = min_dist

    endif





    !write(outunit,'(a,10(1x,g18.8))') 'Average1:',min_dist,max_dist
    !do in = 1,av_group_cnt
    !   write(0,'(a,i10,10(1x,g18.8))') 'AV:',in,av_s(in),av_r(in),av_max_r(in),av_min_r(in)
    !end do

    ! go through the data and remove degeneracies in tangency and wall points

    dir = 1.0   
    last_face = 1


    in = 0
    do while (in.lt.av_group_cnt-1)
       !do in = 1,av_group_cnt-1
       in = in+1
       last_dir = dir
       !write(outunit,'(a,2i8)') 'Index:',in,av_group_cnt
       dir = av_r(in+1)-av_r(in)

       if (dir.eq.0.and.last_face.eq.1) then 
          ! fix wall degeneracy
          call insert_limiter_point(in,last_face)
          last_face = 2
       elseif (dir.eq.0.and.last_face.eq.2) then 
          ! fix tangency degeneracy
          call insert_limiter_point(in,last_face)
          last_face = 1
       elseif (last_dir.eq.0.and.last_face.eq.1) then 
          last_face = 2
       elseif (last_dir.eq.0.and.last_face.eq.2) then 
          last_face = 1
       elseif (dir.lt.0.and.last_dir.gt.0) then
          last_face = 2
       elseif (dir.gt.0.and.last_dir.lt.0) then
          last_face = 1
       elseif (dir.gt.0.and.last_dir.gt.0) then
          last_face = 1
       elseif (dir.lt.0.and.last_dir.lt.0) then
          last_face = 2
       endif

       !write(outunit,'(a,2i6,10(1x,g18.8))') 'FIND:',in,last_face,dir,last_dir,av_r(in),av_s(in)

    end do




    ! Now run through and identify IN, OUT and tangency points - tangency points only occur at the top of a limiter ... not sure if IN/OUT are required


    ! allocate av_type
    if (av_group_cnt.gt.0) then 
       if (allocated(av_type)) deallocate(av_type)
       allocate(av_type(av_group_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_TYPE:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    dir = 1.0   
    av_tan_cnt = 0
    av_wall_cnt = 0

    last_face = 1

    do in = 1,av_group_cnt-1
       last_dir = dir
       dir = av_r(in+1)-av_r(in)

       if (dir.eq.0.and.last_face.eq.1) then 
          av_type(in) = WALL
          av_wall_cnt = av_wall_cnt + 1
          last_face = 2
       elseif (dir.eq.0.and.last_face.eq.2) then 
          av_type(in) = TANGENCY
          av_tan_cnt = av_tan_cnt + 1
          last_face = 1
       elseif (last_dir.eq.0.and.last_face.eq.1) then 
          av_type(in) = WALL
          av_wall_cnt = av_wall_cnt + 1
          last_face = 2
       elseif (last_dir.eq.0.and.last_face.eq.2) then 
          av_type(in) = TANGENCY
          av_tan_cnt = av_tan_cnt + 1
          last_face = 1
       elseif (dir.lt.0.and.last_dir.gt.0) then
          av_type(in) = WALL
          av_wall_cnt = av_wall_cnt + 1
          last_face = 2
       elseif (dir.gt.0.and.last_dir.lt.0) then
          av_type(in) = TANGENCY
          av_tan_cnt = av_tan_cnt + 1
       elseif (dir.gt.0.and.last_dir.gt.0) then
          av_type(in) = SURFACE_START
          last_face = 1
       elseif (dir.lt.0.and.last_dir.lt.0) then
          av_type(in) = SURFACE_END
          last_face = 2
       endif

       !write(0,'(a,3i6,10(1x,g18.8))')       'DIR:',av_tan_cnt,in,last_face,dir,last_dir,av_type(in),av_r(in),av_s(in)
       !write(outunit,'(a,3i6,10(1x,g18.8))') 'DIR:',av_tan_cnt,in,last_face,dir,last_dir,av_type(in),av_r(in),av_s(in)

    end do


    av_type(av_group_cnt) = SURFACE_END

    !
    ! Tangency points
    !

    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_r)) deallocate(av_tan_r)
       allocate(av_tan_r(av_tan_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_TAN_R:IERR =',ierr)
          stop
       endif
       av_tan_r = 0.0
    else
       return
    endif



    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_s)) deallocate(av_tan_s)
       allocate(av_tan_s(av_tan_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_TAN_S:IERR =',ierr)
          stop
       endif
       av_tan_s = 0.0
    else
       return
    endif


    ! Need index to limiter surface vertex as well
    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_ind)) deallocate(av_tan_ind)
       allocate(av_tan_ind(av_tan_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_TAN_IND:IERR =',ierr)
          stop
       endif
       av_tan_ind = 0
    else
       return
    endif

    !
    ! Tangency point radial order
    !

    if (av_tan_cnt.gt.0) then 
       if (allocated(tan_ord_r)) deallocate(tan_ord_r)
       allocate(tan_ord_r(av_tan_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:TAN_ORD_R:IERR =',ierr)
          stop
       endif
       tan_ord_r = 0.0
    else
       return
    endif

    !
    ! wall points
    !

    if (av_wall_cnt.gt.0) then 
       if (allocated(av_wall_s)) deallocate(av_wall_s)
       allocate(av_wall_s(av_wall_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_WALL_S:IERR =',ierr)
          stop
       endif
       av_wall_s = 0.0
    else
       return
    endif

    if (av_wall_cnt.gt.0) then 
       if (allocated(av_wall_r)) deallocate(av_wall_r)
       allocate(av_wall_r(av_wall_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_WALL_R:IERR =',ierr)
          stop
       endif
       av_wall_r = 0.0
    else
       return
    endif


    tmp_tan_cnt = 0
    tmp_wall_cnt = 0

    ! Collect tangency and wall data - check to make sure there are no two consective points with the same R value - i.e. flat tangency point
    double_tan = 0
    double_wall = 0

    do in = 1,av_group_cnt

       if (av_type(in).eq.TANGENCY) then 
          tmp_tan_cnt = tmp_tan_cnt + 1

          if (tmp_tan_cnt.le.av_tan_cnt) then 
             ! record r,s and index into surface array for the tangency point
             av_tan_r(tmp_tan_cnt) = av_r(in)
             av_tan_s(tmp_tan_cnt) = av_s(in)
             av_tan_ind(tmp_tan_cnt) = in
             write(outunit,'(a,i6,10(1x,g18.8))') 'TAN :',tmp_tan_cnt,in,av_r(in),av_s(in)
             if (in.ne.av_group_cnt) then 
                if(av_type(in+1).eq.TANGENCY) then 
                   double_tan = double_tan+1
                endif
             endif
          else
             CALL ERRMSG('Too many tangency points found:',tmp_tan_cnt)
             stop
          endif

       endif

       if (av_type(in).eq.WALL) then 
          tmp_wall_cnt = tmp_wall_cnt + 1
          if (tmp_wall_cnt.le.av_wall_cnt) then 
             av_wall_r(tmp_wall_cnt) = av_r(in)
             av_wall_s(tmp_wall_cnt) = av_s(in)
             write(outunit,'(a,i6,10(1x,g18.8))') 'WALL:',tmp_wall_cnt,in,av_r(in),av_s(in)
             if (in.ne.av_group_cnt) then 
                if(av_type(in+1).eq.WALL) then 
                   double_wall = double_wall+1
                endif
             endif
          else
             CALL ERRMSG('Too many wall points found:',tmp_wall_cnt)
             stop
          endif
       endif

    end do


    write(outunit,*) 'DOUBLES:',double_tan,double_wall

    !write(outunit,'(a,6i10)') 'Tan:',av_tan_cnt, av_wall_cnt,double_tan,double_wall

    !do in = 1,av_wall_cnt
    !   write(0,'(a,i6,10(1x,g18.8))') 'Wall:',in,av_wall_r(in),av_wall_s(in)
    !end do

    !do in = 1,av_tan_cnt
    !   write(0,'(a,i6,10(1x,g18.8))') 'Tan :',in,av_tan_r(in),av_tan_s(in)
    !end do

    !write(outunit,'(a,10(1x,g18.8))') 'Average2:',min_dist,max_dist
    !do in = 1,av_group_cnt
    !   write(0,'(a,i10,10(1x,g18.8))') 'AV:',in,av_s(in),av_r(in),av_max_r(in),av_min_r(in)
    !end do


    ! Find max and min from the reduced intersection data ... loop through revised wall to get this data

    r_limiter_max = min_dist -1.0
    r_limiter_min = max_dist +1.0

    s_limiter_max = min_lc -1.0
    s_limiter_min = max_lc +1.0


    do in = 1,av_group_cnt

       if (av_r(in).lt.min_dist-1.0) then 
          write(0,*) 'PROB:1:',in,min_dist,av_r(in)
       elseif (av_r(in).gt.max_dist+1) then 
          write(0,*) 'PROB:2:',in,max_dist,av_r(in)
       endif

       r_limiter_max = max(r_limiter_max,av_r(in))
       r_limiter_min = min(r_limiter_min,av_r(in))

       s_limiter_max = max(s_limiter_max,av_s(in))
       s_limiter_min = min(s_limiter_min,av_s(in))

    end do

    write(outunit,'(a,10(1x,g18.8))') 'R,S RANGES:',r_limiter_min,r_limiter_max,s_limiter_min,s_limiter_max
    !write(0,'(a,10(1x,g18.8))') 'R,S RANGES:',r_limiter_min,r_limiter_max,s_limiter_min,s_limiter_max



    ! Ok - now need to check for tangency points that are very close together and adjust the wall so that these points lie on the same line. 
    !    - in this case the larger r value will be assigned the smaller ... the minimum separation will be quite small compared to the simulation space 
    !      so this may not have a large impact
    !
    ! What is the minimum radial tangency separation distance?
    ! 
    !    
    ! set separation limit to r_limiter_max - r_limiter_min / 100 ... or just choose an absolute value like 0.001m 
    !
    min_tan_sep = (r_limiter_max-r_limiter_min)/100.0



    ! sort tangency points by r order to move r points that are too close together
    call sort_arrays(0,av_tan_cnt,av_tan_r,av_tan_s,av_tan_ind)


    ! go through and adjust tangency points as needed. 
    ! since wall points are not required to be a polygon vertex this processing isn't required for wall extremities

    do in = 1,av_tan_cnt-1
       if ((av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep) then 
          ! move tangency point in to match
          av_tan_r(in+1) = av_tan_r(in)
          av_r(int(av_tan_ind(in+1))) = av_r(int(av_tan_ind(in)))

          write(outunit,'(a,2i10,l8,15(1x,g18.8))') 'Reduce tan:',in,in+1,(av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep,av_tan_r(in+1),av_tan_r(in),av_tan_r(in+1)-av_tan_r(in),min_tan_sep,av_tan_s(in),av_tan_s(in+1)
          !write(0,'(a,2i10,l8,15(1x,g18.8))') 'Reduce tan:',in,in+1,(av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep,av_tan_r(in+1),av_tan_r(in),av_tan_r(in+1)-av_tan_r(in),min_tan_sep,av_tan_s(in),av_tan_s(in+1)
       endif

    end do


    ! put tangency points back into s order - this should be the same order as before the previous re-order - this is needed for generating cells while the r 
    ! r order is useful for generating rows
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,tan_ord_r,av_tan_ind)

    ! record the S order of the points 

    do in = 1,av_tan_cnt
       tan_ord_r(in) = in
       write(outunit,'(a,i8,10(1x,g18.8))') 'TAN_ORD_R:',in,av_tan_r(in),tan_ord_r(in)
    end do

    ! Sort the S indices by R to get tan_ord_r referencing the R order of the data while sorted by S
    call sort_arrays(0,av_tan_cnt,av_tan_r,av_tan_s,tan_ord_r,av_tan_ind)

    ! Re-sort to S - leaving tan_ord_r unchanged
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,av_tan_ind)


    do in = 1,av_tan_cnt
       write(outunit,'(a,i8,10(1x,g18.8))') 'TANS:',in,av_tan_r(in),av_tan_s(in),av_tan_ind(in)
       !write(0,'(a,i8,10(1x,g18.8))') 'TANS:',in,av_tan_r(in),av_tan_s(in),av_tan_ind(in)
    end do

    do in = 1,av_wall_cnt
       write(outunit,'(a,i8,10(1x,g18.8))') 'WALL:',in,av_wall_r(in),av_wall_s(in)
       !write(0,'(a,i8,10(1x,g18.8))') 'WALL:',in,av_wall_r(in),av_wall_s(in)
    end do


  end subroutine calculate_castem_limiter_surface


  subroutine insert_limiter_point(num,last_face)
    implicit none
    integer :: num
    integer :: last_face
    real, allocatable :: tmp_r(:), tmp_s(:),tmp_min_r(:),tmp_max_r(:)
    integer :: tmp_cnt
    integer :: in,inmod,ierr

    tmp_cnt = av_group_cnt + 1

    ! allocate tmp_r, tmp_s and tmp_type
    if (tmp_cnt.gt.0) then 

       if (allocated(tmp_r)) deallocate(tmp_r)
       allocate(tmp_r(tmp_cnt),stat=ierr)

       if (allocated(tmp_s)) deallocate(tmp_s)
       allocate(tmp_s(tmp_cnt),stat=ierr)

       if (allocated(tmp_min_r)) deallocate(tmp_min_r)
       allocate(tmp_min_r(tmp_cnt),stat=ierr)

       if (allocated(tmp_max_r)) deallocate(tmp_max_r)
       allocate(tmp_max_r(tmp_cnt),stat=ierr)

       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:TMP_ ARRAYS:IERR =',ierr)
          stop
       endif
    else
       return
    endif


    !
    ! Copy data from av arrays 
    !

    !write(0,*) 'Nums:',tmp_cnt,av_group_cnt,last_face,num

    do in = 1, av_group_cnt
       if (in.le.num) then 
          inmod = 0
       else 
          inmod = 1
       endif
       tmp_r(in+inmod) = av_r(in)
       tmp_s(in+inmod) = av_s(in)
       tmp_min_r(in+inmod) = av_min_r(in)
       tmp_max_r(in+inmod) = av_max_r(in)
    end do

    ! determine S value by averaging
    tmp_s(num+1) = (tmp_s(num) + tmp_s(num+2))/2.0

    ! determine new R value ... move outward in R for wall point and inward in R for tangency
    ! assign to average of max or min of intersection range ... 

    if (last_face.eq.1) then 
       ! last limiter facing was heading toward wall ... therefore move new R outward slightly

       if (tmp_r(num).ne.tmp_r(num+2)) then 
          write(0,'(a,i8,10(1x,g18.8))') 'ERROR: non-degenerate wall point being added:',num,tmp_r(num),tmp_r(num+2)
       endif
       if (tmp_max_r(num).ne.tmp_r(num)) then 
          ! use max r to set new point 
          tmp_r(num+1) = (tmp_max_r(num)+tmp_max_r(num+2))/2.0
       else
          ! otherwise choose a very small displacement ... abs(max_dist-min_dist)/1000.0
          tmp_r(num+1) = tmp_r(num) + abs((max_dist-min_dist)/1000.0d0) 
       endif
       tmp_max_r(num+1) = tmp_r(num+1)
       tmp_min_r(num+1) = tmp_r(num+1)
       write(outunit,'(a,3(1x,g18.8),a,3(1x,g18.8))') 'Inserting single wall point:',tmp_r(num),tmp_r(num+1),tmp_r(num+2),':',tmp_s(num),tmp_s(num+1),tmp_s(num+2)

    elseif (last_face.eq.2) then 
       ! last limiter facing was heading toward tip ... therefore move new R inward slightly
       if (tmp_r(num).ne.tmp_r(num+2)) then 
          write(0,'(a,i8,10(1x,g18.8))') 'ERROR: non-degenerate wall point being added:',num,tmp_r(num),tmp_r(num+2)
       endif
       if (tmp_min_r(num).ne.tmp_r(num)) then 
          tmp_r(num+1) = (tmp_min_r(num)+tmp_min_r(num+2))/2.0
       else
          tmp_r(num+1) = tmp_r(num) - abs((max_dist-min_dist)/1000.0d0) 
       endif
       tmp_max_r(num+1) = tmp_r(num+1)
       tmp_min_r(num+1) = tmp_r(num+1)
       write(outunit,'(a,3(1x,g18.8),a,3(1x,g18.8))') 'Inserting single tan point:',tmp_r(num),tmp_r(num+1),tmp_r(num+2),':',tmp_s(num),tmp_s(num+1),tmp_s(num+2)
    else
       call errmsg('INSERT_LIMITER_POINT:Invalid value of last_face =',last_face)
    endif

    ! deallocate and reallocate the av arrays

    av_group_cnt = tmp_cnt

    if (av_group_cnt.gt.0) then 
       if (allocated(av_s)) deallocate(av_s)
       allocate(av_s(av_group_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('RE-ALLOCATION ERROR:AV_S:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (av_group_cnt.gt.0) then 
       if (allocated(av_r)) deallocate(av_r)
       allocate(av_r(av_group_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('RE-ALLOCATION ERROR:AV_R:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (av_group_cnt.gt.0) then 
       if (allocated(av_min_r)) deallocate(av_min_r)
       allocate(av_min_r(av_group_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('RE-ALLOCATION ERROR:AV_MIN_R:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (av_group_cnt.gt.0) then 
       if (allocated(av_max_r)) deallocate(av_max_r)
       allocate(av_max_r(av_group_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('RE-ALLOCATION ERROR:AV_MAX_R:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    ! copy tmp data back to new av arrays and deallocate tmp arrays

    av_s = tmp_s
    av_r = tmp_r
    av_max_r = tmp_max_r
    av_min_r = tmp_min_r

    if (allocated(tmp_s)) deallocate(tmp_s)
    if (allocated(tmp_r)) deallocate(tmp_r)
    if (allocated(tmp_max_r)) deallocate(tmp_max_r)
    if (allocated(tmp_min_r)) deallocate(tmp_min_r)


    !write(0,*) 'EXECUTING INSERT ROUTINE:'
    !write(outunit,'(a,10(1x,g18.8))') 'Average3:',min_dist,max_dist
    !do in = 1,av_group_cnt
    !   write(0,'(a,i10,10(1x,g18.8))') 'AV:',in,av_s(in),av_r(in),av_max_r(in),av_min_r(in)
    !end do


  end subroutine insert_limiter_point



  subroutine calculate_ray_limiter_surface
    !use common_utilities
    implicit none

    integer :: in,if,ierr
    !integer :: node_cnt
    real,allocatable :: surf_sep(:),av_group(:)
    !real,allocatable :: surf_ang(:)

    real :: min_sep
    real,parameter :: max_fact = 5.0
    integer,parameter :: max_block_size = 100
    integer :: it, ix
    real :: tmp_r,tmp_s,tmp_fl,tmp_int
    real :: test_sep,last_sep
    logical :: leave
    integer :: group_id
    integer :: av_cnt

    integer :: tmp_wall_cnt, tmp_tan_cnt,last_face,av_tmp_tan_cnt
    integer :: double_tan,double_wall
    real :: dir, last_dir

    integer :: index_offset

    real*8 :: av_type_org


    ! new variables
    logical :: done
    integer :: tan_cnt,bump
    integer :: sfl, efl, ifl

    !real :: ,last_ang,test_ang
    !real,parameter :: PI=3.141592654
    !real,parameter :: ang_limit = PI/18.0
    !real,parameter :: raddeg=57.29577952 

    !
    ! This routine goes through the intersection data and orders them so that when connected they form a piece-wise continuous description of the 
    ! intersection surface. The basic assumption is that the intersection data will be sortable by S coordinate ... that for a given specified S there is 
    ! only one surface intersection - thus no "islands" in the mesh formed by the intersections. Due to lack of data or less common situations this assumption
    ! may not always be true so the code will adjust the intersections as required. 
    !
    ! In addition, this code will also determine how to connect intersections on different limiters ... these would usually be identified by a "tangency" point 
    ! derivative change in the surface definition forming the bottom side of the tile as opposed to the peaks at limiter tips. Three methods can be used. 
    ! 1) connect the surface intersections
    ! 2) add a surface intersection point at the farthest out location
    ! 3) remove all intersections from the surface until the farther in intersection. 
    !
    !
    ! A second function of this code is to collect, order and create a reference list for all tangency points. If the rho (R) separation of tangency points is less 
    ! than a minimum value these will be moved onto the same tangency line. This data is used in grid generation to guide the choice of rho (R) values for the sides
    ! of the grid polygons
    !
    ! D.Elder 
    !
    ! The RAY data can not be sorted on S ... probably this is true for the CASTEM data as well but it wasn't clear due to other concerns with the CASTEM data. 
    !
    !
    ! Possible surface construction algorithm. 
    ! 1) start at minimum rho, minimum S point in the intersection data (alternatively the corner point if the subset specifiers are being used). 
    ! 2) March outward until the bottom edge of the first/nearest limiter bump
    ! 3) Connect to first bump
    ! 4) Follow bump around until reaching the farthest out point of the next bump
    ! 5) Connect to the next bump in order ... this will be the next indexed intersection on the farthest out rho line of the current bump .. if
    !    there is no further index then move in one field line at a time looking for an intersection to the next limiter. 
    ! 6) Continue until reaching the last limiter where there are no further intersection points. Connect to the field line starting points in this case and continue
    !    to the starting rho value and maximum S value. 
    !
    ! Several intersection data traversal routines are required. 
    ! A) Next inward intersection point - ends with tangency
    ! B) Next outward intersection point - ends at wall
    ! C) Next intersection point on next limiter - may be outward or inward. 
    ! 
    ! Might be worthwhile calculating the bump order first - bump order is defined by tangency points since each bump must have one (I hope).
    ! Need to calculate the bump index of the tangency points. 






    ! run through all intersections - sort node_list based on S values
    !
    ! Alternatively ... go through the intersection points and join each to its nearest neighbour ... this should work given the grid resolution ... then double check 
    ! by verifying the Lc coordinate monotonicity 
    !
    ! Start with an initial point at min_dist, min_lc and end at min_dist, max_lc
    ! This means the number of nodes in the list should be the total number of intersections + 2 ... the first and last nodes will have no cell reference. 
    ! Using just the closest will not work - especially for the big gaps ... need to go back to sorting by S value then check distances to see if there are any fix ups
    ! needed. 



    if (rg_int_win_mins.ne.rg_int_win_maxs) then 
       filter_intersections = .true.
    else
       filter_intersections = .false.
    endif



    ! Initialize averaging from global variable


    ! Generate list of tangency points including S, RHO and bump information

    call find_ray_tangency

    ! map start and end rows of all bumps

    call map_bumps

    ! now we need to use this data to map the limiter surface into the arrays used by the grid generator. Keep in mind that subsets of the data can be specified. 
    ! The tangency points are used to determine which bumps are present in any specific subset. 

    av_group_cnt = tot_n_intsects



    if (av_group_cnt.gt.0) then 
       if (allocated(av_s)) deallocate(av_s)
       allocate(av_s(av_group_cnt),stat=ierr)

       if (allocated(av_r)) deallocate(av_r)
       allocate(av_r(av_group_cnt),stat=ierr)

       if (allocated(av_min_r)) deallocate(av_min_r)
       allocate(av_min_r(av_group_cnt),stat=ierr)

       if (allocated(av_max_r)) deallocate(av_max_r)
       allocate(av_max_r(av_group_cnt),stat=ierr)

       if (allocated(av_angle)) deallocate(av_angle)
       allocate(av_angle(av_group_cnt),stat=ierr)

       if (allocated(av_metric)) deallocate(av_metric)
       allocate(av_metric(av_group_cnt),stat=ierr)

       if (allocated(av_type)) deallocate(av_type)
       allocate(av_type(av_group_cnt),stat=ierr)

       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_ ARRAYS:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    write(0,'(a,i9)') 'AV_GROUP_CNT:',av_group_cnt


    av_cnt = 0

    if (filter_intersections) then 

       ! add end faces - defined in input file
       av_cnt = av_cnt + 1

       av_r(av_cnt) = rg_minr
       av_s(av_cnt) = rg_mins
       av_min_r(av_cnt) = min_dist
       av_max_r(av_cnt) = min_dist

       av_type(av_cnt) = RAY_END
       av_metric(av_cnt) = -1
       av_angle(av_cnt)  = -1


    write(0,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)


       av_cnt = av_cnt + 1

       av_r(av_cnt) = rg_maxr
       av_s(av_cnt) = rg_mins
       av_min_r(av_cnt) = min_dist
       av_max_r(av_cnt) = min_dist

       av_type(av_cnt) = RAY_END
       av_metric(av_cnt) = -1
       av_angle(av_cnt)  = -1

    write(0,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)


!    End point addition is not needed except for intersection subsets
! 
!
!    else
!
!       ! add end points
!       av_cnt = av_cnt + 1
!
!       av_r(av_cnt) = min_dist
!       av_s(av_cnt) = min_lc
!       av_min_r(av_cnt) = min_dist
!       av_max_r(av_cnt) = min_dist
!
!       av_type(av_cnt) = RAY_END
!       av_metric(av_cnt) = -1
!       av_angle(av_cnt)  = -1
!
!    write(0,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
!    write(6,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
!
!
    endif


    ! Map out following the field line start position until reaching the end row of the first bump

    done = .false.

    tan_cnt = 1


    do while (.not.done)

       bump = av_tan_bump(tan_cnt)

       write (0,'(a,10(1x,i8))') 'BUMP_ASSYMETRY:', bump, bump_inf(bump,1),bump_inf(bump,2),bump_inf(bump,3),bump_inf(bump,2)-bump_inf(bump,3)
       write (6,'(a,10(1x,i8))') 'BUMP_ASSYMETRY:', bump, bump_inf(bump,1),bump_inf(bump,2),bump_inf(bump,3),bump_inf(bump,2)-bump_inf(bump,3)

       sfl = bump_inf(bump,1)
       efl = bump_inf(bump,2)

       ! Assign start
       if (tan_cnt.eq.1.and.(.not.filter_intersections)) then 
          do ifl = 1,efl
             call assign_ray_surf(ifl,bump,0,av_cnt)
    write(0,'(a,2i8,10(1x,g18.8))') 'RS:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RS:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)


          end do
       end if


       ! Do first half of bump (first intersections + tangency)
       ! always want RAY_ENTER intersection - no matter what assymmetries
       do ifl = efl,sfl+1,-1
          call assign_ray_surf(ifl,bump,1,av_cnt)

    write(0,'(a,2i8,10(1x,g18.8))') 'RF:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RF:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)

       end do

       ! add tangency point
          call assign_ray_surf(sfl,bump,3,av_cnt)
          ! assign the tan index into the av arrays 
          av_tan_ind(tan_cnt) = av_cnt

    write(0,'(a,2i8,10(1x,g18.8))') 'RT:',sfl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RT:',sfl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)

       ! allow for assymmetric limiters
       efl = bump_inf(bump,3)

       ! Do second half of bump (second intersections)

       do ifl = sfl+1,efl
          ! if beyond first half intersections then request only RAY_EXIT  point for bump on field line
          call assign_ray_surf(ifl,bump,2,av_cnt)

    write(0,'(a,2i8,10(1x,g18.8))') 'RB:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RB:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)

       end do

       ! Transition to next bump or END

       ! Assign end
       if (tan_cnt.eq.av_tan_cnt.and.(.not.filter_intersections)) then 
          do ifl = efl,1,-1
             call assign_ray_surf(ifl,bump,4,av_cnt)

    write(0,'(a,2i8,10(1x,g18.8))') 'RE:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RE:',ifl,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)

          end do

          done = .true.
       else
          tan_cnt = tan_cnt + 1

       end if
    end do



    if (filter_intersections) then 

       ! add end faces - defined in input file
       av_cnt = av_cnt + 1
       av_r(av_cnt) = rg_maxr
       av_s(av_cnt) = rg_maxs
       av_min_r(av_cnt) = min_dist
       av_max_r(av_cnt) = min_dist

       av_type(av_cnt) = RAY_END
       av_metric(av_cnt) = -1
       av_angle(av_cnt)  = -1

    write(0,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)


       av_cnt = av_cnt + 1
       av_r(av_cnt) = rg_minr
       av_s(av_cnt) = rg_maxs
       av_min_r(av_cnt) = min_dist
       av_max_r(av_cnt) = min_dist

       av_type(av_cnt) = RAY_END
       av_metric(av_cnt) = -1
       av_angle(av_cnt)  = -1

    write(0,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
    write(6,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)


!    End point addition is not needed except for intersection subsets
!    else
!
!       ! add end points
!       av_cnt = av_cnt + 1
!       av_r(av_cnt) = min_dist
!       av_s(av_cnt) = max_lc
!       av_min_r(av_cnt) = min_dist
!       av_max_r(av_cnt) = min_dist
!
!       av_type(av_cnt) = RAY_END
!       av_metric(av_cnt) = -1
!       av_angle(av_cnt)  = -1
!
!    write(0,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
!    write(6,'(a,2i8,10(1x,g18.8))') 'RP:',0,av_cnt,av_r(av_cnt),av_s(av_cnt),av_type(av_cnt)
!
!
    endif



    write(0,*) 'AV_CNT, av_group_cnt =',av_cnt, av_group_cnt
    write(6,*) 'AV_CNT, av_group_cnt =',av_cnt, av_group_cnt


    av_group_cnt = av_cnt


    do in = 1,av_group_cnt
       write(6,'(a,i8,10(1x,g18.8))') 'LIM:',in,av_r(in),av_s(in),av_type(in)

    end do

    !write(outunit,'(a,10(1x,g18.8))') 'Average1:',min_dist,max_dist
    !do in = 1,av_group_cnt
    !   write(0,'(a,i10,10(1x,g18.8))') 'AV:',in,av_s(in),av_r(in),av_max_r(in),av_min_r(in)
    !end do

    ! go through the data and remove degeneracies in tangency and wall points
    ! Assume for now that RAY data does not have this problem

    !dir = 1.0   
    !last_face = 1


    !in = 0
    !do while (in.lt.av_group_cnt-1)
    !do in = 1,av_group_cnt-1
    !  in = in+1
    !  last_dir = dir
    !  !write(outunit,'(a,2i8)') 'Index:',in,av_group_cnt
    !  dir = av_r(in+1)-av_r(in)
    !
    !  if (dir.eq.0.and.last_face.eq.1) then 
    !     ! fix wall degeneracy
    !     call insert_limiter_point(in,last_face)
    !     last_face = 2
    !  elseif (dir.eq.0.and.last_face.eq.2) then 
    !     ! fix tangency degeneracy
    !     call insert_limiter_point(in,last_face)
    !     last_face = 1
    !  elseif (last_dir.eq.0.and.last_face.eq.1) then 
    !     last_face = 2
    !  elseif (last_dir.eq.0.and.last_face.eq.2) then 
    !     last_face = 1
    !  elseif (dir.lt.0.and.last_dir.gt.0) then
    !     last_face = 2
    !  elseif (dir.gt.0.and.last_dir.lt.0) then
    !     last_face = 1
    !  elseif (dir.gt.0.and.last_dir.gt.0) then
    !     last_face = 1
    !  elseif (dir.lt.0.and.last_dir.lt.0) then
    !     last_face = 2
    !  endif

    !  !write(outunit,'(a,2i6,10(1x,g18.8))') 'FIND:',in,last_face,dir,last_dir,av_r(in),av_s(in)

    ! end do







    ! Now run through and identify IN, OUT and tangency points - tangency points only occur at the top of a limiter ... not sure if IN/OUT are required


    ! allocate av_type
    !if (av_group_cnt.gt.0) then 
    !   if (allocated(av_type)) deallocate(av_type)
    !   allocate(av_type(av_group_cnt),stat=ierr)
    !   if (ierr.ne.0) then 
    !      call errmsg('ALLOCATION ERROR:AV_TYPE:IERR =',ierr)
    !      stop
    !   endif
    !else
    !   return
    !endif

    dir = 1.0   
    ! Tangency has already been defined from RAY data - however, leave code for now to check this 
    av_tmp_tan_cnt = 0
    av_wall_cnt = 0

    last_face = 1

    do in = 1,av_group_cnt-1
       last_dir = dir
       dir = av_r(in+1)-av_r(in)
       av_type_org = av_type(in)

       if (dir.eq.0.and.last_face.eq.1) then 
          av_type(in) = WALL
          av_wall_cnt = av_wall_cnt + 1
          last_face = 1
       elseif (dir.eq.0.and.last_face.eq.2) then 
          av_type(in) = TANGENCY
          av_tmp_tan_cnt = av_tmp_tan_cnt + 1
          last_face = 2
       elseif (last_dir.eq.0.and.last_face.eq.1) then 
          av_type(in) = WALL
          av_wall_cnt = av_wall_cnt + 1
          if (dir.lt.0) then 
             last_face = 2
          else
             last_face = 1
          endif
       elseif (last_dir.eq.0.and.last_face.eq.2) then 
          av_type(in) = TANGENCY
          av_tmp_tan_cnt = av_tmp_tan_cnt + 1
          if (dir.gt.0) then 
             last_face = 1
          else
             last_face = 2
          endif
       elseif (dir.lt.0.and.last_dir.gt.0) then
          av_type(in) = WALL
          av_wall_cnt = av_wall_cnt + 1
          last_face = 2
       elseif (dir.gt.0.and.last_dir.lt.0) then
          av_type(in) = TANGENCY
          av_tmp_tan_cnt = av_tmp_tan_cnt + 1
          last_face = 1
       elseif (dir.gt.0.and.last_dir.gt.0) then
          av_type(in) = SURFACE_START
          last_face = 1
       elseif (dir.lt.0.and.last_dir.lt.0) then
          av_type(in) = SURFACE_END
          last_face = 2
       endif

       write(0,'(a,3i6,10(1x,g18.8))')       'DIR:',av_tmp_tan_cnt,in,last_face,dir,last_dir,av_type_org,av_type(in),av_r(in),av_s(in)
       write(outunit,'(a,3i6,10(1x,g18.8))') 'DIR:',av_tmp_tan_cnt,in,last_face,dir,last_dir,av_type_org,av_type(in),av_r(in),av_s(in)

    end do


       write(0,'(a,3i6,10(1x,g18.8))')       'END:',av_tmp_tan_cnt,in,last_face,dir,last_dir,av_type_org,av_type(in),av_r(in),av_s(in)
       write(outunit,'(a,3i6,10(1x,g18.8))') 'END:',av_tmp_tan_cnt,in,last_face,dir,last_dir,av_type_org,av_type(in),av_r(in),av_s(in)


    av_type(av_group_cnt) = SURFACE_END





    !
    ! Tangency points
    !

    !if (av_tan_cnt.gt.0) then 
    !   if (allocated(av_tan_r)) deallocate(av_tan_r)
    !   allocate(av_tan_r(av_tan_cnt),stat=ierr)
    !   if (ierr.ne.0) then 
    !      call errmsg('ALLOCATION ERROR:AV_TAN_R:IERR =',ierr)
    !      stop
    !   endif
    !   av_tan_r = 0.0
    !else
    !   return
    !endif



    !if (av_tan_cnt.gt.0) then 
    !   if (allocated(av_tan_s)) deallocate(av_tan_s)
    !   allocate(av_tan_s(av_tan_cnt),stat=ierr)
    !   if (ierr.ne.0) then 
    !      call errmsg('ALLOCATION ERROR:AV_TAN_S:IERR =',ierr)
    !      stop
    !   endif
    !   av_tan_s = 0.0
    !else
    !   return
    !endif


    ! Need index to limiter surface vertex as well
    !if (av_tan_cnt.gt.0) then 
    !   if (allocated(av_tan_ind)) deallocate(av_tan_ind)
    !   allocate(av_tan_ind(av_tan_cnt),stat=ierr)
    !   if (ierr.ne.0) then 
    !      call errmsg('ALLOCATION ERROR:AV_TAN_IND:IERR =',ierr)
    !      stop
    !   endif
    !   av_tan_ind = 0
    !else
    !   return
    !endif

    !
    ! Tangency point radial order
    !

    !if (av_tan_cnt.gt.0) then 
    !   if (allocated(tan_ord_r)) deallocate(tan_ord_r)
    !   allocate(tan_ord_r(av_tan_cnt),stat=ierr)
    !   if (ierr.ne.0) then 
    !      call errmsg('ALLOCATION ERROR:TAN_ORD_R:IERR =',ierr)
    !      stop
    !  endif
    !   tan_ord_r = 0.0
    !else
    !   return
    !endif



    !
    ! wall points
    !

    if (av_wall_cnt.gt.0) then 

       if (allocated(av_wall_s)) deallocate(av_wall_s)
       allocate(av_wall_s(av_wall_cnt),stat=ierr)

       if (allocated(av_wall_r)) deallocate(av_wall_r)
       allocate(av_wall_r(av_wall_cnt),stat=ierr)

       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_WALL ARRAYS:IERR =',ierr)
          stop
       endif
       av_wall_s = 0.0
       av_wall_r = 0.0
    else
       return
    endif


    tmp_tan_cnt = 0
    tmp_wall_cnt = 0

    ! Collect tangency and wall data - check to make sure there are no two consective points with the same R value - i.e. flat tangency point
    double_tan = 0
    double_wall = 0

    do in = 1,av_group_cnt

       if (av_type(in).eq.TANGENCY) then 
          tmp_tan_cnt = tmp_tan_cnt + 1

          if (tmp_tan_cnt.le.av_tan_cnt) then 
             ! record r,s and index into surface array for the tangency point
             !av_tan_r(tmp_tan_cnt) = av_r(in)
             !av_tan_s(tmp_tan_cnt) = av_s(in)
             !av_tan_ind(tmp_tan_cnt) = in
             write(outunit,'(a,i6,10(1x,g18.8))') 'TAN :',tmp_tan_cnt,in,av_r(in),av_s(in),av_tan_r(tmp_tan_cnt),av_tan_s(tmp_tan_cnt)
             if (in.ne.av_group_cnt) then 
                if(av_type(in+1).eq.TANGENCY) then 
                   double_tan = double_tan+1
                endif
             endif
          !else
          !   write(0,'(a,2i10)') 'TANGENCIES:',av_tan_cnt, tmp_tan_cnt
          !   CALL ERRMSG('Too many tangency points found:',tmp_tan_cnt)
          !   stop
          endif

       endif

       if (av_type(in).eq.WALL) then 
          tmp_wall_cnt = tmp_wall_cnt + 1
          if (tmp_wall_cnt.le.av_wall_cnt) then 
             av_wall_r(tmp_wall_cnt) = av_r(in)
             av_wall_s(tmp_wall_cnt) = av_s(in)
             write(outunit,'(a,i6,10(1x,g18.8))') 'WALL:',tmp_wall_cnt,in,av_r(in),av_s(in)
             if (in.ne.av_group_cnt) then 
                if(av_type(in+1).eq.WALL) then 
                   double_wall = double_wall+1
                endif
             endif
          else
             CALL ERRMSG('Too many wall points found:',tmp_wall_cnt)
             stop
          endif
       endif

    end do

    write(0,'(a,3i10)') 'TANGENCIES:',av_tan_cnt, tmp_tan_cnt,av_tmp_tan_cnt
    write(outunit,'(a,3i10)') 'TANGENCIES:',av_tan_cnt, tmp_tan_cnt,av_tmp_tan_cnt



    write(outunit,*) 'DOUBLES:',double_tan,double_wall

    !write(outunit,'(a,6i10)') 'Tan:',av_tan_cnt, av_wall_cnt,double_tan,double_wall

    !do in = 1,av_wall_cnt
    !   write(0,'(a,i6,10(1x,g18.8))') 'Wall:',in,av_wall_r(in),av_wall_s(in)
    !end do

    !do in = 1,av_tan_cnt
    !   write(0,'(a,i6,10(1x,g18.8))') 'Tan :',in,av_tan_r(in),av_tan_s(in)
    !end do

    !write(outunit,'(a,10(1x,g18.8))') 'Average2:',min_dist,max_dist
    !do in = 1,av_group_cnt
    !   write(0,'(a,i10,10(1x,g18.8))') 'AV:',in,av_s(in),av_r(in),av_max_r(in),av_min_r(in)
    !end do


    ! Find max and min from the reduced intersection data ... loop through revised wall to get this data


    r_limiter_max = min_dist -1.0
    r_limiter_min = max_dist +1.0

    s_limiter_max = min_lc -1.0
    s_limiter_min = max_lc +1.0


    do in = 1,av_group_cnt

       if (av_r(in).lt.min_dist-1.0) then 
          write(0,*) 'PROB:1:',in,min_dist,av_r(in)
       elseif (av_r(in).gt.max_dist+1) then 
          write(0,*) 'PROB:2:',in,max_dist,av_r(in)
       endif

       r_limiter_max = max(r_limiter_max,av_r(in))
       r_limiter_min = min(r_limiter_min,av_r(in))

       s_limiter_max = max(s_limiter_max,av_s(in))
       s_limiter_min = min(s_limiter_min,av_s(in))

    end do

    write(outunit,'(a,10(1x,g18.8))') 'R,S RANGES:',r_limiter_min,r_limiter_max,s_limiter_min,s_limiter_max
    !write(0,'(a,10(1x,g18.8))') 'R,S RANGES:',r_limiter_min,r_limiter_max,s_limiter_min,s_limiter_max



    ! Ok - now need to check for tangency points that are very close together and adjust the wall so that these points lie on the same line. 
    !    - in this case the larger r value will be assigned the smaller ... the minimum separation will be quite small compared to the simulation space 
    !      so this may not have a large impact
    !
    ! What is the minimum radial tangency separation distance?
    ! 
    !    
    ! set separation limit to r_limiter_max - r_limiter_min / 100 ... or just choose an absolute value like 0.001m 
    !
    min_tan_sep = (r_limiter_max-r_limiter_min)/100.0



    ! sort tangency points by r order to move r points that are too close together
    call sort_arrays(0,av_tan_cnt,av_tan_r,av_tan_s,av_tan_ind)


    ! go through and adjust tangency points as needed. 
    ! since wall points are not required to be a polygon vertex this processing isn't required for wall extremities

    do in = 1,av_tan_cnt-1
       if ((av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep) then 
          ! move tangency point in to match

          write(outunit,'(a,2i10,l8,15(1x,g18.8))') 'Reduce tan:',in,in+1,(av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep,av_tan_r(in+1),av_tan_r(in),av_tan_r(in+1)-av_tan_r(in),min_tan_sep,av_tan_s(in),av_tan_s(in+1)
          !write(0,'(a,2i10,l8,15(1x,g18.8))') 'Reduce tan:',in,in+1,(av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep,av_tan_r(in+1),av_tan_r(in),av_tan_r(in+1)-av_tan_r(in),min_tan_sep,av_tan_s(in),av_tan_s(in+1)

          av_tan_r(in+1) = av_tan_r(in)
          av_r(int(av_tan_ind(in+1))) = av_r(int(av_tan_ind(in)))

       endif

    end do


    ! put tangency points back into s order - this should be the same order as before the previous re-order - this is needed for generating cells while the r 
    ! r order is useful for generating rows
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,tan_ord_r,av_tan_ind)

    ! record the S order of the points 

    do in = 1,av_tan_cnt
       tan_ord_r(in) = in
       write(outunit,'(a,i8,10(1x,g18.8))') 'TAN_ORD_R2:',in,av_tan_r(in),av_tan_s(in),tan_ord_r(in)
    end do

    ! Sort the S indices by R to get tan_ord_r referencing the R order of the data while sorted by S
    call sort_arrays(0,av_tan_cnt,av_tan_r,av_tan_s,tan_ord_r,av_tan_ind)

    ! Re-sort to S - leaving tan_ord_r unchanged
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,av_tan_ind)


    do in = 1,av_tan_cnt
       write(outunit,'(a,i8,10(1x,g18.8))') 'TANS:',in,av_tan_r(in),av_tan_s(in),av_tan_ind(in)
       !write(0,'(a,i8,10(1x,g18.8))') 'TANS:',in,av_tan_r(in),av_tan_s(in),av_tan_ind(in)
    end do

    do in = 1,av_wall_cnt
       write(outunit,'(a,i8,10(1x,g18.8))') 'WALL:',in,av_wall_r(in),av_wall_s(in)
       !write(0,'(a,i8,10(1x,g18.8))') 'WALL:',in,av_wall_r(in),av_wall_s(in)
    end do


  end subroutine calculate_ray_limiter_surface



  subroutine find_ray_tangency

    implicit none
    integer :: tan_cnt,in,is
    integer :: ierr

    av_tan_cnt = n_tangency


    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_r)) deallocate(av_tan_r)
       allocate(av_tan_r(av_tan_cnt),stat=ierr)

       if (allocated(av_tan_s)) deallocate(av_tan_s)
       allocate(av_tan_s(av_tan_cnt),stat=ierr)

       if (allocated(av_tan_bump)) deallocate(av_tan_bump)
       allocate(av_tan_bump(av_tan_cnt),stat=ierr)

       if (allocated(av_tan_fl)) deallocate(av_tan_fl)
       allocate(av_tan_fl(av_tan_cnt),stat=ierr)

       if (allocated(av_tan_ind)) deallocate(av_tan_ind)
       allocate(av_tan_ind(av_tan_cnt),stat=ierr)

       if (allocated(tan_ord_r)) deallocate(tan_ord_r)
       allocate(tan_ord_r(av_tan_cnt),stat=ierr)

       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:AV_TAN ARRAYS FOR RAY DATA:IERR =',ierr)
          stop
       endif
       av_tan_r = 0.0
       av_tan_s = 0.0
       av_tan_bump = 0.0
       av_tan_fl = 0.0
       av_tan_ind = 0.0
    else
       return
    endif



    tan_cnt = 0

    do in = 1,n_field_lines
       do is = 1,field_line(in)%int_tot
          if (field_line(in)%int_data(is)%int_type.eq.RAY_TAN) then 
             tan_cnt = tan_cnt + 1
             if (tan_cnt.le.av_tan_cnt) then 
                av_tan_r(tan_cnt) = field_line(in)%dist
                av_tan_s(tan_cnt) = field_line(in)%int_data(is)%lc
                av_tan_bump(tan_cnt) = field_line(in)%int_data(is)%bump
                av_tan_fl(tan_cnt) = in
                av_tan_ind(tan_cnt) = is
             else
                call errmsg('FIND_RAY_TANGENCY:TOO MANY TANGENT POINTS FOUND',av_tan_cnt)
                stop 'FIND_RAY_TANGENCY'
             endif

          endif
       end do
    end do

    if (tan_cnt.ne.av_tan_cnt) then
       call errmsg('FIND_RAY_TANGENCY:MISMATCH IN TANGENCY COUNT',tan_cnt)
       stop 'FIND_RAY_TANGENCY'
    endif

    ! sort tangency point arrays by S
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,av_tan_ind,av_tan_bump,av_tan_fl)

    ! write the S ordered array index into another array
    do in = 1,av_tan_cnt
       tan_ord_r(in) = in
       write(outunit,'(a,i8,10(1x,g18.8))') 'TAN_ORD_RA:',in,av_tan_r(in),av_tan_s(in),tan_ord_r(in)
    end do

    ! Sort the S indices by R to get tan_ord_r referencing the R order of the data while sorted by S
    call sort_arrays(0,av_tan_cnt,av_tan_r,av_tan_s,av_tan_ind,av_tan_bump,av_tan_fl,tan_ord_r)

    ! Re-sort to S - leaving tan_ord_r unchanged
    ! The av_tan_bump array should now contain the bump order
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,av_tan_ind,av_tan_bump,av_tan_fl)

    do in = 1,av_tan_cnt
       write(outunit,'(a,i8,10(1x,g18.8))') 'TAN_ORD_RB:',in,av_tan_r(in),av_tan_s(in),tan_ord_r(in)
    end do

  end subroutine find_ray_tangency


  subroutine assign_ray_surf(fl,bump,intn,av_cnt)
    implicit none
    integer :: fl,bump,intn,av_cnt
    integer :: int_id

    ! intn defines which intersection on a specific field line is to be loaded into the 
    ! av_ arrays ... if the specified intn is unavailable the code just continues to execute without inserting 
    ! an intersection point
    !
    ! intn : 0 = field line start position
    !        1 = RAY_ENTER intersection with specific bump (or tangency)
    !        2 = RAY_EXIT intersection with specific bump
    !        3 = RAY_TAN intersection for specific bump
    !        4 = field line end position 
    !

    ! assign all the quantities for the av_ arrays
    ! av_s, av_r, av_min_r, av_max_r, av_metric, av_angle, av_type

    call get_ray_intersection(int_id,fl,bump,intn)


    ! Intersection ID has been assigned
    if (int_id.ge.0) then 

       av_cnt = av_cnt + 1
       av_s(av_cnt)      = field_line(fl)%int_data(int_id)%lc
       av_r(av_cnt)      = field_line(fl)%dist
       av_type(av_cnt)   = field_line(fl)%int_data(int_id)%int_type
       av_metric(av_cnt) = field_line(fl)%int_data(int_id)%metric
       av_angle(av_cnt)  = field_line(fl)%int_data(int_id)%angle
       av_min_r(av_cnt)  = av_r(av_cnt)
       av_max_r(av_cnt)  = av_r(av_cnt)

    endif

  end subroutine assign_ray_surf



  subroutine get_ray_intersection(int_id,fl,bump,intn)
    implicit none
    integer :: int_id,fl,bump,intn
    ! return index to desired intersection point
    ! 
    ! intn : 0 = field line start position
    !        1 = RAY_ENTER intersection with specific bump (or tangency)
    !        2 = RAY_EXIT intersection with specific bump
    !        3 = RAY_TAN intersection for specific bump
    !        4 = field line end position 
    !

    integer :: int_cnt,int_num,icnt
    logical :: done
    integer :: required_intsect

    done = .false.
    int_id = -1

    int_num = field_line(fl)%int_tot

    ! Note: required_intsect is used either for position or type depending on use
    if (intn.eq.0) then 
       required_intsect = 1
    elseif (intn.eq.1) then 
       required_intsect = RAY_ENTER
    elseif (intn.eq.2) then 
       required_intsect = RAY_EXIT
    elseif (intn.eq.3) then 
       required_intsect = RAY_TAN
    elseif (intn.eq.4) then 
       required_intsect = int_num
    endif


    ! return the field line start position - if available - type 5
    if (intn.eq.0.or.intn.eq.4) then 

       if (field_line(fl)%int_data(required_intsect)%int_type.eq.RAY_END) then 
          int_id = required_intsect
       endif

    elseif(intn.eq.1.or.intn.eq.2.or.intn.eq.3) then 

       icnt = 0
       do while (.not.done)
          icnt = icnt + 1

          if (field_line(fl)%int_data(icnt)%bump.eq.bump.and.field_line(fl)%int_data(icnt)%int_type.eq.required_intsect) then 
             if (.not.filter_intersections.or.&
                  (filter_intersections.and. &
                  (field_line(fl)%int_data(icnt)%lc.ge.rg_int_win_mins).and.&
                  (field_line(fl)%int_data(icnt)%lc.le.rg_int_win_maxs))) then 

                int_id = icnt
                done = .true.

             endif

          endif

          if (icnt.eq.int_num) done=.true.

       end do

    endif



  end subroutine get_ray_intersection



  subroutine map_bumps
    implicit none


    integer :: bump_id,last_bump_id
    integer :: ierr,in,int_cnt,is


    nbumps = maxval(av_tan_bump)


    if (nbumps.gt.0) then 
       if (allocated(bump_inf)) deallocate(bump_inf)
       allocate(bump_inf(nbumps,3),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:BUMP_INF:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    bump_inf= 0

    do in = 1,n_field_lines
       last_bump_id = 0
       int_cnt = 0
       do is = 1,field_line(in)%int_tot
          bump_id = field_line(in)%int_data(is)%bump

          ! calculate the start and end rows for each bump
          if (bump_id.ne.0) then 

             if (bump_inf(bump_id,1).eq.0) then 
                bump_inf(bump_id,1) = in
                if (field_line(in)%int_data(is)%int_type.ne.RAY_TAN) then
                   write(0,'(a,4i10,10(1x,g18.8))') 'BUMP START NOT TANGENCY:',bump_id,in,is,field_line(in)%int_data(is)%int_type

                endif
             endif

             ! This allows for assymmetric bumps - the assymmetric sections will only have one intersection (either an entrance or exit depending)
             if (field_line(in)%int_data(is)%int_type.eq.RAY_ENTER) then 
                bump_inf(bump_id,2) = in
             elseif (field_line(in)%int_data(is)%int_type.eq.RAY_EXIT) then
                bump_inf(bump_id,3) = in
             endif

          endif
       end do
    end do

  end subroutine map_bumps










  subroutine generate_grid
    implicit none

    ! This routine generates a grid ... the characteristics may need to be tweaked for each case.
    ! Options can include making the cell vertices match row to row ...
    !
    ! Steps - determine R cell boundary rows ... use tangency points first ... 
    !         then add rows until row separation is at or below the specified value
    !       - maximum radial separation is r_max-rmin/10.0 (?) ... or maybe r_intersection_max - r_min/10.0
    !
    !       - then ... moving along each row ... allocate cell vertices ... tangency points should be pre-allocated
    !       - split each row at limiter surfaces - surface intersection points are the end-points for each cell group 
    !       - collect vertices for each cell and indices - index and order them 
    !

    !
    ! For first grid ... use all tangency points in both r and s ... then refine grid by imposing maximum cell size ... add rows/columns when the cells will be too large.
    !


    ! generate initial S and R positions from tangency information. ... then fill in the data from intersection information.
    !
    !
    ! generate the rows for the grid ... find number of unique tangency R values
    !
    ! Allocate initial r and s arrays at 3 * number of tangency points
    !


    real*8,allocatable :: r_bnds(:)   !,s_bnds(:)
    integer :: init_grid_size
    integer :: in,in1,it,is

    real*8 :: max_s_sep, max_r_sep

    real*8, allocatable :: ints(:,:),int_type(:,:),int_cnt(:)

    real*8, allocatable :: vert_rec1(:),vert_rec2(:),vert_type1(:),vert_type2(:)

    real*8 :: s_start,s_end
    real*8 :: slen, ssep1
    real*8,allocatable :: ssep1a(:)
    integer :: ncells

    integer :: rbnd_cnt,rbnd_add,new_bnds
    !integer :: sbnd_cnt,rbnd_cnt,rbnd_add,new_bnds

    integer :: npts1a, npts1b,npts2a,npts2b
    integer :: vert_cnt1,vert_cnt2,ntot1,ntot2
    real*8 :: r1,r2,r_sep_dist
    integer :: ierr,max_surf_ints
    integer :: min_cells
    logical :: vert_set


    real*8 :: line_length, rfactor, lfactor  ! variables related to ring reduced grids 

    integer :: rings_processed
    integer :: start_count,end_count

    ! initialize module variables to global input values
    grid_option = rg_grid_opt
    max_s_sep = rg_max_s_sep
    max_r_sep = rg_max_r_sep
    min_cells = rg_min_cells

    !write(0,*) '1',rg_grid_opt,rg_block_av,rg_max_s_sep,rg_max_r_sep,rg_min_cells
    !write(0,*) '1a',grid_option, opt_block_av, max_s_sep,max_r_sep,min_cells,r_limiter_max,r_limiter_min
    !write(0,*) '1b',min_dist,max_dist,min_lc,max_lc
    !write(0,*) 'Calc1:',(r_limiter_max-r_limiter_min)/max_r_sep
    !write(0,*) 'Calc1a:',int((r_limiter_max-r_limiter_min)/max_r_sep)
    !write(0,*) 'Calc2:',3 * av_tan_cnt + 1
    !write(0,*) 'Calc3:',max(3 * av_tan_cnt + 1, int((r_limiter_max-r_limiter_min)/max_r_sep))


    init_grid_size = max(3 * av_tan_cnt + 1, int((r_limiter_max-r_limiter_min)/max_r_sep))

    !write(0,*) 'Init_grid_size:',init_grid_size,max_r_sep

    if (init_grid_size.gt.0) then 
       if (allocated(r_bnds)) deallocate(r_bnds)
       allocate(r_bnds(init_grid_size),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:R_BNDS:IERR =',ierr)
          stop
       endif
    else
       return
    endif


    !    if (init_grid_size.gt.0) then 
    !       if (allocated(s_bnds)) deallocate(s_bnds)
    !       allocate(s_bnds(init_grid_size),stat=ierr)
    !    else
    !       return
    !    endif

    rbnd_cnt = 0
    !    sbnd_cnt = 0


    do in = 1,av_tan_cnt
       in1 = tan_ord_r(in)

       if (in.ne.1) then 
          if (av_tan_r(in1).ne.r_bnds(rbnd_cnt)) then 
             rbnd_cnt = rbnd_cnt + 1
             if (rbnd_cnt.gt.size(r_bnds)) call grow_array(r_bnds,rbnd_cnt-1,rbnd_cnt+10)
             r_bnds(rbnd_cnt) = av_tan_r(in1)
             write(outunit,'(a,3i8,10(1x,g18.8))') 'Calc R_bnd:',in,in1,rbnd_cnt,r_bnds(rbnd_cnt),av_tan_r(in1),av_tan_s(in1),av_tan_s(in),r_bnds(rbnd_cnt-1)

          endif
       else
          rbnd_cnt = rbnd_cnt + 1
          if (rbnd_cnt.gt.size(r_bnds)) call grow_array(r_bnds,rbnd_cnt-1,rbnd_cnt+10)
          r_bnds(rbnd_cnt) = av_tan_r(in1)
          write(outunit,'(a,3i8,10(1x,g18.8))') 'Calc R_bnd:',in,in1,rbnd_cnt,r_bnds(rbnd_cnt),av_tan_r(in1),av_tan_s(in1),av_tan_s(in)
       endif

       !       sbnd_cnt = sbnd_cnt + 1
       !       s_bnds(sbnd_cnt) = av_tan_s(in)

    end do


    !
    ! Add locations that are not tangency points  .. then sort the arrays
    !

    !    sbnd_cnt = sbnd_cnt + 1
    !    s_bnds(sbnd_cnt) = s_limiter_min

    !    sbnd_cnt = sbnd_cnt + 1
    !    s_bnds(sbnd_cnt) = s_limiter_max

    rbnd_cnt = rbnd_cnt + 1
    if (rbnd_cnt.gt.size(r_bnds)) call grow_array(r_bnds,rbnd_cnt-1,rbnd_cnt+10)
    r_bnds(rbnd_cnt) = r_limiter_min

    ! put last r_bnd just inside the wall
    rbnd_cnt = rbnd_cnt + 1
    if (rbnd_cnt.gt.size(r_bnds)) call grow_array(r_bnds,rbnd_cnt-1,rbnd_cnt+10)
    r_bnds(rbnd_cnt) = maxval(av_tan_r) + (r_limiter_max - maxval(av_tan_r)) * 0.95



    call sort_arrays(0,rbnd_cnt,r_bnds)
    !call sort_arrays(0,sbnd_cnt,s_bnds)
    !call sort_arrays(1,sbnd_cnt,s_bnds)

    ! run through R bounds and insert extras in situations where the separation is too large. 


    rbnd_add = 0

    do in = 1,rbnd_cnt-1
       new_bnds = int((r_bnds(in+1)-r_bnds(in))/max_r_sep)
       if (new_bnds.ge.1) then 
          r_sep_dist = (r_bnds(in+1)-r_bnds(in))/(new_bnds+1)
          do it = 1,new_bnds
             rbnd_add = rbnd_add + 1
             if (rbnd_cnt+rbnd_add.gt.size(r_bnds)) call grow_array(r_bnds,rbnd_cnt+rbnd_add-1,rbnd_cnt+rbnd_add+10)
             r_bnds(rbnd_cnt+rbnd_add) = r_bnds(in) + r_sep_dist * it
          end do

       endif

    end do

    rbnd_cnt = rbnd_cnt + rbnd_add

    ! re-sort boundaries

    call sort_arrays(0,rbnd_cnt,r_bnds)

    ! write out r_bnds

    do in = 1,rbnd_cnt
       ! slmod - array bounds error on R_BNDS
       write(outunit,'(a,i8,10(1x,g18.8))') 'R_BNDS:',in,r_bnds(in),r_bnds(in)-r_bnds(max(in-1,1))
       !
    end do


    !
    ! Ok - now find intersection points with structure on this new grid for each R cell boundary
    !


    ! Allocate arrays to hold intersection data

    ! max surf ints needs to be 3 X tangency count since each tangency can generate 2 intersections plus the explicit addition of the tangency point
    max_surf_ints = 3 * av_tan_cnt + 2


    !write(0,'(a,2i6)') 'Parameters:',max_surf_ints,av_tan_cnt

    if (rbnd_cnt.gt.0.and.av_tan_cnt.gt.0) then 
       if (allocated(int_cnt)) deallocate(int_cnt)
       allocate(int_cnt(rbnd_cnt),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:INT_CNT:IERR =',ierr)
          stop
       endif
       int_cnt = 0
    else
       return
    endif

    if (rbnd_cnt.gt.0.and.av_tan_cnt.gt.0) then 
       if (allocated(ints)) deallocate(ints)
       allocate(ints(rbnd_cnt,max_surf_ints),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:INTS:IERR =',ierr)
          stop
       endif
       ints = 0
    else
       return
    endif

    if (rbnd_cnt.gt.0.and.av_tan_cnt.gt.0) then 
       if (allocated(int_type)) deallocate(int_type)
       allocate(int_type(rbnd_cnt,max_surf_ints),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:INT_TYPE:IERR =',ierr)
          stop
       endif
       int_type = 0
    else
       return
    endif

    int_cnt = 0.0
    int_type = 0.0
    ints = 0.0



    ! First r_bnd is assumed to only have 2 points - one at max and one at min
    ! These need to be set differently for subset grids

    if (filter_intersections) then 
       int_cnt(1) = 2
       ints(1,1) = rg_mins
       int_type(1,1) = SURFACE_START
       ints(1,2) = rg_maxs
       int_type(1,2) = SURFACE_END
    else
       int_cnt(1) = 2
       ints(1,1) = min_lc
       int_type(1,1) = SURFACE_START
       ints(1,2) = max_lc
       int_type(1,2) = SURFACE_END
    endif

    ! Find all the other intersections
    do in = 2,rbnd_cnt

       call find_field_line_intsects(max_surf_ints,r_bnds(in),int_cnt(in),int_type(in,:),ints(in,:))

    end do

    ! write out intersection points
    do in = 1,rbnd_cnt
       write(outunit,'(a,i8,10(1x,g18.8))')    'INTS  :',in,int_cnt(in)
       !write(0,'(a,i8,10(1x,g18.8))')    'INTS  :',in,int_cnt(in)
       do it = 1,int_cnt(in)
          write(outunit,'(a,2i8,10(1x,g18.8))') '  DATA:', in,it,int_type(in,it),ints(in,it),r_bnds(in)
       end do
    end do

    ! Now I should have an array containing all the intersection data needed for grid generation - whatever grid generation scheme is used. 

    ! Constructing a ring or rings of polygons requires 
    ! 1) the intersection data from two adjacent boundaries
    ! 2) could require any predefined vertices from adjacent rings if we want to match them up
    ! 3) we may want to propagate out the S position of tangency points so that these are already used when the tangency point is reached
    !
    ! Data required depends on the grid/polygon generation scheme
    ! Independent polygons on each ring only needs (1)
    ! Correlating polygons between rings requires more data. 

    ! option 1 ... generate polygons for rings along each pair of R data - ignore cells in other rings
    ! option 2 ... generate polygons for rings along each pair of R data - use cell corners from farther out - add cells where needed 
    ! option 3 ... generate polygons for rings along each pair of R data - tangency points are vertices on ALL rings with R < Rtan

    ! Ok - what base size to allocate for polygons/ring data

    ! rvertp ... zvertp ... nvertp ... nployp 
    ! nrs - number of rings
    ! nks(ir) - number of knots on each ring
    ! korpg(ik,ir) - polygon index for each cell
    ! rs,zs - cell center points
    ! bratio(ik,ir) ... set to 1 for now
    ! 
    ! base estimate of polygons ... can fine tune and reduce storage later
    ! nrings = ntangency * rbnd_cnt 
    ! nknots = nrings * 100 (?) 
    ! need maxnrs and maxnks calculated later

    ! initial estimates of maximum grid size

    max_nrings = (av_tan_cnt+1) * rbnd_cnt
    !max_nknots = max(200,int((max_lc-min_lc)/max_s_sep)*2)
    max_nknots = max(200,max((min_cells*(av_tan_cnt+1))*2,int((max_lc-min_lc)/max_s_sep)*2))
    ! This value is much too large
    !max_npoly = max_nrings * max_nknots
    max_npoly = max(2*rbnd_cnt * max_nknots,max_nrings * min_cells *4)

    write(0,'(a,7i12,3g18.8)') 'Max values:',av_tan_cnt,av_wall_cnt,rbnd_cnt,max_nrings,max_nknots,max_npoly,int((max_lc-min_lc)/max_s_sep),max_lc,min_lc,max_s_sep

    max_wall_segments = max_nrings * 2 + (av_tan_cnt +1) * 2 

    !
    ! Allocate polygon arrays
    ! rvp
    if (max_npoly.gt.0) then 
       if (allocated(rvp)) deallocate(rvp)
       allocate(rvp(4,max_npoly),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:RVP:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    ! zvp
    if (max_npoly.gt.0) then 
       if (allocated(zvp)) deallocate(zvp)
       allocate(zvp(4,max_npoly),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:ZVP:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    !
    ! nvp
    if (max_npoly.gt.0) then 
       if (allocated(nvp)) deallocate(nvp)
       allocate(nvp(max_npoly),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:NVP:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    ! default is all 4 sided polygons
    nvp = 4



    ! allocate ring data arrays
    ! rcen

    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(rcen)) deallocate(rcen)
       allocate(rcen(max_nknots,max_nrings),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:RCEN:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    ! zcen
    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(zcen)) deallocate(zcen)
       allocate(zcen(max_nknots,max_nrings),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:ZCEN:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    ! poly_ref
    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(poly_ref)) deallocate(poly_ref)
       allocate(poly_ref(max_nknots,max_nrings),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:POLY_REF:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    ! nknots
    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(nknots)) deallocate(nknots)
       allocate(nknots(max_nrings),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:NKNOTS:IERR =',ierr)
          stop
       endif
    else
       return
    endif


    ! vertex data for calling gen_ring

    if (max_nknots.gt.0) then 
       if (allocated(vert_rec1)) deallocate(vert_rec1)
       allocate(vert_rec1(max_nknots),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:VERT_REC1:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (max_nknots.gt.0) then 
       if (allocated(vert_rec2)) deallocate(vert_rec2)
       allocate(vert_rec2(max_nknots),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:VERT_REC2:IERR =',ierr)
          stop
       endif
    else
       return
    endif


    if (max_nknots.gt.0) then 
       if (allocated(vert_type1)) deallocate(vert_type1)
       allocate(vert_type1(max_nknots),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:VERT_TYPE1:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    if (max_nknots.gt.0) then 
       if (allocated(vert_type2)) deallocate(vert_type2)
       allocate(vert_type2(max_nknots),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:VERT_TYPE2:IERR =',ierr)
          stop
       endif
    else
       return
    endif


    ! group data into rings based on intersection and then pass the data to the ring generation routine


    ! initialization
    npoly = 0
    nrings = 0
    nknots = 0

    ! Depending on grid generation option - either calculate the initial vertices on the first line or not


    ! grid option 1 - free form - nothing needs to be done
    ! grid option 2 - gridded with tangency points - generate first set of vertices.
    !

    !
    ! Set up the initial vertex arrays
    !

    call init_wall(r_bnds(1),ints(1,1))



    do in = 1,rbnd_cnt -1 

       r1 = r_bnds(in)
       r2 = r_bnds(in+1)

       if (grid_option.eq.0.or.in.eq.1) then 
          do it = 1,int_cnt(in)
             vert_rec1(it) = ints(in,it)
             vert_type1(it) = int_type(in,it)
             !write(0,'(a,i8,10(1x,g18.8))') 'VERT1A:',it,vert_rec1(it),vert_type1(it)
             write(outunit,'(a,i8,10(1x,g18.8))') 'VERT1A:',it,vert_rec1(it),vert_type1(it)
          end do
          vert_cnt1 = int_cnt(in)
       end if


       do it = 1,int_cnt(in+1)
          vert_rec2(it) = ints(in+1,it)
          vert_type2(it) = int_type(in+1,it)
       end do
       vert_cnt2 = int_cnt(in+1)

       ntot1 = vert_cnt1
       ntot2 = vert_cnt2


       !
       ! For grid option 1 - add all tangency data between its end points to each ring as vertices for r values less than the tangency point location
       ! Also for the first row we need to define all vertices subject to the constraints of minimum number of cells between tangency positions and maximum allowed S separation
       !

       if (grid_option.eq.1) then

          if (in.eq.1) then 
             do it = 1,av_tan_cnt
                if (av_tan_r(it).gt.r1) then
                   ntot1 = ntot1+1
                   vert_rec1(ntot1) = av_tan_s(it)
                   vert_type1(ntot1) = FIXED_VERTEX

                   !write(0,'(a,2i8,10(1x,g18.8))') 'VERT1B:',ntot1,it,vert_rec1(ntot1),vert_type1(ntot1)
                   write(outunit,'(a,2i8,10(1x,g18.8))') 'VERT1B:',ntot1,it,vert_rec1(ntot1),vert_type1(ntot1)
                endif
             end do

             vert_cnt1 = ntot1
             call sort_arrays(1,vert_cnt1,vert_rec1,vert_type1)

             do is = 1,vert_cnt1-1
                ! place additional vertices between each tangency/surface using the grid generation contraints
                slen = vert_rec1(is+1) - vert_rec1(is)
                ncells = int(abs(slen)/max_s_sep) + 1

                if (ncells.lt.min_cells) ncells = min_cells

                !ssep1 = abs(slen)/ncells 
                ! remove abs so code will work with ascending or descending orders
                !ssep1 = slen/ncells 

                call gen_cell_spacing(ncells,ssep1a)

                do it = 1,ncells-1
                   !vert_rec1(ntot1+it) = vert_rec1(is)+ ssep1 * it
                   vert_rec1(ntot1+it) = vert_rec1(is)+ ssep1a(it) * slen
                   vert_type1(ntot1+it) = FIXED_VERTEX
                   !write(0,'(a,2i8,10(1x,g18.8))') 'VERT1C:',ntot1+it,is,vert_rec1(ntot1+it),vert_type1(ntot1+it)
                   write(outunit,'(a,2i8,10(1x,g18.8))') 'VERT1C:',ntot1+it,is,vert_rec1(ntot1+it),vert_type1(ntot1+it)
                end do
                ntot1 = ntot1 + ncells -1


             end do

             vert_cnt1 = ntot1
             call sort_arrays(1,vert_cnt1,vert_rec1,vert_type1)


          endif


          ! only add the tangency data to rec1 for the first row - otherwise we get duplication
          do it = 1,av_tan_cnt

             if (av_tan_r(it).gt.r2) then 
                ntot2 = ntot2+1
                vert_rec2(ntot2) = av_tan_s(it)
                vert_type2(ntot2) = FIXED_VERTEX
             endif

          end do

          vert_cnt1 = ntot1
          vert_cnt2 = ntot2

       end if


       call sort_arrays(1,vert_cnt1,vert_rec1,vert_type1)
       call sort_arrays(1,vert_cnt2,vert_rec2,vert_type2)


       !write(0,'(a,2i8,10(1x,g18.8))') 'VERT1:',vert_cnt1,ntot1
       !do it = 1,ntot1
       !   write(0,'(a,i8,10(1x,g18.8))') 'VERT1D:',in,r1,vert_rec1(it),vert_type1(it)
       !end do
       !write(0,'(a,2i8,10(1x,g18.8))') 'VERT2:',vert_cnt2,ntot2
       !do it = 1,ntot2
       !   write(0,'(a,i8,10(1x,g18.8))') 'VERT2D:',in,r2,vert_rec2(it),vert_type2(it)
       !end do


       write(outunit,'(a,2i8,10(1x,g18.8))') 'VERT1:',vert_cnt1,ntot1
       do it = 1,ntot1
          write(outunit,'(a,i8,10(1x,g18.8))') 'VERT1D:',in,r1,vert_rec1(it),vert_type1(it)
       end do
       write(outunit,'(a,2i8,10(1x,g18.8))') 'VERT2:',vert_cnt2,ntot2
       do it = 1,ntot2
          write(outunit,'(a,i8,10(1x,g18.8))') 'VERT2D:',in,r2,vert_rec2(it),vert_type2(it)
       end do

       !
       ! Break the calls to gen_ring based on number of paired surface elements
       !
       ! Copy new vertices to the end of vert_rec2 after call to gen_ring
       !
       ! Implicit assumption is that the first cell boundary is unbroken from one end to the other
       !

       npts1a = 0
       npts1b = 0
       npts2a = 0
       npts2b = 0

       rings_processed = 0

       ! assign start/end indices on each row
       do it = 1,vert_cnt1
          !write(0,'(a,7i8,g18.8)') 'VERT:',in,it,vert_cnt1,npts1a,npts1b,npts2a,npts2b,vert_type1(it)
          write(outunit,'(a,7i8,g18.8)') 'VERT:',in,it,vert_cnt1,npts1a,npts1b,npts2a,npts2b,vert_type1(it)

          ! modify code for grid generation from top down (i.e. max S to min S) ... base organization was reversed to start and end tags are reversed
          !if (vert_type1(it).eq.SURFACE_START) then 
          if (vert_type1(it).eq.SURFACE_END) then 
             npts1a = it
          elseif (vert_type1(it).eq.TANGENCY.or.vert_type1(it).eq.SURFACE_START) then 
             !elseif (vert_type1(it).eq.TANGENCY.or.vert_type1(it).eq.SURFACE_END) then 
             npts1b = it

             ! Now find the matching indices in vert_rec2
             s_start = vert_rec1(npts1a)
             s_end = vert_rec1(npts1b)
             rings_processed = rings_processed + 1


             write(outunit,'(a,6i8,10(1x,g18.8))') 'VERTB:',in,it,vert_cnt1,vert_cnt2,npts1a,npts1b,s_start,s_end

             ! loop through vert_rec2

             npts2a = 0
             npts2b = 0
             vert_set = .false.

             end_count = 0
             start_count = 0

             do is = 1,vert_cnt2
                ! s_start is numerically larger than s_end
                !if (vert_rec2(is).ge.s_start.and.vert_rec2(is).le.s_end) then 
                !if (vert_rec2(is).le.s_start.and.vert_rec2(is).ge.s_end) then 
                   ! this should assign npts2a and npts2b to appropriate values
                !   if (npts2a.eq.0) then 
                !      npts2a = is
                !   endif

                !   npts2b = is
                !   vert_set=.true.

                !elseif (vert_set) then 
                !   exit
                !endif
                !
                ! Need to change algorithm - want the pair of SURFACE_END/SURFACE_START values corresponding to rings_processed

                if (vert_type2(is).eq.SURFACE_END) then
                   end_count = end_count + 1
                   if (end_count.eq.rings_processed) then 
                      npts2a = is
                      write(outunit,'(a,2i10,10(1x,g18.8))') 'NPTS2A:',is,vert_cnt2,vert_rec2(is),s_start,s_end
                   endif
                endif

                if (vert_type2(is).eq.SURFACE_START) then
                   start_count = start_count + 1
                   if (start_count.eq.rings_processed) then 
                      npts2b = is
                      write(outunit,'(a,2i10,10(1x,g18.8))') 'NPTS2B:',is,vert_cnt2,vert_rec2(is)
                      exit
                   endif
                endif

             end do

             line_length = abs(s_start-s_end)
             rfactor = (r_limiter_max-r1)/(r_limiter_max-r_limiter_min)
             lfactor = line_length*rfactor


             if (npts2a.ne.0.and.npts2b.ne.0.and.(lfactor.ge.lcutoff)) then 

                ! record wall segments at each end of ring
                call add_wall_segment(r1,r2,vert_rec1(npts1a),vert_rec2(npts2a))
                call add_wall_segment(r1,r2,vert_rec1(npts1b),vert_rec2(npts2b))


                write(outunit,'(a,6i8,10(1x,g18.8))') 'GEN_RING: CALL: ', in,it,npts1a,npts1b,npts2a,npts2b,r1,r2,s_start,s_end
                call gen_ring(npts1a,npts1b,ntot1,r1,vert_rec1,vert_type1,npts2a,npts2b,ntot2,r2,vert_rec2,vert_type2,max_nknots,max_s_sep,min_cells)

             else

                if (npts1b.gt.npts1a) then 
                   ! Not generating a ring - add a wall segment 
                   call add_wall_segment(r1,r1,vert_rec1(npts1a),vert_rec1(npts1b))
                endif
                
                write(outunit,'(a,6i8,15(1x,g18.8))') 'GEN_RING: NO CALL: ', in,it,npts1a,npts1b,npts2a,npts2b,r1,r2,s_start,s_end,line_length,rfactor,lfactor,lcutoff

             endif
             

             ! If the point is a tangency then it is also the next surface start point
             if (vert_type1(it).eq.TANGENCY) then 
                npts1a = it
             endif

          endif

       end do


       ! if new vertices need to be retained then sort the vert_rec2 array and assign it to vert_rec1 
       ! for grid option 2 - retain vertices as fixed for next row

       if (grid_option.eq.1) then 
          vert_cnt2 = ntot2
          call sort_arrays(1,vert_cnt2,vert_rec2,vert_type2)
          ! assign "2" arrays to "1" arrays for the next row in grid

       endif

       ! move over vertices by one row ... either keeping ones added or not
       vert_cnt1 = vert_cnt2
       vert_rec1 = vert_rec2
       vert_type1 = vert_type2

    end do


    ! Deallocate local storage

    if (allocated(r_bnds)) deallocate(r_bnds)
    !if (allocated(s_bnds)) deallocate(s_bnds)
    if (allocated(int_type)) deallocate(int_type)
    if (allocated(int_cnt)) deallocate(int_cnt)
    if (allocated(ints)) deallocate(ints)
    if (allocated(vert_rec1)) deallocate(vert_rec1)
    if (allocated(vert_type1)) deallocate(vert_type1)
    if (allocated(vert_rec2)) deallocate(vert_rec2)
    if (allocated(vert_type2)) deallocate(vert_type2)


  end subroutine generate_grid


  subroutine init_wall(r_start,s_start)
    implicit none
    ! initialize data structures for collecting wall elements
    real*8 :: r_start,s_start
    integer :: ierr

    n_wall_segments = 0

    wall_start_r = r_start
    wall_start_s = s_start

    if (max_wall_segments.gt.0) then 
       if (allocated(wall_segments)) deallocate(wall_segments)
       allocate(wall_segments(max_wall_segments,5),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:WALL SEGMENTS:IERR =',ierr)
          stop
       endif
    else
       return
    endif

  end subroutine init_wall

  subroutine add_wall_segment(r1,r2,s1,s2) 
    implicit none
    real*8 :: r1,r2,s1,s2

    ! add a wall segment to the list
    n_wall_segments = n_wall_segments + 1

    if (n_wall_segments.gt.max_wall_segments) then 
       call grow_array2(wall_segments,max_wall_segments,max_wall_segments+50,5,5)
       max_wall_segments = max_wall_segments + 50
    endif

    wall_segments(n_wall_segments,1) = r1
    wall_segments(n_wall_segments,2) = s1
    wall_segments(n_wall_segments,3) = r2
    wall_segments(n_wall_segments,4) = s2
    ! flag to mark when the segement gets used
    wall_segments(n_wall_segments,5) = 0

  end subroutine add_wall_segment

  subroutine finalize_wall
    implicit none
    ! take all the accumulated wall segments - connect them together and connect ends to the av_wall data 

    logical :: done
    real*8 :: r_start,s_start,rp,sp,rn,sn
    integer :: rc,walli
    integer :: ierr,in


    ! Allocate storage for wall points

    if (n_wall_segments.gt.0) then 

       if (allocated(wall_r)) deallocate(wall_r)
       allocate(wall_r(n_wall_segments+1),stat=ierr)

       if (allocated(wall_s)) deallocate(wall_s)
       allocate(wall_s(n_wall_segments+1),stat=ierr)

       if (ierr.ne.0) then 
          call errmsg('ALLOCATION ERROR:WALL R,S ARRAYS:IERR =',ierr)
          stop
       endif
    else
       return
    endif

    write(0,'(a,i8,2(1x,g18.8))') 'WALL_SEGMENTS:',n_wall_segments,wall_start_r,wall_start_s
    do in = 1,n_wall_segments
       write(0,'(a,i8,5(1x,g18.8))') 'WS:',in,wall_segments(in,1),wall_segments(in,2),wall_segments(in,3),wall_segments(in,4),wall_segments(in,5)
    end do

    write(6,'(a,i8,2(1x,g18.8))') 'WALL_SEGMENTS:',n_wall_segments,wall_start_r,wall_start_s
    do in = 1,n_wall_segments
       write(6,'(a,i8,5(1x,g18.8))') 'WS:',in,wall_segments(in,1),wall_segments(in,2),wall_segments(in,3),wall_segments(in,4),wall_segments(in,5)
    end do

    done = .false. 

    r_start = wall_start_r
    s_start = wall_start_s

    nwall = 0

    nwall = nwall+1
    wall_r(nwall) = r_start
    wall_s(nwall) = s_start

    rp = r_start
    sp = s_start

    write(0,*) 'START:',rp,sp


    do while (.not.done)

       write(0,*) 'SEARCH:',rp,sp

       call get_next_point(rp,sp,rn,sn,walli,rc)

       write(0,'(a,3i8,6(1x,g18.8))') 'FIND:',nwall,walli,rc,rp,sp,rn,sn

       ! continuing on a set of connected segments
       if (rc.ne.0) then 
          nwall = nwall + 1
          wall_r(nwall) = rn
          wall_s(nwall) = sn
          ! mark this wall segment as used
          wall_segments(walli,5) = 1

          rp = rn
          sp = sn

       else
          ! wall segments have been added along grid edges when a ring is not generated so ... the wall should be contiguous and this should not happen. 
          ! This should only happen at the very end of the wall
          
          call check_wall_segments(ierr)
          
          if (ierr.ne.0) then 

             call errmsg('FINALIZE WALL: ALL WALL SEGMENTS NOT USED',ierr)
             stop 'Finalize_wall'

          endif

          ! Add a connection back to start point 
          nwall = nwall + 1
          wall_r(nwall) = r_start
          wall_s(nwall) = s_start

          done = .true.
       endif


    end do


  end subroutine finalize_wall

  subroutine check_wall_segments(ierr)
    implicit none
    integer :: ierr,in


    ! Loop through the wall segments and verify that all are marked as used - print any that are not

    ierr = 0

    do in = 1,n_wall_segments
       if (wall_segments(in,5).eq.0) then 
          ierr = ierr + 1
       endif 
    end do

  end subroutine check_wall_segments



  subroutine get_next_point(rp,sp,rn,sn,walli,rc)
    implicit none
    real*8 :: rp,sp,rn,sn
    integer :: walli, rc

    integer :: in

    rc = 0

    

    do in = 1,n_wall_segments

       if (wall_segments(in,5).ne.0) cycle

       if (wall_segments(in,1).eq.rp.and.wall_segments(in,2).eq.sp) then 
          rc = 1
          walli = in
          rn = wall_segments(in,3)
          sn = wall_segments(in,4)
       endif

       if (wall_segments(in,3).eq.rp.and.wall_segments(in,4).eq.sp) then 
          rc = 2
          walli = in
          rn = wall_segments(in,1)
          sn = wall_segments(in,2)
       endif


    end do

    return
  end subroutine get_next_point



  subroutine gen_ring(npts1a,npts1b,ntot1,r1,vert_rec1,vert_type1,npts2a,npts2b,ntot2,r2,vert_rec2,vert_type2,maxpts,max_slen,min_cells)
    implicit none
    integer :: npts1a,npts1b,npts2a,npts2b,min_cells,ntot1,ntot2,maxpts
    real*8 :: vert_rec1(maxpts),vert_type1(maxpts),vert_rec2(maxpts),vert_type2(maxpts)
    real*8 :: max_slen
    real*8 :: r1,r2
    real*8 :: s1(maxpts),s2(maxpts),s1_type(maxpts),s2_type(maxpts)

    !
    ! This routine takes in two rows of end points + any tangency or other predefined points then spaces out the corners of the 
    ! polygons to give a reasonable spacing between the two end surfaces
    !
    ! The routine will add to any exiting points without modifying them - if there are no exisitng points a quasi uniform mesh is generated
    !
    ! Code increments the ring count and records knot and polygon information for this ring
    !


    !  r1, sects    r2 , sects  ... end points are considered surfaces and can not be moved

    ! Cases: 

    ! Allocate point type storage


    integer :: npts1,npts2
    integer :: in,it,is,cells_added
    integer :: ncells, ncells2
    integer :: last_tan
    real*8 :: slen,ssep1,ssep2
    real*8 ::  slen1a,slen2a
    real*8,allocatable :: ssep1a(:),ssep2a(:)

    real*8 :: dist(maxpts),cell_dist(maxpts),tmp_dist(maxpts)

    integer :: val_loc(1)


    npts1 = npts1b-npts1a+1
    npts2 = npts2b-npts2a+1

    s1 = 0.0
    s2 = 0.0
    s1_type = 0.0
    s2_type = 0.0

    !
    ! Copy intersection data from vertex arrays into local arrays
    !
    do in = npts1a,npts1b
       s1(in-npts1a+1) = vert_rec1(in)
       s1_type(in-npts1a+1) = vert_type1(in)
       write(outunit,'(a,2i8,10(1x,g18.8))') 'S1:',in,in-npts1a+1,s1(in-npts1a+1),s1_type(in-npts1a+1)
       !write(6,'(a,2i8,10(1x,g18.8))') 'S1:',in,in-npts1a+1,s1(in-npts1a+1),s1_type(in-npts1a+1)
    end do

    do in = npts2a,npts2b
       s2(in-npts2a+1) = vert_rec2(in)
       s2_type(in-npts2a+1) = vert_type2(in)
       write(outunit,'(a,2i8,10(1x,g18.8))') 'S2:',in,in-npts2a+1,s2(in-npts2a+1),s2_type(in-npts2a+1)       
       !write(6,'(a,2i8,10(1x,g18.8))') 'S2:',in,in-npts2a+1,s2(in-npts2a+1),s2_type(in-npts2a+1)       
    end do


    !
    ! Generate polygons and additional vertices as required
    !

    if (npts1.eq.2.and.npts2.eq.2) then 
       ! determine number of vertices req'd

       ! Average slen
       !slen = max(abs(s1(1)-s1(2)),abs(s2(1)-s2(2)))

       slen = s1(2)-s1(1)

       ncells = int(abs(slen)/max_slen) + 1

       if (ncells.lt.min_cells) ncells = min_cells


       !dist1 = s1(2)-s1(1)
       !dist2 = s2(2)-s2(1)

       !call gen_cell_spacing(dist1,dist2,ncells-1,ssep1,ssep2)
       call gen_cell_spacing(ncells,ssep1a)

       slen1a = s1(2)-s1(1)
       slen2a = s2(2)-s2(1)

       !ssep1 = (s1(2)-s1(1))/ncells 
       !ssep2 = (s2(2)-s2(1))/ncells       

       do in = 1,ncells-1
          !s1(npts1+in) = s1(1) + ssep1 * in
          s1(npts1+in) = s1(1) + ssep1a(in) * slen1a
          s1_type(npts1+in) = NEW_VERTEX
          !s2(npts1+in) = s2(1) + ssep2 * in
          s2(npts1+in) = s2(1) + ssep1a(in) * slen2a
          s2_type(npts1+in) = NEW_VERTEX
       end do

       npts1 = npts1 + ncells -1
       npts2 = npts2 + ncells -1

       call sort_arrays(1,npts1,s1,s1_type)
       call sort_arrays(1,npts2,s2,s2_type)

    elseif (npts1.eq.2.and.npts2.gt.2) then 
       ! Case with one or more tangency points on second boundary

       slen = s1(2)-s1(1)

       ncells = int(abs(slen)/max_slen) + 1

       if (ncells.lt.min_cells) ncells = min_cells

       !ssep1 = (s1(2)-s1(1))/ncells 
       slen1a = s1(2)-s1(1)
       call gen_cell_spacing(ncells,ssep1a)


       !write(outunit,*) 'GR:NPTS:',npts1,npts2,ncells,slen

       do in = 1,ncells-1
          !s1(npts1+in) = s1(1) + ssep1 * in
          s1(npts1+in) = s1(1) + ssep1a(in) * slen1a
          s1_type(npts1+in) = NEW_VERTEX
       end do

       ncells2 = ncells - npts2 +2

       slen = (s2(npts2)-s2(1))

       !write(outunit,'(a,3i8,10(1x,g18.8))') 'GR:NPTS2:',npts2,ncells,ncells2,slen

       dist = 0.0
       do in = 1,npts2-1
          dist(in) = abs(s2(in+1)-s2(in))/abs(slen) * ncells2
          !write(outunit,'(a,i8,10(1x,g18.8))') 'GR:DIST:',in,dist(in),s2(in+1),s2(in)
       end do

       ! assign cell dist an initial distribution of one additional vertex / group ... so that division by zero does not occur in cell assignment algorithm
       cell_dist = 1.0
       tmp_dist = 0.0
       do in = 1,ncells2-1
          tmp_dist = dist/cell_dist
          val_loc = maxloc(tmp_dist)
          it = val_loc(1)
          !it = maxloc(dist)
          cell_dist(it) = cell_dist(it)+1.0
          !dist(it) = dist(it)/(cell_dist(it)+1.0)
          !WRITE(outunit,'(a,2i8,10(1x,g18.8))') 'dist2:',in,it,cell_dist(it),dist(it)
       end do

       !write(0,*) 'DIST sum:',sum(cell_dist),ncells2-1

       cells_added = 0
       do in = 1,npts2-1
          ! cell_dist is assigned a base value of 1.0 above so an extra 1.0 is not required here - it must already be non-zero
          !ssep2 = (s2(in+1)-s2(in))/(cell_dist(in)+1.0)
          slen2a = s2(in+1)-s2(in)
          call gen_cell_spacing(int(cell_dist(in)),ssep2a)
          !ssep2 = (s2(in+1)-s2(in))/(cell_dist(in))
          do it = 1,cell_dist(in)-1
             cells_added = cells_added + 1
             !s2(npts2+cells_added) = s2(in) + ssep2 * it
             s2(npts2+cells_added) = s2(in) + ssep2a(it) * slen2a
             s2_type(npts2+cells_added) = NEW_VERTEX
             !write(outunit,'(a,3i8,10(1x,g18.8))') 'Adding:',in,it,cells_added,ssep2,s2(npts2+cells_added),s2_type(npts2+cells_added)
          end do
       end do

       npts1 = npts1 + ncells -1
       npts2 = npts2 + cells_added

       call sort_arrays(1,npts1,s1,s1_type)
       call sort_arrays(1,npts2,s2,s2_type)

    elseif (npts1.gt.2.and.npts2.eq.2) then 
       ! points pre-specified along r1 ... none defined on r2 ... do proportional spacing to match

       slen = (s1(npts1)-s1(1))

       do in = 1,npts1
          dist(in) =   abs((s1(in)-s1(1)))/abs(slen)
       end do

       ! add points with same spacing along r2

       slen = (s2(npts2)-s2(1))

       do in = 2,npts1-1
          s2(npts2+in-1) = s2(1) + dist(in) * slen
          s2_type(npts2+in-1) = NEW_VERTEX
       end do

       npts2 = npts2 + (npts1-2)

       call sort_arrays(1,npts1,s1,s1_type)
       call sort_arrays(1,npts2,s2,s2_type)

    elseif (npts1.gt.2.and.npts2.gt.2) then 
       ! pre defined corner points along R1 PLUS tangency points on R2
       ! In this case, the code expects exact correspondence for the tangency points in R so the problem devolves into proportional 
       ! distribution between tangency points


       cells_added = 0
       last_tan = 1
       do in = 2,npts1-1

          do it = 2,npts2-1

             ! found match for tangency points
             if (s1(in).eq.s2(it)) then 
                do is = last_tan + 1,in-1
                   cells_added = cells_added + 1
                   s2(npts2+cells_added) = abs(s1(is)-s1(last_tan))/abs(s1(in)-s1(last_tan)) * (s2(it)-s2(it-1)) + s2(it-1)
                   s2_type(npts2+cells_added) = NEW_VERTEX
                   !write(outunit,'(a,4i8,10(1x,g18.8))') 'NEW S2:',in,it,is,last_tan,cells_added,s1(is),s1(in),s1(last_tan),s2(it),s2(it-1),s2(npts2+cells_added)
                end do
                last_tan= in
                ! if we have just finished the match for the last tangency point then process the end of the row
                if (it.eq.npts2-1) then 

                   do is = in+1,npts1-1
                      cells_added = cells_added + 1
                      s2(npts2+cells_added) = abs(s1(is)-s1(last_tan))/abs(s1(npts1)-s1(last_tan)) * (s2(npts2)-s2(it)) + s2(it)
                      s2_type(npts2+cells_added) = NEW_VERTEX
                      !write(outunit,'(a,4i8,10(1x,g18.8))') 'NEW S2:',in,it,is,last_tan,cells_added,s1(is),s1(in),s1(last_tan),s2(it),s2(it-1),s2(npts2+cells_added)
                   end do

                endif
                exit
             end if

          end do

       end do

       npts2 = npts2 + cells_added 
       call sort_arrays(1,npts2,s2,s2_type)

       write(outunit,'(a,i8,10(1x,g18.8))') 'NEW s2:',npts2,cells_added

    endif


    ! for grid generation number of points on each side should be the same - matching vertices ... if not an error has occurred
    if (npts1.eq.npts2) then 

       ! order of listing vertices is important as is ordering cells
       ! the base coordinate system for the ITER data is from small R to larger R >0 
       ! and from S < 0 to S > 0 ... the cells are organized wiith DOWN at +S between vertices 1,2 and UP at -S between 3,4
       ! ... so a cell is   1 = s1(in)  2 = s2(in)   3 = s2(in+1)  4 = s1(in+1)
       !
       ! Lets create the ring and cells

       nrings = nrings + 1
       nknots(nrings) = npts1 -1 

       write(outunit,'(a,4i8,10(1x,g18.8))') 'GEN:',npts1,npts2,nrings,nknots(nrings)

       do in = 1,npts1 -1 
          npoly = npoly + 1
          nvp(npoly) = 4
          rvp(1,npoly) = r1
          zvp(1,npoly) = s1(in)
          rvp(2,npoly) = r2
          zvp(2,npoly) = s2(in)
          rvp(3,npoly) = r2
          zvp(3,npoly) = s2(in+1)
          rvp(4,npoly) = r1
          zvp(4,npoly) = s1(in+1)

          poly_ref(in,nrings) = npoly
          rcen(in,nrings) = sum(rvp(1:nvp(npoly),npoly))/nvp(npoly)
          zcen(in,nrings) = sum(zvp(1:nvp(npoly),npoly))/nvp(npoly)
          write(outunit,'(a,3i8,10(1x,g18.8))') 'POLY:',in,npoly,poly_ref(in,nrings),(rvp(it,npoly),zvp(it,npoly),it=1,nvp(npoly)),rcen(in,nrings),zcen(in,nrings)

       end do

    else

       call errmsg('GEN_RING','ERROR in vertex generation')
       write(outunit,'(a,2i10,10(1x,g18.8))') 'ERROR in vertex generation',npts1,npts2
       stop

    endif



    ! copy new vertices and types into vert_rec2 - vert_rec1 was completed on previous iterations

    !write(outunit,'(a,3i8,10(1x,g18.8))') 'Assign vert_rec2:', maxpts,npts2,ntot2

    do in = 1,npts2
       if (s2_type(in).eq.NEW_VERTEX) then 
          ntot2 = ntot2 + 1
          !write(outunit,'(a,3i8,10(1x,g18.8))') 'REC VERT:', in,ntot2,maxpts,s2(in),s2_type(in)

          vert_rec2(ntot2) = s2(in)
          vert_type2(ntot2) = FIXED_VERTEX
       endif
    end do

    ! save resorting vert_rec2 and vert_type2 until after polygon generation for these rings is complete


  end subroutine gen_ring


  subroutine write_grid
    implicit none

    integer :: in,ir,ik
    integer :: iunit
    ! write out the grids polygons
    character*100 :: filename


    ! find free unit number and open the file
    !call find_free_unit_number(iunit)
    !filename = 'grid.out'
    !open(iunit,file=filename,form='formatted')

    iunit = 6

    !write(0,'(a,3i8,10(1x,g18.8))') 'Grid:', npoly,nrings
    write(iunit,'(a,3i8,10(1x,g18.8))') 'Grid:', npoly,nrings

    do in = 1,nrings
       write (iunit,'(a,3i8)') 'Knots:',in,nknots(in)
    enddo

    do ir = 1,nrings

       do ik = 1,nknots(ir)

          write (iunit,'(a,3i8,12(1x,g18.8))') 'Polys:',ik,ir,poly_ref(ik,ir),(rvp(in,poly_ref(ik,ir)),zvp(in,poly_ref(ik,ir)),in=1,nvp(poly_ref(ik,ir))),rcen(ik,ir),zcen(ik,ir)
       end do

    enddo

    close(iunit)


    !filename = 'poly.out'
    !open(iunit,file=filename,form='formatted')

    write(iunit,'(a,3i8,10(1x,g18.8))') 'Polygon vertices:', npoly,nrings

    do ir = 1,nrings

       do ik = 1,nknots(ir)

          do in = 1,nvp(poly_ref(ik,ir))

             write (iunit,'(2(1x,g18.8))') rvp(in,poly_ref(ik,ir)),zvp(in,poly_ref(ik,ir))

          enddo

          write(iunit,'(a)') 

       end do

    enddo

    close(iunit)




  end subroutine write_grid






  subroutine find_field_line_intsects(max_intsects,r_line,int_cnt,int_type,intersects)
    implicit none
    integer :: max_intsects
    real*8 :: r_line,int_cnt,int_type(max_intsects),intersects(max_intsects)

    ! Loop through the limiter surface finding and counting all intersections - record the S value of the intersection - define the type
    ! exclude intersections at tangency points 
    integer :: sect_type

    real*8 :: r_start,s_start,r_end,z_end,s_end,s_int
    real*8 :: r_int,z_int
    real*8 :: int_type_base
    real*8 :: deleted
    real*8 :: last_surface

    real*8, parameter :: eps = 1.0d-6

    integer :: in,it

    int_cnt = 0.0
    int_type_base = SURFACE_START


    do in = 1,av_group_cnt -1

       r_start = av_r(in)
       s_start = av_s(in)
       r_end = av_r(in+1)
       s_end = av_s(in+1)

       call intsect2dp(r_line,s_limiter_min,r_line,s_limiter_max,r_start,s_start,r_end,s_end,r_int,s_int,sect_type)


       if (sect_type.eq.1) then 
          if (.not.(r_int.eq.r_end.and.av_type(in+1).eq.WALL.or.r_int.eq.r_start.and.av_type(in).eq.WALL)) then 
             ! need to do something about any tangency point intersections found - removed later
             ! need to make sure wall point intersections are not included since it results in a zero length polygon side
             ! an intersection point type of 1 indicates that the intersection point lies on both line segments - otherwise ignore the result
             int_cnt = int_cnt + 1.0

             ! use int_type to deterimine the direction from an intersection point in order to determine later whether a degenerate intersection on the wall
             ! or a degenerate intersection at a tangency or wall point has been found. 

             if (r_start-r_end.lt.0.0) then
                int_type(int(int_cnt)) = SURFACE_START
             else
                int_type(int(int_cnt)) = SURFACE_END
             endif

             !if (int_type_base.gt.0.0) then 
             !   int_type(int(int_cnt)) = SURFACE_START
             !else
             !   int_type(int(int_cnt)) = SURFACE_END
             !endif
             !int_type(int(int_cnt)) = int_type_base
             !int_type_base = int_type_base * -1

             intersects(int(int_cnt)) = s_int

             write(outunit,'(a,2i8,12(1x,g18.8))') 'Int    :',in,sect_type,int_cnt,r_int,s_int,int_type_base,r_start,s_start,r_end,s_end,r_int,s_int
             !write(0,'(a,2i8,12(1x,g18.8))') 'Int    :',in,sect_type,int_cnt,r_int,s_int,int_type_base,r_start,s_start,r_end,s_end,r_int,s_int

          endif

       endif

    end do

    ! Add tangency points if any


    do in = 1,av_tan_cnt
       if (av_tan_r(in).eq.r_line) then 
          int_cnt = int_cnt + 1.0

          write(outunit,'(a,2i8,10(1x,g18.8))') 'Int Tan:',in,sect_type,int_cnt,av_tan_r(in),av_tan_s(in)
          !write(0,'(a,2i8,10(1x,g18.8))') 'Int Tan:',in,sect_type,int_cnt,av_tan_r(in),av_tan_s(in)

          int_type(int(int_cnt)) = TANGENCY
          intersects(int(int_cnt)) = av_tan_s(in)
          !write(outunit,'(a,2i8,10(1x,g18.8))') 'Int Tan:',in,sect_type,int_cnt,av_tan_r(in),av_tan_s(in)
       end if
    end do


    ! sort 

    call sort_arrays(0,int(int_cnt),intersects,int_type)

    ! filter out any duplicates (points with S values closer than the minimum allowed - these should be intersections at tangency points ... if given a choice keep
    ! the one identified as a tangency point.


    write(outunit,'(a,10(1x,g18.8))') 'prelim surface intersections output:', int_cnt, r_line
    do in = 1,int_cnt
       write(outunit,'(a,i8,10(1x,g18.8))') 'PINT sect:', in,int_type(in),intersects(in)
    end do


    deleted = 0

    do in = 1,int_cnt-1

       if (int_type(in).ne.DELETE_POINT) then
          do it = in+1,int_cnt
             if (int_type(it).ne.DELETE_POINT) then 

                if (abs(intersects(it)-intersects(in)).lt.eps) then

                   ! consecutive intersections too close together
                   if (int_type(it).eq.TANGENCY.and.int_type(in).ne.TANGENCY) then 
                      int_type(in) = DELETE_POINT
                      deleted = deleted + 1.0
                      write(outunit,'(a,2i8,l4,10(1x,g18.8))') 'Delete1a:',in,it, (abs(intersects(it)-intersects(in)).lt.eps),intersects(it),intersects(in),eps,int_type(in),int_type(it)
                   elseif (int_type(in).eq.TANGENCY.and.int_type(it).ne.TANGENCY) then 
                      int_type(it) = DELETE_POINT
                      deleted = deleted + 1.0
                      write(outunit,'(a,2i8,l4,10(1x,g18.8))') 'Delete1b:',in,it, (abs(intersects(it)-intersects(in)).lt.eps),intersects(it),intersects(in),eps,int_type(in),int_type(it)
                   elseif (int_type(in).eq.int_type(it)) then 
                      ! if the intersection type of the two points is the same .. only delete one of the points
                      int_type(it) = DELETE_POINT
                      deleted = deleted + 1.0
                      write(outunit,'(a,2i8,l4,10(1x,g18.8))') 'Delete1c:',in,it, (abs(intersects(it)-intersects(in)).lt.eps),intersects(it),intersects(in),eps,int_type(in),int_type(it)
                   elseif (int_type(in).ne.int_type(it)) then 
                      ! two points overlapping - should be wall turning point since int_types are different at this stage and tangency dealt with above - delete both
                      int_type(in) = DELETE_POINT
                      int_type(it) = DELETE_POINT
                      deleted = deleted + 2.0
                      write(outunit,'(a,2i8,l4,10(1x,g18.8))') 'Delete2 :',in,it, (abs(intersects(it)-intersects(in)).lt.eps),intersects(it),intersects(in),eps,int_type(in),int_type(it)
                   endif

                   if (int_type(in).eq.DELETE_POINT) exit

                endif
             endif

          end do

       endif
    end do

    ! Filter out the 10.0s

    call sort_arrays(0,int(int_cnt),int_type,intersects)

    int_cnt = int_cnt - deleted

    call sort_arrays(0,int(int_cnt),intersects,int_type)

    ! reset the types appropriately

    int_type(1) = SURFACE_START
    last_surface = SURFACE_START
    do in = 2,int_cnt
       if (int_type(in).ne.TANGENCY) then 
          if (last_surface.eq.SURFACE_START) then 
             int_type(in) = SURFACE_END
          elseif (last_surface.eq.SURFACE_END) then 
             int_type(in) = SURFACE_START
          endif
          last_surface = int_type(in)
       endif
    end do

    if (last_surface.ne.SURFACE_END) then 
       write(outunit,'(a,10(1x,g18.8))') 'surface intersections error:', int_cnt, r_line
       do in = 1,int_cnt
          write(outunit,'(a,i8,10(1x,g18.8))') ' INT sect:', in,int_type(in),intersects(in)
       end do

       !write(0,'(a,10(1x,g18.8))') 'surface intersections error:', int_cnt, r_line
       !do in = 1,int_cnt
       !   write(0,'(a,i8,10(1x,g18.8))') ' INT sect:', in,int_type(in),intersects(in)
       !end do
       call errmsg('SURFACE INTERSECTION ERROR: Last intersection is not a surface on R=',r_line)
       stop
    else   

       write(outunit,'(a,10(1x,g18.8))') 'surface intersections output:', int_cnt, r_line
       do in = 1,int_cnt
          write(outunit,'(a,i8,10(1x,g18.8))') ' INT sect:', in,int_type(in),intersects(in)
       end do

    endif



  end subroutine find_field_line_intsects

  subroutine assign_grid_to_divimp(maxnrs,maxnks,mves,nrs,nks,&
       nves,rves,zves,&
       npolyp,korpg,&
       nvertp,rvertp,zvertp,&
       rs,zs)
    implicit none
    integer :: ierr
    integer :: maxnrs,maxnks,mves

    integer :: nrs,nks(maxnrs),nves,npolyp
    integer :: korpg(maxnks,maxnrs)
    real :: rs(maxnks,maxnrs),zs(maxnks,maxnrs)
    real :: rves(mves),zves(mves)
    integer :: nvertp(maxnrs*maxnks)
    real :: rvertp(5,maxnrs*maxnks),zvertp(5,maxnrs*maxnks)


    ! local declarations

    integer :: in,ik,ir,it,is
    ! Assign grid polygons


    !write(0,*) 'Assign to DIVIMP:'


    ierr = 0

    if (npoly.le.maxnrs*maxnks) then 
       npolyp = npoly

       do in = 1,npoly
          nvertp(in) = nvp(in)
          do it = 1,4
             rvertp(it,in) = rvp(it,in)
             zvertp(it,in) = zvp(it,in)
          end do
       end do
    else
       ierr = 1
       call errmsg('Too many polygons in ribbon grid:',npoly)
    endif


    !write(0,*) 'Assign NRS:',nrs,nrings,maxnrs

    if (nrings.le.maxnrs) then 
       nrs = nrings

       do ir = 1,nrings
          if (nknots(ir).le.maxnks) then 
             nks(ir) = nknots(ir)
             do ik = 1,nknots(ir)
                rs(ik,ir) = rcen(ik,ir)
                zs(ik,ir) = zcen(ik,ir)
                korpg(ik,ir) = poly_ref(ik,ir)
                write(outunit,'(a,3i8,10(1x,g18.8))') 'POLY COPY:',ik,ir,korpg(ik,ir),rs(ik,ir),zs(ik,ir)
             end do
          else
             ierr = 1
             call errmsg('Too many knots on ribbon grid ring#:',ir)
             call errmsg('                      Knot count = :',nknots(ir))
          endif
       end do
    else
       ierr = 1
       call errmsg('Too many rings in ribbon grid:',nrings)
    endif

    if (ierr.ne.0) then 
       call errmsg('Ribbon grid does not fit into DIVIMP storage - increase MAXNRS and/or MAXNKS counts in params : Exiting')
       stop 'Exit in assign_grid_to_divimp'
    endif

    ! calculate rves,zves by running along the cell ends and the limiter surface line. 
    ! one source of these numbers could be the intersects array from field_line_intersections since these
    ! points are used to define the polygon surfaces at the ends of the rings. 

    ! Do I need a connection map to build this or is there an easier way?

    ! for now just assign the averaged surface since all intersection points lie on this surface
    ! ... may need to rebuild later from connection map data

    ! End points should already be included 

    ! This doesn't work because there are too many data points ... need to use the polygon grid ends or at least the end vertices combined
    ! with the av_wall points to define a wall. 

    !nves = av_group_cnt
    !
    !if (nves.le.mves) then 
    !   do in = av_group_cnt,1,-1
    !      rves(in) = av_r(in)
    !      zves(in) = av_s(in)
    !      write(6,'(a,i8,5(1x,g18.8))') 'Vessel:',in,av_r(in),av_s(in)
    !   end do

       ! Add 

    !else
    !   call errmsg('Too many elements in ribbon grid wall specification:',av_group_cnt)
    !   stop "Too many points in ribbon grid wall"
    !endif


    !
    ! Assign the new vessel wall calculated from connecting grid edges
    !

    call finalize_wall


    nves = n_wall_segments

    if (nves.gt.0.and.nves.le.mves) then 
       do in = 1,n_wall_segments
          rves(in) = wall_r(in)
          zves(in) = wall_s(in)
          write(6,'(a,i8,5(1x,g18.8))') 'Vessel:',in,wall_r(in),wall_s(in)
       end do

    else
       call errmsg('Invalid number of elements in ribbon grid wall specification:',n_wall_segments)
       stop "Too many points in ribbon grid wall"
    endif




  end subroutine assign_grid_to_divimp



  subroutine deallocate_castem_storage
    implicit none


    if (allocated(surf_r)) deallocate(surf_r)
    if (allocated(surf_s)) deallocate(surf_s)
    if (allocated(surf_fl)) deallocate(surf_fl)
    if (allocated(surf_int)) deallocate(surf_int)
    if (allocated(surf_sep)) deallocate(surf_sep)

    if (allocated(av_s)) deallocate(av_s)
    if (allocated(av_r)) deallocate(av_r)
    if (allocated(av_type)) deallocate(av_type)
    if (allocated(av_min_r)) deallocate(av_min_r)
    if (allocated(av_max_r)) deallocate(av_max_r)

    if (allocated(av_tan_r)) deallocate(av_tan_r)
    if (allocated(av_tan_s)) deallocate(av_tan_s)
    if (allocated(av_wall_r)) deallocate(av_wall_r)
    if (allocated(av_wall_s)) deallocate(av_wall_s)
    if (allocated(tan_ord_r)) deallocate(tan_ord_r)
    if (allocated(av_tan_ind)) deallocate(av_tan_ind)

    if (allocated(nvp)) deallocate(nvp)
    if (allocated(rvp)) deallocate(rvp)
    if (allocated(zvp)) deallocate(zvp)


    if (allocated(rcen)) deallocate(rcen)
    if (allocated(zcen)) deallocate(zcen)

    if (allocated(nknots)) deallocate(nknots)
    if (allocated(poly_ref)) deallocate(poly_ref)


  end subroutine deallocate_castem_storage

  subroutine grow_array(foo,num1,num2)
    implicit none
    real*8,allocatable :: foo(:)
    integer :: num1,num2,ierr,in

    real*8,allocatable :: tmp_foo(:)

    ierr = 0

    if (allocated(foo)) then 

       allocate(tmp_foo(num2),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW ALLOCATION ERROR:TMP_FOO:IERR =',ierr)
          stop
       endif

       tmp_foo = 0.0

       do in = 1,num1
          tmp_foo(in) = foo(in)
       end do

       deallocate(foo,stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW DE-ALLOCATION ERROR:FOO:IERR =',ierr)
          stop
       endif

       allocate(foo(num2),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW RE-ALLOCATION ERROR:FOO:IERR =',ierr)
          stop
       endif

       foo = tmp_foo

       deallocate(tmp_foo,stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW DE-ALLOCATION ERROR:TMP_FOO:IERR =',ierr)
          stop
       endif

    end if

  end subroutine grow_array

  subroutine grow_array2(foo,num1,num2,numa,numb)
    implicit none
    real*8,allocatable :: foo(:,:)
    integer :: num1,num2,numa,numb,ierr,in,is

    real*8,allocatable :: tmp_foo(:,:)

    ierr = 0

    if (allocated(foo)) then 

       if (allocated(tmp_foo)) deallocate(tmp_foo)
       allocate(tmp_foo(num2,numb),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW2 ALLOCATION ERROR:TMP_FOO:IERR =',ierr)
          stop
       endif

       tmp_foo = 0.0

       do in = 1,num1
          do is = 1,numa
             tmp_foo(in,is) = foo(in,is)
          end do 
       end do

       deallocate(foo,stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW2 DE-ALLOCATION ERROR:FOO:IERR =',ierr)
          stop
       endif

       allocate(foo(num2,numb),stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW2 RE-ALLOCATION ERROR:FOO:IERR =',ierr)
          stop
       endif

       foo = tmp_foo

       deallocate(tmp_foo,stat=ierr)
       if (ierr.ne.0) then 
          call errmsg('GROW2 DE-ALLOCATION ERROR:TMP_FOO:IERR =',ierr)
          stop
       endif

    end if

  end subroutine grow_array2



  subroutine gen_cell_spacing(ncells,ssep)
    implicit none
    integer :: ncells
    integer :: m
    real*8,allocatable :: ssep(:)
    real*8 :: x,dist

    ! This routine generates ncells numbers distributed between 0 and 1
    ! (leaving out the end points)
    ! If the cell_spacing_factor factor is 1.0 these numbers are linearly distributed. Increasing
    ! the exponent enhances grid resolution near the ends of the regions. 


    if (allocated(ssep)) deallocate(ssep)

    allocate(ssep(ncells-1))

    ! cell spacing option 0 - exponential - only supported option so far
    if (cell_spacing_option.eq.0) then 

       do m = 1,ncells-1
          if (m.le.(ncells/2)) then 
             x = real(m)/real(ncells)
             ssep(m) =  x**cell_spacing_factor
          else
             x = real(ncells-m)/real(ncells)

             ssep(m) = 1.0 - x**cell_spacing_factor
          endif

          !write(6,'(a,1x,2i8,10(1x,g18.8))') 'RES:',ncells,m,x,ssep(m)

       end do

    endif


  end subroutine gen_cell_spacing




end module castem_field_line_data

