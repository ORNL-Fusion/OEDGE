module castem_field_line_data

  use error_handling
  !use utilities
  use common_utilities

  implicit none

  private

  save



  public :: read_identifier_data, read_intersection_data,print_field_line_summary,calculate_limiter_surface,generate_grid,write_grid,&
            assign_grid_to_divimp,deallocate_castem_storage



!  type node_struc
!     integer :: fl
!     integer :: in
!     real*8  :: s,r
!  end type node_struc

  type intersection_struc
     integer :: line_id
     integer :: intsect_id
     real*8  :: xi,yi,zi,lc
     integer :: int_type
     logical :: int_used 
!     type (node_struc) :: next
!     type (node_struc) :: last
  end type intersection_struc

  type field_line_struc
     integer :: line_id
     real*8 :: xs,ys,zs
     integer :: int_up, int_down, int_tot
     integer :: int_stored
     type (intersection_struc), allocatable :: int_data(:)
  end type field_line_struc


  type (field_line_struc),allocatable :: field_line(:)

  integer :: n_tangency, n_enter, n_leave, n_field_lines, tot_n_intsects
  real*8 :: max_field_line_len,min_field_line_len
  real*8 :: max_lc,min_lc
  real*8 :: min_r, max_r

  character*256 :: date_castem,date_process

  logical :: header_has_been_loaded

!  type (node_struc) :: init_node

!  type (node_struc),allocatable :: node_list

  integer :: n_nodes
  real*8,allocatable :: surf_r(:),surf_s(:),surf_fl(:),surf_int(:),surf_sep(:)
  real*8,allocatable :: av_s(:),av_r(:),av_type(:),av_min_r(:),av_max_r(:)
  real*8,allocatable :: av_tan_r(:),av_tan_s(:),av_wall_r(:),av_wall_s(:),tan_ord_r(:),av_tan_ind(:)

  real*8 :: r_limiter_max,r_limiter_min,s_limiter_max,s_limiter_min,min_tan_sep

  integer :: av_group_cnt
  integer :: av_wall_cnt, av_tan_cnt

  integer :: intsec_cnt

! options
 
  logical :: opt_block_av= .true.



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


   integer,parameter :: grid_option = 1

contains


  subroutine read_identifier_data(file,ierr)
    implicit none
    character*(*) :: file
    integer :: ierr

    integer :: iunit

    ! local variables
    integer :: id,inter_up,inter_down,inter_tot
    real*8 :: xinit,yinit,zinit
    integer :: in

!    init_node%fl = -1
!    init_node%in = -1

    ! initialize min and max r coordinate for grid
    min_r =  1e24
    max_r = -1e24

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
       allocate(field_line(n_field_lines),stat=ierr)
       if (ierr.ne.0) return
    else
       return
    endif

    tot_n_intsects = n_enter + n_leave + n_tangency
    intsec_cnt = 0

    write(0,'(a,5(1x,i8))') 'Intsects:',tot_n_intsects,n_enter,n_leave,n_tangency


    ! field line array has been allocated ... now read in and allocate the intersection
    ! arrays

    do in = 1,n_field_lines

       read(iunit,*) id,xinit,yinit,zinit,inter_up,inter_down,inter_tot

       call assign_field_line_data(field_line(in),id,xinit,yinit,zinit,inter_up,inter_down,inter_tot,ierr)

    end do


  end subroutine read_identifier_data


  subroutine assign_field_line_data(fl,id,xinit,yinit,zinit,inter_up,inter_down,inter_tot,ierr)
    implicit none
    type(field_line_struc) :: fl
    integer :: id, inter_up,inter_down,inter_tot,ierr
    real*8 :: xinit,yinit,zinit

    min_r = min(min_r,xinit)
    max_r = max(max_r,xinit)

    fl%line_id = id
    fl%xs = xinit
    fl%ys = yinit
    fl%zs = zinit

    fl%int_up = inter_up
    fl%int_down = inter_down 
    fl%int_tot = inter_tot

    fl%int_stored = 0

    if (inter_tot.ne.(inter_up+inter_down)) then 
       call errmsg('Error: Intersections do not add up:',id)
    end if

    if (allocated(fl%int_data)) then 
       deallocate (fl%int_data)
    endif

    ! Allocate space to hold the total number of recorded intersections for this line
    allocate(fl%int_data(inter_tot),stat=ierr)

    write(6,'(a,3(1x,i8))') 'Assign FL:',id,inter_tot,ierr

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

       write(6,'(a)') trim(line)

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




  subroutine read_intersection_data(file,ierr)
    implicit none
    character*(*) :: file
    integer :: ierr


    integer :: iunit
    character*512 :: line

    integer :: line_id,inter_id,itype,iend
    real*8 :: xint,yint,zint,lc
    
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
       write(6,'(a)') trim(line)

       if (iend.eq.0) then 
          
          read(line,*) line_id,inter_id,xint,yint,zint,lc,itype

          call assign_intsect_data(line_id,field_line(line_id),inter_id,xint,yint,zint,lc,itype,ierr)
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

    write(6,'(a,i10)') 'Total intersections read:', sect_cnt


  end subroutine read_intersection_data

  subroutine assign_intsect_data(line_id,fl,inter_id,xint,yint,zint,lc,itype,ierr)
    implicit none
    integer :: line_id
    type(field_line_struc) :: fl
    integer :: inter_id,itype,ierr
    real*8 :: xint,yint,zint,lc
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
       fl%int_data(in)%int_type = itype
       fl%int_data(in)%int_used = .false. 

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

    fname = 'lim.txt'

    ! find free unit number and open the file
    call find_free_unit_number(outunit)

    open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(0,*) 'Writing out limiter surface: ' ,n_nodes

    do in = 1,n_nodes
       write(outunit,'(i8,5(1x,g18.8))') in,surf_r(in),surf_s(in),surf_fl(in),surf_int(in)
    end do

    close(outunit)


    fname = 'lim-av.txt'

    ! find free unit number and open the file
    call find_free_unit_number(outunit)

    open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(0,*) 'Writing out averaged limiter surface: ' ,av_group_cnt

    do in = 1,av_group_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,av_r(in),av_s(in)
    end do

    close(outunit)


    fname = 'tan-av.txt'

    ! find free unit number and open the file
    call find_free_unit_number(outunit)

    open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(0,*) 'Writing out tangency points: ' ,av_tan_cnt

    do in = 1,av_tan_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,av_tan_r(in),av_tan_s(in)
    end do

    close(outunit)

    fname = 'wall-av.txt'

    ! find free unit number and open the file
    call find_free_unit_number(outunit)

    open(outunit,file=trim(fname),form='formatted',iostat=ierr)

    write(0,*) 'Writing out wall points: ' ,av_wall_cnt

    do in = 1,av_wall_cnt
       write(outunit,'(i8,5(1x,g18.8))') in,av_wall_r(in),av_wall_s(in)
    end do

    close(outunit)








  end subroutine print_field_line_summary



  subroutine calculate_limiter_surface
    !use common_utilities
    implicit none

    integer :: in,if,ierr
    integer :: node_cnt
    real,allocatable :: surf_sep(:),av_group(:)
    !real,allocatable :: surf_ang(:)
    integer :: tmp_node_cnt 

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
    ! deivative change in the surface definition forming the bottom side of the tile as opposed to the peaks at limiter tips. Three methods can be used. 
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
    ! Start with an initial point at min_r, min_lc and end at min_r, max_lc
    ! This means the number of nodes in the list should be the total number of intersections + 2 ... the first and last nodes will have no cell reference. 
    ! Using just the closest will not work - especially for the big gaps ... need to go back to sorting by S value then check distances to see if there are any fix ups
    ! needed. 

    n_nodes = tot_n_intsects

    write(0,*) 'n_nodes:',n_nodes,tot_n_intsects

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_s)) deallocate(surf_s)
       allocate(surf_s(n_nodes),stat=ierr)
    else
       return
    endif

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_r)) deallocate(surf_r)
       allocate(surf_r(n_nodes),stat=ierr)
    else
       return
    endif

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_fl)) deallocate(surf_fl)
       allocate(surf_fl(n_nodes),stat=ierr)
    else
       return
    endif

    if (tot_n_intsects.gt.0) then 
       if (allocated(surf_int)) deallocate(surf_int)
       allocate(surf_int(n_nodes),stat=ierr)
    else
       return
    endif


    !surf_fl(1)=1
    !surf_int(1)=-1
    !surf_s(1)=min_lc
    !surf_r(1)=min_r

    !surf_fl(n_nodes)=1
    !surf_int(n_nodes)=-1
    !surf_s(n_nodes)=max_lc
    !surf_r(n_nodes)=min_r


    ! Populate the rest of the node list with intersection data 

    node_cnt = 0

    do if = 1,n_field_lines
       do in = 1,field_line(if)%int_tot
          node_cnt = node_cnt + 1
          surf_fl(node_cnt)=if
          surf_int(node_cnt)=in
          surf_s(node_cnt)=field_line(if)%int_data(in)%lc
          surf_r(node_cnt)= field_line(if)%xs
          !          node_used(node_cnt) =0
       end do
    end do


    write(0,*) 'node_cnt:',node_cnt

    call sort_arrays(0,node_cnt,surf_s,surf_r,surf_fl,surf_int)


    ! if block averaging is on
    ! calculate intersection data separations

    if (opt_block_av) then 

       tmp_node_cnt = 0

       if (tot_n_intsects.gt.0) then 
          if (allocated(surf_sep)) deallocate(surf_sep)
          allocate(surf_sep(n_nodes),stat=ierr)
       else
          return
       endif

       if (tot_n_intsects.gt.0) then 
          if (allocated(av_group)) deallocate(av_group)
          allocate(av_group(n_nodes),stat=ierr)
       else
          return
       endif

       av_group_cnt = 1

       !if (tot_n_intsects.gt.0) then 
       !   if (allocated(surf_ang)) deallocate(surf_ang)
       !   allocate(surf_ang(n_nodes),stat=ierr)
       !else
       !   return
       !endif


       do in = 1,n_nodes-1

          surf_sep(in) = sqrt((surf_s(in+1)-surf_s(in))**2 + (surf_r(in+1)-surf_r(in))**2)
          !surf_ang(in) = atan2c(surf_s(in+1)-surf_s(in),surf_r(in+1)-surf_r(in))

          !             surf_sep(in1,in2) = sqrt((surf_s(in1)-surf_s(in2))**2 + (surf_r(in1)-surf_r(in2))**2)

       end do

       surf_sep(n_nodes) = surf_sep(n_nodes-1)


       ! run through grouping and reordering

       min_sep = minval(surf_sep)

       write(0,*) 'Minimum intersection separation = ', min_sep

       do in = 1,n_nodes
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

          !          write (6,'(a,i9,2l8,10(1x,g18.6))') 'SEP:',in,surf_sep(in).gt.max_fact * last_sep,abs(surf_ang(in)-last_ang).gt.ang_limit,surf_sep(in),max_fact*last_sep,surf_s(in),surf_r(in),surf_ang(in)*raddeg,last_ang*raddeg,ang_limit*raddeg

          write (6,'(a,i9,l8,10(1x,g18.6))') 'SEP:',in,surf_sep(in).gt.max_fact * last_sep,surf_sep(in),max_fact*last_sep,surf_s(in),surf_r(in),av_group(in)

          ! Test to see if next point could be out of series
          !if (surf_sep(in).gt.max_fact * last_sep.or.abs(surf_ang(in)-last_ang).gt.ang_limit) then 
          if (surf_sep(in).gt.max_fact * last_sep) then 
             it = 1
             leave = .false.

             ! Look for additional points that are in series
             do while (.not.leave)
                it = it+1
                if (it.ge.max_block_size.or.(it+in+1).ge.n_nodes) then 
                   leave = .true.
                   av_group_cnt = av_group_cnt + 1
                else

                   test_sep = sqrt((surf_s(in+it)-surf_s(in))**2 + (surf_r(in+it)-surf_r(in))**2)
                   ! found a closer point - swap it with the next one up

                   !test_ang = atan2c(surf_s(in+it)-surf_s(in),surf_r(in+it)-surf_r(in))

                   !                   write(6,'(a,3i8,l8,10(1x,g18.6))') 'test:',in,it,in+it,test_sep.lt.max_fact*last_sep,test_sep,max_fact*last_sep,surf_s(in+it),surf_r(in+it),test_ang*raddeg,last_ang*raddeg, abs(test_ang-last_ang)

                   write(6,'(a,3i8,l8,10(1x,g18.6))') 'test:',in,it,in+it,test_sep.lt.max_fact*last_sep,test_sep,max_fact*last_sep,surf_s(in+it),surf_r(in+it)

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
                         write(6,'(a,2i8,5(1x,g18.6))') 'Moving:',in,ix,in+ix,surf_r(in+ix+1), surf_s(in+ix+1)
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


       do in = 2,n_nodes
          test_sep = sqrt((surf_s(in+1)-surf_s(in))**2 + (surf_r(in+1)-surf_r(in))**2)
          write(6,'(a,i8,l8,10(1x,g18.6))') 'Nodes:', in,surf_sep(in).gt.max_fact*surf_sep(in-1),surf_r(in),surf_s(in),surf_sep(in),test_sep,av_group(in)


       end do

       !
       ! Allocate temp arrays for averaging
       !

       write(0,*) 'Avgroup:',av_group_cnt, av_group(n_nodes)

       av_group_cnt = av_group(n_nodes) +2



       if (av_group_cnt.gt.0) then 
          if (allocated(av_s)) deallocate(av_s)
          allocate(av_s(av_group_cnt),stat=ierr)
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_r)) deallocate(av_r)
          allocate(av_r(av_group_cnt),stat=ierr)
       else
          return
       endif

       ! record max and min r values in each grouping
       if (av_group_cnt.gt.0) then 
          if (allocated(av_min_r)) deallocate(av_min_r)
          allocate(av_min_r(av_group_cnt),stat=ierr)
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_max_r)) deallocate(av_max_r)
          allocate(av_max_r(av_group_cnt),stat=ierr)
       else
          return
       endif

       av_min_r = max_r + 1.0
       av_max_r = min_r - 1.0


       group_id = 0

       do in = 1,n_nodes


          if (av_group(in).ne.group_id) then 
             ! finish up last group
             if (group_id.ne.0) then
                av_s(group_id+1) = av_s(group_id+1)/av_cnt
                av_r(group_id+1) = av_r(group_id+1)/av_cnt
             endif

             ! start next group
             group_id = av_group(in)
             av_cnt = 1
             av_s(group_id+1) = av_s(group_id+1) + surf_s(in)
             av_r(group_id+1) = av_r(group_id+1) + surf_r(in)

             av_min_r(group_id+1) = min(av_min_r(group_id+1),surf_r(in))
             av_max_r(group_id+1) = max(av_max_r(group_id+1),surf_r(in))

          else
             ! continue counting current group
             av_cnt = av_cnt + 1
             av_s(group_id+1) = av_s(group_id+1) + surf_s(in)
             av_r(group_id+1) = av_r(group_id+1) + surf_r(in)

             av_min_r(group_id+1) = min(av_min_r(group_id+1),surf_r(in))
             av_max_r(group_id+1) = max(av_max_r(group_id+1),surf_r(in))

          endif

          !          write(6,'(a,2i10,10(1x,g18.8))') 'Nodes2:',in,group_id,surf_s(in),surf_r(in),av_min_r(group_id+1),av_max_r(group_id+1)

       end do

       ! finish off last group

       av_s(group_id+1) = av_s(group_id+1)/av_cnt
       av_r(group_id+1) = av_r(group_id+1)/av_cnt


       ! add first and last points

       av_r(1) = min_r
       av_s(1) = min_lc
       av_min_r(1) = min_r
       av_max_r(1) = min_r

       av_r(av_group_cnt) = min_r
       av_s(av_group_cnt) = max_lc
       av_min_r(av_group_cnt) = min_r
       av_max_r(av_group_cnt) = min_r


    else

       ! copy over all nodes without block averaging - in case they fix the bugs in the data

       av_group_cnt = n_nodes +2

       if (av_group_cnt.gt.0) then 
          if (allocated(av_s)) deallocate(av_s)
          allocate(av_s(av_group_cnt),stat=ierr)
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_r)) deallocate(av_r)
          allocate(av_r(av_group_cnt),stat=ierr)
       else
          return
       endif
       ! record max and min r values in each grouping
       if (av_group_cnt.gt.0) then 
          if (allocated(av_min_r)) deallocate(av_min_r)
          allocate(av_min_r(av_group_cnt),stat=ierr)
       else
          return
       endif

       if (av_group_cnt.gt.0) then 
          if (allocated(av_max_r)) deallocate(av_max_r)
          allocate(av_max_r(av_group_cnt),stat=ierr)
       else
          return
       endif


       do in = 1,n_nodes
          av_s(in+1) = surf_s(in)
          av_r(in+1) = surf_r(in)
          av_min_r(in+1) = surf_r(in)
          av_max_r(in+1) = surf_r(in)
       end do

       ! add end points

       av_r(1) = min_r
       av_s(1) = min_lc
       av_min_r(1) = min_r
       av_max_r(1) = min_r

       av_r(av_group_cnt) = min_r
       av_s(av_group_cnt) = max_lc
       av_min_r(av_group_cnt) = min_r
       av_max_r(av_group_cnt) = min_r


    endif

    write(0,'(a,10(1x,g18.8))') 'Average:',min_r,max_r
    do in = 1,av_group_cnt
       write(6,'(a,i10,10(1x,g18.8))') 'AV:',in,av_s(in),av_r(in),av_max_r(in),av_min_r(in)
    end do



    ! go through the data and remove degeneracies in tangency and wall points

    dir = 1.0   
    last_face = 1

    do in = 1,av_group_cnt-1
       last_dir = dir
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


    end do


    ! Now run through and identify IN, OUT and tangency points - tangency points only occur at the top of a limiter ... not sure if IN/OUT are required


    ! allocate av_type
    if (av_group_cnt.gt.0) then 
       if (allocated(av_type)) deallocate(av_type)
       allocate(av_type(av_group_cnt),stat=ierr)
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

       write(0,'(a,i6,10(1x,g18.8))') 'DIR:',in,dir,last_dir,av_type(in),av_r(in),av_s(in)
       write(6,'(a,i6,10(1x,g18.8))') 'DIR:',in,dir,last_dir,av_type(in),av_r(in),av_s(in)

    end do


    av_type(av_group_cnt) = SURFACE_END

    !
    ! Tangency points
    !

    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_r)) deallocate(av_tan_r)
       allocate(av_tan_r(av_tan_cnt),stat=ierr)
       av_tan_r = 0.0
    else
       return
    endif



    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_s)) deallocate(av_tan_s)
       allocate(av_tan_s(av_tan_cnt),stat=ierr)
       av_tan_s = 0.0
    else
       return
    endif


    ! Need index to limiter surface vertex as well
    if (av_tan_cnt.gt.0) then 
       if (allocated(av_tan_ind)) deallocate(av_tan_ind)
       allocate(av_tan_ind(av_tan_cnt),stat=ierr)
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
       av_wall_s = 0.0
    else
       return
    endif

    if (av_wall_cnt.gt.0) then 
       if (allocated(av_wall_r)) deallocate(av_wall_r)
       allocate(av_wall_r(av_wall_cnt),stat=ierr)
       av_wall_r = 0.0
    else
       return
    endif


    tmp_tan_cnt = 0
    tmp_wall_cnt = 0

    ! Collect tangency and wall data - check to make sure there are no two consective points with the same R value - i.e. flat tangency point

    do in = 1,av_group_cnt

       if (av_type(in).eq.TANGENCY) then 
          tmp_tan_cnt = tmp_tan_cnt + 1
          ! record r,s and index into surface array for the tangency point
          av_tan_r(tmp_tan_cnt) = av_r(in)
          av_tan_s(tmp_tan_cnt) = av_s(in)
          av_tan_ind(tmp_tan_cnt) = in
          if (in.ne.av_group_cnt) then 
             if(av_type(in+1).eq.TANGENCY) then 
                double_tan = double_tan+1
             endif
          endif
       endif

       if (av_type(in).eq.WALL) then 
          tmp_wall_cnt = tmp_wall_cnt + 1
          av_wall_r(tmp_wall_cnt) = av_r(in)
          av_wall_s(tmp_wall_cnt) = av_s(in)
          if (in.ne.av_group_cnt) then 
             if(av_type(in+1).eq.WALL) then 
                double_wall = double_wall+1
             endif
          endif
       endif

    end do

    write(0,'(a,6i10)') 'Tan:',av_tan_cnt, av_wall_cnt,double_tan,double_wall
    do in = 1,max(av_tan_cnt,av_wall_cnt)

       if(in.le.av_wall_cnt) then 
          write(6,'(a,i6,l8,10(1x,g18.8))') 'Wall:',in,av_wall_s(in).lt.av_tan_s(in),av_wall_r(in),av_wall_s(in)
       endif

       if (in.le.av_tan_cnt) then 
          write(6,'(a,i6,l8,10(1x,g18.8))') 'Tan :',in,av_tan_s(in).lt.av_wall_s(in+1),av_tan_r(in),av_tan_s(in)
       endif
    end do

    ! Find max and min from the reduced intersection data ... loop through revised wall to get this data

    r_limiter_max = min_r -1.0
    r_limiter_min = max_r +1.0

    s_limiter_max = min_lc -1.0
    s_limiter_min = max_lc +1.0


    do in = 1,av_group_cnt

       r_limiter_max = max(r_limiter_max,dble(av_r(in)))
       r_limiter_min = min(r_limiter_min,dble(av_r(in)))

       s_limiter_max = max(s_limiter_max,dble(av_s(in)))
       s_limiter_min = min(s_limiter_min,dble(av_s(in)))

    end do



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

          write(0,'(a,2i10,l8,15(1x,g18.8))') 'Reduce tan:',in,in+1,(av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep,av_tan_r(in+1),av_tan_r(in),av_tan_r(in+1)-av_tan_r(in),min_tan_sep,av_tan_s(in),av_tan_s(in+1)
          write(6,'(a,2i10,l8,15(1x,g18.8))') 'Reduce tan:',in,in+1,(av_tan_r(in+1)-av_tan_r(in)) .lt. min_tan_sep,av_tan_r(in+1),av_tan_r(in),av_tan_r(in+1)-av_tan_r(in),min_tan_sep,av_tan_s(in),av_tan_s(in+1)
       endif
       
    end do


    ! put tangency points back into s order - this should be the same order as before the previous re-order - this is needed for generating cells while the r 
    ! r order is useful for generating rows
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,tan_ord_r,av_tan_ind)

    ! record the S order of the points 

    do in = 1,av_tan_cnt
       tan_ord_r(in) = in
       write(6,'(a,i8,10(1x,g18.8))') 'TAN_ORD_R:',in,av_tan_r(in),tan_ord_r(in)
    end do

    ! Sort the S indices by R to get tan_ord_r referencing the R order of the data while sorted by S
    call sort_arrays(0,av_tan_cnt,av_tan_r,av_tan_s,tan_ord_r,av_tan_ind)

    ! Re-sort to S - leaving tan_ord_r unchanged
    call sort_arrays(0,av_tan_cnt,av_tan_s,av_tan_r,av_tan_ind)

   
    do in = 1,av_tan_cnt
       write(6,'(a,i8,10(1x,g18.8))') 'TANS:',in,av_tan_r(in),av_tan_s(in),av_tan_ind(in)
    end do


  end subroutine calculate_limiter_surface


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
    else
       return
    endif

    if (tmp_cnt.gt.0) then 
       if (allocated(tmp_s)) deallocate(tmp_s)
       allocate(tmp_s(tmp_cnt),stat=ierr)
    else
       return
    endif

    if (tmp_cnt.gt.0) then 
       if (allocated(tmp_min_r)) deallocate(tmp_min_r)
       allocate(tmp_min_r(tmp_cnt),stat=ierr)
    else
       return
    endif

    if (tmp_cnt.gt.0) then 
       if (allocated(tmp_max_r)) deallocate(tmp_max_r)
       allocate(tmp_max_r(tmp_cnt),stat=ierr)
    else
       return
    endif



    !
    ! Copy data from av arrays 
    !


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
       tmp_r(num+1) = (tmp_max_r(num)+tmp_max_r(num+2))/2.0
    elseif (last_face.eq.2) then 
       ! last limiter facing was heading toward tip ... therefore move new R inward slightly
       if (tmp_r(num).ne.tmp_r(num+2)) then 
          write(0,'(a,i8,10(1x,g18.8))') 'ERROR: non-degenerate wall point being added:',num,tmp_r(num),tmp_r(num+2)
       endif
       tmp_r(num+1) = (tmp_min_r(num)+tmp_min_r(num+2))/2.0
    endif

    ! deallocate and reallocate the av arrays

    av_group_cnt = tmp_cnt

    if (av_group_cnt.gt.0) then 
       if (allocated(av_s)) deallocate(av_s)
       allocate(av_s(av_group_cnt),stat=ierr)
    else
       return
    endif

    if (av_group_cnt.gt.0) then 
       if (allocated(av_r)) deallocate(av_r)
       allocate(av_r(av_group_cnt),stat=ierr)
    else
       return
    endif

    if (av_group_cnt.gt.0) then 
       if (allocated(av_min_r)) deallocate(av_min_r)
       allocate(av_min_r(av_group_cnt),stat=ierr)
    else
       return
    endif

    if (av_group_cnt.gt.0) then 
       if (allocated(av_max_r)) deallocate(av_max_r)
       allocate(av_max_r(av_group_cnt),stat=ierr)
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


  end subroutine insert_limiter_point


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
    integer :: min_s_cnt, min_r_cnt

    real*8, allocatable :: ints(:,:),int_type(:,:),int_cnt(:)

    real*8, allocatable :: vert_rec1(:),vert_rec2(:),vert_type1(:),vert_type2(:)

    real*8 :: s_start,s_end
    real*8 :: slen, ssep1
    integer :: ncells

    integer :: rbnd_cnt,rbnd_add,new_bnds
    !integer :: sbnd_cnt,rbnd_cnt,rbnd_add,new_bnds

    integer :: npts1a, npts1b,npts2a,npts2b
    integer :: vert_cnt1,vert_cnt2,ntot1,ntot2
    real*8 :: r1,r2,r_sep_dist
    integer :: ierr,max_surf_ints
    integer :: min_cells

    max_s_sep = 0.5
    max_r_sep = 0.002

    min_cells = 5


    min_s_cnt = 5
    min_r_cnt = 1

    init_grid_size = max(3 * av_tan_cnt + 1, int((r_limiter_max-r_limiter_min)/max_r_sep))


    if (init_grid_size.gt.0) then 
       if (allocated(r_bnds)) deallocate(r_bnds)
       allocate(r_bnds(init_grid_size),stat=ierr)
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
             r_bnds(rbnd_cnt) = av_tan_r(in1)
             write(6,'(a,3i8,10(1x,g18.8))') 'Calc R_bnd:',in,in1,rbnd_cnt,r_bnds(rbnd_cnt),av_tan_r(in1),r_bnds(rbnd_cnt-1),av_tan_s(in1),av_tan_s(in)

          endif
       else
          rbnd_cnt = rbnd_cnt + 1
          r_bnds(rbnd_cnt) = av_tan_r(in1)
          write(6,'(a,3i8,10(1x,g18.8))') 'Calc R_bnd:',in,in1,rbnd_cnt,r_bnds(rbnd_cnt),av_tan_r(in1),av_tan_s(in1),av_tan_s(in)
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
    r_bnds(rbnd_cnt) = r_limiter_min

    ! put last r_bnd just inside the wall
    rbnd_cnt = rbnd_cnt + 1
    r_bnds(rbnd_cnt) = maxval(av_tan_r) + (r_limiter_max - maxval(av_tan_r)) * 0.95

    call sort_arrays(0,rbnd_cnt,r_bnds)
    !call sort_arrays(0,sbnd_cnt,s_bnds)
    !call sort_arrays(1,sbnd_cnt,s_bnds)

    ! run through R bounds and insert extras in situations where the separation is too large. 

    rbnd_add = 0

    do in = 1,rbnd_cnt-1
       new_bnds = int((r_bnds(in+1)-r_bnds(in))/max_r_sep)
       if (new_bnds.gt.1) then 
          r_sep_dist = (r_bnds(in+1)-r_bnds(in))/(new_bnds+1)
          do it = 1,new_bnds
             rbnd_add = rbnd_add + 1
             r_bnds(rbnd_cnt+rbnd_add) = r_bnds(in) + r_sep_dist * it
          end do

       endif

    end do

    rbnd_cnt = rbnd_cnt + rbnd_add

    ! re-sort boundaries

    call sort_arrays(0,rbnd_cnt,r_bnds)

    ! write out r_bnds

    do in = 1,rbnd_cnt
       write(6,'(a,i8,10(1x,g18.8))') 'R_BNDS:',in,r_bnds(in),r_bnds(in)-r_bnds(in-1)
    end do


    !
    ! Ok - now find intersection points with structure on this new grid for each R cell boundary
    !


    ! Allocate arrays to hold intersection data

    max_surf_ints = 2 * av_tan_cnt + 2

    if (rbnd_cnt.gt.0.and.av_tan_cnt.gt.0) then 
       if (allocated(int_cnt)) deallocate(int_cnt)
       allocate(int_cnt(rbnd_cnt),stat=ierr)
       int_cnt = 0
    else
       return
    endif

    if (rbnd_cnt.gt.0.and.av_tan_cnt.gt.0) then 
       if (allocated(ints)) deallocate(ints)
       allocate(ints(rbnd_cnt,max_surf_ints),stat=ierr)
       ints = 0
    else
       return
    endif

    if (rbnd_cnt.gt.0.and.av_tan_cnt.gt.0) then 
       if (allocated(int_type)) deallocate(int_type)
       allocate(int_type(rbnd_cnt,max_surf_ints),stat=ierr)
       int_type = 0
    else
       return
    endif

    int_cnt = 0.0
    int_type = 0.0
    ints = 0.0

    ! First r_bnd is assumed to only have 2 points - one at max and one at min
    int_cnt(1) = 2
    ints(1,1) = min_lc
    int_type(1,1) = SURFACE_START
    ints(1,2) = max_lc
    int_type(1,2) = SURFACE_END

    ! Find all the other intersections
    do in = 2,rbnd_cnt

       call find_field_line_intsects(max_surf_ints,r_bnds(in),int_cnt(in),int_type(in,:),ints(in,:))

    end do

    ! write out intersection points
    do in = 1,rbnd_cnt
       write(6,'(a,i8,10(1x,g18.8))')    'INTS  :',in,int_cnt(in)
       write(0,'(a,i8,10(1x,g18.8))')    'INTS  :',in,int_cnt(in)
       do it = 1,int_cnt(in)
          write(6,'(a,2i8,10(1x,g18.8))') '  DATA:', in,it,int_type(in,it),ints(in,it),r_bnds(in)
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

    max_nrings = av_tan_cnt * rbnd_cnt
    max_nknots = 200
    max_npoly = max_nrings * max_nknots

    ! Allocate polygon arrays
    ! rvp
    if (max_npoly.gt.0) then 
       if (allocated(rvp)) deallocate(rvp)
       allocate(rvp(4,max_npoly),stat=ierr)
    else
       return
    endif

    ! zvp
    if (max_npoly.gt.0) then 
       if (allocated(zvp)) deallocate(zvp)
       allocate(zvp(4,max_npoly),stat=ierr)
    else
       return
    endif

    !
    ! nvp
    if (max_npoly.gt.0) then 
       if (allocated(nvp)) deallocate(nvp)
       allocate(nvp(max_npoly),stat=ierr)
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
    else
       return
    endif

    ! zcen
    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(zcen)) deallocate(zcen)
       allocate(zcen(max_nknots,max_nrings),stat=ierr)
    else
       return
    endif

    ! poly_ref
    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(poly_ref)) deallocate(poly_ref)
       allocate(poly_ref(max_nknots,max_nrings),stat=ierr)
    else
       return
    endif

    ! nknots
    if (max_nrings.gt.0.and.max_nknots.gt.0) then 
       if (allocated(nknots)) deallocate(nknots)
       allocate(nknots(max_nrings),stat=ierr)
    else
       return
    endif


    ! vertex data for calling gen_ring

    if (max_nknots.gt.0) then 
       if (allocated(vert_rec1)) deallocate(vert_rec1)
       allocate(vert_rec1(max_nknots),stat=ierr)
    else
       return
    endif

    if (max_nknots.gt.0) then 
       if (allocated(vert_rec2)) deallocate(vert_rec2)
       allocate(vert_rec2(max_nknots),stat=ierr)
    else
       return
    endif


    if (max_nknots.gt.0) then 
       if (allocated(vert_type1)) deallocate(vert_type1)
       allocate(vert_type1(max_nknots),stat=ierr)
    else
       return
    endif

    if (max_nknots.gt.0) then 
       if (allocated(vert_type2)) deallocate(vert_type2)
       allocate(vert_type2(max_nknots),stat=ierr)
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



    do in = 1,rbnd_cnt -1 

       write(0,*) '5:',in


       r1 = r_bnds(in)
       r2 = r_bnds(in+1)

       if (grid_option.eq.1.or.in.eq.1) then 
          do it = 1,int_cnt(in)
             vert_rec1(it) = ints(in,it)
             vert_type1(it) = int_type(in,it)
             write(6,'(a,i8,10(1x,g18.8))') 'VERT1A:',it,vert_rec1(it),vert_type1(it)
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
       ! For grid option 2 - add all tangency data between its end points to each ring as vertices for r values less than the tangency point location
       ! Also for the first row we need to define all vertices subject to the constraints of minimum number of cells between tangency positions and maximum allowed S separation
       !

       if (grid_option.eq.2) then

          if (in.eq.1) then 
             do it = 1,av_tan_cnt
                if (av_tan_r(it).gt.r1) then
                   ntot1 = ntot1+1
                   vert_rec1(ntot1) = av_tan_s(it)
                   vert_type1(ntot1) = FIXED_VERTEX
                endif
             end do

             vert_cnt1 = ntot1
             call sort_arrays(1,vert_cnt1,vert_rec1,vert_type1)

             do is = 1,vert_cnt1-1
                ! place additional vertices between each tangency/surface using the grid generation contraints
                slen = vert_rec1(is+1) - vert_rec1(is)
                ncells = int(slen/max_s_sep) + 1

                if (ncells.lt.min_cells) ncells = min_cells

                ssep1 = abs(slen)/ncells 

                do it = 1,ncells-1
                   vert_rec1(ntot1+it) = vert_rec1(is)+ ssep1 * it
                   vert_type1(ntot1+it) = FIXED_VERTEX
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
          call sort_arrays(1,vert_cnt1,vert_rec1,vert_type1)
          call sort_arrays(1,vert_cnt2,vert_rec2,vert_type2)

       end if


       write(6,'(a,2i8,10(1x,g18.8))') 'VERT1:',vert_cnt1,ntot1
       do it = 1,ntot1
          write(6,'(a,i8,10(1x,g18.8))') 'VERT1D:',in,r1,vert_rec1(it),vert_type1(it)
       end do
       write(6,'(a,2i8,10(1x,g18.8))') 'VERT2:',vert_cnt2,ntot2
       do it = 1,ntot2
          write(6,'(a,i8,10(1x,g18.8))') 'VERT2D:',in,r2,vert_rec2(it),vert_type2(it)
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

       ! assign start/end indices on each row
       do it = 1,vert_cnt1
          write(0,'(a,10i8)') 'VERT:',in,it,vert_cnt1,npts1a,npts1b,npts2a,npts2b

          if (vert_type1(it).eq.SURFACE_START) then 
             npts1a = it
          elseif (vert_type1(it).eq.TANGENCY.or.vert_type1(it).eq.SURFACE_END) then 
             npts1b = it

             ! Now find the matching indices in vert_rec2
             s_start = vert_rec1(npts1a)
             s_end = vert_rec1(npts1b)


             write(0,'(a,6i8,10(1x,g18.8))') 'VERTB:',in,it,vert_cnt1,vert_cnt2,npts1a,npts1b,s_start,s_end

             ! loop through vert_rec2

             npts2a = 0
             npts2b = 0

             do is = 1,vert_cnt2
                !write(0,*) 'NPTS2A:',is,vert_cnt2,vert_rec2(is)
                if (vert_rec2(is).ge.s_start.and.vert_rec2(is).le.s_end) then 
                   ! this should assign npts2a and npts2b to appropriate values
                   if (npts2a.eq.0) then 
                      npts2a = is
                   endif

                   npts2b = is

                elseif (vert_rec2(is).gt.s_end) then 
                   exit
                endif
                !write(0,*) 'NPTS2B:',is,vert_cnt2,vert_rec2(is)
             end do


             if (npts2a.ne.0.and.npts2b.ne.0) then 

                write(0,'(a,6i8,10(1x,g18.8))') 'GEN_RING: CALL: ', in,it,npts1a,npts1b,npts2a,npts2b,r1,r2,s_start,s_end
                call gen_ring(npts1a,npts1b,ntot1,r1,vert_rec1,vert_type1,npts2a,npts2b,ntot2,r2,vert_rec2,vert_type2,max_nknots,max_s_sep,min_cells)

             else

                write(0,'(a,6i8,10(1x,g18.8))') 'GEN_RING: NO CALL: ', in,it,npts1a,npts1b,npts2a,npts2b,r1,r2,s_start,s_end

             endif

             write (0,*) 'After gen_ring'

             ! If the point is a tangency then it is also the next surface_start point
             if (vert_type1(it).eq.TANGENCY) then 
                npts1a = it
             endif

          endif

          write(0,*) '1:'

       end do



       ! if new vertices need to be retained then sort the vert_rec2 array and assign it to vert_rec1 
       ! for grid option 2 - retain vertices as fixed for next row

       write(0,*) '2:'

       if (grid_option.eq.2) then 
          vert_cnt2 = ntot2
          call sort_arrays(1,vert_cnt2,vert_rec2,vert_type2)
          ! assign "2" arrays to "1" arrays for the next row in grid

       endif

       write(0,*) '3:'

       ! move over vertices by one row ... either keeping ones added or not
       vert_cnt1 = vert_cnt2
       vert_rec1 = vert_rec2
       vert_type1 = vert_type2

       write(0,*) '4:'


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
       write(0,'(a,2i8,10(1x,g18.8))') 'S1:',in,in-npts1a+1,s1(in-npts1a+1),s1_type(in-npts1a+1)
       write(6,'(a,2i8,10(1x,g18.8))') 'S1:',in,in-npts1a+1,s1(in-npts1a+1),s1_type(in-npts1a+1)
    end do

    do in = npts2a,npts2b
       s2(in-npts2a+1) = vert_rec2(in)
       s2_type(in-npts2a+1) = vert_type2(in)
       write(0,'(a,2i8,10(1x,g18.8))') 'S2:',in,in-npts2a+1,s2(in-npts2a+1),s2_type(in-npts2a+1)       
       write(6,'(a,2i8,10(1x,g18.8))') 'S2:',in,in-npts2a+1,s2(in-npts2a+1),s2_type(in-npts2a+1)       
    end do


    !
    ! Generate polygons and additional vertices as required
    !

    if (npts1.eq.2.and.npts2.eq.2) then 
       ! determine number of vertices req'd

       ! Average slen
       !slen = max(abs(s1(1)-s1(2)),abs(s2(1)-s2(2)))

       slen = abs(s1(1)-s1(2))

       ncells = int(slen/max_slen) + 1

       if (ncells.lt.min_cells) ncells = min_cells

       ssep1 = abs(s1(1)-s1(2))/ncells 
       ssep2 = abs(s2(1)-s2(2))/ncells       

       do in = 1,ncells-1
          s1(npts1+in) = s1(1) + ssep1 * in
          s1_type(npts1+in) = NEW_VERTEX
          s2(npts1+in) = s2(1) + ssep2 * in
          s2_type(npts1+in) = NEW_VERTEX
       end do

       npts1 = npts1 + ncells -1
       npts2 = npts2 + ncells -1

       call sort_arrays(1,npts1,s1,s1_type)
       call sort_arrays(1,npts2,s2,s2_type)

    elseif (npts1.eq.2.and.npts2.gt.2) then 
       ! Case with one or more tangency points on second boundary

       slen = abs(s1(1)-s1(2))

       ncells = int(slen/max_slen) + 1

       if (ncells.lt.min_cells) ncells = min_cells

       ssep1 = abs(s1(1)-s1(2))/ncells 

       write(0,*) 'GR:NPTS:',npts1,npts2,ncells,slen

       do in = 1,ncells-1
          s1(npts1+in) = s1(1) + ssep1 * in
          s1_type(npts1+in) = NEW_VERTEX
       end do

       ncells2 = ncells - npts2 +2

       slen = abs(s2(1)-s2(npts2))

       write(0,'(a,3i8,10(1x,g18.8))') 'GR:NPTS2:',npts2,ncells,ncells2,slen

       dist = 0.0
       do in = 1,npts2-1
          dist(in) = abs(s2(in+1)-s2(in))/slen * ncells2
          write(0,'(a,i8,10(1x,g18.8))') 'GR:DIST:',in,dist(in),s2(in+1),s2(in)
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
          WRITE(0,'(a,2i8,10(1x,g18.8))') 'dist2:',in,it,cell_dist(it),dist(it)
       end do

       !write(0,*) 'DIST sum:',sum(cell_dist),ncells2-1

       cells_added = 0
       do in = 1,npts2-1
          ssep2 = abs(s2(in+1)-s2(in))/(cell_dist(in)+1.0)
          do it = 1,cell_dist(in)-1
             cells_added = cells_added + 1
             s2(npts2+cells_added) = s2(in) + ssep2 * it
             s2_type(npts2+cells_added) = NEW_VERTEX
             write(0,'(a,3i8,10(1x,g18.8))') 'Adding:',in,it,cells_added,ssep2,s2(npts2+cells_added),s2_type(npts2+cells_added)
          end do
       end do

       npts1 = npts1 + ncells -1
       npts2 = npts2 + cells_added

       call sort_arrays(1,npts1,s1,s1_type)
       call sort_arrays(1,npts2,s2,s2_type)

    elseif (npts1.gt.2.and.npts2.eq.2) then 
       ! points pre-specified along r1 ... none defined on r2 ... do proportional spacing to match

       slen = abs(s1(npts1)-s1(1))

       do in = 1,npts1
          dist(in) =   abs((s1(in)-s1(1)))/slen
       end do

       ! add points with same spacing along r2


       slen = abs(s2(npts2)-s2(1))

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
                   s2(npts2+cells_added) = abs(s1(is)-s1(last_tan))/abs(s1(in)-s1(last_tan)) * abs(s2(it)-s2(it-1)) + s2(it-1)
                   s2_type(npts2+cells_added) = NEW_VERTEX
                   write(6,'(a,4i8,10(1x,g18.8))') 'NEW S2:',in,it,is,last_tan,cells_added,s1(is),s1(in),s1(last_tan),s2(it),s2(it-1),s2(npts2+cells_added)
                end do
                last_tan= in
                ! if we have just finished the match for the last tangency point then process the end of the row
                if (it.eq.npts2-1) then 

                   do is = in+1,npts1-1
                      cells_added = cells_added + 1
                      s2(npts2+cells_added) = abs(s1(is)-s1(last_tan))/abs(s1(npts1)-s1(last_tan)) * abs(s2(npts2)-s2(it)) + s2(it)
                      s2_type(npts2+cells_added) = NEW_VERTEX
                      write(6,'(a,4i8,10(1x,g18.8))') 'NEW S2:',in,it,is,last_tan,cells_added,s1(is),s1(in),s1(last_tan),s2(it),s2(it-1),s2(npts2+cells_added)
                   end do

                endif
                exit
             end if

          end do

       end do

       npts2 = npts2 + cells_added 
       call sort_arrays(1,npts2,s2,s2_type)

       write(0,'(a,i8,10(1x,g18.8))') 'NEW s2:',npts2,cells_added

    endif


    ! for grid generation number of points on each side should be the same - matching vertices ... if not an error has occurred
    if (npts1.eq.npts2) then 

       ! order of listing vertices is important as is ordering cells
       ! the base coordinate system for the ITER data is from smaell R to larger R >0 
       ! and from S < 0 to S > 0 ... the cells are organized wiith DOWN at -S between vertices 1,2
       ! ... so a cell is   1 = s1(in)  2 = s2(in)   3 = s2(in+1)  4 = s1(in+1)

       ! Lets create the ring and cells

       nrings = nrings + 1
       nknots(nrings) = npts1 -1 

       write(6,'(a,4i8,10(1x,g18.8))') 'GEN:',npts1,npts2,nrings,nknots(nrings)

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
          write(6,'(a,3i8,10(1x,g18.8))') 'POLY:',in,npoly,poly_ref(in,nrings),(rvp(it,npoly),zvp(it,npoly),it=1,nvp(npoly)),rcen(in,nrings),zcen(in,nrings)

       end do

    else

       write(0,'(a,2i10,10(1x,g18.8))') 'ERROR in vertex generation',npts1,npts2
       stop

    endif



    ! copy new vertices and types into vert_rec2 - vert_rec1 was completed on previous iterations

    write(0,'(a,3i8,10(1x,g18.8))') 'Assign vert_rec2:', maxpts,npts2,ntot2

    do in = 1,npts2
       if (s2_type(in).eq.NEW_VERTEX) then 
          ntot2 = ntot2 + 1
          write(0,'(a,3i8,10(1x,g18.8))') 'REC VERT:', in,ntot2,maxpts,s2(in),s2_type(in)

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
    call find_free_unit_number(iunit)

    filename = 'grid.out'
    
    open(iunit,file=filename,form='formatted')

    write(0,'(a,3i8,10(1x,g18.8))') 'Grid:', npoly,nrings
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


    filename = 'poly.out'
    
    open(iunit,file=filename,form='formatted')

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
             if (int_type_base.gt.0.0) then 
                int_type(int(int_cnt)) = SURFACE_START
             else
                int_type(int(int_cnt)) = SURFACE_END
             endif

             int_type(int(int_cnt)) = int_type_base
             int_type_base = int_type_base * -1
             intersects(int(int_cnt)) = s_int

             write(6,'(a,2i8,10(1x,g18.8))') 'Int    :',in,sect_type,int_cnt,r_int,s_int,int_type_base,r_start,s_start,r_end,s_end,r_int,s_int

          endif

       endif

    end do

    ! Add tangency points if any


    do in = 1,av_tan_cnt
       if (av_tan_r(in).eq.r_line) then 
          int_cnt = int_cnt + 1.0
          int_type(int(int_cnt)) = TANGENCY
          intersects(int(int_cnt)) = av_tan_s(in)
          write(6,'(a,2i8,10(1x,g18.8))') 'Int Tan:',in,sect_type,int_cnt,av_tan_r(in),av_tan_s(in)
       end if
    end do


    ! sort 

    call sort_arrays(0,int(int_cnt),intersects,int_type)

    ! filter out any duplicates (points with S values closer than the minimum allowed - these should be intersections at tangency points ... if given a choice keep
    ! the one identified as a tangency point.


    deleted = 0

    do in = 1,int_cnt-1

       if (int_type(in).ne.DELETE_POINT) then
          do it = in+1,int_cnt
             if (int_type(it).ne.DELETE_POINT) then 

                if (abs(intersects(it)-intersects(in)).lt.eps) then

                   ! consecutive intersections too close together
                   if (int_type(it).eq.TANGENCY.and.int_type(in).ne.TANGENCY) then 
                      int_type(in) = DELETE_POINT

                   elseif (int_type(in).eq.TANGENCY.and.int_type(it).ne.TANGENCY) then 
                      int_type(it) = DELETE_POINT
                   else
                      int_type(in) = DELETE_POINT
                   endif

                   deleted = deleted + 1.0
                   write(6,'(a,2i8,l4,10(1x,g18.8))') 'Delete:',in,it, (abs(intersects(it)-intersects(in)).lt.eps),intersects(it),intersects(in),eps,int_type(in),int_type(it)

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
       write(6,'(a,10(1x,g18.8))') 'surface intersections error:', int_cnt, r_line
       do in = 1,int_cnt
          write(6,'(a,i8,10(1x,g18.8))') ' INT sect:', in,int_type(in),intersects(in)
       end do
       stop
    else   

       write(6,'(a,10(1x,g18.8))') 'surface intersections output:', int_cnt, r_line
       do in = 1,int_cnt
          write(6,'(a,i8,10(1x,g18.8))') ' INT sect:', in,int_type(in),intersects(in)
       end do

    endif



  end subroutine find_field_line_intsects

  subroutine assign_grid_to_divimp(maxnrs,maxnks,mves,nrs,nks,&
           nves,rves,zves,&
           npolyp,korpg,&
           nvertp,rvertp,zvertp,&
           rs,zs)
    implicit none
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
       call errmsg('Too many polygons in ribbon grid:',npoly)
    endif



    if (nrs.le.maxnrs) then 
       nrs = nrings
     
       do ir = 1,nrings
          if (nknots(ir).le.maxnks) then 
             nks(ir) = nknots(ir)
             do ik = 1,nknots(ir)
                rs(ik,ir) = rcen(ik,ir)
                zs(ik,ir) = zcen(ik,ir)
                korpg(ik,ir) = poly_ref(ik,ir)
                write(6,'(a,3i8,10(1x,g18.8))') 'POLY COPY:',ik,ir,korpg(ik,ir),rs(ik,ir),zs(ik,ir)
             end do 
          else
             call errmsg('Too many knots on ribbon grid ring#:',ir)
          endif
       end do
    else
       call errmsg('Too many rings in ribbon grid:',nrings)
    endif




    ! calculate rves,zves by running along the cell ends and the limiter surface line. 
    ! one source of these numbers could be the intersects array from field_line_intersections since these
    ! points are used to define the polygon surfaces at the ends of the rings. 

    ! Do I need a connection map to build this or is there an easier way?

    ! for now just assign the averaged surface since all intersection points lie on this surface
    ! ... may need to rebuild later from connection map data
    
    ! End points should already be included 
    if (nves.le.mves) then 
       nves = av_group_cnt
       do in = av_group_cnt,1,-1
          rves(in) = av_r(in)
          zves(in) = av_s(in)
          write(6,'(a,i8,5(1x,g18.8))') 'Vessel:',in,av_r(in),av_s(in)
       end do
       
       ! Add 

    else
       call errmsg('Too many elements in ribbon grid wall specification:',av_group_cnt)
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





end module castem_field_line_data

