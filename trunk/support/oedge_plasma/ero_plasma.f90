module ero_plasma
  use error_handling
  use allocate_arrays
  use utilities
  use oedge_plasma_interface

  integer :: specunit

  integer, parameter :: line_length = 512
  character*6, parameter :: line_form='(a512)'
  character*(line_length) :: buffer

  !
  ! Default plasma and geometry file names
  !
  character*50 ::  erodiv_plasma_fn ='erodiv_plasma.m'
  character*50 ::  erodiv_geometry_fn ='erodiv_geometry.m'
  character*50 ::  out_format = 'ascii_table'

  ! options
  integer :: shift_opt=0,override_offsets=0


  !
  ! Base specifications of divimp plasma and number of ERO blocks in file
  !
  character*(line_length) :: div_plasma_file, div_geo_file
  integer :: n_ero_blocks


  !
  real*8 :: xhat(2),yhat(2) ! vector coordinate transformation used both in the plasma and geometry code - calculated in plasma code

  !
  ! ERO block data
  !
  integer,parameter :: nvert=4
  real*8 :: vert(nvert,2)
  real*8 :: roffset,zoffset
  real*8 :: tor_extent
  integer :: interp_opt,extrap_opt,eroshape_opt
  integer :: remap_bfield

  integer :: divpart_opt   ! 0 off 1 on
  integer :: erodiv_tor_part_mult

  character*(line_length) :: erobg_out_name,eroshape_out_name,divpart_in_name,eropart_out_name

  real*8 :: xmin,xmax,ymin,ymax,zmin,zmax

  !
  ! Grid resolution can be specified by either nx,ny or delx,dely - number of points is rounded off when delx.dely is specified. 
  !
  integer :: nx,ny
  real*8 :: delx,dely

  logical :: divimpbg_loaded = .false.

  ! ero background plasma storage

  real*8,allocatable :: ero_plasma_out(:,:,:)

  ! ERO surface calculation variables

  integer,parameter :: tor_points= 25
  integer :: n_surf,n_pol,n_tor
  real*8 :: x0,y0,z0, alp0, rtor, rz0
  real*8,allocatable :: tpx(:,:),tpy(:,:),tpz(:,:)

  real*8,parameter :: dist_conversion = 1000.0  ! convert m to mm for scale in geometry

  public read_ero_plasma_headers, read_ero_plasma_block, calc_ero_plasma,output_ero_plasma,output_ero_plasma_m,calc_ero_surface,output_ero_surface_m,output_ero_surface,init_ero_options

contains

  subroutine init_ero_options(shiftoption,offset_override)
    implicit none
    integer :: shiftoption,offset_override

    ! option to move the corner points out a bit in the ERO simulation so that particle losses at the edge of hte simulation volume are less of an issue
    shift_opt = shiftoption
    override_offsets = offset_override

  end subroutine init_ero_options



  subroutine read_ero_plasma_headers(erospec_name,nblocks)
    implicit none
    character*(*) :: erospec_name
    integer :: nblocks


    integer :: ios
    logical :: done,read_plasma_file,read_geo_file

    call find_free_unit_number(specunit)

    open(specunit,file=trim(erospec_name),status='old',form='formatted',iostat=ios)

    if (ios.ne.0) then 
       call errmsg('ERROR Opening Specification file: Case = ',trim(erospec_name))
       stop 'Spec file not found'
    endif

    ! set the default output format for the ero plasma file to be ascii table 
    ! the alternative is matlab format. 
    out_format = 'ascii_table'

    ! Loop through specification file looking for tagged data lines

    done = .false.
    read_plasma_file=.false.
    read_geo_file=.false.

    do while (.not.done)  

       read(specunit,line_form,iostat=ios) buffer

       if (ios.ne.0) then 
          call errmsg('ERROR read_ero_spec_headers: Err code = ',ios)
          stop 'Unexpected error reading '//trim(erospec_name)
       endif

       !write(0,'(a,a,a,a)') 'buffer:',trim(buffer),':'

       if (buffer(1:15).eq.'#EROBLOCK') then
          done=.true. 
       elseif (buffer(1:10).eq.'#NUMBLOCKS') then 
          read(buffer(11:),*) n_ero_blocks

          !write(0,'(a,i10)') 'nblocks:',n_ero_blocks
       elseif (buffer(1:14).eq.'#DIVPLASMAFILE') then 
          !read(specunit,line_form,iostat=ios) buffer
          
          read(buffer(15:),*) div_plasma_file
          read_plasma_file = .true.

          !write(0,'(a,a,a,a)') 'bgp:',trim(div_plasma_file),':'
          
          !if (ios.eq.0) then 
          !   div_plasma_file = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_spec_headers: read plasma file name: code = ',ios)
          !   stop 'Unexpected error reading '//trim(erospec_name)
          !endif

       elseif (buffer(1:11).eq.'#DIVGEOFILE') then 
          read(buffer(12:),*) div_geo_file
          read_geo_file = .true.

          !write(0,'(a,a,a,a)') 'grd:',trim(div_geo_file),':'

          !read(specunit,line_form,iostat=ios) buffer

          !if (ios.eq.0) then 
          !   div_geo_file = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_spec_headers: read geo file name: code = ',ios)
          !   stop 'Unexpected error reading '//trim(erospec_name)
          !endif

       elseif (buffer(1:10).eq.'#OUTFORMAT') then 
          read(buffer(11:),*) out_format
       endif

    end do

    !write(0,*) read_geo_file, read_plasma_file

    if ((.not.read_geo_file)) then
       call errmsg('ERROR read_ero_spec_headers','geometry file name not loaded')
       stop 'No geometry file name loaded'
    endif

    if  ((.not.read_plasma_file)) then 
       call errmsg('ERROR read_ero_spec_headers','plasma file name not loaded')
       stop 'No plasma file name loaded'
    endif


    nblocks = n_ero_blocks

    return
  end subroutine read_ero_plasma_headers





  subroutine read_ero_plasma_block
    implicit none

    logical :: done
    integer :: i,ios

    logical :: rv_read,zv_read,erobg_name_read


    ! Initialization and set defaults
    done=.false. 
    rv_read = .false.
    zv_read = .false.
    roffset = 0.0
    zoffset = 0.0
    interp_opt = 0
    extrap_opt = 0
    delx = 0.0
    dely = 0.0
    nx = 100
    ny = 100
    remap_bfield = 1
    eroshape_opt = -1  ! set to invalid option initially as a flag for whether the data was read and mapping requested

    do while (.not.done)  

       read(specunit,line_form,iostat=ios) buffer

       !write(0,'(a,a,a,a)') 'buffer:eb:',trim(buffer),':'

       if (ios.ne.0) then 
          ! need to exit loop if eof reached - not an error since could be last block
          done = .true. 
          cycle
       endif

       if (buffer(1:9).eq.'#EROBLOCK') then
          ! start of next block
          done=.true. 
       elseif (buffer(1:8).eq.'#OFFSETS') then 
          ! read R,Z coordinate offset 
          read(buffer(9:),*) roffset,zoffset
          !write(0,'(a,10(1x,g12.5))') 'offset:',roffset,zoffset

       elseif (buffer(1:10).eq.'#RVERTICES') then 
          ! read R vertices
          read(buffer(11:),*) (vert(i,1),i=1,4)
          rv_read = .true.
          !write(0,'(a,10(1x,g12.5))') 'rvert:',(vert(i,1),i=1,4)

       elseif (buffer(1:10).eq.'#ZVERTICES') then 
          ! read Z vertices
          read(buffer(11:),*) (vert(i,2),i=1,4)
          zv_read = .true.
          !write(0,'(a,10(1x,g12.5))') 'zvert:',(vert(i,2),i=1,4)

       elseif (buffer(1:11).eq.'#TOR_EXTENT') then 
          ! read toroidal extent
          read(buffer(12:),*) tor_extent
          !write(0,'(a,10(1x,g12.5))') 'tor_extent:',tor_extent

       elseif (buffer(1:5).eq.'#NXNY') then 
          ! read nx ny
          read(buffer(6:),*) nx,ny

          !write(0,'(a,2i10,10(1x,g12.5))') 'nxny:',nx,ny

       elseif (buffer(1:8).eq.'#DELTAXY') then 
          ! read delx dely
          read(buffer(9:),*) delx,dely
          !write(0,'(a,10(1x,g12.5))') 'dxdy:',delx,dely

       elseif (buffer(1:12).eq.'#INTERPOLATE') then 
          ! read interpolation option
          read(buffer(13:),*) interp_opt
          !write(0,'(a,2i10,10(1x,g12.5))') 'interp:',interp_opt

       elseif (buffer(1:12).eq.'#EXTRAPOLATE') then 
          ! read interpolation option
          read(buffer(13:),*) extrap_opt
          !write(0,'(a,2i10,10(1x,g12.5))') 'extrap:',extrap_opt


       elseif (buffer(1:13).eq.'#REMAP_BFIELD') then 
          ! read interpolation option
          read(buffer(14:),*) remap_bfield

          !write(0,'(a,2i10,10(1x,g12.5))') 'remap:',remap_bfield

       elseif (buffer(1:10).eq.'#EROBGFILE') then 
          ! read ERO background plasma output file name
          read(buffer(11:),*) erobg_out_name
          !read(specunit,line_form,iostat=ios) buffer
          erobg_name_read = .true.

          !write(0,'(a,a,a,a)') 'erobg:',trim(erobg_out_name),':'

          !if (ios.eq.0) then 
          !   erobg_out_name = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_plasma_block: read erobg output file name: code = ',ios)
          !   stop 'Unexpected error reading erobg_out_name'
          !endif

       elseif (buffer(1:13).eq.'#EROSHAPEFILE') then 
          ! read ERO background plasma output file name
          read(buffer(14:),*) eroshape_out_name
          !read(specunit,line_form,iostat=ios) buffer
          !eroshape_name_read = .true.

          !write(0,'(a,a,a,a)') 'eroshape_fn:',trim(eroshape_out_name),':'

          !if (ios.eq.0) then 
          !   erobg_out_name = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_plasma_block: read erobg output file name: code = ',ios)
          !   stop 'Unexpected error reading erobg_out_name'
          !endif
       elseif (buffer(1:12).eq.'#EROSHAPEOPT') then 
          ! read interpolation option
          read(buffer(13:),*) eroshape_opt
          !write(0,'(a,2i10,10(1x,g12.5))') 'eroshape_opt:',eroshape_opt

       elseif (buffer(1:11).eq.'#DIVPARTOPT') then 
          ! read interpolation option
          read(buffer(12:),*) divpart_opt
          write(0,'(a,2i10,10(1x,g12.5))') 'divpart_opt:',divpart_opt

       elseif (buffer(1:12).eq.'#DIVPARTFILE') then 
          ! read DIVIMP particle output file name
          read(buffer(13:),*) divpart_in_name
          !read(specunit,line_form,iostat=ios) buffer
          !eroshape_name_read = .true.

          !write(0,'(a,a,a,a)') 'divpart_fn:',trim(divpart_in_name),':'

       elseif (buffer(1:12).eq.'#EROPARTFILE') then 
          ! read name for ERO particle input file name
          read(buffer(13:),*) eropart_out_name
          !read(specunit,line_form,iostat=ios) buffer
          !eroshape_name_read = .true.

          !write(0,'(a,a,a,a)') 'eropart_fn:',trim(eropart_out_name),':'

       elseif (buffer(1:15).eq.'#ERODIVPARTMULT') then 
          ! read interpolation option
          read(buffer(16:),*) erodiv_tor_part_mult
          !write(0,'(a,2i10,10(1x,g12.5))') 'erodiv_tor_part_mult:',erodiv_tor_part_mult

       endif

    end do

    if ((.not.rv_read)) then
       call errmsg('ERROR read_ero_plasma_block','R vertices not loaded')
       stop 'R vertices not loaded'
    endif

    if ((.not.zv_read)) then
       call errmsg('ERROR read_ero_plasma_block','Z vertices not loaded')
       stop 'Z vertices not loaded'
    endif

    if ((.not.erobg_name_read)) then
       call errmsg('ERROR read_ero_plasma_block','ERO output plasma file name not read')
       stop 'Output plasma file name not read'
    endif



    return


  end subroutine read_ero_plasma_block




  subroutine calc_ero_plasma(errmsg_unit)
    implicit none

    integer :: ierr
    integer :: errmsg_unit 

    ! algorithm can easily be extended to 3D if required
    real*8 :: xvec(2),yvec(2),bfield(2),xrange,yrange
    real*8 :: xshift(2),yshift(2)
    real*8 :: rt,zt,xt,yt
    real*8 :: rxstep,zxstep,rystep,zystep

    real*8 :: tmp_offset(2),xyoffset(2)

    real*8 :: ne,te,ti,vb,ef,btot,br,bz,bt,psin
    real*8 :: ngrad,tegrad,tigrad
    real*8 :: bx,by
    
    integer :: ix,iy

    !
    ! Check to see if DIVIMP plasma data has been loaded (if not then load it)
    !

    if (.not.divimpbg_loaded) then 
       call load_oedge_data(div_geo_file,div_plasma_file,ierr)

       if (ierr.ne.0) then 
          call errmsg('Calc_ero_plasma: Error loading divimp files:',ierr)
          stop 'DIVIMP files could not be loaded'
       endif

       divimpbg_loaded = .true.

    endif

    ! initialize the oedge plasma loading code


    ! Handle offsets in the local code - not when getting the plasma data from divimp
    call set_oedge_plasma_opts(0.0d0,0.0d0,interp_opt,extrap_opt,errmsg_unit)

    !call set_oedge_plasma_opts(roffset,zoffset,interp_opt,extrap_opt,errmsg_unit)

    !
    !  Figure out how large the ero plasma array needs to be
    !

    !
    ! If delx and dely are both greater than zero then these are used for calculating grid resolution
    !   - nx and ny are over-written
    ! If either delx or dely are less than or equal to zero then nx and ny are used to determine grid resolution
    ! 
    !

    ! Calculate the xvector, yvector, and the xrange and yrange

    ! The vertices are specified as 4 points
    ! - X vector is from point 1 to 2 ... i.e. 2-1
    ! - Y vector is from point 1 to 4 ... i.e. 4-1    
    !
    ! Vertices should be listed as:
    !  1234      
    !       
    ! Spatially these can have any orientation but must be a rectangle with vertices listed either
    ! clockwise or counter-clockwise. 
    !  

    xvec = vert(2,:)-vert(1,:)
    yvec = vert(4,:)-vert(1,:)

    xrange =  sqrt(xvec(1)**2 + xvec(2)**2)
    yrange =  sqrt(yvec(1)**2 + yvec(2)**2)

    xhat = xvec / xrange
    yhat = yvec / yrange

    ! need to move the vertices out slightly along the coordinate vector axes
    ! this is to try to reduce the significance of gyro-orbit losses in ERO for particles injected close to the simulation boundaries
    ! it should not affect the toroidal direction significantly 
    !write(0,*) 'shift_opt:',shift_opt

    if (shift_opt.eq.0) then 
       ! so not shift the corner vertices
       xshift = 0.0
       yshift = 0.0
    elseif (shift_opt.eq.1) then 
       ! fixed shift of 1% for each corner
       xshift = 0.01 * xhat
       yshift = 0.01 * yhat
       !write(0,*) 'xhat:',xhat
       !write(0,*) '0.01*xhat:',0.01*xhat
       !write(0,*) 'xshift:',0.01*xshift

    elseif (shift_opt.eq.2) then 
       ! fixed shift of 5% for each corner
       xshift = 0.05 * xvec
       yshift = 0.05 * yvec
    endif

    ! shift the vertices as follows
    ! 1 moves -x and -y
    ! 2 moves +x and -y
    ! 3 moves +x and +y
    ! 4 moves -x and +y

    vert(1,:) = vert(1,:) - xshift - yshift
    vert(2,:) = vert(2,:) + xshift - yshift
    vert(3,:) = vert(3,:) + xshift + yshift
    vert(4,:) = vert(4,:) - xshift + yshift

    ! After shifting vertices it is now necessary to recalculate xrange and yrange (which requires revising xvec,yvec
    ! The normalized vectors xhat and yhat remain the same 

    xvec = vert(2,:)-vert(1,:)
    yvec = vert(4,:)-vert(1,:)

    xrange =  sqrt(xvec(1)**2 + xvec(2)**2)
    yrange =  sqrt(yvec(1)**2 + yvec(2)**2)


    !write(0,*) 'offset override:',override_offsets,roffset,zoffset
    if (override_offsets.eq.1) then 
       ! calculate the offset automatically as the point midway between vertex 1 and vertex 4.
       ! The ERO geometry should start with 0.0 at the minimum R value with the Z axis centered on this side
       
       roffset = (vert(1,1)+vert(4,1))/2.0
       zoffset = (vert(1,2)+vert(4,2))/2.0
       !write(0,*) 'offsets revised:',override_offsets,roffset,zoffset

    endif



    tmp_offset(1) = -roffset + vert(1,1)
    tmp_offset(2) = -zoffset + vert(1,2)

    xyoffset(1) = dot_product(tmp_offset,xhat)
    xyoffset(2) = dot_product(tmp_offset,yhat)

    write(0,'(a,10(1x,g12.5))') 'xvec :',xvec(1),xvec(2)
    write(0,'(a,10(1x,g12.5))') 'yvec :',yvec(1),yvec(2)
    write(0,'(a,10(1x,g12.5))') 'xhat :',xhat(1),xhat(2)
    write(0,'(a,10(1x,g12.5))') 'yhat :',yhat(1),yhat(2)
    write(0,'(a,10(1x,g12.5))') 'xshif:',xshift(1),xshift(2)
    write(0,'(a,10(1x,g12.5))') 'yshif:',yshift(1),yshift(2)
    write(0,'(a,10(1x,g12.5))') 'vert1:',vert(1,1),vert(1,2)
    write(0,'(a,10(1x,g12.5))') 'vert2:',vert(2,1),vert(2,2)
    write(0,'(a,10(1x,g12.5))') 'vert3:',vert(3,1),vert(3,2)
    write(0,'(a,10(1x,g12.5))') 'vert4:',vert(4,1),vert(4,2)
    write(0,'(a,10(1x,g12.5))') 'tmpof:',tmp_offset(1),tmp_offset(2)
    write(0,'(a,10(1x,g12.5))') 'xyoff:',xyoffset(1),xyoffset(2)


    if (delx.gt.0.0.and.dely.gt.0.0) then 
       ! use delx,dely
       ! data needs to be provided on the boundaries so delx needs to be adjusted so it works out evenly
       nx = int(xrange/delx)
       delx = xrange / real(nx)     

       ny = int(yrange/dely)
       dely = yrange / real(ny)     

    else
       ! use nx,ny

       delx = xrange / real(nx)     
       dely = yrange / real(ny)     

    endif

    rxstep = delx * xhat(1)
    zxstep = delx * xhat(2)

    rystep = dely * yhat(1)
    zystep = dely * yhat(2)

    write(0,'(a,10(1x,g12.5))') 'Steps :',rxstep,zxstep,rystep,zystep
    write(0,'(a,10(1x,g12.5))') 'Offset:', xyoffset(1),xyoffset(2)
    write(0,'(a,10(1x,g12.5))') 'Deltas:', delx,dely
    write(0,'(a,10(1x,g12.5))') 'Vert  :', vert(1,1),vert(1,2)


    ! Allocate storage to hold the ERO output

    call allocate_array(ero_plasma_out,0,nx,0,ny,1,17,'ERO_PLASMA',ierr)

    if (ierr.ne.0) then 
       call errmsg('CALC_ERO_PLASMA: Error allocating ero_plasma_out',ierr)
       stop 'Problem allocating ero_plasma_out'
    endif


    ! loop through the mesh and find the plasma conditions

    do ix = 0,nx

       do iy = 0,ny

          xt = ix * delx + xyoffset(1)
          yt = iy * dely + xyoffset(2)

          rt = ix * rxstep + iy * rystep + vert(1,1)
          zt = ix * zxstep + iy * zystep + vert(1,2)

          !write(0,'(a,10(1x,g12.5))') 'Point:',xt,yt,rt,zt

          call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,psin,btot,br,bz,bt,ngrad,tegrad,tigrad,ierr)

          ero_plasma_out(ix,iy,1) = xt
          ero_plasma_out(ix,iy,2) = yt

          ero_plasma_out(ix,iy,3) = rt
          ero_plasma_out(ix,iy,4) = zt

          ero_plasma_out(ix,iy,5) = ne
          ero_plasma_out(ix,iy,6) = te
          ero_plasma_out(ix,iy,7) = ti
          ero_plasma_out(ix,iy,8) = vb
          ero_plasma_out(ix,iy,9) = ef


          ! Note: Magnetic field vectors probably should be projected onto the 
          !       xhat and yhat coordinate vectors of the mapped ero simulation volume (?)
          !       As long as x->R and y->Z there is no problem but if the coordinates are 
          !       organized otherwise then I am not sure how the magnetic field vectors 
          !       would be interpreted inside ERO.
          !       Btot and Btoroidal would be unchanged but the others may need to be remapped
          !
          ero_plasma_out(ix,iy,10) = btot

          ero_plasma_out(ix,iy,13) = bt

          if (remap_bfield.eq.1) then 
             bfield(1) = br
             bfield(2) = bz

             ! calculate projections of poloidal b-field onto X and Y coordinate vectors
             bx = dot_product(bfield,xhat)
             by = dot_product(bfield,yhat)
             ero_plasma_out(ix,iy,11) = bx
             ero_plasma_out(ix,iy,12) = by
          else
             ero_plasma_out(ix,iy,11) = br
             ero_plasma_out(ix,iy,12) = bz
          endif


          ! gradient data

          ero_plasma_out(ix,iy,14) = ngrad 
          ero_plasma_out(ix,iy,15) = tegrad 
          ero_plasma_out(ix,iy,16) = tigrad 


          ! extrapolation/error flag
          ero_plasma_out(ix,iy,17) = ierr

       end do

    end do


    call output_coordinate_transform('coordinate_transform_data.txt')


  end subroutine calc_ero_plasma



  subroutine calc_ero_surface
    implicit none
    ! This routine calculates the ERO surface based on the DIVIMP wall data

    integer :: nw
    real*8, allocatable :: rw(:),zw(:)
    real*8, allocatable :: surf_flag(:)
    real*8, allocatable :: tmp_array(:)

    integer :: ierr,in_count
    integer :: end_count, surf_ind(2)
    logical :: finished
    integer :: cnt,ind


    integer :: sect
    logical :: found
    real*8 :: rint, zint
    real*8 :: r_int(2),z_int(2)
    real*8 :: r0,thet,thet_extent


    integer :: in,in1,in2,iv,iv2,in_next


    ! output_ero_surface write out the surface shape if required and if wall data is available. 
    !
    ! If the option is specified - write out the 3D ERO limiter shape
    ! - The plasma is defined as being the in the R,Z plane while the limiter extends into 3D along the T axis
    ! 
    ! - In ero .. the limiter shape is specied by Tpx, Tpy, Tpz .
    !   tentative - Tpx toroidal ... Tpz ... Z .... Tpy ... Radial
    ! - X and Z define the surface
    !
    ! Include an option to put in toroidal curvature or not ... 
    !
    ! Need definitions for X0, Y0, Z0, Alp0, RTor, RZ0
    ! 
    ! What do RegInCood and RegTorPRf do? Only matlab? 
    !

    ! X=T
    ! Y=R
    ! Z=Z

    ! initialization
    if (allocated(rw)) deallocate(rw)
    if (allocated(zw)) deallocate(zw)
    if (allocated(surf_flag)) deallocate(surf_flag)



    ! determine segments composing surface in poloidal plane
    ! divide into sufficient number of pieces with straight lines between
    ! extend toroidal surface +/- 25 points in each direction


    if (eroshape_opt.gt.0) then 

       call get_wall_data(nw,rw,zw)

       if (nw.le.0.or.(.not.allocated(rw)).or.(.not.allocated(zw))) then 
          call errmsg('ERROR: CALC_ERO_SURFACE:',' WALL DATA NOT LOADED FROM OEDGE')
          stop 'CALC_ERO_SURFACE: Wall data not loaded'
          !write(0,*) 'Wall:',nw,allocated(rw),allocated(zw)
       endif
       
       !do in = 1,nw
       !   write(6,'(a,i8,2g12.5)') 'WALL:',in,rw(in),zw(in)
       !end do

       allocate(surf_flag(nw),stat=ierr)

       in_count = 0
       surf_flag = 0

       do in = 1,nw
          !write(0,'(a,l4,2i8,20(1x,g12.5))') 'VERT:',inpoly(rw(in),zw(in),nvert,vert(:,1),vert(:,2)),in,nvert,rw(in),zw(in),vert(:,1),vert(:,2)
          if (inpoly(rw(in),zw(in),nvert,vert(:,1),vert(:,2))) then
             surf_flag(in) = 1
             in_count = in_count+1
          endif
       end do

       if (in_count.eq.0) then 
          ! there are no wall vertices inside the ERO sample volume
          ! need to check every wall section for intersections with the ERO sample volume - there should be 2 such intersections if the volume is close to the wall. 
          ! this situation should be uncommon if it occurs at all ... issue error message and stop for now

          call errmsg('calc_ero_surface','no wall sections found in the ERO volume')
          stop 'Calc_ero_surface:no wall points in volume'
       endif


       ! check surf_flag to see if it identifies one contiguous section of surface

       end_count = 0
       surf_ind = 0


       ! find start element

       do in = 1,nw
          if (in.eq.nw) then
             in_next = 1
          else
             in_next = in+1
          endif

          if (surf_flag(in).ne.surf_flag(in_next)) then 
             if (surf_flag(in).eq.1) then 
                surf_ind(2) = in
             else
                surf_ind(1) = in_next
             endif
             end_count = end_count+1
          endif

       end do

       !WRITE(0,*) 'surf_ind:',surf_ind(1),surf_ind(2)


       if (end_count.ne.2) then 
          call errmsg('CALC_ERO_SHAPE: More than 2 surface intersections in ERO volume',end_count)
          stop 'CALC_ERO_SHAPE: More than 2 surface intersections'
       endif


       ! ok - some wall points in volume ... find end-point intersections


       ! cycle over two end points

       do in = 1,2

          if (in.eq.1) then   ! start point

             in1 = surf_ind(in)
             if (in1.eq.1) then 
                in2 = nw
             else
                in2 = in1-1
             endif

          elseif (in.eq.2) then ! end point   
             in1 = surf_ind(in)
             if (in1.eq.nw) then 
                in2 = 1
             else
                in2 = in1+1
             endif

          endif

          sect = 0
          found = .false.

          do iv = 1,nvert

             if (iv.eq.nvert) then 
                iv2=1
             else
                iv2=iv+1
             endif

             !write(0,*) 'in1,in2,iv,iv2:',in1,in2,iv,iv2
             call intsect2dp(vert(iv,1),vert(iv,2),vert(iv2,1),vert(iv2,2),rw(in1),zw(in1),rw(in2),zw(in2),rint,zint,sect)

             if (sect.eq.1) then 

                if (found) then 
                   call errmsg('CALC_ERO_SHAPE:',' Error - Multiple intersections found for same segment')
                else
                   r_int(in) = rint
                   z_int(in) = zint
                   found = .true.
                endif
             endif

          end do

          if (.not.found) then 
             call errmsg('CALC_ERO_SHAPE:','Intersection not found on expected element')
          endif

       end do


       !write(0,*) 'Ints:',r_int,z_int

       ! now we can start calculating the ERO surface - start with only using intersection points and element corners

       n_pol = in_count + 2
       n_tor = 2*tor_points + 1

       !write(0,*) 'n_pol,tor:', n_pol,n_tor

       if (allocated(tpx)) deallocate(tpx)
       if (allocated(tpx)) deallocate(tpy)
       if (allocated(tpx)) deallocate(tpz)


       allocate(tpx(-tor_points:tor_points,n_pol),stat=ierr)  ! T
       allocate(tpy(-tor_points:tor_points,n_pol),stat=ierr)  ! R
       allocate(tpz(-tor_points:tor_points,n_pol),stat=ierr)  ! Z


       ! assign first intersection on mid-point toroidal element 



       cnt = 1
       tpy(0,cnt) = r_int(1)
       tpz(0,cnt) = z_int(1)
       tpx(0,cnt) = 0.0

       ind = surf_ind(1)
       finished = .false.

       do while (.not.finished)

          if (ind.eq.surf_ind(2)) then
             finished = .true.
          endif

          cnt = cnt+1
          tpy(0,cnt) = rw(ind)
          tpz(0,cnt) = zw(ind)
          tpx(0,cnt) = 0.0

          if (ind.eq.nw) then
             ind = 1
          else
             ind = ind+1
          endif

       end do

       ! add last intersection point
       cnt = cnt+1
       tpy(0,cnt) = r_int(2)
       tpz(0,cnt) = z_int(2)
       tpx(0,cnt) = 0.0

       ! ERO requires that the toroidal (tpx) and poloidal (tpz) arrays be in ascending order 
       ! X happens automatically below ... need to check Z and reverse if needed

       if (tpz(0,1) .gt. tpz(0,n_pol)) then 
          ! reorder arrays

          allocate(tmp_array(n_pol),stat=ierr)

          do in = n_pol,1,-1
             tmp_array(n_pol-in+1) = tpz(0,in)
          end do
          tpz(0,:) = tmp_array(:)

          do in = n_pol,1,-1
             tmp_array(n_pol-in+1) = tpy(0,in)
          end do
          tpy(0,:) = tmp_array(:)

          deallocate(tmp_array)

       endif


       ! depending on the shape option - calculate the toroidal projections 


       do cnt = 1,n_pol
          r0 = tpy(0,cnt)
          thet_extent = asin(tor_extent/r0)

          do in = -tor_points,tor_points

             if (eroshape_opt.eq.1) then   !cylindrical .. flat surface extended in toroidal direction
                tpx(in,cnt) = tor_extent * in/tor_points
                tpy(in,cnt) = tpy(0,cnt)
                tpz(in,cnt) = tpz(0,cnt)
             elseif (eroshape_opt.eq.2) then !toroidal ... uses R (tpy) for calculating curvature .. Z is constant 
                tpz(in,cnt) = tpz(0,cnt)

                thet = thet_extent * in/tor_points

                tpy(in,cnt) = r0 * cos(thet)
                tpx(in,cnt) = r0 * sin(thet)

             endif
          end do
       end do

    endif 


    ! Apply coordinate offsets if specified

    do in = -tor_points,tor_points
       do cnt = 1,n_pol
          tpy(in,cnt) = tpy(in,cnt) - roffset
          tpz(in,cnt) = tpz(in,cnt) - zoffset
       end do 
    end do

    ! Remap the offset limiter coordinates into the same coordinate space as the plasma by using the new coordinate basis vectors xhat,yhat
    ! Note that TPx which is the TOROIDAL coordinate is not affected since the coordinate transformation should always be in the poloidal plane

    do in = -tor_points,tor_points
       do cnt = 1,n_pol
          tpy(in,cnt) = tpy(in,cnt)*xhat(1)+tpz(in,cnt)*xhat(2)
          tpz(in,cnt) = tpy(in,cnt)*yhat(1)+tpz(in,cnt)*yhat(2)
       end do 
    end do
    



    !
    ! Set extra parameters
    !
    ! Note: ERO only appears to use Z0 and RTor
    !



    ! clean up

    !write(0,*) 'Deallocating:'

    if (allocated(rw)) deallocate(rw)
    if (allocated(zw)) deallocate(zw)
    if (allocated(surf_flag)) deallocate(surf_flag)

    !write(0,*) 'Deallocated locals:'


    ! calculate max and min of the surface arrays
    xmin = minval(tpx)
    xmax = maxval(tpx)
    ymin = minval(tpy)
    ymax = maxval(tpy)
    zmin = minval(tpz)
    zmax = maxval(tpz)


    ! offsets
    x0 = 0.0
    y0 = 0.0
    z0 = 0.0


    !
    ! Actual midpoint of limiter
    !
    !x0 = (xmax+xmin)/2.0
    !y0 = (ymax+ymin)/2.0
    !z0 = (zmax+zmin)/2.0


    ! midpoint of limiter on its surface
    ! 
    !x0 = tpx(0,npol/2) 
    !y0 = tpy(0,npol/2) 
    !z0 = tpz(0,npol/2) 
    !
    
    rz0 = y0
    
    alp0 = 0.1  ! not sure what this is

    if (eroshape_opt.eq.1) then !cylindrical
       rtor = 1.0e12
    elseif (eroshape_opt.eq.2) then ! toroidal
       rtor = r0
       !rtor = (xmax+xmin)/2.0
       !rtor = x0
    elseif (eroshape_opt.eq.3) then ! toroidal - with Rtor changing sign depending on orientation
       ! this is a test .. may not work .. will only work for coordinate transforms that keep the 1:1 Z,Y->R,Z mapping with the only 
       ! change being the orientation of the "Z" axis relative to R. This will basically set rtor positive for an inner target vertical 
       ! case and negative for an outer target vertical case. Horizontal cases will need more work since the field curvature is across 
       ! the surface - cylindrical approximation might be best for these cases for now. 
       rtor = r0 * xhat(1)
       !rtor = (xmax+xmin)/2.0
       !rtor = x0
    endif


    return

  end subroutine calc_ero_surface

  subroutine output_ero_surface_m
    implicit none

    integer :: outunit,it,ip,ierr

    !
    ! File name and geometry options are read from the erobg input file.
    !
    ! This code writes out a matlab file containing the limiter shape for use in ERO
    !
    ! File name is a parameter but the default is erodiv_shape.m
    ! Coodrdinates in the R,Z plane should match the plasma ... toroidal coordinates can either be flat or curved
    ! based on the toroidal curvature (assumed to have a center at R=0). 
    ! 

    ! sample region is defined by vert(4,2) ... the wall is loaded from the grid file ... intersection points are calculated 
    ! the section of wall within the sample volume is used to define the surface. All turning points are included. Additional 
    ! points may be included if needed. 


    !write(0,*) 'eroshape_opt:',eroshape_opt

    if (eroshape_opt.gt.0) then 
       

       ! write out the ERO limiter shape into a matlab formatted file
       ! Parameters
       !
       ! X0, Y0, Z0, RZ0 (RC0 in ERO), Alp0, RTor
       ! RTor is a toroidal radius quantity - set to 0.0 for cylindrical I think
       !
       ! Format is each profile row on one line for TPx, TPy, TPz .. case matters
       ! Also ... create some additional parameters and catenate them to a parameter template - use jetlim as the template
       !

    call find_free_unit_number(outunit)

    open(outunit,file=trim(eroshape_out_name),form='formatted',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('WRITE_ERO_SHAPE: Problem opening output file:'//trim(eroshape_out_name),ierr)
       stop 'write_ero_plasma: could not open output file'
    endif


    write(outunit,'(a)') '% Data written in MATLAB Compatible format with identifier tags from ERO'
    write(outunit,'(a,2(1x,i8))') '%EROSHAPE:',n_tor,n_pol
    write(outunit,'(a,2i8)') '%  X0 Y0 Z0 RZ0  Alp0  RTor     TPx   TPy    TPz'
    write(outunit,'(a,2i8)') '%  Array data is written out one poloidal profile at a time'
        write(outunit,'(a,2i8)') '%  ERODIV Limiter Shape - dimensions converted to mm from m'


        ! X0
        write(outunit,'(a,g14.6,a)') 'X0 = ',x0*dist_conversion,';'


        ! Y0
        write(outunit,'(a,g14.6,a)') 'Y0 = ',y0*dist_conversion,';'


        ! Z0
        write(outunit,'(a,g14.6,a)') 'Z0 = ',z0*dist_conversion,';'


        ! RZ0
        write(outunit,'(a,g14.6,a)') 'RZ0 = ',rz0*dist_conversion,';'


        ! Alp0
        write(outunit,'(a,g14.6,a)') 'Alp0 = ',alp0,';'


        !RTor
        write(outunit,'(a,g14.6,a)') 'RTor = ',rtor*dist_conversion,';'


        !RegInCood(X0, Y0, Z0, RZ0, Alp0, RTor);
        write(outunit,'(a)') 'RegInCood(X0, Y0, Z0, RZ0, Alp0, RTor);'

        if (n_tor.gt.200) then 
           call errmsg('ERROR: Output_ero_surface : Number of toroidal slices exceeds format statement (200)',n_tor)
           stop 'output_ero_surface: too many toroidal slices'
        endif

        
        !TPx
        write(outunit,'(a)') 'TPx = ['
        do ip = 1,n_pol
           write(outunit,'(200(1x,g15.6))')  (tpx(it,ip)*dist_conversion, it = -tor_points,tor_points)
        end do
        write(outunit,'(a)') ' ];'

        !Tpy

        write(outunit,'(a)') 'TPy = ['
        do ip = 1,n_pol
           write(outunit,'(200(1x,g15.6))')  (tpy(it,ip)*dist_conversion, it = -tor_points,tor_points)
        end do
        write(outunit,'(a)') ' ];'


        !Tpz

        write(outunit,'(a)') 'TPz = ['
        do ip = 1,n_pol
           write(outunit,'(200(1x,g15.6))')  (tpz(it,ip)*dist_conversion, it = -tor_points,tor_points)
        end do
        write(outunit,'(a)') ' ];'

        !RegTorPrf(TPx, TPz, TPy);
        write(outunit,'(a)') 'RegTorPrf(TPx, TPz, TPy);'

        
        ! close 

        close(outunit)

    endif


    write(0,*) 'Surface out: before cleanup'


    call ero_surface_cleanup

    write(0,*) 'Surface out: after cleanup'

  end subroutine output_ero_surface_m



  subroutine output_ero_surface
    implicit none

    integer :: outunit,it,ip,ierr

    !
    ! File name and geometry options are read from the erobg input file.
    !
    ! This code writes out a matlab file containing the limiter shape for use in ERO
    !
    ! File name is a parameter but the default is erodiv_shape.m
    ! Coodrdinates in the R,Z plane should match the plasma ... toroidal coordinates can either be flat or curved
    ! based on the toroidal curvature (assumed to have a center at R=0). 
    ! 

    ! sample region is defined by vert(4,2) ... the wall is loaded from the grid file ... intersection points are calculated 
    ! the section of wall within the sample volume is used to define the surface. All turning points are included. Additional 
    ! points may be included if needed. 


    !write(0,*) 'eroshape_opt:',eroshape_opt

    if (eroshape_opt.gt.0) then 
       

       ! write out the ERO limiter shape into a matlab formatted file
       ! Parameters
       !
       ! X0, Y0, Z0, RZ0 (RC0 in ERO), Alp0, RTor
       ! RTor is a toroidal radius quantity - set to 0.0 for cylindrical I think
       !
       ! Format is each profile row on one line for TPx, TPy, TPz .. case matters
       ! Also ... create some additional parameters and catenate them to a parameter template - use jetlim as the template
       !

    call find_free_unit_number(outunit)

    open(outunit,file=trim(eroshape_out_name),form='formatted',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('WRITE_ERO_SHAPE: Problem opening output file:'//trim(eroshape_out_name),ierr)
       stop 'write_ero_plasma: could not open output file'
    endif


    write(outunit,'(a)') '# Data written in MATLAB Compatible format with identifier tags from ERO'
    write(outunit,'(a,2(1x,i8))') '#EROSHAPE:',n_tor,n_pol
    write(outunit,'(a,2i8)') '#  X0 Y0 Z0 RZ0  Alp0  RTor     TPx   TPy    TPz'
    write(outunit,'(a,2i8)') '#  Array data is written out one poloidal profile at a time'
        write(outunit,'(a,2i8)') '#  ERODIV Limiter Shape - dimensions converted to mm from m'


        ! X0
        write(outunit,'(a,g14.6,a)') '#X0 ',x0*dist_conversion


        ! Y0
        write(outunit,'(a,g14.6,a)') '#Y0 ',y0*dist_conversion


        ! Z0
        write(outunit,'(a,g14.6,a)') '#Z0 ',z0*dist_conversion


        ! RZ0
        write(outunit,'(a,g14.6,a)') '#RZ0 ',rz0*dist_conversion


        ! Alp0
        write(outunit,'(a,g14.6,a)') '#Alp0 ',alp0


        !RTor
        write(outunit,'(a,g14.6,a)') '#RTor ',rtor*dist_conversion


        !RegInCood(X0, Y0, Z0, RZ0, Alp0, RTor);
        !write(outunit,'(a)') 'RegInCood(X0, Y0, Z0, RZ0, Alp0, RTor);'

        if (n_tor.gt.200) then 
           call errmsg('ERROR: Output_ero_surface : Number of toroidal slices exceeds format statement (200)',n_tor)
           stop 'output_ero_surface: too many toroidal slices'
        endif

        
        !TPx
        write(outunit,'(a,i8)') '#TPx ',n_pol
        do ip = 1,n_pol
           write(outunit,'(200(1x,g15.6))')  (tpx(it,ip)*dist_conversion, it = -tor_points,tor_points)
        end do

        !Tpy
        write(outunit,'(a,i8)') '#TPy ',n_pol
        do ip = 1,n_pol
           write(outunit,'(200(1x,g15.6))')  (tpy(it,ip)*dist_conversion, it = -tor_points,tor_points)
        end do

        !Tpz

        write(outunit,'(a,i8)') '#TPz ',n_pol
        do ip = 1,n_pol
           write(outunit,'(200(1x,g15.6))')  (tpz(it,ip)*dist_conversion, it = -tor_points,tor_points)
        end do

        !RegTorPrf(TPx, TPz, TPy);
        !write(outunit,'(a)') 'RegTorPrf(TPx, TPz, TPy);'

        
        ! close 

        close(outunit)

    endif


    !write(0,*) 'Surface out: before cleanup'


    call ero_surface_cleanup

    !write(0,*) 'Surface out: after cleanup'

  end subroutine output_ero_surface




  subroutine write_ero_par_file
    implicit none

    integer :: outunit,rc,ierr
    character*512 :: cmd

    ! this routine contatenates ERODIV options to the end of a ERO template parameter file


    call find_free_unit_number(outunit)

    open(outunit,file='par_erosdiv_add',form='formatted',iostat=ierr)


        write(outunit,'(a)') '# ERODIV input parameters '
        write(outunit,'(a)') 'erodiv_bm_pol_prf'
        write(outunit,'(i8)') n_pol
        write(outunit,'(a)') 'erodiv_bm_tor_prf'
        write(outunit,'(i8)') n_tor
        write(outunit,'(a)') 'erodiv_ranges'
        write(outunit,'(6(1x,g12.5))') xmin*dist_conversion,xmax*dist_conversion,ymin*dist_conversion,ymax*dist_conversion,zmin*dist_conversion,zmax*dist_conversion
        write(outunit,'(a)') 'erodiv_pl_dy'
        write(outunit,'(6(1x,g12.5))') delx
        write(outunit,'(a)') 'erodiv_pl_dz'
        write(outunit,'(6(1x,g12.5))') dely
        write(outunit,'(a)') 'erodiv_pl_fn'
        write(outunit,'(a)') trim(erobg_out_name)
        write(outunit,'(a)') 'erodiv_surface_fn'
        write(outunit,'(a)') trim(eroshape_out_name)
        

    close(outunit)    

    cmd = 'cat par_erodiv_tempate par_erodiv_add > par_erodiv'

    call cissue(cmd,rc)

    if (rc.ne.0) then 
       call errmsg('ERROR: write_ero_par_file: problem executing command',rc)
    endif


  end subroutine write_ero_par_file

  subroutine output_ero_plasma
    implicit none
    ! Added this routine to allow choosing between output file format as an option without
    ! re-writing code
    
    if (trim(out_format).eq.'matlab') then 
       call output_ero_plasma_m
    else
       call output_ero_plasma
    endif

    return
  end subroutine output_ero_plasma


  subroutine output_ero_plasma_m
    implicit none

    ! file name to be saved as is part of the ero input block
    ! data is listed along each x row first

    integer :: outunit,ierr
    integer :: ix,iy,in

    real*8 :: conv_factor

    real*8,allocatable :: grad_xy(:,:)

    character*100 :: data_out_format

    integer, parameter :: write_rz = 1

    data_out_format='(1000(1x,e12.5))'


    call find_free_unit_number(outunit)

    open(outunit,file=trim(erobg_out_name),form='formatted',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('WRITE_ERO_PLASMA: Problem opening output file:'//trim(erobg_out_name),ierr)
       stop 'write_ero_plasma: could not open output file'
    endif


    write(outunit,'(a,2i8)') '%BG DOC  ix iy x y r z ne te ti vpara epara btot br bz btor ngrad tegrad tigrad ierr->(on or off grid)'
    write(outunit,'(a,2i8)') '%EROBG',nx+1,ny+1
    write(outunit,'(a)') '% Data written in MATLAB Compatible format with identifier tags from ERO'

    !do ix = 0,nx
    !   do iy = 0,ny
    !       write(outunit,'(2i8,20(1x,g18.8))') ix,iy,(ero_plasma_out(ix,iy,in),in=1,17)
    !   end do 
    !end do 

    ! ERO plasma files from JET contain the following data ... each line contains space separated floats
    ! Also - the data assume a rectilinear grid ... not sure how to do the coordinate mapping yet
    ! Also need to calculate rectilinear R,Z gradients of Te for some reason

    ! NR = 
    write(outunit,'(a,i8,a)') 'NR = ',nx+1,';'

    ! NZ = 
    write(outunit,'(a,i8,a)') 'NZ = ',ny+1,';'


    ! RC = [ ]
    ! Assume constant axes - write out X data for now
    write(outunit,'(a)') 'RC = ['
    write(outunit,trim(data_out_format)) (ero_plasma_out(ix,0,1)*dist_conversion,ix=0,nx)
    write(outunit,'(a)') '];'

    ! ZC = [ ]
    ! Assume constant axes - write out Y data for now
    write(outunit,'(a)') 'ZC = ['
    write(outunit,trim(data_out_format)) (ero_plasma_out(0,iy,2)*dist_conversion,iy=0,ny)
    write(outunit,'(a)') '];'



    if (write_rz.eq.0) then 
       ! write out the arrays indexed as R then Z
       write(0,*) 'Writing arrays R down the column and Z across the row'

       ! ne = [ ] convert to cm-3
       ! 
       conv_factor = 1.0d-6

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'ne = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,5)*conv_factor,iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! te = [ ]  - ok as eV
       !  
       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'te = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,6),iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! Grad_te = [ ]   - convert to eV/mm
       conv_factor = dist_conversion

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Grad_te = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,15)/conv_factor,iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! Grad_te_RC = [ ] (?)

       allocate(grad_xy(0:nx,0:ny)) 

       grad_xy = 0.0
       ! grad Te on X - convert to eV/mm
       conv_factor = dist_conversion

       call calc_grad(grad_xy,1,6)

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Grad_te_RC = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (grad_xy(ix,iy)/conv_factor,iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! Grad_te_ZC = [ ] (?)

       grad_xy = 0.0
       ! grad Te on Y - convert to eV/mm  
       conv_factor = dist_conversion

       call calc_grad(grad_xy,2,6)

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Grad_te_ZC = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (grad_xy(ix,iy)/conv_factor,iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       if (allocated(grad_xy)) deallocate(grad_xy)

       ! ti = [  - ok as eV

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'ti = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,7),iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! vpar = [   - convert to cm/s  (?) - numbers in sample file seem too small 

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       conv_factor = 100.0

       write(outunit,'(a)') 'vpar = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,8)*conv_factor,iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! Epar = [  - convert to V/mm (?) 
       conv_factor = dist_conversion

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Epar = ['   
       do ix = 0,nx
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,9)*conv_factor,iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! BFld_Tor = [ - tesla should be ok

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'BFld_tor = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) ((ero_plasma_out(ix,iy,10)*ero_plasma_out(ix,iy,13)),iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! BFld_RC = [ - tesla should be ok

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'BFld_RC = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) ((ero_plasma_out(ix,iy,10)*ero_plasma_out(ix,iy,11)),iy=0,ny)
       end do
       write(outunit,'(a)') '];'

       ! BFld_ZC = [ - tesla should be ok

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'BFld_ZC = ['
       do ix = 0,nx
          write(outunit,trim(data_out_format)) ((ero_plasma_out(ix,iy,10)*ero_plasma_out(ix,iy,12)),iy=0,ny)
       end do
       write(outunit,'(a)') '];'

    elseif (write_rz.eq.1) then 
       ! write out the arrays indexed by Z then R

       write(0,*) 'Writing arrays Z down the column and R across the row'

       ! ne = [ ] convert to cm-3
       ! 
       conv_factor = 1.0d-6

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'ne = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,5)*conv_factor,ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! te = [ ]  - ok as eV
       !  
       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'te = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,6),ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! Grad_te = [ ]   - convert to eV/mm
       conv_factor = dist_conversion

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Grad_te = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,15)/conv_factor,ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! Grad_te_RC = [ ] (?)

       allocate(grad_xy(0:nx,0:ny)) 

       grad_xy = 0.0
       ! grad Te on X - convert to eV/mm
       conv_factor = dist_conversion

       call calc_grad(grad_xy,1,6)

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Grad_te_RC = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (grad_xy(ix,iy)/conv_factor,ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! Grad_te_ZC = [ ] (?)

       grad_xy = 0.0
       ! grad Te on Y - convert to eV/mm  
       conv_factor = dist_conversion

       call calc_grad(grad_xy,2,6)

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Grad_te_ZC = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (grad_xy(ix,iy)/conv_factor,ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       if (allocated(grad_xy)) deallocate(grad_xy)

       ! ti = [  - ok as eV

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'ti = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,7),ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! vpar = [   - convert to cm/s  (?) - numbers in sample file seem too small 

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       conv_factor = 100.0

       write(outunit,'(a)') 'vpar = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,8)*conv_factor,ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! Epar = [  - convert to V/mm (?) 
       conv_factor = dist_conversion

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'Epar = ['   
       do iy = 0,ny
          write(outunit,trim(data_out_format)) (ero_plasma_out(ix,iy,9)*conv_factor,ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! BFld_Tor = [ - tesla should be ok

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'BFld_Tor = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) ((ero_plasma_out(ix,iy,10)*ero_plasma_out(ix,iy,13)),ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! BFld_RC = [ - tesla should be ok

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'BFld_RC = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) ((ero_plasma_out(ix,iy,10)*ero_plasma_out(ix,iy,11)),ix=0,nx)
       end do
       write(outunit,'(a)') '];'

       ! BFld_ZC = [ - tesla should be ok

       write(outunit,'(a)') '%'
       write(outunit,'(a)') '% ---------------------------------------------------------------------------------------- '
       write(outunit,'(a)') '%'

       write(outunit,'(a)') 'BFld_ZC = ['
       do iy = 0,ny
          write(outunit,trim(data_out_format)) ((ero_plasma_out(ix,iy,10)*ero_plasma_out(ix,iy,12)),ix=0,nx)
       end do
       write(outunit,'(a)') '];'


    endif

    call ero_plasma_cleanup


  end subroutine output_ero_plasma_m



  subroutine output_ero_plasma_ascii
    implicit none

    ! file name to be saved as is part of the ero input block
    ! data is listed along each x row first

    integer :: outunit,ierr
    integer :: ix,iy,in

    call find_free_unit_number(outunit)

    open(outunit,file=trim(erobg_out_name),form='formatted',iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('WRITE_ERO_PLASMA: Problem opening output file:'//trim(erobg_out_name),ierr)
       stop 'write_ero_plasma: could not open output file'
    endif


    write(outunit,'(a,2i8)') '#BG DOC  ix iy x y r z ne te ti vpara epara btot br bz btor ierr->(on or off grid)'
    write(outunit,'(a,2i8)') '#EROBG',nx+1,ny+1

    do ix = 0,nx
       do iy = 0,ny
           write(outunit,'(2i8,20(1x,g18.8))') ix,iy,(ero_plasma_out(ix,iy,in),in=1,17)
       end do 
    end do 

    call ero_plasma_cleanup


  end subroutine output_ero_plasma_ascii


  subroutine calc_grad(grad_arr,axis,quant)

    real*8,allocatable :: grad_arr(:,:) 
    integer :: axis,quant

    integer :: dist,ixe,ixs,iye,iys

    if (.not.allocated(grad_arr)) return


    do ix = 0,nx

       do iy = 0,ny

          if (axis .eq. 1) then
             ! R gradients - X axis
             iye = iy
             iys = iy
             if (ix .eq.0) then 
                ixe = ix+1
                ixs = ix
             elseif (ix.eq.nx) then     
                ixe = ix
                ixs = ix-1
             else
                ixe = ix+1
                ixs = ix-1
             endif
          elseif (axis.eq.2) then 
             ! Z gradients - Y axis
             ixe = ix
             ixs = ix
             if (iy .eq.0) then 
                iye = iy+1
                iys = iy
             elseif (iy.eq.ny) then     
                iye = iy
                iys = iy-1
             else
                iye = iy+1
                iys = iy-1
             endif

          endif

          dist = ero_plasma_out(ixe,iye,axis) - ero_plasma_out(ixs,iys,axis)

          if (dist.ne.0.0) then 
             grad_arr(ix,iy) = (ero_plasma_out(ixe,iye,quant)-ero_plasma_out(ixs,iys,quant))/dist
          else
             grad_arr(ix,iy) = 0.0
          endif

       end do

    end do

  end subroutine calc_grad




  subroutine ero_plasma_cleanup
    implicit none
    ! remove allocated variables for ERO plasma
    if (allocated(ero_plasma_out)) deallocate(ero_plasma_out)

  end subroutine ero_plasma_cleanup


  subroutine ero_surface_cleanup
    implicit none

    ! delete allocated variables for ERO shape
    if (allocated(tpx)) deallocate(tpx)
    if (allocated(tpy)) deallocate(tpy)
    if (allocated(tpz)) deallocate(tpz)


  end subroutine ero_surface_cleanup


  subroutine output_coordinate_transform(coordinate_filename)
    implicit none
    character*(*) :: coordinate_filename
    integer :: outunit,i,ios


    ! this routine writes out the coordinate transformation data so that 
    ! the particle processing routing can run it in reverse to obtain the 
    ! DIVIMP coordinates from the ERO particle coordinates. 

    call find_free_unit_number(outunit)

    open(outunit,file=trim(coordinate_filename),form='formatted',iostat=ios)
    

    ! Need to write out 
    ! 1) The vertices
    ! 2) The transform vectors
    ! 3) The inverse transform vectors 
    ! 4) The offsets

    write(outunit,'(a,2(1x,g18.8))') '#OFFSETS',roffset,zoffset
    write(outunit,'(a,4(1x,g18.8))') '#RVERTICES',(vert(i,1),i=1,4)
    write(outunit,'(a,4(1x,g18.8))') '#ZVERTICES',(vert(i,2),i=1,4)
    write(outunit,'(a,1(1x,g18.8))') '#TOR_EXTENT',tor_extent
    write(outunit,'(a,4(1x,g18.8))') '#XHAT',(xhat(i),i=1,2)
    write(outunit,'(a,4(1x,g18.8))') '#YHAT',(yhat(i),i=1,2)
    write(outunit,'(a,4(1x,g18.8))') '#XINV',(yhat(i),i=2,1,-1)
    write(outunit,'(a,4(1x,g18.8))') '#YINV',(xhat(i),i=2,1,-1)
    write(outunit,'(a,a)') '#EROPARTFN',' erodiv_ero_part_lost.dat'
    ! Note: DIVIMP reads the file as ero_particle_data.dat ... however, when iterating this will already exist so a different name is needed
    write(outunit,'(a,a)') '#ERODIVPARTFN',' ero_particle_data.out'

    close(outunit)


  end subroutine output_coordinate_transform

  subroutine process_divimp_part_file
    implicit none
    ! this routine reads in the DIVIMP particle track file and makes it compatible with ERO
    ! Note that the toroidal coordinate is set to zero by default to start with. 
    ! Toroidal -> X   Vertical/poloidal -> Y and Radial -> Z

    ! Open the input and output files
    ! Transform coordinates (and units if required) 
    ! Write to output file

    integer :: partin,partout,ios,in
    logical :: done
    integer :: part_count
    real*8 :: absfac, nparticles, scale_fact
    real*8 :: r,z,t,vr,vz,vt,vpara,temp,weight
    integer :: charge

    real*8 :: rprime, zprime, vrprime, vzprime
    real*8 :: tprime,atoms

    real :: rand

    ! Check to see if particle tracking option was on
    if (divpart_opt.eq.0) return


    ! Initialize the random number generator so I can create 
    ! random locations in +/- tor_extent for the particles
    call init_random_seed

    ! open files


    call find_free_unit_number(partin)

    open(partin,file=trim(divpart_in_name),status='old',form='formatted',iostat=ios)

    if (ios.ne.0) then 
       call errmsg('ERROR opening DIVIMP particle input file: Case = ',trim(divpart_in_name))
       return
    endif

    call find_free_unit_number(partout)

    open(partout,file=trim(eropart_out_name),form='formatted',iostat=ios)

    if (ios.ne.0) then 
       call errmsg('ERROR opening DIVIMP-> particle transfer file: Case = ',trim(eropart_out_name))
       return
    endif


    done = .false.

    !
    ! First scan the file and count particles and extract absfac information
    !

    do while (.not.done)

       read(partin,line_form,iostat=ios) buffer

       if (ios.ne.0) then 
          ! need to exit loop if eof reached 
          done = .true. 
          cycle
       endif

       if (buffer(1:1).ne.'#'.and.buffer(1:1).ne.'%') then 
          part_count = part_count+1
       elseif (buffer(1:7).eq.'#ABSFAC') then
          read(buffer(8:),*) absfac,nparticles,scale_fact
          !write(0,'(a,3(1x,g18.8))') 'Scaling: ABSFAC NPART SCALE:',absfac,nparticles,scale_fact
       endif
    end do

    ! Now process the file

    rewind(partin)

    !write(partout,'(a,1x,i12)') '#NPARTICLES',part_count*erodiv_tor_part_mult
    write(partout,'(a,1x,i12)') '#NPARTICLES',part_count*erodiv_tor_part_mult
    write(partout,'(i12)')  part_count*erodiv_tor_part_mult

    done =.false.
    part_count = 0

    do while (.not.done)

       read(partin,line_form,iostat=ios) buffer

       !write(0,'(a,i5,a,a)') 'Re-reading file:',ios,trim(buffer),':'

       if (ios.ne.0) then 
          ! need to exit loop if eof reached 
          done = .true. 
          exit
       endif

       if (buffer(1:1).ne.'#'.and.buffer(1:1).ne.'%') then 

          !# R   Z     CHARGE     VR    VZ     VT    VPARA    TEMP   WEIGHT
          read(buffer,*) r,z,charge,vr,vz,vt,vpara,temp,weight

          ! apply offset and then basis vector transformation
          rprime = r-roffset
          zprime = z-zoffset

          rprime = (rprime * xhat(1) + zprime * xhat(2)) 
          zprime = (rprime * yhat(1) + zprime * yhat(2)) 

          !write(0,'(a,10(1x,g12.5))') 'Transform:',r,z,rprime,zprime,xhat,yhat,roffset,zoffset

          vrprime = vr * xhat(1) + vz * xhat(2)
          vzprime = vr * yhat(1) + vz * yhat(2)

          atoms = weight * scale_fact / real(erodiv_tor_part_mult)

          do in = 1,erodiv_tor_part_mult
             part_count = part_count + 1
             call random_number(rand)
             tprime = (0.5-rand) * tor_extent
             ! T->X  Z -> Y  R-> Z
             write(partout,'(i8,i8,10(1x,e18.8))') part_count,charge,tprime*1000.0,zprime*1000.0,rprime*1000.0, &
                  & vt*100.0,vzprime*100.0,vrprime*100.0,vpara*100.0,temp,atoms
          end do

       endif

    end do

    ! close files
    close(partin)
    close(partout)

  end subroutine process_divimp_part_file


end module ero_plasma
