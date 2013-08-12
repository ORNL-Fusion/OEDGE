module ero_plasma
  use error_handling
  use allocate_arrays
  use oedge_plasma_interface

  integer :: specunit

  integer, parameter :: line_length = 512
  character*6, parameter :: line_form='(a512)'
  character*(line_length) :: buffer


  !
  ! Base specifications of divimp plasma and number of ERO blocks in file
  !
  character*(line_length) :: div_plasma_file, div_geo_file
  integer :: n_ero_blocks


  !
  ! ERO block data
  !
  real*8 :: vert(4,2)
  real*8 :: roffset,zoffset
  integer :: interp_opt,extrap_opt
  character*(line_length) :: erobg_out_name
  !
  ! Grid resolution can be specified by either nx,ny or delx,dely - number of points is rounded off when delx.dely is specified. 
  !
  integer :: nx,ny
  real*8 :: delx,dely

  logical :: divimpbg_loaded = .false.

  ! ero background plasma storage

  real*8,allocatable :: ero_plamsma_out(:,:,:)


  public read_ero_plasma_headers, read_ero_plasma_block, calc_ero_plasma,output_ero_plasma

contains



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

       if (buffer(1:9).eq.'#EROBLOCK') then
          done=.true. 
       elseif (buffer(1:10).eq.'#NUMBLOCKS') then 
          read(buffer(11:),*) n_ero_blocks
       elseif (buffer(1:14).eq.'#DIVPLASMAFILE') then 
          read(specunit,line_form,iostat=ios) buffer
          if (ios.eq.0) then 
             div_plasma_file = trim(buffer)
             read_plasma_file = .true.
          else
             call errmsg('ERROR read_ero_spec_headers: read plasma file name: code = ',ios)
             stop 'Unexpected error reading '//trim(erospec_name)
          endif

       elseif (buffer(1:11).eq.'#DIVGEOFILE') then 
          read(specunit,line_form,iostat=ios) buffer
          if (ios.eq.0) then 
             div_geo_file = trim(buffer)
             read_geo_file = .true.
          else
             call errmsg('ERROR read_ero_spec_headers: read geo file name: code = ',ios)
             stop 'Unexpected error reading '//trim(erospec_name)
          endif

       endif

    end do


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
    integer :: i

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


    do while (.not.done)  

       read(specunit,line_form,iostat=ios) buffer

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
          offsets_read = .true.
       elseif (buffer(1:10).eq.'#RVERTICES') then 
          ! read R vertices
          read(buffer(11:),*) (vert(i,1),i=1,4)
          rv_read = .true.
       elseif (buffer(1:10).eq.'#ZVERTICES') then 
          ! read Z vertices
          read(buffer(11:),*) (vert(i,2),i=1,4)
          zv_read = .true.
       elseif (buffer(1:5).eq.'#NXNY') then 
          ! read nx ny
          read(buffer(6:),*) nx,ny
       elseif (buffer(1:8).eq.'#DELTAXY') then 
          ! read delx dely
          read(buffer(9:),*) delx,dely
       elseif (buffer(1:12).eq.'#INTERPOLATE') then 
          ! read interpolation option
          read(buffer(13:),*) interp_opt
       elseif (buffer(1:12).eq.'#EXTRAPOLATE') then 
          ! read interpolation option
          read(buffer(13:),*) extrap_opt
       elseif (buffer(1:12).eq.'#EROBGFILE') then 
          ! read ERO background plasma output file name
          read(specunit,line_form,iostat=ios) buffer
          if (ios.eq.0) then 
             erobg_out_name = trim(buffer)
             erobg_name_read = .true.
          else
             call errmsg('ERROR read_ero_plasma_block: read erobg output file name: code = ',ios)
             stop 'Unexpected error reading erobg_out_name'
          endif
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




  subroutine calc_ero_plasma
    implicit none

    integer :: ierr
    integer :: errmsg_unit = 0

    real*8 :: xvec(2),yvec(2),xhat(2),yhat(2)
    real*8 :: rt,zt,xt,yt
    real*8 :: rxstep,zxstep,rystep,zystep

    !
    ! Check to see if DIVIMP plasma data has been loaded (if not then load it)
    !

    if (.not.divimpbg_loaded) then 
       call load_oedge_plasma(div_geo_file,div_plasma_file,ierr)

       if (ierr.ne.0) then 
          call errmsg('Calc_ero_plasma: Error loading divimp files:',ierr)
          stop 'DIVIMP files could not be loaded'
       endif

       divimpbg_loaded = .true.

    endif

    ! initialize the oedge plasma loading code

    call set_oedge_plasma_opts(roffset,zoffset,interp_opt,extrap_opt,errmsg_unit)

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

    xvec = vert(2,*)-vert(1,*)
    yvec = vert(4,*)-vert(1,*)

    xrange =  sqrt(xvec(1)**2 + xvec(2)**2)
    yrange =  sqrt(yvec(1)**2 + yvec(2)**2)

    xhat = xvec / xrange
    yhat = yvec / yrange

    offset(1) = roffset - vert(1,1)
    offset(2) = zoffset - vert(1,2)

    xyoffset(1) = dot_product(offset,xhat)
    xyoffset(2) = dot_product(offset,yhat)


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


    ! Allocate storage to hold the ERO output

    call allocate_array(eroplasma,nx,ny,14,'ERO_PLASMA',ierr)

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

          call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

          ero_plasma_out(ix,iy,1) = xt
          ero_plasma_out(ix,iy,2) = yt

          ero_plasma_out(ix,iy,3) = rt
          ero_plasma_out(ix,iy,4) = zt

          ero_plasma_out(ix,iy,5) = ne
          ero_plasma_out(ix,iy,6) = te
          ero_plasma_out(ix,iy,7) = ti
          ero_plasma_out(ix,iy,8) = vb
          ero_plasma_out(ix,iy,9) = ef

          ero_plasma_out(ix,iy,10) = btot
          ero_plasma_out(ix,iy,11) = br
          ero_plasma_out(ix,iy,12) = bz
          ero_plasma_out(ix,iy,13) = bt

          ero_plasma_out(ix,iy,14) = ierr

       end do

    end do


  end subroutine calc_ero_plasma



  subroutine output_ero_plasma
    implicit none

    ! file name to be saved as is part of the ero input block
    ! data is listed along each x row first

    integer :: outunit,ierr
    integer :: ix,iy

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
           write(ounit,'(2i8,20(1x,g18.8))') ix,iy,(ero_plasma_out(ix,iy,in),in=1,14)
       end do 
    end do 

    call ero_plasma_cleanup


  end subroutine output_ero_plasma


  subroutine ero_plasma_cleanup
    implicit none

    if (allocated(ero_plasma_out)) deallocate(ero_plasma_out)


  end subroutine ero_plasma_cleanup



end module ero_plasma
