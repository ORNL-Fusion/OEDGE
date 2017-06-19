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
  integer :: remap_bfield
  character*(line_length) :: erobg_out_name
  !
  ! Grid resolution can be specified by either nx,ny or delx,dely - number of points is rounded off when delx.dely is specified. 
  !
  integer :: nx,ny
  real*8 :: delx,dely

  logical :: divimpbg_loaded = .false.

  ! ero background plasma storage

  real*8,allocatable :: ero_plasma_out(:,:,:)


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

       write(0,'(a,a,a,a)') 'buffer:',trim(buffer),':'

       if (buffer(1:15).eq.'#EROBLOCK') then
          done=.true. 
       elseif (buffer(1:10).eq.'#NUMBLOCKS') then 
          read(buffer(11:),*) n_ero_blocks
          write(0,'(a,i10)') 'nblocks:',n_ero_blocks
       elseif (buffer(1:14).eq.'#DIVPLASMAFILE') then 
          !read(specunit,line_form,iostat=ios) buffer
          
          read(buffer(15:),*) div_plasma_file
          read_plasma_file = .true.
          write(0,'(a,a,a,a)') 'bgp:',trim(div_plasma_file),':'
          
          !if (ios.eq.0) then 
          !   div_plasma_file = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_spec_headers: read plasma file name: code = ',ios)
          !   stop 'Unexpected error reading '//trim(erospec_name)
          !endif

       elseif (buffer(1:11).eq.'#DIVGEOFILE') then 
          read(buffer(12:),*) div_geo_file
          read_geo_file = .true.
          write(0,'(a,a,a,a)') 'grd:',trim(div_geo_file),':'

          !read(specunit,line_form,iostat=ios) buffer

          !if (ios.eq.0) then 
          !   div_geo_file = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_spec_headers: read geo file name: code = ',ios)
          !   stop 'Unexpected error reading '//trim(erospec_name)
          !endif

       endif

    end do

    write(0,*) read_geo_file, read_plasma_file

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
          write(0,'(a,10(1x,g12.5))') 'offset:',roffset,zoffset

       elseif (buffer(1:10).eq.'#RVERTICES') then 
          ! read R vertices
          read(buffer(11:),*) (vert(i,1),i=1,4)
          rv_read = .true.
          write(0,'(a,10(1x,g12.5))') 'rvert:',(vert(i,1),i=1,4)

       elseif (buffer(1:10).eq.'#ZVERTICES') then 
          ! read Z vertices
          read(buffer(11:),*) (vert(i,2),i=1,4)
          zv_read = .true.
          write(0,'(a,10(1x,g12.5))') 'zvert:',(vert(i,2),i=1,4)

       elseif (buffer(1:5).eq.'#NXNY') then 
          ! read nx ny
          read(buffer(6:),*) nx,ny

          write(0,'(a,2i10,10(1x,g12.5))') 'nxny:',nx,ny

       elseif (buffer(1:8).eq.'#DELTAXY') then 
          ! read delx dely
          read(buffer(9:),*) delx,dely
          write(0,'(a,10(1x,g12.5))') 'dxdy:',delx,dely

       elseif (buffer(1:12).eq.'#INTERPOLATE') then 
          ! read interpolation option
          read(buffer(13:),*) interp_opt
          write(0,'(a,2i10,10(1x,g12.5))') 'interp:',interp_opt

       elseif (buffer(1:12).eq.'#EXTRAPOLATE') then 
          ! read interpolation option
          read(buffer(13:),*) extrap_opt
          write(0,'(a,2i10,10(1x,g12.5))') 'extrap:',extrap_opt


       elseif (buffer(1:13).eq.'#REMAP_BFIELD') then 
          ! read interpolation option
          read(buffer(14:),*) remap_bfield

          write(0,'(a,2i10,10(1x,g12.5))') 'remap:',remap_bfield

       elseif (buffer(1:10).eq.'#EROBGFILE') then 
          ! read ERO background plasma output file name
          read(buffer(11:),*) erobg_out_name
          !read(specunit,line_form,iostat=ios) buffer
          erobg_name_read = .true.

          write(0,'(a,a,a,a)') 'erobg:',trim(erobg_out_name),':'

          !if (ios.eq.0) then 
          !   erobg_out_name = trim(buffer)
          !else
          !   call errmsg('ERROR read_ero_plasma_block: read erobg output file name: code = ',ios)
          !   stop 'Unexpected error reading erobg_out_name'
          !endif

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
    real*8 :: xvec(2),yvec(2),xhat(2),yhat(2),bfield(2),xrange,yrange
    real*8 :: rt,zt,xt,yt
    real*8 :: rxstep,zxstep,rystep,zystep

    real*8 :: tmp_offset(2),xyoffset(2)

    real*8 :: ne,te,ti,vb,ef,btot,br,bz,bt,psin
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

    tmp_offset(1) = -roffset + vert(1,1)
    tmp_offset(2) = -zoffset + vert(1,2)

    xyoffset(1) = dot_product(tmp_offset,xhat)
    xyoffset(2) = dot_product(tmp_offset,yhat)

    write(0,'(a,10(1x,g12.5))') 'xvec :',xvec(1),xvec(2)
    write(0,'(a,10(1x,g12.5))') 'yvec :',yvec(1),yvec(2)


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
    write(0,'(a,10(1x,g12.5))') 'xhat  :', xhat(1),xhat(2)
    write(0,'(a,10(1x,g12.5))') 'yhat  :', yhat(1),yhat(2)


    ! Allocate storage to hold the ERO output

    call allocate_array(ero_plasma_out,0,nx,0,ny,1,14,'ERO_PLASMA',ierr)

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

          call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,psin,btot,br,bz,bt,ierr)

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

          ero_plasma_out(ix,iy,14) = ierr

       end do

    end do


  end subroutine calc_ero_plasma



  subroutine output_ero_plasma
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
           write(outunit,'(2i8,20(1x,g18.8))') ix,iy,(ero_plasma_out(ix,iy,in),in=1,14)
       end do 
    end do 

    call ero_plasma_cleanup


  end subroutine output_ero_plasma


  subroutine ero_plasma_cleanup
    implicit none

    if (allocated(ero_plasma_out)) deallocate(ero_plasma_out)


  end subroutine ero_plasma_cleanup



end module ero_plasma
