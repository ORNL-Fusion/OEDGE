module oedge_plasma_interface
  use error_handling

  implicit none

  private

  logical :: debug_code = .false.
  integer :: interpolate_opt
  real*8  :: r_offset,z_offset

  integer :: nrs,mks,irsep,nds,npolyp,nvert

  integer,allocatable :: nks(:)

  !-------------------------------------------------
  !
  ! Geometry information
  !
  ! rs,zs - cell centers
  ! rbnd,zbnd - r,z, coordinates of cell boundary mid-point along field lines
  ! 
  !
  real*8,allocatable :: rs(:,:),zs(:,:),rbnd(:,:),zbnd(:,:)

  ! Connection map for grid contains indices of adjacent cells
  ! ikins - knot index inward
  ! irins - ring index inward
  ! ikouts - knot index outward
  ! irouts - ring index outward
  ! idds(ir,1..2) - index to target elements for each ring - 1 = end of ring, 2=start of ring
  integer,allocatable :: ikins(:,:),irins(:,:),ikouts(:,:),irouts(:,:),idds(:,:)

  ! Index from grid to polygon data
  integer,allocatable :: korpg(:,:)

  ! polygon data 
  ! nvertp - vertex count for cell
  ! rvertp - r coordinates of cell vertices
  ! zvertp - z coordinates of cell vertices
  integer,allocatable :: nvertp(:)
  real*8,allocatable :: rvertp(:,:),zvertp(:,:)

  ! Magnetic field data
  ! btot - magntitude of magnetic field in each cell 
  ! br,bz,bt - unit vector components of the magnetic field direction in each cell
  real*8,allocatable :: btot(:,:),br(:,:),bz(:,:),bt(:,:)


  !-------------------------------------------------
  !
  ! Plasma information
  !
  ! Note: plasma file contains some duplication of geometry information
  !
  ! Plasma quantities
  ! - volume terms
  ! knbs - plasma density (m-3)
  ! ktebs - electron temperature (eV)
  ! ktibs - ion temperature (eV)
  ! kvhs - parallel velocity (m/s) - negative velocities are towards the low S or start of field line -
  !                                - need to check REDEP sign convention
  ! kes - parallel electric field (V/m)
  !
  ! - target quantities
  ! knds  - plasma density (m-3)
  ! kteds - electron temperature (eV)
  ! ktids - ion temperature (eV)
  ! kvds  - parallel flow velocity (m/s)
  ! keds  - electric field (V/m) 
  !

  real*8,allocatable :: knbs(:,:),ktebs(:,:),ktibs(:,:),kvhs(:,:),kes(:,:)
  real*8,allocatable :: knds(:),kteds(:),ktids(:),kvds(:),keds(:)
  real*8,allocatable :: btotd(:),brd(:),bzd(:),btd(:)

  ! The meaning of these will depend on whether cell polygon data or cell center data is being used
  integer :: ik_last,ir_last,iq_last


  public :: get_oedge_plasma,init_oedge_plasma,close_oedge_plasma,find_free_unit_number



contains


  subroutine init_oedge_plasma(filename,r_shift,z_shift,interpolate_option,errmsg_unit,ierr)
    implicit none
    character*(*) :: filename
    integer :: interpolate_option,errmsg_unit
    character*256 :: gridfilename,plasmafilename
    integer :: igridfile,iplasmafile,ierr,ios
    real*8 :: r_shift, z_shift

    ! Set error message unit numbers to standard error only for this application instead of the default standard error and 
    ! standard output. 
    ! errmsg_unit < 0 turns off all error messaging
    ! errmsg_unit = 0 error messages go to standard error
    ! errmsg_unit = 6 error messages go to standard out
    ! errmsg_unit = N error messages go to specified unit number

    call set_errmsg_units(errmsg_unit,-1,-1)

    ! Initialize error return value to 0
    ierr = 0

    ! Assign r_shift, z_shift to r_offset,z_offset internally
    ! This feature allows the coordinate systems between the oedge_plasma_interface code
    ! and the calling routines to be adjusted by a standard offset. Sign conventions remain the
    ! same ... however, the calling routine can then use its own coordinate system if it maps 1:1
    ! other than the location of the origin to the OEDGE coordinate system
    
    r_offset = r_shift
    z_offset = z_shift

    ! Assign chosen interpolation option
    ! 0 - no interpolation - cell value is returned
    ! 1 - nearest neighbours interpolation

    interpolate_opt = interpolate_option

    ! open input files 
    ! read array sizes from the geometry file

    !...  Open geometry file:

    call find_free_unit_number(igridfile)
    gridfilename = trim(filename)//'.grd'

    OPEN(UNIT=igridfile,FILE=trim(gridfilename),STATUS='OLD',IOSTAT=ios)

    if (ios.ne.0) then 

       call errmsg('ERROR OPENING OEDGE GRID FILE:'//trim(gridfilename),ios)

       stop 'ERROR OPENING OEDGE GRID FILE'

    else
       write(0,'(a,i10)') 'GRID FILE NAME:'//trim(gridfilename)//': opened for reading on unit number = ',igridfile 

    endif


    !...  Open plasma file:

    call find_free_unit_number(iplasmafile)
    plasmafilename = trim(filename)//'.bgp'

    OPEN(UNIT=iplasmafile,FILE=trim(plasmafilename),STATUS='OLD',IOSTAT=ios)

    if (ios.ne.0) then 

       call errmsg('ERROR OPENING OEDGE PLASMA FILE:'//trim(plasmafilename),ios)
       stop 'ERROR OPENING OEDGE PLASMA FILE'

    else

       write(0,'(a,i10)') 'PLASMA FILE NAME:'//trim(plasmafilename)//': opened for reading on unit number = ',iplasmafile 

    endif

    ! load geometry data and allocate storage for grid and plasma

    call load_oedge_grid(igridfile,ierr)

    ! load plasma data

    call load_oedge_plasma(iplasmafile,ierr)

    ! combine the target values into one array to make the interpolation routines easier
    ! interpolation is cell center based while the no-interpolation option returns the data in the polygon containing r,z

    call combine_target_data

    ! Initialize start cell to offgrid values
    ik_last = -1
    ir_last = -1
    iq_last = -1


  end subroutine init_oedge_plasma


  subroutine close_oedge_plasma
    implicit none

    ! deallocate any allocated storage

    call deallocate_storage


  end subroutine close_oedge_plasma


  subroutine get_oedge_plasma(rin,zin,ne,te,ti,vb,ef,btoto,bro,bzo,bto,ierr)
    implicit none

    real*8 :: rin,zin,ne,te,ti,vb,ef,btoto,bro,bzo,bto

    real*8 :: r,z

    integer :: ierr

    integer :: ik,ir,iq

    ! Adjust the input coordinates for any origin/coordinate shift or offset specified in the initialization routine
    ! This is useful to map coordinate systems that are otherwise 1:1 with a different origin
    
    r = rin + r_offset
    z = zin + z_offset 

    !write(0,'(a,2(1x,g18.8))') 'R,Z:',r,z,r_offset,z_offset
    !write(6,'(a,2(1x,g18.8))') 'R,Z:',r,z,r_offset,z_offset

    ! Get the OEDGE plasma conditions at the specified R,Z location
    ! Two options - value in cell and interpolated - value in cell is quicker - interpolated is smoother


    ! Get cell 
    call find_poly(r,z,ik,ir,ik_last,ir_last,ierr)


    ! If the R,Z location is not found on grid then exit with non-zero return code
    ! The plasma data returned is undefined (it will contain whatever was passed in)
    if (ierr.ne.0) then 
       return
    endif

    !
    ! If on grid assign last location variables
    !

    ik_last = ik
    ir_last = ir


    ! Interpolate result if required


    if (interpolate_opt.eq.0) then 


       ! no interpolation
       ne = knbs(ik,ir)
       te = ktebs(ik,ir)
       ti = ktibs(ik,ir)
       vb = kvhs(ik,ir)
       ef = kes(ik,ir)
       btoto = btot(ik,ir)
       bro = br(ik,ir)
       bzo = bz(ik,ir)
       bto = bt(ik,ir)

       if (debug_code) then 
          call write_cell_data(ik,ir)
       endif


    elseif (interpolate_opt.eq.1) then
       ! interpolate using Jeff's algorithm
       ! Note: ne,te,ti,vb,ef have target values that are used for interpolation in the first half cell
       !       btot,br,bz,bt do not have target data .. as a result the target values are the same as the first cell center
       ! Note: Target R,Z coordinates are needed for this algorithm.
       !
       ! Jeff's algorithm is based on cell centers and the values at the cell centers. It does not require
       ! the polygon information since it checks to see where the particle is relative to the values at the cell 
       ! centers and then interpolates from there using the perpendicular distance from the point to the lines joining the 
       ! cell centers as the weight factors. 
       ! This makes sense since no matter what the polygon shapes on the actual grid the plasma values are defined and 
       ! essentially calculated for the cell center (this may cause some issues for fluid code quantities like the 
       ! flow velocity that are cell edge based - but does not affect the OEDGE results since everything is cell center 
       ! calculated except at material surfaces)

       ! Find quadrant in cell for interpolation algorithm
       ! This should NOT throw an error - but check just in case
       call find_quadrant(r,z,ik,ir,iq,iq_last,ierr)


       ! If an error is found return uninterpolated results
       if (ierr.ne.0) then 
          ! no interpolation
          ne = knbs(ik,ir)
          te = ktebs(ik,ir)
          ti = ktibs(ik,ir)
          vb = kvhs(ik,ir)
          ef = kes(ik,ir)
          btoto = btot(ik,ir)
          bro = br(ik,ir)
          bzo = bz(ik,ir)
          bto = bt(ik,ir)

          if (debug_code) then 
             call write_cell_data(ik,ir)
          endif

       else

          iq_last = iq

          call interpolate_plasma(r,z,ik,ir,iq,ne,te,ti,vb,ef,btoto,bro,bzo,bto)

       endif

    endif

  end subroutine get_oedge_plasma


  subroutine interpolate_plasma(r,z,ik,ir,iq,ne,te,ti,vb,ef,btoto,bro,bzo,bto)

    real*8 :: r,z,ne,te,ti,vb,ef,btoto,bro,bzo,bto
    integer :: ik,ir,iq,in
    integer :: iks(4),irs(4)


    ! Interpolation is based on cell centers - this causes some issues near X-points for certain quadrants.


    ! Based on the value of iq - the four adjacent cell centers may be identified except at X-points - the problematic conditions can 
    ! be checked by looking for 
    ! 1) irins or irouts of the up or down cells depending on quadrant are not the same
    ! 2) boundary rings can be detected by irins = ir or irouts = ir - these need special treatment as well 

    logical :: xpt, boundary
    real*8 :: elow,ehigh,flow,fhigh,esum,fsum,e1,e2,f1,f2

    xpt = .false.
    boundary = .false.


    ! Get the R,Z coordinates of the 4 corners
    ! Always start at the UP,OUT corner and work around 
    ! Sides between vertices (1,2) and (3,4) are always going across the field lines
    ! Sides between (2,3) and (4,1) are always "parallel" to the field lines
    !
    ! NOTE: IRINS, IROUTS, IKINS, IKOUTS are not set for ik=0 or ik=nks(ir)+1 since these
    !       are the rows next to the target and do not represent actual cells on the grid. 
    !       Special care is taken below when these conditions are encountered. 

    ! IN and DOWN
    if (iq.eq.1) then 

       irs(1) = ir
       iks(1) = ik

       irs(2) = irins(ik,ir)
       iks(2) = ikins(ik,ir)

       if (ik.eq.1) then 
          irs(3) = irins(ik,ir)
          iks(3) = ikins(ik,ir)-1
       else
          irs(3) = irins(ik-1,ir)
          iks(3) = ikins(ik-1,ir)
       endif

       irs(4) = ir
       iks(4) = ik -1 

       ! Check error conditions
       ! Xpoint
       if (irs(2).ne.irs(3)) then 
          xpt = .true.
       endif

       ! boundary ring 
       if (irs(1).eq.irs(2)) then 
          boundary = .true.
       endif

       ! OUT and DOWN
    elseif (iq.eq.2) then

       irs(1) = irouts(ik,ir)
       iks(1) = ikouts(ik,ir)

       irs(2) = ir
       iks(2) = ik

       irs(3) = ir
       iks(3) = ik-1

       if (ik.eq.1) then 
          irs(4) = irouts(ik,ir)
          iks(4) = ikouts(ik,ir)-1
       else
          irs(4) = irouts(ik-1,ir)
          iks(4) = ikouts(ik-1,ir)
       endif

       ! Check error conditions
       ! Xpoint
       if (irs(1).ne.irs(4)) then 
          xpt = .true.

       endif

       ! boundary ring 
       if (irs(1).eq.irs(2)) then 
          boundary = .true.

       endif

       ! OUT and UP
    elseif (iq.eq.3) then

       if (ik.eq.nks(ir)) then 
          irs(1) = irouts(ik,ir)
          iks(1) = ikouts(ik,ir)+1
       else
          irs(1) = irouts(ik+1,ir)
          iks(1) = ikouts(ik+1,ir)
       endif

       irs(2) = ir
       iks(2) = ik +1 

       irs(3) = ir
       iks(3) = ik

       irs(4) = irouts(ik,ir)
       iks(4) = ikouts(ik,ir)

       ! Check error conditions
       ! Xpoint
       if (irs(1).ne.irs(4)) then 
          xpt = .true.
       endif

       ! boundary ring 
       if (irs(3).eq.irs(4)) then 
          boundary = .true.
       endif

       ! IN and UP
    elseif (iq.eq.4) then 

       irs(1) = ir
       iks(1) = ik+1

       if (ik.eq.nks(ir)) then 
          irs(2) = irins(ik,ir)
          iks(2) = ikins(ik,ir)+1
       else
          irs(2) = irins(ik+1,ir)
          iks(2) = ikins(ik+1,ir)
       endif

       irs(3) = irins(ik,ir)
       iks(3) = ikins(ik,ir)

       irs(4) = ir
       iks(4) = ik

       ! Check error conditions
       ! Xpoint
       if (irs(2).ne.irs(3)) then 
          xpt = .true.
       endif

       ! boundary ring 
       if (irs(3).eq.irs(4)) then 
          boundary = .true.
       endif

    endif


    ! If in quadrant adjacent to Xpoint or boundary no interpolation takes place (can be upgraded later) - just return the values in the cell. 

    if (xpt.or.boundary) then 

       ! no interpolation
       ne = knbs(ik,ir)
       te = ktebs(ik,ir)
       ti = ktibs(ik,ir)
       vb = kvhs(ik,ir)
       ef = kes(ik,ir)
       btoto = btot(ik,ir)
       bro = br(ik,ir)
       bzo = bz(ik,ir)
       bto = bt(ik,ir)

       write(6,'(a,2l10,2i10,20(1x,g18.8))') 'X or B:',xpt,boundary,ik,ir,r,z,ne,te,ti
       write(6,'(a,20i8,20(1x,g18.8))') 'X or B:',ik,ir,nks(ir),ikins(ik,ir),irins(ik,ir),ikouts(ik,ir),irouts(ik,ir),((iks(in),irs(in)),in=1,4)


    else


       !    interpolation
       !    f = along field line - constant i moves along a r surface
       !    e = across the field line - constant j moves along a z surface


       ! Distance to side (1,2)
       call distance( r,z, rs(iks(1),irs(1)),zs(iks(1),irs(1)),rs(iks(2),irs(2)),zs(iks(2),irs(2)), flow )
       !call distance( r,z, r(i,j),zz(i,j),r(i,j+1),zz(i,j+1), flow )

       ! Distance to side (3,4)
       call distance( r,z, rs(iks(3),irs(3)),zs(iks(3),irs(3)),rs(iks(4),irs(4)),zs(iks(4),irs(4)), fhigh )
       !call distance( ru,zu, r(i+1,j),zz(i+1,j),r(i+1,j+1),zz(i+1,j+1), fhigh )

       ! Distance to side (2,3)
       call distance( r,z, rs(iks(2),irs(2)),zs(iks(2),irs(2)),rs(iks(3),irs(3)),zs(iks(3),irs(3)), elow )
       !call distance( ru,zu, r(i,j),zz(i,j),r(i+1,j),zz(i+1,j), elow )

       ! Distance to side (4,1)
       call distance( r,z, rs(iks(4),irs(4)),zs(iks(4),irs(4)),rs(iks(1),irs(1)),zs(iks(1),irs(1)), ehigh )
       !call distance( ru,zu, r(i,j+1),zz(i,j+1),r(i+1,j+1),zz(i+1,j+1), ehigh )

       esum = elow+ehigh
       fsum = flow+fhigh

       e1=ehigh/esum
       e2=elow/esum
       f1=fhigh/fsum
       f2=flow/fsum



       call assign_interpolate(irs,iks,e1,e2,f1,f2,knbs,ne)

       if (debug_code) then 
          call debug_interpolate('KNBS',irs,iks,e1,e2,f1,f2,knbs,ne)
       endif

       call assign_interpolate(irs,iks,e1,e2,f1,f2,ktebs,te)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,ktibs,ti)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,kvhs,vb)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,kes,ef)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,btot,btoto)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,br,bro)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,bz,bzo)
       call assign_interpolate(irs,iks,e1,e2,f1,f2,bt,bto)

    endif


  end subroutine interpolate_plasma


  subroutine assign_interpolate(irs,iks,e1,e2,f1,f2,array,val)
    implicit none
    integer :: irs(4),iks(4)
    real*8 :: e1,e2,f1,f2
    real*8,allocatable :: array(:,:)
    real*8 :: val,v12,v34

    ! Array should always be allocated when this routine is called
    if (.not.allocated(array)) then 
       call errmsg('ASSIGN_INTERPOLATE:','CRITICAL ERROR: Array not allocated')
       stop 'Assign_interpolate'
    endif

    v12 = array(iks(1),irs(1)) * e2 + array(iks(2),irs(2)) * e1
    v34 = array(iks(4),irs(4)) * e2 + array(iks(3),irs(3)) * e1

    val = v12 * f1 + v34 * f2



  end subroutine assign_interpolate


  subroutine debug_interpolate(desc,irs,iks,e1,e2,f1,f2,array,val)
    implicit none
    character*(*) :: desc
    integer :: irs(4),iks(4)
    real*8 :: e1,e2,f1,f2
    real*8,allocatable :: array(:,:)
    real*8 :: val,v12,v34

    ! Array should always be allocated when this routine is called
    if (.not.allocated(array)) then 
       call errmsg('DEBUG_INTERPOLATE:','CRITICAL ERROR: Array not allocated')
       stop 'Assign_interpolate'
    endif

    v12 = array(iks(1),irs(1)) * e2 + array(iks(2),irs(2)) * e1
    v34 = array(iks(4),irs(4)) * e2 + array(iks(3),irs(3)) * e1

    !val = v12 * f1 + v34 * f2

    call write_cell_data(iks(1),irs(1))
    call write_cell_data(iks(2),irs(2))
    call write_cell_data(iks(3),irs(3))
    call write_cell_data(iks(4),irs(4))

    write(0,'(a,12(1x,g14.6))') 'DEBUG:'//trim(desc),array(iks(1),irs(1)),array(iks(2),irs(2)),array(iks(4),irs(4)),array(iks(3),irs(3))
    write(0,'(a,12(1x,g14.6))') 'DEBUG:'//trim(desc),array(iks(1),irs(1))*e2,array(iks(2),irs(2)) * e1,array(iks(4),irs(4))*e2,array(iks(3),irs(3)) * e1
    write(0,'(a,12(1x,g14.6))') 'DEBUG:'//trim(desc),e1,e2,f1,f2,v12,v34, v12 * f1 + v34 * f2,val

  end subroutine debug_interpolate




  subroutine find_poly(r,z,ik,ir,ik_last,ir_last,ierr)
    implicit none
    real*8 :: r,z
    integer :: ik,ir,ik_last,ir_last,ierr
    logical :: found,done

    ierr = 0
    found = .false.

    if (ir_last.ge.1.and.ir_last.le.nrs) then 
       if (ik_last.ge.0.and.ik_last.le.nks(ir)+1) then 
          ! Start search in previous cell then surrounding cells
          ik = ik_last
          ir = ir_last
          found = incell(ik,ir,r,z)

          ! previous cell along the field line
          if (.not.found.and.ik_last.ne.1) then 
             ik = ik_last -1
             found = incell(ik,ir,r,z)
          endif

          ! next cell along the field line
          if (.not.found.and.ik_last.ne.nks(ir)) then 
             ik = ik_last +1
             found = incell(ik,ir,r,z)
          endif

          ! inward cell
          if (.not.found) then 
             ik = ikins(ik_last,ir_last)
             ir = irins(ik_last,ir_last)
             found = incell(ik,ir,r,z)
          endif

          ! outward cell
          if (.not.found) then 
             ik = ikouts(ik_last,ir_last)
             ir = irouts(ik_last,ir_last)
             found = incell(ik,ir,r,z)
          endif

       endif
    endif


    ! Start full search

    ! Search outer divertor, SOL, PFZ
    if (.not.found) call search_cells(r,z,ik,ir,irsep,nrs,1,found)
    ! Search inner divertor, SOL, PFZ
    if (.not.found) call search_cells(r,z,ik,ir,irsep,nrs,2,found)
    ! Search confined plasma 
    if (.not.found) call search_cells(r,z,ik,ir,1,irsep-1,0,found)



    if (.not.found) then 
       ! Set error condition - position not found on grid

       call errmsg('FIND_POLY:','Position not found on grid')
       ierr = 1

    endif


  end subroutine find_poly


  subroutine find_quadrant(r,z,ik,ir,iq,iq_last,ierr) 
    implicit none
    real*8 :: r,z
    integer :: ik,ir,iq,iq_last,ierr

    ! local variables
    real*8 :: rvert(4),zvert(4),rside(4),zside(4)
    real*8 :: rv(4),zv(4)
    integer :: nv,iv,ivf,ivlast,ivnext
    logical :: found

    !logical,external :: inpoly

    integer :: cnt,in

    ierr = 0
    iq = 0

    in = korpg(ik,ir)
    nv = nvertp(in)

    if (nv.ne.4) then 
       call errmsg('FIND_QUADRANT:Quadrant code only works for 4 sided polygons',nv)
       ierr = 1
       return
    endif

    rvert = rvertp(:,in)
    zvert = zvertp(:,in)


    !
    !     Split the polygon into 4 pieces and determine which part of the
    !     cell the point R,Z lies in - adjust 
    !     Calculate the midpoints of the sides
    !

    do iv = 1,nv

       ivnext = iv+1
       if (iv.eq.nv) ivnext = 1 

       rside(iv) = (rvert(iv)+rvert(ivnext))/2.0
       zside(iv) = (zvert(iv)+zvert(ivnext))/2.0

    end do

    !rcp = rs(ik,ir)
    !zcp = zs(ik,ir)
    !rcp = (rside(1) + rside(3))/2.0
    !zcp = (zside(1) + zside(3))/2.0

    ! start with quadrant iq_last or 1 if iq_last is not defined

    ivf= 0
    if (iq_last.ge.1.and.iq_last.le.4) then 
       iv = iq_last
    else
       iv = 1
    endif

    found = .false. 
    cnt = 1

    !do while(iv.le.4.and.(.not.found)) 
    do while(cnt.le.4.and.(.not.found)) 
       ivlast = iv-1
       if (iv.eq.1) ivlast = 4 
       !
       !        Determine corners of polygon to check 
       !
       rv(mod(iv-1,4)+1)   = rvert(iv) 
       rv(mod(1+iv-1,4)+1) = rside(iv)

       rv(mod(2+iv-1,4)+1) = rs(ik,ir)
       !rv(mod(2+iv-1,4)+1) = rcp
       rv(mod(3+iv-1,4)+1) = rside(ivlast)

       zv(mod(iv-1,4)+1)   = zvert(iv) 
       zv(mod(1+iv-1,4)+1) = zside(iv)

       zv(mod(2+iv-1,4)+1) = zs(ik,ir)
       !zv(mod(2+iv-1,4)+1) = zcp
       zv(mod(3+iv-1,4)+1) = zside(ivlast)

       found = inpoly(r,z,nv,rv,zv) 

       if (found) ivf = iv

       iv = iv+1
       if (iv.gt.4) iv = 1
       cnt = cnt + 1

    end do

    !
    !     Quadrant is found
    !
    if (ivf.ne.0) then 
       iq = ivf
    else

       !        Error    

       ierr = 1


       call errmsg( 'ERROR in find_quadrant:','Point not found in cell')
       !write(6,'(a,8(1x,g12.5))') 'Last poly:',((rv(iv),zv(iv)),iv=1,4)
       !write(6,'(a,8(1x,g12.5))') 'Point R,Z,S,C:',r,z,s_frac,cross_frac

    endif

  end subroutine find_quadrant




  subroutine search_cells(r,z,ik,ir,ir_start,ir_end,opt,found)
    implicit none
    real*8 :: r,z
    integer :: ik,ir
    integer :: ir_start,ir_end
    integer :: opt
    logical :: found


    integer :: ir_tmp,ik_tmp
    integer :: ik_start,ik_end,ik_step

    found = .false.

    do ir_tmp = ir_start,ir_end

       if (opt.eq.0) then 
          ik_start = 1
          ik_end = nks(ir_tmp)
          ik_step = 1
       elseif (opt.eq.1) then
          ik_start = nks(ir_tmp)
          ik_end = nks(ir_tmp)/2
          ik_step = -1
       elseif (opt.eq.2) then
          ik_start = 1
          ik_end = nks(ir_tmp)/2 -1
          ik_step = 1
       endif

       do ik_tmp = ik_start,ik_end,ik_step

          found = incell(ik_tmp,ir_tmp,r,z)

          if (found) then 
             ik = ik_tmp
             ir = ir_tmp
             return
          endif


       end do

    end do

  end subroutine search_cells

  subroutine load_oedge_grid(infile,ierr)
    implicit none
    integer :: infile,ierr

    integer :: ik,ir,in
    logical :: done

    character*512 :: buffer

    integer :: ios

    !
    !     READ Title line
    !
    read (infile,10) buffer
    !
    if (buffer(1:6).ne.'DIVIMP') then
       call errmsg('NOT A DIVIMP GRID FILE',trim(buffer))
       stop 'Not a valid DIVIMP grid file'
    endif

    ! First line contains overal grid size information


    done = .false.

    ! initialization
    nrs   = 0
    irsep = 0
    nds   = 0
    npolyp= 0
    nvert = 0


    do while (.not.done)

       read(infile,10,iostat=ios,err=1000) buffer

       if (ios.eq.0) then 

          if (buffer(1:4).eq.'NRS:') then
             read (buffer,200) nrs,irsep,nds,npolyp,nvert

             ! Allocate the nks array
             if (allocated(nks)) deallocate(nks)
             if (nrs.gt.0) then 
                allocate(nks(nrs),stat=ierr)
                if (ierr.ne.0) then 
                   call errmsg('ERROR Allocating NKS array:',ierr)
                   stop 'Error allocating NKS array'
                endif
             else
                call errmsg('INVALID GRID FILE: NUMBER OF RINGS: NRS = ',nrs)
                stop 'INVALID GRID FILE'
             endif

          elseif(buffer(1:6).eq.'KNOTS:') then

             if (.not.allocated(nks)) then 
                call errmsg('ERROR NKS Array not yet allocated')
                stop 'NKS not allocated'
             endif

             read (infile,400)  (nks(ir),ir=1,nrs)

             mks = maxval(nks)+1

             ! Allocate the storage for all of the arrays 
             call allocate_storage

          elseif (buffer(1:3).eq.'RS:') then
             !
             !       R cell center coordinate
             !
             read (infile,500) ((rs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:3).eq.'ZS:') then
             !
             !       Z cell center coordinate
             !
             read (infile,500) ((zs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:5).eq.'RBND:') then
             !
             !       R cell bound coordinate
             !
             read (infile,500) ((rbnd(ik,ir),ik=0,nks(ir)),ir=1,nrs)

          elseif (buffer(1:5).eq.'ZBND:') then
             !
             !       Z cell bound coordinate
             !
             read (infile,500) ((zbnd(ik,ir),ik=0,nks(ir)),ir=1,nrs)

          elseif (buffer(1:6).eq.'IRINS:') then
             !
             !       IR Inward connection map
             !
             read (infile,400) ((irins(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:6).eq.'IKINS:') then
             !
             !       IK Inward connection map
             !
             read (infile,400) ((ikins(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:7).eq.'IROUTS:') then
             !
             !       IR outward connection map
             !
             read (infile,400) ((irouts(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:7).eq.'IKOUTS:') then
             !
             !       IK outward connection map
             !
             read (infile,400) ((ikouts(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:5).eq.'IDDS:') then
             !
             !       Target element connection map
             !
             read (infile,400) ((idds(ir,in),ir=1,nrs),in=1,2)

          elseif (buffer(1:6).eq.'KORPG:') then
             !
             !       Cell polygon data index
             !
             read (infile,400) ((korpg(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:7).eq.'NVERTP:') then
             !
             !       Number of vertices in each polygon
             !
             read (infile,400) (nvertp(in),in=1,npolyp)

          elseif (buffer(1:7).eq.'RVERTP:') then
             !
             !    R coordinates of polygon vertices
             !
             read (infile,500) ((rvertp(ik,in),ik=1,nvert),in=1,npolyp)

          elseif (buffer(1:7).eq.'ZVERTP:') then
             !
             !    Z coordinates of polygon vertices
             !
             read (infile,500) ((zvertp(ik,in),ik=1,nvert),in=1,npolyp)

          elseif (buffer(1:5).eq.'BTOT:') then
             !
             !       B-field magnitude in cell
             !
             read (infile,500) ((btot(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:3).eq.'BR:') then
             !
             !       R part of B-field unit vector
             !
             read (infile,500) ((br(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:3).eq.'BZ:') then
             !
             !       Z part of B-field unit vector
             !
             read (infile,500) ((bz(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:3).eq.'BT:') then
             !
             !       T part of B-field unit vector
             !
             read (infile,500) ((bt(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          endif

       else
          done=.true.
       endif
       !
       !     Loop back for continued reading
       !
    end do


    CLOSE (infile)


    ! Assign target data for magnetic field until this is put in transfer file (if needed)
    do ir = irsep,nrs
       btotd(idds(ir,1)) = btot(nks(ir),ir)
       btotd(idds(ir,2)) = btot(1,ir)

       brd(idds(ir,1)) = br(nks(ir),ir)
       brd(idds(ir,2)) = br(1,ir)

       bzd(idds(ir,1)) = bz(nks(ir),ir)
       bzd(idds(ir,2)) = bz(1,ir)

       btd(idds(ir,1)) = bt(nks(ir),ir)
       btd(idds(ir,2)) = bt(1,ir)
    end do


    return

1000 continue

    call errmsg('ERROR READING IN DIVIMP GRID:')
    stop 'ERROR reading in grid file'


    !
    !     Formatting
    !
10  format(a)
100 format(a40)
200 format('NRS:',i5,'IRSEP:',i5,'NDS:',i5,'NPOLYP:',i5,'NVERT:',i5)
400 format(12i6)
500 format(6e18.10)


  end subroutine load_oedge_grid

  subroutine allocate_storage
    use allocate_arrays
    implicit none
    integer :: ierr


    ! This routine uses the values set in nrs, mks, nds, npolyp and nvert to allocate the storage for all of the 
    ! grid and plasma array data.

    ! Allocate basic grid geometry arrays
    call allocate_array(rs,0,mks,1,nrs,'RS',ierr)
    call allocate_array(zs,0,mks,1,nrs,'ZS',ierr)
    call allocate_array(rbnd,0,mks,1,nrs,'RBND',ierr)
    call allocate_array(zbnd,0,mks,1,nrs,'ZBND',ierr)

    ! Allocate connection map arrays
    call allocate_array(irins,0,mks,1,nrs,'IRINS',ierr)
    call allocate_array(ikins,0,mks,1,nrs,'IKINS',ierr)
    call allocate_array(irouts,0,mks,1,nrs,'IROUTS',ierr)
    call allocate_array(ikouts,0,mks,1,nrs,'IKOUTS',ierr)
    call allocate_array(idds,nrs,2,'IDDS',ierr)

    ! Allocate polygon map and cell information
    call allocate_array(korpg,0,mks,1,nrs,'KORPG',ierr)
    call allocate_array(nvertp,npolyp,'NVERTP',ierr)
    call allocate_array(rvertp,nvert,npolyp,'RVERTP',ierr)
    call allocate_array(zvertp,nvert,npolyp,'ZVERTP',ierr)

    ! Allocate magnetic field arrays
    call allocate_array(btot,0,mks,1,nrs,'BTOT',ierr)
    call allocate_array(br,0,mks,1,nrs,'BR',ierr)
    call allocate_array(bz,0,mks,1,nrs,'BZ',ierr)
    call allocate_array(bt,0,mks,1,nrs,'BT',ierr)

    ! Allocate plasma arrays
    ! Volume plasma background
    call allocate_array(knbs,0,mks,1,nrs,'KNBS',ierr)
    call allocate_array(ktebs,0,mks,1,nrs,'KTEBS',ierr)
    call allocate_array(ktibs,0,mks,1,nrs,'KTIBS',ierr)
    call allocate_array(kvhs,0,mks,1,nrs,'KVHS',ierr)
    call allocate_array(kes,0,mks,1,nrs,'KES',ierr)

    ! Target plasma background
    call allocate_array(knds,nds,'KNDS',ierr)
    call allocate_array(kteds,nds,'KTEDS',ierr)
    call allocate_array(ktids,nds,'KTIDS',ierr)
    call allocate_array(kvds,nds,'KVDS',ierr)
    call allocate_array(keds,nds,'KEDS',ierr)

    ! Target magnetic field 
    ! Set to values at first cell center
    call allocate_array(btotd,nds,'BTOTD',ierr)
    call allocate_array(brd,nds,'BTOTD',ierr)
    call allocate_array(bzd,nds,'BTOTD',ierr)
    call allocate_array(btd,nds,'BTOTD',ierr)


  end subroutine allocate_storage

  subroutine deallocate_storage
    implicit none


    ! Deallocate storage used in the OEDGE plasma interface

    if (allocated(nks)) deallocate(nks)

    if (allocated(rs)) deallocate(rs)
    if (allocated(zs)) deallocate(zs)
    if (allocated(rbnd)) deallocate(rbnd)
    if (allocated(zbnd)) deallocate(zbnd)

    if (allocated(irins)) deallocate(irins)
    if (allocated(ikins)) deallocate(ikins)
    if (allocated(irouts)) deallocate(irouts)
    if (allocated(ikouts)) deallocate(ikouts)
    if (allocated(idds)) deallocate(idds)

    if (allocated(korpg)) deallocate(korpg)
    if (allocated(nvertp)) deallocate(nvertp)
    if (allocated(rvertp)) deallocate(rvertp)
    if (allocated(zvertp)) deallocate(zvertp)

    if (allocated(btot)) deallocate(btot)
    if (allocated(br)) deallocate(br)
    if (allocated(bz)) deallocate(bz)
    if (allocated(bt)) deallocate(bt)

    if (allocated(knbs)) deallocate(knbs)
    if (allocated(ktebs)) deallocate(ktebs)
    if (allocated(ktibs)) deallocate(ktibs)
    if (allocated(kvhs)) deallocate(kvhs)
    if (allocated(kes)) deallocate(kes)

    if (allocated(knds)) deallocate(knds)
    if (allocated(kteds)) deallocate(kteds)
    if (allocated(ktids)) deallocate(ktids)
    if (allocated(kvds)) deallocate(kvds)
    if (allocated(keds)) deallocate(keds)

  end subroutine deallocate_storage




  subroutine load_oedge_plasma(infile,ierr)
    implicit none
    integer :: infile,ierr

    !
    !     LOAD_OEDGE_PLASMA:The purpose of this routine is to read in the
    !               DIVIMP background plasma in a DIVIMP specific
    !               format.
    !


    character*256 :: buffer
    integer :: ik,ir,id,ios
    integer :: tmpnrs,tmpnds,tmpirsep
    integer,allocatable :: tmpnks(:)
    logical :: done
    !
    !     READ Title line
    !


    read (infile,10) buffer
    !
    if (buffer(1:6).ne.'DIVIMP') then
       call errmsg('NOT A DIVIMP PLASMA FILE',trim(buffer))
       stop 'Not a valid DIVIMP plasma file'
    endif

    done = .false.

    do while (.not.done) 

       read(infile,10,iostat=ios,err=1000) buffer

       if (ios.eq.0) then 

          if (buffer(1:4).eq.'NRS:') then
             read (buffer,200) tmpnrs,tmpirsep,tmpnds
             !
             !        Check to see if this matches the current grid.
             !
             if (nrs.ne.tmpnrs.or.irsep.ne.tmpirsep.or.nds.ne.tmpnds) then
                !
                !           Grid characteristic mismatch - exit program.
                !
                call errmsg('DIVIMP PLASMA FILE DOES NOT MATCH GRID')
                write (0,*) 'NRS  :',nrs,tmpnrs
                write (0,*) 'IRSEP:',irsep,tmpirsep
                write (0,*) 'NDS  :',nds,tmpnds
                write (0,*) 'PROGRAM EXITING'
                stop 'DIVIMP PLASMA/GRID FILE MISMATCH'

             end if

             ! Allocate tmpnks

             if (allocated(tmpnks)) deallocate(tmpnks)
             allocate(tmpnks(tmpnrs),stat=ierr)
             if (ierr.ne.0) then 
                call errmsg('LOAD_OEDGE_PLASMA: Problem allocating tmpnks',ierr)
                stop 'LOAD_OEDGE_PLASMA:tmpnks'
             endif

          elseif(buffer(1:6).eq.'KNOTS:') then

             read (infile,400)  (tmpnks(ir),ir=1,nrs)
             !
             !        Check to see if knots match
             !
             do ir = 1, nrs

                if (nks(ir).ne.tmpnks(ir)) then

                   call errmsg('DIVIMP PLASMA FILE DOES NOT MATCH GRID')
                   write (0,*) 'IR     :',ir
                   write (0,*) 'NKS(IR):',nks(ir),tmpnks(ir)
                   write (0,*) 'PROGRAM EXITING'
                   stop 'DIVIMP PLASMA/GRID FILE MISMATCH'
                end if
             end do

          elseif (buffer(1:5).eq.'KNBS:') then
             !
             !        Density - volume
             !
             read (infile,500) ((knbs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:5).eq.'KNDS:') then
             !
             !        Density - target
             !
             read (infile,500) (knds(id),id=1,nds)

          elseif (buffer(1:6).eq.'KTEBS:') then
             !
             !        Te - volume
             !
             read (infile,500) ((ktebs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:6).eq.'KTEDS:') then
             !
             !        Te - target
             !
             read (infile,500) (kteds(id),id=1,nds)

          elseif (buffer(1:6).eq.'KTIBS:') then
             !
             !        Ti - volume
             !
             read (infile,500) ((ktibs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:6).eq.'KTIDS:') then
             !
             !        Ti - target
             !
             read (infile,500) (ktids(id),id=1,nds)

          elseif (buffer(1:5).eq.'KVHS:') then
             !
             !        Velocity - volume
             !
             read (infile,500) ((kvhs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:5).eq.'KVDS:') then
             !
             !        Velocity - target
             !
             read (infile,500) (kvds(id),id=1,nds)

          elseif (buffer(1:4).eq.'KES:') then
             !
             !        Electric Field - volume
             !
             read (infile,500) ((kes(ik,ir),ik=1,nks(ir)),ir=1,nrs)

          elseif (buffer(1:5).eq.'KEDS:') then
             !
             !        Electric Field - target
             !
             read (infile,500) (keds(id),id=1,nds)

          endif

       else
          done=.true.
       endif
       !
       !     Loop back for continued reading
       !
    end do

    if (allocated(tmpnks)) deallocate(tmpnks)

    CLOSE (infile)



    return

1000 continue

    call errmsg('ERROR READING IN DIVIMP PLASMA BACKGROUND')
    stop 'ERROR reading in Plasma file'


    !
    !     Formatting
    !
10  format(a)
100 format(a40)
200 format('NRS:',i5,'IRSEP:',i5,'NDS:',i5)
400 format(12i6)
500 format(6e18.10)

  end subroutine load_oedge_plasma


  subroutine combine_target_data
    implicit none
    integer :: in,ir,ik,id

    ! This routine loads the target values into the 0 and nks(ir)+1 elements of each data block

    ! Load the target coordinates first

    do ir = 1,nrs
       ! Note for in < irsep this has little meaning but it won't affect the algorithm
       ! R coordinate
       rs(0,ir) = rbnd(0,ir)
       rs(nks(ir)+1,ir) = rbnd(nks(ir),ir)
       ! Z coordinate
       zs(0,ir) = zbnd(0,ir)
       zs(nks(ir)+1,ir) = zbnd(nks(ir),ir)
    end do

    ! Load magnetic field data - values at the target are not available - first cell is extended to target

    do ir = 1,nrs
       btot(0,ir) = btot(1,ir)
       br(0,ir)   = br(1,ir)
       bz(0,ir)   = bz(1,ir)
       bt(0,ir)   = bt(1,ir)

       btot(nks(ir)+1,ir) = btot(nks(ir),ir)
       br(nks(ir)+1,ir)   = br(nks(ir),ir)
       bz(nks(ir)+1,ir)   = bz(nks(ir),ir)
       bt(nks(ir)+1,ir)   = bt(nks(ir),ir)
    end do

    ! Load target data for plasma quantities

    do ir = irsep,nrs
       id = idds(ir,1)
       knbs(nks(ir)+1,ir)  = knds(id)
       ktebs(nks(ir)+1,ir) = kteds(id)
       ktibs(nks(ir)+1,ir) = ktids(id)
       kvhs(nks(ir)+1,ir)  = kvds(id)
       kes(nks(ir)+1,ir)   = keds(id)

       id = idds(ir,2)
       knbs(0,ir)  = knds(id)
       ktebs(0,ir) = kteds(id)
       ktibs(0,ir) = ktids(id)
       kvhs(0,ir)  = kvds(id)
       kes(0,ir)   = keds(id)
    end do
    
    ! Fill in ends of core plasma rings

    do ir = 1,irsep-1
       knbs(nks(ir)+1,ir)  = knbs(nks(ir),ir) 
       ktebs(nks(ir)+1,ir) = ktebs(nks(ir),ir)
       ktibs(nks(ir)+1,ir) = ktibs(nks(ir),ir)
       kvhs(nks(ir)+1,ir)  = kvhs(nks(ir),ir) 
       kes(nks(ir)+1,ir)   = kes(nks(ir),ir)  

       knbs(0,ir)  = knbs(1,ir) 
       ktebs(0,ir) = ktebs(1,ir)
       ktibs(0,ir) = ktibs(1,ir)
       kvhs(0,ir)  = kvhs(1,ir) 
       kes(0,ir)   = kes(1,ir)  
    end do



    ! If needed storage can be deallocated here for some quantities


  end subroutine combine_target_data




  subroutine find_free_unit_number(unit)
    implicit none
    integer unit
    !
    !     FIND_FREE_UNIT_NUMBER:
    !
    !     This routine scans through unit numbers looking for one that
    !     is not currently in use. This number is returned. This code
    !     is based on the assumption that any unit numbers returned will
    !     be used before this routine is called again asking for another 
    !     number - otherwise it will likely return the previous value.
    !
    integer test_unit
    logical unit_open

    test_unit = 10
    unit_open = .true.

    ! Check for unit number assignment.  
    Do While (Unit_open)
       test_unit=test_unit + 1
       Inquire (Unit = test_unit, Opened = Unit_open)
    End Do

    unit = test_unit

    return
  end subroutine find_free_unit_number




  logical function incell(ik,ir,r,z)
    implicit none
    integer ik,ir
    real*8 r,z
    !      include 'params'
    !      include 'cgeom'
    !
    !     INCELL: This function returns a simple YES/NO decision
    !             about whether the point R,Z is in the cell designated
    !             by IK,IR with a set of vertices defined in an ordered
    !             clockwise fashion. It takes the cross product
    !             between the vector from the vertex to the test point
    !             and the vector from the vertex to the next clockwise
    !             vertex of the polygon. The cross-product must be
    !             the same sign for all vertices - if the
    !             point is outside the polygon it will fail this test
    !             for at least one vertex. (i.e. the cross-product will
    !             be less than zero.) (Suggested solution courtesy
    !             of Ian Youle :-) )
    !
    !             David Elder, Dec 8, 1993
    !
    !             Note: the objectives of the solution method were
    !             simplicity and reasonable computational cost.
    !             This solution avoids the need for square roots
    !             or trigonometric calculations.
    !
    !             Note: in the confined plasma the first and last cells
    !                   are identical. However, S=0 and S=SMAX are at the 
    !                   center of this cell. This causes some inconsistencies
    !                   when calculating particle positions. In order
    !                   to address this - this routine will consder a paricle
    !                   in the first cell on a core ring when it is in the
    !                   second half of the cell and in the last cell of a core
    !                   ring when it is in the first half of the cell. 
    !
    integer :: k,v,nextv,i,nv
    real*8 :: vxr,vxz,vwr,vwz,cp,lastcp
    !
    logical :: res
    !logical,external :: inpoly

    real*8 rc(4),zc(4)

    lastcp = 0.0

    incell = .false.
    k = korpg(ik,ir)
    if (k.eq.0) return
    nv = nvertp(k)
    if (nv.eq.0) return
    do v = 1,nv
       if (v.eq.nv) then
          nextv = 1
       else
          nextv = v+1
       endif
       !
       !        Want the vector cross-product Rx X Rw
       !
       !         vxr = r - rvertp(v,k)
       !         vxz = z - zvertp(v,k)
       !         vwr = rvertp(nextv,k) - rvertp(v,k)
       !         vwz = zvertp(nextv,k) - zvertp(v,k)
       !
       !         cp = vxr*vwz - vxz*vwr
       !
       !         if (cp.lt.0.0)  return
       !
       !          if (   (
       !     >     ( (r-rvertp(v,k)) *
       !     >       (zvertp(nextv,k)-zvertp(v,k)) )
       !     >    -( (z-zvertp(v,k)) *
       !     >       (rvertp(nextv,k)-rvertp(v,k)) )
       !     >           )
       !     >         .lt.0.0) return
       !
       cp =    (((r-rvertp(v,k)) * (zvertp(nextv,k)-zvertp(v,k))) &
            -((z-zvertp(v,k)) * (rvertp(nextv,k)-rvertp(v,k))))

       !
       !         There is a problem for points that should 
       !         lie on the boundary of the cell - i.e. that 
       !         are calculated based on the polygon corners and 
       !         which are mathematically on the polygon surface. 
       !         Numerically, these points can have a cross product
       !         which is close to zero but can vary to either side. 
       !         In order to consider these points in the cell - the 
       !         cross products are set to zero for values less than
       !         a specified limit. In this case the limit is set to 1.0e-7 
       !
       !         This value was determined by examining the range of cross 
       !         product values generated when sampling 50,000 points 
       !         calculated on a polygon with a scale size of 1.0m. 
       !         The maximum error cross product in this case was 6e-8.
       !
       !         D. Elder, Dec 13, 2006
       !
       !         Upon consideration - it might be best to not allow these points
       !         to be considered inside the cell since if they are detected they
       !         can be moved slightly to an appropriate location. 
       !
       !          if (abs(cp).lt.1.0e-7) cp = 0.0 
       !
       if (v.eq.1.and.cp.ne.0.0) lastcp = cp

       if ((lastcp * cp).lt.0.0) return

       if (cp.ne.0.0) lastcp = cp

    end do

    !
    !     Particle has been found in cell
    !
    incell = .true.
    !
    !     Check particles in the first or last cell of core rings
    !     for more accurate assessement.
    !
    if (ir.lt.irsep.and.(ik.eq.1.or.ik.eq.nks(ir))) then 
       !
       !        For a particle found to be in the first or last cell of 
       !       a core ring - need to decide if it is in the first
       !        half or second half and revise incell result accordingly. 
       !
       !        Only the second half of the first cell or the first half
       !        of the last cell should return true
       !
       !
       !        NOTE: Keep in mind that the first and last cells of core
       !              rings are supposed to be identical - if this changes
       !              then this code needs to be modified. 
       !         
       !        Check to see if particle is in the first half of the 
       !        cell - set vertices for call to inpoly.-  
       !         
       !        k was set to cell geometry index at the beginning of this
       !        routine.
       !
       rc(1) = rvertp(1,k)
       rc(2) = rvertp(2,k)
       rc(3) = (rvertp(2,k) + rvertp(3,k)) /2.0
       rc(4) = (rvertp(1,k) + rvertp(4,k)) /2.0

       zc(1) = zvertp(1,k)
       zc(2) = zvertp(2,k)
       zc(3) = (zvertp(2,k) + zvertp(3,k)) /2.0
       zc(4) = (zvertp(1,k) + zvertp(4,k)) /2.0
       !
       !        Check to see of point is in first half of cell
       !
       res = inpoly(r,z,4,rc,zc)
       !
       !        Change value of incell to false for either of the 
       !        invalid cases
       !        First half and in first cell
       !        Last half and in last cell. 
       !

       if ((ik.eq.1.and.res).or.(ik.eq.nks(ir).and.(.not.res))) then 
          incell = .false.
       endif

    endif


    return
  end function incell


  logical function inpoly(r,z,nv,rvert,zvert)
    implicit none
    integer ik,ir
    integer nv
    !      parameter (maxvert=8)
    real*8 r,z,rvert(nv),zvert(nv)
    !
    !     INPOLY: This function returns a simple YES/NO decision
    !             about whether the point R,Z is in the cell designated
    !             by a set of vertices defined in an ordered fashion.
    !             It takes the cross product
    !             between the vector from the vertex to the test point
    !             and the vector from the vertex to the next
    !             vertex of the polygon. The cross-product must be
    !             the same sign for all vertices - if the
    !             point is outside the polygon it will fail this test
    !             for at least one vertex. (i.e. the cross-product will
    !             change sign) (Suggested solution courtesy
    !             of Ian Youle :-) )
    !
    !             David Elder, Dec 8, 1993
    !
    !             Note: the objectives of the solution method were
    !             simplicity and reasonable computational cost.
    !             This solution avoids the need for square roots
    !             or trigonometric calculations.
    !
    integer v,nextv
    real*8 cp,lastcp

    lastcp = 0.0 

    inpoly = .false.

    if (nv.eq.0) return  
    !
    !     Loop through vertices
    !
    do v = 1,nv
       !
       if (v.eq.nv) then
          nextv = 1
       else
          nextv = v+1
       endif
       !
       !        Want the vector cross-product Rx X Rw
       !
       !         vxr = r - rvert(v)
       !         vxz = z - zvert(v)
       !         vwr = rvert(nextv) - rvert(v)
       !         vwz = zvert(nextv) - zvert(v)
       !
       !         cp = vxr*vwz - vxz*vwr
       !
       cp =    (((r-rvert(v)) * (zvert(nextv)-zvert(v))) &
            - ((z-zvert(v)) * (rvert(nextv)-rvert(v))))
       !
       !         There is a problem for points that should 
       !         lie on the boundary of the cell - i.e. that 
       !         are calculated based on the polygon corners and 
       !         which are mathematically on the polygon surface. 
       !         Numerically, these points can have a cross product
       !         which is close to zero but can vary to either side. 
       !         In order to consider these points in the cell - the 
       !         cross products are set to zero for values less than
       !         a specified limit. In this case the limit is set to 1.0e-7 
       !
       !         This value was determined by examining the range of cross 
       !         product values generated when sampling 50,000 points 
       !         calculated on a polygon with a scale size of 1.0m. 
       !         The maximum error cross product in this case was 6e-8.
       !
       !         D. Elder, Dec 13, 2006
       !
       if (abs(cp).lt.1.0e-7) cp = 0.0 

       if (v.eq.1.and.cp.ne.0.0) lastcp = cp
       !
       !         Look for change in sign of cp  
       !
       if ((lastcp * cp).lt.0.0) return 
       !
       if (cp.ne.0.0) lastcp = cp  
       !
    end do
    
    inpoly = .true.
    return
  end function inpoly



  subroutine distance( xstar,zstar,xalpha,zalpha,xbeta,zbeta, f )
    implicit none
    real*8 :: xstar,zstar,xalpha,zalpha,xbeta,zbeta,f
    real*8 :: m,x8,z8      
    !  to compute perp. distance from grid point (xstar,zstar) to line connecting
    !  points (xalpha,zalpha) and (xbeta,zbeta)

    m=(zbeta-zalpha)/(xbeta-xalpha) +1.e-12   ! slope of line between alpha,beta points
    x8=(zstar +xstar/m  -zalpha +m*xalpha)/( m +1/m )
    z8=zstar +(x8-xstar)/m
    f=sqrt( (xstar-x8)**2 +(zstar-z8)**2 )

    !      print 10,xstar,zstar,xalpha,zalpha,xbeta,zbeta
    !      print 10,m,x8,z8,f
10  format( 'm,x8,z8,f', 6g15.5 )

    return
  end subroutine distance


  subroutine fomega( xstar,zstar,xalpha,zalpha,xbeta,zbeta,omega )
    implicit none
    real*8 :: xstar,zstar,xalpha,zalpha,xbeta,zbeta,omega
    real *8 :: m      
    !  to compute cross product magnitude, "omega", of vector connecting
    !  grid point (xalpha,zalpha) to grid point (xbeta,zbeta); with vector connecting particle point (xstar,zstar) 
    !  to grid point (xalpha,zalpha)

    omega=( xbeta-xalpha )*( zstar-zalpha )-(xstar-xalpha )*( zbeta-zalpha )

    !      print 10,xstar,zstar,xalpha,zalpha,xbeta,zbeta,omega
10  format( '...omega', 8g15.5 )

    return
  end subroutine fomega


  subroutine write_cell_data(ik,ir)
    implicit none
    integer :: ik,ir
    ! write out the data for the specific cell for debugging purposes

    write(0,'(a,2i8,10(1x,g14.6))') 'DEBUG:',ik,ir,rs(ik,ir),zs(ik,ir),knbs(ik,ir),ktebs(ik,ir),ktibs(ik,ir),kvhs(ik,ir)

  end subroutine write_cell_data

end module oedge_plasma_interface
