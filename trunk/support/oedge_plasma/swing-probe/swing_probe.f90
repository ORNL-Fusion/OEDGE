program swing_probe
  use error_handling
  use oedge_plasma_interface
  implicit none

  character*512 :: arg
  character*512 :: case_name,gridfilename,plasmafilename,spfilename
  character*512 :: ext_coord_file,line

  integer :: interpolate_option,errmsg_unit,extrapolate_option
  integer :: ierr

  real*8 :: r_offset, z_offset

  real*8 :: ne,te,ti,vb,ef,psin,btot,br,bz,bt
  real*8 :: psin_ext


  integer :: in,nsteps,iunit,ios,extunit

  real*8 :: rcen,zcen,radius,r,z
  real*8 :: start_theta, end_theta, dtheta,theta
  real*8,parameter  :: PI = 3.1415928
  real*8,parameter :: RADDEG=57.29577952,DEGRAD=1.745329252D-02

  logical :: external_coords


  !
  !
  !
  ! Three steps to using the interface
  ! 1) Initialization  -   set_oedge_plasma_opts
  ! 2) Get plasma data -   load_oedge_plasma
  ! 3) Clean up        -   close_oedge_plasma 
  !
  ! Clean up is done when the GET routine will no longer be called - it frees up all allocated arrays
  !
  !
  ! Initialize
  !
  external_coords = .false.
  ios = 0

  ! Initialization for analytic swing probe path
  ! set up the coordinate calculation
  ! do 180 degree sweep - from down to up - use a specifiable degree resolution 
  rcen =  0.977
  zcen = -0.8285
  radius = 0.201
  dtheta = 0.25
  start_theta = -90.0
  end_theta = 90.0



  ! Load command line input - case name, interpolation and error messaging options
  ! 
  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name'   interpolate_option   errmsg_unit
  !
  ! File name is the base case name without the .grd or .bgp file extensions
  !
  call getarg(1,arg)
  case_name = trim(arg)

  if (case_name.eq.'') then 
     ! no arguments specified
     write(0,'(a)') 'No arguments specified on call: This code requires the following input line'
     write(0,'(a)') '<program name>    ''input file case name''   interpolate_option    errmsg_unit'
     stop
  endif


  ! Interpolate option = 0 = off (gives value in cell)
  !                    = 1 = on  (uses interpolation scheme)
  call getarg(2,arg)
  read(arg,*) interpolate_option

  !
  ! errmsg_unit specifies the unit number for error message output ... 
  ! A value less than 0 (e.g. -1 ) will suppress all error messages from the oedge_plasma_interface routines
  !

  call getarg(3,arg)
  read(arg,*) errmsg_unit

  write(0,'(a,a,a,i8,a,i8)') 'SWING PROBE: OEDGE PLASMA: Case = ',trim(case_name),'  Interpolate Option =',interpolate_option,'    Error message output unit =',errmsg_unit


  call getarg(4,arg)

  if (arg(1:2).eq.'-e') then 
     external_coords = .true.
     call getarg(5,arg)
     ext_coord_file = trim(arg)
     write(0,'(a,a)') 'External coordinate file:',trim(ext_coord_file)
  endif


  extrapolate_option = 0

  ! 
  ! R_offset and Z_offset specify a coordinate system shift to be applied to coordinates passed into the OEDGE plasma routine. 
  ! e.g. If R_offset = 1.5 and Z_offset = -1.366 then a coordinate of 0,0 passed into the OEDGE plasma routine will be mapped to 1.5,-1.366 
  ! This makes it easy to deal with different coordinate systems with different origins. 
  !
  ! However, the OEDGE code works in MKS units ... thus distances and geometry are in meters. If the calling code uses a different coordinate scale (e.g. cm) then 
  ! appropriate conversions must be made before the call to get the plasma conditions. 
  !

  r_offset = 0.0
  z_offset = 0.0

  call set_oedge_plasma_opts(r_offset,z_offset,interpolate_option,extrapolate_option,errmsg_unit)

  gridfilename = trim(case_name)//'.grd'
  plasmafilename = trim(case_name)//'.bgp'


  call load_oedge_data(gridfilename,plasmafilename,ierr)

  if (ierr.ne.0) then 
     call errmsg('Problem loading oedge plasma',ierr)
     stop 'plasma not loaded'
  endif


  ! get a unit number and open the output file name
  spfilename = trim(case_name)//'.swing_probe'

  ! open output file
  call find_free_unit_number(iunit)
  open(unit=iunit,file=trim(spfilename),iostat=ios)

  if (external_coords) then 

     write(iunit,'(a,a,3x,3(1x,g18.8))') ' Swing probe results: EXTERNAL_COORDS: ',trim(case_name)

  else

     write(iunit,'(a,a,3x,3(1x,g18.8))') ' Swing probe results: ',trim(case_name),rcen,zcen,radius

  endif

  write(iunit,'(a)') ' IN      THETA(degrees)    R(m)          Z(m)         PSIn        '//&
       & 'Te(eV)     Ti(eV)    ne(m-3)    ef(V/m)    vb(m/s)      Btot(T)     B_r    B_z   B_tor'


  if (external_coords) then 

     call find_free_unit_number(extunit)
     open(unit=extunit,file=trim(ext_coord_file),iostat=ios)

     if (ios.ne.0) then
        call errmsg('ERROR: Problem opening external coordinates file: code =',ios)
        stop 'STOP: Error opening external coordinates'
     endif

     ! Discard initial comment line in external coordinate file
     read(extunit,'(a)') line


     ios = 0
     in = 0
     theta = 0.0

     do while (ios.eq.0) 

        read(extunit,*,iostat=ios) r,z,psin_ext

        if (ios.eq.0) then 

           call get_oedge_plasma(r,z,ne,te,ti,vb,ef,psin,btot,br,bz,bt,ierr)

           if (ierr.eq.0.or.ierr.eq.2) then     ! ierr=0 data found ... ierr=2 extrapolated data returned
              in = in + 1
              write(iunit,'(1i8,20(1x,g18.8))') in,theta*raddeg,r,z,psin_ext,te,ti,ne,ef,vb,btot,br,bz,bt
              !
              !     Off grid data is not an error in this case
              !
              !     elseif (ierr.eq.1) then 
              !        write(error_message_data,'(a,2i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,in,r,z
              !        call errmsg('PLASMA TEST',error_message_data)
           endif



        endif


     end do


  else
     ! Swing probe arc described analytically


     nsteps = int((end_theta-start_theta)/dtheta) +1 

     ! convert to radian

     dtheta = dtheta * degrad
     start_theta = start_theta * degrad
     end_theta = end_theta * degrad


     do in = 0,nsteps

        theta = start_theta + in * dtheta

        r = rcen + radius * cos(theta)
        z = zcen + radius * sin(theta)

        call get_oedge_plasma(r,z,ne,te,ti,vb,ef,psin,btot,br,bz,bt,ierr)

        if (ierr.eq.0.or.ierr.eq.2) then     ! ierr=0 data found ... ierr=2 extrapolated data returned
           write(iunit,'(1i8,20(1x,g18.8))') in,theta*raddeg,r,z,psin,te,ti,ne,ef,vb,btot,br,bz,bt
           !
           !     Off grid data is not an error in this case
           !
           !     elseif (ierr.eq.1) then 
           !        write(error_message_data,'(a,2i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,in,r,z
           !        call errmsg('PLASMA TEST',error_message_data)
        endif

     end do


  endif



  ! Clean up when done calling OEDGE plasma interface

  call close_oedge_plasma

  ! close output file

  close(iunit)

end program swing_probe
