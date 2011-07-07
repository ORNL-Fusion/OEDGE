program plasma_test
  use error_handling
  use oedge_plasma_interface
  implicit none

  character*512 :: arg
  character*256 :: case_name
  integer :: interpolate_option,errmsg_unit
  integer :: ierr
  integer :: ir,iz,nr,nz
  real*8 :: r_offset, z_offset

  real*8 :: rmin,rmax,zmin,zmax 
  real*8 :: dr,dz
  real*8 :: rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt

  !
  !
  !
  ! Three steps to using the interface
  ! 1) Initialization  -   init_oedge_plasma
  ! 2) Get plasma data -   get_oedge_plasma
  ! 3) Clean up        -   close_oedge_plasma 
  !
  ! Clean up is done when the GET routine will no longer be called - it frees up all allocated arrays
  !
  !

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

  write(0,'(a,a,a,i8,a,i8)') 'OEDGE PLASMA: Case = ',trim(case_name),'  Interpolate Option =',interpolate_option,'    Error message output unit =',errmsg_unit


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

  call init_oedge_plasma(case_name,r_offset,z_offset,interpolate_option,errmsg_unit,ierr)

  if (ierr.ne.0) then 
     call errmsg('Problem loading oedge plasma',ierr)
     stop 'plasma not loaded'
  endif


  ! Define sample space - call and report plasma conditions
  ! Values chosen here are arbitrary for testing purposes

  nr = 100
  nz = 100

  rmin = 1.4
  rmax = 1.6

  zmin = -1.37
  zmax = -1.17

  dr = (rmax-rmin)/nr
  dz = (zmax-zmin)/nz

  write(6,*) 
  write(6,'(a,i10)') 'Ordered by R: Interpolation =',interpolate_option
  write(6,*) 

  do ir = 0,nr

     rt = rmin + ir * dr

     do iz = 0,nz

        zt = zmin + iz * dz

        call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

        if (ierr.eq.0) then 
           write(6,'(2i8,20(1x,g18.8))') ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
        else
           write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
           call errmsg('PLASMA TEST',error_message_data)
        endif

     end do
     write(6,*) 

  end do

  write(6,*) 
  write(6,'(a,i10)') 'Ordered by Z',interpolate_option
  write(6,*) 


  do iz = 0,nz

     zt = zmin + iz * dz

     do ir = 0,nr

        rt = rmin + ir * dr

        call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

        if (ierr.eq.0) then 
           write(6,'(2i8,20(1x,g18.8))') ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
        else
           write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
           call errmsg('PLASMA TEST',error_message_data)
        endif

     end do
     write(6,*) 

  end do




  ! Clean up when done calling OEDGE plasma interface

  call close_oedge_plasma



end program plasma_test
