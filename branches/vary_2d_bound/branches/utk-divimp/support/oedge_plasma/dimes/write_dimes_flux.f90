program write_dimes_flux

  use error_handling
  use allocate_arrays
  use utilities
  use oedge_plasma_interface
  implicit none


  character*512 :: arg
  character*256 :: case_name,gridfilename,plasmafilename,fluxfilename
  integer :: interpolate_option,errmsg_unit
  integer :: ierr
  integer :: ir,iz,nr,nz
  integer :: ounit
  character*256 :: outfilename,matlab_outfilename

  real*8 :: r_offset, z_offset

  real*8 :: rmin,rmax,zmin,zmax 
  real*8 :: dr,dz
  real*8 :: rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,psin

  integer :: in,nizs
  real*8, allocatable :: flux(:),energy(:)

  
  integer :: narg
  integer,external :: iargc



  !
  ! This program was developed from the plasma_test sample code. This code writes out a flux file for a grid/wall
  ! covering a specified region of the surface. Jeff can read this file into his code. 
  !

  
  !
  ! Three steps to using the interface
  ! 1) Initialization  -   init_oedge_plasma
  ! 2) Get plasma data -   get_oedge_plasma
  ! 3) Interpolate and output fluxes  
  ! 4) Clean up        -   close_oedge_plasma 
  !

    ! Load command line input - case name, interpolation and error messaging options
  ! 

  narg = iargc()

  if (narg.ne.3) then 
     write(0,'(a)') 'Incorrect number of arguments specified:'
     write(0,'(a)') 'USAGE: program_name <case_name [str]>  <interpolate_option [int]>  <error_message_unit [int] less than 0 is off >'
     stop
  endif

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


  ! Interpolate option = 0 = off (gives value on surface element)
  !                    = 1 = on  (uses interpolation scheme)
  call getarg(2,arg)
  read(arg,*) interpolate_option

  !
  ! errmsg_unit specifies the unit number for error message output ... 
  ! A value less than 0 (e.g. -1 ) will suppress all error messages from the oedge_plasma_interface routines
  !
  call getarg(3,arg)
  read(arg,*) errmsg_unit

  write(0,'(a,a,a,i8,a,i8)') 'OEDGE IMPURITY FLUX: Case = ',trim(case_name),'  Interpolate Option =',interpolate_option,'    Error message output unit =',errmsg_unit


  ! Need to specify offsets to remap flux data to surface near DIMES
  
  ! 170844
  ! DIMES tungsten spot exposures
  !
  ! grid strike point = 1.48063
  ! actual averaged strike point = 1.47903
  ! DIMES center = 1.489
  ! Map DIMES center to grid relative to strike point
  ! 1.489 - 1.47903 = 9.97e-3
  ! R value of DIMES center relative to grid = 1.48063+9.97e-3 = 1.4906
  
  r_offset =  1.4906
  z_offset = -1.25

  ! 170851
  ! DIMES tungsten spot exposures
  !
  ! grid strike point = 1.46925
  ! actual averaged strike point = 1.468576
  ! DIMES center = 1.489
  ! Map DIMES center to grid relative to strike point
  ! 1.489 - 1.468576 = 2.0424e-2
  ! R value of DIMES center relative to grid = 1.46925 + 2.0424e-2 = 1.489674
  
  !r_offset =  1.489674
  !z_offset = -1.25
  
  write(0,'(a,f15.6,a,f15.6)') 'OFFSETS: R_OFFSET=',r_offset, ' Z_OFFSET=',z_offset
  
  gridfilename = trim(case_name)//'.grd'
  plasmafilename = trim(case_name)//'.bgp'
  fluxfilename = trim(case_name)//'.impdep'

  ! offsets are set globally in the oedge_plasma module
  call set_oedge_plasma_opts(r_offset,z_offset,interpolate_option,1,errmsg_unit)

  call load_oedge_data(gridfilename,plasmafilename,ierr)
 
  if (ierr.ne.0) then 
     call errmsg('Problem loading oedge plasma',ierr)
     stop 'plasma not loaded'
  endif

  call load_resolved_deposition_data(fluxfilename,ierr)

  if (ierr.ne.0) then 
     call errmsg('Problem loading oedge flux data',ierr)
     stop 'flux data not loaded'
  endif
  


  ! Define sample space - call and report plasma conditions
  ! Values chosen here are arbitrary for testing purposes

  outfilename = trim(case_name)//'.impurity.fluxes.dat'

  call find_free_unit_number(ounit)

  open(ounit,file=trim(outfilename),status='unknown',form='formatted')

  write(0,'(a)') 'Writing fluxes to file:'//trim(outfilename)

  nr = 120
  !nz = 0
  nz = 60

  rmin = -0.06
  rmax = 0.06

  !zmin = 0.0
  !zmax = 0.0
  zmin = 0.0
  zmax = 0.0

  if ((rmax.eq.rmin).or.nr.eq.0) then 
     dr=0.0
  else
     dr = (rmax-rmin)/nr
  endif
  
  if ((zmax.eq.zmin).or.nz.eq.0) then 
     dz = 0.0
  else
     dz = (zmax-zmin)/nz
  endif

  !call set_debug_code(.true.)


  ! allocate storage for flux and energy arrays
  ! 
  nizs = get_deposition_charge_states()

  call allocate_array(flux,nizs,'FLUX',ierr)
  call allocate_array(energy,nizs,'ENERGY',ierr)



  write(ounit,*) 
  write(ounit,'((1x,a7),2x,2(a10,3x),100(3x,a3,i8,5x,2(7x,a10,2x)))') 'R_INDEX','R(m)','Z(m)',('IZ=',iz,'FLUX','ENERGY',iz=1,nizs)
  write(ounit,*) 

  
  ! This is designed just to work for area near DIMES - generating code for arbitrary sections of wall will require more detailed work

  
  do ir = 0,nr

     rt = rmin + real(ir) * dr
     !write(0,'(a,i8,10(1x,g18.8))') 'R:',ir,rt,rmin,dr,rmin + ir*dr

     ! fixed Z value for DIMES target
     zt = zmin

     ! offsets are factored in inside the oedge_plasma module
     call interpolate_deposition_data(rt,zt,flux,energy,ierr)

     if (ierr.eq.0) then 
        write(ounit,'(i8,2x,,2(1x,g12.5),100(6x,i8,5x,2(1x,g18.8)))') ir,rt,zt,(iz,flux(iz),energy(iz),iz=1,nizs)
     else
        write(ounit,'(i8,2x,2(1x,g12.5),100(1x,i8,2(1x,g18.8)))') ir,rt,zt,(iz,flux(iz),energy(iz),iz=1,nizs)
        write(error_message_data,'(a,2i8,2(1x,g12.5))') 'ERROR GETTING FLUXES:',ierr,ir,rt,zt
        call errmsg('WRITE FLUX',error_message_data)
        !write(6,'(a)') error_message_data

     endif

  end do

  
  if (allocated(flux)) deallocate(flux)
  if (allocated(energy)) deallocate(energy)

  ! Clean up when done calling OEDGE plasma interface

  call close_oedge_plasma

  close(ounit)





end program write_dimes_flux
