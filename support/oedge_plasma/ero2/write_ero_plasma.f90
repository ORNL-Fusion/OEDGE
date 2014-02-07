program write_ero_plasma
  use error_handling
  use ero_plasma
  implicit none

  character*512 :: arg
  character*256 :: case_name,erospec_name
  integer :: nblocks
  character*8,parameter :: erospec_ext = '.erospec'
  integer :: in
  integer :: narg
  integer,external :: iargc
  integer :: errmsg_unit

  ! turn off error messaging
  errmsg_unit=-1

  ! read in the ero plasma specification file
  ! 
  ! it contains the following information
  !
  ! divimp plasma data file name
  ! divimp geometry data file name
  ! number of ero plasas to write out
  ! ero plasma specification blocks
  ! file name for storing the ero solution
  ! 4 vertices expected to be a rectangle
  ! R,Z offset to map the divimp machine coordinates to the ero simulation coordinates
  ! grid resolution in the two directions ... directions are defined by vectors and should be perpendicular
  ! interpolation option to be used 
  ! error message unit number
  ! unit conversion flag (do conversion here or in ERO?)
  !
  ! command line argument is the name of the data file to be read
  !


  ! open ero plasma command file

  narg = iargc()

  if (narg.ne.1) then 
     write(0,'(a)') 'Incorrect number of arguments specified:'
     write(0,'(a)') 'USAGE: program_name <ero_plasma_sepcification_file [str]>'
     stop
  endif

  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name' 
  !
  ! File name is the base case name without the .grd or .bgp file extensions
  !
  call getarg(1,arg)
  case_name = trim(arg)

  if (case_name.eq.'') then 
     ! no arguments specified
     write(0,'(a)') 'No arguments specified on call: This code requires the following input line'
     write(0,'(a)') '<program name>    ''input file name''  '
     stop
  endif

  !erospec_name = trim(case_name) // trim(erospec_ext)
  erospec_name = trim(case_name) 

  write(0,'(a,a,a,i8,a,i8)') 'ERO PLASMA: Case = ',trim(erospec_name)

  ! read data from the file
  ! DIVIMP plasma and geometry file names
  ! Number of ERO data blocks

  call read_ero_plasma_headers(erospec_name,nblocks)

  do in = 1,nblocks

     call read_ero_plasma_block
     call calc_ero_plasma(errmsg_unit)
     call output_ero_plasma 
     call calc_ero_surface
     call output_ero_surface

  end do



end program write_ero_plasma
