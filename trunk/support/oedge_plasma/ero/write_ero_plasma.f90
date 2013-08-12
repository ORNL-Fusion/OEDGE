program write_ero_plasma
  use error_handling
  implicit none

  character*512 :: arg
  character*256 :: case_name,erospec_name
  integer :: nblocks
  character*8,parameter :: erospec_ext = '.erospec'
  integer :: in
  integer :: narg
  integer,external :: iargc



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
  ! Command line is:  'file name'   interpolate_option   errmsg_unit
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


  ! Interpolate option = 0 = off (gives value in cell)
  !                    = 1 = on  (uses interpolation scheme)
  !call getarg(2,arg)
  !read(arg,*) interpolate_option

  !
  ! errmsg_unit specifies the unit number for error message output ... 
  ! A value less than 0 (e.g. -1 ) will suppress all error messages from the oedge_plasma_interface routines
  !
  !call getarg(3,arg)
  !read(arg,*) errmsg_unit

  erospec_name = trim(case_name) // trim(erospec_ext)

  write(0,'(a,a,a,i8,a,i8)') 'ERO PLASMA: Case = ',trim(erospec_name)

  ! read data from the file
  ! DIVIMP plasma and geometry file names
  ! Number of ERO data blocks

  call read_ero_plasma_headers(erospec_name,nblocks)


  do in = 1,nblocks

     call read_ero_plasma_block
     call calc_ero_plasma
     call output_ero_plasma 

  end do



end program write_ero_plasma
