!
! jdemod - Sept 8/2016
!
! The following code is from Sterling Smith. It implements a nuclear approach
! to creating netcdf files rather than the standard work flow. This should 
! be appropriate to adding netcdf output to OEDGE.
! 
! The file name has been changed to f90 from F90 so that CPP is not run. The 
! ifdef block containing the test code has been commented out
!
! Compilation of this in OEDGE requires access to the netcdf and hdf5 libraries.
! It also needs the netcdf.mod file. This file was NOT included in the standard
! netcdf, netcdf-devel, netcdf-fortran and netcdf-fortran-devel packages for
! RHEL ... it should have been in the netcdf-fotran-devel package as far as I 
! could tell but was not. netcdf-fortran 4.2 source had issues with the 
! autoconf and automake functionality on RHEL. 4.4.4 had issues with the base
! netcdf libraries under RHEL being incompatible. So I have added recent versions
! of hdf5, netcdf and netcdf-fortran sources to divimp/libsrc and I have 
! compiled these libraries into ~/divimp/local. Makefiles will need to be 
! adjusted to obtain the required libraries and netcdf.mod from lib and include
! directories respectively.  
!
! The routine uses the error_handling module for issuing error messages
!
! OVERVIEW
!
! This nc_utils module is meant to be an easy interface to netcdf
! 
! While the typical netcdf workflow is
! 
! nc_open
! for all dimensions
!  nc_define_dimension
! endfor
! for all variables
!  nc_define_variable
!  nc_put_attribute1
!  nc_put_attribute2
!  nc_put_attribute3
! endfor 
! nc_switch_to_data_mode
! for all variables
!  nc_put_variable
! endfor
! nc_close
!
! Instead the nc_utils_generic module provides the following functions:
! write_nc, read_nc, open_nc_file, close_nc_file
!
! which can be used to read and write data from/to a netcdf file in a single function call, optionally including
! the indication of dimensions and the long_name and units attributes.
!
! 
! Routines:
!
! <> denote optional arguments
!
!
!  function open_nc_file(filename,<mode>,<check_zero_flag>,<debug>) result (ierr)    
!       filename - character - file to be opened
!       mode - NC_WRITE or NC_READONLY - default NC_READONLY
!       check_zero_flag - flags whether to check for zero values before storing to netcdf file - default TRUE
!       debug - logical - sets the verbose flag - default FALSE
!       ierr     - integer - error return code
!       [Opens netcdf file for access from module routines]
!
!  function close_nc_file result (ierr)
!       ierr - integer - error return code
!       [Closes the netcdf file opened for access in the module]
!
!
!  Routines to read and write various data types
!
!  Character data
!
!  write_nc(var_name,var,<log_name>,<units>) result(ier)
!  read_nc(var_name,var,<log_name>,<units>) result (ier)
!
!  I4, R4 and R8 data  
!  - scalar and arrays 
!
!  Scalars:
!  write_nc(var_name,vars,<long_name>,<units>) result(ier)
!  Arrays: 
!  write_nc(var_name,vars,dim_names,number_of_dimensions,<long_name>,<units>) result(ier)
!
!  Scalars:
!  read_nc(var_name,vars,<long_name>,<units>) result(ier)
!  Arrays: 
!  read_nc(var_name,vars,<long_name>,<units>) result(ier)
!
!  var_name - character - netcdf name of the variable
!  vars - variable to be used to load or store the data
!  dim_names - character array - netcdf identifiers for each of the dimension names
!  number_of_dimensions - integer - dimensionality of the variable to be read or written ... 0=scalar, 1=1d, 2=2d, 3=3d etc 
!  long_name - character - long name of the variable
!  units - character - name for units
!
!  Note: character arrays have yet to be implemented
!
! 
!
! TESTS
! 
! The unit tests for this module are provided at the end of file and can be compiled with
!
! The testing code is used to test both versions of the netcdf nc_utils. 
! 
! $FC -I$NETCDF_DIR/include/ -L$NETCDF_DIR/lib -lnetcdff -lnetcdf -o tnc divimp_types.f90 error_handling.f90 nc_utils_generic_v1.f90 nc_utils_generic.f90 run_nc_test.f90
! 
! which produces the executable 'tnc', which can be run to produce the output
! 
! > ./tnc
! PASS: Write a text string
! PASS: Write another test string
! PASS: Write R8 scalar
! PASS: Not write R8 scalar 0
! PASS: Not write R8 scalar 0
! PASS: Write R4 scalar
! PASS: Not write scalar R4 0
! PASS: Write I4scalar
! PASS: Not write I4 scalar 0
! PASS: Write R8 1d array
! PASS: Not write R8 1d array of zeros
! PASS: Update R8 1d array with value*2
! PASS: Not overwite R8 1d array with wrong dimension name
! PASS: Not write R8 1d array with wrong dimension length
! PASS: Not write R8 1d array with wrong size
! PASS: Write R4 1d array
! PASS: Not write R4 1d array of zeros
! PASS: Update R4 1d array with value*2
! PASS: Not overwite R4 1d array with wrong dimension name
! PASS: Not write R4 1d array with wrong dimension length
! PASS: Not write R4 1d array with wrong size
! PASS: Write I4 1d array
! PASS: Not write I4 1d array of zeros
! PASS: Update I4 1d array with value*2
! PASS: Not overwite I4 1d array with wrong dimension name
! PASS: Not write I4 1d array with wrong dimension length
! PASS: Not write I4 1d array with wrong size
! PASS: Write R8 2d array
! PASS: Not write R8 2d array of zeros
! PASS: Update R8 2d array with value*2
! PASS: Not overwite R8 2d array with wrong dimension name
! PASS: Not overwrite R8 2d array with wrong dimension length
! PASS: Not write R8 2d array with wrong size
! PASS: Write R4 2d array
! PASS: Not write R4 2d array of zeros
! PASS: Update R4 2d array with value*2
! PASS: Not overwite R4 2d array with wrong dimension name
! PASS: Not overwrite R4 2d array with wrong dimension length
! PASS: Not write R4 2d array with wrong size
! PASS: Write I4 2d array
! PASS: Not write I4 2d array of zeros
! PASS: Update I4 2d array with value*2
! PASS: Not overwite I4 2d array with wrong dimension name
! PASS: Not overwrite I4 2d array with wrong dimension length
! PASS: Not write I4 2d array with wrong size
! PASS: Not write R8 3d array of zeros
! PASS: Write R8 3d array
! PASS: Update R8 3d array with value*2
! PASS: Not overwite R8 3d array with wrong dimension name
! PASS: Not write R8 3d array with wrong dimension length
! PASS: Not write R8 3d array with wrong size
! PASS: Not write R4 3d array of zeros
! PASS: Write R4 3d array
! PASS: Update R4 3d array with value*2
! PASS: Not overwite R4 3d array with wrong dimension name
! PASS: Not write R4 3d array with wrong dimension length
! PASS: Not write R4 3d array with wrong size
! PASS: Not write I4 3d array of zeros
! PASS: Write I4 3d array
! PASS: Update I4 3d array with value*2
! PASS: Not overwite I4 3d array with wrong dimension name
! PASS: Not write I4 3d array with wrong dimension length
! PASS: Not write I4 3d array with wrong size
! PASS: Read a text string
! PASS: Read another test string
! PASS: Test string does not exist
! PASS: Test string too short
! PASS: READ R8 scalar
! PASS: READ R4 scalar
! PASS: READ I4 scalar
! PASS: READ scalar not found
! PASS: Read 1d r8 array
! PASS: Assign zeroes to R8 1D array not found
! PASS: Not read 1d r8 array with wrong dimension length
! PASS: Read 1d r4 array
! PASS: Assign zeroes to R4 1D array not found
! PASS: Not read 1d r4 array with wrong dimension length
! PASS: Read 1d i4 array
! PASS: Assign zeroes to I4 1D array not found
! PASS: Not read 1d i4 array with wrong dimension length
! PASS: Read 2d r8 array
! PASS: Assign zeroes to R8 2D array not found
! PASS: Not read array - test_2d/10_arr_r8 - not in database
! PASS: Read 2d r4 array
! PASS: Assign zeroes to R4 2D array not found
! PASS: Not read 2d r4 array with wrong dimension length
! PASS: Read 2d i4 array
! PASS: Assign zeroes to I4 2D array not found
! PASS: Not read 2d i4 array with wrong dimension length
! PASS: Read 3d r8 array
! PASS: Assign zeroes to R8 3D array not found
! PASS: Not read 3d r8 array with wrong dimension length
! PASS: Read 3d r4 array
! PASS: Assign zeroes to R4 3D array not found
! PASS: Not read 3d r4 array with wrong dimension length
! PASS: Read 3d i4 array
! PASS: Assign zeroes to I4 3D array not found
! PASS: Not read 3d i4 array with wrong dimension length
!
! SHORT EXAMPLE
!
!
! use nc_utils_generic
! use netcdf, ONLY: nf90noerr
!
! The netcdf file is opened by the user
! ier=open_nc_file('test.inc',mode=NC_WRITE)
! The desired array is written
! ier = write_nc('var1', 1,'This is var1', 'meter', (/ 'dim1' /), size(var1),v1d=var1)
! OR
! ier = write_nc('var1', 1,long_name='This is var1', units='meter', dim_names=(/ 'dim1' /), nd=size(var1),v1d=var1)
! 
! Test return code for error
! if (ier.ne.0) WRITE (*,*) 'There was a problem with var1'      ! alternatively could add use netcdf, ONLY:nf90noerr  ... and test against that return value
!   
! ier = close_nc_file()
!



MODULE nc_utils_generic

  use error_handling
  use debug_options

  implicit none

  PRIVATE
  PUBLIC ::  write_nc,read_nc,open_nc_file,test_nc_utils,close_nc_file,NC_WRITE,NC_READONLY
  !    public ::  open_nc_file,close_nc_file

  CHARACTER(len=1024) :: err_msg ! Holder of error messages
  INTEGER :: nc_id
  LOGICAL :: verbose = .false.

  !
  ! Define fixed attribute names
  !

  character*9, parameter :: longname='long_name'
  character*5, parameter :: unitsname='units'

  ! will not store NEW variables with a zero value to the database
  ! Note that when reading data ... IF zero_check is true then all variables that 
  ! are not found will be assigned a zero value. Probably do the same without the zero_check flag
  ! but may issue additional error messages
  LOGICAL :: zero_check = .true.


  interface write_nc
     module procedure write_char_data, &
         & write_sc_i4_data, write_1d_i4_data, write_2d_i4_data, write_3d_i4_data, write_4d_i4_data, &
         & write_sc_r4_data, write_1d_r4_data, write_2d_r4_data, write_3d_r4_data, write_4d_r4_data, &
         & write_sc_r8_data, write_1d_r8_data, write_2d_r8_data, write_3d_r8_data, write_4d_r8_data
  end interface write_nc

  interface read_nc
     module procedure read_char_data, &
         & read_sc_i4_data, read_1d_i4_data, read_2d_i4_data, read_3d_i4_data, read_4d_i4_data, &
         & read_sc_r4_data, read_1d_r4_data, read_2d_r4_data, read_3d_r4_data, read_4d_r4_data, &
         & read_sc_r8_data, read_1d_r8_data, read_2d_r8_data, read_3d_r8_data, read_4d_r8_data
  end interface read_nc

  integer, parameter :: NC_READONLY = 0
  integer, parameter :: NC_WRITE = 1



CONTAINS

  function open_nc_file(filename,mode,check_zero_flag,debug) result(ierr)
    ! opens the nc file ... sets nc_id and verbose output/debug
    use error_handling
    use netcdf
    implicit none

    character*(*) :: filename
    logical,optional :: debug,check_zero_flag
    integer,optional :: mode
    integer :: ierr

    integer :: mode_val,create_mode

    verbose = .FALSE.
    zero_check = .TRUE.

    if (present(debug)) verbose = debug
    if (present(check_zero_flag)) zero_check = check_zero_flag

    ! direct error messages only to stderr ... default is both stderr and stdout i.e. fort.0 and fort.6
    call set_errmsg_units(0,-1,-1)

    ! set default mode_val to READONLY/NOWRITE in case the mode passed in isn't specified
    mode_val = NF90_NOWRITE

    ! READ-ONLY
    if (present(mode)) then 
       if (mode.eq.NC_READONLY) then 
          mode_val = NF90_NOWRITE
       elseif (mode.eq.NC_WRITE) then 
          mode_val = NF90_WRITE
       endif
    endif

    ! Try to open the database file ... 
    ierr =  nf90_open(trim(filename),mode_val,nc_id)

    ! If this fails then check read/write request and either issue error or create

    if (ierr .ne. nf90_noerr) then 

       if (mode_val.eq.NF90_WRITE) then
          ! file does not exist but it is to be written ... run nf90_create
          ! Just in case the error was NOT that the file was missing ... specify NOCLOBBER
          ! so we don't overwrite an existing database
          create_mode = NF90_NOCLOBBER          

          ierr = nf90_create(trim(filename),create_mode,nc_id)

          if (ierr.ne.nf90_noerr) then 
             call errmsg('OPEN_NC_FILE: ERROR OPENING:'//trim(filename)//': FILE FOR WRITING',ierr)
             stop 'OPEN_NC_FILE_WRITE_ERROR'
          endif

       elseif (mode_val.eq.NF90_NOWRITE) then 
          ! database was expected but could not be opened ... issue an error message and quit
          call errmsg('OPEN_NC_FILE: ERROR OPENING:'//trim(filename)//': FILE FOR READING',ierr)
       endif

    endif

    ! At this point the netcdf database file has been opened either for reading or writing OR the code exited with an error message
    return
  end function open_nc_file


  function close_nc_file () result (ierr)
    use error_handling
    use netcdf, ONLY : nf90_close,nf90_noerr
    implicit none
    integer :: ierr

    ! this routine closes the nc file opened in this module 

    ierr = nf90_close(nc_id)
    if (ierr.ne.nf90_noerr) then 
       if (verbose) call errmsg('NC_UTILS_GENERIC: CLOSE_NC_FILE: Error closing file = ',ierr)
    endif

    call reset_errmsg_units

  end function close_nc_file


  !-------------------------------------------------
  ! Utility functions - internal

  FUNCTION handle_nf90_error(ier_in) RESULT(ier)
    USE netcdf, ONLY : NF90_NOERR, nf90_close
    use error_handling
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ier_in
    INTEGER :: ier, ier_close
    ! error handler ... just writes an error message if a problem comes up

    ier = ier_in
    IF ( ier .ne. NF90_NOERR ) THEN
       IF ( verbose ) call errmsg("NC_UTILS_GENERIC:"//trim(err_msg),ier)
       !WRITE(*,*) TRIM(err_msg)
    ENDIF
  END FUNCTION handle_nf90_error

  FUNCTION check_existing_variable(var_name, var_type, var_id, dim_names, dim_vals) RESULT(ier)
    USE netcdf, ONLY : nf90_inquire_variable, NF90_MAX_NAME, &
         & NF90_MAX_VAR_DIMS, nf90_noerr, nf90_inquire_dimension, nf90_inq_dimid
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_type, var_id
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), DIMENSION(:),optional ::  dim_vals
    character*(*),intent(in),optional :: dim_names(:)
    CHARACTER(len=NF90_MAX_NAME) :: existing_var_name,dim_name
    integer :: dim_len
    INTEGER ::  existing_type, existing_ndims, ind, ier, dim_id
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: existing_dimids

    ier = nf90_inquire_variable(nc_id, var_id, existing_var_name,existing_type, existing_ndims, existing_dimids)

    !write(0,*) 'checking variable:',var_id

    ! is variable found
    IF ( ier .ne. NF90_NOERR ) THEN
       WRITE(err_msg,*) 'Unable to check existing variable ',var_name
       RETURN
    ENDIF

    ! check if the variable name matches - however, this should always be a match since the 
    ! variable name was used to look up the variable in the first place
    IF ( TRIM(existing_var_name) .ne. TRIM(var_name) ) THEN
       WRITE(err_msg,*) 'Variable name is different ',TRIM(var_name),TRIM(existing_var_name)
       ier = NF90_NOERR + 1
       RETURN
    ENDIF

    ! check if stored type matches expected type
    IF ( existing_type .ne. var_type ) THEN
       WRITE(err_msg,*) 'Type mismatch for variable ',TRIM(var_name)
       ier = NF90_NOERR + 2
       RETURN
    ENDIF

    if (present(dim_names)) then 
       ! check if rank of variable matches the expected rank
       IF ( existing_ndims .ne. SIZE(dim_names) ) THEN
          WRITE(err_msg,*) 'Dimension mismatch for variable ',TRIM(var_name)
          ier = NF90_NOERR + 3
          RETURN
       ENDIF

       ! check if the DIM IDs of the dimensions match  ... NOTE: this is different from checking if the actual sizes match
       DO ind = 1, SIZE(dim_names)
          ! check to see if variable has dimensions with different names than existing
          ier = nf90_inq_dimid(nc_id, trim(dim_names(ind)), dim_id)

          IF ( dim_id .ne. existing_dimids(ind) ) THEN
             WRITE(err_msg,*) 'Attempt to write variable ',TRIM(var_name), ' with existing dimension ',ind,' with ID=',existing_dimids(ind),' different from new dimension ID=',dim_id
             ier = NF90_NOERR + 4
             RETURN
          ENDIF
       END DO
    endif

    if (present(dim_vals)) then 

       ! check if rank of variable matches the expected rank
       IF ( existing_ndims .ne. SIZE(dim_vals) ) THEN
          WRITE(err_msg,*) 'Dimension mismatch for variable ',TRIM(var_name)
          ier = NF90_NOERR + 3
          RETURN
       ENDIF

       ! check the values of the dimensions
       DO ind = 1, SIZE(dim_vals)
          ier = nf90_inquire_dimension(nc_id,existing_dimids(ind),dim_name,dim_len)
          IF ( dim_len .ne. dim_vals(ind) ) THEN
             WRITE(err_msg,*) 'Attempt to read/write variable ',TRIM(var_name), ' with different dimension sizes'
             ier = NF90_NOERR + 4
             RETURN
          ENDIF
       END DO
    endif

    ier = NF90_NOERR
    RETURN
  END FUNCTION check_existing_variable

  FUNCTION check_existing_dimension(var_name,dim_id, dim_name, dim_len) RESULT(ier)
    USE netcdf, ONLY : NF90_MAX_NAME, nf90_noerr, nf90_inquire_dimension
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim_id, dim_len
    CHARACTER(len=*), INTENT(IN) :: dim_name
    character(len=*), intent(in) :: var_name   ! include var_name to make error identification easier
    INTEGER :: existing_len, ier
    CHARACTER(len=NF90_MAX_NAME) :: existing_name

    WRITE(err_msg,*) 'Error inquiring dimension ', TRIM(dim_name)

    ier = handle_nf90_error(nf90_inquire_dimension(nc_id, dim_id, existing_name, existing_len))
    IF ( ier .ne. nf90_noerr ) RETURN

    !write(0,'(5a,i8,3a,i8,a)') 'Check_existing_dimension:',trim(var_name),':',trim(dim_name),':',dim_len,':',trim(existing_name),':',existing_len,':'


    ! check dimension name
    IF (TRIM(dim_name) .ne. TRIM(existing_name) ) THEN
       ier = nf90_noerr + 1
       WRITE(err_msg,*) 'Dim id ',dim_id, ' name mismatch ', TRIM(dim_name), ' ', TRIM(existing_name),' for var=',trim(var_name)
       RETURN
    ENDIF

    ! check value of dimension
    IF (dim_len .ne. existing_len) THEN
       ier = nf90_noerr + 2
       WRITE(err_msg,*) 'Dimension ', TRIM(dim_name), ' is already of size ', existing_len, ' not ', dim_len,' for var=',trim(var_name)
    ENDIF

  END FUNCTION check_existing_dimension

  !  SUBROUTINE convert_mks(units, new_units, fac)
  !    CHARACTER(len=*), INTENT(IN) :: units
  !    CHARACTER(len=MAX_LEN_UNITS), INTENT(OUT) :: new_units
  !    REAL*8, INTENT(OUT) :: fac
  !    CHARACTER(len=LEN_TRIM(units)) :: tunits
  !    tunits = TRIM(units)
  !    IF ( tunits .eq. 'W/cm^3') THEN
  !      fac = 1e6
  !      new_units = 'watts/meter^3'
  !    ELSE IF (tunits .eq.
  !  END SUBROUTINE

  FUNCTION add_units_long_name(var_id, units, long_name, fac) RESULT(ier)
    ! jdemod
    ! fac is part of an unimplemented automated unit conversion utility
    ! will leave it in for now but not really needed
    USE netcdf, ONLY : nf90_put_att, nf90_noerr
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_id
    CHARACTER(len=*), INTENT(IN) :: units, long_name
    REAL*8, INTENT(OUT) :: fac
    INTEGER :: ier
    !    CALL convert_mks(units, new_units, fac)

    ier = switch_to_define_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error adding units attribute',TRIM(units)
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id,unitsname,TRIM(units)))
    IF ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error adding long_name attribute',TRIM(long_name)
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id, longname,TRIM(long_name)))

    RETURN
  END FUNCTION add_units_long_name


  FUNCTION get_units(var_id, units, fac) RESULT(ier)
    ! jdemod
    ! fac is part of an unimplemented automated unit conversion utility
    ! will leave it in for now but not really needed
    ! this routine saves units attribute
    USE netcdf, ONLY : nf90_get_att, nf90_noerr,nf90_inquire_attribute
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_id
    CHARACTER(len=*), INTENT(OUT) :: units
    REAL*8, INTENT(OUT) :: fac
    INTEGER :: ier

    character(:),allocatable :: existing_units
    integer :: unitslen

    ier=nf90_inquire_attribute(nc_id,var_id,unitsname,len=unitslen)
    if ( ier .ne. nf90_noerr ) RETURN

    !    Not implemented
    !    CALL convert_mks(units, new_units, fac)

    ier = switch_to_define_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    if (allocated(existing_units)) deallocate(existing_units)
    allocate(character(unitslen)::existing_units)

    WRITE(err_msg,*) 'Error getting attribute :',unitsname
    ier = handle_nf90_error(nf90_get_att(nc_id, var_id,unitsname,existing_units))

    IF ( ier .eq. nf90_noerr ) then 
       units = existing_units
    endif

    deallocate(existing_units)

    RETURN
  END FUNCTION get_units

  FUNCTION add_units(var_id, units, fac) RESULT(ier)
    ! jdemod
    ! fac is part of an unimplemented automated unit conversion utility
    ! will leave it in for now but not really needed
    ! this routine saves units attribute
    USE netcdf, ONLY : nf90_put_att, nf90_noerr
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_id
    CHARACTER(len=*), INTENT(IN) :: units
    REAL*8, INTENT(OUT) :: fac
    INTEGER :: ier
    !    Not implemented
    !    CALL convert_mks(units, new_units, fac)

    ier = switch_to_define_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error adding units attribute',TRIM(units)
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id,unitsname,TRIM(units)))
    IF ( ier .ne. nf90_noerr ) RETURN

    RETURN
  END FUNCTION add_units


  FUNCTION get_long_name(var_id, long_name) RESULT(ier)
    ! jdemod
    ! add long name attribute
    USE netcdf, ONLY : nf90_get_att, nf90_noerr,nf90_inquire_attribute
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_id
    CHARACTER(len=*), INTENT(OUT) :: long_name
    INTEGER :: ier
    !    CALL convert_mks(units, new_units, fac)

    character(:),allocatable :: existing_long_name
    integer :: namelen

    ier=nf90_inquire_attribute(nc_id,var_id,longname,len=namelen)
    if ( ier .ne. nf90_noerr ) RETURN

    ier = switch_to_define_mode()
    if ( ier .ne. nf90_noerr ) RETURN


    if (allocated(existing_long_name)) deallocate(existing_long_name)
    allocate(character(namelen)::existing_long_name)

    WRITE(err_msg,*) 'Error getting attribute:',longname
    ier = handle_nf90_error(nf90_get_att(nc_id, var_id,longname,existing_long_name))

    ! Note - if long_name is too small then fortran truncates the result
    IF ( ier .eq. nf90_noerr ) then 
       long_name = existing_long_name
    endif
  
    deallocate(existing_long_name)

    RETURN
  END FUNCTION get_long_name


  FUNCTION add_long_name(var_id, long_name) RESULT(ier)
    ! jdemod
    ! add long name attribute
    USE netcdf, ONLY : nf90_put_att, nf90_noerr
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_id
    CHARACTER(len=*), INTENT(IN) :: long_name
    INTEGER :: ier
    !    CALL convert_mks(units, new_units, fac)

    ier = switch_to_define_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error adding long_name attribute',TRIM(long_name)
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id,longname,TRIM(long_name)))

    RETURN
  END FUNCTION add_long_name


  FUNCTION switch_to_define_mode() RESULT(ier)
    ! switching to define mode will not return an error
    USE netcdf, ONLY : nf90_redef, nf90_noerr
    IMPLICIT NONE
    INTEGER :: ier
    ! switch database to define mode .. force to no error
    WRITE(err_msg,*) 'error switching to define mode'
    ier = nf90_redef(nc_id)
    ier = nf90_noerr

    RETURN
  END FUNCTION switch_to_define_mode

  FUNCTION switch_to_data_mode() RESULT(ier)
    ! switching to data mode will not return an error
    USE netcdf, ONLY : nf90_enddef, nf90_noerr
    IMPLICIT NONE
    INTEGER :: ier

    ! switch database to data mode ... forced no error
    WRITE(err_msg,*) 'error switching to data mode'
    ier = nf90_enddef(nc_id)
    ier = nf90_noerr

    RETURN
  END FUNCTION switch_to_data_mode





  !----------------------------------------------------------------
  ! Character Data routines
  ! Only writes a character string

  FUNCTION write_char_data(var_name,var,long_name,units) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_CHAR, nf90_put_var, nf90_max_name
    use error_handling

    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name,  var
    character(len=*), intent(in),optional :: long_name, units
    !integer, intent(in) :: ndims

    !CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    !INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    !integer , optional, INTENT(IN) :: vs
    !integer , DIMENSION(:), optional, INTENT(IN) :: v1d
    !integer , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    !integer , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    !integer , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    !INTEGER, DIMENSION(ndims) :: dim_ids
    !LOGICAL, DIMENSION(ndims) :: create_dim
    logical :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name
    integer :: dim_max,dim_val
    character*2 :: rank_str,type_str

    rank_str='1D'
    type_str='CH'

    ! create a unique dim_name by concatentating the _LEN with the var_name

    dim_name = trim(var_name)//'_LEN'
    dim_max = len(var)
    dim_val = len_trim(var)

    ! Add or check dimension values
    ! Check for non-scalar values

    ! inquire on variable NAME 
    ier = nf90_inq_varid(nc_id, var_name, var_id)


    IF (ier .ne. nf90_noerr) THEN
       ! If the variable is not found 

       ier = nf90_inq_dimid(nc_id, dim_name, dim_id)

       IF (ier .ne. nf90_noerr) THEN

          ier = switch_to_define_mode()
          IF (ier .ne. nf90_noerr) RETURN

          WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
          ! create dimension in netcdf dataset if they are missing
          ! use actual length of string for now ... may want to change to max?
          ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, dim_val, dim_id))
          IF (ier .ne. nf90_noerr) RETURN
       ELSE
          ier = handle_nf90_error(check_existing_dimension(var_name,dim_id, dim_name, dim_val))
          IF (ier .ne. nf90_noerr) RETURN
       ENDIF


       ! switch to define mode and add the variable definition
       ier = switch_to_define_mode()
       IF (ier .ne. nf90_noerr) RETURN

       WRITE(err_msg,*) 'Problem creating variable ', trim(var_name),' with rank '//trim(rank_str)//' and type '//trim(type_str)
       !WRITE(err_msg,*) 'Problem creating variable ', var_name

       ! define character array
       ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_CHAR, dim_id, var_id))
       IF (ier .ne. nf90_noerr) RETURN

    ELSE
       ! variable exists already ... check that it matches ... if it does then data will be over-written
       ! err_msg set by check_existing_variable

       ier = handle_nf90_error(check_existing_variable(var_name, NF90_CHAR, var_id, (/dim_name/) ))
       IF (ier .ne. nf90_noerr) RETURN

    ENDIF

    fac = 1.

    ! add units and long name attribute if present
    if (present(units)) then 
       ier = add_units(var_id,units,fac)
       if ( ier .ne. nf90_noerr ) RETURN
    endif

    if (present(long_name)) then 
       ier = add_long_name(var_id,long_name)
       if ( ier .ne. nf90_noerr ) RETURN
    endif

    !ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN

    ! switch to data mode to store the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Write data to file
    ier = handle_nf90_error(nf90_put_var(nc_id,var_id,var(1:dim_val)))
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION write_char_data


  FUNCTION write_ch_data(var_name,ndims,long_name,units,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! NOTE: This approach does not work ... need a different implementation since character arrays
    !       can not be passed to a higher dimensional argument using an interface. Probably an elegane solution but
    !       no time to find it at the moment
    !       
    ! for writing character arrays - The *len becomes one of the dimensions and must be explicitly assigned
    ! A character variable is treated as a 1D array in this routine. Only a single character is a scalar
    ! 
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_CHAR, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name
    integer, intent(in) :: ndims

    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units

    CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    INTEGER , optional, INTENT(IN) :: v0d
    INTEGER , DIMENSION(:), optional, INTENT(IN) :: v1d
    INTEGER , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    INTEGER , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    INTEGER , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(ndims) :: dim_ids
    !LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str

    type_val = NF90_CHAR
    type_str = 'CH'

    !include 'nc_write_data_generic.f90.inc'

    return

  END FUNCTION write_ch_data





  !----------------------------------------------------------------
  ! INTEGER Data routines
  !
  ! Version 2 - uses rank signature 
  !
  !

  function write_sc_i4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(IN) :: vars
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Dim names and nd are required for everything except scalars
    ! These are still declared since the included code will contain them
    ! However, they are not used and are not optional arguments
    !
    ! Local variables
    !
    !CHARACTER(len=10), DIMENSION(1) :: dim_names
    !INTEGER, DIMENSION(1):: nd

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    !INTEGER, DIMENSION(1) :: dim_ids
    !LOGICAL, DIMENSION(1) :: create_dim
    !CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    integer :: opt
    character*2 :: rank_str,type_str

    ! Set the code option based on the array type
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = 'SC'

    include 'nc_write_data2_generic_scalar.f90.inc'

    return
  end function write_sc_i4_data

  function write_1d_i4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(IN) :: vars(:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '1D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_1d_i4_data

  function write_2d_i4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(IN) :: vars(:,:)
    !
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '2D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_2d_i4_data


  function write_3d_i4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(IN) :: vars(:,:,:)
    !
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '3D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_3d_i4_data


  function write_4d_i4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(IN) :: vars(:,:,:,:)
    !
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '4D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_4d_i4_data



  !----------------------------------------------------------------
  ! REAL*4/FLOAT Data routines
  !
  ! Version 2 - uses rank signature 
  !
  !

  function write_sc_r4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*4 , INTENT(IN) :: vars
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Dim names and nd are required for everything except scalars
    ! These are still declared since the included code will contain them
    ! However, they are not used and are not optional arguments
    !
    ! Local variables
    !
    !CHARACTER(len=10), DIMENSION(1):: dim_names
    !INTEGER, DIMENSION(1) :: nd

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    !INTEGER, DIMENSION(1) :: dim_ids
    !LOGICAL, DIMENSION(1) :: create_dim
    !CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = 'SC'

    include 'nc_write_data2_generic_scalar.f90.inc'

    return
  end function write_sc_r4_data

  function write_1d_r4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*4 , INTENT(IN) :: vars(:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '1D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_1d_r4_data

  function write_2d_r4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*4 , INTENT(IN) :: vars(:,:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '2D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_2d_r4_data

  function write_3d_r4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*4 , INTENT(IN) :: vars(:,:,:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '3D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_3d_r4_data

  function write_4d_r4_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*4 , INTENT(IN) :: vars(:,:,:,:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '4D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_4d_r4_data



  !----------------------------------------------------------------
  ! REAL*8 Data routine
  !
  ! Version 2 - resolve by rank
  !



  function write_sc_r8_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*8 , INTENT(IN) :: vars
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Dim names and nd are required for everything except scalars
    ! These are still declared since the included code will contain them
    ! However, they are not used and are not optional arguments
    !
    ! Local variables
    !
    !CHARACTER(len=10), DIMENSION(1) :: dim_names
    !INTEGER, DIMENSION(1) :: nd

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    !INTEGER, DIMENSION(1) :: dim_ids
    !LOGICAL, DIMENSION(1) :: create_dim
    !CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = 'SC'

    include 'nc_write_data2_generic_scalar.f90.inc'

    return
  end function write_sc_r8_data

  function write_1d_r8_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*8 , INTENT(IN) :: vars(:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '1D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_1d_r8_data

  function write_2d_r8_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*8 , INTENT(IN) :: vars(:,:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '2D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_2d_r8_data

  function write_3d_r8_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*8 , INTENT(IN) :: vars(:,:,:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '3D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_3d_r8_data


  function write_4d_r8_data(var_name,vars,dim_names,nd,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL*8 , INTENT(IN) :: vars(:,:,:,:)
    CHARACTER(len=*), DIMENSION(size(shape(vars))) :: dim_names
    INTEGER, INTENT(IN), DIMENSION(size(shape(vars))) :: nd
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(IN), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(size(shape(vars))) :: dim_ids
    !LOGICAL, DIMENSION(size(shape(vars))) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '4D'

    include 'nc_write_data2_generic.f90.inc'

    return
  end function write_4d_r8_data



  !-------------------------------------------------------------------- 
  ! READ routines


  !----------------------------------------------------------------
  ! Character Data routines
  ! Only reads a character string

  FUNCTION read_char_data(var_name,var,long_name,units) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_inquire_dimension, &
         & nf90_inq_varid, NF90_CHAR, nf90_get_var, nf90_max_name, nf90_max_var_dims,nf90_inquire_variable
    use error_handling

    IMPLICIT NONE

    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name
    CHARACTER(len=*), INTENT(OUT) ::     var
    ! if specified then they are loaded if available
    CHARACTER(len=*), INTENT(OUT), optional :: long_name, units

    !integer, intent(in) :: ndims
    !CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    !INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    !integer , optional, INTENT(IN) :: vs
    !integer , DIMENSION(:), optional, INTENT(IN) :: v1d
    !integer , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    !integer , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    !integer , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    !INTEGER, DIMENSION(ndims) :: dim_ids
    !LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    CHARACTER(len=NF90_MAX_NAME) :: existing_name, existing_dim_name
    INTEGER ::  existing_type, existing_ndims, ind
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: existing_dimids
    integer :: existing_dim_len
    integer :: dim_max


    ! initialize to blank
    var = ''

    dim_name = trim(var_name)//'_LEN'

    dim_max = len(var)
    !dim_val = len_trim(var)

    ! Add or check dimension values
    ! Check for non-scalar values

    ! inquire on variable NAME 
    ! Exit if it is not found

    ier = nf90_inq_varid(nc_id, var_name, var_id)

    if (ier.ne.nf90_noerr) then 
       if (verbose) call errmsg('NC_UTILS_GENERIC: READ_CHAR_DATA: ',trim(var_name)//' NOT IN DATASET')
       ! Because of the non-stored Zero assumption present in some cases - set all returned variables to zero 
       ! - make sure ier is set so calling routine can deal with possible errors
       var = ''
       return
    endif

    ! get the dimension (i.e. length) of the string variable

    ier=handle_nf90_error(nf90_inquire_variable(nc_id,var_id,existing_name,existing_type,existing_ndims,existing_dimids))

    if (existing_type.ne.NF90_CHAR.or.existing_ndims.ne.1) then 
       if (verbose) call errmsg('NC_UTILS_GENERIC: READ TEXT:'//trim(var_name)//' HAS INCORRECT TYPE OR DIMENSIONS',existing_ndims)
       ier = -1
       return
    endif

    ! check to make sure that the max length of the input variable is longer than the string to be read

    ier = handle_nf90_error(nf90_inquire_dimension(nc_id, existing_dimids(1), existing_dim_name, existing_dim_len))
    IF (ier .ne. nf90_noerr) RETURN

    if (dim_max.lt.existing_dim_len) then
       if (verbose) call errmsg('NC_UTILS_GENERIC: INPUT STRING TOO SMALL',existing_dim_len)
       ier = -1
       return
    endif

    ! switch to data mode to load the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error reading ', var_name, ' from file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Read data from file
    ier = handle_nf90_error(nf90_get_var(nc_id,var_id,var(1:existing_dim_len)))

    if ( ier .ne. nf90_noerr ) RETURN

    ! if attributes are requested then load them

    if (present(long_name)) then 
       ier = handle_nf90_error(get_long_name(var_id,long_name))
       if (ier.ne.nf90_noerr) then
          long_name = ''
       endif
    endif

    fac = 1.0

    if (present(units)) then 
       ier = handle_nf90_error(get_units(var_id,units,fac))
       if (ier.ne.nf90_noerr) then
          units = ''
       endif
    endif

  END FUNCTION read_char_data


  FUNCTION read_ch_data(var_name,ndims,long_name,units,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! NOTE: This approach does not work ... need a different implementation since character arrays
    !       can not be passed to a higher dimensional argument using an interface. Probably an elegane solution but
    !       no time to find it at the moment
    !       
    ! import netcdf functions
    ! read character array
    USE netcdf, ONLY : nf90_noerr, nf90_inquire_variable, &
         & NF90_CHAR, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name
    integer, intent(in) :: ndims
    integer :: ier

    character(len=*), optional :: long_name, units

    CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    INTEGER , optional, INTENT(IN) :: v0d
    INTEGER , DIMENSION(:), optional, INTENT(IN) :: v1d
    INTEGER , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    INTEGER , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    INTEGER , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    integer :: var_id,ierr

    !INTEGER :: dim_id, var_id, ier
    !INTEGER, DIMENSION(ndims) :: dim_ids
    !LOGICAL, DIMENSION(ndims) :: create_dim
    !CHARACTER(len=nf90_max_name) :: dim_name

!    integer, dimension(ndims) :: arr_shape

    integer :: type_val
    character*2 :: type_str

    type_val = NF90_CHAR
    type_str = 'CH'

    !include 'nc_read_data_generic.f90.inc'

    return

  END FUNCTION read_ch_data


  !----------------------------------------------------------------
  ! INTEGER Data routines
  !
  ! Version 2 - uses rank signature 
  !

  function read_sc_i4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(OUT) :: vars
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Dim names and nd are required for everything except scalars
    ! These are still declared since the included code will contain them
    ! However, they are not used and are not optional arguments
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(1) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = 'SC'

    include 'nc_read_data2_generic_scalar.f90.inc'

    return
  end function read_sc_i4_data

  function read_1d_i4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(OUT) :: vars(:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    ! opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '1D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_1d_i4_data

  function read_2d_i4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(OUT) :: vars(:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '2D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_2d_i4_data

  function read_3d_i4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(OUT) :: vars(:,:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '3D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_3d_i4_data


  function read_4d_i4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER , INTENT(OUT) :: vars(:,:,:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_INT
    type_str = 'I4'
    rank_str = '4D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_4d_i4_data




  !----------------------------------------------------------------
  ! Real*4 Data routines
  !
  ! Version 2 - uses rank signature 
  !

  function read_sc_r4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*4 , INTENT(OUT) :: vars
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Dim names and nd are required for everything except scalars
    ! These are still declared since the included code will contain them
    ! However, they are not used and are not optional arguments
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(1) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = 'SC'

    include 'nc_read_data2_generic_scalar.f90.inc'

    return
  end function read_sc_r4_data


  function read_1d_r4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*4 , INTENT(out) :: vars(:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '1D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_1d_r4_data



  function read_2d_r4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*4 , INTENT(out) :: vars(:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '2D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_2d_r4_data


  function read_3d_r4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*4 , INTENT(out) :: vars(:,:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '3D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_3d_r4_data


  function read_4d_r4_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*4 , INTENT(out) :: vars(:,:,:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_FLOAT
    type_str = 'R4'
    rank_str = '1D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_4d_r4_data



  !----------------------------------------------------------------
  ! Real*8 Data routines
  !
  ! Version 2 - uses rank signature 
  !

  function read_sc_r8_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*8 , INTENT(OUT) :: vars
    CHARACTER(len=*), INTENT(OUT), optional ::  long_name, units
    !
    ! Dim names and nd are required for everything except scalars
    ! These are still declared since the included code will contain them
    ! However, they are not used and are not optional arguments
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(1) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = 'SC'

    include 'nc_read_data2_generic_scalar.f90.inc'

    return
  end function read_sc_r8_data


  function read_1d_r8_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*8 , INTENT(out) :: vars(:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '1D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_1d_r8_data



  function read_2d_r8_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*8 , INTENT(out) :: vars(:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '2D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_2d_r8_data



  function read_3d_r8_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*8 , INTENT(out) :: vars(:,:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '3D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_3d_r8_data



  function read_4d_r8_data(var_name,vars,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(IN) :: var_name
    real*8 , INTENT(out) :: vars(:,:,:,:)
    !
    ! Optional Arguments
    !
    CHARACTER(len=*), INTENT(out), optional ::  long_name, units
    !
    ! Local variables
    !

    REAL*8 :: fac

    INTEGER :: var_id, ier
!    integer, dimension(size(shape(vars))) :: arr_shape

    integer :: type_val
    character*2 :: type_str
    integer :: opt
    character*2 :: rank_str

    ! Set the code option based on the array type
    opt = size(shape(vars))
    type_val = NF90_DOUBLE
    type_str = 'R8'
    rank_str = '4D'

    include 'nc_read_data2_generic.f90.inc'

    return
  end function read_4d_r8_data



  !-------------------------------------------------------------------------------------------
  ! Test related code



  SUBROUTINE test(ier, str, known_fail)
    USE netcdf, ONLY : nf90_noerr
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ier
    CHARACTER(len=*), INTENT(IN) :: str
    LOGICAL, OPTIONAL :: known_fail
    LOGICAL :: kf
    CHARACTER(len=*), PARAMETER :: fail='FAIL: ', pass='PASS: '
    IF ( .not. PRESENT(known_fail) ) THEN
       kf = .false.
    ELSE
       kf = known_fail
    ENDIF
    IF ( (ier .ne. nf90_noerr .and. kf) .or. &
         &(ier .eq. nf90_noerr .and. .not. kf) ) THEN
       WRITE(*,*) pass//str
    ELSE
       WRITE(*,*) fail//str
    ENDIF
  END SUBROUTINE test



  function test_nc_utils () result(ier)
    USE netcdf, ONLY : nf90_create, nf90_close, nf90_noerr
    !USE nc_utils, ONLY: nc_id, write_1d_array, write_2d_array, test, verbose, write_3d_array, write_scalar
    IMPLICIT NONE
    integer :: ier

    INTEGER,PARAMETER :: N = 100, N1=10, N2=8, N3=8
    integer :: n_r,n1_r
    LOGICAL :: exist
    REAL*8 :: scalar_r8,scalar_r8_z,scalar_r8_r
    REAL*8, DIMENSION(N) :: arr1d_r8, zeroarr1d_r8, arr1d_r8_r       ,zeroarr1d_r8_r
    REAL*8, DIMENSION(N1,N2) :: arr2d_r8, zeroarr2d_r8, arr2d_r8_r   ,zeroarr2d_r8_r
    REAL*8, DIMENSION(N1,N2,N3) :: arr3d_r8, zeroarr3d_r8, arr3d_r8_r,zeroarr3d_r8_r

    REAL*4 :: scalar_r4,scalar_r4_r
    REAL*4, DIMENSION(N) :: arr1d_r4, zeroarr1d_r4, arr1d_r4_r       ,zeroarr1d_r4_r
    REAL*4, DIMENSION(N1,N2) :: arr2d_r4, zeroarr2d_r4, arr2d_r4_r   ,zeroarr2d_r4_r
    REAL*4, DIMENSION(N1,N2,N3) :: arr3d_r4, zeroarr3d_r4, arr3d_r4_r,zeroarr3d_r4_r

    integer :: scalar_i4,scalar_i4_r,scalar_i4_nf_r
    integer, DIMENSION(N) :: arr1d_i4, zeroarr1d_i4, arr1d_i4_r       ,zeroarr1d_i4_r
    integer, DIMENSION(N1,N2) :: arr2d_i4, zeroarr2d_i4, arr2d_i4_r   ,zeroarr2d_i4_r
    integer, DIMENSION(N1,N2,N3) :: arr3d_i4, zeroarr3d_i4, arr3d_i4_r,zeroarr3d_i4_r

    character*40 :: test_string,test_string1,test_string2,test_string_1_r,long_name_r,units_r
    character*30 :: test_string_2_r
    character*6 :: test_string_3_r
    character*10 :: char_array(N1)

    CHARACTER(len=*), PARAMETER :: fn = 'test_nc_utils.nc'

    integer :: i,j,k

    test_string = fn


    char_array(1) = 'A'
    char_array(2) = 'BC'
    char_array(3) = 'DEFGHI'
    char_array(4) = 'jkl'
    char_array(5) = 'mnop'
    char_array(6) = 'qrstuvw'
    char_array(7) = 'Xyz'
    char_array(8) = '012345'
    char_array(9) = '6'
    char_array(10) = '789'

    DO i=1,N
       arr1d_r8(i) = dble(i)
       arr1d_r4(i) = real(i)+1.0
       arr1d_i4(i) = i - 1
    ENDDO
    DO i=1,N1
       DO j=1,N2
          arr2d_r8(i,j) = (i**2+j**2)**0.5
          DO k=1,N3
             arr3d_r8(i,j,k) = (i**2+j**2+k**2)**0.5
          ENDDO
       ENDDO
    ENDDO

    arr2d_r4 = arr2d_r8 + 1.0
    arr2d_i4 = arr2d_r8 - 1

    arr3d_r4 = arr3d_r8 + 1.0
    arr3d_i4 = arr3d_r8 - 1


    zeroarr1d_r8 = 0.0
    zeroarr2d_r8 = 0.0
    zeroarr3d_r8 = 0.0

    zeroarr1d_r4 = 0.0
    zeroarr2d_r4 = 0.0
    zeroarr3d_r4 = 0.0

    zeroarr1d_i4 = 0.0
    zeroarr2d_i4 = 0.0
    zeroarr3d_i4 = 0.0

    scalar_r8 = 10.0
    scalar_r8_z = 0.0
    scalar_r4 = 4.0
    scalar_i4 = 135


    ! Open nc file 
    ier = open_nc_file(fn,NC_WRITE,debug=.false.)

    IF ( ier .ne. nf90_noerr ) THEN
       WRITE(*,*) 'Problem opening file', fn
       STOP
    ENDIF

    ! --------------- character ------------------

    CALL test(write_nc('test_text_string','test_string','This is a text string','T'), 'Write a text string')
    CALL test(write_nc('test_text_string2',test_string,'This is also a text string','T'), 'Write another test string')
    !CALL test(write_nc('test_text_array',size(shape(char_array))+1,'This is a text array','T',&
    !         &(/'char_len','char_dim1'/),(/len(char_array(1)),size(char_array)/),v2d=char_array), 'Write another test string')


    ! -------------- Scalar ---------------------

    ! scalar tests
    CALL test(write_nc('test_scalar_r8',scalar_r8,'This is a scalar','W'), 'Write R8 scalar')
    CALL test(write_nc('zero_scalar_r8',0.0D0,'This is 0','W'), 'Not write R8 scalar 0')
    CALL test(write_nc('zero_scalar_r8',scalar_r8_z,'This is 0','W'), 'Not write R8 scalar 0')
    CALL test(write_nc('test_scalar_r4',scalar_r4,'This is a scalar','W'), 'Write R4 scalar')
    CALL test(write_nc('zero_scalar_r4',0.0,'This is 0','W'), 'Not write scalar R4 0')
    CALL test(write_nc('test_scalar_i4',scalar_i4,'This is a scalar','W'), 'Write I4scalar')
    CALL test(write_nc('zero_scalar_i4',0,'This is 0','W'), 'Not write I4 scalar 0')

    ! Test to see if a variable can be stored with the same name as a dimension

    CALL test(write_nc('dim_arr',N,'This is a scalar','W'), 'Write dim_arr variable')
    

    ! ----------- 1D ---------------------

    ! 1d tests - r8
    CALL test(write_nc('test_1d_arr_r8',arr1d_r8,['dim_arr'],[N],'This is the long name','-'),'Write R8 1d array')
    CALL test(write_nc('zero_1d_arr_r8',zeroarr1d_r8,['dim_arr'],[N],'This array is all 0s','-'),'Not write R8 1d array of zeros')
    CALL test(write_nc('test_1d_arr_r8',arr1d_r8*2.0d0,['dim_arr'],[N],'This is the long name','-'),'Update R8 1d array with value*2')
    ! update arr content for later read comparison
    arr1d_r8 = arr1d_r8 * 2.0
    CALL test(write_nc('test_1d_arr_r8',arr1d_r8*10.0d0,['dim_arr2'],[N],'arr*10','-'),'Not overwite R8 1d array with wrong dimension name', .true.)
    CALL test(write_nc('test_1d/10_arr_r8',arr1d_r8(1:N/10),['dim_arr'],[N/10],'N/10','-'),'Not write R8 1d array with wrong dimension length', .true.)
    CALL test(write_nc('test_1d/10_arr_r8',arr1d_r8(1:N/10),['dim_arr'],[N],'arr(1:N/10)','-'),'Not write R8 1d array with wrong size', .true.)


    ! 1d tests - r4
    CALL test(write_nc('test_1d_arr_r4',arr1d_r4,['dim_arr'],[N],'This is the long name','-'),'Write R4 1d array')
    CALL test(write_nc('zero_1d_arr_r4',zeroarr1d_r4,['dim_arr'],[N],'This array is all 0s','-'),'Not write R4 1d array of zeros')
    CALL test(write_nc('test_1d_arr_r4',arr1d_r4*2.0,['dim_arr'],[N],'This is the long name','-'),'Update R4 1d array with value*2')
    ! update arr content for later read comparison
    arr1d_r4 = arr1d_r4 * 2.0
    CALL test(write_nc('test_1d_arr_r4',arr1d_r4*10.0,['dim_arr2'],[N],'arr*10','-'),'Not overwite R4 1d array with wrong dimension name', .true.)
    CALL test(write_nc('test_1d/10_arr_r4',arr1d_r4(1:N/10),['dim_arr'],[N/10],'N/10','-'),'Not write R4 1d array with wrong dimension length', .true.)
    CALL test(write_nc('test_1d/10_arr_r4',arr1d_r4(1:N/10),['dim_arr'],[N],'arr(1:N/10)','-'),'Not write R4 1d array with wrong size', .true.)

    ! 1d tests - i4
    CALL test(write_nc('test_1d_arr_i4',arr1d_i4,['dim_arr'],[N],'This is the long name','-'),'Write I4 1d array')
    CALL test(write_nc('zero_1d_arr_i4',zeroarr1d_i4,['dim_arr'],[N],'This array is all 0s','-'),'Not write I4 1d array of zeros')
    CALL test(write_nc('test_1d_arr_i4',arr1d_i4*2,['dim_arr'],[N],'This is the long name','-'),'Update I4 1d array with value*2')
    ! update arr content for later read comparison
    arr1d_i4 = arr1d_i4 * 2.0
    CALL test(write_nc('test_1d_arr_i4',arr1d_i4*10,['dim_arr2'],[N],'arr*10','-'),'Not overwite I4 1d array with wrong dimension name', .true.)
    CALL test(write_nc('test_1d/10_arr_i4',arr1d_i4(1:N/10),['dim_arr'],[N/10],'N/10','-'),'Not write I4 1d array with wrong dimension length', .true.)
    CALL test(write_nc('test_1d/10_arr_i4',arr1d_i4(1:N/10),['dim_arr'],[N],'arr(1:N/10)','-'),'Not write I4 1d array with wrong size', .true.)

    ! -------------------- 2D ----------------------

    ! 2d tests - r8
    CALL test(write_nc('test_2d_arr_r8',arr2d_r8, (/'dim_arr1','dim_arr2'/),(/N1,N2/),'This is the long name','-'),'Write R8 2d array')
    CALL test(write_nc('zero_2d_arr_r8',zeroarr2d_r8,(/'zerodim01','zerodim10'/),(/N1,N2/),'This array is all 0s','-'),'Not write R8 2d array of zeros')
    CALL test(write_nc('test_2d_arr_r8',arr2d_r8*2.0d0, (/'dim_arr1','dim_arr2'/),(/N1,N2/),'This is the long name','-'),'Update R8 2d array with value*2')
    ! update arr content for later read comparison
    arr2d_r8 = arr2d_r8 * 2.0
    CALL test(write_nc('test_2d_arr_r8',arr2d_r8*10.0d0,(/'dim_arr3','dim_arr4'/),(/N1,N2/),'arr*10','-'),'Not overwite R8 2d array with wrong dimension name', .true.)
    CALL test(write_nc('test_2d/10_arr_r8',arr2d_r8,(/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/),'N/10','-'),'Not overwrite R8 2d array with wrong dimension length', .true.)
    CALL test(write_nc('test_2d/10_arr_r8',arr2d_r8(1:N1,1:N2/2), (/'dim_arr1','dim_arr2'/),(/N1,N2/),'arr(1:N/10)','-'),'Not write R8 2d array with wrong size', .true.)

    ! 2d tests - r4
    CALL test(write_nc('test_2d_arr_r4',arr2d_r4, (/'dim_arr1','dim_arr2'/),(/N1,N2/),'This is the long name','-'),'Write R4 2d array')
    CALL test(write_nc('zero_2d_arr_r4',zeroarr2d_r4,(/'zerodim01','zerodim10'/),(/N1,N2/),'This array is all 0s','-'),'Not write R4 2d array of zeros')
    CALL test(write_nc('test_2d_arr_r4',arr2d_r4*2.0, (/'dim_arr1','dim_arr2'/),(/N1,N2/),'This is the long name','-'),'Update R4 2d array with value*2')
    ! update arr content for later read comparison
    arr2d_r4 = arr2d_r4 * 2.0
    CALL test(write_nc('test_2d_arr_r4',arr2d_r4*10.0,(/'dim_arr3','dim_arr4'/),(/N1,N2/),'arr*10','-'),'Not overwite R4 2d array with wrong dimension name', .true.)
    CALL test(write_nc('test_2d/10_arr_r4',arr2d_r4,(/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/),'N/10','-'),'Not overwrite R4 2d array with wrong dimension length', .true.)
    CALL test(write_nc('test_2d/10_arr_r4',arr2d_r4(1:N1,1:N2/2), (/'dim_arr1','dim_arr2'/),(/N1,N2/),'arr(1:N/10)','-'),'Not write R4 2d array with wrong size', .true.)


    ! 2d tests - i4
    CALL test(write_nc('test_2d_arr_i4',arr2d_i4, (/'dim_arr1','dim_arr2'/),(/N1,N2/),'This is the long name','-'),'Write I4 2d array')
    CALL test(write_nc('zero_2d_arr_i4',zeroarr2d_i4,(/'zerodim01','zerodim10'/),(/N1,N2/),'This array is all 0s','-'),'Not write I4 2d array of zeros')
    CALL test(write_nc('test_2d_arr_i4',arr2d_i4*2, (/'dim_arr1','dim_arr2'/),(/N1,N2/),'This is the long name','-'),'Update I4 2d array with value*2')
    ! update arr content for later read comparison
    arr2d_i4 = arr2d_i4 * 2.0
    CALL test(write_nc('test_2d_arr_i4',arr2d_i4*10,(/'dim_arr3','dim_ari4'/),(/N1,N2/),'arr*10','-'),'Not overwite I4 2d array with wrong dimension name', .true.)
    CALL test(write_nc('test_2d/10_arr_i4',arr2d_i4,(/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/),'N/10','-'),'Not overwrite I4 2d array with wrong dimension length', .true.)
    CALL test(write_nc('test_2d/10_arr_i4',arr2d_i4(1:N1,1:N2/2), (/'dim_arr1','dim_arr2'/),(/N1,N2/),'arr(1:N/10)','-'),'Not write I4 2d array with wrong size', .true.)


    ! write a variable with the same name as a dimension ... but different value

    CALL test(write_nc('dim_arr2',N1,'This is a scalar','W'), 'Write DIM_ARR2 variable')
    


    !-------------------- 3D -------------------------

    ! 3d tests - r8
    CALL test(write_nc('zero_3d_arr_r8',zeroarr3d_r8, (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/),'This array is all 0s','kg'),'Not write R8 3d array of zeros')
    CALL test(write_nc('test_3d_arr_r8',arr3d_r8,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'This is the long name','kg'),'Write R8 3d array')
    CALL test(write_nc('test_3d_arr_r8',arr3d_r8*2.0d0,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'This is the long name','kg'),'Update R8 3d array with value*2')
    ! update arr content for later read comparison
    arr3d_r8 = arr3d_r8 * 2.0
    CALL test(write_nc('test_3d_arr_r8',arr3d_r8*10.0d0,(/'dim_arr4','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'arr*10','kg'),'Not overwite R8 3d array with wrong dimension name', .true.)
    CALL test(write_nc('test_3d/10_arr_r8',arr3d_r8,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/2/),'N/10','kg'),'Not write R8 3d array with wrong dimension length', .true.)
    CALL test(write_nc('test_3d/10_arr_r8',arr3d_r8(1:N1,1:N2/2,1:N3),(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'arr(1:N/10)','kg'),'Not write R8 3d array with wrong size', .true.)


    ! 3d tests - r4
    CALL test(write_nc('zero_3d_arr_r4',zeroarr3d_r4, (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/),'This array is all 0s','kg'),'Not write R4 3d array of zeros')
    CALL test(write_nc('test_3d_arr_r4',arr3d_r4,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'This is the long name','kg'),'Write R4 3d array')
    CALL test(write_nc('test_3d_arr_r4',arr3d_r4*2.0,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'This is the long name','kg'),'Update R4 3d array with value*2')
    ! update arr content for later read comparison
    arr3d_r4 = arr3d_r4 * 2.0
    CALL test(write_nc('test_3d_arr_r4',arr3d_r4*10.0,(/'dim_arr4','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'arr*10','kg'),'Not overwite R4 3d array with wrong dimension name', .true.)
    CALL test(write_nc('test_3d/10_arr_r4',arr3d_r4,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/2/),'N/10','kg'),'Not write R4 3d array with wrong dimension length', .true.)
    CALL test(write_nc('test_3d/10_arr_r4',arr3d_r4(1:N1,1:N2/2,1:N3),(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'arr(1:N/10)','kg'),'Not write R4 3d array with wrong size', .true.)

    ! 3d tests - i4
    CALL test(write_nc('zero_3d_arr_i4',zeroarr3d_i4, (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/),'This array is all 0s','kg'),'Not write I4 3d array of zeros')
    CALL test(write_nc('test_3d_arr_i4',arr3d_i4,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'This is the long name','kg'),'Write I4 3d array')
    CALL test(write_nc('test_3d_arr_i4',arr3d_i4*2,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'This is the long name','kg'),'Update I4 3d array with value*2')
    ! update arr content for later read comparison
    arr3d_i4 = arr3d_i4 * 2.0
    CALL test(write_nc('test_3d_arr_i4',arr3d_i4*10,(/'dim_ari4','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'arr*10','kg'),'Not overwite I4 3d array with wrong dimension name', .true.)
    CALL test(write_nc('test_3d/10_arr_i4',arr3d_i4,(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/2/),'N/10','kg'),'Not write I4 3d array with wrong dimension length', .true.)
    CALL test(write_nc('test_3d/10_arr_i4',arr3d_i4(1:N1,1:N2/2,1:N3),(/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),'arr(1:N/10)','kg'),'Not write I4 3d array with wrong size', .true.)






    ier = close_nc_file()

    if (ier.ne.nf90_noerr) call errmsg('Problem closing nc file:',ier)


    ! Lets read back and verify the values written

    ! re-open file
    ier = open_nc_file(fn,NC_READONLY,debug=.false.)

    IF ( ier .ne. nf90_noerr ) THEN
       WRITE(*,*) 'Problem opening file for reading:', fn
       return
    ENDIF


    ! Character test

    CALL test(read_nc('test_text_string',test_string_1_r), 'Read a text string')
    if (trim(test_string_1_r).ne.'test_string') then 
       write(0,*) 'READ TEXT FAIL:',trim(test_string_1_r),':','test_string',':'
    endif

    CALL test(read_nc('test_text_string2',test_string_2_r), 'Read another test string')
    if (trim(test_string_2_r).ne.trim(test_string)) then 
       write(0,*) 'READ TEXT FAIL:',trim(test_string_2_r),':',trim(test_string),':'
    endif

    !CALL test(read_nc('test_text_array',size(shape(char_array)), 'Read another test string')
    !if (trim(test_string_2_r.ne.test_string)) then 
    !   write(0,*) 'READ TEXT FAIL:',trim(test_string_2_r),':',trim(test_string),':'
    !endif

    CALL test(read_nc('test_text_string3',test_string_3_r), 'Test string does not exist',.TRUE.)
    CALL test(read_nc('test_text_string2',test_string_3_r), 'Test string too short',.TRUE.)

    ! scalar tests
    CALL test(read_nc('test_scalar_r8',scalar_r8_r), 'READ R8 scalar')
    if (scalar_r8_r.ne.scalar_r8) then 
       write(0,*) 'ERROR: Scalar R8 not read correctly:',scalar_r8,scalar_r8_r
    endif

    CALL test(read_nc('test_scalar_r4',scalar_r4_r), 'READ R4 scalar')
    if (scalar_r4_r.ne.scalar_r4) then 
       write(0,*) 'ERROR: Scalar R4 not read correctly:',scalar_r4,scalar_r4_r
    endif

    CALL test(read_nc('test_scalar_i4',scalar_i4_r), 'READ I4 scalar')
    if (scalar_i4_r.ne.scalar_i4) then 
       write(0,*) 'ERROR: Scalar i4 not read correctly:',scalar_i4,scalar_i4_r
    endif

    CALL test(read_nc('test_scalar_i4_not_found',scalar_i4_nf_r), 'READ scalar not found',.TRUE.)
    if (scalar_r4_r.ne.scalar_r4) then 
       write(0,*) 'ERROR: R8 not read correctly:',scalar_i4,scalar_i4_r
    endif

    ! read back the variables with dimension names



    CALL test(read_nc('dim_arr',N_R), 'Read dim_arr variable')
    if (N_R.ne.N) then 
       write(0,*) 'ERROR: I4 variable dim_arr not read correctly:',N,N_R
    endif


    CALL test(read_nc('dim_arr2',N1_R), 'Read dim_arr2 variable')
    if (N1_R.ne.N1) then 
       write(0,*) 'ERROR: I4 variable dim_arr2 not read correctly:',N1,N1_R
    endif


    ! 1d tests - r8
    CALL test(read_nc('test_1d_arr_r8',arr1d_r8_r,long_name_r,units_r),'Read 1d r8 array')
    if (any(arr1d_r8_r(:).ne.arr1d_r8(:))) call errmsg('FAIL: 1D R8 variable not read back correctly')
    if (verbose) call errmsg('DISPLAY TEST ATTRIBUTE VALUES FOR VAR = '//'test_1d_arr_r8'//':'//trim(long_name_r)//':'//trim(units_r)//':')

    do i = 1,N1
          if (arr1d_r8_r(i).ne.arr1d_r8(i)) then 
             write(0,'(a,i4,2(1x,g12.5))') '1D R8 DATA READ MISMATCH:',i,arr1d_r8_r(i),arr1d_r8(i)
          endif
    end do


    CALL test(read_nc('zero_1d_arr_r8',zeroarr1d_r8_r),'Assign zeroes to R8 1D array not found',.true.)
    if (any(zeroarr1d_r8_r(:).ne.zeroarr1d_r8(:))) call errmsg('FAIL: 1D R8 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_1d_arr_r8',size(shape(arr1d_1d_r8)),v1d=arr1d_r8),'Not overwite 1d array with wrong dimension name', .true.)
    CALL test(read_nc('test_1d/10_arr_r8',arr1d_r8(1:N/10)),'Not read 1d r8 array with wrong dimension length', .true.)



    ! 1d tests - r4
    CALL test(read_nc('test_1d_arr_r4',arr1d_r4_r),'Read 1d r4 array')
    if (any(arr1d_r4_r(:).ne.arr1d_r4(:))) call errmsg('FAIL: 1D R4 variable not read back correctly')
    do i = 1,N1
          if (arr1d_r4_r(i).ne.arr1d_r4(i)) then 
             write(0,'(a,i4,2(1x,g12.5))') '1D R4 DATA READ MISMATCH:',i,arr1d_r4_r(i),arr1d_r4(i)
          endif
    end do



    CALL test(read_nc('zero_1d_arr_r4',zeroarr1d_r4_r),'Assign zeroes to R4 1D array not found',.true.)
    if (any(zeroarr1d_r4_r.ne.zeroarr1d_r4)) call errmsg('FAIL: 1D R4 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_1d_arr_r4',size(shape(arr1d_r4)),v1d=arr1d_r4),'Not overwite 1d array with wrong dimension name', .true.)
    CALL test(read_nc('test_1d/10_arr_r4',arr1d_r4(1:N/10)),'Not read 1d r4 array with wrong dimension length', .true.)

    ! 1d tests - i4
    CALL test(read_nc('test_1d_arr_i4',arr1d_i4_r),'Read 1d i4 array')
    if (any(arr1d_i4_r.ne.arr1d_i4)) call errmsg('FAIL: 1D I4 variable not read back correctly')
    do i = 1,N1
          if (arr1d_i4_r(i).ne.arr1d_i4(i)) then 
             write(0,'(a,i4,2(1x,g12.5))') '1D I4 DATA READ MISMATCH:',i,arr1d_i4_r(i),arr1d_i4(i)
          endif
    end do



    CALL test(read_nc('zero_1d_arr_i4',zeroarr1d_i4_r),'Assign zeroes to I4 1D array not found',.true.)
    if (any(zeroarr1d_i4_r.ne.zeroarr1d_i4)) call errmsg('FAIL: 1D I4 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_1d_arr_i4',size(shape(arr1d_i4)),v1d=arr1d_i4),'Not overwite 1d array with wrong dimension name', .true.)
    CALL test(read_nc('test_1d/10_arr_i4',arr1d_i4(1:N/10)),'Not read 1d i4 array with wrong dimension length', .true.)



    ! 2d tests - r8
    CALL test(read_nc('test_2d_arr_r8',arr2d_r8_r),'Read 2d r8 array')
    if (any(arr2d_r8_r(:,:).ne.arr2d_r8(:,:))) call errmsg('FAIL: 2D R8 variable not read back correctly')
    do i = 1,N1
       do j = 1,N2
          if (arr2d_r8_r(i,j).ne.arr2d_r8(i,j)) then 
             write(0,'(a,2i4,2(1x,g12.5))') '2D R8 DATA READ MISMATCH:',i,j,arr2d_r8_r(i,j),arr2d_r8(i,j)
          endif
       end do
    end do



    CALL test(read_nc('zero_2d_arr_r8',zeroarr2d_r8_r),'Assign zeroes to R8 2D array not found',.true.)
    if (any(zeroarr2d_r8_r(:,:).ne.zeroarr2d_r8(:,:))) call errmsg('FAIL: 2D R8 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_2d_arr_r8',size(shape(arr2d_r8_r)),v2d=arr2d_r8),'Not overwite 2d array with wrong dimension name', .true.)
    CALL test(read_nc('test_2d/10_arr_r8',arr2d_r8_r(1:N1/2,1:N2)),'Not read array - test_2d/10_arr_r8 - not in database', .true.)

    ! 2d tests - r4
    CALL test(read_nc('test_2d_arr_r4',arr2d_r4_r),'Read 2d r4 array')
    if (any(arr2d_r4_r(:,:).ne.arr2d_r4(:,:))) call errmsg('FAIL: 2D R4 variable not read back correctly')
    do i = 1,N1
       do j = 1,N2
          if (arr2d_r4_r(i,j).ne.arr2d_r4(i,j)) then 
             write(0,'(a,2i4,2(1x,g12.5))') '2D R4 DATA READ MISMATCH:',i,j,arr2d_r4_r(i,j),arr2d_r4(i,j)
          endif
       end do
    end do

    CALL test(read_nc('zero_2d_arr_r4',zeroarr2d_r4_r),'Assign zeroes to R4 2D array not found',.true.)
    if (any(zeroarr2d_r4_r(:,:).ne.zeroarr2d_r4(:,:))) call errmsg('FAIL: 2D R4 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_2d_arr_r4',size(shape(arr2d_r4_r)),v2d=arr2d_r4),'Not overwite 2d array with wrong dimension name', .true.)
    CALL test(read_nc('test_2d_arr_r4',arr2d_r4_r(1:N1,1:N2/2)),'Not read 2d r4 array with wrong dimension length', .true.)


    ! 2d tests - i4
    CALL test(read_nc('test_2d_arr_i4',arr2d_i4_r),'Read 2d i4 array')
    if (any(arr2d_i4_r(:,:).ne.arr2d_i4(:,:))) call errmsg('FAIL: 2D I4 variable not read back correctly')
    do i = 1,N1
       do j = 1,N2
          if (arr2d_i4_r(i,j).ne.arr2d_i4(i,j)) then 
             write(0,'(a,2i4,2(1x,g12.5))') '2D I4 DATA READ MISMATCH:',i,j,arr2d_i4_r(i,j),arr2d_i4(i,j)
          endif
       end do
    end do



    CALL test(read_nc('zero_2d_arr_i4',zeroarr2d_i4_r),'Assign zeroes to I4 2D array not found',.true.)
    if (any(zeroarr2d_i4_r(:,:).ne.zeroarr2d_i4(:,:))) call errmsg('FAIL: 2D I4 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_2d_arr_i4',size(shape(arr2d_i4_r)),v2d=arr2d_i4),'Not overwite 2d array with wrong dimension name', .true.)
    CALL test(read_nc('test_2d_arr_i4',arr2d_i4_r(1:N1/2,1:N2/2)),'Not read 2d i4 array with wrong dimension length', .true.)



    ! 3d tests - r8
    CALL test(read_nc('test_3d_arr_r8',arr3d_r8_r),'Read 3d r8 array')
    if (any(arr3d_r8_r(:,:,:).ne.arr3d_r8(:,:,:))) call errmsg('FAIL: 3D R8 variable not read back correctly')
    do i = 1,N1
       do j = 1,N2
          do k = 1,N3
             if (arr3d_r8_r(i,j,k).ne.arr3d_r8(i,j,k)) then 
                write(0,'(a,3i4,2(1x,g12.5))') '3D R8 DATA READ MISMATCH:',i,j,k,arr3d_r8_r(i,j,k),arr3d_r8(i,j,k)
             endif
          end do
       end do
    end do



    CALL test(read_nc('zero_3d_arr_r8',zeroarr3d_r8_r),'Assign zeroes to R8 3D array not found',.true.)
    if (any(zeroarr3d_r8_r(:,:,:).ne.zeroarr3d_r8(:,:,:))) call errmsg('FAIL: 3D R8 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_3d_arr_r8',size(shape(arr3d_r8_r)),v3d=arr3d_r8),'Not overwite 3d array with wrong dimension name', .true.)
    CALL test(read_nc('test_3d_arr_r8',arr3d_r8_r(1:N1,1:N2,1:N3/2)),'Not read 3d r8 array with wrong dimension length', .true.)

    ! 3d tests - r4
    CALL test(read_nc('test_3d_arr_r4',arr3d_r4_r),'Read 3d r4 array')
    if (any(arr3d_r4_r(:,:,:).ne.arr3d_r4(:,:,:))) call errmsg('FAIL: 3D R4 variable not read back correctly')

    CALL test(read_nc('zero_3d_arr_r4',zeroarr3d_r4_r),'Assign zeroes to R4 3D array not found',.true.)
    if (any(zeroarr3d_r4_r(:,:,:).ne.zeroarr3d_r4(:,:,:))) call errmsg('FAIL: 3D R4 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_3d_arr_r4',size(shape(arr3d_r4_r)),v3d=arr3d_r4),'Not overwite 3d array with wrong dimension name', .true.)
    CALL test(read_nc('test_3_arr_r4',arr3d_r4_r(1:N1/2,1:N2,1:N3)),'Not read 3d r4 array with wrong dimension length', .true.)


    ! 3d tests - i4
    CALL test(read_nc('test_3d_arr_i4',arr3d_i4_r),'Read 3d i4 array')
    if (any(arr3d_i4_r(:,:,:).ne.arr3d_i4(:,:,:))) call errmsg('FAIL: 3D I4 variable not read back correctly')

    CALL test(read_nc('zero_3d_arr_i4',zeroarr3d_i4_r),'Assign zeroes to I4 3D array not found',.true.)
    if (any(zeroarr3d_i4_r(:,:,:).ne.zeroarr3d_i4(:,:,:))) call errmsg('FAIL: 3D I4 Zero variable not read back correctly')
    
    !CALL test(read_nc('test_3d_arr_i4',size(shape(arr3d_i4_r)),v3d=arr3d_i4),'Not overwite 3d array with wrong dimension name', .true.)
    CALL test(read_nc('test_3d_arr_i4',arr3d_i4_r(1:N1,1:N2/2,1:N3)),'Not read 3d i4 array with wrong dimension length', .true.)


    ier = close_nc_file()


  END function test_nc_utils


END MODULE nc_utils_generic

 
