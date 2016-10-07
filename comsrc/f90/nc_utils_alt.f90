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
! Instead the nc_utils module provides the functions:
!  write_scalar, write_1d_array, write_2d_array, write_3d_array
! which can be used to write an array to a netcdf file in a single function call, including
! the indication of dimensions and the long_name and units attributes.
!
! TESTS
! 
! The unit tests for this module are provided at the end of file and can be compiled with
! 
! $FC -I$NETCDF_DIR/include/ -L$NETCDF_DIR/lib -lnetcdff -lnetcdf -DTEST_NC_UTILS -r8 -Wl,-rpath,$NETCDF_DIR/lib nc_utils.F90
! 
! which produces the executable (probably a.out), which can be run to produce the output
! 
! > ./a.out
! PASS: Write scalar
! PASS: Not write scalar 0
! PASS: Write 1d array
! PASS: Not write 1d array of zeros
! PASS: Update 1d array with value*2
! PASS: Not write 1d array with wrong dimension
! PASS: Not write 1d array with wrong dimension length
! PASS: Not write 1d array with wrong size
! PASS: Not write 2d array of zeros
! PASS: Write 2d array
! PASS: Not write 2d array with wrong dimension
! PASS: Not write 2d array with wrong dimension length
! PASS: Not write 2d array with wrong size
! PASS: Not write 2d array to 1d existing variable
! PASS: Not write 3d array of zeros
! PASS: Write 3d array
! PASS: Not write 3d array with wrong dimension
! PASS: Not write 3d array with wrong dimension length
! PASS: Not write 3d array with wrong size
!
! SHORT EXAMPLE
!
! USE netcdf, ONLY : nf90_open, nf90_noerr, nf90_close
! USE nc_utils, ONLY: nc_id, write_1d_array
! The netcdf file is opened by the user
! ier = nf90_open('test.nc',0,ncid)
! The netcdf file (or group) ID is (here explicitly) transferred to the nc_utils module
! nc_id = ncid
! The desired array is written
! ier = write_1d_array( var1, 'var1', 'This is var1', 'meter', (/ 'dim1' /), size(var1) )
! IF ( ier .ne. f90_noerr ) WRITE (*,*) 'There was a problem with var1'
! The netcdf file is closed
! ier = nf90_close(ncid)
!
! PUBLIC attributes
!
! nc_id - The integer variable holding the netcdf ID of the file or group to which to write
!           the user must set nc_id after opening the netcdf file or group
! verbose - A logical variable indicating whether to print error messages to STDOUT
! write_scalar - A function for writing a real*8 scalar to the netcdf file
! write_1d_array - A function for writing a 1d real*8 array to the netcdf file
! write_2d_array - A function for writing a 2d real*8 array to the netcdf file
! write_3d_array - A function for writing a 3d real*8 array to the netcdf file
!
! jdemod
!
! Consider removing the "scaling" factor support since it wasn't complete and doesn't really make much 
! sense. The codes should know what the units are and assign them so scaling the results
! doesn't seem very useful. In addition, the automatic scaling routine for units would need
! to be extremely robust to sense what the specific units were and then automatically 
! convert them to MKS units. 
!
!


MODULE nc_utils

  use error_handling

  implicit none

  PRIVATE
  PUBLIC ::  write_nc,open_nc_file,test_nc_utils,close_nc_file
  !    public ::  open_nc_file,close_nc_file

  CHARACTER(len=1024) :: err_msg ! Holder of error messages
  INTEGER :: nc_id
  LOGICAL :: verbose = .true.

  ! will not store NEW variables with a zero value to the database
  ! Note that when reading data ... IF zero_check is true then all variables that 
  ! are not found will be assigned a zero value. Probably do the same without the zero_check flag
  ! but may issue additional error messages
  LOGICAL :: zero_check = .true.


  interface write_nc
     write_char_data, write_i4_data, write_r4_data, write_r8_data
  end interface write_nc


  interface read_nc
     read_char_data, read_i4_data, read_r4_data, read_r8_data
  end interface write_nc

  integer, parameter :: NC_READONLY = 0
  integer, parameter :: NC_WRITE = 1



CONTAINS

  subroutine open_nc_file(filename,ierr,mode,check_zero_flag,debug)
    ! opens the nc file ... sets nc_id and verbose output/debug
    use error_handling
    use netcdf
    implicit none

    character*(*) :: filename
    logical,optional :: debug,check_zero_flag
    integer,optional :: mode
    integer :: ierr

    integer :: mode_val,create_mode

    verbose = .TRUE.
    zero_check = .TRUE.

    if (present(debug)) verbose = debug
    if (present(check_zero_flag) zero_check = check_zero_flag

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
             call errmsg('OPEN_NC_FILE: ERROR OPENING:'\\trim(filename)\\': FILE FOR WRITING',ierr)
             stop 'OPEN_NC_FILE_WRITE_ERROR'
          endif

       elseif (mode_val.eq.NF90_NOWRITE) then 
          ! database was expected but could not be opened ... issue an error message and quit
          call errmsg('OPEN_NC_FILE: ERROR OPENING:'\\trim(filename)\\': FILE FOR READING',ierr)
       endif

    endif

    ! At this point the netcdf database file has been opened either for reading or writing OR the code exited with an error message
    return
  end subroutine open_nc_file


  subroutine close_nc_file
    use error_handling
    use netcdf, ONLY : nf90_close,nf90_noerr
    implicit none
    integer :: ierr

    ! this routine closes the nc file opened in this module 

    ierr = nf90_close(nc_id)
    if (ierr.ne.nf90_noerr) then 
       call errmsg('NC_UTILS_DIVIMP: CLOSE_NC_FILE: Error closing file = ',ierr)
    endif
    
  end subroutine close_nc_file


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
       IF ( verbose ) call errmsg("NC_UTILS_DIVIMP:"//trim(err_msg),ier)
       !WRITE(*,*) TRIM(err_msg)
    ENDIF
  END FUNCTION handle_nf90_error

  FUNCTION check_existing_variable(var_name, var_type, var_id, &
       & dim_ids,dim_vals) RESULT(ier)
    USE netcdf, ONLY : nf90_inquire_variable, NF90_MAX_NAME, &
         & NF90_MAX_VAR_DIMS, nf90_noerr
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_type, var_id
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), DIMENSION(:),optional :: dim_ids, dim_vals
    CHARACTER(len=NF90_MAX_NAME) :: existing_var_name
    INTEGER ::  existing_type, existing_ndims, ind, ier
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: existing_dimids

    ier = nf90_inquire_variable(nc_id, var_id, existing_var_name,existing_type, existing_ndims, existing_dimids)

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


    if (present(dim_ids)) then 

       ! check if rank of variable matches the expected rank
       IF ( existing_ndims .ne. SIZE(dim_ids) ) THEN
          WRITE(err_msg,*) 'Dimension mismatch for variable ',TRIM(var_name)
          ier = NF90_NOERR + 3
          RETURN
       ENDIF

       ! check if the DIM IDs of the dimensions match  ... NOTE: this is different from checking if the actual sizes match
       DO ind = 1, SIZE(dim_ids)
          IF ( dim_ids(ind) .ne. existing_dimids(ind) ) THEN
             WRITE(err_msg,*) 'Attempt to write variable ',TRIM(var_name), ' with existing dimension different than new dimension'
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

  FUNCTION check_existing_dimension(dim_id, dim_name, dim_len) RESULT(ier)
    USE netcdf, ONLY : NF90_MAX_NAME, nf90_noerr, nf90_inquire_dimension
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim_id, dim_len
    CHARACTER(len=*), INTENT(IN) :: dim_name
    INTEGER :: existing_len, ier
    CHARACTER(len=NF90_MAX_NAME) :: existing_name

    WRITE(err_msg,*) 'Error inquiring dimension ', TRIM(dim_name)

    ier = handle_nf90_error(nf90_inquire_dimension(nc_id, dim_id, existing_name, existing_len))
    IF ( ier .ne. nf90_noerr ) RETURN

    ! check dimension name
    IF (TRIM(dim_name) .ne. TRIM(existing_name) ) THEN
       ier = nf90_noerr + 1
       WRITE(err_msg,*) 'Dim id ',dim_id, ' name mismatch ', TRIM(dim_name), ' ', TRIM(existing_name)
       RETURN
    ENDIF

    ! check value of dimension
    IF (dim_len .ne. existing_len) THEN
       ier = nf90_noerr + 2
       WRITE(err_msg,*) 'Dimension ', TRIM(dim_name), ' is already of size ', existing_len, ' not ', dim_len
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
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id, 'units',TRIM(units)))
    IF ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error adding long_name attribute',TRIM(long_name)
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id, 'long_name',TRIM(long_name)))

    RETURN
  END FUNCTION add_units_long_name

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

  FUNCTION write_char_data(var_name,long_name,units,var) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_CHAR, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units, var
    !integer, intent(in) :: ndims

    !CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    !INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    !integer , optional, INTENT(IN) :: vs
    !integer , DIMENSION(:), optional, INTENT(IN) :: v1d
    !integer , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    !integer , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    !integer , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    !REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(ndims) :: dim_ids
    LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    ! create a unique dim_name by concatentating the _LEN with the var_name

    dim_name = trim(var_name)//'_LEN'
    dim_max = len(var)
    dim_val = len_trim(var)

    ! Add or check dimension values
    ! Check for non-scalar values


    ier = nf90_inq_dimid(nc_id, dim_name, dim_id)
    IF (ier .ne. nf90_noerr) THEN
       create_dim = .true.
    ELSE
       create_dim = .false.
       ier = handle_nf90_error(check_existing_dimension(dim_id, dim_name, dim_val))
       IF (ier .ne. nf90_noerr) RETURN
    ENDIF

    ! inquire on variable NAME 
    ier = nf90_inq_varid(nc_id, var_name, var_id)

    IF (ier .ne. nf90_noerr) THEN
       ! If the variable is not found 

       IF (create_dim) THEN
          ! Create any required dimensions that are not already in place for non-scalars
          ier = switch_to_define_mode()
          IF (ier .ne. nf90_noerr) RETURN

          WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
          ! create dimension in netcdf dataset if they are missing
          ! use actual length of string for now ... may want to change to max?
          ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, dim_val, dim_id))
          IF (ier .ne. nf90_noerr) RETURN

       ENDIF

       ! switch to define mode and add the variable definition
       ier = switch_to_define_mode()
       IF (ier .ne. nf90_noerr) RETURN

       WRITE(err_msg,*) 'Problem creating variable ', var_name

       ! define character array
       ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_CHAR, dim_id, var_id))
       IF (ier .ne. nf90_noerr) RETURN

    ELSE
       ! variable exists already ... check that it matches ... if it does then data will be over-written
       ! err_msg set by check_existing_variable

       ier = handle_nf90_error(check_existing_variable(var_name, NF90_CHAR, var_id, dim_id))
       IF (ier .ne. nf90_noerr) RETURN

    ENDIF

    fac = 1.

    ! add units and long name attribute ... may not have anything significant for a character string
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN

    ! switch to data mode to store the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Write data to file
    ier = handle_nf90_error(nf90_put_var(nc_id,var_id,var))
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION write_char_data


  !----------------------------------------------------------------
  ! INTEGER Data routines
  !

  FUNCTION write_i4_data(var_name,long_name,units,ndims,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_INT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units
    integer, intent(in) :: ndims

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
    LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    ! Check to see which rank of data is to be written and set a variable so we do not need 
    ! to check present all the time. NOTE: Non-generic code that depends on the variable rank
    ! must be inside a block using the correct opt identifier

    ! Set variable option to ndims ... then double check that the expected optional argument is set

    if (present(v0d)) then 
       ! scalar
       opt=0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
    elseif (present(v2d)) then 
       ! 1D array
       opt=2
    elseif (present(v3d)) then 
       ! 1D array
       opt=3
    elseif (present(v4d)) then 
       ! 1D array
       opt=4
    endif


    ! check for matching shape and non-zero content
    ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
    !       it becomes impossible to overwrite existing database values with zeroes

    if (opt.eq.1) then 
       ! check matching shape
       IF ( any(SHAPE(v1d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.2) then 
       ! check matching shape
       IF ( any(SHAPE(v2d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.3) then 
       ! check matching shape
       IF ( any(SHAPE(v3d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.4) then 
       ! check matching shape
       IF ( any(SHAPE(v4d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    endif


    ! Add or check dimension values
    ! Check for non-scalar values
    if (opt .ne.0) then 

       DO d = 1,size(nd)
          dim_name = TRIM(dim_names(d))
          ier = nf90_inq_dimid(nc_id, dim_name, dim_ids(d))
          dim_id = dim_ids(d)
          IF (ier .ne. nf90_noerr) THEN
             create_dim(d) = .true.
          ELSE
             create_dim(d) = .false.
             ier = handle_nf90_error(check_existing_dimension(dim_id, dim_name, nd(d)))
             IF (ier .ne. nf90_noerr) RETURN
          ENDIF
       ENDDO

    endif

    ! inquire on variable NAME 
    ier = nf90_inq_varid(nc_id, var_name, var_id)

    IF (ier .ne. nf90_noerr) THEN
       ! If the variable is not found 

       ! check to make sure all values are not zero (if zero_check active) before creating a new variable


       if (zero_check) then 
          if (opt.eq.0) then 
             ! check not all zero
             IF (v0d .eq. 0) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.1) then 
             ! check not all zero
             IF (.not. any(v1d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.2) then 
             ! check not all zero
             IF (.not. any(v2d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.3) then 
             ! check not all zero
             IF (.not. any(v3d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.4) then 
             ! check not all zero
             IF (.not. any(v4d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          endif
       endif




       IF (opt.gt.0.and.any(create_dim) ) THEN
          ! Create any required dimensions that are not already in place for non-scalars
          ier = switch_to_define_mode()
          IF (ier .ne. nf90_noerr) RETURN

          DO d=1,size(nd)
             IF (create_dim(d)) THEN
                ! create dimension
                dim_name = TRIM(dim_names(d))
                WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
                ! create dimensions in netcdf dataset if they are missing
                ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd(d), dim_ids(d)))
                IF (ier .ne. nf90_noerr) RETURN

             ENDIF
          ENDDO
       ENDIF

       ! switch to define mode and add the variable definition
       ier = switch_to_define_mode()
       IF (ier .ne. nf90_noerr) RETURN

       WRITE(err_msg,*) 'Problem creating variable ', var_name

       if (opt.eq.0) then
          ! define scalar
          ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_INT, var_id))
       else
          ! define array
          ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_INT, dim_ids, var_id))
       endif

       IF (ier .ne. nf90_noerr) RETURN

    ELSE
       ! variable exists already ... check that it matches ... if it does then data will be over-written
       ! err_msg set by check_existing_variable

       ! dim_ids should already be a zero length array for scalars
       ier = handle_nf90_error(check_existing_variable(var_name, NF90_INT, var_id, dim_ids))
       IF (ier .ne. nf90_noerr) RETURN

    ENDIF

    fac = 1.

    ! add units and long name attribute
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN

    ! switch to data mode to store the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Write data to file
    if (opt.eq.0) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v0d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v1d(1:nd(1))))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v2d(1:nd(1),1:nd(2))))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v3d(1:nd(1),1:nd(2),1:nd(3))))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v4d(1:nd(1),1:nd(2),1:nd(3),1:nd(4))))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION write_i4_data




  !----------------------------------------------------------------
  ! REAL*4/FLOAT Data routines
  !

  FUNCTION write_r4_data(var_name,long_name,units,ndims,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_FLOAT, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units
    integer, intent(in) :: ndims

    CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    REAL*4 , optional, INTENT(IN) :: v0d
    REAL*4 , DIMENSION(:), optional, INTENT(IN) :: v1d
    REAL*4 , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    REAL*4 , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    REAL*4 , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(ndims) :: dim_ids
    LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    ! Check to see which rank of data is to be written and set a variable so we do not need 
    ! to check present all the time. NOTE: Non-generic code that depends on the variable rank
    ! must be inside a block using the correct opt identifier

    ! Set variable option to ndims ... then double check that the expected optional argument is set

    if (present(v0d)) then 
       ! scalar
       opt=0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
    elseif (present(v2d)) then 
       ! 1D array
       opt=2
    elseif (present(v3d)) then 
       ! 1D array
       opt=3
    elseif (present(v4d)) then 
       ! 1D array
       opt=4
    endif


    ! check for matching shape and non-zero content
    if (opt.eq.1) then 

       ! check matching shape
       IF ( any(SHAPE(v1d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.2) then 

       ! check matching shape
       IF ( any(SHAPE(v2d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.3) then 

       ! check matching shape
       IF ( any(SHAPE(v3d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.4) then 

       ! check matching shape
       IF ( any(SHAPE(v4d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R4_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    endif


    ! Add or check dimension values
    ! Check for non-scalar values
    if (opt .ne.0) then 

       DO d = 1,size(nd)
          dim_name = TRIM(dim_names(d))
          ier = nf90_inq_dimid(nc_id, dim_name, dim_ids(d))
          dim_id = dim_ids(d)
          IF (ier .ne. nf90_noerr) THEN
             create_dim(d) = .true.
          ELSE
             create_dim(d) = .false.
             ier = handle_nf90_error(check_existing_dimension(dim_id, dim_name, nd(d)))
             IF (ier .ne. nf90_noerr) RETURN
          ENDIF
       ENDDO

    endif

    ! inquire on variable NAME 
    ier = nf90_inq_varid(nc_id, var_name, var_id)

    IF (ier .ne. nf90_noerr) THEN
       ! If the variable is not found 

       if (zero_check) then 
          if (opt.eq.0) then 
             ! check not all zero
             IF (v0d .eq. 0) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.1) then 
             ! check not all zero
             IF (.not. any(v1d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.2) then 
             ! check not all zero
             IF (.not. any(v2d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.3) then 
             ! check not all zero
             IF (.not. any(v3d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.4) then 
             ! check not all zero
             IF (.not. any(v4d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          endif
       endif



       IF (opt.gt.0.and.any(create_dim) ) THEN
          ! Create any required dimensions that are not already in place for non-scalars
          ier = switch_to_define_mode()
          IF (ier .ne. nf90_noerr) RETURN

          DO d=1,size(nd)
             IF (create_dim(d)) THEN
                ! create dimension
                dim_name = TRIM(dim_names(d))
                WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
                ! create dimensions in netcdf dataset if they are missing
                ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd(d), dim_ids(d)))
                IF (ier .ne. nf90_noerr) RETURN

             ENDIF
          ENDDO
       ENDIF

       ! switch to define mode and add the variable definition
       ier = switch_to_define_mode()
       IF (ier .ne. nf90_noerr) RETURN

       WRITE(err_msg,*) 'Problem creating variable ', var_name

       if (opt.eq.0) then
          ! define scalar
          ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_FLOAT, var_id))
       else
          ! define array
          ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_FLOAT, dim_ids, var_id))
       endif

       IF (ier .ne. nf90_noerr) RETURN

    ELSE
       ! variable exists already ... check that it matches ... if it does then data will be over-written
       ! err_msg set by check_existing_variable

       ! dim_ids should already be a zero length array for scalars
       ier = handle_nf90_error(check_existing_variable(var_name, NF90_FLOAT, var_id, dim_ids))
       IF (ier .ne. nf90_noerr) RETURN

    ENDIF

    fac = 1.

    ! add units and long name attribute
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN

    ! switch to data mode to store the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Write data to file
    if (opt.eq.0) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v0d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v1d(1:nd(1))))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v2d(1:nd(1),1:nd(2))))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v3d(1:nd(1),1:nd(2),1:nd(3))))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v4d(1:nd(1),1:nd(2),1:nd(3),1:nd(4))))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION write_r4_data


  !----------------------------------------------------------------
  ! REAL*8 Data routines
  !



  FUNCTION write_r8_data(var_name,long_name,units,ndims,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units
    integer, intent(in) :: ndims

    CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    REAL*8 , optional, INTENT(IN) :: v0d
    REAL*8 , DIMENSION(:), optional, INTENT(IN) :: v1d
    REAL*8 , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    REAL*8 , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    REAL*8 , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(ndims) :: dim_ids
    LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    ! Check to see which rank of data is to be written and set a variable so we do not need 
    ! to check present all the time. NOTE: Non-generic code that depends on the variable rank
    ! must be inside a block using the correct opt identifier

    ! Set variable option to ndims ... then double check that the expected optional argument is set

    if (present(v0d)) then 
       ! scalar
       opt=0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
    elseif (present(v2d)) then 
       ! 1D array
       opt=2
    elseif (present(v3d)) then 
       ! 1D array
       opt=3
    elseif (present(v4d)) then 
       ! 1D array
       opt=4
    endif


    ! check for matching shape - check non-zero for new variables only
    if (opt.eq.1) then 
       ! check matching shape
       IF ( any(SHAPE(v1d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R8_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.2) then 
       ! check matching shape
       IF ( any(SHAPE(v2d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R8_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.3) then 
       ! check matching shape
       IF ( any(SHAPE(v3d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R8_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    elseif (opt.eq.4) then 
       ! check matching shape
       IF ( any(SHAPE(v4d) .ne. nd )) THEN
          IF (verbose)  call errmsg('NC_UTILS_DIVIMP: WRITE_R8_DATA: ',var_name//': Variable and dimension shape mismatch')
          ier = nf90_noerr+10
          RETURN
       ENDIF
    endif


    ! Add or check dimension values
    ! Check for non-scalar values
    if (opt .ne.0) then 

       DO d = 1,size(nd)
          dim_name = TRIM(dim_names(d))
          ier = nf90_inq_dimid(nc_id, dim_name, dim_ids(d))
          dim_id = dim_ids(d)
          IF (ier .ne. nf90_noerr) THEN
             create_dim(d) = .true.
          ELSE
             create_dim(d) = .false.
             ier = handle_nf90_error(check_existing_dimension(dim_id, dim_name, nd(d)))
             IF (ier .ne. nf90_noerr) RETURN
          ENDIF
       ENDDO

    endif

    ! inquire on variable NAME 
    ier = nf90_inq_varid(nc_id, var_name, var_id)

    IF (ier .ne. nf90_noerr) THEN
       ! If the variable is not found 

       if (zero_check) then 
          if (opt.eq.0) then 
             ! check not all zero
             IF (v0d .eq. 0) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.1) then 
             ! check not all zero
             IF (.not. any(v1d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.2) then 
             ! check not all zero
             IF (.not. any(v2d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.3) then 
             ! check not all zero
             IF (.not. any(v3d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          elseif (opt.eq.4) then 
             ! check not all zero
             IF (.not. any(v4d .ne. 0)) THEN
                IF (verbose) call errmsg('NC_UTILS_DIVIMP: WRITE_I4_DATA: ',var_name//' not written to netcdf file because it is all 0')
                ier = nf90_noerr
                RETURN
             ENDIF
          endif
       endif

       IF (opt.gt.0.and.any(create_dim) ) THEN
          ! Create any required dimensions that are not already in place for non-scalars
          ier = switch_to_define_mode()
          IF (ier .ne. nf90_noerr) RETURN

          DO d=1,size(nd)
             IF (create_dim(d)) THEN
                ! create dimension
                dim_name = TRIM(dim_names(d))
                WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
                ! create dimensions in netcdf dataset if they are missing
                ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd(d), dim_ids(d)))
                IF (ier .ne. nf90_noerr) RETURN

             ENDIF
          ENDDO
       ENDIF

       ! switch to define mode and add the variable definition
       ier = switch_to_define_mode()
       IF (ier .ne. nf90_noerr) RETURN

       WRITE(err_msg,*) 'Problem creating variable ', var_name

       if (opt.eq.0) then
          ! define scalar
          ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_DOUBLE, var_id))
       else
          ! define array
          ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_DOUBLE, dim_ids, var_id))
       endif

       IF (ier .ne. nf90_noerr) RETURN

    ELSE
       ! variable exists already ... check that it matches ... if it does then data will be over-written
       ! err_msg set by check_existing_variable

       ! dim_ids should already be a zero length array for scalars
       ier = handle_nf90_error(check_existing_variable(var_name, NF90_DOUBLE, var_id, dim_ids))
       IF (ier .ne. nf90_noerr) RETURN

    ENDIF

    fac = 1.

    ! add units and long name attribute
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN

    ! switch to data mode to store the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Write data to file
    if (opt.eq.0) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v0d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v1d(1:nd(1))))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v2d(1:nd(1),1:nd(2))))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v3d(1:nd(1),1:nd(2),1:nd(3))))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v4d(1:nd(1),1:nd(2),1:nd(3),1:nd(4))))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION write_r8_data



  !-------------------------------------------------------------------- 
  ! READ routines


  !----------------------------------------------------------------
  ! Character Data routines
  ! Only reads a character string

  FUNCTION read_char_data(var_name,var,long_name,units) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
         & nf90_inq_varid, nf90_def_var, NF90_CHAR, nf90_put_var, nf90_max_name, nf90_max_var_dims
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name, var
    ! not used for now ... can be used for additional verification if desired later
    CHARACTER(len=*), INTENT(IN), optional :: long_name, units
    !integer, intent(in) :: ndims

    !CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    !INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    !integer , optional, INTENT(IN) :: vs
    !integer , DIMENSION(:), optional, INTENT(IN) :: v1d
    !integer , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    !integer , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    !integer , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    !REAL*8 :: fac

    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(ndims) :: dim_ids
    LOGICAL, DIMENSION(ndims) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name

    CHARACTER(len=NF90_MAX_NAME) :: existing_var_name, existing_dim_name
    INTEGER ::  existing_type, existing_ndims, ind, ier
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: existing_dimids
    integer :: existing_dim_len
    


    ! initialize to blank
    var = ''

    dim_name = trim(var_name)//'_LEN'
    dim_max = len(var)
    dim_val = len_trim(var)

    ! Add or check dimension values
    ! Check for non-scalar values

    ! inquire on variable NAME 
    ! Exit if it is not found

    ier = nf90_inq_varid(nc_id, var_name, var_id)

    if (ier.ne.nf90_noerr) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA: ',trim(var_name)//' NOT IN DATASET')
       ! Because of the non-stored Zero assumption present in some cases - set all returned variables to zero 
       ! - make sure ier is set so calling routine can deal with possible errors
       var = ''
       return
    endif
    
    ! get the dimension (i.e. length) of the string variable

    ier=handle_nf90_error(nf90_inquire_variable(nc_id,var_id,existing_name,existing_type,exisitng_dims,existing_dimids))
    
    if (existing_type.ne_NF90_CHAR.or.existing_dims.ne.1) then 
       call errmsg('NC_UTILS_DIVIMP: READ TEXT:'//trim(var_name)//' HAS INCORRECT TYPE OR DIMENSIONS',existing_dims)
       ier = -1
       return
    endif
    
    ! check to make sure that the max length of the input variable is longer than the string to be read
    
    ier = handle_nf90_error(nf90_inquire_dimension(nc_id, existing_dimids(1), existing_dim_name, existing_dim_len))
    IF (ier .ne. nf90_noerr) RETURN

    if (dim_max.lt.existing_dim_len) then
       call errmsg('NC_UTILS_DIVIMP: INPUT STRING TOO SMALL',existing_dim_len)
       ier -1
       return
    endif

    ! switch to data mode to load the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error reading ', var_name, ' from file'
    ! add data ... note the lower bounds are assumed to be 1. 

    ! Read data from file
    ier = handle_nf90_error(nf90_get_var(nc_id,var_id,var))
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION read_char_data

  ! Integer

  FUNCTION read_i4_data(var_name,ndims,long_name,units,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_noerr, nf90_inquire_variable, &
         & NF90_INT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name,
    character(len=*), intent(in),optional :: long_name, units
    integer, intent(in) :: ndims

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
    
    integer, dimension(ndims) :: arr_shape
    
    ! This routine reads the specified variable
    ! Almost all the specifiers are optional since they are there ONLY for checking ... checks not implemented at the moment
    
    ! The routine performs the basic check that the shape and size of the array specified
    ! matches the data to be read

    ! If the optional variables are specified they will be checked against the stored quantity as an additional
    ! level of checking. 

    ! Check to see which rank of data is to be read and set a variable so we do not need 
    ! to check present all the time. NOTE: Non-generic code that depends on the variable rank
    ! must be inside a block using the correct opt identifier

    ! Set variable option to ndims ... then double check that the expected optional argument is set

    opt = -1

    if (present(v0d)) then 
       ! scalar
       opt=0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
       arr_shape = shape(v1d)
    elseif (present(v2d)) then 
       ! 1D array
       opt=2
       arr_shape = shape(v2d)
    elseif (present(v3d)) then 
       ! 1D array
       opt=3
       arr_shape = shape(v3d)
    elseif (present(v4d)) then 
       ! 1D array
       opt=4
       arr_shape = shape(v4d)
    endif

    if (opt.eq.-1) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA:','NO OUTPUT VARIABLE SPECIFIED')
       ier = -1
       return
    endif

    ! inquire on variable NAME 
    ! Exit if it is not found
    ier = nf90_inq_varid(nc_id, var_name, var_id)

    if (ier.ne.nf90_noerr) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA: ',trim(var_name)//' NOT IN DATASET')
       ! Because of the non-stored Zero assumption present in some cases - set all returned variables to zero 
       ! - make sure ier is set so calling routine can deal with possible errors
       if (opt.eq.0) then 
          v0d = 0
       elseif (opt.eq.1) then 
          v1d = 0
       elseif (opt.eq.2) then 
          v2d = 0
       elseif (opt.eq.3) then 
          v3d = 0
       elseif (opt.eq.4) then 
          v4d = 0
       endif
       return
    endif
   
    ! check for matching shape 
    ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
    !       it becomes impossible to overwrite existing database values with zeroes
    ! This should have checked the name and shape information for var_id

    ier = check_existing_variable(var_name,NF90_INT,var_id,dim_vals=arr_shape)
    if (ier.ne.nf90_noerr) return

    ! switch to data mode to read the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error reading ', var_name, ' to variable'
    ! get data ... note the lower bounds are assumed to be 1. 

    ! Read data from file
    if (opt.eq.0) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v0d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v1d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v2d))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v3d))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v4d))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION read_i4_data

  ! REAL*4

  FUNCTION read_r4_data(var_name,ndims,long_name,units,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_noerr, nf90_inquire_variable, &
         & NF90_FLOAT, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name,
    character(len=*), intent(in),optional :: long_name, units
    integer, intent(in) :: ndims

    CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    REAL*4 , optional, INTENT(IN) :: v0d
    REAL*4 , DIMENSION(:), optional, INTENT(IN) :: v1d
    REAL*4 , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    REAL*4 , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    REAL*4 , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    integer :: var_id,ierr

    !INTEGER :: dim_id, var_id, ier
    !INTEGER, DIMENSION(ndims) :: dim_ids
    !LOGICAL, DIMENSION(ndims) :: create_dim
    !CHARACTER(len=nf90_max_name) :: dim_name
    
    integer, dimension(ndims) :: arr_shape
    
    ! This routine reads the specified variable
    ! Almost all the specifiers are optional since they are there ONLY for checking ... checks not implemented at the moment
    
    ! The routine performs the basic check that the shape and size of the array specified
    ! matches the data to be read

    ! If the optional variables are specified they will be checked against the stored quantity as an additional
    ! level of checking. 

    ! Check to see which rank of data is to be read and set a variable so we do not need 
    ! to check present all the time. NOTE: Non-generic code that depends on the variable rank
    ! must be inside a block using the correct opt identifier

    ! Set variable option to ndims ... then double check that the expected optional argument is set

    opt = -1

    if (present(v0d)) then 
       ! scalar
       opt=0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
       arr_shape = shape(v1d)
    elseif (present(v2d)) then 
       ! 1D array
       opt=2
       arr_shape = shape(v2d)
    elseif (present(v3d)) then 
       ! 1D array
       opt=3
       arr_shape = shape(v3d)
    elseif (present(v4d)) then 
       ! 1D array
       opt=4
       arr_shape = shape(v4d)
    endif

    if (opt.eq.-1) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA:','NO OUTPUT VARIABLE SPECIFIED')
       ier = -1
       return
    endif

    ! inquire on variable NAME 
    ! Exit if it is not found
    ier = nf90_inq_varid(nc_id, var_name, var_id)
    if (ier.ne.nf90_noerr) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA: ',trim(var_name)//' NOT IN DATASET')
       ! Because of the non-stored Zero assumption present in some cases - set all returned variables to zero 
       ! - make sure ier is set so calling routine can deal with possible errors
       if (opt.eq.0) then 
          v0d = 0
       elseif (opt.eq.1) then 
          v1d = 0
       elseif (opt.eq.2) then 
          v2d = 0
       elseif (opt.eq.3) then 
          v3d = 0
       elseif (opt.eq.4) then 
          v4d = 0
       endif
       return
    endif
   
    ! check for matching shape 
    ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
    !       it becomes impossible to overwrite existing database values with zeroes
    ! This should have checked the name and shape information for var_id

    ier = check_existing_variable(var_name,NF90_FLOAT,var_id,dim_vals=arr_shape)
    if (ier.ne.nf90_noerr) return

    ! switch to data mode to read the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error reading ', var_name, ' to variable'
    ! get data ... note the lower bounds are assumed to be 1. 

    ! Read data from file
    if (opt.eq.0) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v0d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v1d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v2d))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v3d))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v4d))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION read_r4_data


  ! REAL*8

  FUNCTION read_r8_data(var_name,ndims,long_name,units,dim_names,nd,v0d,v1d,v2d,v3d,v4d) RESULT(ier)
    ! import netcdf functions
    USE netcdf, ONLY : nf90_noerr, nf90_inquire_variable, &
         & NF90_DOUBLE, nf90_get_var, nf90_max_name
    use error_handling
    IMPLICIT NONE
    ! jdemod - parameterize different shape arrays through optional arguments - should let the one routine
    !          handle any size of array
    ! jdemod - fortran interface distinguishes on TYPE only - not rank - so it is not possible to 
    !          generalize the 1D,2D etc routines since they have identical type signatures

    CHARACTER(len=*), INTENT(IN) :: var_name,
    character(len=*), intent(in),optional :: long_name, units
    integer, intent(in) :: ndims

    CHARACTER(len=*), DIMENSION(ndims), optional :: dim_names
    INTEGER, INTENT(IN), DIMENSION(ndims), optional :: nd

    REAL*8 , optional, INTENT(IN) :: v0d
    REAL*8 , DIMENSION(:), optional, INTENT(IN) :: v1d
    REAL*8 , DIMENSION(:,:), optional, INTENT(IN) :: v2d
    REAL*8 , DIMENSION(:,:,:), optional, INTENT(IN) :: v3d
    REAL*8 , DIMENSION(:,:,:,:), optional, INTENT(IN) :: v4d

    REAL*8 :: fac

    integer :: var_id,ierr

    !INTEGER :: dim_id, var_id, ier
    !INTEGER, DIMENSION(ndims) :: dim_ids
    !LOGICAL, DIMENSION(ndims) :: create_dim
    !CHARACTER(len=nf90_max_name) :: dim_name
    
    integer, dimension(ndims) :: arr_shape
    
    ! This routine reads the specified variable
    ! Almost all the specifiers are optional since they are there ONLY for checking ... checks not implemented at the moment
    
    ! The routine performs the basic check that the shape and size of the array specified
    ! matches the data to be read

    ! If the optional variables are specified they will be checked against the stored quantity as an additional
    ! level of checking. 

    ! Check to see which rank of data is to be read and set a variable so we do not need 
    ! to check present all the time. NOTE: Non-generic code that depends on the variable rank
    ! must be inside a block using the correct opt identifier

    ! Set variable option to ndims ... then double check that the expected optional argument is set

    opt = -1

    if (present(v0d)) then 
       ! scalar
       opt=0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
       arr_shape = shape(v1d)
    elseif (present(v2d)) then 
       ! 1D array
       opt=2
       arr_shape = shape(v2d)
    elseif (present(v3d)) then 
       ! 1D array
       opt=3
       arr_shape = shape(v3d)
    elseif (present(v4d)) then 
       ! 1D array
       opt=4
       arr_shape = shape(v4d)
    endif

    if (opt.eq.-1) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA:','NO OUTPUT VARIABLE SPECIFIED')
       ier = -1
       return
    endif

    ! inquire on variable NAME 
    ! Exit if it is not found
    ier = nf90_inq_varid(nc_id, var_name, var_id)
    if (ier.ne.nf90_noerr) then 
       call errmsg('NC_UTILS_DIVIMP: READ_I4_DATA: ',trim(var_name)//' NOT IN DATASET')
       ! Because of the non-stored Zero assumption present in some cases - set all returned variables to zero 
       ! - make sure ier is set so calling routine can deal with possible errors
       if (opt.eq.0) then 
          v0d = 0
       elseif (opt.eq.1) then 
          v1d = 0
       elseif (opt.eq.2) then 
          v2d = 0
       elseif (opt.eq.3) then 
          v3d = 0
       elseif (opt.eq.4) then 
          v4d = 0
       endif
       return
    endif
   
    ! check for matching shape 
    ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
    !       it becomes impossible to overwrite existing database values with zeroes
    ! This should have checked the name and shape information for var_id

    ier = check_existing_variable(var_name,NF90_DOUBLE,var_id,dim_vals=arr_shape)
    if (ier.ne.nf90_noerr) return

    ! switch to data mode to read the data
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN

    WRITE(err_msg,*) 'Error reading ', var_name, ' to variable'
    ! get data ... note the lower bounds are assumed to be 1. 

    ! Read data from file
    if (opt.eq.0) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v0d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v1d))
    elseif (opt.eq.1) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v2d))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v3d))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v4d))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION read_r8_data




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

  subroutine test_nc_utils
    USE netcdf, ONLY : nf90_create, nf90_close, nf90_noerr
    USE nc_utils, ONLY: nc_id, write_1d_array, write_2d_array, test, verbose, write_3d_array, &
         & write_scalar
    IMPLICIT NONE
    INTEGER,PARAMETER :: N = 100, N1=10, N2=8, N3=8

    LOGICAL :: exist
    REAL*8 :: scalar_r8,scalar_r8_r
    REAL*8, DIMENSION(N) :: arr1d, zeroarr1d, arr1d_r
    REAL*8, DIMENSION(N1,N2) :: arr2d, zeroarr2d, arr2d_r
    REAL*8, DIMENSION(N1,N2,N3) :: arr3d, zeroarr3d, arr3d_r

    REAL*4 :: scalar_r4,scalar_r4_r
    REAL*4, DIMENSION(N) :: arr1d_r4, zeroarr1d_r4, arr1d_r4_r
    REAL*4, DIMENSION(N1,N2) :: arr2d_r4, zeroarr2d_r4, arr2d_r4_r
    REAL*4, DIMENSION(N1,N2,N3) :: arr3d_r4, zeroarr3d_r4, arr2d_r4_r

    integer :: scalar_i4,scalar_i4_r
    integer, DIMENSION(N) :: arr1d_i4, zeroarr1d_i4, arr1d_i4_r
    integer, DIMENSION(N1,N2) :: arr2d_i4, zeroarr2d_i4, arr2d_i4_r
    integer, DIMENSION(N1,N2,N3) :: arr3d_i4, zeroarr3d_i4, arr3d_i4_r

    character*40 :: test_string,test_string1,test_string2

    CHARACTER(len=*), PARAMETER :: fn = 'test_nc_utils.nc'

    test_string = fn

    verbose = .false.
    DO i=1,N
       arr1d(i) = real(i)
    ENDDO
    DO i=1,N1
       DO j=1,N2
          arr2d(i,j) = (i**2+j**2)**0.5
          DO k=1,N3
             arr3d(i,j,k) = (i**2+j**2+k**2)**0.5
          ENDDO
       ENDDO
    ENDDO

    zeroarr1d = 0.0
    zeroarr2d = 0.0
    zeroarr3d = 0.0

    zeroarr1d_r4 = 0.0
    zeroarr2d_r4 = 0.0
    zeroarr3d_r4 = 0.0

    zeroarr1d_i4 = 0.0
    zeroarr2d_i4 = 0.0
    zeroarr3d_i4 = 0.0

    ier = open_nc_file(fn,ier,NC_WRITE,verbose)

    IF ( ier .ne. nf90_noerr ) THEN
       WRITE(*,*) 'Problem opening file', fn
       STOP
    ENDIF

    scalar_r8 = 10.0
    scalar_r4 = 4.0
    scalar_i4 = 135


    ! scalar tests
    CALL test(write_nc('test_scalar_r8','This is a scalar','W',0,v0d=scalar_r8), 'Write scalar')
    CALL test(write_nc('zero_scalar_r8','This is 0','W',0,v0d=0.0D0), 'Not write scalar 0')
    scalar_r8 = 0.0
    CALL test(write_nc('zero_scalar_r8','This is 0','W',0,v0d=scalar_r8), 'Not write scalar 0')
    CALL test(write_nc('test_scalar_r4','This is a scalar','W',0,v0d=scalar_r4), 'Write scalar')
    CALL test(write_nc('zero_scalar_r4','This is 0','W',0,v0d=0.0D0), 'Not write scalar 0')
    CALL test(write_nc('test_scalar_i4','This is a scalar','W',0,v0d=scalar_i4), 'Write scalar')
    CALL test(write_nc('zero_scalar_i4','This is 0','W',0,v0d=0.0D0), 'Not write scalar 0')
    CALL test(write_nc('test_text_string','This is a text string','T','test string'), 'Write a text string')
    CALL test(write_nc('test_text_string2','This is a text string','T',test_string), 'Write another test string')


    ! 1d tests - r8
    CALL test(write_nc('test_1d_arr','This is the long name','-',size(shape(arr1d)),'dim_arr',N,v1d=arr1d),&
         &'Write 1d array')
    CALL test(write_nc('zero_arr','This array is all 0s','-',size(shape(zeroarr1d)),'dim_arr',N,v1d=zeroarr1d),&
         &'Not write 1d array of zeros')
    CALL test(write_nc('test_1d_arr','This is the long name','-',size(shape(arr1d)),'dim_arr',N,v1d=arr1d*2),&
         &'Update 1d array with value*2')
    CALL test(write_nc('test_1d_arr','arr*10','-',size(shape(arr1d)),'dim_arr2',N,v1d=arr1d*10.0),&
         &'Not write 1d array with wrong dimension', .true.)
    CALL test(write_nc('test_1d/10_arr','N/10','-',size(shape(arr1d)),'dim_arr',N/10,v1d=arr1d(1:N/10)),&
         & 'Not write 1d array with wrong dimension length', .true.)
    CALL test(write_nc('test_1d/10_arr','arr(1:N/10)','-',size(shape(arr1d)),'dim_arr',N,v1d=arr1d(1:N/10)),&
         & 'Not write 1d array with wrong size', .true.)

    ! 1d tests - r4
    CALL test(write_nc('test_1d_arr','This is the long name','-',size(shape(arr1d_r4)),'dim_arr',N,v1d=arr1d_r4),&
         &'Write 1d array')
    CALL test(write_nc('zero_arr','This array is all 0s','-',size(shape(zeroarr1d_r4)),'dim_arr',N,v1d=zeroarr1d_r4),&
         &'Not write 1d array of zeros')
    CALL test(write_nc('test_1d_arr','This is the long name','-',size(shape(arr1d_r4)),'dim_arr',N,v1d=arr1d_r4*2),&
         &'Update 1d array with value*2')
    CALL test(write_nc('test_1d_arr','arr*10','-',size(shape(arr1d_r4)),'dim_arr2',N,v1d=arr1d_r4*10.0),&
         &'Not write 1d array with wrong dimension', .true.)
    CALL test(write_nc('test_1d/10_arr','N/10','-',size(shape(arr1d_r4)),'dim_arr',N/10,v1d=arr1d_r4(1:N/10)),&
         & 'Not write 1d array with wrong dimension length', .true.)
    CALL test(write_nc('test_1d/10_arr','arr(1:N/10)','-',size(shape(arr1d_r4)),'dim_arr',N,v1d=arr1d_r4(1:N/10)),&
         & 'Not write 1d array with wrong size', .true.)


    ! 1d tests - i4
    CALL test(write_nc('test_1d_arr','This is the long name','-',size(shape(arr1d_i4)),'dim_arr',N,v1d=arr1d_i4),&
         &'Write 1d array')
    CALL test(write_nc('zero_arr','This array is all 0s','-',size(shape(zeroarr1d_i4)),'dim_arr',N,v1d=zeroarr1d_i4),&
         &'Not write 1d array of zeros')
    CALL test(write_nc('test_1d_arr','This is the long name','-',size(shape(arr1d_i4)),'dim_arr',N,v1d=arr1d_i4*2),&
         &'Update 1d array with value*2')
    CALL test(write_nc('test_1d_arr','arr*10','-',size(shape(arr1d_i4)),'dim_arr2',N,v1d=arr1d_i4*10.0),&
         &'Not write 1d array with wrong dimension', .true.)
    CALL test(write_nc('test_1d/10_arr','N/10','-',size(shape(arr1d_i4)),'dim_arr',N/10,v1d=arr1d_i4(1:N/10)),&
         & 'Not write 1d array with wrong dimension length', .true.)
    CALL test(write_nc('test_1d/10_arr','arr(1:N/10)','-',size(shape(arr1d_i4)),'dim_arr',N,v1d=arr1d_i4(1:N/10)),&
         & 'Not write 1d array with wrong size', .true.)


    ! 2d tests - r8
    CALL test(write_nc('zero_arr_2d','All 0s','-',size(shape(zeroarr2d)),&
         & (/'zerodim01','zerodim10'/),(/N1,N2/),v2d=zeroarr2d)), 'Not write 2d array of zeros')

    CALL test(write_nc('test_2d_arr','This is a 2d array','W',size(shape(arr2d)),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d), 'Write 2d array')

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d)),&
         & (/'dim_arr3','dim_arr4'/),(/N1,N2/),v2d=arr2d), &
         & 'Not write 2d array with wrong dimension', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d)),&
         & (/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/),v2d=arr2d), &
         & 'Not write 2d array with wrong dimension length', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d(1:N1,1:N2/2))),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d(1:N1,1:N2/2) ), &
         & 'Not write 2d array with wrong size', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','W',size(shape(arr2d)),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d), &
         & 'Not write 2d array to 1d existing variable',.true.)




    ! 2d tests - r4
    CALL test(write_nc('zero_arr_2d','All 0s','-',size(shape(zeroarr2d_r4)),&
         & (/'zerodim01','zerodim10'/),(/N1,N2/),v2d=zeroarr2d_r4)), 'Not write 2d array of zeros')

    CALL test(write_nc('test_2d_arr','This is a 2d array','W',size(shape(arr2d_r4)),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d_r4), 'Write 2d array')

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d_r4)),&
         & (/'dim_arr3','dim_arr4'/),(/N1,N2/),v2d=arr2d_r4), &
         & 'Not write 2d array with wrong dimension', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d_r4)),&
         & (/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/),v2d=arr2d_r4), &
         & 'Not write 2d array with wrong dimension length', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d_r4(1:N1,1:N2/2))),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d_r4(1:N1,1:N2/2) ), &
         & 'Not write 2d array with wrong size', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','W',size(shape(arr2d_r4)),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d_r4), &
         & 'Not write 2d array to 1d existing variable',.true.)



    ! 2d tests - i4
    CALL test(write_nc('zero_arr_2d','All 0s','-',size(shape(zeroarr2d_i4)),&
         & (/'zerodim01','zerodim10'/),(/N1,N2/),v2d=zeroarr2d_i4)), 'Not write 2d array of zeros')

    CALL test(write_nc('test_2d_arr','This is a 2d array','W',size(shape(arr2d_i4)),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d_i4), 'Write 2d array')

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d_i4)),&
         & (/'dim_arr3','dim_arr4'/),(/N1,N2/),v2d=arr2d_i4), &
         & 'Not write 2d array with wrong dimension', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d_i4)),&
         & (/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/),v2d=arr2d_i4), &
         & 'Not write 2d array with wrong dimension length', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','m',size(shape(arr2d_i4(1:N1,1:N2/2))),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d_i4(1:N1,1:N2/2) ), &
         & 'Not write 2d array with wrong size', .true.)

    CALL test(write_nc('test_2d_arr','This is a 2d array','W',size(shape(arr2d_i4)),&
         & (/'dim_arr1','dim_arr2'/),(/N1,N2/),v2d=arr2d_i4), &
         & 'Not write 2d array to 1d existing variable',.true.)




    ! 3d tests -r8
    CALL test(write_nc('zero_arr_3d','All 0s','-',size(shape(zeroarr3d)),&
         & (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/),v3d=zeroarr3d), 'Not write 3d array of zeros')

    CALL test(write_nc('test_3d_arr','This is a 3d array','kg',size(shape(arr3d)),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),v3d=arr3d), 'Write 3d array')

    CALL test(write_nc('test_3d_arr','This is a 3d array','kg',size(shape(arr3d)),&
         & (/'dim_arr1','dim_arr3','dim_arr4'/),(/N1,N2,N3/),v3d=arr3d), &
         & 'Not write 3d array with wrong dimension', .true.)

    CALL test(write_nc(arr3d,'test_3d_arr','This is a 3d array','kg',size(shape(arr3d)),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1/2,N2/2,N3/2/),v3d=arr3d), &
         & 'Not write 3d array with wrong dimension length', .true.)

    CALL test(write_nc(arr3d(1:N1,1:N2/2,1:N3),'test_3d_arr','This is a 3d array','kg',size(shape(arr3d(1:N1,1:N2/2,1:N3))),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),v3d=arr3d(1:N1,1:N2/2,1:N3)), &
         & 'Not write 3d array with wrong size', .true.)

    ! 3D tests - R4
    CALL test(write_nc('zero_arr_3d','All 0s','-',size(shape(zeroarr3d_r4)),&
         & (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/),v3d=zeroarr3d_r4), 'Not write 3d array of zeros')

    CALL test(write_nc('test_3d_arr','This is a 3d array','kg',size(shape(arr3d_r4)),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),v3d=arr3d_r4), 'Write 3d array')

    CALL test(write_nc('test_3d_arr','This is a 3d array','kg',size(shape(arr3d_r4)),&
         & (/'dim_arr1','dim_arr3','dim_arr4'/),(/N1,N2,N3/),v3d=arr3d_r4), &
         & 'Not write 3d array with wrong dimension', .true.)

    CALL test(write_nc(arr3d,'test_3d_arr','This is a 3d array','kg',size(shape(arr3d_r4)),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1/2,N2/2,N3/2/),v3d=arr3d_r4), &
         & 'Not write 3d array with wrong dimension length', .true.)

    CALL test(write_nc(arr3d(1:N1,1:N2/2,1:N3),'test_3d_arr','This is a 3d array','kg',size(shape(arr3d_r4(1:N1,1:N2/2,1:N3))),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),v3d=arr3d_r4(1:N1,1:N2/2,1:N3)), &
         & 'Not write 3d array with wrong size', .true.)


    ! 3D tests - I4
    CALL test(write_nc('zero_arr_3d','All 0s','-',size(shape(zeroarr3d_i4)),&
         & (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/),v3d=zeroarr3d_i4), 'Not write 3d array of zeros')

    CALL test(write_nc('test_3d_arr','This is a 3d array','kg',size(shape(arr3d_i4)),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),v3d=arr3d_i4), 'Write 3d array')

    CALL test(write_nc('test_3d_arr','This is a 3d array','kg',size(shape(arr3d_i4)),&
         & (/'dim_arr1','dim_arr3','dim_arr4'/),(/N1,N2,N3/),v3d=arr3d_i4), &
         & 'Not write 3d array with wrong dimension', .true.)

    CALL test(write_nc(arr3d,'test_3d_arr','This is a 3d array','kg',size(shape(arr3d_i4)),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1/2,N2/2,N3/2/),v3d=arr3d_i4), &
         & 'Not write 3d array with wrong dimension length', .true.)

    CALL test(write_nc(arr3d(1:N1,1:N2/2,1:N3),'test_3d_arr','This is a 3d array','kg',size(shape(arr3d_i4(1:N1,1:N2/2,1:N3))),&
         & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/),v3d=arr3d_i4(1:N1,1:N2/2,1:N3)), &
         & 'Not write 3d array with wrong size', .true.)


    call close_nc_file


    ! Lets read back and verify the values written


    ! re-open file
    ier = open_nc_file(fn,ier,NC_READONLY,verbose)

    IF ( ier .ne. nf90_noerr ) THEN
       WRITE(*,*) 'Problem opening file for reading', fn
       STOP
    ENDIF



    ! scalar tests
    CALL test(read_nc('test_scalar_r8',0,'This is a scalar','W',0,v0d=scalar_r8_r), 'READ scalar')
    if (scalar_r8_r.ne.scalar_r8) then 
       write(0,*) 'ERROR: R8 not read correctly:',scalar_r8,scalar_r8_r
    endif

    CALL test(read_nc('test_scalar_r4',0,'This is a scalar','W',0,v0d=scalar_r4_r), 'READ scalar')
    if (scalar_r8_r.ne.scalar_r8) then 
       write(0,*) 'ERROR: R8 not read correctly:',scalar_r8,scalar_r8_r
    endif

    CALL test(read_nc('test_scalar_i4',0,'This is a scalar','W',0,v0d=scalar_i4_r), 'READ scalar')
    if (scalar_r8_r.ne.scalar_r8) then 
       write(0,*) 'ERROR: R8 not read correctly:',scalar_r8,scalar_r8_r
    endif


    CALL test(read_nc('test_scalar_r4',0,'This is a scalar','W',0,v0d=scalar_i4_r), 'READ scalar wrong type',.true.)

    CALL test(read_nc('test_text_string','This is a text string','T',test_string1), 'READ a text string')
    CALL test(read_nc('test_text_string2','This is a text string','T',test_string2), 'READ another test string')











    call close_nc_file



  END subroutine test_nc_utils





END MODULE nc_utils

