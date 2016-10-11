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

MODULE nc_utils

    implicit none

    PRIVATE
    PUBLIC ::  write_scalar, write_1d_array, write_2d_array, write_3d_array, err_msg, test
!    public ::  open_nc_file,close_nc_file

    CHARACTER(len=1024) :: err_msg ! Holder of error messages
    INTEGER :: nc_id
    LOGICAL :: verbose = .true.


!    interface write_scalar
!       module procedure write_i4_scalar,write_r4_scalar,write_r8_scalar
!    end interface write_scalar

!    interface write_1d_array
!       module procedure write_1d_i4_array,write_1d_r4_array,write_1d_r8_array
!    end interface write_1d_array

!    interface write_2d_array
!       module procedure write_2d_i4_array,write_2d_r4_array,write_2d_r8_array
!    end interface write_2d_array

!    interface write_3d_array
!       module procedure write_3d_i4_array,write_3d_r4_array,write_3d_r8_array
!    end interface write_3d_array



CONTAINS

  subroutine open_nc_file(filename,debug,ierr)
    ! opens the nc file ... sets nc_id and verbose output/debug
    use error_handling
    use netcdf, only : nf90_open
    implicit none
    
    character*(*) :: filename
    logical :: debug 
    integer :: ierr
    
    verbose = debug

    ierr =  nf90_open(trim(filename),0,nc_id)
    
    if (ierr.ne.0) then 
       call errmsg('OPEN_NC_FILE:ERROR OPENING NC FILE:'//trim(filename),ierr)
       stop 'OPEN_NC_FILE'
    endif
    return
  end subroutine open_nc_file


  FUNCTION handle_nf90_error(ier_in) RESULT(ier)
    USE netcdf, ONLY : NF90_NOERR, nf90_close
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ier_in
    INTEGER :: ier, ier_close
    ier = ier_in
    IF ( ier .ne. NF90_NOERR ) THEN
      IF ( verbose ) WRITE(*,*) TRIM(err_msg)
    ENDIF
  END FUNCTION handle_nf90_error

  FUNCTION check_existing_variable(var_name, var_type, var_id, &
      & dim_ids) RESULT(ier)
    USE netcdf, ONLY : nf90_inquire_variable, NF90_MAX_NAME, &
      & NF90_MAX_VAR_DIMS, nf90_noerr
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: var_type, var_id
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), DIMENSION(:) :: dim_ids
    CHARACTER(len=NF90_MAX_NAME) :: existing_var_name
    INTEGER ::  existing_type, existing_ndims, ind, ier
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: existing_dimids
    ier = nf90_inquire_variable(nc_id, var_id, existing_var_name, &
      & existing_type, existing_ndims, existing_dimids)
    IF ( ier .ne. NF90_NOERR ) THEN
      WRITE(err_msg,*) 'Unable to check existing variable ',var_name
      RETURN
    ENDIF
    IF ( TRIM(existing_var_name) .ne. TRIM(var_name) ) THEN
      WRITE(err_msg,*) 'Variable name is different ',TRIM(var_name),&
      & TRIM(existing_var_name)
      ier = NF90_NOERR + 1
      RETURN
    ENDIF
    IF ( existing_type .ne. var_type ) THEN
      WRITE(err_msg,*) 'Type mismatch for variable ',TRIM(var_name)
      ier = NF90_NOERR + 2
      RETURN
    ENDIF
    IF ( existing_ndims .ne. SIZE(dim_ids) ) THEN
      WRITE(err_msg,*) 'Dimension mismatch for variable ',TRIM(var_name)
      ier = NF90_NOERR + 3
      RETURN
    ENDIF
    DO ind = 1, SIZE(dim_ids)
      IF ( dim_ids(ind) .ne. existing_dimids(ind) ) THEN
        WRITE(err_msg,*) 'Attempt to write variable ',TRIM(var_name), ' with existing dimension different than new dimension'
        ier = NF90_NOERR + 4
        RETURN
      ENDIF
    END DO
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
    IF (TRIM(dim_name) .ne. TRIM(existing_name) ) THEN
      ier = nf90_noerr + 1
      WRITE(err_msg,*) 'Dim id ',dim_id, ' name mismatch ', TRIM(dim_name), ' ', TRIM(existing_name)
      RETURN
    ENDIF
    IF (dim_len .ne. existing_len) THEN
      ier = nf90_noerr + 2
      WRITE(err_msg,*) 'Dimension ', TRIM(dim_name), ' is already of size ', &
              & existing_len, ' not ', dim_len
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
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id, 'units',&
      & TRIM(units)))
    IF ( ier .ne. nf90_noerr ) RETURN
    WRITE(err_msg,*) 'Error adding long_name attribute',TRIM(long_name)
    ier = handle_nf90_error(nf90_put_att(nc_id, var_id, 'long_name',&
      & TRIM(long_name)))
    RETURN
  END FUNCTION

  FUNCTION switch_to_define_mode() RESULT(ier)
    ! switching to define mode will not return an error
    USE netcdf, ONLY : nf90_redef, nf90_noerr
    IMPLICIT NONE
    INTEGER :: ier
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
    WRITE(err_msg,*) 'error switching to data mode'
    ier = nf90_enddef(nc_id)
    ier = nf90_noerr
    RETURN
  END FUNCTION switch_to_data_mode

  FUNCTION write_scalar(var,var_name,long_name,units) RESULT(ier)
    USE netcdf, ONLY : nf90_noerr, &
        & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var
    IMPLICIT NONE
    REAL*8 , INTENT(IN) :: var
    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units
    REAL*8 :: fac
    INTEGER :: var_id, ier, i
    INTEGER, ALLOCATABLE, DIMENSION(:) :: zero_len_array ! This may not be portable
    ALLOCATE(zero_len_array(0))

    IF (var .eq. 0) THEN
      IF (verbose) WRITE(*,*) var_name, ' not written to netcdf file because it is 0'
      ier = nf90_noerr
      RETURN
    ENDIF
    ier = nf90_inq_varid(nc_id, var_name, var_id)
    IF (ier .ne. nf90_noerr) THEN
      ier = switch_to_define_mode()
      IF (ier .ne. nf90_noerr) RETURN
      WRITE(err_msg,*) 'Problem creating variable ', var_name
      ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_DOUBLE, var_id))
      IF (ier .ne. nf90_noerr) RETURN
    ELSE
      ! err_msg set by check_existing_variable
      ier = handle_nf90_error(check_existing_variable(var_name, &
        & NF90_DOUBLE, var_id, zero_len_array))
      IF (ier .ne. nf90_noerr) RETURN
    ENDIF
    fac = 1.
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN
    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ier = handle_nf90_error(nf90_put_var(nc_id,var_id,var))
    if ( ier .ne. nf90_noerr ) RETURN
    RETURN

  END FUNCTION write_scalar

  FUNCTION write_1d_array(var,var_name,long_name,units,dim_name,nd) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
        & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var
    IMPLICIT NONE
    REAL*8 , DIMENSION(nd), INTENT(IN) :: var
    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units, dim_name
    INTEGER, INTENT(IN) :: nd
    REAL*8 :: fac
    INTEGER :: dim_id, var_id, ier
    LOGICAL :: create_dim

    IF (.not. any(var .ne. 0)) THEN
      IF (verbose) WRITE(*,*) var_name, ' not written to netcdf file because it is all 0'
      ier = nf90_noerr
      RETURN
    ENDIF
    ier = nf90_inq_dimid(nc_id, dim_name, dim_id)
    IF (ier .ne. nf90_noerr) THEN
      create_dim = .true.
    ELSE
      create_dim = .false.
      ier = handle_nf90_error(check_existing_dimension(dim_id, dim_name, nd))
      IF (ier .ne. nf90_noerr) RETURN
    ENDIF
    ier = nf90_inq_varid(nc_id, var_name, var_id)
    IF (ier .ne. nf90_noerr) THEN
      IF ( create_dim ) THEN
        ier = switch_to_define_mode()
        IF (ier .ne. nf90_noerr) RETURN
        WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
        ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd, dim_id))
        IF (ier .ne. nf90_noerr) RETURN
      ENDIF
      ier = switch_to_define_mode()
      IF (ier .ne. nf90_noerr) RETURN
      WRITE(err_msg,*) 'Problem creating variable ', var_name
      ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_DOUBLE, (/dim_id/), &
        & var_id))
      IF (ier .ne. nf90_noerr) RETURN
    ELSE
      ! err_msg set by check_existing_variable
      ier = handle_nf90_error(check_existing_variable(var_name, &
        & NF90_DOUBLE, var_id, (/dim_id/)))
      IF (ier .ne. nf90_noerr) RETURN
    ENDIF
    fac = 1.
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN
    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ier = handle_nf90_error(nf90_put_var(nc_id,var_id,var(1:nd)))
    if ( ier .ne. nf90_noerr ) RETURN

  END FUNCTION write_1d_array

  FUNCTION write_2d_array(var,var_name,long_name,units,dim_names,nd) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
        & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    IMPLICIT NONE
    INTEGER, INTENT(IN), DIMENSION(2) :: nd
    REAL*8 , DIMENSION(:,:), INTENT(IN) :: var
    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units
    CHARACTER(len=*), DIMENSION(2) :: dim_names
    REAL*8 :: fac
    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(2) :: dim_ids
    LOGICAL, DIMENSION(2) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name
    IF ( any(SHAPE(var) .ne. nd )) THEN
      IF (verbose) WRITE(*,*) 'Variable and dimension shape mismatch ',var_name
      ier = nf90_noerr+10
      RETURN
    ENDIF

    IF (.not. any(var .ne. 0)) THEN
      IF (verbose) WRITE(*,*) var_name, ' not written to netcdf file because it is all 0'
      ier = nf90_noerr
      RETURN
    ENDIF
    DO d = 1,2
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
    ier = nf90_inq_varid(nc_id, var_name, var_id)
    IF (ier .ne. nf90_noerr) THEN
      IF ( any(create_dim) ) THEN
        ier = switch_to_define_mode()
        IF (ier .ne. nf90_noerr) RETURN
        DO d=1,size(nd)
          IF (create_dim(d)) THEN
            dim_name = TRIM(dim_names(d))
            WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
            ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd(d), dim_ids(d)))
            IF (ier .ne. nf90_noerr) RETURN
          ENDIF
        ENDDO
      ENDIF
      ier = switch_to_define_mode()
      IF (ier .ne. nf90_noerr) RETURN
      WRITE(err_msg,*) 'Problem creating variable ', var_name
      ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_DOUBLE, dim_ids, &
        & var_id))
      IF (ier .ne. nf90_noerr) RETURN
    ELSE
      ! err_msg set by check_existing_variable
      ier = handle_nf90_error(check_existing_variable(var_name, &
        & NF90_DOUBLE, var_id, dim_ids))
      IF (ier .ne. nf90_noerr) RETURN
    ENDIF
    fac = 1.
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN
    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ier = handle_nf90_error(nf90_put_var(nc_id,var_id,var(1:nd(1),1:nd(2))))
    if ( ier .ne. nf90_noerr ) RETURN
  END FUNCTION write_2d_array
  
  FUNCTION write_3d_array(var,var_name,long_name,units,dim_names,nd) RESULT(ier)
    USE netcdf, ONLY : nf90_inq_dimid, nf90_noerr, nf90_def_dim, &
        & nf90_inq_varid, nf90_def_var, NF90_DOUBLE, nf90_put_var, nf90_max_name
    IMPLICIT NONE
    INTEGER, INTENT(IN), DIMENSION(3) :: nd
    REAL*8 , DIMENSION(:,:,:), INTENT(IN) :: var
    CHARACTER(len=*), INTENT(IN) :: var_name, long_name, units
    CHARACTER(len=*), DIMENSION(3) :: dim_names
    REAL*8 :: fac
    INTEGER :: dim_id, var_id, ier, d, ndum
    INTEGER, DIMENSION(3) :: dim_ids
    LOGICAL, DIMENSION(3) :: create_dim
    CHARACTER(len=nf90_max_name) :: dim_name
    IF ( any(SHAPE(var) .ne. nd )) THEN
      IF (verbose) WRITE(*,*) 'Variable and dimension shape mismatch ',var_name
      ier = nf90_noerr+10
      RETURN
    ENDIF
    IF (.not. any(var .ne. 0)) THEN
      IF (verbose) WRITE(*,*) var_name, ' not written to netcdf file because it is all 0'
      ier = nf90_noerr
      RETURN
    ENDIF
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
    ier = nf90_inq_varid(nc_id, var_name, var_id)
    IF (ier .ne. nf90_noerr) THEN
      IF ( any(create_dim) ) THEN
        ier = switch_to_define_mode()
        IF (ier .ne. nf90_noerr) RETURN
        DO d=1,size(nd)
          IF (create_dim(d)) THEN
            dim_name = TRIM(dim_names(d))
            WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
            ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd(d), dim_ids(d)))
            IF (ier .ne. nf90_noerr) RETURN
          ENDIF
        ENDDO
      ENDIF
      ier = switch_to_define_mode()
      IF (ier .ne. nf90_noerr) RETURN
      WRITE(err_msg,*) 'Problem creating variable ', var_name
      ier = handle_nf90_error(nf90_def_var(nc_id, var_name, NF90_DOUBLE, dim_ids, &
        & var_id))
      IF (ier .ne. nf90_noerr) RETURN
    ELSE
      ! err_msg set by check_existing_variable
      ier = handle_nf90_error(check_existing_variable(var_name, &
        & NF90_DOUBLE, var_id, dim_ids))
      IF (ier .ne. nf90_noerr) RETURN
    ENDIF
    fac = 1.
    ier = add_units_long_name(var_id, units, long_name, fac)
    if ( ier .ne. nf90_noerr ) RETURN
    ier = switch_to_data_mode()
    if ( ier .ne. nf90_noerr ) RETURN
    WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
    ier = handle_nf90_error(nf90_put_var(nc_id,var_id,var(1:nd(1),1:nd(2),1:nd(3))))
    if ( ier .ne. nf90_noerr ) RETURN
  END FUNCTION write_3d_array
  
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
END MODULE nc_utils

!#ifdef TEST_NC_UTILS
!PROGRAM test_nc_utils
!  USE netcdf, ONLY : nf90_create, nf90_close, nf90_noerr
!  USE nc_utils, ONLY: nc_id, write_1d_array, write_2d_array, test, verbose, write_3d_array, &
!    & write_scalar
!  IMPLICIT NONE
!  INTEGER,PARAMETER :: N = 100, N1=10, N2=8, N3=8
!
!  LOGICAL :: exist
!  REAL*8 :: scalar
!  REAL*8, DIMENSION(N) :: arr, zeroarr
!  REAL*8, DIMENSION(N1,N2) :: arr2d, zeroarr2d
!  REAL*8, DIMENSION(N1,N2,N3) :: arr3d, zeroarr3d
!  CHARACTER(len=*), PARAMETER :: fn = 'test_nc_utils.nc'
!  verbose = .false.
!  DO i=1,N
!    arr(i) = real(i)
!  ENDDO
!  DO i=1,N1
!    DO j=1,N2
!      arr2d(i,j) = (i**2+j**2)**0.5
!      DO k=1,N3
!        arr3d(i,j,k) = (i**2+j**2+k**2)**0.5
!      ENDDO 
!    ENDDO
!  ENDDO
!  zeroarr(:) = 0.0
!  zeroarr2d(:,:) = 0.0
!  zeroarr3d(:,:,:) = 0.0
!  ier = nf90_create(fn,0,nc_id)
!  IF ( ier .ne. nf90_noerr ) THEN
!    WRITE(*,*) 'Problem opening file', fn
!    STOP
!  ENDIF
!  ! scalar tests
!  CALL test(write_scalar(10.D0,'test_scalar','This is a scalar','W'), 'Write scalar')
!  CALL test(write_scalar(0.D0,'zero_scalar','This is 0','W'), 'Not write scalar 0')
!  ! 1d tests
!  CALL test(write_1d_array(arr,'test_1d_arr','This is the long name','-','dim_arr',N),&
!                    & 'Write 1d array')
!  CALL test(write_1d_array(zeroarr,'zero_arr','This array is all 0s','-','dim_arr',N),&
!                    & 'Not write 1d array of zeros')
!  CALL test(write_1d_array(arr*2,'test_1d_arr','This is the long name','-','dim_arr',N),&
!                    & 'Update 1d array with value*2')
!  CALL test(write_1d_array(arr*10., 'test_1d_arr','arr*10','-','dim_arr2',N),&
!                    & 'Not write 1d array with wrong dimension', .true.)
!  CALL test(write_1d_array(arr(1:N/10), 'test_1d/10_arr','N/10','-','dim_arr',N/10),&
!                    & 'Not write 1d array with wrong dimension length', .true.)
!  CALL test(write_1d_array(arr(1:N/10), 'test_1d/10_arr','arr(1:N/10)','-','dim_arr',N),&
!                    & 'Not write 1d array with wrong size', .true.)
  ! 2d tests
!  CALL test(write_2d_array(zeroarr2d,'zero_arr_2d','All 0s','-',&
!                    & (/'zerodim01','zerodim10'/),(/N1,N2/)), 'Not write 2d array of zeros')
!  CALL test(write_2d_array(arr2d,'test_2d_arr','This is a 2d array','W',&
!                    & (/'dim_arr1','dim_arr2'/),(/N1,N2/)), 'Write 2d array')
!  CALL test(write_2d_array(arr2d,'test_2d_arr','This is a 2d array','m',&
!                    & (/'dim_arr3','dim_arr4'/),(/N1,N2/)), &
!                    & 'Not write 2d array with wrong dimension', .true.)
!  CALL test(write_2d_array(arr2d,'test_2d_arr','This is a 2d array','m',&
!                    & (/'dim_arr1','dim_arr2'/),(/N1/2,N2/2/)), &
!                    & 'Not write 2d array with wrong dimension length', .true.)
!  CALL test(write_2d_array(arr2d(1:N1,1:N2/2),'test_2d_arr','This is a 2d array','m',&
!                    & (/'dim_arr1','dim_arr2'/),(/N1,N2/)), &
!                    & 'Not write 2d array with wrong size', .true.)
!  CALL test(write_2d_array(arr2d,'test_1d_arr','This is a 2d array','W',&
!                    & (/'dim_arr1','dim_arr2'/),(/N1,N2/)), &
!                    & 'Not write 2d array to 1d existing variable',.true.)
!  ! 3d tests
!  CALL test(write_3d_array(zeroarr3d,'zero_arr_3d','All 0s','-',&
!                    & (/'zerodim01','zerodim10','zerodim03'/),(/N1,N2,N3/)), 'Not write 3d array of zeros')
!  CALL test(write_3d_array(arr3d,'test_3d_arr','This is a 3d array','kg',&
!                    & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/)), 'Write 3d array')
!  CALL test(write_3d_array(arr3d,'test_3d_arr','This is a 3d array','kg',&
!                    & (/'dim_arr1','dim_arr3','dim_arr4'/),(/N1,N2,N3/)), &
!                    & 'Not write 3d array with wrong dimension', .true.)
!  CALL test(write_3d_array(arr3d,'test_3d_arr','This is a 3d array','kg',&
!                    & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1/2,N2/2,N3/2/)), &
!                    & 'Not write 3d array with wrong dimension length', .true.)
!  CALL test(write_3d_array(arr3d(1:N1,1:N2/2,1:N3),'test_3d_arr','This is a 3d array','kg',&
!                    & (/'dim_arr1','dim_arr2','dim_arr3'/),(/N1,N2,N3/)), &
!                    & 'Not write 3d array with wrong size', .true.)
!  ier = nf90_close(nc_id)
!END
!#endif

