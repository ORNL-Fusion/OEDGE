  ! -*-Mode:f90-*-
  

  ! check for matching shape and non-zero content
  ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
  !       it becomes impossible to overwrite existing database values with zeroes

  ! Set opt based on size of argument (dimensionality)
  ! opt = 0 = scalar
  opt = size(shape(vars))

  !check_zero = vars.eq.0

  ! check matching shape
  !if (opt.gt.0) then 
     IF ( any(SHAPE(vars) .ne. nd )) THEN
        IF (verbose)  call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//': Variable and dimension shape mismatch')
        ier = nf90_noerr+10
        RETURN
     ENDIF
  !endif

  ! Add or check dimension values
  ! Check for non-scalar values
  !if (opt.gt.0) then 
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
  !endif

  ! inquire on variable NAME 
  ier = nf90_inq_varid(nc_id, var_name, var_id)

  IF (ier .ne. nf90_noerr) THEN
     ! If the variable is not found 
     ! check to make sure all values are not zero (if zero_check active) before creating a new variable

     if (zero_check) then 
        ! check not all zero
        IF (.not. any(vars.ne.0)) THEN
           IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0')
           ier = nf90_noerr
           RETURN
        ENDIF
     endif

!     IF (opt.gt.0) then 
        if (any(create_dim) ) THEN
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
        endif 
!     ENDIF

     ! switch to define mode and add the variable definition
     ier = switch_to_define_mode()
     IF (ier .ne. nf90_noerr) RETURN

     WRITE(err_msg,*) 'Problem creating variable ', var_name

!     if (opt.eq.0) then
!        ! define scalar
!        ier = handle_nf90_error(nf90_def_var(nc_id, var_name, type_val, var_id))
!     else
        ! define array
        ier = handle_nf90_error(nf90_def_var(nc_id, var_name, type_val, dim_ids, var_id))
!     endif

     IF (ier .ne. nf90_noerr) RETURN

  ELSE
     ! variable exists already ... check that it matches ... if it does then data will be over-written
     ! err_msg set by check_existing_variable

     ! dim_ids should already be a zero length array for scalars
     ier = handle_nf90_error(check_existing_variable(var_name, type_val, var_id, dim_ids))
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

  ! switch to data mode to store the data
  ier = switch_to_data_mode()
  if ( ier .ne. nf90_noerr ) RETURN

  WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
  ! add data ... note the lower bounds are assumed to be 1. 

  ! Write data to file
  ier = handle_nf90_error(nf90_put_var(nc_id,var_id,vars))
  if ( ier .ne. nf90_noerr ) RETURN


