  ! -*-Mode:f90-*-
  
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
     ! 2D array
     opt=2
  elseif (present(v3d)) then 
     ! 3D array
     opt=3
  elseif (present(v4d)) then 
     ! 4D array
     opt=4
  endif


  if (opt.gt.0) then 
     if (.not.(present(nd).and.present(dim_names))) then 
     ! Error - an array specified but nd or dim_names is NOT
        if (verbose) call errmsg ('NC_UTILS_GENERIC:WRITE_'//type_str//'_DATA:','DIMENSIONS (ND and/or DIM_NAMES) NOT PRESENT FOR ARRAY VARIABLE')
        ier = 1
        return
     endif
  endif

  ! check for matching shape and non-zero content
  ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
  !       it becomes impossible to overwrite existing database values with zeroes

  if (opt.eq.1) then 
     ! check matching shape
     IF ( any(SHAPE(v1d) .ne. nd )) THEN
        IF (verbose)  call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//': Variable and dimension shape mismatch')
        ier = nf90_noerr+10
        RETURN
     ENDIF
  elseif (opt.eq.2) then 
     ! check matching shape
     IF ( any(SHAPE(v2d) .ne. nd )) THEN
        IF (verbose)  call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//': Variable and dimension shape mismatch')
        ier = nf90_noerr+10
        RETURN
     ENDIF
  elseif (opt.eq.3) then 
     ! check matching shape
     IF ( any(SHAPE(v3d) .ne. nd )) THEN
        IF (verbose)  call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//': Variable and dimension shape mismatch')
        ier = nf90_noerr+10
        RETURN
     ENDIF
  elseif (opt.eq.4) then 
     ! check matching shape
     IF ( any(SHAPE(v4d) .ne. nd )) THEN
        IF (verbose)  call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//': Variable and dimension shape mismatch')
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
              IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0')
              ier = nf90_noerr
              RETURN
           ENDIF
        elseif (opt.eq.1) then 
           ! check not all zero
           IF (.not. any(v1d .ne. 0)) THEN
              IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0')
              ier = nf90_noerr
              RETURN
           ENDIF
        elseif (opt.eq.2) then 
           ! check not all zero
           IF (.not. any(v2d .ne. 0)) THEN
              IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0')
              ier = nf90_noerr
              RETURN
           ENDIF
        elseif (opt.eq.3) then 
           ! check not all zero
           IF (.not. any(v3d .ne. 0)) THEN
              IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0')
              ier = nf90_noerr
              RETURN
           ENDIF
        elseif (opt.eq.4) then 
           ! check not all zero
           IF (.not. any(v4d .ne. 0)) THEN
              IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0')
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
        ier = handle_nf90_error(nf90_def_var(nc_id, var_name, type_val, var_id))
     else
        ! define array
        ier = handle_nf90_error(nf90_def_var(nc_id, var_name, type_val, dim_ids, var_id))
     endif

     IF (ier .ne. nf90_noerr) RETURN

  ELSE
     ! variable exists already ... check that it matches ... if it does then data will be over-written
     ! err_msg set by check_existing_variable

     ! dim_ids should already be a zero length array for scalars
     ier = handle_nf90_error(check_existing_variable(var_name, type_val, var_id, dim_ids))
     IF (ier .ne. nf90_noerr) RETURN

  ENDIF

  fac = 1.

  if (present(units)) then 
     ier = add_units(var_id,units,fac)
     if ( ier .ne. nf90_noerr ) RETURN
  endif

  if (present(long_name)) then 
     ier = add_long_name(var_id,long_name)
     if ( ier .ne. nf90_noerr ) RETURN
  endif


  ! add units and long name attribute
  ! ier = add_units_long_name(var_id, units, long_name, fac)
  !if ( ier .ne. nf90_noerr ) RETURN

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
  elseif (opt.eq.2) then 
     ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v2d(1:nd(1),1:nd(2))))
  elseif (opt.eq.3) then 
     ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v3d(1:nd(1),1:nd(2),1:nd(3))))
  elseif (opt.eq.4) then 
     ier = handle_nf90_error(nf90_put_var(nc_id,var_id,v4d(1:nd(1),1:nd(2),1:nd(3),1:nd(4))))
  endif
  if ( ier .ne. nf90_noerr ) RETURN


