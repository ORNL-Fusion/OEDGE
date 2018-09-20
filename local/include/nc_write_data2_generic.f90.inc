  ! -*-Mode:f90-*-

  ier = nf90_noerr
  
  !write(0,'(5a,i8,8(3a,i8))')   'nc_write:',rank_str,':',type_str,':',size(shape(vars)),((':',dim_names(i),'=',nd(i)),i=1,size(nd))  
  
  ! check for matching shape and non-zero content
  ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
  !       it becomes impossible to overwrite existing database values with zeroes
  
  ! Set opt based on size of argument (dimensionality)
  ! opt = 0 = scalar
  opt = size(shape(vars))

  if (verbose_debug)  write(0,*) trim(var_name)//':'//type_str//':'//rank_str//' being written'
  
  ! check matching shape
  IF ( any(SHAPE(vars) .ne. nd )) THEN
     if (verbose_debug) write(0,*) trim(var_name)//': variable and dimension shape mismatch ',ier
     IF (verbose)  call errmsg('NC_UTILS_GENERIC: WRITE_'//type_str//'_DATA: ',trim(var_name)//': Variable and dimension shape mismatch')
     ier = nf90_noerr+10
     RETURN
  ENDIF

  if (verbose_debug) write(0,*) trim(var_name)//': before inquire variable name',ier

  ! inquire on variable NAME 
  ier = nf90_inq_varid(nc_id, var_name, var_id)

  if (verbose_debug) write(0,*) trim(var_name)//': inquire variable name',ier

  IF (ier .ne. nf90_noerr) THEN
     ! If the variable is not found 
     ! check to make sure all values are not zero (if zero_check active) before creating a new variable
     if (verbose_debug) write(0,*) trim(var_name)//': variable not found',ier

     if (zero_check) then 
        ! check not all zero
        IF (.not. any(vars.ne.0)) THEN
           if (verbose_debug) write(0,*) trim(var_name)//': ALL ZERO',ier
           IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0',msglvl=1)
           ier = nf90_noerr
           RETURN
        ENDIF
     endif

     ! The create_dims paradigm actually adds more work ... dimensions only need to be checked for new variables ... existing variables can't be 
     ! put in without dimensions. This will also reduce complexity when trying to deal with cases where dimensions are repeated and needs creating.
     ! The first instance will create it and the next will check it. 
     ! Add or check dimension values
     DO d = 1,size(nd)
        dim_name = TRIM(dim_names(d))
        ier = nf90_inq_dimid(nc_id, dim_name, dim_ids(d))
        if (verbose_debug) write(0,*) trim(dim_name)//': inquire dim name',ier

        ! makes more sense to just create missing dimensions inline
        IF (ier .ne. nf90_noerr) THEN
           ! create dimension  entry
           ier = switch_to_define_mode()
           dim_name = TRIM(dim_names(d))
           IF (ier .ne. nf90_noerr) then
              call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',dim_name//' switch to define mode failed')
              RETURN
           endif              
           WRITE(err_msg,*) var_name, ' error creating dimension ', dim_name
           ! create dimensions in netcdf dataset if they are missing
           ier = handle_nf90_error(nf90_def_dim(nc_id, dim_name, nd(d), dim_ids(d)))
           if (verbose_debug) write(0,*) trim(dim_name)//': created dimension',ier,d,nd(d)
           IF (ier .ne. nf90_noerr) RETURN
        ELSE
           ! Check existing dimension matches
           ier = handle_nf90_error(check_existing_dimension(var_name,dim_ids(d), dim_name, nd(d)))
           if (verbose_debug) write(0,*) trim(dim_name)//': check dimension',ier,d,nd(d)
           IF (ier .ne. nf90_noerr) RETURN
        ENDIF
     ENDDO 

     ! switch to define mode and add the variable definition
     ier = switch_to_define_mode()
     if (verbose_debug) write(0,*) trim(var_name)//': switch to define mode',ier


     IF (ier .ne. nf90_noerr) then
        call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',trim(var_name)//' switch to define mode failed')
        RETURN
     endif              

     WRITE(err_msg,*) 'Problem creating variable ', var_name,' with rank '//rank_str//' and type '//type_str
     ier = handle_nf90_error(nf90_def_var(nc_id, var_name, type_val, dim_ids, var_id))
     if (verbose_debug) write(0,*) trim(var_name)//': create variable',ier

     IF (ier .ne. nf90_noerr) RETURN

  ELSE
     ! variable exists already ... check that it matches ... if it does then data will be over-written
     !ier = handle_nf90_error(check_existing_variable(var_name, type_val, var_id, dim_ids))
     ier = handle_nf90_error(check_existing_variable(var_name, type_val, var_id, dim_names, nd))
     if (verbose_debug) write(0,*) trim(var_name)//': variable exists',ier
     IF (ier .ne. nf90_noerr) RETURN

  ENDIF

  if (verbose_debug) write(0,*) trim(var_name)//': added variable',ier


  
  fac = 1.

  
  ! add units and long name attribute if present
  if (present(units)) then 
     ier = add_units(var_id,units,fac)
     IF (ier .ne. nf90_noerr) then
        call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',trim(var_name)//' add_units failed')
        RETURN
     endif              
     if (verbose_debug) write(0,*) trim(var_name)//': added units',ier
  endif



  if (present(long_name)) then 
     ier = add_long_name(var_id,long_name)
     IF (ier .ne. nf90_noerr) then
        call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',trim(var_name)//' add_long_name failed')
        RETURN
     endif              
     if (verbose_debug) write(0,*) trim(var_name)//': added long name',ier
  endif

  ! switch to data mode to store the data
  ier = switch_to_data_mode()
  if (verbose_debug) write(0,*) trim(var_name)//': switch to data mode',ier
  IF (ier .ne. nf90_noerr) then
     call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',trim(var_name)//' switch to data mode failed')
     RETURN
  endif              

  WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
  ! add data ... note the lower bounds are assumed to be 1. 

  ! Write data to file
  ier = handle_nf90_error(nf90_put_var(nc_id,var_id,vars))
  if (verbose_debug) write(0,*) trim(var_name)//': data written',ier

  IF (ier .ne. nf90_noerr) then
     call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',trim(var_name)//' error writing data')
     RETURN
  endif              


