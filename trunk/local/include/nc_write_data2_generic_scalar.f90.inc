  ! -*-Mode:f90-*-
  

  ! check for matching shape and non-zero content
  ! NOTE: Check for non-zero content should only be done for NEW variables ... not for existing ones ... otherwise
  !       it becomes impossible to overwrite existing database values with zeroes

  ! Set opt based on size of argument (dimensionality)
  ! opt = 0 = scalar
  opt = size(shape(vars))

  !check_zero = vars.eq.0
  if (verbose_debug) write(0,*) trim(var_name)//': scalar variable being written'

  ! inquire on variable NAME 
  ier = nf90_inq_varid(nc_id, var_name, var_id)

  if (verbose_debug) write(0,*) trim(var_name)//': inquire variable name',ier


  IF (ier .ne. nf90_noerr) THEN
     ! If the variable is not found 
     ! check to make sure all values are not zero (if zero_check active) before creating a new variable
     if (verbose_debug) write(0,*) trim(var_name)//': variable not found',ier

     if (zero_check) then 
        ! check not all zero
        IF (vars.eq.0) THEN
           if (verbose_debug) write(0,*) trim(var_name)//': VAR=ZERO',ier
           IF (verbose) call errmsg('NC_UTILS_GENERIC: WRITE_'//rank_str//'_'//type_str//'_DATA: ',var_name//' not written to netcdf file because it is all 0',msglvl=1)
           ier = nf90_noerr
           RETURN
        ENDIF
     endif

     ! switch to define mode and add the variable definition
     ier = switch_to_define_mode()
     if (verbose_debug) write(0,*) trim(var_name)//': switch to define mode',ier
     IF (ier .ne. nf90_noerr) RETURN

     WRITE(err_msg,*) 'Problem creating variable ', var_name

     ! define scalar
     ier = handle_nf90_error(nf90_def_var(nc_id, var_name, type_val, var_id))
     if (verbose_debug) write(0,*) trim(var_name)//': create variable',ier

     IF (ier .ne. nf90_noerr) RETURN

  ELSE
     ! variable exists already ... check that it matches ... if it does then data will be over-written
     ! err_msg set by check_existing_variable

     ! dim_ids should already be a zero length array for scalars
     ier = handle_nf90_error(check_existing_variable(var_name, type_val, var_id))
     if (verbose_debug) write(0,*) trim(var_name)//': check variable',ier
     IF (ier .ne. nf90_noerr) RETURN

  ENDIF

  fac = 1.

  ! add units and long name attribute if present
  if (present(units)) then 
     ier = add_units(var_id,units,fac)
     if (verbose_debug) write(0,*) trim(var_name)//': added units',ier
     if ( ier .ne. nf90_noerr ) RETURN
  endif

  if (present(long_name)) then 
     ier = add_long_name(var_id,long_name)
     if (verbose_debug) write(0,*) trim(var_name)//': added long name',ier
     if ( ier .ne. nf90_noerr ) RETURN
  endif

  ! switch to data mode to store the data
  ier = switch_to_data_mode()
  if (verbose_debug) write(0,*) trim(var_name)//': switch to data mode',ier
  IF (ier .ne. nf90_noerr) then
     call errmsg('NC_UTILS_GENERIC: WRITE_SCALAR_DATA: ',trim(var_name)//' switch to data mode failed')
     RETURN
  endif              

  WRITE(err_msg,*) 'Error writing ', var_name, ' to file'
  ! add data ... note the lower bounds are assumed to be 1. 

  ! Write data to file
  ier = handle_nf90_error(nf90_put_var(nc_id,var_id,vars))

  if (verbose_debug) write(0,*) trim(var_name)//': data written?',ier

  IF (ier .ne. nf90_noerr) then
     call errmsg('NC_UTILS_GENERIC: WRITE_SCALAR_DATA: ',trim(var_name)//' error writing data')
     RETURN
  endif              

