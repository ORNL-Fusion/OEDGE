    ! -*-Mode:f90-*-
    
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
       arr_shape = 0
    elseif (present(v1d)) then 
       ! 1D array
       opt=1
       arr_shape = shape(v1d)
    elseif (present(v2d)) then 
       ! 2D array
       opt=2
       arr_shape = shape(v2d)
    elseif (present(v3d)) then 
       ! 3D array
       opt=3
       arr_shape = shape(v3d)
    elseif (present(v4d)) then 
       ! 4D array
       opt=4
       arr_shape = shape(v4d)
    endif

    if (opt.eq.-1) then 
       if (verbose) call errmsg('NC_UTILS_GENERIC: READ_'//type_str//'_DATA:','NO OUTPUT VARIABLE SPECIFIED')
       ier = -1
       return
    endif

    ! inquire on variable NAME 
    ! Exit if it is not found
    ier = nf90_inq_varid(nc_id, var_name, var_id)

    if (ier.ne.nf90_noerr) then 
       if (verbose) call errmsg('NC_UTILS_GENERIC: READ_'//type_str//'_DATA: ',trim(var_name)//' NOT IN DATASET')
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

    ier = check_existing_variable(var_name,type_val,var_id,dim_vals=arr_shape)
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
    elseif (opt.eq.2) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v2d))
    elseif (opt.eq.3) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v3d))
    elseif (opt.eq.4) then 
       ier = handle_nf90_error(nf90_get_var(nc_id,var_id,v4d))
    endif
    if ( ier .ne. nf90_noerr ) RETURN

    ! if attributes are requested then load them

    if (present(long_name)) then 
       ier = handle_nf90_error(get_long_name(var_id,long_name))
       if (ier.ne.nf90_noerr) then
          long_name = ''
       endif
    endif

    if (present(units)) then 
       ier = handle_nf90_error(get_units(var_id,units,fac))
       if (ier.ne.nf90_noerr) then
          units = ''
       endif
    endif

