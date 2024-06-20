module init_test_mod

contains
  subroutine init_test(level,element,temperature,plasma_type)
    use all_variables
    use test_var
    use Mlog_message
    implicit none
    integer*4,intent(in) :: level
    character(len=2),intent(in) :: element
    real*8, intent(in) :: temperature
    integer*4,intent(in) :: plasma_type

    select case(test_configuration) 
    case(1)
       test_theta_shift=0.D0
       test_reg_r=1.D0
       test_reg_theta=1.D0
    case(2)
       test_theta_shift=0.25D0
       test_reg_r=0.2D0
       test_reg_theta=0.2D0
    case(3)
       test_theta_shift=0.25D0
       test_reg_r=1.D0
       test_reg_theta=1.D0
    end select


    call set_test_soledge_parameters(element,temperature)
    call set_test_mesh(level)
    call set_test_geometry()
    call allocate_masks()
    call MD_broadcast_mesh()
    call MD_broadcast_corners_mesh()
    call MD_broadcast_geometry()
    call write_log_message(msg_mesh_loaded)

    call init_super_tridiag()

    call compute_index()

    allocate(global_variables%dt_table(1:global_parameters%N_zones))
    call allocate_arrays_memory()

    call compute_reference_geometry()
    call compute_reference_parameters()

    call compute_transport_coefficients()

    call compute_jacobian()
    call compute_metric_coefficients()
    call compute_sign_metric()

    if(plasma_type.eq.1) then
       call init_test_plasma()
    else
       if(plasma_type.eq.2) then
          call init_test_plasma2()
       else
          call init_test_plasma3()
       end if
    end if

    call compute_implicit_tridiag_coefficients()

    allocate(test_sources(1:global_parameters%N_zones,0:global_parameters%N_ions))
    call allocate_test_sources()
    allocate(global_variables%residuals(0:global_parameters%N_ions))

    call load_amdata()
    call allocate_am_vars()

  end subroutine init_test
end module init_test_mod
