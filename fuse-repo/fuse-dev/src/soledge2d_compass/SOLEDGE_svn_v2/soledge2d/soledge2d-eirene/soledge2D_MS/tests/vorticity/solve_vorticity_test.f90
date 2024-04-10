subroutine solve_vorticity_test()
  use all_variables, only : global_parameters, zones, flags
  use Mvorticity_vars
  use Mdefinitions
  implicit none
  integer*4 :: k

  call compute_vorticity_matrix4(DIRICHLET_TEST)
  call add_vorticity_matrix_BC4(DIRICHLET_TEST) !Dirichlet on the wall

  write(*,*) mat_vort_nnz
  call fill_implicit_pastix(CSC_vort,mat_vort_nnz,vorticity_mat)
  call free_all(vorticity_mat)
  call analyze_implicit_pastix(CSC_vort)

  do k=1,global_parameters%N_Zones
     call MD_broadcast_corners_phi2(zones(k),STEP_OLD)
     call compute_vorticity(zones(k),STEP_OLD,DIRICHLET_TEST)
     call init_RHS_vorticity(zones(k))
     call add_pi_time_plus_one_to_RHS(zones(k),DIRICHLET_TEST)
!     call add_RHS_vorticity_sources(zones(k))
     call add_RHS_pe_Ohm_term(zones(k))
     call add_RHS_test_vorticity_source(zones(k))
  end do

  call compute_vorticity_RHS4(DIRICHLET_TEST) !Dirichlet
  call solve_implicit_pastix(CSC_vort)
  call redistribute_vorticity_sol(CSC_vort)
  call clear_pastix(CSC_vort)

  do k=1,global_parameters%N_Zones
     call add_phi_Ohm_term(zones(k))
     call MD_broadcast_corners_phi2(zones(k),STEP_NEW)
     call compute_vorticity(zones(k),STEP_NEW,DIRICHLET_TEST)
  end do

end subroutine solve_vorticity_test
