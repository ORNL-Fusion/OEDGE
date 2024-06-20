subroutine solve_vorticity()
#include "compile_opt.inc"
  use all_variables, only : global_parameters, zones, flags, drift_flags
  use Mvorticity_vars
  use Mdefinitions
  implicit none
  integer*4 :: k

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)

  !$OMP DO
  do k=1,global_parameters%N_Zones
     call MD_broadcast_plasma(zones(k),STEP_NEW)
     call MD_broadcast_corners2(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER

  !$OMP MASTER
  if((.not.flags%restart).or.(drift_flags%BC_type_phi.eq.1)) then
     call compute_vorticity_matrix4(ZERO_FLUX,LAMBDA_TE)
     call add_vorticity_matrix_BC4(ZERO_FLUX,LAMBDA_TE)
  else
     call compute_vorticity_matrix4(ZERO_FLUX,BOHM_CURRENT)
     call add_vorticity_matrix_BC4(ZERO_FLUX,BOHM_CURRENT)
  end if

  call init_implicit_pastix(CSC_vort)
  call fill_implicit_pastix(CSC_vort,mat_vort_nnz,vorticity_mat)
  call free_all(vorticity_mat)
  call analyze_implicit_pastix(CSC_vort)
  !$OMP END MASTER
  !$OMP BARRIER

  !$OMP DO
  do k=1,global_parameters%N_Zones
     call MD_broadcast_corners_phi2(zones(k),STEP_OLD)
     call compute_vorticity(zones(k),STEP_OLD,ZERO_FLUX)
     call init_RHS_vorticity(zones(k))
#if VORTICITY_PI == 1
     call add_pi_time_plus_one_to_RHS(zones(k),ZERO_FLUX)
#endif
     call add_RHS_vorticity_sources(zones(k))
     if(drift_flags%solve_drift) then				!### Leybros modif ###
        call add_RHS_vorticity_drift(zones(k))
     end if
     call add_RHS_pe_Ohm_term(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP MASTER
  if((.not.flags%restart).or.(drift_flags%BC_type_phi.eq.1)) then
     call compute_vorticity_RHS4(ZERO_FLUX,LAMBDA_TE)
  else
     call compute_vorticity_RHS4(ZERO_FLUX,BOHM_CURRENT)
  end if
  call solve_implicit_pastix(CSC_vort)
  call redistribute_vorticity_sol(CSC_vort)
  call clear_pastix(CSC_vort)
  call free_pastix(CSC_vort)
  !$OMP END MASTER
  !$OMP BARRIER
  !$OMP DO
  do k=1,global_parameters%N_Zones
     call add_phi_Ohm_term(zones(k))
     call MD_broadcast_corners_phi2(zones(k),STEP_NEW)
     call compute_vorticity(zones(k),STEP_NEW,ZERO_FLUX)
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO
  do k=1,global_parameters%N_Zones
     call update_vorticity(zones(k))
     call store_pi_old(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER
  !  !$OMP MASTER
  !  call current_balance_custom2(zones(6),19)
  !  call current_balance_custom2(zones(4),10)
  !  !$OMP END MASTER
  !$OMP END PARALLEL

end subroutine solve_vorticity
