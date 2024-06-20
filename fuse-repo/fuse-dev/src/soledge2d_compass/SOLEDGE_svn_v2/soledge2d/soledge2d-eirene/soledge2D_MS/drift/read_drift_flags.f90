subroutine read_drift_flags
  use all_variables, only : drift_flags, transport_parameters, global_variables
  use Mreaders
  use Msmoothing_vars
  implicit none
  integer :: finput
  finput=100
  open(unit=finput,file='input_drifts.txt',status='unknown')
  call skip_line(finput,5)
  read(finput,1) drift_flags%solve_phi
  call skip_line(finput,2)
  read(finput,2) drift_flags%solve_drift
  call skip_line(finput,2)
  read(finput,3) drift_flags%N_solve_phi
  call skip_line(finput,2)
  read(finput,2) drift_flags%vorticity_restart
  call skip_line(finput,6)
  read(finput,4) drift_flags%use_ExB_radial
  read(finput,1) drift_flags%use_ExB_poloidal
  call skip_line(finput,2)
  read(finput,4) drift_flags%use_gradB_radial
  read(finput,1) drift_flags%use_gradB_poloidal
  call skip_line(finput,5)
  read(finput,5) drift_flags%jdiam
  read(finput,5) drift_flags%jadvW
  read(finput,5) drift_flags%jdiffW
  read(finput,5) drift_flags%jgradBW
  read(finput,5) drift_flags%jExBW
  call skip_line(finput,6)
  read(finput,6) drift_flags%BC_type_phi
  call skip_line(finput,2)
  read(finput,7) transport_parameters%zeta_p
  read(finput,8) transport_parameters%zeta_t
  call skip_line(finput,2)
  read(finput,9)  global_variables%eta_para_smooth
  smoothing_diffusivity=global_variables%eta_para_smooth
  global_variables%eta_para_smooth=0.D0
  call skip_line(finput,2)
  read(finput,10)  drift_flags%eta_perp0
  call skip_line(finput,2)
  read(finput,7)  drift_flags%Mach_lim
  call skip_line(finput,2)
  read(finput,7)  drift_flags%Mach_lim_rad
  call skip_line(finput,2)
  read(finput,11) drift_flags%reverse_B
  close(finput)
1 format(11X,L1)
2 format(14X,L1)
3 format(8X,I6)
4 format(9X,L1)
5 format(15X,L1)
6 format(14X,I1)
7 format(11X,F5.2)
8 format(13X,F5.2)
9 format(14X,F7.4)
10 format(12X,F7.4)
11 format(12X,L1)
end subroutine read_drift_flags

