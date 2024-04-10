program test_advection
#include "compile_opt.inc"
  use test_var
  use all_variables
  use init_test_mod
#if GFORTRAN == 0
  use ifport
#endif
  implicit none
  character(len=32) :: arg
  integer*4 :: level
  !#############################################################
  !the test needs the grid size as argument
  call getarg(1,arg)
  read(arg,'(I2)') level
  call getarg(2,arg)
  read(arg,'(I1)') test_configuration !1: regular ; 2: unregular ; 3: metric
  !the grid size is given by Nx=Nz=2^(1+level)
  !#############################################################
  call init_test(level,'D ',100.D0,1)
 
  allocate(test_metric(1:global_parameters%N_zones))
  call compute_test_metric()

  call save_metric_test(level)

end program test_advection
