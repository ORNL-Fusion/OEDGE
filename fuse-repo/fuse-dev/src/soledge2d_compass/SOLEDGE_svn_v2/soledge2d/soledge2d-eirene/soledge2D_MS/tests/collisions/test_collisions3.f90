program test_collisions3
#include "compile_opt.inc"
  use test_var
  use all_variables
  use init_test_mod
#if GFORTRAN == 0
  use ifport
#endif
  implicit none
  character(len=32) :: arg
  integer*4 :: level,n,k
  !#############################################################
  !the test needs the grid size as argument
  call getarg(1,arg)
  read(arg,'(I2)') level
  call getarg(2,arg)
  read(arg,'(I1)') test_configuration !1: regular ; 2: unregular
  !the grid size is given by Nx=Nz=2^(1+level)
  !#############################################################
  call init_test(level,'He',100.D0,2)

  !set parallel diffusion to zero
  do k=1,global_parameters%N_zones
     do n=0,global_parameters%N_ions
        zones(k)%species(n)%transport_para%kappa=1.d-15
     end do
  end do
 
  call compute_test_collision_sources(3)

  call step_in_soledge_collisions3()

  call save_plasma_test(level) 

end program test_collisions3
