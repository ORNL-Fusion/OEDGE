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
  integer*4 :: k,n
  !#############################################################
  !the test needs the grid size as argument
  call getarg(1,arg)
  read(arg,'(I2)') level
  call getarg(2,arg)
  read(arg,'(I1)') test_configuration !1: regular ; 2: unregular
  !the grid size is given by Nx=Nz=2^(1+level)
  !#############################################################
  call init_test(level,'D ',100.D0,1)

  !set diffusion to zero
  transport_parameters%Dn_p(1)=0.D0
  transport_parameters%Dn_t(1)=0.D0
  transport_parameters%nu_p(1)=0.D0
  transport_parameters%nu_t(1)=0.D0
  transport_parameters%chie_p=0.D0
  transport_parameters%chie_t=0.D0
  transport_parameters%chii_p(1)=0.D0
  transport_parameters%chii_t(1)=0.D0
  transport_parameters%v_pinch(1)=0.D0

  !set parallel diffusion to zero
  do k=1,global_parameters%N_zones
     do n=0,global_parameters%N_ions
        zones(k)%species(n)%transport_para%kappa=1.d-15
     end do
  end do

  call compute_test_advection_sources()

  call step_in_soledge_advection()

  call save_plasma_test(level) 

end program test_advection
