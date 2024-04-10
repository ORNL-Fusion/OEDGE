program test_am
#include "compile_opt.inc"
  use test_var
  use all_variables
  use init_test_mod
#if GFORTRAN == 0
  use ifport
#endif
  implicit none
  character(len=32) :: arg
  character(len=2) :: element
  integer*4 :: k,n,level,n_ite
  real*8 :: temperature
  integer*4 :: n_iterations
  !#############################################################
  !the test needs the grid size as argument
  call getarg(1,arg)
  read(arg,'(A2)') element
  call getarg(2,arg)
  read(arg,'(ES15.7)') temperature
  call getarg(3,arg)
  read(arg,'(I10)') n_iterations
  test_configuration = 1
  !1: regular ; 2: unregular
  level = 2
  !the grid size is given by Nx=Nz=2^(1+level)
  !#############################################################
  call init_test(level,element,temperature,3)

  write(*,*) 'Reference temperature = ',reference_parameters%fields%T0eV,' eV'

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

  global_parameters%n_iterations=n_iterations

  n_ite=0
  zones(1)%species(0)%residuals%resn=1.

  do while ((zones(1)%species(0)%residuals%resn.gt.1.d-10).and.(n_ite.lt.global_parameters%n_iterations))
     call step_in_soledge_am()
     call update_all_fields()
     call save_residuals(n_ite)   
     n_ite=n_ite+1
  end do

  call save_plasma_test(level) 
  call save_residuals_force()

end program test_am
