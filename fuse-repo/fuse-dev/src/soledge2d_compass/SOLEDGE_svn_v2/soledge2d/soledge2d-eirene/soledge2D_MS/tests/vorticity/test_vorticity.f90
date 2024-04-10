program test_vorticity
#include "compile_opt.inc"
  use test_var
  use all_variables
  use init_test_mod
  use Mdefinitions
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
  allocate(test_vort(1:global_parameters%N_Zones))
  
  global_variables%dt=1.d-5
  global_variables%dt_vort=global_variables%dt

  call init_vorticity()

  call compute_index()

  
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

  global_variables%dt=1.d-5
  global_variables%dt_vort=1.d-5

  !set parallel diffusion to zero
  do k=1,global_parameters%N_zones
     call compute_analytic_vorticity(k)
     do n=0,global_parameters%N_ions
        zones(k)%species(n)%transport_para%kappa=1.d-15
     end do
  end do

  !fill ghost cells (still periodic here)
  do k=1,global_parameters%N_zones
     call MD_broadcast_plasma(zones(k),STEP_OLD)
  end do
  do k=1,global_parameters%N_zones
     call MD_broadcast_corners(zones(k),STEP_OLD)
  end do

  !set plasma at time + 1 also
  do k=1,global_parameters%N_zones
     do n=0,global_parameters%N_ions
        zones(k)%species(n)%var(STEP_NEW)%density=zones(k)%species(n)%var(STEP_OLD)%density
        zones(k)%species(n)%var(STEP_NEW)%Gamma=zones(k)%species(n)%var(STEP_OLD)%Gamma
        zones(k)%species(n)%var(STEP_NEW)%temperature=zones(k)%species(n)%var(STEP_OLD)%temperature
     end do
  end do
  !fill ghost cells (still periodic here)
  do k=1,global_parameters%N_zones
     call MD_broadcast_plasma(zones(k),STEP_NEW)
  end do
  do k=1,global_parameters%N_zones
     call MD_broadcast_corners2(zones(k))
  end do

  !set non periodic
  zones(1)%Neighbors(N_East)=-3
  zones(1)%Neighbors(N_West)=-3
  zones(2)%Neighbors(N_East)=-3
  zones(2)%Neighbors(N_West)=-3
  zones(1)%MagNeighbors(N_East)=1
  zones(1)%MagNeighbors(N_West)=1
  zones(2)%MagNeighbors(N_East)=1
  zones(2)%MagNeighbors(N_West)=1
  do k=1,global_parameters%N_zones
     megazones(k)%is_periodic=.false.
     call store_pi_old(zones(k))
  end do
  
  !setting phi = Lambda*Te at time t
  call init_test_phi()
  do k=1,global_parameters%N_zones
     call compute_test_vorticity_sources(k)
  end do

  !now looking for phi at time t+1
  call solve_vorticity_test()

  call save_plasma_test(level)
  call save_vorticity_test(level) 

end program test_vorticity
