module Melectric_field

  implicit none
  
  type :: Telectric_field
     real*8,allocatable :: phi(:,:)             ! electric potential
     real*8,allocatable :: phi_smooth(:,:)             ! electric potential
     real*8,allocatable :: E(:,:)               ! electric field
     real*8,allocatable :: Etheta(:,:)               ! electric field along theta
     real*8,allocatable :: Epsi(:,:)               ! electric field along psi
     real*8,allocatable :: vorticity(:,:) 
     real*8,allocatable :: current(:,:,:) 
     real*8,allocatable :: pi_old(:,:)             ! main ion pressure (not updated at every time step)
     real*8,allocatable :: SW(:,:) ! vorticity source
     real*8,allocatable :: FluxW(:,:,:) ! vorticity source
     real*8,allocatable :: j_parallel(:,:,:) ! Ohm current
     real*8,allocatable :: j_para_adv_W(:,:,:) ! parallel current due to vorticity advection
     real*8,allocatable :: j_perp_adv_W(:,:,:) ! parallel current due to vorticity advection
     real*8,allocatable :: j_diff_W(:,:,:) ! perp current due to vorticity diffusion
     real*8,allocatable :: j_perp(:,:,:) ! perp current due to time variation
     real*8,allocatable :: RHS(:,:) ! vorticity source
     real*8 :: CornersW(2,2,2) ! vorticity corner values
     real*8 :: CornersPhi(2,2,2) ! phi corner values
     real*8 :: CornersPi_old(2,2,2) ! pi corner values
  end type Telectric_field

end module Melectric_field
