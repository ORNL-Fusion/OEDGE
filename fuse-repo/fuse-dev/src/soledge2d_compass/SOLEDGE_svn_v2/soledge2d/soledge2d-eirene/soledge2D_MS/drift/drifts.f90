module Mdrifts

  implicit none
    
  type :: Tdrifts
     real*8,allocatable :: uEp(:,:), uEt(:,:)           ! electric drift velocity in psi and theta
     real*8,allocatable :: udp(:,:), udt(:,:)           ! diamagnetic drift velocity in psi and theta
     real*8,allocatable :: uBp(:,:), uBt(:,:)		! curvature drift in psi and theta
     real*8,allocatable :: jdiam(:,:,:), jExB(:,:,:), jBxDB(:,:,:)		! Current for drift term in vorticity
     real*8,allocatable :: uEps(:,:), uEts(:,:)           ! electric drift velocity in psi and theta save
     real*8,allocatable :: uBps(:,:), uBts(:,:)		! curvature drift in psi and theta save
  end type Tdrifts

  type :: TDriftsExtraVars
     logical,allocatable :: OK_points(:,:)
     integer*4,allocatable :: ref_points(:,:,:,:) !Nx*Nz*3(i,j,k)*Nref_points_max(4)
     integer*4,allocatable :: Nref_points(:,:) !Nx*Nz
  end type TDriftsExtraVars

end module Mdrifts
