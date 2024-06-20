module Mturbulence

  implicit none

  type :: TImplicit_coefficients_ke
     real*8,allocatable :: west_k(:,:)
     real*8,allocatable :: east_k(:,:)
     real*8,allocatable :: west_epsilon(:,:)
     real*8,allocatable :: east_epsilon(:,:)
  end type TImplicit_coefficients_ke
  
  type :: Tkepsilon
     real*8,allocatable :: k(:,:)             ! turbulence intensity
     real*8,allocatable :: epsilon(:,:)       ! turbulence dissipation
     real*8,allocatable :: Sk(:,:)     
     real*8,allocatable :: Sepsilon(:,:)
     real*8,allocatable :: mu_t(:,:)          ! turbulence diffusivity
     real*8,allocatable :: interchange(:,:)   ! interchange grad(P)*grad(B)
     real*8,allocatable :: kh(:,:)            ! kelvin-helmotz grad(v)*grad(v)
     real*8,allocatable :: UE_shear(:,:)      ! ExB vtheta shear (dv_theta/dr)**2
     real*8 :: CornersK(2,2,2)                ! k corner values
     real*8 :: CornersEpsilon(2,2,2)          ! epsilon corner values
     type(TImplicit_coefficients_ke) :: implicit_coefs
  end type Tkepsilon

  type :: Tkepsilon_param
     real*8 :: Cmu, C1e, C2e, C3e
     real*8 :: sigma_k, sigma_epsilon
     real*8 :: sigma_T, sigma_n, sigma_v
     real*8 :: C_interchange
     real*8 :: C_kh
     real*8 :: th_interchange
     real*8 :: mu_max, mu_min
     real*8 :: gradT_weight
     real*8 :: deltaOmega
     real*8 :: tauV
  end type Tkepsilon_param
  
  
end module Mturbulence
