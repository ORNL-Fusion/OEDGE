module Mgeometry
  
  implicit none

  type :: Tmasks
     real*8,allocatable :: chi(:,:)
     real*8,allocatable :: chi1(:,:)
     real*8,allocatable :: chi2(:,:)
     real*8,allocatable :: chi3(:,:)
     real*8,allocatable :: chi4(:,:)
     real*8,allocatable :: chi5(:,:)
     real*8,allocatable :: chi6(:,:)
     integer*4,allocatable :: npts_around_penwall(:,:)
     integer*4,allocatable :: pts_around_penwall(:,:,:,:)
  end type Tmasks
  
  type :: Tmesh
     integer*4 :: Nx
     integer*4 :: Nz
     real*8 :: xmin,xmax
     real*8 :: zmin,zmax
     real*8,allocatable :: z(:,:)
     real*8,allocatable :: x(:,:)
     real*8,allocatable :: z_plus_1half(:,:)
     real*8,allocatable :: z_minus_1half(:,:)
     real*8,allocatable :: x_plus_1half(:,:)
     real*8,allocatable :: x_minus_1half(:,:)
     real*8,allocatable :: Rgeom(:,:)
     real*8,allocatable :: Zgeom(:,:)
     real*8,allocatable :: Rcorner(:,:)
     real*8,allocatable :: Zcorner(:,:)
     real*8,allocatable :: Br(:,:)
     real*8,allocatable :: Bz(:,:)
     real*8,allocatable :: Bphi(:,:)
     real*8,allocatable :: B(:,:)
     real*8 :: cornerx(2,2,2)
     real*8 :: cornerz(2,2,2)
     real*8 :: cornerRg(2,2,2)
     real*8 :: cornerZg(2,2,2)
     real*8 :: cornersB(2,2,2)
     integer*4,allocatable :: index(:,:)
  end type Tmesh

  type :: Tmetric_coefficients
     real*8,allocatable :: cpp(:,:),ctt(:,:),cpt(:,:)
     real*8,allocatable :: c_pp(:,:),c_tt(:,:),c_pt(:,:)
     real*8,allocatable :: Jacobian(:,:),G(:,:)
     real*8,allocatable :: dS_north_PU(:,:),dS_south_PU(:,:)
     real*8,allocatable :: dS_west_PU(:,:),dS_east_PU(:,:)
     real*8,allocatable :: dS_north_DD(:,:),dS_south_DD(:,:)
     real*8,allocatable :: dS_west_DD(:,:),dS_east_DD(:,:)
     real*8,allocatable :: dvol_PU(:,:),dvol_DD(:,:)
     real*8,allocatable :: sinepitch_east(:,:),sinepitch_west(:,:)
     real*8,allocatable :: divb(:,:)
     real*8,allocatable :: delta_r(:,:)
  end type Tmetric_coefficients

  type :: Tjacobian ! T stands for theta, P stands for psi
     real*8,allocatable :: dTdR(:,:)
     real*8,allocatable :: dTdZ(:,:)
     real*8,allocatable :: dPdR(:,:)
     real*8,allocatable :: dPdZ(:,:)
     real*8,allocatable :: dRdT(:,:)
     real*8,allocatable :: dRdP(:,:)
     real*8,allocatable :: dZdT(:,:)
     real*8,allocatable :: dZdP(:,:)
  end type Tjacobian

end module Mgeometry
