module test_var

  type :: Ttest_sources
     real*8,allocatable :: Sn(:,:)
     real*8,allocatable :: SG(:,:)
     real*8,allocatable :: ST(:,:)
     real*8,allocatable :: SW(:,:)
  end type Ttest_sources

  type :: Ttest_metric
     real*8,allocatable :: Jac(:,:)
     real*8,allocatable :: cpp(:,:)
     real*8,allocatable :: ctt(:,:)
     real*8,allocatable :: cpt(:,:)
     real*8,allocatable :: G(:,:)
  end type Ttest_metric

  type :: Ttest_vorticity
     real*8,allocatable :: W(:,:)
     real*8,allocatable :: phi(:,:)
     real*8,allocatable :: r(:,:)
     real*8,allocatable :: theta(:,:)
  end type Ttest_vorticity
  
  real*8,parameter :: test_R0=2.D0
  real*8,parameter :: test_a=0.8D0
  real*8,parameter :: test_psi0=0.4D0
  real*8,parameter :: test_B0=3.D0 
  real*8,parameter :: test_n0=1.d19
  real*8,parameter :: test_Gamma0=1.d21
  real*8,parameter :: test_T0eV=100.
  real*8,parameter :: test_phi0=100.
  real*8,parameter :: test_rmin=0.48D0 ! test_a*0.6
  real*8,parameter :: test_rmax=1.12D0 ! test_a*1.4
  integer*4 :: test_configuration
  real*8 :: test_theta_shift
  real*8 :: test_reg_r
  real*8 :: test_reg_theta
  Type(Ttest_sources),allocatable :: test_sources(:,:)
  Type(Ttest_metric),allocatable :: test_metric(:)
  Type(Ttest_vorticity),allocatable :: test_vort(:)
  
  integer*4,parameter :: TEST_PARA=1
  integer*4,parameter :: TEST_PERP=2

end module test_var
