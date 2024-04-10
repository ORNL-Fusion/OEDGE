module Mflux_surface
  
  type :: TSuper_tridiag
     real*8,allocatable :: a(:)
     real*8,allocatable :: b(:)
     real*8,allocatable :: c(:)
     real*8,allocatable :: S(:)
     real*8,allocatable :: buffer(:)
  end type TSuper_tridiag
  
  type :: TSuper_tridiag_propertites
     integer*4 :: i_psi
     integer*4 :: n_zones
     integer*4,allocatable :: zones(:)
     logical :: is_periodic
  end type TSuper_tridiag_propertites

  type :: TFlux_surface
     type(TSuper_tridiag) :: tridiag
     type(TSuper_tridiag_propertites) :: properties
     integer*4 :: Nz
     integer*4 :: ns_psi
  end type TFlux_surface
  
end module Mflux_surface
