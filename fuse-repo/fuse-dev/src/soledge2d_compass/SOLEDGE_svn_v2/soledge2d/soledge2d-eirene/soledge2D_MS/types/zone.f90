module Mzone

  use Mgeometry
  use Melectric_field
  use Mspecies
  use Mshared_fields
  use Mturbulence
  use Mneutrals
  use Mdrifts
  implicit none

  integer*4,parameter ::  N_North = 1
  integer*4,parameter ::  N_South = 2
  integer*4,parameter ::  N_East = 3
  integer*4,parameter ::  N_West = 4

  type :: TWeno
     real*8 :: C(0:1,0:1),Ctilde(0:1,0:1)
  end type TWeno

  type :: Tzone
     integer*4 :: number
     integer*4 :: Neighbors(4) ! North, South, East and West zone numbers
     integer*4 :: MagNeighbors(4) ! North, South, East and West zone numbers
     type(Tmesh) :: mesh
     type(Tmetric_coefficients) :: metric_coefficients
     type(Tjacobian) :: jacobian
     type(Tmasks) :: masks
     type(Telectric_field) :: electric_fields(2)
     type(Tkepsilon) :: kepsilon(2)
     type(Tshared_fields) :: shared_fields
     type(Tspecies),allocatable :: species(:)
     type(Tneutrals) :: neutrals
     type(TWeno),allocatable :: Weno(:,:)
     type(TDriftsExtraVars) :: DriftsExtra
  end type Tzone

  type :: Tmegazone
     integer*4 :: size ! number of zone in the megazone
     integer*4,allocatable :: zone_number(:)
     logical :: is_periodic
  end type Tmegazone

end module Mzone
