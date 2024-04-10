module Mphysics
  
  implicit none
  
  Type :: TAm_polynom
     real*8 :: Te_min,Te_max
     real*8 :: Ne_min,Ne_max
     real*8,allocatable :: coefficients(:)
  end type TAm_polynom

  Type :: TAmdata
     Type(TAm_polynom) :: recombination_rate_polynom
     Type(TAm_polynom) :: ionization_rate_polynom      
     Type(TAm_polynom) :: line_excitation_polynom
     Type(TAm_polynom) :: line_recombination_polynom !include also Brehmsstrahlung
     real*8 :: Ionization_potential
  end type TAmdata

  Type :: Tcxdata
    character(29):: cx_database1
    character(29):: cx_database2
    character(6) :: ibulk
    character(6) :: iscd1
    character(6) :: iscd2
    character(6) :: iscde
    real*8, dimension(9) :: fit_coeffs
  end type Tcxdata

  Type :: Telement
     character(2) :: symbol
     character(10) :: name
     integer*4 :: Z
     real*8 :: mass    ! set to zero for electrons
     real*8 :: mass2   ! set to proper value for electrons
     integer*4 :: amdata_polynom_degree
     Type(TAmdata),allocatable :: amdatas(:)
     integer*4 :: ncx
     Type(Tcxdata),allocatable :: cxdatas(:)
  end type Telement

  real*8, parameter :: kb = 1.3806d-23
  real*8, parameter :: ev = 1.6022d-19
  real*8, parameter :: pi = 3.14159265359
  real*8, parameter :: m_u = 1.6605d-27
  real*8, parameter :: m_e = 9.109d-31
  real*8, parameter :: epsilon_0 = 8.854187818d-12

end module Mphysics
