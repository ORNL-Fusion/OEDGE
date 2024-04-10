module Mreference_parameters
  
  implicit none

  type :: Tfields_rp
     real*8 :: n0 !reference density (/m3)
     real*8 :: c0 !reference velocity (m/s)
     real*8 :: T0 !reference temperature (K)
     real*8 :: T0eV !reference temperature (eV)
     real*8 :: tau0 !reference time (s)
     real*8 :: phi0 !reference potential (V)
     real*8 :: W0 !reference vorticity 
     real*8 :: k0 !reference k (m2s-2)
     real*8 :: epsilon0 !reference k (m2s-3)
  end type Tfields_rp

  type :: Tgeometry_rp
     real*8 :: R0 !reference major radius (m)
     real*8 :: Rm0 !reference minor radius (m)
     real*8 :: rs0 !reference domain radial thickness (m)
     real*8 :: Bpol0 !reference poloidal field (T)
     real*8 :: Btor0 !reference toroidal field (T)
     real*8 :: qref  !reference safety factor
     real*8 :: Score !Core surface (m2)
     real*8 :: A     !Aspect ratio * (2*pi)
  end type Tgeometry_rp

  type :: Treference_parameters
     type(Tfields_rp) :: fields
     type(Tgeometry_rp) :: geometry
  end type Treference_parameters


end module Mreference_parameters
