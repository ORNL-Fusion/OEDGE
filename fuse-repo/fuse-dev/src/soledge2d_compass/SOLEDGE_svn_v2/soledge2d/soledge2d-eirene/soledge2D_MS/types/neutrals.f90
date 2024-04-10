module Mneutrals

  implicit none

  Type :: Tneutrals
     real*8,allocatable :: density(:,:)
     real*8,allocatable :: Dn(:,:)
     real*8,allocatable :: RHS(:,:)
     real*8,allocatable :: Sn_nn(:,:)
  end type Tneutrals
  
  real*8 :: Rn_fluid
 
end module Mneutrals
