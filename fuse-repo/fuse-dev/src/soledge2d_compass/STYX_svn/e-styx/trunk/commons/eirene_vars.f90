module Meirene_vars
  
  Type :: TEirene_vars
     integer*4 :: feedback
  end type TEirene_vars

  type volume
     real*8, allocatable :: cell(:,:)
  end type volume

  Type(TEirene_vars) :: eirene_vars
  integer*4, parameter :: Ntrimax = 6

end module Meirene_vars
