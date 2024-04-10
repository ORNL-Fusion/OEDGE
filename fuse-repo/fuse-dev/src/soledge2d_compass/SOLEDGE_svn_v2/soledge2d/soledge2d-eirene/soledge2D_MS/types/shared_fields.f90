module Mshared_fields
  
  implicit none

  type :: Tshared_fields
     real*8,allocatable :: log_Lambda(:,:)
     real*8,allocatable :: Zeff(:,:)
  end type Tshared_fields

end module Mshared_fields
