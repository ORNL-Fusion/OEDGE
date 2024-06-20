module Mplots

    implicit none

  ! This module contains the variable for each specie
  ! Specie number 0 is reserved for electrons

  type :: Tplot
     integer*4 :: type
     integer*4 :: nzones
     integer*4,allocatable :: zones(:)
     integer*4 :: coord
     integer*4 :: coord2
  end type Tplot

end module Mplots
