module data_structures

type fit_data
    ! bounds
    real*8 :: Nemin,Nemax,Temin,Temax
    ! degree of fit
    integer :: degree
    ! goodness of fit
    real*8 :: goodness
    ! fit parameters
    real*8, allocatable :: coeffs(:) 
    ! density and temperature arrays (for checks: ability to compare to original data ?)
    integer :: NNe,NTe
    real*8, allocatable :: Ne(:),Te(:) 
  end type fit_data

  type element
     type(fit_data), allocatable :: ionization(:),recombination(:),plt(:),prb(:)
     character(2) :: element
  end type element

  type(element) :: element_

  type :: rate
    real*8, allocatable :: values(:,:)
  end type rate

! number of columns in data files < 10

  integer, parameter :: ncol = 5 


end module data_structures
