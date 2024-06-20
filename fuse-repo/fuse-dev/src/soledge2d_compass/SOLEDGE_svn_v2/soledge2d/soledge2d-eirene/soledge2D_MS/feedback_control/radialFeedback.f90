module MradialFeedback

  implicit none
  
  Type :: RFDset
     real*8, allocatable :: fluxN(:)
     real*8, allocatable :: fluxTe(:)
     real*8, allocatable :: fluxTi(:)
  end type RFDset

  Type :: TradialFeedbackData
     integer*4 :: Nzones
     integer*4,allocatable :: zone_profile(:)
     integer*4 :: Nz
     integer*4 :: Nxtot
     real*8,allocatable :: input_gradN(:)
     real*8,allocatable :: input_gradTe(:)
     real*8,allocatable :: input_gradTi(:)
     real*8,allocatable :: input_N(:)
     real*8,allocatable :: input_Te(:)
     real*8,allocatable :: input_Ti(:)
     real*8,allocatable :: input_D(:)
     real*8,allocatable :: input_Chi(:)
     real*8,allocatable :: x(:)
     real*8,allocatable :: D(:)
     real*8,allocatable :: chie(:)
     real*8,allocatable :: chii(:)
     Type(RFDset),allocatable :: set(:)
     real*8 :: Dmin
     real*8 :: Dmax
     real*8 :: keep
     real*8 :: Gain
     real*8 :: GainG
  end type TradialFeedbackData

  Type(TradialFeedbackData) :: radialFeedbackData

end module MradialFeedback
