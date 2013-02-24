module sol22_input

!
! jdemod - I am creating this data module to start replacing common blocks.
!
! Modules perform the same task but allow for both data sharing and initialization.
! In addition, executable code related to variable setup can also be included. 
!
! Ultimately it would be good for all the unstructured input to be organized and added to files 
! like this one. 
!

  integer,public :: debug_sol22 = .false.   ! 284 - SOL22 debug flag - default 0 = off
  integer,public :: debug_sol22_ir = 1      ! 285 - SOL22 debug - ring for hi res plasma - default = 1
  integer,public ::  debug_sol22_ikopt = 1  ! 286 - SOL22 debug - ikopt (ring end) for detailed plasma - default = 1

contains



end module sol22_input
