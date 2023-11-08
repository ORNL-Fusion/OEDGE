module mod_rundiv_local
  
  implicit none

  ! jdemod
  !
  ! There are a number of local variables declared in rundiv.f that are
  ! passed thorugh various call chains so that they are accessible in various
  ! other routines. This approach is not ideal for the unstructured input code
  ! approach where these values are initialized and then may be updated from 
  ! a tagged input file. As a result, I am moving them to a module so that
  ! they can be used both in the existing call chains without conflicts and
  ! added to the updated unstructured input routines where their values
  ! can be updated without having to add this subset of inputs to a series
  ! of call chains.
  !
  ! From a design perspective these variables should have probably been 
  ! placed in a common block or other centrally accessible storage but when 
  ! the code was initially written that likely did not make sense. 
  !

  

      integer :: nymfs, nimps, nimps2, nizs, niters
      CHARACTER ::     TITLE*174, desc*1024, EQUIL*60
      real :: cpulim
               
  
contains



end module mod_rundiv_local
