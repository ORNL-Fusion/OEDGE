c
c
      function ranset_eirene(ise)
      USE EIRMOD_PRECISION
      implicit none
      integer, intent(in) :: ise
 
      integer :: iseed
      common /cmem/ iseed
      real(dp) :: ranset_eirene, ranset_eirene_reinit
      integer, save :: ifirst=0
 
      if (ise <= 0) then
        if (ifirst == 0) then
         iseed = 9876543
!       else  seed schon gesetzt, uebernimmt dies
        end if
      else
         iseed = abs(ise)
      end if
 
      ifirst = 1
      ranset_eirene=0
      return
 
C     The following ENTRY is for reinitialization of EIRENE
 
      ENTRY ranset_eirene_reinit()
      ifirst = 0
      ranset_eirene_reinit = 0.D0
      return
      end
