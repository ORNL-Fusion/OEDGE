 
c
      FUNCTION EIRENE_SECOND_OWN()
      USE EIRMOD_PRECISION
      implicit none
!pb      real(dp) x05baf,start,reset_second,second_own,time
      real(dp) ::  EIRENE_x05baf, EIRENE_reset_second, EIRENE_second_own
      real(sp) :: start,time
      save start
      data start /0.0/
c x05baf is a NAG routine
!pb      time= EIRENE_x05baf()
      call cpu_time(time)
      EIRENE_second_own=time-start
      RETURN
C
      ENTRY EIRENE_RESET_SECOND()
!pb      start=x05baf()
      call cpu_time(start)
      EIRENE_reset_second=start
      return
      END
