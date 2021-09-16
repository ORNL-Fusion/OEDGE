
c
      FUNCTION SECOND_OWN()
      USE PRECISION
      implicit none
      real(dp) x05baf,start,reset_second,second_own,time
      save start
      data start /0.0/
c x05baf is a NAG routine
      time=x05baf()
!pb      call cpu_time(time)
      second_own=time-start
      RETURN
C
      ENTRY RESET_SECOND()
      start=x05baf()
!pb      call cpu_time(start)
      reset_second=start
      return
      END
