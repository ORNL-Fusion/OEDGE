

      SUBROUTINE EIRENE_BROAD_USR
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsou
      use styx2eirene
      use eirmod_cpes
      IMPLICIT NONE
      if (my_pe .ne. 0) then 
	if (.not.allocated(time_para)) allocate(time_para(nprs))
      	time_para=0
      endif
      

      RETURN
      END
