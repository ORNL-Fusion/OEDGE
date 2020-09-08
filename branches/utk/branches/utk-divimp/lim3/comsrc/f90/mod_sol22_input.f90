module mod_sol22_input


  implicit none

  private
  
  integer,public:: nsol22_opt

  real,allocatable,public:: sol22_regions(:,:)
  character*100,allocatable,public :: sol22_filenames(:)
  


  public :: read_sol22_input


contains



  SUBROUTINE read_sol22_input
      use error_handling
      use mod_reader
      implicit none
      INTEGER   IERR,in,is
!
!  *********************************************************************
!  *                                                                   *
!  *  read_sol22_input: reads input for SOL22 regions                  *
!  *                                                                   *
!  *********************************************************************
!
     CHARACTER COMENT*128,MESAGE*72
!
      Ierr = 0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

      do in = 1,nsol22_opt

 100     IF (IBUF.EQ.0) READ (5,buff_format,ERR=9999,END=9999) BUFFER
         write(9,'(1x,a20,a1,a,1x,a6)') 'SOL22:',':',trim(buffer),'RDSOL22'
         IF (BUFFER(1:1).EQ.'$') GOTO 100

         IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
           CALL ReadUnstructuredInput(BUFFER)
           GOTO 100
         ENDIF

         MESAGE = 'EXPECTING COMMENT  4 reals and a character string'
         IBUF = 0
         READ (BUFFER,*,ERR=9999,END=9999) COMENT,(sol22_regions(in,is),is=1,4),sol22_filenames(in)

       end do
     
      RETURN

 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDSOL22 ERROR',mesage)
         call errmsg('RDSOL22 LINE ',trim(buffer))

      else
         call errmsg('RDSOL22 ERROR',mesage)
         call errmsg('RDSOL22 LINE ',trim(buffer))
      endif

      RETURN
    END SUBROUTINE read_sol22_input


  


end module mod_sol22_input
