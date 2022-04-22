module mod_sol22_input_lim


  implicit none

  private
  
  integer,public:: sol22_opt
  integer,public:: nsol22_opt

  real,allocatable,public:: sol22_regions(:,:)

  character*500,allocatable,public :: sol22_filenames(:)
  character*500,public :: sol22_default_filename
  
  public :: read_sol22_input


contains



  SUBROUTINE read_sol22_input(nsol)
      use error_handling
      use allocate_arrays
      use mod_reader
      implicit none
!
!  *********************************************************************
!  *                                                                   *
!  *  read_sol22_input: reads input for SOL22 regions                  *
!  *                                                                   *
!  *********************************************************************
!
!
      INTEGER   IERR,in,is,nsol
      CHARACTER COMENT*128,MESAGE*72

      if (nsol.le.0) return

      write(0,*) 'read_sol22_input:',nsol,nsol22_opt
      
      call allocate_array(sol22_regions,nsol,5,'SOL22 Region data',ierr)
      if (allocated(sol22_filenames)) deallocate(sol22_filenames)
      allocate(sol22_filenames(nsol))

      write(0,*) 'read_sol22_input:allocated',nsol


      Ierr = 0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

      do in = 1,nsol

 100     IF (IBUF.EQ.0) READ (5,buff_format,ERR=9999,END=9999) BUFFER
         write(9,'(1x,a20,a1,a,1x,a6)') 'SOL22:',':',trim(buffer),'RDSOL22'
         write(0,'(1x,a20,a1,a,1x,a6)') 'SOL22:',':',trim(buffer),'RDSOL22'
         IF (BUFFER(1:1).EQ.'$') GOTO 100

         IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
           CALL ReadUnstructuredInput(BUFFER)
           GOTO 100
         ENDIF

         MESAGE = 'EXPECTING COMMENT  5 reals and a character string'
         IBUF = 0
         READ (BUFFER,*,ERR=9999,END=9999) COMENT,(sol22_regions(in,is),is=1,5),sol22_filenames(in)

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


  


end module mod_sol22_input_lim
