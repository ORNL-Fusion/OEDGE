 
      subroutine EIRENE_remove_char (inchar, outchar ,remchar)
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comprt, only : iunout
 
      implicit none
 
      character(len=*), intent(in) :: inchar, remchar
      character(len=*), intent(out) :: outchar
      integer :: lin, lout, i, io
 
      lout = len(outchar)
      outchar = repeat(' ',lout)
 
      lin=len_trim(inchar)
 
      io=0
 
      do i=1,lin
        if (index(remchar,inchar(i:i)) == 0) then
           io = io + 1
           if (io > lout) then
              write (iunout,*) ' ERROR IN REMOVE_CHAR '
              write (iunout,*) ' OUTPUT STRING IS TOO SHORT'
              write (iunout,*) ' INCHAR = ',inchar
              write (iunout,*) ' REMCHAR = ',remchar
              write (iunout,*) ' STRING SHORTENED TO ',outchar
              return
           end if
           outchar(io:io) = inchar(i:i)
        end if
      end do
 
      return
      end subroutine EIRENE_remove_char
