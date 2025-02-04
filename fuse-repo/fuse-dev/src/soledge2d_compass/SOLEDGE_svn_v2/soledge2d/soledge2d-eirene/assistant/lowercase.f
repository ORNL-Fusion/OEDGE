 
 
      subroutine EIRENE_lowercase (zeile)
 
      IMPLICIT NONE
      character(*), INTENT(INOUT) :: zeile
      INTEGER :: L, I, J, LANF, LEND
      character(26) :: klein, gross
      data klein /'abcdefghijklmnopqrstuvwxyz'/
      data gross /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
 
      LANF=verify(zeile,' ')
      LEND=verify(zeile,' ',.true.)
      if (lanf == 0) return
      l = lend-lanf+1
      zeile(1:l) = zeile(lanf:lend)
      if (lanf > 1) zeile(l+1:lanf-1) = repeat(' ',lanf-1)
      do i=1,l
        j=index(gross,zeile(i:i))
        if (j>0) zeile(i:i)=klein(j:j)
      end do
 
      return
      end subroutine EIRENE_lowercase
 
