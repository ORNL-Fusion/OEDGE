C ===== SOURCE: locstr_usr.f
csw
csw routine to find a specific string in unit fp
csw mar2006
csw s.wiesen@fz-juelich.de
csw
      subroutine eirene_locstr_usr(fp,sstr,ier)
      implicit none
      integer, intent(in) :: fp
      integer, intent(out) :: ier
      character(*), intent(in) :: sstr
      character(256) :: line
      integer :: ierror

      ier=1
      rewind(fp)
      do
         read(fp,'(a256)',iostat=ierror) line
         if(ierror /= 0) return

         if( index(line,sstr) == 1) then
            ier = 0
            return
         endif
      enddo
      end
