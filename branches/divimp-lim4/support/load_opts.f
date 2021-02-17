      subroutine load_options(numopt,opts,maxopts)
      implicit none
      integer numopt,maxopts
      character*(*) opts(maxopts) 
c
c     Locals 
c
      integer in,iargc
      external iargc
c
      numopt=iargc()
c
      if (numopt.le.0) then 
         numopt = 0
         return
      else
c
c        Do not load more than the maximum options
c
         numopt = min(numopt,maxopts)
c
         do in = 1,numopt
c
            call getarg(in,opts(in))
c
         end do
      endif
c
      return
c
      end  
