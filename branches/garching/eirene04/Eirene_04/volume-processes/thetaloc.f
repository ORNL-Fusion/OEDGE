
      

      subroutine thetaloc(xx,a,x,b)
c     ***********************************************************
c     * Ermittlung des Feldindex b mit vorgegebener RAN-Zahl x, *
c     * so dass x zwischen RAN(b) und RAN(b+1)                  *
c     * Modifikation: für x=RAN(bmax) erhält man bmax           *
c     * *********************************************************
      USE PRECISION
      implicit none
      integer a,b
      REAL(DP) :: x,xx(a)
      integer bl,bm,bu

      bl=0
      bu=a+1
      
90    if(bu-bl.gt.1)then
        bm=(bu+bl)*0.5
        if((xx(a).ge.xx(1)).eqv.(x.ge.xx(bm)))then
          bl=bm
        else
          bu=bm
        endif
      goto 90
      endif

      if(x.eq.xx(1))then
        b=1
      else if(x.eq.xx(a))then
        b=a
      else
        b=bl
      endif

      return
      end
