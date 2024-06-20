 
 
 
      subroutine EIRENE_searchbin(xx,n,x,m)
c     ***********************************************************
c     * Ermittlung des Feldindex m mit vorgegebener Zahl x,     *
c     * so dass x zwischen xx(m) und xx(m+1)                    *
c     * Modifikation:   x <=xx(1) --> m=1  x >=xx(n) --> m=n    *
c     * *********************************************************
      USE EIRMOD_PRECISION
      implicit none
      integer n,m
      REAL(DP) :: x,xx(n)
      integer bl,bm,bu
 
      if (x.le.xx(1)) then
        m=1
      else if (x.ge.xx(n)) then
        m=n
      else
 
c  binary search
        bl=0
        bu=n+1
 
80      if (bu-bl.gt.1) then
          bm=(bu+bl)*0.5
          if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(bm))) then
            bl=bm
          else
            bu=bm
          endif
          goto 80
        endif
 
        m=bl
 
      endif
 
      return
      end
