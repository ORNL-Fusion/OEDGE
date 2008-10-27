


c     ***********************************************************
c     --Subroutinen--


      subroutine energyloc(xx,a,x,b)
c     ***********************************************************
c     * Ermittlung des Feldindex b mit vogegebener Energie x,   *
c     * so dass x zwischen Energie(b) und Energie(b+1)          *
c     * Modifikation: Ergebnis: Energien zw. 0.1 eV und 100 eV  *
c     * *********************************************************
      USE PRECISION
      implicit none
      integer a,b
      REAL(DP) :: x,xx(a)
      integer bl,bm,bu

      bl=0
      bu=a+1
      
80    if(bu-bl.gt.1)then
        bm=(bu+bl)*0.5
        if((xx(a).ge.xx(1)).eqv.(x.ge.xx(bm)))then
          bl=bm
        else
          bu=bm
        endif
      goto 80
      endif

      if(x.le.xx(1))then
        b=1
      else if(x.eq.xx(a))then
        b=a
      else
        b=bl
      endif

      return
      end
