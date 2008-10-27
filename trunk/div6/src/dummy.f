c     -*-Fortran-*-
c
c
      REAL FUNCTION PsinToR(psin,z)
      IMPLICIT none
      REAL psin,z
      PsinToR = 0.0
      END


      SUBROUTINE SOL28_OLD(irs,ire,ikopt)
      INTEGER irs,ire,ikopt
      RETURN
      END
 
      SUBROUTINE SetBounds
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'
      INTEGER ik,ir
      DO ir = 1, nrs
        ikbound(ir,IKLO) = 1
        ikbound(ir,IKHI) = nks(ir)
      ENDDO
      RETURN
      END

c      SUBROUTINE BuildMap
c      CALL WN('BuildMap','Calling substitute routine')
c      END


      SUBROUTINE TailorGrid_Old
      implicit none
      CALL WN('BuildMap','Calling substitute routine')
      END

      SUBROUTINE ShapeTarget_Old(d1,d2)
      implicit none
      INCLUDE 'params'
      REAL*8 d1(MAXNKS,MAXNRS),d2(MAXNKS,MAXNRS)
      CALL WN('BuildMap','Calling substitute routine')
      END

      SUBROUTINE StructureGrid_Old(i1)
      implicit none
      integer i1
      CALL WN('StructureGrid','Calling substitute routine')
      END

c      SUBROUTINE UnstructureGrid
c      CALL WN('UnstructureGrid','Calling substitute routine')
c      END

c      SUBROUTINE InsertRing(i1,i2)
c      CALL WN('InsertRing','Calling substitute routine')
c      END

c      SUBROUTINE DeleteRing(i1)
c      CALL WN('DeleteRing','Calling substitute routine')
c      END

c      SUBROUTINE BalanceGrid
c      CALL WN('Balance','Calling substitute routine')
c      END

      SUBROUTINE GridSpace2
      implicit none
c      CALL WN('GridSpace2','Calling substitute routine')
      END

      SUBROUTINE EvaluateDensityPeak(i1,i2,i3)
      implicit none
      integer i1,i2,i3
      CALL WN('EvaluateDensityPeak','Calling substitute routine')
      END

      SUBROUTINE ChopSource(region,ir,source,multi)
      IMPLICIT none
      INCLUDE 'params'
      INTEGER region,ir
      REAL    source(MAXNKS,MAXNRS),multi
      CALL WN('ChopSource','Calling substitute routine')
      END

      SUBROUTINE ApplyPINQeMultiplier(i1,i2)
      implicit none
      integer i1,i2
      CALL WN('ApplyPINQeMultiplier','Calling substitute routine')
      END

      SUBROUTINE ReversePINQeMultiplier(i1,i2)
      implicit none
      integer i1,i2
      CALL WN('ReversePINQeMultiplier','Calling substitute routine')
      END

      REAL FUNCTION GetL1(i1,i2)
      implicit none
      integer i1,i2
      CALL WN('GetL1','Calling substitute routine')
      GetL1 = 0.0
      END

c      SUBROUTINE UpdateTargets(i1)
c      CALL WN('UpdateTargets','Calling substitute routine')
c      END

c      SUBROUTINE LoadPIN
c      END

       SUBROUTINE  AnalyseDensityPeakWidth(idum1)
       implicit none
       integer idum1
       END

      SUBROUTINE CalcLocalMockAdjustment2(ddum1,ddum2,idum1,
     .                                    rdum1)
      implicit none 
      INTEGER idum1
      REAL    rdum1(*)
      REAL*8  ddum1,ddum2
      END
