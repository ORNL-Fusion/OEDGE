c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE dumm1
      IMPLICIT none
      END
c
c ======================================================================
c
c      SUBROUTINE BuildInversionMesh(idum1)
c      IMPLICIT none
c      INTEGER idum1
c      RETURN
c      END
c
c ======================================================================
c
      LOGICAL FUNCTION CheckInversionCell(idum1,rdum1,rdum2)
      IMPLICIT none
      INTEGER idum1
      REAL    rdum1,rdum2
      CheckInversionCell = .FALSE.
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE LoadVesselStructures(idum1)
      IMPLICIT none
      INTEGER idum1
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE ProcessTriangleGrid(idum1)
      IMPLICIT none
      INTEGER idum1
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE ProcessMagneticGrid(idum1)
      IMPLICIT none
      INTEGER idum1
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE ProcessTetrahedronGrid(idum1)
      IMPLICIT none
      INTEGER idum1
      RETURN
      END
c
c ======================================================================
c
c      INTEGER FUNCTION GetNumberOfObjects()
c      IMPLICIT none
c      INTEGER idum1
c      GetNumberOfObjects = -1
c      RETURN
c      END
c
c ======================================================================
c
      SUBROUTINE WN(routine,message)
      IMPLICIT none
      CHARACTER routine*(*),message*(*)
      WRITE(0,'(4A)') ' WARNING ',routine,': ',message
      RETURN 
      END
c
c ======================================================================
c
      SUBROUTINE ER(routine,message,*)
      IMPLICIT none
      CHARACTER routine*(*),message*(*)
      WRITE(0,'(4A)') ' ERROR ',routine,': ',message
      RETURN 1
      END
c
c ======================================================================
c
      SUBROUTINE GetFluidGridEmission(idum1,idum2,idum3,rdum1,rdum2)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none
      INTEGER idum1,idum2,idum3
      REAL    rdum1(idum2,idum3),rdum2
      RETURN
      END
c
c ======================================================================
c
c      RECURSIVE SUBROUTINE LoadTriangleData(flag1,flag2,flag3,normalize,
c     .                                      tdata)
c      IMPLICIT none
c      INTEGER flag1,flag2,flag3,normalize
c      REAL    tdata(*)
c 
c      WRITE(0,*) 'CALL TO LoadTriangleData NOT ALLOWED'
c      STOP
c      END
c
c ====================================================================== 
c
      SUBROUTINE AssignPlasmaQuantities(idum1,idum2,idum3)
      IMPLICIT none
      INTEGER idum1,idum2,idum3
      END
c
c ======================================================================
c
      SUBROUTINE PlotLineShapes(istart,iend,MAXNPIXEL,npixel,pixel)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER istart,iend,MAXNPIXEL,npixel
      TYPE(type_view) :: pixel(MAXNPIXEL)

      RETURN
      END
c
c     ======================================================================
c
      SUBROUTINE NextLine(fp,ntally,icount,rdum,binary)
      IMPLICIT none

      INTEGER   fp,ntally,icount
      LOGICAL   binary
      REAL      rdum(*)   

      WRITE(0,*) 'SHOULD NOT BE IN NextLine'
      STOP
      
      END
c     
c ======================================================================
c
      SUBROUTINE Wrapper_ClearObjects
      IMPLICIT none

      WRITE(0,*) 'SHOULD NOT BE IN Wrapper_ClearObjects'
      STOP
      
      END
c     
c ======================================================================
c     
      SUBROUTINE rayGenerateRibbonGrid
      IMPLICIT none

      WRITE(0,*) 'SHOULD NOT BE IN rayGenerateRibbonGrid'
      STOP
      
      END
c     
c ======================================================================
c     
      SUBROUTINE DumpShinKajita(title)
      IMPLICIT none

      CHARACTER*(*) title
      
      RETURN      
      END
c     
c ======================================================================
c     
