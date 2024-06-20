C
C
      SUBROUTINE EIRENE_STCOOR (X,Y,IFL)
c  instor: zaehler, instor-ter aufruf
c  ifl: flag:=0, first point, .ne.0: else, npl2d(instor)=ifl
c  x,y: co-ordinaten des punktes, stored on xpl2d(instor),ypl2d(instor)
C
c  alle teilstuecke stehen hintereinander auf xpl2d,ypl2d
c  zum entwirren:NUMSUR(inums,..) array
c  inums: zaehler fuer individuelle zu plottende kurvenstuecke
c  jedes bekommt die eigene nummer und arrow
c  daher: ggfls: eine nummer mehrmals im plot vorhanden
c  NUMSUR : nummer des flachestuecks
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CRECH
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: X, Y
      INTEGER, INTENT(IN) :: IFL
      TYPE(PPOINT), POINTER :: POINT
      INTEGER, SAVE :: IFIRST=0
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        INUMS=0
        NULLIFY(FIRST_POINT)
        NULLIFY(LAST_POINT)
      ENDIF
C
      ALLOCATE(POINT)
      POINT%XPL2D = X
      POINT%YPL2D = Y
      POINT%NPL2D = IFL
      POINT%NUMSUR = INOSF
      NULLIFY(POINT%NXTPNT)
 
      IF (ASSOCIATED(LAST_POINT)) THEN
        LAST_POINT%NXTPNT => POINT
        LAST_POINT => POINT
      ELSE
        FIRST_POINT => POINT
        LAST_POINT => FIRST_POINT
      END IF
 
      INUMS = INUMS + 1
      INSTOR = INSTOR + 1
C
      RETURN
 
C     the following ENTRY is for reinitialization of EIRENE
 
      ENTRY EIRENE_STCOOR_REINIT
      IFIRST = 0
      return
 
      END
