C
C
      SUBROUTINE EIRENE_RDCN (ERSETZ,CONST)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
C
      CHARACTER(*), INTENT(IN) :: ERSETZ
      REAL(DP), INTENT(OUT) :: CONST
      INTEGER :: IEXPO, IW, IPUNKT, ILEN, ICON
      CHARACTER(10) :: FORM
C
      ILEN=INDEX(ERSETZ,'>')-2
      FORM='          '
C
C   INTEGER?
C
      IF (INDEX(ERSETZ,'.').EQ.0) THEN
        FORM(1:5)='(I  )'
        IF (ILEN.GT.9) WRITE (FORM(3:4),'(I2)') ILEN
        IF (ILEN.LE.9) WRITE (FORM(3:3),'(I1)') ILEN
        READ (ERSETZ(2:ILEN+1),FORM) ICON
        CONST=FLOAT(ICON)
C
      ELSE
C
C   REAL CONSTANT
C
        IPUNKT=INDEX(ERSETZ,'.')-1
        IEXPO=INDEX(ERSETZ,'E')
        IF (IEXPO.EQ.0) IEXPO=INDEX(ERSETZ,'D')
C
        IF (IEXPO.EQ.0) THEN
C   F-FORMAT
          FORM(1:2)='(F'
          IF (ILEN.GT.9) THEN
            WRITE (FORM(3:4),'(I2)') ILEN
            IW=4
          ELSE
            WRITE (FORM(3:3),'(I1)') ILEN
            IW=3
          ENDIF
          FORM(IW+1:IW+3)='. )'
          WRITE (FORM(IW+2:IW+2),'(I1)') ILEN-IPUNKT
          READ (ERSETZ(2:ILEN+1),FORM) CONST
C
        ELSE
C
C    E-FORMAT
          IEXPO=IEXPO-1
          FORM(1:2)='(E'
          IF (ILEN.GT.9) THEN
            WRITE (FORM(3:4),'(I2)') ILEN
            IW=4
          ELSE
            WRITE (FORM(3:3),'(I1)') ILEN
            IW=3
          ENDIF
          FORM(IW+1:IW+3)='. )'
          WRITE (FORM(IW+2:IW+2),'(I1)') IEXPO-(IPUNKT+1)
          READ (ERSETZ(2:ILEN+1),FORM) CONST
        ENDIF
      ENDIF
C
      RETURN
      END
