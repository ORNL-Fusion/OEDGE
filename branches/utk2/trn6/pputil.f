**==PPCSTR
C
C=======================================================================
      SUBROUTINE PPCSTR(CSTR)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPCSTR
C
C PURPOSE : TO OUTPUT A CHARACTER STRING TO THE POST PROCESSOR FILE
C
C INPUT   : (C *) CSTR  - CHARACTER STRING FOR OUTPUT
C
C OUTPUT  : (C *) CSTR  - CHARACTER STRING FOR OUTPUT
C
C HISTORY : V1.R1.M0 --- 12/03/93 --- CREATION
C
C***********************************************************************
C
C DUMMY ARGUMENTS
C
      CHARACTER*(*) CSTR
C
      CHARACTER CBUFF*133
C
      INTEGER LPP,IFORM,LOUT
      COMMON/CPPOUT/LPP,IFORM
C
      PARAMETER (LOUT = 15)
C
      CBUFF = ' '
      CBUFF = CSTR
      IF(IFORM.EQ.0)THEN
        WRITE(LPP) CBUFF
      ELSE
        WRITE(LPP,'(A)') CBUFF
      ENDIF
C
 9999 RETURN
      END
**++EOF
**==PPRECW
C
C=======================================================================
      SUBROUTINE PPRECW(NAME,DESC,UNITS,RDATA,RSF,TYPE,NDATA,ISTAG)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPRECW
C
C PURPOSE : TO WRITE A DESCRIPTION HEADER + RECORD FOR A GIVEN DATA TYPE
C
C           NB..R*8 DATA IS CONVERTED TO R*4 BEFORE WRITING
C
C INPUT   : (C *) NAME   - DATA VARIABLE NAME
C           (C *) DESC   - DATA DESCRIPTION
C           (C *) UNITS  - DATA UNITS
C           ( - ) RDATA  - DATA (MAY BE I*4, R*4, R*8)
C           ( - ) RSF    - SCALE FACTOR (MAY BE I*4, R*4, R*8)
C           (C *) TYPE   - DATA VARIABLE TYPE ('I4','R4','R8')
C           (I*4) NDATA  - NUMBER OF DATA VALUES
C           (I*4) ISTAG  - STAGGERED VARIABLE FLAG (0-UNSTAGGERED
C                                                   1-STAGGERED)
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 12/03/93 --- CREATION
C           V1.R2.M0 --- 01/06/93 --- L.HORTON MODS: REMOVE MAP OPTION
C                                             & FIX EQUIVALENCING BUGS
C           V1.R3.M0 --- 14/12/94 --- ADD STAGGERED VARIABLE FLAG
C
C***********************************************************************
C
      CHARACTER*(*) NAME,DESC,UNITS,TYPE
      INTEGER       NDATA,ISTAG
      REAL*8        RDATA(NDATA),RSF
C
      INTEGER       LDATA
      PARAMETER    (LDATA =18000)
C
      CHARACTER*132 MSG1
      CHARACTER*80  CBUFF
      CHARACTER*32  FDESC
      CHARACTER*8   FNAME,FUNITS
      CHARACTER*2   FTYPE
      CHARACTER*1   CSTAG
      INTEGER       I
      INTEGER       I4DAT(LDATA),I4SF
      REAL*4        R4DAT(LDATA),R4SF
      REAL*8        R8DAT(LDATA),R8SF
C
      EQUIVALENCE(R8DAT(1),R4DAT(1),I4DAT(1))
      EQUIVALENCE(R8SF    ,R4SF    ,I4SF    )
C
      INTEGER LPP,IFORM,LOUT
      COMMON/CPPOUT/LPP,IFORM
C
      PARAMETER (LOUT = 15)
C
C
C INITIALISE
C ----------
C
      FNAME  = NAME
      FDESC  = DESC
      FUNITS = UNITS
      FTYPE  = TYPE
C
C CHECK INPUT VARIABLES
C ---------------------
C
      IF(FTYPE.NE.'I4' .AND. FTYPE.NE.'R4' .AND. FTYPE.NE.'R8')THEN
        CALL ERRMSS(LOUT,'PPRECW',1,'INVALID TYPE '//FTYPE//
     &             'VARIABLE '//FNAME,' ',' ')
      ENDIF
C
      IF(NDATA.GT.LDATA)THEN
        WRITE(MSG1,*) 'DIMENSION LDATA TOO SMALL'
     &        //' FOR VARIABLE '//FNAME//'  - NDATA =',NDATA
        CALL ERRMSS(LOUT,'PPRECW',1,MSG1,' ',' ')
      ENDIF
C
      IF(ISTAG.EQ.0)THEN
        CSTAG = ' '
      ELSE IF(ISTAG.EQ.1)THEN
        CSTAG = 'S'
      ELSE
        WRITE(MSG1,*) 'INVALID ISTAG VALUE'
     &             //' FOR VARIABLE '//FNAME//'  - ISTAG =',ISTAG
        CALL ERRMSS(LOUT,'PPRECW',1,MSG1,' ',' ')
      ENDIF
C
C SET UP LOCAL 8 BYTE DATA ARRAY
C ------------------------------
C
      R8SF = RSF
      DO 100 I=1,NDATA
        R8DAT(I)  = RDATA(I)
 100  CONTINUE
C
C SCALE THE DATA
C --------------
C
      IF(FTYPE.EQ.'R8')THEN
        DO 110 I=1,NDATA
          R4DAT(I) = R8SF*R8DAT(I)
 110    CONTINUE
      ELSE IF(FTYPE.EQ.'R4')THEN
        DO 120 I=1,NDATA
          R4DAT(I) = R4SF*R4DAT(I)
 120    CONTINUE
      ELSE IF(FTYPE.EQ.'I4')THEN
        DO 130 I=1,NDATA
          I4DAT(I) = I4SF*I4DAT(I)
 130    CONTINUE
      ENDIF
C
C WRITE HEADER AND RECORD
C -----------------------
C
      IF(NDATA.GT.0)THEN
C
        IF(IFORM.EQ.0)THEN
C
          CBUFF  = ' '
          WRITE(CBUFF,'(3A,I8,2A,1X,A)')
     &          '#',FNAME,FTYPE(1:1),NDATA,FDESC,FUNITS,CSTAG
          WRITE(LPP) CBUFF
          WRITE(LPP) (R4DAT(I),I=1,NDATA)
C
        ELSE IF (IFORM.EQ.1)THEN
C
          WRITE(LPP,'(3(A,1X),I8,1X,3(A,1X))')
     &              '#',FNAME,FTYPE(1:1),NDATA,FDESC,FUNITS,CSTAG
          WRITE(LPP,'(10Z8)') (R4DAT(I),I=1,NDATA)
C
        ELSE IF (IFORM.EQ.2)THEN
C
          WRITE(LPP,'(3(A,1X),I8,1X,3(A,1X))')
     &              '#',FNAME,FTYPE(1:1),NDATA,FDESC,FUNITS,CSTAG
          IF(FTYPE(1:1).EQ.'R')THEN
            WRITE(LPP,'(7(1PE11.3))') (R4DAT(I),I=1,NDATA)
          ELSE
            WRITE(LPP,'(7I10)') (I4DAT(I),I=1,NDATA)
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C-----------------------------------------------------------------------
C
 9999 RETURN
      END
**++EOF
