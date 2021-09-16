c     -*Fortran*-
c
c**==PPCSTR
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
c
c     This parameter has to be larger than 5 * the maximum polygons
c
      PARAMETER    (LDATA =40000)
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
c     Reset output channel to be the LIM file so that error 
c     messages are sent there. 
c
c      PARAMETER (LOUT = 15)
c
      PARAMETER (LOUT = 6)
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
c
c          write(0,'(3a,i6)') 'PPRECW:',fname,fdesc,ndata
c          write(6,'(3a,i6)') 'PPRECW:',fname,fdesc,ndata
c          write(6,'(10(1x,e12.5))') (r4dat(i),i=1,ndata)    
c          do i = 1,ndata
c             write(6,'(a,i5,3(1x,g12.5))') fname,i,r4dat(i)
c          end do 
c
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

c
c
c     NOTE: The following are taken from EDGE2D - SUPPORTZ.F
c
c
C
C=======================================================================
      SUBROUTINE CHRLTU( STRIN , STROUT )
      IMPLICIT NONE
C
C++ ....................................................................
C
C VERSION : V1.R1.M0
C
C ROUTINE : CHARACTER LOWER TO UPPER CONVERSION
C           -- -      -     -- -
C
C PURPOSE : TO CONVERT LOWER CASE CHARACTERS IN A STRING TO UPPER CASE.
C
C INPUT   : (C**) STRIN   = INPUT STRING
C
C OUTPUT  : (C**) STROUT  = OUTPUT STRING
C
C PROGRAM : (I*4) I       = INDEX COUNTER
C                 ICODE0  = BEGINNING OF CODE RANGE
C                 ICODE1  = FINISH OF CODE RANGE
c                 ICODE   = CHARACTER CODE
C                 IOFSET  = CODE OFF-SET
C
C METHOD  : A) EACH CHARACTER IN THE INPUT STRING IS CHECKED TO SEE IF
C              IT LIES IN THE LOWER CASE CODE RANGE
C           B) IS SO, THEN AN OFF-SET IS ADDED ONTO THE CODE OF THE
C              CHARACTER TO BRING IT INTO THE CORRESPONDING UPPER CODE
C              RANGE.
C
C ROUTINE : (I*4) IHAT    = 0 --- OUTSIDE INCLUSIVE RANGE
C                         = 1 --- INSIDE  INCLUSIVE RANGE
C
C
C AUTHOR  : J.SPENCE  (K1/0/80)  EXT. 4865
C           JET
C
C HISTORY : V1.R1.M0 --- 01/09/94 --- CREATION
C
C-- ....................................................................
C
C..INPUT
      CHARACTER STRIN*(*)
C
C..OUTPUT
      CHARACTER STROUT*(*)
C
C..PROGRAM
      INTEGER*4 I , ICODE0 , ICODE1 , ICODE , IOFSET , IHAT
C
C-----------------------------------------------------------------------
C                            INITIALISE
C-----------------------------------------------------------------------
C
C..START OF CODE RANGE
      ICODE0 = ICHAR('a')
C
C..END   OF CODE RANGE
      ICODE1 = ICHAR('z')
C
C..OFF-SET
      IOFSET = ICHAR('A') - ICHAR('a')
C
C-----------------------------------------------------------------------
C                    TEST EACH CHARACTER AND CONVERT
C-----------------------------------------------------------------------
C
      DO 100 I           = 1 , MIN0(LEN(STRIN),LEN(STROUT))
C
C..CODE OF I-TH CHARACTER
             ICODE       = ICHAR(STRIN(I:I))
C
C..CONVERT CODE
             ICODE       = ICODE + IOFSET*IHAT(ICODE,ICODE0,ICODE1)
C
C..COPY CONVERTED CHARACTER TO OUTPUT STRING
             STROUT(I:I) = CHAR(ICODE)
C
  100 CONTINUE
C
C-----------------------------------------------------------------------
C                       PAD STROUT WITH BLANKS
C-----------------------------------------------------------------------
C
      DO 200 I           = LEN(STRIN)+1 , LEN(STROUT)
         STROUT(I:I)     = ' '
  200 CONTINUE
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C
C=======================================================================
      SUBROUTINE CHRUTL( STRIN , STROUT )
      IMPLICIT NONE
C
C++ ....................................................................
C
C VERSION : V1.R1.M0
C
C ROUTINE : CHARACTER UPPER TO LOWER CONVERSION
C           -- -      -     -- -
C
C PURPOSE : TO CONVERT UPPER CASE CHARACTERS IN A STRING TO LOWER CASE.
C
C INPUT   : (C**) STRIN   = INPUT STRING
C
C OUTPUT  : (C**) STROUT  = OUTPUT STRING
C
C PROGRAM : (I*4) I       = INDEX COUNTER
C                 ICODE0  = BEGINNING OF CODE RANGE
C                 ICODE1  = FINISH OF CODE RANGE
c                 ICODE   = CHARACTER CODE
C                 IOFSET  = CODE OFF-SET
C
C METHOD  : A) EACH CHARACTER IN THE INPUT STRING IS CHECKED TO SEE IF
C              IT LIES IN THE UPPER CASE CODE RANGE
C           B) IS SO, THEN AN OFF-SET IS ADDED ONTO THE CODE OF THE
C              CHARACTER TO BRING IT INTO THE CORRESPONDING LOWER CODE
C              RANGE.
C
C ROUTINE : (I*4) IHAT    = 0 --- OUTSIDE INCLUSIVE RANGE
C                         = 1 --- INSIDE  INCLUSIVE RANGE
C
C AUTHOR  : J.SPENCE  (K1/0/80)  EXT. 4865
C           JET
C
C HISTORY : V1.R1.M0 --- 01/09/94 --- CREATION
C
C-- ....................................................................
C
C..INPUT
      CHARACTER STRIN*(*)
C
C..OUTPUT
      CHARACTER STROUT*(*)
C
C..PROGRAM
      INTEGER*4 I , ICODE0 , ICODE1 , ICODE , IOFSET , IHAT
C
C-----------------------------------------------------------------------
C                            INITIALISE
C-----------------------------------------------------------------------
C
C..START OF CODE RANGE
      ICODE0 = ICHAR('A')
C
C..END   OF CODE RANGE
      ICODE1 = ICHAR('Z')
C
C..OFF-SET
      IOFSET = ICHAR('a') - ICHAR('A')
C
C-----------------------------------------------------------------------
C                    TEST EACH CHARACTER AND CONVERT
C-----------------------------------------------------------------------
C
      DO 100 I           = 1 , MIN0(LEN(STRIN),LEN(STROUT))
C
C..CODE OF I-TH CHARACTER
             ICODE       = ICHAR(STRIN(I:I))
C
C..CONVERT CODE
             ICODE       = ICODE + IOFSET*IHAT(ICODE,ICODE0,ICODE1)
C
C..COPY CONVERTED CHARACTER TO OUTPUT STRING
             STROUT(I:I) = CHAR(ICODE)
C
  100 CONTINUE
C
C-----------------------------------------------------------------------
C                       PAD STROUT WITH BLANKS
C-----------------------------------------------------------------------
C
      DO 200 I           = LEN(STRIN)+1 , LEN(STROUT)
         STROUT(I:I)     = ' '
  200 CONTINUE
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C
C=======================================================================
      FUNCTION IHAT( ITEST , IMIN , IMAX )
      IMPLICIT NONE
      INTEGER*4 IHAT
C
C++ ....................................................................
C
C VERSION : V1.R1.M0
C
C ROUTINE : INTEGER "TOP HAT" FUNCTION
C           -            ---
C
C PURPOSE : TO RETURN 1 OR 0 DEPENDING ON I BEING IN THE RANGE IMIN TO
C           IMAX (INCLUSIVELY) OR OUTSIDE RESPECTIVELY.
C
C INPUT   : (I*4) ITEST   = INTEGER TO TEST
C           (I*4) IMIN    = MINIMUM VALUE OF RANGE
C           (I*4) IMAXT   = MAXIMUM VALUE OF RANGE
C
C OUTPUT  : (I*4) IHAT    = 0 --- IMIN >  ITEST >  IMAX
C                         = 1 --- IMIN <= ITEST <= IMAX
C
C METHOD  : USE ROUTINE 'ISIGN'.
C
C AUTHOR  : J.SPENCE  (K1/0/80)  EXT. 4865
C           JET
C
C HISTORY : V1.R1.M0 --- 01/09/94 --- CREATION
C
C-- ....................................................................
C
C..INPUT
      INTEGER*4 ITEST , IMIN , IMAX
C
C-----------------------------------------------------------------------
C                            CALCULATION
C             (IMAX+1 MAKES THE UPPER RANGE INCLUSIVE)
C-----------------------------------------------------------------------
C
      IHAT  = ( ISIGN( 1 , ITEST - MIN0(IMIN,IMAX)     )
     &        - ISIGN( 1 , ITEST - MAX0(IMIN,IMAX) - 1 ) ) / 2
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      RETURN
      END
c
c
c     Taken from EDGE2D - util.f
c
C
C=======================================================================
      SUBROUTINE ERRMSS(LUN, MODULE, IER, MSG1, MSG2, MSG3)
      implicit none
C
C***********************************************************************
C
C Prints out error message
C
C INPUT :   LUN    - Unit number for error message (0 --> default)
C           MODULE - Name of calling routine
C           IER    - Error level (0-warning, 1-error+exit, 2-error only)
C           MSGn   - Error message for line n (blank lines not printed)
C
C***********************************************************************
C
      CHARACTER*(*) MODULE,MSG1,MSG2,MSG3
      integer lun, ier
      integer ldef,lout,lenstr
      external lenstr
      DATA LDEF/6/
C
C
      IF (LUN.NE.0) THEN
          LOUT = LUN
      ELSE
          LOUT = LDEF
      END IF
C
      WRITE(LOUT,*)
      IF( IER.NE.0 )THEN
        WRITE(LOUT,*)
        WRITE(LOUT,'(4A)') '*** ERROR(',MODULE(1:LENSTR(MODULE)),') : '
     +                   , MSG1
      ELSE
        WRITE(LOUT,'(4A)')'*** WARNING(',MODULE(1:LENSTR(MODULE)),') : '
     +                   ,MSG1
      ENDIF
      IF (MSG2.NE.' ')
     +    WRITE(LOUT,'(2A)') '                    ',MSG2
      IF (MSG3.NE.' ')
     +    WRITE(LOUT,'(2A)') '                    ',MSG3
C
c      IF( IER.EQ.1 ) CALL EXITX(LOUT)
c
c      For now - if an error is encountered writing out the TRAN file then stop
c
      IF( IER.EQ.1 ) stop
C
C
      RETURN
      END


