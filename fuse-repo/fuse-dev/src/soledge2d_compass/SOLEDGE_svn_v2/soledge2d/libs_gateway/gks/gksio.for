C*
C* Copyright @ 1984 - 1994   Josef Heinen
C*
C* Permission to use, copy, and distribute this software and its
C* documentation for any purpose with or without fee is hereby granted,
C* provided that the above copyright notice appear in all copies and
C* that both that copyright notice and this permission notice appear
C* in supporting documentation.
C*
C* Permission to modify the software is granted, but not the right to
C* distribute the modified code.  Modifications are to be distributed
C* as patches to released version.
C*
C* This software is provided "as is" without express or implied warranty.
C*
C* Send your comments or suggestions to
C*  J.Heinen@kfa-juelich.de.
C*
C*
C* FACILITY:
C*
C*      Graphical Kernel System (GKS)
C*
C* ABSTRACT:
C*
C*      This module contains some I/O routines for GLI GKS.
C*      (Microsoft FORTRAN version)
C*
C* AUTHOR(S):
C*
C*      J. Heinen
C*
C* VERSION:
C*
C*      V1.0
C*
C*

	SUBROUTINE OPFONT
C               open file with character fonts

	INTEGER FONT, CH, BUFFER(256)

	INTEGER UNIT
	PARAMETER (UNIT = 80)

	INTEGER HASHT(0:95)
	INTEGER*1 BUF(256), BUFC(256,0:95)

	INTEGER STATUS
	INTEGER MAP(24), SMAP(24), GRMAP(24), GKMAP(24), GERMAN(11)
	INTEGER*1 ANSI(11), GREEK(14), GMAP(14)

	LOGICAL UMLAUT, SHARPS
	INTEGER CHR, OFFSET, FN

	CHARACTER ENV*80, DRIVE(5)
	INTEGER L, GKS

	SAVE

	DATA DRIVE  / 'C', 'D', 'E', 'A', 'B' /

	DATA MAP    / 1, 18,  1,  6, 12,  3,  8, 11,  4,  7,
     *               10,  2, 13, 14,  5,  9, 15, 16, 17, 20,
     *               21, 19, 22, 23 /
	DATA GRMAP  / 1, 12,  6,  9,  8, 11,  5, 13, 18, 17,
     *               19,  1, 10,  7, 24,  1,  1,  1,  1,  1,
     *                1,  1, 23, 24 /
	DATA GKMAP  / 1,  2,  3,  4,  5,  6,  7, 13, 11, 14,
     *               12, 15, 16, 17, 18, 19, 20, 21,  1,  1,
     *                1,  1, 23, 24 /
	DATA SMAP   / 4,  4,  4,  4,  4,  7,  7,  7, 10, 10,
     *               10,  7,  7,  7,  4,  4,  7,  7,  7,  4,
     *                4,  4,  4,  4 /

	DATA GERMAN / 196, 214, 220, 228, 246, 252,
     *                223, 171, 187, 183, 169 /
	DATA ANSI   / 'A', 'O', 'U', 'a', 'o', 'u',
     *                'b', '<', '>', '.', '@' /
	DATA GREEK  / 'j', 'o', 'q', 'u', 'v', 'w', 'y',
     *                'J', 'O', 'Q', 'U', 'V', 'W', 'Y' /
	DATA GMAP   / ' ', 'w', ' ', 'o', 'y', 'v', 'q',
     *                ' ', 'W', ' ', 'O', 'Y', 'V', 'Q' /

	CALL GETENV('GLI_HOME', ENV)
	IF (ENV(1:1) .EQ. ' ') ENV = 'C:\GKS\'
	L = INDEX(ENV, ' ') - 1

$IF DEFINED(F77L3)
	OPEN (UNIT=UNIT, FILE=ENV(1:L)//'GKSFONT.DAT', 
     *      STATUS='OLD', FORM='UNFORMATTED', ACCESS='TRANSPARENT', 
     *      IOSTAT=STATUS)
	DO 15, I = 1, 5
	    IF (STATUS .NE. 0) OPEN (UNIT=UNIT,
     *		FILE=DRIVE(I)//':\GKS\GKSFONT.DAT', STATUS='OLD',
     *		FORM='UNFORMATTED', ACCESS='TRANSPARENT',
     *		IOSTAT=STATUS)
  15	CONTINUE
$ELSE
	OPEN (UNIT=UNIT, FILE=ENV(1:L)//'GKSFONT.DAT',
     *      STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=256, 
     *      IOSTAT=STATUS)
	DO 15, I = 1, 5
	    IF (STATUS .NE. 0) OPEN (UNIT=UNIT,
     *		FILE=DRIVE(I)//':\GKS\GKSFONT.DAT', STATUS='OLD',
     *		FORM='UNFORMATTED', ACCESS='DIRECT', RECL=256,
     *		IOSTAT=STATUS)
  15	CONTINUE
$ENDIF

	CALL GETENV ('GLIGKS', ENV)
	IF (ENV .EQ. 'GKGKS') THEN
	    GKS = 2
	ELSE IF (ENV .EQ. 'GKSGRAL') THEN
	    GKS = 1
	ELSE
	    GKS = 0
	END IF

        DO 1, I = 1, 256
            HASHT(I) = -1
   1    CONTINUE

	RETURN


	ENTRY LOOKUP (FONT, CH, BUFFER)
C               inquire character stroke

	CHR = CH
	UMLAUT = .FALSE.
	SHARPS = .FALSE.

	IF (CHR .GE. 127) THEN
	    DO 2, I = 1, 11
		IF (CHR. EQ. GERMAN(I)) THEN
		    CHR = ANSI(I)
		    IF (I .LE. 6) THEN
			UMLAUT = .TRUE.
		    ELSE IF (I .EQ. 7) THEN
			SHARPS = .TRUE.
		    END IF
		END IF
   2        CONTINUE
	END IF

	IF (CHR .LT. 32 .OR. CHR .GE. 127) CHR = 32

	FN = MOD(ABS(FONT), 100)
	IF (FN .EQ. 51) THEN
	    FN = 23
	ELSE IF (FN .GT. 23) THEN
	    FN = 1
	END IF

	IF (SHARPS) THEN
	    IF (FN .NE. 23) THEN
		FN = SMAP(FN)
	    ELSE
		CHR = 126
	    END IF
	ELSE IF (GKS .EQ. 1) THEN
            IF (FN .EQ. 13 .OR. FN .EQ. 14) THEN
                DO 3, I = 1, 14
                    IF (CHR .EQ. GREEK(I)) THEN
                        CHR = GMAP(I)
                        GO TO 4
                    END IF
   3            CONTINUE
   4        END IF
	    FN = GRMAP(FN)
	ELSE IF (GKS .EQ. 2) THEN
	    FN = GKMAP(FN)
	END IF

        CHR = CHR - ICHAR(' ')
$IF DEFINED(F77L3)
	OFFSET = 256 * ((MAP(FN) - 1) * 95 + CHR) + 1
$ELSE
	OFFSET = ((MAP(FN) - 1) * 95 + CHR) + 1
$ENDIF
        IF (HASHT(CHR) .NE. OFFSET) THEN
	    READ (UNIT=UNIT, REC=OFFSET) BUF
            DO 5, I = 1, 256
                BUFC(I, CHR) = BUF(I)
   5        CONTINUE
            HASHT(CHR) = OFFSET
        END IF

	DO 6, I = 1, 256
	    BUFFER(I) = BUFC(I, CHR)
   6    CONTINUE

C               append umaut
	IF (UMLAUT .AND. BUFFER(8) .LT. 120-20)
     *      BUFFER(8) = BUFFER(8) + 10

	RETURN


	ENTRY CLFONT
C               close file with character fonts

	CLOSE (UNIT=UNIT)

	RETURN
	END


	SUBROUTINE BUFIN (CONID, STATUS, TERM, NCHARS, CHARS)

	INTEGER CONID, STATUS, TERM, NCHARS
	CHARACTER*(*) CHARS

	IF (TERM .EQ. 0) THEN
	    READ (CONID, '(A)') CHARS(1:NCHARS)
	    STATUS = 1
	ELSE
	    STATUS = 0
	END IF

	RETURN
	END


	SUBROUTINE BUFOUT (CONID, NCHARS, CHARS)

	INTEGER CONID, NCHARS
	CHARACTER*(*) CHARS

	INTEGER N
	CHARACTER BUF*500, NAME*100, FMT*7

	INQUIRE (UNIT=CONID, NAME=NAME)
	IF (NAME .EQ. 'CON') THEN
	    FMT = '(1X, A)'
	ELSE
	    FMT = '(A)'
	END IF

	N = 0
	DO 1, I = 1, NCHARS
	    IF (CHARS(I:I) .EQ. CHAR(10) .OR. N .EQ. 500) THEN
		IF (N .GT. 0) THEN
		    WRITE(CONID, FMT) BUF(1:N)
		    N = 0
		ELSE
		    WRITE(CONID, FMT)
		END IF
	    ELSE
		N = N + 1
		BUF(N:N) = CHARS(I:I)
	    END IF
   1    CONTINUE

	IF (N .GT. 0) WRITE(CONID, FMT) BUF(1:N)
 
	RETURN
	END


	SUBROUTINE BINOUT (CONID, NCHARS, CHARS)

	INTEGER CONID, NCHARS
	CHARACTER*(*) CHARS

	CHARACTER NAME*100, FORM*20
	INTEGER I, J, K, L

	INQUIRE (UNIT=CONID, NAME=NAME, FORM=FORM)
	CLOSE (UNIT=CONID)

$IF DEFINED(F77L3)
	OPEN (UNIT=CONID, FILE=NAME, STATUS='UNKNOWN',
     *	    ACCESS='TRANSPARENT')
$ELSE
	OPEN (UNIT=CONID, FILE=NAME, STATUS='UNKNOWN', FORM='BINARY')
$ENDIF
	J = 1
	K = 0
	L = NCHARS / 8192
	DO 1, I = 1, L
	    K = K + 8192
	    WRITE (CONID) CHARS(J:K)
	    J = K + 1
   1	CONTINUE
	IF (K .LT. NCHARS) WRITE (CONID) CHARS(J:NCHARS)

	RETURN
	END


	SUBROUTINE GKINFO (NCHARS, CHARS)

	INTEGER NCHARS
	CHARACTER*(*) CHARS

	NCHARS = 0
	CHARS = ' '

	RETURN
	END


	SUBROUTINE GKMAGS (MAGS, DPI)

	REAL MAGS
	INTEGER DPI

	MAGS = 0.0
	DPI = 75

	RETURN
	END


	SUBROUTINE LIBSIG (STATUS, ARG)

	INTEGER STATUS
	CHARACTER*(*) ARG
	INTEGER ERRNUM, ERRFIL

	INTEGER ERRLUN, ERRNO
	CHARACTER*6 FCTID

	ERRNO = MOD(STATUS/8, 1024)
	FCTID = ARG

	GOTO (  1,  2,  3,  4,  5,  6,  7,  8,999,999,
     *        999,999,999,999,999,999,999,999,999, 20,
     *         21, 22,999, 24, 25, 26, 27, 28, 29, 30,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999, 50,
     *         51, 52, 53,999,999,999,999,999,999, 60,
     *        999, 62,999, 64, 65, 66,999, 68,999, 70,
     *        999, 72, 73, 74, 75,999,999, 78,999,999,
     *         81,999,999, 84, 85,999,999, 88,999,999,
     *        999,999,999,999,999,999,999,999,999,100), ERRNO
	GOTO 999

   1    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * in the state GKCL', FCTID
	GOTO 999

   2    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * in the state GKOP', FCTID
	GOTO 999

   3    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * in the state WSAC', FCTID
	GOTO 999

   4    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * in the state SGOP', FCTID
	GOTO 999

   5    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * either in the state WSAC or SGOP', FCTID
	GOTO 999

   6    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * either in the state WSOP or WSAC', FCTID
	GOTO 999
      
   7    WRITE (ERRLUN, *) ERRNO, 'GKS not in proper state. GKS must be
     * in one of the states GKOP,WSOP,WSAC,SGOP', FCTID
	GOTO 999

   8    WRITE (ERRLUN, *) ERRNO, 'Specified workstation identifier is
     * invalid', FCTID
	GOTO 999

  20    WRITE (ERRLUN, *) ERRNO, 'Specified workstation identifier is
     * invalid', FCTID
	GOTO 999
 
  21    WRITE (ERRLUN, *) ERRNO, 'Specified connection identifier is
     * invalid', FCTID
	GOTO 999

  22    WRITE (ERRLUN, *) ERRNO, 'Specified workstation type is invalid'
     * , FCTID
	GOTO 999

  24    WRITE (ERRLUN, *) ERRNO, 'Specified workstation is open', FCTID
	GOTO 999

  25    WRITE (ERRLUN, *) ERRNO, 'Specified workstation is not open'
     * , FCTID
	GOTO 999

  26    WRITE (ERRLUN, *) ERRNO, 'Specified workstation cannot be
     * opened', FCTID
	GOTO 999
 
  27    WRITE (ERRLUN, *) ERRNO, 'Workstation Independent Segment
     * Storage is not open', FCTID
	GOTO 999

  28    WRITE (ERRLUN, *) ERRNO, 'Workstation Independent Segment
     * Storage is already open', FCTID
	GOTO 999

  29    WRITE (ERRLUN, *) ERRNO, 'Specified workstation is active'
     * , FCTID
	GOTO 999

  30    WRITE (ERRLUN, *) ERRNO, 'Specified workstation is not active'
     * , FCTID
	GOTO 999

  50    WRITE (ERRLUN, *) ERRNO, 'Transformation number is invalid'
     * , FCTID
	GOTO 999

  51    WRITE (ERRLUN, *) ERRNO, 'Rectangle definition is invalid'
     * , FCTID
	GOTO 999

  52    WRITE (ERRLUN, *) ERRNO, 'Viewport is not within the NDC
     * unit square', FCTID
	GOTO 999

  53    WRITE (ERRLUN, *) ERRNO, 'Workstation window is not within
     * the NDC unit square', FCTID
	GOTO 999

  60    WRITE (ERRLUN, *) ERRNO, 'Polyline index is invalid', FCTID
	GOTO 999

  62    WRITE (ERRLUN, *) ERRNO, 'Linetype is invalid', FCTID
	GOTO 999

  64    WRITE (ERRLUN, *) ERRNO, 'Polymarker index is invalid', FCTID
	GOTO 999

  65    WRITE (ERRLUN, *) ERRNO, 'Colour index is invalid', FCTID
	GOTO 999

  66    WRITE (ERRLUN, *) ERRNO, 'Marker type is invalid', FCTID
	GOTO 999
  68    WRITE (ERRLUN, *) ERRNO, 'Text index is invalid', FCTID
	GOTO 999

  70    WRITE (ERRLUN, *) ERRNO, 'Text font is invalid', FCTID
	GOTO 999

  72    WRITE (ERRLUN, *) ERRNO, 'Character expansion factor is invalid'
     * , FCTID
	GOTO 999

  73    WRITE (ERRLUN, *) ERRNO, 'Character height is invalid', FCTID
	GOTO 999

  74    WRITE (ERRLUN, *) ERRNO, 'Character up vector is invalid'
     * , FCTID
	GOTO 999
 
  75    WRITE (ERRLUN, *) ERRNO, 'Fill area index is invalid', FCTID
	GOTO 999

  78    WRITE (ERRLUN, *) ERRNO, 'Style index is invalid', FCTID
	GOTO 999

  81    WRITE (ERRLUN, *) ERRNO, 'Pattern size value is invalid', FCTID
	GOTO 999

  84    WRITE (ERRLUN, *) ERRNO, 'Dimensions of colour index array are
     * invalid', FCTID
	GOTO 999

  85    WRITE (ERRLUN, *) ERRNO, 'Colour index is invalid', FCTID
	GOTO 999

  88    WRITE (ERRLUN, *) ERRNO, 'Colour is invalid', FCTID
	GOTO 999

 100    WRITE (ERRLUN, *) ERRNO, 'Number of points is invalid', FCTID
	GOTO 999

 200    FORMAT (1X, 'GKS error ', I4, '. ', A, ' in routine ', A)

 999    RETURN


	ENTRY GERSET (ERRNUM, ERRFIL)

	ERRNO = ERRNUM
	ERRLUN = ERRFIL

	RETURN
	END


	SUBROUTINE GKFD (CONID)

	INTEGER CONID

	IF (CONID .LT. 0) STOP 'GKS: invalid connection identifier'

	RETURN
	END


	SUBROUTINE GTWSTY (WSTYPE)

	INTEGER WSTYPE

	WSTYPE = 41

	RETURN
	END


	SUBROUTINE GKLTOI (L, I)

	INTEGER L(2), I(2)

	I(1) = L(1)
	I(2) = L(2)

	RETURN
	END


	SUBROUTINE PUTENV (STRING)
C               change or add value to environment

	CHARACTER*(*) STRING, NAME, VALUE

	INTEGER MAXENV, MAXSTR
	PARAMETER (MAXENV = 20, MAXSTR = 80)

	CHARACTER*80 ENV(MAXENV)
	INTEGER I, J, K, N

	SAVE

	DATA N / 0 /

	I = INDEX(STRING, '=') - 1
	IF (I .GT. 0 .AND. I .LT. MAXSTR) THEN
	    DO 1, J = 1,N
		IF (STRING(1:I) .EQ. ENV(N)(1:I)) THEN
		    K = J
		    GO TO 2
		END IF
   1        CONTINUE
	    IF (N .EQ. MAXENV) STOP 'GKS: environment full'
	    N = N + 1
	    K = N
   2        ENV(K) = STRING
	END IF

	RETURN


	ENTRY GETENV (NAME, VALUE)
C               get environment variable

	VALUE = ' '
	I = LEN(NAME)
	IF (I .GT. 0 .AND. I .LT. MAXSTR) THEN
	    DO 3, J = 1,N
		IF (NAME .EQ. ENV(J)(1:I)) THEN
		    VALUE = ENV(J)(I+2:)
		    GO TO 4
		END IF
   3        CONTINUE
	END IF

   4    RETURN
	END


	SUBROUTINE DPSOP (SIZEX, SIZEY, FMT)
C		initiate Display PostScript environment

	REAL SIZEX, SIZEY
	INTEGER FMT, NCHARS
	CHARACTER*(*) CHARS

	WRITE(*, *) SIZEX, SIZEY, FMT

	STOP 'GKS: Display PostScript not supported on this system'

	RETURN

	ENTRY DPSWR (NCHARS, CHARS)
	
	WRITE(*, *) NCHARS, CHARS

	STOP 'GKS: Display PostScript not supported on this system'

	RETURN

	ENTRY DPSFL

	ENTRY DPSCL

	END


	SUBROUTINE GKTMP (NCHARS, CHARS)

	INTEGER NCHARS
	CHARACTER CHARS*(*)

	NCHARS = 10
	CHARS = 'GLIGKS.BIN'

	RETURN
	END
