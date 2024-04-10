C*
C* Copyright @ 1984 - 1993   Josef Heinen
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

        SUBROUTINE GKCSG (WKID, SEGN, LIA, LR1, LR2, LC, CLEAR)

        INTEGER WKID, SEGN, LIA, LR1, LR2, LC, CLEAR

        INTEGER MAXBUF
        PARAMETER (MAXBUF = 2000)

        INTEGER IA(MAXBUF)
        REAL R1(MAXBUF), R2(MAXBUF)
        CHARACTER*132 CHARS(1)

        IF (LIA .LE. MAXBUF .AND. LR1 .LE. MAXBUF .AND. LR2 .LE. MAXBUF
     *      .AND. LC .LE. 1) THEN
            IF (CLEAR .NE. 0) CALL GCLRWK (WKID, 1)
	    CALL GKDCSG (WKID, SEGN, IA, R1, R2, CHARS)
        ELSE
            STOP 'GKS: not enough segment storage'
        END IF

        RETURN
        END
