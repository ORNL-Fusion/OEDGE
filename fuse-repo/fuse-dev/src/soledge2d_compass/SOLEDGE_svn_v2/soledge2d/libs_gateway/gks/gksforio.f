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
C*

        SUBROUTINE FORTWR (CONID, NCHARS, CHARS, STATUS)

C*  Output a variable length record

        INTEGER CONID, NCHARS, STATUS
        CHARACTER CHARS*(*)

        INTEGER RECL, I, N
        CHARACTER BUF*500

        INQUIRE (UNIT=CONID, RECL=RECL)

        IF (RECL .LT. 500) THEN
            N = 0
            DO 100, I = 1, NCHARS
                IF (CHARS(I:I) .NE. CHAR(10)) THEN
                    N = N + 1
                    BUF(N:N) = CHARS(I:I)
                ELSE
                    IF (N .GT. 0) THEN
                        WRITE(CONID, '(A)', IOSTAT=STATUS) BUF(1:N)
                        N = 0
                    ELSE
                        WRITE(CONID, '(A)', IOSTAT=STATUS)
                    END IF
                END IF
 100        CONTINUE
            IF (N .GT. 0) WRITE(CONID, '(A)', IOSTAT=STATUS) BUF(1:N)
        ELSE
            WRITE (CONID, '(A)', IOSTAT=STATUS) CHARS (1:NCHARS)
        END IF

        IF (STATUS .GT. 0) THEN
            STATUS = -STATUS
        ELSE IF (STATUS .EQ. 0) THEN
            STATUS = NCHARS
        END IF

        RETURN
        END
