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

        INTEGER FUNCTION GORD (SET,LENGTH,ELEMNT)
C               return order of an element

        INTEGER I, LENGTH
        INTEGER SET(LENGTH), ELEMNT

        DO 1 I = 1,LENGTH
          IF (SET(I) .EQ. ELEMNT) THEN
            GORD = I
            RETURN
          END IF
   1    CONTINUE

        GORD = 0

        RETURN
        END


        LOGICAL FUNCTION GANY (SET,LENGTH,ELEMNT)
C               test set for an element

        INTEGER I, LENGTH
        INTEGER SET(LENGTH), ELEMNT

        DO 1 I = 1,LENGTH
          IF (SET(I) .EQ. ELEMNT) THEN
            GANY = .TRUE.
            RETURN
          END IF
   1    CONTINUE

        GANY = .FALSE.

        RETURN
        END


        LOGICAL FUNCTION GALL (SET,LENGTH,ELEMNT)
C               compare all elements of a set

        INTEGER I, LENGTH
        INTEGER SET(LENGTH), ELEMNT

        DO 1 I = 1,LENGTH
          IF (SET(I) .NE. ELEMNT) THEN
            GALL = .FALSE.
            RETURN
          END IF
   1    CONTINUE

        GALL = .TRUE.

        RETURN
        END


        SUBROUTINE GADD (ELEMNT,SET,LENGTH)
C               add an element to a set
 
        INTEGER I, LENGTH
        INTEGER ELEMNT, SET(LENGTH)

        DO 1 I = 1,LENGTH
          IF (SET(I) .EQ. 0) THEN
            SET(I) = ELEMNT
            RETURN
          END IF
   1    CONTINUE

        RETURN
        END


        SUBROUTINE GDEL (ELEMNT,SET,LENGTH)
C               delete an element from a set

        INTEGER I, J, LENGTH
        INTEGER ELEMNT, SET(LENGTH)

        DO 1 I = 1,LENGTH
          IF (SET(I) .EQ. ELEMNT) THEN
            IF (I .LT. LENGTH) THEN
              DO 2 J = I,LENGTH-1
                SET(J) = SET(J+1)
   2          CONTINUE
            END IF
            SET(LENGTH) = 0
            RETURN
          END IF
   1    CONTINUE

        RETURN
        END


        SUBROUTINE GINI (SET,LENGTH)
C               initialize set

        INTEGER I, LENGTH
        INTEGER SET(LENGTH)

        DO 1 I = 1,LENGTH
C               do not remove the following IF-statement; prevents compiler
C               from optimizing DO-loop
          IF (SET(I) .NE. 0) SET(I) = 0
   1    CONTINUE

        RETURN
        END


        SUBROUTINE GNUM (SET,LENGTH,NUM)
C               number of elements

        INTEGER I, LENGTH, NUM
        INTEGER SET(LENGTH)

        NUM = 0
        DO 1 I = 1,LENGTH
          IF (SET(I) .NE. 0) NUM = NUM+1
   1    CONTINUE

        RETURN
        END


        SUBROUTINE GKQPA (PATTRN, PA)
C               inquire predefined pattern array

        INTEGER PATTRN, PA(0:32)

        INTEGER I, J, PATT(0:32,0:119)

        SAVE PATT

        DATA (PATT(0,  I), I =  0,  11) / 12 * 4/
        DATA (PATT(0,  I), I = 12, 119) /108 * 8/

        DATA (PATT(I,  0),I=1,4)/  0,   0,   0,   0/
        DATA (PATT(I,  1),I=1,4)/255, 255, 187, 255/
        DATA (PATT(I,  2),I=1,4)/238, 255, 187, 255/
        DATA (PATT(I,  3),I=1,4)/187, 238, 187, 255/
        DATA (PATT(I,  4),I=1,4)/119, 221, 119, 221/
        DATA (PATT(I,  5),I=1,4)/119, 170, 119, 221/
        DATA (PATT(I,  6),I=1,4)/ 85, 170,  85, 187/
        DATA (PATT(I,  7),I=1,4)/ 85, 170,  85, 170/
        DATA (PATT(I,  8),I=1,4)/136,  85, 136,  34/
        DATA (PATT(I,  9),I=1,4)/136,  34, 136,  34/
        DATA (PATT(I, 10),I=1,4)/  0,  68,   0,  17/
        DATA (PATT(I, 11),I=1,4)/  0,  34,   0,   0/
        DATA (PATT(I, 12),I=1,8)/  0,   0,   0,   0,   0,   0,   0,   0/
        DATA (PATT(I, 13),I=1,8)/128,   0,   8,   0, 128,   0,   8,   0/
        DATA (PATT(I, 14),I=1,8)/ 34, 136,  34, 136,  34, 136,  34, 136/
        DATA (PATT(I, 15),I=1,8)/ 85, 170,  85, 170,  85, 170,  85, 170/
        DATA (PATT(I, 16),I=1,8)/170,   0, 170,   0, 170,   0, 170,   0/
        DATA (PATT(I, 17),I=1,8)/ 85,  85,  85,  85,  85,  85,  85,  85/
        DATA (PATT(I, 18),I=1,8)/ 17,  34,  68, 136,  17,  34,  68, 136/
        DATA (PATT(I, 19),I=1,8)/119, 119, 119, 119, 119, 119, 119, 119/
        DATA (PATT(I, 20),I=1,8)/ 78, 207, 252, 228,  39,  63, 243, 114/
        DATA (PATT(I, 21),I=1,8)/127, 239, 253, 223, 254, 247, 191, 251/
        DATA (PATT(I, 22),I=1,8)/  0, 119, 119, 119,   0, 119, 119, 119/
        DATA (PATT(I, 23),I=1,8)/  0, 127, 127, 127,   0, 247, 247, 247/
        DATA (PATT(I, 24),I=1,8)/127, 255, 255, 255, 255, 255, 255, 255/
        DATA (PATT(I, 25),I=1,8)/127, 191, 223, 255, 253, 251, 247, 255/
        DATA (PATT(I, 26),I=1,8)/125, 187, 198, 187, 125, 254, 254, 254/
        DATA (PATT(I, 27),I=1,8)/  7, 139, 221, 184, 112, 232, 221, 142/
        DATA (PATT(I, 28),I=1,8)/170,  95, 191, 191, 170, 245, 251, 251/
        DATA (PATT(I, 29),I=1,8)/223, 175, 119, 119, 119, 119, 250, 253/
        DATA (PATT(I, 30),I=1,8)/ 64, 255,  64,  64,  79,  79,  79,  79/
        DATA (PATT(I, 31),I=1,8)/127, 255, 247, 255, 127, 255, 247, 255/
        DATA (PATT(I, 32),I=1,8)/119, 255, 221, 255, 119, 255, 221, 255/
        DATA (PATT(I, 33),I=1,8)/119, 221, 119, 221, 119, 221, 119, 221/
        DATA (PATT(I, 34),I=1,8)/ 85, 255,  85, 255,  85, 255,  85, 255/
        DATA (PATT(I, 35),I=1,8)/  0, 255,   0, 255,   0, 255,   0, 255/
        DATA (PATT(I, 36),I=1,8)/238, 221, 187, 119, 238, 221, 187, 119/
        DATA (PATT(I, 37),I=1,8)/  0, 255, 255, 255,   0, 255, 255, 255/
        DATA (PATT(I, 38),I=1,8)/254, 253, 251, 247, 239, 223, 191, 127/
        DATA (PATT(I, 39),I=1,8)/ 85, 255, 127, 255, 119, 255, 127, 255/
        DATA (PATT(I, 40),I=1,8)/  0, 127, 127, 127, 127, 127, 127, 127/
        DATA (PATT(I, 41),I=1,8)/247, 227, 221,  62, 127, 254, 253, 251/
        DATA (PATT(I, 42),I=1,8)/119, 235, 221, 190, 119, 255,  85, 255/
        DATA (PATT(I, 43),I=1,8)/191,  95, 255, 255, 251, 245, 255, 255/
        DATA (PATT(I, 44),I=1,8)/252, 123, 183, 207, 243, 253, 254, 254/
        DATA (PATT(I, 45),I=1,8)/127, 127, 190, 193, 247, 247, 235,  28/
        DATA (PATT(I, 46),I=1,8)/239, 223, 171,  85,   0, 253, 251, 247/
        DATA (PATT(I, 47),I=1,8)/136, 118, 112, 112, 136, 103,   7,   7/
        DATA (PATT(I, 48),I=1,8)/255, 247, 235, 213, 170, 213, 235, 247/
        DATA (PATT(I, 49),I=1,8)/255, 247, 235, 213, 170, 213, 235, 247/
        DATA (PATT(I, 50),I=1,8)/127, 255, 255, 255, 255, 255, 255, 255/
        DATA (PATT(I, 51),I=1,8)/127, 255, 255, 255, 247, 255, 255, 255/
        DATA (PATT(I, 52),I=1,8)/119, 255, 255, 255, 247, 255, 255, 255/
        DATA (PATT(I, 53),I=1,8)/119, 255, 255, 255, 119, 255, 255, 255/
        DATA (PATT(I, 54),I=1,8)/119, 255, 223, 255, 119, 255, 255, 255/
        DATA (PATT(I, 55),I=1,8)/119, 255, 223, 255, 119, 255, 253, 255/
        DATA (PATT(I, 56),I=1,8)/119, 255, 221, 255, 119, 255, 253, 255/
        DATA (PATT(I, 57),I=1,8)/119, 255, 221, 255, 119, 255, 221, 255/
        DATA (PATT(I, 58),I=1,8)/ 87, 255, 221, 255, 119, 255, 221, 255/
        DATA (PATT(I, 59),I=1,8)/ 87, 255, 221, 255, 117, 255, 221, 255/
        DATA (PATT(I, 60),I=1,8)/ 85, 255, 221, 255, 117, 255, 221, 255/
        DATA (PATT(I, 61),I=1,8)/ 85, 255, 221, 255,  85, 255, 221, 255/
        DATA (PATT(I, 62),I=1,8)/ 85, 255,  93, 255,  85, 255, 221, 255/
        DATA (PATT(I, 63),I=1,8)/ 85, 255,  93, 255,  85, 255, 213, 255/
        DATA (PATT(I, 64),I=1,8)/ 85, 255,  85, 255,  85, 255, 213, 255/
        DATA (PATT(I, 65),I=1,8)/ 85, 255,  85, 255,  85, 255,  85, 255/
        DATA (PATT(I, 66),I=1,8)/ 85, 191,  85, 255,  85, 255,  85, 255/
        DATA (PATT(I, 67),I=1,8)/ 85, 191,  85, 255,  85, 251,  85, 255/
        DATA (PATT(I, 68),I=1,8)/ 85, 187,  85, 255,  85, 251,  85, 255/
        DATA (PATT(I, 69),I=1,8)/ 85, 187,  85, 255,  85, 187,  85, 255/
        DATA (PATT(I, 70),I=1,8)/ 85, 187,  85, 239,  85, 187,  85, 255/
        DATA (PATT(I, 71),I=1,8)/ 85, 187,  85, 239,  85, 187,  85, 254/
        DATA (PATT(I, 72),I=1,8)/ 85, 187,  85, 238,  85, 187,  85, 254/
        DATA (PATT(I, 73),I=1,8)/ 85, 187,  85, 238,  85, 187,  85, 238/
        DATA (PATT(I, 74),I=1,8)/ 85, 171,  85, 238,  85, 187,  85, 238/
        DATA (PATT(I, 75),I=1,8)/ 85, 171,  85, 238,  85, 186,  85, 238/
        DATA (PATT(I, 76),I=1,8)/ 85, 170,  85, 238,  85, 186,  85, 238/
        DATA (PATT(I, 77),I=1,8)/ 85, 170,  85, 238,  85, 170,  85, 238/
        DATA (PATT(I, 78),I=1,8)/ 85, 170,  85, 174,  85, 170,  85, 238/
        DATA (PATT(I, 79),I=1,8)/ 85, 170,  85, 174,  85, 170,  85, 234/
        DATA (PATT(I, 80),I=1,8)/ 85, 170,  85, 170,  85, 170,  85, 234/
        DATA (PATT(I, 81),I=1,8)/ 85, 170,  85, 170,  85, 170,  85, 170/
        DATA (PATT(I, 82),I=1,8)/ 21, 170,  85, 170,  85, 170,  85, 170/
        DATA (PATT(I, 83),I=1,8)/ 21, 170,  85, 170,  81, 170,  85, 170/
        DATA (PATT(I, 84),I=1,8)/ 17, 170,  85, 170,  81, 170,  85, 170/
        DATA (PATT(I, 85),I=1,8)/ 17, 170,  85, 170,  17, 170,  85, 170/
        DATA (PATT(I, 86),I=1,8)/ 17, 170,  69, 170,  17, 170,  85, 170/
        DATA (PATT(I, 87),I=1,8)/ 17, 170,  69, 170,  17, 170,  84, 170/
        DATA (PATT(I, 88),I=1,8)/ 17, 170,  68, 170,  17, 170,  84, 170/
        DATA (PATT(I, 89),I=1,8)/ 17, 170,  68, 170,  17, 170,  68, 170/
        DATA (PATT(I, 90),I=1,8)/  1, 170,  68, 170,  17, 170,  68, 170/
        DATA (PATT(I, 91),I=1,8)/  1, 170,  68, 170,  16, 170,  68, 170/
        DATA (PATT(I, 92),I=1,8)/  0, 170,  68, 170,  16, 170,  68, 170/
        DATA (PATT(I, 93),I=1,8)/  0, 170,  68, 170,   0, 170,  68, 170/
        DATA (PATT(I, 94),I=1,8)/  0, 170,   4, 170,   0, 170,  68, 170/
        DATA (PATT(I, 95),I=1,8)/  0, 170,   4, 170,   0, 170,  64, 170/
        DATA (PATT(I, 96),I=1,8)/  0, 170,   0, 170,   0, 170,  64, 170/
        DATA (PATT(I, 97),I=1,8)/  0, 170,   0, 170,   0, 170,   0, 170/
        DATA (PATT(I, 98),I=1,8)/  0,  42,   0, 170,   0, 170,   0, 170/
        DATA (PATT(I, 99),I=1,8)/  0,  42,   0, 170,   0, 162,   0, 170/
        DATA (PATT(I,100),I=1,8)/  0,  34,   0, 170,   0, 162,   0, 170/
        DATA (PATT(I,101),I=1,8)/  0,  34,   0, 170,   0,  34,   0, 170/
        DATA (PATT(I,102),I=1,8)/  0,  34,   0, 138,   0,  34,   0, 170/
        DATA (PATT(I,103),I=1,8)/  0,  34,   0, 138,   0,  34,   0, 168/
        DATA (PATT(I,104),I=1,8)/  0,  34,   0, 136,   0,  34,   0, 168/
        DATA (PATT(I,105),I=1,8)/  0,  34,   0, 136,   0,  34,   0, 136/
        DATA (PATT(I,106),I=1,8)/  0,   2,   0, 136,   0,  34,   0, 136/
        DATA (PATT(I,107),I=1,8)/  0,   2,   0, 136,   0,  32,   0, 136/
        DATA (PATT(I,108),I=1,8)/  0,   0,   0, 136,   0,  32,   0, 136/
        DATA (PATT(I,109),I=1,8)/119, 119, 119, 119, 119, 119, 119, 119/
        DATA (PATT(I,110),I=1,8)/  0, 255, 255, 255,   0, 255, 255, 255/
        DATA (PATT(I,111),I=1,8)/119, 187, 221, 238, 119, 187, 221, 238/
        DATA (PATT(I,112),I=1,8)/238, 221, 187, 119, 238, 221, 187, 119/
        DATA (PATT(I,113),I=1,8)/  0, 119, 119, 119,   0, 119, 119, 119/
        DATA (PATT(I,114),I=1,8)/126, 189, 219, 231, 231, 219, 189, 126/
        DATA (PATT(I,115),I=1,8)/247, 247, 247, 247, 247, 247, 247, 247/
        DATA (PATT(I,116),I=1,8)/255, 255, 255, 255,   0, 255, 255, 255/
        DATA (PATT(I,117),I=1,8)/127, 191, 223, 239, 247, 251, 253, 254/
        DATA (PATT(I,118),I=1,8)/254, 253, 251, 247, 239, 223, 191, 127/
        DATA (PATT(I,119),I=1,8)/247, 247, 247, 247,   0, 247, 247, 247/

        J = PATTRN
        IF (J .GT. 119) THEN
          J = 119
        ELSE IF (J .LT. 0) THEN
          J = 0
        END IF

        K = PATT(0,J)
        DO 1, I = 0, K
          PA(I) = PATT(I,J)
    1   CONTINUE

        RETURN


        ENTRY GKSPA (PATTRN, PA)
C               set pattern array

        IF (PATTRN .GE. 0 .AND. PATTRN .LE. 119) THEN
          K = PA(0)
          IF (K .EQ. 4 .OR. K .EQ. 8 .OR. K .EQ. 32) THEN
            DO 2, I = 0, K
              PATT(I,PATTRN) = PA(I)
    2       CONTINUE
          END IF
        END IF

        RETURN
        END     
