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

      SUBROUTINE GKDTK1 (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS)

      INTEGER LC, LR1, LR2
      INTEGER FCTID, DX, DY, DIMX, IA(3)
      REAL R1(3), R2(3)
      CHARACTER*(*) CHARS(1)
      INTEGER IPTR

      INTEGER L

      ENTRY GKDTK2 (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS)

      ENTRY GKDVT3 (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS)

      WRITE (*, *) FCTID, DX, DY, DIMX, (IA(L), L = 1, DX * DY),
     * LR1, (R1(L), L = 1, LR1), LR2, (R2(L), L=1, LR2),
     * LC, (CHARS(L), L = 1, LC)

      STOP 'GKS: logical device driver not supported on this system'

      RETURN

      ENTRY GKDCGM (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS, IPTR)

      ENTRY GKDMFO (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS, IPTR)

      ENTRY GKDXW (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS, IPTR)

      ENTRY GKDPDF (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS, IPTR)

      ENTRY GKDMFI (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     * LC, CHARS, IPTR)

      WRITE (*, *) FCTID, DX, DY, DIMX, (IA(L), L = 1, DX * DY),
     * LR1, (R1(L), L = 1, LR1), LR2, (R2(L), L=1, LR2),
     * LC, (CHARS(L), L = 1, LC), IPTR

      STOP 'GKS: logical device driver not supported on this system'

      RETURN
      END
