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
C               operating state value
        INTEGER STATE

C *     GKS State List
C               set of open workstations
        INTEGER SOPWS(MNOPWS,6)
C               number of open workstations
        INTEGER NOPWS
C               set of active workstations
        INTEGER SACWS(MNACWS)
C               number of active workstations
        INTEGER NACWS
C               aspect source flags
        INTEGER ASF(13)
C               character height
        REAL CHH
C               character up vector
        REAL CHUP(2)
C               text path
        INTEGER TXP
C               text alignment
        INTEGER TXAL(2)
C               polyline index
        INTEGER LINDEX
C               linetype
        INTEGER LTYPE
C               linewidth scale factor
        REAL LWIDTH
C               polyline colour index
        INTEGER PLCOLI
C               polymarker index
        INTEGER MINDEX
C               marker type
        INTEGER MTYPE
C               marker size and scale factor
        REAL MSZSC
C               polymarker colour index
        INTEGER PMCOLI
C               text index
        INTEGER TINDEX
C               text font and precision
        INTEGER TXFONT,TXPREC
C               character expansion factor
        REAL CHXP
C               character spacing
        REAL CHSP
C               text colour index
        INTEGER TXCOLI
C               fill area index
        INTEGER FINDEX
C               fill area interior style
        INTEGER INTS
C               fill area style index
        INTEGER STYLI
C               fill area colour index
        INTEGER FACOLI
C               current normalization transformation number
        INTEGER CNTNR
C               list of normalization transformation numbers
        INTEGER LNTNR(0:MXNTNR)
C               normalization transformation
        REAL WINDOW(4,0:MXNTNR),VIEWPT(4,0:MXNTNR)
C               clipping indicator
        INTEGER CLIP
C               transformation matrix
        REAL MAT(2,3)
C               name of open segment
        INTEGER OPSG

C *     Miscellaneous stuff
C               list of available workstation types
        INTEGER LWSTY(NWSTY)
C               workstation categories
        INTEGER WSCAT(NWSTY)
C               maximum display surface size (DC)
        REAL GDC(2,NWSTY)
C               maximum display surface size (raster units)
        INTEGER GRU(2,NWSTY)
C               X11 flag
        LOGICAL X11
C               kernel flag
        LOGICAL KERNEL

C               GKS state list structure
        INTEGER GKSL(-3:118), IIA(4)

        EQUIVALENCE (IIA,GKSL(-3))
        EQUIVALENCE (LINDEX,GKSL(1)), (LTYPE,GKSL(2)), (LWIDTH,GKSL(3))
        EQUIVALENCE (PLCOLI,GKSL(4))
        EQUIVALENCE (MINDEX,GKSL(5)), (MTYPE,GKSL(6)), (MSZSC,GKSL(7))
        EQUIVALENCE (PMCOLI,GKSL(8))
        EQUIVALENCE (TINDEX,GKSL(9)), (TXFONT,GKSL(10))
        EQUIVALENCE (TXPREC,GKSL(11)), (CHXP,GKSL(12)), (CHSP,GKSL(13))
        EQUIVALENCE (TXCOLI,GKSL(14)), (CHH,GKSL(15)), (CHUP,GKSL(16))
        EQUIVALENCE (TXP,GKSL(18)), (TXAL,GKSL(19))
        EQUIVALENCE (FINDEX,GKSL(21)), (INTS,GKSL(22)), (STYLI,GKSL(23))
        EQUIVALENCE (FACOLI,GKSL(24))
        EQUIVALENCE (WINDOW,GKSL(25)), (VIEWPT,GKSL(61))
        EQUIVALENCE (CNTNR,GKSL(97)), (CLIP,GKSL(98)), (OPSG,GKSL(99))
        EQUIVALENCE (MAT,GKSL(100)), (ASF,GKSL(106))

        COMMON /GKS/ GKSL, STATE, SOPWS, NOPWS, SACWS, NACWS, LNTNR,
     *               WSCAT, LWSTY, GDC, GRU, X11, KERNEL
CDEC$ PSECT /GKS/ NOSHR

