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
C*      GKS Description Table
C               level of GKS

        INTEGER LEVEL
        PARAMETER (LEVEL=GL0B)
C               number of available workstation types
        INTEGER NWSTY
        PARAMETER (NWSTY=39)
C               max. number of simult. open workstations
        INTEGER MNOPWS
        PARAMETER (MNOPWS=16)
C               max. number of simult. active workstations
        INTEGER MNACWS
        PARAMETER (MNACWS=16)
C               max. number of workstations associated with segment
        INTEGER MNWSAS
        PARAMETER (MNWSAS=1)
C               maximum normalization transformation number
        INTEGER MXNTNR
        PARAMETER (MXNTNR=8)
