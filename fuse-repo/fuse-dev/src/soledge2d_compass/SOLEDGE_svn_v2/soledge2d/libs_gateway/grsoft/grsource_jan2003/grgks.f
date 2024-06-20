C*
C* Copyright @ 1994, 1995   Josef Heinen
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
C*      GR/GR3 Software
C*
C* ABSTRACT:
C*
C*      This module contains the GLIGKS layer for the GR/GR3 software
C*      (Microsoft FORTRAN version).
C*
C* AUTHOR:
C*
C*      Josef Heinen
C*
C* VERSION:
C*
C*      V1.1-00
C*
C*

        subroutine grconf (devtype)

        integer devtype, closewk

        character env*100
        integer wstype

        wstype = -1
        if (devtype .ne. 35) then
            wstype = devtype
        else
            call getenv ('GRSOFT_DEVICE', env)
            if (env .ne. ' ') read (env, '(i6)') wstype
        end if

C*  set GLIGKS environment
        call putenv ('GLI_GKS=GKSGRAL')
        call putenv ('GLI_GKS_CMAP_EXTENT=')

C*  open Workstation Independent Segment Storage (WISS)
        call gopwk (1, 88, 5)
        call gacwk (1)
        call gcrsg (1)

        return

        entry grpan (closewk)
        closewk = closewk
C*  close segment and deactivate WISS
        call gclsg
        call gdawk (1)

C*  copy segment to IBM PC
        call gksgwk (0, 41)

C*  display output panel
        call gkpan (wstype)
        if (wstype .gt. 0) then
            open (unit=4, file='gr.out', status='unknown', recl=500)
            call gksgwk (4, wstype)
            close (unit=4)
        end if

C*  re-activate WISS
        call gacwk (1)
        call gdsg (1)
        call gcrsg (1)

        end

        subroutine gksgwk (conid, wstype)
C*  copy segment to workstation

        integer conid, wstype
        integer err, wkcat, stat, tnr
        real px, py

        call gkopwk (conid, wstype)

        call gcsgwk (2, 1)
        if (conid .eq. 0) then
            call gqwkca (wstype, err, wkcat)
            if (wkcat .eq. 2) call grqlc (2, 1, stat, tnr, px, py)
        end if

        call gkclwk

        end

        subroutine gkopwk (conid, wstype)
C*  open workstation

        integer conid, wstype

        integer errind, devunits, lx, ly
        real rx, ry, xcm, ycm, xratio, yratio
        integer i, coli
        real r(0:7), g(0:7), b(0:7)

        data r /1, 0, 1, 0, 0, 1, 1, 0/
        data g /1, 0, 0, 0, 1, 0, 1, 1/
        data b /1, 0, 0, 1, 0, 1, 0, 1/

        call gopwk (2, conid, wstype)
        call gacwk (2)
        call gqdsp (wstype, errind, devunits, rx, ry, lx, ly)

        call grqscl (xcm, ycm)
        if (rx .gt. xcm .and. ry .gt. ycm) then
            rx = xcm
            ry = ycm
        else
            xratio = xcm / rx
            yratio = ycm / ry
            if (yratio .gt. xratio) then
                rx = xcm / yratio
                ry = ycm / yratio
            else
                rx = xcm / xratio
                ry = ycm / xratio
            end if
        end if

        xmin = 0
        ymin = 0
        xmax = rx
        ymax = ry
C*  setup workstation viewport
        call gswkvp (2, xmin, xmax, ymin, ymax)

        if (xcm .gt. ycm) then
            xmax = 1
            ymax = ycm / xcm
        else
            xmax = xcm / ycm
            ymax = 1
        end if
C*  setup workstation window
        call gswkwn (2, xmin, xmax, ymin, ymax)

        do 1, i = 0, 7
            coli = i
C*  set color representation
            call gscr (2, coli, r(i), g(i), b(i))
    1   continue

        end

        subroutine gkclwk
C*  close workstation

        call gdawk (2)
        call gclwk (2)

        end

        subroutine gkpan (wstype)
C*  display workstation panel

        integer wstype

        write (*, *) '  0  no output'
        write (*, *) ' 53  HP-GL'
        write (*, *) ' 61  PostScript (b/w)'
        write (*, *) ' 62  Color PostScript'
        write (*, *) '104  PBM (Portable BitMap)'
        write (*, *) '?'
        read (*, *) wstype

        if (wstype .ne. 0 .and. wstype .ne. 53 .and. wstype .ne. 61
     * .and. wstype .ne. 62 .and. wstype .ne. 104) wstype = -1

        end
