C*
C* Copyright @ 1984 - 1995   Josef Heinen
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

        SUBROUTINE GKDPS (FCTID,DX,DY,DIMX,IA,LR1,R1,LR2,R2,LC,CHARS)
C*  GKS logical device driver for PostScript printers

        INTEGER LC, LR1, LR2
        INTEGER FCTID,DX,DY,DIMX,IA(3)
        REAL R1(3),R2(3)
        CHARACTER*(*) CHARS

        REAL PI
        PARAMETER (PI = 3.141592)

        EXTERNAL GKPPL,GKPPM,GKPFA,GKPTX

        INTEGER ERRIND,TNR,STYLE,COLI,PATTRN,FONT,PREC,LTYPE
        REAL YRES,WIDTH,SIZE,FACTOR,X,Y,ANGLE
        LOGICAL EMPTY,INIT
        INTEGER LPAGE,PAGES

        CHARACTER PAGE*8

        SAVE

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

C *     Workstation State List
C               connection and type
        INTEGER CONID, WTYPE
C               workstation state
        INTEGER STATE
C               workstation transformation
        REAL WINDOW(4),VIEWPT(4)
        INTEGER HEIGHT


        GOTO (999,999,  2,  3,  4,  5,  6,999,  8,999,
     *        999,999, 12, 13, 14, 15, 16,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999, 48,999,
     *        999,999,999,999, 54, 55,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999) FCTID+1

        IF (LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0)
     *      STOP 'GKS: internal inconsistency'
        GOTO 999

C*  open workstation
    2   CONTINUE
        CONID = IA(2)
        WTYPE = IA(3)

C*  initialize color table
        CALL GKPICT
C*  set up connection
        CALL GKPSCO (CONID,WTYPE)

C*  set default workstation window
        WINDOW(1) = 0.
        WINDOW(2) = 1.
        WINDOW(3) = 0.
        WINDOW(4) = 1.

C*  set default workstation viewport
        VIEWPT(1) = 0.
        IF (WTYPE .EQ. 92) THEN
          VIEWPT(2) = 0.2032
        ELSE
          VIEWPT(2) = 0.19685
        END IF
        VIEWPT(3) = 0.
        VIEWPT(4) = VIEWPT(2)

C*  set up device transformation
        CALL GKPSDT (WINDOW,VIEWPT,HEIGHT,WTYPE)

        PAGES = 0
        INIT = .FALSE.
        EMPTY = .TRUE.
        GOTO 999

C*  close workstation
    3   CONTINUE
        IF (INIT) THEN
          IF (.NOT.EMPTY) THEN
            IF (WTYPE .LT. 63) THEN
              CALL GKPPB ('showpage')
            END IF
          END IF
          CALL GKPPB ('psl restore end % GLI_GKS_dict')
          CALL GKPEP (PAGES)
          CALL GKPPB ('%%Trailer')
          CALL GKPPB ('GLI_GKS_save restore')
        END IF
        IF (PAGES .GT. 0) THEN
          CALL GDEC (PAGES,LPAGE,PAGE)
          CALL GKPPB ('%%Pages: '//PAGE(1:LPAGE))
        ELSE IF (WTYPE .LT. 63) THEN
          CALL GKPHDR
          CALL GKPPB ('%%Trailer')
          CALL GKPPB ('%%Pages: (none)')
        END IF
        CALL GKPU
        IF (WTYPE .GE. 63) CALL DPSCL
        GOTO 999

C*  activate workstation
    4   CONTINUE
        STATE = GACTIV
        GOTO 999

C*  deactivate workstation
    5   CONTINUE
        STATE = GINACT
        GOTO 999

C*  clear workstation
    6   CONTINUE
        IF (INIT) THEN
          IF (.NOT.EMPTY) THEN
            IF (WTYPE .LT. 63) THEN
              CALL GKPPB ('showpage')
            END IF
            EMPTY = .TRUE.
          END IF
          CALL GKPPB ('psl restore end % GLI_GKS_dict')
          CALL GKPEP (PAGES)
          INIT = .FALSE.
        END IF
        GOTO 999

C*  update workstation
    8   CONTINUE
        CALL GKPU
        GOTO 999

C*  polyline
   12   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (.NOT.INIT) THEN
            CALL GKPINI (WTYPE,WINDOW,VIEWPT,HEIGHT,PAGES)
            INIT = .TRUE.
          END IF
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current line type
          CALL GQLN (ERRIND,LTYPE)
          IF (LTYPE .NE. GLSOLI) CALL GKPSLN (LTYPE)
C*  inquire linewidth scale factor
          CALL GQLWSC (ERRIND,WIDTH)
          CALL GKPSLW (WIDTH)
C*  inquire polyline color index
          CALL GQPLCI (ERRIND,COLI)
          CALL GKPSCI (COLI,WTYPE)
C*  polyline
          CALL GKPPL (IA(1),R1,R2,LTYPE,TNR)
          IF (LTYPE .NE. GLSOLI) CALL GKPSLN (GLSOLI)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  polymarker
   13   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (.NOT.INIT) THEN
            CALL GKPINI (WTYPE,WINDOW,VIEWPT,HEIGHT,PAGES)
            INIT = .TRUE.
          END IF
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire marker size scale factor
          CALL GQMKSC (ERRIND,SIZE)
          CALL GKPSMS (SIZE)
C*  set marker rotation angle
          X = 0.0
          Y = 1.0
          CALL GCST (X,Y)
          ANGLE = -ATAN2(X,Y)*180.0/PI
          CALL GKPSMA (ANGLE)
C*  adjust linewidth
          FACTOR = SIZE*1.5
          CALL GKPSLW (FACTOR)
C*  inquire polymarker color index
          CALL GQPMCI (ERRIND,COLI)
          CALL GKPSFG (COLI,WTYPE)
          CALL GSIMPM (IA(1),R1,R2,GKPPM)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  text
   14   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (.NOT.INIT) THEN
            CALL GKPINI (WTYPE,WINDOW,VIEWPT,HEIGHT,PAGES)
            INIT = .TRUE.
          END IF
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire text font and text color index
          CALL GQTXFP (ERRIND,FONT,PREC)
          IF (PREC .NE. GSTRKP) THEN
            CALL GKPSCH (FONT,HEIGHT)
          ELSE
            CALL GKPSLW (1.0)
          END IF
C*  inquire text color index
          CALL GQTXCI (ERRIND,COLI)
          CALL GKPSCI (COLI,WTYPE)
          CALL GKTEXT (R1,R2,CHARS,GKPPL,GKPFA,GKPTX)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  fill area
   15   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (.NOT.INIT) THEN
            CALL GKPINI (WTYPE,WINDOW,VIEWPT,HEIGHT,PAGES)
            INIT = .TRUE.
          END IF
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current fill area interior style
          CALL GQFAIS (ERRIND,STYLE)
C*  inquire fill area color index
          CALL GQFACI (ERRIND,COLI)
          CALL GKPSCI (COLI,WTYPE)
          CALL GKPSLW (1.0)
          IF (STYLE .EQ. GSOLID) THEN
            CALL GKPFA (IA(1),R1,R2,TNR)
          ELSE IF (STYLE .EQ. GPATTR) THEN
            CALL GQFASI (ERRIND,PATTRN)
            CALL GKPPA (IA(1),R1,R2,TNR,PATTRN)
          ELSE
            YRES = 1.0/4650.0
            CALL GFILLA (IA(1),R1,R2,TNR,GKPPL,YRES)
          END IF
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  cell array
   16   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (.NOT.INIT) THEN
            CALL GKPINI (WTYPE,WINDOW,VIEWPT,HEIGHT,PAGES)
            INIT = .TRUE.
          END IF
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GKPCA (R1(1),R1(2),R2(1),R2(2),DX,DY,DIMX,IA,WTYPE)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  set color representation
   48   CONTINUE
        CALL GKPSCR (IA(2),R1(1),R1(2),R1(3))
        GOTO 999

C*  set workstation window
   54   CONTINUE
        WINDOW(1) = R1(1)
        WINDOW(2) = R1(2)
        WINDOW(3) = R2(1)
        WINDOW(4) = R2(2)
        CALL GKPSDT (WINDOW,VIEWPT,HEIGHT,WTYPE)

        IF (INIT) CALL GKPCLP (WINDOW)
        GOTO 999

C*  set workstation viewport
   55   CONTINUE
        VIEWPT(1) = R1(1)
        VIEWPT(2) = R1(2)
        VIEWPT(3) = R2(1)
        VIEWPT(4) = R2(2)
        CALL GKPSDT (WINDOW,VIEWPT,HEIGHT,WTYPE)
        GOTO 999

  999   RETURN
        END


        SUBROUTINE GKPPL (N,PX,PY,LTYPE,TNR)
C*  polyline

        INTEGER LIMIT
        PARAMETER (LIMIT=1000)

        INTEGER LTYPE, N
        REAL PX(N),PY(N)
        INTEGER TNR

        EXTERNAL GKPVEC,GKPCOD

        CALL GKPLIM (LIMIT)
C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,LTYPE,TNR,GKPVEC,GKPCOD)
C*  close stroke
        CALL GKPSK

        RETURN
        END


        SUBROUTINE GKPFA (N,PX,PY,TNR)
C*  fill area

        REAL PX(N),PY(N)
        INTEGER TNR
        INTEGER N

        EXTERNAL GKPVEC,GKPCOD

        CALL GKPLIM (0)
C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,0,TNR,GKPVEC,GKPCOD)
C*  fill
        CALL GKPPB ('fi')

        RETURN
        END


        SUBROUTINE GKPPA (N,PX,PY,TNR,PATTRN)
C*  pattern fill

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER TNR, PATTRN

        INTEGER PA(0:32),I,J,K,L
        CHARACTER STR*64

        EXTERNAL GKPVEC,GKPCOD

        CALL GKPLIM (0)
C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,0,TNR,GKPVEC,GKPCOD)

C*  inquire pattern
        CALL GKQPA (PATTRN,PA)
        K = PA(0)
        IF (K .EQ. 4 .OR. K .EQ. 8) THEN
            J = 2
        ELSE IF (K .EQ. 32) THEN
            J = 1
        ELSE
            STOP 'GKS: invalid pattern'
        END IF
        K = 1
        L = 1
        DO 1 I=1,32
          CALL GKPHEX (PA(K),STR(L:L+1))
          IF (MOD(I,J) .EQ. 0) K = K+1
          IF (K .GT. PA(0)) K = 1
          L = L+2
   1    CONTINUE

        CALL GKPPB ('bp <'//STR(1:64)//'> ep')

        RETURN
        END


        SUBROUTINE GKPTX (X,Y,CHARS)
C*  output hardware characters

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)

        REAL PI
        PARAMETER (PI = 3.141592)

        INTEGER I, J, IDUMMY
        REAL X,Y,YREL
        CHARACTER*(*) CHARS
        INTEGER FONT,HEIGHT

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

        INTEGER ERRIND,PREC,ALH,ALV
        REAL XORG,YORG
        CHARACTER STR*500

        CHARACTER*2 SHOW(GAHNOR:GARITE)
        REAL YFAC(GAVNOR:GABOTT)

        INTEGER CFONT, MAP(32)
        REAL CHH, CHGHT, CAPS(43), YSIZE
        CHARACTER*31 FONTS(43), PITCH, OCHAR
        INTEGER FAMILY, SIZE, FLEN, PLEN, IC 
        REAL UX, UY, PHI
        INTEGER TNR, ANGLE
        CHARACTER BS

        CHARACTER*40 AFONTS(223)
        INTEGER ACAPS(223), AENC(223)
        REAL ACAP

        SAVE

        DATA SHOW /'lj','lj','ct','rj'/
        DATA YFAC /0.,-1.2,-1.0,-0.5,0.,0.2/
        DATA YSIZE /64/

        DATA CFONT /0/, CHGHT /0.0/
        DATA MAP / 22,  9,  5, 14, 18, 26, 13,  1,
     *             24, 11,  7, 16, 20, 28, 13,  3,
     *             23, 10,  6, 15, 19, 27, 13,  2,
     *             25, 12,  8, 17, 21, 29, 13,  4/
        DATA CAPS / 0.662, 0.653, 0.676, 0.669,
     *              0.718, 0.718, 0.718, 0.718,
     *              0.562, 0.562, 0.562, 0.562,
     *              0.667,
     *              0.714, 0.714, 0.714, 0.714,
     *              0.722, 0.722, 0.722, 0.722,
     *              0.740, 0.740, 0.740, 0.740,
     *              0.732, 0.733, 0.732, 0.732,
     *              0.718, 0.718, 0.718, 0.718,
     *              0.681, 0.681, 0.681, 0.681,
     *              0.692, 0.692, 0.681, 0.681,
     *              0.587, 0.587/
        DATA FONTS /
     * 'Times-Roman', 'Times-Italic', 'Times-Bold', 'Times-BoldItalic',
     * 'Helvetica', 'Helvetica-Oblique', 'Helvetica-Bold',
     * 'Helvetica-BoldOblique', 'Courier', 'Courier-Oblique',
     * 'Courier-Bold', 'Courier-BoldOblique', 'Symbol',
     * 'LubalinGraph-Book', 'LubalinGraph-BookOblique',
     * 'LubalinGraph-Demi', 'LubalinGraph-DemiOblique',
     * 'NewCenturySchlbk-Roman', 'NewCenturySchlbk-Italic',
     * 'NewCenturySchlbk-Bold', 'NewCenturySchlbk-BoldItalic',
     * 'AvantGarde-Book', 'AvantGarde-BookOblique', 'AvantGarde-Demi',
     * 'AvantGarde-DemiOblique', 'Souvenir-Light',
     * 'Souvenir-LightItalic', 'Souvenir-Demi', 'Souvenir-DemiItalic',
     * 'Helvetica-Narrow', 'Helvetica-Narrow-Oblique',
     * 'Helvetica-Narrow-Bold', 'Helvetica-Narrow-BoldOblique',
     * 'Bookman-Light', 'Bookman-LightItalic', 'Bookman-Demi',
     * 'Bookman-DemiItalic', 'Palatino-Roman', 'Palatino-Italic',
     * 'Palatino-Bold', 'Palatino-BoldItalic',
     * 'ZapfChancery-MediumItalic', 'ZapfDingbats' /

        DATA (AFONTS(I), I = 1, 43) /
     * 'AGOldFace-BoldOutline', 'AGOldFace-Outline',
     * 'AmericanTypewriter-Bold', 'AmericanTypewriter-BoldA',
     * 'AmericanTypewriter-BoldCond', 'AmericanTypewriter-BoldCondA',
     * 'AmericanTypewriter-Cond', 'AmericanTypewriter-CondA',
     * 'AmericanTypewriter-Light', 'AmericanTypewriter-LightA',
     * 'AmericanTypewriter-LightCond', 'AmericanTypewriter-LightCondA',
     * 'AmericanTypewriter-Medium', 'AmericanTypewriter-MediumA',
     * 'AvantGarde-Book', 'AvantGarde-BookOblique', 'AvantGarde-Demi',
     * 'AvantGarde-DemiOblique', 'Bauhaus-Bold', 'Bauhaus-Demi',
     * 'Bauhaus-Heavy', 'Bauhaus-Light', 'Bauhaus-Medium', 'Bellevue',
     * 'Benguiat-Bold', 'Benguiat-Book', 'BenguiatGothic-Bold',
     * 'BenguiatGothic-BoldOblique', 'BenguiatGothic-Book',
     * 'BenguiatGothic-BookOblique', 'BenguiatGothic-Heavy',
     * 'BenguiatGothic-HeavyOblique', 'BenguiatGothic-Medium',
     * 'BenguiatGothic-MediumOblique', 'Bookman-Bold',
     * 'Bookman-BoldItalic', 'Bookman-Demi', 'Bookman-DemiItalic',
     * 'Bookman-Light', 'Bookman-LightItalic', 'Bookman-Medium',
     * 'Bookman-MediumItalic', 'CaslonTwoTwentyFour-Black' /
        DATA (AFONTS(I), I = 44, 92) /
     * 'CaslonTwoTwentyFour-BlackIt', 'CaslonTwoTwentyFour-Bold',
     * 'CaslonTwoTwentyFour-BoldIt', 'CaslonTwoTwentyFour-Book',
     * 'CaslonTwoTwentyFour-BookIt', 'CaslonTwoTwentyFour-Medium',
     * 'CaslonTwoTwentyFour-MediumIt', 'CastellarMT', 'Cheltenham-Bold',
     * 'Cheltenham-BoldItalic', 'Cheltenham-Book',
     * 'Cheltenham-BookItalic', 'Cheltenham-Light',
     * 'Cheltenham-LightItalic', 'Cheltenham-Ultra',
     * 'Cheltenham-UltraItalic', 'City-Bold', 'City-BoldItalic',
     * 'City-Medium', 'City-MediumItalic', 'Courier', 'Courier-Bold',
     * 'Courier-BoldOblique', 'Courier-Oblique', 'Cushing-Bold',
     * 'Cushing-BoldItalic', 'Cushing-Book', 'Cushing-BookItalic',
     * 'Cushing-Heavy', 'Cushing-HeavyItalic', 'Cushing-Medium',
     * 'Cushing-MediumItalic', 'DorchesterScriptMT', 'Esprit-Black',
     * 'Esprit-BlackItalic', 'Esprit-Bold', 'Esprit-BoldItalic',
     * 'Esprit-Book', 'Esprit-BookItalic', 'Esprit-Medium',
     * 'Esprit-MediumItalic', 'Fenice-Bold', 'Fenice-BoldOblique',
     * 'Fenice-Light', 'Fenice-LightOblique', 'Fenice-Regular',
     * 'Fenice-RegularOblique', 'Fenice-Ultra', 'Fenice-UltraOblique' /
        DATA (AFONTS(I), I = 93, 136) /
     * 'FrizQuadrata', 'FrizQuadrata-Bold',
     * 'Galliard-Black', 'Galliard-BlackItalic', 'Galliard-Bold',
     * 'Galliard-BoldItalic', 'Galliard-Italic', 'Galliard-Roman',
     * 'Galliard-Ultra', 'Galliard-UltraItalic', 'Garamond-Bold',
     * 'Garamond-BoldCondensed', 'Garamond-BoldCondensedItalic',
     * 'Garamond-BoldItalic', 'Garamond-Book', 'Garamond-BookCondensed',
     * 'Garamond-BookCondensedItalic', 'Garamond-BookItalic',
     * 'Garamond-Light', 'Garamond-LightCondensed',
     * 'Garamond-LightCondensedItalic', 'Garamond-LightItalic',
     * 'Garamond-Ultra', 'Garamond-UltraCondensed',
     * 'Garamond-UltraCondensedItalic', 'Garamond-UltraItalic',
     * 'GillSans', 'GillSans-Bold', 'GillSans-BoldItalic',
     * 'GillSans-ExtraBold', 'GillSans-Italic', 'Giovanni-Black',
     * 'Giovanni-BlackItalic', 'Giovanni-Bold', 'Giovanni-BoldItalic',
     * 'Giovanni-Book', 'Giovanni-BookItalic', 'GoudyTextMT',
     * 'Helvetica', 'Helvetica-Bold', 'Helvetica-BoldOblique',
     * 'Helvetica-Oblique', 'Isadora-Bold', 'Isadora-Regular' /
        DATA (AFONTS(I), I = 137, 182) /
     * 'ItcEras-Bold', 'ItcEras-Book', 'ItcEras-Demi', 'ItcEras-Light',
     * 'ItcEras-Medium', 'ItcEras-Ultra', 'ItcKabel-Bold',
     * 'ItcKabel-Book', 'ItcKabel-Demi', 'ItcKabel-Medium',
     * 'ItcKabel-Ultra', 'Korinna-Bold',
     * 'Korinna-KursivBold', 'Korinna-KursivRegular', 'Korinna-Regular',
     * 'Leawood-Black', 'Leawood-BlackItalic', 'Leawood-Bold',
     * 'Leawood-BoldItalic', 'Leawood-Book', 'Leawood-BookItalic',
     * 'Leawood-Medium', 'Leawood-MediumItalic', 'Machine',
     * 'Machine-Bold', 'Madrone', 'NewBaskerville-Bold',
     * 'NewBaskerville-BoldItalic', 'NewBaskerville-Italic',
     * 'NewBaskerville-Roman', 'NuptialScript', 'OfficinaSans-Bold',
     * 'OfficinaSans-BoldItalic', 'OfficinaSans-Book',
     * 'OfficinaSans-BookItalic', 'OfficinaSerif-Bold',
     * 'OfficinaSerif-BoldItalic', 'OfficinaSerif-Book',
     * 'OfficinaSerif-BookItalic', 'PepitaMT', 'Ponderosa', 'Poplar',
     * 'Rosewood-Fill', 'Rosewood-Regular', 'RussellSquare',
     * 'RussellSquare-Oblique' /
        DATA (AFONTS(I), I = 183, 223) /
     * 'Slimbach-Black', 'Slimbach-BlackItalic', 'Slimbach-Bold',
     * 'Slimbach-BoldItalic', 'Slimbach-Book', 'Slimbach-BookItalic',
     * 'Slimbach-Medium', 'Slimbach-MediumItalic',
     * 'Souvenir-Demi', 'Souvenir-DemiItalic', 'Souvenir-Light',
     * 'Souvenir-LightItalic', 'Stencil', 'Tiepolo-Black',
     * 'Tiepolo-BlackItalic', 'Tiepolo-Bold', 'Tiepolo-BoldItalic',
     * 'Tiepolo-Book', 'Tiepolo-BookItalic', 'Times-Bold',
     * 'Times-BoldItalic', 'Times-Italic', 'Times-Roman',
     * 'Usherwood-Black', 'Usherwood-BlackItalic', 'Usherwood-Bold',
     * 'Usherwood-BoldItalic', 'Usherwood-Book', 'Usherwood-BookItalic',
     * 'Usherwood-Medium', 'Usherwood-MediumItalic', 'Veljovic-Black',
     * 'Veljovic-BlackItalic', 'Veljovic-Bold', 'Veljovic-BoldItalic',
     * 'Veljovic-Book', 'Veljovic-BookItalic', 'Veljovic-Medium',
     * 'Veljovic-MediumItalic', 'Willow', 'ZapfChancery-MediumItalic' /

        DATA (ACAPS(I), I = 1, 100) /
     * 706, 706, 686, 686, 683, 683, 688, 688, 683, 683,
     * 683, 683, 686, 686, 740, 740, 740, 740, 700, 700,
     * 700, 700, 700, 706, 658, 658, 691, 691, 691, 691,
     * 691, 691, 691, 691, 681, 681, 681, 681, 681, 681,
     * 681, 681, 680, 680, 680, 680, 680, 680, 680, 680,
     * 715, 703, 703, 703, 703, 703, 703, 703, 703, 706,
     * 706, 706, 706, 562, 562, 562, 562, 714, 714, 714,
     * 714, 714, 714, 714, 714, 589, 666, 666, 666, 666,
     * 666, 666, 666, 666, 692, 692, 692, 692, 692, 692,
     * 692, 692, 658, 658, 682, 682, 680, 680, 680, 680 /
        DATA (ACAPS(I), I = 101, 223) /
     * 682, 682, 623, 630, 630, 623, 623, 630, 630, 627,
     * 623, 630, 630, 624, 622, 630, 630, 623, 682, 682,
     * 682, 682, 682, 660, 660, 660, 660, 660, 660, 739,
     * 718, 718, 718, 718, 695, 700, 667, 667, 667, 667,
     * 667, 667, 702, 702, 702, 702, 702, 714, 714, 714,
     * 714, 709, 709, 709, 709, 709, 709, 709, 709, 717,
     * 717, 750, 660, 660, 660, 660, 604, 685, 685, 685,
     * 685, 685, 685, 685, 685, 567, 850, 750, 630, 648,
     * 720, 720, 670, 670, 670, 670, 670, 670, 670, 670,
     * 732, 732, 732, 732, 748, 614, 614, 614, 614, 614,
     * 614, 676, 669, 653, 662, 627, 627, 627, 627, 627,
     * 627, 627, 627, 626, 626, 626, 626, 626, 626, 626,
     * 626, 750, 708 /

        DATA AENC /
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     * 0, 0, 0 /

        BS = CHAR(92)

C*  device transformation
        CALL GKPDT (X,Y,XORG,YORG)
C*  inquire text font and precision, text alignment
        CALL GQTXFP (ERRIND,IDUMMY,PREC)
        CALL GQTXAL (ERRIND,ALH,ALV)

C*  inquire character-up vector
        CALL GQCHUP (ERRIND,UX,UY)
        CALL GQCNTN (ERRIND,TNR)
        CALL GCNT (UX,UY,TNR)
C*  segment transformation
        CALL GCST (UX,UY)
        ANGLE = NINT(-ATAN2(UX,UY)*180.0/PI)

        IF (PREC .EQ. GSTRP) THEN
C*  align the text
          PHI = ANGLE/180.0*PI
          YREL = YSIZE*YFAC(ALV)
          XORG = XORG-YREL*SIN(PHI)
          YORG = YORG+YREL*COS(PHI)
        END IF

        IF (ANGLE .EQ. 0) THEN
          CALL GKPMOV (XORG,YORG)
        ELSE
          CALL GKPAMT (ANGLE,XORG,YORG)
        END IF

        J = 0
        DO 1 I = 1,LEN(CHARS)
          IC = ICHAR(CHARS(I:I))
          IF (IC. LT. 0) IC = 255+IC
          IF (IC .LE. 127) THEN
            IF (INDEX('()\\',CHARS(I:I)) .GT. 0) THEN
              J = J+1
              STR(J:J) = BS
            END IF
            J = J+1
            STR(J:J) = CHARS(I:I)
          ELSE
            CALL GOCT (IC,OCHAR)
            J = J+1
            STR(J:J) = BS
            J = J+3
            STR(J-2:J) = OCHAR(1:3)
          END IF
   1    CONTINUE

        CALL GKPPB ('('//STR(1:J)//') '//SHOW(ALH))
        IF (ANGLE .NE. 0) CALL GKPPB ('gr')

        RETURN


        ENTRY GKPSCH (FONT,HEIGHT)
C*  set character height

        CALL GSCT
        CALL GCHH (CHH)

        IF (FONT .NE. CFONT .OR. (ABS(CHH-CHGHT) .GT. FEPS)) THEN  
            CFONT = ABS(FONT)
            CHGHT = ABS(CHH)

C*  fonts used by Adobe's Illustrator (contributed by Mika.Heiskanen@fmi.fi)

            IF (CFONT .GE. 10001 .AND. CFONT .LE. 10223) THEN

                FAMILY = CFONT-10000
                FLEN = INDEX(AFONTS(FAMILY),' ') - 1

                YSIZE = CHGHT*HEIGHT
                ACAP = FLOAT(ACAPS(FAMILY))/1000.0
                SIZE = MIN(MAX(INT(YSIZE/ACAP),1),7200)
                CALL GDEC (SIZE,PLEN,PITCH)

                IF (AENC(FAMILY) .EQ. 0) THEN
                    CALL GKPPB ('gsave /'//AFONTS(FAMILY)(1:FLEN)//
     *                  '_ ISOLatin1Encoding')
                    CALL GKPPB ('/'//AFONTS(FAMILY)(1:FLEN)//
     *                  ' encodefont pop grestore')

                    CALL GKPPB ('/'//AFONTS(FAMILY)(1:FLEN)//'_'//
     *                  ' findfont '//PITCH(1:PLEN)//
     *                  ' scalefont setfont')
                ELSE
                    CALL GKPPB ('/'//AFONTS(FAMILY)(1:FLEN)//
     *                  ' findfont '//PITCH(1:PLEN)//
     *                  ' scalefont setfont')
                END IF

            ELSE

                IF (CFONT .GE. 101 .AND. CFONT .LE. 143) THEN
                    FAMILY = CFONT-100
                ELSE IF (CFONT .GE. 1 .AND. CFONT .LE. 32) THEN
                    FAMILY = MAP(CFONT)
                ELSE
                    FAMILY = 9
                END IF
                FLEN = INDEX(FONTS(FAMILY),' ') - 1
               
                YSIZE = CHGHT*HEIGHT
                SIZE = MIN(MAX(INT(YSIZE/CAPS(FAMILY)),1),7200)
                CALL GDEC (SIZE,PLEN,PITCH)

                IF (FAMILY .NE. 13 .AND. FAMILY .NE. 42 .AND.
     *              FAMILY .NE. 43) THEN
                    CALL GKPPB ('gsave /'//FONTS(FAMILY)(1:FLEN)//
     *                  '_ ISOLatin1Encoding')
       
                    CALL GKPPB ('/'//FONTS(FAMILY)(1:FLEN)//
     *                  ' encodefont pop grestore')

                    CALL GKPPB ('/'//FONTS(FAMILY)(1:FLEN)//'_'//
     *                  ' findfont '//PITCH(1:PLEN)//
     *                  ' scalefont setfont')
                ELSE
                    CALL GKPPB ('/'//FONTS(FAMILY)(1:FLEN)//
     *                  ' findfont '//PITCH(1:PLEN)//
     *                  ' scalefont setfont')
                END IF

            END IF

        END IF

        RETURN
        END


        SUBROUTINE GKPCA (XMIN,XMAX,YMIN,YMAX,DX,DY,DIMX,COLIA,WTYPE)
C*  cell array

        REAL XMIN, XMAX, YMIN, YMAX
        INTEGER DX, DY, DIMX, COLIA(*), WTYPE

        INTEGER ERRIND, TNR, CLSW
        REAL CLRT(4)
        REAL X1, Y1, X2, Y2
        INTEGER W, H, X, Y

        CHARACTER*8 SDX, SDY, SX, SY, SW, SH
        INTEGER LSDX, LSDY, LSX, LSY, LSW, LSH
        CHARACTER LINE*78, SCI*6
        INTEGER I, J, II, JJ, CI, LSCI, SWAP

C*  inquire current normalization transformation number
        CALL GQCNTN (ERRIND,TNR)

        CALL GDEC (DX,LSDX,SDX)
        CALL GDEC (DY,LSDY,SDY)

        X1 = XMIN
        Y1 = YMAX
        CALL GNT (X1,Y1,TNR)
        CALL GST (X1,Y1)
        CALL GKPDT (X1,Y1,X1,Y1)

        X2 = XMAX
        Y2 = YMIN
        CALL GNT (X2,Y2,TNR)
        CALL GST (X2,Y2)
        CALL GKPDT (X2,Y2,X2,Y2)

        W = NINT(ABS(X2-X1))
        H = NINT(ABS(Y2-Y1))
        X = NINT(MIN(X1,X2))
        Y = NINT(MIN(Y1,Y2))

        CALL GDEC (X,LSX,SX)
        CALL GDEC (Y,LSY,SY)
        CALL GDEC (W,LSW,SW)
        CALL GDEC (H,LSH,SH)

        CALL GKPPB ('gsave')

C*  inquire current clipping rectangle and set clipping path
        CALL GQCLIP (ERRIND,CLSW,CLRT)
        CALL GKPCLP (CLRT)

        IF (MOD(WTYPE,2) .EQ. 0) THEN
            CALL GKPPB ('/rgbstr '//SDX(1:LSDX)//' string def')
        ELSE
            CALL GKPPB ('/gstr '//SDX(1:LSDX)//' string def')
        END IF
        CALL GKPPB (SX(1:LSX)//' '//SY(1:LSY)//' translate')
        CALL GKPPB (SW(1:LSW)//' '//SH(1:LSH)//' scale')

        SWAP = 0
        IF (X1 .GT. X2) SWAP = 1
        IF (Y1 .GT. Y2) SWAP = SWAP + 2

        IF (SWAP .EQ. 0) THEN
            CALL GKPPB (SDX(1:LSDX)//' '//SDY(1:LSDY)//' 8 ['//
     *          SDX(1:LSDX)//' 0 0 -'//SDY(1:LSDY)//' 0 '//
     *          SDY(1:LSDY)//']')
        ELSE IF (SWAP .EQ. 1) THEN
            CALL GKPPB (SDX(1:LSDX)//' '//SDY(1:LSDY)//' 8 ['//
     *          '-'//SDX(1:LSDX)//' 0 0 -'//SDY(1:LSDY)//' '//
     *          SDX(1:LSDX)//' '//SDY(1:LSDY)//']')
        ELSE IF (SWAP .EQ. 2) THEN
            CALL GKPPB (SDX(1:LSDX)//' '//SDY(1:LSDY)//' 8 ['//
     *          SDX(1:LSDX)//' 0 0 '//SDY(1:LSDY)//' 0 0]')
        ELSE
            CALL GKPPB (SDX(1:LSDX)//' '//SDY(1:LSDY)//' 8 ['//
     *          '-'//SDX(1:LSDX)//' 0 '//SDX(1:LSDX)//' '//
     *          SDY(1:LSDY)//' 0 0]')
        END IF

        IF (MOD(WTYPE,2) .EQ. 0) THEN
            CALL GKPPB ('{currentfile rgbstr readhexstring pop}')
            CALL GKPPB ('false 3 colorimage')
        ELSE
            CALL GKPPB ('{currentfile gstr readhexstring pop}')
            CALL GKPPB ('image')
        END IF

        II = 0
        JJ = 0
        DO 1 J = 1,DY
            DO 2 I = 1,DX
                CI = COLIA(JJ+I)
                IF (CI .GE. 588) THEN
                    CI = 80 + (CI-588)/56 * 12 +
     *                  NINT(MOD(CI-588,56) * 11.0/56.0)
                ELSE IF (CI .GE. 257) THEN
                    CI = 8 + NINT((CI-257)/330.0 * (72-1))
                ELSE IF (CI .LT. 0) THEN
                    CI = 0
                END IF
                CALL GKPQCI (CI, LSCI, SCI, WTYPE)
                II = II+LSCI
                LINE(II-LSCI+1:II) = SCI(1:LSCI)
                IF (II .GE. 78) THEN
                    CALL GKPPB (LINE(1:II))
                    II = 0
                END IF
   2        CONTINUE
            JJ = JJ+DIMX
   1    CONTINUE
        IF (II .GT. 0) CALL GKPPB (LINE(1:II))

        CALL GKPPB ('grestore')

        RETURN
        END


        SUBROUTINE GKPPM (X,Y,MARKER)
C*  output a marker symbol

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

        INTEGER IX, IY, MARKER
        REAL DX, DY, X, Y
        CHARACTER*8 XPOS,YPOS
        INTEGER LXPOS,LYPOS

        CHARACTER*3 MACRO(-20:5)
        DATA MACRO /
     *    'nom','npl','ftr','ftl','tud','fst',' st','fdm','ndm','fhg',
     *    'nhg','fbt','nbt','fsq','nsq','ftd','ntd','ftu','ntu','fci',
     *    ' dt',
     *    ' dt',' pl','fas','nci',' dc'/

C*  device transformation
        CALL GKPDT (X,Y,DX,DY)

        IX = NINT(DX)
        IY = NINT(DY)
        CALL GDEC (IX,LXPOS,XPOS)
        CALL GDEC (IY,LYPOS,YPOS)

        CALL GKPPB (XPOS(1:LXPOS)//' '//YPOS(1:LYPOS)//
     *    ' '//MACRO(MARKER))

        RETURN
        END


        SUBROUTINE GKPCLP (CLRT)
C*  set clipping path

        REAL CLRT(4)

        INTEGER I, J
        REAL CX1, CY1, CX2, CY2
        INTEGER IX1, IX2, IY1, IY2
        CHARACTER*8 X1, X2, Y1, Y2
        INTEGER LX1, LX2, LY1, LY2

        IF (CLRT(1) .LT. CLRT(2)) THEN
            I = 1
        ELSE
            I = 2
        END IF
        IF (CLRT(3) .LT. CLRT(4)) THEN
            J = 3
        ELSE
            J = 4
        END IF

        CALL GKPDTR (CLRT(  I),CLRT(  J),CX1,CY1)
        CALL GKPDTR (CLRT(3-I),CLRT(7-J),CX2,CY2)

        IX1 = INT(CX1)-2
        IY1 = INT(CY1)-2
        IX2 = NINT(CX2)+2
        IY2 = NINT(CY2)+2

        CALL GDEC (IX1,LX1,X1)
        CALL GDEC (IY1,LY1,Y1)
        CALL GDEC (IX2,LX2,X2)
        CALL GDEC (IY2,LY2,Y2)

        CALL GKPPB ('np '//X1(1:LX1)//' '//Y1(1:LY1)//' m '//
     *    X1(1:LX1)//' '//Y2(1:LY2)//' l '//
     *    X2(1:LX2)//' '//Y2(1:LY2)//' l '//
     *    X2(1:LX2)//' '//Y1(1:LY1)//' l cp clip')

        RETURN
        END


        SUBROUTINE GKPSDT (WN,VP,HEIGHT,WTYPE)
C*  set up device transformation

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)

        REAL WN(4),VP(4)
        INTEGER HEIGHT, WTYPE
        REAL MAGS, SIZEX, SIZEY
        LOGICAL LANDSC
        INTEGER YTRANS, LIMIT, ILIMIT

        INTEGER IX, IX1, IX2, IY, IY1, IY2, JX, JY, RX, RY
        INTEGER LX1, LX2, LY1, LY2, IANGLE, NP
        REAL A, B, C, D, E, F, G, H, X, XD, Y, YD
        REAL MW, MH, MAGN

        CHARACTER*8 ANGLE, X1, X2, Y1, Y2
        INTEGER LANGLE
        LOGICAL STROKE
        CHARACTER*1 DC(-1:1,-1:1)

        SAVE

        INCLUDE 'gksdefs.i'

        DATA DC /'F', 'G', 'H', 'E', 'I', 'A', 'D', 'C', 'B'/

        E = (VP(2)-VP(1))/(WN(2)-WN(1))
        IF (WTYPE .EQ. 92) THEN
C*  11in, 180dpi
          F = (1980-1)/0.2794*600/180
        ELSE
C*  11.25in, 600pi
          F = (6750-1)/0.28575
        END IF
        G = (VP(4)-VP(3))/(WN(4)-WN(3))
        IF (WTYPE .EQ. 92) THEN
C*  8in, 180dpi
          H = (1440-1)/0.2032*600/180
        ELSE
C*  7.75in, 600dpi
          H = (4650-1)/0.19685
        END IF

        A = E*F
        B = F*(VP(1)-WN(1)*E)
        C = G*H
        D = H*(VP(3)-WN(3)*G)

        MW = A*(WN(2)-WN(1))
        MH = C*(WN(4)-WN(3))
        HEIGHT = INT(C)

        STROKE = .FALSE.

        RETURN


        ENTRY GKPSIZ (LANDSC, MAGS, SIZEX, SIZEY)
C*  return size information

        IF (ABS(MAGS) .GT. FEPS) THEN
          MAGN = 1.2**MAGS
        ELSE
          MAGN = 1.0
        END IF

        IF (LANDSC) THEN
          SIZEX = MH*MAGN
          SIZEY = MW*MAGN
        ELSE
          SIZEX = MW*MAGN
          SIZEY = MH*MAGN
        END IF

        RETURN


        ENTRY GKPBND (WTYPE, LANDSC, MAGS, YTRANS)
C*  display bounding box

        IF (ABS(MAGS) .GT. FEPS) THEN
          MAGN = 1.2**MAGS
        ELSE
          MAGN = 1.0
        END IF

        IF (WTYPE .LT. 63) THEN
          IF (MOD(WTYPE,2) .EQ. 0) THEN
            IX1 = 21
          ELSE
            IX1 = 18
          END IF
        ELSE
          IX1 = 0
        END IF
        IF (WTYPE .LT. 63) THEN
          IY1 = 15
        ELSE
          IY1 = 0
        END IF
        IF (LANDSC) THEN
          IX2 = IX1+NINT(MH*72/600*MAGN)
          IY2 = IY1+NINT(MW*72/600*MAGN)
          YTRANS = IY2
        ELSE
          IX2 = IX1+NINT(MW*72/600*MAGN)
          IY2 = IY1+NINT(MH*72/600*MAGN)
          YTRANS = IY1
        END IF
        CALL GDEC (IX1,LX1,X1)
        CALL GDEC (IY1,LY1,Y1)
        CALL GDEC (IX2,LX2,X2)
        CALL GDEC (IY2,LY2,Y2)

        CALL GKPPB ('%%BoundingBox: '//
     *    X1(1:LX1)//' '//Y1(1:LY1)//' '//
     *    X2(1:LX2)//' '//Y2(1:LY2))

        RETURN


        ENTRY GKPVEC (X,Y)
C*  re-initialize vector drawing sequence

C*  device transformation
        IX = NINT(A*X+B)
        IY = NINT(C*Y+D)

        IF (STROKE) THEN
          CALL GKPPB ('sk')
          STROKE = .FALSE.
        END IF

        CALL GDEC (IX,LX1,X1)
        CALL GDEC (IY,LY1,Y1)
        CALL GKPPB ('np '//X1(1:LX1)//' '//Y1(1:LY1)//' m')
        NP = 1

        RETURN


        ENTRY GKPCOD (X,Y)

        JX = IX
        JY = IY
C*  device transformation
        IX = NINT(A*X+B)
        IY = NINT(C*Y+D)

        IF (NP .EQ. 1 .OR. IX .NE. JX .OR. IY .NE. JY) THEN
          RX = IX-JX
          RY = IY-JY
          IF (ABS(RX) .GT. 1 .OR. ABS(RY) .GT. 1) THEN
            CALL GDEC (RX,LX1,X1)
            CALL GDEC (RY,LY1,Y1)
            CALL GKPPB (X1(1:LX1)//' '//Y1(1:LY1)//' rl')
          ELSE
            CALL GKPPB (DC(RX,RY))
          END IF
          NP = NP+1

          IF (NP .EQ. LIMIT) THEN
            CALL GKPPB ('sk')
            STROKE = .FALSE.

            CALL GDEC (IX,LX1,X1)
            CALL GDEC (IY,LY1,Y1)
            CALL GKPPB (X1(1:LX1)//' '//Y1(1:LY1)//' m')
            NP = 1
          ELSE
            STROKE = .TRUE.
          END IF
        END IF

        RETURN


        ENTRY GKPSK
C*  stroke

        IF (STROKE) THEN
          CALL GKPPB ('sk')
          STROKE = .FALSE.
        END IF

        RETURN


        ENTRY GKPMOV (X,Y)
C*  move to position

        IX = NINT(X)
        IY = NINT(Y)

        CALL GDEC (IX,LX1,X1)
        CALL GDEC (IY,LY1,Y1)
        CALL GKPPB (X1(1:LX1)//' '//Y1(1:LY1)//' m')

        RETURN


        ENTRY GKPAMT (IANGLE,X,Y)
C*  rotated move

        IX = NINT(X)
        IY = NINT(Y)

        CALL GDEC (IANGLE,LANGLE,ANGLE)
        CALL GDEC (IX,LX1,X1)
        CALL GDEC (IY,LY1,Y1)
        CALL GKPPB (ANGLE(1:LANGLE)//' '//
     *    X1(1:LX1)//' '//Y1(1:LY1)//' am')

        RETURN


        ENTRY GKPDT (X,Y,XD,YD)

C*  device transformation
        XD = A*X+B
        YD = C*Y+D

        RETURN


        ENTRY GKPDTR (X,Y,XD,YD)

C*  device transformation
        XD = A*X+B
        YD = C*Y+D

        RETURN


        ENTRY GKPLIM (ILIMIT)
C*  point limit

        LIMIT = ILIMIT

        RETURN
        END


        SUBROUTINE GKPSLN (LTYPE)
C*  set linetype

        INTEGER LTYPE, N
        CHARACTER*30 PATTRN(-30:4)
        DATA PATTRN /
     *  '[20 15 20 15 20 15 20 30]', '[20 15 20 15 20 30]',
     *  '[20 15 20 30]', '[15 10 15 10 15 10 15 30]',
     *  '[15 10 15 10 15 30]', '[15 10 15 30]',
     *  '[15 10 15 10 15 10 15 20]', '[15 10 15 10 15 20]',
     *  '[15 10 15 20]', '[0 10]', '[0 20]', '[0 40]', '[0 50]',
     *  '[0 15 0 15 0 30]', '[0 15 0 30]', '[30 15 0 15 0 15 0 15]',
     *  '[30 15 0 15 0 15]', '[30 15 0 15]', '[60 20 30 20]',
     *  '[60 20]', '[30 30]', '[30 20]',
     *  '[0 30 0 30 0 60]', '[0 30 0 60]', '[0 60]', '[60 60]',
     *  '[120 40 60 40]', '[120 40]', '[60 30 0 30 0 30 0 30]',
     *  '[60 30 0 30 0 30]', '[]',
     *  '[]', '[60 40]', '[0 30]', '[60 30 0 30]' /

        N = INDEX(PATTRN(LTYPE),']')
        CALL GKPPB (PATTRN(LTYPE)(1:N)//' 0 setdash')

        IF (LTYPE .EQ. 1) THEN
          CALL GKPPB ('0 setlinecap')
        ELSE
          CALL GKPPB ('1 setlinecap')
        END IF

        RETURN
        END


        SUBROUTINE GKPSLW (WIDTH)
C*  set linewidth

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)

        INTEGER N
        REAL WIDTH, CWIDTH
        CHARACTER*12 LWIDTH

        DATA CWIDTH /0.0/

        IF (ABS(WIDTH-CWIDTH) .GT. FEPS) THEN 
            CWIDTH = ABS(WIDTH)
            CALL GFLT (4*CWIDTH,N,LWIDTH)
            CALL GKPPB (LWIDTH(1:N)//' lw')
        END IF

        RETURN
        END


        SUBROUTINE GKPSMS (SIZE)
C*  set marker size scale factor

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)
        
        INTEGER N
        REAL SIZE,CSIZE
        CHARACTER*12 LSIZE

        DATA CSIZE /0.0/

        IF (ABS(SIZE-CSIZE) .GT. FEPS) THEN 
            CSIZE = ABS(SIZE)
            CALL GFLT (CSIZE,N,LSIZE)
            CALL GKPPB (LSIZE(1:N)//' ms')
        END IF

        RETURN
        END


        SUBROUTINE GKPSMA (ANGLE)
C*  set marker rotation angle

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)
        
        INTEGER N
        REAL ANGLE, CANGLE
        CHARACTER*12 LANGLE

        DATA CANGLE /0.0/

        IF (ABS(ANGLE-CANGLE) .GT. FEPS) THEN 
            CANGLE = ABS(ANGLE)
            CALL GFLT (CANGLE,N,LANGLE)
            CALL GKPPB (LANGLE(1:N)//' ma')
        END IF

        RETURN
        END


        SUBROUTINE GKPHDR

        INTEGER NCHARS
        CHARACTER BUFFER*100

        CALL GKINFO (NCHARS, BUFFER)

        CALL GKPPB ('%!PS-Adobe-2.0')
        IF (NCHARS .GT. 0) THEN
            CALL GKPPB ('%%Creator: '//BUFFER(36:NCHARS)//
     *          ', GLI GKS PostScript Device Handler, V4.5')
            CALL GKPPB ('%%+CreationDate: '//BUFFER(1:24))
        ELSE
            CALL GKPPB (
     *          '%%Creator: GLI GKS PostScript Device Handler, V4.5')
        END IF
        CALL GKPPB ('%%+Copyright @ 1993-1995, J.Heinen')
        CALL GKPPB ('%%Pages: (at end)')

        RETURN
        END


        SUBROUTINE GKPINI (WTYPE,WN,VP,HEIGHT,PAGES)
C*  initialize - write PostScript prolog

        INTEGER WTYPE
        REAL WN(4), VP(4)
        INTEGER HEIGHT, PAGES

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)

        INTEGER LPAGE, LMAG, LY
        CHARACTER PAGE*8, MAG*12, X*2, Y*8
        LOGICAL LANDSC
        REAL MAGS, SIZEX, SIZEY
        INTEGER FORM, DPI, YTRANS

        SAVE

C*  GKS symbol definitions
        INCLUDE 'gksdefs.i'
C*  GKS description table
        INCLUDE 'gksdescr.i'
C*  GKS state list
        INCLUDE 'gksstate.i'

        IF (WTYPE .EQ. 92) THEN
          LANDSC = (VP(2)-VP(1) .GT. 0.2032)
        ELSE
          LANDSC = (VP(2)-VP(1) .GT. 0.19685)
        END IF

        IF (PAGES .EQ. 0) THEN
        
        CALL GKMAGS (MAGS, DPI)
        IF (WTYPE .EQ. 92) THEN
          MAGS = ALOG(180.0 / DPI) / ALOG(1.2)
        END IF

        IF (WTYPE .GE. 63) THEN
          CALL GKPSIZ (LANDSC, MAGS, SIZEX, SIZEY)
          IF (WTYPE .EQ. 92) THEN
C*  SIXEL format
            FORM = 0
          ELSE
C*  GIF format
            FORM = 1
          END IF
          CALL DPSOP (SIZEX, SIZEY, FORM)
        END IF

        CALL GKPHDR
        CALL GKPBND (WTYPE, LANDSC, MAGS, YTRANS)
        CALL GKPPB ('%%EndComments')
        CALL GKPPB ('%%BeginProcSet GLI 4.5')

        CALL GKPPB ('save /GLI_GKS_save exch def')
        CALL GKPPB ('/GLI_GKS_dict 150 dict def GLI_GKS_dict begin')

        CALL GKPPB ('/in {72 mul} def')
        CALL GKPPB ('/np {newpath} def')
        CALL GKPPB ('/cp {closepath} def')
        CALL GKPPB ('/m {moveto} def')
        CALL GKPPB ('/l {lineto} def')
        CALL GKPPB ('/A {1 0 rlineto} def')
        CALL GKPPB ('/B {1 1 rlineto} def')
        CALL GKPPB ('/C {0 1 rlineto} def')
        CALL GKPPB ('/D {-1 1 rlineto} def')
        CALL GKPPB ('/E {-1 0 rlineto} def')
        CALL GKPPB ('/F {-1 -1 rlineto} def')
        CALL GKPPB ('/G {0 -1 rlineto} def')
        CALL GKPPB ('/H {1 -1 rlineto} def')
        CALL GKPPB ('/I {1 0 rlineto -1 0 rlineto} def')
        CALL GKPPB ('/am {np gsave translate rotate 0 0 m} def')
        CALL GKPPB ('/gr {grestore} def')
        CALL GKPPB ('/rm {rmoveto} def')
        CALL GKPPB ('/srm {rxy rx s mul ry s mul rm} def')
        CALL GKPPB ('/rl {rlineto} def')
        CALL GKPPB ('/srl {rxy rx s mul ry s mul rl} def')
        CALL GKPPB ('/sk {stroke} def')
        CALL GKPPB ('/csk {closepath stroke} def')
        CALL GKPPB ('/fi {closepath eofill} def')
        CALL GKPPB ('/sg {setgray} def')
        CALL GKPPB ('/sc {setrgbcolor} def')
        IF (WTYPE .GE. 63 .AND. X11) THEN
          CALL GKPPB (
     *      '/sp {dup [/DevicePixel 8] setcolorspace pop setcolor} def')
        END IF
        CALL GKPPB ('/fg {0 sg} def')
        CALL GKPSBG (WTYPE)
        CALL GKPPB ('/lw {setlinewidth} def')
        CALL GKPPB ('/ms {/s exch def} def')
        CALL GKPPB ('/ma {/a exch def} def')
        CALL GKPPB (
     *    '/ct {dup stringwidth pop 2 div neg 0 rmoveto show} def')
        CALL GKPPB ('/rj {dup stringwidth pop neg 0 rmoveto show} def')
        CALL GKPPB ('/lj {show} def')
        CALL GKPPB ('/xy {/y exch def /x exch def} def')
        CALL GKPPB ('/rxy {/ry exch def /rx exch def} def')
        CALL GKPPB ('/sxy {gsave xy x y translate a rotate'//
     *    ' x neg y neg translate} def')
        CALL GKPPB ('/dt {xy np fg x y s 0 360 arc fi} def')
        CALL GKPPB ('/pl {sxy np x y m fg -24 0 srl 48 0 srl'//
     *    ' -24 0 srl 0 24 srl 0 -48 srl sk gr} def')
        CALL GKPPB ('/as {np x y m 0 24 srm 14 -43.4 srl'//
     *    ' -36.8 26.8 srl 45.6 0 srl -36.8 -26.8 srl')
        CALL GKPPB ('14 43.4 srl 14 -43.4 srl} def')
        CALL GKPPB ('/fas {sxy fg as fill fg as csk gr} def')
        CALL GKPPB ('/dc {sxy np x y m fg -24 24 srl 48 -48 srl'//
     *    ' -24 24 srl -24 -24 srl 48 48 srl')
        CALL GKPPB ('sk gr} def')
        CALL GKPPB ('/sq {np x y m 0 24 srm 24 0 srl 0 -48 srl'//
     *    ' -48 0 srl 0 48 srl 24 0 srl} def')
        CALL GKPPB ('/nsq {sxy bg sq fi fg sq csk gr} def')
        CALL GKPPB ('/fsq {sxy fg sq fi fg sq csk gr} def')
        CALL GKPPB ('/ci {np x y 24 s mul 0 360 arc} def')
        CALL GKPPB ('/nci {xy bg ci fi fg ci sk} def')
        CALL GKPPB ('/fci {xy fg ci fi fg ci sk} def')
        CALL GKPPB ('/tu {np x y m 0 28 srm -24 -42 srl'//
     *    ' 48 0 srl -24 42 srl} def')
        CALL GKPPB ('/ntu {sxy bg tu fi fg tu csk gr} def')
        CALL GKPPB ('/ftu {sxy fg tu fi fg tu csk gr} def')
        CALL GKPPB ('/td {np x y m 0 -28 srm -24 42 srl'//
     *    ' 48 0 srl -24 -42 srl} def')
        CALL GKPPB ('/ntd {sxy bg td fi fg td csk gr} def')
        CALL GKPPB ('/ftd {sxy fg td fi fg td csk gr} def')
        CALL GKPPB ('/dm {np x y m 0 24 srm -24 -24 srl'//
     *    ' 24 -24 srl 24 24 srl -24 24 srl} def')
        CALL GKPPB ('/ndm {sxy bg dm fi fg dm csk gr} def')
        CALL GKPPB ('/fdm {sxy fg dm fi fg dm csk gr} def')
        CALL GKPPB ('/bt {np x y m -30 24 srl 0 -48 srl'//
     *    ' 60 48 srl 0 -48 srl -30 24 srl} def')
        CALL GKPPB ('/nbt {sxy bg bt fi fg bt csk gr} def')
        CALL GKPPB ('/fbt {sxy fg bt fi fg bt csk gr} def')
        CALL GKPPB ('/hg {np x y m -24 30 srl 48 0 srl'//
     *    ' -48 -60 srl 48 0 srl -24 30 srl} def')
        CALL GKPPB ('/nhg {sxy bg hg fi fg hg csk gr} def')
        CALL GKPPB ('/fhg {sxy fg hg fi fg hg csk gr} def')
        CALL GKPPB ('/st {sxy bg as fi fg as csk gr} def')
        CALL GKPPB ('/fst {fas} def')
        CALL GKPPB ('/tud {sxy bg tu fi bg td fi fg tu csk'//
     *    ' fg td csk gr} def')
        CALL GKPPB ('/tl {np x y m -14 0 srm 42 -24 srl'//
     *    ' 0 48 srl -42 -24 srl} def')
        CALL GKPPB ('/ftl {sxy fg tl fi fg tl csk gr} def')
        CALL GKPPB ('/tr {np x y m 28 0 srm -42 -24 srl'//
     *    ' 0 48 srl 42 -24 srl} def')
        CALL GKPPB ('/ftr {sxy fg tr fi fg tr csk gr} def')
        CALL GKPPB ('/hpl {np x y m 0 24 srm 8 0 srl'//
     *    ' 0 -16 srl 16 0 srl 0 -16 srl -16 0 srl')
        CALL GKPPB ('0 -16 srl -16 0 srl 0 16 srl -16 0 srl'//
     *    ' 0 16 srl 16 0 srl 0 16 srl 8 0 srl} def')
        CALL GKPPB ('/npl {sxy bg hpl fi fg hpl csk gr} def')
        CALL GKPPB ('/om {np x y m 0 24 srm 16 0 srl'//
     *    ' 8 -8 srl 0 -32 srl -8 -8 srl -32 0 srl')
        CALL GKPPB ('-8 8 srl 0 32 srl 8 8 srl 16 0 srl} def')
        CALL GKPPB ('/nom {sxy bg om fi fg om csk gr} def')
        CALL GKPPB ('/pat1 {/px exch def /pa 16 array def 0 1 15'//
     *    ' {/py exch def /pw 2 string def')
        CALL GKPPB ('pw 0 px py 2 mul 2 getinterval putinterval'//
     *    ' pa py pw put} for} def')
        CALL GKPPB ('/pat2 {/pi exch def /cflag exch def save'//
     *    ' cflag 1 eq {eoclip} {clip}')
        CALL GKPPB ('ifelse newpath {clippath pathbbox} stopped'//
     *    ' not {/ph exch def /pw exch def')
        CALL GKPPB ('/py exch def /px exch def /px px 256 div'//
     *    ' floor 256 mul def')
        CALL GKPPB ('/py py 256 div floor 256 mul def px py translate'//
     *    ' /pw pw px sub 256 div')
        CALL GKPPB ('floor 1 add cvi def /ph ph py sub 256 div'//
     *    ' floor 1 add cvi def')
        CALL GKPPB ('pw 256 mul ph 256 mul scale /pw pw 32 mul def'//
     *    ' /ph ph 32 mul def')
        CALL GKPPB ('/px 0 def /py 0 def pw ph pi'//
     *    ' [pw neg 0 0 ph neg pw ph] {pa py get')
        CALL GKPPB ('/px px 16 add def px pw ge {/px 0 def'//
     *    ' /py py 1 add 16 mod def} if} pi type')
        CALL GKPPB ('/booleantype eq {imagemask} {image} ifelse}'//
     *    ' if restore} def')
        CALL GKPPB ('/bp {closepath gsave} def')
        CALL GKPPB ('/ep {pat1 1 1 pat2 grestore} def')

        CALL GKPPB ('/OF /findfont load def')
        CALL GKPPB ('/findfont {dup GLI_GKS_dict exch known')
        CALL GKPPB ('{GLI_GKS_dict exch get}')
        CALL GKPPB ('if GLI_GKS_dict /OF get exec} def')
        CALL GKPPB ('mark')
        CALL GKPPB ('/ISOLatin1Encoding 8#000 1 8#001'//
     *      ' {StandardEncoding exch get} for')
        CALL GKPPB ('/emdash /endash 8#004 1 8#025'//
     *      ' {StandardEncoding exch get} for')
        CALL GKPPB ('/quotedblleft /quotedblright 8#030 1 8#054'//
     *      ' {StandardEncoding exch get} for')
        CALL GKPPB ('/minus 8#056 1 8#217'//
     *      ' {StandardEncoding exch get} for')
        CALL GKPPB ('/dotlessi 8#301 1 8#317'//
     *      ' {StandardEncoding exch get} for')
        CALL GKPPB ('/space/exclamdown/cent/sterling/currency/yen'//
     *      '/brokenbar/section')
        CALL GKPPB ('/dieresis/copyright/ordfeminine/guillemotleft'//
     *      '/logicalnot/hyphen/registered')
        CALL GKPPB ('/macron/degree/plusminus/twosuperior'//
     *      '/threesuperior/acute/mu/paragraph')
        CALL GKPPB ('/periodcentered/cedilla/onesuperior/ordmasculine'//
     *      '/guillemotright/onequarter')
        CALL GKPPB ('/onehalf/threequarters/questiondown/Agrave'//
     *      '/Aacute/Acircumflex/Atilde')
        CALL GKPPB ('/Adieresis/Aring/AE/Ccedilla/Egrave/Eacute'//
     *      '/Ecircumflex/Edieresis/Igrave')
        CALL GKPPB ('/Iacute/Icircumflex/Idieresis/Eth/Ntilde/Ograve'//
     *      '/Oacute/Ocircumflex/Otilde')
        CALL GKPPB ('/Odieresis/multiply/Oslash/Ugrave/Uacute'//
     *      '/Ucircumflex/Udieresis/Yacute/Thorn')
        CALL GKPPB ('/germandbls/agrave/aacute/acircumflex/atilde'//
     *      '/adieresis/aring/ae/ccedilla')
        CALL GKPPB ('/egrave/eacute/ecircumflex/edieresis/igrave'//
     *      '/iacute/icircumflex/idieresis')
        CALL GKPPB ('/eth/ntilde/ograve/oacute/ocircumflex/otilde'//
     *      '/odieresis/divide/oslash/ugrave')
        CALL GKPPB ('/uacute/ucircumflex/udieresis/yacute/thorn'//
     *      '/ydieresis')
        CALL GKPPB ('256 array astore def cleartomark')
        CALL GKPPB ('/encodefont {findfont dup maxlength dict begin')
        CALL GKPPB ('{1 index /FID ne {def} {pop pop} ifelse} forall')
        CALL GKPPB ('/Encoding exch def dup')
        CALL GKPPB ('/FontName exch def currentdict')
        CALL GKPPB ('definefont end} def')
        CALL GKPPB ('end')
        
        CALL GKPPB ('%%EndProcSet')
        CALL GKPPB ('%%EndProlog')

        END IF
        
        PAGES = PAGES + 1
        CALL GDEC (PAGES,LPAGE,PAGE)
        CALL GKPPB ('%%Page: '//PAGE(1:LPAGE)//' '//PAGE(1:LPAGE))
        
        CALL GKPPB ('%%BeginPageSetup')
        CALL GKPPB ('GLI_GKS_dict begin save /psl exch def')
        
        IF (WTYPE .GE. 63) THEN
          CALL GKPPB ('initgraphics')
          CALL GKPPB ('1 setgray clippath fill')
        END IF

        IF (WTYPE .LT. 63) THEN
            IF (MOD(WTYPE,2) .EQ. 0) THEN
                X = '21'
            ELSE
                X = '18'
            END IF
            IF (LANDSC) THEN
                CALL GDEC (YTRANS,LY,Y)
                CALL GKPPB (X//' '//Y(1:LY)//' translate -90 rotate')
            ELSE
                CALL GKPPB (X//' 15 translate')
            END IF
        ELSE
            IF (LANDSC) THEN
                CALL GDEC (YTRANS,LY,Y)
                CALL GKPPB ('0 '//Y(1:LY)//' translate -90 rotate')
            END IF
        ENDIF

        IF (ABS(MAGS) .GT. FEPS) THEN
            CALL GFLT (1.2**MAGS,LMAG,MAG)
            CALL GKPPB (MAG(1:LMAG)//' 1 in 600 div mul dup scale')
        ELSE
            CALL GKPPB ('1 in 600 div dup scale')
        END IF

        CALL GKPSCI (-1,WTYPE)
        CALL GKPSFG (-1,WTYPE)
        CALL GKPPB ('0 setlinecap 1 setlinejoin')
        CALL GKPSLW (-1.0)
        CALL GKPSMS (-1.0)
        CALL GKPPB ('0 ma')
        CALL GKPSCH (-1,HEIGHT)
C*  set clipping rectangle
        CALL GKPCLP (WN)
        CALL GKPPB ('%%EndPageSetup')
        CALL GKPU

        RETURN


        ENTRY GKPEP (PAGES)
C*  end page

        CALL GDEC (PAGES,LPAGE,PAGE)
        CALL GKPPB ('%%EndPage: '//PAGE(1:LPAGE)//' '//PAGE(1:LPAGE))
        CALL GKPU

        END


        SUBROUTINE GKPICT
C*  initialize color table

        INTEGER ICOLOR, RINT, GINT, BINT, WTYPE
        REAL RED, GREEN, BLUE
        CHARACTER*(*) SCI
        INTEGER LSCI

        INTEGER COLOR, FCOL, PIXEL
        REAL R(0:979), G(0:979), B(0:979)

        REAL GREY
        INTEGER I, J, CI
        CHARACTER*12 SR, SG, SB, PIX
        INTEGER LSR, LSG, LSB, LPIX

        SAVE

C*  GKS symbol definitions
        INCLUDE 'gksdefs.i'
C*  GKS description table
        INCLUDE 'gksdescr.i'
C*  GKS state list
        INCLUDE 'gksstate.i'

        DATA COLOR /1/, FCOL /1/

        DO 1 I = 0,979
           J = I
           CALL GQRGB (J,R(J),G(J),B(J))
   1    CONTINUE

        COLOR = -1

        RETURN


        ENTRY GKPSCR (ICOLOR, RED, GREEN, BLUE)
C*  set color representation

        IF (ICOLOR .GE. 0 .AND. ICOLOR .LT. 980) THEN
           R(ICOLOR) = RED
           G(ICOLOR) = GREEN
           B(ICOLOR) = BLUE
        END IF

        RETURN


        ENTRY GKPSCI (ICOLOR, WTYPE)
C*  set color index

        IF (ICOLOR .LT. 980) THEN
           IF (ICOLOR .NE. COLOR) THEN
              CI = ABS(ICOLOR)
              IF (WTYPE .GE. 63 .AND. X11) THEN
                 CALL GQPIX (CI,PIXEL)
                 CALL GDEC (PIXEL,LPIX,PIX)
                 CALL GKPPB (PIX(1:LPIX)//' sp')
              ELSE
                 IF (MOD(WTYPE,2) .NE. 0) THEN
                    GREY = 0.3*R(CI)+0.59*G(CI)+0.11*B(CI)
                    CALL GFLT (GREY,LSG,SG)
                    CALL GKPPB (SG(1:LSG)//' sg')
                 ELSE
                    CALL GFLT (R(CI),LSR,SR)
                    CALL GFLT (G(CI),LSG,SG)
                    CALL GFLT (B(CI),LSB,SB)
                    CALL GKPPB (SR(1:LSR)//' '//SG(1:LSG)//
     *                 ' '//SB(1:LSB)//' sc')
                 END IF
              END IF
              COLOR = CI
           END IF
        END IF

        RETURN


        ENTRY GKPSFG (ICOLOR, WTYPE)
C*  set foreground color

        IF (ICOLOR .LT. 980) THEN
           IF (ICOLOR .NE. FCOL) THEN
              CI = ABS(ICOLOR)
              IF (WTYPE .GE. 63 .AND. X11) THEN
                 CALL GQPIX (CI,PIXEL)
                 CALL GDEC (PIXEL,LPIX,PIX)
                 CALL GKPPB ('/fg {'//PIX(1:LPIX)//' sp} def')
              ELSE
                 IF (MOD(WTYPE,2) .NE. 0) THEN
                    GREY = 0.3*R(CI)+0.59*G(CI)+0.11*B(CI)
                    CALL GFLT (GREY,LSG,SG)
                    CALL GKPPB ('/fg {'//SG(1:LSG)//' sg} def')
                 ELSE
                    CALL GFLT (R(CI),LSR,SR)
                    CALL GFLT (G(CI),LSG,SG)
                    CALL GFLT (B(CI),LSB,SB)
                    CALL GKPPB ('/fg {'//SR(1:LSR)//' '//SG(1:LSG)//
     *                  ' '//SB(1:LSB)//' sc} def')
                 END IF
              END IF
              FCOL = CI
           END IF
           IF (ICOLOR .NE. COLOR) THEN
              CI = ABS(ICOLOR)
              CALL GKPPB ('fg')
              COLOR = CI
           END IF
        END IF

        RETURN


        ENTRY GKPSBG (WTYPE)
C*  set background color

        IF (WTYPE .GE. 63 .AND. X11) THEN
           CALL GQPIX (0,PIXEL)
           CALL GDEC (PIXEL,LPIX,PIX)
           CALL GKPPB ('/bg {'//PIX(1:LPIX)//' sp} def')
        ELSE
           IF (MOD(WTYPE,2) .NE. 0) THEN
              GREY = 0.3*R(0)+0.59*G(0)+0.11*B(0)
              CALL GFLT (GREY,LSG,SG)
              CALL GKPPB ('/bg {'//SG(1:LSG)//' sg} def')
           ELSE
              CALL GFLT (R(0),LSR,SR)
              CALL GFLT (G(0),LSG,SG)
              CALL GFLT (B(0),LSB,SB)
              CALL GKPPB ('/bg {'//SR(1:LSR)//' '//SG(1:LSG)//
     *           ' '//SB(1:LSB)//' sc} def')
           END IF
        END IF

        RETURN


        ENTRY GKPQCI (ICOLOR, LSCI, SCI, WTYPE)
C*  inquire color index

        IF (MOD(WTYPE,2) .EQ. 0) THEN
            RINT = NINT(R(ICOLOR)*255)
            GINT = NINT(G(ICOLOR)*255)
            BINT = NINT(B(ICOLOR)*255)
            CALL GKPHEX (RINT, SCI(1:2))
            CALL GKPHEX (GINT, SCI(3:4))
            CALL GKPHEX (BINT, SCI(5:6))
            LSCI = 6
        ELSE
            GREY = 0.3*R(ICOLOR)+0.59*G(ICOLOR)+0.11*B(ICOLOR)
            GINT = NINT(GREY*255)
            CALL GKPHEX (GINT, SCI(1:2))
            LSCI = 2
        END IF

        RETURN
        END


        SUBROUTINE GKPHEX (I, STR)
C*  convert integer to hex notation

        INTEGER I
        CHARACTER*(*) STR

        CHARACTER*1 HTAB(0:15)
        DATA HTAB /'0','1','2','3','4','5','6','7','8','9',
     *    'A','B','C','D','E','F'/

        STR(1:1) = HTAB(I/16)
        STR(2:2) = HTAB(MOD(I,16))

        RETURN
        END


        SUBROUTINE GKPSCO (ID,IT)
C*  set up connection

        INTEGER I, L, ID, IPNTR, LPNTR
        INTEGER CONID, WTYPE

        CHARACTER*(*) BUFF
        CHARACTER*500 IOBUFF
        CHARACTER*1 LF

        SAVE

        DATA IPNTR /0/, LPNTR /0/

        LF = CHAR(10)

        CONID = ID
        WTYPE = IT

        RETURN


        ENTRY GKPPB (BUFF)
C*  pack buffer

        L = LEN(BUFF)

        IF (BUFF(1:1) .EQ. '%') THEN
          IF (LPNTR .NE. 0) THEN
            IPNTR = IPNTR+1
            IOBUFF(IPNTR:IPNTR) = LF
            LPNTR = 0
          END IF

        ELSE IF (L .GT. 78-LPNTR) THEN
          IF (IPNTR .NE. 0) THEN
            IPNTR = IPNTR+1
            IOBUFF(IPNTR:IPNTR) = LF
            LPNTR = 0
          END IF
        END IF

        IF (L+2 .GT. 500-IPNTR) THEN
C*  xfer buffer
          IF (WTYPE .GE. 63) THEN
            CALL DPSWR (IPNTR,IOBUFF)
          ELSE
            CALL BUFOUT (CONID,IPNTR,IOBUFF)
          END IF
          IPNTR = 0
        END IF

        IF (LPNTR .NE. 0) THEN
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = ' '
          LPNTR = LPNTR+1
        END IF

        DO 1 I = 1,L
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = BUFF(I:I)
          LPNTR = LPNTR+1
   1    CONTINUE

        IF (BUFF(1:1) .EQ. '%') THEN
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = LF
          LPNTR = 0
        END IF

        RETURN


        ENTRY GKPU
C*  update

        IF (LPNTR .NE. 0) THEN
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = LF
          LPNTR = 0
        END IF
        IF (WTYPE .GE. 63) THEN
          CALL DPSWR (IPNTR,IOBUFF)
          CALL DPSWR (1,LF)
          CALL DPSFL
        ELSE
          CALL BUFOUT (CONID,IPNTR,IOBUFF)
        END IF
        IPNTR = 0

        RETURN
        END
