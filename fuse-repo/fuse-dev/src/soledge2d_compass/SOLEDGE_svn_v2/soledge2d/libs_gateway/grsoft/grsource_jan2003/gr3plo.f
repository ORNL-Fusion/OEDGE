C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C Fuer den Microsoft-Compiler muss das Quellprogramm in wenigstens zwei
C Teile zerlegt werden; die Grenzlinie muss dann vor der Zeile
C C========= ... sein.
C Dieser Compiler ueberprueft sonst die Uebereinstimmung von aktuellen
C Argumenten mit formalen Argumenten zu genau.
C***********************************************************************
C     DEBUG SUBCHK
C     END DEBUG
c update   6. 11. 90 m. busch    save
c update  10.  1. 91 m. busch    common grpic eingefuegt
c                                xcm=39.5, ycm=28.7 durch xmaxcm,ymaxcm
c                                ersetzt
c update  14. 03. 91 G.GRoten    V002,DIAG2
C UPDATE  17. 09. 91 G.Groten    GRSCLC,GRSCLV zurueckgesetzt
C UPDATE  24. 10. 91 G.Groten    GR3SOR, GR3CBC und GR3HDD: IFPT-Index
C UPDATE  29. 10. 91 G.Groten    LIND(43,2) und IF ( ... .EQ.0 )goto 20
C UPDATE  07. 11. 91 G.GROTEN    ISMIN,ISMAX -> GRISMI,GRISMA WEGEN ESSL
C UPDATE  11. 11. 91 G.GROTEN    Kleine Aenderung ( numerisch IBM=CRAY)
C UPDATE  02. 12. 98 G.GROTEN    GR3OUT entfernt
c-----------------------------------------------------------------------
      SUBROUTINE GR3PLO(A,IER,HOW)

c     .. parameters ..
      INTEGER           MFL
      PARAMETER         (MFL=256)
c     ..
c     .. scalar arguments ..
      INTEGER           IER
      CHARACTER*(*)     HOW
c     ..
c     .. array arguments ..
      INTEGER           FLPIC,IF0,IF1,IF2,IF3,IF4,IFL,IL,IL0,IP0,ISY,
     $                  ISZ,IW,KINDA,KOAX,LQ,LW,MAFL,MNF,MNL,MNP,NOTDIM,
     $                  NSCLC,NSCLP,NSCLV,RAHMEN
      REAL              A(46*MNP),EX4(4)
c     ..
c     .. local scalars ..
      REAL              XUR,YUR
      INTEGER           I,KIND,KINDAX
c     ..
c     .. local arrays ..
      REAL              EXT4(4)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3DZN,GR3PFL,GR3PLT,GRNWPN,GRSPTS
c     ..
c     .. intrinsic functions ..
      INTRINSIC         INDEX
c     ..
c     .. common blocks ..
      COMMON            /GR3DRE/IFDREH,XDREH,YDREH,ZDREH,KINDA,KOAX
CDEC$ PSECT /GR3DRE/ NOSHR
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      COMMON            /GR4COM/IFL,JGR,XC,YC,COMA
CDEC$ PSECT /GR4COM/ NOSHR
      COMMON            /GRPIC/FLPIC,NSCLC,NSCLV,NSCLP,RAHMEN,XMXCMQ,
     $                  YMXCMQ,XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON            /GRPP/PP
CDEC$ PSECT /GRPP/ NOSHR
      REAL              XDCPIC,XDREH,XMXCMQ,YDCPIC,YDREH,YMXCMQ,ZDREH
      LOGICAL           IFDREH,WIN,ZENPRO
      REAL              COMA(3,2),GR(3,3),PP(18),RT(3,3),XC(500),YC(500)
      INTEGER           JGR(7,MFL)
c     ..
c     .. save statement ..
      SAVE /GR3DRE/,/GR3SAV/,/GR4COM/,/GRPIC/,/GRPP/
      SAVE
c     ..
c     .. executable statements ..

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3PLO RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 20
*
      END IF

      IF (HOW(1:3).EQ.'ALL') THEN
         KIND   = 1
*
      ELSE IF (HOW(1:3).EQ.'MIX') THEN
         KIND   = 2
*
      ELSE IF (HOW(1:3).EQ.'HID') THEN
         KIND   = 3
*
      ELSE
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3PLO RC=20: MUST BE "ALL" OR "MIX" OR "HID"'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 20
*
      END IF
*
      IF (IFL.EQ.0) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3PLO RC=30: NO SURFACE IS BUILT!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 30
         GO TO 20
*
      END IF
*
      IF (.NOT.WIN) THEN
         EXT4(1) = 0.
         EXT4(3) = 0.
      END IF
*
      IF (INDEX(HOW,',D').GT.0) THEN
         XUR    = 3.
         YUR    = 3.
         IF (XMXCMQ-PP(3)-.000005.GT.PP(1)) XUR    = XMXCMQ - 3.
         IF (YMXCMQ-PP(4)-.000005.GT.PP(2)) YUR    = YMXCMQ - 3.
         IF (RT(3,1).GT.0.) THEN
            CALL GRNWPN(1)
            CALL GRSPTS(20)
*
         ELSE
            CALL GRNWPN(3)
            CALL GRSPTS(18)
         END IF
*
         CALL GR3PFL(XUR,YUR,RT(1,1),RT(2,1),'X')
         IF (RT(3,2).GT.0.) THEN
            CALL GRNWPN(1)
            CALL GRSPTS(20)
*
         ELSE
            CALL GRNWPN(3)
            CALL GRSPTS(18)
         END IF
*
         CALL GR3PFL(XUR,YUR,RT(1,2),RT(2,2),'Y')
         IF (RT(3,3).GT.0.) THEN
            CALL GRNWPN(1)
            CALL GRSPTS(20)
*
         ELSE
            CALL GRNWPN(3)
            CALL GRSPTS(18)
         END IF
*
         CALL GR3PFL(XUR,YUR,RT(1,3),RT(2,3),'Z')
      END IF
*
      IF (INDEX(HOW,',A1').GT.0) THEN
         KINDAX = 1
*
      ELSE IF (INDEX(HOW,',A2').GT.0) THEN
         KINDAX = 2
*
      ELSE IF (INDEX(HOW,',A3').GT.0) THEN
         KINDAX = 3
*
      ELSE IF (INDEX(HOW,',A4').GT.0) THEN
         KINDAX = 4
*
      ELSE IF (INDEX(HOW,',A5').GT.0) THEN
         KINDAX = 5
*
      ELSE IF (INDEX(HOW,',A6').GT.0) THEN
         KINDAX = 6
*
      ELSE
         KINDAX = KINDA
      END IF

      CALL GR3PLT(KIND,A,A(ISY),A(ISZ),IP0,A(IF1),A(IF2),A(IF3),
     $            A(IF4),MNF,IF0,A(IL),A(LQ),MNL,IL0,A(IW),LW,EXT4,
     $            KINDAX,IER)
*
*
      IF (ZENPRO) THEN
c----    call gr3dzn(  x  ,  y   ,   z  ,ip0)
         CALL GR3DZN(A,A(ISY),A(ISZ),IP0)
         ZENPRO = .FALSE.
      END IF
*
      IF (IER.LT.10) IER    = 0

      GO TO 20
c-----------------------------------------------------------------------
      ENTRY             GR3WIN(EX4,IER)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3WIN RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 20
*
      END IF

      IF (EX4(1).GT.EX4(3) .OR. EX4(2).GT.EX4(4)) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3WIN RC=20: XMIN>XMAX    OR   YMIN>YMAX'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 20
*
      END IF

      IF (EX4(1).EQ.EX4(3) .AND. EX4(2).EQ.EX4(4)) THEN
         WIN    = .FALSE.
*
      ELSE
         WIN    = .TRUE.
      END IF
*
      DO 10 I = 1,4
         EXT4(I) = EX4(I)
   10 CONTINUE
      IER    = 0
   20 CONTINUE
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3CEN(A,IER,AB)

c     .. scalar arguments ..
      REAL              AB
      INTEGER           IER
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              A(46*MNP)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3DZN,GR3ZEN
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN
      LOGICAL           ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3CEN RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 10
*
      END IF

      IF (AB.NE.0.) THEN
         CALL GR3ZEN(AB,A,A(ISY),A(ISZ),IP0,ZENPRO)
         ZENPRO = .TRUE.
*
      ELSE
         CALL GR3DZN(A,A(ISY),A(ISZ),IP0)
         ZENPRO = .FALSE.
      END IF
*
      IER    = 0
   10 CONTINUE
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3ZEN(Z1,X,Y,Z,IP0,ZENPRO)
c---- es wird auf die ebene z=0 projeziert.
c---------------------------------------------------------------------

c     .. scalar arguments ..
      REAL              Z1
      LOGICAL           ZENPRO
      INTEGER           IP0
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,ISY,ISZ,IW,LQ,
     $                  LW,MAFL,MNF,MNL,MNP,NOTDIM
      REAL              X(MNP),Y(MNP),Z(MNP)
c     ..
c     .. local scalars ..
      REAL              FAZ,Z0
      INTEGER           I
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MAX,MIN
c     ..
c     .. common blocks ..
      COMMON            /GR3BAC/ICOLBA,CENTR,WIND
CDEC$ PSECT /GR3BAC/ NOSHR
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IPQ,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRQ,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              CENTR,WIN
      REAL              GR(3,3),RT(3,3),WIND(4)
      INTEGER           ICOLBA(8)
c     ..
c     .. save statement ..
      SAVE /GR3BAC/,/GR3SAV/
      SAVE
c     ..
c     .. data statements ..
      DATA              Z0/0./
c     ..
c     .. executable statements ..

c---- normalfall: noch keine zentralprojektion gemacht
      IF (.NOT.ZENPRO) THEN
         WIND(1) = Z1/ (Z1-Z(1))*X(1)
         WIND(2) = Z1/ (Z1-Z(1))*Y(1)
         WIND(3) = WIND(1)
         WIND(4) = WIND(2)
         DO 10 I = 1,IP0
            FAZ    = Z1/ (Z1-Z(I))
            X(I)   = FAZ*X(I)
            Y(I)   = FAZ*Y(I)
            WIND(1) = MIN(WIND(1),X(I))
            WIND(2) = MIN(WIND(2),Y(I))
            WIND(3) = MAX(WIND(3),X(I))
            WIND(4) = MAX(WIND(4),Y(I))
   10    CONTINUE
         Z0     = Z1
c---- sonderfall: schon eine zentralprojektion gemacht
      ELSE
         WIND(1) = Z1* (Z0-Z(1))/ (Z0* (Z1-Z(1)))*X(1)
         WIND(2) = Z1* (Z0-Z(1))/ (Z0* (Z1-Z(1)))*Y(1)
         WIND(3) = WIND(1)
         WIND(4) = WIND(2)
         DO 20 I = 1,IP0
            FAZ    = Z1* (Z0-Z(I))/ (Z0* (Z1-Z(I)))
            X(I)   = FAZ*X(I)
            Y(I)   = FAZ*Y(I)
            WIND(1) = MIN(WIND(1),X(I))
            WIND(2) = MIN(WIND(2),Y(I))
            WIND(3) = MAX(WIND(3),X(I))
            WIND(4) = MAX(WIND(4),Y(I))
   20    CONTINUE
         Z0     = Z1
      END IF

      GO TO 40
c----------------------------------------------------------------------
      ENTRY             GR3DZN(X,Y,Z,IP0)
c---- die zentralprojektion wird rueckgaengig gemacht
      IF (Z0.NE.0.) THEN
         DO 30 I = 1,IP0
            FAZ    = 1. - Z(I)/Z0
            X(I)   = FAZ*X(I)
            Y(I)   = FAZ*Y(I)
   30    CONTINUE
         Z0     = 0.
      END IF

   40 CONTINUE
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3SOR(NR,FMM,INDF,NEXT,INDEXF,IFPT,IDONE,FLAG,MNF,EX4)
c    sorts the faces according to buckets met.  first the faces are
c    sorted using a linked list and then they are rearranged so that
c    the indexes of the faces meeting the rectangle with coordinates
c    ix,iy lie in position ifpt(ix,iy) and ifpt(ix+1,iy) in indf.
c    beware of the fortran trick discribed above.




c     .. scalar arguments ..
      INTEGER           IDONE,INDF,MNF,NR
      LOGICAL           FLAG
c     ..
c     .. array arguments ..
      REAL              EX4(4),FMM(5,INDF)
      INTEGER           IFPT(NR*NR+1),INDEXF(4*MNF),NEXT(4*MNF)
c     ..
c     .. local scalars ..
      REAL              XDIV,YDIV
      INTEGER           I,II,IND,IND1,IND2,INEXT,INEXT1,INTEMP,IPOINT,
     $                  IX,IXHIGH,IXLOW,IY,IYHIGH,IYLOW
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MIN
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      DO 10 I = 1,MNF*4
         NEXT(I) = 0
         INDEXF(I) = 0
   10 CONTINUE

      DO 30 IY = 1,NR
         DO 20 IX = 1,NR
            IFPT(NR*(IY-1)+IX) = 0
   20    CONTINUE
   30 CONTINUE

      XDIV   = NR/ (EX4(3)-EX4(1))
      YDIV   = NR/ (EX4(4)-EX4(2))
      IPOINT = 0

      DO 60 II = 1,INDF
         IF (IPOINT.GT. (4*MNF-NR*NR)) THEN
            FLAG   = .TRUE.
            IDONE  = II - 1
            GO TO 100
*
         END IF

         IXLOW  = (FMM(1,II)-EX4(1))*XDIV + 1.
         IXHIGH = (FMM(2,II)-EX4(1))*XDIV + 1.
         IXHIGH = MIN(NR,IXHIGH)
         IYLOW  = (FMM(3,II)-EX4(2))*YDIV + 1.
         IYHIGH = (FMM(4,II)-EX4(2))*YDIV + 1.
         IYHIGH = MIN(NR,IYHIGH)
         DO 50 IX = IXLOW,IXHIGH
            DO 40 IY = IYLOW,IYHIGH
               IPOINT = IPOINT + 1
               IND    = IFPT(NR*(IY-1)+IX)
               NEXT(IPOINT) = IND
               INDEXF(IPOINT) = II
               IFPT(NR*(IY-1)+IX) = IPOINT
   40       CONTINUE
   50    CONTINUE

   60 CONTINUE

      IPOINT = 0

      DO 90 IY = 1,NR
         DO 80 IX = 1,NR
            INEXT  = IFPT(NR*(IY-1)+IX)
            IFPT(NR*(IY-1)+IX) = IPOINT + 1
c---------- do while (inext.gt.0) --------------------------------------
   70       CONTINUE
            IF (INEXT.GT.0) THEN
               IF (INEXT.LE.IPOINT) THEN
                  INEXT  = NEXT(INEXT)
*
               ELSE
                  IPOINT = IPOINT + 1
                  IF (IPOINT.NE.INEXT) THEN
                     IND1   = INDEXF(INEXT)
                     IND2   = INDEXF(IPOINT)
                     INEXT1 = NEXT(IPOINT)
                     INTEMP = NEXT(INEXT)
                     NEXT(IPOINT) = INEXT
                     NEXT(INEXT) = INEXT1
                     INDEXF(INEXT) = IND2
                     INDEXF(IPOINT) = IND1
                     INEXT  = INTEMP
*
                  ELSE
                     INEXT  = NEXT(INEXT)
                  END IF
*
               END IF
*
               GO TO 70
*
            END IF
*
   80    CONTINUE
   90 CONTINUE

      IFPT(NR*NR+1) = IPOINT + 1

  100 CONTINUE
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3HDD(NR,X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,ILINE,INDF,
     $                  INDL1,LNFC,MNP,MNL,MNF,EX4,FMM,ENRM,IUP,IDOWN,
     $                  DIFF,IFPT,INDEXF,IFLIST,IFREE)

c  x,y,z is an array-group of three dimensional points.
c  iface 1 2 3 4 is an array of integers
c  which indicate the points in x,y,z which belong to the indf
c  faces under consideration.  iline is an array which indicates the
c  endpoints of line segments.  these points are also stored in x,y,z.
c  gr3hdd compares faces to lines and sends visible segments to plotter.
c  if a line is compared to a face on which it lies, it may dissappear
c  due to roundoff error.

c       if lnfc(1,i) = j and lnfc(2,i) = k, the i-th line will not be
c  compared to the j-th or k-th face and the points making up
c  line i are assumed to be among points numbered 1,2, and 4 in face
c  j and 2,3, and 4 in face k.  the last is irrelevant unless the
c  normals of the trianges formed by 1,2,4 and 2,3,4 have z-coordinates
c  of different sign.
c   mnf, and mnl are the maximum numbers of
c   faces, and lines respectively.

******************************************



******************************************






c     next : for backside colors
c     commonblock der standardwerte bzw. der geaenderten tabellenwerte
c     .. parameters ..
      INTEGER           MFL,LXC
      REAL              TEIL
      PARAMETER         (MFL=256,LXC=500,TEIL=.00390625)
c     ..
c     .. scalar arguments ..
      INTEGER           INDF,INDL1,MNF,MNL,MNP,NR
c     ..
c     .. array arguments ..
      REAL              ENRM(2,4,MNF),EX4(4),FMM(5,MNF),X(MNP),Y(MNP),
     $                  Z(MNP)
      INTEGER           IFACE1(MNF),IFACE2(MNF),IFACE3(MNF),IFACE4(MNF),
     $                  IFLIST(MNF),IFPT(NR*NR+1),ILINE(2,MNL),
     $                  INDEXF(4*MNF),LNFC(2,MNL)
      LOGICAL           DIFF(0:MNF),IDOWN(2,MNL),IFREE(MNF),IUP(2,MNL)
c     ..
c     .. local scalars ..
      REAL              D,DET,DIAG2,EZLINE,EZQUAD,FLGROT,LX1,LX2,LY1,
     $                  LY2,LZ1,LZ2,QX1,QX2,QX3,QX4,QY1,QY2,QY3,QY4,R1,
     $                  R2,R3,R4,RMAX,RMIN,T,T1,T2,T3,T4,T5,TDET,TEMP1,
     $                  TEMP2,V1,V2,V3,V4,VV1,VV2,VV3,VV4,VV5,VV6,W0,
     $                  W01,W02,X0,X1,X2,XDIR,XLINE,XMAXL,XMINL,Y0,Y1,
     $                  Y2,YDIR,YLINE,YMAXL,YMINL,Z0,Z1,Z2,ZEIGRO,
     $                  ZEIWIN,ZLINE,ZMINL,V002
      INTEGER           I,I11,I12,I21,I22,I31,I32,I41,I42,IC,ICOL,ICOLO,
     $                  IFA,IFAOLD,II,IJ,IND1,IND2,INDL,INDP,INT1,INT2,
     $                  INT3,INT4,INTLIN,INTOL,INTSYM,IOBER,IP1B,IP2B,
     $                  SIZMRK,ISURF,IUNTER,JJ,KK,LNF1,LNF2,NUM,NX,NY,
     $                  PAVAIL,PSTART
      LOGICAL           DREI1B,DREI2B,FLAG1,FLAG2,LDOWN,LGCL1,LGCL2,
     $                  LGCL3,LGCL4,LGCL5,LGCL6,LN1,LN2,LP1,LP2,LUP,SEE
c     ..
c     .. local arrays ..
      REAL              RM(2,100)
      INTEGER           AVAIL(100),LIND(43,2),PRM(100)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3BUC,GR3CBC,GR3FXL,GR3NRM,GRDRW,GRJMP,GRNWPN,
     $                  GRSPTS
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MAX,MIN
c     ..
c     .. common blocks ..
      COMMON            /GR3BAC/ICOLBA,CENTR,WIND
CDEC$ PSECT /GR3BAC/ NOSHR
      COMMON            /GR4COM/IFL,JGR,XC,YC,COMA
CDEC$ PSECT /GR4COM/ NOSHR
      COMMON            /GRPP/PP
CDEC$ PSECT /GRPP/ NOSHR
      REAL              CENTR
      INTEGER           IFL
      REAL              COMA(3,2),PP(18),WIND(4),XC(LXC),YC(LXC)
      INTEGER           ICOLBA(8),JGR(7,MFL)
c     ..
c     .. equivalences ..
      EQUIVALENCE       (PP(14),ZEIGRO), (PP(15),ZEIWIN),
     $                  (PP(16),INTSYM), (PP(13),INTLIN),
     $                  (PP(17),SIZMRK), (PP(18),FLGROT)
c     ..
c     .. save statement ..
      SAVE /GR3BAC/,/GR4COM/,/GRPP/
      SAVE
c     ..
c     .. executable statements ..
      INDL   = INDL1

      CALL GR3NRM(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,ENRM,INDF,DIFF)
      DO 10 I = 1,INDL
         IUP(1,I) = .TRUE.
         IUP(2,I) = .TRUE.
         IDOWN(1,I) = .TRUE.
         IDOWN(2,I) = .TRUE.
   10 CONTINUE

      DO 20 I = 1,INDL
         IND1   = LNFC(1,I)

         IF (IND1.NE.0) THEN
            IUP(1,I) = .FALSE.
            IDOWN(1,I) = .FALSE.
            IF (DIFF(IND1)) IDOWN(1,INDL) = .TRUE.
         END IF

         IND2   = LNFC(2,I)
         IF (IND2.NE.0) THEN
            IUP(2,I) = .FALSE.
            IDOWN(2,I) = .FALSE.
            IF (DIFF(IND2)) IUP(2,INDL) = .TRUE.
         END IF

   20 CONTINUE

      DO 30 I = 1,MNF
         IFREE(I) = .TRUE.
   30 CONTINUE

      KK     = 0

      DO 40 I = 1,100
         AVAIL(I) = I
   40 CONTINUE

      PAVAIL = 100

      NX     = NR
      NY     = NR

      D      = MAX(EX4(3)-EX4(1),EX4(4)-EX4(2))*0.001

      INTOL  = INTLIN
      IOBER  = 0
      DO 80 ISURF = 1,IFL + KK
         IUNTER = IOBER + 1
         IOBER  = JGR(4,ISURF)
         SEE    = JGR(5,ISURF) .GT. 0
         IF (.NOT.SEE) GO TO 80
         CALL GRSPTS(JGR(5,ISURF))
         IF (JGR(5,ISURF).GT.15) THEN
            V002=0.002
         ELSE
            V002=0.000002
         ENDIF
         DIAG2  =  V002**2* ((PP(7)-PP(5))**2+ (PP(8)-PP(6))**2)
         ICOL   = JGR(3,ISURF)
         ICOLO  = ICOLBA(ICOL)
         IFA    = ICOL
         IFAOLD = IFA
         CALL GRNWPN(IFA)
         IC     = 0
         DO 70 II = IUNTER,IOBER
            PSTART = AVAIL(PAVAIL)
            RM(1,PSTART) = 0.0
            RM(2,PSTART) = 1.0
            PRM(PSTART) = 0
            PAVAIL = PAVAIL - 1

            IND1   = ILINE(1,II)
            IND2   = ILINE(2,II)

            LX1    = X(IND1)
            LY1    = Y(IND1)
            LX2    = X(IND2)
            LY2    = Y(IND2)
            CALL GR3BUC(EX4,LX1,LY1,LX2,LY2,NUM,LIND,NX,NY)
            CALL GR3CBC(NUM,NR,LIND,INDEXF,IFPT,IFLIST,IFREE)

            LZ1    = Z(IND1)
            LZ2    = Z(IND2)

            XDIR   = LX2 - LX1
            YDIR   = LY2 - LY1

            IF (LX1.LE.LX2) THEN
               XMAXL  = LX2
               XMINL  = LX1
*
            ELSE
               XMAXL  = LX1
               XMINL  = LX2
            END IF

            IF (LY1.LE.LY2) THEN
               YMAXL  = LY2
               YMINL  = LY1
*
            ELSE
               YMAXL  = LY1
               YMINL  = LY2
            END IF

            ZMINL  = MIN(LZ1,LZ2)

            DO 50 IJ = 1,NUM
               JJ     = IFLIST(IJ)

               LUP    = .TRUE.
               LDOWN  = .TRUE.

               IF (LNFC(1,II).EQ.JJ) THEN
                  LUP    = IUP(1,II)
                  LDOWN  = IDOWN(1,II)
               END IF

               IF (LNFC(2,II).EQ.JJ) THEN
                  LUP    = IUP(2,II)
                  LDOWN  = IDOWN(2,II)
               END IF

************************************************************************
*     first the quadrilateral and line are compared to see if they     *
*     satisfy equalities which allow them to intersect.  if not there  *
*     is no intersection and the next quadrilateral is considered.     *
************************************************************************
               IF ((LUP.OR.LDOWN) .AND. (ZMINL.LT.FMM(5,JJ)) .AND.
     $             (XMAXL.GE.FMM(1,JJ)) .AND. (XMINL.LE.FMM(2,JJ)) .AND.
     $             (YMAXL.GE.FMM(3,JJ)) .AND. (YMINL.LE.FMM(4,JJ))) THEN

                  INDP   = ABS(IFACE1(JJ))
                  QX1    = X(INDP)
                  QY1    = Y(INDP)

                  INDP   = ABS(IFACE2(JJ))
                  QX2    = X(INDP)
                  QY2    = Y(INDP)

                  INDP   = IFACE3(JJ)
                  QX3    = X(INDP)
                  QY3    = Y(INDP)

                  INDP   = IFACE4(JJ)
                  QX4    = X(INDP)
                  QY4    = Y(INDP)

************************************************************************
*     compares a line segment stored in line to a quadrilateral        *
*     stored in quad to find the two concealed intervals.              *
*     a quadrilateral can not conceal one of its own edges.            *
*     the returned intervals are subintervals of 0,1.                  *
*     if no concealed intervals are found, flag1 and flag2 are false   *
*     after return.  if flag1 is true, the interval r1,r2 is concealed.*
*     if flag2 is true, the interval r3,r4 is concealed                *
*     the parallel projection from z = infinity is assumed.            *
*     the line is parameterized from t = 0 at l1 = line(1,.) to t = 1  *
*     at l2 = line(2,.).  the points in the quadrilateral are assumed  *
*     to lie in a clockwise or counterclockwise order q1 = quad(1,.),  *
*     q2 = quad(2,.), q3 = quad(3,.), q4 = quad(4,.).                  *
************************************************************************

                  FLAG1  = .FALSE.
                  FLAG2  = .FALSE.
************************************************************************
*     we compute the intersection of the line with the triangles.      *
************************************************************************

                  V1     = (QY1-LY1)*XDIR - (QX1-LX1)*YDIR
                  V2     = (QY2-LY1)*XDIR - (QX2-LX1)*YDIR
                  V3     = (QY3-LY1)*XDIR - (QX3-LX1)*YDIR
                  V4     = (QY4-LY1)*XDIR - (QX4-LX1)*YDIR

                  LP1    = .FALSE.
                  LP2    = .FALSE.
                  LN1    = .FALSE.
                  LN2    = .FALSE.

                  IF (V1.GT.0.0) THEN
                     INT1   = 1
                     LP1    = .TRUE.
*
                  ELSE IF (V1.LT.0.0) THEN
                     INT1   = -1
                     LN1    = .TRUE.
*
                  ELSE
                     INT1   = 0
                  END IF

                  IF (V2.GT.0.0) THEN
                     INT2   = 1
                     LP1    = .TRUE.
                     LP2    = .TRUE.
*
                  ELSE IF (V2.LT.0.0) THEN
                     INT2   = -1
                     LN1    = .TRUE.
                     LN2    = .TRUE.
*
                  ELSE
                     INT2   = 0
                  END IF

                  IF (V3.GT.0.0) THEN
                     INT3   = 1
                     LP2    = .TRUE.
*
                  ELSE IF (V3.LT.0.0) THEN
                     INT3   = -1
                     LN2    = .TRUE.
*
                  ELSE
                     INT3   = 0
                  END IF

                  IF (V4.GT.0.0) THEN
                     INT4   = 1
                     LP1    = .TRUE.
                     LP2    = .TRUE.
*
                  ELSE IF (V4.LT.0.0) THEN
                     INT4   = -1
                     LN1    = .TRUE.
                     LN2    = .TRUE.
*
                  ELSE
                     INT4   = 0
                  END IF
*
                  LUP    = LUP .AND. (LP1 .AND. LN1)
                  LDOWN  = LDOWN .AND. (LP2 .AND. LN2)

                  IF (LUP .OR. LDOWN) THEN

********************************************************************** *
*     if int1 - int2 is not zero, line crosses the segment determined  *
*     by the points q(.,1) and quad(.,2), a similar statement holds    *
*     for int2 - int3, int3 - int4, int4 - int1 and int2 - int4.       *
*     in case int1 - int2 is not zero we solve the following system    *
*     for s and t using the resulting x and y equations.               *
*        (1 - t) x1 + t x2 = s q1 + (1 - s) q2                         *
*     this becomes                                                     *
*      t(lx2 - lx1) + s(qx2 -qx1) = qx2 - lx1                          *
*      t(ly2 - ly1) + s(qy2 -qy1) = qy2 - ly1                          *
*     thus if                                                          *
*       det =  (lx2 - lx1)(qy2 - qy1) - (ly2 - ly1)(qx2 - qx1)         *
*       tdet = (qx2 - lx1)(qy2 - qy1) - (qy2 - ly1)(qx2 - qx1)         *
*       sdet = (lx2 - lx1)(qy2 - ly1) - (ly2 - ly1)(qx2 - lx1)         *
*       t = tdet/det and s = sdet/det.                                 *
*     in this case the value of the determinant used as the            *
*     denominator when solving the system is not zero.                 *
********************************************************************** *

                     R1     = 1.0
                     R2     = 0.0
                     R3     = 1.0
                     R4     = 0.0

                     LGCL1  = INT1 .NE. INT2
                     LGCL2  = INT1 .NE. INT4
                     LGCL3  = INT2 .NE. INT4
                     LGCL4  = LGCL3
                     LGCL5  = INT3 .NE. INT4
                     LGCL6  = INT2 .NE. INT3
                     IF (LGCL1 .AND. LGCL2 .AND. LGCL3) THEN
                        VV1    = ABS(V1-V2)
                        VV2    = ABS(V1-V4)
                        VV3    = ABS(V2-V4)

                        IF (VV1.GT.VV2) THEN

                           IF (VV2.GT.VV3) THEN
                              LGCL3  = .FALSE.
*
                           ELSE
                              LGCL2  = .FALSE.
                           END IF

                        ELSE
                           IF (VV1.GT.VV3) THEN
                              LGCL3  = .FALSE.
*
                           ELSE
                              LGCL1  = .FALSE.
                           END IF
*
                        END IF
*
                     END IF

                     IF (LGCL4 .AND. LGCL5 .AND. LGCL6) THEN
                        VV4    = ABS(V2-V4)
                        VV5    = ABS(V4-V3)
                        VV6    = ABS(V2-V3)

                        IF (VV4.GT.VV5) THEN

                           IF (VV5.GT.VV6) THEN
                              LGCL6  = .FALSE.
*
                           ELSE
                              LGCL5  = .FALSE.
                           END IF
*
                        ELSE

                           IF (VV4.GT.VV6) THEN
                              LGCL6  = .FALSE.
*
                           ELSE
                              LGCL4  = .FALSE.
                           END IF
*
                        END IF
*
                     END IF

                     IF (LGCL1) THEN

                        DET    = V2 - V1
                        TDET   = (QX2-LX1)* (QY2-QY1) -
     $                           (QY2-LY1)* (QX2-QX1)
                        T1     = TDET/DET
                        R1     = MIN(R1,T1)
                        R2     = MAX(R2,T1)
                     END IF

                     IF (LGCL2) THEN

                        DET    = V4 - V1
                        TDET   = (QX4-LX1)* (QY4-QY1) -
     $                           (QY4-LY1)* (QX4-QX1)
                        T2     = TDET/DET
                        R1     = MIN(R1,T2)
                        R2     = MAX(R2,T2)
                     END IF

                     IF (LGCL3 .OR. LGCL4) THEN

                        DET    = V4 - V2
                        TDET   = (QX4-LX1)* (QY4-QY2) -
     $                           (QY4-LY1)* (QX4-QX2)
                        T3     = TDET/DET
                        IF (LGCL3) THEN
                           R1     = MIN(R1,T3)
                           R2     = MAX(R2,T3)
                        END IF

                        IF (LGCL4) THEN
                           R3     = MIN(R3,T3)
                           R4     = MAX(R4,T3)
                        END IF
*
                     END IF

                     IF (LGCL5) THEN

                        DET    = V4 - V3
                        TDET   = (QX4-LX1)* (QY4-QY3) -
     $                           (QY4-LY1)* (QX4-QX3)
                        T4     = TDET/DET
                        R3     = MIN(R3,T4)
                        R4     = MAX(R4,T4)

                     END IF

                     IF (LGCL6) THEN
                        DET    = V3 - V2
                        TDET   = (QX3-LX1)* (QY3-QY2) -
     $                           (QY3-LY1)* (QX3-QX2)
                        T5     = TDET/DET
                        R3     = MIN(R3,T5)
                        R4     = MAX(R4,T5)
                     END IF

                     R1     = MAX(R1,0.0)
                     R2     = MIN(R2,1.0)
                     R3     = MAX(R3,0.0)
                     R4     = MIN(R4,1.0)
************************************************************************
*     the calculation of intersections is finished. r1,r2,r3, and      *
*     r4 contain the appropriate values in case the z-coordinate       *
*     of a point in the triangle is greater than the corresponding     *
*     point on the line.                                               *
************************************************************************
                     IF ((R2.GT.R1+.0001) .AND. LUP) THEN
                        T      = (R2+R1)*.5
                        ZLINE  = (1.0-T)*LZ1 + T*LZ2
                        XLINE  = (1.0-T)*LX1 + T*LX2
                        YLINE  = (1.0-T)*LY1 + T*LY2
                        EZQUAD = ENRM(1,4,JJ) - ENRM(1,1,JJ)*XLINE -
     $                           ENRM(1,2,JJ)*YLINE
                        EZLINE = ENRM(1,3,JJ)*ZLINE
                        FLAG1  = EZLINE .LT. EZQUAD
                     END IF

                     IF ((R4.GT.R3+.0001) .AND. LDOWN) THEN
                        T      = (R4+R3)*.5
                        ZLINE  = (1.0-T)*LZ1 + T*LZ2
                        XLINE  = (1.0-T)*LX1 + T*LX2
                        YLINE  = (1.0-T)*LY1 + T*LY2
                        EZQUAD = ENRM(2,4,JJ) - ENRM(2,1,JJ)*XLINE -
     $                           ENRM(2,2,JJ)*YLINE

                        EZLINE = ENRM(2,3,JJ)*ZLINE
                        FLAG2  = EZLINE .LT. EZQUAD
                     END IF
***********************************************************************

                     IF (FLAG1) THEN
                        RMAX   = R2
                        RMIN   = R1

                        FLAG1  = .FALSE.
                        CALL GR3FXL(RMAX,RMIN,PAVAIL,AVAIL,PSTART,PRM,
     $                              RM,FLAG1)

                        IF (FLAG1) WRITE (6,FMT=*
     $                      ) 'TOO MANY INTERVALS ON LINE '

                        IF (PSTART.EQ.0) GO TO 70

                     END IF

                     IF (FLAG2) THEN
                        RMAX   = R4
                        RMIN   = R3

                        FLAG1  = .FALSE.
                        CALL GR3FXL(RMAX,RMIN,PAVAIL,AVAIL,PSTART,PRM,
     $                              RM,FLAG1)

                        IF (FLAG1) WRITE (6,FMT=*
     $                      ) 'TOO MANY INTERVALS ON LINE '
                        IF (PSTART.EQ.0) GO TO 70
                     END IF

                  END IF

               END IF

   50       CONTINUE

            LNF1   = LNFC(1,II)
            LNF2   = LNFC(2,II)
            IF (LNF1.GT.0 .AND. LNF2.GT.0) THEN
               I11    = IFACE1(LNF1)
               I12    = IFACE1(LNF2)
               I21    = IFACE2(LNF1)
               I22    = IFACE2(LNF2)
               I31    = IFACE3(LNF1)
               I41    = IFACE4(LNF1)
               I32    = IFACE3(LNF2)
               I42    = IFACE4(LNF2)
               IF (ABS(I11).EQ.ABS(I21)) I41    = ABS(I11)
               IF (ABS(I21).EQ.I31) I41    = I31
               IF (ABS(I12).EQ.ABS(I22)) I42    = ABS(I12)
               IF (ABS(I22).EQ.I32) I42    = I32
               IF (I11.GT.0 .AND. I12.GT.0 .AND. I21.GT.0 .AND.
     $             I22.GT.0) THEN
                  IFA    = ICOL
*
               ELSE IF (I11.LT.0 .AND. I12.LT.0 .AND. I21.LT.0 .AND.
     $                  I22.LT.0) THEN
                  IFA    = ICOLO
*
               ELSE IF (ABS(I11).EQ.ABS(I21) .OR. ABS(I21).EQ.I31) THEN
*
               ELSE
                  DREI1B = ABS(I12) .NE. I42 .AND.
     $                     ABS(I12) .NE. ABS(I22) .AND.
     $                     ABS(I22) .NE. I42
                  DREI2B = I32 .NE. I42 .AND. I32 .NE. ABS(I22) .AND.
     $                     ABS(I22) .NE. I42
                  IF (ILINE(1,II).EQ.ABS(I12)) THEN
                     IF (ILINE(2,II).EQ.I42) THEN
                        IF (DREI1B) THEN
                           IP1B   = ABS(I22)
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                        IF (DREI2B) THEN
                           IP2B   = I32
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.ABS(I22)) THEN
                        IF (DREI1B) THEN
                           IP1B   = I42
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                        IF (DREI2B) THEN
                           IP2B   = I32
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.I32) THEN
                        IF (DREI1B .OR. DREI2B) THEN
                           IP1B   = ABS(I22)
                           IP2B   = IP1B
*
                        ELSE
                           IP1B   = 0
                           IP2B   = 0
                        END IF
*
                     ELSE
                        IP1B   = 0
                        IP2B   = 0
                     END IF
*
                  ELSE IF (ILINE(1,II).EQ.ABS(I22)) THEN
                     IF (ILINE(2,II).EQ.ABS(I12)) THEN
                        IF (DREI1B) THEN
                           IP2B   = I42
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                        IF (DREI2B) THEN
                           IP1B   = I32
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.I32) THEN
                        IF (DREI1B) THEN
                           IP1B   = ABS(I12)
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                        IF (DREI2B) THEN
                           IP2B   = I42
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.I42) THEN
                        IF (DREI1B) THEN
                           IP1B   = ABS(I12)
                           IP2B   = ABS(I12)
*
                        ELSE IF (DREI2B) THEN
                           IP1B   = I32
                           IP2B   = I32
*
                        ELSE
                           IP1B   = 0
                           IP2B   = 0
                        END IF
*
                     ELSE
                        IP1B   = 0
                        IP2B   = 0
                     END IF
*
                  ELSE IF (ILINE(1,II).EQ.I32) THEN
                     IF (ILINE(2,II).EQ.I42) THEN
                        IF (DREI2B) THEN
                           IP1B   = ABS(I22)
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                        IF (DREI1B) THEN
                           IP2B   = ABS(I12)
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.ABS(I22)) THEN
                        IF (DREI2B) THEN
                           IP1B   = I42
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                        IF (DREI1B) THEN
                           IP2B   = ABS(I12)
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.ABS(I12)) THEN
                        IF (DREI1B .OR. DREI2B) THEN
                           IP1B   = ABS(I22)
                           IP2B   = IP1B
*
                        ELSE
                           IP1B   = 0
                           IP2B   = 0
                        END IF
*
                     ELSE
                        IP1B   = 0
                        IP2B   = 0
                     END IF
*
                  ELSE IF (ILINE(1,II).EQ.I42) THEN
                     IF (ILINE(2,II).EQ.ABS(I12)) THEN
                        IF (DREI1B) THEN
                           IP2B   = ABS(I22)
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                        IF (DREI2B) THEN
                           IP1B   = I32
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.I32) THEN
                        IF (DREI1B) THEN
                           IP1B   = ABS(I12)
*
                        ELSE
                           IP1B   = 0
                        END IF
*
                        IF (DREI2B) THEN
                           IP2B   = ABS(I22)
*
                        ELSE
                           IP2B   = 0
                        END IF
*
                     ELSE IF (ILINE(2,II).EQ.ABS(I22)) THEN
                        IF (DREI1B) THEN
                           IP1B   = ABS(I12)
                           IP2B   = ABS(I12)
*
                        ELSE IF (DREI2B) THEN
                           IP1B   = I32
                           IP2B   = I32
*
                        ELSE
                           IP1B   = 0
                           IP2B   = 0
                        END IF
*
                     ELSE
                        IP1B   = 0
                        IP2B   = 0
                     END IF
*
                  END IF
c-------------------------------------------------------------------
                  IFA    = ICOL
                  IF (IP1B.EQ.0 .AND. IP2B.NE.0) IP1B   = IP2B
                  IF (IP2B.EQ.0 .AND. IP1B.NE.0) IP2B   = IP1B
                  IF (IP1B.NE.0 .AND. IP2B.NE.0) THEN
                     X1     = X(ILINE(1,II)) +
     $                        (X(IP1B)-X(ILINE(1,II)))*TEIL
                     Y1     = Y(ILINE(1,II)) +
     $                        (Y(IP1B)-Y(ILINE(1,II)))*TEIL
                     Z1     = Z(ILINE(1,II)) +
     $                        (Z(IP1B)-Z(ILINE(1,II)))*TEIL
                     X2     = X(ILINE(2,II)) +
     $                        (X(IP2B)-X(ILINE(2,II)))*TEIL
                     Y2     = Y(ILINE(2,II)) +
     $                        (Y(IP2B)-Y(ILINE(2,II)))*TEIL
                     Z2     = Z(ILINE(2,II)) +
     $                        (Z(IP2B)-Z(ILINE(2,II)))*TEIL
                     X0     = (X1+X2)/2
                     Y0     = (Y1+Y2)/2
                     Z0     = (Z1+Z2)/2
                     IF (ENRM(1,3,LNF1).NE.0.) THEN
                        W01    = (ENRM(1,4,LNF1)-ENRM(1,1,LNF1)*X0-
     $                           ENRM(1,2,LNF1)*Y0)/ENRM(1,3,LNF1)
*
                     ELSE
                        W01    = -75.75E20
                     END IF
*
                     IF (ENRM(2,3,LNF1).NE.0.) THEN
                        W02    = (ENRM(2,4,LNF1)-ENRM(2,1,LNF1)*X0-
     $                           ENRM(2,2,LNF1)*Y0)/ENRM(2,3,LNF1)
*
                     ELSE
                        W02    = -75.75E20
                     END IF
*
                     IF (W01.EQ.-75.75E20 .OR. W02.EQ.-75.75E20) THEN
                        W0     = MAX(W01,W02)
*
                     ELSE
                        W0     = MIN(W01,W02)
                     END IF
*
                     IF (W0.GT.Z0) THEN
                        IF (I11.LT.0 .OR. I21.LT.0) THEN
                           IFA    = ICOLO
                        END IF
*
                     ELSE
                        IF (I12.LT.0 .AND. I22.LT.0) THEN
                           IFA    = ICOLO
                        END IF
*
                     END IF
*
                  ELSE
                     IFA    = ICOLO
                  END IF
*
               END IF
*
            ELSE IF (LNF1.GT.0) THEN
               I11    = IFACE1(LNF1)
               I21    = IFACE2(LNF1)
               IF (I11.GT.0 .OR. I21.GT.0) THEN
                  IFA    = ICOL
*
               ELSE IF (I11.LT.0 .AND. I21.LT.0) THEN
                  IFA    = ICOLO
*
               ELSE
                  IFA    = ICOL
               END IF
*
            END IF
c----- do while (pstart.ne.0) ------------------------------------------
   60       CONTINUE
            IF (PSTART.NE.0) THEN
               TEMP1  = RM(1,PSTART)
               X1     = (1.-TEMP1)*LX1 + TEMP1*LX2
               Y1     = (1.-TEMP1)*LY1 + TEMP1*LY2
               TEMP2  = RM(2,PSTART)
               X2     = (1.-TEMP2)*LX1 + TEMP2*LX2
               Y2     = (1.-TEMP2)*LY1 + TEMP2*LY2
               IF (((X2-X1)**2+ (Y2-Y1)**2).GT.DIAG2) THEN
                  IF (IFA.NE.IFAOLD) THEN
                     CALL GRNWPN(IFA)
                     IFAOLD = IFA
                  END IF
*
                  CALL GRJMP(X1,Y1)
                  CALL GRDRW(X2,Y2)
               END IF
c  the visible segment found above is (x1,y1) (x2,y2)
c  locations for visible intervals stored are returned to avail below.

               PAVAIL = PAVAIL + 1
               AVAIL(PAVAIL) = PSTART
               PSTART = PRM(PSTART)
               GO TO 60
*
            END IF

   70    CONTINUE
         SEE    = JGR(5,ISURF) .GT. 0
   80 CONTINUE
      CALL GRSPTS(INTOL)
      END
C=======================================================================
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3PLT(KIND,X,Y,Z,NP,IFACE1,IFACE2,IFACE3,IFACE4,MNF,
     $                  INDF,ILINE,LNFC,MNL,INDL,WORK,LENW,RA4,KINDAX,
     $                  IER)





c     .. parameters ..
      INTEGER           MFL
      PARAMETER         (MFL=256)
c     ..
c     .. scalar arguments ..
      INTEGER           IER,INDF,INDL,KIND,KINDAX,LENW,MNF,MNL,NP
c     ..
c     .. array arguments ..
      REAL              RA4(4),WORK(*),X(NP),Y(NP),Z(NP)
      INTEGER           IFACE1(MNF),IFACE2(MNF),IFACE3(MNF),IFACE4(MNF),
     $                  ILINE(2,MNL),LNFC(2,MNL)
c     ..
c     .. local scalars ..
      REAL              AB,ABMA,D1,D2,D3,FAXL,FAXR,FLGROT,RA41,RA42,
     $                  RA43,RA44,SCX,SCY,SPATEX,SQ,U,WIFA,X0,X1,XADD,
     $                  XFA,XMAS,XMAX,XMIN,XMIT,XN0,XN1,XR,XS,XTEX,Y0,
     $                  Y1,YADD,YFA,YMAS,YMAX,YMIN,YMIT,YN0,YN1,YR,YS,
     $                  YTEX,ZEIGRO,ZEIWI,ZEIWI0,ZEIWIN
      INTEGER           I,I1,I2,IE,IKO,ILI,INTLIN,INTOLD,INTSYM,SIZMRK,
     $                  ISURF,J1,J2,JGI,JGIK,K,K1,K2,KMA,L
      LOGICAL           LINKS,SEE
      CHARACTER*10      C10
c     ..
c     .. local arrays ..
      REAL              EX4(4),RAX(5),RAY(5)
      INTEGER           IAX(3),LCH(3)
c     ..
c     .. external functions ..
      INTEGER           GRISMA,GRISMI
      EXTERNAL          GRISMA,GRISMI
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3ALL,GR3HID,GRAXLIN,GRAXLOG,GRCHRC,GRDRW,
     $                  GRDSH,GRFTOC,GRJMP,GRLN,GRNWPN,GRSCLC,GRSCLV,
     $                  GRSPTS,GSCHH,GSCHSP,GSCHUP,GSCHXP,GTX
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ATAN,ATAN2,MIN,SQRT
c     ..
c     .. common blocks ..
      COMMON            /GR3SAC/CAXS,CORI
CDEC$ PSECT /GR3SAC/ NOSHR
      COMMON            /GR4COM/IFL,JGR,XC,YC,COMA
CDEC$ PSECT /GR4COM/ NOSHR
      COMMON            /GRPIC/FLPIC,NSCLC,NSCLV,NSCLP,RAHMEN,XMXCMQ,
     $                  YMXCMQ,XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON            /GRPP/PP
CDEC$ PSECT /GRPP/ NOSHR

C---- COMMON-Bloecke zu GR3ANT
      INTEGER            NUMAAN, LEMAAN
      PARAMETER          (NUMAAN=32,LEMAAN=32)
      INTEGER            NUA,INDA,LENA,IFOA,ICOA,IBOA
      REAL               SIZA,DS1A,DS2A,DS3A,ANGA,PTXA,PTYA
      CHARACTER*(LEMAAN) CHAA
      COMMON /GRANTN/ NUA,INDA(NUMAAN),LENA(NUMAAN),SIZA(NUMAAN),
     $       IFOA(NUMAAN),DS1A(NUMAAN),DS2A(NUMAAN),DS3A(NUMAAN),
     $       ANGA(NUMAAN),PTXA(NUMAAN),PTYA(NUMAAN),ICOA(NUMAAN),
     $       IBOA(NUMAAN)
CDEC$ PSECT /GRANTN/ NOSHR
      COMMON /GRANTC/ CHAA(NUMAAN)
CDEC$ PSECT /GRANTC/ NOSHR

      REAL              XDCPIC,XMXCMQ,YDCPIC,YMXCMQ
      INTEGER           FLPIC,IFL,NSCLC,NSCLP,NSCLV,RAHMEN
      CHARACTER         CORI
      REAL              COMA(3,2),PP(18),XC(500),YC(500)
      INTEGER           JGR(7,MFL)
      CHARACTER*20      CAXS(3)
c     ..
c     .. equivalences ..
      EQUIVALENCE       (PP(14),ZEIGRO), (PP(15),ZEIWIN),
     $                  (PP(16),INTSYM), (PP(13),INTLIN),
     $                  (PP(17),SIZMRK), (PP(18),FLGROT)
c     ..
c     .. save statement ..
      SAVE /GR3SAC/,/GR4COM/,/GRPIC/, /GRPP/,/GRANTN/ ,/GRANTC/
      SAVE
c     ..
c     .. executable statements ..
c-----------------------------------------------------------------------


c  the minimum and maximum x- and y-coordinates of the projection
c  are found.

      EX4(1) = X(GRISMI(NP,X,1))
      EX4(2) = Y(GRISMI(NP,Y,1))
      EX4(3) = X(GRISMA(NP,X,1))
      EX4(4) = Y(GRISMA(NP,Y,1))
      RA41   = RA4(1)
      RA42   = RA4(2)
      RA43   = RA4(3)
      RA44   = RA4(4)
      IE     = 0
      IF (RA41.EQ.RA43 .AND. RA42.EQ.RA44) THEN
         RA41   = EX4(1)
         RA42   = EX4(2)
         RA43   = EX4(3)
         RA44   = EX4(4)
c------- zooming in y-richtung ?
      ELSE IF (RA42.EQ.RA44) THEN
         YMIT   = RA42
         RA42   = YMIT - (RA43-RA41)*0.3632911
         RA44   = YMIT + (RA43-RA41)*0.3632911
c------- zooming in x-richtung ?
      ELSE IF (RA41.EQ.RA43) THEN
         XMIT   = RA4(1)
         RA41   = XMIT - (RA44-RA42)*0.6881533
         RA43   = XMIT + (RA44-RA42)*0.6881533
*
      ELSE
         IF (EX4(1).LT.RA41 .OR. EX4(2).LT.RA42 .OR. EX4(3).GT.RA43 .OR.
     $       EX4(4).GT.RA44) IE = 40
      END IF

      XADD   = .01* (RA43-RA41)
      YADD   = .01* (RA44-RA42)
      XMAX   = RA43 + XADD
      XMIN   = RA41 - XADD
      YMAX   = RA44 + YADD
      YMIN   = RA42 - YADD

c  hole alte werte, bestimme fenstergroesse

      X0     = PP(1)
      Y0     = PP(2)
      X1     = PP(3)
      Y1     = PP(4)
      PP5 = PP(5)
      PP6 = PP(6)
      PP7 = PP(7)
      PP8 = PP(8)
      ZEIGR0 = ZEIGRO
      ZEIWI0 = ZEIWIN
      CALL GQTXFP(IERR,IFONT,IPREC)
      IF ((YMAX-YMIN)* (X1-X0).GT. (XMAX-XMIN)* (Y1-Y0)) THEN
         YN0    = Y0
         YN1    = Y1
         U      = (XMAX-XMIN)* (Y1-Y0)/ (YMAX-YMIN)
         XN0    = (X1+X0-U)*0.5
         XN1    = XN0 + U
*
      ELSE
         XN0    = X0
         XN1    = X1
         U      = (YMAX-YMIN)* (X1-X0)/ (XMAX-XMIN)
         YN0    = (Y1+Y0-U)*0.5
         YN1    = YN0 + U
      END IF
*
      CALL GRSCLC(XN0,YN0,XN1,YN1)
      CALL GRSCLV(XMIN,YMIN,XMAX,YMAX)
      XFA    = (PP(7)-PP(5))/ (PP(3)-PP(1))
      YFA    = (PP(8)-PP(6))/ (PP(4)-PP(2))
      IF (IE.EQ.40) THEN
         XMAS   = 0.001* (XMAX-XMIN)
         YMAS   = 0.001* (YMAX-YMIN)
         RAX(1) = XMIN + XMAS
         RAY(1) = YMIN + YMAS
         RAX(2) = XMAX - XMAS
         RAY(2) = YMIN + YMAS
         RAX(3) = XMAX - XMAS
         RAY(3) = YMAX - YMAS
         RAX(4) = XMIN + YMAS
         RAY(4) = YMAX - YMAS
         RAX(5) = XMIN + XMAS
         RAY(5) = YMIN + YMAS
         CALL GRNWPN(1)
         CALL GRLN(RAX,RAY,5)
      END IF

c---- soll ein koordinatenkubus gezeichnet werden?

      DO 10 I = 1,IFL
         IF (JGR(1,I).EQ.-32767) GO TO 20
   10 CONTINUE
      GO TO 110

c---- bestimme die beschriftungslaenge 'lch'
   20 CONTINUE
      CALL GRNWPN(JGR(3,I))
      CALL GRSPTS(JGR(5,I))
      INTOLD = JGR(5,I)
      SEE    = INTOLD .GT. 0
      IF (KINDAX.EQ.1 .OR. KINDAX.EQ.3 .OR. KINDAX.EQ.5) JGR(5,I) = 0
      DO 50 K = 1,3
         DO 30 L = 20,2,-1
            IF (CAXS(K) (L:L).NE.' ') GO TO 40
   30    CONTINUE
   40    CONTINUE
         LCH(K) = L
   50 CONTINUE

c---- berechne den schwerpunkt des koordinatenkubus

      SCX    = 0.
      SCY    = 0.
      JGI    = 0
      IF (I.GT.1) JGI    = JGR(4,I-1)
      DO 70 K = 1,4
         JGIK   = JGI + K
         DO 60 L = 1,2
            ILI    = ILINE(L,JGIK)
            SCX    = SCX + X(ILI)
            SCY    = SCY + Y(ILI)
   60    CONTINUE
   70 CONTINUE
      SCX    = SCX*.125
      SCY    = SCY*.125

c---- suche in allen drei richtungen die am weitesten vom zentrum
c---- entfernten achsen.

      DO 90 IKO = 1,3
         ABMA   = 0.
         DO 80 K = 1,4
            JGIK   = JGI + (IKO-1)*4 + K
            I1     = ILINE(1,JGIK)
            I2     = ILINE(2,JGIK)
            AB     = MIN(SQRT((X(I1)-SCX)**2+ (Y(I1)-SCY)**2),
     $               SQRT((X(I2)-SCX)**2+ (Y(I2)-SCY)**2))
C-update--- 11.11.91  1.000002 => 1.000004
            IF (AB.GE.ABMA*1.000004) THEN
               ABMA   = AB
               KMA    = K
            END IF
*
   80    CONTINUE
         IAX(IKO) = KMA
   90 CONTINUE

c---- schreibe die texte aus caxs an diese achsen

      WIFA   = 45./ATAN(1.)
      DO 100 IKO = 1,3
         JGIK   = JGI + (IKO-1)*4 + IAX(IKO)
         I1     = ILINE(1,JGIK)
         I2     = ILINE(2,JGIK)
         IF ((KINDAX.EQ.1.OR.KINDAX.EQ.3.OR.KINDAX.EQ.5) .AND. SEE) THEN
            CALL GRJMP(X(I1),Y(I1))
            CALL GRDRW(X(I2),Y(I2))
         END IF
*
         XR     = X(I2) - X(I1)
         YR     = Y(I2) - Y(I1)
c------- damit keine schrift auf dem kopf steht
         IF (XR.GE.0. .OR. CORI.EQ.'1') THEN
            J1     = I1
            J2     = I2
            K1     = 1
            K2     = 2
*
         ELSE
            J1     = I2
            J2     = I1
            XR     = -XR
            YR     = -YR
            K1     = 2
            K2     = 1
         END IF
*
         IF (XR.NE.0. .OR. YR.NE.0.) THEN
            ZEIWI  = ATAN2(YR,XR)*WIFA
            SQ     = SQRT(XR**2+YR**2)
            XR     = XR/SQ
            YR     = YR/SQ
c---------- liegt der zentralpunkt links oder rechts? ---kreuzprodukt--
            IF (XR* (SCY-Y(J1))-YR* (SCX-X(J1)).LE.0.) THEN
               LINKS  = .TRUE.
               XS     = XR
               YS     = YR
*
            ELSE
               LINKS  = .FALSE.
               XS     = -3.*XR
               YS     = -3.*YR
            END IF
*
            IF (KINDAX.GT.2 .AND. COMA(IKO,1).NE.COMA(IKO,2)) THEN
               IF (IKO.LT.3 .OR. KINDAX.EQ.3 .OR. KINDAX.EQ.4) THEN
                  CALL GRAXLIN(X(J1),Y(J1),X(J2),Y(J2),COMA(IKO,K1),
     $                         COMA(IKO,K2),LINKS,0)
                  FAXL   = 9
                  FAXR   = 7
*
               ELSE
                  CALL GRAXLOG(X(J1),Y(J1),X(J2),Y(J2),COMA(IKO,K1),
     $                         COMA(IKO,K2),LINKS,0)
                  FAXL   = 11
                  FAXR   = 9
               END IF
*
               IF (LINKS) THEN
                  XS     = XS + XR*FAXL
                  YS     = YS + YR*FAXL
*
               ELSE
                  XS     = XS - XR*FAXR
                  YS     = YS - YR*FAXR
               END IF
*
            END IF
c mvs,cray  call grchrc(zeigro,zeiwi,intsym)
c next cms
            CALL GSCHUP(-YR,XR)
            CALL GSCHH(ZEIGRO*XFA)
            CALL GSCHSP(0.)
            CALL GSCHXP(MIN(SQ/ (LCH(IKO)*ZEIGRO*XFA),1.2))
c end cms
            SPATEX = LCH(IKO)*ZEIGRO
            XTEX   = .5* (X(J1)+X(J2)- (SPATEX*XR+ZEIGRO*YS*1.25)*XFA)
            YTEX   = .5* (Y(J1)+Y(J2)- (SPATEX*YR-ZEIGRO*XS*1.25)*YFA)
c mvs cray  if (see .and. lch(iko).gt.0)
c    >         call grtxt(xtex,ytex,lch(iko),caxs(iko))
c next cms
            IF (SEE .AND. LCH(IKO).GT.0) CALL GTX(XTEX,YTEX,
     $          CAXS(IKO) (:LCH(IKO)))
            CALL GSCHXP(1.)
c end cms
            IF (KINDAX.LE.2 .AND. COMA(IKO,1).NE.COMA(IKO,2)) THEN
               CALL GRFTOC(COMA(IKO,K1),C10,L)
               XTEX   = X(J1) + (2.*XR-YS*.5)*ZEIGRO*XFA
               YTEX   = Y(J1) + (2.*YR+XS*.5)*ZEIGRO*YFA
c mvs cray     if (see) call grtxt(xtex,ytex,l,c10)
c next cms
               IF (SEE) CALL GTX(XTEX,YTEX,C10(:L))
c end cms
               CALL GRFTOC(COMA(IKO,K2),C10,L)
               XTEX   = X(J2) - ((L+0)*XR+YS*.5)*XFA*ZEIGRO
               YTEX   = Y(J2) - ((L+0)*YR-XS*.5)*YFA*ZEIGRO
c mvs cray     if (see) call grtxt(xtex,ytex,l,c10)
c next cms
               IF (SEE) CALL GTX(XTEX,YTEX,C10(:L))
c end cms
            END IF
*
         END IF
*
  100 CONTINUE
  110 CONTINUE
      D1     = PP(10)
      D2     = PP(11)
      D3     = PP(12)
c----------------------------------------------------------------------
c     kind = 1 if all lines should be shown.
c            2 if hidden lines should be dotted.
c            3 if hidden lines are not shown.
c-----------------------------------------------------------------------
      CALL GRDSH(1.,0.,1.)
      IF (KIND.EQ.1 .OR. INDF.EQ.0) THEN
         CALL GR3ALL(X,Y,ILINE,EX4)
*
      ELSE
         IF (KIND.EQ.2) THEN
            DO 120 ISURF = 1,IFL
               IF (JGR(5,ISURF).NE.0) JGR(5,ISURF) = JGR(5,ISURF) + 5
  120       CONTINUE
         END IF
*
         CALL GR3HID(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,MNF,INDF,ILINE,
     $               LNFC,MNL,INDL,EX4,WORK,LENW,IER,NP)
         IF (KIND.EQ.2) THEN
            DO 130 ISURF = 1,IFL
               IF (JGR(5,ISURF).NE.0) JGR(5,ISURF) = JGR(5,ISURF) - 5
  130       CONTINUE
            CALL GRDSH(0.1,0.5,0.1)
            CALL GR3ALL(X,Y,ILINE,EX4)
            CALL GRDSH(1.,0.,1.)
         END IF
*
      END IF
*
      JGR(5,I) = INTOLD
*
      DO 200 I=1,NUA
         ICOLA=MOD(ICOA(I),2000)
         CALL GRNWPN(ICOLA)
         CALL GRDSH(DS1A(I),DS2A(I),DS2A(I))
         IDICK=18+(ICOA(I)-ICOLA)/2000*2
         CALL GRSPTS(IDICK)
         CALL GRJMP(X(INDA(I)),Y(INDA(I)))
         XT=X(INDA(I))+PTXA(I)*XFA
         YT=Y(INDA(I))+PTYA(I)*YFA
         CALL GRDRW(XT,YT)
         CALL GRCHRC(SIZA(I),ANGA(I),18)
         IF ( IBOA(I).EQ.2 ) THEN
            IF ( PTYA(I).GE.0. ) THEN
               IBO2=5
            ELSE
               IBO2=1
            ENDIF
         ELSE
            IBO2=3
         ENDIF
         CALL GSTXAL(IBOA(I),IBO2)
         CALL GRFONT(IFOA(I))
         CALL GRTXT(XT,YT,LENA(I),CHAA(I))
  200 CONTINUE
      CALL GRDSH(D1,D2,D3)
      CALL GSTXAL(0,4)
      CALL GRCHRC(ZEIGR0,ZEIWI0,INTSYM)
      CALL GRFONT(IFONT)
      CALL GRSCLC(X0,Y0,X1,Y1)
      CALL GRSCLV(PP5,PP6,PP7,PP8)
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3ALL(X,Y,ILINE,EX4)

c     plot of all lines



c     commonblock der standardwerte bzw. der geaenderten tabellenwerte

c     .. parameters ..
      INTEGER           MFL
      PARAMETER         (MFL=256)
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IFL,IL,IL0,IP0,ISY,ISZ,IW,
     $                  LQ,LW,MAFL,MNF,MNL,MNP,NOTDIM
      REAL              EX4(4),X(MNP),Y(MNP)
      INTEGER           ILINE(2,MNL)
c     ..
c     .. local scalars ..
      REAL              D,FLGROT,X1,X2,Y1,Y2,ZEIGRO,ZEIWIN
      INTEGER           I,II,INTLIN,INTOL,INTSYM,SIZMRK,K,L,L1
      LOGICAL           NWPN,SEE
c     ..
c     .. external subroutines ..
      EXTERNAL          GRLN,GRNWPN,GRSPTS
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MAX
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      COMMON            /GR4COM/IFL,JGR,XC,YC,COMA
CDEC$ PSECT /GR4COM/ NOSHR
      COMMON            /GRPP/PP
CDEC$ PSECT /GRPP/ NOSHR
      REAL              WIN,ZENPRO
      REAL              COMA(3,2),GR(3,3),PP(18),RT(3,3),XC(500),YC(500)
      INTEGER           JGR(7,MFL)
c     ..
c     .. equivalences ..
      EQUIVALENCE       (PP(14),ZEIGRO), (PP(15),ZEIWIN),
     $                  (PP(16),INTSYM), (PP(13),INTLIN),
     $                  (PP(17),SIZMRK), (PP(18),FLGROT)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/,/GR4COM/,/GRPP/
      SAVE
c     ..
c     .. executable statements ..

      D      = MAX(EX4(3)-EX4(1),EX4(4)-EX4(2))*0.0001
      L      = 1
      II     = 0
      X2     = X(ILINE(1,1))
      Y2     = Y(ILINE(1,1))
      CALL GRNWPN(JGR(3,1))
      INTOL  = INTLIN
      SEE    = JGR(5,1) .GT. 0
      DO 20 K = 1,IFL
         CALL GRSPTS(JGR(5,K))
         NWPN   = K .GT. 1
         L1     = JGR(4,K)
         DO 10 I = L,L1
            X1     = X(ILINE(1,I))
            Y1     = Y(ILINE(1,I))
            II     = II + 1
            XC(II) = X2
            YC(II) = Y2
            IF (ABS(X2-X1).GE.D .OR. ABS(Y2-Y1).GE.D .OR. II.EQ.500 .OR.
     $          NWPN) THEN
               IF (II.GT.1 .AND. SEE) CALL GRLN(XC,YC,II)
               IF (NWPN) THEN
                  CALL GRNWPN(JGR(3,K))
                  NWPN   = .FALSE.
                  SEE    = JGR(5,K) .GT. 0
               END IF
*
               II     = 1
               XC(1)  = X1
               YC(1)  = Y1
            END IF
*
            X2     = X(ILINE(2,I))
            Y2     = Y(ILINE(2,I))
   10    CONTINUE
         L      = L1 + 1
   20 CONTINUE
      II     = II + 1
      XC(II) = X2
      YC(II) = Y2
      IF (SEE) CALL GRLN(XC,YC,II)
      CALL GRSPTS(INTOL)
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3HID(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,MNF,INDF,
     $                  ILINE,LNFC,MNL,INDL,EX4,WORK,LENW,IER,MNP)

c  x,y,z           is an array-group of three dimensional points.
c  iface 1 2 3 4   is an array-group of integers which indicate
c                  the points in x,y,z which belong to the indf
c                  faces under consideration.
c  indf            is the current number of faces in iface. if faces are
c                  actually triangles, one point can be entered twice.
c  iline(2,mnl)    is an array which indicates the endpoints of line
c                  segments and mnl is the maximum number of lines
c                  which will be considered.
c  indl            indicates the current number of lines under
c                  consideration. in practice indl should be somewhat
c                  less than mnl since some lines may be added to
c                  improve the picture when faces are non-planar.
c                  the points for iline are in x,y,z.
c  lnfc(2,mnl)     is an array which indicates the indicies of the faces
c                  to which a line belongs.  zero is entered for
c                  nonexistant faces. if lnfc(1,i)=j and lnfc(2,i)=k,
c                  the i-th line will not be compared to the j-th or
c                  k-th face and the points making up line i are assumed
c                  to be among points numbered 1,2, and 4 in face j and
c                  2,3, and 4 in face k. the last is irrelevant unless
c                  the normals of the trianges formed by 1,2,4 and 2,3,4
c                  have z-coordinate of different sign.
c  work(lenw)      working space with length at least 46*mnp



c  hdden3 compares faces to lines and sends visible segments to printer.
c  if a line is compared to a face on which it lies, it may dissappear
c  due to roundoff error.

c     it would be fairly easy to eliminate one half of the
c     faces and the edges in case the surface to be represented
c     is closed.  that would cause a speedup by a factor of 4.


c---- integer iface1(mnf),iface2(mnf),iface3(mnf),iface4(mnf)
c---- integer iline(2,mnl),lnfc(2,mnl)


c     .. scalar arguments ..
      INTEGER           IER,INDF,INDL,LENW,MNF,MNL,MNP
c     ..
c     .. array arguments ..
      REAL              EX4(4),WORK(LENW),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IFACE1(*),IFACE2(*),IFACE3(*),IFACE4(*),
     $                  ILINE(2,*),LNFC(2,*)
c     ..
c     .. local scalars ..
      REAL              XIDONE,XMNF
      INTEGER           IDONE,IW,IW1,IW2,IW3,JW,LW,LW1,LW2,NR
      LOGICAL           FLAG
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3HDD,GR3MMF,GR3SOR
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MAX,SQRT
c     ..
c     .. save statement ..
c     ..
c     .. executable statements ..

      IF (LENW.LT.23*MNF+4*MNL) THEN
         IF (IER.NE.2) WRITE (6,FMT=*)
     $       'GR3HID RC=10: WORK IS NOT BIG ENOUGH'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 20
*
      ELSE
         CALL GR3MMF(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,WORK,INDF)

         XMNF   = INDF
         NR     = MAX(SQRT(XMNF*0.5),1.0)
         IW     = 13*MNF + 1
         IW1    = IW + NR**2 + 2
         IW2    = IW1 + 4*MNF
c---- do until( .not. flag ) ------------------------------------------
   10    CONTINUE
         FLAG   = .FALSE.

         CALL GR3SOR(NR,WORK,INDF,WORK(IW2),WORK(IW1),WORK(IW),IDONE,
     $               FLAG,MNF,EX4)

         IF (FLAG) THEN
            XIDONE = IDONE
            XMNF   = MNF
            NR     = .75*NR*SQRT(XIDONE/XMNF)
            IF (NR.EQ.0) THEN
               WRITE (*,FMT=*) 'NR .EQ. 0 IN GR3HID'
               STOP
*
            END IF
*
            GO TO 10
*
         END IF
c---- end do until -----------------------------------------------------
         IW3    = IW2 + MNF
         LW     = 1 + 22*MNF
         LW1    = LW + 2*MNL
         LW2    = LW1 + 2*MNL
         JW     = 1 + 5*MNF
         CALL GR3HDD(NR,X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,ILINE,INDF,
     $               INDL,LNFC,MNP,MNL,MNF,EX4,WORK,WORK(JW),WORK(LW),
     $               WORK(LW1),WORK(LW2),WORK(IW),WORK(IW1),WORK(IW2),
     $               WORK(IW3))
         IER    = 0
      END IF

   20 CONTINUE
      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3MMF(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,FMM,INDF)

c finds maximum and minimum coordinates for faces.
c the entries in fmm(j,i) are the min x for face i if j = 1,
c the max x for the face if j = 2, the min y if j = 3, the max y
c if j = 4, and the max z if j = 5.



c     .. scalar arguments ..
      INTEGER           INDF
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              FMM(5,MNF),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IFACE1(MNF),IFACE2(MNF),IFACE3(MNF),IFACE4(MNF)
c     ..
c     .. local scalars ..
      REAL              QX1,QX2,QX3,QX4,QY1,QY2,QY3,QY4,QZ1,QZ2,QZ3,QZ4,
     $                  TMAX1,TMAX2,TMIN1,TMIN2
      INTEGER           I,IND1,IND2,IND3,IND4
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MAX,MIN
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      DO 10 I = 1,INDF
         IND1   = ABS(IFACE1(I))
         IND2   = ABS(IFACE2(I))
         IND3   = IFACE3(I)
         IND4   = IFACE4(I)

         QZ1    = Z(IND1)
         QZ2    = Z(IND2)
         QZ3    = Z(IND3)
         QZ4    = Z(IND4)

         FMM(5,I) = MAX(MAX(QZ1,QZ2),MAX(QZ3,QZ4))

         QX1    = X(IND1)
         QX2    = X(IND2)
         QX3    = X(IND3)
         QX4    = X(IND4)

         TMAX1  = QX2
         IF (QX1.GE.QX2) TMAX1  = QX1
         TMIN1  = QX1
         IF (QX1.GE.QX2) TMIN1  = QX2

         TMAX2  = QX4
         IF (QX3.GE.QX4) TMAX2  = QX3
         TMIN2  = QX3
         IF (QX3.GE.QX4) TMIN2  = QX4

         FMM(1,I) = MIN(TMIN1,TMIN2)
         FMM(2,I) = MAX(TMAX1,TMAX2)

         QY1    = Y(IND1)
         QY2    = Y(IND2)
         QY3    = Y(IND3)
         QY4    = Y(IND4)

         TMAX1  = QY2
         IF (QY1.GE.QY2) TMAX1  = QY1
         TMIN1  = QY1
         IF (QY1.GE.QY2) TMIN1  = QY2

         TMAX2  = QY4
         IF (QY3.GE.QY4) TMAX2  = QY3
         TMIN2  = QY3
         IF (QY3.GE.QY4) TMIN2  = QY4

         FMM(4,I) = MAX(TMAX1,TMAX2)
         FMM(3,I) = MIN(TMIN1,TMIN2)
   10 CONTINUE

      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3NRM(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,ENRM,INDF,
     $                  DIFF)

c      finds the normal equations for the two trianges composing
c      each of the faces.

***********************************************



c     .. scalar arguments ..
      INTEGER           INDF
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              ENRM(2,4,INDF),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IFACE1(INDF),IFACE2(INDF),IFACE3(INDF),
     $                  IFACE4(INDF)
      LOGICAL           DIFF(0:INDF)
c     ..
c     .. local scalars ..
      REAL              TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6
      INTEGER           I,IND1,IND2,IND3,IND4
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..
***********************************************
      DO 10 I = 1,INDF
         IFACE1(I) = ABS(IFACE1(I))
         IFACE2(I) = ABS(IFACE2(I))
         IND1   = IFACE1(I)
         IND2   = IFACE2(I)
         IND3   = IFACE3(I)
         IND4   = IFACE4(I)
         TEMP1  = (Y(IND2)-Y(IND1))* (Z(IND4)-Z(IND2)) -
     $            (Y(IND4)-Y(IND2))* (Z(IND2)-Z(IND1))

         TEMP2  = (X(IND2)-X(IND1))* (Z(IND4)-Z(IND2)) -
     $            (X(IND4)-X(IND2))* (Z(IND2)-Z(IND1))

         TEMP3  = (X(IND2)-X(IND1))* (Y(IND4)-Y(IND2)) -
     $            (X(IND4)-X(IND2))* (Y(IND2)-Y(IND1))
         TEMP4  = (Y(IND3)-Y(IND2))* (Z(IND4)-Z(IND2)) -
     $            (Y(IND4)-Y(IND2))* (Z(IND3)-Z(IND2))

         TEMP5  = (X(IND3)-X(IND2))* (Z(IND4)-Z(IND2)) -
     $            (X(IND4)-X(IND2))* (Z(IND3)-Z(IND2))

         TEMP6  = (X(IND3)-X(IND2))* (Y(IND4)-Y(IND2)) -
     $            (X(IND4)-X(IND2))* (Y(IND3)-Y(IND2))

         ENRM(1,1,I) = TEMP1
         IF (TEMP3.LT.0.0) THEN
            ENRM(1,1,I) = -TEMP1
            IFACE1(I) = -IND1
            IF (TEMP6.EQ.0.) IFACE2(I) = -IND2
         END IF
*
         ENRM(1,2,I) = -TEMP2
         IF (TEMP3.LT.0.0) ENRM(1,2,I) = TEMP2
         ENRM(1,3,I) = TEMP3
         IF (TEMP3.LT.0.0) ENRM(1,3,I) = -TEMP3
         ENRM(2,1,I) = TEMP4
         IF (TEMP6.LT.0.0) THEN
            ENRM(2,1,I) = -TEMP4
            IFACE2(I) = -IND2
            IF (TEMP3.EQ.0.) IFACE1(I) = -IND1
         END IF
*
         ENRM(2,2,I) = -TEMP5
         IF (TEMP6.LT.0.0) ENRM(2,2,I) = TEMP5
         ENRM(2,3,I) = TEMP6
         IF (TEMP6.LT.0.0) ENRM(2,3,I) = -TEMP6

         DIFF(I) = TEMP3*TEMP6 .LT. 0.0

         ENRM(1,4,I) = ENRM(1,1,I)*X(IND2) + ENRM(1,2,I)*Y(IND2) +
     $                 ENRM(1,3,I)*Z(IND2)

         ENRM(2,4,I) = ENRM(2,1,I)*X(IND2) + ENRM(2,2,I)*Y(IND2) +
     $                 ENRM(2,3,I)*Z(IND2)

   10 CONTINUE

      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3BUC(EX4,X0,Y0,X1,Y1,NUM,LIND,NX,NY)

c     gets the indexes of the rectangles (buckets) met by the line
c     segment between (x0,y0) and (x1,y1).




c     .. scalar arguments ..
      REAL              X0,X1,Y0,Y1
      INTEGER           NUM,NX,NY
c     ..
c     .. array arguments ..
      REAL              EX4(4)
      INTEGER           LIND(43,2)
c     ..
c     .. local scalars ..
      REAL              TEMP,TI,XDIV,XINEW,XIOLD,XNULL,XNY,XT,XX0,XX1,
     $                  YDIV,YINEW,YIOLD,YNULL,YT,YY0,YY1
      INTEGER           ITEMP,IX,IXBOT,IXT,IXTOP,IY,IYBOT,IYTOP
      LOGICAL           LGCL
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MIN
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      XDIV   = NX/ (EX4(3)-EX4(1))
      YDIV   = NY/ (EX4(4)-EX4(2))
      XNY    = NY

      IF (X0.LT.X1) THEN
         XX0    = X0
         XX1    = X1
         YY0    = Y0
         YY1    = Y1
*
      ELSE
         XX0    = X1
         XX1    = X0
         YY0    = Y1
         YY1    = Y0
      END IF

      XNULL  = (XX0-EX4(1))*XDIV
      XT     = (XX1-EX4(1))*XDIV

      YNULL  = (YY0-EX4(2))*YDIV
      YT     = (YY1-EX4(2))*YDIV

      LGCL   = YNULL .LE. YT

      XIOLD  = XNULL
      YIOLD  = YNULL

      ITEMP  = XNULL
      TEMP   = ITEMP

      IXBOT  = ITEMP + 1

      ITEMP  = XT
      TEMP   = ITEMP

      IF (TEMP.LT.XT) THEN
         IXTOP  = ITEMP
*
      ELSE
         IXTOP  = ITEMP - 1
      END IF

      IXTOP  = MIN(IXTOP,NX-1)

      NUM    = 0

      DO 20 IX = IXBOT,IXTOP
         XINEW  = IX
         TI     = (XINEW-XNULL)/ (XT-XNULL)
         YINEW  = (1.-TI)*YNULL + TI*YT

         IF (LGCL) THEN
            IYBOT  = MIN(YIOLD+1.,XNY)
            IYTOP  = MIN(YINEW+1.,XNY)
*
         ELSE
            IYTOP  = MIN(YIOLD+1.,XNY)
            IYBOT  = MIN(YINEW+1.,XNY)
         END IF

         XIOLD  = XINEW
         YIOLD  = YINEW

         DO 10 IY = IYBOT,IYTOP
            NUM    = NUM + 1
            LIND(NUM,1) = IX
            LIND(NUM,2) = IY
   10    CONTINUE

   20 CONTINUE

      IXT    = XT
      TEMP   = IXT

      IF (TEMP.LT.XT) THEN
         IX     = IXT + 1
*
      ELSE
         IX     = IXT
      END IF

      YINEW  = YT

      IF (LGCL) THEN
         IYBOT  = MIN(YIOLD+1.,XNY)
         IYTOP  = MIN(YINEW+1.,XNY)
*
      ELSE
         IYTOP  = MIN(YIOLD+1.,XNY)
         IYBOT  = MIN(YINEW+1.,XNY)
      END IF

      DO 30 IY = IYBOT,IYTOP
         NUM    = NUM + 1
         LIND(NUM,1) = IX
         LIND(NUM,2) = IY
   30 CONTINUE

      END
C     DEBUG SUBCHK
C     END DEBUG
c----------------------------------------------------------------------
c ibm-version of gr3cbc (must be replaced on cray: gr3cbc cft77 )
      SUBROUTINE GR3CBC(NUM,NR,LIND,INDEXF,IFPT,IFLIST,IFREE)

c cbc-collect buckets
c     the entries in the multibucket sort corresponding to buckets
c     having indicies stored in lind are collected and returned in
c     iflist.  the last face found has index num.  the last entry in
c     lind corresponding to a pair of indicies of a bucket is also num.

c     .. scalar arguments ..
      INTEGER           NR,NUM
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      INTEGER         IFLIST(MNF),IFPT(NR*NR+1),INDEXF(MNF*4),LIND(43,2)
      LOGICAL           IFREE(MNF)
c     ..
c     .. local scalars ..
      INTEGER           II,IND,IX,IY,JJ,NUMF
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      NUMF   = 0

      DO 20 II = 1,NUM
         IX     = LIND(II,1)
         IY     = LIND(II,2)

         IF( NR*(IY-1)+IX .EQ. 0 ) GOTO 20
         DO 10 JJ = IFPT(NR*(IY-1)+IX),IFPT(NR*(IY-1)+IX+1) - 1
            IND    = INDEXF(JJ)

            IF (IFREE(IND)) THEN
               IFREE(IND) = .FALSE.
               NUMF   = NUMF + 1
               IFLIST(NUMF) = IND
            END IF

   10    CONTINUE
   20 CONTINUE

      NUM    = NUMF

      DO 30 II = 1,NUMF
         IFREE(IFLIST(II)) = .TRUE.
   30 CONTINUE

      END
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3FXL(RMAX,RMIN,PAVAIL,AVAIL,PSTART,PRM,RM,FLAG)
*********************************************************************
* the upper and lower bounds of visible intervals are stored in rm. *
* rm(1,*) is the left endpoint and rm(2,*) is the right endpoint.   *
* pstart points at the visible interval with the smallest left end  *
* point.  prm(*) is the index of the visible interval to the right  *
* of the interval stored in position rm(.,*) in rm.  the array      *
* avail indicates the indexes of available locations in rm.  upon   *
* initialization each index of rm should be in avail and pavail     *
* should point at the greatest index of avail.  the indexes of rm   *
* in avail need not be ordered.  pavail is the                      *
* index of the cell in avail with greatest index which contains an  *
* index of an available cell in rm.  upon entry                     *
* rmin and rmax are the endpoints of an invisible interval on the   *
* line segment under consideration.  if pstart is zero on exit,     *
* the segment under consideration is invisible. if flag is true,    *
* there is no more free space in rm.                                *
*********************************************************************




c     .. parameters ..
      INTEGER           NDI
      PARAMETER         (NDI=2*65)
c     ..
c     .. scalar arguments ..
      REAL              RMAX,RMIN
      INTEGER           PAVAIL,PSTART
      LOGICAL           FLAG
c     ..
c     .. array arguments ..
      REAL              RM(2,NDI)
      INTEGER           AVAIL(NDI),PRM(NDI)
c     ..
c     .. local scalars ..
      REAL              R1,R2
      INTEGER           INDEX,LAST,NEXT
      LOGICAL           LGCL1,LGCL2,LGCL3
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      LAST   = 0

   10 CONTINUE

      NEXT   = PSTART

   20 CONTINUE

      R1     = RM(1,NEXT)
      R2     = RM(2,NEXT)

      LGCL3  = R1 .GE. RMAX

      IF (.NOT. (LGCL3.OR.R2.LE.RMIN)) THEN

         LGCL1  = R1 .GE. RMIN
         LGCL2  = R2 .LE. RMAX

         IF (LGCL1 .AND. LGCL2) THEN
            PAVAIL = PAVAIL + 1
            AVAIL(PAVAIL) = NEXT
            IF (LAST.NE.0) THEN
               PRM(LAST) = PRM(NEXT)
*
            ELSE IF (PRM(NEXT).EQ.0) THEN
               PSTART = 0
               GO TO 30
*
            ELSE
               PSTART = PRM(NEXT)
               GO TO 10
*
            END IF
*
            NEXT   = LAST
*
         ELSE IF (.NOT. (LGCL1.OR.LGCL2)) THEN
            INDEX  = AVAIL(PAVAIL)
            PRM(INDEX) = PRM(NEXT)
            PRM(NEXT) = INDEX
            LAST   = NEXT
            PAVAIL = PAVAIL - 1

            IF (PAVAIL.EQ.0) THEN
               FLAG   = .TRUE.
               GO TO 30
*
            END IF

            RM(1,INDEX) = RMAX
            RM(2,INDEX) = R2
            RM(2,NEXT) = RMIN
            NEXT   = INDEX
            LGCL3  = .TRUE.
*
         ELSE IF (.NOT.LGCL1) THEN
            RM(2,NEXT) = RMIN
*
         ELSE
            RM(1,NEXT) = RMAX
         END IF
*
      END IF

      IF (.NOT.LGCL3) THEN
         LAST   = NEXT
         NEXT   = PRM(NEXT)
         IF (NEXT.NE.0) GO TO 20
      END IF

   30 CONTINUE
      END
C     DEBUG SUBCHK
C     END DEBUG
c-----------------------------------------------------------------------
c draw an arrow (pfeil) with x,y- direction from point xu,yu
c-----------------------------------------------------------------------
      SUBROUTINE GR3PFL(XU,YU,X,Y,ART)

c     .. scalar arguments ..
      REAL              X,XU,Y,YU
      CHARACTER         ART
c     ..
c     .. local scalars ..
      REAL              AM,XN1,XN2,YN1,YN2
      INTEGER           I,J
c     ..
c     .. local arrays ..
      REAL              XX(2,7),XY(2,7),XZ(2,7),YX(2,7),YY(2,7),YZ(2,7)
c     ..
c     .. external subroutines ..
      EXTERNAL          GRDRW,GRJMP
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. data statements ..

      DATA              XX/.0,1.2,1.275,1.5,1.725,1.275,1.725,1.5,1.8,
     $                  3.,3.,2.85,3.,2.85/
      DATA              YX/.0,0.0,0.225,0.0,0.225,-.225,-.225,0.0,0.0,
     $                  0.,0.,0.15,0.,-.15/
      DATA              XY/.0,1.2,1.275,1.5,1.500,1.500,1.725,1.5,1.8,
     $                  3.,3.,2.85,3.,2.85/
      DATA              YY/.0,0.0,0.225,0.0,-.225,0.000,0.225,0.0,0.0,
     $                  0.,0.,0.15,0.,-.15/
      DATA              XZ/.0,1.2,1.275,1.725,1.725,1.275,1.275,1.725,
     $                  1.8,3.,3.,2.85,3.,2.85/
      DATA              YZ/.0,0.0,0.225,0.225,0.225,-.225,-.225,-.225,
     $                  0.0,0.,0.,0.15,0.,-.15/
c     ..
c     .. executable statements ..
      IF (ART.EQ.'X') THEN
         DO 10 I = 1,7
            XN1    = XX(1,I)*X - YX(1,I)*Y
            YN1    = XX(1,I)*Y + YX(1,I)*X
            XN2    = XX(2,I)*X - YX(2,I)*Y
            YN2    = XX(2,I)*Y + YX(2,I)*X
            CALL GRJMP(XU+XN1,YU+YN1)
            CALL GRDRW(XU+XN2,YU+YN2)
   10    CONTINUE
*
      ELSE IF (ART.EQ.'Y') THEN
         DO 30 I = 1,7
            IF (I.EQ.2 .AND. XN2.LT.0.) THEN
               DO 20 J = 2,4
                  YY(1,J) = -YY(1,J)
   20          CONTINUE
            END IF
*
            XN1    = XY(1,I)*X - YY(1,I)*Y
            YN1    = XY(1,I)*Y + YY(1,I)*X
            XN2    = XY(2,I)*X - YY(2,I)*Y
            YN2    = XY(2,I)*Y + YY(2,I)*X
            CALL GRJMP(XU+XN1,YU+YN1)
            CALL GRDRW(XU+XN2,YU+YN2)
   30    CONTINUE
         YY(1,2) = ABS(YY(1,2))
         YY(1,3) = -ABS(YY(1,3))
         YY(1,4) = ABS(YY(1,4))
*
      ELSE IF (ART.EQ.'Z') THEN
         AM     = XZ(2,2)
         DO 40 I = 1,7
            IF (I.EQ.2 .AND. ABS(XN2).LT.ABS(YN2)) THEN
               XZ(2,2) = XZ(1,2)
               YZ(2,2) = -YZ(2,2)
               YZ(1,3) = -YZ(1,3)
               YZ(2,3) = -YZ(2,3)
               XZ(1,4) = XZ(2,4)
               YZ(1,4) = -YZ(1,4)
            END IF
*
            XN1    = XZ(1,I)*X - YZ(1,I)*Y
            YN1    = XZ(1,I)*Y + YZ(1,I)*X
            XN2    = XZ(2,I)*X - YZ(2,I)*Y
            YN2    = XZ(2,I)*Y + YZ(2,I)*X
            CALL GRJMP(XU+XN1,YU+YN1)
            CALL GRDRW(XU+XN2,YU+YN2)
   40    CONTINUE
         XZ(2,2) = AM
         YZ(2,2) = ABS(YZ(2,2))
         YZ(1,3) = ABS(YZ(1,3))
         YZ(2,3) = -ABS(YZ(2,3))
         XZ(1,4) = XZ(2,1)
         YZ(1,4) = -ABS(YZ(1,4))
*
      ELSE
         WRITE (*,FMT=*) 'GR3PFL: ACHSE NICHT X,Y ODER Z'
      END IF
*
      END
