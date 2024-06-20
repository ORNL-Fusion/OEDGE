C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C Beim Microsoft-Compiler muss das Programm vor den Zeilen C====== ...
C in wenigstens drei Quellprograamme zerlegt werden, weil dieser
C Compiler zu streng  die Uebereinstimmung von Typen prueft, wenn
C alles im gleichen Quellprogrammdeck steht.
c-----------------------------------------------------------------------
C    Achtung  - M.Busch 10. 4. 91
C    Das Programm wurde mit TOOLPACK bearbeitet
C    TOOLPACK erkennt NUMARG als EXTERNAL - NUMARG ist bei CRAY aber
C    eine INTRINSIC Function
C    Bei Portierung auf andere Rechner evtl. NUMARG rausnehmen
C    und dafuer soregen, dass GR3AXS immer mit 8 Argumenten
C    aufgerufen wird
c-----------------------------------------------------------------------
c -- diese version ist fuer cms auf einer ibm --------------------------
c -- lokalisiere jedes vorkommen von 'ibm','cray','cms','mvs', bevor
c -- du dies fuer die cray oder mvs vorbereitest! -------------------
c -- hier wird nirgends mehr direkt auf gr-software bezug genommen
c -- uebergang zu phigs soll erleichtert werden.
c-----------------------------------------------------------------------
c update 07.05.90: trennung von gr-grundsoftware-programmen: gr3plo
c update 01.08.90: kurvenscharen bei gr3net und gr3nt1
c update  6.11.90: m. busch save
c update 14.03.91: G. GRoten ; bei DRAHT immer LNFC=0
c update 14.03.91: G. GRoten ; bei JCO-4000 Strichstaerke 14
C update 22.07.91: G.Groten  ; COMMON /GRANTN/; NUA=0 annotation Text
C update 19.08.91: M. Busch  ; NUMARG rausgenommen
C update 29.10.91: G.Groten  ; KEN(ILA),PS(3,ILA) => KEN(IPA),PS(3,IPA)
C UPDATE 31.10.91: G.GROTEN  ; FEHLER IN GR3TRI MIT IP0-IPOLD KORRIG.
C Update 07.11.91: G.Groten  , ISMIN,ISMAX -> GRISMI,GRISMA
C Update 20.01.91: G.Groten  , EXTERNAL SUB wegen PC entfernt. === hinzu
C Update 13.10.92: G.Groten  , in GR3QS1 und GR3QS2 Stack vergroessert.
c----------------------------- mit OPT(3),, nogo
c aufgerufene unterprogramme:
c
c    grsclc, grsclv, grdsh, grnwpn, grchrc, grtxt, grftoc, grspts,
c    grln, grnxtf, scan
c
c-----------------------------------------------------------------------
      SUBROUTINE GR3DIM(LWORK,IER)
c
c     nfl ist die anzahl der flaechen in einem bild.
c
c     lwork ist die laenge des hilfspeichers work. es wird, um platz
c     zu sparen, vorausgesetzt, dass real, integer und logical den
c     gleichen speicher belegen duerfen.
c
c     ier als eingabewert bewirkt,
c        wenn 1, dass eine fehlernachricht geschrieben wird,
c                ier den wert eines fehlercode annimmt,
c                und ins aufrufende programm zurueckverzweigt wird;
c        wenn 2, dass ier den wert eines fehlercode annimmt,
c                und ins aufrufende programm zurueckverzweigt wird;
c        sonst , dass eine fehlernachricht geschrieben wird,
c                dass das programm dann mit stop abbricht.
c-----------------------------------------------------------------------







c     .. parameters ..
      INTEGER           MFL
      PARAMETER         (MFL=256)
c     ..
c     .. scalar arguments ..
      REAL              D1,D10,D2,D3,D4,D5,THE1,THE2,THE3
      INTEGER           IER,IFACES,IX,IY,JCO,KINDAX,LWORK,N1,N2,NID,
     $                  NIDIM,NT,NU,NV,NX,NY,I6,I7,I8,I9
      LOGICAL           CHORI
      CHARACTER (len=*) AX1,AX2,AX3
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IFDREH,IFL,IL,IL0,
     $                  IP0,ISY,ISZ,IW,KINDA,KOAX,L,LW,MAFL,MNF,MNL,MNP,
     $                  NOTDIM
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,L,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              A(46*MNP),EX33(3,3),PT(3,NID,N2),TAB(NT,NY),
     $                  TV(3),WERBER(3,2),X(NIDIM,NV),XX(NT),
     $                  Y(NIDIM,NV),YY(NY),Z(NIDIM,NV)
      INTEGER           IGR(2,MFL)
      CHARACTER(len=20) CHAXS(3)
c     ..
c     .. subroutine arguments ..
c     ..
c     .. local scalars ..
      DOUBLE PRECISION  SUX,SUY,SUZ
      INTEGER           I,J,JCOL,K,MIS,NI,NJ
      LOGICAL           LCOL
c     ..
c     .. local arrays ..
      REAL              ROT(3,3),ROT1(3,3),ROT2(3,3),ROT3(3,3)
c     ..
c     .. external functions ..
      INTEGER           GRISMA,GRISMI
      EXTERNAL          GRISMA,GRISMI
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3CN1,GR3CN2,GR3CON,GR3CUB,GR3INV,GR3RTX,GR3TRF
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MAX,MIN,MOD
c     ..
c     .. common blocks ..

C---- COMMON-Bloecke zu GR3ANT
      INTEGER            NUMAAN, LEMAAN
      PARAMETER          (NUMAAN=32,LEMAAN=32)
      INTEGER            NUA,INDA,LENA,IFOA,ICOA,IBOA
      REAL               SIZA,DS1A,DS2A,DS3A,ANGA,PTXA,PTYA
      CHARACTER (len=LEMAAN) CHAA
      COMMON /GRANTN/ NUA,INDA(NUMAAN),LENA(NUMAAN),SIZA(NUMAAN),
     $       IFOA(NUMAAN),DS1A(NUMAAN),DS2A(NUMAAN),DS3A(NUMAAN),
     $       ANGA(NUMAAN),PTXA(NUMAAN),PTYA(NUMAAN),ICOA(NUMAAN),
     $       IBOA(NUMAAN)
CDEC$ PSECT /GRANTN/ NOSHR
      COMMON /GRANTC/ CHAA(NUMAAN)
CDEC$ PSECT /GRANTC/ NOSHR

      COMMON            /GR3BAC/ICOLBA,CENTR,WIND
CDEC$ PSECT /GR3BAC/ NOSHR
      COMMON            /GR3DRE/IFDREH,XDREH,YDREH,ZDREH,KINDA,KOAX
CDEC$ PSECT /GR3DRE/ NOSHR
      COMMON            /GR3SAC/CAXS,CORI
CDEC$ PSECT /GR3SAC/ NOSHR
      COMMON            /GR4COM/IFL,JGR,XC,YC,COMA
CDEC$ PSECT /GR4COM/ NOSHR
      REAL              CENTR,XDREH,YDREH,ZDREH
      LOGICAL           WIN,ZENPRO
      REAL              COMA(3,2),GR(3,3),RT(3,3),WIND(4),XC(500),
     $                  YC(500)
      INTEGER           ICOLBA(8),JGR(7,MFL)
      CHARACTER(len=20) CAXS(3)
      CHARACTER(len=1)  CORI
c     ..
c     .. equivalences ..
      EQUIVALENCE       (LCOL,JCOL)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/,/GRANTN/,/GRANTC/,/GR3BAC/,/GR3DRE/,
     $     /GR3SAC/,/GR4COM/
      SAVE
c     ..
c     .. executable statements ..

c-----------------------------------------------------------------------
      NUA = 0
      KOAX   = 1
      ICOLBA(1) = 1
      ICOLBA(2) = 2
      ICOLBA(3) = 3
      ICOLBA(4) = 4
      ICOLBA(5) = 5
      ICOLBA(6) = 6
      ICOLBA(7) = 7
      ICOLBA(8) = 8
      KINDA  = 4
      DO 20 I = 1,3
         DO 10 J = 1,3
            RT(I,J) = 0.
   10    CONTINUE
         COMA(I,1) = 0.
         COMA(I,2) = 0.
         RT(I,I) = 1.
   20 CONTINUE
      NOTDIM = 1234567890
      WIN    = .FALSE.
      MNP    = LWORK/46
      MNL    = 2*MNP
      MNF    = MNP
      IFL    = 0
      MAFL   = 0
      IL0    = 0
      IF0    = 0
      IP0    = 0
      ISY    = 1 + MNP
      ISZ    = ISY + MNP
      IF1    = ISZ + MNP
      IF2    = IF1 + MNF
      IF3    = IF2 + MNF
      IF4    = IF3 + MNF
      IL     = IF4 + MNF
      L      = IL + MNL + MNL
      IW     = L + MNL + MNL
      LW     = LWORK - IW + 1
      ZENPRO = .FALSE.
      IER    = 0
      GO TO 230

c----------------------------------------------------------------------
      ENTRY             GR3PEN(IGR,IER)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3PEN RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      DO 30 I = 1,MAFL
         JGR(3,I) = IGR(1,I)
         JGR(5,I) = IGR(2,I)
   30 CONTINUE
      IER    = 0
      GO TO 230

c----------------------------------------------------------------------
      ENTRY             GR3EXT(A,IER,EX33)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3PEN RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      EX33(1,1) = A(GRISMI(IP0,A,1))
      EX33(1,3) = A(GRISMA(IP0,A,1))
      EX33(2,1) = A(ISY-1+GRISMI(IP0,A(ISY),1))
      EX33(2,3) = A(ISY-1+GRISMA(IP0,A(ISY),1))
      EX33(3,1) = A(ISZ-1+GRISMI(IP0,A(ISZ),1))
      EX33(3,3) = A(ISZ-1+GRISMA(IP0,A(ISZ),1))
      SUX    = 0.
      SUY    = 0.
      SUZ    = 0.
      DO 40 I = 1,IP0
         SUX    = SUX + A(I)
         SUY    = SUY + A(ISY-1+I)
         SUZ    = SUZ + A(ISZ-1+I)
   40 CONTINUE
      IF (IP0.NE.0) THEN
         EX33(1,2) = SUX/IP0
         EX33(2,2) = SUY/IP0
         EX33(3,2) = SUZ/IP0
      END IF
*
      GO TO 230

      ENTRY             GR3NET(A,IER,NID,PT,N1,IX,N2,IY,IFACES,JCO)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NET RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IFL    = IFL + 1
      MAFL   = IFL

      IF (IFL.GT.MFL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*) 'GR3NET RC=20: TOO MANY SURFACES'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF

      NI     = N1
      JGR(1,IFL) = NI

      IF ((NI-1)*IX+1.GT.NID) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NET RC=30: DIMENSION NID MUST BE > OR = N1'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 30
         GO TO 230
*
      END IF

      NJ     = N2
      JGR(2,IFL) = NJ
      JGR(3,IFL) = MIN(MAX(MOD(JCO+4000,2000),1),8)
      JGR(5,IFL) = MIN(MAX((JCO+2000)/2000,0),8)*2 +
     $             16*MIN(MAX(JCO,0),1)
      IF (JCO.LT.-2000) JGR(5,IFL) = 14

      IF (IP0+NI*NJ.GT.MNP .OR. IF0+ (NI-1)* (NJ-1).GT.MNF .OR.
     $    IL0+ (NI-1)*NJ+NI* (NJ-1).GT.MNL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NET RC=40: Workspace too small'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 40
         GO TO 230
*
      END IF

c---- call gr3con(pt,nid,ni,nj,x,  y   ,  z   ,iface1,iface2,iface3,
c---->            iface4,iline,lnfc,il0,if0,ip0,ifaces)
      CALL GR3CON(PT,NID,NI,IX,NJ,IY,A,A(ISY),A(ISZ),A(IF1),A(IF2),
     $            A(IF3),A(IF4),A(IL),A(L),IL0,IF0,IP0,IFACES)
      JGR(6,IFL) = IP0
      JGR(7,IFL) = IF0
      JGR(4,IFL) = IL0
      IER    = 0
      GO TO 230

c----------------------------------------------------------------------

      ENTRY             GR3NT1(A,IER,NIDIM,X,Y,Z,NU,IX,NV,IY,IFACES,JCO)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NT1 RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IFL    = IFL + 1
      MAFL   = IFL

      IF (IFL.GT.MFL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*) 'GR3NT1 RC=20: TOO MANY SURFACES'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF

      NI     = NU

      IF ((NI-1)*IX+1.GT.NIDIM) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NT1 RC=30: DIMENSION NID MUST BE > OR = N1'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 30
         GO TO 230
*
      END IF
*
      JGR(1,IFL) = NI
      NJ     = NV
      JGR(2,IFL) = NJ
      JGR(3,IFL) = MIN(MAX(MOD(JCO+4000,2000),1),8)
      JGR(5,IFL) = MIN(MAX((JCO+2000)/2000,0),8)*2 +
     $             16*MIN(MAX(JCO,0),1)
      IF (JCO.LT.-2000) JGR(5,IFL) = 14

      IF (IP0+NI*NJ.GT.MNP .OR. IF0+ (NI-1)* (NJ-1).GT.MNF .OR.
     $    IL0+ (NI-1)*NJ+NI* (NJ-1).GT.MNL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NT1 RC=40: Workspace too small'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 40
         GO TO 230
*
      END IF

      CALL GR3CN1(X,Y,Z,NIDIM,NI,IX,NJ,IY,A,A(ISY),A(ISZ),A(IF1),A(IF2),
     $            A(IF3),A(IF4),A(IL),A(L),IL0,IF0,IP0,IFACES)
      JGR(6,IFL) = IP0
      JGR(7,IFL) = IF0
      JGR(4,IFL) = IL0
      IER    = 0
      GO TO 230

c----------------------------------------------------------------------

      ENTRY             GR3SUB(A,IER,SUB,JCO,D1,D2,D3,D4,D5,I6,I7,I8,I9,
     $                  D10)
*
      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3SUB RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IFL    = IFL + 1
      MAFL   = IFL
      IF (IFL.GT.MFL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*) 'GR3SUB RC=20: TOO MANY SURFACES'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF
*
      JGR(3,IFL) = MIN(MAX(MOD(JCO+4000,2000),1),8)
      JGR(5,IFL) = MIN(MAX((JCO+2000)/2000,0),8)*2 +
     $             16*MIN(MAX(JCO,0),1)
      IF (JCO.LT.-2000) JGR(5,IFL) = 14
c     call sub(x,  y   ,  z   ,iface1,iface2,iface3,iface4,iline,lnfc,
c    >         indl,indf,np,d1,d2,d3,d4,d5,I6,I7,I8,I9,d10)
      CALL SUB(A,A(ISY),A(ISZ),A(IF1),A(IF2),A(IF3),A(IF4),A(IL),A(L),
     $         IL0,IF0,IP0,D1,D2,D3,D4,D5,I6,I7,I8,I9,D10)
      JGR(1,IFL) = 0
      JGR(2,IFL) = 0
      JGR(6,IFL) = IP0
      JGR(7,IFL) = IF0
      JGR(4,IFL) = IL0
      IER    = 0

      GO TO 230

c----------------------------------------------------------------------

      ENTRY             GR3FUN(A,IER,NT,TAB,NX,IX,XX,NY,IY,YY,IFACES,
     $                  JCO)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3FUN RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IFL    = IFL + 1

      IF (IFL.GT.MFL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*) 'GR3FUN RC=20: TOO MANY SURFACES'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF
*
      MAFL   = IFL

      NI     = NX

      IF ((NI-1)*IX+1.GT.NT) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3FUN RC=30: DIMENSION NT MUST BE > OR = NX'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 30
         GO TO 230
*
      END IF
*
      JGR(1,IFL) = NI

      NJ     = NY
      JGR(2,IFL) = NJ

      IF (IP0+NI*NJ.GT.MNP .OR. IF0+ (NI-1)* (NJ-1).GT.MNF .OR.
     $    IL0+ (NI-1)*NJ+NI* (NJ-1).GT.MNL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3FUN RC=40: Workspace too small'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 40
         GO TO 230
*
      END IF

      JGR(3,IFL) = MIN(MAX(MOD(JCO+4000,2000),1),8)
      JGR(5,IFL) = MIN(MAX((JCO+2000)/2000,0),8)*2 +
     $             16*MIN(MAX(JCO,0),1)
      IF (JCO.LT.-2000) JGR(5,IFL) = 14
      CALL GR3CN2(TAB,NT,XX,YY,NI,IX,NJ,IY,A,A(ISY),A(ISZ),A(IF1),
     $            A(IF2),A(IF3),A(IF4),A(IL),A(L),IL0,IF0,IP0,IFACES)
      JGR(6,IFL) = IP0
      JGR(7,IFL) = IF0
      JGR(4,IFL) = IL0
      IER    = 0
      GO TO 230

c----------------------------------------------------------------------
c---  die argumente chori und kindax sind spaeter hinzugekommen
c---  darum die mimik mit numarg !
      ENTRY             GR3AXS(A,IER,EX33,WERBER,CHAXS,CHORI,KINDAX,JCO)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3AXS RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IFL    = IFL + 1
      IF (IFL.GT.MFL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*) 'GR3AXS RC=20: TOO MANY SURFACES'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF
*
      IF (IP0+8.GT.MNP .OR. IL0+12.GT.MNL) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3NET RC=30: Workspace too small'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 30
         GO TO 230
*
      END IF
*
      MAFL   = IFL
      JGR(1,IFL) = -32767
      JGR(2,IFL) = 1
C     IF (NUMARG().EQ.6) THEN
c------- lcol ist equivalent jcol
C        LCOL   = CHORI
C        CORI   = '0'
*
C     ELSE
         JCOL   = JCO
         CORI   = '0'
         IF ( CHORI) CORI='1'
         KINDA  = MOD(MAX(ABS(KINDAX)-1,0),6) + 1
C     END IF
*
      JGR(3,IFL) = MIN(MAX(MOD(JCO+4000,2000),1),8)
      JGR(5,IFL) = MIN(MAX((JCO+2000)/2000,0),8)*2 +
     $             16*MIN(MAX(JCO,0),1)

      CALL GR3CUB(A,A(ISY),A(ISZ),A(IL),A(L),IL0,IP0,EX33)
      JGR(6,IFL) = IP0
      JGR(7,IFL) = IF0
      JGR(4,IFL) = IL0

      DO 50 I = 1,3
         CAXS(I) = CHAXS(I)
         COMA(I,1) = WERBER(I,1)
         COMA(I,2) = WERBER(I,2)
         IF (EX33(I,2).EQ.EX33(I,1)) COMA(I,2) = WERBER(I,1)
   50 CONTINUE

      KOAX   = 0
      IER    = 0
      GO TO 230

c----------------------------------------------------------------------

      ENTRY             GR3KIL(IER)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3KIL RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IF (IFL.LE.0) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3KIL RC=20: THERE IS NO SURFACE TO BE KILLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF

      IFL    = IFL - 1
      IF (IFL.GT.0) THEN
         IP0    = JGR(6,IFL)
         IF0    = JGR(7,IFL)
         IL0    = JGR(4,IFL)
*
      ELSE
         IP0    = 0
         IF0    = 0
         IL0    = 0
      END IF

      IER    = 0
      GO TO 230

c----------------------------------------------------------------------

      ENTRY             GR3ROT(A,IER,AX1,THE1,AX2,THE2,AX3,THE3)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ROT RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      MIS    = 0
      IF (AX1.EQ.'Z' .OR. AX1.EQ.'Z ') THEN
*
      ELSE IF (AX1.EQ.'Y' .OR. AX1.EQ.'Y ') THEN
         IF (THE1.NE.0.) THEN
            IF (ZENPRO) GO TO 190
         END IF
*
      ELSE IF (AX1.EQ.'X' .OR. AX1.EQ.'X ') THEN
         IF (THE1.NE.0.) THEN
            IF (ZENPRO) GO TO 190
         END IF
*
      ELSE IF (AX1.EQ.'XU' .OR. AX1.EQ.'YU' .OR. AX1.EQ.'ZU') THEN
         IF (ZENPRO) GO TO 190
         MIS    = MIS + 1
*
      ELSE
         GO TO 180
*
      END IF

      CALL GR3RTX(ROT1,THE1,AX1)

      IF (AX2.EQ.'Z' .OR. AX2.EQ.'Z ') THEN
*
      ELSE IF (AX2.EQ.'Y' .OR. AX2.EQ.'Y ') THEN
         IF (THE2.NE.0.) THEN
            IF (ZENPRO) GO TO 190
         END IF
*
      ELSE IF (AX2.EQ.'X' .OR. AX2.EQ.'X ') THEN
         IF (THE2.NE.0.) THEN
            IF (ZENPRO) GO TO 190
         END IF
*
      ELSE IF (AX2.EQ.'XU' .OR. AX2.EQ.'YU' .OR. AX2.EQ.'ZU') THEN
         IF (ZENPRO) GO TO 190
         MIS    = MIS + 1
*
      ELSE
         GO TO 180
*
      END IF

      CALL GR3RTX(ROT2,THE2,AX2)

      IF (AX3.EQ.'Z' .OR. AX3.EQ.'Z ') THEN
*
      ELSE IF (AX3.EQ.'Y' .OR. AX3.EQ.'Y ') THEN
         IF (THE3.NE.0.) THEN
            IF (ZENPRO) GO TO 190
         END IF
*
      ELSE IF (AX3.EQ.'X' .OR. AX3.EQ.'X ') THEN
         IF (THE3.NE.0.) THEN
            IF (ZENPRO) GO TO 190
         END IF
*
      ELSE IF (AX3.EQ.'XU' .OR. AX3.EQ.'YU' .OR. AX3.EQ.'ZU') THEN
         IF (ZENPRO) GO TO 190
         MIS    = MIS + 1
*
      ELSE
         GO TO 180
*
      END IF

      CALL GR3RTX(ROT3,THE3,AX3)

      IF (MIS.NE.0 .AND. MIS.NE.3) GO TO 180

      DO 80 I = 1,3
         DO 60 K = 1,3
            ROT(I,K) = ROT2(I,1)*ROT1(1,K) + ROT2(I,2)*ROT1(2,K) +
     $                 ROT2(I,3)*ROT1(3,K)
   60    CONTINUE
         DO 70 K = 1,3
            ROT2(I,K) = ROT(I,K)
   70    CONTINUE
   80 CONTINUE

      DO 100 I = 1,3
         DO 90 K = 1,3
            ROT(I,K) = ROT3(I,1)*ROT2(1,K) + ROT3(I,2)*ROT2(2,K) +
     $                 ROT3(I,3)*ROT2(3,K)
   90    CONTINUE
  100 CONTINUE
c---------------------- drehung um ursprungachsen ? -------------------
      IF (MIS.EQ.3) THEN
         CALL GR3INV(GR)
         DO 120 I = 1,3
            DO 110 K = 1,3
               ROT2(I,K) = ROT(I,1)*GR(1,K) + ROT(I,2)*GR(2,K) +
     $                     ROT(I,3)*GR(3,K)
  110       CONTINUE
  120    CONTINUE
         DO 140 I = 1,3
            DO 130 K = 1,3
               ROT(I,K) = RT(I,1)*ROT2(1,K) + RT(I,2)*ROT2(2,K) +
     $                    RT(I,3)*ROT2(3,K)
  130       CONTINUE
  140    CONTINUE
      END IF
c--------------------- rt=gr gesamttransformierte einheitsmatrix -------
      DO 160 I = 1,3
         DO 150 K = 1,3
            GR(I,K) = ROT(I,1)*RT(1,K) + ROT(I,2)*RT(2,K) +
     $                ROT(I,3)*RT(3,K)
  150    CONTINUE
  160 CONTINUE

      DO 170 I = 1,3
         RT(I,1) = GR(I,1)
         RT(I,2) = GR(I,2)
         RT(I,3) = GR(I,3)
  170 CONTINUE

c     call gr3trf(rot,x,  y   ,  z   ,ip0,inc)
      CALL GR3TRF(ROT,A,A(ISY),A(ISZ),IP0,1)
      IER    = 0
      GO TO 230

  180 CONTINUE
      IF (IER.NE.2) WRITE (*,FMT=*)
     $    'GR3ROT RC=20: AXIS NOT ''X'' OR ''Y'' OR ''Z'''
      IF (IER.NE.2 .AND. IER.NE.1) STOP
      IER    = 20
      GO TO 230

  190 CONTINUE
      IF (IER.NE.2) WRITE (*,FMT=*)
     $    'GR3ROT RC=30: AFTER A CENTRAL PROJECTION ONLY ROTATION ''Z'''
      IF (IER.NE.2 .AND. IER.NE.1) STOP
      IER    = 30
      GO TO 230

c----------------------------------------------------------------------

      ENTRY             GR3ORI(A,IER)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ORI RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      IF (ZENPRO) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ORI RC=20: NOT POSSIBLE AFTER GR3CEN!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 230
*
      END IF

      CALL GR3INV(RT)

c     call gr3trf(rt,x,  y   ,  z   ,ip0,inc)
      CALL GR3TRF(RT,A,A(ISY),A(ISZ),IP0,1)

      DO 210 I = 1,3
         DO 200 J = 1,3
            GR(I,J) = 0.
            RT(I,J) = 0.
  200    CONTINUE
         GR(I,I) = 1.
         RT(I,I) = 1.
  210 CONTINUE
      IER    = 0
      GO TO 230
c-----------------------------------------------------------------------
      ENTRY             GR3TRL(A,IER,TV)

      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3TRL RC=10: GR3DIM MUST FIRST BE CALLED!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 230
*
      END IF

      DO 220 I = 1,IP0
         A(I)   = A(I) + TV(1)
         A(ISY-1+I) = A(ISY-1+I) + TV(2)
         A(ISZ-1+I) = A(ISZ-1+I) + TV(3)
  220 CONTINUE
      IER    = 0
  230 CONTINUE
      END
      SUBROUTINE GR3TR1(AR,IER,XYZ,N,RIADJ,RIEND,IFACES,JCO)
c     .. scalar arguments ..
      REAL              XYZ,RIADJ,RIEND
      INTEGER           IER,IFACES,JCO,N
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP)
c     ..
c     .. local scalars ..
      REAL              U,V,W
      INTEGER           I,IPX,IPY,IPZ
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3OBB,GR3OBE,GR3TRA,GR3TRC
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..
   10 CONTINUE
      IF (IP0+N.GT.MNP) THEN
         WRITE (*,FMT=*) 'GR3TR1: zu viele Punkte.'
         STOP
*
      END IF
*
      CALL GR3OBB
      CALL GR3TRA(AR,IER,XYZ,N,RIADJ,RIEND,JCO)
      IPX    = IP0
      IPY    = MNP + IP0
      IPZ    = 2*MNP + IP0
      CALL GR3OBE(AR,IER,IFACES,JCO)
      DO 20 I = 1,N
         CALL GR3TRC(AR(IPX+I),AR(IPY+I),AR(IPZ+I),U,V,W)
         AR(IPX+I) = U
         AR(IPY+I) = V
         AR(IPZ+I) = W
   20 CONTINUE
      END
      SUBROUTINE GR3TRI(AR,IER,XYZ,NG,NPO,NB,LIK,IFACES,JCO)
c     .. scalar arguments ..
      INTEGER           IER,IFACES,JCO,LIK,NG
c     ..
c     .. array arguments ..
      REAL              AR(*),XYZ(3,*)
      INTEGER           NB(NG),NPO(NG)
c     ..
c     .. local scalars ..
      REAL              RV,U,V,W
      INTEGER           I,IE,IN1,IN2,IPX,IPY,IPZ,IS,IS0,IS1,IV,J,JG,K,
     $                  KK,LAST,LASTP,LWK,MAPO
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3OBB,GR3OBE,GR3TRA,GR3TRC,GRDELN,GREDGE,GRMESH
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. equivalences ..
      EQUIVALENCE       (RV,IV)
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      MAPO   = 0
      DO 10 I = 1,NG
         MAPO   = MAPO + NPO(I)
   10 CONTINUE
      IF (IP0+MAPO.GT.MNP) THEN
         WRITE (*,FMT=*) 'GR3TRI: zu viele Punkte.'
         STOP
*
      END IF
*
      CALL GR3OBB
      JG     = 1
      DO 60 J = 1,NG
c------- triangulation ueber xy
         IS     = 46*MNP - 7*NPO(J) + 10
         IS1    = 46*MNP - NPO(J) + 1
         IS0    = IS - 2* (NPO(J)/4)
         IF (IS0.LT.34.375*MNP) THEN
            WRITE (*,FMT=*) 'GR3TRI: Arbeitsspeicher zu klein'
            STOP
*
         END IF
*
         CALL GRMESH(NPO(J),XYZ(1,JG),AR(IS),AR(IS1),IE)
         IF (IE.NE.0) THEN
            IF (IE.EQ.1) THEN
               WRITE (*,FMT=*) 'GR3TRI,GRMESH: weniger als 3 Punkte'
*
            ELSE
               WRITE (*,FMT=*) 'GR3TRI,GRMESH: alle Punkte kollinear'
            END IF
*
            STOP
*
         END IF
*
         DO 20 I = 0,NB(J) - 1
            IF (I.EQ.0) THEN
               IN1    = NB(J)
               IN2    = 1
*
            ELSE
               IN1    = I
               IN2    = I + 1
            END IF
*
            LWK    = NPO(J)/4
            CALL GREDGE(IN1,IN2,XYZ(1,JG),LWK,AR(IS0),AR(IS),AR(IS1),IE)
            IF (IE.NE.0) THEN
               WRITE (*,FMT=*) 'GR3TRI,GREDGE Fehler Nr.',IE,IN1,IN2
               IF (IE.EQ.3) THEN
                  WRITE (*,FMT=*)
     $              'moeglicher Grund: doppelter X,Y-Punkt'
               END IF
*
               STOP
*
            END IF
*
   20    CONTINUE
         KK     = 0
   30    CONTINUE
         K      = 0
         DO 50 I = NB(J),1,-1
   40       CONTINUE
            RV     = AR(IS1+I-1)
            LASTP  = IV
            RV     = AR(IS+LASTP-1)
            LAST   = IV
            IF (LAST.EQ.0) THEN
               RV     = AR(IS+LASTP-2)
               IF ((I.GT.1.AND.IV.NE.I-1) .OR.
     $             (I.EQ.1.AND.IV.NE.NB(J))) THEN
                  CALL GRDELN(NPO(J),I,IV,AR(IS),AR(IS1),IE)
                  IF (IE.NE.0) THEN
                     IF (IE.EQ.2 .OR. IE.EQ.4) GO TO 50
                     WRITE (*,FMT=*) 'GR3TRI: GRDELN,ier=',IE,I,IV
                     STOP
*
                  END IF
*
                  K      = K + 1
                  GO TO 40
*
               END IF
*
            END IF
*
   50    CONTINUE
         KK     = KK + 1
         IF (KK.LE.9 .AND. K.GT.0) GO TO 30
         CALL GR3TRA(AR,IER,XYZ(1,JG),NPO(J),AR(IS),AR(IS1),LIK)
         JG     = JG + NPO(J)
   60 CONTINUE
      IPX    = IP0
      IPY    = MNP + IP0
      IPZ    = 2*MNP + IP0
      IPOLD  = IP0
      CALL GR3OBE(AR,IER,IFACES,JCO)
C---- IP0 ist jetzt geaendert !
      DO 70 I = 1,IP0-IPOLD
         CALL GR3TRC(AR(IPX+I),AR(IPY+I),AR(IPZ+I),U,V,W)
         AR(IPX+I) = U
         AR(IPY+I) = V
         AR(IPZ+I) = W
   70 CONTINUE
      END
      SUBROUTINE GR3OBB

c     .. scalar arguments ..
      INTEGER          IER,IFACES,JCO,M,N
c     ..
c     .. common blocks ..
      REAL             WIN,ZENPRO
      INTEGER          IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,L,LW,
     $                 MAFL,MNF,MNL,MNP,NOTDIM
      REAL             GR(3,3),RT(3,3)
      COMMON           /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                 IF2,IF3,IF4,IL,L,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL             AR(46*MNP),P(3,M,N)
      INTEGER          LIN(M,N)
c     ..
c     .. local scalars ..
      REAL             DIMFAC
      INTEGER          IIA,ILA,ILMAX,IPA,IS,IS1,IS2,IS3,IS4
c     ..
c     .. external subroutines ..
      EXTERNAL         GR3OBJ,GR3OBT,GR3SUB
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..
c     save is,is1,ipa,is2,iia,ila,is3,is4,ilmax

      IS     = 15*MNP + 1
      DIMFAC = (46-15)/8.
      ILMAX  = MNP*DIMFAC
      IS1    = IS + ILMAX*3
      IS2    = IS1 + ILMAX
      IS3    = IS2 + ILMAX
      IS4    = IS3 + ILMAX
      IPA    = 0
      IIA    = 0
      ILA    = 0
      RETURN

      ENTRY            GR3OBP(AR,IER,P,M,N,LIN)
*
      IF (M.LT.3 .OR. M.GT.4) THEN
         WRITE (*,FMT=*) 'GR3OBP: M only 3 or 4'
         STOP
*
      END IF
*
      IF (IPA+N*4.GT.ILMAX) THEN
         WRITE (*,FMT=*) 'GR3OBP: zu viele Teile, AR groesser machen!'
         STOP
*
      END IF
*
      CALL GR3OBT(AR(IS),AR(IS1),AR(IS2),P,M,N,LIN,IPA,IIA,ILA,ILMAX)
      RETURN

      ENTRY            GR3OBE(AR,IER,IFACES,JCO)
*
      CALL GR3SUB(AR,IER,GR3OBJ,JCO,AR(IS),AR(IS1),AR(IS2),AR(IS3),
     $            AR(IS4),IPA,IIA,ILA,IFACES,D4)
      END
C=======================================================================
      SUBROUTINE GR3CON(PT,NID,NI,IX,NJ,IY,X,Y,Z,IFACE1,IFACE2,IFACE3,
     $                  IFACE4,ILINE,LNFC,IL0,IF0,IP0,IFACES)

c     pt is an array of points which determines a surface which is to
c     be plotted.  pt(k,i,j) is the k-th coordinate of the (i,j)-th
c     point.

***********************************************************************



c     .. parameters ..
      REAL              BLIND1,BLIND2
      PARAMETER         (BLIND1=-75.7501E20,BLIND2=-75.7499E20)
      INTEGER           MFL
      PARAMETER         (MFL=256)
c     ..
c     .. scalar arguments ..
      INTEGER           IF0,IFACES,IL0,IP0,IX,IY,NI,NID,NJ,NT
c     ..
c     .. common blocks ..
      INTEGER           IF1,IF2,IF3,IF4,IFL,IFO,IL,ILO,IPO,ISY,ISZ,IW,
     $                  LQ,LW,MAFL,MNF,MNL,MNP,NOTDIM
      REAL              WIN,ZENPRO
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,ILO,IFO,IPO,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              COMA(3,2),GR(3,3),RT(3,3),XC(500),YC(500)
      INTEGER           JGR(7,MFL)
      COMMON            /GR4COM/IFL,JGR,XC,YC,COMA
CDEC$ PSECT /GR4COM/ NOSHR
c     ..
c     .. array arguments ..
      REAL              PT(3,NID, (NJ-1)*IY+1),T(NT, (NJ-1)*IY+1),
     $                  X(MNP),XT((NI-1)*IX+1),XX(NID, (NJ-1)*IY+1),
     $                  Y(MNP),YT((NJ-1)*IY+1),YY(NID, (NJ-1)*IY+1),
     $                  Z(MNP),ZZ(NID, (NJ-1)*IY+1)
      INTEGER           IFACE1(MNF),IFACE2(MNF),IFACE3(MNF),IFACE4(MNF),
     $                  ILINE(2,MNL),LNFC(2,MNL)
c     ..
c     .. local scalars ..
      REAL              SICHTB,XU,YU,ZU
      INTEGER           I,IAL,II,IND1,IND2,INDF,INDL,INDP0,INF,INP,IP,J,
     $                  JJ,KIND,NSI
      LOGICAL           DRAHT,IN1,IN2,KORR
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3TRC
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MAX,MIN,MOD
c     ..
c     .. save statement ..
      SAVE /GR3SAV/, /GR4COM/
      SAVE
c     ..
c     .. executable statements ..

***********************************************************************

c     indl is the number of faces under consideration.
c     indf is the number of lines under consideration.
c     ip is the number of points under consideration.

      KIND   = 1
      GO TO 10
c-----------------------------------------------------------------------
      ENTRY             GR3CN1(XX,YY,ZZ,NID,NI,IX,NJ,IY,X,Y,Z,IFACE1,
     $                  IFACE2,IFACE3,IFACE4,ILINE,LNFC,IL0,IF0,IP0,
     $                  IFACES)
*
      KIND   = 2
      GO TO 10
c-----------------------------------------------------------------------
      ENTRY             GR3CN2(T,NT,XT,YT,NI,IX,NJ,IY,X,Y,Z,IFACE1,
     $                  IFACE2,IFACE3,IFACE4,ILINE,LNFC,IL0,IF0,IP0,
     $                  IFACES)
*
      KIND   = 3
   10 CONTINUE
      INDL   = IL0
      INDF   = IF0
      IP     = IP0
      DRAHT  = MOD(IFACES,2) .EQ. 0
      IN2    = .TRUE.
      IN1    = IFACES .LT. 4
      IN2    = IFACES .LT. 2 .OR. IFACES .EQ. 4 .OR. IFACES .EQ. 5

      SICHTB = 0.
      NSI    = 0
      DO 60 I = 1,NI - 1
         INDP0  = IP
         IF (KIND.EQ.1) THEN
            DO 20 J = 1,NJ - 1
               IP     = IP + 1
               X(IP)  = PT(1, (I-1)*IX+1, (J-1)*IY+1)
               Y(IP)  = PT(2, (I-1)*IX+1, (J-1)*IY+1)
               Z(IP)  = PT(3, (I-1)*IX+1, (J-1)*IY+1)
               IF (Z(IP).LT.BLIND1 .OR. Z(IP).GT.BLIND2) THEN
                  SICHTB = SICHTB + Z(IP)
                  NSI    = NSI + 1
               END IF
*
   20       CONTINUE
*
         ELSE IF (KIND.EQ.2) THEN
            DO 30 J = 1,NJ - 1
               IP     = IP + 1
               X(IP)  = XX((I-1)*IX+1, (J-1)*IY+1)
               Y(IP)  = YY((I-1)*IX+1, (J-1)*IY+1)
               Z(IP)  = ZZ((I-1)*IX+1, (J-1)*IY+1)
               IF (Z(IP).LT.BLIND1 .OR. Z(IP).GT.BLIND2) THEN
                  SICHTB = SICHTB + Z(IP)
                  NSI    = NSI + 1
               END IF
*
   30       CONTINUE
*
         ELSE
            DO 40 J = 1,NJ - 1
               IP     = IP + 1
               IF (T((I-1)*IX+1, (J-1)*IY+1).LT.BLIND1 .OR.
     $             T((I-1)*IX+1, (J-1)*IY+1).GT.BLIND2) THEN
                  CALL GR3TRC(XT((I-1)*IX+1),YT((J-1)*IY+1),
     $                        T((I-1)*IX+1, (J-1)*IY+1),X(IP),Y(IP),
     $                        Z(IP))
                  SICHTB = SICHTB + T((I-1)*IX+1, (J-1)*IY+1)
                  NSI    = NSI + 1
*
               ELSE
                  X(IP)  = XT((I-1)*IX+1)
                  Y(IP)  = YT((J-1)*IY+1)
                  Z(IP)  = T((I-1)*IX+1, (J-1)*IY+1)
               END IF
*
   40       CONTINUE
         END IF
*
         IF (.NOT.DRAHT) THEN
            IP     = INDP0
            DO 50 J = 1,NJ - 1
               IP     = IP + 1
               INDF   = INDF + 1
               IFACE1(INDF) = IP
               IFACE2(INDF) = IP + NJ
               IFACE3(INDF) = IP + NJ + 1
               IFACE4(INDF) = IP + 1
   50       CONTINUE
         END IF
*
         IP     = IP + 1
         IF (KIND.EQ.1) THEN
            X(IP)  = PT(1, (I-1)*IX+1, (NJ-1)*IY+1)
            Y(IP)  = PT(2, (I-1)*IX+1, (NJ-1)*IY+1)
            Z(IP)  = PT(3, (I-1)*IX+1, (NJ-1)*IY+1)
*
         ELSE IF (KIND.EQ.2) THEN
            X(IP)  = XX((I-1)*IX+1, (NJ-1)*IY+1)
            Y(IP)  = YY((I-1)*IX+1, (NJ-1)*IY+1)
            Z(IP)  = ZZ((I-1)*IX+1, (NJ-1)*IY+1)
*
         ELSE
            IF (T((I-1)*IX+1, (NJ-1)*IY+1).LT.BLIND1 .OR.
     $          T((I-1)*IX+1, (NJ-1)*IY+1).GT.BLIND2) THEN
               CALL GR3TRC(XT((I-1)*IX+1),YT((NJ-1)*IY+1),
     $                     T((I-1)*IX+1, (NJ-1)*IY+1),X(IP),Y(IP),Z(IP))
               SICHTB = SICHTB + T((I-1)*IX+1, (NJ-1)*IY+1)
               NSI    = NSI + 1
*
            ELSE
               X(IP)  = XT((I-1)*IX+1)
               Y(IP)  = YT((NJ-1)*IY+1)
               Z(IP)  = T((I-1)*IX+1, (NJ-1)*IY+1)
            END IF
*
         END IF
*
         IF (Z(IP).LT.BLIND1 .OR. Z(IP).GT.BLIND2) THEN
            SICHTB = SICHTB + Z(IP)
            NSI    = NSI + 1
         END IF
*
   60 CONTINUE
      IF (KIND.EQ.1) THEN
         DO 70 J = 1,NJ
            IP     = IP + 1
            X(IP)  = PT(1, (NI-1)*IX+1, (J-1)*IY+1)
            Y(IP)  = PT(2, (NI-1)*IX+1, (J-1)*IY+1)
            Z(IP)  = PT(3, (NI-1)*IX+1, (J-1)*IY+1)
            IF (Z(IP).LT.BLIND1 .OR. Z(IP).GT.BLIND2) THEN
               SICHTB = SICHTB + Z(IP)
               NSI    = NSI + 1
            END IF
*
   70    CONTINUE
*
      ELSE IF (KIND.EQ.2) THEN
         DO 80 J = 1,NJ
            IP     = IP + 1
            X(IP)  = XX((NI-1)*IX+1, (J-1)*IY+1)
            Y(IP)  = YY((NI-1)*IX+1, (J-1)*IY+1)
            Z(IP)  = ZZ((NI-1)*IX+1, (J-1)*IY+1)
            IF (Z(IP).LT.BLIND1 .OR. Z(IP).GT.BLIND2) THEN
               SICHTB = SICHTB + Z(IP)
               NSI    = NSI + 1
            END IF
*
   80    CONTINUE
*
      ELSE
         DO 90 J = 1,NJ
            IP     = IP + 1
            IF (T((NI-1)*IX+1, (J-1)*IY+1).LT.BLIND1 .OR.
     $          T((NI-1)*IX+1, (J-1)*IY+1).GT.BLIND2) THEN
               CALL GR3TRC(XT((NI-1)*IX+1),YT((J-1)*IY+1),
     $                     T((NI-1)*IX+1, (J-1)*IY+1),X(IP),Y(IP),Z(IP))
               SICHTB = SICHTB + T((NI-1)*IX+1, (J-1)*IY+1)
               NSI    = NSI + 1
*
            ELSE
               X(IP)  = XT((NI-1)*IX+1)
               Y(IP)  = YT((J-1)*IY+1)
               Z(IP)  = T((NI-1)*IX+1, (J-1)*IY+1)
            END IF
*
   90    CONTINUE
      END IF

      KORR   = NSI .EQ. (IP-IP0) .OR. DRAHT
      SICHTB = SICHTB/NSI
      IF (.NOT.KORR) THEN
         INP    = IP0
         INF    = IF0
         DO 110 I = 1,NI - 1
            DO 100 J = 1,NJ - 1
               INP    = INP + 1
               INF    = INF + 1
               IF (Z(INP).GE.BLIND1 .AND. Z(INP).LE.BLIND2 .OR.
     $             Z(INP+NJ).GE.BLIND1 .AND. Z(INP+NJ).LE.BLIND2 .OR.
     $             Z(INP+NJ+1).GE.BLIND1 .AND.
     $             Z(INP+NJ+1).LE.BLIND2 .OR.
     $             Z(INP+1).GE.BLIND1 .AND. Z(INP+1).LE.BLIND2) THEN
                  IFACE1(INF) = IP
                  IFACE2(INF) = IP
                  IFACE3(INF) = IP
                  IFACE4(INF) = IP
               END IF
*
  100       CONTINUE
            INP    = INP + 1
  110    CONTINUE
      END IF
*
      IF (JGR(5,IFL).NE.-1000) THEN
         IF (IN1) THEN
            DO 130 JJ = 1,NJ
               DO 120 II = 1,NI - 1
                  IAL    = II*NJ + JJ + IP0
                  IF ((Z(IAL).LT.BLIND1.OR.Z(IAL).GT.BLIND2) .AND.
     $                (Z(IAL-NJ).LT.BLIND1.OR.Z(IAL-NJ).GT.BLIND2)) THEN
                     INDL   = INDL + 1
                     ILINE(1,INDL) = IAL - NJ
                     ILINE(2,INDL) = IAL
                     IF (.NOT.DRAHT) THEN
                        IF (JJ.EQ.NJ) THEN
                           IND1   = 0
*
                        ELSE
                           IND1   = II* (NJ-1) + JJ - NJ + IF0 + 1
                        END IF

                        IF (JJ.EQ.1) THEN
                           IND2   = 0
*
                        ELSE
                           IND2   = II* (NJ-1) + JJ - NJ + IF0
                        END IF
*
                        LNFC(1,INDL) = MAX(IND1,IND2)
                        LNFC(2,INDL) = MIN(IND1,IND2)
                     ELSE
                        LNFC(1,INDL) = 0
                        LNFC(2,INDL) = 0
                     END IF
*
                  END IF
*
  120          CONTINUE
  130       CONTINUE
         END IF
*
         IF (IN2) THEN
            DO 150 II = 1,NI
               DO 140 JJ = 1,NJ - 1
                  IAL    = II*NJ + JJ - NJ + IP0
                  IF ((Z(IAL).LT.BLIND1.OR.Z(IAL).GT.BLIND2) .AND.
     $                (Z(IAL+1).LT.BLIND1.OR.Z(IAL+1).GT.BLIND2)) THEN
                     INDL   = INDL + 1
                     ILINE(1,INDL) = IAL
                     ILINE(2,INDL) = IAL + 1
                     IF (.NOT.DRAHT) THEN
                        IF (II.LT.NI) THEN
                           IND1   = JJ + II* (NJ-1) - NJ + IF0 + 1
*
                        ELSE
                           IND1   = 0
                        END IF

                        IF (II.EQ.1) THEN
                           IND2   = 0
*
                        ELSE
                           IND2   = JJ + (II-1)* (NJ-1) - NJ + IF0 + 1
                        END IF

                        LNFC(1,INDL) = MAX(IND1,IND2)
                        LNFC(2,INDL) = MIN(IND1,IND2)
                     ELSE
                        LNFC(1,INDL) = 0
                        LNFC(2,INDL) = 0
                     END IF
*
                  END IF
*
  140          CONTINUE
  150       CONTINUE
         END IF
*
      END IF
*
      DO 160 I = IP0 + 1,IP
         IF (Z(I).GE.BLIND1 .AND. Z(I).LE.BLIND2) THEN
            Z(I)   = SICHTB
            CALL GR3TRC(X(I),Y(I),Z(I),XU,YU,ZU)
            X(I)   = XU
            Y(I)   = YU
            Z(I)   = ZU
         END IF
*
  160 CONTINUE

      IL0    = INDL
      IF0    = INDF
      IP0    = IP
      END
      SUBROUTINE GR3RTX(ROT,THETA,AX)

c     finds the rotation which rotates a plane through theta
c     degrees. the axis of rotation is specified by ax.


c     .. scalar arguments ..
      REAL              THETA
      CHARACTER         AX
c     ..
c     .. array arguments ..
      REAL              ROT(3,3)
c     ..
c     .. local scalars ..
      REAL              FCOS,FSIN,THETA1,XX
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ATAN,COS,SIN
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

c  xx is used to convert degrees to radians.

      XX     = ATAN(1.)/45.

      THETA1 = XX*THETA

      FSIN   = SIN(THETA1)
      FCOS   = COS(THETA1)

      IF (AX.EQ.'Z') THEN

c     rotation about z coordinate axis.

         ROT(1,1) = FCOS
         ROT(1,2) = -FSIN
         ROT(1,3) = 0.0

         ROT(2,1) = FSIN
         ROT(2,2) = FCOS
         ROT(2,3) = 0.0

         ROT(3,1) = 0.0
         ROT(3,2) = 0.0
         ROT(3,3) = 1.0
*
      ELSE IF (AX.EQ.'Y') THEN

c  rotation about y coordinate axis

         ROT(1,1) = FCOS
         ROT(1,2) = 0.0
         ROT(1,3) = FSIN

         ROT(2,1) = 0.0
         ROT(2,2) = 1.0
         ROT(2,3) = 0.0

         ROT(3,1) = -FSIN
         ROT(3,2) = 0.0
         ROT(3,3) = FCOS
*
      ELSE

c  rotation about x coordinate axis

         ROT(1,1) = 1.0
         ROT(1,2) = 0.0
         ROT(1,3) = 0.0

         ROT(2,1) = 0.0
         ROT(2,2) = FCOS
         ROT(2,3) = -FSIN

         ROT(3,1) = 0.0
         ROT(3,2) = FSIN
         ROT(3,3) = FCOS
      END IF
*
      END
      SUBROUTINE GR3TRF(ROT,X,Y,Z,NP,INC)


c     .. scalar arguments ..
      INTEGER           INC,NP
c     ..
c     .. array arguments ..
      REAL              ROT(3,3),X(NP),Y(NP),Z(NP)
c     ..
c     .. local scalars ..
      REAL              TEMP1,TEMP2,TEMP3
      INTEGER           I
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      DO 10 I = 1,NP,INC
         TEMP1  = ROT(1,1)*X(I) + ROT(1,2)*Y(I) + ROT(1,3)*Z(I)
         TEMP2  = ROT(2,1)*X(I) + ROT(2,2)*Y(I) + ROT(2,3)*Z(I)
         TEMP3  = ROT(3,1)*X(I) + ROT(3,2)*Y(I) + ROT(3,3)*Z(I)
         X(I)   = TEMP1
         Y(I)   = TEMP2
         Z(I)   = TEMP3
   10 CONTINUE

      END
      SUBROUTINE GR3CUB(X,Y,Z,LL,LLV,KL,KP,EX33)

c     .. scalar arguments ..
      INTEGER           KL,KP
c     ..
c     .. common blocks ..
      INTEGER           IF1,IF2,IF3,IF4,IFO,IL,ILO,IPO,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,ILO,IFO,IPO,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
c     ..
c     .. array arguments ..
      REAL              EX33(3,3),X(MNP),Y(MNP),Z(MNP)
      INTEGER           LL(2,MNL),LLV(2,MNL)
c     ..
c     .. local scalars ..
      INTEGER           I,J,K,L
c     ..
c     .. local arrays ..
      INTEGER           N1(12),N2(12)
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. data statements ..
c die reihenfolge ist nur so richtig !  siehe gr3plt !
      DATA              N1/1,2,4,3,1,5,6,2,1,5,7,3/,N2/5,6,8,7,3,7,8,4,
     $                  2,6,8,4/
c     ..
c     .. executable statements ..

c---- berechnung der 8 punkte des kubus

      I      = KP
      DO 30 J = 1,3,2
         DO 20 K = 1,3,2
            DO 10 L = 1,3,2
               I      = I + 1
               X(I)   = EX33(1,J) + (EX33(1,3)-EX33(1,1))*0.001* (J-2)
               Y(I)   = EX33(2,K) + (EX33(2,3)-EX33(2,1))*0.001* (K-2)
               Z(I)   = EX33(3,L) + (EX33(3,3)-EX33(3,1))*0.001* (L-2)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE

c---- berechnung der 12 linien des kubus

      DO 40 I = 1,12
         LL(1,KL+I) = KP + N1(I)
         LL(2,KL+I) = KP + N2(I)
         LLV(1,KL+I) = 0
         LLV(2,KL+I) = 0
   40 CONTINUE

      KP     = KP + 8
      KL     = KL + 12

      END
      SUBROUTINE GR3INV(AS)

c     .. parameters ..
      INTEGER           N
      PARAMETER         (N=3)
c     ..
c     .. array arguments ..
      REAL              AS(N,N)
c     ..
c     .. local scalars ..
      DOUBLE PRECISION  AH,AH1,AMAX,AH1U(4)
      INTEGER           I,IMAX,J,K,J1S
c     ..
c     .. local arrays ..
      DOUBLE PRECISION  A(N,N)
      INTEGER           PIV(N)
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS
c     ..
c     .. executable statements ..

      DO 20 I = 1,N
         DO 10 J = 1,N
            A(I,J) = AS(I,J)
   10    CONTINUE
   20 CONTINUE
      DO 80 K = 1,N
         IMAX   = K
         IF ( K.NE.N ) THEN
            AMAX  = ABS(A(K,K))
            DO 30 I = 1, N - K
               AH1U(I) = ABS(A(K+I,K))
   30       CONTINUE

            J1S = 1
            AH = AH1U(1)
            DO 31 I=2,3-K
               IF ( AH1U(I).GT.AH ) THEN
                  AH = AH1U(I)
                  J1S = I
               ENDIF
   31       CONTINUE
            IF (AH1U(J1S) .GT. AMAX) THEN
               AMAX = AH1U(J1S)
               IMAX = J1S + K
            ENDIF
            DO 40 J = 1,N
               AH     = A(IMAX,J)
               A(IMAX,J) = A(K,J)
               A(K,J) = AH
   40       CONTINUE
         ENDIF
         PIV(K) = IMAX
         AH     = 1E0/A(K,K)
         DO 50 J = 1,N
            A(K,J) = A(K,J)*AH
   50    CONTINUE
         DO 70 I = 1,N
            IF (I.EQ.K) GO TO 70
            AH1    = A(I,K)
            DO 60 J = 1,N
               A(I,J) = A(I,J) - A(K,J)*AH1
   60       CONTINUE
            A(I,K) = -AH1*AH
   70    CONTINUE
         A(K,K) = AH
   80 CONTINUE
      DO 100 J = N - 1,1,-1
         K      = PIV(J)
         DO 90 I = 1,N
            AH     = A(I,K)
            A(I,K) = A(I,J)
            A(I,J) = AH
   90    CONTINUE
  100 CONTINUE
      DO 120 I = 1,N
         DO 110 J = 1,N
            AS(I,J) = A(I,J)
  110    CONTINUE
  120 CONTINUE
      END
      SUBROUTINE GR3TRA(AR,IER,XYZ,N,IADJ,IEND,LINKS)
c     .. scalar arguments ..
      INTEGER           IER,LINKS,N
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP),XYZ(3,MNP)
      INTEGER           IADJ(6*N),IEND(N)
c     ..
c     .. local scalars ..
      INTEGER           I,IA,IB,J,K,KA,KE
c     ..
c     .. local arrays ..
      REAL              POI(3,3)
      INTEGER           LIN(3)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3OBP
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. data statements ..
      DATA              LIN/1,1,1/
c     ..
c     .. executable statements ..

c---- linien und flaechen
      KA     = 1
      DO 40 I = 1,N
         KE     = IEND(I)
         DO 20 J = KA,KE
            IA     = IADJ(J)
            IF (J.LT.KE) THEN
               IB     = IADJ(J+1)
               IF (IB.EQ.0) GO TO 30
*
            ELSE
               IB     = IADJ(KA)
            END IF
*
            IF (IA.GT.I .AND. IB.GT.I) THEN
               DO 10 K = 1,3
                  IF (LINKS.NE.0) THEN
                     POI(K,1) = XYZ(K,I)
                     POI(K,2) = XYZ(K,IA)
                     POI(K,3) = XYZ(K,IB)
*
                  ELSE
                     POI(K,1) = XYZ(K,I)
                     POI(K,2) = XYZ(K,IB)
                     POI(K,3) = XYZ(K,IA)
                  END IF
*
   10          CONTINUE
               CALL GR3OBP(AR,IER,POI,3,1,LIN)
            END IF
*
   20    CONTINUE
   30    CONTINUE
         IADJ(I) = 0
         KA     = KE + 1
   40 CONTINUE
      END
c     algorithm 624 collected algorithms from acm.
c     algorithm appeared in acm-trans. math. software, vol.10, no. 4,
c     dec., 1984, p. 437.
      SUBROUTINE GRMESH(N,Z,IADJ,IEND,IER)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine creates a thiessen triangulation of n
c arbitrarily spaced points in the plane referred to as
c nodes.  the triangulation is optimal in the sense that it
c is as nearly equiangular as possible.  grmesh is part of
c an interpolation package which also provides subroutines
c to reorder the nodes, add a new node, delete an arc, plot
c the mesh, and print the data structure.
c   unless the nodes are already ordered in some reasonable
c fashion, they should be reordered by subroutine reordr for
c increased efficiency before calling grmesh.
c
c input parameters -     n - number of nodes in the mesh.
c                            n .ge. 3.
c
c                        z - n-vectors of coordinates.
c                            (  z(3,i)   defines node i.
c
c                     iadj - vector of length .ge. 6*n-9.
c
c                     iend - vector of length .ge. n.
c
c n, x, and y are not altered by this routine.
c
c output parameters - iadj - adjacency lists of neighbors in
c                            counterclockwise order.  the
c                            list for node i+1 follows that
c                            for node i where x and y define
c                            the order.  the value 0 denotes
c                            the boundary (or a pseudo-node
c                            at infinity) and is always the
c                            last neighbor of a boundary
c                            node.  iadj is unchanged if ier
c                            .ne. 0.
c
c                     iend - pointers to the ends of
c                            adjacency lists (sets of
c                            neighbors) in iadj.  the
c                            neighbors of node 1 begin in
c                            iadj(1).  for k .gt. 1, the
c                            neighbors of node k begin in
c                            iadj(iend(k-1)+1) and k has
c                            iend(k) - iend(k-1) neighbors
c                            including (possibly) the
c                            boundary.  iadj(iend(k)) .eq. 0
c                            iff node k is on the boundary.
c                            iend is unchanged if ier = 1.
c                            if ier = 2 iend contains the
c                            indices of a sequence of n
c                            nodes ordered from left to
c                            right where left and right are
c                            defined by assuming node 1 is
c                            to the left of node 2.
c
c                      ier - error indicator
c                            ier = 0 if no errors were
c                                    encountered.
c                            ier = 1 if n .lt. 3.
c                            ier = 2 if n .ge. 3 and all
c                                    nodes are collinear.
c
c modules referenced by grmesh - grshfd, gradno, grfind,
c                                griadd, grbdya, grswpt,
c                                grswap, grind
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           IER,N
c     ..
c     .. array arguments ..
      REAL              Z(3,N)
      INTEGER           IADJ(6*N-9),IEND(N)
c     ..
c     .. local scalars ..
      REAL              CPROD,DXK,DXR,DYK,DYR,SPROD,XK,XL,XR,YK,YL,YR
      INTEGER           I,IERR,IND,INDX,ITEMP,K,KM1,KM1D2,KMI,KMIN,N0,
     $                  NL,NN,NR
c     ..
c     .. external subroutines ..
      EXTERNAL          GRADNO,GRSHFD
c     ..
c     .. executable statements ..
c
c local parameters -
c
c nn =          local copy of n
c k =           node (index) to be inserted into iend
c km1 =         k-1 - (variable) length of iend
c nl,nr =       iend(1), iend(km1) -- leftmost and rightmost
c                 nodes in iend as viewed from the right of
c                 1-2 when iend contains the initial ordered
c                 set of nodal indices
c xl,yl,xr,yr = x and y coordinates of nl and nr
c dxr,dyr =     xr-xl, yr-yl
c xk,yk =       x and y coordinates of node k
c dxk,dyk =     xk-xl, yk-yl
c cprod =       vector cross product of nl-nr and nl-k --
c                 used to determine the position of node k
c                 with respect to the line defined by the
c                 nodes in iend
c sprod =       scalar product used to determine the
c                 interval containing node k when k is on
c                 the line defined by the nodes in iend
c ind,indx =    indices for iend and iadj, respectively
c n0,itemp =    temporary nodes (indices)
c ierr =        dummy parameter for call to gradno
c km1d2,kmi,i = km1/2, k-i, do-loop index -- used in iend
c                 reordering loop
c kmin =        first node index sent to gradno
c
      NN     = N
      IER    = 1
      IF (NN.LT.3) GO TO 140
      IER    = 0
c
c initialize iend, nl, nr, and k
c
      IEND(1) = 1
      IEND(2) = 2
      XL     = Z(1,1)
      YL     = Z(2,1)
      XR     = Z(1,2)
      YR     = Z(2,2)
      K      = 2
c
c begin loop on nodes 3,4,...
c
   10 CONTINUE
      DXR    = XR - XL
      DYR    = YR - YL
c
c next loop begins here if nl and nr are unchanged
c
   20 CONTINUE
      IF (K.EQ.NN) GO TO 130
      KM1    = K
      K      = KM1 + 1
      XK     = Z(1,K)
      YK     = Z(2,K)
      DXK    = XK - XL
      DYK    = YK - YL
      CPROD  = DXR*DYK - DXK*DYR
      IF (CPROD.GT.0.) GO TO 60
      IF (CPROD.LT.0.) GO TO 80
c
c node k lies on the line containing nodes 1,2,...,k-1.
c   set sprod to (nl-nr,nl-k).
c
      SPROD  = DXR*DXK + DYR*DYK
      IF (SPROD.GT.0.) GO TO 30
c
c node k is to the left of nl.  insert k as the first
c   (leftmost) node in iend and set nl to k.
c
      CALL GRSHFD(1,KM1,1,IEND)
      IEND(1) = K
      XL     = XK
      YL     = YK
      GO TO 10
c
c node k is to the right of nl.  find the leftmost node
c   n0 which lies to the right of k.
c   set sprod to (n0-nl,n0-k).
c
   30 CONTINUE
      DO 40 IND = 2,KM1
         N0     = IEND(IND)
         SPROD  = (XL-Z(1,N0))* (XK-Z(1,N0)) +
     $            (YL-Z(2,N0))* (YK-Z(2,N0))
         IF (SPROD.GE.0.) GO TO 50
   40 CONTINUE
c
c node k is to the right of nr.  insert k as the last
c   (rightmost) node in iend and set nr to k.
c
      IEND(K) = K
      XR     = XK
      YR     = YK
      GO TO 10
c
c node k lies between iend(ind-1) and iend(ind).  insert k
c   in iend.
c
   50 CONTINUE
      CALL GRSHFD(IND,KM1,1,IEND)
      IEND(IND) = K
      GO TO 20
c
c node k is to the left of nl-nr.  reorder iend so that nl
c   is the leftmost node as viewed from k.
c
   60 CONTINUE
      KM1D2  = KM1/2
      DO 70 I = 1,KM1D2
         KMI    = K - I
         ITEMP  = IEND(I)
         IEND(I) = IEND(KMI)
         IEND(KMI) = ITEMP
   70 CONTINUE
c
c node k is to the right of nl-nr.  create a triangulation
c   consisting of nodes 1,2,...,k.
c
   80 CONTINUE
      NL     = IEND(1)
      NR     = IEND(KM1)
c
c create the adjacency lists for the first k-1 nodes.
c   insert neighbors in reverse order.  each node has four
c   neighbors except nl and nr which have three.
c
      DO 90 IND = 1,KM1
         N0     = IEND(IND)
         INDX   = 4*N0
         IF (N0.GE.NL) INDX   = INDX - 1
         IF (N0.GE.NR) INDX   = INDX - 1
         IADJ(INDX) = 0
         INDX   = INDX - 1
         IF (IND.LT.KM1) IADJ(INDX) = IEND(IND+1)
         IF (IND.LT.KM1) INDX   = INDX - 1
         IADJ(INDX) = K
         IF (IND.EQ.1) GO TO 90
         IADJ(INDX-1) = IEND(IND-1)
   90 CONTINUE
c
c create the adjacency list for node k
c
      INDX   = 5*KM1 - 1
      IADJ(INDX) = 0
      DO 100 IND = 1,KM1
         INDX   = INDX - 1
         IADJ(INDX) = IEND(IND)
  100 CONTINUE
c
c replace iend elements with pointers to iadj
c
      INDX   = 0
      DO 110 IND = 1,KM1
         INDX   = INDX + 4
         IF (IND.EQ.NL .OR. IND.EQ.NR) INDX   = INDX - 1
         IEND(IND) = INDX
  110 CONTINUE
      INDX   = INDX + K
      IEND(K) = INDX
c
c add the remaining nodes to the triangulation
c
      IF (K.EQ.NN) GO TO 140
      KMIN   = K + 1
      DO 120 K = KMIN,NN
         CALL GRADNO(K,Z,IADJ,IEND,IERR)
  120 CONTINUE
      GO TO 140
c
c all nodes are collinear
c
  130 CONTINUE
      IER    = 2
  140 CONTINUE
      END
      SUBROUTINE GRSHFD(NFRST,NLAST,KK,IARR)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine shifts a set of contiguous elements of an
c integer array kk positions downward (upward if kk .lt. 0).
c the loops are unrolled in order to increase efficiency.
c
c input parameters - nfrst,nlast - bounds on the portion of
c                                  iarr to be shifted.  all
c                                  elements between and
c                                  including the bounds are
c                                  shifted unless nfrst .gt.
c                                  nlast, in which case no
c                                  shift occurs.
c
c                             kk - number of positions each
c                                  element is to be shifted.
c                                  if kk .lt. 0 shift up.
c                                  if kk .gt. 0 shift down.
c
c                           iarr - integer array of length
c                                  .ge. nlast + max(kk,0).
c
c nfrst, nlast, and kk are not altered by this routine.
c
c output parameter -        iarr - shifted array.
c
c modules referenced by grshfd - none
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           KK,NFRST,NLAST
c     ..
c     .. array arguments ..
      INTEGER           IARR(*)
c     ..
c     .. local scalars ..
      INTEGER           I,IBAK,IMAX,INC,INDX,K,NF,NL,NLP1,NS,NSL
c     ..
c     .. data statements ..
      DATA              INC/5/
c     ..
c     .. executable statements ..
c
c local parameters -
c
c inc =  do-loop increment (unrolling factor) -- if inc is
c          changed, statements must be added to or deleted
c          from the do-loops
c k =    local copy of kk
c nf =   local copy of nfrst
c nl =   local copy of nlast
c nlp1 = nl + 1
c ns =   number of shifts
c nsl =  number of shifts done in unrolled do-loop (multiple
c          of inc)
c i =    do-loop index and index for iarr
c ibak = index for downward shift of iarr
c indx = index for iarr
c imax = bound on do-loop index
c
      K      = KK
      NF     = NFRST
      NL     = NLAST
      IF (NF.GT.NL .OR. K.EQ.0) RETURN
      NLP1   = NL + 1
      NS     = NLP1 - NF
      NSL    = INC* (NS/INC)
      IF (K.LT.0) GO TO 40
c
c shift downward starting from the bottom
c
      IF (NSL.LE.0) GO TO 20
      DO 10 I = 1,NSL,INC
         IBAK   = NLP1 - I
         INDX   = IBAK + K
         IARR(INDX) = IARR(IBAK)
         IARR(INDX-1) = IARR(IBAK-1)
         IARR(INDX-2) = IARR(IBAK-2)
         IARR(INDX-3) = IARR(IBAK-3)
         IARR(INDX-4) = IARR(IBAK-4)
   10 CONTINUE
c
c perform the remaining ns-nsl shifts one at a time
c
   20 CONTINUE
      IBAK   = NLP1 - NSL
   30 CONTINUE
      IF (IBAK.LE.NF) RETURN
      IBAK   = IBAK - 1
      INDX   = IBAK + K
      IARR(INDX) = IARR(IBAK)
      GO TO 30
c
c shift upward starting from the top
c
   40 CONTINUE
      IF (NSL.LE.0) GO TO 60
      IMAX   = NLP1 - INC
      DO 50 I = NF,IMAX,INC
         INDX   = I + K
         IARR(INDX) = IARR(I)
         IARR(INDX+1) = IARR(I+1)
         IARR(INDX+2) = IARR(I+2)
         IARR(INDX+3) = IARR(I+3)
         IARR(INDX+4) = IARR(I+4)
   50 CONTINUE
c
c perform the remaining ns-nsl shifts one at a time
c
   60 CONTINUE
      I      = NSL + NF
   70 CONTINUE
      IF (I.GT.NL) RETURN
      INDX   = I + K
      IARR(INDX) = IARR(I)
      I      = I + 1
      GO TO 70
*
      END
      SUBROUTINE GRADNO(KK,Z,IADJ,IEND,IER)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine adds node kk to a triangulation of a set
c of points in the plane producing a new triangulation.  a
c sequence of edge swaps is then applied to the mesh,
c resulting in an optimal triangulation.  gradno is part
c of an interpolation package which also provides routines
c to initialize the data structure, plot the mesh, and
c delete arcs.
c
c input parameters -   kk - index of the node to be added
c                           to the mesh.  kk .ge. 4.
c
c                     z - vectors of coordinates of the
c                           nodes in the mesh.  (  z(i))
c                           defines node i for i = 1,..,kk.
c
c                    iadj - set of adjacency lists of nodes
c                           1,..,kk-1.
c
c                    iend - pointers to the ends of
c                           adjacency lists in iadj for
c                           each node in the mesh.
c
c iadj and iend may be created by grmesh.
c
c kk, x, and y are not altered by this routine.
c
c output parameters - iadj,iend - updated with the addition
c                                 of node kk as the last
c                                 entry.
c
c                           ier - error indicator
c                                 ier = 0 if no errors
c                                         were encountered.
c                                 ier = 1 if all nodes
c                                         (including kk) are
c                                         collinear.
c
c modules referenced by gradno - grfind, griadd, grbdya,
c                                grshfd, grind, grswpt,
c                                grswap
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           IER,KK
c     ..
c     .. array arguments ..
      REAL              Z(3,KK)
      INTEGER           IADJ(6*KK-9),IEND(KK)
c     ..
c     .. local scalars ..
      REAL              XK,YK
      INTEGER           I1,I2,I3,IN1,IND21,IND2F,INDK1,INDKF,INDKL,IO1,
     $                  IO2,K,KM1,NABOR1
c     ..
c     .. external functions ..
      INTEGER           GRIND
      LOGICAL           GRSWPT
      EXTERNAL          GRIND,GRSWPT
c     ..
c     .. external subroutines ..
      EXTERNAL          GRBDYA,GRFIND,GRIADD,GRSWAP
c     ..
c     .. executable statements ..
c
c local parameters -
c
c k =        local copy of kk
c km1 =      k - 1
c i1,i2,i3 = vertices of a triangle containing k
c indkf =    iadj index of the first neighbor of k
c indkl =    iadj index of the last neighbor of k
c nabor1 =   first neighbor of k before any swaps occur
c io1,io2 =  adjacent neighbors of k defining an arc to
c              be tested for a swap
c in1 =      vertex opposite k -- first neighbor of io2
c              which precedes io1.  in1,io1,io2 are in
c              counterclockwise order.
c indk1 =    index of io1 in the adjacency list for k
c ind2f =    index of the first neighbor of io2
c ind21 =    index of io1 in the adjacency list for io2
c xk,yk =      z(1,k),   z(2,k)
c
      IER    = 0
      K      = KK
c
c initialization
c
      KM1    = K - 1
      XK     = Z(1,K)
      YK     = Z(2,K)
c
c add node k to the mesh
c
      CALL GRFIND(KM1,XK,YK,Z,IADJ,IEND,I1,I2,I3)
      IF (I1.EQ.0) GO TO 50
      IF (I3.EQ.0) CALL GRBDYA(K,I1,I2,IADJ,IEND)
      IF (I3.NE.0) CALL GRIADD(K,I1,I2,I3,IADJ,IEND)
c
c initialize variables for optimization of the mesh
c
      INDKF  = IEND(KM1) + 1
      INDKL  = IEND(K)
      NABOR1 = IADJ(INDKF)
      IO2    = NABOR1
      INDK1  = INDKF + 1
      IO1    = IADJ(INDK1)
c
c begin loop -- find the vertex opposite k
c
   10 CONTINUE
      IND2F  = 1
      IF (IO2.NE.1) IND2F  = IEND(IO2-1) + 1
      IND21  = GRIND(IO2,IO1,IADJ,IEND)
      IF (IND2F.EQ.IND21) GO TO 20
      IN1    = IADJ(IND21-1)
      GO TO 30
c
c in1 is the last neighbor of io2
c
   20 CONTINUE
      IND21  = IEND(IO2)
      IN1    = IADJ(IND21)
      IF (IN1.EQ.0) GO TO 40
c
c swap test -- if a swap occurs, two new arcs are opposite k
c              and must be tested.  indk1 and indkf must be
c              decremented.
c
   30 CONTINUE
      IF (.NOT.GRSWPT(IN1,K,IO1,IO2,Z)) GO TO 40
      CALL GRSWAP(IN1,K,IO1,IO2,IADJ,IEND)
      IO1    = IN1
      INDK1  = INDK1 - 1
      INDKF  = INDKF - 1
      GO TO 10
c
c no swap occurred.  reset io2 and io1, and test for
c   termination.
c
   40 CONTINUE
      IF (IO1.EQ.NABOR1) RETURN
      IO2    = IO1
      INDK1  = INDK1 + 1
      IF (INDK1.GT.INDKL) INDK1  = INDKF
      IO1    = IADJ(INDK1)
      IF (IO1.NE.0) GO TO 10
      RETURN
c
c all nodes are collinear
c
   50 CONTINUE
      IER    = 1
      END
      SUBROUTINE GRFIND(NST,PX,PY,Z,IADJ,IEND,I1,I2,I3)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine locates a point p in a thiessen triangu-
c lation, returning the vertex indices of a triangle which
c contains p.  grfind is part of an interpolation package
c which provides subroutines for creating the mesh.
c
c input parameters -    nst - index of node at which grfind
c                             begins search.  search time
c                             depends on the proximity of
c                             nst to p.
c
c                     px,py - x and y-coordinates of the
c                             point to be located.
c
c                         z - vectors of coordinates of
c                             nodes in the mesh.  (x(i),y(i))
c                             defines node i for i = 1,...,n
c                             where n .ge. 3.
c
c                      iadj - set of adjacency lists of
c                             nodes in the mesh.
c
c                      iend - pointers to the ends of
c                             adjacency lists in iadj for
c                             each node in the mesh.
c
c iadj and iend may be created by grmesh.
c
c input parameters are not altered by this routine.
c
c output parameters - i1,i2,i3 - vertex indices in counter-
c                                clockwise order - vertices
c                                of a triangle containing p
c                                if p is an interior node.
c                                if p is outside of the
c                                boundary of the mesh, i1
c                                and i2 are the first (right
c                                -most) and last (leftmost)
c                                nodes which are visible
c                                from p, and i3 = 0.  if p
c                                and all of the nodes lie on
c                                a single line then i1 = i2
c                                = i3 = 0.
c
c modules referenced by grfind - none
c
c intrinsic function called by grfind - max0
c
c***********************************************************
c
c     .. scalar arguments ..
      IMPLICIT NONE
      REAL              PX,PY
      INTEGER           I1,I2,I3,NST
c     ..
c     .. array arguments ..
      REAL              Z(3,*)
      INTEGER           IADJ(*),IEND(*)
c     ..
c     .. local scalars ..
      REAL              XP,YP
      INTEGER           IND,INDX,N0,N1,N2,N3,N4,NEXT,NF,NL,I
c     ..
c     ..
c     ..
c     .. statement function definitions ..
c
c local parameters -
c
c xp,yp =     local variables containing px and py
c n0,n1,n2 =  nodes in counterclockwise order defining a
c               cone (with vertex n0) containing p
c n3,n4 =     nodes opposite n1-n2 and n2-n1, respectively
c indx,ind =  indices for iadj
c nf,nl =     first and last neighbors of n0 in iadj, or
c               first (rightmost) and last (leftmost) nodes
c               visible from p when p is outside the
c               boundary
c next =      candidate for i1 or i2 when p is outside of
c               the boundary
c left =      statement function which computes the sign of
c               a cross product (z-component).  left(x1,...,
c               y0) = .TRUE. iff (x0,y0) is on or to the
c               left of the vector from (x1,y1) to (x2,y2).
c
c     ..
c     .. executable statements ..
      XP     = PX
      YP     = PY
c
c initialize variables and find a cone containing p
c
      N0     = MAX(NST,1)
      DO 61 i=1,2000000000
      INDX   = IEND(N0)
      NL     = IADJ(INDX)
      INDX   = 1
      IF (N0.NE.1) INDX   = IEND(N0-1) + 1
      NF     = IADJ(INDX)
      N1     = NF
      IF (NL.EQ.0) THEN
c
c        n0 is a boundary node.  set nl to the last nonzero
c        neighbor of n0.
c
         IND    = IEND(N0) - 1
         NL     = IADJ(IND)
         IF (.NOT.LEFT(Z(1,N0),Z(2,N0),Z(1,NF),Z(2,NF),XP,YP)) THEN
c
c           p is outside the boundary
c
            NL     = N0
            GO TO 160
*
         ELSE
            IF (LEFT(Z(1,NL),Z(2,NL),Z(1,N0),Z(2,N0),XP,YP)) GO TO 40
c
c           p is outside the boundary and n0 is the rightmost
c           visible boundary node
c
            I1     = N0
            GO TO 180
         END IF
      END IF
c
c n0 is an interior node.  find n1.
c
      DO 30 INDX=INDX+1,2000000000
         IF (LEFT(Z(1,N0),Z(2,N0),Z(1,N1),Z(2,N1),XP,YP)) GO TO 39
         N1     = IADJ(INDX)
         IF (N1.EQ.NL) GO TO 70
  30  CONTINUE
c
c p is to the left of arc n0-n1.  initialize n2 to the next
c   neighbor of n0.
c
   39 INDX = INDX-1
   40 DO 41 INDX=INDX+1,2000000000
         N2     = IADJ(INDX)
         IF (.NOT.LEFT(Z(1,N0),Z(2,N0),Z(1,N2),Z(2,N2),XP,YP)) GO TO 80
         N1     = N2
         IF (N1.EQ.NL) GO TO 42
   41 CONTINUE
   42 IF (.NOT.LEFT(Z(1,N0),Z(2,N0),Z(1,NF),Z(2,NF),XP,YP)) GO TO 70
      IF (XP.EQ.Z(1,N0) .AND. YP.EQ.Z(2,N0)) GO TO 60
c
c p is left of or on arcs n0-nb for all neighbors nb
c   of n0.
c all points are collinear iff p is left of nb-n0 for
c   all neighbors nb of n0.  search the neighbors of n0
c   in reverse order.  note -- n1 = nl and indx points to
c   nl.
c
   50 CONTINUE
      IF (.NOT.LEFT(Z(1,N1),Z(2,N1),Z(1,N0),Z(2,N0),XP,YP)) GO TO 60
      IF (N1.EQ.NF) GO TO 200
      INDX   = INDX - 1
      N1     = IADJ(INDX)
      GO TO 50
c
c p is to the right of n1-n0, or p=n0.  set n0 to n1 and
c   start over.
c
   60 CONTINUE
      N0     = N1
   61 CONTINUE
c
c p is between arcs n0-n1 and n0-nf
c
   70 CONTINUE
      N2     = NF
c
c p is contained in a cone defined by line segments n0-n1
c   and n0-n2 where n1 is adjacent to n2
c
   80 CONTINUE
      N3     = N0
      DO 121 i=1,2000000000
      IF (LEFT(Z(1,N1),Z(2,N1),Z(1,N2),Z(2,N2),XP,YP)) GO TO 130
c
c set n4 to the first neighbor of n2 following n1
c
      INDX   = IEND(N2)
      IF (IADJ(INDX).EQ.N1) THEN
c
c        n1 is the last neighbor of n2.
c        set n4 to the first neighbor.
c
         INDX   = 1
         IF (N2.NE.1) INDX   = IEND(N2-1) + 1
         N4     = IADJ(INDX)
      ELSE
c
c        n1 is not the last neighbor of n2
c
CGroten  the call is done only to make OPT(3) STRICT possible
         CALL grdum(n2)
         DO 100 INDX = INDX-1,1,-1
            IF (IADJ(INDX).EQ.N1) GO TO 101
  100    CONTINUE
         STOP 'GRFIND'
C 100    CONTINUE
C        INDX   = INDX - 1
C        IF (IADJ(INDX).NE.N1) GO TO 100
  101    N4     = IADJ(INDX+1)
         IF (N4.EQ.0) THEN
c
c           p is outside the boundary
c
            NF     = N2
            NL     = N1
            GO TO 160
         END IF
      END IF
c
c define a new arc n1-n2 which intersects the line
c   segment n0-p
c
      IF (.NOT.LEFT(Z(1,N0),Z(2,N0),Z(1,N4),Z(2,N4),XP,YP)) THEN
         N3     = N2
         N2     = N4
      ELSE
         N3     = N1
         N1     = N4
      ENDIF
  121 CONTINUE
c
c p is in the triangle (n1,n2,n3) and not on n2-n3.  if
c   n3-n1 or n1-n2 is a boundary arc containing p, treat p
c   as exterior.
c
  130 CONTINUE
      INDX   = IEND(N1)
      IF (IADJ(INDX).NE.0) GO TO 150
c
c n1 is a boundary node.  n3-n1 is a boundary arc iff n3
c   is the last nonzero neighbor of n1.
c
      IF (N3.NE.IADJ(INDX-1)) GO TO 140
c
c n3-n1 is a boundary arc
c
      IF (.NOT.LEFT(Z(1,N1),Z(2,N1),Z(1,N3),Z(2,N3),XP,YP)) GO TO 140
c
c p lies on n1-n3
c
      I1     = N1
      I2     = N3
      I3     = 0
      RETURN
c
c n3-n1 is not a boundary arc containing p.  n1-n2 is a
c   boundary arc iff n2 is the first neighbor of n1.
c
  140 CONTINUE
      INDX   = 1
      IF (N1.NE.1) INDX   = IEND(N1-1) + 1
      IF (N2.NE.IADJ(INDX)) GO TO 150
c
c n1-n2 is a boundary arc
c
      IF (.NOT.LEFT(Z(1,N2),Z(2,N2),Z(1,N1),Z(2,N1),XP,YP)) GO TO 150
c
c p lies on n1-n2
c
      I1     = N2
      I2     = N1
      I3     = 0
      RETURN
c
c p does not lie on a boundary arc.
c
  150 CONTINUE
      I1     = N1
      I2     = N2
      I3     = N3
      RETURN
c
c nf and nl are adjacent boundary nodes which are visible
c   from p.  find the first visible boundary node.
c set next to the first neighbor of nf.
c
  160 CONTINUE
      INDX   = 1
      IF (NF.NE.1) INDX   = IEND(NF-1) + 1
      NEXT   = IADJ(INDX)
      IF (LEFT(Z(1,NF),Z(2,NF),Z(1,NEXT),Z(2,NEXT),XP,YP)) GO TO 170
      NF     = NEXT
      GO TO 160
c
c nf is the first (rightmost) visible boundary node
c
  170 CONTINUE
      I1     = NF
c
c find the last visible boundary node.  nl is the first
c   candidate for i2.
c set next to the last neighbor of nl.
c
  180 CONTINUE
      INDX   = IEND(NL) - 1
      NEXT   = IADJ(INDX)
      IF (LEFT(Z(1,NEXT),Z(2,NEXT),Z(1,NL),Z(2,NL),XP,YP)) GO TO 190
      NL     = NEXT
      GO TO 180
c
c nl is the last (leftmost) visible boundary node
c
  190 CONTINUE
      I2     = NL
      I3     = 0
      RETURN
c
c all points are collinear
c
  200 CONTINUE
      I1     = 0
      I2     = 0
      I3     = 0


      contains

      function left(X1,Y1,X2,Y2,X0,Y0) result(erg)
      REAL,intent(in)::   X0,X1,X2,Y0,Y1,Y2
      LOGICAL erg
        
      erg = (X2-X1)* (Y0-Y1) >= (X0-X1)* (Y2-Y1)
      end function left

      END SUBROUTINE GRFIND

      SUBROUTINE GRDUM(I)
C     Nonsense
      I = SIGN(ABS(I),I)
      END
      SUBROUTINE GRIADD(KK,I1,I2,I3,IADJ,IEND)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine adds an interior node to a triangulation
c of a set of kk-1 points in the plane.  iadj and iend are
c updated with the insertion of node kk in the triangle
c whose vertices are i1, i2, and i3.
c
c input parameters -        kk - index of node to be
c                                inserted.  kk .ge. 4.
c
c                     i1,i2,i3 - indices of the vertices of
c                                a triangle containing node
c                                kk -- in counterclockwise
c                                order.
c
c                         iadj - set of adjacency lists
c                                of nodes in the mesh.
c
c                         iend - pointers to the ends of
c                                adjacency lists in iadj for
c                                each node in the mesh.
c
c   iadj and iend may be created by grmesh and must contain
c the vertices i1, i2, and i3.  i1,i2,i3 may be determined
c by grfind.
c
c kk, i1, i2, and i3 are not altered by this routine.
c
c output parameters - iadj,iend - updated with the addition
c                                 of node kk as the last
c                                 entry.  node kk will be
c                                 connected to nodes i1, i2,
c                                 and i3.  no optimization
c                                 of the mesh is performed.
c
c module referenced by griadd - grshfd
c
c intrinsic function called by griadd - mod
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           I1,I2,I3,KK
c     ..
c     .. array arguments ..
      INTEGER           IADJ(6*KK-9),IEND(KK)
c     ..
c     .. local scalars ..
      INTEGER           I,IMAX,IMIN,INDX,IP1,IP2,IP3,ITEMP,K,KM1,N1,N2,
     $                  NF,NL
c     ..
c     .. local arrays ..
      INTEGER           N(3),NFT(3)
c     ..
c     .. external subroutines ..
      EXTERNAL          GRSHFD
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MOD
c     ..
c     .. executable statements ..
c
c local parameters -
c
c k =           local copy of kk
c km1 =         k - 1
c n =           vector containing i1, i2, i3
c nft =         pointers to the tops of the 3 sets of iadj
c                 elements to be shifted downward
c ip1,ip2,ip3 = permutation indices for n and nft
c indx =        index for iadj and n
c nf,nl =       indices of first and last entries in iadj
c                 to be shifted down
c n1,n2 =       first 2 vertices of a new triangle --
c                 (n1,n2,kk)
c imin,imax =   bounds on do-loop index -- first and last
c                 elements of iend to be incremented
c i =           do-loop index
c itemp =       temporary storage location
c
      K      = KK
c
c initialization
c
      N(1)   = I1
      N(2)   = I2
      N(3)   = I3
c
c set up nft
c
      DO 20 I = 1,3
         N1     = N(I)
         INDX   = MOD(I,3) + 1
         N2     = N(INDX)
         INDX   = IEND(N1) + 1
c
c find the index of n2 as a neighbor of n1
c
   10    CONTINUE
         INDX   = INDX - 1
         IF (IADJ(INDX).NE.N2) GO TO 10
         NFT(I) = INDX + 1
   20 CONTINUE
c
c order the vertices by decreasing magnitude.
c   n(ip(i+1)) precedes n(ip(i)) in iend for
c   i = 1,2.
c
      IP1    = 1
      IP2    = 2
      IP3    = 3
      IF (N(2).LE.N(1)) GO TO 30
      IP1    = 2
      IP2    = 1
   30 CONTINUE
      IF (N(3).LE.N(IP1)) GO TO 40
      IP3    = IP1
      IP1    = 3
   40 CONTINUE
      IF (N(IP3).LE.N(IP2)) GO TO 50
      ITEMP  = IP2
      IP2    = IP3
      IP3    = ITEMP
c
c add node k to the adjacency lists of each vertex and
c   update iend.  for each vertex, a set of iadj elements
c   is shifted downward and k is inserted.  shifting starts
c   at the end of the array.
c
   50 CONTINUE
      KM1    = K - 1
      NL     = IEND(KM1)
      NF     = NFT(IP1)
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,3,IADJ)
      IADJ(NF+2) = K
      IMIN   = N(IP1)
      IMAX   = KM1
      DO 60 I = IMIN,IMAX
         IEND(I) = IEND(I) + 3
   60 CONTINUE
c
      NL     = NF - 1
      NF     = NFT(IP2)
      CALL GRSHFD(NF,NL,2,IADJ)
      IADJ(NF+1) = K
      IMAX   = IMIN - 1
      IMIN   = N(IP2)
      DO 70 I = IMIN,IMAX
         IEND(I) = IEND(I) + 2
   70 CONTINUE
c
      NL     = NF - 1
      NF     = NFT(IP3)
      CALL GRSHFD(NF,NL,1,IADJ)
      IADJ(NF) = K
      IMAX   = IMIN - 1
      IMIN   = N(IP3)
      DO 80 I = IMIN,IMAX
         IEND(I) = IEND(I) + 1
   80 CONTINUE
c
c add node k to iend and its neighbors to iadj
c
      INDX   = IEND(KM1)
      IEND(K) = INDX + 3
      DO 90 I = 1,3
         INDX   = INDX + 1
         IADJ(INDX) = N(I)
   90 CONTINUE
      END
      SUBROUTINE GRBDYA(KK,I1,I2,IADJ,IEND)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine adds a boundary node to a triangulation
c of a set of kk-1 points in the plane.  iadj and iend are
c updated with the insertion of node kk.
c
c input parameters -   kk - index of an exterior node to be
c                           added.  kk .ge. 4.
c
c                      i1 - first (rightmost as viewed from
c                           kk) boundary node in the mesh
c                           which is visible from kk - the
c                           line segment kk-i1 intersects
c                           no arcs.
c
c                      i2 - last (leftmost) boundary node
c                           which is visible from kk.
c
c                    iadj - set of adjacency lists of nodes
c                           in the mesh.
c
c                    iend - pointers to the ends of
c                           adjacency lists in iadj for
c                           each node in the mesh.
c
c   iadj and iend may be created by grmesh and must contain
c the vertices i1 and i2.  i1 and i2 may be determined by
c grfind.
c
c kk, i1, and i2 are not altered by this routine.
c
c output parameters - iadj,iend - updated with the addition
c                                 of node kk as the last
c                                 entry.  node kk will be
c                                 connected to i1, i2, and
c                                 all boundary nodes between
c                                 them.  no optimization of
c                                 the mesh is performed.
c
c module referenced by grbdya - grshfd
c
c intrinsic functions called by grbdya - min0, max0
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           I1,I2,KK
c     ..
c     .. array arguments ..
      INTEGER           IADJ(6*KK-9),IEND(KK)
c     ..
c     .. local scalars ..
      INTEGER           I,IMAX,IMIN,INDX,K,KEND,KM1,N1,N2,NEXT,NF,NL,
     $                  NLEFT,NRIGHT
c     ..
c     .. external subroutines ..
      EXTERNAL          GRSHFD
c     ..
c     ..
c     .. executable statements ..
c
c local parameters -
c
c k =            local copy of kk
c km1 =          k - 1
c nright,nleft = local copies of i1, i2
c nf,nl =        indices of iadj bounding the portion of the
c                  array to be shifted
c n1 =           iadj index of the first neighbor of nleft
c n2 =           iadj index of the last neighbor of nright
c i =            do-loop index
c imin,imax =    bounds on do-loop index -- first and last
c                  elements of iend to be incremented
c kend =         pointer to the last neighbor of k in iadj
c next =         next boundary node to be connected to kk
c indx =         index for iadj
c
      K      = KK
      KM1    = K - 1
      NRIGHT = I1
      NLEFT  = I2
c
c initialize variables
c
      NL     = IEND(KM1)
      N1     = 1
      IF (NLEFT.NE.1) N1     = IEND(NLEFT-1) + 1
      N2     = IEND(NRIGHT)
      NF     = MAX(N1,N2)
c
c insert k as a neighbor of max(1,nright,nleft)
c
      CALL GRSHFD(NF,NL,2,IADJ)
      IADJ(NF+1) = K
      IMIN   = MAX(NRIGHT,NLEFT)
      DO 10 I = IMIN,KM1
         IEND(I) = IEND(I) + 2
   10 CONTINUE
c
c initialize kend and insert k as a neighbor of
c   min(nright,nleft)
c
      KEND   = NL + 3
      NL     = NF - 1
      NF     = MIN(N1,N2)
      CALL GRSHFD(NF,NL,1,IADJ)
      IADJ(NF) = K
      IMAX   = IMIN - 1
      IMIN   = MIN(NRIGHT,NLEFT)
      DO 20 I = IMIN,IMAX
         IEND(I) = IEND(I) + 1
   20 CONTINUE
c
c insert nright as the first neighbor of k
c
      IADJ(KEND) = NRIGHT
c
c initialize indx for loop on boundary nodes between nright
c   and nleft
c
      INDX   = IEND(NRIGHT) - 2
   30 CONTINUE
      NEXT   = IADJ(INDX)
      IF (NEXT.EQ.NLEFT) GO TO 40
c
c connect next and k
c
      KEND   = KEND + 1
      IADJ(KEND) = NEXT
      INDX   = IEND(NEXT)
      IADJ(INDX) = K
      INDX   = INDX - 1
      GO TO 30
c
c insert nleft and 0 as the last neighbors of k
c
   40 CONTINUE
      IADJ(KEND+1) = NLEFT
      KEND   = KEND + 2
      IADJ(KEND) = 0
      IEND(K) = KEND
      END
      LOGICAL FUNCTION GRSWPT(IN1,IN2,IO1,IO2,Z)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this function decides whether or not to replace a
c diagonal arc in a quadrilateral with the other diagonal.
c the determination is based on the sizes of the angles
c contained in the 2 triangles defined by the diagonal.
c the diagonal is chosen to maximize the smallest of the
c six angles over the two pairs of triangles.
c
c input parameters -  in1,in2,io1,io2 - node indices of the
c                              four points defining the
c                              quadrilateral.  io1 and io2
c                              are currently connected by a
c                              diagonal arc.  this arc
c                              should be replaced by an arc
c                              connecting in1, in2 if the
c                              decision is made to swap.
c                              in1,io1,io2 must be in
c                              counterclockwise order.
c
c                          z - vectors of nodal coordinates.
c                              (  z( ,i)) are the coord-
c                              inates of node i for i = in1,
c                              in2, io1, or io2.
c
c none of the input parameters are altered by this routine.
c
c output parameter -  grswpt - .TRUE. iff the arc connecting
c                              io1 and io2 is to be replaced
c
c modules referenced by grswpt - none
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER                 IN1,IN2,IO1,IO2
c     ..
c     .. array arguments ..
      REAL                    Z(3,*)
c     ..
c     .. local scalars ..
      REAL                    COS1,COS2,DX11,DX12,DX21,DX22,DY11,DY12,
     $                        DY21,DY22,SIN1,SIN12,SIN2
c     ..
c     .. executable statements ..
c
c local parameters -
c
c dx11,dy11 = x,y coordinates of the vector in1-io1
c dx12,dy12 = x,y coordinates of the vector in1-io2
c dx22,dy22 = x,y coordinates of the vector in2-io2
c dx21,dy21 = x,y coordinates of the vector in2-io1
c sin1 =      cross product of the vectors in1-io1 and
c               in1-io2 -- proportional to sin(t1) where t1
c               is the angle at in1 formed by the vectors
c cos1 =      inner product of the vectors in1-io1 and
c               in1-io2 -- proportional to cos(t1)
c sin2 =      cross product of the vectors in2-io2 and
c               in2-io1 -- proportional to sin(t2) where t2
c               is the angle at in2 formed by the vectors
c cos2 =      inner product of the vectors in2-io2 and
c               in2-io1 -- proportional to cos(t2)
c sin12 =     sin1*cos2 + cos1*sin2 -- proportional to
c               sin(t1+t2)
c
      GRSWPT = .FALSE.
c
c compute the vectors containing the angles t1, t2
c
      DX11   = Z(1,IO1) - Z(1,IN1)
      DX12   = Z(1,IO2) - Z(1,IN1)
      DX22   = Z(1,IO2) - Z(1,IN2)
      DX21   = Z(1,IO1) - Z(1,IN2)
c
      DY11   = Z(2,IO1) - Z(2,IN1)
      DY12   = Z(2,IO2) - Z(2,IN1)
      DY22   = Z(2,IO2) - Z(2,IN2)
      DY21   = Z(2,IO1) - Z(2,IN2)
c
c compute inner products
c
      COS1   = DX11*DX12 + DY11*DY12
      COS2   = DX22*DX21 + DY22*DY21
c
c the diagonals should be swapped iff (t1+t2) .gt. 180
c   degrees.  the following two tests insure numerical
c   stability.
c
      IF (COS1.GE.0. .AND. COS2.GE.0.) RETURN
      IF (COS1.LT.0. .AND. COS2.LT.0.) GO TO 10
c
c compute vector cross products
c
      SIN1   = DX11*DY12 - DX12*DY11
      SIN2   = DX22*DY21 - DX21*DY22
      SIN12  = SIN1*COS2 + COS1*SIN2
      IF (SIN12.GE.0.) RETURN
   10 CONTINUE
      GRSWPT = .TRUE.
      END
      SUBROUTINE GRSWAP(NIN1,NIN2,NOUT1,NOUT2,IADJ,IEND)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this subroutine-swaps the diagonals in a convex quadri-
c lateral.
c
c input parameters -  nin1,nin2,nout1,nout2 - nodal indices
c                            of a pair of adjacent triangles
c                            which form a convex quadrilat-
c                            eral.  nout1 and nout2 are con-
c                            nected by an arc which is to be
c                            replaced by the arc nin1-nin2.
c                            (nin1,nout1,nout2) must be tri-
c                            angle vertices in counterclock-
c                            wise order.
c
c the above parameters are not altered by this routine.
c
c                iadj,iend - triangulation data structure
c                            (see subroutine grmesh).
c
c output parameters - iadj,iend - updated with the arc
c                                 replacement.
c
c modules referenced by grswap - grind, grshfd
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           NIN1,NIN2,NOUT1,NOUT2
c     ..
c     .. array arguments ..
      INTEGER           IADJ(*),IEND(*)
c     ..
c     .. local scalars ..
      INTEGER           I,IMAX,IMIN,IP1,IP2,J,K,NF,NL
c     ..
c     .. local arrays ..
      INTEGER           IN(2),IO(2)
c     ..
c     .. external functions ..
      INTEGER           GRIND
      EXTERNAL          GRIND
c     ..
c     .. external subroutines ..
      EXTERNAL          GRSHFD
c     ..
c     .. executable statements ..
c
c local parameters -
c
c in =        nin1 and nin2 ordered by increasing magnitude
c               (the neighbors of in(1) precede those of
c               in(2) in iadj)
c io =        nout1 and nout2 in increasing order
c ip1,ip2 =   permutation of (1,2) such that io(ip1)
c               precedes io(ip2) as a neighbor of in(1)
c j,k =       permutation of (1,2) used as indices of in
c               and io
c nf,nl =     iadj indices boundary a portion of the array
c               to be shifted
c i =         iend index
c imin,imax = bounds on the portion of iend to be incre-
c               mented or decremented
c
      IN(1)  = NIN1
      IN(2)  = NIN2
      IO(1)  = NOUT1
      IO(2)  = NOUT2
      IP1    = 1
c
c order the indices so that in(1) .lt. in(2) and io(1) .lt.
c   io(2), and choose ip1 and ip2 such that (in(1),io(ip1),
c   io(ip2)) forms a triangle.
c
      IF (IN(1).LT.IN(2)) GO TO 10
      IN(1)  = IN(2)
      IN(2)  = NIN1
      IP1    = 2
   10 CONTINUE
      IF (IO(1).LT.IO(2)) GO TO 20
      IO(1)  = IO(2)
      IO(2)  = NOUT1
      IP1    = 3 - IP1
   20 CONTINUE
      IP2    = 3 - IP1
      IF (IO(2).LT.IN(1)) GO TO 80
      IF (IN(2).LT.IO(1)) GO TO 120
c
c in(1) and io(1) precede in(2) and io(2).  for (j,k) =
c   (1,2) and (2,1), delete io(k) as a neighbor of io(j)
c   by shifting a portion of iadj either up or down and
c   and insert in(k) as a neighbor of in(j).
c
      DO 70 J = 1,2
         K      = 3 - J
         IF (IN(J).GT.IO(J)) GO TO 40
c
c   the neighbors of in(j) precede those of io(j) -- shift
c     down by 1
c
         NF     = 1 + GRIND(IN(J),IO(IP1),IADJ,IEND)
         NL     = -1 + GRIND(IO(J),IO(K),IADJ,IEND)
         IF (NF.LE.NL) CALL GRSHFD(NF,NL,1,IADJ)
         IADJ(NF) = IN(K)
         IMIN   = IN(J)
         IMAX   = IO(J) - 1
         DO 30 I = IMIN,IMAX
            IEND(I) = IEND(I) + 1
   30    CONTINUE
         GO TO 60
c
c   the neighbors of io(j) precede those of in(j) -- shift
c     up by 1
c
   40    CONTINUE
         NF     = 1 + GRIND(IO(J),IO(K),IADJ,IEND)
         NL     = -1 + GRIND(IN(J),IO(IP2),IADJ,IEND)
         IF (NF.LE.NL) CALL GRSHFD(NF,NL,-1,IADJ)
         IADJ(NL) = IN(K)
         IMIN   = IO(J)
         IMAX   = IN(J) - 1
         DO 50 I = IMIN,IMAX
            IEND(I) = IEND(I) - 1
   50    CONTINUE
c
c   reverse (ip1,ip2) for (j,k) = (2,1)
c
   60    CONTINUE
         IP1    = IP2
         IP2    = 3 - IP1
   70 CONTINUE
      RETURN
c
c the vertices are ordered (io(1),io(2),in(1),in(2)).
c   delete io(2) by shifting up by 1
c
   80 CONTINUE
      NF     = 1 + GRIND(IO(1),IO(2),IADJ,IEND)
      NL     = -1 + GRIND(IO(2),IO(1),IADJ,IEND)
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,-1,IADJ)
      IMIN   = IO(1)
      IMAX   = IO(2) - 1
      DO 90 I = IMIN,IMAX
         IEND(I) = IEND(I) - 1
   90 CONTINUE
c
c   delete io(1) by shifting up by 2 and insert in(2)
c
      NF     = NL + 2
      NL     = -1 + GRIND(IN(1),IO(IP2),IADJ,IEND)
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,-2,IADJ)
      IADJ(NL-1) = IN(2)
      IMIN   = IO(2)
      IMAX   = IN(1) - 1
      DO 100 I = IMIN,IMAX
         IEND(I) = IEND(I) - 2
  100 CONTINUE
c
c   shift up by 1 and insert in(1)
c
      NF     = NL + 1
      NL     = -1 + GRIND(IN(2),IO(IP1),IADJ,IEND)
      CALL GRSHFD(NF,NL,-1,IADJ)
      IADJ(NL) = IN(1)
      IMIN   = IN(1)
      IMAX   = IN(2) - 1
      DO 110 I = IMIN,IMAX
         IEND(I) = IEND(I) - 1
  110 CONTINUE
      RETURN
c
c the vertices are ordered (in(1),in(2),io(1),io(2)).
c   delete io(1) by shifting down by 1
c
  120 CONTINUE
      NF     = 1 + GRIND(IO(1),IO(2),IADJ,IEND)
      NL     = -1 + GRIND(IO(2),IO(1),IADJ,IEND)
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,1,IADJ)
      IMIN   = IO(1)
      IMAX   = IO(2) - 1
      DO 130 I = IMIN,IMAX
         IEND(I) = IEND(I) + 1
  130 CONTINUE
c
c   delete io(2) by shifting down by 2 and insert in(1)
c
      NL     = NF - 2
      NF     = 1 + GRIND(IN(2),IO(IP2),IADJ,IEND)
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,2,IADJ)
      IADJ(NF+1) = IN(1)
      IMIN   = IN(2)
      IMAX   = IO(1) - 1
      DO 140 I = IMIN,IMAX
         IEND(I) = IEND(I) + 2
  140 CONTINUE
c
c   shift down by 1 and insert in(2)
c
      NL     = NF - 1
      NF     = 1 + GRIND(IN(1),IO(IP1),IADJ,IEND)
      CALL GRSHFD(NF,NL,1,IADJ)
      IADJ(NF) = IN(2)
      IMIN   = IN(1)
      IMAX   = IN(2) - 1
      DO 150 I = IMIN,IMAX
         IEND(I) = IEND(I) + 1
  150 CONTINUE
      END
      INTEGER FUNCTION GRIND(NVERTX,NABOR,IADJ,IEND)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this function returns the index of nabor in the
c adjacency list for nvertx.
c
c input parameters - nvertx - node whose adjacency list is
c                             to be searched.
c
c                     nabor - node whose index is to be
c                             returned.  nabor must be
c                             connected to nvertx.
c
c                      iadj - set of adjacency lists.
c
c                      iend - pointers to the ends of
c                             adjacency lists in iadj.
c
c input parameters are not altered by this function.
c
c output parameter -  index - iadj(index) = nabor.
c
c modules referenced by grind - none
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER                NABOR,NVERTX
c     ..
c     .. array arguments ..
      INTEGER                IADJ(1),IEND(1)
c     ..
c     .. local scalars ..
      INTEGER                INDX,NB
c     ..
c     .. executable statements ..
c
c local parameters -
c
c nb =   local copy of nabor
c indx = index for iadj
c
      NB     = NABOR
c
c initialization
c
      INDX   = IEND(NVERTX) + 1
c
c search the list of nvertx neighbors for nb
c
   10 CONTINUE
      INDX   = INDX - 1
      IF (IADJ(INDX).NE.NB) GO TO 10
c
      GRIND  = INDX
      END
      SUBROUTINE GRDELN(NN,NOUT1,NOUT2,IADJ,IEND,IER)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   this routine deletes a boundary edge from a triangu-
c lation of a set of points in the plane.  it may be nec-
c essary to force certain edges to be present before call-
c ing delete (see subroutine edge).  note that subroutines
c edge, trfind, and the routines which call trfind (adnode,
c unif, intrc1, and intrc0) should not be called following
c a deletion.
c
c input parameters -    nn - number of nodes in the trian-
c                            gulation.
c
c              nout1,nout2 - pair of adjacent nodes on the
c                            boundary defining the arc to
c                            be removed.  nout2 must be the
c                            last nonzero neighbor of nout1.
c
c the above parameters are not altered by this routine.
c
c                iadj,iend - data structure defining the
c                            triangulation (see subroutine
c                            trmesh).
c
c output parameters - iadj,iend - updated with the removal
c                                 of the arc nout1-nout2
c                                 if ier .eq. 0.
c
c                           ier - error indicator
c                                 ier = 0 if no errors were
c                                         encountered.
c                                 ier = 1 if nout1 or nout2
c                                         is not on the
c                                         boundary.
c                                 ier = 2 if nout1 or nout2
c                                         has only 2 nonzero
c                                         neighbors.
c                                 ier = 3 if nout2 is not
c                                         the last neighbor
c                                         of nout1.
c                                 ier = 4 if a deletion
c                                         would divide the
c                                         mesh into two
c                                         regions.
c
c modules referenced by delete - grshfd, grind
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           IER,NN,NOUT1,NOUT2
c     ..
c     .. array arguments ..
      INTEGER           IADJ(6*NN-9),IEND(NN)
c     ..
c     .. local scalars ..
      INTEGER           I,IMAX,IND12,IND1F,IND1L,IND21,IND2F,IND2L,
     $                  INDFP2,INDLM3,INDN0,INDNF,INDNL,IO1,IO2,IOUT1,
     $                  IOUT2,ITEMP,N,NEWBD,NF,NL
c     ..
c     .. external functions ..
      INTEGER           GRIND
      EXTERNAL          GRIND
c     ..
c     .. external subroutines ..
      EXTERNAL          GRSHFD
c     ..
c     .. executable statements ..
c
c local parameters -
c
c n =           local copy of nn
c iout1,iout2 = local copies of nout1 and nout2
c io1,io2 =     nout1,nout2 in order of increasing magnitude
c ind12 =       index of io2 in the adjacency list for io1
c ind21 =       index of io1 in the adjacency list for io2
c itemp =       temporary storage location for permutations
c ind1f =       iadj index of the first neighbor of io1
c ind1l =       iadj index of the last neighbor of io1
c ind2f =       iadj index of the first neighbor of io2
c ind2l =       iadj index of the last neighbor of io2
c newbd =       the neighbor common to nout1 and nout2
c indnf =       iadj index of the first neighbor of newbd
c indnl =       iadj index of the last neighbor of newbd
c indn0 =       index of 0 in the adjacency list for newbd
c                 before permuting the neighbors
c indfp2 =      indnf + 2
c indlm3 =      indnl - 3
c nf,nl =       bounds on the portion of iadj to be shifted
c i =           do-loop index
c imax =        upper bound on do-loop for shifting iend
c
      N      = NN
      IOUT1  = NOUT1
      IOUT2  = NOUT2
c
c initialize indices
c
      IND1F  = 1
      IF (IOUT1.GT.1) IND1F  = IEND(IOUT1-1) + 1
      IND1L  = IEND(IOUT1)
      IND2F  = 1
      IF (IOUT2.GT.1) IND2F  = IEND(IOUT2-1) + 1
      IND2L  = IEND(IOUT2)
      NEWBD  = IADJ(IND1L-2)
      INDN0  = GRIND(NEWBD,IOUT2,IADJ,IEND)
      INDNL  = IEND(NEWBD)
c
c order vertices such that the neighbors of io1 precede
c   those of io2
c
      IF (IOUT1.GT.IOUT2) GO TO 10
      IO1    = IOUT1
      IO2    = IOUT2
      IND12  = IND1L - 1
      IND21  = IND2F
      GO TO 20
*
   10 CONTINUE
      IO1    = IOUT2
      IO2    = IOUT1
      IND12  = IND2F
      IND21  = IND1L - 1
c
c check for errors
c
   20 CONTINUE
      IF ((IADJ(IND1L).NE.0) .OR. (IADJ(IND2L).NE.0)) GO TO 210
      IF ((IND1L-IND1F.LE.2) .OR. (IND2L-IND2F.LE.2)) GO TO 220
      IF (IADJ(IND1L-1).NE.IOUT2) GO TO 230
      IF (IADJ(INDNL).EQ.0) GO TO 240
c
c delete the edge io1-io2 and make newbd a boundary node
c
      IF (NEWBD.LT.IO1) GO TO 80
      IF (NEWBD.LT.IO2) GO TO 60
c
c the vertices are ordered io1, io2, newbd.
c delete io2 as a neighbor of io1.
c
      NF     = IND12 + 1
      NL     = IND21 - 1
      CALL GRSHFD(NF,NL,-1,IADJ)
      IMAX   = IO2 - 1
      DO 30 I = IO1,IMAX
         IEND(I) = IEND(I) - 1
   30 CONTINUE
c
c delete io1 as a neighbor of io2
c
      NF     = NL + 2
      NL     = INDN0
      CALL GRSHFD(NF,NL,-2,IADJ)
      IMAX   = NEWBD - 1
      DO 40 I = IO2,IMAX
         IEND(I) = IEND(I) - 2
   40 CONTINUE
c
c shift the bottom of iadj up 1 leaving room for 0 as a
c   neighbor of newbd
c
      INDN0  = INDN0 - 1
      NF     = NL + 1
      NL     = IEND(N)
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,-1,IADJ)
      DO 50 I = NEWBD,N
         IEND(I) = IEND(I) - 1
   50 CONTINUE
      GO TO 120
c
c the vertices are ordered io1, newbd, io2.
c delete io2 as a neighbor of io1 leaving room for 0 as a
c   neighbor of newbd.
c
   60 CONTINUE
      NF     = IND12 + 1
      NL     = INDN0
      CALL GRSHFD(NF,NL,-1,IADJ)
      IMAX   = NEWBD - 1
      DO 70 I = IO1,IMAX
         IEND(I) = IEND(I) - 1
   70 CONTINUE
      GO TO 100
c
c the vertices are ordered newbd, io1, io2.
c delete io2 as a neighbor of io1 leaving room for 0 as a
c   neighbor of newbd.
c
   80 CONTINUE
      INDN0  = INDN0 + 1
      NF     = INDN0
      NL     = IND12 - 1
      IF (NF.LE.NL) CALL GRSHFD(NF,NL,1,IADJ)
      IMAX   = IO1 - 1
      DO 90 I = NEWBD,IMAX
         IEND(I) = IEND(I) + 1
   90 CONTINUE
c
c delete io1 as a neighbor of io2
c
  100 CONTINUE
      NF     = IND21 + 1
      NL     = IEND(N)
      CALL GRSHFD(NF,NL,-1,IADJ)
      DO 110 I = IO2,N
         IEND(I) = IEND(I) - 1
  110 CONTINUE
c
c permute the neighbors of newbd with end-around shifts so
c   that 0 is the last neighbor
c
  120 CONTINUE
      INDNF  = 1
      IF (NEWBD.GT.1) INDNF  = IEND(NEWBD-1) + 1
      INDNL  = IEND(NEWBD)
      IF (INDN0-INDNF.GE.INDNL-INDN0) GO TO 160
c
c shift upward
c
      IF (INDN0.GT.INDNF) GO TO 130
      CALL GRSHFD(INDNF+1,INDNL,-1,IADJ)
      GO TO 200
*
  130 CONTINUE
      INDFP2 = INDNF + 2
      IF (INDN0.LT.INDFP2) GO TO 150
      DO 140 I = INDFP2,INDN0
         ITEMP  = IADJ(INDNF)
         CALL GRSHFD(INDNF+1,INDNL,-1,IADJ)
         IADJ(INDNL) = ITEMP
  140 CONTINUE
c
c the last shift is by 2
c
  150 CONTINUE
      ITEMP  = IADJ(INDNF)
      CALL GRSHFD(INDFP2,INDNL,-2,IADJ)
      IADJ(INDNL-1) = ITEMP
      GO TO 200
c
c shift downward
c
  160 CONTINUE
      IF (INDN0.EQ.INDNL) GO TO 200
      IF (INDN0.LT.INDNL-1) GO TO 170
      CALL GRSHFD(INDNF,INDNL-2,1,IADJ)
      IADJ(INDNF) = IADJ(INDNL)
      GO TO 200
*
  170 CONTINUE
      INDLM3 = INDNL - 3
      IF (INDN0.GT.INDLM3) GO TO 190
      DO 180 I = INDN0,INDLM3
         ITEMP  = IADJ(INDNL)
         CALL GRSHFD(INDNF,INDNL-1,1,IADJ)
         IADJ(INDNF) = ITEMP
  180 CONTINUE
c
c the last shift is by 2
c
  190 CONTINUE
      ITEMP  = IADJ(INDNL-1)
      CALL GRSHFD(INDNF,INDLM3,2,IADJ)
      IADJ(INDNF+1) = IADJ(INDNL)
      IADJ(INDNF) = ITEMP
c
c insert 0 as the last neighbor of newbd
c
  200 CONTINUE
      IADJ(INDNL) = 0
      IER    = 0
      RETURN
c
c one of the vertices is not on the boundary
c
  210 CONTINUE
      IER    = 1
      RETURN
c
c one of the vertices has only two nonzero neighbors.  the
c   triangulation would be destroyed by a deletion
c
  220 CONTINUE
      IER    = 2
      RETURN
c
c nout2 is not the last nonzero neighbor of nout1
c
  230 CONTINUE
      IER    = 3
      RETURN
c
c a deletion would divide the mesh into two regions
c   connected at a single node
c
  240 CONTINUE
      IER    = 4
      END
      SUBROUTINE GREDGE(IN1,IN2,Q,LWK,IWK,IADJ,IEND,IER)
c
c***********************************************************
c
c                                               robert renka
c                                       oak ridge natl. lab.
c                                             (615) 576-5139
c
c   given a triangulation of n nodes and a pair of nodal
c indices in1 and in2, this routine swaps arcs as necessary
c to force in1 and in2 to be adjacent.  only arcs which
c intersect in1-in2 are swapped out.  if a thiessen triangu-
c lation is input, the resulting triangulation is as close
c as possible to a thiessen triangulation in the sense that
c all arcs other than in1-in2 are locally optimal.
c   a sequence of calls to edge may be used to force the
c presence of a set of edges defining the boundary of a non-
c convex region.  subsequent deletion of edges outside this
c region (by subroutine delete) results in a nonconvex tri-
c angulation which may serve as a finite element grid.
c (edge should not be called after a call to delete.)  if,
c on the other hand, interpolation is to be performed in the
c nonconvex region, edges must not be deleted, but it is
c still advantageous to have the nonconvex boundary present
c if it is desirable that interpolated values be influenced
c by the geometry.  note that subroutine getnp which is used
c to select the nodes entering into local derivative esti-
c mates will not necessarily return closest nodes if the
c triangulation has been rendered nonoptimal by a call to
c edge.  however, the effect will be merely to further en-
c hance the influence of the nonconvex geometry on interpo-
c lated values.
c
c input parameters - in1,in2 - indices (of q      ) in the
c                              range 1,...,n defining a pair
c                              of nodes to be connected by
c                              an arc.
c
c                        q(3,*) n-vectors containing carte-
c                              sian coordinates of the
c                              nodes.
c
c the above parameters are not altered by this routine.
c
c                        lwk - number of columns reserved
c                              for iwk.  this must be at
c                              least ni -- the number of
c                              arcs which intersect in1-in2.
c                              (ni is bounded by n-3).
c
c                        iwk - integer work array dimension-
c                              ed 2 by lwk (or vector of
c                              length .ge. 2*lwk).
c
c                  iadj,iend - data structure defining the
c                              triangulation.  see subrou-
c                              tine trmesh.
c
c output parameters - lwk - number of iwk columns required
c                           if ier = 0 or ier = 2.  lwk = 0
c                           iff in1 and in2 were adjacent
c                           on input.
c
c                     iwk - contains the indices of the end-
c                           points of the new arcs other
c                           than in1-in2 unless ier .gt. 0
c                           or lwk = 0.  new arcs to the
c                           left of in1->in2 are stored in
c                           the first k-1 columns (left por-
c                           tion of iwk), column k contains
c                           zeros, and new arcs to the right
c                           of in1->in2 occupy columns k+1,
c                           ...,lwk.  (k can be determined
c                           by searching iwk for the zeros.)
c
c               iadj,iend - updated if necessary to reflect
c                           the presence of an arc connect-
c                           ing in1 and in2, unaltered if
c                           ier .ne. 0.
c
c                     ier - error indicator
c                           ier = 0 if no errors were en-
c                                   countered.
c                           ier = 1 if in1 .lt. 1, in2 .lt.
c                                   1, in1 = in2, or lwk
c                                   .lt. 0 on input.
c                           ier = 2 if more space is requir-
c                                   ed in iwk.  see lwk.
c                           ier = 3 if in1 and in2 could not
c                                   be connected due to an
c                                   invalid data structure.
c
c modules referenced by edge - grswap, grind, grshfd, grswpt
c
c***********************************************************
c
c     .. scalar arguments ..
      INTEGER           IER,IN1,IN2,LWK
c     ..
c     .. array arguments ..
      REAL              Q(3,*)
      INTEGER           IADJ(*),IEND(*),IWK(2,LWK)
c     ..
c     .. local scalars ..
      REAL              X0,X1,X2,Y0,Y1,Y2
      INTEGER           I,INDF,INDL,INDX,IO1,IO2,IWC,IWCM1,IWCP1,IWEND,
     $                  IWF,IWL,LFT,N0,N1,N1LST,N2,NEXT,NL,NR
      LOGICAL           SWP
c     ..
c     .. external functions ..
      LOGICAL           GRSWPT
      EXTERNAL          GRSWPT
c     ..
c     .. external subroutines ..
      EXTERNAL          GRSWAP
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
c     ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. statement function definitions ..
c
c local parameters -
c
c n1,n2 =   local copies of in1 and in2 or nodes opposite an
c             arc io1-io2 to be tested for a swap in the
c             optimization loops
c iwend =   input or output value of lwk
c iwl =     iwk (column) index of the last (rightmost) arc
c             which intersects in1->in2
c indf =    iadj index of the first neighbor of in1 or io1
c indx =    iadj index of a neighbor of in1, nl, or io1
c n1lst =   last neighbor of in1
c nl,nr =   endpoints of an arc which intersects in1-in2
c             with nl left in1->in2
c next =    node opposite nl->nr
c iwf =     iwk (column) index of the first (leftmost) arc
c             which intersects in1->in2
c lft =     flag used to determine if a swap results in the
c             new arc intersecting in1-in2 -- lft = 0 iff
c             n0 = in1, lft = -1 implies n0 left in1->in2,
c             and lft = 1 implies n0 left in2->in1
c n0 =      node opposite nr->nl
c iwc =     iwk index between iwf and iwl -- nl->nr is
c             stored in iwk(1,iwc)->iwk(2,iwc)
c iwcp1 =   iwc + 1
c iwcm1 =   iwc - 1
c i =       do-loop index and column index for iwk
c io1,io2 = endpoints of an arc to be tested for a swap in
c             the optimization loops
c indl =    iadj index of the last neighbor of io1
c x1,y1 =   coordinates of in1
c x2,y2 =   coordinates of in2
c x0,y0 =   coordinates of n0
c swp =     flag set to .TRUE. iff a swap occurs in an opti-
c             mization loop
c left =    statement function which returns the value
c             .TRUE. iff (xp,yp) is on or to the left of the
c             vector (xa,ya)->(xb,yb)
c
c     ..
c     .. executable statements ..
c
c store in1, in2, and lwk in local variables and check for
c   errors.
c
      N1     = IN1
      N2     = IN2
      IWEND  = LWK
      IF (N1.LT.1 .OR. N2.LT.1 .OR. N1.EQ.N2 .OR.
     $    IWEND.LT.0) GO TO 350
c
c store the coordinates of n1 and n2 and initialize iwl.
c
      X1     = Q(1,N1)
      Y1     = Q(2,N1)
      X2     = Q(1,N2)
      Y2     = Q(2,N2)
      IWL    = 0
c
c set nr and nl to adjacent neighbors of n1 such that
c   nr left n2->n1 and nl left n1->n2.
c
c   set indf and indx to the indices of the first and last
c     neighbors of n1 and set n1lst to the last neighbor.
c
      INDF   = 1
      IF (N1.GT.1) INDF   = IEND(N1-1) + 1
      INDX   = IEND(N1)
      N1LST  = IADJ(INDX)
      IF (N1LST.EQ.0) INDX   = INDX - 1
      IF (N1LST.EQ.0) GO TO 20
c
c   n1 is an interior node.  loop through the neighbors nl
c     in reverse order until nl left n1->n2.
c
      NL     = N1LST
   10 CONTINUE
      IF (LEFT(X1,Y1,X2,Y2,Q(1,NL),Q(2,NL))) GO TO 20
      INDX   = INDX - 1
      NL     = IADJ(INDX)
      IF (INDX.GT.INDF) GO TO 10
c
c   nl is the first neighbor of n1.  set nr to the last
c     neighbor and test for an arc n1-n2.
c
      NR     = N1LST
      IF (NL.EQ.N2) GO TO 340
      GO TO 40
c
c   nl = iadj(indx) left n1->n2 and indx .gt. indf.  set
c     nr to the preceding neighbor of n1.
c
   20 CONTINUE
      INDX   = INDX - 1
      NR     = IADJ(INDX)
      IF (LEFT(X2,Y2,X1,Y1,Q(1,NR),Q(2,NR))) GO TO 30
      IF (INDX.GT.INDF) GO TO 20
c
c   set nl and nr to the first and last neighbors of n1 and
c     test for an invalid data structure (n1 cannot be a
c     boundary node and cannot be adjacent to n2).
c
      NL     = NR
      NR     = N1LST
      IF (NR.EQ.0 .OR. NR.EQ.N2) GO TO 370
      GO TO 40
c
c   set nl to the neighbor following nr and test for an arc
c     n1-n2.
c
   30 CONTINUE
      NL     = IADJ(INDX+1)
      IF (NL.EQ.N2 .OR. NR.EQ.N2) GO TO 340
c
c store the ordered sequence of intersecting edges nl->nr in
c   iwk(1,iwl)->iwk(2,iwl).
c
   40 CONTINUE
      IWL    = IWL + 1
      IF (IWL.LE.IWEND) IWK(1,IWL) = NL
      IF (IWL.LE.IWEND) IWK(2,IWL) = NR
c
c   set next to the neighbor of nl which follows nr.
c
      INDX   = IEND(NL)
      IF (IADJ(INDX).NE.NR) GO TO 50
c
c   nr is the last neighbor of nl.  set next to the first
c     neighbor.
c
      INDX   = 0
      IF (NL.NE.1) INDX   = IEND(NL-1)
      GO TO 60
c
c   nr is not the last neighbor of nl.  loop through the
c     neighbors in reverse order.
c
   50 CONTINUE
      INDX   = INDX - 1
      IF (IADJ(INDX).NE.NR) GO TO 50
c
c   store next, test for an invalid triangulation (nl->nr
c     cannot be a boundary edge), and test for termination
c     of the loop.
c
   60 CONTINUE
      NEXT   = IADJ(INDX+1)
      IF (NEXT.EQ.0) GO TO 370
      IF (NEXT.EQ.N2) GO TO 80
c
c   set nl or nr to next.
c
      IF (LEFT(X1,Y1,X2,Y2,Q(1,NEXT),Q(2,NEXT))) GO TO 70
      NR     = NEXT
      GO TO 40
*
   70 CONTINUE
      NL     = NEXT
      GO TO 40
c
c iwl is the number of arcs which intersect n1-n2.  store
c   lwk and test for sufficient space.
c
   80 CONTINUE
      LWK    = IWL
      IF (IWL.GT.IWEND) GO TO 360
      IWEND  = IWL
c
c initialize for edge swapping loop -- all possible swaps
c   are applied (even if the new arc again intersects
c   n1-n2), arcs to the left of n1->n2 are stored in the
c   left portion of iwk, and arcs to the right are stored in
c   the right portion.  iwf and iwl index the first and last
c   intersecting arcs.
c
      IER    = 0
      IWF    = 1
c
c top of loop -- set n0 to n1 and nl->nr to the first edge.
c   iwc points to the arc currently being processed.  lft
c   .le. 0 iff n0 left n1->n2.
c
   90 CONTINUE
      LFT    = 0
      N0     = N1
      X0     = X1
      Y0     = Y1
      NL     = IWK(1,IWF)
      NR     = IWK(2,IWF)
      IWC    = IWF
c
c   set next to the node opposite nl->nr unless iwc is the
c     last arc.
c
  100 CONTINUE
      IF (IWC.EQ.IWL) GO TO 210
      IWCP1  = IWC + 1
      NEXT   = IWK(1,IWCP1)
      IF (NEXT.NE.NL) GO TO 150
      NEXT   = IWK(2,IWCP1)
c
c   next right n1->n2 and iwc .lt. iwl.  test for a possible
c     swap.
c
      IF (.NOT.LEFT(X0,Y0,Q(1,NR),Q(2,NR),Q(1,NEXT),Q(2,
     $    NEXT))) GO TO 130
      IF (LFT.GE.0) GO TO 110
      IF (.NOT.LEFT(Q(1,NL),Q(2,NL),X0,Y0,Q(1,NEXT),Q(2,
     $    NEXT))) GO TO 130
c
c   replace nl->nr with n0->next.
c
      CALL GRSWAP(NEXT,N0,NL,NR,IADJ,IEND)
      IWK(1,IWC) = N0
      IWK(2,IWC) = NEXT
      GO TO 140
c
c   swap nl-nr for n0-next, shift columns iwc+1,...,iwl to
c     the left, and store n0-next in the right portion of
c     iwk.
c
  110 CONTINUE
      CALL GRSWAP(NEXT,N0,NL,NR,IADJ,IEND)
      DO 120 I = IWCP1,IWL
         IWK(1,I-1) = IWK(1,I)
         IWK(2,I-1) = IWK(2,I)
  120 CONTINUE
      IWK(1,IWL) = N0
      IWK(2,IWL) = NEXT
      IWL    = IWL - 1
      NR     = NEXT
      GO TO 100
c
c   a swap is not possible.  set n0 to nr.
c
  130 CONTINUE
      N0     = NR
      X0     = Q(1,N0)
      Y0     = Q(2,N0)
      LFT    = 1
c
c   advance to the next arc.
c
  140 CONTINUE
      NR     = NEXT
      IWC    = IWC + 1
      GO TO 100
c
c   next left n1->n2, next .ne. n2, and iwc .lt. iwl.
c     test for a possible swap.
c
  150 CONTINUE
      IF (.NOT.LEFT(Q(1,NL),Q(2,NL),X0,Y0,Q(1,NEXT),Q(2,
     $    NEXT))) GO TO 190
      IF (LFT.LE.0) GO TO 160
      IF (.NOT.LEFT(X0,Y0,Q(1,NR),Q(2,NR),Q(1,NEXT),Q(2,
     $    NEXT))) GO TO 190
c
c   replace nl->nr with next->n0.
c
      CALL GRSWAP(NEXT,N0,NL,NR,IADJ,IEND)
      IWK(1,IWC) = NEXT
      IWK(2,IWC) = N0
      GO TO 200
c
c   swap nl-nr for n0-next, shift columns iwf,...,iwc-1 to
c     the right, and store n0-next in the left portion of
c     iwk.
c
  160 CONTINUE
      CALL GRSWAP(NEXT,N0,NL,NR,IADJ,IEND)
      I      = IWC
  170 CONTINUE
      IF (I.EQ.IWF) GO TO 180
      IWK(1,I) = IWK(1,I-1)
      IWK(2,I) = IWK(2,I-1)
      I      = I - 1
      GO TO 170
*
  180 CONTINUE
      IWK(1,IWF) = N0
      IWK(2,IWF) = NEXT
      IWF    = IWF + 1
      GO TO 200
c
c   a swap is not possible.  set n0 to nl.
c
  190 CONTINUE
      N0     = NL
      X0     = Q(1,N0)
      Y0     = Q(2,N0)
      LFT    = -1
c
c   advance to the next arc.
c
  200 CONTINUE
      NL     = NEXT
      IWC    = IWC + 1
      GO TO 100
c
c   n2 is opposite nl->nr (iwc = iwl).
c
  210 CONTINUE
      IF (N0.EQ.N1) GO TO 240
      IF (LFT.LT.0) GO TO 220
c
c   n0 right n1->n2.  test for a possible swap.
c
      IF (.NOT.LEFT(X0,Y0,Q(1,NR),Q(2,NR),X2,Y2)) GO TO 90
c
c   swap nl-nr for n0-n2 and store n0-n2 in the right
c     portion of iwk.
c
      CALL GRSWAP(N2,N0,NL,NR,IADJ,IEND)
      IWK(1,IWL) = N0
      IWK(2,IWL) = N2
      IWL    = IWL - 1
      GO TO 90
c
c   n0 left n1->n2.  test for a possible swap.
c
  220 CONTINUE
      IF (.NOT.LEFT(Q(1,NL),Q(2,NL),X0,Y0,X2,Y2)) GO TO 90
c
c   swap nl-nr for n0-n2, shift columns iwf,...,iwl-1 to the
c     right, and store n0-n2 in the left portion of iwk.
c
      CALL GRSWAP(N2,N0,NL,NR,IADJ,IEND)
      I      = IWL
  230 CONTINUE
      IWK(1,I) = IWK(1,I-1)
      IWK(2,I) = IWK(2,I-1)
      I      = I - 1
      IF (I.GT.IWF) GO TO 230
      IWK(1,IWF) = N0
      IWK(2,IWF) = N2
      IWF    = IWF + 1
      GO TO 90
c
c iwf = iwc = iwl.  swap out the last arc for n1-n2 and
c   store zeros in iwk.
c
  240 CONTINUE
      CALL GRSWAP(N2,N1,NL,NR,IADJ,IEND)
      IWK(1,IWC) = 0
      IWK(2,IWC) = 0
      IF (IWC.EQ.1) GO TO 290
c
c optimization loops -- optimize the set of new arcs to the
c   left of in1->in2.  the loop is repeated until no swaps
c   are performed.
c
      IWCM1  = IWC - 1
  250 CONTINUE
      SWP    = .FALSE.
      DO 280 I = 1,IWCM1
         IO1    = IWK(1,I)
         IO2    = IWK(2,I)
c
c   set n1 to the neighbor of io1 which follows io2 and set
c     n2 to the neighbor of io1 which precedes io2.
c
         INDF   = 1
         IF (IO1.GT.1) INDF   = IEND(IO1-1) + 1
         INDL   = IEND(IO1)
         INDX   = INDL
         IF (IADJ(INDX).NE.IO2) GO TO 260
c
c   io2 is the last neighbor of io1.
c
         N1     = IADJ(INDF)
         N2     = IADJ(INDX-1)
         GO TO 270
c
c   io2 is not the last neighbor of io1.  loop through the
c     neighbors in reverse order.
c
  260    CONTINUE
         INDX   = INDX - 1
         IF (IADJ(INDX).NE.IO2) GO TO 260
         N1     = IADJ(INDX+1)
         IF (INDX.NE.INDF) N2     = IADJ(INDX-1)
         IF (INDX.EQ.INDF) N2     = IADJ(INDL)
c
c   test io1-io2 for a swap.
c
  270    CONTINUE
         IF (.NOT.GRSWPT(N1,N2,IO1,IO2,Q)) GO TO 280
         SWP    = .TRUE.
         CALL GRSWAP(N1,N2,IO1,IO2,IADJ,IEND)
         IWK(1,I) = N1
         IWK(2,I) = N2
  280 CONTINUE
      IF (SWP) GO TO 250
c
c test for termination.
c
  290 CONTINUE
      IF (IWC.EQ.IWEND) RETURN
      IWCP1  = IWC + 1
c
c optimize the set of new arcs to the right of in1->in2.
c
  300 CONTINUE
      SWP    = .FALSE.
      DO 330 I = IWCP1,IWEND
         IO1    = IWK(1,I)
         IO2    = IWK(2,I)
c
c   set n1 and n2 to the nodes opposite io1->io2 and
c     io2->io1, respectively.
c
         INDF   = 1
         IF (IO1.GT.1) INDF   = IEND(IO1-1) + 1
         INDL   = IEND(IO1)
         INDX   = INDL
         IF (IADJ(INDX).NE.IO2) GO TO 310
c
         N1     = IADJ(INDF)
         N2     = IADJ(INDX-1)
         GO TO 320
c
  310    CONTINUE
         INDX   = INDX - 1
         IF (IADJ(INDX).NE.IO2) GO TO 310
         N1     = IADJ(INDX+1)
         IF (INDX.NE.INDF) N2     = IADJ(INDX-1)
         IF (INDX.EQ.INDF) N2     = IADJ(INDL)
c
  320    CONTINUE
         IF (.NOT.GRSWPT(N1,N2,IO1,IO2,Q)) GO TO 330
         SWP    = .TRUE.
         CALL GRSWAP(N1,N2,IO1,IO2,IADJ,IEND)
         IWK(1,I) = N1
         IWK(2,I) = N2
  330 CONTINUE
      IF (SWP) GO TO 300
      RETURN
c
c in1 and in2 were adjacent on input.
c
  340 CONTINUE
      IER    = 0
      LWK    = 0
      RETURN
c
c parameter out of range
c
  350 CONTINUE
      IER    = 1
      RETURN
c
c insufficient space in iwk
c
  360 CONTINUE
      IER    = 2
      RETURN
c
c invalid triangulation data structure
c
  370 CONTINUE
      IER    = 3
      RETURN
*

      contains

      function left(XA,YA,XB,YB,XP,YP) result(erg)
      REAL,intent(in)::   XA,YA,XB,YB,XP,YP 
      LOGICAL erg
        
      erg = (XB-XA)* (YP-YA) >=  (XP-XA)* (YB-YA)
      end function left

      END SUBROUTINE GREDGE 
    

      SUBROUTINE GR3TBS(X,Y,Z,LV1,LV2,LV3,LV4,LL,LLV,KL,KV,KP,POI,N,M,R,
     $                  CLOSED,IFACES)
c     .. parameters ..
      REAL              ZWOPI
      INTEGER           MMAX
      PARAMETER         (ZWOPI=6.283185,MMAX=128)
c     ..
c     .. scalar arguments ..
      REAL              R
      INTEGER           IFACES,KL,KP,KV,M,N
      LOGICAL           CLOSED
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IQ,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IQ,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              POI(3,N),X(MNP),Y(MNP),Z(MNP)
      INTEGER           LL(2,MNL),LLV(2,MNL),LV1(MNF),LV2(MNF),LV3(MNF),
     $                  LV4(MNF)
c     ..
c     .. local scalars ..
      REAL              CP,FA,PHI,SP,SQE,SQF,SU1,SU2,SUN,SUZ
      INTEGER           I,IV,J,K,L,NN,IL,MOLD
      LOGICAL           CLOS,FIEQLA
c     ..
c     .. local arrays ..
      REAL              E(3),F(3),RE(3,MMAX),SE(3),SE1(3),SE2(3)
c     ..
c     .. intrinsic functions ..
      INTRINSIC         COS,SIN,SQRT
c     ..
      SAVE /GR3SAV/
      SAVE
c     .. executable statements ..

      MOLD = M
      M = ABS(M)
      IF (M.GT.MMAX) THEN
         WRITE (*,FMT=*) 'GR3TUB: M (5. ARGUMENT) IS TO LARGE'
         STOP
*
      END IF

      FIEQLA = .TRUE.
      DO 10 J = 1,3
         FIEQLA = FIEQLA .AND. (POI(J,1).EQ.POI(J,N))
   10 CONTINUE
      CLOS   = CLOSED .AND. .NOT. FIEQLA
      IF (CLOS) THEN
         L      = 2
         IF ( KP+(N+2)*M .GT. MNP ) THEN
           WRITE(*,*) 'GR3TUB: LAR must be greater.'
           STOP
         END IF
         DO 20 K = 1,M
            X(KP+K) = POI(1,1)
            Y(KP+K) = POI(2,1)
            Z(KP+K) = POI(3,1)
            X(KP+ (N+1)*M+K) = POI(1,N)
            Y(KP+ (N+1)*M+K) = POI(2,N)
            Z(KP+ (N+1)*M+K) = POI(3,N)
   20    CONTINUE
*
      ELSE
         L      = 1
      END IF

c---- ************ 2 zur startkurvenrichtung senkrechte einheitsvektoren

      E(1)   = POI(2,2) - POI(2,1) - POI(3,2) + POI(3,1)
      E(2)   = POI(3,2) - POI(3,1) - POI(1,2) + POI(1,1)
      E(3)   = POI(1,2) - POI(1,1) - POI(2,2) + POI(2,1)
      SQE    = SQRT(E(1)**2+E(2)**2+E(3)**2)
      E(1)   = E(1)/SQE
      E(2)   = E(2)/SQE
      E(3)   = E(3)/SQE
      F(1)   = (POI(2,2)-POI(2,1))*E(3) - (POI(3,2)-POI(3,1))*E(2)
      F(2)   = (POI(3,2)-POI(3,1))*E(1) - (POI(1,2)-POI(1,1))*E(3)
      F(3)   = (POI(1,2)-POI(1,1))*E(2) - (POI(2,2)-POI(2,1))*E(1)
      SQF    = SQRT(F(1)**2+F(2)**2+F(3)**2)
      F(1)   = F(1)/SQF
      F(2)   = F(2)/SQF
      F(3)   = F(3)/SQF

c---- ************ m zur startkurvenrichtung senkrechte r-vektoren
c---- ************ erster punktekranz

      DO 30 K = 1,M - 1
         PHI    = (K-1)*ZWOPI/ (M-1)
         CP     = COS(PHI)
         SP     = SIN(PHI)
         RE(1,K) = R* (CP*E(1)+SP*F(1))
         RE(2,K) = R* (CP*E(2)+SP*F(2))
         RE(3,K) = R* (CP*E(3)+SP*F(3))
         X(KP+ (L-1)*M+K) = RE(1,K) + POI(1,1)
         Y(KP+ (L-1)*M+K) = RE(2,K) + POI(2,1)
         Z(KP+ (L-1)*M+K) = RE(3,K) + POI(3,1)
   30 CONTINUE
      RE(1,M) = RE(1,1)
      RE(2,M) = RE(2,1)
      RE(3,M) = RE(3,1)
      X(KP+L*M) = X(KP+ (L-1)*M+1)
      Y(KP+L*M) = Y(KP+ (L-1)*M+1)
      Z(KP+L*M) = Z(KP+ (L-1)*M+1)

      SU1    = 0.
      DO 40 J = 1,3
         SE1(J) = POI(J,2) - POI(J,1)
         SU1    = SU1 + SE1(J)**2
   40 CONTINUE

c---- ************ schleife ueber alle inneren punkte

      DO 140 I = 2,N - 1
         L      = L + 1
c------- ************ berechnung der normalen der knickebene

         SU1    = SQRT(SU1)
         SU2    = 0.
         DO 50 J = 1,3
            SE2(J) = POI(J,I+1) - POI(J,I)
            SE1(J) = SE1(J)/SU1
            SU2    = SU2 + SE2(J)**2
   50    CONTINUE
         SU2    = SQRT(SU2)
         DO 60 J = 1,3
            SE2(J) = SE2(J)/SU2
            SE(J)  = SE1(J) + SE2(J)
   60    CONTINUE

c------- ************ schleife ueber den punktekranz am i-knick

         IF ( KP+L*M .GT. MNP ) THEN
           WRITE(*,*) 'GR3TUB: LAR must be greater.'
           STOP
         END IF
         DO 80 K = 1,M - 1
            SUZ    = 0.
            SUN    = 0.
            DO 70 J = 1,3
               SUZ    = SUZ + RE(J,K)*SE(J)
               SUN    = SUN + SE(J)*SE1(J)
   70       CONTINUE
            FA     = -SUZ/SUN
            X(KP+ (L-1)*M+K) = POI(1,I) + RE(1,K) + SE1(1)*FA
            Y(KP+ (L-1)*M+K) = POI(2,I) + RE(2,K) + SE1(2)*FA
            Z(KP+ (L-1)*M+K) = POI(3,I) + RE(3,K) + SE1(3)*FA
   80    CONTINUE
         X(KP+L*M) = X(KP+ (L-1)*M+1)
         Y(KP+L*M) = Y(KP+ (L-1)*M+1)
         Z(KP+L*M) = Z(KP+ (L-1)*M+1)

c------- ************ 2 neue zur kurvenrichtung senkrechte einheitsvekt.

         SU1    = 0.
         SUZ    = 0.
         SE2(1) = X(KP+ (L-1)*M+1) - POI(1,I)
         SE2(2) = Y(KP+ (L-1)*M+1) - POI(2,I)
         SE2(3) = Z(KP+ (L-1)*M+1) - POI(3,I)
         DO 90 J = 1,3
            SE1(J) = POI(J,I+1) - POI(J,I)
            SU1    = SU1 + SE1(J)**2
            SUZ    = SUZ + SE1(J)*SE2(J)
   90    CONTINUE
         FA     = SUZ/SU1
         SUN    = 0.
         DO 100 J = 1,3
            E(J)   = SE2(J) - FA*SE1(J)
            SUN    = SUN + E(J)**2
  100    CONTINUE
         SUN    = SQRT(SUN)
         E(1)   = E(1)/SUN
         E(2)   = E(2)/SUN
         E(3)   = E(3)/SUN
         F(1)   = SE1(2)*E(3) - SE1(3)*E(2)
         F(2)   = SE1(3)*E(1) - SE1(1)*E(3)
         F(3)   = SE1(1)*E(2) - SE1(2)*E(1)
         SQF    = SQRT(F(1)**2+F(2)**2+F(3)**2)
         F(1)   = F(1)/SQF
         F(2)   = F(2)/SQF
         F(3)   = F(3)/SQF

c---- ************ m zur neuen kurvenrichtung senkrechte r-vektoren

         DO 120 K = 1,M - 1
            PHI    = (K-1)*ZWOPI/ (M-1)
            CP     = COS(PHI)
            SP     = SIN(PHI)
            DO 110 J = 1,3
               RE(J,K) = R* (CP*E(J)+SP*F(J))
  110       CONTINUE
  120    CONTINUE
         DO 130 J = 1,3
            RE(J,M) = RE(J,1)
  130    CONTINUE
  140 CONTINUE

c---- ******************* letzter punktekranz

      IF (FIEQLA) THEN

         IF ( KP+N*M .GT. MNP ) THEN
           WRITE(*,*) 'GR3TUB: LAR must be greater.'
           STOP
         END IF

c------- ************ berechnung der normalen der knickebene

         SU1    = SQRT(SU1)
         SU2    = 0.
         DO 150 J = 1,3
            SE2(J) = POI(J,2) - POI(J,1)
            SE1(J) = SE1(J)/SU1
            SU2    = SU2 + SE2(J)**2
  150    CONTINUE
         SU2    = SQRT(SU2)
         DO 160 J = 1,3
            SE2(J) = SE2(J)/SU2
            SE(J)  = SE1(J) + SE2(J)
  160    CONTINUE

c------- ************ schleife ueber den punktekranz am i-knick

         DO 180 K = 1,M - 1
            SUZ    = 0.
            SUN    = 0.
            DO 170 J = 1,3
               SUZ    = SUZ + RE(J,K)*SE(J)
               SUN    = SUN + SE(J)*SE1(J)
  170       CONTINUE
            FA     = -SUZ/SUN
            X(KP+ (N-1)*M+K) = POI(1,N) + RE(1,K) + SE1(1)*FA
            Y(KP+ (N-1)*M+K) = POI(2,N) + RE(2,K) + SE1(2)*FA
            Z(KP+ (N-1)*M+K) = POI(3,N) + RE(3,K) + SE1(3)*FA
            X(KP+K) = X(KP+ (N-1)*M+K)
            Y(KP+K) = Y(KP+ (N-1)*M+K)
            Z(KP+K) = Z(KP+ (N-1)*M+K)
  180    CONTINUE
         X(KP+N*M) = X((N-1)*M+1)
         Y(KP+N*M) = Y((N-1)*M+1)
         Z(KP+N*M) = Z((N-1)*M+1)
         X(KP+M) = X(KP+1)
         Y(KP+M) = Y(KP+1)
         Z(KP+M) = Z(KP+1)
*
      ELSE
         L      = L + 1
         IF ( KP+L*M .GT. MNP ) THEN
           WRITE(*,*) 'GR3TUB: LAR must be greater.'
           STOP
         END IF
         DO 190 K = 1,M - 1
            X(KP+ (L-1)*M+K) = RE(1,K) + POI(1,N)
            Y(KP+ (L-1)*M+K) = RE(2,K) + POI(2,N)
            Z(KP+ (L-1)*M+K) = RE(3,K) + POI(3,N)
  190    CONTINUE
         X(KP+L*M) = X(KP+ (L-1)*M+1)
         Y(KP+L*M) = Y(KP+ (L-1)*M+1)
         Z(KP+L*M) = Z(KP+ (L-1)*M+1)
      END IF
*
      NN     = N
      IF (CLOS) NN     = N + 2
c---- erzeugen der faces ----------------------------------------------
      IV     = KV
      IF (IFACES.NE.0) THEN
         DO 210 I = 1,NN - 1
           IF (MOLD .GT. 0) THEN
            DO 200 J = 1,M - 1
               IV     = IV + 1
               LV1(IV) = KP + (I-1)*M + J
               LV2(IV) = KP + (I-1)*M + J + 1
               LV3(IV) = KP + I*M + J + 1
               LV4(IV) = KP + I*M + J
  200       CONTINUE
           ELSE
            DO 201 J = 1,M - 1
               IV     = IV + 1
               LV4(IV) = KP + (I-1)*M + J
               LV3(IV) = KP + (I-1)*M + J + 1
               LV2(IV) = KP + I*M + J + 1
               LV1(IV) = KP + I*M + J
  201       CONTINUE
           END IF
  210    CONTINUE
      END IF
c---- erzeugen der linien
      IL     = KL
      DO 230 J = 2,M - 1
         DO 220 I = 1,NN - 1
            IL     = IL + 1
            LL(1,IL) = KP + (I-1)*M + J
            LL(2,IL) = KP + I*M + J
            IF (IFACES.NE.0) THEN
               LLV(1,IL) = KV + (I-1)* (M-1) + J - 1
               LLV(2,IL) = KV + (I-1)* (M-1) + J
*
            ELSE
               LLV(1,IL) = 0
               LLV(2,IL) = 0
            END IF
*
  220    CONTINUE
  230 CONTINUE
      DO 240 I = 1,NN - 1
         IL     = IL + 1
         LL(1,IL) = KP + (I-1)*M + 1
         LL(2,IL) = KP + I*M + 1
         IF (IFACES.NE.0) THEN
            LLV(1,IL) = KV + I* (M-1)
            LLV(2,IL) = KV + (I-1)* (M-1) + 1
*
         ELSE
            LLV(1,IL) = 0
            LLV(2,IL) = 0
         END IF
*
  240 CONTINUE
      DO 260 I = 2,NN - 1
         DO 250 J = 1,M - 1
            IL     = IL + 1
            LL(1,IL) = KP + (I-1)*M + J
            LL(2,IL) = KP + (I-1)*M + J + 1
            IF (IFACES.NE.0) THEN
               LLV(1,IL) = KV + (I-1)* (M-1) + J
               LLV(2,IL) = KV + (I-2)* (M-1) + J
*
            ELSE
               LLV(1,IL) = 0
               LLV(2,IL) = 0
            END IF
*
  250    CONTINUE
  260 CONTINUE
      DO 270 J = 1,M - 1
         IL     = IL + 1
         LL(1,IL) = KP + J
         LL(2,IL) = KP + J + 1
         IF (IFACES.NE.0) THEN
            LLV(1,IL) = KV + J
            LLV(2,IL) = 0
*
         ELSE
            LLV(1,IL) = 0
            LLV(2,IL) = 0
         END IF
*
  270 CONTINUE
      DO 280 J = 1,M - 1
         IL     = IL + 1
         LL(1,IL) = KP + (NN-1)*M + J
         LL(2,IL) = KP + (NN-1)*M + J + 1
         IF (IFACES.NE.0) THEN
            LLV(1,IL) = KV + (NN-2)* (M-1) + J
            LLV(2,IL) = 0
*
         ELSE
            LLV(1,IL) = 0
            LLV(2,IL) = 0
         END IF
*
  280 CONTINUE
      KP     = KP + M*NN
      KL     = IL
      KV     = IV
      M = MOLD
      END
      SUBROUTINE GR3ARR(AR,IER,XYZ,RICLAN,N,KIND,RMAX,IFACES,JCO)
c---- kind=+ punkt am fuss, kind=0 in der mitte, kind=- an der spitze
c     .. scalar arguments ..
      REAL              AR,RMAX
      INTEGER           IER,IFACES,JCO,KIND,N
c     ..
c     .. array arguments ..
      REAL              RICLAN(3,N),XYZ(3,N)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3ARW,GR3SUB
c     ..
c     .. executable statements ..
      CALL GR3SUB(AR,IER,GR3ARW,JCO,XYZ(1,1),RICLAN(1,1),N,KIND,RMAX,
     $            IFACES,D1,D2,D3,D4)
      END
      SUBROUTINE GR3ARW(X,Y,Z,LV1,LV2,LV3,LV4,LL,LLV,KL,KV,KP,XYZ,RIC,N,
     $                  KIND,RMAX,IFACES)
c     .. parameters ..
      REAL              A,B,ZWOPI
      PARAMETER         (A=0.1,B=0.5,ZWOPI=6.283185)
c     ..
c     .. scalar arguments ..
      REAL              RMAX
      INTEGER           IFACES,KIND,KL,KP,KV,N
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IQ,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IQ,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              RIC(3,N),X(MNP),XYZ(3,N),Y(MNP),Z(MNP)
      INTEGER           LL(2,MNL),LLV(2,MNL),LV1(MNF),LV2(MNF),LV3(MNF),
     $                  LV4(MNF)
c     ..
c     .. local scalars ..
      REAL              AX,CP,FA1,FA2,FAC,PHI,R,SP,SQE,SQF,WU
      INTEGER           I,IP,IV,J,K,KIN,L,IL
c     ..
c     .. local arrays ..
      REAL              E(3),F(3),POI(3),RE(3),RICNOR(3)
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,COS,MOD,SIN,SQRT
c     ..
c     .. executable statements ..

      IL     = KL
      DO 90 L = 1,N
         IP     = KP
         IV     = KV
         WU     = SQRT(RIC(1,L)**2+RIC(2,L)**2+RIC(3,L)**2)
         IF (WU.LT.1E-16) THEN
            WRITE (*,FMT=*) 'GR3ARW: no direction given'
            STOP
*
         END IF
*
         RICNOR(1) = RIC(1,L)/WU
         RICNOR(2) = RIC(2,L)/WU
         RICNOR(3) = RIC(3,L)/WU
         AX     = A*WU
         R      = RMAX*AX/ (AX+RMAX)
         KIN    = KIND
         IF (KIN.NE.0) KIN    = KIN/ABS(KIN)
         FAC    = (KIN-1)*.5
         POI(1) = XYZ(1,L) + WU*RICNOR(1)*FAC
         POI(2) = XYZ(2,L) + WU*RICNOR(2)*FAC
         POI(3) = XYZ(3,L) + WU*RICNOR(3)*FAC
         IP     = IP + 1
         X(IP)  = POI(1)
         Y(IP)  = POI(2)
         Z(IP)  = POI(3)

c---- ************ 2 zum richtungsvektor senkrechte einheitsvektoren

         E(1)   = RICNOR(2) - RICNOR(3)
         E(2)   = RICNOR(3) - RICNOR(1)
         E(3)   = RICNOR(1) - RICNOR(2)
         SQE    = SQRT(E(1)**2+E(2)**2+E(3)**2)
         E(1)   = E(1)/SQE
         E(2)   = E(2)/SQE
         E(3)   = E(3)/SQE
         F(1)   = RICNOR(2)*E(3) - RICNOR(3)*E(2)
         F(2)   = RICNOR(3)*E(1) - RICNOR(1)*E(3)
         F(3)   = RICNOR(1)*E(2) - RICNOR(2)*E(1)
         SQF    = SQRT(F(1)**2+F(2)**2+F(3)**2)
         F(1)   = F(1)/SQF
         F(2)   = F(2)/SQF
         F(3)   = F(3)/SQF

c---- ************ 3 zum richtungsvektor senkrechte r-vektoren
c---- ************ erster punktekranz

         FA1    = WU - R/B*1.5
         FA2    = WU - R/B*2.0
         DO 10 K = 1,3
            PHI    = (K-1)*ZWOPI/3
            CP     = COS(PHI)
            SP     = SIN(PHI)
            RE(1)  = R* (CP*E(1)+SP*F(1))
            RE(2)  = R* (CP*E(2)+SP*F(2))
            RE(3)  = R* (CP*E(3)+SP*F(3))
            IP     = IP + 1
            X(IP)  = RE(1) + POI(1)
            Y(IP)  = RE(2) + POI(2)
            Z(IP)  = RE(3) + POI(3)
            IP     = IP + 1
            X(IP)  = X(IP-1) + RICNOR(1)*FA1
            Y(IP)  = Y(IP-1) + RICNOR(2)*FA1
            Z(IP)  = Z(IP-1) + RICNOR(3)*FA1
            IP     = IP + 1
            X(IP)  = X(IP-2) + RICNOR(1)*FA2 + RE(1)
            Y(IP)  = Y(IP-2) + RICNOR(2)*FA2 + RE(2)
            Z(IP)  = Z(IP-2) + RICNOR(3)*FA2 + RE(3)
   10    CONTINUE
         IP     = IP + 1
         X(IP)  = POI(1) + RICNOR(1)*WU
         Y(IP)  = POI(2) + RICNOR(2)*WU
         Z(IP)  = POI(3) + RICNOR(3)*WU
c---- die flaechen
         IF (IFACES.NE.0) THEN
            DO 20 I = 1,3
               IV     = IV + 1
               LV1(IV) = KP + 1
               LV2(IV) = KP + 2 + MOD(I,3)*3
               LV3(IV) = KP + 2 + (I-1)*3
               LV4(IV) = LV3(IV)
   20       CONTINUE
            DO 30 I = 1,3
               IV     = IV + 1
               LV1(IV) = KP + 2 + MOD(I,3)*3
               LV2(IV) = KP + 3 + MOD(I,3)*3
               LV3(IV) = KP + 3 + (I-1)*3
               LV4(IV) = KP + 2 + (I-1)*3
   30       CONTINUE
            DO 40 I = 1,3
               IV     = IV + 1
               LV1(IV) = KP + 3 + MOD(I,3)*3
               LV2(IV) = KP + 4 + MOD(I,3)*3
               LV3(IV) = KP + 4 + (I-1)*3
               LV4(IV) = KP + 3 + (I-1)*3
   40       CONTINUE
            DO 50 I = 1,3
               IV     = IV + 1
               LV1(IV) = KP + 11
               LV2(IV) = KP + 4 + (I-1)*3
               LV3(IV) = KP + 4 + MOD(I,3)*3
               LV4(IV) = LV3(IV)
   50       CONTINUE
         END IF
c---- die linien
         DO 60 I = 1,3
            IL     = IL + 1
            LL(1,IL) = KP + 1
            LL(2,IL) = KP + 2 + (I-1)*3
            IF (IFACES.NE.0) THEN
               LLV(1,IL) = KV + I
               LLV(2,IL) = KV + MOD(I+1,3) + 1
*
            ELSE
               LLV(1,IL) = 0
               LLV(2,IL) = 0
            END IF
*
            IL     = IL + 1
            LL(1,IL) = KP + 2 + (I-1)*3
            LL(2,IL) = KP + 3 + (I-1)*3
            IF (IFACES.NE.0) THEN
               LLV(1,IL) = KV + 3 + I
               LLV(2,IL) = KV + MOD(I+1,3) + 4
*
            ELSE
               LLV(1,IL) = 0
               LLV(2,IL) = 0
            END IF
*
            IL     = IL + 1
            LL(1,IL) = KP + 3 + (I-1)*3
            LL(2,IL) = KP + 4 + (I-1)*3
            IF (IFACES.NE.0) THEN
               LLV(1,IL) = KV + 6 + I
               LLV(2,IL) = KV + MOD(I+1,3) + 7
*
            ELSE
               LLV(1,IL) = 0
               LLV(2,IL) = 0
            END IF
*
            IL     = IL + 1
            LL(1,IL) = KP + 4 + (I-1)*3
            LL(2,IL) = KP + 11
            IF (IFACES.NE.0) THEN
               LLV(1,IL) = KV + 9 + I
               LLV(2,IL) = KV + MOD(I+1,3) + 10
*
            ELSE
               LLV(1,IL) = 0
               LLV(2,IL) = 0
            END IF
*
   60    CONTINUE
         DO 80 J = 1,3
            DO 70 I = 1,3
               IL     = IL + 1
               LL(1,IL) = KP + (I-1)*3 + 1 + J
               LL(2,IL) = KP + MOD(I,3)*3 + 1 + J
               IF (IFACES.NE.0) THEN
                  LLV(1,IL) = KV + 3* (J-1) + I
                  LLV(2,IL) = KV + 3*J + I
*
               ELSE
                  LLV(1,IL) = 0
                  LLV(2,IL) = 0
               END IF
*
   70       CONTINUE
   80    CONTINUE
         KP     = IP
         KV     = IV
   90 CONTINUE
      KL     = IL
      END
      SUBROUTINE GR3OBT(PS,KEN,LL,P,M,N,L,IPA,IIA,ILA,ILMAX)

c     .. scalar arguments ..
      INTEGER           IIA,ILA,ILMAX,IPA,M,N
c     ..
c     .. array arguments ..
      REAL              P(3,M,N),PS(3,ILMAX)
      INTEGER           KEN(ILMAX),L(M,N),LL(ILMAX)
c     ..
c     .. local scalars ..
      INTEGER           I,II,IL,IP,J
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      IP     = IPA
      II     = IIA
      IL     = ILA
      DO 20 I = 1,N
         II     = II + 1
         DO 10 J = 1,M
            IP     = IP + 1
            PS(1,IP) = P(1,J,I)
            PS(2,IP) = P(2,J,I)
            PS(3,IP) = P(3,J,I)
            KEN(IP) = J + 8*II
            IF (L(J,I).NE.0) THEN
               IL     = IL + 1
               IF (IL.GT.ILMAX) THEN
                  WRITE (*,FMT=*)
     $              'GR3OBT: zu viele Linien, AR zu klein!'
                  STOP
*
               END IF
*
               LL(IL) = J + 8*II
            END IF
*
   10    CONTINUE
         IF (M.EQ.3) THEN
            IP     = IP + 1
            PS(1,IP) = PS(1,IP-3)
            PS(2,IP) = PS(2,IP-3)
            PS(3,IP) = PS(3,IP-3)
            KEN(IP) = 4 + 8*II
         END IF
*
   20 CONTINUE
      IIA    = II
      IPA    = IP
      ILA    = IL
      END
      SUBROUTINE GR3OBJ(X,Y,Z,LV1,LV2,LV3,LV4,LL,LLV,KL,KV,KP,PS,KEN,
     $                  LIN,IND,ISO,IPA,IIA,ILA,IFACES)
c----- ipa : anzahl punkte
c----- iia : anzahl vierecke
c----- ila : anzahl linien

c     .. parameters ..
      REAL              EPS
      PARAMETER         (EPS=1E-4)
c     ..
c     .. scalar arguments ..
      INTEGER           IFACES,IIA,ILA,IPA,KL,KP,KV
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              PS(3,IPA),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IND(IPA),ISO(2,ILA),KEN(IPA),LIN(ILA),LL(2,MNL),
     $                  LLV(2,MNL),LV1(MNF),LV2(MNF),LV3(MNF),LV4(MNF)
c     ..
c     .. local scalars ..
      REAL              PX,PY,PZ,UX,UY,UZ
      INTEGER           I,IA1,IA2,IP,IP1,IP2,IPH,IZ,IZI,IZV,JL,LLV1,
     $                  LLV2,MAXIN,MININ
      LOGICAL           LOX,LOY,LOZ
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3QS1,GR3QS2
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MAX,MIN
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      CALL GR3QS1(PS,IND,IPA)
      IP     = KP
      PX     = (PS(1,IND(1))+EPS)*.99
      PY     = (PS(2,IND(1))+EPS)*.99
      PZ     = (PS(3,IND(1))+EPS)*.99
      DO 10 I = 1,IPA
         UX     = PS(1,IND(I))
         UY     = PS(2,IND(I))
         UZ     = PS(3,IND(I))
         LOX    = ABS(UX-PX) .GT. EPS* (ABS(UX)+ABS(PX))
         LOY    = ABS(UY-PY) .GT. EPS* (ABS(UY)+ABS(PY))
         LOZ    = ABS(UZ-PZ) .GT. EPS* (ABS(UZ)+ABS(PZ))
         IF (LOX .OR. LOY .OR. LOZ) THEN
            PX     = PS(1,IND(I))
            PY     = PS(2,IND(I))
            PZ     = PS(3,IND(I))
            IP     = IP + 1
            X(IP)  = PX
            Y(IP)  = PY
            Z(IP)  = PZ
         END IF
*
         IZ     = KEN(IND(I))
         IZV    = IZ/8
         IZI    = IZ - IZV*8
         IF (IZI.LT.3) THEN
            IF (IZI.EQ.1) THEN
               LV1(KV+IZV) = IP
*
            ELSE
               LV2(KV+IZV) = IP
            END IF
*
         ELSE
            IF (IZI.EQ.3) THEN
               LV3(KV+IZV) = IP
*
            ELSE
               LV4(KV+IZV) = IP
            END IF
*
         END IF
*
   10 CONTINUE
      DO 20 I = 1,ILA
         IZ     = LIN(I)
         IZV    = IZ/8
         IZI    = IZ - IZV*8
         IF (IZI.LT.3) THEN
            IF (IZI.EQ.1) THEN
               IP1    = LV1(KV+IZV)
               IP2    = LV2(KV+IZV)
*
            ELSE
               IP1    = LV2(KV+IZV)
               IP2    = LV3(KV+IZV)
            END IF
*
         ELSE
            IF (IZI.EQ.3) THEN
               IP1    = LV3(KV+IZV)
               IP2    = LV4(KV+IZV)
*
            ELSE
               IP1    = LV4(KV+IZV)
               IP2    = LV1(KV+IZV)
            END IF
*
         END IF
*
         IF (IP1.GT.IP2) THEN
            IPH    = IP1
            IP1    = IP2
            IP2    = IPH
            LIN(I) = -LIN(I)
         END IF
*
         ISO(1,I) = IP1
         ISO(2,I) = IP2
   20 CONTINUE
      CALL GR3QS2(ISO,IND,ILA)
      IA1    = -1
      IA2    = -1
      DO 30 I = 1,ILA
         IF (ISO(1,IND(I)).NE.IA1 .OR. ISO(2,IND(I)).NE.IA2) THEN
            IA1    = ISO(1,IND(I))
            IA2    = ISO(2,IND(I))
            KEN(IND(I)) = 0
*
         ELSE
            MININ  = MIN(IND(I),IND(I-1))
            MAXIN  = MAX(IND(I),IND(I-1))
            ISO(1,MININ) = 0
            KEN(MAXIN) = ABS(LIN(MININ))/8
         END IF
*
   30 CONTINUE
      JL     = KL
      DO 40 I = 1,ILA
         IF (ISO(1,I).NE.0) THEN
            JL     = JL + 1
            IZ     = LIN(I)
            IF (IZ.GT.0) THEN
               LL(1,JL) = ISO(1,I)
               LL(2,JL) = ISO(2,I)
*
            ELSE
               LL(1,JL) = ISO(2,I)
               LL(2,JL) = ISO(1,I)
            END IF
*
            IZV    = ABS(IZ)/8
            IF (IFACES.NE.0) THEN
               LLV(1,JL) = IZV + KV
               IF (KEN(I).NE.0) THEN
                  LLV(2,JL) = KEN(I) + KV
*
               ELSE
                  LLV(2,JL) = 0
               END IF
*
            ELSE
               LLV(1,JL) = 0
               LLV(2,JL) = 0
            END IF
*
            LLV1   = LLV(1,JL)
            LLV2   = LLV(2,JL)
         END IF
*
   40 CONTINUE
      KP     = IP
      IF (IFACES.NE.0) KV     = KV + IIA
      KL     = JL
      END
      SUBROUTINE GR3QS1(A,IND,N)
c---- sortiere punkte ('a' enthaelt n-mal drei koordinaten)



c     .. parameters ..
      INTEGER           NS,NMAX,K,MI,MIH
      PARAMETER         (NS=22,NMAX=2**NS,K=4,MI=2**K-1,MIH=MI/2)
c     ..
c     .. scalar arguments ..
      INTEGER           N
c     ..
c     .. array arguments ..
      REAL              A(3,N)
      INTEGER           IND(N)
c     ..
c     .. local scalars ..
      REAL              AHX,AHY,AHZ,GX,GY,GZ,ZX1,ZX2,ZY1,ZY2,ZZ1,ZZ2
      INTEGER           I,IH,IMI,INI,IPL,J,O,P,R1,R2,S
c     ..
c     .. local arrays ..
      INTEGER           STACO(NS),STACP(NS)
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      P      = N
      IF (P.GT.NMAX) THEN
         WRITE (*,FMT=*) 'GR3QS1: zu viele Elemente'
         STOP
*
      END IF
*
      S      = 0
      O      = 1
      DO 10 I = 1,N
         IND(I) = I
   10 CONTINUE
      GO TO 30

   20 CONTINUE
      O      = STACO(S)
      P      = STACP(S)
      S      = S - 1
   30 CONTINUE
      IF (P-O.LE.MI) THEN
c------- sortiere mit einfacherem verfahren (ersetzen)
         DO 50 I = O + 1,P
            AHX    = A(1,IND(I))
            AHY    = A(2,IND(I))
            AHZ    = A(3,IND(I))
            J      = I
            INI    = IND(I)
   40       CONTINUE
            IF (J.GT.1) THEN
               IF (AHX.LT.A(1,IND(J-1)) .OR.
     $             AHX.EQ.A(1,IND(J-1)) .AND. (AHY.LT.A(2,
     $             IND(J-1)).OR.AHY.EQ.A(2,IND(J-1)).AND.AHZ.LT.A(3,
     $             IND(J-1)))) THEN
                  IND(J) = IND(J-1)
                  J      = J - 1
                  GO TO 40
*
               ELSE
                  IND(J) = INI
               END IF
*
            ELSE
               IND(J) = INI
            END IF
*
   50    CONTINUE
         IF (S.GT.0) GO TO 20
         GO TO 80
*
      ELSE
c------- sortiere mit quicksort
         I      = (O+P)/2
         IMI    = IND(I-MIH)
         IPL    = IND(I+MIH)
         IF (A(1,IPL).GT.A(1,IMI) .OR. A(1,IPL).EQ.A(1,IMI) .AND.
     $       (A(2,IPL).GT.A(2,IMI).OR.A(2,IPL).EQ.A(2,IMI).AND.A(3,
     $       IPL).GT.A(3,IMI))) THEN
            ZX1    = A(1,IMI)
            ZY1    = A(2,IMI)
            ZZ1    = A(3,IMI)
            ZX2    = A(1,IPL)
            ZY2    = A(2,IPL)
            ZZ2    = A(3,IPL)
*
         ELSE
            ZX2    = A(1,IMI)
            ZY2    = A(2,IMI)
            ZZ2    = A(3,IMI)
            ZX1    = A(1,IPL)
            ZY1    = A(2,IPL)
            ZZ1    = A(3,IPL)
         END IF
*
         GX     = A(1,IND(I))
         GY     = A(2,IND(I))
         GZ     = A(3,IND(I))
         IF (GX.LT.ZX1 .OR. GX.EQ.ZX1 .AND.
     $       (GY.LT.ZY1.OR.GY.EQ.ZY1.AND.GZ.LT.ZZ1)) THEN
            GX     = ZX1
            GY     = ZY1
            GZ     = ZZ1
*
         ELSE IF (GX.GT.ZX2 .OR. GX.EQ.ZX2 .AND.
     $            (GY.GT.ZY2.OR.GY.EQ.ZY2.AND.GZ.GT.ZZ2)) THEN
            GX     = ZX2
            GY     = ZY2
            GZ     = ZZ2
         END IF
*
         I      = O
         J      = P
   60    CONTINUE
         IF (A(1,IND(I)).LT.GX .OR. A(1,IND(I)).EQ.GX .AND.
     $       (A(2,IND(I)).LT.GY.OR.A(2,IND(I)).EQ.GY.AND.A(3,
     $       IND(I)).LT.GZ)) THEN
            I      = I + 1
            GO TO 60
*
         ELSE
   70       CONTINUE
            IF (A(1,IND(J)).GT.GX .OR. A(1,IND(J)).EQ.GX .AND.
     $          (A(2,IND(J)).GT.GY.OR.A(2,IND(J)).EQ.GY.AND.A(3,
     $          IND(J)).GT.GZ)) THEN
               J      = J - 1
               GO TO 70
*
            END IF
*
         END IF
*
         IF (J.GE.I) THEN
            IH     = IND(I)
            IND(I) = IND(J)
            IND(J) = IH
            I      = I + 1
            J      = J - 1
            IF (J.GE.I) GO TO 60
         END IF
c------- eintragen in den stack
         IF (J-O.GE.P-I) THEN
            R1     = S + 1
            R2     = S + 2
*
         ELSE
            R1     = S + 2
            R2     = S + 1
         END IF
*
         STACO(R1) = O
         STACO(R2) = I
         STACP(R1) = J
         STACP(R2) = P
         S      = S + 2
      END IF
*
      GO TO 20
*
   80 CONTINUE
      END
      SUBROUTINE GR3QS2(A,IND,N)
c---- sortiere linien ('a' enthaelt die indizes von je 2 punkten)



c     .. parameters ..
      INTEGER           NS,NMAX,K,MI,MIH
      PARAMETER         (NS=22,NMAX=2**NS,K=4,MI=2**K-1,MIH=MI/2)
c     ..
c     .. scalar arguments ..
      INTEGER           N
c     ..
c     .. array arguments ..
      INTEGER           A(2,N),IND(N)
c     ..
c     .. local scalars ..
      INTEGER           AHX,AHY,GX,GY,I,IH,IMI,INI,IPL,J,O,P,R1,R2,S,
     $                  ZX1,ZX2,ZY1,ZY2
c     ..
c     .. local arrays ..
      INTEGER           STACO(NS),STACP(NS)
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      P      = N
      IF (P.GT.NMAX) THEN
         WRITE (*,FMT=*) 'GR3QS2: zu viele Elemente'
         STOP
*
      END IF
*
      S      = 0
      O      = 1
      DO 10 I = 1,N
         IND(I) = I
   10 CONTINUE
      GO TO 30

   20 CONTINUE
      O      = STACO(S)
      P      = STACP(S)
      S      = S - 1
   30 CONTINUE
      IF (P-O.LE.MI) THEN
c------- sortiere mit einfacherem verfahren (ersetzen)
         DO 50 I = O + 1,P
            AHX    = A(1,IND(I))
            AHY    = A(2,IND(I))
            J      = I
            INI    = IND(I)
   40       CONTINUE
            IF (J.GT.1) THEN
               IF (AHX.LT.A(1,IND(J-1)) .OR.
     $             AHX.EQ.A(1,IND(J-1)) .AND. AHY.LT.A(2,IND(J-1))) THEN
                  IND(J) = IND(J-1)
                  J      = J - 1
                  GO TO 40
*
               ELSE
                  IND(J) = INI
               END IF
*
            ELSE
               IND(J) = INI
            END IF
*
   50    CONTINUE
         IF (S.GT.0) GO TO 20
         GO TO 80
*
      ELSE
c------- sortiere mit quicksort
         I      = (O+P)/2
         IMI    = IND(I-MIH)
         IPL    = IND(I+MIH)
         IF (A(1,IPL).GT.A(1,IMI) .OR. A(1,IPL).EQ.A(1,IMI) .AND.
     $       A(2,IPL).GT.A(2,IMI)) THEN
            ZX1    = A(1,IMI)
            ZY1    = A(2,IMI)
            ZX2    = A(1,IPL)
            ZY2    = A(2,IPL)
*
         ELSE
            ZX2    = A(1,IMI)
            ZY2    = A(2,IMI)
            ZX1    = A(1,IPL)
            ZY1    = A(2,IPL)
         END IF
*
         GX     = A(1,IND(I))
         GY     = A(2,IND(I))
         IF (GX.LT.ZX1 .OR. GX.EQ.ZX1 .AND. GY.LT.ZY1) THEN
            GX     = ZX1
            GY     = ZY1
*
         ELSE IF (GX.GT.ZX2 .OR. GX.EQ.ZX2 .AND. GY.GT.ZY2) THEN
            GX     = ZX2
            GY     = ZY2
         END IF
*
         I      = O
         J      = P
   60    CONTINUE
         IF (A(1,IND(I)).LT.GX .OR. A(1,IND(I)).EQ.GX .AND.
     $       A(2,IND(I)).LT.GY) THEN
            I      = I + 1
            GO TO 60
*
         ELSE
   70       CONTINUE
            IF (A(1,IND(J)).GT.GX .OR. A(1,IND(J)).EQ.GX .AND.
     $          A(2,IND(J)).GT.GY) THEN
               J      = J - 1
               GO TO 70
*
            END IF
*
         END IF
*
         IF (J.GE.I) THEN
            IH     = IND(I)
            IND(I) = IND(J)
            IND(J) = IH
            I      = I + 1
            J      = J - 1
            IF (J.GE.I) GO TO 60
         END IF
c------- eintragen in den stack
         IF (J-O.GE.P-I) THEN
            R1     = S + 1
            R2     = S + 2
*
         ELSE
            R1     = S + 2
            R2     = S + 1
         END IF
*
         STACO(R1) = O
         STACO(R2) = I
         STACP(R1) = J
         STACP(R2) = P
         S      = S + 2
      END IF
*
      GO TO 20
*
   80 CONTINUE
      END
      SUBROUTINE GR3BLK(AR,IER,ZENT,BDH,N,IFACES,JCO)

c     .. scalar arguments ..
      REAL              AR
      INTEGER           IER,IFACES,JCO,N
c     ..
c     .. array arguments ..
      REAL              BDH(3,N),ZENT(3,N)
c     ..
c     .. local scalars ..
      INTEGER           I,J,K,L
c     ..
c     .. local arrays ..
      REAL              P(3,4),PNT(3,2,2,2)
      INTEGER           IP(3,4,6),LIN(4)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3OBB,GR3OBE,GR3OBP
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. data statements ..

      DATA              IP/1,1,1,2,1,1,2,2,1,1,2,1,1,1,1,1,1,2,2,1,2,2,
     $                  1,1,2,1,1,2,1,2,2,2,2,2,2,1,1,2,1,2,2,1,2,2,2,1,
     $                  2,2,1,1,1,1,2,1,1,2,2,1,1,2,1,1,2,1,2,2,2,2,2,2,
     $                  1,2/
      DATA              LIN/1,1,1,1/
c     ..
c     .. executable statements ..

      CALL GR3OBB

      DO 70 L = 1,N
         DO 30 K = 1,2
            DO 20 J = 1,2
               DO 10 I = 1,2
                  PNT(1,I,J,K) = ZENT(1,L) + BDH(1,L)* ((2*I-3)*0.5)
                  PNT(2,I,J,K) = ZENT(2,L) + BDH(2,L)* ((2*J-3)*0.5)
                  PNT(3,I,J,K) = ZENT(3,L) + BDH(3,L)* (K-1)
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE

         DO 60 K = 1,6
            DO 50 J = 1,4
               DO 40 I = 1,3
                  P(I,J) = PNT(I,IP(1,5-J,K),IP(2,5-J,K),IP(3,5-J,K))
   40          CONTINUE
   50       CONTINUE
            CALL GR3OBP(AR,IER,P,4,1,LIN)
   60    CONTINUE
   70 CONTINUE
      CALL GR3OBE(AR,IER,IFACES,JCO)

      END
      SUBROUTINE GR3IMP(AR,IE,N1,N2,F,NX,IX,X,NY,IY,Y,NZ,IZ,Z,CF,IF,JCO)

c     .. scalar arguments ..
      REAL              CF,F,X,Y,Z
      INTEGER           IE,IF,IX,IY,IZ,JCO,N1,N2,NX,NY,NZ
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP)
c     ..
c     .. local scalars ..
      INTEGER           K
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3IMF
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      K      = 46*MNP - 6*NX*NY + 1
      IF (K.LE.34.375*MNP) THEN
         WRITE (*,FMT=*) 'GR3IMP: Arbeitsspeicher zu klein'
         STOP 'GR3IMP'
*
      END IF
*
      CALL GR3IMF(AR,IE,N1,N2,F,NX,IX,X,NY,IY,Y,NZ,IZ,Z,CF,IF,JCO,AR(K))
      END
      SUBROUTINE GR3IMF(AR,IE,N1,N2,F,NX,IX,X,NY,IY,Y,NZ,IZ,Z,C,IF,JC,Q)

c     .. parameters ..
      REAL              BLIND1,BLIND2
      PARAMETER         (BLIND1=-75.7501E20,BLIND2=-75.7499E20)
c     ..
c     .. scalar arguments ..
      REAL              C
      INTEGER           IE,IF,IX,IY,IZ,JC,N1,N2,NX,NY,NZ
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP),F(N1,N2, (NZ-1)*IZ+1),Q(3,NX,NY,0:1),
     $                  X(N1),Y(N2),Z((NZ-1)*IZ+1)
c     ..
c     .. local scalars ..
      REAL              AF,C11,C12,C21,C22,C31,C32,F1,F2,FA
      INTEGER           I,I1,I2,IP,J,J1,J2,JN,K,K1,K2,KK,L,M
c     ..
c     .. local arrays ..
      REAL              W(2,0:11),XYZ(3,7)
      INTEGER           IPX(2,0:11),IPY(2,0:11),IPZ(2,0:11),LI1(3,0:11),
     $                  LI2(3,0:11),LIN(5)
      LOGICAL           NOTYET(0:11),OB0(0:11)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3OBB,GR3OBE,GR3OBP,GR3TRC
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. data statements ..

      DATA              IPX/0,1,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,
     $                  0,0/
      DATA              IPY/0,0,0,1,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,
     $                  0,1/
      DATA              IPZ/0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,1,
     $                  1,1/
      DATA              LI1/3,2,1,0,3,2,7,10,6,4,11,7,0,5,8,1,6,9,2,7,
     $                  10,3,4,11,4,0,5,5,1,6,11,8,9,8,9,10/
      DATA              LI2/5,8,4,6,9,5,1,0,3,2,1,0,11,7,3,8,4,0,9,5,1,
     $                  10,6,2,9,10,11,10,11,8,6,2,7,7,3,4/
      DATA              LIN/1,1,1,1,0/
c     ..
c     .. executable statements ..

      CALL GR3OBB

c------- transformation der eckpunkte eines kubus von den ursprueng-
c------- lichen koordinaten in die kartesichen koordinaten, in denen
c------- auch der hidden line/surface algorithmus arbeitet.
c------- es ist nicht gut, schon in den urspruenglichen koordinaten
c------- zu interpolieren und nur die interpolierten punkte zu trans-
c------- formieren. ineinanderliegende flaechen schneiden sich dann
c------- gelegentlich.

c---- zuweisen ebene 2 der transformationsmatrix:

      DO 20 J = 1,NY
         DO 10 I = 1,NX
            CALL GR3TRC(X((I-1)*IX+1),Y((J-1)*IY+1),Z(1),Q(1,I,J,1),
     $                  Q(2,I,J,1),Q(3,I,J,1))
   10    CONTINUE
   20 CONTINUE

      DO 110 K = 1,NZ - 1

c---- zuweisen ebene 2 => ebene1; neuberechnen der ebene 2 transform.

         DO 40 J = 1,NY
            DO 30 I = 1,NX
               Q(1,I,J,0) = Q(1,I,J,1)
               Q(2,I,J,0) = Q(2,I,J,1)
               Q(3,I,J,0) = Q(3,I,J,1)
               CALL GR3TRC(X((I-1)*IX+1),Y((J-1)*IY+1),Z(K*IZ+1),
     $                     Q(1,I,J,1),Q(2,I,J,1),Q(3,I,J,1))
   30       CONTINUE
   40    CONTINUE

         DO 100 J = 1,NY - 1
            DO 90 I = 1,NX - 1
               DO 50 L = 0,11
                  NOTYET(L) = .TRUE.
                  I1     = (I+IPX(1,L)-1)*IX + 1
                  I2     = (I+IPX(2,L)-1)*IX + 1
                  J1     = (J+IPY(1,L)-1)*IY + 1
                  J2     = (J+IPY(2,L)-1)*IY + 1
                  K1     = (K+IPZ(1,L)-1)*IZ + 1
                  K2     = (K+IPZ(2,L)-1)*IZ + 1
                  F1     = F(I1,J1,K1)
                  F2     = F(I2,J2,K2)
                  IF ((F1.LE.BLIND1.OR.F1.GE.BLIND2) .AND.
     $                (F2.LE.BLIND1.OR.F2.GE.BLIND2)) THEN
                     W(1,L) = F1 - C
                     W(2,L) = F2 - C
                     OB0(L) = (W(1,L).LT.0. .AND. W(2,L).GE.0.) .OR.
     $                        (W(2,L).LT.0. .AND. W(1,L).GE.0.)
*
                  ELSE
                     OB0(L) = .FALSE.
                  END IF
*
   50          CONTINUE
               DO 80 L = 0,9
                  IF (OB0(L) .AND. NOTYET(L)) THEN
                     M      = L
                     IP     = 0
   60                CONTINUE
                     IF (NOTYET(M)) THEN
                        IP     = IP + 1
                        FA     = -W(1,M)/ (W(2,M)-W(1,M))
                        AF     = 1. - FA
                        C11    = Q(1,I+IPX(1,M),J+IPY(1,M),IPZ(1,M))
                        C12    = Q(1,I+IPX(2,M),J+IPY(2,M),IPZ(2,M))
                        C21    = Q(2,I+IPX(1,M),J+IPY(1,M),IPZ(1,M))
                        C22    = Q(2,I+IPX(2,M),J+IPY(2,M),IPZ(2,M))
                        C31    = Q(3,I+IPX(1,M),J+IPY(1,M),IPZ(1,M))
                        C32    = Q(3,I+IPX(2,M),J+IPY(2,M),IPZ(2,M))
                        XYZ(1,IP) = AF*C11 + FA*C12
                        XYZ(2,IP) = AF*C21 + FA*C22
                        XYZ(3,IP) = AF*C31 + FA*C32
                        NOTYET(M) = .FALSE.

                        DO 70 KK = 1,3

                           IF (W(1,M).LT.0.) THEN
                              JN     = LI1(KK,M)
*
                           ELSE
                              JN     = LI2(KK,M)
                           END IF

                           IF (OB0(JN)) THEN
                              M      = JN
                              GO TO 60
*
                           END IF
*
   70                   CONTINUE
                     END IF
c                    ein stueck zu ende
                     IF (IP.EQ.3) THEN
                        CALL GR3OBP(AR,IE,XYZ,3,1,LIN)
*
                     ELSE IF (IP.EQ.4) THEN
                        CALL GR3OBP(AR,IE,XYZ,4,1,LIN)
*
                     ELSE IF (IP.EQ.5) THEN
                        CALL GR3OBP(AR,IE,XYZ,4,1,LIN(2))
                        XYZ(1,6) = XYZ(1,1)
                        XYZ(2,6) = XYZ(2,1)
                        XYZ(3,6) = XYZ(3,1)
                        CALL GR3OBP(AR,IE,XYZ(1,4),3,1,LIN(3))
*
                     ELSE IF (IP.EQ.6) THEN
                        CALL GR3OBP(AR,IE,XYZ,4,1,LIN(2))
                        XYZ(1,7) = XYZ(1,1)
                        XYZ(2,7) = XYZ(2,1)
                        XYZ(3,7) = XYZ(3,1)
                        CALL GR3OBP(AR,IE,XYZ(1,4),4,1,LIN(2))
                     END IF
*
                  END IF

   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      CALL GR3OBE(AR,IE,IF,JC)
      END
      SUBROUTINE GR3HHL(AR,IER,NX,IX,X,NY,IY,Y,NW,WE,N1,F,IG,IFA,JCO)


c     .. scalar arguments ..
      REAL              F,WE,X,Y
      INTEGER           IER,IFA,IG,IX,IY,JCO,N1,NW,NX,NY
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP)
c     ..
c     .. local scalars ..
      INTEGER           I0
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3ISL
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..
      I0     = 46*MNP - 6*NX*NY + 1
      IF (I0.LE.34.375*MNP) THEN
         WRITE (*,FMT=*) 'GR3HHL: Arbeitsspeicher zu klein.'
         STOP 'GR3HHL'
*
      END IF
*
      CALL GR3ISL(AR,IER,NX,IX,X,NY,IY,Y,NW,WE,N1,F,IG,IFA,JCO,AR(I0))
      END
      SUBROUTINE GR3ISL(AR,IER,NX,IX,X,NY,IY,Y,NW,WE,N1,F,IG,IFA,JCO,Q)
c
c     es hat sich durch test ergeben, dass es besser ist, wie hier
c     erst die koordinateneckpunkte zu transformieren und dann zu
c     interpolieren, als umgekehrt.


c     .. parameters ..
      REAL              BLIND1,BLIND2
      PARAMETER         (BLIND1=-75.7501E20,BLIND2=-75.7499E20)
c     ..
c     .. scalar arguments ..
      INTEGER           IER,IFA,IG,IX,IY,JCO,N1,NW,NX,NY
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP),F(N1, (NY-1)*IY+1),Q(3,NX,NY,0:1),
     $                  WE(NW),X(N1),Y((NY-1)*IY+1)
c     ..
c     .. local scalars ..
      REAL              AF,C11,C12,C21,C22,C31,C32,F1,F2,FA
      INTEGER           I,IP,J,JN,K,KK,L,M
c     ..
c     .. local arrays ..
      REAL              W(2,0:11),XYZ(3,7)
      INTEGER           IPX(2,0:11),IPY(2,0:11),IPZ(2,0:11),LI1(3,0:11),
     $                  LI2(3,0:11),LIN(4),MRK(7)
      LOGICAL           NOTYET(0:11),OB0(0:11)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3ISH,GR3OBB,GR3OBE,GR3OBP,GR3TRC
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. data statements ..

      DATA              IPX/0,1,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,
     $                  0,0/
      DATA              IPY/0,0,0,1,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,
     $                  0,1/
      DATA              IPZ/0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,1,
     $                  1,1/
      DATA              LI1/3,2,1,0,3,2,7,10,6,4,11,7,0,5,8,1,6,9,2,7,
     $                  10,3,4,11,4,0,5,5,1,6,11,8,9,8,9,10/
      DATA              LI2/5,8,4,6,9,5,1,0,3,2,1,0,11,7,3,8,4,0,9,5,1,
     $                  10,6,2,9,10,11,10,11,8,6,2,7,7,3,4/
c     ..
c     .. executable statements ..

      CALL GR3OBB

c---- zuweisen ebene 2

      DO 20 J = 1,NY
         DO 10 I = 1,NX
            CALL GR3TRC(X((I-1)*IX+1),Y((J-1)*IY+1),WE(1),Q(1,I,J,1),
     $                  Q(2,I,J,1),Q(3,I,J,1))
   10    CONTINUE
   20 CONTINUE

      DO 110 K = 1,NW - 1

c---- zuweisen ebene 2 => ebene1; neues berechnen der ebene 2
         DO 40 J = 1,NY
            DO 30 I = 1,NX
               Q(1,I,J,0) = Q(1,I,J,1)
               Q(2,I,J,0) = Q(2,I,J,1)
               Q(3,I,J,0) = Q(3,I,J,1)
               CALL GR3TRC(X((I-1)*IX+1),Y((J-1)*IY+1),WE(K+1),
     $                     Q(1,I,J,1),Q(2,I,J,1),Q(3,I,J,1))
   30       CONTINUE
   40    CONTINUE

         DO 100 J = 1,NY - 1
            DO 90 I = 1,NX - 1
               DO 50 L = 0,11
                  NOTYET(L) = .TRUE.
                  F1     = F((I+IPX(1,L)-1)*IX+1, (J+IPY(1,L)-1)*IY+1)
                  F2     = F((I+IPX(2,L)-1)*IX+1, (J+IPY(2,L)-1)*IY+1)
                  IF ((F1.LE.BLIND1.OR.F1.GE.BLIND2) .AND.
     $                (F2.LE.BLIND1.OR.F2.GE.BLIND2)) THEN
                     W(1,L) = F1 - WE(K+IPZ(1,L))
                     W(2,L) = F2 - WE(K+IPZ(2,L))
                     OB0(L) = (W(1,L).LT.0. .AND. W(2,L).GE.0.) .OR.
     $                        (W(2,L).LT.0. .AND. W(1,L).GE.0.)
*
                  ELSE
                     OB0(L) = .FALSE.
                  END IF
*
   50          CONTINUE
               DO 80 L = 0,9
                  IF (OB0(L) .AND. NOTYET(L)) THEN
                     M      = L
                     IP     = 0
   60                CONTINUE
                     IF (NOTYET(M)) THEN
                        IP     = IP + 1
                        FA     = -W(1,M)/ (W(2,M)-W(1,M))
                        AF     = 1. - FA
                        C11    = Q(1,I+IPX(1,M),J+IPY(1,M),IPZ(1,M))
                        C12    = Q(1,I+IPX(2,M),J+IPY(2,M),IPZ(2,M))
                        C21    = Q(2,I+IPX(1,M),J+IPY(1,M),IPZ(1,M))
                        C22    = Q(2,I+IPX(2,M),J+IPY(2,M),IPZ(2,M))
                        C31    = Q(3,I+IPX(1,M),J+IPY(1,M),IPZ(1,M))
                        C32    = Q(3,I+IPX(2,M),J+IPY(2,M),IPZ(2,M))
                        XYZ(1,IP) = AF*C11 + FA*C12
                        XYZ(2,IP) = AF*C21 + FA*C22
                        XYZ(3,IP) = AF*C31 + FA*C32
                        MRK(IP) = IPZ(1,M) + IPZ(2,M)
                        NOTYET(M) = .FALSE.
                        DO 70 KK = 1,3

                           IF (W(1,M).LT.0.) THEN
                              JN     = LI1(KK,M)
*
                           ELSE
                              JN     = LI2(KK,M)
                           END IF

                           IF (OB0(JN)) THEN
                              M      = JN
                              GO TO 60
*
                           END IF
*
   70                   CONTINUE
                     END IF
c                    ein stueck zu ende
                     IF (IP.EQ.3) THEN
                        MRK(4) = MRK(1)
                        CALL GR3ISH(MRK,IG,LIN,3)
                        CALL GR3OBP(AR,IER,XYZ,3,1,LIN)
*
                     ELSE IF (IP.EQ.4) THEN
                        MRK(5) = MRK(1)
                        CALL GR3ISH(MRK,IG,LIN,4)
                        CALL GR3OBP(AR,IER,XYZ,4,1,LIN)
*
                     ELSE IF (IP.EQ.5) THEN
                        CALL GR3ISH(MRK,IG,LIN,3)
                        LIN(4) = 0
                        CALL GR3OBP(AR,IER,XYZ,4,1,LIN)
                        XYZ(1,6) = XYZ(1,1)
                        XYZ(2,6) = XYZ(2,1)
                        XYZ(3,6) = XYZ(3,1)
                        MRK(6) = MRK(1)
                        CALL GR3ISH(MRK(4),IG,LIN,2)
                        LIN(3) = 0
                        CALL GR3OBP(AR,IER,XYZ(1,4),3,1,LIN)
*
                     ELSE IF (IP.EQ.6) THEN
                        CALL GR3ISH(MRK,IG,LIN,3)
                        LIN(4) = 0
                        CALL GR3OBP(AR,IER,XYZ,4,1,LIN)
                        XYZ(1,7) = XYZ(1,1)
                        XYZ(2,7) = XYZ(2,1)
                        XYZ(3,7) = XYZ(3,1)
                        MRK(7) = MRK(1)
                        CALL GR3ISH(MRK(4),IG,LIN,3)
                        CALL GR3OBP(AR,IER,XYZ(1,4),4,1,LIN)
                     END IF
*
                  END IF

   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      CALL GR3OBE(AR,IER,IFA,JCO)
      END
      SUBROUTINE GR3ISH(MRK,IG,LIN,N)

c     .. scalar arguments ..
      INTEGER           IG,N
c     ..
c     .. array arguments ..
      INTEGER           LIN(N),MRK(N+1)
c     ..
c     .. local scalars ..
      INTEGER           I
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      IF (IG.EQ.0) THEN
         DO 10 I = 1,N
            LIN(I) = 1
   10    CONTINUE
*
      ELSE IF (IG.LT.0) THEN
         DO 20 I = 1,N
            IF (MRK(I).EQ.1 .OR. MRK(I+1).EQ.1 .OR.
     $          MRK(I).NE.MRK(I+1)) THEN
               LIN(I) = 1
*
            ELSE
               LIN(I) = 0
            END IF
*
   20    CONTINUE
*
      ELSE
         DO 30 I = 1,N
            IF (MRK(I).NE.1 .AND. MRK(I+1).NE.1 .AND.
     $          MRK(I).EQ.MRK(I+1)) THEN
               LIN(I) = 1
*
            ELSE
               LIN(I) = 0
            END IF
*
   30    CONTINUE
      END IF

      END
      SUBROUTINE GR3MRK(AR,IER,ZENT,R,KIND,N,IFACES,JCO)


c     .. scalar arguments ..
      INTEGER           IER,IFACES,JCO,N
c     ..
c     .. common blocks ..
      REAL              WIN,ZENPRO
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              GR(3,3),RT(3,3)
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
c     ..
c     .. array arguments ..
      REAL              AR(46*MNP),R(N),ZENT(3,N)
      INTEGER           KIND(N)
c     ..
c     .. local scalars ..
      REAL              SUX,SUY,SUZ,X,Y,Z
      INTEGER           I,J,K,KI,L,NFACES
c     ..
c     .. local arrays ..
      REAL              DODA(3,20),HEXA(3,8),OKTA(3,6),P(3,6),Q(3,20),
     $                  TETRA(3,4)
      INTEGER           IDODA(5,12),IHEXA(4,6),IOKTA(3,8),ITETRA(3,4),
     $                  LIN(5),NN(3,20)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3OBB,GR3OBE,GR3OBP,GR3TRC
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MOD
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. data statements ..

c     data    itetra      /3,2,1, 3,1,4, 2,4,1, 2,3,4/
      DATA              DODA/0.4911235,0.3568221,-0.7946545,-0.1875925,
     $                  0.5773503,-0.7946545,-0.607062,0.0,-0.7946545,
     $                  -0.1875925,-0.5773503,-0.7946545,0.4911235,
     $                  -0.3568221,-0.7946545,0.7946545,0.5773503,
     $                  -0.1875925,-0.303531,0.9341724,-0.1875925,
     $                  -0.982247,0.0,-0.1875925,-0.303531,-0.934172,
     $                  -0.1875925,0.7946545,-0.577350,-0.1875925,
     $                  0.982247,0.0,0.1875925,0.303531,0.934172,
     $                  0.1875925,-0.7946545,0.57735,0.1875925,
     $                  -0.7946545,-0.57735,0.1875925,0.303531,
     $                  -0.934172,0.1875925,0.607062,0.0,0.7946545,
     $                  0.1875925,0.577350,0.7946545,-0.4911235,
     $                  0.356822,0.7946545,-0.4911235,-0.356822,
     $                  0.7946545,0.1875925,-0.577350,0.7946545/
      DATA              IDODA/1,2,3,4,5,2,7,13,8,3,3,8,14,9,4,4,9,15,10,
     $                  5,5,10,11,6,1,1,6,12,7,2,7,12,17,18,13,8,13,18,
     $                  19,14,9,14,19,20,15,10,15,20,16,11,6,11,16,17,
     $                  12,20,19,18,17,16/
      DATA              NN/1,4,5,1,3,4,4,10,5,1,5,6,1,6,2,1,2,3,4,3,9,4,
     $                  9,10,5,10,11,5,11,6,6,11,7,2,6,7,2,8,3,3,8,9,9,
     $                  12,10,10,12,11,2,7,8,9,8,12,11,12,7,12,8,7/
      DATA              TETRA/0.,0.,1.,.842809,0.,-.333333,-.4714045,
     $                  .816497,-.333333,-.4714045,-.816497,-.333333/
      DATA              ITETRA/1,2,3,4,1,3,1,4,2,4,3,2/
      DATA              HEXA/1.,1.,1.,-1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,1.,
     $                  1.,-1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,-1./
      DATA              IHEXA/1,4,3,2,1,5,8,4,1,2,6,5,2,3,7,6,3,4,8,7,5,
     $                  6,7,8/
      DATA              OKTA/0.,0.,1.,0.,0.,-1.,1.,0.,0.,0.,1.,0.,-1.,
     $                  0.,0.,0.,-1.,0./
      DATA              IOKTA/1,4,3,5,4,1,5,1,6,1,3,6,5,6,2,4,5,2,3,4,2,
     $                  6,3,2/
      DATA              LIN/1,1,1,1,0/
c     ..
c     .. executable statements ..
      CALL GR3OBB
      DO 230 L = 1,N
         CALL GR3TRC(ZENT(1,L),ZENT(2,L),ZENT(3,L),X,Y,Z)
         KI     = MOD(ABS(KIND(L)-1),5) + 1
         IF (KI.EQ.1) THEN
c------- tetraeder -----------------------
            NFACES = 4
            DO 10 I = 1,4
               Q(1,I) = X + TETRA(1,I)*R(L)
               Q(2,I) = Y + TETRA(2,I)*R(L)
               Q(3,I) = Z + TETRA(3,I)*R(L)
   10       CONTINUE
            DO 40 I = 1,NFACES
               DO 30 J = 1,3
                  DO 20 K = 1,3
                     P(K,J) = Q(K,ITETRA(J,I))
   20             CONTINUE
   30          CONTINUE
               CALL GR3OBP(AR,IER,P,3,1,LIN)
   40       CONTINUE
*
         ELSE IF (KI.EQ.2) THEN
c------- hexaeder (wuerfel ) -------------
            NFACES = 6
            DO 50 I = 1,8
               Q(1,I) = X + .57735*HEXA(1,I)*R(L)
               Q(2,I) = Y + .57735*HEXA(2,I)*R(L)
               Q(3,I) = Z + .57735*HEXA(3,I)*R(L)
   50       CONTINUE
            DO 80 I = 1,NFACES
               DO 70 J = 1,4
                  DO 60 K = 1,3
                     P(K,J) = Q(K,IHEXA(5-J,I))
   60             CONTINUE
   70          CONTINUE
               CALL GR3OBP(AR,IER,P,4,1,LIN)
   80       CONTINUE
*
         ELSE IF (KI.EQ.3) THEN
c------- oktaeder ------------------------
            NFACES = 8
            DO 90 I = 1,6
               Q(1,I) = X + OKTA(1,I)*R(L)
               Q(2,I) = Y + OKTA(2,I)*R(L)
               Q(3,I) = Z + OKTA(3,I)*R(L)
   90       CONTINUE
            DO 120 I = 1,NFACES
               DO 110 J = 1,3
                  DO 100 K = 1,3
                     P(K,J) = Q(K,IOKTA(4-J,I))
  100             CONTINUE
  110          CONTINUE
               CALL GR3OBP(AR,IER,P,3,1,LIN)
  120       CONTINUE
*
         ELSE IF (KI.EQ.4) THEN
c------- dodekaeder ----------------------
            NFACES = 12
            DO 130 I = 1,20
               Q(1,I) = X + DODA(1,I)*R(L)
               Q(2,I) = Y + DODA(2,I)*R(L)
               Q(3,I) = Z + DODA(3,I)*R(L)
  130       CONTINUE
            DO 170 I = 1,NFACES
               DO 150 J = 1,5
                  DO 140 K = 1,3
                     P(K,J) = Q(K,IDODA(6-J,I))
  140             CONTINUE
  150          CONTINUE
               DO 160 K = 1,3
                  P(K,6) = Q(K,IDODA(5,I))
  160          CONTINUE
               CALL GR3OBP(AR,IER,P,4,1,LIN(2))
               CALL GR3OBP(AR,IER,P(1,4),3,1,LIN(3))
  170       CONTINUE
*
         ELSE
c------- ikosaeder -----------------------
            NFACES = 20
            DO 190 I = 1,12
               SUX    = 0.
               SUY    = 0.
               SUZ    = 0.
               DO 180 J = 1,5
                  SUX    = SUX + DODA(1,IDODA(J,I))
                  SUY    = SUY + DODA(2,IDODA(J,I))
                  SUZ    = SUZ + DODA(3,IDODA(J,I))
  180          CONTINUE
               SUX    = SUX*0.2*1.25841*R(L)
               SUY    = SUY*0.2*1.25841*R(L)
               SUZ    = SUZ*0.2*1.25841*R(L)
               Q(1,I) = X + SUX
               Q(2,I) = Y + SUY
               Q(3,I) = Z + SUZ
  190       CONTINUE
            DO 220 I = 1,NFACES
               DO 210 J = 1,3
                  DO 200 K = 1,3
                     P(K,J) = Q(K,NN(4-J,I))
  200             CONTINUE
  210          CONTINUE
               CALL GR3OBP(AR,IER,P,3,1,LIN)
  220       CONTINUE
         END IF
*
  230 CONTINUE
      CALL GR3OBE(AR,IER,IFACES,JCO)
      END
      SUBROUTINE GR3ANT(AR,IER,PTXT,ANTXT,SIZTXT,PTRLEN,PTRANG,IDIREC,
     $                  IBOUND,IFONT,ICOL)
C---- GR3ANT : GR3 ANnotation Text
C---- AR     : Arbeitsspeicher der GR3-Routinen
C---- IER    : Fehlercode
C---- PTXT(3): Punkt im 3D des Textes
C---- ANTXT  : der Text selbst
C---- SIZTXT : Groesse in CM der Buchstaben, bezogen auf GRSCLC
C---- PTRLEN : Laenge des Zeigestriches in CM.
C---- PTRANG : Winkel in Grad des Zeigestriches mit der Waagerechten
C---- IDIREC : Richtung des Textes (0     waagerecht,
C----                               sonst senkrecht zum Zeigestrich)
C---- IBOUND : Lage des Textes zum Ende des Zeigestrichs
C----          0       : zentral
C----          negativ : rechtsbuendig
C----          positiv : linksbuendig
C---- IFONT  : Schriftart
C---- ICOL   : Farbe

C---- Strichlierung wie vorher mit GRDASH gesetzt

      INTEGER IER, IDIREC, IBOUND, ICOL
      REAL AR(*), PTXT(3), SIZTXT, PTRLEN, PTRANG
      CHARACTER(*) ANTXT

      PARAMETER (PIFAC=3.141593/180.)

C---- COMMON-Bloecke zu GR3ANT
      INTEGER            NUMAAN, LEMAAN
      PARAMETER          (NUMAAN=32,LEMAAN=32)
      INTEGER            NUA,INDA,LENA,IFOA,ICOA,IBOA
      REAL               SIZA,DS1A,DS2A,DS3A,ANGA,PTXA,PTYA
      CHARACTER(len=LEMAAN) CHAA
      COMMON /GRANTN/ NUA,INDA(NUMAAN),LENA(NUMAAN),SIZA(NUMAAN),
     $       IFOA(NUMAAN),DS1A(NUMAAN),DS2A(NUMAAN),DS3A(NUMAAN),
     $       ANGA(NUMAAN),PTXA(NUMAAN),PTYA(NUMAAN),ICOA(NUMAAN),
     $       IBOA(NUMAAN)
CDEC$ PSECT /GRANTN/ NOSHR
      COMMON /GRANTC/ CHAA(NUMAAN)
CDEC$ PSECT /GRANTC/ NOSHR

C---- COMMON mit Zeigern auf AR

      LOGICAL WIN,ZENPRO
      REAL RT(3,3),GR(3,3)
      COMMON /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LL,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR

C---- COMMON mit Plotparametern

      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRANTN/, /GRANTC/, /GR3SAV/, /GRPP/
      SAVE

C-----------------------------------------------------------------------
      IF (NOTDIM.NE.1234567890) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ANT RC=10: GR3DIM must first be called!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 999
*
      END IF

      NUA = NUA+1

      IF (NUA.GT.NUMAAN) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ANT RC=20: Zu viele Annotation-Texte! Max:',NUMAAN
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 20
         GO TO 999
*
      END IF

      IP0 = IP0+1

      IF (IP0.GT.MNP) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ANT RC=30: Zu viele Punkte in AR! LAR vergroessern!'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 30
         GO TO 999
*
      END IF

      INDA(NUA) = IP0
      AR(IP0      ) = PTXT(1)
      AR(IP0+MNP  ) = PTXT(2)
      AR(IP0+2*MNP) = PTXT(3)

      IF ( LEN(ANTXT).GT.LEMAAN ) THEN
         IF (IER.NE.2) WRITE (*,FMT=*)
     $       'GR3ANT RC=40: Annotation-Text zu lang! Max:',LEMAAN
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 40
         GO TO 999
*
      END IF

      LENA(NUA) = MIN(LEN(ANTXT),LEMAAN)
      CHAA(NUA) = ANTXT
      SIZA(NUA) = SIZTXT
      IFOA(NUA) = IFONT
      DS1A(NUA) = PP(10)
      DS2A(NUA) = PP(11)
      DS3A(NUA) = PP(12)
      IF ( IDIREC.EQ.0 ) THEN
         ANGA(NUA) = 0.
      ELSE
         ANGA(NUA) = PTRANG-90.
      ENDIF
      PTXA(NUA) = ABS(PTRLEN)*COS(PTRANG*PIFAC)
      PTYA(NUA) = ABS(PTRLEN)*SIN(PTRANG*PIFAC)
      ICOA(NUA) = ICOL
      IF ( IBOUND.EQ.0 ) THEN
         IBOA(NUA) = 2
      ELSE
         IF ( IBOUND.LT.0 ) THEN
            IBOA(NUA) = 3
         ELSE
            IBOA(NUA) = 1
         ENDIF
      ENDIF
  999 END
      SUBROUTINE GR3BKS(I1,I2,I3,I4,I5,I6,I7,I8)
c---- farben der rueckseite (backside)

c     .. scalar arguments ..
      INTEGER           I1,I2,I3,I4,I5,I6,I7,I8
c     ..
c     .. common blocks ..
      COMMON            /GR3BAC/ICOLBA,CENTR,WIND
CDEC$ PSECT /GR3BAC/ NOSHR
      REAL              CENTR
      REAL              WIND(4)
      INTEGER           ICOLBA(8)
c     ..
c     .. save statement ..
      SAVE /GR3BAC/
c     ..
c     .. executable statements ..

      IF (I1.GE.1 .AND. I1.LE.8) ICOLBA(1) = I1
      IF (I2.GE.1 .AND. I2.LE.8) ICOLBA(2) = I2
      IF (I3.GE.1 .AND. I3.LE.8) ICOLBA(3) = I3
      IF (I4.GE.1 .AND. I4.LE.8) ICOLBA(4) = I4
      IF (I5.GE.1 .AND. I5.LE.8) ICOLBA(5) = I5
      IF (I6.GE.1 .AND. I6.LE.8) ICOLBA(6) = I6
      IF (I7.GE.1 .AND. I7.LE.8) ICOLBA(7) = I7
      IF (I8.GE.1 .AND. I8.LE.8) ICOLBA(8) = I8
      END
c-----------------------------------------------------------------------
c only for ibm (must be deleted on cray)
      INTEGER FUNCTION GRISMI(N,V,IV)
c     .. scalar arguments ..
      INTEGER                IV,N
c     ..
c     .. array arguments ..
      REAL                   V(1+ (N-1)*IV)
c     ..
c     .. local scalars ..
      REAL                   VMIN
      INTEGER                I,J
c     ..
c     .. executable statements ..
      J      = 1
      IF (IV.LT.0) J      = 1 - (N-1)*IV
      VMIN   = V(J)
      GRISMI  = 1
      DO 10 I = 2,N
         J      = J + IV
         IF (V(J).LT.VMIN) THEN
            VMIN   = V(J)
            GRISMI  = I
         END IF
*
   10 CONTINUE
      END
c--------------------------------------------------------------------
c   only for ibm (must be deleted on cray)
      INTEGER FUNCTION GRISMA(N,V,IV)
c     .. scalar arguments ..
      INTEGER                IV,N
c     ..
c     .. array arguments ..
c     REAL                   V(1+ (N-1)*IV)
      REAL                   V(*)
c     ..
c     .. local scalars ..
      REAL                   VMAX
      INTEGER                I,J
c     ..
c     .. executable statements ..
      J      = 1
      IF (IV.LT.0) J      = 1 - (N-1)*IV
      VMAX   = V(J)
      GRISMA  = 1
      DO 10 I = 2,N
         J      = J + IV
         IF (V(J).GT.VMAX) THEN
            VMAX   = V(J)
            GRISMA  = I
         END IF
*
   10 CONTINUE
      END
C=======================================================================
      SUBROUTINE GR3TUB(AR,IER,POI,N,M,R,CLOSED,IFACES,JCO)
c     .. scalar arguments ..
      REAL              AR,R
      INTEGER           IER,IFACES,JCO,N,M
      LOGICAL CLOSED
c     ..
c     .. array arguments ..
      REAL              POI(3,N)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3SUB,GR3TBS
c     ..
c     .. executable statements ..

      CALL GR3SUB(AR,IER,GR3TBS,JCO,POI(1,1),N,M,R,CLOSED,IFACES,
     $            I1,I2,I3,R4)
      END
