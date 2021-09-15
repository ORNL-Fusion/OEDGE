C EIRENE07 COMPILATION
C ===== SOURCE: algebr.f
C
C
C-----------------------------------------------------------------------
      SUBROUTINE ALGEBR (TERM,OPER,IZIF,CONST,NOP)
C-----------------------------------------------------------------------
C
C     AUTOR:           ST. HUBER
C     IHK-KENNZIFFER:  121
C     DATUM:           25-NOV-1988
C
C     FUNKTION:
C
C     DAS PROGRAMM LIEST ARITHMETISCHE AUSDRUECKE EIN,
C     RUFT DAS UNTERPROGRAMM ZERLEG AUF
C     UND GIBT ENTWEDER DIE EINZELNEN ZERLEGUNGEN ODER DIE
C     REGELVERLETZUNG AUS
C
C-----------------------------------------------------------------------
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
C
C     ARGUMENT:
C
         CHARACTER(*), INTENT(INOUT) :: TERM
C           : EINZULESENDER AUSDRUCK

         CHARACTER(1), INTENT(OUT) :: OPER(*)
         INTEGER, INTENT(OUT) :: IZIF(4,*)
         REAL(DP), INTENT(INOUT) :: CONST(*)
         INTEGER, INTENT(OUT) :: NOP
C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: ZMAX = 20
C           : ANZAHL DER MAXIMALEN ZERLEGUNGEN

         INTEGER, PARAMETER :: MAXLEN = 72
C           : MAXIMALE STRINGLAENGE

C
C     LOKALE VARIABLEN :
C
         INTEGER ::  LAENGE
C           : AKTUELLE LAENGE VON TERM

         CHARACTER(MAXLEN) :: HLFTERM
C           : HILFSSTRING ZUM UMSPEICHERN

         CHARACTER(MAXLEN+2):: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         INTEGER :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         INTEGER :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

         CHARACTER(MAXLEN) :: PART(ZMAX)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER :: IPART(ZMAX)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(MAXLEN) :: ARITH(ZMAX)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER :: IARITH(ZMAX)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

         CHARACTER(MAXLEN) :: HILFE
C           : ARBEITSSPEICHER FUER UNTERPROGRAMM ZERLEG

         INTEGER :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

CHR
CHR      VARIABLEN ZUR MODIFIKATION DES PROGRAMMES
         INTEGER :: NR, ANFANG, ENDE, FELDIND, IK, IKM, IKP
         CHARACTER(4) :: FO
         CHARACTER(12) :: ERSETZ(10)
         character(10) :: buchst
chr
C
C     HILFSVARIABLEN :
C
         INTEGER :: MAXI, I
chr
chr   string, der die neuen variablennamen enthaelt
      buchst='ABCDEFGHIJ'
chr
C
C        LESE TERM UND WERTE AUS
C
         IF (TERM .NE. ' ') THEN
C
C           VERARBEITUNG, FUER DEN FALL, DASS KEINE LEERZEILE
C           EINGELESEN WURDE
C
chr         vorbereiten des terms fuer die weitere verarbeitung, d.h.
chr         bringen der operanden in die vom programm verlangte
chr         zweistellige alphabetische form
            nr=0
101         anfang=index(term,'<')
            if (anfang.ne.0) then
               nr=nr+1
               ende=index(term,'>')
chr            abspeichern der ersetzten operanden
chr            der i-te operand wird hierbei in der i-ten dimension
chr            des feldes ersetz abgelegt und durch den i-ten
chr            buchstaben im alphabet ersetzt. gleichheit von
chr            operanden wird hierbei nicht beruecksichtigt
               ersetz(nr)=term(anfang:ende)
               if (anfang.gt.1) then
               hlfterm=term(1:anfang-1)//buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
               else
               hlfterm=buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
               endif
               term=hlfterm
               goto 101
            endif
chr
C
C           ERMITTELN DER LAENGE VON TERM
C
            LAENGE=LEN(TERM)
20          IF (TERM(LAENGE:LAENGE) .EQ. ' ') THEN
               LAENGE=LAENGE-1
               GOTO 20
            ENDIF

            AUSDRU=TERM
            AKTLEN=LAENGE

            CALL ZERLEG(AUSDRU, AKTLEN, IPART, PART, IARITH, ARITH,
     >                  TEIL, HILFE, ERROR)

            IF (ERROR .EQ. 0) THEN
C
C              AUSGABE DER ZERLEGUNG
C
               MAXI=0
               DO 45, I=1,TEIL
                   MAXI=MAX(IPART(I),MAXI)
45                 CONTINUE

               NOP=TEIL
               DO 30, I=1,TEIL
chr               ausgabe der zerlegung in der form:
chr                    operator ziffer1 ziffer2 ziffer3 ziffer4
chr               wobei jeweils die 1. und 2. sowie die 3. und 4.
chr               ziffer eine einheit bilden.
chr               entweder stellt eine solche einheit einen operanden
chr               dar oder ein zwischenergebnis mit der 1. ziffer als
chr               nummer und der 2. als 0 zur kennzeichnung des paares
chr               als zwischenergebnis
                  if (part(i)(7:7).ne.'Z') then
                     OPER(I)=PART(I)(9:9)
                     feldind=index(buchst,part(i)(7:7))
                     IK=INDEX(ERSETZ(FELDIND),',')
                     IF (IK.EQ.0) THEN
                       IZIF(1,I)=-I
                       IZIF(2,I)=0
                       CALL RDCN (ERSETZ(FELDIND),CONST(I))
                     ELSE
                       IKM=IK-1
                       IKP=IK+1
                       FO(1:4)='(I )'
                       WRITE (FO(3:3),'(I1)') IK-2
                       READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(1,I)
                       IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                         IZIF(2,I)
                       ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                         IZIF(2,I)
                       ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                         IZIF(2,I)
                       ENDIF
                     ENDIF

                     IF (PART(I)(10:10).NE.'Z') THEN
                       FELDIND=INDEX(BUCHST,PART(I)(10:10))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IZIF(3,I)=-I
                         IZIF(4,I)=0
                         CALL RDCN (ERSETZ(FELDIND),CONST(I))
                       ELSE
                         IKM=IK-1
                         IKP=IK+1
                         FO(1:4)='(I )'
                         WRITE (FO(3:3),'(I1)') IK-2
                         READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(3,I)
                         IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                          IZIF(4,I)
                         ENDIF
                       ENDIF
                     ELSE
                        READ(PART(I)(12:12),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     ENDIF
                  ELSE
                     OPER(I)=PART(I)(10:10)
                     READ(PART(I)(9:9),'(I1)') IZIF(1,I)
                     IZIF(2,I)=0
                     if (part(i)(11:11).ne.'Z') then
                       FELDIND=INDEX(BUCHST,PART(I)(11:11))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IZIF(3,I)=-I
                         IZIF(4,I)=0
                         CALL RDCN (ERSETZ(FELDIND),CONST(I))
                       ELSE
                         IKM=IK-1
                         IKP=IK+1
                         FO(1:4)='(I )'
                         WRITE (FO(3:3),'(I1)') IK-2
                         READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(3,I)
                         IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                          IZIF(4,I)
                         ENDIF
                       ENDIF
                     else
                        READ(PART(I)(13:13),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     endif
                  endif
30                CONTINUE
            ELSE
C
C              AUSGABE DER FEHLERMELDUNG
C
               WRITE(iunout,'(2A)') '0FOLGENDE REGELVERLETZUNG ',
     >                        'WURDE ERKANNT:'
               CALL MECKER(ERROR)
               NOP=0
            ENDIF
         ENDIF
C
C     ENDE VON ALGEBR
C
      RETURN
      END
C ===== SOURCE: arellp.f
C
C
*DK ARELLP
      SUBROUTINE ARELLP
C  INPUT:
     >                 (EPS1,EPS2,ELL1,ELL2,TR1,TR2,
     >                  XHALB1,XHALB2,ALPHA,BETA,IFLAG,
C  OUTPUT:
     >                  AELL,SX,SY,
     >                  X1ALPHA,Y1ALPHA,X2ALPHA,Y2ALPHA,
     >                  X1BETA,Y1BETA,X2BETA,Y2BETA)
C
C  ELLIPSE, MIT TRIANGULARITAET, LT. L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,
C                                    PHYS.FLUIDS, 24,8,1981,P1431
C  EQS.(110/111) AND (A9)
C
C  X=EPS+     R*COS(T)+TR*(COS(2T)
C  Y=    ELL*(R*SIN(T)-TR*(SIN(2T))
C
C VARIABLENBESCHREIBUNG
C
C UEBERGABEPARAMETER
C
C INPUT:
C
C EPS1   X-KOORDINATE DES MITTELPUNKTES DER GROESSEREN ELLIPSE
C EPS2   X-KOORDINATE DES MITTELPUNKTES DER KLEINEREN ELLIPSE
C        DIE Y-KOORDINATEN DER MITTELPUNKTE MUESSEN =0 SEIN
C ELL1   ELLIPTIZITAET YHALB1 / XHALB1
C ELL2   ELLIPTIZITAET YHALB2 / XHALB2
C TR1    TRIANGULARITAET
C TR2    TRIANGULARITAET
C XHALB1 HALBACHSE DER GROESSEREN ELLIPSE IN X-RICHTUNG
C XHALB2 HALBACHSE DER KLEINEREN ELLIPSE IN X-RICHTUNG
C ALPHA  GROESSERER DER WINKEL, ZWISCHEN 0 UND 2*PI, VON +X-ACHSE
C BETA   KLEINERER DER WINKEL, ZWISCHEN 0 UND 2*PI, VON +X-ACHSE
C
C OUTPUT:  IFLAG=1 OR IFLAG=0
C
C AELL   FLAECHENINHALT DER ZU BERECHNENDEN FLAECHE
C SX     X-WERT DES SCHWERPUNKTES DER ZU BERECHNENDEN FLAECHE
C SY     Y-WERT DES SCHWERPUNKTES DER ZU BERECHNENDEN FLAECHE
C
C OUTPUT:  IFLAG=2 OR IFLAG=0
C
C X1ALPHA X-WERT DES SCHNITTPUNKTES DER GROESSEREN ELLIPSE MIT DER
C         GERADEN, DIE DURCH DEN GROESSEREN WINKEL GEGEBEN IST.
C Y1ALPHA Y-WERT DES SCHNITTPUNKTES DER GROESSEREN ELLIPSE MIT DER
C         GERADEN, DIE DURCH DEN GROESSEREN WINKEL GEGEBEN IST.
C X2ALPHA X-WERT DES SCHNITTPUNKTES DER KLEINEREN ELLIPSE MIT DER
C         GERADEN, DIE DURCH DEN GROESSEREN WINKEL GEGEBEN IST.
C Y2ALPHA Y-WERT DES SCHNITTPUNKTES DER KLEINEREN ELLIPSE MIT DER
C         GERADEN, DIE DURCH DEN GROESSEREN WINKEL GEGEBEN IST.
C X1BETA X-WERT DES SCHNITTPUNKTES DER GROESSEREN ELLIPSE MIT DER
C        GERADEN, DIE DURCH DEN KLEINEREN WINKEL GEGEBEN IST.
C Y1BETA Y-WERT DES SCHNITTPUNKTES DER GROESSEREN ELLIPSE MIT DER
C        GERADEN, DIE DURCH DEN KLEINEREN WINKEL GEGEBEN IST.
C X2BETA X-WERT DES SCHNITTPUNKTES DER KLEINEREN ELLIPSE MIT DER
C        GERADEN, DIE DURCH DEN KLEINEREN WINKEL GEGEBEN IST.
C Y2BETA Y-WERT DES SCHNITTPUNKTES DER KLEINEREN ELLIPSE MIT DER
C        GERADEN, DIE DURCH DEN KLEINEREN WINKEL GEGEBEN IST.
C
C
      USE PRECISION
      IMPLICIT NONE
C  INPUT:
      REAL(DP), INTENT(IN) :: EPS1, EPS2, ELL1, ELL2, TR1, TR2,
     >                      XHALB1, XHALB2, ALPHA, BETA
      INTEGER, INTENT (IN) :: IFLAG
C  OUTPUT:
      REAL(DP), INTENT(OUT) :: AELL, SX, SY,
     >                       X1ALPHA, Y1ALPHA, X2ALPHA, Y2ALPHA,
     >                       X1BETA, Y1BETA, X2BETA, Y2BETA

      REAL(DP) :: INTSIN, INTCOS, PI, SIN2, COS2, COS3, COS2T, COSSIN,
     .          COS2SIN, T1, T2, T3, T4, T5, T6, T7, T8, T9, TRYA, TRXS,
     .          TRXA, TRYS, EPPA, YHALB1, YHALB2, TRX1, SIN2TSIN,
     .          COS22T, SIN22T, COS2TCOS, ELLS, ELLA, EPPS, DENOM, TRX2,
     .          TRY1, TRY2
      REAL(DP) ::
     .          T1458, T1464, T1465, T1457, T1441, T1444, T1448, T1509,
     .          T1516, T1525, T1502, T1469, T1472, T1482, T1440, T1234,
     .          T1235, T1257, T1227, T1206, T1220, T1223, T1432, T1433,
     .          T1436, T1338, T1276, T1294, T1335, T1528, T1531, T1003,
     .          T1007, T1008, T1000, T1129, T1136, T1140, T1120, T1111,
     .          T1114, T1115, T1192, T1196, T1199, T1176, T1163, T1167,
     .          T1168, T1106, T1049, T1060, T1071, T1038, T1022, T1026,
     .          T1031, T1101, T1102, T1105, T1098, T1081, T1085, T1090,
     .          T1473
      REAL(DP) :: T68, T72, T76, T64, T54, T59, T60, T91, T87, T79, T83,
     .          T86, T49, T19, T21, T25, T11, T10, T40, T41, T45, T37,
     .          T27, T28, T31, T58, T61, T62, T57, T51, T53, T56, T78,
     .          T81, T85, T75, T63, T67, T70, T48, T20, T23, T24, T16,
     .          T13, T14, T15, T38, T42, T44, T35, T26, T30, T33, T96,
     .          T99, T95, T88, T89, T94, T12
      REAL(DP) :: T103, T107, T938, T939, T943, T934, T914, T917, T933,
     .          T947, T973, T996, T910, T892, T893, T894, T889, T877,
     .          T878, T888, T905, T906, T909, T904, T898, T899, T900,
     .          T632, T633, T634, T619, T602, T603, T613, T655, T659,
     .          T667, T650, T638, T641, T647, T598, T571, T575, T580,
     .          T568, T564, T565, T566, T595, T596, T597, T592, T584,
     .          T585, T588, T821, T824, T825, T786, T766, T767, T770,
     .          T911, T993, T994, T895, T828, T842, T879, T759, T684,
     .          T687, T693, T683, T670, T671, T678, T743, T747, T754,
     .          T742, T700, T708, T709, T221, T225, T230, T218, T209,
     .          T213, T214, T274, T275, T279, T267, T231, T242, T261,
     .          T206, T164, T168, T169, T162, T108, T130, T146, T187,
     .          T188, T196, T184, T173, T178, T183, T440, T450, T463,
     .          T439, T387, T402, T559, T560, T561, T555, T499, T553,
     .          T554, T386, T329, T334, T337, T322, T296, T301, T303,
     .          T374, T377, T379, T367, T338, T348, T361, T134, T140
      REAL(DP) :: T133, T125, T126, T131, T154, T155, T156, T148, T144,
     .          T145, T147, T124, T100, T118, T119, T123, T114, T105,
     .          T109, T113, T692, T701, T704, T660, T651, T652, T656,
     .          T731, T734, T737, T730, T722, T728, T729, T637, T607,
     .          T610, T611, T606, T562, T574, T601, T629, T630, T631,
     .          T626, T618, T622, T623, T834, T845, T838, T830, T809,
     .          T815, T816, T867, T871, T872, T862, T846, T858, T859,
     .          T806, T775, T779, T780, T774, T752, T761, T764, T798,
     .          T801, T805, T797, T784, T785, T259, T260, T262, T252,
     .          T247, T248, T251, T286, T294, T298, T271, T264, T265,
     .          T269, T243, T192, T200, T201, T167, T159, T160, T163,
     .          T224, T228, T229, T220, T211, T215, T219, T435, T438,
     .          T441, T427, T394, T416, T426, T481, T508, T520, T452,
     .          T442, T446, T449, T383, T345, T346, T349, T341, T305,
     .          T315, T332, T372, T375, T376, T371, T358, T359, T362,
     .          T410, T835, T141, T793

C     BERECHNUNG VON PI
      PI = 4.*ATAN(1.)
C
      INTSIN = COS(BETA)-COS(ALPHA)
      INTCOS = SIN(ALPHA) - SIN(BETA)
      COSSIN = 0.5*(SIN(ALPHA)*SIN(ALPHA) - SIN(BETA)*SIN(BETA))
      COS2SIN = 1./3.*(-COS(ALPHA)*COS(ALPHA)*COS(ALPHA) +
     .                 COS(BETA)*COS(BETA)*COS(BETA))
      COS2 = 0.5*(ALPHA+SIN(ALPHA)*COS(ALPHA)-BETA-SIN(BETA)*COS(BETA))
      SIN2 = 0.5*(ALPHA-SIN(ALPHA)*COS(ALPHA)-BETA+SIN(BETA)*COS(BETA))
      COS3 = SIN(ALPHA) - 1./3.*SIN(ALPHA)*SIN(ALPHA)*SIN(ALPHA) -
     .       SIN(BETA) + 1./3.*SIN(BETA)*SIN(BETA)*SIN(BETA)
      COS2T = 0.5*(SIN(2.*ALPHA) - SIN(2.*BETA))
      COS22T = 0.5*(ALPHA+SIN(4.*ALPHA)/4.-BETA-SIN(4.*BETA)/4.)
      SIN22T = 0.5*(ALPHA-SIN(4.*ALPHA)/4.-BETA+SIN(4.*BETA)/4.)
      COS2TCOS = 1./3.*(2.*COS(ALPHA)*SIN(2.*ALPHA) -
     .                  SIN(ALPHA)*COS(2.*ALPHA) -
     .                  2.*COS(BETA)*SIN(2.*BETA) +
     .                  SIN(BETA)*COS(2.*BETA))
      SIN2TSIN = 1./3.*(COS(ALPHA)*SIN(2.*ALPHA) -
     .                  2.*SIN(ALPHA)*COS(2.*ALPHA) -
     .                  COS(BETA)*SIN(2.*BETA) +
     .                  2.*SIN(BETA)*COS(2.*BETA))
C
C  BERECHNUNG DER HALBACHSEN IN Y-RICHTUNG
      YHALB1 = ELL1*XHALB1
      YHALB2 = ELL2*XHALB2
      TRX1   = TR1
      TRX2   = TR2
      TRY1   = -ELL1*TR1
      TRY2   = -ELL2*TR2
C
      DENOM=(XHALB1 - XHALB2)+1.D-30
      ELLS = (YHALB1 - YHALB2) / DENOM
      ELLA = (YHALB2*XHALB1 - YHALB1*XHALB2) / DENOM
      EPPS = (EPS1 - EPS2) / DENOM
      EPPA = (EPS2*XHALB1 - EPS1*XHALB2) / DENOM
      TRXS = (TRX1 - TRX2) / DENOM
      TRXA = (TRX2*XHALB1 - TRX1*XHALB2) / DENOM
      TRYS = (TRY1 - TRY2) / DENOM
      TRYA = (TRY2*XHALB1 - TRY1*XHALB2) / DENOM
C
      IF (IFLAG.EQ.2) GOTO 1000
C
C     IF (NLTRI) THEN
C
C  BERECHNUNG DER FLAECHE
C
          AELL = 0.5*(XHALB1**2-XHALB2**2)*(EPPS*ELLS*INTCOS +
     .           2.*EPPS*TRYS*COS2T + ELLS*COS2 +
     .           (2.*TRYS+TRXS*ELLS)*COS2TCOS + 2.*TRYS*TRXS*COS22T+
     .           ELLS*SIN2 + (2.*ELLS*TRXS+TRYS)*SIN2TSIN +
     .           2.*TRYS*TRXS*SIN22T) +
     .           (XHALB1-XHALB2)*(EPPS*ELLA*INTCOS+2.*EPPS*TRYA*COS2T+
     .           ELLA*COS2 + (2.*TRYA+TRXS*ELLA)*COS2TCOS +
     .           2.*TRXS*TRYA*COS22T +
     .           2.*ELLS*TRXA*SIN2TSIN + 2.*TRYS*TRXA*SIN22T)
C
C  BERECHNUNG DES SCHWERPUNKTES
C
      SX = 0.
      SY = 0.
      t1 = 2*alpha
      t2 = sin(t1)
      t3 = t2**2
      t4 = t3*t2
      t5 = trxs**2
      t6 = xhalb1**2
      t7 = t6*xhalb1
      t8 = trys*t7
      t9 = t5*t8
      t12 = epps**2
      t13 = xhalb2**2
      t14 = t13*xhalb2
      t15 = trys*t14
      t16 = t12*t15
      t20 = trys*epps*t14
      t23 = sin(alpha)
      t24 = eppa*xhalb1
      t26 = epps*ella*t24
      t30 = trxs*trys*trxa*t13
      t33 = eppa*xhalb2
      t35 = epps*ella*t33
      t38 = epps*t8
      t42 = trxa**2
      t44 = trys*t42*xhalb1
      t48 = t12*trya*t13
      t51 = eppa*t13
      t53 = epps*ells*t51
      t56 = eppa*t6
      t57 = trys*t56
      t58 = epps*t57
      t61 = t23**2
      t62 = t61*t23
      t63 = ells*t7
      t67 = epps*trya*t24
      t70 = epps*ells*t56
      t75 = t12*t8
      t78 = t12*t63
      t81 = ells*t14
      t85 = trys*t42*xhalb2
      t88 = t13*t2
      t89 = trxa*t88
      t94 = sin(3*alpha)
      t95 = xhalb1*t94
      t96 = trxa*t95
      t99 = xhalb1*t23
      t100 = eppa*t99
      t105 = eppa*t95
      t109 = t7*t2
      t113 = t6*t23
      t114 = trxs*t113
      t118 = t6*t94
      t119 = trxs*t118
      t123 = cos(t1)
      t124 = t123**2
      t125 = t124*t2
      t126 = t13*t125
      t131 = trya*t88
      t133 = t14*t23
      t134 = ells*t133
      t140 = sin(4*alpha)
      t141 = t7*t140
      t144 = t6*t2
      t145 = trya*t144
      t147 = cos(alpha)
      t148 = t23*t147
      t154 = sin(5*alpha)
      t155 = xhalb1*t154
      t156 = trxa*t155
      t159 = t13*t23
      t160 = trxs*t159
      t163 = eppa*t113
      t167 = eppa*t118
      t192 = trxa*xhalb1*t2
      t200 = t6*t154
      t201 = trxs*t200
      t211 = t14*t2
      t215 = t14*t94
      t219 = t2*t123
      t220 = t6*t219
      t224 = trxa*t113
      t228 = t147**2
      t229 = t228*t23
      t243 = t13*t154
      t247 = t13*t94
      t248 = trxs*t247
      t251 = trxa*t159
      t252 = ells*t251
      t259 = t7*t23
      t260 = ells*t259
      t262 = trys*t259
      t264 = t7*t94
      t265 = trys*t264
      t269 = ella*t159
      t271 = t14*t140
      t286 = xhalb1*t219
      t294 = epps*t133
      t298 = t6*t125
      t305 = trxa*t144
      t315 = trxa*t118
      t332 = trys*trxs*t14*t219
      t341 = t7*t154
      t345 = trxa*t99
      t346 = ella*t345
      t349 = ella*t96
      t358 = trxa*t247
      t359 = ells*t358
      t362 = t14*t154
      t371 = xhalb2*t154
      t372 = trxa*t371
      t375 = xhalb2*t23
      t376 = trxa*t375
      t383 = t6*t140
      t394 = trxa*t243
      t416 = trxa*xhalb2*t2
      t426 = xhalb2*t94
      t427 = trxa*t426
      t435 = ella*t113
      t438 = ella*t427
      t441 = xhalb2*t219
      t442 = eppa*t441
      t446 = eppa*t375
      t449 = eppa*t426
      t452 = t7*t219
      t481 = ella*t376
      t508 = eppa*t159
      t520 = eppa*t247
      t562 = eppa*t286
      t574 = t13*t140
      t601 = t13*t219
      t606 = t6*alpha
      t607 = trxa*t606
      t610 = t13*alpha
      t611 = trxa*t610
      t618 = t7*alpha
      t622 = xhalb1*alpha
      t623 = trxa*t622
      t626 = eppa*t622
      t629 = t14*alpha
      t630 = trxs*t629
      t631 = trys*t630
      t637 = eppa*t606
      t651 = xhalb2*alpha
      t652 = eppa*t651
      t656 = trxs*t610
      t660 = trxs*t606
      t692 = trxa*t651
      t701 = trys*t618
      t704 = eppa*t610
      t722 = sin(beta)
      t728 = 2*beta
      t729 = sin(t728)
      t730 = t729**2
      t731 = t730*t729
      t734 = trys*trxa*trxs*t6
      t737 = trys*t51
      t752 = t12*trya*t6
      t761 = epps*t737
      t764 = t5*t15
      t774 = sin(3*beta)
      t775 = xhalb2*t774
      t779 = t13*t722
      t780 = eppa*t779
      t784 = t14*t722
      t785 = epps*t784
      t793 = ella*t779
      t797 = sin(4*beta)
      t798 = t13*t797
      t801 = ells*t784
      t805 = t722**2
      t806 = t805*t722
      t809 = t14*t729
      t815 = cos(beta)
      t816 = t722*t815
      t830 = eppa*t775
      t834 = cos(t728)
      t835 = t729*t834
      t838 = trys*trxs*t14*t835
      t846 = t12*ella*t13
      t858 = xhalb2*t722
      t859 = eppa*t858
      t862 = xhalb2*t835
      t867 = t13*t729
      t871 = t834**2
      t872 = t871*t729
      t877 = xhalb1*t835
      t878 = eppa*t877
      t888 = trxa*t858
      t889 = ella*t888
      t892 = xhalb1*t774
      t893 = trxa*t892
      t894 = ella*t893
      t898 = sin(5*beta)
      t899 = xhalb2*t898
      t900 = trxa*t899
      t904 = xhalb1*t722
      t905 = trxa*t904
      t906 = ella*t905
      t909 = trxa*t775
      t910 = ella*t909
      t914 = epps*trya*t33
      t917 = trxa*xhalb2*t729
      t933 = t6*t774
      t934 = trxs*t933
      t938 = t6*t722
      t939 = trxs*t938
      t943 = eppa*t904
      t947 = eppa*t892
      t973 = t6*t797
      t996 = xhalb1*t898
      t1000 = trxa*t996
      t1003 = trxs*t779
      t1007 = t13*t774
      t1008 = trxs*t1007
      t1022 = t6*t729
      t1026 = t6*t835
      t1031 = t14*t797
      t1038 = eppa*t1007
      t1049 = trxa*t867
      t1060 = trxa*xhalb1*t729
      t1071 = t13*t835
      t1081 = eppa*t862
      t1085 = t7*t722
      t1090 = t14*t774
      t1098 = trxa*t1022
      t1101 = trxa*t779
      t1102 = ells*t1101
      t1105 = trxa*t1007
      t1106 = ells*t1105
      t1111 = trys*t1085
      t1114 = t7*t774
      t1115 = trys*t1114
      t1120 = t13*t872
      t1129 = t7*t835
      t1136 = t13*t898
      t1140 = t14*t898
      t1163 = t7*t898
      t1167 = t6*t898
      t1168 = trxs*t1167
      t1176 = trya*t867
      t1192 = t6*t872
      t1196 = trya*t1022
      t1199 = trxa*t933
      t1206 = trxa*t938
      t1220 = ells*t1085
      t1223 = trxa*t1136
      t1227 = t7*t797
      t1234 = t815**2
      t1235 = t1234*t722
      t1257 = eppa*t933
      t1276 = t7*t729
      t1294 = eppa*t938
      t1335 = t12*t81
      t1338 = ella*t938
      t1432 = xhalb2*beta
      t1433 = eppa*t1432
      t1436 = t14*beta
      t1440 = t13*beta
      t1441 = trxs*t1440
      t1444 = trxa*t1432
      t1448 = t12*ella*t6
      t1457 = trxs*t1436
      t1458 = trys*t1457
      t1464 = xhalb1*beta
      t1465 = trxa*t1464
      t1469 = eppa*t1464
      t1473 = t6*beta
      t1482 = trxs*t1473
      t1502 = trxa*t1440
      t1509 = eppa*t1440
      t1516 = t7*beta
      t1525 = eppa*t1473
      t1528 = trxa*t1473
      t1531 = trys*t1516
      SX = SX -TRYS*TRXA*T1469-TRXS*TRYS*T1525+TRXS*TRYA*T562/2-EPPS*ELL
     +s*t1516/2+epps*ella*t1440/2-ells*trxs*t1516/4-epps*ells*t7*t816/6+
     +epps*ella*t13*t816/2+trxs*trys*trxa*t1120/6-3.0/40.0*trxs*ells*t12
     +23+trys*trxs*t1090/36+ells*t14*t1235/9+3.0/40.0*trys*trxa*t1167-tr
     +ys*trxs*t1140/60-2.0/3.0*trxs*trya*t1060-trys*trxa*t1026/4+trxs*tr
     +ys*t1049/3-ells*trxs*t1276/4-trxs*trys*t704+ells*trxs*t1257/12-trx
     +s*trys*t1098/3-t5*ella*t933/24+ells*trxs*t1227/48+trxs*ella*t859/2
     +-ella*trxa*xhalb2*t140/16+epps*trya*trxa*t286/2-t5*ella*t247/24-tr
     +xs*trys*trxa*t126/6-epps*trya*trxa*t877/2+trxs*trya*trxa*xhalb1*t1
     +25/3-ells*epps*t14*t148/6+trya*t1440/4+trxs*ella*t830/6-ells*t42*t
     +892/6+ells*t42*t996/10+epps*ella*t1003/2+epps*ella*t1008/6-ells*t1
     +525/2-ella*trxs*t973/16-t5*ells*t215/12+t1458/3+t23*t1448/2-trxs*t
     +rys*epps*t1129/6-epps*trya*t933/3+t5*t145/3-epps*trya*t938-ella*tr
     +xs*t574/16+ella*trxa*xhalb2*t797/16-trxs*ella*t947/6+t731*t764/9-t
     +2*t761/2+t729*t16/3-trxs*ella*t943/2-epps*ella*t939/2-ella*t1482/4
     +-epps*ella*t934/6-t1338/3-t5*t1338/4+t2*t752/2+ells*trxs*t809/4-t5
     +*ells*t1140/60+t5*ells*t1090/12+trya*trxs*t1136/10-ells*trxa*t943+
     +trys*trxa*t878/2+trxs*trys*t1163/60+ells*trxs*t271/48+t5*ells*t362
     +/60-trxs*trys*t1129/6-t5*trya*t1192/6
      SX = SX+EPPS*ELLS*T1199/12+2.0/3.0*T
     +rys*trxs*t784-3.0/4.0*epps*ells*t1206+t5*ella*t1007/24+t5*ells*t11
     +63/60-t5*ells*t1114/12+t5*ella*t1136/40+2.0/3.0*trxs*trya*t917-ell
     +s*trxs*t211/4+t2*t57/2-t729*t57/2-t5*trys*t14*t125/9-t722*t78/3-tr
     +ya*t372/10+epps*t359/12-t731*t44/3-trys*t1294/4+trya*t345-trxs*try
     +a*t878/2-t722*t1448/2+t2*t75/3+ella*trxs*t867/4+t23*t78/3-t729*t75
     +2/2+t23*t70/2+t2*t67+t62*t63/9+ells*trxa*t859+t729*t914+trxs*ella*
     +t900/20+ella*t1433/2-t722*t26+epps*t910/6+trya*t830/3+ella*t6*t229
     +/6+epps*t346/2-ells*trxa*t830/3+ella*trxs*t798/16-ella*t416/4+2.0/
     +9.0*t5*trys*t809+2.0/3.0*ells*trxs*epps*t259+2.0/3.0*ells*trxs*t78
     +5+trya*t888-t731*t9/9+trys*t1101/2+3.0/4.0*trxs*ells*t780+trxs*t26
     +5/36+ells*t42*t775/6+ella*t1441/4+trys*t784/3-t722*t70/2-ella*t146
     +9/2-trys*trxa*t601/4-trys*trxa*t652-epps*trya*t692-t729*t75/3+trya
     +*t1003+trys*t1257/12+ells*t1457/4-2.0/3.0*ells*trxs*epps*t1085-ell
     +a*trxa*xhalb1*t797/16-ells*epps*t629/2-t1531/12+t729*t737/2+trya*t
     +900/10-t23*t846/2+t4*t734/3+t435/3+ella*t1444/4+ells*trxs*t618/4+e
     +lls*t305/4+trxs*t910/12+trys*t358/24-trys*t133/3+trys*trxa*t626-ep
     +ps*t906/2+ells*t1502/4+epps*trya*t660+trxs*t1102/4-trys*t1038/12+t
     +838/6-trxs*t1531/3-epps*trya*t656-trxs*trya*t652
      SX = SX+EPPS*ELLA*T606/2+
     +trxs*trya*t626+epps*t889/2-trys*t1031/48-t731*t734/3+trxs*trys*t63
     +7+trys*t1227/48-epps*t438/6-3.0/4.0*epps*t252+epps*trya*t623+t722*
     +t53/2+epps*ells*t618/2+trys*trxa*epps*t606+trya*t1008/6-epps*ella*
     +t610/2-t5*trya*t126/6-t729*t67-t2*t48/2-t23*t53/2-trxs*ella*t372/2
     +0+t4*t44/3-trys*t251/2+trya*t606/4+ells*trxa*t574/16+t2*t38/3-ells
     +*trxa*t105/3-t23*t35-trya*t610/4-t4*t764/9-t4*t30/3-trys*trxa*t562
     +/2-ella*t692/4-3.0/4.0*trxs*ells*t508+trys*trxa*epps*t1440-trys*tr
     +xa*epps*t610-ella*t652/2+t729*t761/2+ells*t1509/2-trxs*trys*epps*t
     +1516-t5*t269/4+5.0/24.0*trxs*t1106+t701/12-trxs*trya*t1469-epps*el
     +la*t1473/2-epps*trya*t1482-t5*ella*t243/40+epps*trya*t1441+trxs*tr
     +ya*t1433+trxs*trys*t1509+3.0/40.0*trxs*ells*t394+ells*trxa*t100-tr
     +ys*epps*t215/18+ella*trxs*t144/4+t729*t48/2+t722*t846/2+ella*trxs*
     +t383/16-epps*t631-2.0/9.0*t5*trys*t211+trys*t1502/2-ella*t656/4-tr
     +xs*trya*t442/2-2.0/9.0*t5*trys*t1276-trys*trxa*t1081/2-ells*t42*t8
     +99/10-ells*trxs*t141/48+trys*trxa*t1433+ells*epps*t1436/2+epps*try
     +a*t1444-trxs*trys*t341/60-epps*trya*t1465+trxs*ells*t520/12+trxs*t
     +rys*epps*t618-ella*trxs*t88/4-ells*t42*t426/6-t5*trys*t7*t872/9-ep
     +ps*ella*t160/2+trys*t1090/9-t5*t134/6
      SX = SX+T731*T85/3-TRYS*T1528/2-ELLS
     +*trxa*t446+trxs*trys*trxa*t298/6-epps*t481/2-ells*t7*t1235/9+trys*
     +trxa*t1071/4+ella*t13*t1235/6-trya*t973/16+ells*trxa*t449/3+ells*t
     +42*t371/10-trxs*ella*t449/6-trya*t427/6-trys*t1206/2-trxs*ella*t44
     +6/2+epps*t1458-trya*t446-5.0/24.0*ells*trxa*t934-ells*t630/4-trya*
     +t893/6+ella*trxa*xhalb1*t140/16-t5*ella*t1167/40-3.0/4.0*ells*trxs
     +*t1294+3.0/40.0*ells*trxa*t1168-ella*t6*t1235/6+epps*ells*t7*t148/
     +6+t5*trys*t7*t125/9+trxs*trys*t452/6-epps*ella*t13*t148/2+trys*trx
     +a*t442/2-trxs*trys*trxa*t1192/6+trxs*t701/3+trys*t1199/24-ella*t13
     +*t229/6+t4*t9/9-3.0/40.0*trys*t1223+epps*ella*t6*t148/2-t5*t1220/6
     ++2.0/3.0*trxs*trya*t192+epps*t265/18-t2*t16/3+ells*t7*t229/9+3.0/4
     +0.0*trys*t394+t729*t20/3-2.0/9.0*t1220+epps*t262/2+ella*eppa*xhalb
     +1*t148/2+t731*t30/3-2.0/3.0*trxs*trya*t416-ells*t611/4+trxs*ella*t
     +156/20-ells*t14*t229/9-trys*t629/12-ells*trxa*t383/16-ells*t704/2-
     +trxs*trya*trxa*xhalb1*t872/3+epps*t349/6-epps*ella*t248/6-t2*t20/3
     ++3.0/4.0*epps*ells*t224-trya*t449/3+trys*trxs*t362/60+t145/4+ells*
     +t42*t95/6-epps*trya*trxa*t441/2+trys*t271/48-trya*t1168/10-t729*t5
     +8/2-t5*ells*t341/60-2.0/9.0*t134+t5*ells*t264/12-trys*trxa*epps*t1
     +473-epps*ells*t315/12+t5*trya*t298/6-t332/6
      SX = SX+trya*t383/16-t5*t1196/
     +3-t131/4-3.0/40.0*trys*trxa*t200-t1196/4+trya*t909/6-epps*trya*t24
     +7/3+trys*t780/4+3.0/4.0*ells*trxs*t163-trya*t248/6-epps*trya*trxs*
     +t601/2-ells*trxs*t167/12
      SX = SX+ELLS*T637/2-ELLA*T1465/4-TRYA*T947/3-TRYA
     +*t1000/10+epps*trya*trxs*t1071/2+ells*trxa*t947/3+trxs*trya*t1081/
     +2+trys*epps*t1090/18+epps*trya*t779+t5*trya*t1120/6-trxs*ella*t100
     +0/20+trya*t100+epps*trya*trxs*t220/2-2.0/3.0*trys*trxs*t133+5.0/24
     +.0*ells*trxa*t119-ells*t1528/4-trxs*trya*trxa*xhalb2*t125/3+trys*t
     +rxa*t220/4+t5*ella*t200/40+t5*ella*t118/24-t631/3-2.0/3.0*trxs*t11
     +11+ella*t626/2-ells*trxa*t798/16-ells*trxa*t939/4+ells*trxa*t973/1
     +6+trya*t798/16-trya*trxs*t243/10-trys*trxs*t215/36+ells*t607/4-trx
     +s*t252/4-trxs*t894/12-ells*t42*t155/10+trys*t1436/12-epps*trya*trx
     +s*t1026/2-ella*trxs*t1022/4-ells*trxs*t1031/48+epps*trya*t1007/3-t
     +rxs*ells*t1038/12+epps*trya*t113+epps*trya*t118/3-2.0/3.0*ells*trx
     +s*t294+trxs*trys*t305/3+ella*t660/4+t1176/4+t5*t260/6-trya*t939+tr
     +xs*trya*trxa*xhalb2*t872/3+ella*t623/4-trya*t934/6+trxs*t349/12+t5
     +*t1176/3-epps*trya*t159-trys*t315/24-t2*t737/2-t23*t1335/3-epps*t8
     +94/6-trys*t611/2-t806*t63/9-epps*t1115/18+trya*t859+trys*t224/2+t8
     +06*t81/9+trxs*trys*epps*t452/6-epps*t1111/2+t722*t1335/3-trys*t167
     +/12+trys*t163/4-epps*t1106/12+trys*t607/2+trya*t201/10-trya*t1473/
     +4-ella*eppa*xhalb1*t816/2+t793/3+2.0/9.0*t801-trxs*t1115/36+3.0/4.
     +0*epps*t1102-t5*t131/3+trya*t105/3-epps*ella*t6*t816/2
      SX = SX+t5*trys*t14
     +*t872/9-ells*t1098/4+t5*t801/6+t5*t435/4-trxs*t438/12-trys*t215/9+
     +ella*t192/4-t729*t38/3-t4*t85/3+epps*trya*trxa*t862/2-trxs*trys*t8
     +9/3+ella*t917/4+t722*t35+trxs*ella*t100/2-trys*t141/48+t23*t26+trx
     +s*ella*t105/6+ella*eppa*xhalb2*t816/2
      SX = SX+ELLS*TRXS*T109/4-ELLS*T89/4+
     +2.0/9.0*t260-ella*eppa*xhalb2*t148/2-trya*t376-trxs*t481/2+trxs*t8
     +89/2+t262/3-t2*t914+epps*ella*t114/2-trya*t943-trya*t905+trya*t96/
     +6+trys*t520/12+ells*epps*t14*t816/6+t265/9+t2*t58/2+epps*ella*t119
     +/6+2.0/3.0*trxs*t262-trya*t574/16+trya*t156/10-t1111/3-t62*t81/9+e
     +pps*t838/6-3.0/40.0*ells*trxa*t201+ells*trxa*t114/4-trys*t1105/24+
     +trya*t114+trys*t785/2-trys*t294/2+t5*t793/4-trys*t508/4+2.0/9.0*t5
     +*trys*t109-5.0/24.0*trxs*t359-t1115/9-trxs*t906/2-ella*t1060/4+trx
     +s*t346/2+trya*t119/6-t269/3-epps*t332/6-trya*t160+ells*t1049/4
      SX = SX/(AELL+1.D-30)
C
      t1 = 2*alpha
      t2 = cos(t1)
      t3 = t2**2
      t4 = xhalb1**2
      t6 = trya*trys*t4
      t7 = epps*t6
      t10 = t3*t2
      t11 = trys**2
      t12 = xhalb2**2
      t13 = t12*xhalb2
      t14 = t11*t13
      t15 = trxs*t14
      t19 = cos(4*alpha)
      t21 = ells*trya*t4
      t24 = cos(alpha)
      t25 = t24**2
      t26 = ells**2
      t27 = t26*t13
      t28 = epps*t27
      t31 = t25*t24
      t33 = ells*ella*t12
      t37 = ella*trya*xhalb2
      t40 = t4*xhalb1
      t41 = t11*t40
      t42 = trxs*t41
      t45 = trxs*t6
      t48 = t26*t40
      t49 = epps*t48
      t53 = ells*ella*t4
      t54 = epps*t53
      t58 = ella**2
      t59 = t58*xhalb1
      t60 = epps*t59
      t64 = ella*trya*xhalb1
      t67 = t58*xhalb2
      t68 = trxs*t67
      t72 = ells*trya*t12
      t76 = t26*trxs*t13
      t79 = epps*t14
      t83 = t26*trxs*t40
      t86 = trya**2
      t87 = t86*xhalb1
      t88 = trxs*t87
      t91 = epps*t67
      t94 = t86*xhalb2
      t95 = epps*t94
      t100 = trxs*t94
      t103 = epps*t87
      t107 = trya*trys*t12
      t108 = epps*t107
      t119 = epps*t41
      t123 = cos(3*alpha)
      t124 = t12*t123
      t130 = sin(t1)**2
      t131 = t130*t2
      t141 = t4*t2
      t146 = xhalb2*t2
      t156 = trya*t146
      t162 = cos(5*alpha)
      t163 = xhalb1*t162
      t164 = trya*t163
      t168 = xhalb1*t24
      t169 = trya*t168
      t173 = xhalb1*t123
      t178 = xhalb2*t162
      t183 = xhalb2*t24
      t184 = ella*t183
      t187 = xhalb2*t123
      t188 = ella*t187
      t196 = t12*t131
      t201 = trxs*t107
      t206 = trxs*t59
      t209 = trya*t124
      t213 = t12*t24
      t214 = trya*t213
      t218 = trya*t188
      t221 = trya*t184
      t224 = t4*t24
      t225 = trya*t224
      t230 = sin(alpha)**2
      t231 = t230*t24
      t242 = t12*t2
      t243 = trya*t242
      t251 = t40*t2
      t261 = trys*t124
      t267 = t4*t162
      t274 = t4*t123
      t275 = trya*t274
      t279 = t4*t131
      t286 = trxa*t242
      t296 = t13*t123
      t298 = ells*trys*t296
      t301 = t13*t24
      t303 = ells*trys*t301
      t322 = trys*t274
      t329 = trya*t141
      t334 = ella*t169
      t337 = trys*t213
      t338 = ella*t337
      t348 = xhalb1*t2
      t349 = trya*t348
      t361 = t12*t162
      t362 = trys*t361
      t367 = t40*t162
      t372 = t40*t24
      t374 = ells*trys*t372
      t377 = t40*t123
      t379 = ells*trys*t377
      t386 = trys*t224
      t387 = ella*t386
      t402 = ella*t261
      t410 = t13*t2
      t427 = ella*t322
      t439 = ella*t173
      t440 = trya*t439
      t450 = trys*t267
      t463 = t13*t162
      t499 = trya*t178
      t553 = cos(3*beta)
      t554 = t12*t553
      t555 = trya*t554
      t559 = cos(beta)
      t560 = t4*t559
      t561 = trys*t560
      t564 = xhalb1*t559
      t565 = trya*t564
      t566 = ella*t565
      t568 = t40*t553
      t571 = t40*t559
      t575 = t13*t559
      t580 = t13*t553
      t584 = cos(5*beta)
      t585 = t13*t584
      t588 = xhalb2*t553
      t592 = xhalb2*t584
      t595 = 2*beta
      t596 = cos(t595)
      t597 = xhalb2*t596
      t598 = trya*t597
      t602 = trys*t554
      t603 = ella*t602
      t607 = cos(4*beta)
      t613 = t4*t584
      t618 = t12*t559
      t619 = trya*t618
      t632 = xhalb2*t559
      t633 = trya*t632
      t634 = ella*t633
      t637 = trya*t588
      t638 = ella*t637
      t641 = trya*t592
      t647 = ells*trys*t580
      t650 = t12*t596
      t655 = trya*t650
      t659 = ells*trys*t575
      t667 = trxa*t141
      t670 = t4*t553
      t671 = trya*t670
      t678 = trya*t560
      t683 = trys*t670
      t684 = ella*t683
      t687 = ella*t561
      t692 = t596**2
      t693 = t692*t596
      t700 = t559**2
      t701 = t700*t559
      t704 = xhalb1*t584
      t708 = trys*t618
      t709 = ella*t708
      t729 = t4*t596
      t742 = t12*t584
      t743 = trys*t742
      t747 = t40*t596
      t754 = ells*trys*t571
      t759 = ells*trys*t568
      t766 = xhalb1*t596
      t767 = trya*t766
      t770 = t40*t584
      t786 = trys*t613
      t793 = trxa*t729
      t797 = xhalb1*t553
      t821 = epps*t33
      t824 = trya*t797
      t825 = ella*t824
      t828 = trya*t704
      t835 = trxa*t650
      t842 = trya*t729
      t879 = t13*t596
      t893 = sin(t595)**2
      t894 = t893*t596
      t895 = t4*t894
      t911 = t12*t894
      t993 = sin(beta)**2
      t994 = t993*t559
      SY = SY+TRYS*TRXA*ELLA*T163/10-TRXS*TRYS*ELLS*T463/60-TRYS*TRXA*EL
     +la*t588/6-ells*trxs*ella*t650/4+ells*trxs*trya*t361/40+ells*trxa*t
     +rya*t183+ells*trxa*trya*t187/6-ells*trys*t13*t19/48-trys*trxa*trya
     +*xhalb1*t131/3-ells*ella*t4*t231/6-trys*trxa*ells*t274/6+trys*trxa
     +*ells*t267/10-trys*trxa*ells*t224+ella*trys*t242/4+trys*trxa*ells*
     +t213+trys*trxa*ella*t592/10-trxs*t11*t40*t131/9+ells*ella*t12*t231
     +/6-2.0/9.0*trxs*t11*t251-ells*trxa*ella*t597/2+ells*trxa*ella*xhal
     +b2*t607/8-t647/9-trya*t743/40-trys*trxa*ella*t632-t31*t48/9-t3*t11
     +9/6+t3*t107/4+trys*trxa*trya*xhalb1*t894/3-t11*t377/36+trya*t362/4
     +0+t86*t797/6-trya*t708/4-t700*t821/2-ella*t349/4-2.0/9.0*t26*t575-
     +t692*t14/6-t86*t588/6+t693*t88/3-t693*t201/3+ells*ella*t4*t994/6-e
     +lla*t598/4+ells*trys*t13*t607/48+epps*t825/2+t700*t54/2-trxs*t659-
     +t3*t7/2+t86*t704/10+t701*t48/9-t25*t54/2+t25*t28/6-3.0/16.0*t19*t6
     +4-trxs*t379/12-t25*t49/6+t3*t79/6-trxs*t374-t607*t68/16-t25*t60/2-
     +t10*t45/3-epps*t303/6-trys*trxa*ells*t361/10-t700*t91/2+trxs*trys*
     +t842/3+ells*ella*t213/3+trxs*ella*t786/40+t700*t60/2+ells*trys*t41
     +0/4-trxs*t334+ells*trxa*t565+ells*trxa*t824/6-ells*ella*t618/3-try
     +a*t322/8-trxs*trys*trya*t911/6-t693*t15/9+t31*t27/9-t10*t88/3-ells
     +*trxa*t828/10-epps*ells*t678/4+t3*t95/2-t3*t103/2+t3*t108/2-t338/2
     ++2.0/3.0*trys*trxa*t156+2.0/9.0*trxs*t11*t747
      SY = SY+T692*T119/6-TRYS*TRX
     +a*ella*t178/10+t11*trxa*t895/6-ells*trxs*t619+t607*t206/16-ells*t6
     +55/4-t11*t835/3+t692*t103/2-ells*trxs*ella*t141/4+ells*trxa*ella*t
     +146/2-t26*trxs*t879/6+2.0/3.0*trys*trxa*t767-t607*t83/48-ells*trys
     +*t879/4+ells*ella*t560/3-trxs*t634+t709/2+epps*t709/4-ells*trxa*tr
     +ya*t173/6-epps*t379/6+t26*t40*t994/9-t26*t13*t994/9+ells*trxa*t641
     +/10-t693*t100/3-epps*t687/4-ells*trxa*t633+ells*trys*t747/4-t11*t6
     +67/3-ells*trxs*trya*t742/40+ells*trxs*trya*t613/40+epps*t659/6-t68
     +7/2+3.0/16.0*t607*t64+trxs*t638/12+t25*t91/2+t3*t14/6+t659/3+ells*
     +t243/4-trys*trxa*trya*xhalb2*t894/3+t86*t178/10+trxs*t11*t13*t131/
     +9+3.0/20.0*trxs*ella*t828-trxs*trys*t655/3+epps*t387/4+t692*t41/6+
     +trys*trxa*ells*t124/6-ella*trys*t650/4+t298/9+t11*t585/60-trxs*t64
     +7/12-ells*trxa*t637/6+t26*t286/4+t607*t76/48-t31*t59/3-trys*trxa*e
     +lls*t613/10+t86*t187/6+t701*t59/3-t603/6+epps*t684/4+t11*t286/3-tr
     +xs*t427/24-t11*trxa*t911/6+2.0/9.0*trxs*t11*t410-trxs*t603/24+2.0/
     +9.0*t26*t301+t26*t793/4-t31*t53/3+trya*t261/8+t19*t83/48-t10*t42/9
     +-epps*t647/6-2.0/9.0*trxs*t11*t879-ells*t329/4-trya*t450/40-trxs*t
     +218/12+trys*trxa*ells*t560+epps*t759/6-epps*t754/6+trxs*t221+trxs*
     +trys*ells*t585/60+t10*t100/3+t11*t793/3+trxs*t440/12-trxs*t11*t13*
     +t894/9+t693*t45/3+ells*trxs*t678+t387/2-t19*t21/16
      SY = SY+trya*t337/4-t60
     +7*t72/16+trxs*t298/12+trys*trxa*ells*t670/6+epps*ells*t671/4-epps*
     +t566/2+ells*trxs*t671/24
      SY = SY+ELLA*T156/4+EPPS*T334/2-T11*T575/6-EPPS*T
     +338/4-trya*t386/4+epps*t298/6-t11*t580/36+trya*t561/4+epps*t634/2-
     +t692*t108/2+trxs*t566+ells*trxs*ella*t729/4-t11*t770/60-3.0/20.0*t
     +rxs*ella*t641-epps*ells*t555/4+epps*ells*t619/4+t607*t21/16-t11*t3
     +72/6-ells*trxa*ella*xhalb1*t607/8-t86*t592/10-2.0/3.0*trys*trxa*t5
     +98-t692*t79/6-t26*trxa*t12*t19/16-2.0/9.0*t26*t372-ells*trxs*t555/
     +24+ella*t767/4-t26*t835/4+t11*t301/6-trys*trxa*t439/6-trys*trxa*el
     +la*t168-ells*trxa*t499/10+t31*t67/3-t700*t28/6-t26*t667/4+ells*trx
     +a*ella*t766/2-t11*t463/60+trxs*trys*trya*t196/6-t3*t41/6+trxs*trys
     +*t243/3-trya*t602/8+ells*trxa*ella*xhalb1*t19/8-ells*ella*t12*t994
     +/6+ells*trxs*ella*t242/4+t684/6+ells*trxs*t209/24+ells*trxs*t214+3
     +.0/20.0*trxs*ella*t499+ells*trys*t40*t19/48+trxs*trys*trya*t895/6+
     +t26*t13*t231/9-trxs*ella*t450/40+3.0/16.0*t19*t37+trxs*trys*ells*t
     +367/60+trxs*t759/12-3.0/20.0*trxs*ella*t164-ells*trxa*ella*t348/2-
     +trxs*t825/12-t701*t33/3+trxs*t754-t19*t76/48+t11*t571/6+t25*t821/2
     +-trxs*trys*trya*t279/6+epps*t402/4+trxs*t687+epps*ells*t225/4+t700
     +*t49/6-epps*t440/2-t303/3-t26*trxa*t4*t607/16-epps*ells*t275/4+try
     +s*trxa*ells*t742/10-ells*trxs*trya*t267/40+t11*t568/36+t26*trxs*t4
     +10/6-trys*trxa*ells*t618-trys*trxa*ells*t554/6
      SY = SY+trys*trxa*trya*xhal
     +b2*t131/3+t26*trxa*t12*t607/16+t26*trxa*t4*t19/16-epps*t427/4-t427
     +/6-epps*t221/2+trxs*ella*t362/40-ells*trys*t251/4+epps*t218/2-epps
     +*t638/2+trxs*t11*t40*t894/9-epps*ells*t214/4-3.0/16.0*t607*t37-t69
     +2*t107/4+trxs*t402/24+t402/6
      SY = SY+T701*T53/3-T26*T40*T231/9-ELLS*TRXA*E
     +lla*xhalb2*t19/8-2.0/3.0*trys*trxa*t349+trxs*t338-trys*trxa*ella*t
     +704/10+trys*trxa*ella*t564+trys*trxa*ella*t797/6+t19*t68/16+t693*t
     +42/9-t19*t206/16-t701*t67/3+ells*t842/4-ella*trys*t141/4-t3*t6/4+t
     +10*t201/3+trxs*t684/24-trxs*t709+t11*t367/60-trxs*trys*t329/3+t692
     +*t6/4+epps*ells*t209/4-ells*trxa*t169-t692*t95/2+t26*trxs*t747/6-t
     +rxs*t387+t692*t7/2-trxs*trys*ells*t770/60-t86*t163/10-ells*trxs*t2
     +75/24-ells*trys*t40*t607/48-t11*trxa*t279/6+trys*trxa*t184+trys*tr
     +xa*t188/6-ells*trxs*t225-ells*ella*t224/3+trya*t786/40-t26*trxs*t2
     +51/6+trya*t683/8-t86*t173/6+ells*trxa*t164/10+trxs*t303+t11*trxa*t
     +196/6-t379/9+epps*t374/6-trxs*ella*t743/40+ella*trys*t729/4+t374/3
     +-epps*t603/4+t759/9-t754/3+t11*t296/36+2.0/9.0*t26*t571+t19*t72/16
     ++t31*t33/3+t10*t15/9-t701*t27/9
      SY = SY/(AELL+1.D-30)
C
C     ELSE
C
C  VEREINFACHTE BERECHNUNG DER FLAECHE UND DES SCHWERPUNKTES,
C  FALLS KEINE TRIANGULARITAET (ALLE TRI... GLEICH 0).
C
C         A1 = 0.5*ELLS*(XHALB1*XHALB1-XHALB2*XHALB2)*(ALPHA-BETA)
C         A2 = ELLA*(XHALB1-XHALB2)*(0.5*(ALPHA+SIN(ALPHA)*
C    .        COS(ALPHA))-0.5*(BETA+SIN(BETA)*COS(BETA)))
C         A3 = EPPS*ELLA*(XHALB1-XHALB2)*(SIN(ALPHA)-SIN(BETA))
C         A4 = EPPS*ELLS*0.5*(XHALB1*XHALB1-XHALB2*XHALB2)*
C    .         (SIN(ALPHA)-SIN(BETA))
C         AELL = A1 + A2 + A3 + A4
C
C  BERECHNUNG DES SCHWERPUNKTES
C
C         SX = 1./AELL*(1./3.*(XHALB1**3-XHALB2**3)*
C    .         ELLS*((EPPS**2+1.)*INTCOS+EPPS*COS2+EPPS*(ALPHA-BETA))+
C    .         0.5*(XHALB1**2-XHALB2**2)*(ELLA*COS3+2.*EPPS*ELLA*
C    .         COS2+(EPPS**2*ELLA+EPPA*EPPS*ELLS)*INTCOS+EPPA*ELLS*
C    .         (ALPHA-BETA))+(XHALB1-XHALB2)*EPPA*ELLA*
C    .         (COS2+EPPS*INTCOS))
C         SY = 1./AELL*(1./3.*(XHALB1**3-XHALB2**3)*
C    .         ELLS**2*(INTSIN+EPPS*COSSIN)+0.5*(XHALB1**2-
C    .         XHALB2**2)*ELLA*ELLS*(COS2SIN+INTSIN+2.*EPPS*COSSIN)+
C    .         (XHALB1-XHALB2)*ELLA**2*(COS2SIN+EPPS*COSSIN))
C
C  BERECHNUNG DER SCHNITTPUNKTE
C
C         X2BETA = EPPS*XHALB2+EPPA + XHALB2 * COS(BETA)
C         Y2BETA = (ELLS*XHALB2+ELLA)*SIN(BETA)
C         X1BETA = EPPS*XHALB1+EPPA + XHALB1 * COS(BETA)
C         Y1BETA = (ELLS*XHALB1+ELLA)*SIN(BETA)
C         X2ALPHA = EPPS*XHALB2+EPPA + XHALB2 * COS(ALPHA)
C         Y2ALPHA = (ELLS*XHALB2+ELLA)*SIN(ALPHA)
C         X1ALPHA = EPPS*XHALB1+EPPA + XHALB1 * COS(ALPHA)
C         Y1ALPHA = (ELLS*XHALB1+ELLA)*SIN(ALPHA)
C     ENDIF
C
C  BERECHNUNG DER SCHNITTPUNKTE
C
1000  CONTINUE
      IF (IFLAG.EQ.1) GOTO 2000
C
          X2BETA = EPPS*XHALB2+EPPA+XHALB2*COS(BETA)+
     .             (TRXS*XHALB2+TRXA)*COS(2.*BETA)
          Y2BETA = (ELLS*XHALB2+ELLA)*SIN(BETA)+
     .             (TRYS*XHALB2+TRYA)*SIN(2.*BETA)
          X1BETA = EPPS*XHALB1+EPPA + XHALB1 * COS(BETA)+
     .             (TRXS*XHALB1+TRXA)*COS(2.*BETA)
          Y1BETA = (ELLS*XHALB1+ELLA)*SIN(BETA)+
     .             (TRYS*XHALB1+TRYA)*SIN(2.*BETA)
          X2ALPHA = EPPS*XHALB2+EPPA + XHALB2 * COS(ALPHA)+
     .              (TRXS*XHALB2+TRXA)*COS(2.*ALPHA)
          Y2ALPHA = (ELLS*XHALB2+ELLA)*SIN(ALPHA)+
     .              (TRYS*XHALB2+TRYA)*SIN(2.*ALPHA)
          X1ALPHA = EPPS*XHALB1+EPPA + XHALB1 * COS(ALPHA)+
     .              (TRXS*XHALB1+TRXA)*COS(2.*ALPHA)
          Y1ALPHA = (ELLS*XHALB1+ELLA)*SIN(ALPHA)+
     .              (TRYS*XHALB1+TRYA)*SIN(2.*ALPHA)
C
2000  CONTINUE
      RETURN
      END
C ===== SOURCE: arpoly.f
C
C
      SUBROUTINE ARPOLY(X,Y,N,M1,MI,ME,AREA2,XCOM,YCOM)
C
C  THIS SUBROUTINE EVALUATES SOME QUANTITIES OF M POLYGONS OF LENGTH N
C  GIVEN BY THE X COORDINATES X(J,I)
C  AND      THE Y COORDINATES Y(J,I), J=MI,ME, I=1,N
C  M1 IS THE LEADING DIMENSION IN THE ARRAYS X AND Y  AS SPECIFIED IN
C     THE CALLING PROGRAM.
C  OUTPUT:
C  AREA2(J,I) = AREA OF POLYGON SEGMENT J,I,J+1,I+1
C  XCOM(J,I)  = X COORDINATE OF CENTER OF MASS OF SEGMENT J,I,J+1,I+1
C  YCOM(J,I)  = Y COORDINATE OF CENTER OF MASS OF SEGMENT J,I,J+1,I+1
C
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, M1, MI, ME
      REAL(DP), INTENT(IN) :: X(M1,*), Y(M1,*)
      REAL(DP), INTENT(OUT) :: AREA2(M1,*), XCOM(*), YCOM(*)
      REAL(DP) :: A1, B1, C1, A2, B2, C2, RHO1, RHO2, AREAD1, AREAD2,
     .          XCD1, YCD1, XCD2, YCD2
      INTEGER :: IP1, J, JP1, I, NM1, MEM1, IN

      MEM1=ME-1
      NM1=N-1
      DO 10 J=MI,MEM1
        JP1=J+1
        DO 20 I=1,NM1
          IP1=I+1
          IN=J+(I-1)*ME
C  AREA OF CELLS
          A1=SQRT((X(J,I)-X(J,IP1))**2+(Y(J,I)-Y(J,IP1))**2)
          B1=SQRT((X(JP1,IP1)-X(J,IP1))**2+(Y(JP1,IP1)-Y(J,IP1))**2)
          C1=SQRT((X(J,I)-X(JP1,IP1))**2+(Y(J,I)-Y(JP1,IP1))**2)
          RHO1=0.5*(A1+B1+C1)
          A2=SQRT((X(J,I)-X(JP1,I))**2+(Y(J,I)-Y(JP1,I))**2)
          B2=SQRT((X(JP1,IP1)-X(JP1,I))**2+(Y(JP1,IP1)-Y(JP1,I))**2)
          C2=C1
          RHO2=0.5*(A2+B2+C2)
          AREAD1=SQRT(MAX(0._DP,RHO1*(RHO1-A1)*(RHO1-B1)*(RHO1-C1)))
          AREAD2=SQRT(MAX(0._DP,RHO2*(RHO2-A2)*(RHO2-B2)*(RHO2-C2)))
          AREA2(J,I)=AREAD1+AREAD2
C  CENTER OF MASS
          XCD1=(X(J,I)+X(JP1,IP1)+X(J,IP1))/3.
          YCD1=(Y(J,I)+Y(JP1,IP1)+Y(J,IP1))/3.
          XCD2=(X(J,I)+X(JP1,I)+X(JP1,IP1))/3.
          YCD2=(Y(J,I)+Y(JP1,I)+Y(JP1,IP1))/3.
          XCOM(IN)=(XCD1*AREAD1 + XCD2*AREAD2) / (AREA2(J,I)+1.E-60_DP)
          YCOM(IN)=(YCD1*AREAD1 + YCD2*AREAD2) / (AREA2(J,I)+1.E-60_DP)
20      CONTINUE
10    CONTINUE
      RETURN
      END
C ===== SOURCE: arquad.f
C
C
      FUNCTION ARQUAD(X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
C  THIS FUNCTION EVALUATES THE AREA OF A QUADRANGLE GIVEN BY
C  THE COORDINATES (X1,Y1), (X2,Y2), (X3,Y3) AND (X4,Y4)
C  THIS SUBROUTINE IN A REDUCED FORM OF SUBR. APOLYG (FORMERLY: POLYGM)
C
C  AREA OF CELLS
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, Y1, X2, Y2, X3, Y3, X4, Y4
      REAL(DP) :: A1, B1, C1, A2, B2, C2, RHO1, RHO2, AREAD1, AREAD2,
     .          ARQUAD

      A1=SQRT((X1-X4)**2+(Y1-Y4)**2)
      B1=SQRT((X3-X4)**2+(Y3-Y4)**2)
      C1=SQRT((X1-X3)**2+(Y1-Y3)**2)
      RHO1=0.5*(A1+B1+C1)
      A2=SQRT((X1-X2)**2+(Y1-Y2)**2)
      B2=SQRT((X3-X2)**2+(Y3-Y2)**2)
      C2=C1
      RHO2=0.5*(A2+B2+C2)
      AREAD1=SQRT(MAX(0._DP,RHO1*(RHO1-A1)*(RHO1-B1)*(RHO1-C1)))
      AREAD2=SQRT(MAX(0._DP,RHO2*(RHO2-A2)*(RHO2-B2)*(RHO2-C2)))
      ARQUAD=AREAD1+AREAD2
      RETURN
      END
C ===== SOURCE: artria.f
C
C
      FUNCTION ARTRIA (X1,Y1,X2,Y2,X3,Y3)
C
C  THIS FUNCTION EVALUATES THE AREA OF A TRIANGLE,
C  WITH A SIGN DETERMINED BY ORIENTATION
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, Y1, X2, Y2, X3, Y3
      REAL(DP) :: ARTRIA
      ARTRIA = 0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
      END
C ===== SOURCE: b_proj.f
C
C
      subroutine B_PROJ(B,EB,V,V_PARALLEL,V_PERP,PHI)
C
C  PROJECTS A VECTOR V ONTO ANOTHER VECTOR B
C  AT ENTRY: B_PROJ : NORMALIZE B TO UNIT VECTOR EB FIRST
C  AT ENTRY: B_PROJN: EB IS KNOWN TO BE THE UNIT VECTOR OF B
C  INPUT : B(3), V(3)
C  OUTPUT: V_PARALLEL, V_PERP   (V_PARALLEL CAN BE POSITIVE OR NEGATIVE)
C
      USE PRECISION
      implicit none
      REAL(DP), INTENT(INOUT) :: B(3), V(3)
      REAL(DP), INTENT(OUT) :: EB(3), V_PARALLEL, V_PERP, PHI
      REAL(DP) :: BNI, VQ

      BNI=1./SQRT(SUM(B*B)+1.D-30)
      EB=B*BNI
C
      ENTRY B_PROJN(B,EB,V,V_PARALLEL,V_PERP,PHI)
      V_PARALLEL=SUM(V*EB)
C
      VQ=SUM(V*V)
      V_PERP    =SQRT(VQ-V_PARALLEL**2)
C
C     PHI= ?
	PHI=0._DP   ! for the time being
      return
      end
C ===== SOURCE: b_proji.f


      subroutine B_PROJI(B,EB,V,V_PARALLEL,V_PERP,PHI)
C
C  INVERS OF B_PROJ: GIVEN B(3), V_PARALLEL, V_PERP AND PHI
C  CALCULATE V(3)
C  ENTRY B_PROJI : FIRST NORMALIZE B TO EB.
C  ENTRY B_PROJIN: EB IS KNOWN TO BE THE UNIT VECTOR OF B
C  INPUT:  B(3), V_PARALLEL, V_PERP, PHI
C          SIGN OF V_PARALLEL IS WITH RESPECT TO B(3)
C  OUTPUT: V(3)
C
      USE PRECISION
      implicit none
      REAL(DP), INTENT(IN) :: B(3), V_PARALLEL, V_PERP, PHI
      REAL(DP), INTENT(INOUT) :: EB(3),V(3)
      REAL(DP) :: V_P(3), E1(3), E2(3), E3(3)
      REAL(DP) :: VQ, BNI, VNI, B12I, B12

      BNI=1./SQRT(SUM(B*B)+1.D-30)
      EB=B*BNI
C
      ENTRY B_PROJIN(B,EB,V,V_PARALLEL,V_PERP,PHI)
C
      E1=EB
C
      B12=B(1)**2+B(2)**2
      IF (B12.GT.1.D-30) THEN
        B12I=1./SQRT(B12)
        E2(1)=-B(2)*B12I
        E2(2)= B(1)*B12I
        E2(3)= 0.
      ELSE
        E2(1)=1.
        E2(2)=0.
        E2(3)=0.
      ENDIF
C  E3 = E1 X E2
      E3(1)= E1(2)*E2(3)-E1(3)*E2(2)
      E3(2)=-E1(1)*E2(3)+E1(3)*E2(1)
      E3(3)= E1(1)*E2(2)-E1(2)*E2(1)
C
      V_P(1)=V_PARALLEL
      V_P(2)=COS(PHI)*V_PERP
      V_P(3)=SIN(PHI)*V_PERP
C  V = [E1,E2,E3] * V_PER
      V(1)=V_P(1)*E1(1)+V_P(2)*E2(1)+V_P(3)*E3(1)
      V(2)=V_P(1)*E1(2)+V_P(2)*E2(2)+V_P(3)*E3(2)
      V(3)=V_P(1)*E1(3)+V_P(2)*E2(3)+V_P(3)*E3(3)
      VNI=1./SQRT(SUM(V**2)+1.D-30)
      V = V * VNI
      return
      end
C ===== SOURCE: bilinear_int.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****s* INTERPOLATION/BILINEAR_INT
C NAME
C     BILINEAR_INT (D, RX, RY, Z)
C DESCRIPTION
C     Bilinear interpolation of data
C INPUTS
C     REAL*8, DIMENSION(2,2):: D  ! data to interpolate
C     REAL*8    :: RX, RY ! relative position at which interpolations should be done
C OUTPUT
C     REAL*8    :: Z  ! interpolated data
C     ******
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BILINEAR_INT (D, RX, RY, Z)
      USE PRECISION
      IMPLICIT NONE

!---- input variables
      REAL(DP), DIMENSION(2,2), INTENT(IN):: D
      REAL(DP), INTENT(IN)  :: RX, RY
!---- output variables
      REAL(DP), INTENT(OUT)  :: Z

!---- local variables
      REAL(DP)  :: Z1, Z2

      Z1 = D(1,1)*RX + D(2,1)*(1-RX)
      Z2 = D(1,2)*RX + D(2,2)*(1-RX)
      Z  = Z1*RY + Z2*(1-RY)

      RETURN
      END SUBROUTINE BILINEAR_INT
C ===== SOURCE: binsearch.f
      subroutine binsearch(xx,n,x,i)
c     ***********************************************************
c     * ermittlung des Feldindex i mit vorgegebener zahl x,     *
c     * so dass x zwischen xx(i) und xx(i+1)                    *
c     * ferner: falls x<=xx(1) i=1, und falls x>=xx(n) i=n      *
c     * *********************************************************
      use PRECISION
      implicit none

      integer, intent(in) :: n
      integer, intent(out) :: i
      real(dp), intent(in) :: xx(n), x
      integer :: bl, bm, bu

c  savest version:
c  all tests included
      entry binsearch_0(xx,n,x,i)

      bl=0
      bu=n+1

      if (xx(n).ge.xx(1)) then
! monoton increasing
        if(x.le.xx(1))then
          i=1
        else if(x.ge.xx(n))then
          i=n
        else
c  binary search
          do while (bu-bl.gt.1)
            bm=(bu+bl)*0.5
            if(x.ge.xx(bm)) then
              bl=bm
            else
              bu=bm
            endif
          enddo
          i=bl
        end if
      else
! monoton decreasing
        if(x.ge.xx(1))then
          i=1
        else if(x.le.xx(n))then
          i=n
        else
c  binary search
          do while (bu-bl.gt.1)
            bm=(bu+bl)*0.5
            if(x.le.xx(bm)) then
              bl=bm
            else
              bu=bm
            endif
          enddo
          i=bl
        end if
      end if
c
      return

c  fast version:
c  we already know: a)  xx is monotonically increasing (not decreasing)
c                   b)  x  lies between xx(1) and xx(n)
      entry binsearch_2(xx,n,x,i)

      bl=0
      bu=n+1

c  binary search
      do while (bu-bl.gt.1)
        bm=(bu+bl)*0.5
        if(x.ge.xx(bm))then
          bl=bm
        else
          bu=bm
        endif
      end do
c
      i=bl

      return
      end
C ===== SOURCE: bitget.f
C
C
      LOGICAL FUNCTION BITGET(IBIT,N1LOW,N1UP,IROW,JCOL,NBITS)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N1LOW, N1UP, IROW, JCOL, NBITS
      INTEGER, INTENT(INOUT) :: IBIT(N1LOW:N1UP,*)
      INTEGER :: JELEM, JB, ISET

      JELEM=JCOL/NBITS
      JB=MOD(JCOL,NBITS)
      IF (JB == 0) THEN
        JB=NBITS-1
      ELSE
        JELEM=JELEM+1
        JB=JB-1
      END IF

      BITGET = BTEST(IBIT(IROW,JELEM),JB)

      RETURN
      END
C ===== SOURCE: bitset.f
C
C
      SUBROUTINE BITSET(IBIT,N1LOW,N1UP,IROW,JCOL,ISET,NBITS)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N1LOW, N1UP, IROW, JCOL, ISET, NBITS
      INTEGER, INTENT(INOUT) :: IBIT(N1LOW:N1UP,*)
      INTEGER :: JELEM, JB

      JELEM=JCOL/NBITS
      JB=MOD(JCOL,NBITS)
      IF (JB == 0) THEN
        JB=NBITS-1
      ELSE
        JELEM=JELEM+1
        JB=JB-1
      END IF

      IF (ISET == 0) IBIT(IROW,JELEM) = IBCLR(IBIT(IROW,JELEM),JB)
      IF (ISET == 1) IBIT(IROW,JELEM) = IBSET(IBIT(IROW,JELEM),JB)

      RETURN
      END
C ===== SOURCE: dbl_poly.f
      subroutine dbl_poly (cf, al, pl, cou, dum, ji, je,
     .                     rcmin, rcmax, fp, ifexmn, ifexmx )

      use precision
      
      implicit none

      real(dp), intent(in) :: cf(9,9), fp(6)
      real(dp), intent(in) :: al, pl, rcmin, rcmax
      real(dp), intent(out) :: cou, dum(9) 
      integer, intent(in) :: ji, je, ifexmn, ifexmx
      real(dp) :: s01, s02, ds12, expo1, expo2, ccxm1, ccxm2,
     .            fpar1, fpar2, fpar3, extrap
      integer :: kk, jj, j, i, ii, ifex

      dum = 0._dp

      if ((ifexmn .ne. 0) .and. (al < rcmin)) then
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(<S*V>)
        S01=RCMIN
        S02=LOG(2._DP)+RCMIN
        DS12=S02-S01
        EXPO1=0.
        EXPO2=0.
        DO 1 J=1,9
          JJ=J-1
          DO 1 I=1,9
            II=I-1
            EXPO1=EXPO1+S01**II*PL**JJ*CF(I,J)
            EXPO2=EXPO2+S02**II*PL**JJ*CF(I,J)
 1      CONTINUE
        CCXM1=EXPO1
        CCXM2=EXPO2
        FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
        FPAR2=      (CCXM2-CCXM1)/DS12
        FPAR3=0.D0
C
        IFEX=5
        COU=EXTRAP(AL,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)

      else

        do jj = je, ji, -1
          dum(jj) = cf(9,jj)
          do kk = 8, 1, -1
            dum(jj) = dum(jj) * al + cf(kk,jj)
          end do 
        end do

        cou = dum(9)

        do jj = 8, 1, -1
          cou = cou * pl + dum(jj)
        end do

      end if

      return
      end subroutine dbl_poly
C ===== SOURCE: dekey.f
C
C
      SUBROUTINE DEKEY (A,NLFLD,NFIRST0,NFIRST1,NFIX,NDIM)
C
C   MODIFY A INTEGER FIELD NLFLD ACCORDING TO A CHARACTER STRING A
C
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NFIRST0, NFIRST1, NFIX, NDIM
      INTEGER, INTENT(INOUT) :: NLFLD(NFIRST0:NFIRST1,NDIM)
      CHARACTER(*), INTENT(IN) :: A
      INTEGER :: LOG, NMIN, IN, IN2, J, NMAX, IANF, LA
      CHARACTER(5) :: FORM

      DATA FORM/'(I  )'/
C
      LA=LEN(A)
      IANF=1
C
    1 IN=INDEX(A(IANF:LA),'/')
      IF (IN.EQ.0) RETURN
      WRITE (FORM(3:4),'(I2)') IN-1
      READ (A(IANF:LA),FORM) NMIN
      IN2=INDEX(A(IANF+IN:LA),' ')
      WRITE (FORM(3:4),'(I2)') IN2-1
      READ (A(IANF+IN:LA),FORM) NMAX
      IANF=IANF+IN+IN2
      LOG=0
      IF (NMIN.GT.0) LOG=1
      DO 2 J=ABS(NMIN),NMAX
2        NLFLD(NFIX,J)=LOG
      IF (IANF.GE.LA) RETURN
      GOTO 1
C
      END
C ===== SOURCE: dekeyb.f
C
C
      SUBROUTINE DEKEYB (A,NLFLD,NFIRST0,NFIRST1,NFIX,NDIM,NBITS)
C
C   MODIFY A LOGICAL FIELD NLFLD ACCORDING TO A CHARACTER STRING A
C
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NFIRST0, NFIRST1, NFIX, NDIM, NBITS
      INTEGER, INTENT(INOUT) :: NLFLD(NFIRST0:NFIRST1,NDIM)
      CHARACTER(*), INTENT(IN) ::  A
      CHARACTER(5) :: FORM
      INTEGER :: LOG, IN, NMIN, IANF, J, NMAX, IN2, LA
      DATA FORM/'(I  )'/
C
      LA=LEN(A)
      IANF=1
C
    1 IN=INDEX(A(IANF:LA),'/')
      IF (IN.EQ.0) RETURN
      WRITE (FORM(3:4),'(I2)') IN-1
      READ (A(IANF:LA),FORM) NMIN
      IN2=INDEX(A(IANF+IN:LA),' ')
      WRITE (FORM(3:4),'(I2)') IN2-1
      READ (A(IANF+IN:LA),FORM) NMAX
      IANF=IANF+IN+IN2
      LOG=0
      IF (NMIN.GT.0) LOG=1
      DO J=ABS(NMIN),NMAX
        CALL BITSET(NLFLD,NFIRST0,NFIRST1,NFIX,J,LOG,NBITS)
      ENDDO
      IF (IANF.GE.LA) RETURN
      GOTO 1
C
      END
C ===== SOURCE: deter.f
C
C
C*DK DETER
      FUNCTION DETER(A11,A21,A31,A12,A22,A32,A13,A23,A33)
C
C  DETERMINAT OF (3,3)-MATRIX A
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A11, A21, A31, A12, A22, A32, A13, A23,
     .                        A33
      REAL(DP) :: S, DETER

      S=A11*(A22*A33-A23*A32)
      S=S-A12*(A21*A33-A23*A31)
      S=S+A13*(A21*A32-A22*A31)
      DETER=S
      RETURN
      END
C ===== SOURCE: deter4x4.f
      function deter4x4 (a)
      use precision
      implicit none

      real(dp), intent(in) :: a(4,4)
      real(dp) :: deter4x4, d11, d12, d13, d14, deter

      d11 = deter(a(2,2),a(3,2),a(4,2),
     .            a(2,3),a(3,3),a(4,3),
     .            a(2,4),a(3,4),a(4,4))

      d12 = deter(a(2,1),a(3,1),a(4,1),
     .            a(2,3),a(3,3),a(4,3),
     .            a(2,4),a(3,4),a(4,4))

      d13 = deter(a(2,1),a(3,1),a(4,1),
     .            a(2,2),a(3,2),a(4,2),
     .            a(2,4),a(3,4),a(4,4))

      d14 = deter(a(2,1),a(3,1),a(4,1),
     .            a(2,2),a(3,2),a(4,2),
     .            a(2,3),a(3,3),a(4,3))

      deter4x4 = a(1,1)*d11 - a(1,2)*d12 + a(1,3)*d13 - a(1,4)*d14

      return
      end





C ===== SOURCE: fehler.f



C-----------------------------------------------------------------------
              SUBROUTINE FEHLER(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DIESES UNTERPROGRAMM UEBERPRUEFT MOEGLICHE REGELVERLETZUNGEN
C     VON AUSDRU
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: MAXOPT=15
C           : MAXIMAL ZULAESSIGE ANZAHL VON OPERATOREN IN AUSDRU

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     LOKALE VARIABLEN :
C
         INTEGER :: OANDEN
C           : ANZAHL DER OPERANDEN IN AUSDRU

         INTEGER :: OTOREN
C           : ANZAHL DER OPERATOREN IN AUSDRU


      ERROR=0
C
C     UEBERPRUEFUNG AUF GUELTIGE ZEICHEN
C
      CALL SIGNOK(AUSDRU,AKTLEN,ERROR)
      IF (ERROR .EQ. 0) THEN
C
C        UEBERPRUEFUNG AUF KORREKTE KLAMMERUNG
C
         CALL KLAMME(AUSDRU,AKTLEN,ERROR)
         IF (ERROR .EQ. 0) THEN
C
C           UEBERPRUEFUNG AUF KORREKTHEIT DER OPERANDEN
C
            CALL OPRAND(AUSDRU,AKTLEN,OANDEN,ERROR)
            IF (ERROR .EQ. 0) THEN
C
C           UEBERPRUEFUNG AUF KORREKTHEIT DER OPERATOREN
C
               CALL OPERAT(AUSDRU,AKTLEN,OTOREN,ERROR)
               IF (ERROR .EQ. 0) THEN
                  IF (OTOREN .GT. MAXOPT) THEN
C
C                    AUSDRU ENTHAELT MEHR ALS DIE MAXIMAL ZULAESSIGE
C                    OPERATORENANZAHL
C
                     ERROR=12
                  ELSEIF (OTOREN .NE. OANDEN-1) THEN
C
C                    DIE ANZAHL DER OPERATOREN IST VERSCHIEDEN VON
C                    DER ANZAHL DER OPERANDEN-1
C
                     ERROR=13
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
C     ENDE VON FEHLER
C
      END
C ===== SOURCE: ftcre.f
C
C
      SUBROUTINE FTCRE (F,C)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: F
      CHARACTER(10), INTENT(OUT) :: C
      WRITE(C,'(1P,E10.3)') F
      RETURN
      END
C ===== SOURCE: ftcri.f
C
C
      SUBROUTINE FTCRI (I,C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      CHARACTER(6), INTENT(OUT) :: C
      WRITE(C,'(1P,I6)') I
      RETURN
      END
C ===== SOURCE: headng.f
C
C
C*DK HEADNG
      SUBROUTINE HEADNG (A,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      INTEGER :: I
      CHARACTER(1) :: U(72)
      DATA U/72*'='/
      WRITE (iunout,'(1X,A)') A
      WRITE (iunout,60) (U(I),I=1,N)
60    FORMAT (1X,72A1)
      RETURN
      END
C ===== SOURCE: idez.f
C
C
      FUNCTION IDEZ(I,J,IZIF)
C
C   I IS AN INTEGER WITH  N DEZIMALS, THE FUNCTION RETURNS
C   IDEZ, THE J.TH DEZIMAL
C   J MUST BE .LE. IZIF
C      IF N .LT. IZIF AND J .GT. N  , IDEZ=0
C      IF N .GT. IZIF AND J .EQ. IZIF  , IDEZ=PART OF I ON THE LEFT OF
C      THE DEZIMAL J, INCLUDING THE J TH DEZIMAL
C      E.G.: I=1234,IZIF=3, J=1: IDEZ=4
C                           J=2: IDEZ=3
C                           J=3: IDEZ=12
C                           J=4: ERROR ,J GT IZIF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, J, IZIF
      INTEGER :: IDEZ, K, IZ, IQ, IT, IDIF

      IDIF=IZIF-J+1
      IF (IDIF.LE.0) THEN
        CALL MASAGE ('ERROR IN FUNCTION IDEZ                         ')
        CALL EXIT_OWN(1)
      ENDIF
      IZ=I
      DO 1 K=1,IDIF
        IT=10**(IZIF-K)
        IQ=IZ/IT
        IZ=IZ-IQ*IT
1     CONTINUE
      IDEZ=IQ
      RETURN
      END
C ===== SOURCE: iexp10.f
C
C
C----------------------------------------------------------------------*
C FUNCTION IEXP10
C----------------------------------------------------------------------*
      FUNCTION IEXP10(ZAHL)
C
C  FUNKTION, WELCHE DEN 10-ER-EXPONENTEN EINER ZAHL
C  BERECHNET. (AUSNAHME: IEXP10(0)=0)
C
C  EINGABE: ZAHL             REAL
C  AUSGABE: IEXP10(ZAHL)     INTEGER
C
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ZAHL
      REAL(DP) :: ZHL2
      INTEGER :: I, IEXP10

      I=0
      IF( ZAHL.EQ.0) THEN
         IEXP10=I
         RETURN
      ENDIF
      ZHL2=ABS(ZAHL)
  10  IF (ZHL2.LT.1) THEN
         I=I-1
         ZHL2=ZHL2*10
         GOTO 10
      ENDIF
  20  IF (ZHL2.GE.10) THEN
         I=I+1
         ZHL2=ZHL2/10
         GOTO 20
      ENDIF
      IEXP10=I
      RETURN
      END
C ===== SOURCE: illz.f
C
c
      FUNCTION ILLZ(N,V,IV)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: V(*)
      INTEGER, INTENT(IN) :: N, IV
      INTEGER :: J, ILLZ
      J=1
      IF (IV.LT.0) J=1-(N-1)*IV
      DO 17,ILLZ=1,N
CODER    IF (V(J).LT.0) GOTO 18    BEI REAL ODER INTEGER
CPB      IF (.NOT.V(J)) GOTO 18
         IF (V(J)) GOTO 18
   17    J=J+IV
   18 ILLZ=ILLZ-1
      END
C ===== SOURCE: indtal.f
C
C
C*DK INDTAL
      SUBROUTINE INDTAL (IND,M,NX,NY,NZ,NB)
C
C   called from subr. STATIS
C               subr. STATIS_BGK
C               subr. STATIS_COP
C   SIMILAR TO SUBR. INTTAL:
C   PROVIDE ARRAY IND(J,K), J=1,NRAD, K=1,8
C   THE VALUES IN = IND(J,K) , FOR EACH CELL J,
C   ARE THE INDICES OF THOSE BINS TO WHICH CELL J CONTRIBUTES
C   FOR EACH EIRENE CELL JE: IND(JE,1)=JE,
C
C   2ND INDEX .LE. 8, BECAUSE J MAY CONTRIBUTE TO
C
C    1     J              (FUNCTION OF X,Y,Z)
C    2     X      AVERAGE (FUNCTION OF Y,Z)
C    3     Y      AVERAGE (FUNCTION OF X,Z)
C    4     Z      AVERAGE (FUNCTION OF X,Y)
C    5     X-Y    AVERAGE (FUNCTION OF Z)
C    6     X-Z    AVERAGE (FUNCTION OF Y)
C    7     Y-Z    AVERAGE (FUNCTION OF X)
C    8     X-Y-Z  AVERAGE   --
C
C   M: LEADING DIMENSION OF IND IN CALLING PROGRAM
C  FURTHER (TO SPEED UP SUBR.STATIS)
C
C
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: M, NX, NY, NZ, NB
      INTEGER, INTENT(OUT) :: IND(M,8)
      INTEGER :: IX, IY, IZ, NIR, NIRT, NIRTP, IADD, K, IR, IJ, IK,
     .           NIT, NIPR, NIP, NITP, NS, I, NXM, NYM, NZM, IB, J,
     .           N1DEL, N2DEL
C
C NS: NUMBER OF STANDARD MESH CELLS
      NS=NX*NY*NZ*NB
C
      DO 10 I=1,M
        DO 10 J=1,8
          IND(I,J)=0
10    CONTINUE
C
C  additional cells contribute only to themselves
C
      DO 11 I=NS+1,M
        IND(I,1)=I
11    CONTINUE
C
      IF (NS.EQ.0) RETURN
C
      N1DEL=0
      IF (NY.GT.1.OR.NZ.GT.1) N1DEL=NX
      N2DEL=0
      IF (NZ.GT.1) N2DEL=NY
C
      NXM=MAX(1,NX-1)
      NYM=MAX(1,NY-1)
      NZM=MAX(1,NZ-1)

C  LOOP OVER STANDARD MESH BLOCKS
      DO 1000 IB=1,NB
      IADD=(IB-1)*NX*NY*NZ
C
      NIRTP=NX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
      DO 100 IY=1,NYM
        NIRT=NX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 101 IZ=1,NZM
          NIR=NX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 102 IX=1,NXM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            IND(K,1)=K
            IND(K,2)=NIR
            IND(K,6)=NIRT
            IND(K,8)=NIRTP
102       CONTINUE
101     CONTINUE
100   CONTINUE
C
      DO 200 IZ=1,NZM
        NIPR=NX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
        DO 201 IX=1,NXM
          NIP=IX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 202 IY=1,NYM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            IND(K,1)=K
            IND(K,3)=NIP
            IND(K,5)=NIPR
202       CONTINUE
201     CONTINUE
200   CONTINUE
C
      DO 300 IX=1,NXM
        NITP=IX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 301 IY=1,NYM
          NIT=IX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
          DO 302 IZ=1,NZM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            IND(K,1)=K
            IND(K,4)=NIT
            IND(K,7)=NITP
302       CONTINUE
301     CONTINUE
300   CONTINUE
C
1000  CONTINUE
C  teste: gibt es ein ir, fuer verschiedene j aber gleiche ind(ir,j)
C         dann wuerde naemlich in statis doppelt gezaehlt
C         falls weniger als 3 dimensionen, kann sowas vorkommen.
      DO 2000 J=1,8
2000  continue
      do 2001 ir=1,ns
        do j=1,8
          ij=ind(ir,j)
          do k=1,8
            ik=ind(ir,k)
            if (ik.ne.0.and.ik.eq.ij.and.k.ne.j) then
              if (k.gt.j) ind(ir,k)=0
              if (j.gt.k) ind(ir,j)=0
            endif
          enddo
        enddo
2001  continue
C
      do 2002 ir=1,ns
        do j=1,8
          ij=ind(ir,j)
        enddo
2002  continue

C
      RETURN
      END
C ===== SOURCE: inttal.f
C
C
C*DK INTTAL
      SUBROUTINE INTTAL (A,VOL,J,M,N,YINT,NX,NY,NZ,NB)
C
C   INTEGRATE A(J,K), K=1,N, RESULT IS YINT
C   USE VOL(K) AS WEIGHTING
C   J FIXED, (SPECIES INDEX)
C   M: LEADING DIMENSION OF A IN CALLING PROGRAM
C   A IS A VOLUME AVERAGED TALLY
C
C  IN CASE NX*NY*NZ*NB NOT ZERO, A STANDARD MESH HAS BEEN DEFINED.
C                                THEN THE VARIOUS PROJECTIONS
C                                ONTO THE 1ST,2ND AND 3RD MESH ARE
C                                COMPUTED
C                                IX=1,NX-1, IY=1,NY, IZ=1,NZ-1
C
      USE PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: J, M, N, NX, NY, NZ, NB
      REAL(DP), INTENT(INOUT) :: A(M,*)
      REAL(DP), INTENT(IN) :: VOL(*)
      REAL(DP), INTENT(OUT) :: YINT
      REAL(DP) :: VR, YR, YPR, VRTP, VRT, YRT, VTP, YTP, VT, YT, YP,
     .            VPR, A0, VP, YRTP
      INTEGER :: IX, IY, IZ, NIR, NIRTP, NIRT, NITP, NIT, NIPR, NIP,
     .           N1DEL, N2DEL, KX, KY, KZ, IB, K, NXM, NYM, NZM, IADD,
     .           KB, NS
C
      N1DEL=0
      IF (NY.GT.1.OR.NZ.GT.1) N1DEL=NX
      N2DEL=0
      IF (NZ.GT.1) N2DEL=NY
      NXM=MAX(1,NX-1)
      NYM=MAX(1,NY-1)
      NZM=MAX(1,NZ-1)
      NS=NX*NY*NZ*NB
C
C  LOOP OVER ALL CELLS
      YINT=0.
      DO 5 KB=1,NB
        IADD=(KB-1)*NX*NY*NZ
        DO 5 KX=1,NXM
        DO 5 KY=1,NYM
        DO 5 KZ=1,NZM
          K=KX+((KY-1)+(KZ-1)*N2DEL)*N1DEL+IADD
          YINT=YINT+VOL(K)*A(J,K)
5     CONTINUE
      DO 6 K=NS+1,N
        YINT=YINT+VOL(K)*A(J,K)
6     CONTINUE
C
C  LOOP OVER STANDARD MESH BLOCKS
      IF (NS.EQ.0) RETURN
C
      DO 1000 IB=1,NB
      IADD=(IB-1)*NX*NY*NZ
C
C  INTEGRATE OVER RADIAL CO-ORDINATE: YR
C  INTEGRATE OVER RADIAL AND TOROIDAL CO-ORDINATE: YRT
C  INTEGRATE OVER ALL THREE CO-ORDINATES: YRTP
      YRTP=0.
      VRTP=1.D-60
      NIRTP=NX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
      DO 100 IY=1,NYM
        YRT=0.
        VRT=1.D-60
        NIRT=NX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 101 IZ=1,NZM
          YR=0.
          VR=1.D-60
          NIR=NX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 102 IX=1,NXM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YR=YR+VOL(K)*A(J,K)
            YRT=YRT+VOL(K)*A(J,K)
            YRTP=YRTP+VOL(K)*A(J,K)
            VR=VR+VOL(K)
            VRT=VRT+VOL(K)
            VRTP=VRTP+VOL(K)
102       CONTINUE
          A(J,NIR)=YR/VR
101     CONTINUE
        A(J,NIRT)=YRT/VRT
100   CONTINUE
      A(J,NIRTP)=YRTP/VRTP
C  INTEGRATE OVER POLOIDAL CO-ORDINATE: YP
C  INTEGRATE OVER POLOIDAL AND RADIAL CO-ORDINATE: YPR
      DO 200 IZ=1,NZM
        YPR=0.
        VPR=1.D-60
        NIPR=NX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
        DO 201 IX=1,NXM
          YP=0.
          VP=1.D-60
          NIP=IX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 202 IY=1,NYM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YP=YP+VOL(K)*A(J,K)
            YPR=YPR+VOL(K)*A(J,K)
            VP=VP+VOL(K)
            VPR=VPR+VOL(K)
202       CONTINUE
          A(J,NIP)=YP/VP
201     CONTINUE
        A0=A(J,NIPR)
        A(J,NIPR)=YPR/VPR
200   CONTINUE
C  INTEGRATE OVER TOROIDAL CO-ORDINATE: YT
C  INTEGRATE OVER TOROIDAL AND POLOIDAL CO-ORDINATE: YTP
      DO 300 IX=1,NXM
        YTP=0.
        VTP=1.D-60
        NITP=IX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 301 IY=1,NYM
          YT=0.
          VT=1.D-60
          NIT=IX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
          DO 302 IZ=1,NZM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YT=YT+VOL(K)*A(J,K)
            YTP=YTP+VOL(K)*A(J,K)
            VT=VT+VOL(K)
            VTP=VTP+VOL(K)
302       CONTINUE
          A(J,NIT)=YT/VT
301     CONTINUE
        A(J,NITP)=YTP/VTP
300   CONTINUE
1000  CONTINUE
      RETURN
      END
C ===== SOURCE: intvol.f
C
C
C*DK INTVOL
      SUBROUTINE INTVOL (A,J,M,N,YINT,NX,NY,NZ,NB)
C
C   SIMILAR TO INTTAL
C   INTEGRATE A(J,K), K=1,N, RESULT IS YINT
C   USE 1. AS WEIGHTING. PUT TOTALS RATHER THAN MEANS ON AVERAGE-TALLY
C   LOCATIONS
C   J FIXED, (SPECIES INDEX)
C   M: LEADING DIMENSION OF A IN CALLING PROGRAM
C   A IS A TOTAL TALLY (EG. CELL VOLUME)
C
C  IN CASE NX*NY*NZ*NB NOT ZERO, A STANDARD MESH HAS BEEN DEFINED.
C                                THEN THE VARIOUS PROJECTIONS
C                                ONTO THE 1ST,2ND AND 3RD MESH ARE
C                                COMPUTED
C                                IX=1,NX-1, IY=1,NY, IZ=1,NZ-1
C
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: J, M, N, NX, NY, NZ, NB
      REAL(DP), INTENT(INOUT) :: A(M,*)
      REAL(DP), INTENT(OUT) :: YINT
      REAL(DP) :: YR, YRT, YRTP, YT, YTP, YP, YPR
      INTEGER :: IY, IZ, IB, K, NIR, NIRTP, NIT, NITP, NIP, IX, NIPR,
     .           KZ, N1DEL, N2DEL, IADD, KB, KX, KY, NXM, NIRT, NS,
     .           NYM, NZM
C
      N1DEL=0
      IF (NY.GT.1.OR.NZ.GT.1) N1DEL=NX
      N2DEL=0
      IF (NZ.GT.1) N2DEL=NY
      NXM=MAX(1,NX-1)
      NYM=MAX(1,NY-1)
      NZM=MAX(1,NZ-1)
      NS=NX*NY*NZ*NB
C
C  LOOP OVER ALL CELLS
      YINT=0.
      DO 5 KB=1,NB
        IADD=(KB-1)*NX*NY*NZ
        DO 5 KX=1,NXM
        DO 5 KY=1,NYM
        DO 5 KZ=1,NZM
          K=KX+((KY-1)+(KZ-1)*N2DEL)*N1DEL+IADD
          YINT=YINT+A(J,K)
5     CONTINUE
      DO 6 K=NS+1,N
        YINT=YINT+A(J,K)
6     CONTINUE
C
C  LOOP OVER STANDARD MESH BLOCKS
      IF (NS.EQ.0) RETURN
C
      DO 1000 IB=1,NB
      IADD=(IB-1)*NX*NY*NZ
C
C  INTEGRATE OVER RADIAL CO-ORDINATE: YR
C  INTEGRATE OVER RADIAL AND TOROIDAL CO-ORDINATE: YRT
C  INTEGRATE OVER ALL THREE CO-ORDINATES: YRTP
      YRTP=0.
      NIRTP=NX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
      DO 100 IY=1,NYM
        YRT=0.
        NIRT=NX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 101 IZ=1,NZM
          YR=0.
          NIR=NX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 102 IX=1,NXM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YR=YR+A(J,K)
            YRT=YRT+A(J,K)
            YRTP=YRTP+A(J,K)
102       CONTINUE
          A(J,NIR)=YR
101     CONTINUE
        A(J,NIRT)=YRT
100   CONTINUE
      A(J,NIRTP)=YRTP
C  INTEGRATE OVER POLOIDAL CO-ORDINATE: YP
C  INTEGRATE OVER POLOIDAL AND RADIAL CO-ORDINATE: YPR
      DO 200 IZ=1,NZM
        YPR=0.
        NIPR=NX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
        DO 201 IX=1,NXM
          YP=0.
          NIP=IX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 202 IY=1,NYM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YP=YP+A(J,K)
            YPR=YPR+A(J,K)
202       CONTINUE
          A(J,NIP)=YP
201     CONTINUE
        A(J,NIPR)=YPR
200   CONTINUE
C  INTEGRATE OVER TOROIDAL CO-ORDINATE: YT
C  INTEGRATE OVER TOROIDAL AND POLOIDAL CO-ORDINATE: YTP
      DO 300 IX=1,NXM
        YTP=0.
        NITP=IX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 301 IY=1,NYM
          YT=0.
          NIT=IX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
          DO 302 IZ=1,NZM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YT=YT+A(J,K)
            YTP=YTP+A(J,K)
302       CONTINUE
          A(J,NIT)=YT
301     CONTINUE
        A(J,NITP)=YTP
300   CONTINUE
1000  CONTINUE
      RETURN
      END
C ===== SOURCE: invert.f
C
C
      SUBROUTINE INVERT(A,AI)

      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(3,3)
      REAL(DP), INTENT(OUT) :: AI(3,3)
      REAL(DP) :: DET
C
      DET = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1)
     .    + A(1,3) * A(2,1) * A(3,2) - A(1,3) * A(2,2) * A(3,1)
     .    - A(1,2) * A(2,1) * A(3,3) - A(1,1) * A(2,3) * A(3,2)
C
C     WRITE (*,*) DET
C
      AI(1,1) = (A(2,2) * A(3,3) - A(2,3) * A(3,2)) / DET
      AI(1,2) = (A(1,3) * A(3,2) - A(1,2) * A(3,3)) / DET
      AI(1,3) = (A(1,2) * A(2,3) - A(2,2) * A(1,3)) / DET
      AI(2,1) = (A(2,3) * A(3,1) - A(2,1) * A(3,3)) / DET
      AI(2,2) = (A(1,1) * A(3,3) - A(1,3) * A(3,1)) / DET
      AI(2,3) = (A(1,3) * A(2,1) - A(1,1) * A(2,3)) / DET
      AI(3,1) = (A(2,1) * A(3,2) - A(2,2) * A(3,1)) / DET
      AI(3,2) = (A(1,2) * A(3,1) - A(1,1) * A(3,2)) / DET
      AI(3,3) = (A(1,1) * A(2,2) - A(1,2) * A(2,1)) / DET
C
      RETURN
      END
C ===== SOURCE: klamme.f



C-----------------------------------------------------------------------
             SUBROUTINE KLAMME(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFUNG AUF KORREKTE KLAMMERUNG
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: KLAMAX=3
C           : MAXIMAL ZULAESSIGE ANZAHL VON KLAMMERSCHACHTELUNGEN

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     LOKALE VARIABLEN :
C
         INTEGER :: ANZAHL
C           : VARIABLE, DIE JE NACH IHREM WERT ANGIBT
C             = 0  ANZAHL DER OEFFNENDEN KLAMMERN IST GLEICH DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN
C             > 0  ANZAHL DER OEFFNENDEN KLAMMERN IST GROESSER DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN
C             < 0  ANZAHL DER OEFFNENDEN KLAMMERN IST KLEINER DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN

         INTEGER :: MINKLA
C           : MINIMUM VON ANZAHL

         INTEGER :: MAXKLA
C           : MAXIMUM VON ANZAHL

         INTEGER :: KLAMON
C           : POSITION DER ZULETZ GEFUNDENEN OEFFNENDEN KLAMMER

         INTEGER :: KLAMOF
C           : POSITION DER ZULETZ GEFUNDENEN SCHLIESSENDEN KLAMMER

C
C     HILFSVARIABLE :
C
         INTEGER :: I


      ANZAHL=0
      MAXKLA=0
      MINKLA=0
      KLAMON=0
      KLAMOF=-1

      DO 10, I=1,AKTLEN
         IF (AUSDRU(I:I) .EQ. '(' .AND. ERROR .EQ. 0 ) THEN
            ANZAHL=ANZAHL+1
            MAXKLA=MAX(MAXKLA,ANZAHL)
            KLAMON=I

            IF (KLAMON .EQ. KLAMOF+1) THEN
C
C              ZWISCHEN ZWEI KLAMMERAUSDRUECKEN STEHT KEIN OPERATOR
C
               ERROR=2
            ENDIF
         ELSEIF (AUSDRU(I:I) .EQ. ')' .AND.  ERROR .EQ. 0 ) THEN
            ANZAHL=ANZAHL-1
            MINKLA=MIN(MINKLA,ANZAHL)
            KLAMOF=I

            IF (KLAMOF .EQ. KLAMON+1) THEN
C
C              DER KLAMMERAUSDRUCK IST LEER
C
               ERROR=3
            ENDIF
         ENDIF
10       CONTINUE

      IF ((ANZAHL .NE. 0  .OR. MINKLA .LT. 0) .AND. ERROR .EQ. 0) THEN
C
C        FALSCHE KLAMMERSETZUNG
C
         ERROR=4
      ELSEIF (MAXKLA .GT. KLAMAX .AND. ERROR .EQ. 0) THEN
C
C        ES WURDEN MEHR ALS KLAMAX KLAMMERN GESCHACHTELT
C
         ERROR=5
      ENDIF
C
C     ENDE VON KLAMME
C
      END
C ===== SOURCE: learca.f
C
C
C*DK LEARCA
      FUNCTION LEARCA (X,R,N1,N,NS,TEXT)
C
C   THIS FUNCTION COMPUTES THE INDEX OF THE SMALLER MESHPOINT OF THE
C   INTERVALL CONTAINING THE POINT X IN THE MESH R(NS,I),I=1,N
C   N1 IS THE LEADING DIMENSION OF THE FIELD R, AS SPECIFIED IN THE
C   CALLING PROGRAM AND NS IS A FIXED INDEX.
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: TEXT
      INTEGER, INTENT(IN) :: N1, N, NS
      REAL(DP), INTENT(IN) :: R(N1,*)
      REAL(DP), INTENT(IN) :: X
      INTEGER :: NNN, I, J, LEARCA

      NNN=1
      IF (X.LT.R(NS,1)-1.D-12) GOTO 20
13    DO 10 J=2,N
        I=J
        IF (X-R(NS,J).LE.0.0) GOTO 15
10    CONTINUE
      NNN=N
      IF (X.GT.R(NS,N)+1.D-12) GOTO 20
15    LEARCA=I-1
      RETURN
C
20    WRITE (iunout,*) 'X OUT OF RANGE IN LEARCA'
      WRITE (iunout,*)  X,NNN,R(NS,NNN)
      WRITE (iunout,*) 'LEARCA= ',NNN,' RETURNED TO SUBR. ',TEXT
      LEARCA=NNN
C
      RETURN
      END
C ===== SOURCE: leer.f
C
C
C*DK LEER
      SUBROUTINE LEER (N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      INTEGER :: I
      DO 1 I=1,N
         WRITE (iunout,60)
1        CONTINUE
60    FORMAT ('      ')
      RETURN
      END
C ===== SOURCE: lowercase.f


      subroutine lowercase (zeile)

      IMPLICIT NONE
      character(*), INTENT(INOUT) :: zeile
      INTEGER :: L, I, J, LANF, LEND
      character(26) :: klein, gross
      data klein /'abcdefghijklmnopqrstuvwxyz'/
      data gross /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      LANF=verify(zeile,' ')
      LEND=verify(zeile,' ',.true.)
      if (lanf == 0) return
      l = lend-lanf+1
      zeile(1:l) = zeile(lanf:lend)
      zeile(l+1:lend) = repeat(' ',lanf)
      do i=1,l
        j=index(gross,zeile(i:i))
        if (j>0) zeile(i:i)=klein(j:j)
      end do

      return
      end subroutine lowercase

C ===== SOURCE: margin.f
C
C
      SUBROUTINE MARGIN (NLFLD,NDIM,NMIN,NMAX,LOG)
C
C  FIND MINIMAL AND MAXIMAL LOGICAL OF VALUE LOG IN THE
C  LOGICAL FIELD NLFLD
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      INTEGER, INTENT(OUT) :: NMIN, NMAX
      LOGICAL, INTENT(IN) :: NLFLD(NDIM)
      LOGICAL, INTENT(IN) :: LOG
      INTEGER :: J
C
      IF (LOG) THEN
C
         DO 1 J=1,NDIM
            NMIN=J
            IF (NLFLD(J)) GOTO 2
    1    CONTINUE
         NMIN=NDIM+1
C
    2    NMAX=NDIM
         DO 3 J=NDIM,NMIN,-1
            NMAX=J
            IF (NLFLD(J)) GOTO 4
    3    CONTINUE
C
      ELSE
C
         DO 5 J=1,NDIM
            NMIN=J
            IF (.NOT.NLFLD(J)) GOTO 6
    5    CONTINUE
         NMIN=NDIM+1
C
    6    NMAX=NDIM
         DO 7 J=NDIM,NMIN,-1
            NMAX=J
            IF (.NOT.NLFLD(J)) GOTO 4
    7    CONTINUE
C
      ENDIF
C
    4 RETURN
      END
C ===== SOURCE: masage.f
C
C
C*DK MASAGE
      SUBROUTINE MASAGE (A)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      WRITE (iunout,60) A
 60   FORMAT (1X,A)
      RETURN
      END
C ===== SOURCE: masaj1.f


C*DK MASAJ1
      SUBROUTINE MASAJ1 (A,M,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(13), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N, M(N)
      INTEGER :: J
      WRITE (iunout,60) A
60    FORMAT (1X,A13)
      DO 62 J=1,N
62    WRITE (iunout,61) J,M(J)
61    FORMAT(1X,I2,1X,I6)
      RETURN
      END
C ===== SOURCE: masajr.f


C*DK MASAJR
      SUBROUTINE MASAJR (A,B,N,R)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(24), INTENT(IN) :: A
      CHARACTER(8), INTENT(IN) :: B
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: R
      WRITE (iunout,60) A,B,N,R
60    FORMAT (1X,A24,1X,A8,1X,I2,1X,1PE12.4)
      RETURN
      END
C ===== SOURCE: masal1.f
C
C
C*DK MASAL1
      SUBROUTINE MASAL1 (A,NL,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(6), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: NL(N)
      INTEGER :: I
      WRITE (iunout,60) A,(NL(I),I=1,N)
60    FORMAT (1X,A6/(20L3))
      RETURN
      END
C ===== SOURCE: masbox.f
C
C
C*DK MASBOX
      SUBROUTINE MASBOX (STRING)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: STRING
      CHARACTER(80) :: AST(5)
      INTEGER :: I, J, ILEN

      DO I=1,5
        DO J=1,80
          AST(I)(J:J)=' '
        ENDDO
      ENDDO

      ILEN=LEN(STRING)
      DO I=1,ILEN+6
        AST(1)(I:I)='*'
        AST(5)(I:I)='*'
      ENDDO
      AST(2)(1:1)='*'
      AST(2)(ILEN+6:ILEN+6)='*'
      AST(3)(1:1)='*'
      AST(3)(ILEN+6:ILEN+6)='*'
      AST(4)(1:1)='*'
      AST(4)(ILEN+6:ILEN+6)='*'
C
      DO I=1,ILEN
        AST(3)(I+3:I+3)=STRING(I:I)
      ENDDO
C
      WRITE (iunout,'(//1X)')
      DO I=1,5
        WRITE (iunout,'(A80)') AST(I)
      ENDDO
      WRITE (iunout,'(//1X)')
C
      RETURN
      END
C ===== SOURCE: masbr2.f
C
C
      SUBROUTINE MASBR2 (A,NLFLD,NFIRST0,NFIRST1,N0,N1,M,NBITS)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: NFIRST0, NFIRST1, N0, N1, M, NBITS
      INTEGER, INTENT(IN) :: NLFLD(NFIRST0:NFIRST1,*)
      INTEGER :: I, J, K, JA, IZ, JE, NCH
      LOGICAL BITGET
C
      NCH=110
      IZ=M/NCH
      IF (MOD(M,NCH).NE.0) IZ=IZ+1
C
      WRITE (iunout,'(///1X,A)') A
      DO 1 I=1,IZ
         JA=(I-1)*NCH+1
         JE=MIN(I*NCH,M)
         WRITE (iunout,'(/6X,A1,12(1X,10I1))') 'I',(MOD(J,10),J=JA,JE)
         WRITE (iunout,'(1X,130A1)') ('-',J=JA,JE+16)
         DO 2 K=N0,N1
            WRITE (iunout,'(1X,I4,1X,A1,12(1X,10L1))') K,'I',
     .       (BITGET(NLFLD,NFIRST0,NFIRST1,K,J,NBITS),J=JA,JE)
2        CONTINUE
1     CONTINUE
      RETURN
      END
C ===== SOURCE: masir2.f
C
C
      SUBROUTINE MASIR2 (A,NLFLD,NFIRST0,NFIRST1,N0,N1,M)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: NFIRST0, NFIRST1, N0, N1, M
      INTEGER, INTENT(IN) :: NLFLD(NFIRST0:NFIRST1,*)
      INTEGER :: I, J, K, JA, IZ, JE, NCH
C
      NCH=110
      IZ=M/NCH
      IF (MOD(M,NCH).NE.0) IZ=IZ+1
C
      WRITE (iunout,'(///1X,A)') A
      DO 1 I=1,IZ
         JA=(I-1)*NCH+1
         JE=MIN(I*NCH,M)
         WRITE (iunout,'(/6X,A1,12(1X,10I1))') 'I',(MOD(J,10),J=JA,JE)
         WRITE (iunout,'(1X,130A1)') ('-',J=JA,JE+16)
         DO 2 K=N0,N1
            WRITE (iunout,'(1X,I4,1X,A1,12(1X,10I1))') K,'I',
     .            (NLFLD(K,J),J=JA,JE)
2        CONTINUE
1     CONTINUE
      RETURN
      END
C ===== SOURCE: masj1.f
C
C*DK MASJ1
      SUBROUTINE MASJ1 (A,I)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I
      REAL(DP) :: FI
      IF (I.LT.1000000) THEN
        WRITE (IUNOUT,60) A,I
      ELSEIF (I.LT.100000000) THEN
        WRITE (IUNOUT,70) A,I
      ELSE
        FI=I
        WRITE (IUNOUT,80) A,FI
      ENDIF
60    FORMAT (1X,A8,3X,I6)
70    FORMAT (1X,A8,3X,I8)
80    FORMAT (1X,A8,3X,1PE12.4)
      RETURN
      END
C ===== SOURCE: masj1r.f
C
C*DK MASJ1R
      SUBROUTINE MASJ1R (A,I,B)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(16), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B
      INTEGER, INTENT(IN) :: I
      WRITE (iunout,60) A,I,B
60    FORMAT (1X,A16,3X,I6,3X,1PE12.4)
      RETURN
      END
C ===== SOURCE: masj2.f
C
C*DK MASJ2
      SUBROUTINE MASJ2 (A,I,J)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(16), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I, J
      REAL(DP) :: FJ
      IF (J.LT.1000000) THEN
        WRITE (iunout,60) A,I,J
      ELSEIF (J.LT.100000000) THEN
        WRITE (iunout,70) A,I,J
      ELSE
        fj=j
        WRITE (iunout,80) A,I,fj
      ENDIF
60    FORMAT (1X,A16,2(3X,I6))
70    FORMAT (1X,A16,(3X,I6,3X,I8))
80    FORMAT (1X,A16,(3X,I6,3X,1pe12.4))
      RETURN
      END
C ===== SOURCE: masj2r.f
C
C*DK MASJ2R
      SUBROUTINE MASJ2R (A,I,J,B)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(24), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B
      INTEGER, INTENT(IN) :: I, J
      WRITE (iunout,60) A,I,J,B
60    FORMAT (1X,A24,3X,I6,3X,I6,3X,1PE12.4)
      RETURN
      END
C ===== SOURCE: masj3.f
C
C*DK MASJ3
      SUBROUTINE MASJ3 (A,I,J,K)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(24), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I, J, K
      WRITE (iunout,60) A,I,J,K
60    FORMAT (1X,A24,1X,3(I6,3X))
      RETURN
      END
C ===== SOURCE: masj4.f
C
C*DK MASJ4
      SUBROUTINE MASJ4 (A,I,J,K,L)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(32), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I, J, K, L
      WRITE (iunout,60) A,I,J,K,L
60    FORMAT (1X,A32/1X,4(I6,3X))
      RETURN
      END
C ===== SOURCE: masj5.f
C
C
C*DK MASJ5
      SUBROUTINE MASJ5 (A,I,J,K,L,M)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(40), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I, J, K, L, M
      WRITE (iunout,60) A,I,J,K,L,M
60    FORMAT (1X,A40/1X,5(I6,3X))
      RETURN
      END
C ===== SOURCE: masj6.f
C
C
C*DK MASJ6
      SUBROUTINE MASJ6 (A,I,J,K,L,M,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(48), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I, J, K, L, M, N
      WRITE (iunout,60) A,I,J,K,L,M,N
60    FORMAT (1X,A48/1X,6(I6,3X))
      RETURN
      END
C ===== SOURCE: masjr2.f
C
C*DK MASJR2
      SUBROUTINE MASJR2 (A,I,B,C)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(24), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B, C
      INTEGER, INTENT(IN) :: I
      WRITE (iunout,60) A,I,B,C
60    FORMAT (1X,A24,3X,I6,3X,2(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masl1.f
C
C*DK MASL1
      SUBROUTINE MASL1 (A,B)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B
      WRITE (iunout,60) A,B
60    FORMAT (1X,A,5X,L3)
      RETURN
      END
C ===== SOURCE: masl2.f
C
C*DK MASL2
      SUBROUTINE MASL2 (A,B,C)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B, C
      WRITE (iunout,60) A,B,C
60    FORMAT (1X,A,5X,2L3)
      RETURN
      END
C ===== SOURCE: masl3.f
C
C
C*DK MASL3
      SUBROUTINE MASL3 (A,B,C,D)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B, C, D
      WRITE (iunout,60) A,B,C,D
60    FORMAT (1X,A,5X,3L3)
      RETURN
      END
C ===== SOURCE: masl4.f
C
C
C*DK MASL4
      SUBROUTINE MASL4 (A,B,C,D,E)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B, C, D, E
      WRITE (iunout,60) A,B,C,D,E
60    FORMAT (1X,A,5X,4L3)
      RETURN
      END
C ===== SOURCE: masl5.f
C
C*DK MASL5
      SUBROUTINE MASL5 (A,B,C,D,E,F)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B, C, D, E, F
      WRITE (iunout,60) A,B,C,D,E,F
60    FORMAT (1X,A,5X,5L3)
      RETURN
      END
C ===== SOURCE: masl6.f
C
C*DK MASL6
      SUBROUTINE MASL6 (A,B,C,D,E,F,G)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(36), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B, C, D, E, F, G
      WRITE (iunout,60) A,B,C,D,E,F,G
60    FORMAT (1X,A36,5X,6L3)
      RETURN
      END
C ===== SOURCE: maslr2.f
C
C
      SUBROUTINE MASLR2 (A,NLFLD,NFIRST0,NFIRST1,N0,N1,M)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: NFIRST0, NFIRST1, N0, N1, M
      LOGICAL, INTENT(IN) :: NLFLD(NFIRST0:NFIRST1,*)
      INTEGER :: I, J, K, JA, IZ, JE, NCH
C
      NCH=110
      IZ=M/NCH
      IF (MOD(M,NCH).NE.0) IZ=IZ+1
C
      WRITE (iunout,'(///1X,A)') A
      DO 1 I=1,IZ
         JA=(I-1)*NCH+1
         JE=MIN(I*NCH,M)
         WRITE (iunout,'(/6X,A1,12(1X,10I1))') 'I',(MOD(J,10),J=JA,JE)
         WRITE (iunout,'(1X,130A1)') ('-',J=JA,JE+16)
         DO 2 K=N0,N1
            WRITE (iunout,'(1X,I4,1X,A1,12(1X,10L1))') K,'I',
     .            (NLFLD(K,J),J=JA,JE)
2        CONTINUE
1     CONTINUE
      RETURN
      END
C ===== SOURCE: masprm.f
C
C
C*DK MARPRM
      SUBROUTINE MASPRM (A,NA,MA,B,NB,MB,IERR)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A, B
      INTEGER, INTENT(IN) :: NA, MA, NB, MB
      INTEGER, INTENT(INOUT) :: IERR
      IERR=IERR+1
      WRITE (iunout,*) 'PARAMETER ERROR DETECTED '
      WRITE (iunout,*) A(1:NA),' MUST BE >= ',B(1:NB)
      WRITE (iunout,*) A(1:NA),' = ',MA
      WRITE (iunout,*) B(1:NB),' = ',MB
      RETURN
      END
C ===== SOURCE: masr1.f
C
C
C*DK MASR1
      SUBROUTINE MASR1 (A,B)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B
      WRITE (iunout,60) A,B
60    FORMAT (1X,A8,3X,1PE12.4)
      RETURN
      END
C ===== SOURCE: masr2.f
C
C
C
C*DK MASR2
      SUBROUTINE MASR2 (A,B,C)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(16), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B, C
      WRITE (iunout,60) A,B,C
60    FORMAT (1X,A16,2(3X,1PE12.4))
      RETURN
      END
C ===== SOURCE: masr3.f
C
C*DK MASR3
      SUBROUTINE MASR3 (A,B,C,D)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(24), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B, C, D
      WRITE (iunout,60) A,B,C,D
60    FORMAT (1X,A24,1X,3(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masr4.f
C
C*DK MASR4
      SUBROUTINE MASR4 (A,B,C,D,E)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(32), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B, C, D, E
      WRITE (iunout,60) A,B,C,D,E
60    FORMAT (1X,A32/1X,4(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masr5.f
C
C*DK MASR5
      SUBROUTINE MASR5 (A,B,C,D,E,F)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(40), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B, C, D, E, F
      WRITE (iunout,60) A,B,C,D,E,F
60    FORMAT (1X,A40/1X,5(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masr6.f
C
C*DK MASR6
      SUBROUTINE MASR6 (A,B,C,D,E,F,G)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(48), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B, C, D, E, F, G
      WRITE (iunout,60) A,B,C,D,E,F,G
60    FORMAT (1X,A48/1X,6(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masrr1.f
C
C*DK MASRR1
      SUBROUTINE MASRR1 (A,B,N,IS)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
C  IS IST DIE ANZAHL DER SPALTEN (MAXIMAL: 6)
      CHARACTER(11), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B(*)
      INTEGER, INTENT(IN) :: N, IS
      INTEGER :: IE, I, J, IA
      WRITE (iunout,60) A
60    FORMAT (1X,A11)
      IA=1
      IE=MIN0(N,IA-1+IS)
      DO 62 J=1,N
      IF (IA.GT.IE) GOTO 63
      WRITE (iunout,61) (I,B(I),I=IA,IE)
      IA=IE+1
      IE=MIN0(N,IE+IS)
62    CONTINUE
63    CONTINUE
61    FORMAT(6(1X,I4,1X,1PE12.4))
      RETURN
      END
C ===== SOURCE: masrr2.f
C
C*DK MASRR2
      SUBROUTINE MASRR2 (A,B,C,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(22), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: B(N), C(N)
      INTEGER :: J
      WRITE (iunout,60) A
60    FORMAT (1X,A22)
      DO 5 J=1,N
5        WRITE (iunout,61) J,B(J),C(J)
61    FORMAT (1X,I4,1X,2(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masrr3.f
C
C
C
C*DK MASRR3
      SUBROUTINE MASRR3 (A,B,C,D,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(22), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: B(N), C(N), D(N)
      INTEGER :: J
      WRITE (iunout,60) A
60    FORMAT (1X,A22)
      DO 5 J=1,N
5        WRITE (iunout,61) J,B(J),C(J),D(J)
61    FORMAT (1X,I4,1X,3(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masrr4.f
C
C
C*DK MASRR4
      SUBROUTINE MASRR4 (A,B,C,D,E,N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(22), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: B(N), C(N), D(N), E(N)
      INTEGER :: J
      WRITE (iunout,60) A
60    FORMAT (1X,A22)
      DO 5 J=1,N
5        WRITE (iunout,61) J,B(J),C(J),D(J),E(J)
61    FORMAT (1X,I4,1X,4(1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masxr1.f
C
C
C*DK MASXR1
      SUBROUTINE MASXR1 (A,B,KI,KE,KST,TX)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(14), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B(*)
      INTEGER, INTENT(IN) :: KI, KE, KST
      INTEGER :: K
      IF (KE.LT.KI) RETURN
      WRITE (iunout,60) A
60    FORMAT (1X,A14)
         WRITE (iunout,61) (TX(K),B(K),K=KI,KE,KST)
61    FORMAT (1X,6(A8,1PE12.4,1X))
      RETURN
      END
C ===== SOURCE: masxr2.f
C
C
C*DK MASXR2
      SUBROUTINE MASXR2 (A,B,C,KI,KE,KST,N,M,TX)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(14), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: KI, KE, KST, N, M
      REAL(DP), INTENT(IN) :: B(M),C(N,M)
      INTEGER :: J, K
      IF (KE.LT.KI) RETURN
      WRITE (iunout,60) A
60    FORMAT (1X,A14)
      DO 5 J=1,M
5        WRITE (iunout,61) J,B(J),(TX(K),C(K,J),K=KI,KE,KST)
61    FORMAT (1X,I4,1X,1PE12.4,2X,6(A8,1PE12.4,1X))
      RETURN
      END
C ===== SOURCE: masxr3.f
C
C
C*DK MASXR3
      SUBROUTINE MASXR3 (A,B,C,D,KI,KE,KST,N,M,TX)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(14), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: KI, KE, KST, N, M
      REAL(DP), INTENT(IN) :: B(M),C(N,M), D(N,M)
      INTEGER :: J, K
      IF (KE.LT.KI) RETURN
      WRITE (iunout,60) A
60    FORMAT (1X,A14)
      DO 5 J=1,M
5        WRITE (iunout,61) B(J),(TX(K),C(K,J),D(K,J),K=KI,KE,KST)
61    FORMAT (1X,1PE12.4,2X,6(A8,1X,2(1PE12.4,2X)))
      RETURN
      END
C ===== SOURCE: masxr4.f
C
C
C*DK MASXR4
      SUBROUTINE MASXR4 (A,B,C,D,KI,KE,KST,N,M,TX)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(16), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: KI, KE, KST, N, M
      REAL(DP), INTENT(IN) :: B(M),C(M), D(N,M)
      INTEGER :: J, K
      IF (KE.LT.KI) RETURN
      WRITE (iunout,60) A
60    FORMAT (1X,A16)
      DO 5 J=1,M
5        WRITE (iunout,61) B(J),C(J),(TX(K),D(K,J),K=KI,KE,KST)
61    FORMAT (1X,2(1PE12.4,3X),6(A8,2X,1PE12.4,3X))
      RETURN
      END
C ===== SOURCE: masyr1.f
C
C
C*DK MASYR1
      SUBROUTINE MASYR1 (A,B,NL,I1,N0,N1,M0,M1,TX)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(9), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I1, N0, N1, M0, M1
      REAL(DP), INTENT(IN) :: B(N0:N1,M0:M1)
      LOGICAL, INTENT(IN) :: NL(N0:N1,M0:M1)
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(8) :: KK(100)
      INTEGER :: KKK(100)
      INTEGER :: J, KE, K
      KE=0
      DO 10 J=1,N1
        IF (.NOT.NL(J,I1)) GOTO 10
        KE=KE+1
        KK(KE)=TX(J)
        KKK(KE)=J
10    CONTINUE
      IF (KE.EQ.0) RETURN
      WRITE (iunout,60) A
60    FORMAT (1X,A9)
5        WRITE (iunout,61) (KK(K),B(KKK(K),I1),K=1,KE)
61    FORMAT (1X,6(1A8,1PE12.4,2X))
      RETURN
      END
C ===== SOURCE: maxmn2.f
C
C
C*DK MAXMN2
      SUBROUTINE MAXMN2 (F,ID,IMIN,IMAX,JMIN,JMAX,XMIN,XMAX)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ID, IMIN, IMAX, JMIN, JMAX
      REAL(DP), INTENT(IN) :: F(ID,*)
      REAL(DP), INTENT(OUT) :: XMIN, XMAX

      XMIN = MINVAL(F(IMIN:IMAX,JMIN:JMAX))
      XMAX = MAXVAL(F(IMIN:IMAX,JMIN:JMAX))

      RETURN
      END
C ===== SOURCE: mecker.f



C-----------------------------------------------------------------------
                SUBROUTINE MECKER(ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     AUSGABE DER REGELVERLETZUNGEN
C
C-----------------------------------------------------------------------
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN


      IF     (ERROR .EQ.  1) THEN
        WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT EIN UNGUELTIGES ZEICHEN.'
      ELSEIF (ERROR .EQ.  2) THEN
        WRITE(iunout,*) 
     >  'ZWISCHEN ZWEI KLAMMERAUSDRUECKEN BEFINDET SICH',
     >  ' KEIN OPERATOR.'
      ELSEIF (ERROR .EQ.  3) THEN
         WRITE(iunout,*) 'EIN KLAMMERAUSDRUCK IST LEER.'
      ELSEIF (ERROR .EQ.  4) THEN
         WRITE(iunout,*) 'DER AUSDRUCK IST FALSCH GEKLAMMERT.'
      ELSEIF (ERROR .EQ.  5) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 3 INEINANDER',
     >              'GESCHACHTELTE KLAMMERN.'
      ELSEIF (ERROR .EQ.  6) THEN
        WRITE(iunout,*) 
     >  'EIN OPERAND BESTEHT AUS MEHR ALS ZWEI BUCHSTABEN.'
      ELSEIF (ERROR .EQ.  7) THEN
         WRITE(iunout,*) 'ZWEI OPERATOREN STEHEN NEBENEINANDER.'
      ELSEIF (ERROR .EQ.  8) THEN
         WRITE(iunout,*) 'NACH EINEM OPERATOR FOLGT EINE SCHLIESSENDE',
     >              ' KLAMMER.'
      ELSEIF (ERROR .EQ.  9) THEN
        WRITE(iunout,*) 
     >  'DAS ERSTES ZEICHEN DES AUSDRUCKS IST OPERATOR, ',
     >  'UND KEIN PRAEFIX.'
      ELSEIF (ERROR .EQ. 10) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENDET MIT EINEM OPERATOR.'
      ELSEIF (ERROR .EQ. 11) THEN
         WRITE(iunout,*) 'DAS ERSTE ZEICHEN IN EINEM KLAMMERAUSDRUCK',
     >              ' IST EIN OPERATOR.'
      ELSEIF (ERROR .EQ. 12) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 15 OPERATOREN.'
      ELSEIF (ERROR .EQ. 13) THEN
         WRITE(iunout,*) 'DIE (ANZAHL DER OPERANDEN)-1 IST UNGLEICH ',
     >              ' DER (ANZAHL DER OPERATOREN).'
      ENDIF
C
C     ENDE VON MECKER
C
      END
C ===== SOURCE: ncelln.f
C
C
      SUBROUTINE NCELLN(NCELL,NR,NP,NT,NA,NB,NR1,NP2,NT3,NBL,LR,LP,LT)
C  FIND STANDARD 3D MESH CELL NUMBERS NR,NP,NT FROM 1D ARRAY CELL NUMBER NCELL
C
C  INPUT:
C     LR : NLRAD  LOGICAL: TRUE IF THERE IS A GRID IN X (RADIAL) DIRECTION
C                 NR1ST > 1 THEN
C     LP : NLPOL  SAME, FOR Y (POLOIDAL) GRID
C     LT : NLTOR  SAME, FOR Z (TOROIDAL) GRID
C     NR1: NR1ST  (INPUT BLOCK 2, NR1ST-1 IS THE NUMBER OF X-CELLS)
C     NP2: NP2ND  (INPUT BLOCK 2, NP2ND-1 IS THE NUMBER OF Y-CELLS)
C     NT3: NT3RD  (INPUT BLOCK 2, NT3RD-1 IS THE NUMBER OF Z-CELLS)
C     NBL: NBMLT  DEFAULT: =1  (GT 1 ONLY IF THERE ARE MORE THAN BLOCKS)
C  OUTPUT:
C     NR : X (RADIAL) CELL NUMBER.
C                     IF NR=NR1ST: THIS CELL CONTAINS RADIAL AVERAGE
C     NP : Y (POLOIDAL) CELL NUMBER.
C                     IF NP=NP2ND: THIS CELL CONTAINS POLOIDAL AVERAGE
C     NT : Z (TOROIDAL) CELL NUMBER.
C                     IF NT=NT3RD: THIS CELL CONTAINS TOROIDAL AVERAGE
C
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCELL, NR1, NP2, NT3, NBL
      INTEGER, INTENT(OUT) :: NR, NP, NT, NA, NB
      LOGICAL, INTENT(IN) :: LR,LP,LT
      INTEGER :: NCELL0, NBLCKA, NP2T3, NSTRD, NR1P2, N12, NCELL1,
     .           NCELL2, NCELL3, N23

      NSTRD=NR1*NP2*NT3
      IF (NCELL.GT.NSTRD*NBL) THEN
C  ADDITIONAL CELL OUTSIDE THE STANDARD GRID. RETURN: NA:=ADDITIONAL CELL NO.
        NA=NCELL-NSTRD*NBL
        NB=NBL+1
        NR=0
        NP=1
        NT=1
        RETURN
      ELSEIF (NCELL.EQ.0) THEN
        NA=0
        NB=1
        NR=0
        NP=1
        NT=1
        RETURN
      ENDIF
C
C  AT THIS POINT: 1 <= NCELL <= NSURF, HENCE: NA=0
C
      NR1P2=0
      IF (LP.OR.LT) NR1P2=NR1
      NP2T3=0
      IF (LT) NP2T3=NP2
C
      NCELL0=NCELL
      NA=0
      NB=(NCELL0-1)/NSTRD+1
      NBLCKA=NSTRD*(NB-1)+NA
C
      NCELL1=NCELL0-NBLCKA
      IF (LT) THEN
        N23=(NCELL1-1)/NR1P2+1
        NT=(N23-1)/NP2T3+1
      ELSE
        NT=1
      ENDIF
C
      NCELL2=NCELL1-(NT-1)*NP2T3*NR1P2
      IF (LP) THEN
        N12=NCELL2
        NP=(N12-1)/NR1P2+1
      ELSE
        NP=1
      ENDIF
C
      NCELL3=NCELL2-(NP-1)*NR1P2
      NR=NCELL3
      RETURN
      END
C ===== SOURCE: operat.f



C-----------------------------------------------------------------------
          SUBROUTINE OPERAT(AUSDRU,AKTLEN,OTOREN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFUNG EINER REGELVERLETZUNG BEI OPERATOREN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION
C
         CHARACTER(5), PARAMETER :: FAKTOR='^*/+-'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: OTOREN
C           : ANZAHL DER OPERATOREN IN AUSDRU

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS

      OTOREN=0

      I=1
C
C     REPEAT
C
10    CONTINUE
         POS= INDEX (FAKTOR, AUSDRU(I:I) )
         IF (POS .GT. 0  .AND.  ERROR .EQ. 0 ) THEN
            IF (INDEX(FAKTOR, AUSDRU(I+1:I+1)) .GT. 0) THEN
C
C              2 OPERATOREN STEHEN NEBENEINANDER
C
               ERROR=7
            ELSEIF (INDEX(')',AUSDRU(I+1:I+1)) .GT. 0) THEN
C
C              NACH EINEM OPERATOR FOLGT EINE SCHLIESSENDE KLAMMER
C
               ERROR=8
            ENDIF
            OTOREN=OTOREN+1
         ENDIF
         I=I+1
      IF (I .LT. AKTLEN  .AND.  ERROR .EQ. 0 ) GOTO 10
C
C     UNTIL : BIS AUSDRU VOLLSTAENDIG DURCHLAUFEN
C             ODER FEHLER AUFGETRETEN
C
C     UEBERPRUEFEN, OB DER 1. OPERAND IN AUSDRU EIN PRAEFIX BESITZT
C
      POS=INDEX(FAKTOR,AUSDRU(1:1))
      IF (POS .GT. 0 ) THEN
         IF ( POS .LE. 3) THEN
C
C           AUSDRU BEGINN MIT EINEM OPERATOR
C
            ERROR=9
         ELSE
            OTOREN=OTOREN-1
         ENDIF
      ENDIF

      IF (INDEX(FAKTOR,AUSDRU(AKTLEN:AKTLEN)).GT. 0  .AND.
     >     ERROR .EQ. 0) THEN
C
C        AUSDRU ENDET MIT EINEM OPERATOR
C
         ERROR=10
      ENDIF

      DO 20, I=1,AKTLEN-1
C
C        DIE ENDGUELTIGE ANZAHL DER OPERATOREN ERGIBT SICH AUS
C        DER BISHERIGEN ANZAHL DER OPERATOREN ABZUEGLICH DER
C        ANZAHL DER PRAEFIXE
C
         IF (AUSDRU(I:I) .EQ. '(' .AND.
     >        INDEX(FAKTOR,AUSDRU(I+1:I+1)) .GT. 0 ) THEN
            IF (INDEX(FAKTOR, AUSDRU(I+1:I+1) ) .GT. 3) THEN
               OTOREN=OTOREN-1
            ELSE
C
C              EINER GEOEFFNETEN KLAMMER FOLGT EIN OPERATOR
C
               ERROR=11
            ENDIF
         ENDIF
20    CONTINUE
C
C     ENDE VON OPRATO
C
      END
C ===== SOURCE: oprand.f



C-----------------------------------------------------------------------
           SUBROUTINE OPRAND(AUSDRU,AKTLEN,OANDEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFEN AUF ZULAESSIGE OPERANDEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(26), PARAMETER ::
     P                  BUCHST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: OANDEN
C           : ANZAHL DER OPERANDEN IN AUSDRU

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS


      OANDEN=0
      AUSDRU(AKTLEN+1:AKTLEN+2)='  '

      I=1
C
C     WHILE-1 : SOLANGE AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C               UND KEIN FEHLER AUFGETRETEN
C
20    IF (I .LE. AKTLEN  .AND.  ERROR .EQ. 0 ) THEN
         POS= INDEX (BUCHST, AUSDRU(I:I) )
C
C        WHILE-2 : SOLANGE ANFANG VON EINEM OPERANDEN GEFUNDEN
C                  UND AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C                  UND KEIN FEHLER AUFGETRETEN
C
10       IF (POS .GT. 0  .AND. I .LE. AKTLEN .AND. ERROR .EQ. 0 ) THEN
            IF (INDEX(BUCHST, AUSDRU(I+1:I+1)) .GT. 0  .AND.
     >          INDEX(BUCHST, AUSDRU(I+2:I+2)) .GT. 0 ) THEN
C
C              OPERAND BESTEHT AUS MEHR ALS 2 BUCHSTABEN
C               => ZU LANG
C
               ERROR=6
            ELSE
C
C              GEHE WEITER IM STRING
C
               I=I+2
               POS=INDEX( BUCHST, AUSDRU(I:I) )
               OANDEN=OANDEN+1
            ENDIF
            GOTO 10
         ENDIF
C
C        ENDWHILE-2
C
C        KEIN ANFANG VON EINEM OPERANDEN GEFUNDEN, GEHE WEITER IM STRING
C
         I=I+1
         GOTO 20
      ENDIF
C
C     ENDWHILE-1
C
C     ENDE VON OPRAND
C
      END
C ===== SOURCE: page.f
C
C
C*DK PAGE
      SUBROUTINE PAGE
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      WRITE (iunout,60)
60    FORMAT ('1         ')
      RETURN
      END
C ===== SOURCE: prttal.f
C
C
C
      SUBROUTINE PRTTAL(T1,T2,T3,PROF,NR,NP,NT,NB,NTT,IFLAG,IFILE)
C
C  PRINT VOLUME AVERAGED TALLIES
C  IFLAG=-1:  ONLY HEADER IS PRINTED
C  IFLAG= 0:  ONLY MEAN VALUES IN EACH BLOCK
C  IFLAG= 1:  ADDITIONALLY: 1D PROFILES (THESE MAY BE AVERAGES)
C  IFLAG= 2:  ADDITIONALLY: 2D PROFILES (THESE MAY BE AVERAGES)
C  IFLAG= 3:  ADDITIONALLY: 3D PROFILES
C  IFLAG> 3:  ONLY FULL PROFILES, NO AVERAGES
C
C  IFILE> 0:  WRITE FULL TALLY ONTO STREAM FORT.IFILE
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: T1, T2, T3
      REAL(DP), INTENT(IN) :: PROF(*)
      INTEGER, INTENT(IN) :: NR, NP, NT, NB, NTT, IFLAG, IFILE
      REAL(DP) :: H(6)
      INTEGER :: K(6)
      INTEGER :: JR, JP, JT, IJ, N1DEL, N2DEL, IADD, JA, IA, IB, I,
     .           IC, IT, IP, NRM, NS, NTM, NPM, IRAD, IST, NCOL, IR
      CHARACTER(1) :: TL(72)

      DATA TL/72*'='/
      SAVE

      CALL LEER(3)
      WRITE (iunout,*) TL
      WRITE (iunout,*) TL
      WRITE (iunout,*) 'TALLY:   ',T1
      WRITE (iunout,*) 'SPECIES: ',T2
      WRITE (iunout,*) 'UNITS:   ',T3
      WRITE (iunout,*) TL
      WRITE (iunout,*) TL
      CALL LEER(1)
      IF (IFILE.GT.0) THEN
        OPEN (UNIT=IFILE,position='APPEND')
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) T1
        WRITE (IFILE,*) T2
        WRITE (IFILE,*) T3
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) NR,NP,NT,NB,NTT
        DO IRAD=1,NTT,5
          WRITE (IFILE,*) (PROF(IR),IR=IRAD,MIN(IRAD+4,NTT))
        ENDDO
        CLOSE (UNIT=IFILE)
      ENDIF

11111 IF (IFLAG.LT.0) RETURN
C  NCOL: NUMBER OF PRINTED DATA PER LINE, .LE.6
      NCOL=5
C
      NS=NR*NP*NT*NB
      NRM=MAX(1,NR-1)
      NPM=MAX(1,NP-1)
      NTM=MAX(1,NT-1)
      N1DEL=0
      IF (NP.GT.1.OR.NT.GT.1) N1DEL=NR
      N2DEL=0
      IF (NT.GT.1) N2DEL=NP
C
C LOOP OVER STANDARD MESH BLOCKS
C
      IF (NS.EQ.0) GOTO 50000
      DO 10000 IB=1,NB
      IF (NB.GT.1) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,777) IB
        WRITE (iunout,*) TL
      ENDIF
      IADD=(IB-1)*NR*NP*NT
C
      IF (IFLAG.LE.2) GOTO 1000
      IF (NR.LE.1.OR.NP.LE.1.OR.NT.LE.1) GOTO 1000
C
C  3 D PROFILES
C
      WRITE (iunout,*) TL
      WRITE (iunout,81)
      WRITE (iunout,*) TL
      CALL LEER(1)
      DO 1 JT=1,NTM
        WRITE (iunout,77) JT
        CALL LEER(1)
        DO 11 JP=1,NPM
          IJ=1
          IR=0
          WRITE (iunout,7) JP
110       DO 111 JR=IJ,NRM
            IC=JR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 112
111       CONTINUE
112       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 110
C  NEXT SEGMENT
          CALL LEER(2)
11      CONTINUE
        WRITE (iunout,*) TL
1     CONTINUE
      WRITE (iunout,*) TL
      IF (IFLAG.GT.3) GOTO 10000
C
C  2 D PROFILES
C
1000  CONTINUE
C
      IF (IFLAG.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NP.LE.1) GOTO 2000
      IF (NP.LE.1.AND.NT.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NT.LE.1) GOTO 2000
C
C  POLOIDAL AND TOROIDAL PROFILE, RADIALLY AVERAGED
C
      IF (NTM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1) WRITE (iunout,881)
        IF (NR.GT.1) CALL LEER(1)
        DO 2 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IP=0
220       DO 222 JP=IJ,NPM
            IC=NR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IP=IP+1
            IJ=IJ+1
            K(IP)=JP
            H(IP)=PROF(IC)
            IF (IP.GE.NCOL) GOTO 223
222       CONTINUE
223       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IP)
          IP=0
          IF (IJ.LE.NPM) GOTO 220
C  NEXT SEGMENT
          CALL LEER(2)
2       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND POLOIDAL PROFILE, TOROIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,81)
        IF (NT.GT.1) WRITE (iunout,883)
        IF (NT.GT.1) CALL LEER(1)
        DO 3 JP=1,NPM
          WRITE (iunout,7) JP
          IJ=1
          IR=0
330       DO 333 JR=IJ,NRM
            IC=JR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 334
333       CONTINUE
334       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 330
C  NEXT SEGMENT
          CALL LEER(2)
3       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND TOROIDAL PROFILE, POLOIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NTM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1) WRITE (iunout,882)
        IF (NP.GT.1) CALL LEER(1)
        DO 4 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IR=0
440       DO 444 JR=IJ,NRM
            IC=JR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 445
444       CONTINUE
445       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 440
C  NEXT SEGMENT
          CALL LEER(2)
4       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
C  1 D PROFILES
C
2000  CONTINUE
C
      IF (IFLAG.LE.0) GOTO 3000
C  RADIAL PROFILE, POLOIDALLY AND TOROIDALLY AVERAGED
C
      IF (NRM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1.AND.NT.EQ.1) WRITE (iunout,882)
        IF (NP.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NP.GT.1.AND.NT.GT.1) WRITE (iunout,8883)
        IJ=1
        IR=0
1110    DO 1111 JR=IJ,NRM
          IC=JR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IR=IR+1
          IJ=IJ+1
          K(IR)=JR
          H(IR)=PROF(IC)
          IF (IR.GE.NCOL) GOTO 1112
1111    CONTINUE
1112    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IR)
        IR=0
        IF (IJ.LE.NRM) GOTO 1110
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  POLOIDAL PROFILE, RADIALLY AND TOROIDALLY AVERAGED
C
      IF (NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1.AND.NT.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NR.GT.1.AND.NT.GT.1) WRITE (iunout,8882)
        IJ=1
        IP=0
1220    DO 1222 JP=IJ,NPM
          IC=NR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IP=IP+1
          IJ=IJ+1
          K(IP)=JP
          H(IP)=PROF(IC)
          IF (IP.GE.NCOL) GOTO 1223
1222    CONTINUE
1223    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IP)
        IP=0
        IF (IJ.LE.NPM) GOTO 1220
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  TOROIDAL PROFILE, RADIALLY AND POLOIDALLY AVERAGED
C
      IF (NTM.GT.1) THEN
        IJ=1
        IT=0
        WRITE (iunout,83)
        IF (NR.GT.1.AND.NP.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NP.GT.1) WRITE (iunout,882)
        IF (NR.GT.1.AND.NP.GT.1) WRITE (iunout,8881)
1330    DO 1333 JT=IJ,NTM
          IC=NR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
          IT=IT+1
          IJ=IJ+1
          K(IT)=JT
          H(IT)=PROF(IC)
          IF (IT.GE.NCOL) GOTO 1334
1333    CONTINUE
1334    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IT)
        IT=0
        IF (IJ.LE.NTM) GOTO 1330
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
3000  CONTINUE
      IC=NR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
      WRITE (iunout,8888) PROF(IC)
      WRITE (iunout,*) TL
      CALL LEER(4)
C
10000 CONTINUE
C
50000 CONTINUE
C  ADDITIONAL CELLS
      IF (NTT.GT.NS) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,7777)
        WRITE (iunout,*) TL
      ENDIF
      IJ=NS+1
      IA=0
550   DO 555 JA=IJ,NTT
        IC=JA
        IA=IA+1
        IJ=IJ+1
        K(IA)=JA-NS
        H(IA)=PROF(IC)
        IF (IA.GE.NCOL) GOTO 556
555   CONTINUE
556   CONTINUE
      WRITE (iunout,6) (K(IC),H(IC),IC=1,IA)
      IA=0
      IF (IJ.LE.NTT) GOTO 550
      CALL LEER(2)
C
6     FORMAT (1X,6(I4,2X,1PE12.4,2X))
7     FORMAT (1X,'Y- OR POLOIDAL SEGMENT NUMBER ',I4)
77    FORMAT (1X,'Z- OR TOROIDAL SEGMENT NUMBER ',I4)
777   FORMAT (1X,'STANDARD MESH BLOCK NUMBER ',I4)
7777  FORMAT (1X,'ADDITIONAL CELLS ')
81    FORMAT (1X,'X- OR RADIAL PROFILE ')
82    FORMAT (1X,'Y- OR POLOIDAL PROFILE ')
83    FORMAT (1X,'Z- OR TOROIDAL PROFILE ')
881   FORMAT (1X,'X- OR RADIAL AVERAGE ')
882   FORMAT (1X,'Y- OR POLOIDAL AVERAGE ')
883   FORMAT (1X,'Z- OR TOROIDAL AVERAGE ')
8881  FORMAT (1X,'X- OR RAD. AND Y- OR POL. AVERAGE ',1PE12.4)
8882  FORMAT (1X,'X- OR RAD. AND Z- OR TOR. AVERAGE ',1PE12.4)
8883  FORMAT (1X,'Y- OR POL. AND Z- OR TOR. AVERAGE ',1PE12.4)
8888  FORMAT (1X,'BLOCK AVERAGE ',1PE12.4)
      RETURN
      END
C ===== SOURCE: prttls.f
C
C
      SUBROUTINE PRTTLS(T1,T2,T3,PROF,NR,NP,NT,NB,NTT,IFLAG,IFILE,
     .                  IR1,IR2,IP1,IP2,IT1,IT2)
C
C  SIMILAR TO PRTTAL, BUT FOR TOTAL "SURFACE TALLIES" (FLUXES, AREAS)
C  IFLAG=-1:  ONLY HEADER IS PRINTED
C  IFLAG= 0:  ONLY MEAN VALUES IN EACH BLOCK
C  IFLAG= 1:  ADDITIONALLY: 1D AVERAGES
C  IFLAG= 2:  ADDITIONALLY: 2D AVERAGES
C  IFLAG= 3:  ADDITIONALLY: 3D PROFILES
C  IFLAG> 3:  ONLY FULL PROFILES, NO AVERAGES
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: T1, T2, T3
      REAL(DP), INTENT(IN) :: PROF(*)
      INTEGER, INTENT(IN) :: NR, NP, NT, NB, NTT, IFLAG, IFILE,
     .                       IR1, IR2, IP1, IP2, IT1, IT2
      INTEGER, PARAMETER :: NSTREAM=15
      REAL(DP) :: H(6)
      INTEGER :: K(6), ISTREAM(NSTREAM)
      INTEGER :: JR, JP, JT, IJ, N1DEL, N2DEL, IADD, JA, IA, IB, I,
     .           IC, IT, IP, NRM, NS, NTM, NPM, IRAD, IST, NCOL, IR,
     .           IRM, IPM, ITM
      CHARACTER(1) :: TL(72)

      DATA TL/72*'='/
      DATA ISTREAM/6,50,20,21,29,30,31,32,33,10,11,12,13,14,15/
      SAVE

      CALL LEER(3)
      WRITE (iunout,*) TL
      WRITE (iunout,*) TL
      WRITE (iunout,*) 'TALLY:   ',T1
      WRITE (iunout,*) 'SPECIES: ',T2
      WRITE (iunout,*) 'UNITS:   ',T3
      WRITE (iunout,*) TL
      WRITE (iunout,*) TL
      CALL LEER(1)
      IF (IFILE.GT.0) THEN
        DO IST=1,NSTREAM
          IF (IFILE.EQ.ISTREAM(IST)) GOTO 11111
        ENDDO
        OPEN(UNIT=IFILE,POSITION='APPEND')
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) T1
        WRITE (IFILE,*) T2
        WRITE (IFILE,*) T3
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) NR,NP,NT,NB,NTT
        DO IRAD=1,NTT,5
          WRITE (IFILE,*) (PROF(IR),IR=IRAD,MIN(IRAD+4,NTT))
        ENDDO
        CLOSE (UNIT=IFILE)
      ENDIF
C
11111 IF (IFLAG.LT.0) RETURN
C  NCOL: NUMBER OF PRINTED DATA PER LINE, .LE.6
      NCOL=5
C
      NS=NR*NP*NT*NB
      NRM=MAX(1,NR-1)
      NPM=MAX(1,NP-1)
      NTM=MAX(1,NT-1)
      IRM=MAX(1,IR2-1)
      IPM=MAX(1,IP2-1)
      ITM=MAX(1,IT2-1)
      N1DEL=0
      IF (NP.GT.1.OR.NT.GT.1) N1DEL=NR
      N2DEL=0
      IF (NT.GT.1) N2DEL=NP
C
C LOOP OVER STANDARD MESH BLOCKS
C
      IF (NS.EQ.0) GOTO 50000
      DO 10000 IB=1,NB
      IF (NB.GT.1) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,777) IB
        WRITE (iunout,*) TL
      ENDIF
      IADD=(IB-1)*NR*NP*NT
C
      IF (IFLAG.LE.2) GOTO 1000
      IF (NR.LE.1.OR.NP.LE.1.OR.NT.LE.1) GOTO 1000
C
C  3 D PROFILES: OUT, BECAUSE SURFACE TALLIES
C
C
C  2 D PROFILES
C
1000  CONTINUE
C
      IF (IFLAG.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NP.LE.1) GOTO 2000
      IF (NP.LE.1.AND.NT.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NT.LE.1) GOTO 2000
C
C  POLOIDAL AND TOROIDAL PROFILE, RADIALLY AVERAGED
C
      IF (NTM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1) WRITE (iunout,881)
        IF (NR.GT.1) CALL LEER(1)
        DO 2 JT=IT1,ITM
          WRITE (iunout,77) JT
          IJ=IP1
          IP=0
220       DO 222 JP=IJ,IPM
            IC=NR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IP=IP+1
            IJ=IJ+1
            K(IP)=JP
            H(IP)=PROF(IC)
            IF (IP.GE.NCOL) GOTO 223
222       CONTINUE
223       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IP)
          IP=0
          IF (IJ.LE.IPM) GOTO 220
C  NEXT SEGMENT
          CALL LEER(2)
2       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND POLOIDAL PROFILE, TOROIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,81)
        IF (NT.GT.1) WRITE (iunout,883)
        IF (NT.GT.1) CALL LEER(1)
        DO 3 JP=IP1,IPM
          WRITE (iunout,7) JP
          IJ=IR1
          IR=0
330       DO 333 JR=IJ,IRM
            IC=JR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 334
333       CONTINUE
334       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.IRM) GOTO 330
C  NEXT SEGMENT
          CALL LEER(2)
3       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND TOROIDAL PROFILE, POLOIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NTM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1) WRITE (iunout,882)
        IF (NP.GT.1) CALL LEER(1)
          DO 4 JT=IT1,ITM
          WRITE (iunout,77) JT
          IJ=IR1
          IR=0
440       DO 444 JR=IJ,IRM
            IC=JR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 445
444       CONTINUE
445       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.IRM) GOTO 440
C  NEXT SEGMENT
          CALL LEER(2)
4       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
C  1 D PROFILES
C
2000  CONTINUE
C
      IF (IFLAG.LE.0) GOTO 3000
C  RADIAL PROFILE, POLOIDALLY AND TOROIDALLY AVERAGED
C
      IF (NRM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1.AND.NT.EQ.1) WRITE (iunout,882)
        IF (NP.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NP.GT.1.AND.NT.GT.1) WRITE (iunout,8883)
        IJ=IR1
        IR=0
1110    DO 1111 JR=IJ,IRM
          IC=JR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IR=IR+1
          IJ=IJ+1
          K(IR)=JR
          H(IR)=PROF(IC)
          IF (IR.GE.NCOL) GOTO 1112
1111    CONTINUE
1112    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IR)
        IR=0
        IF (IJ.LE.IRM) GOTO 1110
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  POLOIDAL PROFILE, RADIALLY AND TOROIDALLY AVERAGED
C
      IF (NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1.AND.NT.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NR.GT.1.AND.NT.GT.1) WRITE (iunout,8882)
        IJ=IP1
        IP=0
1220    DO 1222 JP=IJ,IPM
          IC=NR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IP=IP+1
          IJ=IJ+1
          K(IP)=JP
          H(IP)=PROF(IC)
          IF (IP.GE.NCOL) GOTO 1223
1222    CONTINUE
1223    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IP)
        IP=0
        IF (IJ.LE.IPM) GOTO 1220
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  TOROIDAL PROFILE, RADIALLY AND POLOIDALLY AVERAGED
C
      IF (NTM.GT.1) THEN
        IJ=IT1
        IT=0
        WRITE (iunout,83)
        IF (NR.GT.1.AND.NP.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NP.GT.1) WRITE (iunout,882)
        IF (NR.GT.1.AND.NP.GT.1) WRITE (iunout,8881)
1330    DO 1333 JT=IJ,ITM
          IC=NR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
          IT=IT+1
          IJ=IJ+1
          K(IT)=JT
          H(IT)=PROF(IC)
          IF (IT.GE.NCOL) GOTO 1334
1333    CONTINUE
1334    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IT)
        IT=0
        IF (IJ.LE.ITM) GOTO 1330
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
3000  CONTINUE
      IC=NR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
      WRITE (iunout,8888) PROF(IC)
      WRITE (iunout,*) TL
      CALL LEER(4)
C
10000 CONTINUE
C
50000 CONTINUE
      CALL LEER(2)
C
6     FORMAT (1X,6(I4,2X,1PE12.4,2X))
7     FORMAT (1X,'Y- OR POLOIDAL SEGMENT NUMBER ',I4)
77    FORMAT (1X,'Z- OR TOROIDAL SEGMENT NUMBER ',I4)
777   FORMAT (1X,'STANDARD MESH BLOCK NUMBER ',I4)
7777  FORMAT (1X,'ADDITIONAL CELLS ')
81    FORMAT (1X,'X- OR RADIAL PROFILE ')
82    FORMAT (1X,'Y- OR POLOIDAL PROFILE ')
83    FORMAT (1X,'Z- OR TOROIDAL PROFILE ')
881   FORMAT (1X,'X- OR RADIAL TOTAL ')
882   FORMAT (1X,'Y- OR POLOIDAL TOTAL ')
883   FORMAT (1X,'Z- OR TOROIDAL TOTAL ')
8881  FORMAT (1X,'X- OR RAD. AND Y- OR POL. TOTAL ',1PE12.4)
8882  FORMAT (1X,'X- OR RAD. AND Z- OR TOR. TOTAL ',1PE12.4)
8883  FORMAT (1X,'Y- OR POL. AND Z- OR TOR. TOTAL ',1PE12.4)
8888  FORMAT (1X,'BLOCK TOTAL ',1PE12.4)
      RETURN
      END
C ===== SOURCE: prtvol.f


      SUBROUTINE PRTVOL(T1,T2,T3,PROF,NR,NP,NT,NB,NTT,IFLAG,IFILE)
C
C  SIMILAR TO PRTTAL, BUT FOR "TOTAL TALLIES" SUCH AS CELL VOLUMES
C  IFLAG=-1:  ONLY HEADER IS PRINTED
C  IFLAG= 0:  ONLY MEAN VALUES IN EACH BLOCK
C  IFLAG= 1:  ADDITIONALLY: 1D AVERAGES
C  IFLAG= 2:  ADDITIONALLY: 2D AVERAGES
C  IFLAG= 3:  ADDITIONALLY: 3D PROFILES
C  IFLAG> 3:  ONLY FULL PROFILES, NO AVERAGES
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: T1, T2, T3
      REAL(DP), INTENT(IN) :: PROF(*)
      INTEGER, INTENT(IN) :: NR, NP, NT, NB, NTT, IFLAG, IFILE
      INTEGER, PARAMETER :: NSTREAM=15
      REAL(DP) :: H(6)
      INTEGER :: K(6), ISTREAM(NSTREAM)
      INTEGER :: JR, JP, JT, IJ, N1DEL, N2DEL, IADD, JA, IA, IB, I,
     .           IC, IT, IP, NRM, NS, NTM, NPM, IRAD, IST, NCOL, IR
      CHARACTER(1) :: TL(72)

      DATA TL/72*'='/
      DATA ISTREAM/6,50,20,21,29,30,31,32,33,10,11,12,13,14,15/
      SAVE
      CALL LEER(3)
      WRITE (iunout,*) TL
      WRITE (iunout,*) TL
      WRITE (iunout,*) 'TALLY:   ',T1
      WRITE (iunout,*) 'SPECIES: ',T2
      WRITE (iunout,*) 'UNITS:   ',T3
      WRITE (iunout,*) TL
      WRITE (iunout,*) TL
      CALL LEER(1)
      IF (IFILE.GT.0) THEN
        DO IST=1,NSTREAM
          IF (IFILE.EQ.ISTREAM(IST)) GOTO 11111
        ENDDO
        OPEN (UNIT=IFILE,POSITION='APPEND')
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) T1
        WRITE (IFILE,*) T2
        WRITE (IFILE,*) T3
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) NR,NP,NT,NB,NTT
        DO IRAD=1,NTT,5
          WRITE (IFILE,*) (PROF(IR),IR=IRAD,MIN(IRAD+4,NTT))
        ENDDO
        CLOSE (UNIT=IFILE)
      ENDIF
C
11111 IF (IFLAG.LT.0) RETURN
C  NCOL: NUMBER OF PRINTED DATA PER LINE, .LE.6
      NCOL=5
C
      NS=NR*NP*NT*NB
      NRM=MAX(1,NR-1)
      NPM=MAX(1,NP-1)
      NTM=MAX(1,NT-1)
      N1DEL=0
      IF (NP.GT.1.OR.NT.GT.1) N1DEL=NR
      N2DEL=0
      IF (NT.GT.1) N2DEL=NP
C
C LOOP OVER STANDARD MESH BLOCKS
C
      IF (NS.EQ.0) GOTO 50000
      DO 10000 IB=1,NB
      IF (NB.GT.1) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,777) IB
        WRITE (iunout,*) TL
      ENDIF
      IADD=(IB-1)*NR*NP*NT
C
      IF (IFLAG.LE.2) GOTO 1000
      IF (NR.LE.1.OR.NP.LE.1.OR.NT.LE.1) GOTO 1000
C
C  3 D PROFILES
C
      WRITE (iunout,*) TL
      WRITE (iunout,81)
      WRITE (iunout,*) TL
      CALL LEER(1)
      DO 1 JT=1,NTM
        WRITE (iunout,77) JT
        CALL LEER(1)
        DO 11 JP=1,NPM
          IJ=1
          IR=0
          WRITE (iunout,7) JP
110       DO 111 JR=IJ,NRM
            IC=JR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 112
111       CONTINUE
112       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 110
C  NEXT SEGMENT
          CALL LEER(2)
11      CONTINUE
        WRITE (iunout,*) TL
1     CONTINUE
      WRITE (iunout,*) TL
      IF (IFLAG.GT.3) GOTO 10000
C
C  2 D PROFILES
C
1000  CONTINUE
C
      IF (IFLAG.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NP.LE.1) GOTO 2000
      IF (NP.LE.1.AND.NT.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NT.LE.1) GOTO 2000
C
C  POLOIDAL AND TOROIDAL PROFILE, RADIALLY AVERAGED
C
      IF (NTM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1) WRITE (iunout,881)
        IF (NR.GT.1) CALL LEER(1)
        DO 2 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IP=0
220       DO 222 JP=IJ,NPM
            IC=NR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IP=IP+1
            IJ=IJ+1
            K(IP)=JP
            H(IP)=PROF(IC)
            IF (IP.GE.NCOL) GOTO 223
222       CONTINUE
223       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IP)
          IP=0
          IF (IJ.LE.NPM) GOTO 220
C  NEXT SEGMENT
          CALL LEER(2)
2       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND POLOIDAL PROFILE, TOROIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,81)
        IF (NT.GT.1) WRITE (iunout,883)
        IF (NT.GT.1) CALL LEER(1)
        DO 3 JP=1,NPM
          WRITE (iunout,7) JP
          IJ=1
          IR=0
330       DO 333 JR=IJ,NRM
            IC=JR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 334
333       CONTINUE
334       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 330
C  NEXT SEGMENT
          CALL LEER(2)
3       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND TOROIDAL PROFILE, POLOIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NTM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1) WRITE (iunout,882)
        IF (NP.GT.1) CALL LEER(1)
        DO 4 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IR=0
440       DO 444 JR=IJ,NRM
            IC=JR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 445
444       CONTINUE
445       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 440
C  NEXT SEGMENT
          CALL LEER(2)
4       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
C  1 D PROFILES
C
2000  CONTINUE
C
      IF (IFLAG.LE.0) GOTO 3000
C  RADIAL PROFILE, POLOIDALLY AND TOROIDALLY AVERAGED
C
      IF (NRM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1.AND.NT.EQ.1) WRITE (iunout,882)
        IF (NP.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NP.GT.1.AND.NT.GT.1) WRITE (iunout,8883)
        IJ=1
        IR=0
1110    DO 1111 JR=IJ,NRM
          IC=JR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IR=IR+1
          IJ=IJ+1
          K(IR)=JR
          H(IR)=PROF(IC)
          IF (IR.GE.NCOL) GOTO 1112
1111    CONTINUE
1112    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IR)
        IR=0
        IF (IJ.LE.NRM) GOTO 1110
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  POLOIDAL PROFILE, RADIALLY AND TOROIDALLY AVERAGED
C
      IF (NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1.AND.NT.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NR.GT.1.AND.NT.GT.1) WRITE (iunout,8882)
        IJ=1
        IP=0
1220    DO 1222 JP=IJ,NPM
          IC=NR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IP=IP+1
          IJ=IJ+1
          K(IP)=JP
          H(IP)=PROF(IC)
          IF (IP.GE.NCOL) GOTO 1223
1222    CONTINUE
1223    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IP)
        IP=0
        IF (IJ.LE.NPM) GOTO 1220
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  TOROIDAL PROFILE, RADIALLY AND POLOIDALLY AVERAGED
C
      IF (NTM.GT.1) THEN
        IJ=1
        IT=0
        WRITE (iunout,83)
        IF (NR.GT.1.AND.NP.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NP.GT.1) WRITE (iunout,882)
        IF (NR.GT.1.AND.NP.GT.1) WRITE (iunout,8881)
1330    DO 1333 JT=IJ,NTM
          IC=NR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
          IT=IT+1
          IJ=IJ+1
          K(IT)=JT
          H(IT)=PROF(IC)
          IF (IT.GE.NCOL) GOTO 1334
1333    CONTINUE
1334    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IT)
        IT=0
        IF (IJ.LE.NTM) GOTO 1330
        CALL LEER(1)
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
3000  CONTINUE
      IC=NR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
      WRITE (iunout,8888) PROF(IC)
      WRITE (iunout,*) TL
      CALL LEER(4)
C
10000 CONTINUE
C
50000 CONTINUE
C  ADDITIONAL CELLS
      IF (NTT.GT.NS) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,7777)
        WRITE (iunout,*) TL
      ENDIF
      IJ=NS+1
      IA=0
550   DO 555 JA=IJ,NTT
        IC=JA
        IA=IA+1
        IJ=IJ+1
        K(IA)=JA-NS
        H(IA)=PROF(IC)
        IF (IA.GE.NCOL) GOTO 556
555   CONTINUE
556   CONTINUE
      WRITE (iunout,6) (K(IC),H(IC),IC=1,IA)
      IA=0
      IF (IJ.LE.NTT) GOTO 550
      CALL LEER(2)
C
6     FORMAT (1X,6(I4,2X,1PE12.4,2X))
7     FORMAT (1X,'Y- OR POLOIDAL SEGMENT NUMBER ',I4)
77    FORMAT (1X,'Z- OR TOROIDAL SEGMENT NUMBER ',I4)
777   FORMAT (1X,'STANDARD MESH BLOCK NUMBER ',I4)
7777  FORMAT (1X,'ADDITIONAL CELLS ')
81    FORMAT (1X,'X- OR RADIAL PROFILE ')
82    FORMAT (1X,'Y- OR POLOIDAL PROFILE ')
83    FORMAT (1X,'Z- OR TOROIDAL PROFILE ')
881   FORMAT (1X,'X- OR RADIAL TOTAL ')
882   FORMAT (1X,'Y- OR POLOIDAL TOTAL ')
883   FORMAT (1X,'Z- OR TOROIDAL TOTAL ')
8881  FORMAT (1X,'X- OR RAD. AND Y- OR POL. TOTAL ',1PE12.4)
8882  FORMAT (1X,'X- OR RAD. AND Z- OR TOR. TOTAL ',1PE12.4)
8883  FORMAT (1X,'Y- OR POL. AND Z- OR TOR. TOTAL ',1PE12.4)
8888  FORMAT (1X,'BLOCK TOTAL ',1PE12.4)
      RETURN
      END
C ===== SOURCE: rdcn.f
C
C
      SUBROUTINE RDCN (ERSETZ,CONST)
      USE PRECISION
      IMPLICIT NONE
C
      CHARACTER(*), INTENT(IN) :: ERSETZ
      REAL(DP), INTENT(OUT) :: CONST
      INTEGER :: IEXPO, IW, IPUNKT, ILEN, ICON
      CHARACTER(10) :: FORM
C
      ILEN=INDEX(ERSETZ,'>')-2
      FORM='          '
C
C   INTEGER?
C
      IF (INDEX(ERSETZ,'.').EQ.0) THEN
        FORM(1:5)='(I  )'
        IF (ILEN.GT.9) WRITE (FORM(3:4),'(I2)') ILEN
        IF (ILEN.LE.9) WRITE (FORM(3:3),'(I1)') ILEN
        READ (ERSETZ(2:ILEN+1),FORM) ICON
        CONST=FLOAT(ICON)
C
      ELSE
C
C   REAL CONSTANT
C
        IPUNKT=INDEX(ERSETZ,'.')-1
        IEXPO=INDEX(ERSETZ,'E')
        IF (IEXPO.EQ.0) IEXPO=INDEX(ERSETZ,'D')
C
        IF (IEXPO.EQ.0) THEN
C   F-FORMAT
          FORM(1:2)='(F'
          IF (ILEN.GT.9) THEN
            WRITE (FORM(3:4),'(I2)') ILEN
            IW=4
          ELSE
            WRITE (FORM(3:3),'(I1)') ILEN
            IW=3
          ENDIF
          FORM(IW+1:IW+3)='. )'
          WRITE (FORM(IW+2:IW+2),'(I1)') ILEN-IPUNKT
          READ (ERSETZ(2:ILEN+1),FORM) CONST
C
        ELSE
C
C    E-FORMAT
          IEXPO=IEXPO-1
          FORM(1:2)='(E'
          IF (ILEN.GT.9) THEN
            WRITE (FORM(3:4),'(I2)') ILEN
            IW=4
          ELSE
            WRITE (FORM(3:3),'(I1)') ILEN
            IW=3
          ENDIF
          FORM(IW+1:IW+3)='. )'
          WRITE (FORM(IW+2:IW+2),'(I1)') IEXPO-(IPUNKT+1)
          READ (ERSETZ(2:ILEN+1),FORM) CONST
        ENDIF
      ENDIF
C
      RETURN
      END
C ===== SOURCE: readtl.f
C
C
      SUBROUTINE READTL(T1,T2,T3,PROF,NR,NP,NT,NB,NTT,IFLAG,IFILE)
C
C  READ VOLUME AVERAGED TALLIES FROM STREAM IFILE
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(72), INTENT(OUT) :: T1
!pb      CHARACTER(72), INTENT(OUT) :: T2, T3
      CHARACTER(24), INTENT(OUT) :: T2, T3
      REAL(DP), INTENT(OUT) :: PROF(*)
      INTEGER, INTENT(IN) :: IFILE, IFLAG
      INTEGER, INTENT(OUT) :: NR, NP, NT, NB, NTT
      INTEGER, PARAMETER :: NSTREAM=15
      REAL(DP) :: H(6)
      INTEGER :: K(6), ISTREAM(NSTREAM)
      INTEGER :: I, IRAD, IST, IR
!PB      CHARACTER(72) :: TL(72)
      CHARACTER(72) :: TL

C     DATA TL/72*'='/
      DATA ISTREAM/6,50,20,21,29,30,31,32,33,10,11,12,13,14,15/
      SAVE
      IF (IFILE.GT.0) THEN
        DO IST=1,NSTREAM
          IF (IFILE.EQ.ISTREAM(IST)) THEN
            WRITE (iunout,*) 'ERROR IN INPUT BLOCK 5, INDPRO '
            WRITE (iunout,*) 'STREAM NO. ',IFILE,' IS NOT AVAILABLE'
            GOTO 11111
          ENDIF
        ENDDO
        READ (IFILE,'(A72)') TL
        READ (IFILE,'(A72)') TL
        READ (IFILE,'(A72)') T1
        READ (IFILE,'(A24)') T2
        READ (IFILE,'(A24)') T3
        READ (IFILE,'(A72)') TL
        READ (IFILE,'(A72)') TL
        READ (IFILE,*) NR,NP,NT,NB,NTT
        DO IRAD=1,NTT,5
          READ (IFILE,*) (PROF(IR),IR=IRAD,IRAD+4)
        ENDDO
      ENDIF
11111 CONTINUE
      RETURN
      END
C ===== SOURCE: resetp.f
C
C
C*DK RESETP
      SUBROUTINE RESETP (CRV,DEN,JC,ISEP,M1,N)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JC, ISEP, M1,N
      REAL(DP), INTENT(OUT) :: CRV(M1,N)
      REAL(DP), INTENT(IN) :: DEN(N)
      INTEGER :: J
      DO 5 J=ISEP,N
5        CRV(JC,J-ISEP+1)=DEN(J)
      RETURN
      END
C ===== SOURCE: rotate.f

C
C*DK ROTATE
      SUBROUTINE ROTATE(VLABX,VLABY,VLABZ,VLOCX,VLOCY,VLOCZ,
     .                  CX,CY,CZ,CS)
C
C   VLAB IST DIE RICHTUNG DES EINFLIEGENDEN TEILCHENS IM LABORSYSTEM
C   VLOC IST DIE RICHTUNG DES REFLEKTIERTEN TEILCHENS IM LOKALEN SYSTEM
C   C IST DIE POSITIVE X RICHTUNG IM LOKALEN SYSTEM
C   ES WIRD VLOC INS LABORSYSTEM ZURUECKTRANSFORMIERT UND ALS VLAB
C   ZURUECKGEGEBEN
C   DAS LOKALE SYSTEM WIRD SO BESTIMMT, DASS DAS TEILCHEN IN SEINER
C   X-Z-EBENE EINFLIEGT, MIT POSITIVER Z-GESCHWINDIGKEIT.
C   ES WIRD VORAUSGESETZT, DASS CS = COS(C,VLAB) POSITIV
C   UND DASS VLOCX NEGATIV (D.H. BEIM OUTPUT COS(C,VLAB).LE.0.)
C
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: VLABX, VLABY, VLABZ, VLOCX, VLOCY,
     .                           VLOCZ
      REAL(DP), INTENT(IN) :: CX, CY, CZ, CS
      REAL(DP) :: SS, SSI, A2, A3, B2, B3, C2, C3

      SS=SQRT(1.-CS*CS)
      SSI=1./SS
C
C     A1=CX
C     B1=CY
C     C1=CZ
C
      A3=(VLABX-CS*CX)*SSI
      B3=(VLABY-CS*CY)*SSI
      C3=(VLABZ-CS*CZ)*SSI
C   (A2,B2,C2)=(A3,B3,C3) KREUZ (A1,B1,C1)
C
      A2=B3*CZ-C3*CY
      B2=C3*CX-A3*CZ
      C2=A3*CY-B3*CX
C
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLABX=CX*VLOCX+A2*VLOCY+A3*VLOCZ
      VLABY=CY*VLOCX+B2*VLOCY+B3*VLOCZ
      VLABZ=CZ*VLOCX+C2*VLOCY+C3*VLOCZ
      RETURN
C
      ENTRY ROTATF(VLABX,VLABY,VLABZ,VLOCX,VLOCY,VLOCZ,
     .             CX,CY,CZ)
C  HIER IST ENTWEDER CS=1., D.H. SENKRECHTER EINFLUG, ODER
C  DIE ORIENTIERUNG DER Y-Z-ACHSEN IM LOCALEN SYSTEM SPIELT
C  WEGEN DER SYMMETRIE DER VERTEILUNG DES REFLEXIONSWINKELS
C  KEINE ROLLE, IST INSB. UNABHAENGIG VON VLAB WAEHLBAR
C
C  1. FALL:  ABS(CZ).NE.1.
C
      IF (ABS(CZ).GE.0.99999) GOTO 1
C
      SS=SQRT(CY*CY+CX*CX)
      SSI=1./SS
C
      A2=-CY*SSI
      B2=CX*SSI
C     C2=0.
C
      A3=-CZ*B2
      B3=CZ*A2
C     C3=SS
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLABX=CX*VLOCX+A2*VLOCY+A3*VLOCZ
      VLABY=CY*VLOCX+B2*VLOCY+B3*VLOCZ
      VLABZ=CZ*VLOCX+         SS*VLOCZ
      RETURN
C
C  2. FALL: CZ=1. ODER CZ=-1., D.H. CX=CY=0.
C
 1    CONTINUE
C     A2=0.
C     B2=-CZ
C     C2=0.
C
C     A3=1.=CZ*CZ
C     B3=0.
C     C3=0.
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLABX=                     VLOCZ
      VLABY=        -CZ*VLOCY
      VLABZ=CZ*VLOCX
      RETURN
C
      ENTRY ROTATI(VLABX,VLABY,VLABZ,VLOCX,VLOCY,VLOCZ,
     .             CX,CY,CZ)
C  WIE BEI ENTRY ROTATF, ABER ES WIRD MIT INVERSER
C  (=TRANSPONIERTER) MATRIX GEDREHT.
C
C  1. FALL:  ABS(CZ).NE.1.
C
      IF (ABS(CZ).GE.0.99999) GOTO 2
C
      SS=SQRT(CY*CY+CX*CX)
      SSI=1./SS
C
      A2=-CY*SSI
      B2=CX*SSI
C     C2=0.
C
      A3=-CZ*B2
      B3=CZ*A2
C     C3=SS
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLOCX=CX*VLABX+CY*VLABY+CZ*VLABZ
      VLOCY=A2*VLABX+B2*VLABY
      VLOCZ=A3*VLABX+B3*VLABY+SS*VLABZ
      RETURN
C
C  2. FALL: CZ=1. ODER CZ=-1., D.H. CX=CY=0.
C
 2    CONTINUE
C     A2=0.
C     B2=-CZ
C     C2=0.
C
C     A3=1.=CZ*CZ
C     B3=0.
C     C3=0.
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLOCX=                  CZ*VLABZ
      VLOCY=        -CZ*VLABY
      VLOCZ=   VLABX
      RETURN
      END
C ===== SOURCE: ruksub.f



C-----------------------------------------------------------------------
        SUBROUTINE RUKSUB(TEIL, IPART, PART, IARITH, ARITH, HILFE)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ERSETZEN VON '^' DURCH '**' UND EINFUEGEN VON 'Z' FUER
C     ZERLEGUNG
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(10), PARAMETER :: ZIFFER='0123456789'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

C
C     EIN/AUSGABEPARAMETER :
C
         CHARACTER(*), INTENT(INOUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER, INTENT(INOUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(*), INTENT(INOUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER, INTENT(INOUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

         CHARACTER(*), INTENT(INOUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, J

      DO 10, J=1,TEIL
C
C        ERSETZEN VON '^' DURCH '**' IN PART(J)
C
         I=2
20       IF (I .LT. IPART(J)) THEN
            IF (PART(J)(I:I) .EQ. '^') THEN
CPB            HILFE=PART(J)(1:I-1)//'**'//PART(J)(I+1:IPART(J))//' '
CPB            I=I+1
CPB            IPART(J)=IPART(J)+1
CPB            PART(J)=HILFE
            ENDIF
            I=I+1
            GOTO 20
         ENDIF
C
C        ERSETZEN VON '^' DURCH '**' IN ARITH(J)
C
         I=2
30       IF (I .LT. IARITH(J)) THEN
            IF (ARITH(J)(I:I) .EQ. '^') THEN
CPB            HILFE=ARITH(J)(1:I-1)//'**'//ARITH(J)(I+1:IARITH(J))//' '
CPB            I=I+1
CPB            IARITH(J)=IARITH(J)+1
CPB            ARITH(J)=HILFE
            ENDIF
            I=I+1
            GOTO 30
         ENDIF
C
C        EINFUEGEN  VON 'Z' FUER ZERLEGUNG IN PART(J)
C
         I=1
40       IF (I .LT. IPART(J)) THEN
            IF (INDEX(ZIFFER,PART(J)(I:I)) .GT. 0) THEN
               IF (I .EQ. 1) THEN
                  HILFE='Z'//PART(J)(I:IPART(J))//' '
               ELSE
                  HILFE=PART(J)(1:I-1)//'Z'//
     >                   PART(J)(I:IPART(J))//' '
               ENDIF

               I=I+2
               IPART(J)=IPART(J)+1
               PART(J)=HILFE
            ENDIF
            I=I+1
            GOTO 40
         ENDIF
C
C        EINFUEGEN VON 'Z' FUER ZERLEGUNG IN ARITH(J)
C
         I=1
50       IF (I .LT. IARITH(J)) THEN
            IF (INDEX(ZIFFER,ARITH(J)(I:I)) .GT. 0) THEN
               IF (I .EQ. 1) THEN
                  HILFE='Z'//ARITH(J)(I:IARITH(J))//' '
               ELSE
                  HILFE=ARITH(J)(1:I-1)//'Z'//
     >                   ARITH(J)(I:IARITH(J))//' '
               ENDIF

               I=I+2
               IARITH(J)=IARITH(J)+1
               ARITH(J)=HILFE
            ENDIF
            I=I+1
            GOTO 50
         ENDIF

10       CONTINUE
C
C     ENDE VON RUKSUB
C
      END
C ===== SOURCE: schrit.f



C-----------------------------------------------------------------------
      SUBROUTINE SCHRIT(AUSDRU, AKTLEN, ALPHA, OMEGA, TEIL,
     >                  IPART, PART, IARITH, ARITH)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ZERLEGUNG EINES TEILAUSDRUCKS,
C     DIESER AUSDRUCK ENTHAELT KEINE KLAMMERN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(5), PARAMETER :: FAKTOR='^*/+-'

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         INTEGER, INTENT(INOUT) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

         INTEGER, INTENT(INOUT) :: ALPHA
C           : POSITION DES ERSTEN ZEICHENS VOM TEILAUSDRUCK

         INTEGER, INTENT(INOUT) :: OMEGA
C           : POSITION DES LETZTEN ZEICHENS VOM TEILAUSDRUCK

         CHARACTER(*), INTENT(INOUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER, INTENT(INOUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(*), INTENT(INOUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER, INTENT(INOUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

C
C     LOKALE VARIABLEN :
C
         INTEGER :: BEGINN
C           : POSITION IN AUSDRU, BEI DER DIE ZERLEGUNG BEGINNT

         INTEGER :: ENDE
C           : POSITION IN AUSDRU, BEI DER DIE ZERLEGUNG BEENDET WIRD

         INTEGER :: OPER1
C           : POSITION DES 1. OPERATORS IN AUSDRU

         INTEGER :: OPER2
C           : POSITION DES 2. OPERATORS IN AUSDRU

         INTEGER :: RANGE1
C           : RANG DES 1. OPERATORS

         INTEGER :: RANGE2
C           : RANG DES 2. OPERATORS

         LOGICAL :: ERFOLG
C           : GIBT AN, OB EINE ZERLEGUNG ERFOLGEN DARF

         LOGICAL :: PREFIX
C           : GIBT AN, OB 1.OPERAND IN AUSDRUCK EIN PRAEFIX BESITZT

         CHARACTER(2) :: TEILCH
C           : INHALT VON TEIL

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS


      BEGINN=ALPHA
      ENDE=OMEGA
      ERFOLG=.FALSE.
      I=BEGINN-1
C
C     SUCHE NACH DEM 1. OPERATOR IM TEILSTRING
C
10       CONTINUE
         I=I+1
         POS=INDEX(FAKTOR, ausdru(I:I))
      IF (POS .LE. 0 .AND. I .LT. omega) GOTO 10

      IF (I .EQ. BEGINN) THEN
         PREFIX=.TRUE.
      ELSE
         PREFIX=.FALSE.
      ENDIF
C
C     WHILE : SOLANGE TEILSTRING EINEN OPERATOR ENTHAELT
C
60    IF (POS .GT. 0) THEN
C
C        POSITION UND RANGBESTIMMUNG DES OPERATORS
C
         OPER1=I
         RANGE1=POS/2

C
C        REPEAT-1
C
30          CONTINUE

C
C           SUCHE NACH DEM 2. OPERATOR IM TEILSTRING
C           REPEAT-2
C
20             CONTINUE
               I=I+1
               POS=INDEX(FAKTOR,ausdru(I:I))
            IF (POS .LE. 0  .AND.  I .LT. omega) GOTO 20
C
C           UNTIL-2 : BIS TEILSTRING VOLLSTAENDIG DURCHLAUFEN ODER
C                     WEITEREN OPERATOR IM TEILSTRING GEFUNDEN
C
            IF (POS .GT. 0) THEN
C
C              WEITEREN OPERATOR IN TEILSTRING GEFUNDEN,
C              POSITION UND RANGBESTIMMUNG DES OPERATORS
C
               OPER2=I
               RANGE2=POS/2
C
C              UEBERPRUEFEN, OB NOCH EIN OPERATOR IM TEILSTRING
C              GESUCHT WERDEN MUSS, ODER EINE ZERLEGUNG ERFOLGEN KANN
C
               IF (RANGE2 .LT. RANGE1 .AND. .NOT. PREFIX) THEN
                  BEGINN=OPER1+1
                  OPER1=OPER2
                  RANGE1=RANGE2
               ELSEIF (RANGE2 .EQ. RANGE1  .AND.  RANGE1 .EQ. 0
     >                  .AND. .NOT. PREFIX) THEN
                  BEGINN=OPER1+1
                  OPER1=OPER2
                  RANGE1=RANGE2
               ELSE
                  ENDE=OPER2-1
                  ERFOLG=.TRUE.
C                 zerlegung von beginn bis ende
               ENDIF

            ELSE
C
C              ES IST KEIN WEITERER OPERATOR IN TEILSTRING VORHANDEN,
C              DIE POSITION DES LETZTEN OPERATORS IST OPER1
C
               ERFOLG=.TRUE.
            ENDIF
         IF ( .NOT. ERFOLG) GOTO 30
C
C        UNTIL-1 : BIS EINE ZERLEGUNG ERFOLGEN DARF
C
         TEIL=TEIL+1
C
C        UMWANDELN DER ZAHL TEIL IN EINE ZEICHENKETTE TEILCH
C
         WRITE(TEILCH,'(I2)') TEIL
         IF (TEIL .LT. 10) THEN
            TEILCH(1:1)='0'
         ENDIF
C
C        ABSPEICHERN DER TEIL-TEN ELEMENTARZERLEGUNG
C
         IPART(TEIL)= 5 + ENDE+1-BEGINN
         PART(TEIL)=TEILCH //' = '//AUSDRU(BEGINN:ENDE) // ' '
C
C        ABSPEICHERN DER TEIL-TEN GENERATION VON AUSDRU
C
         IARITH(TEIL)=AKTLEN - (ENDE+1-BEGINN) +2
         IF (BEGINN .EQ. 1) then
            if (ENDE .EQ. AKTLEN) THEN
               ARITH(TEIL)=TEILCH
            ELSE
               ARITH(TEIL)=TEILCH//AUSDRU(ENDE+1:AKTLEN)
            ENDIF
         ELSEIF (ENDE .EQ. AKTLEN) THEN
               ARITH(TEIL)=AUSDRU(1:BEGINN-1)//TEILCH
         ELSE
            ARITH(TEIL)=AUSDRU(1:BEGINN-1)//TEILCH
     >                   //AUSDRU(ENDE+1:AKTLEN)
         ENDIF
         AKTLEN=IARITH(TEIL)
         AUSDRU=ARITH(TEIL)
C
C        WERTE ZURUECKSETZEN
C
         OMEGA= OMEGA - (ENDE+1-BEGINN) +2
         ERFOLG=.FALSE.
         BEGINN=alpha
         ENDE=omega
         I=BEGINN-1
         PREFIX=.FALSE.
C
C        SUCHE NACH DEM 1. OPERATOR IM TEILSTRING
C
50          CONTINUE
            I=I+1
            POS=INDEX(FAKTOR, AUSDRU(I:I))
         IF (POS .LE. 0 .AND. I .LT. OMEGA) GOTO 50

      GOTO 60
      ENDIF
C
C     ENDWHILE
C     ZERLEGUNG BEENDET, DA TEILSTRING KEIN OPERATOR MEHR ENTHAELT
C
C     ENDE VON SCHRIT
C
      END
C ===== SOURCE: second_own.f

c
      FUNCTION SECOND_OWN()
      USE PRECISION
      implicit none
      real(dp) x05baf,start,reset_second,second_own,time
      save start
      data start /0.0/
c x05baf is a NAG routine
!pb      time=x05baf()
      call cpu_time(time)
      second_own=time-start
      RETURN
C
      ENTRY RESET_SECOND()
!pb      start=x05baf()
      call cpu_time(start)
      reset_second=start
      return
      END
C ===== SOURCE: setref.f
C
C
      SUBROUTINE SETREF(AFF,AFFI,IFLAG,CC1,CC2,CC3)
C
C  SET REFLECTION MATRIX AFF AND INVERSE REFLECTION MATRIX AFFI = AFF
C  INPUT:
C  IFLAG=1:  REFLECTION HYPERPLANE NORMAL VECTOR C1,C2,C3
C            I.E. REFLECTION AT PLANE X*C1+Y*C2+Z*C3+0 (!!!)=0
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: AFF(3,3),AFFI(3,3)
      REAL(DP), INTENT(IN) :: CC1, CC2, CC3
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: C, C1, C2, C3, CN
      INTEGER :: I, J

      C1=CC1
      C2=CC2
      C3=CC3
      IF (IFLAG.EQ.1) THEN
C  NORMALIZE ROTATION AXIS
        C=C1*C1+C2*C2+C3*C3
        IF (C.LE.0.) THEN
          WRITE (iunout,*) 
     .    'WARNING: INVALID REFLEC. PLANE IN SUBR. SETREF'
          WRITE (iunout,*) 'NO REFLECTION CARRIED OUT'
          DO 1 J=1,3
            DO 1 I=1,3
              AFF(I,J)=0.
              AFFI(I,J)=0.
1         CONTINUE
          DO 2 J=1,3
            AFF(J,J)=1.
            AFFI(J,J)=1.
2         CONTINUE
          RETURN
        ENDIF
        CN=SQRT(C)
        C1=C1/CN
        C2=C2/CN
        C3=C3/CN
        AFF(1,1)=1.-2.*C1*C1
        AFF(2,1)=  -2.*C2*C1
        AFF(3,1)=  -2.*C3*C1
        AFF(1,2)=  -2.*C1*C2
        AFF(2,2)=1.-2.*C2*C2
        AFF(3,2)=  -2.*C3*C2
        AFF(1,3)=  -2.*C1*C3
        AFF(2,3)=  -2.*C2*C3
        AFF(3,3)=1.-2.*C3*C3
        DO 10 I=1,3
          DO 10 J=1,3
            AFFI(I,J)=AFF(I,J)
10      CONTINUE
      ELSE
      ENDIF
      RETURN
      END
C ===== SOURCE: setrot.f
C
C
      SUBROUTINE SETROT(AFF,AFFI,IFLAG,CC1,CC2,CC3,CC4)
C
C  SET ROTATION MATRIX AFF AND INVERSE ROTATION MATRIX AFFI
C  INPUT:
C  IFLAG=1:  ROTATION AXIS C1,C2,C3 AND ROTATION ANGLE C4 DEGREES
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: AFF(3,3),AFFI(3,3)
      REAL(DP), INTENT(IN) :: CC1, CC2, CC3, CC4
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: C, C1, C2, C3, C4, CAL, SAL, CN, ANG, PI
      INTEGER :: I, J

      DATA PI/3.141592654/

      C1=CC1
      C2=CC2
      C3=CC3
      C4=CC4
      IF (IFLAG.EQ.1) THEN
C  NORMALIZE ROTATION AXIS
        C=C1*C1+C2*C2+C3*C3
        IF (C.LE.0.) THEN
          WRITE (iunout,*) 
     .    'WARNING: INVALID ROTATION AXIS IN SUBR. SETROT'
          WRITE (iunout,*) 'NO ROTATION CARRIED OUT'
          DO 1 J=1,3
            DO 1 I=1,3
              AFF(I,J)=0.
              AFFI(I,J)=0.
1         CONTINUE
          DO 2 J=1,3
            AFF(J,J)=1.
            AFFI(J,J)=1.
2         CONTINUE
          RETURN
        ENDIF
        CN=SQRT(C)
        C1=C1/CN
        C2=C2/CN
        C3=C3/CN
        ANG=C4*PI/180.
        CAL=COS(ANG)
        SAL=SIN(ANG)
        AFF(1,1)=CAL+(1-CAL)*C1*C1
        AFF(2,1)=    (1-CAL)*C2*C1+SAL*C3
        AFF(3,1)=    (1-CAL)*C3*C1-SAL*C2
        AFF(1,2)=    (1-CAL)*C1*C2-SAL*C3
        AFF(2,2)=CAL+(1-CAL)*C2*C2
        AFF(3,2)=    (1-CAL)*C3*C2+SAL*C1
        AFF(1,3)=    (1-CAL)*C1*C3+SAL*C2
        AFF(2,3)=    (1-CAL)*C2*C3-SAL*C1
        AFF(3,3)=CAL+(1-CAL)*C3*C3
        CAL=COS(-ANG)
        SAL=SIN(-ANG)
        AFFI(1,1)=CAL+(1-CAL)*C1*C1
        AFFI(2,1)=    (1-CAL)*C2*C1+SAL*C3
        AFFI(3,1)=    (1-CAL)*C3*C1-SAL*C2
        AFFI(1,2)=    (1-CAL)*C1*C2-SAL*C3
        AFFI(2,2)=CAL+(1-CAL)*C2*C2
        AFFI(3,2)=    (1-CAL)*C3*C2+SAL*C1
        AFFI(1,3)=    (1-CAL)*C1*C3+SAL*C2
        AFFI(2,3)=    (1-CAL)*C2*C3-SAL*C1
        AFFI(3,3)=CAL+(1-CAL)*C3*C3
      ELSE
      ENDIF
      RETURN
      END
C ===== SOURCE: signok.f



C-----------------------------------------------------------------------
            SUBROUTINE SIGNOK(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DAS UNTERPROGRAMM UEBERPRUEFT AUSDRU AUF ZULAESSIGE ZEICHEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(33), PARAMETER ::
     P                  GULTIG='ABCDEFGHIJKLMNOPQRSTUVWXYZ+-*/^()'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     HILFSVARIABLE :
C
         INTEGER :: I


      I=1
C
C     WHILE : SOLANGE KEIN FEHLER AUFGETRETEN UND STRING NOCH NICHT
C             VOLLSTAENDIG DURCHLAUFEN
C
10    IF (ERROR .EQ. 0  .AND.  I .LE. AKTLEN) THEN
         IF (INDEX(GULTIG,AUSDRU(I:I)) .GT. 0) THEN
            I=I+1
         ELSE
C
C           UNGUELTIGES ZEICHEN IN AUSDRU GEFUNDEN
C
            ERROR=1
         ENDIF
         GOTO 10
      ENDIF
C
C     ENDWHILE
C
C     ENDE VON SIGNOK
C
      END
C ===== SOURCE: smach.f
C
      FUNCTION SMACH(I)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP)  SMACH
      IF(I.EQ.1) THEN
        smach=epsilon(1.d0)
      ELSEIF(I.EQ.2) THEN
        smach=tiny(1.d0)
      ELSEIF(I.EQ.3) THEN
        smach=huge(1.d0)
      ELSE
        WRITE(iunout,'(A)')
     @   ' SMACH CALLED WITH PARAMETER JOB NOT EQUAL TO 1,2 OR 3.'
      ENDIF
      RETURN
      END
C ===== SOURCE: sngl_poly.f
      function sngl_poly (cf, al, rcmin, rcmax, fpp, ifexmn, ifexmx) 
     .                   result(cou)

      use precision
      
      implicit none

      real(dp), intent(in) :: cf(9), fpp(6)
      real(dp), intent(in) :: al, rcmin, rcmax
      integer, intent(in) :: ifexmn, ifexmx
      real(dp) :: cou, fp(6), s01, s02, ds12, expo1, expo2, ccxm1, 
     .            ccxm2, extrap 
      integer :: ii, if8, ifex

      if ((ifexmn .ne. 0) .and. (al < rcmin)) then

C  ELAB BELOW MINIMUM ENERGY FOR FIT:

        FP = FPP
        ifex = ifexmn

C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN
        IF (IFEXMN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(SIGMA)
          S01=RCMIN
          S02=LOG(2._DP)+RCMIN
          DS12=S02-S01
          EXPO1=CF(9)
          EXPO2=CF(9)
          DO 1 II=1,8
            IF8=9-II
            EXPO1=EXPO1*S01+CF(IF8)
            EXPO2=EXPO2*S02+CF(IF8)
 1        CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FP(1)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FP(2)=      (CCXM2-CCXM1)/DS12
          FP(3)=0.D0
C
          IFEX=5
        ENDIF

        COU=EXTRAP(AL,IFEX,FP(1),FP(2),FP(3))
        cou = log(cou)

      elseif ((ifexmx .ne. 0) .and. (al > rcmax)) then

C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:

C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
        FP = FPP
        COU=EXTRAP(AL,IFEXMX,FP(4),FP(5),FP(6))
        cou = log(cou)

      else

C  ELAB WITHIN ENERGY RANGE OF FIT:

        cou = cf(9)

        do ii = 8, 1, -1
          cou = cou * al + cf(ii)
        end do

      end if

      return
      end function sngl_poly
C ===== SOURCE: spoint.f
C
C
      SUBROUTINE SPOINT (A0,A1,A2,A3,X,Y,Z,XR,YR,ZR,P,T)
C                                                     X+ LAMBDA * XR
C  SCHNITTPUNKT DER EBENE A0,A1,A2,A3 MIT DER GERADEN Y+ LAMBDA * YR
C                                                     Z+ LAMBDA * ZR
C  OUTPUT:  P(3): SCHNITTPUNKT
C           T:    LAMBDA
C
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A0, A1, A2, A3, X, Y, Z, XR, YR, ZR
      REAL(DP), INTENT(OUT) :: P(3), T
      REAL(DP) :: XNEN

      XNEN=A1*XR+A2*YR+A3*ZR
      IF (ABS(XNEN).LT.1.D-30) THEN
        P(1)=1.D60
        P(2)=1.D60
        P(3)=1.D60
        T=1.D60
      ELSE
        T=-(A0+A1*X+A2*Y+A3*Z)/XNEN
        P(1)=X+T*XR
        P(2)=Y+T*YR
        P(3)=Z+T*ZR
      ENDIF
      RETURN
      END
C ===== SOURCE: subtit.f



C-----------------------------------------------------------------------
               SUBROUTINE SUBTIT(AUSDRU, AKTLEN, HILFE)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ELIMINATION VON BLANKS UND SUBTITUTION VON '**' DURCH '^'
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         CHARACTER(*), INTENT(INOUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS

C
C     ELIMINATION VON BLANKS
C
      I=1
C
C     WHILE-1 : SOLANGE AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C
30    IF (I .LT. AKTLEN) THEN
C
C        ENTFERNEN VON BLANKS
C
         IF (AUSDRU(I:I) .EQ. ' ') THEN
            IF ( I .EQ. 1) THEN
               HILFE=AUSDRU(I+1:AKTLEN) // ' '
            ELSE
               HILFE=AUSDRU(1:I-1) // AUSDRU(I+1:AKTLEN) // ' '
            ENDIF
            AUSDRU=HILFE
            AKTLEN=AKTLEN-1
         ELSE
            I=I+1
         ENDIF
         GOTO 30
      ENDIF
C
C     ENDWHILE-1
C
C
C     SUBTITUTION VON '**' DURCH '^'
C
      POS=INDEX(AUSDRU,'**')
C
C     WHILE-2 : SOLANGE AUSDRU NOCH '**' ENTHAELT
C
10    IF ( POS .GT. 1 .AND. POS .LT. AKTLEN-1) THEN
C
C        ERSETZEN VON '**' DURCH '^'
C
         HILFE=AUSDRU(1:POS-1) // '^' // AUSDRU(POS+2:AKTLEN)
         AUSDRU=HILFE
         AKTLEN=AKTLEN-1
         POS=INDEX(AUSDRU,'**')
         GOTO 10
      ENDIF
C
C     ENDWHILE-2
C
C     ENDE VON SUBTIT
C
      END
C ===== SOURCE: swordr.f
C
C
      SUBROUTINE SWORDR (X,NX)
C
C  INVERTIEREN DER REIHENFOLGE DER WERTE AUF EINEM ARRAY:
C  X(1)-->X(NX)
C  X(2)-->X(NX-1)
C       .
C       .
C       .
C
      USE PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NX
      REAL(DP), INTENT(INOUT) :: X(NX)
      REAL(DP) :: H
      INTEGER :: IS, I

      IS=NX/2
      DO 1 I=1,IS
      H=X(I)
      X(I)=X(NX-I+1)
1     X(NX-I+1)=H
      RETURN
      END
C ===== SOURCE: swrdrm.f
C
C
      SUBROUTINE SWRDRM (F,NROW,IX,IY,ISW)

      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NROW, IX, IY, ISW
      REAL(SP), INTENT(INOUT) :: F(NROW,*)
      REAL(SP) :: H
      INTEGER :: I, IS, J

      IF (MOD(ISW,2).EQ.1) THEN
        IS=IX/2
        DO 1 I=1,IS
        DO 1 J=1,IY
          H=F(I,J)
          F(I,J)=F(IX-I+1,J)
1         F(IX-I+1,J)=H
      ENDIF
      IF (ISW.GE.2) THEN
        IS=IY/2
        DO 2 J=1,IS
        DO 2 I=1,IX
          H=F(I,J)
          F(I,J)=F(I,IY-J+1)
2         F(I,IY-J+1)=H
      ENDIF
C
      RETURN
      END
C ===== SOURCE: symet.f
C
C
      SUBROUTINE SYMET(ESTIM,NTAL,NRAD,NR1ST,NP2ND,NT3RD,NAD,N1,
     .                 LP,LT)
C
C  SYMMETRISE TALLIES
C
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: ESTIM(*)
      INTEGER, INTENT(IN) :: NAD(*),N1(*)
      INTEGER, INTENT(IN) :: NTAL, NRAD, NR1ST, NP2ND, NT3RD
      LOGICAL, INTENT(IN) :: LP,LT
      REAL(DP) :: SAV
      INTEGER :: IR, IP, IT, IND, I1, INDEX1, INDEX2, J1, J2, ITAL,
     .           NSYM, NSYH

      IF (LP) THEN
        NSYM=NP2ND
        NSYH=(NSYM-1)/2
        DO 30 ITAL=1,NTAL
          IND=NAD(ITAL)*NRAD
          DO 5 I1=1,N1(ITAL)
            DO 10 IR=1,NR1ST
            DO 10 IT=1,NT3RD
            DO 10 IP=1,NSYH
                  J1=IR+((IT-1)*NP2ND+IP-1)*NR1ST
                  J2=IR+((IT-1)*NP2ND+NSYM-IP-1)*NR1ST
                  INDEX1=IND+(J1-1)*N1(ITAL)+I1
                  INDEX2=IND+(J2-1)*N1(ITAL)+I1
                  SAV=(ESTIM(INDEX1)+ESTIM(INDEX2))*0.5
                  ESTIM(INDEX1)=SAV
                  ESTIM(INDEX2)=SAV
10          CONTINUE
5         CONTINUE
30      CONTINUE
      ENDIF
      IF (LT) THEN
        NSYM=NT3RD
        NSYH=(NSYM-1)/2
        DO 130 ITAL=1,NTAL
          IND=NAD(ITAL)*NRAD
          DO 105 I1=1,N1(ITAL)
            DO 110 IR=1,NR1ST
            DO 110 IP=1,NP2ND
            DO 110 IT=1,NSYH
                  J1=IR+((IT-1)*NP2ND+IP-1)*NR1ST
                  J2=IR+((NSYM-IT-1)*NP2ND+IP-1)*NR1ST
                  INDEX1=IND+(J1-1)*N1(ITAL)+I1
                  INDEX2=IND+(J2-1)*N1(ITAL)+I1
                  SAV=(ESTIM(INDEX1)+ESTIM(INDEX2))*0.5
                  ESTIM(INDEX1)=SAV
                  ESTIM(INDEX2)=SAV
110         CONTINUE
105        CONTINUE
130     CONTINUE
      ENDIF
      RETURN
      END
C ===== SOURCE: trmain.f
C
C
      SUBROUTINE TRMAIN(X,NSECND)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: X
      INTEGER, INTENT(IN) :: NSECND
      REAL(DP) :: S, SECOND_OWN
      S=SECOND_OWN( )
      X=FLOAT(NSECND)-S
      RETURN
      END
C ===== SOURCE: uppercase.f


      subroutine uppercase (zeile)

      IMPLICIT NONE
      character(*), INTENT(INOUT) :: zeile
      INTEGER :: L, I, J, LANF, LEND
      character(26) :: klein, gross
      data klein /'abcdefghijklmnopqrstuvwxyz'/
      data gross /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      LANF=verify(zeile,' ')
      LEND=verify(zeile,' ',.true.)
      if (lanf == 0) return
      l = lend-lanf+1
      if (l < lend) then
        zeile(1:l) = zeile(lanf:lend)
        zeile(l+1:lend) = repeat(' ',lanf)
      end if
      do i=1,l
        j=index(klein,zeile(i:i))
        if (j>0) zeile(i:i)=gross(j:j)
      end do

      return
      end subroutine uppercase

C ===== SOURCE: wrstrt.f
!pb  27.11.06: open and close statements for fort.10 moved here
C
C
      SUBROUTINE WRSTRT(IG,NSTRAI,IESTM1,IESTM2,IESTM3,
     .                  TALLYV,TALLYS,TALLYL,
     .                  ISDVI1,STAT1,ISDVI2,STAT2,
     .                  ISDVC1,SIGC,ISDVC2,SIGCS,
     .                  IBGKI,SIG_BGK,JBGKI,SIGS_BGK,
     .                  ICOPI,SIG_COP,JCOPI,SIGS_COP,
     .                  ISPCI,TRCFLE)

      USE PRECISION
      USE PARMMOD, ONLY: EIRENE_SPECTRUM, SPECT_ARRAY
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      TYPE(SPECT_ARRAY), INTENT(INOUT) :: TALLYL(*)
      REAL(DP), INTENT(INOUT) :: TALLYV(*), TALLYS(*),
     .                         STAT1(*), SIG_BGK(*), SIG_COP(*)
      REAL(DP), INTENT(INOUT) :: STAT2(*), SIGC(*), SIGCS(*)
      REAL(DP), INTENT(INOUT) :: SIGS_BGK(*), SIGS_COP(*)
      INTEGER, INTENT(IN) :: IG, NSTRAI, IESTM1, IESTM2, ISDVI1, ISDVI2,
     .                       ISDVC1, ISDVC2, IBGKI, JBGKI, ICOPI, JCOPI,
     .                       IESTM3, ISPCI
      LOGICAL, INTENT(IN) :: TRCFLE

      INTEGER :: IMAX11, IMAX12, IMAX21, IMAX22, IMAX23, IMAX24, IMAX2,
     .           IMAX31, IMAX32, IMAX41, IMAX42, NRECL, IRC, ISTRA,
     .           JINI, J, JEND, IMAX, ISPC, IMAXS
!      REAL(DP), DIMENSION(:), POINTER :: PSGM
C
C  WRITE DATA FOR SINGLE STRATA ON TEMP. FILE FT10
C
      NRECL=1500
      IMAX11=IESTM1/NRECL+1
      IMAX12=IESTM2/NRECL+1
      IMAX21=ISDVI1/NRECL+1
      IMAX22=ISDVI2/NRECL+1
      IMAX23=ISDVC1/NRECL+1
      IMAX24=ISDVC2/NRECL+1
      IMAX2=IMAX21+IMAX22+IMAX23+IMAX24
      IMAX31=IBGKI/NRECL+1
      IMAX32=JBGKI/NRECL+1
      IMAX41=ICOPI/NRECL+1
      IMAX42=JCOPI/NRECL+1
      IMAXS=0
      DO ISPC=1,IESTM3
        IMAXS=IMAXS+1
        IMAXS=IMAXS+TALLYL(ISPC)%PSPC%NSPC/NRECL+1
        IF (ISPCI.NE.0) THEN
          IMAXS=IMAXS+2*(TALLYL(ISPC)%PSPC%NSPC/NRECL+1)
        END IF
      END DO
      IMAX=IMAX11+IMAX12+IMAX2+IMAX31+IMAX32+IMAX41+IMAX42+IMAXS
      ISTRA=IG
      IRC=ISTRA*IMAX+1
      IF (TRCFLE.AND.IG.NE.0) WRITE (iunout,*) 'WRITE STRATUM NO. ',IG
      IF (TRCFLE.AND.IG.EQ.0) WRITE (iunout,*) 'WRITE SUM OVER STRATA '
C
      OPEN (UNIT=10,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NRECL,
     .      STATUS='UNKNOWN')

      JINI=1
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMV'
1     JEND=MIN0(JINI-1+NRECL,IESTM1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM1)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (TALLYV(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM1) GOTO 12
      JINI=JEND+1
      IRC=IRC+1
      GOTO 1

12    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMS'
      IRC=IRC+1
      JINI=1
11    JEND=MIN0(JINI-1+NRECL,IESTM2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM2)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (TALLYS(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM2) GOTO 2
      JINI=JEND+1
      IRC=IRC+1
      GOTO 11
C
2     CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 1'
      IRC=IRC+1
      JINI=1
3     JEND=MIN0(JINI-1+NRECL,ISDVI1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI1)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (STAT1(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI1) GOTO 21
      JINI=JEND+1
      IRC=IRC+1
      GOTO 3
C
21    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 2'
      IRC=IRC+1
      JINI=1
22    JEND=MIN0(JINI-1+NRECL,ISDVI2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI2)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (STAT2(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI2) GOTO 23
      JINI=JEND+1
      IRC=IRC+1
      GOTO 22
C
23    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 3'
      IRC=IRC+1
      JINI=1
24    JEND=MIN0(JINI-1+NRECL,ISDVC1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC1)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (SIGC(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC1) GOTO 25
      JINI=JEND+1
      IRC=IRC+1
      GOTO 24
C
25    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 4'
      IRC=IRC+1
      JINI=1
26    JEND=MIN0(JINI-1+NRECL,ISDVC2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC2)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (SIGCS(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC2) GOTO 4
      JINI=JEND+1
      IRC=IRC+1
      GOTO 26
C
4     CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS_BGK'
      IRC=IRC+1
      JINI=1
5     JEND=MIN0(JINI-1+NRECL,IBGKI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IBGKI)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (SIG_BGK(J),J=JINI,JEND)
      IF (JEND.EQ.IBGKI) GOTO 6
      JINI=JEND+1
      IRC=IRC+1
      GOTO 5
C
6     CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'SUM STATIS_BGK'
      IRC=IRC+1
      JINI=1
61    JEND=MIN0(JINI-1+NRECL,JBGKI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.JBGKI)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (SIGS_BGK(J),J=JINI,JEND)
      IF (JEND.EQ.JBGKI) GOTO 62
      JINI=JEND+1
      IRC=IRC+1
      GOTO 61
C
62    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS_COP'
      IRC=IRC+1
      JINI=1
7     JEND=MIN0(JINI-1+NRECL,ICOPI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ICOPI)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (SIG_COP(J),J=JINI,JEND)
      IF (JEND.EQ.ICOPI) GOTO 8
      JINI=JEND+1
      IRC=IRC+1
      GOTO 7
C
8     CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'SUM STATIS_COP'
      IRC=IRC+1
      JINI=1
9     JEND=MIN0(JINI-1+NRECL,JCOPI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.JCOPI)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10,REC=IRC) (SIGS_COP(J),J=JINI,JEND)
      IF (JEND.EQ.JCOPI) GOTO 13
      JINI=JEND+1
      IRC=IRC+1
      GOTO 9
C
13    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'SPECTRA'
      DO ISPC=1,IESTM3
        IRC=IRC+1
        WRITE (10,REC=IRC) TALLYL(ISPC)%PSPC%SPCMIN,
     .                     TALLYL(ISPC)%PSPC%SPCMAX,
     .                     TALLYL(ISPC)%PSPC%SPCDEL,
     .                     TALLYL(ISPC)%PSPC%SPCDELI,
     .                     TALLYL(ISPC)%PSPC%SPCINT,
     .                     TALLYL(ISPC)%PSPC%SGMS,
     .                     TALLYL(ISPC)%PSPC%STVS,
     .                     TALLYL(ISPC)%PSPC%EES,
     .                     TALLYL(ISPC)%PSPC%NSPC,
     .                     TALLYL(ISPC)%PSPC%ISPCTYP,
     .                     TALLYL(ISPC)%PSPC%ISPCSRF,
     .                     TALLYL(ISPC)%PSPC%IPRTYP,
     .                     TALLYL(ISPC)%PSPC%IPRSP,
     .                     TALLYL(ISPC)%PSPC%IMETSP
        DO JINI=1,TALLYL(ISPC)%PSPC%NSPC,NRECL
          IRC=IRC+1
          JEND=MIN(TALLYL(ISPC)%PSPC%NSPC, JINI+NRECL-1)
          WRITE (10,REC=IRC) (TALLYL(ISPC)%PSPC%SPC(J),J=JINI,JEND)
        END DO
        IF (ISPCI.NE.0) THEN
          DO JINI=1,TALLYL(ISPC)%PSPC%NSPC,NRECL
            IRC=IRC+1
            JEND=MIN(TALLYL(ISPC)%PSPC%NSPC, JINI+NRECL-1)
            WRITE (10,REC=IRC) (TALLYL(ISPC)%PSPC%SGM(J),J=JINI,JEND)
          END DO
          DO JINI=1,TALLYL(ISPC)%PSPC%NSPC,NRECL
            IRC=IRC+1
            JEND=MIN(TALLYL(ISPC)%PSPC%NSPC, JINI+NRECL-1)
            WRITE (10,REC=IRC) (TALLYL(ISPC)%PSPC%SDV(J),J=JINI,JEND)
          END DO
        END IF
      END DO

      CLOSE (UNIT=10)
C
      RETURN
C
      ENTRY RSTRT(IG,NSTRAI,IESTM1,IESTM2,IESTM3,
     .            TALLYV,TALLYS,TALLYL,
     .            ISDVI1,STAT1,ISDVI2,STAT2,
     .            ISDVC1,SIGC,ISDVC2,SIGCS,
     .            IBGKI,SIG_BGK,JBGKI,SIGS_BGK,
     .            ICOPI,SIG_COP,JCOPI,SIGS_COP,
     .            ISPCI,TRCFLE)
C
C  READ DATA FOR SINGLE STRATA OR SUM OVER STRATA FROM TEMP. FILE 10
C
      NRECL=1500
      IMAX11=IESTM1/NRECL+1
      IMAX12=IESTM2/NRECL+1
      IMAX21=ISDVI1/NRECL+1
      IMAX22=ISDVI2/NRECL+1
      IMAX23=ISDVC1/NRECL+1
      IMAX24=ISDVC2/NRECL+1
      IMAX2=IMAX21+IMAX22+IMAX23+IMAX24
      IMAX31=IBGKI/NRECL+1
      IMAX32=JBGKI/NRECL+1
      IMAX41=ICOPI/NRECL+1
      IMAX42=JCOPI/NRECL+1
      IMAXS=0
      DO ISPC=1,IESTM3
        IMAXS=IMAXS+1
        IMAXS=IMAXS+TALLYL(ISPC)%PSPC%NSPC/NRECL+1
        IF (ISPCI.NE.0) THEN
          IMAXS=IMAXS+2*(TALLYL(ISPC)%PSPC%NSPC/NRECL+1)
        END IF
      END DO
      IMAX=IMAX11+IMAX12+IMAX2+IMAX31+IMAX32+IMAX41+IMAX42+IMAXS
      ISTRA=IG
      IRC=ISTRA*IMAX+1
      IF (TRCFLE.AND.IG.NE.0) WRITE (iunout,*) 'READ STRATUM NO. ',IG
      IF (TRCFLE.AND.IG.EQ.0) WRITE (iunout,*) 'READ SUM OVER STRATA '

      OPEN (UNIT=10,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NRECL,
     .      STATUS='OLD')

C
      JINI=1
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMV'
10    JEND=MIN0(JINI-1+NRECL,IESTM1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM1)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (TALLYV(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM1) GOTO 15
      JINI=JEND+1
      IRC=IRC+1
      GOTO 10
C
15    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMS'
      JINI=1
      IRC=IRC+1
16    JEND=MIN0(JINI-1+NRECL,IESTM2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM2)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (TALLYS(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM2) GOTO 20
      JINI=JEND+1
      IRC=IRC+1
      GOTO 16
C
20    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 1'
      IRC=IRC+1
      JINI=1
30    JEND=MIN0(JINI-1+NRECL,ISDVI1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI1)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (STAT1(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI1) GOTO 31
      JINI=JEND+1
      IRC=IRC+1
      GOTO 30
C
31    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 2'
      IRC=IRC+1
      JINI=1
32    JEND=MIN0(JINI-1+NRECL,ISDVI2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI2)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                                      IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (STAT2(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI2) GOTO 33
      JINI=JEND+1
      IRC=IRC+1
      GOTO 32
C
33    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 3'
      IRC=IRC+1
      JINI=1
34    JEND=MIN0(JINI-1+NRECL,ISDVC1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC1)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (SIGC(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC1) GOTO 35
      JINI=JEND+1
      IRC=IRC+1
      GOTO 34
C
35    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 4'
      IRC=IRC+1
      JINI=1
36    JEND=MIN0(JINI-1+NRECL,ISDVC2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC2)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (SIGCS(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC2) GOTO 40
      JINI=JEND+1
      IRC=IRC+1
      GOTO 36
C
40    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS_BGK'
      IRC=IRC+1
      JINI=1
50    JEND=MIN0(JINI-1+NRECL,IBGKI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IBGKI)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (SIG_BGK(J),J=JINI,JEND)
      IF (JEND.EQ.IBGKI) GOTO 60
      JINI=JEND+1
      IRC=IRC+1
      GOTO 50
C
60    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'SUM STATIS_BGK'
      IRC=IRC+1
      JINI=1
65    JEND=MIN0(JINI-1+NRECL,JBGKI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.JBGKI)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (SIGS_BGK(J),J=JINI,JEND)
      IF (JEND.EQ.JBGKI) GOTO 66
      JINI=JEND+1
      IRC=IRC+1
      GOTO 65
C
66    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS_COP'
      IRC=IRC+1
      JINI=1
70    JEND=MIN0(JINI-1+NRECL,ICOPI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ICOPI)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (SIG_COP(J),J=JINI,JEND)
      IF (JEND.EQ.ICOPI) GOTO 80
      JINI=JEND+1
      IRC=IRC+1
      GOTO 70
C
80    CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'SUM STATIS_COP'
      IRC=IRC+1
      JINI=1
90    JEND=MIN0(JINI-1+NRECL,JCOPI)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.JCOPI)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10,REC=IRC) (SIGS_COP(J),J=JINI,JEND)
      IF (JEND.EQ.JCOPI) GOTO 100
      JINI=JEND+1
      IRC=IRC+1
      GOTO 90
C
100   CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'SPECTRA'
      DO ISPC=1,IESTM3
        IRC=IRC+1
        READ (10,REC=IRC) TALLYL(ISPC)%PSPC%SPCMIN,
     .                    TALLYL(ISPC)%PSPC%SPCMAX,
     .                    TALLYL(ISPC)%PSPC%SPCDEL,
     .                    TALLYL(ISPC)%PSPC%SPCDELI,
     .                    TALLYL(ISPC)%PSPC%SPCINT,
     .                    TALLYL(ISPC)%PSPC%SGMS,
     .                    TALLYL(ISPC)%PSPC%STVS,
     .                    TALLYL(ISPC)%PSPC%EES,
     .                    TALLYL(ISPC)%PSPC%NSPC,
     .                    TALLYL(ISPC)%PSPC%ISPCTYP,
     .                    TALLYL(ISPC)%PSPC%ISPCSRF,
     .                    TALLYL(ISPC)%PSPC%IPRTYP,
     .                    TALLYL(ISPC)%PSPC%IPRSP,
     .                    TALLYL(ISPC)%PSPC%IMETSP
        DO JINI=1,TALLYL(ISPC)%PSPC%NSPC,NRECL
          IRC=IRC+1
          JEND=MIN(TALLYL(ISPC)%PSPC%NSPC, JINI+NRECL-1)
          READ (10,REC=IRC) (TALLYL(ISPC)%PSPC%SPC(J),J=JINI,JEND)
        END DO
        IF (ISPCI.NE.0) THEN
          DO JINI=1,TALLYL(ISPC)%PSPC%NSPC,NRECL
            IRC=IRC+1
            JEND=MIN(TALLYL(ISPC)%PSPC%NSPC, JINI+NRECL-1)
            READ (10,REC=IRC) (TALLYL(ISPC)%PSPC%SGM(J),J=JINI,JEND)
          END DO
          DO JINI=1,TALLYL(ISPC)%PSPC%NSPC,NRECL
            IRC=IRC+1
            JEND=MIN(TALLYL(ISPC)%PSPC%NSPC, JINI+NRECL-1)
            READ (10,REC=IRC) (TALLYL(ISPC)%PSPC%SDV(J),J=JINI,JEND)
          END DO
        END IF
      END DO

      CLOSE (UNIT=10)
C
      RETURN
      END
C ===== SOURCE: xcordf.f
C
C
C*DK XCORDF
      SUBROUTINE XCORDF(RQ,E,ELQ,AI,BI,X,XS,Y,YS)
C
C   INTERSECTION OF LINE: AI*Y+BI*X=0 WITH ELLIPSE (X-E)**2+Y*Y/ELQ=RQ
C   (X,Y)=INTERSECTION-POINT IS RETURNED
C   IF NO INTERSECTION, X IS SET TO 1.D30
C   THE SOLUTION WITH LARGER X IS (X,Y), THE OTHER IS (XS,YS)
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: RQ, E, ELQ, AI, BI
      REAL(DP), INTENT(OUT) :: X, XS, Y, YS
      REAL(DP) :: T, A, AA, AAA, B

       IF (AI.EQ.0.) GOTO 2
         T=-BI/AI
         A=T*T/ELQ
         AA=A+1.
         AAA=1./AA
         B=RQ*AA-E*E*A
         IF (B.LE.0.) GOTO 3
1        B=SQRT(B)
         X=AAA*(E+B)
         Y=T*X
         XS=AAA*(E-B)
         YS=T*XS
         RETURN
2      X=0.
       XS=0.
       IF (RQ.LE.E*E) GOTO 3
       B=SQRT((RQ-E*E)*ELQ)
       Y=B
       YS=-B
       RETURN
C    NO INTERSECTION
3      X=1.D30
       XS=1.D30
       Y=0.
       YS=0.
       RETURN
       END
C ===== SOURCE: zerleg.f



C-----------------------------------------------------------------------
      SUBROUTINE ZERLEG(AUSDRU, AKTLEN, IPART, PART, IARITH, ARITH,
     >                  TEIL, HILFE, ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DAS UNTERPROGRAMM UBERPRUEFT DEN AUSDRUCK AUF REGELVERLETZUNG
C     UND ZERLEGT IHN IN EINZELOPERATIONEN, SO DASS DER AUSDRUCK
C     BEI SEQUENTIELLER ABARBEITUNG RICHTIG AUSGEWERTET WIRD
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     EIN/AUSGABEPARAMETER :
C
         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

         CHARACTER(*), INTENT(OUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER, INTENT(OUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(*), INTENT(OUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER, INTENT(OUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

         CHARACTER(*), INTENT(OUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

         INTEGER, INTENT(OUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     LOKALE VARIABLEN :
C
         INTEGER :: ALPHA
C           : POSITION DES ERSTEN ZEICHENS VOM TEILAUSDRUCK

         INTEGER :: OMEGA
C           : POSITION DES LEZTEN ZEICHENS VOM TEILAUSDRUCK

C
C     HILFSVARIABLEN :
C
      INTEGER :: POS


C
C     ELIMINATION VON BLANKS UND ERSETZEN VON '**' DURCH '^'
C
      CALL SUBTIT(AUSDRU, AKTLEN, HILFE)
C
C     UEBERPRUEFUNG VERSCHIEDENER REGELVERLETZUNGEN
C
      CALL FEHLER(AUSDRU, AKTLEN, ERROR)
      IF (ERROR .EQ. 0) THEN
C
C        BEGINN DER ZERLEGUNG
C
         TEIL=0

         ALPHA=1
         OMEGA=AKTLEN
         POS=INDEX( AUSDRU(1:AKTLEN), ')' )
C
C        WHILE : SOLANGE NOCH KLAMMERN IN AUSDRU VORHANDEN
C
22       IF (POS .GT. 0) THEN
C
C           BESTIMME LETZTES ZEICHEN IM ERSTEN INNERSTEN KLAMMERAUSDRUCK
C
            OMEGA=POS-1
C
C           BESTIMME ERSTES ZEICHEN IM ERSTEN INNERSTEN KLAMMERAUSDRUCK
C
            POS=0
11             POS=INDEX( AUSDRU(POS+1:OMEGA), '(' ) +POS
            IF ( INDEX ( AUSDRU(POS+1:OMEGA), '(' ) .GT. 0 ) GOTO 11
            ALPHA=POS+1
C
C           ZERLEGUNG DES KLAMMERAUSDRUCKS
C
            CALL SCHRIT(AUSDRU, AKTLEN, ALPHA, OMEGA, TEIL,
     >                  IPART, PART, IARITH, ARITH)
C
C           UEBERPRUEFUNG,OB DIE KLAMMERN WEGFALLEN KOENNEN
C
            IF (OMEGA+1-ALPHA .LE. 2
     >          .AND. INDEX('+-',AUSDRU(ALPHA:ALPHA)) .EQ. 0) THEN
C
C              ENTFERNUNG DER KLAMMERN
C
               IF (ALPHA-1 .EQ. 1) THEN
                  IF (OMEGA+1 .EQ. AKTLEN) THEN
                     HILFE=AUSDRU(ALPHA:OMEGA)
                  ELSE
                     HILFE=AUSDRU(ALPHA:OMEGA)//AUSDRU(OMEGA+2:AKTLEN)
                  ENDIF
               ELSEIF (OMEGA+1 .EQ. AKTLEN) THEN
                  HILFE=AUSDRU(1:ALPHA-2)//AUSDRU(ALPHA:OMEGA)
               ELSE
                  HILFE=AUSDRU(1:ALPHA-2)//AUSDRU(ALPHA:OMEGA)
     >                   //AUSDRU(OMEGA+2:AKTLEN)
               ENDIF
               AKTLEN=AKTLEN-2
               AUSDRU=HILFE
               IARITH(TEIL)=AKTLEN
               ARITH(TEIL)=AUSDRU
            ENDIF

            ALPHA=1
            OMEGA=AKTLEN
            POS=INDEX( AUSDRU(1:AKTLEN), ')' )
            GOTO 22
         ENDIF
C
C        ENDWHILE
C
C        AUSDRUCK ENTHAELT KEINE KLAMMERAUSDRUECKE MEHR
C          => ZERLEGUNG DES GESAMMTEN AUSDRUCKS
C
         CALL SCHRIT(AUSDRU, AKTLEN, ALPHA, OMEGA, TEIL,
     >               IPART, PART, IARITH, ARITH)

         IF (TEIL .EQ. 0) THEN
C
C           ES SIND KEINE ZERLEGUNGEN GEMACHT WORDEN, DA
C           DER AUSDRUCK NUR AUS EINEM OPERAND BESTEHT
C
            TEIL = 1
            PART(1) = 'Z01 = ' // AUSDRU
            IPART(1) = 6 + AKTLEN
            ARITH(1) = 'Z01'
            IARITH(1) = 3
         ELSE
C
C           ERSETZEN VON '^' DURCH '**' UND EINFUEGEN VON 'Z' FUER
C           ZERLEGUNG
C
            CALL RUKSUB(TEIL, IPART, PART, IARITH, ARITH, HILFE)
         ENDIF

      ENDIF
C
C     ENDE VON ZERLEG
C
      END
