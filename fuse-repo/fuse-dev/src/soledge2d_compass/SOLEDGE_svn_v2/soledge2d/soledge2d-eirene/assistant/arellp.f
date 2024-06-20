C
C
*DK ARELLP
      SUBROUTINE EIRENE_ARELLP
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
      USE EIRMOD_PRECISION
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
