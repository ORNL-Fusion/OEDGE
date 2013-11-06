c     -*-Fortran-*-
C@PROCESS OPT(1),VECTOR(LEV(0))
C
C
c slmod begin
      SUBROUTINE IZTAU (CRMI,CRMB,CION,RIZB,CIOPTA,cprint,nizs)
c
c      SUBROUTINE IZTAU (CRMI,CRMB,CION,RIZB,CIOPTA,cprint)
c slmod end
      IMPLICIT  none
C     INCLUDE   "PARAMS"
      include    'params'
C     INCLUDE   "CGEOM"
      include    'cgeom'
C     INCLUDE   "CADAS"
      include    'cadas'
C     INCLUDE   "CNOCO"
      include    'cnoco'
C     INCLUDE   "CIONIZ"
      include    'cioniz'
c
c     Temporary includes 
c
      include 'cedge2d'
      include 'temp' 
      integer   printiz

      INTEGER   CION,CIOPTA,cprint
c slmod begin
      INTEGER   NIZS
c slmod end
      REAL      CRMI,RIZB,crmb
C
C***********************************************************************
C
C      THIS VERSION BY C.M.FARRELL  APRIL 1989
C      ---------------------------------------
C         THIS SUBROUTINE CALCULATES EXPECTED IONISATION & RECOMBINATION
C      TIMES TO THE IONISATION STATES 0 TO CION
C      OF A REQUESTED ELEMENT.
C      THE IONISATION TIMES ARE CALCULATED FOR EACH OF A GIVEN
C      SET OF TEMPERATURES AND DENSITIES.
C      THEY ARE CALCULATED FROM THE EQUATIONS
C      T = 1/(NE.EXP(-KIS/TB).(TB/KIS)**(1/2).SUM)
C      SUM = SUM(ANS.(LOG10(TB/KIS))**N), N=0,5
C      NE IS THE BACKGROUND ELECTRON DENSITY SO NE = RNB * RIZB
C      KIS IS A FUNCTION OF ATOMIC NUMBER AND IONISATION STATE AND
C      ANS IS A FUNCTION OF ATOMIC NUMBER, IONISATION STATE AND N.
C      VALUES OF BOTH FOR A SELECTED SET OF ELEMENTS
C      ARE STORED IN THE ROUTINE.
C      THE RECOMBINATION TIMES ARE FOR ELECTRON-ION RECOMBINATION
C      AND ARE EXTRACTED FROM THE ABELS VAN MAANEN PACKAGE.
C
C         PARAMETERS -
C     CRMI   : IMPURITY MASS
c     crmb   : background plasma ion mass 
C     CION   : ELEMENT (BY ATOMIC NUMBER)
C     RIZB   : BACKGROUND ION IONISATION STATE ZB
C     CIOPTA : IONISATION OPTION CAN BE 0,1,2,3,4,5,6
C
C     MAXIZS : HIGHEST ELECTRON IONISATION TIMES CAN BE CALCULATED FOR
C
C         ERRORS  -
C      CION MUST BE AN AVAILABLE ELEMENT:
C                    WILL BE SET TO 6 (IE CARBON)
C
C***********************************************************************
C
      INTEGER   N,IN,IZ,IK,IR,ntemp
      PARAMETER (N=6)
      REAL      KIS(12,10),ANS(N,12,10),SZ(12,10)
      REAL      KTI,LOGKTI,SUM,DUMMY,DENOM
c
      real      rion,rrec,rcxr
      real      tmpte,tmpti,tmpne 
C
      INTEGER   ICODE
      INTEGER   IONH,IONHE,IONLI,IONBE,IONB,IONC,IONO,IONCR,IONFE,IONNI
      PARAMETER(IONH=1, IONHE=2, IONLI=3, IONBE=4, IONB=5)
      PARAMETER(IONC=6, IONO=7, IONCR=8, IONFE=9, IONNI=10)
C
      DATA  KIS(1,IONH)/
     H  13.60/
      DATA  (KIS(IZ,IONHE), IZ = 1, 2)/
     H  24.60,  54.42/
      DATA  (KIS(IZ,IONLI), IZ = 1, 3)/
     L   5.39,  75.64, 122.45/
      DATA  (KIS(IZ,IONBE), IZ = 1, 4)/
     B   9.32,  18.21, 153.89, 217.71/
      DATA  (KIS(IZ,IONB), IZ = 1, 5)/
     B   8.30,  25.15,  37.93, 259.37, 340.22/
      DATA  (KIS(IZ,IONC), IZ = 1, 6)/
     C  11.26,  24.38,  41.38,  64.49, 392.08, 489.98/
      DATA  (KIS(IZ,IONO), IZ = 1, 8)/
     O  13.62,  35.12,  54.93,  68.60, 103.68, 138.12, 739.32, 871.39/
      DATA  (KIS(IZ,IONCR), IZ = 1, 12)/
     C   6.80,  16.50,  31.00,  49.10,  69.30,  90.60, 161.10, 184.70,
     R                                 209.30, 244.40, 270.80, 298.00/
      DATA  (KIS(IZ,IONFE), IZ = 1, 12)/
     F   7.90,  16.20,  30.60,  54.80,  75.00,  99.00, 125.00, 151.10,
     E                                 233.60, 214.00, 290.30, 330.80/
      DATA  (KIS(IZ,IONNI), IZ = 1, 12)/
     N   7.60,  18.20,  36.20,  54.90,  76.10, 108.00, 133.00, 162.00,
     I                                 193.00, 224.60, 321.00, 352.00/
C
C    >   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
C    >                                   0.00,   0.00,   0.00,   0.00,
C
      DATA  (ANS(IN,1,IONH), IN = 1, N)/
     H   2.3742E-08, -3.6866E-09, -1.0366E-08,
     1               -3.8010E-09,  3.4159E-09,  1.6834E-09/
C
      DATA  (ANS(IN,1,IONHE), IN = 1, N)/
     H   1.4999E-08,  5.6656E-10, -6.0821E-09,
     1               -3.5894E-09,  1.5529E-09,  1.3207E-09/
      DATA  (ANS(IN,2,IONHE), IN = 1, N)/
     H   3.4356E-09, -1.6865E-09, -6.9236E-10,
     2                9.7863E-11,  1.5591E-10,  6.2236E-11/
C
      DATA  (ANS(IN,1,IONLI), IN = 1, N)/
     L   9.9655E-08, -5.5941E-08, -5.5228E-08,
     1                4.0589E-08,  1.4800E-08, -1.3120E-08/
      DATA  (ANS(IN,2,IONLI), IN = 1, N)/
     L   3.4023E-09, -7.6588E-10, -8.6078E-10,
     2               -8.9748E-10,  4.1661E-10,  3.3188E-10/
      DATA  (ANS(IN,3,IONLI), IN = 1, N)/
     L   1.1786E-09, -8.7637E-10, -9.3373E-11,
     3                2.1137E-10,  1.9017E-11, -4.0679E-11/
C
      DATA  (ANS(IN,1,IONBE), IN = 1, N)/ 7.4206E-08, -1.5520E-08,
     B      -3.9403E-08,  7.2155E-09,  1.1098E-08, -2.5501E-09/
      DATA  (ANS(IN,2,IONBE), IN = 1, N)/ 1.7136E-08, -1.3997E-08,
     B      -2.9656E-11,  3.0777E-09, -1.1676E-10, -5.1930E-10/
      DATA  (ANS(IN,3,IONBE), IN = 1, N)/ 1.6039E-09, -6.4336E-10,
     B      -7.7804E-10,  3.3527E-10,  2.1889E-10, -1.0600E-10/
      DATA  (ANS(IN,4,IONBE), IN = 1, N)/ 4.9587E-10, -3.6870E-10,
     B      -3.9284E-11,  8.8928E-11,  8.0007E-12, -1.7114E-11/
C
      DATA  (ANS(IN,1,IONB), IN = 1, N)/ 5.8365E-08, 1.0047E-08,
     B      -3.6230E-08, -7.3448E-09,  1.0220E-08,  1.6951E-09/
      DATA  (ANS(IN,2,IONB), IN = 1, N)/ 2.0590E-08, -9.8899E-09,
     B      -6.0949E-09,  2.7762E-09,  1.6499E-09, -6.7692E-10/
      DATA  (ANS(IN,3,IONB), IN = 1, N)/ 5.7023E-09, -4.6578E-09,
     B      -9.8688E-12,  1.0242E-09, -3.8855E-11, -1.7281E-10/
      DATA  (ANS(IN,4,IONB), IN = 1, N)/ 7.3539E-10, -2.9498E-10,
     B      -3.5672E-10,  1.5372E-10,  1.0036E-10, -4.8602E-11/
      DATA  (ANS(IN,5,IONB), IN = 1, N)/ 2.5458E-10, -1.8930E-10,
     B      -2.0169E-11,  4.5657E-11,  4.1076E-12, -8.7867E-12/
C
      DATA  (ANS(IN,1,IONC), IN = 1, N)/ 5.9848E-08, 1.1903E-08,
     C      -3.0140E-08, -1.3693E-08,  8.3748E-09,  4.0150E-09/
      DATA  (ANS(IN,2,IONC), IN = 1, N)/ 2.8395E-08, -1.6698E-08,
     C      -2.3557E-09,  3.2161E-10,  9.6016E-10,  5.2713E-10/
      DATA  (ANS(IN,3,IONC), IN = 1, N)/ 8.0945E-09, -3.6568E-09,
     C      -3.9572E-09,  2.3820E-09,  1.0515E-09, -7.9301E-10/
      DATA  (ANS(IN,4,IONC), IN = 1, N)/ 2.7464E-09, -2.0070E-09,
     C      -2.3595E-11,  4.2011E-10, -8.1600E-11, -3.9729E-11/
      DATA  (ANS(IN,5,IONC), IN = 1, N)/ 3.9495E-10, -1.5842E-10,
     C      -1.9158E-10,  8.2555E-11,  5.3899E-11, -2.6102E-11/
      DATA  (ANS(IN,6,IONC), IN = 1, N)/ 1.4715E-10, -1.0941E-10,
     C      -1.1657E-11,  2.6389E-11,  2.3742E-12, -5.0786E-12/
C
      DATA  (ANS(IN,1,IONO), IN = 1, N)/ 3.3559E-08,   1.3449E-08,
     O      -6.7112E-09, -1.9976E-08,  1.6214E-09,  6.5852E-09/
      DATA  (ANS(IN,2,IONO), IN = 1, N)/ 2.4476E-08, -5.3141E-09,
     O      -7.3316E-09, -4.4515E-09,  2.4257E-09,  1.9791E-09/
      DATA  (ANS(IN,3,IONO), IN = 1, N)/ 1.4741E-08, -8.7905E-09,
     O      -8.6099E-10, -2.4143E-10,  1.2598E-10,  6.4901E-10/
      DATA  (ANS(IN,4,IONO), IN = 1, N)/ 5.7525E-09, -2.3510E-09,
     O      -2.3387E-09,  1.4063E-09,  1.2451E-10, -1.9423E-10/
      DATA  (ANS(IN,5,IONO), IN = 1, N)/ 3.1263E-09, -1.6820E-09,
     O      -9.4803E-10,  8.4311E-11,  5.0188E-10, -4.1921E-11/
      DATA  (ANS(IN,6,IONO), IN = 1, N)/ 1.0099E-09, -6.5165E-10,
     O       2.8863E-12,  3.0336E-11, -1.4065E-11,  4.5106E-11/
      DATA  (ANS(IN,7,IONO), IN = 1, N)/ 1.5258E-10, -6.1203E-11,
     O      -7.4014E-11,  3.1894E-11,  2.0823E-11, -1.0084E-11/
      DATA  (ANS(IN,8,IONO), IN = 1, N)/ 6.2090E-11, -4.6167E-11,
     O      -4.9189E-12,  1.1135E-11,  1.0018E-12, -2.1430E-12/
C
      DATA  (ANS(IN,1,IONCR), IN = 1, N)/ 1.7476E-07, -9.4892E-08,
     C      -3.9137E-08,  2.2035E-08,  1.0542E-08, -4.5512E-09/
      DATA  (ANS(IN,2,IONCR), IN = 1, N)/ 1.6346E-07, -2.0193E-08,
     C      -3.3862E-08, -3.4535E-08, -7.6987E-09,  4.3482E-09/
      DATA  (ANS(IN,3,IONCR), IN = 1, N)/ 6.4591E-08, -1.3383E-08,
     C      -1.8749E-08, -2.5181E-09,  2.3497E-09, -3.5688E-09/
      DATA  (ANS(IN,4,IONCR), IN = 1, N)/ 2.6807E-08, -1.1048E-08,
     C      -1.0254E-08, -6.4356E-11,  6.9254E-09, -1.7988E-09/
      DATA  (ANS(IN,5,IONCR), IN = 1, N)/ 1.2173E-08, -6.2551E-09,
     C      -2.0154E-09,  2.6449E-09, -6.3679E-10, -3.5022E-10/
      DATA  (ANS(IN,6,IONCR), IN = 1, N)/ 9.5663E-09, -7.0338E-09,
     C       1.8930E-09, -1.1051E-09, -5.6458E-10,  8.6219E-10/
      DATA  (ANS(IN,7,IONCR), IN = 1, N)/ 6.0363E-09, -4.6514E-09,
     C      -4.8253E-10,  1.2458E-09,  8.7311E-11, -2.5985E-10/
      DATA  (ANS(IN,8,IONCR), IN = 1, N)/ 4.3538E-09, -3.4216E-09,
     C      -3.1101E-10,  9.2875E-10,  5.0125E-11, -1.9485E-10/
      DATA  (ANS(IN,9,IONCR), IN = 1, N)/ 3.0908E-09, -2.5123E-09,
     C      -1.3583E-10,  6.5966E-10,  1.1133E-11, -1.3365E-10/
      DATA  (ANS(IN,10,IONCR), IN = 1, N)/ 2.0484E-09, -1.6875E-09,
     C      -2.4637E-10,  6.1647E-10,  3.7796E-11, -1.5391E-10/
      DATA  (ANS(IN,11,IONCR), IN = 1, N)/ 1.4481E-09, -1.1805E-09,
     C      -1.0586E-10,  3.5501E-10,  1.3576E-11, -7.9488E-11/
      DATA  (ANS(IN,12,IONCR), IN = 1, N)/ 1.0909E-09, -8.5926E-10,
     C      -5.9555E-11,  2.1767E-10,  8.1838E-12, -4.3076E-11/
C
      DATA  (ANS(IN,1,IONFE), IN = 1, N)/    1.4438E-07, -8.0018E-08,
     F      -6.1752E-08,  5.6502E-08,  1.2350E-08, -1.7668E-08/
      DATA  (ANS(IN,2,IONFE), IN = 1, N)/    4.6220E-08, -2.4962E-08,
     F       2.8052E-09, -2.8845E-09,  4.2137E-09, -2.1823E-09/
      DATA  (ANS(IN,3,IONFE), IN = 1, N)/    3.2093E-08, -3.7148E-09,
     F      -1.5464E-08,  6.2761E-09,  1.0119E-09, -1.4006E-09/
      DATA  (ANS(IN,4,IONFE), IN = 1, N)/    2.6806E-08, -3.6843E-09,
     F      -5.1535E-09, -5.3097E-09, -1.5243E-09,  5.4680E-10/
      DATA  (ANS(IN,5,IONFE), IN = 1, N)/    1.7044E-08, -3.5184E-09,
     F      -4.8036E-09, -8.1005E-10,  5.3660E-10, -8.6145E-10/
      DATA  (ANS(IN,6,IONFE), IN = 1, N)/    1.6784E-08, -1.1368E-08,
     F      -1.0575E-09,  1.3537E-09,  3.8511E-10, -3.5680E-11/
      DATA  (ANS(IN,7,IONFE), IN = 1, N)/    1.1951E-08, -9.0125E-09,
     F       7.8484E-10,  4.2869E-10, -1.3924E-10,  2.2839E-10/
      DATA  (ANS(IN,8,IONFE), IN = 1, N)/    7.0219E-09, -5.1840E-09,
     F       6.7389E-10, -5.0196E-11, -2.2030E-10,  3.1117E-10/
      DATA  (ANS(IN,9,IONFE), IN = 1, N)/    3.4559E-09, -2.6630E-09,
     F      -2.7626E-10,  7.1325E-10,  4.9986E-11, -1.4878E-10/
      DATA  (ANS(IN,10,IONFE), IN = 1, N)/   3.2349E-09, -1.6688E-09,
     F      -9.1742E-10,  2.3160E-10,  2.0134E-10, -1.1015E-10/
      DATA  (ANS(IN,11,IONFE), IN = 1, N)/   1.8841E-09, -1.5320E-09,
     F      -8.2512E-11,  4.0239E-10,  6.6853E-12, -8.1537E-11/
      DATA  (ANS(IN,12,IONFE), IN = 1, N)/   1.3008E-09, -1.0716E-09,
     F      -1.5646E-10,  3.9149E-10,  2.4002E-11, -9.7738E-11/
C
      DATA  (ANS(IN,1,IONNI), IN = 1, N)/   1.1655E-07, -3.9394E-08,
     N      -6.9896E-08,  3.4807E-08,  1.7806E-08, -1.1861E-08/
      DATA  (ANS(IN,2,IONNI), IN = 1, N)/   4.6772E-08, -1.9826E-08,
     N      -6.6507E-09,  1.2882E-08, -6.5867E-09, -1.1212E-09/
      DATA  (ANS(IN,3,IONNI), IN = 1, N)/   3.1226E-08, -1.2061E-08,
     N      -5.7631E-09, -1.2886E-09,  1.7137E-09,  8.7919E-10/
      DATA  (ANS(IN,4,IONNI), IN = 1, N)/   2.4688E-08, -7.2247E-09,
     N      -2.3266E-09, -6.6474E-10, -1.3530E-09, -1.0559E-09/
      DATA  (ANS(IN,5,IONNI), IN = 1, N)/   8.1992E-09, -9.4906E-10,
     N      -3.9507E-09,  1.6034E-09,  2.5852E-10, -3.5782E-10/
      DATA  (ANS(IN,6,IONNI), IN = 1, N)/   9.4899E-09, -1.5693E-09,
     N      -1.4282E-09, -1.7224E-09, -7.9560E-10,  1.3403E-10/
      DATA  (ANS(IN,7,IONNI), IN = 1, N)/   7.9357E-09, -1.3089E-09,
     N      -3.5502E-09,  1.2384E-10,  1.0690E-09, -7.2708E-10/
      DATA  (ANS(IN,8,IONNI), IN = 1, N)/   4.4730E-09, -1.8453E-09,
     N      -1.7110E-09, -1.0738E-11,  1.1556E-09, -3.0014E-10/
      DATA  (ANS(IN,9,IONNI), IN = 1, N)/   2.6192E-09, -1.3458E-09,
     N      -4.3365E-10,  5.6908E-10, -1.3701E-10, -7.5354E-11/
      DATA  (ANS(IN,10,IONNI), IN = 1, N)/ 2.4493E-09, -1.8009E-09,
     N       4.8466E-10, -2.8294E-10, -1.4455E-10,  2.2075E-10/
      DATA  (ANS(IN,11,IONNI), IN = 1, N)/ 2.1454E-09, -1.6532E-09,
     N      -1.7150E-10,  4.4279E-10,  3.1031E-11, -9.2362E-11/
      DATA  (ANS(IN,12,IONNI), IN = 1, N)/ 1.6548E-09, -1.3005E-09,
     N      -1.1821E-10,  3.5301E-10,  1.9052E-11, -7.4059E-11/
C
C     DATA  (ANS(IN,1,IONNI), IN = 1, N)/ 0.0000E-00,   0.0000E-00,
C    O       0.0000E-00,  0.0000E-00,  0.0000E-00,  0.0000E-00/
C
       DATA (SZ(IN,IONH ), IN=1, 1) / 3.142E-8 /
       DATA (SZ(IN,IONHE), IN=1, 2) / 2.547E-8, 3.924E-9 /
       DATA (SZ(IN,IONLI), IN=1, 3) / 9.002E-8, 4.679E-9, 1.136E-9 /
       DATA (SZ(IN,IONBE), IN=1, 4) / 1.015E-7, 1.586E-8, 1.795E-9,
     1      4.796E-10 /
       DATA (SZ(IN,IONB ), IN=1, 5) / 1.111E-7, 2.402E-8, 5.281E-9,
     1      8.196E-10, 2.457E-10 /
       DATA (SZ(IN,IONC ), IN=1, 6) / 1.244E-7, 3.201E-8, 8.936E-9,
     1      2.912E-9, 4.444E-10, 1.424E-10 /
       DATA (SZ(IN,IONO ), IN=1, 8) / 9.343E-8, 3.387E-8, 1.607E-8,
     1      6.867E-9, 3.024E-9, 1.220E-9, 1.704E-10, 5.996E-11 /
       DATA (SZ(IN,IONCR), IN=1,12) / 1.970E-7, 2.040E-7, 8.240E-8,
     1      3.020E-8, 1.630E-8, 1.040E-8, 5.610E-9, 4.010E-9,
     2      2.800E-9, 1.770E-9, 1.290E-9, 1.000E-9 /
       DATA (SZ(IN,IONFE), IN=1,12) / 1.540E-7, 7.210E-8, 5.780E-8,
     1      3.340E-8, 2.170E-8, 1.730E-8, 1.200E-8, 7.290E-9,
     2      3.230E-9, 2.790E-9, 1.700E-9, 1.120E-9 /
       DATA (SZ(IN,IONNI), IN=1,12) / 1.380E-7, 7.850E-8, 4.310E-8,
     1      3.490E-8, 1.480E-8, 1.190E-8, 1.020E-8, 5.050E-9,
     2      3.480E-9, 2.640E-9, 1.990E-9, 1.520E-9 /
C
C-----------------------------------------------------------------------
C     USE ADAS/NOCORONA PACKAGE TO CALCULATE RATES IF CIOPTA=3,4,5,6
C-----------------------------------------------------------------------
C
       IF (CIOPTA.EQ.3.OR.CIOPTA.EQ.4.OR.CIOPTA.EQ.5.OR.CIOPTA.EQ.6)
     >   GOTO 600
C
C-----------------------------------------------------------------------
C      FOR CIOPTA=1,2 IDENTIFY REQUESTED ELEMENT AND SET ICODE
C-----------------------------------------------------------------------
C
      IF      (CION .EQ. 1) THEN
         ICODE = IONH
      ELSE IF (CION .EQ. 2) THEN
         ICODE = IONHE
      ELSE IF (CION .EQ. 3) THEN
         ICODE = IONLI
      ELSE IF (CION .EQ. 4) THEN
         ICODE = IONBE
      ELSE IF (CION .EQ. 5) THEN
         ICODE = IONB
      ELSE IF (CION .EQ. 6) THEN
         ICODE = IONC
      ELSE IF (CION .EQ. 8) THEN
         ICODE = IONO
      ELSE IF (CION .EQ. 24) THEN
         ICODE = IONCR
      ELSE IF (CION .EQ. 26) THEN
         ICODE = IONFE
      ELSE IF (CION .EQ. 28) THEN
         ICODE = IONNI
      ELSE
         WRITE(6, *)  ' ERROR - SUBROUTINE IZTAU '
         WRITE(6, *)  ' IONISATION TIMES REQUESTED'
         WRITE(6, *)  ' OF ELEMENT ATOMIC NUMBER ', CION
         WRITE(6, *)  ' THESE ARE NOT AVAILABLE, CARBON ASSUMED'
         CION = 6
         ICODE = IONC
      ENDIF
550   CONTINUE
C
C-----------------------------------------------------------------------
C     USE MAXIMUM VALUES FOR IONISATION RATES IF CIOPTA=2.
C-----------------------------------------------------------------------
C
       IF (CIOPTA.EQ.2) GOTO 300
C
C-----------------------------------------------------------------------
C  STANDARD CASE:     FOR EACH IONISATION LEVEL
C-----------------------------------------------------------------------
C
      DO 200 IZ = 0, CION-1
       DO 100 IR = 1, NRS
         DO 100 IK = 1, NKS(IR)
          KTI = KTEBS(IK,IR) / KIS(IZ+1,ICODE)
C
C-------- IONISATION TIME VERY LARGE ...
C
          IF (KTI .LT. 0.02)  THEN
               KFIZS(IK,IR,IZ) = HI
C
C-------- SUM SERIES IN RIZ EQUATION
C-------- (LOOP IS EXPANDED OUT FOR OPTIMIZATION)
C
          ELSE
            LOGKTI = LOG10(KTI)
            SUM    = ANS(1,IZ+1,ICODE)
            DUMMY  = LOGKTI
            SUM    = SUM   + ANS(2,IZ+1,ICODE) * DUMMY
            DUMMY  = DUMMY * LOGKTI
            SUM    = SUM   + ANS(3,IZ+1,ICODE) * DUMMY
            DUMMY  = DUMMY * LOGKTI
            SUM    = SUM   + ANS(4,IZ+1,ICODE) * DUMMY
            DUMMY  = DUMMY * LOGKTI
            SUM    = SUM   + ANS(5,IZ+1,ICODE) * DUMMY
            DUMMY  = DUMMY * LOGKTI
            SUM    = SUM   + ANS(6,IZ+1,ICODE) * DUMMY
C
C---------- CALCULATE EXPECTED IONISATION TIME
C---------- THE 1.E6 FACTOR IS BECAUSE ANS IS IN CM**3/S
C
            DENOM = SQRT(KTI) * SUM * KNBS(IK,IR) * RIZB
            IF (DENOM.NE.0.0) THEN
              KFIZS(IK,IR,IZ) = 1.E6 * EXP(1.0/KTI) / DENOM
            ELSE
              KFIZS(IK,IR,IZ) = HI
            ENDIF
          ENDIF
          KFRCS(IK,IR,IZ) = 0.0
  100   CONTINUE
  200 CONTINUE
C
      GOTO 999
C
C-----------------------------------------------------------------------
C         CIOPTA = 2   USE MAXIMUM VALUES
C     CALCULATE IONISATION TIMES T = 1/(NE * SZ)
C-----------------------------------------------------------------------
C
  300  CONTINUE
C
C---- NE AND NOT NB IS USED IN THIS EXPRESSION, WHERE NE = IZ * NB
C
      DO 410 IZ = 0, CION-1
       DO 400 IR = 1, NRS
         DO 400 IK = 1, NKS(IR)
          KFIZS(IK,IR,IZ) = 1.0E6/(KNBS(IK,IR)*RIZB*SZ(IZ+1,ICODE))
          KFRCS(IK,IR,IZ) = 0.0
  400   CONTINUE
  410 CONTINUE
C
      GOTO 999
C
C-----------------------------------------------------------------------
C    NOCORONA INTERFACE.  EXTRACT IONISATION AND RECOMBINATION TIMES
C    FROM THE PACKAGE HELD IN JETPPS.ATOMIC.FORT(NOCORONA)
C    CIOPTA=3  IONISATION AND RECOMBINATION
C    CIOPTA=4  IONISATION ONLY
C    CIOPTA=5  AS 3, WITH IONISATION DISABLED AFTER GIVEN STATE
C    CIOPTA=6  AS 4, WITH IONISATION DISABLED AFTER GIVEN STATE
C    (OPTIONS 5,6 USEFUL FOR FOLLOWING NICKEL - DON'T NEED TO STORE 28
C     IONISATION STATES WORTH OF DATA)
C-----------------------------------------------------------------------
C
 
  600 CONTINUE
c
      CALL RZERO(KFIZS, MAXNKS*MAXNRS*(MAXIZS+1))
      CALL RZERO(KFRCS, MAXNKS*MAXNRS*(MAXIZS+1))
c
c     NOCORONA package specified.
c
      if (cdatopt.eq.0) then
c
      DO 700 IR = 1, NRS
C
C---- THE INITIALISATION OF NOCORONA IS DONE IN RUNLIM3 ...
C---- CALCULATE ELECTRON DENSITIES IN CM-3 FROM ION DENSITIES IN M-3.
C
      DO 610 IK = 1, NKS(IR)
        PNES(IK) = KNBS(IK,IR) * 1.E-6 * RIZB
        PTES(IK) = KTEBS(IK,IR)
  610 CONTINUE
C
C---- CALCULATE IONISATION AND RECOMBINATION RATES ...
C---- VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY GOES
C---- OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.
C
c      CALL RRATES (PTES, PNES, PRATES, NKS(IR))
c
      CALL NCRRATES ( NKS(IR))
c
C
      DO 630 IZ = 0, CION
        DO 620 IK = 1, NKS(IR)
          IF (IZ.EQ.CION) THEN
            KFIZS(IK,IR,IZ) = 0.0
          ELSE
            KFIZS(IK,IR,IZ) = 1.E6 /
     >        (PRATES(1,IZ+1,1,IK)*KNBS(IK,IR)*RIZB)
          ENDIF
          IF ((CIOPTA.EQ.3.OR.CIOPTA.EQ.5).AND.IZ.GT.0) THEN
            KFRCS(IK,IR,IZ) = 1.E6 /
     >        (PRATES(2,IZ,1,IK)*KNBS(IK,IR)*RIZB)
          ELSE
            KFRCS(IK,IR,IZ) = 0.0
          ENDIF
  620   CONTINUE
  630 CONTINUE
C
  700 CONTINUE
C
 9100 FORMAT(//1X,'ELEMENT WITH MASS ',G11.4,' IS NOT INCLUDED',
     >    ' IN THE NOCORONA PACKAGE.',/)
c
c     ADAS package specified
c
      elseif (cdatopt.eq.1) then
c slmod begin
        if (sloutput) then
          WRITE(0,*) 
          WRITE(0,*) '-------------------------------------------------'
          WRITE(0,*) 'WARNING IZTAU: CION->NIZS in places -SL, 23/11/11'
          WRITE(0,*) '-------------------------------------------------'
          WRITE(0,*) 
        endif
c slmod end
c
C
C---- DO ONE RING AT A TIME
C
      DO 800 IR = 1, NRS
C
      DO 710 IK = 1, NKS(IR)
        PNESA(IK) = KNBS(IK,IR) * RIZB
        PTESA(IK) = KTEBS(IK,IR)
  710 CONTINUE
C
C---- CALCULATE IONISATION AND RECOMBINATION RATES ...
C
      write(year,'(i2.2)') iyearz
      call xxuid(useridz)            
c      YEAR = '89'
c      YEARDF = '89'
      ICLASS = 2
C
c slmod begin
      DO 730 IZ = 0, NIZS-1
c
c      DO 730 IZ = 0, CION-1
c slmod end
        CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,   
     >              PCOEF(1,IZ+1))                                      
c
c        CALL ADASRD(YEAR,YEARDF,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
c     >              PCOEF(1,IZ+1))
        DO 720 IK = 1, NKS(IR)
          IF (PCOEF(IK,IZ+1).GT.0.0) THEN
            KFIZS(IK,IR,IZ) = 1.0 / 
     >                 (adas_iz_rate_mult*PCOEF(IK,IZ+1)*PNESA(IK))
          ELSE
            KFIZS(IK,IR,IZ) = HI
          ENDIF
c
  720   CONTINUE
  730 CONTINUE
C
      IF (CIOPTA.EQ.3 .OR. CIOPTA.EQ.5) THEN
        ICLASS = 1
C
c slmod begin
        DO 735 IZ = 1, NIZS
c
c        DO 735 IZ = 1, CION
c slmod end
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,   
     >                PCOEF(1,IZ))                                      
c          CALL ADASRD(YEAR,YEARDF,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
c     >                PCOEF(1,IZ))
          DO 725 IK = 1, NKS(IR)
            IF (PCOEF(IK,IZ).GT.0.0) THEN
              KFRCS(IK,IR,IZ) = 1.0 / 
     >                 (adas_rec_rate_mult*PCOEF(IK,IZ)*PNESA(IK))
            ELSE
              KFRCS(IK,IR,IZ) = HI
            ENDIF
c
c
  725     CONTINUE
  735   CONTINUE
      ENDIF
C
  800 CONTINUE
c
c     ADPAK atomic physics style data - b2frates tables
c
      elseif (cdatopt.eq.2) then
c
c        Loop through background plasma
c
         do ir = 1,nrs
c    
            do ik = 1,nks(ir)
c       
               do iz = 0,cion
c
c                 Scale the ion temperature by mass in AMU 
c
                  tmpne = knbs(ik,ir) * rizb                  
                  tmpti = ktibs(ik,ir) / crmb
                  tmpte = ktebs(ik,ir)  
c                  
                  call mcrates(tmpte,tmpti,tmpne,
     >                         iz,cion,cion,rion,rrec,rcxr,0)
c
c                 Check rates
c
                  if (check_rates) then 

                     if (iz.eq.1) 
     >                     rate_calc(ik,ir,1) = 
     >                    (rrec * e2dnes(ik,ir)) * e2dnzs(ik,ir,1) 
     >                           * karea2(ik,ir)

                     write(6,'(a,3i4,1p,10(1x,g10.4))') 'Rate1:',ik,ir,
     >                 iz, e2dnzs(ik,ir,iz),e2dnes(ik,ir),
     >                 rion,rrec,
     >                 (rion*e2dnes(ik,ir))*e2dnzs(ik,ir,iz),
     >                 (rrec*e2dnes(ik,ir))*e2dnzs(ik,ir,iz)



                  endif




c
c                  if (cprint.eq.3.or.cprint.eq.9) then 
c                     write (6,*) 'MCRATES:',ik,ir,iz,rion,rrec,rcxr
c                  endif
c
c                 Ionization
c
                  IF (rion.GT.0.0.and.tmpne.gt.0.0) THEN
                     KFIZS(IK,IR,IZ) = 1.0/(rion*tmpne)
                  ELSE
                     KFIZS(IK,IR,IZ) = HI
                  ENDIF
c
c                 Recombination
c
                  if (ciopta.eq.3.or.ciopta.eq.5) then 
c                   
                     IF (rrec.GT.0.0.and.tmpne.gt.0.0) THEN
                        KFRCS(IK,IR,IZ)=1.0/(rrec*tmpne)
                     ELSE
                        KFRCS(IK,IR,IZ) = HI
                     ENDIF
c
                  endif
c
               end do

            end do
       
         end do                   
c
c     INEL atomic physics data 
c
      elseif (cdatopt.eq.3) then
c
c        Loop through background plasma
c
         do ir = 1,nrs
c    
            do ik = 1,nks(ir)
c       
               do iz = 0,cion
c
c                 Scale the ion temperature by mass in AMU 
c
                  tmpne = knbs(ik,ir) * rizb                  
                  tmpti = ktibs(ik,ir) / crmb
                  tmpte = ktebs(ik,ir)  
c                  
c                 Call for Te-based quantities 
c
                  call imprates(tmpte,iz,cion,rion,rrec,rcxr)
c
c
c                 Check rates
c
                  if (check_rates) then 

                     if (iz.eq.1) 
     >                     rate_calc(ik,ir,1) = 
     >                    (rrec * e2dnes(ik,ir)) * e2dnzs(ik,ir,1) 
     >                           * karea2(ik,ir)

                     write(6,'(a,3i4,1p,10(1x,g10.4))') 'Rate1:',ik,ir,
     >                 iz, e2dnzs(ik,ir,iz),e2dnes(ik,ir),
     >                 rion,rrec,
     >                 (rion*e2dnes(ik,ir))*e2dnzs(ik,ir,iz),
     >                 (rrec*e2dnes(ik,ir))*e2dnzs(ik,ir,iz)



                  endif




c
c                  if (cprint.eq.3.or.cprint.eq.9) then 
c                     write (6,*) 'INEL:',ik,ir,iz,rion,rrec,rcxr
c                  endif
c
c                 Ionization
c
                  IF (rion.GT.0.0.and.tmpne.gt.0.0) THEN
                     KFIZS(IK,IR,IZ) = 1.0/(rion*tmpne)
                  ELSE
                     KFIZS(IK,IR,IZ) = HI
                  ENDIF
c
c                 Recombination
c
                  if (ciopta.eq.3.or.ciopta.eq.5) then 
c                   
                     IF (rrec.GT.0.0.and.tmpne.gt.0.0) THEN
                        KFRCS(IK,IR,IZ)=1.0/(rrec*tmpne)
                     ELSE
                        KFRCS(IK,IR,IZ) = HI
                     ENDIF
c
                  endif
c
               end do

            end do
       
         end do                   
c
c     End of CDATOPT BLOCK - NOCORONA/ADAS/ADPAK
c
      endif
c
C
C-----------------------------------------------------------------------
C         FINISHED: VALIDATE IONISATION DATA BEFORE LEAVING
C-----------------------------------------------------------------------
C
  999 CONTINUE
c slmod begin
      DO 650 IZ = 0, NIZS
c
c      DO 650 IZ = 0, CION
c slmod end
       DO 640 IR = 1, NRS
        DO 640 IK = 1, NKS(IR)
          IF (KFIZS(IK,IR,IZ).LT.0.0 .OR. KFIZS(IK,IR,IZ).GT.HI)
     >        KFIZS(IK,IR,IZ) = HI
          IF (KFRCS(IK,IR,IZ).LT.0.0 .OR. KFRCS(IK,IR,IZ).GT.HI)
     >        KFRCS(IK,IR,IZ) = HI
  640   CONTINUE
  650 CONTINUE
C
c     Print out a sampling of the atomic physics cross-sections if
c     the print option has been selected.
c
      if (cprint.eq.3.or.cprint.eq.9) then 

c
      write (6,*) 'PRINT RATE DATA: ',cdatopt
      write (6,*)

      
c      printiz = cion-1

      printiz = 1


      if (cdatopt.eq.1) then 

         ntemp = 30

         PTESA(1) = 0.1
         PTESA(2) = 0.25
         PTESA(3) = 0.50
         PTESA(4) = 0.75
         PTESA(5) = 1.0
         PTESA(6) = 1.25
         PTESA(7) = 1.50
         PTESA(8) = 1.75
         PTESA(9) = 2.0
         PTESA(10) = 2.25
         PTESA(11) = 2.5
         PTESA(12) = 2.75
         PTESA(13) = 3.0
         PTESA(14) = 4.0
         PTESA(15) = 5.0
         PTESA(16) = 6.0
         PTESA(17) = 7.0
         PTESA(18) = 8.0
         PTESA(19) = 9.0 
         PTESA(20) = 10.0
         PTESA(21) = 20.0
         PTESA(22) = 30.0
         PTESA(23) = 40.0
         PTESA(24) = 50.0
         PTESA(25) = 60.0
         PTESA(26) = 70.0
         PTESA(27) = 80.0
         PTESA(28) = 90.0
         PTESA(29) = 100.0
         PTESA(30) = 150.0

         do in = 1,4

            tmpne = 1.0e17 * 10**(in) 

            do ik = 1,ntemp
 
               pnesa(ik) = tmpne

            end do 
C
C---- CALCULATE IONISATION AND RECOMBINATION RATES ...
C
            write(year,'(i2.2)') iyearz
            call xxuid(useridz)            
c
            ICLASS = 2
C
            DO IZ = 0, min(printiz,cion-1)
               CALL ADASRD(YEAR,CION,IZ+1,ICLASS,Ntemp,PTESA,PNESA,   
     >              PCOEF(1,IZ+1))                                      

            end do
c
            write (6,*) 'ADAS Ionization rate coefficients:'//
     >                  ' Density = ', tmpne
c
            write(6,1001) (iz,iz=0,min(printiz,cion-1))
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=0,min(printiz,cion-1))
            end do
c
C
c           RECOMBINATION
c
            iclass = 1
c
            DO IZ = 1, min(printiz+1,cion)
               CALL ADASRD(YEAR,CION,IZ,ICLASS,Ntemp,PTESA,PNESA,   
     >              PCOEF(1,IZ))                                      

            end do
c
            write (6,*) 'ADAS Recombination rate coefficients:'//
     >                  ' Density = ', tmpne
c
            write(6,1001)  (iz,iz=1,min(printiz+1,cion))
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz),iz=1,min(printiz+1,cion))
            end do
c
c        End density loop
c
         end do

c
c     ADPAK atomic physics style data - b2frates tables
c
      elseif (cdatopt.eq.2) then

         ntemp = 30

         PTESA(1) = 0.1
         PTESA(2) = 0.25
         PTESA(3) = 0.50
         PTESA(4) = 0.75
         PTESA(5) = 1.0
         PTESA(6) = 1.25
         PTESA(7) = 1.50
         PTESA(8) = 1.75
         PTESA(9) = 2.0
         PTESA(10) = 2.25
         PTESA(11) = 2.5
         PTESA(12) = 2.75
         PTESA(13) = 3.0
         PTESA(14) = 4.0
         PTESA(15) = 5.0
         PTESA(16) = 6.0
         PTESA(17) = 7.0
         PTESA(18) = 8.0
         PTESA(19) = 9.0 
         PTESA(20) = 10.0
         PTESA(21) = 20.0
         PTESA(22) = 30.0
         PTESA(23) = 40.0
         PTESA(24) = 50.0
         PTESA(25) = 60.0
         PTESA(26) = 70.0
         PTESA(27) = 80.0
         PTESA(28) = 90.0
         PTESA(29) = 100.0
         PTESA(30) = 150.0
c
         do in = 1,4 
c
            tmpne = 1.0e17 * 10**(in) 
c
c            tmpne = 1.0e18

C
C---- CALCULATE IONISATION RATES ...
C
C
c            DO IZ = 0, CION-1
            DO IZ = 0, min(printiz,cion-1)
               do ik = 1,ntemp
                  call mcrates(ptesa(ik),ptesa(ik),tmpne,
     >                         iz,cion,cion,rion,rrec,rcxr,0)
                  pcoef(ik,iz+1) = rion 
                end do
            end do
c
            write (6,*) 'B2-FRATES Ionization rate coefficients:'
     >                  //' Density = ', tmpne
c
            write(6,1001)  (iz,iz=0,min(printiz,cion-1))
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=0,min(printiz,cion-1))
            end do

C
c           RECOMBINATION
c
            DO IZ = 1, min(printiz+1,cion)
               do ik = 1,ntemp
                  call mcrates(ptesa(ik),ptesa(ik),tmpne,
     >                         iz,cion,cion,rion,rrec,rcxr,0)
                  pcoef(ik,iz+1) = rrec
                end do
            end do
c
            write (6,*) 'B2-FRATES Recombination rate coefficients:'
     >                  //' Density = ', tmpne
c
            write(6,1001)  (iz,iz=1,min(printiz+1,cion))
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=1,min(printiz+1,cion))
            end do

C
c           Charge Exchange
c
            DO IZ = 1, min(printiz+1,cion)
               do ik = 1,ntemp
                  call mcrates(ptesa(ik),ptesa(ik),tmpne,
     >                         iz,cion,cion,rion,rrec,rcxr,0)
                  pcoef(ik,iz+1) = rcxr
                end do
            end do
c
            write (6,*) 'B2-FRATES CX Recombination rate coefficients:'
     >                  //' Density = ', tmpne
c
            write(6,1001)  (iz,iz=1,min(printiz+1,cion))
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=1,min(printiz+1,cion))
            end do

c
c        End density loop
c
         end do
c
c
c     INEL atomic physics data
c
      elseif (cdatopt.eq.3) then

         ntemp = 30

         PTESA(1) = 0.1
         PTESA(2) = 0.25
         PTESA(3) = 0.50
         PTESA(4) = 0.75
         PTESA(5) = 1.0
         PTESA(6) = 1.25
         PTESA(7) = 1.50
         PTESA(8) = 1.75
         PTESA(9) = 2.0
         PTESA(10) = 2.25
         PTESA(11) = 2.5
         PTESA(12) = 2.75
         PTESA(13) = 3.0
         PTESA(14) = 4.0
         PTESA(15) = 5.0
         PTESA(16) = 6.0
         PTESA(17) = 7.0
         PTESA(18) = 8.0
         PTESA(19) = 9.0 
         PTESA(20) = 10.0
         PTESA(21) = 20.0
         PTESA(22) = 30.0
         PTESA(23) = 40.0
         PTESA(24) = 50.0
         PTESA(25) = 60.0
         PTESA(26) = 70.0
         PTESA(27) = 80.0
         PTESA(28) = 90.0
         PTESA(29) = 100.0
         PTESA(30) = 150.0
c
C
C---- CALCULATE IONISATION RATES ...
C
C
            DO IZ = 0, CION-1
               do ik = 1,ntemp
                  call imprates(ptesa(ik),iz,cion,rion,rrec,rcxr)
                  pcoef(ik,iz+1) = rion 
                end do
            end do
c
            write (6,*) 'INEL Ionization rate coefficients:'
c
            write(6,1001)  (iz,iz=0,cion-1)
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=0,cion-1)
            end do
C
c           RECOMBINATION
c
            DO IZ = 1, CION
               do ik = 1,ntemp
                  call imprates(ptesa(ik),iz,cion,rion,rrec,rcxr)
                  pcoef(ik,iz+1) = rrec
                end do
            end do
c
            write (6,*) 'INEL Recombination rate coefficients:'
c
            write(6,1001)  (iz,iz=1,cion)
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=1,cion)
            end do

C
c           Charge Exchange
c
            DO IZ = 1, CION
               do ik = 1,ntemp
                  call imprates(ptesa(ik),iz,cion,rion,rrec,rcxr)
                  pcoef(ik,iz+1) = rcxr
                end do
            end do
c
            write (6,*) 'INEL CX Recombination rate coefficients:'
c
            write(6,1001)  (iz,iz=1,cion)
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz+1),iz=1,cion)
            end do

c
c     End of CDATOPT BLOCK - ADAS/ADPAK/INEL - DATA PRINT OUT
c
      endif
c
      write (6,*) 'End of RATE DATA'
      write (6,*)
c
c     End of CPRINT if statement       
c
      endif
c
c     Formatting
c

1001  format(3x,'N',5x,'Tev',2x,10(2x,'IZ =',i2,2x)) 
1000  format(i5,1x,f7.2,10(1x,1p,g9.3))

c
      RETURN
      END
