      SUBROUTINE IZTAU (CRMI,NXS,NYS,CION,CIZB,CIOPTA)                          
      use mod_params
      use mod_comt2
      use mod_cnoco
      use mod_cadas
      use mod_slcom
      IMPLICIT  NONE
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c      INCLUDE   'cnoco'                                                         
C     INCLUDE   (CNOCO)                                                         
c slmod begin
c      INCLUDE   'slcom'                                                         
c slmod end
c
c      include   'cadas'
c
      INTEGER   NXS,NYS,CION,CIZB,CIOPTA                                        
      REAL      CRMI                                                            
C                                                                               
C***********************************************************************        
C                                                                               
C      THIS VERSION BY C.M.FARRELL   SEPT 1987                                  
C      ---------------------------------------                                  
C         THIS SUBROUTINE CALCULATES EXPECTED IONISATION & RECOMBINATIONION     
C      TIMES TO THE IONISATION STATES 0 TO CION                                 
C      OF A REQUESTED ELEMENT.                                                  
C      THE IONISATION TIMES ARE CALCULATED FOR EACH OF A GIVEN                  
C      SET OF TEMPERATURES AND DENSITIES.                                       
C      THEY ARE CALCULATED FROM THE EQUATIONS                                   
C      T = 1/(NE.EXP(-KIS/TB).(TB/KIS)**(1/2).SUM)                              
C      SUM = SUM(ANS.(LOG10(TB/KIS))**N), N=0,5                                 
C      NE IS THE BACKGROUND ELECTRON DENSITY SO NE = RNB * CIZB                 
C      KIS IS A FUNCTION OF ATOMIC NUMBER AND IONISATION STATE AND              
C      ANS IS A FUNCTION OF ATOMIC NUMBER, IONISATION STATE AND N.              
C      THE EQUATION ABOVE IS VALID FOR BACKGROUND TEMPERATURES IN 
C      THE RANGE  I/10 < TB < 10*I. WHERE I IS THE IONIZTAION ENERGY
C      OF THE SPECIFIC STATE OF THE PARTICULAR ION. FOR TB > 10*I THE 
C      FOLLOWING EQUATION IS USED.
C      T = 1/((TB/KIS)**(-1/2).{BNS(1).LN(TB/KIS) + SUM1})
C      SUM1 = SUM(BNS(N+2) . (TB/KIS)**(-N))  N = 0,2 
C      VALUES OF BOTH FOR A SELECTED SET OF ELEMENTS                            
C      ARE STORED IN THE ROUTINE.                                               
C      THE RECOMBINATION TIMES ARE FOR ELECTRON-ION RECOMBINATION               
C      AND ARE EXTRACTED FROM THE ABELS VAN MAANEN PACKAGE.                     
C                                                                               
C         PARAMETERS -                                                          
C     CFIZS : ARRAY EXPECTED IONISATION TIMES TO BE RETURNED IN (S)             
C     CFRCS : ARRAY EXPECTED RECOMBINATIONN TIMES TO BE RETURNED IN (S)         
C     CRMI   : IMPURITY MASS                                                    
C     NQXSO  : NUMBER OF SETS OF IONISATION TIMES OUTBOARD REQUIRED             
C     NQXSI  : NUMBER OF SETS OF IONISATION TIMES INBOARD REQUIRED              
C     CION   : ELEMENT (BY ATOMIC NUMBER)                                       
C     CTEMBS : SET OF BACKGROUND TEMPERATURES TB  (EV)                          
C     CRNBS  : SET OF BACKGROUND ION DENSITIES NB (M**-3)                       
C     CIZB   : BACKGROUND ION IONISATION STATE ZB                               
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
      INTEGER   N,IN,IZ,IX,IY,M,ik,ir
      PARAMETER (N=6,M=4)                                                       
      REAL      KIS(12,10),ANS(N,12,10),BNS(M,12,10),SZ(12,10)               
      REAL      KTI,LOGKTI,SUM,TEMP,RIZB,DENOM                                  
C                                                                               
      ! jdemod - adding 3D poloidal zone support
      integer :: pz
c

      INTEGER   ICODE 
      INTEGER   IONH,IONHE,IONLI,IONBE,IONB,IONC,IONO,IONCR,IONFE,IONNI         
c slmod begin
     +          ,IONN,NSTEP,TAG
      REAL      V,E,SIGMA,SIGMAV,BETA,X,X1,X2,XSTEP,WEIGHT

      PARAMETER(IONN=11)
c slmod end
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
C     THE FOLLOWING GIVE THE COEFFICIENT FOR THE HIGH BACKGROUND   
C     TEMPERATURE CASES
C
C                                                                               
      DATA  (BNS(IN,1,IONH), IN = 1, M)/                                        
     H   2.4617E-08, 9.5987E-08, -9.2464E-07, 3.9974E-06/             
C                                                                               
      DATA  (BNS(IN,1,IONHE), IN = 1, M)/                                       
     H   3.1373E-08, 4.7894E-08, -7.7361E-07, 3.7367E-06/
      DATA  (BNS(IN,2,IONHE), IN = 1, M)/                                       
     H   3.0755E-09, 1.1892E-08, -1.1505E-07, 5.0451E-07/
C                                                                               
      DATA  (BNS(IN,1,IONLI), IN = 1, M)/                                       
     L   4.5456E-08, 2.7800E-07, -1.5830E-06, 5.4652E-06/
      DATA  (BNS(IN,2,IONLI), IN = 1, M)/                                       
     L   7.3446E-09, 5.3440E-10, -5.6346E-08, 2.9555E-07/
      DATA  (BNS(IN,3,IONLI), IN = 1, M)/                                       
     L   1.9755E-09, -1.0918E-09, 8.8579E-10, 6.0762E-09/
C                                                                               
      DATA  (BNS(IN,1,IONBE), IN = 1, M)/ 
     B   2.1732E-07, -2.1648E-07, 2.8110E-07, 5.3070E-07/
      DATA  (BNS(IN,2,IONBE), IN = 1, M)/ 
     B   6.4891E-08, -1.1198E-07, 6.4469E-07, -1.9617E-06/
      DATA  (BNS(IN,3,IONBE), IN = 1, M)/ 
     B   2.7903E-09, -2.4350E-10, -9.8088E-09, 5.2155E-08/
      DATA  (BNS(IN,4,IONBE), IN = 1, M)/ 
     B   8.3329E-10, -4.6056E-10, 3.7364E-10, 2.5630E-09/
C                                                                               
      DATA  (BNS(IN,1,IONB), IN = 1, M)/ 
     B   3.0952E-07, -4.9239E-07, 1.3750E-06, -2.5382E-06/
      DATA  (BNS(IN,2,IONB), IN = 1, M)/ 
     B   4.8123E-08, -4.1483E-08, 5.5408E-08, 1.0022E-07/
      DATA  (BNS(IN,3,IONB), IN = 1, M)/ 
     B   1.1270E-08, -9.8068E-09, 9.3591E-09, 4.9053E-08/
      DATA  (BNS(IN,4,IONB), IN = 1, M)/ 
     B   1.2752E-09, -1.1129E-10, -4.4828E-09, 2.3836E-08/
      DATA  (BNS(IN,5,IONB), IN = 1, M)/ 
     B   4.2656E-10, -2.3576E-10, 1.9126E-10, 1.3120E-09/
C                                                                               
      DATA  (BNS(IN,1,IONC), IN = 1, M)/ 
     C   3.7442E-07, -6.5826E-07, 2.05201E-06, -4.4688E-06/
      DATA  (BNS(IN,2,IONC), IN = 1, M)/ 
     C   6.0150E-08, -4.0215E-08, -2.7928E-08, 5.5510E-07/
      DATA  (BNS(IN,3,IONC), IN = 1,M)/ 
     C   1.7372E-08, -1.5019E-08, 4.1269E-08, -8.8550E-08/
      DATA  (BNS(IN,4,IONC), IN = 1, M)/ 
     C   5.8147E-09, -5.5184E-09, -4.9521E-09, 8.9187E-08/
      DATA  (BNS(IN,5,IONC), IN = 1, M)/ 
     C   6.8613E-10, -5.9877E-11, -2.4120E-09, 1.2825E-08/
      DATA  (BNS(IN,6,IONC), IN = 1, M)/ 
     C   2.4680E-10, -1.3641E-10, 1.1066E-10, 7.5911E-10/
C                                                                               
      DATA  (BNS(IN,1,IONO), IN = 1, M)/ 
     O   3.2684E-07, -6.7409E-07, 2.3910E-06, -6.0366E-06/
      DATA  (BNS(IN,2,IONO), IN = 1, M)/ 
     O   4.9066E-08, 2.1777E-08, -5.4684E-07, 2.7249E-06/
      DATA  (BNS(IN,3,IONO), IN = 1, M)/ 
     O   1.7523E-08, 2.7812E-08, -3.1460E-07, 1.4290E-06/
      DATA  (BNS(IN,4,IONO), IN = 1, M)/ 
     O   6.6074E-09, 1.6482E-08, -1.7855E-07, 8.0360E-07/
      DATA  (BNS(IN,5,IONO), IN = 1, M)/ 
     O   2.0855E-09, 7.7560E-09, -4.8817E-08, 1.7611E-07/
      DATA  (BNS(IN,6,IONO), IN = 1, M)/ 
     O   9.3640E-10, 4.2792E-10, 1.3924E-09, -5.9247E-09/
      DATA  (BNS(IN,7,IONO), IN = 1, M)/ 
     O   2.6498E-10, -2.3124E-11, -9.3150E-10, 4.9530E-09/
      DATA  (BNS(IN,8,IONO), IN = 1, M)/ 
     O   1.0406E-10, -5.7515E-11, 4.6660E-11, 3.2007E-10/
C                                                                               
      DATA  (BNS(IN,1,IONCR), IN = 1, M)/
     C   3.9118E-07, -3.3201E-07, 4.6427E-07, 7.3151E-07/
      DATA  (BNS(IN,2,IONCR), IN = 1, M)/ 
     C   3.5012E-07, -1.7316E-07, -5.9974E-08, 2.0628E-06/
      DATA  (BNS(IN,3,IONCR), IN = 1, M)/ 
     C   1.6474E-07, -1.4478E-07, 2.2817E-07, 2.7072E-07/
      DATA  (BNS(IN,4,IONCR), IN = 1, M)/ 
     C   5.6173E-08, -3.8160E-08, 2.9052E-08, 1.5319E-07/
      DATA  (BNS(IN,5,IONCR), IN = 1, M)/ 
     C   3.3365E-08, -2.4577E-08, -9.9264E-08, 8.2255E-07/
      DATA  (BNS(IN,6,IONCR), IN = 1, M)/ 
     C   1.8396E-08, -8.6349E-09, -4.1371E-08, 3.2309E-07/
      DATA  (BNS(IN,7,IONCR), IN = 1, M)/ 
     C   9.3502E-09, -4.1167E-09, 1.4319E-09, 3.2742E-08/
      DATA  (BNS(IN,8,IONCR), IN = 1, M)/ 
     C   6.5503E-09, -2.6522E-09, 6.1972E-10, 2.2761E-08/
      DATA  (BNS(IN,9,IONCR), IN = 1, M)/ 
     C   4.5173E-09, -1.7327E-09, 5.6334E-10, 1.4330E-08/
      DATA  (BNS(IN,10,IONCR), IN = 1, M)/
     C   2.5340E-09, -2.3010E-10, -1.1865E-09, 7.8059E-09/
      DATA  (BNS(IN,11,IONCR), IN = 1, M)/ 
     C   2.0060E-09, -5.9170E-10, -8.4996E-11, 6.1521E-09/
      DATA  (BNS(IN,12,IONCR), IN = 1, M)/ 
     C   1.6949E-09, -7.9867E-10, 6.0553E-10, 4.7274E-09/
C                                                                               
      DATA  (BNS(IN,1,IONFE), IN = 1, M)/ 
     F   3.4615E-07, -4.3391E-07, 1.7618E-06, -5.3387E-06/
      DATA  (BNS(IN,2,IONFE), IN = 1, M)/ 
     F   2.1839E-07, -3.8290E-07, 1.0867E-06, -1.8547E-06/
      DATA  (BNS(IN,3,IONFE), IN = 1, M)/ 
     F   1.7440E-07, -3.0837E-07, 9.8490E-07, -2.1688E-06/
      DATA  (BNS(IN,4,IONFE), IN = 1, M)/
     F   5.7171E-08, -2.8314E-08, -8.8711E-09, 3.3389E-07/
      DATA  (BNS(IN,5,IONFE), IN = 1, M)/ 
     F   4.3518E-07, -3.8256E-08, 6.0497E-08, 7.0196E-08/
      DATA  (BNS(IN,6,IONFE), IN = 1, M)/ 
     F   2.9874E-08, -1.3733E-08, -3.2756E-08, 3.1443E-07/
      DATA  (BNS(IN,7,IONFE), IN = 1, M)/  
     F   2.2082E-08, -1.3775E-08, -2.1429E-09, 1.5872E-07/
      DATA  (BNS(IN,8,IONFE), IN = 1, M)/  
     F   1.2868E-08, -6.3972E-09, -1.7289E-08, 1.6339E-07/
      DATA  (BNS(IN,9,IONFE), IN = 1, M)/ 
     F   5.3531E-09, -2.3566E-09, 8.1922E-10, 1.8746E-08/
      DATA  (BNS(IN,10,IONFE), IN = 1, M)/  
     F   3.0753E-09, 2.1304E-09, -4.5673E-09, 3.2858E-09/
      DATA  (BNS(IN,11,IONFE), IN = 1, M)/  
     F   2.7519E-09, -1.0533E-09, 3.3880E-10, 8.7309E-09/
      DATA  (BNS(IN,12,IONFE), IN = 1, M)/  
     F   1.6092E-09, -1.4612E-10, -7.5350E-10, 4.9571E-09/
C                                                                               
      DATA  (BNS(IN,1,IONNI), IN = 1, M)/ 
     N   2.6366E-07, 2.0322E-07, 2.7985E-07, 2.6778E-07/
      DATA  (BNS(IN,2,IONNI), IN = 1, M)/  
     N   2.3090E-07, -3.8813E-07, 9.9775E-07, -1.2563E-06/
      DATA  (BNS(IN,3,IONNI), IN = 1, M)/ 
     N   1.0627E-07, -1.4265E-07, 3.6229E-07, -5.1833E-07/
      DATA  (BNS(IN,4,IONNI), IN = 1, M)/ 
     N   6.9061E-08, -4.9755E-08, -1.0610E-07, 1.0766E-06/
      DATA  (BNS(IN,5,IONNI), IN = 1, M)/ 
     N   4.4555E-08, -7.8783E-08, 2.5162E-07, -5.5408E-07/
      DATA  (BNS(IN,6,IONNI), IN = 1, M)/ 
     N   2.0100E-08, -9.9498E-09, -3.0294E-09, 1.1756E-07/
      DATA  (BNS(IN,7,IONNI), IN = 1, M)/  
     N   2.0347E-08, -1.7900E-08, 2.8098E-08, 3.1384E-08/
      DATA  (BNS(IN,8,IONNI), IN = 1, M)/ 
     N   9.3730E-09, -6.3674E-09, 4.8476E-09, 2.5561E-08/
      DATA  (BNS(IN,9,IONNI), IN = 1, M)/
     N   7.1788E-09, -5.2880E-09, -2.1358E-08, 1.7698E-07/
      DATA  (BNS(IN,10,IONNI), IN = 1, M)/ 
     N   4.7098E-09, -2.2108E-09, -1.0592E-08, 8.2721E-08/
      DATA  (BNS(IN,11,IONNI), IN = 1, M)/
     N   3.3232E-09, -1.4630E-09, 5.0857E-10, 1.1638E-08/
      DATA  (BNS(IN,12,IONNI), IN = 1, M)/ 
     N   2.4897E-09, -1.0081E-09, 2.3555E-10, 8.6511E-09/
C
C
C
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
C                     IDENTIFY REQUESTED ELEMENT AND SET ICODE                  
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
c slmod begin
      ELSE IF (CION.EQ.7) THEN
        IF (N2OPT.EQ.1) 
     +    WRITE(0,*) 'Warning (IZTAU): Nitrogen data has been patched. '
        ICODE = IONN
c slmod end
      ELSE IF (CION .EQ. 8) THEN                                                
         ICODE = IONO                                                           
      ELSE IF (CION .EQ. 10) THEN 
         WRITE(6,*) '1:', CION,CIOPTA
         IF (CIOPTA.GT.2) GOTO 550
C
C        IF NEON IS SPECIFIED BUT THE IONIZATION OPTION IS NOT
C        FROM THE NOCORONA PACKAGE THEN AN ERROR HAS OCCURRED 
C        ASSUME CARBON (AS BELOW) AND CONTINUE
C
         WRITE(6, *)  ' ERROR - SUBROUTINE IZTAU '                              
         WRITE(6, *)  ' IONISATION TIMES REQUESTED'                             
         WRITE(6, *)  ' OF ELEMENT ATOMIC NUMBER ', CION                        
         WRITE(6, *)  ' THESE ARE NOT AVAILABLE, CARBON ASSUMED'                
         CION = 6                                                               
         ICODE = IONC                                                           
      ELSE IF (CION.EQ.17) THEN 
         IF (CIOPTA.GT.2) GOTO 550
C
C        IF CHLORINE IS SPECIFIED BUT THE IONIZATION OPTION IS NOT
C        FROM THE NOCORONA PACKAGE THEN AN ERROR HAS OCCURRED 
C        ASSUME CARBON (AS BELOW) AND CONTINUE
C
         WRITE(6, *)  ' ERROR - SUBROUTINE IZTAU '                              
         WRITE(6, *)  ' IONISATION TIMES REQUESTED'                             
         WRITE(6, *)  ' OF ELEMENT ATOMIC NUMBER ', CION                        
         WRITE(6, *)  ' THESE ARE NOT AVAILABLE, CARBON ASSUMED'                
         CION = 6                                                               
         ICODE = IONC                                                           
      ELSE IF (CION.EQ.18) THEN 
         IF (CIOPTA.GT.2) GOTO 550
C
C        IF ARGON IS SPECIFIED BUT THE IONIZATION OPTION IS NOT
C        FROM THE NOCORONA PACKAGE THEN AN ERROR HAS OCCURRED 
C        ASSUME CARBON (AS BELOW) AND CONTINUE
C
         WRITE(6, *)  ' ERROR - SUBROUTINE IZTAU '                              
         WRITE(6, *)  ' IONISATION TIMES REQUESTED'                             
         WRITE(6, *)  ' OF ELEMENT ATOMIC NUMBER ', CION                        
         WRITE(6, *)  ' THESE ARE NOT AVAILABLE, CARBON ASSUMED'                
         CION = 6                                                               
         ICODE = IONC                                                           
      ELSE IF (CION .EQ. 24) THEN                                               
         ICODE = IONCR                                                          
      ELSE IF (CION .EQ. 26) THEN                                               
         ICODE = IONFE                                                          
      ELSE IF (CION .EQ. 28) THEN                                               
         ICODE = IONNI                                                          
      else if (cion.eq.74) then 
         ! tungsten specified
         IF (CIOPTA.GT.2) GOTO 550
         WRITE(6, *)  ' ERROR - SUBROUTINE IZTAU '                              
         WRITE(6, *)  ' IONISATION TIMES REQUESTED'                             
         WRITE(6, *)  ' OF ELEMENT ATOMIC NUMBER ', CION                        
         WRITE(6, *)  ' THESE ARE NOT AVAILABLE, CARBON ASSUMED'                
         CION = 6                                                               
         ICODE = IONC                                                           
      ELSE                                                                      
         WRITE(6, *)  ' ERROR - SUBROUTINE IZTAU '                              
         WRITE(6, *)  ' IONISATION TIMES REQUESTED'                             
         WRITE(6, *)  ' OF ELEMENT ATOMIC NUMBER ', CION                        
         WRITE(6, *)  ' THESE ARE NOT AVAILABLE, CARBON ASSUMED'                
         CION = 6                                                               
         ICODE = IONC                                                           
      ENDIF                                                                     
C
C                                                                               
550   CONTINUE 
      RIZB = REAL (CIZB)                                                        
C                                                                               
C-----------------------------------------------------------------------        
C     USE MAXIMUM VALUES FOR IONISATION RATES IF CIOPTA=2.                      
C-----------------------------------------------------------------------        
C                                                                               
       IF (CIOPTA.EQ.2) GOTO 300                                                
C                                                                               
C-----------------------------------------------------------------------        
C     USE NOCORONA PACKAGE TO CALCULATE RATES IF CIOPTA=3,4,5,6                 
C-----------------------------------------------------------------------        
C                                                                               
       IF (CIOPTA.EQ.3.OR.CIOPTA.EQ.4.OR.CIOPTA.EQ.5.OR.CIOPTA.EQ.6)            
     >   GOTO 600                                                               
C                                                                               
C-----------------------------------------------------------------------        
C  STANDARD CASE:     FOR EACH IONISATION LEVEL                                 
C-----------------------------------------------------------------------        
C                                                                               
       DO pz = 1,maxpzone
       DO 200  IZ = 0, CION-1                                                    
       DO 100 IY = -NYS, NYS                                                    
        DO 100 IX = 1, NXS                                                    
          KTI = CTEMBS(IX,IY,pz) / KIS(IZ+1,ICODE)                                 
C
C         FOR A TEMPERATURE RATIO LESS THAN 10 USE THE STANDARD
C         FORMULA OTHERWISE USE THE HIGH TEMPERATURE FORMULA  
C
        IF (KTI.LE.10.0) THEN 
C                                                                               
C-------- IONISATION TIME VERY LARGE ...                                        
C                                                                               
          IF (KTI .LT. 0.02)  THEN                                              
               CFIZS(IX,IY,IZ,pz) = 1.E20                                          
C                                                                               
C-------- SUM SERIES IN RIZ EQUATION                                            
C-------- (LOOP IS EXPANDED OUT FOR OPTIMIZATION)                               
C                                                                               
          ELSE                                                                  
            LOGKTI = LOG10(KTI)                                                 
            SUM  = ANS(1,IZ+1,ICODE)                                            
            TEMP = LOGKTI                                                       
            SUM  = SUM  + ANS(2,IZ+1,ICODE) * TEMP                              
            TEMP = TEMP * LOGKTI                                                
            SUM  = SUM  + ANS(3,IZ+1,ICODE) * TEMP                              
            TEMP = TEMP * LOGKTI                                                
            SUM  = SUM  + ANS(4,IZ+1,ICODE) * TEMP                              
            TEMP = TEMP * LOGKTI                                                
            SUM  = SUM  + ANS(5,IZ+1,ICODE) * TEMP                              
            TEMP = TEMP * LOGKTI                                                
            SUM  = SUM  + ANS(6,IZ+1,ICODE) * TEMP                              
C                                                                               
C---------- CALCULATE EXPECTED IONISATION TIME                                  
C---------- THE 1.E6 FACTOR IS BECAUSE ANS IS IN CM**3/S                        
C                                                                               
            DENOM = SQRT(KTI) * SUM * CRNBS(IX,IY,pz) * RIZB                       
            IF (DENOM.NE.0.0) THEN                                              
              CFIZS(IX,IY,IZ,pz) = 1.E6 * EXP(1.0/KTI) / DENOM                     
            ELSE                                                                
              CFIZS(IX,IY,IZ,pz) = 1.E20                                           
            ENDIF                                                               
          ENDIF                                                                 
C
C         HIGH RELATIVE BACKGROUND TEMPERATURE
C           
        ELSE
C                                                                               
C-------- SUM SERIES IN RIZ EQUATION                                            
C-------- (LOOP IS EXPANDED OUT FOR OPTIMIZATION)                               
C                                                                               
            SUM  = BNS(1,IZ+1,ICODE)*LOG(KTI)                          
            SUM  = SUM  + BNS(2,IZ+1,ICODE) 
            TEMP = 1.0/KTI                                                   
            SUM  = SUM  + BNS(3,IZ+1,ICODE) * TEMP                              
            TEMP = TEMP /KTI                                                
            SUM  = SUM  + BNS(4,IZ+1,ICODE) * TEMP                              
C                                                                               
C---------- CALCULATE EXPECTED IONISATION TIME                                  
C---------- THE 1.E6 FACTOR IS BECAUSE ANS IS IN CM**3/S                        
C                                                                               
            DENOM = (1.0/SQRT(KTI)) * SUM * CRNBS(IX,IY,pz) * RIZB                
            IF (DENOM.NE.0.0) THEN                                              
              CFIZS(IX,IY,IZ,pz) = 1.E6 / DENOM                     
            ELSE                                                                
              CFIZS(IX,IY,IZ,pz) = 1.E20                                           
            ENDIF                                                               
         ENDIF
         CFRCS(IX,IY,IZ,pz) = 0.0                                                 
  100   CONTINUE                                                                
  200 CONTINUE                                                                  

! END OF PZ LOOP
      
      end do
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
       do pz=1,maxpzone
       DO 410 IZ = 0, CION-1                                                     
       DO 400 IY = -NYS, NYS                                                    
        DO 400 IX = 1, NXS                                                      
          CFIZS(IX,IY,IZ,pz) = 1.0E6 /                                             
     >      (CRNBS(IX,IY,pz) * RIZB * SZ(IZ+1,ICODE))                              
          CFRCS(IX,IY,IZ,pz) = 0.0                                                 
  400   CONTINUE                                                                
  410 CONTINUE                                                                  
      end do
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
      cfizs = 0.0
      cfrcs = 0.0
      !CALL RZERO(CFIZS, MAXNXS*MAXNYS*(MAXIZS+1))
      !CALL RZERO(CFRCS, MAXNXS*MAXNYS*(MAXIZS+1))
c
c     NOCORONA package specified.
c
      if (cdatopt.eq.0) then

c
         do pz = 1,maxpzone
         DO 700 IY = -NYS, NYS                                                     
C                                                                               
C---- THE INITIALISATION OF NOCORONA IS DONE IN RUNLIM3 ...                     
C---- CALCULATE ELECTRON DENSITIES IN CM-3 FROM ION DENSITIES IN M-3.           
C                                                                               
      DO 610 IX = 1, NXS                                                        
        PNES(IX) = CRNBS(IX,IY,pz) * 1.E-6 * RIZB                                  
        ptes(ix) = ctembs(ix,iy,pz)
 610  CONTINUE                                                                  
C                                                                               
C---- CALCULATE IONISATION AND RECOMBINATION RATES ...                          
C---- VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY GOES           
C---- OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.                
C                                                                               
      CALL RRATES (PTES, PNES, PRATES, MAXNXS)                          
C                                                                               
      DO 630 IZ = 0, CION                                                       
        DO 620 IX = 1, NXS                                                      
          IF (IZ.EQ.CION) THEN                                                  
            CFIZS(IX,IY,IZ,PZ) = 0.0                                               
          ELSE                                                                  
            CFIZS(IX,IY,IZ,pz) = 1.E6 /                                            
     >        (PRATES(1,IZ+1,1,IX)*CRNBS(IX,IY,pz)*RIZB)                           
          ENDIF                                                                 
          IF ((CIOPTA.EQ.3.OR.CIOPTA.EQ.5).AND.IZ.GT.0) THEN                    
            CFRCS(IX,IY,IZ,pz) = 1.E6 /                                            
     >        (PRATES(2,IZ,1,IX)*CRNBS(IX,IY,pz)*RIZB)                             
          ELSE                                                                  
            CFRCS(IX,IY,IZ,pz) = 0.0                                               
          ENDIF                                                                 
  620   CONTINUE                                                                
  630 CONTINUE                                                                  
C                                                                               
  700 CONTINUE                                                                  
      ! end of pz loop
      end do

c     slmod begin
      IF (N2OPT.EQ.1) CALL GetN2Rate(rizb)
c slmod end
C                                                                               
 9100 FORMAT(//1X,'ELEMENT WITH MASS ',G11.4,' IS NOT INCLUDED',                
     >    ' IN THE NOCORONA PACKAGE.',/)                                        

c
c     ADAS package specified
c
c
c     *** NOTE ***   IY = IR   ... IX = IK
c
c     Code has been patched from the DIVIMP (4.04) implementation
c

      elseif (cdatopt.eq.1) then
c
C
C---- DO ONE RING AT A TIME
C
         do pz = 1,maxpzone
         DO 800 IR = -NYS, NYS
C
      DO 710 IK = 1, NXS
        PNESA(IK) = CRNBS(IK,IR,pz) * RIZB
        PTESA(IK) = CTEMBS(IK,IR,pz)
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
      DO 730 IZ = 0, CION-1
        CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NXS,PTESA,PNESA,   
     >              PCOEF(1,IZ+1))                                      
c        CALL ADASRD(YEAR,YEARDF,CION,IZ+1,ICLASS,NXS,PTESA,PNESA,
c     >              PCOEF(1,IZ+1))
        DO 720 IK = 1, NXS

          IF (PCOEF(IK,IZ+1).GT.0.0) THEN
            CFIZS(IK,IR,IZ,pz) = 1.0 / (PCOEF(IK,IZ+1)*PNESA(IK))
          ELSE
            CFIZS(IK,IR,IZ,pz) = HI
          ENDIF
  720   CONTINUE
  730 CONTINUE
C
      IF (CIOPTA.EQ.3 .OR. CIOPTA.EQ.5) THEN
        ICLASS = 1
C
        DO 735 IZ = 1, CION
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NXS,PTESA,PNESA,   
     >                PCOEF(1,IZ))                                      
c          CALL ADASRD(YEAR,YEARDF,CION,IZ,ICLASS,NXS,PTESA,PNESA,
c     >                PCOEF(1,IZ))
          DO 725 IK = 1, NXS
            IF (PCOEF(IK,IZ).GT.0.0) THEN
              CFRCS(IK,IR,IZ,pz) = 1.0 / (PCOEF(IK,IZ)*PNESA(IK))
            ELSE
              CFRCS(IK,IR,IZ,pz) = HI
            ENDIF
  725     CONTINUE
  735   CONTINUE
      ENDIF
C
  800 CONTINUE
      ! end of pz loop
      end do
c     
c     End of ADAS/Nocorona block
c
      endif

C                                                                               
C-----------------------------------------------------------------------        
C         FINISHED: VALIDATE IONISATION DATA BEFORE LEAVING                     
C         PRINT TABLE OF RESULTS IF REQUIRED                                    
C-----------------------------------------------------------------------        
C                                                                               
  999 CONTINUE                                                                  
      do pz = 1,maxpzone
      DO  IZ = 0, CION                                                       
       WRITE (6,9000) IZ                                                        
       DO IY = -NYS, NYS                                                    
         DO  IX = 1, NXS                                                   
          IF (CFIZS(IX,IY,IZ,pz).LT.0.0.OR.CFIZS(IX,IY,IZ,pz).GT.1.E20)             
     >        CFIZS(IX,IY,IZ,pz) = 1.E20                                           
          IF (CFRCS(IX,IY,IZ,pz).LT.0.0.OR.CFRCS(IX,IY,IZ,pz).GT.1.E20)             
     >        CFRCS(IX,IY,IZ,pz) = 1.E20                                    
          IF (5*(IX/5).EQ.IX.AND.IY.EQ.INT(NYS/2))          
     >      WRITE (6,9001) CTEMBS(IX,IY,pz),CRNBS(IX,IY,pz),
     >                     CFIZS(IX,IY,IZ,pz),CFRCS(IX,IY,IZ,pz)        
         end do
       end do
      end do
      end do
C     
 9000 FORMAT(//1X,'IZTAU: SAMPLE IONISATION AND RECOMBINATION TIMES',           
     >  ' FOR IONISATION STATE',I4,                                             
     > //1X,'TEMPERATURE       DENSITY     ',
     >      'IONISATION TIME   RECOMB. TIME', 
     >  /1X,'------------------------------',
     >      '------------------------------') 
 9001 FORMAT(1X,F10.3,5X,G11.4,5X,G11.4,7X,G11.4)                           
      RETURN                                                                    
      END                                                                       
C
C
C
      SUBROUTINE ADASRD(YEAR,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)     
      use mod_params
      use error_handling
      use mod_cadas2
C     
C  READ THE REQUESTED RATE COEFFICIENT FROM THE ADAS MASTER ELEMENT
C  FILES:
C        ICLASS = 1: RECOMBINATION RATE COEFFICIENT
C                 2: IONISATION RATE COEFFICIENT
C                 3: CHARGE EXCHANGE RECOMBINATION COEFFICIENT
C                 4: POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLU
C                 5: POWER COEFFICIENT FOR LINE RADIATION
C                 6: POWER COEFFICIENT FOR CHARGE EXCHANGE
C  THIS ROUTINE USES THE STANDARD ADAS EXTRACTION ROUTINE D2DATA AND
C  REALLY ONLY PROVIDES A 'CLEAN' INTERFACE, TAKING CARE OF CHANGES
C  IN UNITS AND IN PRECISION OF VARIABLES.  IF THE REQUESTED DATA
C  DOESN'T EXIST (IFAIL=1 RETURNED FROM D2DATA) THE PROGRAM IS STOPPED.
C
      IMPLICIT NONE
c      INCLUDE   'params'
c
c      include    'params'
C     INCLUDE   "CADAS2"
c
c      include    'cadas2'
C
      CHARACTER*2 YEAR
      INTEGER IZ0, IZ1, ICLASS, NPTS
      REAL TE(NPTS), NE(NPTS), COEF(NPTS)
c
c      logical lintrp(maxpts)
C
      INTEGER I, J
C
      IEVCUT = 0
C
      DO J = 1, NPTS
        DTEV(J) = DBLE(ALOG10(TE(J)))
        DDENS(J) = DBLE(ALOG10(NE(J)*1.0E-6))
      ENDDO
      CALL D2DATA(YEAR, YEAR, TITLF, IFAIL,
     >            IZ0, IZ1, ICLASS, NPTS, IEVCUT,
     >            MAXADS, ITMAXD, IDMAXD, IZMAXD,
     >            DTEV, DDENS,
     >            DTEVD, DDENSD, DRCOFD, ZDATA,
     >            DRCOFI
c     >            , LINTRP
     >            )
      IF (IFAIL.EQ.1) THEN
        !WRITE(6,1000) IZ0, IZ1,YEAR
        !WRITE(7,1000) IZ0, IZ1,YEAR
        !WRITE(0,1000) IZ0, IZ1,YEAR

        write(error_message_data,1000) IZ0, IZ1,YEAR

        call errmsg('ADASRD:',trim(error_message_data))

        STOP
      ENDIF
C
C  EXTRAPOLATED VALUES ARE RETURNED AS ZERO!
C
      DO J = 1, NPTS
        IF (DRCOFI(J).NE.0.0) THEN
          COEF(J) = 10.**SNGL(DRCOFI(J)) * 1.0E-6
        ELSE
          COEF(J) = 0.0
        ENDIF
      ENDDO
C
 1000 FORMAT(' ERROR READING REQUESTED ATOMIC DATA!',
     >       ' MASTER ELEMENT FILE FOR NUCLEAR CHARGE ',I2,
     >       ' AND ION CHARGE ',I2,
     >       ' WAS NOT FOUND IN YEAR ',A2)
C
      RETURN
      END
 

c
c slmod begin
c
c ======================================================================
c
c subroutine: GetN2Rate
c
c ======================================================================
c
      SUBROUTINE GetN2Rate(rizb)
      use mod_params
      use mod_comt2
      use mod_cnoco
      use mod_comxyt
      use mod_slcom
      implicit none

      REAL    rizb

c      INCLUDE   'params'                                                        
c      INCLUDE   'comt2'                                                         
c      INCLUDE   'comxyt'                                                         
c      INCLUDE   'cnoco'                                                         
c      INCLUDE   'slcom'                                                         

      INTEGER   IONN,NSTEP,TAG
      REAL      V,E,SIGMA,SIGMAV,BETA,X,X1,X2,XSTEP,WEIGHT,SigmaI,SigmaN
      REAL      RES,OLDSIG
      integer   ix,iy,iz
      integer pz


      ! fix pz = 1 - this feature will not be compatible with 3D poloidal plasma slices
      ! development can be done later if this is desirable

      pz = 1
      
      WRITE(63,*) ' '
      WRITE(63,*) 'Molecular nitrogen ionisation data:'
      WRITE(63,*) ' '

      IY = 1
      IZ = 0

      DO IX = 1, NXS
c
c Integrate over the Maxwellian:
c
        RES     = 0.0
        OLDSIG = 0.0

        SIGMAV = 0.0
        
        BETA = 9.1095E-31 / (2.0 * 1.6E-19 * CTEMBS(IX,IY,pz))
       
        NSTEP = 10001

        X1 = SQRT(2.0 *    0.1 * 1.6E-19 / 9.1095E-31) 
        X2 = SQRT(2.0 * 1000.0 * 1.6E-19 / 9.1095E-31) 
 
        XSTEP = (X2 - X1) / NSTEP 
c
c       Integrate using Simpson's method (neglecting end points):
c
        TAG = 1
        DO X = X1+XSTEP, X2-XSTEP, XSTEP
       
          IF (TAG.EQ.1) THEN
            WEIGHT = 4.0 / 3.0
            TAG    = 0
          ELSE
            WEIGHT = 2.0 / 3.0
            TAG    = 1
          ENDIF
         
          V      = X
          E      = 0.5 * 9.1095E-31 * V**2.0 / 1.6E-19
          SIGMA  = SigmaI(E)

          OLDSIG = SIGMAV
          SIGMAV = SIGMAV + XSTEP * V**3.0 * SIGMA *
     +                      EXP(-BETA * V**2.0) * WEIGHT
          RES    = 100.0 * ABS(OLDSIG - SIGMAV) / SIGMAV

          IF (RES<0.001) GOTO 10
        ENDDO

        WRITE(0,*) 'WARNING: SIGMAV for NIRATE not converged',RES

10      CONTINUE

        SIGMAV = SIGMAV * 4.0 * PI * (BETA / PI)**1.5

        NIRATE(IX) = 1.0 / (SIGMAV * CRNBS(IX,IY,pz) * RIZB)
      ENDDO
c
c Find the rate co-efficient for dissociation:
c
      DO IX = 1, NXS
        RES    = 0.0
        OLDSIG = 0.0

        SIGMAV = 0.0
        
        BETA = 9.1095E-31 / (2.0 * 1.6E-19 * CTEMBS(IX,IY,pz))
       
        NSTEP = 10001

        X1 = SQRT(2.0 *    0.1 * 1.6E-19 / 9.1095E-31) 
        X2 = SQRT(2.0 * 1000.0 * 1.6E-19 / 9.1095E-31) 

        XSTEP = (X2 - X1) / NSTEP 

        TAG = 1
        DO X = X1+XSTEP, X2-XSTEP, XSTEP
       
          IF (TAG.EQ.1) THEN
            WEIGHT = 4.0 / 3.0
            TAG    = 0
          ELSE
            WEIGHT = 2.0 / 3.0
            TAG    = 1
          ENDIF
         
          V      = X
          E      = 0.5 * 9.1095E-31 * V**2.0 / 1.6E-19
          SIGMA  = SigmaN(E)

          OLDSIG = SIGMAV
          SIGMAV = SIGMAV + XSTEP * V**3.0 * SIGMA *
     +                      EXP(-BETA * V**2.0) * WEIGHT
          RES    = 100.0 * ABS(OLDSIG - SIGMAV) / SIGMAV

c          IF (IX.EQ.28) THEN
c            WRITE(63,'(2E12.4,F10.1,E12.4,E16.6)') 
c     +        X,V,E,SIGMA,SIGMAV 
c          ENDIF

          IF (RES<0.001) GOTO 20
        ENDDO

        WRITE(0,*) 'WARNING: SIGMAV for NNRATE not converged',RES

20      CONTINUE

        SIGMAV = SIGMAV * 4.0 * PI * (BETA / PI)**1.5

        NNRATE(IX) = 1.0 / (SIGMAV * CRNBS(IX,IY,pz) * RIZB)

        N2RATE(IX) = (NIRATE(IX) * NNRATE(IX)) / 
     +               (NIRATE(IX) + NNRATE(IX))
      ENDDO
c
c Output:
c
      IY = 1
      IZ = 0

      WRITE(63,*) ' '
      WRITE(63,'(2A4,A8,A6,8A10)') 
     +  'ix','iz','xs','te','nb',
     +  'nirate','nnrate','n2rate',
     +  'cfizs0',
     +  's-v I','s-v N','s-v N2'
  
      DO IX = 1, NXS
        WRITE(63,'(2I4,F8.3,F6.1,8E10.3)') 
     +    IX,IZ,xs(ix),
     +    ctembs(ix,iy,pz),crnbs(ix,iy,pz),
     +    nirate(ix),nnrate(ix),n2rate(ix),
     +    CFIZS(IX,IY,IZ,pz),
     +    1.0/(NiRATE(IX) * (CRNBS(IX,IY,pz) * RIZB)),
     +    1.0/(NnRATE(IX) * (CRNBS(IX,IY,pz) * RIZB)),
     +    1.0/(N2RATE(IX) * (CRNBS(IX,IY,pz) * RIZB))
      ENDDO

      WRITE(63,*) ' '

      RETURN
      END
c
c END OF SUBROUTINE ====================================================
c
c
c ======================================================================
c
c function: SigmaN
c
c ======================================================================
c      
      REAL FUNCTION SigmaI(E)

      REAL E

      SigmaI  = 5.580 * LOG(E / 15.5) 
      SigmaI  = SigmaI + (-5.300) * (1.0 - 15.5 / E) ** 1.0
      SigmaI  = SigmaI + (-2.836) * (1.0 - 15.5 / E) ** 2.0
      SigmaI  = SigmaI / (15.5 * E) * 1.0E-17 

      IF (SigmaI.LT.0.0) SigmaI = 0.0

      RETURN 
      END
c
c ======================================================================
c
c
c ======================================================================
c
c function: SigmaN
c
c ======================================================================
c      
      REAL FUNCTION SigmaN(E)

      REAL E

      IF (E.LT.12.0) THEN
        SigmaN = 0.0

      ELSEIF (E.LT.14.0) THEN
        SigmaN = (0.04 - 0.01) / (14.0 - 12.0) * (E - 12.0) + 0.01

      ELSEIF (E.LT.16.0) THEN
        SigmaN = (0.20 - 0.04) / (16.0 - 14.0) * (E - 14.0) + 0.04

      ELSEIF (E.LT.18.0) THEN
        SigmaN = (0.36 - 0.20) / (18.0 - 16.0) * (E - 16.0) + 0.20

      ELSEIF (E.LT.20.0) THEN
        SigmaN = (0.52 - 0.36) / (20.0 - 18.0) * (E - 18.0) + 0.36

      ELSEIF (E.LT.25.0) THEN
        SigmaN = (0.87 - 0.52) / (25.0 - 20.0) * (E - 20.0) + 0.52

      ELSEIF (E.LT.30.0) THEN
        SigmaN = (1.04 - 0.52) / (30.0 - 25.0) * (E - 25.0) + 0.87

      ELSEIF (E.LT.50.0) THEN
        SigmaN = (1.23 - 1.04) / (50.0 - 30.0) * (E - 30.0) + 1.04

      ELSE
        SigmaN = (0.95 - 1.23) / (200.0 - 50.0) * (E - 50.0) + 1.23
      ENDIF

      SigmaN = SigmaN * 1.0E-16 * 1.0E-4

      IF (SigmaN.LT.0.0) SigmaN = 0.0

      RETURN 
      END
c
c ======================================================================
c
c slmod end
c
