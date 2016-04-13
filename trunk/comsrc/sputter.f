      SUBROUTINE SYIELD (MATT,MATP,CNEUTD,
     >                   CBOMBF,CBOMBZ,cbomb_frac,CION,CIZB,CRMB,cebd)
      IMPLICIT NONE
      INTEGER MATT,MATP,CNEUTD,CBOMBF,CBOMBZ,CION,CIZB
      REAL    CRMB,cebd,cbomb_frac
C
C  *********************************************************************
C  *                                                                   *
C  *  SYIELD:  SETS UP SPUTTERING YIELD DATA IN COMMON AREA CYIELD.    *
C  *  THE DATA IS TAKEN FROM BODANSKY  NUCLEAR FUSION  1984            *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  MATP : CODE INTEGER GIVING BACKGROUND MATERIAL                   *
C  *  CNEUTD : NEUT SPUTTER OPTION                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include    'params'
C     INCLUDE "CYIELD"
      include    'cyield'
      REAL ETH(7,12), ETF(7,12), Q(7,12) , ebd(12)
      LOGICAL IDATA(7,12)
      INTEGER I,J,NSPEC
      CHARACTER*18 TARMAT(19)
      CHARACTER*6  PLAMAT(7)
C
      NSPEC=7
      NTARS=12
C
C  NSPEC = NUMBER OF IMPURITY SPECIES IN PLASMA.
C  NTARS = NUMBER OF TARGET MATERIALS.
C   CETH = THRESHOLD ENERGY FOR TARGET-ION INTERACTION CONSIDERED (EV)
C   CETF = THOMAS-FERMI INTERACTION POTENTIAL (EV)
C     CQ = YIELD FACTOR (ATOMS/ION)
C CIDATA = LOGICAL FLAG INDICATING WHETHER DATA IS AVAILABLE (T OR F)
C   CEBD = TARGET MATERIAL SURFACE BINDING ENERGY (EV)
C
C  SEE NOTES 33, 86, 113, 301, 303 FOR SPUTTERING CONSTANTS.
C
 
      DATA TARMAT/
     &  ' ALUMINIUM       ',' BERYLLIUM       ',' COPPER          ',
     &  ' GRAPHITE        ',' TITANIUM        ',' IRON            ',
     &  ' NICKEL          ',' MOLYBDENUM      ',' TUNGSTEN        ',
     &  ' BORON           ',' TITANIUM CARBIDE',' SILICON CARBIDE ',
     &  ' "DEUTERIUM"     ',' "HELIUM"        ',' "NEON"          ',
     &  ' "ARGON"         ',' "OXYGEN"        ',' "CHLORINE"      ',
     &  ' "NITROGEN"      ' /
 
      DATA PLAMAT/
     &  ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    '/
 
      DATA  ETH/
     &             53.0, 24.0, 22.0, 20.0,  0.0, 25.0,  0.0,
c
c     jdemod - july 2009
c
c     Beryllium data.The Eth, Etf and Q values here have been 
c     adjusted to the values used in 1990 LIM runs which matched
c     the JET Be limiter experimental results.
c     Numbers have been adjusted for D->Be and Be->Be
c
c     &            20.0,  9.0, 21.0, 30.0, 40.0, 25.0, 70.0,
     &             20.0, 10.0, 21.0, 30.0, 40.0, 35.0, 70.0,
     &             55.0, 34.0, 30.0, 19.0,  0.0, 36.0,  0.0,
     &             35.0, 30.0, 30.0, 29.0, 42.0, 42.0,  0.0,
CN    &             35.0, 30.0, 30.0, 29.0, 42.0, 42.0,  0.0,
CO    &             35.0, 28.0, 30.0, 29.0, 44.0, 44.0,  0.0,
     &             80.0, 50.0, 40.0, 25.0, 35.0, 40.0,  0.0,
     &             64.0, 44.0, 40.0, 33.0, 35.0, 40.0,  0.0,
     &             50.0, 31.0, 25.0, 22.0, 28.0, 34.4, 48.0,
     &            199.0, 90.0, 70.0, 46.0, 55.0, 64.0,  0.0,
     &            443.0,220.0,140.0,110.0, 80.0, 65.0, 40.0,
     &             20.0, 20.0, 21.0, 30.0, 40.0, 35.0,  0.0,
     &             65.0, 31.0, 29.0, 18.0,  0.0,  0.0, 90.0,
     &             38.0, 25.0, 22.0, 18.0,  0.0,  0.0, 90.0/
 
      DATA ETF/
     &            1059.0,1097.0,1135.0,2448.0,10295.0,34545.0,15718.0,
c
c     jdemod - july 2009
c
c     Beryllium data.The Eth, Etf and Q values here have been 
c     adjusted to the values used in 1990 LIM runs which matched
c     the JET Be limiter experimental results.
c     Numbers have been adjusted for D->Be and Be->Be
c
c     &           256.0, 282.0, 308.0, 780.0, 4152.0, 2208.0, 6970.0,
     &            256.0, 280.0, 308.0, 780.0, 4152.0, 2210.0, 6970.0,
     &            2926.0,2971.0,3017.0,6292.0,22696.0,244619.0,32723.0,
     &            415.0,447.0,479.0,1087.0,5687.0,5687.0,9298.0,
CN    &            415.0,447.0,479.0,1087.0,5687.0,5687.0,9298.0,
CO    &            415.0,446.0,479.0,1087.0,5680.0,5680.0,9298.0,
     &            2054.0,2096.0,2138.0,4502.0,16947.0,117898.0,24843.0,
     &            2544.0,2589.0,2634.0,5514.0,20247.0,174096.0,29839.0,
     &            2799.0,2846.0,2893.0,6044.0,22011.0,206960.0,31856.0,
     &            4718.0,4767.0,4816.0,9944.0,34183.0,533048.0,48322.0,
     &          9870.0,9923.0,9977.0,20373.0,66507.0,1998599.0,91979.0,
     &             333.0, 361.0, 389.0, 895.0,4857.0,3716.0,8023.0,
     &            1161.0,1198.0,1236.0,2653.0,10917.0,0.0,16540.0,
     &            767.0,803.0,840.0,1841.0,8311.0,0.0,12997.0/
 
      DATA Q/
     &            0.043 ,0.093 ,0.2   ,0.34  ,0.0   ,5.4   ,0.0  ,
c
c     jdemod - july 2009
c
c     Beryllium data.The Eth, Etf and Q values here have been 
c     adjusted to the values used in 1990 LIM runs which matched
c     the JET Be limiter experimental results.
c     Numbers have been adjusted for D->Be and Be->Be
c
c     &           0.1   ,0.3   ,0.24  ,0.59  ,1.6   ,1.4   ,1.3  ,
     &            0.1   ,0.27  ,0.24  ,0.59  ,1.6   ,0.94  ,1.3  ,
     &            0.1   ,0.23  ,0.35  ,0.81  ,0.0   ,17.0  ,0.0  ,
     &            0.035 ,0.1   ,0.20  ,0.32  ,1.5   ,1.5   ,0.0  ,
CN    &            0.035 ,0.1   ,0.20  ,0.32  ,1.5   ,1.5   ,0.0  ,
CO    &            0.035 ,0.14  ,0.20  ,0.32  ,1.9   ,1.9   ,0.0  ,
     &            0.017 ,0.055 ,0.1   ,0.125 ,1.0   ,3.7   ,0.0  ,
     &            0.042 ,0.13  ,0.21  ,0.44  ,3.2   ,13.0  ,0.0  ,
     &            0.041 ,0.115 ,0.22  ,0.45  ,2.5   ,11.7  ,1.55 ,
     &            0.007 ,0.023 ,0.045 ,0.12  ,0.93  ,18.0  ,0.0  ,
     &            0.007 ,0.019 ,0.038 ,0.106 ,0.93  ,20.0  ,2.2  ,
     &            0.1   ,0.15  ,0.24  ,0.59  ,1.6   ,0.94  ,0.0  ,
     &            0.025 ,0.04  ,0.096 ,0.36  ,0.0   ,0.0   ,1.3  ,
     &            0.039 ,0.11  ,0.18  ,0.36  ,0.0   ,0.0   ,1.3  /
 
      DATA IDATA/
     &            4*.TRUE.,.FALSE.,.TRUE.,.FALSE.,
     &            7*.TRUE.,
     &            4*.TRUE.,.FALSE.,.TRUE.,.FALSE.,
     &            6*.TRUE.,.FALSE.,
     &            6*.TRUE.,.FALSE.,
     &            6*.TRUE.,.FALSE.,
     &            7*.TRUE.,
     &            6*.TRUE.,.FALSE.,
     &            7*.TRUE.,
     &            6*.TRUE.,.FALSE.,
     &            4*.TRUE.,2*.FALSE.,.TRUE.,
     &            4*.TRUE.,2*.FALSE.,.TRUE./
C
C  TABLE OF BINDING ENERGIES TO BE USED AS DEFAULT WHEN ZERO IS
C  SPECIFIED IN THE INPUT FILE.  FOR THE TWO COMPOUNDS AND FOR
C  THE GASEOUS IMPURITIES, I HAVE SET EBD = 0
C
      DATA EBD/3.36,3.38,3.52,7.42,4.89,4.34,
     &         4.46,6.83,8.68,5.73,0.00,0.00/
C
C-----------------------------------------------------------------------
C INITIALISE COMMON BLOCK CYIELD.
C-----------------------------------------------------------------------
C
      DO 10 I=1,NTARS
        DO 20 J=1,NSPEC
         CETH(J,I)   = ETH(J,I)
         CETF(J,I)   = ETF(J,I)
         CQ(J,I)     = Q(J,I)
         CIDATA(J,I) = IDATA(J,I)
   20   CONTINUE
   10 CONTINUE
C
C-----------------------------------------------------------------------
C  ASSIGN TARGET AND BACKGROUND ION MATERIALS.
C  AT PRESENT TARGET MATERIALS 11,12 ARE NOT USED.
C  MATERIAL 13 IS A SPECIAL "DEUTERIUM" CASE, WHERE YIELD IS ALWAYS 1.0
C-----------------------------------------------------------------------
C
      MATT = 4
      IF (CION.EQ.13) MATT = 1
      IF (CION.EQ.4)  MATT = 2
      IF (CION.EQ.29) MATT = 3
      IF (CION.EQ.6)  MATT = 4
      IF (CION.EQ.22) MATT = 5
      IF (CION.EQ.26) MATT = 6
      IF (CION.EQ.28) MATT = 7
      IF (CION.EQ.42) MATT = 8
      IF (CION.EQ.74) MATT = 9
      IF (CION.EQ.5)  MATT = 10
      IF (CION.EQ.1)  MATT = 13
      IF (CION.EQ.2)  MATT = 14
      IF (CION.EQ.10) MATT = 15
      IF (CION.EQ.18) MATT = 16
      IF (CION.EQ.8)  MATT = 17
      IF (CION.EQ.17) MATT = 18
      if (cion.eq.7)  MATT = 19
C
C---- TARGET BINDING ENERGY FROM INTERNAL DATA IF EXTERNAL INPUT ZERO
C
      IF (CEBD.EQ.0.0 .AND. MATT.LE.12) CEBD = EBD(MATT)
C
C---- PLASMA MATERIAL.  CAN BE SET EXPLICITLY FROM INPUT DATA
C
      MATP = 6
      IF (NINT(CRMB).LE.4) MATP = NINT (CRMB)
      IF (CIZB.EQ.6) MATP = 5
      IF (CIZB.EQ.8) MATP = 7

      IF (CNEUTD.EQ.1) then
         MATP = CBOMBF
         flux_frac = cbomb_frac
      else
         flux_frac = 1.0
      endif

c
      call prb
      call prc('TARGET MATERIAL IS     '//TARMAT(MATT))
      call prc('BOMBARDING IONS ARE    '//PLAMAT(MATP))
c
      IF (CNEUTD.EQ.1) then
         CALL PRI ('         WITH ZIMP', CBOMBZ)
         CALL PRR ('         ZIMP FLUX FRACTION',flux_frac)
      endif

      RETURN
      END
C
C
C     function to implement customized yields for W; Krieger IPP/97
C
      real function yldtung(tempe, tempi)

      integer i
      real tempd(8), yield(8), tempe, tempi

*     data from Thoma, PPCF 97
      data tempd /2.5, 5., 10., 15., 20., 30., 50., 70./
      data yield /1.e-7, 2.8e-5, 2.3e-4, 4.5e-4, 7.5e-4, 1.e-3,
     >            1.3e-3, 2.e-3/

      if (tempi.lt.tempd(1)) then
c
c       cannot set it to 0 -> code crashes if no sputtered atoms
c       generated   Krieger IPP/97
c        yldtung = 1.e-7
c
        yldtung = 0.0  
      else if (tempi.ge.tempd(6)) then
        yldtung = yield(6)
      else
        do 10 i=1,5
          if ((tempi.ge.tempd(i)).and.(tempi.lt.tempd(i+1))) then
            yldtung=yield(i)+(yield(i+1)-yield(i))/
     >                       (tempd(i+1)-tempd(i))*
     >                       (tempi     -tempd(i))
            return
          endif
  10    continue
      endif

      return
      end
C
C
c      New: pass Te, Ti to yield function; Krieger IPP/97
c
       FUNCTION YIELD(MATP,MATT,ENERGY,Te,ti)
       use eckstein_2002_yield_data
       use eckstein_2007_yield_data
       IMPLICIT none
       REAL YIELD,ENERGY,X1,X12,X2,te,ti
       INTEGER MATP,MATT
C
C  *********************************************************************
C  *                                                                   *
C  *  YIELD:   CALCULATES YIELD OF MATERIAL "MATP" HITTING TARGET      *
C  *  MADE OF MATERIAL "MATT" WITH AN ENERGY "ENERGY".                 *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATP : CODE INTEGER GIVING MATERIAL OF IMPACTING PARTICLE        *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  ENERGY : ENERGY OF IMPACTING PARTICLE                            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "CYIELD"
      include    'cyield'
c
      include 'params'
      include 'comtor'     
c
      real yld93,yld96,yldtung,yld_be_2002,yld_c_2002
      external yld93,yld96,yldtung,yld_be_2002,yld_c_2002
C
      if (csputopt.eq.2) then 
         yield = yld93(MATP,MATT,ENERGY) * flux_frac
         return
      elseif (csputopt.eq.3) then 
         yield = yld96(MATP,MATT,ENERGY) * flux_frac
         return
      elseif (csputopt.eq.4) then 
         yield = const_yield * flux_frac
         return
      elseif (csputopt.eq.5) then
c
c        Krieger IPP/97 - custom W yield 
c
c        special sputter formula for tungsten
c
c        Carbon sputtering? 
c
c         if (matp.eq.5.and.matt.eq.9.and.csputopt.eq.4) then
c
         if (matt.eq.9.and.matp.ne.6.and.ti.gt.0.0) then
           yield = yldtung(Te,Ti) * flux_frac
c
c
c
         elseif ((matt.eq.2.or.matt.eq.4).and.
     >           (matp.eq.2.or.matp.eq.6)) then 
c
c           Add special option for Beryllium or Carbon sputtering here 
c       
c           Eckstein 2002
c
c           The second parameter can only have one value so the incident angle  
c           data for hydrogenic sputtering must match the 
c           self sputtering.
c
c           Be -> Be
c           D -> Be  
c
c           The second parameter is an option to allow selection 
c           of different data sets as required. 
c     
c           It specifies the incident angle : a value of -1 is used to load
c           averaged data. 
c
            yield = yld2002(energy,matp,matt,extra_sputter_angle) 
     >              * flux_frac
c
c
c
         else
           yield = yld96(MATP,MATT,ENERGY) * flux_frac   
         endif
         return
c
      elseif (csputopt.eq.6) then 
c
c        This option calls yield routines based on Eckstein's
c        2007 tabulation/parameterization of the sputtering yield data
c
c        This data has no explicit angular dependence the numbers are
c        for normal incidence. 
c
c        The data in this routine is taken from 
c
! Topics in Applied Physics 110
!
! Behrisch and Eckstein (eds)
!
! "Sputtering by Particle Bombardment, Experiments and Computer Calculations from Threshold to MeV Energies"
! 
! Chapter "Sputtering Yields" by W. Eckstein, Springer 2007, pgs 33 - 186
c
c        At the present time only certain bombarding and target materials are supported - if an
c        unsupported combination is specified then the data defaults to '96 
c
         if (eckstein2007_data_available) then 

            yield = yield_2007(matp,matt,energy) * flux_frac

         else

            yield = yld96(MATP,MATT,ENERGY)  * flux_frac 

         endif
c
         return

c        
      endif
c
      IF (MATT.EQ.13.OR.MATT.EQ.14.OR.MATT.EQ.15
     >      .OR.MATT.EQ.16.OR.MATT.EQ.17.or.matt.eq.18
     >      .or.matt.eq.19) THEN
        YIELD = 1.0 * flux_frac
        RETURN
      ENDIF
C
      IF(ENERGY.LE.0.0) GOTO 100
      IF(CIDATA(MATP,MATT)) THEN
      IF((CETH(MATP,MATT)/ENERGY).GT.1.0E0) GOTO 100
         X1=ENERGY/CETF(MATP,MATT)
         X12=SQRT(X1)
         X2=CETH(MATP,MATT)/ENERGY
c
         YIELD=CQ(MATP,MATT)*(3.441*X12*LOG(X1+2.718))/(1.0+6.355*
     1          X12+X1*(6.882*X12-1.708))*(1-X2)*(1-X2)*(1.0-X2**
     2          (2.0/3.0)) * flux_frac
c
c     The 1984 Bohdansky paper has a typo - either equation 6.13
c     or equation 6.16 is wrong. If eq 6.16 was correct then the 
c     expression could be further simplified so I think 6.13 must
c     be correct and the above equation for the yield is thus 
c     correct. 
c
c     The following equation matches equation
c     6.16 from page 64 of the Nuclear Fusion 1984 article
c     by Langley, Bohdansky, Eckstein ... 
c
c     In the paper the expresion is (1-X2)**(2/3) not (1-X2**(2/3))
c     However, this disagrees with eq 6.13 in the same paper so there
c     is an inconsistency. 
c
c
c         YIELD=CQ(MATP,MATT)*(3.441*X12*LOG(X1+2.718))/(1.0+6.355*
c     1          X12+X1*(6.882*X12-1.708))*(1-X2)*(1-X2)*(1.0-X2)**
c     2          (2.0/3.0)
c
      ELSE
         YIELD=0.0
      ENDIF
      RETURN
C  ERROR TRAPPING, OUTSIDE RANGE SET YIELD=0,   E=0 OR E/E0 > 1
  100 YIELD=0.0
      END




C
C
C
      SUBROUTINE SYLD93(MATT,MATP,CNEUTD,
     >                  CBOMBF,CBOMBZ,cbomb_frac,CION,CIZB,CRMB,CEBD)
      IMPLICIT NONE
      INTEGER MATT,MATP,CNEUTD,CBOMBF,CBOMBZ,CION,CIZB
      REAL    CRMB,CEBD,cbomb_frac
C
C  *********************************************************************
C  *                                                                   *
C  *  SYLD93:  SETS UP SPUTTERING YIELD DATA IN COMMON AREA CYIELD.    *
C  *  THE DATA IS TAKEN FROM ECKSTEIN IPP 9/82 (FEB 1993)              *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  MATP : CODE INTEGER GIVING BACKGROUND MATERIAL                   *
C  *  CNEUTD : NEUT SPUTTER OPTION                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include    'params'
C     INCLUDE "CYIELD"
      include    'cyield'
      REAL ETH(7,12), ETF(7,12), Q(7,12), EBD(12)
      LOGICAL IDATA(7,12)
      INTEGER I,J,NSPEC
      CHARACTER*18 TARMAT(19)
      CHARACTER*6  PLAMAT(7)
C
      NSPEC=7
      NTARS=12
C
C  NSPEC = NUMBER OF IMPURITY SPECIES IN PLASMA.
C  NTARS = NUMBER OF TARGET MATERIALS.
C   CETH = THRESHOLD ENERGY FOR TARGET-ION INTERACTION CONSIDERED (EV)
C   CETF = THOMAS-FERMI INTERACTION POTENTIAL (EV)
C     CQ = YIELD FACTOR (ATOMS/ION)
C CIDATA = LOGICAL FLAG INDICATING WHETHER DATA IS AVAILABLE (T OR F)
C
C     DATA FROM ECKSTEIN'S LATEST REPORT IPP 9/82 HAS BEEN USED.
C     THE TWO COMPOUNDS (TITANIUM CARBIDE AND SILICON CARBIDE)
C     IN THE ORIGINAL REFERENCE ARE NOT AVAILABLE IN THE LATEST
C     REPORT AND HAVE BEEN REPLACED WITH LITHIUM AND CHROMIUM.
C     NOTE ALSO THAT ECKSTEIN HAS CHANGED THE DEFINITION OF HIS
C     NUCLEAR STOPPING CROSS SECTION.  HE HAS ALSO RECOMMENDED
C     THAT WE USE FITS TO EXPERIMENTAL DATA, WHERE AVAILABLE, IN
C     PREFERENCE TO THE CALCULATIONS BASED ON EMPIRICAL FORMULAE
C     (PP 335-338, WHICH ARE AVAILABLE FOR A LARGE RANGE OF
C     PROJECTILE-TARGET COMBINATIONS).  FOR THE TIME BEING I HAVE
C     REPLACED THE GENERAL TABLES WITH FITS TO EXPERIMENTAL DATA
C     ONLY FOR H, D, AND HE ON BE.  THIS RESULTS IN HIGHER YIELDS
C     (APPROX 2X).
C                       LORNE HORTON MAY 93
c
c     changed parameters for boron, D, He and B projectiles;
c     Krieger, IPP 95
c     the new values are from Eckstein IPP 9/82, exp. data p.28
C
 
      DATA TARMAT/
     &  ' ALUMINIUM       ',' BERYLLIUM       ',' COPPER          ',
     &  ' GRAPHITE        ',' TITANIUM        ',' IRON            ',
     &  ' NICKEL          ',' MOLYBDENUM      ',' TUNGSTEN        ',
     &  ' BORON           ',' LITHIUM         ',' CHROMIUM        ',
     &  ' "DEUTERIUM"     ',' "HELIUM"        ',' "NEON"          ',
     &  ' "ARGON"         ',' "OXYGEN"        ',' "CHLORINE"      ',
     &  ' "NITROGEN"      ' /
 
      DATA PLAMAT/
     &  ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    '/
 
      DATA  ETH/
     &     23.87, 14.86, 12.91, 12.51, 16.32, 24.02, 18.55,
C    &     12.99, 13.09, 14.69, 16.40, 28.08, 24.17, 32.71,
     &     12.2 , 10.0 , 14.69, 13.9 , 28.08, 24.17, 32.71,
     &     57.25, 28.90, 20.64, 17.07, 13.27, 25.17, 14.01,
     &     31.11, 27.64, 29.48, 32.15, 52.98, 52.98, 61.54,
     &     59.49, 31.51, 23.71, 20.56, 19.45, 34.96, 21.23,
     &     61.39, 31.63, 23.12, 19.54, 16.70, 31.03, 17.95,
     &     66.80, 34.12, 24.69, 20.67, 17.00, 31.89, 18.14,
     &    172.36, 83.30, 56.47, 44.28, 25.75, 48.83, 25.47,
     &    447.02,209.37,136.26,102.07, 41.20, 62.06, 35.92,
*          Krieger - changed to Boron experimental data. 
*    &     23.14, 21.56, 23.46, 25.83, 43.25, 40.97, 50.30,
     &     23.14, 13.00, 23.46, 20.00, 43.25, 65.00, 50.30,             
     &      6.22,  6.92,  8.03,  9.10, 15.94, 11.94, 18.61,
     &     54.47, 28.39, 21.01, 17.96, 16.07, 29.46, 17.40/
C
      DATA ETF/
     &    1059.,   1097.,   1135.,   2448.,  10297.,  34550.,  15720.,
     &     256.,    282.,    308.,    720.,   4153.,   2208.,   6971.,
     &    2926.,   2972.,   3017.,   6293.,  22701., 224652.,  32727.,
     &     415.,    447.,    479.,   1087.,   5688.,   5688.,   9298.,
     &    2054.,   2097.,   2139.,   4503.,  16949., 117915.,  24846.,
     &    2544.,   2590.,   2635.,   5517.,  20270., 174122.,  29437.,
     &    2799.,   2846.,   2893.,   6045.,  22014., 206991.,  31860.,
     &    4719.,   4768.,   4817.,   9945.,  34188., 533127.,  48329.,
     &    9871.,   9925.,   9978.,  20376.,  66517.,1998893.,  91993.,
     &     333.,    361.,    389.,    894.,   4856.,   3717.,   8021.,
     &     185.,    209.,    232.,    557.,   3506.,   1129.,   6014.,
     &    2296.,   2340.,   2383.,   5002.,  18577., 144458.,  27091./
C
      DATA Q/
     &     0.08,  0.14,  0.19,  0.37,  1.65,  4.21,  2.36,
C    &     0.07,  0.11,  0.14,  0.28,  1.00,  0.67,  1.35,
     &     0.128, 0.220, 0.14,  0.707, 1.00,  0.67,  1.35,
     &     0.08,  0.14,  0.19,  0.38,  1.83, 14.23,  2.73,
     &     0.05,  0.10,  0.10,  0.20,  0.75,  0.75,  1.02,
CD   &     0.05,  0.08,  0.10,  0.20,  0.75,  0.75,  1.02,
     &     0.06,  0.11,  0.15,  0.30,  1.41,  7.44,  2.07,
     &     0.07,  0.12,  0.16,  0.33,  1.59, 10.44,  2.36,
     &     0.07,  0.12,  0.16,  0.33,  1.60, 11.51,  2.38,
     &     0.05,  0.09,  0.12,  0.24,  1.20, 16.27,  1.81,
     &     0.04,  0.07,  0.10,  0.20,  1.02, 33.47,  1.55,
*          Krieger - changed to Boron experimental data. 
*    &     0.05,  0.08,  0.11,  0.21,  0.80,  0.67,  1.08,
     &     0.05,  .136,  0.11,  0.31,  0.80,  2.00,  1.08,              
     &     0.10,  0.16,  0.21,  0.40,  1.37,  0.69,  1.82,
     &     0.07,  0.12,  0.17,  0.34,  1.61,  9.54,  2.38/
C
      DATA IDATA/
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE./
C
C  TABLE OF BINDING ENERGIES TO BE USED AS DEFAULT WHEN ZERO IS
C  SPECIFIED IN THE INPUT FILE.  FOR THE GASEOUS IMPURITIES I HAVE
C  SET EBD = 0
C
      DATA EBD/3.36,3.38,3.52,7.42,4.89,4.34,
     &         4.46,6.83,8.68,5.73,1.67,4.12/
C
C-----------------------------------------------------------------------
C INITIALISE COMMON BLOCK CYIELD.
C-----------------------------------------------------------------------
C
      DO 10 I=1,NTARS
        DO 20 J=1,NSPEC
         CETH(J,I)   = ETH(J,I)
         CETF(J,I)   = ETF(J,I)
         CQ(J,I)     = Q(J,I)
         CIDATA(J,I) = IDATA(J,I)
   20   CONTINUE
   10 CONTINUE
C
C-----------------------------------------------------------------------
C  ASSIGN TARGET AND BACKGROUND ION MATERIALS.
C  MATERIALS 13-18 ARE SPECIAL "GAS" TARGET CASES, WHERE THE YIELD IS
C  ALWAYS 1.0
C-----------------------------------------------------------------------
C
      MATT = 4
      IF (CION.EQ.13) MATT = 1
      IF (CION.EQ.4)  MATT = 2
      IF (CION.EQ.29) MATT = 3
      IF (CION.EQ.6)  MATT = 4
      IF (CION.EQ.22) MATT = 5
      IF (CION.EQ.26) MATT = 6
      IF (CION.EQ.28) MATT = 7
      IF (CION.EQ.42) MATT = 8
      IF (CION.EQ.74) MATT = 9
      IF (CION.EQ.5)  MATT = 10
      IF (CION.EQ.3)  MATT = 11
      IF (CION.EQ.24) MATT = 12
      IF (CION.EQ.1)  MATT = 13
      IF (CION.EQ.2)  MATT = 14
      IF (CION.EQ.10) MATT = 15
      IF (CION.EQ.18) MATT = 16
      IF (CION.EQ.8)  MATT = 17
      IF (CION.EQ.17) MATT = 18
      IF (CION.EQ.7)  MATT = 19
C
C---- TARGET BINDING ENERGY FROM INTERNAL DATA IF EXTERNAL INPUT ZERO
C
      IF (CEBD.EQ.0.0 .AND. MATT.LE.12) CEBD = EBD(MATT)
C
C---- PLASMA MATERIAL.  CAN BE SET EXPLICITLY FROM INPUT DATA
C
      MATP = 6
      IF (NINT(CRMB).LE.4) MATP = NINT (CRMB)
      IF (CIZB.EQ.6) MATP = 5
      IF (CIZB.EQ.8) MATP = 7

      IF (CNEUTD.EQ.1) then
         MATP = CBOMBF
         flux_frac = cbomb_frac
      else
         flux_frac = 1.0
      endif

      WRITE (7,*) 'TARGET MATERIAL IS     ' , TARMAT(MATT)
      WRITE (7,*) 'BOMBARDING IONS ARE    ' , PLAMAT(MATP)

      IF (CNEUTD.EQ.1) then
         CALL PRI ('         WITH ZIMP', CBOMBZ)
         CALL PRR ('         ZIMP FLUX FRACTION',flux_frac)
      endif

      RETURN
      END
C
C
C
       REAL FUNCTION YLD93(MATP,MATT,ENERGY)
       IMPLICIT NONE
       REAL ENERGY,X1,X12,X2
       INTEGER MATP,MATT
C
C  *********************************************************************
C  *                                                                   *
C  *  YLD93:   CALCULATES YIELD OF MATERIAL "MATP" HITTING TARGET      *
C  *  MADE OF MATERIAL "MATT" WITH AN ENERGY "ENERGY".                 *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATP : CODE INTEGER GIVING MATERIAL OF IMPACTING PARTICLE        *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  ENERGY : ENERGY OF IMPACTING PARTICLE                            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "CYIELD"
      include    'cyield'
C
      IF (MATT.EQ.13.OR.MATT.EQ.14.OR.MATT.EQ.15
     >      .OR.MATT.EQ.16.OR.MATT.EQ.17.OR.MATT.EQ.18
     >      .or.matt.eq.19) THEN
        YLD93 = 1.0
        RETURN
      ENDIF
C
      IF(ENERGY.LE.0.0) GOTO 100
      IF(CIDATA(MATP,MATT)) THEN
      IF((CETH(MATP,MATT)/ENERGY).GT.1.0E0) GOTO 100
         X1=ENERGY/CETF(MATP,MATT)
         X12=SQRT(X1)
         X2=CETH(MATP,MATT)/ENERGY
         YLD93=CQ(MATP,MATT)*(0.5*LOG(1.0+1.2288*X1))
     1         /(X1+0.1728*X12+0.008*X1**0.1504)
     2         *(1-X2)*(1-X2)*(1.0-X2**(2.0/3.0))
      ELSE
         YLD93=0.0
      ENDIF
      RETURN
C  ERROR TRAPPING, OUTSIDE RANGE SET YIELD=0,   E=0 OR E/E0 > 1
  100 YLD93=0.0
      END
c
c
c
      SUBROUTINE SYLD96(MATT,MATP,CNEUTD,
     >                  CBOMBF,CBOMBZ,cbomb_frac,CION,CIZB,CRMB,CEBD)
      IMPLICIT none
      INTEGER MATT,MATP,CNEUTD,CBOMBF,CBOMBZ,CION,CIZB
      REAL    CRMB,CEBD,cbomb_frac
C
C  *********************************************************************
C  *                                                                   *
C  *  SYLD96:  SETS UP SPUTTERING YIELD DATA IN COMMON AREA CYIELD.    *
C  *  THE DATA IS TAKEN FROM ECKSTEIN IPP 9/82 (FEB 1993)              *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  MATP : CODE INTEGER GIVING BACKGROUND MATERIAL                   *
C  *  CNEUTD : NEUT SPUTTER OPTION                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include    'params'
C     INCLUDE "CYIELD"
      include    'cyield'
      REAL ETH(7,12), ETF(7,12), Q(7,12), EBD(12)
      LOGICAL IDATA(7,12)
      INTEGER I,J,NSPEC
      CHARACTER*18 TARMAT(19)
      CHARACTER*6  PLAMAT(7)
C
      NSPEC=7
      NTARS=12
C
C  NSPEC = NUMBER OF IMPURITY SPECIES IN PLASMA.
C  NTARS = NUMBER OF TARGET MATERIALS.
C   CETH = THRESHOLD ENERGY FOR TARGET-ION INTERACTION CONSIDERED (EV)
C   CETF = THOMAS-FERMI INTERACTION POTENTIAL (EV)
C     CQ = YIELD FACTOR (ATOMS/ION)
C CIDATA = LOGICAL FLAG INDICATING WHETHER DATA IS AVAILABLE (T OR F)
C
C     DATA FROM ECKSTEIN'S LATEST REPORT IPP 9/82 HAS BEEN USED.
C     THE TWO COMPOUNDS (TITANIUM CARBIDE AND SILICON CARBIDE)
C     IN THE ORIGINAL REFERENCE ARE NOT AVAILABLE IN THE LATEST
C     REPORT AND HAVE BEEN REPLACED WITH LITHIUM AND CHROMIUM.
C     NOTE ALSO THAT ECKSTEIN HAS CHANGED THE DEFINITION OF HIS
C     NUCLEAR STOPPING CROSS SECTION.  HE HAS ALSO RECOMMENDED
C     THAT WE USE FITS TO EXPERIMENTAL DATA, WHERE AVAILABLE, IN
C     PREFERENCE TO THE CALCULATIONS BASED ON EMPIRICAL FORMULAE
C     (PP 335-338, WHICH ARE AVAILABLE FOR A LARGE RANGE OF
C     PROJECTILE-TARGET COMBINATIONS).  FOR THE TIME BEING I HAVE
C     REPLACED THE GENERAL TABLES WITH FITS TO EXPERIMENTAL DATA
C     ONLY FOR H, D, AND HE ON BE.  THIS RESULTS IN HIGHER YIELDS
C     (APPROX 2X).
C                       LORNE HORTON MAY 93
C
 
      DATA TARMAT/
     &  ' ALUMINIUM       ',' BERYLLIUM       ',' COPPER          ',
     &  ' GRAPHITE        ',' TITANIUM        ',' IRON            ',
     &  ' NICKEL          ',' MOLYBDENUM      ',' TUNGSTEN        ',
     &  ' BORON           ',' LITHIUM         ',' CHROMIUM        ',
     &  ' "DEUTERIUM"     ',' "HELIUM"        ',' "NEON"          ',
     &  ' "ARGON"         ',' "OXYGEN"        ',' "CHLORINE"      ',
     &  ' "NITROGEN"      ' /
 
      DATA PLAMAT/
     &  ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    '/
c slmod begin
c      PLAMAT(8) is reserved for Ne in eckstein_2007_yield_data.f90. -SL, 10/05/12
c slmod end
      DATA  ETH/
     &     23.87, 14.86, 12.91, 12.51, 16.32, 24.02, 18.55,
c     &     12.99, 13.09, 14.69, 16.40, 28.08, 24.17, 32.71,
     &     12.2 , 10.0 , 14.69, 13.9 , 28.08, 24.17, 32.71,
     &     57.25, 28.90, 20.64, 17.07, 13.27, 25.17, 14.01,
c
c    NOTE: For H,D,T on C - set the Eth so that the threshold
c          energy is compatible with the value of Ebd when the 
c          truncated Thompson distribution is in use.
c
c          Eth = Ebd/(Gamma * (1.0 - Gamma)) 
c          Gamma H = 0.2840  Eth(min) = 36.50
c          Gamma D = 0.4898           = 29.70
c          Gamma T = 0.6400           = 35.62
c
c
c     &     31.00, 27.00, 29.00, 32.15, 52.98, 52.98, 61.54,
c
c
c     jdemod - july 2009
c              Change the data for He->C to match the Garching experimental
c              fits from Eckstein 2003 pg 46
c              Eth = 30.2, Etf=1087, Q=0.387
c
c     &     36.50, 29.70, 35.62, 32.15, 52.98, 52.98, 61.54,
c
c
     &     36.50, 29.70, 35.62, 30.2, 52.98, 52.98, 61.54,
c     &     31.11, 27.64, 29.48, 32.15, 52.98, 52.98, 61.54,
     &     59.49, 31.51, 23.71, 20.56, 19.45, 34.96, 21.23,
     &     61.39, 31.63, 23.12, 19.54, 16.70, 31.03, 17.95,
     &     66.80, 34.12, 24.69, 20.67, 17.00, 31.89, 18.14,
     &    172.36, 83.30, 56.47, 44.28, 25.75, 48.83, 25.47,
*     &    447.02,209.37,136.26,102.07, 41.20, 62.06, 35.92,
*         W->D corr. Krieger IPP/97
     &    447.02,201.00,136.26,102.07, 41.20, 62.06, 35.92,
     &     23.14, 21.56, 23.46, 25.83, 43.25, 40.97, 50.30,
     &      6.22,  6.92,  8.03,  9.10, 15.94, 11.94, 18.61,
     &     54.47, 28.39, 21.01, 17.96, 16.07, 29.46, 17.40/
C
      DATA ETF/
     &    1059.,   1097.,   1135.,   2448.,  10297.,  34550.,  15720.,
     &     256.,    282.,    308.,    720.,   4153.,   2208.,   6971.,
     &    2926.,   2972.,   3017.,   6293.,  22701., 224652.,  32727.,
     &     415.,    447.,    479.,   1087.,   5688.,   5688.,   9298.,
     &    2054.,   2097.,   2139.,   4503.,  16949., 117915.,  24846.,
     &    2544.,   2590.,   2635.,   5517.,  20270., 174122.,  29437.,
     &    2799.,   2846.,   2893.,   6045.,  22014., 206991.,  31860.,
     &    4719.,   4768.,   4817.,   9945.,  34188., 533127.,  48329.,
     &    9871.,   9925.,   9978.,  20376.,  66517.,1998893.,  91993.,
     &     333.,    361.,    389.,    894.,   4856.,   3717.,   8021.,
     &     185.,    209.,    232.,    557.,   3506.,   1129.,   6014.,
     &    2296.,   2340.,   2383.,   5002.,  18577., 144458.,  27091./
C
      DATA Q/
     &     0.08,  0.14,  0.19,  0.37,  1.65,  4.21,  2.36,
C    &     0.07,  0.11,  0.14,  0.28,  1.00,  0.67,  1.35,
     &     0.128, 0.220, 0.14,  0.707, 1.00,  0.67,  1.35,
     &     0.08,  0.14,  0.19,  0.38,  1.83, 14.23,  2.73,
c     jdemod - july 2009
c              Change the data for He->C to match the Garching experimental
c              fits from Eckstein 2003 pg 46
c              Eth = 30.2, Etf=1087, Q=0.387
c
c     &     0.035, 0.10,  0.12,  0.20,  0.75,  0.75,  1.02,
c
     &     0.035, 0.10,  0.12,  0.387,  0.75,  0.75,  1.02,
C93  &     0.05,  0.10,  0.10,  0.20,  0.75,  0.75,  1.02,
CD   &     0.05,  0.08,  0.10,  0.20,  0.75,  0.75,  1.02,
     &     0.06,  0.11,  0.15,  0.30,  1.41,  7.44,  2.07,
     &     0.07,  0.12,  0.16,  0.33,  1.59, 10.44,  2.36,
     &     0.07,  0.12,  0.16,  0.33,  1.60, 11.51,  2.38,
     &     0.05,  0.09,  0.12,  0.24,  1.20, 16.27,  1.81,
*     &     0.04,  0.07,  0.10,  0.20,  1.02, 33.47,  1.55,
*          W->D corr. Krieger IPP/97
     &     0.04,  0.035, 0.10,  0.20,  1.02, 33.47,  1.55,
     &     0.05,  0.08,  0.11,  0.21,  0.80,  0.67,  1.08,
     &     0.10,  0.16,  0.21,  0.40,  1.37,  0.69,  1.82,
     &     0.07,  0.12,  0.17,  0.34,  1.61,  9.54,  2.38/
C
      DATA IDATA/
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE.,
     &            7*.TRUE./
C
C  TABLE OF BINDING ENERGIES TO BE USED AS DEFAULT WHEN ZERO IS
C  SPECIFIED IN THE INPUT FILE.  FOR THE GASEOUS IMPURITIES I HAVE
C  SET EBD = 0
C
      DATA EBD/3.36,3.38,3.52,7.42,4.89,4.34,
     &         4.46,6.83,8.68,5.73,1.67,4.12/
C
C-----------------------------------------------------------------------
C INITIALISE COMMON BLOCK CYIELD.
C-----------------------------------------------------------------------
C
      DO 10 I=1,NTARS
        DO 20 J=1,NSPEC
         CETH(J,I)   = ETH(J,I)
         CETF(J,I)   = ETF(J,I)
         CQ(J,I)     = Q(J,I)
         CIDATA(J,I) = IDATA(J,I)
   20   CONTINUE
   10 CONTINUE
C
C-----------------------------------------------------------------------
C  ASSIGN TARGET AND BACKGROUND ION MATERIALS.
C  MATERIALS 13-18 ARE SPECIAL "GAS" TARGET CASES, WHERE THE YIELD IS
C  ALWAYS 1.0
C-----------------------------------------------------------------------
C
      MATT = 4
      IF (CION.EQ.13) MATT = 1
      IF (CION.EQ.4)  MATT = 2
      IF (CION.EQ.29) MATT = 3
      IF (CION.EQ.6)  MATT = 4
      IF (CION.EQ.22) MATT = 5
      IF (CION.EQ.26) MATT = 6
      IF (CION.EQ.28) MATT = 7
      IF (CION.EQ.42) MATT = 8
      IF (CION.EQ.74) MATT = 9
      IF (CION.EQ.5)  MATT = 10
      IF (CION.EQ.3)  MATT = 11
      IF (CION.EQ.24) MATT = 12
      IF (CION.EQ.1)  MATT = 13
      IF (CION.EQ.2)  MATT = 14
      IF (CION.EQ.10) MATT = 15
      IF (CION.EQ.18) MATT = 16
      IF (CION.EQ.8)  MATT = 17
      IF (CION.EQ.17) MATT = 18
      IF (CION.EQ.7)  MATT = 19
C
C---- TARGET BINDING ENERGY FROM INTERNAL DATA IF EXTERNAL INPUT ZERO
C
      IF (CEBD.EQ.0.0 .AND. MATT.LE.12) CEBD = EBD(MATT)
C
C---- PLASMA MATERIAL.  CAN BE SET EXPLICITLY FROM INPUT DATA
C
      MATP = 6
      IF (NINT(CRMB).LE.4) MATP = NINT (CRMB)
      IF (CIZB.EQ.6) MATP = 5
      IF (CIZB.EQ.8) MATP = 7

      IF (CNEUTD.EQ.1) then
         MATP = CBOMBF
         flux_frac = cbomb_frac
      else
         flux_frac = 1.0
      endif

      WRITE (7,*) 'TARGET MATERIAL IS     ' , TARMAT(MATT)
      WRITE (7,*) 'BOMBARDING IONS ARE    ' , PLAMAT(MATP)
      

      IF (CNEUTD.EQ.1) then
         CALL PRI ('         WITH ZIMP', CBOMBZ)
         CALL PRR ('         ZIMP FLUX FRACTION',flux_frac)
      endif

      RETURN
      END
C
C
C
       REAL FUNCTION YLD96(MATP,MATT,ENERGY)
       IMPLICIT none
       REAL ENERGY,X1,X12,X2
       INTEGER MATP,MATT
C
C  *********************************************************************
C  *                                                                   *
C  *  YLD96:   CALCULATES YIELD OF MATERIAL "MATP" HITTING TARGET      *
C  *  MADE OF MATERIAL "MATT" WITH AN ENERGY "ENERGY".                 *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATP : CODE INTEGER GIVING MATERIAL OF IMPACTING PARTICLE        *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  ENERGY : ENERGY OF IMPACTING PARTICLE                            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "CYIELD"
      include    'cyield'
C
      IF (MATT.EQ.13.OR.MATT.EQ.14.OR.MATT.EQ.15
     >      .OR.MATT.EQ.16.OR.MATT.EQ.17.OR.MATT.EQ.18
     >      .or.matt.eq.19) THEN
        YLD96 = 1.0
        RETURN
      ENDIF
C
      IF(ENERGY.LE.0.0) GOTO 100
      IF(CIDATA(MATP,MATT)) THEN
      IF((CETH(MATP,MATT)/ENERGY).GT.1.0E0) GOTO 100
         X1=ENERGY/CETF(MATP,MATT)
         X12=SQRT(X1)
         X2=CETH(MATP,MATT)/ENERGY
         YLD96=CQ(MATP,MATT)*(0.5*LOG(1.0+1.2288*X1))
     1         /(X1+0.1728*X12+0.008*X1**0.1504)
     2         *(1-X2)*(1-X2)*(1.0-X2**(2.0/3.0))
      ELSE
         YLD96=0.0
      ENDIF
      RETURN
C  ERROR TRAPPING, OUTSIDE RANGE SET YIELD=0,   E=0 OR E/E0 > 1
  100 YLD96=0.0
      END
c
c
c
      real function yldchem96 (e0,tmpdflux,matp,matt,tsurf)
      implicit none
      real e0,tmpdflux,tsurf
      integer matp,matt
c
      include 'params'
      include 'comtor'
c
c
c     YLDCHEM96: This function calculates the chemical sputtering
c              yield for Deuterium on Graphite based on the formulae
c              of  C. Garcia-Rosales, J. Roth. (1996)
c
c              The yield is composed of a temperature dependent
c              component and an impact energy dependent component.
c
c              David Elder, 1996 Apr 19
c
      real ysurf,ytherm,yphys,dflux,yphys2,eth2,csp3,c_const
      real d_const, etf, q, x1,x12,x2
      real tev,t
      real yield
      external yield
c
c     Test for negative dflux - convert to positive and issue 
c     warning.      
c
      if (tmpdflux.lt.0.0) then 
         dflux = abs(tmpdflux) 
         write(6,*) 'WARNING: Negative flux found in YLDCHEM2 :',
     >               tmpdflux  
      else
         dflux = tmpdflux
      endif 
c
      t = tsurf
c
      if (matt.ne.4) then
         write (6,*) 'Error: Chemical sputtering called'//
     >               ' for non-carbon target'
         yldchem96 = 0.0
         return
      endif
c
      if (matp.eq.1) then 
         eth2 = 2.0
         d_const = 250.0
         etf = 415.0
         q = 0.035
      elseif (matp.eq.2) then 
         eth2 = 1.0
         d_const = 125.0
         etf = 447.0
         q = 0.1  
      elseif (matp.eq.3) then 
         eth2 = 1.0
         d_const = 83.0
         etf = 479.0 
         q = 0.12
      endif
c
      tev = t  * ( 1.38e-23 / 1.6e-19 )
c
      yphys = yield(matp,matt,e0,0.0,0.0)
c
c     Calculate yphys2 - same as Yphys except eth is different 
c
      IF (E0.LE.0.0.or.(eth2/e0).GT.1.0E0) then 
        yphys2 = 0.0
      else
        X1=E0/etf
        X12=SQRT(X1)
        X2=eth2/e0
        Yphys2= q * (0.5*LOG(1.0+1.2288*X1))
     >        /(X1+0.1728*X12+0.008*X1**0.1504)
     >        *(1-X2)*(1-X2)*(1.0-X2**(2.0/3.0))
      endif  
c     
c     calculate c_const
c     
      c_const = 1.0 / ( 1.0 + 1.0e13 * exp(-2.45/tev)) 
c     
c     calculate csp3
c     
      csp3 = c_const * (2.0e-32 * dflux + exp(-1.7/tev))/
     >      (2.0e-32 * dflux 
     >       + (1.0 + 2e29/dflux * exp(-1.8/tev))*exp(-1.7/tev))      
c     
c     Calculate Ysurf
c     
      ysurf = csp3 * yphys2 / (1.0 + exp((e0-90)/50))
c
c     Calculate Ytherm
c      
      ytherm = csp3 * 0.033 * exp(-1.7/tev) /
     >         (2.0e-32 * dflux + exp(-1.7/tev)) 
c
c     Calculate Ychem
c
      yldchem96 = ysurf + ytherm * (1.0 + d_const * yphys)
c
      return
      end
c
c
c
      subroutine test_phys_yld(matt,matp)
      implicit none
      integer matt,matp

      real YIELD
      external yield

      real yself,yphys
      real energy
      integer in

      ! loop energy from 5 to 500 eV in 5eV increments and calculate
      ! both the self sputtering and matp-> matt

      write(6,'(a)') 'CALCULATION OF TEST YIELDS'//
     >               ' FOR CASE MATERIALS:'
      write(6,'(a)') '  ENERGY (eV)      Yphys        Yself'

      do in = 1,200

         energy = in * 5.0

         ! assume te=ti=0.0 for these print outs
         yphys = yield(matp,matt,energy,0.0,0.0) 
         yself = yield(6,matt,energy,0.0,0.0)

         write(6,'(a,3(1x,g12.5))') 'YIELD:',energy,yphys,yself
         
      end do
      return
      end
c
c
c
      real function yldchem (e0,tmpdflux,matp,matt,tsurf)
      implicit none
      real e0,tmpdflux,tsurf
      integer matp,matt
c
      include 'params'
c
      include 'comtor'
c
c     YLDCHEM is being used to access all of the different 
c     chemical sputtering yield formulae - so that changes
c     throughout the DIVIMP code will not be required. 
c
c     YLDCHEM: This function calculates the chemical sputtering
c              yield for Deuterium on Graphite based on the formulae
c              in the publication "Revised Formula for Chemical
c              Sputtering of Carbon", C. Garcia-Rosales, J. Roth.
c
c              The yield is composed of a temperature dependent
c              component and an impact energy dependent component.
c
c              David Elder, 1994 Oct 27
c
      real ychth,ychath,yphys,dflux
      real tev,t
      real yield,ychem,yldchem96
      external yield,yldchem96
c
c     Test for negative dflux - convert to positive and issue 
c     warning.      
c
      if (tmpdflux.lt.0.0) then 
         dflux = abs(tmpdflux) 
         write(6,*) 'WARNING: Negative flux found in YLDCHEM :',
     >               tmpdflux  
      else
         dflux = tmpdflux
      endif 
c
c     Needed quantities
c
      t = tsurf
c
      if (matt.ne.4) then
         write (6,*) 'Error: Chemical sputtering called'//
     >               ' for non-carbon target'
         yldchem = 0.0
         return
      endif
c
c     YLDCHEM is being used to access all of the different 
c     chemical sputtering yield formulae - so that changes
c     throughout the DIVIMP code will not be required. 
c
c
c     Alternate Garcia-Rosales/Roth (1996) - DIVIMP
c
      if (cchemopt.eq.2) then  
         yldchem = yldchem96(e0,tmpdflux,matp,matt,t)
         return
c
c     Houyang's Chemical Sputtering Options. 
c
      elseif (cchemopt.gt.2.and.cchemopt.le.8) then 
c
c        Houyang's options are numbered the same - except
c        ABOVE the two DIVIMP options. 
c
c        Convert fluxes to cm-2/sec
c
         CALL SPUTCHEM(cchemopt-2,E0,t,dFLUX/1.0e4,YCHEM,crmb) 
         yldchem = ychem
         return
c
c     Constant Yield   
c
      elseif (cchemopt.eq.9) then 
         yldchem = const_yield
         return
c
c     Modification to one of Houyang's options by Tom Rognlien
c
      elseif (cchemopt.eq.10.or.cchemopt.eq.11) then 
c
c        Options in SPUTCHEM are numbered the same - except
c        ABOVE the two DIVIMP options. 
c
c        Convert fluxes to cm-2 s-1
c
         CALL SPUTCHEM(cchemopt-2,E0,t,dFLUX/1.0e4,YCHEM,crmb) 
         yldchem = ychem
         return
      endif 
c
c
c     Original Chemical Sputtering Option continues here.
c
      tev = t  * ( 1.38e-23 / 1.6e-19 )
c
      yphys = yield (matp,matt,e0,0.0,0.0)
c
c     Thermal chemical yield
c
c     Note that the formulae from Garcia-Rosales and Roth use the
c     deuterium flux in cm-2s-1 and so it is necessary to convert the
c     DIVIMP calculated dflux which is in terms of m-2s-1 by dividing
c     it by 1.0e4.
c
c
      ychth = (6.0e19 * exp(-1.0/tev)) /
     >        (1.0e15 + 3.0e27 * exp(-2.0/tev)) *
     >        (2.0 + 200.0 * yphys) * (1.0e16/(dflux/1.0e4))**(0.1)
c
c     Athermal Chemical Yield
c
      ychath = 0.05 * exp( e0 * 1.0e-3 * (20.0 - 1.0/tev)) /
     >   ((1.0+exp((e0-150.0)/25.0)) * (1.0+exp((t-740.0)/25.0)))
c
c     Total chemical yield
c
      yldchem = ychth + ychath
c
      return
      end
c
c 
      SUBROUTINE SPUTCHEM(IOPTCHEM,E0,TEMP,FLUX,YCHEM,mb)

      IMPLICIT NONE
      INTEGER  IOPTCHEM
      REAL     E0,TEMP,FLUX,YCHEM,YGARCIA,YHAASZ,YROTH96,
     >         YHAASZ97,yhaasz97m,yhaasz97fm,mb
      INTRINSIC MAX
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FOR D --> C                                  *
C  *                                                                   *  
C  *  IOPTCHEM       -  Options for chemical sputtering:               *
C  *         1       -  Garcia-Rosales' formula (EPS94)                *
C  *         2       -  according to Pospieszczyk (EPS95)              *
C  *         3       -  Vietzke (in Physical processes of the inter-   *
C  *                    action of Fusion Plasmas with Solids)          *
C  *         4       -  Haasz (Submitted to J.Nucl.Mater.,Dec. 1995)   *
C  *         5       -  Roth & Garcia-Rosales (Submitted to Nucl.      *
C  *                    Fusion, March 1996)                            *
C  *         6       -  Haasz 1997 (Brian Mech's PhD Thesis)           *
c  *         8       -  Modified Haasz 1997 - projected to low energy  *
c  *
C  *                                                                   *
C  *  E0  (eV)       -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *  FLUX (cm-2s-1) -  Ion or neutral flux                            *
C  *  YCHEM          -  Chemical Sputtering yield                      *
c  *  MB             -  Background ion mass - AMU                      *
C  *                                                                   *
C  *********************************************************************
C
 
      IF      (IOPTCHEM.EQ.1) THEN
        YCHEM = YGARCIA(E0,TEMP,FLUX)
      ELSE IF (IOPTCHEM.EQ.2) THEN
        YCHEM = 0.04254*(MAX(5E18,FLUX)/5E18)**(-0.477)
      ELSE IF (IOPTCHEM.EQ.3) THEN
        YCHEM = 0.0215*(MAX(1E14,FLUX)/1E16)**(-0.1)
      ELSE IF (IOPTCHEM.EQ.4) THEN
        YCHEM = YHAASZ(E0,TEMP)
      ELSE IF (IOPTCHEM.EQ.5) THEN
c
c       Perform unit conversion to m-2 s-1 in the call so that the
c       value of "flux" in the calling routine is not changed.  
c
c       FLUX  = FLUX*1E4                        
c
        YCHEM = YROTH96(E0,TEMP,FLUX*1.0e4)
      ELSE IF (IOPTCHEM.EQ.6) THEN
        YCHEM = YHAASZ97(E0,TEMP,mb)
      ELSE IF (IOPTCHEM.EQ.8) THEN
        YCHEM = YHAASZ97M(E0,TEMP,mb)
      ELSE IF (IOPTCHEM.EQ.9) THEN
        YCHEM = YHAASZ97FM(E0,TEMP,mb,flux*1.0e4)
      END IF

      RETURN
      END

C -------------

      FUNCTION YROTH96(E0,TEMP,FLUX)

      IMPLICIT NONE
      REAL    E0,TEMP,FLUX
      REAL    ETHC,ETFC,QC,SNC
      REAL    CSURF,CSP3
      REAL    YPHYS,YSURF,YTHERM,YROTH96
      INTRINSIC MIN
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING CALCULATED BY Garcia-Rosales' FORMULA        *
C  *                                                                   *  
C  *  ETHC (eV)  -  Threshold energy for D -> C physical sputtering    *
C  *  ETFC (eV)  -  Thomas-Fermi energy                                *
C  *  SNC        -  Stopping power                                     *
C  *  QC         -  Fitting parameters                                 *
C  *                                                                   *
C  *  CSURF      -  Fitting parameters                                 *
C  *  CSP3       -  Carbon at surface                                  *
C  *                                                                   *
C  *  YPHYS      -  Physical sputtering yield                          *
C  *  YSURF      -  Sputtering due to SURFACE process                  *
C  *  YTHERM     -  Sputtering due to THERMAL process                  *
C  *  YROTH96    -  Total CHEMICAL sputtering yield                    *
C  *                                                                   *
C  *********************************************************************

C ---------------------------------------------------------
C Total Chemical Sputtering Yield:
C            Ychem = Ysurf+Ytherm*(1+D*Yphys)
C ---------------------------------------------------------

C 
C 1> PHYSICAL SPUTTERING YIELD
C
      ETHC = 27.0
      ETFC = 447.0      
      QC   = 0.1
C
C  - Stopping Power      
C
      SNC = 0.5*LOG(1.+1.2288*E0/ETFC)/(E0/ETFC
     >    + 0.1728*SQRT(E0/ETFC)
     >    + 0.008*(E0/ETFC)**0.1504)
C
C  - Physical Sputtering Yield
C
      IF (E0.GT.ETHC) THEN
         YPHYS = QC*SNC*(1-(ETHC/E0)**(2./3.))*(1-ETHC/E0)**2
      ELSE
         YPHYS = 0.0
      ENDIF
C
C 2> YSURF: Surface Process 
C
      CSURF  = 1/(1.+1E13*EXP(-2.45*11604/TEMP))
      CSP3   = CSURF*(2E-32*FLUX+EXP(-1.7*11604/TEMP))
     >         /(2E-32*FLUX+(1+2E29/FLUX*EXP(-1.8*11604/TEMP))
     >         *EXP(-1.7*11604/TEMP))

      IF (E0.GT.1.) THEN
         YSURF = CSP3*QC*SNC*(1-(1./E0)**(2./3.))*(1-1./E0)**2
     >           /(1.+EXP((MIN(90.0,E0)-90.)/50.))
      ELSE
         YSURF = 0.0
      ENDIF
C
C 3> YTHERM: Thermak Activated Process
C
      YTHERM = CSP3*0.033*EXP(-1.7*11604/TEMP)
     >         /(2E-32*FLUX+EXP(-1.7*11604/TEMP))
C
C 4> YCHEM: Total Chemical Sputtering Yield
C
      YROTH96 = YSURF + YTHERM * (1 + 125 * YPHYS)

CW    WRITE(6,*) 'YROTH96 = ',YROTH96

      RETURN
      END
      
      
C -------------

      FUNCTION YGARCIA(E0,TEMP,FLUX)

      IMPLICIT NONE
      REAL    E0,TEMP,FLUX
      REAL    ETHC,ETFC,QC,SNC
      REAL    YPHYS,YCHEM_TH,YCHEM_ATH,YGARCIA
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING CALCULATED BY Garcia-Rosales' FORMULA        *
C  *                                                                   *  
C  *  ETHC (eV)  -  Threshold energy for D -> C physical sputtering    *
C  *  ETFC (eV)  -  Thomas-Fermi energy                                *
C  *  SNC        -  Stopping power                                     *
C  *  QC         -  Fitting parameters                                 *
C  *                                                                   *
C  *  YPHYS      -  Physical sputtering yield                          *
C  *  YCHEM_TH   -  Thermal activated mechanism                        *
C  *  YCHEM_ATH  -  Athermal mechanism                                 *
C  *                                                                   *
C  *********************************************************************
C            

      ETHC = 27.0
      ETFC = 447.0      
      QC   = 0.1
C
C Check for impact energies below threshold
C
      IF (E0.GT.ETHC) THEN
C
C Stopping Power      
C
      SNC = 0.5*LOG(1.+1.2288*E0/ETFC)/(E0/ETFC
     >    + 0.1728*SQRT(E0/ETFC)
     >    + 0.008*(E0/ETFC)**0.1504)
C 
C Physical Sputtering Yield
C
         YPHYS = QC*SNC*(1-(ETHC/E0)**(2./3.))*(1-ETHC/E0)**2
C
      ELSE
         YPHYS = 0.0
      ENDIF
C
C Chemical Sputtering Yield
C
      YCHEM_TH = 6E19*EXP(-1.*11604/TEMP)
     >         /(1E15+3E27*EXP(-2.*11604/TEMP))
     >         * (2.+200*YPHYS)*(MAX(1E16,FLUX)/1E16)**(-0.1)

      YCHEM_ATH = 0.05*EXP(E0*1E-3*(20.-1*11604./TEMP))
     >          / ((1.+EXP((E0-150.)/25.))*(1.+EXP((TEMP-740.)/25.)))

      YGARCIA = YCHEM_TH + YCHEM_ATH

CW    WRITE(6,*) 'YTH = ',YCHEM_TH,'YATH = ',YCHEM_ATH

      RETURN
      END
      
      
C -------------


      FUNCTION YHAASZ(E0,TEMP)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (Dec. 1995)            *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE

      REAL    E0,TEMP
      REAL    FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >        FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >        FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4),
     >        FITC1050(4),FITC1100(4)
      REAL    POLY_C(4),YFIT,FITE0
      REAL    YHAASZ
      INTEGER I
C
C     Poly. fit c. /       a0,      a1,      a2,      a3
C
      DATA FITC300 / -0.01789, 0.02309, 0.00089,-0.00315/
      DATA FITC350 / -0.01691, 0.02020, 0.00451,-0.00407/
      DATA FITC400 / -0.01128, 0.01230, 0.00922,-0.00483/
      DATA FITC450 / -0.00401, 0.00453, 0.01226,-0.00493/
      DATA FITC500 / -0.01000, 0.02097,-0.00032,-0.00153/
      DATA FITC550 / -0.02022, 0.04019,-0.01430, 0.00253/
      DATA FITC600 /  0.00047,-0.00319, 0.00950,-0.00025/
      DATA FITC650 /  0.02921,-0.05657, 0.03467,-0.00226/
      DATA FITC700 /  0.00561,-0.00081,-0.01044, 0.00939/
      DATA FITC750 /  0.00225, 0.00205,-0.00949, 0.00800/
      DATA FITC800 /  0.00900,-0.02109, 0.01366, 0.00048/
      DATA FITC850 /  0.00483,-0.01691, 0.01513,-0.00152/
      DATA FITC900 /  0.00569,-0.02211, 0.02185,-0.00427/
      DATA FITC950 /  0.00317,-0.01827, 0.02081,-0.00482/
      DATA FITC1000/  0.00436,-0.02075, 0.02290,-0.00574/
      DATA FITC1050/  0.00463,-0.02082, 0.02285,-0.00601/
      DATA FITC1100/  0.00920,-0.02942, 0.02802,-0.00723/
C
C Find right polynomial fit coefficients for a given temperature
C
      IF      (TEMP.LE.300.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
         ENDDO
      ELSE IF (TEMP.GT.300.0 .AND. TEMP.LE.350.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
         ENDDO
      ELSE IF (TEMP.GT.350.0 .AND. TEMP.LE.400.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
         ENDDO
      ELSE IF (TEMP.GT.400.0 .AND. TEMP.LE.450.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
         ENDDO
      ELSE IF (TEMP.GT.450.0 .AND. TEMP.LE.500.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
         ENDDO
      ELSE IF (TEMP.GT.500.0 .AND. TEMP.LE.550.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
         ENDDO
      ELSE IF (TEMP.GT.550.0 .AND. TEMP.LE.600.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
         ENDDO
      ELSE IF (TEMP.GT.600.0 .AND. TEMP.LE.650.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
         ENDDO
      ELSE IF (TEMP.GT.650.0 .AND. TEMP.LE.700.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
         ENDDO
      ELSE IF (TEMP.GT.700.0 .AND. TEMP.LE.750.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
         ENDDO
      ELSE IF (TEMP.GT.750.0 .AND. TEMP.LE.800.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
         ENDDO
      ELSE IF (TEMP.GT.800.0 .AND. TEMP.LE.850.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
         ENDDO
      ELSE IF (TEMP.GT.850.0 .AND. TEMP.LE.900.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
         ENDDO
      ELSE IF (TEMP.GT.900.0 .AND. TEMP.LE.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
         ENDDO
      ELSE IF (TEMP.GT.950.0 .AND. TEMP.LE.1000.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1000(I)
         ENDDO
      ELSE IF (TEMP.GT.1000.0 .AND. TEMP.LE.1050.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1050(I)
         ENDDO
      ELSE
         DO I = 1,4
           POLY_C(I) = FITC1100(I)
         ENDDO
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.0)  THEN
         FITE0 = 10.
      ELSE IF (E0.GT.200.0) THEN
         FITE0 = 200.
      ELSE
         FITE0 = E0
      ENDIF 
C
      YFIT = 0.0
      DO I = 1,4
        YFIT = YFIT + POLY_C(I)*ALOG10(FITE0)**(I-1)
      ENDDO

      YHAASZ = YFIT

CW    WRITE(6,*) 'YHAASZ = ',YHAASZ

      RETURN
      END
c
c
c
C -------------


      FUNCTION YHAASZ97(E0,TEMP,mb)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE
c
c     Routine Originally written by Houyang Guo 
c     - fits contained here are for total C/D+ yields as taken from 
c       Brian Mech's PhD thesis data. The experimental data is for
c       a flux of 1e18 /m2
c     - The mass dependence adjustment to obtain approximate fits for
c       C/H+ was suggested by Jim Davis and implemented by David Elder.
c       The mass correction factor is:
c       E0 <= 50eV  - Y = Y * sqrt (mH/mD)
c       E0 >= 100eV - Y = Y * mH/mD
c       50eV < E0 < 100eV - Y= Y * (sqrt(mH/mD) + (mH/mD - sqrt(mH/mD)) * (E0-50)/50)
c


      REAL    E0,TEMP,mb
      REAL    FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >        FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >        FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4)
      REAL    POLY_C(4),YFIT,FITE0
      REAL    YHAASZ97,mass_cor
      INTEGER I
C
C     Poly. fit c. /       a0,      a1,      a2,      a3
C
      DATA FITC300 / -0.03882, 0.07432,-0.03470, 0.00486/
      DATA FITC350 / -0.05185, 0.10126,-0.05065, 0.00797/
      DATA FITC400 / -0.06089, 0.12186,-0.06240, 0.01017/
      DATA FITC450 / -0.08065, 0.16884,-0.09224, 0.01625/
      DATA FITC500 / -0.08872, 0.19424,-0.10858, 0.01988/
      DATA FITC550 / -0.08728, 0.20002,-0.11420, 0.02230/
      DATA FITC600 / -0.05106, 0.13146,-0.07514, 0.01706/
      DATA FITC650 /  0.07373,-0.13263, 0.09571,-0.01672/
      DATA FITC700 /  0.02722,-0.03599, 0.02064, 0.00282/
      DATA FITC750 /  0.09052,-0.18253, 0.12362,-0.02109/
      DATA FITC800 /  0.02604,-0.05480, 0.04025,-0.00484/
      DATA FITC850 /  0.03478,-0.08537, 0.06883,-0.01404/
      DATA FITC900 /  0.02173,-0.06399, 0.05862,-0.01380/
      DATA FITC950 / -0.00086,-0.01858, 0.02897,-0.00829/
      DATA FITC1000/ -0.01551, 0.01359, 0.00600,-0.00353/
C
C Find right polynomial fit coefficients for a given temperature
C
      IF      (TEMP.LE.300.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
         ENDDO
      ELSE IF (TEMP.GT.300.0 .AND. TEMP.LE.350.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
         ENDDO
      ELSE IF (TEMP.GT.350.0 .AND. TEMP.LE.400.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
         ENDDO
      ELSE IF (TEMP.GT.400.0 .AND. TEMP.LE.450.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
         ENDDO
      ELSE IF (TEMP.GT.450.0 .AND. TEMP.LE.500.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
         ENDDO
      ELSE IF (TEMP.GT.500.0 .AND. TEMP.LE.550.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
         ENDDO
      ELSE IF (TEMP.GT.550.0 .AND. TEMP.LE.600.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
         ENDDO
      ELSE IF (TEMP.GT.600.0 .AND. TEMP.LE.650.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
         ENDDO
      ELSE IF (TEMP.GT.650.0 .AND. TEMP.LE.700.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
         ENDDO
      ELSE IF (TEMP.GT.700.0 .AND. TEMP.LE.750.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
         ENDDO
      ELSE IF (TEMP.GT.750.0 .AND. TEMP.LE.800.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
         ENDDO
      ELSE IF (TEMP.GT.800.0 .AND. TEMP.LE.850.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
         ENDDO
      ELSE IF (TEMP.GT.850.0 .AND. TEMP.LE.900.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
         ENDDO
      ELSE IF (TEMP.GT.900.0 .AND. TEMP.LE.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
         ENDDO
      ELSE IF (TEMP.GT.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1000(I)
         ENDDO
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.0)  THEN
         FITE0 = 10.
      ELSE IF (E0.GT.200.0) THEN
         FITE0 = 200.
      ELSE
         FITE0 = E0
      ENDIF 
C
      YFIT = 0.0
      DO I = 1,4
        YFIT = YFIT + POLY_C(I)*ALOG10(FITE0)**(I-1)
      ENDDO
c
c     Mass correction factor 
c
      if (mb.eq.2.0) then 
         mass_cor = 1.0
      elseif (e0.le.50.0) then  
         mass_cor = sqrt(mb/2.0) 
      elseif(e0.gt.50.0.and.e0.lt.100.0) then 
         mass_cor = sqrt(mb/2.0) + (mb/2.0-sqrt(mb/2.0))
     >                   * (e0-50.0)/50.0
      elseif(e0.ge.100.0) then 
         mass_cor = mb/2.0
      endif
c

      YHAASZ97 = YFIT * mass_cor

CW    WRITE(6,*) 'YHAASZ97 = ',YHAASZ97

      RETURN
      END
c
c
c
      FUNCTION YHAASZ97M(E0,TEMP,mb)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  with the addition of a new fit below 10 eV as suggested by       *
C  *  J.Davis and parameterized by G. Porter; now interpolates between *
C  *  5 and 10 eV to lower value (YDAVIS98), and is fixed below 5 eV   *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE

      real E0, TEMP,mb
      real YHAASZ97M, YDAVIS98, YHAASZ97
      real mass_cor
      real m1,m2,m3,reducf,FRAC

      DATA m1/602.39/, m2/202.24/, m3/43.561/, reducf/0.2/
c
c     Mass correction is left out of the modifications below 10eV since it 
c     isn't clear that the behaviour of atmomic H and D differ substantially.
c     So the routine here will not mass correct below 10eV and will just
c     go linearly from the mass corrected value at 10eV to the 
c     constant (atomic) level at 5eV.
c

      mass_cor = 1.0

c
c     Mass correction factor 
c
c      if (mb.eq.2.0) then 
c         mass_cor = 1.0
c      elseif (e0.le.50.0) then  
c         mass_cor = sqrt(mb/2.0) 
c      elseif(e0.gt.50.0.and.e0.lt.100.0) then 
c         mass_cor = sqrt(mb/2.0) + (mb/2.0-sqrt(mb/2.0))
c     >                   * (e0-50.0)/50.0
c      elseif(e0.ge.100.0) then 
c         mass_cor = mb/2.0
c      endif
c
      IF (E0 .GE. 10) THEN
         YHAASZ97M = YHAASZ97(E0,TEMP,mb)
      ELSEIF (E0 .LT. 10. .AND. E0 .GE. 5.) THEN
         FRAC = (E0-5.)/5.
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1)**2 + m3) * mass_cor
         YHAASZ97M = FRAC*YHAASZ97(E0,TEMP,mb)+ (1.-FRAC)*YDAVIS98
      ELSEIF (E0 .LT. 5.) THEN
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1)**2 + m3) * mass_cor
         YHAASZ97M = YDAVIS98   
      ENDIF

      RETURN
      END
c
c
c
      FUNCTION YHAASZ97FM(E0,TEMP,mb,flux)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  with the addition of a new fit below 10 eV as suggested by       *
C  *  J.Davis and parameterized by G. Porter; now interpolates between *
C  *  5 and 10 eV to lower value (YDAVIS98), and is fixed below 5 eV   *
c  *                                                                   *
c  *  YHAASZ97FM - this routine has also been modified to include      *
c  *               a flux dependence. Flux is in units of m-2 s-1.     *
c  *               Flux dependence involves reducing the effective     *
c  *               surface temperature and thus shifting the profiles  *
C  *                                                                   *
C  *********************************************************************
c

c
C            
      IMPLICIT NONE

      real E0, TEMP,mb,flux
      real YHAASZ97FM, YDAVIS98, YHAASZ97
      real mass_cor 
      real temp_eff
      real m1,m2,m3,reducf,FRAC

      DATA m1/602.39/, m2/202.24/, m3/43.561/, reducf/0.2/
c
c     The flux dependence in this case involves just shifting the 
c     effective temperature at which the yield curves are 
c     appllied depending on the flux. A maximum shift of 100K
c     is applied for fluxes over 1e20. A zero shift is applied for
c     fluxes under 1e19 and linear between these.     
c
      if (flux.le.1.0e19) then 

         temp_eff = temp

      elseif (flux.ge.1.0e20) then 

         temp_eff = temp - 100.0

      else

         temp_eff = temp - 100.0 * (flux-1.0e19)/9.0e19

      endif
c
c     Mass correction is left out of the modifications below 10eV since it 
c     isn't clear that the behaviour of atmomic H and D differ substantially.
c     So the routine here will not mass correct below 10eV and will just
c     go linearly from the mass corrected value at 10eV to the 
c     constant (atomic) level at 5eV.
c

      mass_cor = 1.0
c
c     Mass correction factor 
c
c      if (mb.eq.2.0) then 
c         mass_cor = 1.0
c      elseif (e0.le.50.0) then  
c         mass_cor = sqrt(mb/2.0) 
c      elseif(e0.gt.50.0.and.e0.lt.100.0) then 
c         mass_cor = sqrt(mb/2.0) + (mb/2.0-sqrt(mb/2.0))
c     >                   * (e0-50.0)/50.0
c      elseif(e0.ge.100.0) then 
c         mass_cor = mb/2.0
c      endif
c
c
c
      IF (E0 .GE. 10) THEN
         YHAASZ97FM = YHAASZ97(E0,TEMP_eff,mb)
      ELSEIF (E0 .LT. 10. .AND. E0 .GE. 5.) THEN
         FRAC = (E0-5.)/5.
         YDAVIS98 = reducf/(m2*((TEMP_eff/m1)**2 - 1)**2 + m3)
     >              * mass_cor
         YHAASZ97FM = FRAC*YHAASZ97(E0,TEMP_eff,mb)+ (1.-FRAC)*YDAVIS98
      ELSEIF (E0 .LT. 5.) THEN
         YDAVIS98 = reducf/(m2*((TEMP_eff/m1)**2 - 1)**2 + m3) 
     >              * mass_cor
         YHAASZ97FM = YDAVIS98   
      ENDIF


      RETURN
      END
