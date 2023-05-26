module mod_cx

  implicit none

      !     WF'96: THIS CLUSTER OF ROUTINES WAS ADAPTED FROM NIMBUS,
      !     WITH THE AIM OF CALCULATING THE ION-NEUTRAL
      !     CHARGE EXCHANGE COLLISION RATE. CALLED VIA

      !     pmomloss() =>  estscx() => scx_s() => cxsig()

      !     CHARGE EXCHANGE   H + H+ -> H+ + H

      ! I   IZ       = ATOMIC NUMBER
      ! I   AN       = MASS OF NEUTRAL (AMU)
      ! I   EN       = NEUTRAL ENERGY (EV)
      ! I   WNX/Y/Z  = NEUTRAL DIRECTION COSINES
      ! I   TI       = ION TEMPERATURE (EV)
      ! I   vp       = velocity of PLASMA FLOW
      ! I   AI       = ION MASS (AMU)
      ! I   WPX/Y/Z  = FLOW DIRECTION COSINES
      ! I   IXTYPE   = C.X. CROSS SECTION MODEL
      ! O   SIGCX    = C.X. CROSS SECTION (EFFECTIVE)
      ! O   SVCX     = C.X. REACTION RATE

      !     RETURNS MICROSCOPIC CHARGE EXCHANGE X-SECTION  BETWEEN H-ISOTOPES
      !     AVERAGED OVER ALL TARGET VELOCITIES (MAXWELL + FLOW) FOR NEUTRALS
      !     IN A PLASMA REGION CONTAINING H+,D+,DT+,T+.

      !     NB:  originally in cgs;  converted to mks / eV , both in/output.

  ! moved common block variables to module variables
  real*8 :: tis,c1,c2,c3,c4
  integer :: idqh

contains


    subroutine CXSIG (IZ,AN,EN,WNX,WNY,WNZ,TI,VP,AI,WPX,WPY,WPZ,IXTYPE,SIGCX,SVCX)

      implicit none
      integer ip,ixtype
      real*8 IZ,AN,EN,WNX,WNY,WNZ,TI,VP,AI,WPX,WPY,WPZ,SIGCX,SVCX

      real*8 vn,vnx,vny,vnz,vpx,vpy,vpz,vr2,ES,ENS,ESTS,C,RRATE,ER,RM
      !COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      !real*8 tis,c1,c2,c3,c4
      !integer idqh

      !real*8 xtfct
      !external xtfct

      !     VELOCITY OF THE NEUTRAL

      IF(IZ.ge.2.0) GOTO 300
      vn = SQRT(EN/AN) * 1.3841E+06 * 0.01
      vnx = vn * WNX
      vny = vn * WNY

      !      write(6,*) 'wnxy',wnx,wny,wnz
      !      write(6,*) 'vn,xyz',vn,vnx,vny,vnz

      !     VELOCITY OF THE PLASMA

      !      vp = SQRT((TE+TI)/AI) * 0.978E+06 * MACH

      vnz = vn * WNZ
      vp  = abs (vp)
      vpx = vp * WPX
      vpy = vp * WPY

      !      write(6,*) 'Wxyz',WPX,WPY,WPZ
      !      write(6,*) 'vp,xyz',vp,vpx,vpy,vpz

      !     MEAN RELATIVE VELOCITY BETWEEN PLASMA AND NEUTRAL

      vpz = vp * WPZ
      IF(VP.LE.0.0) THEN
         vr2 = vn ** 2
      ELSE
         vr2 = (vpx-vnx)**2 + (vpy-vny)**2 + (vpz-vnz)**2

         !     SCALED SHIFTED NEUTRAL ENERGY (EV)

      ENDIF

      !     SHIFTED NEUTRAL ENERGY

      ES = 0.52197E-12 * vr2 * 1.0E+4

      !     SCALED ION TEMPERATURE

      ENS = ES * AN

      !     RATIO OF SHIFTED NEUTRAL ENERGY TO SCALED ION TEMP.

      TIS = TI / AI

      ESTS = ES / TIS
      IF (ES.GE.4.0E+4) GO TO 20
      IF (ES.GE.5.0E+3) GO TO 10
      IF (ESTS.GE.5.0) GO TO 80
      IP = 1
      GO TO 30
10    IF(ESTS.GE.26.0) GO TO 80
      IP = 2
      GO TO 30
20    IF(ESTS.GE.100.0) GO TO 80
      IP = 3
30    C = SQRT(AI/AN)
      C1 = SQRT(ENS/TI)
      C2 = C*C1
      C3 = 0.0
      IF (C2**2.LT.174.0) C3 = EXP(-C2**2)
      C4 = 2.0*C2
      IDQH = 1
      GO TO ( 40 , 50 , 60 ),IP
40    CALL DQH04P(XTFCT,RM,ixtype)
      GO TO 70
50    CALL DQH12P(XTFCT,RM,ixtype)
      GO TO 70
60    CALL DQH32P(XTFCT,RM,ixtype)

70    RRATE = 0.780939E+6 * RM * TI/AI * SQRT(AN/ENS)* 1.0E-6
      if (vn.gt.0.0) then
         SIGCX = RRATE / vn
      else
         !         write (6,*) 'ERROR: SOLASCV - Momentum Option 9 :'//
         !     >               ' Vn =< 0.0 ',vn,sigcx
         SIGCX = 0.0

      endif

      !     MEAN RELATIVE NEUTRAL/ION ENERGY

      GOTO 90

80    ER = 1.5 * AN * TIS + ENS

      !     CONSERVE THE REACTION RATE

      CALL HXHP(ER/AN,ixtype,SIGCX)
      SIGCX = SIGCX * SQRT(ER/EN) * 1.0E-4
90    CONTINUE
      SVCX = SIGCX * vn

      RETURN
300   WRITE(6,*) 'ERROR: NON-HYDROGENIC CALL FROM CXSIG'
      write(6,*) 'iz=', iz
      write(6,*) 'var...', ip,ixtype,IZ,AN,EN,WNX,WNY,WNZ,TI,VP,AI,WPX,WPY,WPZ,SIGCX,SVCX,vn,vnx,vny,vnz,vpx,vpy,vpz,&
           vr2,ES,ENS,ESTS,C,RRATE,ER,RM
      STOP

      RETURN

    END subroutine CXSIG


    SUBROUTINE HXHP(ER,IXTYP,SIGMA)

      !     RETURNS MICROSCOPIC CROSS SECTION SIGMA FOR CHARGE-EXCHANGE
      !     REACTION BETWEEN HYDROGEN ATOMS/IONS

      !     ER=RELATIVE ENERGY

      !                    MICROSCOPIC CROSS SECTION OF REACTION
      !                    (H+)+H-->H+(H+) IS GIVEN

      !     ICXTYP=1:JANEV - 'ELEMENTARY PROCESSES IN HYDROGEN-HELIUM PLASMA',
      !                       SPRINGER (1987))



      implicit none
      real*8 ER,SIGMA,ESCALE
      integer ixtyp

      !external PNFIT

      !        DIMENSION R318(9)
      real*8 R318(9)

      DATA R318/-3.274123792568E+01,-8.916456579806E-02,-3.016990732025E-02,9.205482406462E-03,2.400266568315E-03,&
           -1.927122311323E-03,3.654750340106E-04,-2.788866460622E-05,7.422296363524E-07/
      !                                   JANEV (1987)
      IF( IXTYP.EQ.1 ) THEN
         ESCALE = ER
         IF(ESCALE.LT.0.1) ESCALE = 0.1
         CALL PNFIT(9,R318,ESCALE,SIGMA)
      ELSE
         ESCALE = ER
         IF(ESCALE.GT.200.0) GO TO 10
         !                                   GREENLAND,1984 (E.LE.200 EV)

         !         Caution: Alog10 was used previouly

         IF(ESCALE.LT.0.001) ESCALE=0.001
         SIGMA = 4.549E-15 - 1.033E-15 * LOG10(ESCALE)
         !                                   RIVIERE  (E.GT.200 EV)
         RETURN
10       SIGMA = 6.937E-15 * (1.0 - 0.155*LOG10(ESCALE))**2
         IF(ESCALE.GE. 15000.0) THEN
            SIGMA = SIGMA/(1.0 + 0.1112E-14 * (ESCALE**3.3))
         ENDIF

      END IF
      RETURN

    END SUBROUTINE HXHP



    SUBROUTINE DQH04P(FCT,Y,ixtype)

      !     ..................................................................

      !        SUBROUTINE DQH04P

      !        PURPOSE
      !           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
      !                               0 TO +INFINITY).

      !        USAGE
      !           CALL DQH04P (FCT,Y)
      !           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT

      !        DESCRIPTION OF PARAMETERS
      !           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
      !                    SUBPROGRAM USED.
      !           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.

      !        REMARKS
      !           NONE

      !        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
      !           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
      !           MUST BE FURNISHED BY THE USER.

      !        METHOD
      !           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSSIAN-HERMITE
      !           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
      !           FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.
      !           FOR REFERENCE, SEE
      !           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
      !           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
      !           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
      !           TR00.1100 (MARCH 1964), PP.213-214.

      !     ..................................................................


      implicit none
      REAL*8  X,Y,FCT

      integer ixtype
      !COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      !real*8 tis,c1,c2,c3,c4
      !integer idqh

      EXTERNAL FCT
      X=.29306374202572440D1
      Y=.19960407221136762D-3*FCT(X,ixtype)
      X=.19816567566958429D1
      Y=Y+.17077983007413475D-1*FCT(X,ixtype)
      X=.11571937124467802D1
      Y=Y+.20780232581489188D0*FCT(X,ixtype)
      X=.38118699020732212D0
      Y=Y+.66114701255824129D0*FCT(X,ixtype)
      RETURN

    END SUBROUTINE DQH04P


    SUBROUTINE DQH12P(FCT,Y,ixtype)
      !.......................................................................
      !.......................................................................


      !        SUBROUTINE DQH12P

      !        PURPOSE
      !           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
      !                               0 TO +INFINITY).

      !        USAGE
      !           CALL DQH12P (FCT,Y)
      !           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT

      !        DESCRIPTION OF PARAMETERS
      !           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
      !                    SUBPROGRAM USED.
      !           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.

      !        REMARKS
      !           NONE

      !        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
      !           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
      !           MUST BE FURNISHED BY THE USER.

      !        METHOD
      !           EVALUATION IS DONE BY MEANS OF 12-POINT GAUSSIAN-HERMITE
      !           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
      !           FCT(X) IS A POLYNOMIAL UP TO DEGREE 47.
      !           FOR REFERENCE, SEE
      !           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
      !           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
      !           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
      !           TR00.1100 (MARCH 1964), PP.213-214.

      !     ..................................................................

      implicit none
      REAL*8 X,Y,FCT

      integer ixtype
      !COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      !real*8 tis,c1,c2,c3,c4
      !integer idqh

      EXTERNAL FCT
      X=.60159255614257397D1
      Y=.16643684964891089D-15*FCT(X,IXTYPE)
      X=.52593829276680444D1
      Y=Y+.65846202430781701D-12*FCT(X,IXTYPE)
      X=.46256627564237873D1
      Y=Y+.30462542699875639D-9*FCT(X,IXTYPE)
      X=.40536644024481495D1
      Y=Y+.40189711749414297D-7*FCT(X,IXTYPE)
      X=.35200068130345247D1
      Y=Y+.21582457049023336D-5*FCT(X,IXTYPE)
      X=.30125461375655648D1
      Y=Y+.56886916364043798D-4*FCT(X,IXTYPE)
      X=.25238810170114270D1
      Y=Y+.8236924826884175D-3*FCT(X,IXTYPE)
      X=.20490035736616989D1
      Y=Y+.70483558100726710D-2*FCT(X,IXTYPE)
      X=.15842500109616941D1
      Y=Y+.37445470503230746D-1*FCT(X,IXTYPE)
      X=.11267608176112451D1
      Y=Y+.12773962178455916D0*FCT(X,IXTYPE)
      X=.67417110703721224D0
      Y=Y+.28617953534644302D0*FCT(X,IXTYPE)
      X=.22441454747251559D0
      Y=Y+.42693116386869925D0*FCT(X,IXTYPE)
      RETURN

    END SUBROUTINE DQH12P


    SUBROUTINE DQH32P(FCT,Y,ixtype)
      !.......................................................................
      !.......................................................................


      !        SUBROUTINE DQH32P

      !        PURPOSE
      !           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
      !                               0 TO +INFINITY).

      !        USAGE
      !           CALL DQH32P (FCT,Y)
      !           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT

      !        DESCRIPTION OF PARAMETERS
      !           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
      !                    SUBPROGRAM USED.
      !           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.

      !        REMARKS
      !           NONE

      !        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
      !           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
      !           MUST BE FURNISHED BY THE USER.

      !        METHOD
      !           EVALUATION IS DONE BY MEANS OF 32-POINT GAUSSIAN-HERMITE
      !           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
      !           FCT(X) IS A POLYNOMIAL UP TO DEGREE 127.
      !           FOR REFERENCE, SEE
      !           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
      !           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
      !           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
      !           TR00.1100 (MARCH 1964), PP.213-214.

      !     ..................................................................


      implicit none
      REAL*8 X,Y,FCT

      integer ixtype
      !COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      !real*8 tis,c1,c2,c3,c4
      !integer idqh

      EXTERNAL FCT
      X=.10526123167960546D2
      Y=.55357065358569428D-48*FCT(X,IXTYPE)
      X=.9895287586829539D1
      Y=Y+.16797479901081592D-42*FCT(X,IXTYPE)
      X=.9373159549646721D1
      Y=Y+.34211380112557405D-38*FCT(X,IXTYPE)
      X=.8907249099964770D1
      Y=Y+.15573906246297638D-34*FCT(X,IXTYPE)
      X=.8477529083379863D1
      Y=Y+.25496608991129993D-31*FCT(X,IXTYPE)
      X=.8073687285010225D1
      Y=Y+.19291035954649669D-28*FCT(X,IXTYPE)
      X=.7689540164040497D1
      Y=Y+.7861797788925910D-26*FCT(X,IXTYPE)
      X=.7321013032780949D1
      Y=Y+.19117068833006428D-23*FCT(X,IXTYPE)
      X=.69652411205511075D1
      Y=Y+.29828627842798512D-21*FCT(X,IXTYPE)
      X=.66201122626360274D1
      Y=Y+.31522545665037814D-19*FCT(X,IXTYPE)
      X=.62840112287748282D1
      Y=Y+.23518847106758191D-17*FCT(X,IXTYPE)
      X=.59556663267994860D1
      Y=Y+.12800933913224380D-15*FCT(X,IXTYPE)
      X=.56340521643499721D1
      Y=Y+.52186237265908475D-14*FCT(X,IXTYPE)
      X=.53183252246332709D1
      Y=Y+.16283407307097204D-12*FCT(X,IXTYPE)
      X=.50077796021987682D1
      Y=Y+.39591777669477239D-11*FCT(X,IXTYPE)
      X=.47018156474074998D1
      Y=Y+.7615217250145451D-10*FCT(X,IXTYPE)
      X=.43999171682281376D1
      Y=Y+.11736167423215493D-8*FCT(X,IXTYPE)
      X=.41016344745666567D1
      Y=Y+.14651253164761094D-7*FCT(X,IXTYPE)
      X=.38065715139453605D1
      Y=Y+.14955329367272471D-6*FCT(X,IXTYPE)
      X=.35143759357409062D1
      Y=Y+.12583402510311846D-5*FCT(X,IXTYPE)
      X=.32247312919920357D1
      Y=Y+.8788499230850359D-5*FCT(X,IXTYPE)
      X=.29373508230046218D1
      Y=Y+.51259291357862747D-4*FCT(X,IXTYPE)
      X=.26519724354306350D1
      Y=Y+.25098369851306249D-3*FCT(X,IXTYPE)
      X=.23683545886324014D1
      Y=Y+.10363290995075777D-2*FCT(X,IXTYPE)
      X=.20862728798817620D1
      Y=Y+.36225869785344588D-2*FCT(X,IXTYPE)
      X=.18055171714655449D1
      Y=Y+.10756040509879137D-1*FCT(X,IXTYPE)
      X=.15258891402098637D1
      Y=Y+.27203128953688918D-1*FCT(X,IXTYPE)
      X=.12472001569431179D1
      Y=Y+.58739981964099435D-1*FCT(X,IXTYPE)
      X=.9692694230711780D0
      Y=Y+.10849834930618684D0*FCT(X,IXTYPE)
      X=.69192230581004458D0
      Y=Y+.17168584234908370D0*FCT(X,IXTYPE)
      X=.41498882412107868D0
      Y=Y+.23299478606267805D0*FCT(X,IXTYPE)
      X=.13830224498700972D0
      Y=Y+.27137742494130398D0*FCT(X,IXTYPE)
      RETURN

    END SUBROUTINE DQH32P


    SUBROUTINE PNFIT(NJ,A,T,SVJ)


      !     POLYNOMIAL FORM USED FOR MOLECULAR REACTIONS WITH ELECTRONS
      !     AND JANEV C.X. X.S.

      !     SVJ=RATE COEFFICIENT (CM**3/SEC)
      !     T=TEMPERATURE (EV)
      !     ALOG(SVJ)=A(1)+A(2)*ALOG(T)+...+A(NJ)*ALOG(T)**(NJ-1)
      !              =A(1)+X*(A(2)+X*(A(3)+X*(...+X*(A(N-1)+A(NJ)*X)...)))

      !      DIMENSION A(NJ)
      !      REAL*8 TT

      implicit none
      integer nj,i

      real*8 a(*),t,svj,tt,sv,x
      TT = T
      X = LOG(TT)
      SV = A(NJ)
      DO  I = NJ-1, 1, -1
         SV = SV*X + A(I)
      end do
      SVJ = EXP(SV)
      RETURN

    END SUBROUTINE PNFIT


    REAL*8 FUNCTION XTFCT(X,IXTYPE)

      implicit none
      integer ixtype


      REAL*8 X,D,D1,D2,F,ER,SIGMA
      !COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      !real*8 tis,c1,c2,c3,c4
      !integer idqh

      !      GO TO ( 10 , 20 , 30 , 40 ), IDQH

      ER = TIS * X * X

      !      GO TO 50
      !   20 CALL HXHEP(ER,SIGMA)
      !      GO TO 50
      !   30 CALL XIPP(ER,SIGMA)
      !      GO TO 50
      !   40 CALL ELXS(KES,ER*AR,SIGMA)

10    CALL HXHP(ER,ixtype,SIGMA)
50    IF(C3.EQ.0.0) GO TO 60
      D = C4*X
      IF(D.GT.174.0) GO TO 60
      D = EXP(D)
      F = (D - 1.0D0/D)*C3
      GO TO 70
60    D1 = C2*(2.0D0*X-C2)
      D2 = -C2*(2.0D0*X+C2)
      F = EXP(D1) - EXP(D2)
70    XTFCT = X*X*F*SIGMA
      RETURN

    END FUNCTION XTFCT





end module mod_cx
