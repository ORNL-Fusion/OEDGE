C EIRENE07 COMPILATION
C ===== SOURCE: cdef.f
C
C
      SUBROUTINE CDEF(AL,JI,JE,K,COU,NTE,CF,LEXP,LTEST,LSUM,IFCLF)
C
C  EIRENE ATOMIC DATA , DEFAULT OR FROM FILE: POLYNOM FIT FORMAT
C
C  INPUT:
C    AL(J),J=1,NTE
C  OUTPUT:
C    MAXWELLIAN RATES, AL=LN(KT), KT IN (EV)
C  K
C       K>0: A&M DATA FROM FILES HYDHEL, METHANE OR AMJUEL
C           K:  NUMBER OF REACTION IN EIRENE "CREAC"-ARRAY
C
C       K<0: HARD WIRED EIRENE DEFAULT ATOMIC AND MOLECULAR DATA PACKAGE
C       EACH NUMBER ABS(K) CORRESPONDS TO ONE SPECIFIC REACTION DATA FIT
C       (SEE COMMENTS BELOW)
C       COEFFICIENTS FOR THESE DEFAULT REACTIONS ARE SPECIFIED IN
C       SUBROUTINE SETUP_DEFAULT_REACTIONS
C
C  JI,JE
C  FOR JI<=J<=JE RETURN:
C       1<=J<=9:  JTH ENERGY COEFFICIENT OF TWO PARAM. FITS
C                 AT TEMPERATUR KT (EV)
C                 (E.G.: RATES FOR HEAVY PARTICLE INTERACTIONS)
C                 IF ONLY ONE FIT AVAILABLE, ITS INDEX IS J=1
C                 E.Q. FOR ELECT. IMP. RATES AS FUNCTION OF TE
C  IFCLF 
C       FLAG DESCRIBING WHICH SORT OF REACTION DATA IS TO USED
C       = 1    POTENTIAL
C       = 2    CROSS SECTION
C       = 3    RATE COEFFICIENT
C       = 4    MOMENTUM-WEIGHTED RATE COEFFICIENT
C       = 5    ENERGY-WEIGHTES RATE COEFFICIENT
C       = 6    OTHER 
C
C  FIT FROM JANEV ET AL, PPPL-TM-368, 1985  (PREPRINT) OR:
C           SPRINGER SERIES ON ATOMS AND PLASMAS, VOL 4, 1987
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COMXS
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL(*)
      REAL(DP), INTENT(OUT) :: COU(0:9,*),CF(9,0:9)
      INTEGER, INTENT(IN) :: JI, JE, K, NTE, IFCLF
      LOGICAL, INTENT(IN) :: LEXP, LTEST, LSUM
      REAL(DP) :: CTEST
      INTEGER :: ICELL, J, II
C
C
C  DATA FROM A&M DATA FILES
      SELECT CASE (IFCLF)
      CASE (1)
        IF (.NOT.REACDAT(K)%LPOT) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR POTENTIAL',K
          CALL EXIT_OWN(1)
        END IF
        CF(1:9,1) = REACDAT(K)%POT%POLY%DBLPOL(1:9,1)
      CASE (2)
        IF (.NOT.REACDAT(K)%LCRS) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR CROSS SECTION',K
          CALL EXIT_OWN(1)
        END IF
        CF(1:9,1) = REACDAT(K)%CRS%POLY%DBLPOL(1:9,1)
      CASE (3)
        IF (.NOT.REACDAT(K)%LRTC) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR RATE COEFFICIENT',K
          CALL EXIT_OWN(1)
        END IF
        CF(1:9,JI:JE)=REACDAT(K)%RTC%POLY%DBLPOL(1:9,JI:JE)
      CASE (4)
        IF (.NOT.REACDAT(K)%LRTCMW) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' MOMENTUM-WEIGHTED RATE COEFFICIENT',K
          CALL EXIT_OWN(1)
        END IF
        CF(1:9,JI:JE)=REACDAT(K)%RTCMW%POLY%DBLPOL(1:9,JI:JE)
      CASE (5)
        IF (.NOT.REACDAT(K)%LRTCEW) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' ENERGY-WEIGHTED RATE COEFFICIENT',K
          CALL EXIT_OWN(1)
        END IF
        CF(1:9,JI:JE)=REACDAT(K)%RTCEW%POLY%DBLPOL(1:9,JI:JE)
      CASE (6)
        IF (.NOT.REACDAT(K)%LOTH) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' OTHER REACTION',K
          CALL EXIT_OWN(1)
        END IF
        CF(1:9,JI:JE)=REACDAT(K)%OTH%POLY%DBLPOL(1:9,JI:JE)
      CASE DEFAULT
        WRITE (IUNOUT,*) ' WRONG FLAG IFCLF SPECIFIED IN CDEF'
        WRITE (IUNOUT,*) ' IFCLF = ', IFCLF 
        CALL EXIT_OWN(1)
      END SELECT

      IF (LTEST) THEN
        DO 13 J=JI,JE
          CTEST=0.
          DO 14 II=1,9
            CTEST=CTEST+ABS(CF(II,J))
14        CONTINUE
          IF (CTEST.LE.EPS60) GOTO 990
13      CONTINUE
      ENDIF
C
      IF (LSUM) THEN
        DO J=JI,JE
          DO ICELL=1,NTE
            COU(J,ICELL)=CF(9,J)
          END DO
        END DO
C
        DO 21 J=JI,JE
          DO 22 II=8,1,-1
            DO 23 ICELL=1,NTE
              COU(J,ICELL)=COU(J,ICELL)*AL(ICELL)+CF(II,J)
23          CONTINUE
22        CONTINUE
21      CONTINUE
      ENDIF
C
      IF (LEXP) THEN
        DO 24 ICELL=1,NTE
          COU(JE,ICELL)=EXP(MAX(-100._DP,COU(JE,ICELL)))
24      CONTINUE
      ENDIF
C
      RETURN
C
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBROUTINE CDEF: ZERO FIT COEFFICIENTS'
      WRITE (iunout,*) 'J,K = ',J,K,'  EXIT CALLED!'
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: cdefn.f
C
C
      SUBROUTINE CDEFN(AL,PL,K,COU,NTE,CF,LEXP,LTEST,LSUM)
C
C  EIRENE DEFAULT ATOMIC DATA FOR INTERACTION WITH HYDROGEN
C  SAME AS CDEF, BUT FOR 2 PARAMETER FITTING EXPRESSIONS
C
C  INPUT:
C    AL(J),J=1,NTE, PL(I),I=1,NTE
C  OUTPUT:
C    MAXWELLIAN RATES, AL=LN(KT), KT IN (EV)
C                      PL=LN(NE), NE IN (1/CM**3)
C  K
C       K>0: A&M DATA FROM FILES HYDHEL, METHANE OR AMJUEL
C           K:  NUMBER OF REACTION IN EIRENE "CREAC"-ARRAY
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COMXS
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL(*), PL(*)
      REAL(DP), INTENT(OUT) :: COU(0:9,*), CF(9,0:9)
      INTEGER, INTENT(IN) :: K, NTE
      LOGICAL, INTENT(IN) :: LEXP, LTEST, LSUM
      REAL(DP) :: DUMMP(9)
      REAL(DP) :: CCXM1, CCXM2, EXPO1, EXPO2, FPAR1, FPAR2, FPAR3, S01,
     .          S02, DS12, CTEST, EXTRAP
      INTEGER :: I, JJ, KK, IFEX, J, II, ICELL
C
C  K>0:
C  DATA FROM ARRAY CREAC(9,0:9,K)
C
      IF (K.GT.0) THEN
C  DATA FROM A&M DATA FILES
        IF (.NOT.REACDAT(K)%LRTC) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR RATE COEFFICIENT',K
          CALL EXIT_OWN(1)
        END IF
        IF (REACDAT(K)%RTC%IFIT == 1) THEN
          WRITE (IUNOUT,*) ' ONLY 1D FIT AVAILABLE FOR RATE COEFF.',K
          CALL EXIT_OWN(1)
        END IF
        DO 11 J=1,9
          DO 12 II=1,9
            CF(II,J)=REACDAT(K)%RTC%POLY%DBLPOL(II,J)
12        CONTINUE
11      CONTINUE
        IF (LTEST) THEN
          DO 13 J=1,9
            CTEST=0.
            DO 14 II=1,9
              CTEST=CTEST+ABS(CF(II,J))
14          CONTINUE
            IF (CTEST.LE.EPS60) GOTO 990
13        CONTINUE
        ENDIF
      ELSE
        GOTO 990
      ENDIF
C
      IF (LSUM) THEN
        DO 20 ICELL=1,NTE
          COU(1,ICELL)=0.
20      CONTINUE
C
        DO 25 ICELL=1,NTE
          IF (AL(ICELL).LT.REACDAT(K)%RTC%POLY%RCMN) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(<S*V>)
            S01=REACDAT(K)%RTC%POLY%RCMN
            S02=LOG(2._DP)+REACDAT(K)%RTC%POLY%RCMN
            DS12=S02-S01
            EXPO1=0.
            EXPO2=0.
            DO 1 J=1,9
              JJ=J-1
              DO 1 I=1,9
                II=I-1
                EXPO1=EXPO1+S01**II*PL(ICELL)**JJ*CF(I,J)
                EXPO2=EXPO2+S02**II*PL(ICELL)**JJ*CF(I,J)
1           CONTINUE
            CCXM1=EXPO1
            CCXM2=EXPO2
            FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
            FPAR2=      (CCXM2-CCXM1)/DS12
            FPAR3=0.D0
C
            IFEX=5
            COU(1,ICELL)=EXTRAP(AL(ICELL),IFEX,FPAR1,FPAR2,FPAR3)
            if (.not.lexp) cou(1,icell)=log(cou(1,icell))
C
          ELSE
C
            DO 22 JJ=9,1,-1
              DUMMP(JJ)=CF(9,JJ)
              DO 23 KK=8,1,-1
                DUMMP(JJ)=DUMMP(JJ)*AL(ICELL)+CF(kk,jj)
23            CONTINUE
22          CONTINUE
            cou(1,icell)=dummp(9)
            DO 24 JJ=8,1,-1
              cou(1,icell)=cou(1,icell)*PL(icell)+DUMMP(JJ)
24          CONTINUE
C
C           DO 22 J=1,9
C             JJ=J-1
C             DO 22 I=1,9
C               II=I-1
C               COU(1,ICELL)=COU(1,ICELL)+AL(ICELL)**II*
C    .                                    PL(ICELL)**JJ*CF(I,J)
C22          CONTINUE
C
            if (lexp) COU(1,ICELL)=EXP(MAX(-100._DP,COU(1,ICELL)))
          ENDIF
25      CONTINUE
      ENDIF
C
C
      RETURN
C
C
990   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN SUBROUTINE CDEFN: ZERO FIT COEFFICIENTS'
      WRITE (iunout,*) 'J,K = ',J,K,'  EXIT CALLED!'
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: condense.f
C
C
      SUBROUTINE CONDENSE
C
C  CONDENSE COLLISION KERNEL, IF SOME SECONDARIES ARE NOT FOLLOWED
C  BY EIRENE, I.E., IF NFOLA(IATM), NFOLM(IMOL), NFOLI(IION) LT 0
C  FOR SOME TEST PARTICLE SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: ISP, IION, ICOL, IATM, IMOL, IRDS, IMDS

      DO 10 IATM=1,NATMI
C  NRCA=0 ?
        DO 100 ICOL=1,NRCA(IATM)
100     CONTINUE
10    CONTINUE
C
C
      DO 20 IMOL=1,NMOLI
C  currently: only electron impact collisions
        DO 200 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 220 IION=1,NIONI
            ISP=NSPAM+IION
            IF (PIODS(IRDS,IION).GT.0) THEN
              IF (NFOLI(IION).LT.0) THEN
                WRITE (iunout,*) 'TEST ION ',TEXTS(ISP),
     .                           ' CAN BE CONDENSED'
              ENDIF
            ENDIF
220       CONTINUE
200     CONTINUE
20    CONTINUE
C
      RETURN
      END
C ===== SOURCE: cross.f
C 0406: default resonant cx for He in He+/He++ plasma added:
C       Janev (HYDHEL) ,1987, reactions 5.3.1 and 6.3.1
C
      FUNCTION CROSS(AL,K,IR,TEXT)
C
C  CROSS SECTION
C    AL=LN(ELAB), ELAB IN (EV)
C    RETURN CROSS SECTION IN CM**2
C
C  K>0 :  DATA FROM ARRAY CREAC(9,0:9,K)
C
C  K<0 :  DEFAULT MODEL DEFINED IN SETUP_DEFAULT_REACTIONS
C
C  K=-1:  H + H+ --> H+ + H   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 3.1.8
C
C  K=-2:  He + He+ --> He+ + He   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 5.3.1
C
C  K=-3:  He + He++ --> H++ + He   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 6.3.1
C
      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL
      INTEGER, INTENT(IN) :: K, IR
      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      REAL(DP) :: B(8), FP(6)
      REAL(DP) :: S01, S02, DS12, EXPO1, EXPO2, CROSS,
     .            CCXM1, CCXM2, EXPO, EXTRAP, E, XI, SNGL_POLY
      INTEGER :: IF8, II, I
      type(poly_data), pointer :: rp

C
      IF ((K >= -10) .AND. (K <= NREAC)) THEN
        IF (IFTFLG(K,1) == 0) THEN

          RP => REACDAT(K)%CRS%POLY
          EXPO = SNGL_POLY(RP%DBLPOL,AL,RP%RCMN,RP%RCMX,RP%FPARM,
     .                     RP%IFEXMN,RP%IFEXMX)
          CROSS = EXP(MAX(-100._DP,EXPO))

        ELSE IF (IFTFLG(K,1) == 3) THEN
C  default extrapolation ifexmn=-1 not yet available
C  ELAB BELOW MINIMUM ENERGY FOR FIT:
          IF (AL.LT.REACDAT(K)%CRS%POLY%RCMN) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
            FP = REACDAT(K)%CRS%POLY%FPARM
            CROSS=EXTRAP(AL,REACDAT(K)%CRS%POLY%IFEXMN,
     .                   FP(1),FP(2),FP(3))
C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
          ELSEIF (AL.GT.REACDAT(K)%CRS%POLY%RCMX) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
            FP = REACDAT(K)%CRS%POLY%FPARM
            CROSS=EXTRAP(AL,REACDAT(K)%CRS%POLY%IFEXMX,
     .                   FP(4),FP(5),FP(6))
          ELSE
            E = EXP(AL)
            XI = REACDAT(K)%CRS%POLY%DBLPOL(1,1)
            B(1:8) = REACDAT(K)%CRS%POLY%DBLPOL(2:9,1)
            CROSS = B(1)*LOG(E/XI)
            DO I=1,7
              CROSS = CROSS + B(I+1)*(1.D0-XI/E)**I
            END DO
            CROSS = CROSS * 1.D-13/(XI*E)
          ENDIF
        ELSE
          WRITE (iunout,*) ' WRONG FITTING FLAG IN CROSS '
          WRITE (iunout,*) ' K = ',K,' IFTFLG = ',IFTFLG(K,1)
          WRITE (iunout,*) 'REACTION NO. ',IR
          CALL EXIT_OWN(1)
        END IF
      ELSE
        WRITE (iunout,*) 'ERROR IN CROSS: K= ',K
        WRITE (iunout,*) 'CALLED FROM ',TEXT
        WRITE (iunout,*) 'REACTION NO. ',IR
        CALL EXIT_OWN(1)
      ENDIF
      RETURN
      END
C ===== SOURCE: energy_rate_coeff.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only 
!pb  30.11.06: divide energy rate coefficient by ELCHA to get correct units

      function energy_rate_coeff (ir, p1, p2, lexp, iprshft) 
     .                           result (rate)

      use precision
      use parmmod
      use comxs
      use ccona
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp) :: rate, sngl_poly, dum(9), rcmin, rcmax, fp(6), q1, q2
      real(dp), save :: xlog10e, xln10, dsub, xlnelch
      integer :: jfexmn, jfexmx
      integer, save :: ifirst=0, ifsub=0 

      interface
        function intp_adas (ad,p1,p2) result(res)
          use precision
          use comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function intp_adas
      end interface

      if (.not.reacdat(ir)%lrtcew) then
        write (iunout,*) ' no data for energy weighted rate',
     .                   ' coefficient available for reaction ',ir
        call exit_own(1)
      end if

      rate = 0._dp

      if ((reacdat(ir)%rtcew%ifit == 1) .or. 
     .    (reacdat(ir)%rtcew%ifit == 2)) then
        rcmin  = reacdat(ir)%rtcew%poly%rcmn
        rcmax  = reacdat(ir)%rtcew%poly%rcmx
        fp     = reacdat(ir)%rtcew%poly%fparm
        jfexmn = reacdat(ir)%rtcew%poly%ifexmn
        jfexmx = reacdat(ir)%rtcew%poly%ifexmx
      end if

      if (mod(iftflg(ir,2),100) == 10) then
       
        rate = reacdat(ir)%rtcew%poly%dblpol(1,1)

      elseif (reacdat(ir)%rtcew%ifit == 1) then

        rate = sngl_poly(reacdat(ir)%rtcew%poly%dblpol(1:9,1),p1,
     .                   rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))

      else if (reacdat(ir)%rtcew%ifit == 2) then
         
        if (ifsub == 0) then
          ifsub = 1
          dsub = log(1.e8_dp)
        end if

        q2 = p2
        if (iprshft > 0) q2 = q2 - dsub
         
        call dbl_poly(reacdat(ir)%rtcew%poly%dblpol,p1,q2,rate,dum,1,9,
     .                rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
        
      else if (reacdat(ir)%rtcew%ifit == 3) then
         
        if (ifirst == 0) then
          ifirst = 1
          xln10 = log(10._dp)
          xlog10e = 1._dp/xln10
          xlnelch = log(elcha)
        end if

        q1 = xlog10e*p1
        q2 = xlog10e*p2
        rate = intp_adas(reacdat(ir)%rtcew%adas,q1,q2)

        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate
        end if
        rate = rate - xlnelch
        
      end if
      

      return

      end function energy_rate_coeff
     
C ===== SOURCE: extrap.f
C
C
      FUNCTION EXTRAP(ELAB,IFLAG,FP1,FP2,FP3)
C
C  NOTE:
C  INPUT:  ELAB IS LOG OF RELATIVE ENERGY, OR LOG OF TEMP
C  OUTPUT: EXTRAP IS NOT LOG, BUT THE TRUE VALUE
C
C  FUNCTION FOR EXTRAPOLATING SINGLE PARAMETER FITS BEYOND THEIR
C  RANGE OF VALIDITY
C  TYPE  IFLAG=1--4: JANEV ET AL. , SPRINGER, 1987, P13
C  TYPE  IFLAG=5  BACHMANN ET AL., IPP-REPORT, .....ELASTIC
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ELAB, FP1, FP2, FP3
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: X, EL, EXTRAP

      IF (IFLAG.EQ.1) THEN
C  NON ZERO THRESHOLD
        EXTRAP=0.
      ELSEIF (IFLAG.EQ.2) THEN
C  EXTRAPOLATION AT HIGH ENERGY END FOR REACTIONS WITH NON ZERO THRESHOLD
C  FP1 SHOULD BE = E_THRESHOLD (EV)
        EL=EXP(ELAB)
        X=EL/FP1
        EXTRAP=FP2*X**FP3*LOG(X)
      ELSEIF (IFLAG.EQ.3) THEN
        EXTRAP=EXP(FP1+FP2*ELAB)
      ELSEIF (IFLAG.EQ.4) THEN
C
C  OUT
C       EXTRAP=EXP((FP1+FP2*ELAB)**2)
C
      ELSEIF (IFLAG.EQ.5) THEN
C  LINEAR OR QUADRATIC EXTRAPOLATION IN LN(SIGMA)
        EXTRAP=EXP(FP1+FP2*ELAB+FP3*ELAB**2)
      ELSE
        GOTO 999
      ENDIF
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'ERROR IN EXTRAP. EXIT CALLED '
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: feelei1.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  30.11.06: DELPOT introduced

      FUNCTION FEELEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: ELEIC(9), FEELEI1, PLS, DEIMIN, EE, FTABEI1, DSUB,
     .            ELEI, ENERGY_RATE_COEFF, DELE
      INTEGER :: J, I, KK, II

      FEELEI1=0.D0
      KK=NELREI(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEELEI1=-EIONHE
        CASE (-2)
            FEELEI1=EELDS1(IREI,1)
        CASE (-3)
            FEELEI1=-1.5*TEIN(K)
        CASE (-4)
            FEELEI1=-EIONH
        CASE (-5)
            FEELEI1=-10.5
        CASE (-6)
            FEELEI1=-25.0
        CASE (-7)
            FEELEI1=EELDS1(IREI,1)
        CASE (-8)
            FEELEI1=-10.5
        CASE (-9)
            FEELEI1=-15.5
        CASE (-10)
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            FEELEI1=-0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        IF (JELREI(IREI) == 1) THEN
          ELEI = ENERGY_RATE_COEFF(KK,TEINL(K),0._DP,.TRUE.,0)
          FEELEI1=-ELEI*DEIN(K)*FACREA(KK)/(FTABEI1(IREI,K)+EPS60)
        ELSE
!pb          DSUB=LOG(1.D8)
          DEIMIN=LOG(1.D8)
!pb          PLS=MAX(DEIMIN,DEINL(K))-DSUB
          PLS=MAX(DEIMIN,DEINL(K))
          ELEI = ENERGY_RATE_COEFF(KK,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELEI+FACREA(KK)+DEINL(K))
          FEELEI1=-EXP(EE)/(FTABEI1(IREI,K)+EPS60)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          FEELEI1=FEELEI1+DELE
        END IF
      END IF

      RETURN
      END
C ===== SOURCE: feelrc1.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  19.12.06: bremsstrahlung added


      FUNCTION FEELRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: ELRC1(9), PLS, DELE, EE, FEELRC1, FTABRC1, DSUB,
     .            DEIMIN, ELRC, ENERGY_RATE_COEFF, BREMS, Z, ngffmh 
      INTEGER :: J, I, KK, II
      LOGICAL :: LADAS

      FEELRC1=0.D0
      KK=NELRRC(IRRC)

      IF (KK < 0) THEN

        SELECT CASE (KK)
        CASE (-1)
            FEELRC1=-1.5*TEIN(K)*FTABRC1(IRRC,K)
        CASE (-2)
            FEELRC1=EELRC1(IRRC,1)*FTABRC1(IRRC,K)
        CASE (-3)
            FEELRC1=-1.5*TEIN(K)*FTABRC1(IRRC,K)
        END SELECT

      ELSE IF (KK > 0) THEN

        IF (JELRRC(IRRC) == 1) THEN
          ELRC = ENERGY_RATE_COEFF(KK,TEINL(K),0._DP,.FALSE.,0)
          ELRC=ELRC+FACREA(KK)
          ELRC=EXP(MAX(-100._DP,ELRC))
          FEELRC1=-ELRC*DEIN(K)
        ELSE
!pb          DSUB=LOG(1.D8)
          DEIMIN=LOG(1.D8)
!pb          PLS=MAX(DEIMIN,DEINL(K))-DSUB
          PLS=MAX(DEIMIN,DEINL(K))
          ELRC = ENERGY_RATE_COEFF(KK,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELRC+DEINL(K)+FACREA(KK))
          FEELRC1=-EXP(EE)
        END IF

        LADAS = IS_RTCEW_ADAS(KK)
        IF (LADAS.AND.(NCHRGP(IPLS) /= 0)) THEN
          Z = NCHRGP(IPLS)
          BREMS = 1.54E-32_DP * TEIN(J)**0.5 * Z**2 *
     .            ngffmh(Z**2 * 13.6_DP/TEIN(J))*
     .            DEIN(K)*EXP(FACREA(KK))/ELCHA
          FEELRC1=FEELRC1 + BREMS
        END IF

        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          FEELRC1=FEELRC1+DELE*FTABRC1(IRRC,K)
        END IF
      END IF

      RETURN
      END
C ===== SOURCE: fehvds1.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FEHVDS1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: FEHVDS1, EHVDS, FTABEI1, RATE_COEFF
      INTEGER :: II, KK

      FEHVDS1=0.D0
      KK=NREAHV(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEHVDS1=EHVDS1(IREI,1)
        CASE (-2)
            FEHVDS1=6.
        CASE (-3)
            FEHVDS1=10.0
        CASE (-4)
            FEHVDS1=8.6
        CASE (-5)
            FEHVDS1=0.5
        CASE (-6)
            FEHVDS1=0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        EHVDS = RATE_COEFF(KK,TEINL(K),0._DP,.FALSE.,0)
        EHVDS=EXP(MAX(-100._DP,EHVDS+FACREA(KK)))
        FEHVDS1=EHVDS*DEIN(K)/(FTABEI1(IREI,K)+EPS60)
      END IF

      RETURN
      END
C ===== SOURCE: feplcx3.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FEPLCX3 (IRCX,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRCX, K
      REAL(DP) :: PLS, ADD, EPCX, FEPLCX3, RATE_COEFF
      INTEGER :: II, KK, IPLSTI

      FEPLCX3=0.D0
      KK=NELRCX(IRCX)
      IPLSTI=MPLSTI(IPLS)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
C  DEFAULT CX MODEL
            FEPLCX3=1.5*TIIN(IPLSTI,K)+EDRIFT(IPLS,K)
        CASE (-2)
C  MEAN ENERGY FROM DRIFTING MONOENERGETIC
            FEPLCX3=EPLCX3(IRCX,1,1)+EDRIFT(IPLS,K)
        CASE (-3)
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
            FEPLCX3=1.5*TIIN(IPLSTI,K)+EDRIFT(IPLS,K)
        END SELECT
      ELSE
C  MEAN ENERGY FROM SINGLE PARAMETER FIT KK
        PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
        EPCX = RATE_COEFF(KK,PLS,0._DP,.FALSE.,0)
        ADD=EPLCX3(IRCX,1,1)
        FEPLCX3=EPCX*DIIN(IPLS,K)*ADD
      END IF

      RETURN
      END
C ===== SOURCE: feplel3.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FEPLEL3 (IREL,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREL, K
      REAL(DP) :: PLS, ADD, EPEL, FEPLEL3, RATE_COEFF
      INTEGER :: II, KK, IPLSTI

      FEPLEL3=0.D0
      KK=NELREL(IREL)
      IPLSTI=MPLSTI(IPLS)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
C  DEFAULT EL MODEL
C  OUT
        CASE (-2)
C  MEAN ENERGY FROM DRIFTING MONOENERGETIC
            FEPLEL3=EPLEL3(IREL,1,1)+EDRIFT(IPLS,K)
        CASE (-3)
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
            FEPLEL3=1.5*TIIN(IPLSTI,K)+EDRIFT(IPLS,K)
        END SELECT
      ELSE
C  MEAN ENERGY FROM SINGLE PARAMETER FIT KK
        PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
        EPEL = RATE_COEFF(KK,PLS,0._DP,.FALSE.,0)
        ADD=EPLEL3(IREL,1,1)
        FEPLEL3=EPEL*DIIN(IPLS,K)*ADD
      END IF

      RETURN
      END
C ===== SOURCE: feplot3.f


      FUNCTION FEPLOT3 (IROT,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS
      USE PHOTON

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IROT, K
      REAL(DP) :: PLS, ADD, ELCX, FEPLOT3
      INTEGER :: II, KK, IPLSTI

      FEPLOT3=0.D0
      KK=NELROT(IROT)
      IPLSTI=MPLSTI(IPLS)
      SELECT CASE (KK)
      CASE (-2)
         FEPLOT3=EPLOT3(IROT,1,1)+EDRIFT(IPLS,K)
      CASE (-3)
         FEPLOT3=1.5*TIIN(IPLSTI,K)+EDRIFT(IPLS,K)
      CASE DEFAULT
        WRITE (iunout,*) ' FEPLOT3 NOT READY '
        WRITE (iunout,*) ' CALCULATION ABANDONNED '
        CALL EXIT_OWN(1)
      END SELECT

      RETURN
      END
C ===== SOURCE: feplpi3.f


      FUNCTION FEPLPI3 (IRPI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: FEPLPI3
      INTEGER :: KK, IPLSTI

      FEPLPI3=0.D0
      KK=NELRPI(IRPI)
      IPLSTI=MPLSTI(IPLS)
      SELECT CASE (KK)
      CASE (-1)
        FEPLPI3=EPLPI3(IRPI,1,1)
      CASE (-2)
        FEPLPI3=1.5*TIIN(IPLSTI,K)+EDRIFT(IPLS,K)
      END SELECT

      RETURN
      END
C ===== SOURCE: fi.f
C
C
      FUNCTION FI(R,ER,B,IFLAG,P,DFI)

C  EVALUATE EFFECTIVE POTENTIAL FUNCTION FI AT R
C  EVALUATE DFI(R)/DR AT R
C  RETURN FI=FI(R), DFI=DFI(R)/DR
C     --------------
C  IFLAG=1:  H+ + H
C  IFLAG=2:  H+ + NOBLE GASES,  H+ + H2,  HE+ + HE
C

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: R, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: DFI
      INTEGER, INTENT(IN) :: IFLAG

      REAL(DP) :: DSS, U, DU, SS, RI, RIQ, G1, REFF, RLIM, RR1, RR,
     .            EMRR, G2, EX2, B2, V, DV, EX, RLOW, R2, R3, FI, G,
     .            RMI, EPS

      B2=B*B
C
      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        R2=R*R
        R3=R2*R
C  FIND V=V(R) and DV=DV(R)/DR
        IF (R.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSEIF (R.LT.RLOW) THEN
          R2=RLOW*RLOW
          R3=R2*RLOW
          EX=EXP(-RLOW)
          EX2=EX*EX
          SS=(1.+RLOW+R2/3.)*EX-1.
          RI=1./RLOW
          V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          DV=0.D0
        ELSE
          EX=EXP(-R)
          EX2=EX*EX
          SS=(1.+R+R2/3.)*EX-1.
          RI=1./R
          V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          RIQ=RI*RI
          DSS=-R/3.*(1.+R)*EX
          U=RI-(1.+RI)*EX2-(1.+R)*EX
          DU=-RIQ+(RIQ+2.*RI+2.)*EX2+R*EX
          DV=27.211*((DU*SS-U*DSS)/(SS*SS)-RIQ)
        ENDIF
C  FIND FI=FI(R) AND DFI=DFI(R)/DR
        FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSEIF (IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        R2=R*R
        R3=R2*R
C
        RR=R*RMI
        EMRR=1.-RR
C
C       G2=1.00+(1.0-G2)*MAX(0.D0,-EMRR)/EMRR
C       REFF=-G1*G2*EMRR
        IF (RR.LT.1.0) THEN
          G=G1
        ELSE
          G=G1*G2
        ENDIF
        REFF=-G*EMRR
C
        IF (REFF.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSE
          EX=EXP(-REFF)
          EX2=EX*EX
          V=EPS*(EX2-(EX+EX))
          DV=-2.*EPS*RMI*G*(EX2-EX)
        ENDIF
        FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSE
        WRITE (iunout,*) 'ERROR IN FUNCTION FI. IFLAG INVALID.'
        WRITE (iunout,*) 'IFLAG = ',IFLAG
        CALL EXIT_OWN(1)
      ENDIF
C
      RETURN
      END
C ===== SOURCE: fivec.f
C
C
      SUBROUTINE FIVEC(ER,B,IFLAG,P)

C  VECTORIZED VERSION OF FUNCTION FI FOR GAUSS MEHLER QUADRATURE
C  EVALUATE EFFECTIVE POTENTIAL FUNCTION AT AR(I),I=1,NFI
C  NOTE: NFI.LE.128 IS NOT CHECKED, BUT USED
C  RETURN FI(AR(I)) IN THE ARRAY AFI(I),I=1,NFI
C     --------------
C  IFLAG=1:  H+ + H
C  IFLAG=2:  H+ + NOBLE GASES , H+ + H2, HE+ + HE
C

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: RI, G1, SS, EX, EX2, RR, RLIM, REFF, EMRR, G2, V,
     .          R, R2, B2, RLOW, RR1, G, RMI, EPS
      INTEGER :: IFI

      REAL(DP) :: AR(128), AFI(128)
      INTEGER :: NFI
      COMMON /CFI/ AR,AFI,NFI

      B2=B*B
C
      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        DO 1 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
C  FIND V=V(R)
          IF (R.GT.160.D0) THEN
            V=0.D0
          ELSEIF (R.LT.RLOW) THEN
            R2=RLOW*RLOW
            EX=EXP(-RLOW)
            EX2=EX*EX
            SS=(1.+RLOW+R2/3.)*EX-1.
            RI=1./RLOW
            V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          ELSE
            EX=EXP(-R)
            EX2=EX*EX
            SS=(1.+R+R2/3.)*EX-1.
            RI=1./R
            V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
1       CONTINUE
C
      ELSEIF(IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES, (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        DO 2 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
          RR=R*RMI
          EMRR=1.-RR
          IF (RR.LT.1) THEN
            G=G1
          ELSE
            G=G1*G2
          ENDIF
          REFF=-G*EMRR
C
          IF (REFF.GT.160.D0) THEN
            V=0.D0
          ELSE
            EX=EXP(-REFF)
            EX2=EX*EX
            V=EPS*(EX2-(EX+EX))
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
2       CONTINUE
C
      ELSE
        WRITE (iunout,*) 'ERROR IN FUNCTION FIVEC. IFLAG INVALID.'
        WRITE (iunout,*) 'IFLAG = ',IFLAG
        CALL EXIT_OWN(1)
      ENDIF
C
      RETURN
      END
C ===== SOURCE: fpatha.f
C   29.11.05:   comments,  empty lines, etc... to syncronize with
C               other fpath.. routines
C               use CSPEI removed
C               IRDS --> IREI
C               SIGMAX NOW SET ONLY FOR ACTIVE REACTIONS
C               added: jcou,ncou
!pb  30.08.06:  data structure for reaction data redefined
!pb  12.10.06:  modcol revised
!pb  22.11.06:  flag for shift of first parameter to rate_coeff introduced
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
C
      FUNCTION FPATHA (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR 
C   "BEAM OF NEUTRAL ATOMS" OF VELOCITY VEL
C   IN DRIFTING MAXWELLIAN PLASMA-BACKGROUND
C
C   INPUT:
C   K         :  CURRENT GRID CELL
C   JCOU, NCOU:  THERE WILL BE NCOU CALLS TO FPATH, FOR SAME TEST PARTICLE
C                COORDINATES. THIS CURRENT CALL IS CALL NO. JCOU.

C   OUTPUT: COMMON COMLCA
C           CFLAG: FLAG FOR SAMPLING OF POST COLLISION STATES
C           CFLAG(1,...): EI
C           CFLAG(3,...): CX
C           CFLAG(4,...): II
C           CFLAG(5,...): EL
C           CFLAG(6,...): RC
C
C  CFLAG(...,1):
C      =0:   VI: DELTA COLLISION IN VELOCITY SPACE (BUT DIFFERENT
C                                                   SPECIES ALLOWED)
C      =1:   VI: MONOENERGETIC AND ISOTROPIC IN CENTER OF MASS SYSTEM
C      =2:   VI: MAXWELL PLUS DRIFT
C      =3:   VI: SIGMA-V-WEIGHTED MAXWELLIAN PLUS DRIFT
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS
      USE CESTIM

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(7,3)
      INTEGER, INTENT(IN) :: K,JCOU,NCOU

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: TBPI3(9), TBCX3(9), TBEL3(9), FP(6)
      REAL(DP) :: EPPI3(9), EPCX3(9), EPEL3(9)
      REAL(DP) :: CEL, RMN, RLMS, ER, RMI, RMSI, CXS, VEFFQ, TBCX, VEFF,
     .          FEPLCX3, SIG, FEPLEL3, TBEL, RATE_COEFF, SNGL_POLY,
     .          ELTHDUM, CTCHDUM, SIGMAX, FTABEI1, EHEAVY, FEELEI1,
     .          FEHVDS1, DENEL, FPATHA, VX, VY, VZ, PVELQ0, ELAB,
     .          VRELQ, VREL, FEPLPI3, CII, CROSS, ELB, PLS, TBPI, EXPO,
     .          RCMIN, RCMAX, ENERGY_RATE_COEFF
      INTEGER :: IBGK, IRCX, IAEL, IREL, IPL, IAT, IAEI, IREI, IAPI,
     .           IRPI, KK, IREAC, IACX, II, IF8, JAN, J, I1, I2, IPLSTI,
     .           IPLSV, IROT
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb      IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb      ENDIF
      FPATHA=1.D10
      SIGMAX=0.D0
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      PVELQ0=VEL*VEL

      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
C
      DO 3 IPLS=1,NPLSV
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C  ELECTRON IMPACT COLLISION - RATE - COEFFICIENT
C  NO MASS SCALING NEEDED FOR BULK ELECTRONS
C
20    IF (LGAEI(IATM,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 30
      DO 10 IAEI=1,NAEII(IATM)
        IREI=LGAEI(IATM,IAEI)
        IF (MODCOL(1,2,IREI).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IREI)=TABDS1(IREI,K)
          ELSE
            SIGVEI(IREI)=FTABEI1(IREI,K)
          END IF
        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IREI,5)=EELDS1(IREI,K)
          EHEAVY=EHVDS1(IREI,K)
        ELSE
          ESIGEI(IREI,5)=FEELEI1(IREI,K)
          EHEAVY=FEHVDS1(IREI,K)
        ENDIF
C
        ESIGEI(IREI,1)=EATDS(IREI,0,1)*E0+EATDS(IREI,0,2)*EHEAVY
        ESIGEI(IREI,2)=EMLDS(IREI,0,1)*E0+EMLDS(IREI,0,2)*EHEAVY
        ESIGEI(IREI,3)=EIODS(IREI,0,1)*E0+EIODS(IREI,0,2)*EHEAVY
        ESIGEI(IREI,4)=EPLDS(IREI,  1)*E0+EPLDS(IREI,  2)*EHEAVY

        SIGMAX=MAX(SIGMAX,SIGVEI(IREI))
        SIGEIT=SIGEIT+SIGVEI(IREI)
10    CONTINUE
C
C  GENERAL ION IMPACT ON ATOM IATM, BULK ION SPEZIES IPLS=1,NPLSI
C  30--->40
C
30    IF (LGAPI(IATM,0,0).EQ.0) GOTO 40
      DO 36 IAPI=1,NAPII(IATM)
        IRPI=LGAPI(IATM,IAPI,0)
        IPLS=LGAPI(IATM,IAPI,1)
        IPLSV=MPLSV(IPLS)
        IPLSTI=MPLSTI(IPLS)
        IF (LGVAC(K,IPLS)) GOTO 36
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(4,2,IRPI).EQ.1) THEN
C  MAXWELL
          IF (NSTORDR >= NRAD) THEN
            SIGVPI(IRPI)=TABPI3(IRPI,K,1)
          ELSE
CDR  SIGVPI(IRPI)=FTABPI3 : NOT READY
            KK=NREAPI(IRPI)
            PLS=TIINL(IPLSTI,K)+ADDPI(IRPI,IPLS)
            TBPI = RATE_COEFF(KK,PLS,0._DP,.TRUE.,0)*DIIN(IPLS,K)
            SIGVPI(IRPI)=TBPI
          END IF
        ELSEIF (MODCOL(4,2,IRPI).EQ.2) THEN
C  BEAM - MAXWELL
C
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFPI(IRPI))
          IF (NSTORDR >= NRAD) THEN
            TBPI3(1:NSTORDT) = TABPI3(IRPI,K,1:NSTORDT)
            FP = 0._DP
            RCMIN = -HUGE(1._DP)
            RCMAX = HUGE(1._DP)
            EXPO = SNGL_POLY(TBPI3,ELB,RCMIN,RCMAX,FP,0,0)
          ELSE
! CALCULATE RATE-COEFFICIENT 
            KK=NREAPI(IRPI)
            PLS=TIINL(IPLSTI,K)+ADDPI(IRPI,IPLS)
            EXPO = RATE_COEFF(KK,PLS,ELB,.FALSE.,0) + DIINL(IPLS,K)
     .             + FACREA(KK)
          END IF
          SIGVPI(IRPI)=EXP(EXPO)
        ELSEIF (MODCOL(4,2,IRPI).EQ.3) THEN
C  BEAM - BEAM
          VRELQ=ZTI(IPLS)+PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFPI(IRPI)
          IREAC=MODCOL(4,1,IRPI)
          CII=CROSS(ELAB,IREAC,IRPI,'FPATHA II')
          SIGVPI(IRPI)=CII*VREL*DENIO(IPLS)
        ELSE
          GOTO 991
        ENDIF
C
C  2.A ELECTRON ENERGY LOSS / COLLISION (EV)
C
        ESIGPI(IRPI,2)=0.D0
C
        SIGMAX=MAX(SIGMAX,SIGVPI(IRPI))
        SIGPIT=SIGPIT+SIGVPI(IRPI)
C
C  2.B BULK ION ENERGY LOSS / COLLISION (EV)
C
        IF (NSTORDR >= NRAD) THEN
          ESIGPI(IRPI,1)=EPLPI3(IRPI,K,1)
        ELSE
          ESIGPI(IRPI,1)=FEPLPI3(IRPI,K)
        END IF
        CFLAG(4,1)=2
36    CONTINUE
C
C  CHARGE EXCHANGE RATE COEFFICIENT OF ATOMS IATM
C  WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
C  40--->50
C
40    CONTINUE
      IF (LGACX(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
      DO 41 IACX=1,NACXI(IATM)
        IRCX=LGACX(IATM,IACX,0)
        IPLS=LGACX(IATM,IACX,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF (LGVAC(K,IPLS)) GOTO 41
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(3,2,IRCX).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE, IGNORE NEUTRAL VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
CDR  SIGVCX(IRCX)=FTABCX3 : NOT READY
            KK=NREACX(IRCX)
            PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
            TBCX = RATE_COEFF(KK,PLS,0._DP,.TRUE.,0)*DIIN(IPLS,K)
            SIGVCX(IRCX)=TBCX
          END IF
        ELSEIF (MODCOL(3,2,IRCX).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
C     HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
C           HENCE: USE BEAM-BEAM RATE INSTEAD.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,IRCX)
            CXS=CROSS(ELAB,IREAC,IRCX,'FPATHA CX1')
            SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBCX3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT 
              KK=NREACX(IRCX)
              PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = RATE_COEFF(KK,PLS,ELB,.FALSE.,0) + DIINL(IPLS,K)
     .               + FACREA(KK)
            END IF
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,IRCX).EQ.3) THEN
C  MODEL 3:
C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLSV)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHA CX2')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSEIF (MODCOL(3,2,IRCX).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHA CX3')
          SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
C
        SIGMAX=MAX(SIGMAX,SIGVCX(IRCX))
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(3,4,IRCX).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (LEA.AND.(IESTCX(IRCX,3).EQ.0)) THEN  ! for tracklength estimator only
            IF (NSTORDR >= NRAD) THEN
              ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
            ELSE
              ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
            END IF
          END IF  ! this was for tracklength estimator only
          CFLAG(3,1)=2
        ELSEIF (MODCOL(3,4,IRCX).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (LEA.AND.(IESTCX(IRCX,3).EQ.0)) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
          IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
            EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
            FP = 0._DP
            RCMIN = -HUGE(1._DP)
            RCMAX = HUGE(1._DP)
            EXPO = SNGL_POLY(EPCX3,ELB,RCMIN,RCMAX,FP,0,0)
          ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT 
            KK=NELRCX(IRCX)
            PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
            EXPO = ENERGY_RATE_COEFF(KK,PLS,ELB,.FALSE.,0) 
     .             + DIINL(IPLS,K) + FACREA(KK)
          END IF
          ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
          ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(3,1)=3
        ELSEIF (MODCOL(3,4,IRCX).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (LEA.AND.(IESTCX(IRCX,3).EQ.0)) THEN  ! for tracklength estimator only
            IF (NSTORDR >= NRAD) THEN
              ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
            ELSE
              ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
            END IF
          END IF  ! this was for tracklength estimator only
          CFLAG(3,1)=1
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C  ELASTIC COLLISIONS OF ATOMS IATM  WITH IONS OF SPEZIES IPLS=1,NPLSI
C  50--->60
C
50    CONTINUE
      IF (LGAEL(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 60
      DO 51 IAEL=1,NAELI(IATM)
        IREL=LGAEL(IATM,IAEL,0)
        IPLS=LGAEL(IATM,IAEL,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IBGK=NPBGKP(IPLS,1)
        IF (LGVAC(K,IPLS)) GOTO 51
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(5,2,IREL).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE, IGNORE NEUTRAL VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVEL(IREL)=TABEL3(IREL,K,1)
          ELSE
            KK=NREAEL(IREL)
            PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
            TBEL = RATE_COEFF(KK,PLS,0._DP,.TRUE.,0)*DIIN(IPLS,K)
            SIGVEL(IREL)=TBEL
          END IF
        ELSEIF (MODCOL(5,2,IREL).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
C  TEMPERATURE TOO LOW, USE: BEAM_ATOM - BEAM_DRIFT RATECOEFF.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFEL(IREL)
            IREAC=MODCOL(5,1,IREL)
            CEL=CROSS(ELAB,IREAC,IREL,'FPATHA EL1')
            SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBEL3(1:NSTORDT) = TABEL3(IREL,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT 
              KK=NREAEL(IREL)
              PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
              EXPO = RATE_COEFF(KK,PLS,ELB,.FALSE.,0) + DIINL(IPLS,K)
     .               + FACREA(KK)
            END IF
            SIGVEL(IREL)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(5,2,IREL).EQ.3) THEN
C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLSV)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,IREL)
C  FIND SIGMA FROM OAK RIDGE "ELASTIC" DATA TABLES
          IF (LHABER) THEN
            RMN=RMASSA(IATM)
            RMI=RMASSP(IPLS)
            RMSI=1./(RMN+RMI)
            RLMS=RMN*RMI*RMSI
            ER=RLMS*VEFFQ*CVELI2
cdr  flag -1.0_DP: only sigma(ER), but no scattering angle evaluated
            CALL SCATANG (ER,-1.0_DP,ELTHDUM,CTCHDUM,SIG)
            CEL= SIG*AU_TO_CM2
          ELSE
C  FIND SIGMA FROM AMJUEL DATA TABLES (BACHMANN ET AL.)
            CEL=CROSS(ELAB,IREAC,IREL,'FPATHA EL2')
          END IF
          SIGVEL(IREL)=CEL*VEFF*DENIO(IPLS)
        ELSEIF (MODCOL(5,2,IREL).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,IREL)
          CEL=CROSS(ELAB,IREAC,IREL,'FPATHA EL3')
          SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
        ELSE
          GOTO 995
        ENDIF

        SIGMAX=MAX(SIGMAX,SIGVEL(IREL))
        SIGELT=SIGELT+SIGVEL(IREL)

        IF (IBGK.NE.0) SIGBGK=SIGBGK+SIGVEL(IREL)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(5,4,IREL).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=FEPLEL3(IREL,K)
          END IF
          CFLAG(5,1)=2
        ELSEIF (MODCOL(5,4,IREL).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (IESTEL(IREL,3).EQ.0) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
          IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
            EPEL3(1:NSTORDT) = EPLEL3(IREL,K,1:NSTORDT)
            FP = 0._DP
            RCMIN = -HUGE(1._DP)
            RCMAX = HUGE(1._DP)
            EXPO = SNGL_POLY(EPEL3,ELB,RCMIN,RCMAX,FP,0,0)
          ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT 
            KK=NELREL(IREL)
            PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
            EXPO = ENERGY_RATE_COEFF(KK,PLS,ELB,.FALSE.,0) 
     .             + DIINL(IPLS,K) + FACREA(KK)
          END IF
          ESIGEL(IREL,1)=EXP(EXPO)/SIGVEL(IREL)
          ESIGEL(IREL,1)=ESIGEL(IREL,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(5,1)=3
        ELSEIF (MODCOL(5,4,IREL).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=FEPLEL3(IREL,K)
          END IF
          CFLAG(5,1)=1
        ELSE
          GOTO 995
        ENDIF
51    CONTINUE
C
60    CONTINUE
C
C     TOTAL
C
100   CONTINUE
C
      IF (SIGEIT.GT.0._DP) THEN
        DO IREI=1,NRDS
          IF (SIGVEI(IREI) .LE. SIGMAX*1.D-10) THEN
            SIGEIT=SIGEIT-SIGVEI(IREI)
            SIGVEI(IREI) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGPIT.GT.0._DP) THEN
        DO IRPI=1,NRPI
          IF (SIGVPI(IRPI) .LE. SIGMAX*1.D-10) THEN
            SIGPIT=SIGPIT-SIGVPI(IRPI)
            SIGVPI(IRPI) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGCXT.GT.0._DP) THEN
        DO IRCX=1,NRCX
          IF (SIGVCX(IRCX) .LE. SIGMAX*1.D-10) THEN
            SIGCXT=SIGCXT-SIGVCX(IRCX)
            SIGVCX(IRCX) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGELT.GT.0._DP) THEN
        DO IREL=1,NREL
          IF (SIGVEL(IREL) .LE. SIGMAX*1.D-10) THEN
            SIGELT=SIGELT-SIGVEL(IREL)
            SIGVEL(IREL) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGOTT.GT.0._DP) THEN
        DO IROT=1,NROT
          IF (SIGVOT(IROT) .LE. SIGMAX*1.D-10) THEN
            SIGOTT=SIGOTT-SIGVOT(IROT)
            SIGVOT(IROT) = 0.D0
          END IF
        END DO
      END IF

C
      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT+SIGOTT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHA=VEL/SIGTOT
        ZMFPI=1./FPATHA
      ENDIF
C
      RETURN
990   CONTINUE
      WRITE (iunout,*) 'IATM,IREI,MODCOL(1,J,IREI) '
      WRITE (iunout,*) IATM,IREI,(MODCOL(1,J,IREI),J=1,4)
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN FPATHA: INCONSISTENT ION IMP. DATA'
      WRITE (iunout,*) 'IATM,IRPI,MODCOL(4,J,IRPI) '
      WRITE (iunout,*) IATM,IRPI,(MODCOL(4,J,IRPI),J=1,4)
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN FPATHA: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (iunout,*) 'IATM,IRCX,MODCOL(3,J=1,4,IRCX) '
      WRITE (iunout,*) IATM,IRCX,(MODCOL(3,J,IRCX),J=1,4)
      CALL EXIT_OWN(1)
995   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN FPATHA: INCONSISTENT ELASTIC COLL. DATA'
      WRITE (iunout,*) 'IATM,IREL,MODCOL(5,J,IREL) '
      WRITE (iunout,*) IATM,IREL,(MODCOL(5,J,IREL),J=1,4)
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: fpathi.f
c  25.11.05: option modcol(3,4...)=3 added 
c            (adopted from fpatha)
c            cx rate option 4 added (adopted from fpatha)
c            syncronized with fpatha
c  still missing: el (and bgk) and pi reactions
C               added: jcou,ncou
!pb  30.08.06:  data structure for reaction data redefined
!pb  12.10.06:  modcol revised
!pb  22.11.06:  flag for shift of first parameter to rate_coeff introduced
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
C
      FUNCTION FPATHI (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR
C   "BEAM TEST IONS" OF VELOCITY VEL 
C   IN DRIFTING MAXWELLIAN PLASMA-BACKGROUND
C
C   INPUT:
C   K         :  CURRENT GRID CELL
C   JCOU, NCOU:  THERE WILL BE NCOU CALLS TO FPATH, FOR SAME TEST PARTICLE
C                COORDINATES. THIS CURRENT CALL IS CALL NO. JCOU.
 
C   OUTPUT: COMMON COMLCA
C           CFLAG: FLAG FOR SAMPLING OF POST COLLISION STATES
C                  SAME AS IN FPATHA.F
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(7,3)
      INTEGER, INTENT(IN) :: K,JCOU,NCOU

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: TBPI3(9), TBCX3(9), TBEL3(9), FP(6)
      REAL(DP) :: EPPI3(9), EPCX3(9), EPEL3(9)
      REAL(DP) :: ELB, EXPO, FEPLCX3, ELAB, VEFF, CROSS, CXS, SIGMAX,
     .          VX, VY, VZ, PVELQ0, FPATHI, DENEL, FTABEI1, PLS,
     .          TBCX, VEFFQ, FEELEI1, EHEAVY, FEHVDS1, RATE_COEFF,
     .          VREL,VRELQ, ENERGY_RATE_COEFF, RCMIN, RCMAX, SNGL_POLY 
      INTEGER :: J, IF8, JAN, IREAC, IPL, IIO, IIEI, KK, II,
     .           IBGK, IRCX, IIEL, IICX, IIPI, I1, I2, IPLSTI, IPLSV,
     .           IREI, IRPI
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb      IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb      ENDIF
      FPATHI=1.D10
      SIGMAX=0.D0
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      PVELQ0=VEL*VEL

      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
C
      DO 3 IPLS=1,NPLSV
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C  ELECTRON IMPACT COLLISION - RATE - COEFFICIENT
C  NO MASS SCALING NEEDED FOR BULK ELECTRONS
C
20    IF (LGIEI(IION,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 30
      DO 10 IIEI=1,NIDSI(IION)
        IREI=LGIEI(IION,IIEI)
        IF (MODCOL(1,2,IREI).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IREI)=TABDS1(IREI,K)
          ELSE
            SIGVEI(IREI)=FTABEI1(IREI,K)
          END IF
        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IREI,5)=EELDS1(IREI,K)
          EHEAVY=EHVDS1(IREI,K)
        ELSE
          ESIGEI(IREI,5)=FEELEI1(IREI,K)
          EHEAVY=FEHVDS1(IREI,K)
        ENDIF
C
        ESIGEI(IREI,1)=EATDS(IREI,0,1)*E0+EATDS(IREI,0,2)*EHEAVY
        ESIGEI(IREI,2)=EMLDS(IREI,0,1)*E0+EMLDS(IREI,0,2)*EHEAVY
        ESIGEI(IREI,3)=EIODS(IREI,0,1)*E0+EIODS(IREI,0,2)*EHEAVY
        ESIGEI(IREI,4)=EPLDS(IREI,  1)*E0+EPLDS(IREI,  2)*EHEAVY

        SIGMAX=MAX(SIGMAX,SIGVEI(IREI))
        SIGEIT=SIGEIT+SIGVEI(IREI)
10    CONTINUE
C
C  GENERAL ION IMPACT ON TEST ION IION, BULK ION SPEZIES IPLS=1,NPLSI
C  30--->40
C
30    IF (LGIPI(IION,0,0).EQ.0) GOTO 40
      DO 36 IIPI=1,NIPII(IION)
        SIGVPI(IIPI)=0.
        SIGPIT=SIGPIT+SIGVPI(IIPI)
36    CONTINUE
C
C  CHARGE EXCHANGE RATE COEFFICIENT OF TEST ION IION
C  WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
C  40--->50
C
40    CONTINUE
      IF (LGICX(IION,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
      DO 41 IICX=1,NICXI(IION)
        IRCX=LGICX(IION,IICX,0)
        IPLS=LGICX(IION,IICX,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF (LGVAC(K,IPLS)) GOTO 41
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(3,2,IRCX).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE, IGNORE NEUTRAL VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
CDR  SIGVCX(IRCX)=FTABCX3 : NOT READY
            KK=NREACX(IRCX)
            PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
            TBCX = RATE_COEFF(KK,PLS,0._DP,.TRUE.,0)*DIIN(IPLS,K)
            SIGVCX(IRCX)=TBCX
          END IF
        ELSEIF (MODCOL(3,2,IRCX).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
C     HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
C           HENCE: USE BEAM-BEAM RATE INSTEAD.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,IRCX)
            CXS=CROSS(ELAB,IREAC,IRCX,'FPATHI CX1')
            SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBCX3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT 
              KK=NREACX(IRCX)
              PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = RATE_COEFF(KK,PLS,ELB,.FALSE.,0) + DIINL(IPLS,K)
     .               + FACREA(KK)
            END IF
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,IRCX).EQ.3) THEN
C  MODEL 3:
C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLSV)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHI CX2')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSEIF (MODCOL(3,2,IRCX).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHI CX3')
          SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
C
        SIGMAX=MAX(SIGMAX,SIGVCX(IRCX))
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(3,4,IRCX).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (IESTCX(IRCX,3).EQ.0) THEN  ! for tracklength estimator only
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=2
          ENDIF  ! this was for tracklength estimator only
        ELSEIF (MODCOL(3,4,IRCX).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (IESTCX(IRCX,3).EQ.0) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT 
              KK=NELRCX(IRCX)
              PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = ENERGY_RATE_COEFF(KK,PLS,ELB,.FALSE.,0) 
     .               + DIINL(IPLS,K) + FACREA(KK)
            END IF
            ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
            ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(3,1)=3
        ELSEIF (MODCOL(3,4,IRCX).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (IESTCX(IRCX,3).EQ.0) THEN  ! for tracklength estimator only
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          ENDIF  ! this was for tracklength estimator only
          CFLAG(3,1)=1
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C  ELASTIC COLLISIONS OF TEST-ION IION  WITH IONS OF SPEZIES IPLS=1,NPLSI 
C  50--->60
C
50    CONTINUE
C
60    CONTINUE
C
C     TOTAL
C
100   CONTINUE
C
      IF (SIGEIT.GT.0._DP) THEN
        DO IREI=1,NRDS
          IF (SIGVEI(IREI) .LE. SIGMAX*1.D-10) THEN
            SIGEIT=SIGEIT-SIGVEI(IREI)
            SIGVEI(IREI) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGPIT.GT.0._DP) THEN
        DO IRPI=1,NRPI
          IF (SIGVPI(IRPI) .LE. SIGMAX*1.D-10) THEN
            SIGPIT=SIGPIT-SIGVPI(IRPI)
            SIGVPI(IRPI) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGCXT.GT.0._DP) THEN
        DO IRCX=1,NRCX
          IF (SIGVCX(IRCX) .LE. SIGMAX*1.D-10) THEN
            SIGCXT=SIGCXT-SIGVCX(IRCX)
            SIGVCX(IRCX) = 0.D0
          END IF
        END DO
      END IF
C
      SIGTOT=SIGPIT+SIGCXT+SIGEIT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHI=VEL/SIGTOT
        ZMFPI=1./FPATHI
      ENDIF
C
      RETURN
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN FPATHI: INCONSISTENT ELEC. IMP. DATA'
      WRITE (iunout,*) 'IION,IREI,MODCOL(1,J=1,4,IREI) '
      WRITE (iunout,*) IION,IREI,(MODCOL(1,J,IREI),J=1,4)
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN FPATHI: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (iunout,*) 'IION,IRCX,MODCOL(3,J=1,4,IRCX) '
      DO 994 IPL=1,NPLSI
        WRITE (iunout,*) IION,IRCX,(MODCOL(3,J,IRCX),J=1,4)
994   CONTINUE
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: fpathm.f
!pb  30.08.06:  data structure for reaction data redefined
!pb  12.10.06:  modcol revised
!pb  22.11.06:  flag for shift of first parameter to rate_coeff introduced
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
C
C
      FUNCTION FPATHM (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR NEUTRAL
C   "BEAM MOLECULES" OF VELOCITY VEL IN DRIFTING MAXWELLIAN PLASMA-BACKGROUND
C   IN CELL K
C   CFLAG:  AS IN FUNCTION FPATHA
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS
      USE CTRCEI

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(7,3)
      INTEGER, INTENT(IN) :: K, JCOU,NCOU

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: TBPI3(9), TBCX3(9), TBEL3(9), FP(6)
      REAL(DP) :: EPPI3(9), EPCX3(9), EPEL3(9)
      REAL(DP) :: VREL, VRELQ, TBEL, FEPLCX3, CROSS, EXPO, ELB, CEL,
     .          SIGMAX, FEPLEL3, VX, VY, VZ, EHEAVY, FTABEI1, FPATHM,
     .          PVELQ0, DENEL, FEELEI1, VEFF, VEFFQ, CXS, ELAB,
     .          TBCX, FEHVDS1, PLS, RATE_COEFF, SNGL_POLY, RCMIN, RCMAX,
     .          ENERGY_RATE_COEFF
      INTEGER :: IBGK, IREL, IMEL, J, JAN, IF8, IML, IPL, IMDS, IRDS,
     .           II, IREAC, IMCX, IMPI, KK, IRCX, I1, I2, IPLSTI, IPLSV,
     .           IREI, IRPI
C
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb      IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb      ENDIF
      FPATHM=1.D10
      SIGMAX=0.D0
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSV
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C  ELECTRON IMPACT RATE COEFFICIENT
C
      IF (LGMEI(IMOL,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 25
      DO 10 IMDS=1,NMDSI(IMOL)
        IRDS=LGMEI(IMOL,IMDS)
        IF (MODCOL(1,2,IRDS).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IRDS)=TABDS1(IRDS,K)
          ELSE
            SIGVEI(IRDS)=FTABEI1(IRDS,K)
          END IF
        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IRDS,5)=EELDS1(IRDS,K)
          EHEAVY=EHVDS1(IRDS,K)
        ELSE
          ESIGEI(IRDS,5)=FEELEI1(IRDS,K)
          EHEAVY=FEHVDS1(IRDS,K)
        ENDIF
C
        ESIGEI(IRDS,1)=EATDS(IRDS,0,1)*E0+EATDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,2)=EMLDS(IRDS,0,1)*E0+EMLDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,3)=EIODS(IRDS,0,1)*E0+EIODS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,4)=EPLDS(IRDS,  1)*E0+EPLDS(IRDS,  2)*EHEAVY
C
        SIGMAX=MAX(SIGMAX,SIGVEI(IRDS))
        SIGEIT=SIGEIT+SIGVEI(IRDS)
10    CONTINUE
C
C
C  IONIZATION OF MOLECULE BY ION IMPACT, ION SPEZIES IPLS=1,NPLSI
C
25    CONTINUE
      DO 30 IMPI=1,NMPII(IMOL)
        SIGVPI(IMPI)=0.
        SIGPIT=SIGPIT+SIGVPI(IMPI)
30    CONTINUE
C
C  CHARGE EXCHANGE RATE COEFFICIENT FOR MOLECULE IMOL
C  WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
C  40--->50
C
40    CONTINUE
      IF (LGMCX(IMOL,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
      DO 41 IMCX=1,NMCXI(IMOL)
        IRCX=LGMCX(IMOL,IMCX,0)
        IPLS=LGMCX(IMOL,IMCX,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF (LGVAC(K,IPLS)) GOTO 41
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(3,2,IRCX).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE, IGNORE NEUTRAL VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
CDR  SIGVCX(IRCX)=FTABCX3 : NOT READY
            KK=NREACX(IRCX)
            PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
            TBCX = RATE_COEFF(KK,PLS,0._DP,.TRUE.,0)*DIIN(IPLS,K)
            SIGVCX(IRCX)=TBCX
          END IF
        ELSEIF (MODCOL(3,2,IRCX).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
C     HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
C           HENCE: USE BEAM-BEAM RATE INSTEAD.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,IRCX)
            CXS=CROSS(ELAB,IREAC,IRCX,'FPATHM CX')
            SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBCX3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT 
              KK=NREACX(IRCX)
              PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = RATE_COEFF(KK,PLS,ELB,.FALSE.,0) + DIINL(IPLS,K)
     .               + FACREA(KK)
            END IF
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,IRCX).EQ.3) THEN
C  MODEL 3:
C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLSV)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHM CX')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
C
        SIGMAX=MAX(SIGMAX,SIGVCX(IRCX))
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(3,4,IRCX).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=2
        ELSEIF (MODCOL(3,4,IRCX).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (IESTCX(IRCX,3).EQ.0) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT 
              KK=NELRCX(IRCX)
              PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = ENERGY_RATE_COEFF(KK,PLS,ELB,.FALSE.,0) 
     .               + DIINL(IPLS,K) + FACREA(KK)
            END IF
            ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
            ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(3,1)=3
        ELSEIF (MODCOL(3,4,IRCX).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=1
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C  ELASTIC COLLISIONS OF MOLECULE IMOL WITH IONS OF SPEZIES IPLS=1,NPLSI
C  50--->60
C
50    CONTINUE
      IF (LGMEL(IMOL,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 60
      DO 51 IMEL=1,NMELI(IMOL)
        IREL=LGMEL(IMOL,IMEL,0)
        IPLS=LGMEL(IMOL,IMEL,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IBGK=NPBGKP(IPLS,1)
        IF (LGVAC(K,IPLS)) GOTO 51
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(5,2,IREL).EQ.1) THEN
C  MAXWELL
          IF (NSTORDR >= NRAD) THEN
            SIGVEL(IREL)=TABEL3(IREL,K,1)
          ELSE
            KK=NREAEL(IREL)
            PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
            TBEL = RATE_COEFF(KK,PLS,0._DP,.TRUE.,0)*DIIN(IPLS,K)
            SIGVEL(IREL)=TBEL
          END IF
        ELSEIF (MODCOL(5,2,IREL).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFEL(IREL)
            IREAC=MODCOL(5,1,IREL)
            CEL=CROSS(ELAB,IREAC,IREL,'FPATHM EL')
            SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBEL3(1:NSTORDT) = TABEL3(IREL,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(TBEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT 
              KK=NREAEL(IREL)
              PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
              EXPO = RATE_COEFF(KK,PLS,ELB,.FALSE.,0) + DIINL(IPLS,K)
     .               + FACREA(KK)
            END IF
            SIGVEL(IREL)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(5,2,IREL).EQ.3) THEN
C  BEAM - BEAM
          VRELQ=ZTI(IPLS)+PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,IREL)
          CEL=CROSS(ELAB,IREAC,IREL,'FPATHM EL')
          SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
        ELSE
          GOTO 995
        ENDIF

        SIGMAX=MAX(SIGMAX,SIGVEL(IREL))
        SIGELT=SIGELT+SIGVEL(IREL)
        IF (IBGK.NE.0) SIGBGK=SIGBGK+SIGVEL(IREL)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(5,4,IREL).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=FEPLEL3(IREL,K)
          END IF
          CFLAG(5,1)=2
        ELSEIF (MODCOL(5,4,IREL).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (IESTEL(IREL,3).EQ.0) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY 
! PRECALCULATING TEMPERATURE DEPENDENCIES
              EPEL3(1:NSTORDT) = EPLEL3(IREL,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = SNGL_POLY(EPEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT 
              KK=NELREL(IREL)
              PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
              EXPO = ENERGY_RATE_COEFF(KK,PLS,ELB,.FALSE.,0)
     .               + DIINL(IPLS,K) + FACREA(KK)
            END IF
            ESIGEL(IREL,1)=EXP(EXPO)/SIGVEL(IREL)
            ESIGEL(IREL,1)=ESIGEL(IREL,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(5,1)=3
C       ELSEIF (MODCOL(5,4,IREL).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  NOT READY
        ELSE
          GOTO 995
        ENDIF
51    CONTINUE
C
60    CONTINUE
C
C     TOTAL
C
      IF (SIGEIT.GT.0._DP) THEN
        DO IREI=1,NRDS
          IF (SIGVEI(IREI) .LE. SIGMAX*1.D-10) THEN
            SIGEIT=SIGEIT-SIGVEI(IREI)
            SIGVEI(IREI) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGPIT.GT.0._DP) THEN
        DO IRPI=1,NRPI
          IF (SIGVPI(IRPI) .LE. SIGMAX*1.D-10) THEN
            SIGPIT=SIGPIT-SIGVPI(IRPI)
            SIGVPI(IRPI) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGCXT.GT.0._DP) THEN
        DO IRCX=1,NRCX
          IF (SIGVCX(IRCX) .LE. SIGMAX*1.D-10) THEN
            SIGCXT=SIGCXT-SIGVCX(IRCX)
            SIGVCX(IRCX) = 0.D0
          END IF
        END DO
      END IF

      IF (SIGELT.GT.0._DP) THEN
        DO IREL=1,NREL
          IF (SIGVEL(IREL) .LE. SIGMAX*1.D-10) THEN
            SIGELT=SIGELT-SIGVEL(IREL)
            SIGVEL(IREL) = 0.D0
          END IF
        END DO
      END IF
C
100   CONTINUE
      SIGTOT=SIGPIT+SIGCXT+SIGEIT+SIGELT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHM=VEL/SIGTOT
        ZMFPI=1./FPATHM
      ENDIF
C
      RETURN
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN FPATHM: INCONSISTENT ELEC. IMP. DATA'
      WRITE (iunout,*) 'IMOL,IRDS,MODCOL(1,J,IRDS) '
      WRITE (iunout,*) IMOL,IRDS,(MODCOL(1,J,IRDS),J=1,4)
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN FPATHM: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (iunout,*) 'IMOL,IRCX,(MODCOL(3,J,IRCX),J=1,4) '
      WRITE (iunout,*) IMOL,IRCX,(MODCOL(3,J,IRCX),J=1,4)
      CALL EXIT_OWN(1)
995   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN FPATHM: INCONSISTENT ELASTIC COLL. DATA'
      WRITE (iunout,*) 'IMOL,IREL,(MODCOL(5,J,IREL),J=1,4) '
      WRITE (iunout,*) IMOL,IREL,(MODCOL(5,J,IREL),J=1,4)
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: fpathph.f
C  27.6.05: iadd removed
C  02.01.06:    SIGMAX NOW SET ONLY FOR ACTIVE REACTIONS
C               added: jcou,ncou
!pb  12.10.06:  modcol revised
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
C
      FUNCTION FPATHPH (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR PHOTON
C   "BEAM" OF VELOCITY (E0,VEL_X,Y,Z) IN DRIFTING MAXWELLIAN BACKGROUND
C   IN CELL K
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS
      USE CSPEI
      USE PHOTON

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(7,3)
      INTEGER, INTENT(IN) :: K, JCOU,NCOU

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: FPATHPH, sigmax, sigv, feplot3, vx, vy, vz,
     .            DENEL, PVELQ0, fac
      integer :: il, kk, ipl, irot, ipot, j, iph, i1, i2
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb      IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb      ENDIF
      FPATHPH = 1.E10_DP
      SIGMAX=0.D0

      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
C
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSV
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
csw
csw OT processes (photonic reactions)
csw
60    CONTINUE
      if(phv_lgphot(iphot,0,0) == 0) goto 70
      do 61 ipot=1,phv_nphoti(iphot)
        irot=phv_lgphot(iphot,ipot,0)
        ipls =phv_lgphot(iphot,ipot,1)
        il   =phv_lgphot(iphot,ipot,2)
        kk   =phv_lgphot(iphot,ipot,3)
        IF (LGVAC(K,IPLS)) GOTO 61
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(7,2,   IROT).EQ.1) THEN
          GOTO 997
        ELSEIF (MODCOL(7,2,   IROT).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE. FULL ACCOUNT FOR DOPPLER SHIFT
cdr       kk   = nreaot(irot)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr  getcoeff liefert nun maxw. average ueber Ti(ipls), z.b. voigt, ....
          call PH_GETCOEFF(kk,iphot,0,k,ipls,fac,sigv)
          sigv=sigv*diin(ipls,k)
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpat
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVOT(irot)=sigv
          GOTO 997
        ELSEIF (MODCOL(7,2,   IROT).EQ.4) THEN
C  MODEL 4:
C  BEAM-BEAM RATE. IGNORE DOPPLER SHIFT DUE TO THERMAL MOTION,
C                  INCLUDE DOPPLER SHIFT DUE TO DIRECTED MOTION
cdr       kk   = nreaot(irot)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr       ireac=modcol(7,1,irot)
cdr  ireac entspricht "typ" in Getcoeff - cross section (lorentz, vdw, ...)
cdr  allerdings kann hier der "querschnitt" von hintergrundparametern abhaengen
          call PH_GETCOEFF(kk,iphot,0,k,ipls,fac,sigv)
          sigv=sigv*diin(ipls,k)
cdr  besser: sig vc. E0_effective, ohne v = vrel = c, dann
cdr  dann:   sigv=sig * c * diin
cdr  denn:   sigv enthaelt hier keine faltung ueber background maxw. vel-verteil
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpat
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVOT(irot)=sigv
cdr       ESIGOT(irot,1)=e0*sigv   ziemlich sicher falsch
        ELSE
          GOTO 997
        ENDIF

        SIGMAX=MAX(SIGMAX,SIGVOT(IROT))
        SIGOTT=SIGOTT+SIGVOT(IROT)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(7,4,     IROT).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGOT(IROT,1)=EPLOT3(IROT,K,1)
          ELSE
            ESIGOT(IROT,1)=FEPLOT3(IROT,K)
          END IF
          CFLAG(7,1)=2
        ELSEIF (MODCOL(7,4,     IROT).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGOT(IROT,1)=EPLOT3(IROT,K,1)
          ELSE
            ESIGOT(IROT,1)=FEPLOT3(IROT,K)
          END IF
          CFLAG(7,1)=1
        ELSE
          GOTO 997
        ENDIF
61    CONTINUE

70    CONTINUE
c
C     TOTAL
C
100   CONTINUE

C
      IF (SIGOTT.GT.0._DP) THEN
        DO IROT=1,NROT
          IF (SIGVOT(IROT) .LE. SIGMAX*1.D-10) THEN
            SIGOTT=SIGOTT-SIGVOT(IROT)
            SIGVOT(IROT) = 0.D0
          END IF
        END DO
      END IF

      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT+SIGOTT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHPH=VEL/SIGTOT
        ZMFPI=1./FPATHPH
      ENDIF
C
      RETURN
997   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN FPATHPH: INCONSISTENT PHOTON COLL. DATA.'
      WRITE (iunout,*) 'IPHOT,IROT,MODCOL(7,J,IROT) '
      WRITE (iunout,*) IPHOT,IROT,(MODCOL(7,J,IROT),J=1,4)
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: ftabei1.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FTABEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: TBEIC(9), FTABEI1, DEIMIN, DSUB, PLS, TBEI, RATE_COEFF
      INTEGER :: J, I, II, KK

      TBEI=0.D0
      KK = NREAEI(IREI)

!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
!pb      PLS=MAX(DEIMIN,DEINL(K))-DSUB
      PLS=MAX(DEIMIN,DEINL(K))

      TBEI = RATE_COEFF(KK,TEINL(K),PLS,.TRUE.,1)*FACREA(KK)
      IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)

      FTABEI1 = TBEI

      RETURN
      END

C ===== SOURCE: ftabrc1.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FTABRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: DEIMIN, DSUB, FTABRC1, ZX, TBRC, PL, RATE_COEFF
      INTEGER :: II, KK

      TBRC=0.D0
      KK = NREARC(IRRC)

      IF (KK == 0) THEN
        ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
        TBRC=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(K)

      ELSE

!pb        DSUB=LOG(1.D8)
        DEIMIN=LOG(1.D8)
!pb        PL=MAX(DEIMIN,DEINL(K))-DSUB
        PL=MAX(DEIMIN,DEINL(K))

        TBRC = RATE_COEFF(KK,TEINL(K),PL,.TRUE.,1)*FACREA(KK)
        IF (IFTFLG(KK,2) < 100) TBRC=TBRC*DEIN(K)
        
      ENDIF

      FTABRC1 = TBRC

      RETURN
      END

C ===== SOURCE: gaumeh.f
C
C
      SUBROUTINE GAUMEH(RS,ER,B,IFLAG,P,N,IGAUS,RESULT)
C
C  IGAUS=1
C  GAUSS MEHLER QUADRATURE, N=5,10,20, W(X)=1./SQRT(1+X)/SQRT(1-X)
C                           INTEGRATION FROM A=-1 TO B=1
C                           TRANSFORMED TO A=0,B=1, AND
C                           W(X)=1./SQRT(X)/SQRT(1-X)
C

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: RS, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: RESULT
      INTEGER, INTENT(IN) :: IFLAG, N, IGAUS

      REAL(DP) :: AR(128), AFI(128)
      INTEGER :: NFI
      COMMON /CFI/ AR,AFI,NFI
      REAL(DP) :: X5A(5),X5B(5),X10A(10),X10B(10),X20A(20),X20B(20)
      REAL(DP) :: X5(5),X10(10),X20(20)
C
      REAL(DP) :: XG10(10),WG10(10),XG10A(10),XG10B(10)
      REAL(DP) :: SUM, F, X, W5, W10, W20
      INTEGER :: I, IFIRST, IFI
C
      SAVE
      DATA X5 /0.95105654E+00,0.58778542E+00,0.31391647E-06,
     .       - 0.58778363E+00,-0.95105600E+00/,
     .W5 /0.62831837E+00/
      DATA X10/9.876884E-01, 8.910065E-01,7.071068E-01,4.539905E-01,
     .         1.564344E-01,-1.564345E-01,-4.539905E-01,-7.071068E-01,
     .        -8.910065E-01,-9.876884E-01/,
     .W10/3.141593E-01/
      DATA X20/9.969173E-01,9.723699E-01,9.238795E-01,8.526402E-01,
     .         7.604060E-01,6.494480E-01,5.224985E-01,3.826834E-01,
     .         2.334453E-01,
     .         7.845905E-02,-7.845914E-02,-2.334454E-01,-3.826835E-01,
     .        -5.224986E-01,-6.494481E-01,-7.604060E-01,-8.526402E-01,
     .        -9.238796E-01,-9.723700E-01,-9.969174E-01/,
     .W20/1.570796E-01/
C
C  IGAUS=2
C  GAUSS MEHLER QUADRATURE, N=10, W(X)=1./SQRT(1-X)
C                           INTEGRATION FROM A=0 TO B=1
      DATA XG10/
     .    .013695585480651072296,
     .    .070758123420104485088,
     .    .16782834791297652828,
     .    .29588270759990962383,
     .    .44298868539955669668,
     .    .59543571523425259033,
     .    .73901490631777356542,
     .    .86034375925702257262,
     .    .94811360601967719799,
     .    .99414369156320387287/
      DATA WG10/
     .    .035228014278304236584,
     .    .081202859600773882686,
     .    .12534409666821812627,
     .    .16655348315340949082,
     .    .20386023963448093937,
     .    .23638906392303662412,
     .    .26337727689835336360,
     .    .28419221863676482994,
     .    .29834597294520589473,
     .    .30550677426145261144/
      DATA IFIRST /0/
      IF (IFIRST.EQ.0) THEN
C SET ROOTS AND WEIGHTS FOR GAUSS QUADRATURE RULES
        IFIRST=1
        DO 5 I=1,5
          X=0.5*(1.+X5(I))
          X5A(I)=SQRT(X*(1.-X))
          X5B(I)=1./X
5       CONTINUE
        DO 10 I=1,10
          X=0.5*(1.+X10(I))
          X10A(I)=SQRT(X*(1.-X))
          X10B(I)=1./X
10      CONTINUE
        DO 11 I=1,10
          X=XG10(I)
          XG10A(I)=SQRT(1.-X)
          XG10B(I)=1./X
11      CONTINUE
        DO 20 I=1,20
          X=0.5*(1.+X20(I))
          X20A(I)=SQRT(X*(1.-X))
          X20B(I)=1./X
20      CONTINUE
      ENDIF
C
      NFI=N
      IF (N.EQ.5) THEN
        DO 50 IFI=1,5
          AR(IFI)=RS*X5B(IFI)
50      CONTINUE
C
        CALL FIVEC(ER,B,IFLAG,P)
        SUM=0.D0
        DO 51 IFI=1,5
          F=AFI(IFI)
          IF (F.GT.0.D0)
     .    SUM=SUM+X5A(IFI)/SQRT(F)
51      CONTINUE
        RESULT=W5*SUM
C
      ELSEIF (N.EQ.10) THEN
        IF (IGAUS.EQ.1) THEN
          DO 100 IFI=1,10
            AR(IFI)=RS*X10B(IFI)
100       CONTINUE
C
          CALL FIVEC(ER,B,IFLAG,P)
          SUM=0.D0
          DO 101 IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+X10A(IFI)/SQRT(F)
101       CONTINUE
          RESULT=W10*SUM
        ELSEIF (IGAUS.EQ.2) THEN
          DO 110 IFI=1,10
            AR(IFI)=RS*XG10B(IFI)
110       CONTINUE
C
          CALL FIVEC(ER,B,IFLAG,P)
          SUM=0.D0
          DO 111 IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+XG10A(IFI)/SQRT(F)*WG10(IFI)
111       CONTINUE
          RESULT=SUM
        ENDIF
C
      ELSEIF (N.EQ.20) THEN
        DO 200 IFI=1,20
          AR(IFI)=RS*X20B(IFI)
200     CONTINUE
C
        CALL FIVEC(ER,B,IFLAG,P)
        SUM=0.D0
        DO 201 IFI=1,20
          F=AFI(IFI)
          IF (F.GT.0.D0)
     .    SUM=SUM+X20A(IFI)/SQRT(F)
201     CONTINUE
        RESULT=W20*SUM
C
      ELSE
        WRITE (iunout,*) 'ERROR IN GAUMEH, WRONG PARAMETER N'
        CALL EXIT_OWN(1)
      ENDIF
      RETURN
      END
C ===== SOURCE: h_colrad.f
C*
C*     COLLISIONAL-RADIATIVE MODEL OF
C*
C*     ATOMIC HYDROGEN
C*
C*
C   ASSUME: SLOW SPECIES: H,H+
C   OUTPUT:
C   R0(..) : TRAIN OF H* TRAVELING WITH H+
C   R1(..) : TRAIN OF H* TRAVELING WITH H
C   R2(..) : TRAIN OF H* TRAVELING WITH external source, e.g. radiation trap.
C
C
c   C: stoesse nach oben
c   F: stoesse nach unten
c   A: spontan nach unten
c   S: ionisation
c
c   ALPHA(N): THREEBODY recombination from H+
c   BETA(N) : radiative rec. from H+
c   C(1,N)  : excitation from ground state
C   Q2(N)   : ???   ->  H*(N)  (external source)
C
c   reduced pop coeff r0,r1  are per electron. hence: times "densel"
c   for pop0,pop1 - arrays of reduced population coefficients
C*
C***********************************************************************
      SUBROUTINE H_COLRAD (TEMP, DENSEL, Q2, POP0, POP1, POP2,
     .                     ALPCR, SCR, SCRRAD,
     .                     E_ALPCR, E_SCR, E_SCRRAD,
     .                     E_ALPCR_T, E_SCR_T,E_SCRRAD_T)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

C--------- ATOMIC PARAMETER ------------------------------------------
      REAL(DP), INTENT(IN) :: TEMP, DENSEL
      REAL(DP), INTENT(IN) :: Q2(40)
      REAL(DP), INTENT(OUT) ::   ALPCR,    SCR,     SCRRAD
      REAL(DP), INTENT(OUT) :: E_ALPCR,  E_SCR,   E_SCRRAD
      REAL(DP), INTENT(OUT) :: E_ALPCR_T,E_SCR_T, E_SCRRAD_T
      REAL(DP), INTENT(OUT) :: POP0(40), POP1(40), POP2(40)

      REAL(DP), SAVE :: A(40,40), E_AT(40), OSC(40,40)
      REAL(DP), SAVE :: A21SAVE, POP_ESC

      REAL(DP) :: C(40,40),S(40),F(40,40)
     &           ,SAHA(40),BETA(40),ALPHA(40)
      REAL(DP) :: R1(40),R0(40),R2(40)
C
cdr   integer, parameter :: lupa=34, lima=40
      integer, parameter :: lupa=34, lima=34
      INTEGER, SAVE :: IFRST=0
      INTEGER :: IP
c
      IF (LIMA.GT.40.OR.LUPA.GT.LIMA) THEN
        WRITE (iunout,*) 'LIMA, LUPA ??? ',LIMA,LUPA
        CALL EXIT_OWN(1)
      ENDIF
C
C ATOM
c  pop_esc   : lyman alpha population escape factor
c  pop_esc= 1: lyman alpha opt. thin
c  pop_esc= 0: lyman alpha opt. thick
c
c  only once and for all !!
c
      IF (IFRST == 0) THEN
        POP_ESC=1.
        CALL EINSTN(OSC,A,E_AT,40,A21SAVE,POP_ESC)
        IFRST = 1
      END IF

      E_ALPCR_T=0._DP
c
c
C ATOM
      CALL CLSAHA(TEMP,SAHA)
      CALL RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA)

C
C***********************************************************************
C
C ATOMIC HYDROGEN
C
      CALL POPCOF(DENSEL,SAHA,C,F,S,A,ALPHA,BETA,LUPA,LIMA,R0,R1,
     &            R2,Q2)

C TRAINS OF ELECTRONICALLY EXCITED H
      DO IP=2,LIMA
C
CDR COUPLING TO H+ IONS
        POP0(IP)=R0(IP)*DENSEL
CDR COUPLING TO H ATOMS
        POP1(IP)=R1(IP)*DENSEL
CDR COUPLING TO EXTERNAL SOURCE Q
        POP2(IP)=R2(IP)

      END DO
C  EFFECTIVE IONISATION RATES
      CALL IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUPA,LIMA,
     &              ALPCR,SCR,F,R2,Q2,SCRRAD)
C***********************************************************************
C  EFFECTIVE ELECTRON COOLING RATES
      CALL E_IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUPA,LIMA,
     &              ALPCR,SCR,F,R2,Q2,SCRRAD,
     &              E_ALPCR,E_SCR,E_SCRRAD,E_AT,
     &              E_SCR_T,E_SCRRAD_T)
C

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      RETURN
      END SUBROUTINE H_COLRAD

C**********************************************************************
      SUBROUTINE EINSTN(F,A,E_AT,LIM,A21SAVE,POP_ESC)
C
C     CALCULATION OF OSCILLATOR STRENGTH AND EINSTEIN COEFFICIENT
C     FOR ATOMIC HYDROGEN
C     E_AT are the energy levels. Ground state is E_AT(1)=0.
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(40,40)
      DIMENSION A(40,40)
      DIMENSION E_AT(40)

      UH=13.595
      DO 100 I=1,LIM
        P=I
        E_AT(I)=UH*(1.0-1.0/P**2)
100   CONTINUE

      DO 101 I=1,LIM-1
      DO 102 J=I+1,LIM
      AI=I
      AJ=J
      X=1.-(AI/AJ)**2
      G=(.9935+.2328/AI-.1296/AI**2)
     *-(.6282-.5598/AI+.5299/AI**2)/(AI*X)
     *+(.3387-1.181/AI+1.470/AI**2)/(AI*X)**2
      IF(I.NE.1.AND.I.NE.2) GO TO 300
      G=1.0785-.2319/X+.02947/X**2
  200 IF(I.NE.1) GO TO 300
      G=1.1330-.4059/X+.07014/X**2
  300 F(I,J)=2.**6/(3.*SQRT(3.)*3.1416)*(AI/AJ)**3/(2.*AI**2)*G/X**3
      A(J,I)=8.03E9*AI**2/AJ**2*(AI**(-2)-AJ**(-2))**2*F(I,J)
  102 CONTINUE
  101 CONTINUE
cdr  nur a(2-->1) rausnehmen
      A21SAVE=A(2,1)
      A(2,1)=A(2,1)*pop_esc
      RETURN
      END

C***********************************************************************
      SUBROUTINE CLSAHA(TEMP,SAHA)
C
C     SAHA-BOLTZMANN COEFFICIENT FOR ATOMIC HYDROGEN
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SAHA(40)

      TE=TEMP*1.1605E4

      DO 101 I=1,40
      P=I
      UION=13.595/TEMP/P**2
c     if (uion.gt.65) then
c       saha(i)=1.e30
c     else
        SAHA(I)=P**2*EXP(UION)/2.414D15/SQRT(TE**3)
c     endif
  101 CONTINUE
      RETURN
      END

C***********************************************************************
      SUBROUTINE RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA)
C
C     RATE COEFFICIENT FOR ATOMIC HYDROGEN
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION OSC(40,40),C(40,40),F(40,40),U(40,40)
      DIMENSION SAHA(40),S(40),ALPHA(40),BETA(40),UION(40)

C     INITIALIZATION
      DO 1 I=1,40
      S(I)=0.0
      ALPHA(I)=0.0
      BETA(I)=0.0
      UION(I)=0.0
      DO 1 J=1,40
      C(I,J)=0.0
      F(I,J)=0.0
      U(I,J)=0.0
    1 CONTINUE

      TE=TEMP*1.1605E4
      UH=13.595

      DO 101 I=1,40
      P=I
  101 UION(I)=13.595/TEMP/P**2
      DO 102 I=1,40
      DO 102 J=1,40
  102 U(I,J)=UION(I)-UION(J)

      IF(TE.GT.5.0E3) THEN
        CALL EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)

        DO 105 I=1,40
  105     ALPHA(I)=S(I)*SAHA(I)

      ELSE
        CALL EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)
        S(1)=ALPHA(1)/SAHA(1)
        DO 19 I=2,40
   19     ALPHA(I)=S(I)*SAHA(I)
      END IF
c
      DO 602 I=1,40
      P=I
      XP=UH/TEMP/P**2
      CALL CLBETA(XP,P,XS)
  602 BETA(I)=5.197D-14*(UH/TEMP)**.5/P*XS
c

      RETURN
      END

C***********************************************************************
      SUBROUTINE EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)
C
C     EXCITATION RATE COEFFICIENT
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(40,40),OSC(40,40),C(40,40),F(40,40)
      DIMENSION S(40),ALPHA(40)
C
      DO 1 I=1,40
      S(I)=0.0
      ALPHA(I)=0.0
      DO 1 J=1,40
      C(I,J)=0.0
    1 F(I,J)=0.0

      TE=TEMP*1.1605E4

C*********  1 -> J
      I=1
      P=I
      DO 100 J=2,40
      Q=J
      CALL COF1N(U(I,J),OSC(I,J),TE,F1,I,J)
      F(J,1)=F1
  100 C(1,J)=Q**2/P**2/EXP(U(I,J))*F(J,I)

C*********  2-10 -> J
      DO 110 I=2,10
      P=I
      DO 110 J=I+1,40
      Q=J
      CALL COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
      CALL COFJO(U(I,J),OSC(I,J),TE,CJ,I,J)

      GG=((P-2.)/8.)**0.25
      C(I,J)=(1.-GG)*CJ+GG*CV
110   F(J,I)=P**2/Q**2*EXP(U(I,J))*C(I,J)


C*********  I(>11) -> J

      DO 120 I=11,39
      P=I
      DO 120 J=I+1,40
      Q=J
      CALL COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
      C(I,J)=CV
  120 F(J,I)=P**2/Q**2*EXP(U(I,J))*C(I,J)

C*********  S  1 ->
      I=1

      IF(TE.GT.5.0E3) THEN
      CALL COFJS(TE,S1,I)
      S(1)=S1

      ELSE
      CALL COFJS2(TE,AL,I)
      ALPHA(1)=AL
      END IF
C
C*********  S  2-10 ->
      DO 210 I=2,10
      P=I

      CALL COFJS(TE,SJ,I)
      CALL COFVS(TEMP,SV,I)

      GGG=((P-2.)/8.)**0.25
  210 S(I)=(1.-GGG)*SJ+GGG*SV

C*********  S  I(>11) ->
      DO 220 I=11,40
      CALL COFVS(TEMP,SV,I)
  220 S(I)=SV
      RETURN
      END

C***********************************************************************
      SUBROUTINE COF1N(U,OSC,TE,F,I,J)
C
C     C(1,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN 1 -> J
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION GINT
      EXTERNAL GINT

      P=I
      Q=J
      X=1-P**2/Q**2
      Y=U
C
      R=0.45*X
      IF(J.EQ.2) THEN
      RX=-0.95
      ELSE IF(J.EQ.3) THEN
      RX=-0.95
      ELSE IF(J.EQ.4) THEN
      RX=-0.95
      ELSE IF(J.EQ.5) THEN
      RX=-0.95
      ELSE IF(J.EQ.6) THEN
      RX=-0.94
      ELSE IF(J.EQ.7) THEN
      RX=-0.93
      ELSE IF(J.EQ.8) THEN
      RX=-0.92
      ELSE IF(J.EQ.9) THEN
      RX=-0.91
      ELSE
      RX=-0.9
      END IF
C
      Z=R+Y
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X-0.603/X**2)
      C1=2*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z

      IF(Y.GT.50.0) THEN
      E1Y=EXP(-Y)*GINT(Y)/Y
      E1Z=EXP(-Z)*GINT(Z)/Z
      E2Y=EXP(-Y)*(1.0-GINT(Y))
      E2Z=EXP(-Z)*(1.0-GINT(Z))
      E1=(1/Y+0.5)*E1Y+RX*(1/Z+0.5)*E1Z
      E2=E2Y/Y+RX*E2Z/Z
      F=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)*P**2/Q**2*EXP(Y)
      ELSE

      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)

      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y+RX*(1/Z+0.5)*E1Z
      E2=E2Y/Y+RX*E2Z/Z

      F=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)*P**2/Q**2*EXP(Y)
      END IF
      RETURN
      END

C***********************************************************************
      SUBROUTINE COFJO(U,OSC,TE,C,I,J)
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
C
C
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT REAL*8(A-H,O-Z)
      P=I
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      Q=J
      X=1-P**2/Q**2
      Y=U
      R=1.94/P**1.57*X
      Z=R+Y
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      C1=2*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)
      IF (ICON.NE.0) GOTO 1000
      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y-(1/Z+0.5)*E1Z
      E2=E2Y/Y-E2Z/Z

      C=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)

      RETURN
 1000 WRITE(iunout,*) 'ERROR IN COFJO        ICON = ',ICON
      STOP
      END

C***********************************************************************
      SUBROUTINE COFVR(U,OSC,TEMP,C,I,J)
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.VRIENS, A.H.M.SMEETS, PHYS. REV. A 22, 940 (1980).
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      UH=13.595
      P=I
      BN=1.4/P*LOG(P)-0.7/P-0.51/P**2+1.16/P**3-0.55/P**4
      Q=J
      S1=Q-P
      X=1-P**2/Q**2
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      DELTA=EXP(-B/A)+0.06*S1**2/P**2/Q
      G1=1+TEMP/UH*P**3
      G2=3+11*S1**2/P**2
      G3=6+1.6*S1*Q+0.3/S1**2+0.8*Q**1.5/S1**0.5*(S1-0.6)
      GAMMA=UH*LOG(G1)*G2/G3
      C1=1.6D-7*TEMP**0.5/(TEMP+GAMMA)*EXP(-U)
      C2=0.3*TEMP/UH+DELTA

      C=C1*(A*LOG(C2)+B)

      RETURN
      END

C***********************************************************************
      SUBROUTINE COFJS(TE,S,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(0:2,40)

      G(0,1)=1.1330
      G(1,1)=-0.4059
      G(2,1)=0.07014
      G(0,2)=1.0785
      G(1,2)=-0.2319
      G(2,2)=0.02947
      DO 350 N=3,40
      G(0,N)=0.9935+0.2328/N-0.1296/N**2
      G(1,N)=-0.6282/N+0.5598/N**2-0.5299/N**3
  350 G(2,N)=0.3887/N**2-1.181/N**3+1.470/N**4

      IF (I.EQ.1) THEN

      P=1.0

      Y=1.57770E5/TE
      R=0.45
      Z=R+Y
      RX=-0.59
      A=0.0
      DO 223 K=0,2
  223 A=A+G(K,1)/(K+3)
      A=A*1.9603*P

      B=0.66667*P**2*(5-0.603)
      C1=2*P**2
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)

      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=1/Y*E1Y+RX*1/Z*E1Z
      E0Y=EXP(-Y)/Y
      E0Z=EXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY+RX*EGZ

      S=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)

      ELSE
      P=I
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      Y=1.57770E5/TE/P**2
      R=1.94/P**1.57
      Z=R+Y
      A=0.0
      DO 222 K=0,2
  222 A=A+G(K,I)/(K+3)
      A=A*1.9603*P
      B=0.66667*P**2*(5+BN)
      C1=2*P**2
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)
      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=1/Y*E1Y-1/Z*E1Z
      E0Y=EXP(-Y)/Y
      E0Z=EXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY-EGZ

      S=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)
      END IF
      RETURN
      END

C***********************************************************************
      SUBROUTINE COFJS2(TE,W,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C     FOR LOW TEMPERATURE
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(0:2,40)
      DOUBLE PRECISION GINT
      EXTERNAL GINT

      G(0,1)= 1.1330
      G(1,1)=-0.4059
      G(2,1)= 0.07014
      G(0,2)= 1.0785
      G(1,2)=-0.2319
      G(2,2)= 0.02947
      DO 350 N=3,40
      G(0,N)=0.9935+0.2328/N-0.1296/N**2
      G(1,N)=-0.6282/N+0.5598/N**2-0.5299/N**3
  350 G(2,N)=0.3887/N**2-1.181/N**3+1.470/N**4

C     IF (I.EQ.1) THEN
      PP=I
      P=1.0
      Y=1.57770E5/TE
      R=0.45
      Z=R+Y
      RX=-0.59

      A=0.0
      DO 223 K=0,2
  223 A=A+G(K,1)/(K+3)
      A=A*1.9603*P

      B=0.66667*P**2*(5-0.603)
      C1=2*P**2
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z

      E1Y=EXP(-Y)*GINT(Y)/Y
      E1Z=EXP(-Z)*GINT(Z)/Z
      E2Y=EXP(-Y)*(1.0-GINT(Y))
      E2Z=EXP(-Z)*(1.0-GINT(Z))
      E1=1/Y*E1Y+RX*1/Z*E1Z
      E0Y=EXP(-Y)/Y
      E0Z=EXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY+RX*EGZ

      W=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)*EXP(13.595/TE*1.1605E4)
     */2.414D15/SQRT(TE**3)

      RETURN
      END

C***********************************************************************
      SUBROUTINE COFVS(TEMP,S,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C
C
C     L.VRIENS, A.H.M.SMEETS, PHYS. REV. A 22, 940 (1980).
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)

      P=I
      UI=13.595/TEMP/P**2
      UIZ=UI**2.33+4.38*UI**1.72+1.32*UI
      S=9.56D-6/TEMP**1.5*EXP(-UI)/UIZ

      RETURN
      END

cdr
      double precision FUNCTION GINT(xX)
      REAL*8 Xx
      double precision X,GG(8)
cdr
      DATA GG/0.2677737343,8.6347608925,18.0590169730,8.5733287401,
     *        3.9584960228,21.0996530827,25.6329561486,9.5733223454/
cdr
      x=dble(xx)
cdr
      GINT=(GG(1)+GG(2)*X+GG(3)*X**2+GG(4)*X**3+X**4)/(GG(5)+GG(6)
     *      *X+GG(7)*X**2+GG(8)*X**3+X**4)
      RETURN
      END


C***********************************************************************
      SUBROUTINE CLBETA(XP,P,S)
C
C     RADIATIVE RECOMBINATION RATE COEFFICIENT
C
C
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 GAUNT3,PP,XPP,A,B,EPSA,EPSR
      COMMON PP,XPP
      EXTERNAL GAUNT3

      II=INT(P)
      PP=P
      XPP=XP

      A=0.0
      B=20.0
cdr   EPSA=1.0D-5
      EPSA=1.0D-4
      EPSR=1.0D-5
      NMIN=15
      NMAX=511

      CALL  AQC8(A,B,GAUNT3,EPSA,EPSR,NMIN,NMAX,S,ERR,N,ICON)
C     IF (ICON.NE.0) GOTO 1000

C     WRITE (iunout,10) II,S,ERR,N,ICON
C  10 FORMAT(1H ,1I3,2(2X,1PD10.3),2X,I4,2X,I5,/)

      RETURN
C1000 WRITE(iunout,*) 'ERROR IN CLBETA       ICON = ',ICON
C     STOP
      END

C***********************************************************************

      REAL*8 FUNCTION GAUNT3(X)
      REAL*8 PP,XPP,U,B,X
      COMMON PP,XPP
      U=X/XPP
      B=PP
      GAUNT3=(1./(U+1.)+0.1728*(U-1.)/B**(2./3.)/(U+1.)**(5./3.)-0.0496*
     *(U**2+4./3.*U+1.)/B**(4./3.)/(U+1.)**(7./3.))*EXP(-X)
      RETURN
      END

C***********************************************************************
      SUBROUTINE POPCOF(DENSEL,SAHA,C,F,S,A,ALPHA,BETA,LUP,LIM,R0,R1,
     &      R2,Q2)
C
C     SOLUTION OF RATE EQUATION FOR ATOMIC HYDROGEN
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8      C(40,40),F(40,40),A(40,40),W(40,40)
     &         ,SAHA(40),S(40),ALPHA(40),BETA(40),R0(40),R1(40)
     &         ,       Q2(40),R2(40)
     &         ,BLAX(40),VW(40),WA(40,40)
      dimension ip(40)

      DO 201 K=2,LUP-1

        DO 202 L=2,K
cdr stoss bevoelkerung von k von unten
  202     W(K,L)=C(L,K)*DENSEL

cc diagonale
cc entvoelkerung durch stoesse nach unten
        SUMF=0.
        DO 301 I=1,K-1
  301     SUMF=SUMF+F(K,I)
cc entvoelkerung durch stoesse nach oben
        SUMC=0.
        DO 302 I=K+1,LIM
  302     SUMC=SUMC+C(K,I)
cc  spontan nach unten
        SUMA=0.
        DO 303 I=1,K-1
  303     SUMA=SUMA+A(K,I)
cdr entvoelkerung von k: stoesse nach unten, nach oben, ionis, spontan
cdr                      nach unten
        W(K,K)=-(DENSEL*(SUMF+SUMC+S(K))+SUMA)

cc diagonale fertig

        DO 203 L=K+1,LUP
cdr bevoelkerung durch: stoesse von oben, spontan von oben
  203     W(K,L)=DENSEL*F(L,K)+A(L,K)

  201 CONTINUE

cdr k loop finished
cdr: jetzt: dito fuer k=lup, d.h. bevoelkerung von oben entfaellt
      DO 211 L=2,LUP-1

  211 W(LUP,L)=C(L,LUP)*DENSEL

      SUMF=0.
      DO 311 I=1,LUP-1
  311 SUMF=SUMF+F(LUP,I)
      SUMC=0.0
      DO 313 I=LUP+1,LIM
  313 SUMC=SUMC+C(LUP,I)
      SUMA=0.
      DO 312 I=1,LUP-1
  312 SUMA=SUMA+A(LUP,I)

      W(LUP,LUP)=-(DENSEL*(SUMF+SUMC+S(LUP))+SUMA)

C  RECHTE SEITEN:

c  vorbereiten fuer recombination
      DO 550 K=2,LUP
        SUMF=0.0
        DO 500 I=LUP+1,LIM
  500     SUMF=SUMF+F(I,K)*SAHA(I)
        SUMAS=0.0
        DO 501 I=LUP+1,LIM
  501     SUMAS=SUMAS+SAHA(I)*A(I,K)
c
c  matrixelemente: 1/s  (densel*rate-coeff. )
c  rechte seiten : cm**3/s, nicht: 1/s, also fuer elektronendichte=1
c                                      (bzw: stosspartnerdichte =1)
c  geht wg. linearitaet.
c  recombination e + H+ --> H*
        W(K,LUP+1)=-(DENSEL*SUMF+SUMAS+(DENSEL*ALPHA(K)+BETA(K)))
c  ionisation e + H --> H*
        W(K,LUP+2)=-C(1,K)
c  external source: Q2
        W(K,LUP+3)=-Q2(K)
  550 CONTINUE

cdr w besetzt fuer w(i,j) i=2,lup,j=2,lup+3
cdr geht gut, solange lup<38
      DO 3001 II=LUP,LUP+2
cdr reduziere w indices um 1: auf wa: i=1,lup-1,j=1,(lup-1)+3
        DO 402 I=1,LUP-1
        DO 402 J=1,LUP-1+3
  402     WA(I,J)=W(I+1,J+1)
        DO 3000 J=1,LUP-1
          BLAX(J)=WA(J,II)
 3000   CONTINUE

        CALL LAX(WA,40,LUP-1,  BLAX,0.0,1,IS,VW,IP,ICON)
c

        DO 3010 J=1,LUP-1
          IF(II.EQ.LUP) THEN
coupling to H+
            R0(J+1)=BLAX(J)
          ELSE IF(II.EQ.LUP+1) THEN
coupling to H-groundstate
            R1(J+1)=BLAX(J)
          ELSE IF(II.EQ.LUP+2) THEN
coupling to Q2
            R2(J+1)=BLAX(J)
          END IF
 3010   CONTINUE
 3001 CONTINUE


      RETURN
      END

C***********************************************************************
      SUBROUTINE IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUP,LIM,
     &ALPCR,SCR,F,R2,Q2,SCRRAD)
C
C     EFFECTIVE IONIZATION AND RECOMBINATION RATE COEFFICIENTS
C     FOR ATOMIC HYDROGEN
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(40,40),S(40),SAHA(40),A(40,40),ALPHA(40),BETA(40),
     &          R0(40),R1(40),R2(40),Q2(40),F(40,40)
C
      DO 5000 I=LUP+1,LIM
        R0(I)=1.0*SAHA(I)
        R1(I)=0.0
        R2(I)=0.0
 5000 CONTINUE
C
      SCR=S(1)
      DO 5001 I=2,LUP
        SUSCR=C(1,I)-R1(I)*(F(I,1)*DENSEL+A(I,1))
 5001 SCR=SCR+SUSCR

      SCRRAD=0.00
      DO 5002 I=2,LUP
        SUSRAD=Q2(I)-R2(I)*(F(I,1)*DENSEL+A(I,1))
 5002 SCRRAD=SCRRAD+SUSRAD

      ALPCR1=DENSEL*ALPHA(1)+BETA(1)

      ALPCR2=0.0
      DO 5003 I=2,LIM

 5003 ALPCR2=ALPCR2+R0(I)*(DENSEL*F(I,1)+A(I,1))

      ALPCR=ALPCR1+ALPCR2
      RETURN
      END

C***********************************************************************
      SUBROUTINE E_IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUP,LIM,
     &                    ALPCR,SCR,F,R2,Q2,SCRRAD,
     &                    E_ALPCR,E_SCR,E_SCRRAD,E_AT,
     &                    E_SCR_T,E_SCRRAD_T)
C
C     EFFECTIVE ELECTRON ENERGY LOSS IONIZATION AND RECOMBINATION
C     RATE COEFFICIENTS FOR ATOMIC HYDROGEN
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(40,40),S(40),SAHA(40),A(40,40),ALPHA(40),BETA(40),
     &          R0(40),R1(40),R2(40),Q2(40),F(40,40),E_AT(40)
C
      DO 5000 I=LUP+1,LIM
        R0(I)=1.0*SAHA(I)
        R1(I)=0.0
        R2(I)=0.0
 5000 CONTINUE
C
      E_SCR=0.0
      E_SCRRAD=0.0
      E_ALPCR=0.0
C
C  FOR TESTING
      E_SCR_T=0.0
      E_SCRRAD_T=0.0
C
C  EFFECTIV ELECTRON COOLING CORRESPONDING TO
C  "ORDINARY" COUPLING TO GROUND STATE S(I)
C
C  PART I
C  1--> inf.  R1(1)=1
      UH=13.595
      DE=(E_AT(1)-UH)
      E_SCR  =S(1)*DE
C  PART II
C  i1--> inf.  R1(i1)=...
      do i1=2,lup
        DE=(E_AT(i1)-UH)
        E_SCR  =E_SCR+r1(i1)*densel*S(I1)*DE
      enddo
C  PART III
C  1--> I2.  R1(1)=1
      DO I2=2,LUP
        DE=(E_AT(1)-E_AT(I2))
        SUSCR=C(1,I2)*DE
        E_SCR=E_SCR+SUSCR
      ENDDO
C  PART IV
C  i1--> i2. i1<i2,  R1(i1)=...
      DO  I1=2,LUP-1
        DO  I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=r1(i1)*densel*C(I1,I2)*DE
          E_SCR=E_SCR+SUSCR
        ENDDO
      ENDDO
C  PART V
C  i2--> i1.  i2>i1, R1(i2)=...
C            (includes i1=1, i.e., invers to PART III)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))
          SUSCR  =r1(i2)*densel*F(I2,I1)*DE
          E_SCR=E_SCR+SUSCR
C  separte treatment of radiation losses alone, only for testing
          SUSCR_T=r1(i2)*A(I2,I1)*(-1.)*DE
          E_SCR_T=E_SCR_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only.  evaluate second formular for E_SCR,
C                  using radiation loss E_SCR_T and effective rate SCR
      UH=13.595
      DE=(E_AT(1)-UH)
      E_SCR_T=E_SCR_T+SCR *DE

C  contribution for "ordinary" coupling to ground state done

c  next: contribution for coupling to external source Q2
c        use same codes as for E_SCR, but with Q2(K) <-- C(1,K)
c                                     and S(1)=0
c        and, of course, r2(k) instead of r1(k)
C        furthermore: no electron energy losses associated with Q2: E_Q2=0
C
C
C  PART I
C  zero, because S(1)=0
C  PART II
C  i1--> inf.  R1(i1)=...
      do i1=2,lup
        DE=(E_AT(i1)-UH)
        E_SCRRAD  =E_SCRRAD+r2(i1)*densel*S(I1)*DE
      enddo
C  PART III
C  1--> I2.  R2(1)=1
      DO I2=2,LUP
C  NO ELECTRON ENERGY LOSS IS ASSOCIATED WITH Q2
C  STRICTLY: TOGETHER WITH Q2 THERE SHOULD ALSO BE AN E_Q2
C  FOR EXTERNAL SOURCES.
        DE=0
        SUSCR=Q2(I2)*DE
        E_SCRRAD=E_SCRRAD+SUSCR
      ENDDO
C  PART IV
C  i1--> i2. i1<i2,  R2(i1)=...
      DO I1=2,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=r2(i1)*densel*C(I1,I2)*DE
          E_SCRRAD=E_SCRRAD+SUSCR
        ENDDO
      ENDDO
C  PART V
C  i2--> i1.  i2>i1, R2(i2)=...
C            (includes i1=1, i.e., invers to PART III)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))
          SUSCR  =r2(i2)*densel*F(I2,I1)*DE
          E_SCRRAD=E_SCRRAD+SUSCR
C  separte treatment of radiation losses alone, only for testing
C         SUSCR_T=r2(i2)*A(I2,I1)*(-1.)*DE
C         E_SCRRAD_T=E_SCRRAD_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only.  evaluate second formular for E_SCRRAD,
C                  using radiation loss E_SCRRAD_T and effective rate SCRRAD
C     UH=13.595
C     DE=(E_AT(1)-UH)
C     E_SCRRAD_T=E_SCRRAD_T+SCRRAD *DE

c  this test only works if in PART III one would set:
C       DE=(E_AT(1)-E_AT(I2))
C  But for external sources this is not necessarily the case
C  test done, 11.01.05,  o.k.,  then test switched off.

c  next: contribution for coupling to H+

C  TO BE WRITTEN  E_ALPRC=....

      RETURN
      END

C********************************************************************
C            C(I,J)  FROM  JOHNSON
C *******************************************************************
      SUBROUTINE JOHN(TEMP,CJ,OSC)
C
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 OSC(40,40),CJ(40,40)
      TE=TEMP*1.1605E4
C
      DO 1 I=1,40
      DO 2 J=1,40
      P=I
      Q=J
      IF (I.EQ.1) THEN
      BN=-0.603
      ELSE
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      END IF
      EP=13.6/P**2
      EQ=13.6/Q**2
      U=(EP-EQ)/TEMP
C
      IF (U.GT.0) THEN
      X=1-P**2/Q**2
      Y=U
      R=1.94/P**1.57*X
      Z=R+Y
      A=2*P**2*OSC(I,J)/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      C1=2*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)
      E1Y=-E1Y

      E1Z=-E1Z

      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y-(1/Z+0.5)*E1Z
      E2=E2Y/Y-E2Z/Z

      CJ(I,J)=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)
      ELSE
      CJ(I,J)=0.0

      END IF
C
    2 CONTINUE
    1 CONTINUE
C
C
      RETURN
      END
c
      subroutine LAX(A,N1,N,B,eps,ifl,is,vw,ip,icon)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      double precision a(n1,n1),B(*),vw(*),d
      dimension ip(*)
      dimension iw(100),indx(100)
      if (n1.gt.100) then
        write (iunout,*) 'error in lax'
      stop
      endif
      call galpd(a,n1,n,b,iw,ier)
      if (ier.eq.1) then
        write (iunout,*) 'error in lax, matrix ist singulaer'
      endif
      return
      end
c
      SUBROUTINE GALPD(A,NA,NG,B,IW,IER)
C
C***********************************************************************
C*   GAUSS-ALGORITHMUS ZUR LOESUNG LINEARER GLEICHUNGS-SYSTEME MIT     *
C*   PIVOTIERUNG.                                                      *
C*   GENAUIGKEIT:   DOUBLD-PRECISION                   (01.07.1991)    *
C***********************************************************************
C    A(NA,NA): KOEFFIZIENTEN-MATRIX DES GLEICHUNGS-SYSTEMS
C    NA      : DIMENSION VON A WIE IM AUFRUFENDEN PROGRAMM ANGEGEBEN
C    NG      : ANZAHL DER UNBEKANNTEN   (NG <= NA)
C    B(NG)   : ELEMENTE DER RECHTEN SEITE DES GLEICHUNGS-SYSTEMS
C    IW(NG)  : INTEGER-HILFS-ARRAY FUER EINE MOEGLICHE PROGRAMM-
C              INTERNE UMNUMERIERUNG DER GLEICHUNGEN
C    IER     : ERROR-INDEX (IER = 1: MATRIX SINGULAER)
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NA,NA),B(NG),IW(NG)
      DATA ZERO /1.D-71/
      IER=0
C
C     ******************************************************************
C     DER FALL:   NG = 2
C     ******************************************************************
C
      IF(NG.EQ.2) THEN
                  H=A(1,1)*A(2,2)-A(2,1)*A(1,2)
                  AI=B(1)*A(2,2)-B(2)*A(1,2)
                  AK=A(1,1)*B(2)-A(2,1)*B(1)
                  B(1)=AI/H
                  B(2)=AK/H
                  RETURN
                  ENDIF
C
C     ******************************************************************
C     NUMMERN DER UNBEKANNTEN AUF IW ABSPEICHERN.
C     ******************************************************************
C
      DO 1 K=1,NG
    1    IW(K)=K
C
C     ******************************************************************
C     DIE A-MATRIX AUF DREIECKS-FORM BRINGEN.
C     ******************************************************************
C
      DO 10 I=1,NG
C
C        ===============================================================
C        Pivot-Element suchen  (  Zeile IZ,  Spalte KS )
C        ===============================================================
C
         AP=0
         IZ=0
         KS=0
         DO 3 M=I,NG
            DO 2 N=I,NG
               AMN=DABS(A(M,N))
               IF(AMN.GT.AP) THEN
                             AP=AMN

                             IZ=M

                             KS=N

                             ENDIF
    2          CONTINUE
    3      CONTINUE
C
C        ============================
C        Zeilen umordnen, wenn IZ > I
C        ============================
C
         IF(IZ.GT.I) THEN
                     DO 4 N=I,NG
                        H=A(I,N)
                        A(I,N)=A(IZ,N)
    4                   A(IZ,N)=H
                     H=B(I)
                     B(I)=B(IZ)
                     B(IZ)=H
                     ENDIF
C
C        ===============================================
C        Spalten umordnen und Unbekannte neu numerieren
C        ===============================================
C
         IF(KS.GT.I) THEN
                     IH=IW(I)
                     IW(I)=IW(KS)
                     IW(KS)=IH
                     DO 5 M=1,NG
                        AK=A(M,KS)
                        A(M,KS)=A(M,I)
    5                   A(M,I)=AK
                     ENDIF
C
C        ===============================================================
C        Total-Pivotierung. Spalte i,  Zeile 1 ... i-1, i+1 ... NG
C        zu Null machen.
C        ===============================================================
C
         AP=A(I,I)
         IF(DABS(AP).LT.ZERO) GOTO 15
         AP=1/AP
         DO 8 M=1,NG
            IF(M.EQ.I) GOTO 8
            IF(DABS(A(M,I)).GT.ZERO) THEN
                                     Q=A(M,I)*AP
                                     DO 7 N=I,NG
    7                                   A(M,N)=A(M,N)-A(I,N)*Q
                                     B(M)=B(M)-B(I)*Q
                                     ENDIF
    8       CONTINUE
   10    CONTINUE
C
C     ******************************************************************
C     WENN A(N,N) ^= 0, KOENNEN DIE UNBEKANNTEN BERECHNET WERDEN.
C     -:  SIE WERDEN ZUNAECHST AUF A(N,I), I=1,...,N GESETZT.
C     -:  DANN DIE BERECHNETEN UNBEKANNTEN IN DER RICHTIGEN ANORDNUNG
C         AUF B(I) SCHREIBEN UND AN DAS AUFRUFENDE PROGRAMM ZURUECKGEBEN
C     ******************************************************************
C
      DO 12 M=1,NG
   12    A(M,M)=B(M)/A(M,M)
      DO 14 M=1,NG
         II=IW(M)
   14    B(II)=A(M,M)
C
       RETURN

C     ******************************************************************
C     MATRIX IST SINGULAER.
C       -: IER = 1 SETZTEN
C       -: RETURN
C     ******************************************************************
C
   15 IER=1
      RETURN
      END


      subroutine expi(x,ei,icon)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c  exponential integral
c  -int exp(-x)/x, von -x nach unendlich     x<0,  identisch mit
c  +int exp(x)/x,  von -unendl. bis x
c
c   PV +int exp(-x)/x von unendl bis -x     x>0 ,  identisch mit
c   PV +int exp(x)/x von -unendl bis x
      double precision mmdei,dei,xx
cdr   if (x.gt.0) then
        xx=x
        dei=-mmdei(xx,ier)
        icon=ier
        ei=dei
cdr   else
cdr     write (*,*) 'argument in expi lt.0, call exit'
cdr     stop
cdr   endif
      return
      end
c
      subroutine aqc8(a,b,f,epsa,epsr,nmin,nmax,S,err,n,icon)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c  S=integral von a bis b, der function f(x) (external).
c  epsa,epsr : absolute and relative errors, input
c   nmin,nmax  min u max anzahl der functionsaufrufe
c   (nmax<511, nmin>15)
c
c output
c S
c err: estim absolut error
c n  anzahl der functionsaufrufe
c icon: error code
      external f
      external midpnt
      call qromo(f,a,b,s,midpnt,epsa)
      return
      end



c
      DOUBLE PRECISION FUNCTION MMDEI (S,IER)
C  exponential integral function
c  imsl routine, dort: mmdei(ipot=2,s,ier), also:
c  s muss gt.0, mmdei ist dann: integral (s bis unendlich) von
c               exp(-t)/t dt
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      X=S
      Y=ABS(X)
      Z=0.25D+0*Y
      IF(Z-1.0D0)11,11,12
   11 VALUE=((((((((((((((((((((-.483702D-8*Z+.2685377D-7)*Z-.11703642D-
     1 6)*Z+.585911692D-6)*Z-.2843937873D-5)*Z+.1284394756D-4)*Z-.547380
     1 648948D-4)*Z+.21992775413732D-3)*Z-.8290046678016D-3)*Z+.00291878
     1 85699858)*Z-.0095523774824170)*Z+.028895941356815)*Z-.08026651003
     1 2735)*Z+.20317460364863)*Z-.46439909294939)*Z+.948148
     1 1480886)*Z-1.706666666669
     1          )*Z+2.6666666666702)*Z-3.5555555555546)*Z+3.999999999999
     1 4)*Z-3.9999999999996)*Z+.57721566490143+LOG(Y)
      VALUE=-VALUE
      GO TO 13
   12 Z=1.0D0/Z
      VALUE=(((((((((((((((((((-.77769007188383D-3*Z+.84295295893998D-2
     1 )*Z-.04272827083185)*Z+.13473041266261)*Z-.29671551148)*Z+.486188
     1 06480858)*Z-.61743468824936)*Z+.62656052544291)*Z-.52203502518727
     1 )*Z+.36771959879483)*Z-.22727998629908)*Z+.12965447884319)*Z
     1 -.72886262704894D-1)*Z+.043418637381012)*Z-.29244589614262D-1
     1 )*Z+.23433749581188D-1)*Z-.023437317333964)*Z+.03124999448124
     1 )*Z-.062499999910765)*Z+.24999999999935)*Z-.20933859981036D-14
c     if (y.gt.160.D0) then
c       value=0
c     else
        VALUE=VALUE*EXP(-Y)
c     endif
   13 VAL=VALUE
      MMDEI=VAL
      ier=0
      RETURN
      END

      SUBROUTINE QROMO(FUNC,A,B,SS,CHOOSE,eps)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (JMAX=14,JMAXP=JMAX+1,KM=4,K=KM+1)
      DIMENSION S(JMAXP),H(JMAXP)
      external choose,func
      H(1)=1.0D0
      DO 11 J=1,JMAX

        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0d0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) then
            RETURN
          endif
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      PAUSE 'Too many steps.'
      write(iunout,*) '(W) Too many steps.'
      END

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N

        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I

          DIF=DIFT

        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
ctk       IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)


        ELSE
          DY=D(NS)


          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
C
      SUBROUTINE MIDPNT(FUNC,A,B,S,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      external func
      save
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END
c
C ===== SOURCE: intp_adas.f
      function intp_adas(ad,p1,p2) result(res)

      use precision
      use comxs, only: adas_data

      implicit none

      type(adas_data), pointer :: ad
      real(dp), intent(in) :: p1, p2
      real(dp) :: res, rx, ry
      integer :: ide, ite
      
      if (p1 <= ad%temp(1)) then
        ite = 1
      else if (p1 >= ad%temp(ad%ntemp)) then
        ite = ad%ntemp-1
      else 
        call binsearch_2 (ad%temp, ad%ntemp, p1, ite)
      end if  

      if (p2 <= ad%dens(1)) then
        ide = 1
      else if (p2 >= ad%dens(ad%ndens)) then
        ide = ad%ndens-1
      else 
        call binsearch_2 (ad%dens, ad%ndens, p2, ide)
      end if  
      
      rx = (ad%temp(ite+1) - p1) * ad%dte(ite)
      ry = (ad%dens(ide+1) - p2) * ad%dde(ide)
      call bilinear_int (ad%fit(ite:ite+1,ide:ide+1), rx, ry, res)

      return
      end function intp_adas
C ===== SOURCE: mom_rate_coeff.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only 

      function mom_rate_coeff (ir, p1, p2, lexp, iprshft) result (rate)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp) :: rate, sngl_poly, dum(9), rcmin, rcmax, fp(6), q1, q2
      real(dp), save :: xlog10e, xln10, dsub
      integer :: jfexmn, jfexmx
      integer, save :: ifirst=0, ifsub=0  

      interface
        function intp_adas (ad,p1,p2) result(res)
          use precision
          use comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function intp_adas
      end interface


      if (.not.reacdat(ir)%lrtcmw) then
        write (iunout,*) ' no data for momentum weighted rate',
     .                   ' coefficient available for reaction ',ir
        call exit_own(1)
      end if

      rate = 0._dp

      if ((reacdat(ir)%rtcmw%ifit == 1) .or. 
     .    (reacdat(ir)%rtcmw%ifit == 2)) then
        rcmin  = reacdat(ir)%rtcmw%poly%rcmn
        rcmax  = reacdat(ir)%rtcmw%poly%rcmx
        fp     = reacdat(ir)%rtcmw%poly%fparm
        jfexmn = reacdat(ir)%rtcmw%poly%ifexmn
        jfexmx = reacdat(ir)%rtcmw%poly%ifexmx
      end if

      if (mod(iftflg(ir,2),100) == 10) then
       
        rate = reacdat(ir)%rtcmw%poly%dblpol(1,1)

      elseif (reacdat(ir)%rtcmw%ifit == 1) then

        rate = sngl_poly(reacdat(ir)%rtcmw%poly%dblpol(1:9,1),p1,
     .                   rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))

      else if (reacdat(ir)%rtcmw%ifit == 2) then
         
        if (ifsub == 0) then
          ifsub = 1
          dsub = log(1.e8_dp)
        end if

        q2 = p2
        if (iprshft > 0) q2 = q2 - dsub
         
        call dbl_poly(reacdat(ir)%rtcmw%poly%dblpol,p1,q2,rate,dum,1,9,
     .                rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
        
      else if (reacdat(ir)%rtcmw%ifit == 3) then
         
        if (ifirst == 0) then
          ifirst = 1
          xln10 = log(10._dp)
          xlog10e = 1._dp/xln10
        end if

        q1 = xlog10e*p1
        q2 = xlog10e*p2
        rate = intp_adas(reacdat(ir)%rtcmw%adas,q1,q2)

        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate
        end if
        
      end if
      

      return

      end function mom_rate_coeff
     
C ===== SOURCE: ngffmh.f
CX UNIX PORT - SCCS info: Module @(#)$Header: /home/boerner/Eirene-Repository/Eirene/volume-processes/ngffmh.f,v 1.1 2007/01/16 15:13:32 boerner Exp $ Date $Date: 2007/01/16 15:13:32 $
CX
!pb       REAL*8 FUNCTION NGFFMH(GAM2)
       FUNCTION NGFFMH(GAM2)
!pb       IMPLICIT REAL*8(A-H,O-Z)
       USE PRECISION
       IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ************** FORTRAN77 SUBROUTINE: NGFFMH  ************************
C
C  VERSION:  1.0
C
C  PURPOSE:
C
C  EVALUATES ELECTRON TEMPERATURE AND FREQUENCY AVERAGED HYDROGENIC
C  FREE FREE GAUNT FACTOR.
C  OBTAINED FROM INTERPOLATION OF KARZAS & LATTER (1959) FIG.6
C  FOR -3<LOG10(Z0*Z0*IH/KTE)<1. OUTSIDE THIS RANGE A VERY APPROXIMATE
C  EXTRAPOLATION IS PERFORMED WITH GFFMH=1 IN THE INFINITE LIMITS.
C
C  INPUT:
C       GAM2=Z0*Z0*IH/KTE
C  OUTPUT:
C       NGFFMH=MAXWELL AND FREQUENCY AVERAGED FREE-FREE GAUNT FACTOR.
C  ********* H.P.SUMMERS, JET         12 JAN 1987    ******************
C-----------------------------------------------------------------------
C
C UNIX-IDL CONVERSION:
C
C VERSION: 1.1                          DATE: 22-08-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST CONVERTED. NO CHANGES.
C
C VERSION: 1.2                          DATE: 13-09-99
C MODIFIED: Martin O'Mullane
C               - Define function name as real*8.
C
C-----------------------------------------------------------------------
       REAL(DP), INTENT(IN) :: GAM2
       REAL(DP) :: NGFFMH, GAM2L
       INTEGER :: K

!PB       DIMENSION GAM2LA(17),GA(17)
       REAL(DP) :: GAM2LA(17),GA(17)
       DATA GAM2LA/-3.0D0,-2.75D0,-2.50D0,-2.25D0,-2.00D0,-1.75D0,
     &-1.50D0,-1.25D0,-1.00D0,-0.75D0,-0.50D0,-0.25D0,0.00D0,0.25D0,
     &0.50D0,0.75D0,1.00D0/
       DATA GA/1.139D0,1.151D0,1.167D0,1.189D0,1.215D0,1.248D0,1.283D0,
     &1.326D0,1.370D0,1.411D0,1.431D0,1.436D0,1.433D0,1.415D0,1.379D0,
     &1.338D0,1.296D0/

       GAM2L=DLOG10(GAM2)
       IF(GAM2L.LE.-3.0D0)GO TO 30
       IF(GAM2L.GE.1.0D0)GO TO 40
C  INTERPOLATION REGION, LOCATE GAM2L IN GAM2LA ARRAY
       K=0
   20  K=K+1
       IF(GAM2L.GT.GAM2LA(K)) GO TO 20
       K=K-1
       IF(K.EQ.16)K=15
C  QUADRATIC INTERPOLATION
       NGFFMH=(GAM2L-GAM2LA(K+1))*(GAM2L-GAM2LA(K+2))/
     &((GAM2LA(K)-GAM2LA(K+1))*(GAM2LA(K)-GAM2LA(K+2)))*GA(K)+
     &(GAM2L-GAM2LA(K))*(GAM2L-GAM2LA(K+2))/
     &((GAM2LA(K+1)-GAM2LA(K))*(GAM2LA(K+1)-GAM2LA(K+2)))*GA(K+1)+
     &(GAM2L-GAM2LA(K))*(GAM2L-GAM2LA(K+1))/
     &((GAM2LA(K+2)-GAM2LA(K))*(GAM2LA(K+2)-GAM2LA(K+1)))*GA(K+2)
       RETURN
C  EXTRAPOLATION FOR LOW GAM2
   30  NGFFMH=1.0D0-0.417D0/GAM2L
       RETURN
   40  NGFFMH=1.0D0+0.296D0/GAM2L
       RETURN
       END

C ===== SOURCE: prep_rtcs.f
      subroutine prep_rtcs (ir, iflg, ji, je, p1, cf)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, iflg, ji, je
      real(dp), intent(in) :: p1
      real(dp), intent(out) :: cf(9)
      real(dp) :: dum
      type(poly_data), pointer :: rp

      select case (iflg)
      case (3)
        if (.not.reacdat(ir)%lrtc) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR RATE COEFFICIENT',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%rtc%poly

      case (4)
        if (.not.reacdat(ir)%lrtcmw) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' MOMENTUM-WEIGHTED RATE COEFFICIENT',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%rtcmw%poly

      case (5)
        if (.not.reacdat(ir)%lrtcew) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' ENERGY-WEIGHTED RATE COEFFICIENT',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%rtcew%poly

      case (6)
        if (.not.reacdat(ir)%loth) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' OTHER REACTION',ir
          CALL EXIT_OWN(1)
        END IF
        rp => reacdat(ir)%oth%poly

      case default
        write (iunout,*) ' call to prep_rtcs with wrong flag, iflg = ',
     .                     iflg
        write (iunout,*) ' 3 <= iflg <= 6 assumed '
        call exit_own(1)
      end select

      call dbl_poly (rp%dblpol,p1,0._dp,dum,cf,ji,je,
     .               rp%rcmn, rp%rcmx, rp%fparm, rp%ifexmn, rp%ifexmx)

      return
      end subroutine prep_rtcs
C ===== SOURCE: rate_coeff.f
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only 
!pb  07.12.06: double declaration of dsub removed

      function rate_coeff (ir, p1, p2, lexp, iprshft) result (rate)

!   lexp:     rate coefficient in cm**3/sec
!  .not.lexp: ln(rate coefficient in cm**3/sec)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp) :: rate, sngl_poly, dum(9), rcmin, rcmax, fp(6), q1, q2
      real(dp), save :: xlog10e, xln10, dsub
      integer :: jfexmn, jfexmx
      integer, save :: ifirst=0, ifsub=0 

      interface
        function intp_adas (ad,p1,p2) result(res)
          use precision
          use comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function intp_adas
      end interface


      if (.not.reacdat(ir)%lrtc) then
        write (iunout,*) ' no data for rate coefficient available',
     .                    ' for reaction ',ir
        call exit_own(1)
      end if

      rate = 0._dp

      if ((reacdat(ir)%rtc%ifit == 1) .or. 
     .    (reacdat(ir)%rtc%ifit == 2)) then
        rcmin  = reacdat(ir)%rtc%poly%rcmn
        rcmax  = reacdat(ir)%rtc%poly%rcmx
        fp     = reacdat(ir)%rtc%poly%fparm
        jfexmn = reacdat(ir)%rtc%poly%ifexmn
        jfexmx = reacdat(ir)%rtc%poly%ifexmx
      end if
        
      if (mod(iftflg(ir,2),100) == 10) then
       
        rate = reacdat(ir)%rtc%poly%dblpol(1,1)

      elseif (reacdat(ir)%rtc%ifit == 1) then

        rate = sngl_poly(reacdat(ir)%rtc%poly%dblpol(1:9,1),p1,
     .                   rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))

      else if (reacdat(ir)%rtc%ifit == 2) then
         
        if (ifsub == 0) then
          ifsub = 1
          dsub = log(1.e8_dp)
        end if

        q2 = p2
        if (iprshft > 0) q2 = q2 - dsub

        call dbl_poly(reacdat(ir)%rtc%poly%dblpol,p1,q2,rate,dum,1,9,
     .                rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
        
      else if (reacdat(ir)%rtc%ifit == 3) then
         
        if (ifirst == 0) then
          ifirst = 1
          xln10 = log(10._dp)
          xlog10e = 1._dp/xln10
        end if

        q1 = xlog10e*p1
        q2 = xlog10e*p2
        rate = intp_adas(reacdat(ir)%rtc%adas,q1,q2)

        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate
        end if
        
      end if
      
      return

      end function rate_coeff
     
C ===== SOURCE: rstern.f
C
C
      FUNCTION RSTERN (ER,B,IFLAG,P)
C
C     IFLAG=1:  H+  + H, PURELY REPULSIVE POTENTIAL
C     IFLAG=2:  MORSE POTENTIAL, He+ + He, H+ + Noble Gases, H+ + H2
C
C               PARAMETERS: P(1):
C                           P(2):
C                           P(3):
C
C  RSTERN(ER,B) IS THE LARGEST ROOT OF THE EQUATION:
C
C     FI(R):=1.-V(R)/ER-(B/R)**2=0.
C
C   HERE:  ER COLLISION ENERGY (EV)
C          B  IMPACT PARAMETER
C          V  INTERACTION POTENTIAL (EV)
C
C     DIMENSION RV0(3),RVM(3),RVW(3),VM(3),VW(3),VSW(3)

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG

      REAL(DP) :: FIW, RV, RTSAF, BQ, TOL, RSTERN, FI, DFI, FITEST,
     .          rup, rlw, vw, rw, r0
      SAVE
C     DATA IFIRST /0/
C     DATA RV0/0.,0.99699,2.18039/
C     DATA RVM/0.,1.4556 ,2.835539/
C     DATA RVW/0.,1.99515,3.490687/
C     DATA VM /0.,-2.     ,-2.7/
C     DATA VW /0.,-1.5    ,-2.025/
C     DATA VSW/0.,1.284688,1.4283/
C
C
      IF (IFLAG.EQ.1) THEN
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAF
C  FOR IFLAG=1, RLW=B IS A LOWER BOUND AND RUP=INF IS AN UPPER BOUND.
C  FIRST TRY TO FIND BETTER BOUNDS:
        RUP=B
C
10      RLW=RUP
        RUP=RLW*2.
        FITEST=FI(RUP,ER,B,IFLAG,P,DFI)
        IF (FITEST.LT.0.D0) GOTO 10
C
        TOL=1.D-7*(RUP-RLW)
        RSTERN=RTSAF(RLW,RUP,TOL,ER,B,IFLAG,P)
C
      ELSEIF (IFLAG.EQ.2) THEN
C  GENERALISED MORSE POTENTIAL
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAFE
C  FOR IFLAG=2 AND IFLAG=3, THE FOLLOWING INTERVAL CAN BE USED:
C  R0 IS THE ROOT OF THE INTERACTION POTENTIAL V(R): V(R0)=0.
C  RM IS THE RADIUS OF THE MINIMUM OF V(R)           V(RM)= VM = -EPS
C  RW IS THE RADIUS OF THE POINT OF INFLECTION       V(RM)= VW
C        RM=P(4)
         R0=P(5)
         RW=P(6)
         VW=P(8)
C
        BQ=B*B
        RV=R0
C
C  CASE 1: B < R0,
C
        IF (RV.GT.B) THEN
          RLW=B
          RUP=RV
        ENDIF
C
C  CASE 2: R0 < B,
C
        IF (RV.LT.B) THEN
          RLW=RV
          RUP=B
C  TRY TO IMPROVE LOWER BOUND RLW
          FIW=1.-VW/ER-BQ/(RW**2)
          IF (FIW.LT.0.D0) RLW=RW
        ENDIF
C
        TOL=1.D-7*(RUP-RLW)
        RSTERN=RTSAF(RLW,RUP,TOL,ER,B,IFLAG,P)
      ELSE
        GOTO 990
      ENDIF
C
      RETURN
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN FUNCTION RSTERN '
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: rtsaf.f
C
C
      FUNCTION RTSAF(X1,X2,XACC,ER,B,IFLAG,P)

C  TAKEN FROM: NUMERICAL RECIPES, W.H.PRESS ET AL.,
C              CAMBRIDGE UNIV. PRESS, 1989, P258
C  MODIFIED, TO FIND LARGEST ROOT, IN CASE MORE THAN ONE ROOTS
C  AND TO SPEED UP

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, X2, XACC, ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG
      INTEGER, PARAMETER :: MAXIT=100
      REAL(DP) :: DF, XDIST, F, FI, TEMP, DX, DXOLD, XH, RTSAF, XL
      INTEGER :: J

      XL=X1
      XH=X2
C
      F=FI(XH,ER,B,IFLAG,P,DF)
      XDIST=(XH-XL)
      DX=MIN(XDIST,F/DF)
      DXOLD=DX
      RTSAF=XH-DX
C
      F=FI(RTSAF,ER,B,IFLAG,P,DF)

      IF(F.LT.0.D0) THEN
        XL=RTSAF
      ELSE
        XH=RTSAF
      ENDIF
C
      DO 11 J=1,MAXIT
        IF(((RTSAF-XH)*DF-F)*((RTSAF-XL)*DF-F).GT.0.
     *      .OR. ABS(F+F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          RTSAF=XL+DX
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=RTSAF
          RTSAF=RTSAF-DX
        ENDIF
C
        IF(ABS(DX).LT.XACC) RETURN
C
        F=FI(RTSAF,ER,B,IFLAG,P,DF)
        IF(F.LT.0.D0) THEN
          XL=RTSAF
        ELSE
          XH=RTSAF
        ENDIF
11    CONTINUE
      WRITE (iunout,*) 'RTSAF EXCEEDING MAXIMUM ITERATIONS'
      RETURN
      END
C ===== SOURCE: scatang.f
C  6.12.05  COMMENTS, NAMES MODIFIED
C
      SUBROUTINE SCATANG (EREL,RAN,ELTHETA,CTTHETA,SIGM)
c     ***********************************************************
c     *                                                         *
c     *Programm zur Berechnung der Streuwinkel bei vorgegebener *
c     *Energie und Zufallszahl (Monte-Carlo-Simulation) für eine*
c     *atomare Reaktion.                                        *
c     *Die verschiedenen Reaktionen werden mittels der Daten von*
c     *P.S. Krstic und D.R. Schultz (Atomic and Plasma-Material *
c     *Interaction Data for Fusion), Differentielle Querschnitte*
c     *mit zugehörigen Winkeln, berechnet.                      *
c     *Das Programm berechnet Reaktionen zwischen Energien von  *
c     *0.1 eV und 100 eV (CM). Dabei sind 31 Energiewerte aus   *
c     *den o.g.Daten vorgegeben.                                *
c     *Input:    Dateien der jeweiligen Reaktion                *
c     *          REL. Reaktionsenergy: EREL                     *
c     *          Zufallszahl/Flag    : RAN                      *
c     *Output:   Verwendete Energie: energy                     *
c     *          Streuwinkel: eltheta (elastisch)               *
c     *                       cttheta (ladungsaustausch)        *
c     *Totale Streuquerschnitte + Momente für alle Energien:    *
c     *  sigma(i), gamma(i), moment(i), viscos(i)               *
c     *mit dem Feldindex i der Energie (zugeh. Vorgabe i=m).    *
c     *                                                         *
c     *Mögliche Energien: Feld energy(i)                        *
c     *Winkel der diff.Querschnitte abh.von Energie: theta(k,i) *
c     *mit k=1..768                                             *
c     *                                                         *
c     *(19.9.2000)                          Torsten Haberscheidt*
c     *                                                         *
c     *                                                         *
c     *modifications: 24.3.03:                                  *
c     *   remove "ct"-scattering angle evaluation               *
c     *   remove "el"-scattering angle evaluation for ran.lt.0  *
c     *           in this case: only SIGM =SIGM(EREL)           *
c     *           is returned                                   *
c     ***********************************************************

c     --Parameter--


      USE PRECISION
      USE CCONA

      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: EREL,RAN
      REAL(DP), INTENT(OUT) :: SIGM

      integer i,j,k,l,m,n,dat,IFIRST,IUN
      REAL(DP) :: z,w
      REAL(DP) :: energy(31)
c     REAL(DP) :: theta(768,31),el(770),ct(770),dtheta(770)
      REAL(DP) :: theta(768,31),el(770),dtheta(770)
      REAL(DP) :: sig,sigma(31)
c     REAL(DP) :: gam,gamma(31)
      REAL(DP) :: mom,moment(31),vis,viscos(31)
      REAL(DP) :: normel(770)
c     REAL(DP) :: normct(770)
      REAL(DP) :: elangle(768,31)
c     REAL(DP) :: ctangle(768,31)
      REAL(DP) :: eltheta, cttheta
      CHARACTER(10) :: FILNAM
      DATA IFIRST /0/
      SAVE

c     ***********************************************************
c     -- Einlesen der Daten --

      IF (IFIRST == 0) THEN
        IFIRST = 1

        iun=23

        do 10 i=1,31
          IF (I <= 10) THEN
            WRITE (FILNAM,'(a3,i1,a4)') 'el-',i-1,'.dat'
          ELSE
            WRITE (FILNAM,'(a3,i2,a4)') 'el-',i-1,'.dat'
          END IF

          open (iun,file=filnam)

          read(iUN,*)
          read(iUN,100) energy(i)
100       format(T19,E9.4)

          do 20 j=1,13
            read(iUN,*)
20        continue

          dat=0

          do 30 k=1,768
c           read(iun,200,end=5) theta(k,i), el(k), ct(k)
            read(iun,200,end=5) theta(k,i), el(k)
            dat= dat+1
30        continue

200       format(T5, E11.6, T20, E11.6, T35, E11.6)
5         continue
          close(iUN)

c       *******************************************************
c       --Verarbeitung der Daten--

c       Winkelintervalle

          dtheta(1) = 0.5*(theta(1,i) + theta(2,i))

          do 40 k=2,dat-1
            dtheta(k) = 0.5*(theta(k+1,i) - theta(k-1,i))
40        continue

          dtheta(dat) = pia - 0.5*(theta(dat-1,i) + theta(dat,i))

c       Streuquerschnitte

          sig=0
c         gam=0
          mom=0
          vis=0

          do 50 k=1,dat
            sig = sig + (dtheta(k) * el(k))
c           gam = gam + (dtheta(k) * ct(k))
            mom = mom + (dtheta(k) * el(k))*(1-cos(theta(k,i)))
            vis = vis + (dtheta(k) * el(k))*(sin(theta(k,i)))**2
50        continue

          sigma(i) = sig
c         gamma(i) = gam
          moment(i)= mom
          viscos(i)= vis

c       Normierung

          do 60 k=1,dat
            normel(k) = (dtheta(k) * el(k)) / sig
c           normct(k) = (dtheta(k) * ct(k)) / gam
60        continue

c       Generieren von Zahlen [0;1]

          elangle(1,i) = normel(1)
c         ctangle(1,i) = normct(1)

          do 70 k=2,dat
            elangle(k,i) = elangle(k-1,i) + normel(k)
c           ctangle(k,i) = ctangle(k-1,i) + normct(k)
70        continue



10      continue

      END IF   ! END OF IFIRST-BLOCK

c     ***********************************************************

c   Berechnung des Monte-Carlo-Winkels mit linearer Interpolation
c   still missing: interpolation in energy. currently: use fixed energy array

c  binary search in energy array. Find index m
c  to be done: use equidist. grid on log scale (see ORNL elastic web page)
      call searchbin(energy,31,erel,m)

c  elastic (both with and without charge transfer, IP-model)

      SIGM=SIGMA(M)

      if (ran.lt.0._DP) return

c  binary search in angle array, given the energy energy(m)
      if(ran.lt.elangle(1,m)) then
        eltheta = ran*theta(1,m)/elangle(1,m)
      else if (ran.eq.elangle(dat,m)) then
        eltheta = pia
      else
        call searchbin(elangle(1:dat,m),dat,ran,n)
        z=(ran - elangle(n,m))/(elangle(n+1,m) - elangle(n,m))+n

        if(z.gt.(dat-1))then
          eltheta = (z-n)*(pia - theta(n,m)) + theta(n,m)
        else
          eltheta = (z-n)*(theta(n+1,m) - theta(n,m))+theta(n,m)
        endif
      endif


C CHARGE TRANSFER COMPONENT ONLY
c currently not used

c     if(ran.lt.ctangle(1,m))then
c       cttheta = ran*theta(1,m)/ctangle(1,m)
c     else if(ran.eq.ctangle(dat,m)) then
c       cttheta = pia
c     else
c       call searchbin(ctangle(1:dat,m),dat,ran,l)
c       w=(ran - ctangle(l,m))/(ctangle(l+1,m) - ctangle(l,m))+l
c
c       if(w.gt.(dat-1)) then
c         cttheta = (w-l)*(pia - theta(l,m)) + theta(l,m)
c       else
c         cttheta = (w-l)*(theta(l+1,m) - theta(l,m))+theta(l,m)
c       endif
c     endif



      RETURN
      end
C ===== SOURCE: searchbin.f



      subroutine searchbin(xx,n,x,m)
c     ***********************************************************
c     * Ermittlung des Feldindex m mit vorgegebener Zahl x,     *
c     * so dass x zwischen xx(m) und xx(m+1)                    *
c     * Modifikation:   x <=xx(1) --> m=1  x >=xx(n) --> m=n    *
c     * *********************************************************
      USE PRECISION
      implicit none
      integer n,m
      REAL(DP) :: x,xx(n)
      integer bl,bm,bu

      if (x.le.xx(1)) then
        m=1
      else if (x.ge.xx(n)) then
        m=n
      else

c  binary search
        bl=0
        bu=n+1

80      if (bu-bl.gt.1) then
          bm=(bu+bl)*0.5
          if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(bm))) then
            bl=bm
          else
            bu=bm
          endif
          goto 80
        endif

        m=bl

      endif

      return
      end
C ===== SOURCE: setamd.f
C 27.6.05:  PHV_NROTA, PHV_NROTPH REMOVED
C
      SUBROUTINE SETAMD(ICAL)
C
C  SET ATOMIC AND MOLECULAR DATA: DRIVER
C
      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMSOU
      USE CZT1
      USE PHOTON

      IMPLICIT NONE

      real(dp) :: tpb1, tpb2, second_own
      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, IRPI, IRDS

!pb      tpb1 = second_own()

      IF (ICAL == 0) THEN
        NRCX=0
        NREL=0
        NRPI=0
        NRDS=0
        NREC=0
        NBGV=0
        NROT=0
        CALL XSECTA_PARAM
        CALL XSECTM_PARAM
        CALL XSECTI_PARAM
        CALL XSECTP_PARAM
        CALL XSECTPH_PARAM

        NRCX=MAX(1,NRCX)
        NREL=MAX(1,NREL)
        NRPI=MAX(1,NRPI)
        NRDS=MAX(1,NRDS)
        NREC=MAX(1,NREC)
        NBGV=MAX(1,NBGV)
        NROT=MAX(1,NROT)

        CALL SET_PARMMOD(2)
        CALL ALLOC_COMXS(2)
        CALL ALLOC_COMSOU(2)
        CALL ALLOC_CZT1(2)

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for setamd(0) ',tpb2-tpb1
!pb        tpb1 = tpb2

        RETURN
      ELSE
        CALL INIT_CMDTA(2)

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for init_cmdta ',tpb2-tpb1
!pb        tpb1 = tpb2

      END IF

      NRCXI=0
      NRELI=0
      NRPII=0
      NREII=0
      NRRCI=0
      NRBGI=0
      CALL XSECTA

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsecta ',tpb2-tpb1
!pb        tpb1 = tpb2

      CALL XSECTM

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectm ',tpb2-tpb1
!pb        tpb1 = tpb2

      CALL XSECTI

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsecti ',tpb2-tpb1
!pb        tpb1 = tpb2

      CALL XSECTP

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectp ',tpb2-tpb1
!pb        tpb1 = tpb2

      CALL XSECTPH

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectph ',tpb2-tpb1
!pb        tpb1 = tpb2

      CALL CONDENSE

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for condense ',tpb2-tpb1
!pb        tpb1 = tpb2


      IPATDS = 0
      IPMLDS = 0
      IPIODS = 0
      IPPLDS = 0
      DO IRDS=1,NRDS
        ipatds(IRDS,0)=COUNT(PATDS(IRDS,1:) > 0)
        IF (ipatds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03
             IPATDS(IRDS,1:ipatds(IRDS,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATDS(IRDS,1:) > 0)
        END IF
        ipmlds(IRDS,0)=COUNT(PMLDS(IRDS,1:) > 0)
        IF (ipmlds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03
             IPMLDS(IRDS,1:ipmlds(IRDS,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLDS(IRDS,1:) > 0)
        END IF
        ipiods(IRDS,0)=COUNT(PIODS(IRDS,1:) > 0)
        IF (ipiods(IRDS,0).GT.0) THEN         ! inserted by Derek Harting 26.03.
             IPIODS(IRDS,1:ipiods(IRDS,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIODS(IRDS,1:) > 0)
        END IF
        ipplds(IRDS,0)=COUNT(PPLDS(IRDS,1:) > 0)
        IF (ipplds(IRDS,0).GT.0) THEN         ! inserted by Derek Harting 26.03.
             IPPLDS(IRDS,1:ipplds(IRDS,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLDS(IRDS,1:) > 0)
        END IF
      END DO

      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
      IPPLPI = 0
      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IF (ipatpi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATPI(IRPI,1:) > 0)
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IF (ipmlpi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLPI(IRPI,1:) > 0)
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IF (ipiopi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIOPI(IRPI,1:) > 0)
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IF (ipplpi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLPI(IRPI,1:) > 0)
      END DO

!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for packs und counts ',tpb2-tpb1
!pb        tpb1 = tpb2


      RETURN
      END
C ===== SOURCE: setup_default_reactions.f
      subroutine setup_default_reactions

      use precision
      use parmmod
      use comxs

      implicit none

      integer :: ir


!  SPECIFY DEFAULT MODEL FOR RATE COEFFCIENTS

C
C -K=1:   E + HE --> 2E + HE+
C  RATE COEFFICIENT, JANEV, 2.3.9
      IR = -1
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
!pb      REACDAT(IR)%RTC%DBLPOL(1:9,1) =
!pb     . (/-4.409865E+01,2.391597E+01,-1.075323E+01,3.058039,
!pb     .   -5.685119E-01,6.795391E-02,-5.009056E-03,2.067236E-04,
!pb     .   -3.649161E-06/)
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     .  (/-4.409864886561d+01, 2.391596563469d+01,-1.075323019821d+01,
     .     3.058038757198d+00,-5.685118909884d-01, 6.795391233790d-02,
     .    -5.009056101857d-03, 2.067236157507d-04,-3.649161410833d-06/) 

C -K=2:   FREE
C -K=3:   FREE

C -K=4:   E + H --> H+ + 2E
C  RATE COEFFICIENT, JANEV, 2.1.5
      IR = -4
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-3.271397E+01,1.353656E+01,-5.739329E 00,1.563155E 00,
     .   -2.877056E-01,3.482560E-02,-2.631976E-03,1.119544E-04,
     .   -2.039150E-06/)

C -K=5:  E + H2 --> H + H + E
C  RATE COEFFICIENT, JANEV, 2.2.5, PREPRINT (CORRECT), NOT "BOOK"
      IR = -5
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-2.7872175E+01,1.0522527E+01,-4.9732123E+00,
     .    1.4511982E+00,-3.0627906E-01,4.4333795E-02,
     .   -4.0963442E-03, 2.1596703E-04,-4.9285453E-06/)

C -K=6:  E + H2 --> H+ + H + 2E
C  RATE COEFFICIENT, JANEV, 2.2.10
      IR = -6
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-3.834597E+01,1.426322E+01,-5.826467E+00,
     .    1.727941E+00,-3.598121E-01,4.822199E-02,
     .   -3.909403E-03,1.738777E-04,-3.252845E-06/)

C -K=7: E + H2 --> H2+(VIB) + 2E
C  RATE COEFFICIENT, JANEV, 2.2.9
      IR = -7
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-3.568640E+01,1.733469E+01,-7.767469E+00,
     .    2.211579E+00,-4.169840E-01,5.088290E-02,
     .   -3.832738E-03,1.612863E-04,-2.893392E-06/)

C -K=8: E + H2+(VIB) --> H + H+ + E
C  RATE COEFFICIENT, JANEV, 2.2.12
      IR = -8
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-1.781416E+01,2.277799E+00,-1.266868E+00,
     .    4.296170E-01,-9.609908E-02,1.387958E-02,
     .   -1.231349E-03,6.042383E-05,-1.247521E-06/)

C -K=9: E + H2+(VIB) --> H+ + H+ + 2E
C  RATE COEFFICIENT, JANEV, 2.2.11
      IR = -9
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-3.746192E+01,1.559355E+01,-6.693238E+00,
     .    1.981700E+00,-4.044820E-01,5.352392E-02,
     .   -4.317452E-03,1.918499E-04,-3.591779E-06/)

C -K=10: E + H2+(VIB) --> H + H(N)
C  RATE COEFFICIENT, JANEV, 2.2.14
      IR = -10
      ALLOCATE(REACDAT(IR)%RTC)
      ALLOCATE(REACDAT(IR)%RTC%POLY)
      ALLOCATE(REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LRTC = .TRUE.
      REACDAT(IR)%RTC%IFIT = 1
      REACDAT(IR)%RTC%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%RTC%POLY%FPARM = 0._DP
      REACDAT(IR)%RTC%POLY%IFEXMN = 0
      REACDAT(IR)%RTC%POLY%IFEXMX = 0
      REACDAT(IR)%RTC%POLY%DBLPOL(1:9,1) =
     . (/-1.670436E+01,-6.035645E-01,-1.942746E-08,
     .   -2.005952E-07,2.962996E-08,2.134293E-08,
     .   -6.353973E-09,6.152557E-10,-2.025362E-11/)

C -K=11:  FREE
C

      


!  SPECIFY DEFAULT MODEL FOR CROSS SECTION



C  K=-1:  H + H+ --> H+ + H   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 3.1.8
      IR = -1
      ALLOCATE(REACDAT(IR)%CRS)
      ALLOCATE(REACDAT(IR)%CRS%POLY)
      ALLOCATE(REACDAT(IR)%CRS%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LCRS = .TRUE.
      REACDAT(IR)%CRS%IFIT = 1
      REACDAT(IR)%CRS%POLY%DBLPOL(1:9,1) =
     . (/-3.274124E+01,-8.916457E-02,-3.016991E-02,
     .    9.205482E-03, 2.400267E-03,-1.927122E-03,
     .    3.654750E-04,-2.788866E-05, 7.422296E-07/)
      IFTFLG(IR,1) = 0
      REACDAT(IR)%CRS%POLY%RCMN = -2.3025851D+00
      REACDAT(IR)%CRS%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%CRS%POLY%FPARM(1) = -3.2945896D+01
      REACDAT(IR)%CRS%POLY%FPARM(2) = -1.713112D-01
      REACDAT(IR)%CRS%POLY%FPARM(3) = 0._DP
C  USE ASYMPTOTIC EXPRESSION NO. IFMN=5
      REACDAT(IR)%CRS%POLY%IFEXMN = 5
      REACDAT(IR)%CRS%POLY%IFEXMX = 5

C
C  K=-2:  He + He+ --> He+ + He   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 5.3.1
      IR = -2
      ALLOCATE(REACDAT(IR)%CRS)
      ALLOCATE(REACDAT(IR)%CRS%POLY)
      ALLOCATE(REACDAT(IR)%CRS%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LCRS = .TRUE.
      REACDAT(IR)%CRS%IFIT = 1
      REACDAT(IR)%CRS%POLY%DBLPOL(1:9,1) =
     . (/ -3.369296194290d+01,-8.324653178943d-02, 6.660151719388d-03,
     .    -3.592504363592d-03,-1.745382918016d-04, 1.497204460315d-04,
     .    -2.152122621503d-05, 1.473684503283d-06,-4.401831552698d-08/)
      IFTFLG(IR,1) = 0
      REACDAT(IR)%CRS%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%CRS%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%CRS%POLY%FPARM(1) = 0._DP
      REACDAT(IR)%CRS%POLY%FPARM(2) = 0._DP
      REACDAT(IR)%CRS%POLY%FPARM(3) = 0._DP
C  USE ASYMPTOTIC EXPRESSION NO. IFMN=5
      REACDAT(IR)%CRS%POLY%IFEXMN = 5
      REACDAT(IR)%CRS%POLY%IFEXMX = 5

C
C  K=-3:  He + He++ --> H++ + He   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 6.3.1
      IR = -3
      ALLOCATE(REACDAT(IR)%CRS)
      ALLOCATE(REACDAT(IR)%CRS%POLY)
      ALLOCATE(REACDAT(IR)%CRS%POLY%DBLPOL(1:9,1))
      REACDAT(IR)%LCRS = .TRUE.
      REACDAT(IR)%CRS%IFIT = 1
      REACDAT(IR)%CRS%POLY%DBLPOL(1:9,1) =
     . (/ -3.459818117569d+01,-8.748942423786d-02,-2.445604128495d-02,
     .     2.392295193337d-03, 9.876388162277d-04,-2.282012750308d-04,
     .     3.598361283629d-06, 1.940270105613d-06,-1.105794797036d-07/)
      IFTFLG(IR,1) = 0
      REACDAT(IR)%CRS%POLY%RCMN = -HUGE(1._DP)
      REACDAT(IR)%CRS%POLY%RCMX = HUGE(1._DP)
      REACDAT(IR)%CRS%POLY%FPARM(1) = 0._DP
      REACDAT(IR)%CRS%POLY%FPARM(2) = 0._DP
      REACDAT(IR)%CRS%POLY%FPARM(3) = 0._DP
C  USE ASYMPTOTIC EXPRESSION NO. IFMN=5
      REACDAT(IR)%CRS%POLY%IFEXMN = 5
      REACDAT(IR)%CRS%POLY%IFEXMX = 5



      end subroutine setup_default_reactions
C ===== SOURCE: xsecta.f
c 24.11.05: chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
C  6.12.05: comments changed: default cx only for H on p. No He default cx
C  2.5.06:  default resonant cx added for He on He+ and He on He++
C           also modified: cross.f, xsecta_param.f
! 30.08.06: array PLS and COUN changed to allocatable arrays
! 30.08.06: data structure for reaction data redefined          
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
C
      SUBROUTINE XSECTA
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR ATOMIC SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMSOU
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: CF(9,0:9)
      REAL(DP), ALLOCATABLE :: PLS(:)
      REAL(DP) :: FACTKK, CHRDF0, EELEC, RMASS, DSUB, DEIMIN, EHEAVY,
     .            EBULK, TMASS, PMASS, COU, RATE_COEFF
      INTEGER :: II, IML, IM, IIO, NTE, ISTORE, ISCND, ISCDE, IFRST,
     .           IAT, IREI, IATM, IDSC1, J, IPLS1, IPLS, IION1, NRC,
     .           KK, ISPZB, IAEL, ITYPB, IREL, IBGK, IA, ISP, IP,
     .           IAPI, IRPI, IACX, IDSC, IPL, IAEI, IESTM, IRCX, IPLSTI,
     .           ISTORE_MDCL
      INTEGER, EXTERNAL :: IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2
c slmod begin
      REAL(DP), ALLOCATABLE :: DHELP(:)
      LOGICAL :: OVERWRITE = .FALSE.
c slmod end
C
C
C   ELECTRON IMPACT COLLISIONS:
C
C  FIND SPECIES INDEX OF ION AFTER IONIZATION EVENT FOR THE DEFAULT
C  ELECTRON IMPACT IONIZATION MODELS FROM INPUT MASS AND
C  AND CHARGE NUMBER
C
C
      ALLOCATE (PLS(NSTORDR))

!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 70 J=1,NSBOX
!pb          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
          PLS(J)=MAX(DEIMIN,DEINL(J))
70      CONTINUE
      END IF
C
      DO 100 IATM=1,NATMI
        IDSC1=0
        LGAEI(IATM,0)=0
C
        DO NRC=1,NRCA(IATM)
          KK=IREACA(IATM,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
        IF (NRCA(IATM).EQ.0.AND.NCHARA(IATM).LE.2) THEN
C
C  DEFAULT H,D,T OR HE ELEC. IMP. IONIZATION MODEL
C
          IION1=0
          IPLS1=0
          DO 52 IPLS=1,NPLSI
            IF (NCHARP(IPLS).EQ.NCHARA(IATM).AND.
     .          NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .          NCHRGP(IPLS).EQ.1) THEN
              IPLS1=IPLS
C
              IDSC1=IDSC1+1
              NREII=NREII+1
              IREI=NREII
              LGAEI(IATM,IDSC1)=IREI
C
              PELDS(IREI)=1.
              PPLDS(IREI,IPLS1)=1.
              EPLDS(IREI,1)=1.D0
              EPLDS(IREI,2)=0.D0
              GOTO 50
            ENDIF
52        CONTINUE
          GOTO 100
C
50        CONTINUE
          NTE=NSBOX
          IF (NSTORDR < NRAD) NTE=1
          IF (NCHARA(IATM).EQ.1) THEN
c  hydrogenic atoms
            ISTORE=-4
            EELEC=-EIONH
          ELSEIF (NCHARA(IATM).EQ.2) THEN
c  helium atoms
            ISTORE=-1
            EELEC=-EIONHE
          ENDIF
C
          IF (NSTORDR >= NRAD) THEN
            DO 80 J=1,NSBOX
              COU = RATE_COEFF(ISTORE,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
80          CONTINUE
C  NO RADIATION LOSS INCLUDED
            EELDS1(IREI,1:NSBOX)=EELEC
            NREAEI(IREI) = ISTORE
            JEREAEI(IREI) = 1
            NELREI(IREI) = ISTORE
          ELSE
            FACREA(ISTORE) = 0._DP
            NREAEI(IREI) = ISTORE
            JEREAEI(IREI) = 1
            NELREI(IREI) = ISTORE
          ENDIF

          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
C
C  TRACKLENGTH ESTIMATOR FOR ALL COLLISION RATE CONTRIBUTIONS
C
          IESTEI(IREI,1)=0
          IESTEI(IREI,2)=0
          IESTEI(IREI,3)=0
C
          NAEII(IATM)=IDSC1
C
C  NON DEFAULT ELEC. IMP. COLLISION MODEL,
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 90 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.D0
            IAT=NSPH+IATM
            RMASS=RMASSA(IATM)
            IFRST=ISCD1A(IATM,NRC)
            ISCND=ISCD2A(IATM,NRC)
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EHEAVY=ESCD1A(IATM,NRC)+ESCD2A(IATM,NRC)
            EELEC=EELECA(IATM,NRC)
c slmod begin
c...        Just nasty, but using EELECA=999 to overwrite ionisation
c           data with data from previous run:
            IF (EELEC.EQ.999.0) THEN
              OVERWRITE = .TRUE.
              EELEC = 2.0         ! Correct I believe...
            ENDIF
c slmod end
            IDSC1=IDSC1+1
            NREII=NREII+1
            IREI=NREII
            LGAEI(IATM,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IAT,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
c slmod begin
            IF (OVERWRITE) THEN
              WRITE(0,*) 'HERE TO CALL PROUSR:',nrc,eeleca(iatm,nrc)
c...          Overwrite TABDS1(2,*) with data passed to the PROUSR routine:
              ALLOCATE(DHELP(NSURF))
              CALL PROUSR (DHELP,-100,0.0D0,0.0D0,0.0D0,0.0D0,
     .                                0.0D0,0.0D0,.FALSE.,NSURF)
              TABDS1(2,1:NSURF)=DHELP(1:NSURF)                  ! Index assumption (the '2')... 
              DEALLOCATE(DHELP)
            ENDIF
c slmod end
90        CONTINUE
          NAEII(IATM)=IDSC1
       ENDIF
C
        NAEIIM(IATM)=NAEII(IATM)-1
        LGAEI(IATM,0)=NAEII(IATM)
C
        DO 95 IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          CALL XSTEI_1(IREI)
95      CONTINUE
C
100   CONTINUE
C
C
C   CHARGE EXCHANGE:
C
C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IATM AND ANY IPLS
      DO 200 IATM=1,NATMI
        IDSC=0
        LGACX(IATM,0,0)=0
        LGACX(IATM,0,1)=0
C
C   DEFAULT MODEL 100 --- 129: RESONANT CX FOR H  + P,
C                 130 --- 139: RESONANT CX FOR HE + HE+,
C                 140 --- 149: RESONANT CX FOR HE + HE++,
C
        IF (NRCA(IATM).EQ.0) THEN
          DO 155 IPLS=1,NPLSI
C CHECK: "ATOMIC" COLLISION PARTNERS ONLY
            IF (NPRT(NSPAMI+IPLS).NE.1.OR.NPRT(NSPH+IATM).NE.1) GOTO 155
C
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1.AND.
     .          NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
C
C  FIND BULK SECONDARIES
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 155
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 125
124           CONTINUE
              GOTO 155
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
C  PROJECTILE MASS IS 1.
C  TARGET     MASS IS 1.
              PMASS=1.*PMASSA
              TMASS=1.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-1
              ISTORE_MDCL = -1
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE+ ION
C
C  FIND BULK SECONDARIES
              DO 131 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 133
131           CONTINUE
              GOTO 155
133           DO 134 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 135
134           CONTINUE
              GOTO 155

C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
135           CONTINUE
C  PROJECTILE MASS IS 4.
C  TARGET     MASS IS 4.
              PMASS=4.*PMASSA
              TMASS=4.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-2
              ISTORE_MDCL = -2
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.2) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE++ ION
C
C  FIND BULK SECONDARIES
              DO 141 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.2)
     .          GOTO 143
141           CONTINUE
              GOTO 155
143           DO 144 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 145
144           CONTINUE
              GOTO 155

C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
145           CONTINUE
C  PROJECTILE MASS IS 4.
C  TARGET     MASS IS 4.
              PMASS=4.*PMASSA
              TMASS=4.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-3
              ISTORE_MDCL = -3
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
            ELSE
              GOTO 155
            ENDIF


            IDSC=IDSC+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGACX(IATM,IDSC,0)=IRCX
            LGACX(IATM,IDSC,1)=IPLS
            N1STX(IRCX,1)=1
            N1STX(IRCX,2)=IAT
            N1STX(IRCX,3)=1
            N2NDX(IRCX,1)=4
            N2NDX(IRCX,2)=IPL
            N2NDX(IRCX,3)=1
            MODCOL(3,1,IRCX)=ISTORE_MDCL

            DEFCX(IRCX)=LOG(CVELI2*PMASS)
            EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  TRACKLENGTH ESTIMATOR FOR ALL COLLISION RATE CONTRIBUTIONS
C
            IESTCX(IRCX,1)=0
            IESTCX(IRCX,2)=0
            IESTCX(IRCX,3)=0
C
C  DEFAULT BULK ION ENERGY LOSS RATE = 1.5*TI+EDRIFT PER COLLISION
C
            IF (NSTORDR >= NRAD) THEN
              IPLSTI=MPLSTI(IPLS)
              DO 150 J=1,NSBOX
                EPLCX3(IRCX,J,1)=1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J)
150           CONTINUE
              NELRCX(IRCX) = -1
            ELSE
              NELRCX(IRCX) = -1
            END IF
C
            MODCOL(3,2,IRCX)=3
            MODCOL(3,4,IRCX)=3
C
155       CONTINUE
C
          NACXI(IATM)=IDSC
C
C  NON DEFAULT CX MODEL:
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 160 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.3) GOTO 160
C
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.D0
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IDSC=IDSC+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGACX(IATM,IDSC,0)=IRCX
            LGACX(IATM,IDSC,1)=IPLS
            FDLMCX(IRCX)=FLDLMA(IATM,NRC)
            IAT=NSPH+IATM
            IPL=IPLS
            RMASS=RMASSA(IATM)
            IFRST=ISCD1A(IATM,NRC)
            ISCND=ISCD2A(IATM,NRC)
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EBULK=EBULKA(IATM,NRC)
            CALL XSTCX(RMASS,IRCX,IAT,IPL,
     .                 IFRST,ISCND,EBULK,CHRDF0,
     .                 ISCDE,IESTM,KK,FACTKK)
C
160       CONTINUE
C
          NACXI(IATM)=IDSC
C  NO CX MODEL DEFINED
        ELSE
          NACXI(IATM)=0
        ENDIF
C
        NACXIM(IATM)=NACXI(IATM)-1
C
        LGACX(IATM,0,0)=0.
        DO 180 IACX=1,NACXI(IATM)
          LGACX(IATM,0,0)=LGACX(IATM,0,0)+LGACX(IATM,IACX,0)
180     CONTINUE
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IATM=1,NATMI
        IDSC=0
        LGAEL(IATM,0,0)=0
        LGAEL(IATM,0,1)=0
C
C   AT PRESENT NO DEFAULT MODEL
C
        IF (NRCA(IATM).EQ.0) THEN
          NAELI(IATM)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 230 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGAEL(IATM,IDSC,0)=IREL
            LGAEL(IATM,IDSC,1)=IPLS
C
C  SPECIAL TREATMENT: BGK COLLISIONS AMONGST TESTPARTICLES
            IF (IBGKA(IATM,NRC).NE.0) THEN
              IF (NPBGKA(IATM).EQ.0) THEN
                NRBGI=NRBGI+3
                IBGK=NRBGI/3
                NPBGKA(IATM)=IBGK
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKA(IATM)
              ELSE
                GOTO 999
              ENDIF
C  SELF OR CROSS COLLISION?
              ITYPB=IDEZ(IBGKA(IATM,NRC),1,3)
              ISPZB=IDEZ(IBGKA(IATM,NRC),3,3)
              IF (ITYPB.NE.1.OR.ISPZB.NE.IATM) THEN
C  CROSS COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKA(IATM,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C  BGK-COLLISION PARAMETERS DONE
C
            IAT=NSPH+IATM
            IPL=IPLS
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EBULK=EBULKA(IATM,NRC)
            CALL XSTEL(IREL,IAT,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE

          NAELI(IATM)=IDSC
        ENDIF
C
        NAELIM(IATM)=NAELI(IATM)-1
C
        LGAEL(IATM,0,0)=0.
        DO 280 IAEL=1,NAELI(IATM)
          LGAEL(IATM,0,0)=LGAEL(IATM,0,0)+LGAEL(IATM,IAEL,0)
280     CONTINUE
C
300   CONTINUE
C
C   GENERAL HEAVY PARTICLE IMPACT COLLISIONS
C
      CALL XSTAPI(PLS)
C
C
      DO 1000 IATM=1,NATMI
C
        IF (TRCAMD) THEN
          CALL MASBOX ('ATOMIC SPECIES IATM = '//TEXTS(NSPH+IATM))
          CALL LEER(1)
C
          IF (LGAEI(IATM,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 870 IAEI=1,NAEII(IATM)
              IREI=LGAEI(IATM,IAEI)
              CALL XSTEI_2(IREI)
870         CONTINUE
          ENDIF
C
C
          CALL LEER(2)
          IF (LGACX(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 890 IACX=1,NACXI(IATM)
              IRCX=LGACX(IATM,IACX,0)
              IPL =LGACX(IATM,IACX,1)
              CALL XSTCX_2(IRCX,IPL)
890         CONTINUE
          ENDIF
          IF (LGAEL(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 895 IAEL=1,NAELI(IATM)
              IREL=LGAEL(IATM,IAEL,0)
              IPL =LGAEL(IATM,IAEL,1)
              CALL XSTEL_2(IREL,IPL)
895         CONTINUE
          ENDIF
          CALL LEER(2)
          IF (LGAPI(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO GENERAL ION IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 885 IAPI=1,NAPII(IATM)
              IRPI=LGAPI(IATM,IAPI,0)
              IPLS=LGAPI(IATM,IAPI,1)
              CALL LEER(1)
              WRITE (iunout,*) 'GENERAL ION IMPACT REACTION NO. IRPI= ',
     .                          IRPI
              CALL LEER(1)
              WRITE (iunout,*) 'INCIDENT BULK ION: IPLS:'
              WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPLS)
              CALL LEER(1)
              WRITE (iunout,*) 'ELECTRONS: PELPI,EELPI'
              WRITE (iunout,*) 'EL      ', PELPI(IRPI),0.D0
              CALL LEER(1)
              WRITE (iunout,*) 'BULK ION SECONDARIES:'
              IF (IPPLPI(IRPI,0).GT.0) THEN
                WRITE (iunout,*) 'BULK IONS: PPLPI,EPLPI '
                DO IP=1,IPPLPI(IRPI,0)
                  IPL=IPPLPI(IRPI,IP)
                  WRITE (iunout,*) TEXTS(NSPAMI+IPL),PPLPI(IRPI,IPL),
     .                                          EPLPI(IRPI,IPL)
                ENDDO
              ELSE
                WRITE (iunout,*) 'NONE'
              ENDIF
              CALL LEER(1)
              WRITE (iunout,*) 'TEST PARTICLE SECONDARIES:'
              IF (P2NPI(IRPI).EQ.0.D0) THEN
                WRITE (iunout,*) 'NONE'
              ENDIF
              DO IA=1,IPATPI(IRPI,0)
                IAT=IPATPI(IRPI,IA)
                ISP=NSPH+IAT
                WRITE (iunout,*) 'ATOM     IATM= ',
     .                      TEXTS(ISP),PATPI(IRPI,IAT)
              ENDDO
              DO IM=1,IPMLPI(IRPI,0)
                IML=IPMLPI(IRPI,IM)
                ISP=NSPA+IML
                WRITE (iunout,*) 'MOLECULE IMOL= ',
     .                      TEXTS(ISP),PMLPI(IRPI,IML)
              ENDDO
              DO II=1,IPIOPI(IRPI,0)
                IIO=IPIOPI(IRPI,II)
                ISP=NSPAM+IIO
                WRITE (iunout,*) 'TEST ION IION= ',
     .                      TEXTS(ISP),PIOPI(IRPI,IIO)
              ENDDO
            CALL LEER(1)
885         CONTINUE
          ENDIF
        ENDIF
1000  CONTINUE

      DEALLOCATE (PLS)
C
      RETURN
C
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTA.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IATM = ',IATM
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IATM,IREL '
      WRITE (iunout,*) IATM,IREL,IPLS
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: xsecta_param.f
C  may 2006:  default helium resonant cx added
C
      SUBROUTINE XSECTA_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMSOU
      USE CTEXT
      USE COMXS
      USE CSPEI
      USE PHOTON

      IMPLICIT NONE

      INTEGER :: IPLS, NRC, IATM, IAT, KK, IPL, NNROT
C
      DO 100 IATM=1,NATMI
C
        IF (NRCA(IATM).EQ.0.AND.NCHARA(IATM).LE.2) THEN
C
C  DEFAULT H,D,T OR HE ELEC. IMP. IONIZATION MODEL
C
          DO 52 IPLS=1,NPLSI
            IF (NCHARP(IPLS).EQ.NCHARA(IATM).AND.
     .          NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .          NCHRGP(IPLS).EQ.1) THEN
              NRDS=NRDS+1
              GOTO 50
            ENDIF
52        CONTINUE
          GOTO 100
C
50        CONTINUE
C
C  NON DEFAULT ELEC. IMP. COLLISION MODEL,
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 90 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
            NRDS=NRDS+1
90        CONTINUE
        ENDIF
C
100   CONTINUE
C
C
C   CHARGE EXCHANGE:
C
      DO 200 IATM=1,NATMI
C
C   HYDROGENIC AND HELIUM DEFAULT MODEL 100 --- 140
C
        IF (NRCA(IATM).EQ.0) THEN
          DO 155 IPLS=1,NPLSI
C CHECK: "ATOMIC" COLLISION PARTNERS ONLY
            IF (NPRT(NSPAMI+IPLS).NE.1.OR.NPRT(NSPH+IATM).NE.1) GOTO 155
C
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1.AND.
     .          NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
C
C  FIND BULK SECONDARIES
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 155
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 125
124           CONTINUE
              GOTO 155
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
              NRCX=NRCX+1
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE+ ION
C
C  FIND BULK SECONDARIES
              DO 131 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 133
131           CONTINUE
              GOTO 155
133           DO 134 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 135
134           CONTINUE
              GOTO 155

C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
135           CONTINUE
              NRCX=NRCX+1
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.2) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE++ ION
C
C  FIND BULK SECONDARIES
              DO 141 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.2)
     .          GOTO 143
141           CONTINUE
              GOTO 155
143           DO 144 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 145
144           CONTINUE
              GOTO 155

C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
145           CONTINUE
              NRCX=NRCX+1

            ENDIF
155       CONTINUE
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            NRCX=NRCX+1
C
130       CONTINUE
C
C  NO CX MODEL DEFINED
        ELSE
        ENDIF
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IATM=1,NATMI
C
C   AT PRESENT NO DEFAULT MODEL
C
        IF (NRCA(IATM).EQ.0) THEN
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 230 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
            NREL=NREL+1
C
C  SPECIAL TREATMENT: BGK COLLISIONS AMONGST TESTPARTICLES
            IF (IBGKA(IATM,NRC).NE.0) THEN
              IF (NPBGKA(IATM).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE

        ENDIF
C
300   CONTINUE
csw
csw  COLLISIONS OF ATOMS WITH PHOTON BACKGROUND, OT - type
csw
cdr  not active
c      IF (.FALSE.) THEN
c      nnrot=0
c      do iatm=1,natmi
c         if(nrca(iatm) > 0) then
c            do nrc=1,nrca(iatm)
c               kk=ireaca(iatm,nrc)
c               if(iswr(kk) == 7) then
c                  nNROT=nNROT+1
c               endif
c            enddo
c         endif
c      enddo
c      call PH_ALLOC_XSECTA(nnrot)
c      END IF
csw
C
      CALL XSTAPI_PARAM
C
      RETURN
C
      END





C ===== SOURCE: xsecti.f
C  27.6.05  irds --> irei
c 24.11.05 use nprt(ispz) to check if iion is a molecular ion
c          otherwise he+ atomic ions could be confused with d2+ molecular
c          ions and then get assigned the wrong default collision model
c 24.11.05 chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
C
      SUBROUTINE XSECTI
C
C  TABLE FOR REACTION RATES FOR TEST IONS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: FACTKK, CHRDF0, RMASS, ACCMAS, ACCINV, DSUB, DEIMIN,
     .          EHEAVY, EELEC, EBULK, COU, RATE_COEFF
      INTEGER :: ICOUNT, IA1, IP2, IPLS, ITEST, IIO, IION, IDSC1,
     .           NRC, J, IPLS1, IPLS2, IATM, KK, IATM1, IATM2, ITYPB,
     .           ISPZB, III, IDSC, IREL, IBGK, IIDS, IERR, IMOL, IIEL,
     .           IIEI, IREI, IESTM, IFRST, ISCND, ISCDE, IPL, IICX,
     .           IDSC2, IRCX
      INTEGER, EXTERNAL :: IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2
C
!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
!pb          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
          PLS(J)=MAX(DEIMIN,DEINL(J))
10      CONTINUE
      END IF
C
C
C  SET TEST IONIC SPECIES ATOMIC AND MOLECULAR DATA;
C
C  STORE "DEFAULT DISSOCIATION MODEL" DATA
C  IN EACH CELL.
C  FOR HYDROGENIC MOLECULE IONS ONLY THIS MEANS: ZERO MFP,
C  INSTANTANOUS DECAY INTO ATOMS OR BULK IONS
C  FOR ALL OTHER SPECIES: INFINITE MFP, I.E. NO COLLISIONS
C
C
      DO 100 IION=1,NIONI
        IDSC1=0
        LGIEI(IION,0)=0
C
        DO NRC=1,NRCI(IION)
          KK=IREACI(IION,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
C   FIRST: DEAL WITH EI (ELECTRON IMPACT) COLLISIONS
C 
        IF (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2.AND.
     .      NPRT(NSPAM+IION).GT.1) THEN
C  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
C  USE NPRT.GT.1 TO DISTINGUISH ATOMIC FROM MOLECULAR IONS

C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
C  H2+:
          IF (NMASSI(IION).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
C  HD+:
          ELSEIF (NMASSI(IION).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
C  D2+:
          ELSEIF (NMASSI(IION).EQ.4) THEN
C  TEST: D2+ OR HT+, USE TEXTS(IION)
            IF (INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
C  D2+ TEST ION IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0) THEN
C  HT+ TEST ION IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (iunout,*) 'TEST ION NO ',IION,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT+:
          ELSEIF (NMASSI(IION).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
C  T2+:
          ELSEIF (NMASSI(IION).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS (MAY BE SPLITTED INTO 1A AND 1B)
          IF (IATM1.NE.IATM2) THEN
            FACTKK=0.5
            ICOUNT=1
          ELSE
            FACTKK=1.D0
            ICOUNT=2
          ENDIF
C
          IA1=IATM1
          IP2=IPLS2
7000      ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PATDS(IREI,IA1)=PATDS(IREI,IA1)+1.
          PPLDS(IREI,IP2)=PPLDS(IREI,IP2)+1.
          ACCMAS=ACCMAS+RMASSA(IA1)
          ACCMAS=ACCMAS+RMASSP(IP2)
          ACCINV=ACCINV+1./RMASSA(IA1)
          ACCINV=ACCINV+1./RMASSP(IP2)
          P2ND(IREI,NSPH+IA1)=P2ND(IREI,NSPH+IA1)+1.
          EATDS(IREI,IA1,1)=RMASSA(IA1)/ACCMAS
          EATDS(IREI,IA1,2)=1./RMASSA(IA1)/ACCINV
          EPLDS(IREI,      1)=RMASSP(IP2)/ACCMAS
          EPLDS(IREI,      2)=1./RMASSP(IP2)/ACCINV
          EATDS(IREI,0,    1)=EATDS(IREI,IA1,1)
          EATDS(IREI,0,    2)=EATDS(IREI,IA1,2)
          PELDS(IREI)=0.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 73 J=1,NSBOX
              COU = RATE_COEFF(-8,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)*FACTKK
73          CONTINUE
            EELDS1(IREI,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 8.6 EV
            EHVDS1(IREI,1:NSBOX)=8.6
            NREAEI(IREI) = -8
            JEREAEI(IREI) = 1
            NELREI(IREI) = -8
            NREAHV(IREI) = -4
          ELSE
            FACREA(-8) = FACTKK
            NREAEI(IREI) = -8
            JEREAEI(IREI) = 1
            NELREI(IREI) = -8
            NREAHV(IREI) = -4
          END IF
          IF (ICOUNT.EQ.1) THEN
            IA1=IATM2
            IP2=IPLS1
            ICOUNT=2
            GOTO 7000
          ENDIF
C  SECOND PROCESS
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PPLDS(IREI,IPLS1)=PPLDS(IREI,IPLS1)+1.
          PPLDS(IREI,IPLS2)=PPLDS(IREI,IPLS2)+1.
          EPLDS(IREI,1)=1.0
          EPLDS(IREI,2)=1.0
          PELDS(IREI)=1.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              COU = RATE_COEFF(-9,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
71          CONTINUE
            EELDS1(IREI,1:NSBOX)=-15.5
C  TRANSFERRED KINETIC ENERGY: 0.5 EV
            EHVDS1(IREI,1:NSBOX)=0.5
            NREAEI(IREI) = -9
            JEREAEI(IREI) = 1
            NELREI(IREI) = -9
            NREAHV(IREI) = -5
          ELSE
            FACREA(-9) = 1._DP
            NREAEI(IREI) = -9
            JEREAEI(IREI) = 1
            NELREI(IREI) = -9
            NREAHV(IREI) = -5
          END IF
C  THIRD PROCESS
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PATDS(IREI,IATM1)=PATDS(IREI,IATM1)+1.
          PATDS(IREI,IATM2)=PATDS(IREI,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IREI,NSPH+IATM1)=P2ND(IREI,NSPH+IATM1)+1.
          P2ND(IREI,NSPH+IATM2)=P2ND(IREI,NSPH+IATM2)+1.
          EATDS(IREI,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IREI,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IREI,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IREI,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IREI,0,    1)=EATDS(IREI,IATM1,1)+EATDS(IREI,IATM2,1)
          EATDS(IREI,0,    2)=EATDS(IREI,IATM1,2)+EATDS(IREI,IATM2,2)
          PELDS(IREI)=-1.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              COU = RATE_COEFF(-10,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
72          CONTINUE
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            EELDS1(IREI,1:NSBOX)=-0.88*TEIN(1:NSBOX)
C  TRANSFERRED KINETIC ENERGY: INGOING ELECTRON ENERGY
            EHVDS1(IREI,1:NSBOX)=0.88*TEIN(1:NSBOX)
            NREAEI(IREI) = -10
            JEREAEI(IREI) = 1
            NELREI(IREI) = -10
            NREAHV(IREI) = -6
          ELSE
            FACREA(-10) = 1._DP
            NREAEI(IREI) = -10
            JEREAEI(IREI) = 1
            NELREI(IREI) = -10
            NREAHV(IREI) = -6
          END IF
C
76        CONTINUE
C
          NIDSI(IION)=IDSC1
C
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=-NCHRGI(IION)
            IIO=NSPAM+IION
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EHEAVY=ESCD1I(IION,NRC)+ESCD2I(IION,NRC)
            EELEC=EELECI(IION,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IREI=NREII
            LGIEI(IION,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IIO,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
90        CONTINUE
C
          NIDSI(IION)=IDSC1
C
        ENDIF
C
        NIDSIM(IION)=NIDSI(IION)-1
        LGIEI(IION,0)=NIDSI(IION)
C
100   CONTINUE

C   SECONDLY: DEAL WITH CX (CHARGE EXCHANGE) COLLISIONS

      DO 200 IION=1,NIONI
        IDSC2=0
        LGICX(IION,0,0)=0
        LGICX(IION,0,1)=0
C
C  THERE ARE CURRENTLY NO DEFAULT CX RATES FOR TEST IONS
C
        IF (NRCI(IION).EQ.0) THEN
          NICXI(IION)=0
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 130 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=-NCHRGI(IION)
            IPLS=IDEZ(IBULKI(IION,NRC),3,3)
            IDSC2=IDSC2+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGICX(IION,IDSC2,0)=IRCX
            LGICX(IION,IDSC2,1)=IPLS
            IIO=NSPAM+IION
            IPL=IPLS
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            CALL XSTCX(RMASS,IRCX,IIO,IPL,IFRST,ISCND,EBULK,
     .                       CHRDF0,ISCDE,IESTM,KK,FACTKK)
C
130       CONTINUE
C
          NICXI(IION)=IDSC2
        ENDIF
C
        NICXIM(IION)=NICXI(IION)-1
        LGICX(IION,0,0)=0
        DO IICX=1,NICXI(IION)
          LGICX(IION,0,0)=LGICX(IION,0,0)+LGICX(IION,IICX,0)
        ENDDO
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IION=1,NIONI
        IDSC=0
        LGIEL(IION,0,0)=0
        LGIEL(IION,0,1)=0
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCI(IION).EQ.0) THEN
          NIELI(IION)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 230 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKI(IION,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGIEL(IION,IDSC,0)=IREL
            LGIEL(IION,IDSC,1)=IPLS

C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKI(IION,NRC).NE.0) THEN
              IF (NPBGKI(IION).EQ.0) THEN
                NRBGI=NRBGI+3
                IBGK=NRBGI/3
                NPBGKI(IION)=IBGK
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKI(IION)
              ELSE
                GOTO 999
              ENDIF
C  SELF OR CROSS COLLISION?
              ITYPB=IDEZ(IBGKI(IION,NRC),1,3)
              ISPZB=IDEZ(IBGKI(IION,NRC),3,3)
              IF (ITYPB.NE.3.OR.ISPZB.NE.IION) THEN
C  CROSS COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKI(IION,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C
C
            III=NSPAM+IION
            IPL=IPLS
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            CALL XSTEL(IREL,III,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE
C
          NIELI(IION)=IDSC
        ENDIF
C
        NIELIM(IION)=NIELI(IION)-1
C
        LGIEL(IION,0,0)=0.
        DO 280 IIEL=1,NIELI(IION)
          LGIEL(IION,0,0)=LGIEL(IION,0,0)+LGIEL(IION,IIEL,0)
280     CONTINUE
C
300   CONTINUE
C
      DO 1000 IION=1,NIONI
C
        DO 500 IIEI=1,NIDSI(IION)
          IREI=LGIEI(IION,IIEI)
          CALL XSTEI_1(IREI)
500     CONTINUE
C
        IF (TRCAMD) THEN
          CALL MASBOX ('TEST ION SPECIES IION = '//TEXTS(NSPAM+IION))
          CALL LEER(1)
          IF (LGICX(IION,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 215 IICX=1,NICXI(IION)
              IRCX=LGICX(IION,IICX,0)
              IPL =LGICX(IION,IICX,1)
              CALL XSTCX_2(IRCX,IPL)
215         CONTINUE
          ENDIF
C
          IF (LGIEI(IION,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISION'
            CALL LEER(1)
          ELSE
            DO 210 IIDS=1,NIDSI(IION)
              IREI=LGIEI(IION,IIDS)
              CALL XSTEI_2(IREI)
210         CONTINUE
          ENDIF
          CALL LEER(1)
          IF (LGIEL(IION,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 815 IIEL=1,NIELI(IION)
              IREL=LGIEL(IION,IIEL,0)
              IPL =LGIEL(IION,IIEL,1)
              CALL XSTEL_2(IREL,IPL)
815         CONTINUE
          ENDIF
        ENDIF
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IION,IPLS ',KK,IION,IPLS
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTI.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IION = ',IION
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (iunout,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (iunout,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IION,IREL '
      WRITE (iunout,*) IION,IREL,IPLS
      CALL EXIT_OWN(1)
      RETURN
      END
C ===== SOURCE: xsecti_param.f
C
C
      SUBROUTINE XSECTI_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: ICOUNT, NRC, KK, IION, IPLS, ITEST, IATM, IATM1,
     .           IATM2, IPLS1, IPLS2
C
      DO 100 IION=1,NIONI
C
        IF (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
C  H2+:
          IF (NMASSI(IION).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
C  HD+:
          ELSEIF (NMASSI(IION).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
C  D2+:
          ELSEIF (NMASSI(IION).EQ.4) THEN
C  TEST: D2+ OR HT+, USE TEXTS(IION)
            IF (INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
C  D2+ TEST ION IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0) THEN
C  HT+ TEST ION IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (iunout,*) 'TEST ION NO ',IION,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT+:
          ELSEIF (NMASSI(IION).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
C  T2+:
          ELSEIF (NMASSI(IION).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS (MAY BE SPLITTED INTO 1A AND 1B)
          IF (IATM1.NE.IATM2) THEN
            ICOUNT=1
          ELSE
            ICOUNT=2
          ENDIF
C
7000      CONTINUE
          NRDS=NRDS+1
          IF (ICOUNT.EQ.1) THEN
            ICOUNT=2
            GOTO 7000
          ENDIF
C  SECOND PROCESS
          NRDS=NRDS+1
C  THIRD PROCESS
          NRDS=NRDS+1
C
76        CONTINUE
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            NRDS=NRDS+1
90        CONTINUE
C
        ENDIF
C
100   CONTINUE

      DO 200 IION=1,NIONI
C
C  NO DEFAULT CX RATES
C
        IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 130 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            NRCX=NRCX+1
C
130       CONTINUE
C
        ENDIF
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IION=1,NIONI
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 230 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            NREL=NREL+1

C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKI(IION,NRC).NE.0) THEN
              IF (NPBGKI(IION).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE
C
        ENDIF
C
300   CONTINUE
C
      RETURN
C
      END
C ===== SOURCE: xsectm.f
C 27.6.05  irds --> irei
c 24.11.05 use nprt(ispz) to check if imol is really a molecule.
c          otherwise He atoms could be confused with d2 molecules,
c          if they accidentally are specified in the molecule block 4b
c 24.11.05 chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
C
      SUBROUTINE XSECTM
C
C  TABLE FOR REACTION RATES FOR MOLECULES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: ACCMAS, ACCINV, FACTKK, DEIMIN, EELEC, CHRDF0, DSUB,
     .          RMASS, EBULK, EHEAVY, COU, RATE_COEFF
      INTEGER :: ITEST, IATM, IPLS, IION, IA1, IP2, ION, ICOUNT,
     .           IION3, IDSC1, NRC, KK, J, IMOL, IPLS1, IPLS2, IPLS3,
     .           IATM1, IATM2, ITYPB, ISPZB, IMEL, IDSC, IREL, IBGK,
     .           IMDS, IERR, IMCX, ISCND, ISCDE, IESTM, IML, IFRST,
     .           IRCX, IREI, IPL, IDSC2
      INTEGER, EXTERNAL :: IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2
C
!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
!pb          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
          PLS(J)=MAX(DEIMIN,DEINL(J))
10      CONTINUE
      END IF
C
C
C  SET MOLECULAR DATA
C
C  STORE "DEFAULT DISSOCIATION MODEL" DATA
C  IN EACH CELL. THIS DEFAULT MODEL (JANEV/LANGER, PPPL) MAY BE USED
C  FOR HYDROGENIC MOLECULES ONLY.
C
C
      DO 100 IMOL=1,NMOLI
        IDSC1=0
        LGMEI(IMOL,0)=0
C
        DO NRC=1,NRCM(IMOL)
          KK=IREACM(IMOL,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
C  CHECK IF THIS REALLY IS A  MOLECULE: USE NPRT(ISPZ).GT.1?

        IF (NPRT(NSPA+IMOL).LE.1) THEN
          WRITE (IUNOUT,*) 'SEVERE INPUT ERROR DETECTED IN XSECTM: '
          WRITE (IUNOUT,*) 'IMOL= ',IMOL,' CARRIES ONLY ONE FLUX UNIT'
          WRITE (IUNOUT,*) 'EXIT CALLED FROM XSECTM '
        ENDIF

C  YES, "IMOL" IS A MOLECULE !

      
        IF (NRCM(IMOL).EQ.0.AND.NCHARM(IMOL).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2 DISSOCIATION AND IONIZATION
C  TO ALL HYDROGENIC MOLECULES
C
          EELEC=-EIONH2
C
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
          IPLS3=0
          IION3=0
C  H2:
          IF (NMASSM(IMOL).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
            DO 25 IION=1,NIONI
              IF (NMASSI(IION).EQ.2.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
25          CONTINUE
C  HD:
          ELSEIF (NMASSM(IMOL).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
            DO 35 IION=1,NIONI
              IF (NMASSI(IION).EQ.3.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
35          CONTINUE
C  D2:  (OR HT ?)
          ELSEIF (NMASSM(IMOL).EQ.4) THEN
C  TEST: D2 OR HT, USE TEXTS(IMOL)
            IF (INDEX(TEXTS(NSPA+IMOL),'D').NE.0) THEN
C  D2 MOLECULE IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
              DO 45 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .              INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
                  IION3=IION
                ENDIF
45            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPA+IMOL),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPA+IMOL),'T').NE.0) THEN
C  HT MOLECULE IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
              DO 48 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .             (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0)) THEN
                  IION3=IION
                ENDIF
48            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (iunout,*) 'MOLECULE NO ',IMOL,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT:
          ELSEIF (NMASSM(IMOL).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
            DO 55 IION=1,NIONI
              IF (NMASSI(IION).EQ.5.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
55          CONTINUE
C  T2:
          ELSEIF (NMASSM(IMOL).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
            DO 65 IION=1,NIONI
              IF (NMASSI(IION).EQ.6.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
65          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2*IION3
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
          PATDS(IREI,IATM1)=PATDS(IREI,IATM1)+1.
          PATDS(IREI,IATM2)=PATDS(IREI,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IREI,NSPH+IATM1)=P2ND(IREI,NSPH+IATM1)+1.
          P2ND(IREI,NSPH+IATM2)=P2ND(IREI,NSPH+IATM2)+1.
          EATDS(IREI,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IREI,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IREI,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IREI,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IREI,0,    1)=EATDS(IREI,IATM1,1)+EATDS(IREI,IATM2,1)
          EATDS(IREI,0,    2)=EATDS(IREI,IATM1,2)+EATDS(IREI,IATM2,2)
          PELDS(IREI)=0.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 70 J=1,NSBOX
              COU = RATE_COEFF(-5,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
70          CONTINUE
            EELDS1(IREI,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 6 EV
            EHVDS1(IREI,1:NSBOX)=6.
            NREAEI(IREI)=-5
            JEREAEI(IREI)=1
            NELREI(IREI)=-5
            NREAHV(IREI)=-2
          ELSE
            FACREA(-5) = 0._DP
            NREAEI(IREI)=-5
            JEREAEI(IREI)=1
            NELREI(IREI)=-5
            NREAHV(IREI)=-2
          END IF
C  SECOND PROCESS (MAY BE SPLITTED INTO 2A AND 2B)
          IF (IATM1.NE.IATM2) THEN
            FACTKK=0.5
            ICOUNT=1
          ELSE
            FACTKK=1.D0
            ICOUNT=2
          ENDIF
C
          IA1=IATM1
          IP2=IPLS2
73        ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
          PATDS(IREI,IA1)=PATDS(IREI,IA1)+1.
          PPLDS(IREI,IP2)=PPLDS(IREI,IP2)+1.
          ACCMAS=ACCMAS+RMASSA(IA1)
          ACCMAS=ACCMAS+RMASSP(IP2)
          ACCINV=ACCINV+1./RMASSA(IA1)
          ACCINV=ACCINV+1./RMASSP(IP2)
          P2ND(IREI,NSPH+IA1)=P2ND(IREI,NSPH+IA1)+1.
          EATDS(IREI,IA1,1)=RMASSA(IA1)/ACCMAS
          EATDS(IREI,IA1,2)=1./RMASSA(IA1)/ACCINV
          EPLDS(IREI,      1)=RMASSP(IP2)/ACCMAS
          EPLDS(IREI,      2)=1./RMASSP(IP2)/ACCINV
          EATDS(IREI,0,    1)=EATDS(IREI,IA1,1)
          EATDS(IREI,0,    2)=EATDS(IREI,IA1,2)
          PELDS(IREI)=1.0
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              COU = RATE_COEFF(-6,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)*FACTKK
71          CONTINUE
            EELDS1(IREI,1:NSBOX)=-25.0
C  TRANSFERRED KINETIC ENERGY: 10 EV
            EHVDS1(IREI,1:NSBOX)=10.0
            NREAEI(IREI) = -6
            JEREAEI(IREI) = 1
            NELREI(IREI) = -6
            NREAHV(IREI) = -3
          ELSE
            FACREA(-6) = LOG(FACTKK)
            NREAEI(IREI) = -6
            JEREAEI(IREI) = 1
            NELREI(IREI) = -6
            NREAHV(IREI) = -3
          END IF
          IF (ICOUNT.EQ.1) THEN
            IA1=IATM2
            IP2=IPLS1
            ICOUNT=2
            GOTO 73
          ENDIF
C
C  THIRD PROCESS
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
          ION=NSPAM+IION3
          PIODS(IREI,IION3)=PIODS(IREI,IION3)+1.
          P2ND(IREI,ION)=P2ND(IREI,ION)+1.
          EIODS(IREI,IION3,1)=1.
          EIODS(IREI,IION3,2)=0.
          EIODS(IREI,0,1)=1.
          EIODS(IREI,0,2)=0.
          PELDS(IREI)=1.0
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              COU = RATE_COEFF(-7,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
72          CONTINUE
            EELDS1(IREI,1:NSBOX)=EELEC
            NREAEI(IREI) = -7
            JEREAEI(IREI) = 1
            NELREI(IREI) = -7
          ELSE
            FACREA(-7) = 0._DP
            NREAEI(IREI) = -7
            JEREAEI(IREI) = 1
            NELREI(IREI) = -7
            EELDS1(IREI,1)=EELEC
          END IF
C
76        CONTINUE
          NMDSI(IMOL)=IDSC1
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 90 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.
            IML=NSPA+IMOL
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EHEAVY=ESCD1M(IMOL,NRC)+ESCD2M(IMOL,NRC)
            EELEC=EELECM(IMOL,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IREI=NREII
            LGMEI(IMOL,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IML,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
90        CONTINUE
          NMDSI(IMOL)=IDSC1
        ENDIF
C
        NMDSIM(IMOL)=NMDSI(IMOL)-1
        LGMEI(IMOL,0)=NMDSI(IMOL)
100   CONTINUE
C
C
C   CHARGE EXCHANGE:

      DO 200 IMOL=1,NMOLI
        IDSC2=0
        LGMCX(IMOL,0,0)=0
        LGMCX(IMOL,0,1)=0
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMCXI(IMOL)=0
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 130 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.
            IPLS=IDEZ(IBULKM(IMOL,NRC),3,3)
            IDSC2=IDSC2+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGMCX(IMOL,IDSC2,0)=IRCX
            LGMCX(IMOL,IDSC2,1)=IPLS
            IML=NSPA+IMOL
            IPL=IPLS
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL XSTCX(RMASS,IRCX,IML,IPL,
     .                 IFRST,ISCND,EBULK,CHRDF0,
     .                 ISCDE,IESTM,KK,FACTKK)
C
130       CONTINUE
C
          NMCXI(IMOL)=IDSC2
        ENDIF
C
        NMCXIM(IMOL)=NMCXI(IMOL)-1
C
        LGMCX(IMOL,0,0)=0
        DO 161 IMCX=1,NMCXI(IMOL)
          LGMCX(IMOL,0,0)=LGMCX(IMOL,0,0)+LGMCX(IMOL,IMCX,0)
161     CONTINUE
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IMOL=1,NMOLI
        IDSC=0
        LGMEL(IMOL,0,0)=0
        LGMEL(IMOL,0,1)=0
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMELI(IMOL)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 230 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKM(IMOL,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGMEL(IMOL,IDSC,0)=IREL
            LGMEL(IMOL,IDSC,1)=IPLS
C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKM(IMOL,NRC).NE.0) THEN
              IF (NPBGKM(IMOL).EQ.0) THEN
                NRBGI=NRBGI+3
                IBGK=NRBGI/3
                NPBGKM(IMOL)=IBGK
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKM(IMOL)
              ELSE
                GOTO 999
              ENDIF
C  SELF OR CROSS COLLISION?
              ITYPB=IDEZ(IBGKM(IMOL,NRC),1,3)
              ISPZB=IDEZ(IBGKM(IMOL,NRC),3,3)
              IF (ITYPB.NE.2.OR.ISPZB.NE.IMOL) THEN
C  CROSS COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKM(IMOL,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C
            IML=NSPA+IMOL
            IPL=IPLS
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL XSTEL(IREL,IML,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE
C
          NMELI(IMOL)=IDSC
        ENDIF
C
        NMELIM(IMOL)=NMELI(IMOL)-1
C
        LGMEL(IMOL,0,0)=0.
        DO 280 IMEL=1,NMELI(IMOL)
          LGMEL(IMOL,0,0)=LGMEL(IMOL,0,0)+LGMEL(IMOL,IMEL,0)
280     CONTINUE
C
C
300   CONTINUE
C
      DO 1000 IMOL=1,NMOLI
C
        DO 500 IMDS=1,NMDSI(IMOL)
          IREI=LGMEI(IMOL,IMDS)
          CALL XSTEI_1(IREI)
500     CONTINUE
C
        IF (TRCAMD) THEN
          CALL MASBOX ('MOLECULAR SPECIES IMOL = '//TEXTS(NSPA+IMOL))
          CALL LEER(1)
          IF (LGMCX(IMOL,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 880 IMCX=1,NMCXI(IMOL)
              IRCX=LGMCX(IMOL,IMCX,0)
              IPL =LGMCX(IMOL,IMCX,1)
              CALL XSTCX_2(IRCX,IPL)
880         CONTINUE
          ENDIF
C
          IF (LGMEI(IMOL,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 870 IMDS=1,NMDSI(IMOL)
              IREI=LGMEI(IMOL,IMDS)
              CALL XSTEI_2(IREI)
870         CONTINUE
          ENDIF
C
          IF (LGMEL(IMOL,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 895 IMEL=1,NMELI(IMOL)
              IREL=LGMEL(IMOL,IMEL,0)
              IPL =LGMEL(IMOL,IMEL,1)
              CALL XSTEL_2(IREL,IPL)
895         CONTINUE
          ENDIF
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE '
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IMOL,IPLS ',KK,IMOL,IPLS
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTM.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IMOL = ',IMOL
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (iunout,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (iunout,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IMOL,IREL '
      WRITE (iunout,*) IMOL,IREL,IPLS
      CALL EXIT_OWN(1)
      RETURN
C
      END
C ===== SOURCE: xsectm_param.f
C
C
      SUBROUTINE XSECTM_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: ITEST, ICOUNT, NRC, IION, IMOL, IATM, IPLS, IATM1,
     .           IATM2, IPLS1, IPLS2, IPLS3, IION3, KK
C
      DO 100 IMOL=1,NMOLI
C
        IF (NRCM(IMOL).EQ.0.AND.NCHARM(IMOL).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2 DISSOCIATION AND IONIZATION
C  TO ALL HYDROGENIC MOLECULES
C
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
          IPLS3=0
          IION3=0
C  H2:
          IF (NMASSM(IMOL).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
            DO 25 IION=1,NIONI
              IF (NMASSI(IION).EQ.2.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
25          CONTINUE
C  HD:
          ELSEIF (NMASSM(IMOL).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
            DO 35 IION=1,NIONI
              IF (NMASSI(IION).EQ.3.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
35          CONTINUE
C  D2:  (OR HT ?)
          ELSEIF (NMASSM(IMOL).EQ.4) THEN
C  TEST: D2 OR HT, USE TEXTS(IMOL)
            IF (INDEX(TEXTS(NSPA+IMOL),'D').NE.0) THEN
C  D2 MOLECULE IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
              DO 45 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .              INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
                  IION3=IION
                ENDIF
45            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPA+IMOL),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPA+IMOL),'T').NE.0) THEN
C  HT MOLECULE IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
              DO 48 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .             (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0)) THEN
                  IION3=IION
                ENDIF
48            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (iunout,*) 'MOLECULE NO ',IMOL,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT:
          ELSEIF (NMASSM(IMOL).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
            DO 55 IION=1,NIONI
              IF (NMASSI(IION).EQ.5.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
55          CONTINUE
C  T2:
          ELSEIF (NMASSM(IMOL).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
            DO 65 IION=1,NIONI
              IF (NMASSI(IION).EQ.6.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
65          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2*IION3
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS
          NRDS=NRDS+1
C  SECOND PROCESS (MAY BE SPLITTED INTO 2A AND 2B)
          IF (IATM1.NE.IATM2) THEN
            ICOUNT=1
          ELSE
            ICOUNT=2
          ENDIF
C
73        CONTINUE
          NRDS=NRDS+1
          IF (ICOUNT.EQ.1) THEN
            ICOUNT=2
            GOTO 73
          ENDIF
C
C  THIRD PROCESS
          NRDS=NRDS+1
C
76        CONTINUE
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 90 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            NRDS=NRDS+1
90        CONTINUE
        ENDIF
C
100   CONTINUE
C
C
C   CHARGE EXCHANGE:

      DO 200 IMOL=1,NMOLI
C
        IF (NRCM(IMOL).EQ.0) THEN
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 130 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            NRCX=NRCX+1
C
130       CONTINUE
C
        ENDIF
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IMOL=1,NMOLI
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMELI(IMOL)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 230 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            NREL=NREL+1
C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKM(IMOL,NRC).NE.0) THEN
              IF (NPBGKM(IMOL).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE
C
        ENDIF
C
C
C
300   CONTINUE
C
      RETURN
C
      END
C ===== SOURCE: xsectp.f
C  aug. 05:  corrected electron energy loss rate for default rec. rate
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
C
      SUBROUTINE XSECTP
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR BULK ION SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI
      USE PHOTON

      IMPLICIT NONE
C
      REAL(DP) :: PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: DELE, FCTKKL, EEMX, ZX, DSUB, DEIMIN, RMASS2, FACTKK,
     .            RMASS2_2, CORSUM, COU, RATE_COEFF, ENERGY_RATE_COEFF,
     .            BREMS, TOT_BREMS, Z, ngffmh
      INTEGER :: IIRC, IION3, IPLS3, IATM3, IMOL3, KK, NRC, IATM,
     .           IRRC, J, IDSC, IPLS, NSERC5, KREAD, I, MODC, IATM1,
     .           ITYP, ISPZ, ITYP2, ISPZ2, IPHOT3
      INTEGER, EXTERNAL :: IDEZ
      LOGICAL :: LEXP, LADAS
      SAVE
C
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 70 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))
70      CONTINUE
      END IF
C
C   RECOMBINATION
C
      DO 1000 IPLS=1,NPLSI
C
        IDSC=0
        LGPRC(IPLS,0)=0
C
        DO NRC=1,NRCP(IPLS)
          KK=IREACP(IPLS,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.7) GOTO 994
        ENDDO
C
        IF (NRCP(IPLS).EQ.0) THEN
C
          IF (NCHARP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
C
C  DEFAULT HYDROGENIC RECOMBINATION MODEL
C  HYDR. RECOMBINATION RATE-COEFFICIENT (1/S/CCM) E + H+ --> H + RAD.
C  GORDEEV ET. AL., PIS'MA ZH. EHKSP. TEOR. FIZ. 25 (1977) 223.
C
            DO 52 IATM=1,NATMI
              IF (NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .                            NCHRGP(IPLS).EQ.1) THEN
C
                IDSC=IDSC+1
                NRRCI=NRRCI+1
                IF (NRRCI.GT.NREC) GOTO 992
                IRRC=NRRCI
                LGPRC(IPLS,IDSC)=IRRC
                IF (NSTORDR >= NRAD) THEN
                  DO 51 J=1,NSBOX
                    ZX=EIONH/MAX(1.E-5_DP,TEIN(J))
C  rate = rate coeff: <sig v> times electr. density
                    TABRC1(IRRC,J)=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(J)
C  maxw. electron energy loss rate due to recombination
c                   corsum=0._dp  !  old default: 1.5*Te
C  correction due to energy dependence in rec. cross section
C  corsum=d(ln<sig v>)/d(ln Te)
c  corsum approx -0.5 for Te --> 0
c  corsum approx  0.0 for Te approx 11.5
c  corsum approx +1.0 for Te --> infty
                    corsum=(-0.5_dp*zx+0.59)/(zx+0.59)
                    EELRC1(IRRC,J)=-(1.5+CORSUM)*TEIN(J)*TABRC1(IRRC,J)
51                CONTINUE
                  NREARC(IRRC) = 0
                  JEREARC(IRRC) = 0
                  NELRRC(IRRC) = -1
                ELSE
                  NREARC(IRRC) = 0
                  JEREARC(IRRC) = 0
                  NELRRC(IRRC) = -1
                END IF
                IATM1=IATM
                NATPRC(IRRC)=IATM1
                NIOPRC(IRRC)=0
                NPLPRC(IRRC)=0
                NMLPRC(IRRC)=0
C
                MODCOL(6,2,IRRC)=1
                MODCOL(6,4,IRRC)=1
              ENDIF
52          CONTINUE
C
            NPRCI(IPLS)=IDSC
          ENDIF
C
C  NON DEFAULT MODEL:  240--
C
        ELSEIF (NRCP(IPLS).GT.0) THEN
          DO 82 NRC=1,NRCP(IPLS)
            KK=IREACP(IPLS,NRC)
csw check photonic process
            if(iswr(kk)==7) then
               idsc=idsc+1
               nrrci=nrrci+1
               IF (NRRCI.GT.NREC) GOTO 992
               call XSTRC(ipls,nrc,idsc,nrrci)
               cycle
csw end branch
            ELSEIF (ISWR(KK).EQ.6) THEN
C
              FACTKK=FREACP(IPLS,NRC)
              IF (FACTKK.EQ.0.D0) FACTKK=1.
C  RECOMBINATION MODEL FOR BULK IONS
              IDSC=IDSC+1
              NRRCI=NRRCI+1
              IF (NRRCI.GT.NREC) GOTO 992
              IRRC=NRRCI
              LGPRC(IPLS,IDSC)=IRRC
C
              ITYP=IDEZ(ISCD1P(IPLS,NRC),1,3)
              ISPZ=IDEZ(ISCD1P(IPLS,NRC),3,3)
              IF (ITYP.EQ.3) THEN
                NIOPRC(IRRC)=ISPZ
                RMASS2=RMASSI(ISPZ)
              ELSEIF (ITYP.EQ.4) THEN
                NPLPRC(IRRC)=ISPZ
                RMASS2=RMASSP(ISPZ)
              ELSEIF (ITYP.EQ.1) THEN
                NATPRC(IRRC)=ISPZ
                RMASS2=RMASSA(ISPZ)
              ELSEIF (ITYP.EQ.2) THEN
                NMLPRC(IRRC)=ISPZ
                RMASS2=RMASSM(ISPZ)
              ELSEIF (ITYP.EQ.0) THEN
                NPHPRC(IRRC)=ISPZ
                RMASS2=0.
              ENDIF

              ITYP2=IDEZ(ISCD2P(IPLS,NRC),1,3)
              ISPZ2=IDEZ(ISCD2P(IPLS,NRC),3,3)
              IF (ITYP2.EQ.4) THEN
                NPLPRC_2(IRRC)=ISPZ2
                RMASS2_2      =RMASSP(ISPZ2)
              ELSE
                RMASS2_2=0._DP
              ENDIF
C  CHECK MASS CONSERVATION
              IF (RMASSP(IPLS).NE.(RMASS2+RMASS2_2)) GOTO 993
C
C  1.) CROSS SECTION(TE)
C           NOT NEEDED
C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3) --> RATE (1/S)
C
C  2.A) RATE COEFFICIENT = CONST.
C           TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
              IF (IDEZ(MODCLF(KK),3,5).EQ.1) THEN
                IF (NSTORDR >= NRAD) THEN
                  LEXP = .NOT. (MOD(IFTFLG(KK,2),100) == 10)
                  DO J=1,NSBOX
                    COU = RATE_COEFF(KK,TEINL(J),0._DP,LEXP,0)
                    TABRC1(IRRC,J)=COU*FACTKK
                    IF (IFTFLG(KK,2) < 100)
     .                TABRC1(IRRC,J)=TABRC1(IRRC,J)*DEIN(J)
                  END DO
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 1
                ELSE
C  DON'T STORE DATA, BUT COMPUTE THEM THEN NEEDED
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 1
                  FACREA(KK) = FACTKK
                END IF
                MODCOL(6,2,IRRC)=1
C             ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT
              ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
                IF (NSTORDR >= NRAD) THEN
                  DO J=1,NSBOX
                    COU = RATE_COEFF(KK,TEINL(J),PLS(J),.TRUE.,1)
                    TABRC1(IRRC,J)=COU*FACTKK
                    IF (IFTFLG(KK,2) < 100)
     .                TABRC1(IRRC,J)=TABRC1(IRRC,J)*DEIN(J)
                  END DO

                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                ELSE
C  DON'T STORE DATA, BUT COMPUTE THEM WHEN NEEDED
                  FACREA(KK) = FACTKK
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                END IF
                MODCOL(6,2,IRRC)=1
              ENDIF
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ELECTRON ENERGY LOSS RATE
C
              NSERC5=IDEZ(ISCDEP(IPLS,NRC),5,5)
              IF (NSERC5.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATECOEFF.
                IF (NSTORDR >= NRAD) THEN
                  DO 101 J=1,NSBOX
                    EELRC1(IRRC,J)=EELECP(IPLS,NRC)*TABRC1(IRRC,J)
101               CONTINUE
                  NELRRC(IRRC) = -2
                ELSE
                  NELRRC(IRRC) = -2
                  EELRC1(IRRC,1)=EELECP(IPLS,NRC)
                END IF
                MODCOL(6,4,IRRC)=1
              ELSEIF (NSERC5.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ELECTRON = -1.5*TE*RATECOEFF.
                IF (NSTORDR >= NRAD) THEN
                  DO 102 J=1,NSBOX
                    EELRC1(IRRC,J)=-1.5*TEIN(J)*TABRC1(IRRC,J)
102               CONTINUE
                  NELRRC(IRRC) = -3
                ELSE
                  NELRRC(IRRC) = -3
                END IF
                MODCOL(6,4,IRRC)=1
              ELSEIF (NSERC5.EQ.3) THEN
C  4.C)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELECP(IPLS,NRC)
                MODC=IDEZ(MODCLF(KREAD),5,5)
                LADAS = IS_RTCEW_ADAS(KREAD)
                Z = NCHRGP(IPLS)
                TOT_BREMS = 0._DP
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >= NRAD) THEN

                    DO J = 1, NSBOX
                      IF (LGVAC(J,NPLS+1).OR.(NCHRGP(IPLS)==0)) THEN
                        BREMS = 0._DP
                      ELSE
                        BREMS = 1.54E-32_DP * TEIN(J)**0.5 * Z**2 *
     .                          ngffmh(Z**2 * 13.6_DP/TEIN(J)) *
     .                          DEIN(J)*FACTKK/ELCHA
                      END IF
                      TOT_BREMS = TOT_BREMS + BREMS*DIIN(IPLS,J)
                      EELRC1(IRRC,J)=ENERGY_RATE_COEFF(KREAD,TEINL(J),
     .                               0._DP,.TRUE.,0)*DEIN(J)*FACTKK
                      IF (LADAS) EELRC1(IRRC,J) = EELRC1(IRRC,J) + BREMS
                    END DO

                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                    FACREA(KREAD) = LOG(FACTKK)
                  END IF
                  MODCOL(6,4,IRRC)=1
C  4.D)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C               ELSEIF (MODC.EQ.2) THEN
C        IRRELEVANT
C                 MODCOL(6,4,IRRC)=2
C  4.E)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    FCTKKL=LOG(FACTKK)
                    DO J = 1, NSBOX
                      IF (LGVAC(J,NPLS+1).OR.(NCHRGP(IPLS)==0)) THEN
                        BREMS = 0._DP
                      ELSE
                        BREMS = 1.54E-32_DP * TEIN(J)**0.5 * Z**2 *
     .                          ngffmh(Z**2 * 13.6_DP/TEIN(J)) *
     .                          DEIN(J)*FACTKK/ELCHA
                      END IF
                      TOT_BREMS = TOT_BREMS + BREMS*DIIN(IPLS,J)
                      EELRC1(IRRC,J)=ENERGY_RATE_COEFF(KREAD,TEINL(J),
     .                               PLS(J),.FALSE.,1)
                      EEMX=MAX(-100._DP,EELRC1(IRRC,J)+FCTKKL+DEINL(J))
                      EELRC1(IRRC,J)=-EXP(EEMX)
                      IF (LADAS) EELRC1(IRRC,J) = EELRC1(IRRC,J) + BREMS
                    END DO

                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                    FACREA(KREAD) = LOG(FACTKK)
                  END IF
                  MODCOL(6,4,IRRC)=1
                ENDIF
                WRITE (IUNOUT,*) ' BREMSSTRAHLUNG = ', TOT_BREMS
                IF (DELPOT(KREAD).NE.0.D0) THEN
                  DELE=DELPOT(KREAD)
                  IF (NSTORDR >= NRAD) THEN
                    DO 110 J=1,NSBOX
                      EELRC1(IRRC,J)=EELRC1(IRRC,J)+
     .                               DELE*TABRC1(IRRC,J)
 110                CONTINUE
                  END IF
                ENDIF
              ENDIF
            ENDIF
C
82        CONTINUE
          NPRCI(IPLS)=IDSC
C
C  NO MODEL DEFINED
        ELSE
          NPRCI(IPLS)=0
        ENDIF
C
        NPRCIM(IPLS)=NPRCI(IPLS)-1
        LGPRC(IPLS,0)=NPRCI(IPLS)
C
        IF (TRCAMD) THEN
          CALL MASBOX ('BULK ION SPECIES IPLS = '//TEXTS(NSPAMI+IPLS))
          CALL LEER(1)
          IF (LGPRC(IPLS,0).EQ.0) THEN
            WRITE (iunout,*) 'NO RECOMBINATION '
          ELSE
            DO 220 IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              WRITE (iunout,*) 'RECOMBINATION NO. IRRC= ',IRRC
              WRITE (iunout,*) 'RECOMBINATION INTO SPECIES:'
              IION3=NIOPRC(IRRC)
              IF (IION3.NE.0) WRITE (iunout,*) 'TEST ION IION= ',
     .                                     TEXTS(NSPAM+IION3)
              IPLS3=NPLPRC(IRRC)
              IF (IPLS3.NE.0) WRITE (iunout,*) 'BULK ION IPLS= ',
     .                                     TEXTS(NSPAMI+IPLS3)
              IATM3=NATPRC(IRRC)
              IF (IATM3.NE.0) WRITE (iunout,*) 'ATOM     IATM= ',
     .                                     TEXTS(NSPH+IATM3)
              IMOL3=NMLPRC(IRRC)
              IF (IMOL3.NE.0) WRITE (iunout,*) 'MOLECULE IMOL= ',
     .                                     TEXTS(NSPA+IMOL3)
              IPHOT3=NPHPRC(IRRC)
              IF (IPHOT3.NE.0) WRITE (iunout,*) 'PHOTON  IPHOT= ',
     .                                     TEXTS(IPHOT3)
C  and, possibly, a second secondary
              IION3=NIOPRC_2(IRRC)
              IF (IION3.NE.0) WRITE (iunout,*) 'TEST ION IION= ',
     .                                     TEXTS(NSPAM+IION3)
              IPLS3=NPLPRC_2(IRRC)
              IF (IPLS3.NE.0) WRITE (iunout,*) 'BULK ION IPLS= ',
     .                                     TEXTS(NSPAMI+IPLS3)
              IATM3=NATPRC_2(IRRC)
              IF (IATM3.NE.0) WRITE (iunout,*) 'ATOM     IATM= ',
     .                                     TEXTS(NSPH+IATM3)
              IMOL3=NMLPRC_2(IRRC)
              IF (IMOL3.NE.0) WRITE (iunout,*) 'MOLECULE IMOL= ',
     .                                     TEXTS(NSPA+IMOL3)
              IPHOT3=NPHPRC_2(IRRC)
              IF (IPHOT3.NE.0) WRITE (iunout,*) 'PHOTON  IPHOT= ',
     .                                     TEXTS(IPHOT3)
C
C             WRITE (iunout,*) 'ELECTRONS: PELPRC,EELRC1'
C             IF (NSTORDR >= NRAD) THEN
C               WRITE (iunout,*) 'EL      ',1.,EELRC1(IRRC,1)
C             ELSE
C               WRITE (iunout,*) 'EL      ',1.,FEELRC1(IRRC,1)
C             END IF
220         CONTINUE
          ENDIF
          CALL LEER(1)
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR RECOMBINATION'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (iunout,*) 'NREC TOO SMALL, CHECK PARAMETER STATEMENTS'
      CALL EXIT_OWN(1)
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (iunout,*) 'MASS CONSERVATION VIOLATED, IPLS,IRRC ',
     .                  IPLS,IRRC
      CALL EXIT_OWN(1)
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTP.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IPLS = ',IPLS
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
C
      END
C ===== SOURCE: xsectp_param.f
C
C
      SUBROUTINE XSECTP_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: NRC, IATM, IPLS, KK
C
C
C   RECOMBINATION
C
      DO 1000 IPLS=1,NPLSI
C
        IF (NRCP(IPLS).EQ.0) THEN
C
          IF (NCHARP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
C
C  DEFAULT HYDROGENIC RECOMBINATION MODEL
C  HYDR. RECOMBINATION RATE-COEFFICIENT (1/S/CCM) E + H+ --> H + RAD.
C  GORDEEV ET. AL., PIS'MA ZH. EHKSP. TEOR. FIZ. 25 (1977) 223.
C
            DO 52 IATM=1,NATMI
              IF (NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .                            NCHRGP(IPLS).EQ.1) THEN
C
                NREC=NREC+1
              ENDIF
52          CONTINUE
C
          ENDIF
C
C  NON DEFAULT MODEL:  240--
C
        ELSEIF (NRCP(IPLS).GT.0) THEN
          DO 82 NRC=1,NRCP(IPLS)
            KK=IREACP(IPLS,NRC)
            IF ((ISWR(KK).NE.6) .AND. (ISWR(KK).NE.7)) GOTO 82
C
C  RECOMBINATION MODEL FOR BULK IONS
            NREC=NREC+1
C
82        CONTINUE
C
C  NO MODEL DEFINED
        ELSE
        ENDIF
C
C
1000  CONTINUE
C
      RETURN
C
      END
C ===== SOURCE: xsectph.f
C
C
      SUBROUTINE XSECTPH
C
C  TABLE FOR REACTION RATES FOR PHOTONS
C
      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CTRCEI
      USE PHOTON
      IMPLICIT NONE
csw
csw   PHOTON COLLISIONS, OT - type
csw
      integer :: kk,iphot,idsc,nrc,ipl0,ipl1,ipl2,ityp1,ityp2,ifnd,
     .    updf,mode, idot


      idot=0

      DO IPHOT=1,NPHOTI
        IDSC=0
        PHV_LGPHOT(IPHOT,0,0)=0
        PHV_LGPHOT(IPHOT,0,1)=0
C
C   AT PRESENT NO DEFAULT MODEL
C
        IF (NRCPH(IPHOT).EQ.0) THEN
          PHV_NPHOTI(IPHOT)=0
C
C  NON DEFAULT "OT" MODEL:
C
        ELSEIF(NRCPH(IPHOT) > 0) THEN
          DO NRC=1,NRCPH(IPHOT)
            KK=IREACPH(IPHOT,NRC)
            IF (ISWR(KK).NE.7) CYCLE
            IDSC=IDSC+1
            IDOT=IDOT+1
            NREAOT(IDOT) = KK
            CALL PH_XSECTPH (IPHOT,NRC,IDSC)
          ENDDO
          PHV_NPHOTI(IPHOT)=IDSC
C  NO "OT" MODEL DEFINED
        ELSE
          PHV_NPHOTI(IPHOT)=0
        ENDIF

CDR     PHV_NPHOTIM(IPHOT)=PHV_NPHOTI(IPHOT)-1

        PHV_LGPHOT(IPHOT,0,0)=IDSC

      ENDDO
csw
csw output:
csw
      DO IPHOT=1,NPHOTI
C
        IF (TRCAMD) THEN
          CALL MASBOX ('PHOTON SPECIES IPHOT = '//TEXTS(IPHOT))
          CALL LEER(1)
C
          IF(PHV_NPHOTI(iphot).eq.0) then
            CALL LEER(1)
            WRITE (iunout,*) 'NO "OT"-REACTION FOR THIS PHOTON'
            CALL LEER(1)
          ELSE
            DO IDSC=1,PHV_NPHOTI(IPHOT)
              CALL LEER(1)
              WRITE (iunout,*) '(OTHER) REACTION NO. IROT= ',IDSC
              CALL LEER(1)

                  ipl0=PHV_LGPHOT(iphot,idsc,1)
                  ifnd=PHV_LGPHOT(iphot,idsc,2)
                  kk=PHV_LGPHOT(iphot,idsc,3)
                  updf=PHV_LGPHOT(iphot,idsc,4)
                  mode=PHV_LGPHOT(iphot,idsc,5)

                  ityp1=PHV_N1STOTph(iphot,idsc,1)
                  ipl1= PHV_N1STOTph(iphot,idsc,2)
                  ityp2=PHV_N2NDOTph(iphot,idsc,1)
                  ipl2= PHV_N2NDOTph(iphot,idsc,2)

                  write (iunout,*) 'irot,ipl0,il,kk,updf,mode'
                  write (iunout,*)  idsc,ipl0,ifnd,kk,updf,mode
                  write (iunout,*) 'ityp1,ipl1,ityp2,ipl2'
                  write (iunout,*)  ityp1,ipl1,ityp2,ipl2
                  call leer(1)
               enddo
            endif
         endif
      enddo

      RETURN
      END SUBROUTINE XSECTPH
C ===== SOURCE: xsectph_param.f
C  28.6.05
C nnrot --> nrot, for subr. find_param, setamd, after removing phv_nrota,
C                 and phv_nrotph
C
      SUBROUTINE XSECTPH_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMUSR
      USE PHOTON
      IMPLICIT NONE
csw
csw  PHOTON COLLISION, OT - type
csw
      integer :: iphot,nrc,kk

      do iphot=1,nphoti
         if(nrcph(iphot) > 0) then
            do nrc=1,nrcph(iphot)
               kk=ireacph(iphot,nrc)
               if(iswr(kk) == 7) then
                 NROT=NROT+1
               endif
            enddo
         endif
      enddo
cdr
c  this call is still necessary because not all OT-process data
c  have already been moved to module COMXS.
C  Still some clean-up work to be done
cdr
      call PH_ALLOC_XSECTPH(nrot)

      RETURN

      END SUBROUTINE XSECTPH_PARAM
C ===== SOURCE: xstapi.f
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
C
C
      SUBROUTINE XSTAPI(PLS)
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR ATOMIC SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMPRT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: PLS(*)
      REAL(DP) :: CF(9)
      REAL(DP) :: FCTKKL, P2N, FACTKK, TMASS, ADDT, ADDTL, PMASS, 
     .            CHRDIF, COU, RATE_COEFF
      INTEGER :: NSEPI4, NSEPI5, IAPI, I, NEND, J, IAT, IO, ION, IA,
     .           IML, IM, MODC, KK, NRC, IDSC, IRPI, IIO, IPL, IPLSTI
      INTEGER, EXTERNAL :: IDEZ
      SAVE
C
C   ION IMPACT COLLISIONS
C
      DO 1000 IATM=1,NATMI
        IDSC=0
        LGAPI(IATM,0,0)=0
        LGAPI(IATM,0,1)=0
C
C  NO DEFAULT MODEL
C

        IF (NRCA(IATM).EQ.0) THEN
          NAPII(IATM)=0
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.4) GOTO 130
C
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
C  INCIDENT BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IPLSTI=MPLSTI(IPLS)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 990
            IDSC=IDSC+1
            NRPII=NRPII+1
            IF (NRPII.GT.NRPI) GOTO 999
            IRPI=NRPII
            NREAPI(IRPI) = KK
            LGAPI(IATM,IDSC,0)=IRPI
            LGAPI(IATM,IDSC,1)=IPLS
            PPLPI(IRPI,IPLS)=PPLPI(IRPI,IPLS)-1.D0
C  SECONDARY INDEX, FIRST SECONDARY
            ITYP=IDEZ(ISCD1A(IATM,NRC),1,3)
            ISPZ=IDEZ(ISCD1A(IATM,NRC),3,3)
            IF (ITYP.EQ.3) PIOPI(IRPI,ISPZ)=PIOPI(IRPI,ISPZ)+1.D0
            IF (ITYP.EQ.4) PPLPI(IRPI,ISPZ)=PPLPI(IRPI,ISPZ)+1.D0
C  SECONDARY INDEX, SECOND SECONDARY
            ITYP=IDEZ(ISCD2A(IATM,NRC),1,3)
            ISPZ=IDEZ(ISCD2A(IATM,NRC),3,3)
            IF (ITYP.EQ.3) PIOPI(IRPI,ISPZ)=PIOPI(IRPI,ISPZ)+1.D0
            IF (ITYP.EQ.4) PPLPI(IRPI,ISPZ)=PPLPI(IRPI,ISPZ)+1.D0
C
            CHRDIF=0.
            DO 133 IIO=1,NIONI
              CHRDIF=CHRDIF+PIOPI(IRPI,IIO)*NCHRGI(IIO)
133         CONTINUE
            DO 134 IPL=1,NPLSI
              CHRDIF=CHRDIF+PPLPI(IRPI,IPL)*NCHRGP(IPL)
134         CONTINUE
            PELPI(IRPI)=PELPI(IRPI)+CHRDIF
C
            PPLPI(IRPI,IPLS)=PPLPI(IRPI,IPLS)+1.D0
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
            PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT: TARGET AT REST)
            TMASS=MASST(KK)*PMASSA
C
            ADDT=PMASS/RMASSP(IPLS)
            ADDTL=LOG(ADDT)
            ADDPI(IRPI,IPLS) = ADDTL
C
C CROSS SECTION (E-LAB)
            IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
              MODCOL(4,1,IRPI)=KK
              MODCOL(4,2,IRPI)=3
              IF (FACTKK.NE.1.D0)
     .        WRITE (iunout,*) 
     .          'FREACA NOT READY FOR CROSS SECTION IN XSTAPI'
            ENDIF
C
C RATE COEFFICIENT
            MODC=IDEZ(MODCLF(KK),3,5)
            IF (MODC.GE.1.AND.MODC.LE.2) THEN
              MODCOL(4,2,IRPI)=MODC
              IF (MODC.EQ.1) NEND=1
              IF (MODC.EQ.2) NEND=NSTORDT
              IF (NSTORDR >= NRAD) THEN
                DO 142 J=1,NSBOX
                  PLS(J)=TIINL(IPLSTI,J)+ADDTL
142             CONTINUE
                IF (MODC.EQ.1) THEN
                  DO 145 J=1,NSBOX
                    COU = RATE_COEFF(KK,PLS(J),0._DP,.TRUE.,0)
                    TABPI3(IRPI,J,1)=COU*DIIN(IPLS,J)*FACTKK
145               CONTINUE
                ELSEIF (MODC.EQ.2) THEN
                  FCTKKL=LOG(FACTKK)
                  DO J=1,NSBOX
                    CALL PREP_RTCS (KK,3,1,NEND,PLS(J),CF)
                    TABPI3(IRPI,J,1:NEND)=CF(1:NEND)
                    TABPI3(IRPI,J,1)=TABPI3(IRPI,J,1)+
     .                               DIINL(IPLS,J)+FCTKKL
                  END DO
                END IF
              ELSE
                FACREA(KK) = LOG(FACTKK)
              END IF
            ENDIF
C
            DEFPI(IRPI)=LOG(CVELI2*PMASS)
            EEFPI(IRPI)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4A. BULK ION ENERGY LOSS RATE
C
C  SET ENERGY LOSS RATE OF IMPACTING ION
            NSEPI4=IDEZ(ISCDEA(IATM,NRC),4,5)
            IF (NSEPI4.EQ.0) THEN
              IF (NSTORDR >= NRAD) THEN
                DO 151 J=1,NSBOX
                  EPLPI3(IRPI,J,1)=EBULKA(IATM,NRC)
151             CONTINUE
                NELRPI(IRPI)=-1
              ELSE
                NELRPI(IRPI)=-1
                EPLPI3(IRPI,1,1)=EBULKA(IATM,NRC)
              END IF
              MODCOL(4,4,IRPI)=1
            ELSEIF (NSEPI4.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 252 J=1,NSBOX
                  EPLPI3(IRPI,J,1)=1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J)
252             CONTINUE
                NELRPI(IRPI) = -2
              ELSE
                NELRPI(IRPI) = -2
              END IF
              MODCOL(4,4,IRPI)=1
            ELSE
              WRITE (iunout,*) 'NSEPI4 ILL DEFINED IN XSTAPI '
              CALL EXIT_OWN(1)
            ENDIF
C
C  4B. BULK ELECTRON ENERGY LOSS RATE
C
C  SET NET ENERGY LOSS RATE OF ELECTRON (IF ANY INVOLVED)
            NSEPI5=IDEZ(ISCDEA(IATM,NRC),5,5)
            IF (NSEPI5.EQ.0) THEN
C             MODCOL(4,4,IRPI)=1
            ELSE
              WRITE (iunout,*) 'NSEPI5 ILL DEFINED IN XSTAPI '
              CALL EXIT_OWN(1)
            ENDIF
C
130       CONTINUE
C
          NAPII(IATM)=IDSC
C  NO MODEL DEFINED
        ELSE
          NAPII(IATM)=0
        ENDIF
C
        NAPIIM(IATM)=NAPII(IATM)-1
C
        LGAPI(IATM,0,0)=0
        DO 180 IAPI=1,NAPII(IATM)
          LGAPI(IATM,0,0)=LGAPI(IATM,0,0)+LGAPI(IATM,IAPI,0)
180     CONTINUE
C
        DO 500 IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          DO 510 IAT=1,NATMI
            IA=IAT
            PATPI(IRPI,0)=PATPI(IRPI,0)+
     +                          PATPI(IRPI,IAT)
510       CONTINUE
          DO 520 IML=1,NMOLI
            IM=NATMI+IML
            PMLPI(IRPI,0)=PMLPI(IRPI,0)+
     +                          PMLPI(IRPI,IML)
520       CONTINUE
          DO 530 ION=1,NIONI
            IO=NSPAM+ION
            PIOPI(IRPI,0)=PIOPI(IRPI,0)+
     +                          PIOPI(IRPI,ION)
530       CONTINUE
          DO 540 IPL=1,NPLSI
            PPLPI(IRPI,0)=PPLPI(IRPI,0)+
     +                          PPLPI(IRPI,IPL)
540       CONTINUE
C
C
          P2NPI(IRPI)=PATPI(IRPI,0)+PMLPI(IRPI,0)+
     .                   PIOPI(IRPI,0)
          P2N=P2NP(IRPI,NSPAMI)
          DO 550 ISPZ=1,NSPAMI
            IF (P2N.GT.0.D0)
     .      P2NP(IRPI,ISPZ)=P2NP(IRPI,ISPZ)/P2N
550       CONTINUE
500     CONTINUE
1000  CONTINUE
C
      RETURN
C
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTAPI: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ION IMPACT COLLISION'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTAPI: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'INSUFFICIENT STORAGE FOR PI: NRPI=',NRPI
      CALL EXIT_OWN(1)
      RETURN
C
      END
C ===== SOURCE: xstapi_param.f
C
C
      SUBROUTINE XSTAPI_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMPRT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: KK, NRC
C
C   ION IMPACT COLLISIONS
C
      DO 1000 IATM=1,NATMI
C
C  NO DEFAULT MODEL
C
       IF (NRCA(IATM).EQ.0) THEN
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.4) GOTO 130
            NRPI=NRPI+1
C
130       CONTINUE
C
C  NO MODEL DEFINED
        ELSE
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
      END
C ===== SOURCE: xstcx.f
C 24.11.05: chrdf0 introduced (called from xsecta, xsectm, xsecti)
C 08.08.06: error exit 991 introduced: charge conservation violation
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
! 08.01.07: pls = 0.dp, twice, preset.
C
      SUBROUTINE XSTCX(RMASS,IRCX,ISP,IPL,ISCD1,ISCD2,EBULK,
     .                 CHRDF0,ISCDE,IESTM,KK,FACTKK)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EBULK, FACTKK, CHRDF0
      INTEGER, INTENT(IN) :: IRCX, ISP, IPL, ISCD1, ISCD2, ISCDE,
     .                       IESTM, KK
      REAL(DP) :: PLS(NSTORDR), CF(9,0:9), CFF(9)
      REAL(DP) :: ADD, ADDL, RMTEST, RMBULK, FCTKKL, ADDTL, CHRDIF,
     .            ADDT, TMASS, PMASS, COU, RATE_COEFF, ENERGY_RATE_COEFF
      INTEGER :: ITYP1, ITYP2, ISPZ1, IERR, ISPZ2, IATM, IPLS, KREAD,
     .           J, NEND, MODC, NSECX4, I, IPL2, IIO2, IPLTI
      INTEGER, EXTERNAL :: IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2

      SAVE
C
C  SET NON DEFAULT CHARGE EXCHANGE COLLISION PROCESS NO. IRCX
C
      IF (IPL.LE.0.OR.IPL.GT.NPLSI) GOTO 990
      IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
      RMBULK=RMASSP(IPL)
      RMTEST=RMASS
      IPLTI=MPLSTI(IPL)
C
C  1ST SECONDARY INDEX
      N1STX(IRCX,1)=IDEZ(ISCD1,1,3)
      N1STX(IRCX,2)=IDEZ(ISCD1,3,3)
      N1STX(IRCX,3)=0
      IF (N1STX(IRCX,1).LT.4) N1STX(IRCX,3)=1
C
      IF (N1STX(IRCX,1).EQ.1) THEN
        IF (RMBULK.NE.RMASSA(N1STX(IRCX,2))) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.2) THEN
        IF (RMBULK.NE.RMASSM(N1STX(IRCX,2))) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.3) THEN
        IF (RMBULK.NE.RMASSI(N1STX(IRCX,2))) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.4) THEN
        IF (RMBULK.NE.RMASSP(N1STX(IRCX,2))) GOTO 992
      ENDIF
C  2ND SECONDARY INDEX
      N2NDX(IRCX,1)=IDEZ(ISCD2,1,3)
      N2NDX(IRCX,2)=IDEZ(ISCD2,3,3)
      N2NDX(IRCX,3)=N1STX(IRCX,3)
      IF (N2NDX(IRCX,1).LT.4) N2NDX(IRCX,3)=N2NDX(IRCX,3)+1
C
      IF (N2NDX(IRCX,1).EQ.1) THEN
        IF (RMTEST.NE.RMASSA(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
        IF (RMTEST.NE.RMASSM(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
        IF (RMTEST.NE.RMASSI(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
        IF (RMTEST.NE.RMASSP(N2NDX(IRCX,2))) GOTO 992
      ENDIF
C
      CHRDIF=CHRDF0-NCHRGP(IPL)
      IF (N1STX(IRCX,1).EQ.3) THEN
        IIO2=N1STX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGI(IIO2)
      ENDIF
      IF (N1STX(IRCX,1).EQ.4) THEN
        IPL2=N1STX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGP(IPL2)
      ENDIF
      IF (N2NDX(IRCX,1).EQ.3) THEN
        IIO2=N2NDX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGI(IIO2)
      ENDIF
      IF (N2NDX(IRCX,1).EQ.4) THEN
        IPL2=N2NDX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGP(IPL2)
      ENDIF
      IF (CHRDIF.NE.0) GOTO 991
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
      PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT; TARGET AT REST)
      TMASS=MASST(KK)*PMASSA
      ADDT=PMASS/RMASSP(IPL)
      ADDTL=LOG(ADDT)
      NREACX(IRCX) = KK
      ADDCX(IRCX,IPL) = ADDTL
C
C CROSS SECTION (E-LAB)
      IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(3,1,IRCX)=KK
        MODCOL(3,2,IRCX)=3
        IF (FACTKK.NE.1.D0)
     .  WRITE (iunout,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTCX'
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),3,5)
      IF (MODC.GE.1.AND.MODC.LE.2) THEN
        MODCOL(3,2,IRCX)=MODC
        IF (MODC.EQ.1) NEND=1
        IF (MODC.EQ.2) NEND=NSTORDT
        IF (NSTORDR >= NRAD) THEN
          PLS=0._DP
          DO 242 J=1,NSBOX
            PLS(J)=TIINL(IPLTI,J)+ADDTL
242       CONTINUE
          IF (MODC.EQ.1) THEN
            DO 245 J=1,NSBOX
              COU = RATE_COEFF(KK,PLS(J),0._DP,.TRUE.,0)
              TABCX3(IRCX,J,1)=COU*DIIN(IPL,J)*FACTKK
245         CONTINUE
          ELSEIF (MODC.EQ.2) THEN
            FCTKKL=LOG(FACTKK)
            DO J=1,NSBOX
              CALL PREP_RTCS (KK,3,1,NEND,PLS(J),CFF)
              TABCX3(IRCX,J,1:NEND) = CFF(1:NEND)
              TABCX3(IRCX,J,1)=TABCX3(IRCX,J,1)+DIINL(IPL,J)+FCTKKL
            END DO
          END IF
        ELSE
          FACREA(KK) = LOG(FACTKK)
        END IF
      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS SECTION AT LEAST?
        IF (MODCOL(3,2,IRCX).NE.3) GOTO 996
      ENDIF
      DEFCX(IRCX)=LOG(CVELI2*PMASS)
      EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  3. BULK PARTICLE MOMENTUM LOSS RATE
C
C
C  4. BULK PARTICLE ENERGY LOSS RATE
C
      NSECX4=IDEZ(ISCDE,4,5)
      IF (NSECX4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK PARTICLE = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLCX3(IRCX,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
            ENDDO
            NELRCX(IRCX) = -3
          ELSE
            NELRCX(IRCX) = -3
          END IF
        ELSE
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLCX3(IRCX,J,1)=EBULK+EDRIFT(IPL,J)
251         CONTINUE
            NELRCX(IRCX) = -2
          ELSE
            NELRCX(IRCX) = -2
            EPLCX3(IRCX,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(3,4,IRCX)=3
      ELSEIF (NSECX4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLCX3(IRCX,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
252         CONTINUE
            NELRCX(IRCX) = -3
          ELSE
            NELRCX(IRCX) = -3
          END IF
        ELSE
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTCX '
          WRITE (iunout,*) 'MODIFIED TREATMENT OF CHARGE EXCHANGE '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLCX3(IRCX,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
            NELRCX(IRCX) = -2
          ELSE
            NELRCX(IRCX) = -2
            EPLCX3(IRCX,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(3,4,IRCX)=1
C     ELSEIF (NSECX4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSECX4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (iunout,*) 
     .        'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (iunout,*) 'IN CX COLLISION IRCX= ',IRCX
            WRITE (iunout,*) 'BECAUSE NO ENERGY WEIGHTED RATE AVAILABLE'
cdr         WRITE (iunout,*) 'OOPS: '
cdr         WRITE (iunout,*) 'COLLISION ESTIMATOR NOT READY '
cdr         WRITE (iunout,*) 'CALL EXIT '
cdr         CALL EXIT_OWN(1)
          ENDIF
          IESTCX(IRCX,3)=1
          MODCOL(3,4,IRCX)=2
        ELSE
C  ION ENERGY AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELRCX(IRCX) = KREAD
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(3,4,IRCX)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          IF (NSTORDR >= NRAD) THEN
            PLS=0._DP
            DO 253 J=1,NSBOX
              PLS(J)=TIINL(IPLTI,J)+ADDTL
253         CONTINUE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              DO 254 J=1,NSBOX
                EPLCX3(IRCX,J,1)=ENERGY_RATE_COEFF(KREAD,PLS(J),0._DP,
     .                           .FALSE.,0)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              DO 257 J=1,NSBOX
                CALL PREP_RTCS (KREAD,5,1,NEND,PLS(J),CFF)
                EPLCX3(IRCX,J,1:NEND) = CFF(1:NEND)
                EPLCX3(IRCX,J,1) = EPLCX3(IRCX,J,1)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF
          ELSE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLCX3(IRCX,1,1)=ADD
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREA(KREAD) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        IERR=5
        GOTO 996
      ENDIF

C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTCX(IRCX,1)=IDEZ(IESTM,1,3)
      IESTCX(IRCX,2)=IDEZ(IESTM,2,3)
      IF (IESTCX(IRCX,3).EQ.0) IESTCX(IRCX,3)=IDEZ(IESTM,3,3)
C
      ITYP1=N1STX(IRCX,1)
      ITYP2=N2NDX(IRCX,1)
      IF (IESTCX(IRCX,1).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (iunout,*) 'IRCX = ',IRCX
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,1)=0
      ENDIF
      IF (IESTCX(IRCX,2).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (iunout,*) 'IRCX = ',IRCX
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,2)=0
      ENDIF
      IF (IESTCX(IRCX,3).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR EN.-BALANCE '
        WRITE (iunout,*) 'IRCX = ',IRCX
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,3)=0
      ENDIF
      RETURN
C
      ENTRY XSTCX_2(IRCX,IPL)
C
      CALL LEER(1)
      WRITE (iunout,*) 'CHARGE EXCHANGE REACTION NO. IRCX= ',IRCX
      CALL LEER(1)
      WRITE (iunout,*) 'CHARGE EXCHANGE WITH BULK IONS IPLS:'
      WRITE (iunout,*) '1ST AND 2ND NEXT GEN. SPECIES I2ND1, I2ND2:'
      ITYP1=N1STX(IRCX,1)
      ITYP2=N2NDX(IRCX,1)
      ISPZ1=N1STX(IRCX,2)
      ISPZ2=N2NDX(IRCX,2)
      IF (ITYP1.EQ.1) TEXTS1=TEXTS(NSPH+ISPZ1)
      IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
      IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
      IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
      IF (ITYP2.EQ.1) TEXTS2=TEXTS(NSPH+ISPZ2)
      IF (ITYP2.EQ.2) TEXTS2=TEXTS(NSPA+ISPZ2)
      IF (ITYP2.EQ.3) TEXTS2=TEXTS(NSPAM+ISPZ2)
      IF (ITYP2.EQ.4) TEXTS2=TEXTS(NSPAMI+ISPZ2)
      WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPL),'I2ND1= ',TEXTS1,
     .                    'I2ND2= ',TEXTS2
      CALL LEER(1)
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE '
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED  '
      WRITE (iunout,*) 'CHARGE CONSERVATION VIOLATED '
      WRITE (iunout,*) 'IRCX, TEST-SPECIES, BULK SPECIES ',IRCX,
     .                  TEXTS(ISP),TEXTS(NSPAMI+IPL)
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK ',KK
      CALL EXIT_OWN(1)
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'EBULK_ION .LE.0, BUT MONOENERGETIC DISTRIBUTION?'
      WRITE (iunout,*) 'CHECK ENERGY FLAG ISCDEA'
      WRITE (iunout,*) 'KK,ISCDEA ',KK,ISCDEA
      CALL EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED  '
      WRITE (iunout,*) 'NO CROSS SECTION AVAILABLE FOR NON DEFAULT CX'
      WRITE (iunout,*) 'KK ',KK
      WRITE (iunout,*) 'EITHER PROVIDE CROSS SECTION OR USE DIFFERENT '
      WRITE (iunout,*) 'POST COLLISION SAMPLING FLAG ISCDEA'
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: xstei.f
!pb  28.06.06: bug fix for NSTORAM=0 and MODC=1
!pb            FACREA=FACTKK instead of FACREA=log(FACTKK) as single 
!pb            polynomial fit is linear in it parameter
!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  22.11.06: DELPOT introduced
C
C
      SUBROUTINE XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
c   isp:  incident test particle species identifier



C
C  SET NON DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EHEAVY, CHRDF0, EELEC, FACTKK
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IREI, ISP, IFRST, ISCND, ISCDE, IESTM, KK
      REAL(DP) :: CF(9,0:9)
      REAL(DP) :: EFLAG, CHRDIF, FCTKKL, FEHVDS1, EE, TB, FEELEI1, EN,
     .          P2N, EA, EI, ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMAS,
     .          ACCMSA, ACCINA, ACCINM, ACCMSP, ACCINV, COU, RATE_COEFF, 
     .          ENERGY_RATE_COEFF, DELE
      INTEGER :: MODC, KREAD, IM, IA, IERR, J, IPL, I, IP, IRAD, IO,
     .           ION, ISPZ, III, INUM, ITYP, ISPE, ICOUNT, IAT,
     .           IMM, IIO, IAA, IML, IMIN, IMAX
      INTEGER, EXTERNAL :: IDEZ

      ITYP=IDEZ(IFRST,1,3)
      INUM=IDEZ(IFRST,2,3)
      ISPE=IDEZ(IFRST,3,3)
C
C ACCUMULATED MASS OF SECONDARIES: ACCMAS (AMU)
      ACCMAS=0.D0
      ACCMSA=0.D0
      ACCMSM=0.D0
      ACCMSI=0.D0
      ACCMSP=0.D0
      ACCINV=0.D0
      ACCINA=0.D0
      ACCINM=0.D0
      ACCINI=0.D0
      ACCINP=0.D0
C
      ICOUNT=1
 85   CONTINUE
      IF (ITYP.EQ.1) THEN
        IAT=ISPE
        IAA=NSPH+IAT
        PATDS(IREI,IAT)=PATDS(IREI,IAT)+INUM
        P2ND(IREI,IAA)=P2ND(IREI,IAA)+INUM
        ACCMAS=ACCMAS+INUM*RMASSA(IAT)
        ACCMSA=ACCMSA+INUM*RMASSA(IAT)
        ACCINV=ACCINV+INUM/RMASSA(IAT)
        ACCINA=ACCINA+INUM/RMASSA(IAT)
        EATDS(IREI,IAT,1)=RMASSA(IAT)
        EATDS(IREI,IAT,2)=1./RMASSA(IAT)
      ELSEIF (ITYP.EQ.2) THEN
        IML=ISPE
        IMM=NSPA+IML
        PMLDS(IREI,IML)=PMLDS(IREI,IML)+INUM
        P2ND(IREI,IMM)=P2ND(IREI,IMM)+INUM
        ACCMAS=ACCMAS+INUM*RMASSM(IML)
        ACCMSM=ACCMSM+INUM*RMASSM(IML)
        ACCINV=ACCINV+INUM/RMASSM(IML)
        ACCINM=ACCINM+INUM/RMASSM(IML)
        EMLDS(IREI,IML,1)=RMASSM(IML)
        EMLDS(IREI,IML,2)=1./RMASSM(IML)
      ELSEIF (ITYP.EQ.3) THEN
        IIO=ISPE
        III=NSPAM+IIO
        PIODS(IREI,IIO)=PIODS(IREI,IIO)+INUM
        P2ND(IREI,III)=P2ND(IREI,III)+INUM
        ACCMAS=ACCMAS+INUM*RMASSI(IIO)
        ACCMSI=ACCMSI+INUM*RMASSI(IIO)
        ACCINV=ACCINV+INUM/RMASSI(IIO)
        ACCINI=ACCINI+INUM/RMASSI(IIO)
        EIODS(IREI,IIO,1)=RMASSI(IIO)
        EIODS(IREI,IIO,2)=1./RMASSI(IIO)
      ELSEIF (ITYP.EQ.4) THEN
        IPL=ISPE
        PPLDS(IREI,IPL)=PPLDS(IREI,IPL)+INUM
        ACCMAS=ACCMAS+INUM*RMASSP(IPL)
        ACCMSP=ACCMSP+INUM*RMASSP(IPL)
        ACCINV=ACCINV+INUM/RMASSP(IPL)
        ACCINP=ACCINP+INUM/RMASSP(IPL)
      ENDIF
C
      IF (ISCND.NE.0.AND.ICOUNT.EQ.1) THEN
        ITYP=IDEZ(ISCND,1,3)
        INUM=IDEZ(ISCND,2,3)
        ISPE=IDEZ(ISCND,3,3)
        ICOUNT=2
        GOTO 85
      ENDIF
C
      IF (ABS(ACCMAS-RMASS).GT.1.D-10) THEN
        WRITE (IUNOUT,*) 'MESSAGE FROM XSTEI.F: '
        WRITE (IUNOUT,*) 'FOR INCIDENT SPECIES ',TEXTS(ISP)
        WRITE (iunout,*) 'MASS CONSERVATION VIOLATED FOR REACT. KK= ',KK
c slmod begin - debug
        WRITE (iunout,*) '  ITYP =',ityp
        WRITE (iunout,*) '  IAT  =',iat
        WRITE (iunout,*) '  IML  =',iml
        WRITE (iunout,*) '  IIO  =',iio
        WRITE (iunout,*) '  IPL  =',ipl
c slmod end
        CALL EXIT_OWN(1)
      ENDIF
C
      DO IAT=1,NATMI
        EATDS(IREI,IAT,1)=EATDS(IREI,IAT,1)/ACCMAS
        EATDS(IREI,IAT,2)=EATDS(IREI,IAT,2)/ACCINV
      ENDDO
      EATDS(IREI,0,1)=ACCMSA/ACCMAS
      EATDS(IREI,0,2)=ACCINA/ACCINV
      DO IML=1,NMOLI
        EMLDS(IREI,IML,1)=EMLDS(IREI,IML,1)/ACCMAS
        EMLDS(IREI,IML,2)=EMLDS(IREI,IML,2)/ACCINV
      ENDDO
      EMLDS(IREI,0,1)=ACCMSM/ACCMAS
      EMLDS(IREI,0,2)=ACCINM/ACCINV
      DO IIO=1,NIONI
        EIODS(IREI,IIO,1)=EIODS(IREI,IIO,1)/ACCMAS
        EIODS(IREI,IIO,2)=EIODS(IREI,IIO,2)/ACCINV
      ENDDO
      EIODS(IREI,0,1)=ACCMSI/ACCMAS
      EIODS(IREI,0,2)=ACCINI/ACCINV
      EPLDS(IREI,1)=ACCMSP/ACCMAS
      EPLDS(IREI,2)=ACCINP/ACCINV
C
      CHRDIF=CHRDF0
      DO 83 IIO=1,NIONI
        CHRDIF=CHRDIF+PIODS(IREI,IIO)*NCHRGI(IIO)
83    CONTINUE
      DO 84 IPL=1,NPLSI
        CHRDIF=CHRDIF+PPLDS(IREI,IPL)*NCHRGP(IPL)
84    CONTINUE
      PELDS(IREI)=PELDS(IREI)+CHRDIF
C
C
C  1.) CROSS SECTION(TE) : NOT NEEDED
C
C
C  2.) RATE COEFFICIENT (CM**3/S) * ELECTRON DENSITY (CM**-3)
C
C
C  2.A) RATE COEFFICIENT = CONST.
C     TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
      IF (IDEZ(MODCLF(KK),3,5).EQ.1) THEN
        IF (NSTORDR >= NRAD) THEN
C  RATE:  (1/S)
C  RATE COEFFICIENT: (CM^3/S)
          DO J=1,NSBOX
            COU = RATE_COEFF(KK,TEINL(J),0._DP,.TRUE.,0)
            TABDS1(IREI,J)=COU*FACTKK
            IF (IFTFLG(KK,2) < 100) 
     .        TABDS1(IREI,J)=TABDS1(IREI,J)*DEIN(J)
          END DO
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ELSE
          FACREA(KK) = FACTKK
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ENDIF
        MODCOL(1,2,IREI)=1
C     ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C  TO BE WRITTEN
C       MODCOL(1,2,IREI)=2
      ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        IF (NSTORDR >= NRAD) THEN
          FCTKKL=LOG(FACTKK)
          DO J=1,NSBOX
            COU = RATE_COEFF(KK,TEINL(J),PLS(J),.FALSE.,1)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DEINL(J)
            TB=MAX(-100._DP,TB)
            TABDS1(IREI,J)=EXP(TB)
          END DO
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ELSE
          FACREA(KK) = FACTKK
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ENDIF
        MODCOL(1,2,IREI)=1
      ELSE
        IERR=1
        GOTO 996
      ENDIF
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ENERGY RATES
C
C  4.A: ELECTRON ENERGY LOSS RATES
C
      EFLAG=IDEZ(ISCDE,5,5)
      IF (EFLAG.EQ.0) THEN
C  4.A1) ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 101 J=1,NSBOX
                  EELDS1(IREI,J)=EELEC
101             CONTINUE
                NELREI(IREI)=-2
              ELSE
                NELREI(IREI)=-2
                EELDS1(IREI,1)=EELEC
              END IF
              MODCOL(1,4,IREI)=1
      ELSEIF (EFLAG.EQ.1) THEN
              IF (NSTORDR >= NRAD) THEN
                DO 103 J=1,NSBOX
                  EELDS1(IREI,J)=-1.5*TEIN(J)
103             CONTINUE
                NELREI(IREI)=-3
              ELSE
                NELREI(IREI)=-3
              END IF
              MODCOL(1,4,IREI)=1
      ELSEIF (EFLAG.EQ.3) THEN
C  4.A2) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELEC
                MODC=IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    DO 102 J=1,NSBOX
                      EELDS1(IREI,J)=-ENERGY_RATE_COEFF(KREAD,TEINL(J),
     .                                0._DP,.TRUE.,0)*DEIN(J)*FACTKK/
     .                                (TABDS1(IREI,J)+EPS60) 
102                 CONTINUE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                    FACREA(KREAD)=FACTKK
                  ENDIF
                  MODCOL(1,4,IREI)=1
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    FCTKKL=LOG(FACTKK)
                    DO J = 1, NSBOX
                      EE = ENERGY_RATE_COEFF(KREAD,TEINL(J),PLS(J),
     .                                       .FALSE.,1)
                      EE = MAX(-100._DP,EE+FCTKKL+DEINL(J))
                      EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
                    END DO
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                    FACREA(KREAD)=LOG(FACTKK)
                  END IF
                  MODCOL(1,4,IREI)=1
                ENDIF
                IF (DELPOT(KREAD).NE.0.D0) THEN
                  DELE=DELPOT(KREAD)
                  IF (NSTORDR >= NRAD) THEN
                    DO J=1,NSBOX
                      EELDS1(IREI,J)=EELDS1(IREI,J)+DELE
                    END DO  
                  END IF
                ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  4.B: HEAVY PARTICLE ENERGY GAIN RATE
C
      EFLAG=IDEZ(ISCDE,3,5)
      IF (EFLAG.EQ.0) THEN
C  4.B1)  RATE = CONST.*RATECOEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 201 J=1,NSBOX
            EHVDS1(IREI,J)=EHEAVY
201       CONTINUE
          NREAHV(IREI)=-1
        ELSE
          NREAHV(IREI)=-1
          EHVDS1(IREI,1)=EHEAVY
        END IF
C     ELSEIF (EFLAG.EQ.1) THEN
C        NOT A VALID OPTION
      ELSEIF (EFLAG.EQ.3) THEN
C  4.B3)  ENERGY RATE = EN.WEIGHTED RATE(TE)
        KREAD=EHEAVY
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 202 J=1,NSBOX
              EHVDS1(IREI,J)=ENERGY_RATE_COEFF(KREAD,TEINL(J),0._DP,
     .                .TRUE.,0)*DEIN(J)*FACTKK/(TABDS1(IREI,J)+EPS60)
202         CONTINUE
            NREAHV(IREI)=KREAD
          ELSE
            NREAHV(IREI)=KREAD
            FACREA(KREAD)=LOG(FACTKK)
          END IF
        ELSE
          WRITE (iunout,*) 'INVALID OPTION IN XSTEI '
          CALL EXIT_OWN(1)
        ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEI(IREI,1)=IDEZ(IESTM,1,3)
      IESTEI(IREI,2)=IDEZ(IESTM,2,3)
      IESTEI(IREI,3)=IDEZ(IESTM,3,3)
C
      IF (IESTEI(IREI,1).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,1)=0
      ENDIF
      IF (IESTEI(IREI,2).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,2)=0
      ENDIF
      CALL LEER(1)
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTEI_1(IREI)
C
          DO 510 IAT=1,NATMI
            IA=NSPH+IAT
            PATDS(IREI,0)=PATDS(IREI,0)+
     +                          PATDS(IREI,IAT)
            P2ND(IREI,IA)=P2ND(IREI,IA-1)+
     +                          P2ND(IREI,IA)
510       CONTINUE
          DO 520 IML=1,NMOLI
            IM=NSPA+IML
            PMLDS(IREI,0)=PMLDS(IREI,0)+
     +                          PMLDS(IREI,IML)
            P2ND(IREI,IM)=P2ND(IREI,IM-1)+
     +                          P2ND(IREI,IM)
520       CONTINUE
          DO 530 ION=1,NIONI
            IO=NSPAM+ION
            PIODS(IREI,0)=PIODS(IREI,0)+
     +                          PIODS(IREI,ION)
            P2ND(IREI,IO)=P2ND(IREI,IO-1)+
     +                          P2ND(IREI,IO)
530       CONTINUE
          DO 540 IPL=1,NPLSI
            PPLDS(IREI,0)=PPLDS(IREI,0)+
     +                          PPLDS(IREI,IPL)
540       CONTINUE
C
          P2NDS(IREI)=PATDS(IREI,0)+PMLDS(IREI,0)+
     .                   PIODS(IREI,0)
          P2N=P2ND(IREI,NSPAMI)
          DO 550 ISPZ=NSPH+1,NSPAMI
            IF (P2N.GT.0.D0)
     .      P2ND(IREI,ISPZ)=P2ND(IREI,ISPZ)/P2N
550       CONTINUE
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTEI_2(IREI)
C
      CALL LEER(1)
      WRITE (iunout,*) 'ELEC. IMPACT REACTION NO. IREI= ',IREI
      CALL LEER(1)
      EI=1.D30
      EA=-1.D30
      imin=0
      imax=0
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELDS1(IREI,IRAD)
        ELSE
          EN=FEELEI1(IREI,IRAD)
        END IF
        if (en < ei) imin=irad
        if (en > ea) imax=irad
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
875   CONTINUE
      WRITE (iunout,*) 'BACKGROUND SECONDARIES:'
      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (iunout,*) 'ELECTRONS: PELDS, CONSTANT ENERGY: EEL'
        WRITE (iunout,'(1X,A8,2(1PE12.4))') 'EL      ',PELDS(IREI),EI
      ELSE
        WRITE (iunout,*) 
     .    'ELECTRONS: PELDS, ENERGY RANGE: EEL_MIN,EEL_MAX'
        WRITE (iunout,'(1X,A8,3(1PE12.4))') 'EL      ',PELDS(IREI),EI,EA
      ENDIF
      write (iunout,*) ' imin = ', imin, ' imax = ',imax
C
      EI=1.D30
      EA=-1.D30
      DO 876 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 876
        IF (NSTORDR >= NRAD) THEN
          EN=EHVDS1(IREI,IRAD)
        ELSE
          EN=FEHVDS1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
876   CONTINUE
      IF (PPLDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'BULK IONS: PPLDS '
        DO 874 IPL=1,NPLSI
          IP=NSPAMI+IPL
          IF (PPLDS(IREI,IPL).NE.0.D0)
     .      WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLDS(IREI,IPL)
874     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EPLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EPLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2),
     .                                 ' * EHEAVY '
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (iunout,*) 'BULK IONS: NONE'
      ENDIF
      CALL LEER(1)
C
      WRITE (iunout,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NDS(IREI).EQ.0.D0) THEN
        WRITE (iunout,*) 'NONE'
        CALL LEER(1)
        RETURN
      ENDIF
C
      IF (PATDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'ATOMS    : PATDS '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATDS(IREI,IAT).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IA),PATDS(IREI,IAT)
871     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EATDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EATDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'MOLECULES: PMLDS '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLDS(IREI,IML).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IM),PMLDS(IREI,IML)
872     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EMLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EMLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIODS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'TEST IONS: PIODS '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIODS(IREI,IIO).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IO),PIODS(IREI,IIO)
873     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EIODS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EIODS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      RETURN
C
C
C-----------------------------------------------------------------------
C
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI, MODCLF(KK) ',MODCLF(KK)
      WRITE (iunout,*) IREI,KK
      CALL EXIT_OWN(1)
997   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: ISCDE FLAG'
      WRITE (iunout,*) IREI
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: xstel.f
!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
cdr  05.01.07:  write(6,...) --> write(iunout,...) in one place
C
C
      SUBROUTINE XSTEL(IREL,ISP,IPL,
     .                 EBULK,ISCDE,IESTM,
     .                 KK,FACTKK)
C
C  RETURNS:
C    MODCOL(5,...)
C    TABEL3(IREL,NCELL,...)
C    EPLEL3(IREL,NCELL,...)
C    DEFEL(IREL)
C    EEFEL(IREL)
C    IESTEL(IREL,...)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: EBULK, FACTKK
      INTEGER, INTENT(IN) :: IREL, ISP, IPL, ISCDE, IESTM, KK
      REAL(DP) :: PLS(NSTORDR), CF(9,0:9), CFF(9)
      REAL(DP) :: FCTKKL, ADD, ADDL, ADDT, ADDTL, PMASS, TMASS, COU,
     .            RATE_COEFF, ENERGY_RATE_COEFF
      INTEGER :: I, NSEEL4, NEND, J, KREAD, MODC, IDEZ, IERR, IPLTI

      SAVE
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
      PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT; TARGET AT REST)
      TMASS=MASST(KK)*PMASSA
      ADDT=PMASS/RMASSP(IPL)
      ADDTL=LOG(ADDT)
      ADDEL(IREL,IPL) = ADDTL
      NREAEL(IREL) = KK
      IPLTI = MPLSTI(IPL)
C
C POTENTIAL
      IF (IDEZ(MODCLF(KK),1,5).EQ.1) THEN
        MODCOL(5,0,IREL)=KK
      ENDIF
C
C CROSS SECTION (E-LAB), IN FUNCTION CROSS, K=KK
      IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(5,1,IREL)=KK
        MODCOL(5,2,IREL)=3
        IF (FACTKK.NE.1.D0)
     .    WRITE (iunout,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTEL'
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),3,5)
      IF (MODC.GE.1.AND.MODC.LE.2) THEN
        MODCOL(5,2,IREL)=MODC
        IF (MODC.EQ.1) NEND=1
        IF (MODC.EQ.2) NEND=NSTORDT
        IF (NSTORDR >= NRAD) THEN
          DO 242 J=1,NSBOX
            PLS(J)=TIINL(IPLTI,J)+ADDTL
242       CONTINUE
          IF (MODC.EQ.1) THEN
            DO 245 J=1,NSBOX
              COU = RATE_COEFF(KK,PLS(J),0._DP,.TRUE.,0)
              TABEL3(IREL,J,1)=COU*DIIN(IPL,J)*FACTKK
245         CONTINUE
          ELSEIF (MODC.EQ.2) THEN
            FCTKKL=LOG(FACTKK)
            DO J=1,NSBOX
              CALL PREP_RTCS(KK,3,1,NEND,PLS(J),CFF)
              TABEL3(IREL,J,1:NEND) = CFF(1:NEND) 
              TABEL3(IREL,J,1)=TABEL3(IREL,J,1)+DIINL(IPL,J)+FCTKKL
            END DO
          END IF
        ELSE
          FACREA(KK) = LOG(FACTKK)
        END IF
      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS SECTION AT LEAST?
        IF (MODCOL(5,2,IREL).NE.3) GOTO 993
      ENDIF
      DEFEL(IREL)=LOG(CVELI2*PMASS)
      EEFEL(IREL)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4. BULK ION ENERGY LOSS RATE
C
      NSEEL4=IDEZ(ISCDE,4,5)
      IF (NSEEL4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK ION = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
            ENDDO
            NELREL(IREL) = -3
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK
251         CONTINUE
            NELREL(IREL) = -1
          ELSE
            NELREL(IREL) = -1
            EPLEL3(IREL,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(5,4,IREL)=3
      ELSEIF (NSEEL4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
252         CONTINUE
            NELREL(IREL) = -3
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTEL '
          WRITE (iunout,*) 'MODIFIED TREATMENT OF ELASTIC COLLISIONS '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
            NELREL(IREL) = -2
          ELSE
            NELREL(IREL) = -2
            EPLEL3(IREL,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(5,4,IREL)=1
C     ELSEIF (NSEEL4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEEL4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (iunout,*) 
     .        'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (iunout,*) 'IN ELASTIC COLLISION IREL= ',IREL
            WRITE (iunout,*) 'BECAUSE NO ENERGY WEIGHTED RATE AVAILABLE'
          ENDIF
          IESTEL(IREL,3)=1
          MODCOL(5,4,IREL)=2
        ELSE  ! ION ENERGY AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELREL(IREL) = KREAD
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(5,4,IREL)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          IF (NSTORDR >= NRAD) THEN
            DO 253 J=1,NSBOX
              PLS(J)=TIINL(IPLTI,J)+ADDTL
253         CONTINUE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              DO 254 J=1,NSBOX
                EPLEL3(IREL,J,1)=ENERGY_RATE_COEFF(KREAD,PLS(J),0._DP,
     .                           .FALSE.,0)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              DO 257 J=1,NSBOX
                CALL PREP_RTCS(KK,5,1,NEND,PLS(J),CFF)
                EPLEL3(IREL,J,1:NEND) = CFF(1:NEND) 
                EPLEL3(IREL,J,1) = EPLEL3(IREL,J,1)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF
          ELSE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLEL3(IREL,1,1)=ADD
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREA(KREAD) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        IERR=5
        GOTO 993
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEL(IREL,1)=IDEZ(IESTM,1,3)
      IESTEL(IREL,2)=IDEZ(IESTM,2,3)
      IF (IESTEL(IREL,3).EQ.0) IESTEL(IREL,3)=IDEZ(IESTM,3,3)
C
C
C
      IF (IESTEL(IREL,2).EQ.0.AND.NPBGKP(IPL,1).EQ.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: TR.L.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (iunout,*) 'IREL = ',IREL
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO COLLISION ESTIMATOR '
        IESTEL(IREL,2)=1
      ENDIF
      IF (IESTEL(IREL,3).EQ.0.AND.NPBGKP(IPL,1).EQ.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: TR.L.EST NOT AVAILABLE FOR EN.-BALANCE '
        WRITE (iunout,*) 'IREL = ',IREL
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO COLLISION ESTIMATOR '
        IESTEL(IREL,3)=1
      ENDIF
      RETURN
C
      ENTRY XSTEL_2(IREL,IPL)
C
      CALL LEER(1)
      WRITE (iunout,*) 'ELASTIC COLLISION NO. IREL= ',IREL
      CALL LEER(1)
      WRITE (iunout,*) 'ELASTIC COLLISION WITH BULK IONS IPLS:'
      WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPL)
C
      CALL LEER(1)
      RETURN
C
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEL, SPECIES ISP: '
      WRITE (iunout,*) ISP,IREL
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: xstrc.f
!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised

      SUBROUTINE XSTRC(ipls,nrc,idsc,irrc)
cdr
cdr  to replace ph_xsectp, as called from XSECTP.F,
cdr  prepare volume recombination processes: bulk (+ bulk)--> test (+ bulk)
cdr  i.e.                                    H+    +  e   --> H    (+ rad.)
cdr  i.e.                                    H(n=2)       --> Ly-alpha (+H(n=1))
cdr
c    ipls: incident bulk
c    nrc : index of reaction in list of all reactions for IPLS
c    idsc: index of RC reaction for species ipls
c    irrc: index of RC reaction in NRRC arrays
c
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE COMXS
      USE PHOTON
      IMPLICIT NONE
      integer, intent(in) :: ipls,nrc,idsc,irrc
      integer :: kk,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,
     .           j
      integer, external :: idez
      real(dp) :: factkk,aik

c  fetch data for process nrc

      kk = IREACP(ipls,nrc)
      if (kk /= idreac) call get_reaction(kk)

      FACTKK=FREACP(IPLS,NRC)
      IF (FACTKK.EQ.0.D0) FACTKK=1.
      aik=reaction%aik


      IPL0 =IDEZ(IBULKP(ipls,nrc),3,3)
      IPL1 =IDEZ(ISCD1P(ipls,nrc),3,3)
      IPL2 =IDEZ(ISCD2P(ipls,nrc),3,3)
      ITYP0=IDEZ(IBULKP(ipls,nrc),1,3)
      ITYP1=IDEZ(ISCD1P(ipls,nrc),1,3)
      ITYP2=IDEZ(ISCD2P(ipls,nrc),1,3)

      LGPRC(IPLS,IDSC)=IRRC

      facrea(kk) = factkk
      NREARC(irrc) = kk
      do j=1,nrad
         tabrc1(irrc,j)=aik*factkk
      enddo
      modcol(6,2,irrc)=1

      select case(ityp1)
      case(0)
         NPHPRC(IRRC)=IPL1
      case(1)
         NATPRC(IRRC)=IPL1
      case(2)
         NMLPRC(IRRC)=IPL1
      case(3)
         NIOPRC(IRRC)=IPL1
      case(4)
         NPLPRC(IRRC)=IPL1
      case default
         write(iunout,*) 'volume-processes: xstrc.f:'
         write(iunout,*) '   ityp1=',ityp1,' not allowed'
         call exit(1)
      end select
      select case(ityp2)
      case(0)
         NPHPRC_2(IRRC)=IPL2
      case(1)
         NATPRC_2(IRRC)=IPL2
      case(2)
         NMLPRC_2(IRRC)=IPL2
      case(3)
         NIOPRC_2(IRRC)=IPL2
      case(4)
         NPLPRC_2(IRRC)=IPL2
      case default
         write(iunout,*) 'volume-processes: xstrc.f:'
         write(iunout,*) '   ityp2=',ityp2,' not allowed'
         call exit(1)
      end select
      return
      END SUBROUTINE XSTRC
