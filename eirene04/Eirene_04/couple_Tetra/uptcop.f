C
C
      SUBROUTINE UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CSPEZ
      USE CZT1
      USE CGEOM
      USE COMPRT
      USE CSDVI
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                      XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: P, WTRSIG, EION, V0_PARB, PARMOM_0, DIST, WTR
      INTEGER :: IAEL, IREL, IPL2, IAEI, IRDS, IBGK, IICX, IIEI, IIEL,
     .           IMEL, IPL1, I, IPL, IIO, IRD, IP, IR, IML, IAT, IFIRST,
     .           IRCX, IADD, ICOU, IACX, IRDD, IMCX, IMEI, IPLSTI,
     .           IPLSV, IPLV
      INTEGER, SAVE :: NMTSP
      REAL(DP), ALLOCATABLE, SAVE ::
     . CNDYNA(:), CNDYNM(:), CNDYNI(:)
CDR
      REAL(DP), ALLOCATABLE, SAVE ::
     . VPX(:),    VPY(:),    VRX(:),    VRY(:)
CDR
      DATA IFIRST/0/
      SAVE
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNA(NATM))
        ALLOCATE (CNDYNM(NMOL))
        ALLOCATE (CNDYNI(NION))
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        DO 11 IAT=1,NATMI
11        CNDYNA(IAT)=AMUA*RMASSA(IAT)
        DO 12 IML=1,NMOLI
12        CNDYNM(IML)=AMUA*RMASSM(IML)
        DO 13 IIO=1,NIONI
13        CNDYNI(IIO)=AMUA*RMASSI(IIO)
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
C
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        DO 1 I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
1       CONTINUE
        DO 2 IR=1,NR1STM
          DO 2 IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
2       CONTINUE
C
        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI
C
      ENDIF
C
C  WV=WEIGHT/VEL
C
C  ATOMS
      IF (ITYP.EQ.1) THEN
        DO 20 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 20
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C  1,NPLSI:
C              PARTICLE CHARGE EXCHANGE RATE DUE TO IPLS: #/S
C              WITH ATOM SPECIES IATM=1,NATMI, PER ION
C  EACH RATE IS WEIGHTED WITH THE FACTOR (E0/EI-1), E0 BEING
C  THE NEUTRAL PARTICLE ENERGY, EI THE MEAN PLASMA ION ENERGY
C  THESE RATES ARE SCALED IN THE SHORT CYCLE WITH EI*NI
C
C
          IF (NCPVI.LT.NPLSI) GOTO 20
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 51
          DO 52 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 52
            IPLSTI = MPLSTI(IPLS)
            EION=1.5*TIIN(IPLSTI,IRD)+EDRIFT(IPLS,IRD)
            WTRSIG=WTR*SIGVCX(IRCX)/DIIN(IPLS,IRD)
            COPV(IPLS,IRDD)=COPV(IPLS,IRDD)+WTRSIG*(E0/EION-1.)
            LMETSP(NMTSP+IPLS)=.TRUE.
52        CONTINUE
51        CONTINUE
C
C.........................................
C
C   MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C.........................................
C
C
C  CONTRIBUTIONS FROM ATOMS
C  NPLSI+1, 2*NPLSI:
C
          IF (NCPVI.LT.2*NPLSI) GOTO 20
C
          IADD=NPLSI
          V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD))
          PARMOM_0=V0_PARB*CNDYNA(IATM)
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 59
          DO 56 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 56
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 56
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHA
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM(IPL,IRD)
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                        SIGN(1._DP,BVIN(IPLV,IRD))
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
56        CONTINUE
59        CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 61 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 62 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IRDS)*P
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  IPLV=MPLSV(IPL)
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                            SIGN(1._DP,BVIN(IPLV,IRD))
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
62            CONTINUE
            ENDIF
61        CONTINUE
C
C  ION IMPACT IONIZATION CONTRIBUTION: NOT INCLUDED
C
C
C  ELASTIC CONTRIBUTION FROM ATOMS
C
C
          IF (LGAEL(IATM,0,0).EQ.0) GOTO 80
C  DEFAULT TRACKLENGTH ESTIMATOR (BGK APPROXIMATION)
          DO 81  IAEL=1,NAELI(IATM)
            IREL=LGAEL(IATM,IAEL,0)
            IPLS=LGAEL(IATM,IAEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 81
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 81
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB*
     .                      SIGN(1._DP,BVIN(IPLSV,IRD))
            LMETSP(NMTSP+IPL2)=.TRUE.
81        CONTINUE
80      CONTINUE
C
20      CONTINUE
C
C  MOLECULES
      ELSEIF (ITYP.EQ.2) THEN
C
        DO 200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 200
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C
C  CONTRIBUTIONS FROM MOLECULES
C  2*NPLSI+1, 3*NPLSI:
C
          IF (NCPVI.LT.3*NPLSI) GOTO 200
C
          IADD=2*NPLSI
          V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD))
          PARMOM_0=V0_PARB*CNDYNM(IMOL)
C
          IF (LGMCX(IMOL,0,0).EQ.0) GOTO 590
          DO 560 IMCX=1,NMCXI(IMOL)
            IRCX=LGMCX(IMOL,IMCX,0)
            IPLS=LGMCX(IMOL,IMCX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 560
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 560
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHM
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM(IPL,IRD)
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                        SIGN(1._DP,BVIN(IPLV,IRD))
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
560       CONTINUE
590       CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 610 IMEI=1,NMDSI(IMOL)
            IRDS=LGMEI(IMOL,IMEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 620 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IRDS)*P
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  IPLV=MPLSV(IPL)
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                            SIGN(1._DP,BVIN(IPLV,IRD))
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
620           CONTINUE
            ENDIF
610       CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM MOLECULES
C
C
          IF (LGMEL(IMOL,0,0).EQ.0) GOTO 800
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 810 IMEL=1,NMELI(IMOL)
            IREL=LGMEL(IMOL,IMEL,0)
            IPLS=LGMEL(IMOL,IMEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 810
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 810
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLSV,IRD))
            LMETSP(NMTSP+IPL2)=.TRUE.
810       CONTINUE
800     CONTINUE
C
C
200     CONTINUE
C
C  TEST IONS
C
      ELSEIF (ITYP.EQ.3) THEN
C
        DO 2000 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 2000
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C
C  CONTRIBUTIONS FROM TEST IONS
C  3*NPLSI+1, 4*NPLSI:
C
          IF (NCPVI.LT.4*NPLSI) GOTO 2000
C
          IADD=3*NPLSI
          V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD))
          PARMOM_0=V0_PARB*CNDYNI(IION)
C
          IF (LGICX(IION,0,0).EQ.0) GOTO 5900
          DO 5600 IICX=1,NICXI(IION)
            IRCX=LGICX(IION,IICX,0)
            IPLS=LGICX(IION,IICX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 5600
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 5600
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHI
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
C
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM(IPL,IRD)
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                        SIGN(1._DP,BVIN(IPLV,IRD))
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
5600      CONTINUE
5900      CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 6100 IIEI=1,NIDSI(IION)
            IRDS=LGIEI(IION,IIEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 6200 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IRDS)*P
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  IPLV=MPLSV(IPL)
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                            SIGN(1._DP,BVIN(IPLV,IRD))
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
6200          CONTINUE
            ENDIF
6100      CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM TEST IONS
C
          IF (LGIEL(IION,0,0).EQ.0) GOTO 8000
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 8100 IIEL=1,NIELI(IION)
            IREL=LGIEL(IION,IIEL,0)
            IPLS=LGIEL(IION,IIEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 8100
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 8100
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLSV,IRD))
            LMETSP(NMTSP+IPL2)=.TRUE.
8100      CONTINUE
8000    CONTINUE
C
C
2000    CONTINUE
C
C
      ENDIF
C
      RETURN
      END
