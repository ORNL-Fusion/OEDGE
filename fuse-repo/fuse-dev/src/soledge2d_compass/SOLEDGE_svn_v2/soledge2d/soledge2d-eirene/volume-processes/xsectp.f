C  aug. 05:  corrected electron energy loss rate for default rec. rate
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
! 25.03.07: check of mass conservation only for up to two secondaries
C
      SUBROUTINE EIRENE_XSECTP
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR BULK ION SPECIES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_PHOTON
 
      IMPLICIT NONE
C
      REAL(DP) :: PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: DELE, FCTKKL, EEMX, ZX, DSUB, DEIMIN, RMASS2, FACTKK,
     .            RMASS2_2, CORSUM, COU, EIRENE_RATE_COEFF, 
     .            EIRENE_ENERGY_RATE_COEFF,
     .            BREMS, Z, eirene_ngffmh, ERATE
      INTEGER :: IIRC, IION3, IPLS3, IATM3, IMOL3, KK, NRC, IATM,
     .           IRRC, J, IDSC, IPLS, NSERC5, KREAD, I, MODC, IATM1,
     .           ITYP, ISPZ, ITYP2, ISPZ2, IPHOT3
      INTEGER, EXTERNAL :: EIRENE_IDEZ
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
c test yannick
          	  NREARC(IRRC) = IRRC        	
c                  NREARC(IRRC) = 0
                  JEREARC(IRRC) = 0
                  NELRRC(IRRC) = -1
                ELSE
c test yannick
c                  NREARC(IRRC) = 0
        	  NREARC(IRRC) = IRRC
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
            if(iswr(kk)==7) then    ! OT Processes
               idsc=idsc+1
               nrrci=nrrci+1
               IF (NRRCI.GT.NREC) GOTO 992
               call EIRENE_XSTRC(ipls,nrc,idsc,nrrci)
               cycle
csw end branch
            ELSEIF (ISWR(KK).EQ.6) THEN  !  RC Processes
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
              ITYP=EIRENE_IDEZ(ISCD1P(IPLS,NRC),1,3)
              ISPZ=EIRENE_IDEZ(ISCD1P(IPLS,NRC),3,3)
              IF ((ISPZ < 1) .OR. (ISPZ > MAXSPC(ITYP))) GOTO 995
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
 
              ITYP2=EIRENE_IDEZ(ISCD2P(IPLS,NRC),1,3)
              ISPZ2=EIRENE_IDEZ(ISCD2P(IPLS,NRC),3,3)
              IF ((ISCD2P(IPLS,NRC) /= 0) .AND.
     .           ((ISPZ2 < 1) .OR. (ISPZ2 > MAXSPC(ITYP2)))) GOTO 995
              IF (ITYP2.EQ.4) THEN
                NPLPRC_2(IRRC)=ISPZ2
                RMASS2_2      =RMASSP(ISPZ2)
              ELSE
                RMASS2_2=0._DP
              ENDIF
C  CHECK MASS CONSERVATION
              IF (REACDAT(KK)%NOSEC < 3) THEN
                IF (RMASSP(IPLS).NE.(RMASS2+RMASS2_2)) GOTO 993
              END IF
C
C  1.) CROSS SECTION(TE)
C           NOT NEEDED
C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3) --> RATE (1/S)
C
C  2.A) RATE COEFFICIENT = CONST.
C           TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
              IF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.1) THEN
                IF (NSTORDR >= NRAD) THEN
                  LEXP = .NOT. (MOD(IFTFLG(KK,2),100) == 10)
                  DO J=1,NSBOX
!pb                    IF (LGVAC(J,IPLS)) CYCLE
                    IF (LGVAC(J,NPLS+1).AND.IFTFLG(KK,2) < 100) CYCLE
                    COU = EIRENE_RATE_COEFF(KK,TEINL(J),0._DP,
     .                    LEXP,0,ERATE)
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
                END IF
                MODCOL(6,2,IRRC)=1
C             ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT
              ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
                IF (NSTORDR >= NRAD) THEN
                  DO J=1,NSBOX
!pb                    IF (LGVAC(J,IPLS)) CYCLE
                    IF (LGVAC(J,NPLS+1).AND.IFTFLG(KK,2) < 100) CYCLE
                    COU = EIRENE_RATE_COEFF(KK,TEINL(J),PLS(J),
     .                   .TRUE.,1,ERATE)
                    TABRC1(IRRC,J)=COU*FACTKK
                    IF (IFTFLG(KK,2) < 100)
     .                TABRC1(IRRC,J)=TABRC1(IRRC,J)*DEIN(J)
                  END DO
 
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                ELSE
C  DON'T STORE DATA, BUT COMPUTE THEM WHEN NEEDED
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                END IF
                MODCOL(6,2,IRRC)=1
              ENDIF
              FACREA(KK,1) = FACTKK
              FACREA(KK,2) = LOG(FACTKK)
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ELECTRON ENERGY LOSS RATE
C
              NSERC5=EIRENE_IDEZ(ISCDEP(IPLS,NRC),5,5)
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
                IF ((KREAD < 1) .OR. (KREAD > NREACI)) GOTO 996
                MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
                LADAS = EIRENE_IS_RTCEW_ADAS(KREAD)
                Z = NCHRGP(IPLS)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >= NRAD) THEN
 
                    DO J = 1, NSBOX
!pb                      IF (LGVAC(J,IPLS)) CYCLE
                      IF (LGVAC(J,NPLS+1).AND.IFTFLG(KK,2) < 100) CYCLE
                      IF (LGVAC(J,NPLS+1).OR.(NCHRGP(IPLS)==0)) THEN
                        BREMS = 0._DP
                      ELSE
                        BREMS = 1.54E-32_DP * TEIN(J)**0.5 * Z**2 *
     .                          eirene_ngffmh(Z**2 * 13.6_DP/TEIN(J)) *
     .                          DEIN(J)*FACTKK/ELCHA
                      END IF
                      EELRC1(IRRC,J)=EIRENE_ENERGY_RATE_COEFF(KREAD,
     .                               TEINL(J),
     .                               0._DP,.TRUE.,0)*DEIN(J)*FACTKK
C  SUBTRACT BREMSTRAHLUNG FROM ADAS PRB RATE
                      IF (LADAS) EELRC1(IRRC,J) = EELRC1(IRRC,J) + BREMS
                    END DO
 
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
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
!pb                      IF (LGVAC(J,IPLS)) CYCLE
                      IF (LGVAC(J,NPLS+1).AND.IFTFLG(KK,2) < 100) CYCLE
                      IF (LGVAC(J,NPLS+1).OR.(NCHRGP(IPLS)==0)) THEN
                        BREMS = 0._DP
                      ELSE
                        BREMS = 1.54E-32_DP * TEIN(J)**0.5 * Z**2 *
     .                          eirene_ngffmh(Z**2 * 13.6_DP/TEIN(J)) *
     .                          DEIN(J)*FACTKK/ELCHA
                      END IF
                      EELRC1(IRRC,J)=EIRENE_ENERGY_RATE_COEFF(KREAD,
     .                               TEINL(J),
     .                               PLS(J),.FALSE.,1)
                      EEMX=MAX(-100._DP,EELRC1(IRRC,J)+FCTKKL+DEINL(J))
                      EELRC1(IRRC,J)=-EXP(EEMX)
C  SUBTRACT BREMSTRAHLUNG FROM ADAS PRB RATE
                      IF (LADAS) EELRC1(IRRC,J) = EELRC1(IRRC,J) + BREMS
                    END DO
 
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                  END IF
                  MODCOL(6,4,IRRC)=1
                ENDIF
 
                FACREA(KREAD,1) = FACTKK
                FACREA(KREAD,2) = LOG(FACTKK)
 
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
          CALL EIRENE_MASBOX
     .  ('BULK ION SPECIES IPLS = '//TEXTS(NSPAMI+IPLS))
          CALL EIRENE_LEER(1)
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
          CALL EIRENE_LEER(1)
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR RECOMBINATION'
      CALL EIRENE_EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED '
      WRITE (iunout,*) 'NREC TOO SMALL, CHECK PARAMETER STATEMENTS'
      CALL EIRENE_EXIT_OWN(1)
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED '
      WRITE (iunout,*) 'MASS CONSERVATION VIOLATED, IPLS,IRRC ',
     .                  IPLS,IRRC
      CALL EIRENE_EXIT_OWN(1)
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTP.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IPLS = ',IPLS
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      EIRENE_'
      CALL EIRENE_EXIT_OWN(1)
995   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED '
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED '
      WRITE (iunout,*)
     .  'WRONG REACTION INDEX SPECIFIED FOR KREAD IN REACTION KK'
      WRITE (iunout,*) 'KK ',KK
      WRITE (IUNOUT,*) 'KREAD ',KREAD
      CALL EIRENE_EXIT_OWN(1)
C
      END
