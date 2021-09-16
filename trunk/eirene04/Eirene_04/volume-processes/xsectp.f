C
C
      SUBROUTINE XSECTP
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR BULK ION SPECIES
C
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
      USE PHOTON

      IMPLICIT NONE
C
      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: DELE, FCTKKL, EEMX, ZX, DSUB, DEIMIN, RMASS2, FACTKK, 
     .            RMASS2_2
      INTEGER :: IIRC, IION3, IPLS3, IATM3, IMOL3, KK, NRC, IATM,
     .           IRRC, J, IDSC, IPLS, NSERC5, KREAD, I, MODC, IATM1,
     .           ITYP, ISPZ, ITYP2, ISPZ2, IPHOT3
      INTEGER, EXTERNAL :: IDEZ
      SAVE
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 70 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
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
                    TABRC1(IRRC,J)=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(J)
                    EELRC1(IRRC,J)=-1.5*TEIN(J)*TABRC1(IRRC,J)
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
              ENDIF
52          CONTINUE
C
            MODCOL(6,2,NSPAMI+IPLS,1)=1
            MODCOL(6,4,NSPAMI+IPLS,1)=1
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
               call PH_XSECTP(ipls,nrc,idsc,nrrci)
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
csw added branch
              ELSEIF (ITYP.EQ.0) THEN
                NPHPRC(IRRC)=ISPZ
                RMASS2=0.
              ENDIF

              ITYP2=IDEZ(ISCD2P(IPLS,NRC),1,3)
              ISPZ2=IDEZ(ISCD2P(IPLS,NRC),3,3)
              IF (ITYP2.EQ.4) THEN
                NPLPRC_2=ISPZ2
                RMASS2_2=RMASSP(ISPZ2)
              ELSE
                RMASS2_2=0._DP
              ENDIF
C  CHECK MASS CONSERVATION
              IF (RMASSP(IPLS).NE.(RMASS2+RMASS2_2)) GOTO 993
C
C  1.) CROSS SECTION(TE)
C           NOT NEEDED
C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3)
C
C  2.A) RATE COEFFICIENT = CONST.
C           TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
              IF (IDEZ(MODCLF(KK),3,5).EQ.1) THEN
                IF (NSTORDR >= NRAD) THEN
                  IF (MOD(IFTFLG(KK,2),100) == 10) THEN
C  RATE:  (1/S)
                    COUN(1,1:NSBOX)=CREAC(1,1,KK)
                  ELSE
C  RATE COEFFICIENT: (CM^3/S)
                    CALL CDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,
     .                         .FALSE.,.TRUE.)
                  END IF
                  IF (IFTFLG(KK,2) < 100) THEN
                    DO J=1,NSBOX
                      TABRC1(IRRC,J)=TABRC1(IRRC,J)+
     +                               COUN(1,J)*DEIN(J)*FACTKK
                    END DO
                  ELSE
                    DO J=1,NSBOX
                      TABRC1(IRRC,J)=TABRC1(IRRC,J)+
     +                               COUN(1,J)*FACTKK
                    END DO
                  ENDIF
!pb
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 1
                ELSE
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 1
                  FACREA(KK) = LOG(FACTKK)
                END IF
                MODCOL(6,2,NSPAMI+IPLS,1)=1
C             ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT
              ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
                IF (NSTORDR >= NRAD) THEN
                  CALL CDEFN(TEINL,PLS,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     ,                       .TRUE.)
                  IF (IFTFLG(KK,2) < 100) THEN
                    DO 93 J=1,NSBOX
                      TABRC1(IRRC,J)=COUN(1,J)*DEIN(J)*FACTKK
93                  CONTINUE
                  ELSE
                    DO J=1,NSBOX
                      TABRC1(IRRC,J)=COUN(1,J)*FACTKK
                    END DO
                  ENDIF
!pb
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                ELSE
                  FACREA(KK) = FACTKK
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                END IF
                MODCOL(6,2,NSPAMI+IPLS,1)=1
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
                MODCOL(6,4,NSPAMI+IPLS,1)=1
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
                MODCOL(6,4,NSPAMI+IPLS,1)=1
              ELSEIF (NSERC5.EQ.3) THEN
C  4.C)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELECP(IPLS,NRC)
                MODC=IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >= NRAD) THEN
                    CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                         .FALSE.,.TRUE.)
                    DO 104 J=1,NSBOX
                      EELRC1(IRRC,J)=-COUN(1,J)*DEIN(J)*FACTKK
104                 CONTINUE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                    FACREA(KREAD) = LOG(FACTKK)
                  END IF
                  MODCOL(6,4,NSPAMI+IPLS,1)=1
C  4.D)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C               ELSEIF (MODC.EQ.2) THEN
C        IRRELEVANT
C                 MODCOL(6,4,NSPAMI+IPLS,1)=2
C  4.E)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    CALL CDEF (TEINL,1,9,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                         .FALSE.,.TRUE.)
                    DO 106 J=1,NSBOX
                      EELRC1(IRRC,J)=COUN(9,J)
106                 CONTINUE
                    DO 108 I=8,1,-1
                      DO 107 J=1,NSBOX
                        EELRC1(IRRC,J)=EELRC1(IRRC,J)*PLS(J)+
     .                                      COUN(I,J)
107                   CONTINUE
108                 CONTINUE
                    FCTKKL=LOG(FACTKK)
                    DO 109 J=1,NSBOX
                      EEMX=MAX(-100._DP,EELRC1(IRRC,J)+FCTKKL+DEINL(J))
                      EELRC1(IRRC,J)=-EXP(EEMX)
109                 CONTINUE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                    FACREA(KREAD) = LOG(FACTKK)
                  END IF
                  MODCOL(6,4,NSPAMI+IPLS,1)=1
                ENDIF
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
            WRITE (6,*) 'NO RECOMBINATION '
          ELSE
            DO 220 IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              WRITE (6,*) 'RECOMBINATION NO. IRRC= ',IRRC
              WRITE (6,*) 'RECOMBINATION INTO SPECIES:'
              IION3=NIOPRC(IRRC)
              IF (IION3.NE.0) WRITE (6,*) 'TEST ION IION= ',
     .                                     TEXTS(NSPAM+IION3)
              IPLS3=NPLPRC(IRRC)
              IF (IPLS3.NE.0) WRITE (6,*) 'BULK ION IPLS= ',
     .                                     TEXTS(NSPAMI+IPLS3)
              IATM3=NATPRC(IRRC)
              IF (IATM3.NE.0) WRITE (6,*) 'ATOM     IATM= ',
     .                                     TEXTS(NSPH+IATM3)
              IMOL3=NMLPRC(IRRC)
              IF (IMOL3.NE.0) WRITE (6,*) 'MOLECULE IMOL= ',
     .                                     TEXTS(NSPA+IMOL3)
              IPHOT3=NPHPRC(IRRC)
              IF (IPHOT3.NE.0) WRITE (6,*) 'PHOTON  IPHOT= ',
     .                                     TEXTS(IPHOT3)
C             WRITE (6,*) 'ELECTRONS: PELPRC,EELRC1'
C             IF (NSTORDR >= NRAD) THEN
C               WRITE (6,*) 'EL      ',1.,EELRC1(IRRC,1)
C             ELSE
C               WRITE (6,*) 'EL      ',1.,FEELRC1(IRRC,1)
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
      WRITE (6,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR RECOMBINATION'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (6,*) 'NREC TOO SMALL, CHECK PARAMETER STATEMENTS'
      CALL EXIT_OWN(1)
993   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (6,*) 'MASS CONSERVATION VIOLATED, IPLS,IRRC ',IPLS,IRRC
      CALL EXIT_OWN(1)
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTP.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IPLS = ',IPLS
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
C
      END
