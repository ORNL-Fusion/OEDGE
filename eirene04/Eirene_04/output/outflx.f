C
C
C 22.10.03; ein falsches write statement rausgenommen (bei reemitted
C           from incident test ions, totals.
C 27.03.04; iliin=-3 option (only net fluxes on transp. surfaces) re-enforced
C           simultaneous changes in escape.f.
C           not active for bulk particle fluxes updated in subr. locate  
C           not yet in repository
C 05.05.04:  printout of spectrum: text improved
C
C
      SUBROUTINE OUTFLX(A,ISTRA)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CGRID
      USE CSPEZ
      USE CTRCEI
      USE CTEXT
      USE CSDVI
      USE CLGIN
      USE COUTAU
      USE CTRIG

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      CHARACTER(32), INTENT(IN) :: A
      REAL(DP) :: SUMA1(0:NATM,0:NSTRA), VARA1(0:NATM,0:NSTRA),
     .          SUMM1(0:NMOL,0:NSTRA), VARM1(0:NMOL,0:NSTRA),
     .          SUMI1(0:NION,0:NSTRA), VARI1(0:NION,0:NSTRA),
     .          SUMPH1(0:NPHOT,0:NSTRA), VARPH1(0:NPHOT,0:NSTRA),
     .          SUMP1(0:NPLS,0:NSTRA), VARP1(0:NPLS,0:NSTRA),
     .          SUMA2(0:NATM,0:NSTRA), VARA2(0:NATM,0:NSTRA),
     .          SUMM2(0:NMOL,0:NSTRA), VARM2(0:NMOL,0:NSTRA),
     .          SUMI2(0:NION,0:NSTRA), VARI2(0:NION,0:NSTRA),
     .          SUMPH2(0:NPHOT,0:NSTRA), VARPH2(0:NPHOT,0:NSTRA),
     .          SUMP2(0:NPLS,0:NSTRA), VARP2(0:NPLS,0:NSTRA),
     .          SUMS(0:NADS,0:NSTRA),  VARS(0:NADS,0:NSTRA),
     .          SUML(0:NALS,0:NSTRA),  VARL(0:NALS,0:NSTRA)
      REAL(DP) :: HELP(NRAD), HELPP(NLMPGS)
      REAL(DP) :: PGAINE, SUM1, DUMMY, SUMMT, SUMMEM, SUMMTI, SUMMEI,
     .          SUMMTA, SUMMTM, SUMMA, SUMMS, SUMMM, SUMMI, SUMMP,
     .          SUMA, SUMM, SUMI, SUMME, SUMMTP, SUMMEP, SUMP, TTTT,
     .          SUMMEA, SUMMEPH, SUMMTPH, SUMMPH, SUMPH, EN
      INTEGER :: NR, NP, NT, MSURFG, J, NCELL, NTOTAL, N1, N2, N3, I,
     .           ITRII, IPLGN, ISPR, ITEXT, ISTS, NFTI, NFTE, NTCO,
     .           NF, K, ITALS, I0, IADS, IALS, IION, IPLS, IATM, N,
     .           IMOL, IPHOT, ISPC, IOUT, ISF, IE, IT
      INTEGER :: IADTYP(0:4)
      LOGICAL :: LGVRA1(0:NATM,0:NSTRA),
     .           LGVRM1(0:NMOL,0:NSTRA),
     .           LGVRI1(0:NION,0:NSTRA),
     .           LGVRPH1(0:NPHOT,0:NSTRA),
     .           LGVRP1(0:NPLS,0:NSTRA)
      LOGICAL :: LGVRA2(0:NATM,0:NSTRA),
     .           LGVRM2(0:NMOL,0:NSTRA),
     .           LGVRI2(0:NION,0:NSTRA),
     .           LGVRPH2(0:NPHOT,0:NSTRA),
     .           LGVRP2(0:NPLS,0:NSTRA)
      LOGICAL :: LGVARS(0:NADS,0:NSTRA),
     .           LGVARL(0:NALS,0:NSTRA)
      LOGICAL :: LOGADS(0:NADS,0:NSTRA),
     .           LOGALS(0:NALS,0:NSTRA)
      LOGICAL :: PRINTED(NLIMPS)
      CHARACTER(8) :: TEXTA(NADS), TEXTL(NALS)
      CHARACTER(24) :: TXT24
      CHARACTER(72) :: TXT72
      CHARACTER(10) :: TEXTYP(0:4)

      PRINTED = .FALSE.
      CALL LEER(1)
      WRITE (6,9999) A
C
C  SURFACE LOOP
C
      DO 10000 ISPR=1,NSURPR
C
C
        ITEXT=0
        I=NPRSRF(ISPR)
        CALL LEER(2)
        IF (IGJUM0(I).NE.0) THEN
C  CHECK, IF SURFACE "I" IS STILL THERE AS ONE SIDE OF A TRIANGLE
          IF (LEVGEO.EQ.4) THEN
            DO ITRII=1,NTRII
            DO IPLGN=1,3
              ISTS=ABS(INMTI(IPLGN,ITRII))
              IF (ISTS.EQ.I) GOTO 5
            ENDDO
            ENDDO
          ENDIF
          WRITE (6,*) ' SURFACE NO. ',I,' : OUT '
          PRINTED(I) = .TRUE.
          GOTO 10000
        ENDIF
5       CONTINUE
C  PRINT SURFACE AREA (NOT FOR "TIME SURFACE")
        CALL MASBOX(TXTSFL(I))
        IF (NTIME.GE.1.AND.I.EQ.NLIM+NSTSI) GOTO 1
        IF (SAREA(I).NE.666.) THEN
          WRITE (6,'(A22,1P,1E12.4)') ' SURFACE AREA (CM**2) ',SAREA(I)
        ELSE
          WRITE (6,*) 'SURFACE AREA (CM**2) ','?'
        ENDIF
1       CONTINUE
C
C
        if (i.gt.nlim.and.levgeo.le.4.and.nlmpgs.ne.nlimps) then
C
C  SPATIAL RESOLUTION ON NON DEFAULT STANDARD SURFACE?
          HELP=0.D0
          ISTS=I-NLIM
C
          ITALS=NPRTLS(ISPR)
          IF (ITALS.LE.0.OR.ITALS.GT.NTALS) GOTO 11
          I0=0
          IF (NFRTWI(ITALS).GT.1) I0=1
          NFTI=1
          NFTE=NFSTWI(ITALS)
          IF (NSPEZS(ISPR,1).GT.0) THEN
            NFTI=NSPEZS(ISPR,1)
            NFTE=MAX(NFTI,NSPEZS(ISPR,2))
          ENDIF
          NF=NFRSTW(ITALS)
          DO 10 K=NFTI,NFTE
            IF (K.GT.NFSTWI(ITALS)) THEN
              CALL LEER(1)
              WRITE (6,*) 'SPECIES INDEX OUT OF RANGE IN OUTFLX'
              WRITE (6,*) 'ISTS,ITALS, K, ',ISTS,ITALS,K
              CALL LEER(1)
              GOTO 10
            ENDIF
C
            DO J=1,NLMPGS
              HELPP(J)=ESTIMS(NADDW(ITALS)+K,J)
            END DO
C
            IF (LEVGEO.LE.3) THEN
              IF (INUMP(ISTS,2).NE.0) then
C  POLOIDAL SURFACE
                sum1=0
                np=1
                do nr=1,nr1st
                do nt=1,nt3rd
                  MSURFG=NR+(NT-1)*NR1P2
                  MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                  NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                  HELP(ncell)=HELPP(MSURFG)
                  sum1=sum1+HELPP(msurfg)
                ENDDO
                ENDDO
                N1=NR1ST             
                N2=1
                N3=NT3RD
              ELSEIF (INUMP(ISTS,1).NE.0) then
C  RADIAL SURFACE
                sum1=0
                NR=1
                do np=1,np2nd
                do nt=1,nt3rd
                  MSURFG=NP+(NT-1)*NP2T3
                  MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                  NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                  HELP(ncell)=HELPP(MSURFG)
                  sum1=sum1+HELPP(msurfg)
                ENDDO
                ENDDO
                N1=1
                N2=NP2ND
                N3=NT3RD
              ELSEIF (INUMP(ISTS,3).NE.0) then
C  TOROIDAL SURFACE
                sum1=0
                nt=1
                do nr=1,nr1st
                do np=1,np2nd
                  MSURFG=Nr+(Np-1)*Nr1p2
                  MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                  NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                  HELP(ncell)=HELPP(MSURFG)
                  sum1=sum1+HELPP(msurfg)
                ENDDO
                ENDDO
                N1=NR1ST
                N2=NP2ND
                N3=1
              ENDIF
            ELSE IF (LEVGEO.EQ.4) THEN
              sum1=0.D0
              ntco=0
              DO NP=1,3
                DO NR=1,NTRII
                  IF (INMTI(NP,NR) == NLIM+ISTS) THEN
                    MSURFG=NLIM+NSTS+INSPAT(NP,NR)
                    NTCO=NTCO+1
                    HELP(NTCO)=HELPP(MSURFG)
                    SUM1=SUM1+HELPP(MSURFG)
                  END IF
                END DO
              END DO
              N1=NTCO+1
              N2=1
              N3=1
            END IF
            write (6,*) 'test ',sum1
            NTOTAL=N1*N2*N3
            CALL INTVOL (HELP,1,1,NTOTAL,DUMMY,N1,N2,N3,1)
            IF (ABS(DUMMY) > EPS60) THEN
              CALL PRTTLS(TXTTLW(K,ITALS),TXTSPW(K,ITALS),
     .                  TXTUNW(K,ITALS),
     .                  HELP,N1,N2,N3,1,NTOTAL,NFLAGS(ISPR),
     .                  NTLSFL(ISPR),
     .                  IRPTA(ISTS,1),IRPTE(ISTS,1),IRPTA(ISTS,2),
     .                  IRPTE(ISTS,2),1,1)
            ELSE
              CALL PRTTLS(TXTTLW(K,ITALS),TXTSPW(K,ITALS),
     .                  TXTUNW(K,ITALS),
     .                  HELP,N1,N2,N3,1,NTOTAL,-1,NTLSFL(ISPR),
     .                  IRPTA(ISTS,1),IRPTE(ISTS,1),IRPTA(ISTS,2),
     .                  IRPTE(ISTS,2),1,1)
              CALL MASAGE
     .            ('IDENTICAL ZERO, NOT PRINTED                  ')
              CALL LEER(2)
            END IF
10        CONTINUE
11        CONTINUE
        ENDIF


C  SPECTRA
        IF (NTLSFL(ISPR) > 0) THEN
          ISF=I
          IF (I < 0) ISF=ABS(I)+NLIM
          IOUT = NTLSFL(ISPR)

          TEXTYP(0) = 'PHOTONS   '
          TEXTYP(1) = 'ATOMS     '
          TEXTYP(2) = 'MOLECULES '
          TEXTYP(3) = 'TEST IONS '
          TEXTYP(4) = 'BULK IONS '
          IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

          DO ISPC=1,NADSPC
            IF (ESTIML(ISPC)%PSPC%ISPCSRF == ISF) THEN
              WRITE (IOUT,*)
              WRITE (IOUT,*)
              WRITE (IOUT,*) ' SPECTRUM CALCULATED FOR SURFACE ',I
              IT = ESTIML(ISPC)%PSPC%ISPCTYP
              IF (IT == 1) THEN
                WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                            'PARTICLE FLUX IN AMP'
              ELSE
                WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                            'ENERGY FLUX IN WATT '
              END IF
              WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
     .              TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
              IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
                WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
              ELSE
                WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .                 TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .                       ESTIML(ISPC)%PSPC%IPRSP)
              END IF
              WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ', 
     .               ESTIML(ISPC)%PSPC%SPCMIN
              WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ', 
     .               ESTIML(ISPC)%PSPC%SPCMAX
              WRITE (IOUT,'(A16,4x,I6)') ' NUMBER OF BINS ', 
     .               ESTIML(ISPC)%PSPC%NSPC
              WRITE (IOUT,*)
              IF (ESTIML(ISPC)%PSPC%SPCINT > EPS60) THEN
                IF (NSIGI_SPC == 0) THEN
                  DO IE=1, ESTIML(ISPC)%PSPC%NSPC
                    EN = ESTIML(ISPC)%PSPC%SPCMIN + 
     .                   (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
                    WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
     .                 ESTIML(ISPC)%PSPC%SPC(IE)
                  END DO
                ELSE
                  DO IE=1, ESTIML(ISPC)%PSPC%NSPC
                    EN = ESTIML(ISPC)%PSPC%SPCMIN + 
     .                   (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
                    WRITE (IOUT,'(I6,3ES12.4)') IE,EN,
     .                   ESTIML(ISPC)%PSPC%SPC(IE),
     .                   ESTIML(ISPC)%PSPC%SDV(IE)
                  END DO
                END IF
              ELSE
                WRITE (IOUT,*) ' SPECTRUM IDENTICAL 0 '
              END IF
              WRITE (IOUT,*)
              WRITE (IOUT,*) ' INTEGRAL OF SPECTRUM ',
     .               ESTIML(ISPC)%PSPC%SPCINT
              IF (NSIGI_SPC > 0)
     .          WRITE (IOUT,*) ' STANDARD DEVIATION  ',
     .               ESTIML(ISPC)%PSPC%SGMS
            END IF
          END DO
        END IF

        IF (PRINTED(I)) CYCLE
C
C  *****************************************************
C   INCIDENT FLUXES, POSITIVE PARTIAL FLUXES, NET FLUXES
C  *****************************************************
C
C
C   SURFACE AVERAGED TALLY NO.1 AND NO.26
C
        SUMMT=0.
        SUMME=0.
        SUMMTP=0.
        SUMMEP=0.
        SUMA=0.
        SUMM=0.
        SUMI=0.
        SUMP=0.
        SUMPH=0.
        DO 20 IATM=1,NATMI
          LGVRA1(IATM,ISTRA)=.FALSE.
          LGVRA2(IATM,ISTRA)=.FALSE.
          SUMA1(IATM,ISTRA)=POTAT(IATM,I)
          SUMMT=SUMMT+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
          SUMA=SUMA+SUMA1(IATM,ISTRA)
          SUMA2(IATM,ISTRA)=EOTAT(IATM,I)
          SUMME=SUMME+SUMA2(IATM,ISTRA)
20      CONTINUE
C
        DO 21 N=1,NSIGSI
          IF (IIHW(N).EQ.1) THEN
            DO 22 IATM=1,NATMI
              IF (IGHW(N).NE.IATM) GOTO 22
              VARA1(IATM,ISTRA)=SIGMAW(N,I)
              LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
22          CONTINUE
          ELSEIF (IIHW(N).EQ.26) THEN
            DO 522 IATM=1,NATMI
              IF (IGHW(N).NE.IATM) GOTO 522
              VARA2(IATM,ISTRA)=SIGMAW(N,I)
              LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
522         CONTINUE
          ENDIF
21      CONTINUE
C
C   SURFACE AVERAGED TALLY NO.7 AND NO.32
C
        DO 30 IMOL=1,NMOLI
          LGVRM1(IMOL,ISTRA)=.FALSE.
          LGVRM2(IMOL,ISTRA)=.FALSE.
          SUMM1(IMOL,ISTRA)=POTML(IMOL,I)
          SUMMT=SUMMT+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
          SUMM=SUMM+SUMM1(IMOL,ISTRA)
          SUMM2(IMOL,ISTRA)=EOTML(IMOL,I)
          SUMME=SUMME+SUMM2(IMOL,ISTRA)
30      CONTINUE
C
        DO 23 N=1,NSIGSI
          IF (IIHW(N).EQ.7) THEN
            DO 24 IMOL=1,NMOLI
              IF (IGHW(N).NE.IMOL) GOTO 24
              VARM1(IMOL,ISTRA)=SIGMAW(N,I)
              LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
24          CONTINUE
          ELSEIF (IIHW(N).EQ.32) THEN
            DO 524 IMOL=1,NMOLI
              IF (IGHW(N).NE.IMOL) GOTO 524
              VARM2(IMOL,ISTRA)=SIGMAW(N,I)
              LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
524         CONTINUE
          ENDIF
23      CONTINUE
C
C   SURFACE AVERAGED TALLY NO.13 AND NO.38
C
      DO 40 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=POTIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        SUMI2(IION,ISTRA)=EOTIO(IION,I)
        SUMME=SUMME+SUMI2(IION,ISTRA)
40    CONTINUE
C
      DO 25 N=1,NSIGSI
        IF (IIHW(N).EQ.13) THEN
          DO 26 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 26
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
26        CONTINUE
        ELSEIF (IIHW(N).EQ.38) THEN
          DO 526 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 526
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
526       CONTINUE
        ENDIF
25    CONTINUE
C
C   SURFACE AVERAGED TALLY NO 19 and NO 44
C
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=POTPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        SUMPH2(IPHOT,ISTRA)=EOTPHT(IPHOT,I)
        SUMME=SUMME+SUMPH2(IPHOT,ISTRA)
      ENDDO

      DO N=1,NSIGSI
        IF (IIHW(N).EQ.19) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.44) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE AVERAGED TALLY NO.25 AND NO.50
C
      DO 50 IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        LGVRP2(IPLS,ISTRA)=.FALSE.
        SUMP1(IPLS,ISTRA)=POTPL(IPLS,I)
        SUMMTP=SUMMTP+SUMP1(IPLS,ISTRA)*NPRT(NSPAMI+IPLS)
        SUMP=SUMP+SUMP1(IPLS,ISTRA)
        SUMP2(IPLS,ISTRA)=EOTPL(IPLS,I)
        SUMMEP=SUMMEP+SUMP2(IPLS,ISTRA)
50    CONTINUE
C
      DO 27 N=1,NSIGSI
        IF (IIHW(N).EQ.25) THEN
          DO 28 IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) GOTO 28
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
28        CONTINUE
        ELSEIF (IIHW(N).EQ.50) THEN
          DO 528 IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) GOTO 528
            VARP2(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP2(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
528       CONTINUE
        ENDIF
27    CONTINUE
C
C
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMP)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO FLUXES INCIDENT ON THIS SURFACE             ')
        CALL LEER(1)
      ELSE
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
          WRITE (6,*) 'FLUX INCIDENT ON SURFACE:'
C  SURFACE AVERAGED TALLY NO. 1
          IF (SUMA.NE.0.D0) THEN
            WRITE (6,*) 'INCIDENT: ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 26
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 7
          IF (SUMM.NE.0.D0) THEN
            WRITE (6,*) 'INCIDENT: MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 32
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 13
          IF (SUMI.NE.0.D0) THEN
            WRITE (6,*) 'INCIDENT: TEST IONS'
            CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 38
            CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0.D0) THEN
            WRITE (6,*) 'INCIDENT: PHOTONS'
            CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 44
            CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 25
          IF (SUMP.NE.0.D0) THEN
            WRITE (6,*) 'INCIDENT: BULK IONS (RECYCLING SOURCE)'
            CALL MASYR1('P-FLUX:  ',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
C  SURFACE AVERAGED TALLY NO. 50
            CALL MASYR1('E-FLUX:  ',SUMP2,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP2,LGVRP2,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
          ENDIF
          CALL LEER (1)
          WRITE (6,*) 'TOTAL INCIDENT "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (6,*) '(EXCLUDING BULK IONS (RECYCLING SOURCE) '
          CALL MASR1 ('TOT.PFLX',SUMMT)
          CALL MASR1 ('TOT.EFLX',SUMME)
C
        ELSEIF (ILIIN(I).LT.0.AND.ILIIN(I).NE.-3) THEN
          IF (SUMA.NE.0.D0.OR.SUMM.NE.0.D0.OR.SUMI.NE.0.D0
     .                    .OR.SUMPH.NE.0.D0)
     .    WRITE (6,*) 'PARTIAL PARTICLE AND ENERGY CURRENTS, POSITIVE '
C  SURFACE AVERAGED TALLY NO. 1
          IF (SUMA.NE.0.D0) THEN
            WRITE (6,*) 'ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 26
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 7
          IF (SUMM.NE.0.D0) THEN
            WRITE (6,*) 'MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 32
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 13
          IF (SUMI.NE.0.D0) THEN
            WRITE (6,*) 'TEST IONS'
            CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 38
            CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0.D0) THEN
            WRITE (6,*) 'PHOTONS'
            CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 44
            CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 25
          IF (SUMP.NE.0.D0) THEN
            WRITE (6,*) 'FLUX INCIDENT ON SURFACE (RECYCLING SOURCE):'
            WRITE (6,*) 'BULK IONS'
            CALL MASYR1('P-FLUX:  ',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
C  SURFACE AVERAGED TALLY NO. 50
            CALL MASYR1('E-FLUX:  ',SUMP2,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP2,LGVRP2,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
          ENDIF
          CALL LEER (1)
          WRITE (6,*) 'TOTAL POSITIVE "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (6,*) '(EXCLUDING BULK IONS (RECYCLING SOURCE) '
          CALL MASR1 ('POS.PFLX',SUMMT)
          CALL MASR1 ('POS.EFLX',SUMME)
C
        ELSEIF (ILIIN(I).EQ.-3) THEN
          IF (SUMA.NE.0.D0.OR.SUMM.NE.0.D0.OR.SUMI.NE.0.D0
     .                    .OR.SUMPH.NE.0.D0)
     .    WRITE (6,*) 'NET PARTICLE AND ENERGY CURRENTS '
C  SURFACE AVERAGED TALLY NO. 1
          IF (SUMA.NE.0.D0) THEN
            WRITE (6,*) 'ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 26
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 7
          IF (SUMM.NE.0.D0) THEN
            WRITE (6,*) 'MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 32
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 13
          IF (SUMI.NE.0.D0) THEN
            WRITE (6,*) 'TEST IONS'
            CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 38
            CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0.D0) THEN
            WRITE (6,*) 'PHOTONS'
            CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 44
            CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 25
          IF (SUMP.NE.0.D0) THEN
            WRITE (6,*) 'FLUX INCIDENT ON SURFACE (RECYCLING SOURCE):'
            WRITE (6,*) 'BULK IONS'
            CALL MASYR1('P-FLUX:  ',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
C  SURFACE AVERAGED TALLY NO. 50
            CALL MASYR1('E-FLUX:  ',SUMP2,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP2,LGVRP2,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
          ENDIF
          CALL LEER (1)
          WRITE (6,*) 'TOTAL NET "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (6,*) '(EXCLUDING BULK IONS (RECYCLING SOURCE) '
          CALL MASR1 ('NET PFLX',SUMMT)
          CALL MASR1 ('NET EFLX',SUMME)


        ENDIF
C
C  INDEPENDENT OF VALUE AND SIGN OF ILIIN:
C
        IF (SUMMTP.GT.0.D0) THEN
          CALL LEER (1)
          WRITE (6,*) 'TOTAL INCIDENT RECYCLING SOURCE "ATOMIC" FLUXES'
          WRITE (6,*) 'BULK IONS'
          CALL MASR1 ('SRC.PFLX',SUMMTP)
          CALL MASR1 ('SRC.EFLX',SUMMEP)
        ENDIF
      ENDIF
C
C  ******************************************
C   REEMITTED FLUXES, NEGATIVE PARTIAL FLUXES
C  ******************************************
C
C   FIRST: FROM INCIDENT ATOMS
C
      SUMMTA=0.
      SUMMEA=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.2 AND NO.27
C
      DO 102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        SUMA1(IATM,ISTRA)=PRFAAT(IATM,I)
        SUMMTA=SUMMTA+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        SUMA2(IATM,ISTRA)=ERFAAT(IATM,I)
        SUMMEA=SUMMEA+SUMA2(IATM,ISTRA)
102   CONTINUE
C
      DO 121 N=1,NSIGSI
        IF (IIHW(N).EQ.2) THEN
          DO 122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
122       CONTINUE
        ELSEIF (IIHW(N).EQ.27) THEN
          DO 622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
622       CONTINUE
        ENDIF
121   CONTINUE
C
C
C   SURFACE AVERAGED TALLY NO.8 AND NO.33
C
      DO 103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        SUMM1(IMOL,ISTRA)=PRFAML(IMOL,I)
        SUMMTA=SUMMTA+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        SUMM2(IMOL,ISTRA)=ERFAML(IMOL,I)
        SUMMEA=SUMMEA+SUMM2(IMOL,ISTRA)
103   CONTINUE
C
      DO 123 N=1,NSIGSI
        IF (IIHW(N).EQ.8) THEN
          DO 124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
124       CONTINUE
        ELSEIF (IIHW(N).EQ.33) THEN
          DO 624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
624       CONTINUE
        ENDIF
123   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.14 AND NO.39
C
      DO 104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=PRFAIO(IION,I)
        SUMMTA=SUMMTA+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        SUMI2(IION,ISTRA)=ERFAIO(IION,I)
        SUMMEA=SUMMEA+SUMI2(IION,ISTRA)
104   CONTINUE
C
      DO 125 N=1,NSIGSI
        IF (IIHW(N).EQ.14) THEN
          DO 126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
126       CONTINUE
        ELSEIF (IIHW(N).EQ.39) THEN
          DO 626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
626       CONTINUE
        ENDIF
125   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.20 AND NO.45
C
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=PRFAPHT(IPHOT,I)
        SUMMTA=SUMMTA+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        SUMPH2(IPHOT,ISTRA)=ERFAPHT(IPHOT,I)
        SUMMEA=SUMMEA+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.20) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.45) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'NO FLUXES REEMITTED FROM INCIDENT ATOMS '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (6,*) 'FLUX REEMITTED FROM INCIDENT ATOMS:'
C  SURFACE AVERAGED TALLY NO. 2
        IF (SUMA.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 27
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 8
        IF (SUMM.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 33
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 14
        IF (SUMI.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 39
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 20
        IF (SUMPH.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 45
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTA)
        CALL MASR1 ('TOT.EFLX',SUMMEA)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 2
          IF (SUMA.NE.0.D0) THEN
            WRITE (6,*) 'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
            ITEXT=1
            WRITE (6,*) 'ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 27
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL LEER (1)
            CALL MASR1 ('NEG.PFLX',SUMMTA)
            CALL MASR1 ('NEG.EFLX',SUMMEA)
          ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT MOLECULES
      SUMMTM=0.
      SUMMEM=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.3 AND NO.28
C
      DO 1102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        SUMA1(IATM,ISTRA)=PRFMAT(IATM,I)
        SUMMTM=SUMMTM+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        SUMA2(IATM,ISTRA)=ERFMAT(IATM,I)
        SUMMEM=SUMMEM+SUMA2(IATM,ISTRA)
1102  CONTINUE
C
      DO 1121 N=1,NSIGSI
        IF (IIHW(N).EQ.3) THEN
          DO 1122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 1122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
1122      CONTINUE
        ELSEIF (IIHW(N).EQ.28) THEN
          DO 1622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 1622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
1622      CONTINUE
        ENDIF
1121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.9 AND NO.34
C
      DO 1103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        SUMM1(IMOL,ISTRA)=PRFMML(IMOL,I)
        SUMMTM=SUMMTM+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        SUMM2(IMOL,ISTRA)=ERFMML(IMOL,I)
        SUMMEM=SUMMEM+SUMM2(IMOL,ISTRA)
1103  CONTINUE
C
      DO 1123 N=1,NSIGSI
        IF (IIHW(N).EQ.9) THEN
          DO 1124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 1124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
1124      CONTINUE
        ELSEIF (IIHW(N).EQ.34) THEN
          DO 1624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 1624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
1624      CONTINUE
        ENDIF
1123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.15 AND NO.40
C
      DO 1104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=PRFMIO(IION,I)
        SUMMTM=SUMMTM+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        SUMI2(IION,ISTRA)=ERFMIO(IION,I)
        SUMMEM=SUMMEM+SUMI2(IION,ISTRA)
1104  CONTINUE
C
      DO 1125 N=1,NSIGSI
        IF (IIHW(N).EQ.15) THEN
          DO 1126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 1126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
1126      CONTINUE
        ELSEIF (IIHW(N).EQ.40) THEN
          DO 1626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 1626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
1626      CONTINUE
        ENDIF
1125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.21 AND NO.46
C
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=PRFMPHT(IPHOT,I)
        SUMMTM=SUMMTM+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        SUMPH2(IPHOT,ISTRA)=ERFMPHT(IPHOT,I)
        SUMMEM=SUMMEM+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.21) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.46) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'NO FLUXES REEMITTED FROM INCIDENT MOLECULES '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (6,*) 'FLUX REEMITTED FROM INCIDENT MOLECULES:'
C  SURFACE AVERAGED TALLY NO. 3
        IF (SUMA.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 28
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 9
        IF (SUMM.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 34
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 15
        IF (SUMI.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 40
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 21
        IF (SUMPH.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 46
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTM)
        CALL MASR1 ('TOT.EFLX',SUMMEM)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 9
          IF (SUMM.NE.0.D0) THEN
            IF (ITEXT.EQ.0)
     .      WRITE (6,*) 'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
            ITEXT=1
            WRITE (6,*) 'MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 34
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL LEER (1)
            CALL MASR1 ('NEG.PFLX',SUMMTM)
            CALL MASR1 ('NEG.EFLX',SUMMEM)
          ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT TEST-IONS
      SUMMTI=0.
      SUMMEI=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.4 AND NO.29
C
      DO 2102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        SUMA1(IATM,ISTRA)=PRFIAT(IATM,I)
        SUMMTI=SUMMTI+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        SUMA2(IATM,ISTRA)=ERFIAT(IATM,I)
        SUMMEI=SUMMEI+SUMA2(IATM,ISTRA)
2102  CONTINUE
C
      DO 2121 N=1,NSIGSI
        IF (IIHW(N).EQ.4) THEN
          DO 2122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 2122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
2122      CONTINUE
        ELSEIF (IIHW(N).EQ.29) THEN
          DO 2622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 2622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
2622      CONTINUE
        ENDIF
2121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.10 AND NO.35
C
      DO 2103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        SUMM1(IMOL,ISTRA)=PRFIML(IMOL,I)
        SUMMTI=SUMMTI+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        SUMM2(IMOL,ISTRA)=ERFIML(IMOL,I)
        SUMMEI=SUMMEI+SUMM2(IMOL,ISTRA)
2103  CONTINUE
C
      DO 2123 N=1,NSIGSI
        IF (IIHW(N).EQ.10) THEN
          DO 2124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 2124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
2124      CONTINUE
        ELSEIF (IIHW(N).EQ.35) THEN
          DO 2624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 2624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
2624      CONTINUE
        ENDIF
2123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.16 AND NO.41
C
      DO 2104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=PRFIIO(IION,I)
        SUMMTI=SUMMTI+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        SUMI2(IION,ISTRA)=ERFIIO(IION,I)
        SUMMEI=SUMMEI+SUMI2(IION,ISTRA)
2104  CONTINUE
C
      DO 2125 N=1,NSIGSI
        IF (IIHW(N).EQ.16) THEN
          DO 2126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 2126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
2126      CONTINUE
        ELSEIF (IIHW(N).EQ.41) THEN
          DO 2626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 2626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
2626      CONTINUE
        ENDIF
2125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO 22 AND NO.47
C
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=PRFIPHT(IPHOT,I)
        SUMMTI=SUMMTI+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        SUMPH2(IPHOT,ISTRA)=ERFIPHT(IPHOT,I)
        SUMMEI=SUMMEI+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.22) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.47) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'NO FLUXES REEMITTED FROM INCIDENT TEST IONS '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (6,*) 'FLUX REEMITTED FROM INCIDENT TEST IONS:'
C  SURFACE AVERAGED TALLY NO. 4
        IF (SUMA.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 29
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 10
        IF (SUMM.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 35
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 16
        IF (SUMI.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 41
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 22
        IF (SUMPH.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 47
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTI)
        CALL MASR1 ('TOT.EFLX',SUMMEI)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 16
        IF (SUMI.NE.0.D0) THEN
        IF (ITEXT.EQ.0)
     .  WRITE (6,*) 'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE '
        ITEXT=1
        WRITE (6,*) 'TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 41
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL LEER (1)
        CALL MASR1 ('NEG.PFLX',SUMMTI)
        CALL MASR1 ('NEG.EFLX',SUMMEI)
        ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT PHOTONS
      SUMMTPH=0.
      SUMMEPH=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.5 AND NO.30
C
      DO 4102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        SUMA1(IATM,ISTRA)=PRFPHAT(IATM,I)
        SUMMTPH=SUMMTPH+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        SUMA2(IATM,ISTRA)=ERFPHAT(IATM,I)
        SUMMEPH=SUMMEPH+SUMA2(IATM,ISTRA)
4102  CONTINUE
C
      DO 4121 N=1,NSIGSI
        IF (IIHW(N).EQ.5) THEN
          DO 4122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 4122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
4122      CONTINUE
        ELSEIF (IIHW(N).EQ.30) THEN
          DO 4622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 4622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
4622      CONTINUE
        ENDIF
4121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.11 AND NO.36
C
      DO 4103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        SUMM1(IMOL,ISTRA)=PRFPHML(IMOL,I)
        SUMMTPH=SUMMTPH+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        SUMM2(IMOL,ISTRA)=ERFPHML(IMOL,I)
        SUMMEPH=SUMMEPH+SUMM2(IMOL,ISTRA)
4103  CONTINUE
C
      DO 4123 N=1,NSIGSI
        IF (IIHW(N).EQ.11) THEN
          DO 4124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 4124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
4124      CONTINUE
        ELSEIF (IIHW(N).EQ.36) THEN
          DO 4624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 4624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
4624      CONTINUE
        ENDIF
4123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.17 AND NO.42
C
      DO 4104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=PRFPHIO(IION,I)
        SUMMTPH=SUMMTPH+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        SUMI2(IION,ISTRA)=ERFPHIO(IION,I)
        SUMMEPH=SUMMEPH+SUMI2(IION,ISTRA)
4104   CONTINUE
C
      DO 4125 N=1,NSIGSI
        IF (IIHW(N).EQ.17) THEN
          DO 4126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 4126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
4126      CONTINUE
        ELSEIF (IIHW(N).EQ.42) THEN
          DO 4626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 4626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
4626      CONTINUE
        ENDIF
4125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO 23 AND NO.48
C
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=PRFPHPHT(IPHOT,I)
        SUMMTI=SUMMTI+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        SUMPH2(IPHOT,ISTRA)=ERFPHPHT(IPHOT,I)
        SUMMEI=SUMMEI+SUMPH2(IPHOT,ISTRA)
      ENDDO

      DO N=1,NSIGSI
        IF (IIHW(N).EQ.23) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.48) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'NO FLUXES REEMITTED FROM INCIDENT PHOTONS   '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (6,*) 'FLUX REEMITTED FROM INCIDENT PHOTONS:'
C  SURFACE AVERAGED TALLY NO. 5
        IF (SUMA.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 30
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 11
        IF (SUMM.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 36
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 17
        IF (SUMI.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 42
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 23
        IF (SUMPH.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 48
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTPH)
        CALL MASR1 ('TOT.EFLX',SUMMEPH)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 23
        IF (SUMPH.NE.0.D0) THEN
        IF (ITEXT.EQ.0)
     .  WRITE (6,*) 'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE '
        ITEXT=1
        WRITE (6,*) 'PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 48
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL LEER (1)
        CALL MASR1 ('NEG.PFLX',SUMMTPH)
        CALL MASR1 ('NEG.EFLX',SUMMEPH)
        ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT BULK-IONS
      SUMMTP=0.
      SUMMEP=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.6 AND NO.31
C
      DO 3102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        SUMA1(IATM,ISTRA)=PRFPAT(IATM,I)
        SUMMTP=SUMMTP+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        SUMA2(IATM,ISTRA)=ERFPAT(IATM,I)
        SUMMEP=SUMMEP+SUMA2(IATM,ISTRA)
3102  CONTINUE
C
      DO 3121 N=1,NSIGSI
        IF (IIHW(N).EQ.6) THEN
          DO 3122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 3122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
3122      CONTINUE
        ELSEIF (IIHW(N).EQ.31) THEN
          DO 3622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 3622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
3622      CONTINUE
        ENDIF
3121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.12 AND NO.37
C
      DO 3103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        SUMM1(IMOL,ISTRA)=PRFPML(IMOL,I)
        SUMMTP=SUMMTP+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        SUMM2(IMOL,ISTRA)=ERFPML(IMOL,I)
        SUMMEP=SUMMEP+SUMM2(IMOL,ISTRA)
3103  CONTINUE
C
      DO 3123 N=1,NSIGSI
        IF (IIHW(N).EQ.12) THEN
          DO 3124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 3124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
3124      CONTINUE
        ELSEIF (IIHW(N).EQ.37) THEN
          DO 3624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 3624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
3624      CONTINUE
        ENDIF
3123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.18 AND NO.43
C
      DO 3104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=PRFPIO(IION,I)
        SUMMTP=SUMMTP+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        SUMI2(IION,ISTRA)=ERFPIO(IION,I)
        SUMMEP=SUMMEP+SUMI2(IION,ISTRA)
3104  CONTINUE
C
      DO 3125 N=1,NSIGSI
        IF (IIHW(N).EQ.18) THEN
          DO 3126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 3126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
3126      CONTINUE
        ELSEIF (IIHW(N).EQ.43) THEN
          DO 3626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 3626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
3626      CONTINUE
        ENDIF
3125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.24 AND NO.49
C
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=PRFPPHT(IPHOT,I)
        SUMMTP=SUMMTP+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        SUMPH2(IPHOT,ISTRA)=ERFPPHT(IPHOT,I)
        SUMMEP=SUMMEP+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.24) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.49) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0) THEN
        CALL LEER(1)
        WRITE (6,*) 'NO FLUXES RE-EMITTED FROM INCIDENT BULK IONS '
        CALL LEER(1)
      ELSE
        CALL LEER(1)
C
        WRITE (6,*) 'FLUX RE-EMITTED FROM INCIDENT BULK IONS:'
        WRITE (6,*) '(RECYCLING SOURCE) '
C  SURFACE AVERAGED TALLY NO. 6
        IF (SUMA.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 31
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA, 
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 12
        IF (SUMM.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 37
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 18
        IF (SUMI.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 43
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 24
        IF (SUMPH.NE.0.D0) THEN
        WRITE (6,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 49
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        WRITE (6,*) 'TOTAL REEMITTED "ATOMIC" FLUXES, AMPERE AND WATT'
        WRITE (6,*) 'CONTRIB. FROM INCIDENT BULK IONS (REC. SOURCE)'
        CALL MASR1 ('TOT.PFLX',SUMMTP)
        CALL MASR1 ('TOT.EFLX',SUMMEP)
C
      ENDIF
C
      IF (ABS(SUMMTA)+ABS(SUMMTM)+ABS(SUMMTI)+
     .    ABS(SUMMTPH).NE.0.D0) THEN
        CALL LEER (1)
C
        IF (ILIIN(I).GT.0) THEN
          WRITE (6,*) 'TOTAL REEMITTED "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (6,*) '(EXCLUDING CONTRIB. FROM INCIDENT BULK IONS)'
          CALL MASR1 ('TOT.PFLX',SUMMTA+SUMMTM+SUMMTI+SUMMTPH)
          CALL MASR1 ('TOT.EFLX',SUMMEA+SUMMEM+SUMMEI+SUMMEPH)
          CALL LEER (1)
C
        ELSEIF (ILIIN(I).LT.0.AND.ILIIN(I).NE.-3) THEN
          WRITE (6,*) 'TOTAL NEGATIVE "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (6,*) '(EXCLUDING CONTRIB. FROM INCIDENT BULK IONS)'
          CALL MASR1 ('TOT.PFLX',SUMMTA+SUMMTM+SUMMTI+SUMMTPH)
          CALL MASR1 ('TOT.EFLX',SUMMEA+SUMMEM+SUMMEI+SUMMEPH)
        ENDIF
      ENDIF
C
C
C   SPUTTERED FLUXES
C
C   SURFACE AVERAGED TALLY NO.51
C
      SUMMT=0.
      SUMMA=0.
      DO 202 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        SUMA1(IATM,ISTRA)=SPTAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
202     CONTINUE
C
      DO 221 N=1,NSIGSI
        IF (IIHW(N).EQ.51) THEN
          DO 222 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 222
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
222       CONTINUE
        ENDIF
221   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.52
C
      SUMMM=0.
      DO 203 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        SUMM1(IMOL,ISTRA)=SPTML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
203   CONTINUE
C
      DO 223 N=1,NSIGSI
        IF (IIHW(N).EQ.52) THEN
          DO 224 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 224
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
224       CONTINUE
        ENDIF
223   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.53
C
      SUMMI=0.
      DO 204 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        SUMI1(IION,ISTRA)=SPTIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
204   CONTINUE
C
      DO 225 N=1,NSIGSI
        IF (IIHW(N).EQ.53) THEN
          DO 226 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 226
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
226       CONTINUE
        ENDIF
225   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.54
C
      SUMMPH=0.
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        SUMPH1(IPHOT,ISTRA)=SPTPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMI=SUMMI+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.54) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE AVERAGED TALLY NO.55
C
      SUMMP=0.
      DO 205 IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        SUMP1(IPLS,ISTRA)=SPTPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
205   CONTINUE
C
      DO 227 N=1,NSIGSI
        IF (IIHW(N).EQ.55) THEN
          DO 228 IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) GOTO 228
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
228       CONTINUE
        ENDIF
227   CONTINUE
C
      IF (SUMMT.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO FLUXES SPUTTERED FROM THIS SURFACE         ')
        CALL LEER(1)
        GOTO 206
      ENDIF
      CALL LEER(1)
      WRITE (6,*) 'FLUX SPUTTERED FROM SURFACE:'
      IF (SUMMA.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 51
        CALL MASYR1('ATOMS    ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
      ENDIF
      IF (SUMMM.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 52
        CALL MASYR1('MOLECULES',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
      ENDIF
      IF (SUMMI.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 53
        CALL MASYR1('TEST IONS',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
      ENDIF
      IF (SUMMPH.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 54
        CALL MASYR1('PHOTONS',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
      ENDIF
      IF (SUMMP.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 55
        CALL MASYR1('BULK IONS',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .               TEXTS(NSPAMI+1))
        CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .               TEXTS(NSPAMI+1))
      ENDIF
      CALL LEER (1)
      CALL MASR1 ('TOT. FLX',SUMMT)

206   CONTINUE
C
C   ADDITIONAL SURFACE AVERAGED TALLIES
C
C   SURFACE AVERAGED TALLY NO. 57
C
      SUMMS=0.
      DO 302 IADS=1,NADSI
        TEXTA(IADS)=TXTSPW(IADS,NTLSA)
        LGVARS(IADS,ISTRA)=.FALSE.
        SUMS(IADS,ISTRA)=ADDS(IADS,I)
        LOGADS(IADS,ISTRA)=ADDS(IADS,I).NE.0.
        SUMMS=SUMMS+SUMS(IADS,ISTRA)
302   CONTINUE
C
      DO 321 N=1,NSIGSI
        IF (IIHW(N).NE.57) GOTO 321
        IF (IGHW(N).EQ.0) GOTO 321
        DO 322 IADS=1,NADSI
          IF (IGHW(N).NE.IADS) GOTO 322
          VARS(IADS,ISTRA)=SIGMAW(N,I)
          LGVARS(IADS,ISTRA)=LOGADS(IADS,ISTRA)
322     CONTINUE
321   CONTINUE
C
      IF (SUMMS.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO ADDITIONAL SURFACE TALLIES AT THIS SURFACE ')
        CALL LEER(1)
        GOTO 305
      ENDIF
      CALL LEER(1)
      WRITE (6,*) 'ADDITIONAL SURFACE TALLIES:'
      IF (SUMMS.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 57
        CALL MASYR1('ADD.TALLY',SUMS,LOGADS,ISTRA,0,NADS,0,NSTRA,TEXTA)
      ENDIF
      CALL LEER (1)
C
305   CONTINUE
C
C
C   ALGEBRAIC SURFACE AVERAGED TALLIES
C
C   SURFACE AVERAGED TALLY NO. 58
C
      SUMMS=0.
      DO 402 IALS=1,NALSI
        TEXTL(IALS)=TXTSPW(IALS,NTLSR)
        LGVARL(IALS,ISTRA)=.FALSE.
        SUML(IALS,ISTRA)=ALGS(IALS,I)
        LOGALS(IALS,ISTRA)=ALGS(IALS,I).NE.0.
        SUMMS=SUMMS+SUML(IALS,ISTRA)
402   CONTINUE
C
      DO 421 N=1,NSIGSI
        IF (IIHW(N).NE.58) GOTO 421
        IF (IGHW(N).EQ.0) GOTO 421
        DO 422 IALS=1,NALSI
          IF (IGHW(N).NE.IALS) GOTO 422
          VARL(IALS,ISTRA)=SIGMAW(N,I)
          LGVARL(IALS,ISTRA)=LOGALS(IALS,ISTRA)
422     CONTINUE
421   CONTINUE
C
      IF (SUMMS.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO ALGEBRAIC SURFACE TALLIES AT THIS SURFACE ')
        CALL LEER(1)
        GOTO 405
      ENDIF
      CALL LEER(1)
      WRITE (6,*) 'ALGEBRAIC SURFACE TALLIES:'
      DO IALS=1,NALSI
        WRITE (6,*) 'NO ',IALS,'  ',TXTTLW(IALS,NTLSR)
      ENDDO
      IF (SUMMS.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 58
        CALL MASYR1('         ',SUML,LOGALS,ISTRA,0,NALS,0,NSTRA,TEXTL)
      ENDIF
      CALL LEER (1)
C
405   CONTINUE
C

C  SPECTRA
      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO ISPC=1,NADSPC
        IF (ESTIML(ISPC)%PSPC%ISPCSRF == I) THEN
          CALL LEER (1)
          WRITE (6,'(A33)') ' SPECTRUM CALCULATED FOR SURFACE '
          IT = ESTIML(ISPC)%PSPC%ISPCTYP
          IF (IT == 1) THEN
            WRITE (6,'(A20,A40)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (6,'(A20,A40)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
          END IF
          WRITE (6,'(A20,A8)') ' TYPE OF PARTICLE : ',
     .                       TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
          IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
            WRITE (6,'(A10,10X,A16)') ' SPECIES :','SUM OVER SPECIES'
          ELSE
            WRITE (6,'(A10,10X,A8)') ' SPECIES :',
     .             TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .                   ESTIML(ISPC)%PSPC%IPRSP)
          END IF
          WRITE (6,'(A22,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .           ESTIML(ISPC)%PSPC%SPCINT
          IF (NSIGI_SPC > 0)
     .      WRITE (6,'(A22,ES12.4)') ' STANDARD DEVIATION   ',
     .           ESTIML(ISPC)%PSPC%SGMS
        END IF
      END DO

      PRINTED(I) = .TRUE.
10000 CONTINUE
      RETURN
9999  FORMAT (1X,A32)
      END
