*DK HALPHA
      SUBROUTINE HALPHA (IST,IAD1,IAD2,IAD3,IAD4,IAD5,IADS,IADF)
c
C  SUBROUTINE FOR H-ALPHA EMISSIVITY. (BALMER SERIES)
C  CALLED FROM EIRENE, SECTION DIAGNO, SUBR. SIGAL
C  THE H-ALPHA EMISSIVITY PROFILE (PHOTONS/S/CM**3) IS COMPUTED
C  AND WRITTEN ONTO TALLIES ADDV(IAD1,...),... FOR STRATUM NO. IST
C  IAD1: CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  IAD2: CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  IAD3: CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  IAD4: CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  IAD5: CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  IADF: FULCHER
C  IADS: SUM OVER ALL CONTRINUTIONS
c
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE COMSIG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMPRT
      USE COMSOU
      USE CLGIN
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IAD1, IAD2, IAD3, IAD4, IAD5, IADS, IST,
     .                       IADF
      REAL(DP) :: DA31(0:8,0:8),
     .          DP31(0:8,0:8),
     .          DM31(0:8,0:8),
     .          DI31(0:8,0:8),
     .          DN31(0:8,0:8),
     .          DNFUL(0:8,0:8)
      REAL(DP) :: RHMH2(0:8), RH2PH2(0:8,0:8)
      REAL(DP) :: DATM3, DION3, DMOL3, DPLS3, DNML3, RATIO7, TEI, DEJ,
     .          SIGADD1, SIGADD2, SIGADD3, SIGADD4, SIGADD5, TEF, DEF,
     .          RHMP2, DA, DM, FAC43, FAC53, POWALF1, DNFU3, POWALFF,
     .          DPP, DI, DN, RATIO2, SIGADD, FAC21, FAC31, FAC41, FAC51,
     .          FAC32, FAC42, FAC52, POWALF, POWALF2, POWALF3, POWALF4,
     .          POWALF5, DE, TE, FACFUL, SIGADDF
      INTEGER :: IRC, IFIRST, NCELC, IERROR, IR, I, J
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      CHARACTER(8) :: FILNAM
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(6) :: CISTRA
C
      SAVE
C
      DATA IFIRST/0/
C  RADIATIVE TRANSITION RATES (1/S)
C  BALMER ALPHA
      FAC32=4.410E7
C  BALMER BETA
      FAC42=8.419E6
C  BALMER GAMMA
      FAC52=2.530E6
C
C  LYMAN ALPHA
      FAC21=4.699E8
C  LYMAN BETA
      FAC31=5.575E7
C  LYMAN GAMMA
      FAC41=1.278E7
C  LYMAN DELTA
      FAC51=4.125E6
C
C  PASCHEN ALPHA
      FAC43=8.986E6
C  PASCHEN BETA
      FAC53=2.201E6
C
C  FULCHER
      FACFUL=2.53E7      
C
C  INITIALIZE ATOMIC DATA ARRAYS
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        IERROR=0
C
C  READ REDUCED POPULATION COEFFICIENT FOR HYDR. ATOMS FROM FILE AMJUEL
C  AND PUT THEM FROM CREAC(..,..,IR) ONTO DA,DPP,DM,DI, AND DN ARRAY
C
        IR=NREACI+1
        IF (IR.GT.NREAC) THEN
          WRITE (6,*) 'FROM SUBROUTINE HALPHA: '
          CALL MASPRM('NREAC',5,NREAC,'IR',2,IR,IERROR)
          CALL EXIT_OWN(1)
        ENDIF
C
        FILNAM='AMJUEL  '
        H123='H.12'
        CRC='OT '
C
C  H(n=3)/H(n=1)
        REAC='2.1.5a   '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO J=1,9
          DO I=1,9
            DA31(J-1,I-1)=CREAC(J,I,NREACI+1)
          ENDDO
        ENDDO
C  H(n=3)/H+
        REAC='2.1.8a   '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO J=1,9
          DO I=1,9
            DP31(J-1,I-1)=CREAC(J,I,NREACI+1)
          ENDDO
        ENDDO
C  H(n=3)/H2(g)
        REAC='2.2.5a   '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO J=1,9
          DO I=1,9
            DM31(J-1,I-1)=CREAC(J,I,NREACI+1)
          ENDDO
        ENDDO
C  H(n=3)/H2+(g)
        REAC='2.2.14a  '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO J=1,9
          DO I=1,9
            DI31(J-1,I-1)=CREAC(J,I,NREACI+1)
          ENDDO
        ENDDO
C  H(n=3)/H-
        REAC='7.2a     '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO J=1,9
          DO I=1,9
            DN31(J-1,I-1)=CREAC(J,I,NREACI+1)
          ENDDO
        ENDDO
C  H2(n=3,D)/nH2
        REAC='2.2.5fu  '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO J=1,9
          DO I=1,9
            DNFUL(J-1,I-1)=CREAC(J,I,NREACI+1)
          ENDDO
        ENDDO
C
C  NOW READ RATIO OF DENSITIES:
C
C  FIRST: H-/H2
        FILNAM='AMJUEL  '
        H123='H.11'
        REAC='7.0c     '
        CRC='OT '
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO I=1,9
          RHMH2(I-1)=CREAC(I,1,NREACI+1)
        ENDDO

C  NEXT : H2+/H2
        FILNAM='AMJUEL  '
        H123='H.12'
        REAC='2.0c     '
        CRC='OT '
C  2.0c INCLUDES ION CONVERION (CX) ON H2(V)
C  OLD VERSION (WITHOUT THIS CX) SHOULD BE RECOVERED BY
C  READING 2.0b INSTEAD, AND OMITTING THE H- CHANNEL 5.
C
CDR: IF CX ON H2(V) IS NOT INCLUDED IN A SPECIFIC NEUTRAL TRANSPORT EQUATION
CDR  (SEE INPUT BLOCK 4), THEN IT SHOULD NOT BE INCLUDED HERE EITHER.
c slmod begin
c...    In Dgamma as well:
        WRITE(0,*) 'NOTE: REAC=2.0b IS OFF, ASK DETLEV'
c        REAC='2.0b    '
c slmod end
CDR
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
        DO I=1,9
          DO J=1,9
            RH2PH2(I-1,J-1)=CREAC(I,J,NREACI+1)
          ENDDO
        ENDDO
      ENDIF
C
C  END OF INITIALIZATION
C
      IF (IESTR.EQ.IST) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (6,*) 'ERROR IN HALPHA: DATA FOR STRATUM ISTRA= ',IST
        WRITE (6,*) 'ARE NOT AVAILABLE. H-ALPHA ABANDONNED'
        RETURN
      ENDIF
C
C  LOOP OVER 2D MESH
C
      POWALF=0.
      POWALF1=0.
      POWALF2=0.
      POWALF3=0.
      POWALF4=0.
      POWALF5=0.
      POWALFF=0.

      IF (MAX(IAD1,IAD2,IAD3,IAD4,IAD5,IADS,IADF) > NADV) GOTO 999
      ADDV(IAD1,1:NRTAL) = 0.D0
      ADDV(IAD2,1:NRTAL) = 0.D0
      ADDV(IAD3,1:NRTAL) = 0.D0
      ADDV(IAD4,1:NRTAL) = 0.D0
      ADDV(IAD5,1:NRTAL) = 0.D0
      ADDV(IADF,1:NRTAL) = 0.D0
      ADDV(IADS,1:NRTAL) = 0.D0

      DO 1000 NCELL=1,NSBOX
C
C  LOCAL PLASMA DATA
C
        IF (LGVAC(NCELL,NPLS+1)) THEN
          TE=TVAC
          DE=DVAC
          NCELC=NCLTAL(NCELL)
        ELSE
          TE=TEIN(NCELL)
          DE=DEIN(NCELL)
          NCELC=NCLTAL(NCELL)
        ENDIF
        SIGADD1=0.
        SIGADD2=0.
        SIGADD3=0.
        SIGADD4=0.
        SIGADD5=0.
        SIGADDF=0.
!pb        IF (LGVAC(NCELL,0)) GOTO 500
        IF (LGVAC(NCELL,NPLS+1)) GOTO 500
        DEF=LOG(DE*1.D-8)
        TEF=LOG(TE)
        DATM3=0.
        DPLS3=0.
        DMOL3=0.
        DION3=0.
        DNML3=0.
        DNFU3=0.
        DO 150 J=0,8
          DEJ=DEF**J
          DO 150 I=0,8
            TEI=TEF**I
            DATM3=DATM3+DA31(I,J)*TEI*DEJ
            DPLS3=DPLS3+DP31(I,J)*TEI*DEJ
            DMOL3=DMOL3+DM31(I,J)*TEI*DEJ
            DION3=DION3+DI31(I,J)*TEI*DEJ
            DNML3=DNML3+DN31(I,J)*TEI*DEJ
            DNFU3=DNFU3+DNFUL(I,J)*TEI*DEJ
150     CONTINUE
        DATM3=EXP(DATM3)
        DPLS3=EXP(DPLS3)
        DMOL3=EXP(DMOL3)
        DION3=EXP(DION3)
        DNML3=EXP(DNML3)
        DNFU3=EXP(DNFU3)

C  RATIO OF DENSITIES: H- TO H2, COLL. EQUIL. IN VIBRATION
C  (ONLY TE-DEPENDENT)

        RATIO7=0
        DO 160 I=0,8
          TEI=TEF**I
          RATIO7=RATIO7+RHMH2(I)*TEI
160     CONTINUE
        RATIO7=EXP(RATIO7)

C  RATIO OF DENSITIES: H2+ TO H2, INCL. ION CONVERSION
C  RATIO OF DENSITIES: H2+ TO H2, EXCL. ION CONVERSION

        RATIO2=0
        DO 170 J=0,8
          DEJ=DEF**J
          DO 170 I=0,8
            TEI=TEF**I
            RATIO2=RATIO2+RH2PH2(I,J)*TEI*DEJ
170     CONTINUE
        RATIO2=EXP(RATIO2)

C
C  CHANNEL 1
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENA (IONIZATION)
C
        DO 200 IATM=1,NATMI
          IF (NCHARA(IATM).NE.1) GOTO 200
          DA=DATM3*PDENA(IATM,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD1=SIGADD1+DA*FAC32
200     CONTINUE
C
C  CHANNEL 2
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN DIIN (RECOMBINATION)
C
        DO 205 IPLS=1,NPLSI
          IF (NCHARP(IPLS).NE.1.OR.NCHRGP(IPLS).NE.1) GOTO 205
          DPP=DPLS3*DIIN(IPLS,NCELL)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD2=SIGADD2+DPP*FAC32
205     CONTINUE
C
C  CHANNEL 3
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENM (DISSOCIATION OF H2)
C
        DO 210 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 210
          DM=DMOL3*PDENM(IMOL,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD3=SIGADD3+DM*FAC32
210     CONTINUE
C
C  CHANNEL 4
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENI: (DISSOCIATION OF H2+)
C
C  REVISED: USE (PDENM * DENSITY RATIO H2+/H2) NOW, INSTEAD OF PDENI
C
C       DO 215 IION=1,NIONI
C         IF (NCHARI(IION).NE.2) GOTO 215
C         DI=DION3*PDENI(IION,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
C         SIGADD4=SIGADD4+DI*FAC32
C215     CONTINUE
        DO 215 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 215
          DI=DION3*PDENM(IMOL,NCELC)*RATIO2
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD4=SIGADD4+DI*FAC32
215     CONTINUE
C
C  CHANNEL 5
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN H- DENSITY (CHARGE EXCHANGE RECOMBINATION)
C
C  REVISED: USE (PDENM * DENSITY RATIO H-/H2) NOW, INSTEAD OF PDENI
C
        DO 220 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 220
          DN=DNML3*PDENM(IMOL,NCELC)*RATIO7
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD5=SIGADD5+DN*FAC32
220     CONTINUE
C
C
C  CHANNEL 6
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENM (FULCHER)
C
        DO IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) CYCLE
          DM=DNFU3*PDENM(IMOL,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADDF=SIGADDF+DM*FACFUL
        END DO  
C
500     CONTINUE
C
        SIGADD=SIGADD1+SIGADD2+SIGADD3+SIGADD4+SIGADD5
C
        ADDV(IAD1,NCELC)=ADDV(IAD1,NCELC)+SIGADD1*VOL(NCELL)
        ADDV(IAD2,NCELC)=ADDV(IAD2,NCELC)+SIGADD2*VOL(NCELL)
        ADDV(IAD3,NCELC)=ADDV(IAD3,NCELC)+SIGADD3*VOL(NCELL)
        ADDV(IAD4,NCELC)=ADDV(IAD4,NCELC)+SIGADD4*VOL(NCELL)
        ADDV(IAD5,NCELC)=ADDV(IAD5,NCELC)+SIGADD5*VOL(NCELL)
        ADDV(IADF,NCELC)=ADDV(IADF,NCELC)+SIGADDF*VOL(NCELL)
        ADDV(IADS,NCELC)=ADDV(IADS,NCELC)+SIGADD*VOL(NCELL)
C
C
        POWALF1=POWALF1+SIGADD1*3.028E-19*VOL(NCELL)
        POWALF2=POWALF2+SIGADD2*3.028E-19*VOL(NCELL)
        POWALF3=POWALF3+SIGADD3*3.028E-19*VOL(NCELL)
        POWALF4=POWALF4+SIGADD4*3.028E-19*VOL(NCELL)
        POWALF5=POWALF5+SIGADD5*3.028E-19*VOL(NCELL)
        POWALFF=POWALFF+SIGADDF*3.028E-19*VOL(NCELL)
C
        POWALF =POWALF +SIGADD *3.028E-19*VOL(NCELL)
C
1000  CONTINUE

      ADDV(IAD1,1:NSBOX_TAL)=ADDV(IAD1,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD2,1:NSBOX_TAL)=ADDV(IAD2,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD3,1:NSBOX_TAL)=ADDV(IAD3,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD4,1:NSBOX_TAL)=ADDV(IAD4,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD5,1:NSBOX_TAL)=ADDV(IAD5,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IADF,1:NSBOX_TAL)=ADDV(IADF,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IADS,1:NSBOX_TAL)=ADDV(IADS,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)

      CALL LEER(2)
      CALL FTCRI(IST,CISTRA)
      IF (IST.GT.0)
     .CALL MASBOX('SUBR. HALPHA CALLED, FOR STRATUM NO. '//CISTRA)
      IF (IST.EQ.0)
     .CALL MASBOX('SUBR. HALPHA CALLED, FOR SUM OVER STRATA')
      CALL LEER(1)
      WRITE (6,*) ' RADIATED POWER BY HALPHA:',POWALF
      WRITE (6,*) ' COUPL. TO GROUNDSTATE   :',POWALF1
      WRITE (6,*) ' COUPLING TO CONTINUUM   :',POWALF2
      WRITE (6,*) ' COUPLING TO MOLECULES   :',POWALF3
      WRITE (6,*) ' COUPLING TO MOL.IONS    :',POWALF4
      WRITE (6,*) ' COUPLING TO NEG.IONS    :',POWALF5
      WRITE (6,*) ' COUPLING TO FULCHER     :',POWALFF
      CALL LEER(2)

      CALL INTTAL (ADDV,VOLTAL,IAD1,NADV,NSBOX_TAL,ADDVI(IAD1,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL INTTAL (ADDV,VOLTAL,IAD2,NADV,NSBOX_TAL,ADDVI(IAD2,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL INTTAL (ADDV,VOLTAL,IAD3,NADV,NSBOX_TAL,ADDVI(IAD3,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL INTTAL (ADDV,VOLTAL,IAD4,NADV,NSBOX_TAL,ADDVI(IAD4,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL INTTAL (ADDV,VOLTAL,IAD5,NADV,NSBOX_TAL,ADDVI(IAD5,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL INTTAL (ADDV,VOLTAL,IADF,NADV,NSBOX_TAL,ADDVI(IADF,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL INTTAL (ADDV,VOLTAL,IADS,NADV,NSBOX_TAL,ADDVI(IADS,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
C
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU)
        WRITE (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (6,*) 'WRITE 11  IRC= ',IRC
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU)
        WRITE (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (6,*) 'WRITE 11  IRC= ',IRC
      ENDIF
C
      RETURN
999   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. HALPHA '
      WRITE (6,*) 'NO STORAGE AVAILBALE ON ADDITIONAL TALLY ADDV '
      WRITE (6,*) 'STORAGE REQUESTED FOR IADV= ',IAD1,IAD2,IAD3,IAD4,
     .             IAD5,IADS,IADF
      WRITE (6,*) 'CHECK INPUT BLOCK 10A '
      CALL EXIT_OWN(1)
      END