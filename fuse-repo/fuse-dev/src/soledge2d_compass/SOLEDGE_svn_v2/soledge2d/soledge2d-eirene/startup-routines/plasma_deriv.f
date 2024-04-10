c  new in 2004:
c  density models to contruct backgound data from other given data :
c      Saha, Boltzmann, Corona, Colrad, File (fort.13,or: fort.10)
c
c  presently:  "File" and "Boltzmann": may affect electron density.
c              hence: done prior to electron density, etc...
c              "Corona", "Colrad", "Saha": need electron density as
c                            input, or, at least, do not affect n_e
c                            hence: done after electron density, etc...
C  may05
c  1) additional density model only for neutrals?  removed
c  2) boltzmann factor only if Ti gt tvac
c     to be checked: correct low T limit: everything in lower level?
c  3) new ti only if nlmlti
c  4) new vi only if nlmlv (nlmlv is new, in cinit.f, and set in input.f)
c  june05
c     bug fix: in colrad model: vnew=max(vvac,vold) nonsense, because
c          vold<0 possible. replaced by vnew=vold
c
c  march 06
c     new option: icall > 0, and call base_density
c       allows to use output tallies and special "density model" to
c       construct new input tallies (densities, temperatures, drift velocities)
c       e.g. for post processing (diagno), or for iterations.
c
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  06.03.07: new density models 'CONSTANT' and 'MULTIPLY' introduced
!pb            'CONSTANT' sets constant plasma profiles
!pb            'MULTIPLY' creates a new bulkdensity by multiplying an
!pb            existing plasma density with a factor specified in input block 5
 
c
      SUBROUTINE EIRENE_PLASMA_DERIV (ICALL)
c  input:
c    icall=0
c      called prior to Monte Carlo Loop
c      in this call all density models referring to output tallies
c      are ignored.
c    icall=1
c      called after to Monte Carlo Loop and sum over strata
c        this allows to put output tallies from a run onto the
c        background for a next iteration or post processing.
c        in this call all density models referring to input tallies  are
c        ignored.
c      write fort.13 after all density models are done.
c  set derived plasma parameters:
c   DEIN             : electron density (from quasineutrality)
c   DEINL            : log electron density (with cutoffs)
c   TEINL            : log electron temperature (with cutoffs)
c   LGVAC(...,NPLS+1): electron vacuum flag
c   LGVAC(...,IPLS)  : bulk ion vacuum flag
c   LGVAC(...,0)     : background vacuum flag
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CTRCEI
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEI
      USE EIRMOD_COMXS
      USE EIRMOD_CESTIM
      use EIRMOD_csdvi
      use EIRMOD_csdvi_bgk
      use EIRMOD_csdvi_cop
      use EIRMOD_comsou
      use EIRMOD_cspei
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ICALL
      REAL(DP) :: ZTII, ZTNI, FCT2, FCRG, FCT1, EIRENE_VDION, ZTEI, 
     .            ZTNE,EMPLS, FCT0, TEPLS, DEPLS, DIPLS, AM1, TEF, DEF,
     .            TEI, DEJ, TVACL, DVACL, BOLTZFAC, RCORONA, RCOLRAD,
     .            TEIDEJ, EIRENE_RATE_COEFF, RCMIN, RCMAX, ERATE
      REAL(DP) :: tpb1, tpb2, EIRENE_second_own
      REAL(DP) :: COEF(0:8), COEF2D(0:8,0:8), FP(6)
      REAL(DP), ALLOCATABLE :: DEINTF(:), SUMNI(:), SUMMNI(:),
     .                         BASE_DENSITY(:), BASE_TEMP(:)
      INTEGER :: I, IR, IN, IP, IPM, J, IPLS, IOLD, ISW, IRE, I1,
     .           IO, IPLSTI, IPLSV, IOLDTI, IOLDV, IBS, JFEXMN, JFEXMX
 
      TYPE(EIRENE_SPECTRUM), POINTER :: SPEC
      LOGICAL :: FOUND
 
      FP = 0._DP
      RCMIN = -HUGE(1._DP)
      RCMAX =  HUGE(1._DP)
      JFEXMN = 0
      JFEXMX = 0
 
      tpb1 = EIRENE_second_own()
 
      IBS = 0
      DO IPLS=1,NPLSI
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF (INDEX(CDENMODEL(IPLS),'FORT.13') > 0) THEN
          CALL EIRENE_ALLOC_BCKGRND
          ALLOCATE(DEINTF(NRAD))
          OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
          REWIND 13
          READ (13,IOSTAT=IO) TEINTF,TIINTF,DEINTF,DIINTF,
     .                        VXINTF,VYINTF,VZINTF
          IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMUSR, IO= ',IO
          CLOSE (UNIT=13)
          IF (IO.EQ.0) THEN
            IOLD=TDMPAR(IPLS)%TDM%ISP(1)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)
            IF (NLMLTI) TIIN(IPLSTI,:)=TIINTF(IOLDTI,:)
            DIIN(IPLS,:)=DIINTF(IOLD,:)
            IF (NLMLV) THEN
              VXIN(IPLSV,:)=VXINTF(IOLDV,:)
              VYIN(IPLSV,:)=VYINTF(IOLDV,:)
              VZIN(IPLSV,:)=VZINTF(IOLDV,:)
            END IF
          ENDIF
          DEALLOCATE(DEINTF)
        ELSEIF (INDEX(CDENMODEL(IPLS),'FORT.10') > 0) THEN
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)
 
          ALLOCATE (BASE_DENSITY(NRAD))
          ALLOCATE (BASE_TEMP(NRAD))
          CALL EIRENE_GET_BASE_DENSITY(1)
 
          IF (NLMLTI) TIIN(IPLSTI,:)=BASE_TEMP(:)
          IF (NLMLV) THEN
            VXIN(IPLSV,:)=VXIN(IOLDV,:)
            VYIN(IPLSV,:)=VYIN(IOLDV,:)
            VZIN(IPLSV,:)=VZIN(IOLDV,:)
          END IF
 
          DO IR=1,NSBOX
            DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR))
          ENDDO
          DEALLOCATE (BASE_DENSITY)
          DEALLOCATE (BASE_TEMP)
 
          IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
            FOUND = .TRUE.
            DO IR=1,NSBOX
              IF (LSPCCLL(IR)) THEN
                IF (FOUND) ALLOCATE (SPEC)
                CALL EIRENE_GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  BACK_SPEC(IBS)%PSPC => SPEC
                END IF
              END IF
            ENDDO
            IF (.NOT.FOUND) DEALLOCATE(SPEC)
          END IF
 
        ELSEIF (INDEX(CDENMODEL(IPLS),'CONSTANT') > 0) THEN
 
          IF (NLMLTI) TIIN(IPLSTI,:)=MAX(TVAC,TDMPAR(IPLS)%TDM%TVAL)
          DIIN(IPLS,:)=MAX(DVAC,TDMPAR(IPLS)%TDM%DVAL)
          IF (NLMLV) THEN
            VXIN(IPLSV,:)=TDMPAR(IPLS)%TDM%VXVAL
            VYIN(IPLSV,:)=TDMPAR(IPLS)%TDM%VYVAL
            VZIN(IPLSV,:)=TDMPAR(IPLS)%TDM%VZVAL
          END IF
 
        ELSEIF (INDEX(CDENMODEL(IPLS),'MULTIPLY') > 0) THEN
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)
 
          ALLOCATE (BASE_DENSITY(NRAD))
          ALLOCATE (BASE_TEMP(NRAD))
          CALL EIRENE_GET_BASE_DENSITY(1)
 
          DO IR=1,NSBOX
            IF (NLMLTI) TIIN(IPLSTI,IR)=MAX(TVAC,BASE_TEMP(IR)*
     .                             TDMPAR(IPLS)%TDM%TFACTOR)
            DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR)*
     .                             TDMPAR(IPLS)%TDM%DFACTOR)
            IF (NLMLV) THEN
              VXIN(IPLSV,IR)=VXIN(IOLDV,IR)*TDMPAR(IPLS)%TDM%VFACTOR
              VYIN(IPLSV,IR)=VYIN(IOLDV,IR)*TDMPAR(IPLS)%TDM%VFACTOR
              VZIN(IPLSV,IR)=VZIN(IOLDV,IR)*TDMPAR(IPLS)%TDM%VFACTOR
            END IF
          ENDDO
 
          DEALLOCATE (BASE_DENSITY)
          DEALLOCATE (BASE_TEMP)
 
        ELSEIF (INDEX(CDENMODEL(IPLS),'BOLTZMANN') > 0) THEN
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)
 
          ALLOCATE (BASE_DENSITY(NRAD))
          ALLOCATE (BASE_TEMP(NRAD))
          CALL EIRENE_GET_BASE_DENSITY(1)
 
          IF (NLMLTI) TIIN(IPLSTI,:)=BASE_TEMP(:)
          IF (NLMLV) THEN
            VXIN(IPLSV,:)=VXIN(IOLDV,:)
            VYIN(IPLSV,:)=VYIN(IOLDV,:)
            VZIN(IPLSV,:)=VZIN(IOLDV,:)
          END IF
 
          FOUND = .TRUE.
          DO IR=1,NSBOX
            BOLTZFAC=0._DP
            IF (TIIN(IOLDTI,IR).GT.TVAC)
     .        BOLTZFAC=TDMPAR(IPLS)%TDM%G_BOLTZ*
     .                 EXP(-TDMPAR(IPLS)%TDM%DELTAE/BASE_TEMP(IR))
            DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR))*BOLTZFAC
 
            IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
              IF (LSPCCLL(IR)) THEN
                IF (FOUND) ALLOCATE (SPEC)
                CALL EIRENE_GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%SPC = SPEC%SPC * BOLTZFAC
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  BACK_SPEC(IBS)%PSPC => SPEC
                END IF
              END IF
            END IF
 
          ENDDO
          IF (.NOT.FOUND) DEALLOCATE(SPEC)
          DEALLOCATE (BASE_DENSITY)
          DEALLOCATE (BASE_TEMP)
        END IF
      END DO
 
cdr   tpb2 = EIRENE_second_own()
cdr   write (6,*) ' cputime for fort.13, fort.10, boltzmann ',tpb2-tpb1
cdr   tpb1 = tpb2
 
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
C  SET ELECTRON-DENSITY FROM QUASI-NEUTRALITY, AND DRIFT ENERGY (EV)
      TVACL=LOG(TVAC)
      DVACL=LOG(DVAC*1.E-8_DP)
      LGVAC=.TRUE.
      DO 5102 J=1,NSBOX
        DEIN(J)=0.
        DO IPLS=1,NPLSI
          DEIN(J)=DEIN(J)+DBLE(NCHRGP(IPLS))*DIIN(IPLS,J)
        END DO
C  SET 'LOG OF TEMPERATURE AND DENSITY' ARRAYS
        ZTEI=MAX(TVAC,MIN(TEIN(J),1.E10_DP))
        TEINL(J)=LOG(ZTEI)
        ZTNE=MAX(DVAC,MIN(DEIN(J),1.E20_DP))
        DEINL(J)=LOG(ZTNE)
        TEPLS=TEIN(J)
        DEPLS=DEIN(J)
        LGVAC(J,NPLS+1)=TEPLS.LE.TVAC.OR.DEPLS.LE.DVAC
        LGVAC(J,0)     =LGVAC(J,0).AND.LGVAC(J,NPLS+1)
5102  CONTINUE
 
cdr      tpb2 = EIRENE_second_own()
cdr      write (6,*) ' cputime for log values ',tpb2-tpb1
cdr      tpb1 = tpb2
 
C
C   SPECIAL DENSITY MODELS
C
      ALLOCATE (BASE_DENSITY(NRAD))
      ALLOCATE (BASE_TEMP(NRAD))
      DO IPLS=1,NPLSI
 
        IF (LEN_TRIM(CDENMODEL(IPLS)) == 0) CYCLE
 
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
 
        SELECT CASE (CDENMODEL(IPLS))
 
        CASE ('SAHA      ')
!PB   TO BE WRITTEN
          WRITE (iunout,*) ' DENSITY PROFILE ACCORDING TO SAHA IS NOT ',
     .                'AVAILABLE '
          WRITE (iunout,*) ' PLEASE CHOOSE DIFFERENT OPTION FOR ION ',
     .                'DENSITY ',IPLS
          CALL EIRENE_EXIT_OWN(1)
 
        CASE ('CORONA    ')
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)
 
          CALL EIRENE_GET_BASE_DENSITY(1)
 
          IF (NLMLTI) TIIN(IPLSTI,:)=BASE_TEMP(:)
          IF (NLMLV) THEN
            VXIN(IPLSV,:)=VXIN(IOLDV,:)
            VYIN(IPLSV,:)=VYIN(IOLDV,:)
            VZIN(IPLSV,:)=VZIN(IOLDV,:)
          END IF
 
          REACDAT(NREACI+1)%LRTC = .FALSE.
          CALL EIRENE_SLREAC (NREACI+1,TDMPAR(IPLS)%TDM%FNAME(1),
     .                 TDMPAR(IPLS)%TDM%H2(1),
     .                 TDMPAR(IPLS)%TDM%REACTION(1),
     .                 TDMPAR(IPLS)%TDM%CR(1),
     .                 RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
          AM1=1._DP/TDMPAR(IPLS)%TDM%A_CORONA
c  now COEF contains the excitation rate,
c  and AMI is the radiative decay rate
c  compute new density of species ipls from equilibrium between
c  these two processes  for the "ground state" density BASE_DENSITY
          FOUND = .TRUE.
          DO IR=1,NSURF
            IF (LGVAC(IR,NPLS+1)) THEN
              TEF=TVACL
            ELSE
              TEF=TEINL(IR)
            END IF
            RCORONA = EIRENE_RATE_COEFF(NREACI+1,TEF,0._DP,.TRUE.,
     .                                  0,ERATE)
            DIIN(IPLS,IR)=BASE_DENSITY(IR)*RCORONA*DEIN(IR)*AM1
            DIIN(IPLS,IR)=MAX(DVAC,DIIN(IPLS,IR))
 
            IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
              IF (LSPCCLL(IR)) THEN
                IF (FOUND) ALLOCATE (SPEC)
                CALL EIRENE_GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  SPEC%SPC = SPEC%SPC * RCORONA*DEIN(IR)*AM1
                  BACK_SPEC(IBS)%PSPC => SPEC
                END IF
              END IF
            END IF
 
          END DO
          IF (.NOT.FOUND) DEALLOCATE(SPEC)
 
cdr      tpb2 = EIRENE_second_own()
cdr      write (6,*) ' cputime for corona ',ipls,tpb2-tpb1
cdr      tpb1 = tpb2
 
        CASE ('COLRAD    ')
          IF (.NOT.ALLOCATED(SUMNI)) THEN
            ALLOCATE (SUMNI(NRAD))
            ALLOCATE (SUMMNI(NRAD))
          END IF
          SUMNI = 0._DP
          SUMMNI = 0._DP
          DIIN(IPLS,:)=0._DP
C  ARE THERE MULTIPLE ION TEMPERATURES?
          if (nlmlti) then
            TIIN(IPLSTI,:)=0._DP
          endif
C  ARE THERE MULTIPLE ION DRIFT VELOCITIES?
          if (nlmlv) then
            VXIN(IPLSV,:)=0._DP
            VYIN(IPLSV,:)=0._DP
            VZIN(IPLSV,:)=0._DP
          endif
          FOUND = .TRUE.
          DO IRE=1,TDMPAR(IPLS)%TDM%NRE
            IOLD=TDMPAR(IPLS)%TDM%ISP(IRE)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)
            CALL EIRENE_GET_BASE_DENSITY(IRE)
            REACDAT(NREACI+1)%LOTH = .FALSE.
            CALL EIRENE_SLREAC (NREACI+1,TDMPAR(IPLS)%TDM%FNAME(IRE),
     .                   TDMPAR(IPLS)%TDM%H2(IRE),
     .                   TDMPAR(IPLS)%TDM%REACTION(IRE),
     .                   TDMPAR(IPLS)%TDM%CR(IRE),
     .                   RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
            I1=INDEX(TDMPAR(IPLS)%TDM%H2(IRE),'.')
            READ (TDMPAR(IPLS)%TDM%H2(IRE)(I1+1:),*) ISW
            SELECT CASE (ISW)
            CASE (11)
c  only temperature dependence in reduced population coefficient
              COEF(0:8)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(1:9,1)
              DO IR=1,NSURF
                IF (LGVAC(IR,NPLS+1)) CYCLE
                TEF=TEINL(IR)
                RCOLRAD=0._DP
                DO I=0,8
                  TEI=TEF**I
                  RCOLRAD=RCOLRAD+COEF(I)*TEI
                END DO
                RCOLRAD=EXP(RCOLRAD)
                DIIN(IPLS,IR)=DIIN(IPLS,IR)+BASE_DENSITY(IR)*RCOLRAD
                if (nlmlti) then
                  TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+
     .                            BASE_DENSITY(IR)*BASE_TEMP(IR)
                endif
                if (nlmlv) then
                  VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VXIN(IOLDV,IR)
                  VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VYIN(IOLDV,IR)
                  VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VZIN(IOLDV,IR)
                endif
                SUMNI(IR) = SUMNI(IR) + BASE_DENSITY(IR)
                SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*BASE_DENSITY(IR)
 
                IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
                  IF (LSPCCLL(IR)) THEN
                    IF (FOUND) ALLOCATE (SPEC)
                    CALL EIRENE_GET_SPECTRUM (IR,IRE,SPEC,FOUND)
                    IF (FOUND) THEN
                      IF (IRE == 1) THEN
                        IBS = IBS + 1
                        SPEC%IPRTYP = 4
                        SPEC%IPRSP = IPLS
                        SPEC%SPC = SPEC%SPC * RCOLRAD
                        BACK_SPEC(IBS)%PSPC => SPEC
                      ELSE
                        IF (  (BACK_SPEC(IBS)%PSPC%NSPC == SPEC%NSPC)
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMIN==SPEC%SPCMIN)
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMAX==SPEC%SPCMAX))
     .                  THEN
                          SPEC%SPC = SPEC%SPC * RCOLRAD
                          BACK_SPEC(IBS)%PSPC%SPC =
     .                         BACK_SPEC(IBS)%PSPC%SPC + SPEC%SPC
                          DEALLOCATE(SPEC)
                        ELSE
                          WRITE (IUNOUT,*) ' ERROR IN PLASMA_DERIV,',
     .                       ' DENSITY MODEL COLRAD '
                          WRITE (IUNOUT,*) ' TWO SPECTRA CONTRIBTING ',
     .                       ' TO THE SAME BACKGROUND SPECTRUM DO NOT',
     .                       ' MATCH '
                          WRITE (IUNOUT,*) ' IPLS, IRE ',IPLS, IRE
                          CALL EIRENE_EXIT_OWN(1)
                        END IF
                      END IF
                    END IF
                  END IF
                END IF
 
              END DO
              IF (.NOT.FOUND) DEALLOCATE(SPEC)
            CASE (12)
c  temperature and density dependence in reduced population coefficient
              COEF2D(0:8,0:8)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(1:9,1:9)
              DO IR=1,NSURF
                IF (LGVAC(IR,NPLS+1)) CYCLE
                TEF=TEINL(IR)
                DEF=LOG(DEIN(IR)*1.E-8_DP)
                RCOLRAD=0._DP
                DEJ=1._DP
                DO J=0,8
                  TEIDEJ=DEJ
                  DO I=0,8
                    RCOLRAD=RCOLRAD+COEF2D(I,J)*TEIDEJ
                    TEIDEJ=TEIDEJ*TEF
                  END DO
                  DEJ=DEJ*DEF
                END DO
                RCOLRAD=EXP(RCOLRAD)
                DIIN(IPLS,IR)=DIIN(IPLS,IR)+BASE_DENSITY(IR)*RCOLRAD
                if (nlmlti) then
                  TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+
     .                            BASE_DENSITY(IR)*BASE_TEMP(IR)
                endif
                if (nlmlv) then
                  VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VXIN(IOLDV,IR)
                  VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VYIN(IOLDV,IR)
                  VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VZIN(IOLDV,IR)
                endif
                SUMNI(IR) = SUMNI(IR) + BASE_DENSITY(IR)
                SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*BASE_DENSITY(IR)
 
                IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
                  IF (LSPCCLL(IR)) THEN
                    IF (FOUND) ALLOCATE (SPEC)
                    CALL EIRENE_GET_SPECTRUM (IR,IRE,SPEC,FOUND)
                    IF (FOUND) THEN
                      IF (IRE == 1) THEN
                        IBS = IBS + 1
                        SPEC%IPRTYP = 4
                        SPEC%IPRSP = IPLS
                        SPEC%SPC = SPEC%SPC * RCOLRAD
                        BACK_SPEC(IBS)%PSPC => SPEC
                      ELSE
                        IF (  (BACK_SPEC(IBS)%PSPC%NSPC == SPEC%NSPC)
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMIN==SPEC%SPCMIN)
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMAX==SPEC%SPCMAX))
     .                  THEN
                          SPEC%SPC = SPEC%SPC * RCOLRAD
                          BACK_SPEC(IBS)%PSPC%SPC =
     .                         BACK_SPEC(IBS)%PSPC%SPC + SPEC%SPC
                          DEALLOCATE(SPEC)
                        ELSE
                          WRITE (IUNOUT,*) ' ERROR IN PLASMA_DERIV,',
     .                       ' DENSITY MODEL COLRAD '
                          WRITE (IUNOUT,*) ' TWO SPECTRA CONTRIBTING ',
     .                       ' TO THE SAME BACKGROUND SPECTRUM DO NOT',
     .                       ' MATCH '
                          WRITE (IUNOUT,*) ' IPLS, IRE ',IPLS, IRE
                          CALL EIRENE_EXIT_OWN(1)
                        END IF
                      END IF
                    END IF
                  END IF
                END IF
 
              END DO  ! IR
              IF (.NOT.FOUND) DEALLOCATE(SPEC)
            CASE DEFAULT
              WRITE (iunout,*) ' H.',ISW,
     .             ' NOT FORESEEN IN COLRAD DENSITY MODEL '
            END SELECT
          END DO   ! IRE
 
          DIIN(IPLS,:)=MAX(DVAC,DIIN(IPLS,:))
          if (nlmlti) then
            TIIN(IPLSTI,:)=MAX(TVAC,TIIN(IPLSTI,:)/(SUMNI(:)+eps60))
          endif
          if (nlmlv) then
            VXIN(IPLSV,:)=VXIN(IPLSV,:)/(SUMMNI(:)+eps60)
            VYIN(IPLSV,:)=VYIN(IPLSV,:)/(SUMMNI(:)+eps60)
            VZIN(IPLSV,:)=VZIN(IPLSV,:)/(SUMMNI(:)+eps60)
          endif
 
        CASE DEFAULT
!  NOTHING TO BE DONE HERE, ALREADY COMPLETED
        END SELECT ! density model
 
cdr      tpb2 = EIRENE_second_own()
cdr      write (6,*) ' cputime for colrad ',ipls,tpb2-tpb1
cdr      tpb1 = tpb2
 
      END DO
      IF (ALLOCATED(SUMNI)) THEN
        DEALLOCATE (SUMNI)
        DEALLOCATE (SUMMNI)
      END IF
      DEALLOCATE (BASE_DENSITY)
      DEALLOCATE (BASE_TEMP)
 
      NBACK_SPEC = IBS
 
cdr      tpb2 = EIRENE_second_own()
cdr      write (6,*) ' cputime for density models ',tpb2-tpb1
cdr      tpb1 = tpb2
 
C
C  special density models done
 
C  SET DRIFT ENERGY (EV)
      DO J=1,NSBOX
        DO IPLS=1,NPLSI
          IPLSV=MPLSV(IPLS)
          IF (NLDRFT) THEN
            IF (INDPRO(4) == 8) THEN
              IF(IPLS.EQ.1) THEN
                EDRIFT(IPLS,J)=CVRSSP(IPLS)*EIRENE_VDION(J)**2
              ELSE
                WRITE(iunout,*)'WARNING! IPLS>1 NO DRIFT!'
                EDRIFT(IPLS,J)=0.D0
              END IF
            ELSE
              EDRIFT(IPLS,J)=CVRSSP(IPLS)*
     .              (VXIN(IPLSV,J)**2+VYIN(IPLSV,J)**2+VZIN(IPLSV,J)**2)
            END IF
          ELSE
            EDRIFT(IPLS,J)=0.D0
          ENDIF
        END DO
      END DO
C
C  SET B_PERP
C
      DO J=1,NSBOX
        IF (ABS(BXIN(J)) > EPS10) THEN
           BYPERP(J) = 1._DP
           BXPERP(J) = -BYIN(J)/BXIN(J)
        ELSEIF (ABS(BYIN(J)) > EPS10) THEN
           BXPERP(J) = 0._DP
           BYPERP(J) = -BXIN(J)/BYIN(J)
        ELSE
           BXPERP(J) = 1._DP
           BYPERP(J) = 0._DP
        END IF
C  CHECK ORIENTATION
        IF (BXIN(J)*BYPERP(J)-BXPERP(J)*BYIN(J) < 0._DP) THEN
           BXPERP(J) = -BXPERP(J)
           BYPERP(J) = -BYPERP(J)
        END IF
      END DO
 
      DO 5103 J=1,NSBOX
C  SET 'VACUUM REGION FLAGS'
C  LGVAC(...,0)  VACUUM, NO REACTION RATES AT ALL
C  LGVAC(...,IPLS)       NO REACTION RATES FOR BACKGROUND SPECIES IPLS
C  LGVAC(...,NPLS+1)     NO REACTION RATES FOR BACKGROUND ELECTRONS
C                        BUT PERHAPS FOR NEUTRAL BACKGROUND
        DO 5106 IPLS=1,NPLSI
          IPLSTI=MPLSTI(IPLS)
          EMPLS=1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J)
          DIPLS=DIIN(IPLS,J)
          LGVAC(J,IPLS)=EMPLS.LE.TVAC.OR.DIPLS.LE.DVAC
          LGVAC(J,0)   =LGVAC(J,0).AND.LGVAC(J,IPLS)
5106    CONTINUE
5103  CONTINUE
C
      DO 5105 IPLS=1,NPLSI
C  FACTOR FOR MOST PROBABLE SPEED
        FCT0=1./RMASSP(IPLS)*2.*CVEL2A*CVEL2A
C  FACTOR FOR MEAN SPEED
        FCT1=1./RMASSP(IPLS)*8./PIA*CVEL2A*CVEL2A
C  FACTOR FOR ROOT MEAN SQUARE SPEED
        FCT2=1./RMASSP(IPLS)*3.*CVEL2A*CVEL2A
        FCRG=CVEL2A/SQRT(RMASSP(IPLS))
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        BVIN(IPLSV,:)=0._DP
        PARMOM(IPLS,:)=0._DP
        DO 5105 J=1,NSBOX
          ZTII=MAX(TVAC,MIN(TIIN(IPLSTI,J),1.E10_DP))
          TIINL(IPLSTI,J)=LOG(ZTII)
          IF (NLDRFT) THEN
            BVIN(IPLSV,J)=BXIN(J)*VXIN(IPLSV,J)+
     .                    BYIN(J)*VYIN(IPLSV,J)+
     .                    BZIN(J)*VZIN(IPLSV,J)
            PARMOM(IPLS,J)=BVIN(IPLSV,J)*SIGN(1._DP,BVIN(IPLSV,J))*
     .                     AMUA*RMASSP(IPLS)
          ENDIF
C
C  ZT1: FOR "EFFECTIVE" PLASMA PARTICLE VELOCITY IN CROSS SECTIONS
C       FOR HEAVY PARTICLE INTERACTIONS
C       SQRT(ZT1) IS THE MEAN VELOCITY AT TI=ZTII, TAKEN AS
C       ROOT MEAN SQUARE SPEED
C
          ZT1(IPLS,J)=FCT2*ZTII
C
C  ZRGQ: VARIANCE FOR SAMPLING FROM MAXWELLIAN VELOCITY DISTRIBUTION
C  ZRG=SQRT(ZRGQ) = STANDARD DEVIATION
C
          ZRG(IPLS,J)=FCRG*SQRT(ZTII)
C
          ZTNI=MAX(DVAC,MIN(DIIN(IPLS,J),1.E20_DP))
          DIINL(IPLS,J)=LOG(ZTNI)
5105  CONTINUE
C
      IF (LEVGEO.EQ.3) THEN
        DO 5161 I=1,NPPLG-1
          DO 5162 IP=NPOINT(2,I),NPOINT(1,I+1)-1
            IPM=IP-1
            DO 5163 IPLS=0,NPLS+1
              DO 5163 IR=1,NR1STM
                IN=IR+IPM*NR1ST
                LGVAC(IN,IPLS)=.TRUE.
5163        CONTINUE
5162      CONTINUE
5161    CONTINUE
      ENDIF
 
!
 
C
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C
 
cdr      tpb2 = EIRENE_second_own()
cdr      write (6,*) ' cputime for edrift, b_perp, etc. ',tpb2-tpb1
cdr      tpb1 = tpb2
 
      IF ((NFILEL >=1) .AND. (NFILEL <=5)) THEN
         NFILEL=3
         CALL EIRENE_WRPLAM(TRCFLE,0)
!      ELSE
      ELSE IF (NFILEL > 5) THEN
         NFILEL=9
         CALL EIRENE_WRPLAM_XDR(TRCFLE,0)
      END IF
 
cdr      tpb2 = EIRENE_second_own()
cdr      write (6,*) ' cputime for wrplam ',tpb2-tpb1
cdr      tpb1 = tpb2
 
 
      RETURN
 
      CONTAINS
 
      SUBROUTINE EIRENE_GET_BASE_DENSITY(IRE)
c  input: ire, number of reaction/density that contributes to the
c              evaluation of the expression for the selected species
c              ipls with special density/temperature option
      INTEGER, INTENT(IN) :: IRE
      INTEGER :: IG, IT, ISTRA, ITYP
 
      BASE_DENSITY = 0._DP
      BASE_TEMP = 0._DP
 
      IF (ICALL > 0) THEN
 
C FOR CALLS AFTER PARTICLE TRACING
 
        ISTRA = TDMPAR(IPLS)%TDM%ISTR(IRE)
        ITYP  = TDMPAR(IPLS)%TDM%ITP(IRE)
        IF (ISTRA.EQ.IESTR.OR.ITYP.EQ.4) THEN
C  NOTHING TO BE DONE
        ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
          IESTR=ISTRA
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .                 NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
          IESTR=ISTRA
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .                 NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSE
          WRITE (6,*)
     .      'ERROR IN PLASMA_DERIV: DATA FOR STRATUM ISTRA= ',ISTRA
          WRITE (6,*) 'ARE NOT AVAILABLE. '
          RETURN
        ENDIF
      END IF   ! ICALL > 0
 
      SELECT CASE (TDMPAR(IPLS)%TDM%ITP(IRE))
      CASE(0)
        IF (ASSOCIATED(PDENPH)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0)  THEN
              BASE_DENSITY(IG) = PDENPH(IOLD,IT)
              BASE_TEMP(IG) = EDENPH(IOLD,IT)/(PDENPH(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(1)
        IF (ASSOCIATED(PDENA)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0) THEN
              BASE_DENSITY(IG) = PDENA(IOLD,IT)
              BASE_TEMP(IG) = EDENA(IOLD,IT)/(PDENA(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(2)
        IF (ASSOCIATED(PDENM)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0)  THEN
              BASE_DENSITY(IG) = PDENM(IOLD,IT)
              BASE_TEMP(IG) = EDENM(IOLD,IT)/(PDENM(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(3)
        IF (ASSOCIATED(PDENI)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0)  THEN
              BASE_DENSITY(IG) = PDENI(IOLD,IT)
              BASE_TEMP(IG) = EDENI(IOLD,IT)/(PDENI(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(4)
        BASE_DENSITY(1:NRAD) = DIIN(IOLD,1:NRAD)
        BASE_TEMP(1:NRAD) = TIIN(IOLDTI,1:NRAD)
      CASE DEFAULT
        WRITE (IUNOUT,*) 'BULK SPECIES ', IPLS,
     .                   ' DENSITY MODEL: ',CDENMODEL(IPLS)
        WRITE (IUNOUT,*) ' WRONG PARTICLE TYPE SPECIFIED '
        WRITE (IUNOUT,*) ' CHECK INPUT AND RERUN CASE '
        CALL EIRENE_EXIT_OWN(1)
      END SELECT
 
      END SUBROUTINE EIRENE_GET_BASE_DENSITY
 
 
 
      SUBROUTINE EIRENE_GET_SPECTRUM (ICELL,IRE,SPEC,FOUND)
 
      INTEGER, INTENT(IN) :: ICELL, IRE
      TYPE(EIRENE_SPECTRUM), INTENT(OUT) :: SPEC
      LOGICAL, INTENT(OUT) :: FOUND
      INTEGER :: ISPC
 
      FOUND = .FALSE.
 
      DO ISPC = 1, NADSPC
        IF ((ESTIML(ISPC)%PSPC%ISRFCLL == 2) .AND.
     .      (ESTIML(ISPC)%PSPC%ISPCSRF == ICELL) .AND.
     .      (ESTIML(ISPC)%PSPC%IPRTYP == TDMPAR(IPLS)%TDM%ITP(IRE)).AND.
     .      (ESTIML(ISPC)%PSPC%IPRSP == TDMPAR(IPLS)%TDM%ISP(IRE))) THEN
          ALLOCATE (SPEC%SPC(0:ESTIML(ISPC)%PSPC%NSPC+1))
          ALLOCATE (SPEC%SDV(0:ESTIML(ISPC)%PSPC%NSPC+1))
          ALLOCATE (SPEC%SGM(0:ESTIML(ISPC)%PSPC%NSPC+1))
          SPEC = ESTIML(ISPC)%PSPC
          FOUND = .TRUE.
        END IF
      END DO
 
      END SUBROUTINE EIRENE_GET_SPECTRUM
 
 
      END SUBROUTINE EIRENE_PLASMA_DERIV
 
