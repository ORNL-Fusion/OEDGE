! 23.08.06: VPX, VPY, VRX, VRY changed to ALLOCATABLE, SAVE to speed up
!           subroutine call (save time in storage allocation)
C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CZT1
      USE EIRMOD_CTRIG
      USE EIRMOD_CGEOM
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNP(:)
      REAL(DP), ALLOCATABLE :: PPPL_COP(:,:), CPPV(:,:),
     .                         EPPL_COP(:), EPEL(:)
      REAL(DP) :: RECTOT, SUMN, SUMM, SUMEI, SUMEE, RECADD, EEADD,
     .            EIRENE_FTABRC1, EIRENE_FEELRC1, PIADD, EIADD
      INTEGER :: IPL, JPLS, ISR, ISTEP, IIRC, IU, IPLSTI, IRRC, IN, INC
      INTEGER, SAVE :: IFIRST, ISTROLD=-1
      DATA IFIRST/0/
 
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNP(NPLS))
        DO IPL=1,NPLSI
          CNDYNP(IPL)=AMUA*RMASSP(IPL)
        END DO
      ENDIF

      IF (ISTRA /= ISTROLD) THEN
C
C
C  ADD CONTRIBUTIONS FROM VOLUME RECOMBINATION SOURCE
C
        ALLOCATE (PPPL_COP(NPLS,NRAD))
        ALLOCATE (CPPV(NPLS,NRAD))
        ALLOCATE (EPPL_COP(NRAD))
        ALLOCATE (EPEL(NRAD))
        PPPL_COP =0.D0
        CPPV = 0.D0
        EPPL_COP = 0.D0
        EPEL = 0.D0

        IF (NLVOL(ISTRA)) THEN
C
          RECTOT = 0._DP
          DO 7473 IPLS=1,NPLSI
            JPLS = NSPEZ(ISTRA)
            IF ((JPLS > 0) .AND. (JPLS <= NPLSI) .AND.
     .          (IPLS /= JPLS)) CYCLE
            IPLSTI= MPLSTI(IPLS)
            DO ISR=1, NSRFSI(ISTRA)
              ISTEP = SORIND(ISR,ISTRA) 
              DO 7472 IIRC=1,NPRCI(IPLS)
                IRRC=LGPRC(IPLS,IIRC)
                IF ((ISTEP > 0) .AND. (ISTEP /= IRRC)) CYCLE
                SUMN=0.0
                SUMM=0.0
                SUMEI=0.0
                SUMEE=0.0
                DO IN=1,NTRII
                  INC=NCLTAL(IN)
                  IF (NSTORDR >= NRAD) THEN
                    RECADD=-TABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                    EEADD=  EELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  ELSE
                    RECADD=-EIRENE_FTABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                    EEADD=  EIRENE_FEELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  END IF
                  PPPL_COP(IPLS,INC)=PPPL_COP(IPLS,INC)+RECADD
                  SUMN=SUMN+RECADD*VOL(IN)
                  PIADD=PARMOM(IPLS,IN)*RECADD
                  CPPV(IPLS,INC)=CPPV(IPLS,INC)+PIADD
                  SUMM=SUMM+PIADD*VOL(IN)
                  EIADD=(1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))*RECADD
                  EPPL_COP(INC)=EPPL_COP(INC)+EIADD
                  SUMEI=SUMEI+EIADD*VOL(IN)
                  EPEL(INC)=EPEL(INC)+EEADD
                  SUMEE=SUMEE+EEADD*VOL(IN)
                END DO
                RECTOT = RECTOT + SUMN
                WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
                CALL EIRENE_MASR4('SUMN, SUMM, SUMEI, SUMEE        ',
     .                      SUMN,SUMM,SUMEI,SUMEE)
7472          CONTINUE
            END DO
7473      CONTINUE
        END IF

        IF (LCOPV) THEN
          IU = UBOUND(COPV,1)

          IF (IU >= NPLSI) 
     .      COPV(1:NPLSI,:) = PPPL_COP(1:NPLSI,:)
          IF (IU >= 2*NPLSI) 
     .      COPV(NPLSI+1:2*NPLSI,:) = CPPV(1:NPLSI,:)
          IF (IU >= 2*NPLSI+1) 
     .      COPV(2*NPLSI+1,:) = EPPL_COP(:)
          IF (IU >= 2*NPLSI+2) 
     .      COPV(2*NPLSI+2,:) = EPEL(:)
        END IF

        DEALLOCATE (PPPL_COP)
        DEALLOCATE (CPPV)
        DEALLOCATE (EPPL_COP)
        DEALLOCATE (EPEL)

        ISTROLD = ISTRA

      END IF

      RETURN
      END
 
 
 
