!pb  121206  check if tallies are available before using them

      SUBROUTINE SAVE_TALLIES (ISTRAI)
C
C  SAVE EIRENE TALLIES, SCALE PER UNIT FLUX (AMP), ON COMMON BRASCL
C  WTOTP IS NEGATIVE IN EIRENE (SINK FOR IONS)
C  ALL STRATA WHICH ARE NOT SPECIFIED BY INPUT BLOCK 14 (FROM
C  PLASMA CODE DATA) ARE NOT RESCALED HERE
C

      USE PRECISION
      USE PARMMOD
      USE BRASPOI
      USE CCOUPL
      USE COUTAU
      USE COMUSR
      USE CGRID
      USE CESTIM

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRAI
      REAL(DP) :: FLXI
      INTEGER :: IN, IATM, IMOL, IPLS, IION, ICPV

      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL

      IF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).NE.0.) THEN
         FLXI=-1._DP/WTOTP(0,ISTRAI)
      ELSEIF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).EQ.0.) THEN
         RETURN
      ELSEIF (ISTRAI.GT.NTARGI) THEN
         FLXI=1._DP
      ENDIF

      CALL FREE_SIMARR(ISTRAI)
      CALL FREE_MULARR(ISTRAI)

      DO IPLS=1,NPLSI
        DO IN=1,NSBOX_TAL
          IF (LPAPL) THEN 
	  IF (PAPL(IPLS,IN) .NE. 0.D0) THEN
!pb            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = PAPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => PAPLS(ISTRAI)%PMUL
            PAPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LPMPL) THEN 
          IF (PMPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = PMPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => PMPLS(ISTRAI)%PMUL
            PMPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LPIPL) THEN 
          IF (PIPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = PIPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => PIPLS(ISTRAI)%PMUL
            PIPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMAPL) THEN 
          IF (MAPL(IPLS,IN) .NE. 0.D0) THEN
!pb            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MAPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MAPLS(ISTRAI)%PMUL
            MAPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMMPL) THEN 
          IF (MMPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MMPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MMPLS(ISTRAI)%PMUL
            MMPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMIPL) THEN 
          IF (MIPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MIPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MIPLS(ISTRAI)%PMUL
            MIPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMPHPL) THEN 
          IF (MPHPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MPHPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MPHPLS(ISTRAI)%PMUL
            MPHPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO IN=1,NSBOX_TAL
	IF (LEAEL) THEN
        IF (EAEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EAEL(IN)*FLXI
          CPSIM%NXTSIM => EAELS(ISTRAI)%PSIM
          EAELS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEMEL) THEN 
        IF (EMEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EMEL(IN)*FLXI
          CPSIM%NXTSIM => EMELS(ISTRAI)%PSIM
          EMELS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEIEL) THEN 
        IF (EIEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EIEL(IN)*FLXI
          CPSIM%NXTSIM => EIELS(ISTRAI)%PSIM
          EIELS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEAPL) THEN 
        IF (EAPL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EAPL(IN)*FLXI
          CPSIM%NXTSIM => EAPLS(ISTRAI)%PSIM
          EAPLS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEMPL) THEN 
        IF (EMPL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EMPL(IN)*FLXI
          CPSIM%NXTSIM => EMPLS(ISTRAI)%PSIM
          EMPLS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEIPL) THEN 
        IF (EIPL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EIPL(IN)*FLXI
          CPSIM%NXTSIM => EIPLS(ISTRAI)%PSIM
          EIPLS(ISTRAI)%PSIM => CPSIM
        ENDIF
	ENDIF
      ENDDO

      DO IATM=1,NATMI
        DO IN=1,NSBOX_TAL
	  IF (LPDENA) THEN
          IF (PDENA(IATM,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IATM
            CPMUL%ICM = IN
            CPMUL%VALUEM = PDENA(IATM,IN)*FLXI
            CPMUL%NXTMUL => PDENAS(ISTRAI)%PMUL
            PDENAS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LEDENA) THEN 
          IF (EDENA(IATM,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IATM
            CPMUL%ICM = IN
            CPMUL%VALUEM = EDENA(IATM,IN)*FLXI
            CPMUL%NXTMUL => EDENAS(ISTRAI)%PMUL
            EDENAS(ISTRAI)%PMUL => CPMUL
          ENDIF
	  ENDIF
        ENDDO
      ENDDO

      DO IMOL=1,NMOLI
        DO IN=1,NSBOX_TAL
	  IF (LPDENM) THEN
          IF (PDENM(IMOL,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IMOL
            CPMUL%ICM = IN
            CPMUL%VALUEM = PDENM(IMOL,IN)*FLXI
            CPMUL%NXTMUL => PDENMS(ISTRAI)%PMUL
            PDENMS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO IION=1,NIONI
        DO IN=1,NSBOX_TAL
	  IF (LPDENM) THEN
          IF (PDENI(IION,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = IION
            CPMUL%ICM = IN
            CPMUL%VALUEM = PDENI(IION,IN)*FLXI
            CPMUL%NXTMUL => PDENIS(ISTRAI)%PMUL
            PDENIS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO ICPV=1,NCPVI
        DO IN=1,NSBOX_TAL
	  IF (LCOPV) THEN
          IF (COPV(ICPV,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => NEW_MULARR()
            CPMUL%IART = ICPV
            CPMUL%ICM = IN
            CPMUL%VALUEM = COPV(ICPV,IN)*FLXI
            CPMUL%NXTMUL => COPVS(ISTRAI)%PMUL
            COPVS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

