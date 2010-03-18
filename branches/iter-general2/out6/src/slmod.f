c
c ======================================================================
c
      SUBROUTINE LDADAS2(CZ,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                   TENUM,TEVALS,NENUM,NEVALS,MODE,
     .                   CVALS,WAVE,IRCODE)
c
c      SUBROUTINE LDADAS(CZ,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
c     >                  CVALS,WAVE,IRCODE)
c
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDADAS:  CODE TO LOAD THE REQUESTED LINE EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE LINE IS SPECIFIED BY THE    *
C  *           NUCLEAR CHARGE OF THE EMITTING ION, THE                 *
C  *           CHARGE OF THE EMITTING IONISATION STATE, AN ID FLAG     *
C  *           WHICH LOCATES THE INPUT FILE, AND THREE BLOCK SELECTOR  *
C  *           NUMBERS, ONE EACH FOR EMISSION BY ELECTRON EXCITATION,  *
C  *           RECOMBINATION FROM THE NEXT HIGHER IONISATION STATE,    *
C  *           AND CHARGE EXCHANGE FROM THE NEXT HIGHER IONISATION     *
C  *           STATE.  IN ADDITION TO THE EMISSION PROFILE, THE        *
C  *           ROUTINE RETURNS THE WAVELENGTH OF THE TRANSITION AND    *
C  *           AN ERROR CODE FROM THE ADAS EXTRACTION ROUTINE, SPEC.   *
C  *           NOTE THAT THERE IS NO CHECKING OF THE SELECTOR NUMBERS! *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *                                                                   *
C  *********************************************************************
C
      include 'params'
      INCLUDE 'cgeom'
      include 'pindata'
      include 'dynam2'
      include 'comtor'
C
      REAL      tevals(MAXGXS,MAXNGS),nevals(MAXGXS,MAXNGS)
      INTEGER   tenum,nenum,mode

      INTEGER   CZ,IZ,ADASYR,ISELE,ISELR,ISELX,IRCODE
      REAL      WAVE,CVALS(MAXGXS,MAXNGS)
c
c      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
c
      CHARACTER ADASID*(*),ADASEX*(*)
C
      INTEGER   IR,IK,IADAS,NPAIRS,IKK
      REAL*8    TADAS(20),DADAS(20)
      REAL*8    WLNGTH,PECAE(20),PECAR(20),PECAX(20)
      LOGICAL*4 LTRNG(20),LDRNG(20)
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2
C
      WAVE = 0.0
      IRCODE = 0
      CALL XXUID(ADASID)
      IF (ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
      CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
      DO IR = 1,NENUM
        DO IK = 1,TENUM,20
          NPAIRS = MIN0(20,TENUM-(IK-1))
c
c      DO IR = 1,NRS
c        DO IK = 1,NKS(IR),20
c
c          NPAIRS = MIN0(20,NKS(IR)-(IK-1))
c
          DO IADAS = 1,NPAIRS
            TADAS(IADAS) = DBLE(TEVALS(IK+(IADAS-1),IR))
            DADAS(IADAS) = DBLE(1.E-6*RIZB*NEVALS(IK+(IADAS-1),IR))
c
c            TADAS(IADAS) = DBLE(KTEBS(IK+(IADAS-1),IR))
c            DADAS(IADAS) = DBLE(1.E-6*RIZB*KNBS(IK+(IADAS-1),IR))
c
          ENDDO
C
          CALL DZERO(PECAE,NPAIRS)
          IF (ISELE.GT.0) THEN
            CALL SPEC(ISELE,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELE.EQ.-1) THEN
C
C  JUST LOAD EMISSION MEASURE FOR ISEL = -1
C    - SINCE THIS VALUE CAN EXCEED THE UNIX SINGLE PRECISION LIMIT,
C      WORK IN DENSITY UNITS OF 10**18
C
            CALL DINIT(PECAE,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          CALL DZERO(PECAR,NPAIRS)
          IF (ISELR.GT.0) THEN
            CALL SPEC(ISELR,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELR.EQ.-1) THEN
            CALL DINIT(PECAR,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
C  FOR IMPURITIES USE THE THIRD SWITCH FOR CX, FOR HYDROGEN
C  ADD IN A MOLECULAR CONTRIBUTION INSTEAD
C
          CALL DZERO(PECAX,NPAIRS)
CLDH - USE ION TEMPERATURE FOR CX RATE
          IF (CZ.GT.1.0) THEN
            DO IADAS = 1,NPAIRS
              TADAS(IADAS) = DBLE(KTIBS(IK+(IADAS-1),IR))
            ENDDO
          ENDIF
CLDH
          IF (ISELX.GT.0) THEN
            CALL SPEC(ISELX,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELX.EQ.-1) THEN
            CALL DINIT(PECAX,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          DO IADAS = 1,NPAIRS
            IKK = IK + (IADAS-1)
            IF (CZ.GT.1.0) THEN
              CVALS(IKK,IR) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ)
     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ+1)
     >        +PECAX(IADAS) * PINATOM(IKK,IR) * SDLIMS(IKK,IR,IZ+1))
            ELSE
C
C---- HYDROGEN DENSITIES ARE IN DIFFERENT ARRAYS
C
              IF (mode.EQ.1) THEN
                CVALS(IKK,IR) = PECAE(IADAS) * 1.0D-6
              ELSEIF (mode.EQ.2) THEN
                CVALS(IKK,IR) = PECAR(IADAS) * 1.0D-6
              ELSEIF (mode.EQ.3) THEN
                CVALS(IKK,IR) = PECAX(IADAS) * 1.0D-6
              ENDIF

c
c              CVALS(IKK,IR) = 1.D-6*
c     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR) * PINATOM(IKK,IR)
c     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR) * KNBS(IKK,IR)
c     >        +PECAX(IADAS) * KNBS(IKK,IR)   * PINMOL(IKK,IR))
c
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      WAVE = WLNGTH
C
      RETURN
      END
