C
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
          DO 122 IPLS=1,NPLSI
C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IATM AND IPLS
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 122
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS).AND.NCHRGP(IPLS).EQ.1)
     .          GOTO 125
124           CONTINUE
              GOTO 122
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
              NRCX=NRCX+1
C
            ENDIF
122       CONTINUE
C
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
csw  PHOTON COLLISION, OT - type
csw
      nnrot=0
      do iatm=1,natmi
         if(nrca(iatm) > 0) then
            do nrc=1,nrca(iatm)
               kk=ireaca(iatm,nrc)
               if(iswr(kk) == 7) then
                  nNROT=nNROT+1
               endif
            enddo
         endif
      enddo
      call PH_ALLOC_XSECTA(nnrot)
csw
C
      CALL XSTAPI_PARAM
C
      RETURN
C
      END





