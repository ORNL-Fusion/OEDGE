C  may 2006:  default helium resonant cx added
C  april 07:  setting of NRPI adapted (ion impact collisions)
C
      SUBROUTINE EIRENE_XSECTA_PARAM
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_COMSOU
      USE EIRMOD_CTEXT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_PHOTON
 
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
          DO 155 IPLS=1,NPLSI
C CHECK: "ATOMIC" COLLISION PARTNERS ONLY
            IF (NPRT(NSPAMI+IPLS).NE.1.OR.NPRT(NSPH+IATM).NE.1) GOTO 155
C
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1.AND.
     .          NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
C
C  FIND BULK SECONDARIES
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 155
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 125
124           CONTINUE
              GOTO 155
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
              NRCX=NRCX+1
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE+ ION
C
C  FIND BULK SECONDARIES
              DO 131 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 133
131           CONTINUE
              GOTO 155
133           DO 134 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 135
134           CONTINUE
              GOTO 155
 
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
135           CONTINUE
              NRCX=NRCX+1
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.2) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE++ ION
C
C  FIND BULK SECONDARIES
              DO 141 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.2)
     .          GOTO 143
141           CONTINUE
              GOTO 155
143           DO 144 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 145
144           CONTINUE
              GOTO 155
 
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
145           CONTINUE
              NRCX=NRCX+1
 
            ENDIF
155       CONTINUE
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
csw  COLLISIONS OF ATOMS WITH PHOTON BACKGROUND, OT - type
csw
cdr  not active
c      IF (.FALSE.) THEN
c      nnrot=0
c      do iatm=1,natmi
c         if(nrca(iatm) > 0) then
c            do nrc=1,nrca(iatm)
c               kk=ireaca(iatm,nrc)
c               if(iswr(kk) == 7) then
c                  nNROT=nNROT+1
c               endif
c            enddo
c         endif
c      enddo
c      call PH_ALLOC_XSECTA(nnrot)
c      END IF
csw
C
C
C   ION IMPACT COLLISIONS
C
 
      DO 1000 IATM=1,NATMI
C
C  NO DEFAULT MODEL
C
       IF (NRCA(IATM).EQ.0) THEN
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 150 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.4) GOTO 150
            NRPI=NRPI+1
C
150       CONTINUE
C
C  NO MODEL DEFINED
        ELSE
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
      END
 
 
 
 
 
