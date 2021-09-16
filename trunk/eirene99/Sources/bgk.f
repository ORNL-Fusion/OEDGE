c  code segment: bgk
c
c  only needed, if some test species are labeled as bgk-species
c               with non-linear self interactions
c               this segment contains a routine which updates the tallies
c               required for iteration (UPTBGK),
C               and a routine (MODBGK) for doing the iterations.
c               the standard deviations for the "bgk-tallies" are
c               computed in subroutine STATIS_BGK
C
c
c
      SUBROUTINE UPTBGK(WV,NPBGK)
C  UPDATE BGK-SPECIFIC TALLIES, TRACKLENGTH ESTIMATORS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CESTIM'
      INCLUDE 'COMPRT'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'CUPD'
      INCLUDE 'CGRID'
      INCLUDE 'CTEXT'
      INCLUDE 'CSDVI'
      CHARACTER*8 TXT
      DATA IFIRST/0/
      SAVE
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        NSBGK=NRBGI/3
        DO IBGK=1,NSBGK
          ITP=0
          DO ISP=1,NATMI
            IF (NPBGKA(ISP).EQ.IBGK) THEN
              ITP=1
              IAT=ISP
              TXT=TEXTS(IAT)
              GOTO 1
            ENDIF
          ENDDO
          DO ISP=1,NMOLI
            IF (NPBGKM(ISP).EQ.IBGK) THEN
              ITP=2
              IML=ISP
              TXT=TEXTS(NSPA+IML)
              GOTO 1
            ENDIF
          ENDDO
          DO ISP=1,NIONI
c slmod begin - not tr (fixed in 02)
c...bug?
            IF (NPBGKI(ISP).EQ.IBGK) THEN
c
c            IF (NPBGKI(ISP).EQ.IBG) THEN
c slmod end
              ITP=3
              IIO=ISP
              TXT=TEXTS(NSPAM+IIO)
              GOTO 1
            ENDIF
          ENDDO
          WRITE (6,*) 'SPECIES ERROR IN UPTBGK'
          CALL EXIT
1         CONTINUE
C
C  BGK-SPECIES NO. IBGK
          IUPD1=(IBGK-1)*3+1
          IUPD2=(IBGK-1)*3+2
          IUPD3=(IBGK-1)*3+3
          TXTTAL(IUPD1,NTALB)='BGK TALLY: FLUX DENSITY IN X DIRECTION '
          TXTTAL(IUPD2,NTALB)='BGK TALLY: FLUX DENSITY IN Y DIRECTION '
          TXTTAL(IUPD3,NTALB)='BGK TALLY: FLUX DENSITY IN Z DIRECTION '
          TXTUNT(IUPD1,NTALB)='#/CM**2/S               '
          TXTUNT(IUPD2,NTALB)='#/CM**2/S               '
          TXTUNT(IUPD3,NTALB)='#/CM**2/S               '
          TXTSPC(IUPD1,NTALB)=TXT
          TXTSPC(IUPD2,NTALB)=TXT
          TXTSPC(IUPD3,NTALB)=TXT
          IBGVE(IUPD1)=1
          IBGVE(IUPD2)=1
          IBGVE(IUPD3)=1
          IBGRC(IUPD1)=ITP
          IBGRC(IUPD2)=ITP
          IBGRC(IUPD3)=ITP
        ENDDO

        NMTSP=NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI+NCPVI
C
C  END OF IFIRST BLOCK
      ENDIF
C
C  UPDATE BGK TALLIES
C  PRESENTLY: UPDATE TRANSPORT FLUX VECTOR ON BGKV-TALLY
C
      IBGK=NPBGK
      IUPD1=(IBGK-1)*3+1
      IUPD2=(IBGK-1)*3+2
      IUPD3=(IBGK-1)*3+3
      LMETSP(NMTSP+IUPD1)=.TRUE.
      LMETSP(NMTSP+IUPD2)=.TRUE.
      LMETSP(NMTSP+IUPD3)=.TRUE.
      DO 51 I=1,NCOU
        DIST=CLPD(I)
        WTRV=WV*DIST*VEL
        WTRVX=WTRV*VELX
        WTRVY=WTRV*VELY
        WTRVZ=WTRV*VELZ
        IRD=NRCELL+NUPC(I)*NR1P2+NBLCKA
        BGKV(IUPD1,IRD)=BGKV(IUPD1,IRD)+WTRVX
        BGKV(IUPD2,IRD)=BGKV(IUPD2,IRD)+WTRVY
        BGKV(IUPD3,IRD)=BGKV(IUPD3,IRD)+WTRVZ
51    CONTINUE

      RETURN
      END
C
      SUBROUTINE MODBGK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'COMPRT'
      INCLUDE 'COMNNL'
      INCLUDE 'COMSOU'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CINIT'
      INCLUDE 'COMXS'
      INCLUDE 'CSPEI'
      INCLUDE 'CTRCEI'
      INCLUDE 'CCONA'
      INCLUDE 'CESTIM'
      INCLUDE 'CSDVI'
      INCLUDE 'CSDVI_BGK'
      INCLUDE 'CSDVI_COP'
      INCLUDE 'COUTAU'
      INCLUDE 'CLOGAU'
      INCLUDE 'CZT1'
C
CPB   DIMENSION PDEN(NRAD),EDEN(NRAD),PDEN2(NRAD),EDEN2(NRAD),
CPB  .          ENERGY(NPLS,NRAD)
      DOUBLE PRECISION, ALLOCATABLE :: PDEN(:), EDEN(:),
     .          PDEN2(:), EDEN2(:), ENERGY(:,:)
      DIMENSION ITYP1(NPLS),ITYP2(NPLS),ISPZ1(NPLS),ISPZ2(NPLS),
     .          IREL1(NPLS),INRC1(NPLS),RATM(3)
      LOGICAL LMARK(NPLS)
      LOGICAL TRCSAV
C
      DIMENSION OUTAU(NOUTAU),ESTIM(NESTIM),SDVI(NSDVI),SDVI_BGK(NSBGK),
     .          SDVI_COP(NSCOP)
      EQUIVALENCE
     .  (OUTAU(1),PDENAI(0,0)),
     .  (ESTIM(1),PDENA(1,1)),
     .  (SDVI(1),SIGMA(1,1)),
     .  (SDVI_BGK(1),SIGMA_BGK(1,1)),
     .  (SDVI_COP(1),SIGMA_COP(1,1))
      DIMENSION RCMDTA(NMDTA),ICMDTA(MMDTA)
      EQUIVALENCE (RCMDTA(1),TABDS1(1,1)),(ICMDTA(1),MODCOL(1,1,1,0))
c slmod begin - tr (modified)
      COMMON /RESCOM/ EIRRES
      REAL*8          EIRRES(7,NPLS)

c      COMMON /RESCOM/ RRN,RRE,RRM,RATN,RATE,RATM,RESN,RESE,RESM

      DIMENSION ASDSAV(NPLS,NRAD,5)
      DOUBLE PRECISION SLTIME

      CALL DSET(ASDSAV,NPLS*NRAD*5,0.0D0)

c      CALL Clock_1970(SLTIME)

      IF (output)
     .WRITE (0 ,*) 'WRITING EXCHANGE RATE DATA'
c
c      jdemod - since assignment to SLTIME is commented out above
c               also comment out this write statement
c
c      WRITE (97,*) 'PARTICLE, MOMENTUM AND ENERGY EXCHANGE RATES ',
c     .             SLTIME
c
      WRITE (97,*)
c slmod end

      CALL LEER(3)
      WRITE (6,*) 'MODBGK CALLED AFTER ITERATION IITER= ',IITER
      CALL LEER(3)
C
      A=1.D0/DFLOAT(IITER)
      ISTRA=0
      IF (NSTRAI.EQ.1.AND.IESTR.EQ.1) ISTRA=1
      IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.NE.0) THEN
        IESTR=0
        CALL RSTRT(0,NSTRAI,NESTIM,NSDVI,ESTIM,SDVI,
     .             NSBGK,SDVI_BGK,NSCOP,SDVI_COP,TRCFLE)
      ELSE
        WRITE (6,*) 'ERROR IN MODBGK: DATA FOR STRATUM ISTRA= ',ISTRA
        WRITE (6,*) 'ARE NOT AVAILABLE. EXIT CALLED'
        CALL EXIT
      ENDIF
C
      if (nbmlt.gt.1) then
        write (6,*) 'option not ready in modbgk '
        write (0,*) 'MARK: *** option not ready in modbgk, but not '//
     .              'exiting ***'
c        call exit
      endif
C
      NBLCKA=0
      ALLOCATE (PDEN(NRAD))
      ALLOCATE (EDEN(NRAD))
      ALLOCATE (PDEN2(NRAD))
      ALLOCATE (EDEN2(NRAD))
      ALLOCATE (ENERGY(NPLS,NRAD))
c
C  LOOP OVER THOSE BACKGROUND ION SPECIES, WHICH ARE ARTIFICIAL
C  SPECIES FOR (NON-LINEAR) ITERATIONS
C .........................................................................
      DO 1000 IPLS=1,NPLSI
C .........................................................................
C
C  IS IPLS AN ARTIFICIAL "BGK BACKGROUND SPECIES"?
C
        ITYP1(IPLS)=0
        ISPZ1(IPLS)=0
C
        IF (NPBGKP(IPLS,1).EQ.0) GOTO 1000
c slmod begin - tr
        WRITE (97,*) 'IPLS= ',IPLS
c slmod end
C
C  YES. FIND INCIDENT TEST PARTICLE COLLISION PARTNER: ITYP, ISPZ, IREL

        IBGK1=NPBGKP(IPLS,1)
        IUP1=(IBGK1-1)*3+1
        IUP2=(IBGK1-1)*3+2
        IUP3=(IBGK1-1)*3+3
C
C  TRY ATOMS
        DO IATM=1,NATMI
          IF (NPBGKA(IATM).EQ.IBGK1) THEN
            ITYP1(IPLS)=1
            ISPZ1(IPLS)=IATM
            FACT1=CVRSSA(IATM)
            RMAS1=RMASSA(IATM)
            DO IRAD=1,NRAD
              PDEN(IRAD)=PDENA(IATM,IRAD)
              EDEN(IRAD)=EDENA(IATM,IRAD)
            ENDDO
C  FIND INDEX  NRC
            DO NRC=1,NRCA(IATM)
              IP=IDEZ(IBULKA(IATM,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 10
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
10          DO IAEL=1,NAELI(IATM)
              IF  (LGAEL(IATM,IAEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGAEL(IATM,IAEL,0)
                GOTO 100
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY MOLECULES
        DO IMOL=1,NMOLI
          IF (NPBGKM(IMOL).EQ.IBGK1) THEN
            ITYP1(IPLS)=2
            ISPZ1(IPLS)=IMOL
            FACT1=CVRSSM(IMOL)
            RMAS1=RMASSM(IMOL)
            DO IRAD=1,NRAD
              PDEN(IRAD)=PDENM(IMOL,IRAD)
              EDEN(IRAD)=EDENM(IMOL,IRAD)
            ENDDO
C  FIND INDEX  NRC
            DO NRC=1,NRCM(IMOL)
              IP=IDEZ(IBULKM(IMOL,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 20
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
20          DO IMEL=1,NMELI(IMOL)
              IF  (LGMEL(IMOL,IMEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGMEL(IMOL,IMEL,0)
                GOTO 100
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY TEST IONS
        DO IION=1,NIONI
          IF (NPBGKI(IION).EQ.IBGK1) THEN
            ITYP1(IPLS)=3
            ISPZ1(IPLS)=IION
            FACT1=CVRSSI(IION)
            RMAS1=RMASSI(IION)
            DO IRAD=1,NRAD
              PDEN(IRAD)=PDENI(IION,IRAD)
              EDEN(IRAD)=EDENI(IION,IRAD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCI(IION)
              IP=IDEZ(IBULKI(IION,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 30
              ENDIF
            ENDDO
C  FIND INDEX IREL
30          DO IIEL=1,NIELI(IION)
              IF  (LGIEL(IION,IIEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGIEL(IION,IIEL,0)
                GOTO 100
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
        GOTO 995
C
C  SELF COLLISION OR CROSS COLLISION
C
100     CONTINUE
C
        IF (NPBGKP(IPLS,2).EQ.0) THEN
          ITYP2(IPLS)=0
          ISPZ2(IPLS)=0
C
        ELSEIF (NPBGKP(IPLS,2).NE.0) THEN
C
C  CROSS COLLISION, FIND SECOND COLLISION PARTNER
C  THIS IS NOT THE INGOING COLLIDING TESTPARTICLE,
C  (DETERMINING, E.G., MASS AND DENSITY OF BACKGROUND PARTICLE)
C  BUT THE SECOND, 'ARTIFICIAL', PARTNER
C
          ITYP2(IPLS)=IDEZ(NPBGKP(IPLS,2),1,3)
          ISPZ2(IPLS)=IDEZ(NPBGKP(IPLS,2),3,3)
C
          IF (ITYP2(IPLS).EQ.1) THEN
            IATM2=ISPZ2(IPLS)
            FACT2=CVRSSA(IATM2)
            RMAS2=RMASSA(IATM2)
            DO IRAD=1,NRAD
              PDEN2(IRAD)=PDENA(IATM2,IRAD)
              EDEN2(IRAD)=EDENA(IATM2,IRAD)
            ENDDO
            IBGK2=NPBGKA(IATM2)
          ELSEIF (ITYP2(IPLS).EQ.2) THEN
            IMOL2=ISPZ2(IPLS)
            FACT2=CVRSSM(IMOL2)
            RMAS2=RMASSM(IMOL2)
            DO IRAD=1,NRAD
              PDEN2(IRAD)=PDENM(IMOL2,IRAD)
              EDEN2(IRAD)=EDENM(IMOL2,IRAD)
            ENDDO
            IBGK2=NPBGKM(IMOL2)
          ELSEIF (ITYP2(IPLS).EQ.3) THEN
            IION2=ISPZ2(IPLS)
            FACT2=CVRSSI(IION2)
            RMAS2=RMASSI(IION2)
            DO IRAD=1,NRAD
              PDEN2(IRAD)=PDENI(IION2,IRAD)
              EDEN2(IRAD)=EDENI(IION2,IRAD)
            ENDDO
            IBGK2=NPBGKI(IION2)
          ENDIF
          IUP12=(IBGK2-1)*3+1
          IUP22=(IBGK2-1)*3+2
          IUP32=(IBGK2-1)*3+3
C
        ENDIF
C
        IF (RMASSP(IPLS).NE.RMAS1) THEN
          WRITE (6,*) 'INCONSISTENT MASS FOR IPLS= ',IPLS
          WRITE (6,*) 'RMASSP(IPLS),RMAS1 ',RMASSP(IPLS),RMAS1
          CALL EXIT
        ENDIF
C
        CNDYN=AMUA*RMAS1
C
        NXM=MAX(1,NR1STM)
        NYM=MAX(1,NP2NDM)
        NZM=MAX(1,NT3RDM)
c slmod begin - tr
        IF (NLMLT) THEN
          NZM=NBMLT
        ENDIF
c slmod end
C
        RRN=0.
        RRE=0.
        RRM=0.
        RATN=0.
        RATE=0.
        RATM=0.
        RESN=0.
        RESE=0.
        RESM=0.
C
        IREL=IREL1(IPLS)
C
        IF (ITYP2(IPLS).EQ.0) THEN
C
        WRITE (6,*) 'MODBGK: SELF COLLISION, IPLS',IPLS
        WRITE (6,*) 'ITYP,ISPZ,IBGK,IREL ',ITYP1(IPLS),ISPZ1(IPLS),
     .                                     IBGK1,IREL1(IPLS)

            DO 80 IT=1,NZM

        DO 80 IR=1,NXM
          DO 80 IP=1,NYM

c slmod begin - tr
              IF (NLMLT) THEN
                IRAD=IR+(IP-1)*NR1P2+(IT-1)*NSTRD

c                WRITE(0,*) 'MARK: IRAD=',irad,
c     .                      IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA,
c     .                      NSTRD,NSURF
              ELSE
                IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
              ENDIF
c
c              IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
c slmod end
C
              TBEL=0.
              IF (LGVAC(IRAD,IPLS)) GOTO 81
              IF (NSTORDR >= NRAD) THEN
                TBEL = TABEL3(IREL,IRAD,1)
              ELSE
                KK=NREAEL(IREL)
                PLS=TIINL(IPLS,IRAD)+ADDEL(IREL,IPLS)
                TBEL = CREAC(9,1,KK)
                DO II=8,1,-1
                  TBEL = TBEL*PLS + CREAC(II,1,KK)
                END DO
                TBEL=EXP(MAX(-100.D0,TBEL))*DIIN(IPLS,IRAD)
              END IF
81            CONTINUE
C DELTA_N
              DOLD=DIIN(IPLS,IRAD)
              DEL=DOLD-PDEN(IRAD)
              RATN=RATN+TBEL*DEL*VOL(IRAD)
c              IF (TRCMOD) THEN
c                WRITE (6,*) 'IR,T,TBEL,VOL,RATN ',IRAD,TIIN(IPLS,IRAD),
c     .                               TBEL,VOL(IRAD),TBEL*DEL*VOL(IRAD)
c                WRITE (6,*) 'DOLD,DNEW,DEL      ',DOLD,PDEN(IRAD),DEL
c              ENDIF
              RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
C DELTA_E
c             EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
c    .             DIIN(IPLS,IRAD)
              EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
     .             PDEN(IRAD)
              DEL=EOLD-EDEN(IRAD)
              RATE=RATE+TBEL*DEL*VOL(IRAD)
              RESE=RESE+TBEL*ABS(DEL)*VOL(IRAD)
C DELTA_V
              DELX=BGKV(IUP1,IRAD)-VXIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
              DELY=BGKV(IUP2,IRAD)-VYIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
              DELZ=BGKV(IUP3,IRAD)-VZIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
              RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
              RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
              RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C NEW T, NEW V
              VX=BGKV(IUP1,IRAD)/(PDEN(IRAD)+EPS60)
              VY=BGKV(IUP2,IRAD)/(PDEN(IRAD)+EPS60)
              VZ=BGKV(IUP3,IRAD)/(PDEN(IRAD)+EPS60)
              VXIN(IPLS,IRAD)=VX
              VYIN(IPLS,IRAD)=VY
              VZIN(IPLS,IRAD)=VZ
              ED=(VX**2+VY**2+VZ**2)*FACT1
              TIIN(IPLS,IRAD)=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED)/1.5
C NEW N
              DIIN(IPLS,IRAD)=PDEN(IRAD)
              RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C             RRM=?
              RRE=RRE+EDEN(IRAD)*VOL(IRAD)
80      CONTINUE
c
c  same as do 80 loop , for additional cell region
c
        DO 90 IRAD=NSURF+1,NSURF+NRADD
C
          TBEL=0.
          IF (LGVAC(IRAD,IPLS)) GOTO 91
          IF (NSTORDR >= NRAD) THEN
            TBEL = TABEL3(IREL,IRAD,1)
          ELSE
            KK=NREAEL(IREL)
            PLS=TIINL(IPLS,IRAD)+ADDEL(IREL,IPLS)
            TBEL = CREAC(9,1,KK)
            DO II=8,1,-1
              TBEL = TBEL*PLS + CREAC(II,1,KK)
            END DO
            TBEL=EXP(MAX(-100.D0,TBEL))*DIIN(IPLS,IRAD)
          END IF
91        CONTINUE
C DELTA_N
          DOLD=DIIN(IPLS,IRAD)
          DEL=DOLD-PDEN(IRAD)
          RATN=RATN+TBEL*DEL*VOL(IRAD)
c          IF (TRCMOD) THEN
c            WRITE (6,*) 'IR,T,TBEL,VOL,RATN ',IRAD,TIIN(IPLS,IRAD),
c     .                           TBEL,VOL(IRAD),TBEL*DEL*VOL(IRAD)
c            WRITE (6,*) 'DOLD,DNEW,DEL      ',DOLD,PDEN(IRAD),DEL
c          ENDIF
          RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
C DELTA_E
c         EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
c    .          DIIN(IPLS,IRAD)
          EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
     .          PDEN(IRAD)
          DEL=EOLD-EDEN(IRAD)
          RATE=RATE+TBEL*DEL*VOL(IRAD)
          RESE=RESE+TBEL*ABS(DEL)*VOL(IRAD)
C DELTA_V
          DELX=BGKV(IUP1,IRAD)-VXIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
          DELY=BGKV(IUP2,IRAD)-VYIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
          DELZ=BGKV(IUP3,IRAD)-VZIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
          RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
          RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
          RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C NEW T, NEW V
          VX=BGKV(IUP1,IRAD)/(PDEN(IRAD)+EPS60)
          VY=BGKV(IUP2,IRAD)/(PDEN(IRAD)+EPS60)
          VZ=BGKV(IUP3,IRAD)/(PDEN(IRAD)+EPS60)
          VXIN(IPLS,IRAD)=VX
          VYIN(IPLS,IRAD)=VY
          VZIN(IPLS,IRAD)=VZ
          ED=(VX**2+VY**2+VZ**2)*FACT1
c slmod begin - tr
c...temp: Ignore 1st additional cell (leave as vacuum)
          IF (IRAD.EQ.NSURF+1) CYCLE
c slmod end
          TIIN(IPLS,IRAD)=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED)/1.5
C NEW N
          DIIN(IPLS,IRAD)=PDEN(IRAD)
          RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C         RRM=?
          RRE=RRE+EDEN(IRAD)*VOL(IRAD)
90      CONTINUE
C
        ELSE
C
        WRITE (6,*) 'MODBGK: CROSS COLLISION, IPLS ',IPLS
        WRITE (6,*) 'ITYP1,ISPZ1,IBGK1,IREL1 ',ITYP1(IPLS),ISPZ1(IPLS),
     .                                         IBGK1,IREL1(IPLS)
        WRITE (6,*) 'ITYP2,ISPZ2,IBGK2       ',ITYP2(IPLS),ISPZ2(IPLS),
     .                                         IBGK2
        DO 180 IR=1,NXM
          DO 180 IP=1,NYM
            DO 180 IT=1,NZM
c slmod begin - tr
              IF (NLMLT) THEN
                IRAD=IR + (IP-1)*NR1P2 + (IT-1)*NSTRD
              ELSE
                IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
              ENDIF
c
c              IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
c slmod end
C
              TBEL=0
              IF (LGVAC(IRAD,IPLS)) GOTO 181
              IF (NSTORDR >= NRAD) THEN
                TBEL = TABEL3(IREL,IRAD,1)
              ELSE
                KK=NREAEL(IREL)
                PLS=TIINL(IPLS,IRAD)+ADDEL(IREL,IPLS)
                TBEL = CREAC(9,1,KK)
                DO II=8,1,-1
                  TBEL = TBEL*PLS + CREAC(II,1,KK)
                END DO
                TBEL=EXP(MAX(-100.D0,TBEL))*DIIN(IPLS,IRAD)
              END IF
181           CONTINUE
c             EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
c    .              DIIN(IPLS,IRAD)
              EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
     .              PDEN(IRAD)
              DOLD=DIIN(IPLS,IRAD)
              DEL=EOLD-EDEN(IRAD)
              RATE=RATE+TBEL*DEL*VOL(IRAD)
              DELX=BGKV(IUP1,IRAD)-VXIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
              DELY=BGKV(IUP2,IRAD)-VYIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
              DELZ=BGKV(IUP3,IRAD)-VZIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
              RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
              RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
              RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C
              VXIN1=BGKV(IUP1 ,IRAD)/(PDEN (IRAD)+EPS60)
              VXIN2=BGKV(IUP12,IRAD)/(PDEN2(IRAD)+EPS60)
              VXMIX=(RMAS1*VXIN1+RMAS2*VXIN2)/(RMAS1+RMAS2)
              VYIN1=BGKV(IUP2 ,IRAD)/(PDEN (IRAD)+EPS60)
              VYIN2=BGKV(IUP22,IRAD)/(PDEN2(IRAD)+EPS60)
              VYMIX=(RMAS1*VYIN1+RMAS2*VYIN2)/(RMAS1+RMAS2)
              VZIN1=BGKV(IUP3 ,IRAD)/(PDEN (IRAD)+EPS60)
              VZIN2=BGKV(IUP32,IRAD)/(PDEN2(IRAD)+EPS60)
              VZMIX=(RMAS1*VZIN1+RMAS2*VZIN2)/(RMAS1+RMAS2)
              ED1=(VXIN1**2+VYIN1**2+VZIN1**2)*FACT1
              ED2=(VXIN2**2+VYIN2**2+VZIN2**2)*FACT2
              VXIN(IPLS,IRAD)=VXMIX
              VYIN(IPLS,IRAD)=VYMIX
              VZIN(IPLS,IRAD)=VZMIX
              T1=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED1)/1.5
              T2=(EDEN2(IRAD)/(PDEN2(IRAD)+EPS60)-ED2)/1.5
              RM=2.D0*RMAS1*RMAS2/(RMAS1+RMAS2)**2
              TMIX=T1+RM*(T2-T1+
     .             FACT2/3.D0*((VXIN1-VXIN2)**2+(VYIN1-VYIN2)**2+
     .                         (VZIN1-VZIN2)**2))
              TIIN(IPLS,IRAD)=TMIX
              DEL=DIIN(IPLS,IRAD)-PDEN(IRAD)
              RATN=RATN+TBEL*DEL*VOL(IRAD)
              RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
              DIIN(IPLS,IRAD)=PDEN(IRAD)
              RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C             RRM=?
              RRE=RRE+EDEN(IRAD)*VOL(IRAD)
180     CONTINUE
c
c  same as do 180 loop, for additional cell region
c
        DO 190 IRAD=NSURF+1,NSURF+NRADD
C
c slmod begin - tr
c...temp: Ignore 1st additional cell (leave as vacuum)
          IF (IRAD.EQ.NSURF+1) CYCLE
c slmod end
          TBEL=0
          IF (LGVAC(IRAD,IPLS)) GOTO 191
          IF (NSTORDR >= NRAD) THEN
            TBEL = TABEL3(IREL,IRAD,1)
          ELSE
            KK=NREAEL(IREL)
            PLS=TIINL(IPLS,IRAD)+ADDEL(IREL,IPLS)
            TBEL = CREAC(9,1,KK)
            DO II=8,1,-1
              TBEL = TBEL*PLS + CREAC(II,1,KK)
            END DO
            TBEL=EXP(MAX(-100.D0,TBEL))*DIIN(IPLS,IRAD)
          END IF
191       CONTINUE
C         EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
C    .          DIIN(IPLS,IRAD)
          EOLD=(1.5*TIIN(IPLS,IRAD)+EDRIFT(IPLS,IRAD))*
     .          PDEN(IRAD)
          DOLD=DIIN(IPLS,IRAD)
          DEL=EOLD-EDEN(IRAD)
          RATE=RATE+TBEL*DEL*VOL(IRAD)
          DELX=BGKV(IUP1,IRAD)-VXIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
          DELY=BGKV(IUP2,IRAD)-VYIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
          DELZ=BGKV(IUP3,IRAD)-VZIN(IPLS,IRAD)*DIIN(IPLS,IRAD)
          RATM(1)=RATM(1)+TBEL*DELX*VOL(IRAD)
          RATM(2)=RATM(2)+TBEL*DELY*VOL(IRAD)
          RATM(3)=RATM(3)+TBEL*DELZ*VOL(IRAD)
C
          VXIN1=BGKV(IUP1 ,IRAD)/(PDEN (IRAD)+EPS60)
          VXIN2=BGKV(IUP12,IRAD)/(PDEN2(IRAD)+EPS60)
          VXMIX=(RMAS1*VXIN1+RMAS2*VXIN2)/(RMAS1+RMAS2)
          VYIN1=BGKV(IUP2 ,IRAD)/(PDEN (IRAD)+EPS60)
          VYIN2=BGKV(IUP22,IRAD)/(PDEN2(IRAD)+EPS60)
          VYMIX=(RMAS1*VYIN1+RMAS2*VYIN2)/(RMAS1+RMAS2)
          VZIN1=BGKV(IUP3 ,IRAD)/(PDEN (IRAD)+EPS60)
          VZIN2=BGKV(IUP32,IRAD)/(PDEN2(IRAD)+EPS60)
          VZMIX=(RMAS1*VZIN1+RMAS2*VZIN2)/(RMAS1+RMAS2)
          ED1=(VXIN1**2+VYIN1**2+VZIN1**2)*FACT1
          ED2=(VXIN2**2+VYIN2**2+VZIN2**2)*FACT2
          VXIN(IPLS,IRAD)=VXMIX
          VYIN(IPLS,IRAD)=VYMIX
          VZIN(IPLS,IRAD)=VZMIX
          T1=(EDEN(IRAD)/(PDEN(IRAD)+EPS60)-ED1)/1.5
          T2=(EDEN2(IRAD)/(PDEN2(IRAD)+EPS60)-ED2)/1.5
          RM=2.D0*RMAS1*RMAS2/(RMAS1+RMAS2)**2
          TMIX=T1+RM*(T2-T1+
     .         FACT2/3.D0*((VXIN1-VXIN2)**2+(VYIN1-VYIN2)**2+
     .                     (VZIN1-VZIN2)**2))
          TIIN(IPLS,IRAD)=TMIX
          DEL=DIIN(IPLS,IRAD)-PDEN(IRAD)
          RATN=RATN+TBEL*DEL*VOL(IRAD)
          RESN=RESN+TBEL*ABS(DEL)*VOL(IRAD)
          DIIN(IPLS,IRAD)=PDEN(IRAD)
          RRN=RRN+PDEN(IRAD)*VOL(IRAD)
C         RRM=?
          RRE=RRE+EDEN(IRAD)*VOL(IRAD)
C
190     CONTINUE
c
        ENDIF
c
        CALL LEER(2)
        WRITE (6,*) 'PARTICLE, MOMENTUM AND ENERGY EXCHANGE RATES '
        CALL LEER(1)
        WRITE (6,'(1X,A8,1X,1PE12.4)') 'RATN=   ',RATN*ELCHA
        WRITE (6,'(1X,A8,1X,3(1PE12.4))') 'RATM=   ',
     .                    RATM(1)*ELCHA*CNDYN,
     .                    RATM(2)*ELCHA*CNDYN,RATM(3)*ELCHA*CNDYN
        WRITE (6,'(1X,A8,1X,1PE12.4)') 'RATE=   ',RATE*ELCHA
        CALL LEER(2)
        WRITE (6,*) 'RESIDUA (1/SEC)'
        CALL LEER(1)
        CALL MASR1('RESN=   ',RESN/(RRN+EPS60))
C       CALL MASR3('RESM=                   ',RATM(1)/(RRM+EPS60),
C   .                    ,RATM(2)/(RRM+EPS60),RATM(3)/(RRM+EPS60))
        CALL MASR1('RESE=   ',RESE/(RRE+EPS60))
        CALL LEER(2)
c slmod begin - tr
        EIRRES(1,IPLS) = RATN*ELCHA
        EIRRES(2,IPLS) = RATM(1)*ELCHA*CNDYN
        EIRRES(3,IPLS) = RATM(2)*ELCHA*CNDYN
        EIRRES(4,IPLS) = RATM(3)*ELCHA*CNDYN
        EIRRES(5,IPLS) = RATE*ELCHA
        EIRRES(6,IPLS) = RESN/(RRN+EPS60)
        EIRRES(7,IPLS) = RESE/(RRE+EPS60)

        WRITE (97,*)
        WRITE (97,'(1X,A8,1X,1PE12.4)') 'RATN=   ',EIRRES(1,IPLS)
        WRITE (97,'(1X,A8,1X,3(1PE12.4))') 'RATM=   ',
     .                                     EIRRES(2,IPLS),
     .                      EIRRES(3,IPLS),EIRRES(4,IPLS)
        WRITE (97,'(1X,A8,1X,1PE12.4)') 'RATE=   ',EIRRES(5,IPLS)
        WRITE (97,*)
        WRITE (97,*) 'RESIDUA (1/SEC)'
        WRITE (97,*)
        WRITE (97,*) 'RESN=   ',RESN/(RRN+EPS60)
        WRITE (97,*) 'RESE=   ',RESE/(RRE+EPS60)
        WRITE (97,*)
c slmod end
C
C.................................................................
1000  CONTINUE
C.................................................................
c slmod begin - not tr
c...  Trying to get BGK rates set properly?
c      ENTRY MODBG1
c slmod end
      CALL LEER(2)
C
C  SAVE OVERHEAD, IF GEOMETRY DATA ALREADY AVAILABLE ON FILE
C
      IF (NFILEM.EQ.1) NFILEM=2
C
C  SET INDPRO=7, AND
C  WRITE PLASMA DATA ONTO RWK FOR CALL TO SUBR. PLASMA BELOW
C  TI,NI AND (VX,VY,VZ) FOR IPLS=1,NPLSI
C  PLAY SAVE: WRITE WHOLE RWK ARRAY.
C
      DO 500 I=1,5
        INDPRO(I)=7
500   CONTINUE
      NRWK1=(6+5*NPLS+NAIN)*NRAD
      IF (NID1 < NRWK1) THEN
        WRITE (6,*) ' RWK-ARRAY IS TOO SMALL TO HOLD PLASMA-DATA '
        WRITE (6,*) ' CHECK PARAMETER NSMSTRA '
        CALL EXIT
      END IF
      DO 510 I=1,NRWK1
        RWK(I)=0.
510   CONTINUE
      DO 550 IR=1,NXM
        DO 550 IP=1,NYM
          DO 550 IT=1,NZM
c slmod begin - tr
            IF (NLMLT) THEN
              IRAD=IR + (IP-1)*NR1P2 + (IT-1)*NSTRD
            ELSE
              IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
            ENDIF
c	    
c            IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
c slmod end
            RWK  ((0+0*NPLS)*NRAD+1   +(IRAD-1)*1   )= TEIN(IRAD)
            DO 520 IPLS=1,NPLSI
              RWK((1+0*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= TIIN(IPLS,IRAD)
              RWK((1+1*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= DIIN(IPLS,IRAD)
              RWK((1+2*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VXIN(IPLS,IRAD)
              RWK((1+3*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VYIN(IPLS,IRAD)
              RWK((1+4*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VZIN(IPLS,IRAD)
520         CONTINUE
            RWK  ((1+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BXIN(IRAD)
            RWK  ((2+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BYIN(IRAD)
            RWK  ((3+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BZIN(IRAD)
            RWK  ((5+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= VOL(IRAD)
            DO 530 IAIN=1,NAINI
              RWK((6+5*NPLS)*NRAD+IAIN+(IRAD-1)*NAIN)= ADIN(IAIN,IRAD)
530         CONTINUE
550   CONTINUE
C
c  same as do 550 loop , for additional cell region
c
      DO 570 IRAD=NSURF+1,NSURF+NRADD
        RWK  ((0+0*NPLS)*NRAD+1   +(IRAD-1)*1   )= TEIN(IRAD)
        DO 560 IPLS=1,NPLSI
          RWK((1+0*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= TIIN(IPLS,IRAD)
          RWK((1+1*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= DIIN(IPLS,IRAD)
          RWK((1+2*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VXIN(IPLS,IRAD)
          RWK((1+3*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VYIN(IPLS,IRAD)
          RWK((1+4*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VZIN(IPLS,IRAD)
560     CONTINUE
        RWK  ((1+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BXIN(IRAD)
        RWK  ((2+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BYIN(IRAD)
        RWK  ((3+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BZIN(IRAD)
        RWK  ((5+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= VOL(IRAD)
        DO 565 IAIN=1,NAINI
          RWK((6+5*NPLS)*NRAD+IAIN+(IRAD-1)*NAIN)= ADIN(IAIN,IRAD)
565     CONTINUE
570   CONTINUE
c slmod begin - tr
c...Also save additional cell data:
      DO IPLS=1,NPLSI
        DO IRAD=NSURF+1,NSURF+NRADD
          ASDSAV(IPLS,IRAD,1)=DIIN(IPLS,IRAD)
          ASDSAV(IPLS,IRAD,2)=TIIN(IPLS,IRAD)
          ASDSAV(IPLS,IRAD,3)=VXIN(IPLS,IRAD)
          ASDSAV(IPLS,IRAD,4)=VYIN(IPLS,IRAD)
          ASDSAV(IPLS,IRAD,5)=VZIN(IPLS,IRAD)
        ENDDO
      ENDDO
c slmod end
C .........................................................................
C  NOW: NEW COLLISION RATES MUST BE SET FOR THE NEXT ITERATION
C .........................................................................
C
C  IN CASE OF CROSS COLLISION, SOME MODIFICATIONS ON THE
C  BACKGROUND PARAMETERS ARE REQUIRED TEMPORARYLY TO ENFORCE
C  A SPECIFIC RELATION BETWEEN TAU_1,2 AND TAU_2,1
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE PROFILES
C  E.G.: EDRIFT, NEEDED FOR EPLEL3
C
      CALL PLASMA_DERIV
C
      TRCSAV=TRCAMD
      TRCAMD=.FALSE.
C
      CALL LEER(1)
      DO IPLS=1,NPLSI
        LMARK(IPLS)=.FALSE.
      ENDDO
      DO IPLS1=1,NPLSI
        IF (NPBGKP(IPLS1,2).NE.0) THEN
C  IPLS1 IS A CROSS COLLISION TALLY
C  FIND CORRESPONDING 2ND CROSS COLLISION TALLY
          IPLS2=0
          DO IPLS=1,NPLSI
            IF (ITYP1(IPLS).EQ.ITYP2(IPLS1).AND.
     .          ISPZ1(IPLS).EQ.ISPZ2(IPLS1).AND.
     .          ITYP2(IPLS).EQ.ITYP1(IPLS1).AND.
     .          ISPZ2(IPLS).EQ.ISPZ1(IPLS1)) IPLS2=IPLS
          ENDDO
          IF (IPLS2.EQ.0) GOTO 300
          CALL LEER(1)
          WRITE (6,*) 'CORRESPONDING CROSS COLLISION SPECIES '
          WRITE (6,*) 'IPLS1,IPLS2 ',IPLS1,IPLS2
          IF (LMARK(IPLS1).OR.LMARK(IPLS2)) GOTO 300
C  IPLS2 IS THE SECOND CROSS COLLISION TALLY
          WRITE (6,*) 'MODIFY PARAMETERS FOR CROSS COLLISIONALITIES '
          WRITE (6,*) 'IPLS1,IPLS2 ',IPLS1,IPLS2
          CALL LEER(1)
          LMARK(IPLS1)=.TRUE.
          LMARK(IPLS2)=.TRUE.
          DO IRAD=1,NSBOX
            DS1=DIIN(IPLS1,IRAD)
            DIIN(IPLS1,IRAD)=DIIN(IPLS2,IRAD)
            DIIN(IPLS2,IRAD)=DS1
C
            ENERGY(IPLS1,IRAD)=1.5*TIIN(IPLS1,IRAD)+EDRIFT(IPLS1,IRAD)
            ENERGY(IPLS2,IRAD)=1.5*TIIN(IPLS2,IRAD)+EDRIFT(IPLS2,IRAD)
C
            TS1=0.5*(TIIN(IPLS1,IRAD)+TIIN(IPLS2,IRAD))
            TIIN(IPLS1,IRAD)=TS1
            TIIN(IPLS2,IRAD)=TS1
          ENDDO
300       CONTINUE
        ENDIF
      ENDDO
      CALL LEER(2)
C
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE MODIFIED PROFILES
C
      CALL PLASMA_DERIV
C
C  RESET BGK-ATOMIC AND MOLECULAR DATA ARRAYS
C
      DO IPLS=1,NPLSI
        IF (NPBGKP(IPLS,1).NE.0) THEN
          ITYP=ITYP1(IPLS)
          ISPZ=ISPZ1(IPLS)
          IREL=IREL1(IPLS)
          NRC=INRC1(IPLS)
          IF (ITYP.EQ.1) THEN
            ISP=ISPZ
            KK=IREACA(ISPZ,NRC)
            EBULK=EBULKA(ISPZ,NRC)
            ISCDE=ISCDEA(ISPZ,NRC)
            IESTM=IESTMA(ISPZ,NRC)
            FACTKK=FREACA(ISPZ,NRC)
          ELSEIF (ITYP.EQ.2) THEN
            ISP=NATMI+ISPZ
            KK=IREACM(ISPZ,NRC)
            EBULK=EBULKM(ISPZ,NRC)
            ISCDE=ISCDEM(ISPZ,NRC)
            IESTM=IESTMM(ISPZ,NRC)
            FACTKK=FREACM(ISPZ,NRC)
          ELSEIF (ITYP.EQ.3) THEN
            ISP=NATMI+NMOLI+ISPZ
            KK=IREACI(ISPZ,NRC)
            EBULK=EBULKI(ISPZ,NRC)
            ISCDE=ISCDEI(ISPZ,NRC)
            IESTM=IESTMI(ISPZ,NRC)
            FACTKK=FREACI(ISPZ,NRC)
          ENDIF
          IF (FACTKK.EQ.0.D0) FACTKK=1.D0
C  BGK COLLISION, RESET TABEL3, EPLEL3
          CALL XSTEL(IREL,ISP,IPLS,EBULK,ISCDE,IESTM,KK,FACTKK)
          IF (NPBGKP(IPLS,2).NE.0) THEN
C  CROSS COLLISION, RESET EPLEL3 FOR TRACKLENGTH ESTIMATOR
            IF (NSTORDR >= NRAD) THEN
              DO J=1,NSBOX
                EPLEL3(IREL,J,1)=ENERGY(IPLS,J)
              ENDDO
            ELSE
              NELREL(IREL)=-3
            END IF
          ENDIF
        ENDIF
      ENDDO
C
      TRCAMD=TRCSAV
C
C
C  RESTORE PLASMA DATA FROM RWK ARRAY
C
      CALL PLASMA
c slmod begin - tr
c...This is sort of unnecessary, since INDPRO=7 will restore the data
c   on additional cells using the RWK array, but presently I am unsure
c   of changing from INDPRO=6.
      DO IPLS=1,NPLSI
        DO IRAD=NSURF+1,NSURF+NRADD
          DIIN(IPLS,IRAD)=ASDSAV(IPLS,IRAD,1)
          TIIN(IPLS,IRAD)=ASDSAV(IPLS,IRAD,2)
          VXIN(IPLS,IRAD)=ASDSAV(IPLS,IRAD,3)
          VYIN(IPLS,IRAD)=ASDSAV(IPLS,IRAD,4)
          VZIN(IPLS,IRAD)=ASDSAV(IPLS,IRAD,5)
        ENDDO
      ENDDO
c slmod end
      CALL PLASMA_DERIV
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C
      NFILEL=3

      CALL WRPLAM(TRCFLE)
C
      DEALLOCATE (PDEN)
      DEALLOCATE (EDEN)
      DEALLOCATE (PDEN2)
      DEALLOCATE (EDEN2)
      DEALLOCATE (ENERGY)
C
      RETURN
C
995   CONTINUE
      WRITE (6,*) 'SPECIES ERROR IN MODBGK'
      CALL EXIT
C
999   CONTINUE
      WRITE (6,*) 'ERROR IN MODBGK. IPLS,IBGK= ',IPLS,IBGK1,IBGK2
      CALL EXIT
      END
c
      SUBROUTINE STATIS_BGK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CCONA'
      INCLUDE 'COMUSR'
      INCLUDE 'COUTAU'
      INCLUDE 'CGRID'
      INCLUDE 'CESTIM'
      INCLUDE 'CSDVI'
      INCLUDE 'CSDVI_BGK'
      LOGICAL LP,LT
C
      SAVE
C
      ENTRY STATS0_BGK
C
      RETURN
C
      ENTRY STATS1_BGK(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
C
      IF (NBGVI.EQ.0) RETURN
      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
C
      DO 1012 IBGV=1,NBGVI
C
        IF (LMETSP(NSPAN(NTALB)+IBGV-1)) THEN
          SD1S=0.
          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)
            SD1=BGKV(IBGV,IR)-SDVIA_BGK(IBGV,IR)
            SD1S=SD1S+SD1
            SDVIA_BGK(IBGV,IR)=BGKV(IBGV,IR)
            SIGMA_BGK(IBGV,IR)=SIGMA_BGK(IBGV,IR)+SD1*SD1
          END DO
          SGMS_BGK(IBGV)=SGMS_BGK(IBGV)+SD1S*SD1S
        END IF
1012  CONTINUE

C  STATISTICS FOR PDENA AND EDENA
      IBGV = NBGVI
      DO IAT=1,NATMI
C
        IF (LMETSP(NSPAN(1)+IAT-1)) THEN
          IBGV=NBGVI+IAT
          SD1S=0.
          SD2S=0.
          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)
            SD1=PDENA(IAT,IR)-SDVIA_BGK(IBGV,IR)
            SD1S=SD1S+SD1
            SDVIA_BGK(IBGV,IR)=PDENA(IAT,IR)
            SIGMA_BGK(IBGV,IR)=SIGMA_BGK(IBGV,IR)+SD1*SD1

            SD2=EDENA(IAT,IR)-SDVIA_BGK(IBGV+NATMI,IR)
            SD2S=SD2S+SD2
            SDVIA_BGK(IBGV+NATMI,IR)=EDENA(IAT,IR)
            SIGMA_BGK(IBGV+NATMI,IR)=SIGMA_BGK(IBGV+NATMI,IR)+SD2*SD2
          END DO
          SGMS_BGK(IBGV)=SGMS_BGK(IBGV)+SD1S*SD1S
          SGMS_BGK(IBGV+NATMI)=SGMS_BGK(IBGV+NATMI)+SD2S*SD2S
        END IF
      END DO
      IBGV = NBGVI+2*NATMI

C  STATISTICS FOR PDENM AND EDENM
      DO IMO=1,NMOLI
C
        IF (LMETSP(NSPAN(2)+IMO-1)) THEN
          IBGV=NBGVI+2*NATMI+IMO
          SD1S=0.
          SD2S=0.
          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)
            SD1=PDENM(IMO,IR)-SDVIA_BGK(IBGV,IR)
            SD1S=SD1S+SD1
            SDVIA_BGK(IBGV,IR)=PDENM(IMO,IR)
            SIGMA_BGK(IBGV,IR)=SIGMA_BGK(IBGV,IR)+SD1*SD1

            SD2=EDENM(IMO,IR)-SDVIA_BGK(IBGV+NMOLI,IR)
            SD2S=SD2S+SD2
            SDVIA_BGK(IBGV+NMOLI,IR)=EDENM(IMO,IR)
            SIGMA_BGK(IBGV+NMOLI,IR)=SIGMA_BGK(IBGV+NMOLI,IR)+SD2*SD2
          END DO
          SGMS_BGK(IBGV)=SGMS_BGK(IBGV)+SD1S*SD1S
          SGMS_BGK(IBGV+NMOLI)=SGMS_BGK(IBGV+NMOLI)+SD2S*SD2S
        END IF
      END DO
      IBGV = IBGV+NMOLI
C
C
1020  CONTINUE
      RETURN
C
      ENTRY STATS2_BGK(XN,FSIG,ZFLUX)
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1  MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NBGVI.EQ.0) GOTO 2200
C
      DO 2112 IBGV=1,NBGVI
        DS=SUM(BGKV(IBGV,1:NSB))
        DO 2111 IR=1,NSB
          D=BGKV(IBGV,IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0.D0,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
2111    CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0.D0,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
2112  CONTINUE
C
C   STATISTICS FOR PDENA
      IBGV=NBGVI
      DO IAT=1,NATMI
        IBGV=IBGV+1
        DS=SUM(PDENA(IAT,1:NSB))
        DO IR=1,NSB
          D=PDENA(IAT,IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0.D0,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0.D0,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
C   STATISTICS FOR EDENA
      DO IAT=1,NATMI
        IBGV=IBGV+1
        DS=SUM(EDENA(IAT,1:NSB))
        DO IR=1,NSB
          D=EDENA(IAT,IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0.D0,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0.D0,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
C   STATISTICS FOR PDENM
      DO IMO=1,NMOLI
        IBGV=IBGV+1
        DS=SUM(PDENM(IMO,1:NSB))
        DO IR=1,NSB
          D=PDENM(IMO,IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0.D0,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0.D0,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
C   STATISTICS FOR EDENM
      DO IMO=1,NMOLI
        IBGV=IBGV+1
        DS=SUM(EDENM(IMO,1:NSB))
        DO IR=1,NSB
          D=EDENM(IMO,IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0.D0,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0.D0,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
2200  CONTINUE
      RETURN
      END
