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
      SUBROUTINE UPTBGK(XSTOR2,XSTORV2,WV,NPBGK,IFLAG)
C
C  UPDATE BGK-SPECIFIC TALLIES, TRACKLENGTH ESTIMATORS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CUPD
      USE CGRID
      USE COMPRT
      USE CTEXT
      USE CSDVI
      USE COMXS
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                      XSTORV2(NSTORV,N2ND+N3RD)
      REAL(DP), INTENT(IN) :: WV
      INTEGER, INTENT(IN) :: NPBGK, IFLAG
      REAL(DP) :: DIST, WTRV, WTRVX, WTRVY, WTRVZ
      INTEGER :: I, NMTSP, IUPD2, IUPD3, IRD, IFIRST, NSBGK, IBGK,
     .           IML, IIO, IUPD1, ITP, ISP, IAT
      CHARACTER(8) :: TXT
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
              TXT=TEXTS(NSPH+IAT)
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
            IF (NPBGKI(ISP).EQ.IBGK) THEN
              ITP=3
              IIO=ISP
              TXT=TEXTS(NSPAM+IIO)
              GOTO 1
            ENDIF
          ENDDO
          WRITE (6,*) 'SPECIES ERROR IN UPTBGK'
          CALL EXIT_OWN(1)
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

        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI+NCPVI
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
