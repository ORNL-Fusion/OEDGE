cdr  28.4.04:  nhsts(ispz) option connected (to select species
cdr            for trajectory plot. see modification to input.f, 28.4.04
      SUBROUTINE PLT2D

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CRECH
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CPL3D
      USE CPLMSK
      USE CPLOT
      USE CINIT
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CTRCEI
      USE CGEOM
      USE CTETRA
      USE COMPRT
      USE COMSOU
      USE COMSPL
      USE CTEXT
      USE CLGIN
      USE CTRIG

      IMPLICIT NONE
C
      INTEGER,PARAMETER :: NTXHST=18

      REAL(DP) :: XX(101), YY(101)
      REAL(DP) :: DSD(3), AFF(3,3), AFFI(3,3)
      REAL(DP) :: TET(4,3), EBENE(4), CTPNTS(4,3)
      REAL(DP) :: XNP05, YYIA, XN1, YN, FX, FY, XN2, Z1, XN3, YNP, XN, 
     .          R, SLT, ZW1, ZW2, YWO, XWO, XR, THET, PPHI, XPL, RWN,
     .          YPL, ZPL, ZWN, YW1, YW2, SCLFCY, YMI2D, YMA2D, SCLFCX,
     .          X, Y, Z, XW2, DM, RR, XW1, ABSMAX, ORDMAX, XPLO, YPLO,
     .          ZPLO, XNULL, XMI2D, XMA2D, YNULL, XWN, YWN, YT2, TESTN,
     .          XT2, XTIP, P, YTIP, TR, RS, EP, EL, XT, YT, XTN, YTN 
      REAL(SP) :: XPS(5), YPS(5)
      INTEGER :: ISPL(NTXHST),ICLR(2*NSTS+1),IDSH(2*NSTS+1),
     .           ISWC(2*NSTS+1)
      INTEGER, ALLOCATABLE :: IFARB(:,:), IDASH(:,:)
      INTEGER :: ICP, ISTR, IC, IC1, IC2, NCTPNT, IDUMMY, ICT, IEN, 
     .           ITH, ITHPL, IAN, NTDUM, NTT, IECKE2,
     .           LEARCA, NT, ICOLOR, NU, J, ISYM, IFLAG, IERR, IWRIT,
     .           IR, IP, ISTS, IN, IY, IB, IT, IA, NRET, I, NSW, ISW,
     .           ISP, IHELP, K, IFL, ISYM_ERR, IRA, IRE, IPA, IPE, 
     .           ITA, ITE
      LOGICAL :: PLSAV1, PLSAV2, LSTORE
      CHARACTER(20) :: TXTHST(NTXHST)
      CHARACTER(10) :: CX, CY, CX0, CY0, CZ0
      CHARACTER(4) :: CH

      SAVE
      DATA ABSMAX,ORDMAX/21.,21./
      DATA XNULL,YNULL/9.,4./,XWN,YWN/0.,0./
      DATA IWRIT/0/,ISPL/2,3,12,13,4,5,6,14,15,16,17,0,0,0,7,8,9,48/
      DATA TXTHST
     .           /'LOCATE(1)           ',
     .            'ELECTR. IMPACT(2)   ',
     .            'HEAVY PAR. IMPACT(3)',
     .            'PHOTON IMPACT(4)    ',
     .            'ELASTIC COLL.(5)    ',
     .            'CHARGE EXCHANGE(6)  ',
     .            'FOKKER PLANCK(7)    ',
     .            'SURFACE(8)          ',
     .            'SPLITTING(9)        ',
     .            'RUSSIAN ROULETTE(10)',
     .            'PERIODICITY(11)     ',
     .            'RESTART:A. SPLT.(12)',
     .            'SAVE:COND. EXP.(13) ',
     .            'RESTART:COND EXP(14)',
     .            'TIME LIMIT(15)      ',
     .            'GENERATION LIMIT(16)',
     .            'FLUID LIMIT(17)     ',
     .            'ERROR DETECTED      '/
C
C  SYMBOL FOR PARTICLE TRACING ERROR
      ISYM_ERR=NTXHST
C
C  PREPARE PLOT OF STANDARD-MESH AND ADDITIONAL SURFACES
C
      IF (.NOT.NLPL2D) GOTO 300
C
      CALL FTCRE(CH2MX,CX)
      CALL FTCRE(CH2MY,CY)
      CALL FTCRE(CH2X0,CX0)
      CALL FTCRE(CH2Y0,CY0)
      CALL FTCRE(CH2Z0,CZ0)
C
C  CO-ORDINATES FOR TEXT FOR HISTORIES ON PLOT, TOP LEFT,
C  AND SCALING FACTORS
C  THESE 4 NUMBERS MAY BE MODIFIED IN SUBR. PLT3D
      XN2D=-8.
      YN2D=12.
      FX2D=1.
      FY2D=1.
C
C
C  NULLPUNKT AUF DEM PAPIER
      X0PL=XNULL
      Y0PL=YNULL
C  ACHSENLAENGEN
      LENX=ABSMAX
      LENY=ORDMAX
C  ACHSENUNTERTEILUNG VORGEGEBEN?
      STPSZX=0.2
      STPSZY=0.2
      INTNRX=10
      INTNRY=10
C  ACHSE LOGARITHMISCH?
      LOGX=.FALSE.
      LOGY=.FALSE.
C  ZEICHNE NETZLINIEN EIN
      GRIDX=.TRUE.
      GRIDY=.TRUE.
C  MACHE GRADE GRENZEN
      FITX=.FALSE.
      FITY=.FALSE.

      MINX=-1.
      MAXX=1.
      MINY=-1.
      MAXY=1.
C
C  PLOT X AND Y AXIS
C
      CALL GRNXTB(1)
      CALL PLTMSK(IERR)
C
      CALL GRSPTS (16)
      CALL GRCHRC (0.3,0.,16)
      CALL GRSCLV (0.,0.,REAL(ABSMAX,KIND(1.E0)),
     .                   REAL(ORDMAX,KIND(1.E0)))
C
      CALL GRTXT (-8.,24.,72,TXTRUN)
      CALL GRTXT (-8.,22.,15,'SCALING FACTORS')
      CALL GRTXT (-8.,21.,7,'FACT-X=')
      CALL GRTXT (-8.,20.,7,'FACT-Y=')
      CALL GRTXT (-5.3,21.,10,CX)
      CALL GRTXT (-5.3,20.,10,CY)
      CALL GRTXT (-8.,18.,15,'ORIGIN         ')
      CALL GRTXT (-8.,17.,7,'CH2X0= ')
      CALL GRTXT (-8.,16.,7,'CH2Y0= ')
      CALL GRTXT (-5.3,17.,10,CX0)
      CALL GRTXT (-5.3,16.,10,CY0)
      CALL GRTXT (-8.,14.,10,'PLOTTED AT')
      CALL GRTXT (-8.,13.,4,'Z = ')
      CALL GRTXT (-5.3,13.,10,CZ0)
      XMI2D=CH2X0-CH2MX
      XMA2D=CH2X0+CH2MX
      YMI2D=CH2Y0-CH2MY
      YMA2D=CH2Y0+CH2MY
      CALL GRSCLV(REAL(XMI2D,KIND(1.E0)),REAL(YMI2D,KIND(1.E0)),
     .            REAL(XMA2D,KIND(1.E0)),REAL(YMA2D,KIND(1.E0)))
C
C  SCALE FACTORS: USER CO-ORDINATES TO CM:
C  X-DIRECTION:
      SCLFCX=ABSMAX/(XMA2d-XMI2d)
C  Y-DIRECTION:
      SCLFCY=ORDMAX/(YMA2d-YMI2d)
C
C   PLOT R-GRID
C
C  PLOT RADIAL GRID
C
      IF (.NOT.PL1ST.OR.NR1ST.LE.1.OR.(PLCUT(1).AND.(LEVGEO.NE.5)))
     .    GOTO 170
C
      IF (LEVGEO.EQ.1) THEN
C  X-Y-PLANE
        IF (PLCUT(3)) THEN
          ALLOCATE (IFARB(N1ST,N2ND))
          ALLOCATE (IDASH(N1ST,N2ND))
          IFARB = 1
          IDASH = 0
          DO J=1,NSTSI
            NU = INUMP(J,1)
            IF (NU <= 0) CYCLE 
            IPA = IRPTA(J,2)
            IPE = IRPTE(J,2)-1
            IFARB(NU,IPA:IPE) = ILCOL(NLIM+J)
            IDASH(NU,IPA:IPE) = ILIIN(NLIM+J)
            IF (NLSPLT(NU)) IFARB(NU,:) = 2
          END DO
C  X-Z-PLANE
        ELSEIF (PLCUT(2)) THEN
          ALLOCATE (IFARB(N1ST,N3RD))
          ALLOCATE (IDASH(N1ST,N3RD))
          IFARB = 1
          IDASH = 0
          DO J=1,NSTSI
            NU = INUMP(J,1)
            IF (NU <= 0) CYCLE 
            ITA = IRPTA(J,3)
            ITE = IRPTE(J,3)-1
            IFARB(NU,ITA:ITE) = ILCOL(NLIM+J)
            IDASH(NU,ITA:ITE) = ILIIN(NLIM+J)
            IF (NLSPLT(NU)) IFARB(NU,:) = 2
          END DO
        ENDIF
        DO NU=NPLINR,NPLOTR,NPLDLR
          LSTORE = .FALSE.
          XW1=RSURF(NU)
          IF (PLCUT(3)) THEN
            IF (NLTRA) XW1=XW1+RMTOR
            IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
              IF (NLPOL) THEN
                DO IP=NPLINP, MIN(NPLOTP-1,NP2NDM), NPLDLP
                  IF (IFARB(NU,IP) == 666) CYCLE
                  XX(1)=XW1
                  XX(2)=XW1
                  YY(1)=PSURF(IP)
                  YY(2)=PSURF(IP+1)
                  CALL GRNWPN(IFARB(NU,IP))
                  IF (IDASH(NU,IP) <= 0) THEN
                    CALL GRDSH(0.2,0.5,0.2)
                  ELSE
                    CALL GRDSH(1.,0.,1.)
                  END IF
                  CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
                END DO
              ELSE 
                XX(1)=XW1
                XX(2)=XW1
                YY(1)=YIA
                YY(2)=YAA
                CALL GRNWPN(IFARB(NU,1))
                IF (IDASH(NU,1) <= 0) THEN
                  CALL GRDSH(0.2,0.5,0.2)
                ELSE
                  CALL GRDSH(1.,0.,1.)
                END IF
                CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
              END IF  
            ENDIF
          ELSEIF (PLCUT(2)) THEN
            IF (NLTRZ) THEN
              IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
                IF (NLTOR) THEN
                  DO IT=NPLINT, MIN(NPLOTT-1,NT3RDM), NPLDLT
                    IF (IFARB(NU,IT) == 666) CYCLE
                    XX(1)=XW1
                    XX(2)=XW1
                    YY(1)=ZSURF(IT)
                    YY(2)=ZSURF(IT+1)
                    CALL GRNWPN(IFARB(NU,IT))
                    IF (IDASH(NU,IT) <= 0) THEN
                      CALL GRDSH(0.2,0.5,0.2)
                    ELSE
                      CALL GRDSH(1.,0.,1.)
                    END IF
                    CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
                  END DO
                ELSE ! .NOT.NLTOR
                  XX(1)=XW1
                  XX(2)=XW1
                  YY(1)=ZIA
                  YY(2)=ZAA
                  CALL GRNWPN(IFARB(NU,1))
                  IF (IDASH(NU,1) <= 0) THEN
                    CALL GRDSH(0.2,0.5,0.2)
                  ELSE
                    CALL GRDSH(1.,0.,1.)
                  END IF
                  CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
                END IF  
              ENDIF
            ELSEIF (NLTRA) THEN
              X=RMTOR+XW1
              Z=TANAL*X
              RR=SQRT(X*X+Z*Z)
              IF (NTTRA.GT.100) THEN
                WRITE (6,*) 'ERROR IN PLT2D '
                CALL EXIT_OWN(1)
              ENDIF
              CALL GRNWPN(IFARB(NU,1))
              IF (IDASH(NU,1) <= 0) THEN
                 CALL GRDSH(0.2,0.5,0.2)
              ELSE
                 CALL GRDSH(1.,0.,1.)
              END IF
              DO J=1,NTTRA+1
                XX(J)=RR*COS((J-1)*2.*ALPHA)
                YY(J)=RR*SIN((J-1)*2.*ALPHA)
              ENDDO
              CALL PLTLNE(NTTRA+1,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
            ENDIF
          ENDIF
        END DO

        call grnwpn(1)
        do j=1,nsurf
          if (nstgrd(j).eq.1) then
            call ncelln(j,ir,ip,it,ia,ib,nr1st,np2nd,nt3rd,nbmlt,
     .                  nlrad,nlpol,nltor)
            IF (IR >= NR1ST) CYCLE
            IF (PLCUT(3)) THEN
              IF (IP >= NP2ND) CYCLE
              XPS(1)=RSURF(IR)
              YPS(1)=PSURF(IP)
              XPS(2)=RSURF(IR)
              YPS(2)=PSURF(IP+1)
              XPS(3)=RSURF(IR+1)
              YPS(3)=PSURF(IP+1)
              XPS(4)=RSURF(IR+1)
              YPS(4)=PSURF(IP)
              XPS(5)=RSURF(IR)
              YPS(5)=PSURF(IP)
            ELSE IF (PLCUT(2)) THEN
              IF (IT >= NT3RD) CYCLE
              XPS(1)=RSURF(IR)
              YPS(1)=ZSURF(IT)
              XPS(2)=RSURF(IR)
              YPS(2)=ZSURF(IT+1)
              XPS(3)=RSURF(IR+1)
              YPS(3)=ZSURF(IT+1)
              XPS(4)=RSURF(IR+1)
              YPS(4)=ZSURF(IT)
              XPS(5)=RSURF(IR)
              YPS(5)=ZSURF(IT)
            END IF
            IF (NLTRA) XPS(1:5)=XPS(1:5)+RMTOR
            CALL GRFILL(5,XPS,YPS,1,1)
          ENDIF
        enddo
        
        IF (ALLOCATED(IFARB)) THEN
          DEALLOCATE (IFARB)
          DEALLOCATE (IDASH)
        END IF
        
!        DO 144 NU=NPLINR,NPLOTR,NPLDLR
!          LSTORE = .FALSE.
!          DO 145 J=1,NSTSI
!            IF (NU.EQ.INUMP(J,1)) THEN
!              IF (ILCOL(NLIM+J) == 666) GOTO 144
!              CALL GRNWPN(ILCOL(NLIM+J))
!              LSTORE = PLSTOR
!              INOSF=NLIM+J
!              IF (ILIIN(NLIM+J).LE.0) GOTO 146
!              CALL GRDSH(1.,0.,1.)
!              GOTO 147
!            ENDIF
!145       CONTINUE
!146       CALL GRDSH(0.2,0.5,0.2)
!147       CONTINUE
!          IF (NLSPLT(NU)) CALL GRNWPN(2)
!          XW1=RSURF(NU)
!          IF (PLCUT(3)) THEN
!            IF (NLTRA) XW1=XW1+RMTOR
!            IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
!              XX(1)=XW1
!              XX(2)=XW1
!              YY(1)=YW1
!              YY(2)=YW2
!              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
!            ENDIF
!          ELSEIF (PLCUT(2)) THEN
!            IF (NLTRZ) THEN
!              IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
!                XX(1)=XW1
!                XX(2)=XW1
!                YY(1)=YW1
!                YY(2)=YW2
!                CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
!              ENDIF
!            ELSEIF (NLTRA) THEN
!              X=RMTOR+XW1
!              Z=TANAL*X
!              RR=SQRT(X*X+Z*Z)
!              IF (NTTRA.GT.100) THEN
!                WRITE (6,*) 'ERROR IN PLT2D '
!                CALL EXIT_OWN(1)
!              ENDIF
!              DO J=1,NTTRA+1
!                XX(J)=RR*COS((J-1)*2.*ALPHA)
!                YY(J)=RR*SIN((J-1)*2.*ALPHA)
!              ENDDO
!              CALL PLTLNE(NTTRA+1,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
!            ENDIF
!          ENDIF
!          CALL GRNWPN(1)
!144     CONTINUE

        CALL GRDSH(1.,0.,1.)
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        IF (PLCUT(2)) THEN
C
C  X-Z-PLANE
          Y=CH2Z0
          IF (NLTOR) THEN
            IF (NLTRZ) THEN
              YW1=ZSURF(NPLINT)
              YW2=ZSURF(NPLOTT)
            ELSEIF (NLTRA) THEN
              YW1=ZSURF(NPLINT)
              YW2=ZSURF(NPLOTT)
            ELSEIF (NLTRT) THEN
              GOTO 990
            ENDIF
          ELSEIF (.NOT.NLTOR) THEN
            IF (NLTRZ) THEN
              YW1=ZIA
              YW2=ZAA
            ELSEIF (NLTRA) THEN
              YW1=ZIA*DEGRAD
              YW2=ZAA*DEGRAD
            ELSEIF (NLTRT) THEN
              GOTO 990
            ENDIF
          ENDIF
          DO 134 NU=NPLINR,NPLOTR,NPLDLR
            LSTORE=.FALSE.
            DO 135 J=1,NSTSI
              IF (NU.EQ.INUMP(J,1)) THEN
                IF (ILCOL(NLIM+J) == 666) GOTO 134
                CALL GRNWPN(ILCOL(NLIM+J))
                LSTORE = PLSTOR
                INOSF=NLIM+J
                IF (ILIIN(NLIM+J).LE.0) GOTO 136
                CALL GRDSH(1.,0.,1.)
                GOTO 137
              ENDIF
135         CONTINUE
136         CALL GRDSH(0.2,0.5,0.2)
137         CONTINUE
            IF (NLSPLT(NU)) CALL GRNWPN(2)
C  PLOT AT Y=CH2Z0 , TO BE WRITTEN. PRESENTLY AT Y=0
            IF (CH2Z0.NE.0.) THEN
              WRITE (6,*) 'PLOTOPTION CH2Z0.NE.0 NOT READY. EXIT '
              CALL EXIT_OWN(1)
            ENDIF
C  FIND X= CONST. LINES AT Y=CH2Z0: XW1, XW2
            XW1=RSURF(NU)+EP1(NU)
            XW2=-RSURF(NU)+EP1(NU)
C
            IF (NLTRZ) THEN
              XX(1)=XW1
              XX(2)=XW1
              YY(1)=YW1
              YY(2)=YW2
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
              XX(1)=XW2
              XX(2)=XW2
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
            ELSEIF (NLTRA) THEN
              X=RMTOR+XW1
              Z=TANAL*X
              RR=SQRT(X*X+Z*Z)
              IF (NTTRA.GT.100) THEN
                WRITE (6,*) 'ERROR IN PLT2D '
                CALL EXIT_OWN(1)
              ENDIF
              DO J=1,NTTRA
                XX(J)=RR*COS(ZSURF(J))
                YY(J)=RR*SIN(ZSURF(J))
              ENDDO
              CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
              X=RMTOR+XW2
              Z=TANAL*X
              RR=SQRT(X*X+Z*Z)
              IF (NTTRA.GT.100) THEN
                WRITE (6,*) 'ERROR IN PLT2D '
                CALL EXIT_OWN(1)
              ENDIF
              DO J=1,NTTRA
                XX(J)=RR*COS(ZSURF(J))
                YY(J)=RR*SIN(ZSURF(J))
              ENDDO
              CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
            ENDIF
            CALL GRNWPN(1)
134       CONTINUE
          CALL GRDSH(1.,0.,1.)
C
        ELSEIF (PLCUT(3)) THEN
C
C    X-Y PLOT
          DO 140 NU=NPLINR,NPLOTR,NPLDLR
            LSTORE = .FALSE.
            DO 141 J=1,NSTSI
              IF (NU.EQ.INUMP(J,1)) THEN
                IF (ILCOL(NLIM+J) == 666) GOTO 140
                CALL GRNWPN(ILCOL(NLIM+J))
                LSTORE = PLSTOR
                INOSF=NLIM+J
                IF (ILIIN(NLIM+J).LE.0) GOTO 142
                CALL GRDSH(1.,0.,1.)
                GOTO 143
              ENDIF
141         CONTINUE
142         CALL GRDSH(0.2,0.5,0.2)
143         CONTINUE
            IF (NLSPLT(NU)) CALL GRNWPN(2)
            DM=0.1
            RS=RSURF(NU)
            EP=EP1(NU)
            EL=ELL(NU)
            TR=TRI(NU)
            CALL PLGELR (RS,EP,EL,TR,DM,100,XX,YY,NRET,PSURF,NP2ND)
            IF (NLTRA) THEN
              DO I=1,NRET
                XX(I)=XX(I)+RMTOR
              ENDDO
            ENDIF
            CALL PLTLNE(NRET,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
            CALL GRNWPN(1)
140       CONTINUE
          CALL GRDSH(1.,0.,1.)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
        IF (PLCUT(3)) THEN
C
          DO 158 NU=NPLINR,NPLOTR,NPLDLR
            IF (NLSPLT(NU)) CALL GRNWPN(2)
C  NSW/2 VERSCHIEDENE STUECKE ZU FLAECHE NU, NSW GERADE
            NSW=0
            DO 1561 J=1,NSTSI
C
              IF (NU.EQ.INUMP(J,1)) THEN
                IF (ILCOL(NLIM+J) == 666) CYCLE
                NSW=NSW+1
                ICLR(NSW)=ILCOL(NLIM+J)
                IDSH(NSW)=1
                IF (ILIIN(NLIM+J).GT.0) IDSH(NSW)=0
                ISWC(NSW)=MAX(1,IRPTA(J,2))
                ICLR(NSW+1)=1
                IDSH(NSW+1)=1
                ISWC(NSW+1)=MIN(NPOINT(2,NPPLG),IRPTE(J,2))
                NSW=NSW+1
              ENDIF
1561        CONTINUE
C
C  SORTIEREN, NACH "POLOIDALEM BEGIN" < "POLOIDALEM ENDE"
            IF (NSW.EQ.0) GOTO 1597
157         ISW=0
            DO 159 I=1,NSW-3,2
              IF (ISWC(I).GT.ISWC(I+2)) THEN
                ISW=ISW+1

                IHELP=ISWC(I)
                ISWC(I)=ISWC(I+2)
                ISWC(I+2)=IHELP
                IHELP=ISWC(I+1)
                ISWC(I+1)=ISWC(I+3)
                ISWC(I+3)=IHELP

                IHELP=ICLR(I)
                ICLR(I)=ICLR(I+2)
                ICLR(I+2)=IHELP
                IHELP=ICLR(I+1)
                ICLR(I+1)=ICLR(I+3)
                ICLR(I+3)=IHELP

                IHELP=IDSH(I)
                IDSH(I)=IDSH(I+2)
                IDSH(I+2)=IHELP
                IHELP=IDSH(I+1)
                IDSH(I+1)=IDSH(I+3)
                IDSH(I+3)=IHELP
              ENDIF
159         CONTINUE
            IF (ISW.GT.0) GOTO 157
C  SORTIEREN FERTIG
C
            I=0
1581        I=I+1
1585        IF (ISWC(I).EQ.ISWC(I+1)) THEN
              ICLR(I)=ICLR(I+1)
              IDSH(I)=IDSH(I+1)
              DO 1596 J=I+2,NSW
                ISWC(J-1)=ISWC(J)
                ICLR(J-1)=ICLR(J)
                IDSH(J-1)=IDSH(J)
1596          CONTINUE
              NSW=NSW-1
              IF (I.LT.NSW) GOTO 1585
            ENDIF
            IF (I.LT.NSW-1) GOTO 1581
C
1597        CONTINUE
            NSW=NSW+1
            ISWC(NSW)=100000
            ICLR(NSW)=1
            IDSH(NSW)=1
C
            ISW=1
            CALL GRDSH (0.2,0.5,0.2)
            DO 150 I=1,NPPLG
              LSTORE = .FALSE.
              DO 150 K=NPOINT(1,I),NPOINT(2,I)
                IFL=K-NPOINT(1,I)
                XTN=XPOL(NU,K)
                IF (NLTRA) XTN=XTN+RMTOR
C               IF (NLTRT) XTN=XTN+RMTOR
                YTN=YPOL(NU,K)
                CALL TSTCHM(IFL,XTN,YTN,XT,YT,IN,TESTN,
     .                      XMI2D,XMA2D,YMI2D,YMA2D,XT2,YT2)
C
                IF (IFL.EQ.0.AND.K.EQ.ISWC(ISW)) THEN
                  LSTORE = PLSTOR
                  INOSF=0
                  IF (.NOT.NLSPLT(NU)) CALL GRNWPN (ICLR(ISW))
                  IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
                  IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
                  ISW=ISW+1
                ENDIF
C
                IF (IN.EQ.0) IN=4
                GOTO (151,152,153,154,155),IN
151               CONTINUE
                    CALL GRDRW (REAL(XTN,KIND(1.E0)),
     .                          REAL(YTN,KIND(1.E0)))
                    IF (LSTORE) CALL STCOOR(XTN,YTN,1)
                    GOTO 156
152               CONTINUE
                    CALL GRDRW(REAL(XT,KIND(1.E0)),
     .                         REAL(YT,KIND(1.E0)))
                    CALL GRJMP(REAL(XTN,KIND(1.E0)),
     .                         REAL(YTN,KIND(1.E0)))
                    IF (LSTORE) THEN
                      CALL STCOOR (XT,YT,0)
                      CALL STCOOR (XTN,YTN,1)
                    END IF
                    GOTO 156
153               CONTINUE
                    CALL GRJMP(REAL(XT,KIND(1.E0)),
     .                         REAL(YT,KIND(1.E0)))
                    CALL GRDRW(REAL(XTN,KIND(1.E0)),
     .                         REAL(YTN,KIND(1.E0)))
                    IF (LSTORE) THEN
                      CALL STCOOR (XT,YT,0)
                      CALL STCOOR (XTN,YTN,1)
                    END IF
                    GOTO 156
154               CONTINUE
                    CALL GRJMP(REAL(XTN,KIND(1.E0)),
     .                         REAL(YTN,KIND(1.E0)))
                    IF (LSTORE) CALL STCOOR(XTN,YTN,0)
                    GOTO 156
155               CONTINUE
                    CALL GRJMP(REAL(XT,KIND(1.E0)),
     .                         REAL(YT,KIND(1.E0)))
                    CALL GRDRW(REAL(XT2,KIND(1.E0)),
     .                         REAL(YT2,KIND(1.E0)))
                    IF (LSTORE) THEN
                      CALL STCOOR (XT,YT,0)
                      CALL STCOOR (XT2,YT2,1)
                    END IF
                    GOTO 156
156             CONTINUE
C
                IF (IFL.GT.0.AND.K.EQ.ISWC(ISW)) THEN
                  IF (.NOT.NLSPLT(NU)) CALL GRNWPN (ICLR(ISW))
                  IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
                  IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
                  ISW=ISW+1
                  CALL GRJMP (REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
                  LSTORE = PLSTOR
                  INOSF=0
                  IF (LSTORE) CALL STCOOR (XTN,YTN,0)
                ENDIF
C
C
150         CONTINUE
            CALL GRNWPN(1)
158       CONTINUE
          CALL GRDSH(1.,0.,1.)

c  fill in isolated cells (defined in couple_... or any other user
c  segment): all cells j with nstgrd(j)=1
          do j=1,nsurf
            if (nstgrd(j).eq.1) then
              call ncelln(j,ir,ip,it,ia,ib,nr1st,np2nd,nt3rd,nbmlt,
     .                                     nlrad,nlpol,nltor)
              XPS(1)=XPOL(IR,IP)
              YPS(1)=YPOL(IR,IP)
              XPS(2)=XPOL(IR,IP+1)
              YPS(2)=YPOL(IR,IP+1)
              XPS(3)=XPOL(IR+1,IP+1)
              YPS(3)=YPOL(IR+1,IP+1)
              XPS(4)=XPOL(IR+1,IP)
              YPS(4)=YPOL(IR+1,IP)
              XPS(5)=XPOL(IR,IP)
              YPS(5)=YPOL(IR,IP)
              IF (NLTRA) XPS(1:5)=XPS(1:5)+RMTOR
              CALL GRFILL(5,XPS,YPS,1,1)
            ENDIF
          enddo
          IF (PLARR) THEN
            CALL GRNWPN(2)
            do ir=1,nr1st
              do ip=1,np2nd-1
                if (inmp1i(ir,ip,0).ne.0) then
                  ists=inmp1i(ir,ip,0)
                  xtn=0.5d0*(xpol(ir,ip)+xpol(ir,ip+1))
                  IF (NLTRA) XTN=XTN+RMTOR
                  ytn=0.5d0*(ypol(ir,ip)+ypol(ir,ip+1))
                  xtip=xtn+10.*plnx(ir,ip)*isign(1,ilside(nlim+ists))
                  ytip=ytn+10.*plny(ir,ip)*isign(1,ilside(nlim+ists))
                  CALL TSTCHM(1,XTN,YTN,XTip,YTip,IN,TESTN,
     .                      XMI2D,XMA2D,YMI2D,YMA2D,XT2,YT2)
                  if (testn .ne. 2)
     .            call grarrw (REAL(xtn,KIND(1.E0)),REAL(ytn,KIND(1.E0))
     .                        ,REAL(xtip,KIND(1.E0)),
     .                         REAL(ytip,KIND(1.E0)),0.4,0.4,0)
                endif
              enddo
            enddo
            CALL GRNWPN(1)
          ENDIF
C
        ELSEIF (PLCUT(2)) THEN
C
C  X-Z PLOT
          Y=CH2Z0
          IF (NLTOR) THEN
            IF (NLTRZ) THEN
              YW1=ZSURF(NPLINT)
              YW2=ZSURF(NPLOTT)
            ELSEIF (NLTRA) THEN
              YW1=ZSURF(NPLINT)
              YW2=ZSURF(NPLOTT)
            ELSEIF (NLTRT) THEN
              GOTO 990
            ENDIF
          ELSEIF (.NOT.NLTOR) THEN
            IF (NLTRZ) THEN
              YW1=ZIA
              YW2=ZAA
            ELSEIF (NLTRA) THEN
              YW1=0.
              YW2=PI2A
            ELSEIF (NLTRT) THEN
              GOTO 990
            ENDIF
          ENDIF
C
          DO 1544 NU=NPLINR,NPLOTR,NPLDLR
            LSTORE = .FALSE.
            IF (NLSPLT(NU)) CALL GRNWPN(2)
C  FIND IY SUCH THAT YPOL(NU,IY-1) < Y < YPOL(NU,IY), IF POSSIBLE
            DO 1549 IY=2,NP2ND
              IP=0
              IF (YPOL(NU,IY-1).LE.YPOL(NU,IY)) THEN
                IF (YPOL(NU,IY-1).LE.Y.AND.Y.LE.YPOL(NU,IY)) IP=IY-1
              ELSE
                IF (YPOL(NU,IY).LE.Y.AND.Y.LE.YPOL(NU,IY-1)) IP=IY
              ENDIF
              IF (IP.EQ.0) GOTO 1549
C  IP FOUND. NOW CHECK, IF THAT WAS A NON DEFAULT SURFACE SEGMENT
              DO 1545 J=1,NSTSI
                IF (NU.EQ.INUMP(J,1).AND.IP.GE.IRPTA(J,2).AND.
     .                                   IP.LE.IRPTE(J,2)) THEN
                  IF (ILCOL(NLIM+J) == 666) GOTO 1549
                  CALL GRNWPN(ILCOL(NLIM+J))
                  LSTORE = PLSTOR
                  INOSF=NLIM+J
                  IF (ILIIN(NLIM+J).LE.0) GOTO 1546
                  CALL GRDSH(1.,0.,1.)
                  GOTO 1547
                ENDIF
1545          CONTINUE
1546          CALL GRDSH(0.2,0.5,0.2)
1547          CONTINUE
C  NOW PLOT THAT RADIAL SURFACE
              XW1=XPOL(NU,IP)
              IF (NLTRZ) THEN
                IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
                  CALL GRJMP (REAL(XW1,KIND(1.E0)),REAL(YW1,KIND(1.E0)))
                  CALL GRDRW (REAL(XW1,KIND(1.E0)),REAL(YW2,KIND(1.E0)))
                  IF (LSTORE) THEN
                    CALL STCOOR (XW1,YW1,0)
                    CALL STCOOR (XW1,YW2,1)
                  END IF
                ENDIF
              ELSEIF (NLTRA) THEN
                X=RMTOR+XW1
                Z=TANAL*X
                RR=SQRT(X*X+Z*Z)
                IF (NTTRA.GT.100) THEN
                  WRITE (6,*) 'ERROR IN PLT2D '
                  CALL EXIT_OWN(1)
                ENDIF
                DO J=1,NTTRA
                  XX(J)=RR*COS(ZSURF(J))
                  YY(J)=RR*SIN(ZSURF(J))
                ENDDO
                CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
C             ELSEIF (NLTRT) THEN
              ENDIF
1549        CONTINUE
            CALL GRNWPN(1)
1544      CONTINUE
          CALL GRDSH(1.,0.,1.)
C
        ENDIF
C
      ELSEIF (LEVGEO.EQ.4.AND.PLCUT(3)) THEN
C
        INSTOR=0
        DO 1010,I=1,NTRII
          P=0.
          DO 1020,J=1,3
            IECKE2 = J+1
            IF (IECKE2 .EQ. 4) IECKE2 = 1
C           SEITE J GEHOERT ZUM RAND
            ISTS=ABS(INMTI(J,I))
            IF (ISTS .GT. 0 .AND.
     .          ISTS .LE. NLIM+NSTSI) THEN
C   SEITE J HAT BESONDERE EIGENSHAFTEN (NON-DEFAULT STD.FLAECHE)
              IF (ILCOL(ISTS) == 666) CYCLE
              CALL GRDSH (1.,0.,1.)
              CALL GRNWPN (ILCOL(ISTS))
              LSTORE = PLSTOR
              INOSF=ISTS
            ELSE
              CALL GRDSH (0.2,0.5,0.2)
              CALL GRNWPN (1)
              LSTORE = .FALSE.
              INOSF=0
            ENDIF
            XX(1)=XTRIAN(NECKE(J,I))
            YY(1)=YTRIAN(NECKE(J,I))
            XX(2)=XTRIAN(NECKE(IECKE2,I))
            YY(2)=YTRIAN(NECKE(IECKE2,I))
            IF (NLTRA) THEN
              XX(1)=XX(1)+RMTOR
              XX(2)=XX(2)+RMTOR
            ENDIF
            dsd(j)=((xx(1)-xx(2))**2+(yy(1)-yy(2))**2)**0.5
            p=p+dsd(j)
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
1020      CONTINUE
          CALL GRNWPN(1)
          IF (PLNUMV) THEN
C  Write number of triangle
c  p: halber umfang des dreiecks
c  r: radius des inkreises des dreiecks
            p=p*0.5
            r=((p-dsd(1))*(p-dsd(2))*(p-dsd(3))/p)**0.5
            WRITE (CH,'(I4)') I
            slt=r/2.
            X=XCOM(I)
            IF (I.GE.1000) X=X-SLT
            IF (I.GE.100) X=X-SLT
            IF (I.GE.10) X=X-SLT
            Y=YCOM(I)
            if (x.lt.xmi2d.or.x.gt.xma2d.or.y.lt.ymi2d.or.y.gt.yma2d)
     .         goto 1010
            slt=slt*sclfcx
            call grchrc(REAL(slt,KIND(1.E0)),0.,idummy)
            CALL GRTXT (REAL(X,KIND(1.E0)),REAL(Y,KIND(1.E0)),4,CH)
            call grchrc(0.3,0.,idummy)
          ENDIF
1010    CONTINUE
        CALL GRDSH (1.,0.,1.)
        CALL GRNWPN (1)
C
      ELSEIF (LEVGEO.EQ.5) THEN
        CALL GRCLP(1)
        IF (PLCUT(1)) THEN
          EBENE = (/-CH2Z0,1._DP,0._DP,0._DP/)
          IC1=2
          IC2=3
        ELSEIF (PLCUT(2)) THEN
          EBENE = (/-CH2Z0,0._DP,1._DP,0._DP/)
          IC1=1
          IC2=3
        ELSEIF (PLCUT(3)) THEN
          EBENE = (/-CH2Z0,0._DP,0._DP,1._DP/)
          IC1=1
          IC2=2
        END IF
        NU_LOOP:DO NU=NPLINR,NPLOTR,NPLDLR
          IF (VOL(NU) < EPS30) CYCLE NU_LOOP
          DO J=1,NSTSI
            IF (NU.EQ.INUMP(J,1)) THEN
              IF (ILCOL(NLIM+J) == 666) CYCLE NU_LOOP
              CALL GRNWPN(ILCOL(NLIM+J))
              IF (ILIIN(NLIM+J).LE.0) GOTO 1646
              CALL GRDSH(1.,0.,1.)
              GOTO 1647
            ENDIF
          ENDDO
1646      CALL GRDSH(0.2,0.5,0.2)
1647      CONTINUE
          TET(1,1:3) = (/ XTETRA(NTECK(1,NU)),
     .                    YTETRA(NTECK(1,NU)),
     .                    ZTETRA(NTECK(1,NU)) /)
          TET(2,1:3) = (/ XTETRA(NTECK(2,NU)),
     .                    YTETRA(NTECK(2,NU)),
     .                    ZTETRA(NTECK(2,NU)) /)
          TET(3,1:3) = (/ XTETRA(NTECK(3,NU)),
     .                    YTETRA(NTECK(3,NU)),
     .                    ZTETRA(NTECK(3,NU)) /)
          TET(4,1:3) = (/ XTETRA(NTECK(4,NU)),
     .                    YTETRA(NTECK(4,NU)),
     .                    ZTETRA(NTECK(4,NU)) /)
          CALL TETRA_SCHNITT (EBENE, TET, NCTPNT, CTPNTS)
          IF (NCTPNT > 0) THEN
            CALL GRJMP (REAL(CTPNTS(1,IC1),KIND(1.E0)),
     .                  REAL(CTPNTS(1,IC2),KIND(1.E0)))
            DO ICT = 2, NCTPNT
              CALL GRDRW (REAL(CTPNTS(ICT,IC1),KIND(1.E0)),
     .                    REAL(CTPNTS(ICT,IC2),KIND(1.E0)))
            END DO
            CALL GRDRW (REAL(CTPNTS(1,IC1),KIND(1.E0)),
     .                  REAL(CTPNTS(1,IC2),KIND(1.E0)))
          END IF
          CALL GRNWPN(1)
        ENDDO NU_LOOP
        CALL GRDSH(1.,0.,1.)
        CALL GRCLP(0)
C
      ELSE
C
C  OPTION (PL1ST.AND.LEVGEO.EQ.4.AND.PLCUT(2)) IS STILL MISSING
C
        GOTO 990
C
      ENDIF
C
170   IF (.NOT.PL2ND.OR..NOT.NLPOL.OR.PLCUT(2))   GOTO 200
C
C  PLOT POLOIDAL GRID
C
      IF (LEVGEO.EQ.1) THEN
C
        CALL GRDSH(0.2,0.5,0.2)
        IF (PLCUT(3)) THEN
C X-Y-PLANE
          ALLOCATE (IFARB(N1ST,N2ND))
          ALLOCATE (IDASH(N1ST,N2ND))
          IFARB = 1
          IDASH = 0
          DO J=1,NSTSI
            NU = INUMP(J,2)
            IF (NU <= 0) CYCLE 
            IRA = IRPTA(J,1)
            IRE = IRPTE(J,1)-1
            IFARB(IRA:IRE,NU) = ILCOL(NLIM+J)
            IDASH(IRA:IRE,NU) = ILIIN(NLIM+J)
            IF (NLSPLT(N1ST+NU)) IFARB(NU,:) = 2
          END DO
          XW1=RSURF(NPLINR)
          XW2=RSURF(NPLOTR)
          IF (NLTRA) XW1=XW1+RMTOR
          IF (NLTRA) XW2=XW2+RMTOR
        ELSEIF (PLCUT(1)) THEN
C Y-Z-PLANE
          ALLOCATE (IFARB(N2ND,N3RD))
          ALLOCATE (IDASH(N2ND,N3RD))
          IFARB = 1
          IDASH = 0
          DO J=1,NSTSI
            NU = INUMP(J,2)
            IF (NU <= 0) CYCLE 
            ITA = IRPTA(J,3)
            ITE = IRPTE(J,3)-1
            IFARB(NU,ITA:ITE) = ILCOL(NLIM+J)
            IDASH(NU,ITA:ITE) = ILIIN(NLIM+J)
            IF (NLSPLT(N1ST+NU)) IFARB(NU,:) = 2
          END DO
          IF (NLTOR) THEN
            IF (NLTRZ) THEN
              ZW1=ZSURF(NPLINT)
              ZW2=ZSURF(NPLOTT)
            ELSEIF (NLTRA) THEN
              ZW1=ZSURF(NPLINT)
              ZW2=ZSURF(NPLOTT)
            ELSEIF (NLTRT) THEN
              GOTO 990
            ENDIF
          ELSEIF (.NOT.NLTOR) THEN
            IF (NLTRZ) THEN
              ZW1=ZIA
              ZW2=ZAA
            ELSEIF (NLTRA) THEN
              ZW1=0.
              ZW2=PI2A
            ELSEIF (NLTRT) THEN
              GOTO 990
            ENDIF
          ENDIF
        ENDIF
        DO NU=NPLINP,NPLOTP,NPLDLP
          LSTORE = .FALSE.
          YW1=PSURF(NU)
          IF (PLCUT(3)) THEN
            DO IR=NPLINR, MIN(NPLOTR-1,NR1STM), NPLDLR
              IF (IFARB(IR,NU) == 666) CYCLE
              XX(1)=RSURF(IR)
              XX(2)=RSURF(IR+1)
              YY(1)=YW1
              YY(2)=YW1
              CALL GRNWPN(IFARB(IR,NU))
              IF (IDASH(IR,NU) <= 0) THEN
                CALL GRDSH(0.2,0.5,0.2)
              ELSE
                CALL GRDSH(1.,0.,1.)
              END IF
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
            END DO
          ELSEIF (PLCUT(1)) THEN
            DO IT=NPLINT, MIN(NPLOTT-1,NT3RDM), NPLDLT
              IF (IFARB(NU,IT) == 666) CYCLE
              XX(1)=XW1
              XX(2)=XW1
              YY(1)=ZSURF(IT)
              YY(2)=ZSURF(IT+1)
              CALL GRNWPN(IFARB(NU,IT))
              IF (IDASH(NU,IT) <= 0) THEN
                CALL GRDSH(0.2,0.5,0.2)
              ELSE
                CALL GRDSH(1.,0.,1.)
              END IF
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
            END DO
          ENDIF
        END DO

        IF (ALLOCATED(IFARB)) THEN
          DEALLOCATE (IFARB)
          DEALLOCATE (IDASH)
        END IF

!        DO 1171 NU=NPLINP,NPLOTP,NPLDLP
!          LSTORE = .FALSE.
!          DO 1172 J=1,NSTSI
!            IF (NU.EQ.INUMP(J,2)) THEN
!              IF (ILCOL(NLIM+J) == 666) GOTO 1171
!              CALL GRNWPN(ILCOL(NLIM+J))
!              LSTORE = PLSTOR
!              INOSF=NLIM+J
!              IF (ILIIN(NLIM+J).LE.0) GOTO 1173
!              CALL GRDSH(1.,0.,1.)
!              GOTO 1174
!            ENDIF
!1172      CONTINUE
!1173      CALL GRDSH(0.2,0.5,0.2)
!1174      CONTINUE
!          IF (NLSPLT(N1ST+NU)) CALL GRNWPN(2)
!          YW1=PSURF(NU)
!          IF (PLCUT(3)) THEN
!            XX(1)=XW1
!            XX(2)=XW2
!            YY(1)=YW1
!            YY(2)=YW1
!            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
!          ELSEIF (PLCUT(1)) THEN
!            XX(1)=ZW1
!            XX(2)=ZW2
!            YY(1)=YW1
!            YY(2)=YW1
!            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
!          ENDIF
!          CALL GRNWPN(1)
!1171    CONTINUE
        CALL GRDSH(1.,0.,1.)
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.PLCUT(3)) THEN
C
        DO 176 NU=NPLINP,NPLOTP,NPLDLP
          IF (NLSPLT(N1ST+NU)) CALL GRNWPN(2)
          IAN=MAX(1,NPLINR)
          IEN=MIN(NR1ST,NPLOTR)
          NSW=0
          DO 177 J=1,NSTSI
C
            IF (NU.EQ.INUMP(J,2)) THEN
              IF (ILCOL(NLIM+J) == 666) CYCLE
              NSW=NSW+1
              ICLR(NSW)=ILCOL(NLIM+J)
              IDSH(NSW)=1
              IF (ILIIN(NLIM+J).GT.0) IDSH(NSW)=0
              ISWC(NSW)=MAX(IAN,IRPTA(J,1))
              ICLR(NSW+1)=1
              IDSH(NSW+1)=1
              ISWC(NSW+1)=MIN(IEN,IRPTE(J,1))
              NSW=NSW+1
            ENDIF
177       CONTINUE
C
C  SORTIEREN, NACH "RADIALEM BEGIN" < "RADIALEM ENDE"
          IF (NSW.EQ.0) GOTO 1797
179       ISW=0
          DO 178 I=1,NSW-3,2
            IF (ISWC(I).GT.ISWC(I+2)) THEN
              ISW=ISW+1

              IHELP=ISWC(I)
              ISWC(I)=ISWC(I+2)
              ISWC(I+2)=IHELP
              IHELP=ISWC(I+1)
              ISWC(I+1)=ISWC(I+3)
              ISWC(I+3)=IHELP

              IHELP=ICLR(I)
              ICLR(I)=ICLR(I+2)
              ICLR(I+2)=IHELP
              IHELP=ICLR(I+1)
              ICLR(I+1)=ICLR(I+3)
              ICLR(I+3)=IHELP

              IHELP=IDSH(I)
              IDSH(I)=IDSH(I+2)
              IDSH(I+2)=IHELP
              IHELP=IDSH(I+1)
              IDSH(I+1)=IDSH(I+3)
              IDSH(I+3)=IHELP
            ENDIF
178       CONTINUE
          IF (ISW.GT.0) GOTO 179
C
          I=0
1781      I=I+1
1785      IF (ISWC(I).EQ.ISWC(I+1)) THEN
            ICLR(I)=ICLR(I+1)
            IDSH(I)=IDSH(I+1)
            DO 1796 J=I+2,NSW
              ISWC(J-1)=ISWC(J)
              ICLR(J-1)=ICLR(J)
              IDSH(J-1)=IDSH(J)
1796        CONTINUE
            NSW=NSW-1
            IF (I.LT.NSW) GOTO 1785
          ENDIF
          IF (I.LT.NSW-1) GOTO 1781
C
1797      CONTINUE
          NSW=NSW+1
          ISWC(NSW)=100000
          ICLR(NSW)=1
          IDSH(NSW)=1
          ISW=1
          CALL GRDSH (0.2,0.5,0.2)
          LSTORE = .FALSE.
          DO 1751 I=IAN,IEN
            IFL=I-IAN
            XTN=XPOL(I,NU)
            IF (NLTRA) XTN=XTN+RMTOR
            YTN=YPOL(I,NU)
            CALL TSTCHM(IFL,XTN,YTN,XT,YT,IN,TESTN,
     .                  XMI2D,XMA2D,YMI2D,YMA2D,XT2,YT2)
C
            IF (IFL.EQ.0.AND.I.EQ.ISWC(ISW)) THEN
              IF (.NOT.NLSPLT(N1ST+NU)) CALL GRNWPN (ICLR(ISW))
              LSTORE = PLSTOR
              INOSF=0
              IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
              IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
              ISW=ISW+1
            ENDIF
C
            IF (IN.EQ.0) IN=4
            GOTO (171,172,173,174,1752),IN
171           CONTINUE
                CALL GRDRW (REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
                IF (LSTORE) CALL STCOOR(XTN,YTN,1)
                GOTO 175
172           CONTINUE
                CALL GRDRW(REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
                CALL GRJMP(REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
                IF (LSTORE) THEN
                  CALL STCOOR (XT,YT,0)
                  CALL STCOOR (XTN,YTN,1)
                END IF
                GOTO 175
173           CONTINUE
                CALL GRJMP(REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
                CALL GRDRW(REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
                IF (LSTORE) THEN
                  CALL STCOOR (XT,YT,0)
                  CALL STCOOR (XTN,YTN,1)
                END IF
                GOTO 175
174           CONTINUE
                CALL GRJMP(REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
                IF (LSTORE) CALL STCOOR(XTN,YTN,0)
                GOTO 175
1752          CONTINUE
                CALL GRJMP(REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
                CALL GRDRW(REAL(XT2,KIND(1.E0)),REAL(YT2,KIND(1.E0)))
                IF (LSTORE) THEN
                  CALL STCOOR (XT,YT,0)
                  CALL STCOOR (XT2,YT2,1)
                END IF
                GOTO 175
175       CONTINUE
C
          IF (IFL.GT.0.AND.I.EQ.ISWC(ISW)) THEN
            IF (.NOT.NLSPLT(N1ST+NU)) CALL GRNWPN (ICLR(ISW))
            IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
            IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
            ISW=ISW+1
            LSTORE = PLSTOR
            INOSF=0
            CALL GRJMP (REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
            IF (LSTORE) CALL STCOOR(XTN,YTN,0)
          ENDIF
C
1751      CONTINUE
          CALL GRNWPN(1)
176     CONTINUE
        CALL GRDSH(1.,0.,1.)
C
        if (PLARR) then
          call grnwpn(2)
          do ip=1,np2nd
            do ir=1,nr1st-1
              if (inmp2i(ir,ip,0).ne.0) then
                ists=inmp2i(ir,ip,0)
                xtn=0.5d0*(xpol(ir,ip)+xpol(ir+1,ip))
                IF (NLTRA) XTN=XTN+RMTOR
                ytn=0.5d0*(ypol(ir,ip)+ypol(ir+1,ip))
                xtip=xtn+10.*pplnx(ir,ip)*isign(1,ilside(nlim+ists))
                ytip=ytn+10.*pplny(ir,ip)*isign(1,ilside(nlim+ists))
                CALL TSTCHM(1,XTN,YTN,XTip,YTip,IN,TESTN,
     .                      XMI2D,XMA2D,YMI2D,YMA2D,XT2,YT2)
                if (testn .ne. 2)
     .          call grarrw (REAL(xtn,KIND(1.E0)),REAL(ytn,KIND(1.E0)),
     .                       REAL(xtip,KIND(1.E0)),
     .                       REAL(ytip,KIND(1.E0)),
     .                       0.4,0.4,0)
              endif
            enddo
          enddo
          call grnwpn(1)
        endif
C
      ELSE
C
C  OPTION (PL2ND.AND.(LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.PLCUT(1)) IS STILL MISSING
        GOTO 991
C
      ENDIF
C
C  PLOT 3RD GRID (Z OR TOROIDAL)
C
200   IF (.NOT.PL3RD.OR..NOT.NLTOR.OR.PLCUT(3))  GOTO 220
C
      IF (NLTRZ.AND.LEVGEO.EQ.1) THEN
C Y-Z-PLANE
        IF (PLCUT(1)) THEN
          XW1=YIA
          XW2=YAA
          IF (NLPOL)THEN
            XW1=PSURF(NPLINP)
            XW2=PSURF(NPLOTP)
          ENDIF
C X-Z-PLANE
        ELSEIF (PLCUT(2)) THEN
          XW1=RIA
          XW2=RAA
          IF (NLRAD) THEN
            XW1=RSURF(NPLINR)
            XW2=RSURF(NPLOTR)
          ENDIF
        ENDIF
        DO 204 NU=NPLINT,NPLOTT,NPLDLT
          LSTORE = .FALSE.
          DO 205 J=1,NSTSI
            IF (NU.EQ.INUMP(J,3)) THEN
              IF (ILCOL(NLIM+J) == 666) GOTO 204
              CALL GRNWPN(ILCOL(NLIM+J))
              LSTORE = PLSTOR
              INOSF=NLIM+J
              IF (ILIIN(NLIM+J).LE.0) GOTO 206
              CALL GRDSH(1.,0.,1.)
              GOTO 207
            ENDIF
205       CONTINUE
206       CALL GRDSH(0.2,0.5,0.2)
207       CONTINUE
          IF (NLSPLT(N1ST+N2ND+NU)) CALL GRNWPN(2)
          ZW1=ZSURF(NU)
          IF (PLCUT(1)) THEN
            XX(1)=ZW1
            XX(2)=ZW1
            YY(1)=XW1
            YY(2)=XW2
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
          ELSEIF (PLCUT(2)) THEN
            XX(1)=XW1
            XX(2)=XW2
            YY(1)=ZW1
            YY(2)=ZW1
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)
          ENDIF
          CALL GRNWPN(1)
204     CONTINUE
        CALL GRDSH(1.,0.,1.)
C
      ELSE
C
        GOTO 992
C
      ENDIF
C
C
220   CONTINUE
C
C   PLOT  ADDITIONAL SURFACES
C
      IF (.NOT.PLADD) GOTO 250
C
        LZR=.TRUE.
C
        IF (NLTRA) THEN
          CALL XSHADD(RMTOR,1,NLIMI)
          DO I=1,NLIMI
            IF (ILCOL(I) == 666) CYCLE
            IF (ILTOR(I).GT.0) THEN
              CALL SETROT (AFF,AFFI,1,0.0_DP,1.0_DP,0.0_DP,
     .                    -ZZONE(ILTOR(I))*180._DP/PIA)
              CALL ROTADD (AFF,AFFI,I,I)
            ENDIF
C
            CALL PLTADD(I,I)
C
            IF (ILTOR(I).GT.0) THEN
              CALL SETROT (AFF,AFFI,1,0.0_DP,1.0_DP,0.0_DP,
     .                     ZZONE(ILTOR(I))*180._DP/PIA)
              CALL ROTADD (AFF,AFFI,I,I)
            ENDIF
          ENDDO
          CALL XSHADD(-RMTOR,1,NLIMI)
        ELSEIF (.NOT.NLTRA) THEN
          CALL PLTADD(1,NLIMI)
        ENDIF
C
250   CONTINUE
C
300   CONTINUE
C
      IF (.NOT.NLPL3D) GOTO 400
C
      PLSAV1=PLNUMS
      PLSAV2=PLSTOR
      PLNUMS=.FALSE.
      PLSTOR=.FALSE.
      CALL PLT3D (XN2D,YN2D,FX2D,FY2D,ITH,XMI2D,XMA2D,YMI2D,YMA2D)
      PLNUMS=PLSAV1
      PLSTOR=PLSAV2
C  IF NLTRA: 3D PLOTTING IS DONE IN THE LOCAL CO-ORD. SYSTEM NO. ITHPL
      ITHPL=ITH
C
400   CONTINUE
      RETURN
C
C  PLOT PARTICLE HISTORIES IN GEOMETRY-PLOT
C
      ENTRY CHCTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)
c slmod begin
      IF (.TRUE.) THEN
        CALL USRTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)
      ENDIF 
c slmod end
C
C  IFLAG.EQ.0 ONLY SYMBOL AT XPLO,YPLO,ZPLO
C  IFLAG.NE.0 TRACK FROM LAST POSITION (PREVIOUS CALL) TO XPLO,YPLO,ZPLO
C  ISYM       NUMBER OF SYMBOL FOR THE CURRENT EVENT
C
C  TEXT FOR PARTICLE-HISTORIES PLOT, ONLY AT FIRST CALL TO THIS ENTRY
C
      XN=XN2D
      YN=YN2D
      FX=FX2D
      FY=FY2D
      IF (IWRIT.EQ.0.AND.PLHST) THEN
        IWRIT=1
        IF (.NOT.NLPL3D) CALL GRSCLV(0.,0.,REAL(ABSMAX,KIND(1.E0)),
     .                                     REAL(ORDMAX,KIND(1.E0)))
        XNP05=XN+0.5/FX
        DO IA=1,NTXHST
          YYIA=YN-(0.75*(IA-1))/FY
          CALL GRJMPS (REAL(XN,KIND(1.E0)),REAL(YYIA+0.15,KIND(1.E0)),
     .                 ISPL(IA))
          CALL GRTXT (REAL(XNP05,KIND(1.E0)),REAL(YYIA,KIND(1.E0)),20,
     .                TXTHST(IA))
        ENDDO
C
        IF (NLPL3D) THEN
          XN1=XN+0./FX
          XN2=XN1+2./FX
          XN3=XN2+1.5/FX
          YNP=YYIA-1.5/FY
        ELSE
C  FX=FY=1.
          XN1=ABSMAX+4./FX
          XN2=XN1+2./FX
          XN3=XN2+1.5/FX
          YNP=ORDMAX+0.65/FY
        ENDIF
        IC=1
        DO 510 I=1,NATMI
          ISP=NSPH+I
          IF (NHSTS(ISP).EQ.-1) GOTO 510
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
510     CONTINUE
        DO 511 I=1,NMOLI
          ISP=NSPA+I
          IF (NHSTS(ISP).EQ.-1) GOTO 511
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
511     CONTINUE
        DO 512 I=1,NIONI
          ISP=NSPAM+I
          IF (NHSTS(ISP).EQ.-1) GOTO 512
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
512     CONTINUE
        DO 513 I=1,NPLSI
          ISP=NSPAMI+I
          IF (NHSTS(ISP).EQ.-1) GOTO 513
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
513     CONTINUE
        DO 514 I=1,NPHOTI
          ISP=I
          IF (NHSTS(ISP).EQ.-1) GOTO 514
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
514     CONTINUE
        IF (.NOT.NLPL3D)
     .     CALL GRSCLV(REAL(XMI2D,KIND(1.E0)),REAL(YMI2D,KIND(1.E0)),
     .                 REAL(XMA2D,KIND(1.E0)),REAL(YMA2D,KIND(1.E0)))
      ENDIF
C
C  WRITE TRACK DATA
C
      IF (TRCHST) THEN
        CALL LEER(1)
        WRITE (6,*) TXTHST(ISYM)
        IF (ISYM.EQ.1.OR..NOT.NLTRC)  CALL MASJ1('NPANU   ',NPANU)
        IF (ISPZ.GT.0.AND.ISPZ.LE.NSPZ) THEN
          WRITE (6,'(1X,A8)') TEXTS(ISPZ)
        ELSE
          WRITE (6,'(1X,A13)') 'LINE OF SIGHT'
        ENDIF
        CALL MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
        CALL MASR6 ('VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .               VELX,VELY,VELZ,VEL,E0,WEIGHT)
        CALL MASR1 ('TIME    ',TIME)
        CALL MASR1 ('XGENER  ',XGENER)
        CALL MASJ3 ('NRCELL,NACELL,NBLOCK    ',NRCELL,NACELL,NBLOCK)
        IF (NLPLG) CALL MASJ1('IPOLG   ',IPOLG)
        IF (NLTRA) THEN
          CALL MASJ1('IPERID  ',IPERID)
          CALL MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI/DEGRAD)
        ENDIF
        IF (NLPOL) THEN
          CALL MASJ1 ('NPCELL  ',NPCELL)
        ENDIF
        IF (NLTOR) THEN
          CALL MASJ1 ('NTCELL  ',NTCELL)
        ENDIF
        IF (NLTET.AND.(IPOLGN>0).AND.(MRSURF>0)) THEN
          CALL MASJ1 ('INMTIT  ',INMTIT(IPOLGN,MRSURF))
        END IF
C  CALLED FROM DIAGNO, WITH ISTRA=0?
        IF (ISTRA.EQ.0) THEN
          ISTR=1
        ELSE
          ISTR=ISTRA
        ENDIF
        IF ((ISYM.GE.6.AND.ISYM.LE.10).OR.
     .      (ISYM.EQ.1.AND.NLSRF(ISTR))) THEN
          CALL MASJ4 ('MRSURF,MPSURF,MTSURF,MASURF     ',
     .                 MRSURF,MPSURF,MTSURF,MASURF)
          CALL MASR1 ('SCOS    ',SCOS)
        ENDIF
      ENDIF
C
C  PLOT HISTORY
C
      IF (.NOT.PLHST) GOTO 410
C
C  PREPARE PLOT DATA
C
      IF (NLPL3D) THEN
C 3D PLOT: NEW CO-ORDINATES: XPL,YPL,ZPL, PROJECTED TO XWN,YWN
        IF (NLTRA) THEN
          IF (.NOT.NLTOR) THEN
            Z1=ZSURF(1)
            PPHI=MOD(PHI+PI2A-Z1,PI2A)+Z1
            NTT=LEARCA(PPHI,ZSURF,1,NTTRA,1,'CHCTRC ')
            CALL FZRTOR(XPLO,ZPLO,NTT,XR,THET,NTDUM,.FALSE.,0)
          ELSEIF (NLTOR) THEN
            CALL FZRTOR(XPLO,ZPLO,NTCELL,XR,THET,NTDUM,.FALSE.,0)
          ENDIF
          CALL FZRTRI(XPL,ZPL,ITHPL,XR,THET,NTDUM)
          YPL=YPLO
        ELSEIF (NLTRZ) THEN
          XPL=XPLO
          YPL=YPLO
          ZPL=ZPLO
        ENDIF
        CALL PL3D(XPL,YPL,ZPL,XWN,YWN)
      ELSEIF (NLPL2D) THEN
C 2D PLOT: NEW CO-ORDINATES: XWN,YWN
        IF (PLCUT(3)) THEN
          XWN=XPLO
          YWN=YPLO
          IF (NLTRA) THEN
            XWN=XPLO+RMTOR
          ENDIF
        ELSEIF (PLCUT(2)) THEN
          IF (NLTRA) THEN
            XWN=XPLO
            CALL FZRTRA(XWN,ZWN,PHI,NT)
            XWN=XWN+RMTOR
            RWN=SQRT(XWN**2+ZWN**2)
            XWN=RWN*COS(PHI)
            YWN=RWN*SIN(PHI)
          ELSE
            XWN=XPLO
            YWN=ZPLO
          ENDIF
        ELSEIF (PLCUT(1)) THEN
          XWN=ZPLO
          YWN=YPLO
        ENDIF
      ENDIF
C
C  NOW PLOT TRACK FROM XWO,YWO  TO  XWN,YWN
C
      ICOLOR=MOD(ISPZ+1,7)
      IF (ICOLOR.EQ.0) ICOLOR=7
      CALL GRNWPN(ICOLOR)
      IF (ICOL.EQ.1) CALL GRDSH(0.2,0.5,0.2)
C
      CALL TSTCHM(IFLAG,XWN,YWN,XT,YT,IN,TESTN,XMI2D,XMA2D,YMI2D,YMA2D,
     .            XT2,YT2)
      IF (IN.EQ.0) IN=4
      IF (IN.LE.2) 
     .   CALL GRJMP (REAL(XWO,KIND(1.E0)),REAL(YWO,KIND(1.E0)))
      IF (ILINIE.EQ.0.OR.NHSTS(ISPZ).EQ.-1) GOTO 404
      GOTO (401,402,403,404,405),IN
401     CALL GRDRW (REAL(XWN,KIND(1.E0)),REAL(YWN,KIND(1.E0)))
        GOTO 406
402     CALL GRDRW (REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
        CALL GRJMP (REAL(XWN,KIND(1.E0)),REAL(YWN,KIND(1.E0)))
        GOTO 406
403     CALL GRJMP (REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
        CALL GRDRW (REAL(XWN,KIND(1.E0)),REAL(YWN,KIND(1.E0)))
        GOTO 406
404     CALL GRJMP (REAL(XWN,KIND(1.E0)),REAL(YWN,KIND(1.E0)))
        GOTO 406
405     CALL GRJMP (REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
        CALL GRDRW (REAL(XT2,KIND(1.E0)),REAL(YT2,KIND(1.E0)))
        GOTO 406
406   CONTINUE
C
C  PLOT SYMBOL
C
      IF (TESTN.EQ.0..AND.NHSTS(ISPZ).NE.-1) THEN
        DO 408 J=1,8
          IF (ISYM.EQ.ISYPLT(J).OR.ISYM.EQ.ISYM_ERR) THEN
            CALL GRJMPS(REAL(XWN,KIND(1.E0)),REAL(YWN,KIND(1.E0)),
     .                  ISPL(ISYM))
            GOTO 409
          ENDIF
408     CONTINUE
409     CONTINUE
      ENDIF
      IF (NLTRC.AND.TRCHST) THEN
        WRITE (6,*) 'ICOLOR,IN,ISYM,IFLAG,NHSTS ',
     .               ICOLOR,IN,ISYM,IFLAG,NHSTS(ISPZ)
        CALL LEER(1)
      ENDIF
C
      IF (ICOL.EQ.1) CALL GRDSH(1.,0.,1.)
      CALL GRNWPN(1)
C
410   CONTINUE
      XWO=XWN
      YWO=YWN
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'OPTION IN PLT2D NOT READY: X- OR RADIAL GRID'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (6,*) 'OPTION IN PLT2D NOT READY: Y- OR POLOIDAL GRID'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (6,*) 'OPTION IN PLT2D NOT READY: Z- OR TOROIDAL GRID'
      CALL EXIT_OWN(1)
      END
