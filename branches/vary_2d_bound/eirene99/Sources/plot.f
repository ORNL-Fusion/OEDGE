      SUBROUTINE PLT2D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CLGIN'
      INCLUDE 'CINIT'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CADGEO'
      INCLUDE 'CGEOM'
      INCLUDE 'CPLMSK'
      INCLUDE 'CPLOT'
      INCLUDE 'CPL3D'
      INCLUDE 'CGRID'
      INCLUDE 'COMUSR'
      INCLUDE 'CLOGAU'
      INCLUDE 'CTRCEI'
      INCLUDE 'CCONA'
      INCLUDE 'CPOLYG'
      INCLUDE 'COMSOU'
      INCLUDE 'CTEXT'
      INCLUDE 'CTRIG'
      INCLUDE 'CRECH'
      INCLUDE 'COMSPL'
c slmod begin
      DATA count /0.0/
c slmod end
C
      DIMENSION XX(101),YY(101)
      CHARACTER*10 CX,CY,CX0,CY0
      PARAMETER (NTXHST=15)
      CHARACTER*20 TXTHST(NTXHST)
      LOGICAL PLSAV1,PLSAV2
      REAL*4 XPS(5),YPS(5)
      DIMENSION ISPL(15),ICLR(2*NSTS+1),IDSH(2*NSTS+1),ISWC(2*NSTS+1)
      DIMENSION DSD(3),AFF(3,3),AFFI(3,3)
      CHARACTER*4 CH
      DIMENSION XCOMTR(NTRIS),YCOMTR(NTRIS)
      EQUIVALENCE (XCOM(1,1),XCOMTR(1)),
     .            (YCOM(1,1),YCOMTR(1))
      SAVE
      DATA ABSMAX,ORDMAX/21.,21./
      DATA XNULL,YNULL/9.,4./,XWN,YWN/0.,0./
      DATA IWRIT/0/,ISPL/2,3,4,5,6,14,15,16,17,0,0,0,7,8,48/
      DATA TXTHST
     .           /'LOCATE(1)           ',
     .            'ELECTR. IMPACT(2)   ',
     .            'ION IMPACT(3)       ',
     .            'CHARGE EXCHANGE(4)  ',
     .            'FOKKER PLANCK(5)    ',
     .            'SURFACE(6)          ',
     .            'SPLITTING(7)        ',
     .            'RUSSIAN ROULETTE(8) ',
     .            'PERIODICITY(9)      ',
     .            'RESTART:A. SPLT.(10)',
     .            'SAVE:COND. EXP.(11) ',
     .            'RESTART:COND EXP(12)',
     .            'TIME LIMIT(13)      ',
     .            'GENERATION LIMIT(14)',
     .            'ERROR DETECTED      '/
C
C  PREPARE PLOT OF STANDARD-MESH AND ADDITIONAL SURFACES
C
      IF (.NOT.NLPL2D) GOTO 300
C
      CALL FTCRE(CH2MX,CX)
      CALL FTCRE(CH2MY,CY)
      CALL FTCRE(CH2X0,CX0)
      CALL FTCRE(CH2Y0,CY0)
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
      CALL GRSCLV (0.,0.,sngl(ABSMAX),sngl(ORDMAX))
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
      XMI2D=CH2X0-CH2MX
      XMA2D=CH2X0+CH2MX
      YMI2D=CH2Y0-CH2MY
      YMA2D=CH2Y0+CH2MY
      CALL GRSCLV(sngl(XMI2D),sngl(YMI2D),sngl(XMA2D),sngl(YMA2D))
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
      IF (.NOT.PL1ST.OR.NR1ST.LE.1.OR.PLCUT(1)) GOTO 170
C
      IF (LEVGEO.EQ.1) THEN
C  X-Y-PLANE
        IF (PLCUT(3)) THEN
          YW1=YIA
          YW2=YAA
          IF (NLPOL) THEN
            YW1=PSURF(NPLINP)
            YW2=PSURF(NPLOTP)
          ENDIF
C  X-Z-PLANE
        ELSEIF (PLCUT(2)) THEN
          IF (NLTOR) THEN
            YW1=ZSURF(NPLINT)
            YW2=ZSURF(NPLOTT)
          ELSEIF (.NOT.NLTOR) THEN
            YW1=ZIA
            YW2=ZAA
          ENDIF
        ENDIF
        DO 144 NU=NPLINR,NPLOTR,NPLDLR
          DO 145 J=1,NSTSI
            IF (NU.EQ.INUMP(J,1)) THEN
              CALL GRNWPN(ILCOL(NLIM+J))
              IF (ILIIN(NLIM+J).LE.0) GOTO 146
              CALL GRDSH(1.,0.,1.)
              GOTO 147
            ENDIF
145       CONTINUE
146       CALL GRDSH(0.2,0.5,0.2)
147       CONTINUE
          IF (NLSPLT(NU)) CALL GRNWPN(2)
          XW1=RSURF(NU)
          IF (PLCUT(3)) THEN
            IF (NLTRA) XW1=XW1+RMTOR
            IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
              XX(1)=XW1
              XX(2)=XW1
              YY(1)=YW1
              YY(2)=YW2
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
            ENDIF
          ELSEIF (PLCUT(2)) THEN
            IF (NLTRZ) THEN
              IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
                XX(1)=XW1
                XX(2)=XW1
                YY(1)=YW1
                YY(2)=YW2
                CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
              ENDIF
            ELSEIF (NLTRA) THEN
              X=RMTOR+XW1
              Z=TANAL*X
              RR=SQRT(X*X+Z*Z)
              IF (NTTRA.GT.100) THEN
                WRITE (6,*) 'ERROR IN PLT2D '
                CALL EXIT
              ENDIF
              DO J=1,NTTRA+1
                XX(J)=RR*COS((J-1)*2.*ALPHA)
                YY(J)=RR*SIN((J-1)*2.*ALPHA)
              ENDDO
              CALL PLTLNE(NTTRA+1,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
            ENDIF
          ENDIF
          CALL GRNWPN(1)
144     CONTINUE
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
            DO 135 J=1,NSTSI
              IF (NU.EQ.INUMP(J,1)) THEN
                CALL GRNWPN(ILCOL(NLIM+J))
                IF (ILIIN(NLIM+J).LE.0) GOTO 136
                CALL GRDSH(1.,0.,1.)
                GOTO 137
              ENDIF
135         CONTINUE
136         CALL GRDSH(0.2,0.5,0.2)
137         CONTINUE
            IF (NLSPLT(NU)) CALL GRNWPN(2)
C  PLOT AT Y=CH2Z0 , TO BE WRITTEN. PRESENTLY AT Y=0
            XW1=RSURF(NU)+EP1(NU)
            XW2=-RSURF(NU)+EP1(NU)
            IF (NLTRZ) THEN
              XX(1)=XW1
              XX(2)=XW1
              YY(1)=YW1
              YY(2)=YW2
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
              XX(1)=XW2
              XX(2)=XW2
              CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
            ELSEIF (NLTRA) THEN
              X=RMTOR+XW1
              Z=TANAL*X
              RR=SQRT(X*X+Z*Z)
              IF (NTTRA.GT.100) THEN
                WRITE (6,*) 'ERROR IN PLT2D '
                CALL EXIT
              ENDIF
              DO J=1,NTTRA
                XX(J)=RR*COS(ZSURF(J))
                YY(J)=RR*SIN(ZSURF(J))
              ENDDO
              CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
              X=RMTOR+XW2
              Z=TANAL*X
              RR=SQRT(X*X+Z*Z)
              IF (NTTRA.GT.100) THEN
                WRITE (6,*) 'ERROR IN PLT2D '
                CALL EXIT
              ENDIF
              DO J=1,NTTRA
                XX(J)=RR*COS(ZSURF(J))
                YY(J)=RR*SIN(ZSURF(J))
              ENDDO
              CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
            ENDIF
            CALL GRNWPN(1)
134       CONTINUE
          CALL GRDSH(1.,0.,1.)
C
        ELSEIF (PLCUT(3)) THEN
C
C    X-Y PLOT
          DO 140 NU=NPLINR,NPLOTR,NPLDLR
            DO 141 J=1,NSTSI
              IF (NU.EQ.INUMP(J,1)) THEN
                CALL GRNWPN(ILCOL(NLIM+J))
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
            CALL PLTLNE(NRET,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
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
                  IF (.NOT.NLSPLT(NU)) CALL GRNWPN (ICLR(ISW))
                  IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
                  IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
                  ISW=ISW+1
                ENDIF
C
                IF (IN.EQ.0) IN=4
                GOTO (151,152,153,154,155),IN
151               CONTINUE
                    CALL GRDRW (SNGL(XTN),SNGL(YTN))
                    GOTO 156
152               CONTINUE
                    CALL GRDRW(SNGL(XT),SNGL(YT))
                    CALL GRJMP(SNGL(XTN),SNGL(YTN))
                    GOTO 156
153               CONTINUE
                    CALL GRJMP(SNGL(XT),SNGL(YT))
                    CALL GRDRW(SNGL(XTN),SNGL(YTN))
                    GOTO 156
154               CONTINUE
                    CALL GRJMP(SNGL(XTN),SNGL(YTN))
                    GOTO 156
155               CONTINUE
                    CALL GRJMP(SNGL(XT),SNGL(YT))
                    CALL GRDRW(SNGL(XT2),SNGL(YT2))
                    GOTO 156
156             CONTINUE
C
                IF (IFL.GT.0.AND.K.EQ.ISWC(ISW)) THEN
                  IF (.NOT.NLSPLT(NU)) CALL GRNWPN (ICLR(ISW))
                  IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
                  IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
                  ISW=ISW+1
                  CALL GRJMP (SNGL(XTN),SNGL(YTN))
                ENDIF
C
C
150         CONTINUE
            CALL GRNWPN(1)
158       CONTINUE
          CALL GRDSH(1.,0.,1.)
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
     .            call grarrw (sngl(xtn),sngl(ytn),sngl(xtip),
     .                         sngl(ytip),0.4,0.4,0)
                endif
              enddo
            enddo
            CALL GRNWPN(1)
          ENDIF
c  fill in isolated cells (defined in couple_... or any other user
c  segment): all cells j with nstgrd(j)=1
          do j=1,nsurf
            if (nstgrd(j).eq.1) then
              call ncelln(j,ir,ip,it,ia,ib,nr1st,np2nd,nt3rd,nbmlt,
     .                                     nlrad,nlpol,nltor)
              XPS(1)=XPOL(IR,IP)
              IF (NLTRA) XPS(1)=XPS(1)+RMTOR
              YPS(1)=YPOL(IR,IP)
              XPS(2)=XPOL(IR,IP+1)
              IF (NLTRA) XPS(2)=XPS(2)+RMTOR
              YPS(2)=YPOL(IR,IP+1)
              XPS(3)=XPOL(IR+1,IP+1)
              IF (NLTRA) XPS(3)=XPS(3)+RMTOR
              YPS(3)=YPOL(IR+1,IP+1)
              XPS(4)=XPOL(IR+1,IP)
              IF (NLTRA) XPS(4)=XPS(4)+RMTOR
              YPS(4)=YPOL(IR+1,IP)
              XPS(5)=XPOL(IR,IP)
              IF (NLTRA) XPS(5)=XPS(5)+RMTOR
              YPS(5)=YPOL(IR,IP)
              CALL GRFILL(5,XPS,YPS,1,1)
            ENDIF
          enddo
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
            DO 1545 J=1,NSTSI
              IF (NU.EQ.INUMP(J,1)) THEN
                CALL GRNWPN(ILCOL(NLIM+J))
                IF (ILIIN(NLIM+J).LE.0) GOTO 1546
                CALL GRDSH(1.,0.,1.)
                GOTO 1547
              ENDIF
1545        CONTINUE
1546        CALL GRDSH(0.2,0.5,0.2)
1547        CONTINUE
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
              XW1=XPOL(NU,IP)
              IF (NLTRZ) THEN
                IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
                  CALL GRJMP (SNGL(XW1),SNGL(YW1))
                  CALL GRDRW (SNGL(XW1),SNGL(YW2))
                ENDIF
              ELSEIF (NLTRA) THEN
                X=RMTOR+XW1
                Z=TANAL*X
                RR=SQRT(X*X+Z*Z)
                IF (NTTRA.GT.100) THEN
                  WRITE (6,*) 'ERROR IN PLT2D '
                  CALL EXIT
                ENDIF
                DO J=1,NTTRA
                  XX(J)=RR*COS(ZSURF(J))
                  YY(J)=RR*SIN(ZSURF(J))
                ENDDO
                CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
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
              CALL GRDSH (1.,0.,1.)
              CALL GRNWPN (ILCOL(ISTS))
            ELSE
              CALL GRDSH (0.2,0.5,0.2)
              CALL GRNWPN (1)
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
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
1020      CONTINUE
          IF (PLNUMV) THEN
C  Write number of triangle
c  p: halber umfang des dreiecks
c  r: radius des inkreises des dreiecks
            p=p*0.5
            r=((p-dsd(1))*(p-dsd(2))*(p-dsd(3))/p)**0.5
            WRITE (CH,'(I4)') I
            slt=r/2.
            X=XCOMTR(I)
            IF (I.GE.1000) X=X-SLT
            IF (I.GE.100) X=X-SLT
            IF (I.GE.10) X=X-SLT
            Y=YCOMTR(I)
            if (x.lt.xmi2d.or.x.gt.xma2d.or.y.lt.ymi2d.or.y.gt.yma2d)
     .         goto 1010
            slt=slt*sclfcx
            call grchrc(sngl(slt),0.,idummy)
            CALL GRTXT (sngl(X),sngl(Y),4,CH)
            call grchrc(0.3,0.,idummy)
          ENDIF
1010    CONTINUE
        CALL GRDSH (1.,0.,1.)
        CALL GRNWPN (1)
C
      ELSEIF (LEVGEO.EQ.5) THEN
        DO NU=NPLINR,NPLOTR,NPLDLR
          DO J=1,NSTSI
            IF (NU.EQ.INUMP(J,1)) THEN
              CALL GRNWPN(ILCOL(NLIM+J))
              IF (ILIIN(NLIM+J).LE.0) GOTO 1646
              CALL GRDSH(1.,0.,1.)
              GOTO 1647
            ENDIF
          ENDDO
1646      CALL GRDSH(0.2,0.5,0.2)
1647      CONTINUE
          IF (NLSPLT(NU)) CALL GRNWPN(2)
          CALL PLTUSR(.TRUE.,NU)
          CALL GRNWPN(1)
        ENDDO
        CALL GRDSH(1.,0.,1.)
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
          XW1=RSURF(NPLINR)
          XW2=RSURF(NPLOTR)
          IF (NLTRA) XW1=XW1+RMTOR
          IF (NLTRA) XW2=XW2+RMTOR
        ELSEIF (PLCUT(1)) THEN
C Y-Z-PLANE
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
        DO 1171 NU=NPLINP,NPLOTP,NPLDLP
          DO 1172 J=1,NSTSI
            IF (NU.EQ.INUMP(J,2)) THEN
              CALL GRNWPN(ILCOL(NLIM+J))
              IF (ILIIN(NLIM+J).LE.0) GOTO 1173
              CALL GRDSH(1.,0.,1.)
              GOTO 1174
            ENDIF
1172      CONTINUE
1173      CALL GRDSH(0.2,0.5,0.2)
1174      CONTINUE
          IF (NLSPLT(N1ST+NU)) CALL GRNWPN(2)
          YW1=PSURF(NU)
          IF (PLCUT(3)) THEN
            XX(1)=XW1
            XX(2)=XW2
            YY(1)=YW1
            YY(2)=YW1
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
          ELSEIF (PLCUT(1)) THEN
            XX(1)=ZW1
            XX(2)=ZW2
            YY(1)=YW1
            YY(2)=YW1
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
          ENDIF
          CALL GRNWPN(1)
1171    CONTINUE
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
              IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
              IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
              ISW=ISW+1
            ENDIF
C
            IF (IN.EQ.0) IN=4
            GOTO (171,172,173,174,1752),IN
171           CONTINUE
                CALL GRDRW (sngl(XTN),sngl(YTN))
                GOTO 175
172           CONTINUE
                CALL GRDRW(sngl(XT),sngl(YT))
                CALL GRJMP(sngl(XTN),sngl(YTN))
                GOTO 175
173           CONTINUE
                CALL GRJMP(sngl(XT),sngl(YT))
                CALL GRDRW(sngl(XTN),sngl(YTN))
                GOTO 175
174           CONTINUE
                CALL GRJMP(sngl(XTN),sngl(YTN))
                GOTO 175
1752          CONTINUE
                CALL GRJMP(sngl(XT),sngl(YT))
                CALL GRDRW(sngl(XT2),sngl(YT2))
                GOTO 175
175       CONTINUE
C
          IF (IFL.GT.0.AND.I.EQ.ISWC(ISW)) THEN
            IF (.NOT.NLSPLT(N1ST+NU)) CALL GRNWPN (ICLR(ISW))
            IF (IDSH(ISW).EQ.0) CALL GRDSH (1.,0.,1.)
            IF (IDSH(ISW).EQ.1) CALL GRDSH (0.2,0.5,0.2)
            ISW=ISW+1
            CALL GRJMP (sngl(XTN),sngl(YTN))
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
     .          call grarrw (sngl(xtn),sngl(ytn),sngl(xtip),sngl(ytip),
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
c slmod begin - f90
c...bug?
c XAA and XIA are not defined in common blocks, but YAA and YIA are defined
c in CGRID.  This is the only reference to XAA and YAA in the code.
          WRITE(0,*) 'WARNING: Assigning zero XW1 and XW2 in PLOT.F'

          XW1=0
          XW2=0
c
c          XW1=XIA
c          XW2=XAA
c slmod end
          IF (NLRAD) THEN
            XW1=RSURF(NPLINR)
            XW2=RSURF(NPLOTR)
          ENDIF
        ENDIF
        DO 204 NU=NPLINT,NPLOTT,NPLDLT
          DO 205 J=1,NSTSI
            IF (NU.EQ.INUMP(J,3)) THEN
              CALL GRNWPN(ILCOL(NLIM+J))
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
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
          ELSEIF (PLCUT(2)) THEN
            XX(1)=XW1
            XX(2)=XW2
            YY(1)=ZW1
            YY(2)=ZW1
            CALL PLTLNE(2,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
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
            IF (ILTOR(I).GT.0) THEN
              CALL SETROT (AFF,AFFI,1,0.0D0,1.0D0,0.0D0,
     .                    -ZZONE(ILTOR(I))*180.D0/PIA)
              CALL ROTADD (AFF,AFFI,I,I)
            ENDIF
          ENDDO
        ENDIF
C
        CALL PLTADD(1,NLIMI)
C
        IF (NLTRA) THEN
          DO I=1,NLIMI
            IF (ILTOR(I).GT.0) THEN
              CALL SETROT (AFF,AFFI,1,0.0D0,1.0D0,0.0D0,
     .                     ZZONE(ILTOR(I))*180.D0/PIA)
              CALL ROTADD (AFF,AFFI,I,I)
            ENDIF
          ENDDO
          CALL XSHADD(-RMTOR,1,NLIMI)
        ENDIF
C
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
C
C  IFLAG.EQ.0 ONLY SYMBOL AT XPLO,YPLO,ZPLO
C  IFLAG.NE.0 TRACK FROM LAST POSITION (PREVIOUS CALL) TO XPLO,YPLO,ZPLO
C  ISYM       NUMBER OF SYMBOL FOR THE CURRENT EVENT
C
C  TEXT FOR PARTICLE-HISTORIES PLOT, ONLY AT FIRST CALL TO THIS ENTRY
C
c slmod begin
      IF (iflag.EQ.0) count = count + 1.0

      WRITE(80,'(I8,3F14.6,3I8,1P,1E10.2,0P,F8.4,4I5)') 
     .  NINT(count),XPLO,YPLO,ZPLO,IFLAG,ISYM,npanu,e0,weight,nacell,
     .  ifpath,iupdte,ityp
c slmod end
      XN=XN2D
      YN=YN2D
      FX=FX2D
      FY=FY2D
      IF (IWRIT.EQ.0.AND.PLHST) THEN
        IWRIT=1
        IF (.NOT.NLPL3D) CALL GRSCLV(0.,0.,sngl(ABSMAX),sngl(ORDMAX))
        XNP05=XN+0.5/FX
        DO IA=1,NTXHST
          YYIA=YN-(0.75*(IA-1))/FY
          CALL GRJMPS (sngl(XN),sngl(YYIA+0.15),ISPL(IA))
          CALL GRTXT (sngl(XNP05),sngl(YYIA),20,TXTHST(IA))
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
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (sngl(XN1),sngl(YNP))
          CALL GRDRW (sngl(XN2),sngl(YNP))
          CALL GRTXT (sngl(XN2),sngl(YNP-0.15),8,TEXTS(I))
510     CONTINUE
        DO 511 I=1,NMOLI
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (sngl(XN1),sngl(YNP))
          CALL GRDRW (sngl(XN2),sngl(YNP))
          CALL GRTXT (sngl(XN2),sngl(YNP-0.15),8,TEXTS(NSPA+I))
511     CONTINUE
        DO 512 I=1,NIONI
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (sngl(XN1),sngl(YNP))
          CALL GRDRW (sngl(XN2),sngl(YNP))
          CALL GRTXT (sngl(XN2),sngl(YNP-0.15),8,TEXTS(NSPAM+I))
512     CONTINUE
        DO 513 I=1,NPLSI
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          CALL GRNWPN (ICP)
          CALL GRJMP (sngl(XN1),sngl(YNP))
          CALL GRDRW (sngl(XN2),sngl(YNP))
          CALL GRTXT (sngl(XN2),sngl(YNP-0.15),8,TEXTS(NSPAMI+I))
513     CONTINUE
        IF (.NOT.NLPL3D)
     .     CALL GRSCLV(sngl(XMI2D),sngl(YMI2D),sngl(XMA2D),sngl(YMA2D))
      ENDIF
C
C  WRITE TRACK DATA
C
      IF (TRCHST) THEN
        CALL LEER(1)
        WRITE (6,*) TXTHST(ISYM)
        IF (ISYM.EQ.1)  CALL MASJ1('NPANU   ',NPANU)
        WRITE (6,'(1X,A8)') TEXTS(ISPZ)
        CALL MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
        CALL MASR6 ('VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .               VELX,VELY,VELZ,VEL,E0,WEIGHT)
        CALL MASR1 ('TIME    ',TIME)
        CALL MASJ1 ('IGENER  ',IGENER)
        CALL MASJ4 ('NRCELL,IPOLG,NACELL,NBLOCK      ',
     .               NRCELL,IPOLG,NACELL,NBLOCK)
        IF (NLTOR) THEN
          CALL MASJ1 ('NTCELL  ',NTCELL)
        ENDIF
        IF (NLTRA) THEN
          CALL MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI)
        ENDIF
        IF (NLPOL) THEN
          CALL MASJ1 ('NPCELL  ',NPCELL)
        ENDIF
        IF ((ISYM.GE.6.AND.ISYM.LE.10).OR.
     .      (ISYM.EQ.1.AND.NLSRF(ISTRA))) THEN
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
      IF (IN.LE.2) CALL GRJMP (sngl(XWO),sngl(YWO))
      IF (ILINIE.EQ.0) GOTO 404
      GOTO (401,402,403,404,405),IN
401     CALL GRDRW (sngl(XWN),sngl(YWN))
        GOTO 406
402     CALL GRDRW (sngl(XT),sngl(YT))
        CALL GRJMP (sngl(XWN),sngl(YWN))
        GOTO 406
403     CALL GRJMP (sngl(XT),sngl(YT))
        CALL GRDRW (sngl(XWN),sngl(YWN))
        GOTO 406
404     CALL GRJMP (sngl(XWN),sngl(YWN))
        GOTO 406
405     CALL GRJMP (sngl(XT),sngl(YT))
        CALL GRDRW (sngl(XT2),sngl(YT2))
        GOTO 406
406   CONTINUE
C
C  PLOT SYMBOL
C
      IF (TESTN.EQ.0.) THEN
        DO 408 J=1,8
          IF (ISYM.EQ.ISYPLT(J).OR.ISYM.EQ.15) THEN
            CALL GRJMPS(sngl(XWN),sngl(YWN),ISPL(ISYM))
            GOTO 409
          ENDIF
408     CONTINUE
409     CONTINUE
      ENDIF
      IF (NLTRC.AND.TRCHST)
     .WRITE (6,*) 'ICOLOR,IN,ISYM,IFLAG ',ICOLOR,IN,ISYM,IFLAG
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
      CALL EXIT
991   CONTINUE
      WRITE (6,*) 'OPTION IN PLT2D NOT READY: Y- OR POLOIDAL GRID'
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'OPTION IN PLT2D NOT READY: Z- OR TOROIDAL GRID'
      CALL EXIT
      END
C
      SUBROUTINE PLGELL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  AT ENTRY PLGELR:
C  THIS SUBROUTINE FITS A POLYGON OF LENGTH N
C  TO A RADIAL SURFACE (E.G.: ELLIPSE, CENTERED AT X=EP,Y=0.,...)
C  THE TWO SEMIAXES MEASURED FROM X=EP ARE:
C  AHALB: X DIRECTION, BHALB: Y DIRECTION
C  DM IS THE MINIMAL DISTANCE OF TWO POINTS FOR PLOTTING. E.G:
C  DM=0.01*MIN(CHXM,CHYM), IF PLGELL IS CALLED FROM SUBR. PLT2D
C  DM=0.0, IF PLGELL IS CALLED FORM SUBR. PLT3D
C  AT ENTRY PLGELP:
C  THIS SUBROUTINE FITS A POLYGON OF LENGTH N
C  TO A POLOIDAL SURFACE (E.G.: THETA=CONST, IF NLCRC,...)
C
C  OUTPUT: XX(J),YY(J),J=1,NRET (ACTUAL NUMBER OF POINTS ON POLYGON,
C                                PROBABLY LESS THEN N)
C
      INCLUDE 'CCONA'
      DIMENSION XX(*),YY(*),THETA(*)
C
      ENTRY PLGELR (AHALB,EP,ELL,TRI,DM,N,XX,YY,NRET,THETA,NP)
C
      AHB=AHALB
      BHB=AHALB*ELL
      ATR=TRI
      BTR=-TRI*ELL
      IF (NP.GT.1) THEN
        DX=(THETA(NP)-THETA(1))/DBLE(N-1)
      ELSE
        DX=PI2A/DBLE(N-1)
      ENDIF
C  START AT THE RIGHTMOST POINT THETA(1)=0
      XOLD=EP+AHB*COS(THETA(1))+ATR*COS(THETA(1)*2.)
      YOLD=   BHB*SIN(THETA(1))+BTR*SIN(THETA(1)*2.)
      XX(1)=XOLD
      YY(1)=YOLD
      ICOU=1
C
      DO 100 I=2,N
        XAD=DBLE(I-1)*DX
        XH=THETA(1)+XAD
        XNEW=EP+AHB*COS(XH)+ATR*COS(2.*XH)
        YNEW=   BHB*SIN(XH)+BTR*SIN(2.*XH)
C
C   SKIP THIS PART?
        DMM=SQRT((XNEW-XOLD)**2+(YNEW-YOLD)**2)
        IF (DMM.LT.DM.AND.I.LT.N) GOTO 100
C   NO!
        ICOU=ICOU+1
        XX(ICOU)=XNEW
        YY(ICOU)=YNEW
        XOLD=XX(ICOU)
        YOLD=YY(ICOU)
100   CONTINUE
C
      NRET=ICOU
      RETURN
C
      ENTRY PLGELP
      RETURN
      END
C
      SUBROUTINE TSTCHM(I1,XTN,YTN,XT,YT,IINDEX,TESTN,XMI,XMA,YMI,YMA,
     .                     XT2,YT2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C INPUT:
C I1=0:
C       FIRST CALL FOR THIS TRACK. INITIALIZE CHECK DATA AND JUMP ONLY
C       SET XTO,YTO=XTN,YTN
C    ELSE:
C       CARRY OUT CHECK FOR TRACK XTO,YTO  TO  XTN,YTN AND RETURN DATA
C       SET XTO,YTO=XTN,YTN
C OUTPUT:
C     TESTN: FLAG FOR NEW POINT XTN,YTN
C     TESTN= 0:  COMPLETELY WITHIN CHAMBER
C     TESTN= 1:  EITHER WITHIN X RANGE BUT OUTSIDE Y RANGE,
C                OR     WITHIN Y RANGE BUT OUTSIDE X RANGE
C     TESTN= 2:  BOTH: OUTSIDE X RANGE AND OUTSIDE Y RANGE
C
C     IINDEX: FLAG FOR TRACK XTO,YTO TO  XTN,YTN
C     IINDEX= 0:  INITIAL CALL (I1=0), NO PLOT
C     IINDEX= 1:  COMPLETE TRACK INSIDE
C     IINDEX= 2:  OLD POINT INSIDE, NEW POINT OUTSIDE
C     IINDEX= 3:  NEW POINT INSIDE, OLD POINT OUTSIDE
C     IINDEX= 4:  BOTH POINTS OUTSIDE, AND NO INTERSECTIONS. NO PLOT
C     IINDEX= 5:  BOTH POINTS OUTSIDE, BUT PART OF TRACK INSIDE.  PLOT
C     XT,YT: INTERSETION POINT ON BOUNDARY  (IF ANY)
C     XT2,YT2: SECOND INTERSECTION POINT ON BOUNDARY (IF ANY)
C
      DIMENSION XS(2),YS(2)
C
      SAVE TESTO,TBDXO,TBDYO,XTO,YTO
C
      TSEC(ZXO,ZYO,ZXN,ZYN,ZBX1,ZBY1,ZBX2,ZBY2)=
     .     ((ZXO-ZBX1)*(ZYN-ZYO)-(ZYO-ZBY1)*(ZXN-ZXO))/
     .     ((ZBX2-ZBX1)*(ZYN-ZYO)-(ZBY2-ZBY1)*(ZXN-ZXO)+1.D-30)
C
      IF (I1.EQ.0) THEN
        TBDXN=0.
        TBDYN=0.
        IF (XTN.LT.XMI) TBDXN=-1.
        IF (XTN.GT.XMA) TBDXN=1.
        IF (YTN.LT.YMI) TBDYN=-1.
        IF (YTN.GT.YMA) TBDYN=1.
        TESTN=ABS(TBDXN)+ABS(TBDYN)
        IINDEX=0
        GOTO 100
      ENDIF
C
      TBDXN=0.
      TBDYN=0.
      IF (XTN.LT.XMI) TBDXN=-1.
      IF (XTN.GT.XMA) TBDXN=1.
      IF (YTN.LT.YMI) TBDYN=-1.
      IF (YTN.GT.YMA) TBDYN=1.
      TESTN=ABS(TBDXN)+ABS(TBDYN)
      TEST=TESTN+TESTO
      IF (TEST.EQ.0.) THEN
        IINDEX=1
      ELSEIF (TESTO.EQ.0..AND.TESTN.NE.0.) THEN
        IF (TBDXN.LT.0.) THEN
          XT=XMI
          TIME=(XT-XTO)/(XTN-XTO)
          YT=YTO+TIME*(YTN-YTO)
          IF (YT.LE.YMA.AND.YT.GE.YMI) GOTO 146
        ELSEIF (TBDXN.GT.0.) THEN
          XT=XMA
          TIME=(XT-XTO)/(XTN-XTO)
          YT=YTO+TIME*(YTN-YTO)
          IF (YT.LE.YMA.AND.YT.GE.YMI) GOTO 146
        ENDIF
        IF (TBDYN.LT.0.) THEN
          YT=YMI
          TIME=(YT-YTO)/(YTN-YTO)
          XT=XTO+TIME*(XTN-XTO)
        ELSEIF (TBDYN.GT.0.) THEN
          YT=YMA
          TIME=(YT-YTO)/(YTN-YTO)
          XT=XTO+TIME*(XTN-XTO)
        ENDIF
146     CONTINUE
        IINDEX=2
      ELSEIF (TESTN.EQ.0..AND.TESTO.NE.0.) THEN
        IF (TBDXO.LT.0.) THEN
          XT=XMI
          TIME=(XT-XTN)/(XTN-XTO)
          YT=YTN+TIME*(YTN-YTO)
          IF (YT.LE.YMA.AND.YT.GE.YMI) GOTO 148
        ELSEIF (TBDXO.GT.0.) THEN
          XT=XMA
          TIME=(XT-XTN)/(XTN-XTO)
          YT=YTN+TIME*(YTN-YTO)
          IF (YT.LE.YMA.AND.YT.GE.YMI) GOTO 148
        ENDIF
        IF (TBDYO.LT.0.) THEN
          YT=YMI
          TIME=(YT-YTN)/(YTN-YTO)
          XT=XTN+TIME*(XTN-XTO)
          ELSEIF (TBDYO.GT.0.) THEN
          YT=YMA
          TIME=(YT-YTN)/(YTN-YTO)
          XT=XTN+TIME*(XTN-XTO)
        ENDIF
148     CONTINUE
        IINDEX=3
      ELSE
        TESTX=ABS(TBDXO+TBDXN)
        TESTY=ABS(TBDYO+TBDYN)
        IF (TESTX.EQ.2..OR.TESTY.EQ.2.) THEN
          IINDEX=4
        ELSE
C   TEST ALL BOUNDARIES
          IT=0
C   LOWER BOUNDARY
          T=TSEC(XTO,YTO,XTN,YTN,XMI,YMI,XMA,YMI)
          IF (T.GE.0..AND.T.LE.1.) THEN
            IT=IT+1
            XS(IT)=XMI+T*(XMA-XMI)
            YS(IT)=YMI
          ENDIF
C   UPPER BOUNDARY
          T=TSEC(XTO,YTO,XTN,YTN,XMI,YMA,XMA,YMA)
          IF (T.GE.0..AND.T.LE.1.) THEN
            IT=IT+1
            XS(IT)=XMI+T*(XMA-XMI)
            YS(IT)=YMA
          ENDIF
C   LEFT BOUNDARY
          T=TSEC(XTO,YTO,XTN,YTN,XMI,YMI,XMI,YMA)
          IF (T.GE.0..AND.T.LE.1.) THEN
            IT=IT+1
            XS(IT)=XMI
            YS(IT)=YMI+T*(YMA-YMI)
          ENDIF
C   RIGHT BOUNDARY
          T=TSEC(XTO,YTO,XTN,YTN,XMA,YMI,XMA,YMA)
          IF (T.GE.0..AND.T.LE.1.) THEN
            IT=IT+1
            XS(IT)=XMA
            YS(IT)=YMI+T*(YMA-YMI)
          ENDIF
C
          IF (IT.EQ.0) THEN
            IINDEX=4
          ELSEIF (IT.EQ.2) THEN
            IINDEX=5
            XT=XS(1)
            YT=YS(1)
            XT2=XS(2)
            YT2=YS(2)
          ELSE
            WRITE (6,*) ' NONSENSE IN TSTCHM '
            WRITE (6,*) ' XTO,YTO ',XTO,YTO
            WRITE (6,*) ' XTN,YTN ',XTN,YTN
            IINDEX=4
          ENDIF
        ENDIF
      ENDIF
100   TESTO=TESTN
      TBDXO=TBDXN
      TBDYO=TBDYN
      XTO=XTN
      YTO=YTN
      RETURN
      END
C
      SUBROUTINE PLTLNE(NRET,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(*),YY(*)
c
      DO 148 I=1,NRET
        IFL=I-1
        XTN=XX(I)
        YTN=YY(I)
        CALL TSTCHM(IFL,XTN,YTN,XT,YT,IN,TESTN,
     .              XMI2D,XMA2D,YMI2D,YMA2D,XT2,YT2)
C
        IF (IN.EQ.0) IN=4
        GOTO (1481,1482,1483,1484,14852),IN
1481      CONTINUE
            CALL GRDRW (sngl(XTN),sngl(YTN))
            GOTO 1485
1482      CONTINUE
            CALL GRDRW(sngl(XT),sngl(YT))
            CALL GRJMP(sngl(XTN),sngl(YTN))
            GOTO 1485
1483      CONTINUE
            CALL GRJMP(sngl(XT),sngl(YT))
            CALL GRDRW(sngl(XTN),sngl(YTN))
            GOTO 1485
1484      CONTINUE
            CALL GRJMP(sngl(XTN),sngl(YTN))
            GOTO 1485
14852     CONTINUE
            CALL GRJMP(sngl(XT),sngl(YT))
            CALL GRDRW(sngl(XT2),sngl(YT2))
            GOTO 1485
1485    CONTINUE
C
148   CONTINUE
      return
      end
C
      SUBROUTINE PLTADD (MANF,MEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  PLOT ADDITIONAL SURFACES ISURF=MANF,MEND, PROJECTED INTO A PLANE
C  IF LZR,     PLOT IS DONE IN THIS ROUTINE
C  IF.NOT.LZR, PLOT DATA ARE COLLECTED IN COMMON CRECH
C
      DIMENSION  G(4,2),R(4,2),AFF(3,3),AFFI(3,3)
      INCLUDE 'PARMMOD'
      INCLUDE 'CLMSUR'
      INCLUDE 'CTRCEI'
      INCLUDE 'CPLOT'
      INCLUDE 'CADGEO'
      INCLUDE 'CLGIN'
      INCLUDE 'CCONA'
      INCLUDE 'CTEXT'
      INCLUDE 'CRECH'
      CHARACTER*1 CH1
      CHARACTER*2 CH2
      CHARACTER*3 CH3
      EXTERNAL ELLO,ELLU,PARA1,PARA2O,PARA2U,HYPP,HYPM,
     .         GERADY,GERADX
      LOGICAL L1,L2,L3,L4,L5
      LOGICAL LBOX
      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
C
      XCHL=CH2X0-CH2MX
      XCHR=CH2X0+CH2MX
      YCHL=CH2Y0-CH2MY
      YCHR=CH2Y0+CH2MY
      DCHX=2.*CH2MX
      DCHY=2.*CH2MY
      G(1,1)=XCHL
      G(1,2)=YCHL
      G(2,1)=XCHL
      G(2,2)=YCHR
      G(3,1)=XCHR
      G(3,2)=YCHL
      G(4,1)=XCHL
      G(4,2)=YCHL
      R(1,1)=0.
      R(1,2)=DCHY
      R(2,1)=DCHX
      R(2,2)=0.
      R(3,1)=0.
      R(3,2)=DCHY
      R(4,1)=DCHX
      R(4,2)=0.
C
      DO 2000 ICUT=1,3
      IF (.NOT.PLCUT(ICUT)) GOTO 2000
      IF (ICUT.EQ.1) THEN
        CALL ZEROA2(AFF,3,3)
        AFF(1,3)=1.
        AFF(2,2)=1.
        AFF(3,1)=-1.
        CALL ZEROA2(AFFI,3,3)
        AFFI(1,3)=-1.
        AFFI(2,2)=1.
        AFFI(3,1)=1.
        CALL ROTADD(AFF,AFFI,MANF,MEND)
      ELSEIF (ICUT.EQ.2) THEN
        CALL ZEROA2(AFF,3,3)
        AFF(1,1)=1.
        AFF(2,3)=1.
        AFF(3,2)=-1.
        CALL ZEROA2(AFFI,3,3)
        AFFI(1,1)=1.
        AFFI(2,3)=-1.
        AFFI(3,2)=1.
        CALL ROTADD(AFF,AFFI,MANF,MEND)
      ELSEIF (ICUT.EQ.3) THEN
C  NO ROTATION NEEDED
      ENDIF
C
      IF (CH2Z0.NE.0.) CALL ZSHADD(-CH2Z0,MANF,MEND)
C
C   PLOT SURFACES
      ZPLT=0.
      INSTOR=0
C
      DO 200 ISURF=MANF,MEND
        J=ISURF
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) INOSF=ISURF
C
C  PLOT SURFACE NO. ISURF
C
        IF (IGJUM0(ISURF).NE.0) GOTO 200
        IF (TRCPLT) THEN
          CALL LEER(1)
          WRITE (6,*) TXTSFL(ISURF)
        ENDIF
        LBOX=ILBOX(ISURF).GT.0.AND.LZR
        IF (LZR) THEN
          ICOLOR=ILCOL(ISURF)
          CALL GRNWPN(ICOLOR)
          IF (ILIIN(ISURF).LE.0) THEN
            CALL GRDSH(0.2,0.5,0.2)
          ELSE
            CALL GRDSH(1.,0.,1.)
          ENDIF
        ENDIF
C
C   PLOT SCHNITT MIT ALLGEMEINEM 3, 4 ODER 5 ECK
C
        IF (RLB(ISURF).GE.2.) THEN
          AMXPX=MAX(P1(1,ISURF),P2(1,ISURF))
          AMNPX=MIN(P1(1,ISURF),P2(1,ISURF))
          AMXPY=MAX(P1(2,ISURF),P2(2,ISURF))
          AMNPY=MIN(P1(2,ISURF),P2(2,ISURF))
          AMXPZ=MAX(P1(3,ISURF),P2(3,ISURF))
          AMNPZ=MIN(P1(3,ISURF),P2(3,ISURF))
          IF (RLB(ISURF).GE.3.) THEN
            AMXPX=MAX(AMXPX,P3(1,ISURF))
            AMNPX=MIN(AMNPX,P3(1,ISURF))
            AMXPY=MAX(AMXPY,P3(2,ISURF))
            AMNPY=MIN(AMNPY,P3(2,ISURF))
            AMXPZ=MAX(AMXPZ,P3(3,ISURF))
            AMNPZ=MIN(AMNPZ,P3(3,ISURF))
          ENDIF
          IF (RLB(ISURF).GE.4.) THEN
            AMXPX=MAX(AMXPX,P4(1,ISURF))
            AMNPX=MIN(AMNPX,P4(1,ISURF))
            AMXPY=MAX(AMXPY,P4(2,ISURF))
            AMNPY=MIN(AMNPY,P4(2,ISURF))
            AMXPZ=MAX(AMXPZ,P4(3,ISURF))
            AMNPZ=MIN(AMNPZ,P4(3,ISURF))
          ENDIF
          IF (RLB(ISURF).GE.5.) THEN
            AMXPX=MAX(AMXPX,P5(1,ISURF))
            AMNPX=MIN(AMNPX,P5(1,ISURF))
            AMXPY=MAX(AMXPY,P5(2,ISURF))
            AMNPY=MIN(AMNPY,P5(2,ISURF))
            AMXPZ=MAX(AMXPZ,P5(3,ISURF))
            AMNPZ=MIN(AMNPZ,P5(3,ISURF))
          ENDIF
          IF(AMXPX.LT.XCHL.OR.AMNPX.GT.XCHR.OR.
     .       AMXPY.LT.YCHL.OR.AMNPY.GT.YCHR.OR.
     .       AMXPZ.LT.ZPLT.OR.AMNPZ.GT.ZPLT) THEN
            IF (TRCPLT) WRITE (6,6660)
            GOTO 200
          ENDIF
          IF (TRCPLT) THEN
            WRITE (6,*) ' AMNPX,AMNPY,AMNPZ ',AMNPX,AMNPY,AMNPZ
            WRITE (6,*) ' AMXPX,AMXPY,AMXPZ ',AMXPX,AMXPY,AMXPZ
            WRITE (6,*) ' A0,A1,A2,A3 ',
     .                   A0LM(ISURF),A1LM(ISURF),A2LM(ISURF),A3LM(ISURF)
          ENDIF
c
          IF (ABS(A3LM(ISURF)).GE.1.-EPS10) THEN
C   N ECK PARALLEL ZU Z=ZPLT
            IF (TRCPLT) WRITE (6,*) 'N ECK PARALLEL ZU Z=ZPLT'
            IF (ABS(A0LM(ISURF)-ZPLT).GT.EPS10) GOTO 200
C
            L1=P1(1,ISURF).GE.XCHL.AND.P1(1,ISURF).LE.XCHR.AND.
     .         P1(2,ISURF).GE.YCHL.AND.P1(2,ISURF).LE.YCHR
            L2=P2(1,ISURF).GE.XCHL.AND.P2(1,ISURF).LE.XCHR.AND.
     .         P2(2,ISURF).GE.YCHL.AND.P2(2,ISURF).LE.YCHR
            L3=.FALSE.
            IF (RLB(ISURF).GE.3.) THEN
              L3=P3(1,ISURF).GE.XCHL.AND.P3(1,ISURF).LE.XCHR.AND.
     .           P3(2,ISURF).GE.YCHL.AND.P3(2,ISURF).LE.YCHR
            ENDIF
            L4=.FALSE.
            IF (RLB(ISURF).GE.4.) THEN
              L4=P4(1,ISURF).GE.XCHL.AND.P4(1,ISURF).LE.XCHR.AND.
     .           P4(2,ISURF).GE.YCHL.AND.P4(2,ISURF).LE.YCHR
            ENDIF
            L5=.FALSE.
            IF (RLB(ISURF).GE.5.) THEN
              L5=P5(1,ISURF).GE.XCHL.AND.P5(1,ISURF).LE.XCHR.AND.
     .           P5(2,ISURF).GE.YCHL.AND.P5(2,ISURF).LE.YCHR
            ENDIF
            IF (TRCPLT) WRITE (6,*) 'L1,L2,L3,L4,L5 ',L1,L2,L3,L4,L5
C   PLOT P1, P2
            IF (L1.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P1(1,ISURF)),sngl(P1(2,ISURF)))
                CALL GRDRW (sngl(P2(1,ISURF)),sngl(P2(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P1(1,ISURF),P1(2,ISURF),0)
                CALL STCOOR (P2(1,ISURF),P2(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P1(1,ISURF),P1(2,ISURF),P2(1,ISURF),
     .                    P2(2,ISURF),L1,L2,EPS10)
            ENDIF
C   PLOT P1, P3
            IF (L1.AND.L3) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P1(1,ISURF)),sngl(P1(2,ISURF)))
                CALL GRDRW (sngl(P3(1,ISURF)),sngl(P3(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P1(1,ISURF),P1(2,ISURF),0)
                CALL STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P1(1,ISURF),P1(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L1,L3,EPS10)
            ENDIF
C   PLOT P2, P3
            IF (L2.AND.L3.AND..NOT.L4) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P2(1,ISURF)),sngl(P2(2,ISURF)))
                CALL GRDRW (sngl(P3(1,ISURF)),sngl(P3(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P2(1,ISURF),P2(2,ISURF),0)
                CALL STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P2(1,ISURF),P2(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L2,L3,EPS10)
            ENDIF
C   PLOT P4, P2
            IF (L4.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P4(1,ISURF)),sngl(P4(2,ISURF)))
                CALL GRDRW (sngl(P2(1,ISURF)),sngl(P2(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P4(1,ISURF),P4(2,ISURF),0)
                CALL STCOOR (P2(1,ISURF),P2(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P4(1,ISURF),P4(2,ISURF),P2(1,ISURF),
     .                    P2(2,ISURF),L4,L2,EPS10)
            ENDIF
C   PLOT P4, P3
            IF (L4.AND.L3.AND..NOT.L5) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P4(1,ISURF)),sngl(P4(2,ISURF)))
                CALL GRDRW (sngl(P3(1,ISURF)),sngl(P3(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P4(1,ISURF),P4(2,ISURF),0)
                CALL STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P4(1,ISURF),P4(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L4,L3,EPS10)
            ENDIF
C   PLOT P4, P5
            IF (L4.AND.L5) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P4(1,ISURF)),sngl(P4(2,ISURF)))
                CALL GRDRW (sngl(P5(1,ISURF)),sngl(P5(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P4(1,ISURF),P4(2,ISURF),0)
                CALL STCOOR (P5(1,ISURF),P5(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P4(1,ISURF),P4(2,ISURF),P5(1,ISURF),
     .                    P5(2,ISURF),L4,L5,EPS10)
            ENDIF
C   PLOT P5, P3
            IF (L5.AND.L3) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(P5(1,ISURF)),sngl(P5(2,ISURF)))
                CALL GRDRW (sngl(P3(1,ISURF)),sngl(P3(2,ISURF)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (P5(1,ISURF),P5(2,ISURF),0)
                CALL STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,P5(1,ISURF),P5(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L5,L3,EPS10)
            ENDIF
C   BERECHNE SCHNITTGERADE
          ELSE
            VR1=-A2LM(ISURF)
            VR2=A1LM(ISURF)
            VR3=0.
            V3=ZPLT
            IF (TRCPLT) WRITE (6,*) ' VR1,VR2,VR3 ',VR1,VR2,VR3
            IF (ABS(A1LM(ISURF)).GT.ABS(A2LM(ISURF))) THEN
              V1=(-A0LM(ISURF)-A3LM(ISURF)*V3)/A1LM(ISURF)
              V2=0.
            ELSE
              V1=0.
              V2=(-A0LM(ISURF)-A3LM(ISURF)*V3)/A2LM(ISURF)
            ENDIF
            IS=0
C
            XMU=FMU(V1,V2,V3,VR1,VR2,VR3,
     .              P1(1,ISURF),P1(2,ISURF),P1(3,ISURF),
     .              P2(1,ISURF)-P1(1,ISURF),P2(2,ISURF)-P1(2,ISURF),
     .              P2(3,ISURF)-P1(3,ISURF),EPS10)
            IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
              IS=1
              XP1=P1(1,ISURF)+XMU*(P2(1,ISURF)-P1(1,ISURF))
              YP1=P1(2,ISURF)+XMU*(P2(2,ISURF)-P1(2,ISURF))
            ENDIF
            IF (RLB(ISURF).GE.3.) THEN
              XMU=FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P1(1,ISURF),P1(2,ISURF),P1(3,ISURF),
     .                P3(1,ISURF)-P1(1,ISURF),P3(2,ISURF)-P1(2,ISURF),
     .                P3(3,ISURF)-P1(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P1(1,ISURF)+XMU*(P3(1,ISURF)-P1(1,ISURF))
                  YP1=P1(2,ISURF)+XMU*(P3(2,ISURF)-P1(2,ISURF))
                ELSE
                  XP2=P1(1,ISURF)+XMU*(P3(1,ISURF)-P1(1,ISURF))
                  YP2=P1(2,ISURF)+XMU*(P3(2,ISURF)-P1(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
            ENDIF
            IF (RLB(ISURF).GE.4.) THEN
              XMU=FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P4(1,ISURF),P4(2,ISURF),P4(3,ISURF),
     .                P2(1,ISURF)-P4(1,ISURF),P2(2,ISURF)-P4(2,ISURF),
     .                P2(3,ISURF)-P4(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P4(1,ISURF)+XMU*(P2(1,ISURF)-P4(1,ISURF))
                  YP1=P4(2,ISURF)+XMU*(P2(2,ISURF)-P4(2,ISURF))
                ELSE
                  XP2=P4(1,ISURF)+XMU*(P2(1,ISURF)-P4(1,ISURF))
                  YP2=P4(2,ISURF)+XMU*(P2(2,ISURF)-P4(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
              XMU=FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P4(1,ISURF),P4(2,ISURF),P4(3,ISURF),
     .                P3(1,ISURF)-P4(1,ISURF),P3(2,ISURF)-P4(2,ISURF),
     .                P3(3,ISURF)-P4(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P4(1,ISURF)+XMU*(P3(1,ISURF)-P4(1,ISURF))
                  YP1=P4(2,ISURF)+XMU*(P3(2,ISURF)-P4(2,ISURF))
                ELSE
                  XP2=P4(1,ISURF)+XMU*(P3(1,ISURF)-P4(1,ISURF))
                  YP2=P4(2,ISURF)+XMU*(P3(2,ISURF)-P4(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
            ENDIF
            IF (RLB(ISURF).GE.5) THEN
              XMU=FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P4(1,ISURF),P4(2,ISURF),P4(3,ISURF),
     .                P5(1,ISURF)-P4(1,ISURF),P5(2,ISURF)-P4(2,ISURF),
     .                P5(3,ISURF)-P4(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P4(1,ISURF)+XMU*(P5(1,ISURF)-P4(1,ISURF))
                  YP1=P4(2,ISURF)+XMU*(P5(2,ISURF)-P4(2,ISURF))
                ELSE
                  XP2=P4(1,ISURF)+XMU*(P5(1,ISURF)-P4(1,ISURF))
                  YP2=P4(2,ISURF)+XMU*(P5(2,ISURF)-P4(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
              XMU=FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P5(1,ISURF),P5(2,ISURF),P5(3,ISURF),
     .                P3(1,ISURF)-P5(1,ISURF),P3(2,ISURF)-P5(2,ISURF),
     .                P3(3,ISURF)-P5(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P5(1,ISURF)+XMU*(P3(1,ISURF)-P5(1,ISURF))
                  YP1=P5(2,ISURF)+XMU*(P3(2,ISURF)-P5(2,ISURF))
                ELSE
                  XP2=P5(1,ISURF)+XMU*(P3(1,ISURF)-P5(1,ISURF))
                  YP2=P5(2,ISURF)+XMU*(P3(2,ISURF)-P5(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
            ENDIF
C
            IF (IS.EQ.0) THEN
              WRITE (6,*) ' SCHNITTGERADE AUSSERHALB DES VIERECKS'
            ELSEIF (IS.EQ.1.AND.
     .              XP1.GE.XCHL.AND.XP1.LE.XCHR.AND.YP1.GE.YCHL.AND.
     .              YP1.LE.YCHR) THEN
              IF (LZR) THEN
                CALL GRSPTS (30)
                CALL GRJMP (sngl(XP1),sngl(YP1))
                CALL GRDRW (sngl(XP1),sngl(YP1))
                CALL GRSPTS(16)
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP1,YP1,0)
            ENDIF
            GOTO 200
C
90          CONTINUE
            L1=XP1.GE.XCHL.AND.XP1.LE.XCHR.AND.
     .         YP1.GE.YCHL.AND.YP1.LE.YCHR
            L2=XP2.GE.XCHL.AND.XP2.LE.XCHR.AND.
     .         YP2.GE.YCHL.AND.YP2.LE.YCHR
            IF (L1.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP (sngl(XP1),sngl(YP1))
                CALL GRDRW (sngl(XP2),sngl(YP2))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL STCOOR (XP1,YP1,0)
                CALL STCOOR (XP2,YP2,1)
              ENDIF
            ELSE
              CALL PSIDE (G,R,XP1,YP1,XP2,YP2,L1,L2,EPS10)
            ENDIF
          ENDIF
          GOTO 200
        ENDIF
C
C   ALLE ANDEREN FALLE: RLB.LT.2.
C
        IF (RLB(ISURF).EQ.1.) THEN
          IF (ZLIMS1(ISURF,1).GT.ZPLT.OR.ZLIMS2(ISURF,1).LT.ZPLT) THEN
            IF (TRCPLT) THEN
              WRITE (6,*) 'ZLIMS1,ZPLT,ZLIMS2 ',
     .                     ZLIMS1(ISURF,1),ZPLT,ZLIMS2(ISURF,1)
            ENDIF
            GOTO 200
          ENDIF
        ENDIF
        XM=0.
        YM=0.
        AK=A4LM(ISURF)
        BK=A7LM(ISURF)
        CK=A5LM(ISURF)
        DK=A1LM(ISURF)+A8LM(ISURF)*ZPLT
        EK=A2LM(ISURF)+A9LM(ISURF)*ZPLT
        FK=A0LM(ISURF)+(A3LM(ISURF)+A6LM(ISURF)*ZPLT)*ZPLT
        IF (TRCPLT) THEN
          WRITE (6,*) 'COEFFICIENTS OF CURVE-EQUATION'
          WRITE (6,*) 'AK*X**2+BK*XY+CK*Y**2+DK*X+EK*Y+FK=0.    '
          WRITE (6,6662)  AK,BK,CK,DK,EK,FK
        ENDIF
C
C   DIE KURVENGLEICHUNG LAUTET
C   AK*X**2+BK*XY+CK*Y**2+DK*X+EK*Y+FK=0.    '
C   TRANSFORMIERE BK AUF 0., ABER ERHALTE INVARIANTEN
C
        DKCK=DK*CK
        DKCK2=DKCK+DKCK
        AKEK=AK*EK
        AKEK2=AKEK+AKEK
C  ALTE INVARIANTEN:
        S=AK+CK
        DEL=AK*CK-BK*BK*0.25
        DELS=DEL
        IF (ABS(DEL).LE.EPS10) DEL=0.
        DET=FK*DEL+(DK*(BK*EK-DKCK2)-EK*(AKEK2-DK*BK))*0.125
        DETS=DET
        IF (ABS(DET).LE.EPS10) DET=0.
        IF (TRCPLT) THEN
          WRITE (6,*) 'S,DELS,DETS,DEL,DET '
          WRITE (6,*) S,DELS,DETS,DEL,DET
        ENDIF
C
        ALF=0.
        IF (BK.EQ.0.) THEN
          COSA=1.
          SINA=0.
          TANA=0.
          A=AK
          C=CK
          D=DK
          E=EK
          F=FK
        ELSE
          IF (AK.NE.CK) ALF=ATAN(BK/(AK-CK))*0.5
          IF (AK.EQ.CK) ALF=45.*PIA/180.
          SINA=SIN(ALF)
          IF (SIGN(1.D0,SINA).NE.SIGN(1.D0,BK)) THEN
            ALF=ALF+PIA
            SINA=SIN(ALF)
          ENDIF
          COSA=COS(ALF)
          TANA=TAN(ALF)
          SINAQ=SINA*SINA
          COSAQ=1.-SINAQ
          SCA=SINA*COSA
          A=AK*COSAQ+BK*SCA+CK*SINAQ
          IF (ABS(A).LT.EPS10) A=0.
          C=S-A
          IF (ABS(C).LT.EPS10) C=0.
          F=FK
          IF (ABS(F).LT.EPS10) F=0.
          D=DK*COSA+EK*SINA
          IF (ABS(D).LT.EPS10) D=0.
          E=EK*COSA-DK*SINA
          IF (ABS(E).LT.EPS10) E=0.
          DH=D*0.5
          EH=E*0.5
          DET=DETER(A,0.D0,DH,0.D0,C,EH,DH,EH,F)
          IF (ABS(DET).LE.EPS10) DET=0.
        ENDIF
C
        IF (TRCPLT) THEN
          WRITE (6,*) 'AFTER TRANSFORMATION: B=0., NEW COEFFICIENTS'
          WRITE (6,6661) A,C,D,E,F
          WRITE (6,*) 'DETERMINANT DET= ',DET
        ENDIF
C
        IF (RLB(ISURF).EQ.1.) THEN
          XL1=MAX(XCHL,XLIMS1(ISURF,1))
          XL1=MIN(XCHR,XL1)
          XL2=MAX(XCHL,XLIMS2(ISURF,1))
          XL2=MIN(XCHR,XL2)
          YL1=MAX(YCHL,YLIMS1(ISURF,1))
          YL1=MIN(YCHR,YL1)
          YL2=MAX(YCHL,YLIMS2(ISURF,1))
          YL2=MIN(YCHR,YL2)
        ELSE
          XL1=XCHL
          XL2=XCHR
          YL1=YCHL
          YL2=YCHR
        ENDIF
C
        X13=YL1*SINA+XL1*COSA
        X14=YL2*SINA+XL1*COSA
        X23=YL1*SINA+XL2*COSA
        X24=YL2*SINA+XL2*COSA
C
        XANF=MIN(X13,X14,X23,X24)
        XEND=MAX(X13,X14,X23,X24)
C
        Y13=-XL1*SINA+YL1*COSA
        Y14=-XL2*SINA+YL1*COSA
        Y23=-XL1*SINA+YL2*COSA
        Y24=-XL2*SINA+YL2*COSA
C
        YANF=MIN(Y13,Y14,Y23,Y24)
        YEND=MAX(Y13,Y14,Y23,Y24)
C
        IF (RLB(ISURF).EQ.1.5) THEN
          XC1=XLIMS1(ISURF,1)
          XC2=XLIMS2(ISURF,1)
          YC1=YLIMS1(ISURF,1)
          YC2=YLIMS2(ISURF,1)
        ENDIF
C
        IF (RLB(ISURF).LT.0.) THEN
           DO 401 I=1,ILIN(ISURF)
              ALIN(I)=ALIMS(ISURF,I)
              XLIN(I)=XLIMS(ISURF,I)
              YLIN(I)=YLIMS(ISURF,I)
              ZLIN(I)=ZLIMS(ISURF,I)
401        CONTINUE
           DO 402 I=1,ISCN(ISURF)
              A0S(I)=ALIMS0(ISURF,I)
              A1S(I)=XLIMS1(ISURF,I)
              A2S(I)=YLIMS1(ISURF,I)
              A3S(I)=ZLIMS1(ISURF,I)
              A4S(I)=XLIMS2(ISURF,I)
              A5S(I)=YLIMS2(ISURF,I)
              A6S(I)=ZLIMS2(ISURF,I)
              A7S(I)=XLIMS3(ISURF,I)
              A8S(I)=YLIMS3(ISURF,I)
              A9S(I)=ZLIMS3(ISURF,I)
402        CONTINUE
           MLIN=ILIN(ISURF)
           MSCN=ISCN(ISURF)
        ENDIF
        IF (TRCPLT) THEN
          WRITE (6,*) 'PLOT REGION:'
          WRITE (6,*) 'XL1,XL2,YL1,YL2 ',XL1,XL2,YL1,YL2
          WRITE (6,*) 'PLOT REGION AFTER TRANSFORMATION'
          WRITE (6,*) 'XANF,XEND,YANF,YEND ',XANF,XEND,YANF,YEND
        ENDIF
C
C
C     FALLUNTERSCHEIDUNGEN
C
        IF (ABS(DEL).LE.EPS10) GOTO 1000
        XN=-DET/DEL
        IF (ABS(XN).LE.EPS10) XN=0.
        IF (XN) 405,500,400
405     IF (DEL.LT.-EPS10) GOTO 450
        IF (A.LT.-EPS10.AND.C.LT.-EPS10) GOTO 410
        GOTO 407
400     IF (DEL.LT.-EPS10) GOTO 450
        IF (A.GT.EPS10.AND.C.GT.EPS10) GOTO 410
C
C     KEINE REELLE LOESUNG
C
407     IF (TRCPLT) WRITE (6,6664)
        GOTO 200
C
C     ELLIPSE
C
410     AHALB=SQRT(ABS(XN/A))
        BHALB=SQRT(ABS(XN/C))
        XM=-D/(2.*A)
        YM=-E/(2.*C)
        IF (RLB(ISURF).EQ.1.5) THEN
           XA=-AHALB
           XE=AHALB
        ELSE
C   PLOTGEBIET ENTHAELT GEDREHTE BEGRENZUNGSBOX
C   BEGRENZUNGSBOX WIRD ALS PLOTGRENZE NACH RUECKTRANSFORMATION GENUTZT
           XA=MAX(XM-AHALB,XANF+1.D-8)-XM
           XE=MIN(XM+AHALB,XEND-1.D-8)-XM
        ENDIF
        CALL PLTIN (RLB(ISURF),ELLO,XA,XE,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'ELLO'
        CALL PLTIN (RLB(ISURF),ELLU,XA,XE,INN2,LBOX)
        IF (TRCPLT) WRITE (6,*) 'ELLU'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (6,6665)
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (6,6665)
        GOTO 200
C
C     HYPERBEL
C
450     CONTINUE
        RAD=MAX(0.D0,XN/A)
        X1=-D/A/2.+SQRT(RAD)
        X2=-D/A/2.-SQRT(RAD)
        INN1=0
        INN2=0
        INN3=0
        INN4=0
        IF (X1.GT.XL2) GOTO  451
        XA=MAX(X1,XANF)+1.D-6
        XE=XEND
        CALL PLTIN (RLB(ISURF),HYPP,XA,XE,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'HYPP'
        CALL PLTIN (RLB(ISURF),HYPM,XA,XE,INN2,LBOX)
        IF (TRCPLT) WRITE (6,*) 'HYPM'
        IF (INN1+INN2.EQ.0.AND.TRCPLT) WRITE (6,6666)
451     IF (X2.LT.XL1) GOTO  452
        XA=XANF
        XE=MIN(X2,XEND)-1.D-6
        CALL PLTIN (RLB(ISURF),HYPP,XA,XE,INN3,LBOX)
        IF (TRCPLT) WRITE (6,*) 'HYPP'
        CALL PLTIN (RLB(ISURF),HYPM,XA,XE,INN4,LBOX)
        IF (TRCPLT) WRITE (6,*) 'HYPM'
452     IF (INN3+INN4.EQ.0.AND.TRCPLT) WRITE (6,6666)
        GOTO 200
C
C
C     EIN PUNKT
C
500     CONTINUE
        IF (DEL.LT.-EPS10) GOTO 510
        X=-D/(2.*A)
        Y=-E/(2.*C)
        XP=XTRAN(X,Y)
        YP=YTRAN(X,Y)
        IF (TRCPLT) WRITE (6,*) 'POINT'
        IF (XP.LT.XL1.OR.XP.GT.XL2.OR.YP.LT.YL1.OR.YP.GT.YL2) THEN
          IF (TRCPLT) WRITE (6,6667)
          GOTO 200
        ENDIF
        IF (LZR) THEN
          CALL GRCHRC (0.1,0.0,16)
          CALL GRJMPS (sngl(XP),sngl(YP),3)
          CALL GRCHRC (0.3,0.0,16)
        ENDIF
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,0)
        GOTO 200
C
C     ZWEI SICH SCHNEIDENDE GERADEN
C
510     CONTINUE
        A0=-SQRT(ABS(A/C))
        A1=-SQRT(ABS(A/C))*D/(2.*A)-E/(2.*C)
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADY'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (6,6670)
        A0=SQRT(ABS(A/C))
        A1=SQRT(ABS(A/C))*D/(2.*A)-E/(2.*C)
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN2,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADY'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (6,6670)
        GOTO 200
C
1000    CONTINUE
        IF (ABS(A).GT.EPS10) GOTO 1300
        IF (ABS(C).GT.EPS10) GOTO 1400
        IF (ABS(E).GT.EPS10.AND.ABS(E).GE.ABS(D)) GOTO 1200
        IF (ABS(D).GT.EPS10.AND.ABS(D).GE.ABS(E)) GOTO 1100
        IF (ABS(F).GT.EPS10) GOTO 1002
        IF (TRCPLT) THEN
        WRITE (6,*) 'GLEICHUNG DER FORM 0.=0.'
        ENDIF
        GOTO 200
1002    IF (TRCPLT) WRITE (6,1001) F
1001    FORMAT (//1X,'KEIN PLOTT, DENN GLEICHUNG DER FORM F=',
     .            1PE12.4,' = 0')
        GOTO 200
C
C     GERADE;  D*X + E*Y + F = 0, D UNGLEICH 0., E=0. MOEGLICH
C
1100    CONTINUE
        A0=-E/D
        A1=-F/D
        CALL PLTIN (RLB(ISURF),GERADX,YANF,YEND,INN,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADX'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (6,6669)
        GOTO 200
C
C     GERADE; D*X + E*Y + F = 0, E UNGLEICH 0., D=0. MOEGLICH
C
1200    CONTINUE
        A0=-D/E
        A1=-F/E
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADY'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (6,6669)
        GOTO 200
C
C     PARABEL; A*X**2 + D*X + E*Y + F = 0
C
1300    CONTINUE
        IF (ABS(DET).LE.EPS10) GOTO 1310
        CALL PLTIN (RLB(ISURF),PARA1,XANF,XEND,INN,LBOX)
        IF (TRCPLT) WRITE (6,*) 'PARA1'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (6,6668)
        GOTO 200
C
C     PAAR PARALLELER GERADEN; A*X**2 + D*X + F = 0
C
1310    SQ=(D*D/(4.*A)-F)/A
        SQR=0.
        IF (SQ.GT.0.) SQR=SQRT(SQ)
        IF (SQ.LT.0.) THEN
          IF (TRCPLT) WRITE (6,6664)
          GOTO 200
        ENDIF
        A0=0.
        A1=-D/(2.*A)+SQR
        CALL PLTIN (RLB(ISURF),GERADX,YANF,YEND,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADX'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (6,6670)
        A1=-D/(2.*A)-SQR
        CALL PLTIN (RLB(ISURF),GERADX,YANF,YEND,INN2,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADX'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (6,6670)
        GOTO 200
C
C     PARABEL; C*Y**2 + D*X + E*Y + F = 0
C
1400    CONTINUE
        IF (ABS(DET).LE.EPS10) GOTO 1500
        XS=(E*E/(4.*C)-F)/D
        YS=PARA2O(XS+1.)
        IF (YS.LT.1.D50) THEN
           XAN=MIN(MAX(XANF,XS),XEND)
           XEN=XEND
        ELSE
           XEN=MAX(MIN(XEND,XS),XANF)
           XAN=XANF
        ENDIF
        IF (ABS(XEN-XAN).LT.EPS10) THEN
           IF (TRCPLT) WRITE (6,6668)
           GOTO 200
        ENDIF
        CALL PLTIN (RLB(ISURF),PARA2O,XAN,XEN,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'PARA2O'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (6,6671)
        CALL PLTIN (RLB(ISURF),PARA2U,XAN,XEN,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'PARA2U'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (6,6671)
        GOTO 200
C
C     PAAR PARALLELER GERADEN; C*Y**2 + E*Y + F = 0
C
1500    CONTINUE
        A0=0.
        SQ=(E*E/(4.*C)-F)/C
        SQR=0
        IF (SQ.GT.0.) SQR=SQRT(SQ)
        IF (SQ.LT.0.) THEN
           IF (TRCPLT) WRITE (6,6664)
           GOTO 200
        ENDIF
        A1=-E/(2.*C)+SQR
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN1,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADY'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (6,6670)
        A1=-E/(2.*C)-SQR
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN2,LBOX)
        IF (TRCPLT) WRITE (6,*) 'GERADY'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (6,6670)
C
C  SURFACE NO ISURF DONE
C
200   CONTINUE
C
C  END OF DO LOOP OVER SURFACE NUMBERS ISURF
C
      IF (LZR) THEN
        CALL GRDSH(1.,0.,1.)
        CALL GRNWPN(1)
      ENDIF
C
      IF (CH2Z0.NE.0.) CALL ZSHADD(CH2Z0,MANF,MEND)
      IF (ICUT.EQ.1) THEN
        CALL ROTADD(AFFI,AFF,MANF,MEND)
      ELSEIF (ICUT.EQ.2) THEN
        CALL ROTADD(AFFI,AFF,MANF,MEND)
      ELSEIF (ICUT.EQ.3) THEN
C  NO INVERS ROTATION NEEDED
      ENDIF
C
2000  CONTINUE
C
      IF (TRCPLT) WRITE (6,*) 'INSTOR= ',INSTOR
C
      IF (PLNUMS.OR.PLARR) THEN
        DO 3000 I=1,INUMS
          IPA=NUMSUR(I,2)
          IPE=NUMSUR(I,3)
          ARC=0.
          DO 3001 IARC=IPA+1,IPE
            ARC=ARC+SQRT((XPL2D(IARC)-XPL2D(IARC-1))**2+
     .                   (YPL2D(IARC)-YPL2D(IARC-1))**2)
3001      CONTINUE
          ARC05=ARC*0.5
          ARC=0.
          DO 3002 IARC=IPA+1,IPE
            ARC=ARC+SQRT((XPL2D(IARC)-XPL2D(IARC-1))**2+
     .                   (YPL2D(IARC)-YPL2D(IARC-1))**2)
            IF (ARC.GT.ARC05) GOTO 3003
3002      CONTINUE
C  POINT BETWEEN IARC-1 AND IARC
3003      CONTINUE
          IF (IARC.GT.1) THEN
            XINI=(XPL2D(IARC-1)+XPL2D(IARC))*0.5
            YINI=(YPL2D(IARC-1)+YPL2D(IARC))*0.5
          ENDIF
C
C  PLOT ARROWS: SURFACE NORMAL
          IF (PLARR) THEN
c           xlst=
c           ylst=
c           alen=
c           awid=
c           icode=
c           call grarrw(sngl(XINI),sngl(YINI),xlst,ylst,
c    .                  sngl(alen),sngl(awid),isngl(icode))
c
          ENDIF
C  PLOT SURFACE NUMBERS
          IF (PLNUMS) THEN
            IF (NUMSUR(I,1).LT.10) THEN
              WRITE (CH1,'(I1)') NUMSUR(I,1)
              CALL GRTXT (sngl(XINI),sngl(YINI),1,CH1)
            ELSEIF (NUMSUR(I,1).LT.100) THEN
              WRITE (CH2,'(I2)') NUMSUR(I,1)
              CALL GRTXT (sngl(XINI),sngl(YINI),2,CH2)
            ELSEIF (NUMSUR(I,1).LT.1000) THEN
              WRITE (CH3,'(I3)') NUMSUR(I,1)
              CALL GRTXT (sngl(XINI),sngl(YINI),3,CH3)
            ENDIF
          ENDIF
3000    CONTINUE
      ENDIF
C
      RETURN
6660  FORMAT (//1X,' CLOSED POLYGON OUTSIDE PLOTREGION')
6661  FORMAT (/1X,' A,C,D,E,F',/1X,1P,5E12.4)
6662  FORMAT (/1X,' AK,BK,CK,DK,EK,FK',/1X,1P,6E12.4)
6664  FORMAT (//1X,'NO REELL SOLUTION')
6665  FORMAT (//1X,'HALF ELLIPSE OUTSIDE PLOTREGION')
6666  FORMAT (//1X,'HYPERBEL OUTSIDE PLOTREGION')
6667  FORMAT (//1X,'SINGLE POINT OUTSIDE PLOTREGION')
6668  FORMAT (//1X,'PARABEL OUTSIDE PLOTREGION')
6669  FORMAT (//1X,'STRAIGHT LINE OUTSIDE PLOTREGION')
6670  FORMAT (//1X,'ONE OF THE TWO STRAIGHT LINES OUTSIDE PLOTREGION')
6671  FORMAT (//1X,'HALF PARABEL OUTSIDE PLOTREGION')
      END
C
      FUNCTION PARA1 (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      PARA1=-A/E*X*X-D/E*X-F/E
      RETURN
      END
C
C
      FUNCTION PARA2O (XIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      DH=-E/(2.*C)
      RAD=(E*E)/(4.*C*C)-(D*XIN+F)/C
      IF (RAD.GE.0.) GOTO 1
      PARA2O=1.D50
      RETURN
1     Y=SQRT(RAD)
      PARA2O=Y-DH
      RETURN
      END
C
C
      FUNCTION PARA2U (XIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      DH=-E/(2.*C)
      RAD=(E*E)/(4.*C*C)-(D*XIN+F)/C
      IF (RAD.GE.0.) GOTO 1
      PARA2U=1.D50
      RETURN
1     Y=-SQRT(RAD)
      PARA2U=Y-DH
      RETURN
      END
C
C
      FUNCTION ELLO (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      ELLO=BHALB*SQRT(MAX(0.D0,1.D0-X*X/(AHALB*AHALB)))
      RETURN
      END
C
C
      FUNCTION ELLU (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      ELLU=-BHALB*SQRT(MAX(0.D0,1.D0-X*X/(AHALB*AHALB)))
      RETURN
      END
C
C
      FUNCTION HYPP(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      AH=A*X*X
      DH=D*X
      RAD=E*E/(4.*C*C)-(F+AH+DH)/C
      IF (RAD.GE.0.) GOTO 1
      HYPP=1.D50
      RETURN
1     HYPP=-E/(2.*C)+SQRT(RAD)
      RETURN
      END
C
C
      FUNCTION HYPM(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      AH=A*X*X
      DH=D*X
      RAD=E*E/(4.*C*C)-(F+AH+DH)/C
      IF (RAD.GE.0.) GOTO 1
      HYPM=1.D50
      RETURN
1     HYPM=-E/(2.*C)-SQRT(RAD)
      RETURN
      END
C
C
      FUNCTION GERADY(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      GERADY=A0*X+A1
      RETURN
      END
C
C
      FUNCTION GERADX(Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CLMSUR'
      GERADX=Y
      Y=A0*Y+A1
      RETURN
      END
C
      SUBROUTINE PLTIN (RIB,FCN,ANF,END,INN,LBOX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CTRCEI'
      INCLUDE 'CLMSUR'
      LOGICAL LBOX
      EXTERNAL FCN
      IF (RIB.EQ.1..OR.RIB.EQ.0.)
     .                CALL PLTKI (FCN,ANF,END,INN,TRCPLT,LBOX)
      IF (RIB.EQ.1.5) CALL PLTKA (FCN,ANF,END,INN,TRCPLT,LBOX)
      IF (RIB.LT.0.)  CALL PLTKU (FCN,ANF,END,INN,TRCPLT,LBOX)
      RETURN
      END
C
C
C
C
      SUBROUTINE PLTKI (FCN,XANF,XEND,INN,TRCPLT,LBOX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CLMSUR'
      INCLUDE 'CCONA'
      INCLUDE 'CRECH'
      INCLUDE 'CPLOT'
      LOGICAL TRCPLT,LBOX
      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      GERAY(XG)=(YY-YYO)/(XX-XXO)*XG+YY-XX*(YY-YYO)/(XX-XXO)
      GERAX(YG)=(YG-YY)*(XX-XXO)/(YY-YYO)+XX
      IF (TRCPLT) WRITE (6,*) 'PLTKI'
      IF (LBOX) THEN
        CALL GRJMP(sngl(XL1),sngl(YL1))
        CALL GRDRW(sngl(XL2),sngl(YL1))
        CALL GRDRW(sngl(XL2),sngl(YL2))
        CALL GRDRW(sngl(XL1),sngl(YL2))
        CALL GRDRW(sngl(XL1),sngl(YL1))
      ENDIF
      INC=100
      XINC=DBLE(INC)
      DX=(XEND-XANF)/XINC
      IFLAG=0
      ISIDE=0
      INN=0
      XX=0.
      YY=0.
      DO 431 I=1,INC+1
         X=XANF+(I-1)*DX
         Y=FCN(X)
         XXO=XX
         YYO=YY
         XX=XTRAN(X+XM,Y+YM)
         YY=YTRAN(X+XM,Y+YM)
         IJUMP=ISIDE
         ISIDE=0
         IF (XX.LE.XL1+EPS10) ISIDE=3
         IF (XX.GE.XL2-EPS10) ISIDE=2
         IF (YY.LE.YL1+EPS10) ISIDE=4
         IF (YY.GE.YL2-EPS10) ISIDE=1
         IF (ISIDE.NE.0) THEN
             IF (IFLAG.EQ.0) GOTO 431
             IJUMP=ISIDE
             IF (IFLAG.NE.0) GOTO 432
         ENDIF
C        ISIDE=0
         IF (IFLAG.EQ.0.AND.I.NE.1) GOTO 432
439      IF (I.NE.1) IFLAG=24
         IF (IFLAG.EQ.0) THEN
            IF (LZR) CALL GRJMP (sngl(XX),sngl(YY))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,0)
         ENDIF
         IF (IFLAG.NE.0) THEN
            IF (LZR) CALL GRDRW (sngl(XX),sngl(YY))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,1)
         ENDIF
         INN=MAX0(INN,IFLAG)
         GOTO 431
432      GOTO (451,452,453,454),IJUMP
         GOTO 439
451      XP=GERAX(YL2)
         YP=YL2
         GOTO 429
452      XP=XL2
CPB      YP=YTRAN(XP,FCN(XP-XM)+YM)
         YP=GERAY(XL2)
         GOTO 429
453      XP=XL1
CPB      YP=YTRAN(XP,FCN(XP-XM)+YM)
         YP=GERAY(XL1)
         GOTO 429
454      XP=GERAX(YL1)
         YP=YL1
         GOTO 429
C
429      CONTINUE
         IF (IFLAG.EQ.0) THEN
            IF (LZR) CALL GRJMP (sngl(XP),sngl(YP))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,0)
         ENDIF
         IF (IFLAG.NE.0) THEN
            IF (LZR) CALL GRDRW (sngl(XP),sngl(YP))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
         ENDIF
         INN=MAX0(INN,IFLAG)
         IF (IFLAG.EQ.0) GOTO 439
         IFLAG=0
C
431   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE PLTKA (FCN,XANF,XEND,INN,TRCPLT,LBOX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CLMSUR'
      INCLUDE 'CRECH'
      INCLUDE 'CPLOT'
      LOGICAL TRCPLT,LBX,LBXO,LPA,L1,L2,L3,L4,LBOX
      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      GERAY(XG)=(YY-YYO)/(XX-XXO)*XG+YY-XX*(YY-YYO)/(XX-XXO)
      GERAX(YG)=(YG-YY)*(XX-XXO)/(YY-YYO)+XX
      IF (TRCPLT) WRITE (6,*) 'PLTKA'
      IF (LBOX) THEN
        CALL GRJMP(sngl(XL1),sngl(YL1))
        CALL GRDRW(sngl(XL2),sngl(YL1))
        CALL GRDRW(sngl(XL2),sngl(YL2))
        CALL GRDRW(sngl(XL1),sngl(YL2))
        CALL GRDRW(sngl(XL1),sngl(YL1))
      ENDIF
      INC=100
      XINC=DBLE(INC)
      DX=(XEND-XANF)/XINC
      IFLAG=0
      INN=0
      X=XANF
      Y=FCN(X)
      XX=XTRAN(X+XM,Y+YM)
      YY=YTRAN(X+XM,Y+YM)
C*****LPA=.T. IF POINT INSIDE PLOTAREA
      LPA=XX.GE.XL1.AND.XX.LE.XL2.AND.YY.GE.YL1.AND.YY.LE.YL2
C*****LBX=.T. IF POINT INSIDE BOX
      LBX=XX.GE.XC1.AND.XX.LE.XC2.AND.YY.GE.YC1.AND.YY.LE.YC2
      IF (LPA.AND..NOT.LBX) THEN
        IF (LZR) CALL GRJMP (sngl(XX),sngl(YY))
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,0)
        IFLAG=24
        INN=24
      ENDIF
      DO 1 I=1,INC
         X=XANF+I*DX
         Y=FCN(X)
         XXO=XX
         YYO=YY
         LBXO=LBX
         XX=XTRAN(X+XM,Y+YM)
         YY=YTRAN(X+XM,Y+YM)
         LPA=XX.GE.XL1.AND.XX.LE.XL2.AND.YY.GE.YL1.AND.YY.LE.YL2
         LBX=XX.GE.XC1.AND.XX.LE.XC2.AND.YY.GE.YC1.AND.YY.LE.YC2
         IF (IFLAG.EQ.0.AND.(LBX.OR..NOT.LPA)) GOTO 1
         IF (IFLAG.NE.0.AND.(LPA.AND..NOT.LBX)) THEN
            IF (LZR) CALL GRDRW (sngl(XX),sngl(YY))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,1)
         ELSE
            IF (LBX.OR.LBXO) THEN
               XB1=XC1
               XB2=XC2
               YB1=YC1
               YB2=YC2
            ELSE
               XB1=XL1
               XB2=XL2
               YB1=YL1
               YB2=YL2
            ENDIF
            L1=XX.LE.XB1.OR.XXO.LE.XB1
            L2=XX.GE.XB2.OR.XXO.GE.XB2
            L3=YY.LE.YB1.OR.YYO.LE.YB1
            L4=YY.GE.YB2.OR.YYO.GE.YB2
            IF (L1) XP=XB1
            IF (L2) XP=XB2
            IF (L3) YP=YB1
            IF (L4) YP=YB2
CPB   IM FALLE EINER "SCHIEFEN" ELLIPSE FEHLER MOEGLICH
            IF (L1.OR.L2) YP=FCN(XP-XM)+YM
            IF (L3.OR.L4) XP=GERAX(YP)
            IF (IFLAG.EQ.0) THEN
               IF (LZR) THEN
                  CALL GRJMP (sngl(XP),sngl(YP))
                  CALL GRDRW (sngl(XX),sngl(YY))
               ENDIF
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                  CALL STCOOR (XP,YP,0)
                  CALL STCOOR (XX,YY,1)
               ENDIF
               IFLAG=24
               INN=24
            ELSE
               IF (LZR) CALL GRDRW (sngl(XP),sngl(YP))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
               IFLAG=0
            ENDIF
         ENDIF
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE PLTKU (FCN,XANF,XEND,INN,TRCPLT,LBOX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CLMSUR'
      INCLUDE 'CRECH'
      INCLUDE 'CPLOT'
      LOGICAL TRCPLT,LIN,LINO,LTST,LPLA,LPLAO,LPL,LGL,L1,L2,L3,L4
      LOGICAL LBOX
      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      IF (TRCPLT) WRITE (6,*) 'PLTKU'
      IF (LBOX) THEN
        CALL GRJMP(sngl(XL1),sngl(YL1))
        CALL GRDRW(sngl(XL2),sngl(YL1))
        CALL GRDRW(sngl(XL2),sngl(YL2))
        CALL GRDRW(sngl(XL1),sngl(YL2))
        CALL GRDRW(sngl(XL1),sngl(YL1))
      ENDIF
C
      INC=100
100   CONTINUE
      XINC=DBLE(INC)
      DX=(XEND-XANF)/XINC
      INN=0
      IFLAG=0
C
      DO 1 I=1,MLIN
1        ALIN(I)=ALIN(I)+ZLIN(I)*ZPLT
      DO 2 I=1,MSCN
         A0S(I)=A0S(I)+(A3S(I)+A6S(I)*ZPLT)*ZPLT
         A1S(I)=A1S(I)+A8S(I)*ZPLT
         A2S(I)=A2S(I)+A9S(I)*ZPLT
2     CONTINUE
      DO 10 I=1,INC+1
         X=XANF+(I-1)*DX
         Y=FCN(X)
         XX=XTRAN(X+XM,Y+YM)
         YY=YTRAN(X+XM,Y+YM)
         LPLA=XX.GE.XL1.AND.XX.LE.XL2.AND.YY.GE.YL1.AND.YY.LE.YL2
         LIN=.TRUE.
         DO 11 J=1,MLIN
            GL=ALIN(J)+XLIN(J)*XX+YLIN(J)*YY
            LIN=LIN.AND.GL.LE.0.
11       CONTINUE
         DO 12 J=1,MSCN
            GS=A0S(J)+(A1S(J)+A4S(J)*XX)*XX+(A2S(J)+
     .         A5S(J)*YY)*YY+A7S(J)*XX*YY
            LIN=LIN.AND.GS.LE.0.
12       CONTINUE
C
         IF (I.EQ.1) THEN
C*****FIRST POINT OF LINE
            IF (LIN.AND.LPLA) THEN
               IF (LZR) CALL GRJMP (sngl(XX),sngl(YY))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,0)
               INN=24
               IFLAG=24
            ENDIF
         ELSE
C*****ALL OTHER POINTS
C*****POINT OUT OF PLOT AREA
            IF (.NOT.LPLA.AND..NOT.LPLAO) GOTO 16
C*****POINT OUT OF CONFIGURATION
            IF (.NOT.LIN.AND..NOT.LINO) GOTO 16
C*****POINT INSIDE PLOTAREA AND INSIDE THE CONFIGURATION
            IF (LIN.AND.LINO.AND.LPLA.AND.LPLAO) THEN
               IF (LZR) CALL GRDRW (sngl(XX),sngl(YY))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,1)
            ELSE
C*****LINE CROSSES BOUNDARY
C***** ... OF CONFIGURATION
               LGL=LIN.AND.LINO
C***** ... OF PLOT AREA
               LPL=LPLA.AND.LPLAO
               IF (.NOT.LPL) THEN
                  L1=XX.LE.XL1.OR.XXO.LE.XL1
                  L2=XX.GE.XL2.OR.XXO.GE.XL2
                  L3=YY.LE.YL1.OR.YYO.LE.YL1
                  L4=YY.GE.YL2.OR.YYO.GE.YL2
                  IF (LPLA) THEN
                     XOK=XX
                     YOK=YY
                     XOUT=XXO
                     YOUT=YYO
                  ELSE
                     XOK=XXO
                     YOK=YYO
                     XOUT=XX
                     YOUT=YY
                  ENDIF
               ENDIF
C*****INTERPOLATE POINT ON CONFIGURATION BOUNDARY
               IF (.NOT.LGL) THEN
                  IF (LIN) THEN
                     XS=XANF+(I-1)*DX
                     DXX=-DX
                  ELSE
                    XS=XANF+(I-2)*DX
                    DXX=DX
                  ENDIF
C
13                DXX=DXX*0.5
                  XT=XS+DXX
                  YT=FCN(XT)
                  XTT=XTRAN(XT+XM,YT+YM)
                  YTT=YTRAN(XT+XM,YT+YM)
                  LTST=.TRUE.
                  DO 14 J=1,MLIN
                     GL=ALIN(J)+XLIN(J)*XTT+YLIN(J)*YTT
                     LTST=LTST.AND.GL.LE.0.
14                CONTINUE
                  DO 15 J=1,MSCN
                     GS=A0S(J)+(A1S(J)+A4S(J)*XTT)*XTT+
     .                  (A2S(J)+A5S(J)*YTT)*YTT+A7S(J)*XTT*YTT
                     LTST=LTST.AND.GS.LE.0.
15                CONTINUE
                  IF (LTST) XS=XS+DXX
                  IF (ABS(DXX).GT.1.E-3) GOTO 13
C
                  YS=FCN(XS)
                  XP=XTRAN(XS+XM,YS+YM)
                  YP=YTRAN(XS+XM,YS+YM)
                  IF (XP.GE.XL1.AND.XP.LE.XL2.AND.
     .                YP.GE.YL1.AND.YP.LE.YL2) THEN
C*****INTERPOLATED POINT INSIDE PLOT AREA, CAN BE PLOTTED
                     IF (LIN) THEN
                        IF (LZR) THEN
                           CALL GRJMP (sngl(XP),sngl(YP))
                           CALL GRDRW (sngl(XX),sngl(YY))
                        ENDIF
                        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                           CALL STCOOR (XP,YP,0)
                           CALL STCOOR (XX,YY,1)
                        ENDIF
                        INN=24
                        IFLAG=24
                     ELSE
                        IF (LZR) CALL GRDRW (sngl(XP),sngl(YP))
                        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                           CALL STCOOR (XP,YP,1)
                        ENDIF
                        IFLAG=0
                     ENDIF
                  ELSE
C*****INTERPOLATED POINT OUTSIDE PLOTAREA
                     L1=XP.LE.XL1
                     L2=XP.GE.XL2
                     L3=YP.LE.YL1
                     L4=YP.GE.YL2
                     XOUT=XP
                     YOUT=YP
                     IF (LPL) THEN
                        IF (LIN) THEN
                           XOK=XX
                           YOK=YY
                        ELSE
                           XOK=XXO
                           YOK=YYO
                        ENDIF
                     ENDIF
                     LPL=.FALSE.
                  ENDIF
               ENDIF
               IF (.NOT.LPL) THEN
C*****INTERPOLATE POINT ON BOUNDARY OF PLOTAREA
                  IF (L1) XP=XL1
                  IF (L2) XP=XL2
                  IF (L3) YP=YL1
                  IF (L4) YP=YL2
                  IF (L1.OR.L2) YP=YOK-(XOK-XP)*(YOK-YOUT)/(XOK-XOUT)
                  IF (L3.OR.L4) XP=XOK-(YOK-YP)*(XOK-XOUT)/(YOK-YOUT)
                  IF (IFLAG.EQ.0) THEN
                     IF (LZR) THEN
                        CALL GRJMP (sngl(XP),sngl(YP))
                        CALL GRDRW (sngl(XX),sngl(YY))
                     ENDIF
                     IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                        CALL STCOOR (XP,YP,0)
                        CALL STCOOR (XX,YY,1)
                     ENDIF
                     INN=24
                     IFLAG=24
                  ELSE
                     IF (LZR) CALL GRDRW (sngl(XP),sngl(YP))
                     IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                       CALL STCOOR (XP,YP,1)
                     ENDIF
                     IFLAG=0
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
16       XXO=XX
         YYO=YY
         LINO=LIN
         LPLAO=LPLA
10    CONTINUE
      IF (INN.GT.0.OR.INC.GE.400) RETURN
      INC=INC*2
      GOTO 100
      END
C
      FUNCTION FMU (V1,V2,V3,VR1,VR2,VR3,P1,P2,P3,PR1,PR2,PR3,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      XNEN=PR1*VR2-PR2*VR1
      IF (ABS(XNEN).GT.EPS) THEN
        FMU=(VR2*(V1-P1)-VR1*(V2-P2))/XNEN
        RETURN
      ENDIF
      XNEN=PR1*VR3-PR3*VR1
      IF (ABS(XNEN).GT.EPS) THEN
        FMU=(VR3*(V1-P1)-VR1*(V3-P3))/XNEN
        RETURN
      ENDIF
      XNEN=PR2*VR3-PR3*VR2
      IF (ABS(XNEN).GT.EPS) THEN
        FMU=(VR3*(V2-P2)-VR2*(V3-P3))/XNEN
        RETURN
      ENDIF
C     WRITE (6,*) ' SCHNITTGERADE PARALLEL ZU VIERECKSEITE '
      FMU=10.
      RETURN
      END
C
      SUBROUTINE GSP (P1,P2,R1,R2,Q1,Q2,S1,S2,XLA,XMU,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      S3=0.
      R3=0.
      PK1=R2*S3-R3*S2
      PK2=R3*S1-R1*S3
      PK3=R1*S2-R2*S1
      XNO=SQRT(PK1*PK1+PK2*PK2+PK3*PK3)
      IF (XNO.LT.EPS) THEN
       XLA=10.
       XMU=10.
      ELSE IF (ABS(R1).LT.EPS) THEN
       XLA=(P1-Q1)/S1
       XMU=-(P2-Q2-S2*XLA)/R2
      ELSE
       XLA=(P2-Q2)/S2
       XMU=-(P1-Q1-S1*XLA)/R1
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSIDE (G,R,P11,P12,P21,P22,L1,L2,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CRECH'
      INCLUDE 'CPLOT'
      DIMENSION G(4,2),R(4,2)
      LOGICAL L1,L2
      IS=0
      DO 100 I=1,4
      CALL GSP (G(I,1),G(I,2),R(I,1),R(I,2),P11,P12,P21-P11,P22-P12,
     .          XLA,XMU,EPS)
      IF (XLA.GE.0..AND.XLA.LE.1..AND.XMU.GE.0..AND.XMU.LE.1.) THEN
        XP=G(I,1)+XMU*R(I,1)
        YP=G(I,2)+XMU*R(I,2)
        IF (L1) THEN
          IF (LZR) THEN
            CALL GRJMP (sngl(P11),sngl(P12))
            CALL GRDRW (sngl(XP),sngl(YP))
          ENDIF
          IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
            CALL STCOOR (P11,P12,0)
            CALL STCOOR (XP,YP,1)
          ENDIF
          RETURN
        ELSE IF (L2) THEN
          IF (LZR) THEN
            CALL GRJMP (sngl(XP),sngl(YP))
            CALL GRDRW (sngl(P21),sngl(P22))
          ENDIF
          IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
            CALL STCOOR (XP,YP,0)
            CALL STCOOR (P21,P22,1)
          ENDIF
          RETURN
        ELSE
          IS=IS+1
          IF (IS.EQ.1) THEN
            IF (LZR) CALL GRJMP (sngl(XP),sngl(YP))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,0)
          ELSEIF (IS.EQ.2) THEN
            IF (LZR) CALL GRDRW (sngl(XP),sngl(YP))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
            RETURN
          ENDIF
        ENDIF
      ENDIF
100   CONTINUE
      RETURN
      END
C
      SUBROUTINE GRNXTB(K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  K=1 CALL VOR EINEM BILD AUS EIGENER GRSOFTWARE
C  K=2 CALL VOR EINEM BILD MIT GRBLD
      DATA IFIRST/0/,IFRST2/0/
      WRITE (6,*) 'GRNXTB K,IFIRST,IFRST2 ',K,IFIRST,IFRST2
      GOTO (1,2),K
1     IF (IFIRST.EQ.1) THEN
        CALL GRNXTF
        RETURN
      ELSE
        IFIRST=1
      ENDIF
      RETURN
C
2     IF (IFIRST.EQ.0) THEN
        CALL GRSCLC (0.,0.0,39.5,28.7)
        IFIRST=1
        IFRST2=1
      ELSE IF (IFRST2.EQ.0) THEN
        CALL GRNXTF
        CALL GRSCLC (0.,0.,39.5,28.7)
        IFRST2=1
      ELSE
        CALL GRSCLC (3.,3.5,39.5,28.7)
      ENDIF
      RETURN
      END
C
      SUBROUTINE EXIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL GREND
      STOP
      END
C
      SUBROUTINE PLT3D (XR,YR,FAKX,FAKY,ITH,ABSMIN,ABSMAX,ORDMIN,ORDMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGEOM'
      INCLUDE 'CTRCEI'
      INCLUDE 'COMSOU'
      INCLUDE 'CGRID'
      INCLUDE 'CINIT'
      INCLUDE 'CADGEO'
      INCLUDE 'CPL3D'
      INCLUDE 'CLGIN'
      INCLUDE 'CCONA'
      INCLUDE 'CLOGAU'
      INCLUDE 'CTEXT'
      INCLUDE 'CUPD'
      INCLUDE 'CPOLYG'
      INCLUDE 'CPLOT'
      INCLUDE 'COMUSR'
      INCLUDE 'CRECH'
      PARAMETER (NPLY=501)
      DIMENSION AL(10),AR(10),XP(NPLY),YP(NPLY),ZPLOT(N3RD+NTOR),
     .          XSAVE(NPLY,N3RD+NTOR),YSAVE(NPLY,N3RD+NTOR),
     .          PHIAN(9),PHIEN(9)
      DIMENSION ILT(N3RD+NTOR)
      REAL*4 XPS(NPLY),YPS(NPLY)
      DIMENSION XX(101),YY(101)
      LOGICAL PLABLE(NLIM),LPERID(NLIM),LSYMET(NLIM),
     .        LERR1,LERR2,LSAVE,PLT1,PLT2,PLT3
C
C  DEFAULT TOROIDAL GRID PLOT OPTION
C  IN CASE THAT NO TOROIDAL GRID IS DEFINED
C
      IF (.NOT.NLTOR.AND..NOT.NLTRA) THEN
        IPLTS(3)=1
        IPLAS(3,1)=1
        IPLES(3,1)=2
      ENDIF
C
      IF (NLTRA) THEN
        ITH=ANGLE3
        IF (ITH.LE.0.OR.ITH.GE.NTTRA) ITH=1
        WIN=ZZONE(ITH)
        RMT=RMTOR
        WINJ=WIN
      ELSEIF (NLTRZ) THEN
        WIN=0.
        RMT=0.
        WINJ=WIN
      ENDIF
C
      CALL ZEROA1(XP,NPLY)
      CALL ZEROA1(YP,NPLY)
      CALL ZEROL(PLABLE,1,NLIM)
C
      NZAD=5
      NINNE=6
      NIN=20
C
      XMI3D=CH3X0-CH3MX
      XMA3D=CH3X0+CH3MX
      YMI3D=CH3Y0-CH3MY
      YMA3D=CH3Y0+CH3MY
      ZMI3D=CH3Z0-CH3MZ
      ZMA3D=CH3Z0+CH3MZ
      ABSMAX=-1.E30
      ABSMIN=1.D30
      ORDMAX=-1.E30
      ORDMIN=1.D30
C
C  FIND MAXIMA AND MINIMA FOR SCALING OF 3D PLOT WINDOW
C
      LSAVE=PLBOX
      PLBOX=.FALSE.
      CALL PL3D (XMI3D,YMI3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMI3D,YMI3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMI3D,YMA3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMI3D,YMA3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMI3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMI3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMA3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMA3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
C
C
      XNULL=8.
      YNULL=2.
      XCENT=24.
      YCENT=24.
      DELX=ABSMAX-ABSMIN
      DELY=ORDMAX-ORDMIN
C  VERIFY THAT XCENT/YCENT=DELX/DELY
C  IN ORDER TO AVOID RESCALING DURING PLOTTING
C
      YN=XCENT*DELY/DELX
      XN=XCENT
      IF (YN.GT.YCENT) THEN
        XN=YCENT/YN*XCENT
        YN=YCENT
      ENDIF
C
C
      PLBOX=LSAVE
      CALL GRNXTB(1)
      CALL GRSCLC(sngl(XNULL),sngl(YNULL),sngl(XN+XNULL),
     .            sngl(YN+YNULL))
      CALL GRSCLV(sngl(ABSMIN),sngl(ORDMIN),sngl(ABSMAX),sngl(ORDMAX))
      FAKX=XN/(ABSMAX-ABSMIN)
      FAKY=YN/(ORDMAX-ORDMIN)
C
C  PLOT ADDITIONAL SURFACES
C
C  IF NLTRA, ASSUME THAT THIS SURFACE IS GIVEN IN LOCAL TOROIDAL SYSTEM
C            NO. ILTOR. PLOT IS DONE IN CO-ORDINATE SYSTEM OF CELL ITH
C            WHICH WAS SELECTED BY THE INPUT FLAG ANGLE3
C
      DO 100 I=1,5
        IF (.NOT.PL3A(I)) GOTO 100
        DO 10 IP=1,IPLTA(I)
        DO 10 J=IPLAA(I,IP),IPLEA(I,IP)
          IF (J.GT.NLIMI) GOTO 10
          IF (IGJUM0(J).NE.0) THEN
            IF (TRCPLT) THEN
              WRITE (6,*) 'SURFACE NO. ',J,' OUT'
            ENDIF
            GOTO 10
          ELSE
            IF (NLTRA.AND.ILTOR(J).LE.0) THEN
C  STILL TO BE WRITTEN: BETTER WAY OF IDENTIFYING TOROIDALLY SYMMETRIC
C                       SURFACE
              IF (P3(3,J).GT.1.D50.AND.ZLIMS1(J,1).LT.-1.D15.AND.
     .                                 ZLIMS2(J,1).GT.1.D15) THEN
                LSYMET(J)=.TRUE.
                IF (TRCPLT) THEN
                  WRITE (6,*) 'SURFACE NO. ',J,' TOROIDALLY SYMMETRIC'
                  WRITE (6,*) 'PLOT LATER INTO STANDARD MESH '
                ENDIF
                GOTO 10
              ELSE
                LPERID(J)=.TRUE.
                NJZ=0
                DO 50 IPZ=1,IPLTS(3)
                  IF (IPLAS(3,IPZ).LT.IPLES(3,IPZ)) THEN
                    DO 51 JJ=IPLAS(3,IPZ),IPLES(3,IPZ)-1
                      NJZ=NJZ+1
                      ILT(NJZ)=JJ
51                  CONTINUE
                  ELSE
                    DO 52 JJ=IPLAS(3,IPZ),NTTRA-1
                      NJZ=NJZ+1
                      ILT(NJZ)=JJ
52                  CONTINUE
                    DO 53 JJ=1,IPLES(3,IPZ)-1
                      NJZ=NJZ+1
                      ILT(NJZ)=JJ
53                  CONTINUE
                  ENDIF
50              CONTINUE
                IF (TRCPLT) THEN
                  WRITE (6,*) 'SURFACE NO. ',J,' TOROIDALLY PERODIC'
                  WRITE (6,*) 'PLOT ADD. SURFACE NO. ',J
                  WRITE (6,*) 'INTO ',NJZ, ' TOROIDAL SEGMENTS'
                ENDIF
              ENDIF
            ELSEIF (.NOT.NLTRA.OR.ILTOR(J).GT.0) THEN
              NJZ=1
              ILT(1)=ILTOR(J)
              LSYMET(J)=.FALSE.
              LPERID(J)=.FALSE.
              IF (TRCPLT) THEN
                WRITE (6,*) 'PLOT ADD. SURFACE NO. ',J
              ENDIF
            ENDIF
          ENDIF
C

          DO IJZ=1,NJZ
          IF (NLTRA) WINJ=ZZONE(ILT(IJZ))
C
C** 1 <= RLB < 2 ?
C
          IF (RLB(J).EQ.1..OR.RLB(J).EQ.1.5) THEN
C
C**GEKRUEMMTE FLAECHE ODER  EBENENPAAR ?
            IF (JUMLIM(J).EQ.0) THEN
              CALL FL2O (A0LM(J),A1LM(J),A2LM(J),
     .                   A3LM(J),A4LM(J),A5LM(J),
     .                   A6LM(J),A7LM(J),A8LM(J),
     .                   A9LM(J),MERK,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                   RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
              IF (TRCPLT) THEN
                WRITE (6,*) 'FL2O CALLED '
                WRITE (6,*) 'MERK= ',MERK
              ENDIF
C**ZYLINDER: FINDE ACHSE
              IF (MERK.EQ.4) THEN
                IF (TRCPLT) WRITE (6,*) 'CX,CY,CZ,RZYL ',CX,CY,CZ,RZYL
                IF (ABS(CX).GT.1.D-6.AND.
     .              MAX(ABS(CY),ABS(CZ)).LE.1.D-6) THEN
                  T1=(XLIMS1(J,1)-ZX0)/CX
                  T2=(XLIMS2(J,1)-ZX0)/CX
                ELSEIF (ABS(CY).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CZ)).LE.1.D-6) THEN
                  T1=(YLIMS1(J,1)-ZY0)/CY
                  T2=(YLIMS2(J,1)-ZY0)/CY
                ELSEIF (ABS(CZ).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CY)).LE.1.D-6) THEN
                  T1=(ZLIMS1(J,1)-ZZ0)/CZ
                  T2=(ZLIMS2(J,1)-ZZ0)/CZ
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
C  ZYLINDER: GGFLS MEHRERE TEILSTUECKE
                CALL CTQUA (A0LM(J),A1LM(J),A2LM(J),A3LM(J),A4LM(J),
     .                      A5LM(J),A6LM(J),A7LM(J),A8LM(J),A9LM(J),
     .                      XLIMS1(J,1),XLIMS2(J,1),YLIMS1(J,1),
     .                      YLIMS2(J,1),ZLIMS1(J,1),ZLIMS2(J,1),
     .                      RLB(J),RZYL,ZX0,ZY0,ZZ0,PHIAN,PHIEN,IPRT)
                IF (TRCPLT) THEN
                  WRITE (6,*) ' T1,T2 ',T1,T2
                  WRITE (6,*) ' IPRT ',IPRT
                  WRITE (6,*) ' PHIAN,PHIEN ',(PHIAN(JP),PHIEN(JP),
     .                                         JP=1,IPRT)
                ENDIF
                DO 109 IPR=1,IPRT
                  PHA=PHIAN(IPR)
                  PHE=PHIEN(IPR)
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                         RZYL,NZAD,NINNE,NIN,ILCOL(J),
     .                         IGFIL(J).NE.0,
     .                         J,0,AL,0,AR,PHA,PHE)
109             CONTINUE
C**KEGEL: BISLANG NUR EIN STUECK MOEGLICH. FINDE ACHSE
              ELSEIF (MERK.EQ.8) THEN
                IF (TRCPLT) WRITE (6,*) 'CX,CY,CZ,RZYL ',CX,CY,CZ,RZYL
                IF (ABS(CX).GT.1.D-6.AND.
     .              MAX(ABS(CY),ABS(CZ)).LE.1.D-6) THEN
                  T1=(XLIMS1(J,1)-ZX0)/CX
                  T2=(XLIMS2(J,1)-ZX0)/CX
                ELSEIF (ABS(CY).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CZ)).LE.1.D-6) THEN
                  T1=(YLIMS1(J,1)-ZY0)/CY
                  T2=(YLIMS2(J,1)-ZY0)/CY
                ELSEIF (ABS(CZ).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CY)).LE.1.D-6) THEN
                  T1=(ZLIMS1(J,1)-ZZ0)/CZ
                  T2=(ZLIMS2(J,1)-ZZ0)/CZ
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
                CALL CONE (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                     RZYL,NZAD,NINNE,NIN,ILCOL(J),
     .                     IGFIL(J).NE.0,J,0,AL,0,AR)
C**KUGEL,ELLIPSOID: BISLANG NUR EIN STUECK MOEGLICH
C             ELSEIF (MERK.EQ.8) CALL SPHERE (ZX0,ZY0,ZZ0,HX,HY,HZ,
C    .                                   NZAD,NINNE,NIN,ILCOL(J),
c    .                                   IGFIL(J).NE.0,J,0,AL,0,AR)
C**KUGEL, ELLIPSOID
              ELSEIF (MERK == 13) THEN
                CALL ELLIPSOID (ZX0,ZY0,ZZ0,CX,CY,CZ,XLIMS1(J,1),
     .               YLIMS1(J,1),ZLIMS1(J,1),XLIMS2(J,1),YLIMS2(J,1),
     .               ZLIMS2(J,1),RLB(J),ILCOL(J),5,5,5)
C**PAAR VON EBENEN (ODER EINE DOPPELEBENE)
              ELSEIF (MERK.EQ.1.OR.MERK.EQ.2.OR.MERK.EQ.3) THEN
                IF (TRCPLT) WRITE (6,*) 'B0,B1,B2,B3 ',B0,B1,B2,B3
                CALL PLANE (B0,B1,B2,B3,RLB(J),NLIM,EPS10,
     .                      ALIMS,XLIMS,YLIMS,ZLIMS,
     .                      ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                             XLIMS2,YLIMS2,ZLIMS2,
     .                             XLIMS3,YLIMS3,ZLIMS3,
     .                      ILCOL(J),IGFIL(J).NE.0,J)
                IF (MERK.NE.1) THEN
                  IF (TRCPLT) WRITE (6,*) 'F0,F1,F2,F3 ',F0,F1,F2,F3
                  CALL PLANE (F0,F1,F2,F3,RLB(J),NLIM,EPS10,
     .                        ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                ENDIF
              ELSE
                PLABLE(J)=.TRUE.
                GOTO 10
              ENDIF
C**EINE EBENE
            ELSEIF (JUMLIM(J).NE.0) THEN
              CALL PLANE (A0LM(J),A1LM(J),A2LM(J),A3LM(J),RLB(J),NLIM,
     .              EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                    ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                           XLIMS2,YLIMS2,ZLIMS2,
     .                           XLIMS3,YLIMS3,ZLIMS3,
     .                    ILCOL(J),IGFIL(J).NE.0,J)
            ENDIF
C
C**RLB >= 3 ? EIN EBENENSTUECK, DURCH POLYGON BEGRENZT
C
          ELSEIF (RLB(J).GT.2.) THEN
C
            CALL PRLLO(P1(1,J),P2(1,J),P3(1,J),P4(1,J),P5(1,J),
     .                 ILCOL(J),IGFIL(J).NE.0)
C
C**RLB < 0 ? ERST EINIGE OPTIONEN VORHANDEN, REST: CALL PLTUSR
C
          ELSEIF (RLB(J).LT.0) THEN
C
            IF (JUMLIM(J).EQ.0) THEN
              CALL FL2O (A0LM(J),A1LM(J),A2LM(J),
     .                   A3LM(J),A4LM(J),A5LM(J),
     .                   A6LM(J),A7LM(J),A8LM(J),
     .                   A9LM(J),MERK,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                   RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
              IF (TRCPLT) THEN
                WRITE (6,*) 'FL2O CALLED '
                WRITE (6,*) 'MERK= ',MERK
              ENDIF
C
C**PAAR VON EBENEN ODER DOPPELEBENE ?
              IF (MERK.LE.3) THEN
                IF (ISCN(J).EQ.0) THEN
                  CALL PLANE (B0,B1,B2,B3,RLB(J),
     .                        NLIM,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                  CALL PLANE (F0,F1,F2,F3,RLB(J),
     .                        NLIM,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
C**ZYLINDER ?
              ELSEIF (MERK.EQ.4) THEN
C**ZYLINDER BEGRENZT DURCH MAXIMAL 9 EBENEN
               IF (ISCN(J).EQ.0) THEN
                 CALL ZYLPLN (ZX0,ZY0,ZZ0,CX,CY,CZ,RZYL,J,NZAD,NINNE,
     .                        NIN)
C**ZYLINDER BEGRENZT DURCH MAXIMAL EINE FLAECHE ZWEITER ORDNUNG
               ELSEIF (ILIN(J).EQ.0.AND.ISCN(J).EQ.1) THEN
                 IB=1
                 CALL FL2O (ALIMS0(J,IB),XLIMS1(J,IB),YLIMS1(J,IB),
     .                      ZLIMS1(J,IB),XLIMS2(J,IB),YLIMS2(J,IB),
     .                      ZLIMS2(J,IB),XLIMS3(J,IB),YLIMS3(J,IB),
     .                      ZLIMS3(J,IB),MERK2,ZX0B,ZY0B,ZZ0B,CXB,CYB,
     .                      CZB,RZYLB,B0B,B1B,B2B,B3B,F0B,F1B,F2B,F3B,
     .                      EPS10,NMACH)
                 IF (TRCPLT) THEN
                   WRITE (6,*) 'FL2O CALLED '
                   WRITE (6,*) 'MERK2= ',MERK2
                 ENDIF
C**ZYLINDER BEGRENZT VON 2 EBENEN
                 IF (MERK2.LE.3) THEN
                   AL(1)=B0B
                   AL(2)=B1B
                   AL(3)=B2B
                   AL(4)=B3B
                   AR(1)=F0B
                   AR(2)=F1B
                   AR(3)=F2B
                   AR(4)=F3B
                   CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,4,TA,TD,LERR1)
                   CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AR,4,TB,TD,LERR2)
                   IF (LERR1.OR.LERR2) THEN
                     IF (TRCPLT)
     .               WRITE (6,*)' FEHLER IN BERANDUNG VON FLAECHE ',J
                     PLABLE(J)=.TRUE.
                     GOTO 10
                   ENDIF
                   IF (TA.LT.TB) THEN
                     T1=TA-2.*RZYL
                     T2=TB+2.*RZYL
                   ELSE
                     T1=TB-2.*RZYL
                     T2=TA+2.*RZYL
                     DO 15 II=1,4
                       TD=AL(II)
                       AL(II)=AR(II)
                       AR(II)=TD
15                   CONTINUE
                   ENDIF
                   CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                          RZYL,NZAD,NINNE,NIN,
     .                          ILCOL(J),IGFIL(J).NE.0,J,4,AL,4,AR,
     .                          0.D0,360.D0)
C**ZYLINDER BEGRENZT VON ECHT GEKRUEMMTEN FLAECHE 2TER ORDNUNG
                ELSEIF (MERK2.GE.4) THEN
                  IB=1
                  AL(1)=ALIMS0(J,IB)
                  AL(2)=XLIMS1(J,IB)
                  AL(3)=YLIMS1(J,IB)
                  AL(4)=ZLIMS1(J,IB)
                  AL(5)=XLIMS2(J,IB)
                  AL(6)=YLIMS2(J,IB)
                  AL(7)=ZLIMS2(J,IB)
                  AL(8)=XLIMS3(J,IB)
                  AL(9)=YLIMS3(J,IB)
                  AL(10)=ZLIMS3(J,IB)
                  CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,10,TA,TB,LERR1)
                  IF (LERR1) THEN
                    WRITE (6,*)' FEHLER IN DER BERANDUNG VON FLAECHE',J
                    PLABLE(J)=.TRUE.
                    GOTO 10
                  ENDIF
                  IF (TA.LT.TB) THEN
                    T1=TA-2.*RZYL
                    T2=TB+2.*RZYL
                  ELSE
                    T1=TB-2.*RZYL
                    T2=TA+2.*RZYL
                  ENDIF
                  DO 18 K=1,10
                    AR(K)=AL(K)
18                CONTINUE
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                      RZYL,NZAD,NINNE,NIN,
     .                      ILCOL(J),IGFIL(J).NE.0,
     .                      J,10,AL,10,AR,0.D0,360.D0)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
               ENDIF
C**KUGEL, ELLIPSOID
              ELSEIF (MERK == 13) THEN
                CALL ELLIPSOID (ZX0,ZY0,ZZ0,CX,CY,CZ,XLIMS1(J,1),
     .               YLIMS1(J,1),ZLIMS1(J,1),XLIMS2(J,1),YLIMS2(J,1),
     .               ZLIMS2(J,1),RLB(J),ILCOL(J),5,5,5)
              ELSE
                PLABLE(J)=.TRUE.
                GOTO 10
              ENDIF
C
C**EBENE MIT RLB.LT.0 OPTION
C
            ELSEIF (JUMLIM(J).NE.0) THEN
C
C**EBENE BEGRENZT DURCH ANDERE EBENEN
              IF (ISCN(J).EQ.0) THEN
                CALL PLANE (A0LM(J),A1LM(J),A2LM(J),A3LM(J),RLB(J),
     .                      NLIM,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                      ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                             XLIMS2,YLIMS2,ZLIMS2,
     .                             XLIMS3,YLIMS3,ZLIMS3,
     .                      ILCOL(J),IGFIL(J).NE.0,J)
              ELSEIF (ILIN(J).EQ.0) THEN
C**EBENE BEGRENZT DURCH EINEN ODER MEHRERE ZYLINDER?
                IB=0
20              IB=IB+1
                CALL FL2O (ALIMS0(J,IB),XLIMS1(J,IB),YLIMS1(J,IB),
     .                     ZLIMS1(J,IB),XLIMS2(J,IB),YLIMS2(J,IB),
     .                     ZLIMS2(J,IB),XLIMS3(J,IB),YLIMS3(J,IB),
     .                     ZLIMS3(J,IB),MERK2,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                     RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
                IF (TRCPLT) THEN
                  WRITE (6,*) 'FL2O CALLED '
                  WRITE (6,*) 'MERK2= ',MERK2
                ENDIF
                IF (MERK2.EQ.4) THEN
                  AL(1)=A0LM(J)
                  AL(2)=A1LM(J)
                  AL(3)=A2LM(J)
                  AL(4)=A3LM(J)
                  CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,4,TA,TD,LERR1)
                  IF (LERR1) THEN
                    WRITE (6,*)' FEHLER IN DER BERANDUNG VON FLAECHE',J
                    PLABLE(J)=.TRUE.
                    GOTO 10
                  ENDIF
                  T1=TA-2.*RZYL
                  T2=TA+4.*RZYL
                  NP=1
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                         RZYL,NP,NINNE,NIN,
     .                         ILCOL(J),IGFIL(J).NE.0,J,4,AL,0,AR,
     .                         0.D0,360.D0)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
                IF (IB.LT.ISCN(J)) GOTO 20
C**EBENE BEGRENZT DURCH ALLE ANDERE OPTIONEN
              ELSE
                PLABLE(J)=.TRUE.
                GOTO 10
              ENDIF
C
            ENDIF
          ENDIF
C
C END NJZ LOOP
          ENDDO
C
10      CONTINUE
100   CONTINUE
C
C
C  PLOTTE DIEJENIGEN FLAECHEN, DIE NICHT AUTOMATISCH
C  MOEGLICH WAREN, IN DER USER-SUPPLIED ROUTINE PLTUSR
      I1=0
      DO 200 J=1,NLIMI
        IF (PLABLE(J)) THEN
          IF (NLTRA) WINJ=ZZONE(ILTOR(J))
          CALL PLTUSR(PLABLE(J),J)
        ENDIF
C
        IF (PLABLE(J)) THEN
          IF (I1.EQ.0) THEN
            WRITE (6,*) 'MESSAGE FROM SUBR. PLT3D :'
            WRITE (6,*) 'SURFACE NO. J COULD NOT BE PLOTTED'
            I1=1
          ENDIF
          WRITE (6,*) 'J= ',J
        ENDIF
200   CONTINUE
C
C  PLOT SURFACES OF STANDARD MESH
C
500   IF (.NOT.(PL3S(1).OR.PL3S(2).OR.PL3S(3))) GOTO 10000
C
C  TOROIDAL GRID
C
      DO 3000 IPZ=1,IPLTS(3)
C
        NJZ=0
        IF (IPLAS(3,IPZ).LT.IPLES(3,IPZ)) THEN
          DO 3001 J=IPLAS(3,IPZ),IPLES(3,IPZ)
            NJZ=NJZ+1
            IF (J.EQ.1.OR.(NLTOR.OR.NLTRA)) THEN
              ZPLOT(NJZ)=ZSURF(J)
            ELSEIF (J.EQ.2) THEN
              ZPLOT(NJZ)=ZAA
            ENDIF
3001      CONTINUE
        ELSE
          DO 3002 J=IPLAS(3,IPZ),NTTRA
            NJZ=NJZ+1
            ZPLOT(NJZ)=ZSURF(J)
3002      CONTINUE
          DO 3003 J=2,IPLES(3,IPZ)
            NJZ=NJZ+1
            ZPLOT(NJZ)=ZSURF(J)
3003      CONTINUE
        ENDIF
C
        DO 3100 IZ=1,NJZ
          CALL GRNWPN(2)
          PHI=ZPLOT(IZ)
C
C  PHI = CONST , PLOT POLOIDAL CROSS SECTION AT TOROIDAL POSITION PHI
C
          IST=5
          IS=0
C
          DO 1000 IBR=1,IPLTS(1)
C
            IS=0
            IST=10
C
            DO 1100 IR=IPLAS(1,IBR),IPLES(1,IBR)
              IF (TRCPLT) WRITE (6,*) 'PLOT RAD. STAND. SURFACE NO. ',IR
C
C NR: POINTS TO BE PLOTTED ON RADIAL SURFACE IR
C
              IF (NLCRC.OR.NLELL.OR.NLTRI) THEN
                DM=0.
                RS=RSURF(IR)
                EP=EP1(IR)
                EL=ELL(IR)
                TR=TRI(IR)
                CALL PLGELR(RS,EP,EL,TR,DM,100,XX,YY,NR,PSURF,NP2ND)
C
                DO 1120 J=1,NR
                  Y=YY(J)
                  IF (NLTRA) THEN
                    RR=XX(J)+RMTOR
                    X=RR*COS(PHI)
                    Z=RR*SIN(PHI)
                    CALL TORLOC(WIN,RMT,X,Z)
                    WINJ=WIN
                  ELSEIF (NLTRZ) THEN
                    X=XX(J)
                    Z=PHI
                  ENDIF
C
                  CALL PL3D(X,Y,Z,XP(J),YP(J))
1120            CONTINUE
                do 1130 jj=1,nr
                  xps(jj)=xp(jj)
                  yps(jj)=yp(jj)
1130            continue
                CALL GRLN(XPS,YPS,NR)
C
C
              ELSEIF (NLPLG) THEN
C
                NR=0
                DO 1160 K=1,NPPLG
                  IAN=NPOINT(1,K)
                  IEN=NPOINT(2,K)
                  KIN=NR+1
                  DO 1165 J=IAN,IEN
                    IF (NR.GE.NPLY) THEN
                      WRITE (6,*) 'FROM PLT3D: NOT ENOUGH STORAGE   '
                      WRITE (6,*) 'INCREASE PARAMETER NPLY '
                      GOTO 1165
                    ENDIF
                    NR=NR+1
                    X=XPOL(IR,J)
                    Y=YPOL(IR,J)
                    IF (NLTRA) THEN
                      RR=X+RMTOR
                      X=RR*COS(PHI)
                      Z=RR*SIN(PHI)
                      CALL TORLOC(WIN,RMT,X,Z)
                      WINJ=WIN
                    ELSEIF (NLTRZ) THEN
                      Z=PHI
                    ENDIF
                    CALL PL3D(X,Y,Z,XP(NR),YP(NR))
1165              CONTINUE
                  do 1162 jj=1,nr
                    xps(jj)=xp(jj)
                    yps(jj)=yp(jj)
1162              continue
                  CALL GRLN (XPS,YPS,NR)
1160            CONTINUE
C
C
              ELSE
C  TO BE WRITTEN
              ENDIF
C
C  RADIAL SURFACE IR PLOTTED, NR POINTS
C  NEXT: SAVE CO-ORDINATES
C
              IF (IS+NR/IST+1.LE.NPLY) THEN
                DO 1200 J=1,NR,IST
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(J)
                  YSAVE(IS,IZ)=YP(J)
1200            CONTINUE
C  LAST POINT:
                IF (MOD(NR,IST).NE.0) THEN
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(NR)
                  YSAVE(IS,IZ)=YP(NR)
                ENDIF
              ELSE
                WRITE(6,*) ' STORAGE EXCEEDED IN PLT3D, IR= ',IR
              ENDIF
C
1100        CONTINUE
C
C  ALL RADIAL SURFACES IN BLOCK IBR DONE
C
1000      CONTINUE
C
C  ALL BLOCKS FOR RADIAL SURFACES DONE
C
C  ARE THERE TOROIDALLY (OR Z) SYMMETRIC ADDITIONAL SURFACES
C
          ISSTD=IS
          CALL GRNWPN(1)
          LZR=.FALSE.
C
          CH2X0S=CH2X0
          CH2Y0S=CH2Y0
          CH2MXS=CH2MX
          CH2MYS=CH2MY
          CH2X0=CH3X0
          CH2Y0=CH3Y0
          CH2MX=CH3MX
          CH2MY=CH3MY
          DO 1300 IBA=1,5
            IF (.NOT.PL3A(IBA)) GOTO 1300
            DO 1310 IA=1,IPLTA(IBA)
              DO 1320 J=IPLAA(IBA,IA),IPLEA(IBA,IA)
                IF (IGJUM0(J).NE.0) GOTO 1320
                IF (.NOT.LSYMET(J).OR.J.GT.NLIMI) GOTO 1320
                PLT1=PLCUT(1)
                PLT2=PLCUT(2)
                PLT3=PLCUT(3)
                PLCUT(1)=.FALSE.
                PLCUT(2)=.FALSE.
                PLCUT(3)=.TRUE.
                CALL PLTADD(J,J)
                PLCUT(1)=PLT1
                PLCUT(2)=PLT2
                PLCUT(3)=PLT3
                NA=0
C  AT PRESENT: ONLY FIRST AND LAST POINT, IF STRAIGHT LINE
C              OR 10 POINTS, IF CURVED LINE
                ISTP=MAX(1,(INSTOR-1)/10)
                IF (JUMLIM(J).GT.0) ISTP=INSTOR-1
                IF (INSTOR.LT.2) GOTO 1320
                DO 1340 ID=1,INSTOR,ISTP
                  NA=NA+1
                  X=XPL2D(ID)
                  Y=YPL2D(ID)
                  IF (NLTRA) THEN
                    RR=X+RMTOR
                    X=RR*COS(PHI)
                    Z=RR*SIN(PHI)
                    CALL TORLOC(WIN,RMT,X,Z)
                    WINJ=WIN
                  ELSEIF (NLTRZ) THEN
                    Z=PHI
                  ENDIF
                  CALL PL3D(X,Y,Z,XP(NA),YP(NA))
                  IF (NPL2D(ID).EQ.0) CALL GRJMP (sngl(XP(NA)),
     .               sngl(YP(NA)))
                  IF (NPL2D(ID).EQ.1) CALL GRDRW (sngl(XP(NA)),
     .               sngl(YP(NA)))
C
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(NA)
                  YSAVE(IS,IZ)=YP(NA)
1340            CONTINUE
1320          CONTINUE
1310        CONTINUE
1300      CONTINUE
          CH2X0=CH2X0S
          CH2Y0=CH2Y0S
          CH2MX=CH2MXS
          CH2MY=CH2MYS
C
3100    CONTINUE
C
C   LINES OF CONSTANT POLOIDAL POSITION
C
        CALL GRDSH(0.2,0.5,0.2)
        CALL GRNWPN(2)
        DO 2000 J=1,IS
          IF (J.GT.ISSTD) CALL GRNWPN(1)
          CALL GRJMP (sngl(XSAVE(J,1)),sngl(YSAVE(J,1)))
          DO 2000 IZ=2,NJZ
            CALL GRDRW(sngl(XSAVE(J,IZ)),sngl(YSAVE(J,IZ)))
2000    CONTINUE
        CALL GRDSH(1.,0.,1.)
        CALL GRNWPN(1)
C
3000  CONTINUE
C
10000 CONTINUE
C
C  BESCHRIFTUNG
C
      XH=(ABSMIN+ABSMAX)/2.
      YH=ORDMAX+2./FAKY
      CALL GRTXT (sngl(XH),sngl(YH),27,'CHECK OF GEOMETRICAL INPUT:')
      YH=YH-0.5/FAKY
      CALL GRTXT (sngl(XH),sngl(YH),72,TXTRUN)
      YH=YH-0.75/FAKY
      DO 11000 J=1,5
        IF (PL3A(J)) CALL GRTXT(sngl(XH),sngl(YH),16,TEXTLA(J))
        IF (PL3A(J)) YH=YH-0.5/FAKY
11000 CONTINUE
C
C  RETURN CO-ORDINATES FOR PLOTS OF PARTICLE TRACKS
      XR=ABSMIN-6./FAKX
      YR=ORDMAX-0./FAKY
      RETURN
      END
C
      SUBROUTINE TORLOC(W3,RM,X1,Z1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C TRANSFORM FROM TORUS- TO LOCAL SYSTEM AT PHI=W3
      SAVE WO,S,C
      DATA WO,S,C/0.,0.,1./
      IF (WO.EQ.W3) GOTO 1
      S=SIN(W3)
      C=COS(W3)
      WO=W3
1     CONTINUE
C
      XS=X1
      X1=C*XS+S*Z1
      Z1=-S*XS+C*Z1
      X1=X1-RM
      RETURN
      END
C
      SUBROUTINE LOCTOR(W3,RM,X1,Z1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C TRANSFORM FROM LOCAL SYSTEM AT PHI=W3 TO TORUS SYSTEM
      SAVE WO,S,C
      DATA WO,S,C/0.,0.,1./
      IF (WO.EQ.W3) GOTO 1
      S=SIN(-W3)
      C=COS(-W3)
      WO=W3
1     CONTINUE
C
      X1=X1+RM
      XS=X1
      X1=C*XS+S*Z1
      Z1=-S*XS+C*Z1
      RETURN
      END
C
      SUBROUTINE PL3Q(CORD,N,IO,NF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XP(101),YP(101),CORD(*)
      REAL*4 XPS(101),YPS(101)
      LOGICAL NF
      IF (N.GT.100)  WRITE (6,*) 'WARNING FROM PL3Q, N= ',N
      N3=3*N
      I=0
      DO 100 J=1,N3,3
        I=I+1
        CALL PL3D (CORD(J),CORD(J+1),CORD(J+2),XP(I),YP(I))
100   CONTINUE
      XP(N+1)=XP(1)
      YP(N+1)=YP(1)
C
      IF (IO.GE.2) CALL GRNWPN(IO)
      do 200 jj=1,n+1
        xps(jj)=xp(jj)
        yps(jj)=yp(jj)
200   continue
      CALL GRLN (XPS,YPS,N+1)
C  AUSFUELLEN DES KURVENZUGES MIT FARBE NO. IO
      IF (NF) CALL GRFILL(N+1,XPS,YPS,1,1)
      IF (IO.GE.2) CALL GRNWPN(1)
      RETURN
      END
C
C
      SUBROUTINE PLANE (A0,A1,A2,A3,RL,N1,EPS,
     .  AL,XL,YL,ZL,AL0,XL1,YL1,ZL1,XL2,YL2,ZL2,XL3,YL3,ZL3,IO,NF,NUM)
C  PLOT PLANE, SECTION INSIDE A BOX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(3,36),XYZG(3,36),ANGLE(36),CORD(108),PS(3)
      DIMENSION AL(N1,*),XL(N1,*),YL(N1,*),ZL(N1,*),
     .          AL0(N1,*),XL1(N1,*),YL1(N1,*),ZL1(N1,*),
     .                    XL2(N1,*),YL2(N1,*),ZL2(N1,*),
     .                    XL3(N1,*),YL3(N1,*),ZL3(N1,*)
      LOGICAL NF
C
      IS=0
      IF (RL.GT.0) THEN
        X1=XL1(NUM,1)
        X2=XL2(NUM,1)
        Y1=YL1(NUM,1)
        Y2=YL2(NUM,1)
        Z1=ZL1(NUM,1)
        Z2=ZL2(NUM,1)
        DX=XL2(NUM,1)-XL1(NUM,1)
        DY=YL2(NUM,1)-YL1(NUM,1)
        DZ=ZL2(NUM,1)-ZL1(NUM,1)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z1,DX,0.D0,0.D0,P(1,1),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z1,0.D0,0.D0,DZ,P(1,2),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z2,DX,0.D0,0.D0,P(1,3),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y1,Z2,0.D0,0.D0,DZ,P(1,4),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y1,Z1,0.D0,DY,0.D0,P(1,5),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y2,Z1,DX,0.D0,0.D0,P(1,6),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y2,Z1,0.D0,0.D0,DZ,P(1,7),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y2,Z2,DX,0.D0,0.D0,P(1,8),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y2,Z2,0.D0,0.D0,DZ,P(1,9),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y1,Z2,0.D0,DY,0.D0,P(1,10),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z1,0.D0,DY,0.D0,P(1,11),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z2,0.D0,DY,0.D0,P(1,12),T)
        IPOINT=12
        DO 10 I=1,IPOINT
          IF ((P(1,I)-X1.GE.-EPS.AND.X2-P(1,I).GE.-EPS).AND.
     .       (P(2,I)-Y1.GE.-EPS.AND.Y2-P(2,I).GE.-EPS).AND.
     .       (P(3,I)-Z1.GE.-EPS.AND.Z2-P(3,I).GE.-EPS)) THEN
            XYZG(1,IS+1)=P(1,I)
            XYZG(2,IS+1)=P(2,I)
            XYZG(3,IS+1)=P(3,I)
            IS=IS+1
          ENDIF
10      CONTINUE
      ELSEIF (RL.GT.-10.) THEN
        ILN=-RL
        DO 20 I=1,ILN
          IP=I+1
          DO 21 J=IP,ILN
C  SCHNITTPUNKT EBENE I, J, UND A0,A1,A2,A3
            DET=DETER(A1,XL(NUM,I),XL(NUM,J),A2,YL(NUM,I),YL(NUM,J),
     .                A3,ZL(NUM,I),ZL(NUM,J))
            IF (ABS(DET).LE.EPS) GOTO 21
            B1=-A0
            B2=-AL(NUM,I)
            B3=-AL(NUM,J)
            PS(1)=DETER(B1,B2,B3,A2,YL(NUM,I),YL(NUM,J),
     .                  A3,ZL(NUM,I),ZL(NUM,J))/DET
            PS(2)=DETER(A1,XL(NUM,I),XL(NUM,J),B1,B2,B3,
     .                  A3,ZL(NUM,I),ZL(NUM,J))/DET
            PS(3)=DETER(A1,XL(NUM,I),XL(NUM,J),A2,YL(NUM,I),YL(NUM,J),
     .                  B1,B2,B3)/DET
C  CHECKE, OB ALLE ANDEREN LINEAREN UNGLEICHUNGEN ERFUELLT SIND
            DO 40 K=1,ILN
              IF (K.EQ.I.OR.K.EQ.J) GOTO 40
              TEST=PS(1)*XL(NUM,K)+PS(2)*YL(NUM,K)+
     .             PS(3)*ZL(NUM,K)+AL(NUM,K)
              IF (TEST.GT.EPS) GOTO 21
40          CONTINUE
            IS=IS+1
            XYZG(1,IS)=PS(1)
            XYZG(2,IS)=PS(2)
            XYZG(3,IS)=PS(3)
21        CONTINUE
20      CONTINUE
      ENDIF
C
      IF (IS.LE.2) THEN
        WRITE (6,*) 'WENIGER ALS 3 SCHNITTPUNKTE ',
     .              'DER EBENE MIT DEM QUADER GEFUNDEN'
        WRITE (6,*) 'NUMMER DER FLAECHE: ',NUM
        RETURN
      ENDIF
C
      XMIT=0.
      YMIT=0.
      ZMIT=0.
      DO 200 I=1,IS
         XMIT=XMIT+XYZG(1,I)
         YMIT=YMIT+XYZG(2,I)
200      ZMIT=ZMIT+XYZG(3,I)
      XMIT=XMIT/DBLE(IS)
      YMIT=YMIT/DBLE(IS)
      ZMIT=ZMIT/DBLE(IS)
C
C     BESTIMME DIE WINKEL
      PI=4.*ATAN(1.)
      ICHECK=0
      ICOUNT=0
150   IF ((ABS(A2).LT.1.E-6.AND.ABS(A3).LT.1.E-6).OR.ICHECK.GT.1) THEN
        DO 300 I=1,IS
          ANGLE(I)=ATAN2((XYZG(3,I)-ZMIT),(XYZG(2,I)-YMIT))/PI*180.
300     CONTINUE
        ICHECK=1
      ELSEIF ((ABS(A1).LT.1.E-6.AND.ABS(A3).LT.1.E-6).OR.
     .        MOD(ICHECK,2).EQ.1) THEN
        DO 400 I=1,IS
          ANGLE(I)=ATAN2((XYZG(3,I)-ZMIT),(XYZG(1,I)-XMIT))/PI*180.
400     CONTINUE
        ICHECK=2
      ELSE
        DO 500 I=1,IS
          ANGLE(I)=ATAN2((XYZG(2,I)-YMIT),(XYZG(1,I)-XMIT))/PI*180.
500     CONTINUE
        ICHECK=3
      ENDIF
C
C     SORTIERE NACH WINKELN
      DO 600 I=1,IS-1
        ANG=ANGLE(I)
        ISORT=I
        DO 610 J=I+1,IS
          IF (ANGLE(J).LT.ANG) THEN
            ISORT=J
            ANG=ANGLE(J)
          ENDIF
610     CONTINUE
        DO 620 J=1,3
          HELP=XYZG(J,I)
          XYZG(J,I)=XYZG(J,ISORT)
620       XYZG(J,ISORT)=HELP
        HELP=ANGLE(I)
        ANGLE(I)=ANGLE(ISORT)
        ANGLE(ISORT)=HELP
600   CONTINUE
C
      ICOUNT=ICOUNT+1
      DO 700 I=1,IS-1
        DO 700 J=I+1,IS
           IF (ABS(ANGLE(J)-ANGLE(I)).LT.1.E-6.AND.
     .        (XYZG(1,I)-XYZG(1,J))**2+(XYZG(2,I)-XYZG(2,J))**2+
     .        (XYZG(3,I)-XYZG(3,J))**2.GT.1.E-6.AND.
     .         ICOUNT.LT.3) GOTO 150
700   CONTINUE
C
      II=0
      DO 1000 I=1,IS
        DO 1100 J=1,3
          II=II+1
1100      CORD(II)=XYZG(J,I)
1000  CONTINUE
      CALL PL3Q(CORD,IS,IO,NF)
C
      RETURN
      END
C
      SUBROUTINE PRLLO (P1,P2,P3,P4,P5,IO,NF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P1(3),P2(3),P3(3),P4(3),P5(3),CORD(15)
      LOGICAL NF
      DO 1 I=1,3
      CORD(I)=P1(I)
      CORD(3+I)=P2(I)
      CORD(6+I)=P4(I)
      CORD(9+I)=P5(I)
1     CORD(12+I)=P3(I)
      CALL PL3Q (CORD,5,IO,NF)
      RETURN
      END
C
      SUBROUTINE ZYLIND (X0,Y0,Z0,VX,VY,VZ,T1,T2,RAD,NK,NP,NA,IO,NF,NUM,
     .                   ILEFT,AL,IRIGHT,AR,PHIAN,PHIEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  ZYLINDERACHSE IST GERADE X+T*V, T1<T<T2, RADIUS RAD.
C  NK KREISE, NA STUETZSTELLEN AUF KREIS (POLYGON, NA-ECK)
C  NP ANZAHL DER LINIEN FUER PHI=CONST.
C  WENN ILEFT (IRIGHT) .NE. 0, DANN WIRD DER ERSTE (I=1, D.H.T=T1)
C  BZW DER LETZTE (I=NK, T=T2) KREIS DES ZYLINDERS DURCH DIE SCHNITT-
C  FLAECHE DIESES ZYLINDERS MIT DER GLEICHUNG
C  A(1)+A(2)*X+A(3)*Y+....+A(10)*Y*Z=0.
C  ERSETZT. T1 UND T2 SIND SO EINZUGEBEN, DASS DIE SCHNITTFLAECHEN
C  AUSSERHALB DIESES PARAMETERBEREICHES LIEGEN.
C  DABEI SIND NUR DIE ERSTEN ILEFT (IRIGHT) KOEFFIZIENTEN .NE.0
C  D.H. ILEFT (IRIGHT) <= 4 ENTSPRICHT DEM SCHNITT MIT EINER EBENE.
      DIMENSION P(3,101),XP(101),YP(101),AL(10),AR(10)
      REAL*4 XPS(101),YPS(101)
      LOGICAL NF
C
      NA=MIN0(NA,100)
      NP=MIN0(NP,100)
      IF (NK.EQ.1) THEN
        DT=0.
      ELSE
        DT=(T2-T1)/DBLE(NK-1)
      ENDIF
      PI=4.*ATAN(1.)
C
      AX=VX
      AY=VY
      AZ=VZ
      IF (ABS(AZ).GT.1.E-20) THEN
        BX=1.
        BY=1.
        BZ=-(AX+AY)/AZ
        PHID=ACOS(BX/SQRT(2.+BZ*BZ))
      ELSE IF (ABS(AX).GT.1.E-20) THEN
        BY=1.
        BZ=1.
        BX=-(AY+AZ)/AX
        PHID=ACOS(BY/SQRT(2.+BX*BX))
      ELSE IF (ABS(AY).GT.1.E-20) THEN
        BX=1.
        BZ=1.
        BY=-(AX+AZ)/AY
        PHID=ACOS(BZ/SQRT(2.+BY*BY))
      ELSE
        WRITE (6,*) 'FEHLER IN DER EINGABE VON VX,VY,VZ ',NUM,VX,VY,VZ
        CALL EXIT
      ENDIF
      CX=AY*BZ-AZ*BY
      CY=AZ*BX-AX*BZ
      CZ=AX*BY-AY*BX
      BA=SQRT(AX*AX+AY*AY+AZ*AZ)
      BB=SQRT(BX*BX+BY*BY+BZ*BZ)
      BC=SQRT(CX*CX+CY*CY+CZ*CZ)
      AX=AX/BA
      AY=AY/BA
      AZ=AZ/BA
      BX=BX/BB
      BY=BY/BB
      BZ=BZ/BB
      CX=CX/BC
      CY=CY/BC
      CZ=CZ/BC
C
      DANG=(PHIEN-PHIAN)/DBLE(NA)*PI/180.
      PHAN=PHIAN*PI/180.-PHID
C  BERECHNE EINEN KREIS IN RICHTUNG DER ZYLINDERACHSE, MIT
C  DEM 0-PUNKT ALS MITTELPUNKT UND RADIUS RAD
C  VON PHIAN BIS PHIEN
      DO 1 J=1,NA+1
        PHI=PHAN+(J-1)*DANG
        XK=RAD*COS(PHI)
        YK=RAD*SIN(PHI)
        P(1,J)=XK*BX+YK*CX
        P(2,J)=XK*BY+YK*CY
        P(3,J)=XK*BZ+YK*CZ
1     CONTINUE
C  SETZTE EINEN VERSCHIEBUNGSVEKTOR AUF DER ZYLINDERACHSE
C  INNERHALB DES BEREICHES T1----T2, FUER SHNITT-OPTION
      PXS=X0+(T1+T2)/2.*VX
      PYS=Y0+(T1+T2)/2.*VY
      PZS=Z0+(T1+T2)/2.*VZ
C PLOTTE DIE KREISSTUECKE, NK STUECK
      DO 2 I=1,NK
        T=T1+(I-1)*DT
        PX=X0+T*VX
        PY=Y0+T*VY
        PZ=Z0+T*VZ
        IF (I.EQ.1.AND.ILEFT.NE.0) THEN
          CALL SHNITT(P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,1,NA+1,1)
        ELSEIF (I.EQ.NK.AND.IRIGHT.NE.0) THEN
          CALL SHNITT(P,PXS,PYS,PZS,VX,VY,VZ,AR,IRIGHT,XP,YP,1,NA+1,1)
        ELSE
          DO 3 J=1,NA+1
            PXX=P(1,J)+PX
            PYY=P(2,J)+PY
            PZZ=P(3,J)+PZ
3           CALL PL3D (PXX,PYY,PZZ,XP(J),YP(J))
        ENDIF
        IF (IO.GE.2) CALL GRNWPN(IO)
        do 7 jj=1,na+1
          xps(jj)=xp(jj)
          yps(jj)=yp(jj)
7       continue
        CALL GRLN (XPS,YPS,NA+1)
C  FAERBE DIE ENDEN DES ZYLINDERS EIN
        IF ((I.EQ.1.OR.I.EQ.NK).AND.NF) CALL GRFILL(NA+1,XPS,YPS,1,1)
        IF (IO.GE.2) CALL GRNWPN(1)
2     CONTINUE
C
C  PLOTTE PHI=CONST LINIEN, INSGESAMT NP STUECK
C
C  SETZE NEUEN KREIS UM 0-PUNKT, MIT NP STUETZSTELLEN
      DANG=(PHIEN-PHIAN)/DBLE(NP-1)*PI/180.
      PHAN=PHIAN*PI/180.-PHID
      DO 6 J=1,NP
        PHI=PHAN+(J-1)*DANG
        XK=RAD*COS(PHI)
        YK=RAD*SIN(PHI)
        P(1,J)=XK*BX+YK*CX
        P(2,J)=XK*BY+YK*CY
        P(3,J)=XK*BZ+YK*CZ
6     CONTINUE
C
      IF (NK.EQ.1) RETURN
      IA=1
      IE=NK
      IF (ILEFT.NE.0) IA=2
      IF (IRIGHT.NE.0) IE=NK-1
      DO 5 J=1,NP
        IF (ILEFT.NE.0) THEN
          CALL SHNITT (P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,J,J,1)
        ENDIF
        DO 4 I=IA,IE
          T=T1+(I-1)*DT
          PX=X0+T*VX
          PY=Y0+T*VY
          PZ=Z0+T*VZ
          CALL PL3D (P(1,J)+PX,P(2,J)+PY,P(3,J)+PZ,XP(I),YP(I))
4       CONTINUE
        IF (IRIGHT.NE.0) THEN
          CALL SHNITT (P,PXS,PYS,PZS,VX,VY,VZ,AR,IRIGHT,XP,YP,J,J,NK)
        ENDIF
        do 9 jj=1,nk
          xps(jj)=xp(jj)
          yps(jj)=yp(jj)
9       continue
        CALL GRLN (XPS,YPS,NK)
5     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CONE (X0,Y0,Z0,VX,VY,VZ,T1,T2,ALF,NK,NP,NA,IO,NF,NUM,
     .                 ILEFT,AL,IRIGHT,AR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  KEGELACHSE IST GERADE X+T*V, T1<T<T2, OEFFNUNGSWINKEL ALF
C  NK KREISE, NA STUETZSTELLEN AUF KREIS (POLYGON, NA-ECK)
C  NP ANZAHL DER LINIEN FUER PHI=CONST.
C  WENN ILEFT (IRIGHT) .NE. 0, DANN WIRD DER ERSTE (I=1, D.H.T=T1)
C  BZW DER LETZTE (I=NK, T=T2) KREIS DES KEGELS DURCH DIE SCHNITT-
C  FLAECHE DIESES KEGELS MIT DER GLEICHUNG
C  A(1)+A(2)*X+A(3)*Y+....+A(10)*Y*Z=0.
C  ERSETZT. T1 UND T2 SIND SO EINZUGEBEN, DASS DIE SCHNITTFLAECHEN
C  AUSSERHALB DIESES PARAMETERBEREICHES LIEGEN.
C  DABEI SIND NUR DIE ERSTEN ILEFT (IRIGHT) KOEFFIZIENTEN .NE.0
C  D.H. ILEFT (IRIGHT) <= 4 ENTSPRICHT DEM SCHNITT MIT EINER EBENE.
      DIMENSION XP(101),YP(101),AL(10),AR(10)
      REAL*4 XPS(101),YPS(101)
      LOGICAL NF
C
      NA=MIN0(NA,100)
      NP=MIN0(NP,100)
      IF (NK.EQ.1) THEN
        DT=0.
      ELSE
        DT=(T2-T1)/DBLE(NK-1)
      ENDIF
      PI=4.*ATAN(1.)
C
      AX=VX
      AY=VY
      AZ=VZ
      IF (ABS(AZ).GT.1.E-20) THEN
        BX=1.
        BY=1.
        BZ=-(AX+AY)/AZ
      ELSE IF (ABS(AX).GT.1.E-20) THEN
        BY=1.
        BZ=1.
        BX=-(AY+AZ)/AX
      ELSE IF (ABS(AY).GT.1.E-20) THEN
        BX=1.
        BZ=1.
        BY=-(AX+AZ)/AY
      ELSE
        WRITE (6,*) 'FEHLER IN DER EINGABE VON VX,VY,VZ ',NUM,VX,VY,VZ
        CALL EXIT
      ENDIF
      CX=AY*BZ-AZ*BY
      CY=AZ*BX-AX*BZ
      CZ=AX*BY-AY*BX
      BA=SQRT(AX*AX+AY*AY+AZ*AZ)
      BB=SQRT(BX*BX+BY*BY+BZ*BZ)
      BC=SQRT(CX*CX+CY*CY+CZ*CZ)
      AX=AX/BA
      AY=AY/BA
      AZ=AZ/BA
      BX=BX/BB
      BY=BY/BB
      BZ=BZ/BB
      CX=CX/BC
      CY=CY/BC
      CZ=CZ/BC
C
      DANG=2.*PI/DBLE(NA)
C  BERECHNE EINEN KREIS IN RICHTUNG DER ZYLINDERACHSE, MIT
C  DEM 0-PUNKT ALS MITTELPUNKT UND RADIUS RAD
C  SETZTE EINEN VERSCHIEBUNGSVEKTOR AUF DER ZYLINDERACHSE
C  INNERHALB DES BEREICHES T1----T2, FUER SHNITT-OPTION
      TH=(T1+T2)/2.
C PLOTTE DIE KREISE
      DO 2 I=1,NK
        T=T1+(I-1)*DT
        IF (I.EQ.1.AND.ILEFT.NE.0) THEN
          CALL SCCONE(X0,Y0,Z0,-VX,-VY,-VZ,ALF,TH,T1,BX,BY,BZ,CX,CY,CZ,
     .                DANG,AL,ILEFT,XP,YP,1,NA+1,1)
        ELSEIF (I.EQ.NK.AND.IRIGHT.NE.0) THEN
          CALL SCCONE(X0,Y0,Z0,VX,VY,VZ,ALF,TH,T2,BX,BY,BZ,CX,CY,CZ,
     .                DANG,AR,IRIGHT,XP,YP,1,NA+1,1)
        ELSE
        PX=X0+T*VX
        PY=Y0+T*VY
        PZ=Z0+T*VZ
        RAD=T*TAN(ALF)
        DO 3 J=1,NA+1
          PHI=(J-1)*DANG
          XK=RAD*COS(PHI)
          YK=RAD*SIN(PHI)
          PXX=XK*BX+YK*CX+PX
          PYY=XK*BY+YK*CY+PY
          PZZ=XK*BZ+YK*CZ+PZ
3         CALL PL3D (PXX,PYY,PZZ,XP(J),YP(J))
        ENDIF
        IF (IO.GE.2) CALL GRNWPN(IO)
        do 7 jj=1,na+1
          xps(jj)=xp(jj)
          yps(jj)=yp(jj)
7       continue
        CALL GRLN (XPS,YPS,NA+1)
C  FAERBE DIE ENDEN DES CONES EIN
        IF ((I.EQ.1.OR.I.EQ.NK).AND.NF) CALL GRFILL(NA+1,XPS,YPS,1,1)
        IF (IO.GE.2) CALL GRNWPN(1)
2     CONTINUE
C
C  SETZE NEUEN KREIS UM 0-PUNKT, MIT NP STUETZSTELLEN
      DANG=2.*PI/DBLE(NP)
C
C  PLOTTE PHI=CONST LINIEN, INSGESAMT NP STUECK
      IA=1
      IE=NK
      IF (ILEFT.NE.0) IA=2
      IF (IRIGHT.NE.0) IE=NK-1
      DO 5 J=1,NP
        IF (ILEFT.NE.0) THEN
          CALL SCCONE(X0,Y0,Z0,-VX,-VY,-VZ,ALF,TH,T1,BX,BY,BZ,CX,CY,CZ,
     .                DANG,AL,ILEFT,XP,YP,J,J,1)
        ENDIF
        DO 4 I=IA,IE
          T=T1+(I-1)*DT
          PX=X0+T*VX
          PY=Y0+T*VY
          PZ=Z0+T*VZ
          RAD=T*TAN(ALF)
          PHI=(J-1)*DANG
          XK=RAD*COS(PHI)
          YK=RAD*SIN(PHI)
          PXX=XK*BX+YK*CX
          PYY=XK*BY+YK*CY
          PZZ=XK*BZ+YK*CZ
          CALL PL3D (PXX+PX,PYY+PY,PZZ+PZ,XP(I),YP(I))
4       CONTINUE
        IF (IRIGHT.NE.0) THEN
          CALL SCCONE(X0,Y0,Z0,VX,VY,VZ,ALF,TH,T2,BX,BY,BZ,CX,CY,CZ,
     .                DANG,AR,IRIGHT,XP,YP,J,J,NK)
        ENDIF
        do 9 jj=1,nk
          xps(jj)=xp(jj)
          yps(jj)=yp(jj)
9       continue
        CALL GRLN (XPS,YPS,NK)
5     CONTINUE
      RETURN
      END
C
      SUBROUTINE SHNITT(P,PX,PY,PZ,VX,VY,VZ,A,I,XP,YP,JA,JE,IXS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CALLED FROM ZYLIND
      DIMENSION P(3,*),A(*),XP(*),YP(*)
      LOGICAL LERR
      DATA EPS12 /1.E-12/
      LERR=.FALSE.
      IX=IXS-1
C     WRITE (6,*) ' SHNITT  I = ',I
      IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
C       WRITE (6,*) ' VX,VY,VZ ',VX,VY,VZ
C       WRITE (6,*) ' A ',A(1),A(2),A(3),A(4)
        XN=VX*A(2)+VY*A(3)+VZ*A(4)
C       WRITE (6,*) ' XN ',XN
        DO 100 J=JA,JE
          XX=P(1,J)+PX
          YY=P(2,J)+PY
          ZZ=P(3,J)+PZ
          IF (LERR) THEN
            ALAMDA=0.
            GOTO 101
          ENDIF
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
          IF (ABS(XN).LT.EPS12) THEN
            WRITE (6,*) 'ERROR IN SUBR. SHNITT. SET ALAMDA=0.'
            WRITE (6,*) 'NO INTERSECTION FOUND WITH PLANE'
            ALAMDA=0.
            LERR=.TRUE.
          ELSE
            ALAMDA=(-A(1)-(A(2)*XX+A(3)*YY+A(4)*ZZ))/XN
            IF (ALAMDA.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SHNITT. SET ALAMDA=0.'
              WRITE (6,*) 'NO INTERSECTION IN POSITIV DIRECTION'
              WRITE (6,*) 'WITH PLANE '
              ALAMDA=0.
              LERR=.TRUE.
            ENDIF
          ENDIF
C
101       XX=XX+ALAMDA*VX
          YY=YY+ALAMDA*VY
          ZZ=ZZ+ALAMDA*VZ
          IX=IX+1
          CALL PL3D(XX,YY,ZZ,XP(IX),YP(IX))
100     CONTINUE
      ELSEIF (I.GT.4) THEN
C  SCHNITTKURVE MIT VOLLER GLEICHUNG 2TER ORDNUNG
        AA=(A(5)*VX+A(8)*VY+A(9)*VZ)*VX+(A(6)*VY+A(10)*VZ)*VY+
     .      A(7)*VZ*VZ
        DO 200 J=JA,JE
          XX=P(1,J)+PX
          YY=P(2,J)+PY
          ZZ=P(3,J)+PZ
          IF (LERR) THEN
            ALAMDA=0.
            GOTO 201
          ENDIF
          BB=(A(2)+2.*A(5)*XX+A(8)*YY+A(9)*ZZ)*VX+
     .       (A(3)+2.*A(6)*YY+A(8)*XX+A(10)*ZZ)*VY+
     .       (A(4)+2.*A(7)*ZZ+A(9)*XX+A(10)*YY)*VZ
          CC=A(1)+(A(2)+A(5)*XX+A(8)*YY+A(9)*ZZ)*XX+
     .       (A(3)+A(6)*YY+A(10)*ZZ)*YY+(A(4)+A(7)*ZZ)*ZZ
C
          IF (ABS(AA).LT.EPS12.AND.ABS(BB).LT.EPS12) THEN
            LERR=.TRUE.
            ALAMDA=0.
            WRITE (6,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'AA=BB=0.'
            GOTO 201
          ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
            ALAM1=-CC/BB
            ALAM2=-CC/BB
          ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
            DET=-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
              WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (6,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 201
            ENDIF
            ALAM1=SQRT(DET)
            ALAM2=-SQRT(DET)
          ELSE IF (ABS(CC).LT.EPS12) THEN
            ALAM1=0.
            ALAM2=-BB/AA
          ELSE
            DET=BB*BB/(4*AA*AA)-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
              WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (6,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 201
            ENDIF
            ROOT=SQRT(DET)
            ALAM1=-BB/(2.*AA)+ROOT
            ALAM2=-BB/(2.*AA)-ROOT
          ENDIF
          IF (LERR) GOTO 201
C  DECIDE, WHICH ONE OF THE 2 SOLUTIONS TO TAKE
          IF (ALAM1*ALAM2.LT.0.) THEN
            ALAMDA=MAX(ALAM1,ALAM2)
          ELSE IF (ABS(ALAM1).GT.ABS(ALAM2)) THEN
            ALAMDA=ALAM2
          ELSE
            ALAMDA=ALAM1
          ENDIF
          IF (ALAMDA.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0.'
            WRITE (6,*) 'INTERSECTION IN WRONG DIRECTION'
            ALAMDA=0.
            LERR=.TRUE.
            GOTO 201
          ENDIF
201       XX=XX+ALAMDA*VX
          YY=YY+ALAMDA*VY
          ZZ=ZZ+ALAMDA*VZ
          IX=IX+1
          CALL PL3D(XX,YY,ZZ,XP(IX),YP(IX))
200     CONTINUE
      ENDIF
      RETURN
      END
C
      SUBROUTINE SCCONE(X0,Y0,Z0,VX,VY,VZ,ALF,T1,T2,BX,BY,BZ,CX,CY,CZ,
     .                  DANG,A,I,XP,YP,JA,JE,IXS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CALLED FROM CONE
      DIMENSION A(*),XP(*),YP(*)
      LOGICAL LERR
      DATA EPS12 /1.E-12/
      LERR=.FALSE.
      IX=IXS-1
      RAD1=T1*TAN(ALF)
      RAD2=T2*TAN(ALF)
      PX1=X0+T1*VX
      PY1=Y0+T1*VY
      PZ1=Z0+T1*VZ
      PX2=X0+T2*VX
      PY2=Y0+T2*VY
      PZ2=Z0+T2*VZ
      DO 100 J=JA,JE
        PHI=(J-1)*DANG
        XK1=RAD1*COS(PHI)
        YK1=RAD1*SIN(PHI)
        PXX1=XK1*BX+YK1*CX
        PYY1=XK1*BY+YK1*CY
        PZZ1=XK1*BZ+YK1*CZ
        XK2=RAD2*COS(PHI)
        YK2=RAD2*SIN(PHI)
        PXX2=XK2*BX+YK2*CX
        PYY2=XK2*BY+YK2*CY
        PZZ2=XK2*BZ+YK2*CZ
        PX21=PXX2-PXX1
        PY21=PYY2-PYY1
        PZ21=PZZ2-PZZ1
        IF (LERR) THEN
          ALAMDA=0.
          GOTO 101
        ENDIF
        IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
          XN=PX21*A(2)+PY21*A(3)+PZ21*A(4)
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
          IF (ABS(XN).LT.EPS12) THEN
            WRITE (6,*) 'ERROR IN SUBR. SCCONE. SET ALAMDA=0.'
            WRITE (6,*) 'NO INTERSECTION FOUND WITH PLANE'
            ALAMDA=0.
            LERR=.TRUE.
          ELSE
            ALAMDA=(-A(1)-(A(2)*PXX1+A(3)*PYY1+A(4)*PZZ1))/XN
            IF (ALAMDA.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SCCONE. SET ALAMDA=0.'
              WRITE (6,*) 'NO INTERSECTION IN POSITIV DIRECTION'
              WRITE (6,*) 'WITH PLANE '
              ALAMDA=0.
              LERR=.TRUE.
            ENDIF
          ENDIF
C
        ELSEIF (I.GT.4) THEN
C  SCHNITTKURVE MIT VOLLER GLEICHUNG 2TER ORDNUNG
          AA=(A(5)*PX21+A(8)*PY21+A(9)*PZ21)*PX21+
     .       (A(6)*PY21+A(10)*PZ21)*PY21+A(7)*PZ21*PZ21
          BB=(A(2)+2.*A(5)*PXX1+A(8)*PYY1+A(9)*PZZ1)*PX21+
     .       (A(3)+2.*A(6)*PYY1+A(8)*PXX1+A(10)*PZZ1)*PY21+
     .       (A(4)+2.*A(7)*PZZ1+A(9)*PXX1+A(10)*PYY1)*PZ21
          CC=A(1)+(A(2)+A(5)*PXX1+A(8)*PYY1+A(9)*PZZ1)*PXX1+
     .       (A(3)+A(6)*PYY1+A(10)*PZZ1)*PYY1+(A(4)+A(7)*PZZ1)*PZZ1
C
          IF (ABS(AA).LT.EPS12.AND.ABS(BB).LT.EPS12) THEN
            LERR=.TRUE.
            ALAMDA=0.
            WRITE (6,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'AA=BB=0.'
            GOTO 101
          ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
            ALAM1=-CC/BB
            ALAM2=-CC/BB
          ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
            DET=-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
              WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (6,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 101
            ENDIF
            ALAM1=SQRT(DET)
            ALAM2=-SQRT(DET)
          ELSE IF (ABS(CC).LT.EPS12) THEN
            ALAM1=0.
            ALAM2=-BB/AA
          ELSE
            DET=BB*BB/(4*AA*AA)-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
              WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (6,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 101
            ENDIF
            ROOT=SQRT(DET)
            ALAM1=-BB/(2.*AA)+ROOT
            ALAM2=-BB/(2.*AA)-ROOT
          ENDIF
          IF (LERR) GOTO 101
C  DECIDE, WHICH ONE OF THE 2 SOLUTIONS TO TAKE
          IF (ALAM1*ALAM2.LT.0.) THEN
            ALAMDA=MAX(ALAM1,ALAM2)
          ELSE IF (ABS(ALAM1).GT.ABS(ALAM2)) THEN
            ALAMDA=ALAM2
          ELSE
            ALAMDA=ALAM1
          ENDIF
          IF (ALAMDA.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0.'
            WRITE (6,*) 'INTERSECTION IN WRONG DIRECTION'
            ALAMDA=0.
            LERR=.TRUE.
            GOTO 101
          ENDIF
        ENDIF
C
101     XX=PXX1+ALAMDA*PX21
        YY=PYY1+ALAMDA*PY21
        ZZ=PZZ1+ALAMDA*PZ21
        IX=IX+1
        CALL PL3D(XX,YY,ZZ,XP(IX),YP(IX))
100   CONTINUE
      RETURN
      END
C
      SUBROUTINE SECQUA(XX,YY,ZZ,VX,VY,VZ,A,I,ALAM1,ALAM2,LERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C MORE GENERAL THAN SHNITT OR CSCONE, CALLED FROM PL3D ITSELF
      DIMENSION A(*)
      LOGICAL LERR
      DATA EPS12 /1.E-12/
      LERR=.FALSE.
      IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
        XN=VX*A(2)+VY*A(3)+VZ*A(4)
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
        IF (ABS(XN).LT.EPS12) THEN
          WRITE (6,*) 'ERROR IN SUBR. SECQUA. SET ALAM1=ALAM2=0.'
          WRITE (6,*) 'NO INTERSECTION FOUND WITH PLANE'
          ALAM1=0.
          ALAM2=0.
          LERR=.TRUE.
        ELSE
          ALAM1=(-A(1)-(A(2)*XX+A(3)*YY+A(4)*ZZ))/XN
          ALAM2=ALAM1
        ENDIF
C
      ELSEIF (I.GT.4) THEN
C  SCHNITTKURVE MIT VOLLER GLEICHUNG 2TER ORDNUNG
        AA=(A(5)*VX+A(8)*VY+A(9)*VZ)*VX+(A(6)*VY+A(10)*VZ)*VY+
     .      A(7)*VZ*VZ
        BB=(A(2)+2.*A(5)*XX+A(8)*YY+A(9)*ZZ)*VX+
     .     (A(3)+2.*A(6)*YY+A(8)*XX+A(10)*ZZ)*VY+
     .     (A(4)+2.*A(7)*ZZ+A(9)*XX+A(10)*YY)*VZ
        CC=A(1)+(A(2)+A(5)*XX+A(8)*YY+A(9)*ZZ)*XX+
     .     (A(3)+A(6)*YY+A(10)*ZZ)*YY+(A(4)+A(7)*ZZ)*ZZ
C
        IF (ABS(AA).LT.EPS12.AND.ABS(BB).LT.EPS12) THEN
          LERR=.TRUE.
          ALAM1=0.
          ALAM2=0.
          WRITE (6,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
          WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
          WRITE (6,*) 'AA=BB=0.'
        ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
          ALAM1=-CC/BB
          ALAM2=-CC/BB
        ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
          DET=-CC/AA
          IF (DET.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'DETERMINANT DET= ',DET
            ALAM1=0.
            ALAM2=0.
            LERR=.TRUE.
          ELSE
            ALAM1=SQRT(DET)
            ALAM2=-SQRT(DET)
          ENDIF
        ELSE IF (ABS(CC).LT.EPS12) THEN
          ALAM1=0.
          ALAM2=-BB/AA
        ELSE
          DET=BB*BB/(4*AA*AA)-CC/AA
          IF (DET.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'DETERMINANT DET= ',DET
            ALAM1=0.
            ALAM2=0.
            LERR=.TRUE.
          ELSE
            ROOT=SQRT(DET)
            ALAM1=-BB/(2.*AA)+ROOT
            ALAM2=-BB/(2.*AA)-ROOT
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE CTQUA (A0LM,A1LM,A2LM,A3LM,A4LM,A5LM,A6LM,A7LM,A8LM,
     .                  A9LM,XLIMS1,XLIMS2,YLIMS1,YLIMS2,ZLIMS1,ZLIMS2,
     .                  RLB,RAD,CX,CY,CZ,PHIAN,PHIEN,IPART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PHIAN(*),PHIEN(*),PHIANG(10)
      LOGICAL LX0,LY0,LZ0,LCTX1,LCTX2,LCTY1,LCTY2
C
      EPS10=1.E-10
      PI=4.*ATAN(1.)
      PI180=PI/180.
C
      LX0=A1LM**2+A4LM**2+A7LM**2+A8LM**2.LT.EPS10
      LY0=A2LM**2+A5LM**2+A7LM**2+A9LM**2.LT.EPS10
      LZ0=A3LM**2+A6LM**2+A8LM**2+A9LM**2.LT.EPS10
C
      IF (.NOT.(LX0.OR.LY0.OR.LZ0)) THEN
C  SCHIEFER ZYLINDER
        IPART=1
        PHIAN(1)=0.
        PHIEN(1)=360.
        RETURN
      ENDIF
C
      IF ((LX0.AND.LY0).OR.(LX0.AND.LZ0).OR.(LY0.AND.LZ0)) THEN
        WRITE (6,*) ' ERROR IN CTQUA '
        WRITE (6,*) ' EQUATION DOES NOT DESCRIBE A CYLINDER '
        WRITE (6,*) ' IT DESCRIBES TWO PARALLEL PLANES '
        WRITE (6,*) ' NO PLOT IS DONE '
        IPART=0
        RETURN
      ENDIF
C
      IF (LX0) THEN
C  CYLINDER PARALLEL TO X
        A0=A0LM
        A1=A2LM
        A2=A3LM
        A3=A5LM
        A4=A6LM
        XL1=YLIMS1
        XL2=YLIMS2
        YL1=ZLIMS1
        YL2=ZLIMS2
        XM=CY
        YM=CZ
      ELSEIF (LY0) THEN
C  CYLINDER PARALLEL TO Y
        A0=A0LM
        A1=A1LM
        A2=A3LM
        A3=A4LM
        A4=A6LM
        XL1=XLIMS1
        XL2=XLIMS2
        YL1=ZLIMS1
        YL2=ZLIMS2
        XM=CX
        YM=CZ
      ELSE
C  CYLINDER PARALLEL TO Z
        A0=A0LM
        A1=A1LM
        A2=A2LM
        A3=A4LM
        A4=A5LM
        XL1=XLIMS1
        XL2=XLIMS2
        YL1=YLIMS1
        YL2=YLIMS2
        XM=CX
        YM=CY
      ENDIF
C
      LCTX1=XM-RAD.LE.XL1
      LCTX2=XM+RAD.GE.XL2
      LCTY1=YM-RAD.LE.YL1
      LCTY2=YM+RAD.GE.YL2
C
      IF (.NOT.(LCTX1.OR.LCTX2.OR.LCTY1.OR.LCTY2)) THEN
C  ZYLINDER LIEGT KOMPLETT IM QUADER
        IPART=1
        PHIAN(1)=0.
        PHIEN(1)=360.
        RETURN
      ENDIF
C
      IPHI=0
      IF (LCTX1) CALL CTCIRC (A0,A1,A2,A3,A4,XL1,YL1,YL2,XM,YM,
     .                        PHIANG,IPHI,0)
      IF (LCTX2) CALL CTCIRC (A0,A1,A2,A3,A4,XL2,YL1,YL2,XM,YM,
     .                        PHIANG,IPHI,0)
      IF (LCTY1) CALL CTCIRC (A0,A2,A1,A4,A3,YL1,XL1,XL2,YM,XM,
     .                        PHIANG,IPHI,1)
      IF (LCTY2) CALL CTCIRC (A0,A2,A1,A4,A3,YL2,XL1,XL2,YM,XM,
     .                        PHIANG,IPHI,1)
C
C  SORTIERE WINKEL
10    ISORT=0
      DO 15 I=1,IPHI-1
        IF (PHIANG(I+1).LT.PHIANG(I)) THEN
          PHIH=PHIANG(I)
          PHIANG(I)=PHIANG(I+1)
          PHIANG(I+1)=PHIH
          ISORT=ISORT+1
        ENDIF
15    CONTINUE
      IF (ISORT.GT.0) GOTO 10
C
      IPHI=IPHI+1
      PHIANG(IPHI)=PHIANG(1)
C
      IPART=0
      DO 20 I=1,IPHI-1
        IF (ABS(PHIANG(I+1)-PHIANG(I)).LT.EPS10) GOTO 20
        IF (PHIANG(I+1).LT.PHIANG(I)) PHIANG(I+1)=PHIANG(I+1)+360.
        PHI=0.5*(PHIANG(I)+PHIANG(I+1))*PI180
        XPHI=XM+RAD*COS(PHI)
        YPHI=YM+RAD*SIN(PHI)
        IF (RLB.LT.1.5) THEN
          IF (XPHI.GE.XL1.AND.XPHI.LE.XL2.AND.
     .        YPHI.GE.YL1.AND.YPHI.LE.YL2) THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ELSE
          IF (XPHI.LE.XL1.OR.XPHI.GE.XL2.OR.YPHI.LE.YL1.OR.YPHI.GE.YL2)
     .       THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ENDIF
20    CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE CTCIRC (A0,A1,A2,A3,A4,XP,YL1,YL2,XM,YM,
     .                   PHIANG,IPHI,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PHIANG(*)
C
C   SOLVE  A0+A1*X+A2*Y+A3*X**2+A4*Y**2=0 FOR X=XP
C
      PI=4.*ATAN(1.)
      RD=0.25*A2**2-(A0+A1*XP+A3*XP*XP)*A4
      IF (RD.LT.0.) THEN
        WRITE (6,*) ' WURZEL NEGATIV IN CTCIRC '
        RETURN
      ENDIF
      YP1=(-0.5*A2+SQRT(RD))/A4
      YP2=(-0.5*A2-SQRT(RD))/A4
C
      IF (ISW.EQ.0) THEN
        PHI1=ATAN2((YP1-YM),(XP-XM))/PI*180.
        PHI2=ATAN2((YP2-YM),(XP-XM))/PI*180.
      ELSE
        PHI1=ATAN2((XP-XM),(YP1-YM))/PI*180.
        PHI2=ATAN2((XP-XM),(YP2-YM))/PI*180.
      ENDIF
C
      PHI1=MOD(PHI1+360.D0,360.D0)
      PHI2=MOD(PHI2+360.D0,360.D0)
C
      IF (YP1.GE.YL1.AND.YP1.LE.YL2) THEN
        IPHI=IPHI+1
        PHIANG(IPHI)=PHI1
      ENDIF
      IF (YP2.GE.YL1.AND.YP2.LE.YL2) THEN
        IPHI=IPHI+1
        PHIANG(IPHI)=PHI2
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PL3D(PX,PY,PZ,PP1,PP2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CPL3D'
      INCLUDE 'CPLOT'
      SAVE
      DATA IFIRST,IPBOX/0,0/

      FXI(X,Y) = (F10+F11*X)+F12*Y
      FZETA(X,Y,Z) = ((F00+F01*X)+F02*Y)+F03*Z

      X=PX
      Y=PY
      Z=PZ
      XMIN = CH3X0-CH3MX
      YMIN = CH3Y0-CH3MY
      ZMIN = CH3Z0-CH3MZ
      XMAX = CH3X0+CH3MX
      YMAX = CH3Y0+CH3MY
      ZMAX = CH3Z0+CH3MZ
      BREITE = 24.
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        FA = ATAN(1.)/45.
        W3=SQRT(3.)
        XO=.5*(XMIN+XMAX)
        DX=XMAX-XMIN
        YO=.5*(YMIN+YMAX)
        DY=YMAX-YMIN
        ZO=.5*(ZMIN+ZMAX)
        DZ=ZMAX-ZMIN
        TH=FA*ANGLE2
        CH=FA*ANGLE1
        SIT=SIN(TH)
        SIP=SIN(CH)
        COT=COS(TH)
        COP=COS(CH)
        BREI=BREITE/W3
        F00=(.5+((XO/DX*SIP-YO/DY*COP)*SIT-ZO/DZ*COT)/W3)*BREITE
        F01=-SIP*SIT/DX*BREI
        F02=COP*SIT/DY*BREI
        F03=COT/DZ*BREI
        F10=(.5-(YO/DY*SIP+XO/DX*COP)/W3)*BREITE
        F11=COP/DX*BREI
        F12=SIP/DY*BREI
      ENDIF
      IF (PLBOX.AND.IPBOX.EQ.0) THEN
        IPBOX=1
        PP1 = FXI(CH3X0,CH3Z0)
        PP2 = FZETA(CH3X0,CH3Z0,CH3Y0)
        CALL GRMRKS(0.4)
        CALL GRNWPN(2)
        CALL GRJMPS(sngl(PP1),sngl(PP2),203)
        CALL GRMRKS(0.2)
        PPX11= FXI(XMIN,ZMIN)
        PPY11= FZETA(XMIN,ZMIN,YMIN)
        CALL GRJMP(sngl(PPX11),sngl(PPY11))
        PPX12 = FXI(XMAX,ZMIN)
        PPY12 = FZETA(XMAX,ZMIN,YMIN)
        CALL GRDRW(sngl(PPX12),sngl(PPY12))
        PPX13 = FXI(XMAX,ZMIN)
        PPY13 = FZETA(XMAX,ZMIN,YMAX)
        CALL GRDRW(sngl(PPX13),sngl(PPY13))
        PPX14 = FXI(XMIN,ZMIN)
        PPY14 = FZETA(XMIN,ZMIN,YMAX)
        CALL GRDRW(sngl(PPX14),sngl(PPY14))
        CALL GRDRW(sngl(PPX11),sngl(PPY11))
        PPX21= FXI(XMIN,ZMAX)
        PPY21= FZETA(XMIN,ZMAX,YMIN)
        CALL GRDRW(sngl(PPX21),sngl(PPY21))
        PPX22 = FXI(XMAX,ZMAX)
        PPY22 = FZETA(XMAX,ZMAX,YMIN)
        CALL GRDRW(sngl(PPX22),sngl(PPY22))
        PPX23 = FXI(XMAX,ZMAX)
        PPY23 = FZETA(XMAX,ZMAX,YMAX)
        CALL GRDRW(sngl(PPX23),sngl(PPY23))
        PPX24 = FXI(XMIN,ZMAX)
        PPY24 = FZETA(XMIN,ZMAX,YMAX)
        CALL GRDRW(sngl(PPX24),sngl(PPY24))
        CALL GRDRW(sngl(PPX21),sngl(PPY21))
        CALL GRJMP(sngl(PPX12),sngl(PPY12))
        CALL GRDRW(sngl(PPX22),sngl(PPY22))
        CALL GRJMP(sngl(PPX13),sngl(PPY13))
        CALL GRDRW(sngl(PPX23),sngl(PPY23))
        CALL GRJMP(sngl(PPX14),sngl(PPY14))
        CALL GRDRW(sngl(PPX24),sngl(PPY24))
        CALL GRNWPN(1)
      ENDIF
C
      IF (RMT.EQ.0.OR.WIN.EQ.WINJ) GOTO 1
C
C  TOROIDAL GEOMETRY, PLOT IN LOCAL CO-ORDINATE SYSTEM AT WIN, BUT
C  X,Y,Z ARE GIVEN IN LOCAL CO-ORDIANTE SYSTEM AT WINJ
      CALL LOCTOR(WINJ,RMT,X,Z)
      CALL TORLOC(WIN,RMT,X,Z)
C
1     CONTINUE
C
      PP1 = FXI(X,Z)
      PP2 = FZETA(X,Z,Y)
      END
C
C
      SUBROUTINE PLTEIR (ISTRA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA
C  PLOT PLASMA TALLIES ONLY ONCE.
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'COMSOU'
      INCLUDE 'CTRCEI'
      INCLUDE 'COUTAU'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CESTIM'
      INCLUDE 'CPLOT'
      INCLUDE 'CGRPTL'
      INCLUDE 'CCONA'
      INCLUDE 'CSPEI'
      INCLUDE 'CPLMSK'
      INCLUDE 'CTEXT'
      INCLUDE 'CSDVI'
      INCLUDE 'CSDVI_BGK'
      INCLUDE 'CSDVI_COP'
      INCLUDE 'CLOGAU'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCOUPL'
C
      DIMENSION ESTIM(NESTIM),OUTAU(NOUTAU),SDVI(NSDVI),PLPRM(NPLPRM),
     .          SDVI_BGK(NSBGK),SDVI_COP(NSCOP)
      DIMENSION VECTOR(NRAD,NPLT),VECSAV(NRAD,NPLT),
     .          VSDVI(NRAD,NPLT),
     .          YMN2(NPLT),YMX2(NPLT),YMNLG2(NPLT),YMXLG2(NPLT)
      DIMENSION IR1(NPLT),IR2(NPLT),IRS(NPLT)
      LOGICAL LPLOT2(NPLT),LSDVI(NPLT),LPLTT2,LINLOG
      CHARACTER*24 TXUNIT(NPLT),TXSPEC(NPLT)
      CHARACTER*24 TXUNT1,TXSPC1
      CHARACTER*72 TXTALL(NPLT)
      CHARACTER*72 TXTLL1
      CHARACTER*72 HEAD,HEAD0,HEAD1,HEAD2,HEAD3,HEAD4,HEAD5,HEAD6,HEAD7,
     .             HEAD8,TXHEAD
C
      EQUIVALENCE
     .  (RWK(NID3P),VECTOR(1,1)),
     .  (RWK(NID3P+NRAD*NPLT),VSDVI(1,1)),
     .  (RWK(NID3P+NRAD*2*NPLT),VECSAV(1,1))
C
C  LAST INDEX: NID3+NRAD*(3*NPLT)
C
      EQUIVALENCE
     .  (ESTIM(1),PDENA(1,1)),
     .  (SDVI(1),SIGMA(1,1)),
     .  (SDVI_BGK(1),SIGMA_BGK(1,1)),
     .  (SDVI_COP(1),SIGMA_COP(1,1)),
     .  (OUTAU(1),PDENAI(0,0)),
     .  (PLPRM(1),TEIN(1))
      SAVE
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        ISAVE=ISTRA
        IFIRST=1
      ENDIF
C
      IF (TRCPLT)
     .    WRITE (6,*) 'PLTEIR CALLED, ISTRA, XMCP: ',ISTRA,XMCP(ISTRA)
C
C  NULLPUNKT AUF DEM PAPIER
      X0PL=10.
      Y0PL=3.
C  ACHSENLAENGEN
      LENX=25.
      LENY=20.
C  ACHSENUNTERTEILUNG VORGEGEBEN?
C  NEIN!
      STPSZX=0.
      STPSZY=0.
      INTNRX=0
      INTNRY=0
C  ACHSE LOGARITHMISCH?
      LOGX=.FALSE.
C     LOGY VIA INPUT
C  LOG. ACHSE MIN
      MINLY=0
C  LOG. ACHSE MAX
C     MAXLY WERDEN BERECHNET IN ANPSGL
C  ZEICHNE NETZLINIEN EIN
      GRIDX=.TRUE.
      GRIDY=.TRUE.
C  MACHE GRADE GRENZEN, X-ACHSE (Y ACHSE, NUR WENN TALZMI=TALZMA=666.)
      FITX=.TRUE.
C
C
      IF ((NFILEN.EQ.1.OR.NFILEN.EQ.2).AND.IESTR.NE.ISTRA) THEN
        IF (XMCP(ISTRA).GT.1.) THEN
          IESTR=ISTRA
          CALL RSTRT(ISTRA,NSTRAI,NESTIM,NSDVI,ESTIM,SDVI,
     .               NSBGK,SDVI_BGK,NSCOP,SDVI_COP,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL SYMET(ESTIM,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .                 NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IF (XMCP(ISTRA).GT.1.) THEN
          IESTR=ISTRA
          CALL RSTRT(ISTRA,NSTRAI,NESTIM,NSDVI,ESTIM,SDVI,
     .               NSBGK,SDVI_BGK,NSCOP,SDVI_COP,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL SYMET(ESTIM,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .                 NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ENDIF
      ELSEIF (NFILEN.EQ.0.AND.ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
      ELSE
        WRITE (6,*) 'ERROR IN PLTEIR: DATA FOR STRATUM ISTRA= ',ISTRA
        WRITE (6,*) 'ARE NOT AVAILABLE. PLOTS ABANDONNED'
        RETURN
      ENDIF
C
      IF (ISTRA.EQ.0)
     .HEAD='SUM OVER STRATA
     .          '
      IF (ISTRA.NE.0) THEN
      HEAD='STRATUM NO.
     .          '
      WRITE (HEAD(13:15),'(I3)') ISTRA
      ENDIF
C
      HEAD0='VOLUME AVERAGED BACKGROUND TALLY, INPUT
     .           '
      HEAD1='DEFAULT VOLUME AVERAGED TALLY, TRACKLENGTH ESTIMATED
     .           '
      HEAD2='ADDITIONAL VOLUME AVERAGED TALLY, TRACKLENGTH ESTIMATED
     .           '
      HEAD3='ADDITIONAL VOLUME AVERAGED TALLY, COLLISION ESTIMATED
     .           '
      HEAD4='VOLUME AVERAGED TALLY, SNAPSHOT ESTIMATED
     .           '
      HEAD5='VOLUME AVERAGED TALLY, FOR COUPLING TO PLASMA CODE
     .           '
      HEAD6='BGK TALLY
     .           '
      HEAD7='ALGEBRAIC FUNCTION OF VOLUME AVERAGED TALLIES
     .           '
      HEAD8='RELATIVE STANDARD DEVIATION
     .           '
C
      IALG=0
C
C  .......................................
C
C   LOOP OVER NVOLPL
C  .......................................
C
      DO 10000 IBLD=1,NVOLPL
C
        IF (PLTL2D(IBLD).OR.PLTL3D(IBLD)) THEN
C
          DO 110 ICURV=1,NSPTAL(IBLD)
            JTAL=NPTALI(IBLD,ICURV)
C  REDO ALGEBRAIC TALLY IN CASE NFILEN=2 OR NFILEN=7
            IF (JTAL.EQ.NTALR.AND.IALG.EQ.0.AND.
     .          (NFILEN.EQ.2.OR.NFILEN.EQ.7)) THEN
              CALL ALGTAL
              IALG=1
              DO 105 IALV=1,NALVI
                CALL INTTAL (ALGV,VOL,IALV,NALV,NSBOX,ALGVI(IALV,ISTRA),
     .                       NR1ST,NP2ND,NT3RD,NBMLT)
105           CONTINUE
            ENDIF
            ITL=IABS(JTAL)
C  PLOT OUTPUT TALLIES ONLY FOR STRATA WITH TWO OR MORE HISTORIES
            IF (JTAL.GT.0.AND.XMCP(ISTRA).LE.1) GOTO 10000
C  PLOT INPUT TALLIES ONLY ONCE PER ITERATION
            IF (JTAL.LT.0.AND.ISAVE.NE.ISTRA) GOTO 10000
            TXHEAD=HEAD0
            IF (JTAL.GT.0) TXHEAD=HEAD1
            IF (JTAL.EQ.NTALA) TXHEAD=HEAD2
            IF (JTAL.EQ.NTALC) TXHEAD=HEAD3
            IF (JTAL.EQ.NTALT) TXHEAD=HEAD4
            IF (JTAL.EQ.NTALM) TXHEAD=HEAD5
            IF (JTAL.EQ.NTALB) TXHEAD=HEAD6
            IF (JTAL.EQ.NTALR) TXHEAD=HEAD7
            IF (TRCPLT) THEN
              CALL LEER(1)
              WRITE (6,*) 'PLOT REQUESTED FOR TALLY NO. ',JTAL
            ENDIF
C
            LPLTT2=.FALSE.
C
C .............................................
C
C  PUT TALLY ONTO ARRAY: VECTOR
C .............................................
C
C
            LSDVI(ICURV)=.FALSE.
            LPLOT2(ICURV)=.FALSE.
            ISPZ=ISPTAL(IBLD,ICURV)
            IF (JTAL.LT.0.) THEN
              NF=NFRSTP(ITL)
              IF (ISPZ.EQ.0) THEN
                DO 111 I=1,NRAD
111               VECTOR(I,ICURV)=0.
                DO 112 K=1,NF
                DO 112 I=1,NRAD
                  IINDEX=NADDP(ITL)*NRAD+(I-1)*NF+K
                  VECTOR(I,ICURV)=VECTOR(I,ICURV)+PLPRM(IINDEX)
112             CONTINUE
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NF) THEN
                DO 113 I=1,NRAD
                  IINDEX=NADDP(ITL)*NRAD+(I-1)*NF+ISPZ
                  VECTOR(I,ICURV)=PLPRM(IINDEX)
113             CONTINUE
              ELSE
                IF (TRCPLT) THEN
                  WRITE (6,*) 'SPECIES INDEX OUT OF RANGE '
                  WRITE (6,*) 'ICURV,ISPTAL(IBLD,ICURV) ',ICURV,ISPZ
                  WRITE (6,*) 'ALL PLOTS FOR THIS TALLY TURNED OFF '
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF
            ELSEIF (JTAL.GE.0) THEN
              NFT=NFSTVI(ITL)
              NF=NFIRST(ITL)
              IF (ISPZ.EQ.0) THEN
                DO 121 I=1,NRAD
121               VECTOR(I,ICURV)=0.
                DO 122 K=1,NFT
                  DO 122 I=1,NRAD
                    IINDEX=NADDV(ITL)*NRAD+(I-1)*NF+K
                    VECTOR(I,ICURV)=VECTOR(I,ICURV)+ESTIM(IINDEX)
122             CONTINUE
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NFT) THEN
                DO 125 I=1,NRAD
                  IINDEX=NADDV(ITL)*NRAD+(I-1)*NF+ISPZ
                  VECTOR(I,ICURV)=ESTIM(IINDEX)
125             CONTINUE
              ELSE
                IF (TRCPLT) THEN
                  WRITE (6,*) 'SPECIES INDEX OUT OF RANGE '
                  WRITE (6,*) 'ICURV,ISPTAL(IBLD,ICURV) ',ICURV,ISPZ
                  WRITE (6,*) 'ALL PLOTS FOR THIS TALLY TURNED OFF '
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF
C
              IF (PLTLER(IBLD)) THEN
C  CHECK IF STANDARD DEVIATION IS AVAILABLE FOR THIS TALLY
                DO 126 N=1,NSIGVI
                  IF (IIH(N).NE.JTAL) GOTO 126
                  IF (IGH(N).NE.ISPZ.AND.IGH(N).NE.0) GOTO 126
                  LSDVI(ICURV)=.TRUE.
                  DO 127 I=1,NRAD
                    VSDVI(I,ICURV)=SIGMA(N,I)
127               CONTINUE
126             CONTINUE
              ENDIF
C
            ENDIF
C
            DO 129 I=1,NRAD
              VECSAV(I,ICURV)=VECTOR(I,ICURV)
129         CONTINUE
C
110       CONTINUE
C
C ...................................
C                                   .
C    VECTOR(IC,ICURV) IS SET NOW    .
C ...................................
C
C
          LOGY=PLTLLG(IBLD)
C
          IF (PLTL2D(IBLD)) THEN
C
C  ............................
C
C   SET ABSCISSA FOR 2D PLOT
C  ............................
C
            IXSET2=0
            IF (LEVGEO.EQ.1.OR.LEVGEO.EQ.2) THEN
C   USE RADIAL SURFACE CENTERED GRID "RHOSRF" FOR EACH POLOIDAL POSITION
              IXTL2=NSURF
              DO 130 I=1,NR1ST
                XXP2D(I)=RHOSRF(I)
130           CONTINUE
              DO 131 J=2,NP2ND*NT3RD*NBMLT
                DO 131 I=1,NR1ST
                  XXP2D(I+(J-1)*NR1ST)=XXP2D(I)
131            CONTINUE
              DO 138 I=NSURF+1,NRAD
138             XXP2D(I)=0.
              IXSET2=1
            ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
C   USE PERPEND. ARCLENGTH "BGLP" IN CASE OF POLYGON GRID AND NLPOL
              IXTL2=NSURF
              DO 133 I=1,NR1ST
                DO 133 J=1,NP2ND
                  DO 133 K=1,NT3RD
                    IRAD=I+((J-1)+(K-1)*NP2T3)*NR1P2
                    XXP2D(IRAD)=BGLP(I,J)
133           CONTINUE
              DO 136 I=NSURF+1,NRAD
136             XXP2D(I)=0.
              IXSET2=1
            ELSE
C   NO 2D PLOTOPTIONS AVAILABLE
            ENDIF
C
            IF (IXSET2.NE.1) THEN
              WRITE (6,*) ' NO GRID SET FOR 2D PLOTTING '
              WRITE (6,*) ' IXSET2 = ',IXSET2
              WRITE (6,*) ' NO 2D PLOTTING IS DONE '
              GOTO 1000
            ENDIF
C
            XMI=XXP2D(NPLIN2(IBLD,1))*(1.+1.E-6)
            XMA=XXP2D(NPLOT2(IBLD,1))/(1.+1.E-6)
C
C  IN CASE OF LSMOT2, SET ZONE CENTERED ABSCISSA
C  GRID FROM SURFACE CENTERED  GRID  "X"
C
            IF (LSMOT2(IBLD)) THEN
              DO 137 J=1,NRAD-1
                XXP2D(J)=(XXP2D(J)+XXP2D(J+1))*0.5
137           CONTINUE
            ENDIF
C
            DO 140 ICURV=1,NSPTAL(IBLD)
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IR1(ICURV)=NPLIN2(IBLD,ICURV)
              IR2(ICURV)=NPLOT2(IBLD,ICURV)
              IRS(ICURV)=NPLDL2(IBLD,ICURV)
C
              YMNLG2(ICURV)=1.D60
              YMXLG2(ICURV)=-1.D60
              IF (JTAL.GT.0.) THEN
                I0=0
                IF (NFRSTI(ITL).GT.1) I0=1
                INDX=NADDI(ITL)*NSTRAP+NFRSTI(ITL)*ISTRA+ISPZ+I0
                IF (OUTAU(INDX).EQ.0.) THEN
                  IF (TRCPLT) THEN
                    WRITE (6,*) 'TALLY NO. ',JTAL,' CURVE NO. ',ICURV
                    WRITE (6,*) 'NOT PLOTTED BECAUSE'
                    WRITE (6,*) 'ZERO INTEGRAL (OUTAU(INDX)=0.) '
                    WRITE (6,*) 'INDX,NADDI(JTAL),NFRSTI(JTAL),I0'
                    WRITE (6,*)  INDX,NADDI(ITL),NFRSTI(ITL),I0
                  ENDIF
                  YMN2(ICURV)=0.
                  YMX2(ICURV)=0.
                  YMNLG2(ICURV)=0.
                  YMXLG2(ICURV)=0.
                  GOTO 140
                ENDIF
              ENDIF
              LPLOT2(ICURV)=.TRUE.
              LPLTT2=.TRUE.
              I1=IR1(ICURV)
              I2=IR2(ICURV)
              IF (LSMOT2(IBLD)) I2=I2-1
              IS=IRS(ICURV)
C
C YMNLG2, YMXLG2: REAL MAX/MIN, FOR LEGENDE ON 2D PLOT ONLY
              DO 141 I=I1,I2,IS
                YMNLG2(ICURV)=MIN(YMNLG2(ICURV),VECTOR(I,ICURV))
141           CONTINUE
              DO 142 I=I1,I2,IS
                YMXLG2(ICURV)=MAX(YMXLG2(ICURV),VECTOR(I,ICURV))
142           CONTINUE
C
C YMN2, YMX2: FOR AXIS
              FITY=.TRUE.
              IF (TALZMI(IBLD).NE.666.) THEN
                IF (.NOT.LOGY) FITY=.FALSE.
                YMN2(ICURV)=TALZMI(IBLD)
                DO 143 I=1,NRAD
                  VECTOR(I,ICURV)=MAX(YMN2(ICURV),VECTOR(I,ICURV))
143             CONTINUE
                IF (LOGY) YMN2(ICURV)=YMN2(ICURV)*(1.+1.E-6)
              ELSE
                YMN2(ICURV)=YMNLG2(ICURV)
              ENDIF
C
              IF (TALZMA(IBLD).NE.666.) THEN
                IF (.NOT.LOGY) FITY=.FALSE.
                YMX2(ICURV)=TALZMA(IBLD)
                DO 144 I=1,NRAD
                  VECTOR(I,ICURV)=MIN(YMX2(ICURV),VECTOR(I,ICURV))
144             CONTINUE
                IF (LOGY) YMX2(ICURV)=YMX2(ICURV)/(1.+1.E-6)
              ELSE
                YMX2(ICURV)=YMXLG2(ICURV)
              ENDIF
140         CONTINUE
C
C  PLOT ALL CURVES REQUESTED FROM THIS TALLY INTO ONE PICTURE
            DO 150 ICURV=1,NSPTAL(IBLD)
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IF (ISPZ.EQ.0) THEN
                TXSPEC(ICURV)='SUM OVER SPECIES        '
                IF (JTAL.LT.0) TXUNIT(ICURV)=TXTPUN(1,ITL)
                IF (JTAL.LT.0) TXTALL(ICURV)=TXTPLS(1,ITL)
                IF (JTAL.GE.0) TXUNIT(ICURV)=TXTUNT(1,ITL)
                IF (JTAL.GE.0) TXTALL(ICURV)=TXTTAL(1,ITL)
              ELSE
                IF (JTAL.LT.0) THEN
                  TXTALL(ICURV)=TXTPLS(ISPZ,ITL)
                  TXSPEC(ICURV)=TXTPSP(ISPZ,ITL)
                  TXUNIT(ICURV)=TXTPUN(ISPZ,ITL)
                ELSE
                  TXTALL(ICURV)=TXTTAL(ISPZ,ITL)
                  TXSPEC(ICURV)=TXTSPC(ISPZ,ITL)
                  TXUNIT(ICURV)=TXTUNT(ISPZ,ITL)
                ENDIF
              ENDIF
150         CONTINUE
            IERR=0
            CALL PLTTLY (XXP2D,VECTOR,VSDVI,YMN2,YMX2,
     .             IR1,IR2,IRS,
     .             NSPTAL(IBLD),TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .             LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,LHIST2(IBLD),IERR)
            IF (IERR.GT.0) GOTO 1000
            IF (TRCPLT) THEN
              WRITE (6,*) '2D PLOT FOR TALLY NO. ',JTAL,' DONE'
              WRITE (6,*) 'XMIN= ',XMI,' XMAX= ',XMA
              DO 160 ICURV=1,NSPTAL(IBLD)
                IF (LPLOT2(ICURV))
     .            WRITE (6,*) 'ICURV= ',ICURV,
     .                        ' YMIN= ',YMNLG2(ICURV),
     .                        ' YMAX= ',YMXLG2(ICURV),
     .                        ' LSDVI= ',LSDVI(ICURV)
160           CONTINUE
            ENDIF
C
          ENDIF
C
1000      CONTINUE
C
C   3D PLOT GRID
C
          IF (PLTL3D(IBLD)) THEN
C
            DO 1040 ICURV=1,NSPTAL(IBLD)
              DO 1035 I=1,NRAD
                 VECTOR(I,ICURV)=VECSAV(I,ICURV)
1035          CONTINUE
C  SYMMETRY CONDITION AT POLAR ANGLE THETA=YIA AND THETA=2*PI+YIA
              IF (LEVGEO.EQ.2.AND.IYTL3.EQ.NP2ND) THEN
                DO 1036 I=1,IXTL3
1036              VECTOR(I+NP2NDM*NR1ST,ICURV)=VECTOR(I,ICURV)
              ENDIF
1040        CONTINUE
C
C SET QUASIRECTANGULAR PLOT GRIDS XXP3D (IX), IX=1,IXTL3
C                             AND YYP3D (IY), IY=1,IYTL3
C FOR EACH 3D PICTURE
C
            IXSET3=0
            IYSET3=0
            IF (LEVGEO.EQ.1) THEN
              IXTL3=NR1ST
              DO 218 I=1,IXTL3
218             XXP3D(I)=RHOSRF(I)
              IXSET3=1
              IF (NLTOR.AND.NLTRZ) THEN
                IYTL3=NT3RD
                DO 220 I=1,IYTL3
220               YYP3D(I)=ZSURF(I)
                DO I=1,NR1ST
                  DO J=1,NT3RD
                    XPOL(I,J)=RHOSRF(I)
                    YPOL(I,J)=ZSURF(J)
                  ENDDO
                ENDDO
                IYSET3=1
              ELSEIF (NLPOL) THEN
                IYTL3=NP2ND
                DO 221 I=1,IYTL3
221               YYP3D(I)=PSURF(I)
                DO I=1,NR1ST
                  DO J=1,NP2ND
                    XPOL(I,J)=RHOSRF(I)
                    YPOL(I,J)=PSURF(J)
                  ENDDO
                ENDDO
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.2) THEN
C
              IXTL3=NR1ST
              DO 223 I=1,IXTL3
223             XXP3D(I)=RHOSRF(I)
              IXSET3=1
              IF (NLTOR.AND.NLTRZ) THEN
                IYTL3=NT3RD
                DO 226 I=1,IYTL3
226               YYP3D(I)=ZSURF(I)
                IYSET3=1
              ELSEIF (NLPOL) THEN
                IYTL3=NP2ND
                DO 225 I=1,IYTL3-1
225               YYP3D(I)=0.5*(PSURF(I+1)+PSURF(I))
                YYP3D(NP2ND)=PSURF(1)+PI2A
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.3) THEN
C
              IXTL3=NR1ST
              DO 228 IX=1,IXTL3
228             XXP3D(IX)=IX
              IXSET3=1
              IYTL3=NP2ND
              DO 230 IX=1,IYTL3
230             YYP3D(IX)=IX
              IYSET3=1
C
            ELSEIF (LEVGEO.EQ.4) THEN
C
              IXSET3=1
              IYSET3=1
            ELSE
              WRITE (6,*) '3D PLOT OPTION TO BE WRITTEN, LEVGEO '
              CALL EXIT
            ENDIF
C
C  LOOP ICURV=1,....
C
            ICINC=1
            IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) ICINC=2
            DO 1160 ICURV=1,NSPTAL(IBLD),ICINC
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IF (ISPZ.EQ.0) THEN
                TXSPEC(1)='SUM OVER SPECIES        '
                IF (JTAL.LT.0) TXUNT1=TXTPUN(1,ITL)
                IF (JTAL.LT.0) TXTLL1=TXTPLS(1,ITL)
                IF (JTAL.GE.0) TXUNT1=TXTUNT(1,ITL)
                IF (JTAL.GE.0) TXTLL1=TXTTAL(1,ITL)
              ELSE
                IF (JTAL.LT.0) THEN
                  TXTLL1=TXTPLS(ISPZ,ITL)
                  TXSPC1=TXTPSP(ISPZ,ITL)
                  TXUNT1=TXTPUN(ISPZ,ITL)
                ELSE
                  TXTLL1=TXTTAL(ISPZ,ITL)
                  TXSPC1=TXTSPC(ISPZ,ITL)
                  TXUNT1=TXTUNT(ISPZ,ITL)
                ENDIF
              ENDIF
C
              IF (IXSET3+IYSET3.LT.2) THEN
                WRITE (6,*) ' NO GRIDS SET FOR 3D PLOTTING '
                WRITE (6,*) ' IXSET3,IYSET3 = ',IXSET3,IYSET3
                WRITE (6,*) ' NO 3D PLOTTING IS DONE '
                GOTO 10000
              ENDIF
C
              LINLOG=PLTLLG(IBLD)
              TMIN=TALZMI(IBLD)
              TMAX=TALZMA(IBLD)
C
1200          CONTINUE
C
C
C  CONTOUR PLOTS
              IF (LCNTR3(IBLD)) THEN
                CALL ISOLNE (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS PLOTS
              ELSEIF (LRAPS3(IBLD)) THEN
                CALL RPSCOL (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  3D HISTOGRAM
              ELSEIF (LHIST3(IBLD)) THEN
                CALL PL3DPG (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,TALW1(IBLD),TALW2(IBLD),
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  3D SURFACE PLOTS, IN CUBE
              ELSEIF (LSMOT3(IBLD)) THEN
                CALL PLOT3D (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,TALW1(IBLD),TALW2(IBLD),
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  VECTOR FIELD PLOT
              ELSEIF (LVECT3(IBLD)) THEN
                CALL VECLNE (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS VECTOR PLOTS
              ELSEIF (LRPVC3(IBLD)) THEN
                CALL RPSVEC (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
              ELSE
                WRITE (6,*) 'NO 3D PLOT OPTION FOR IBLD= ',IBLD
                GOTO 10000
              ENDIF
C
C   PLOT STANDARD DEVIATION PROFILE FOR THIS TALLY, IF REQUESTED
C
              IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) GOTO 1160
              IF (PLTLER(IBLD).AND.LSDVI(ICURV)) THEN
                TXUNT1='%                       '
                TXHEAD=HEAD8
                DO 1222 I=1,NRAD
                  VECTOR(I,ICURV)=VSDVI(I,ICURV)
1222            CONTINUE
                LINLOG=.FALSE.
                TMIN=0.
                TMAX=100.
                LSDVI(ICURV)=.FALSE.
                GOTO 1200
              ENDIF
C
1160        CONTINUE
C  LOOP ICURV FINISHED
          ENDIF
C
        ELSEIF (NSPTAL(IBLD).GT.0) THEN
          WRITE (6,*) 'PLOT REQUEST FOR TALLY NO. ',JTAL,' BUT NEITHER'
          WRITE (6,*) 'PLTL2D NOR PLTL3D TRUE. NO PLOT FOR THIS TALLY'
        ENDIF
C
10000 CONTINUE
C  LOOP IBLD FINISHED
C
      RETURN
      END
C
      SUBROUTINE PLTTLY (X,Y,VBAR,YMN,YMX,IR1,IR2,IRS,NKURV,TXTTAL,
     .                   TXTSPC,TXTUNT,TXTRUN,TXHEAD,
     .                   LBAR,XMI,XMA,YMNLG,YMXLG,LPLOT,LHIST,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CCONA'
      INCLUDE 'CPLMSK'
C
      DIMENSION X(*),Y(NRAD,*),VBAR(NRAD,*),YMN(*),YMX(*),YMNLG(*),
     .          YMXLG(*)
      DIMENSION IR1(*),IR2(*),IRS(*)
      CHARACTER*(*) TXTTAL(*),TXTSPC(*),TXTUNT(*),TXTRUN,TXHEAD
      CHARACTER*10 CHR
      LOGICAL LBAR(*),LPLOT(*),LHIST
C
      IKURV=0
      DO 1 I=1,NKURV
        IF (LPLOT(I)) IKURV=IKURV+1
1     CONTINUE
      IF (IKURV.EQ.0) RETURN
C
      CALL GRNXTB(1)
      CALL GRSCLC (0.,0.,39.5,28.7)
      CALL GRSCLV (0.,0.,39.5,28.7)
      IT=LEN(TXTRUN)
      CALL GRTXT (1.,27.5,IT,TXTRUN)
      IT=LEN(TXHEAD)
      CALL GRTXT (1.,26.75,IT,TXHEAD)
      CALL GRTXT (1.,26.,15,'MIN. ABSZISSA =')
      WRITE (CHR,'(1P,E10.3)') XMI
      CALL GRTXTC (10,CHR)
      CALL GRTXT (1.,25.5,15,'MAX. ABSZISSA =')
      WRITE (CHR,'(1P,E10.3)') XMA
      CALL GRTXTC (10,CHR)
      YA=24.75
      IFI=0
      DO 2 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 2
        IF (IFI.EQ.0) THEN
          IFI=1
          IT=LEN(TXTTAL(I))
          CALL GRTXT (1.5,24.75,IT,TXTTAL(I))
        ENDIF
        YA=YA-1.
        IT=LEN(TXTSPC(I))
        ISY=I+1
        CALL GRJMPS (0.5,sngl(YA),ISY)
        CALL GRTXT (1.5,sngl(YA),IT,TXTSPC(I))
        YA=YA-0.5
        IT=LEN(TXTUNT(I))
        CALL GRTXT (1.5,sngl(YA),IT,TXTUNT(I))
        YA=YA-0.5
        CALL GRTXT (1.5,sngl(YA),11,'MAX. VALUE =')
        WRITE (CHR,'(1P,E10.3)') YMXLG(I)
        CALL GRTXTC (10,CHR)
        YA=YA-0.5
        CALL GRTXT (1.5,sngl(YA),11,'MIN. VALUE =')
        WRITE (CHR,'(1P,E10.3)') YMNLG(I)
        CALL GRTXTC (10,CHR)
2     CONTINUE
C
      MINX=XMI
      MAXX=XMA
      MINY=1.D30
      MAXY=-1.D30
      DO 3 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 3
        MINY=MIN(MINY,real(YMN(I)))
        MAXY=MAX(MAXY,real(YMX(I)))
3     CONTINUE
C
      IF (LOGY) THEN
        MAXY=MAX(sngl(1.D-48),MAXY)
        MINY=MAX(MINY,MAXY/sngl(1.D12))
      ENDIF
C
C  PLOT X AND Y AXIS
C
      CALL PLTMSK(IERR)
      IF (IERR.GT.0) RETURN
C
C  PREPARE ARRAYS FOR PLOTTING ON LOG. SCALE
C
      IF (LOGY) THEN
        YMINY=MAX(10.D0**MINLY,10.D0**(MAXLY-12))
        DO 20 ICURV=1,NKURV
          IF (.NOT.LPLOT(ICURV)) GOTO 20
          DO 21 I=IR1(ICURV),IR2(ICURV)-1,IRS(ICURV)
21          Y(I,ICURV)=LOG10(MAX(YMINY,Y(I,ICURV)))
20      CONTINUE
      ENDIF
C
C  PLOT !
C
      DO 50 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 50
        I1=IR1(I)
        IS=IRS(I)
        I2=IR2(I)-IS
C
        IF (LHIST) THEN
C  PLOT LINES
          IF (LOGY) THEN
            YMY=MINLY
          ELSE
            YMY=MINY
          ENDIF
          CALL GRJMP(sngl(X(I1)),sngl(YMY))
          CALL GRDRW(sngl(X(I1)),sngl(Y(I1,I)))
          CALL GRDRW(sngl(X(I1+IS)),sngl(Y(I1,I)))
          DO 30 J=I1+IS,I2,IS
            CALL GRDRW (sngl(X(J)),sngl(Y(J,I)))
30          CALL GRDRW (sngl(X(J+IS)),sngl(Y(J,I)))
          CALL GRDRW(sngl(X(I2+IS)),sngl(YMY))
C  PLOT SYMBOLS
          ISY=I+1
          NP=(I2-I1+1)/IS
          NPS=MAX0(NP/7,1)
          DO 31 J=I1,I2,IS*NPS
31          CALL GRJMPS (sngl(0.5*(X(J)+X(J+IS))),sngl(Y(J,I)),ISY)
C
        ELSEIF (.NOT.LHIST) THEN
C  PLOT LINES
          CALL GRJMP(sngl(X(I1)),sngl(Y(I1,I)))
          DO 33 J=I1+IS,I2,IS
33          CALL GRDRW (sngl(X(J)),sngl(Y(J,I)))
C  PLOT SYMBOLS
          ISY=I+1
          NP=(I2-I1+1)/IS
          NPS=MAX0(NP/7,1)
          DO 35 J=I1,I2,IS*NPS
35          CALL GRJMPS (sngl(X(J)),sngl(Y(J,I)),ISY)
        ENDIF
C
C  PLOT ERROR BARS
          IF (LBAR(I)) THEN
            DO 40 J=I1,I2,IS
              FP=1.+VBAR(J,I)*0.01
              FM=1.-VBAR(J,I)*0.01
              IF (LOGY) THEN
                AA=10.**Y(J,I)
                ST1=LOG10(MAX(1.D-30,AA*FM))
                ST2=LOG10(MAX(1.D-30,AA*FP))
                ST1=MAX(DBLE(MINLY),MIN(DBLE(MAXLY),ST1))
                ST2=MAX(DBLE(MINLY),MIN(DBLE(MAXLY),ST2))
              ELSE
                ST1=Y(J,I)*FM
                ST2=Y(J,I)*FP
                ST1=MAX(dble(MINY),MIN(dble(MAXY),ST1))
                ST2=MAX(dble(MINY),MIN(dble(MAXY),ST2))
              ENDIF
            IF (LHIST) THEN
              CALL GRJMP (sngl(0.5*(X(J)+X(J+IS))),sngl(ST1))
              CALL GRDRW (sngl(0.5*(X(J)+X(J+IS))),sngl(ST2))
            ELSE
              CALL GRJMP (sngl(X(J)),sngl(ST1))
              CALL GRDRW (sngl(X(J)),sngl(ST2))
            ENDIF
40        CONTINUE
        ENDIF
50    CONTINUE
      CALL GRCHRC (0.3,0.,16)
C
      RETURN
      END
C
C----------------------------------------------------------------------*
C SUBROUTINE PLTMSK                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE PLTMSK(IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  PLOTPROGRAMM ZUR ERSTELLUNG EINES ZEICHENFELDES.
C  DIE EINGABEWERTE GEBEN AN WO SICH DAS FELD IM BILD BEFINDEN SOLL UND
C  WELCHE ART VON X- UND Y-ACHSE GEWUENSCHT WIRD.
C  DAS HEISST OB MIT LOGARITHMISCHER ODER LINEARER ACHSENEINTEILUNG,
C  MIT GITTERLINIEN ODER NUR MIT KURZEN MARKIERUNGEN VERSEHEN
C  ODER BESCHRIFTUNG VOM PROGRAMM ANGEPASST ODER SELBST VORGEGEBEN.
C
C  EINGABEDATEN:  COMMONBLOECKE ORIGIN,XAXES,YAXES
C
      INCLUDE 'CPLMSK'
C
      CALL GRSCLC(sngl(X0PL),sngl(Y0PL),sngl(X0PL+LENX),
     .            sngl(Y0PL+LENY))
C  PLOT X AXIS
      IAX=1
      CALL PLTAXI(IERR,IAX)
C  PLOT Y AXIS
      IAX=2
      CALL PLTAXI(IERR,IAX)
      IF (IERR.GT.0) RETURN
C
      IF (LOGX) THEN
         XMIN=REAL(MINLX)
         XMAX=REAL(MAXLX)
      ELSE
         XMIN=MINX
         XMAX=MAXX
      ENDIF
C
      IF (LOGY) THEN
         YMIN=REAL(MINLY)
         YMAX=REAL(MAXLY)
      ELSE
         YMIN=MINY
         YMAX=MAXY
      ENDIF
C
      CALL GRSCLV(sngl(XMIN),sngl(YMIN),sngl(XMAX),sngl(YMAX))
C
      RETURN
      END
C----------------------------------------------------------------------*
C SUBROUTINE PLTAXI                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE PLTAXI(IERR,IAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  PLOTPROGRAMM ZUM ZEICHNEN EINER X- ODER Y- ACHSE.
C  DIE EINGABEWERTE GEBEN AN, UM WELCHE ACHSE ES SICH HANDELT,
C  OB MIT LOGARITHMISCHER ODER LINEARER ACHSENEINTEILUNG,
C  MIT GITTERLINIEN ODER NUR MIT KURZEN MARKIERUNGEN VERSEHEN
C  ODER BESCHRIFTUNG VOM PROGRAMM ANGEPASST ODER SELBST VORGEGEBEN
C  UND WIE LANG DIE ANDERE ACHSE IST.
C
C  EINGABE :  COMMON XAXES,YAXES
C  IAX=1: PLOTTE X ACHSE
C  IAX=2: PLOTTE Y ACHSE
C
      INCLUDE 'CPLMSK'
      INTEGER       FCTR
      CHARACTER     CPARAM*7,CFCTR*3,CZHNLB*3,CEINLB*1
C
C  ZEICHNEN EINER X-ACHSE
C
      IF (IAX.EQ.1) THEN
C
C  LINEARE X-ACHSENEINTEILUNG
C
         IF (.NOT.LOGX) THEN
            IERR=0
            IF (FITX) CALL ANPSG(MINX,MAXX,INTNRX,STPSZX,IERR)
            IF (IERR.GT.0) RETURN
            CALL GRSPTS(20)
            CALL GRDSH(1.,0.,1.)
            CALL GRSCLV(MINX,0.,MAXX,LENY)
C
C  GITTERLINIEN BZW. MARKIERUNGEN AN DER X-ACHSE
C
            DO 5 J=0,INTNRX
               T=MINX+J*STPSZX
               CALL GRJMP(sngl(T),-0.1)
               IF (GRIDX) THEN
                  IF (J.NE.0.AND.J.NE.INTNRX) THEN
                    CALL GRSPTS(16)
                    CALL GRDSH(0.2,0.5,0.2)
                  ENDIF
                  CALL GRDRW(sngl(T),LENY)
                  CALL GRSPTS(20)
                  CALL GRDSH(1.,0.,1.)
               ELSE
                  CALL GRDRW(sngl(T),0.)
               ENDIF
    5       CONTINUE
            IF (.NOT.GRIDX) THEN
               CALL GRJMP(MINX,0.)
               CALL GRDRW(MINX,LENY)
               CALL GRJMP(MAXX,0.)
               CALL GRDRW(MAXX,LENY)
            ENDIF
C
C  LABELS AN DER X-ACHSE
C
            CALL GRCHRC(0.3,0.,16)
            DO 10 J=2,INTNRX+1,2
               FCTR=IEXP10(DBLE(MINX+MAXX)*0.5)
               PARAM=(MINX+(J-1)*STPSZX)/10.**FCTR
               PARAM1=PARAM*100000.
               IPARAM=NINT(PARAM1)
               PARAM=REAL(IPARAM)/100000.
               WRITE(CPARAM,'(F7.3)') PARAM
               IL=7
               T=MINX+(J-1)*STPSZX-(MAXX-MINX)*0.8/LENX
               CALL GRTXT(sngl(T),-0.5,IL,CPARAM)
   10       CONTINUE
            T=MAXX-1.6*(MAXX-MINX)/LENX
            IF (ABS(FCTR).GE.10) THEN
              WRITE(CFCTR,'(I3)') FCTR
            ELSE
              WRITE(CFCTR,'(I2)') FCTR
              CFCTR(3:3)=' '
            ENDIF
            IF (FCTR.LT.0.) CALL GRTXT(sngl(T),-1.0,8,'*10**'//CFCTR)
            IF (FCTR.GT.0.) CALL GRTXT(sngl(T),-1.0,7,
     .                                 '*10**'//CFCTR(2:3))
C
C  LOGARITHMISCHE X-ACHSENEINTEILUNG
C
          ELSE
            IERR=0
            IF (FITX) CALL ANPSGL(MINX,MAXX,MINLX,MAXLX,IERR)
            IF (IERR.GT.0) RETURN
            CALL GRSCLV(REAL(MINLX),0.,REAL(MAXLX),LENY)
C
C  GITTERLINIEN BZW. MARKIERUNGEN AN DER X-ACHSE
C
            EXPR=LENX/(MAXLX-MINLX)
            IF (EXPR.GE.6) THEN
               DO 15 RI=MINLX,MAXLX-1
                  DO 15 RJ=1,9
                     ARG=RJ*10.**RI
                     T=LOG10(ARG)
                     IF (RJ.LT.2.) THEN
                        CALL GRSPTS(20)
                     ELSE
                        CALL GRSPTS(16)
                     ENDIF
                     CALL GRJMP(sngl(T),-0.1)
                     IF (GRIDX) THEN
                        IF (RI.NE.MINLX.OR.RJ.NE.1) THEN
                          CALL GRDSH(0.2,0.5,0.2)
                          CALL GRSPTS(16)
                        ENDIF
                        CALL GRDRW(sngl(T),LENY)
                        CALL GRDSH(1.,0.,1.)
                     ELSE
                        CALL GRDRW(sngl(T),0.)
                     ENDIF
   15          CONTINUE
            ELSE IF(EXPR.GE.4) THEN
               DO 20 RI=MINLX,MAXLX-1
                  CALL GRSPTS(20)
                  CALL GRJMP(sngl(RI),-0.1)
                  IF (GRIDX) THEN
                     IF (RI.NE.MINLX) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                     CALL GRDRW(sngl(RI),LENY)
                     CALL GRDSH(1.,0.,1.)
                     CALL GRSPTS(20)
                  ELSE
                     CALL GRDRW(sngl(RI),0.)
                  ENDIF
                  DO 20 RJ=2,8,2
                     CALL GRSPTS(16)
                     ARG=RJ*10.**RI
                     T=LOG10(ARG)
                     CALL GRJMP(sngl(T),-0.1)
                     IF (GRIDX) THEN
                        CALL GRDSH(0.2,0.5,0.2)
                        CALL GRDRW(sngl(T),LENY)
                        CALL GRDSH(1.,0.,1.)
                     ELSE
                        CALL GRDRW(sngl(T),0.)
                     ENDIF
   20          CONTINUE
            ELSE IF (EXPR.GE.1) THEN
               DO 25 RI=MINLX,MAXLX-1
                  CALL GRSPTS(20)
                  CALL GRJMP(sngl(RI),-0.1)
                  IF (GRIDX) THEN
                     IF (RI.NE.MINLX) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                     CALL GRDRW(sngl(RI),LENY)
                     CALL GRDSH(1.,0.,1.)
                     CALL GRSPTS(20)
                  ELSE
                     CALL GRDRW(sngl(RI),0.)
                  ENDIF
                  DO 25 RJ=2,5,3
                     CALL GRSPTS(16)
                     ARG=RJ*10.**RI
                     T=LOG10(ARG)
                     CALL GRJMP(sngl(T),-0.1)
                     IF (GRIDX) THEN
                        CALL GRDSH(0.2,0.5,0.2)
                        CALL GRDRW(sngl(T),LENY)
                        CALL GRDSH(1.,0.,1.)
                     ELSE
                        CALL GRDRW(sngl(T),0.)
                     ENDIF
   25          CONTINUE
            ELSE IF (EXPR.GE.0) THEN
               DO 30 RI=MINLX,MAXLX-1
                  CALL GRSPTS(20)
                  CALL GRJMP(sngl(RI),-0.1)
                  IF (GRIDX) THEN
                     IF (RI.NE.MINLX) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                     CALL GRDRW(sngl(RI),LENY)
                     CALL GRDSH(1.,0.,1.)
                     CALL GRSPTS(20)
                  ELSE
                       CALL GRDRW(sngl(RI),0.)
                  ENDIF
   30          CONTINUE
            ENDIF
            CALL GRJMP(REAL(MAXLX),-0.1)
C
            IF (GRIDX) THEN
               CALL GRDRW(REAL(MAXLX),LENY)
            ELSE
               CALL GRDRW(REAL(MAXLX),0.)
            ENDIF
C
            IF (.NOT.GRIDX) THEN
               CALL GRJMP(REAL(MINLX),0.)
               CALL GRDRW(REAL(MINLX),LENY)
               CALL GRJMP(REAL(MAXLX),0.)
               CALL GRDRW(REAL(MAXLX),LENY)
            ENDIF
C
C  10-ER LABELS AN DER X-ACHSE
C
            CALL GRCHRC(0.3,0.,20)
            DO 35 I=MINLX+1,MAXLX
               T=I-0.3*(MAXLX-MINLX)/LENX
               CALL GRTXT(sngl(T),-0.90,2,'10')
  35        CONTINUE
C
            DO 40 I=MINLX+1,MAXLX
               WRITE(CZHNLB,'(I3)') I
               IF (ABS(I).LT.10) THEN
                  CALL GRTXT(REAL(I),-0.6,2,CZHNLB(2:3))
               ELSE
                  CALL GRTXT(REAL(I),-0.60,3,CZHNLB)
               ENDIF
  40        CONTINUE
C
C  EINFACHE LABELS AN DER X-ACHSE
C
            CALL GRCHRC(0.3,0.,18)
            EXPR=LENX/(MAXLX-MINLX)
            IF (EXPR.GE.8) THEN
               DO 45 RI=MINLX,MAXLX-1
                  DO 45 J=2,9
                     WRITE(CEINLB,'(I1)') J
                     ARG=REAL(J)*10.**RI
                     T=LOG10(ARG)-0.05*(MAXLX-MINLX)/LENX
                     CALL GRTXT(sngl(T),-0.5,1,CEINLB)
   45          CONTINUE
            ELSE IF (EXPR.GE.4) THEN
               DO 50 RI=MINLX,MAXLX-1
                  DO 50 J=2,8,2
                     WRITE(CEINLB,'(I1)') J
                     ARG=REAL(J)*10.**RI
                     T=LOG10(ARG)-0.05*(MAXLX-MINLX)/LENX
                     CALL GRTXT(sngl(T),-0.5,1,CEINLB)
   50          CONTINUE
            ELSE IF (EXPR.GE.2) THEN
               DO 55 RI=MINLX,MAXLX-1
                  DO 55 J=2,5,3
                     WRITE(CEINLB,'(I1)') J
                     ARG=REAL(J)*10.**RI
                     T=LOG10(ARG)-0.05*(MAXLX-MINLX)/LENX
                     CALL GRTXT(sngl(T),-0.5,1,CEINLB)
   55          CONTINUE
            ENDIF
C
C KEINE LABELS FUER
C       LENX/(MAXLX-MINLX)>=0
C
         ENDIF
C
C  ZEICHNEN EINER Y-ACHSE
C
      ELSEIF (IAX.EQ.2) THEN
C
C  LINEARE  Y-ACHSENEINTEILUNG
C
         IF (.NOT.LOGY) THEN
            IERR=0
            IF (FITY) THEN
              CALL ANPSG(MINY,MAXY,INTNRY,STPSZY,IERR)
            ELSE
              INTNRY=10
              STPSZY=(MAXY-MINY)/10.
            ENDIF
            IF (IERR.GT.0) RETURN
            CALL GRSPTS(20)
            CALL GRDSH(1.,0.,1.)
            CALL GRSCLV(0.,MINY,LENX,MAXY)
C
C   GITTERLINIEN BZW. MARKIERUNGEN AN DER Y-ACHSE
C
            DO 110 J=0,INTNRY
               T=MINY+J*STPSZY
               CALL GRJMP(-0.1,sngl(T))
               IF (GRIDY) THEN
                  IF (J.NE.0.AND.J.NE.INTNRY) THEN
                    CALL GRSPTS(16)
                    CALL GRDSH(0.2,0.5,0.2)
                  ENDIF
                  CALL GRDRW(LENX,sngl(T))
                  CALL GRSPTS(20)
                  CALL GRDSH(1.,0.,1.)
               ELSE
                  CALL GRDRW(0.,sngl(T))
               ENDIF
  110       CONTINUE
            IF (.NOT.GRIDY) THEN
               CALL GRJMP(0.,MINY)
               CALL GRDRW(LENX,MINY)
               CALL GRJMP(0.,MAXY)
               CALL GRDRW(LENX,MAXY)
            ENDIF
C
C  LABELS AN DER Y-ACHSE
C
            CALL GRCHRC(0.3,90.,18)
            DO 115 J=2,INTNRY+1,2
               FCTR=IEXP10(DBLE(MINY+MAXY)*0.5)
               PARAM=(MINY+(J-1)*STPSZY)/10.**FCTR
               PARAM1=PARAM*100000.
               IPARAM=NINT(PARAM1)
               PARAM=REAL(IPARAM)/100000.
               WRITE(CPARAM,'(F7.3)') PARAM
               IL=7
               T=MINY+(J-1)*STPSZY-(MAXY-MINY)*0.8/LENY
               CALL GRTXT(-0.45,sngl(T),IL,CPARAM)
  115      CONTINUE
           T=MAXY-1.6*(MAXY-MINY)/LENY
           IF (ABS(FCTR).GE.10) THEN
             WRITE(CFCTR,'(I3)') FCTR
           ELSE
             WRITE(CFCTR,'(I2)') FCTR
             CFCTR(3:3)=' '
           ENDIF
           IF (FCTR.LT.0.) CALL GRTXT(-0.85,sngl(T),8,'*10**'//CFCTR)
           IF (FCTR.GT.0.) CALL GRTXT(-0.85,sngl(T),7,
     .                                '*10**'//CFCTR(2:3))
C
C  LOGARITHMISCHE Y-ACHSENEINTEILUNG
C
         ELSE
           IERR=0
           IF (FITY) CALL ANPSGL(MINY,MAXY,MINLY,MAXLY,IERR)
           IF (IERR.GT.0) RETURN
           CALL GRSCLV(0.,REAL(MINLY),LENX,REAL(MAXLY))
C
C  GITTERLINIEN BZW. MARKIERUNGEN AN DER Y-ACHSE
C
           EXPR=LENY/(MAXLY-MINLY)
           IF (EXPR.GE.6) THEN
              DO 120 RI=MINLY,MAXLY-1
                 DO 120 RJ=1,9
                    ARG=RJ*10.**RI
                    T=LOG10(ARG)
                    IF (RJ.LT.2.) THEN
                       CALL GRSPTS(20)
                    ELSE
                       CALL GRSPTS(16)
                    ENDIF
                    CALL GRJMP(-0.1,sngl(T))
                    IF (GRIDY) THEN
                       IF (RI.NE.MINLY.OR.RJ.NE.1) THEN
                         CALL GRDSH(0.2,0.5,0.2)
                         CALL GRSPTS(16)
                       ENDIF
                       CALL GRDRW(LENX,sngl(T))
                       CALL GRDSH(1.,0.,1.)
                    ELSE
                       CALL GRDRW(0.,sngl(T))
                    ENDIF
  120         CONTINUE
           ELSE IF (EXPR.GE.4) THEN
              DO 125 RI=MINLY,MAXLY-1
                 CALL GRSPTS(20)
                 CALL GRJMP(-0.1,sngl(RI))
                 IF (GRIDY) THEN
                     IF (RI.NE.MINLY) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                    CALL GRDRW(LENX,sngl(RI))
                    CALL GRSPTS(20)
                    CALL GRDSH(1.,0.,1.)
                 ELSE
                    CALL GRDRW(0.,sngl(RI))
                 ENDIF
                 DO 125 RJ=2,8,2
                    CALL GRSPTS(16)
                    ARG=RJ*10.**RI
                    T=LOG10(ARG)
                    CALL GRJMP(-0.1,sngl(T))
                    IF (GRIDY) THEN
                       CALL GRDSH(0.2,0.5,0.2)
                       CALL GRDRW(LENX,sngl(T))
                       CALL GRDSH(1.,0.,1.)
                    ELSE
                       CALL GRDRW(0.,sngl(T))
                    ENDIF
  125         CONTINUE
           ELSE IF (EXPR.GE.1) THEN
              DO 130 RI=MINLY,MAXLY-1
                 CALL GRSPTS(20)
                 CALL GRJMP(-0.1,sngl(RI))
                 IF (GRIDY) THEN
                     IF (RI.NE.MINLY) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                    CALL GRDRW(LENX,sngl(RI))
                    CALL GRSPTS(20)
                    CALL GRDSH(1.,0.,1.)
                 ELSE
                    CALL GRDRW(0.,sngl(RI))
                 ENDIF
                 DO 130 RJ=2,5,3
                    CALL GRSPTS(18)
                    ARG=RJ*10.**RI
                    T=LOG10(ARG)
                    CALL GRJMP(-0.1,sngl(T))
                    IF (GRIDY) THEN
                       CALL GRDSH(0.2,0.5,0.2)
                       CALL GRDRW(LENX,sngl(T))
                       CALL GRDSH(1.,0.,1.)
                    ELSE
                       CALL GRDRW(0.,sngl(T))
                    ENDIF
  130         CONTINUE
           ELSE IF (EXPR.GE.0) THEN
              DO 135 RI=MINLY,MAXLY-1
                 CALL GRSPTS(20)
                 CALL GRJMP(-0.1,sngl(RI))
                 IF (GRIDY) THEN
                    IF (RI.NE.MINLY) THEN
                      CALL GRSPTS(16)
                      CALL GRDSH(0.2,0.5,0.2)
                    ENDIF
                    CALL GRDRW(LENX,sngl(RI))
                    CALL GRSPTS(20)
                    CALL GRDSH(1.,0.,1.)
                 ELSE
                    CALL GRDRW(0.,sngl(RI))
                 ENDIF
  135         CONTINUE
           ENDIF
C
           CALL GRJMP(-0.1,REAL(MAXLY))
           IF (GRIDY) THEN
              CALL GRDRW(LENX,REAL(MAXLY))
           ELSE
              CALL GRDRW(0.,REAL(MAXLY))
           ENDIF
           IF (.NOT.GRIDY) THEN
              CALL GRJMP(0.,REAL(MINLY))
              CALL GRDRW(LENX,REAL(MINLY))
              CALL GRJMP(0.,REAL(MAXLY))
              CALL GRDRW(LENX,REAL(MAXLY))
           ENDIF
C
C  10-ER LABELS AN DER Y-ACHSE
C
           CALL GRCHRC(0.3,0.,20)
           DO 140 I=MINLY+1,MAXLY
              T=I-0.1*(MAXLY-MINLY)/LENY
              CALL GRTXT(-1.20,sngl(T),2,'10')
  140      CONTINUE
C
           DO 145 I=MINLY+1,MAXLY
              WRITE(CZHNLB,'(I3)') I
              IF (ABS(I).LT.10) THEN
                 CALL GRTXT(-0.85,REAL(I)+0.05,2,CZHNLB(2:3))
              ELSE
                 CALL GRTXT(-0.85,REAL(I)+0.05,3,CZHNLB)
              ENDIF
  145      CONTINUE
C
C  EINFACHE LABELS AN DER Y-ACHSE
C
           CALL GRCHRC(0.3,0.,18)
           EXPR=LENY/(MAXLY-MINLY)
           IF (EXPR.GE.8) THEN
              DO 150 RI=MINLY,MAXLY-1
                 DO 150 J=2,9
                    WRITE(CEINLB,'(I1)') J
                    ARG=REAL(J)*10.**RI
                    T=LOG10(ARG)-0.1*(MAXLY-MINLY)/LENY
                    CALL GRTXT(-0.4,sngl(T),1,CEINLB)
  150         CONTINUE
           ELSE IF (EXPR.GE.4) THEN
              DO 155 RI=MINLY,MAXLY-1
                 DO 155 J=2,8,2
                    WRITE(CEINLB,'(I1)') J
                    ARG=REAL(J)*10.**RI
                    T=LOG10(ARG)-0.1*(MAXLY-MINLY)/LENY
                    CALL GRTXT(-0.4,sngl(T),1,CEINLB)
  155         CONTINUE
           ELSE IF (EXPR.GE.2) THEN
              DO 160 RI=MINLY,MAXLY-1
                 DO 160 J=2,5,3
                    WRITE(CEINLB,'(I1)') J
                    ARG=REAL(J)*10.**RI
                    T=LOG10(ARG)-0.1*(MAXLY-MINLY)/LENY
                    CALL GRTXT(-0.4,sngl(T),1,CEINLB)
  160         CONTINUE
           ENDIF
C
C KEINE LABELS FUER
C   LENY/(MAXLY-MINLY) >= 0
C
         ENDIF
      ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------*
C SUBROUTINE ANPSGL                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE ANPSGL(MIN,MAX,MINL,MAXL,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  DIESES PROGRAMM ERHAELT ZWEI INTERVALLGRENZEN (MIN,MAX) UND
C  BERECHNET WELCHE ZEHNERPOTENZ ALS LINKE INTERVALLGRENZE UND
C  WELCHE ALS RECHT IN FRAGE KOMMT (MINL,MAXL).
C  DAS INTERVALL (MINL,MAXL) EIGNET SICH DANN ZUR LOGARTHMISCHEN
C  SKALIERUNG EINER KOORDINATENACHSE.
C
C  EINGABE :  MIN   DEC DBLE(6)
C             MAX   DEC DBLE(6)
C  AUSGABE :  MINL  BIN FIXED(15)
C             MAXL  BIN FIXED(15)                                  *
C
      real    MIN,MAX
C
      IF (MIN.LE.0.OR.MAX.LE.0) THEN
         WRITE(6,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(6,*)  '***    MIN<=0  ODER  MAX<=0    ***'
         WRITE(6,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MIN.GT.MAX) THEN
         WRITE(6,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(6,*)  '***        MIN  >  MAX         ***'
         WRITE(6,*)  'MIN= ',MIN,' MAX= ',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MAX-MIN.GE.0) THEN
         MINL=IEXP10(DBLE(MIN))
         MAXL=IEXP10(DBLE(MAX))+1
      ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------*
C SUBROUTINE ANPSG                                                     *
C----------------------------------------------------------------------*
      SUBROUTINE ANPSG(MIN,MAX,INTNR,STPSZ,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  DIESES PROGRAMM ERHAELT ZWEI INTERVALLGRENZEN (MIN,MAX) UND
C  FORMT DIESE MEIST "KRUMMEN" ZAHLEN UM IN SOLCHE, DIE SICH
C  ZUR ACHSENBESCHRIFTUNG EIGNEN.
C  INTNR GIBT GIBT ANZAHL DER TEILINTERVALLE ZWISCHEN MIN UND MAX
C  AN, DIE SICH AUS DER UMFORMUNG ERGEBEN.
C  STPSZ IST DIE LAENGE EINES SOLCHEN TEILINTERVALLES.
C
C  EINGABE :  MIN   DEC DBLE(6)
C             MAX   DEC DBLE(6)
C  AUSGABE :  MIN   DEC DBLE(6)
C             MAX   DEC DBLE(6)
C             INTNR BIN FIXED(15)
C             STPSZ DEC DBLE(6)
C
      REAL    MAX,MIN,STPSZ
      DOUBLE PRECISION    LKS,RTS,MDDL,STPS2,DIFF
C
      DIFF=MAX-MIN
      IF (DIFF.EQ.0) THEN
         IPT=IEXP10(DBLE(MAX))
         MIN=MIN-10.**(IPT-1)
         MAX=MAX+10.**(IPT-1)
         DIFF=MAX-MIN
      ENDIF
      IF (DIFF.GT.0) THEN
         STPS2=DIFF/10
         IPT=IEXP10(STPS2)
         STPS2=STPS2/10.**IPT
         STPS2=AINT(STPS2+0.5)*10.**IPT
         IPT=IEXP10(DIFF)
         MDDL=(MAX+MIN)*0.5/10.**(IPT-1)
         MDDL=AINT(MDDL)*10.**(IPT-1)
         LKS=MDDL
         RTS=MDDL
         INTNR=0
    5    IF (LKS.GT.MIN) THEN
            INTNR=INTNR+1
            LKS=LKS-STPS2
            GOTO 5
         ENDIF
   10    IF (RTS.LT.MAX) THEN
            INTNR=INTNR+1
            RTS=RTS+STPS2
            GOTO 10
         ENDIF
         IPT=IEXP10(LKS)
         MIN=LKS+SIGN(1.D0,LKS)*10.D0**(-14+IPT)
         IPT=IEXP10(RTS)
         MAX=RTS+SIGN(1.D0,RTS)*10.D0**(-14+IPT)
         IPT=IEXP10(STPS2)
         STPSZ=STPS2+10.**(-14+IPT)
      ELSE
         WRITE(6,*)  '------------------------------------'
         WRITE(6,*)  'PARAMETERFEHLER IN ANPSG: MIN > MAX'
         WRITE(6,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE PLOT3D (ARR,IBLD,ICURV,
     .                   IX,IY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,LOGL,
     .                   ZMA,ZMI,W1,W2,
     .                   HEAD,RUNID,TXHEAD,TRC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ON INPUT:  LOGL: USE LOG SCALE FOR ORDINATE
C             ARR(IJ),I=1,IX-1,J=1,IY-1 ,IJ=I+(J-1)*NR1ST
C                          ARRAY TO BE PLOTTED
C             XX(I),I=1,IX X-GRID BOUNDARIES
C             YY(J),J=1,IY Y GRID BOUNDARIES
C  ARR(IJ) IS THEN SET ONTO 2D ARRAY FALT(I,J)
C
C
C  LPOLAR: R-THETA CO-ORDINATES
C  LKARTH: X-Y     CO-ORDINATES
C
C  FALT-->XZY (3,...)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CPLOT'
      PARAMETER (NYDIM=N2NDPLG)
      DIMENSION XX(*),YY(*),ARR(*)
      REAL*4
     .   FALT(N1ST,NYDIM),X(N1ST,NYDIM),Y(N1ST,NYDIM),Z2(N1ST,NYDIM)
      PARAMETER (LAR=46*N1ST*NYDIM)
      REAL*4 AR(LAR),EXT(3,3),VALU(3,2),DCM,YH
      REAL*4 YHLF,XHLF,FMIN,FMAX,REMIN,REMAX,XMIN,XMAX,YMIN,YMAX
      LOGICAL LOGL,TRC
      CHARACTER*72 TEXT1
      CHARACTER*24 TEXT2,TEXT3
      CHARACTER*72 HEAD,RUNID,TXHEAD
      CHARACTER*17 CH
      CHARACTER*20 CHAXS(3)
C
C
      ier=1
      IF (LEVGEO.LE.3.AND.LEVGEO.GT.1.AND.LPTOR3(IBLD)) THEN
C
        IXM=IX-1
        IYM=IY-1
        IZ=IYM*NR1ST
C
        DO 1 I=1,N1ST
          DO 1 J=1,NYDIM
            FALT(I,J)=-75.75E20
1       CONTINUE
C
        IF (LOGL) THEN
          DO 3 J=1,IZ
            ARR(J)=LOG10(MAX(1.D-48,ARR(J)))
3         CONTINUE
        ENDIF
C
C  SET ONTO 2D ARRAY FOR PLOTTING
C
        DO 20 I=1,IX
          DO 20 J=1,IYM
            DXXX=ARR(I+(J-1)*NR1ST)
            FALT(I,J)=DXXX
            X(I,J)=XPOL(I,J)
            Y(I,J)=YPOL(I,J)
20      CONTINUE
C
        IXM=IX-1
        IYM=IY-1
C
C     SEARCH FOR MINIMUM AND MAXIMUM AND REPLACE, IF REQUIRED
C
        IZN=IXM*IYM
        FMIN=FALT(1,1)
        FMAX=FALT(1,1)
        XMIN=X(1,1)
        XMAX=X(1,1)
        YMIN=Y(1,1)
        YMAX=Y(1,1)
        DO 25 J=1,IXM
        DO 25 I=1,IYM
          FMIN=MIN(FMIN,FALT(J,I))
          FMAX=MAX(FMAX,FALT(J,I))
          XMIN=MIN(XMIN,X(J,I))
          XMAX=MAX(XMAX,X(J,I))
          YMIN=MIN(YMIN,Y(J,I))
          YMAX=MAX(YMAX,Y(J,I))
25      CONTINUE
C
        REMIN=ZMI
        REMAX=ZMA
        IF (LOGL) THEN
          REMIN=LOG10(MAX(1.D-48,ZMI))
          REMAX=LOG10(MAX(1.D-48,ZMA))
        ENDIF
        IF (ZMI.EQ.666.) REMIN=FMIN
        IF (ZMA.EQ.666.) REMAX=FMAX
        IF (REMIN.GE.REMAX) THEN
          REMIN=REMIN-1.
          REMAX=REMAX+1.
        ENDIF
C
        DO 30 J=1,IXM
        DO 30 I=1,IYM
          FALT(J,I)=MIN(REMAX,FALT(J,I))
          FALT(J,I)=MAX(REMIN,FALT(J,I))
30      CONTINUE
C
C  PLOT SMOOTH SURFACE
C
        DCM=MAX(ABS(XMAX-XMIN),ABS(YMAX-YMIN))
        XHLF=0.5*(XMAX+XMIN)
        YHLF=0.5*(YMAX+YMIN)
        XMINN=XHLF-0.5*DCM
        XMAXN=XHLF+0.5*DCM
        YMINN=YHLF-0.5*DCM
        YMAXN=YHLF+0.5*DCM
C
C     NORMIEREN DER WERTE
C
        DO IX=1,IXM
          DO IY=1,IYM
            X(IX,IY) = (X(IX,IY)-XMINN)/(XMAXN-XMINN)
            Y(IX,IY) = (Y(IX,IY)-YMINN)/(YMAXN-YMINN)
            FALT(IX,IY) = (FALT(IX,IY)-REMIN)/(REMAX-REMIN)
            Z2(IX,IY) = 0.
          ENDDO
        ENDDO
C
        CALL GRNXTB (1)
        CALL GRSCLC(10.,3.,34.,27.)
        CALL GRSCLV(2.,2.,26.,26.)
        CALL GR3DIM(LAR,IER)
        CALL GR3NT1(AR,IER,n1st,X,Y,Z2,IXM,1,IYM,1,1,1)
        CALL GR3NT1(AR,IER,n1st,X,Y,FALT,IXM,1,IYM,1,1,2)
        CALL GR3EXT(AR,IER,EXT)
        VALU(1,1)=XMIN
        VALU(1,2)=XMAX
        VALU(2,1)=YMIN
        VALU(2,2)=YMAX
        VALU(3,1)=REMIN
        VALU(3,2)=REMAX
        CHAXS(1) = ' '
        CHAXS(2) = ' '
        CHAXS(3) = ' '
        CALL GR3AXS(AR,IER,EXT,VALU,CHAXS,.FALSE.,4,1)
        CALL GR3ROT(AR,IER,'Z',SNGL(W1),'X',SNGL(W2),'Y',0.0)
        CALL GR3PLO(AR,IER,'HID')
C
      ELSE
        WRITE (6,*) 'INVALID OPTION IN PLOT3D  '
        WRITE (6,*) '3D-PLOT ABANDONED  '
        RETURN
      ENDIF
C
C     WRITE TEXT ONTO THE PLOT
C
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,YH,72,RUNID)
      YH=26.75
      CALL GRTXT (1.,YH,72,HEAD)
      YH=26.00
      CALL GRTXT (1.,YH,72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,YH,10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,YH-0.5,10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,YH-1.,10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,YH-2.,10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') FMAX
      CALL GRTXT (1.,YH-2.5,10,CH)
      CALL GRTXT (1.,YH-3.,10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') FMIN
      CALL GRTXT (1.,YH-3.5,10,CH)
C
      RETURN
      END
C
      SUBROUTINE PL3DPG (ARR,IBLD,ICURV,
     .                   IX,IY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,LOGL,
     .                   ZMA,ZMI,W1,W2,
     .                   HEAD,RUNID,TXHEAD,TRC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CGEOM'
      INCLUDE 'CGRID'
      INCLUDE 'CPOLYG'
      INCLUDE 'CPLOT'
      DIMENSION XX(*),YY(*),ARR(*)
      PARAMETER (NYDIM=N2NDPLG)
      PARAMETER (LAR=46*128*128)
      REAL*4 AR(LAR),EXT(3,3),VALU(3,2)
      REAL*4 XYZ(3,128,128)
      real*4 yh
      LOGICAL LOGL,TRC
      CHARACTER*72 TEXT1,HEAD,RUNID,TXHEAD
      CHARACTER*24 TEXT2,TEXT3
      CHARACTER*17 CH
      CHARACTER*20 CHAXS(3)
C
      DO 1 I=1,3
      DO 1 J=1,128
      DO 1 K=1,128
1       XYZ(I,J,K)=-75.75E20
C
      IF (LEVGEO.LE.1.OR.LEVGEO.GT.3.OR..NOT.LPTOR3(IBLD)) THEN
        WRITE (6,*) 'PLOTOPTION NOT READY. RETURN FROM PL3DPG '
        RETURN
      ENDIF
C
      IPAN=1
      IPEN=NP2ND
C
C  SEARCH MINIMA AND MAXIMA OF INDEPENDENT VARIABLES X AND Y
C
      XMI=1.D60
      XMA=-1.D60
      YMI=1.D60
      YMA=-1.D60
      RMI=1.D60
      RMA=-1.D60
      DO 10 I=1,NR1ST
        DO 10 J=1,NPPLG
          DO 10 K=NPOINT(1,J),NPOINT(2,J)
            XMI=MIN(XMI,XPOL(I,K))
            XMA=MAX(XMA,XPOL(I,K))
            YMI=MIN(YMI,YPOL(I,K))
            YMA=MAX(YMA,YPOL(I,K))
C  SEARCH MINIMA AND MAXIMA OF DEPENDENT VARIABLE Z=ARR
            IF (K.LT.NPOINT(2,J).AND.I.LT.NR1ST) THEN
              IR=I+(K-1)*NR1ST
              IF (LOGL) THEN
                RMI=MIN(RMI,MAX(1.D-48,ARR(IR)))
                RMA=MAX(RMA,MAX(1.D-48,ARR(IR)))
                ARR(IR)=LOG10(MAX(1.D-48,ARR(IR)))
              ELSE
                RMI=MIN(RMI,ARR(IR))
                RMA=MAX(RMA,ARR(IR))
              ENDIF
            ENDIF
10    CONTINUE
C
      REMIN=ZMI
      REMAX=ZMA
      IF (LOGL) THEN
        IF (ZMI.EQ.666.) THEN
          REMIN=LOG10(RMI)
        ELSE
          REMIN=LOG10(MAX(1.D-48,ZMI))
        ENDIF
        IF (ZMA.EQ.666.) THEN
          REMAX=LOG10(RMA)
        ELSE
          REMAX=LOG10(MAX(1.D-48,ZMA))
        ENDIF
      ELSE
        IF (ZMI.EQ.666.) REMIN=RMI
        IF (ZMA.EQ.666.) REMAX=RMA
      ENDIF
C
      IF (TRC) THEN
        WRITE (6,*) 'PL3DPG:  ,IPAN,IPEN,XMI,XMA,YMI,YMA,REMIN,REMAX'
        WRITE (6,*) '        ',IPAN,IPEN,XMI,XMA,YMI,YMA,REMIN,REMAX
      ENDIF
C
      IF (ABS(REMAX-REMIN)/MAX(REMAX,1.D-30).LT.0.01)
     .    REMAX=REMIN+0.01*REMIN*SIGN(1.D0,REMIN)
C
C
      DX=MAX(ABS(XMA-XMI),ABS(YMA-YMI))*0.5
      XMT=0.5*(XMA+XMI)
      YMT=0.5*(YMA+YMI)
      XMINN=XMT-DX
      XMAXN=XMT+DX
      YMINN=YMT-DX
      YMAXN=YMT+DX
C
C     SET PLOT PARAMETER
C
C  DRDMPA(15) FOR EACH VALID PART, SEE BELOW
      CALL GRNXTB (1)
      CALL GRSCLC(10.,3.,34.,27.)
      CALL GRSCLV(2.,2.,26.,26.)
      CALL GR3DIM(LAR,IER)
C
      DO 30 K=1,NPPLG
        IPX=1
        DO 25 I=1,NR1ST-1
          IPY=1
          DO 20 J=NPOINT(1,K),NPOINT(2,K)-1
            IF (IPX+2.GT.128) GOTO 999
            IF (IPY+2.GT.128) GOTO 999
            XYZ(1,IPX+1,IPY+1)=XPOL(I,J)
            XYZ(1,IPX+1,IPY+2)=XPOL(I,J+1)
            XYZ(1,IPX+2,IPY+1)=XPOL(I+1,J)
            XYZ(1,IPX+2,IPY+2)=XPOL(I+1,J+1)
C
            XYZ(2,IPX+1,IPY+1)=YPOL(I,J)
            XYZ(2,IPX+1,IPY+2)=YPOL(I,J+1)
            XYZ(2,IPX+2,IPY+1)=YPOL(I+1,J)
            XYZ(2,IPX+2,IPY+2)=YPOL(I+1,J+1)
C
            IR=I+(J-1)*NR1ST
            AAR=ARR(IR)
            AAR=MAX(MIN(REMAX,AAR),REMIN)
            XYZ(3,IPX+1,IPY+1)=AAR
            XYZ(3,IPX+1,IPY+2)=AAR
            XYZ(3,IPX+2,IPY+1)=AAR
            XYZ(3,IPX+2,IPY+2)=AAR
            IPY=IPY+2
20        CONTINUE
          IPX=IPX+2
25      CONTINUE
C
        IF (IPX+1.GT.128) GOTO 999
        IF (IPY+1.GT.128) GOTO 999
        DO 26 I=2,IPY
          XYZ(1,1,I)=XYZ(1,2,I)
          XYZ(2,1,I)=XYZ(2,2,I)
          XYZ(3,1,I)=REMIN
          XYZ(1,IPX+1,I)=XYZ(1,IPX,I)
          XYZ(2,IPX+1,I)=XYZ(2,IPX,I)
          XYZ(3,IPX+1,I)=REMIN
26      CONTINUE
        DO 27 I=2,IPX
          XYZ(1,I,1)=XYZ(1,I,2)
          XYZ(2,I,1)=XYZ(2,I,2)
          XYZ(3,I,1)=REMIN
          XYZ(1,I,IPY+1)=XYZ(1,I,IPY)
          XYZ(2,I,IPY+1)=XYZ(2,I,IPY)
          XYZ(3,I,IPY+1)=REMIN
27      CONTINUE
        XYZ(1,1,1)=XYZ(1,2,2)
        XYZ(2,1,1)=XYZ(2,2,2)
        XYZ(3,1,1)=REMIN
        XYZ(1,IPX+1,1)=XYZ(1,IPX,2)
        XYZ(2,IPX+1,1)=XYZ(2,IPX,2)
        XYZ(3,IPX+1,1)=REMIN
        XYZ(1,IPX+1,IPY+1)=XYZ(1,IPX,IPY)
        XYZ(2,IPX+1,IPY+1)=XYZ(2,IPX,IPY)
        XYZ(3,IPX+1,IPY+1)=REMIN
        XYZ(1,1,IPY+1)=XYZ(1,2,IPY)
        XYZ(2,1,IPY+1)=XYZ(2,2,IPY)
        XYZ(3,1,IPY+1)=REMIN
        DO I=1,IPX+1
          DO J=1,IPY+1
            XYZ(1,I,J)=(XYZ(1,I,J)-XMINN)/(XMAXN-XMINN)
            XYZ(2,I,J)=(XYZ(2,I,J)-YMINN)/(YMAXN-YMINN)
            XYZ(3,I,J)=(XYZ(3,I,J)-REMIN)/(REMAX-REMIN)
          ENDDO
        ENDDO
        CALL GR3NET(AR,IER,128,XYZ,IPX+1,1,IPY+1,1,1,2)
30    CONTINUE
      CALL GR3EXT(AR,IER,EXT)
      VALU(1,1)=XMI
      VALU(1,2)=XMA
      VALU(2,1)=YMI
      VALU(2,2)=YMA
      VALU(3,1)=REMIN
      VALU(3,2)=REMAX
      CHAXS(1) = ' '
      CHAXS(2) = ' '
      CHAXS(3) = ' '
      CALL GR3AXS(AR,IER,EXT,VALU,CHAXS,.FALSE.,4,1)
      CALL GR3ROT(AR,IER,'Z',SNGL(W1),'X',SNGL(W2),'Y',0.0)
      CALL GR3PLO(AR,IER,'HID')
C
C     WRITE TEXT AND MEAN VALUE ONTO THE PLOT
C
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,YH,72,RUNID)
      YH=26.75
      CALL GRTXT (1.,YH,72,HEAD)
      YH=26.00
      CALL GRTXT (1.,YH,72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,YH,10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,YH-0.5,10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,YH-1.,10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,YH-2.,10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') RMA
      CALL GRTXT (1.,YH-2.5,10,CH)
      CALL GRTXT (1.,YH-3.,10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') RMI
      CALL GRTXT (1.,YH-3.5,10,CH)
C
      RETURN
999   CONTINUE
      WRITE (6,*) 'NOT ENOUGH STORAGE FOR 3D HISTOGRAM PLOT'
      WRITE (6,*) 'REDUCE PLOT AREA '
      WRITE (6,*) 'PLOT ABANDONNED'
      RETURN
      END
C
C
      SUBROUTINE ZYLPLN (ZX0,ZY0,ZZ0,ZVX,ZVY,ZVZ,RZYL,JS,NZAD,NINNE,NIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
      INCLUDE 'CLGIN'
      INCLUDE 'CCONA'
C
      DIMENSION AFF(3,3),AFFI(3,3),PHIAN(20,9),PHIEN(20,9),ANGLE(19),
     .          IPAR(20),TAR(20),AL(4),AR(4)
C
      PHID=180.-ACOS(1./SQRT(2.))*RADDEG
C     PHID=0.
      ITEST=0
      WRITE (6,*) ' ZYLPLN CALLED FOR JS = ',JS
      WRITE (6,*) ' ZX0,ZY0,ZZ0 ',ZX0,ZY0,ZZ0
      WRITE (6,*) ' ZVX,ZVY,ZVZ ',ZVX,ZVY,ZVZ
C  WINKEL ZWISCHEN (ZVX,ZVY,ZVZ) UND (0,0,1)
      COSD=ZVZ/SQRT(ZVX*ZVX+ZVY*ZVY+ZVZ*ZVZ)
      SIND=SQRT(1.-COSD*COSD)
      WRITE (6,*) ' COSD,SIND ',COSD,SIND
C  DREHACHSE (0,0,1)X(ZVX,ZVY,ZVZ)
      BETC=SQRT(ZVX*ZVX+ZVY*ZVY)
      C1=-ZVY/(BETC+EPS60)
      C2=ZVX/(BETC+EPS60)
      C3=0.
C  DREHMATRIX NACH KORN & KORN
4711  CONTINUE
      EC=1.-COSD
      AFF(1,1)=COSD+EC*C1*C1
      AFF(1,2)=EC*C1*C2-C3*SIND
      AFF(1,3)=EC*C1*C3+C2*SIND
      AFF(2,1)=EC*C1*C2+C3*SIND
      AFF(2,2)=COSD+EC*C2*C2
      AFF(2,3)=EC*C2*C3-C1*SIND
      AFF(3,1)=EC*C3*C1-C2*SIND
      AFF(3,2)=EC*C3*C2+C1*SIND
      AFF(3,3)=COSD+EC*C3*C3
C
C  BERECHNE AFF**-1
      CALL INVERT(AFF,AFFI)
C
      V1=ZVX*AFFI(1,1)+ZVY*AFFI(2,1)+ZVZ*AFFI(3,1)
      V2=ZVX*AFFI(1,2)+ZVY*AFFI(2,2)+ZVZ*AFFI(3,2)
      V3=ZVX*AFFI(1,3)+ZVY*AFFI(2,3)+ZVZ*AFFI(3,3)
      WRITE (6,*) ' V1,V2,V3 ',V1,V2,V3
      IF (ABS(V1)+ABS(V2).GT.EPS10.OR.V3.LT.0.) THEN
        IF (ITEST.EQ.0) THEN
          ITEST=ITEST+1
          C1=-C1
          C2=-C2
          C3=-C3
          GOTO 4711
        ELSE
          WRITE (6,*) ' NONSENSE IN ZYLPLN   JS = ',JS
          WRITE (6,*) ' PLOT OF THIS SURFACE ABANDONNED '
          RETURN
        ENDIF
      ENDIF
C
      WRITE (6,*) ' VOR ROTADD ALIMS,XLIMS,YLIMS,ZLIMS '
      WRITE (6,'(1X,1P,4E12.4)') (ALIMS(JS,I),XLIMS(JS,I),YLIMS(JS,I),
     .                          ZLIMS(JS,I),I=1,ILIN(JS))
      CALL XSHADD (-ZX0,JS,JS)
      CALL YSHADD (-ZY0,JS,JS)
      CALL ZSHADD (-ZZ0,JS,JS)
      CALL ROTADD (AFF,AFFI,JS,JS)
      WRITE (6,*) ' NACH ROTADD ALIMS,XLIMS,YLIMS,ZLIMS '
      WRITE (6,'(1X,1P,4E12.4)') (ALIMS(JS,I),XLIMS(JS,I),YLIMS(JS,I),
     .                          ZLIMS(JS,I),I=1,ILIN(JS))
C
C  SUCHE T-INTERVALL
C
      TMAX=1.D60
      TMIN=-1.D60
      DO 100 I=1,ILIN(JS)
        IF (ABS(XLIMS(JS,I)).LT.EPS10.AND.ABS(YLIMS(JS,I)).LT.EPS10.AND.
     .      ABS(ZLIMS(JS,I)).LT.EPS10) THEN
          WRITE (6,*) ' NONSENSE IN ZYLPLN WITH SURFACE ',JS
          WRITE (6,*) ' PLOT OF SURFACE ABORTED '
          CALL ROTADD (AFFI,AFF,JS,JS)
          CALL XSHADD (ZX0,JS,JS)
          CALL YSHADD (ZY0,JS,JS)
          CALL ZSHADD (ZZ0,JS,JS)
          RETURN
        ENDIF
C
C
C  WINKEL DER EBENE MIT DER Z=0 EBENE
        BPLN=SQRT(XLIMS(JS,I)**2+YLIMS(JS,I)**2+ZLIMS(JS,I)**2)
        CPHI=ZLIMS(JS,I)/BPLN
C  EBENE PARALLEL ZUM ZYLINDER?
        IF (ABS(CPHI).LT.EPS10) GOTO 100
        HYP=RZYL/CPHI
        DLT=SQRT(HYP*HYP-RZYL*RZYL)
C  SCHNITT DER EBENE MIT (0,0,1)
        T=-ALIMS(JS,I)/ZLIMS(JS,I)
        TTEST=T-1.
C  FESTSTELLEN, OB (0,0,TTEST) EIN GUELTIGER PUNKT IST
        E=ZLIMS(JS,I)*TTEST+ALIMS(JS,I)
        IF (E.LE.0.) THEN
C  EBENE IST OBERE BEGRENZUNG
          IF (T+DLT.LT.TMAX) THEN
            TMAX=T+DLT
            TMXD=T-DLT
            IMAX=I
          ENDIF
        ELSE
C  EBENE IST UNTERE BEGRENZUNG
          IF (T-DLT.GT.TMIN) THEN
            TMIN=T-DLT
            TMND=T+DLT
            IMIN=I
          ENDIF
        ENDIF
C
100   CONTINUE
C
      IF (TMIN*TMAX.GT.0.) THEN
        IF (ABS(TMIN).GT.ABS(TMAX)) THEN
           HILF=TMIN
           HILFD=TMND
           IHILF=IMIN
           TMIN=TMAX
           TMND=TMXD
           IMIN=IMAX
           TMAX=HILF
           TMXD=HILFD
           IMAX=IHILF
         ENDIF
      ENDIF
C
      WRITE (6,*) ' TMIN,TMAX ',TMIN,TMAX
      WRITE (6,*) ' TMND,TMXD ',TMND,TMXD
      WRITE (6,*) ' IMIN,IMAX ',IMIN,IMAX
      DLT=(TMAX-TMIN)/(NZAD+1)
C
      DO 200 IZ=1,NZAD+2
        T=TMIN+(IZ-1)*DLT
        IF (IZ.EQ.1) T=TMND
        IF (IZ.EQ.NZAD+2) T=TMXD
        IANG=0
C
C  SUCHE DIE WINKELBEREICHE
C
        DO 110 I=1,ILIN(JS)
          B0=ALIMS(JS,I)+ZLIMS(JS,I)*T
          B1=XLIMS(JS,I)
          B2=YLIMS(JS,I)
          IANG=IANG+2
          CALL SECANG (B0,B1,B2,RZYL,ANGLE(IANG-1),ANGLE(IANG))
110     CONTINUE
        WRITE (6,*) ' ANGLE VOR SORT ',(ANGLE(I),I=1,IANG)
C
C  SORTIERE WINKEL
120     ISORT=0
        DO 115 I=1,IANG-1
          IF (ANGLE(I+1).LT.ANGLE(I)) THEN
            PHIH=ANGLE(I)
            ANGLE(I)=ANGLE(I+1)
            ANGLE(I+1)=PHIH
            ISORT=ISORT+1
          ENDIF
115     CONTINUE
        IF (ISORT.GT.0) GOTO 120
        WRITE (6,*) ' ANGLE NACH SORT ',(ANGLE(I),I=1,IANG)
C
        IANG=IANG+1
        ANGLE(IANG)=ANGLE(1)
C
        IPART=0
        DO 130 I=1,IANG-1
          IF (ANGLE(I+1).LT.ANGLE(I)) ANGLE(I+1)=ANGLE(I+1)+360.
          IF (ABS(ANGLE(I+1)-ANGLE(I)).LT.EPS10) GOTO 130
          PHI=0.5*(ANGLE(I)+ANGLE(I+1))*DEGRAD
          X=RZYL*COS(PHI)
          Y=RZYL*SIN(PHI)
          Z=T
C         WRITE (6,*) ' ANGLE(I),(I+1) ',ANGLE(I),ANGLE(I+1)
C         WRITE (6,*) ' PHI,X,Y,Z, ',PHI,X,Y,Z
          DO 125 IE=1,ILIN(JS)
            TEST=ALIMS(JS,IE)+XLIMS(JS,IE)*X+YLIMS(JS,IE)*Y+
     .           ZLIMS(JS,IE)*Z
C           WRITE (6,*) ' IE,TEST ',IE,TEST
            IF (TEST.GT.0.) GOTO 130
125       CONTINUE
          IPART=IPART+1
          PHIAN(IZ,IPART)=ANGLE(I)
          PHIEN(IZ,IPART)=ANGLE(I+1)
130     CONTINUE
        TAR(IZ)=T
        IPAR(IZ)=IPART
C
        WRITE (6,*) ' IZ,IPAR,TAR ',IZ,IPAR(IZ),TAR(IZ)
        WRITE (6,'(1X,1P,2E12.4)')(PHIAN(IZ,JP),PHIEN(IZ,JP),JP=1,IPART)
C
200   CONTINUE
C
C
      CALL ROTADD (AFFI,AFF,JS,JS)
      CALL XSHADD (ZX0,JS,JS)
      CALL YSHADD (ZY0,JS,JS)
      CALL ZSHADD (ZZ0,JS,JS)
C
C
      AL(1)=ALIMS(JS,IMIN)
      AL(2)=XLIMS(JS,IMIN)
      AL(3)=YLIMS(JS,IMIN)
      AL(4)=ZLIMS(JS,IMIN)
C
      AR(1)=ALIMS(JS,IMAX)
      AR(2)=XLIMS(JS,IMAX)
      AR(3)=YLIMS(JS,IMAX)
      AR(4)=ZLIMS(JS,IMAX)
C
C
      CALL ZYLND2 (ZX0,ZY0,ZZ0,ZVX,ZVY,ZVZ,TAR,RZYL,NZAD+2,NINNE,NIN,
     .     ILCOL(JS),IGFIL(JS).NE.0,JS,4,AL,4,AR,PHIAN,PHIEN,20,IPAR)
C
      RETURN
      END
C
C
      SUBROUTINE SECANG (B0,B1,B2,RAD,ANG1,ANG2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'CCONA'
C
      ANG1=0.
      ANG2=360.
C
      IF (ABS(B1).LT.EPS10.AND.ABS(B2).LT.EPS10) RETURN
C
      IF (ABS(B1).LT.EPS10) THEN
        Y1=-B0/B2
        IF (ABS(Y1).GE.RAD) RETURN
        X1=SQRT(RAD*RAD-Y1*Y1)
        X2=-X1
        Y2=Y1
      ELSEIF (ABS(B2).LT.EPS10) THEN
        X1=-B0/B1
        IF (ABS(X1).GE.RAD) RETURN
        Y1=SQRT(RAD*RAD-X1*X1)
        X2=X1
        Y2=-Y1
      ELSE
        B1Q=B1*B1
        B2Q=B2*B2
        B0Q=B0*B0
        XNEN=B1Q+B2Q
        B0B1=B0*B1
        RD=B0B1*B0B1/XNEN/XNEN-(B0Q-B2Q*RAD*RAD)/XNEN
        IF (RD.LT.0.) RETURN
        X1=-B0B1/XNEN+SQRT(RD)
        X2=-B0B1/XNEN-SQRT(RD)
        Y1=-B0/B2-B1/B2*X1
        Y2=-B0/B2-B1/B2*X2
      ENDIF
C
      ANG1=ATAN2(Y1,X1)*RADDEG
      IF (ANG1.LT.0.) ANG1=ANG1+360.
      ANG2=ATAN2(Y2,X2)*RADDEG
      IF (ANG2.LT.0.) ANG2=ANG2+360.
C
      RETURN
      END
C
      SUBROUTINE ZYLND2(X0,Y0,Z0,VX,VY,VZ,TAR,RAD,NK,NP,NA,IO,NF,NUM,
     .                   ILEFT,AL,IRIGHT,AR,PHIAN,PHIEN,NDP,IPART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  ZYLINDERACHSE IST GERADE X+T*V, T1<T<T2, RADIUS RAD.
C  NK KREISE, NA STUETZSTELLEN AUF KREIS (POLYGON, NA-ECK)
C  NP ANZAHL DER LINIEN FUER PHI=CONST.
C  WENN ILEFT (IRIGHT) .NE. 0, DANN WIRD DER ERSTE (I=1, D.H.T=T1)
C  BZW DER LETZTE (I=NK, T=T2) KREIS DES ZYLINDERS DURCH DIE SCHNITT-
C  FLAECHE DIESES ZYLINDERS MIT DER GLEICHUNG
C  A(1)+A(2)*X+A(3)*Y+....+A(10)*Y*Z=0.
C  ERSETZT. T1 UND T2 SIND SO EINZUGEBEN, DASS DIE SCHNITTFLAECHEN
C  AUSSERHALB DIESES PARAMETERBEREICHES LIEGEN.
C  DABEI SIND NUR DIE ERSTEN ILEFT (IRIGHT) KOEFFIZIENTEN .NE.0
C  D.H. ILEFT (IRIGHT) <= 4 ENTSPRICHT DEM SCHNITT MIT EINER EBENE.
C
C
      INCLUDE 'CCONA'
C
      DIMENSION P(3,101),XP(101),YP(101),AL(*),AR(*),TAR(NK),
     .          IPART(NK),PHIAN(NDP,*),PHIEN(NDP,*)
      REAL*4 XPS(101),YPS(101)
      LOGICAL NF
C
      NA=MIN0(NA,100)
      NP=MIN0(NP,100)
      PI=4.*ATAN(1.)
C
      AX=VX
      AY=VY
      AZ=VZ
      IF (ABS(AZ).GT.1.E-20) THEN
        BX=1.
        BY=1.
        BZ=-(AX+AY)/AZ
        PHID=ACOS(BX/SQRT(2.+BZ*BZ))
      ELSE IF (ABS(AX).GT.1.E-20) THEN
        BY=1.
        BZ=1.
        BX=-(AY+AZ)/AX
        PHID=ACOS(BY/SQRT(2.+BX*BX))
      ELSE IF (ABS(AY).GT.1.E-20) THEN
        BX=1.
        BZ=1.
        BY=-(AX+AZ)/AY
        PHID=ACOS(BZ/SQRT(2.+BY*BY))
      ELSE
        WRITE (6,*) 'FEHLER IN DER EINGABE VON VX,VY,VZ ',NUM,VX,VY,VZ
        CALL EXIT
      ENDIF
      WRITE (6,*) ' PHID ',PHID*RADDEG
      WRITE (6,*) ' VX,VY,VZ ',VX,VY,VZ
      CX=AY*BZ-AZ*BY
      CY=AZ*BX-AX*BZ
      CZ=AX*BY-AY*BX
      BA=SQRT(AX*AX+AY*AY+AZ*AZ)
      BB=SQRT(BX*BX+BY*BY+BZ*BZ)
      BC=SQRT(CX*CX+CY*CY+CZ*CZ)
      AX=AX/BA
      AY=AY/BA
      AZ=AZ/BA
      BX=BX/BB
      BY=BY/BB
      BZ=BZ/BB
      CX=CX/BC
      CY=CY/BC
      CZ=CZ/BC
C
C PLOTTE DIE KREISSTUECKE, NK STUECK
      DO 2 I=1,NK
        T=TAR(I)
        DO 10 IP=1,IPART(I)
          DANG=(PHIEN(I,IP)-PHIAN(I,IP))/DBLE(NA)*PI/180.
          PHAN=PHIAN(I,IP)*PI/180.-PHID
C  BERECHNE EINEN KREIS IN RICHTUNG DER ZYLINDERACHSE, MIT
C  DEM 0-PUNKT ALS MITTELPUNKT UND RADIUS RAD
C  VON PHIAN BIS PHIEN
          DO 1 J=1,NA+1
            PHI=PHAN+(J-1)*DANG
            XK=RAD*COS(PHI)
            YK=RAD*SIN(PHI)
            P(1,J)=XK*BX+YK*CX
            P(2,J)=XK*BY+YK*CY
            P(3,J)=XK*BZ+YK*CZ
    1     CONTINUE
C  SETZTE EINEN VERSCHIEBUNGSVEKTOR AUF DER ZYLINDERACHSE
C  INNERHALB DES BEREICHES T1----T2, FUER SHNITT-OPTION
          PXS=X0+(TAR(1)+TAR(NK))/2.*VX
          PYS=Y0+(TAR(1)+TAR(NK))/2.*VY
          PZS=Z0+(TAR(1)+TAR(NK))/2.*VZ
          PX=X0+T*VX
          PY=Y0+T*VY
          PZ=Z0+T*VZ
          IF (I.EQ.1.AND.ILEFT.NE.0) THEN
          WRITE (6,*) ' LEFT END OF ZYLINDER '
          WRITE (6,*) (AL(ILFT),ILFT=1,ILEFT)
          CALL SHNITT(P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,1,NA+1,1)
          ELSEIF (I.EQ.NK.AND.IRIGHT.NE.0) THEN
          WRITE (6,*) ' RIGHT END OF ZYLINDER '
          WRITE (6,*) (AR(IRGHT),IRGHT=1,IRIGHT)
          CALL SHNITT(P,PXS,PYS,PZS,VX,VY,VZ,AR,IRIGHT,XP,YP,1,NA+1,1)
          ELSE
            DO 3 J=1,NA+1
              PXX=P(1,J)+PX
              PYY=P(2,J)+PY
              PZZ=P(3,J)+PZ
3             CALL PL3D (PXX,PYY,PZZ,XP(J),YP(J))
          ENDIF
          IF (IO.GE.2) CALL GRNWPN(IO)
          do 7 jj=1,na+1
            xps(jj)=xp(jj)
            yps(jj)=yp(jj)
7         continue
          CALL GRLN (XPS,YPS,NA+1)
C  FAERBE DIE ENDEN DES ZYLINDERS EIN
          IF ((I.EQ.1.OR.I.EQ.NK).AND.NF) CALL GRFILL(NA+1,XPS,YPS,1,1)
          IF (IO.GE.2) CALL GRNWPN(1)
10      CONTINUE
2     CONTINUE
C
C  PLOTTE PHI=CONST LINIEN, INSGESAMT NP STUECK
C
C  SETZE NEUEN KREIS UM 0-PUNKT, MIT NP STUETZSTELLEN
      DANG=2.*PI/DBLE(NP-1)
      PHAN=-PHID
      DO 6 J=1,NP
        PHI=PHAN+(J-1)*DANG
        XK=RAD*COS(PHI)
        YK=RAD*SIN(PHI)
        P(1,J)=XK*BX+YK*CX
        P(2,J)=XK*BY+YK*CY
        P(3,J)=XK*BZ+YK*CZ
6     CONTINUE
C
      IA=1
      IE=NK
      IF (ILEFT.NE.0) IA=2
      IF (IRIGHT.NE.0) IE=NK-1
      DO 5 J=1,NP
        PHI=PHAN+(J-1)*DANG
        PHIDEG=PHI*RADDEG
        JP=0
        IF (ILEFT.NE.0) THEN
          DO 11 IP=1,IPART(1)
            IF (PHIDEG.LT.PHIAN(1,IP).OR.PHIDEG.GT.PHIEN(1,IP)) GOTO 11
            JP=JP+1
            CALL SHNITT(P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,J,J,JP)
11        CONTINUE
        ENDIF
        DO 4 I=IA,IE
          DO 12 IP=1,IPART(I)
            IF (PHIDEG.LT.PHIAN(I,IP).OR.PHIDEG.GT.PHIEN(I,IP)) GOTO 12
            JP=JP+1
            T=TAR(I)
            PX=X0+T*VX
            PY=Y0+T*VY
            PZ=Z0+T*VZ
            CALL PL3D (P(1,J)+PX,P(2,J)+PY,P(3,J)+PZ,XP(JP),YP(JP))
            GOTO 4
12        CONTINUE
          GOTO 14
4       CONTINUE
        IF (IRIGHT.NE.0) THEN
          DO 13 IP=1,IPART(NK)
            IF (PHIDEG.LT.PHIAN(NK,IP).OR.PHIDEG.GT.PHIEN(NK,IP))
     .          GOTO 13
            JP=JP+1
            CALL SHNITT (P,PXS,PYS,PZS,VX,VY,VZ,AR,IRIGHT,XP,YP,J,J,JP)
13        CONTINUE
        ENDIF
14      CONTINUE
        IF (JP.GT.1) THEN
          do 9 jj=1,jp
            xps(jj)=xp(jj)
            yps(jj)=yp(jj)
9         continue
          CALL GRLN (XPS,YPS,JP)
        endif
5     CONTINUE
      RETURN
      END

      SUBROUTINE ISOLNE (AORIG,IBLD,ICURV,
     .                   IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE CARRIES OUT A CONTOUR PLOT
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ON VERTICES IS PERFORMED
C
C  IN CASE (LEVGEO=1 OR LEVGEO=2).AND.NLTOR, THE RHOSRF AND ZSURF GRIDS
C                                            ARE USED.
C  IN CASE (LEVGEO=1            ).AND.NLPOL, THE RHOSRF AND PSURF GRIDS
C                                            ARE USED.
C  IN CASE (LEVGEO=3 OR LEVGEO=2).AND.NLPOL, THE FULL POLYGON GRIDS
C                                            ARE USED.
C  IN CASE LEVGEO=4,                         THE FULL TRIANGULAR MESH
C                                            IS USED
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CPOLYG'
      INCLUDE 'CGRID'
      INCLUDE 'CTRIG'
      INCLUDE 'CGEOM'
      INCLUDE 'CLOGAU'
      INCLUDE 'CPLOT'
C
      DIMENSION AORIG(*)
      DIMENSION A(N1ST,N2ND+N3RD),AA(NRAD,1)
      EQUIVALENCE (A(1,1),AA(1,1))
      DIMENSION XX(*),YY(*)
      REAL*4 XY(8000)
      REAL*4 YH
C
      LOGICAL LOGL,TRC
      CHARACTER*72 TEXT1,HEAD,RUNID,TXHEAD
      CHARACTER*24 TEXT2,TEXT3
      CHARACTER*17 CH
C
      ZINT(X1,X2,AA1,AA2,A1)=X1+(A1-AA1)/(AA2-AA1+1.D-30)*(X2-X1)
C
C     PLOT 18 CONTOURS, WITH 6 DIFFERENT COLOURS
      NISO=18
      IISO=3
C
C  SEARCH FOR MAXIMUM AND MINIMUM RMI AND RMA
C
      RMI=1.D60
      RMA=-1.D60
C
      IF (LEVGEO .LE. 2) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 21 IR=1,IXX-1
          DO 21 IP=1,IYY-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
21      CONTINUE
      ELSEIF (LEVGEO.EQ.3) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 20 IR=1,NR1ST-1
        DO 20 IPART=1,NPPLG
          DO 20 IP=NPOINT(1,IPART),NPOINT(2,IPART)-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
20      CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO 22 I=1,NTRII
          RMI=MIN(RMI,AORIG(I))
          RMA=MAX(RMA,AORIG(I))
22      CONTINUE
      ENDIF
C
C  INTERPOLATE: AORIG --> A
C
      IF (LEVGEO.EQ.4) THEN
        CALL CELINT(AORIG,AA,LOGL,IBLD,ICURV,NRAD)
      ELSEIF (LEVGEO.LE.3) THEN
        CALL CELINT(AORIG,A,LOGL,IBLD,ICURV,N1ST)
      ENDIF
C
C  SEARCH FOR XMIN,XMAX,YMIN,YMAX
C
      XMIN=1.D60
      YMIN=1.D60
      XMAX=-1.D60
      YMAX=-1.D60
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = ZSURF(1)
        YMAX = ZSURF(NT3RD)
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = PSURF(1)
        YMAX = PSURF(NP2ND)
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
C  SUFFICIENT TO SEARCH ON OUTERMOST RADIAL SURFACE (BECAUSE: CONVEX)
        DO 5 IP = 1,NP2ND
          XMIN = MIN(XMIN,XPOL(NR1ST,IP))
          XMAX = MAX(XMAX,XPOL(NR1ST,IP))
          YMIN = MIN(YMIN,YPOL(NR1ST,IP))
          YMAX = MAX(YMAX,YPOL(NR1ST,IP))
5       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
C  SEARCH ON WHOLE MESH
        DO 1 IR=1,NR1ST
          NP1=NPOINT(1,1)
          NP2=NPOINT(2,NPPLG)
          XMIN=MIN(XMIN,XPOL(IR,NP1),XPOL(IR,NP2))
          YMIN=MIN(YMIN,YPOL(IR,NP1),YPOL(IR,NP2))
          XMAX=MAX(XMAX,XPOL(IR,NP1),XPOL(IR,NP2))
          YMAX=MAX(YMAX,YPOL(IR,NP1),YPOL(IR,NP2))
1       CONTINUE
C
        DO 4 IPART=1,NPPLG
          DO 2 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            XMIN=MIN(XMIN,XPOL(1,IP))
            YMIN=MIN(YMIN,YPOL(1,IP))
            XMAX=MAX(XMAX,XPOL(1,IP))
            YMAX=MAX(YMAX,YPOL(1,IP))
2         CONTINUE
          DO 3 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            XMIN=MIN(XMIN,XPOL(NR1ST,IP))
            YMIN=MIN(YMIN,YPOL(NR1ST,IP))
            XMAX=MAX(XMAX,XPOL(NR1ST,IP))
            YMAX=MAX(YMAX,YPOL(NR1ST,IP))
3         CONTINUE
4       CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
C  SEARCH ON WHOLE MESH
        DO I=1,NRKNOT
          XMIN=MIN(XMIN,XTRIAN(I))
          YMIN=MIN(YMIN,YTRIAN(I))
          XMAX=MAX(XMAX,XTRIAN(I))
          YMAX=MAX(YMAX,YTRIAN(I))
        ENDDO
      ENDIF
C
C
      CM=20.
      DX=(XMAX-XMIN)*FCABS1(IBLD)
      DY=(YMAX-YMIN)*FCABS2(IBLD)
      FAK=CM/MAX(DX,DY)
C
C  PLOT FRAME
C
      CALL GRNXTB (1)
      CALL GRSCLC (10.,4.,sngl(10.+DX*FAK),sngl(4.+DY*FAK))
      CALL GRSCLV (sngl(XMIN),sngl(YMIN),sngl(XMAX),sngl(YMAX))
      CALL GRAXS (7,'X=3,Y=3',6,'R (CM)',6,'Z (CM)')
C
C  SCALE FACTORS: USER CO-ORDINATES TO CM:
C  X-DIRECTION:
      SCLFCX=((10.+DX*FAK)-10.)/(XMAX-XMIN)
C  Y-DIRECTION:
      SCLFCY=((4.+DY*FAK)-4.)/(YMAX-YMIN)
C
C  PLOT BOUNDARY OF MESH
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR) THEN
        CALL GRJMP(sngl(XMIN),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMIN))
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        CALL GRJMP(sngl(XMIN),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMIN))
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
        DO 7 IR=1,NR1ST,NR1STM
          CALL GRJMP(sngl(XPOL(IR,1)),sngl(YPOL(IR,1)))
          DO 9 IP = 2,NP2ND
9           CALL GRDRW(sngl(XPOL(IR,IP)),sngl(YPOL(IR,IP)))
7       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
        CALL GRJMP(sngl(XPOL(1,NPOINT(1,1))),sngl(YPOL(1,NPOINT(1,1))))
          DO 10 IR=2,NR1ST
10          CALL GRDRW (sngl(XPOL(IR,NPOINT(1,1))),
     .                  sngl(YPOL(IR,NPOINT(1,1))))
C
        CALL GRJMP (sngl(XPOL(1,NPOINT(2,NPPLG))),
     .              sngl(YPOL(1,NPOINT(2,NPPLG))))
        DO 11 IR=2,NR1ST
          NP=NPOINT(2,NPPLG)
11        CALL GRDRW (sngl(XPOL(IR,NP)),sngl(YPOL(IR,NP)))
        DO 15 I=1,NPPLG
          CALL GRJMP (sngl(XPOL(1,NPOINT(1,I))),
     .                sngl(YPOL(1,NPOINT(1,I))))
          DO 12 IP=NPOINT(1,I),NPOINT(2,I)
12          CALL GRDRW (sngl(XPOL(1,IP)),sngl(YPOL(1,IP)))
          CALL GRJMP (sngl(XPOL(NR1ST,NPOINT(1,I))),
     .                sngl(YPOL(NR1ST,NPOINT(1,I))))
          DO 13 IP=NPOINT(1,I),NPOINT(2,I)
13          CALL GRDRW (sngl(XPOL(NR1ST,IP)),sngl(YPOL(NR1ST,IP)))
15      CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO ITR=1,NTRII
          IF (NCHBAR(1,ITR) .EQ. 0) THEN
            CALL GRJMP(SNGL(XTRIAN(NECKE(1,ITR))),
     .                 SNGL(YTRIAN(NECKE(1,ITR))))
            CALL GRDRW(SNGL(XTRIAN(NECKE(2,ITR))),
     .                 SNGL(YTRIAN(NECKE(2,ITR))))
          ENDIF
          IF (NCHBAR(2,ITR) .EQ. 0) THEN
            CALL GRJMP(SNGL(XTRIAN(NECKE(3,ITR))),
     .                 SNGL(YTRIAN(NECKE(3,ITR))))
            CALL GRDRW(SNGL(XTRIAN(NECKE(2,ITR))),
     .                 SNGL(YTRIAN(NECKE(2,ITR))))
          ENDIF
          IF (NCHBAR(3,ITR) .EQ. 0) THEN
            CALL GRJMP(SNGL(XTRIAN(NECKE(1,ITR))),
     .                 SNGL(YTRIAN(NECKE(1,ITR))))
            CALL GRDRW(SNGL(XTRIAN(NECKE(3,ITR))),
     .                 SNGL(YTRIAN(NECKE(3,ITR))))
          ENDIF
        ENDDO
      ENDIF
C
      IF (LOGL) THEN
        IF (ZMI.EQ.666.) THEN
          RAMIN=LOG10(MAX(1.D-48,RMI))
        ELSE
          RAMIN=LOG10(MAX(1.D-48,ZMI))
        ENDIF
        IF (ZMA.EQ.666.) THEN
          RAMAX=LOG10(MAX(1.D-48,RMA))
        ELSE
          RAMAX=LOG10(MAX(1.D-48,ZMA))
        ENDIF
      ELSE
        RAMIN=ZMI
        IF (ZMI.EQ.666.) RAMIN=RMI
        RAMAX=ZMA
        IF (ZMA.EQ.666.) RAMAX=RMA
      ENDIF
C
      IF (ABS(RAMAX-RAMIN)/MAX(RAMAX,1.D-30).LT.0.01)
     .    RAMAX=RAMIN+0.01*RAMIN*SIGN(1.D0,RAMIN)
C
C
      DA=(RAMAX-RAMIN)/DBLE(NISO-1)
C
      ICOLOR=1
      IF (((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR).OR.
     .     (LEVGEO.EQ.1.AND.NLPOL)) THEN
        DO 1000 IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO 1100 IR=1,IXX-1
          DO 1100 IP=1,IYY-1
            A1=A(IR,IP)
            A2=A(IR+1,IP)
            A3=A(IR+1,IP+1)
            A4=A(IR,IP+1)
            ACMIN=MIN(A1,A2,A3,A4)
            ACMAX=MAX(A1,A2,A3,A4)
            IT=0
            IF (ACONT.GE.ACMIN.AND.ACONT.LE.ACMAX) THEN
              IF (ACONT.GE.MIN(A1,A2).AND.
     .            ACONT.LE.MAX(A1,A2)) THEN
                XY(IC+1)=ZINT(XX(IR),XX(IR+1),A1,A2,ACONT)
                XY(IC+2)=YY(IP)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A2,A3).AND.
     .            ACONT.LE.MAX(A2,A3)) THEN
                XY(IC+1)=XX(IR+1)
                XY(IC+2)=ZINT(YY(IP),YY(IP+1),A2,A3,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A3,A4).AND.
     .            ACONT.LE.MAX(A3,A4)) THEN
                XY(IC+1)=ZINT(XX(IR+1),XX(IR),A3,A4,ACONT)
                XY(IC+2)=YY(IP+1)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A4,A1).AND.
     .            ACONT.LE.MAX(A4,A1)) THEN
                XY(IC+1)=XX(IR)
                XY(IC+2)=ZINT(YY(IP+1),YY(IP),A4,A1,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (MOD(IT,2).EQ.1) THEN
C               WRITE (6,*) ' ACHTUNG !!!! '
C               WRITE (6,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (6,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
                XY(IC+1)=XY(IC-IT*2+1)
                XY(IC+2)=XY(IC-IT*2+2)
                IC=IC+2
              ENDIF
              IF (IC+8.GT.800) THEN
                CALL XYPLOT (XY,IC)
                IC=0
              ENDIF
            ENDIF
1100      CONTINUE
          IF (IC.GT.0) CALL XYPLOT (XY,IC)
1000    CONTINUE
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL) THEN
        DO 100 IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO 110 IR=1,NR1STM
          DO 110 IPART=1,NPPLG
          DO 110 IP=NPOINT(1,IPART),NPOINT(2,IPART)-1
            A1=A(IR,IP)
            A2=A(IR+1,IP)
            A3=A(IR+1,IP+1)
            A4=A(IR,IP+1)
            ACMIN=MIN(A1,A2,A3,A4)
            ACMAX=MAX(A1,A2,A3,A4)
            IT=0
            IF (ACONT.GE.ACMIN.AND.ACONT.LE.ACMAX) THEN
              IF (ACONT.GE.MIN(A1,A2).AND.
     .            ACONT.LE.MAX(A1,A2)) THEN
                XY(IC+1)=ZINT(XPOL(IR,IP),XPOL(IR+1,IP),A1,A2,ACONT)
                XY(IC+2)=ZINT(YPOL(IR,IP),YPOL(IR+1,IP),A1,A2,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A2,A3).AND.
     .            ACONT.LE.MAX(A2,A3)) THEN
                XY(IC+1)=ZINT(XPOL(IR+1,IP),XPOL(IR+1,IP+1),A2,A3,ACONT)
                XY(IC+2)=ZINT(YPOL(IR+1,IP),YPOL(IR+1,IP+1),A2,A3,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A3,A4).AND.
     .            ACONT.LE.MAX(A3,A4)) THEN
                XY(IC+1)=ZINT(XPOL(IR+1,IP+1),XPOL(IR,IP+1),A3,A4,ACONT)
                XY(IC+2)=ZINT(YPOL(IR+1,IP+1),YPOL(IR,IP+1),A3,A4,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A4,A1).AND.
     .            ACONT.LE.MAX(A4,A1)) THEN
                XY(IC+1)=ZINT(XPOL(IR,IP+1),XPOL(IR,IP),A4,A1,ACONT)
                XY(IC+2)=ZINT(YPOL(IR,IP+1),YPOL(IR,IP),A4,A1,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (MOD(IT,2).EQ.1) THEN
C               WRITE (6,*) ' ACHTUNG !!!! '
C               WRITE (6,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (6,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
                XY(IC+1)=XY(IC-IT*2+1)
                XY(IC+2)=XY(IC-IT*2+2)
                IC=IC+2
              ENDIF
              IF (IC+8.GT.800) THEN
                CALL XYPLOT (XY,IC)
                IC=0
              ENDIF
            ENDIF
110       CONTINUE
          IF (IC.GT.0) CALL XYPLOT (XY,IC)
100     CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO  ITR=1,NTRII
            A1=AA(NECKE(1,ITR),1)
            A2=AA(NECKE(2,ITR),1)
            A3=AA(NECKE(3,ITR),1)
            ACMIN=MIN(A1,A2,A3)
            ACMAX=MAX(A1,A2,A3)
            IT=0
            IF (ACONT.GE.ACMIN.AND.ACONT.LE.ACMAX) THEN
              IF (ACONT.GE.MIN(A1,A2).AND.
     .            ACONT.LE.MAX(A1,A2)) THEN
                XY(IC+1)=ZINT(XTRIAN(NECKE(1,ITR)),
     .                        XTRIAN(NECKE(2,ITR)),A1,A2,ACONT)
                XY(IC+2)=ZINT(YTRIAN(NECKE(1,ITR)),
     .                        YTRIAN(NECKE(2,ITR)),A1,A2,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A2,A3).AND.
     .            ACONT.LE.MAX(A2,A3)) THEN
                XY(IC+1)=ZINT(XTRIAN(NECKE(2,ITR)),
     .                        XTRIAN(NECKE(3,ITR)),A2,A3,ACONT)
                XY(IC+2)=ZINT(YTRIAN(NECKE(2,ITR)),
     .                        YTRIAN(NECKE(3,ITR)),A2,A3,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A3,A1).AND.
     .            ACONT.LE.MAX(A3,A1)) THEN
                XY(IC+1)=ZINT(XTRIAN(NECKE(3,ITR)),
     .                        XTRIAN(NECKE(1,ITR)),A3,A1,ACONT)
                XY(IC+2)=ZINT(YTRIAN(NECKE(3,ITR)),
     .                        YTRIAN(NECKE(1,ITR)),A3,A1,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (MOD(IT,2).EQ.1) THEN
C               WRITE (6,*) ' ACHTUNG !!!! '
C               WRITE (6,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (6,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
                XY(IC+1)=XY(IC-IT*2+1)
                XY(IC+2)=XY(IC-IT*2+2)
                IC=IC+2
              ENDIF
              IF (IC+8.GT.800) THEN
                CALL XYPLOT (XY,IC)
                IC=0
              ENDIF
            ENDIF
          ENDDO
          IF (IC.GT.0) CALL XYPLOT (XY,IC)
        ENDDO
      ENDIF
C
C
C     WRITE TEXT AND MEAN VALUE ONTO THE PLOT
C
      CALL GRNWPN (1)
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,YH,72,RUNID)
      YH=26.75
      CALL GRTXT (1.,YH,72,HEAD)
      YH=26.00
      CALL GRTXT (1.,YH,72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,YH,10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,YH-0.5,10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,YH-1.,10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,YH-2.,10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') RMA
      CALL GRTXT (1.,YH-2.5,10,CH)
      CALL GRTXT (1.,YH-3.,10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') RMI
      CALL GRTXT (1.,YH-3.5,10,CH)
C
      ICOLOR=1
      YH=YH-4.
      DO 200 IS=1,NISO
        ACONT=RAMIN+(IS-1)*DA
        IF (LOGL) ACONT=10.**ACONT
        IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
        CALL GRNWPN(ICOLOR)
        YH=YH-0.5
        CALL GRJMP (1.,YH+0.25)
        CALL GRDRW (2.5,YH+0.25)
        CALL GRNWPN (1)
        WRITE (CH,'(1P,E10.3)') ACONT
        CALL GRTXT (3.,YH,10,CH)
200   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE VECLNE (AORIG,BORIG,IBLD,ICURV,IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A VECTOR FIELD PLOT
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ON VERTICES IS PERFORMED
C
C  IN CASE (LEVGEO=1 OR LEVGEO=2).AND.NLTOR, THE RHOSRF AND ZSURF GRIDS
C                                            ARE USED.
C  IN CASE (LEVGEO=1            ).AND.NLPOL, THE RHOSRF AND PSURF GRIDS
C                                            ARE USED.
C  IN CASE (LEVGEO=3 OR LEVGEO=2).AND.NLPOL, THE FULL POLYGON GRIDS
C                                            ARE USED.
C  IN CASE LEVGEO=4,                         THE FULL TRIANGULAR MESH
C                                            IS USED
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CPOLYG'
      INCLUDE 'CGRID'
      INCLUDE 'CTRIG'
      INCLUDE 'CGEOM'
      INCLUDE 'CLOGAU'
      INCLUDE 'CESTIM'
      INCLUDE 'CPLOT'
C
      DIMENSION XCOMTR(NTRIS),YCOMTR(NTRIS)
      EQUIVALENCE (XCOM(1,1),XCOMTR(1)),
     .            (YCOM(1,1),YCOMTR(1))
      DIMENSION AORIG(*),BORIG(*)
      DIMENSION XX(*),YY(*)
      REAL*4 XY(800)
      REAL*4 YH
C
      LOGICAL LOGL,TRC
      CHARACTER*72 TEXT1,HEAD,RUNID,TXHEAD
      CHARACTER*24 TEXT2,TEXT3
      CHARACTER*17 CH
C
C  SEARCH FOR XMIN,XMAX,YMIN,YMAX
C
      XMIN=1.D60
      YMIN=1.D60
      XMAX=-1.D60
      YMAX=-1.D60
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = ZSURF(1)
        YMAX = ZSURF(NT3RD)
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = PSURF(1)
        YMAX = PSURF(NP2ND)
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
C  SUFFICIENT TO SEARCH ON OUTERMOST RADIAL SURFACE (BECAUSE: CONVEX)
        DO 5 IP = 1,NP2ND
          XMIN = MIN(XMIN,XPOL(NR1ST,IP))
          XMAX = MAX(XMAX,XPOL(NR1ST,IP))
          YMIN = MIN(YMIN,YPOL(NR1ST,IP))
          YMAX = MAX(YMAX,YPOL(NR1ST,IP))
5       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
C  SEARCH ON WHOLE MESH
        DO 1 IR=1,NR1ST
          NP1=NPOINT(1,1)
          NP2=NPOINT(2,NPPLG)
          XMIN=MIN(XMIN,XPOL(IR,NP1),XPOL(IR,NP2))
          YMIN=MIN(YMIN,YPOL(IR,NP1),YPOL(IR,NP2))
          XMAX=MAX(XMAX,XPOL(IR,NP1),XPOL(IR,NP2))
          YMAX=MAX(YMAX,YPOL(IR,NP1),YPOL(IR,NP2))
1       CONTINUE
C
        DO 4 IPART=1,NPPLG
          DO 2 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            XMIN=MIN(XMIN,XPOL(1,IP))
            YMIN=MIN(YMIN,YPOL(1,IP))
            XMAX=MAX(XMAX,XPOL(1,IP))
            YMAX=MAX(YMAX,YPOL(1,IP))
2         CONTINUE
          DO 3 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            XMIN=MIN(XMIN,XPOL(NR1ST,IP))
            YMIN=MIN(YMIN,YPOL(NR1ST,IP))
            XMAX=MAX(XMAX,XPOL(NR1ST,IP))
            YMAX=MAX(YMAX,YPOL(NR1ST,IP))
3         CONTINUE
4       CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
C  SEARCH ON WHOLE MESH
        DO I=1,NRKNOT
          XMIN=MIN(XMIN,XTRIAN(I))
          YMIN=MIN(YMIN,YTRIAN(I))
          XMAX=MAX(XMAX,XTRIAN(I))
          YMAX=MAX(YMAX,YTRIAN(I))
        ENDDO
      ENDIF
C
C
      CM=20.
      DX=FCABS1(IBLD) * (XMAX-XMIN)
      DY=FCABS2(IBLD) * (YMAX-YMIN)
      FAK=CM/MAX(DX,DY)
C
C  PLOT FRAME
C
      CALL GRNXTB (1)
      CALL GRSCLC (10.,4.,sngl(10.+DX*FAK),sngl(4.+DY*FAK))
      CALL GRSCLV (sngl(XMIN),sngl(YMIN),sngl(XMAX),sngl(YMAX))
      CALL GRAXS (7,'X=3,Y=3',6,'R (CM)',6,'Z (CM)')
C
C  SCALE FACTORS: USER CO-ORDINATES TO CM:
C  X-DIRECTION:
      SCLFCX=((10.+DX*FAK)-10.)/(XMAX-XMIN)
C  Y-DIRECTION:
      SCLFCY=((4.+DY*FAK)-4.)/(YMAX-YMIN)
C
C  PLOT BOUNDARY OF MESH
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR) THEN
        CALL GRJMP(sngl(XMIN),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMIN))
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        CALL GRJMP(sngl(XMIN),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMAX))
        CALL GRDRW(sngl(XMAX),sngl(YMIN))
        CALL GRDRW(sngl(XMIN),sngl(YMIN))
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
        DO 7 IR=1,NR1ST,NR1STM
          CALL GRJMP(sngl(XPOL(IR,1)),sngl(YPOL(IR,1)))
          DO 9 IP = 2,NP2ND
9           CALL GRDRW(sngl(XPOL(IR,IP)),sngl(YPOL(IR,IP)))
7       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
        CALL GRJMP(sngl(XPOL(1,NPOINT(1,1))),sngl(YPOL(1,NPOINT(1,1))))
          DO 10 IR=2,NR1ST
10          CALL GRDRW (sngl(XPOL(IR,NPOINT(1,1))),
     .                  sngl(YPOL(IR,NPOINT(1,1))))
C
        CALL GRJMP (sngl(XPOL(1,NPOINT(2,NPPLG))),
     .              sngl(YPOL(1,NPOINT(2,NPPLG))))
        DO 11 IR=2,NR1ST
          NP=NPOINT(2,NPPLG)
11        CALL GRDRW (sngl(XPOL(IR,NP)),sngl(YPOL(IR,NP)))
        DO 15 I=1,NPPLG
          CALL GRJMP (sngl(XPOL(1,NPOINT(1,I))),
     .                sngl(YPOL(1,NPOINT(1,I))))
          DO 12 IP=NPOINT(1,I),NPOINT(2,I)
12          CALL GRDRW (sngl(XPOL(1,IP)),sngl(YPOL(1,IP)))
          CALL GRJMP (sngl(XPOL(NR1ST,NPOINT(1,I))),
     .                sngl(YPOL(NR1ST,NPOINT(1,I))))
          DO 13 IP=NPOINT(1,I),NPOINT(2,I)
13          CALL GRDRW (sngl(XPOL(NR1ST,IP)),sngl(YPOL(NR1ST,IP)))
15      CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO ITR=1,NTRII
          IF (NCHBAR(1,ITR) .EQ. 0) THEN
            CALL GRJMP(SNGL(XTRIAN(NECKE(1,ITR))),
     .                 SNGL(YTRIAN(NECKE(1,ITR))))
            CALL GRDRW(SNGL(XTRIAN(NECKE(2,ITR))),
     .                 SNGL(YTRIAN(NECKE(2,ITR))))
          ENDIF
          IF (NCHBAR(2,ITR) .EQ. 0) THEN
            CALL GRJMP(SNGL(XTRIAN(NECKE(3,ITR))),
     .                 SNGL(YTRIAN(NECKE(3,ITR))))
            CALL GRDRW(SNGL(XTRIAN(NECKE(2,ITR))),
     .                 SNGL(YTRIAN(NECKE(2,ITR))))
          ENDIF
          IF (NCHBAR(3,ITR) .EQ. 0) THEN
            CALL GRJMP(SNGL(XTRIAN(NECKE(1,ITR))),
     .                 SNGL(YTRIAN(NECKE(1,ITR))))
            CALL GRDRW(SNGL(XTRIAN(NECKE(3,ITR))),
     .                 SNGL(YTRIAN(NECKE(3,ITR))))
          ENDIF
        ENDDO
      ENDIF
C
C PLOT 2D VECTOR FIELD
C
      IF (((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR).OR.
     .     (LEVGEO.EQ.1.AND.NLPOL)) THEN
C
C 1. SEARCH FOR MINIMA AND MAXIMA
          VMIN=1.D60
          VMAX=-1.D60
          DO 900 IR=1,IXX-1
          DO 900 IP=1,IYY-1
            IRAD = IR + (IP-1)*NR1ST
            VABS = SQRT(AORIG(IRAD)**2 + BORIG(IRAD)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
900       CONTINUE
C 2. SCALING
          DO 1100 IR=1,IXX-1
          DO 1100 IP=1,IYY-1
            IRAD = IR + (IP-1)*NR1ST
            VX =  5 * (AORIG(IRAD) / VMAX)
            VY =  5 * (BORIG(IRAD) / VMAX)
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOM(IR,IP)-VX/2
            YM=YCOM(IR,IP)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(sngl(xm),sngl(ym),sngl(xm+vx),sngl(ym+vy),
     .                  SNGL(PLFL),SNGL(BRFL),1)
1100      CONTINUE
C
      ELSEIF ((LEVGEO.EQ.2.AND.NLPOL).OR.
     .         LEVGEO.EQ.3) THEN
C
        WRITE (6,*) 'VECTOR PLOTS FOR (LEVGEO.EQ.2.AND.NLPOL).OR.'
        WRITE (6,*) 'LEVGEO.EQ.3 NOT READY'
        CALL EXIT
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C PLOT 2D VECTOR FIELD
C
C 1. SEARCH FOR MINIMA AND MAXIMA
          VMIN=1.D60
          VMAX=-1.D60
          DO 1500 IR=1,NTRII
            VABS = SQRT(AORIG(IR)**2 + BORIG(IR)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
1500      CONTINUE
C 2. SCALING
          DO 1700 IR=1,NTRII
            VX = 10 * (AORIG(IR) / VMAX)
            VY = 10 * (BORIG(IR) / VMAX)
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOMTR(IR)-VX/2
            YM=YCOMTR(IR)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(sngl(xm),sngl(ym),sngl(xm+vx),sngl(ym+vy),
     .                  SNGL(PLFL),SNGL(BRFL),1)
1700      CONTINUE
      ENDIF
C
2000  CONTINUE
C
C     WRITE TEXT, MAXIMUM AND MINIMUM VALUE ONTO THE PLOT
C
      TEXT1='2D VECTOR FIELD'
      CALL GRNWPN (1)
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,YH,72,RUNID)
      YH=26.75
      CALL GRTXT (1.,YH,72,HEAD)
      YH=26.00
      CALL GRTXT (1.,YH,72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,YH,10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,YH-0.5,10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,YH-1.,10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,YH-2.,10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') VMAX
      CALL GRTXT (1.,YH-2.5,10,CH)
      CALL GRTXT (1.,YH-3.,10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') VMIN
      CALL GRTXT (1.,YH-3.5,10,CH)
C
      RETURN
      END
C
      SUBROUTINE RPSVEC (AORIG,BORIG,IBLD,ICURV,
     .                   IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A DATASET USED BY RAPS TO PRODUCE
C  A VECTOR PLOT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CPOLYG'
      INCLUDE 'CGRID'
      INCLUDE 'CTRIG'
      INCLUDE 'CGEOM'
      INCLUDE 'CLOGAU'
      INCLUDE 'CESTIM'
      INCLUDE 'CPLOT'
C
      DIMENSION AORIG(*),BORIG(*)
      DIMENSION YWERT(N1ST,N2ND+N3RD),ywert1(nrad,1)
      equivalence (ywert(1,1),ywert1(1,1))
      DIMENSION ZWERT(N1ST,N2ND+N3RD),zwert1(nrad,1)
      equivalence (zwert(1,1),zwert1(1,1))
      INTEGER ZUORD(NKNOT,0:20,2)
      DIMENSION XX(*),YY(*)
      REAL*4 XY(800)
      REAL*4 YH
C
      LOGICAL LOGL,TRC
      CHARACTER*72 TEXT1,HEAD,RUNID,TXHEAD
      CHARACTER*24 TEXT2,TEXT3
      CHARACTER*17 CH
c
      data nvplot/0/
C
      WRITE (60,*) RUNID
      WRITE (60,*) TXHEAD
      WRITE (60,*) HEAD
      WRITE (60,*) TEXT1
      WRITE (60,*) TEXT2
      WRITE (60,*) TEXT3
      WRITE (60,*)
c
      nraps2=80
      nvplot=nvplot+1
      ch(1:7)='vecplot'
      if (nvplot.lt.10) then
        write (ch(8:8),'(i1)') nvplot
        lench=8
      else
        write (ch(8:9),'(i2)') nvplot
        lench=9
      endif
C
      NRAPS2=NRAPS2+nvplot
      nraps=nraps+1
      IRAPS=IRAPS+1
C
C  WRITE VALUE OF THE VECTOR TO RAPS-FILE IN ORDER TO HAVE A
C  SHADED PLOT
      OPEN (UNIT=NRAPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS
C
      IF (LEVGEO.EQ.4) THEN
        CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD)
        CALL CELINT(BORIG,ZWERT1,LOGL,IBLD,ICURV,NRAD)
      ELSEIF (LEVGEO.LE.3) THEN
        CALL CELINT(AORIG,YWERT,LOGL,IBLD,ICURV,N1ST)
        CALL CELINT(BORIG,ZWERT,LOGL,IBLD,ICURV,N1ST)
      ENDIF
C
      IF (LEVGEO.LE.2.AND.NLTOR) THEN
        DO 1100 IR=1,NR1ST
          DO 3100 IT=1,NT3RD
            IF (ZMI.NE.666.) YWERT(IR,IT)=MAX(YWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IT)=MIN(YWERT(IR,IT),ZMA)
            IF (ZMI.NE.666.) ZWERT(IR,IT)=MAX(ZWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) ZWERT(IR,IT)=MIN(ZWERT(IR,IT),ZMA)
            BETRAG = SQRT(YWERT(IR,IT)**2+ZWERT(IR,IT)**2)
            WRITE (NRAPS,*) BETRAG
            YWERT(IR,IT)=YWERT(IR,IT)/(BETRAG+1.D-20)
            ZWERT(IR,IT)=ZWERT(IR,IT)/(BETRAG+1.D-20)
3100      CONTINUE
1100    CONTINUE
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL) THEN
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
              IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
              IF (ZMI.NE.666.) ZWERT(IR,IP)=MAX(ZWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) ZWERT(IR,IP)=MIN(ZWERT(IR,IP),ZMA)
              BETRAG = SQRT(YWERT(IR,IP)**2+ZWERT(IR,IP)**2)
              WRITE (NRAPS,*) BETRAG
              YWERT(IR,IP)=YWERT(IR,IP)/(BETRAG+1.D-20)
              ZWERT(IR,IP)=ZWERT(IR,IP)/(BETRAG+1.D-20)
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
      ELSEIF (LEVGEO.EQ.4) THEN
        DO 60 I=1,NRKNOT
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          IF (ZMI.NE.666.) ZWERT1(I,1)=MAX(ZWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) ZWERT1(I,1)=MIN(ZWERT1(I,1),ZMA)
          BETRAG = SQRT(YWERT1(I,1)**2+ZWERT1(I,1)**2)
          WRITE (NRAPS,*) BETRAG
          YWERT1(I,1)=YWERT1(I,1)/(BETRAG+1.D-20)
          ZWERT1(I,1)=ZWERT1(I,1)/(BETRAG+1.D-20)
60      CONTINUE
C
      ELSE
        WRITE (6,*) 'UNWRITTEN OPTION IN RPSCOL: PLOT ABANDONNED '
        RETURN
      ENDIF
C
      CLOSE (UNIT=NRAPS)
C
C  WRITE VECTOR-COMPONENTS TO RAPS-FILE IN ORDER TO FORM A VECTORPLOT
C
C
      OPEN (UNIT=NRAPS2,file=ch(1:lench),
     .                  ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS2

      WRITE(NRAPS2,'(1X,A5,8X,A4,50(11X,I1))') '-1111',
     .'PFEI',1,3,1,1,1
C
      IF (LEVGEO.LE.2.AND.NLTOR) THEN
        WRITE (6,*) 'UNWRITTEN OPTION IN RPSVEC: PLOT ABANDONNED '

      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL) THEN
C
C       IR+1,IP+1          IR,IP+1        IR-1,IP+1_
C           +----------------+----------------+
C           |        5       |        4       |
C           | 6              |               3|
C           |                |                |
C       IR+1,IP ---------- IR,IP -------- IR-1,IP
C           |                |                |
C           | 7              |              2 |
C           |        8       |        1       |
C           + ---------------+----------------+
C       IR+1,IP-1          IR,IP-1        IR-1,IP-1
C
        DO IR=1,NR1ST
          DO IPART=1,NPPLG
            IPLOOP: DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
C  RIGHT HEMISPHERE
              IF (IR.GT.1) THEN
C  LOWER RIGHT QUADRANT
                IF (IP.GT.NPOINT(1,IPART)) THEN
C  SIDE 1
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR,IP-1),YPOL(IR,IP-1),
     .                         XPOL(IR-1,IP-1),YPOL(IR-1,IP-1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 2
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR-1,IP-1),YPOL(IR-1,IP-1),
     .                         XPOL(IR-1,IP),YPOL(IR-1,IP),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
C  UPPER RIGHT QUADRANT
                IF (IP.LT.NPOINT(2,IPART)) THEN
C  SIDE 3
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR-1,IP),YPOL(IR-1,IP),
     .                         XPOL(IR-1,IP+1),YPOL(IR-1,IP+1),XL,XM)
                  if (xm .ge. 0. .and. xm .le. 1.) GOTO 100
C  SIDE 4
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR,IP+1),YPOL(IR,IP+1),
     .                         XPOL(IR-1,IP+1),YPOL(IR-1,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
              ENDIF
C  LEFT HEMISPHERE
              IF (IR.LT.NR1ST) THEN
C  UPPER LEFT QUADRANT
                IF (IP.LT.NPOINT(2,IPART)) THEN
C  SIDE 5
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP+1),YPOL(IR+1,IP+1),
     .                         XPOL(IR,IP+1),YPOL(IR,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 6
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP),YPOL(IR+1,IP),
     .                         XPOL(IR+1,IP+1),YPOL(IR+1,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
C  LOWER LEFT QUADRANT
                IF (IP.GT.NPOINT(1,IPART)) THEN
C  SIDE 7
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP),YPOL(IR+1,IP),
     .                         XPOL(IR+1,IP-1),YPOL(IR+1,IP-1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 8
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP-1),YPOL(IR+1,IP-1),
     .                         XPOL(IR,IP-1),YPOL(IR,IP-1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
              ENDIF
              CYCLE IPLOOP
100           CONTINUE
              ywert(IR,IP)=ywert(IR,IP)*xl*0.9
              zwert(IR,IP)=zwert(IR,IP)*xl*0.9
            ENDDO IPLOOP
          ENDDO
        ENDDO
        I=0
        DO IR=1,NR1ST
          DO IPART=1,NPPLG
            DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
              I=I+1
              BETRAG=SQRT(YWERT(IR,IP)**2+ZWERT(IR,IP)**2)
              IF (BETRAG .GT. 1.E-5)
     .        WRITE(nraps2,'(I6,1P,5E12.4)')
     .             I,YWERT(IR,IP),zwert(IR,IP),0.,0.,0.
            enddo
          enddo
        enddo

      ELSEIF (LEVGEO.EQ.4) THEN
        DO 41 I=1,NRKNOT
          DO 51 J=0,20
            DO 51 K=1,2
            ZUORD(I,J,K) = 0
51        CONTINUE
41      CONTINUE
        DO 40 J=1,NTRII
          DO 50 I=1,3
            ZUORD(NECKE(I,J),0,1) = ZUORD(NECKE(I,J),0,1) + 1
c zuord darf maximal 20 werden
            if (zuord(necke(i,j),0,1).gt.20) then
              write (6,*) 'error in rpsvec: zuord'
              call exit
            endif
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0,1),1) = J
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0,1),2) = I
50        CONTINUE
40      CONTINUE
        DO 61 I=1,NRKNOT
          IF (ZUORD(I,0,1).LT.1) THEN
            WRITE (6,*) 'ERROR IN RPSVEC: POINT ',I,' NOT IN MESH'
          ENDIF
          DO 70 J=1,ZUORD(I,0,1)
*  K: NUMBER OF TRIANGLE CONTAINING NODE I
            K=ZUORD(I,J,1)
*  IA: NUMBER OF NODE I IN TRIANGLE K
            ia=ZUORD(I,J,2)
c  ib, ic are the two other nodes in triangle k
            ib=ia+1
            if (ib.eq.4) ib=1
            ic=ib+1
            if (ic.eq.4) ic=1
c
            call muelam (xtrian(i),ytrian(i),ywert1(i,1),zwert1(i,1),
     .                   xtrian(necke(ib,k)),ytrian(necke(ib,k)),
     .                   xtrian(necke(ic,k)),ytrian(necke(ic,k)),xl,xm)
            if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) then
              ywert1(i,1)=ywert1(i,1)*xl*0.9
              zwert1(i,1)=zwert1(i,1)*xl*0.9
              goto 61
            endif
70        CONTINUE
c   no intersection found.
c   point I is on a boundary, and the vector is pointing outside
c   the computational volume. don't plot it.
          YWERT1(I,1) = 0.
          YWERT1(I,1) = 0.
61      continue

        DO I=1,NRKNOT
          BETRAG=SQRT(YWERT1(I,1)**2+ZWERT1(I,1)**2)
          IF (BETRAG .GT. 1.D-5)
     .    WRITE(nraps2,'(I6,1P,5E12.4)') I,YWERT1(I,1),ZWERT1(I,1),
     .                                   0.,0.,0.
        enddo
      ENDIF
C
      WRITE(NRAPS2,'(1X,A5,8X,A3,50(11X,I1))') '-9999',
     .           'FIN',0,0,0
      CLOSE (UNIT=NRAPS2)
C
      RETURN
      END
c
c
      subroutine muelam (x1,y1,vx,vy,x2,y2,x3,y3,xl,xm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      d=vx*(y2-y3)-vy*(x2-x3)+1.e-20
      xl=((y2-y3)*(x2-x1)-(x2-x3)*(y2-y1))/d
      xm=(-vy*(x2-x1)+vx*(y2-y1))/d

      return
      end
C
C
      SUBROUTINE XYPLOT (XY,NC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL*4 XY(*)
C
      NBEG=1
      IANF=3
1     CONTINUE
      DO 2 I=IANF+2,NC,2
        DXY=ABS(XY(IANF)-XY(I))+ABS(XY(IANF+1)-XY(I+1))
        IF (DXY.LT.1.E-3) THEN
          IF (I-IANF.GT.2) THEN
            HILF=XY(IANF+2)
            XY(IANF+2)=XY(I)
            XY(I)=HILF
            HILF=XY(IANF+3)
            XY(IANF+3)=XY(I+1)
            XY(I+1)=HILF
            IF (I-IANF.GT.4) THEN
              IF (MOD((I-IANF)/2,2).EQ.1) THEN
                J=I+2
              ELSE
                J=I-2
              ENDIF
              HILF=XY(IANF+4)
              XY(IANF+4)=XY(J)
              XY(J)=HILF
              HILF=XY(IANF+5)
              XY(IANF+5)=XY(J+1)
              XY(J+1)=HILF
            ENDIF
          ENDIF
          IANF=IANF+4
          GOTO 1
        ENDIF
2     CONTINUE
C
      CALL GRJMP (XY(NBEG),XY(NBEG+1))
      DO 5 I=NBEG+2,IANF,4
5       CALL GRDRW (XY(I),XY(I+1))
C
      IF (IANF.LT.NC-1) THEN
        NBEG=IANF+2
        IANF=IANF+4
        GOTO 1
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE CELINT (AORIG,YWERT,LOGL,IBLD,ICURV,N1DIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  INTERPOLATION FROM CELL CENTERS (AORIG)
C  TO QUANTITIES AT CELL VERTICES
C  INPUT : AORIG
C  OUTPUT: YWERT
C  IN CASE LOGL=.TRUE.: YWERT=LOG10(YWERT)
C
C    ******M******
C    *     *     *
C    *  D  I  E  *
C    *     *     *   |
C    L**H**A**F**J   |
C    *     *     *
C    *  C  G  B  *   IP
C    *     *     *
C    ******K******
C
C        <--- IR
C

      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CPOLYG'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CPLOT'
      INCLUDE 'CTRIG'
      INCLUDE 'CLOGAU'
      INCLUDE 'CCONA'
      INCLUDE 'CLGIN'
C
      DIMENSION AORIG(*)
      DIMENSION TEILA(4),TEILWERT(4),YWERT(N1DIM,*)
      INTEGER ZUORD(NKNOT,0:20)
      DOUBLE PRECISION IX,IY,JX,JY,KX,KY,LX,LY,MX,MY
C
C
      LOGICAL LOGL
      INTEGER AGEF,EGEF
      INTEGER AGLEICH(4),EGLEICH(4)
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.LPPOL3(IBLD)) THEN
        IP=1
        IF (NLPOL) IP=IPROJ3(IBLD,ICURV)
        IF (IP.LE.0.OR.IP.GT.NP2ND) IP=1
        DO 1100 IR=1,NR1ST
          DO 3100 IT=1,NT3RD
C   WERTEBEARBEITUNG
            DO 3101,J=1,4
              TEILA(J) = 0.
              TEILWERT(J) = 0.
3101        CONTINUE
C  UNTEN RECHTS
            IF ((IR .NE. 1) .AND. (IT .NE. 1)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM 1. POLYGON
C             UND IST NICHT ANFANG EINES POLYGONS
              AX = RHOSRF(IR)
              AY = ZSURF(IT)
              JX = RHOSRF(IR-1)
              KY = ZSURF(IT-1)
              GY = 0.5 * (AY + KY)
              FX = 0.5 * (AX + JX)
              IRD=IR-1+((IP-1)+(IT-2)*NP2T3)*NR1P2
              TEILWERT(1) = AORIG(IRD)
              TEILA(1) = ABS(AX-FX)*ABS(AY-GY)
            ENDIF
C  UNTEN LINKS
            IF ((IR .NE. NR1ST) .AND. (IT .NE. 1)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM LETZTEN
C             POLYGON UND IST NICHT ANFANGSPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = ZSURF(IT)
              LX = RHOSRF(IR+1)
              KY = ZSURF(IT-1)
              GY = 0.5 * (AY + KY)
              HX = 0.5 * (AX + LX)
              IRD=IR+((IP-1)+(IT-2)*NP2T3)*NR1P2
              TEILWERT(2) = AORIG(IRD)
              TEILA(2) = ABS(AX-HX)*ABS(AY-GY)
            ENDIF
C  OBEN LINKS
            IF ((IR .NE. NR1ST) .AND. (IT .NE. NT3RD)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM LETZTEN POLYGON
C             UND IST NICHT DER ENDPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = ZSURF(IT)
              LX = RHOSRF(IR+1)
              MY = ZSURF(IT+1)
              IY = 0.5 * (AY + MY)
              HX = 0.5 * (AX + LX)
              IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
              TEILWERT(3) = AORIG(IRD)
              TEILA(3) = ABS(AX-HX)*ABS(AY-IY)
            ENDIF
C  OBEN RECHTS
            IF ((IR .NE. 1) .AND. (IT .NE. NT3RD)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM 1. POLYGON
C             UND IST NICHT DER ENDPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = ZSURF(IT)
              JX = RHOSRF(IR-1)
              MY = ZSURF(IT+1)
              IY = 0.5 * (AY + MY)
              FX = 0.5 * (AX + JX)
              IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
              TEILWERT(4) = AORIG(IRD)
              TEILA(4) = ABS(AX-FX)*ABS(AY-IY)
            ENDIF
C
            AGES = TEILA(1) + TEILA(2) + TEILA(3) + TEILA(4)
            YWERT(IR,IT) = 0.
            DO 3103,J=1,4
              YWERT(IR,IT) = YWERT(IR,IT) + TEILA(J)/AGES*TEILWERT(J)
3103        CONTINUE
            IF (LOGL) YWERT(IR,IT)=LOG10(MAX(1.D-48,YWERT(IR,IT)))
3100      CONTINUE
1100    CONTINUE
C
      ELSEIF (LEVGEO.EQ.1.AND.LPTOR3(IBLD)) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 1110 IR=1,NR1ST
          DO 3110 IP=1,NP2ND
C   WERTEBEARBEITUNG
            DO 3111,J=1,4
              TEILA(J) = 0.
              TEILWERT(J) = 0.
3111        CONTINUE
C           UNTEN RECHTS
            IF ((IR .NE. 1) .AND. (IP .NE. 1)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM 1. POLYGON
C             UND IST NICHT ANFANG EINES POLYGONS
              AX = RHOSRF(IR)
              AY = PSURF(IP)
              JX = RHOSRF(IR-1)
              KY = PSURF(IP-1)
              GY = 0.5 * (AY + KY)
              FX = 0.5 * (AX + JX)
              IRD=IR-1+((IP-2)+(IT-1)*NP2T3)*NR1P2
              TEILWERT(1) = AORIG(IRD)
              TEILA(1) = ABS(AX-FX)*ABS(AY-GY)
            ENDIF
C  UNTEN LINKS
            IF ((IR .NE. NR1ST) .AND. (IP .NE. 1)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM LETZTEN
C             POLYGON UND IST NICHT ANFANGSPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = PSURF(IP)
              LX = RHOSRF(IR+1)
              KY = PSURF(IP-1)
              GY = 0.5 * (AY + KY)
              HX = 0.5 * (AX + LX)
              IRD=IR+((IP-2)+(IT-1)*NP2T3)*NR1P2
              TEILWERT(2) = AORIG(IRD)
              TEILA(2) = ABS(AX-HX)*ABS(AY-GY)
            ENDIF
C  OBEN LINKS
            IF ((IR .NE. NR1ST) .AND. (IP .NE. NP2ND)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM LETZTEN POLYGON
C             UND IST NICHT DER ENDPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = PSURF(IP)
              LX = RHOSRF(IR+1)
              MY = PSURF(IP+1)
              IY = 0.5 * (AY + MY)
              HX = 0.5 * (AX + LX)
              IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
              TEILWERT(3) = AORIG(IRD)
              TEILA(3) = ABS(AX-HX)*ABS(AY-IY)
            ENDIF
C  OBEN RECHTS
            IF ((IR .NE. 1) .AND. (IP .NE. NP2ND)) THEN
C             AKTUELLER PUNKT LIEGT NICHT AUF DEM 1. POLYGON
C             UND IST NICHT DER ENDPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = PSURF(IP)
              JX = RHOSRF(IR-1)
              MY = PSURF(IP+1)
              IY = 0.5 * (AY + MY)
              FX = 0.5 * (AX + JX)
              IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
              TEILWERT(4) = AORIG(IRD)
              TEILA(4) = ABS(AX-FX)*ABS(AY-IY)
            ENDIF
C
            AGES = TEILA(1) + TEILA(2) + TEILA(3) + TEILA(4)
            YWERT(IR,IP) = 0.
            DO 3113,J=1,4
              YWERT(IR,IP) = YWERT(IR,IP) + TEILA(J)/AGES*TEILWERT(J)
3113        CONTINUE
            IF (LOGL) YWERT(IR,IP)=LOG10(MAX(1.D-48,YWERT(IR,IP)))
3110      CONTINUE
1110    CONTINUE
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTOR3(IBLD)) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
C           RANDPUNKTE
            IRA=1
            IRE=NR1ST
            IPA=NPOINT(1,IPART)
            IPE=NPOINT(2,IPART)
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
C             WERTEBEARBEITUNG
              DO 31,J=1,4
                TEILA(J) = 0.
                TEILWERT(J) = 0.
31            CONTINUE
C             UNTEN RECHTS
              IF (IR .NE. IRA) THEN
C               AKTUELLER PUNKT LIEGT NICHT AUF DEM INNERSTEN POLYGON
                AX = XPOL(IR,IP)
                AY = YPOL(IR,IP)
                JX = XPOL(IR-1,IP)
                JY = YPOL(IR-1,IP)
                IF (IP .NE. IPA) THEN
C                 DER AKTUELLE PUNKT IST NICHT ANFANG EINES TEILPOLYGONS
                  KX = XPOL(IR,IP-1)
                  KY = YPOL(IR,IP-1)
                  BX = XCOM(IR-1,IP-1)
                  BY = YCOM(IR-1,IP-1)
                  IF (IR .NE. IRE) THEN
C                   DER AKTUELLE PUNKT LIEGT NICHT AUF DEM
C                   AEUSSERSTEN POLYGON
                    CX = XCOM(IR,IP-1)
                    CY = YCOM(IR,IP-1)
                    CALL SCHNITP(AX,AY,KX,KY,BX,BY,CX,CY,GX,GY)
                  ELSE
C                   DER AKTUELLE PUNKT LIEGT AUF DEM
C                   AEUSSERSTEN POLYGON
                    GX = 0.5 * (AX + KX)
                    GY = 0.5 * (AY + KY)
                  ENDIF
                  IF (IP .NE. IPE) THEN
C                   DER AKTUELLE PUNKT IST NICHT ENDE EINES TEILPOLYGONS
                    EX = XCOM(IR-1,IP)
                    EY = YCOM(IR-1,IP)
                    CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                    IRD=IR-1+((IP-2)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(1) = AORIG(IRD)
                    TEILA(1) = FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                  ELSE
C                   DER AKTUELLE PUNKT IST ENDE EINES TEILPOLYGON
C                   ANGRENZENDE TEILPOLYGONE SUCHEN
C                   ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    AGEF = 0
                    DO 1111 JPART=IPART+1,NPPLG
                      JP = NPOINT(1,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        AGEF = AGEF + 1
                        AGLEICH(AGEF) = JP
                      ENDIF
1111                CONTINUE
                    DO 2111 JPART=1,IPART
                      JP = NPOINT(1,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        AGEF = AGEF + 1
                        AGLEICH(AGEF) = JP
                      ENDIF
2111                CONTINUE
C                   ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    EGEF = 0
                    DO 1112 JPART = 1,NPPLG
                      IF (IPART .NE. JPART) THEN
                        JP = NPOINT(2,JPART)
                        IF ((AX - XPOL(IR,JP))**2 +
     >                      (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                          EGEF = EGEF + 1
                          EGLEICH(EGEF) = JP
                        ENDIF
                      ENDIF
1112                CONTINUE

                    IF (AGEF + EGEF .EQ. 3) THEN
C                     XPUNKT
C                     ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                     PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                     ENDPUNKTE
                      EX = XCOM(IR-1,AGLEICH(AGEF))
                      EY = YCOM(IR-1,AGLEICH(AGEF))
                      CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                      MX = XPOL(IR,AGLEICH(AGEF)+1)
                      MY = YPOL(IR,AGLEICH(AGEF)+1)
                      DX = XCOM(IR,AGLEICH(AGEF))
                      DY = YCOM(IR,AGLEICH(AGEF))
                      CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                      TEILA(1) = FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
     >                         + FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)
                      IRD=IR-1+((IP-2)+(IT-1)*NP2T3)*NR1P2
                      IRG=IR-1+((AGLEICH(AGEF)-1)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(1) = (AORIG(IRD)+AORIG(IRG))/2.
                    ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                     GENAU EIN ANGRENZENDES TEILPOLYGON
C                     GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                      EX = XCOM(IR-1,AGLEICH(1))
                      EY = YCOM(IR-1,AGLEICH(1))
                      CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                      IRD=IR-1+((IP-2)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(1)=AORIG(IRD)
                      TEILA(1) = FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                      JTEST = INMP2I(IR,IP,0)
                      IF (JTEST .NE. 0) THEN
                        JTEST = JTEST + NLIM
                        IF (ILIIN(JTEST) .GT. 0) THEN
C                         NICHT INTERPOLIEREN
                          FX = 0.5 * (AX + JX)
                          FY = 0.5 * (AY + JY)
                          IRD=IR-1+((IP-2)+(IT-1)*NP2T3)*NR1P2
                          TEILWERT(1)=AORIG(IRD)
                          TEILA(1)=FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                        ENDIF
                      ENDIF
                    ELSEIF (AGEF + EGEF .EQ. 0) THEN
C                     KEIN ANGRENZENDES TEILPOLYGON GEFUNDEN
                      FX = 0.5 * (AX + JX)
                      FY = 0.5 * (AY + JY)
                      IRD=IR-1+((IP-2)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(1) = AORIG(IRD)
                      TEILA(1) = FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                    ELSE
C  MECKERN
                    ENDIF
                  ENDIF
                ELSE
C                 DER AKTUELLE PUNKT IST ANFANG EINES TEILPOLYGONS
                  EX = XCOM(IR-1,IP)
                  EY = YCOM(IR-1,IP)
C                 ANGRENZENDE TEILPOLYGONE SUCHEN
C                 ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  AGEF = 0
                  DO 1120 JPART=1,NPPLG
                    IF (IPART .NE. JPART) THEN
                       JP = NPOINT(1,JPART)
                       IF ((AX - XPOL(IR,JP))**2 +
     >                     (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                         AGEF = AGEF + 1
                         AGLEICH(AGEF) = JP
                       ENDIF
                    ENDIF
1120              CONTINUE
C                 ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  EGEF = 0
                  DO 1121 JPART = IPART,NPPLG
                    JP = NPOINT(2,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      EGEF = EGEF + 1
                      EGLEICH(EGEF) = JP
                    ENDIF
1121              CONTINUE
                  DO 2121 JPART = 1,IPART-1
                    JP = NPOINT(2,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      EGEF = EGEF + 1
                      EGLEICH(EGEF) = JP
                    ENDIF
2121              CONTINUE
                  IF (AGEF + EGEF .EQ. 3) THEN
C                   XPUNKT
C                   ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                   PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                   ENDPUNKTE
                    BX = XCOM(IR-1,EGLEICH(1)-1)
                    BY = YCOM(IR-1,EGLEICH(1)-1)
                    KX = XPOL(IR,EGLEICH(1)-1)
                    KY = YPOL(IR,EGLEICH(1)-1)
                    CX = XCOM(IR,EGLEICH(1)-1)
                    CY = YCOM(IR,EGLEICH(1)-1)
                    LX = XPOL(IR+1,EGLEICH(1))
                    LY = YPOL(IR+1,EGLEICH(1))
                    DX = XCOM(IR,AGLEICH(1))
                    DY = YCOM(IR,AGLEICH(1))
                    MX = XPOL(IR,AGLEICH(1)+1)
                    MY = YPOL(IR,AGLEICH(1)+1)
                    EX = XCOM(IR-1,AGLEICH(1))
                    EY = YCOM(IR-1,AGLEICH(1))
                    CALL SCHNITP(AX,AY,KX,KY,BX,BY,CX,CY,GX,GY)
                    CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                    CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                    TEILA(1) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)+
     >                         FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                    IRD=IR+((EGLEICH(1)-2)+(IT-1)*NP2T3)*NR1P2
                    IRG=IR+((AGLEICH(1)-1)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(1) = (AORIG(IRD)+AORIG(IRG))/2.
                  ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                   GENAU EIN ANGRENZENDES TEILPOLYGON
C                   GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                    BX = XCOM(IR-1,EGLEICH(1)-1)
                    BY = YCOM(IR-1,EGLEICH(1)-1)
                    KX = XPOL(IR,EGLEICH(1)-1)
                    KY = YPOL(IR,EGLEICH(1)-1)
                    IF (IR .NE. IRE) THEN
C                     AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                     AEUSSERSTEN POLYGON
                      CX = XCOM(IR,EGLEICH(1)-1)
                      CY = YCOM(IR,EGLEICH(1)-1)
                      CALL SCHNITP(AX,AY,KX,KY,CX,CY,BX,BY,GX,GY)
                    ELSE
C                     AKTUELLER PUNKT LIEGT AUF DEM
C                     AEUSSERSTEN POLYGON
                      GX = 0.5 * (AX + KX)
                      GY = 0.5 * (AY + KY)
                    ENDIF
                    CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                    TEILA(1) = FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                    IRD=IR-1+((EGLEICH(1)-2)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(1) = AORIG(IRD)
                    JTEST = INMP2I(IR,IP,0)
                    IF (JTEST .NE. 0) THEN
                      JTEST = JTEST + NLIM
                      IF (ILIIN(JTEST) .GT. 0) THEN
C                       NICHT INTERPOLIEREN
                        TEILWERT(1) = 0.
                        TEILA(1) = 0.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
C  UNTEN LINKS
              IF (IR .NE. IRE) THEN
C               AKTUELLER PUNKT LIEGT NICHT AUF DEM AEUSSERSTEN
C               POLYGONS
                AX = XPOL(IR,IP)
                AY = YPOL(IR,IP)
                LX = XPOL(IR+1,IP)
                LY = YPOL(IR+1,IP)
                IF (IP .NE. IPA) THEN
C                 AKTUELLER PUNKT IST NICHT ANFANGSPUNKT EINES
C                 TEILPOLYGONS
                  KX = XPOL(IR,IP-1)
                  KY = YPOL(IR,IP-1)
                  CX = XCOM(IR,IP-1)
                  CY = YCOM(IR,IP-1)
                  IF (IR .NE. IRA) THEN
C                   AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                   INNERSTEN POLYGON
                    BX = XCOM(IR-1,IP-1)
                    BY = YCOM(IR-1,IP-1)
                    CALL SCHNITP(AX,AY,KX,KY,BX,BY,CX,CY,GX,GY)
                  ELSE
C                   AKTUELLER PUNKT LIEGT AUF DEM
C                   INNERSTEN POLYGON
                    GX = 0.5 * (AX + KX)
                    GY = 0.5 * (AY + KY)
                  ENDIF
                  IF (IP .NE. IPE) THEN
C                   AKTUELLER PUNKT IST NICHT ENDPUNKT EINES
C                   TEILPOLYGONS
                    DX = XCOM(IR,IP)
                    DY = YCOM(IR,IP)
                    CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                    IRD=IR+((IP-2)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(2) = AORIG(IRD)
                    TEILA(2) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)
                  ELSE
C                   AKTUELLER PUNKT IST ENDPUNKT EINES
C                   TEILPOLYGONS
C                   ANGRENZENDE TEILPOLYGONE SUCHEN
C                   ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    AGEF = 0
                    DO 110 JPART=IPART+1,NPPLG
                      JP = NPOINT(1,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        AGEF = AGEF + 1
                        AGLEICH(AGEF) = JP
                      ENDIF
110                 CONTINUE
                    DO 210 JPART=1,IPART
                      JP = NPOINT(1,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        AGEF = AGEF + 1
                        AGLEICH(AGEF) = JP
                      ENDIF
210                 CONTINUE
C                   ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    EGEF = 0
                    DO 111 JPART = 1,NPPLG
                      IF (IPART .NE. JPART) THEN
                        JP = NPOINT(2,JPART)
                        IF ((AX - XPOL(IR,JP))**2 +
     >                      (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                          EGEF = EGEF + 1
                          EGLEICH(EGEF) = JP
                        ENDIF
                      ENDIF
 111                CONTINUE
                    IF (AGEF + EGEF .EQ. 3) THEN
C                     XPUNKT
C                     ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                     PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                     ENDPUNKTE
                      DX = XCOM(IR,AGLEICH(1))
                      DY = YCOM(IR,AGLEICH(1))
                      CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                      MX = XPOL(IR,AGLEICH(1)+1)
                      MY = YPOL(IR,AGLEICH(1)+1)
                      EX = XCOM(IR-1,AGLEICH(1))
                      EY = YCOM(IR-1,AGLEICH(1))
                      CALL SCHNITP(AX,AY,MX,MY,DX,DY,EX,EY,IX,IY)
                      TEILA(2) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)+
     >                           FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                      IRD=IR+((IP-2)+(IT-1)*NP2T3)*NR1P2
                      IRG=IR+((AGLEICH(1)-1)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(2) = (AORIG(IRD)+AORIG(IRG))/2.
                    ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                     GENAU EIN ANGRENZENDES TEILPOLYGON
C                     GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                      DX = XCOM(IR,AGLEICH(1))
                      DY = YCOM(IR,AGLEICH(1))
                      CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                      IRD=IR+((IP-2)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(2)=AORIG(IRD)
                      TEILA(2) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)
                      JTEST = INMP2I(IR,IP,0)
                      IF (JTEST .NE. 0) THEN
                        JTEST = JTEST + NLIM
                        IF (ILIIN(JTEST) .GT. 0) THEN
C                         NICHT INTERPOLIEREN
                          HX = 0.5 * (AX + LX)
                          HY = 0.5 * (AY + LY)
                          IRD=IR+((IP-2)+(IT-1)*NP2T3)*NR1P2
                          TEILWERT(2)=AORIG(IRD)
                          TEILA(2)=FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)
                        ENDIF
                      ENDIF
                    ELSEIF (AGEF + EGEF .EQ. 0) THEN
C                     KEIN ANGRENZENDES TEILPOLYGON GEFUNDEN
                      HX = 0.5 * (AX + LX)
                      HY = 0.5 * (AY + LY)
                      IRD=IR+((IP-2)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(2) = AORIG(IRD)
                      TEILA(2) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)
                    ELSE
C   MECKERN
                    ENDIF
                  ENDIF
                ELSE
C                 AKTUELLER PUNKT IST ANFANGSPUNKT EINES
C                 TEILPOLYGONS
C                 ANGRENZENDE TEILPOLYGONE SUCHEN
C                 ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  DX = XCOM(IR,IP)
                  DY = YCOM(IR,IP)
                  AGEF = 0
                  DO 120 JPART=1,NPPLG
                    IF (IPART .NE. JPART) THEN
                      JP = NPOINT(1,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        AGEF = AGEF + 1
                        AGLEICH(AGEF) = JP
                      ENDIF
                    ENDIF
120               CONTINUE
C                 ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  EGEF = 0
                  DO 121 JPART = IPART,NPPLG
                    JP = NPOINT(2,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                 (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      EGEF = EGEF + 1
                      EGLEICH(EGEF) = JP
                    ENDIF
 121              CONTINUE
                  DO 221 JPART = 1,IPART-1
                    JP = NPOINT(2,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      EGEF = EGEF + 1
                      EGLEICH(EGEF) = JP
                    ENDIF
 221              CONTINUE
                  IF (AGEF + EGEF .EQ. 3) THEN
C                   XPUNKT
C                   ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                   PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                   ENDPUNKTE
                    CX = XCOM(IR,EGLEICH(EGEF)-1)
                    CY = YCOM(IR,EGLEICH(EGEF)-1)
                    KX = XPOL(IR,EGLEICH(EGEF)-1)
                    KY = YPOL(IR,EGLEICH(EGEF)-1)
                    BX = XCOM(IR-1,EGLEICH(EGEF)-1)
                    BY = YCOM(IR-1,EGLEICH(EGEF)-1)
                    JX = XPOL(IR-1,EGLEICH(EGEF))
                    JY = YPOL(IR-1,EGLEICH(EGEF))
                    EX = XCOM(IR-1,AGLEICH(1))
                    EY = YCOM(IR-1,AGLEICH(1))
                    MX = XPOL(IR,AGLEICH(1)+1)
                    MY = YPOL(IR,AGLEICH(1)+1)
                    DX = XCOM(IR,AGLEICH(1))
                    DY = YCOM(IR,AGLEICH(1))
                    CALL SCHNITP(AX,AY,KX,KY,CX,CY,BX,BY,GX,GY)
                    CALL SCHNITP(AX,AY,JX,JY,EX,EY,BX,BY,FX,FY)
                    CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                    TEILA(2) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)+
     >                         FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                    IRD=IR-1+((AGLEICH(1)-1)+(IT-1)*NP2T3)*NR1P2
                    IRG=IR-1+((EGLEICH(EGEF)-2)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(2)=(AORIG(IRD)+AORIG(IRG))/2.
                  ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                   GENAU EIN ANGRENZENDES TEILPOLYGON
C                   GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                    CX = XCOM(IR,EGLEICH(1)-1)
                    CY = YCOM(IR,EGLEICH(1)-1)
                    KX = XPOL(IR,EGLEICH(1)-1)
                    KY = YPOL(IR,EGLEICH(1)-1)
                    IF (IR .NE. IRA) THEN
C                     AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                     INNERSTEN POLYGON
                      BX = XCOM(IR-1,EGLEICH(1)-1)
                      BY = YCOM(IR-1,EGLEICH(1)-1)
                      CALL SCHNITP(AX,AY,KX,KY,CX,CY,BX,BY,GX,GY)
                    ELSE
C                     AKTUELLER PUNKT LIEGT AUF DEM
C                     INNERSTEN POLYGON
                      GX = 0.5 * (AX + KX)
                      GY = 0.5 * (AY + KY)
                    ENDIF
                    CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                    TEILA(2) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)
                    IRD=IR+((EGLEICH(1)-2)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(2) = AORIG(IRD)
                    JTEST = INMP2I(IR,IP,0)
                    IF (JTEST .NE. 0) THEN
                      JTEST = JTEST + NLIM
                      IF (ILIIN(JTEST) .GT. 0) THEN
C                         NICHT INTERPOLIEREN
                        TEILWERT(2) = 0.
                        TEILA(2) = 0.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
C  OBEN LINKS
              IF (IR .NE. IRE) THEN
C               AKTUELLER PUNKT LIEGT NICHT AUF DEM
C               AEUSSERSTEN POLYGON
                AX = XPOL(IR,IP)
                AY = YPOL(IR,IP)
                LX = XPOL(IR+1,IP)
                LY = YPOL(IR+1,IP)
                IF (IP .NE. IPE) THEN
C                 AKTUELLER PUNKT IST NICHT ENDPUNKT EINES
C                 TEILPOLYGONS
                  MX = XPOL(IR,IP+1)
                  MY = YPOL(IR,IP+1)
                  DX = XCOM(IR,IP)
                  DY = YCOM(IR,IP)
                  IF (IR .NE. IRA) THEN
C                   AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                   INNERSTEN POLYGON
                    EX = XCOM(IR-1,IP)
                    EY = YCOM(IR-1,IP)
                    CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                  ELSE
C                   AKTUELLER PUNKT LIEGT AUF DEM
C                   INNERSTEN POLYGON
                    IX = 0.5 * (AX + MX)
                    IY = 0.5 * (AY + MY)
                  ENDIF
                  IF (IP .NE. IPA) THEN
C                   AKTUELLER PUNKT IST NICHT ANFANGSPUNKT EINES
C                   TEILPOLYGONS
                    CX = XCOM(IR,IP-1)
                    CY = YCOM(IR,IP-1)
                    CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                    IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(3) = AORIG(IRD)
                    TEILA(3) = FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                  ELSE
C                   AKTUELLER PUNKT IST ANFANGSPUNKT EINES
C                   TEILPOLYGONS
C                   ANGRENZENDE TEILPOLYGONE SUCHEN
C                   ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    AGEF = 0
                    DO 140 JPART=1,NPPLG
                      IF (IPART .NE. JPART) THEN
                        JP = NPOINT(1,JPART)
                        IF ((AX - XPOL(IR,JP))**2 +
     >                     (AY-YPOL(IR,JP))**2.LT.EPS10)THEN
                          AGEF = AGEF + 1
                          AGLEICH(AGEF) = JP
                        ENDIF
                      ENDIF
140                 CONTINUE
C                   ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    EGEF = 0
                    DO 141 JPART = IPART,NPPLG
                      JP = NPOINT(2,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                         EGEF = EGEF + 1
                         EGLEICH(EGEF) = JP
                       ENDIF
 141                CONTINUE
                    DO 241 JPART = 1,IPART-1
                      JP = NPOINT(2,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        EGEF = EGEF + 1
                        EGLEICH(EGEF) = JP
                      ENDIF
 241                CONTINUE
                    IF (AGEF + EGEF .EQ. 3) THEN
C                     XPUNKT
C                     ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                     PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                     ENDPUNKTE
                      CX = XCOM(IR,EGLEICH(EGEF)-1)
                      CY = YCOM(IR,EGLEICH(EGEF)-1)
                      KX = XPOL(IR,EGLEICH(EGEF)-1)
                      KY = YPOL(IR,EGLEICH(EGEF)-1)
                      BX = XCOM(IR-1,EGLEICH(EGEF)-1)
                      BY = YCOM(IR-1,EGLEICH(EGEF)-1)
                      CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                      CALL SCHNITP(AX,AY,KX,KY,CX,CY,BX,BY,GX,GY)
                      TEILA(3) = FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)+
     >                           FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)
                      IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      IRG=IR+((EGLEICH(EGEF)-2)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(3)=(AORIG(IRD)+AORIG(IRG))/2.
                    ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                     GENAU EIN ANGRENZENDES TEILPOLYGON
C                     GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                      CX = XCOM(IR,EGLEICH(1)-1)
                      CY = YCOM(IR,EGLEICH(1)-1)
                      CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                      IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(3)=AORIG(IRD)
                      TEILA(3) = FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                      JTEST = INMP2I(IR,IP,0)
                      IF (JTEST .NE. 0) THEN
                        JTEST = JTEST + NLIM
                        IF (ILIIN(JTEST) .GT. 0) THEN
C                         NICHT INTERPOLIEREN
                          HX = 0.5 * (AX + LX)
                          HY = 0.5 * (AY + LY)
                          IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                          TEILWERT(3) = AORIG(IRD)
                          TEILA(3) = FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                        ENDIF
                      ENDIF
                    ELSEIF (AGEF + EGEF .EQ. 0) THEN
C                     KEIN ANGRENZENDES TEILPOLYGON GEFUNDEN
                      HX = 0.5 * (AX + LX)
                      HY = 0.5 * (AY + LY)
                      IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(3) = AORIG(IRD)
                      TEILA(3) = FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                    ELSE
C  MECKERN
                    ENDIF
                  ENDIF
                ELSE
C                 AKTUELLER PUNKT IST ENDPUNKT EINES
C                 TEILPOLYGONS
                  CX = XCOM(IR,IP-1)
                  CY = YCOM(IR,IP-1)
C                 ANGRENZENDE TEILPOLYGONE SUCHEN
C                 ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  AGEF = 0
                  DO 150 JPART=IPART+1,NPPLG
                    JP = NPOINT(1,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      AGEF = AGEF + 1
                      AGLEICH(AGEF) = JP
                    ENDIF
150               CONTINUE
                  DO 250 JPART=1,IPART
                    JP = NPOINT(1,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      AGEF = AGEF + 1
                      AGLEICH(AGEF) = JP
                    ENDIF
250               CONTINUE
C                 ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  EGEF = 0
                  DO 151 JPART = 1,NPPLG
                    IF (IPART .NE. JPART) THEN
                      JP = NPOINT(2,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        EGEF = EGEF + 1
                        EGLEICH(EGEF) = JP
                      ENDIF
                    ENDIF
 151              CONTINUE
                  IF (AGEF + EGEF .EQ. 3) THEN
C                   XPUNKT
C                   ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                   PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                   ENDPUNKTE
                    DX = XCOM(IR,AGLEICH(1))
                    DY = YCOM(IR,AGLEICH(1))
                    MX = XPOL(IR,AGLEICH(1)+1)
                    MY = YPOL(IR,AGLEICH(1)+1)
                    EX = XCOM(IR-1,AGLEICH(1))
                    EY = YCOM(IR-1,AGLEICH(1))
                    JX = XPOL(IR-1,AGLEICH(1))
                    JY = YPOL(IR-1,AGLEICH(1))
                    BX = XCOM(IR-1,EGLEICH(1)-1)
                    BY = YCOM(IR-1,EGLEICH(1)-1)
                    KX = XPOL(IR,EGLEICH(1)-1)
                    KY = YPOL(IR,EGLEICH(1)-1)
                    CX = XCOM(IR,EGLEICH(1)-1)
                    CY = YCOM(IR,EGLEICH(1)-1)
                    CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                    CALL SCHNITP(AX,AY,JX,JY,EX,EY,BX,BY,FX,FY)
                    CALL SCHNITP(AX,AY,KX,KY,BX,BY,CX,CY,GX,GY)
                    TEILA(3) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)+
     >                         FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                    IRD=IR-1+((AGLEICH(1)-1)+(IT-1)*NP2T3)*NR1P2
                    IRG=IR-1+((EGLEICH(1)-2)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(3)=(AORIG(IRD)+AORIG(IRG))/2.
                  ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                   GENAU EIN ANGRENZENDES TEILPOLYGON
C                   GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                    DX = XCOM(IR,AGLEICH(1))
                    DY = YCOM(IR,AGLEICH(1))
                    MX = XPOL(IR,AGLEICH(1)+1)
                    MY = YPOL(IR,AGLEICH(1)+1)
                    IF (IR .NE. IRA) THEN
C                     AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                     INNERSTEN POLYGON
                      EX = XCOM(IR-1,AGLEICH(1))
                      EY = YCOM(IR-1,AGLEICH(1))
                      CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                    ELSE
C                     AKTUELLER PUNKT LIEGT AUF DEM
C                     INNERSTEN POLYGON
                      IX = 0.5 * (AX + MX)
                      IY = 0.5 * (AY + MY)
                    ENDIF
                    CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                    TEILA(3) = FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                    IRD=IR+((AGLEICH(1)-1)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(3) = AORIG(IRD)
                    JTEST = INMP2I(IR,IP,0)
                    IF (JTEST .NE. 0) THEN
                      JTEST = JTEST + NLIM
                      IF (ILIIN(JTEST) .GT. 0) THEN
C                       NICHT INTERPOLIEREN
                        TEILWERT(3) = 0.
                        TEILA(3) = 0.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
C  OBEN RECHTS
              IF (IR .NE. IRA) THEN
C               AKTUELLER PUNKT LIEGT NICHT AUF DEM INNERSTEN POLYGON
                AX = XPOL(IR,IP)
                AY = YPOL(IR,IP)
                JX = XPOL(IR-1,IP)
                JY = YPOL(IR-1,IP)
                IF (IP .NE. IPE) THEN
C                 AKTUELLER PUNKT IST NICHT ENDPUNKT EINES TEILPOLYGONS
                  MX = XPOL(IR,IP+1)
                  MY = YPOL(IR,IP+1)
                  EX = XCOM(IR-1,IP)
                  EY = YCOM(IR-1,IP)
                  IF (IR .NE. IRE) THEN
C                   AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                   AEUSSERSTEN POLYGON
                    DX = XCOM(IR,IP)
                    DY = YCOM(IR,IP)
                    CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                  ELSE
C                   AKTUELLER PUNKT LIEGT AUF DEM
C                   AEUSSERSTEN POLYGON
                    IX = 0.5 * (AX + MX)
                    IY = 0.5 * (AY + MY)
                  ENDIF
                  IF (IP .NE. IPA) THEN
C                   AKTUELLER PUNKT IST NICHT ANFANGSPUNKT EINES
C                   TEILPOLYGONS
                    BX = XCOM(IR-1,IP-1)
                    BY = YCOM(IR-1,IP-1)
                    CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                    IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(4) = AORIG(IRD)
                    TEILA(4) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)
                  ELSE
C                   AKTUELLER PUNKT IST ANFANGSPUNKT EINES
C                   TEILPOLYGONS
C                   ANGRENZENDE TEILPOLYGONE SUCHEN
C                   ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    AGEF = 0
                    DO 160 JPART=1,NPPLG
                      IF (IPART .NE. JPART) THEN
                        JP = NPOINT(1,JPART)
                        IF ((AX - XPOL(IR,JP))**2 +
     >                     (AY-YPOL(IR,JP))**2.LT.EPS10)THEN
                          AGEF = AGEF + 1
                          AGLEICH(AGEF) = JP
                        ENDIF
                      ENDIF
160                 CONTINUE
C                   ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                   MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                    EGEF = 0
                    DO 161 JPART = IPART,NPPLG
                      JP = NPOINT(2,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        EGEF = EGEF + 1
                        EGLEICH(EGEF) = JP
                      ENDIF
 161                CONTINUE
                    DO 261 JPART = 1,IPART-1
                      JP = NPOINT(2,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        EGEF = EGEF + 1
                        EGLEICH(EGEF) = JP
                      ENDIF
 261                CONTINUE
                    IF (AGEF + EGEF .EQ. 3) THEN
C                     XPUNKT
C                     ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                     PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                     ENDPUNKTE
                      BX = XCOM(IR-1,EGLEICH(1)-1)
                      BY = YCOM(IR-1,EGLEICH(1)-1)
                      KX = XPOL(IR,EGLEICH(1)-1)
                      KY = YPOL(IR,EGLEICH(1)-1)
                      CX = XCOM(IR,EGLEICH(1)-1)
                      CY = YCOM(IR,EGLEICH(1)-1)
                      CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                      CALL SCHNITP(AX,AY,KX,KY,BX,BY,CX,CY,GX,GY)
                      TEILA(4) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)+
     >                           FLAECH(AX,AY,FX,FY,BX,BY,GX,GY)
                      IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      IRG=IR-1+((EGLEICH(1)-2)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(4)=(AORIG(IRD)+AORIG(IRG))/2.
                    ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                     GENAU EIN ANGRENZENDES TEILPOLYGON
C                     GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                      BX = XCOM(IR-1,EGLEICH(1)-1)
                      BY = YCOM(IR-1,EGLEICH(1)-1)
                      CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                      IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(4)=AORIG(IRD)
                      TEILA(4) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)
                      JTEST = INMP2I(IR,IP,0)
                      IF (JTEST .NE. 0) THEN
                        JTEST = JTEST + NLIM
                        IF (ILIIN(JTEST) .GT. 0) THEN
C                         NICHT INTERPOLIEREN
                          FX = 0.5 * (AX + JX)
                          FY = 0.5 * (AY + JY)
                          IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                          TEILWERT(4)=AORIG(IRD)
                          TEILA(4)=FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)
                        ENDIF
                      ENDIF
                    ELSEIF (AGEF + EGEF .EQ. 0) THEN
C                     KEIN ANGRENZENDES TEILPOLYGON GEFUNDEN
                      FX = 0.5 * (AX + JX)
                      FY = 0.5 * (AY + JY)
                      IRD=IR-1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      TEILWERT(4) = AORIG(IRD)
                      TEILA(4) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)
                    ELSE
C    MECKERN
                    ENDIF
                  ENDIF
                ELSE
C                 AKTUELLER PUNKT IST ENDPUNKT EINES
C                 TEILPOLYGONS
                  BX = XCOM(IR-1,IP-1)
                  BY = YCOM(IR-1,IP-1)
C                 ANGRENZENDE TEILPOLYGONE SUCHEN
C                 ANFANGSPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  AGEF = 0
                  DO 170 JPART=IPART+1,NPPLG
                    JP = NPOINT(1,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      AGEF = AGEF + 1
                      AGLEICH(AGEF) = JP
                    ENDIF
170               CONTINUE
                  DO 270 JPART=1,IPART
                    JP = NPOINT(1,JPART)
                    IF ((AX - XPOL(IR,JP))**2 +
     >                  (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                      AGEF = AGEF + 1
                      AGLEICH(AGEF) = JP
                    ENDIF
270               CONTINUE
C                 ENDPUNKTE DER TEILPOLYGONE AUF GLEICHHEIT
C                 MIT DEM AKTUELLEN PUNKT UEBERPRUFEN
                  EGEF = 0
                  DO 171 JPART = 1,NPPLG
                    IF (IPART .NE. JPART) THEN
                      JP = NPOINT(2,JPART)
                      IF ((AX - XPOL(IR,JP))**2 +
     >                    (AY - YPOL(IR,JP))**2 .LT. EPS10) THEN
                        EGEF = EGEF + 1
                        EGLEICH(EGEF) = JP
                      ENDIF
                    ENDIF
 171              CONTINUE
                  IF (AGEF + EGEF .EQ. 3) THEN
C                   XPUNKT
C                   ES EXISTIEREN GENAU 3 MIT DEM AKTUELLEN
C                   PUNKT UEBEREINSTIMMENDE ANFANGS- BZW.
C                   ENDPUNKTE
                    EX = XCOM(IR-1,AGLEICH(AGEF))
                    EY = YCOM(IR-1,AGLEICH(AGEF))
                    MX = XPOL(IR,AGLEICH(AGEF)+1)
                    MY = YPOL(IR,AGLEICH(AGEF)+1)
                    DX = XCOM(IR,AGLEICH(AGEF))
                    DY = YCOM(IR,AGLEICH(AGEF))
                    LX = XPOL(IR+1,AGLEICH(AGEF))
                    LY = YPOL(IR+1,AGLEICH(AGEF))
                    CX = XCOM(IR,EGLEICH(1)-1)
                    CY = YCOM(IR,EGLEICH(1)-1)
                    BX = XCOM(IR-1,EGLEICH(1)-1)
                    BY = YCOM(IR-1,EGLEICH(1)-1)
                    KX = XPOL(IR,EGLEICH(1)-1)
                    KY = YPOL(IR,EGLEICH(1)-1)
                    CALL SCHNITP(AX,AY,KX,KY,CX,CY,BX,BY,GX,GY)
                    CALL SCHNITP(AX,AY,LX,LY,CX,CY,DX,DY,HX,HY)
                    CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                    TEILA(4) = FLAECH(AX,AY,GX,GY,CX,CY,HX,HY)+
     >                         FLAECH(AX,AY,HX,HY,DX,DY,IX,IY)
                    IRD=IR+((EGLEICH(1)-2)+(IT-1)*NP2T3)*NR1P2
                    IRG=IR+((AGLEICH(AGEF)-1)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(4)=(AORIG(IRD)+AORIG(IRG))/2.
                  ELSEIF (AGEF + EGEF .EQ. 1) THEN
C                   GENAU EIN ANGRENZENDES TEILPOLYGON
C                   GEFUNDEN, UEBER GRENZE INTERPOLIEREN
                    EX = XCOM(IR-1,AGLEICH(1))
                    EY = YCOM(IR-1,AGLEICH(1))
                    MX = XPOL(IR,AGLEICH(1)+1)
                    MY = YPOL(IR,AGLEICH(1)+1)
                    IF (IR .NE. IRE) THEN
C                     AKTUELLER PUNKT LIEGT NICHT AUF DEM
C                     AEUSSERSTEN POLYGON
                      DX = XCOM(IR,AGLEICH(1))
                      DY = YCOM(IR,AGLEICH(1))
                      CALL SCHNITP(AX,AY,MX,MY,EX,EY,DX,DY,IX,IY)
                    ELSE
C                     AKTUELLER PUNKT LIEGT AUF DEM
C                     AEUSSERSTEN POLYGON
                      IX = 0.5 * (AX + MX)
                      IY = 0.5 * (AY + MY)
                    ENDIF
                    CALL SCHNITP(AX,AY,JX,JY,BX,BY,EX,EY,FX,FY)
                    TEILA(4) = FLAECH(AX,AY,IX,IY,EX,EY,FX,FY)
                    IRD=IR-1+((AGLEICH(1)-1)+(IT-1)*NP2T3)*NR1P2
                    TEILWERT(4) = AORIG(IRD)
                    JTEST = INMP2I(IR,IP,0)
                    IF (JTEST .NE. 0) THEN
                      JTEST = JTEST + NLIM
                      IF (ILIIN(JTEST) .GT. 0) THEN
C                       NICHT INTERPOLIEREN
                        TEILWERT(4) = 0.
                        TEILA(4) = 0.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
C
              AGES = TEILA(1) + TEILA(2) + TEILA(3) + TEILA(4)
              YWERT(IR,IP) = 0.
              DO 33,J=1,4
                YWERT(IR,IP) = YWERT(IR,IP) + TEILA(J)/AGES*TEILWERT(J)
33            CONTINUE
              IF (LOGL) YWERT(IR,IP)=LOG10(MAX(1.D-48,YWERT(IR,IP)))
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
C
        DO 41 I=1,NRKNOT
          YWERT(I,1) = 0.
          DO 51 J=0,20
            ZUORD(I,J) = 0
51        CONTINUE
41      CONTINUE
        DO 40 J=1,NTRII
          DO 50 I=1,3
            ZUORD(NECKE(I,J),0) = ZUORD(NECKE(I,J),0) + 1
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0)) = J
50        CONTINUE
40      CONTINUE
        DO 60 I=1,NRKNOT
          GESA = 0
          DO 70 J=1,ZUORD(I,0)
            K=ZUORD(I,J)
            GESA = GESA + VOL(K)
            YWERT(I,1) = YWERT(I,1) + VOL(K) * AORIG(K)
70        CONTINUE
          IF (GESA .NE. 0) THEN
            YWERT(I,1) = YWERT(I,1)/GESA
          ELSE
            WRITE (6,*) 'ERROR IN CELINT: POINT ',I,' NOT IN MESH'
            YWERT(I,1) = YWERT(I,1)
          ENDIF
          IF (LOGL) YWERT(I,1)=LOG10(MAX(1.D-48,YWERT(I,1)))
60      CONTINUE
C
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE RPSCOL (AORIG,IBLD,ICURV,
     .                   I1,I2,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C    ******M******
C    *     *     *
C    *  D  I  E  *
C    *     *     *   |
C    L**H**A**F**J   |
C    *     *     *
C    *  C  G  B  *   IP
C    *     *     *
C    ******K******
C
C        <--- IR
C

      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CPOLYG'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CPLOT'
      INCLUDE 'CTRIG'
      INCLUDE 'CLOGAU'
C
      DIMENSION AORIG(*)
      DIMENSION YWERT(N1ST,N2ND+N3RD),ywert1(nrad,1)
      equivalence (ywert(1,1),ywert1(1,1))
      DIMENSION XX(*),YY(*)
C
      LOGICAL LOGL,TRC
C
      CHARACTER*72 TEXT1,HEAD,RUNID,TXHEAD
      CHARACTER*24 TEXT2,TEXT3
C
      WRITE (60,*) RUNID
      WRITE (60,*) TXHEAD
      WRITE (60,*) HEAD
      WRITE (60,*) TEXT1
      WRITE (60,*) TEXT2
      WRITE (60,*) TEXT3
      WRITE (60,*)
C
      NRAPS=NRAPS+1
      IRAPS=IRAPS+1
C
      OPEN (UNIT=NRAPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS
C
      IF (LEVGEO.EQ.4) THEN
        CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD)
      ELSEIF (LEVGEO.LE.3) THEN
        CALL CELINT(AORIG,YWERT,LOGL,IBLD,ICURV,N1ST)
      ENDIF
C
      IF (LEVGEO.LE.2.AND.LPPOL3(IBLD)) THEN
        LPPOLR=.TRUE.
        IF (.NOT.NLRAD.OR..NOT.NLTOR) THEN
          WRITE (6,*) 'ERROR IN RPSCOL. LPPOL3? '
          RETURN
        ENDIF
        DO 1100 IR=1,NR1ST
          DO 3100 IT=1,NT3RD
            IF (ZMI.NE.666.) YWERT(IR,IT)=MAX(YWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IT)=MIN(YWERT(IR,IT),ZMA)
            WRITE (NRAPS,*) YWERT(IR,IT)
3100      CONTINUE
1100    CONTINUE
C
      ELSEIF (LEVGEO.LE.2.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        IF (.NOT.NLRAD.OR..NOT.NLTOR) THEN
          WRITE (6,*) 'ERROR IN RPSCOL. LPTOR3? '
          RETURN
        ENDIF
        DO 1 IR=1,NR1ST
          DO 3 IP=1,NP2ND
            IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
            WRITE (NRAPS,*) YWERT(IR,IP)
3         CONTINUE
1       CONTINUE
C
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        IF (.NOT.NLRAD.OR..NOT.NLPOL) THEN
          WRITE (6,*) 'ERROR IN RPSCOL. LPTOR3? '
          RETURN
        ENDIF
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
              IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
              WRITE (NRAPS,*) YWERT(IR,IP)
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        DO 60 I=1,NRKNOT
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          WRITE(NRAPS,*) YWERT1(I,1)
60      CONTINUE
C
      ELSE
        WRITE (6,*) 'UNWRITTEN OPTION IN RPSCOL: PLOT ABANDONNED '
      ENDIF
C
      CLOSE (UNIT=NRAPS)
C
      RETURN
      END
C
      SUBROUTINE RPSOUT
C
C  ANZ: NUMBER OF CELLS
C  WRITE (17,...) LABELED CO-ORDINATES OF VERTICES
C  WRITE (18,...) LABELING NUMBER OF VERTICES FOR EACH CELL
C  WRITE (19,...) VALUE OF TALLY AT VERTEX (BY LABEL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGEOM'
      INCLUDE 'CPOLYG'
      INCLUDE 'CRECH'
      INCLUDE 'COMUSR'
      INCLUDE 'CPLOT'
      INCLUDE 'CGRID'
      INCLUDE 'CTRIG'
      INCLUDE 'CLOGAU'
      INCLUDE 'CLGIN'
      DIMENSION YWERT(2*NPTAL)
      INTEGER ANZ
C
      I=0
      NSTAB=0
2     CONTINUE
      I=I+1
      IF (I+1.LE.INSTOR) THEN
        IF (NPL2D(I).EQ.0) THEN
C  START OF A NEW LINE
          NSTAB=NSTAB+1
        ELSEIF (NPL2D(I+1).EQ.1) THEN
C  CONTINUATION OF A LINE
          NSTAB=NSTAB+1
        ENDIF
        GOTO 2
      ENDIF
C
      IF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
        NST=NSTSI
        IF (NTIME.GT.0) NST=NST-1
        DO ISTS=1,NST
          IF (ILIIN(NLIM+ISTS) > 0) THEN
            IDIMP=TRANSFER(MAXLOC(INUMP(ISTS,1:3)),1)
            IF (IDIMP == 1) THEN
              IA=IRPTA(ISTS,2)
              IE=IRPTE(ISTS,2)
            ELSEIF (IDIMP == 2) THEN
              IA=IRPTA(ISTS,1)
              IE=IRPTE(ISTS,1)
            END IF
            IF ((IDIMP>0) .AND. (IDIMP<=2))  NSTAB=NSTAB+(IE-IA)
          END IF
        END DO
      END IF
C
      DO 5 IF=1,IRAPS
        NRPS=60+IF
        OPEN (UNIT=NRPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
        REWIND NRPS
5     CONTINUE
C
      WRITE(17,'(1X,A5,8X,A4,11X,A1,11X,A1,11X,A1,11X,A1)') '-1111',
     .'NPCO','1','2','1','1'
      WRITE(19,'(1X,A5,8X,A4,50(11X,I1))') '-1111',
     .'NPST',1,IRAPS,1,(1,I=1,IRAPS)
C
      IF (LEVGEO.EQ.1.AND.LPTORR) THEN
        ANZ = NR1STM*NP2NDM
C
        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'QUAM4        1',ANZ,'4'
        I=0
        DO 100 IR=1,NR1ST
           DO 110 IP=1,NP2ND
             IF (IP .NE. NP2ND) THEN
               I = I + 1
               IF (IR .LT. NR1ST) THEN
                 WRITE(18,'(1X,A1,4I6)') '0',
     .                     (IR-1)*NP2ND+IP,
     .                     (IR-1)*NP2ND+IP+1,
     .                     IR*NP2ND+IP+1,
     .                     IR*NP2ND+IP
               ENDIF
               WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR),PSURF(IP)
             ENDIF
             DO 105 IF=1,IRAPS
               READ (60+IF,*) YWERT(IF)
105          CONTINUE
             IF (IP .NE. NP2ND) THEN
               WRITE(19,'(I6,1P,50E12.4)') I,  (YWERT(IF),IF=1,IRAPS)
             ELSE
               WRITE(19,'(I6,1P,50E12.4)') I+1,(YWERT(IF),IF=1,IRAPS)
             ENDIF
110        CONTINUE
           I = I + 1
           WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR),PSURF(NP2ND)
100     CONTINUE
        NCO=I
C
      ELSEIF (LEVGEO.LE.2.AND.LPPOLR) THEN
        ANZ = NR1STM*NT3RDM
C
        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'QUAM4        1',ANZ,'4'
        I=0
        DO 1100 IR=1,NR1ST
           DO 2100 IT=1,NT3RD
             IF (IT .NE. NT3RD) THEN
               I = I + 1
               IF (IR .LT. NR1ST) THEN
                 WRITE(18,'(1X,A1,4I6)') '0',
     .                     (IR-1)*NT3RD+IT,
     .                     (IR-1)*NT3RD+IT+1,
     .                     IR*NT3RD+IT+1,
     .                     IR*NT3RD+IT
               ENDIF
               WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR),ZSURF(IT)
             ENDIF
             DO 2105 IF=1,IRAPS
               READ (60+IF,*) YWERT(IF)
2105         CONTINUE
             IF (IT .NE. NT3RD) THEN
               WRITE(19,'(I6,1P,50E12.4)') I,  (YWERT(IF),IF=1,IRAPS)
             ELSE
               WRITE(19,'(I6,1P,50E12.4)') I+1,(YWERT(IF),IF=1,IRAPS)
             ENDIF
2100       CONTINUE
           I = I + 1
           WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR),ZSURF(NT3RD)
1100    CONTINUE
        NCO=I
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
        ANZ = 0
        DO 1 IR=1,NR1ST-1
          DO 1 IPPLG=1,NPPLG
             DO 1 IP=NPOINT(1,IPPLG),NPOINT(2,IPPLG)-1
               ANZ = ANZ + 1
1       CONTINUE

        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'QUAM4        1',ANZ,'4'
        I=0
        DO 10 IR=1,NR1ST
          NPUNKT = 0
          DO 11 IPPLG=1,NPPLG
            NPUNKT = NPUNKT+NPOINT(2,IPPLG)-NPOINT(1,IPPLG)+1
11        CONTINUE
          DO 20 IPPLG=1,NPPLG
            DO 30 IP=NPOINT(1,IPPLG),NPOINT(2,IPPLG)
              IF (IP .NE. NPOINT(2,IPPLG)) THEN
                I = I + 1
                IF (IR .LT. NR1ST) THEN
                  WRITE(18,'(1X,A1,4I6)') '0',
     .                      (IR-1)*NPUNKT+IP,
     .                      (IR-1)*NPUNKT+IP+1,
     .                      IR*NPUNKT+IP+1,
     .                      IR*NPUNKT+IP
                ENDIF
                WRITE(17,'(I6,1P,2E12.4)') I,XPOL(IR,IP),YPOL(IR,IP)
              ENDIF
              DO 25 IF=1,IRAPS
                READ (60+IF,*) YWERT(IF)
25            CONTINUE
              IF (IP .NE. NPOINT(2,IPPLG)) THEN
                WRITE(19,'(I6,1P,50E12.4)') I,(YWERT(IF),IF=1,IRAPS)
              ELSE
                WRITE(19,'(I6,1P,50E12.4)') I+1,(YWERT(IF),IF=1,IRAPS)
              ENDIF
30          CONTINUE
            I = I + 1
            WRITE(17,'(I6,1P,2E12.4)') I,XPOL(IR,NPOINT(2,IPPLG)),
     .                                   YPOL(IR,NPOINT(2,IPPLG))
20        CONTINUE
10      CONTINUE
        NCO=I
C
      ELSEIF (LEVGEO.EQ.4.AND.LPTORR) THEN
        ANZ=NTRII
        DO 40 I=1,NRKNOT
          DO 50 IF=1,IRAPS
            READ(60+IF,*) YWERT(IF)
50        CONTINUE
          WRITE(19,'(I6,1P,50E12.4)') I,(YWERT(IF),IF=1,IRAPS)
          WRITE(17,'(I6,1P,2E12.4)') I,XTRIAN(I),YTRIAN(I)
40      CONTINUE
        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'TRIM3        1',ANZ,'3'

        DO 60 I=1,NTRII
          WRITE(18,'(1X,A1,3I6)') '0',NECKE(1,I),NECKE(2,I),NECKE(3,I)
60      CONTINUE
        NCO=NRKNOT
      ELSE
      ENDIF
c
      IGR=1
      IF (NCO > 0) IGR=2
      WRITE(18,'(A8,2I6,5X,A1)') 'FLA2    ',IGR,NSTAB,'2'
c
      I=0
3     CONTINUE
      I=I+1
      IF (I+1.LE.INSTOR) THEN
        IF (NPL2D(I).EQ.0) THEN
C  START OF A NEW LINE
          NCO=NCO+1
          WRITE (17,'(I6,1P,2E12.4)') NCO,XPL2D(I),YPL2D(I)
          NCO=NCO+1
          WRITE (17,'(I6,1P,2E12.4)') NCO,XPL2D(I+1),YPL2D(I+1)
          WRITE (18,'(1X,A1,2I6)') '0',NCO-1,NCO
        ELSEIF (NPL2D(I+1).EQ.1) THEN
C  CONTINUATION OF A LINE
          NCO=NCO+1
          WRITE (17,'(I6,1P,2E12.4)') NCO,XPL2D(I),YPL2D(I)
          NCO=NCO+1
          WRITE (17,'(I6,1P,2E12.4)') NCO,XPL2D(I+1),YPL2D(I+1)
          WRITE (18,'(1X,A1,2I6)') '0',NCO-1,NCO
        ENDIF
        GOTO 3
      ENDIF
C
      IF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
        NST=NSTSI
        IF (NTIME.GT.0) NST=NST-1
        DO ISTS=1,NST
          IF (ILIIN(NLIM+ISTS) > 0) THEN
          IDIMP=TRANSFER(MAXLOC(INUMP(ISTS,1:3)),1)

          IF (IDIMP == 1) THEN
            IR=INUMP(ISTS,1)
            IA=IRPTA(ISTS,2)
            IE=IRPTE(ISTS,2)
            DO IP=IA,IE-1
              NCO=NCO+1
              WRITE (17,'(I6,1P,2E12.4)') NCO,XPOL(IR,IP),YPOL(IR,IP)
              NCO=NCO+1
              WRITE (17,'(I6,1P,2E12.4)')
     .                               NCO,XPOL(IR,IP+1),YPOL(IR,IP+1)
              WRITE (18,'(1X,A1,2I6)') '0',NCO-1,NCO
            END DO

          ELSEIF (IDIMP ==2) THEN
            IP=INUMP(ISTS,IDIMP)
            IA=IRPTA(ISTS,1)
            IE=IRPTE(ISTS,1)
            DO IR=IA,IE-1
              NCO=NCO+1
              WRITE (17,'(I6,1P,2E12.4)') NCO,XPOL(IR,IP),YPOL(IR,IP)
              NCO=NCO+1
              WRITE (17,'(I6,1P,2E12.4)')
     .                               NCO,XPOL(IR+1,IP),YPOL(IR+1,IP)
              WRITE (18,'(1X,A1,2I6)') '0',NCO-1,NCO
            END DO
          END IF
          END IF
        END DO
      END IF
c
      WRITE(17,'(1X,A5,8X,A3,12X,A1,11X,A1,11X,A1,11X,A1)') '-9999',
     .           'FIN','0','0','0','0'
      WRITE(19,'(1X,A5,8X,A3,50(11X,I1))') '-9999',
     .           'FIN',(0,IF=1,IRAPS)
      RETURN
      END

      SUBROUTINE SCHNITP(X1,Y1,X2,Y2,X3,Y3,X4,Y4,EX,EY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION MUE

      MUE = ((Y2-Y4)*(X3-X4)+(X4-X2)*(Y3-Y4))/
     .      ((X1-X2)*(Y3-Y4)-(Y1-Y2)*(X3-X4)+1.D-20)
      EX = X1 + MUE * (X1-X2)
      EY = Y1 + MUE * (Y1-Y2)
      END


      FUNCTION FLAECH (X1,Y1,X2,Y2,X3,Y3,X4,Y4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FLAECH1 = 0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
      FLAECH2 = 0.5*(X1*(Y3-Y4)+X3*(Y4-Y1)+X4*(Y1-Y3))
      FLAECH = ABS(FLAECH1) + ABS(FLAECH2)
      END
c
C
      SUBROUTINE STCOOR (X,Y,IFL)
c  instor: zaehler, instor-ter aufruf
c  ifl: flag:=0, first point, .ne.0: else, npl2d(instor)=ifl
c  x,y: co-ordinaten des punktes, stored on xpl2d(instor),ypl2d(instor)
c
c  alle teilstuecke stehen hintereinander auf xpl2d,ypl2d
c  zum entwirren:NUMSUR(inums,..) array
c  inums: zaehler fuer individuelle zu plottende kurvenstuecke
c  jedes bekommt die eigene nummer und arrow
c  daher: ggfls: eine nummer mehrmals im plot vorhanden
c  NUMSUR(1,inums) : nummer des flachestuecks
c  NUMSUR(2,inums) : index des anfangspunktes in xpl2d/ypl2d arrays
c  NUMSUR(3,inums) : index des endpunktes in xpl2d/ypl2d arrays
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CRECH'
      DATA IFIRST/0/
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        INUMS=0
      ENDIF
C
      IF (INSTOR.LT.NRECH) THEN
        INSTOR=INSTOR+1
        XPL2D(INSTOR)=X
        YPL2D(INSTOR)=Y
        NPL2D(INSTOR)=IFL
        IF (INSTOR.EQ.NRECH) THEN
          WRITE (6,*) ' MAXIMUM NUMBER OF CO-ORDINATES STORED IN STCOOR'
          WRITE (6,*) ' IS REACHED.'
          WRITE (6,*) ' STORING OF FURTHER CO-ORDINATES IS STOPPED '
          WRITE (6,*) ' INCREASE PARAMETER NRECH AND RECOMPILE '
        ENDIF
C
        IF (IFL.EQ.0) THEN
          INUMS=INUMS+1
          NUMSUR(INUMS,1)=INOSF
          NUMSUR(INUMS,2)=INSTOR
        ENDIF
        NUMSUR(INUMS,3)=INSTOR
      ENDIF
C
      RETURN
      END

      SUBROUTINE ELLIPSOID (X0,Y0,Z0,CX,CY,CZ,XLIMS1,YLIMS1,ZLIMS1,
     .                      XLIMS2,YLIMS2,ZLIMS2,RLB,ILCOL,NX,NY,NZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'CCONA'
      DIMENSION XP(2*101), YP(2*101), PHIAN(9), PHIEN(9)
      INTEGER IP(2*101)

      XA = MAX(X0-CX,XLIMS1)
      XE = MIN(X0+CX,XLIMS2)
      IF (NX > 1) THEN
        DX = (XE-XA)/FLOAT(NX-1)
      ELSE
        DX = 0.D0
      END IF

      YA = MAX(Y0-CY,YLIMS1)
      YE = MIN(Y0+CY,YLIMS2)
      IF (NY > 1) THEN
        DY = (YE-YA)/FLOAT(NY-1)
      ELSE
        DY = 0.D0
      END IF

      ZA = MAX(Z0-CZ,ZLIMS1)
      ZE = MIN(Z0+CZ,ZLIMS2)
      IF (NZ > 1) THEN
        DZ = (ZE-ZA)/FLOAT(NZ-1)
      ELSE
        DZ = 0.D0
      END IF

      N = 101

!  PLOT Z-ISOLINES
      CALL GRNWPN(ILCOL)

      DO IZ = 1, NZ
        Z = ZA + (IZ-1)*DZ
        RAD = 1.D0 - (Z-Z0)**2/CZ**2
        IF (RAD < 0.D0) CYCLE
        RX = CX * SQRT(RAD)
        RY = CY * SQRT(RAD)

        CALL CTELL (X0,Y0,RX,RY,XLIMS1,YLIMS1,XLIMS2,YLIMS2,RLB,
     .              PHIAN,PHIEN,IPART)

        DO IPRT=1,IPART
          DPH = (PHIEN(IPRT)-PHIAN(IPRT)) / FLOAT(N-1)
          X = X0 + RX*COS(PHIAN(IPRT))
          Y = Y0 + RY*SIN(PHIAN(IPRT))
          CALL PL3D (X,Y,Z,XS,YS)
          CALL GRJMP (SNGL(XS),SNGL(YS))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            X = X0 + RX*COS(PHI)
            Y = Y0 + RY*SIN(PHI)
            CALL PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (SNGL(XS),SNGL(YS))
          END DO
        END DO

      END DO ! IZ

!  PLOT Y-ISOLINES

      DO IY = 1, NY
        Y = YA + (IY-1)*DY
        RAD = 1.D0 - (Y-Y0)**2/CY**2
        IF (RAD < 0.D0) CYCLE
        RX = CX * SQRT(RAD)
        RZ = CZ * SQRT(RAD)

        CALL CTELL (X0,Z0,RX,RZ,XLIMS1,ZLIMS1,XLIMS2,ZLIMS2,RLB,
     .              PHIAN,PHIEN,IPART)

        DO IPRT=1,IPART
          DPH = (PHIEN(IPRT)-PHIAN(IPRT)) / FLOAT(N-1)
          X = X0 + RX*COS(PHIAN(IPRT))
          Z = Z0 + RZ*SIN(PHIAN(IPRT))
          CALL PL3D (X,Y,Z,XS,YS)
          CALL GRJMP (SNGL(XS),SNGL(YS))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            X = X0 + RX*COS(PHI)
            Z = Z0 + RZ*SIN(PHI)
            CALL PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (SNGL(XS),SNGL(YS))
          END DO
        END DO
      END DO ! IY

!  PLOT X-ISOLINES

      DO IX = 1, NX
        X = XA + (IX-1)*DX
        RAD = 1.D0 - (X-X0)**2/CX**2
        IF (RAD < 0.D0) CYCLE
        RY = CY * SQRT(RAD)
        RZ = CZ * SQRT(RAD)

        CALL CTELL (Y0,Z0,RY,RZ,YLIMS1,ZLIMS1,YLIMS2,ZLIMS2,RLB,
     .              PHIAN,PHIEN,IPART)

        DO IPRT=1,IPART
          DPH = (PHIEN(IPRT)-PHIAN(IPRT)) / FLOAT(N-1)
          Y = Y0 + RY*COS(PHIAN(IPRT))
          Z = Z0 + RZ*SIN(PHIAN(IPRT))
          CALL PL3D (X,Y,Z,XS,YS)
          CALL GRJMP (SNGL(XS),SNGL(YS))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            Y = Y0 + RY*COS(PHI)
            Z = Z0 + RZ*SIN(PHI)
            CALL PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (SNGL(XS),SNGL(YS))
          END DO
        END DO
      END DO ! IY

      CALL GRNWPN(1)

      RETURN
      END

      SUBROUTINE CTELL (X0,Y0,RX,RY,XLIMS1,YLIMS1,XLIMS2,YLIMS2,RLB,
     .                  PHIAN,PHIEN,IPART)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'CCONA'
      DIMENSION PHIAN(*),PHIEN(*),PHIANG(10)

      IPHI = 0

!  DETERMINE INTERSECTION POINTS OF ELLIPSE WITH BOX

!  INTERSECTION WITH X=XLIMS1
      QX1 = (XLIMS1 - X0) / (RX+EPS30)
      IF (ABS(QX1) <= 1.D0) THEN
        PHI = ACOS(QX1)
        Y1 = Y0 + RY*SIN(PHI)
        IF ((Y1 >= YLIMS1) .AND. (Y1 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
        PHI = PI2A - PHI
        Y2 = Y0 + RY*SIN(PHI)
        IF ((Y2 >= YLIMS1) .AND. (Y2 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
      END IF

!  INTERSECTION WITH X=XLIMS2
      QX2 = (XLIMS2 - X0) / (RX+EPS30)
      IF (ABS(QX2) <= 1.D0) THEN
        PHI = ACOS(QX2)
        Y1 = Y0 + RY*SIN(PHI)
        IF ((Y1 >= YLIMS1) .AND. (Y1 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
        PHI = PI2A - PHI
        Y2 = Y0 + RY*SIN(PHI)
        IF ((Y2 >= YLIMS1) .AND. (Y2 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
      END IF

!  INTERSECTION WITH Y=YLIMS1
      QY1 = (YLIMS1 - Y0) / (RY+EPS30)
      IF (ABS(QY1) <= 1.D0) THEN
        PHI1 = ASIN(QY1)
        PHI2 = PIA - PHI1
        IF (PHI1 < 0.D0) PHI1 = PHI1 + PI2A
        X1 = X0 + RX*COS(PHI1)
        IF ((X1 >= XLIMS1) .AND. (X1 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI1
        END IF
        X2 = X0 + RX*COS(PHI2)
        IF ((X2 >= XLIMS1) .AND. (X2 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI2
        END IF
      END IF

!  INTERSECTION WITH Y=YLIMS2
      QY1 = (YLIMS2 - Y0) / (RY+EPS30)
      IF (ABS(QY1) <= 1.D0) THEN
        PHI1 = ASIN(QY1)
        PHI2 = PIA - PHI1
        IF (PHI1 < 0.D0) PHI1 = PHI1 + PI2A
        X1 = X0 + RX*COS(PHI1)
        IF ((X1 >= XLIMS1) .AND. (X1 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI1
        END IF
        X2 = X0 + RX*COS(PHI2)
        IF ((X2 >= XLIMS1) .AND. (X2 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI2
        END IF
      END IF

      IF (IPHI == 0) THEN
        IPART = 1
        PHIAN(1) = 0.D0
        PHIEN(1) = PI2A
        RETURN
      END IF
C
C  SORT PHIANG
      ISORT = 1
      DO WHILE (ISORT > 0)
        ISORT=0
        DO I=1,IPHI-1
          IF (PHIANG(I+1).LT.PHIANG(I)) THEN
            PHIH=PHIANG(I)
            PHIANG(I)=PHIANG(I+1)
            PHIANG(I+1)=PHIH
            ISORT=ISORT+1
          ENDIF
        END DO
      END DO
C
      IPHI=IPHI+1
      PHIANG(IPHI)=PHIANG(1)
C
      IPART=0
      DO I=1,IPHI-1
        IF (ABS(PHIANG(I+1)-PHIANG(I)).LT.EPS10) CYCLE
        IF (PHIANG(I+1).LT.PHIANG(I)) PHIANG(I+1)=PHIANG(I+1)+PI2A
        PHI=0.5*(PHIANG(I)+PHIANG(I+1))
        XPHI=X0+RX*COS(PHI)
        YPHI=Y0+RY*SIN(PHI)
        IF (RLB.LT.1.5) THEN
          IF (XPHI.GE.XLIMS1.AND.XPHI.LE.XLIMS2.AND.
     .        YPHI.GE.YLIMS1.AND.YPHI.LE.YLIMS2) THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ELSE
          IF (XPHI.LE.XLIMS1.OR.XPHI.GE.XLIMS2.OR.
     .        YPHI.LE.YLIMS1.OR.YPHI.GE.YLIMS2)THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ENDIF
      END DO
C
      RETURN
      END
