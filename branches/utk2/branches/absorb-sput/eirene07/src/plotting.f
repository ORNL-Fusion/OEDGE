C EIRENE07 COMPILATION
C ===== SOURCE: anpsg.f
C
C
C----------------------------------------------------------------------*
C SUBROUTINE ANPSG                                                     *
C----------------------------------------------------------------------*
      SUBROUTINE ANPSG(MIN,MAX,INTNR,STPSZ,IERR)
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
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      REAL(SP), INTENT(INOUT) :: MAX, MIN
      REAL(SP), INTENT(OUT) :: STPSZ
      INTEGER, INTENT(INOUT) :: IERR
      INTEGER, INTENT(OUT) :: INTNR
      REAL(DP) :: LKS, RTS, MDDL, STPS2, DIFF
      INTEGER :: IPT, IEXP10
C
      DIFF=MAX-MIN
      IF (DIFF.EQ.0) THEN
         IPT=IEXP10(REAL(MAX,KIND=DP))
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
         MIN=LKS+SIGN(1._DP,LKS)*10.D0**(-14+IPT)
         IPT=IEXP10(RTS)
         MAX=RTS+SIGN(1._DP,RTS)*10.D0**(-14+IPT)
         IPT=IEXP10(STPS2)
         STPSZ=STPS2+10.**(-14+IPT)
      ELSE
         WRITE(iunout,*)  '------------------------------------'
         WRITE(iunout,*)  'PARAMETERFEHLER IN ANPSG: MIN > MAX'
         WRITE(iunout,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
C
      RETURN
      END
C ===== SOURCE: anpsgl.f
C
C
C----------------------------------------------------------------------*
C SUBROUTINE ANPSGL                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE ANPSGL(MIN,MAX,MINL,MAXL,IERR)
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
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(SP), INTENT(IN) ::  MIN, MAX
      INTEGER, INTENT(INOUT) :: IERR
      INTEGER, INTENT(OUT) :: MINL, MAXL
      INTEGER :: IEXP10
C
      IF (MIN.LE.0.OR.MAX.LE.0) THEN
         WRITE(iunout,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(iunout,*)  '***    MIN<=0  ODER  MAX<=0    ***'
         WRITE(iunout,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MIN.GT.MAX) THEN
         WRITE(iunout,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(iunout,*)  '***        MIN  >  MAX         ***'
         WRITE(iunout,*)  'MIN= ',MIN,' MAX= ',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MAX-MIN.GE.0) THEN
         MINL=IEXP10(REAL(MIN,KIND=DP))
         MAXL=IEXP10(REAL(MAX,KIND=DP))+1
      ENDIF
C
      RETURN
      END
C ===== SOURCE: celint.f
C
!pb  19.12.06:  initialise YWERT
C
      SUBROUTINE CELINT (AORIG,YWERT,LOGL,IBLD,ICURV,N1DIM,IERR)
C
C  INTERPOLATION FROM CELL CENTERS (AORIG)
C  TO QUANTITIES AT CELL VERTICES
C  INPUT : AORIG
C  OUTPUT: YWERT
C
C  IN CASE LOGL=.TRUE.: YWERT=LOG10(YWERT)
C
C  THE FOLLOWING GRID OPTIONS ARE AVAILBLE CURRENTLY:
C
C  IN CASE (LEVGEO=1 OR LEVGEO=2).AND.LPPOL3, THE RHOSRF AND ZSURF GRIDS
C                                             ARE USED.
C  IN CASE (LEVGEO=1            ).AND.LPTOR3, THE RHOSRF AND PSURF GRIDS
C                                             ARE USED.
C  IN CASE (LEVGEO=1            ).AND.LPRAD3, TO BE WRITTEN: PSURF AND ZSURF
C                                                
C  IN CASE (LEVGEO=3 OR LEVGEO=2).AND.LPTOR3, THE FULL POLYGON GRIDS
C                                             ARE USED.
C  IN CASE LEVGEO=4,             .AND.LPTOR3  THE FULL TRIANGULAR MESH
C                                             IS USED
C  IN CASE LEVGEO=5,             .AND.LRPSCUT A FULL TRIANGULAR MESH
C                                             IS USED AS OPTAINED FROM RPSCUT
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
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CLGIN
      USE CTRIG
      USE CTETRA

      IMPLICIT NONE
C
      TYPE(CELL_ELEM), POINTER :: CUR

      INTEGER, INTENT(IN) :: IBLD, ICURV, N1DIM
      INTEGER, INTENT(OUT) :: IERR
      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(OUT) :: YWERT(N1DIM,*)
      LOGICAL, INTENT(IN) :: LOGL

      REAL(DP) :: TEILA(4), TEILWERT(4), VOLSUM(NCOORD)
      REAL(DP) :: IX, IY, JX, JY, KX, KY, LX, LY, MX, MY
      REAL(DP) :: FLAECH, AX, AY, BX, BY, CX, CY, DX, DY, EX, EY, FX,
     .            FY, GX, GY, HX, HY, GESA, AGES, DIST
      REAL(DP) :: dxcom(nr1st-1), dycom(np2nd-1),
     .            dzcom(nt3rd-1), vorig(nr1st-1, np2nd-1, nt3rd-1),
     .            value(nr1st, np2nd, nt3rd),
     .            dist1, dist2, dist3, dist4, dist5, dist6, dist7,
     .            dist8, summedist
      REAL(DP), ALLOCATABLE, SAVE :: XSTGRD(:)
      INTEGER :: ZUORD(NKNOT,0:50)
      INTEGER :: AGEF, EGEF
      INTEGER :: AGLEICH(4), EGLEICH(4)
      INTEGER :: JPART, JP, K, JTEST, IRG, I, J, IT, IR, IP, IPA, IPE,
     .           IRA, IRE, IRD, IPART, IRMEGM, IRMAG, IREGM, IRAG,
     .           IRMIP, IRIP, IRMIPM, IRIPM, IC, IN, IRM1, IPM1, ITM1,
     .           IN1, IN2, IN3, IN4, IN5, IN6, IN7, IN8
C

      IF ((LEVGEO <= 3) .AND. .NOT.ALLOCATED(XSTGRD)) THEN
        ALLOCATE (XSTGRD(NRAD))
        WHERE (NSTGRD > 0)
          XSTGRD = 0._DP
        ELSEWHERE
          XSTGRD = 1._DP
        END WHERE
      END IF

      IERR=0
C   X-Z PLOT ON Y=CONST PLANE
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.LPPOL3(IBLD)) THEN
        IP=1
        IF (NLPOL) IP=IPROJ3(IBLD,ICURV)
        IF (IP.LE.0.OR.IP.GT.NP2ND) IP=1
        YWERT(1:N1ST,1:N2ND+N3RD) = 0._DP
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
              TEILWERT(1) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(1) = ABS(AX-FX)*ABS(AY-GY)*XSTGRD(IRD)
            ENDIF
C  UNTEN LINKS
            IF ((IR .NE. NR1ST) .AND. (IT .NE. 1)) THEN
C  AKTUELLER PUNKT LIEGT NICHT AUF DEM LETZTEN
C  POLYGON UND IST NICHT ANFANGSPUNKT EINES POLYGONS
              AX = RHOSRF(IR)
              AY = ZSURF(IT)
              LX = RHOSRF(IR+1)
              KY = ZSURF(IT-1)
              GY = 0.5 * (AY + KY)
              HX = 0.5 * (AX + LX)
              IRD=IR+((IP-1)+(IT-2)*NP2T3)*NR1P2
              TEILWERT(2) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(2) = ABS(AX-HX)*ABS(AY-GY)*XSTGRD(IRD)
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
              TEILWERT(3) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(3) = ABS(AX-HX)*ABS(AY-IY)*XSTGRD(IRD)
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
              TEILWERT(4) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(4) = ABS(AX-FX)*ABS(AY-IY)*XSTGRD(IRD)
            ENDIF
C
            AGES = TEILA(1) + TEILA(2) + TEILA(3) + TEILA(4) + EPS60
            YWERT(IR,IT) = 0.
            DO 3103,J=1,4
              YWERT(IR,IT) = YWERT(IR,IT) + TEILA(J)/AGES*TEILWERT(J)
3103        CONTINUE
            IF (LOGL) YWERT(IR,IT)=LOG10(MAX(1.E-48_DP,YWERT(IR,IT)))
3100      CONTINUE
1100    CONTINUE
C
C   X-Y PLOT ON Z=CONST PLANE
      ELSEIF (LEVGEO.EQ.1.AND.LPTOR3(IBLD)) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        YWERT(1:N1ST,1:N2ND+N3RD) = 0._DP
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
              TEILWERT(1) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(1) = ABS(AX-FX)*ABS(AY-GY)*XSTGRD(IRD)
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
              TEILWERT(2) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(2) = ABS(AX-HX)*ABS(AY-GY)*XSTGRD(IRD)
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
              TEILWERT(3) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(3) = ABS(AX-HX)*ABS(AY-IY)*XSTGRD(IRD)
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
              TEILWERT(4) = AORIG(IRD)*XSTGRD(IRD)
              TEILA(4) = ABS(AX-FX)*ABS(AY-IY)*XSTGRD(IRD)
            ENDIF
C
            AGES = TEILA(1) + TEILA(2) + TEILA(3) + TEILA(4) + EPS60
            YWERT(IR,IP) = 0.
            DO 3113,J=1,4
              YWERT(IR,IP) = YWERT(IR,IP) + TEILA(J)/AGES*TEILWERT(J)
3113        CONTINUE
            IF (LOGL) YWERT(IR,IP)=LOG10(MAX(1.E-48_DP,YWERT(IR,IP)))
3110      CONTINUE
1110    CONTINUE
C
C  X-Y-Z PLOT (CUBE)
      ELSEIF (LEVGEO.EQ.1.AND.NLTRZ
     .   .AND.NLRAD.AND.NLPOL.AND.NLTOR
     .   .AND..NOT.(LPPOL3(IBLD).OR.LPTOR3(IBLD).OR.LPRAD3(IBLD))) THEN
         do ir=1,nr1st-1
            dxcom(ir) = ((rsurf(ir)-rsurf(ir+1))/2.)**2
         enddo
         do ip=1,np2nd-1
            dycom(ip) = ((psurf(ip)-psurf(ip+1))/2.)**2
         enddo
         do it=1,nt3rd-1
            dzcom(it) = ((zsurf(it)-zsurf(it+1))/2.)**2
         enddo
         do ir=1,nr1st-1
            do ip=1,np2nd-1
              do it=1,nt3rd-1
                 IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                 vorig(ir,ip,it) = aorig(ird)
              enddo
            enddo
         enddo
c berechnung im inneren des quaders
         do ir=2,nr1st-1
            do ip=2,np2nd-1
              do it=2,nt3rd-1
                 irm1=ir-1
                 ipm1=ip-1
                 itm1=it-1
                 in1=irm1+((IPM1-1)+(ITM1-1)*NP2T3)*NR1P2
                 in2=irm1+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
                 in3=irm1+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
                 in4=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                 in5=ir+((IPM1-1)+(ITM1-1)*NP2T3)*NR1P2
                 in6=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
                 in7=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
                 in8=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
                 dist1 = sqrt(dxcom(ir-1)+dycom(ip-1)+dzcom(it-1))*
     .                   xstgrd(in1)
                 dist2 = sqrt(dxcom(ir-1)+dycom(ip-1)+dzcom(it))*
     .                   xstgrd(in2)
                 dist3 = sqrt(dxcom(ir-1)+dycom(ip)+dzcom(it-1))*
     .                   xstgrd(in3)
                 dist4 = sqrt(dxcom(ir-1)+dycom(ip)+dzcom(it))*
     .                   xstgrd(in4)
                 dist5 = sqrt(dxcom(ir)+dycom(ip-1)+dzcom(it-1))*
     .                   xstgrd(in5)
                 dist6 = sqrt(dxcom(ir)+dycom(ip-1)+dzcom(it))*
     .                   xstgrd(in6)
                 dist7 = sqrt(dxcom(ir)+dycom(ip)+dzcom(it-1))*
     .                   xstgrd(in7)
                 dist8 = sqrt(dxcom(ir)+dycom(ip)+dzcom(it))*
     .                   xstgrd(in8)
                 summedist = dist1+dist2+dist3+dist4+dist5+dist6+
     .                       dist7+dist8+eps60
                 value(ir,ip,it) = (dist1*vorig(ir-1,ip-1,it-1) +
     .                              dist2*vorig(ir-1,ip-1,it)   +
     .                              dist3*vorig(ir-1,ip,it-1)   +
     .                              dist4*vorig(ir-1,ip,it)     +
     .                              dist5*vorig(ir,ip-1,it-1)   +
     .                              dist6*vorig(ir,ip-1,it)     +
     .                              dist7*vorig(ir,ip,it-1)     +
     .                              dist8*vorig(ir,ip,it))/summedist
              enddo
            enddo
         enddo
c berechnung der seiten
         do ir=2,nr1st-1
            do ip=2,np2nd-1
               irm1=ir-1
               ipm1=ip-1
               it=1
               in1=irm1+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
               in2=irm1+((IP-1)  +(IT-1)*NP2T3)*NR1P2
               in3=ir  +((IPM1-1)+(IT-1)*NP2T3)*NR1P2
               in4=ir  +((IP-1)  +(IT-1)*NP2T3)*NR1P2
               dist1 = sqrt(dxcom(ir-1)+dycom(ip-1)+dzcom(1))*
     .                 xstgrd(in1)
               dist2 = sqrt(dxcom(ir-1)+dycom(ip)+dzcom(1))*
     .                 xstgrd(in2)
               dist3 = sqrt(dxcom(ir)+dycom(ip-1)+dzcom(1))*
     .                 xstgrd(in3)
               dist4 = sqrt(dxcom(ir)+dycom(ip)+dzcom(1))*
     .                 xstgrd(in4)
               summedist = dist1+dist2+dist3+dist4+eps60
               value(ir,ip,1) = (dist1*vorig(ir-1,ip-1,1) +
     .              dist2*vorig(ir-1,ip,1)   +
     .              dist3*vorig(ir,ip-1,1)   +
     .              dist4*vorig(ir,ip,1))/summedist

               irm1=ir-1
               ipm1=ip-1
               it=nt3rd-1
               in1=irm1+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
               in2=irm1+((IP-1)  +(IT-1)*NP2T3)*NR1P2
               in3=ir  +((IPM1-1)+(IT-1)*NP2T3)*NR1P2
               in4=ir  +((IP-1)  +(IT-1)*NP2T3)*NR1P2
               dist1 = sqrt(dxcom(ir-1)+dycom(ip-1)+dzcom(nt3rd-1))*
     .                 xstgrd(in1)
               dist2 = sqrt(dxcom(ir-1)+dycom(ip)+dzcom(nt3rd-1))*
     .                 xstgrd(in2)
               dist3 = sqrt(dxcom(ir)+dycom(ip-1)+dzcom(nt3rd-1))*
     .                 xstgrd(in3)
               dist4 = sqrt(dxcom(ir)+dycom(ip)+dzcom(nt3rd-1))*
     .                 xstgrd(in4)
               summedist = dist1+dist2+dist3+dist4+eps60
               value(ir,ip,nt3rd) = (dist1*vorig(ir-1,ip-1,nt3rd-1) +
     .              dist2*vorig(ir-1,ip,nt3rd-1)   +
     .              dist3*vorig(ir,ip-1,nt3rd-1)   +
     .              dist4*vorig(ir,ip,nt3rd-1))/summedist
            enddo
         enddo
         do ir=2,nr1st-1
            do it=2,nt3rd-1
               irm1=ir-1
               ip=1
               itm1=it-1
               in1=irm1+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
               in2=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
               in3=ir  +((IP-1)+(ITM1-1)*NP2T3)*NR1P2
               in4=ir  +((IP-1)+(IT-1)*NP2T3)*NR1P2
               dist1 = sqrt(dxcom(ir-1)+dycom(1)+dzcom(it-1))*
     .                 xstgrd(in1)
               dist2 = sqrt(dxcom(ir-1)+dycom(1)+dzcom(it))*
     .                 xstgrd(in2)
               dist3 = sqrt(dxcom(ir)+dycom(1)+dzcom(it-1))*
     .                 xstgrd(in3)
               dist4 = sqrt(dxcom(ir)+dycom(1)+dzcom(it))*
     .                 xstgrd(in4)
               summedist = dist1+dist2+dist3+dist4+eps60
               value(ir,1,it) = (dist1*vorig(ir-1,1,it-1) +
     .              dist2*vorig(ir-1,1,it)   +
     .              dist3*vorig(ir,1,it-1)   +
     .              dist4*vorig(ir,1,it))/summedist

               irm1=ir-1
               ip=np2nd-1
               itm1=it-1
               in1=irm1+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
               in2=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
               in3=ir  +((IP-1)+(ITM1-1)*NP2T3)*NR1P2
               in4=ir  +((IP-1)+(IT-1)*NP2T3)*NR1P2
               dist1 = sqrt(dxcom(ir-1)+dycom(np2nd-1)+dzcom(it-1))*
     .                 xstgrd(in1)
               dist2 = sqrt(dxcom(ir-1)+dycom(np2nd-1)+dzcom(it))*
     .                 xstgrd(in2)
               dist3 = sqrt(dxcom(ir)+dycom(np2nd-1)+dzcom(it-1))*
     .                 xstgrd(in3)
               dist4 = sqrt(dxcom(ir)+dycom(np2nd-1)+dzcom(it))*
     .                 xstgrd(in4)
               summedist = dist1+dist2+dist3+dist4+eps60
               value(ir,np2nd,it) = (dist1*vorig(ir-1,np2nd-1,it-1) +
     .              dist2*vorig(ir-1,np2nd-1,it)   +
     .              dist3*vorig(ir,np2nd-1,it-1)   +
     .              dist4*vorig(ir,np2nd-1,it))/summedist
            enddo
         enddo
         do ip=2,np2nd-1
            do it=2,nt3rd-1
               ir=1
               ipm1=ip-1
               itm1=it-1
               in1=ir+((IPM1-1)+(ITM1-1)*NP2T3)*NR1P2
               in2=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
               in3=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
               in4=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
               dist1 = sqrt(dxcom(1)+dycom(ip-1)+dzcom(it-1))*
     .                 xstgrd(in1)
               dist2 = sqrt(dxcom(1)+dycom(ip-1)+dzcom(it))*
     .                 xstgrd(in2)
               dist3 = sqrt(dxcom(1)+dycom(ip)+dzcom(it-1))*
     .                 xstgrd(in3)
               dist4 = sqrt(dxcom(1)+dycom(ip)+dzcom(it))*
     .                 xstgrd(in4)
               summedist = dist1+dist2+dist3+dist4+eps60
               value(1,ip,it) = (dist1*vorig(1,ip-1,it-1) +
     .              dist2*vorig(1,ip-1,it)   +
     .              dist3*vorig(1,ip,it-1)   +
     .              dist4*vorig(1,ip,it))/summedist

               ir=nr1st-1
               ipm1=ip-1
               itm1=it-1
               in1=ir+((IPM1-1)+(ITM1-1)*NP2T3)*NR1P2
               in2=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
               in3=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
               in4=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
               dist1 = sqrt(dxcom(nr1st-1)+dycom(ip-1)+dzcom(it-1))*
     .                 xstgrd(in1)
               dist2 = sqrt(dxcom(nr1st-1)+dycom(ip-1)+dzcom(it))*
     .                 xstgrd(in2)
               dist3 = sqrt(dxcom(nr1st-1)+dycom(ip)+dzcom(it-1))*
     .                 xstgrd(in3)
               dist4 = sqrt(dxcom(nr1st-1)+dycom(ip)+dzcom(it))*
     .                 xstgrd(in4)
               summedist = dist1+dist2+dist3+dist4+eps60
               value(nr1st,ip,it) = (dist1*vorig(nr1st-1,ip-1,it-1) +
     .              dist2*vorig(nr1st-1,ip-1,it)   +
     .              dist3*vorig(nr1st-1,ip,it-1)   +
     .              dist4*vorig(nr1st-1,ip,it))/summedist
            enddo
         enddo
c berechnung der kanten
         do it=2,nt3rd-1
            ir=1
            ip=1
            itm1=it-1
            in1=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(1)+dycom(1)+dzcom(it-1))*xstgrd(in1)
            dist2= sqrt(dxcom(1)+dycom(1)+dzcom(it))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(1,1,it) = (dist1*vorig(1,1,it-1) +
     .           dist2*vorig(1,1,it))/summedist

            ir=nr1st-1
            ip=1
            itm1=it-1
            in1=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(nr1st-1)+dycom(1)+dzcom(it-1))*xstgrd(in1)
            dist2= sqrt(dxcom(nr1st-1)+dycom(1)+dzcom(it))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(nr1st,1,it) = (dist1*vorig(nr1st-1,1,it-1) +
     .           dist2*vorig(nr1st-1,1,it))/summedist

            ir=1
            ip=np2nd-1
            itm1=it-1
            in1=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(1)+dycom(np2nd-1)+dzcom(it-1))*xstgrd(in1)
            dist2= sqrt(dxcom(1)+dycom(np2nd-1)+dzcom(it))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(1,np2nd,it) = (dist1*vorig(1,np2nd-1,it-1) +
     .           dist2*vorig(1,np2nd-1,it))/summedist

            ir=nr1st-1
            ip=np2nd-1
            itm1=it-1
            in1=ir+((IP-1)+(ITM1-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1 = sqrt(dxcom(nr1st-1)+dycom(np2nd-1)+dzcom(it-1))*
     .              xstgrd(in1)
            dist2 = sqrt(dxcom(nr1st-1)+dycom(np2nd-1)+dzcom(it))*
     .              xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(nr1st,np2nd,it) = (dist1*vorig(nr1st-1,np2nd-1,it-1) +
     .           dist2*vorig(nr1st-1,np2nd-1,it))/summedist

         enddo

         do ir=2,nr1st-1
            irm1=ir-1
            ip=1
            it=1
            in1=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir  +((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1 = sqrt(dxcom(ir-1)+dycom(1)+dzcom(1))*xstgrd(in1)
            dist2 = sqrt(dxcom(ir)+dycom(1)+dzcom(1))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(ir,1,1) = (dist1*vorig(ir-1,1,1) +
     .           dist2*vorig(ir,1,1))/summedist

            irm1=ir-1
            ip=1
            it=nt3rd-1
            in1=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir  +((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(ir-1)+dycom(1)+dzcom(nt3rd-1))*xstgrd(in1)
            dist2= sqrt(dxcom(ir)+dycom(1)+dzcom(nt3rd-1))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(ir,1,nt3rd) = (dist1*vorig(ir-1,1,nt3rd-1) +
     .           dist2*vorig(ir,1,nt3rd-1))/summedist

            irm1=ir-1
            ip=np2nd-1
            it=1
            in1=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir  +((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(ir-1)+dycom(np2nd-1)+dzcom(1))*xstgrd(in1)
            dist2= sqrt(dxcom(ir)+dycom(np2nd-1)+dzcom(1))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(ir,np2nd,1) = (dist1*vorig(ir-1,np2nd-1,1) +
     .           dist2*vorig(ir,np2nd-1,1))/summedist

            irm1=ir-1
            ip=np2nd-1
            it=nt3rd-1
            in1=irm1+((IP-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir  +((IP-1)+(IT-1)*NP2T3)*NR1P2
            dist1 = sqrt(dxcom(ir-1)+dycom(np2nd-1)+dzcom(nt3rd-1))*
     .              xstgrd(in1)
            dist2 = sqrt(dxcom(ir)+dycom(np2nd-1)+dzcom(nt3rd-1))*
     .              xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(ir,np2nd,nt3rd) = (dist1*vorig(ir-1,np2nd-1,nt3rd-1) +
     .           dist2*vorig(ir,np2nd-1,nt3rd-1))/summedist

         enddo

         do ip=2,np2nd-1
            ir=1
            ipm1=ip-1
            it=1
            in1=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)  +(IT-1)*NP2T3)*NR1P2
            dist1 = sqrt(dxcom(1)+dycom(ip-1)+dzcom(1))*xstgrd(in1)
            dist2 = sqrt(dxcom(1)+dycom(ip)+dzcom(1))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(1,ip,1) = (dist1*vorig(1,ip-1,1) +
     .           dist2*vorig(1,ip,1))/summedist

            ir=1
            ipm1=ip-1
            it=nt3rd-1
            in1=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)  +(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(1)+dycom(ip-1)+dzcom(nt3rd-1))*xstgrd(in1)
            dist2= sqrt(dxcom(1)+dycom(ip)+dzcom(nt3rd-1))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(1,ip,nt3rd) = (dist1*vorig(1,ip-1,nt3rd-1) +
     .           dist2*vorig(1,ip,nt3rd-1))/summedist

            ir=nr1st-1
            ipm1=ip-1
            it=1
            in1=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)  +(IT-1)*NP2T3)*NR1P2
            dist1= sqrt(dxcom(nr1st-1)+dycom(ip-1)+dzcom(1))*xstgrd(in1)
            dist2= sqrt(dxcom(nr1st-1)+dycom(ip)+dzcom(1))*xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(nr1st,ip,1) = (dist1*vorig(nr1st-1,ip-1,1) +
     .           dist2*vorig(nr1st-1,ip,1))/summedist

            ir=nr1st-1
            ipm1=ip-1
            it=nt3rd-1
            in1=ir+((IPM1-1)+(IT-1)*NP2T3)*NR1P2
            in2=ir+((IP-1)  +(IT-1)*NP2T3)*NR1P2
            dist1 = sqrt(dxcom(nr1st-1)+dycom(ip-1)+dzcom(nt3rd-1))*
     .              xstgrd(in1)
            dist2 = sqrt(dxcom(nr1st-1)+dycom(ip)+dzcom(nt3rd-1))*
     .              xstgrd(in2)
            summedist = dist1+dist2+eps60
            value(nr1st,ip,nt3rd) = (dist1*vorig(nr1st-1,ip-1,nt3rd-1) +
     .           dist2*vorig(nr1st-1,ip,nt3rd-1))/summedist

         enddo
c berechnung der ecken
         in1=1+((1-1)+(1-1)*NP2T3)*NR1P2
         value(1,1,1) = vorig(1,1,1)*xstgrd(in1)
         in2=1+((1-1)+(nt3rd-1)*NP2T3)*NR1P2
         value(1,1,nt3rd) = vorig(1,1,nt3rd-1)*xstgrd(in2)
         in3=1+((np2nd-1)+(1-1)*NP2T3)*NR1P2
         value(1,np2nd,1) = vorig(1,np2nd-1,1)*xstgrd(in3)
         in4=1+((np2nd-1)+(nt3rd-1)*NP2T3)*NR1P2
         value(1,np2nd,nt3rd) = vorig(1,np2nd-1,nt3rd-1)*xstgrd(in4)
         in5=nr1st+((1-1)+(1-1)*NP2T3)*NR1P2
         value(nr1st,1,1) = vorig(nr1st-1,1,1)*xstgrd(in5)
         in6=nr1st+((1-1)+(nt3rd-1)*NP2T3)*NR1P2
         value(nr1st,1,nt3rd) = vorig(nr1st-1,1,nt3rd-1)*xstgrd(in6)
         in7=nr1st+((np2nd-1)+(1-1)*NP2T3)*NR1P2
         value(nr1st,np2nd,1) = vorig(nr1st-1,np2nd-1,1)*xstgrd(in7)
         in8=nr1st+((np2nd-1)+(nt3rd-1)*NP2T3)*NR1P2
         value(nr1st,np2nd,nt3rd) = vorig(nr1st-1,np2nd-1,nt3rd-1)*
     .                              xstgrd(in8)
c  fertig, zurueck auf ywert
         YWERT(1:NRAD,1) = 0._DP
         do ir=1,nr1st
            do ip=1,np2nd
               do it=1,nt3rd
                  IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  ywert(ird,1) = value(ir,ip,it)
               enddo
            enddo
         enddo

      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTOR3(IBLD)) THEN
C
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        YWERT(1:N1ST,1:N2ND+N3RD) = 0._DP
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
              IC=INDPOINT(IR,IP)
              AGES = 0._DP
              YWERT(IR,IP) = 0.
              CUR => COORCELL(IC)%PCELL
              DO WHILE (ASSOCIATED(CUR))
                IN=CUR%NOCELL
                DIST1=1._DP/SQRT((XCOM(IN)-XPOL(IR,IP))**2+
     .                           (YCOM(IN)-YPOL(IR,IP))**2)
                AGES = AGES + DIST1*XSTGRD(IN)
                YWERT(IR,IP) = YWERT(IR,IP) + AORIG(IN)*DIST1
                CUR => CUR%NEXT_CELL
              END DO
              YWERT(IR,IP) = YWERT(IR,IP) / (AGES+EPS60)
              IF (LOGL) YWERT(IR,IP)=LOG10(MAX(1.E-48_DP,YWERT(IR,IP)))
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
      ELSEIF ((LEVGEO.EQ.4.AND.LPTOR3(IBLD)) .OR.
     .        (LEVGEO.EQ.5.AND.LRPSCUT)) THEN
C
        YWERT(1:NRAD,1) = 0._DP
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
            DIST=1._DP/SQRT((XTRIAN(I)-XCOM(K))**2+
     .                      (YTRIAN(I)-YCOM(K))**2)
            GESA = GESA + DIST
            YWERT(I,1) = YWERT(I,1) + DIST * AORIG(K)
70        CONTINUE
          IF (GESA .NE. 0) THEN
            YWERT(I,1) = YWERT(I,1)/GESA
          ELSE
            WRITE (iunout,*) 'ERROR IN CELINT: POINT ',I,' NOT IN MESH'
            YWERT(I,1) = YWERT(I,1)
          ENDIF
          IF (LOGL) YWERT(I,1)=LOG10(MAX(1.E-48_DP,YWERT(I,1)))
60      CONTINUE
C
C
      ELSEIF (LEVGEO.EQ.5.AND..NOT.LRPSCUT) THEN
         VOLSUM = 0.
         YWERT(1:NRAD,1) = 0.
         DO I=1,NTET
              DO J=1,4
                 DIST1=SQRT((XTETRA(NTECK(J,I))-XTCEN(I))**2 +
     .                      (YTETRA(NTECK(J,I))-YTCEN(I))**2 +
     .                      (ZTETRA(NTECK(J,I))-ZTCEN(I))**2)
                 VOLSUM(NTECK(J,I))=VOLSUM(NTECK(J,I))+1._DP/DIST1
                 YWERT(NTECK(J,I),1) = YWERT(NTECK(J,I),1)+
     .                                 AORIG(I)/DIST1
            ENDDO
         ENDDO
         WHERE (ABS(VOLSUM(1:NCOOR)) > 1.D-10)
           YWERT(1:NCOOR,1) = YWERT(1:NCOOR,1)/VOLSUM(1:NCOOR)
         END WHERE
         IF (LOGL)
     .    YWERT(1:NCOOR,1)=LOG10(MAX(1.E-48_DP,YWERT(1:NCOOR,1)))
      ELSE
        WRITE (iunout,*) 'CELINT CALLED WITH INVALID OPTIONS '
        WRITE (iunout,*) 'PLOT ABANDONNED '
        IERR=1
      ENDIF
C
      RETURN
      END
C ===== SOURCE: cone.f
C
C
      SUBROUTINE CONE (X0,Y0,Z0,VX,VY,VZ,T1,T2,ALF,NK,NP,NA,IO,NF,NUM,
     .                 ILEFT,AL,IRIGHT,AR)
C
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
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X0, Y0, Z0, VX, VY, VZ, T1, T2
      REAL(DP), INTENT(IN) :: AL(10), AR(10)
      INTEGER, INTENT(IN) :: NK, NP, NA, IO, NUM, ILEFT, IRIGHT
      LOGICAL, INTENT(IN) :: NF
      REAL(DP), ALLOCATABLE :: XP(:), YP(:)
      REAL(DP) :: T, TH, PX, PY, PZ, DANG, CX, CY, CZ, BA, BB, BC, PXX,
     .          PYY, PZZ, PHI, RAD, XK, YK, ALF, AX, AY, AZ, BX, BY, BZ,
     .          PI, DT
      REAL(SP), ALLOCATABLE :: XPS(:), YPS(:)
      INTEGER :: I, IA, JJ, NAPK, IE, J
C
      NAPK = MAX(NA,NP,NK) + 1
      ALLOCATE (XP(NAPK))
      ALLOCATE (YP(NAPK))
      ALLOCATE (XPS(NAPK))
      ALLOCATE (YPS(NAPK))
C
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
        WRITE (iunout,*) 'FEHLER IN DER EINGABE VON VX,VY,VZ ',
     .                    NUM,VX,VY,VZ
        CALL EXIT_OWN(1)
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

      DEALLOCATE (XP)
      DEALLOCATE (YP)
      DEALLOCATE (XPS)
      DEALLOCATE (YPS)

      RETURN
      END
C ===== SOURCE: ctcirc.f
C
C
      SUBROUTINE CTCIRC (A0,A1,A2,A3,A4,XP,YL1,YL2,XM,YM,
     .                   PHIANG,IPHI,ISW)
C
C   SOLVE  A0+A1*X+A2*Y+A3*X**2+A4*Y**2=0 FOR X=XP
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A0, A1, A2, A3, A4, XP, YL1, YL2, XM, YM
      REAL(DP), INTENT(OUT) :: PHIANG(*)
      INTEGER, INTENT(IN) :: ISW
      INTEGER, INTENT(INOUT) :: IPHI
      REAL(DP) :: PI, RD, PHI1, PHI2, YP1, YP2

      PI=4.*ATAN(1.)
      RD=0.25*A2**2-(A0+A1*XP+A3*XP*XP)*A4
      IF (RD.LT.0.) THEN
        WRITE (iunout,*) ' WURZEL NEGATIV IN CTCIRC '
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
C ===== SOURCE: ctell.f


      SUBROUTINE CTELL (X0,Y0,RX,RY,XLIMS1,YLIMS1,XLIMS2,YLIMS2,RLB,
     .                  PHIAN,PHIEN,IPART)

      USE PRECISION
      USE CCONA

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X0, Y0, RX, RY, XLIMS1, YLIMS1, XLIMS2,
     .                      YLIMS2, RLB
      REAL(DP), INTENT(OUT) :: PHIAN(*), PHIEN(*)
      INTEGER, INTENT(OUT) :: IPART
      REAL(DP) :: PHIANG(10)
      REAL(DP) :: PHI1, PHI2, X1, X2, Y1, Y2, XPHI, YPHI, QY1, QX2,
     .            PHIH, PHI, QX1
      INTEGER :: IPHI, ISORT, I

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
C ===== SOURCE: ctqua.f
C
C
      SUBROUTINE CTQUA (A0LM,A1LM,A2LM,A3LM,A4LM,A5LM,A6LM,A7LM,A8LM,
     .                  A9LM,XLIMS1,XLIMS2,YLIMS1,YLIMS2,ZLIMS1,ZLIMS2,
     .                  RLB,RAD,CX,CY,CZ,PHIAN,PHIEN,IPART)

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: A0LM, A1LM, A2LM, A3LM, A4LM, A5LM, A6LM,
     .                      A7LM, A8LM, A9LM, XLIMS1, XLIMS2, YLIMS1,
     .                      YLIMS2, ZLIMS1, ZLIMS2, RLB, RAD, CX, CY, CZ
      REAL(DP), INTENT(OUT) :: PHIAN(*), PHIEN(*)
      INTEGER, INTENT(OUT) :: IPART
      REAL(DP) :: PHIANG(10)
      REAL(DP) :: A0, A1, A2, A3, A4, XL1, XL2, EPS10, PI, PI180, YL1,
     .          YL2, PHI, PHIH, XPHI, YPHI, XM, YM
      INTEGER :: I, IPHI, ISORT
      LOGICAL :: LX0, LY0, LZ0, LCTX1, LCTX2, LCTY1, LCTY2
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
        WRITE (iunout,*) ' ERROR IN CTQUA '
        WRITE (iunout,*) ' EQUATION DOES NOT DESCRIBE A CYLINDER '
        WRITE (iunout,*) ' IT DESCRIBES TWO PARALLEL PLANES '
        WRITE (iunout,*) ' NO PLOT IS DONE '
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
C ===== SOURCE: ellipsoid.f


      SUBROUTINE ELLIPSOID (X0,Y0,Z0,CX,CY,CZ,XLIMS1,YLIMS1,ZLIMS1,
     .                      XLIMS2,YLIMS2,ZLIMS2,RLB,ILCOL,NX,NY,NZ)

      USE PRECISION
      USE CCONA

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X0, Y0, Z0, CX, CY, CZ, XLIMS1, YLIMS1,
     .                      ZLIMS1, XLIMS2, YLIMS2, ZLIMS2, RLB
      INTEGER, INTENT(IN) :: ILCOL, NX, NY, NZ
      REAL(DP) :: XP(2*101), YP(2*101), PHIAN(9), PHIEN(9)
      REAL(DP) :: RX, RY, RZ, DPH, DZ, X, Y, Z, ZE, RAD, XS, YS, PHI,
     .          DX, DY, XA, XE, YA, YE, ZA
      INTEGER :: IPART, IPRT, N, IZ, IX, IY, I
      INTEGER :: IP(2*101)

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
          CALL GRJMP (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            X = X0 + RX*COS(PHI)
            Y = Y0 + RY*SIN(PHI)
            CALL PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
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
          CALL GRJMP (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            X = X0 + RX*COS(PHI)
            Z = Z0 + RZ*SIN(PHI)
            CALL PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
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
          CALL GRJMP (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            Y = Y0 + RY*COS(PHI)
            Z = Z0 + RZ*SIN(PHI)
            CALL PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          END DO
        END DO
      END DO ! IY

      CALL GRNWPN(1)

      RETURN
      END
C ===== SOURCE: ello.f
C
C
      FUNCTION ELLO (X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: ELLO
      ELLO=BHALB*SQRT(MAX(0.D0,1.D0-X*X/(AHALB*AHALB)))
      RETURN
      END
C ===== SOURCE: ellu.f
C
C
      FUNCTION ELLU (X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: ELLU
      ELLU=-BHALB*SQRT(MAX(0.D0,1.D0-X*X/(AHALB*AHALB)))
      RETURN
      END
C ===== SOURCE: exit.f
C
C
      SUBROUTINE EXIT(ICC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ICC
      CALL GREND
      STOP
      END
C ===== SOURCE: exit_own.f
C
C
      SUBROUTINE EXIT_OWN (ICC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ICC
      CALL GREND
      STOP
      END SUBROUTINE EXIT_OWN
C ===== SOURCE: flaech.f


      FUNCTION FLAECH (X1,Y1,X2,Y2,X3,Y3,X4,Y4)

      USE PRECISION

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, Y1, X2, Y2, X3, Y3, X4, Y4
      REAL(DP) :: FLAECH1, FLAECH2, FLAECH

      FLAECH1 = 0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
      FLAECH2 = 0.5*(X1*(Y3-Y4)+X3*(Y4-Y1)+X4*(Y1-Y3))
      FLAECH = ABS(FLAECH1) + ABS(FLAECH2)
      END
C ===== SOURCE: fmu.f
C
C
      FUNCTION FMU (V1,V2,V3,VR1,VR2,VR3,P1,P2,P3,PR1,PR2,PR3,EPS)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: V1, V2, V3, VR1, VR2, VR3,
     .                      P1, P2, P3, PR1, PR2, PR3, EPS
      REAL(DP) :: XNEN, FMU

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
C     WRITE (iunout,*) ' SCHNITTGERADE PARALLEL ZU VIERECKSEITE '
      FMU=10.
      RETURN
      END
C ===== SOURCE: geradx.f
C
C
      FUNCTION GERADX(Y)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: Y
      REAL(DP) :: GERADX
      GERADX=Y
      Y=A0*Y+A1
      RETURN
      END
C ===== SOURCE: gerady.f
C
C
      FUNCTION GERADY(X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: GERADY
      GERADY=A0*X+A1
      RETURN
      END
C ===== SOURCE: grnxtb.f
C  june 05:  new option: K=3 (wg. iteravitve mode)
C  june 05:  name: calling routine
C
      SUBROUTINE GRNXTB(K,name)
C
C  K=1 CALL VOR EINEM BILD AUS EIGENER GRSOFTWARE
C  K=2 CALL VOR EINEM BILD MIT GRBLD  (KURVEF,....)
C  K=3 CALL NACH EINER VOLLEN ITERATION. RESET IFRST1,IFRST2
C                                        CLOSE PICTURES
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: K
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      INTEGER, SAVE :: IFRST1=0, IFRST2=0

      WRITE (iunout,*) 'GRNXTB CALLED FROM ',NAME
      WRITE (iunout,*) 'K,IFRST1,IFRST2 ',K,IFRST1,IFRST2

      GOTO (1,2,3),K
1     IF (IFRST1.EQ.1) THEN
        CALL GRNXTF
        RETURN
      ELSE
        IFRST1=1
      ENDIF
      RETURN
C
2     IF (IFRST1.EQ.0) THEN
        CALL GRSCLC (0.,0.0,39.5,28.7)
        CALL GRSCLC (3.,3.5,39.5,28.7)
        IFRST1=1
        IFRST2=1
      ELSE IF (IFRST2.EQ.0) THEN
        CALL GRNXTF
        CALL GRSCLC (0.,0.,39.5,28.7)
        CALL GRSCLC (3.,3.5,39.5,28.7)
        IFRST2=1
      ELSE
        CALL GRSCLC (3.,3.5,39.5,28.7)
      ENDIF
      RETURN

3     CONTINUE
      IF (IFRST1.EQ.1) CALL GRNXTF
      IFRST1=0
      IF (IFRST2.EQ.1) THEN
        CALL GRSCLC (0.,0.,39.5,28.7)
        CALL GRSCLC (3.,3.5,39.5,28.7)
      ENDIF
      IFRST2=0
      RETURN
      END
C ===== SOURCE: gsp.f
C
C
      SUBROUTINE GSP (P1,P2,R1,R2,Q1,Q2,S1,S2,XLA,XMU,EPS)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: P1, P2, R1, R2, Q1, Q2, S1, S2, EPS
      REAL(DP), INTENT(OUT) :: XLA, XMU
      REAL(DP) :: S3, R3, PK1, PK2, PK3, XNO

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
C ===== SOURCE: hypm.f
C
C
      FUNCTION HYPM(X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: HYPM, DH, RAD, AH
      AH=A*X*X
      DH=D*X
      RAD=E*E/(4.*C*C)-(F+AH+DH)/C
      IF (RAD.GE.0.) GOTO 1
      HYPM=1.D50
      RETURN
1     HYPM=-E/(2.*C)-SQRT(RAD)
      RETURN
      END
C ===== SOURCE: hypp.f
C
C
      FUNCTION HYPP(X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: HYPP, DH, RAD, AH
      AH=A*X*X
      DH=D*X
      RAD=E*E/(4.*C*C)-(F+AH+DH)/C
      IF (RAD.GE.0.) GOTO 1
      HYPP=1.D50
      RETURN
1     HYPP=-E/(2.*C)+SQRT(RAD)
      RETURN
      END
C ===== SOURCE: isolne.f


      SUBROUTINE ISOLNE (AORIG,IBLD,ICURV,
     .                   IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE CARRIES OUT A CONTOUR PLOT
C
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ONTO VERTICES IS PERFORMED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTRIG

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV, IXX, IYY
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: ZINT, SCLFCX, SCLFCY, RAMIN, RAMAX, FAK, CM, DX, DY,
     .          A1, A2, A3, A4, ACMIN, ACMAX, DA, ACONT, RMI, XMIN,
     .          XMAX, YMIN, YMAX, RMA, X1, X2, AA1, AA2
      REAL(DP), ALLOCATABLE :: A(:,:),AA(:,:)
      REAL(SP) :: XY(8000)
      REAL(SP) :: YH
      INTEGER :: NP, ITR, NP1, NP2, ICOLOR, IS, IC, IISO, NISO, IERR,
     .           IPART, IT, IR, IP, I
      CHARACTER(17) :: CH
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
      IF (LEVGEO .LE. 2.AND.LPTOR3(IBLD)) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 21 IR=1,IXX-1
          DO 21 IP=1,IYY-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
21      CONTINUE
      ELSEIF (LEVGEO .LE. 2.AND.LPPOL3(IBLD)) THEN
        IP=1
        IF (NLPOL) IP=IPROJ3(IBLD,ICURV)
        IF (IP.LE.0.OR.IP.GT.NP2ND) IP=1
        DO 23 IR=1,IXX-1
          DO 23 IT=1,IYY-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
23      CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
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
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
        DO 22 I=1,NTRII
          RMI=MIN(RMI,AORIG(I))
          RMA=MAX(RMA,AORIG(I))
22      CONTINUE
      ELSE
        WRITE (iunout,*) 'MISSING OPTION IN ISOLNE: RMI,RMA '
        WRITE (iunout,*) 'PLOT ABANDONNED '
        RETURN
      ENDIF
C
C  INTERPOLATE: AORIG --> A
C
      IF (LEVGEO.EQ.4) THEN
        ALLOCATE (AA(NRAD,1))
        CALL CELINT(AORIG,AA,LOGL,IBLD,ICURV,NRAD,IERR)
      ELSEIF (LEVGEO.LE.3) THEN
         ALLOCATE (A(N1ST,MAX(N2ND,N3RD)))
        CALL CELINT(AORIG,A,LOGL,IBLD,ICURV,N1ST,IERR)
      ENDIF
      IF (IERR.GT.0) RETURN
C
C  SEARCH FOR XMIN,XMAX,YMIN,YMAX
C
      XMIN=1.D60
      YMIN=1.D60
      XMAX=-1.D60
      YMAX=-1.D60
C
      IF ((LEVGEO.LE.2).AND.LPPOL3(IBLD)) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = ZSURF(1)
        YMAX = ZSURF(NT3RD)
      ELSEIF (LEVGEO.EQ.1.AND.LPTOR3(IBLD)) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = PSURF(1)
        YMAX = PSURF(NP2ND)
      ELSEIF (LEVGEO.EQ.2.AND.LPTOR3(IBLD)) THEN
C  SUFFICIENT TO SEARCH ON OUTERMOST RADIAL SURFACE (BECAUSE: CONVEX)
        DO 5 IP = 1,NP2ND
          XMIN = MIN(XMIN,XPOL(NR1ST,IP))
          XMAX = MAX(XMAX,XPOL(NR1ST,IP))
          YMIN = MIN(YMIN,YPOL(NR1ST,IP))
          YMAX = MAX(YMAX,YPOL(NR1ST,IP))
5       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
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
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
C  SEARCH ON WHOLE MESH
        DO I=1,NRKNOT
          XMIN=MIN(XMIN,XTRIAN(I))
          YMIN=MIN(YMIN,YTRIAN(I))
          XMAX=MAX(XMAX,XTRIAN(I))
          YMAX=MAX(YMAX,YTRIAN(I))
        ENDDO
      ELSE
        WRITE (iunout,*) 
     .    'MISSING OPTION IN ISOLNE: XMIN,XMAX,YMIN,YMAX '
        WRITE (iunout,*) 'PLOT ABANDONNED '
        RETURN
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
      CALL GRNXTB (1,'ISOLNE.F')
      CALL GRSCLC (10.,4.,REAL(10.+DX*FAK,KIND(1.E0)),
     .                    REAL(4.+DY*FAK,KIND(1.E0)))
      CALL GRSCLV (REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)),
     .             REAL(XMAX,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
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
      IF (LEVGEO.EQ.1) THEN
        CALL GRJMP(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
      ELSEIF (LEVGEO.EQ.2.AND.LPPOL3(IBLD)) THEN
        CALL GRJMP(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
      ELSEIF (LEVGEO.EQ.2.AND.LPTOR3(IBLD)) THEN
        DO 7 IR=1,NR1ST,NR1STM
          CALL GRJMP(real(XPOL(IR,1),KIND(1.E0)),
     .               real(YPOL(IR,1),KIND(1.E0)))
          DO 9 IP = 2,NP2ND
9           CALL GRDRW(real(XPOL(IR,IP),KIND(1.E0)),
     .                 real(YPOL(IR,IP),KIND(1.E0)))
7       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
        CALL GRJMP(real(XPOL(1,NPOINT(1,1)),KIND(1.E0)),
     .             real(YPOL(1,NPOINT(1,1)),KIND(1.E0)))
          DO 10 IR=2,NR1ST
10          CALL GRDRW (real(XPOL(IR,NPOINT(1,1)),KIND(1.E0)),
     .                  real(YPOL(IR,NPOINT(1,1)),KIND(1.E0)))
C
        CALL GRJMP (real(XPOL(1,NPOINT(2,NPPLG)),KIND(1.E0)),
     .              real(YPOL(1,NPOINT(2,NPPLG)),KIND(1.E0)))
        DO 11 IR=2,NR1ST
          NP=NPOINT(2,NPPLG)
11        CALL GRDRW (real(XPOL(IR,NP),KIND(1.E0)),
     .                real(YPOL(IR,NP),KIND(1.E0)))
        DO 15 I=1,NPPLG
          CALL GRJMP (real(XPOL(1,NPOINT(1,I)),KIND(1.E0)),
     .                real(YPOL(1,NPOINT(1,I)),KIND(1.E0)))
          DO 12 IP=NPOINT(1,I),NPOINT(2,I)
12          CALL GRDRW (real(XPOL(1,IP),KIND(1.E0)),
     .                  real(YPOL(1,IP),KIND(1.E0)))
          CALL GRJMP (real(XPOL(NR1ST,NPOINT(1,I)),KIND(1.E0)),
     .                real(YPOL(NR1ST,NPOINT(1,I)),KIND(1.E0)))
          DO 13 IP=NPOINT(1,I),NPOINT(2,I)
13          CALL GRDRW (real(XPOL(NR1ST,IP),KIND(1.E0)),
     .                  real(YPOL(NR1ST,IP),KIND(1.E0)))
15      CONTINUE
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
        DO ITR=1,NTRII
          IF (NCHBAR(1,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(1,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(2,ITR)),KIND(1.E0)))
          ENDIF
          IF (NCHBAR(2,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(3,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(3,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(2,ITR)),KIND(1.E0)))
          ENDIF
          IF (NCHBAR(3,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(1,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(3,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(3,ITR)),KIND(1.E0)))
          ENDIF
        ENDDO
      ELSE
        WRITE (iunout,*) 'MISSING OPTION IN ISOLNE: PLOT GRID BOUNDARY '
        WRITE (iunout,*) 'PLOT ABANDONNED '
        IF (ALLOCATED(A)) DEALLOCATE (A)
        IF (ALLOCATED(AA)) DEALLOCATE (AA)
        RETURN
      ENDIF
C
      IF (LOGL) THEN
        IF (ZMI.EQ.666.) THEN
          RAMIN=LOG10(MAX(1.E-48_DP,RMI))
        ELSE
          RAMIN=LOG10(MAX(1.E-48_DP,ZMI))
        ENDIF
        IF (ZMA.EQ.666.) THEN
          RAMAX=LOG10(MAX(1.E-48_DP,RMA))
        ELSE
          RAMAX=LOG10(MAX(1.E-48_DP,ZMA))
        ENDIF
      ELSE
        RAMIN=ZMI
        IF (ZMI.EQ.666.) RAMIN=RMI
        RAMAX=ZMA
        IF (ZMA.EQ.666.) RAMAX=RMA
      ENDIF
C
      IF (ABS(RAMAX-RAMIN)/MAX(RAMAX,1.E-30_DP).LT.0.01)
     .    RAMAX=RAMIN+0.01*RAMIN*SIGN(1._DP,RAMIN)
C
C
      DA=(RAMAX-RAMIN)/DBLE(NISO-1)
C
      ICOLOR=1
      IF ((LEVGEO.EQ.1).OR.(LEVGEO.EQ.2.AND.LPPOL3(IBLD))) THEN
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
C               WRITE (iunout,*) ' ACHTUNG !!!! '
C               WRITE (iunout,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (iunout,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
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
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTOR3(IBLD)) THEN
        DO 100 IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO 110 IR=1,NR1ST-1
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
C               WRITE (iunout,*) ' ACHTUNG !!!! '
C               WRITE (iunout,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (iunout,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
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
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
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
C               WRITE (iunout,*) ' ACHTUNG !!!! '
C               WRITE (iunout,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (iunout,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
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
      ELSE
        WRITE (iunout,*) 'MISSING OPTION IN ISOLNE: PLOT CONTOURS '
        WRITE (iunout,*) 'PLOT ABANDONNED '
        IF (ALLOCATED(A)) DEALLOCATE (A)
        IF (ALLOCATED(AA)) DEALLOCATE (AA)
        RETURN
      ENDIF
C
C
C     WRITE TEXT AND MEAN VALUE ONTO THE PLOT
C
      CALL GRNWPN (1)
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,RUNID)
      YH=26.75
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,HEAD)
      YH=26.00
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,REAL(YH-0.5,KIND(1.E0)),10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,REAL(YH-1.,KIND(1.E0)),10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,REAL(YH-2.,KIND(1.E0)),10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') RMA
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') RMI
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)
C
      ICOLOR=1
      YH=YH-4.
      DO 200 IS=1,NISO
        ACONT=RAMIN+(IS-1)*DA
        IF (LOGL) ACONT=10.**ACONT
        IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
        CALL GRNWPN(ICOLOR)
        YH=YH-0.5
        CALL GRJMP (1.,REAL(YH+0.25,KIND(1.E0)))
        CALL GRDRW (2.5,REAL(YH+0.25,KIND(1.E0)))
        CALL GRNWPN (1)
        WRITE (CH,'(1P,E10.3)') ACONT
        CALL GRTXT (3.,REAL(YH,KIND(1.E0)),10,CH)
200   CONTINUE
C
      IF (ALLOCATED(A)) DEALLOCATE (A)
      IF (ALLOCATED(AA)) DEALLOCATE (AA)

      RETURN
      END
C ===== SOURCE: loctor.f
C
C
      SUBROUTINE LOCTOR(W3,RM,X1,Z1)
C
C TRANSFORM FROM LOCAL SYSTEM AT PHI=W3 TO TORUS SYSTEM
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: W3, RM
      REAL(DP), INTENT(INOUT) :: X1, Z1
      REAL(DP), SAVE :: WO=0.D0, S=0.D0, C=1.D0
      REAL(DP) :: XS

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
C ===== SOURCE: muelam.f
C
C
      SUBROUTINE MUELAM (X1,Y1,VX,VY,X2,Y2,X3,Y3,XL,XM)

      USE PRECISION

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, Y1, VX, VY, X2, Y2, X3, Y3
      REAL(DP), INTENT(OUT) :: XL, XM
      REAL(DP) :: D

      D=VX*(Y2-Y3)-VY*(X2-X3)+1.E-20
      XL=((Y2-Y3)*(X2-X1)-(X2-X3)*(Y2-Y1))/D
      XM=(-VY*(X2-X1)+VX*(Y2-Y1))/D

      RETURN
      END
C ===== SOURCE: para1.f
C
C
      FUNCTION PARA1 (X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: PARA1
      PARA1=-A/E*X*X-D/E*X-F/E
      RETURN
      END
C ===== SOURCE: para2o.f
C
C
      FUNCTION PARA2O (XIN)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XIN
      REAL(DP) :: RAD, Y, DH, PARA2O
      DH=-E/(2.*C)
      RAD=(E*E)/(4.*C*C)-(D*XIN+F)/C
      IF (RAD.GE.0.) GOTO 1
      PARA2O=1.D50
      RETURN
1     Y=SQRT(RAD)
      PARA2O=Y-DH
      RETURN
      END
C ===== SOURCE: para2u.f
C
C
      FUNCTION PARA2U (XIN)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XIN
      REAL(DP) :: RAD, Y, DH, PARA2U
      DH=-E/(2.*C)
      RAD=(E*E)/(4.*C*C)-(D*XIN+F)/C
      IF (RAD.GE.0.) GOTO 1
      PARA2U=1.D50
      RETURN
1     Y=-SQRT(RAD)
      PARA2U=Y-DH
      RETURN
      END
C ===== SOURCE: pl3d.f
C
C
      SUBROUTINE PL3D(PX,PY,PZ,PP1,PP2)

      USE PRECISION
      USE PARMMOD
      USE CPL3D
      USE CPLOT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PX, PY, PZ
      REAL(DP), INTENT(OUT) :: PP1, PP2
      REAL(DP) :: FXI, FZETA, X, Y, Z, COP, COT, SIP, SIT, PPY11, PPX11,
     .          YO, ZO, DY, DZ, CH, TH, PPX12, PPX13, PPX14, PPX21,
     .          PPX22, PPX23, PPX24, PPY12, PPY13, PPY14, PPY21, PPY22,
     .          PPY23, PPY24, DX, F00, F01, F02, F03, F10, F11, F12,
     .          BREITE, XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, XO, W3, FA,
     .          BREI
      INTEGER :: IFIRST, IPBOX
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
        CALL GRJMPS(REAL(PP1),REAL(PP2),203)
        CALL GRMRKS(0.2)
        PPX11= FXI(XMIN,ZMIN)
        PPY11= FZETA(XMIN,ZMIN,YMIN)
        CALL GRJMP(REAL(PPX11,KIND(1.E0)),REAL(PPY11,KIND(1.E0)))
        PPX12 = FXI(XMAX,ZMIN)
        PPY12 = FZETA(XMAX,ZMIN,YMIN)
        CALL GRDRW(REAL(PPX12,KIND(1.E0)),REAL(PPY12,KIND(1.E0)))
        PPX13 = FXI(XMAX,ZMIN)
        PPY13 = FZETA(XMAX,ZMIN,YMAX)
        CALL GRDRW(REAL(PPX13,KIND(1.E0)),REAL(PPY13,KIND(1.E0)))
        PPX14 = FXI(XMIN,ZMIN)
        PPY14 = FZETA(XMIN,ZMIN,YMAX)
        CALL GRDRW(REAL(PPX14,KIND(1.E0)),REAL(PPY14,KIND(1.E0)))
        CALL GRDRW(REAL(PPX11,KIND(1.E0)),REAL(PPY11,KIND(1.E0)))
        PPX21= FXI(XMIN,ZMAX)
        PPY21= FZETA(XMIN,ZMAX,YMIN)
        CALL GRDRW(REAL(PPX21,KIND(1.E0)),REAL(PPY21,KIND(1.E0)))
        PPX22 = FXI(XMAX,ZMAX)
        PPY22 = FZETA(XMAX,ZMAX,YMIN)
        CALL GRDRW(REAL(PPX22,KIND(1.E0)),REAL(PPY22,KIND(1.E0)))
        PPX23 = FXI(XMAX,ZMAX)
        PPY23 = FZETA(XMAX,ZMAX,YMAX)
        CALL GRDRW(REAL(PPX23,KIND(1.E0)),REAL(PPY23,KIND(1.E0)))
        PPX24 = FXI(XMIN,ZMAX)
        PPY24 = FZETA(XMIN,ZMAX,YMAX)
        CALL GRDRW(REAL(PPX24,KIND(1.E0)),REAL(PPY24,KIND(1.E0)))
        CALL GRDRW(REAL(PPX21,KIND(1.E0)),REAL(PPY21,KIND(1.E0)))
        CALL GRJMP(REAL(PPX12,KIND(1.E0)),REAL(PPY12,KIND(1.E0)))
        CALL GRDRW(REAL(PPX22,KIND(1.E0)),REAL(PPY22,KIND(1.E0)))
        CALL GRJMP(REAL(PPX13,KIND(1.E0)),REAL(PPY13,KIND(1.E0)))
        CALL GRDRW(REAL(PPX23,KIND(1.E0)),REAL(PPY23,KIND(1.E0)))
        CALL GRJMP(REAL(PPX14,KIND(1.E0)),REAL(PPY14,KIND(1.E0)))
        CALL GRDRW(REAL(PPX24,KIND(1.E0)),REAL(PPY24,KIND(1.E0)))
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

C     the following ENTRY is for reintialization of EIRE (DMH)
      
      ENTRY PL3D_REINIT
      IFIRST = 0
      IPBOX = 0
      return

      END
C ===== SOURCE: pl3dpg.f
C  3D HISTOGRAM PLOT
C
      SUBROUTINE PL3DPG (ARR,IBLD,ICURV,
     .                   IX,IY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,LOGL,
     .                   ZMA,ZMI,W1,W2,
     .                   HEAD,RUNID,TXHEAD,TRC)

      USE PRECISION
      USE PARMMOD
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      INTEGER, PARAMETER :: LAR=46*128*128

      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(INOUT) :: ARR(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI, W1, W2
      INTEGER, INTENT(IN) :: IBLD, ICURV, IX, IY
      LOGICAL, INTENT(IN) :: LOGL,TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: REMIN, REMAX, XMT, DX, YMI, YMA, RMI, RMA, AAR, XMINN,
     .          XMAXN, YMINN, YMAXN, YMT, XMI, XMA
!pb      REAL(SP) :: AR(LAR), EXT(3,3), VALU(3,2)
!pb      REAL(SP) :: XYZ(3,128,128)
!pb      real(sp) :: yh
      REAL :: AR(LAR), EXT(3,3), VALU(3,2)
      REAL :: XYZ(3,128,128)
      real :: yh
      INTEGER :: IR, IPX, IPY, IER, IPAN, IPEN, K, I, J
      CHARACTER(17) :: CH
      CHARACTER(20) :: CHAXS(3)
C
      DO 1 I=1,3
      DO 1 J=1,128
      DO 1 K=1,128
1       XYZ(I,J,K)=-75.75E20
C
      IF (LEVGEO.LE.1.OR.LEVGEO.GT.3.OR..NOT.LPTOR3(IBLD)) THEN
        WRITE (iunout,*) 'PLOTOPTION NOT READY. RETURN FROM PL3DPG '
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
                RMI=MIN(RMI,MAX(1.E-48_DP,ARR(IR)))
                RMA=MAX(RMA,MAX(1.E-48_DP,ARR(IR)))
                ARR(IR)=LOG10(MAX(1.E-48_DP,ARR(IR)))
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
          REMIN=LOG10(MAX(1.E-48_DP,ZMI))
        ENDIF
        IF (ZMA.EQ.666.) THEN
          REMAX=LOG10(RMA)
        ELSE
          REMAX=LOG10(MAX(1.E-48_DP,ZMA))
        ENDIF
      ELSE
        IF (ZMI.EQ.666.) REMIN=RMI
        IF (ZMA.EQ.666.) REMAX=RMA
      ENDIF
C
      IF (TRC) THEN
        WRITE (iunout,*) 
     .    'PL3DPG:  ,IPAN,IPEN,XMI,XMA,YMI,YMA,REMIN,REMAX'
        WRITE (iunout,*) 
     .    '        ',IPAN,IPEN,XMI,XMA,YMI,YMA,REMIN,REMAX
      ENDIF
C
      IF (ABS(REMAX-REMIN)/MAX(REMAX,1.E-30_DP).LT.0.01)
     .    REMAX=REMIN+0.01*REMIN*SIGN(1._DP,REMIN)
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
      CALL GRNXTB (1,'PL3DPG.F')
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
      CALL GR3ROT(AR,IER,'Z',REAL(W1,KIND(1.E0)),
     .            'X',REAL(W2,KIND(1.E0)),'Y',0.0)
      CALL GR3PLO(AR,IER,'HID')
C
C     WRITE TEXT AND MEAN VALUE ONTO THE PLOT
C
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,RUNID)
      YH=26.75
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,HEAD)
      YH=26.00
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,REAL(YH-0.5,KIND(1.E0)),10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,REAL(YH-1.,KIND(1.E0)),10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,REAL(YH-2.,KIND(1.E0)),10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') RMA
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') RMI
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'NOT ENOUGH STORAGE FOR 3D HISTOGRAM PLOT'
      WRITE (iunout,*) 'REDUCE PLOT AREA '
      WRITE (iunout,*) 'PLOT ABANDONNED'
      RETURN
      END
C ===== SOURCE: pl3q.f
C
C
      SUBROUTINE PL3Q(CORD,N,IO,NF)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: CORD(*)
      INTEGER, INTENT(IN) :: N, IO
      LOGICAL, INTENT(IN) :: NF
      REAL(DP) :: XP(N+1),YP(N+1)
      REAL(SP) :: XPS(N+1),YPS(N+1)
      INTEGER :: J, JJ, I, N3

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
C ===== SOURCE: plane.f
C
C
      SUBROUTINE PLANE (A0,A1,A2,A3,RL,N1,EPS,
     .  AL,XL,YL,ZL,AL0,XL1,YL1,ZL1,XL2,YL2,ZL2,XL3,YL3,ZL3,IO,NF,NUM)
C
C  PLOT PLANE, SECTION INSIDE A BOX
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A0, A1, A2, A3, RL, EPS
      INTEGER, INTENT(IN) :: N1, IO, NUM
      REAL(DP), INTENT(IN) :: AL(N1,*),  XL(N1,*),  YL(N1,*),  ZL(N1,*),
     .                      AL0(N1,*), XL1(N1,*), YL1(N1,*), ZL1(N1,*),
     .                      XL2(N1,*), YL2(N1,*), ZL2(N1,*),
     .                      XL3(N1,*), YL3(N1,*), ZL3(N1,*)
      LOGICAL NF
      REAL(DP) :: P(3,36), XYZG(3,36), ANGLE(36), CORD(108), PS(3)
      REAL(DP) :: B1, B2, B3, DET, DETER, TEST, T, DX, DY, DZ, HELP,
     .          XMIT, YMIT, ZMIT, PI, ANG, X1, X2, Y1, Y2, Z1, Z2
      INTEGER :: K, J, IPOINT, IP, ILN, I, II, ISORT, ICOUNT, ICHECK,
     .           IS
C
      IS=0
      IF (RL.GT.0) THEN
        X1=XL1(1,NUM)
        X2=XL2(1,NUM)
        Y1=YL1(1,NUM)
        Y2=YL2(1,NUM)
        Z1=ZL1(1,NUM)
        Z2=ZL2(1,NUM)
        DX=XL2(1,NUM)-XL1(1,NUM)
        DY=YL2(1,NUM)-YL1(1,NUM)
        DZ=ZL2(1,NUM)-ZL1(1,NUM)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z1,DX,0._DP,0._DP,P(1,1),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z1,0._DP,0._DP,DZ,P(1,2),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z2,DX,0._DP,0._DP,P(1,3),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y1,Z2,0._DP,0._DP,DZ,P(1,4),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y1,Z1,0._DP,DY,0._DP,P(1,5),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y2,Z1,DX,0._DP,0._DP,P(1,6),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y2,Z1,0._DP,0._DP,DZ,P(1,7),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y2,Z2,DX,0._DP,0._DP,P(1,8),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y2,Z2,0._DP,0._DP,DZ,P(1,9),T)
        CALL SPOINT (A0,A1,A2,A3,X2,Y1,Z2,0._DP,DY,0._DP,P(1,10),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z1,0._DP,DY,0._DP,P(1,11),T)
        CALL SPOINT (A0,A1,A2,A3,X1,Y1,Z2,0._DP,DY,0._DP,P(1,12),T)
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
            DET=DETER(A1,XL(I,NUM),XL(J,NUM),A2,YL(I,NUM),YL(J,NUM),
     .                A3,ZL(I,NUM),ZL(J,NUM))
            IF (ABS(DET).LE.EPS) GOTO 21
            B1=-A0
            B2=-AL(I,NUM)
            B3=-AL(J,NUM)
            PS(1)=DETER(B1,B2,B3,A2,YL(I,NUM),YL(J,NUM),
     .                  A3,ZL(I,NUM),ZL(J,NUM))/DET
            PS(2)=DETER(A1,XL(I,NUM),XL(J,NUM),B1,B2,B3,
     .                  A3,ZL(I,NUM),ZL(J,NUM))/DET
            PS(3)=DETER(A1,XL(I,NUM),XL(J,NUM),A2,YL(I,NUM),YL(J,NUM),
     .                  B1,B2,B3)/DET
C  CHECKE, OB ALLE ANDEREN LINEAREN UNGLEICHUNGEN ERFUELLT SIND
            DO 40 K=1,ILN
              IF (K.EQ.I.OR.K.EQ.J) GOTO 40
              TEST=PS(1)*XL(K,NUM)+PS(2)*YL(K,NUM)+
     .             PS(3)*ZL(K,NUM)+AL(K,NUM)
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
        WRITE (iunout,*) 'WENIGER ALS 3 SCHNITTPUNKTE ',
     .              'DER EBENE MIT DEM QUADER GEFUNDEN'
        WRITE (iunout,*) 'NUMMER DER FLAECHE: ',NUM
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
C ===== SOURCE: plgell.f
C
C
      SUBROUTINE PLGELL
C
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
      USE PRECISION
      USE CCONA

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: THETA(*)
      REAL(DP), INTENT(OUT) :: XX(*), YY(*)
      REAL(DP), INTENT(IN) :: AHALB, EP, ELL, TRI, DM
      INTEGER, INTENT(IN) :: N, NP
      INTEGER, INTENT(OUT) :: NRET
      REAL(DP) :: XOLD, YOLD, ATR, BTR, DX, XNEW, YNEW, DMM, XAD, XH,
     .          BHB, AHB
      INTEGER :: ICOU, I
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
C ===== SOURCE: plot3d.f
C  3D SMOOTH SURFACE PLOT
C
      SUBROUTINE PLOT3D (ARR,IBLD,ICURV,
     .                   NX,NY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,LOGL,
     .                   ZMA,ZMI,W1,W2,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  ON INPUT:  LOGL: USE LOG SCALE FOR ORDINATE
C             ARR(IJ),I=1,NX-1,J=1,NY-1 ,IJ=I+(J-1)*NR1ST
C                          ARRAY TO BE PLOTTED
C             XX(I),I=1,NX X-GRID BOUNDARIES
C             YY(J),J=1,NY Y GRID BOUNDARIES
C  ARR(IJ) IS THEN SET ONTO 2D ARRAY FALT(I,J)
C
C
C  LPOLAR: R-THETA CO-ORDINATES
C  LKARTH: X-Y     CO-ORDINATES
C
C  FALT-->XZY (3,...)
C
      USE PRECISION
      USE PARMMOD
      USE CPLOT
      USE CGRID
      USE CGEOM
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(INOUT) :: ARR(*)
      REAL(DP), INTENT(IN) :: ZMI, ZMA, W1, W2
      INTEGER, INTENT(IN) :: NX, NY, IBLD, ICURV
      LOGICAL, INTENT(IN) :: LOGL,TRC
      CHARACTER(72), INTENT(IN) :: TEXT1
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3
      CHARACTER(72), INTENT(IN) :: HEAD, RUNID, TXHEAD

      REAL(DP) :: XMINN, XMAXN, DXXX, YMINN, YMAXN
!pb      REAL(SP) :: FALT(N1ST,N2NDPLG),
!pb     .          X(N1ST,N2NDPLG), Y(N1ST,N2NDPLG), Z2(N1ST,N2NDPLG)
!pb      REAL(SP) :: EXT(3,3), VALU(3,2), DCM, YH
!pb      REAL(SP) :: YHLF, XHLF, FMIN, FMAX, REMIN, REMAX,
!pb     .          XMIN, XMAX, YMIN, YMAX
!pb      REAL(SP), ALLOCATABLE :: AR(:)
      REAL :: FALT(N1ST,N2NDPLG),
     .        X(N1ST,N2NDPLG), Y(N1ST,N2NDPLG), Z2(N1ST,N2NDPLG)
      REAL :: EXT(3,3), VALU(3,2), DCM, YH
      REAL :: YHLF, XHLF, FMIN, FMAX, REMIN, REMAX,
     .        XMIN, XMAX, YMIN, YMAX
      REAL, ALLOCATABLE :: AR(:)
      INTEGER :: I, J, IXM, IYM, IZ, IZN, IER, LAR, IX, IY
      CHARACTER(17) :: CH
      CHARACTER(20) :: CHAXS(3)
C
C
      LAR=46*N1ST*N2NDPLG
      ALLOCATE (AR(LAR))
      ier=1
      IF (LEVGEO.LE.3.AND.LEVGEO.GT.1.AND.LPTOR3(IBLD)) THEN
C
        IXM=NX-1
        IYM=NY-1
        IZ=IYM*NR1ST
C
        DO 1 I=1,N1ST
          DO 1 J=1,N2NDPLG
            FALT(I,J)=-75.75E20
1       CONTINUE
C
        IF (LOGL) THEN
          DO 3 J=1,IZ
            ARR(J)=LOG10(MAX(1.E-48_DP,ARR(J)))
3         CONTINUE
        ENDIF
C
C  SET ONTO 2D ARRAY FOR PLOTTING
C
        DO 20 I=1,NX
          DO 20 J=1,IYM
            DXXX=ARR(I+(J-1)*NR1ST)
            FALT(I,J)=DXXX
            X(I,J)=XPOL(I,J)
            Y(I,J)=YPOL(I,J)
20      CONTINUE
C
        IXM=NX-1
        IYM=NY-1
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
          REMIN=LOG10(MAX(1.E-48_DP,ZMI))
          REMAX=LOG10(MAX(1.E-48_DP,ZMA))
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
        CALL GRNXTB (1,'PLOT3D.F')
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
        CALL GR3ROT(AR,IER,'Z',REAL(W1,KIND(1.E0)),'X',
     .              REAL(W2,KIND(1.E0)),'Y',0.0)
        CALL GR3PLO(AR,IER,'HID')
C
      ELSE
        WRITE (iunout,*) 'INVALID OPTION IN PLOT3D  '
        WRITE (iunout,*) '3D-PLOT ABANDONED  '
        RETURN
      ENDIF
C
C     WRITE TEXT ONTO THE PLOT
C
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,RUNID)
      YH=26.75
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,HEAD)
      YH=26.00
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,REAL(YH-0.5,KIND(1.E0)),10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,REAL(YH-1.,KIND(1.E0)),10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,REAL(YH-2.,KIND(1.E0)),10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') FMAX
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') FMIN
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)

      DEALLOCATE (AR)
C
      RETURN
      END
C ===== SOURCE: plt2d.f
cdr  28.4.04:  nhsts(ispz) option connected (to select species
cdr            for trajectory plot. see modification to input.f, 28.4.04
cdr  24.8.06:  plot symbols corrected to more recent GR  software standards
!pb  5.10.06:  plot for triangle geometry in x-z plane added

C   2D GEOMETRY (AND TRAJECTORY) PLOT

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
      USE CFPLK
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
     .          XT2, XTIP, P, YTIP, TR, RS, EP, EL, XT, YT, XTN, YTN,
     .          DXX, DYY, A, B
      REAL(SP) :: XPS(5), YPS(5)
      INTEGER :: ISPL(NTXHST),ICLR(2*NSTS+1),IDSH(2*NSTS+1),
     .           ISWC(2*NSTS+1)
      INTEGER, ALLOCATABLE :: IFARB(:,:), IDASH(:,:), ICPSPZ(:)
      INTEGER :: ICP, ISTR, IC, IC1, IC2, NCTPNT, IDUMMY, ICT, IEN,
     .           ITH, ITHPL, IAN, NTDUM, NTT, IECKE2,
     .           LEARCA, NT, ICOLOR, NU, J, ISYM, IFLAG, IERR, IWRIT,
     .           IR, IP, ISTS, IN, IY, IB, IT, IA, NRET, I, NSW, ISW,
     .           ISP, IHELP, K, IFL, ISYM_ERR, IRA, IRE, IPA, IPE,
     .           ITA, ITE, JJ
      LOGICAL :: PLSAV1, PLSAV2, LSTORE
      CHARACTER(20) :: TXTHST(NTXHST)
      CHARACTER(10) :: CX, CY, CX0, CY0, CZ0
      CHARACTER(4) :: CH

      SAVE
      DATA ABSMAX,ORDMAX/21.,21./
      DATA XNULL,YNULL/9.,4./,XWN,YWN/0.,0./
      DATA IWRIT/0/,ISPL/2,101,103,205,100,206,208,104,105,
     .                   106,107,108,200,201,202,204,207,4/
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
      IF (.NOT.ALLOCATED(ICPSPZ)) ALLOCATE (ICPSPZ(NSPZ))
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
      CALL GRNXTB(1,'PLT2D.F')
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
            IPE = MIN(IRPTE(J,2)-1,N2ND)
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
            ITE = MIN(IRPTE(J,3)-1,N3RD)
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
                WRITE (iunout,*) 'ERROR IN PLT2D '
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
!                WRITE (iunout,*) 'ERROR IN PLT2D '
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
              WRITE (iunout,*) 'PLOTOPTION CH2Z0.NE.0 NOT READY. EXIT '
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
                WRITE (iunout,*) 'ERROR IN PLT2D '
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
                WRITE (iunout,*) 'ERROR IN PLT2D '
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
                  WRITE (iunout,*) 'ERROR IN PLT2D '
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
      ELSEIF (LEVGEO.EQ.4) THEN
C
        IF (PLCUT(3)) THEN
          INSTOR=0
          DO 1010,I=1,NTRII
            P=0.
            DO 1020,J=1,3
              IECKE2 = J+1
              IF (IECKE2 .EQ. 4) IECKE2 = 1
C             SEITE J GEHOERT ZUM RAND
              ISTS=ABS(INMTI(J,I))
              IF (ISTS .GT. 0 .AND.
     .            ISTS .LE. NLIM+NSTSI) THEN
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
1020        CONTINUE
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
     .           goto 1010
              slt=slt*sclfcx
              call grchrc(REAL(slt,KIND(1.E0)),0.,idummy)
              CALL GRTXT (REAL(X,KIND(1.E0)),REAL(Y,KIND(1.E0)),4,CH)
              call grchrc(0.3,0.,idummy)
            ENDIF
1010      CONTINUE
          CALL GRDSH (1.,0.,1.)
          CALL GRNWPN (1)

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

          DO I=1,NTRII
            DO J=1,3
              IECKE2 = J+1
              IF (IECKE2 .EQ. 4) IECKE2 = 1
              XX(1)=XTRIAN(NECKE(J,I))
              YY(1)=YTRIAN(NECKE(J,I))
              XX(2)=XTRIAN(NECKE(IECKE2,I))
              YY(2)=YTRIAN(NECKE(IECKE2,I))

C  CHECK IF Y IS ENCLOSED BY YY(1) AND YY(2)
C  Y=CH2Z0 INTERSECTS WITH TRIANGLE
              IF ((MIN(YY(1),YY(2)) <= Y) .AND. 
     .            (MAX(YY(1),YY(2)) >= Y)) THEN  
                DXX = XX(1) - XX(2)  
                DYY = YY(1) - YY(2)  
C  FIND INTERSECTION POINT WITH TRIANGLE SIDE                
                IF (ABS(DXX) < EPS10) THEN
C  X=CONST
                  XW1 = XX(1)
                ELSE IF (ABS(DYY) < EPS10) THEN
C  Y=CONST
                  XW1 = XX(1)
                ELSE
                  A = DYY/DXX
                  B = YY(1) - A*XX(1)
                  XW1 = (Y-B)/A
                END IF
                ISTS=ABS(INMTI(J,I))
                IF (ISTS .GT. 0 .AND.
     .              ISTS .LE. NLIM+NSTSI) THEN
C   SEITE J HAT BESONDERE EIGENSCHAFTEN (NON-DEFAULT STD.FLAECHE)
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
C  NOW PLOT THAT RADIAL SURFACE
                IF (NLTRZ) THEN
                  IF (XW1.GE.XMI2D.AND.XW1.LE.XMA2D) THEN
                    CALL GRJMP (REAL(XW1,KIND(1.E0)),
     .                          REAL(YW1,KIND(1.E0)))
                    CALL GRDRW (REAL(XW1,KIND(1.E0)),
     .                          REAL(YW2,KIND(1.E0)))
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
                    WRITE (iunout,*) 'ERROR IN PLT2D '
                    CALL EXIT_OWN(1)
                  ENDIF
                  DO JJ=1,NTTRA
                    XX(JJ)=RR*COS(ZSURF(JJ))
                    YY(JJ)=RR*SIN(ZSURF(JJ))
                  ENDDO
                  CALL PLTLNE(NTTRA,XX,YY,XMI2D,XMA2D,
     .                        YMI2D,YMA2D,LSTORE)
                  CALL GRDSH (1.,0.,1.)
                  CALL GRNWPN (1)
C               ELSEIF (NLTRT) THEN
                ENDIF
              ENDIF
            END DO  ! J
          END DO ! I 
C
        END IF
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
        ICPSPZ=0
        DO 505 I=1,NPHOTI
          ISP=I
          IF (NHSTS(ISP).EQ.-1) GOTO 505
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          ICPSPZ(ISP)=ICP
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
505     CONTINUE
        DO 510 I=1,NATMI
          ISP=NSPH+I
          IF (NHSTS(ISP).EQ.-1) GOTO 510
          IC=IC+1
          ICP=MOD(IC,7)
          IF (ICP.EQ.0) ICP=7
          YNP=YNP-0.75/FY
          ICPSPZ(ISP)=ICP
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
          ICPSPZ(ISP)=ICP
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
          ICPSPZ(ISP)=ICP
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
          ICPSPZ(ISP)=ICP
          CALL GRNWPN (ICP)
          CALL GRJMP (REAL(XN1,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRDRW (REAL(XN2,KIND(1.E0)),REAL(YNP,KIND(1.E0)))
          CALL GRTXT (REAL(XN2,KIND(1.E0)),REAL(YNP-0.15,KIND(1.E0)),8,
     .                TEXTS(ISP))
513     CONTINUE
        IF (.NOT.NLPL3D)
     .     CALL GRSCLV(REAL(XMI2D,KIND(1.E0)),REAL(YMI2D,KIND(1.E0)),
     .                 REAL(XMA2D,KIND(1.E0)),REAL(YMA2D,KIND(1.E0)))
      ENDIF
C
C  WRITE TRACK DATA
C
      IF (TRCHST) THEN
        CALL LEER(1)
        WRITE (iunout,*) TXTHST(ISYM)
        IF (ISPZ.GT.0.AND.ISPZ.LE.NSPZ) THEN
          IF (ISYM.EQ.1.OR..NOT.NLTRC)  CALL MASJ1('NPANU   ',NPANU)
          WRITE (iunout,'(1X,A8)') TEXTS(ISPZ)
        ELSE
          WRITE (iunout,'(1X,A14,1X,I6)') 'LINE OF SIGHT ',NPANU
        ENDIF
        CALL MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
        IF (ITYP.EQ.3) THEN
          IF (LCART) THEN
          CALL MASR6 
     .         ('VELX,VELY,VELZ,VEL,E0,WEIGHT                     ',
     .           VELX,VELY,VELZ,VEL,E0,WEIGHT)
          ELSE
          CALL MASR6 
     .         ('VLXPAR,VLYPAR,VLZPAR,VELPAR,E0,WEIGHT            ',
     .           VLXPAR,VLYPAR,VLZPAR,VELPAR,E0,WEIGHT)
          ENDIF
        ELSE
          CALL MASR6 
     .         ('VELX,VELY,VELZ,VEL,E0,WEIGHT                     ',
     .                 VELX,VELY,VELZ,VEL,E0,WEIGHT)
        ENDIF
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
      ICOLOR=ICPSPZ(ISPZ)
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
        WRITE (iunout,*) 'ICOLOR,IN,ISYM,IFLAG,NHSTS ',
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
      WRITE (iunout,*) 'OPTION IN PLT2D NOT READY: X- OR RADIAL GRID'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'OPTION IN PLT2D NOT READY: Y- OR POLOIDAL GRID'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'OPTION IN PLT2D NOT READY: Z- OR TOROIDAL GRID'
      CALL EXIT_OWN(1)

C     following ENTRY is for reinitialization of EIRENE (DMH)
      
      ENTRY PLT2D_REINIT
      ABSMAX = 21.
      ORDMAX = 21
      XNULL = 9.
      YNULL = 4.
      XWN = 0.
      YWN =0.
      IWRIT = 0
      ISPL = (/2,101,103,205,100,206,208,104,105,
     .         106,107,108,200,201,202,204,207,4/)
      return 

      END
C ===== SOURCE: plt3d.f
C  3D GEOMETRY (AND TRAJECTORY) PLOT
C
      SUBROUTINE PLT3D (XR,YR,FAKX,FAKY,ITH,ABSMIN,ABSMAX,ORDMIN,ORDMAX)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CRECH
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CPL3D
      USE CPLOT
      USE CINIT
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CTRCEI
      USE CGEOM
      USE COMSOU
      USE CTEXT
      USE CLGIN
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: XR, YR, FAKX, FAKY,
     .                       ABSMIN, ABSMAX, ORDMIN, ORDMAX
      INTEGER, INTENT(OUT) :: ITH

      INTEGER, PARAMETER :: NPLY=501
      REAL(DP) :: AL(10), AR(10), XP(NPLY), YP(NPLY), ZPLOT(N3RD+NTOR),
     .          XSAVE(NPLY,N3RD+NTOR), YSAVE(NPLY,N3RD+NTOR),
     .          PHIAN(9),PHIEN(9)
      REAL(DP) :: XX(101),YY(101)
      REAL(DP) :: DM, RS, Y, TR, EP, EL, PHI, TA, TD, F1B, F2B, F3B,
     .          TB, RR, X, Z, CH2MXS, CH2MYS,
     .          CH2X0S, CH2Y0S, XCENT, YCENT, XNULL, YNULL, DELX, DELY,
     .          XH, YH, XMA3D, YMA3D, ZMA3D, XMI3D, YMI3D, ZMI3D, XN,
     .          YN, ZX0B, ZY0B, ZZ0B, PHA, PHE, CXB, B1B, B2B, B3B,
     .          CYB, CZB, F0B, B0B, RZYLB, CX, CY, CZ, ZX0, ZY0, ZZ0,
     .          RZYL, T1, T2, F0, F1, F2, F3
      REAL(SP) :: XPS(NPLY),YPS(NPLY)
      INTEGER:: ILT(N3RD+NTOR)
      INTEGER :: IR, IBR, IST, IS, NR, I1, IZ, NP, II, K, ID, NA, ISTP,
     .           IA, IAN, IEN, KIN, ISSTD, IBA, NJZ, J, JJ, IPZ,
     .           IP, I, NINNE, NZAD, NIN, MERK2, IPR, IB, MERK,
     .           IJZ, JP, IPRT
      LOGICAL :: PLABLE(NLIM), LPERID(NLIM), LSYMET(NLIM),
     .           LERR1, LERR2, LSAVE, PLT1, PLT2, PLT3
      TYPE(PPOINT), POINTER :: CUR
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
        ITH=1
        WIN=0.
        RMT=0.
        WINJ=WIN
      ENDIF
C
      XP=0.D0
      YP=0.D0
      PLABLE(1:NLIM)=.FALSE.
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
      CALL GRNXTB(1,'PLT3D.F')
      CALL GRSCLC(REAL(XNULL,KIND(1.E0)),REAL(YNULL,KIND(1.E0)),
     .            REAL(XN+XNULL,KIND(1.E0)),
     .            REAL(YN+YNULL,KIND(1.E0)))
      CALL GRSCLV(REAL(ABSMIN,KIND(1.E0)),REAL(ORDMIN,KIND(1.E0)),
     .            REAL(ABSMAX,KIND(1.E0)),REAL(ORDMAX,KIND(1.E0)))
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
              WRITE (iunout,*) 'SURFACE NO. ',J,' OUT'
            ENDIF
            GOTO 10
          ELSE
            IF (NLTRA.AND.ILTOR(J).LE.0) THEN
C  STILL TO BE WRITTEN: BETTER WAY OF IDENTIFYING TOROIDALLY SYMMETRIC
C                       SURFACE
              IF (P3(3,J).GT.1.D50.AND.ZLIMS1(1,J).LT.-1.D15.AND.
     .                                 ZLIMS2(1,J).GT.1.D15) THEN
                LSYMET(J)=.TRUE.
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SURFACE NO. ',J,
     .              ' TOROIDALLY SYMMETRIC'
                  WRITE (iunout,*) 'PLOT LATER INTO STANDARD MESH '
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
                  WRITE (iunout,*) 'SURFACE NO. ',J,
     .              ' TOROIDALLY PERODIC'
                  WRITE (iunout,*) 'PLOT ADD. SURFACE NO. ',J
                  WRITE (iunout,*) 'INTO ',NJZ, ' TOROIDAL SEGMENTS'
                ENDIF
              ENDIF
            ELSEIF (.NOT.NLTRA.OR.ILTOR(J).GT.0) THEN
              NJZ=1
              ILT(1)=ILTOR(J)
              LSYMET(J)=.FALSE.
              LPERID(J)=.FALSE.
              IF (TRCPLT) THEN
                WRITE (iunout,*) 'PLOT ADD. SURFACE NO. ',J
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
                WRITE (iunout,*) 'FL2O CALLED '
                WRITE (iunout,*) 'MERK= ',MERK
              ENDIF
C**ZYLINDER: FINDE ACHSE
              IF (MERK.EQ.4) THEN
                IF (TRCPLT) WRITE (iunout,*) 'CX,CY,CZ,RZYL ',
     .                                        CX,CY,CZ,RZYL
                IF (ABS(CX).GT.1.D-6.AND.
     .              MAX(ABS(CY),ABS(CZ)).LE.1.D-6) THEN
                  T1=(XLIMS1(1,J)-ZX0)/CX
                  T2=(XLIMS2(1,J)-ZX0)/CX
                ELSEIF (ABS(CY).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CZ)).LE.1.D-6) THEN
                  T1=(YLIMS1(1,J)-ZY0)/CY
                  T2=(YLIMS2(1,J)-ZY0)/CY
                ELSEIF (ABS(CZ).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CY)).LE.1.D-6) THEN
                  T1=(ZLIMS1(1,J)-ZZ0)/CZ
                  T2=(ZLIMS2(1,J)-ZZ0)/CZ
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
C  ZYLINDER: GGFLS MEHRERE TEILSTUECKE
                CALL CTQUA (A0LM(J),A1LM(J),A2LM(J),A3LM(J),A4LM(J),
     .                      A5LM(J),A6LM(J),A7LM(J),A8LM(J),A9LM(J),
     .                      XLIMS1(1,J),XLIMS2(1,J),YLIMS1(1,J),
     .                      YLIMS2(1,J),ZLIMS1(1,J),ZLIMS2(1,J),
     .                      RLB(J),RZYL,ZX0,ZY0,ZZ0,PHIAN,PHIEN,IPRT)
                IF (TRCPLT) THEN
                  WRITE (iunout,*) ' T1,T2 ',T1,T2
                  WRITE (iunout,*) ' IPRT ',IPRT
                  WRITE (iunout,*) ' PHIAN,PHIEN ',(PHIAN(JP),PHIEN(JP),
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
                IF (TRCPLT) WRITE (iunout,*) 'CX,CY,CZ,RZYL ',
     .                                        CX,CY,CZ,RZYL
                IF (ABS(CX).GT.1.D-6.AND.
     .              MAX(ABS(CY),ABS(CZ)).LE.1.D-6) THEN
                  T1=(XLIMS1(1,J)-ZX0)/CX
                  T2=(XLIMS2(1,J)-ZX0)/CX
                ELSEIF (ABS(CY).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CZ)).LE.1.D-6) THEN
                  T1=(YLIMS1(1,J)-ZY0)/CY
                  T2=(YLIMS2(1,J)-ZY0)/CY
                ELSEIF (ABS(CZ).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CY)).LE.1.D-6) THEN
                  T1=(ZLIMS1(1,J)-ZZ0)/CZ
                  T2=(ZLIMS2(1,J)-ZZ0)/CZ
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
                CALL ELLIPSOID (ZX0,ZY0,ZZ0,CX,CY,CZ,XLIMS1(1,J),
     .               YLIMS1(1,J),ZLIMS1(1,J),XLIMS2(1,J),YLIMS2(1,J),
     .               ZLIMS2(1,J),RLB(J),ILCOL(J),5,5,5)
C**PAAR VON EBENEN (ODER EINE DOPPELEBENE)
              ELSEIF (MERK.EQ.1.OR.MERK.EQ.2.OR.MERK.EQ.3) THEN
                IF (TRCPLT) WRITE (iunout,*) 'B0,B1,B2,B3 ',B0,B1,B2,B3
                CALL PLANE (B0,B1,B2,B3,RLB(J),9,EPS10,
     .                      ALIMS,XLIMS,YLIMS,ZLIMS,
     .                      ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                             XLIMS2,YLIMS2,ZLIMS2,
     .                             XLIMS3,YLIMS3,ZLIMS3,
     .                      ILCOL(J),IGFIL(J).NE.0,J)
                IF (MERK.NE.1) THEN
                  IF (TRCPLT) WRITE (iunout,*) 'F0,F1,F2,F3 ',
     .                                          F0,F1,F2,F3
                  CALL PLANE (F0,F1,F2,F3,RLB(J),9,EPS10,
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
              CALL PLANE (A0LM(J),A1LM(J),A2LM(J),A3LM(J),RLB(J),9,
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
                WRITE (iunout,*) 'FL2O CALLED '
                WRITE (iunout,*) 'MERK= ',MERK
              ENDIF
C
C**PAAR VON EBENEN ODER DOPPELEBENE ?
              IF (MERK.LE.3) THEN
                IF (ISCN(J).EQ.0) THEN
                  CALL PLANE (B0,B1,B2,B3,RLB(J),
     .                        9,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                  CALL PLANE (F0,F1,F2,F3,RLB(J),
     .                        9,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
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
                 CALL FL2O (ALIMS0(IB,J),XLIMS1(IB,J),YLIMS1(IB,J),
     .                      ZLIMS1(IB,J),XLIMS2(IB,J),YLIMS2(IB,J),
     .                      ZLIMS2(IB,J),XLIMS3(IB,J),YLIMS3(IB,J),
     .                      ZLIMS3(IB,J),MERK2,ZX0B,ZY0B,ZZ0B,CXB,CYB,
     .                      CZB,RZYLB,B0B,B1B,B2B,B3B,F0B,F1B,F2B,F3B,
     .                      EPS10,NMACH)
                 IF (TRCPLT) THEN
                   WRITE (iunout,*) 'FL2O CALLED '
                   WRITE (iunout,*) 'MERK2= ',MERK2
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
                     IF (TRCPLT) WRITE (iunout,*)
     .                  ' FEHLER IN BERANDUNG VON FLAECHE ',J
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
     .                          0._DP,360._DP)
C**ZYLINDER BEGRENZT VON ECHT GEKRUEMMTEN FLAECHE 2TER ORDNUNG
                ELSEIF (MERK2.GE.4) THEN
                  IB=1
                  AL(1)=ALIMS0(IB,J)
                  AL(2)=XLIMS1(IB,J)
                  AL(3)=YLIMS1(IB,J)
                  AL(4)=ZLIMS1(IB,J)
                  AL(5)=XLIMS2(IB,J)
                  AL(6)=YLIMS2(IB,J)
                  AL(7)=ZLIMS2(IB,J)
                  AL(8)=XLIMS3(IB,J)
                  AL(9)=YLIMS3(IB,J)
                  AL(10)=ZLIMS3(IB,J)
                  CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,10,TA,TB,LERR1)
                  IF (LERR1) THEN
                    WRITE (iunout,*)
     .                ' FEHLER IN DER BERANDUNG VON FLAECHE',J
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
     .                      J,10,AL,10,AR,0._DP,360._DP)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
               ENDIF
C**KUGEL, ELLIPSOID
              ELSEIF (MERK == 13) THEN
                CALL ELLIPSOID (ZX0,ZY0,ZZ0,CX,CY,CZ,XLIMS1(1,J),
     .               YLIMS1(1,J),ZLIMS1(1,J),XLIMS2(1,J),YLIMS2(1,J),
     .               ZLIMS2(1,J),RLB(J),ILCOL(J),5,5,5)
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
     .                      9,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                      ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                             XLIMS2,YLIMS2,ZLIMS2,
     .                             XLIMS3,YLIMS3,ZLIMS3,
     .                      ILCOL(J),IGFIL(J).NE.0,J)
              ELSEIF (ILIN(J).EQ.0) THEN
C**EBENE BEGRENZT DURCH EINEN ODER MEHRERE ZYLINDER?
                IB=0
20              IB=IB+1
                CALL FL2O (ALIMS0(IB,J),XLIMS1(IB,J),YLIMS1(IB,J),
     .                     ZLIMS1(IB,J),XLIMS2(IB,J),YLIMS2(IB,J),
     .                     ZLIMS2(IB,J),XLIMS3(IB,J),YLIMS3(IB,J),
     .                     ZLIMS3(IB,J),MERK2,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                     RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'FL2O CALLED '
                  WRITE (iunout,*) 'MERK2= ',MERK2
                ENDIF
                IF (MERK2.EQ.4) THEN
                  AL(1)=A0LM(J)
                  AL(2)=A1LM(J)
                  AL(3)=A2LM(J)
                  AL(4)=A3LM(J)
                  CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,4,TA,TD,LERR1)
                  IF (LERR1) THEN
                    WRITE (iunout,*)
     .                ' FEHLER IN DER BERANDUNG VON FLAECHE',J
                    PLABLE(J)=.TRUE.
                    GOTO 10
                  ENDIF
                  T1=TA-2.*RZYL
                  T2=TA+4.*RZYL
                  NP=1
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                         RZYL,NP,NINNE,NIN,
     .                         ILCOL(J),IGFIL(J).NE.0,J,4,AL,0,AR,
     .                         0._DP,360._DP)
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
            WRITE (iunout,*) 'MESSAGE FROM SUBR. PLT3D :'
            WRITE (iunout,*) 'SURFACE NO. J COULD NOT BE PLOTTED'
            I1=1
          ENDIF
          WRITE (iunout,*) 'J= ',J
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
              IF (TRCPLT) WRITE (iunout,*) 
     .          'PLOT RAD. STAND. SURFACE NO. ',IR
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
                      WRITE (iunout,*) 'FROM PLT3D: NOT ENOUGH STORAGE '
                      WRITE (iunout,*) 'INCREASE PARAMETER NPLY '
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
                WRITE(iunout,*) ' STORAGE EXCEEDED IN PLT3D, IR= ',IR
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
                CUR => FIRST_POINT
                DO WHILE(ASSOCIATED(CUR))
                  NA=NA+1
                  X=CUR%XPL2D
                  Y=CUR%YPL2D
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
                  IF (CUR%NPL2D.EQ.0) CALL GRJMP (
     .               REAL(XP(NA),KIND(1.E0)),
     .               REAL(YP(NA),KIND(1.E0)))
                  IF (CUR%NPL2D.EQ.1) CALL GRDRW (
     .               REAL(XP(NA),KIND(1.E0)),
     .               REAL(YP(NA),KIND(1.E0)))
C
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(NA)
                  YSAVE(IS,IZ)=YP(NA)
                  DO ID=1,ISTP
                    IF (ASSOCIATED(CUR)) CUR => CUR%NXTPNT
                  END DO
                END DO
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
          CALL GRJMP (REAL(XSAVE(J,1),KIND(1.E0)),
     .                REAL(YSAVE(J,1),KIND(1.E0)))
          DO 2000 IZ=2,NJZ
            CALL GRDRW(REAL(XSAVE(J,IZ),KIND(1.E0)),
     .                 REAL(YSAVE(J,IZ),KIND(1.E0)))
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
      CALL GRTXT (REAL(XH,KIND(1.E0)),REAL(YH,KIND(1.E0)),27,
     .            'CHECK OF GEOMETRICAL INPUT:')
      YH=YH-0.5/FAKY
      CALL GRTXT (REAL(XH,KIND(1.E0)),REAL(YH,KIND(1.E0)),72,TXTRUN)
      YH=YH-0.75/FAKY
      DO 11000 J=1,5
        IF (PL3A(J)) CALL GRTXT(REAL(XH,KIND(1.E0)),REAL(YH,KIND(1.E0)),
     .                          16,TEXTLA(J))
        IF (PL3A(J)) YH=YH-0.5/FAKY
11000 CONTINUE
C
C  RETURN CO-ORDINATES FOR PLOTS OF PARTICLE TRACKS
      XR=ABSMIN-6./FAKX
      YR=ORDMAX-0./FAKY
      RETURN
      END
C ===== SOURCE: pltadd.f
C
C
      SUBROUTINE PLTADD (MANF,MEND)
C
C  PLOT ADDITIONAL SURFACES ISURF=MANF,MEND, PROJECTED INTO A PLANE
C  IF LZR,     PLOT IS DONE IN THIS ROUTINE
C  IF.NOT.LZR, PLOT DATA ARE COLLECTED IN MODULE CRECH
C
      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CADGEO
      USE CCONA
      USE CLMSUR
      USE CPLOT
      USE CTRCEI
      USE CTEXT
      USE CLGIN
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MANF, MEND
      REAL(DP) :: G(4,2), R(4,2), AFF(3,3), AFFI(3,3)
      REAL(DP) :: ELLO, ELLU, PARA1, PARA2O, PARA2U,
     .          HYPP, HYPM, GERADY, GERADX
      REAL(DP) :: XTRAN, YTRAN, XI, ETA, XA, XE, XN, X1, X2,
     .          XANF, XEND, RAD, X13, Y13, X14, Y14, X23,
     .          Y23, X24, Y24, YANF, YEND, ARC05, XINI, XEN,
     .          ARC, YINI, XAN, XP, YP, X, Y, XS, YS, SQ, SQR,
     .          VR1, VR2, VR3, AMNPY, AMNPZ, AMXPY, AMXPZ, XMU,
     .          FMU, XP1, V1, V2, V3, XCHL, XCHR, YCHL, YCHR,
     .          DCHX, AMXPX, AMNPX, DCHY, YP1, TANA, SINAQ,
     .          COSAQ, ALF, DELS, DET, DETS, DETER,
     .          SCA, DH, EH, CK, DK, EK, BK, XP2, YP2, AK, AKEK2,
     .          S, DEL, AKEK, FK, DKCK, DKCK2
      INTEGER :: INN1, INN2, INN3, I, INN4, INN, IS, ICOLOR, J,
     .           ICUT, ISURF
      LOGICAL :: L1, L2, L3, L4, L5
      LOGICAL :: LBOX, LLISTE
      CHARACTER(1) :: CH1
      CHARACTER(2) :: CH2
      CHARACTER(3) :: CH3
      TYPE(PPOINT), POINTER :: CUR, SURFAN, SURFEN

      EXTERNAL ELLO,ELLU,PARA1,PARA2O,PARA2U,HYPP,HYPM,
     .         GERADY,GERADX

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
C  PLOTTING PLANE IS X-Y, AT Z=CH2Z0
      DO 2000 ICUT=1,3
      IF (.NOT.PLCUT(ICUT)) GOTO 2000
      IF (ICUT.EQ.1) THEN
C  PLOT Y-Z AT X=CH2Z0. HENCE MAP:
C  Y --> Y
C  Z --> X
C  X --> Z
        IF (CH2Z0.NE.0.) CALL XSHADD(-CH2Z0,MANF,MEND)
        AFF=0.D0
        AFF(1,3)=1.
        AFF(2,2)=1.
        AFF(3,1)=-1.
        AFFI=0.D0
        AFFI(1,3)=-1.
        AFFI(2,2)=1.
        AFFI(3,1)=1.
        CALL ROTADD(AFF,AFFI,MANF,MEND)
      ELSEIF (ICUT.EQ.2) THEN
C  PLOT X-Z AT Y=CH2Z0. HENCE MAP:
C  X --> X
C  Z --> Y
C  Y --> -Z
        IF (CH2Z0.NE.0.) CALL YSHADD(-CH2Z0,MANF,MEND)
        AFF=0.D0
        AFF(1,1)=1.
        AFF(2,3)=1.
        AFF(3,2)=-1.
        AFFI=0.D0
        AFFI(1,1)=1.
        AFFI(2,3)=-1.
        AFFI(3,2)=1.
        CALL ROTADD(AFF,AFFI,MANF,MEND)
      ELSEIF (ICUT.EQ.3) THEN
C  PLOT X-Y AT Z=CH2Z0. HENCE MAP:
        IF (CH2Z0.NE.0.) CALL ZSHADD(-CH2Z0,MANF,MEND)
C  NO ROTATION NEEDED
      ENDIF
C
C
C   PLOT SURFACES
      ZPLT=0.
      INSTOR=0
C
      DO 200 ISURF=MANF,MEND
        IF (ILCOL(ISURF) == 666) GOTO 200
        J=ISURF
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) INOSF=ISURF
C
C  PLOT SURFACE NO. ISURF
C
        IF (IGJUM0(ISURF).NE.0) GOTO 200
        IF (TRCPLT) THEN
          CALL LEER(1)
          WRITE (iunout,*) TXTSFL(ISURF)
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
            IF (TRCPLT) WRITE (iunout,6660)
            GOTO 200
          ENDIF
          IF (TRCPLT) THEN
            WRITE (iunout,*) ' AMNPX,AMNPY,AMNPZ ',AMNPX,AMNPY,AMNPZ
            WRITE (iunout,*) ' AMXPX,AMXPY,AMXPZ ',AMXPX,AMXPY,AMXPZ
            WRITE (iunout,*) ' A0,A1,A2,A3 ',
     .                   A0LM(ISURF),A1LM(ISURF),A2LM(ISURF),A3LM(ISURF)
          ENDIF
c
          IF (ABS(A3LM(ISURF)).GE.1.-EPS10) THEN
C   N ECK PARALLEL ZU Z=ZPLT
            IF (TRCPLT) WRITE (iunout,*) 'N ECK PARALLEL ZU Z=ZPLT'
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
            IF (TRCPLT) WRITE (iunout,*) 'L1,L2,L3,L4,L5 ',
     .                                    L1,L2,L3,L4,L5
C   PLOT P1, P2
            IF (L1.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P1(1,ISURF),KIND(1.E0)),
     .                      REAL(P1(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P2(1,ISURF),KIND(1.E0)),
     .                      REAL(P2(2,ISURF),KIND(1.E0)))
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
                CALL GRJMP (REAL(P1(1,ISURF),KIND(1.E0)),
     .                      REAL(P1(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
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
                CALL GRJMP (REAL(P2(1,ISURF),KIND(1.E0)),
     .                      REAL(P2(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
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
                CALL GRJMP (REAL(P4(1,ISURF),KIND(1.E0)),
     .                      REAL(P4(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P2(1,ISURF),KIND(1.E0)),
     .                      REAL(P2(2,ISURF),KIND(1.E0)))
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
                CALL GRJMP (REAL(P4(1,ISURF),KIND(1.E0)),
     .                      REAL(P4(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
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
                CALL GRJMP (REAL(P4(1,ISURF),KIND(1.E0)),
     .                      REAL(P4(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P5(1,ISURF),KIND(1.E0)),
     .                      REAL(P5(2,ISURF),KIND(1.E0)))
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
                CALL GRJMP (REAL(P5(1,ISURF),KIND(1.E0)),
     .                      REAL(P5(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
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
            IF (TRCPLT) WRITE (iunout,*) ' VR1,VR2,VR3 ',VR1,VR2,VR3
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
              WRITE (iunout,*) ' SCHNITTGERADE AUSSERHALB DES VIERECKS'
            ELSEIF (IS.EQ.1.AND.
     .              XP1.GE.XCHL.AND.XP1.LE.XCHR.AND.YP1.GE.YCHL.AND.
     .              YP1.LE.YCHR) THEN
              IF (LZR) THEN
                CALL GRSPTS (30)
                CALL GRJMP (REAL(XP1,KIND(1.E0)),REAL(YP1,KIND(1.E0)))
                CALL GRDRW (REAL(XP1,KIND(1.E0)),REAL(YP1,KIND(1.E0)))
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
                CALL GRJMP (REAL(XP1,KIND(1.E0)),REAL(YP1,KIND(1.E0)))
                CALL GRDRW (REAL(XP2,KIND(1.E0)),REAL(YP2,KIND(1.E0)))
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
          IF (ZLIMS1(1,ISURF).GT.ZPLT.OR.ZLIMS2(1,ISURF).LT.ZPLT) THEN
            IF (TRCPLT) THEN
              WRITE (iunout,*) 'ZLIMS1,ZPLT,ZLIMS2 ',
     .                     ZLIMS1(1,ISURF),ZPLT,ZLIMS2(1,ISURF)
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
          WRITE (iunout,*) 'COEFFICIENTS OF CURVE-EQUATION'
          WRITE (iunout,*) 'AK*X**2+BK*XY+CK*Y**2+DK*X+EK*Y+FK=0.    '
          WRITE (iunout,6662)  AK,BK,CK,DK,EK,FK
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
          WRITE (iunout,*) 'S,DELS,DETS,DEL,DET '
          WRITE (iunout,*) S,DELS,DETS,DEL,DET
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
          IF (SIGN(1._DP,SINA).NE.SIGN(1._DP,BK)) THEN
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
          DET=DETER(A,0._DP,DH,0._DP,C,EH,DH,EH,F)
          IF (ABS(DET).LE.EPS10) DET=0.
        ENDIF
C
        IF (TRCPLT) THEN
          WRITE (iunout,*) 
     .      'AFTER TRANSFORMATION: B=0., NEW COEFFICIENTS'
          WRITE (iunout,6661) A,C,D,E,F
          WRITE (iunout,*) 'DETERMINANT DET= ',DET
        ENDIF
C
        IF (RLB(ISURF).EQ.1.) THEN
          XL1=MAX(XCHL,XLIMS1(1,ISURF))
          XL1=MIN(XCHR,XL1)
          XL2=MAX(XCHL,XLIMS2(1,ISURF))
          XL2=MIN(XCHR,XL2)
          YL1=MAX(YCHL,YLIMS1(1,ISURF))
          YL1=MIN(YCHR,YL1)
          YL2=MAX(YCHL,YLIMS2(1,ISURF))
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
          XC1=XLIMS1(1,ISURF)
          XC2=XLIMS2(1,ISURF)
          YC1=YLIMS1(1,ISURF)
          YC2=YLIMS2(1,ISURF)
        ENDIF
C
        IF (RLB(ISURF).LT.0.) THEN
           DO 401 I=1,ILIN(ISURF)
              ALIN(I)=ALIMS(I,ISURF)
              XLIN(I)=XLIMS(I,ISURF)
              YLIN(I)=YLIMS(I,ISURF)
              ZLIN(I)=ZLIMS(I,ISURF)
401        CONTINUE
           DO 402 I=1,ISCN(ISURF)
              A0S(I)=ALIMS0(I,ISURF)
              A1S(I)=XLIMS1(I,ISURF)
              A2S(I)=YLIMS1(I,ISURF)
              A3S(I)=ZLIMS1(I,ISURF)
              A4S(I)=XLIMS2(I,ISURF)
              A5S(I)=YLIMS2(I,ISURF)
              A6S(I)=ZLIMS2(I,ISURF)
              A7S(I)=XLIMS3(I,ISURF)
              A8S(I)=YLIMS3(I,ISURF)
              A9S(I)=ZLIMS3(I,ISURF)
402        CONTINUE
           MLIN=ILIN(ISURF)
           MSCN=ISCN(ISURF)
        ENDIF
        IF (TRCPLT) THEN
          WRITE (iunout,*) 'PLOT REGION:'
          WRITE (iunout,*) 'XL1,XL2,YL1,YL2 ',XL1,XL2,YL1,YL2
          WRITE (iunout,*) 'PLOT REGION AFTER TRANSFORMATION'
          WRITE (iunout,*) 'XANF,XEND,YANF,YEND ',XANF,XEND,YANF,YEND
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
407     IF (TRCPLT) WRITE (iunout,6664)
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
           XA=MAX(XM-AHALB,XANF+1.E-8_DP)-XM
           XE=MIN(XM+AHALB,XEND-1.E-8_DP)-XM
        ENDIF
        CALL PLTIN (RLB(ISURF),ELLO,XA,XE,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'ELLO'
        CALL PLTIN (RLB(ISURF),ELLU,XA,XE,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'ELLU'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6665)
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6665)
        GOTO 200
C
C     HYPERBEL
C
450     CONTINUE
        RAD=MAX(0._DP,XN/A)
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
        IF (TRCPLT) WRITE (iunout,*) 'HYPP'
        CALL PLTIN (RLB(ISURF),HYPM,XA,XE,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPM'
        IF (INN1+INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6666)
451     IF (X2.LT.XL1) GOTO  452
        XA=XANF
        XE=MIN(X2,XEND)-1.D-6
        CALL PLTIN (RLB(ISURF),HYPP,XA,XE,INN3,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPP'
        CALL PLTIN (RLB(ISURF),HYPM,XA,XE,INN4,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPM'
452     IF (INN3+INN4.EQ.0.AND.TRCPLT) WRITE (iunout,6666)
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
        IF (TRCPLT) WRITE (iunout,*) 'POINT'
        IF (XP.LT.XL1.OR.XP.GT.XL2.OR.YP.LT.YL1.OR.YP.GT.YL2) THEN
          IF (TRCPLT) WRITE (iunout,6667)
          GOTO 200
        ENDIF
        IF (LZR) THEN
          CALL GRCHRC (0.1,0.0,16)
          CALL GRJMPS (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)),3)
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
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        A0=SQRT(ABS(A/C))
        A1=SQRT(ABS(A/C))*D/(2.*A)-E/(2.*C)
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        GOTO 200
C
1000    CONTINUE
        IF (ABS(A).GT.EPS10) GOTO 1300
        IF (ABS(C).GT.EPS10) GOTO 1400
        IF (ABS(E).GT.EPS10.AND.ABS(E).GE.ABS(D)) GOTO 1200
        IF (ABS(D).GT.EPS10.AND.ABS(D).GE.ABS(E)) GOTO 1100
        IF (ABS(F).GT.EPS10) GOTO 1002
        IF (TRCPLT) THEN
        WRITE (iunout,*) 'GLEICHUNG DER FORM 0.=0.'
        ENDIF
        GOTO 200
1002    IF (TRCPLT) WRITE (iunout,1001) F
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
        IF (TRCPLT) WRITE (iunout,*) 'GERADX'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (iunout,6669)
        GOTO 200
C
C     GERADE; D*X + E*Y + F = 0, E UNGLEICH 0., D=0. MOEGLICH
C
1200    CONTINUE
        A0=-D/E
        A1=-F/E
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (iunout,6669)
        GOTO 200
C
C     PARABEL; A*X**2 + D*X + E*Y + F = 0
C
1300    CONTINUE
        IF (ABS(DET).LE.EPS10) GOTO 1310
        CALL PLTIN (RLB(ISURF),PARA1,XANF,XEND,INN,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'PARA1'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (iunout,6668)
        GOTO 200
C
C     PAAR PARALLELER GERADEN; A*X**2 + D*X + F = 0
C
1310    SQ=(D*D/(4.*A)-F)/A
        SQR=0.
        IF (SQ.GT.0.) SQR=SQRT(SQ)
        IF (SQ.LT.0.) THEN
          IF (TRCPLT) WRITE (iunout,6664)
          GOTO 200
        ENDIF
        A0=0.
        A1=-D/(2.*A)+SQR
        CALL PLTIN (RLB(ISURF),GERADX,YANF,YEND,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADX'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        A1=-D/(2.*A)-SQR
        CALL PLTIN (RLB(ISURF),GERADX,YANF,YEND,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADX'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
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
           IF (TRCPLT) WRITE (iunout,6668)
           GOTO 200
        ENDIF
        CALL PLTIN (RLB(ISURF),PARA2O,XAN,XEN,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'PARA2O'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6671)
        CALL PLTIN (RLB(ISURF),PARA2U,XAN,XEN,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'PARA2U'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6671)
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
           IF (TRCPLT) WRITE (iunout,6664)
           GOTO 200
        ENDIF
        A1=-E/(2.*C)+SQR
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        A1=-E/(2.*C)-SQR
        CALL PLTIN (RLB(ISURF),GERADY,XANF,XEND,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
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
      IF (ICUT.EQ.1) THEN
        CALL ROTADD(AFFI,AFF,MANF,MEND)
        IF (CH2Z0.NE.0.) CALL XSHADD(CH2Z0,MANF,MEND)
      ELSEIF (ICUT.EQ.2) THEN
        CALL ROTADD(AFFI,AFF,MANF,MEND)
        IF (CH2Z0.NE.0.) CALL YSHADD(CH2Z0,MANF,MEND)
      ELSEIF (ICUT.EQ.3) THEN
        IF (CH2Z0.NE.0.) CALL ZSHADD(CH2Z0,MANF,MEND)
C  NO INVERS ROTATION NEEDED
      ENDIF
C
2000  CONTINUE
C
      IF (TRCPLT) WRITE (iunout,*) 'INSTOR= ',INSTOR
C
      IF (PLNUMS.OR.PLARR) THEN
        CUR => FIRST_POINT
        DO WHILE (ASSOCIATED(CUR))
          IF (CUR%NPL2D == 0) SURFAN => CUR
          ARC=0.
          lliste = ASSOCIATED(CUR%NXTPNT)
          if (lliste) lliste = lliste .and. (CUR%NXTPNT%NPL2D == 1)
          DO WHILE (lliste)
            ARC=ARC+SQRT((CUR%XPL2D - CUR%NXTPNT%XPL2D)**2 +
     .                   (CUR%YPL2D - CUR%NXTPNT%YPL2D)**2)
            CUR => CUR%NXTPNT
            lliste = lliste .and. ASSOCIATED(CUR%NXTPNT)
            if (lliste) lliste = lliste .and. (CUR%NXTPNT%NPL2D == 1)
          END DO
          SURFEN => CUR
          ARC05=ARC*0.5
          ARC=0.
          CUR => SURFAN
          DO WHILE (ARC < ARC05)
            ARC=ARC+SQRT((CUR%XPL2D - CUR%NXTPNT%XPL2D)**2 +
     .                   (CUR%YPL2D - CUR%NXTPNT%YPL2D)**2)
            IF (ARC < ARC05) CUR => CUR%NXTPNT
          END DO
C  POINT BETWEEN CUR AND CUR%NXTPNT
          XINI=(CUR%XPL2D+CUR%NXTPNT%XPL2D)*0.5
          YINI=(CUR%YPL2D+CUR%NXTPNT%YPL2D)*0.5
C
C  PLOT ARROWS: SURFACE NORMAL
          IF (PLARR) THEN
c           xlst=
c           ylst=
c           alen=
c           awid=
c           icode=
c           call grarrw(REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),xlst,ylst,
c    .                  REAL(alen,KIND(1.E0)),REAL(awid,KIND(1.E0)),icode)
c
          ENDIF
C  PLOT SURFACE NUMBERS
          IF (PLNUMS) THEN
            IF (CUR%NUMSUR.LT.10) THEN
              WRITE (CH1,'(I1)') CUR%NUMSUR
              CALL GRTXT (REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),
     .                    1,CH1)
            ELSEIF (CUR%NUMSUR.LT.100) THEN
              WRITE (CH2,'(I2)') CUR%NUMSUR
              CALL GRTXT (REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),
     .                    2,CH2)
            ELSEIF (CUR%NUMSUR.LT.1000) THEN
              WRITE (CH3,'(I3)') CUR%NUMSUR
              CALL GRTXT (REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),
     .                    3,CH3)
            ENDIF
          ENDIF
          CUR => SURFEN%NXTPNT
        END DO
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
C ===== SOURCE: pltaxi.f
C
C
C----------------------------------------------------------------------*
C SUBROUTINE PLTAXI                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE PLTAXI(IERR,IAX)
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
      USE CPLMSK

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IAX
      INTEGER, INTENT(OUT) :: IERR
      REAL(DP) :: EXPR, RI, ARG, RJ, T, PARAM1, PARAM,CENTRAL
      INTEGER :: FCTR, IL, IPARAM, I, J, IEXP10
      CHARACTER(7) :: CPARAM
      CHARACTER(3) :: CFCTR, CZHNLB
      CHARACTER(1) :: CEINLB
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
            CALL GRSCLV(REAL(MINX,KIND(1.E0)),0.,REAL(MAXX,KIND(1.E0)),
     .                  REAL(LENY,KIND(1.E0)))
C
C  GITTERLINIEN BZW. MARKIERUNGEN AN DER X-ACHSE
C
            DO 5 J=0,INTNRX
               T=MINX+J*STPSZX
               CALL GRJMP(REAL(T,KIND(1.E0)),-0.1)
               IF (GRIDX) THEN
                  IF (J.NE.0.AND.J.NE.INTNRX) THEN
                    CALL GRSPTS(16)
                    CALL GRDSH(0.2,0.5,0.2)
                  ENDIF
                  CALL GRDRW(REAL(T,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
                  CALL GRSPTS(20)
                  CALL GRDSH(1.,0.,1.)
               ELSE
                  CALL GRDRW(REAL(T,KIND(1.E0)),0.)
               ENDIF
    5       CONTINUE
            IF (.NOT.GRIDX) THEN
               CALL GRJMP(REAL(MINX,KIND(1.E0)),0.)
               CALL GRDRW(REAL(MINX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
               CALL GRJMP(REAL(MAXX,KIND(1.E0)),0.)
               CALL GRDRW(REAL(MAXX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
            ENDIF
C
C  LABELS AN DER X-ACHSE
C
            CALL GRCHRC(0.3,0.,16)
            DO 10 J=2,INTNRX+1,2
               CENTRAL=REAL((MINX+MAXX)*0.5,DP)
               IF (ABS(CENTRAL).LT.1.E-30)
     .           CENTRAL=MAX(ABS(MINX),ABS(MAXX))
               FCTR=IEXP10(CENTRAL)
               PARAM=(MINX+(J-1)*STPSZX)/10.**FCTR
               PARAM1=PARAM*100000.
               IPARAM=NINT(PARAM1)
               PARAM=REAL(IPARAM,KIND(1.E0))/100000.
               if (Abs(param).lt.100.) then
                 WRITE(CPARAM,'(F7.3)') PARAM
               elseif (Abs(param).lt.1000.) then
                 WRITE(CPARAM,'(F7.2)') PARAM
               else
                 WRITE(CPARAM,'(F7.1)') PARAM
               endif
               IL=7
               T=MINX+(J-1)*STPSZX-(MAXX-MINX)*0.8/LENX
               CALL GRTXT(REAL(T,KIND(1.E0)),-0.5,IL,CPARAM)
   10       CONTINUE
            T=MAXX-1.6*(MAXX-MINX)/LENX
            IF (ABS(FCTR).GE.10) THEN
              WRITE(CFCTR,'(I3)') FCTR
            ELSE
              WRITE(CFCTR,'(I2)') FCTR
              CFCTR(3:3)=' '
            ENDIF
            IF (FCTR.LT.0.) CALL GRTXT(REAL(T,KIND(1.E0)),
     .                                 -1.0,8,'*10**'//CFCTR)
            IF (FCTR.GT.0.) CALL GRTXT(REAL(T,KIND(1.E0)),-1.0,7,
     .                                 '*10**'//CFCTR(2:3))
C
C  LOGARITHMISCHE X-ACHSENEINTEILUNG
C
          ELSE
            IERR=0
            IF (FITX) CALL ANPSGL(MINX,MAXX,MINLX,MAXLX,IERR)
            IF (IERR.GT.0) RETURN
            CALL GRSCLV(REAL(MINLX,KIND(1.E0)),0.,
     .                  REAL(MAXLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
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
                     CALL GRJMP(REAL(T,KIND(1.E0)),-0.1)
                     IF (GRIDX) THEN
                        IF (RI.NE.MINLX.OR.RJ.NE.1) THEN
                          CALL GRDSH(0.2,0.5,0.2)
                          CALL GRSPTS(16)
                        ENDIF
                        CALL GRDRW(REAL(T,KIND(1.E0)),
     .                             REAL(LENY,KIND(1.E0)))
                        CALL GRDSH(1.,0.,1.)
                     ELSE
                        CALL GRDRW(REAL(T,KIND(1.E0)),0.)
                     ENDIF
   15          CONTINUE
            ELSE IF(EXPR.GE.4) THEN
               DO 20 RI=MINLX,MAXLX-1
                  CALL GRSPTS(20)
                  CALL GRJMP(REAL(RI,KIND(1.E0)),-0.1)
                  IF (GRIDX) THEN
                     IF (RI.NE.MINLX) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                     CALL GRDRW(REAL(RI,KIND(1.E0)),
     .                          REAL(LENY,KIND(1.E0)))
                     CALL GRDSH(1.,0.,1.)
                     CALL GRSPTS(20)
                  ELSE
                     CALL GRDRW(REAL(RI,KIND(1.E0)),0.)
                  ENDIF
                  DO 20 RJ=2,8,2
                     CALL GRSPTS(16)
                     ARG=RJ*10.**RI
                     T=LOG10(ARG)
                     CALL GRJMP(REAL(T,KIND(1.E0)),-0.1)
                     IF (GRIDX) THEN
                        CALL GRDSH(0.2,0.5,0.2)
                        CALL GRDRW(REAL(T,KIND(1.E0)),
     .                             REAL(LENY,KIND(1.E0)))
                        CALL GRDSH(1.,0.,1.)
                     ELSE
                        CALL GRDRW(REAL(T,KIND(1.E0)),0.)
                     ENDIF
   20          CONTINUE
            ELSE IF (EXPR.GE.1) THEN
               DO 25 RI=MINLX,MAXLX-1
                  CALL GRSPTS(20)
                  CALL GRJMP(REAL(RI,KIND(1.E0)),-0.1)
                  IF (GRIDX) THEN
                     IF (RI.NE.MINLX) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                     CALL GRDRW(REAL(RI,KIND(1.E0)),
     .                          REAL(LENY,KIND(1.E0)))
                     CALL GRDSH(1.,0.,1.)
                     CALL GRSPTS(20)
                  ELSE
                     CALL GRDRW(REAL(RI,KIND(1.E0)),0.)
                  ENDIF
                  DO 25 RJ=2,5,3
                     CALL GRSPTS(16)
                     ARG=RJ*10.**RI
                     T=LOG10(ARG)
                     CALL GRJMP(REAL(T,KIND(1.E0)),-0.1)
                     IF (GRIDX) THEN
                        CALL GRDSH(0.2,0.5,0.2)
                        CALL GRDRW(REAL(T,KIND(1.E0)),
     .                             REAL(LENY,KIND(1.E0)))
                        CALL GRDSH(1.,0.,1.)
                     ELSE
                        CALL GRDRW(REAL(T,KIND(1.E0)),0.)
                     ENDIF
   25          CONTINUE
            ELSE IF (EXPR.GE.0) THEN
               DO 30 RI=MINLX,MAXLX-1
                  CALL GRSPTS(20)
                  CALL GRJMP(REAL(RI,KIND(1.E0)),-0.1)
                  IF (GRIDX) THEN
                     IF (RI.NE.MINLX) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                     CALL GRDRW(REAL(RI,KIND(1.E0)),
     .                          REAL(LENY,KIND(1.E0)))
                     CALL GRDSH(1.,0.,1.)
                     CALL GRSPTS(20)
                  ELSE
                       CALL GRDRW(REAL(RI,KIND(1.E0)),0.)
                  ENDIF
   30          CONTINUE
            ENDIF
            CALL GRJMP(REAL(MAXLX,KIND(1.E0)),-0.1)
C
            IF (GRIDX) THEN
               CALL GRDRW(REAL(MAXLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
            ELSE
               CALL GRDRW(REAL(MAXLX,KIND(1.E0)),0.)
            ENDIF
C
            IF (.NOT.GRIDX) THEN
               CALL GRJMP(REAL(MINLX,KIND(1.E0)),0.)
               CALL GRDRW(REAL(MINLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
               CALL GRJMP(REAL(MAXLX,KIND(1.E0)),0.)
               CALL GRDRW(REAL(MAXLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
            ENDIF
C
C  10-ER LABELS AN DER X-ACHSE
C
            CALL GRCHRC(0.3,0.,20)
            DO 35 I=MINLX+1,MAXLX
               T=I-0.3*(MAXLX-MINLX)/LENX
               CALL GRTXT(REAL(T,KIND(1.E0)),-0.90,2,'10')
  35        CONTINUE
C
            DO 40 I=MINLX+1,MAXLX
               WRITE(CZHNLB,'(I3)') I
               IF (ABS(I).LT.10) THEN
                  CALL GRTXT(REAL(I,KIND(1.E0)),-0.6,2,CZHNLB(2:3))
               ELSE
                  CALL GRTXT(REAL(I,KIND(1.E0)),-0.60,3,CZHNLB)
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
                     ARG=REAL(J,KIND(1.E0))*10.**RI
                     T=LOG10(ARG)-0.05*(MAXLX-MINLX)/LENX
                     CALL GRTXT(REAL(T,KIND(1.E0)),-0.5,1,CEINLB)
   45          CONTINUE
            ELSE IF (EXPR.GE.4) THEN
               DO 50 RI=MINLX,MAXLX-1
                  DO 50 J=2,8,2
                     WRITE(CEINLB,'(I1)') J
                     ARG=REAL(J,KIND(1.E0))*10.**RI
                     T=LOG10(ARG)-0.05*(MAXLX-MINLX)/LENX
                     CALL GRTXT(REAL(T,KIND(1.E0)),-0.5,1,CEINLB)
   50          CONTINUE
            ELSE IF (EXPR.GE.2) THEN
               DO 55 RI=MINLX,MAXLX-1
                  DO 55 J=2,5,3
                     WRITE(CEINLB,'(I1)') J
                     ARG=REAL(J,KIND(1.E0))*10.**RI
                     T=LOG10(ARG)-0.05*(MAXLX-MINLX)/LENX
                     CALL GRTXT(REAL(T,KIND(1.E0)),-0.5,1,CEINLB)
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
            CALL GRSCLV(0.,REAL(MINY,KIND(1.E0)),REAL(LENX,KIND(1.E0)),
     .                  REAL(MAXY,KIND(1.E0)))
C
C   GITTERLINIEN BZW. MARKIERUNGEN AN DER Y-ACHSE
C
            DO 110 J=0,INTNRY
               T=MINY+J*STPSZY
               CALL GRJMP(-0.1,REAL(T,KIND(1.E0)))
               IF (GRIDY) THEN
                  IF (J.NE.0.AND.J.NE.INTNRY) THEN
                    CALL GRSPTS(16)
                    CALL GRDSH(0.2,0.5,0.2)
                  ENDIF
                  CALL GRDRW(REAL(LENX,KIND(1.E0)),REAL(T,KIND(1.E0)))
                  CALL GRSPTS(20)
                  CALL GRDSH(1.,0.,1.)
               ELSE
                  CALL GRDRW(0.,REAL(T,KIND(1.E0)))
               ENDIF
  110       CONTINUE
            IF (.NOT.GRIDY) THEN
               CALL GRJMP(0.,REAL(MINY,KIND(1.E0)))
               CALL GRDRW(REAL(LENX,KIND(1.E0)),REAL(MINY,KIND(1.E0)))
               CALL GRJMP(0.,REAL(MAXY,KIND(1.E0)))
               CALL GRDRW(REAL(LENX,KIND(1.E0)),REAL(MAXY,KIND(1.E0)))
            ENDIF
C
C  LABELS AN DER Y-ACHSE
C
            CALL GRCHRC(0.3,90.,18)
            DO 115 J=2,INTNRY+1,2
               FCTR=IEXP10(REAL((MINY+MAXY)*0.5,DP))
               PARAM=(MINY+(J-1)*STPSZY)/10.**FCTR
               PARAM1=PARAM*100000.
               IPARAM=NINT(PARAM1)
               PARAM=REAL(IPARAM,KIND(1.E0))/100000.
               WRITE(CPARAM,'(F7.3)') PARAM
               IL=7
               T=MINY+(J-1)*STPSZY-(MAXY-MINY)*0.8/LENY
               CALL GRTXT(-0.45,REAL(T,KIND(1.E0)),IL,CPARAM)
  115      CONTINUE
           T=MAXY-1.6*(MAXY-MINY)/LENY
           IF (ABS(FCTR).GE.10) THEN
             WRITE(CFCTR,'(I3)') FCTR
           ELSE
             WRITE(CFCTR,'(I2)') FCTR
             CFCTR(3:3)=' '
           ENDIF
           IF (FCTR.LT.0.) CALL GRTXT(-0.85,REAL(T,KIND(1.E0)),8,
     .                                '*10**'//CFCTR)
           IF (FCTR.GT.0.) CALL GRTXT(-0.85,REAL(T,KIND(1.E0)),7,
     .                                '*10**'//CFCTR(2:3))
C
C  LOGARITHMISCHE Y-ACHSENEINTEILUNG
C
         ELSE
           IERR=0
           IF (FITY) CALL ANPSGL(MINY,MAXY,MINLY,MAXLY,IERR)
           IF (IERR.GT.0) RETURN
           CALL GRSCLV(0.,REAL(MINLY,KIND(1.E0)),REAL(LENX,KIND(1.E0)),
     .                 REAL(MAXLY,KIND(1.E0)))
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
                    CALL GRJMP(-0.1,REAL(T,KIND(1.E0)))
                    IF (GRIDY) THEN
                       IF (RI.NE.MINLY.OR.RJ.NE.1) THEN
                         CALL GRDSH(0.2,0.5,0.2)
                         CALL GRSPTS(16)
                       ENDIF
                       CALL GRDRW(REAL(LENX,KIND(1.E0)),
     .                            REAL(T,KIND(1.E0)))
                       CALL GRDSH(1.,0.,1.)
                    ELSE
                       CALL GRDRW(0.,REAL(T,KIND(1.E0)))
                    ENDIF
  120         CONTINUE
           ELSE IF (EXPR.GE.4) THEN
              DO 125 RI=MINLY,MAXLY-1
                 CALL GRSPTS(20)
                 CALL GRJMP(-0.1,REAL(RI,KIND(1.E0)))
                 IF (GRIDY) THEN
                     IF (RI.NE.MINLY) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                    CALL GRDRW(REAL(LENX,KIND(1.E0)),
     .                         REAL(RI,KIND(1.E0)))
                    CALL GRSPTS(20)
                    CALL GRDSH(1.,0.,1.)
                 ELSE
                    CALL GRDRW(0.,REAL(RI,KIND(1.E0)))
                 ENDIF
                 DO 125 RJ=2,8,2
                    CALL GRSPTS(16)
                    ARG=RJ*10.**RI
                    T=LOG10(ARG)
                    CALL GRJMP(-0.1,REAL(T,KIND(1.E0)))
                    IF (GRIDY) THEN
                       CALL GRDSH(0.2,0.5,0.2)
                       CALL GRDRW(REAL(LENX,KIND(1.E0)),
     .                            REAL(T,KIND(1.E0)))
                       CALL GRDSH(1.,0.,1.)
                    ELSE
                       CALL GRDRW(0.,REAL(T,KIND(1.E0)))
                    ENDIF
  125         CONTINUE
           ELSE IF (EXPR.GE.1) THEN
              DO 130 RI=MINLY,MAXLY-1
                 CALL GRSPTS(20)
                 CALL GRJMP(-0.1,REAL(RI,KIND(1.E0)))
                 IF (GRIDY) THEN
                     IF (RI.NE.MINLY) THEN
                       CALL GRSPTS(16)
                       CALL GRDSH(0.2,0.5,0.2)
                     ENDIF
                    CALL GRDRW(REAL(LENX,KIND(1.E0)),
     .                         REAL(RI,KIND(1.E0)))
                    CALL GRSPTS(20)
                    CALL GRDSH(1.,0.,1.)
                 ELSE
                    CALL GRDRW(0.,REAL(RI,KIND(1.E0)))
                 ENDIF
                 DO 130 RJ=2,5,3
                    CALL GRSPTS(18)
                    ARG=RJ*10.**RI
                    T=LOG10(ARG)
                    CALL GRJMP(-0.1,REAL(T,KIND(1.E0)))
                    IF (GRIDY) THEN
                       CALL GRDSH(0.2,0.5,0.2)
                       CALL GRDRW(REAL(LENX,KIND(1.E0)),
     .                            REAL(T,KIND(1.E0)))
                       CALL GRDSH(1.,0.,1.)
                    ELSE
                       CALL GRDRW(0.,REAL(T,KIND(1.E0)))
                    ENDIF
  130         CONTINUE
           ELSE IF (EXPR.GE.0) THEN
              DO 135 RI=MINLY,MAXLY-1
                 CALL GRSPTS(20)
                 CALL GRJMP(-0.1,REAL(RI,KIND(1.E0)))
                 IF (GRIDY) THEN
                    IF (RI.NE.MINLY) THEN
                      CALL GRSPTS(16)
                      CALL GRDSH(0.2,0.5,0.2)
                    ENDIF
                    CALL GRDRW(REAL(LENX,KIND(1.E0)),
     .                         REAL(RI,KIND(1.E0)))
                    CALL GRSPTS(20)
                    CALL GRDSH(1.,0.,1.)
                 ELSE
                    CALL GRDRW(0.,REAL(RI,KIND(1.E0)))
                 ENDIF
  135         CONTINUE
           ENDIF
C
           CALL GRJMP(-0.1,REAL(MAXLY,KIND(1.E0)))
           IF (GRIDY) THEN
              CALL GRDRW(REAL(LENX,KIND(1.E0)),REAL(MAXLY,KIND(1.E0)))
           ELSE
              CALL GRDRW(0.,REAL(MAXLY,KIND(1.E0)))
           ENDIF
           IF (.NOT.GRIDY) THEN
              CALL GRJMP(0.,REAL(MINLY,KIND(1.E0)))
              CALL GRDRW(REAL(LENX,KIND(1.E0)),REAL(MINLY,KIND(1.E0)))
              CALL GRJMP(0.,REAL(MAXLY,KIND(1.E0)))
              CALL GRDRW(REAL(LENX,KIND(1.E0)),REAL(MAXLY,KIND(1.E0)))
           ENDIF
C
C  10-ER LABELS AN DER Y-ACHSE
C
           CALL GRCHRC(0.3,0.,20)
           DO 140 I=MINLY+1,MAXLY
              T=I-0.1*(MAXLY-MINLY)/LENY
              CALL GRTXT(-1.20,REAL(T,KIND(1.E0)),2,'10')
  140      CONTINUE
C
           DO 145 I=MINLY+1,MAXLY
              WRITE(CZHNLB,'(I3)') I
              IF (ABS(I).LT.10) THEN
                 CALL GRTXT(-0.85,REAL(I,KIND(1.E0))+0.05,2,CZHNLB(2:3))
              ELSE
                 CALL GRTXT(-0.85,REAL(I,KIND(1.E0))+0.05,3,CZHNLB)
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
                    ARG=REAL(J,KIND(1.E0))*10.**RI
                    T=LOG10(ARG)-0.1*(MAXLY-MINLY)/LENY
                    CALL GRTXT(-0.4,REAL(T,KIND(1.E0)),1,CEINLB)
  150         CONTINUE
           ELSE IF (EXPR.GE.4) THEN
              DO 155 RI=MINLY,MAXLY-1
                 DO 155 J=2,8,2
                    WRITE(CEINLB,'(I1)') J
                    ARG=REAL(J,KIND(1.E0))*10.**RI
                    T=LOG10(ARG)-0.1*(MAXLY-MINLY)/LENY
                    CALL GRTXT(-0.4,REAL(T,KIND(1.E0)),1,CEINLB)
  155         CONTINUE
           ELSE IF (EXPR.GE.2) THEN
              DO 160 RI=MINLY,MAXLY-1
                 DO 160 J=2,5,3
                    WRITE(CEINLB,'(I1)') J
                    ARG=REAL(J,KIND(1.E0))*10.**RI
                    T=LOG10(ARG)-0.1*(MAXLY-MINLY)/LENY
                    CALL GRTXT(-0.4,REAL(T,KIND(1.E0)),1,CEINLB)
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
C ===== SOURCE: plteir.f
cdr  30.4.04:  call plttly for spectra corrected.
cdr            first bin (no.0) and last bin (no. nsts+1) contain
cdr            the fluxes outside the range of spectra.
cpb  30.7.04:  deal with switched off tallies
cdr  10.6.05:  further modifications of plot for spectra (text,
c              total, plot vs. wavelength, plot 2 spectra into same frame
c    7.12.06:  in call to rstrt: one argument was wrong: sgms_cop--> sgms_bgk
!pb  18.12.06: general checking of XMCP removed to allow plots of 
!              input tallies even is no Monte Carlo particle has been followed
!    10.01.07: ENTRY PLTEIR_REINIT added for reinitialization of EIRENE
C
C
      SUBROUTINE PLTEIR (ISTRA)
C
C  ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA
C  PLOT PLASMA TALLIES ONLY ONCE.
C
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CESTIM
      USE CCONA
      USE CGRPTL
      USE CLOGAU
      USE CPLMSK
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CTRCEI
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMSOU
      USE CTEXT
      USE COUTAU
      USE CSPEI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA

      REAL(DP), ALLOCATABLE :: VECTOR(:,:),VECSAV(:,:),VSDVI(:,:)
      REAL(DP), ALLOCATABLE :: XSPEC(:),YSPEC(:,:),VSPEC(:,:),
     .          WLSPEC(:),YSPECWL(:,:),VSPECWL(:,:)
      REAL(DP) :: XXP3D_DUM(1), YYP3D_DUM(1)
      REAL(DP) :: YMN2(NPLT), YMX2(NPLT), YMNLG2(NPLT), YMXLG2(NPLT)
      REAL(DP) :: XMI, XMA, TMIN, TMAX, XI, XE, DEL, OUTAUI,
     .            SPCAN, SPC00,WL00,DE, DW
      INTEGER :: IR1(NPLT), IR2(NPLT), IRS(NPLT)
      INTEGER :: IXXE, IXXI, IYYE, IYYI, K, IINDEX, ISPC, NSPS, INULL,
     .           NF, NFT, I, IA, N, IXSET2, ISPZ, IALG, N1SDVI, ISAVE,
     .           IFIRST, IALV, ITL, JTAL, IBLD, ICURV, IE, IXSET3, IS,
     .           IERR, ICINC, IYSET3, IX, I2M, J, IRAD, I0, I1, I2, IT,
     .           INDX
      LOGICAL :: LPLOT2(NPLT), LSDVI(NPLT), LPLTT2, LINLOG, L_SAME
      CHARACTER(24) :: TXUNIT(NPLT), TXSPEC(NPLT)
      CHARACTER(24) :: TXUNT1, TXSPC1
      CHARACTER(72) :: TXTALL(NPLT)
      CHARACTER(72) :: TXTLL1
      CHARACTER(72) :: HEAD,  HEAD0, HEAD1, HEAD2, HEAD3, HEAD4,
     .                 HEAD5, HEAD6, HEAD7, HEAD8, HEAD9, HEAD10, TXHEAD
C
C
      SAVE
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        ISAVE=ISTRA
        IFIRST=1
      ENDIF
C
      IF (TRCPLT)
     .    WRITE (iunout,*) 'PLTEIR CALLED, ISTRA, XMCP: ',
     .                      ISTRA,XMCP(ISTRA)
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
C  NEW FRAME FOR EACH PICTURE IN PLTTLY
      L_SAME=.FALSE.
C
C
      IF (IESTR.EQ.ISTRA) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IF (XMCP(ISTRA).GT.1.) THEN
        IESTR=ISTRA
        CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL SYMET(ESTIMV,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IF (XMCP(ISTRA).GT.1.) THEN
        IESTR=ISTRA
        CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL SYMET(ESTIMV,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN PLTEIR: DATA FOR STRATUM ISTRA= ',
     .                    ISTRA
        WRITE (iunout,*) 'ARE NOT AVAILABLE. PLOTS ABANDONNED'
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
      HEAD9='SPECTRUM (VS. ENERGY, EV)
     .           '
      HEAD10='SPECTRUM (VS. WAVELENGTH, NM)
     .            '
C
      IALG=0
C
C  .......................................
C
C   LOOP OVER NVOLPL
C  .......................................
C
      IF (NVOLPL > 0) THEN
        ALLOCATE (VECTOR(NRAD,NPLT))
        IF (ANY(PLTL2D(1:NVOLPL).AND.PLTL3D(1:NVOLPL)))
     .     ALLOCATE (VECSAV(NRAD,NPLT))
        IF (ANY(PLTLER(1:NVOLPL))) THEN
          ALLOCATE (VSDVI(NRAD,NPLT))
          N1SDVI = NRAD
        ELSE
          ALLOCATE (VSDVI(1,1))
          N1SDVI = 1
        END IF
      END IF

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
                CALL INTTAL (ALGV,VOLTAL,IALV,NALV,NSBOX_TAL,
     .                       ALGVI(IALV,ISTRA),
     .                       NR1TAL,NP2TAL,NT3TAL,NBMLT)
105           CONTINUE
            ENDIF
            ITL=IABS(JTAL)
C  PLOT OUTPUT TALLIES ONLY FOR STRATA WITH TWO OR MORE HISTORIES
            IF (JTAL.GT.0.AND.XMCP(ISTRA).LE.1) GOTO 10000
C  PLOT INPUT TALLIES ONLY ONCE PER ITERATION
            IF (JTAL.LT.0.AND.ISAVE.NE.ISTRA) GOTO 10000
            TXHEAD=HEAD0
            IF (JTAL.GT.0)     TXHEAD=HEAD1
            IF (JTAL.EQ.NTALA) TXHEAD=HEAD2
            IF (JTAL.EQ.NTALC) TXHEAD=HEAD3
            IF (JTAL.EQ.NTALT) TXHEAD=HEAD4
            IF (JTAL.EQ.NTALM) TXHEAD=HEAD5
            IF (JTAL.EQ.NTALB) TXHEAD=HEAD6
            IF (JTAL.EQ.NTALR) TXHEAD=HEAD7
            IF (TRCPLT) THEN
              CALL LEER(1)
              WRITE (iunout,*) 'PLOT REQUESTED FOR TALLY NO. ',JTAL
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
              DO 111 I=1,NRAD
111             VECTOR(I,ICURV)=0.
              IF (ISPZ.EQ.0) THEN
                SELECT CASE (ITL)
                CASE (1)
                  VECTOR(1:NSBOX,ICURV) = TEIN(1:NSBOX)
                CASE (2)
                  VECTOR(1:NSBOX,ICURV) = SUM(TIIN(1:NF,1:NSBOX),1)
                CASE (3)
                  VECTOR(1:NSBOX,ICURV) = DEIN(1:NSBOX)
                CASE (4)
                  VECTOR(1:NSBOX,ICURV) = SUM(DIIN(1:NF,1:NSBOX),1)
                CASE (5)
                  VECTOR(1:NSBOX,ICURV) = SUM(VXIN(1:NF,1:NSBOX),1)
                CASE (6)
                  VECTOR(1:NSBOX,ICURV) = SUM(VYIN(1:NF,1:NSBOX),1)
                CASE (7)
                  VECTOR(1:NSBOX,ICURV) = SUM(VZIN(1:NF,1:NSBOX),1)
                CASE (8)
                  VECTOR(1:NSBOX,ICURV) = BXIN(1:NSBOX)
                CASE (9)
                  VECTOR(1:NSBOX,ICURV) = BYIN(1:NSBOX)
                CASE (10)
                  VECTOR(1:NSBOX,ICURV) = BZIN(1:NSBOX)
                CASE (11)
                  VECTOR(1:NSBOX,ICURV) = BFIN(1:NSBOX)
                CASE (12)
                  VECTOR(1:NSBOX,ICURV) = SUM(ADIN(1:NF,1:NSBOX),1)
                CASE (13)
                  VECTOR(1:NSBOX,ICURV) = SUM(EDRIFT(1:NF,1:NSBOX),1)
                CASE (14)
                  VECTOR(1:NSBOX,ICURV) = VOL(1:NSBOX)
                CASE (15)
                  VECTOR(1:NSBOX,ICURV) = SUM(WGHT(1:NF,1:NSBOX),1)
                CASE (16)
                  VECTOR(1:NSBOX,ICURV) = BXPERP(1:NSBOX)
                CASE (17)
                  VECTOR(1:NSBOX,ICURV) = BYPERP(1:NSBOX)
                CASE DEFAULT
                  WRITE (iunout,*) ' WRONG TALLY NUMBER IN PLTEIR',
     .                        ' JTAL = ',JTAL
                  WRITE (iunout,*) ' NO PLOT PERFORMED '
                  CALL LEER(1)
                  GOTO 10000
                END SELECT
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NF) THEN
                SELECT CASE (ITL)
                CASE (1)
                  VECTOR(1:NSBOX,ICURV) = TEIN(1:NSBOX)
                CASE (2)
                  VECTOR(1:NSBOX,ICURV) = TIIN(MPLSTI(ISPZ),1:NSBOX)
                CASE (3)
                  VECTOR(1:NSBOX,ICURV) = DEIN(1:NSBOX)
                CASE (4)
                  VECTOR(1:NSBOX,ICURV) = DIIN(ISPZ,1:NSBOX)
                CASE (5)
                  VECTOR(1:NSBOX,ICURV) = VXIN(MPLSV(ISPZ),1:NSBOX)
                CASE (6)
                  VECTOR(1:NSBOX,ICURV) = VYIN(MPLSV(ISPZ),1:NSBOX)
                CASE (7)
                  VECTOR(1:NSBOX,ICURV) = VZIN(MPLSV(ISPZ),1:NSBOX)
                CASE (8)
                  VECTOR(1:NSBOX,ICURV) = BXIN(1:NSBOX)
                CASE (9)
                  VECTOR(1:NSBOX,ICURV) = BYIN(1:NSBOX)
                CASE (10)
                  VECTOR(1:NSBOX,ICURV) = BZIN(1:NSBOX)
                CASE (11)
                  VECTOR(1:NSBOX,ICURV) = BFIN(1:NSBOX)
                CASE (12)
                  VECTOR(1:NSBOX,ICURV) = ADIN(ISPZ,1:NSBOX)
                CASE (13)
                  VECTOR(1:NSBOX,ICURV) = EDRIFT(ISPZ,1:NSBOX)
                CASE (14)
                  VECTOR(1:NSBOX,ICURV) = VOL(1:NSBOX)
                CASE (15)
                  VECTOR(1:NSBOX,ICURV) = WGHT(ISPZ,1:NSBOX)
                CASE (16)
                  VECTOR(1:NSBOX,ICURV) = BXPERP(1:NSBOX)
                CASE (17)
                  VECTOR(1:NSBOX,ICURV) = BYPERP(1:NSBOX)
                CASE DEFAULT
                  WRITE (iunout,*) ' WRONG TALLY NUMBER IN PLTEIR',
     .                        ' JTAL = ',JTAL
                  WRITE (iunout,*) ' NO PLOT PERFORMED '
                  CALL LEER(1)
                  GOTO 10000
                END SELECT
              ELSE
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE '
                  WRITE (iunout,*) 'ICURV,ISPTAL(IBLD,ICURV) ',
     .                              ICURV,ISPZ
                  WRITE (iunout,*) 
     .              'ALL PLOTS FOR THIS TALLY TURNED OFF '
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF
            ELSEIF (JTAL.GE.0) THEN
              IF (.NOT.LIVTALV(JTAL)) THEN
                WRITE (iunout,*) TXTTAL(1,JTAL)
                WRITE (iunout,*) 'TALLY SWITCHED OFF '
                WRITE (iunout,*) 'ALL PLOTS FOR THIS TALLY TURNED OFF '
                GOTO 10000
              END IF
              NFT=NFSTVI(ITL)
              NF=NFIRST(ITL)
              IF (ISPZ.EQ.0) THEN
                DO 121 I=1,NRAD
121               VECTOR(I,ICURV)=0.
                DO 122 K=1,NFT
                  DO 122 I=1,NRAD
                    VECTOR(I,ICURV)=VECTOR(I,ICURV)+
     .                              ESTIMV(NADDV(ITL)+K,NCLTAL(I))
122             CONTINUE
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NFT) THEN
                DO 125 I=1,NRAD
                  VECTOR(I,ICURV)=ESTIMV(NADDV(ITL)+ISPZ,NCLTAL(I))
125             CONTINUE
              ELSE
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE '
                  WRITE (iunout,*) 'ICURV,ISPTAL(IBLD,ICURV) ',
     .                              ICURV,ISPZ
                  WRITE (iunout,*) 
     .              'ALL PLOTS FOR THIS TALLY TURNED OFF '
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
                    VSDVI(I,ICURV)=SIGMA(N,NCLTAL(I))
127               CONTINUE
126             CONTINUE
              ENDIF
C
            ENDIF
C
            IF (PLTL2D(IBLD) .AND. PLTL3D(IBLD)) THEN
              DO 129 I=1,NRAD
                VECSAV(I,ICURV)=VECTOR(I,ICURV)
129           CONTINUE
            END IF
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
C  SET ABSCISSA FROM GRID DATA, OR USE ONE OF THE INPUT OPTIONS:
            IXSET2=0
            IF (TALXMI(IBLD).NE.0..OR.TALXMA(IBLD).NE.0.) THEN
              XMI=TALXMI(IBLD)
              XMA=TALXMA(IBLD)
              IA= 100000000
              IE=-100000000
              DO ICURV=1,NSPTAL(IBLD)
                IA=MIN(IA,NPLIN2(IBLD,ICURV))
                IE=MAX(IE,NPLOT2(IBLD,ICURV))
              ENDDO
              DEL=IE-IA
              IF (XMI.LT.XMA) THEN
C  EQUIDISTANT IN LIN SCALE
                DO I=IA,IE
                  XXP2D(I)=XMI+(I-IA)/DEL*(XMA-XMI)
                ENDDO
              ELSEIF (XMI.GT.XMA.AND.XMI.GT.0.AND.XMA.GT.0) THEN
C  EQUIDISTANT IN LOG SCALE
                XI=LOG(XMA)
                XE=LOG(XMI)
                DO I=IA,IE
                  XXP2D(I)=EXP(XI+(I-IA)/DEL*(XE-XI))
                ENDDO
C  USER DEFINED ABSCISSA, XXP2D_USR
              ELSEIF (XMI.GT.XMA.AND.(XMI.LE.0.OR.XMA.LE.0)) THEN
                DO I=IA,IE
                  XXP2D(I)=XXP2D_USR(I,IBLD)
                ENDDO
              ENDIF
              IXSET2=1
              XMI=XXP2D(IA)*(1.+1.E-6)
              XMA=XXP2D(IE)/(1.+1.E-6)
              GOTO 139
            ENDIF
C
C  TRY DEFAULT OPTION TO SET PLOT GRID FROM 1.ST (RADIAL) GRID
C
            IXSET2=0
            IF (LEVGEO.EQ.1.OR.LEVGEO.EQ.2) THEN
C   USE RADIAL SURFACE CENTERED GRID "RHOSRF" FOR EACH POLOIDAL POSITION
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
            XMI=XXP2D(NPLIN2(IBLD,1))*(1.+1.E-6)
            XMA=XXP2D(NPLOT2(IBLD,1))/(1.+1.E-6)
C
139         CONTINUE
C
            IF (IXSET2.NE.1) THEN
              WRITE (iunout,*) ' NO GRID SET FOR 2D PLOTTING '
              WRITE (iunout,*) ' IXSET2 = ',IXSET2
              WRITE (iunout,*) ' NO 2D PLOTTING IS DONE '
              GOTO 1000
            ENDIF
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
                CALL FETCH_OUTAU (OUTAUI,JTAL,ISPZ,ISTRA,IUNOUT)
                IF (OUTAUI.EQ.0.) THEN
                  IF (TRCPLT) THEN
                    WRITE (iunout,*) 'TALLY NO. ',JTAL,
     .                               ' CURVE NO. ',ICURV
                    WRITE (iunout,*) 'NOT PLOTTED BECAUSE'
                    WRITE (iunout,*) 'ZERO INTEGRAL (OUTAU(INDX)=0.) '
                    WRITE (iunout,*) 'INDX,NADDI(JTAL),NFRSTI(JTAL),I0'
                    WRITE (iunout,*)  INDX,NADDI(ITL),NFRSTI(ITL),I0
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
              I2M=I2-1
              IS=IRS(ICURV)
C
C YMNLG2, YMXLG2: REAL MAX/MIN, FOR LEGENDE ON 2D PLOT ONLY
              DO 141 I=I1,I2M,IS
                YMNLG2(ICURV)=MIN(YMNLG2(ICURV),VECTOR(I,ICURV))
141           CONTINUE
              DO 142 I=I1,I2M,IS
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
            L_SAME=.FALSE.
            CALL PLTTLY (XXP2D,VECTOR,VSDVI,YMN2,YMX2,
     .             IR1,IR2,IRS,
     .             NSPTAL(IBLD),TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .             LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,LHIST2(IBLD),IERR,
     .             N1SDVI,NRAD,L_SAME)
            IF (TRCPLT) THEN
              IF (IERR.GT.0) THEN
                WRITE (iunout,*) '2D PLOT FOR TALLY NO. ',JTAL,
     .                           ' ABANDONED'
                WRITE (iunout,*) 'ERROR CODE FROM SUBR. PLTTLY: ',IERR
                WRITE (iunout,*) 'XMI,XMA ',XMI,XMA
                GOTO 1000
              ENDIF
              WRITE (iunout,*) '2D PLOT FOR TALLY NO. ',JTAL,' DONE'
              WRITE (iunout,*) 'XMIN= ',XMI,' XMAX= ',XMA
              DO 160 ICURV=1,NSPTAL(IBLD)
                IF (LPLOT2(ICURV))
     .            WRITE (iunout,*) 'ICURV= ',ICURV,
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
              IF (PLTL2D(IBLD)) THEN
                DO 1035 I=1,NRAD
                   VECTOR(I,ICURV)=VECSAV(I,ICURV)
1035            CONTINUE
              END IF
C  SYMMETRY CONDITION AT POLAR ANGLE THETA=YIA AND THETA=2*PI+YIA
C  NOT READY: IXTL3 NOT DEFINED HERE. ENFORCE SYMMETRY AUTOMATICALLY EARLIER
C             IF (LEVGEO.EQ.2.AND.IYTL3.EQ.NP2ND) THEN
C               DO 1036 I=1,IXTL3
C1036             VECTOR(I+NP2NDM*NR1ST,ICURV)=VECTOR(I,ICURV)
C             ENDIF
1040        CONTINUE
C
C SET QUASIRECTANGULAR PLOT GRIDS XXP3D (IX), IX=1,IXTL3
C                             AND YYP3D (IY), IY=1,IYTL3
C FOR EACH 3D PICTURE
C
            IXSET3=0
            IYSET3=0
            IF (LEVGEO.EQ.1) THEN
              IF (NLRAD.AND..NOT.LPRAD3(IBLD)) THEN
                IXTL3=NR1ST
                DO 218 I=1,IXTL3
218               XXP3D(I)=RHOSRF(I)
                IXSET3=1
              ENDIF
              IF (NLTOR.AND.NLTRZ.AND..NOT.LPTOR3(IBLD)) THEN
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
              ENDIF
              IF (NLPOL.AND..NOT.LPPOL3(IBLD)) THEN
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
              IF (NLRAD.AND..NOT.LPRAD3(IBLD)) THEN
                IXTL3=NR1ST
                DO 223 I=1,IXTL3
223               XXP3D(I)=RHOSRF(I)
                IXSET3=1
              ENDIF
              IF (NLTOR.AND.NLTRZ.AND..NOT.LPTOR3(IBLD)) THEN
                IYTL3=NT3RD
                DO 226 I=1,IYTL3
226               YYP3D(I)=ZSURF(I)
                IYSET3=1
              ENDIF
              IF (NLPOL.AND..NOT.LPPOL3(IBLD)) THEN
                IYTL3=NP2ND
                DO 225 I=1,IYTL3-1
225               YYP3D(I)=0.5*(PSURF(I+1)+PSURF(I))
                YYP3D(NP2ND)=PSURF(1)+PI2A
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.3) THEN
C
              IF (LPTOR3(IBLD)) THEN
                IXTL3=NR1ST
                DO 228 IX=1,IXTL3
228               XXP3D(IX)=IX
                IXSET3=1
                IYTL3=NP2ND
                DO 230 IX=1,IYTL3
230               YYP3D(IX)=IX
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.4) THEN
C
              IF (LPTOR3(IBLD)) THEN
                IXSET3=1
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.5) THEN
C
!              IF (LPTOR3(IBLD)) THEN
                IXSET3=1
                IYSET3=1
!              ENDIF
            ELSE
              WRITE (iunout,*) '3D PLOT OPTION TO BE WRITTEN, LEVGEO '
              CALL EXIT_OWN(1)
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
                WRITE (iunout,*) ' NO GRIDS SET FOR 3D PLOTTING '
                WRITE (iunout,*) ' IXSET3,IYSET3 = ',IXSET3,IYSET3
                WRITE (iunout,*) ' NO 3D PLOTTING IS DONE '
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
     .                       IXTL3,IYTL3,XXP3D_DUM,YYP3D_DUM,
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
                IXXI=NPLI13(IBLD,ICURV)
                IXXE=NPLO13(IBLD,ICURV)
                IYYI=NPLI23(IBLD,ICURV)
                IYYE=NPLO23(IBLD,ICURV)
                CALL VECLNE (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXXI,IXXE,IYYI,IYYE,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS VECTOR PLOTS
              ELSEIF (LRPVC3(IBLD)) THEN
                CALL RPSVEC (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D_DUM,YYP3D_DUM,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
              ELSE
                WRITE (iunout,*) 'NO 3D PLOT OPTION FOR IBLD= ',IBLD
                GOTO 10000
              ENDIF
C
C   PLOT STANDARD DEVIATION PROFILE FOR THIS TALLY, IF REQUESTED
C
              IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) GOTO 1160
              IF (PLTLER(IBLD).AND.LSDVI(ICURV)) THEN
                TXUNT1='%                       '
                TXHEAD=HEAD8
                INULL=0
                DO 1222 I=1,NRAD
                  VECTOR(I,ICURV)=VSDVI(I,ICURV)
                  IF (LRAPS3(IBLD).AND.
     .               (ABS(VECTOR(I,ICURV)) < EPS30)) THEN
                    VECTOR(I,ICURV)=101._DP
                    INULL = INULL + 1
                  END IF
1222            CONTINUE
                LINLOG=.FALSE.
                TMIN=0.
                TMAX=100.
                IF (LRAPS3(IBLD).AND.INULL.GT.0) THEN
                  TMAX=101.
                  WRITE (IUNOUT,*) 'RAPS GRAPHICS FOR STD. DEVIATION:  '
                  WRITE (IUNOUT,*) 'IBLD, ICURV ',IBLD, ICURV
                  WRITE (IUNOUT,*)  INULL, ' CELLS WITH 0 HISTORIES '
                  WRITE (IUNOUT,*) 'STD. DEV. SET = 101% IN THESE CELLS'
                  WRITE (IUNOUT,*) 'TO PERMIT SPECIAL CHOICE OF COLOUR '
                END IF
                LSDVI(ICURV)=.FALSE.
                GOTO 1200
              ENDIF
C
1160        CONTINUE
C  LOOP ICURV FINISHED
          ENDIF
C
        ELSEIF (NSPTAL(IBLD).GT.0) THEN
          WRITE (iunout,*) 'PLOT REQUEST FOR TALLY NO. ',JTAL,
     .                     ' BUT NEITHER'
          WRITE (iunout,*) 
     .       'PLTL2D NOR PLTL3D TRUE. NO PLOT FOR THIS TALLY'
        ENDIF
C
10000 CONTINUE
C  LOOP IBLD FINISHED
C

      DO ISPC=1,NADSPC
C  THERE ARE NSPS BINS, AND NSPS+1 ENERGY BIN BOUNDARIES
        NSPS=ESTIML(ISPC)%PSPC%NSPC
        ALLOCATE (XSPEC(NSPS+1))
        ALLOCATE (YSPEC(NSPS+1,1))
        ALLOCATE (VSPEC(NSPS+1,1))
        SPCAN=ESTIML(ISPC)%PSPC%SPCMIN
        SPC00=ESTIML(ISPC)%PSPC%ESP_00
C  x axis: cell faces
        DO I=1,NSPS+1
          XSPEC(I)=SPCAN+(I-1)*ESTIML(ISPC)%PSPC%SPCDEL-SPC00
        END DO
C  y axis: cell averages (approx: cell centres)
        DO I=1,NSPS
          YSPEC(I,1)=ESTIML(ISPC)%PSPC%SPC(I)
          IF (NSIGI_SPC > 0) VSPEC(I,1)=ESTIML(ISPC)%PSPC%SGM(I)
        END DO

        YMN2(1)=MINVAL(YSPEC(1:NSPS,1))
        YMX2(1)=MAXVAL(YSPEC(1:NSPS,1))
        IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
        YMNLG2(1)=YMN2(1)
        YMXLG2(1)=YMX2(1)
        LSDVI(1)=NSIGI_SPC > 0
        LPLOT2(1)=.TRUE.
        IR1(1)=1
        IR2(1)=NSPS+1
        IRS(1)=1
        XMI=XSPEC(1)
        XMA=XSPEC(NSPS+1)
        LOGY=.FALSE.
        FITY=.FALSE.
        IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0) THEN
         TXTALL(1)='SPECTRUM FOR SURFACE        PARTICLE TYPE        '//
     .            'SPECIES                '
        ELSE
         TXTALL(1)='SPECTRUM FOR CELL           PARTICLE TYPE        '//
     .            'SPECIES                '
        ENDIF
        WRITE (TXTALL(1)(22:27),'(I6)') ESTIML(ISPC)%PSPC%ISPCSRF
        WRITE (TXTALL(1)(43:48),'(I6)') ESTIML(ISPC)%PSPC%IPRTYP
        WRITE (TXTALL(1)(58:63),'(I6)') ESTIML(ISPC)%PSPC%IPRSP
        IT = ESTIML(ISPC)%PSPC%ISPCTYP
        TXSPEC=REPEAT(' ',24)
        TXUNIT=REPEAT(' ',24)
        IF (IT == 1) TXUNIT='AMP/BIN(EV)             '
        IF (IT == 2) TXUNIT='WATT/BIN(EV)            '
        TXHEAD=REPEAT(' ',72)
        TXHEAD(1:30)=HEAD9(1:30)
        TXHEAD(32:42)='INTEGRAL: '
        WRITE (TXHEAD(43:55),'(ES12.4)') ESTIML(ISPC)%PSPC%SPCINT
        IERR=0
        L_SAME=.TRUE.
        IF (ISPC.EQ.1) L_SAME=.FALSE.
        L_SAME=ESTIML(ISPC)%PSPC%SPC_SAME .NE. 1.D0
        CALL PLTTLY (XSPEC,YSPEC,VSPEC,YMN2,YMX2,
     .       IR1,IR2,IRS,
     .       1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .       LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .       NSPS+1,NSPS+1,L_SAME)
        DEALLOCATE (XSPEC)
        DEALLOCATE (YSPEC)
        DEALLOCATE (VSPEC)

      END DO

C  NOW REPEAT SAME PLOTS, BUT VS. WAVELENGTH

      DO ISPC=1,NADSPC
C  THERE ARE NSPS BINS, AND NSPS+1 ENERGY BIN BOUNDARIES
        NSPS=ESTIML(ISPC)%PSPC%NSPC
        ALLOCATE (XSPEC(NSPS+1))
        ALLOCATE (YSPEC(NSPS+1,1))
        ALLOCATE (VSPEC(NSPS+1,1))
        IF (NPHOTI > 0) THEN
          ALLOCATE (WLSPEC(NSPS+1))
          ALLOCATE (YSPECWL(NSPS+1,1))
          ALLOCATE (VSPECWL(NSPS+1,1))
        END IF
        SPCAN=ESTIML(ISPC)%PSPC%SPCMIN
        SPC00=ESTIML(ISPC)%PSPC%ESP_00
C  x axis: cell faces
        DO I=1,NSPS+1
          XSPEC(I)=SPCAN+(I-1)*ESTIML(ISPC)%PSPC%SPCDEL
        END DO
C  y axis: cell averages (approx: cell centres)
        DO I=1,NSPS
          YSPEC(I,1)=ESTIML(ISPC)%PSPC%SPC(I)
          IF (NSIGI_SPC > 0) VSPEC(I,1)=ESTIML(ISPC)%PSPC%SGM(I)
        END DO

        IF (NPHOTI > 0) THEN
C  PLOT ALSO VS. WAVELENGTH (NM)
          WL00            =HPCL/MAX(1.E-6_DP,SPC00)*1.E7_DP
C  x axis: cell faces
          DO I=1,NSPS+1
            WLSPEC(NSPS-I+1+1)=HPCL/MAX(1.E-6_DP,XSPEC(I))*1.E7_DP
            WLSPEC(NSPS-I+1+1)=WLSPEC(NSPS-I+1+1)-WL00
          END DO
C  y axis: cell averages (approx: cell centres)
          DO I=1,NSPS
            YSPECWL(NSPS-I+1,1)=YSPEC(I,1)
            IF (NSIGI_SPC > 0) VSPECWL(NSPS-I+1,1)=VSPEC(I,1)
          END DO
C  rescaling:  flux/ev to flux/nm
          DO I=1,NSPS
            DE=XSPEC(NSPS-I+1+1)-XSPEC(NSPS-I+1)
            DW=WLSPEC(I+1)-WLSPEC(I)
            YSPECWL(I,1) = YSPECWL(I,1)*DE/DW
          END DO
        END IF

        DEALLOCATE (XSPEC)
        DEALLOCATE (YSPEC)
        DEALLOCATE (VSPEC)

        IF (NPHOTI > 0) THEN
          YMN2(1)=MINVAL(YSPECWL(1:NSPS,1))
          YMX2(1)=MAXVAL(YSPECWL(1:NSPS,1))
          IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
          YMNLG2(1)=YMN2(1)
          YMXLG2(1)=YMX2(1)
          LSDVI(1)=NSIGI_SPC > 0
          LPLOT2(1)=.TRUE.
          IR1(1)=1
          IR2(1)=NSPS+1
          IRS(1)=1
          XMI=WLSPEC(1)
          XMA=WLSPEC(NSPS+1)
          LOGY=.TRUE.
          FITY=.TRUE.
          IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0) THEN
            TXTALL(1)=
     .      'SPECTRUM FOR SURFACE        PARTICLE TYPE        '//
     .      'SPECIES                '
          ELSE
            TXTALL(1)=
     .      'SPECTRUM FOR CELL           PARTICLE TYPE        '//
     .      'SPECIES                '
          END IF
          WRITE (TXTALL(1)(22:27),'(I6)') ESTIML(ISPC)%PSPC%ISPCSRF
          WRITE (TXTALL(1)(43:48),'(I6)') ESTIML(ISPC)%PSPC%IPRTYP
          WRITE (TXTALL(1)(58:63),'(I6)') ESTIML(ISPC)%PSPC%IPRSP
          TXSPEC=REPEAT(' ',24)
          TXUNIT=REPEAT(' ',24)
          IF (IT == 1) TXUNIT='AMP/BIN(NM),            '
          IF (IT == 2) TXUNIT='WATT/BIN(NM),           '
          TXHEAD=REPEAT(' ',72)
          TXHEAD(1:30)=HEAD10(1:30)
          TXHEAD(32:42)='INTEGRAL: '
          WRITE (TXHEAD(43:55),'(ES12.4)') ESTIML(ISPC)%PSPC%SPCINT
          IERR=0
          L_SAME=.TRUE.
          IF (ISPC.EQ.1) L_SAME=.FALSE.
          L_SAME=ESTIML(ISPC)%PSPC%SPC_SAME .NE. 1.D0
          CALL PLTTLY (WLSPEC,YSPECWL,VSPECWL,YMN2,YMX2,
     .       IR1,IR2,IRS,
     .       1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .       LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .       NSPS+1,NSPS+1,L_SAME)
          DEALLOCATE (WLSPEC)
          DEALLOCATE (YSPECWL)
          DEALLOCATE (VSPECWL)
        END IF
      END DO


      IF (ALLOCATED(VECTOR)) DEALLOCATE(VECTOR)
      IF (ALLOCATED(VECSAV)) DEALLOCATE(VECSAV)
      IF (ALLOCATED(VSDVI))  DEALLOCATE(VSDVI)
      RETURN

C     the following ENTRY is for reinitialization of EIRENE (DMH)
      
      ENTRY PLTEIR_REINIT
      IFIRST = 0
      return
      END
C ===== SOURCE: pltin.f
C
      SUBROUTINE PLTIN (RIB,FCN,ANF,END,INN,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CLMSUR
      USE CTRCEI

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: RIB, ANF, END
      REAL(DP) :: FCN
      INTEGER, INTENT(OUT) :: INN
      LOGICAL :: LBOX
      EXTERNAL FCN

      IF (RIB.EQ.1..OR.RIB.EQ.0.)
     .                CALL PLTKI (FCN,ANF,END,INN,TRCPLT,LBOX)
      IF (RIB.EQ.1.5) CALL PLTKA (FCN,ANF,END,INN,TRCPLT,LBOX)
      IF (RIB.LT.0.)  CALL PLTKU (FCN,ANF,END,INN,TRCPLT,LBOX)
      RETURN
      END
C ===== SOURCE: pltka.f
C
C
      SUBROUTINE PLTKA (FCN,XANF,XEND,INN,TRCPLT,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CLMSUR
      USE CPLOT
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP) :: FCN
      REAL(DP), INTENT(IN) :: XANF, XEND
      INTEGER, INTENT(OUT) :: INN
      LOGICAL, INTENT(IN) ::  TRCPLT, LBOX

      REAL(DP) :: XTRAN, YTRAN, XI, ETA, GERAX, XG, YG, XINC, DX,
     .          X, Y, XB1, YB1, XP, YP, XB2, YB2, XX, YY, XXO, YYO
      INTEGER :: I, IFLAG, INC
      LOGICAL :: LBX, LBXO, LPA, L1, L2, L3, L4

      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
!pb      GERAY(XG)=(YY-YYO)/(XX-XXO)*XG+YY-XX*(YY-YYO)/(XX-XXO)
      GERAX(YG)=(YG-YY)*(XX-XXO)/(YY-YYO)+XX
      IF (TRCPLT) WRITE (iunout,*) 'PLTKA'
      IF (LBOX) THEN
        CALL GRJMP(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
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
        IF (LZR) CALL GRJMP (REAL(XX,KIND(1.E0)),REAL(YY,KIND(1.E0)))
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
            IF (LZR) CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                           REAL(YY,KIND(1.E0)))
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
                  CALL GRJMP (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
                  CALL GRDRW (REAL(XX,KIND(1.E0)),REAL(YY,KIND(1.E0)))
               ENDIF
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                  CALL STCOOR (XP,YP,0)
                  CALL STCOOR (XX,YY,1)
               ENDIF
               IFLAG=24
               INN=24
            ELSE
               IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                              REAL(YP,KIND(1.E0)))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
               IFLAG=0
            ENDIF
         ENDIF
1     CONTINUE
      RETURN
      END
C ===== SOURCE: pltki.f
C
C
      SUBROUTINE PLTKI (FCN,XANF,XEND,INN,TRCPLT,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CCONA
      USE CLMSUR
      USE CPLOT
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XANF, XEND
      REAL(DP) :: FCN
      INTEGER, INTENT(OUT) :: INN
      LOGICAL, INTENT(IN) :: TRCPLT,LBOX

      REAL(DP) :: XTRAN, YTRAN, XI, ETA, GERAX, GERAY, XG, YG, X, Y, DX,
     .          XP, YP, XX, YY, XXO, YYO, XINC
      INTEGER :: I, ISIDE, IFLAG, IJUMP, INC

      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      GERAY(XG)=(YY-YYO)/(XX-XXO)*XG+YY-XX*(YY-YYO)/(XX-XXO)
      GERAX(YG)=(YG-YY)*(XX-XXO)/(YY-YYO)+XX
      IF (TRCPLT) WRITE (iunout,*) 'PLTKI'
      IF (LBOX) THEN
        CALL GRJMP(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
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
            IF (LZR) CALL GRJMP (REAL(XX,KIND(1.E0)),
     .                           REAL(YY,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,0)
         ENDIF
         IF (IFLAG.NE.0) THEN
         IF (LZR) CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                        REAL(YY,KIND(1.E0)))
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
            IF (LZR) CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,0)
         ENDIF
         IF (IFLAG.NE.0) THEN
            IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
         ENDIF
         INN=MAX0(INN,IFLAG)
         IF (IFLAG.EQ.0) GOTO 439
         IFLAG=0
C
431   CONTINUE
      RETURN
      END
C ===== SOURCE: pltku.f
C
C
      SUBROUTINE PLTKU (FCN,XANF,XEND,INN,TRCPLT,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CLMSUR
      USE CPLOT
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP) :: FCN
      REAL(DP), INTENT(IN) :: XANF, XEND
      INTEGER, INTENT(OUT) :: INN
      LOGICAL, INTENT(IN) :: TRCPLT, LBOX

      REAL(DP) :: XTRAN, YTRAN, XI, ETA, DX, XINC, XOUT, YOUT, XS, YS,
     .          XOK, YOK, XXO, YYO, XP, YP, XTT, YTT, XT, YT, GL, GS,
     .          X, Y, XX, YY, DXX
      INTEGER :: INC, IFLAG, I, J
      LOGICAL :: LIN, LINO, LTST, LPLA, LPLAO, LPL, LGL, L1, L2, L3, L4

      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      IF (TRCPLT) WRITE (iunout,*) 'PLTKU'
      IF (LBOX) THEN
        CALL GRJMP(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
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
               IF (LZR) CALL GRJMP (REAL(XX,KIND(1.E0)),
     .                              REAL(YY,KIND(1.E0)))
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
               IF (LZR) CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                              REAL(YY,KIND(1.E0)))
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
                           CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                                 REAL(YP,KIND(1.E0)))
                           CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                                 REAL(YY,KIND(1.E0)))
                        ENDIF
                        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                           CALL STCOOR (XP,YP,0)
                           CALL STCOOR (XX,YY,1)
                        ENDIF
                        INN=24
                        IFLAG=24
                     ELSE
                        IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                                       REAL(YP,KIND(1.E0)))
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
                        CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                              REAL(YP,KIND(1.E0)))
                        CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                              REAL(YY,KIND(1.E0)))
                     ENDIF
                     IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                        CALL STCOOR (XP,YP,0)
                        CALL STCOOR (XX,YY,1)
                     ENDIF
                     INN=24
                     IFLAG=24
                  ELSE
                     IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                                    REAL(YP,KIND(1.E0)))
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
C ===== SOURCE: pltlne.f
C
C
      SUBROUTINE PLTLNE(NRET,XX,YY,XMI2D,XMA2D,YMI2D,YMA2D,LSTORE)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: XX(*), YY(*)
      REAL(DP), INTENT(IN) :: XMI2D, XMA2D, YMI2D, YMA2D
      INTEGER, INTENT(IN) :: NRET
      LOGICAL, INTENT(IN) :: LSTORE
      REAL(DP) :: XT, YT, XTN, YTN, XT2, YT2, TESTN
      INTEGER :: IN, IFL, I
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
            CALL GRDRW (REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
            IF (LSTORE) CALL STCOOR(XTN,YTN,1)
            GOTO 1485
1482      CONTINUE
            CALL GRDRW(REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
            CALL GRJMP(REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
            IF (LSTORE) THEN
              CALL STCOOR(XT,YT,0)
              CALL STCOOR(XTN,YTN,1)
            END IF
            GOTO 1485
1483      CONTINUE
            CALL GRJMP(REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
            CALL GRDRW(REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
            IF (LSTORE) THEN
              CALL STCOOR(XT,YT,0)
              CALL STCOOR(XTN,YTN,1)
            END IF
            GOTO 1485
1484      CONTINUE
            CALL GRJMP(REAL(XTN,KIND(1.E0)),REAL(YTN,KIND(1.E0)))
            IF (LSTORE) CALL STCOOR(XTN,YTN,0)
            GOTO 1485
14852     CONTINUE
            CALL GRJMP(REAL(XT,KIND(1.E0)),REAL(YT,KIND(1.E0)))
            CALL GRDRW(REAL(XT2,KIND(1.E0)),REAL(YT2,KIND(1.E0)))
            IF (LSTORE) THEN
              CALL STCOOR(XT,YT,0)
              CALL STCOOR(XT2,YT2,1)
            END IF
            GOTO 1485
1485    CONTINUE
C
148   CONTINUE
      return
      end
C ===== SOURCE: pltmsk.f
C
C
C----------------------------------------------------------------------*
C SUBROUTINE PLTMSK                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE PLTMSK(IERR)
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
C
      USE CPLMSK

      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: IERR
      REAL(DP) :: XMIN, XMAX, YMIN, YMAX
      INTEGER :: IAX
      CALL GRSCLC(REAL(X0PL,KIND(1.E0)),REAL(Y0PL,KIND(1.E0)),
     .            REAL(X0PL+LENX,KIND(1.E0)),
     .            REAL(Y0PL+LENY,KIND(1.E0)))
C  PLOT X AXIS
      IAX=1
      CALL PLTAXI(IERR,IAX)
C  PLOT Y AXIS
      IAX=2
      CALL PLTAXI(IERR,IAX)
      IF (IERR.GT.0) RETURN
C
      IF (LOGX) THEN
         XMIN=REAL(MINLX,KIND(1.E0))
         XMAX=REAL(MAXLX,KIND(1.E0))
      ELSE
         XMIN=MINX
         XMAX=MAXX
      ENDIF
C
      IF (LOGY) THEN
         YMIN=REAL(MINLY,KIND(1.E0))
         YMAX=REAL(MAXLY,KIND(1.E0))
      ELSE
         YMIN=MINY
         YMAX=MAXY
      ENDIF
C
      CALL GRSCLV(REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)),
     .            REAL(XMAX,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
C
      RETURN
      END
C ===== SOURCE: plttly.f
C  10.6.05:  L_SAME:  USE SAME FRAME AS IN PREVIOUS CALL
C  8.8.06 :  GRPP taken out
C
      SUBROUTINE PLTTLY (X,Y,VBAR,YMN,YMX,IR1,IR2,IRS,NKURV,TXTTAL,
     .                   TXTSPC,TXTUNT,TXTRUN,TXHEAD,
     .                   LBAR,XMI,XMA,YMNLG,YMXLG,LPLOT,LHIST,IERR,
     .                   N1BAR,N1DIM,L_SAME)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CPLMSK

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N1BAR, N1DIM
      REAL(DP), INTENT(IN) :: X(*), VBAR(N1BAR,*),
     .                      YMN(*), YMX(*), YMNLG(*), YMXLG(*)
      REAL(DP), INTENT(INOUT) :: Y(N1DIM,*)
      INTEGER, INTENT(IN) :: IR1(*), IR2(*), IRS(*), NKURV
      LOGICAL, INTENT(IN) :: LBAR(*), LPLOT(*), L_SAME
      LOGICAL, INTENT(IN) :: LHIST
      CHARACTER(LEN=*), INTENT(IN) :: TXTTAL(*),TXTSPC(*),TXTUNT(*),
     .                                TXTRUN, TXHEAD

      REAL(DP) :: YA, YMINY, YMY, AA, FM, ST1, ST2, FP, DMINY, DMAXY,
     .          XMI, XMA
      REAL(DP) :: XMIN, XMAX, YMIN, YMAX
      REAL(SP), SAVE :: PRMSAVE(8)
      INTEGER :: I1, I2, IS, ICURV, IT, ISY, NP, J, NPS, IKURV,
     .           I, IERR, IPEN1, IPEN2
      CHARACTER(10) :: CHR
      CHARACTER(12) :: CHR12
      SAVE YA,IPEN1,IPEN2
C
      IKURV=0
      DO 1 I=1,NKURV
        IF (LPLOT(I)) IKURV=IKURV+1
1     CONTINUE
      IF (IKURV.EQ.0) RETURN
C
      MINX=XMI
      MAXX=XMA
      MINY=1.D30
      MAXY=-1.D30
      DO 3 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 3
        MINY=MIN(MINY,REAL(YMN(I),SP))
        MAXY=MAX(MAXY,REAL(YMX(I),SP))
3     CONTINUE
C
      IF (LOGY) THEN
        DMAXY=MAXY
        DMINY=MINY
        DMAXY=MAX(1.E-36_DP,DMAXY)
        MAXY=DMAXY
        DMINY=MAX(DMINY,DMAXY/1.E12_DP)
        MINY=DMINY
      ENDIF
C
      IF (MINY.GE.MAXY) THEN
        IERR=1
        RETURN
      ENDIF
C
      IF (.NOT.L_SAME) THEN
C  NEW FRAME
        IPEN1=0
        IPEN2=0
        CALL GRNXTB(1,'PLTTLY.F')
        CALL GRSCLC (0.,0.,39.5,28.7)
        CALL GRSCLV (0.,0.,39.5,28.7)
        IT=LEN(TXTRUN)
        CALL GRTXT (1.,27.5,IT,TXTRUN)
        IT=LEN(TXHEAD)
        CALL GRTXT (1.,26.75,IT,TXHEAD)
        CALL GRTXT (1.,26.,15,'MIN. ABSZISSA =')
        WRITE (CHR12,'(1P,E12.4)') XMI
        CALL GRTXTC (12,CHR12)
        CALL GRTXT (1.,25.5,15,'MAX. ABSZISSA =')
        WRITE (CHR12,'(1P,E12.4)') XMA
        CALL GRTXTC (12,CHR12)
        YA=24.75
      ELSE         !   IF (L_SAME) THEN
! SET SCALING FOR TEXT
        CALL GRSCLC (0.,0.,39.5,28.7)
        CALL GRSCLV (0.,0.,39.5,28.7)
      END IF

      DO 2 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 2
        IPEN1=IPEN1+1
        CALL GRNWPN(IPEN1)
        ISY=IPEN1+1
        CALL GRJMPS (0.5,REAL(YA,KIND(1.E0)),ISY)
        IT=LEN(TXTTAL(I))
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),IT,TXTTAL(I))
        YA=YA-0.5
        IT=LEN(TXTSPC(I))
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),IT,TXTSPC(I))
        YA=YA-0.5
        IT=LEN(TXTUNT(I))
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),IT,TXTUNT(I))
        YA=YA-0.5
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),11,'MAX. VALUE =')
        WRITE (CHR,'(1P,E10.3)') YMXLG(I)
        CALL GRTXTC (10,CHR)
        YA=YA-0.5
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),11,'MIN. VALUE =')
        WRITE (CHR,'(1P,E10.3)') YMNLG(I)
        CALL GRTXTC (10,CHR)
        YA=YA-1.0
2     CONTINUE
      CALL GRNWPN(1)

      IF (L_SAME) THEN
! RESTORE SCALING
        CALL GRSCLC (PRMSAVE(1),PRMSAVE(2),PRMSAVE(3),PRMSAVE(4))
        CALL GRSCLV (PRMSAVE(5),PRMSAVE(6),PRMSAVE(7),PRMSAVE(8))
      END IF
C
C
C  PLOT X AND Y AXIS
C
      IF (.NOT.L_SAME) THEN
        CALL PLTMSK(IERR)
        IF (LOGX) THEN
          XMIN=REAL(MINLX,KIND(1.E0))
          XMAX=REAL(MAXLX,KIND(1.E0))
        ELSE
          XMIN=MINX
          XMAX=MAXX
        ENDIF
C
        IF (LOGY) THEN
          YMIN=REAL(MINLY,KIND(1.E0))
          YMAX=REAL(MAXLY,KIND(1.E0))
        ELSE
          YMIN=MINY
          YMAX=MAXY
        ENDIF
        IF (ABS((YMAX-YMIN)/(YMAX+EPS60)) < EPS6) THEN
          YMIN = YMIN - 1.
          YMAX = YMAX + 1.
        END IF
        PRMSAVE(1) = X0PL
        PRMSAVE(2) = Y0PL
        PRMSAVE(3) = X0PL+LENX
        PRMSAVE(4) = Y0PL+LENY
        PRMSAVE(5) = XMIN
        PRMSAVE(6) = YMIN
        PRMSAVE(7) = XMAX
        PRMSAVE(8) = YMAX
        IF (IERR.GT.0) RETURN
      ENDIF
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
        IPEN2=IPEN2+1
        CALL GRNWPN(IPEN2)
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
          CALL GRJMP(REAL(X(I1),KIND(1.E0)),REAL(YMY,KIND(1.E0)))
          CALL GRDRW(REAL(X(I1),KIND(1.E0)),REAL(Y(I1,I),KIND(1.E0)))
          CALL GRDRW(REAL(X(I1+IS),KIND(1.E0)),REAL(Y(I1,I),KIND(1.E0)))
          DO J=I1+IS,I2,IS
            CALL GRDRW (REAL(X(J),KIND(1.E0)),REAL(Y(J,I),KIND(1.E0)))
            CALL GRDRW (REAL(X(J+IS),KIND(1.E0)),
     .                  REAL(Y(J,I),KIND(1.E0)))
          END DO
          CALL GRDRW(REAL(X(I2+IS),KIND(1.E0)),REAL(YMY,KIND(1.E0)))
C  PLOT SYMBOLS
          ISY=IPEN2+1
          NP=(I2-I1+1)/IS
          NPS=MAX0(NP/7,1)
          DO J=I1,I2,IS*NPS
            CALL GRJMPS (REAL(0.5*(X(J)+X(J+IS)),KIND(1.E0)),
     .                   REAL(Y(J,I),KIND(1.E0)),ISY)
          END DO
C
        ELSEIF (.NOT.LHIST) THEN
C  PLOT LINES
          CALL GRJMP(REAL(X(I1),KIND(1.E0)),REAL(Y(I1,I),KIND(1.E0)))
          DO 33 J=I1+IS,I2,IS
33          CALL GRDRW (REAL(X(J),KIND(1.E0)),REAL(Y(J,I),KIND(1.E0)))
C  PLOT SYMBOLS
          ISY=IPEN2+1
          NP=(I2-I1+1)/IS
          NPS=MAX0(NP/7,1)
          DO J=I1,I2,IS*NPS
            CALL GRJMPS (REAL(X(J),KIND(1.E0)),REAL(Y(J,I),KIND(1.E0)),
     .                   ISY)
          END DO
        ENDIF
C
C  PLOT ERROR BARS
          IF (LBAR(I)) THEN
            DO 40 J=I1,I2,IS
              FP=1.+VBAR(J,I)*0.01
              FM=1.-VBAR(J,I)*0.01
              IF (LOGY) THEN
                AA=10.**Y(J,I)
                ST1=LOG10(MAX(1.E-30_DP,AA*FM))
                ST2=LOG10(MAX(1.E-30_DP,AA*FP))
                ST1=MAX(REAL(MINLY,DP),MIN(REAL(MAXLY,DP),ST1))
                ST2=MAX(REAL(MINLY,DP),MIN(REAL(MAXLY,DP),ST2))
              ELSE
                ST1=Y(J,I)*FM
                ST2=Y(J,I)*FP
                ST1=MAX(REAL(MINY,DP),MIN(REAL(MAXY,DP),ST1))
                ST2=MAX(REAL(MINY,DP),MIN(REAL(MAXY,DP),ST2))
              ENDIF
            IF (LHIST) THEN
              CALL GRJMP (REAL(0.5*(X(J)+X(J+IS)),KIND(1.E0)),
     .                    REAL(ST1,KIND(1.E0)))
              CALL GRDRW (REAL(0.5*(X(J)+X(J+IS)),KIND(1.E0)),
     .                    REAL(ST2,KIND(1.E0)))
            ELSE
              CALL GRJMP (REAL(X(J),KIND(1.E0)),REAL(ST1,KIND(1.E0)))
              CALL GRDRW (REAL(X(J),KIND(1.E0)),REAL(ST2,KIND(1.E0)))
            ENDIF
40        CONTINUE
        ENDIF
50    CONTINUE
      CALL GRNWPN(1)

      CALL GRCHRC (0.3,0.,16)
C
      RETURN
      END
C ===== SOURCE: prllo.f
C
C
      SUBROUTINE PRLLO (P1,P2,P3,P4,P5,IO,NF)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3), P5(3)
      INTEGER, INTENT(IN) :: IO
      LOGICAL, INTENT(IN) :: NF
      REAL(DP) :: CORD(15)
      INTEGER :: I

      DO 1 I=1,3
      CORD(I)=P1(I)
      CORD(3+I)=P2(I)
      CORD(6+I)=P4(I)
      CORD(9+I)=P5(I)
1     CORD(12+I)=P3(I)
      CALL PL3Q (CORD,5,IO,NF)
      RETURN
      END
C ===== SOURCE: pside.f
C
C
      SUBROUTINE PSIDE (G,R,P11,P12,P21,P22,L1,L2,EPS)

      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CPLOT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: G(4,2), R(4,2)
      REAL(DP), INTENT(IN) :: P11, P12, P21, P22, EPS
      LOGICAL, INTENT(IN) :: L1, L2
      REAL(DP) :: XLA, XP, YP, XMU
      INTEGER :: I, IS

      IS=0
      DO 100 I=1,4
      CALL GSP (G(I,1),G(I,2),R(I,1),R(I,2),P11,P12,P21-P11,P22-P12,
     .          XLA,XMU,EPS)
      IF (XLA.GE.0..AND.XLA.LE.1..AND.XMU.GE.0..AND.XMU.LE.1.) THEN
        XP=G(I,1)+XMU*R(I,1)
        YP=G(I,2)+XMU*R(I,2)
        IF (L1) THEN
          IF (LZR) THEN
            CALL GRJMP (REAL(P11,KIND(1.E0)),REAL(P12,KIND(1.E0)))
            CALL GRDRW (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
          ENDIF
          IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
            CALL STCOOR (P11,P12,0)
            CALL STCOOR (XP,YP,1)
          ENDIF
          RETURN
        ELSE IF (L2) THEN
          IF (LZR) THEN
            CALL GRJMP (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
            CALL GRDRW (REAL(P21,KIND(1.E0)),REAL(P22,KIND(1.E0)))
          ENDIF
          IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
            CALL STCOOR (XP,YP,0)
            CALL STCOOR (P21,P22,1)
          ENDIF
          RETURN
        ELSE
          IS=IS+1
          IF (IS.EQ.1) THEN
            IF (LZR) CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,0)
          ELSEIF (IS.EQ.2) THEN
            IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
            RETURN
          ENDIF
        ENDIF
      ENDIF
100   CONTINUE
      RETURN
      END
C ===== SOURCE: rpscol.f
!pb 19.12.06: write minimum and maximum value onto RAPS log-file
C
      SUBROUTINE RPSCOL (AORIG,IBLD,ICURV,
     .                   I1,I2,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A PLOTFILE FOR THE RAPS GRAPHICS SYSTEM
C
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ONTO VERTICES IS PERFORMED
C
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
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTRIG
      USE COMPRT, ONLY: IUNOUT
      USE CCONA

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV, I1, I2
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP), ALLOCATABLE :: YWERT(:,:),ywert1(:,:), BORIG(:)
      REAL(DP) :: WMIN, WMAX
      INTEGER :: IR, IERR, IT, I, IPART, IP, ICASE, IRD
C
      WRITE (60,*) RUNID
      WRITE (60,*) TXHEAD
      WRITE (60,*) HEAD
      WRITE (60,*) TEXT1
      WRITE (60,*) TEXT2
      WRITE (60,*) TEXT3
C
      NRAPS=NRAPS+1
      IRAPS=IRAPS+1
C
      OPEN (UNIT=NRAPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS
C
      ICASE=0
      IERR=0
C  2D (X,Y) FEM MESH
      IF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
         ALLOCATE (YWERT1(NRAD,1))
         CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         ICASE=4
         WMIN=MINVAL(YWERT1(1:NRAD,1))
         WMAX=MAXVAL(YWERT1(1:NRAD,1))
C  2D MESH, EITHER X,Y  OR X,Z OR Y,Z PROJECTION
      ELSEIF (LEVGEO.LE.3.AND.
     .       (LPPOL3(IBLD).OR.LPTOR3(IBLD).OR.LPRAD3(IBLD))) THEN
         ALLOCATE (YWERT(N1ST,N2ND+N3RD))
         CALL CELINT(AORIG,YWERT,LOGL,IBLD,ICURV,N1ST,IERR)
         ICASE=3
         WMIN=MINVAL(YWERT(1:N1ST,1:N2ND+N3RD))
         WMAX=MAXVAL(YWERT(1:N1ST,1:N2ND+N3RD))
C  3D MESH, X,Y,Z CUBE,  NO PROJECTION
      ELSEIF (LEVGEO.EQ.1.AND.NLTRZ
     .   .AND.NLRAD.AND.NLPOL.AND.NLTOR
     .   .AND..NOT.(LPPOL3(IBLD).OR.LPTOR3(IBLD).OR.LPRAD3(IBLD))) THEN
         ALLOCATE (YWERT1(NRAD,1))
         CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         ICASE=1
         WMIN=MINVAL(YWERT1(1:NRAD,1))
         WMAX=MAXVAL(YWERT1(1:NRAD,1))
C  3D MESH, TETRAHEDONS, NO PROJECTION
      ELSEIF ((LEVGEO.EQ.5).AND..NOT.LRPSCUT) THEN
         ALLOCATE (YWERT1(NRAD,1))
         CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         ICASE=5
         WMIN=MINVAL(YWERT1(1:NRAD,1))
         WMAX=MAXVAL(YWERT1(1:NRAD,1))
C  3D MESH, TETRAHEDONS, PROJECTION INTO PLANE
      ELSEIF ((LEVGEO.EQ.5).AND.LRPSCUT) THEN
         ALLOCATE (YWERT1(NRAD,1))
         ALLOCATE (BORIG(NRAD))
         CALL RPSCUT (AORIG,BORIG)
         CALL CELINT(BORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         DEALLOCATE (BORIG)
         ICASE=4
         WMIN=MINVAL(YWERT1(1:NRAD,1))
         WMAX=MAXVAL(YWERT1(1:NRAD,1))
      ENDIF
      IF (ZMI .NE. 666.) WMIN = MIN(ZMI, WMIN)
      IF (ZMA .NE. 666.) WMAX = MAX(ZMA, WMAX)
      WRITE (60,*) WMIN, WMAX
      WRITE (60,*)

      IF (IERR.GT.0) THEN
        IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
        IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
        RETURN
      END IF
C
      IF (ICASE.EQ.1) THEN

         do ir=1,nr1st
            do ip=1,np2nd
               do it=1,nt3rd
                  IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IF (ZMI.NE.666.) YWERT1(IRD,1)=
     .                             MAX(YWERT1(IRD,1),ZMI)
                  IF (ZMA.NE.666.) YWERT1(IRD,1)=
     .                             MIN(YWERT1(IRD,1),ZMA)
                  IF (ABS(YWERT1(IRD,1)) < EPS30) YWERT1(IRD,1)=0._DP
                  WRITE (NRAPS,*) YWERT1(IRD,1)
               enddo
            enddo
         enddo


      ELSEIF (LEVGEO.LE.2.AND.LPPOL3(IBLD)) THEN
        LPPOLR=.TRUE.
        IPPOLR=MAX(1,IPROJ3(IBLD,ICURV))
        IF (.NOT.NLRAD.OR..NOT.NLTOR) THEN
          WRITE (iunout,*) 'ERROR IN RPSCOL. LPPOL3? '
          IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
          IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
          RETURN
        ENDIF
        DO 1100 IR=1,NR1ST
          DO 3100 IT=1,NT3RD
            IF (ZMI.NE.666.) YWERT(IR,IT)=MAX(YWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IT)=MIN(YWERT(IR,IT),ZMA)
            IF (ABS(YWERT(IR,IT)) < EPS30) YWERT(IR,IT)=0._DP
            WRITE (NRAPS,*) YWERT(IR,IT)
3100      CONTINUE
1100    CONTINUE
C
      ELSEIF (LEVGEO.LE.2.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        IPTORR=MAX(1,IPROJ3(IBLD,ICURV))
        IF (.NOT.NLRAD.OR..NOT.NLPOL) THEN
          WRITE (iunout,*) 'ERROR IN RPSCOL. LPTOR3? '
          IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
          IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
          RETURN
        ENDIF
        DO 1 IR=1,NR1ST
          DO 3 IP=1,NP2ND
            IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
            IF (ABS(YWERT(IR,IP)) < EPS30) YWERT(IR,IP)=0._DP
            WRITE (NRAPS,*) YWERT(IR,IP)
3         CONTINUE
1       CONTINUE
C
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        IPTORR=MAX(1,IPROJ3(IBLD,ICURV))
        IF (.NOT.NLRAD.OR..NOT.NLPOL) THEN
          WRITE (iunout,*) 'ERROR IN RPSCOL. LPTOR3? '
          IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
          IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
          RETURN
        ENDIF
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
              IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
              IF (ABS(YWERT(IR,IP)) < EPS30) YWERT(IR,IP)=0._DP
              WRITE (NRAPS,*) YWERT(IR,IP)
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
C  icase=4
!PB      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
      ELSEIF (ICASE == 4) THEN
        LPTORR=.TRUE.
        IPTORR=MAX(1,IPROJ3(IBLD,ICURV))
        DO 60 I=1,NRKNOT
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          IF (ABS(YWERT1(I,1)) < EPS30) YWERT1(I,1)=0._DP
          WRITE(NRAPS,*) YWERT1(I,1)
60      CONTINUE
C
C  icase=5
      ELSEIF ((LEVGEO.EQ.5).AND..NOT.LRPSCUT) THEN
        DO I=1,NCOORD
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          IF (ABS(YWERT1(I,1)) < EPS30) YWERT1(I,1)=0._DP
          WRITE(NRAPS,*) YWERT1(I,1)
        enddo
C
      ELSE
        WRITE (iunout,*) 'UNWRITTEN OPTION IN RPSCOL: PLOT ABANDONNED '
      ENDIF
C
      CLOSE (UNIT=NRAPS)
C
      IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
      IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)

      RETURN
      END



C ===== SOURCE: rpscut.f
      subroutine rpscut (aorig,values)

      use precision
      use parmmod
      use ctetra
      use ctrig
      use cgeom
      use ccona
      use cplot
      use module_avltree

      implicit none

      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(OUT) :: VALUES(*)

      integer, allocatable, save :: tetra_kanten(:,:), kanten(:,:), 
     .             iedge(:), icono(:), itetno(:)
      integer, save :: kanten_nummer(4,4) = reshape(
     .    (/ -1,  1,  2,  3,
     .        1, -1,  4,  5,
     .        2,  4, -1,  6,
     .        3,  5,  6, -1 /), (/ 4, 4 /) )
      integer, save :: angrenzende_seiten(2,6) = reshape (
     .     (/ 1, 2,
     .        1, 4,
     .        2, 4,
     .        1, 3,
     .        2, 3,
     .        3, 4 /), (/ 2, 6 /) )
      integer, save :: nkanten, ifirst=0
      real(dp) :: timi, timen, second_own
      real(dp), allocatable, save :: spar(:)
      integer :: i, itri
 
      
!pb   if (.not.allocated(tetra_kanten)) then
      if (ifirst == 0) then
        ifirst = 1
        allocate (tetra_kanten(6,ntet))
        allocate (kanten(2,3*ntet))
        tetra_kanten=0
        kanten=0
        nkanten=0
!        call suche_kanten
!      timen = second_own()
!      write (0,*) ' cpu time spend in suche_kanten ',timen-timi,' sec'
!        kanten=0
!        nkanten=0
        timi = second_own()
        call suche_kanten2
        timen = second_own()
        write (0,*) ' cpu time spend in suche_kanten2 ',
     .                timen-timi,' sec'

        call schneide_kanten 
        call berechne_koordinaten
        call bilde_dreiecke 

        deallocate (tetra_kanten)
        deallocate (kanten)
        deallocate (iedge)
        deallocate (icono)
        deallocate (spar)
      end if

      do i=1,ntrii
        values(i) = aorig(itetno(i))
      end do

      return


      contains


      subroutine suche_kanten

      implicit none
      type(tet_elem), pointer :: cur
      integer :: itet, j, nxt_side, akt_tet, akt_edge, isi, ip1
      integer :: ic, i, nump1, nump2, noedge, jc, no_side, nxt_tet

      nkanten = 0
! für jede Koordinate
      do ic=1,ncoord

! durchsuche die Kanten der angrenzenden Tetraeder
        cur => coortet(ic)%ptet
        do while (associated(cur))
          itet = cur%notet

! Koordinate ic ist Eckpunkt des Tetraeders itet, 
! bestimme die Nummer der Ecke
          do i=1,4
            if (nteck(i,itet) == ic) ip1 = i
        end do  

! bearbeite alle Kanten, die von Eckpunkt nump1 ausgehen
! d.h. die Kanten zu den anderen Eckpunkten 
        do i=1,4
          nump1 = ip1
          if (i == nump1) cycle
          nump2 = i
          noedge = kanten_nummer(nump1,nump2)
! die Kante ist schon bei einer anderen Koordinate gefunden worden 
            if (tetra_kanten(noedge,itet) /= 0) cycle
! die Kante ist neu
          jc = nteck(nump2,itet)
            nkanten = nkanten+1
          kanten(1,nkanten) = ic
          kanten(2,nkanten) = jc
          tetra_kanten(noedge,itet) = nkanten

! durchlaufe alle Nachbartetraeder und markiere die gemeinsame Kante
        
          akt_tet = itet
          akt_edge = noedge
            no_side = angrenzende_seiten(1,akt_edge)
          isi = 0
          do
            nxt_tet = ntbar(no_side,akt_tet)
          if (nxt_tet == 0) then
            isi = isi+1
! an beiden Seiten bis zum Rand gelaufen
            if (isi > 1) exit
                no_side = angrenzende_seiten(2,akt_edge)
            nxt_tet = ntbar(no_side,itet)
            if (nxt_tet == 0) exit
            akt_tet = itet
          end if
            nxt_side = ntseite(no_side,akt_tet)
          if (nxt_tet == itet) exit
! bestimme die Kantennummer auf dem neuen Tetraeder
          do j=1,4
            if (nteck(j,nxt_tet) == ic) nump1 = j
            if (nteck(j,nxt_tet) == jc) nump2 = j
          end do
          noedge = kanten_nummer(nump1,nump2)
            tetra_kanten(noedge,nxt_tet) = nkanten
          no_side = angrenzende_seiten(1,noedge) + 
     .              angrenzende_seiten(2,noedge) - nxt_side
          akt_tet = nxt_tet
          end do  
                
        end do  

          cur => cur%next_tet

        end do  
      
      end do

      write (0,*) 'min(tetra_kanten) ',minval(tetra_kanten(:,1:ntet))
      return
      end subroutine suche_kanten
      

      subroutine suche_kanten2

      implicit none
      integer :: itet, j, nxt_side, akt_tet, akt_edge, isi, ik
      integer :: ic, i, nump1, nump2, noedge, jc, no_side, nxt_tet

      integer :: kpunkte(2,6) = reshape (
     .           (/ 1,2,  1,3,  1,4,  2,3,  2,4,  3,4 /), (/ 2, 6 /) )

      nkanten = 0
! für jeden Tetraeder
      do itet=1,ntet

! durchsuche die Kanten des Tetraeders
        do ik=1,6
          
! die Kante ist schon bei einem anderen Tetraeder gefunden worden 
          if (tetra_kanten(ik,itet) /= 0) cycle

! die Kante ist neu
          nump1 = kpunkte(1,ik)
          nump2 = kpunkte(2,ik)
          ic = nteck(nump1,itet)
          jc = nteck(nump2,itet)

          nkanten = nkanten+1
          kanten(1,nkanten) = ic
        kanten(2,nkanten) = jc
        tetra_kanten(ik,itet) = nkanten
      
! durchlaufe alle Nachbartetraeder und markiere die gemeinsame Kante
        
          noedge = ik
        akt_tet = itet
        akt_edge = noedge
          no_side = angrenzende_seiten(1,akt_edge)
        isi = 0
        do
          nxt_tet = ntbar(no_side,akt_tet)
          if (nxt_tet == 0) then
            isi = isi+1
! an beiden Seiten bis zum Rand gelaufen
          if (isi > 1) exit
              no_side = angrenzende_seiten(2,akt_edge)
            nxt_tet = ntbar(no_side,itet)
            if (nxt_tet == 0) exit
            akt_tet = itet
          end if
          nxt_side = ntseite(no_side,akt_tet)
          if (nxt_tet == itet) exit
! bestimme die Kantennummer auf dem neuen Tetraeder
          do j=1,4
            if (nteck(j,nxt_tet) == ic) nump1 = j
            if (nteck(j,nxt_tet) == jc) nump2 = j
          end do
          noedge = kanten_nummer(nump1,nump2)
          tetra_kanten(noedge,nxt_tet) = nkanten
          no_side = angrenzende_seiten(1,noedge) + 
     .              angrenzende_seiten(2,noedge) - nxt_side
          akt_tet = nxt_tet
        end do  

        end do  ! ik
      
      end do  ! itet
      return
      end subroutine suche_kanten2


      subroutine schneide_kanten
      implicit none
      integer :: ik, k1, k2, ic
      real(dp) :: p1(3), p2(3), v(3), a(3)
      real(dp) :: d, xnen, t, x, y, z, dst
      logical :: inserted

      type(TAVLTree), pointer :: baum

      allocate(spar(nkanten))
      allocate(iedge(nkanten))
      allocate(icono(nkanten))

      spar = 1.E30_dp
      iedge = 0
      icono = 0
      nrknot = 0
      baum => NewTree()

      a(1:3) = cutplane(2:4)
      d = cutplane(1)

      do ik=1,nkanten

        k1 = kanten(1,ik)
        k2 = kanten(2,ik)
      
        p1(1) = xtetra(k1)
        p1(2) = ytetra(k1)
        p1(3) = ztetra(k1)

        p2(1) = xtetra(k2)
        p2(2) = ytetra(k2)
        p2(3) = ztetra(k2)

        v = p2 - p1 

        xnen = sum(a*v)
        if (abs(xnen) < eps10) cycle

        t = -(d + sum(a*p1)) / xnen
        
        if((t > -eps6) .and. (t < 1._dp+eps6)) then  

          x = xtetra(k1) + t*(xtetra(k2)-xtetra(k1)) 
          y = ytetra(k1) + t*(ytetra(k2)-ytetra(k1)) 
          z = ztetra(k1) + t*(ztetra(k2)-ztetra(k1)) 
          dst = sqrt( (xtetra(k2)-xtetra(k1))**2 +
     .                (ytetra(k2)-ytetra(k1))**2 +
     .                (ztetra(k2)-ztetra(k1))**2 )

          ic = nrknot + 1
          inserted=.false.
          call insert (baum, x, y, z, dst, ic, inserted)

          if (inserted) then
            nrknot = nrknot + 1
            iedge(nrknot) = ik
          end if

          spar(ik) = t
          icono(ik) = ic
        end if  

      end do

      call DestroyTree(baum)

      end subroutine schneide_kanten


      subroutine berechne_koordinaten

      implicit none
      real(dp) :: a(3), b(3), c(3), p(3), orig(3)
      real(dp) :: am(2,2), amm1(2,2), rhs(2), x(3)
      real(dp) :: spat, bnorm, cnorm, detam
      integer :: i, ik, k1, k2

! a ist die richtung der normalen zur schnittebene
      a(1:3) = cutplane(2:4)

! berechne vektor b, b senkrecht zu a
! d.h.  a1*b1 + a2*b2 + a3*b3 = 0

      if (abs(a(1)) > eps10) then
        b(1) = -(a(2)+a(3))/a(1)
        b(2) = 1._dp
        b(3) = 1._dp
      else if (abs(a(2)) > eps10) then
        b(1) = 1._dp
        b(2) = -(a(1)+a(3))/a(2)
        b(3) = 1._dp
      else if (abs(a(3)) > eps10) then
        b(1) = 1._dp
        b(2) = 1._dp
        b(3) = -(a(1)+a(2))/a(3)
      else 
        write (6,*) ' Problem in berechne_koordinaten '
        write (6,*) ' die Koeffizienten der Schnittebene sind 0 '
        write (6,*) cutplane
        call exit (1) 
      end if 

! berechne vektor c = a x b ==> c ist senkrecht zu a und b

      c(1) = a(2)*b(3) - b(2)*a(3)
      c(2) = a(3)*b(1) - b(3)*a(1)
      c(3) = a(1)*b(2) - b(1)*a(2)
      
! berechne spatprodukt, (a x b)c > 0 ==> rechtsystem, sonst c=-c

      spat =  a(1)*b(2)*b(3) + a(2)*b(3)*c(1) + a(3)*b(1)*c(2)
     .      - c(1)*b(2)*a(3) - c(2)*b(3)*a(1) - c(3)*b(1)*a(2)
      if (spat < 0._dp) c = -c
      
! normiere b und c
      bnorm = sqrt(sum(b*b))
      cnorm = sqrt(sum(c*c))  
      
      b = b / bnorm
      c = c / cnorm
      
! b und c sind eine orthonormalbasis der schnittebene 

              
      nknot = nrknot
      nknots = nknot
      ntri = ntet / 2
      ntris = ntri

      call dealloc_ctrig
      call alloc_ctrig

      ik = iedge(1)
      k1 = kanten(1,ik)
      k2 = kanten(2,ik)
      orig(1) = xtetra(k1) + spar(ik)*(xtetra(k2)-xtetra(k1))
      orig(2) = ytetra(k1) + spar(ik)*(ytetra(k2)-ytetra(k1))
      orig(3) = ztetra(k1) + spar(ik)*(ztetra(k2)-ztetra(k1))

      am(1,1) = sum(b*b)
      am(1,2) = sum(b*c)
      am(2,1) = sum(b*c)
      am(2,2) = sum(c*c)

      detam = am(1,1)*am(2,2) - am(1,2)*am(2,1)

      if (abs(detam) < eps10) then
        write (6,*) ' problem in berechne_koordinaten '
        write (6,*) ' gls zur koordinatentransformation nicht loesbar'
        call exit (1)
      end if

      amm1(1,1) = am(2,2) / detam
      amm1(1,2) = -am(1,2) / detam
      amm1(2,1) = -am(2,1) / detam
      amm1(2,2) = am(1,1) / detam

      do i=1,nrknot

        ik = iedge(i)
        k1 = kanten(1,ik)
        k2 = kanten(2,ik)

        x(1) = xtetra(k1) + spar(ik)*(xtetra(k2)-xtetra(k1)) 
        x(2) = ytetra(k1) + spar(ik)*(ytetra(k2)-ytetra(k1)) 
        x(3) = ztetra(k1) + spar(ik)*(ztetra(k2)-ztetra(k1)) 

        rhs(1) = sum(b*(x-orig))
        rhs(2) = sum(c*(x-orig))

        xtrian(i) = amm1(1,1)*rhs(1) + amm1(1,2)*rhs(2)
        ytrian(i) = amm1(2,1)*rhs(1) + amm1(2,2)*rhs(2)
      end do

      end subroutine berechne_koordinaten



      subroutine sort_ueberpruefen(ipunkt)
  
      implicit none

      integer, intent(inout), dimension(4) :: ipunkt
      real(dp), dimension(4,2)             :: punkt
      real(dp), dimension(2,2)             :: a
      real(dp), dimension(2)               :: b, t, cen
      real(dp)                             :: x1, x2
      real(dp)                             :: d1, d2, d3, d4, d, deta
      logical                              :: lsw1, lsw2
      integer                              :: ih, i
         
      punkt(1:4,1) = (/ (xtrian(ipunkt(i)), i=1,4) /)
      punkt(1:4,2) = (/ (ytrian(ipunkt(i)), i=1,4) /)

      cen(1)=sum(punkt(1:4,1))*0.25_dp
      cen(2)=sum(punkt(1:4,2))*0.25_dp

      d1 = sqrt(sum((punkt(1,:)-cen)**2))
      d2 = sqrt(sum((punkt(2,:)-cen)**2))
      d3 = sqrt(sum((punkt(3,:)-cen)**2))
      d4 = sqrt(sum((punkt(4,:)-cen)**2))
      d = 1._dp/max(min(d1,d2,d3,d4),eps10)

      lsw1=.true.
      lsw2=.true.
      do while (lsw1 .or. lsw2)

!  test p1-p2 and p3-p4
        a(1,1) = punkt(2,1) - punkt(1,1)
        a(2,1) = punkt(2,2) - punkt(1,2)
        a(1,2) = punkt(3,1) - punkt(4,1)
        a(2,2) = punkt(3,2) - punkt(4,2)

        b(1) = punkt(3,1) - punkt(1,1)
        b(2) = punkt(3,2) - punkt(1,2)
            
        a = a*d
        b = b*d
            
        deta = a(1,1)*a(2,2) - a(2,1)*a(1,2)
        if (abs(deta) < eps10) then
! p1-p2 und p3-p4 sind parallel
        lsw1 = .false.
        else
        x1 = (b(1)*a(2,2) - b(2)*a(1,2))/deta
        x2 = (b(2)*a(1,1) - b(1)*a(2,1))/deta
          if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .        (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 2 and 3
            t(1:2) = punkt(2,1:2)
            punkt(2,1:2) = punkt(3,1:2)
            punkt(3,1:2) = t(1:2)
          ih = ipunkt(2)
            ipunkt(2) = ipunkt(3)
            ipunkt(3) = ih
            lsw1 = .true.
          else
            lsw1 = .false.
          end if
        end if

!  test p1-p4 and p2-p3
        a(1,1) = punkt(4,1) - punkt(1,1)
        a(2,1) = punkt(4,2) - punkt(1,2)
        a(1,2) = punkt(2,1) - punkt(3,1)
        a(2,2) = punkt(2,2) - punkt(3,2)

        b(1) = punkt(2,1) - punkt(1,1)
        b(2) = punkt(2,2) - punkt(1,2)
            
        a = a*d
        b = b*d
            
        deta = a(1,1)*a(2,2) - a(2,1)*a(1,2)
        if (abs(deta) < eps10) then
! p1-p4 und p2-p3 sind parallel
        lsw2 = .false.
        else
        x1 = (b(1)*a(2,2) - b(2)*a(1,2))/deta
        x2 = (b(2)*a(1,1) - b(1)*a(2,1))/deta
      
          if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .        (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 3 and 4
            t(1:2) = punkt(3,1:2)
            punkt(3,1:2) = punkt(4,1:2)
            punkt(4,1:2) = t(1:2)
            ih = ipunkt(3)
            ipunkt(3) = ipunkt(4)
            ipunkt(4) = ih
            lsw2 = .true.
          else
            lsw2 = .false.
          end if
        end if  

      end do                 ! do while

      end subroutine sort_ueberpruefen



      subroutine bilde_dreiecke

      implicit none
      integer :: iknot(6), ipunkt(4)
      integer :: itet, icount, j, ih, i1, i2, i3, ico
      real(dp) :: ar, artri3

      if (.not.allocated(itetno)) allocate (itetno(ntri))
      itetno = 0

      ntrii = 0
      do itet=1,ntet
!        iknot(1:6) = (/ (icono(tetra_kanten(j,itet)),j=1,6) /)
!        icount = count(iknot > 0)
        icount=0
        iloop: do i=1,6
          ico = icono(tetra_kanten(i,itet))
          if (ico > 0) then
            do j=1,icount
              if (ico == iknot(j)) cycle iloop
            end do
            icount = icount + 1
            iknot(icount) = ico
          end if
        end do iloop

        if (icount == 3) then
! schnittgebilde ist ein Dreieck
          ntrii = ntrii + 1
          necke(1:3,ntrii) = iknot(1:3)
          i1=necke(1,ntrii)
          i2=necke(2,ntrii)
          i3=necke(3,ntrii)
!          isum=i1+i2+i3
!          imn=min(i1,i2,i3)
!          imx=max(i1,i2,i3)
!          imid=isum-imn-imx
!          if ((imn == imid) .or. (imx == imid)) then
!            ntrii = ntrii - 1
!            cycle
!          end if
          itetno(ntrii) = itet
          ar = artri3(xtrian(i1),ytrian(i1),0._dp,
     .                xtrian(i2),ytrian(i2),0._dp,
     .                xtrian(i3),ytrian(i3),0._dp)
          if (ar < 0._dp) then
! orientierung des dreiecks stimmt nicht, drehe die punkte 2 und 3
            ih = necke(2,ntrii)
            necke(2,ntrii) = necke(3,ntrii)
            necke(3,ntrii) = ih
          end if

        else if (icount ==4) then
! schnittgebilde ist ein Viereck
          ipunkt = iknot(1:4)

! sorge dafuer, dass die Linie p1,p2,p3,p4 ein Viereck bildet und
! sich nicht schneidet
          call sort_ueberpruefen (ipunkt)
! bilde erstes Dreieck p1, p2, p3
          ntrii = ntrii + 1
          necke(1:3,ntrii) = (/ ipunkt(1), ipunkt(2), ipunkt(3) /)   
          i1=necke(1,ntrii)
          i2=necke(2,ntrii)
          i3=necke(3,ntrii)
!          isum=i1+i2+i3
!          imn=min(i1,i2,i3)
!          imx=max(i1,i2,i3)
!          imid=isum-imn-imx
!          if ((imn < imid) .and. (imx > imid)) then
            itetno(ntrii) = itet
            ar = artri3(xtrian(i1),ytrian(i1),0._dp,
     .                  xtrian(i2),ytrian(i2),0._dp,
     .                  xtrian(i3),ytrian(i3),0._dp)
            if (ar < 0._dp) then
! orientierung des dreiecks stimmt nicht, drehe die punkte 2 und 3
              ih = necke(2,ntrii)
              necke(2,ntrii) = necke(3,ntrii)
              necke(3,ntrii) = ih
            end if
!          else
!            ntrii = ntrii-1
!          end if

! bilde zweites Dreieck p1, p3, p4
          ntrii = ntrii + 1
          necke(1:3,ntrii) = (/ ipunkt(1), ipunkt(3), ipunkt(4) /)   
          i1=necke(1,ntrii)
          i2=necke(2,ntrii)
          i3=necke(3,ntrii)
!          isum=i1+i2+i3
!          imn=min(i1,i2,i3)
!          imx=max(i1,i2,i3)
!          imid=isum-imn-imx
!          if ((imn == imid) .or. (imx == imid)) then
!            ntrii = ntrii - 1
!            cycle
!          end if
          itetno(ntrii) = itet
          ar = artri3(xtrian(i1),ytrian(i1),0._dp,
     .                xtrian(i2),ytrian(i2),0._dp,
     .                xtrian(i3),ytrian(i3),0._dp)
          if (ar < 0._dp) then
! orientierung des dreiecks stimmt nicht, drehe die punkte 2 und 3
            ih = necke(2,ntrii)
            necke(2,ntrii) = necke(3,ntrii)
            necke(3,ntrii) = ih
        end if
         
        end if       
      end do

      do itri=1,ntrii
        xcom(itri) = (xtrian(necke(1,itri)) + xtrian(necke(2,itri)) +
     .                xtrian(necke(3,itri))) / 3._dp
        ycom(itri) = (ytrian(necke(1,itri)) + ytrian(necke(2,itri)) +
     .                ytrian(necke(3,itri))) / 3._dp
      end do

      end subroutine bilde_dreiecke


      end subroutine rpscut
C ===== SOURCE: rpsout.f
C
C
      SUBROUTINE RPSOUT
C
C  ANZ: NUMBER OF CELLS
C  WRITE (17,...) LABELED CO-ORDINATES OF VERTICES
C  WRITE (18,...) LABELING NUMBER OF VERTICES FOR EACH CELL
C  WRITE (19,...) VALUE OF TALLY AT VERTEX (BY LABEL)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CRECH
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CLGIN
      USE CTRIG
      USE CTETRA
      USE CGRPTL
      USE CCONA

      IMPLICIT NONE

      REAL(DP) :: YWERT(2*NPTAL)
      INTEGER :: NCELL, NCO, IR, IP, IGR, IPPLG, NPUNKT, ISTS, IDIMP,
     .           NST, I, NSTAB, NRPS, IT, IF, IA, IE, ANZ
      integer :: ird1,ird2,ird3,ird4,ird5,ird6,ird7,ird8
      integer :: igroups, ipoints, ipl
      integer :: is, is1, anz2
      real(DP), ALLOCATABLE :: phelp(:,:),valcont(:), xcont(:), ycont(:)
      real(DP) :: xm, ym, length_m_p, xgeomin, xgeomax,
     .            ygeomin, ygeomax, epsrel, ahelp, vecax, vecay, vecbx,
     .            vecby,length_veca, length_vecb, del_eps,
     .            x1,x2,x3,x4,y1,y2,y3,y4,atri1,atri2
      integer :: icont, ipoint, iloop,ncont, ip_start, idel
      logical :: del_point

      TYPE(PPOINT), POINTER :: CUR
C
      allocate(valcont(iraps))
      valcont = huge(1._DP)
      igroups = 0
      I=0
      NSTAB=0
      CUR => FIRST_POINT
      IF (ASSOCIATED(CUR)) THEN
        DO WHILE (ASSOCIATED(CUR%NXTPNT))
          IF (CUR%NPL2D.EQ.0) THEN
C  START OF A NEW LINE
            NSTAB=NSTAB+1
          ELSEIF (CUR%NXTPNT%NPL2D.EQ.1) THEN
C  CONTINUATION OF A LINE
            NSTAB=NSTAB+1
          ENDIF
          CUR => CUR%NXTPNT
        END DO
      END IF
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
C  3D
      if ((levgeo.eq.1.and.nlrad.and.nlpol.and.nltor.and.nltrz) .or.
     .    (levgeo.eq.5.and..not.lrpscut).or.lraps3d) then
         WRITE(17,'(1X,A5,8X,A4,11X,A1,11X,A1,11X,A1,11X,A1)') '-1111',
     .        'NPCO','1','3','1','1'
C  2D
      else
         WRITE(17,'(1X,A5,8X,A4,11X,A1,11X,A1,11X,A1,11X,A1)') '-1111',
     .        'NPCO','1','2','1','1'
      endif
C
      WRITE(19,'(1X,A5,8X,A4,50(9X,I3))') '-1111',
     .        'NPST',1,IRAPS,1,(1,I=1,IRAPS)
C
      if (levgeo.eq.1.and.nlrad.and.nlpol.and.nltor.and.nltrz) then
         ANZ=(nr1st-1)*(np2nd-1)*(nt3rd-1)
         i=0
         do ir=1,nr1st
            do ip=1,np2nd
               do it=1,nt3rd
                  DO IF=1,IRAPS
                     READ(60+IF,*) YWERT(IF)
                  enddo
                  i = i+1
                  WRITE(19,'(I6,1P,50E12.4)') I,(YWERT(IF),IF=1,IRAPS)
                  WRITE(17,'(I6,1P,3E12.4)')
     .                       I,rsurf(ir),psurf(ip),zsurf(it)
               enddo
            enddo
         enddo

         WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
         WRITE(18,'(A14,I6,5X,A1)') 'HEXE8        1',ANZ,'8'

         do ir=1,nr1st-1
            do ip=1,np2nd-1
               do it=1,nt3rd-1
                  IRD1=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IRD2=IR+1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IRD3=IR+1+((IP)+(IT-1)*NP2T3)*NR1P2
                  IRD4=IR+((IP)+(IT-1)*NP2T3)*NR1P2
                  IRD5=IR+((IP-1)+(IT)*NP2T3)*NR1P2
                  IRD6=IR+1+((IP-1)+(IT)*NP2T3)*NR1P2
                  IRD7=IR+1+((IP)+(IT)*NP2T3)*NR1P2
                  IRD8=IR+((IP)+(IT)*NP2T3)*NR1P2
                  WRITE(18,'(1X,A1,8I10)') '0',ird1,ird2,ird3,ird4,
     .                                        ird5,ird6,ird7,ird8
               enddo
            enddo
         enddo

      ELSEIF (LEVGEO.EQ.1.AND.LPTORR) THEN
        ANZ = 0
        IT=IPTORR
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            IF (NSTGRD(NCELL).EQ.0) ANZ=ANZ+1
          ENDDO
        ENDDO
C
        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'QUAM4        1',ANZ,'4'

C  DO 100, DO 110: RETAIN SAME SEQUENCE FOR READING FROM FORT(60+IF)
C                  AS IT WAS THE CASE FOR WRITING (RPSCOL,RPSVEC)

        I=0
        IT=IPTORR
        DO 100 IR=1,NR1ST
           DO 110 IP=1,NP2ND

C  FORT 60+IF WAS WRITTEN IN RPSCOL OR RPSVEC IN SAME DO LOOPS
             DO 105 IF=1,IRAPS
               READ (60+IF,*) YWERT(IF)
105          CONTINUE

             IF (IP .NE. NP2ND) THEN
               I = I + 1
               WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR)*FCABS1(1),
     .                                      PSURF(IP)*FCABS2(1)
C  EXCLUDE DEAD CELLS ON FORT.18
               NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
               IF (IR.LT.NR1ST.AND.NSTGRD(NCELL).EQ.0) THEN
                 WRITE(18,'(1X,A1,4I10)') '0',
     .                     (IR-1)*NP2ND+IP,
     .                     (IR-1)*NP2ND+IP+1,
     .                     IR*NP2ND+IP+1,
     .                     IR*NP2ND+IP
               ENDIF
             ENDIF
             IF (IP .NE. NP2ND) THEN
               WRITE(19,'(I6,1P,50E12.4)') I,  (YWERT(IF),IF=1,IRAPS)
             ELSE
               WRITE(19,'(I6,1P,50E12.4)') I+1,(YWERT(IF),IF=1,IRAPS)
             ENDIF
110        CONTINUE
           I = I + 1
           WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR)*FCABS1(1),
     .                                  PSURF(NP2ND)*FCABS2(1)
100     CONTINUE
        NCO=I
C
      ELSEIF (LEVGEO.LE.2.AND.LPPOLR) THEN
        ANZ = 0
        IP=IPPOLR
        DO IR=1,NR1STM
          DO IT=1,NT3RDM
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            IF (NSTGRD(NCELL).EQ.0) ANZ=ANZ+1
          ENDDO
        ENDDO
C
        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'QUAM4        1',ANZ,'4'

        I=0
        IP=IPPOLR
        DO 1100 IR=1,NR1ST
          DO 2100 IT=1,NT3RD

C  FORT 60+IF WAS WRITTEN IN RPSCOL OR RPSVEC IN SAME DO LOOPS
            DO 2105 IF=1,IRAPS
              READ (60+IF,*) YWERT(IF)
2105        CONTINUE

            IF (IT .NE. NT3RD) THEN
              I = I + 1
              WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR),ZSURF(IT)
C  EXCLUDE DEAD CELLS ON FORT.18
              NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
              IF (IR.LT.NR1ST.AND.NSTGRD(NCELL).EQ.0) THEN
                WRITE(18,'(1X,A1,4I10)') '0',
     .                    (IR-1)*NT3RD+IT,
     .                    (IR-1)*NT3RD+IT+1,
     .                     IR*NT3RD+IT+1,
     .                     IR*NT3RD+IT
              ENDIF
            ENDIF
            IF (IT .NE. NT3RD) THEN
              WRITE(19,'(I6,1P,50E12.4)') I,  (YWERT(IF),IF=1,IRAPS)
            ELSE
              WRITE(19,'(I6,1P,50E12.4)') I+1,(YWERT(IF),IF=1,IRAPS)
            ENDIF
2100      CONTINUE
          I = I + 1
          WRITE(17,'(I6,1P,2E12.4)') I,RSURF(IR),ZSURF(NT3RD)
1100    CONTINUE
        NCO=I
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
        ANZ = 0
        IT=IPTORR
        DO  IR=1,NR1ST-1
          DO  IPPLG=1,NPPLG
            DO  IP=NPOINT(1,IPPLG),NPOINT(2,IPPLG)-1
              NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
              IF (NSTGRD(NCELL).EQ.0) ANZ=ANZ+1
            ENDDO
          ENDDO
        ENDDO

        if (lraps3d.and.lr3dcon) then
          WRITE(18,'(A3,2X,A14,45X,I8)')
     .         'PSS','1BEISPIELDATEN',(ANZ+nstab)*(iplane-1)
        else
          WRITE(18,'(A3,2X,A14,45X,I8)')'PSS','1BEISPIELDATEN',
     .                                   (ANZ+NSTAB)*iplane
        endif

        I=0
        ipoints=0
        do ipl=0,iplane-1
          DO IF=1,IRAPS
            open (60+IF)
          enddo

          if (lraps3d.and.lr3dcon) then
            if (ipl < iplane-1)
     .        WRITE(18,'(A8,2I6,5X,A1)') 'HEXE8        ',ipl+1,ANZ,'8'
            igroups = iplane-1
          else
            WRITE(18,'(A8,2I6,5X,A1)') 'QUAM4        ',ipl+1,ANZ,'4'
            igroups = iplane
          endif

          IT=IPTORR
          xgeomin = huge(1._dp)
          xgeomax = -xgeomin
          ygeomin = huge(1._dp)
          ygeomax = -ygeomin
          DO 10 IR=1,NR1ST
            NPUNKT = 0
            DO 11 IPPLG=1,NPPLG
              NPUNKT = NPUNKT+NPOINT(2,IPPLG)-NPOINT(1,IPPLG)+1
 11         CONTINUE
            DO 20 IPPLG=1,NPPLG
              xgeomin = min(xgeomin,minval(
     .                  xpol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              xgeomax = max(xgeomax,maxval(
     .                  xpol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              ygeomin = min(ygeomin,minval(
     .                  ypol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              ygeomax = max(ygeomax,maxval( 
     .                  ypol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              DO 30 IP=NPOINT(1,IPPLG),NPOINT(2,IPPLG)

C  FORT 60+IF WAS WRITTEN IN RPSCOL OR RPSVEC IN SAME DO LOOPS
                DO 25 IF=1,IRAPS
                  READ (60+IF,*) YWERT(IF)
 25             CONTINUE

                IF (IP .NE. NPOINT(2,IPPLG)) THEN
                  I = I + 1
                  if (lraps3d.and.nltra) then
                    WRITE(17,'(I6,1P,3E12.4)')
     .                    I,XPOL(IR,IP)*cos(ipl*rapsdel),YPOL(IR,IP),
     .                      xpol(ir,ip)*sin(ipl*rapsdel)
                  elseif (lraps3d.and.nltrz) then
                    WRITE(17,'(I6,1P,3E12.4)')
     .                    I,XPOL(IR,IP),YPOL(IR,IP),ipl*rapsdel
                  else
                    WRITE(17,'(I6,1P,2E12.4)') I,XPOL(IR,IP),YPOL(IR,IP)
                  endif
C  EXCLUDE DEAD CELLS ON FORT.18
                  NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IF (IR.LT.NR1ST.AND.NSTGRD(NCELL).EQ.0) THEN
                    if (lraps3d.and.lr3dcon) then
                      if (ipl < iplane-1) then
                        WRITE(18,'(1X,A1,8I10)') '0',
     .                       (IR-1)*NPUNKT+IP+ipl*ipoints,
     .                       (IR-1)*NPUNKT+IP+1+ipl*ipoints,
     .                        IR*NPUNKT+IP+1+ipl*ipoints,
     .                        IR*NPUNKT+IP+ipl*ipoints,
     .                       (IR-1)*NPUNKT+IP+(ipl+1)*ipoints,
     .                       (IR-1)*NPUNKT+IP+1+(ipl+1)*ipoints,
     .                        IR*NPUNKT+IP+1+(ipl+1)*ipoints,
     .                        IR*NPUNKT+IP+(ipl+1)*ipoints
                      endif
                    else
                      WRITE(18,'(1X,A1,4I10)') '0',
     .                     (IR-1)*NPUNKT+IP+ipl*ipoints,
     .                     (IR-1)*NPUNKT+IP+1+ipl*ipoints,
     .                      IR*NPUNKT+IP+1+ipl*ipoints,
     .                      IR*NPUNKT+IP+ipl*ipoints
                    endif
                  ENDIF
                ENDIF
                IF (IP .NE. NPOINT(2,IPPLG)) THEN
                  WRITE(19,'(I6,1P,50E12.4)') I,(YWERT(IF),IF=1,IRAPS)
                ELSE
                  WRITE(19,'(I6,1P,50E12.4)') I+1,(YWERT(IF),IF=1,IRAPS)
                ENDIF
30            CONTINUE
              I = I + 1
              if (lraps3d.and.nltra) then
                WRITE(17,'(I6,1P,3E12.4)')
     .                I,XPOL(IR,NPOINT(2,IPPLG))*cos(ipl*rapsdel),
     .                  YPOL(IR,NPOINT(2,IPPLG)),
     .                  xpol(ir,NPOINT(2,IPPLG))*sin(ipl*rapsdel)
              elseif (lraps3d.and.nltrz) then
                WRITE(17,'(I6,1P,3E12.4)')
     .                I,XPOL(IR,NPOINT(2,IPPLG)),
     .                  YPOL(IR,NPOINT(2,IPPLG)),ipl*rapsdel
              else
                WRITE(17,'(I6,1P,2E12.4)') I,XPOL(IR,NPOINT(2,IPPLG)),
     .                                       YPOL(IR,NPOINT(2,IPPLG))
              endif
20          CONTINUE
10        CONTINUE
          if (ipl.eq.0) ipoints=i
          DO IF=1,IRAPS
            close (60+IF)
          enddo
        enddo
        NCO=I
C
      ELSEIF ((LEVGEO.EQ.4.AND.LPTORR) .OR.
     .        (LEVGEO.EQ.5.AND.LRPSCUT)) THEN
C TO BE DONE: NSTGRD.NE.0 AUSBLENDEN, ANZ NEU BERECHENEN.
        IT=IPTORR
        ANZ=NTRII
        xgeomin = minval(xtrian(1:nrknot))
        xgeomax = maxval(xtrian(1:nrknot))
        ygeomin = minval(ytrian(1:nrknot))
        ygeomax = maxval(ytrian(1:nrknot))
        if (lraps3d.and.lr3dcon) then
          WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',
     .                                   (ANZ+NSTAB)*(iplane-1)
        else
          WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',
     .                                   (ANZ+NSTAB)*iplane
        endif
        do ipl=0,iplane-1
          DO IF=1,IRAPS
            open(60+IF)
          enddo
          DO 40 I=1,NRKNOT
            DO 50 IF=1,IRAPS
              READ(60+IF,*) YWERT(IF)
              if (ywert(if) < valcont(if)) valcont(if) = ywert(if)
50          CONTINUE
            WRITE(19,'(I6,1P,50E12.4)') I+ipl*nrknot,
     .                                  (YWERT(IF),IF=1,IRAPS)
            if (lraps3d.and.nltra) then
              WRITE(17,'(I6,1P,3E12.4)') I+ipl*nrknot,
     .              XTRIAN(I)*cos(ipl*rapsdel),
     .              YTRIAN(I),
     .              XTRIAN(I)*sin(ipl*rapsdel)
            elseif (lraps3d.and.nltrz) then
              WRITE(17,'(I6,1P,3E12.4)') I+ipl*nrknot,
     .              XTRIAN(I),YTRIAN(I),ipl*rapsdel
            else
              WRITE(17,'(I6,1P,2E12.4)') I,XTRIAN(I),YTRIAN(I)
            endif
40        CONTINUE
          DO IF=1,IRAPS
            close(60+IF)
          enddo
          if (lraps3d.and.lr3dcon) then
            if (ipl < iplane-1) then
              WRITE(18,'(A8,2I6,5X,A1)') 'PENTA6  ',ipl+1,ANZ,'6'

              DO I=1,NTRII
                WRITE(18,'(1X,A1,6I10)') '0',NECKE(1,I)+ipl*nrknot,
     .                                      NECKE(2,I)+ipl*nrknot,
     .                                      NECKE(3,I)+ipl*nrknot,
     .                                      NECKE(1,I)+(ipl+1)*nrknot,
     .                                      NECKE(2,I)+(ipl+1)*nrknot,
     .                                      NECKE(3,I)+(ipl+1)*nrknot
              enddo
            endif
            igroups = iplane-1
          else
            WRITE(18,'(A8,2I6,5X,A1)') 'TRIM3   ',ipl+1,ANZ,'3'

            DO 60 I=1,NTRII
              WRITE(18,'(1X,A1,3I10)') '0',NECKE(1,I)+ipl*nrknot,
     .                                    NECKE(2,I)+ipl*nrknot,
     .                                    NECKE(3,I)+ipl*nrknot
60          CONTINUE
            igroups = iplane
          endif
        enddo
        NCO=NRKNOT*iplane
C
      ELSEIF (LEVGEO.EQ.5.AND..NOT.LRPSCUT) THEN
C TO BE DONE: NSTGRD.NE.0 AUSBLENDEN, ANZ NEU BERECHENEN.
        anz=ntet-ntet_collaps
        do i=1,ncoor
          do if=1,iraps
            read(60+if,*) ywert(if)
          enddo
          WRITE(19,'(I6,1P,50E12.4)') I,(YWERT(IF),IF=1,IRAPS)
          WRITE(17,'(I6,1P,3E12.4)') I,XTETRA(I),YTETRA(I),ZTETRA(I)
        enddo
        WRITE(18,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18,'(A14,I6,5X,A1)') 'TET4         1',ANZ,'4'
        do i=1,ntet
          if (ntbar(1,i) >= 0)
     .      WRITE(18,'(1X,A1,4I10)') '0',NTECK(1,I),NTECK(2,I),
     .                                  NTECK(3,I),NTECK(4,I)
        enddo
      ELSE
      ENDIF
C
      if (nstab > 0) then
c ist das jemals getestet worden ???
        IGR=igroups+1
        do ipl=0,iplane-1
          if (lraps3d.and.lr3dcon) then
             if (ipl < iplane-1) then
                 WRITE(18,'(A8,2I6,5X,A1)') 'QUAM4   ',IGR+ipl,NSTAB,'4'
                 igroups = igroups+1
              endif
          else
            WRITE(18,'(A8,2I6,5X,A1)') 'FLA2    ',IGR+ipl,NSTAB,'2'
            igroups = igroups+1
          endif
C
          CUR => FIRST_POINT
          IF (ASSOCIATED(CUR)) THEN
            DO WHILE (ASSOCIATED(CUR%NXTPNT))
              IF (CUR%NPL2D.EQ.0) THEN
C  START OF A NEW LINE
                if (lraps3d.and.nltra) then
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                   CUR%XPL2D*cos(ipl*rapsdel),
     .                   CUR%YPL2D,
     .                   CUR%XPL2D*sin(ipl*rapsdel)
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                   CUR%NXTPNT%XPL2D*cos(ipl*rapsdel),
     .                   CUR%NXTPNT%YPL2D,
     .                   CUR%NXTPNT%XPL2D*sin(ipl*rapsdel)
                elseif (lraps3d.and.nltrz) then
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,CUR%XPL2D,CUR%YPL2D,
     .                                        ipl*rapsdel
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                   CUR%NXTPNT%XPL2D,
     .                   CUR%NXTPNT%YPL2D,
     .                   ipl*rapsdel
                else
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,2E12.4)') NCO,CUR%XPL2D,CUR%YPL2D
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,2E12.4)')
     .                   NCO,CUR%NXTPNT%XPL2D,CUR%NXTPNT%YPL2D
                endif
                if (lraps3d.and.lr3dcon) then
                  if (ipl < iplane-1)
     .              WRITE (18,'(1X,A1,4I10)') '0',
     .                     NCO-1,NCO,nco+ipl*2*nstab,
     .                     nco-1+ipl*2*nstab
                else
                  WRITE (18,'(1X,A1,2I10)') '0',NCO-1,NCO
                endif
              ELSEIF (CUR%NXTPNT%NPL2D.EQ.1) THEN
C  CONTINUATION OF A LINE
                if (lraps3d.and.nltra) then
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                   CUR%XPL2D*cos(ipl*rapsdel),
     .                   CUR%YPL2D,
     .                   CUR%XPL2D*sin(ipl*rapsdel)
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                   CUR%NXTPNT%XPL2D*cos(ipl*rapsdel),
     .                   CUR%NXTPNT%YPL2D,
     .                   CUR%NXTPNT%XPL2D*sin(ipl*rapsdel)
                elseif (lraps3d.and.nltrz) then
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,CUR%XPL2D,
     .                   CUR%YPL2D,
     .                   ipl*rapsdel
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                   CUR%NXTPNT%XPL2D,
     .                   CUR%NXTPNT%YPL2D,
     .                   ipl*rapsdel
                else
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,2E12.4)') NCO,CUR%XPL2D,CUR%YPL2D
                  NCO=NCO+1
                  WRITE (17,'(I6,1P,2E12.4)')
     .                   NCO,CUR%NXTPNT%XPL2D,CUR%NXTPNT%YPL2D
                endif
                if (lraps3d.and.lr3dcon) then
                  if (ipl < iplane-1)
     .               WRITE (18,'(1X,A1,4I10)') '0',
     .                      NCO-1,NCO,nco+ipl*2*nstab,
     .                      nco-1+ipl*2*nstab
                else
                  WRITE (18,'(1X,A1,2I10)') '0',NCO-1,NCO
                endif
              ENDIF
              CUR => CUR%NXTPNT
            END DO
          END IF
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
                    if (lraps3d.and.nltra) then
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP)*cos(ipl*rapsdel),
     .                       YPOL(IR,IP),
     .                       XPOL(IR,IP)*sin(ipl*rapsdel)
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP+1)*cos(ipl*rapsdel),
     .                       YPOL(IR,IP+1),
     .                       XPOL(IR,IP+1)*sin(ipl*rapsdel)
                    elseif (lraps3d.and.nltrz) then
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP),
     .                       YPOL(IR,IP),
     .                       ipl*rapsdel
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP+1),
     .                       YPOL(IR,IP+1),
     .                       ipl*rapsdel
                    else
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR,IP),YPOL(IR,IP)
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR,IP+1),YPOL(IR,IP+1)
                    endif
                    if (lraps3d.and.lr3dcon) then
                      if (ipl < iplane-1)
     .                  WRITE (18,'(1X,A1,4I10)') '0',
     .                         NCO-1,NCO,nco+ipl*2*nstab,
     .                         nco-1+ipl*2*nstab
                    else
                      WRITE (18,'(1X,A1,2I10)') '0',NCO-1,NCO
                    endif
                  END DO

                ELSEIF (IDIMP ==2) THEN
                  IP=INUMP(ISTS,IDIMP)
                  IA=IRPTA(ISTS,1)
                  IE=IRPTE(ISTS,1)
                  DO IR=IA,IE-1
                    if (lraps3d.and.nltra) then
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP)*cos(ipl*rapsdel),
     .                       YPOL(IR,IP),
     .                       XPOL(IR,IP)*sin(ipl*rapsdel)
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR+1,IP)*cos(ipl*rapsdel),
     .                       YPOL(IR+1,IP),
     .                       XPOL(IR+1,IP)*sin(ipl*rapsdel)
                    elseif (lraps3d.and.nltrz) then
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP),
     .                       YPOL(IR,IP),
     .                       ipl*rapsdel
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR+1,IP),
     .                       YPOL(IR+1,IP),
     .                       ipl*rapsdel
                    else
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR,IP),YPOL(IR,IP)
                      NCO=NCO+1
                      WRITE (17,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR+1,IP),YPOL(IR+1,IP)
                    endif
                    if (lraps3d.and.lr3dcon) then
                      if (ipl < iplane-1)
     .                  WRITE (18,'(1X,A1,4I10)') '0',
     .                         NCO-1,NCO,nco+ipl*2*nstab,
     .                         nco-1+ipl*2*nstab
                    else
                      WRITE (18,'(1X,A1,2I6)') '0',NCO-1,NCO
                    endif
                  END DO
                END IF
              END IF
            END DO
          END IF
        enddo
      endif
C
      if (lraps3d  .and. ((levgeo == 3) .or. (levgeo == 4)) ) then
        if (ncontour > 0) then
          ALLOCATE (phelp(maxval(nconpoint(1:ncontour))-1,2))
          ALLOCATE (xcont(maxval(nconpoint(1:ncontour))))
          ALLOCATE (ycont(maxval(nconpoint(1:ncontour))))
          igr=igroups+1
          epsrel = (xgeomax-xgeomin+ygeomax-ygeomin)/2.*eps5
          do icont=1,ncontour
            write(65,*) 'xcontour ycontour icont: ',icont
            do ipoint = 1,nconpoint(icont)
               write(65,'(2es12.4)') xcontour(ipoint,icont),
     .              ycontour(ipoint,icont)
            enddo
            write(65,*)
c           bereinigte contour erstellen
            xcont(1:nconpoint(icont))=xcontour(1:nconpoint(icont),icont)
            ycont(1:nconpoint(icont))=ycontour(1:nconpoint(icont),icont)
            ncont = nconpoint(icont)
            del_point = .true.
            ip_start = 1
            do while (del_point)
               do ipoint=ip_start,ncont-1
                  vecax = xcont(ipoint)-xcont(ipoint+1)
                  vecay = ycont(ipoint)-ycont(ipoint+1)
                  length_veca = sqrt(vecax**2+vecay**2)
                  if (length_veca <
     .                (100./10**ceiling(log10(abs(xcont(ipoint)))) +
     .                 100./10**ceiling(log10(abs(xcont(ipoint)))))/2.)
     .            then
                     if (ipoint == 1) then
                        x1 = xcont(ncont-1)
                        y1 = ycont(ncont-1)
                     else
                        x1 = xcont(ipoint-1)
                        y1 = ycont(ipoint-1)
                     endif
                     x2 = xcont(ipoint)
                     y2 = ycont(ipoint)
                     x3 = xcont(ipoint+1)
                     y3 = ycont(ipoint+1)
                     if (ipoint == (ncont-1)) then
                        x4 = xcont(ipoint+2)
                        y4 = ycont(ipoint+2)
                     else
                        x4 = xcont(2)
                        y4 = ycont(2)
                     endif
                     atri1 = x1*(y3-y2)+x3*(y2-y1)+x2*(y1-y3)
                     atri2 = x4*(y3-y2)+x3*(y2-y4)+x2*(y4-y3)
                     if ((atri1 < -eps5) .and. (atri2 < -eps5)) then
                        ncont = ncont - 1
                        xcont(ipoint) = (xcont(ipoint)+
     .                                   xcont(ipoint+1))/2.
                        ycont(ipoint) = (ycont(ipoint)+
     .                                   ycont(ipoint+1))/2.
                        do idel=ipoint+1, ncont
                           xcont(idel) = xcont(idel+1)
                           ycont(idel) = ycont(idel+1)
                        enddo
                        xcont(ncont+1) = 0.
                        ycont(ncont+1) = 0.
                        ip_start=max(1,ipoint-1)
                        exit
                     elseif (((atri1 > 0) .and. (atri2 < 0)) .or.
     .                       (abs(atri1) < eps5)) then
                        ncont = ncont - 1
                        do idel=ipoint, ncont
                           xcont(idel) = xcont(idel+1)
                           ycont(idel) = ycont(idel+1)
                        enddo
                        xcont(ncont+1) = 0.
                        ycont(ncont+1) = 0.
                        ip_start=max(1,ipoint-1)
                        exit
                     elseif (((atri1 < 0) .and. (atri2 > 0)) .or.
     .                       (abs(atri2) < eps5)) then
                        ncont = ncont - 1
                        do idel=ipoint+1, ncont
                           xcont(idel) = xcont(idel+1)
                           ycont(idel) = ycont(idel+1)
                        enddo
                        xcont(ncont+1) = 0.
                        ycont(ncont+1) = 0.
                        ip_start=ipoint
                        exit
                     endif
                  endif
                  if (ipoint == (ncont-1)) then
                     del_point = .false.
                  endif
               enddo
             enddo
             do ipoint=1,NCONT-1
c     punkt berechnen
               if (ipoint == 1) then
                  vecax=xcont(ncont-1)-xcont(ipoint)
                  vecay=ycont(ncont-1)-ycont(ipoint)
               else
                  vecax=xcont(ipoint-1)-xcont(ipoint)
                  vecay=ycont(ipoint-1)-ycont(ipoint)
               endif
               vecbx = xcont(ipoint+1)-xcont(ipoint)
               vecby = ycont(ipoint+1)-ycont(ipoint)
               length_veca = sqrt(vecax**2+vecay**2)
               length_vecb = sqrt(vecbx**2+vecby**2)
               vecax = vecax/length_veca
               vecay = vecay/length_veca
               vecbx = vecbx/length_vecb
               vecby = vecby/length_vecb
               del_eps = min(epsrel,length_veca/10.,length_vecb/10.)
               xm = xcont(ipoint) + vecax + vecbx
               ym = ycont(ipoint) + vecay + vecby
               length_m_p = sqrt((xm-xcont(ipoint))**2+
     .              (ym-ycont(ipoint))**2)
               ahelp = xm*(ycont(ipoint+1)-ycont(ipoint)) +
     .              xcont(ipoint+1) * (ycont(ipoint)-ym) +
     .              xcont(ipoint) * (ym - ycont(ipoint+1))
c     punkt in richtung m verschieben
               if ((length_m_p < eps10) .or. (abs(ahelp) < eps10)) then
                  xm = -1.*ycont(ipoint)+ycont(ipoint+1)
                  ym = xcont(ipoint)-xcont(ipoint+1)
                  length_m_p = sqrt(xm**2+ym**2)
                  phelp(ipoint,1) = xcont(ipoint)+xm/length_m_p*del_eps
                  phelp(ipoint,2) = ycont(ipoint)+ym/length_m_p*del_eps
               iloop = 1
               do while (
     .              (abs(xcont(ipoint)-phelp(ipoint,1))<
     .              100./10**ceiling(log10(abs(xcont(ipoint)))))
     .              .and.
     .              (abs(ycont(ipoint)-phelp(ipoint,2))<
     .              100./10**ceiling(log10(abs(ycont(ipoint))))))
                  iloop = iloop + 1
                  phelp(ipoint,1) = xcont(ipoint)+xm/
     .                 length_m_p*del_eps*iloop
                  phelp(ipoint,2) = ycont(ipoint)+ym/
     .                 length_m_p*del_eps*iloop
               enddo
               else
               phelp(ipoint,1) = xcont(ipoint)+
     .              ahelp/abs(ahelp+eps30) *
     .              (xm-xcont(ipoint))/length_m_p*del_eps
               phelp(ipoint,2) = ycont(ipoint)+
     .              ahelp/abs(ahelp+eps30) *
     .              (ym-ycont(ipoint))/length_m_p*del_eps
               iloop = 1
               do while (
     .              (abs(xcont(ipoint)-phelp(ipoint,1))<
     .              100./10**ceiling(log10(abs(xcont(ipoint)))))
     .              .and.
     .              (abs(ycont(ipoint)-phelp(ipoint,2))<
     .              100./10**ceiling(log10(abs(ycont(ipoint))))))
                  iloop = iloop + 1
                  phelp(ipoint,1) = xcont(ipoint)+
     .                 ahelp/abs(ahelp+eps30) *
     .                 (xm-xcont(ipoint))/length_m_p*
     .                 del_eps*iloop
                  phelp(ipoint,2) = ycont(ipoint)+
     .                 ahelp/abs(ahelp+eps30) *
     .                 (ym-ycont(ipoint))/length_m_p*
     .                 del_eps*iloop
               enddo
             endif
            enddo
            do ipl=0,iplane-1
               write(65,*) 'origx origy neux neuy'
               do ipoint=1,ncont-1
                  write(65,'(4es12.4)')  xcont(ipoint), ycont(ipoint),
     .                 phelp(ipoint,1),phelp(ipoint,2)
c     punkte schreiben
                  if (nltra) then
                     WRITE(17,'(I6,1P,3E12.4)') nco+ipoint,
     .                    phelp(Ipoint,1)*cos(ipl*rapsdel),
     .                    phelp(Ipoint,2),
     .                    phelp(Ipoint,1)*sin(ipl*rapsdel)
                     WRITE(19,'(I6,1P,50E12.4)') nco+ipoint,
     .                                          (valcont(IF),IF=1,IRAPS)
                  elseif (nltrz) then
                     WRITE(17,'(I6,1P,3E12.4)') nco+ipoint,
     .                    phelp(Ipoint,1),
     .                    phelp(Ipoint,2),
     .                    ipl*rapsdel
                     WRITE(19,'(I6,1P,50E12.4)') nco+ipoint,
     .                                          (valcont(IF),IF=1,IRAPS)
                  endif
               enddo
               if (ipl < iplane-1) then
c     elementgruppe anfangen
                  WRITE(18,'(A8,2I6,5X,A1)')
     .                 'QUAM4   ',IGR+ipl,ncont-1,'4'
                  do ipoint=1,ncont-1
c     element schreiben
                     if (ipoint < ncont-1) then
                        WRITE(18,'(1X,A1,4I10)') '0',
     .                       nco+ipoint,nco+ipoint+1,
     .                       nco+ipoint+1+ncont-1,
     .                       nco+ipoint+ncont-1
                     else
                        WRITE(18,'(1X,A1,4I10)') '0',
     .                       nco+ipoint,nco+1,
     .                       nco+1+ncont-1,
     .                       nco+ipoint+ncont-1
                     endif
                  enddo
               endif
               nco = nco + ncont-1
            enddo
            igr=igr+iplane-1
          enddo
          deallocate(phelp)
          deallocate(xcont)
          deallocate(ycont)
        endif
      endif

      WRITE(17,'(1X,A5,8X,A3,12X,A1,11X,A1,11X,A1,11X,A1)') '-9999',
     .           'FIN','0','0','0','0'
      WRITE(19,'(1X,A5,8X,A3,50(11X,I1))') '-9999',
     .           'FIN',(0,IF=1,IRAPS)
      deallocate(valcont)

      RETURN
      END
C ===== SOURCE: rpsvec.f
C  april 2006:  levgeo=1 option and nlpol: added
C
      SUBROUTINE RPSVEC (AORIG,BORIG,IBLD,ICURV,
     .                   IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A DATASET USED BY RAPS TO PRODUCE
C  A VECTORFIELD PLOT
C
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ONTO VERTICES IS PERFORMED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTRIG
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: AORIG(*), BORIG(*)
      REAL(DP), INTENT(IN) :: XX(*),YY(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV, IXX, IYY
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: XL, XM, BETRAG
      REAL(DP), ALLOCATABLE :: YWERT(:,:), YWERT1(:,:),
     .                       ZWERT(:,:), ZWERT1(:,:)
      INTEGER :: I, IP, IPART, IA, IB, IC, J, K, IT, LENCH, IERR, IR,
     .           NRAPS2, NVPLOT
      INTEGER :: ZUORD(NKNOT,0:20,2)
      REAL(SP) :: XY(800)
      REAL(SP) :: YH
      CHARACTER(17) :: CH
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
        ALLOCATE (YWERT1(NRAD,1))
        ALLOCATE (ZWERT1(NRAD,1))
        CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
        CALL CELINT(BORIG,ZWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
      ELSEIF (LEVGEO.LE.3) THEN
        ALLOCATE (YWERT(N1ST,N2ND+N3RD))
        ALLOCATE (ZWERT(N1ST,N2ND+N3RD))
        CALL CELINT(AORIG,YWERT,LOGL,IBLD,ICURV,N1ST,IERR)
        CALL CELINT(BORIG,ZWERT,LOGL,IBLD,ICURV,N1ST,IERR)
      ELSE
        WRITE (IUNOUT,*) 'UNWRITTEN OPTION IN RPSVEC '
        WRITE (IUNOUT,*) 'CELINT CANNOT BE CARRIED OUT'
        IERR=1
      ENDIF
      IF (IERR.GT.0) THEN
        IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
        IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
        IF (ALLOCATED(ZWERT)) DEALLOCATE (ZWERT)
        IF (ALLOCATED(ZWERT1)) DEALLOCATE (ZWERT1)
        RETURN
      END IF
C
C  SCALE VECTORS TO UNIT LENGTH
C
C  PROJECTION INTO X,Z PLANE
      IF (LEVGEO.LE.2.AND.NLTOR
     .     .AND.LPPOL3(IBLD)) THEN
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
C
C  PROJECTION INTO X,Y PLANE,  CARTHESIAN
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL
     .     .AND.LPTOR3(IBLD)) THEN
        DO 1101 IR=1,NR1ST
          DO 3101 IP=1,NP2ND
            IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
            IF (ZMI.NE.666.) ZWERT(IR,IP)=MAX(ZWERT(IR,IP),ZMI)
            IF (ZMA.NE.666.) ZWERT(IR,IP)=MIN(ZWERT(IR,IP),ZMA)
            BETRAG = SQRT(YWERT(IR,IP)**2+ZWERT(IR,IP)**2)
            WRITE (NRAPS,*) BETRAG
            YWERT(IR,IP)=YWERT(IR,IP)/(BETRAG+1.D-20)
            ZWERT(IR,IP)=ZWERT(IR,IP)/(BETRAG+1.D-20)
3101      CONTINUE
1101    CONTINUE
C
C  PROJECTION INTO X,Y PLANE, POLAR OR GENERAL CURVILINEAR (POLYGON)
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL
     .        .AND.LPTOR3(IBLD)) THEN
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
        WRITE (iunout,*) 'UNWRITTEN OPTION IN RPSVEC: PLOT ABANDONNED '
        WRITE (iunout,*) 'SCALING OF VECTORS FAILED  '
        IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
        IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
        IF (ALLOCATED(ZWERT)) DEALLOCATE (ZWERT)
        IF (ALLOCATED(ZWERT1)) DEALLOCATE (ZWERT1)
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

      IF ((LEVGEO.EQ.1).AND.NLPOL
     .     .AND.LPTOR3(IBLD)) THEN
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
          IPLOOP1: DO IP=1,NP2ND
C  RIGHT HEMISPHERE
            IF (IR.GT.1) THEN
C  LOWER RIGHT QUADRANT
              IF (IP.GT.1) THEN
C  SIDE 1
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR),PSURF(IP-1),
     .                       RHOSRF(IR-1),PSURF(IP-1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
C  SIDE 2
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR-1),PSURF(IP-1),
     .                       RHOSRF(IR-1),PSURF(IP),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
              ENDIF
C  UPPER RIGHT QUADRANT
              IF (IP.LT.NP2ND) THEN
C  SIDE 3
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR-1),PSURF(IP),
     .                       RHOSRF(IR-1),PSURF(IP+1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
C  SIDE 4
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR),PSURF(IP+1),
     .                       RHOSRF(IR-1),PSURF(IP+1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
              ENDIF
            ENDIF
C  LEFT HEMISPHERE
            IF (IR.LT.NR1ST) THEN
C  UPPER LEFT QUADRANT
              IF (IP.LT.NP2ND) THEN
C  SIDE 5
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR+1),PSURF(IP+1),
     .                       RHOSRF(IR),PSURF(IP+1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
C  SIDE 6
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR+1),PSURF(IP),
     .                       RHOSRF(IR+1),PSURF(IP+1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
              ENDIF
C  LOWER LEFT QUADRANT
              IF (IP.GT.1) THEN
C  SIDE 7
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR+1),PSURF(IP),
     .                       RHOSRF(IR+1),PSURF(IP-1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
C  SIDE 8
                CALL MUELAM (RHOSRF(IR),PSURF(IP),
     .                       YWERT(IR,IP),ZWERT(IR,IP),
     .                       RHOSRF(IR+1),PSURF(IP-1),
     .                       RHOSRF(IR),PSURF(IP-1),XL,XM)
                if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 101
              ENDIF
            ENDIF
            CYCLE IPLOOP1
101         CONTINUE
            ywert(IR,IP)=ywert(IR,IP)*xl*0.9
            zwert(IR,IP)=zwert(IR,IP)*xl*0.9
          ENDDO IPLOOP1
        ENDDO

        I=0
        DO IR=1,NR1ST
          DO IP=1,NP2ND
            I=I+1
            BETRAG=SQRT(YWERT(IR,IP)**2+ZWERT(IR,IP)**2)
            IF (BETRAG .GT. 1.E-5)
     .      WRITE(nraps2,'(I6,1P,5E12.4)')
     .           I,YWERT(IR,IP),zwert(IR,IP),0.,0.,0.
          enddo
        enddo

      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL
     .        .AND.LPTOR3(IBLD)) THEN
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
            IPLOOP2: DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
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
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
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
              CYCLE IPLOOP2
100           CONTINUE
              ywert(IR,IP)=ywert(IR,IP)*xl*0.9
              zwert(IR,IP)=zwert(IR,IP)*xl*0.9
            ENDDO IPLOOP2
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
              write (iunout,*) 'error in rpsvec: zuord'
              call exit_own(1)
            endif
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0,1),1) = J
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0,1),2) = I
50        CONTINUE
40      CONTINUE
        DO 61 I=1,NRKNOT
          IF (ZUORD(I,0,1).LT.1) THEN
            WRITE (iunout,*) 'ERROR IN RPSVEC: POINT ',I,' NOT IN MESH'
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
      ELSE
        WRITE (iunout,*) 'UNWRITTEN OPTION IN RPSVEC: PLOT ABANDONNED '
        WRITE (iunout,*) 'PRINTING OF RAPS FILE FAILED '
      ENDIF
C
      WRITE(NRAPS2,'(1X,A5,8X,A3,50(11X,I1))') '-9999',
     .           'FIN',0,0,0
      CLOSE (UNIT=NRAPS2)
C
      IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
      IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
      IF (ALLOCATED(ZWERT)) DEALLOCATE (ZWERT)
      IF (ALLOCATED(ZWERT1)) DEALLOCATE (ZWERT1)
      RETURN
      END
C ===== SOURCE: sccone.f
C
C
      SUBROUTINE SCCONE(X0,Y0,Z0,VX,VY,VZ,ALF,T1,T2,BX,BY,BZ,CX,CY,CZ,
     .                  DANG,A,I,XP,YP,JA,JE,IXS)
C
C  CALLED FROM CONE
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(*)
      REAL(DP), INTENT(OUT) :: XP(*),YP(*)
      REAL(DP), INTENT(IN) :: X0, Y0, Z0, VX, VY, VZ, ALF, T1, T2,
     .                      BX, BY, BZ, CX, CY, CZ, DANG
      INTEGER, INTENT(IN) :: I, JA, JE, IXS
      REAL(DP) :: EPS12, XK2, YK2, PXX1, PYY1, PZZ1, PXX2, PYY2, PZZ2,
     .          PX2, PY2, PZ2, XK1, YK1, PHI, DET, ALAM1, ALAM2, ROOT,
     .          XX, YY, ZZ, CC, PX21, PY21, PZ21, ALAMDA, AA, BB, XN,
     .          RAD1, RAD2, PX1, PY1, PZ1
      INTEGER :: J, IX
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
            WRITE (iunout,*) 'ERROR IN SUBR. SCCONE. SET ALAMDA=0.'
            WRITE (iunout,*) 'NO INTERSECTION FOUND WITH PLANE'
            ALAMDA=0.
            LERR=.TRUE.
          ELSE
            ALAMDA=(-A(1)-(A(2)*PXX1+A(3)*PYY1+A(4)*PZZ1))/XN
            IF (ALAMDA.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SCCONE. SET ALAMDA=0.'
              WRITE (iunout,*) 'NO INTERSECTION IN POSITIV DIRECTION'
              WRITE (iunout,*) 'WITH PLANE '
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
            WRITE (iunout,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
            WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (iunout,*) 'AA=BB=0.'
            GOTO 101
          ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
            ALAM1=-CC/BB
            ALAM2=-CC/BB
          ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
            DET=-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
              WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (iunout,*) 'DETERMINANT DET= ',DET
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
              WRITE (iunout,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
              WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (iunout,*) 'DETERMINANT DET= ',DET
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
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0.'
            WRITE (iunout,*) 'INTERSECTION IN WRONG DIRECTION'
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
C ===== SOURCE: schnitp.f


      SUBROUTINE SCHNITP(X1,Y1,X2,Y2,X3,Y3,X4,Y4,EX,EY)
C   INTERSECTION POINT E OF 2 STRAIGHT LINES G1 AND G2
C   G1 IS DEFINED BY POINTS 1 AND 2
C   G2 IS DEFINED BY POINTS 3 AND 4

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, Y1, X2, Y2, X3, Y3, X4, Y4
      REAL(DP), INTENT(OUT) :: EX, EY
      REAL(DP) :: MUE

      MUE = ((Y2-Y4)*(X3-X4)+(X4-X2)*(Y3-Y4))/
     .      ((X1-X2)*(Y3-Y4)-(Y1-Y2)*(X3-X4)+1.D-20)
      EX = X2 + MUE * (X1-X2)
      EY = Y2 + MUE * (Y1-Y2)
      END
C ===== SOURCE: secang.f
C
C
      SUBROUTINE SECANG (B0,B1,B2,RAD,ANG1,ANG2)

      USE PRECISION
      USE CCONA
      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: B0, B1, B2, RAD
      REAL(DP), INTENT(OUT) :: ANG1, ANG2
      REAL(DP) :: B0Q, B1Q, B2Q, X2, Y2, B0B1, RD, XNEN, X1, Y1
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
C ===== SOURCE: secqua.f
C
C
      SUBROUTINE SECQUA(XX,YY,ZZ,VX,VY,VZ,A,I,ALAM1,ALAM2,LERR)
C
C MORE GENERAL THAN SHNITT OR CSCONE, CALLED FROM PL3D ITSELF
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(*)
      REAL(DP), INTENT(IN) :: XX, YY, ZZ, VX, VY, VZ
      REAL(DP), INTENT(OUT) :: ALAM1, ALAM2
      INTEGER, INTENT(IN) :: I
      LOGICAL, INTENT(OUT) :: LERR
      REAL(DP) :: EPS12, XN, AA, DET, ROOT, BB, CC

      DATA EPS12 /1.E-12/

      LERR=.FALSE.
      IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
        XN=VX*A(2)+VY*A(3)+VZ*A(4)
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
        IF (ABS(XN).LT.EPS12) THEN
          WRITE (iunout,*) 'ERROR IN SUBR. SECQUA. SET ALAM1=ALAM2=0.'
          WRITE (iunout,*) 'NO INTERSECTION FOUND WITH PLANE'
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
          WRITE (iunout,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
          WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
          WRITE (iunout,*) 'AA=BB=0.'
        ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
          ALAM1=-CC/BB
          ALAM2=-CC/BB
        ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
          DET=-CC/AA
          IF (DET.LT.0.) THEN
            WRITE (iunout,*) 
     .        'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
            WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (iunout,*) 'DETERMINANT DET= ',DET
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
            WRITE (iunout,*) 
     .        'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
            WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (iunout,*) 'DETERMINANT DET= ',DET
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
C ===== SOURCE: shnitt.f
C
C
      SUBROUTINE SHNITT(P,PX,PY,PZ,VX,VY,VZ,A,I,XP,YP,JA,JE,IXS)
C
C  CALLED FROM ZYLIND
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: P(3,*), A(*)
      REAL(DP), INTENT(OUT) :: XP(*), YP(*)
      REAL(DP), INTENT(IN) :: PX, PY, PZ, VX, VY, VZ
      INTEGER, INTENT(IN) :: I, JA, JE, IXS

      REAL(DP) :: EPS12, XX, YY, ZZ, XN, ALAMDA, ALAM2, DET, ROOT,
     .            ALAM1, AA, BB, CC
      INTEGER :: J, IX
      LOGICAL :: LERR

      DATA EPS12 /1.E-12_DP/

      LERR=.FALSE.
      IX=IXS-1
C     WRITE (iunout,*) ' SHNITT  I = ',I
      IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
C       WRITE (iunout,*) ' VX,VY,VZ ',VX,VY,VZ
C       WRITE (iunout,*) ' A ',A(1),A(2),A(3),A(4)
        XN=VX*A(2)+VY*A(3)+VZ*A(4)
C       WRITE (iunout,*) ' XN ',XN
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
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT. SET ALAMDA=0.'
            WRITE (iunout,*) 'NO INTERSECTION FOUND WITH PLANE'
            ALAMDA=0.
            LERR=.TRUE.
          ELSE
            ALAMDA=(-A(1)-(A(2)*XX+A(3)*YY+A(4)*ZZ))/XN
            IF (ALAMDA.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SHNITT. SET ALAMDA=0.'
              WRITE (iunout,*) 'NO INTERSECTION IN POSITIV DIRECTION'
              WRITE (iunout,*) 'WITH PLANE '
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
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
            WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (iunout,*) 'AA=BB=0.'
            GOTO 201
          ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
            ALAM1=-CC/BB
            ALAM2=-CC/BB
          ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
            DET=-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
              WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (iunout,*) 'DETERMINANT DET= ',DET
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
              WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
              WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (iunout,*) 'DETERMINANT DET= ',DET
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
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0.'
            WRITE (iunout,*) 'INTERSECTION IN WRONG DIRECTION'
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
C ===== SOURCE: stcoor.f
C
C
      SUBROUTINE STCOOR (X,Y,IFL)
c  instor: zaehler, instor-ter aufruf
c  ifl: flag:=0, first point, .ne.0: else, npl2d(instor)=ifl
c  x,y: co-ordinaten des punktes, stored on xpl2d(instor),ypl2d(instor)
C
c  alle teilstuecke stehen hintereinander auf xpl2d,ypl2d
c  zum entwirren:NUMSUR(inums,..) array
c  inums: zaehler fuer individuelle zu plottende kurvenstuecke
c  jedes bekommt die eigene nummer und arrow
c  daher: ggfls: eine nummer mehrmals im plot vorhanden
c  NUMSUR : nummer des flachestuecks
C
      USE PRECISION
      USE PARMMOD
      USE CRECH

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Y
      INTEGER, INTENT(IN) :: IFL
      TYPE(PPOINT), POINTER :: POINT
      INTEGER, SAVE :: IFIRST=0
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        INUMS=0
        NULLIFY(FIRST_POINT)
        NULLIFY(LAST_POINT)
      ENDIF
C
      ALLOCATE(POINT)
      POINT%XPL2D = X
      POINT%YPL2D = Y
      POINT%NPL2D = IFL
      POINT%NUMSUR = INOSF
      NULLIFY(POINT%NXTPNT)

      IF (ASSOCIATED(LAST_POINT)) THEN
        LAST_POINT%NXTPNT => POINT
        LAST_POINT => POINT
      ELSE
        FIRST_POINT => POINT
        LAST_POINT => FIRST_POINT
      END IF

      INUMS = INUMS + 1
      INSTOR = INSTOR + 1
C
      RETURN

C     the following ENTRY is for reinitialization of EIRENE
     
      ENTRY STCOOR_REINIT
      IFIRST = 0
      return

      END
C ===== SOURCE: torloc.f
C
C
      SUBROUTINE TORLOC(W3,RM,X1,Z1)
C
C TRANSFORM FROM TORUS- TO LOCAL SYSTEM AT PHI=W3
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: W3, RM
      REAL(DP), INTENT(INOUT) :: X1, Z1
      REAL(DP), SAVE :: WO=0.D0, S=0.D0, C=1.D0
      REAL(DP) :: XS

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
C ===== SOURCE: tstchm.f
C
C
      SUBROUTINE TSTCHM(I1,XTN,YTN,XT,YT,IINDEX,TESTN,XMI,XMA,YMI,YMA,
     .                     XT2,YT2)
C
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
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: XTN, YTN, XT, YT, XT2, YT2
      REAL(DP), INTENT(IN) :: XMI, XMA, YMI, YMA
      REAL(DP), INTENT(OUT) :: TESTN
      INTEGER, INTENT(IN) :: I1
      INTEGER, INTENT(OUT) :: IINDEX
      REAL(DP) :: XS(2),YS(2)
      REAL(DP) :: TBDXN, TBDYN, TEST, XTO, YTO, TBDXO, TBDYO, T,
     .          TIME, TESTX, TESTY, TESTO, TSEC, ZXO, ZYO, ZXN,
     .          ZYN, ZBX1, ZBY1, ZBX2, ZBY2
      INTEGER :: IT

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
            WRITE (iunout,*) ' NONSENSE IN TSTCHM '
            WRITE (iunout,*) ' XTO,YTO ',XTO,YTO
            WRITE (iunout,*) ' XTN,YTN ',XTN,YTN
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
C ===== SOURCE: veclne.f
C
C
      SUBROUTINE VECLNE (AORIG,BORIG,IBLD,ICURV,
     .                   IXXI,IXXE,IYYI,IYYE,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A VECTOR FIELD PLOT
C  MODIFIED 160600: LEVGEO=3 OPTION ADDED
C                   SCALING OF ARROWS AUTOMATICALLY
C  ARGUMENTS: XX,YY OUT, IXX,IYY --> IXXI,IXXE,IYYI,IYYE
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTRIG
      IMPLICIT NONE
C
C
      REAL(DP), INTENT(IN) :: AORIG(*), BORIG(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV
      INTEGER, INTENT(INOUT) :: IXXI, IXXE, IYYI, IYYE
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: HEAD, RUNID, TXHEAD
      CHARACTER(72), INTENT(OUT) :: TEXT1
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: SCLFCX, SCLFCY, VMIN, VMAX, DX, DY, CM, FAK, BRFL,
     .          D, DIAMETER, XM, YM, VX, VY, VABS, XMIN, XMAX, YMIN,
     .          YMAX, PLFL
      REAL(SP) XY(800)
      REAL(SP) YH
      INTEGER :: ITR, IRAD, I, IR, IPART, NP1, NP2, IP
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
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NT3RD)
        IF (IYYE.LE.0) IYYE=NT3RD
        XMIN = RHOSRF(IXXI)
        XMAX = RHOSRF(IXXE)
        YMIN = ZSURF(IYYI)
        YMAX = ZSURF(IYYE)
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NP2ND)
        IF (IYYE.LE.0) IYYE=NP2ND
        XMIN = RHOSRF(IXXI)
        XMAX = RHOSRF(IXXE)
        YMIN = PSURF(IYYI)
        YMAX = PSURF(IYYE)
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
C  SUFFICIENT TO SEARCH ON OUTERMOST RADIAL SURFACE (BECAUSE: CONVEX)
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NP2ND)
        IF (IYYE.LE.0) IYYE=NP2ND
        DO 5 IP = IYYI,IYYE
          XMIN = MIN(XMIN,XPOL(IXXE,IP))
          XMAX = MAX(XMAX,XPOL(IXXE,IP))
          YMIN = MIN(YMIN,YPOL(IXXE,IP))
          YMAX = MAX(YMAX,YPOL(IXXE,IP))
5       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
C  SEARCH ON WHOLE MESH
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NP2ND)
        IF (IYYE.LE.0) IYYE=NP2ND
        DO 1 IR=IXXI,IXXE
          NP1=MAX(IYYI,NPOINT(1,1))
          NP2=MIN(IYYE,NPOINT(2,NPPLG))
          XMIN=MIN(XMIN,XPOL(IR,NP1),XPOL(IR,NP2))
          YMIN=MIN(YMIN,YPOL(IR,NP1),YPOL(IR,NP2))
          XMAX=MAX(XMAX,XPOL(IR,NP1),XPOL(IR,NP2))
          YMAX=MAX(YMAX,YPOL(IR,NP1),YPOL(IR,NP2))
1       CONTINUE
C
        DO 4 IPART=1,NPPLG
          DO 2 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            IF (IP.LT.IYYI.OR.IP.GT.IYYE) GOTO 2
            XMIN=MIN(XMIN,XPOL(IXXI,IP))
            YMIN=MIN(YMIN,YPOL(IXXI,IP))
            XMAX=MAX(XMAX,XPOL(IXXI,IP))
            YMAX=MAX(YMAX,YPOL(IXXI,IP))
2         CONTINUE
          DO 3 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            IF (IP.LT.IYYI.OR.IP.GT.IYYE) GOTO 3
            XMIN=MIN(XMIN,XPOL(IXXE,IP))
            YMIN=MIN(YMIN,YPOL(IXXE,IP))
            XMAX=MAX(XMAX,XPOL(IXXE,IP))
            YMAX=MAX(YMAX,YPOL(IXXE,IP))
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
      CALL GRNXTB (1,'VECLNE.F')
      CALL GRSCLC (10.,4.,REAL(10.+DX*FAK,KIND(1.E0)),
     .             REAL(4.+DY*FAK,KIND(1.E0)))
      CALL GRSCLV (REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)),
     .             REAL(XMAX,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
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
        CALL GRJMP(REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)))
        CALL GRDRW(REAL(XMIN,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
        CALL GRDRW(REAL(XMAX,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
        CALL GRDRW(REAL(XMAX,KIND(1.E0)),REAL(YMIN,KIND(1.E0)))
        CALL GRDRW(REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)))
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        CALL GRJMP(REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)))
        CALL GRDRW(REAL(XMIN,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
        CALL GRDRW(REAL(XMAX,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
        CALL GRDRW(REAL(XMAX,KIND(1.E0)),REAL(YMIN,KIND(1.E0)))
        CALL GRDRW(REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)))
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
        DO 7 IR=IXXI,IXXE,IXXE-IXXI
          CALL GRJMP(REAL(XPOL(IR,IYYI),KIND(1.E0)),
     .               REAL(YPOL(IR,IYYI),KIND(1.E0)))
          DO IP = IYYI+1,IYYE
            CALL GRDRW(REAL(XPOL(IR,IP),KIND(1.E0)),
     .                 REAL(YPOL(IR,IP),KIND(1.E0)))
          END DO
7       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
        NP1=MAX(IYYI,NPOINT(1,1))
        NP2=MIN(IYYE,NPOINT(2,NPPLG))
        CALL GRJMP(REAL(XPOL(IXXI,NP1),KIND(1.E0)),
     .             REAL(YPOL(IXXI,NP1),KIND(1.E0)))
          DO 10 IR=IXXI+1,IXXE
10          CALL GRDRW (REAL(XPOL(IR,NP1),KIND(1.E0)),
     .                  REAL(YPOL(IR,NP1),KIND(1.E0)))
C
        CALL GRJMP (REAL(XPOL(IXXI,NP2),KIND(1.E0)),
     .              REAL(YPOL(IXXI,NP2),KIND(1.E0)))
        DO 11 IR=IXXI+1,IXXE
11        CALL GRDRW (REAL(XPOL(IR,NP2),KIND(1.E0)),
     .                REAL(YPOL(IR,NP2),KIND(1.E0)))
        DO 15 I=1,NPPLG
          NP1=MAX(IYYI,NPOINT(1,I))
          NP2=MIN(IYYE,NPOINT(2,I))
          CALL GRJMP (REAL(XPOL(IXXI,NP1),KIND(1.E0)),
     .                REAL(YPOL(IXXI,NP1),KIND(1.E0)))
          DO 12 IP=NP1,NP2
12          CALL GRDRW (REAL(XPOL(IXXI,IP),KIND(1.E0)),
     .                  REAL(YPOL(IXXI,IP),KIND(1.E0)))
          CALL GRJMP (REAL(XPOL(IXXE,NPOINT(1,I)),KIND(1.E0)),
     .                REAL(YPOL(IXXE,NPOINT(1,I)),KIND(1.E0)))
          DO 13 IP=NP1,NP2
13          CALL GRDRW (REAL(XPOL(IXXE,IP),KIND(1.E0)),
     .                  REAL(YPOL(IXXE,IP),KIND(1.E0)))
15      CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO ITR=1,NTRII
          IF (NCHBAR(1,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(1,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(2,ITR)),KIND(1.E0)))
          ENDIF
          IF (NCHBAR(2,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(3,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(3,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(2,ITR)),KIND(1.E0)))
          ENDIF
          IF (NCHBAR(3,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(1,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(3,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(3,ITR)),KIND(1.E0)))
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
          DO 900 IR=IXXI,IXXE-1
          DO 900 IP=IYYI,IYYE-1
            IRAD = IR + (IP-1)*NR1ST
            VABS = SQRT(AORIG(IRAD)**2 + BORIG(IRAD)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
900       CONTINUE
C 2. SCALING
          DO 1100 IR=IXXI,IXXE-1
          DO 1100 IP=IXXI,IYYE-1
            IRAD = IR + (IP-1)*NR1ST
            VX = (AORIG(IRAD) / VMAX)
            VY = (BORIG(IRAD) / VMAX)
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOM(IRAD)-VX/2
            YM=YCOM(IRAD)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(REAL(xm,KIND(1.E0)),REAL(ym,KIND(1.E0)),
     .                  REAL(xm+vx,KIND(1.E0)),REAL(ym+vy,KIND(1.E0)),
     .                  REAL(PLFL,KIND(1.E0)),REAL(BRFL,KIND(1.E0)),1)
1100      CONTINUE
C
      ELSEIF ((LEVGEO.EQ.2.AND.NLPOL).OR.
     .         LEVGEO.EQ.3) THEN
C
C 1. SEARCH FOR MINIMA AND MAXIMA
          VMIN= 1.D60
          VMAX=-1.D60
          DIAMETER=1.D60
          DO 901 IR=IXXI,IXXE-1
          DO 901 IP=IXXI,IYYE-1
            IRAD = IR + (IP-1)*NR1ST
            D=SQRT((XPOL(IR,IP)-XPOL(IR+1,IP+1))**2+
     .             (YPOL(IR,IP)-YPOL(IR+1,IP+1))**2)
            DIAMETER=MIN(DIAMETER,D)
            VABS = SQRT(AORIG(IRAD)**2 + BORIG(IRAD)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
C 2. SCALING
            VX =  AORIG(IRAD) / VABS * D * 0.5
            VY =  BORIG(IRAD) / VABS * D * 0.5
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOM(IRAD)-VX/2
            YM=YCOM(IRAD)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(REAL(xm,KIND(1.E0)),REAL(ym,KIND(1.E0)),
     .                  REAL(xm+vx,KIND(1.E0)),REAL(ym+vy,KIND(1.E0)),
     .                  REAL(PLFL,KIND(1.E0)),REAL(BRFL,KIND(1.E0)),1)
901       CONTINUE
C
      ELSEIF (LEVGEO.EQ.4) THEN
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
            XM=XCOM(IR)-VX/2
            YM=YCOM(IR)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(REAL(xm,KIND(1.E0)),REAL(ym,KIND(1.E0)),
     .                  REAL(xm+vx,KIND(1.E0)),REAL(ym+vy,KIND(1.E0)),
     .                  REAL(PLFL,KIND(1.E0)),REAL(BRFL,KIND(1.E0)),1)
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
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,RUNID)
      YH=26.75
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,HEAD)
      YH=26.00
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,REAL(YH-0.5,KIND(1.E0)),10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,REAL(YH-1.,KIND(1.E0)),10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,REAL(YH-2.,KIND(1.E0)),10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') VMAX
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') VMIN
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)
C
      RETURN
      END
C ===== SOURCE: xyplot.f
C
C
      SUBROUTINE XYPLOT (XY,NC)

      USE PRECISION

      IMPLICIT NONE
C
      REAL(SP), INTENT(INOUT) :: XY(*)
      INTEGER, INTENT(IN) :: NC
      REAL(DP) :: DXY, HILF
      INTEGER :: I, IANF, NBEG, J
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
      CALL GRJMP (REAL(XY(NBEG),KIND(1.E0)),REAL(XY(NBEG+1),KIND(1.E0)))
      DO 5 I=NBEG+2,IANF,4
5       CALL GRDRW (REAL(XY(I),KIND(1.E0)),REAL(XY(I+1),KIND(1.E0)))
C
      IF (IANF.LT.NC-1) THEN
        NBEG=IANF+2
        IANF=IANF+4
        GOTO 1
      ENDIF
C
      RETURN
      END
C ===== SOURCE: zylind.f
C
C
      SUBROUTINE ZYLIND (X0,Y0,Z0,VX,VY,VZ,T1,T2,RAD,NK,NP,NA,IO,NF,NUM,
     .                   ILEFT,AL,IRIGHT,AR,PHIAN,PHIEN)
C
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
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X0, Y0, Z0, VX, VY, VZ, T1, T2, RAD,
     .                      PHIAN, PHIEN
      REAL(DP), INTENT(IN) :: AL(10), AR(10)
      INTEGER, INTENT(IN) :: NK, NP, NA, IO, NUM, ILEFT, IRIGHT
      LOGICAL, INTENT(IN) :: NF
      REAL(DP), ALLOCATABLE :: P(:,:), XP(:), YP(:)
      REAL(DP) :: PHI, PHAN, XK, YK, PXS, PYS, CX, CY, CZ, BA, BB, BC,
     .          DANG, PZS, PSS, PX, T, PY, PYY, PXX, PZ, AX, AY, AZ,
     .          PHID, BX, BY, BZ, PI, DT, PZZ
      REAL(SP), ALLOCATABLE :: XPS(:), YPS(:)
      INTEGER :: J, IA, JJ, IE, I, NAPK
C
      NAPK = MAX(NA,NP,NK) + 1
      ALLOCATE (P(3,NAPK))
      ALLOCATE (XP(NAPK))
      ALLOCATE (YP(NAPK))
      ALLOCATE (XPS(NAPK))
      ALLOCATE (YPS(NAPK))

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
        WRITE (iunout,*) 'FEHLER IN DER EINGABE VON VX,VY,VZ ',
     .    NUM,VX,VY,VZ
        CALL EXIT_OWN(1)
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

      DEALLOCATE (P)
      DEALLOCATE (XP)
      DEALLOCATE (YP)
      DEALLOCATE (XPS)
      DEALLOCATE (YPS)

      RETURN
      END
C ===== SOURCE: zylnd2.f
C
C
      SUBROUTINE ZYLND2(X0,Y0,Z0,VX,VY,VZ,TAR,RAD,NK,NP,NA,IO,NF,NUM,
     .                   ILEFT,AL,IRIGHT,AR,PHIAN,PHIEN,NDP,IPART)
C
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
      USE PRECISION
      USE CCONA
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NK, NP, NA, IO, NUM, ILEFT, IRIGHT, NDP
      REAL(DP), INTENT(IN) :: AL(*), AR(*), PHIAN(NDP,*), PHIEN(NDP,*),
     .                      TAR(NK)
      REAL(DP), INTENT(IN) :: X0, Y0, Z0, VX, VY, VZ, RAD
      INTEGER, INTENT(IN) :: IPART(NK)
      LOGICAL, INTENT(IN) :: NF

      REAL(DP) :: XK, YK, PHI, PI, PXS, PYS, PZS, PHAN, BA, BB, BC,
     .          T, PX, PHIDEG, PY, PZ, PXX, PYY, PZZ, CZ, BX, BY, BZ,
     .          AX, AY, AZ, CX, CY, PHID, DANG
      INTEGER :: J, I, IP, IA, IE, JP, JJ, ILFT, IRGHT, NAPK

      REAL(DP), ALLOCATABLE :: P(:,:), XP(:), YP(:)
      REAL(SP), ALLOCATABLE :: XPS(:), YPS(:)
C
      NAPK = MAX(NA,NP,NK) + 1
      ALLOCATE (P(3,NAPK))
      ALLOCATE (XP(NAPK))
      ALLOCATE (YP(NAPK))
      ALLOCATE (XPS(NAPK))
      ALLOCATE (YPS(NAPK))
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
        WRITE (iunout,*) 'FEHLER IN DER EINGABE VON VX,VY,VZ ',
     .    NUM,VX,VY,VZ
        CALL EXIT_OWN(1)
      ENDIF
      WRITE (iunout,*) ' PHID ',PHID*RADDEG
      WRITE (iunout,*) ' VX,VY,VZ ',VX,VY,VZ
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
          WRITE (iunout,*) ' LEFT END OF ZYLINDER '
          WRITE (iunout,*) (AL(ILFT),ILFT=1,ILEFT)
          CALL SHNITT(P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,1,NA+1,1)
          ELSEIF (I.EQ.NK.AND.IRIGHT.NE.0) THEN
          WRITE (iunout,*) ' RIGHT END OF ZYLINDER '
          WRITE (iunout,*) (AR(IRGHT),IRGHT=1,IRIGHT)
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

      DEALLOCATE (P)
      DEALLOCATE (XP)
      DEALLOCATE (YP)
      DEALLOCATE (XPS)
      DEALLOCATE (YPS)

      RETURN
      END
C ===== SOURCE: zylpln.f
C
C
      SUBROUTINE ZYLPLN (ZX0,ZY0,ZZ0,ZVX,ZVY,ZVZ,RZYL,JS,NZAD,NINNE,NIN)

      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CCONA
      USE CLGIN
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: ZX0, ZY0, ZZ0, ZVX, ZVY, ZVZ, RZYL
      INTEGER, INTENT(IN) :: JS, NZAD, NINNE, NIN
C
      REAL(DP) :: AFF(3,3),AFFI(3,3),PHIAN(20,9),PHIEN(20,9),ANGLE(19),
     .          TAR(20),AL(4),AR(4)
      REAL(DP) :: B0, B1, B2, PHI, PHIH, HILF, TMND, HILFD, X, Y, Z,
     .          PHID, BETC, C1, COSD, SIND, HYP, DLT, BPLN, CPHI, E, T,
     .          TMXD, TTEST, TMIN, EC, V1, C2, C3, TMAX, V2, V3,
     .          TEST
      INTEGER :: IPAR(20)
      INTEGER :: ISORT, IPART, IMIN, IANG, IZ, IE, JP, IMAX, ITEST,
     .           I, IHILF

C
      PHID=180.-ACOS(1./SQRT(2.))*RADDEG
C     PHID=0.
      ITEST=0
      WRITE (iunout,*) ' ZYLPLN CALLED FOR JS = ',JS
      WRITE (iunout,*) ' ZX0,ZY0,ZZ0 ',ZX0,ZY0,ZZ0
      WRITE (iunout,*) ' ZVX,ZVY,ZVZ ',ZVX,ZVY,ZVZ
C  WINKEL ZWISCHEN (ZVX,ZVY,ZVZ) UND (0,0,1)
      COSD=ZVZ/SQRT(ZVX*ZVX+ZVY*ZVY+ZVZ*ZVZ)
      SIND=SQRT(1.-COSD*COSD)
      WRITE (iunout,*) ' COSD,SIND ',COSD,SIND
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
      WRITE (iunout,*) ' V1,V2,V3 ',V1,V2,V3
      IF (ABS(V1)+ABS(V2).GT.EPS10.OR.V3.LT.0.) THEN
        IF (ITEST.EQ.0) THEN
          ITEST=ITEST+1
          C1=-C1
          C2=-C2
          C3=-C3
          GOTO 4711
        ELSE
          WRITE (iunout,*) ' NONSENSE IN ZYLPLN   JS = ',JS
          WRITE (iunout,*) ' PLOT OF THIS SURFACE ABANDONNED '
          RETURN
        ENDIF
      ENDIF
C
      WRITE (iunout,*) ' VOR ROTADD ALIMS,XLIMS,YLIMS,ZLIMS '
      WRITE (iunout,'(1X,1P,4E12.4)') 
     .      (ALIMS(I,JS),XLIMS(I,JS),YLIMS(I,JS),
     .                   ZLIMS(I,JS),I=1,ILIN(JS))
      CALL XSHADD (-ZX0,JS,JS)
      CALL YSHADD (-ZY0,JS,JS)
      CALL ZSHADD (-ZZ0,JS,JS)
      CALL ROTADD (AFF,AFFI,JS,JS)
      WRITE (iunout,*) ' NACH ROTADD ALIMS,XLIMS,YLIMS,ZLIMS '
      WRITE (iunout,'(1X,1P,4E12.4)') 
     .             (ALIMS(I,JS),XLIMS(I,JS),YLIMS(I,JS),
     .                          ZLIMS(I,JS),I=1,ILIN(JS))
C
C  SUCHE T-INTERVALL
C
      TMAX=1.D60
      TMIN=-1.D60
      DO 100 I=1,ILIN(JS)
        IF (ABS(XLIMS(I,JS)).LT.EPS10.AND.ABS(YLIMS(I,JS)).LT.EPS10.AND.
     .      ABS(ZLIMS(I,JS)).LT.EPS10) THEN
          WRITE (iunout,*) ' NONSENSE IN ZYLPLN WITH SURFACE ',JS
          WRITE (iunout,*) ' PLOT OF SURFACE ABORTED '
          CALL ROTADD (AFFI,AFF,JS,JS)
          CALL XSHADD (ZX0,JS,JS)
          CALL YSHADD (ZY0,JS,JS)
          CALL ZSHADD (ZZ0,JS,JS)
          RETURN
        ENDIF
C
C
C  WINKEL DER EBENE MIT DER Z=0 EBENE
        BPLN=SQRT(XLIMS(I,JS)**2+YLIMS(I,JS)**2+ZLIMS(I,JS)**2)
        CPHI=ZLIMS(I,JS)/BPLN
C  EBENE PARALLEL ZUM ZYLINDER?
        IF (ABS(CPHI).LT.EPS10) GOTO 100
        HYP=RZYL/CPHI
        DLT=SQRT(HYP*HYP-RZYL*RZYL)
C  SCHNITT DER EBENE MIT (0,0,1)
        T=-ALIMS(I,JS)/ZLIMS(I,JS)
        TTEST=T-1.
C  FESTSTELLEN, OB (0,0,TTEST) EIN GUELTIGER PUNKT IST
        E=ZLIMS(I,JS)*TTEST+ALIMS(I,JS)
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
      WRITE (iunout,*) ' TMIN,TMAX ',TMIN,TMAX
      WRITE (iunout,*) ' TMND,TMXD ',TMND,TMXD
      WRITE (iunout,*) ' IMIN,IMAX ',IMIN,IMAX
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
          B0=ALIMS(I,JS)+ZLIMS(I,JS)*T
          B1=XLIMS(I,JS)
          B2=YLIMS(I,JS)
          IANG=IANG+2
          CALL SECANG (B0,B1,B2,RZYL,ANGLE(IANG-1),ANGLE(IANG))
110     CONTINUE
        WRITE (iunout,*) ' ANGLE VOR SORT ',(ANGLE(I),I=1,IANG)
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
        WRITE (iunout,*) ' ANGLE NACH SORT ',(ANGLE(I),I=1,IANG)
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
C         WRITE (iunout,*) ' ANGLE(I),(I+1) ',ANGLE(I),ANGLE(I+1)
C         WRITE (iunout,*) ' PHI,X,Y,Z, ',PHI,X,Y,Z
          DO 125 IE=1,ILIN(JS)
            TEST=ALIMS(IE,JS)+XLIMS(IE,JS)*X+YLIMS(IE,JS)*Y+
     .           ZLIMS(IE,JS)*Z
C           WRITE (iunout,*) ' IE,TEST ',IE,TEST
            IF (TEST.GT.0.) GOTO 130
125       CONTINUE
          IPART=IPART+1
          PHIAN(IZ,IPART)=ANGLE(I)
          PHIEN(IZ,IPART)=ANGLE(I+1)
130     CONTINUE
        TAR(IZ)=T
        IPAR(IZ)=IPART
C
        WRITE (iunout,*) ' IZ,IPAR,TAR ',IZ,IPAR(IZ),TAR(IZ)
        WRITE (iunout,'(1X,1P,2E12.4)')
     .    (PHIAN(IZ,JP),PHIEN(IZ,JP),JP=1,IPART)
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
      AL(1)=ALIMS(IMIN,JS)
      AL(2)=XLIMS(IMIN,JS)
      AL(3)=YLIMS(IMIN,JS)
      AL(4)=ZLIMS(IMIN,JS)
C
      AR(1)=ALIMS(IMAX,JS)
      AR(2)=XLIMS(IMAX,JS)
      AR(3)=YLIMS(IMAX,JS)
      AR(4)=ZLIMS(IMAX,JS)
C
C
      CALL ZYLND2 (ZX0,ZY0,ZZ0,ZVX,ZVY,ZVZ,TAR,RZYL,NZAD+2,NINNE,NIN,
     .     ILCOL(JS),IGFIL(JS).NE.0,JS,4,AL,4,AR,PHIAN,PHIEN,20,IPAR)
C
      RETURN
      END
