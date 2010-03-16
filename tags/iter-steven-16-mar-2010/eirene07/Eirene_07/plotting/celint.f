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
