C
      SUBROUTINE VOLUME (IND)

C  CALCULATE VOLUME-ELEMENTS FOR VOLUME AVERAGED TALLIES
C  THE CELL VOLUMES VOL MUST BE THOSE SEEN BY THE TESTPARTICLES
C  I.E. NOT NECESSARLY THE TRUE ONES.
C  ONE COMMON FACTOR (LENGTH OF THE CELL IN IGNORABLE
C  DIMENSION) ACTS LIKE A SCALING FACTOR FOR THESE TALLIES.
C
C  IN CASE OF THE NLTRA OPTION: IF (NLTOR):
C                               VOL = TAN(ALPHA)*XCOM*AREA*2.
C                               VOL = VOLUME OF ONE OF THE NTTRAM
C                                     CYLINDRICAL SEGMENTS
C                               AREA= AREA OF THE CELL
C                               ALPHA = 0.5*(2*PI / NTTRAM)
C                               XCOM= X - CENTER OF MASS OF AREA
C                               (I.E., NTTRAM=PI/ALPHA)
C                               IF (.NOT.NLTOR):
C                               VOL= NTTRAM TIMES THE VOLUME GIVEN
C                                    GIVEN ABOVE, IE:
C                               VOL= TAN(ALPHA)/ALPHA*2*PI*XCOM*AREA
C  IN CASE OF THE NLTRT OPTION: VOL = 2*PI*XCOM*AREA
C  (REGARDLESS OF NLTRZ,...)    VOL = VOLUME OF THE CELL, NO TOROIDAL
C                                     RESOLUTION
C                               AREA= AREA OF THE CELL
C                               XCOM= X - CENTER OF MASS OF AREA
C  NOTE: NLTRT=TRUE INTRODUCES INCONCISTENCY, SINCE PARTICLES
C        MAY SEE CYLINDER, BUT VOLUME IS COMPUTED FOR TORUS
C
C  NOTE: XCOM IS CENTER OF MASS IN TORUS SYSTEM, I.E.,
C        THE LARGE RADIUS RMTOR MUST BE ADDED TO THE
C        XCOM EVALUATED IN LOCAL SYSTEMS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTETRA
      USE CTRIG

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IND
      REAL(DP), ALLOCATABLE, SAVE :: AREAP(:,:)
      REAL(DP) :: AREA1(0:N1ST)
      REAL(DP) :: PC1(3), PC2(3), PC3(3), PC4(3)
      REAL(DP) :: AREAR, VOLSR, CAL_VOL, TWOTHIRD, VSAVE, FAC2, FAC3,
     .          PI2AT, AELL, DONE, DNULL, SY, X1, X2, Y1, XNULL, SX,
     .          ARTRIA, AR, XC, Y2, X3, Y3, X4, Y4
      INTEGER :: NCELL1, I, KP, IC4, IC, ITET, IC1, IC2, IC3, NCELLK,
     .           IT, NCELLJ, NCELL, IPP, I1ST, IRAD, IR, IFLAG, J,
     .           IP, JP, IN, K, IRP
!pb      SAVE

C     IND=1: 1-ST GRID, RAD. RESOLUTION
C     IND=2: 2-ND GRID, POL. RESOLUTION
C     IND=3: 3-RD GRID, TOR. RESOLUTION
C     IND=4: ADDITIONAL CELL REGION
C
      GOTO(100,200,300,400),IND
C
100   CONTINUE
C
      ALLOCATE (AREAP(N1STS,N2NDPLG))

      DO 101 IRAD=1,NRAD
        VOL(IRAD)=0.
101   CONTINUE
      AREA1(0)=0.D0
      DO 102 I1ST=1,N1ST
        AREA1(I1ST)=0.D0
102   CONTINUE
      DO 103 IRAD=1,NRAD
        AREA(IRAD)=0.D0
103   CONTINUE
C
      IF (LEVGEO.EQ.1) THEN
C
C 1D SLAB-MODEL, DY = YDF, DZ = ZDF
C
        DO 110 IR=1,NR1STM
          AREA1(IR)=(RSURF(IR+1)-RSURF(IR))*YDF
          SX=(RSURF(IR+1)+RSURF(IR))*0.5
          IF (NLTRZ) THEN
            VOL(IR)=AREA1(IR)*ZDF
          ELSEIF (NLTRA) THEN
            VOL(IR)=AREA1(IR)*(SX+RMTOR)*TANAL/ALPHA*PI2A
          ELSE
            WRITE (6,*) 'INVALID OPTION IN SUBR. VOLUME, EXIT CALLED'
            CALL EXIT_OWN(1)
          ENDIF
110     CONTINUE
        GOTO 190
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
C 1D GRID OF CIRCLES, ELLIPSES OR TRIANGULAR FLUXSURFACES
C
        XNULL=0.
        IFLAG=1
        DNULL=0.
        DONE=1.
        CALL ARELLP(EP1(1),DNULL,ELL(1),DONE,TRI(1),DNULL,
     .              RSURF(1),DNULL,PI2A,XNULL,IFLAG,
     .              AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
        AREA1(0)=AELL
        DO 121 IR=1,NR1STM
C  AREA, CENTER OF GRAVITY
          IRP=IR+1
          CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                TRI(IRP),TRI(IR),
     .                RSURF(IRP),RSURF(IR),PI2A,XNULL,IFLAG,
     .                AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
          AREA1(IR)=AELL
          IF (NLTRT) THEN
            VOL(IR)=AREA1(IR)*(SX+RMTOR)*PI2A
          ELSEIF (NLTRZ) THEN
            VOL(IR)=AREA1(IR)*ZDF
          ELSEIF (NLTRA) THEN
            VOL(IR)=AREA1(IR)*(SX+RMTOR)*TANAL/ALPHA*PI2A
          ENDIF
121     CONTINUE
        GOTO 190
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
C   1D GRID OF POLYGONS
C
        DO 139 K=1,NPPLG
          DO 139 J=NPOINT(1,K),NPOINT(2,K)-1
          AR=ARTRIA(0._DP,0._DP,XPOL(1,J),YPOL(1,J),
     .                          XPOL(1,J+1),YPOL(1,J+1))
          AREA1(0)=AREA1(0)+AR
139     CONTINUE
        AREA1(0)=ABS(AREA1(0))
C
        CALL ARPOLY(XPOL,YPOL,NRPLG,N1STS,1,NR1ST,AREAP,XCOM,YCOM)
        DO IR=1,NR1STM
          DO IP=1,NRPLG-1
            IN = IR + (IP-1)*NR1ST
            AREA(IN) = AREAP(IR,IP)
          END DO
        END DO
C
        DO 131 IR=1,NR1STM
          DO 132 K=1,NPPLG
            DO 132 J=NPOINT(1,K),NPOINT(2,K)-1
              AREA1(IR)=AREA1(IR)+AREAP(IR,J)
132       CONTINUE
131     CONTINUE
C
        IF (NLTRT) THEN
C   PARTICLES SEE A TORUS
          DO 134 IR=1,NR1STM
            XC=0.
            AR=0.
            DO 135 JP=1,NPPLG
              DO 135 J=NPOINT(1,JP),NPOINT(2,JP)-1
                IN = IR + (J-1)*NR1ST
                XC=XC+AREAP(IR,J)*XCOM(IN)
                AR=AR+AREAP(IR,J)
135         CONTINUE
            XC=XC/(AR+EPS60)
            VOL(IR)=AREA1(IR)*(XC+RMTOR)*PI2A
134       CONTINUE
        ELSEIF (NLTRZ) THEN
C   PARTICLES SEE A CYLINDER OF LENGTH DZ = ZDF
          DO 133 IR=1,NR1STM
            VOL(IR)=AREA1(IR)*ZDF
133       CONTINUE
C   PARTICLES SEE A TORUS APPROXIMATED BY NTTRAM STRAIGHT CYLINDERS
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 136 IR=1,NR1STM
            XC=0.
            AR=0.
            DO 137 JP=1,NPPLG
              DO 137 J=NPOINT(1,JP),NPOINT(2,JP)-1
                IN = IR + (J-1)*NR1ST
                XC=XC+AREAP(IR,J)*XCOM(IN)
                AR=AR+AREAP(IR,J)
137         CONTINUE
            XC=XC/(AR+EPS60)
            VOL(IR)=AREA1(IR)*(XC+RMTOR)*PI2AT
136       CONTINUE
        ENDIF
        GOTO 190
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  GRID DEFINED BY FINITE ELEMENTS
C
        DO 150 IR=1,NTRII
          XCOM(IR) = (XTRIAN(NECKE(1,IR))+XTRIAN(NECKE(2,IR))+
     .                  XTRIAN(NECKE(3,IR)))/3.
          YCOM(IR) = (YTRIAN(NECKE(1,IR))+YTRIAN(NECKE(2,IR))+
     .                  YTRIAN(NECKE(3,IR)))/3.
          AR=0.5*(XTRIAN(NECKE(2,IR))*(YTRIAN(NECKE(3,IR))
     >           -YTRIAN(NECKE(1,IR)))+XTRIAN(NECKE(3,IR))*
     >           (YTRIAN(NECKE(1,IR))
     >           -YTRIAN(NECKE(2,IR)))+XTRIAN(NECKE(1,IR))*
     >           (YTRIAN(NECKE(2,IR))-YTRIAN(NECKE(3,IR))))
          AREA(IR) = AR
150     CONTINUE
C
C
C   PARTICLES SEE A TORUS
        IF (NLTRT) THEN
          DO 154 IR=1,NTRII
            VOL(IR)=AREA(IR)*(XCOM(IR)+RMTOR)*PI2A
154       CONTINUE
C   PARTICLES SEE A CYLINDER OF LENGTH DZ = ZDF
        ELSEIF (NLTRZ) THEN
          DO 153 IR=1,NTRII
            VOL(IR)=AREA(IR)*ZDF
153       CONTINUE
C   PARTICLES SEE A TORUS APPROXIMATED BY NTTRAM STRAIGHT CYLINDERS
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 155 IR=1,NTRII
            VOL(IR)=AREA(IR)*(XCOM(IR)+RMTOR)*PI2AT
155       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.5) THEN
        TWOTHIRD=2.0D0/3.0D0
        DO ITET=1,NTET
C  CENTER OF MASS
          XTCEN(ITET)=0.D0
          YTCEN(ITET)=0.D0
          ZTCEN(ITET)=0.D0
          DO J=1,4
            IC=NTECK(J,ITET)
            XTCEN(ITET)=XTCEN(ITET) + XTETRA(IC)
            YTCEN(ITET)=YTCEN(ITET) + YTETRA(IC)
            ZTCEN(ITET)=ZTCEN(ITET) + ZTETRA(IC)
          END DO
          XTCEN(ITET)=XTCEN(ITET)*0.25D0
          YTCEN(ITET)=YTCEN(ITET)*0.25D0
          ZTCEN(ITET)=ZTCEN(ITET)*0.25D0

          IC1 = NTECK(1,ITET)
          IC2 = NTECK(2,ITET)
          IC3 = NTECK(3,ITET)
          IC4 = NTECK(4,ITET)
C  CALCULATE VOLUMES
          PC1(1:3)= (/ XTETRA(IC1), YTETRA(IC1), ZTETRA(IC1) /)
          PC2(1:3)= (/ XTETRA(IC2), YTETRA(IC2), ZTETRA(IC2) /)
          PC3(1:3)= (/ XTETRA(IC3), YTETRA(IC3), ZTETRA(IC3) /)
          PC4(1:3)= (/ XTETRA(IC4), YTETRA(IC4), ZTETRA(IC4) /)
          VOL(ITET) = CAL_VOL(PC1,PC2,PC3,PC4)
          IF (VOL(ITET) < -EPS10) THEN
            WRITE (6,*) ' WARNING ! '
            WRITE (6,*) ' VOL(',ITET,') < 0 VOL= ',VOL(ITET)
          END IF
          IF (SUM(NTBAR(1:4,ITET)) < 0) VOL(ITET) = 0._DP ! COLLAPSED TET
          VOL(ITET) = MAX(VOL(ITET),0._DP)
          AREA(ITET)=VOL(ITET)**TWOTHIRD
        END DO
C
      ELSEIF (LEVGEO.EQ.6) THEN
C
C  GENERAL GEOMETRY OPTION: PROVIDE CELL VOLUMES (CM**3)
C                      ON ARRAY VOL(IC),IC=1,NSURFM
C                     (ALSO PROVIDE CENTER OF CELL:
C                      XCOM(IC),YCOM(IC))
C
        CALL VOLUSR(NR1ST,VOL)
C
      ENDIF
C
190   CONTINUE
C
C
C
C  SET RADIAL SURFACE LABELING MESHES RHOSRF AND RHOZNE
C
      VOLSR=0.
      AREAR=0.
      IF (LEVGEO.EQ.1) THEN
C  RHOSRF(J) = RSURF(J)
        DO 195 J=1,NR1ST
          RHOSRF(J)=RSURF(J)
195     CONTINUE
      ELSEIF (LEVGEO.EQ.2) THEN
C  PI*RHOSRF(J)**2 = AREA ENCLOSED BY ORIGIN AND SURFACE J
        DO 196 J=1,NR1ST
          AREAR=AREAR+AREA1(J-1)
          RHOSRF(J)=SQRT(AREAR/PIA)
196     CONTINUE
      ELSEIF (LEVGEO.EQ.3) THEN
C  PI*RHOSRF(J)**2 = AREA ENCLOSED BY ORIGIN AND SURFACE J
        DO 197 J=1,NR1ST
          AREAR=AREAR+AREA1(J-1)
          RHOSRF(J)=SQRT(AREAR/PIA)
197     CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
C  RHOSRF AND RHOZNE ARE NOT DEFINED FOR FEM-OPTION
        RETURN
      ELSEIF (LEVGEO.EQ.5) THEN
C  RHOSRF AND RHOZNE ARE NOT DEFINED FOR TETRAHEDRON-OPTION
        RETURN
      ELSEIF (LEVGEO.EQ.6) THEN
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
        RETURN
      ENDIF
C
      DO 199 J=1,NR1STM
        RHOZNE(J)=0.5*(RHOSRF(J)+RHOSRF(J+1))
199   CONTINUE
C  MIRROR POINT
      RHOZNE(NR1ST)=RHOSRF(NR1ST)+0.5*(RHOSRF(NR1ST)-RHOSRF(NR1STM))
C
      IF (LEVGEO.LE.3) THEN
        CALL LEER(1)
        WRITE (6,*) 'FLUX-SURFACE LABELING GRIDS'
        CALL MASRR2('  N, RHOSRF,RHOZNE    ',
     .                    RHOSRF,RHOZNE,NR1ST)
        CALL LEER(2)
      ENDIF
C
      RETURN
C
C  2D (R-THETA OR X-Y) VOLUME ELEMENTS
C
200   CONTINUE
C
      IF (LEVGEO.EQ.1) THEN
C
        KP=1
        DO 220 I=1,NR1STM
          NCELL1 = I+(      (KP-1)*NP2T3)*NR1P2
          VSAVE = VOL(NCELL1)
          DO 220 J=1,NP2NDM
            FAC2=(PSURF(J+1)-PSURF(J))/YDF
            NCELLJ=I+((J-1)+(KP-1)*NP2T3)*NR1P2
            XCOM(NCELLJ)=RHOZNE(I)
            YCOM(NCELLJ)=PHZONE(J)
            VOL(NCELLJ)=VSAVE*FAC2
220     CONTINUE
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        IFLAG=1
        IT=1
        DO 240 IR=1,NR1STM
          IRP=IR+1
          DO 250 IP=1,NP2NDM
            IPP=IP+1
            CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                  TRI(IRP),TRI(IR),
     .                  RSURF(IRP),RSURF(IR),PSURF(IPP),PSURF(IP),IFLAG,
     .                  AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            AREA(NCELL)=AELL
            XCOM(NCELL)=SX
            YCOM(NCELL)=SY
250       CONTINUE
240     CONTINUE
C
        IF (NLTRT) THEN
          DO 262 I=1,NR1STM
            DO 262 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA(NCELL)*(XCOM(NCELL)+RMTOR)*PI2A
              IF (VOL(NCELL).GE.0.D0) GOTO 262
              WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT_OWN(1)
262       CONTINUE
        ELSEIF (NLTRZ) THEN
          DO 260 I=1,NR1STM
            DO 260 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA(NCELL)*ZDF
              IF (VOL(NCELL).GE.0.D0) GOTO 260
              WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT_OWN(1)
260       CONTINUE
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 261 I=1,NR1STM
            DO 261 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA(NCELL)*(XCOM(NCELL)+RMTOR)*PI2AT
              IF (VOL(NCELL).GE.0.D0) GOTO 261
              WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT_OWN(1)
261       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
        IF (NLTRT) THEN
          DO 265 I=1,NR1STM
            DO 265 JP=1,NPPLG
              DO 265 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*(XCOM(NCELL)+RMTOR)*PI2A
                IF (VOL(NCELL).GE.0.D0) GOTO 265
                WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
                CALL MASJ2('J,I             ',I,J)
C               CALL EXIT_OWN(1)
265       CONTINUE
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 267 I=1,NR1STM
            DO 267 JP=1,NPPLG
              DO 267 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*(XCOM(NCELL)+RMTOR)*PI2AT
                IF (VOL(NCELL).GE.0.D0) GOTO 267
                WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
                CALL MASJ2('J,I             ',I,J)
C               CALL EXIT_OWN(1)
267       CONTINUE
        ELSEIF (NLTRZ) THEN
          DO 268 I=1,NR1STM
            DO 268 JP=1,NPPLG
              DO 268 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*ZDF
268       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  FINITE ELEMENT OPTION: NOTHING TO BE DONE HERE
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
C  TETRAHEDRON OPTION: NOTHING TO BE DONE HERE
C
      ELSEIF (LEVGEO.EQ.6) THEN
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
      ENDIF
C
      RETURN
C
C    2D (R-Z), (R-PHI) OR (X-Z) VOLUME ELEMENTS
C
300   CONTINUE
C
      IF (NLTRZ.OR.NLTRA.OR.NLTRT) THEN
C
        IT=1
        DO 320 J=1,NP2ND
        DO 320 I=1,NR1ST
          NCELL1 = I+((J-1)            )*NR1P2
          VSAVE=VOL(NCELL1)
          DO 320 K=1,NT3RDM
            FAC3=(ZSURF(K+1)-ZSURF(K))/ZDF
            NCELLK=I+((J-1)+(K-1)*NP2T3)*NR1P2
            VOL(NCELLK)=VSAVE*FAC3
320     CONTINUE
C
      ENDIF
C
      RETURN
C
C   ADDITIONAL CELL VOLUMES
C
400   CONTINUE
C
C   ADDITIONAL CELL VOLUMES ARE DEFAULTED TO 1. AT PRESENT
C
      DO 410 J=NSURF+1,NSBOX
        VOL(J)=1.
410   CONTINUE
      IF (NLADD) THEN
        DO 411 J=NSURF+1,NSBOX
          IF (VOLADD(J-NSURF).GT.0.D0) VOL(J)=VOLADD(J-NSURF)
411     CONTINUE
      ENDIF

      DEALLOCATE (AREAP)
C
      RETURN
C
999   CONTINUE
      WRITE (6,*) 'UNWRITTEN OPTION CALLED IN SUBR. VOLUME '
      CALL EXIT_OWN(1)
      END