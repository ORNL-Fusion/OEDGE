C
      SUBROUTINE PLASMA_DERIV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CGRID'
      INCLUDE 'CLOGAU'
      INCLUDE 'CZT1'
      INCLUDE 'CGEOM'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCONA'
      DIMENSION ZTII(NPLS,NRAD)
c slmod begin - not tr
      COMMON /DENICOM/ DEINTF 
      DIMENSION        DEINTF(NPLS,NRAD)
c slmod end
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
C  SET ELECTRON-DENSITY FROM QUASI-NEUTRALITY, AND DRIFT ENERGY (EV)
      DO 5102 J=1,NSBOX
        DEIN(J)=0.
        DO 5102 IPLS=1,NPLSI
c slmod begin - not tr
c...      Electron density on the standard grid is passed from DIVIMP:
c          IF (J.LE.NSURF) THEN
c            DEIN(J)=DEIN(J)+DEINTF(IPLS,J)
c          ELSE
            DEIN(J)=DEIN(J)+DBLE(NCHRGP(IPLS))*DIIN(IPLS,J)
c          ENDIF
c
c          DEIN(J)=DEIN(J)+DBLE(NCHRGP(IPLS))*DIIN(IPLS,J)
c slmod end
          IF (NLDRFT) THEN
            EDRIFT(IPLS,J)=CVRSSP(IPLS)*
     .              (VXIN(IPLS,J)**2+VYIN(IPLS,J)**2+VZIN(IPLS,J)**2)
          ELSE
            EDRIFT(IPLS,J)=0.D0
          ENDIF
5102  CONTINUE
C
      DO 5103 J=1,NSBOX
C  SET 'LOG OF TEMPERATURE AND DENSITY' ARRAYS
        ZTEI=MAX(TVAC,MIN(TEIN(J),1.D10))
        TEINL(J)=LOG(ZTEI)
        ZTNE=MAX(DVAC,MIN(DEIN(J),1.D20))
        DEINL(J)=LOG(ZTNE)
C  SET 'VACUUM REGION FLAGS'
        LGVAC(J,0)=.TRUE.
        TMPLS=TEIN(J)
        DNPLS=DEIN(J)
        LGVAC(J,NPLS+1)=TMPLS.LE.TVAC.OR.DNPLS.LE.DVAC
        LGVAC(J,0)     =LGVAC(J,0).AND.LGVAC(J,NPLS+1)
        DO 5106 IPLS=1,NPLSI
          EMPLS=1.5*TIIN(IPLS,J)+EDRIFT(IPLS,J)
          DNPLS=DIIN(IPLS,J)
          LGVAC(J,IPLS)=EMPLS.LE.TVAC.OR.DNPLS.LE.DVAC
          LGVAC(J,0)   =LGVAC(J,0).AND.LGVAC(J,IPLS)
5106    CONTINUE
5103  CONTINUE
C
      DO 5105 IPLS=1,NPLSI
        DO 5105 J=1,NSBOX
          ZTII(IPLS,J)=MAX(TVAC,MIN(TIIN(IPLS,J),1.D10))
          TIINL(IPLS,J)=LOG(ZTII(IPLS,J))
          ZTNI=MAX(DVAC,MIN(DIIN(IPLS,J),1.D20))
          DIINL(IPLS,J)=LOG(ZTNI)
5105  CONTINUE
C
      IF (LEVGEO.EQ.3) THEN
        DO 5161 I=1,NPPLG-1
          DO 5162 IP=NPOINT(2,I),NPOINT(1,I+1)-1
            IPM=IP-1
            DO 5163 IPLS=0,NPLS+1
              DO 5163 IR=1,NR1STM
                IN=IR+IPM*NR1ST
                LGVAC(IN,IPLS)=.TRUE.
5163        CONTINUE
5162      CONTINUE
5161    CONTINUE
      ENDIF
C
C
C  ZT1: FOR "EFFECTIVE" PLASMA PARTICLE VELOCITY IN CROSS SECTIONS
C       AND SINGLE PARAMETER RATE COEF. FOR HEAVY PARTICLE INTERACTIONS
C       SQRT(ZT1) IS THE MEAN VELOCITY AT TI=ZTII
C
      DO 1210 IPLS=1,NPLSI
        FC=1./RMASSP(IPLS)*8./PIA*CVEL2A*CVEL2A
        DO 1210 J=1,NSBOX
          ZT1(IPLS,J)=FC*ZTII(IPLS,J)
1210  CONTINUE
C
C  ZRGQ: VARIANCE FOR SAMPLING FROM MAXWELLIAN VELOCITY DISTRIBUTION
C  ZRG=SQRT(ZRGQ) = STANDARD DEVIATION
C
      DO 1215 IPLS=1,NPLSI
        FC=CVEL2A/SQRT(RMASSP(IPLS))
        DO 1215 J=1,NSBOX
          ZRG(IPLS,J)=FC*SQRT(ZTII(IPLS,J))
1215  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PROFN(PRO,PRO0,PROS,P,Q,E,SEP,PROVAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      DIMENSION PRO(*)
      IF (LEVGEO.LE.2) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
        IF (LEVGEO.EQ.2) THEN
          JM1=LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFN       ')
          AR=AREAA (SEP,JM1,ARCA,YR,EP1R)
          RHOSEP=SQRT(AR*PIAI)
        ELSEIF (LEVGEO.EQ.1) THEN
          RHOSEP=SEP
        ENDIF
      ELSE
        WRITE (6,*) 'WARNING: SUBR. PROFN CALLED WITH LEVGEO.GT.2 '
        WRITE (6,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
      RDI=1./(RHOSEP-RHOSRF(1))
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).LE.RHOSEP) THEN
          ZR=RHOZNE(J)-RHOSRF(1)
          RORP=ZR*RDI
          IF (RORP.LE.0.D0) THEN
            PRO(J)=PRO0
          ELSE
            FACT=(1.0-RORP**P)**Q
            PRO(J)=(PRO0-PROS)*FACT+PROS
          ENDIF
        ELSE
          PRO(J)=PROS*EXP((RHOSEP-RHOZNE(J))/E)
        ENDIF
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
C
C
C
C
      SUBROUTINE PROFE(PRO,PRO0,RINP,A1,E,SEP,PROVAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      DIMENSION PRO(*)
      RHO=SEP
      RHOINP=RINP
      IF (LEVGEO.LE.2) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
        IF (LEVGEO.EQ.2) THEN
          JM1=LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFE       ')
          AR=AREAA (SEP,JM1,ARCA,YR,EP1R)
          RHO=SQRT(AR*PIAI)
          JM1=LEARCA(RINP,RSURF,1,NLOCAL,1,'PROFE       ')
          AR=AREAA (RINP,JM1,ARCA,YR,EP1R)
          RHOINP=SQRT(AR*PIAI)
        ENDIF
      ELSE
        WRITE (6,*) 'WARNING: PROFE CALLED WITH LEVGEO.GT.2 '
        WRITE (6,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).LE.RHO) THEN
          ZR=RHOZNE(J)-RHOINP
          PRO(J)=PRO0*EXP(-ZR/A1)
        ELSE
          DIST=RHO-RHOINP
          PROS=PRO0*EXP(-DIST/A1)
          PRO(J)=PROS*EXP((RHO-RHOZNE(J))/E)
        ENDIF
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
C
      SUBROUTINE PROFS(PRO,PRO0,PRO1,SEP,PROVAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      DIMENSION PRO(*)
      RHOSEP=SEP
      IF (LEVGEO.LE.2) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND RADIAL SURFACE LABELING CO-ORDINATE AT "SEP"
C
        IF (LEVGEO.EQ.2) THEN
          JM1=LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFS       ')
          AR=AREAA (SEP,JM1,ARCA,YR,EP1R)
          RHOSEP=SQRT(AR*PIAI)
        ELSEIF (LEVGEO.EQ.1) THEN
          RHOSEP=SEP
        ENDIF
      ELSE
        WRITE (6,*) 'WARNING: PROFS CALLED WITH LEVGEO.GT.2 '
        WRITE (6,*) 'CONSTANT PLASMA PARAMETERS RETURNED'
        NLOCAL=NR1ST
        DO 25 J=1,NLOCAL-1
          PRO(J)=PRO0
25      CONTINUE
        PRO(NR1ST)=0.
        RETURN
      ENDIF
C
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).GT.RHOSEP) GOTO 15
        PRO(J)=PRO0
        GOTO 20
15      PRO(J)=PRO1
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
C
C
      SUBROUTINE PROFR (PRO,IINDEX,NSPZI,NSPZ1,NDAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CSPEI'
      DIMENSION PRO(NSPZ1,*)
      DO 20 ISPZ=1,NSPZI
        DO 20 J=1,NDAT
          PRO(ISPZ,J)=RWK(IINDEX*NRAD+ISPZ+(J-1)*NSPZ1)
20    CONTINUE
      RETURN
      END
C
      SUBROUTINE VOLUME (IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'CPOLYG'
      INCLUDE 'CGEOM'
      INCLUDE 'CCONA'
      INCLUDE 'CLOGAU'
      INCLUDE 'COMUSR'
      INCLUDE 'CTRIG'
      INCLUDE 'CUPD'
      DIMENSION VSAVE(N1ST*N2ND)
      DIMENSION AREA2(N1ST,N2ND),AREAP(N1STS,N2NDPLG),AREAT(NTRI),
     .          AREA1(0:N1ST)
      DIMENSION XPH(N1ST,N2ND),YPH(N1ST,N2ND)
      DIMENSION XCOMTR(NTRIS),YCOMTR(NTRIS)
      EQUIVALENCE (AREA(1),AREA2(1,1))
!PB   EQUIVALENCE (AREA(1),AREAP(1,1))
      EQUIVALENCE (AREA(1),AREAT(1))
      EQUIVALENCE (XCOM(1,1),XCOMTR(1)),
     .            (YCOM(1,1),YCOMTR(1))
c slmod begin - tr
      DIMENSION XRAD(NRAD)
c slmod end
      SAVE

C     IND=1: 1-ST GRID, RAD. RESOLUTION
C     IND=2: 2-ND GRID, POL. RESOLUTION
C     IND=3: 3-RD GRID, TOR. RESOLUTION
C     IND=4: ADDITIONAL CELL REGION
C
      IF (output) WRITE(0,*) 'MARK: IND=',ind
      GOTO(100,200,300,400),IND
C
100   CONTINUE
C
      IF (output) WRITE(0,*) 'MARK: RADIAL GRID'

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
            CALL EXIT
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
c slmod begin - grid - tr
            IF (GRIDOPT.EQ.1) THEN
              AR=ARTRIA(0.D0,0.D0,XVERT(1,J  ,1),YVERT(1,J  ,1),
     .                            XVERT(1,J+1,1),YVERT(1,J+1,1))
            ELSE
              AR=ARTRIA(0.D0,0.D0,XPOL(1,J),YPOL(1,J),
     .                            XPOL(1,J+1),YPOL(1,J+1))
            ENDIF
c
c          AR=ARTRIA(0.D0,0.D0,XPOL(1,J),YPOL(1,J),
c    .                        XPOL(1,J+1),YPOL(1,J+1))
c slmod end
          AREA1(0)=AREA1(0)+AR
139     CONTINUE
        AREA1(0)=ABS(AREA1(0))
C
c slmod begin - tr
        CALL ARPOLY(XPOL,YPOL,XVERT,YVERT,N2ND,NRPLG,N1STS,1,
     .              NR1ST,AREAP,XCOM,YCOM)
c
c        CALL ARPOLY(XPOL,YPOL,NRPLG,N1STS,1,NR1ST,AREAP,XCOM,YCOM)
c slmod end
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
                XC=XC+AREAP(IR,J)*XCOM(IR,J)
                AR=AR+AREAP(IR,J)
135         CONTINUE
            XC=XC/(AR+EPS60)
            VOL(IR)=AREA1(IR)*(XC+RMTOR)*PI2A
134       CONTINUE
        ELSEIF (NLTRZ) THEN
C   PARTICLES SEE A CYLINDER OF LENGTH DZ = ZDF
          IF (output)
     .    WRITE(6,*) 'MARK: CYLINDER',zdf,nr1stm,nsbox
c          WRITE(6,*) 'MARK:- CYLINDER',zdf,nr1stm,nsbox
          DO 133 IR=1,NR1STM
            VOL(IR)=AREA1(IR)*ZDF

c              WRITE(6,'(A,I6,1P,2E12.4,0P)')
c     .          'MARK: AREA1 ',ir,area1(ir),vol(ir)
     
c            WRITE(6,'(A,I6,1P,2E12.4,0P)')
c     .        'MARK:- AREA1 ',ir,area1(ir),vol(ir)
133       CONTINUE
C   PARTICLES SEE A TORUS APPROXIMATED BY NTTRAM STRAIGHT CYLINDERS
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 136 IR=1,NR1STM
            XC=0.
            AR=0.
            DO 137 JP=1,NPPLG
              DO 137 J=NPOINT(1,JP),NPOINT(2,JP)-1
                XC=XC+AREAP(IR,J)*XCOM(IR,J)
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
          XCOMTR(IR) = (XTRIAN(NECKE(1,IR))+XTRIAN(NECKE(2,IR))+
     .                  XTRIAN(NECKE(3,IR)))/3.
          YCOMTR(IR) = (YTRIAN(NECKE(1,IR))+YTRIAN(NECKE(2,IR))+
     .                  YTRIAN(NECKE(3,IR)))/3.
          AR=0.5*(XTRIAN(NECKE(2,IR))*(YTRIAN(NECKE(3,IR))
     >           -YTRIAN(NECKE(1,IR)))+XTRIAN(NECKE(3,IR))*
     >           (YTRIAN(NECKE(1,IR))
     >           -YTRIAN(NECKE(2,IR)))+XTRIAN(NECKE(1,IR))*
     >           (YTRIAN(NECKE(2,IR))-YTRIAN(NECKE(3,IR))))
          AREAT(IR) = AR
150     CONTINUE
C
C
C   PARTICLES SEE A TORUS
        IF (NLTRT) THEN
          DO 154 IR=1,NTRII
            VOL(IR)=AREAT(IR)*(XCOMTR(IR)+RMTOR)*PI2A
154       CONTINUE
C   PARTICLES SEE A CYLINDER OF LENGTH DZ = ZDF
        ELSEIF (NLTRZ) THEN
          DO 153 IR=1,NTRII
            VOL(IR)=AREAT(IR)*ZDF
153       CONTINUE
C   PARTICLES SEE A TORUS APPROXIMATED BY NTTRAM STRAIGHT CYLINDERS
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 155 IR=1,NTRII
            VOL(IR)=AREAT(IR)*(XCOMTR(IR)+RMTOR)*PI2AT
155       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
C  GENERAL GEOMETRY OPTION: PROVIDE CELL VOLUMES (CM**3)
C                      ON ARRAY VOL(IC),IC=1,NSURFM
C                     (ALSO PROVIDE CENTER OF CELL:
C                      XCOM(IC,1),YCOM(IC,1))
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
      IT=1
      DO 201 IR=1,NR1ST
        DO 201 IP=1,NP2ND
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          VSAVE(NCELL)=VOL(NCELL)
201   CONTINUE
      DO 202 ISURF=1,NSURF
        VOL(ISURF)=0.
202   CONTINUE
C
      IF (LEVGEO.EQ.1) THEN
C
        KP=1
        DO 220 J=1,NP2NDM
          FAC2=(PSURF(J+1)-PSURF(J))/YDF
          DO 220 I=1,NR1STM
            XCOM(I,J)=RHOZNE(I)
            YCOM(I,J)=PHZONE(J)
            NCELLJ=I+((J-1)+(KP-1)*NP2T3)*NR1P2
            NCELL1=I+(      (KP-1)*NP2T3)*NR1P2
            VOL(NCELLJ)=VSAVE(NCELL1)*FAC2
220     CONTINUE
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        IFLAG=1
        DO 240 IR=1,NR1STM
          IRP=IR+1
          DO 250 IP=1,NP2NDM
            IPP=IP+1
            CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                  TRI(IRP),TRI(IR),
     .                  RSURF(IRP),RSURF(IR),PSURF(IPP),PSURF(IP),IFLAG,
     .                  AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
            AREA2(IR,IP)=AELL
            XCOM(IR,IP)=SX
            YCOM(IR,IP)=SY
250       CONTINUE
240     CONTINUE
C
        IF (NLTRT) THEN
          DO 262 I=1,NR1STM
            DO 262 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA2(I,J)*(XCOM(I,J)+RMTOR)*PI2A
              IF (VOL(NCELL).GE.0.D0) GOTO 262
              WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT
262       CONTINUE
        ELSEIF (NLTRZ) THEN
          DO 260 I=1,NR1STM
            DO 260 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA2(I,J)*ZDF
              IF (VOL(NCELL).GE.0.D0) GOTO 260
              WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT
260       CONTINUE
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 261 I=1,NR1STM
            DO 261 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA2(I,J)*(XCOM(I,J)+RMTOR)*PI2AT
              IF (VOL(NCELL).GE.0.D0) GOTO 261
              WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT
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
                VOL(NCELL)=ABS(AREAP(I,J))*(XCOM(I,J)+RMTOR)*PI2A
                IF (VOL(NCELL).GE.0.D0) GOTO 265
                WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
                CALL MASJ2('J,I             ',I,J)
C               CALL EXIT
265       CONTINUE
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 267 I=1,NR1STM
            DO 267 JP=1,NPPLG
              DO 267 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*(XCOM(I,J)+RMTOR)*PI2AT
c slmod begin - tr
                VOL(NCELL)=VOL(NCELL)*VOLCOR2
c slmod end
                IF (VOL(NCELL).GE.0.D0) GOTO 267
                WRITE (6,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
                CALL MASJ2('J,I             ',I,J)
C               CALL EXIT
267       CONTINUE
        ELSEIF (NLTRZ) THEN
          DO 268 I=1,NR1STM
            DO 268 JP=1,NPPLG
              DO 268 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*ZDF
c                WRITE(6,'(A,2I4,2X,I5,1P,2E12.4,0P)')
c     .            'MARK:- CYL2: ',i,j,ncell,ABS(AREAP(I,J)),vol(ncell)
268       CONTINUE
c slmod begin - tr
c...      Calculate major radius value for the center-points of the
c         cells on the standard grid:
          CALL DSET(XRAD,NRAD,0.0D0)
          DO IR=1,NR1STM
            DO IP=1,NRPLG-1
              IN=IR+(IP-1)*NR1ST
              IF (GRIDOPT.EQ.1) THEN
                XRAD(IN)=0.25D0*(XVERT(IR,IP  ,1)+XVERT(IR,IP+1,1)+
     .                           XVERT(IR,IP+1,2)+XVERT(IR,IP  ,2))
              ELSE
                XRAD(IN)=0.25D0*(XPOL(IR  ,IP  )+XPOL(IR  ,IP+1)+
     .                           XPOL(IR+1,IP+1)+XPOL(IR+1,IP  ))
              ENDIF
            ENDDO
          ENDDO
c...      Calculation of toroidal volume for volume recombination:
          DO J=1,NSBOX
            VOL2(J)=2.0D0*3.14D0*XRAD(J)*AREA(J)*EIRTORFRAC
          ENDDO
c slmod end
        ENDIF
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  FINITE ELEMENT OPTION: NOTHING TO BE DONE HERE
C
      ELSEIF (LEVGEO.EQ.5) THEN
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
      IT=1
      DO 301 IR=1,NR1ST
        DO 301 IP=1,NP2ND
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          VSAVE(NCELL)=VOL(NCELL)
301   CONTINUE
      DO 302 ISURF=1,NSURF
        VOL(ISURF)=0.
302   CONTINUE
C
      IF (NLTRZ.OR.NLTRA.OR.NLTRT) THEN
C
        DO 320 K=1,NT3RDM
          FAC3=(ZSURF(K+1)-ZSURF(K))/ZDF
          DO 320 J=1,NP2ND
          DO 320 I=1,NR1ST
            NCELLK=I+((J-1)+(K-1)*NP2T3)*NR1P2
            NCELL1=I+((J-1)            )*NR1P2
            VOL(NCELLK)=VSAVE(NCELL1)*FAC3
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
      IF (output) WRITE(0,*) 'MARK: ADDITIONAL CELL VOLUMES'
      DO 410 J=NSURF+1,NSBOX
        VOL(J)=1.
410   CONTINUE
      IF (NLADD) THEN
        DO 411 J=NSURF+1,NSBOX
          IF (VOLADD(J-NSURF).GT.0.D0) VOL(J)=VOLADD(J-NSURF)
411     CONTINUE
      ENDIF
C
      RETURN
C
999   CONTINUE
      WRITE (6,*) 'UNWRITTEN OPTION CALLED IN SUBR. VOLUME '
      CALL EXIT
      END
C
      SUBROUTINE MCARLO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  MONTE CARLO CALCULATION
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CLGIN'
      INCLUDE 'CRAND'
      INCLUDE 'COMSOU'
      INCLUDE 'CSPEI'
      INCLUDE 'CSDVI'
      INCLUDE 'CSDVI_BGK'
      INCLUDE 'CSDVI_COP'
      INCLUDE 'CSPEZ'
      INCLUDE 'CESTIM'
      INCLUDE 'CGEOM'
      INCLUDE 'CGRID'
      INCLUDE 'COMSPL'
      INCLUDE 'COMUSR'
      INCLUDE 'COUTAU'
      INCLUDE 'CLOGAU'
      INCLUDE 'CTRCEI'
      INCLUDE 'CCONA'
      INCLUDE 'CZT1'
      INCLUDE 'COMNNL'
      INCLUDE 'CCOUPL'
c slmod begin - temp - not tr
c...  For the call to USRBGK
c      INCLUDE 'COMXS'
c slmod end
C
      CHARACTER*6 CIS
      CHARACTER*10 CDATE,CTIME
      LOGICAL NLPOLS,NLTORS
      DIMENSION RPST(NPARTC),IPST(MPARTC)
      EQUIVALENCE (RPST(1),X0),(IPST(1),IPOLG)
C
      DIMENSION OUTAU(NOUTAU),ESTIM(NESTIM),SDVI(NSDVI),
     .          SDVI_BGK(NSBGK),SDVI_COP(NSCOP)
      LOGICAL SPEZ(LSPEZ)
C
      DIMENSION
     . STV(NSD,NRAD),STVW(NSDW,NLIMPS),STVC(0:2,NCV,NRAD),
     . STVS(NSD),STVWS(NSDW),STVCS(0:2,NCV),
     . COVS(NCV,NHST)
C
      DIMENSION
     . EE(NSD,NRAD),FF(NSDW,NLIMPS),
     . EES(NSD),FFS(NSDW),
     . SDVIA(NSD,NRAD),SDVIAW(NSDW,NLIMPS),SDVIAC(2,NCV,NRAD)
      DIMENSION
     . ZVOLIN(NRAD),ZVOLIW(NRAD),
     . XTIM(0:NSTRA),SCLTAL(N1MX,NTALV)
C
C   EQUIVALENCE: COMMON: COUTAU,CESTIM,CSDVI,CSPEZ   ONTO 1 ARRAY EACH
C
      EQUIVALENCE
     .  (OUTAU(1),PDENAI(0,0)),
     .  (ESTIM(1),PDENA(1,1)),
     .  (SDVI(1),SIGMA(1,1)),
     .  (SDVI_BGK(1),SIGMA_BGK(1,1)),
     .  (SDVI_COP(1),SIGMA_COP(1,1)),
     .  (SPEZ(1),LOGATM(0,0))
C
C  EQUIVALENCE: SUM OVER STRATA, RWK(1:NID1)
C
C  RWK(1)        UNTIL RWK(NESTM1)        : SUM OVER VOLUME TALLIES
C  RWK(NESTM1+1) UNTIL RWK(NESTIM)        : SUM OVER SURFACE TALLIES
C  RWK(NESTIM+1) UNTIL RWK(NESTIM+NSDVI)  : SUM OVER STANDARD DEVIATIONS
C
      EQUIVALENCE
     .(RWK(1+NESTIM*NSMSTRA                          ),STV(1,1)),
     .(RWK(1+NESTIM*NSMSTRA+NSD* NRAD                ),STVW(1,1)),
     .(RWK(1+NESTIM*NSMSTRA+NSD* NRAD   +NSDW* NLIMPS),STVS(1)),
     .(RWK(1+NESTIM*NSMSTRA+NSD*(NRAD+1)+NSDW* NLIMPS),STVWS(1)),
     .(RWK(1+NESTIM*NSMSTRA+NSD*(NRAD+1)+NSDW*(NLIMPS+1)),STVC(0,1,1)),
     .(RWK(1+NESTIM*NSMSTRA+NSD*(NRAD+1)+NSDW*(NLIMPS+1)+3*NCV*NRAD),
     .                                               STVCS(0,1))
C  LAST INDEX: NID1
C
C   EQUIVALENCE: WORK ARRAYS FOR STANDARD DEVIATION
C TO BE SET TO ZERO ONLY ONCE
      EQUIVALENCE
     . (RWK(NID3P                                        ),EE(1,1)),
     . (RWK(NID3P+NSD* NRAD                              ),FF(1,1)),
     . (RWK(NID3P+NSD* NRAD   +NSDW* NLIMPS   +2*NCV*NRAD),EES(1)),
     . (RWK(NID3P+NSD*(NRAD+1)+NSDW* NLIMPS   +2*NCV*NRAD),FFS(1))
C  LAST INDEX:  NNID3
      PARAMETER(NNID3=
     .          NID3+NSD*(NRAD+1)+NSDW*(NLIMPS+1)+2*NCV*(NRAD+1))
C
C
C TO BE SET TO ZERO FOR EACH STRATUM
C
      EQUIVALENCE
     . (RWK(NNID3+1                     ),SDVIA(1,1)),
     . (RWK(NNID3+1+NSD*NRAD            ),SDVIAW(1,1)),
     . (RWK(NNID3+1+NSD*NRAD+NSDW*NLIMPS),SDVIAC(1,1,1))
C  LAST INDEX:  NNID4
      PARAMETER(NNID4=NNID3+NSD*NRAD+NSDW*NLIMPS+2*NCV*NRAD)
C
C
C   EQUIVALENCE: LOCAL IN SUBR. MCARLO
C
      EQUIVALENCE
     . (RWK(NNID4+1),ZVOLIN(1)),
     . (RWK(NNID4+1+1*NRAD),ZVOLIW(1))
C  LAST INDEX:  NNID4+2*NRAD
C
      LOGICAL LGSTOP
      SAVE XTIM,XX1,SECND
      DATA N2/2/
c slmod begin - tr
      COMMON /TRASHCOM/ nlost
      INTEGER           nlost

      COMMON /STATSCOM/ npar,xpos,ypos,zpos,xvel,yvel,zvel
      INTEGER           npar(9)
      REAL              xpos(9),ypos(9),zpos(9),
     .                  xvel(9),yvel(9),zvel(9)

      COMMON /SLTEMP2/ CXCNT
      INTEGER          CXCNT

      DOUBLE PRECISION slltime,slntime

      INTEGER parcnt(5)
      CHARACTER*12 partag(5)

      partag(1) = 'ATOM'
      partag(2) = 'MOLECULE'
      partag(3) = 'TEST ION'
      partag(4) = ' ??? '
      partag(5) = ' ??? '

      DO i = 1, 9
        npar(i) = 0
        xpos(i) = 0.0
        ypos(i) = 0.0
        zpos(i) = 0.0
        xvel(i) = 0.0
        yvel(i) = 0.0
        zvel(i) = 0.0

      ENDDO

      DO i = 1, 5
        parcnt(i) = 0
      ENDDO

      nlost = 0


c      DIVSUR = NR1ST/2
      DIVSUR = 1


      WRITE(6,*) 
      WRITE(6,*) '********************'
      WRITE(6,*) 
      WRITE(6,*) 'DIVSUR=',DIVSUR
      WRITE(6,*) 
      WRITE(6,*) '********************'
      WRITE(6,*) 

c slmod end
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      IF (NFILEN.NE.0) THEN
        NREC10=1500
        OPEN (UNIT=10,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC10)
        NREC11=NOUTAU
        OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
      ENDIF
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C-------------------------------------------------------------------
C
C** INITIALIZE SOME DATA AND SUBROUTINES (ONCE FOR ALL STRATA) *****
C
C  SCLTAL: FLAG FOR SCALING OF VOLUME AVERAGED TALLY
C  SCLTAL =0  1.
C         =1  ZVOLIN(ICELL)
C         =2  ZW
C         =3  ZVOLIW(ICELL)
C         =4  ZWW
C
      CALL ZEROA2(SCLTAL,N1MX,NTALV)
      SCLTAL(1,1)=1
      SCLTAL(1,2)=1
      SCLTAL(1,3)=1
      SCLTAL(1,4)=1
      SCLTAL(1,5)=1
      SCLTAL(1,6)=1
      SCLTAL(1,7)=3
      SCLTAL(1,8)=3
      SCLTAL(1,9)=3
      SCLTAL(1,10)=3
      SCLTAL(1,11)=3
      SCLTAL(1,12)=3
      SCLTAL(1,13)=3
      SCLTAL(1,14)=3
      SCLTAL(1,15)=3
      SCLTAL(1,16)=3
      SCLTAL(1,17)=3
      SCLTAL(1,18)=3
      SCLTAL(1,19)=3
      SCLTAL(1,20)=3
      SCLTAL(1,21)=3
      SCLTAL(1,22)=3
      SCLTAL(1,23)=3
      SCLTAL(1,24)=3
      SCLTAL(1,25)=3
      SCLTAL(1,26)=3
      SCLTAL(1,27)=3
      SCLTAL(1,28)=3
      SCLTAL(1,29)=3
      SCLTAL(1,30)=3
      SCLTAL(1,31)=3
      SCLTAL(1,32)=3
      SCLTAL(1,33)=3
      SCLTAL(1,34)=3
      SCLTAL(1,35)=3
      SCLTAL(1,36)=3
C
C  INITIALIZE RANDOM NUMBER ARRAYS
      INIV1=0
      INIV2=0
      INIV3=0
      INIV4=0
C  DETERMINE MAXIMAL INTEGER (DEPENDING ON MACHINE)
      IF (NLCRR) THEN
        INTMAX=ILARG(NMACH)
      ENDIF
C
C  IRNDVC MUST BE EVEN AND NOT LARGER THEN 64 (COMMON CRAND)
C  IRNDVC IS THE NUMBER OF RANDOM VECTORS PRODUCED IN ONE CALL TO
C  TO RANDOM SAMPLING ROUTINES
      IF (NLCRR) THEN
        IRNDVC=2
      ELSE
        IRNDVC=64
      ENDIF
      IRNDVH=IRNDVC/2
C
C  INITIALIZE SUBR. STATIS
      CALL STATS0
      CALL STATS0_BGK
      CALL STATS0_COP
C  INITIALIZE SUBR. COVAR
      CALL COVAR0
C  INITIALIZE SUBR. REFLEC AND SPUTER
      CALL REFLC0
      CALL SPUTR0
C  INITIALIZE SUBR. SAMVOL
      CALL SAMVL0
c slmod begin - not tr
      IF (.FALSE..AND.NITER.EQ.1) THEN
c
c        WRITE(0,*) 'TRYING TO RESTART BGK COLLISIONS'
c        WRITE(6,*) 'TRYING TO RESTART BGK COLLISIONS'
c
c        CALL RSTRT(0,NSTRAI,NESTIM,NSDVI,ESTIM,SDVI,
c     .             NSBGK,SDVI_BGK,NSCOP,SDVI_COP,TRCFLE)
c
c        CALL MODBGK
c
      ENDIF
c slmod end
C
      IESTR=-1
      IF (NFILEN.EQ.2.OR.NFILEN.EQ.7) GOTO 2000
C
C**** CLEAR WORK AREA FOR SUM OVER STRATA ****************************
C
      DO J=1,NID1
        RWK(J)=0.
      ENDDO
C  ARRAYS: EE,FF,GG,....
      NI=NID3P
      NE=NNID3
      DO J=NI,NE
        RWK(J)=0.
      ENDDO
C  BGK-ARRAYS:
      IF (NSIGI_BGK.GT.0) THEN
        DO IB=1,NBGVI_STAT
          STVS_BGK(IB)=0.D0
          EES_BGK(IB)=0.D0
          DO IR=1,NRAD
            STV_BGK(IB,IR)=0.D0
            EE_BGK(IB,IR)=0.D0
          ENDDO
        ENDDO
      ENDIF
C  COP-ARRAYS:
      IF (NSIGI_COP.GT.0) THEN
        DO IC=1,NCPVI_STAT
          STVS_COP(IC)=0.D0
          EES_COP(IC)=0.D0
          DO IR=1,NRAD
            STV_COP(IC,IR)=0.D0
            EE_COP(IC,IR)=0.D0
          ENDDO
        ENDDO
      ENDIF
C
C**** INITIALIZE COMMONS COUTAU AND CSPEZ
C
      DO 3 J=1,NOUTAU
3       OUTAU(J)=0.
      FASCL(0)=1.
      FMSCL(0)=1.
      FISCL(0)=1.
C
      DO 4 J=1,LSPEZ
4       SPEZ(J)=.FALSE.
C
C
C   MAXIMAL CALCULATION TIME ALLOWED FOR EACH STRATUM,
C   PROPORTIONAL TO NPTS(ISTRA), OR FLUX(ISTRA) (INPUT)
C   OR LINEAR COMBINATION THEREOF
C   THEREFORE NUMBER OF TEST PARTICLES MAY BE LESS THAN NPTS
C   BUT DO AT LEAST 2 PARTICLES, IN CASE NPTS(ISTRA).GE.2
C
      CALL TRMAIN(XX,NTCPU)
      XTIM(0)=SECOND_OWN()
      SECND=XTIM(0)
C  REMAINING CPU TIME
      XX1=XX-N2
      XPT=0.
      XFL=0.
      DO 7 ISTRA=1,NSTRAI
        IF (NPTS(ISTRA).LE.0.AND.FLUX(ISTRA).GT.0.D0) THEN
          FLUX(ISTRA)=0.D0
          WRITE (6,*) 'STRATUM ISTRA= ',ISTRA,' TURNED OFF. ZERO NPTS'
          CALL LEER(1)
        ENDIF
        IF (NPTS(ISTRA).GT.0.AND.FLUX(ISTRA).LE.0.D0) THEN
          NPTS(ISTRA)=0
          WRITE (6,*) 'STRATUM ISTRA= ',ISTRA,' TURNED OFF. ZERO FLUX'
          CALL LEER(1)
        ENDIF
        XPT=XPT+FLOAT(NPTS(ISTRA))
        XFL=XFL+FLUX(ISTRA)
7     CONTINUE
      XPT1=0.
      XFL1=0.
      DO 8 ISTRA=1,NSTRAI
        XPT1=XPT1+NPTS(ISTRA)
        XFL1=XFL1+FLUX(ISTRA)
        XTIM(ISTRA)=XTIM(0)+XX1*((1.-ALLOC)*XPT1/(XPT+EPS60)+
     +                           (   ALLOC)*XFL1/(XFL+EPS60))
8     CONTINUE
c slmod begin
      DO ISTRA=1,NSTRAI
        IF (XTIM(ISTRA)-XTIM(ISTRA-1).LT.1.0D0) THEN
          DO ISTRA2 = NSTRAI, ISTRA, -1
            XTIM(ISTRA2) = XTIM(ISTRA2) + 
     .                     (1.0D0 - (XTIM(ISTRA)-XTIM(ISTRA-1)))
          ENDDO
        ENDIF
      ENDDO
c slmod end   
C
c slmod begin - tr
      WRITE(6,*) 'NTIME       =',ntime
      WRITE(6,*) 'XX XX1      =',xx,xx1
      WRITE(6,*) 'N2 NSTRAI   =',n2,nstrai
      WRITE(6,*) 'XLPT1 XFL1  =',xpt1,xfl1
      WRITE(6,*) 'ALLOC       =',alloc
      WRITE(6,*) 'XPT XFL     =',xpt,xfl
      WRITE(6,*) 'NPTS12      =',npts(1),npts(2)
      WRITE(6,*) 'FLUX12      =',flux(1),flux(2)

      CALL CLOCK(slltime)
c slmod end
      CALL MASAGE ('LOOP OVER STRATA STARTS AT CPU TIME(SEC):    ')
      CALL MASR1 ('STARTTIM',XTIM(0))
      CALL MASAGE ('CPU TIME ASSIGNED TO STRATA (SEC) :          ')
      IF (ALLOC.EQ.0.D0) THEN
        CALL MASAGE ('PROPORTIONAL NPTS(ISTRA)                     ')
      ELSEIF (ALLOC.EQ.1.) THEN
        CALL MASAGE ('PROPORTIONAL FLUX(ISTRA)                     ')
      ELSE
        CALL MASAGE ('WEIGHTED ALLOCATION BETWEEN NPTS AND FLUX    ')
      ENDIF
      DO 9 ISTRA=1,NSTRAI
        DELT=XTIM(ISTRA)-XTIM(ISTRA-1)
9       CALL MASJ1R ('STRATUM,TIME    ',ISTRA,DELT)
      CALL LEER(2)
C
C  ASSIGN NUMBER OF PARTICLES TO BE STORED ON CENSUS, PROPORTIONAL
C  TO CPU TIME ASSIGNED TO EACH STRATUM
C
      IF (NPRNLI.GT.0) THEN
        WRITE(6,*) 'MAXIMUM NUMBER OF PARTICLES THAT WILL BE SAVED '
        WRITE(6,*) 'FOR SNAPSHOT ESTIMATORS: PROPORTIONAL TO CPU-'
        WRITE(6,*) 'TIME ALLOCATED FOR EACH STRATUM'
        DO 10 ISTRA=1,NSTRAI
c slmod begin
c...      
c
c
          NPRNLS(ISTRA)=NPRNLI*(XTIM(ISTRA)-XTIM(ISTRA-1))/XX1
          NPRNLS(ISTRA)=MAX(2,NPRNLS(ISTRA))
c
c          NPRNLS(ISTRA)=NPRNLI*(XTIM(ISTRA)-XTIM(ISTRA-1))/XX1
c slmod end 
          WRITE(*,*) 'ISTRA, NUMBER = ',ISTRA,NPRNLS(ISTRA)
10      CONTINUE
      ENDIF
c slmod begin - tr
c...note: This should keep stratum 1 from begin lost on short cases.
      rdum1 = reset_second()
c slmod end
C
C**** STRATA LOOP ****************************************************
C
      NPANU=0
      DO 1000 ISTRA=1,NSTRAI
c slmod begin - tr - new
        DO i = 1, 5
          parcnt(i) = 0
        ENDDO
c slmod end
        IF (output)
     .  WRITE(0,*) 'MARK: MCARLO: START OF STRATA LOOP'
        CALL LEER(2)
        WRITE (6,*) 'BEGIN TO WORK ON STRATUM NO. ',ISTRA
        CALL LEER(2)
        IPANU=0
        XMCP(ISTRA)=0.
C
C  INITIALIZE RANDOM NUMBER GENERATOR FOR STRATUM ISTRA
c
c  NOTE ********************************* 
c
c  jdemod - the problem with EIRENE reproducibility does not appear
c           to be a problem with the random number generator but rather 
c           with the time allocation which is system load dependent
c           among other factors - so whent the number of particles/strata
c           depends on timing that is different for each case, then the 
c           results are also different. 
c
c       *********************************        
c     
c     A random number seed value of 2 sets the seed to 2 AND also sets the number of particles/strata to 1000
c     The objective is to obtain repeatable EIRENE runs for use in test cases        
c     A random number seed value of 3 sets the seed to 3 AND also sets the number of particles/strata to 10000
c     The objective is to obtain repeatable EIRENE runs for use in test cases        
c        
        if (ninitl(istra).eq.2) then
           npts(istra) = 1000
           write(0,*) 'EIRENE99: SPECIAL RANDOM NUMBER SEED = 2 :'//
     >               'STRATA=',istra,' :PARTICLES SET TO 1000/STRATA'
     >          //': DEBUG USE ONLY'          
        elseif (ninitl(istra).eq.3) then
           npts(istra) = 10000
           write(0,*) 'EIRENE99: SPECIAL RANDOM NUMBER SEED = 3 :'// 
     >               'STRATA=',istra,' :PARTICLES SET TO 10000/STRATA'
     >          //': DEBUG USE ONLY'          
        endif
c
        IF (NINITL(ISTRA).GT.0) THEN
          NINIST=NINITL(ISTRA)
          IF (output) THEN
          WRITE(0,*) 'MARK: SEED=',NINIST
          ENDIF
          WRITE(6,*) 'MARK: SEED=',NINIST
          CALL RANSET(NINIST)
          CALL RANGET(ISEED)
          INIV1=0
          INIV2=0
          INIV3=0
          INIV4=0
        ELSEIF (NINITL(ISTRA).LT.0) THEN
          CALL DATE_AND_TIME(CDATE,CTIME)
          READ(CTIME(1:6),*) NINITL(ISTRA)
          WRITE (6,*) 'MARK: NINITL(ISTRA) SET TO ',NINITL(ISTRA)
          IF (output)
     .    WRITE (0,*) 'MARK: NINITL(ISTRA) SET TO ',NINITL(ISTRA)
          NINIST=NINITL(ISTRA)
          CALL RANSET(NINIST)
          CALL RANGET(ISEED)
          INIV1=0
          INIV2=0
          INIV3=0
          INIV4=0
C       ELSEIF (NINITL(ISTRA).EQ.0) THEN
C  DON'T INITALIZE FOR THIS STRATUM, NOTHING TO BE DONE HERE
        ENDIF
C
        FASCL(ISTRA)=1.
        FMSCL(ISTRA)=1.
        FISCL(ISTRA)=1.
C
C  CLEAR WORK AREA FOR THIS STRATUM
C
        IESTR=-1
        DO 31 J=1,NESTIM
31        ESTIM(J)=0.

c...sltmp
        DO i1 = 0, NMOMCHA
          DO i2 = 1, NRAD
            DO i3 = 1, NCPV
              copv2(i3,i2,i1) = 0.0
            ENDDO
          ENDDO
        ENDDO
        SUM1 = 0.0

        DO 32 J=1,NSDVI
32        SDVI(J)=0.
C  ARRAYS SDVIA,SDVIAW,SDVIAC
        NI=NNID3+1
        NE=NNID4
        DO 33 J=NI,NE
33        RWK(J)=0.
C  BGK-ARRAYS
        SGMS_BGK=0.D0
        SIGMA_BGK=0.D0
        SDVIA_BGK=0.D0
C  COP-ARRAYS
        SGMS_COP=0.D0
        SIGMA_COP=0.D0
        SDVIA_COP=0.D0
C
C  ENFORCE TOROIDAL OR POLOIDAL SYMMETRY FOR THIS STRATUM
C
        IF (NLAVRP(ISTRA)) THEN
          NLPOLS=NLPOL
          NLPOL=.FALSE.
        ENDIF
C
        IF (NLAVRT(ISTRA)) THEN
          NLTORS=NLTOR
          NLTOR=.FALSE.
        ENDIF
C
        IPRNLS=0
C
C  INITIALIZE SUBR. LOCATE
        CALL LOCAT0
C  INITIALIZE SUBR. SAMPLE
        IF (NLSRF(ISTRA)) CALL SAMSF0
C
        IF (NPTS(ISTRA).LE.0) GOTO 1000
C
C  LOCATE AND FOLLOW MC-PARTICLES
C
        CALL FTCRI(ISTRA,CIS)
        CALL MASBOX ('LAUNCH PARTICLES FOR STRATUM NUMBER ISTRA='//CIS)
        OVER=SECOND_OWN()-SECND
        CALL MASR1 ('OVERHEAD',OVER)
        DO 99 IS=ISTRA,NSTRAI
          XTIM(IS)=XTIM(IS)+OVER
99      CONTINUE
        WRITE (6,*) 'XTIM(ISTRA)= ',XTIM(ISTRA)
C
        LGLAST=.FALSE.
        LGSTOP=.FALSE.
        DO 100 IPTSI=1,NPTS(ISTRA)
          IF (LGLAST.AND.LGSTOP) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO FURTHER COMP.TIME AVAIL. FOR THIS STRATUM'
            WRITE (6,*) 'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR'
            WRITE (6,*) 'THIS STRATUM'
            CALL MASJ2 ('ISTRA,IPANU=    ',ISTRA,IPANU)
            IF (NPRNLI.GT.0) THEN
              WRITE (6,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
              CALL MASJ1 ('IPRNLS= ',IPRNLS)
            ENDIF
            GOTO 101
          ELSEIF (LGLAST.AND..NOT.LGSTOP) THEN
            CALL LEER(1)
            WRITE (6,*) 'CENSUS ARRAYS FILLED FOR THIS STRATUM'
            WRITE (6,*) 'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR'
            WRITE (6,*) 'THIS STRATUM'
            CALL MASJ2 ('ISTRA,IPANU=    ',ISTRA,IPANU)
            WRITE (6,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
            CALL MASJ1 ('IPRNLS= ',IPRNLS)
            GOTO 101
          ENDIF
          SECND1=SECOND_OWN()
          LGLAST = IPTSI.EQ.NPTS(ISTRA)
          LGLAST = LGLAST.OR.(SECND1.GT.XTIM(ISTRA).AND.IPTSI.GE.2)
          LGSTOP = LGLAST
C  NEXT MONTE CARLO HISTORY
          IF (NLCRR) THEN
C  INITIALIZE RANDOM NUMBERS FOR EACH PARTICLE, TO GENERATE CORRELATION
c slmod begin - temp - not tr
            STOP 'ERROR: ATTEMPTING TO CORRELATE PARTICLES'
c slmod end
            CALL RANSET(ISEED)
            DUMMY=RANF( )
            CALL RANGET(ISEED)
            ISEED=INTMAX-ISEED
            INIV1=0
            INIV2=0
            INIV3=0
            INIV4=0
          ENDIF
          XMCP(ISTRA)=XMCP(ISTRA)+1.
          NPANU=NPANU+1
          IPANU=IPANU+1
          NLEVEL=0
c slmod begin - tr
c...      For now, assume that all particles start out in the first
c         segment, which isn't really true. (Better test if I can 
c         use NLTOR=.TRUE., but I don't really want it.)
          IF (NLTRA.AND..NOT.NLTOR) THEN
c THIS NEEDS TO BE SET WHEN THE PARTICLE LAUNCHES, NTCELL IS PROBABLY SET THERE
c            WRITE(0,*) 'NTRSEG=1' 
            NTRSEG=1
          ENDIF
c slmod end
C  LOCATE TEST PARTICLE: SET INITIAL CO-ORDINATES
          CALL LOCAT1(IPANU)
C  IS BIRTH PROCESS SURVIVED?
          IF (.NOT.LGPART) GOTO 110
          IF (output)
     .    WRITE(6,'(A,3I5,A)') ' MARK: MCARLO: BIRTH SURVIVED  ITYP= ',
     .             ITYP,NPANU,IPANU,' '//partag(ityp)
C
102       CONTINUE
C  FOLLOW NEUTRAL PARTICLE
          IF (ITYP.EQ.1.OR.ITYP.EQ.2) THEN
c slmod begin - tr
            cxcnt = 0

            parcnt(ityp) = parcnt(ityp) + 1

            IF (.FALSE.) THEN
              CALL CLOCK(slntime)
              sldtime = slntime - slltime
              slltime = slntime

              WRITE(0,'(5X,A,I6,F12.6)') 'GO FOLNEUT  NPANU,DTIME =',
     .                                   npanu,sldtime
              WRITE(6,'(5X,A,I6,F12.6)') 'GO FOLNEUT  NPANU,DTIME =',
     .                                   npanu,sldtime
            ENDIF

c slmod end
            CALL FOLNEUT
C  FOLLOW TEST ION
          ELSEIF (ITYP.EQ.3) THEN
            CALL FOLION
          ENDIF
C  NEXT GENERATION ?
c          IF (output)
c     .    WRITE(6,*) 'MARK: MCARLO: NEXT GENERATION? LGPART= ',LGPART
          IF (LGPART) GOTO 102
C
110       CONTINUE
C  NUMBER OF REMAINING NODES AND NUMBER OF LEVELS AT NEXT NODE
          IF (NLEVEL.GT.0) THEN
104         INODES=NODES(NLEVEL)-1
            NODES(NLEVEL)=INODES
            IF(INODES.LE.0) GO TO 103
C  RESTORE VARIABLES AND START NEW TRACK
            DO 105 J=1,NPARTC
              RPST(J)=RSPLST(NLEVEL,J)
105         CONTINUE
            DO 106 J=1,MPARTC
              IPST(J)=ISPLST(NLEVEL,J)
106         CONTINUE
            ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
            NLSRFX=MRSURF.GT.0
            NLSRFY=MPSURF.GT.0
            NLSRFZ=MTSURF.GT.0
            NLSRFA=MASURF.GT.0
            IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,10)
            IF (NLSTOR) CALL STORE(200)
            GOTO 102
C  RETURN TO PREVIOUS LEVEL
103         CONTINUE
            NLEVEL=NLEVEL-1
            IF(NLEVEL.GT.0) GOTO 104
          ENDIF
C  HISTORY HAS ENDED
C
C   MEAN SQUARE
          IF (NSIGI.GT.0) CALL STATS1     (NSBOX,NR1ST,NP2ND,NT3RD,
     .                                     NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
          IF (NSIGI_BGK.GT.0) CALL STATS1_BGK (NSBOX,NR1ST,NP2ND,NT3RD,
     .                                     NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
          IF (NSIGI_COP.GT.0) CALL STATS1_COP (NSBOX,NR1ST,NP2ND,NT3RD,
     .                                     NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))

C  RESET INDEX-ARRAY
          NCLMT = 0
          IMETCL = 0
          LMETSP=.FALSE.
C   COVARIANCES
          IF (NCOVI.GT.0)  CALL COVAR1
C
          IF (TRCTIM) THEN
            SECND2=SECOND_OWN( )
            SECDEL=SECND2-SECND1
            CALL MASJ1R('PART., CPU TIME ',NPANU,SECDEL)
          ENDIF
100     CONTINUE
        CALL LEER(1)
        WRITE (6,*) 'ALL REQUESTED TRAJECTORIES COMPLETED'
        WRITE (6,*) 'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR'
        WRITE (6,*) 'THIS STRATUM'
        CALL MASJ2 ('ISTRA,IPANU=    ',ISTRA,IPANU)
        IF (NPRNLI.GT.0) THEN
          WRITE (6,*) 'M.C. HISTORIES THAT SURVIVED TO CENSUS'
          CALL MASJ1 ('IPRNLS= ',IPRNLS)
        ENDIF
C       GOTO 101
101     CONTINUE
        SECND=SECOND_OWN()
c slmod begin - tr
c
c         Some elementary statisitcs:
c
        WRITE(6,*)
        WRITE(6,'(A,I5)') 'STRATUM = ',istra
        WRITE(6,'(A,I5)') 'NTRAC   = ',npar(istra)

        IF (npar(istra).GT.0) THEN
          WRITE(6,'(A,3F12.5)') 'AVG POS = ',
     .      xpos(istra)/npar(istra),ypos(istra)/npar(istra),
     .      zpos(istra)/npar(istra)
          WRITE(6,'(A,3F12.5)') 'AVG VEL = ',
     .      xvel(istra)/npar(istra),yvel(istra)/npar(istra),
     .      zvel(istra)/npar(istra)
        ENDIF

        WRITE(0,'(A,I5,A)') ' STATUS: Time for stratum = ',INT(secnd),
     .                      ' s'
        IF (output) THEN

          WRITE(0,'(A,5I5)') ' STATUS: PARCNT= ',(parcnt(i),i=1,5)
        ENDIF

c slmod end
C
C**** PARTICLE TRACING FOR THIS STRATUM FINISHED **********************
C
C
C  UPDATE AND CHECK LOGICALS FOR TALLIES
C
      DO 120  IMOL=1,NMOLI
        LOGMOL(0,ISTRA)=LOGMOL(0,ISTRA).OR.LOGMOL(IMOL,ISTRA)
        LOGMOL(IMOL,0)=LOGMOL(IMOL,0).OR.LOGMOL(IMOL,ISTRA)
120   CONTINUE
      DO 130  IATM=1,NATMI
        LOGATM(IATM,0)=LOGATM(IATM,0).OR.LOGATM(IATM,ISTRA)
        LOGATM(0,ISTRA)=LOGATM(0,ISTRA).OR.LOGATM(IATM,ISTRA)
130   CONTINUE
      DO 133  IION=1,NIONI
        LOGION(IION,0)=LOGION(IION,0).OR.LOGION(IION,ISTRA)
        LOGION(0,ISTRA)=LOGION(0,ISTRA).OR.LOGION(IION,ISTRA)
133   CONTINUE
      DO 135  IPLS=1,NPLSI
        LOGPLS(IPLS,0)=LOGPLS(IPLS,0).OR.LOGPLS(IPLS,ISTRA)
        LOGPLS(0,ISTRA)=LOGPLS(0,ISTRA).OR.LOGPLS(IPLS,ISTRA)
135   CONTINUE
      LOGMOL(0,0)=LOGMOL(0,0).OR.LOGMOL(0,ISTRA)
      LOGION(0,0)=LOGION(0,0).OR.LOGION(0,ISTRA)
      LOGATM(0,0)=LOGATM(0,0).OR.LOGATM(0,ISTRA)
      LOGPLS(0,0)=LOGPLS(0,0).OR.LOGPLS(0,ISTRA)
C
C  NUMBER OF LOCATED M.C. HISTORIES FOR THIS STRATUM: XMCP(ISTRA)
C
      IF(XMCP(ISTRA).LT.1.) GOTO 1111
C
      WTT=0.
      DO 200 IATM=1,NATMI
        WTOTA(0,ISTRA)=WTOTA(0,ISTRA)+WTOTA(IATM,ISTRA)
200     WTT=WTT+WTOTA(IATM,ISTRA)*NPRT(IATM)
      DO 201 IMOL=1,NMOLI
        WTOTM(0,ISTRA)=WTOTM(0,ISTRA)+WTOTM(IMOL,ISTRA)
201     WTT=WTT+WTOTM(IMOL,ISTRA)*NPRT(NSPA+IMOL)
      DO 202 IION=1,NIONI
        WTOTI(0,ISTRA)=WTOTI(0,ISTRA)+WTOTI(IION,ISTRA)
202     WTT=WTT+WTOTI(IION,ISTRA)*NPRT(NSPAM+IION)
      DO 203 IPLS=1,NPLSI
        WTOTP(0,ISTRA)=WTOTP(0,ISTRA)+WTOTP(IPLS,ISTRA)
203     WTT=WTT-WTOTP(IPLS,ISTRA)*NPRT(NSPAMI+IPLS)
      CALL LEER(2)
      WRITE (6,*) 'TOTAL WEIGHT OF PRIMARY SOURCE PARTICLES '
      WRITE (6,*) 'BULK IONS, ATOMS, MOLECULES, TEST IONS '
      CALL MASR4 ('WTPLS,WTATM,WTMOL,WTION         ',
     .     WTOTP(0,ISTRA),WTOTA(0,ISTRA),WTOTM(0,ISTRA),WTOTI(0,ISTRA))
      WRITE (6,*) 'TOTAL NUMBER OF MONTE CARLO HISTORIES'
      IMCP=XMCP(ISTRA)
      CALL MASJ1 ('NPART   ',IMCP)
C
C
C  SET SOME SCALING CONSTANTS
C
C  FACTOR FOR FLUXES (AMP) (INPUT FLUX "FLUXT" IS IN AMP)
      FLXFAC(ISTRA)=0.
      IF (SCALV(ISTRA).NE.0.D0) THEN
C  NON DEFAULT SCALING OPTION
        IF (IVLSF(ISTRA).EQ.1) THEN
C  SCALE TO ENFORCE CERTAIN VALUE OF VOLUME TALLY
          IS=ISCLS(ISTRA)
          IT=ISCLT(ISTRA)
          IC=ISCL1(ISTRA)
          IF (ISCL2(ISTRA).GT.0) THEN
            I1=ISCL1(ISTRA)
            I2=ISCL2(ISTRA)
            I3=ISCL3(ISTRA)
            IB=ISCLB(ISTRA)
            IA=ISCLA(ISTRA)
            NBLCKA=NSTRD*(IB-1)+IA
            IC=I1+((I2-1)+(I3-1)*NP2T3)*NR1P2+NBLCKA
          ENDIF
          IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 207
          IF (IS.LT.0.OR.IS.GT.NFSTVI(IT)) GOTO 207
          IF (IC.LT.0.OR.IC.GT.NSBOX) GOTO 207
          IADD=NADDV(IT)*NRAD
          IGFF=NFIRST(IT)
          INDX=IADD+(IC-1)*IGFF+IS
          IF (SCLTAL(1,IT).EQ.1) THEN
            VALUE=ESTIM(INDX)/VOL(IC)/ELCHA
          ELSEIF (SCLTAL(1,IT).EQ.2) THEN
            VALUE=ESTIM(INDX)/ELCHA
          ELSEIF (SCLTAL(1,IT).EQ.3) THEN
            VALUE=ESTIM(INDX)/VOL(IC)
          ELSEIF (SCLTAL(1,IT).EQ.4) THEN
            VALUE=ESTIM(INDX)
          ENDIF
          IF (ABS(VALUE).LE.EPS60) GOTO 207
          VAL=SCALV(ISTRA)
          FLXFAC(ISTRA)=VAL/VALUE
          FLX=FLXFAC(ISTRA)*WTT
          FLUXT(ISTRA)=FLX
        ELSEIF (IVLSF(ISTRA).EQ.2) THEN
C  SCALE TO ENFORCE CERTAIN VALUE OF SURFACE TALLY
          GOTO 207
        ELSE
          GOTO 207
        ENDIF
        GOTO 205
207     WRITE (6,*) 'INCONSISTENT INPUT FOR SCALING OF STRATUM ISTRA '
        WRITE (6,*) 'ISTRA ',ISTRA,IS,IT,IC,VALUE
        WRITE (6,*) 'USE DEFAULT SCALING (FLUX(ISTRA)) '
        FLUXT(ISTRA)=FLUX(ISTRA)
        IF (WTT.NE.0.D0) FLXFAC(ISTRA)=FLUXT(ISTRA)/WTT
205     CONTINUE
      ELSE
C  DEFAULT SCALING OPTION: USE FLUX(ISTRA)
        FLUXT(ISTRA)=FLUX(ISTRA)
        IF (WTT.NE.0.D0) FLXFAC(ISTRA)=FLUXT(ISTRA)/WTT
      ENDIF
C
C  TOTAL TEST PARTICLE FLUX (AMP)
      WRITE (6,*) 'TOTAL SOURCE STRENGTH FOR TEST PARTICLE SPECIES'
      CALL MASR1 ('FLUXT=  ',FLUXT(ISTRA))
      CALL LEER(2)
C
C  ZONE IN-DEPENDENT SCALING FACTORS
      ZWW=FLXFAC(ISTRA)
      ZW=FLXFAC(ISTRA)/ELCHA
C  ZONE DEPENDENT SCALING FACTORS
      DO 206 IC=1,NSBOX
        ZVOLIN(IC)=0.
        ZVOLIW(IC)=0.
        IF (VOL(IC).NE.0.D0) THEN
          ZVOLIN(IC)=ZW /VOL(IC)
          ZVOLIW(IC)=ZWW/VOL(IC)
        ENDIF
206   CONTINUE
      ZVOLNT=ZW /VOLTOT
      ZVOLWT=ZWW/VOLTOT
C
C   STATISTICS , IF REQUIRED
C
      IF (XMCP(ISTRA).LE.1.) GOTO 219
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: STATISTICS'
C
C  FACTORS FOR STANDARD DEVIATION
      ZFLUX=FLXFAC(ISTRA)*XMCP(ISTRA)
      FSIG=SQRT(XMCP(ISTRA)/(XMCP(ISTRA)-1.))
C
      IF (NSIGI.GT.0) THEN
        CALL STATS2(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        DO 210 J=1,NSDVI1
          SDVI(J)=SDVI(J)*100.D0
210     CONTINUE
      ENDIF
      IF (NSIGI_BGK.GT.0) THEN
        CALL STATS2_BGK(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        DO 211 IB=1,NBGVI_STAT
          SGMS_BGK(IB)=SGMS_BGK(IB)*100.D0
          DO 212 J=1,NSBOX
            SIGMA_BGK(IB,J)=SIGMA_BGK(IB,J)*100.D0
212       CONTINUE
211     CONTINUE
      ENDIF
      IF (NSIGI_COP.GT.0) THEN
        CALL STATS2_COP(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        DO 213 IC=1,NCPVI_STAT
          SGMS_COP(IC)=SGMS_COP(IC)*100.D0
          DO 214 J=1,NSBOX
            SIGMA_COP(IC,J)=SIGMA_COP(IC,J)*100.D0
214       CONTINUE
213     CONTINUE
      ENDIF
C
219   CONTINUE
C
C
C*****VOLUME AVERAGED TALLIES  220 - 239
C
C  ATOMIC PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: VOLUME AVERAGED TALLIES'
      DO 220 IATM=1,NATMI
        IF (LOGATM(IATM,ISTRA)) THEN
          DO 221 J=1,NSBOX
            PDENA(IATM,J)=PDENA(IATM,J)*ZVOLIN(J)
            EDENA(IATM,J)=EDENA(IATM,J)*ZVOLIN(J)
            PAAT(IATM,J)=PAAT(IATM,J)*ZVOLIW(J)
            PMAT(IATM,J)=PMAT(IATM,J)*ZVOLIW(J)
            PIAT(IATM,J)=PIAT(IATM,J)*ZVOLIW(J)
            PGENA(IATM,J)=PGENA(IATM,J)*ZVOLIW(J)
            EGENA(IATM,J)=EGENA(IATM,J)*ZVOLIW(J)
            VGENA(IATM,J)=VGENA(IATM,J)*ZVOLIW(J)
221       CONTINUE
        ENDIF
        SCLTAL(IATM,1)=1
        SCLTAL(IATM,4)=1
        SCLTAL(IATM,8)=3
        SCLTAL(IATM,13)=3
        SCLTAL(IATM,18)=3
        SCLTAL(IATM,NTALV-8)=3
        SCLTAL(IATM,NTALV-5)=3
        SCLTAL(IATM,NTALV-2)=3
220   CONTINUE
C
C  MOLECULAR PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO 222 IMOL=1,NMOLI
        IF (LOGMOL(IMOL,ISTRA)) THEN
          DO 223 J=1,NSBOX
            PDENM(IMOL,J)=PDENM(IMOL,J)*ZVOLIN(J)
            EDENM(IMOL,J)=EDENM(IMOL,J)*ZVOLIN(J)
            PAML(IMOL,J)=PAML(IMOL,J)*ZVOLIW(J)
            PMML(IMOL,J)=PMML(IMOL,J)*ZVOLIW(J)
            PIML(IMOL,J)=PIML(IMOL,J)*ZVOLIW(J)
            PGENM(IMOL,J)=PGENM(IMOL,J)*ZVOLIW(J)
            EGENM(IMOL,J)=EGENM(IMOL,J)*ZVOLIW(J)
            VGENM(IMOL,J)=VGENM(IMOL,J)*ZVOLIW(J)
223       CONTINUE
        ENDIF
        SCLTAL(IMOL,2)=1
        SCLTAL(IMOL,5)=1
        SCLTAL(IMOL,9)=3
        SCLTAL(IMOL,14)=3
        SCLTAL(IMOL,19)=3
        SCLTAL(IMOL,NTALV-7)=3
        SCLTAL(IMOL,NTALV-4)=3
        SCLTAL(IMOL,NTALV-1)=3
222   CONTINUE
C
C  TEST ION PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO 225 IION=1,NIONI
        IF (LOGION(IION,ISTRA)) THEN
          DO 226 J=1,NSBOX
            PDENI(IION,J)=PDENI(IION,J)*ZVOLIN(J)
            EDENI(IION,J)=EDENI(IION,J)*ZVOLIN(J)
            PAIO(IION,J)=PAIO(IION,J)*ZVOLIW(J)
            PMIO(IION,J)=PMIO(IION,J)*ZVOLIW(J)
            PIIO(IION,J)=PIIO(IION,J)*ZVOLIW(J)
            PGENI(IION,J)=PGENI(IION,J)*ZVOLIW(J)
            EGENI(IION,J)=EGENI(IION,J)*ZVOLIW(J)
            VGENI(IION,J)=VGENI(IION,J)*ZVOLIW(J)
226       CONTINUE
        ENDIF
        SCLTAL(IION,3)=1
        SCLTAL(IION,6)=1
        SCLTAL(IION,10)=3
        SCLTAL(IION,15)=3
        SCLTAL(IION,20)=3
        SCLTAL(IION,NTALV-6)=3
        SCLTAL(IION,NTALV-3)=3
        SCLTAL(IION,NTALV-0)=3
225   CONTINUE
C
C  BULK ION PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO 227 IPLS=1,NPLSI
        IF (LOGPLS(IPLS,ISTRA)) THEN
          DO 228 J=1,NSBOX
            PAPL(IPLS,J)=PAPL(IPLS,J)*ZVOLIW(J)
            PMPL(IPLS,J)=PMPL(IPLS,J)*ZVOLIW(J)
            PIPL(IPLS,J)=PIPL(IPLS,J)*ZVOLIW(J)
228       CONTINUE
        SCLTAL(IPLS,11)=3
        SCLTAL(IPLS,16)=3
        SCLTAL(IPLS,21)=3
        ENDIF
227   CONTINUE
C
C  ADDITIONAL TRACKLENGTH ESTIMATED TALLIES FOR THE STRATUM ISTRA
C  TALLY NO. NTALA
C
      DO 230 IADV=1,NADVI
        IF (IADVE(IADV).EQ.1) THEN
C  SCALE # PER VOLUME
          DO 231 J=1,NSBOX
            ADDV(IADV,J)=ADDV(IADV,J)*ZVOLIN(J)
231       CONTINUE
          SCLTAL(IADV,NTALA)=1
        ELSEIF (IADVE(IADV).EQ.2) THEN
C  SCALE # PER CELL
          DO 232 J=1,NSBOX
            ADDV(IADV,J)=ADDV(IADV,J)*ZW
232       CONTINUE
          SCLTAL(IADV,NTALA)=2
C  SCALE AMP/S PER VOLUME
        ELSEIF (IADVE(IADV).EQ.3) THEN
          DO 233 J=1,NSBOX
            ADDV(IADV,J)=ADDV(IADV,J)*ZVOLIW(J)
233       CONTINUE
          SCLTAL(IADV,NTALA)=3
C  SCALE AMP/S PER CELL
        ELSEIF (IADVE(IADV).EQ.4) THEN
          DO 234 J=1,NSBOX
            ADDV(IADV,J)=ADDV(IADV,J)*ZWW
234       CONTINUE
          SCLTAL(IADV,NTALA)=4
C       ELSE
C  DON'T SCALE AT ALL
          SCLTAL(IADV,NTALA)=0
        ENDIF
230   CONTINUE
C
C  ADDITIONAL COLLISION ESTIMATED TALLIES FOR THE STRATUM ISTRA
C  TALLY NO. NTALC
C
      DO 235 ICLV=1,NCLVI
        IF (ICLVE(ICLV).EQ.1) THEN
C  SCALE # PER VOLUME
          DO 236 J=1,NSBOX
            COLV(ICLV,J)=COLV(ICLV,J)*ZVOLIN(J)
236       CONTINUE
          SCLTAL(ICLV,NTALC)=1
        ELSEIF (ICLVE(ICLV).EQ.2) THEN
C  SCALE # PER CELL
          DO 237 J=1,NSBOX
            COLV(ICLV,J)=COLV(ICLV,J)*ZW
237       CONTINUE
          SCLTAL(ICLV,NTALC)=2
        ELSEIF (ICLVE(ICLV).EQ.3) THEN
C  SCALE AMP/S PER VOLUME
          DO 238 J=1,NSBOX
            COLV(ICLV,J)=COLV(ICLV,J)*ZVOLIW(J)
238       CONTINUE
          SCLTAL(ICLV,NTALC)=3
        ELSEIF (ICLVE(ICLV).EQ.4) THEN
C  SCALE AMP/S PER CELL
          DO 239 J=1,NSBOX
            COLV(ICLV,J)=COLV(ICLV,J)*ZWW
239       CONTINUE
          SCLTAL(ICLV,NTALC)=4
C       ELSE
C  DON'T SCALE AT ALL
          SCLTAL(ICLV,NTALC)=0
        ENDIF
235   CONTINUE
C
C  ADDITIONAL SNAPSHOT ESTIMATED TALLIES FOR THE STRATUM ISTRA
C  TALLY NO. NTALT, FIRST: SNAPV=SNAPV*DTIMV, THEN: SCALING
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: SNAP SHOT TALLIES'
      DO 245 ISNV=1,NSNVI
        IF (ISNVE(ISNV).EQ.1) THEN
C  SCALE # PER VOLUME
          DO 246 J=1,NSBOX
            SNAPV(ISNV,J)=SNAPV(ISNV,J)*DTIMV*ZVOLIN(J)
246       CONTINUE
          SCLTAL(ISNV,NTALT)=1
        ELSEIF (ISNVE(ISNV).EQ.2) THEN
C  SCALE # PER CELL
          DO 247 J=1,NSBOX
            SNAPV(ISNV,J)=SNAPV(ISNV,J)*DTIMV*ZW
247       CONTINUE
          SCLTAL(ISNV,NTALT)=2
        ELSEIF (ISNVE(ISNV).EQ.3) THEN
C  SCALE AMP/S PER VOLUME
          DO 248 J=1,NSBOX
            SNAPV(ISNV,J)=SNAPV(ISNV,J)*DTIMV*ZVOLIW(J)
248       CONTINUE
          SCLTAL(ISNV,NTALT)=3
        ELSEIF (ISNVE(ISNV).EQ.4) THEN
C  SCALE AMP/S PER CELL
          DO 249 J=1,NSBOX
            SNAPV(ISNV,J)=SNAPV(ISNV,J)*DTIMV*ZWW
249       CONTINUE
          SCLTAL(ISNV,NTALT)=4
C       ELSE
C  DON'T SCALE AT ALL
          DO J=1,NSBOX
            SNAPV(ISNV,J)=SNAPV(ISNV,J)*DTIMV
          ENDDO
          SCLTAL(ISNV,NTALT)=0
        ENDIF
245   CONTINUE
C
C  TALLIES FOR FLUID CODE COUPLING FOR THE STRATUM ISTRA
C  TALLY NO. NTALM
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: FLUID CODE TALLIES'
      DO 255 ICPV=1,NCPVI
        IF (ICPVE(ICPV).EQ.1) THEN
C  SCALE # PER VOLUME
          DO 256 J=1,NSBOX
            COPV(ICPV,J)=COPV(ICPV,J)*ZVOLIN(J)
256       CONTINUE
c slmod begin - tr
          IF (.TRUE..OR.ISTRA.EQ.NSTRAI) THEN
c            WRITE(0,*) 'MARK: SCALING # PER VOLUME ',icpv
            DO J=1,NSBOX            
              DO I1 = 1, NMOMCHA
                COPV2(ICPV,J,I1)=COPV2(ICPV,J,I1)*ZVOLIN(J)
              ENDDO
            ENDDO
          ENDIF
c slmod end
          SCLTAL(ICPV,NTALM)=1
        ELSEIF (ICPVE(ICPV).EQ.2) THEN
C  SCALE # PER CELL
c...sltmp
          STOP 'STOP: COPV BUSINESS 10'
          DO 257 J=1,NSBOX
            COPV(ICPV,J)=COPV(ICPV,J)*ZW
257       CONTINUE
          SCLTAL(ICPV,NTALM)=2
        ELSEIF (ICPVE(ICPV).EQ.3) THEN
C  SCALE AMP/S PER VOLUME
          DO 258 J=1,NSBOX
            COPV(ICPV,J)=COPV(ICPV,J)*ZVOLIW(J)
258       CONTINUE
c slmod begin - tr
          DO J=1,NSBOX
            DO I1 = 1, NMOMCHA
              COPV2(ICPV,J,I1)=COPV2(ICPV,J,I1)*ZVOLIW(J)
            ENDDO
          ENDDO
c slmod end
          SCLTAL(ICPV,NTALM)=3
        ELSEIF (ICPVE(ICPV).EQ.4) THEN
C  SCALE AMP/S PER CELL
          DO 259 J=1,NSBOX
c...sltmp
            STOP 'STOP: COPV BUSINESS 09'
            COPV(ICPV,J)=COPV(ICPV,J)*ZWW
259       CONTINUE
          SCLTAL(ICPV,NTALM)=4
C       ELSE
C  DON'T SCALE AT ALL
          SCLTAL(ICPV,NTALM)=0
        ENDIF
255   CONTINUE
C
C  TALLIES FOR FLUID CODE COUPLING FOR THE STRATUM ISTRA
C  TALLY NO. NTALB
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: NTALB TALLIES'
      DO 265 IBGV=1,NBGVI
        IF (IBGVE(IBGV).EQ.1) THEN
C  SCALE # PER VOLUME
          DO 266 J=1,NSBOX
            BGKV(IBGV,J)=BGKV(IBGV,J)*ZVOLIN(J)
266       CONTINUE
          SCLTAL(IBGV,NTALB)=1
        ELSEIF (IBGVE(IBGV).EQ.2) THEN
C  SCALE # PER CELL
          DO 267 J=1,NSBOX
            BGKV(IBGV,J)=BGKV(IBGV,J)*ZW
267       CONTINUE
          SCLTAL(IBGV,NTALB)=2
        ELSEIF (IBGVE(IBGV).EQ.3) THEN
C  SCALE AMP/S PER VOLUME
          DO 268 J=1,NSBOX
            BGKV(IBGV,J)=BGKV(IBGV,J)*ZVOLIW(J)
268       CONTINUE
          SCLTAL(IBGV,NTALB)=3
        ELSEIF (IBGVE(IBGV).EQ.4) THEN
C  SCALE AMP/S PER CELL
          DO 269 J=1,NSBOX
            BGKV(IBGV,J)=BGKV(IBGV,J)*ZWW
269       CONTINUE
          SCLTAL(IBGV,NTALB)=4
C       ELSE
C  DON'T SCALE AT ALL
          SCLTAL(IBGV,NTALB)=0
        ENDIF
265   CONTINUE
C
C  OTHER TALLIES ESTIMATED FROM HISTORIES
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: OTHER TALLIES'
      DO 270 J=1,NSBOX
        PAEL(J)=PAEL(J)*ZVOLIW(J)
        EAEL(J)=EAEL(J)*ZVOLIW(J)
        EAAT(J)=EAAT(J)*ZVOLIW(J)
        EAML(J)=EAML(J)*ZVOLIW(J)
        EAPL(J)=EAPL(J)*ZVOLIW(J)
        EAIO(J)=EAIO(J)*ZVOLIW(J)
C
        PMEL(J)=PMEL(J)*ZVOLIW(J)
        EMEL(J)=EMEL(J)*ZVOLIW(J)
        EMAT(J)=EMAT(J)*ZVOLIW(J)
        EMML(J)=EMML(J)*ZVOLIW(J)
        EMIO(J)=EMIO(J)*ZVOLIW(J)
        EMPL(J)=EMPL(J)*ZVOLIW(J)
C
        PIEL(J)=PIEL(J)*ZVOLIW(J)
        EIEL(J)=EIEL(J)*ZVOLIW(J)
        EIAT(J)=EIAT(J)*ZVOLIW(J)
        EIML(J)=EIML(J)*ZVOLIW(J)
        EIIO(J)=EIIO(J)*ZVOLIW(J)
        EIPL(J)=EIPL(J)*ZVOLIW(J)
270   CONTINUE
      SCLTAL(1,7)=3
      SCLTAL(1,12)=3
      SCLTAL(1,17)=3
      SCLTAL(1,22)=3
      SCLTAL(1,23)=3
      SCLTAL(1,24)=3
      SCLTAL(1,25)=3
      SCLTAL(1,26)=3
      SCLTAL(1,27)=3
      SCLTAL(1,28)=3
      SCLTAL(1,29)=3
      SCLTAL(1,30)=3
      SCLTAL(1,31)=3
      SCLTAL(1,32)=3
      SCLTAL(1,33)=3
      SCLTAL(1,34)=3
      SCLTAL(1,35)=3
      SCLTAL(1,36)=3
C
C   REPLACE DEFAULT TALLIES BY USER SUPPLIED
C   COLLISION ESTIMATED TALLIES
C   THIS IS DONE BEFORE VOLUME INTEGRATION THUS THERE IS THE RISK TO
C   DESTROY TERMS NEEDED FOR GLOBAL BALANCES
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: USER SUPPLIED TALLIES'
      DO 285 ICLV=1,NCLVI
        IS=ICLVS(ICLV)
        IT=ICLVT(ICLV)
        IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 285
        IGFFT=NFSTVI(IT)
        IGFF=NFIRST(IT)
        IF (IS.LE.0.OR.IS.GT.IGFFT) GOTO 285
        IADD=NADDV(IT)*NRAD
        DO 286 J=1,NSBOX
          INDX=IADD+(J-1)*IGFF+IS
          ESTIM(INDX)=COLV(ICLV,J)
286     CONTINUE
285   CONTINUE
C
C   REPLACE DEFAULT TALLIES BY USER SUPPLIED
C   TRACKLENGTH ESTIMATED TALLIES
C   THIS IS DONE BEFORE VOLUME INTEGRATION THUS THERE IS THE RISK TO
C   DESTROY TERMS NEEDED FOR GLOBAL BALANCES
C
      DO 290 IADV=1,NADVI
        IS=IADVS(IADV)
        IT=IADVT(IADV)
        IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 290
        IGFFT=NFSTVI(IT)
        IGFF=NFIRST(IT)
        IF (IS.LE.0.OR.IS.GT.IGFFT) GOTO 290
        IADD=NADDV(IT)*NRAD
        DO 295 J=1,NSBOX
          INDX=IADD+(J-1)*IGFF+IS
          ESTIM(INDX)=ADDV(IADV,J)
295     CONTINUE
290   CONTINUE
C
C   INTEGRATE VOLUME AVERAGED PROFILES   450 --- 459
C
      IF (output)
     .WRITE(0,*) 'MARK: MCARLO: INTEGRATING PROFILES'
      CALL INTTAL (PAEL,VOL,1,1,NSBOX,PAELI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (PMEL,VOL,1,1,NSBOX,PMELI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (PIEL,VOL,1,1,NSBOX,PIELI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EAEL,VOL,1,1,NSBOX,EAELI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EAAT,VOL,1,1,NSBOX,EAATI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EAML,VOL,1,1,NSBOX,EAMLI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EAIO,VOL,1,1,NSBOX,EAIOI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EAPL,VOL,1,1,NSBOX,EAPLI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EMEL,VOL,1,1,NSBOX,EMELI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EMAT,VOL,1,1,NSBOX,EMATI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EMML,VOL,1,1,NSBOX,EMMLI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EMIO,VOL,1,1,NSBOX,EMIOI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EMPL,VOL,1,1,NSBOX,EMPLI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EIEL,VOL,1,1,NSBOX,EIELI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EIAT,VOL,1,1,NSBOX,EIATI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EIML,VOL,1,1,NSBOX,EIMLI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EIIO,VOL,1,1,NSBOX,EIIOI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      CALL INTTAL (EIPL,VOL,1,1,NSBOX,EIPLI(ISTRA),NR1ST,NP2ND,NT3RD,
     .                                             NBMLT)
      DO 450 IATM=1,NATMI
        IF (.NOT.LOGATM(IATM,ISTRA)) GOTO 450
          CALL INTTAL (PDENA,VOL,IATM,NATM,NSBOX,PDENAI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (EDENA,VOL,IATM,NATM,NSBOX,EDENAI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PAAT,VOL,IATM,NATM,NSBOX,PAATI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PMAT,VOL,IATM,NATM,NSBOX,PMATI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PIAT,VOL,IATM,NATM,NSBOX,PIATI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PGENA,VOL,IATM,NATM,NSBOX,PGENAI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (EGENA,VOL,IATM,NATM,NSBOX,EGENAI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (VGENA,VOL,IATM,NATM,NSBOX,VGENAI(IATM,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
450   CONTINUE
      DO 451 IMOL=1,NMOLI
        IF (.NOT.LOGMOL(IMOL,ISTRA)) GOTO 451
          CALL INTTAL (PDENM,VOL,IMOL,NMOL,NSBOX,PDENMI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (EDENM,VOL,IMOL,NMOL,NSBOX,EDENMI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PAML,VOL,IMOL,NMOL,NSBOX,PAMLI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PMML,VOL,IMOL,NMOL,NSBOX,PMMLI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PIML,VOL,IMOL,NMOL,NSBOX,PIMLI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PGENM,VOL,IMOL,NMOL,NSBOX,PGENMI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (EGENM,VOL,IMOL,NMOL,NSBOX,EGENMI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (VGENM,VOL,IMOL,NMOL,NSBOX,VGENMI(IMOL,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
451   CONTINUE
      DO 452 IION=1,NIONI
        IF (.NOT.LOGION(IION,ISTRA)) GOTO 452
          CALL INTTAL (PDENI,VOL,IION,NION,NSBOX,PDENII(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (EDENI,VOL,IION,NION,NSBOX,EDENII(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PAIO,VOL,IION,NION,NSBOX,PAIOI(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PMIO,VOL,IION,NION,NSBOX,PMIOI(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PIIO,VOL,IION,NION,NSBOX,PIIOI(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (PGENI,VOL,IION,NION,NSBOX,PGENII(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (EGENI,VOL,IION,NION,NSBOX,EGENII(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
          CALL INTTAL (VGENI,VOL,IION,NION,NSBOX,VGENII(IION,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
452   CONTINUE
      DO 454 IPLS=1,NPLSI
        CALL INTTAL (PAPL,VOL,IPLS,NPLS,NSBOX,PAPLI(IPLS,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (PMPL,VOL,IPLS,NPLS,NSBOX,PMPLI(IPLS,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (PIPL,VOL,IPLS,NPLS,NSBOX,PIPLI(IPLS,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
454   CONTINUE
      DO 455 IADV=1,NADVI
        IF (IADVE(IADV).NE.2.AND.IADVE(IADV).NE.4)
     .  CALL INTTAL (ADDV,VOL,IADV,NADV,NSBOX,ADDVI(IADV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        IF (IADVE(IADV).EQ.2.OR.IADVE(IADV).EQ.4)
     .  CALL INTVOL (ADDV,     IADV,NADV,NSBOX,ADDVI(IADV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
455   CONTINUE
      DO 456 ICLV=1,NCLVI
        IF (ICLVE(ICLV).NE.2.AND.ICLVE(ICLV).NE.4)
     .  CALL INTTAL (COLV,VOL,ICLV,NCLV,NSBOX,COLVI(ICLV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        IF (ICLVE(ICLV).EQ.2.OR.ICLVE(ICLV).EQ.4)
     .  CALL INTVOL (COLV,     ICLV,NCLV,NSBOX,COLVI(ICLV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
456   CONTINUE
      DO 457 ISNV=1,NSNVI
        IF (ISNVE(ISNV).NE.2.AND.ISNVE(ISNV).NE.4)
     .  CALL INTTAL (SNAPV,VOL,ISNV,NSNV,NSBOX,SNAPVI(ISNV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        IF (ISNVE(ISNV).EQ.2.OR.ISNVE(ISNV).EQ.4)
     .  CALL INTVOL (SNAPV,    ISNV,NSNV,NSBOX,SNAPVI(ISNV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
457   CONTINUE
      DO 458 ICPV=1,NCPVI
        CALL INTTAL (COPV,VOL,ICPV,NCPV,NSBOX,COPVI(ICPV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
c slmod begin - tr
        DO I1 = 1, NMOMCHA
          CALL INTTAL (COPV2(1,1,I1),
     .                 VOL,ICPV,NCPV,NSBOX,RDUM1,
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
        ENDDO
c slmod end
458   CONTINUE
      DO 459 IBGV=1,NBGVI
        CALL INTTAL (BGKV,VOL,IBGV,NBGV,NSBOX,BGKVI(IBGV,ISTRA),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
459   CONTINUE
C
C   SYMMETRISE VOLUME AVERAGED TALLIES?
      IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
        CALL SYMET(ESTIM,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .             NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
      ENDIF
C
C  WORK WITH VOLUME AVERAGED TALLIES FOR THIS STRATUM FINISHED
C
C  SCALE SURFACE AVERAGED ESTIMATORS AND OTHER FLUXES 600 - 630
C
      DO 600 J=NESTM1+1,NESTM1+NESTM2
600     ESTIM(J)=ESTIM(J)*FLXFAC(ISTRA)
C
      EPATI(ISTRA)=EPATI(ISTRA)*FLXFAC(ISTRA)
      EPMLI(ISTRA)=EPMLI(ISTRA)*FLXFAC(ISTRA)
      EPIOI(ISTRA)=EPIOI(ISTRA)*FLXFAC(ISTRA)
C
      ETOTA(ISTRA)=ETOTA(ISTRA)*FLXFAC(ISTRA)
      ETOTM(ISTRA)=ETOTM(ISTRA)*FLXFAC(ISTRA)
      ETOTI(ISTRA)=ETOTI(ISTRA)*FLXFAC(ISTRA)
      ETOTP(ISTRA)=ETOTP(ISTRA)*FLXFAC(ISTRA)
      PTRASH(ISTRA)=PTRASH(ISTRA)*FLXFAC(ISTRA)
      ETRASH(ISTRA)=ETRASH(ISTRA)*FLXFAC(ISTRA)
      DO 610 IATM=0,NATM
        PPATI(IATM,ISTRA)=PPATI(IATM,ISTRA)*FLXFAC(ISTRA)
610     WTOTA(IATM,ISTRA)=WTOTA(IATM,ISTRA)*FLXFAC(ISTRA)
      DO 615 IMOL=0,NMOL
        PPMLI(IMOL,ISTRA)=PPMLI(IMOL,ISTRA)*FLXFAC(ISTRA)
615     WTOTM(IMOL,ISTRA)=WTOTM(IMOL,ISTRA)*FLXFAC(ISTRA)
      DO 620 IION=0,NION
        EELFI(IION,ISTRA)=EELFI(IION,ISTRA)*FLXFAC(ISTRA)
        PPIOI(IION,ISTRA)=PPIOI(IION,ISTRA)*FLXFAC(ISTRA)
620     WTOTI(IION,ISTRA)=WTOTI(IION,ISTRA)*FLXFAC(ISTRA)
      DO 625 IPLS=0,NPLS
625     WTOTP(IPLS,ISTRA)=WTOTP(IPLS,ISTRA)*FLXFAC(ISTRA)
C
C   SUM OVER SURFACE INDEX
C   IN THE SURFACE AVERAGED ESTIMATORS
C
      DO 632 IMOL=1,NMOLI
        IF (.NOT.LOGMOL(IMOL,ISTRA)) GOTO 632
        DO 633 J=1,NLIMPS
          IF (ILIIN(J).LE.0) GOTO 633
          PRFAMI(IMOL,ISTRA)=PRFAMI(IMOL,ISTRA)+PRFAML(IMOL,J)
          PRFMMI(IMOL,ISTRA)=PRFMMI(IMOL,ISTRA)+PRFMML(IMOL,J)
          PRFIMI(IMOL,ISTRA)=PRFIMI(IMOL,ISTRA)+PRFIML(IMOL,J)
          PRFPMI(IMOL,ISTRA)=PRFPMI(IMOL,ISTRA)+PRFPML(IMOL,J)
          POTMLI(IMOL,ISTRA)=POTMLI(IMOL,ISTRA)-POTML(IMOL,J)
          ERFAMI(IMOL,ISTRA)=ERFAMI(IMOL,ISTRA)+ERFAML(IMOL,J)
          ERFMMI(IMOL,ISTRA)=ERFMMI(IMOL,ISTRA)+ERFMML(IMOL,J)
          ERFIMI(IMOL,ISTRA)=ERFIMI(IMOL,ISTRA)+ERFIML(IMOL,J)
          ERFPMI(IMOL,ISTRA)=ERFPMI(IMOL,ISTRA)+ERFPML(IMOL,J)
          EOTMLI(IMOL,ISTRA)=EOTMLI(IMOL,ISTRA)-EOTML(IMOL,J)
          SPTMLI(IMOL,ISTRA)=SPTMLI(IMOL,ISTRA)+SPTML(IMOL,J)
633     CONTINUE
632   CONTINUE
C
      DO 635 IATM=1,NATMI
        IF (.NOT.LOGATM(IATM,ISTRA)) GOTO 635
        DO 636 J=1,NLIMPS
          IF (ILIIN(J).LE.0) GOTO 636
          PRFAAI(IATM,ISTRA)=PRFAAI(IATM,ISTRA)+PRFAAT(IATM,J)
          PRFMAI(IATM,ISTRA)=PRFMAI(IATM,ISTRA)+PRFMAT(IATM,J)
          PRFIAI(IATM,ISTRA)=PRFIAI(IATM,ISTRA)+PRFIAT(IATM,J)
          PRFPAI(IATM,ISTRA)=PRFPAI(IATM,ISTRA)+PRFPAT(IATM,J)
          POTATI(IATM,ISTRA)=POTATI(IATM,ISTRA)-POTAT(IATM,J)
          ERFAAI(IATM,ISTRA)=ERFAAI(IATM,ISTRA)+ERFAAT(IATM,J)
          ERFMAI(IATM,ISTRA)=ERFMAI(IATM,ISTRA)+ERFMAT(IATM,J)
          ERFIAI(IATM,ISTRA)=ERFIAI(IATM,ISTRA)+ERFIAT(IATM,J)
          ERFPAI(IATM,ISTRA)=ERFPAI(IATM,ISTRA)+ERFPAT(IATM,J)
          EOTATI(IATM,ISTRA)=EOTATI(IATM,ISTRA)-EOTAT(IATM,J)
          SPTATI(IATM,ISTRA)=SPTATI(IATM,ISTRA)+SPTAT(IATM,J)
636     CONTINUE
635   CONTINUE
C
      DO 638 IION=1,NIONI
        IF (.NOT.LOGION(IION,ISTRA)) GOTO 638
        DO 639 J=1,NLIMPS
          IF (ILIIN(J).LE.0) GOTO 639
          PRFAII(IION,ISTRA)=PRFAII(IION,ISTRA)+PRFAIO(IION,J)
          PRFMII(IION,ISTRA)=PRFMII(IION,ISTRA)+PRFMIO(IION,J)
          PRFIII(IION,ISTRA)=PRFIII(IION,ISTRA)+PRFIIO(IION,J)
          PRFPII(IION,ISTRA)=PRFPII(IION,ISTRA)+PRFPIO(IION,J)
          POTIOI(IION,ISTRA)=POTIOI(IION,ISTRA)-POTIO(IION,J)
          ERFAII(IION,ISTRA)=ERFAII(IION,ISTRA)+ERFAIO(IION,J)
          ERFMII(IION,ISTRA)=ERFMII(IION,ISTRA)+ERFMIO(IION,J)
          ERFIII(IION,ISTRA)=ERFIII(IION,ISTRA)+ERFIIO(IION,J)
          ERFPII(IION,ISTRA)=ERFPII(IION,ISTRA)+ERFPIO(IION,J)
          EOTIOI(IION,ISTRA)=EOTIOI(IION,ISTRA)-EOTIO(IION,J)
          SPTIOI(IION,ISTRA)=SPTIOI(IION,ISTRA)+SPTIO(IION,J)
639     CONTINUE
638   CONTINUE
C
      DO 641 IPLS=1,NPLSI
        IF (.NOT.LOGPLS(IPLS,ISTRA)) GOTO 641
        DO 642 J=1,NLIMPS
          IF (ILIIN(J).LE.0) GOTO 642
          POTPLI(IPLS,ISTRA)=POTPLI(IPLS,ISTRA)-POTPL(IPLS,J)
          EOTPLI(IPLS,ISTRA)=EOTPLI(IPLS,ISTRA)-EOTPL(IPLS,J)
          SPTPLI(IPLS,ISTRA)=SPTPLI(IPLS,ISTRA)+SPTPL(IPLS,J)
642     CONTINUE
641   CONTINUE
C
      DO 651 ISPZ=1,NSPTOT
        DO 652 J=1,NLIMPS
          IF (ILIIN(J).LE.0) GOTO 652
            SPUMPI(ISPZ,ISTRA)=SPUMPI(ISPZ,ISTRA)+SPUMP(ISPZ,J)
652     CONTINUE
651   CONTINUE
C
C  SUM OVER SPECIES INDEX FOR INTEGRATED VOLUME AVERAGED TALLIES
C                         AND INTEGRATED SURFACE AVERAGED TALLIES
C
      DO 661 IPLS=1,NPLSI
        PAPLI(0,ISTRA)=PAPLI(0,ISTRA)+PAPLI(IPLS,ISTRA)
        PMPLI(0,ISTRA)=PMPLI(0,ISTRA)+PMPLI(IPLS,ISTRA)
        PIPLI(0,ISTRA)=PIPLI(0,ISTRA)+PIPLI(IPLS,ISTRA)
        POTPLI(0,ISTRA)=POTPLI(0,ISTRA)+POTPLI(IPLS,ISTRA)
        EOTPLI(0,ISTRA)=EOTPLI(0,ISTRA)+EOTPLI(IPLS,ISTRA)
        SPTPLI(0,ISTRA)=SPTPLI(0,ISTRA)+SPTPLI(IPLS,ISTRA)
661   CONTINUE
      DO 662 IION=1,NIONI
        PDENII(0,ISTRA)=PDENII(0,ISTRA)+PDENII(IION,ISTRA)
        EDENII(0,ISTRA)=EDENII(0,ISTRA)+EDENII(IION,ISTRA)
        PPIOI (0,ISTRA)=PPIOI (0,ISTRA)+PPIOI (IION,ISTRA)
        PAIOI (0,ISTRA)=PAIOI (0,ISTRA)+PAIOI (IION,ISTRA)
        PMIOI (0,ISTRA)=PMIOI (0,ISTRA)+PMIOI (IION,ISTRA)
        PIIOI (0,ISTRA)=PIIOI (0,ISTRA)+PIIOI (IION,ISTRA)
        POTIOI(0,ISTRA)=POTIOI(0,ISTRA)+POTIOI(IION,ISTRA)
        PRFAII(0,ISTRA)=PRFAII(0,ISTRA)+PRFAII(IION,ISTRA)
        PRFMII(0,ISTRA)=PRFMII(0,ISTRA)+PRFMII(IION,ISTRA)
        PRFIII(0,ISTRA)=PRFIII(0,ISTRA)+PRFIII(IION,ISTRA)
        PRFPII(0,ISTRA)=PRFPII(0,ISTRA)+PRFPII(IION,ISTRA)
        EOTIOI(0,ISTRA)=EOTIOI(0,ISTRA)+EOTIOI(IION,ISTRA)
        ERFAII(0,ISTRA)=ERFAII(0,ISTRA)+ERFAII(IION,ISTRA)
        ERFMII(0,ISTRA)=ERFMII(0,ISTRA)+ERFMII(IION,ISTRA)
        ERFIII(0,ISTRA)=ERFIII(0,ISTRA)+ERFIII(IION,ISTRA)
        ERFPII(0,ISTRA)=ERFPII(0,ISTRA)+ERFPII(IION,ISTRA)
        SPTIOI(0,ISTRA)=SPTIOI(0,ISTRA)+SPTIOI(IION,ISTRA)
        PGENII(0,ISTRA)=PGENII(0,ISTRA)+PGENII(IION,ISTRA)
        EGENII(0,ISTRA)=EGENII(0,ISTRA)+EGENII(IION,ISTRA)
        VGENII(0,ISTRA)=VGENII(0,ISTRA)+VGENII(IION,ISTRA)
        EELFI (0,ISTRA)=EELFI (0,ISTRA)+EELFI (IION,ISTRA)
662   CONTINUE
      DO 663 IMOL=1,NMOLI
        PDENMI(0,ISTRA)=PDENMI(0,ISTRA)+PDENMI(IMOL,ISTRA)
        EDENMI(0,ISTRA)=EDENMI(0,ISTRA)+EDENMI(IMOL,ISTRA)
        PPMLI (0,ISTRA)=PPMLI (0,ISTRA)+PPMLI (IMOL,ISTRA)
        PAMLI (0,ISTRA)=PAMLI (0,ISTRA)+PAMLI (IMOL,ISTRA)
        PMMLI (0,ISTRA)=PMMLI (0,ISTRA)+PMMLI (IMOL,ISTRA)
        PIMLI (0,ISTRA)=PIMLI (0,ISTRA)+PIMLI (IMOL,ISTRA)
        POTMLI(0,ISTRA)=POTMLI(0,ISTRA)+POTMLI(IMOL,ISTRA)
        PRFAMI(0,ISTRA)=PRFAMI(0,ISTRA)+PRFAMI(IMOL,ISTRA)
        PRFMMI(0,ISTRA)=PRFMMI(0,ISTRA)+PRFMMI(IMOL,ISTRA)
        PRFIMI(0,ISTRA)=PRFIMI(0,ISTRA)+PRFIMI(IMOL,ISTRA)
        PRFPMI(0,ISTRA)=PRFPMI(0,ISTRA)+PRFPMI(IMOL,ISTRA)
        EOTMLI(0,ISTRA)=EOTMLI(0,ISTRA)+EOTMLI(IMOL,ISTRA)
        ERFAMI(0,ISTRA)=ERFAMI(0,ISTRA)+ERFAMI(IMOL,ISTRA)
        ERFMMI(0,ISTRA)=ERFMMI(0,ISTRA)+ERFMMI(IMOL,ISTRA)
        ERFIMI(0,ISTRA)=ERFIMI(0,ISTRA)+ERFIMI(IMOL,ISTRA)
        ERFPMI(0,ISTRA)=ERFPMI(0,ISTRA)+ERFPMI(IMOL,ISTRA)
        SPTMLI(0,ISTRA)=SPTMLI(0,ISTRA)+SPTMLI(IMOL,ISTRA)
        PGENMI(0,ISTRA)=PGENMI(0,ISTRA)+PGENMI(IMOL,ISTRA)
        EGENMI(0,ISTRA)=EGENMI(0,ISTRA)+EGENMI(IMOL,ISTRA)
        VGENMI(0,ISTRA)=VGENMI(0,ISTRA)+VGENMI(IMOL,ISTRA)
663   CONTINUE
      DO 664 IATM=1,NATMI
        PDENAI(0,ISTRA)=PDENAI(0,ISTRA)+PDENAI(IATM,ISTRA)
        EDENAI(0,ISTRA)=EDENAI(0,ISTRA)+EDENAI(IATM,ISTRA)
        PPATI (0,ISTRA)=PPATI (0,ISTRA)+PPATI (IATM,ISTRA)
        PAATI (0,ISTRA)=PAATI (0,ISTRA)+PAATI (IATM,ISTRA)
        PMATI (0,ISTRA)=PMATI (0,ISTRA)+PMATI (IATM,ISTRA)
        PIATI (0,ISTRA)=PIATI (0,ISTRA)+PIATI (IATM,ISTRA)
        POTATI(0,ISTRA)=POTATI(0,ISTRA)+POTATI(IATM,ISTRA)
        PRFAAI(0,ISTRA)=PRFAAI(0,ISTRA)+PRFAAI(IATM,ISTRA)
        PRFMAI(0,ISTRA)=PRFMAI(0,ISTRA)+PRFMAI(IATM,ISTRA)
        PRFIAI(0,ISTRA)=PRFIAI(0,ISTRA)+PRFIAI(IATM,ISTRA)
        PRFPAI(0,ISTRA)=PRFPAI(0,ISTRA)+PRFPAI(IATM,ISTRA)
        EOTATI(0,ISTRA)=EOTATI(0,ISTRA)+EOTATI(IATM,ISTRA)
        ERFAAI(0,ISTRA)=ERFAAI(0,ISTRA)+ERFAAI(IATM,ISTRA)
        ERFMAI(0,ISTRA)=ERFMAI(0,ISTRA)+ERFMAI(IATM,ISTRA)
        ERFIAI(0,ISTRA)=ERFIAI(0,ISTRA)+ERFIAI(IATM,ISTRA)
        ERFPAI(0,ISTRA)=ERFPAI(0,ISTRA)+ERFPAI(IATM,ISTRA)
        SPTATI(0,ISTRA)=SPTATI(0,ISTRA)+SPTATI(IATM,ISTRA)
        PGENAI(0,ISTRA)=PGENAI(0,ISTRA)+PGENAI(IATM,ISTRA)
        EGENAI(0,ISTRA)=EGENAI(0,ISTRA)+EGENAI(IATM,ISTRA)
        VGENAI(0,ISTRA)=VGENAI(0,ISTRA)+VGENAI(IATM,ISTRA)
664   CONTINUE
      DO 665 IADV=1,NADVI
        ADDVI(0,ISTRA)=ADDVI(0,ISTRA)+ADDVI(IADV,ISTRA)
665   CONTINUE
      DO 666 ICLV=1,NCLVI
        COLVI(0,ISTRA)=COLVI(0,ISTRA)+COLVI(ICLV,ISTRA)
666   CONTINUE
      DO 667 ISNV=1,NSNVI
        SNAPVI(0,ISTRA)=SNAPVI(0,ISTRA)+SNAPVI(ISNV,ISTRA)
667   CONTINUE
      DO 668 ICPV=1,NCPVI
        COPVI(0,ISTRA)=COPVI(0,ISTRA)+COPVI(ICPV,ISTRA)
668   CONTINUE
      DO 669 IBGV=1,NBGVI
        BGKVI(0,ISTRA)=BGKVI(0,ISTRA)+BGKVI(IBGV,ISTRA)
669   CONTINUE
C
      CALL GETSCL (ISTRA,FATM,FMOL,FION)
C
      IF (.NOT.NLSCL) THEN
        CALL LEER(1)
        WRITE (6,*) 'NO RESCALING DONE (NLSCL=FALSE)'
        CALL LEER(2)
C
C  RESCALE TRACKLENGTH ESTIMATED VOLUME AVERAGED TALLIES TO ENSURE
C  PERFECT PARTICLE BALANCE
C
      ELSEIF (NLSCL) THEN
        FASCL(ISTRA)=FATM
        FMSCL(ISTRA)=FMOL
        FISCL(ISTRA)=FION
C
C  CARRY OUT SCALING OF VOLUME AND SURFACE TALLIES, RESP.
C
C  ATOM TALLIES
C
        DO 2101 IATM=1,NATMI
          DO 111 J=1,NSBOX
            PDENA(IATM,J)=PDENA(IATM,J)*FATM
            EDENA(IATM,J)=EDENA(IATM,J)*FATM
            PAAT(IATM,J)=PAAT(IATM,J)*FATM
            PMAT(IATM,J)=PMAT(IATM,J)*FMOL
            PIAT(IATM,J)=PIAT(IATM,J)*FION
            PGENA(IATM,J)=PGENA(IATM,J)*FATM
            EGENA(IATM,J)=EGENA(IATM,J)*FATM
            VGENA(IATM,J)=VGENA(IATM,J)*FATM
111       CONTINUE
          DO 310 J=1,NLIMPS
            POTAT(IATM,J)=POTAT(IATM,J)*FATM
            PRFAAT(IATM,J)=PRFAAT(IATM,J)*FATM
            PRFMAT(IATM,J)=PRFMAT(IATM,J)*FMOL
            PRFIAT(IATM,J)=PRFIAT(IATM,J)*FION
            EOTAT(IATM,J)=EOTAT(IATM,J)*FATM
            ERFAAT(IATM,J)=ERFAAT(IATM,J)*FATM
            ERFMAT(IATM,J)=ERFMAT(IATM,J)*FMOL
            ERFIAT(IATM,J)=ERFIAT(IATM,J)*FION
            SPTAT(IATM,J)=SPTAT(IATM,J)*FATM
310       CONTINUE
2101    CONTINUE
        DO 2111 IATM=0,NATMI
          PDENAI(IATM,ISTRA)=PDENAI(IATM,ISTRA)*FATM
          EDENAI(IATM,ISTRA)=EDENAI(IATM,ISTRA)*FATM
          PAATI(IATM,ISTRA)=PAATI(IATM,ISTRA)*FATM
          PMATI(IATM,ISTRA)=PMATI(IATM,ISTRA)*FMOL
          PIATI(IATM,ISTRA)=PIATI(IATM,ISTRA)*FION
          POTATI(IATM,ISTRA)=POTATI(IATM,ISTRA)*FATM
          PRFAAI(IATM,ISTRA)=PRFAAI(IATM,ISTRA)*FATM
          PRFMAI(IATM,ISTRA)=PRFMAI(IATM,ISTRA)*FMOL
          PRFIAI(IATM,ISTRA)=PRFIAI(IATM,ISTRA)*FION
          EOTATI(IATM,ISTRA)=EOTATI(IATM,ISTRA)*FATM
          ERFAAI(IATM,ISTRA)=ERFAAI(IATM,ISTRA)*FATM
          ERFMAI(IATM,ISTRA)=ERFMAI(IATM,ISTRA)*FMOL
          ERFIAI(IATM,ISTRA)=ERFIAI(IATM,ISTRA)*FION
          SPTATI(IATM,ISTRA)=SPTATI(IATM,ISTRA)*FATM
          PGENAI(IATM,ISTRA)=PGENAI(IATM,ISTRA)*FATM
          EGENAI(IATM,ISTRA)=EGENAI(IATM,ISTRA)*FATM
          VGENAI(IATM,ISTRA)=VGENAI(IATM,ISTRA)*FATM
2111    CONTINUE
        DO 2112 J=1,NSBOX
          EAAT(J)=EAAT(J)*FATM
          EMAT(J)=EMAT(J)*FMOL
          EIAT(J)=EIAT(J)*FION
2112    CONTINUE
        EAATI(ISTRA)=EAATI(ISTRA)*FATM
        EMATI(ISTRA)=EMATI(ISTRA)*FMOL
        EIATI(ISTRA)=EIATI(ISTRA)*FION
C
C  MOLECULE TALLIES
C
        DO 2115 IMOL=1,NMOLI
          DO 115 J=1,NSBOX
            PDENM(IMOL,J)=PDENM(IMOL,J)*FMOL
            EDENM(IMOL,J)=EDENM(IMOL,J)*FMOL
            PAML(IMOL,J)=PAML(IMOL,J)*FATM
            PMML(IMOL,J)=PMML(IMOL,J)*FMOL
            PIML(IMOL,J)=PIML(IMOL,J)*FION
            PGENM(IMOL,J)=PGENM(IMOL,J)*FMOL
            EGENM(IMOL,J)=EGENM(IMOL,J)*FMOL
            VGENM(IMOL,J)=VGENM(IMOL,J)*FMOL
115       CONTINUE
          DO 315 J=1,NLIMPS
            POTML(IMOL,J)=POTML(IMOL,J)*FMOL
            PRFAML(IMOL,J)=PRFAML(IMOL,J)*FATM
            PRFMML(IMOL,J)=PRFMML(IMOL,J)*FMOL
            PRFIML(IMOL,J)=PRFIML(IMOL,J)*FION
            EOTML(IMOL,J)=EOTML(IMOL,J)*FMOL
            ERFAML(IMOL,J)=ERFAML(IMOL,J)*FATM
            ERFMML(IMOL,J)=ERFMML(IMOL,J)*FMOL
            ERFIML(IMOL,J)=ERFIML(IMOL,J)*FION
            SPTML(IMOL,J)=SPTML(IMOL,J)*FMOL
315       CONTINUE
2115     CONTINUE
        DO 2116 IMOL=0,NMOLI
          PDENMI(IMOL,ISTRA)=PDENMI(IMOL,ISTRA)*FMOL
          EDENMI(IMOL,ISTRA)=EDENMI(IMOL,ISTRA)*FMOL
          PAMLI(IMOL,ISTRA)=PAMLI(IMOL,ISTRA)*FATM
          PMMLI(IMOL,ISTRA)=PMMLI(IMOL,ISTRA)*FMOL
          PIMLI(IMOL,ISTRA)=PIMLI(IMOL,ISTRA)*FION
          POTMLI(IMOL,ISTRA)=POTMLI(IMOL,ISTRA)*FMOL
          PRFAMI(IMOL,ISTRA)=PRFAMI(IMOL,ISTRA)*FATM
          PRFMMI(IMOL,ISTRA)=PRFMMI(IMOL,ISTRA)*FMOL
          PRFIMI(IMOL,ISTRA)=PRFIMI(IMOL,ISTRA)*FION
          EOTMLI(IMOL,ISTRA)=EOTMLI(IMOL,ISTRA)*FMOL
          ERFAMI(IMOL,ISTRA)=ERFAMI(IMOL,ISTRA)*FATM
          ERFMMI(IMOL,ISTRA)=ERFMMI(IMOL,ISTRA)*FMOL
          ERFIMI(IMOL,ISTRA)=ERFIMI(IMOL,ISTRA)*FION
          SPTMLI(IMOL,ISTRA)=SPTMLI(IMOL,ISTRA)*FMOL
          PGENMI(IMOL,ISTRA)=PGENMI(IMOL,ISTRA)*FMOL
          EGENMI(IMOL,ISTRA)=EGENMI(IMOL,ISTRA)*FMOL
          VGENMI(IMOL,ISTRA)=VGENMI(IMOL,ISTRA)*FMOL
2116    CONTINUE
        DO 2117 J=1,NSBOX
          EAML(J)=EAML(J)*FATM
          EMML(J)=EMML(J)*FMOL
          EIML(J)=EIML(J)*FION
2117    CONTINUE
        EAMLI(ISTRA)=EAMLI(ISTRA)*FATM
        EMMLI(ISTRA)=EMMLI(ISTRA)*FMOL
        EIMLI(ISTRA)=EIMLI(ISTRA)*FION
C
C  TEST ION TALLIES
C
        DO 420 IION=1,NIONI
          DO 421 J=1,NSBOX
            PDENI(IION,J)=PDENI(IION,J)*FION
            EDENI(IION,J)=EDENI(IION,J)*FION
            PAIO(IION,J)=PAIO(IION,J)*FATM
            PMIO(IION,J)=PMIO(IION,J)*FMOL
            PIIO(IION,J)=PIIO(IION,J)*FION
            PGENI(IION,J)=PGENI(IION,J)*FION
            EGENI(IION,J)=EGENI(IION,J)*FION
            VGENI(IION,J)=VGENI(IION,J)*FION
421       CONTINUE
          DO 422 J=1,NLIMPS
            POTIO(IION,J)=POTIO(IION,J)*FION
            PRFAIO(IION,J)=PRFAIO(IION,J)*FATM
            PRFMIO(IION,J)=PRFMIO(IION,J)*FMOL
            PRFIIO(IION,J)=PRFIIO(IION,J)*FION
            EOTIO(IION,J)=EOTIO(IION,J)*FION
            ERFAIO(IION,J)=ERFAIO(IION,J)*FATM
            ERFMIO(IION,J)=ERFMIO(IION,J)*FMOL
            ERFIIO(IION,J)=ERFIIO(IION,J)*FION
            SPTIO(IION,J)=SPTIO(IION,J)*FION
422       CONTINUE
420     CONTINUE
        DO 431 IION=0,NIONI
          PDENII(IION,ISTRA)=PDENII(IION,ISTRA)*FION
          EDENII(IION,ISTRA)=EDENII(IION,ISTRA)*FION
          PAIOI(IION,ISTRA)=PAIOI(IION,ISTRA)*FATM
          PMIOI(IION,ISTRA)=PMIOI(IION,ISTRA)*FMOL
          PIIOI(IION,ISTRA)=PIIOI(IION,ISTRA)*FION
          POTIOI(IION,ISTRA)=POTIOI(IION,ISTRA)*FION
          PRFAII(IION,ISTRA)=PRFAII(IION,ISTRA)*FATM
          PRFMII(IION,ISTRA)=PRFMII(IION,ISTRA)*FMOL
          PRFIII(IION,ISTRA)=PRFIII(IION,ISTRA)*FION
          EOTIOI(IION,ISTRA)=EOTIOI(IION,ISTRA)*FION
          ERFAII(IION,ISTRA)=ERFAII(IION,ISTRA)*FATM
          ERFMII(IION,ISTRA)=ERFMII(IION,ISTRA)*FMOL
          ERFIII(IION,ISTRA)=ERFIII(IION,ISTRA)*FION
          SPTIOI(IION,ISTRA)=SPTIOI(IION,ISTRA)*FION
          PGENII(IION,ISTRA)=PGENII(IION,ISTRA)*FION
          EGENII(IION,ISTRA)=EGENII(IION,ISTRA)*FION
          VGENII(IION,ISTRA)=VGENII(IION,ISTRA)*FION
431     CONTINUE
        DO 432 J=1,NSBOX
          EAIO(J)=EAIO(J)*FATM
          EMIO(J)=EMIO(J)*FMOL
          EIIO(J)=EIIO(J)*FION
432     CONTINUE
        EAIOI(ISTRA)=EAIOI(ISTRA)*FATM
        EMIOI(ISTRA)=EMIOI(ISTRA)*FMOL
        EIIOI(ISTRA)=EIIOI(ISTRA)*FION
C
C  ADDITIONAL TRACKLENGTH ESTIMATED TALLIES
C
        DO 423 IADV=1,NADVI
          IF (IADRC(IADV).EQ.1) THEN
            FADD=FATM
          ELSEIF (IADRC(IADV).EQ.2) THEN
            FADD=FMOL
          ELSEIF (IADRC(IADV).EQ.3) THEN
            FADD=FION
          ELSE
            GOTO 423
          ENDIF
          DO 424 J=1,NSBOX
            ADDV(IADV,J)=ADDV(IADV,J)*FADD
424       CONTINUE
          ADDVI(IADV,ISTRA)=ADDVI(IADV,ISTRA)*FADD
423     CONTINUE
C
C  ADDITIONAL COLLISION ESTIMATED TALLIES
C
        DO 426 ICLV=1,NCLVI
          IF (ICLRC(ICLV).EQ.1) THEN
            FADD=FATM
          ELSEIF (ICLRC(ICLV).EQ.2) THEN
            FADD=FMOL
          ELSEIF (ICLRC(ICLV).EQ.3) THEN
            FADD=FION
          ELSE
            GOTO 426
          ENDIF
          DO 427 J=1,NSBOX
            COLV(ICLV,J)=COLV(ICLV,J)*FADD
427       CONTINUE
          COLVI(ICLV,ISTRA)=COLVI(ICLV,ISTRA)*FADD
426     CONTINUE
C
C  SNAPSHOT TALLIES
C
C  TO BE WRITTEN
C
C  TALLIES FOR COUPLING TO FLUID PLASMA CODE
C
        DO 435 ICPV=1,NCPVI
          IF (ICPRC(ICPV).EQ.1) THEN
            FADD=FATM
          ELSEIF (ICPRC(ICPV).EQ.2) THEN
            FADD=FMOL
          ELSEIF (ICPRC(ICPV).EQ.3) THEN
            FADD=FION
          ELSE
            GOTO 435
          ENDIF
          DO 436 J=1,NSBOX
c...sltmp
            STOP 'STOP: COPV BUSINESS Z'
            COPV(ICPV,J)=COPV(ICPV,J)*FADD
436       CONTINUE
          COPVI(ICPV,ISTRA)=COPVI(ICPV,ISTRA)*FADD
435     CONTINUE
C
C  TALLIES FOR BGK SELF COLLISION ITERATIONS
C
        DO 437 IBGV=1,NBGVI
          IF (IBGRC(IBGV).EQ.1) THEN
            FADD=FATM
          ELSEIF (IBGRC(IBGV).EQ.2) THEN
            FADD=FMOL
          ELSEIF (IBGRC(IBGV).EQ.3) THEN
            FADD=FION
          ELSE
            GOTO 437
          ENDIF
          DO 438 J=1,NSBOX
            BGKV(IBGV,J)=BGKV(IBGV,J)*FADD
438       CONTINUE
          BGKVI(IBGV,ISTRA)=BGKVI(IBGV,ISTRA)*FADD
437     CONTINUE
C
C  BULK ION TALLIES
C
        DO 447 IPLS=1,NPLSI
          DO 447 J=1,NSBOX
            PAPL(IPLS,J)=PAPL(IPLS,J)*FATM
            PMPL(IPLS,J)=PMPL(IPLS,J)*FMOL
            PIPL(IPLS,J)=PIPL(IPLS,J)*FION
447     CONTINUE
        DO 448 IPLS=0,NPLSI
          PAPLI(IPLS,ISTRA)=PAPLI(IPLS,ISTRA)*FATM
          PMPLI(IPLS,ISTRA)=PMPLI(IPLS,ISTRA)*FMOL
          PIPLI(IPLS,ISTRA)=PIPLI(IPLS,ISTRA)*FION
448     CONTINUE
        DO 449 J=1,NSBOX
          EAPL(J)=EAPL(J)*FATM
          EMPL(J)=EMPL(J)*FMOL
          EIPL(J)=EIPL(J)*FION
449     CONTINUE
        EAPLI(ISTRA)=EAPLI(ISTRA)*FATM
        EMPLI(ISTRA)=EMPLI(ISTRA)*FMOL
        EIPLI(ISTRA)=EIPLI(ISTRA)*FION
C
C  ELECTRON TALLIES
C
        DO 551 J=1,NSBOX
          PAEL(J)=PAEL(J)*FATM
          PMEL(J)=PMEL(J)*FMOL
          PIEL(J)=PIEL(J)*FION
551     CONTINUE
        PAELI(ISTRA)=PAELI(ISTRA)*FATM
        PMELI(ISTRA)=PMELI(ISTRA)*FMOL
        PIELI(ISTRA)=PIELI(ISTRA)*FION
        DO 552 J=1,NSBOX
          EAEL(J)=EAEL(J)*FATM
          EMEL(J)=EMEL(J)*FMOL
          EIEL(J)=EIEL(J)*FION
552     CONTINUE
        EAELI(ISTRA)=EAELI(ISTRA)*FATM
        EMELI(ISTRA)=EMELI(ISTRA)*FMOL
        EIELI(ISTRA)=EIELI(ISTRA)*FION
C
        CALL LEER(1)
        WRITE (6,*) ('RESCALING OF TRACKLENGTH TALLIES COMPLETED')
        WRITE (6,*) ('RESCALING FACTORS:                        ')
        CALL MASR3 ('FATM,FMOL,FION          ',FATM, FMOL, FION)
        CALL LEER(2)
C
      ENDIF
C
C   ALGEBRAIC EXPRESSION IN TALLIES 801--900
C
      IF (NALVI.GT.0.OR.NALSI.GT.0) THEN
C
        CALL ALGTAL
C
        DO 830 IALV=1,NALVI
          CALL INTTAL (ALGV,VOL,IALV,NALV,NSBOX,ALGVI(IALV,ISTRA),
     .                 NR1ST,NP2ND,NT3RD,NBMLT)
830     CONTINUE
C
        DO 832 IALS=1,NALSI
          ALGSI(IALS,ISTRA)=0.
          DO 831 J=1,NLIMPS
            ALGSI(IALS,ISTRA)=ALGSI(IALS,ISTRA)+ALGS(IALS,J)
831       CONTINUE
832     CONTINUE
C
      ENDIF
C
C  SCALE STANDARD DEVIATIONS, WHICH ARE NOT GIVEN IN % REL.ERROR
C  1/XMCP IS INCLUDED IN ZVOLIN,ZW,ZWW,... FOR TALLY AVERAGING
C  THEREFORE IT MUST BE MULTIPLIED HERE BECAUSE ONLY FLUX SCALING
C
      IF (XMCP(ISTRA).LE.1.D0) GOTO 950
C
      DO 900 ISDV=1,NSIGCI
        DO 901 ITAL=1,NTALR
          IF (IIHC(1,ISDV).NE.ITAL) GOTO 901
          ISPZ=MAX(1,IGHC(1,ISDV))
          IF (SCLTAL(ISPZ,ITAL).EQ.1) THEN
            DO 902 ICELL=1,NSBOX
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTRA)
902         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZVOLNT*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLNT*XMCP(ISTRA)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.2) THEN
            DO 903 ICELL=1,NSBOX
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZW*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZW*
     .                             XMCP(ISTRA)
903         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZW*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZW*XMCP(ISTRA)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.3) THEN
            DO 904 ICELL=1,NSBOX
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTRA)
904         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZVOLWT*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLWT*XMCP(ISTRA)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.4) THEN
            DO 905 ICELL=1,NSBOX
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTRA)
905         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZWW*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZWW*XMCP(ISTRA)
          ENDIF
901     CONTINUE
C
        DO 911 ITAL=1,NTALR
          IF (IIHC(2,ISDV).NE.ITAL) GOTO 911
          ISPZ=MAX(1,IGHC(2,ISDV))
          IF (SCLTAL(ISPZ,ITAL).EQ.1) THEN
            DO 912 ICELL=1,NSBOX
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTRA)
912         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZVOLNT*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLNT*XMCP(ISTRA)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.2) THEN
            DO 913 ICELL=1,NSBOX
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZW*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZW*
     .                             XMCP(ISTRA)
913         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZW*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZW*XMCP(ISTRA)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.3) THEN
            DO 914 ICELL=1,NSBOX
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTRA)
914         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZVOLWT*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLWT*XMCP(ISTRA)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.4) THEN
            DO 915 ICELL=1,NSBOX
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTRA)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTRA)
915         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZWW*XMCP(ISTRA)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZWW*XMCP(ISTRA)
          ENDIF
911     CONTINUE
900   CONTINUE
C
950   CONTINUE
C
C  CALL INTERFACE TO OTHER CODES TO RETURN DATA. STRATUM ISTRA
C
      IF (NMODE.GT.0) THEN
        IESTR=ISTRA
        ISTRAA=ISTRA
        ISTRAE=ISTRA
        CALL IF3COP(ISTRAA,ISTRAE)
      ENDIF
C
C  WRITE RESULTS FOR THIS STRATUM ON TEMP. FILE
C
      IESTR=ISTRA
      IF (NFILEN.EQ.1) THEN
        CALL WRSTRT(ISTRA,NSTRAI,NESTIM,NSDVI,ESTIM,SDVI,
     .              NSBGK,SDVI_BGK,NSCOP,SDVI_COP,TRCFLE)
      ENDIF
C
C  UPDATE TALLIES FOR  "SUM OVER STRATA"
C
      IF (NSTRAI.EQ.1) GOTO 1111
C
      DO 1101 IMOL=0,NMOLI
        PPMLI (IMOL,0)=PPMLI (IMOL,0)+PPMLI (IMOL,ISTRA)
        WTOTM (IMOL,0)=WTOTM (IMOL,0)+WTOTM (IMOL,ISTRA)
        PDENMI(IMOL,0)=PDENMI(IMOL,0)+PDENMI(IMOL,ISTRA)
        EDENMI(IMOL,0)=EDENMI(IMOL,0)+EDENMI(IMOL,ISTRA)
        PAMLI (IMOL,0)=PAMLI (IMOL,0)+PAMLI (IMOL,ISTRA)
        PMMLI (IMOL,0)=PMMLI (IMOL,0)+PMMLI (IMOL,ISTRA)
        PIMLI (IMOL,0)=PIMLI (IMOL,0)+PIMLI (IMOL,ISTRA)
        POTMLI(IMOL,0)=POTMLI(IMOL,0)+POTMLI(IMOL,ISTRA)
        PRFAMI(IMOL,0)=PRFAMI(IMOL,0)+PRFAMI(IMOL,ISTRA)
        PRFMMI(IMOL,0)=PRFMMI(IMOL,0)+PRFMMI(IMOL,ISTRA)
        PRFIMI(IMOL,0)=PRFIMI(IMOL,0)+PRFIMI(IMOL,ISTRA)
        PRFPMI(IMOL,0)=PRFPMI(IMOL,0)+PRFPMI(IMOL,ISTRA)
        EOTMLI(IMOL,0)=EOTMLI(IMOL,0)+EOTMLI(IMOL,ISTRA)
        ERFAMI(IMOL,0)=ERFAMI(IMOL,0)+ERFAMI(IMOL,ISTRA)
        ERFMMI(IMOL,0)=ERFMMI(IMOL,0)+ERFMMI(IMOL,ISTRA)
        ERFIMI(IMOL,0)=ERFIMI(IMOL,0)+ERFIMI(IMOL,ISTRA)
        ERFPMI(IMOL,0)=ERFPMI(IMOL,0)+ERFPMI(IMOL,ISTRA)
        SPTMLI(IMOL,0)=SPTMLI(IMOL,0)+SPTMLI(IMOL,ISTRA)
        PGENMI(IMOL,0)=PGENMI(IMOL,0)+PGENMI(IMOL,ISTRA)
        EGENMI(IMOL,0)=EGENMI(IMOL,0)+EGENMI(IMOL,ISTRA)
        VGENMI(IMOL,0)=VGENMI(IMOL,0)+VGENMI(IMOL,ISTRA)
1101  CONTINUE
      DO 1102 IATM=0,NATMI
        PPATI (IATM,0)=PPATI (IATM,0)+PPATI (IATM,ISTRA)
        WTOTA (IATM,0)=WTOTA (IATM,0)+WTOTA (IATM,ISTRA)
        PDENAI(IATM,0)=PDENAI(IATM,0)+PDENAI(IATM,ISTRA)
        EDENAI(IATM,0)=EDENAI(IATM,0)+EDENAI(IATM,ISTRA)
        PAATI (IATM,0)=PAATI (IATM,0)+PAATI (IATM,ISTRA)
        PMATI (IATM,0)=PMATI (IATM,0)+PMATI (IATM,ISTRA)
        PIATI (IATM,0)=PIATI (IATM,0)+PIATI (IATM,ISTRA)
        POTATI(IATM,0)=POTATI(IATM,0)+POTATI(IATM,ISTRA)
        PRFAAI(IATM,0)=PRFAAI(IATM,0)+PRFAAI(IATM,ISTRA)
        PRFMAI(IATM,0)=PRFMAI(IATM,0)+PRFMAI(IATM,ISTRA)
        PRFIAI(IATM,0)=PRFIAI(IATM,0)+PRFIAI(IATM,ISTRA)
        PRFPAI(IATM,0)=PRFPAI(IATM,0)+PRFPAI(IATM,ISTRA)
        EOTATI(IATM,0)=EOTATI(IATM,0)+EOTATI(IATM,ISTRA)
        ERFAAI(IATM,0)=ERFAAI(IATM,0)+ERFAAI(IATM,ISTRA)
        ERFMAI(IATM,0)=ERFMAI(IATM,0)+ERFMAI(IATM,ISTRA)
        ERFIAI(IATM,0)=ERFIAI(IATM,0)+ERFIAI(IATM,ISTRA)
        ERFPAI(IATM,0)=ERFPAI(IATM,0)+ERFPAI(IATM,ISTRA)
        SPTATI(IATM,0)=SPTATI(IATM,0)+SPTATI(IATM,ISTRA)
        PGENAI(IATM,0)=PGENAI(IATM,0)+PGENAI(IATM,ISTRA)
        EGENAI(IATM,0)=EGENAI(IATM,0)+EGENAI(IATM,ISTRA)
        VGENAI(IATM,0)=VGENAI(IATM,0)+VGENAI(IATM,ISTRA)
1102  CONTINUE
      DO 1103 IION=0,NIONI
        PPIOI (IION,0)=PPIOI (IION,0)+PPIOI (IION,ISTRA)
        EELFI (IION,0)=EELFI (IION,0)+EELFI (IION,ISTRA)
        WTOTI (IION,0)=WTOTI (IION,0)+WTOTI (IION,ISTRA)
        PDENII(IION,0)=PDENII(IION,0)+PDENII(IION,ISTRA)
        EDENII(IION,0)=EDENII(IION,0)+EDENII(IION,ISTRA)
        PAIOI (IION,0)=PAIOI (IION,0)+PAIOI (IION,ISTRA)
        PMIOI (IION,0)=PMIOI (IION,0)+PMIOI (IION,ISTRA)
        PIIOI (IION,0)=PIIOI (IION,0)+PIIOI (IION,ISTRA)
        POTIOI(IION,0)=POTIOI(IION,0)+POTIOI(IION,ISTRA)
        PRFAII(IION,0)=PRFAII(IION,0)+PRFAII(IION,ISTRA)
        PRFMII(IION,0)=PRFMII(IION,0)+PRFMII(IION,ISTRA)
        PRFIII(IION,0)=PRFIII(IION,0)+PRFIII(IION,ISTRA)
        PRFPII(IION,0)=PRFPII(IION,0)+PRFPII(IION,ISTRA)
        EOTIOI(IION,0)=EOTIOI(IION,0)+EOTIOI(IION,ISTRA)
        ERFAII(IION,0)=ERFAII(IION,0)+ERFAII(IION,ISTRA)
        ERFMII(IION,0)=ERFMII(IION,0)+ERFMII(IION,ISTRA)
        ERFIII(IION,0)=ERFIII(IION,0)+ERFIII(IION,ISTRA)
        ERFPII(IION,0)=ERFPII(IION,0)+ERFPII(IION,ISTRA)
        SPTIOI(IION,0)=SPTIOI(IION,0)+SPTIOI(IION,ISTRA)
        PGENII(IION,0)=PGENII(IION,0)+PGENII(IION,ISTRA)
        EGENII(IION,0)=EGENII(IION,0)+EGENII(IION,ISTRA)
        VGENII(IION,0)=VGENII(IION,0)+VGENII(IION,ISTRA)
1103  CONTINUE
      DO 1104 IPLS=0,NPLSI
        WTOTP (IPLS,0)=WTOTP (IPLS,0)+WTOTP (IPLS,ISTRA)
        PAPLI (IPLS,0)=PAPLI (IPLS,0)+PAPLI (IPLS,ISTRA)
        PMPLI (IPLS,0)=PMPLI (IPLS,0)+PMPLI (IPLS,ISTRA)
        PIPLI (IPLS,0)=PIPLI (IPLS,0)+PIPLI (IPLS,ISTRA)
        POTPLI(IPLS,0)=POTPLI(IPLS,0)+POTPLI(IPLS,ISTRA)
        EOTPLI(IPLS,0)=EOTPLI(IPLS,0)+EOTPLI(IPLS,ISTRA)
        SPTPLI(IPLS,0)=SPTPLI(IPLS,0)+SPTPLI(IPLS,ISTRA)
1104  CONTINUE
      DO 1105 IADV=0,NADVI
        ADDVI(IADV,0)=ADDVI(IADV,0)+ADDVI(IADV,ISTRA)
1105  CONTINUE
      DO 1106 ICLV=0,NCLVI
        COLVI(ICLV,0)=COLVI(ICLV,0)+COLVI(ICLV,ISTRA)
1106  CONTINUE
      DO 1107 ISNV=0,NSNVI
        SNAPVI(ISNV,0)=SNAPVI(ISNV,0)+SNAPVI(ISNV,ISTRA)
1107  CONTINUE
      DO 1108 ICPV=0,NCPVI
        COPVI(ICPV,0)=COPVI(ICPV,0)+COPVI(ICPV,ISTRA)
1108  CONTINUE
      DO 1109 IBGV=0,NBGVI
        BGKVI(IBGV,0)=BGKVI(IBGV,0)+BGKVI(IBGV,ISTRA)
1109  CONTINUE
      PAELI(0)=PAELI(0)+PAELI(ISTRA)
      PMELI(0)=PMELI(0)+PMELI(ISTRA)
      PIELI(0)=PIELI(0)+PIELI(ISTRA)
C
      EAELI(0)=EAELI(0)+EAELI(ISTRA)
      EAATI(0)=EAATI(0)+EAATI(ISTRA)
      EAMLI(0)=EAMLI(0)+EAMLI(ISTRA)
      EAIOI(0)=EAIOI(0)+EAIOI(ISTRA)
      EAPLI(0)=EAPLI(0)+EAPLI(ISTRA)
C
      EMELI(0)=EMELI(0)+EMELI(ISTRA)
      EMATI(0)=EMATI(0)+EMATI(ISTRA)
      EMMLI(0)=EMMLI(0)+EMMLI(ISTRA)
      EMIOI(0)=EMIOI(0)+EMIOI(ISTRA)
      EMPLI(0)=EMPLI(0)+EMPLI(ISTRA)
C
      EIELI(0)=EIELI(0)+EIELI(ISTRA)
      EIATI(0)=EIATI(0)+EIATI(ISTRA)
      EIMLI(0)=EIMLI(0)+EIMLI(ISTRA)
      EIIOI(0)=EIIOI(0)+EIIOI(ISTRA)
      EIPLI(0)=EIPLI(0)+EIPLI(ISTRA)
C
      EPATI(0)=EPATI(0)+EPATI(ISTRA)
      EPMLI(0)=EPMLI(0)+EPMLI(ISTRA)
      EPIOI(0)=EPIOI(0)+EPIOI(ISTRA)
C
C
      FLUXT(0)=FLUXT(0)+FLUXT(ISTRA)
      XMCP(0)=XMCP(0)+XMCP(ISTRA)
      PTRASH(0)=PTRASH(0)+PTRASH(ISTRA)
      ETRASH(0)=ETRASH(0)+ETRASH(ISTRA)
      ETOTA(0)=ETOTA(0)+ETOTA(ISTRA)
      ETOTM(0)=ETOTM(0)+ETOTM(ISTRA)
      ETOTI(0)=ETOTI(0)+ETOTI(ISTRA)
      ETOTP(0)=ETOTP(0)+ETOTP(ISTRA)
C
C
      DO 1160 J=1,(NESTM1+NESTM2)*NSMSTRA
        RWK(J)=RWK(J)+ESTIM(J)
1160  CONTINUE
C
      DO 1170 ISDV=1,NSIGCI
        DO 1172 ICELL=1,NSBOX
          STVC(0,ISDV,ICELL)=STVC(0,ISDV,ICELL)+SIGMAC(0,ISDV,ICELL)
          STVC(1,ISDV,ICELL)=STVC(1,ISDV,ICELL)+SIGMAC(1,ISDV,ICELL)**2
          STVC(2,ISDV,ICELL)=STVC(2,ISDV,ICELL)+SIGMAC(2,ISDV,ICELL)**2
1172    CONTINUE
        STVCS(0,ISDV)=STVCS(0,ISDV)+SGMCS(0,ISDV)
        STVCS(1,ISDV)=STVCS(1,ISDV)+SGMCS(1,ISDV)**2
        STVCS(2,ISDV)=STVCS(2,ISDV)+SGMCS(2,ISDV)**2
1170  CONTINUE
C
1111  CONTINUE
        WRITE(6,*) 'CPU TIME USED UNTIL END OF STRATUM ISTRA '
        WRITE(6,*) 'ISTRA, CPU(S) ',ISTRA,SECOND_OWN()
        CALL LEER(2)
c slmod begin - tr
c...    Sum momentum loss components for this stratum:
        DO i1 = 1, NRAD
          DO i2 = 1, NMOMCHA
            IF (i2.EQ.14) CYCLE
            copv2(3,i1,0) = copv2(3,i1,0) + copv2(3,i1,i2)
          ENDDO
        ENDDO
        DO I1 = 1, NCPV
          DO I2 = 1, NRAD
            DO I3 = 0, NMOMCHA
              COPV3(I1,I2,I3) = COPV3(I1,I2,I3) + COPV2(I1,I2,I3) 
            ENDDO
          ENDDO
        ENDDO

c...    Store stratum tallies for specified additional cells:
        DO ISTR=1,NSTRDAT
          DO IATM=1,NATMI
            PSTRDATA(IATM,ISTR,ISTRA)=PDENA(IATM,NSURF+STRDAT(ISTR))
            ESTRDATA(IATM,ISTR,ISTRA)=EDENA(IATM,NSURF+STRDAT(ISTR))
          ENDDO                                                    
          DO IMOL=1,NMOLI                                          
            PSTRDATM(IMOL,ISTR,ISTRA)=PDENM(IMOL,NSURF+STRDAT(ISTR))
            ESTRDATM(IMOL,ISTR,ISTRA)=EDENM(IMOL,NSURF+STRDAT(ISTR))
          ENDDO
        ENDDO
c slmod end
1000  CONTINUE
C
C*** STRATA LOOP FINISHED *******************************************
C
C    STATISTICS, SUM OVER STRATA
C
      IF (NSTRAI.EQ.1.OR.XMCP(0).LE.1) GOTO 2000
C
      DO 1207 K=1,NSIGVI
        DO 1208 I=1,NSBOX
          ST=MAX(0.D0,STV(K,I))
          STV(K,I)=SQRT(ST)/(ABS(EE(K,I))+EPS60)
1208    CONTINUE
        ST=MAX(0.D0,STVS(K))
        STVS(K)=SQRT(ST)/(ABS(EES(K))+EPS60)
1207  CONTINUE
C
      IF (NSIGI_BGK.GT.0) THEN
        DO 1217 K=1,NBGVI_STAT
          DO 1218 I=1,NSBOX
            ST=MAX(0.D0,STV_BGK(K,I))
            STV_BGK(K,I)=SQRT(ST)/(ABS(EE_BGK(K,I))+EPS60)
1218      CONTINUE
          ST=MAX(0.D0,STVS_BGK(K))
          STVS_BGK(K)=SQRT(ST)/(ABS(EES_BGK(K))+EPS60)
1217    CONTINUE
      ENDIF
C
      IF (NSIGI_COP.GT.0) THEN
        DO K=1,NCPVI_STAT
          DO I=1,NSBOX
            ST=MAX(0.D0,STV_COP(K,I))
            STV_COP(K,I)=SQRT(ST)/(ABS(EE_COP(K,I))+EPS60)
          END DO
          ST=MAX(0.D0,STVS_COP(K))
          STVS_COP(K)=SQRT(ST)/(ABS(EES_COP(K))+EPS60)
        END DO
      ENDIF
C
      DO 1221 K=1,NSIGSI
        DO 1222 J=1,NLIMPS
          FFF=FF(K,J)*FF(K,J)
          ST=MAX(0.D0,STVW(K,J))
          STVW(K,J)=SQRT(ST)/(ABS(FF(K,J))+EPS60)
1222    CONTINUE
        ST=MAX(0.D0,STVWS(K))
        STVWS(K)=SQRT(ST)/(ABS(FFS(K))+EPS60)
1221  CONTINUE
C
      DO 1240 ISDV=1,NSIGCI
        DO 1242 ICELL=1,NSBOX
          STVC(1,ISDV,ICELL)=SQRT(MAX(0.D0,STVC(1,ISDV,ICELL)))
          STVC(2,ISDV,ICELL)=SQRT(MAX(0.D0,STVC(2,ISDV,ICELL)))
1242    CONTINUE
        STVCS(1,ISDV)=SQRT(MAX(0.D0,STVCS(1,ISDV)))
        STVCS(2,ISDV)=SQRT(MAX(0.D0,STVCS(2,ISDV)))
1240  CONTINUE
C
C  CONVERT RELATIVE ERRORS TO %-ERRORS
C
      DO 1250 J=NESTIM*NSMSTRA+1,NESTIM*NSMSTRA+NSDVI1
        RWK(J)=RWK(J)*100.D0
1250  CONTINUE
C
      IF (NSIGI_BGK.GT.0) THEN
        DO 1251 IB=1,NBGVI_STAT
          STVS_BGK(IB)=STVS_BGK(IB)*100.D0
          DO 1252 J=1,NSBOX
            STV_BGK(IB,J)=STV_BGK(IB,J)*100.D0
1252      CONTINUE
1251    CONTINUE
      ENDIF
C
      IF (NSIGI_COP.GT.0) THEN
        DO IB=1,NCPVI_STAT
          STVS_COP(IB)=STVS_COP(IB)*100.D0
          DO J=1,NSBOX
            STV_COP(IB,J)=STV_COP(IB,J)*100.D0
          END DO
        END DO
      ENDIF
C
C  PUT SUM OVER STRATA BACK ONTO CESTIM, CSDVI
C
      IF (NSMSTRA == 1) THEN
        DO 1260 J=1,NESTIM
          ESTIM(J)=RWK(J)
1260    CONTINUE
        DO 1270 J=1,NSDVI
          SDVI(J)=RWK(NESTIM+J)
1270    CONTINUE
        IF (NSIGI_BGK.GT.0) THEN
          DO 1271 IB=1,NBGVI_STAT
            SGMS_BGK(IB)=STVS_BGK(IB)
            DO 1272 J=1,NSBOX
              SIGMA_BGK(IB,J)=STV_BGK(IB,J)
1272        CONTINUE
1271      CONTINUE
        ENDIF
        IF (NSIGI_COP.GT.0) THEN
          DO 1273 IC=1,NCPVI_STAT
            SGMS_COP(IC)=STVS_COP(IC)
            DO 1274 J=1,NSBOX
              SIGMA_COP(IC,J)=STV_COP(IC,J)
1274        CONTINUE
1273      CONTINUE
        ENDIF
C
C   ALGEBRAIC EXPRESSION IN TALLIES, SUM OVER STRATA  1271--1279
C
        IF (NALVI.GT.0.OR.NALSI.GT.0) THEN
C
          CALL ALGTAL
C
          DO 1571 IALV=1,NALVI
            CALL INTTAL (ALGV,VOL,IALV,NALV,NSBOX,ALGVI(IALV,0),
     .                   NR1ST,NP2ND,NT3RD,NBMLT)
1571      CONTINUE
C
          DO 1572 IALS=1,NALSI
            ALGSI(IALS,0)=0.
            DO 1573 J=1,NLIMPS
              ALGSI(IALS,0)=ALGSI(IALS,0)+ALGS(IALS,J)
1573        CONTINUE
1572      CONTINUE
C
        ENDIF
C
C  WRITE RESULTS FOR SUM OVER STRATA ON TEMP. FILE
C
        IESTR=0
        IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN

          WRITE(0,*) 'MARK: WRITING SUM TO TEMP FILE' 

          CALL WRSTRT(0,NSTRAI,NESTIM,NSDVI,ESTIM,SDVI,
     .                NSBGK,SDVI_BGK,NSCOP,SDVI_COP,TRCFLE)
        ENDIF
      ENDIF
C
C
2000  CONTINUE
C
C  SAVE OR RESTORE SOME DATA FOR "EIRENE RECALL OPTION NFILE.NE.0"
C  FROM FILE "FT11"
C  NOTE: RECORD IRC=3 MAY BE USED IN INTERFACING ROUTINE INFCOP
C
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
        IF (TRCFLE) WRITE (6,*) 'WRITE DATA FOR RECALL OPTION '
        IRC=1
        WRITE (11,REC=IRC) SPEZ
        IF (TRCFLE)   WRITE (6,*) 'WRITE 11  IRC= ',IRC
        IRC=2
        WRITE (11,REC=IRC) OUTAU
        IF (TRCFLE)   WRITE (6,*) 'WRITE 11  IRC= ',IRC
      ELSEIF (NFILEN.EQ.2.OR.NFILEN.EQ.7) THEN
        IF (TRCFLE) WRITE (6,*) 'READ DATA FOR RECALL OPTION'
        IRC=1
        READ (11,REC=IRC) SPEZ
        IF (TRCFLE)   WRITE (6,*) 'READ 11  IRC= ',IRC
        IRC=2
        READ (11,REC=IRC) OUTAU
        IF (TRCFLE)   WRITE (6,*) 'READ 11  IRC= ',IRC
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE STORE(IFLAG)
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        OPEN (UNIT=16,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      ENDIF
      GOTO (10,20,30,40,50,60,70,80,90),IFLAG
      WRITE (6,*) 'IFLAG OUT OF RANGE IN SUBR. STORE '
      WRITE (6,*) 'EXIT CALLED '
      CALL EXIT
C  LOCATE
10    CONTINUE
20    CONTINUE
C  IONIZATION
30    CONTINUE
40    CONTINUE
50    CONTINUE
C  SURFACE
60    CONTINUE
70    CONTINUE
80    CONTINUE
C  NEW CELL
90    CONTINUE
      RETURN
      END
C
      SUBROUTINE WRGEOM(TRCFLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGEOM'
      INCLUDE 'CGRID'
      INCLUDE 'CPOLYG'
      INCLUDE 'CTRIG'
      INCLUDE 'CLGIN'
      INCLUDE 'CADGEO'
      LOGICAL TRCFLE
      DIMENSION RCGM(NCGM),RCGRD(NCGRD),RCPLYG(NCPLYG),RCTRIG(NCTRIG),
     .          RCLGN(NLGIN),RCADG(NADGEO)
      DIMENSION ICGM(MCGM),ICGRD(MCGRD),ICPLYG(MCPLYG),ICTRIG(MCTRIG),
     .          ICLGN(MLGIN),ICADG(MADGEO)
      DIMENSION LCLGN(LLGIN)
      EQUIVALENCE (RCGM(1),VOLADD(1)),
     .            (ICGM(1),NPOINT(1,1)),
     .            (RCGRD(1),RSURF(1)),
     .            (ICGRD(1),LEVGEO),
     .            (RCPLYG(1),VPLX(1,1)),
     .            (ICPLYG(1),NRPLG),
     .            (RCTRIG(1),XTRIAN(1)),
     .            (ICTRIG(1),NECKE(1,1)),
     .            (RCLGN(1),RLWMN(0)),
     .            (ICLGN(1),ILSWCH(0)),
     .            (LCLGN(1),IGFIL(0)),
     .            (RCADG(1),A0LM(1)),
     .            (ICADG(1),NLIMI)
C
      OPEN (UNIT=12,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 12
      WRITE (12) RCGM,ICGM
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCGM,ICGM '
      WRITE (12) RCGRD,ICGRD
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCGRD,ICGRD'
      WRITE (12) RCPLYG,ICPLYG
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCPLYG,ICPLYG'
      WRITE (12) RCTRIG,ICTRIG
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCTRIG,ICTRIG'
      WRITE (12) RCLGN,ICLGN,LCLGN
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCLGN,ICLGN,LCLGN'
      WRITE (12) RCADG,ICADG
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCADG,ICADG'
      CLOSE (UNIT=12)
      RETURN
C
      ENTRY RGEOM(TRCFLE)
      OPEN (UNIT=12,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 12
      READ (12) RCGM,ICGM
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCGM,ICGM '
      READ (12) RCGRD,ICGRD
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCGRD,ICGRD'
      READ (12) RCPLYG,ICPLYG
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCPLYG,ICPLYG'
      READ (12) RCTRIG,ICTRIG
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCTRIG,ICTRIG'
      READ (12) RCLGN,ICLGN,LCLGN
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCLGN,ICLGN,LCLGN'
      READ (12) RCADG,ICADG
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCADG,ICADG'
      CLOSE (UNIT=12)
      RETURN
      END
C
C
      SUBROUTINE WRPLAM(TRCFLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'COMSOU'
      INCLUDE 'COMXS'
      INCLUDE 'CZT1'
      INCLUDE 'CSTEP'
      LOGICAL TRCFLE
      DIMENSION RCMUSR(NUSR),RCMDTA(NMDTA),RCMAMF(NAMF),RCZT1(NZT1),
     .          RCMSOU(NOMSOU),RCSTEP(NSTPP)
      DIMENSION ICMUSR(MUSR),ICMDTA(MMDTA),ICMAMF(MAMF),
     .          ICMSOU(MOMSOU),ICSTEP(MSTPP)
      LOGICAL   LCMUSR(LUSR),LCMSOU(LOMSOU)
      EQUIVALENCE (RCMUSR(1),TEIN(1)),
     .            (RCMDTA(1),TABDS1(1,1)),
     .            (RCMAMF(1),CREAC(1,0,1)),
     .            (RCMSOU(1),FLUX(1)),
     .            (RCSTEP(1),FLSTEP(0,1,1))
      EQUIVALENCE (ICMUSR(1),NSPA),
     .            (ICMDTA(1),MODCOL(1,1,1,0)),
     .            (ICMAMF(1),NREACI),
     .            (ICMSOU(1),IVLSF(1)),
     .            (ICSTEP(1),IRSTEP(1,1))
      EQUIVALENCE (LCMUSR(1),LGVAC(1,0)),
     .            (LCMSOU(1),NLPNT(1))
      EQUIVALENCE (RCZT1(1),RSQDVI(1))
C
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      WRITE (13) RCMUSR,ICMUSR,LCMUSR
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMUSR,ICMUSR,LCMUSR '
      WRITE (13) RCMDTA,ICMDTA
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMDTA,ICMDTA'
      WRITE (13) RCMAMF,ICMAMF
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMAMF,ICMAMF'
      WRITE (13) RCZT1
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCZT1'
      WRITE (13) RCMSOU,ICMSOU,LCMSOU
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'
      WRITE (13) RCSTEP,ICSTEP
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCSTEP,ICSTEP'
      CLOSE (UNIT=13)
      RETURN
C
      ENTRY RPLAM(TRCFLE)
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      READ (13) RCMUSR,ICMUSR,LCMUSR
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMUSR,ICMUSR,LCMUSR '
      READ (13) RCMDTA,ICMDTA
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMDTA,ICMDTA'
      READ (13) RCMAMF,ICMAMF
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMAMF,ICMAMF'
      READ (13) RCZT1
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCZT1'
      READ (13) RCMSOU,ICMSOU,LCMSOU
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMSOU,ICMSOU,LCMSOU'
      READ (13) RCSTEP,ICSTEP
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCSTEP,ICSTEP'
      CLOSE (UNIT=13)
      RETURN
      END
C
C
      SUBROUTINE UPDATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C ESTIMATORS ARE UPDATED FOR EACH TRACK TAKING T/VEL SEC.
C T (CM) IS STORED ON CLPD ARRAY FOR ONE OR MORE CELLS, THAT HAVE
C BEEN CROSSED WITHOUT COLLISION.
C
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'COMUSR'
      INCLUDE 'CESTIM'
      INCLUDE 'CGRID'
      INCLUDE 'COMXS'
      INCLUDE 'CSPEZ'
      INCLUDE 'CSDVI'
      DIMENSION XSTOR(NSTOR),XSTOR2(NSTOR,N2ND+N3RD)
      EQUIVALENCE (XSTOR(1),SIGVCX(1))
C
C  ESTIMATORS FOR ATOMS
C
      ENTRY UPDATM (XSTOR2)
C
      WV=WEIGHT/VEL
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,WV)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,WV)
      IF (NPBGKA(IATM).GT.0) THEN
                      NPBGK=NPBGKA(IATM)
                      CALL UPTBGK(WV,NPBGK)
      ENDIF
C
      VELQ=VEL*VEL
C
      DO 51 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        IRD=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        EDENA(IATM,IRD)=EDENA(IATM,IRD)+WTRE0
        PDENA(IATM,IRD)=PDENA(IATM,IRD)+WTR
        LMETSP(IATM)=.TRUE.
C
C    ESTIMATORS FOR SOURCES AND SINKS
C    NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C    POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRD,0)) GOTO 51
C
        DO 32 K=1,NSTOR
          XSTOR(K)=XSTOR2(K,I)
32      CONTINUE
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        PAAT(IATM,IRD)=PAAT(IATM,IRD)-WTRSIG
        EAAT(IRD)     =EAAT(IRD)     -WTRSIG*E0
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGACX(IATM,0,0).EQ.0) GOTO 43
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 44  IACX=1,NACXI(IATM)
          IRCX=LGACX(IATM,IACX,0)
          IPLS=LGACX(IATM,IACX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
          IF (IESTCX(IRCX,1).NE.0) THEN
            PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
            LMETSP(NSPAMI+IPLS)=.TRUE.
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              PAAT(IAT1,IRD)= PAAT(IAT1,IRD)+WTRSIG
              LMETSP(IAT1)=.TRUE.
C           ELSEIF (N1STX(IRCX,1).EQ.2) THEN
C             IML1=N1STX(IRCX,2)
C             LOGMOL(IML1,ISTRA)=.TRUE.
C             PAML(IML1,IRD)= PAML(IML1,IRD)+WTRSIG
C             LMETSP(NATMI+IML1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              PAIO(IIO1,IRD)= PAIO(IIO1,IRD)+WTRSIG
              LMETSP(NSPAM+IIO1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              PAPL(IPL1,IRD)= PAPL(IPL1,IRD)+WTRSIG
              LMETSP(NSPAMI+IPL1)=.TRUE.
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              PAAT(IAT2,IRD)= PAAT(IAT2,IRD)+WTRSIG
              LMETSP(IAT2)=.TRUE.
C           ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
C             IML2=N2NDX(IRCX,2)
C             LOGMOL(IML2,ISTRA)=.TRUE.
C             PAML(IML2,IRD)= PAML(IML2,IRD)+WTRSIG
C             LMETSP(NATMI+IML2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              PAIO(IIO2,IRD)= PAIO(IIO2,IRD)+WTRSIG
              LMETSP(NSPAM+IIO2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              PAPL(IPL2,IRD)= PAPL(IPL2,IRD)+WTRSIG
              LMETSP(NSPAMI+IPL2)=.TRUE.
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTCX(IRCX,3).NE.0) THEN
            EAAT(IRD)=EAAT(IRD)          +WTRSIG*E0
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            EAPL(IRD)     =EAPL(IRD)     -WTRSIG*ESIGCX(IRCX,1)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              EAAT(IRD)     = EAAT(IRD)     +WTRSIG*ESIGCX(IRCX,1)
C           ELSEIF (N1STX(IRCX,1).EQ.2) THEN
C             IML1=N1STX(IRCX,2)
C             LOGMOL(IML1,ISTRA)=.TRUE.
C             EAML(IRD)     = EAML(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              EAIO(IRD)     = EAIO(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              EAPL(IRD)     = EAPL(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              EAAT(IRD)     = EAAT(IRD)     +WTRSIG*E0
C           ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
C             IML2=N2NDX(IRCX,2)
C             LOGMOL(IML2,ISTRA)=.TRUE.
C             EAML(IRD)     = EAML(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              EAIO(IRD)     = EAIO(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              EAPL(IRD)     = EAPL(IRD)     +WTRSIG*E0
            ENDIF
          ENDIF
C
44      CONTINUE
43      CONTINUE
C
C  ELASTIC NEUTRAL BULK-ION COLLISION CONTRIBUTION
C
        IF (LGAEL(IATM,0,0).EQ.0) GOTO 60
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 61  IAEL=1,NAELI(IATM)
          IREL=LGAEL(IATM,IAEL,0)
          IPLS=LGAEL(IATM,IAEL,1)
C  DO NOT UPDATE BGK TALLIES HERE
          IBGK=NPBGKP(IPLS,1)
          IF (IBGK.NE.0) GOTO 61
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVEL(IREL)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTEL(IREL,1).NE.0) THEN
            PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C  UPDATE TRACKLENGTH ESTIMATOR
C           PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
C           PAPL(IPLS,IRD)=PAPL(IPLS,IRD)+WTRSIG
C           LMETSP(NSPAMI+IPLS)=.TRUE.
            PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
            LMETSP(IATM)=.TRUE.
          ENDIF
C
          IF (IESTEL(IREL,3).NE.0) THEN
C  COLLISION ESTIMATOR IN SUBR. COLLIDE?
C  COMPENSATE PRE COLLISION RATES HERE
            EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            EAPL(IRD)=EAPL(IRD)-WTRSIG*ESIGEL(IREL,1)
C
C  FIRST SECONDARY: = INCIDENT ION. REMAINS SAME PARTICLE BY DEFAULT
            EAPL(IRD)=EAPL(IRD)+WTRSIG*E0
C  SECOND SECONDARY: = INCIDENT ATOM. REMAINS SAME PARTICLE BY DEFAULT
            EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGEL(IREL,1)
          ENDIF
C
61      CONTINUE
60      CONTINUE
C
C  ELECTRON IMPACT COLLISION CONTRIBUTION
C
        IF (LGAEI(IATM,0).EQ.0) GOTO 57
C
        DO 55 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          IF (SIGVEI(IRDS).LE.0.D0) GOTO 55
C
          WTRSIG=WTR*SIGVEI(IRDS)
C
C  COLLISION ESTIMATOR:
C  NOT AVAILABLE
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
          PAEL(IRD)=PAEL(IRD)+WTRSIG*PELDS(IRDS)
C
C  POST COLLISION CONTRIBUTIONS
          IF (PATDS(IRDS,0).GT.0.D0) THEN
            DO 52 IAT=1,NATMI
              P=PATDS(IRDS,IAT)
              IF (P.GT.0.D0) THEN
                LOGATM(IAT,ISTRA)=.TRUE.
                PAAT(IAT,IRD)=PAAT(IAT,IRD)+P*WTRSIG
                LMETSP(IAT)=.TRUE.
              ENDIF
52          CONTINUE
          ENDIF
          IF (PMLDS(IRDS,0).GT.0.D0) THEN
            DO 53 IML=1,NMOLI
              P=PMLDS(IRDS,IML)
              IF (P.GT.0.D0) THEN
                LOGMOL(IML,ISTRA)=.TRUE.
                PAML(IML,IRD)=PAML(IML,IRD)+P*WTRSIG
                LMETSP(NATMI+IML)=.TRUE.
              ENDIF
53          CONTINUE
          ENDIF
          IF (PIODS(IRDS,0).GT.0.D0) THEN
            DO 54 IIO=1,NIONI
              P=PIODS(IRDS,IIO)
              IF (P.GT.0.D0) THEN
                LOGION(IIO,ISTRA)=.TRUE.
                PAIO(IIO,IRD)=PAIO(IIO,IRD)+P*WTRSIG
                LMETSP(NSPAM+IIO)=.TRUE.
              ENDIF
54          CONTINUE
          ENDIF
          IF (PPLDS(IRDS,0).NE.0.D0) THEN
            DO 56 IPL=1,NPLSI
              P=PPLDS(IRDS,IPL)
              IF (P.NE.0.D0) THEN
                LOGPLS(IPL,ISTRA)=.TRUE.
                PAPL(IPL,IRD)=PAPL(IPL,IRD)+P*WTRSIG
                LMETSP(NSPAMI+IPL)=.TRUE.
              ENDIF
56          CONTINUE
          ENDIF
C
          IF (IESTEI(IRDS,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
            EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
C
          ELSE
C
            EAEL(IRD)=EAEL(IRD)+WTRSIG*ESIGEI(IRDS,0)
            EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGEI(IRDS,1)
            EAML(IRD)=EAML(IRD)+WTRSIG*ESIGEI(IRDS,2)
            EAIO(IRD)=EAIO(IRD)+WTRSIG*ESIGEI(IRDS,3)
            EAPL(IRD)=EAPL(IRD)+WTRSIG*ESIGEI(IRDS,4)
C
          ENDIF
55      CONTINUE
57      CONTINUE
C
C  PLASMA ION IMPACT CONTRIBUTION
C
        IF (LGAPI(IATM,0,0).EQ.0) GOTO 59
        DO 58  IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          IPLS=LGAPI(IATM,IAPI,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.

          WTRSIG=WTR*SIGVPI(IRPI)
C
C  PRE COLLISION BULK ION CONTRIBUTION, ASSUME: INCIDENT ION IS LOST
C
          PAPL(IPLS,IRD)=PAPL(IPLS,IRD)+WTRSIG*(-1.D0)
          LMETSP(NSPAMI+IPLS)=.TRUE.
          EAPL(IRD)     =EAPL(IRD)     +WTRSIG*ESIGPI(IRPI,1)
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION
C
          EAEL(IRD)=EAEL(IRD)+WTRSIG*ESIGPI(IRPI,0)
          PAEL(IRD)=PAEL(IRD)+WTRSIG*PELPI(IRPI)
C
          IF (PIOPI(IRPI,0).GT.0) THEN
C SO NICHT  EAIO(IRD)     = EAIO(IRD)     +WTRSIG*E0
            DO IIO=1,NIONI
              P=PIOPI(IRPI,IIO)
              IF (P.GT.0) THEN
                LOGION(IIO,ISTRA)=.TRUE.
                PAIO(IIO,IRD)= PAIO(IIO,IRD)+WTRSIG*P
                LMETSP(NSPAM+IIO)=.TRUE.
              ENDIF
            ENDDO
          ENDIF
          IF (PPLPI(IRPI,0).GT.0) THEN
C SO NICHT  EAPL(IRD)     = EAPL(IRD)     +WTRSIG*E0
            DO IPL=1,NPLSI
              P=PPLPI(IRPI,IPL)
              IF (P.GT.0) THEN
                LOGPLS(IPL,ISTRA)=.TRUE.
                PAPL(IPL,IRD)= PAPL(IPL,IRD)+WTRSIG*P
                LMETSP(NSPAMI+IPL)=.TRUE.
              ENDIF
            ENDDO
          ENDIF
58      CONTINUE
59      CONTINUE
51    CONTINUE
      RETURN
C
C  ESTIMATORS FOR MOLECULES
C
      ENTRY UPDMOL (XSTOR2)
C
      WV=WEIGHT/VEL
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,WV)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,WV)
      IF (NPBGKM(IMOL).GT.0) THEN
                      NPBGK=NPBGKM(IMOL)
                      CALL UPTBGK(WV,NPBGK)
      ENDIF
C
      VELQ=VEL*VEL
C
      DO 71 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        IRD=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        EDENM(IMOL,IRD)=EDENM(IMOL,IRD)+WTRE0
        PDENM(IMOL,IRD)=PDENM(IMOL,IRD)+WTR
        LMETSP(NATMI+IMOL)=.TRUE.
C
C    ESTIMATORS FOR SOURCES AND SINKS
C    NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C    POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRD,0)) GOTO 71
C
        DO 75 K=1,NSTOR
          XSTOR(K)=XSTOR2(K,I)
75      CONTINUE
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        PMML(IMOL,IRD)=PMML(IMOL,IRD)-WTRSIG
        EMML(IRD)     =EMML(IRD)     -WTRSIG*E0
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGMCX(IMOL,0,0).EQ.0) GOTO 79
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 76  IMCX=1,NMCXI(IMOL)
          IRCX=LGMCX(IMOL,IMCX,0)
          IPLS=LGMCX(IMOL,IMCX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTCX(IRCX,1).NE.0) THEN
            PMML(IMOL,IRD)=PMML(IMOL,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            PMPL(IPLS,IRD)=PMPL(IPLS,IRD)-WTRSIG
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              PMAT(IAT1,IRD)= PMAT(IAT1,IRD)+WTRSIG
              LMETSP(IAT1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              PMML(IML1,IRD)= PMML(IML1,IRD)+WTRSIG
              LMETSP(NATMI+IML1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              PMIO(IIO1,IRD)= PMIO(IIO1,IRD)+WTRSIG
              LMETSP(NSPAM+IIO1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              PMPL(IPL1,IRD)= PMPL(IPL1,IRD)+WTRSIG
              LMETSP(NSPAMI+IPL1)=.TRUE.
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              PMAT(IAT2,IRD)= PMAT(IAT2,IRD)+WTRSIG
              LMETSP(IAT2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              PMML(IML2,IRD)= PMML(IML2,IRD)+WTRSIG
              LMETSP(NATMI+IML2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              PMIO(IIO2,IRD)= PMIO(IIO2,IRD)+WTRSIG
              LMETSP(NSPAM+IIO2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              PMPL(IPL2,IRD)= PMPL(IPL2,IRD)+WTRSIG
              LMETSP(NSPAMI+IPL2)=.TRUE.
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
          IF (IESTCX(IRCX,3).NE.0) THEN
            EMML(IRD)=EMML(IRD)          +WTRSIG*E0
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            EMPL(IRD)     =EMPL(IRD)     -WTRSIG*ESIGCX(IRCX,1)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              EMAT(IRD)     = EMAT(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              EMML(IRD)     = EMML(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              EMIO(IRD)     = EMIO(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              EMPL(IRD)     = EMPL(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ENDIF
C  SECOND SECONDARY: PREVIOUS MOLECULE IMOL
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              EMAT(IRD)     = EMAT(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              EMML(IRD)     = EMML(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              EMIO(IRD)     = EMIO(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              EMPL(IRD)     = EMPL(IRD)     +WTRSIG*E0
            ENDIF
          ENDIF
C
76      CONTINUE
79      CONTINUE
C
C  ELASTIC NEUTRAL BULK-ION COLLISION CONTRIBUTION
C
        IF (LGMEL(IMOL,0,0).EQ.0) GOTO 80
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 81  IMEL=1,NMELI(IMOL)
          IREL=LGMEL(IMOL,IMEL,0)
          IPLS=LGMEL(IMOL,IMEL,1)
C  DO NOT UPDATE BGK TALLIES HERE
          IBGK=NPBGKP(IPLS,1)
          IF (IBGK.NE.0) GOTO 81
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVEL(IREL)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
          IF (IESTEL(IREL,1).NE.0) THEN
            PMML(IMOL,IRD)=PMML(IMOL,IRD)+WTRSIG
          ELSE
C  UPDATE TRACKLENGTH ESTIMATOR
C           PMPL(IPLS,IRD)=PMPL(IPLS,IRD)-WTRSIG
C           PMPL(IPLS,IRD)=PMPL(IPLS,IRD)+WTRSIG
C           LMETSP(NSPAMI+IPLS)=.TRUE.
            PMML(IMOL,IRD)=PMML(IMOL,IRD)+WTRSIG
            LMETSP(NATMI+IMOL)=.TRUE.
          ENDIF
          IF (IESTEL(IREL,3).NE.0) THEN
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
            EMML(IRD)=EMML(IRD)+WTRSIG*E0
          ELSE
            EMPL(IRD)=EMPL(IRD)-WTRSIG*ESIGEL(IREL,1)
            EMPL(IRD)=EMPL(IRD)+WTRSIG*E0
            EMML(IRD)=EMML(IRD)+WTRSIG*ESIGEL(IREL,1)
          ENDIF
C
81      CONTINUE
80      CONTINUE
C
C  ELECTRON IMPACT COLLISION CONTRIBUTION
C
        IF (LGMEI(IMOL,0).EQ.0) GOTO 100
C
        DO 90 IMEI=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMEI)
          IF (SIGVEI(IRDS).LE.0.D0) GOTO 90
C
          WTRSIG=WTR*SIGVEI(IRDS)
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
          PMEL(IRD)=PMEL(IRD)+WTRSIG*PELDS(IRDS)
C
C  POST COLLISION CONTRIBUTIONS
          IF (PATDS(IRDS,0).GT.0.D0) THEN
            DO 92 IAT=1,NATMI
              P=PATDS(IRDS,IAT)
              IF (P.GT.0.D0) THEN
                LOGATM(IAT,ISTRA)=.TRUE.
                PMAT(IAT,IRD)=PMAT(IAT,IRD)+P*WTRSIG
                LMETSP(IAT)=.TRUE.
              ENDIF
92          CONTINUE
          ENDIF
          IF (PMLDS(IRDS,0).GT.0.D0) THEN
            DO 93 IML=1,NMOLI
              P=PMLDS(IRDS,IML)
              IF (P.GT.0.D0) THEN
                LOGMOL(IML,ISTRA)=.TRUE.
                PMML(IML,IRD)=PMML(IML,IRD)+P*WTRSIG
                LMETSP(NATMI+IML)=.TRUE.
              ENDIF
93          CONTINUE
          ENDIF
          IF (PIODS(IRDS,0).GT.0.D0) THEN
            DO 94 IIO=1,NIONI
              P=PIODS(IRDS,IIO)
              IF (P.GT.0.D0) THEN
                LOGION(IIO,ISTRA)=.TRUE.
                PMIO(IIO,IRD)=PMIO(IIO,IRD)+P*WTRSIG
                LMETSP(NSPAM+IIO)=.TRUE.
              ENDIF
94          CONTINUE
          ENDIF
          IF (PPLDS(IRDS,0).NE.0.D0) THEN
            DO 96 IPL=1,NPLSI
              P=PPLDS(IRDS,IPL)
              IF (P.NE.0.D0) THEN
                LOGPLS(IPL,ISTRA)=.TRUE.
                PMPL(IPL,IRD)=PMPL(IPL,IRD)+P*WTRSIG
                LMETSP(NSPAMI+IPL)=.TRUE.
              ENDIF
96          CONTINUE
          ENDIF
C
          IF (IESTEI(IRDS,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
            EMML(IRD)=EMML(IRD)+WTRSIG*E0
C
          ELSE
C
            EMEL(IRD)=EMEL(IRD)+WTRSIG*ESIGEI(IRDS,0)
            EMAT(IRD)=EMAT(IRD)+WTRSIG*ESIGEI(IRDS,1)
            EMML(IRD)=EMML(IRD)+WTRSIG*ESIGEI(IRDS,2)
            EMIO(IRD)=EMIO(IRD)+WTRSIG*ESIGEI(IRDS,3)
            EMPL(IRD)=EMPL(IRD)+WTRSIG*ESIGEI(IRDS,4)
C
          ENDIF
90      CONTINUE
100     CONTINUE
C
71    CONTINUE
      RETURN
C
C
C  ESTIMATORS FOR TEST IONS
C
      ENTRY UPDION (XSTOR2)
C
      WV=WEIGHT/VEL
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,WV)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,WV)
      IF (NPBGKI(IION).GT.0) THEN
                      NPBGK=NPBGKI(IION)
                      CALL UPTBGK(WV,NPBGK)
      ENDIF
C
      VELQ=VEL*VEL
C
      DO 111 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        IRD=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        EDENI(IION,IRD)=EDENI(IION,IRD)+WTRE0
        PDENI(IION,IRD)=PDENI(IION,IRD)+WTR
        LMETSP(NSPAM+IION)=.TRUE.
C
C    ESTIMATORS FOR SOURCES AND SINKS
C    NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C    POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRD,0)) GOTO 111
C
        DO 115 K=1,NSTOR
          XSTOR(K)=XSTOR2(K,I)
115     CONTINUE
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        PIIO(IION,IRD)=PIIO(IION,IRD)-WTRSIG
        EIIO(IRD)     =EIIO(IRD)     -WTRSIG*E0
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGICX(IION,0,0).EQ.0) GOTO 119
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 116  IICX=1,NICXI(IION)
          IRCX=LGICX(IION,IICX,0)
          IPLS=LGICX(IION,IICX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
          IF (IESTCX(IRCX,1).NE.0) THEN
            PIIO(IION,IRD)=PIIO(IION,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            PIPL(IPLS,IRD)=PIPL(IPLS,IRD)-WTRSIG
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              PIAT(IAT1,IRD)= PIAT(IAT1,IRD)+WTRSIG
              LMETSP(IAT1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              PIML(IML1,IRD)= PIML(IML1,IRD)+WTRSIG
              LMETSP(NATMI+IML1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              PIIO(IIO1,IRD)= PIIO(IIO1,IRD)+WTRSIG
              LMETSP(NSPAM+IIO1)=.TRUE.
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              PIPL(IPL1,IRD)= PIPL(IPL1,IRD)+WTRSIG
              LMETSP(NSPAMI+IPL1)=.TRUE.
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              PIAT(IAT2,IRD)= PIAT(IAT2,IRD)+WTRSIG
              LMETSP(IAT2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              PIML(IML2,IRD)= PIML(IML2,IRD)+WTRSIG
              LMETSP(NATMI+IML2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              PIIO(IIO2,IRD)= PIIO(IIO2,IRD)+WTRSIG
              LMETSP(NSPAM+IIO2)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              PIPL(IPL2,IRD)= PIPL(IPL2,IRD)+WTRSIG
              LMETSP(NSPAMI+IPL2)=.TRUE.
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
          IF (IESTCX(IRCX,3).NE.0) THEN
            EIIO(IRD)=EIIO(IRD)          +WTRSIG*E0
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            EIPL(IRD)     =EIPL(IRD)     -WTRSIG*ESIGCX(IRCX,1)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              EIAT(IRD)     = EIAT(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              EIML(IRD)     = EIML(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              EIIO(IRD)     = EIIO(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              EIPL(IRD)     = EIPL(IRD)     +WTRSIG*ESIGCX(IRCX,1)
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              EIAT(IRD)     = EIAT(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              EIML(IRD)     = EIML(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              EIIO(IRD)     = EIIO(IRD)     +WTRSIG*E0
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              EIPL(IRD)     = EIPL(IRD)     +WTRSIG*E0
            ENDIF
          ENDIF
C
116     CONTINUE
119     CONTINUE
C
C  ELECTRON IMPACT COLLISION CONTRIBUTION
C
        IF (LGIEI(IION,0).EQ.0) GOTO 130
C
        DO 120 IIEI=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIEI)
          IF (SIGVEI(IRDS).LE.0.D0) GOTO 120
C
          WTRSIG=WTR*SIGVEI(IRDS)
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
          PIEL(IRD)=PIEL(IRD)+WTRSIG*PELDS(IRDS)
C
C  POST COLLISION CONTRIBUTIONS
          IF (PATDS(IRDS,0).GT.0.D0) THEN
            DO 122 IAT=1,NATMI
              P=PATDS(IRDS,IAT)
              IF (P.GT.0.D0) THEN
                LOGATM(IAT,ISTRA)=.TRUE.
                PIAT(IAT,IRD)=PIAT(IAT,IRD)+P*WTRSIG
                LMETSP(IAT)=.TRUE.
              ENDIF
122         CONTINUE
          ENDIF
          IF (PMLDS(IRDS,0).GT.0.D0) THEN
            DO 123 IML=1,NMOLI
              P=PMLDS(IRDS,IML)
              IF (P.GT.0.D0) THEN
                LOGMOL(IML,ISTRA)=.TRUE.
                PIML(IML,IRD)=PIML(IML,IRD)+P*WTRSIG
                LMETSP(NATMI+IML)=.TRUE.
              ENDIF
123         CONTINUE
          ENDIF
          IF (PIODS(IRDS,0).GT.0.D0) THEN
            DO 124 IIO=1,NIONI
              P=PIODS(IRDS,IIO)
              IF (P.GT.0.D0) THEN
                LOGION(IIO,ISTRA)=.TRUE.
                PIIO(IIO,IRD)=PIIO(IIO,IRD)+P*WTRSIG
                LMETSP(NSPAM+IIO)=.TRUE.
              ENDIF
124         CONTINUE
          ENDIF
          IF (PPLDS(IRDS,0).NE.0.D0) THEN
            DO 126 IPL=1,NPLSI
              P=PPLDS(IRDS,IPL)
              IF (P.NE.0.D0) THEN
                LOGPLS(IPL,ISTRA)=.TRUE.
                PIPL(IPL,IRD)=PIPL(IPL,IRD)+P*WTRSIG
                LMETSP(NSPAMI+IPL)=.TRUE.
              ENDIF
126         CONTINUE
          ENDIF
C
          IF (IESTEI(IRDS,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
            EIIO(IRD)=EIIO(IRD)+WTRSIG*E0
C
          ELSE
C
            EIEL(IRD)=EIEL(IRD)+WTRSIG*ESIGEI(IRDS,0)
            EIAT(IRD)=EIAT(IRD)+WTRSIG*ESIGEI(IRDS,1)
            EIML(IRD)=EIML(IRD)+WTRSIG*ESIGEI(IRDS,2)
            EIIO(IRD)=EIIO(IRD)+WTRSIG*ESIGEI(IRDS,3)
            EIPL(IRD)=EIPL(IRD)+WTRSIG*ESIGEI(IRDS,4)
C
          ENDIF
120     CONTINUE
130     CONTINUE
C
111   CONTINUE
      RETURN
C
      END
C
      SUBROUTINE LOCATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  LOCATE MONTE-CARLO PARTICLE
C
C  CALLED AT ENTRY LOCAT0 AT INITIALISATION FOR EACH STRATUM ISTRA
C     PURPOSE: PRECOMPUTING SOME QUANTITIES TO SPEED UP RANDOM SAMPLING
C              DURING PARTICLE TRACING
C
C  CALLED AT ENTRY LOCAT1 FOR EACH NEW LAUNCHED MONTE CARLO TRAJECTORY
C  FROM PARTICLE LOOP IN SUBR. MCARLO
C     PURPOSE: SET INITIAL TEST FLIGHT STATE, DEFINED BY THE VARIABLES
C              NO. 1 ... TO NPARTC+MPARTC OF COMMON BLOCK "COMPRT"
C              I.E.,
C                  X0... TO IUPDTE
C  CALLED PROGRAMS: SAMPNT (POINT SOURCE)
C                   SAMLNE (LINE SOURCE)  (NOT READY)
C                   SAMSRF (SURFACE SOURCE)
C                   SAMVOL (VOLUME SOURCE)
C  LOCAL VARIABLES: TEWL,TIWL(IPLS),DIWL(IPLS),
C                   VXWL(IPLS),VYWL(IPLS),VZWL(IPLS):
C
C                   THESE ARE BACKGROUND PARAMETERS USED FOR SAMPLING
C                   IN VELOCITY SPACE, IN CASE NLPLS, I.E., IF THE
C                   TEST FLIGHT STARTS AS BACKGROUND PARTICLE, THEN
C                   "RECOMBINING" INTO AT TEST PARTICLE
C                   EG. AT A SURFACE (NLSRF) OR IN THE VOLUME (NLVOL)
C                   IN THE OPPOSITE CASE (.NOT.NLPLS) PARAMETERS
C                   FOR THE SAMPLING DISTRIBUTION ARE SPECIFIED
C                   BY INPUT PARAMETERS IN BLOCK 7.,EG. SORENE,SORENI
C                   SORVDX,SORVDY,SORVDZ AND APPROPRIATE NEMOD2 AND
C                   NEMOD3 FLAGS
C
C                   WEISPZ(ISPZ):
C
C                   ANALOG SPECIES SAMPLING DISTRIBUTION
C                   SPECIES SAMPLING MAY ALSO BE DONE BY BIASED SOURCE
C                   SAMPLING, USING THE DATM,DMOL,DION OR DPLS DISTRIB.
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'COMNNL'
      INCLUDE 'CZT1'
      INCLUDE 'COMSOU'
      INCLUDE 'CTRCEI'
      INCLUDE 'COMXS'
      INCLUDE 'CESTIM'
      INCLUDE 'CLGIN'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CTRIG'
      INCLUDE 'CSPEZ'
      INCLUDE 'COMUSR'
      INCLUDE 'COUTAU'
      INCLUDE 'CLOGAU'
      INCLUDE 'CCONA'
      INCLUDE 'COMSPL'
      INCLUDE 'CPOLYG'
      DIMENSION RPST(NPARTC),IPST(MPARTC)
      DIMENSION RPSTT(NPARTT),IPSTT(MPARTT)
      DIMENSION DUMT(3),DUMV(3)
      DIMENSION WMM(NSRFS),IICSOR(NSRFS),
     .          ITISOR(NSRFS),IUPSOR(NSRFS),IFPSOR(NSRFS),
     .          VXWL(NPLS),VYWL(NPLS),VZWL(NPLS),VPWL(NPLS),
     .          TIWL(NPLS),DIWL(NPLS),
     .          WEISPZ(NSPZ)
      EQUIVALENCE (RPST(1),X0),(IPST(1),IPOLG)
      EQUIVALENCE (RPSTT(1),X0),(IPSTT(1),NPANU)
      LOGICAL NLSPUT,NLTST
      SAVE
C
      ENTRY LOCAT0
C
C  PREPARE DATA FOR SAMPLING SUBSTRATA FOR STRATUM ISTRA: 1--10
C
      DO 1 ISPZ=1,NSPZ
        WEISPZ(ISPZ)=-1.
1     CONTINUE
C
      SUM=0.
      DO 2 ISRFS=1,NSRFSI(ISTRA)
2       SUM=SUM+SORWGT(ISRFS,ISTRA)
c slmod begin - bug - not tr (fixed)
c...Redundant:
c      DO ISRFS=1,NSRFSI(ISTRA)
c        SUM=SUM+SORWGT(ISRFS,ISTRA)
c      ENDDO
c slmod end
      IF (SUM.LE.0.D0) THEN
        WRITE (6,*) 'NO SOURCE MODEL FOR STRATUM NO'
        WRITE (6,*) 'ISTRA=',ISTRA,' BECAUSE THE SUM OF THE FLUXES'
        WRITE (6,*) 'FROM THE SUBSTRATA DEFINED BY'
        WRITE (6,*) 'SORWGT(SUBSTRATUM,STRATUM) IS .LE. ZERO'
        WRITE (6,*) 'THIS STRATUM IS TURNED OF !!'
        NPTS(ISTRA)=0
        RETURN
      ENDIF
      SUM1=0.
      NLIMSQ=NSRFSI(ISTRA)
      DO 4 ISOUR=1,NSRFSI(ISTRA)
        SUM1=SUM1+SORWGT(ISOUR,ISTRA)
        WMM(ISOUR)=SUM1/SUM
4     CONTINUE
C
C  PREPARE SOME DATA FOR ENERGY SAMPLING AND HISTORY INITIALIZATION
C
      NEMOD1=IDEZ(NEMODS(ISTRA),1,4)
      NEMOD2=IDEZ(NEMODS(ISTRA),2,4)
      NEMOD3=IDEZ(NEMODS(ISTRA),3,4)
      NEMDSP=IDEZ(NEMODS(ISTRA),4,4)
C
      DO 5 ISRFS=1,NSRFSI(ISTRA)
        IF (SORIFL(ISRFS,ISTRA).NE.0) THEN
          IDUMM=SORIFL(ISRFS,ISTRA)
          ITISOR(ISRFS)=IDEZ(IDUMM,1,4)
          IF (ITISOR(ISRFS).EQ.2) ITISOR(ISRFS)=-1
          IFPSOR(ISRFS)=IDEZ(IDUMM,2,4)
          IF (IFPSOR(ISRFS).EQ.2) IFPSOR(ISRFS)=-1
          IUPSOR(ISRFS)=IDEZ(IDUMM,3,4)
          IF (IUPSOR(ISRFS).EQ.2) IUPSOR(ISRFS)=-1
          IICSOR(ISRFS)=IDEZ(IDUMM,4,4)
          IF (IICSOR(ISRFS).EQ.2) IICSOR(ISRFS)=-1
        ELSE
          ITISOR(ISRFS)=0
          IFPSOR(ISRFS)=0
          IUPSOR(ISRFS)=0
          IICSOR(ISRFS)=0
        ENDIF
5     CONTINUE
C
      SNORM=SQRT(SORCTX(ISTRA)**2+SORCTY(ISTRA)**2+SORCTZ(ISTRA)**2)
      IF (SNORM.GT.EPS10) THEN
        SORCTX(ISTRA)=SORCTX(ISTRA)/SNORM
        SORCTY(ISTRA)=SORCTY(ISTRA)/SNORM
        SORCTZ(ISTRA)=SORCTZ(ISTRA)/SNORM
      ENDIF
C
      IF (TRCSOU) THEN
        WRITE (6,*) 'NEMOD1,NEMOD2,NEMOD3 ',NEMOD1,NEMOD2,NEMOD3
        WRITE (6,*) 'SNORM  ',SNORM
        WRITE (6,*) 'ISRFS,IICSOR(I),ITISOR(I),IFPSOR(I),IUPSOR(I)'
        DO 6 I=1,NSRFSI(ISTRA)
          WRITE (6,*) I,IICSOR(I),ITISOR(I),IFPSOR(I),IUPSOR(I)
6       CONTINUE
      ENDIF
C
C  PREPARE SOME DATA FOR SPECIES SAMPLING
C

      RETURN
C
      ENTRY LOCAT1(IPANU)
C
C  TENTATIVELY ASSUME: A NEXT GENERATION PARTICLE WILL BE BORN
      LGPART=.TRUE.
C
C   SET SOME DEFAULT DATA TO INITIALIZE THIS HISTORY
C
      WEIGHT=1.0
      IATM=0
      IMOL=0
      IION=0
      IPLS=0
C
      ITIME=1
      IFPATH=1
      IUPDTE=1
C
      NCELL=0
      NBLOCK=1
      NACELL=0
      NBLCKA=0
      NRCELL=0
      NPCELL=1
      NTCELL=1
      IPOLG=1
      IPOLGN=1
      ICOL=0
C
C  DETAILED PRINTOUT OF TRAJECTORY FOR THIS PARTICLE?
C
      NLTRC=NPANU.GE.I1TRC.AND.NPANU.LE.I2TRC
C
C  =====================================================
C  =SAMPLE STARTING POINT FOR  ATOMS, MOLECULES OR IONS=
C  =====================================================
C
      LGTIME=NPRNLI.GT.0
C  DISTANCE TO "TIME-SURFACE"
      IF (.NOT.LGTIME) THEN
        DTIMVI=1.D30
      ELSEIF (LGTIME) THEN
        DTIMVI=TIME0+DTIMV
      ENDIF
C
C   SOURCE DUE TO TIME DEP. MODE, READ PARTICLES FROM CENSUS: RPARTC,IPARTC
      IF (NTIME.GT.0) WRITE(0,*) 'NLCNS',NLCNS(ISTRA),ISTRA,NSTRAI
      IF (NLCNS(ISTRA).AND.ISTRA.EQ.NSTRAI) THEN
C   LABELS  11---20
C   AT PRESENT: ONLY ONE SUBSTRATUM
        ISECT=1
        NLSTOR=IPANU.LE.ISTOR(ISECT,ISTRA)
C
C   BINARY SEARCH IN RPARTW ARRAY
c slmod begin - new
c
c This binary search seems to suck.  Why is it done this way?  Why not just 
c assign an integer from the range of A?      
c
c slmod end 
        A=RANF()*RPARTW(IPRNL)
        I1=0
        I2=IPRNL
9       IM=(I1+I2)/2
        IF(A.LT.RPARTW(IM)) THEN
          I2=IM
          GOTO 9
        ELSEIF(A.GT.RPARTW(IM+1)) THEN
          I1=IM
          GOTO 9
        ENDIF
C
c slmod being - new - bug?
c...ISTRA gets over-written...hard to believe this never got noticed,
c   so likely it appeared after the last use of this code
c
        ISTRA2=ISTRA
c slmod end
        IPTSI=IM+1
        WRITE(0,'(A,3I10)') 
     .    'TIMES:',IPTSI,0,IPRNL
        DO 11 J=1,NPARTT
          RPSTT(J)=RPARTC(IPTSI,J)
11      CONTINUE
        NPANUO=NPANU
        DO 12 J=1,MPARTT
          IPSTT(J)=IPARTC(IPTSI,J)
12      CONTINUE
c slmod begin - new - bug?
        ISTRA=ISTRA2
c slmod end
        ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
        ITYP=ISPEZI(ISPZ,0)
        IATM=ISPEZI(ISPZ,1)
        IMOL=ISPEZI(ISPZ,2)
        IION=ISPEZI(ISPZ,3)
        IPLS=ISPEZI(ISPZ,4)
        IF (output)
     .  WRITE(6,*) 'MARK: LOCATE: SETTING ITYP?  ',
     .             ITYP,IATM,IMOL,IION,IPLS
        NPANU=NPANUO
        NLSRFX=.FALSE.
        NLSRFY=.FALSE.
        NLSRFZ=.FALSE.
c slmod begin - bug? - not tr (fixed)
c...NSTST not declared anywhere? Likely should be NSTSI.
        MSURF=NLIM+NSTSI
c
c        MSURF=NLIM+NSTST
c slmod end
C
        WEIGHT=1.D0
C
        IF (ITYP.EQ.1) THEN
          WTOTA(IATM,ISTRA)=WTOTA(IATM,ISTRA)+WEIGHT
          ETOTA(ISTRA)=ETOTA(ISTRA)+E0*WEIGHT
        ELSEIF (ITYP.EQ.2) THEN
          WTOTM(IMOL,ISTRA)=WTOTM(IMOL,ISTRA)+WEIGHT
          ETOTM(ISTRA)=ETOTM(ISTRA)+E0*WEIGHT
        ELSEIF (ITYP.EQ.3) THEN
          WTOTI(IION,ISTRA)=WTOTI(IION,ISTRA)+WEIGHT
          ETOTI(ISTRA)=ETOTI(ISTRA)+E0*WEIGHT
        ELSE
          WRITE (6,*) 'ERROR IN LOCATE, CALL EXIT '
          CALL EXIT
        ENDIF
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,1)
        IF (NLSTOR) CALL STORE(1)
C
        GOTO 5000
C
C  POINT SOURCE MODEL  21---30
C
      ELSEIF (NLPNT(ISTRA)) THEN
C
C   FIRSTLY FIND POINT NUMBER IPOINT
        IPOINT=1
        IF (NLIMSQ.GT.1) THEN
          ZV=RANF( )
          DO 21 IPOINT=1,NLIMSQ
            IF (ZV.LT.WMM(IPOINT)) GOTO 22
21        CONTINUE
22        CONTINUE
        ENDIF
        ISECT=IPOINT
        NLSTOR=IPANU.LE.ISTOR(ISECT,ISTRA)
C
C   NEXT FIND CO-ORDINATES AND CELL INDICES,
C   LOCAL BACKGROUND TEMPERATURES TIWL AND TEWL, AND
C   LOCAL PLASMA DRIFT VELOCITIES VXWL,VYWL,VZWL FOR EACH BULK
C   ION SPECIES IPLS=1,NPLSI
C
C   NLPT=POINT INDEX IN (NSRFS) SOURCE ARRAYS
        CALL SAMPNT (IPOINT,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
        IF (.NOT.LGPART) RETURN
C
        IF (ITISOR(IPOINT).NE.0) THEN
          ITIME=ITISOR(IPOINT)
        ENDIF
        IF (IFPSOR(IPOINT).NE.0) THEN
          IFPATH=IFPSOR(IPOINT)
        ENDIF
        IF (IUPSOR(IPOINT).NE.0) THEN
          IUPATH=IUPSOR(IPOINT)
        ENDIF
        MSURF=0
C
C   LINE SOURCE  31---50
C
      ELSEIF (NLLNE(ISTRA)) THEN
        ILINE=1
        ISECT=ILINE
        MSURF=0
        WRITE (6,*) 'LINE SOURCE OPTION STILL TO BE WRITTEN. EXIT'
        CALL EXIT
C
C   SURFACE SOURCE MODEL  51---70
C
      ELSEIF (NLSRF(ISTRA)) THEN
        IF (output) THEN
          WRITE(6,*) 'MARK: '
          WRITE(6,*) 'MARK: LOCATE: LAUNCHING SURFACE : NPANU= ',NPANU
        ENDIF
C
C   FIRST FIND SOURCE-SURFACE NUMBER ISURF
        ISURF=1
        IF (NLIMSQ.GT.1) THEN
          ZV=RANF( )
          DO 51 ISURF=1,NLIMSQ
            IF (ZV.LT.WMM(ISURF)) GOTO 52
51        CONTINUE
52        CONTINUE
        ENDIF
        ISECT=ISURF
        NLSTOR=IPANU.LE.ISTOR(ISECT,ISTRA)
C
C   NEXT FIND POSITION ON THIS SOURCE SURFACE, AS WELL AS
C   CELL INDICES, LOCAL TEMPERATURES TIWL AND TEWL, AND
C   LOCAL PLASMA DRIFT VELOCITIES VXWL,VYWL,VZWL FOR EACH BULK
C   ION SPECIES IPLS=1,NPLSI
C
        CALL SAMSF1 (ISURF,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
        IF (.NOT.LGPART) RETURN
C
C   MSURF: NUMBER OF NON-DEFAULT (OR ADDITIONAL) SURFACE
C   MSURF=0 MEANS: SOURCE NOT ON ANY KNOWN SURFACE. DEFAULT REFLECTION MODEL
        MSURF=0
        IF (LEVGEO.EQ.4) THEN
          MSURF=ABS(INMTI(IPOLG,NRCELL))
        ELSE
          IF (MASURF.GT.0) THEN
            MSURF=MASURF
            ITRSF=0
          ELSEIF (MRSURF.GT.0) THEN
            ITRSF=INMP1I(MRSURF,NPCELL,NTCELL)
          ELSEIF (MPSURF.GT.0) THEN
            ITRSF=INMP2I(NRCELL,MPSURF,NTCELL)
          ELSEIF (MTSURF.GT.0) THEN
            ITRSF=INMP3I(NRCELL,NPCELL,MTSURF)
          ENDIF
          IF (ITRSF.GT.0) MSURF=NLIM+ITRSF
        ENDIF
C
C  SET ICOS AND SCOS SUCH AS IF THE SOURCE PARTICLE HAD ARRIVED
C  AT THE SURFACE FROM THE CORRECT SIDE AND IS NOW REFLECTED
C  (NOTE: THE FLAG "IWEI" USED IN SUBR. STDCOL AND ADDCOL
C  WILL ALWAYS BE POSITIVE WITH THIS DEFINITION OF SCOS)
C  THIS DEFAULT SETTING MAY BE OVERRULED BY SORIFL FLAG
C
        IF (IICSOR(ISURF).NE.0) THEN
          ICOS=IICSOR(ISURF)
        ELSEIF (ILSIDE(MSURF).NE.0) THEN
          ICOS=ISIGN(1,ILSIDE(MSURF))
        ELSE
          GOTO 990
        ENDIF
C
        SCOS=ICOS
C
        IF (ITISOR(ISURF).NE.0) THEN
          ITIME=ITISOR(ISURF)
        ELSEIF (ISWICH(1,MSURF).NE.0) THEN
          ITIME=ISWICH(1,MSURF)*ICOS
        ENDIF
        IF (IFPSOR(ISURF).NE.0) THEN
          IFPATH=IFPSOR(ISURF)
        ELSEIF (ISWICH(2,MSURF).NE.0) THEN
          IFPATH=ISWICH(2,MSURF)*ICOS
        ENDIF
        IF (IUPSOR(ISURF).NE.0) THEN
          IUPDTE=IUPSOR(ISURF)
        ELSEIF (ISWICH(3,MSURF).NE.0) THEN
          IUPDTE=ISWICH(3,MSURF)*ICOS
        ENDIF
C
C  FIND SURFACE NORMAL AT PLACE OF BIRTH
C
        IF (INDIM(ISURF,ISTRA).EQ.0) THEN
          CALL ADDNOR(X0,Y0,Z0,SCOS,MSURF,*55,*55)
        ELSEIF (INDIM(ISURF,ISTRA).GT.0) THEN
          CALL STDNOR (X0,Y0,Z0,INDIM(ISURF,ISTRA),SCOS,MSURF,*55,*55)
        ENDIF
55      CONTINUE
C
C  VOLUME SOURCE MODEL  71---90
C
      ELSEIF (NLVOL(ISTRA)) THEN
        IF (output) THEN
          WRITE(6,*) 'MARK: '
          WRITE(6,*) 'MARK: LOCATE: LAUNCHING VOLUME : NPANU= ',NPANU
        ENDIF
C  SUBSTRATA OF VOLUME SOURCE: IVOLM
        IVOLM=1
        IF (NLIMSQ.GT.1) THEN
          ZV=RANF( )
          DO 71 IVOLM=1,NLIMSQ
            IF (ZV.LT.WMM(IVOLM)) GOTO 72
71        CONTINUE
72        CONTINUE
        ENDIF
        ISECT=IVOLM
        CALL SAMVL1(IVOLM,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
        IF (.NOT.LGPART) RETURN
        NLSTOR=IPANU.LE.ISTOR(ISECT,ISTRA)
        MSURF=0
      ENDIF
C
      IRCELL=NRCELL
      IPCELL=NPCELL
      ITCELL=NTCELL
C
C  SAMPLE STARTING TIME
C
      ITMSTP=1
      IF (.NOT.LGTIME) THEN
        TIME=0.
      ELSEIF (LGTIME) THEN
        ISOR=SORLIM(ISECT,ISTRA)
        INDTEC=IDEZ(ISOR,4,4)
        IF (INDTEC.EQ.0) INDTEC=2
        IF (INDTEC.LE.1) TIME=TIME0
        IF (INDTEC.EQ.2) TIME=TIME0+RANF()*DTIMV
      ENDIF
C
C  INITIAL POSITION OF PARTICLE IS DEFINED NOW, FURTHERMORE:
C    NRCELL,NPCELL,NTCELL,IPOLG,NBLOCK,NACELL,
C    AND THE LOCAL BACKGROUND PARAMETERS
C    TEWL,(TIWL(IPLS),VXWL(IPLS),VYWL(IPLS),VZWL(IPLS),IPLS=1,NPLSI)
C
C    PLUS: WEISPZ FOR SOURCE SPECIES SAMPLING
C          WEISPZ IS THE ANALOG SAMPLING DISTRIBUTION
C          DPLS,DATM,DMOL,DION ARE THE NONANALOG SAMPLING DISTRIBUTIONS
C
C    PLUS: CRTX,CRTY,CRTZ,SCOS
C
C .........................................................................
C
C  FIND TYPE AND SPECIES INDEX AND RELATED CONSTANTS 100---199
C .........................................................................
C
      IF (NLATM(ISTRA)) THEN
        IF (output)
     .  WRITE(6,*) 'MARK: LOCATE: ITYP=1, WEIGHT=',WEIGHT
        ITYP=1
        IF (NSPEZ(ISTRA).LT.0) THEN
C  CHECK RADON-NIKODYM CONDITION FOR NON-ANALOG SAMPLING
          DO IATM=1,NATMI
            IF (DATD(IATM).LE.0.D0.AND.WEISPZ(IATM).GT.0.D0) THEN
              GOTO 992
            ENDIF
          ENDDO
        ENDIF
C  FIXED SPECIES INDEX
        IATM=NSPEZ(ISTRA)
        IF (IATM.LT.0.OR.IATM.GT.NATMI) THEN
C  SPECIES SAMPLING FROM DATM
          FR=RANF( )
          DO 102 I=1,NATMIM
            IATM=I
            IF (FR.LE.DATM(IATM)) GOTO 101
102       CONTINUE
          IATM=NATMI
101       CONTINUE
          IF (NSPEZ(ISTRA).LT.0) THEN
C  WEIGHT CORRECTION
            DAT=DATD(IATM)
            IF (WEISPZ(IATM).LT.0.D0) GOTO 999
            WEIGHT=WEIGHT*WEISPZ(IATM)/DAT
          ENDIF
        ELSEIF (IATM.EQ.0) THEN
C  ANALOG SPECIES SAMPLING FROM WEISPZ
          FR=RANF( )
          SUM=0.
          DO 112 I=1,NATMIM
            IATM=I
            IF (WEISPZ(IATM).LT.0.D0) GOTO 999
            SUM=SUM+WEISPZ(IATM)
            IF (FR.LE.SUM) GOTO 111
112       CONTINUE
          IATM=NATMI
111       CONTINUE
        ENDIF
        RSQDV=RSQDVA(IATM)*SQ2I
      ELSEIF (NLMOL(ISTRA)) THEN
        ITYP=2
        IF (NSPEZ(ISTRA).LT.0) THEN
C  CHECK RADON-NIKODYM CONDITION FOR NON-ANALOG SAMPLING
          DO IMOL=1,NMOLI
            IF (DMLD(IMOL).LE.0.D0.AND.WEISPZ(IMOL).GT.0.D0) THEN
              GOTO 992
            ENDIF
          ENDDO
        ENDIF
C  FIXED SPECIES INDEX
        IMOL=NSPEZ(ISTRA)
        IF (IMOL.LT.0.OR.IMOL.GT.NMOLI) THEN
C  NONANALOG SPECIES SAMPLING
          FR=RANF( )
          DO 104 I=1,NMOLIM
            IMOL=I
            IF (FR.LE.DMOL(IMOL)) GOTO 103
104       CONTINUE
          IMOL=NMOLI
103       CONTINUE
C  WEIGHT CORRECTION
          IF (NSPEZ(ISTRA).LT.0) THEN
            DML=DMLD(IMOL)
            IF (WEISPZ(IMOL).LT.0.D0) GOTO 999
            WEIGHT=WEIGHT*WEISPZ(IMOL)/DML
          ENDIF
        ELSEIF (IMOL.EQ.0) THEN
C  ANALOG SPECIES SAMPLING
          FR=RANF( )
          SUM=0.
          DO 114 I=1,NMOLIM
            IMOL=I
            IF (WEISPZ(IMOL).LT.0.D0) GOTO 999
            SUM=SUM+WEISPZ(IMOL)
            IF (FR.LE.SUM) GOTO 113
114       CONTINUE
          IMOL=NMOLI
113       CONTINUE
        ENDIF
        RSQDV=RSQDVM(IMOL)*SQ2I
      ELSEIF (NLION(ISTRA)) THEN
        ITYP=3
        IF (NSPEZ(ISTRA).LT.0) THEN
C  CHECK RADON-NIKODYM CONDITION FOR NON-ANALOG SAMPLING
          DO IION=1,NIONI
            IF (DIOD(IION).LE.0.D0.AND.WEISPZ(IION).GT.0.D0) THEN
              GOTO 992
            ENDIF
          ENDDO
        ENDIF
C  FIXED SPECIES INDEX
        IION=NSPEZ(ISTRA)
        IF (IION.LT.0.OR.IION.GT.NIONI) THEN
C  NONANALOG SPECIES SAMPLING
          FR=RANF( )
          DO 106 I=1,NIONIM
            IION=I
            IF (FR.LE.DION(IION)) GOTO 105
106       CONTINUE
          IION=NIONI
105       CONTINUE
C  WEIGHT CORRECTION
          IF (NSPEZ(ISTRA).LT.0) THEN
            DIO=DIOD(IION)
            IF (WEISPZ(IION).LT.0.D0) GOTO 999
            WEIGHT=WEIGHT*WEISPZ(IION)/DIO
          ENDIF
        ELSEIF (IION.EQ.0) THEN
C  ANALOG SPECIES SAMPLING
          FR=RANF( )
          SUM=0.
          DO 116 I=1,NIONIM
            IION=I
            IF (WEISPZ(IION).LT.0.D0) GOTO 999
            SUM=SUM+WEISPZ(IION)
            IF (FR.LE.SUM) GOTO 115
116       CONTINUE
          IION=NIONI
115       CONTINUE
        ENDIF
        RSQDV=RSQDVI(IION)*SQ2I
      ELSEIF (NLPLS(ISTRA)) THEN
        IF (output)
     .  WRITE(6,*) 'MARK: LOCATE: ITYP=4, WEIGHT=',WEIGHT
        ITYP=4
        IF (NSPEZ(ISTRA).LT.0) THEN
C  CHECK RADON-NIKODYM CONDITION FOR NON-ANALOG SAMPLING
          DO IPLS=1,NPLSI
            IF (DPLD(IPLS).LE.0.D0.AND.WEISPZ(IPLS).GT.0.D0) THEN
              GOTO 992
            ENDIF
          ENDDO
        ENDIF
C  FIXED SPECIES INDEX
        IPLS=NSPEZ(ISTRA)
        IF (IPLS.LT.0.OR.IPLS.GT.NPLSI) THEN
C  NONANALOG SPECIES SAMPLING
          IF (output)
     .    WRITE(6,*) 'MARK: LOCATE: NONANALOGUE SAMPLING  WEIGHT= ',
     .               WEIGHT
          FR=RANF( )
          IF (output)
     .    WRITE(6,*) 'MARK: LOCATE: NPLSIM= ',NPLSIM
          DO 108 I=1,NPLSIM
            IPLS=I
            IF (output)
     .      WRITE(6,'(A,F8.4,I4,F8.4)')
     .        ' MARK: LOCATE: FR,IPLS,DPLS(IPLS)= ',
     .        FR,IPLS,DPLS(IPLS)
            IF (FR.LE.DPLS(IPLS)) GOTO 107
108       CONTINUE
          IPLS=NPLSI
107       CONTINUE
C  WEIGHT CORRECTION
          IF (NSPEZ(ISTRA).LT.0) THEN
            DPL=DPLD(IPLS)
            IF (WEISPZ(IPLS).LT.0.D0) GOTO 999
            WEIGHT=WEIGHT*WEISPZ(IPLS)/DPL
            IF (output)
     .      WRITE(6,'(A,3I6,1P,3E10.2,0P)')
     .        ' MARK: LOCATE: NSPEZ= ',
     .        istra,nspez(istra),ipls,
     .        dpld(ipls),weispz(ipls),weight
          ENDIF
        ELSEIF (IPLS.EQ.0) THEN
C  ANALOG SPECIES SAMPLING
          FR=RANF( )
          SUM=0.
          DO 118 I=1,NPLSIM
            IPLS=I
            IF (WEISPZ(IPLS).LT.0.D0) GOTO 999
            SUM=SUM+WEISPZ(IPLS)
            IF (FR.LE.SUM) GOTO 117
118       CONTINUE
          IPLS=NPLSI
117       CONTINUE
        ENDIF
        RSQDV=RSQDVP(IPLS)*SQ2I
      ENDIF
C
      ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
      IF (output)
     .WRITE(6,*) 'MARK: LOCATE: SPECIES SAMPLING DONE  ITYP= ',ITYP
C  .............................................................
C
C  SPECIES SAMPLING DONE
C  .............................................................
C
C  MAKE SURE NOT TO WASTE TIME IN PARTICLES WITH ZERO WEIGHT
C
      IF (output)
     .WRITE(6,*) 'MARK: LOCATE: WEIGHT= ',WEIGHT
      LGPART=WEIGHT.GT.0.D0
      IF (.NOT.LGPART) RETURN
C
C  PARAMETERS FOR VELOCITY SAMPLING DISTRIBUTION:
C  TEWD,TIWD,VXWD,VYWD,VZDW
C
      IF (NEMOD2.EQ.1) THEN
C  SET SAMPLING TEMPERATURES FROM FIXED INPUT DATA
        TIWD=ABS(SORENI(ISTRA))
        TEWD=ABS(SORENE(ISTRA))
      ELSEIF (NEMOD2.EQ.2) THEN
C  NOT IN USE
      ELSEIF (NEMOD2.EQ.3) THEN
C  SET SAMPLING TEMPERATURES FROM LOCAL PLASMA DATA FOR SPECIES IPLTI
        IPLTI=NEMDSP
        IF (IPLTI.LT.1.OR.IPLTI.GT.NPLSI) GOTO 999
        TIWD=TIWL(IPLTI)
        TEWD=TEWL
      ELSE
C  DEFAULT: ONLY FOR NLPLS=TRUE:
C  SET SAMPLING TEMPERATURES FROM LOCAL PLASMA DATA FOR SPECIES IPLS
        TEWD=TEWL
        IF (NLPLS(ISTRA)) THEN
          IPL=IPLS
          TIWD=TIWL(IPL)
        ELSE
C  SET SAMPLING ION-TEMPERATURE TO ZERO
          TIWD=0.
        ENDIF
      ENDIF
C
      IF (NEMOD3.EQ.1) THEN
C  SET SAMPLING DRIFT VELOCITIES FROM INPUT DATA FOR DRIFT VELOCITY
        VXWD=SORVDX(ISTRA)
        VYWD=SORVDY(ISTRA)
        VZWD=SORVDZ(ISTRA)
      ELSEIF (NEMOD3.EQ.2) THEN
C  SET SAMPLING DRIFT VELOCITIES FROM INPUT DATA FOR MACH NUMBER
        CS=SQRT(1.*TIWD+TEWD)*RSQDV
        VXWD=SORVDX(ISTRA)*CS
        VYWD=SORVDY(ISTRA)*CS
        VZWD=SORVDZ(ISTRA)*CS
      ELSEIF (NEMOD3.EQ.3) THEN
        IPLV=NEMDSP
        IF (IPLV.LT.1.OR.IPLV.GT.NPLSI) GOTO 999
        VXWD=VXWL(IPLV)
        VYWD=VYWL(IPLV)
        VZWD=VZWL(IPLV)
      ELSE
C  DEFAULT: ONLY FOR NLPLS=TRUE:
C  SET SAMPLING DRIFT VELOCITIES FROM BACKGROUND DATA FOR SPECIES IPL
        IF (NLPLS(ISTRA)) THEN
          IPL=IPLS
          VXWD=VXWL(IPL)
          VYWD=VYWL(IPL)
          VZWD=VZWL(IPL)
        ELSE
          VXWD=0.
          VYWD=0.
          VZWD=0.
        ENDIF
      ENDIF
      IF (output)
     .WRITE(6,*) 'MARK: LOCATE: FINDING VELOCITY VECTOR  ITYP= ',ITYP
C
C  .....................................
C
C  FIND VELOCITY VECTOR NEXT
C  .....................................
C
C  PURELY ATOMIC SOURCE?  200 --- 299
C
      IF (NLATM(ISTRA)) THEN
        IF (NEMOD1.EQ.1) THEN
          EMAX=SORENI(ISTRA)
        ELSEIF (NEMOD1.EQ.6) THEN
          EMAX=0.
        ELSE
          GOTO 998
        ENDIF
        IF (EMAX.GT.0) THEN
          E0=EMAX
          VEL=SQRT(E0)*RSQDVA(IATM)
C
C  COSINE LIKE OR GAUSSIAN ANGLE DISTRIBUTION
C
C  IN CASE (CRTX,CRTY,CRTZ) NE (0.,0.,0.)
C  ASSUME NORMAL INCIDENCE AND USE REFLECTION MODEL ANGULAR DISTRIBUTION
          VELX=CRTX
          VELY=CRTY
          VELZ=CRTZ
          CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
C         VEL_MEAN=VEL
C         E0_MEAN=E0
        ELSEIF (EMAX.LE.0..AND.TIWD.GT.0..AND..NOT.NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM SHIFTED TRUNCATED MAXWELLIAN FLUX
C              AROUND INNER (!) NORMAL AT TEMP. TW (EV) = TIWD
          VWD=SQRT(VXWD**2+VYWD**2+VZWD**2)
          CALL VELOCS (TIWD,0.D0,VWD,VXWD,VYWD,VZWD,RSQDVA(IATM),
     .                 CVRSSA(IATM),
     .                 -CRTX,-CRTY,-CRTZ,E0,VELX,VELY,VELZ,VEL)
C  MODIFY ANGULAR DISTRIBUTION IN CASE SORCOS .NE. 0.5 (I.E., IN CASE
C  A NON-COSINE DISTRIBUTION IS REQUESTED
          IF (ABS(SORCOS(ISTRA)-0.5).GT.1.D-5) THEN
            VELX=CRTX
            VELY=CRTY
            VELZ=CRTZ
            CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                  SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
C           VEL_MEAN=VEL
C           E0_MEAN=E0
          ENDIF
        ELSEIF (EMAX.LE.0..AND.TIWD.GT.0..AND.NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM MAXWELLIAN AT TEMP. TW (EV) =TIWD
C
          NFLAG=2
          IDUM=1
          DUMT(1)=SQRT(TIWD/RMASSA(IATM))*CVEL2A
          DUMT(2)=DUMT(1)
          DUMT(3)=DUMT(1)
          DUMV(1)=0
          DUMV(2)=0
          DUMV(3)=0
          CALL VELOCX(0,VXO,VYO,VZO,VO,IO,NO,VELQ,NFLAG,
     .                IDUM,DUMT,DUMV)
C         E0=VELQ*CVRSSA(IATM)
C         E0_MEAN=1.5*TIWD+0.
        ELSE
          GOTO 998
        ENDIF
C
        WTOTA(IATM,ISTRA)=WTOTA(IATM,ISTRA)+WEIGHT
        ETOTA(ISTRA)=ETOTA(ISTRA)+E0*WEIGHT
        IF (NADSI.GE.1.AND.NLSRF(ISTRA)) CALL UPSUSR(WEIGHT,2)
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,1)
        IF (NLSTOR) CALL STORE(1)
C
C  PURELY MOLECULAR SOURCE?  300 --- 399
C
      ELSEIF (NLMOL(ISTRA)) THEN
C
        IF (NEMOD1.EQ.1) THEN
          EMAX=SORENI(ISTRA)
        ELSEIF (NEMOD1.EQ.6) THEN
          EMAX=0.
        ELSE
          GOTO 998
        ENDIF
C
        IF (EMAX.GT.0.D0) THEN
C  MONOENERGETIC DISTRIBUTION
          E0=EMAX
          VEL=RSQDVM(IMOL)*SQRT(E0)
C
C  COSINE LIKE OR GAUSSIAN ANGLE DISTRIBUTION
C
C  IN CASE (CRTX,CRTY,CRTZ) NE (0.,0.,0.)
C  ASSUME NORMAL INCIDENCE AND USE REFLECTION MODEL ANGULAR DISTRIBUTION
          VELX=CRTX
          VELY=CRTY
          VELZ=CRTZ
          CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
C         VEL_MEAN=VEL
C         E0_MEAN=E0
        ELSEIF (EMAX.LE.0..AND.TIWD.GT.0..AND..NOT.NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM SHIFTED TRUNCATED MAXWELLIAN FLUX
C              AROUND INNER (!) NORMAL AT TEMP. TIWL
C
          VWD=SQRT(VXWD**2+VYWD**2+VZWD**2)
          CALL VELOCS (TIWD,0.D0,VWD,VXWD,VYWD,VZWD,RSQDVM(IMOL),
     .                 CVRSSM(IMOL),
     .                 -CRTX,-CRTY,-CRTZ,E0,VELX,VELY,VELZ,VEL)
C  MODIFY ANGULAR DISTRIBUTION IN CASE SORCOS .NE. 0.5 (I.E., IN CASE
C  A NON-COSINE DISTRIBUTION IS REQUESTED
          IF (ABS(SORCOS(ISTRA)-0.5).GT.1.D-5) THEN
            VELX=CRTX
            VELY=CRTY
            VELZ=CRTZ
            CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                  SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
            VEL_MEAN=VEL
            E0_MEAN=E0
          ENDIF
        ELSEIF (EMAX.LE.0..AND.TIWD.GT.0..AND.NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM MAXWELLIAN AT TEMP. TW (EV) =TIWD
C
          NFLAG=2
          IDUM=1
          DUMT(1)=SQRT(TIWD/RMASSM(IMOL))*CVEL2A
          DUMT(2)=DUMT(1)
          DUMT(3)=DUMT(1)
          DUMV(1)=0
          DUMV(2)=0
          DUMV(3)=0
          CALL VELOCX(0,VXO,VYO,VZO,VO,IO,NO,VELQ,NFLAG,
     .                IDUM,DUMT,DUMV)
C         E0=VELQ*CVRSSM(IMOL)
C         E0_MEAN=1.5*TIWD+0.
        ELSE
          GOTO 998
        ENDIF
C
        WTOTM(IMOL,ISTRA)=WTOTM(IMOL,ISTRA)+WEIGHT
        ETOTM(ISTRA)=ETOTM(ISTRA)+WEIGHT*E0
        IF (NADSI.GE.1) CALL UPSUSR(WEIGHT,2)
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,1)
        IF (NLSTOR) CALL STORE(1)
C
C  PURELY TEST IONIC SOURCE?  400 --- 499
C
      ELSEIF (NLION(ISTRA)) THEN
C
        IF (NEMOD1.EQ.1) THEN
          EMAX=SORENI(ISTRA)
        ELSEIF (NEMOD1.EQ.2.OR.NEMOD1.EQ.3) THEN
          EMAX=SORENI(ISTRA)*TIWD+SORENE(ISTRA)*TEWD
        ELSEIF (NEMOD1.EQ.4.OR.NEMOD1.EQ.5) THEN
          VPERP=VXWD*CRTX+VYWD*CRTY+VZWD*CRTZ
          IF (VPERP.GT.0.D0) GOTO 996
          VPARX=VXWD-VPERP*CRTX
          VPARY=VYWD-VPERP*CRTY
          VPARZ=VZWD-VPERP*CRTZ
          VPAR=SQRT(VPARX**2+VPARY**2+VPARZ**2)
          VTERM=SQRT(TIWD/RMASSI(IION))*CVELAA
          VPERP=VPERP/VTERM
          VPAR=VPAR/VTERM
          EMAX=EMAXW(TIWD,VPERP,VPAR)
        ELSEIF (NEMOD1.EQ.6.OR.NEMOD1.EQ.7) THEN
          EMAX=0.
        ELSE
          GOTO 998
        ENDIF
C
        IF (EMAX.GT.0.D0) THEN
          E0=EMAX
          VEL=SQRT(E0)*RSQDVI(IION)
C
C  COSINE LIKE OR GAUSSIAN ANGLE DISTRIBUTION
C
C  IN CASE (CRTX,CRTY,CRTZ) NE (0.,0.,0.D0)
C  ASSUME NORMAL INCIDENCE AND USE REFLECTION MODEL ANGULAR DISTRIBUTION
          VELX=CRTX
          VELY=CRTY
          VELZ=CRTZ
          CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
          VEL_MEAN=VEL
          E0_MEAN=E0
        ELSEIF (EMAX.LE.0..AND.TIWD.GT.0..AND..NOT.NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM SHIFTED TRUNCATED MAXWELLIAN FLUX
C              AROUND INNER (!) NORMAL AT TEMP. TW (EV)
          VWD=SQRT(VXWD**2+VYWD**2+VZWD**2)
          CALL VELOCS (TIWD,0.D0,VWD,VXWD,VYWD,VZWD,RSQDVI(IION),
     .                  CVRSSI(IION),
     .                 -CRTX,-CRTY,-CRTZ,E0,VELX,VELY,VELZ,VEL)
C  MODIFY ANGULAR DISTRIBUTION IN CASE SORCOS .NE. 0.5 (I.E., IN CASE
C  A NON-COSINE DISTRIBUTION IS REQUESTED
          IF (ABS(SORCOS(ISTRA)-0.5).GT.EPS10) THEN
            VELX=CRTX
            VELY=CRTY
            VELZ=CRTZ
            CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                  SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
            VEL_MEAN=VEL
            E0_MEAN=E0
          ENDIF
        ELSEIF (EMAX.LE.0..AND.TIWD.GT.0..AND.NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM MAXWELLIAN AT TEMP. TW (EV) =TIWD
C
          NFLAG=2
          IDUM=1
          DUMT(1)=SQRT(TIWD/RMASSI(IION))*CVEL2A
          DUMT(2)=DUMT(1)
          DUMT(3)=DUMT(1)
          DUMV(1)=0
          DUMV(2)=0
          DUMV(3)=0
          CALL VELOCX(0,VXO,VYO,VZO,VO,IO,NO,VELQ,NFLAG,
     .                IDUM,DUMT,DUMV)
          E0=VELQ*CVRSSI(IION)
          E0_MEAN=1.5*TIWD+0.
        ELSE
          GOTO 998
        ENDIF
C
        WTOTI(IION,ISTRA)=WTOTI(IION,ISTRA)+WEIGHT
        ETOTI(ISTRA)=ETOTI(ISTRA)+E0*WEIGHT
        IF (NADSI.GE.1) CALL UPSUSR(WEIGHT,2)
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,1)
        IF (NLSTOR) CALL STORE(1)
C
C  PURELY BULK IONIC SOURCE?   500  ---  599
C
C  SOURCE DEFINED BY PRE COLLISION RATE OF BULK PARTICLES
C  THE RESULTING TEST PARTICLES MAY BE EITHER ATOMS, MOLECULES OR TEST
C  IONS. IN THIS CASE NOT THE TOTAL TEST PARTICLE FLUX BUT THE
C  THE TOTAL BULK ION FLUX IS SCALED TO A PRESCRIBED VALUE
C
C  SET ENERGY OF THE INCIDENT BULK ION : EMAX
C  IF EMAX.EQ.0., SAMPLE FROM SHIFTED TRUNCATED MAXWELLIAN
C  (ADD SHEATH CONTRIBUTION ESHET IF REQUIRED)
C
      ELSEIF (NLPLS(ISTRA)) THEN
        IF (output)
     .  WRITE(6,*) 'MARK: LOCATE: PROCESSING BULK SOURCE  ITYP= ',ITYP
C
        IF (NLSRF(ISTRA)) THEN
C
          IF (NEMOD1.EQ.1) THEN
            EMAX=SORENI(ISTRA)
          ELSEIF (NEMOD1.EQ.2.OR.NEMOD1.EQ.3) THEN
            EMAX=SORENI(ISTRA)*TIWD+SORENE(ISTRA)*TEWD
          ELSEIF (NEMOD1.EQ.4.OR.NEMOD1.EQ.5) THEN
            VPERP=VXWD*CRTX+VYWD*CRTY+VZWD*CRTZ
            IF (VPERP.LT.0.D0) GOTO 996
            VPARX=VXWD-VPERP*CRTX
            VPARY=VYWD-VPERP*CRTY
            VPARZ=VZWD-VPERP*CRTZ
            VPAR=SQRT(VPARX**2+VPARY**2+VPARZ**2)
            VTERM=SQRT(TIWD/RMASSP(IPLS))*CVELAA
            VPERP=VPERP/VTERM
            VPAR=VPAR/VTERM
            EMAX=EMAXW(TIWD,VPERP,VPAR)
          ELSEIF (NEMOD1.EQ.6.OR.NEMOD1.EQ.7) THEN
            EMAX=0.
          ELSEIF (NEMOD1.GT.7) THEN
            GOTO 998
          ENDIF
C
          IF (NEMOD1.EQ.3.OR.NEMOD1.EQ.5.OR.NEMOD1.EQ.7) THEN
            IF (FSHEAT(MSURF).LE.0.D0) THEN
              GAMMA=0.
              CUR=0.
              DO 550 IP=1,NPLSI
                VPWL(IP)=SQRT(VXWL(IP)**2+VYWL(IP)**2+VZWL(IP)**2)
C               DIWL(IP)=DIWL(IP)
550           CONTINUE
              ESHET=NCHRGP(IPLS)*SHEATH(TEWL,DIWL,VPWL,
     .                                  NCHRGP,GAMMA,CUR,NPLSI,MSURF)
            ELSE
              ESHET=NCHRGP(IPLS)*FSHEAT(MSURF)*TEWL
            ENDIF
          ELSE
            ESHET=0.
          ENDIF
C
          CRTX=-CRTX
          CRTY=-CRTY
          CRTZ=-CRTZ
C
          LOGPLS(IPLS,ISTRA)=.TRUE.
          IF (EMAX.GT.0.D0) THEN
C  CONSTANT VELOCITY NORMAL ONTO THE WALL
            E0=EMAX+ESHET
            VEL=SQRT(E0)*RSQDVP(IPLS)
            VELX=CRTX
            VELY=CRTY
            VELZ=CRTZ
            CALL REFANG(SORCOS(ISTRA),SORMAX(ISTRA),SORCTX(ISTRA),
     .                  SORCTY(ISTRA),SORCTZ(ISTRA),NAMODS(ISTRA),SNORM)
            E0_MEAN=E0
            VEL_MEAN=VEL
          ELSEIF (EMAX.LE.0.D0.AND.TIWD.GT.0.D0) THEN
C  SAMPLE FROM SHIFTED TRUNCATED MAXWELLIAN FLUX AND ACCELERATE IN SHEATH
            VWD=SQRT(VXWD**2+VYWD**2+VZWD**2)
            CALL VELOCS(TIWD,ESHET,VWD,VXWD,VYWD,VZWD,RSQDVP(IPLS),
     .                  CVRSSP(IPLS),
     .                  -CRTX,-CRTY,-CRTZ,E0,VELX,VELY,VELZ,VEL)
          ENDIF
C
          CRTX=-CRTX
          CRTY=-CRTY
          CRTZ=-CRTZ
          IF (output)
     .    WRITE(6,*) 'MARK: LOCATE: ION CREATED  ITYP= ',ITYP
C
C  A BULK ION, HITTING A SURFACE, HAS BEEN CREATED.
C
C  UPDATE PARTICLE EFFLUX  ONTO SURFACE MSURF
C  UPDATE ENERGY FLUX ONTO SURFACE MSURF
C
C  SPATIAL RESOLUTION ON NON DEFAULT STANDARD SURFACE?
c slmod begin - juelich - not tr (fixed in new version)
          IF (MSURF.GT.NLIM.AND.LEVGEO.LT.4.AND.NLMPGS.GT.NLIMPS) THEN
c
c          IF (MSURF.GT.NLIM.AND.LEVGEO.LT.4) THEN
c slmod end
            ISTS=MSURF-NLIM
            IF (INUMP(ISTS,1).NE.0) MSURFG=NPCELL+(NTCELL-1)*NP2T3
            IF (INUMP(ISTS,2).NE.0) MSURFG=NRCELL+(NTCELL-1)*NR1P2
            IF (INUMP(ISTS,3).NE.0) MSURFG=NRCELL+(NPCELL-1)*NR1P2
            MSURFG=NLIM+NSTSI+MSURFG+(ISTS-1)*NGITT
            FLX=FLXOUT(MSURFG)
c slmod begin - debug - not tr
c            IF (DEBUGOPT.NE.0) THEN
c              WRITE(0,*) 'ESCAPE: ASSIGNING MSURFG 01 =',msurfg
c              WRITE(6,*) 'ESCAPE: ASSIGNING MSURFG 01 =',msurfg
c            ENDIF
c slmod end
          ELSE
            MSURFG=0
            FLX=FLXOUT(MSURF)
          ENDIF
C
          WTOTP(IPLS,ISTRA)=WTOTP(IPLS,ISTRA)-WEIGHT
          ETOTP(ISTRA)=ETOTP(ISTRA)-E0*WEIGHT
          POTPL(IPLS,MSURF)=POTPL(IPLS,MSURF)+WEIGHT
          EOTPL(IPLS,MSURF)=EOTPL(IPLS,MSURF)+WEIGHT*E0
          IF (MSURFG.GT.0) THEN
            POTPL(IPLS,MSURFG)=POTPL(IPLS,MSURFG)+WEIGHT
            EOTPL(IPLS,MSURFG)=EOTPL(IPLS,MSURFG)+E0*WEIGHT
          ENDIF
          IF (NADSI.GE.1) CALL UPSUSR(-WEIGHT,1)
C
          IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,1)
          IF (NLSTOR) CALL STORE(2)
C
C  REFLECT THIS ION AS TEST PARTICLE FROM SURFACE NO. MSURF
C
C  BUT FIRST: CALL SPUTTER MODEL IF REQUIRED
C
          FMASS=DBLE(NMASSP(IPLS))
          FCHAR=DBLE(NCHARP(IPLS))
C
          WGHTSP=0.
          WGHTSC=0.
          YIELD1=0.
          YIELD2=0.
C
          NLSPUT=.FALSE.
          IF (ILSPT(MSURF).NE.0) THEN
C  SAVE INCIDENT PARTICLE'S SPEED AND ENERGY
            E0S=E0
            WEIGHS=WEIGHT
            VELS=VEL
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            ISPZS=ISPZ
C
            CALL SPUTR1(WMINS,FMASS,FCHAR,FLX,
     .                  ISRS(ISPZ,MSURF),
     .                  YIELD1,
     .                  ISSPTP,ESPTP,VSPTP,VXSPTP,VYSPTP,VZSPTP,
     .                  ISRC(ISPZ,MSURF),
     .                  YIELD2,
     .                  ISSPTC,ESPTC,VSPTC,VXSPTC,VYSPTC,VZSPTC)
            NLSPUT=YIELD1.GT.0..OR.YIELD2.GT.0.
            WGHTSP=WEIGHT*YIELD1
            WGHTSC=WEIGHT*YIELD2
C
C  UPDATE SPUTTER SURFACE TALLIES
C
            SPTPL(IPLS,MSURF)=SPTPL(IPLS,MSURF)+WGHTSP+WGHTSS
          ENDIF
C
C  PHYSICAL SPUTTERING
C
          IF (YIELD1.GT.0..AND.ISRS(ISPZ,MSURF).NE.0) THEN
C  FOLLOW SPUTTERED PARTICLES LATER. PUT THEM INTO STATISTICAL CELLAR
C
            ISPZ=ISSPTP
            ITYP=ISPEZI(ISPZ,0)
            IATM=ISPEZI(ISPZ,1)
            IMOL=ISPEZI(ISPZ,2)
            IION=ISPEZI(ISPZ,3)
            IPLS=ISPEZI(ISPZ,4)
            E0=ESPTP
            WEIGHT=WGHTSP
            VEL=VSPTP
            VELX=VXSPTP
            VELY=VYSPTP
            VELZ=VZSPTP
C
C.....................................................................
C  SPLITTING
C
            NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
            DO 533 J=1,NPARTC
              RSPLST(NLEVEL,J)=RPST(J)
533         CONTINUE
            DO 534 J=1,MPARTC
              ISPLST(NLEVEL,J)=IPST(J)
534         CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
            NODES(NLEVEL)=2
C
C  SPLITTING DONE. NEXT: SURFACE TALLIES
C.....................................................................
C
            IF (NLTRC.AND.TRCHST) THEN
              WRITE (6,*) 'AFTER SUBR. SPUTER: PHYS. SPUTTERING'
              WRITE (6,'(1X,A8)') TEXTS(ISPZ)
              CALL MASR1('YIELDP  ',YIELD1)
              CALL MASR6 (
     .           'VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .            VELX,VELY,VELZ,VEL,E0,WEIGHT)
            ENDIF
C
            IF (ITYP.EQ.1) THEN
              PPATI(IATM,ISTRA)=PPATI(IATM,ISTRA)+WEIGHT
              EPATI(ISTRA)=EPATI(ISTRA)+E0*WEIGHT
              PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WEIGHT
              ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WEIGHT
                ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITYP.EQ.2) THEN
              PPMLI(IMOL,ISTRA)=PPMLI(IMOL,ISTRA)+WEIGHT
              EPMLI(ISTRA)=EPMLI(ISTRA)+E0*WEIGHT
              PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WEIGHT
              ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WEIGHT
                ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITYP.EQ.3) THEN
              PPIOI(IION,ISTRA)=PPIOI(IION,ISTRA)+WEIGHT
              EPIOI(ISTRA)=EPIOI(ISTRA)+E0*WEIGHT
              PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WEIGHT
              ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WEIGHT
                ERFPIO(IION,MSURFG)=ERFPIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
            IF (NADSI.GE.1) CALL UPSUSR(WEIGHT,2)
          ENDIF
C
C  CHEMICAL SPUTTERING
C
          IF (YIELD2.GT.0..AND.ISRC(ISPZ,MSURF).NE.0.) THEN
C  FOLLOW SPUTTERED PARTICLES LATER. PUT THEM INTO STATISTICAL CELLAR
            ISPZ=ISSPTC
            ITYP=ISPEZI(ISPZ,0)
            IATM=ISPEZI(ISPZ,1)
            IMOL=ISPEZI(ISPZ,2)
            IION=ISPEZI(ISPZ,3)
            IPLS=ISPEZI(ISPZ,4)
            E0=ESPTC
            WEIGHT=WGHTSC
            VEL=VSPTC
            VELX=VXSPTC
            VELY=VYSPTC
            VELZ=VZSPTC
C
C.....................................................................
C  SPLITTING
C
            NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
            DO 535 J=1,NPARTC
              RSPLST(NLEVEL,J)=RPST(J)
535         CONTINUE
            DO 536 J=1,MPARTC
              ISPLST(NLEVEL,J)=IPST(J)
536         CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
            NODES(NLEVEL)=2
C
C  SPLITTING DONE. NEXT: SURFACE TALLIES
C.....................................................................
C
            IF (NLTRC.AND.TRCHST) THEN
              WRITE (6,*) 'AFTER SUBR. SPUTER: CHEM. SPUTTERING'
              WRITE (6,'(1X,A8)') TEXTS(ISPZ)
              CALL MASR1('YIELDC  ',YIELD2)
              CALL MASR6 (
     .           'VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .            VELX,VELY,VELZ,VEL,E0,WEIGHT)
            ENDIF
C
            IF (ITYP.EQ.1) THEN
              PPATI(IATM,ISTRA)=PPATI(IATM,ISTRA)+WEIGHT
              EPATI(ISTRA)=EPATI(ISTRA)+E0*WEIGHT
              PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WEIGHT
              ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WEIGHT
                ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITYP.EQ.2) THEN
              PPMLI(IMOL,ISTRA)=PPMLI(IMOL,ISTRA)+WEIGHT
              EPMLI(ISTRA)=EPMLI(ISTRA)+E0*WEIGHT
              PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WEIGHT
              ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WEIGHT
                ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITYP.EQ.3) THEN
              PPIOI(IION,ISTRA)=PPIOI(IION,ISTRA)+WEIGHT
              EPIOI(ISTRA)=EPIOI(ISTRA)+E0*WEIGHT
              PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WEIGHT
              ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WEIGHT
                ERFPIO(IION,MSURFG)=ERFPIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
            IF (NADSI.GE.1) CALL UPSUSR(WEIGHT,2)
          ENDIF
C
C  RESTORE INCIDENT PARTICLE, FOR SURFACE REFLECTION ROUTINE
C
          IF (ILSPT(MSURF).NE.0) THEN
            E0=E0S
            WEIGHT=WEIGHS
            VEL=VELS
            VELX=VELXS
            VELY=VELYS
            VELZ=VELZS
            ISPZ=ISPZS
            LGPART=.FALSE.
          ENDIF
C
C
C  NEXT: CALL REFLECTION MODEL
C
540       CONTINUE
          IF (output)
     .    WRITE(6,*) 'MARK: LOCATE: REFLECTION MODEL  ITYP= ',ITYP
c          WRITE(0,*) 'debug reflc1:',ityp
          CALL REFLC1 (WMINS,FMASS,FCHAR,NPRT(ISPZ),
     .                 ISRF(ISPZ,MSURF),ISRT(ISPZ,MSURF))
          ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
C
          IF (NLTRC.AND.TRCHST) THEN
            IF (LGPART) THEN
              WRITE (6,*) 'AFTER SUBR. REFLEC SL: '
              WRITE (6,'(1X,A8)') TEXTS(ISPZ)
              CALL MASR6 (
     .           'VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .            VELX,VELY,VELZ,VEL,E0,WEIGHT)
            ELSE
              WRITE (6,*) 'ABSORBED IN SUBR. REFLEC'
            ENDIF
          ENDIF
C
          IF (ITYP.EQ.1) THEN
            PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WEIGHT
            ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+E0*WEIGHT
            IF (MSURFG.GT.0) THEN
              PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WEIGHT
              ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+E0*WEIGHT
            ENDIF
          ELSEIF (ITYP.EQ.2) THEN
            PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WEIGHT
            ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+E0*WEIGHT
            IF (MSURFG.GT.0) THEN
              PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WEIGHT
              ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+E0*WEIGHT
            ENDIF
          ELSEIF (ITYP.EQ.3) THEN
            PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WEIGHT
            ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+E0*WEIGHT
            IF (MSURFG.GT.0) THEN
              PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WEIGHT
              ERFPIO(IION,MSURFG)=ERFPIO(IION,MSURFG)+E0*WEIGHT
            ENDIF
          ENDIF
C
          IF (NLSTOR) CALL STORE(1)
          IF (NADSI.GE.1) CALL UPSUSR(WEIGHT,2)
C
        ELSEIF (NLVOL(ISTRA)) THEN
C
C  SAMPLE FROM MAXWELLIAN AT LOCAL PLASMA PARAMETERS TIIN AND (VXIN,VYIN,VZIN)
C  IN CELL ICELL=NCELL
C
          NFLAG=2
          IDUM=1
          CALL VELOCX(NCELL,VXO,VYO,VZO,VO,IO,NO,VELQ,NFLAG,
     .                IDUM,DUMT,DUMV)
          E0=VELQ*CVRSSP(IPLS)
          WTOTP(IPLS,ISTRA)=WTOTP(IPLS,ISTRA)-WEIGHT
          ETOTP(ISTRA)=ETOTP(ISTRA)-E0*WEIGHT
          IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,1)
          IF (NLSTOR) CALL STORE(2)
C  RECOMBINING BULK ION (IPLS,E0,WEIGHT,...) IS NOW IDENTIFIED
C  FIND TYPE AND SPECIES OF NEW TEST PARTICLE FROM RECOMB. PROCESS: IREC
          IREC=SORIND(IVOLM,ISTRA)
          IF (NATPRC(IREC).GT.0) THEN
            ITYP=1
            IATM=NATPRC(IREC)
            IF (IATM.GT.NATMI) GOTO 999
          ELSEIF (NMLPRC(IREC).GT.0) THEN
            ITYP=2
            IMOL=NMLPRC(IREC)
            IF (IMOL.GT.NMOLI) GOTO 999
          ELSEIF (NIOPRC(IREC).GT.0) THEN
            ITYP=3
            IION=NIOPRC(IREC)
            IF (IION.GT.NIONI) GOTO 999
          ELSE
            GOTO 999
          ENDIF
          IF (NLTRC.AND.TRCHST) THEN
            WRITE (6,*) 'AFTER RECOMBINATION: '
            CALL MASJ5 ('ITYP,IATM,IMOL,IION,IPLS                ',
     .                   ITYP,IATM,IMOL,IION,IPLS)
          ENDIF
C
          IF (NLSTOR) CALL STORE(1)
C
        ELSEIF (NLLNE(ISTRA)) THEN
          WRITE (6,*) 'BULK ION LINE SOURCE NOT READY, EXIT CALLED '
          CALL EXIT
C
        ELSEIF (NLPNT(ISTRA)) THEN
          WRITE (6,*) 'BULK ION POINT SOURCE NOT READY, EXIT CALLED '
          CALL EXIT
        ENDIF
C
        IF (ITYP.EQ.1) THEN
          PPATI(IATM,ISTRA)=PPATI(IATM,ISTRA)+WEIGHT
          EPATI(ISTRA)=EPATI(ISTRA)+E0*WEIGHT
        ELSEIF (ITYP.EQ.2) THEN
          PPMLI(IMOL,ISTRA)=PPMLI(IMOL,ISTRA)+WEIGHT
          EPMLI(ISTRA)=EPMLI(ISTRA)+E0*WEIGHT
        ELSEIF (ITYP.EQ.3) THEN
          PPIOI(IION,ISTRA)=PPIOI(IION,ISTRA)+WEIGHT
          EPIOI(ISTRA)=EPIOI(ISTRA)+E0*WEIGHT
        ENDIF
C
        IF (output)
     .  WRITE(6,*) 'MARK: LOCATE: 5000 ITYP= ',ITYP
      ENDIF
C
5000  CONTINUE
C
C  HAS THE SOURCE PARTICLE BEEN ABSORBED IN SUBR. REFLEC OR SPUTER?
C
      IF (.NOT.LGPART) RETURN
C
C  IS THE PARTICLE LAUNCHED OUTSIDE THE COMPUTATIONAL BOX?
C
C  TEST FOR CORRECT CELL NUMBER AT BIRTH POINT
C  KILL PARTICLE, IF WRONG CELL INDICES
C
      IF (NLTEST) THEN
        CALL CLLTST(*997)
      ELSE
        NLTST=.FALSE.
        NLTST=NLTST.OR.(NLRAD.AND.(NRCELL.GT.NR1ST.OR.NRCELL.LT.0))
        NLTST=NLTST.OR.(NLPOL.AND.(NPCELL.GT.NP2ND.OR.NPCELL.LT.1))
        NLTST=NLTST.OR.(NLTOR.AND.(NTCELL.GT.NT3RD.OR.NTCELL.LT.1))
        NLTST=NLTST.OR.(NRCELL.EQ.0.AND.
     .                            (NACELL.GT.NRADD.OR.NACELL.LT.1))
        IF (NLTST) GOTO 995
      ENDIF
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN LOCATE: ILSIDE OF SOURCE SURFACE IS 0.'
      WRITE (6,*) 'THUS NO OUTER NORMAL CAN BE DEFINED. EXIT CALLED'
      WRITE (6,*) 'SET EITHER ILSIDE NE 0 OR USE SORIFL FLAG '
      WRITE (6,*) 'MSURF,ISTSF,NRCELL,NPCELL,NTCELL '
      WRITE (6,*)  MSURF,ITRSF,NRCELL,NPCELL,NTCELL
      CALL EXIT
991   CONTINUE
      WRITE (6,*) 'ERROR IN LOCATE: INCONSISTENT INPUT FLAGS   '
      WRITE (6,*) 'MSURF = ',MSURF
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN LOCATE: RADON-NIKODYM CONDITION    '
      WRITE (6,*) 'VIOLATED FOR NON-ANALOG SOURCE SPECIES SAMPLING'
      WRITE (6,*) 'CHECK DATM,DMOL,DION OR DPLS ARRAYS (BLOCK) 6 '
      CALL EXIT
995   CONTINUE
      WRITE (6,*) 'PARTICLE LAUNCHED OUTSIDE THE COMPUTATIONAL BOX'
      WRITE (6,*) 'OR WITH INVALID CELL INDICES'
      WRITE (6,*) 'NPANU,X0,Y0,Z0 ',NPANU,X0,Y0,Z0
      WRITE (6,*) 'NRCELL,NPCELL,NTCELL,NBLOCK,NACELL ',
     .             NRCELL,NPCELL,NTCELL,NBLOCK,NACELL
c slmod begin - tr
      PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      NLOST=NLOST+1     
      RETURN
c
c      CALL EXIT
c slmod end
996   CONTINUE
      WRITE (6,*) 'BULK ION LAUNCHED IN WRONG DIRECTION'
      WRITE (6,*) 'NPANU,VXWL,VYWL,VZWL ',NPANU,VXWL(IPLS),VYWL(IPLS),
     .                                    VZWL(IPLS)
      WRITE (6,*) '      CRTX,CRTY,CRTZ ',CRTX,CRTY,CRTZ
      CALL EXIT
997   CONTINUE
      WRITE (6,*) 'TEST PARTICLE LAUNCHED WITH INVALID CELL INDICES'
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,15)
      IF (NLSTOR) CALL STORE(100)
      WEIGHT=0.
      LGPART=.FALSE.
      RETURN
998   CONTINUE
      WRITE (6,*) 'ERROR IN LOCATE: NEMODS,ITYP= ',NEMODS(ISTRA),ITYP
      WRITE (6,*) 'INVALID OPTION. TIWD= ',TIWD
      CALL EXIT
999   CONTINUE
      WRITE (6,*) 'ERROR IN LOCATE: TYP OR SPECIES OUT OF RANGE'
      CALL EXIT
      END
C
