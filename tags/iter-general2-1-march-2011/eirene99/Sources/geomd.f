*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
c slmod begin - not tr
      SUBROUTINE GEOMD(NDXA,NDYA,XPLG,YPLG,NPLP,NPNT,NR1ST,
     .                 PUX,PUY,PVX,PVY,NBMLT)
c
c      SUBROUTINE GEOMD(NDXA,NDYA,XPLG,YPLG,NPLP,NPNT,NR1ST,
c     .                 PUX,PUY,PVX,PVY,NLMLT)
c slmod end
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'PARMMOD'
C
C
      DIMENSION XPLG(N1STS,*),YPLG(N1STS,*),NPNT(2,*),
     .          PUX(*),PUY(*),PVX(*),PVY(*)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      COMMON /LINEIR/
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
C
      CHARACTER*80 LINE,LINE2
C   DIMENSIONIERUNG FUER GITTER
      INTEGER DIMXH,DIMYH,NNCUT,
     1 NXCUT1(10),NXCUT2(10),NYCUT1(10),NYCUT2(10)

      DIMENSION DUMMI(3)
      REAL*8 MERK(NDY)
c slmod begin - grid - not tr (moved to INFCOP)
      IF (GEOMOPT.EQ.1) THEN
        CALL GEODIV (NDXA,NDYA,XPLG,YPLG,NPLP,NPNT,NR1ST,
     .               PUX,PUY,PVX,PVY,NBMLT)
        RETURN
      ENDIF
c slmod end
C  ACTUAL MESH USED IN THIS RUN
C
C      EINLESEROUTINE ANGEPASST AUF BRAAMS-OUTPUT
C   GEAENDERTE DIMENSIONIERUNG BZW. CUT-POSITION
C        MUSS PER HAND ANGEPASST WERDEN:
C        PARAMETER DIMXH,DIMYH                        RFS 14.5.1991
       REWIND 30
       READ (30,*)
       READ (30,*)
3366  FORMAT(/)
      read(30,*) dimxh,dimyh,nncut
      if(nncut.gt.10) stop 'Increase array sizes for cut'
      read(30,*) (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut)
      read(30,*)
      DO 10 IX = 1, DIMXH
       IF (IX.LE.nxcut1(1)-1) THEN
        DO 12 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,3333) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),XPLG(IY,IX)
          READ (30,3333) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),YPLG(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),XPLG(dimyh+1,IX),XPLG(IY,IX)
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),YPLG(dimyh+1,IX),YPLG(IY,IX)
         ENDIF
12      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(1)) THEN
        DO 14 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,3333) XPLG(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),XPLG(IY,IX)
         READ (30,3333) YPLG(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),YPLG(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) XPLG(IY,nxcut2(2)),XPLG(dimyh+1,nxcut2(2)),
     .                                XPLG(dimyh+1,IX),XPLG(IY,IX)
          READ (30,3333) YPLG(IY,nxcut2(2)),YPLG(dimyh+1,nxcut2(2)),
     .                                YPLG(dimyh+1,IX),YPLG(IY,IX)
         ENDIF
14      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(2)).AND.(IX.LE.nxcut1(2)-1)) THEN
        DO 16 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),XPLG(IY,IX+1)
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),YPLG(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   XPLG(dimyh+1,IX+1),XPLG(IY,IX+1)
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   YPLG(dimyh+1,IX+1),YPLG(IY,IX+1)
         ENDIF
16      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(2)) THEN
        DO 18 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,3333) XPLG(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),XPLG(IY,IX+1)
         READ (30,3333) YPLG(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),YPLG(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) XPLG(IY,nxcut2(1)+1),XPLG(dimyh+1,nxcut2(1)+1),
     .                                XPLG(dimyh+1,IX+1),XPLG(IY,IX+1)
          READ (30,3333) YPLG(IY,nxcut2(1)+1),YPLG(dimyh+1,nxcut2(1)+1),
     .                                YPLG(dimyh+1,IX+1),YPLG(IY,IX+1)
         ENDIF
18      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(1)).AND.(IX.LE.dimxh-1)) THEN
        DO 22 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,3333) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   XPLG(IY,IX+2)
          READ (30,3333) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   YPLG(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   XPLG(dimyh+1,IX+2),XPLG(IY,IX+2)
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   YPLG(dimyh+1,IX+2),YPLG(IY,IX+2)
         ENDIF
22      CONTINUE
       ENDIF
       IF (IX.EQ.dimxh) THEN
        DO 24 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,3333) XPLG(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  XPLG(IY,IX+2)
         READ (30,3333) YPLG(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  YPLG(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) XPLG(IY,dimxh+3),XPLG(dimyh+1,dimxh+3),
     .                                XPLG(dimyh+1,IX+2),XPLG(IY,IX+2)
          READ (30,3333) YPLG(IY,dimxh+3),YPLG(dimyh+1,dimxh+3),
     .                                YPLG(dimyh+1,IX+2),YPLG(IY,IX+2)
         ENDIF
24      CONTINUE
       ENDIF
10    CONTINUE

3333  FORMAT(4E15.7)

C    ANZAHL DER TEILSTUECKE PRO POLYGON
       NPLP = 3
C   ANFANGSPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
       NPNT(1,1)=1
C   ENDPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
       NPNT(2,1)=nxcut2(2)
C   ANFANGSPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
       NPNT(1,2)=nxcut2(2)+1
C   ENDPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
       NPNT(2,2)=nxcut2(1)+1
C   ANFANGSPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
       NPNT(1,3)=nxcut2(1)+2
C   ENDPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
       NPNT(2,3)=dimxh+3
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=XPLG(IY,IX)
          Y1(IX)=YPLG(IY,IX)
          X2(IX)=XPLG(IY,IX+1)
          Y2(IX)=YPLG(IY,IX+1)
          X3(IX)=XPLG(IY+1,IX)
          Y3(IX)=YPLG(IY+1,IX)
          X4(IX)=XPLG(IY+1,IX+1)
          Y4(IX)=YPLG(IY+1,IX+1)
1014    CONTINUE
C
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
C
      NP=NPNT(2,NPLP)
      DO 1020 J=1,NDYA+1
        DO 1020 I=1,NP
          XPLG(J,I)=XPLG(J,I)*100.
          YPLG(J,I)=YPLG(J,I)*100.
          IF (ABS(XPLG(J,I)).LT.5.D-5) XPLG(J,I)=0.
          IF (ABS(YPLG(J,I)).LT.5.D-5) YPLG(J,I)=0.
1020  CONTINUE
c slmod begin - not tr
      WRITE(61,*) 'GEOMETRY DATA - GEOMD'
      WRITE(61,*)
      WRITE(61,*) 'ndxa   ',ndxa
      WRITE(61,*) 'ndya   ',ndya
      WRITE(61,*) 'nplp   ',nplp
      WRITE(61,*)
      WRITE(61,*) '1: ',npnt(1,1),npnt(2,1)
      WRITE(61,*) '2: ',npnt(1,2),npnt(2,2)
      WRITE(61,*) '3: ',npnt(1,3),npnt(2,3)

      DO ix = 1, npnt(2,nplp)
         WRITE(61,*)
         DO iy = 1, ndya+1
          WRITE(61,'(2F10.4,2I6)') xplg(iy,ix),yplg(iy,ix),iy,ix
        ENDDO
      ENDDO
c slmod end
      RETURN
*//END GEOMD//
      END
