C
C
      SUBROUTINE
     .  EIRENE_ZYLND2(X0,Y0,Z0,VX,VY,VZ,TAR,RAD,NK,NP,NA,IO,NF,NUM,
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
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
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
        CALL EIRENE_EXIT_OWN(1)
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
          CALL
     .  EIRENE_SHNITT(P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,1,NA+1,1)
          ELSEIF (I.EQ.NK.AND.IRIGHT.NE.0) THEN
          WRITE (iunout,*) ' RIGHT END OF ZYLINDER '
          WRITE (iunout,*) (AR(IRGHT),IRGHT=1,IRIGHT)
          CALL
     .  EIRENE_SHNITT(P,PXS,PYS,PZS,VX,VY,VZ,AR,IRIGHT,XP,YP,1,NA+1,1)
          ELSE
            DO 3 J=1,NA+1
              PXX=P(1,J)+PX
              PYY=P(2,J)+PY
              PZZ=P(3,J)+PZ
3             CALL EIRENE_PL3D (PXX,PYY,PZZ,XP(J),YP(J))
          ENDIF
          IF (IO.GE.2) CALL GRNWPN(IO)
          do 7 jj=1,na+1
            xps(jj)=xp(jj)
            yps(jj)=yp(jj)
7         continue
          CALL GRLN (XPS,YPS,NA+1)
C  FAERBE DIE ENDEN DES ZYLINDERS EIN
          IF ((I.EQ.1.OR.I.EQ.NK).AND.NF) CALL
     .  GRFILL(NA+1,XPS,YPS,1,1)
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
            CALL
     .  EIRENE_SHNITT(P,PXS,PYS,PZS,-VX,-VY,-VZ,AL,ILEFT,XP,YP,J,J,JP)
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
            CALL EIRENE_PL3D
     .  (P(1,J)+PX,P(2,J)+PY,P(3,J)+PZ,XP(JP),YP(JP))
            GOTO 4
12        CONTINUE
          GOTO 14
4       CONTINUE
        IF (IRIGHT.NE.0) THEN
          DO 13 IP=1,IPART(NK)
            IF (PHIDEG.LT.PHIAN(NK,IP).OR.PHIDEG.GT.PHIEN(NK,IP))
     .          GOTO 13
            JP=JP+1
            CALL EIRENE_SHNITT
     .  (P,PXS,PYS,PZS,VX,VY,VZ,AR,IRIGHT,XP,YP,J,J,JP)
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
