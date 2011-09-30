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
