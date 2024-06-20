C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRDRLG(drdmpa,TEXTX,TEXTY,TEXTZ,IOPT,XX,YY,ZZ,IX,IY,IZ)
C     LETZTES UPDATE 29.9.89 GROTEN (GRAXLOG, LETZTER PARAMETER INTEGER)
C     LETZTES UPDATE 25.9.90 GROTEN
C     LETZTES UPDATE 02.5.96 GROTEN
C     BEREICHE FUER IOPT=1,2:
C     IOPT=+1,+2: ZUSAETZLICH KOORDINATENACHSEN!
C                 DANN WERDEN IX,IY,IZ BENOETIGT:
C                 >0: LINEARE ACHSE
C                 <0: LOGARITHMISCHE ACHSE
C                 =0: KEINE ACHSE
C***********************************************************************
      CHARACTER(len=*) TEXTX,TEXTY,TEXTZ
      DIMENSION drdmpa(11),XR(4),YR(4)
      LOGICAL OKX,OKY,OKZ,LINKS
      DIMENSION Z(2,4),ZM(2,3),ZL(2),ZO(2),ZR(2),ZU(2)
      REAL LMAX,LTEXT,LPFEIL,lmas,XA1,YA1,XA2,YA2
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRPP/
      EQUIVALENCE (PP(16),INTSYM)
      EQUIVALENCE (ZM(1,1),ZL(1)), (ZM(1,2),ZO(1)), (ZM(1,3),ZR(1)),
     *            (ZO(1),ZU(1))
C          EPSWIN:   ZUR GUELTIGKEITSABFRAGE FUER DIE ACHSENBESCHRIFTUNG
C          GRDFAK:   UMRECHNUNG VON BOGENMASS IN ALTGRAD
C          BOGFAK:   UMRECHNUNG VON ALTGRAD IN BOGENMASS
C          SPACE:    ABSTAND DES TEXTES ZUR LINIE BZGL. DER TEXTGROESSE
C          SPACE2:   SPACE FUER IOPT=1
      DATA EPSWIN/1./, GRDFAK/57.2958/, BOGFAK/0.0174534/,
     *     SPACE/0.5/, SPACE2/0.4/

 
C     query number of Workstations
C     falls der Dummy Treiber mit WK Typ 4711 existiert ist man im
C     CMS und IWK=2 wird gesetzt, da query text extent auf die dummy
C     wk gemacht wird
C     sonst wird angenommen, das GKS auf der CRAY oder einem anderen
C     GKS laeuft und nur eine WK mit WKID 1 offen ist
C
      CALL GQEWK(1,IERR,NUMBER,IWKT)
      DO 40 I=1,NUMBER
      CALL GQEWK(I,IERR,NUMBER,IWKT)
c     testweise , da Aerger mit Versatec Routinen
      IF (IWKT.GT.16000) GOTO 40
      IF ( IWKT.EQ.4711 .and. ierr.eq.0 ) THEN
         IWKID=2
         goto  50
      ENDIF
 40   CONTINUE
      IWKID=1
 50   CONTINUE

c 14.11.91
C      query wk connection and type
      CALL GQWKC(IWKID,IERR,ICon,IWT)
      CALL GQWKCA (IWT,IERR, ICAT)

C
C     drdmpa(2)=PHI     |PHI|<=360 GRAD
C     drdmpa(3)=THETA   |THETA|<=90 GRAD
C     drdmpa(4)=BREITE
C     drdmpa(5)=RET.CODE AUS GR-PROG.
C
      IF (drdmpa(5) .GT. 0.) then
         write(*,*) 'GRDRLG: drdmpa(5) > 0  at start'
         return
      end if
      IF (ABS(IOPT).LT.1 .OR. ABS(IOPT).GT.2) THEN
         drdmpa(5)=100.
         write(*,*) 'GRDRLG: drdmpa(5) = 100 (IOPT not ok)'
         RETURN
      ELSE
         drdmpa(5)=0.
      END IF
************************************************************************
*     WELCHE ACHSEN DUERFEN BESCHRIFTET WERDEN?                        *
*     KEINE X-ACHSENBESCHRIFTUNG, FALLS                                *
*           THETA=0 & |PHI|=90 ODER 270                                *
*     KEINE Y-ACHSENBESCHRIFTUNG, FALLS                                *
*           THETA=0 & |PHI|=180 ODER 360                               *
*     KEINE Z-ACHSENBESCHRIFTUNG, FALLS                                *
*           |THETA|=90                                                 *
************************************************************************
      PHI=drdmpa(2)
      THETA=drdmpa(3)
      ABSP=ABS(PHI)
      ABST=ABS(THETA)
C
      OKX=ABST.GT.EPSWIN .AND. ABS(ABSP-90.).GT.EPSWIN
     *                   .AND. ABS(ABSP-270.).GT.EPSWIN
      OKY=ABST.GT.EPSWIN .AND. ABSP.GT.EPSWIN
     *                   .AND. ABS(ABSP-180.).GT.EPSWIN
     *                   .AND. ABS(ABSP-360.).GT.EPSWIN
      OKZ=ABS(ABST-90.) .GT. EPSWIN
C
      SINP=SIN(PHI*BOGFAK)
      COSP=COS(PHI*BOGFAK)
      SINT=SIN(THETA*BOGFAK)
      COST=COS(THETA*BOGFAK)
      PP14=PP(14)
      PP15=PP(15)
      BREITE=drdmpa(4)
      X0=BREITE*0.5
      Y0=X0
c 14.11.91
c     nicht gqtxx , wenn NICHT GTSGRAL-GKS und MO- wktyp
      IF (iwkid.eq.1.and. ICAT.EQ.4)  goto 9999
C
      IF (OKX) THEN
        DO 1 L=LEN(TEXTX),1,-1
          LX=L
          IF (TEXTX(L:L) .NE. ' ') GOTO 2
    1   continue
        LX=0
    2   PHIX=ATAN2(-SINP*SINT,COSP)
        PHIXG=PHIX*GRDFAK
        CALL GRCHRC(PP14,PHIXG,INTSYM)
        XA1=XC(X0,PHIX,4.)
        YA1=YC(Y0,PHIX,4.)
        IF (LX.GT.0) THEN
          CALL GQTXX(IWKID,XA1,YA1,TEXTX(:LX),IRC,XE1,YE1,XR,YR)
        ELSE
          XE1=XA1
          YE1=YA1
        END IF
C       STEHT DER TEXT AUF DEM KOPF?
        IF (COS(PHIX) .GE. 0.) THEN
          PHIXT=PHIXG
          XA1=XA1-X0
          YA1=YA1-Y0
        ELSE
          PHIXT=PHIXG+180.
          XA1=XE1-X0
          YA1=YE1-Y0
        END IF

        PHIX90=(PHIXT+90.)*BOGFAK
        XA1=XC(XA1,PHIX90,SPACE)
        YA1=YC(YA1,PHIX90,SPACE)
        XR(3)=XC(XR(2),PHIX90,1.+SPACE+SPACE)
        YR(3)=YC(YR(2),PHIX90,1.+SPACE+SPACE)
        XR(4)=XC(XR(1),PHIX90,1.+SPACE+SPACE)
        YR(4)=YC(YR(1),PHIX90,1.+SPACE+SPACE)
      ELSE
      END IF
      IF (OKY) THEN
        DO 3 L=LEN(TEXTY),1,-1
          LY=L
          IF (TEXTY(L:L) .NE. ' ') GOTO 4
    3   continue
        LY=0
    4   PHIY=ATAN2(COSP*SINT,SINP)
        PHIYG=PHIY*GRDFAK
        CALL GRCHRC(PP14,PHIYG,INTSYM)
        XA2=XC(X0,PHIY,4.)
        YA2=YC(Y0,PHIY,4.)
        IF (LY.GT.0) THEN
          CALL GQTXX(IWKID,XA2,YA2,TEXTY(:LY),IRC,XE2,YE2,XR,YR)
        ELSE
          XE2=XA2
          YE2=YA2
        END IF
        XL=XC(XE2,PHIY,4.)
        YL=YC(YE2,PHIY,4.)

C       STEHT DER TEXT AUF DEN KOPF?
        IF (COS(PHIY) .GE. 0.) THEN
          PHIYT=PHIYG
          XA2=XA2-X0
          YA2=YA2-Y0
        ELSE
          PHIYT=PHIYG+180.
          XA2=XE2-X0
          YA2=YE2-Y0
        END IF

        PHIY90=(PHIYT+90.)*BOGFAK
        XA2=XC(XA2,PHIY90,SPACE)
        YA2=YC(YA2,PHIY90,SPACE)
        XR(3)=XC(XR(2),PHIY90,1.+SPACE+SPACE)
        YR(3)=YC(YR(2),PHIY90,1.+SPACE+SPACE)
        XR(4)=XC(XR(1),PHIY90,1.+SPACE+SPACE)
        YR(4)=YC(YR(1),PHIY90,1.+SPACE+SPACE)
      END IF

      IF (OKZ) THEN
        DO 5 L=LEN(TEXTZ),1,-1
          LZ=L
          IF (TEXTZ(L:L) .NE. ' ') GOTO 6
    5   continue
        LZ=0
        YE=Y0+4.*PP14
    6   CALL GRCHRC(PP14,90.,16)
        IF (LZ .GT. 0)
     *     CALL GQTXX(IWKID,X0,Y0+4.*PP14,TEXTZ(:LZ),IRC,XE,YE,XR,YR)
C       WEGEN PHIZ=90 GRAD VEREINFACHT SICH DIE MINIMAX-ABFRAGE
      END IF
C
************************************************************************
*     Z(*,1) ENTHAELT (-1,-1,-1)  |   IN DER X-Z-EBENE,                *
*     Z(*,2) ENTHAELT (-1,1,-1)   |   NACH DEN DREHUNGEN               *
*     Z(*,3) ENTHAELT (1,-1,-1)   |   UM PHI UND THETA                 *
*     Z(*,4) ENTHAELT (1,1,-1)    |                                    *
*     ZM(*,1) ENTHAELT DIE LINKE  (UNTERE) ECKE   (ZL)                 *
*     ZM(*,2) ENTHAELT DIE UNTERE          ECKE   (ZU)                 *
*     ZM(*,3) ENTHAELT DIE RECHTE (UNTERE) ECKE   (ZR)                 *
************************************************************************
  300 Z(1,1)=-COSP-SINP
      Z(2,1)=SINP*SINT-COSP*SINT-COST
      Z(1,2)=-COSP+SINP
      Z(2,2)=SINP*SINT+COSP*SINT-COST
      Z(1,3)=COSP-SINP
      Z(2,3)=-SINP*SINT-COSP*SINT-COST
      Z(1,4)=COSP+SINP
      Z(2,4)=-SINP*SINT+COSP*SINT-COST
      DO 311 K=1,3
        DO 310 I=1,2
          ZM(I,K)=Z(I,1)
 310    CONTINUE
 311  CONTINUE
      DO 320 K=2,4
        IF (Z(2,K) .LT. ZU(2)) THEN
          ZU(1)=Z(1,K)
          ZU(2)=Z(2,K)
          END IF
        IF (Z(1,K) .LT. ZL(1)) THEN
          ZL(1)=Z(1,K)
          ZL(2)=Z(2,K)
        ELSE IF (Z(1,K) .GT. ZR(1)) THEN
          ZR(1)=Z(1,K)
          ZR(2)=Z(2,K)
          END IF
  320   CONTINUE
C
C     UMRECHNEN IN DIE AKTUELLEN EINHEITEN
C
      WLEN=BREITE/SQRT(3.)
      LMAS=WLEN*ABS(COS(THETA*BOGFAK))
      DO 331 K=1,3
        DO 330 I=1,2
          ZM(I,K)=BREITE*0.5+WLEN*0.5*ZM(I,K)
 330    CONTINUE
 331  CONTINUE
C
C     X-ACHSE
      IF (OKX) THEN
        IF (Z(2,1) .LT. Z(2,2)) THEN
          IA=1
        ELSE
          IA=2
          END IF
        IE=IA+2
        if (abs(iopt).eq.2) ia=3-ia
        IE=IA+2
C       UMRECHNEN IN DIE AKTUELLEN EINHEITEN
        XA=BREITE*0.5+WLEN*0.5*Z(1,IA)
        YA=BREITE*0.5+WLEN*0.5*Z(2,IA)
        XE=BREITE*0.5+WLEN*0.5*Z(1,IE)
        YE=BREITE*0.5+WLEN*0.5*Z(2,IE)
        if (abs(iopt).eq.2) then
           Ya=Ya+lmas
           Ye=Ye+lmas
           fa=-1.
        else
           fa=1.
        endif
        IF (IOPT .GT.  0) THEN
C         ACHSENBESCHRIFTUNG
          LINKS=XE .LT. XA
          if (abs(iopt).eq.2) links=.not.links
          IF (IX .GT. 0) THEN
            CALL GRAXLIN(XA,YA,XE,YE,drdmpa(6),drdmpa(7),LINKS,0)
          ELSE IF (IX .LT. 0) THEN
            SMIN=10.**drdmpa(6)
            SMAX=10.**drdmpa(7)
            CALL GRAXLOG(XA,YA,XE,YE,SMIN,SMAX,LINKS,0)
            END IF
          END IF
        IF (LX .GT. 0) THEN
          LMAX=WLEN*ABS(COS(PHIX))*0.9
          CALL GRCHRC(PP14,PHIXT,INTSYM)
          CALL GQTXX(IWKID,XA,YA,TEXTX(:LX),IRC,XE1,YE1,XR,YR)
          LTEXT=SQRT((XE1-XA)**2+(YE1-YA)**2)
          LPFEIL=LTEXT+6.*PP14
          IF (LPFEIL .LE. LMAX) THEN
            SIZE=PP14
          ELSE
            SIZE=LMAX/(LTEXT/PP14+6.)
            CALL GRCHRC(SIZE,PHIXT,INTSYM)
            CALL GQTXX(IWKID,XA,YA,TEXTX(:LX),IRC,XE1,YE1,XR,YR)
            LPFEIL=LMAX
            END IF
          XE1=XC(XE1,PHIXT*BOGFAK,6.)
          YE1=YC(YE1,PHIXT*BOGFAK,6.)
          IF (ABS(PHIXT-PHIXG) .LT. 1.) THEN
C           PFEIL UND TEXT HABEN DIE GLEICHE RICHTUNG
            XA=XA+0.5*(XE-XE1)
            YA=YA+0.5*(YE-YE1)
            TFAK=2.
          ELSE
C           PFEIL UND TEXT HABEN VERSCHIEDENE RICHTUNGEN
            XA=XE-0.5*(XE-XE1)
            YA=YE-0.5*(YE-YE1)
            TFAK=LTEXT/SIZE+2.
            END IF
          XA=XC(XA,PHIX90,-fa*XX*PP14/SIZE)
          YA=YC(YA,PHIX90,-fa*XX*PP14/SIZE)
          XT=XC(XA,PHIX,TFAK)
          YT=YC(YA,PHIX,TFAK)
          XT=XC(XT,PHIX90,-fa*(1.+SPACE2))
          YT=YC(YT,PHIX90,-fa*(1.+SPACE2))
          XXL=XA+LPFEIL*COS(PHIX)
          XYL=YA+LPFEIL*SIN(PHIX)
          CALL GRARRW(XA,YA,XXL,XYL,2.*SIZE,.3*SIZE,1)
          IF (LX .GT. 0) CALL GRTXT(XT,YT,LX,TEXTX)
          END IF
        END IF
C
C     Y-ACHSE
      IF (OKY) THEN
        IF (Z(2,1) .LT. Z(2,3)) THEN
          IA=1
        ELSE
          IA=3
        END IF
        if (abs(iopt).eq.2) ia=4-ia
        IE=IA+1
C       UMRECHNEN IN DIE AKTUELLEN EINHEITEN
        XA=BREITE*0.5+WLEN*0.5*Z(1,IA)
        YA=BREITE*0.5+WLEN*0.5*Z(2,IA)
        XE=BREITE*0.5+WLEN*0.5*Z(1,IE)
        YE=BREITE*0.5+WLEN*0.5*Z(2,IE)
        if (abs(iopt).eq.2) then
           YA=YA+lmas
           YE=Ye+lmas
           fa=-1.
        else
           fa=1.
        endif
        IF (IOPT .GT.  0) THEN
C         ACHSENBESCHRIFTUNG
          LINKS=XE .LT. XA
          if (abs(iopt).eq.2) links=.not.links
          IF (IY .GT. 0) THEN
            CALL GRAXLIN(XA,YA,XE,YE,drdmpa(8),drdmpa(9),LINKS,0)
          ELSE IF (IY .LT. 0) THEN
            SMIN=10.**drdmpa(8)
            SMAX=10.**drdmpa(9)
            CALL GRAXLOG(XA,YA,XE,YE,SMIN,SMAX,LINKS,0)
            END IF
          END IF
        IF (LY .GT. 0) THEN
          LMAX=WLEN*ABS(COS(PHIY))*0.9
          CALL GRCHRC(PP14,PHIYT,INTSYM)
          CALL GQTXX(IWKID,XA,YA,TEXTY(:LY),IRC,XE1,YE1,XR,YR)
          LTEXT=SQRT((XE1-XA)**2+(YE1-YA)**2)
          LPFEIL=LTEXT+6.*PP14
          IF (LPFEIL .LE. LMAX) THEN
            SIZE=PP14
          ELSE
            SIZE=LMAX/(LTEXT/PP14+6.)
            CALL GRCHRC(SIZE,PHIYT,INTSYM)
            CALL GQTXX(IWKID,XA,YA,TEXTY(:LY),IRC,XE1,YE1,XR,YR)
            LPFEIL=LMAX
            END IF
          XE1=XC(XE1,PHIYT*BOGFAK,6.)
          YE1=YC(YE1,PHIYT*BOGFAK,6.)
          IF (ABS(PHIYT-PHIYG) .LT. 1.) THEN
C           PFEIL UND TEXT HABEN DIE GLEICHE RICHTUNG
            XA=XA+0.5*(XE-XE1)
            YA=YA+0.5*(YE-YE1)
            TFAK=2.
          ELSE
C           PFEIL UND TEXT HABEN VERSCHIEDENE RICHTUNGEN
            XA=XE-0.5*(XE-XE1)
            YA=YE-0.5*(YE-YE1)
            TFAK=LTEXT/SIZE+2.
            END IF
          XA=XC(XA,PHIY90,-fa*YY*PP14/SIZE)
          YA=YC(YA,PHIY90,-fa*YY*PP14/SIZE)
          XT=XC(XA,PHIY,TFAK)
          YT=YC(YA,PHIY,TFAK)
          XT=XC(XT,PHIY90,-fA*(1.+SPACE2))
          YT=YC(YT,PHIY90,-fa*(1.+SPACE2))
          YXL=XA+LPFEIL*COS(PHIY)
          YYL=YA+LPFEIL*SIN(PHIY)
          CALL GRARRW(XA,YA,YXL,YYL,2.*SIZE,.3*SIZE,1)
          IF (LY .GT. 0) CALL GRTXT(XT,YT,LY,TEXTY)
          END IF
        END IF
C
C     Z-ACHSE
      IF (OKZ) THEN
        LMAX=lmas*0.9
        IF (LZ .GT. 0) THEN
          CALL GRCHRC(PP14,90.,INTSYM)
          CALL GQTXX(IWKID,0.,0.,TEXTZ(:LZ),IRC,XE1,YE1,XR,YR)
          LPFEIL=YE1+6.*PP14
          IF (LPFEIL .LE. LMAX) THEN
            SIZE=PP14
          ELSE
            SIZE=LMAX/(YE1/PP14+6.)
            CALL GRCHRC(SIZE,90.,INTSYM)
            CALL GQTXX(IWKID,0.,0.,TEXTZ(:LZ),IRC,XE1,YE1,XR,YR)
            LPFEIL=LMAX
            END IF
C         PFEIL RECHTS
          XA=ZR(1)+ZZ*PP14
          YA=ZR(2)+0.5*(LMAX-LPFEIL)
          CALL GRARRW(XA,YA,XA,YA+LPFEIL,2.*SIZE,.3*SIZE,1)
          XT=XA+2.*PP14
          YT=YA+(1.+SPACE)*SIZE
          IF (LZ .GT. 0) CALL GRTXT(XT,YT,LZ,TEXTZ)
C         PFEIL LINKS
          XA=ZL(1)-ZZ*PP14
          YA=ZL(2)+0.5*(LMAX-LPFEIL)
          CALL GRARRW(XA,YA,XA,YA+LPFEIL,2.*SIZE,.3*SIZE,1)
          XT=XA-PP14
          YT=YA+(1.+SPACE)*SIZE
          IF (LZ .GT. 0) CALL GRTXT(XT,YT,LZ,TEXTZ)
          END IF
        IF (IOPT .GT.  0) THEN
C         ACHSENBESCHRIFTUNG
          IF (IZ .GT. 0) THEN
            CALL GRAXLIN(ZL(1),ZL(2),ZL(1),ZL(2)+LMAs,drdmpa(10),
     *      drdmpa(11),.TRUE.,0)
            CALL GRAXLIN(ZR(1),ZR(2),ZR(1),ZR(2)+LMAs,drdmpa(10),
     *      drdmpa(11),.FALSE.,0)
          ELSE IF (IZ .LT. 0) THEN
            SMIN=10.**drdmpa(10)
            SMAX=10.**drdmpa(11)
            CALL GRAXLOG(ZL(1),ZL(2),ZL(1),ZL(2)+LMAs,SMIN,
     *      SMAX,.TRUE.,0)
            CALL GRAXLOG(ZR(1),ZR(2),ZR(1),ZR(2)+LMAs,SMIN,
     *      SMAX,.FALSE.,0)
            END IF
          END IF
        END IF
C
  900 CALL GRCHRC(PP14,PP15,INTSYM)



      contains 


      function XC(XSTART,ALPHA,F)  result(erg)
      REAL,intent(in) :: XSTART,ALPHA,F
      REAL erg
      
      erg = XSTART+COS(ALPHA)*F*PP(14)
      END function xc    
     

      function YC(YSTART,ALPHA,F)  result(erg)
      REAL,intent(in) :: YSTART,ALPHA,F
      REAL erg
      
      erg = YSTART+SIN(ALPHA)*F*PP(14)
      END function yc    
     
    
     
 9999 END subroutine grdrlg
