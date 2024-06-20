C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRAXSL(IX0,IY0,I10)
************************************************************************
*                                                                      *
*     LOGARITHMISCHES ACHSENPROGRAMM                                   *
*     AUTOR: GERD GROTEN, KFA JUELICH, ZAM                             *
*                    UPDATE 29.9.89 GROTEN (LETTER)                    *
*                    UPDATE 18.6.91 GROTEN:GSTXAL ->GTSGRAL            *
*                    UPDATE 15.7.91 GROTEN:Markierung innen,Netzlinien *
*                    UPDATE  5.5.92 GROTEN:mehr Unterteilungen,Netz log*
*           LETZTES UPDATE 10.10.95 GROTEN:PP integer :Fehler Power WS *
*                                                                      *
*     IX0: STEUERGROESSE FUER DIE X-ACHSE:                             *
*          +/-1: X-ACHSE UNTEN                                         *
*          +/-2: X-ACHSE OBEN                                          *
*          +/-3: X-ACHSE UNTEN UND OBEN                                *
*          +/-4: X-ACHSE UNTEN           nur Decaden                   *
*          +/-5: X-ACHSE OBEN            nur Decaden                   *
*          +/-6: X-ACHSE UNTEN UND OBEN  nur Decaden                   *
*             0: KEINE X-ACHSE                                         *
*            >0: BESCHRIFTUNG MIT ZAHLEN                               *
*            <0: KEINE BESCHRIFTUNG                                    *
*          absoluter Zahlenwert um 10 erhoeht: Markierungsstriche innen*
*          absoluter Zahlenwert um 100 erhoeht: Netzlinien Dekaden     *
*          absoluter Zahlenwert um 1000 erhoeht: Netzl. auch Zwischenw.*
*     IY0: STEUERGROESSE FUER DIE Y-ACHSE:                             *
*          +/-1: Y-ACHSE LINKS                                         *
*          +/-2: Y-ACHSE RECHTS                                        *
*          +/-3: Y-ACHSE LINKS UND RECHTS                              *
*          +/-4: Y-ACHSE LINKS            nur Dekaden                  *
*          +/-5: Y-ACHSE RECHTS           nur Dekaden                  *
*          +/-6: Y-ACHSE LINKS UND RECHTS nur Dekaden                  *
*             0: KEINE Y-ACHSE                                         *
*            >0: BESCHRIFTUNG MIT ZAHLEN                               *
*            <0: KEINE BESCHRIFTUNG                                    *
*          absoluter Zahlenwert um 10 erhoeht: Markierungsstriche innen*
*          absoluter Zahlenwert um 100 erhoeht: Netzlinien Dekaden     *
*          absoluter Zahlenwert um 1000 erhoeht: Netzl. auch Zwischenw.*
*     I10: FALLS I10>0, DANN WERDEN DIE ACHSEN IN                      *
*          VOLLE DEKADEN UNTErTEILT.                                   *
*                                                                      *
*     DIE RICHTUNG, IN DIE EINE ACHSE GEZEICHNET WERDEN                *
*     SOLL, WIRD LOGARITHMISCH MIT GRSCLV SKALIERT.                    *
*                                                                      *
************************************************************************
      CHARACTER (len=3) ZAHL
      CHARACTER (len=5) ZAHLAN
      COMMON /GRLOGY/ ISLOGY,YINNEN,YNETZ,YL,DSH1,DSH2,DSH3
CDEC$ PSECT /GRLOGY/ NOSHR

      INTEGER PP0(18), PP(18)
      COMMON /GRPP/ PP0
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRLOGY/ , /GRPP/
      REAL FL(8)
      LOGICAL X2X5,INTERX,INTERY,XINNEN,YINNEN,XNETZ,YNETZ,XL,YL
      CHARACTER(len=8) ZL
      LOGICAL FLGROT
      EQUIVALENCE (PP(1),PP1), (PP(2),PP2), (PP(3),PP3), (PP(4),PP4),
     $ (PP(5),PP5), (PP(6),PP6), (PP(7),PP7), (PP(8),PP8),
     $ (PP(9),IFONT), (PP(10),PP10), (PP(11),PP11), (PP(13),INTLIN),
     $ (PP(14),CSIZE), (PP(16),ICOLOR), (PP(17),PP17), (PP(18),FLGROT)
      DATA FL/.30103004,.698970079,.47712125,.602059991,.77815125,
     $        .845098040,.90308999,.954242509/
      DATA ZL/'25346789'/

C
C     PP5, PP6, PP7, PP8 SIND NICHT-LOGARITHMIERTE WERTE
C
      CALL GQTXAL(IER,ITX0,ITY0)
      IX=IABS(IX0)
      IF ( IX.GE.1000 ) THEN
         IX=MOD(IX,1000)
         XL=.TRUE.
         XNETZ=.TRUE.
      ELSE
         XL=.FALSE.
         XNETZ=.FALSE.
      ENDIF
      IF ( IX.GE.100 ) THEN
         IX=MOD(IX,100)
         XNETZ=.TRUE.
      ENDIF
      IF ( IX.GE.10 ) THEN
         IX=MOD(IX,10)
         XINNEN=.TRUE.
      ELSE
         XINNEN=.FALSE.
      ENDIF
      IF (IX.GT.3) THEN
         INTERX=.FALSE.
         IX=IX-3
      ELSE
         INTERX=.TRUE.
      ENDIF
      IY=IABS(IY0)
      IF ( IY.GE.1000 ) THEN
         IY=MOD(IY,1000)
         YL=.TRUE.
         YNETZ=.TRUE.
      ELSE
         YL=.FALSE.
         YNETZ=.FALSE.
      ENDIF
      IF ( IY.GE.100 ) THEN
         IY=MOD(IY,100)
         YNETZ=.TRUE.
      ENDIF
      IF ( IY.GE.10 ) THEN
         IY=MOD(IY,10)
         YINNEN=.TRUE.
      ELSE
         YINNEN=.FALSE.
      ENDIF
      LETTER=1
      IF (IY0.LT.0) LETTER=2
      IF (IY.GT.3) THEN
         INTERY=.FALSE.
         IY=IY-3
      ELSE
         INTERY=.TRUE.
      ENDIF
      DO 5 I=1,18
        PP(I)=PP0(I)
    5 CONTINUE
      CALL GRDSH(1.,0.,1.)
      YDIF=PP4-PP2
      XDIF=PP3-PP1
      IF (IX .GT. 0) THEN
************************************************************************
*       X-ACHSE                                                        *
************************************************************************
        XMINL=ALOG10(PP5)
        IXMIN=NINT(XMINL)
        IF (10.**(IXMIN-0.001) .GT. PP5) IXMIN=IXMIN-1
        IF (I10 .GT. 0) XMINL=REAL(IXMIN)
        XMAXL=ALOG10(PP7)
        IXMAX=NINT(XMAXL)
        IF (10.**(IXMAX+0.001) .LT. PP7) IXMAX=IXMAX+1
        IF (I10 .GT. 0) XMAXL=REAL(IXMAX)
        CALL GRSCLV(XMINL,0.,XMAXL,YDIF)
C       ABSTANDSKONTROLLE
        ABST=XDIF/(XMAXL-XMINL)
        X2X5= (ABST .GE. 4.*CSIZE) .AND. INTERX
        ISTEP=5.*CSIZE/ABST
        IF (ISTEP .EQ. 0) ISTEP=1
        IXMIN=(IXMIN-(ISTEP-1))/ISTEP*ISTEP
        ML = 2
        IF ( REAL(IXMAX-IXMIN)/ISTEP .LE. 2. ) ML = 8
        XM1=CSIZE*0.5
        XM2=CSIZE
        IF (XINNEN) THEN
           XM1=-XM1
           XM2=-XM2
           TM1=CSIZE*0.6
           TM2=CSIZE*0.6
        ELSE
           TM1=CSIZE*0.6
           TM2=CSIZE*1.1
        ENDIF
        IF (IX .NE. 2) THEN
C
C         X-ACHSE UNTEN
C
          IF (XINNEN) THEN
             Y1=0.
             Y2=-XM2
             IF ( IX.EQ.3 ) THEN
                Y3=YDIF+XM2
             ELSE
                Y3=YDIF
             ENDIF
             TY1=-CSIZE*.65
             TY2=-CSIZE*.5
             TY3=-CSIZE*.35
          ELSE
             Y1=-XM2
             Y2=0.
             Y3=YDIF
             TY1=-CSIZE*.65-XM2
             TY2=-CSIZE*.5-XM1
             TY3=-CSIZE*.35-XM2
          ENDIF
          CALL GRJMP(XMINL,0.)
          CALL GRDRW(XMAXL,0.)
          CALL GRCHRC(CSIZE,0.,16)
          CALL GSTXAL(2,1)
          DO 10 I=IXMIN,IXMAX,ISTEP
             X=REAL(I)
             IF (.NOT. INTER(XMINL,X,XMAXL)) GOTO 10
             CALL GRJMP(X,Y1)
             CALL GRDRW(X,Y2)
             IF ( XNETZ .AND. X.NE.XMINL .AND. X.NE.XMAXL ) THEN
                CALL GRDSH(PP10,PP11,PP12)
                CALL GRJMP(X,Y2)
                CALL GRDRW(X,Y3)
                CALL GRDSH(1.,0.,1.)
             ENDIF
             IF (IX0 .GT. 0) THEN
                CALL GRTXT(X,TY1,2,'10')
             ENDIF
   10     CONTINUE
          CALL GRCHRC(CSIZE*0.67,0.,16)
          DO 20 I=IXMIN,IXMAX,ISTEP
             X=REAL(I)
             IF (X2X5) THEN
                CALL GSTXAL(2,1)
                DO 55 J=1,ML
                   IF (INTER(XMINL,X+FL(J),XMAXL)) THEN
                      CALL GRJMP(X+FL(J),0.)
                      CALL GRDRW(X+FL(J),-XM1)
                      IF (XL) THEN
                         CALL GRDSH(PP10,PP11,PP12)
                         CALL GRJMP(X+FL(J),Y2)
                         CALL GRDRW(X+FL(J),Y3)
                         CALL GRDSH(1.,0.,1.)
                      ENDIF
                      IF (IX0 .GT. 0) THEN
                         CALL GRTXT(X+FL(J),TY2,1,ZL(J:J))
                      ENDIF
                   ENDIF
   55           CONTINUE
             ENDIF
             IF (IX0.GT.0 .AND. INTER(XMINL,X,XMAXL)) THEN
                WRITE (ZAHL,'(I3)') I
                DO 88 IPOS=1,3
                   IF ( ZAHL(IPOS:IPOS) .NE.' ' ) GOTO 89
   88           CONTINUE
   89           ZAHLAN='  '//ZAHL(IPOS:)
                CALL GSTXAL(1,1)
                CALL GRTXT(X,TY3,-1,ZAHLAN)
             ENDIF
   20     CONTINUE
          XNETZ=.FALSE.
          XL=.FALSE.
       ENDIF
       IF (IX .GT. 1) THEN
C
C         X-ACHSE OBEN
C
          IF (XINNEN) THEN
             Y1=YDIF
             Y2=YDIF+XM2
             TY1=YDIF+CSIZE*.25
             TY2=YDIF+CSIZE*.4
             TY3=YDIF+CSIZE*.9
          ELSE
             Y1=YDIF+XM2
             Y2=YDIF
             TY1=YDIF+CSIZE*.25+XM2
             TY2=YDIF+CSIZE*.4+XM1
             TY3=YDIF+CSIZE*.9+XM2
          ENDIF
          Y3=0.
          CALL GRJMP(XMINL,YDIF)
          CALL GRDRW(XMAXL,YDIF)
          CALL GRCHRC(CSIZE,0.,16)
          CALL GSTXAL(2,5)
          DO 11 I=IXMIN,IXMAX,ISTEP
             X=REAL(I)
             IF (.NOT. INTER(XMINL,X,XMAXL)) GOTO 11
             CALL GRJMP(X,Y1)
             CALL GRDRW(X,Y2)
             IF ( XNETZ .AND. X.NE.XMINL .AND. X.NE.XMAXL ) THEN
                CALL GRDSH(PP10,PP11,PP12)
                CALL GRJMP(X,Y2)
                CALL GRDRW(X,Y3)
                CALL GRDSH(1.,0.,1.)
             ENDIF
             IF (IX0 .GT. 0) THEN
                CALL GRTXT(X,TY1,2,'10')
             ENDIF
   11     CONTINUE
          CALL GRCHRC(CSIZE*0.67,0.,16)
          DO 21 I=IXMIN,IXMAX,ISTEP
             X=REAL(I)
             IF (X2X5) THEN
                CALL GSTXAL(2,5)
                DO 66 J=1,ML
                   IF (INTER(XMINL,X+FL(J),XMAXL)) THEN
                      CALL GRJMP(X+FL(J),YDIF)
                      CALL GRDRW(X+FL(J),YDIF+XM1)
                      IF (XL) THEN
                         CALL GRDSH(PP10,PP11,PP12)
                         CALL GRJMP(X+FL(J),Y2)
                         CALL GRDRW(X+FL(J),Y3)
                         CALL GRDSH(1.,0.,1.)
                      ENDIF
                      IF (IX0 .GT. 0) THEN
                         CALL GRTXT(X+FL(J),TY2,1,ZL(J:J))
                      ENDIF
                   ENDIF
   66           CONTINUE
             ENDIF
             IF (IX0.GT.0 .AND. INTER(XMINL,X,XMAXL)) THEN
                WRITE (ZAHL,'(I3)') I
                DO 78 IPOS=1,3
                   IF ( ZAHL(IPOS:IPOS) .NE.' ' ) GOTO 79
   78           CONTINUE
   79           ZAHLAN='  '//ZAHL(IPOS:)
                CALL GSTXAL(1,5)
                CALL GRTXT(X,TY3,-1,ZAHLAN)
             ENDIF
   21     CONTINUE
        ENDIF
      ENDIF
      IF (IY .GT. 0) THEN
************************************************************************
*       Y-ACHSE                                                        *
************************************************************************
        CALL GRCHRC(CSIZE,0.,16)
        YMINL=ALOG10(PP6)
        IYMIN=NINT(YMINL)
        IF (10.**(IYMIN-0.001) .GT. PP6) IYMIN=IYMIN-1
        IF (I10 .GT. 0) YMINL=REAL(IYMIN)
        YMAXL=ALOG10(PP8)
        IYMAX=NINT(YMAXL)
        IF (10.**(IYMAX+0.001) .LT. PP8) IYMAX=IYMAX+1
        IF (I10 .GT. 0) YMAXL=REAL(IYMAX)
        SMIN=10.**YMINL
        SMAX=10.**YMAXL
        IF (.NOT. INTERY) THEN
           SMIN=-SMIN
           SMAX=-SMAX
        ENDIF
        CALL GRSCLV(PP1,PP2,PP3,PP4)
        ISLOGY=1234567890
        DSH1=PP10
        DSH2=PP11
        DSH3=PP12
C       Y-ACHSE LINKS
        IF (IY .NE. 2) THEN
           CALL GRAXLOG(PP1,PP2,PP1,PP4,SMIN,SMAX,.TRUE.,LETTER)
           YNETZ=.FALSE.
           YL=.FALSE.
        ENDIF
C       Y-ACHSE RECHTS
        IF (IY .GT. 1)
     *    CALL GRAXLOG(PP3,PP2,PP3,PP4,SMIN,SMAX,.FALSE.,LETTER)
        ENDIF
        ISLOGY=0
************************************************************************
      IF (IX .GT. 0) THEN
         IF (IY .GT. 0) THEN
            CALL GRSCLV(XMINL,YMINL,XMAXL,YMAXL)
         ELSE
            CALL GRSCLV(XMINL,PP6,XMAXL,PP8)
         ENDIF
         CALL GRCHRC(CSIZE,PP15,0)
      ELSE
         IF (IY .GT. 0) THEN
            CALL GRSCLV(PP5,YMINL,PP7,YMAXL)
            CALL GRCHRC(CSIZE,PP15,0)
         ENDIF
      ENDIF
      CALL GSTXAL(ITX0,ITY0)
      CALL GRDSH(PP10,PP11,PP12)

      contains 
      
      function inter(X1,X,X2) result(erg)
      REAL,intent(in) :: X1,X,X2
      LOGICAL erg
      
      erg = X1-0.0001 <= X .AND. X <= X2+0.0001
      END function inter    
     
      END subroutine graxsl
