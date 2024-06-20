C
C
C----------------------------------------------------------------------*
C SUBROUTINE PLTAXI                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE EIRENE_PLTAXI(IERR,IAX)
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
      USE EIRMOD_CPLMSK
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IAX
      INTEGER, INTENT(OUT) :: IERR
      REAL(DP) :: EXPR, RI, ARG, RJ, T, PARAM1, PARAM,CENTRAL
      INTEGER :: FCTR, IL, IPARAM, I, J, EIRENE_IEXP10
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
            IF (FITX) CALL EIRENE_ANPSG(MINX,MAXX,INTNRX,STPSZX,IERR)
            IF (IERR.GT.0) RETURN
            CALL GRSPTS(20)
            CALL GRDSH(1.,0.,1.)
            CALL
     .  GRSCLV(REAL(MINX,KIND(1.E0)),0.,REAL(MAXX,KIND(1.E0)),
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
                  CALL
     .  GRDRW(REAL(T,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
                  CALL GRSPTS(20)
                  CALL GRDSH(1.,0.,1.)
               ELSE
                  CALL GRDRW(REAL(T,KIND(1.E0)),0.)
               ENDIF
    5       CONTINUE
            IF (.NOT.GRIDX) THEN
               CALL GRJMP(REAL(MINX,KIND(1.E0)),0.)
               CALL
     .  GRDRW(REAL(MINX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
               CALL GRJMP(REAL(MAXX,KIND(1.E0)),0.)
               CALL
     .  GRDRW(REAL(MAXX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
            ENDIF
C
C  LABELS AN DER X-ACHSE
C
            CALL GRCHRC(0.3,0.,16)
            DO 10 J=2,INTNRX+1,2
               CENTRAL=REAL((MINX+MAXX)*0.5,DP)
               IF (ABS(CENTRAL).LT.1.E-30)
     .           CENTRAL=MAX(ABS(MINX),ABS(MAXX))
               FCTR=EIRENE_IEXP10(CENTRAL)
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
            IF (FITX) CALL EIRENE_ANPSGL(MINX,MAXX,MINLX,MAXLX,IERR)
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
               CALL
     .  GRDRW(REAL(MAXLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
            ELSE
               CALL GRDRW(REAL(MAXLX,KIND(1.E0)),0.)
            ENDIF
C
            IF (.NOT.GRIDX) THEN
               CALL GRJMP(REAL(MINLX,KIND(1.E0)),0.)
               CALL
     .  GRDRW(REAL(MINLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
               CALL GRJMP(REAL(MAXLX,KIND(1.E0)),0.)
               CALL
     .  GRDRW(REAL(MAXLX,KIND(1.E0)),REAL(LENY,KIND(1.E0)))
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
                  CALL
     .  GRTXT(REAL(I,KIND(1.E0)),-0.6,2,CZHNLB(2:3))
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
              CALL EIRENE_ANPSG(MINY,MAXY,INTNRY,STPSZY,IERR)
            ELSE
              INTNRY=10
              STPSZY=(MAXY-MINY)/10.
            ENDIF
            IF (IERR.GT.0) RETURN
            CALL GRSPTS(20)
            CALL GRDSH(1.,0.,1.)
            CALL
     .  GRSCLV(0.,REAL(MINY,KIND(1.E0)),REAL(LENX,KIND(1.E0)),
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
                  CALL
     .  GRDRW(REAL(LENX,KIND(1.E0)),REAL(T,KIND(1.E0)))
                  CALL GRSPTS(20)
                  CALL GRDSH(1.,0.,1.)
               ELSE
                  CALL GRDRW(0.,REAL(T,KIND(1.E0)))
               ENDIF
  110       CONTINUE
            IF (.NOT.GRIDY) THEN
               CALL GRJMP(0.,REAL(MINY,KIND(1.E0)))
               CALL
     .  GRDRW(REAL(LENX,KIND(1.E0)),REAL(MINY,KIND(1.E0)))
               CALL GRJMP(0.,REAL(MAXY,KIND(1.E0)))
               CALL
     .  GRDRW(REAL(LENX,KIND(1.E0)),REAL(MAXY,KIND(1.E0)))
            ENDIF
C
C  LABELS AN DER Y-ACHSE
C
            CALL GRCHRC(0.3,90.,18)
            DO 115 J=2,INTNRY+1,2
               FCTR=EIRENE_IEXP10(REAL((MINY+MAXY)*0.5,DP))
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
           IF (FITY) CALL EIRENE_ANPSGL(MINY,MAXY,MINLY,MAXLY,IERR)
           IF (IERR.GT.0) RETURN
           CALL
     .  GRSCLV(0.,REAL(MINLY,KIND(1.E0)),REAL(LENX,KIND(1.E0)),
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
              CALL
     .  GRDRW(REAL(LENX,KIND(1.E0)),REAL(MAXLY,KIND(1.E0)))
           ELSE
              CALL GRDRW(0.,REAL(MAXLY,KIND(1.E0)))
           ENDIF
           IF (.NOT.GRIDY) THEN
              CALL GRJMP(0.,REAL(MINLY,KIND(1.E0)))
              CALL
     .  GRDRW(REAL(LENX,KIND(1.E0)),REAL(MINLY,KIND(1.E0)))
              CALL GRJMP(0.,REAL(MAXLY,KIND(1.E0)))
              CALL
     .  GRDRW(REAL(LENX,KIND(1.E0)),REAL(MAXLY,KIND(1.E0)))
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
                 CALL
     .  GRTXT(-0.85,REAL(I,KIND(1.E0))+0.05,2,CZHNLB(2:3))
              ELSE
                 CALL
     .  GRTXT(-0.85,REAL(I,KIND(1.E0))+0.05,3,CZHNLB)
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
