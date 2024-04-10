C@PROCESS OPT(3) IL(DIM) NOGOSTMT NOSDUMP
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C----------------------------------------------------------------------
C UPDATE 27.9.90 M. BUSCH    Daten werden formatfrei geschrieben
C                            frueher (20A4) d.h. CRAY unvertraeglich
C UPDATE 16.12.90 Groten     Abfragen auf ungueltige Daten entfernt
C                            Formate zu den I/O Anweisungen getan
C                            PARAMETER,SAVE,IFTHENELSE  modern FORTRAN
C UPDATE 29. 9.91 Groten
C----------------------------------------------------------------------
      SUBROUTINE KURVEF(X,Y,IST,ISY)
C     AUTOR: GERD GROTEN
C       X     BEREICH DER ABSZISSEN
C       Y     BEREICH DER ORDINATEN
C       IST   ANZAHL DER STUETZSTELLEN
C       ISY   SYMBOL UND FARBE DES ZEICHENSTIFTES
C
      INTEGER IST,IO,ISY,LB,KB,JJ,JJJ,ISATZ,J,JSY
      PARAMETER (IO=91)
      REAL   X(IST),Y(IST),ZW(100),ZU(100),CILBER,XMIN,XMAX,YMIN,YMAX
      COMMON /GRCIL/ CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
      SAVE ZW, /GRCIL/
      DATA   ZW /100*0./

      IF (IST.LT.1) THEN
         WRITE(*,*) 'KURVEF:  IST < 1'
         GOTO 9999
      ENDIF

      JSY=MOD(ISY,2000)
      IF (JSY.GT.1318                     .OR.
     *    JSY.GT. 118 .AND. JSY.LT.200    .OR.
     *    JSY.GT. 318 .AND. JSY.LT.400    .OR.
     *    JSY.GT. 518 .AND. JSY.LT.600    .OR.
     *    JSY.GT. 718 .AND. JSY.LT.800    .OR.
     *    JSY.GT. 918 .AND. JSY.LT.1000   .OR.
     *    JSY.GT.1118 .AND. JSY.LT.1200   .OR.
     *    ISY.LT.-32                        ) THEN
         WRITE(*,*) 'KURVEF:  ISY OUT OF RANGE ',ISY
         GOTO 9999
      ENDIF

      IF (CILBER(3).EQ.0) THEN
         XMIN=X(1)
         XMAX=XMIN
         YMIN=Y(1)
         YMAX=YMIN
      ELSE
         XMIN=CILBER(4)
         XMAX=CILBER(5)
         YMIN=CILBER(6)
         YMAX=CILBER(7)
      ENDIF
      CILBER(3)=CILBER(3)+1

      DO 30 J=1,IST
         XMIN=AMIN1(XMIN,X(J))
         XMAX=AMAX1(XMAX,X(J))
         YMIN=AMIN1(YMIN,Y(J))
         YMAX=AMAX1(YMAX,Y(J))
   30 CONTINUE
      CILBER(4)=XMIN
      CILBER(5)=XMAX
      CILBER(6)=YMIN
      CILBER(7)=YMAX

      ISATZ=(IST-1)/100+1
      LB=1
      KB=MIN(100,IST)
      DO 50 J=1,ISATZ
         JJJ=1
         DO 40 JJ=LB,KB
            ZW(JJJ)=X(JJ)
            ZU(JJJ)=Y(JJ)
            JJJ=JJJ+1
   40    CONTINUE
         WRITE(IO)IST,ISY,ZW,ZU
         LB=LB+100
         KB=MIN(KB+100,IST)
   50 CONTINUE

 9999 END
