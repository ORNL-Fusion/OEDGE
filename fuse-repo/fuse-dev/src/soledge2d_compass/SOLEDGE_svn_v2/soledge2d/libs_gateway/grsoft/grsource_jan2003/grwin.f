C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRWIN(WX1,WY1,WX2,WY2,IDUMMY)
CCC   D. Bartel, Ausschnittbildung fuer ADMGDF, 3820,4250
CCC    9.12.87
CCC   11.7.91 M. Busch Modifikation fuer GTSGRAL-GKS
CCC           Luft wird nicht mehr abgeschnitten, wie es frueher war
CCC           Ausschnitt wird auf Ausgabegearet (GDDM WK,CGM, PLotter,
CCC           GKSM ) abgebildet
CCC           Drehen vorher mit GR90DG - sonst spaeter mit PLOT EXEC
CCC           Einziehen in GML - mit PLOT EXEC und 'Frame Output'
CCC
CCC           Achtung - laeuft nicht auf CRAY in GRLIB ( COMMON wird dort
CCC                      benutzt
CCC   UPDATE: 19.7.91 M. BUSCH
CCC   UPDATE: 15.4.92 M. BUSCH fuer GRWINC, GRWINV: fur Berechnung von
CCC                            WID und HEI wird nachgesehen, welches
CCC                            die laengere Seite ist

CCC   UPDATE: 27.5.92 M. BUSCH GRWIN ignorieren wenn GR90DG aktiv
CCC   UPDATE: 23.11.92 M.BUSCH GRWINC, GRWINV korrigiert
      INTEGER  FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $         IGRSHW , ERRUN, IGR3 , SHOWPR, IGR3PL
      REAL XMAXCM,YMAXCM, XDCPIC,YDCPIC
      LOGICAL FLGROT
      COMMON /GRWIND/ X1,X2,Y1,Y2,IWIN
CDEC$ PSECT /GRWIND/ NOSHR
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON / GRCOM1/ IGRSHW , ERRUN, IGR3 , SHOWPR, IGR3PL
CDEC$ PSECT / GRCOM1/ NOSHR

      SAVE /GRWIND/ ,/GRPP/,/GRPIC/,  / GRCOM1/
      EQUIVALENCE (PP(18),FLGROT)
      idummy=idummy

      IF (FLGROT) goto 99
      TX1=WX1 * XDCPIC
      TX2=WX2 * XDCPIC
      TY1=WY1 * YDCPIC
      TY2=WY2 * YDCPIC

      GOTO 100
      ENTRY GRWINC(WX1,WY1,WX2,WY2,IDUMMY)
      IF (FLGROT) goto 99
c 23.11.92 WID immer xmaxcm, HEI immer ymaxcm
c     IF (XMAXCM.LT.YMAXCM) THEN
cc    IF (.NOT.FLGROT ) THEN
         WID=XMAXCM
         HEI=YMAXCM
c     ELSE
c        WID=YMAXCM
c        HEI=XMAXCM
c     ENDIF
c 23.11.92 Neu: multipl. mit Bildverhaeltnis xdcpic bzw. ydcpic
      TX1=xdcpic*WX1/WID
      TX2=xdcpic*WX2/WID
      TY1=ydcpic*WY1/HEI
      TY2=ydcpic*WY2/HEI
      GOTO 100
      ENTRY GRWINV(WX1,WY1,WX2,WY2,IDUMMY )
      IF (FLGROT) goto 99
      TX1=PP(1)+(PP(3)-PP(1))*(WX1-PP(5))/(PP(7)-PP(5))
      TY1=PP(2)+(PP(4)-PP(2))*(WY1-PP(6))/(PP(8)-PP(6))
      TX2=PP(1)+(PP(3)-PP(1))*(WX2-PP(5))/(PP(7)-PP(5))
      TY2=PP(2)+(PP(4)-PP(2))*(WY2-PP(6))/(PP(8)-PP(6))
c 23.11.92 WID immer xmaxcm, HEI immer ymaxcm
cc    IF (.NOT.FLGROT ) THEN
c     IF (XMAXCM.LT.YMAXCM) THEN
         WID=XMAXCM
         HEI=YMAXCM
c     ELSE
c        WID=YMAXCM
c        HEI=XMAXCM
c     ENDIF
c 23.11.92 Neu: multipl. mit Bildverhaeltnis xdcpic bzw. ydcpic
      TX1=xdcpic*TX1/WID
      TX2=xdcpic*TX2/WID
      TY1=ydcpic*TY1/HEI
      TY2=ydcpic*TY2/HEI
100   CONTINUE
      IF(  TX1 .LT. 0. .OR. TX1 .GT. 1.
     1.OR. TX2 .LT. 0. .OR. TX2 .GT. 1.
     2.OR. TY1 .LT. 0. .OR. TY1 .GT. 1.
     3.OR. TY2 .LT. 0. .OR. TY2 .GT. 1.  ) THEN
         WRITE(6,*)' GRWIN: Auschnittsparameter ausserhalb'//
     1   ' Geltungsbereich!'
         WRITE(6,*)'        Es werden die Maximalwerte eingesetzt.'
C        CALL OBEY('CP SLEEP 4 SEC')
        TX1  =0.
        TX2  =1.
        TY1  =0.
        TY2  =1.
      ENDIF
      X1=TX1
      X2=TX2
      Y1=TY1
      Y2=TY2
      IWIN=1234567890
      goto 999
99    write(errun,*) 'GRWIN Aufrufe werden ignoriert, wenn das Bild'
      write(errun,*) 'mit GR90DG gdreht wurde'
      write(errun,*)
     $        'Ein ADMGDF Graphikfile wird auf dem Standard GDDM '
      write(errun,*) 'Fenster erstellt'
      write(errun,*)
     $        'Zum Einziehen in GML kann man evtl. Leerflaechen am Rand'
      write(errun,*)
     $        'des Bildes mit Hilfe der PLOT Exec abschneiden'
999   RETURN
      END
