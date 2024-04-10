C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
c     Update: 21. 2.1991 Busch
      SUBROUTINE GRSHOW
      WRITE(6,*)'******************************************************'
      WRITE(6,*)'***           GR-Software  Warning                 ***'
      WRITE(6,*)'*** GRSHOW steht interaktiv nicht zur Verfuegung   ***'
      WRITE(6,*)'******************************************************'
      RETURN
      ENTRY GRDEL
      WRITE(6,*)'******************************************************'
      WRITE(6,*)'***           GR-Software  Warning                 ***'
      WRITE(6,*)'*** GRDEL steht interaktiv nicht zur Verfuegung    ***'
      WRITE(6,*)'******************************************************'
      END
      subroutine grdcur(np,xp,yp)

C     Unix Version
C     Stand: 17.2.94

      real xp(np),yp(np)

      real wi(4), vi(4)
      REAL XLN(300), YLN(300), YMR(300), XMR(300)
      REAL MOUT(2,3),pp(18)
      CHARACTER GKSTYP*7
      LOGICAL FLGROT
CCC
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON /grpp/ PP
CDEC$ PSECT /grpp/ NOSHR
      COMMON /grgks/ gkstyp
CDEC$ PSECT /grgks/ NOSHR
      COMMON /SCALE/ XMAXDC, XUNITS, YUNITS
CDEC$ PSECT /SCALE/ NOSHR
      COMMON /GRREST/ MAXPKT, NRLN, XCURR, YCURR, XLN, YLN, NRMR, XMR,
     1       YMR
CDEC$ PSECT /GRREST/ NOSHR
      EQUIVALENCE (PP(16),ICOLOR),(PP(17),SIZMRK),(PP(18),FLGROT)

      SAVE /GRPIC/,/grpp/ , /grgks/ ,/SCALE/,/GRREST/
      if (GKSTYP.EQ.'XGKS'.or.GKSTYP.EQ.'GRALGKS') THEN

      WRITE(6,*)'******************************************************'
      WRITE(6,*)'***           GR-Software  Warning                 ***'
      WRITE(6,*)'***           XGKS, GRALGKS                        ***'
      WRITE(6,*)'*** GRSHOW steht interaktiv nicht zur Verfuegung   ***'
      WRITE(6,*)'******************************************************'
      goto 999
      endif

CCC   AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN
CCC
      IF (NRLN .GT. 1) THEN
        CALL GPL(NRLN, XLN, YLN)
        XCURR = XLN(NRLN)
        YCURR = YLN(NRLN)
        NRLN = 0
      END IF
CCC   AUSGABE EVT. VORHANDENER GRJMPS-DATEN
CCC
      IF(NRMR.GT.0) THEN
        CALL GPM(NRMR,XMR,YMR)
        NRMR=0
      ENDIF


CCC     WENN FLGROT='Y', DANN MUSS DAS BILD GEDREHT WERDEN:
      IF (FLGROT) CALL GEVTM(0.5, 0.5, 0., 0., 3.141592*0.5, 1., 1., 1,
     1   MOUT)

       CALL GQMK(IERR,MTYPE)

C 17. 2.94
C--- offene wk's
      CALL GQOPWK(1,ierr,nopn,iwk)
      do 1 ii=1,nopn
         CALL GQOPWK(ii,ierr,nopn,iwk)
         CALL GQWKC ( IWK, IERR,ICON, IWT)
         CALL GQWKCA(IWT,IERR,ICAT)
C---nur vom Typ in/out
         if ( icat.ne.2)  goto 1
         goto 111
1     continue
      goto 999

CCC   UMRECHNUNG DES WINDOW AUF DEN GESAMTEN VIEWPORTBEREICH
C Weltkoordinaten
111    XLI=PP(5)-PP(1)*XUNITS
       YLI=PP(6)-PP(2)*YUNITS
       XRE=PP(7)+(XMAXCM-PP(3))*XUNITS
       YRE=PP(8)+(YMAXCM-PP(4))*YUNITS
C Bildflaeche auf gesammten Viewport bezogen umrechnen ( NT 1)
       call gqnt(1,ierr,wi,vi)
       xdiff = vi(2)-vi(1)
       ydiff = vi(4)-vi(3)
       if ( xdiff .gt. ydiff) then
          xdc=1
          ydc= ydiff / xdiff
       else
          ydc=1
          xdc= xdiff / ydiff
       endif
       CALL GSMK(3)
C Setzen von Markergroesse und -farbe fuer graphischen Cursor
       SIZE=SIZMRK
       ICOL=ICOLOR
       CALL GRMRKS(0.2)
       CALL GRNWPN(1)

       do 100 i=1,np

C beim GLIGKS kommt hier immer Wert mit Transformation 0 zurueck
C GSVPIP nicht realisiert- daher hier selbst umrechnen

       call grqlc(1,1,istat,it,xx,yy)

        XIN=XX
        YIN=YY
c
       if ( xdiff .gt. ydiff ) then
          yy=yy/ydc
       else
          xx=xx/xdc
       endif
       xx=XLI+(XRE-XLI)*xx
       yy=YLI+(YRE-YLI)*yy

       xp(i)=xx
       yp(i)=yy


C     Modifikation von Herrn Kleefeld, AVR REW085 eingebaut
C     Modifikation: jetzt GRDCUR auch bei GR90DG moeglich
c     Busch 25.10.91
         IF ( FLGROT ) THEN
            x Steig = ( PP(7) - PP(5) ) / ( PP(3) - PP(1) )
            y Steig = ( PP(8) - PP(6) ) / ( PP(4) - PP(2) )
            xcm           = ( XP(I)-PP(5) ) / x Steig + PP(1)
            ycm           = ( YP(I)-PP(6) ) / y Steig + PP(2)
            x1cm          = ycm
            y1cm          = XMAXCM - xcm
            XP(I) = ( x1cm - PP(1) ) * x Steig + PP(5)
            YP(I) = ( y1cm - PP(2) ) * y Steig + PP(6)
            CALL GPM(1,XP(I),YP(I))
c           CALL GPM(1,XIN ,YIN  )
            call guwk(1,1)
         ELSE
            CALL GPM(1,XP(I),YP(I))
            call guwk(1,1)
         END IF

100    continue
       CALL GSMK(MTYPE)
C Setzen von Markergroesse und -farbe fuer graphischen Cursor
      CALL GRMRKS(SIZE)
      CALL GRNWPN(ICOL)
999    end
