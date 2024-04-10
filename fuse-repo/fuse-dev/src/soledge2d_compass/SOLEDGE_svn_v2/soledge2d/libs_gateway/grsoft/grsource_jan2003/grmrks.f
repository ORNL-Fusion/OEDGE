C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRMRKS(HOEHE)
      REAL    HOEHE

*KFA 30.7.90 HG CHARACTERGROESSE MARKER ZURUECKGESETZT (2.0)
*KFA  2.8.90 BUSCH MARKERSIZE  AUF 1. GESETZT, MUSS NOCH
*Kfa  2) nominale Markersize ist auf allen Geraeten 2mm
*        Der Marker Size Scale Faktor wird in Abhaengigkeit
*        der HOEHE gesetzt.
* UPDATE 15.5.91 BUSCH
C UPDATE 10. 3.1992 Busch , GKS Treiber erlaubt max. 300 Marker pro
C                   Aufruf
C                   Ausgabe evtl. vorhandener Draw + Jump Daten
C----------------------------------------------------------------------

      COMMON /SCALE/ XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR
C---- COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /MKRTP/ MKTYP
CDEC$ PSECT /MKRTP/ NOSHR
      SAVE /SCALE/ ,/GRPP/,/GRREST/,/MKRTP/
C------- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

         IF (NRLN.GT.1) THEN
            CALL Gpl (NRLN,XLN,YLN)
            XCURR=0.
            YCURR=0.
            NRLN=0
         ENDIF

C------- AUSGABE EVT. VORHANDENER GRJMPS,.. -DATEN

         IF (NRMR.GT.0) THEN
            CALL GPM (NRMR,XMR,YMR)
            NRMR=0
         ENDIF

      PP(17) =HOEHE

C---- SET MARKER SIZE SCALE FACTOR
C Busch
C     Marker nominal 0.02 cm auf allen Ausgabegeraeten
C     Division durch 2 ergibt eine Makersize HOEHE auf allen WK's

      SIZEMK =  HOEHE / 0.2
      CALL GSMKSC( SIZEMK )


      END
