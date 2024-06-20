C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM) FIPS(F)
C ------------------ MIT OPT(3) NOSDUMP NOGOSTMT UEBERSETZEN ! ---------
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C ----------------------------------------------------------------------
C ----G R D R A X-------------------------------------------------------
C     UPDATE: 10. 1. 91 M. BUSCH FESTE CM ANGABE (39.3) DURCH WERT AUS
C                       COMMON GRPIC ERSETZT
C     UPDATE  27.10.94  G.Groten : bei GRFRBN farbe+2000*intensitaet
C ----------------------------------------------------------------------
      SUBROUTINE GRDRAX(LX,TXTX,LY,TXTY,LZ,TXTZ,BRTXT)

      INTEGER RAHMEN,FLPIC,ifa,idi
      CHARACTER(len=*) TXTX,TXTY,TXTZ
      CHARACTER (len=40) TEXT
C --- PP:= TABELLE DER PLOTPARAMETER (SIEHE GR-HANDBUCH)
      EQUIVALENCE (PP(16),INTSYM),(PP(13),INTLIN),(PP(16),ICOLOR)
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRDRA/ SIP,SIT,COP,FA,IRIFA
CDEC$ PSECT /GRDRA/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMXCM9,YMXCM9, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      SAVE /GRPP/,/GRDRA/,/GRPIC/
C --- SIP:= SIN(PHI)
C --- COP:= COS(PHI)
C --- SIT:= SIN(THETA)
C --- FA := PI/180
C --- IRIFA:= FARBE DES RICHTUNGSDREIBEINS

      XCM = XMXCM9 - XMXCM9/200.
      XTXT=(XCM-PP(1)+PP(3)-PP(1))*.5
C     XTXT=(39.3-PP(1)+PP(3)-PP(1))*.5
      YTXT=XTXT-(PP(3)-PP(1))-PP(2)
      AX=ATAN2(-SIP*SIT,COP)/FA
      AY=ATAN2( COP*SIT,SIP)/FA
      AZ=90.
      BRUR=PP(14)
      AUR=PP(15)
      INUR=INTSYM
      CALL GRCHRC(BRTXT,AX,18)
      idi = irifa/2000
      ifa = irifa-idi*2000
      jfa = icolor
      jdi = intlin
      CALL GRNWPN( IFA )
      CALL GRSPTS( idi)
      TEXT(1:2)='  '
      LTEXT=LX+2
      TEXT(3:LTEXT)=TXTX(1:LX)
      CALL GRTXT(XTXT,YTXT,LTEXT,TEXT)
      CALL GRCHRC(BRTXT,AY,18)
      LTEXT=LY+2
      TEXT(3:LTEXT)=TXTY(1:LY)
      CALL GRTXT(XTXT,YTXT,LTEXT,TEXT)
      CALL GRCHRC(BRTXT,AZ,18)
      LTEXT=LZ+2
      TEXT(3:LTEXT)=TXTZ(1:LZ)
      CALL GRTXT(XTXT,YTXT,LTEXT,TEXT)
      CALL GRCHRC(BRUR,AUR,INUR)
      CALL GRNWPN( jFA )
      CALL GRSPTS( jdi)
      END
