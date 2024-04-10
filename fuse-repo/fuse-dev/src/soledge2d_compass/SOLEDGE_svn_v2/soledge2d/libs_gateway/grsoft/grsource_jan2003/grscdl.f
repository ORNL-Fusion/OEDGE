C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C------------------------------------ SCRIPTED DOUBLE LINE
C N     - ANZAHL DER KURVENPUNKTE
C X(N)  - X-WERTE DER KURVEN
C Y(N)  - Y-WERTE DER KURVEN
C IDASH - 0 ODER 1, 2,3,4     DURCHGEZ., GESTRIC., PUNKT., STRICHPUNKT.
C       - +I*8, I=1,2,3,4     DOPPELLINIE, KOMBINIERT
C CHA   - CHARACTERSTRING IN DIE KURVE, WENN NICHT ' '
C NH    - LAENGE DER HILFSVEKTOREN, N+20 REICHT IMMER
C XH    - Hilfsvektor
C YH    - Hilfsvektor
C----------------------------------------------------------------------
C UPDATE 7.8.89 ; "ABS(DET) .LT.0.5E-2"    STATT "DET. EQ. 0."
C UPDATE 28.6.90; BUSCH, GSTXAL RICHTIG GESEZT FUER GTSGRAL GKS
C UPDATE  7.8.90; Groten, "IF (K.LT.2) GOTO 99"
C UPDATE  7.8.90; Groten, deutsche Umlaute mit GRTXTSCN
C UPDATE 10.12.90; BUSCH , PP(9) DURCH COMMON GRPIC ERSETZT
C UPDATE 10.9.91; Groten, alten CHARACTER_UP_Vector merken & ruecksetze
C UPDATE 20. 2.1991 Busch , GKS Treiber erlaubt max. 300 Werte pro Kurve
C                   GPL Aufruf durch GRLN Aufruf ersetzt

      SUBROUTINE GRSCDL(N,X,Y,IDASH,CHA,NH,XH,YH)
      REAL XH(NH),YH(NH),X(N),Y(N)
      CHARACTER(len=*) CHA

      PARAMETER (DEL=0.05,APRO=1.14265)
      DOUBLE PRECISION SUKR
      INTEGER IB(2),IT(2)
      INTEGER FLPIC , RAHMEN
      CHARACTER (len=256) TEXT1
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
       SAVE /GRPP/ ,/GRPIC/

      IF (N.LT.2) THEN
         WRITE(*,*) 'GRSCDL: N MUST BE GREATER THAN 1'
         STOP 10
      ENDIF
      K=1
      DO 1 I=2,N
         IF ( X(I).NE.X(I-1) .OR. Y(I).NE.Y(I-1) ) THEN
            K=K+1
            X(K)=X(I)
            Y(K)=Y(I)
         ENDIF
    1 CONTINUE
      IF (K.LT.2) GOTO 99
      IF (K.GT.NH) THEN
         WRITE(*,*) 'GRSCDL: NH MUST BE GREATER THAN OR EQUAL TO N'
         STOP 30
      ENDIF
      IDAS=MAX(IDASH,1)
      FAX=(PP(3)-PP(1))/(PP(7)-PP(5))
      FAY=(PP(4)-PP(2))/(PP(8)-PP(6))
      FAXY=FAX*FAY
      DELX=DEL/FAX
      DELY=DEL/FAY
      IF (CHA.EQ.' ') THEN
         NT=1
         IB(1)=1
         IT(1)=K
      ELSE
         ZALE=LEN(CHA)*PP(14)
C------  BERECHNE NAEHERUNGSWEISE DIE LAENGE EINES KURVENSTUECKS IN CM
C------  MAXIMALER RELATIVER FEHLER IST ETWA 1 %
         DL=0.
         DO 13 I=2,K
            A1=ABS(X(I)-X(I-1))*FAX
            B1=ABS(Y(I)-Y(I-1))*FAY
            DL=DL+A1+B1-APRO*A1*B1/(A1+B1)
   13    CONTINUE
C------- DL MITTLERE LAENGE EINES ABSCHNITTS IN CM
         DL=DL/(K-1)
         MAZA=ZALE/DL+0.95
         IF (MOD(MAZA,2).NE.0) MAZA=MAZA+1
         if (2*maza.gt.n) then
            NT=1
            IB(1)=1
            IT(1)=K
            goto 777
         endif
         LIRE=MAZA/2
         MAZAM1=MAZA-1
         EPSX=DL*.025/FAX
         EPSY=DL*.025/FAY
         IF (ABS(X(K)-X(1)).LE.EPSX.AND. ABS(Y(K)-Y(1)).LE.EPSY) THEN
C---------- ZYKLISCH
            IF (K+MAZAM1.GT.NH) THEN
               WRITE(*,*) 'GRSCDL: NH IS TO LOW, CAUSE P(1)=P(N)'
               WRITE(*,*) '        THE PICTURE IS OK IN MOST CASES'
               MAZAM1=NH-K
            ENDIF
            DO 4 I=2,MAZAM1
               XH(K+I-1)=X(I)
               YH(K+I-1)=Y(I)
    4       CONTINUE
            KS=K+MAZAM1
         ELSE
C---------- NICHT ZYKLISCH
            KS=K
         ENDIF
         DO 5 I=1,K
            XH(I)=X(I)
            YH(I)=Y(I)
    5    CONTINUE
C------- BERECHNUNG DER WINKELABWEICHUNG VON DER GERADEN
         X2=(XH(2)-XH(1))*FAX
         Y2=(YH(2)-YH(1))*FAY
         DO 3 I=2,KS-1
            X1=X2
            X2=(XH(I+1)-XH(I))*FAX
            Y1=Y2
            Y2=(YH(I+1)-YH(I))*FAY
            XI=-X2*Y1+Y2*X1
            YI= X2*X1+Y2*Y1
            IF (xi.ne.0. .or. yi.ne.0.) then
               XH(I)=ABS(ATAN2(XI,YI))
            else
               xh(i)=xh(i-1)
            endif
    3    CONTINUE
C------- BERECHNE MINIMUM DER MITTLEREN WINKELABWEICHUNG (POSITION)
         IMI=LIRE+1
         SUKR=0.
         DO 7 I=2,MAZA
            SUKR=SUKR+XH(I)
    7    CONTINUE
         SUMIN=SUKR
         DO 8 I=LIRE+2,KS-LIRE
            SUKR=SUKR-XH(I-LIRE)+XH(I+LIRE-1)
            IF (SUKR.LT.SUMIN) THEN
               IMI=I
               SUMIN=SUKR
            ENDIF
    8    CONTINUE
         ILI=IMI-LIRE
         IRE=IMI+LIRE
         IF (ILI.EQ.1) THEN
            IB(1)=IRE
            IT(1)=K
            NT=1
            IF (K-IRE .LT. 1) NT=0
            IA=1
            IE=IRE
         ELSEIF (IRE.GT.K) THEN
            IB(1)=IRE-K+1
            IT(1)=ILI
            NT=1
            IF (IT(1)-IB(1) .LT. 1) NT=0
            IA=IB(1)
            IE=ILI
         ELSEIF (IRE.EQ.K) THEN
            IB(1)=1
            IT(1)=ILI
            NT=1
            IF (IT(1)-IB(1) .LT. 1) NT=0
            IA=ILI
            IE=K
         ELSE
            IB(1)=1
            IT(1)=ILI
            NT=2
            IF (ILI .LT. 2) NT=1
            IB(NT)=IRE
            IT(NT)=K
            IF (NT.EQ.1 .AND. K-IRE .LT. 1) NT=0
            IA=ILI
            IE=IRE
         ENDIF
         IF (NT.NE.0 .AND. (X(IA).NE.X(IE) .OR. Y(IA).NE.Y(IE))) THEN
            XUP=(Y(IA)-Y(IE))*FAY
            YUP=(X(IE)-X(IA))*FAX
            IF (YUP.LT.0.) THEN
               XUP=-XUP
               YUP=-YUP
            ENDIF
            SPACE=SQRT((FAX*(X(IE)-X(IA)))**2+(FAY*(Y(IE)-Y(IA)))**2)
C---------- HOLE ALTES TEXTALIGNMENT ZUM SICHERN
            CALL GQTXAL(J,IALH,IALL)
C---------- SETZE TEXTALIGNMENT ZENTRIERT, HALB
C---------- CALL GSTXAL(2,7) FEHLER BEI TU BERLIN
            CALL GSTXAL(2,3)
C---------- HOLE ALTE TEXTHOEHE ZUM SICHERN
            CALL GQCHH(J,CHH)
C---------- SETZE TEXTHOEHE NEU
            CALL GSCHH(PP(14))
C---------- HOLE ALTEN CHARACTER_EXPANSION_FACTOR
            CALL GQCHXP(J,CXP)
C---------- SETZE CHARACTER_EXPANSION_FACTOR
            CALL GSCHXP(SPACE/(ZALE+1.))
C---------- HOLE ALTEn Character-up-Vektor
            CALL GQCHUP(J,CUX,CUY)
C---------- SETZE CHARACTER_UP_VECTOR
            CALL GSCHUP(XUP,YUP)
C--------   PP(9) = FCTR = 10.5 : XMAXCM=39.5, YMAXCM=28.7
C BUSCH 10.12.90
C           XMAXCM=3.7619048*PP(9)
C           YMAXCM=2.7333333*PP(9)
C---------  WINDOW AUF NDC (0#0,1#0.7266) ABBILDEN
            CALL GSWN(1,0.,XMAXCM,0.,YMAXCM)
            XM=(X(IA)+X(IE))*.5
            YM=(Y(IA)+Y(IE))*.5
            XM=PP(1)+(XM-PP(5))*FAX
            YM=PP(2)+(YM-PP(6))*FAY
            CALL GRTXSCN ( CHA, LEN(CHA), TEXT1 ,LTE1 )
            CALL GTX(XM,YM,TEXT1(:LTE1))
C-------    UMRECHNUNG DES WINDOW AUF DEN GESAMTEN VIEWPORTBEREICH
            XLI=PP(5)-PP(1)/FAX
            YLI=PP(6)-PP(2)/FAY
            XRE=PP(7)+(XMAXCM-PP(3))/FAX
            YRE=PP(8)+(YMAXCM-PP(4))/FAY
C---------- WINDOW AUF NDC (0#0,1#0.7266) ABBILDEN
            CALL GSWN(1,XLI,XRE,YLI,YRE)
C---------- SETZE ALTEN CHARACTER_UP_VECTOR zurueck
            CALL GSCHUP(CUX,CUY)
C---------- SETZE ALTEN CHARACTER_EXPANSION_FACTOR ZURUECK
            CALL GSCHXP(CXP)
C---------- SETZE ALTE TEXTHOEHE ZURUECK
            CALL GSCHH(CHH)
C---------- SETZE ALTES TEXTALIGNMENT ZURUECK
            CALL GSTXAL(IALH,IALL)
         ENDIF
      ENDIF
 777  DO 9 L=1,NT
         N1=IB(L)
         N2=IT(L)
         NN=N2-N1+1
         IF (IDAS.LT.8) THEN
            IF (IDAS.GT.4) IDAS=1
            CALL GSLN(IDAS)
            IF (NN.GT.1) THEN
               CALL GRLN(X(N1),Y(N1),NN)
            ENDIF
         ELSE
            A2=(Y(N1)-Y(N1+1))*FAY
            B2=(X(N1+1)-X(N1))*FAX
            A=ABS(A2)
            B=ABS(B2)
C---------  NAEHERUNG FUER SQRT(A**2+B**2)
C------  MAXIMALER RELATIVER FEHLER IST ETWA 1 %
            D2=A+B-APRO*A*B/(A+B)
            C2P=(Y(N1+1)*X(N1)-X(N1+1)*Y(N1))*FAXY+DEL*D2
            C2M=(Y(N1+1)*X(N1)-X(N1+1)*Y(N1))*FAXY-DEL*D2
            XH(N1)=X(N1)-DELX*A2/D2
            YH(N1)=Y(N1)-DELY*B2/D2
            X(N1)=X(N1)+DELX*A2/D2
            Y(N1)=Y(N1)+DELY*B2/D2
            DO 2 I=N1+1,N2-1
               A1=A2
               B1=B2
               C1P=C2P
               C1M=C2M
               D1=D2
               A2=(Y(I)-Y(I+1))*FAY
               B2=(X(I+1)-X(I))*FAX
               A=ABS(A2)
               B=ABS(B2)
C------------  NAEHERUNG FUER SQRT(A**2+B**2)
C------  MAXIMALER RELATIVER FEHLER IST ETWA 1 %
               D2=A+B-APRO*A*B/(A+B)
               C2P=(Y(I+1)*X(I)-X(I+1)*Y(I))*FAXY+DEL*D2
               C2M=(Y(I+1)*X(I)-X(I+1)*Y(I))*FAXY-DEL*D2
               DET=A1*B2-B1*A2
               IF ( ABS(DET).LT. 0.5E-2 ) THEN
                  XH(I)=X(I)-DELX*A2/D2
                  YH(I)=Y(I)-DELY*B2/D2
                  X(I) =X(I)+DELX*A2/D2
                  Y(I) =Y(I)+DELY*B2/D2
               ELSE
                  XH(I)=(B1*C2P-B2*C1P)/(DET*FAX)
                  YH(I)=(C1P*A2-A1*C2P)/(DET*FAY)
                  X(I) =(B1*C2M-B2*C1M)/(DET*FAX)
                  Y(I) =(C1M*A2-A1*C2M)/(DET*FAY)
               ENDIF
   2        CONTINUE
            XH(N2)=X(N2)-DELX*A2/D2
            YH(N2)=Y(N2)-DELY*B2/D2
            X(N2)=X(N2)+DELX*A2/D2
            Y(N2)=Y(N2)+DELY*B2/D2
            IDA=IDAS/8
            IF (IDA.GT.4) IDA=1
            CALL GSLN(IDA)
            CALL GRLN(X(N1),Y(N1),NN)
            IDA=MAX(MOD(IDAS,8),1)
            IF (IDA.GT.4) IDA=1
            CALL GSLN(IDA)
            CALL GRLN(XH(N1),YH(N1),NN)
         ENDIF
    9 CONTINUE
   99 END
