C@PROCESS NOGOSTMT NOSDUMP OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRSHD (X1,Y1,X2,Y2,D,TH,N1,N2)
      REAL X1(*), Y1(*), X2(*), Y2(*), D, TH
      INTEGER N1,N2

C---- 6.10.87, HE. GERLACH
C---- X1,Y1  X UND Y KOORDINATEN POLYGON 1
C---- X2,Y2  X UND Y KOORDINATEN POLYGON 2
C---- D      ABSTAND DER LINIEN (WERT BEZOGEN AUF GRSCLC, ALSO IN CM)
C---- TH     WINKEL DER LINIEN IN GRAD
C---- N1     ANZAHL DER PUNKTE POLYGON 1
C---- N2     ANZAHL DER PUNKTE POLYGON 2
C---- UPDATE 8.8.1990

      PARAMETER (PIFAC=6.2831853072/360)
      REAL TS(20)
      COMMON/GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRPP/
      DATA M1/0/

C---- BERECHNUNG DER KOEFFIZIENTEN FUER DREHUNG UND STAUCHUNG
      T1=TH*PIFAC
      C=COS(T1)
      S=SIN(T1)
      FAX=(PP(3)-PP(1))/(PP(7)-PP(5))
      DX=PP(1)-PP(5)*FAX
      FAY=(PP(4)-PP(2))/(PP(8)-PP(6))
      DY=PP(2)-PP(6)*FAY
      CFAX=C*FAX
      SFAY=S*FAY
      SFAX=S*FAX
      CFAY=C*FAY
      ADX=C*DX+S*DY
      ADY=C*DY-S*DX

C---- VON NUN AN CM-KOORDINATEN
      PP5=PP(5)
      PP6=PP(6)
      PP7=PP(7)
      PP8=PP(8)
      CALL GRSCLV(PP(1),PP(2),PP(3),PP(4))

C---- BERECHNE YMIN UND YMAX IN DEN UM TH GRAD GEDREHTEN CM-KOORDINATEN
      RYMIN=RY(X1(1),Y1(1))
      RYMAX=RYMIN
      DO 2 I=2,N1
         T1=RY(X1(I),Y1(I))
         RYMIN=MIN(T1,RYMIN)
         RYMAX=MAX(T1,RYMAX)
    2 CONTINUE

      DO 9 I=1,N2
         T1=RY(X2(I),Y2(I))
         RYMIN=MIN(T1,RYMIN)
         RYMAX=MAX(T1,RYMAX)
    9 CONTINUE

C---- IN DIESEN KOORDINATEN LAEUFT RYMIN AUFWAERTS BIS RYMAX
C---- DO WHILE( RYMIN.LE.RYMAX )
    3    M=0
         IP=0
         JU=99

         RXD=RX(X1(N1),Y1(N1))
         XD= RX(X2(N2),Y2(N2))
         RYD=RY(X1(N1),Y1(N1))
         YD= RY(X2(N2),Y2(N2))
         CALL GRSTCK (XD,YD,RXD,RYD,RYMIN,TS,M)

         RXD=RX(X1(1),Y1(1))
         XD= RX(X2(1),Y2(1))
         RYD=RY(X1(1),Y1(1))
         YD= RY(X2(1),Y2(1))
         CALL GRSTCK (XD,YD,RXD,RYD,RYMIN,TS,M)

         DO 5 I=2,N1
            RXD=RX(X1(I),Y1(I))
            XD= RX(X1(I-1),Y1(I-1))
            RYD=RY(X1(I),Y1(I))
            YD= RY(X1(I-1),Y1(I-1))
            CALL GRSTCK (XD,YD,RXD,RYD,RYMIN,TS,M)
    5    CONTINUE

         DO 8 I=2,N2
            RXD=RX(X2(I),Y2(I))
            XD= RX(X2(I-1),Y2(I-1))
            RYD=RY(X2(I),Y2(I))
            YD= RY(X2(I-1),Y2(I-1))
            CALL GRSTCK (XD,YD,RXD,RYD,RYMIN,TS,M)
    8    CONTINUE

         M=M/2*2
         IF (M.GT.0) THEN
            IF (M1.EQ.0) THEN
               M1=M+1
            ELSE
               M1=0
            ENDIF

C           ZURUECKDREHEN INS CM-KOORDINATENSYSTEM VON GRSCLC, ZEICHNEN
            DO 6 J=1,M
               I=IABS(M1-J)
               XO=C*TS(I)-S*RYMIN
               YO=S*TS(I)+C*RYMIN
               IF(IP.EQ.0) THEN
                 CALL GRJMP(XO,YO)
               ELSE IF(IP.EQ.99) THEN
                 CALL GRDRW(XO,YO)
               ENDIF
               IP=IP+JU
               JU=-JU
    6       CONTINUE

         ENDIF

         RYMIN=RYMIN+D
         IF (RYMIN.LT.RYMAX) GOTO 3
C---- CONTINUE ZU DO WHILE

      CALL GRSCLV(PP5,PP6,PP7,PP8)

      contains 
      
      function  RX(C1,C2) result(erg)
      REAL,intent(in) :: C1,C2
      REAL  erg
      
      erg =  C1*CFAX+C2*SFAY+ADX
      END function RX  


  
      function  RY(C1,C2) result(erg)
      REAL,intent(in) :: C1,C2
      REAL  erg
      
      erg =  -C1*SFAX+C2*CFAY+ADY
      END function RY   
          
 

      END subroutine grshd
      
C@PROCESS NOGOSTMT NOSDUMP OPT(3)
C.....SUBROUTINE   GRSTCK     V060      06/01/70  PRODUCT NUMBER 94342
C.....COPYRIGHT 1968 CALIFORNIA COMPUTER PRODUCTS
C.....CHANGED 6.11.87 AT KFA JUELICH BY G.GROTEN

      SUBROUTINE GRSTCK (XD,YD,RXD,RYD,RYMIN,TS,M)

      DIMENSION TS(*)

      XD=XD-RXD
      YD=YD-RYD

      IF (YD.NE.0E0) THEN
         RYD=RYMIN-RYD
    1    TEST=RYD/YD

         IF (TEST.GT.0E0 .AND. TEST.LT.1E0) THEN
            T=TEST*XD+RXD
            DO 2 J=1,M
               IF (TS(J).GE.T) THEN
                  A=TS(J)
                  TS(J)=T
                  T=A
               ENDIF
    2       CONTINUE
            M=M+1
            TS(M)=T
         ELSE IF (TEST.LT.0E0 .OR. TEST.GT.1E0) THEN
         ELSE
            RYMIN=RYMIN+.0001
            RYD=RYD+.0001
            GO TO 1
         ENDIF

      ENDIF


      END
