C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C----------------------------------------------------------------------
C---- PROGRAMM ZUM ZEICHNEN VON HINTEREINANDER ANGEORDNETEN SCHEIBEN
C---- ODER SPAEREN (KUGELN)
C----------------------------------------------------------------------
C---- AUTOR: GROTEN, KFA JUELICH, ZAM
C---- DATUM:  2.7.1987
C---- UPDATE  8.8.1990
C---- UPDATE  8.1.1992 GROTEN EPSR
C----------------------------------------------------------------------
C NS<0  : VERHINDERT SORTIEREN
C ZP=0. : KEINE ZENTRALPROJEKTION
      SUBROUTINE GRSPHR(NS,NR,VS,INFA,CHI,PSI,ZP,HS,IPO,MDM,IER)
      PARAMETER( PI=3.1415926536,MAXLIN=1500,INT2=32768,ZWEIPI=2.*PI)
      PARAMETER( PID2=PI*0.5, PID180=PI/180.)
      PARAMETER (MAXM=20,C666=3.5,IZH30=1024*1024*1024)
      INTEGER INFA(NS),IPO(*)
      REAL VS(4,*),HS(*)
      REAL SINUS(MAXM),TANG(MAXM),COSIN(MAXM),ECKE(3,2,2,2),XX(2),YY(2)
      LOGICAL INLI,INNEN1,INNEN2,INNEN,DURCH,KASTEN
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRPP/
C
C KASTEN BEI ZENTRALPROJEKTION  ?
C
      IF (ZP.LT.0.) THEN
         KASTEN=.TRUE.
         AZP=-ZP
      ELSE
         KASTEN=.FALSE.
         AZP=ZP
      ENDIF
C
C FEHLER BEI DER DIMENSIONIERUNG ?
C
      IF (NR.LT.1 .OR. NR.GT.20) THEN
         IF ( IER.NE.2 ) WRITE(*,*)
     >   'GRSPHR: NR MUST BE POSITIVE AND LESS THAN OR EQUAL TO 20'
         IF (IER.NE.1 .AND. IER.NE.2) STOP 'GRSPHR 30'
         IER=30
         RETURN
      ENDIF
      NN=ABS(NS)
      N=NN*NR
      M=4*N
      IF (MDM.LT.N) THEN
         IF ( IER.NE.2 ) WRITE(*,*)
     >   'GRSPHR: MDM IS TOO SMALL FOR ALL CASES'
         IF (IER.NE.1 .AND. IER.NE.2) STOP 'GRSPHR 10'
         IER=10
         RETURN
      ENDIF
C
C ERGAENZE DIE INNEREN KREISE FUER KUGELN, BERECHNE MITTLERE X,Y
C
      IF (CHI.NE.90.) THEN
         TGCHI=TAN(CHI*PID180)
      ELSE
C------- EPSR  : kleinste positive reelle Zahl X mit '1.+X > 1.'
         EPSR  = .000015258791
         EPSI = EPSR*.5
   2     IF ( 1.+EPSI .GT. 1. ) THEN
            EPSR = EPSI
            EPSI = EPSR*.5
            GOTO 2
         ENDIF
         TGCHI=tan(pid2*(1.-epsr))
      ENDIF
      XM=0.
      YM=0.
      IF (NN.GE.1) THEN
         XMIN= VS(1,1)-VS(4,1)
         XMAX= VS(1,1)+VS(4,1)
         YMIN= VS(2,1)-VS(4,1)
         YMAX= VS(2,1)+VS(4,1)
         ZMIN= VS(3,1)-VS(4,1)
         ZMAX= VS(3,1)+VS(4,1)
      ENDIF
      DO 11 I=1,NN
         XM=XM+VS(1,I)
         YM=YM+VS(2,I)
         XMIN=MIN(XMIN,VS(1,I)-VS(4,I))
         XMAX=MAX(XMAX,VS(1,I)+VS(4,I))
         YMIN=MIN(YMIN,VS(2,I)-VS(4,I))
         YMAX=MAX(YMAX,VS(2,I)+VS(4,I))
         ZMIN=MIN(ZMIN,VS(3,I)-VS(4,I))
         ZMAX=MAX(ZMAX,VS(3,I)+VS(4,I))
         HS(-3+4*NR*I)=VS(1,I)
         HS(-2+4*NR*I)=VS(2,I)
         HS(-1+4*NR*I)=VS(3,I)
         HS(4*NR*I)=VS(4,I)
         HS(M-1+2*NR*I)=-PI
         HS(M+2*NR*I)= PI
  11  CONTINUE
      XM=XM/NN
      YM=YM/NN
C
C SORTIERE NACH Z ( =VS(3,*) )
C
      IF (NS.GT.0) CALL GRSSRT(HS(-3+4*NR),INFA,NN,NR)

      DO 12 I=1,NR-1
         SINUS(I)=SIN(PID2/NR*I)
         COSIN(I)=COS(PID2/NR*I)
         TANG(I)=SINUS(I)/COSIN(I)
  12  CONTINUE

      ALPHA=PSI*PID180
      IPOS=N+1
      IF (IPOS.GT.MDM) GOTO 4711
      DO 14 I=1,NN
         INTENS=MOD(INFA(I),INT2)
         IWI=I*NR
         IPO(IWI)=IWI
         IPO(N+IWI)=0
         DO 13 J=1,NR-1
            IWJ=(I-1)*NR+J
            HS(-3+4*IWJ)=HS(-3+4*NR*I)
            HS(-2+4*IWJ)=HS(-2+4*NR*I)
            HS(-1+4*IWJ)=HS(-1+4*NR*I)+HS(4*NR*I)*COSIN(J)
            HS(4*IWJ)=HS(4*NR*I)*SINUS(J)
            ACO=TGCHI*TANG(J)
            IF (ACO .LT. 1. AND. ACO.GE.0.)  THEN
               HS(M-1+2*IWJ)= -PI
               HS(M+2*IWJ)=  PI
               IPO(IWJ)=IWJ
               IPO(N+IWJ)=0
            ELSE IF (ACO .GT. -1. AND. ACO.LT.0.)  THEN
               IPO(IWJ)=0
               IPO(N+IWJ)=0
               IF (INTENS.GT.16) THEN
                  HS(M-1+2*IWJ)= -PI
                  HS(M+2*IWJ)=  PI
                  IPO(IWJ)=IWJ+IZH30
                  IPO(N+IWJ)=0
               ENDIF
            ELSE
               ACO=ACOS(-1./ACO)
               HS(M-1+2*IWJ)=-ACO+ALPHA
               HS(M+2*IWJ)= ACO+ALPHA
               IPO(IWJ)=IWJ
               IPO(N+IWJ)=0
               IF (INTENS.GT.16) THEN
                  IF ( HS(M-1+2*IWJ)+ZWEIPI.GT. ZWEIPI) THEN
                     WI1=HS(M+2*IWJ)-ZWEIPI
                     WI2=HS(M-1+2*IWJ)
                  ELSE
                     WI1=HS(M+2*IWJ)
                     WI2=HS(M-1+2*IWJ)+ZWEIPI
                  ENDIF
                  IPO(N+IWJ)=IPOS+IZH30
                  HS(M-1+2*IPOS)=WI1
                  HS(M+2*IPOS)=WI2
                  IPO(N+IPOS)=0
                  IPOS=IPOS+1
                  IF (IPOS.GT.MDM) GOTO 4711
               ENDIF
            ENDIF
  13     CONTINUE
         IF (AZP.NE.0.) THEN
            IF (TGCHI.EQ.0.)  THEN
               HS(M-1+2*IWI)= -PI
               HS(M+2*IWI)=  PI
            ELSE
               HS(M-1+2*IWI)= ALPHA - PID2
               HS(M+2*IWI)= ALPHA + PID2
               IF (INTENS.GT.16) THEN
                  IF ( HS(M-1+2*IWI)+ZWEIPI .GT. ZWEIPI) THEN
                     WI1=HS(M+2*IWI)-ZWEIPI
                     WI2=HS(M-1+2*IWI)
                  ELSE
                     WI1=HS(M+2*IWI)
                     WI2=HS(M-1+2*IWI)+ZWEIPI
                  ENDIF
                  IPO(N+IWI)=IPOS+IZH30
                  HS(M-1+2*IPOS)=WI1
                  HS(M+2*IPOS)=WI2
                  IPO(N+IPOS)=0
                  IPOS=IPOS+1
                  IF (IPOS.GT.MDM) GOTO 4711
               ENDIF
            ENDIF
         ENDIF
  14  CONTINUE
C
C WENN GEWUENSCHT ZENTRALPROJEKTION
C
      IF (AZP.NE.0.) THEN
         IF (AZP.NE.666.) THEN
            ZPP=AZP
         ELSE
            IF (N.GE.1) THEN
               ZPP=C666*(HS(1)-XM)+HS(3)-ZMIN
            ENDIF
            DO 20 I=1,N
               ZPP=MAX(ZPP,C666*(HS(-3+4*I)-XM)+HS(-1+4*I)-ZMIN)
               ZPP=MAX(ZPP,C666*(HS(-2+4*I)-YM)+HS(-1+4*I)-ZMIN)
   20       CONTINUE
         ENDIF
         ZSU=ZPP+ZMIN
         DO 19 I=1,N
            ZFA=ZPP/(ZSU-HS(-1+4*I))
            HS(-3+4*I)=(HS(-3+4*I)-XM)*ZFA+XM
            HS(-2+4*I)=(HS(-2+4*I)-YM)*ZFA+YM
            HS(4*I)=HS(4*I)*ZFA
  19     CONTINUE
C
C  BESTIMMUNG DES MASSTABS
C
         XMI=(XMIN-XM)*ZPP/(ZSU-ZMAX)+XM
         XMA=(XMAX-XM)*ZPP/(ZSU-ZMAX)+XM
         YMI=(YMIN-YM)*ZPP/(ZSU-ZMAX)+YM
         YMA=(YMAX-YM)*ZPP/(ZSU-ZMAX)+YM
      ELSE
         XMI=XMIN
         XMA=XMAX
         YMI=YMIN
         YMA=YMAX
      ENDIF
      IF( (YMA-YMI)*(PP(3)-PP(1)).LT.(XMA-XMI)*(PP(4)-PP(2)) ) THEN
         YQ=(YMA+YMI)*.5
         DY=(PP(4)-PP(2))/(PP(3)-PP(1))*(XMA-XMI)*.5
         YMI=YQ-DY
         YMA=YQ+DY
      ELSE
         XQ=(XMA+XMI)*.5
         DX=(PP(3)-PP(1))/(PP(4)-PP(2))*(YMA-YMI)*.5
         XMI=XQ-DX
         XMA=XQ+DX
      ENDIF
      CALL GRSCLV(XMI,YMI,XMA,YMA)
C
C ZEICHNE DEN KASTEN
C
      IF (KASTEN) THEN
         DO 16 I=1,2
            DO 15 J=1,2
               ECKE(3,I,J,1)=ZMAX
               ECKE(3,I,J,2)=ZMIN
               ECKE(1,1,I,J)=XMAX
               ECKE(1,2,I,J)=XMIN
               ECKE(2,I,1,J)=YMAX
               ECKE(2,I,2,J)=YMIN
 15         CONTINUE
 16      CONTINUE
         CALL GRDSH(1.,0.,1.)
         DO 18 I=1,2
            DO 17 J=1,2
               XX(1)=(ECKE(1,I,J,1)-XM)*ZPP/(ZSU-ECKE(3,I,J,1))+XM
               XX(2)=(ECKE(1,I,J,2)-XM)*ZPP/(ZSU-ECKE(3,I,J,2))+XM
               YY(1)=(ECKE(2,I,J,1)-YM)*ZPP/(ZSU-ECKE(3,I,J,1))+YM
               YY(2)=(ECKE(2,I,J,2)-YM)*ZPP/(ZSU-ECKE(3,I,J,2))+YM
               CALL GPL(2,XX,YY)
               XX(1)=(ECKE(1,1,I,J)-XM)*ZPP/(ZSU-ECKE(3,1,I,J))+XM
               XX(2)=(ECKE(1,2,I,J)-XM)*ZPP/(ZSU-ECKE(3,2,I,J))+XM
               YY(1)=(ECKE(2,1,I,J)-YM)*ZPP/(ZSU-ECKE(3,1,I,J))+YM
               YY(2)=(ECKE(2,2,I,J)-YM)*ZPP/(ZSU-ECKE(3,2,I,J))+YM
               CALL GPL(2,XX,YY)
               XX(1)=(ECKE(1,I,1,J)-XM)*ZPP/(ZSU-ECKE(3,I,1,J))+XM
               XX(2)=(ECKE(1,I,2,J)-XM)*ZPP/(ZSU-ECKE(3,I,2,J))+XM
               YY(1)=(ECKE(2,I,1,J)-YM)*ZPP/(ZSU-ECKE(3,I,1,J))+YM
               YY(2)=(ECKE(2,I,2,J)-YM)*ZPP/(ZSU-ECKE(3,I,2,J))+YM
               CALL GPL(2,XX,YY)
   17    CONTINUE
 18     CONTINUE
      ENDIF
C
C UEBERDECKTE KREISSTUECKE WEGLASSEN
C
      DO 40 I=2,N
         DO 35 J=1,I-1
            IF (IPO(J).LE.0) GOTO 35
            ABS2 =(HS(-3+4*J)-HS(-3+4*I))**2+(HS(-2+4*J)-HS(-2+4*I))**2
C
C FALLUNTERSCHEIDUNG:
C
C  A) 2 GETRENNTE SCHEIBEN -> KEINE WEITEREN BERECHNUNGEN
C
            IF ((HS(4*J)+HS(4*I))**2 .LE. ABS2) GOTO 35
C
C  B) 2 SCHEIBEN, DIE TEILWEISE UEBEREINANDERLIEGEN
C
            IF ((HS(4*J)-HS(4*I))**2 .LT. ABS2) THEN
               PHI = ACOS( (HS(4*I)**2 + ABS2 - HS(4*J)**2) /
     >                   (2. * SQRT(ABS2) * HS(4*I)) )
               ALPHA=ATAN2(HS(-2+4*J)-HS(-2+4*I),HS(-3+4*J)-HS(-3+4*I))
C
C              BERECHNUNG VON ANFANGS- UND ENDWINKEL DES BOGENS
C
               WI1 = ALPHA+PHI
               WI2 = ALPHA-PHI
               IF (ALPHA.GE.0.) THEN
                  WI1=WI1-ZWEIPI
               ELSE
                  WI2=WI2+ZWEIPI
               ENDIF
               IF (IPO(I).GT.IZH30) THEN
                  JP=IPO(I)-IZH30
                  DURCH=.FALSE.
               ELSE
                  JP=IPO(I)
                  DURCH=.TRUE.
               ENDIF
               INLI=.TRUE.
               IF (JP.EQ.0) GOTO 40
               DO 30 K=1,2 000 000
                  CALL GRQARC(HS(M-1+2*JP),HS(M+2*JP),WI1,ZKP1,INNEN1)
                  CALL GRQARC(HS(M-1+2*JP),HS(M+2*JP),WI2,ZKP2,INNEN2)
                  IF (INNEN1) THEN
                     IF (INNEN2) THEN
                        IF (ZKP1.EQ.ZKP2 .OR.
     >                   (HS(M-1+2*JP).EQ.-PI.AND.HS(M+2*JP).EQ.PI))THEN
                           HS(M-1+2*JP)=WI1
                           HS(M+2*JP)=WI2
                        ELSE
                           WI1A=HS(M-1+2*JP)
                           WI2A=HS(M+2*JP)
                           HS(M-1+2*JP)=WI1+ZKP1
                           HS(M+2*JP)=WI2A
                           JPN=IPO(N+JP)
                           IF (DURCH) THEN
                              IPO(N+JP)=IPOS
                           ELSE
                              IPO(N+JP)=IPOS+IZH30
                           ENDIF
                           JP=IPOS
                           IPOS=IPOS+1
                           IF (IPOS.GT.MDM) GOTO 4711
                           IPO(N+JP)=JPN
                           HS(M-1+2*JP)=WI1A
                           HS(M+2*JP)=WI2+ZKP2
                           INLI=.FALSE.
                        ENDIF
                     ELSE
                        HS(M-1+2*JP)=WI1+ZKP1
                     ENDIF
                  ELSE
                     IF (INNEN2) THEN
                        HS(M+2*JP)=WI2+ZKP2
                     ELSE
                        CALL GRQARC(WI1,WI2,HS(M-1+2*JP),ZKP1,INNEN)
                        IF (INNEN) THEN
                        ELSE
                           IF (INLI) THEN
                              IPO(I)=IPO(N+JP)
                              IF (IPO(I).EQ.0) GOTO 40
                              GOTO 33
                           ELSE
                              IPO(N+JPO)=IPO(N+JP)
C                             JP=JPO
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  JPO=JP
                  INLI=.FALSE.
                  IF (IPO(N+JP).EQ.0) GOTO 35
 33               IF (IPO(N+JP).GT.IZH30) THEN
                     JP=IPO(N+JP)-IZH30
                     DURCH=.FALSE.
                  ELSE
                     JP=IPO(N+JP)
                     DURCH=.TRUE.
                  ENDIF
  30           CONTINUE
C
C  C) 2 SCHEIBEN, DIE GANZ UEBEREINANDERLIEGEN (UNTERE KLEINER)
C
            ELSE IF( HS(4*I).LE.HS(4*J) ) THEN
               IPO(I)=0
               GOTO 40
C
C  D) 2 SCHEIBEN, DIE GANZ UEBEREINANDERLIEGEN (UNTERE GROESSER)
C
            ELSE
               IPO(J)=-IPO(J)
            ENDIF
  35     CONTINUE
  40  CONTINUE
C
C ZEICHNEN ALLER KREISBOEGEN
C
      DO 51 I=1,N
         IF (ABS(IPO(I)).GE.IZH30) THEN
            IP=ABS(IPO(I))-IZH30
            DURCH=.FALSE.
         ELSE
            IP=ABS(IPO(I))
            DURCH=.TRUE.
         ENDIF
         IF (IP.NE.0) THEN
            INTENS=MOD(INFA((I-1)/NR+1),INT2)
            CALL GRSPTS(INTENS)
            ICOLOR=(INFA((I-1)/NR+1)-INTENS)/INT2
            CALL GRNWPN(ICOLOR)
         ENDIF
         DO 49 K=1,2 000 000
            IF (IP.NE.0) THEN
               IF (DURCH) THEN
                  CALL GRDSH(1.,0.,1.)
                  CALL GRCRCL (HS(-3+4*I),HS(-2+4*I),HS(4*I),
     >               HS(M-1+2*IP),HS(M+2*IP))
               ELSE
                  CALL GRDSH(.1,.5,.1)
                  CALL GRCRCL (HS(-3+4*I),HS(-2+4*I),HS(4*I),
     >               HS(M-1+2*IP),HS(M+2*IP))
               ENDIF
               IF (IPO(N+IP).GT.IZH30) THEN
                  IP=IPO(N+IP)-IZH30
                  DURCH=.FALSE.
               ELSE
                  IP=IPO(N+IP)
                  DURCH=.TRUE.
               ENDIF
            ELSE
               GOTO 50
            ENDIF
   49    CONTINUE
   50 CONTINUE
   51 CONTINUE
      CALL GRDSH(1.,0.,1.)
      IER=0
      RETURN
 4711 IF ( IER.NE.2 ) WRITE(*,*)
     >   'GRSPHR: MDM IS TOO SMALL FOR THIS CASE'
      IF (IER.NE.1 .AND. IER.NE.2) STOP 'GRSPHR 20'
      IER=20
      END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
      SUBROUTINE GRQARC(PH1,PH2,WI,W,OBIN)
      REAL ZP(3)
      LOGICAL OBIN
      DATA ZP/0.,-6.2831853072,+6.2831853072/
      DO 1 I=1,3
         W=ZP(I)
         IF (WI+W.GT.PH1 .AND. WI+W.LT.PH2) THEN
            OBIN=.TRUE.
            GOTO 2
         ENDIF
  1   CONTINUE
      OBIN=.FALSE.
  2   END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C
C   NAME                -   GRSSRT
C-----------------------------------------------------------------------
C
C   PURPOSE             -   SORTING OF AN ARRAY A(4,*) IN DESCENDING
C                           ORDER OF A(3,*). IFA IS CHANGED IN SAME WAY.
C
C   COMPILER            -   VS FORTRAN
C
C   USAGE               -   CALL GRSSRT(A,IFA,N)
C
C   ARGUMENTS    A      -   ARRAY (4,*) OF LENGTH N.
C                           ON INPUT, A CONTAINS THE ARRAY TO BE
C                           SORTED.
C                           ON OUTPUT, A CONTAINS THE SORTED VALUES
C                           OF THE ARRAY.                 (INPUT/OUTPUT)
C
C                IFA    -   INTEGER ARRAY TO BE CHANGED IN SAME ORDER
C
C                N      -   INTEGER VARIABLE CONTAINING THE
C                           NUMBER OF ELEMENTS IN THE ARRAY TO
C                           BE SORTED.                          (INPUT)
C
C   REQD. ROUTINES      -   NONE
C
C   REFERENCE           -   KNUTH D.E., THE ART OF COMPUTER PRO-
C                           GRAMMING VOL.3, SORTING AND SEARCHING,
C                           ADDISON-WESLEY, 1973
C
C   REMARKS             -   THERE IS NO CRAY VERSION
C
C   ALGORITHM           -   QUICKSORT IS USED
C
      SUBROUTINE GRSSRT(A,IFA,JJ,INCA)
      INTEGER   IFA(JJ),IT(20),LT(20)
      REAL   A(4,*),X,T1,T2,T3,T4
C                                  FIRST EXECUTABLE STATEMENT
      J=JJ
      JA=(J-1)*INCA+1
      I=1
      IA=1
      M=1
   10 IF (J-I.LE.1) GOTO 80
      L=(J+I)/2
      LA=(L-1)*INCA+1
      T1=A(1,LA)
      T2=A(2,LA)
      T3=A(3,LA)
      T4=A(4,LA)
      IFT=IFA(L)
      A(1,LA)=A(1,IA)
      A(2,LA)=A(2,IA)
      A(3,LA)=A(3,IA)
      A(4,LA)=A(4,IA)
      IFA(L)=IFA(I)
      N=J
      NA=(N-1)*INCA+1
      K=I+1
      KA=(K-1)*INCA+1
   20 IF (K.GT.N) GOTO 50
      IF (A(3,KA).GE.T3) GOTO 40
      NE=N+1
      NEA=(NE-1)*INCA+1
      DO 30 NN=K,N
         NE=NE-1
         NEA=NEA-INCA
         IF (A(3,NEA).GT.T3) GOTO 45
   30 CONTINUE
      N=K-1
      NA=(N-1)*INCA+1
      GOTO 50
   45 X=A(1,KA)
      A(1,KA)=A(1,NEA)
      A(1,NEA)=X
      X=A(2,KA)
      A(2,KA)=A(2,NEA)
      A(2,NEA)=X
      X=A(3,KA)
      A(3,KA)=A(3,NEA)
      A(3,NEA)=X
      X=A(4,KA)
      A(4,KA)=A(4,NEA)
      A(4,NEA)=X
      IX=IFA(K)
      IFA(K)=IFA(NE)
      IFA(NE)=IX
      N=NE-1
      NA=(N-1)*INCA+1
   40 K=K+1
      KA=(K-1)*INCA+1
      GOTO 20
   50 A(1,IA)=A(1,NA)
      A(2,IA)=A(2,NA)
      A(3,IA)=A(3,NA)
      A(4,IA)=A(4,NA)
      IFA(I)=IFA(N)
      A(1,NA)=T1
      A(2,NA)=T2
      A(3,NA)=T3
      A(4,NA)=T4
      IFA(N)=IFT
      IF (N+N.LE.I+J) GOTO 60
      LT(M)=I
      IT(M)=N-1
      I=N+1
      IA=(I-1)*INCA+1
      GOTO 70
   60 LT(M)=N+1
      IT(M)=J
      J=N-1
      JA=(J-1)*INCA+1
   70 M=M+1
      GOTO 10
   80 IF (I.GE.J) GOTO 90
      IF (A(3,IA).GE.A(3,JA)) GOTO 90
      X=A(1,IA)
      A(1,IA)=A(1,JA)
      A(1,JA)=X
      X=A(2,IA)
      A(2,IA)=A(2,JA)
      A(2,JA)=X
      X=A(3,IA)
      A(3,IA)=A(3,JA)
      A(3,JA)=X
      X=A(4,IA)
      A(4,IA)=A(4,JA)
      A(4,JA)=X
      IX=IFA(I)
      IFA(I)=IFA(J)
      IFA(J)=IX
   90 M=M-1
      IF (M.LE.0) RETURN
      I=LT(M)
      IA=(I-1)*INCA+1
      J=IT(M)
      JA=(J-1)*INCA+1
      GOTO 10
      END
