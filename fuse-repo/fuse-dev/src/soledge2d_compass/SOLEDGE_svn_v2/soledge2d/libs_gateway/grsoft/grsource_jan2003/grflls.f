C@PROCESS NOSDUMP NOGOSTMT OPT(3) IL(DIM)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C---- WRITTEN 18.11.91 GROTEN
C UPDATE 20. 2.1992 Busch , GKS Treiber erlaubt max. 300 Werte pro Kurve
C                   GPL Aufruf durch GRLN Aufruf ersetzt
C UPDATE 12. 3.1992 Busch , GRDSH darf nicht auf KFA Marker wirken,
C                   pp(10-12) muss 1.,0.,1. sein

      SUBROUTINE GRFLLS(N,X,Y,NS)
C---- IMPLICIT NONE
      INTEGER N,NS
      REAL X(N),Y(N)

      REAL PP(8),DSH1,DSH2,DSH3,ZEIGRO,ZEIWI,SYGRO
      INTEGER IFONT,INTLIN,ICOL
      LOGICAL FLGROT
      COMMON /GRPP/ PP,IFONT,DSH1,DSH2,DSH3,INTLIN,ZEIGRO,ZEIWI,ICOL,
     $              SYGRO,FLGROT
CDEC$ PSECT /GRPP/ NOSHR
      REAL XMAXDC,XUNITS,YUNITS
      COMMON /SCALE/ XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR
      INTEGER NSCLC,NSCLV,NSCLP
      REAL FLPIC,XMAXCM,YMAXCM,XDCPIC,YDCPIC,RAHMEN
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR

      SAVE  /GRPP/, /SCALE/ ,/GRPIC/
      INTEGER IS,KIND,I,J,IERR,NUMTR
      REAL    XX(25),YY(25),SCLV(4),SCLC(4),U,V

      REAL X0(3),Y0(3)
      REAL X1(4),Y1(4)
      REAL X2(6),Y2(6)
      REAL X3(24),Y3(24)
      REAL X4(12),Y4(12)
      REAL X5(9),Y5(9)
      REAL X6(12),Y6(12)
      REAL X7(12),Y7(12)
      REAL X8(10),Y8(10)
      REAL X9(5),Y9(5)

      DATA X0/.4330127,0.,-.4330127/
      DATA Y0/    -.25,.5,     -.25/

      DATA X1/.3535534,-.3535534,-.3535534, .3535534/
      DATA Y1/.3535534, .3535534,-.3535534,-.3535534/
      DATA X2/.5,.25     ,-.25    ,-.5,-.25     ,.25       /
      DATA Y2/.0,.4330127,.4330127, .0,-.4330127, -.4330127/

      DATA X3/
     $   -0.500000,-0.482963,-0.433012,-0.353553,-0.249999,-0.129409,
     $    0.000001, 0.129410, 0.250000, 0.353554, 0.433013, 0.482963,
     $    0.500000, 0.482963, 0.433013, 0.353554, 0.250000, 0.129410,
     $    0.000001,-0.129409,-0.249999,-0.353553,-0.433012,-0.482963/
      DATA Y3/
     $    0.000000,-0.129411,-0.250001,-0.353554,-0.433013,-0.482963,
     $   -0.500000,-0.482963,-0.433013,-0.353553,-0.250000,-0.129410,
     $    0.000000, 0.129410, 0.250000, 0.353553, 0.433013, 0.482963,
     $    0.500000, 0.482963, 0.433013, 0.353554, 0.250001, 0.129411/

      DATA X4/.47,.16,.16,-.16,-.16,-.47,-.47,-.16,-.16,.16, .16,.47 /
      DATA Y4/.16,.16,.47, .47, .16,.16,-.16,-.16,-.47,-.47,-.16,-.16/

      DATA X5/
     $     0.080000, 0.080000, 0.473013, 0.393013, 0.000000,-0.393013,
     $    -0.473013,-0.080000,-0.080000/
      DATA Y5/
     $    -0.500000,-0.069282, 0.180718, 0.319282, 0.092376, 0.319282,
     $     0.180718,-0.046188,-0.500000/

      DATA X6/.5,.25,.1,-.1,-.25,-.5,-.2,-.1,-.25,.25,.1,.2/
      DATA Y6/0,.433013,.173205,.173205,.433013,0,0,-.173205,
     $         -.433013,-.433013,-.173205,0/

      DATA X7/ .218508,.353553,.135045,-.135045,-.353553,-.218508,
     $       -.218508,-.353553,-.135045, .135045, .353553, .218508/
      DATA Y7/.135045,.353553,.218508,.218508,.353553,.135045,
     $       -.135045,-.353553,-.218508,-.218508,-.353553,-0.135045/

      DATA X8/0,.112257,.475528,.181636,.293893,0,
     $        -.293893,-.181636,-.475528,-.112257/
      DATA Y8/.5,.154508,.154508,-.059017,-.404508,-.190983,
     $        -.404508,-.059017,.154508,.154508/

      DATA X9/ 0,.293893,-.475528,.475528,-.293893/
      DATA Y9/ .5,-.404508,.154508,.154508,-.404508/


c     Default : durchgezogene Linie
      dsh1s=dsh1
      dsh2s=dsh2
      dsh3s=dsh3
      call grdsh(1.,0.,1.)

C---- MARKERAUSGABE NUR IM GEAENDERTEN WINDOW MOEGLICH
C---- WINDOW WIRD AUF CM ZURUECKGESETZT
      CALL GQCNTN(IERR,NUMTR)
      CALL GQNT(NUMTR,IERR,SCLV,SCLC)
C---- NEUES WINDOW IN CM , NUMTR=2 Clipping
      IF (NUMTR.EQ.1) THEN
         CALL GSWN(NUMTR,0.,XMAXCM,0.,YMAXCM)
      ELSE
         CALL GSWN(NUMTR,PP(1),PP(3),PP(2),PP(4))
      ENDIF

      IS = MOD(NS,100)+1
      KIND = MOD((NS-100)/100,2)
     

      select case(is)

      case(1)

      DO 12 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 11 I=1,3
            XX(I) = U+X0(I)*SYGRO
            YY(I) = V+Y0(I)*SYGRO
   11    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(3,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL (I,XX,YY)
c           CALL GRLN(XX,YY,I)
         ENDIF
   12 CONTINUE
     

      case(2)

      DO 22 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 21 I=1,4
            XX(I) = U+X1(I)*SYGRO
            YY(I) = V+Y1(I)*SYGRO
   21    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(4,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   22 CONTINUE
     

      case(3)

      DO 32 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 31 I=1,6
            XX(I) = U+X2(I)*SYGRO
            YY(I) = V+Y2(I)*SYGRO
   31    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(6,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   32 CONTINUE
     

      case(4)

      DO 42 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 41 I=1,24
            XX(I) = U+X3(I)*SYGRO
            YY(I) = V+Y3(I)*SYGRO
   41    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(24,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   42 CONTINUE
   

      case(5)

      DO 52 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 51 I=1,12
            XX(I) = U+X4(I)*SYGRO
            YY(I) = V+Y4(I)*SYGRO
   51    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(12,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   52 CONTINUE
      

      case(6)

      DO 62 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 61 I=1,9
            XX(I) = U+X5(I)*SYGRO
            YY(I) = V+Y5(I)*SYGRO
   61    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(9,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   62 CONTINUE
      

      case(7)

      DO 72 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 71 I=1,12
            XX(I) = U+X6(I)*SYGRO
            YY(I) = V+Y6(I)*SYGRO
   71    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(12,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   72 CONTINUE
      

      case(8)

      DO 82 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 81 I=1,12
            XX(I) = U+X7(I)*SYGRO
            YY(I) = V+Y7(I)*SYGRO
   81    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(12,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   82 CONTINUE
      

      case(9)

      DO 92 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 91 I=1,10
            XX(I) = U+X8(I)*SYGRO
            YY(I) = V+Y8(I)*SYGRO
   91    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(10,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
   92 CONTINUE
      

      case(10)

      DO 102 J=1,N
         U = PP(1)+(X(J)-PP(5))/XUNITS
         V = PP(2)+(Y(J)-PP(6))/YUNITS
         DO 101 I=1,5
            XX(I) = U+X9(I)*SYGRO
            YY(I) = V+Y9(I)*SYGRO
  101    CONTINUE
         IF (KIND.EQ.1) THEN
            CALL GRFILL(5,XX,YY,KIND,0)
         ELSE
            XX(I) = XX(1)
            YY(I) = YY(1)
            CALL GPL(I,XX,YY)
         ENDIF
  102 CONTINUE 

      end select

C---- DAS ORIGINALWINDOW MUSS WIEDER ETABLIERT WERDEN
      CALL GSWN(NUMTR,SCLV(1),SCLV(2),SCLV(3),SCLV(4))

c     Dashline wieder restaurieren
      call GRDSH(DSH1S,DSH2S,DSH3S)
      dsh1=dsh1s
      dsh2=dsh2s
      dsh3=dsh3s

      END
