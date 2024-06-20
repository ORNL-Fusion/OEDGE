C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 4.8.92  Groten: Dimension von XX,YY auf 7 erhoeht.  SUBCHK
C UPDATE 3.2.93  Groten: an Ausblendung korrigiert; II-Schleife
C UPDATE 2.12.89 Groten: Argumente X und Y unused, entfernt.

      SUBROUTINE GR2ISL(NX,IX,NY,IY,WE,N1,F,Q)
c     Wird von GRHHFL und von CRMENU aufgerufen; nicht fuer Benutzer!

c     .. parameters ..
      REAL              BLIND1,BLIND2
      PARAMETER         (BLIND1=-75.7501E20,BLIND2=-75.7499E20)
c     ..
c     .. scalar arguments ..
      INTEGER           IX,IY,N1,NX,NY
c     ..
c     .. array arguments ..
      REAL              F(N1, (NY-1)*IY+1),Q(2,NX,NY), WE(2)
c     ..
c     .. local scalars ..
      REAL              AF,C11,C12,C21,C22,F1,F2,FA
      INTEGER           I,IP,J,JN,KK,L,M
c     ..
c     .. local arrays ..
      REAL              W(2,0:11),XX(7),YY(7)
      INTEGER           IPX(2,0:11),IPY(2,0:11),LI1(3,0:11),LI2(3,0:11),
     $                  IPZ(2,0:11)
      LOGICAL           NOTYET(0:11),OB0(0:11)
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. data statements ..

      DATA              IPX/0,1,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,
     $                  0,0/
      DATA              IPY/0,0,0,1,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,
     $                  0,1/
      DATA              IPZ/0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,1,
     $                  1,1/
      DATA              LI1/3,2,1,0,3,2,7,10,6,4,11,7,0,5,8,1,6,9,2,7,
     $                  10,3,4,11,4,0,5,5,1,6,11,8,9,8,9,10/
      DATA              LI2/5,8,4,6,9,5,1,0,3,2,1,0,11,7,3,8,4,0,9,5,1,
     $                  10,6,2,9,10,11,10,11,8,6,2,7,7,3,4/
c     ..
c     .. executable statements ..

         DO 100 J = 1,NY - 1
            DO 90 I = 1,NX - 1
               DO 50 L = 0,11
                  NOTYET(L) = .TRUE.
                  F1     = F((I+IPX(1,L)-1)*IX+1, (J+IPY(1,L)-1)*IY+1)
                  F2     = F((I+IPX(2,L)-1)*IX+1, (J+IPY(2,L)-1)*IY+1)
                  IF ((F1.LE.BLIND1.OR.F1.GE.BLIND2) .AND.
     $                (F2.LE.BLIND1.OR.F2.GE.BLIND2)) THEN
                     W(1,L) = F1 - WE(1+IPZ(1,L))
                     W(2,L) = F2 - WE(1+IPZ(2,L))
                     OB0(L) = (W(1,L).LT.0. .AND. W(2,L).GE.0.) .OR.
     $                        (W(2,L).LT.0. .AND. W(1,L).GE.0.)
*
                  ELSE
                     DO 49 II=0,11
                        OB0(II) = .FALSE.
   49                CONTINUE
                     GOTO 51
                  END IF
*
   50          CONTINUE
   51          DO 80 L = 0,9
                  IF (OB0(L) .AND. NOTYET(L)) THEN
                     M      = L
                     IP     = 0
   60                CONTINUE
                     IF (NOTYET(M)) THEN
                        IP     = IP + 1
                        FA     = -W(1,M)/ (W(2,M)-W(1,M))
                        AF     = 1. - FA
                        C11    = Q(1,I+IPX(1,M),J+IPY(1,M))
                        C12    = Q(1,I+IPX(2,M),J+IPY(2,M))
                        C21    = Q(2,I+IPX(1,M),J+IPY(1,M))
                        C22    = Q(2,I+IPX(2,M),J+IPY(2,M))
                        XX(IP) = AF*C11 + FA*C12
                        YY(IP) = AF*C21 + FA*C22
                        NOTYET(M) = .FALSE.
                        DO 70 KK = 1,3

                           IF (W(1,M).LT.0.) THEN
                              JN     = LI1(KK,M)
*
                           ELSE
                              JN     = LI2(KK,M)
                           END IF

                           IF (OB0(JN)) THEN
                              M      = JN
                              GO TO 60
*
                           END IF
*
   70                   CONTINUE
                     END IF

c                    ein stueck zu ende

                     IF (IP.GT.2) CALL GRFILL(IP,XX,YY,1,0)
                  END IF

   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
      END
