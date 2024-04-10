C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRFTOC(F,C,L)
CCC LETZTE AENDERUNG 1. 2.85 GROTEN
CCC UPDATE 25. 9.1990  GROTEN
CCC UPDATE 18. 9.1991  GROTEN wegen RS 6000
      CHARACTER (len=6) FO
      CHARACTER (len=10) C,H
      CHARACTER (len=12) D
      INTEGER KK(8)
      DATA KK/9,9,9,9,9,8,8,7/
C
      IF(ABS(F).GT.999999E0.OR.(ABS(F).LT. .01 .AND. F.NE.0.)) THEN
         WRITE(C,'(1P,E10.3)') F
         IF (C(7:8).EQ.'E+') THEN
            H(1:2)=C(9:10)
            C(8:9)=C(9:10)
            L=9
            M=3
         ELSE
            M=4
            L=10
         ENDIF
   10    IF (C(L-M:L-M).EQ.'0') THEN
            H(1:M)=C(L-M+1:L)
            C(L-M:L)=H(1:M)
            L=L-1
            GOTO 10
         ENDIF
         IF (F.GE.0.) THEN
            H=C(2:L)
            C=H
            L=L-1
         ENDIF
         IF (C(L-1:L-1).EQ.'0') THEN
            H(1:1)=C(L:L)
            C(L-1:L-1)=H(1:1)
            L=L-1
         ENDIF
      ELSE IF (F.EQ.AINT(F)) THEN
         WRITE(C,'(I7)') INT(F)
         DO 1 I=1,7
            IF(C(I:I).NE.' ') GOTO 2
    1       CONTINUE
    2    H=C(I:7)
         C=H
         L=7-I+1
      ELSE
         WRITE(D,'(F12.4)') F
         DO 3 I=1,8
            IF(D(I:I).NE.' ') GOTO 4
    3       CONTINUE
    4    L=KK(I)
         K=MIN(I-1,4)
         IF( F.GE.0. .AND. I.LE.5 ) THEN
            L=L-1
            K=K-1
         ENDIF
         WRITE(FO,'(''(F'',I1,''.'',I1,'')'')') L,K
         WRITE(C,FO) F
   20    IF (C(L:L).EQ.'0') THEN
            L=L-1
            GOTO 20
         ENDIF
   30    IF( C(1:1).EQ.' ' .OR. C(1:1).EQ.'0' ) THEN
            H=C(2:L)
            C=H
            L=L-1
            GOTO 30
         ENDIF
         IF ( C(1:2).EQ.'-0') THEN
            H=C(3:)
            C(2:)=H
            L=L-1
         ENDIF
      ENDIF
      END
