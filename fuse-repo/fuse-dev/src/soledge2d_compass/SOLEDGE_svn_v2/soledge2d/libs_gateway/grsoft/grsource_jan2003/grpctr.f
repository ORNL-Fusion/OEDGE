C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 16. 8.1990
      SUBROUTINE GRPCTR(N,LP,Z)
      INTEGER N,LP
      CHARACTER (len=8) Z

      CHARACTER (len=25) C
      WRITE(C,'(I25)') N
      DO 1 I=25-7,25
         IF (C(I:I).NE.' ') GOTO 2
    1 CONTINUE
    2 Z=C(I:)
      LP=26-I
      END
