C
C
      SUBROUTINE GRID_1(C,NT,N2,N3,C1,C2,C3,CT,IL)
C
C  MAKE A 1D GRID C(I),I=1,NT
C  CONSISTS OF UP TO 3 SECTIONS OF EQUIDISTANT GRID POINTS C
C  C(1) =C1
C  C(N2)=C2
C  C(N3)=C3
C  C(NT)=CT
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: C1, C2, C3, CT
      REAL(DP), INTENT(OUT) :: C(*)
      INTEGER, INTENT(IN) :: NT, IL
      INTEGER, INTENT(INOUT) :: N2, N3
      INTEGER :: ND1, ND2, J

      IF (C2.LT.C1) GOTO 991
      IF (C3.LT.C2) GOTO 992

C  IS THERE A 3RD SECTION ?

      IF (CT.GT.C3.AND.N3.GT.1.AND.N3.LT.NT) THEN
C YES. KNOTS IN THIS SECTION:
        ND2=NT-N3+1
      ELSE
C NO.
        ND2=1
        CT=C3
        N3=NT
      ENDIF
C
C  IS THERE A 2ND SECTION
C
      IF (C3.GT.C2.AND.N2.GT.1.AND.N2.LT.N3) THEN
C YES. KNOTS IN THIS SECTION:
        ND1=N3-N2+1
      ELSE
C NO.
        ND1=1
        C3=C2
        N2=N3
      ENDIF

      DO 101 J=1,N2
        C(J)=C1+DBLE(J-1)/DBLE(N2-1)*(C2-C1)
101   CONTINUE
      DO 102 J=N2+1,N3
        C(J)=C2+DBLE(J-N2)/DBLE(N3-N2)*(C3-C2)
102   CONTINUE
      DO 103 J=N3+1,NT
        C(J)=C3+DBLE(J-N3)/DBLE(NT-N3)*(CT-C3)
103   CONTINUE
C
      RETURN
C
991   CONTINUE
      IF (IL.EQ.1) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 1ST GRID.  RGA > RIA ?'
        WRITE (iunout,*) 'RI1,RGA = ',C1,C2
      ELSEIF (IL.EQ.2) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 2ND GRID.  YGA > YIA ?'
        WRITE (iunout,*) 'YIA,YGA = ',C1,C2
      ELSEIF (IL.EQ.3) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 3RD GRID.  ZGA > ZIA ?'
        WRITE (iunout,*) 'ZIA,ZGA = ',C1,C2
      ELSE
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: IN GRID_1   C2 > C1  ?'
        WRITE (iunout,*) 'C1,C2 = ',C1,C2
      ENDIF
      CALL EXIT_OWN(1)
992   CONTINUE
      IF (IL.EQ.1) THEN
        WRITE (iunout,*) 
     .   'GRID DATA INCONSISTENCY: 1ST GRID.  RAA > RGA ?'
        WRITE (iunout,*) 'RGA,RAA = ',C2,C3
      ELSEIF (IL.EQ.2) THEN
        WRITE (iunout,*) 
     .   'GRID DATA INCONSISTENCY: 2ND GRID.  YAA > YGA ?'
        WRITE (iunout,*) 'YGA,YAA = ',C2,C3
      ELSEIF (IL.EQ.3) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 3RD GRID.  ZAA > ZGA ?'
        WRITE (iunout,*) 'ZGA,ZAA = ',C2,C3
      ELSE
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: IN GRID_1   C3 > C2  ?'
        WRITE (iunout,*) 'C2,C3 = ',C2,C3
      ENDIF
      CALL EXIT_OWN(1)
      RETURN
      END
