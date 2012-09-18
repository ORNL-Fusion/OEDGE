c     Krieger IPP 2012 - fixed INTEL name conflict: ranf->ranf_eirene
      function ranf_eirene ()
      implicit double precision (a-h,o-z)
      ra=h1rn(dummy)
      ranf_eirene=ra
      return
      end
c
      subroutine ranset (iseed)
      implicit double precision (a-h,o-z)
      iseed1=iseed
      iseed2=2000000-iseed
      call h1rnin(iseed1,iseed2)
      return
      end
c
      subroutine ranget(iseed)
      implicit double precision (a-h,o-z)
      data ifirst/0/
      save large
      if (ifirst.eq.0) then
        large=ilarg(2)
        ifirst=1
      endif
      ran=ranf_eirene()
      iseed=ran*large
      return
      end
C
      FUNCTION ILLZ(N,V,IV)
      INTEGER N,IV,J,ILLZ
CODER REAL*8  V(*)
CODER REAL*4  V(*)
CODER INTEGER V(*)
      LOGICAL V(*)
      J=1
      IF (IV.LT.0) J=1-(N-1)*IV
      DO 17,ILLZ=1,N
CODER    IF (V(J).LT.0) GOTO 18    BEI REAL ODER INTEGER
CPB      IF (.NOT.V(J)) GOTO 18
         IF (V(J)) GOTO 18
   17    J=J+IV
   18 ILLZ=ILLZ-1
      END

c
      double precision FUNCTION SECOND_OWN()
      real*8 x05baf,start,reset_second
      save start
      data start /0.0/
c x05baf is a NAG routine
      SECOND_OWN=x05baf()-start
      RETURN
C
      ENTRY RESET_SECOND()
      start=x05baf()
      reset_second=start
      return
      END


      SUBROUTINE DEVCSF(N,A,LDA,EVAL,EVEC,LDEVEC)
      INTEGER N,LDA,LDEVEC,IERR,J,I,IAA
      PARAMETER (IAA=50)
      DOUBLE PRECISION A(LDA,N),EVAL(N),EVAL1(IAA),EVEC(LDEVEC,N),
     .       FV1(IAA),FV2(IAA)
      DOUBLE PRECISION EVEC1(IAA,IAA),TVEC(IAA*IAA)
      EQUIVALENCE (EVEC1(1,1),TVEC(1))

      entry evcsf (n,a,lda,eval,evec,ldevec)

      IF (N .GT. IAA) THEN
         WRITE(6,*) 'MATRIX ZU KLEIN IN EVCSF'
         CALL EXIT
      ENDIF
      CALL RS(LDA,N,A,EVAL1,1,EVEC1,FV1,FV2,IERR)
      DO 10,I=1,N
         EVAL(I) = EVAL1(N-I+1)
         DO 20,J=1,LDEVEC
            EVEC(J,I) = TVEC(((N-I+1)-1)*LDA+J)
20       CONTINUE
10    CONTINUE
      END


C   NAME                -   RS
C-----------------------------------------------------------------------
C
C   PURPOSE             -   CALCULATION OF EIGENVALUES AND OPTIONALLY
C                           EIGENVECTORS OF A REAL SYMMETRIC MATRIX
C
C   COMPILER            -   VS FORTRAN
C
C   PRECISION           -   IBM DOUBLE
C
C   USAGE               -   CALL RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
C   ARGUMENTS    NM     -   MUST BE SET TO THE ROW DIMENSION OF THE
C                           TWO-DIMENSIONAL ARRAY PARAMETERS AS DECLARED
C                           IN THE CALLING PROGRAM DIMENSION STATEMENT.
C                                                                (INPUT)
C
C                N      -   ORDER OF THE MATRIX  A.              (INPUT)
C
C                A      -   MATRIX OF DIMENSION N BY N,
C                           CONTAINS THE REAL SYMMETRIC MATRIX   (INPUT)
C
C                W      -   VECTOR OF LENGTH N,
C                           CONTAINS THE EIGENVALUES IN ASCENDING ORDER
C                                                               (OUTPUT)
C
C                MATZ   -   INTEGER VARIABLE SET EQUAL TO ZERO IF
C                           ONLY EIGENVALUES ARE DESIRED,  OTHERWISE
C                           SET TO ANY NON-ZERO INTEGER FOR BOTH
C                           EIGENVALUES AND EIGENVECTORS.        (INPUT)
C
C                Z      -   MATRIX OF DIMENSION N BY N,
C                           CONTAINS THE EIGENVECTORS IF MATZ IS
C                           NOT ZERO.                           (OUTPUT)
C
C                FV1    -   VECTOR OF LENGTH N. WORK AREA.
C
C                FV2    -   VECTOR OF LENGTH N. WORK AREA.
C
C                IERR   -   INTEGER OUTPUT VARIABLE SET EQUAL TO
C                           AN ERROR COMPLETION CODE DESCRIBED IN
C                           SECTION 2.3 OF THE REFERENCE. THE
C                           NORMAL COMPLETION CODE IS ZERO.     (OUTPUT)
C
C   REQD. ROUTINES      -   TRED1, TQLRAT AND
C                           TRED2, TQL2 ARE PROVIDED.
C
C   REFERENCES          -   B.T. SMITH ET.AL.
C                           MATRIX EIGENSYSTEM ROUTINES - EISPACK GUIDE
C                           SECOND EDITION 1976
C                           LECTURE NOTES IN COMPUTER SCIENCES 6
C
C   REMARKS             -   THIS SUBROUTINE CALLS THE RECOMMENDED
C                           SEQUENCE OF SUBROUTINES FROM THE EIGENSYSTEM
C                           SUBROUTINE PACKAGE (EISPACK) TO FIND THE
C                           EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C                           OF A REAL SYMMETRIC MATRIX.
C
C                           THE CORRESPONDING CRAY-VERSION IS A MEMBER
C                           OF THE $SCILIB.
C
C   ALGORITHM           -   THE MATRIX A IS REDUCED TO A SYMMETRIC
C                           TRIDIAGONAL MATRIX BY ORTHOGONAL TRANSFOR-
C                           MATIONS. THE EIGENVALUE PROBLEM FOR THE RE-
C                           DUCED MATRIX IS SOLVED BY THE QL ALGORITHM.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ------------------------------------------------------------------
C
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
CCCCC
CCCCC     AENDERUNG:    REAL ---> REAL*8
CCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER N,NM,IERR,MATZ
C     REAL A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT-
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT,
C
C        N  IS THE ORDER OF THE MATRIX  A,
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX,
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT-
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     MODIFIED FOR CFT BY MIKE ESS, CRI, JULY 1980
C
C     ------------------------------------------------------------------
C
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
C     ********** LAST CARD OF RS **********
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TQLRAT(N,D,E2,IERR)
CCCCC
CCCCC     AENDERUNG:    REAL ---> REAL*8
CCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
C     REAL D(N),E2(N)
      REAL*8 D(N),E2(N)
C     REAL B,C,F,G,H,P,R,S,MACHEP
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
C     REAL SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     MODIFIED FOR CFT BY MIKE ESS, CRI, JULY 1980
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
C
C     MACHEP IS THE SMALLEST MACHINE REPRESENTABLE REAL NUMBER SUCH
C      THAT 1.+MACHEP .NE. 1.
C
      MACHEP = SMACH( 1 )
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0
      B = 0.0
      E2(N) = 0.0
C
      DO 290 L = 1, N
         J = 0
C        H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
C        S = SQRT(E2(L))
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0 * S)
C        R = SQRT(P*P+1.0)
         R = DSQRT(P*P+1.0)
C        D(L) = S / (P + SIGN(R,P))
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** RATIONAL QL TRANSFORMATION **********
         G = D(M)
         IF (G .EQ. 0.0) G = B
         H = G
         S = 0.0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H .EQ. 0.0) GO TO 210
C        IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0) GO TO 130
  210    P = D(L) + F
C     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQLRAT **********
      END
C
      FUNCTION SMACH(I)
      REAL*8  SMACH,X02AJF,X02AKF,X02ALF
      IF(I.EQ.1) THEN
        SMACH=X02AJF()
      ELSEIF(I.EQ.2) THEN
        SMACH=X02AKF()
      ELSEIF(I.EQ.3) THEN
        SMACH=X02ALF()
      ELSE
        WRITE(6,'(A)')
     @   ' SMACH CALLED WITH PARAMETER JOB NOT EQUAL TO 1,2 OR 3.'
      ENDIF
      RETURN
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
c slmod begin - not tr
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,MACHEP,SMACH
      DOUBLE PRECISION DSQRT,DABS,DSIGN
c      REAL*8 D(N),E(N),Z(NM,N)
c      REAL*8 B,C,F,G,H,P,R,S,MACHEP,SMACH
c      REAL*8 DSQRT,DABS,DSIGN
c slmod end
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT:
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1;
C
C        E HAS BEEN DESTROYED;
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
      MACHEP = SMACH ( 1 )
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
c slmod begin - f90
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
c
c         H = MACHEP * (DABS(D(L)) + DABS(E(L)))
c slmod end
         IF (B .LT. H) B = H
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
         DO 110 M = L, N
c slmod begin - f90
            IF (ABS(E(M)) .LE. B) GO TO 120
c
c            IF (DABS(E(M)) .LE. B) GO TO 120
c slmod end
C     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
         P = (D(L+1) - D(L)) / (2.0D0 * E(L))
c slmod begin - f90
         R = SQRT(P*P+1.0D0)
         H = D(L) - E(L) / (P + SIGN(R,P))
c
c         R = DSQRT(P*P+1.0D0)
c         H = D(L) - E(L) / (P + DSIGN(R,P))
c slmod end
C
         DO 140 I = L, N
  140    D(I) = D(I) - H
C
         F = F + H
C     :::::::::: QL TRANSFORMATION ::::::::::
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         MML = M - L
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
c slmod begin - f90
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0D0)
c
c            IF (DABS(P) .LT. DABS(E(I))) GO TO 150
c            C = E(I) / P
c            R = DSQRT(C*C+1.0D0)
c slmod end
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
c slmod begin - f90
            R = SQRT(C*C+1.0D0)
c
c            R = DSQRT(C*C+1.0D0)
c slmod end
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     :::::::::: FORM VECTOR ::::::::::
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
c slmod begin - f90
         IF (ABS(E(L)) .GT. B) GO TO 130
c
c         IF (DABS(E(L)) .GT. B) GO TO 130
c slmod end
  220    D(L) = D(L) + F
  240 CONTINUE
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
C     :::::::::: LAST CARD OF TQL2 ::::::::::
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
c slmod begin - f90
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)
      DOUBLE PRECISION F,G,H,SCALE
      DOUBLE PRECISION DSQRT,DABS,DSIGN
c
c      REAL*8 A(NM,N),D(N),E(N),E2(N)
c      REAL*8 F,G,H,SCALE
c      REAL*8 DSQRT,DABS,DSIGN
c slmod end
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT:
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
  100 D(I) = A(I,I)
C     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- ::::::::::
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 120 K = 1, L
c slmod begin - f90
  120    SCALE = SCALE + ABS(A(I,K))
c
c  120    SCALE = SCALE + DABS(A(I,K))
c slmod end
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
c slmod begin - f90
         G = -SIGN(SQRT(H),F)
c
c         G = -DSIGN(DSQRT(H),F)
c slmod end
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0D0
C
         DO 240 J = 1, L
            G = 0.0D0
C     :::::::::: FORM ELEMENT OF A*U ::::::::::
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
C     :::::::::: FORM ELEMENT OF P ::::::::::
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
C
         H = F / (H + H)
C     :::::::::: FORM REDUCED A ::::::::::
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
C
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
C
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
C
      RETURN
C     :::::::::: LAST CARD OF TRED1 ::::::::::
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C
      INTEGER I,J,K,L,N,II,NM,JP1
c slmod begin - f90
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
      DOUBLE PRECISION DSQRT,DABS,DSIGN
c
c      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
c      REAL*8 F,G,H,HH,SCALE
c      REAL*8 DSQRT,DABS,DSIGN
c slmod end
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT:
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C          PRODUCED IN THE REDUCTION;
C
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 120 K = 1, L
c slmod begin- f90
  120    SCALE = SCALE + ABS(Z(I,K))
c
c  120    SCALE = SCALE + DABS(Z(I,K))
c slmod end
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
c slmod begin - f90
         G = -SIGN(SQRT(H),F)
c
c         G = -DSIGN(DSQRT(H),F)
c slmod end
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0D0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.0D0
C     :::::::::: FORM ELEMENT OF A*U ::::::::::
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     :::::::::: FORM ELEMENT OF P ::::::::::
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     :::::::::: FORM REDUCED A ::::::::::
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
         DO 280 K = 1, L
  280    Z(I,K) = SCALE * Z(I,K)
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0D0
      E(1) = 0.0D0
C     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0D0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
C     :::::::::: LAST CARD OF TRED2 ::::::::::
      END





      SUBROUTINE DLSBRR(NRA,NCA,A,LDA,B,TOL,X,RES,KBASIS)
      INTEGER NRA,NCA,LDA,KBASIS,I,J,IAA
      PARAMETER (IAA=50)
      REAL*8 A(LDA,NCA),B(NRA),X(NCA),RES(NRA),TOL,Q(IAA+2,IAA+2),R(IAA)
      INTEGER S(IAA)

      IF ((NRA.GT.IAA).OR.(NCA.GT.IAA)) THEN
         WRITE (6,*) 'MATRIX IN LSBRR ZU KLEIN'
         CALL EXIT
      ENDIF
      DO 10,I=1,NRA
         DO 20,J=1,NCA
            Q(I,J) = A(I,J)
20       CONTINUE
         Q(I,NCA+1) = B(I)
10    CONTINUE
      CALL MA20A(Q,RES,X,R,S,IAA+2,NRA,NCA,TOL)
      KBASIS = Q(NRA+1,NCA+2)
      END

      SUBROUTINE LSBRR(NRA,NCA,A,LDA,B,TOL,X,RES,KBASIS)
      INTEGER NRA,NCA,LDA,KBASIS,I,J,IAA
      PARAMETER (IAA=50)
      REAL*8 A(LDA,NCA),B(NRA),X(NCA),RES(NRA),TOL,Q(IAA+2,IAA+2),R(IAA)
      INTEGER S(IAA)
      write (6,*) ' lsbrr is called '
      call exit
      end

c
      double precision function e1 (x)
      double precision x
      ifail = 0
      e1 = s13aaf(x, ifail)
      return
      end
C
*
      FUNCTION H1RN(DUMMY)
*
*#**********************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN FSU-SCRI-87-50 AND MODIFIED BY F. JAMES 1988 TO *
*# PRODUCE VECTOR OF NUMBERS.                                          *
*# ENTRIES ARE:                                                        *
*#     FUNCTION    H1RN(DUMMY)     SINGLE RANDOM NUMBER                *
*#     SUBROUTINE  H1RNV(VEC,LEN)  VECTOR OF RANDOM NUMBERS            *
*#     SUBROUTINE  H1RNIN(IJ,KL)   INITIALISE WITH SEEDS               *
*#     SUBROUTINE  H1RNIV(VEC)     INITIALISE/RESTART WITH SEED ARRAY  *
*#     SUBROUTINE  H1RNSV(VEC)     SAVE SEED ARRAY VEC(100)            *
*#                                                                     *
*# NOTE: -H1RNIN OR H1RNIV MUST BE CALLED BEFORE GENERATING ANY        *
*#        RANDOM NUMBER(S).                                            *
*#       -H1RNSV SAVES SEED ARRAY INTO VEC(100) ONLY. THE USER HAS TO  *
*#        OUTPUT IT.                                                   *
*#                                                                     *
*# CHANGED BY: G. GRINDHAMMER AT: 90/03/14                             *
*# REASON :                                                            *
*#**********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*16    FLAG,CHECK
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN
cpb         WRITE(6,*) ' H1RN (RANMAR): INITIALIZED WITH DEFAULT SEED'
            ISEED1      = 12345
            ISEED2      = 98765
            CALL H1RNIN(ISEED1,ISEED2)
            FIRST = .FALSE.
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
  100 CONTINUE
      H1RN = U(I)-U(J)
      IF(H1RN .LT. 0.) H1RN = H1RN + 1.
      U(I) = H1RN
      I = I - 1
      IF( I .EQ. 0) I=97
      J = J - 1
      IF( J .EQ. 0) J=97
      C = C - CD
      IF( C .LT. 0) C = C + CM
      H1RN = H1RN-C
      IF(H1RN .LE. 0.) H1RN = H1RN + 1.
      IF(H1RN .GE. 1.) GOTO 100
      RETURN
      END
*
      SUBROUTINE H1RNV(RVEC,LEN)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8          RVEC(1)
      CHARACTER*16    FLAG,CHECK
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN
            WRITE(6,*) ' H1RNV (RANMAR): CALL H1RNIN OR H1RNIV BEFORE',
     >                 ' CALLING H1RN.'
            STOP
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
      DO 200 IVEC = 1,LEN
  190    CONTINUE
         UNI = U(I)-U(J)
         IF(UNI .LT. 0.) UNI = UNI + 1.
         U(I) = UNI
         I = I - 1
         IF( I .EQ. 0) I=97
         J = J - 1
         IF( J .EQ. 0) J=97
         C = C - CD
         IF( C .LT. 0) C = C + CM
         UNI = UNI-C
         IF(UNI .LE. 0.) UNI = UNI + 1.
         IF(UNI .GE. 1.) GOTO 190
         RVEC(IVEC) = UNI
  200 CONTINUE
      RETURN
      END
*
      SUBROUTINE H1RNIN(IJ,KL)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*16    FLAG
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
      IJ = IABS(IJ)
      KL = IABS(KL)
      IJ = MOD(IJ,31329)
      KL = MOD(KL,30082)
      I  = MOD(IJ/177, 177) + 2
      J  = MOD(IJ, 177)     + 2
      K  = MOD(KL/169, 178) + 1
      L  = MOD(KL, 169)
*
      DO 300 II= 1, 97
         S= 0.
         T= 0.5
         DO 250 JJ= 1,24
            M = MOD(MOD(I*J,179)*K, 179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1, 169)
            IF ( MOD(L*M,64) .GE. 32) S = S + T
  250       T = 0.5*T
  300 U(II) = S
      C  =   362436./16777216.
      CD =  7654321./16777216.
      CM = 16777213./16777216.
      IP = 97
      JP = 33
*
      FLAG = 'H1RN INITIALISED'
cpb   WRITE(6,610) IJ,KL,I,J,K,L
  610 FORMAT(' H1RNIN: H1RN (RANMAR) INITIALISED WITH 2 SEEDS: ',
     >       2I7,4I4)
*
      RETURN
      END
*
      SUBROUTINE H1RNIV(VEC)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8          VEC(100)
      CHARACTER*16    FLAG
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
      DO 400 IC = 1, 97
  400 U(IC) = VEC(IC)
      C  = VEC(98)
      CD =  7654321./16777216.
      CM = 16777213./16777216.
      IP = NINT(VEC(99))
      JP = NINT(VEC(100))
*
      FLAG = 'H1RN INITIALISED'
      WRITE(6,*) ' H1RNIV: H1RN (RANMAR) INITIALISED/RESTARTED WITH',
     >           ' SEED ARRAY VEC(100)'
*
      RETURN
      END
*
      SUBROUTINE H1RNSV(VEC)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8          VEC(100)
      COMMON /RASET1/ U(97),C,CD,CM,I,J
*
      DO 10 IC = 1, 97
   10 VEC(IC) = U(IC)
      VEC(98) = C
      VEC(99) = REAL(I)
      VEC(100)= REAL(J)
      RETURN
      END
c
c
      function erf (x)
      implicit real*8 (a-h,o-z)
c  s15aef is a nag routine: error function
c  erf(x)=2/sqrt(pi)*int^x_0 dt exp(-(t**2))
      erf=s15aef(x,ifail)
      return
      end
