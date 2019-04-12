c     -*-Fortran-*-
c from outlos3D.f
      subroutine transform_vect(m,a)
      implicit none
      real*8 m(3,3),a(3)
c
c     Transform_VECT: rotate the vector to get its value in the new
c                     coordinate system - in which the rotation to map 
c                     the camera view will be entirely around the rotated
c                     X-axis.
c
      integer i,j
      real*8 tmpa(3) 
c
c      write(6,'(a,3(1x,g12.5))') 'START VECT:',(a(i),i=1,3)
c
      do i = 1,3 
c      
         tmpa(i) = m(i,1)*a(1) + m(i,2)*a(2)+m(i,3)*a(3)
c
      end do
c
c     Assign result to a
c
      do i = 1,3
c
         a(i) = tmpa(i)
c
      end do 
c
c      write(6,'(a,3(1x,g12.5))') 'END VECT  :',(a(i),i=1,3)
c
c      write(6,'(a,12(1x,g12.5))') 'M-TRANS:',
c     >                     ((m(i,j),j=1,3),i=1,3)
c
      return
      end
c
c
c
c
c ======================================================================
c (from comsrc/utility_com.f -SL, 09.01.07)
c 
      SUBROUTINE FITTER (N1,X1,F1,N2,X2,F2,METHOD)
      implicit none
      INTEGER N1,N2
      REAL    X1(N1),F1(N1),X2(N2),F2(N2)
      CHARACTER*(*) METHOD
C
C  *********************************************************************
C  *                                                                   *
C  *  FITTER:   THIS ROUTINE INTERPOLATES BETWEEN A SET OF GIVEN DATA  *
C  *  POINTS AND EXTRAPOLATES OFF THE ENDS OF THE RANGE.  SEVERAL      *
C  *  METHODS ARE USED ACCORDING TO THE VALUE OF "METHOD":-            *
C  *                                                                   *
C  *  "SPLINE"  THIS OPTION USES THE HARWELL CUBIC SPLINE GENERATORS   *
C  *  TO MODEL A SYSTEM OF A FEW DATA POINTS (N1,X1,F1), FROM WHICH WE *
C  *  INTERPOLATE / EXTRAPOLATE TO THE REQUIRED NEW SYSTEM OF MORE     *
C  *  DATAPOINTS (N2,X2,F2).  ALL ARGUMENTS ARE INPUT EXCEPT F2, WHICH *
C  *  IS CALCULATED IN THIS ROUTINE.                                   *
C  *                                                                   *
C  *  "LINEAR"  LINEAR INTERPOLATION                                   *
C  *                                                                   *
C  *  CHRIS FARRELL  (HUNTERSKIL)  SEPT 88                             *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER I1,I2
      REAL FDASH1(1001),WORK(3003),TG01B
C-----------------------------------------------------------------------
      IF     (METHOD.EQ.'SPLINE') THEN
        IF (N1.GT.1001) GOTO 1001
        IF (N1.GT.1) THEN
          CALL TB04A (N1,X1,F1,FDASH1,WORK)
          IF (WORK(1).GT.1.E-6) GOTO 1001
        ENDIF
        DO 100 I2 = 1, N2
          IF     (X2(I2).LE.X1(1)) THEN
            F2(I2) = F1(1)
          ELSEIF (X2(I2).GE.X1(N1)) THEN
            F2(I2) = F1(N1)
          ELSE
            F2(I2) = TG01B (-1,N1,X1,F1,FDASH1,X2(I2))
          ENDIF
  100   CONTINUE
C-----------------------------------------------------------------------
      ELSEIF (METHOD.EQ.'LINEAR') THEN
        DO 210 I2 = 1, N2
          IF     (X2(I2).LE.X1(1)) THEN
            F2(I2) = F1(1)
C           WRITE (6,'(1X,F7.4,'' ='',F7.4)') F2(I2),F1(1)
          ELSEIF (X2(I2).GE.X1(N1)) THEN
            F2(I2) = F1(N1)
C           WRITE (6,'(1X,F7.4,'' ='',F7.4)') F2(I2),F1(N1)
          ELSE
            I1 = 2
  200       IF (X2(I2).GT.X1(I1)) THEN
              I1 = I1 + 1
              GOTO 200
            ENDIF
            F2(I2) = F1(I1) - (X1(I1)-X2(I2)) * (F1(I1)-F1(I1-1)) /
     >                                          (X1(I1)-X1(I1-1))
C           WRITE (6,9001) F2(I2),F1(I1),X1(I1),X2(I2),F1(I1),F1(I1-1),
C    >                                                 X1(I1),X1(I1-1)
          ENDIF
  210   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
      RETURN
C
 1001 CALL PRC ('FITTER:  CUBIC SPLINE ERROR.  CHECK N1 <= 1001')
      STOP
C
 9001 FORMAT(1X,F7.4,' =',F7.4,' - (',F7.4,' -',F7.4,') * (',F7.4,
     >  ' -',F7.4,') / (',F7.4,' -',F7.4,')')
      END
c
c ======================================================================
c (from comsrc/harw.f -SL, 09.01.07)
C
C     TB04A
C     =====
C     IBM  : HARWELL LIBRARY ROUTINE TO CREATE A CUBIC SPLINE
C     CRAY : USE HSL ROUTINE WITH DOUBLE PRECISION NAME TB04AD.
C
      SUBROUTINE TB04A (N,X,F,FDASH,WORK)
      implicit none
      INTEGER N,I
      REAL X(N),F(N),FDASH(N),WORK(3*N)
C
C     AUTOMATIC X2,F2,FDASH2,WORK2
      DOUBLE PRECISION X2(1001),F2(1001),FDASH2(1001),WORK2(3003)
C
C     IN ORDER TO INTERFACE THESE ROUTINES PROPERLY, SINCE
C     FORTRAN DOES NOT DO ARGUMENT TYPE CHECKING IT IS
C     NECESSARY TO CONVERT THE SINGLE PRECISION VARIABLES TO
C     DOUBLE PRECISION. THIS IS COMMENTED OUT FOR CRAY USAGE
C     SINCE ONLY REALS ARE USED IN THIS ENVIRONMENT (64 BIT)
C
      DO 10 I = 1,N
         X2(I) = DBLE(X(I))
         F2(I) = DBLE(F(I))
         FDASH2(I) = DBLE(FDASH(I))
         WORK2(I) = DBLE(WORK(I))
         WORK2(N+I) = DBLE(WORK(N+I))
         WORK2(2*N+I) = DBLE(WORK(2*N+I))
10    CONTINUE

      CALL TB04AD(N,X2,F2,FDASH2,WORK2)

      DO 20 I = 1,N
         X(I) = SNGL(X2(I))
         F(I) = SNGL(F2(I))
         FDASH(I) = SNGL(FDASH2(I))
         WORK(I) = SNGL(WORK2(I))
         WORK(N+I) = SNGL(WORK2(N+I))
         WORK(2*N+I) = SNGL(WORK2(2*N+I))
20    CONTINUE

      RETURN
      END
c
c ======================================================================
c  (from comsrc/harw.f -SL, 09.01.07)
c
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS TB04AD
      SUBROUTINE TB04AD(N,X,F,D,A)
      implicit none
      integer n
C STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION A,D,F,H1,H2,P,X
      DIMENSION X(N),F(N),D(N),A(*)
      integer np
      DATA NP/6/
c
c     Local variables
c
      integer i,j,k
c
C F(I) ARE THE FUNCTION VALUES AT THE POINTS X(I) FOR I=1,N AND
C THE SPLINE DERIVATIVES D(I) ARE FOUND.  THE DIMENSION OF A MUST
C NOT BE LESS THAN 3*N. PERIPHERAL NP MUST BE AN OUTPUT MEDIUM.
      DO 5 I=2,N
      IF(X(I)-X(I-1))1,1,5
1     WRITE(NP,3)I
3      FORMAT(29H RETURN FROM TB04AD BECAUSE X,I3,13H OUT OF ORDER)
      A(1)=1.0D0
      RETURN
5     CONTINUE
      DO 30 I=1,N
      J=2
      IF(I-1)6,10,6
6     J=N-1
      IF(I.EQ.N)GO TO 10
      H1=1.0D0/(X(I)-X(I-1))
      H2=1.0D0/(X(I+1)-X(I))
      A(3*I-2)=H1
      A(3*I-1)=2.0D0*(H1+H2)
      A(3*I)=H2
      D(I)=3.0D0*(F(I+1)*H2*H2+F(I)*(H1*H1-H2*H2)-F(I-1)*H1*H1)
      GO TO 30
10    H1=1.0D0/(X(J)-X(J-1))
      H2=1.0D0/(X(J+1)-X(J))
      A(3*I-2)=H1*H1
      A(3*I-1)=H1*H1-H2*H2
      A(3*I)=-H2*H2
      D(I)=2.D0*(H1*H1*H1*(F(J)-F(J-1))+H2*H2*H2*(F(J)-F(J+1)))
30    CONTINUE
      P=A(4)/A(1)
      A(5)=A(5)-P*A(2)
      A(6)=A(6)-P*A(3)
      D(2)=D(2)-P*D(1)
      DO 50 I=3,N
      K=3*I-4
      P=A(K+2)/A(K)
      A(K+3)=A(K+3)-P*A(K+1)
      D(I)=D(I)-P*D(I-1)
      IF(I.NE.N-1)GO TO 50
      P=A(K+5)/A(K)
      A(K+5)=A(K+6)-P*A(K+1)
      A(K+6)=A(K+7)
      D(N)=D(N)-P*D(N-2)
50    CONTINUE
      D(N)=D(N)/A(3*N-1)
      DO 60 I=3,N
      J=N+2-I
60    D(J)=(D(J)-A(3*J)*D(J+1))/A(3*J-1)
      D(1)=(D(1)-D(2)*A(2)-D(3)*A(3))/A(1)
      A(1)=0.0D0
      RETURN
      END
c
c ======================================================================
c (from comsrc/harw.f -SL, 09.01.07)
C
C***********************************************************************
C
C     TG01B
C     =====
C     IBM  : HARWELL LIBRARY FUNCTION TO EXTRACT VALUE ALONG SPLINE.
C     CRAY : USE HSL FUNCTION WITH DOUBLE PRECISION NAME TG01BD.
C
      REAL FUNCTION TG01B (IFLAG,N,X,F,FDASH,X0)
      implicit none
      integer n,iflag
      REAL X(N),F(N),FDASH(N),X0
C     AUTOMATIC X2,F2,FDASH2,X02
      DOUBLE PRECISION X2(1001),F2(1001),FDASH2(1001),X02,TG01BD
      EXTERNAL TG01BD
      INTEGER I
C
C     THIS IS AN INTERFACE BETWEEN A REAL AND DOUBLE PRECISION
C     FUNCTION. THE CONVERSIONS WERE NOT REQUIRED ON A CRAY WHERE
C     A REAL WAS EQUIVALENT TO DOUBLE PRECISION ON A SMALLER MACHINE
C     AND THE CRAY DOUBLE PRECISION WAS DEFAULTED TO REAL. HOWEVER,
C     THIS DIFFERENCE IS SIGNIFICANT ON A WORKSTATION.
C
C     D.ELDER   1991 APRIL 24
C
      X02 = DBLE(X0)
      DO 10 I = 1,N
         X2(I) = DBLE(X(I))
         F2(I) = DBLE(F(I))
         FDASH2(I) = DBLE(FDASH(I))
10    CONTINUE

      TG01B = SNGL(TG01BD (IFLAG,N,X2,F2,FDASH2,X02))

      X0 = SNGL(X02)
      DO 20 I = 1,N
         X(I) = SNGL(X2(I))
         F(I) = SNGL(F2(I))
         FDASH(I) = SNGL(FDASH2(I))
20    CONTINUE

      RETURN
      END
C
c
c ======================================================================
c (from comsrc/harw.f -SL, 09.01.07)
c
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS TG01BD
      DOUBLE PRECISION FUNCTION TG01BD(IX,N,U,S,D,X)
      implicit none
C
C**********************************************************************
C
C      TG01BD - FUNCTION ROUTINE TO EVALUATE A CUBIC SPLINE GIVEN SPLINE
C     VALUES AND FIRST DERIVATIVE VALUES AT THE GIVEN KNOTS.
C
C     THE SPLINE VALUE IS DEFINED AS ZERO OUTSIDE THE KNOT RANGE,WHICH
C     IS EXTENDED BY A ROUNDING ERROR FOR THE PURPOSE.
C
C                 F = TG01BD(IX,N,U,S,D,X)
C
C       IX    ALLOWS CALLER TO TAKE ADVANTAGE OF SPLINE PARAMETERS SET
C             ON A PREVIOUS CALL IN CASES WHEN X POINT FOLLOWS PREVIOUS
C             X POINT. IF IX < 0 THE WHOLE RANGE IS SEARCHED FOR KNOT
C             INTERVAL; IF IX > 0 IT IS ASSUMED THAT X IS GREATER THAN
C             THE X OF THE PREVIOUS CALL AND SEARCH STARTED FROM THERE.
C       N     NUMBER OF KNOTS.
C       U     THE KNOTS.
C       S     THE SPLINE VALUES.
C       D     THE FIRST DERIVATIVE VALUES OF THE SPLINE AT THE KNOTS.
C       X     THE POINT AT WHICH THE SPLINE VALUE IS REQUIRED.
C       F     THE VALUE OF THE SPLINE AT THE POINT X.
C
C                                      MODIFIED JULY 1970
C
C**********************************************************************
C
      integer iflg,ieps,ix,n,j
      DOUBLE PRECISION A,B,D,H,Q1,Q2,S,SS,U,X,Z
C
C     ALLOWABLE ROUNDING ERROR ON POINTS AT EXTREAMS OF KNOT RANGE
C     IS 2**IEPS*MAX(!U(1)!,!U(N)!).
c slmod begin
c...  How did this ever work? -SL 19.08.06
      DIMENSION U(N),S(N),D(N)
c
c      DIMENSION U(1),S(1),D(1)
c slmod end
      DATA IFLG/0/,IEPS/-50/
C
C       TEST WETHER POINT IN RANGE.
      IF(X.LT.U(1)) GO TO 990
      IF(X.GT.U(N)) GO TO 991
C
C       JUMP IF KNOT INTERVAL REQUIRES RANDOM SEARCH.
      IF(IX.LT.0.OR.IFLG.EQ.0) GO TO 12
C       JUMP IF KNOT INTERVAL SAME AS LAST TIME.
      IF(X.LE.U(J+1)) GO TO 8
C       LOOP TILL INTERVAL FOUND.
    1 J=J+1
   11 IF(X.GT.U(J+1)) GO TO 1
      GO TO 7
C
C       ESTIMATE KNOT INTERVAL BY ASSUMING EQUALLY SPACED KNOTS.
   12 J=DABS(X-U(1))/(U(N)-U(1))*(N-1)+1
C       ENSURE CASE X=U(N) GIVES J=N-1.
      J=MIN0(J,N-1)
C       INDICATE THAT KNOT INTERVAL INSIDE RANGE HAS BEEN USED.
      IFLG=1
C       SEARCH FOR KNOT INTERVAL CONTAINING X.
      IF(X.GE.U(J)) GO TO 11
    2 J=J-1
      IF(X.LT.U(J)) GO TO 2
C
C       CALCULATE SPLINE PARAMETERS FOR JTH INTERVAL.
    7 H=U(J+1)-U(J)
      Q1=H*D(J)
      Q2=H*D(J+1)
      SS=S(J+1)-S(J)
      B=3D0*SS-2D0*Q1-Q2
      A=Q1+Q2-2D0*SS
C
C       CALCULATE SPLINE VALUE.
    8 Z=(X-U(J))/H
      TG01BD=((A*Z+B)*Z+Q1)*Z+S(J)
      RETURN
C       TEST IF X WITHIN ROUNDING ERROR OF U(1).
  990 IF(X.LE.U(1)-2D0**IEPS*DMAX1(DABS(U(1)),DABS(U(N)))) GO TO 99
      J=1
      GO TO 7
C       TEST IF X WITHIN ROUNDING ERROR OF U(N).
  991 IF(X.GE.U(N)+2D0**IEPS*DMAX1(DABS(U(1)),DABS(U(N)))) GO TO 99
      J=N-1
      GO TO 7
   99 IFLG=0
C       FUNCTION VALUE SET TO ZERO FOR POINTS OUTSIDE THE RANGE.
      TG01BD=0D0
      RETURN
      END
c
c ======================================================================
c (based on PRC in comsrc/utility_com.f -SL, 09.01.07)
C
C  *********************************************************************
C  *  PRC:  PRINTS A CHARACTER STRING                                  *
C  *********************************************************************
C
      SUBROUTINE PRC(STRING)
      implicit none
c
      integer stringlen,lenstr,datunit
      external lenstr
      CHARACTER STRING*(*)
c
      datunit = 0

      WRITE (datunit,'(1X,A)') STRING(1:LEN_TRIM(string))
      RETURN
      END
c
c ======================================================================
c
      INTEGER FUNCTION CH1(STRING)
      IMPLICIT none
      CHARACTER STRING*(*)
      INTEGER   I1
      CH1=1
      DO I1 = 1, LEN_TRIM(STRING)
        IF (STRING(I1:I1).NE.' ') EXIT
        CH1=CH1+1
      ENDDO
      RETURN
      END
c
c ======================================================================
c
c ======================================================================
c
c function: CalcPerp2
c
c Calculate the perpendicular distance from a point to a line, but never
c mind the check if the intersection point is between the end points (combine
c with CalcPerp in the future to avoid code repetition).
c
      REAL*8 FUNCTION CalcPerp2(a,b,c,t) 

      IMPLICIT none

      REAL*8, INTENT(IN)  :: a(3),b(3),c(3)
      REAL*8, INTENT(OUT) :: t

      REAL*8 p(3),delta(3),dist

      REAL*8     DTOL
      PARAMETER (DTOL = 1.0D-7)


      CalcPerp2 = -1.0D0

      delta = b - a

      IF (DABS(delta(1)).GT.DTOL.OR.DABS(delta(2)).GT.DTOL.OR.
     .    DABS(delta(3)).GT.DTOL) THEN

        t = ((c(1) - a(1)) * delta(1) + (c(2) - a(2)) * delta(2) + 
     .       (c(3) - a(3)) * delta(3)) /
     .      (delta(1)**2 + delta(2)**2 + delta(3)**2)

        p(1:3) = a(1:3) + t * delta(1:3)

        IF (.TRUE.) THEN
c       IF ((t+DTOL).GE.0.0D0.AND.(t-DTOL).LE.1.0D0) THEN

          dist = DSQRT((p(1) - c(1))**2 + (p(2) - c(2))**2 +
     .                 (p(3) - c(3))**2)

c          WRITE(0,*) 'A=',a
c          WRITE(0,*) 'B=',b
c          WRITE(0,*) 'C=',c
c          WRITE(0,*) 'P=',p
c          WRITE(0,*) 'DELTA=',delta

          IF (dist.LT.DTOL*10.0D0) THEN
c           Point C is on the line AB:
            CalcPerp2 = 0.0D0
          ELSE
c           Point of perpendicular intersection is displaced from
c           the point C:
            CalcPerp2 = dist
          ENDIF
        ELSE
          CalcPerp2 = -1.0D0
        ENDIF
      ELSE
c       If the points are all identicle, then return a positive result,
c       otherwise indicate that the problem was ill-posed:
        IF (DABS(a(1) - c(1)).LT.DTOL.AND.DABS(a(2) - c(2)).LT.DTOL.AND.
     .      DABS(a(3) - c(3)).LT.DTOL) THEN
          CalcPerp2 =  0.0D0
        ELSE
          CalcPerp2 = -1.0D0
        ENDIF
      ENDIF

      RETURN
      END

