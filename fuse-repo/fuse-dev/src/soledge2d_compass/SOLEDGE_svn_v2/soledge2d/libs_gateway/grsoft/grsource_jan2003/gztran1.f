************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
       SUBROUTINE gzdotran(a, b, m, type)
*  AUTHOR   Z. ZOWIRUCHA
*  CHANGED. M. BUSCH    12.4.91
*                       INTEGER    --> INTEGER
*                       REAL   *4  --> REAL
************************************************************************
*
       IMPLICIT NONE
************************************************************************
*  declare dummy arguments
*
       REAL       a(1:3)
       REAL       b(1:3)
       REAL       m(1:4,1:4)
       INTEGER    type
*-----------------------------------------------------------------------
*
       INTEGER    i
*
*-----------------------------------------------------------------------
       IF(type.EQ.1) THEN
*
       DO 100 i=1, 3

                b(i)= m(1,i)*a(1) +
     +                m(2,i)*a(2) +
     +                m(3,i)*a(3) + m(4,i)
 100   CONTINUE
       ELSE
       DO 102 i=1, 3
*
                b(i)= m(1,i)*a(1) +
     +                m(2,i)*a(2) +
     +                m(3,i)*a(3)
*
 102   CONTINUE
       ENDIF
*
*-------------------------------------------------------   end gzdotran
       END
*
*
*
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzrot(matrix, axis, alfa)
************************************************************************
*
*  dummy arguments:
*    matrix - rotation transformation matrix;
*    axis   - axis of rotation (1- x, 2- y, 3 -z);
*    alfa   - gzangle of rotation;
*
*  function:  this subroutine computes a rotation matrix (in three
*             dimentions) with gzangle alfa about x, y or z axis;
*
************************************************************************
       IMPLICIT NONE
*  declare dummy arguments
       REAL       matrix(1:4, 1:4)
       INTEGER    axis
       REAL       alfa
*
*  declare local variables
       REAL       s, c
       INTEGER    i, j
       INTEGER    m, m1, m2
*
*-----------------------------------------------------------------------
       DO 100 i=1, 4
         DO 110 j=1, 4
            matrix(j, i)=0.0
 110     CONTINUE
 100   CONTINUE
       m=axis
*
       matrix(4, 4)=1.0
       matrix(m, m)=1.0
       m1 = MOD(m,3) + 1
       m2 = MOD(m1,3) + 1
*
       c = COS(alfa)
       s = SIN(alfa)
       matrix(m1, m1)=c
       matrix(m2, m2)=c
       matrix(m1, m2)=s
       matrix(m2, m1)=-s
* ---------------------------------------------------   end gzrot
       END
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzshift(matrix, shift1, shift2, shift3)
*
*
       REAL   matrix(1:4,1:4), shift1, shift2, shift3
*
*
       DO 100 i=1,4
         DO 110 j=1,4
            IF(j.EQ.i) THEN
               matrix(j,i)=1.0
            ELSE
               matrix(j,i)=0.0
            ENDIF
*
  110    CONTINUE
  100  CONTINUE
*
       matrix(4,1)=shift1
       matrix(4,2)=shift2
       matrix(4,3)=shift3
*
       RETURN
       END
*
*
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzaxis(matrix, p1, p2, theta)
*
*
       REAL     matrix(1:4,1:4), p1(1:3), p2(1:3), theta
*
*
       REAL     B(4,4), C(4,4), D(4,4), E(4,4), BI(4,4), CI(4,4)
       REAL     DI(4,4),   alfa,beta
       REAL     x1, y1, z1,  x2, y2, z2
*
       x1=p1(1)
       y1=p1(2)
       z1=p1(3)
*
       x2=p2(1)
       y2=p2(2)
       z2=p2(3)
*
       CALL gzshift(B, -x1,-y1,-z1)
       beta= -gzangle(x2-x1, y2-y1)
*
       CALL gzrot(C, 3, beta)
       alfa= -gzangle( SQRT( (x2-x1)**2 + (y2-y1)**2), z1-z2)
       CALL gzrot(D, 2, alfa)
       CALL gzrot(E, 3, theta)
*
       CALL gzshift(BI,  x1, y1, z1)
       CALL gzrot(CI, 3, -beta)
       CALL gzrot(DI,2,  -alfa)
*
       CALL gzmult(C,B,C)
       CALL gzmult(D,C,D)
       CALL gzmult(E,D,E)
       CALL gzmult(DI,E,DI)
       CALL gzmult(CI,DI,CI)
       CALL gzmult(BI,CI,matrix)
*
       RETURN
       END
*
*
*
*
*
       REAL FUNCTION gzangle(A,B)
*
*
       REAL     A, B
*
*
*
*
       REAL    PI
       DATA    PI/3.141592653/
*
       IF( ABS(A).NE.0.0) THEN
          gzangle= ATAN(B/A)


       IF(A.LT.0.0)  gzangle=gzangle+pi
       ELSE
         IF(B.LT.0.0) THEN
            gzangle= pi + pi/2.0
         ELSE
            gzangle=pi/2.0
         ENDIF
         IF(ABS(B).EQ.0.0) gzangle = 0.0
       ENDIF
*
*
*
       RETURN
       END
