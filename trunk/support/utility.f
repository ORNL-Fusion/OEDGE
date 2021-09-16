      INTEGER FUNCTION LENSTR (ASTR)
      implicit none
      CHARACTER*(*) ASTR
C
C  *********************************************************************
C  *                                                                   *
C  *  LENSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *
C  *          ANY TRAILING BLANKS.                                     *
C  *                                                                   *
C  *********************************************************************
C
      integer i
c
      DO 10 I = LEN(ASTR),1,-1
         IF (ASTR(I:I) .NE. ' ') THEN
            LENSTR = I
            RETURN
         ENDIF
   10 CONTINUE
      LENSTR = 1
      RETURN
      END
C
C
C
      INTEGER FUNCTION EXTSTR (ASTR,start)
      implicit none
      CHARACTER*(*) ASTR
      integer start,i
c
      extstr = 0
C
C  *********************************************************************
C  *                                                                   *
C  *  EXTSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *
C  *          ANY TRAILING BLANKS AND THE EFFECTIVE STARTING POINT BY  *
c  *          IGNORING LEADING BLANKS.                                 *
C  *                                                                   *
C  *********************************************************************
C
      DO 10 I = LEN(ASTR),1,-1
         IF (ASTR(I:I) .NE. ' ') THEN
            extSTR = I
            goto 20
         ENDIF
   10 CONTINUE

 20   if (extstr.gt.1) then
         do i = 1, extstr
            if (astr(I:I).ne.' ') then
               start = i
               return
            endif
         end do
      else
         extSTR = 1
         start  = 1
      endif
c
      RETURN
      END
C
C  *********************************************************************
C  *  RZERO:  ZEROES A REAL ARRAY  ...     C.M.FARRELL   NOV 1987      *
C  *********************************************************************
C
      SUBROUTINE RZERO (RARRAY, N)
      implicit none
      INTEGER I,N
      REAL RARRAY(N)
      DO 100 I = 1, N
        RARRAY(I) = 0.0
  100 CONTINUE
      RETURN
      END
