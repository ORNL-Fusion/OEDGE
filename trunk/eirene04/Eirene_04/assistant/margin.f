C
C
      SUBROUTINE MARGIN (NLFLD,NDIM,NMIN,NMAX,LOG)
C
C  FIND MINIMAL AND MAXIMAL LOGICAL OF VALUE LOG IN THE
C  LOGICAL FIELD NLFLD
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      INTEGER, INTENT(OUT) :: NMIN, NMAX
      LOGICAL, INTENT(IN) :: NLFLD(NDIM)
      LOGICAL, INTENT(IN) :: LOG
      INTEGER :: J
C
      IF (LOG) THEN
C
         DO 1 J=1,NDIM
            NMIN=J
            IF (NLFLD(J)) GOTO 2
    1    CONTINUE
         NMIN=NDIM+1
C
    2    NMAX=NDIM
         DO 3 J=NDIM,NMIN,-1
            NMAX=J
            IF (NLFLD(J)) GOTO 4
    3    CONTINUE
C
      ELSE
C
         DO 5 J=1,NDIM
            NMIN=J
            IF (.NOT.NLFLD(J)) GOTO 6
    5    CONTINUE
         NMIN=NDIM+1
C
    6    NMAX=NDIM
         DO 7 J=NDIM,NMIN,-1
            NMAX=J
            IF (.NOT.NLFLD(J)) GOTO 4
    7    CONTINUE
C
      ENDIF
C
    4 RETURN
      END
