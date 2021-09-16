C
C
      SUBROUTINE PL3Q(CORD,N,IO,NF)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: CORD(*)
      INTEGER, INTENT(IN) :: N, IO
      LOGICAL, INTENT(IN) :: NF
      REAL(DP) :: XP(N+1),YP(N+1)
      REAL(SP) :: XPS(N+1),YPS(N+1)
      INTEGER :: J, JJ, I, N3

      N3=3*N
      I=0
      DO 100 J=1,N3,3
        I=I+1
        CALL PL3D (CORD(J),CORD(J+1),CORD(J+2),XP(I),YP(I))
100   CONTINUE
      XP(N+1)=XP(1)
      YP(N+1)=YP(1)
C
      IF (IO.GE.2) CALL GRNWPN(IO)
      do 200 jj=1,n+1
        xps(jj)=xp(jj)
        yps(jj)=yp(jj)
200   continue
      CALL GRLN (XPS,YPS,N+1)
C  AUSFUELLEN DES KURVENZUGES MIT FARBE NO. IO
      IF (NF) CALL GRFILL(N+1,XPS,YPS,1,1)
      IF (IO.GE.2) CALL GRNWPN(1)
      RETURN
      END
