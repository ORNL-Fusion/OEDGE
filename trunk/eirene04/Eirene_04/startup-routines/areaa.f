C  OCT 03: return best guess, even if R out of range, rather then exit
C
      FUNCTION AREAA (R,N,ARCA,Y,EP1R,ELLR)
C
C   FOR LEVGEO=2 OPTION
C
C   THIS FUNCTION INTERPOLATES LINEAR STANDARD MESH PARAMETERS ON POINTS
C   LYING BETWEEN THE SURFACES OF THIS MESH
C   AND EVALUATES THE ENCLOSED (CROSS SECTIONAL) AREA AND THE ARCLENGTH
C
C  INPUT:
C     R
C     N MUST BE GIVEN SUCH THAT RSURF(N)<=R<=RSURF(N+1)
C  OUTPUT:
C     AREAA(R)=AREA WITHIN EIRENE-SURFACE LABELED BY R (RSURF)
C                   I.E. CROSS SECTIONAL AREA (CM**2)
C     ARCA(R)=ARCLENGTH AT R
C     Y(R)
C     EP1R(R)
C     ELLR(R)
C     TRIA(R) (TO BE WRITTEN)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CGRID

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: R 
      REAL(DP), INTENT(OUT) :: ARCA, Y, EP1R, ELLR
      INTEGER, INTENT(IN) :: N
      REAL(DP) :: RRN, APB, RRI, RRD, XL1QQ, Q1, XL1, XL1Q, AREAA
      INTEGER :: NP
C
      NP=N+1
      IF (N.LE.0.OR.N.GE.NR1STM) GOTO 998
      IF (R.LT.RSURF(N).OR.R.GT.RSURF(NP)) GOTO 999
      RRI=RSURF(N)
      RRD=RSURF(NP)-RRI
      RRN=(R-RRI)/RRD
C
      ELLR=ELL(N)+RRN*(ELL(NP)-ELL(N))
      EP1R=EP1(N)+RRN*(EP1(NP)-EP1(N))
 100  Y=ELLR*R
C
      APB=R+Y+EPS60
      XL1=(R-Y)/APB
      XL1Q=XL1*XL1
      XL1QQ=XL1Q*XL1Q
      Q1=(16.-0.75*XL1QQ)/(64.-16.*XL1Q)
C
      AREAA=PIA*Y*R
      ARCA=APB*Q1*PIA*4.
C
      RETURN
C
998   WRITE (6,*) 'ERROR IN FUNCTION AREAA: N, NR1STM ',N,NR1STM
      CALL EXIT_OWN(1)
999   WRITE (6,*) 'ERROR IN FUNCTION AREAA: R,N ',R,N
      IF (R.LT.RSURF(N)) THEN
        ELLR=ELL(N)
        EP1R=EP1(N)
      ELSE
        ELLR=ELL(NP)
        EP1R=EP1(NP)
      ENDIF
      GOTO 100
      END