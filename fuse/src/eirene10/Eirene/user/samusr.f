C
C
      SUBROUTINE SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE COMUSR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6
      REAL(DP), INTENT(OUT) :: X0,Y0,Z0,TEWL,TIWL(*),DIWL(*),
     .                         VXWL(*),VYWL(*),VZWL(*),
     .                         EFWL(*), SHWL, WEISPZ(*)
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IER

      entry sm0usr (is1,is2,sorad1,sorad2,sorad3,sorad4,sorad5,sorad6)
      return

      entry SM1USR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)

      x0 = 0._dp
      y0 = 0._dp
      z0 = 0._dp

      irusr = 0
      ipusr = 0
      itusr = 0
      iausr = 0
      ibusr = 0

      tiwl(1:nplsti) = 0._dp
      tewl = 0._dp
      diwl(1:nplsi) = 0._dp
      vxwl(1:nplsv) = 0._dp
      vywl(1:nplsv) = 0._dp
      vzwl(1:nplsv) = 1._dp
      efwl(1:nplsi) = 0._dp
      shwl = 0._dp
      weispz(1:nspz) = 0._dp

      RETURN
      END
