

      FUNCTION FEPLOT3 (IROT,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS
      USE PHOTON

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IROT, K
      REAL(DP) :: PLS, ADD, ELCX, FEPLOT3
      INTEGER :: II, KK, IPLSTI

      FEPLOT3=0.D0
      KK=NELROT(IROT)
      IPLSTI=MPLSTI(IPLS)
      SELECT CASE (KK)
      CASE (-2)
         FEPLOT3=EPLOT3(IROT,1,1)+EDRIFT(IPLS,K)
      CASE (-3)
         FEPLOT3=1.5*TIIN(IPLSTI,K)+EDRIFT(IPLS,K)
      CASE DEFAULT
        WRITE (6,*) ' FEPLOT3 NOT READY '
        WRITE (6,*) ' CALCULATION ABANDONNED '
        CALL EXIT_OWN(1)
      END SELECT

      RETURN
      END