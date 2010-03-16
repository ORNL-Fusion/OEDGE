C
C
      SUBROUTINE REFLEC_photon
C
C  REFLECT ESCAPING PHOTONS
C  INPUT:
C       ILREF = 0  PERFECT ABSORPTION
C       ILREF = 1  HOLLMAN DATABASE
C
C       ITYP  = 0  INCIDENT PHOTON
C  OUTPUT:
C     LGPART= TRUE AND:
C       ITYP = 0  PHOTON IPHOT IS RETURNED TO CALLING PROGRAM
C     LGPART= FALSE  NO PARTICLE IS RETURNED (ABSORBTION)
C       ITYP = -1
C
      USE PRECISION
      use comprt

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: WMIN, XMP, XCP
      INTEGER, INTENT(IN) :: NPRIN
      INTEGER, INTENT(INOUT) :: IGASF, IGAST
C
C---------------------------------------------------------------------
C

      ENTRY REFLC0_PHOTON
C
      write (iunout,*) ' REFLC0_PHOTON: no calculation performed '
      write (iunout,*) ' photonic reflection not yet available '
C
      RETURN
C
      ENTRY REFLC1_PHOTON (WMIN,XMP,XCP,NPRIN,IGASF,IGAST)
C
      write (iunout,*) ' REFLC1_PHOTON: no calculation performed '
      write (iunout,*) ' photonic reflection not yet available '
C
      RETURN
C
      END
