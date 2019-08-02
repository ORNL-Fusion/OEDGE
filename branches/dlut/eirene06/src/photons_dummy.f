C EIRENE06 COMPILATION
C ===== SOURCE: phot_reflec.f
      module phot_reflec
      use precision
      use comprt
      implicit none

      private

      public :: init_refl_hlm, reflect_hollmann

      contains

      subroutine init_refl_hlm()
      integer :: i

      write (iunout,*) ' init_refl_hlm: nothing done '
      write (iunout,*) ' photonic reflection not yet available '
      return

      end subroutine init_refl_hlm


! reflect
! Die Routine entscheidet, ob ein Photon mit Einfallwinkel theta_i und
! Wellenlaenge lambda_i reflektiert wird. Im Falle einer
! Reflektion ist flag 1, ansonsten 0.
      subroutine reflect_hollmann(theta_i, lambda_i, mat, flag,
     .                   theta_out, alpha_out, rprob)
      real(dp), intent(in) :: theta_i, lambda_i
      integer, intent(in) :: mat
      real(dp), intent(out) :: theta_out, alpha_out, rprob
      integer, intent(out) :: flag

      theta_out = 0._dp
      alpha_out = 0._dp
      rprob = 0._dp
      flag = 0

      write (iunout,*) ' reflect_hollmann: nothing done '
      write (iunout,*) ' photonic reflection not yet available '
      return

      end subroutine reflect_hollmann

      end module phot_reflec




C ===== SOURCE: reflec_photon.f
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
