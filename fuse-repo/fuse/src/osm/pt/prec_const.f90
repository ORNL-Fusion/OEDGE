!---------------------------------------------
! File : Prec_Const.f90
! Date : 08/07/2004
! Constantes
!---------------------------------------------

module prec_const
  implicit none
! Precision des reels
  integer, parameter :: float = SELECTED_REAL_KIND(13,99)
  
! Constantes reelles predefinies
  REAL(float),PARAMETER :: ZE     = 0.0_float
  REAL(float),PARAMETER :: HF     = 0.5_float
  REAL(float),PARAMETER :: ON     = 1.0_float
  REAL(float),PARAMETER :: TW     = 2.0_float
  REAL(float),PARAMETER :: TH     = 3.0_float
  REAL(float),PARAMETER :: FO     = 4.0_float
  REAL(float),PARAMETER :: FI     = 5.0_float
  REAL(float),PARAMETER :: PI     = 3.141592653589793238462643383279502884197_float
end module prec_const

