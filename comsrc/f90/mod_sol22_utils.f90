module mod_sol22_utils

  !use mod_sol22_sources

  implicit none


  public :: getcs_sol22,getcs_sol22_dbl


contains




  REAL FUNCTION GetCs_sol22(te,ti)
    use mod_solparams
    use mod_solcommon
    !use mod_params
    !use mod_comtor
    IMPLICIT none
    !     INCLUDE 'params'
    !     INCLUDE 'comtor'
    REAL, INTENT(IN) :: te,ti

    !     Te,i in eV
    !     a    in amu
    !     result is m s-1

    !GetCs_sol22 = 9.78817E+03 * SQRT(0.5 * (1.0 + zb) * (te + ti) / mb)
    GetCs_sol22 = 9.78817E+03 * SQRT(0.5 * (1.0 + zb) * (te + ti) / mb * (econv/mconv))

    !write(0,*) 'GETCS:',getcs_sol22,te,ti,mb,zb,9.78817E+03,sqrt(econv/mconv),SQRT(0.5 * (1.0 + zb) * (te + ti)*econv / (mb*mconv))

    RETURN
  END FUNCTION GetCs_sol22

  REAL*8 FUNCTION GetCs_sol22_dbl(te,ti)
    use mod_solparams
    use mod_solcommon
    !use mod_params
    !use mod_comtor
    IMPLICIT none
    !     INCLUDE 'params'
    !     INCLUDE 'comtor'
    REAL*8, INTENT(IN) :: te,ti

    !     Te,i in eV
    !     a    in amu
    !     result is m s-1

    !GetCs_sol22_dbl = 9.78817D+03 * SQRT(0.5D0 * (1.0D0 + zb) * (te + ti) / mb)
    GetCs_sol22_dbl = SQRT(0.5D0 * (1.0D0 + zb) * (te + ti) / mb * (econv/mconv))

    !write(0,*) 'GETCS_DBL:',getcs_sol22_dbl,te,ti,mb,zb,9.78817E+03,sqrt(econv/mconv),SQRT(0.5 * (1.0 + zb) * (te + ti)*econv / (mb*mconv))

    RETURN
  END FUNCTION GetCs_sol22_dbl
    
end module mod_sol22_utils
