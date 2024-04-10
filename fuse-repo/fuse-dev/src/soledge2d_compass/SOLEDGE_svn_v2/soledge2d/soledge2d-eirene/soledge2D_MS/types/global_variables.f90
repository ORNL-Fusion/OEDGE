module Mglobal_variables

  type :: Tresiduals
     real*8 :: resn
     real*8 :: resG
     real*8 :: resT
  end type Tresiduals
  
  type :: Tglobal_variables
     real*8,allocatable :: dt_table(:)
     real*8 :: dt
     real*8 :: dt_vort
     real*8 :: dt_vort_old
     real*8 :: sign_metric
     real*8 :: min_density
     real*8 :: min_temperature
     real*8 :: Teps
     real*8 :: tempus
     real*8 :: eta_para_smooth
     real*8 :: epsilon_weno
     type(Tresiduals),allocatable :: residuals(:)
  end type Tglobal_variables

  type :: Tpenalisation_parameters
     real*8 :: eta
     real*8 :: eta2
     real*8 :: dump_pen
     real*8 :: gain
     integer*4 :: keep
  end type Tpenalisation_parameters

  type :: Tintegrals
     real*8,allocatable :: flux_tot_out_ac(:)
     real*8,allocatable :: flux_tot_in_ac(:)
     real*8,allocatable :: flux_totE_out_ac(:)
     real*8,allocatable :: flux_totE_in_ac(:)
     real*8,allocatable :: source_n(:)
     real*8,allocatable :: source_E(:)
     real*8,allocatable :: variation_E(:)
     real*8,allocatable :: stored_E(:)
     real*8,allocatable :: variation_N(:)
     real*8,allocatable :: stored_N(:)
     real*8,allocatable :: source_ionz_tot(:)
     real*8,allocatable :: source_E_ionz(:)
     real*8,allocatable :: tot_N(:),tot_E(:)
     real*8,allocatable :: total_radiation(:)
  end type Tintegrals

end module Mglobal_variables
