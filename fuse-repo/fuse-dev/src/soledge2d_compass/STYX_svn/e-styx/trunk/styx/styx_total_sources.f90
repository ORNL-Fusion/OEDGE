subroutine styx_total_sources(data_in,data_intg)
use eirmod_precision
  use styx2eirene
  implicit none
 
  real(dp), dimension(Neir_cells), intent(in) :: data_in
  real(dp), intent(out) :: data_intg
  integer :: ITRI

  data_intg=0._dp
 
  do ITRI=1,Neir_cells
  	data_intg = data_intg + vol_tri_eirene(ITRI)*data_in(ITRI)
  enddo


end subroutine
