subroutine styx_vertices_to_centre_minT(data_in,data_out)
  use all_variables, only : reference_parameters
  use eirmod_precision
  use styx2eirene
  use eirmod_comprt, only : iunout
  implicit none

  real(dp), dimension(1:NRKNOT_styx,Ntor_cells+1), intent(in) :: data_in
  real(dp), dimension(Neir_cells), intent(out) :: data_out
  real(dp) temp(3)
  real(dp) T0

  integer :: ITRI,itor,icell,itor2

  T0=reference_parameters%fields%T0eV

  do itor=1,Ntor_cells
    do ITRI=1,Ntri_styx
     icell = itri + (itor-1)*(Ntri_styx+1)
     temp(1)=data_in(NVERT(1,ITRI),itor)
     temp(2)=data_in(NVERT(2,ITRI),itor)
     temp(3)=data_in(NVERT(3,ITRI),itor)
     data_out(icell)=max(minval(temp),3.d-1/T0)
    enddo
  enddo

  if (is_3D) then
    do itor=1,Ntor_cells
      do ITRI=1,Ntri_styx
       icell = itri + (itor-1)*(Ntri_styx+1)
       itor2 = itor+1
       temp(1)=data_in(NVERT(1,ITRI),itor2)
       temp(2)=data_in(NVERT(2,ITRI),itor2)
       temp(3)=data_in(NVERT(3,ITRI),itor2)     
       data_out(icell)=0.5d0*(data_out(icell)+max(minval(temp),3.d-1/T0))
      enddo
    enddo
  endif

end subroutine styx_vertices_to_centre_minT 
