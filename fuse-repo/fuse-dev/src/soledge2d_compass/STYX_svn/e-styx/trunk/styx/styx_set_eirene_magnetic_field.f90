subroutine styx_set_eirene_magnetic_field
  use all_variables, only : interp_data2
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  implicit none

  integer*4 :: itor
  real*8, dimension(:,:), allocatable :: BX3D, BY3D, BZ3D


  allocate(BX_tri(Neir_cells))
  allocate(BY_tri(Neir_cells))
  allocate(BZ_tri(Neir_cells))
  allocate(BF_tri(Neir_cells))

  allocate(BX3D(NRKNOT_styx,Ntor_cells+1),BY3D(NRKNOT_styx,Ntor_cells+1), &
		BZ3D(NRKNOT_styx,Ntor_cells+1))
  	
  do itor=1,Ntor_cells+1
    BX3D(:,itor)=Interp_data2%Knots_Br
    BY3D(:,itor)=Interp_data2%Knots_Bz
    BZ3D(:,itor)=Interp_data2%Knots_Bphi
  enddo

  call styx_vertices_to_centre(BX3D,BX_tri)
  call styx_vertices_to_centre(BY3D,BY_tri)
  call styx_vertices_to_centre(BZ3D,BZ_tri)

  BF_tri=sqrt(BX_tri**2+BY_tri**2+BZ_tri**2)
  
  where(BF_tri>0._dp)
    BX_tri=BX_tri/BF_tri
    BY_tri=BY_tri/BF_tri
    BZ_tri=BZ_tri/BF_tri
  elsewhere
    BX_tri=1._dp/sqrt(3._dp)
    BY_tri=1._dp/sqrt(3._dp)
    BZ_tri=1._dp/sqrt(3._dp)
  end where

  deallocate(BX3D,BY3D,BZ3D)

end subroutine styx_set_eirene_magnetic_field
