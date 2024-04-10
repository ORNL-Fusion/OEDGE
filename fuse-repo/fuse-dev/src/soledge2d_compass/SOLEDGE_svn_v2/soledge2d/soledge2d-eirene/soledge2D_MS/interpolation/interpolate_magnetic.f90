subroutine interpolate_magnetic()
  use all_variables, only : interp_data2, zones, global_parameters
  use Minterpolation_types
  implicit none
  Type(Interpolation_temp) :: temp
  integer*4 :: k,n,nfields
  integer*4 :: Nx,Nz,ind
  nfields=4
  allocate(temp%knots_val(interp_data2%N_knots,nfields))
  allocate(temp%zones(global_parameters%N_Zones))
  do k=1,global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(temp%zones(k)%val(1:Nx,1:Nz,1:nfields))
     ind=0
     ind=ind+1
     temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%mesh%Br(1:Nx,1:Nz)
     ind=ind+1
     temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%mesh%Bz(1:Nx,1:Nz)
     ind=ind+1
     temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%mesh%Bphi(1:Nx,1:Nz)
     ind=ind+1
     temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%mesh%B(1:Nx,1:Nz)
  end do
  call interpolate(temp,nfields)
  do k=1,Interp_Data2%N_knots
     ind=0
     ind=ind+1
     Interp_Data2%knots_Br(k)=temp%knots_val(k,ind)
     ind=ind+1
     Interp_Data2%knots_Bz(k)=temp%knots_val(k,ind)
     ind=ind+1
     Interp_Data2%knots_Bphi(k)=temp%knots_val(k,ind)
     ind=ind+1
     Interp_Data2%knots_B(k)=temp%knots_val(k,ind)
  end do
  do k=1,global_parameters%N_Zones
     deallocate(temp%zones(k)%val)
  end do
end subroutine interpolate_magnetic
