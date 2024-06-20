subroutine interpolate_extra_fields()
#include "compile_opt.inc"
  use all_variables, only : interp_data2, zones, global_parameters, global_variables, flags
  use Minterpolation_types
  implicit none
  Type(Interpolation_temp) :: temp
  integer*4 :: k,n,nfields
  integer*4 :: Nx,Nz,ind
  nfields=(global_parameters%N_ions+1)
  allocate(temp%knots_val(interp_data2%N_knots,nfields))
  allocate(temp%zones(global_parameters%N_Zones))
  do k=1,global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(temp%zones(k)%val(1:Nx,1:Nz,1:nfields))
     ind=0
     do n=0,global_parameters%N_ions
        ind=ind+1
        temp%zones(k)%val(1:Nx,1:Nz,ind)=zones(k)%species(n)%sources%rad(1:Nx,1:Nz)
     end do
  end do
  call interpolate(temp,nfields)
  do k=1,Interp_Data2%N_knots
     ind=0
     do n=0,global_parameters%N_ions
        ind=ind+1
        Interp_Data2%knots_radiation(k,n)=temp%knots_val(k,ind)
     end do
  end do
  do k=1,global_parameters%N_Zones
     deallocate(temp%zones(k)%val)
  end do
  deallocate(temp%knots_val)
#if VORTICITY_PASTIX == 1
  !interpolate phi
  nfields=3
  allocate(temp%knots_val(interp_data2%N_knots,nfields))
  do k=1,global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(temp%zones(k)%val(1:Nx,1:Nz,1:nfields))
     temp%zones(k)%val(1:Nx,1:Nz,1)=zones(k)%electric_fields(1)%phi(1:Nx,1:Nz)
     temp%zones(k)%val(1:Nx,1:Nz,2)=zones(k)%electric_fields(1)%vorticity(1:Nx,1:Nz)
     temp%zones(k)%val(1:Nx,1:Nz,3)=zones(k)%electric_fields(1)%E(1:Nx,1:Nz) 
  end do
  call interpolate(temp,nfields)
  do k=1,Interp_Data2%N_knots
     Interp_Data2%knots_phi(k)=temp%knots_val(k,1)
     Interp_Data2%knots_vorticity(k)=temp%knots_val(k,2)
     Interp_Data2%knots_Epara(k)=temp%knots_val(k,3)
  end do
  do k=1,global_parameters%N_Zones
     deallocate(temp%zones(k)%val)
  end do
  deallocate(temp%knots_val)
#endif
  if(flags%turbulence_model.eq.1) then
     !interpolate k - epsilon
     nfields=2
     allocate(temp%knots_val(interp_data2%N_knots,nfields))
     do k=1,global_parameters%N_zones
        Nx=Zones(k)%mesh%Nx
        Nz=Zones(k)%mesh%Nz
        allocate(temp%zones(k)%val(1:Nx,1:Nz,1:nfields))
        temp%zones(k)%val(1:Nx,1:Nz,1)=zones(k)%kepsilon(1)%k(1:Nx,1:Nz)
        temp%zones(k)%val(1:Nx,1:Nz,2)=zones(k)%kepsilon(1)%epsilon(1:Nx,1:Nz)
     end do
     call interpolate(temp,nfields)
     do k=1,Interp_Data2%N_knots
        Interp_Data2%knots_k(k)=temp%knots_val(k,1)
        Interp_Data2%knots_epsilon(k)=temp%knots_val(k,2)
     end do
     do k=1,global_parameters%N_Zones
        deallocate(temp%zones(k)%val)
     end do
     deallocate(temp%knots_val)
  end if
end subroutine interpolate_extra_fields
