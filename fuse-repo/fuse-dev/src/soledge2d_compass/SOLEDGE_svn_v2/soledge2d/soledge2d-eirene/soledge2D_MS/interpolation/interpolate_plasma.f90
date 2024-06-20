subroutine interpolate_plasma()
  use all_variables, only : interp_data2, zones, global_parameters, global_variables
  use Minterpolation_types
  implicit none
  Type(Interpolation_temp) :: temp
  integer*4 :: k,n,nfields
  integer*4 :: Nx,Nz,ind
  nfields=(global_parameters%N_ions+1)*4
  allocate(temp%knots_val(interp_data2%N_knots,nfields))
  allocate(temp%zones(global_parameters%N_Zones))
  do k=1,global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(temp%zones(k)%val(1:Nx,1:Nz,1:nfields))
     ind=0
     do n=0,global_parameters%N_ions
        ind=ind+1
        temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%species(n)%var(1)%density(1:Nx,1:Nz)
        ind=ind+1
        temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%species(n)%var(1)%velocity(1:Nx,1:Nz)
        ind=ind+1
        temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%species(n)%var(1)%temperature(1:Nx,1:Nz)
        ind=ind+1
        temp%zones(k)%val(1:Nx,1:Nz,ind)=Zones(k)%species(n)%sources%rad(1:Nx,1:Nz)
     end do
  end do
  call interpolate(temp,nfields)
  do k=1,Interp_Data2%N_knots
     ind=0
     do n=0,global_parameters%N_ions
        ind=ind+1
        Interp_Data2%knots_density(k,n,1)=max(temp%knots_val(k,ind),global_variables%min_density)
        ind=ind+1
        Interp_Data2%knots_velocity(k,n,1)=temp%knots_val(k,ind)
        ind=ind+1
        Interp_Data2%knots_temperature(k,n,1)=max(temp%knots_val(k,ind),global_variables%min_temperature)
        ind=ind+1
        Interp_Data2%knots_radiation(k,n)=temp%knots_val(k,ind)
     end do
  end do
  do k=1,global_parameters%N_Zones
     deallocate(temp%zones(k)%val)
  end do
end subroutine interpolate_plasma
