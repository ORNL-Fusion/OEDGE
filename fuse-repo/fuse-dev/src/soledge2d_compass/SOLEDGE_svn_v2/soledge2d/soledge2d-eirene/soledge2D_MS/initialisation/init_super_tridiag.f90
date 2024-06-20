subroutine init_super_tridiag
  use all_variables, only : global_parameters, megazones,&
       zones, flux_surfaces
  implicit none
  integer*4 :: k,izone,n
  integer*4 :: ind,length
  integer*4 :: Nx,i
  integer*4 :: n_zones
  !counting flux surfaces
  global_parameters%N_psi_surfaces=0
  do k=1,global_parameters%N_megazones
     izone=megazones(k)%zone_number(1)
     Nx=zones(izone)%mesh%Nx
     global_parameters%N_psi_surfaces=global_parameters%N_psi_surfaces+Nx
  end do
  allocate(flux_surfaces(1:global_parameters%N_psi_surfaces))
  !filling info about flux surfaces
  ind=1
  do k=1,global_parameters%N_megazones
     length=0
     Nx=zones(megazones(k)%zone_number(1))%mesh%Nx
     do n=1,megazones(k)%size
        izone=megazones(k)%zone_number(n)
        length=length+zones(izone)%mesh%Nz
     end do
     do i=1,Nx
        allocate(flux_surfaces(ind)%tridiag%a(1:length))
        allocate(flux_surfaces(ind)%tridiag%b(1:length))
        allocate(flux_surfaces(ind)%tridiag%c(1:length))
        allocate(flux_surfaces(ind)%tridiag%S(1:length))
        allocate(flux_surfaces(ind)%tridiag%buffer(0:length+1))
        flux_surfaces(ind)%tridiag%buffer=0.d0
        flux_surfaces(ind)%properties%is_periodic=megazones(k)%is_periodic
        flux_surfaces(ind)%ns_psi=ind
        flux_surfaces(ind)%Nz=length
        flux_surfaces(ind)%properties%i_psi=i
        flux_surfaces(ind)%properties%n_zones=megazones(k)%size
        n_zones=flux_surfaces(ind)%properties%n_zones
        allocate(flux_surfaces(ind)%properties%zones(1:n_zones))
        do n=1,n_zones
           flux_surfaces(ind)%properties%zones(n)=megazones(k)%zone_number(n)
        end do
        ind=ind+1
     end do
  end do
end subroutine init_super_tridiag
