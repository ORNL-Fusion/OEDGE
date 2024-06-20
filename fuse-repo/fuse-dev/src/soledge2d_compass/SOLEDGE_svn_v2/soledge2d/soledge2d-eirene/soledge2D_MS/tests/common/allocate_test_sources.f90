subroutine allocate_test_sources()
  use test_var
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: n,k,Nx,Nz
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do n=0,global_parameters%N_ions
        allocate(test_sources(k,n)%Sn(1:Nx,1:Nz))
        allocate(test_sources(k,n)%SG(1:Nx,1:Nz))
        allocate(test_sources(k,n)%ST(1:Nx,1:Nz))
        allocate(test_sources(k,n)%SW(1:Nx,1:Nz))
     end do
  end do
end subroutine allocate_test_sources
