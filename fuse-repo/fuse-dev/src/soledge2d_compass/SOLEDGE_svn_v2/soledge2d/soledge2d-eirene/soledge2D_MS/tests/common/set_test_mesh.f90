subroutine set_test_mesh(level)
  use all_variables, only : global_parameters,zones,megazones
  use Mzone
  implicit none
  integer*4,intent(in) :: level
  integer*4 :: Nx,Nz
  integer*4 :: i,j,k
  global_parameters%N_zones=2
  allocate(zones(1:global_parameters%N_zones))
  global_parameters%N_megazones=2
  allocate(megazones(1:global_parameters%N_megazones))
  ! set mesh
  Nx=2**(level+1)
  Nz=2**(level+1)
  do k=1,2
     zones(k)%mesh%Nx=Nx/2
     zones(k)%mesh%Nz=Nz
     call allocate_zone_mesh(zones(k)%mesh)
     do i=1,Nx/2
        do j=1,Nz
           zones(k)%mesh%x_minus_1half(i,j)=0.5*(k-1)+real(i-1)/real(Nx/2)*0.5d0
           zones(k)%mesh%x_plus_1half(i,j)=0.5*(k-1)+real(i)/real(Nx/2)*0.5d0
           zones(k)%mesh%z_minus_1half(i,j)=real(j-1)/real(Nz)
           zones(k)%mesh%z_plus_1half(i,j)=real(j)/real(Nz)
           zones(k)%mesh%x(i,j)=(zones(k)%mesh%x_minus_1half(i,j)+zones(k)%mesh%x_plus_1half(i,j))*0.5D0
           zones(k)%mesh%z(i,j)=(zones(k)%mesh%z_minus_1half(i,j)+zones(k)%mesh%z_plus_1half(i,j))*0.5D0
           zones(k)%mesh%xmin=0.5d0*(k-1)
           zones(k)%mesh%xmax=0.5d0*k
           zones(k)%mesh%zmin=0.D0
           zones(k)%mesh%zmax=1.D0
        end do
     end do
  end do
  zones(1)%Neighbors(N_North)=2
  zones(1)%Neighbors(N_South)=-1
  zones(1)%Neighbors(N_East)=1
  zones(1)%Neighbors(N_West)=1
  zones(2)%Neighbors(N_North)=-1
  zones(2)%Neighbors(N_South)=1
  zones(2)%Neighbors(N_East)=2
  zones(2)%Neighbors(N_West)=2
  zones(1)%MagNeighbors(N_North)=0
  zones(1)%MagNeighbors(N_South)=1
  zones(1)%MagNeighbors(N_East)=0
  zones(1)%MagNeighbors(N_West)=0
  zones(2)%MagNeighbors(N_North)=1
  zones(2)%MagNeighbors(N_South)=0
  zones(2)%MagNeighbors(N_East)=0
  zones(2)%MagNeighbors(N_West)=0
  do k=1,2
     megazones(k)%size=1
     allocate(megazones(k)%zone_number(1:megazones(k)%size))
     megazones(k)%zone_number(1)=k
     megazones(k)%is_periodic=.true.
  end do
end subroutine set_test_mesh
