subroutine compute_index()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 i,j,k
  integer*4 Nx,Nz
  integer*4 index_
  index_=1
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=0,Nx+1
        do j=0,Nz+1
           Zones(k)%mesh%index(i,j)=index_
           index_=index_+1
        end do
     end do
  end do
end subroutine compute_index
