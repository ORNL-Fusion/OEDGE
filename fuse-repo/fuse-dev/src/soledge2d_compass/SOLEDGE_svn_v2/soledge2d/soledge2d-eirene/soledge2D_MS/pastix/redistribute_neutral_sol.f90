subroutine redistribute_neutral_sol(CSC_)
  use all_variables, only : global_parameters, zones
  use Mpastix_solve
  implicit none
  Type(CSC) :: CSC_
  integer*4 i,j,k
  integer*4 Nx,Nz
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=0,Nx+1
        do j=0,Nz+1
           Zones(k)%neutrals%density(i,j)=CSC_%b(Zones(k)%mesh%index(i,j))
        end do
     end do
  end do
end subroutine redistribute_neutral_sol
