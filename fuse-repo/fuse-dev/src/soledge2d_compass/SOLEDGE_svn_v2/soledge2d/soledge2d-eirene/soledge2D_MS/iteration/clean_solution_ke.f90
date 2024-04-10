subroutine clean_solution_ke(zone)
  use all_variables, only : global_parameters, global_variables
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=0,Nx+1
     do j=0,Nz+1
        if(zone%kepsilon(2)%k(i,j).lt.1.d-6) then
           zone%kepsilon(2)%k(i,j)=1.d-6
        end if
        if(zone%kepsilon(2)%epsilon(i,j).lt.1.d-6) then
           zone%kepsilon(2)%epsilon(i,j)=1.d-6
        end if
     end do
  end do
end subroutine clean_solution_ke
