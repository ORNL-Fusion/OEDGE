subroutine clean_solution2(zone)
  use all_variables, only : global_parameters, global_variables
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     do i=0,Nx+1
        do j=0,Nz+1
           if(zone%species(n)%var(2)%temperature(i,j).lt.global_variables%min_temperature) then
              zone%species(n)%var(2)%temperature(i,j)=global_variables%min_temperature
           end if
        end do
     end do
  end do
end subroutine clean_solution2
