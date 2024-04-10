subroutine update_fields(zone)
  use all_variables, only : global_parameters, flags
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=1,global_parameters%n_ions
     zone%species(n)%var(1)%density=zone%species(n)%var(2)%density
     zone%species(n)%var(1)%Gamma=zone%species(n)%var(2)%Gamma
  end do
  if(flags%solve_temperature) then
     do n=0,global_parameters%n_ions
        zone%species(n)%var(1)%temperature=zone%species(n)%var(2)%temperature
     end do
  end if
  if(flags%turbulence_model.eq.1) then
     zone%kepsilon(1)%k=max(zone%kepsilon(2)%k,1.d-6)
     zone%kepsilon(1)%epsilon=zone%kepsilon(2)%epsilon
  end if
end subroutine update_fields
