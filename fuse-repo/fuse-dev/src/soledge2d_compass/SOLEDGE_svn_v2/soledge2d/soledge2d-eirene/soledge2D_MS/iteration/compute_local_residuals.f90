subroutine compute_local_residuals(zone)
  use all_variables, only : global_parameters,flags
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%n_ions
     zone%species(n)%residuals%resn=maxval(abs(zone%species(n)%var(2)%density(1:Nx,1:Nz)&
          -zone%species(n)%var(1)%density(1:Nx,1:Nz))*(1.D0-zone%masks%chi2(1:Nx,1:Nz)))
     zone%species(n)%residuals%resG=maxval(abs(zone%species(n)%var(2)%Gamma(1:Nx,1:Nz)&
          -zone%species(n)%var(1)%Gamma(1:Nx,1:Nz))*(1.D0-zone%masks%chi2(1:Nx,1:Nz)))
     if(flags%solve_temperature) then
        zone%species(n)%residuals%resT=maxval(abs(zone%species(n)%var(2)%temperature(1:Nx,1:Nz)&
             -zone%species(n)%var(1)%temperature(1:Nx,1:Nz))*(1.D0-zone%masks%chi2(1:Nx,1:Nz)))
     else
        zone%species(n)%residuals%resT=0.D0
     end if
  end do
end subroutine compute_local_residuals
