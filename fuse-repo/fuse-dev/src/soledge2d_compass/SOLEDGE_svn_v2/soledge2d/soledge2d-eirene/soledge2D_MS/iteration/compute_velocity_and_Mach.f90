subroutine compute_velocity_and_Mach(zone,STEP)
  use all_variables, only : global_parameters,global_variables
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n,Nx,Nz,i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=1,global_parameters%N_ions
     do i=0,Nx+1
        do j=0,Nz+1
           if(zone%species(n)%var(STEP)%density(i,j).lt.global_variables%min_density*0.1d0) then
              zone%species(n)%var(STEP)%velocity(i,j)=0.d0
              zone%species(n)%var(STEP)%Mach(i,j)=0.d0
           else
              zone%species(n)%var(STEP)%velocity(i,j)=zone%species(n)%var(STEP)%Gamma(i,j)&
                   /zone%species(n)%var(STEP)%density(i,j)
              zone%species(n)%var(STEP)%Mach(i,j)=zone%species(n)%var(STEP)%velocity(i,j)/&
                   sqrt((zone%species(n)%var(STEP)%temperature(i,j)+zone%species(n)%charge*zone%species(0)%var(STEP)%temperature(i,j))/&
                   zone%species(n)%element%mass)
           end if
        end do
     end do
  end do
end subroutine compute_velocity_and_Mach
