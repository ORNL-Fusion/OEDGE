subroutine init_feedback_plasma()
  use all_variables, only : reference_parameters, global_parameters, zones
  use MradialFeedback
  implicit none
  integer*4 :: i,j,k,n
  integer*4 :: Nx,Nz
  real*8 :: x,D,keep
  real*8 :: interpolate_feedback
  real*8 :: dens,temp
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do i=1,Nx
        x=zones(k)%mesh%x(i,1)
        dens=interpolate_feedback(x,4) ! n
        do n=0,global_parameters%N_ions
           zones(k)%species(n)%var(1)%density(i,:)=dens/reference_parameters%fields%n0
        end do
        temp=interpolate_feedback(x,6) ! Ti
        do n=1,global_parameters%N_ions
           zones(k)%species(n)%var(1)%temperature(i,:)=temp/reference_parameters%fields%T0eV
        end do
        temp=interpolate_feedback(x,5) ! Te
        zones(k)%species(0)%var(1)%temperature(i,:)=temp/reference_parameters%fields%T0eV
     end do
  end do
end subroutine init_feedback_plasma
