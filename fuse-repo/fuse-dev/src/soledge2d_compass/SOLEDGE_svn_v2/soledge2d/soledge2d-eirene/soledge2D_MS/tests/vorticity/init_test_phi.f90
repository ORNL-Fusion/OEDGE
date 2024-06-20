subroutine init_test_phi()
  use test_var
  use all_variables, only : zones, global_parameters,reference_parameters
  use Mphysics
  implicit none
  integer*4 :: k,n
  real*8 :: r,theta,test_theta,test_r
  integer*4 :: i,j,Nx,Nz
  real*8 :: Lambda
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do i=0,Nx+1
        do j=0,Nz+1
           !phi = lambda * Te
           Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+Zones(k)%species(1)%var(1)%temperature(i,j)&
                /Zones(k)%species(0)%var(1)%temperature(i,j)))
           zones(k)%electric_fields(1)%phi(i,j)=Lambda*zones(k)%species(0)%var(1)%temperature(i,j)
        end do
     end do
  end do
end subroutine init_test_phi
