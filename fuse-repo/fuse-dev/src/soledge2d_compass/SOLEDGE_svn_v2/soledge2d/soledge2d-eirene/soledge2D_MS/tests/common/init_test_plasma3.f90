subroutine init_test_plasma3()
  use test_var
  use all_variables, only : zones, global_parameters,reference_parameters
  use Mphysics
  implicit none
  integer*4 :: k,n
  integer*4 :: nion
  do k=1,global_parameters%N_Zones
     do nion=1,global_parameters%N_ions
        !for all ions
        zones(k)%species(nion)%var(1)%density=1.d-10
        zones(k)%species(nion)%var(1)%Gamma=0.D0
        zones(k)%species(nion)%var(1)%temperature=1.D0
     end do
     ! singly charged ion
     zones(k)%species(1)%var(1)%density=1.D0
     !electrons
     zones(k)%species(0)%var(1)%temperature=1.D0
  end do
end subroutine init_test_plasma3
