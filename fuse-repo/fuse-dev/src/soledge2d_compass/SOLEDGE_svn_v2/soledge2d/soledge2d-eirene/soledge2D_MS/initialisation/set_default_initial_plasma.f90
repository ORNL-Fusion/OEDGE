subroutine set_default_initial_plasma()
#include "compile_opt.inc"
  use all_variables, only : zones, global_parameters
  implicit none
  integer*4 :: k,n
  do k=1,global_parameters%N_Zones
     !electrons
     zones(k)%species(0)%var(1)%temperature=1.D0
     !main ions
     zones(k)%species(1)%var(1)%density=1.D0
     zones(k)%species(1)%var(1)%Gamma=0.D0
     zones(k)%species(1)%var(1)%temperature=1.D0
     zones(k)%species(1)%penalisation_memories%alpham=0.D0
     zones(k)%species(1)%penalisation_memories%alphap=0.D0
     !other ions
     do n=2,global_parameters%N_ions
        if(zones(k)%species(n)%charge.eq.zones(k)%species(n)%element%Z) then
           zones(k)%species(n)%var(1)%density=1.D-8
           zones(k)%species(n)%var(1)%Gamma=0.D0
           zones(k)%species(n)%var(1)%temperature=1.D0
        else
           zones(k)%species(n)%var(1)%density=1.D-8
           zones(k)%species(n)%var(1)%Gamma=0.D0
           zones(k)%species(n)%var(1)%temperature=1.D0
        end if
        zones(k)%species(n)%penalisation_memories%alpham=0.D0
        zones(k)%species(n)%penalisation_memories%alphap=0.D0
     end do
  end do
end subroutine set_default_initial_plasma
