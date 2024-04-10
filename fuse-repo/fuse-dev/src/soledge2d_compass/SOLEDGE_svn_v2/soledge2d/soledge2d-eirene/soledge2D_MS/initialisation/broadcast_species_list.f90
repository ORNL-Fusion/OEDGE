subroutine broadcast_species_list(zone)
  use MZone
  use all_variables, only : global_parameters
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 n
  zone%species(0)%element%symbol='e'
  zone%species(0)%element%name='electron'
  zone%species(0)%element%mass=0.
  zone%species(0)%element%mass2=m_e/m_u
  zone%species(0)%charge=-1.
  do n=1,global_parameters%N_ions
     zone%species(n)%element%symbol=global_parameters%element_list(global_parameters%ions_list(n,1))%symbol
     zone%species(n)%element%name=global_parameters%element_list(global_parameters%ions_list(n,1))%name
     zone%species(n)%element%mass=global_parameters%element_list(global_parameters%ions_list(n,1))%mass
     zone%species(n)%element%mass2=global_parameters%element_list(global_parameters%ions_list(n,1))%mass2
     zone%species(n)%element%Z=global_parameters%element_list(global_parameters%ions_list(n,1))%Z
     zone%species(n)%element_index=global_parameters%ions_list(n,1)
     zone%species(n)%charge=global_parameters%ions_list(n,2)
     if(zone%species(n)%charge.eq.zone%species(n)%element%Z) then
        zone%species(n)%compute_ionization=.false.
        ! ... because this species does not ionize
     else
        zone%species(n)%compute_ionization=.true.
     end if
     if(zone%species(n)%charge.eq.1) then
        zone%species(n)%compute_recombination=.false.
        ! ... because this is computed by eirene
     else
        zone%species(n)%compute_recombination=.true.
     end if

  end do
end subroutine broadcast_species_list
