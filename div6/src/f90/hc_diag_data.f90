module hc_diag_data

  private

  integer :: num_states
    real,allocatable :: total_entered_state(:)
    real,allocatable :: ave_kin_energy(:)
    real,allocatable :: ave_energy_released(:)
    real,allocatable :: ave_tot_energy_gained(:)
    real,allocatable :: ave_transitions(:)

    public::init_hc_diag_data,end_hc_diag_data,print_hc_diag_data
    public::record_state_energy_diag_data

contains


  subroutine init_hc_diag_data(number_hc_species)
    implicit none
    integer :: number_hc_species,flag

    num_states = number_hc_species

    call allocate_diag_data


  end subroutine init_hc_diag_data

  subroutine end_hc_diag_data
    implicit none

    call deallocate_diag_data

  end subroutine end_hc_diag_data



  subroutine record_state_energy_diag_data(cur_hc_spec,hc_v,sput_weight,e_hc_gained,kin_energy_added,cur_hc_mass,&
                                           &transition_count)
    use hc_velocity_type
    implicit none
    integer :: cur_hc_spec
    type(hc_velocity_type1) :: hc_v
    real :: sput_weight
    real :: e_hc_gained,kin_energy_added,cur_hc_mass
    integer :: transition_count

    total_entered_state(cur_hc_spec) = total_entered_state(cur_hc_spec) + sput_weight
    !
    ! average kinetic energy in eV 
    !
    ! [1/2 (m/e)] ** (-1/2) = [0.5 * 1.67e-27/1.6e-19] ** (-1/2) = 1.38e4
    !
    ave_kin_energy(cur_hc_spec) = ave_kin_energy(cur_hc_spec) + cur_hc_mass * (hc_v%vtot/1.38e4)**2 * sput_weight
    ave_energy_released(cur_hc_spec) = ave_energy_released(cur_hc_spec) + e_hc_gained * sput_weight
    ave_tot_energy_gained(cur_hc_spec) = ave_tot_energy_gained(cur_hc_spec) + kin_energy_added * sput_weight
    ave_transitions(cur_hc_spec) = ave_transitions(cur_hc_spec) + real(transition_count) * sput_weight

    write(6,'(a,i6,10(1x,g12.5))') 'RECORD:',cur_hc_spec,sput_weight,cur_hc_mass,hc_v%vtot,mag_v(hc_v%v),&
        &   cur_hc_mass * (hc_v%vtot/1.38e4)**2 ,e_hc_gained,kin_energy_added,real(transition_count),ave_transitions(cur_hc_spec)


  end subroutine record_state_energy_diag_data



  subroutine print_hc_diag_data(ounit)
    implicit none
    integer :: ounit
    
    !
    ! Print out hc diagnostic data
    !
    
    call print_state_energy_diag_data(ounit)


  end subroutine print_hc_diag_data


  subroutine print_state_energy_diag_data(ounit)
    use hc_init_lib_data ! hc interaction and characteristic data
    implicit none
    integer :: ounit
    !
    ! Local variables
    !
    integer hc_species
    ! 
    ! Print out the state energy diagnostic data
    !

    !
    ! Calculate state average data
    !
    call calculate_state_energy_averages

    !
    ! Print out the data
    ! 

    Write (Ounit,'(1X,A,47X,A)') "Species specific kinetic diagnostic data for all launches:","Species Produced"
    Write (Ounit,'(42X,A9,2X,A4,2X,100A11)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,HC_Species &
                        &= 1,Num_states,1)
    Write (Ounit,62) "Number of paricles enetering state","#",&
         & (total_entered_state(hc_species), HC_Species = 1, Num_states, 1)
    Write (Ounit,62) "Average kinetic energy of particles reaching state","eV", &
         & (ave_kin_energy(HC_Species), HC_Species = 1, Num_states, 1)
    Write (Ounit,62) "Average energy released reaching state","eV",&
         & (ave_energy_released(HC_species), HC_Species = 1, Num_states, 1)
    Write (Ounit,62) "Average total kinetic energy gained reaching state","eV",&
         & (ave_tot_energy_gained(Hc_species), HC_Species = 1, Num_states, 1)
    Write (Ounit,62) "Average number of transitions to reach state","#",&
         & (ave_transitions(Hc_species), HC_Species = 1, Num_states, 1)
    Write (Ounit,*) ""


    Write (6,'(1X,A,47X,A)') "Species specific kinetic diagnostic data for all launches:","Species Produced"
    Write (6,65) "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,HC_Species &
                        &= 1,Num_states,1)
    Write (6,64) "Number of paricles enetering state","#",&
         & (total_entered_state(hc_species), HC_Species = 1, Num_states, 1)
    Write (6,64) "Average kinetic energy of particles reaching state","eV", &
         & (ave_kin_energy(HC_Species), HC_Species = 1, Num_states, 1)
    Write (6,64) "Average energy released reaching state","eV",&
         & (ave_energy_released(HC_species), HC_Species = 1, Num_states, 1)
    Write (6,64) "Average total kinetic energy gained reaching state","eV",&
         & (ave_tot_energy_gained(Hc_species), HC_Species = 1, Num_states, 1)
    Write (6,64) "Average number of transitions to reach state","#",&
         & (ave_transitions(Hc_species), HC_Species = 1, Num_states, 1)
    Write (6,*) ""

    

60     Format (1X,A50,2X,A4,2X,100F11.8) ! Good for fractions/times up to 1.00000000
605    Format (1X,A50,2X,A4,2X,100F11.7) ! Good for fractions/times up to 10.0000000
61     Format (1X,A50,2X,A4,2X,100F11.6) ! Good for locations up to 9999.999999
62     Format (1X,A50,2X,A4,2X,100F11.3) ! Good for counts up to 99,999.99
63     Format (1X,A50,2X,A4,2X,100F11.1) ! Good for counts up to 999,999.9

64     format (1x,a50,2x,a4,2x,100(1x,g12.5))
65     format (42X,A9,2X,A4,2X,100(1x,a12))

  end subroutine print_state_energy_diag_data


  subroutine calculate_state_energy_averages
    implicit none
    integer :: in

    do in = 1,num_states
       if (total_entered_state(in).gt.0.0) then
          ave_kin_energy(in) = ave_kin_energy(in)/total_entered_state(in)
          ave_energy_released(in) = ave_energy_released(in)/total_entered_state(in)
          ave_tot_energy_gained(in) = ave_tot_energy_gained(in)/total_entered_state(in)
          ave_transitions(in) = ave_transitions(in)/total_entered_state(in)
       else
          ave_kin_energy(in) = 0.0
          ave_energy_released(in) = 0.0
          ave_tot_energy_gained(in) = 0.0
          ave_transitions(in) = 0.0
       endif
    end do

  end subroutine calculate_state_energy_averages






subroutine allocate_diag_data
  use allocate_arrays
  implicit none
  integer :: flag

  ! jdemod  - replace direct calls to allocate with the allocate arrays module that handles both error checking and initialization
  
    !allocate(total_entered_state(num_states),stat=flag)
    call allocate_array(total_entered_state,num_states,'total_entered_state',flag)

    !allocate(ave_kin_energy(num_states),stat=flag)
    call allocate_array(ave_kin_energy,num_states,'ave_kin_energy',flag)

    !allocate(ave_energy_released(num_states),stat=flag)
    call allocate_array(ave_energy_released,num_states,'ave_energy_released',flag)

    !allocate(ave_tot_energy_gained(num_states),stat=flag)
    call allocate_array(ave_tot_energy_gained,num_states,'ave_tot_energy_gained',flag)

    !allocate(ave_transitions(num_states),stat=flag)
    call allocate_array(ave_transitions,num_states,'ave_transitions',flag)

    !      total_entered_state = 0.0
    !      ave_kin_energy = 0.0
    !      ave_energy_released = 0.0
    !      ave_tot_energy_gained = 0.0
    !      ! jdemod - initialization forgotten for this one
    !      ave_transitions = 0.0
          
end subroutine allocate_diag_data

subroutine deallocate_diag_data
  implicit none

    deallocate(total_entered_state)
    deallocate(ave_kin_energy)
    deallocate(ave_energy_released)
    deallocate(ave_tot_energy_gained)
    deallocate(ave_transitions)

end subroutine deallocate_diag_data




end module hc_diag_data
