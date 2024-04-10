subroutine affect_element_and_charge_to_species()
  use all_variables, only : global_parameters, element_variables
  implicit none
  integer*4 :: n,n_element,n_charge,Z,n_ion
  allocate(global_parameters%ions_list(1:global_parameters%N_ions,1:2))
  n=1
  do n_element=1,global_parameters%N_species
     do n_charge=1,global_parameters%element_list(n_element)%Z
        global_parameters%ions_list(n,1)=n_element
        global_parameters%ions_list(n,2)=n_charge
        n=n+1
     end do
     allocate(element_variables(n_element)%core_outflux(0:global_parameters%element_list(n_element)%Z))
  end do
  allocate(global_parameters%ind_ion_0(1:global_parameters%N_species))
  n=1
  do n_ion=1,global_parameters%N_ions
     if(global_parameters%ions_list(n_ion,2).eq.1) then
        global_parameters%ind_ion_0(n)=n_ion
        n=n+1
     end if
  end do
end subroutine affect_element_and_charge_to_species
