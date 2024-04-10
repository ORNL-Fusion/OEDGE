subroutine compute_ions_core_flux()
  use all_variables, only : zones, global_parameters,&
       reference_parameters, element_variables
  implicit none
  integer*4 :: j,k,n
  integer*4 :: Nx,Nz
  integer*4 :: n_element,n_charge
  do n=1,global_parameters%N_ions
     n_element=global_parameters%ions_list(n,1)
     n_charge=global_parameters%ions_list(n,2)
     element_variables(n_element)%core_outflux(n_charge)=0.
     do k=1,global_parameters%N_zones
        Nx=zones(k)%mesh%Nx
        Nz=zones(k)%mesh%Nz
        if(zones(k)%Neighbors(1).eq.-1) then
           do j=1,Nz
              element_variables(n_element)%core_outflux(n_charge)=&
                   element_variables(n_element)%core_outflux(n_charge)&
                   +zones(k)%species(n)%fluxes%fluxn(Nx,j,1)*&
                   zones(k)%metric_coefficients%dS_North_PU(Nx,j)
           end do
        end if
        if(zones(k)%Neighbors(2).eq.-1) then
           do j=1,Nz
              element_variables(n_element)%core_outflux(n_charge)=&
                   element_variables(n_element)%core_outflux(n_charge)&
                   -zones(k)%species(n)%fluxes%fluxn(1,j,2)*&
                   zones(k)%metric_coefficients%dS_South_PU(1,j)
           end do
        end if
     end do
  end do
  ! For the n_charge = 0 ==> neutral flux from Eirene
  do n_element=1,global_parameters%N_species
     element_variables(n_element)%total_flux_core_nonionized=0.
     !let us sum over the charge state between 0 and Z-1
     do n_charge=0,global_parameters%element_list(n_element)%Z-1
        element_variables(n_element)%total_flux_core_nonionized=&
             element_variables(n_element)%total_flux_core_nonionized+&
             element_variables(n_element)%core_outflux(n_charge)
     end do
  end do
end subroutine compute_ions_core_flux
