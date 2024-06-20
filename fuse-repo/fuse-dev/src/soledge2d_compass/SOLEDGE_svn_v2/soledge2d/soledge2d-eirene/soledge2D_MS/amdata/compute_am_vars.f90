subroutine compute_am_vars(zone)
  use all_variables, only : global_parameters,reference_parameters
  use MZone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  real*8,allocatable :: pow_logNe(:,:,:)
  real*8,allocatable :: pow_logTe(:,:,:)
  integer*4,parameter :: max_degree=8
  integer*4 :: Nx,Nz
  integer*4 :: degree
  integer*4 :: n,m,l,ind,k
  integer*4 :: element_index,charge
  real*8 :: Ne_min,Ne_max,Te_min,Te_max
  integer*4 :: i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(pow_logNe(0:max_degree-1,1:Nx,1:Nz))
  allocate(pow_logTe(0:max_degree-1,1:Nx,1:Nz))
  do n=1,global_parameters%N_ions
     element_index=zone%species(n)%element_index
     charge=zone%species(n)%charge
     degree=global_parameters%element_list(element_index)%amdata_polynom_degree

     !reset
     zone%species(n)%am_vars%ionization_rate_coefficient=0.D0
     zone%species(n)%am_vars%recombination_rate_coefficient=0.D0
     zone%species(n)%am_vars%radiation_function_excitation=0.D0
     zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung=0.D0

     !ionization
     if(zone%species(n)%compute_ionization) then
        Te_min=global_parameters%element_list(element_index)%amdatas(charge+1)%ionization_rate_polynom%Te_min
        Te_max=global_parameters%element_list(element_index)%amdatas(charge+1)%ionization_rate_polynom%Te_max
        Ne_min=global_parameters%element_list(element_index)%amdatas(charge+1)%ionization_rate_polynom%Ne_min
        Ne_max=global_parameters%element_list(element_index)%amdatas(charge+1)%ionization_rate_polynom%Ne_max
        call compute_power_logs(zone,pow_logNe,pow_logTe,max_degree,Nx,Nz,Te_min,Te_max,Ne_min,Ne_max)
        ind=0
        do m=0,degree-1
           do l=0,degree-1-m
              ind=ind+1
              zone%species(n)%am_vars%ionization_rate_coefficient=zone%species(n)%am_vars%ionization_rate_coefficient&
                   +global_parameters%element_list(element_index)%amdatas(charge+1)%ionization_rate_polynom%coefficients(ind)&
                   *pow_logNe(m,:,:)*pow_logTe(l,:,:)
           end do
        end do
        zone%species(n)%am_vars%ionization_rate_coefficient=exp(log(10.D0)*&
             zone%species(n)%am_vars%ionization_rate_coefficient)
     else
        zone%species(n)%am_vars%ionization_rate_coefficient=0.D0 ! no ionization
     end if
     !de-dimensionalization
     zone%species(n)%am_vars%ionization_rate_coefficient=zone%species(n)%am_vars%ionization_rate_coefficient&
          *(reference_parameters%fields%n0*reference_parameters%fields%tau0)

     !recombination
     if(zone%species(n)%compute_recombination) then
        Te_min=global_parameters%element_list(element_index)%amdatas(charge)%recombination_rate_polynom%Te_min
        Te_max=global_parameters%element_list(element_index)%amdatas(charge)%recombination_rate_polynom%Te_max
        Ne_min=global_parameters%element_list(element_index)%amdatas(charge)%recombination_rate_polynom%Ne_min
        Ne_max=global_parameters%element_list(element_index)%amdatas(charge)%recombination_rate_polynom%Ne_max
        call compute_power_logs(zone,pow_logNe,pow_logTe,max_degree,Nx,Nz,Te_min,Te_max,Ne_min,Ne_max)
        ind=0
        do m=0,degree-1
           do l=0,degree-1-m
              ind=ind+1
              zone%species(n)%am_vars%recombination_rate_coefficient=zone%species(n)%am_vars%recombination_rate_coefficient&
                   +global_parameters%element_list(element_index)%amdatas(charge)%recombination_rate_polynom%coefficients(ind)&
                   *pow_logNe(m,:,:)*pow_logTe(l,:,:)
           end do
        end do
        zone%species(n)%am_vars%recombination_rate_coefficient=exp(log(10.D0)*&
             zone%species(n)%am_vars%recombination_rate_coefficient)
     else
        zone%species(n)%am_vars%recombination_rate_coefficient=0.D0 ! recombination to H0 computed with eirene
     end if
     !de-dimensionalization
     zone%species(n)%am_vars%recombination_rate_coefficient=zone%species(n)%am_vars%recombination_rate_coefficient&
          *(reference_parameters%fields%n0*reference_parameters%fields%tau0)

     !line excitation
     if(zone%species(n)%compute_ionization) then
        Te_min=global_parameters%element_list(element_index)%amdatas(charge+1)%line_excitation_polynom%Te_min
        Te_max=global_parameters%element_list(element_index)%amdatas(charge+1)%line_excitation_polynom%Te_max
        Ne_min=global_parameters%element_list(element_index)%amdatas(charge+1)%line_excitation_polynom%Ne_min
        Ne_max=global_parameters%element_list(element_index)%amdatas(charge+1)%line_excitation_polynom%Ne_max
        call compute_power_logs(zone,pow_logNe,pow_logTe,max_degree,Nx,Nz,Te_min,Te_max,Ne_min,Ne_max)
        ind=0
        do m=0,degree-1
           do l=0,degree-1-m
              ind=ind+1
              zone%species(n)%am_vars%radiation_function_excitation=zone%species(n)%am_vars%radiation_function_excitation&
                   +global_parameters%element_list(element_index)%amdatas(charge+1)%line_excitation_polynom%coefficients(ind)&
                   *pow_logNe(m,:,:)*pow_logTe(l,:,:)
           end do
        end do
        zone%species(n)%am_vars%radiation_function_excitation=exp(log(10.D0)*&
             zone%species(n)%am_vars%radiation_function_excitation)
     else
        zone%species(n)%am_vars%radiation_function_excitation=0.D0 ! recombination to H0 computed with eirene
     end if
     !de-dimensionalization
     zone%species(n)%am_vars%radiation_function_excitation=zone%species(n)%am_vars%radiation_function_excitation&
          *(reference_parameters%fields%n0*reference_parameters%fields%tau0/(kb*reference_parameters%fields%T0))
     

     !line recombination + BML
     if(zone%species(n)%compute_recombination) then
        Te_min=global_parameters%element_list(element_index)%amdatas(charge)%line_recombination_polynom%Te_min
        Te_max=global_parameters%element_list(element_index)%amdatas(charge)%line_recombination_polynom%Te_max
        Ne_min=global_parameters%element_list(element_index)%amdatas(charge)%line_recombination_polynom%Ne_min
        Ne_max=global_parameters%element_list(element_index)%amdatas(charge)%line_recombination_polynom%Ne_max
        call compute_power_logs(zone,pow_logNe,pow_logTe,max_degree,Nx,Nz,Te_min,Te_max,Ne_min,Ne_max)
        ind=0
        do m=0,degree-1
           do l=0,degree-1-m
              ind=ind+1
              zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung=zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung&
                   +global_parameters%element_list(element_index)%amdatas(charge)%line_recombination_polynom%coefficients(ind)&
                   *pow_logNe(m,:,:)*pow_logTe(l,:,:)
           end do
        end do
        zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung=exp(log(10.D0)*&
             zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung)
     else
        zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung=0.D0 ! recombination to H0 computed with eirene
     end if
     !de-dimensionalization
     zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung=zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung&
          *(reference_parameters%fields%n0*reference_parameters%fields%tau0/(kb*reference_parameters%fields%T0))


  end do
  deallocate(pow_logNe,pow_logTe)
end subroutine compute_am_vars
