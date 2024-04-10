subroutine compute_am_sources(zone)
  use all_variables, only : global_parameters, global_variables, reference_parameters
  use MZone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: n,i,j
  integer*4 :: Nx,Nz
  real*8,allocatable :: Source_r(:,:),Source_i(:,:),Source(:,:)
  real*8 :: max_rate
  integer*4 :: element_index,charge
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  max_rate=1.d0/global_variables%dt*0.5d0
  allocate(Source_r(1:Nx,1:Nz),Source_i(1:Nx,1:Nz),Source(1:Nx,1:Nz))
  do n=1,global_parameters%N_ions
     Source_r=0.D0
     Source_i=0.D0
     Source=0.D0
     if(zone%species(n)%compute_recombination) then
        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%recombination_rate_coefficient,max_rate)
        Source_r=Source_r+zone%species(n-1)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n-1)%am_vars%ionization_rate_coefficient,max_rate)
!!$        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%recombination_rate_coefficient
!!$        Source_r=Source_r+zone%species(n-1)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n-1)%am_vars%ionization_rate_coefficient
     end if
     if(zone%species(n)%compute_ionization) then
        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%ionization_rate_coefficient,max_rate)
        Source_i=Source_i+zone%species(n+1)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n+1)%am_vars%recombination_rate_coefficient,max_rate)
!!$        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%ionization_rate_coefficient
!!$        Source_i=Source_i+zone%species(n+1)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n+1)%am_vars%recombination_rate_coefficient
     end if
     Source=Source_r+Source_i
     do i=1,Nx
        do j=1,Nz
           if(zone%masks%chi2(i,j).eq.1) then
              Source(i,j)=0.d0
           end if
        end do
     end do
     zone%species(n)%sources%Sn_am=Source

     Source_r=0.D0
     Source_i=0.D0
     Source=0.D0
     if(zone%species(n)%compute_recombination) then
        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%recombination_rate_coefficient,max_rate)&
             *zone%species(n)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n)%element%mass
        Source_r=Source_r+zone%species(n-1)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n-1)%am_vars%ionization_rate_coefficient,max_rate)&
             *zone%species(n-1)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n-1)%element%mass
!!$        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%recombination_rate_coefficient&
!!$             *zone%species(n)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n)%element%mass
!!$        Source_r=Source_r+zone%species(n-1)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n-1)%am_vars%ionization_rate_coefficient&
!!$             *zone%species(n-1)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n-1)%element%mass
     end if
     if(zone%species(n)%compute_ionization) then
        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%ionization_rate_coefficient,max_rate)&
             *zone%species(n)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n)%element%mass
        Source_i=Source_i+zone%species(n+1)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n+1)%am_vars%recombination_rate_coefficient,max_rate)&
             *zone%species(n+1)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n+1)%element%mass
!!$        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%ionization_rate_coefficient&
!!$             *zone%species(n)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n)%element%mass
!!$        Source_i=Source_i+zone%species(n+1)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n+1)%am_vars%recombination_rate_coefficient&
!!$             *zone%species(n+1)%var(1)%velocity(1:Nx,1:Nz)*zone%species(n+1)%element%mass
     end if
     Source=Source_r+Source_i
     do i=1,Nx
        do j=1,Nz
           if(zone%masks%chi2(i,j).eq.1) then
              Source(i,j)=0.d0
           end if
        end do
     end do
     zone%species(n)%sources%SG_am=Source/zone%species(n)%element%mass

     Source_r=0.D0
     Source_i=0.D0
     Source=0.D0
     if(zone%species(n)%compute_recombination) then
        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%recombination_rate_coefficient,max_rate)&
             *(0.5d0*zone%species(n)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n)%element%mass&
             +1.5d0*zone%species(n)%var(1)%temperature(1:Nx,1:Nz))
        Source_r=Source_r+zone%species(n-1)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n-1)%am_vars%ionization_rate_coefficient,max_rate)&
             *(0.5d0*zone%species(n-1)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n-1)%element%mass&
             +1.5d0*zone%species(n-1)%var(1)%temperature(1:Nx,1:Nz))
!!$        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%recombination_rate_coefficient&
!!$             *(0.5d0*zone%species(n)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n)%element%mass&
!!$             +1.5d0*zone%species(n)%var(1)%temperature(1:Nx,1:Nz))
!!$        Source_r=Source_r+zone%species(n-1)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n-1)%am_vars%ionization_rate_coefficient&
!!$             *(0.5d0*zone%species(n-1)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n-1)%element%mass&
!!$             +1.5d0*zone%species(n-1)%var(1)%temperature(1:Nx,1:Nz))
     end if
     if(zone%species(n)%compute_ionization) then
        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%ionization_rate_coefficient,max_rate)&
             *(0.5d0*zone%species(n)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n)%element%mass&
             +1.5d0*zone%species(n)%var(1)%temperature(1:Nx,1:Nz))
        Source_i=Source_i+zone%species(n+1)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n+1)%am_vars%recombination_rate_coefficient,max_rate)&
             *(0.5d0*zone%species(n+1)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n+1)%element%mass&
             +1.5d0*zone%species(n+1)%var(1)%temperature(1:Nx,1:Nz))
!!$        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%ionization_rate_coefficient&
!!$             *(0.5d0*zone%species(n)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n)%element%mass&
!!$             +1.5d0*zone%species(n)%var(1)%temperature(1:Nx,1:Nz))
!!$        Source_i=Source_i+zone%species(n+1)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n+1)%am_vars%recombination_rate_coefficient&
!!$             *(0.5d0*zone%species(n+1)%var(1)%velocity(1:Nx,1:Nz)**2*zone%species(n+1)%element%mass&
!!$             +1.5d0*zone%species(n+1)%var(1)%temperature(1:Nx,1:Nz))
     end if
     Source=Source_r+Source_i
     do i=1,Nx
        do j=1,Nz
           if(zone%masks%chi2(i,j).eq.1) then
              Source(i,j)=0.d0
           end if
        end do
     end do
     zone%species(n)%sources%SE_am=Source
  end do

  !electrons
  zone%species(0)%sources%rad=0.D0
  Source_r=0.D0
  Source_i=0.D0
  Source=0.D0
  do n=1,global_parameters%N_ions
     zone%species(n)%sources%rad=0.D0
     element_index=zone%species(n)%element_index
     charge=zone%species(n)%charge
     if(zone%species(n)%compute_recombination) then
        Source_r=Source_r+zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%recombination_rate_coefficient,max_rate)*&
             global_parameters%element_list(element_index)%amdatas(charge)%Ionization_potential/reference_parameters%fields%T0eV
        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung
        zone%species(n)%sources%rad=zone%species(n)%sources%rad&
             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung
        zone%species(0)%sources%rad=zone%species(0)%sources%rad&
             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung
!!$        Source_r=Source_r+zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%recombination_rate_coefficient*&
!!$             global_parameters%element_list(element_index)%amdatas(charge)%Ionization_potential/reference_parameters%fields%T0eV
!!$        Source_r=Source_r-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung
!!$        zone%species(n)%sources%rad=zone%species(n)%sources%rad&
!!$             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung
!!$        zone%species(0)%sources%rad=zone%species(0)%sources%rad&
!!$             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%radiation_function_recombination_bremsstrahlung
     end if
     if(zone%species(n)%compute_ionization) then
        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*min(zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%ionization_rate_coefficient,max_rate)*&
             global_parameters%element_list(element_index)%amdatas(charge+1)%Ionization_potential/reference_parameters%fields%T0eV
        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%radiation_function_excitation
        zone%species(n)%sources%rad=zone%species(n)%sources%rad&
             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%radiation_function_excitation
        zone%species(0)%sources%rad=zone%species(0)%sources%rad&
             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
             *zone%species(n)%am_vars%radiation_function_excitation
!!$        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%ionization_rate_coefficient*&
!!$             global_parameters%element_list(element_index)%amdatas(charge+1)%Ionization_potential/reference_parameters%fields%T0eV
!!$        Source_i=Source_i-zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%radiation_function_excitation
!!$        zone%species(n)%sources%rad=zone%species(n)%sources%rad&
!!$             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%radiation_function_excitation
!!$        zone%species(0)%sources%rad=zone%species(0)%sources%rad&
!!$             +zone%species(n)%var(1)%density(1:Nx,1:Nz)*zone%species(0)%var(1)%density(1:Nx,1:Nz)&
!!$             *zone%species(n)%am_vars%radiation_function_excitation
     end if
  end do
  Source=Source_r+Source_i
  do i=1,Nx
     do j=1,Nz
        if(zone%masks%chi2(i,j).eq.1) then
           Source(i,j)=0.d0
        end if
     end do
  end do
  zone%species(0)%sources%SE_am=Source

end subroutine compute_am_sources
