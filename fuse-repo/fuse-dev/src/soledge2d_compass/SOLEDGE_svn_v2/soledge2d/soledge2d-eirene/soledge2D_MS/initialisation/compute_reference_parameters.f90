subroutine compute_reference_parameters()
  use all_variables, only : reference_parameters, boundary_conditions, flags
  implicit none
  call compute_reference_density()
  call compute_reference_temperature()
  call compute_reference_electric_fields()
  call compute_reference_velocity_and_time()
  call compute_reference_turbulence()
  call save_references()
end subroutine compute_reference_parameters

subroutine compute_reference_geometry()
  use all_variables, only : reference_parameters, flags
  implicit none
  if(.not.flags%is_SLAB) then
     ! in the "slab" case, the reference geometry is set by the slab input file
     ! in the "non-slab" case, the following arbitrary values are used
     reference_parameters%geometry%R0=2.
     reference_parameters%geometry%Rm0=0.5
     reference_parameters%geometry%rs0=0.1
     reference_parameters%geometry%Bpol0=0.2
     reference_parameters%geometry%Btor0=2.
  end if
end subroutine compute_reference_geometry

subroutine compute_reference_density()
  use all_variables, only : reference_parameters, boundary_conditions, flags, transport_parameters
  implicit none
  if(reference_parameters%fields%n0.eq.0.) then 
! n0 has not been set in the input file so it is computed from BC
     if(boundary_conditions%BCn_model(1).eq.0) then
        ! density value BC
        reference_parameters%fields%n0=boundary_conditions%BCn(1)
     else
        ! flux boundary condition
        reference_parameters%fields%n0=max(boundary_conditions%BCn(1)*reference_parameters%geometry%rs0&
             /transport_parameters%Dn_p(1),1.d19)
     end if
  end if
end subroutine compute_reference_density


subroutine compute_reference_temperature()
  use all_variables, only : reference_parameters, boundary_conditions, flags, transport_parameters, pi
  use Mphysics
  implicit none
  real*8 :: electron_core_energy_flux
  if(reference_parameters%fields%T0eV.eq.0.) then 
! T0eV has not been set in the input file so it is computed from BC
     if(boundary_conditions%BCT_model(1).eq.0) then
        ! temperature value BC
        reference_parameters%fields%T0eV=boundary_conditions%BCTe
        reference_parameters%fields%T0=reference_parameters%fields%T0eV*eV/kb
     else
        ! flux boundary condition
        electron_core_energy_flux=boundary_conditions%BCTe/reference_parameters%geometry%Score
        reference_parameters%fields%T0=electron_core_energy_flux*reference_parameters%geometry%rs0&
             /(1.5d0*kb*reference_parameters%fields%n0*transport_parameters%chie_p)
        reference_parameters%fields%T0eV=reference_parameters%fields%T0*kb/eV
     end if
  else
     reference_parameters%fields%T0=reference_parameters%fields%T0eV*eV/kb
  end if
end subroutine compute_reference_temperature


subroutine compute_reference_velocity_and_time()
  use all_variables, only : reference_parameters, boundary_conditions, flags, global_parameters
  use Mphysics
  implicit none
  reference_parameters%fields%c0=sqrt(kb*reference_parameters%fields%T0&
       /(m_u))
  reference_parameters%fields%tau0=2.D0*pi*reference_parameters%geometry%R0&
       /reference_parameters%fields%c0
end subroutine compute_reference_velocity_and_time

subroutine compute_reference_turbulence()
  use all_variables, only : reference_parameters, boundary_conditions, flags, global_parameters
  use Mphysics
  implicit none
  reference_parameters%fields%k0=1./(reference_parameters%fields%tau0*0.09)
  reference_parameters%fields%epsilon0=reference_parameters%fields%k0/reference_parameters%fields%tau0
  ! 1 for 1m2/s
  ! 0.09 for Cmu=0.09
end subroutine compute_reference_turbulence

subroutine compute_reference_electric_fields()
  use all_variables, only : reference_parameters, boundary_conditions, flags, global_parameters
  use Mphysics
  implicit none
  reference_parameters%fields%phi0=reference_parameters%fields%T0eV
  reference_parameters%fields%W0=(reference_parameters%fields%phi0*reference_parameters%fields%n0)&
       /(reference_parameters%geometry%rs0**2)*m_u
end subroutine compute_reference_electric_fields
