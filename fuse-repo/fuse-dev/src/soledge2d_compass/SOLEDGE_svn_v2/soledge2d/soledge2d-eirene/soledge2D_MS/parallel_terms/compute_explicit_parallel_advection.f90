subroutine compute_explicit_parallel_advection(zone)
#include "compile_opt.inc"
  use all_variables, only : global_parameters, flags, drift_flags
  use Mzone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  real*8,allocatable :: Fluxes_n(:,:,:)
  real*8,allocatable :: Fluxes_G(:,:,:)
  real*8,allocatable :: Fluxes_E(:,:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  integer*4 :: i,j
  character(50) :: filename
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes_n(1:Nx,1:Nz,1:4))
  allocate(Fluxes_G(1:Nx,1:Nz,1:4))
  allocate(Fluxes_E(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  do n=1,global_parameters%n_ions
     call compute_source_para(zone,Fluxes_n,Fluxes_G,Fluxes_E,Nx,Nz,n)
     !density
     zone%species(n)%fluxes%fluxn=zone%species(n)%fluxes%fluxn+Fluxes_n
!!$     !$OMP CRITICAL
!!$     if(n.eq.1) then
!!$        write(filename,"(A6,I0)") "FluxE_",zone%number
!!$        open(unit=10,file=trim(filename),status='unknown')
!!$        do i=1,Nx
!!$           write(10,100) (Fluxes_n(i,j,3),j=1,Nz)
!!$        end do
!!$        close(10)
!!$        write(filename,"(A6,I0)") "FluxW_",zone%number
!!$        open(unit=10,file=trim(filename),status='unknown')
!!$        do i=1,Nx
!!$           write(10,100) (Fluxes_n(i,j,4),j=1,Nz)
!!$        end do
!!$        close(10)
!!$     end if
!!$100  format(512es15.7)
!!$     !$OMP END CRITICAL
    Source=divergence(zone,Fluxes_n,Nx,Nz)
     zone%species(n)%sources%Sn=zone%species(n)%sources%Sn-Source
     !velocity
     zone%species(n)%fluxes%fluxG=zone%species(n)%fluxes%fluxG+Fluxes_G
     Source=divergence(zone,Fluxes_G,Nx,Nz)
     zone%species(n)%sources%SG=zone%species(n)%sources%SG-Source
     !temperature
     zone%species(n)%fluxes%fluxE=zone%species(n)%fluxes%fluxE+Fluxes_E
     Source=divergence(zone,Fluxes_E,Nx,Nz)
!!$     !$OMP CRITICAL
!!$     if(n.eq.1) then
!!$        write(filename,"(A5,I0)") "Spol_",zone%number
!!$        open(unit=10,file=trim(filename),status='unknown')
!!$        do i=1,Nx
!!$           write(10,100) (Source(i,j),j=1,Nz)
!!$        end do
!!$        close(10)
!!$     end if
!!$     !$OMP END CRITICAL
     zone%species(n)%sources%SE=zone%species(n)%sources%SE-Source
  end do
  ! electrons
  call compute_source_para_e(zone,Fluxes_E,Nx,Nz)
  !temperature
  zone%species(0)%fluxes%fluxE=zone%species(0)%fluxes%fluxE+Fluxes_E
  Source=divergence(zone,Fluxes_E,Nx,Nz)
  zone%species(0)%sources%SE=zone%species(0)%sources%SE-Source
#if VORTICITY_PASTIX == 1
  ! vorticity
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  call compute_source_para_passive_scalar(zone,zone%electric_fields(1)%vorticity,Fluxes,Nx,Nz)
  !  call set_perp_flux_to_zero_on_the_wall(zone,Fluxes_W,Nx,Nz)
  if(drift_flags%jadvW) then
     zone%electric_fields(1)%j_para_adv_W=-Fluxes
  else
     zone%electric_fields(1)%j_para_adv_W=0.d0
  end if
  deallocate(Fluxes)
#endif
  if(flags%turbulence_model.eq.1) then
     ! k
     allocate(Fluxes(1:Nx,1:Nz,1:4))
     Fluxes=0.
     call compute_source_para_passive_scalar(zone,zone%kepsilon(1)%k*zone%species(0)%var(1)%density,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%kepsilon(1)%Sk=zone%kepsilon(1)%Sk-Source
     ! epsilon
     Fluxes=0.
     call compute_source_para_passive_scalar(zone,zone%kepsilon(1)%epsilon*zone%species(0)%var(1)%density,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%kepsilon(1)%Sepsilon=zone%kepsilon(1)%Sepsilon-Source
     deallocate(Fluxes)
  end if
  deallocate(Fluxes_n,Fluxes_G,Fluxes_E)
  deallocate(Source)
end subroutine compute_explicit_parallel_advection
