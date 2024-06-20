subroutine compute_pdivv_term(zone)
  use all_variables, only : reference_parameters, global_variables, flags, global_parameters
  use Mzone
  use Mphysics
  use Moperator
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: Nx,Nz,n_ion
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  integer*4 :: i,j,n
  real*8 :: R0,rs0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  do n=0,global_parameters%N_ions
     Fluxes=0.D0
     do i=1,Nx
        do j=1,Nz
!!$           Fluxes(i,j,3)=0.5D0*(zone%species(n)%var(1)%density(i,j)*zone%species(n)%var(1)%temperature(i,j)&
!!$                +zone%species(n)%var(1)%density(i,j+1)*zone%species(n)%var(1)%temperature(i,j+1))
!!$           Fluxes(i,j,3)=Fluxes(i,j,3)*0.5d0*(zone%species(n)%drifts%uEt(i,j)+zone%species(n)%drifts%uBt(i,j)+&
!!$                zone%species(n)%drifts%uEt(i,j+1)+zone%species(n)%drifts%uBt(i,j+1))*(2.d0*pi*R0/rs0)
!!$           Fluxes(i,j,4)=0.5D0*(zone%species(n)%var(1)%density(i,j)*zone%species(n)%var(1)%temperature(i,j)&
!!$                +zone%species(n)%var(1)%density(i,j-1)*zone%species(n)%var(1)%temperature(i,j-1))
!!$           Fluxes(i,j,4)=Fluxes(i,j,4)*0.5d0*(zone%species(n)%drifts%uEt(i,j)+zone%species(n)%drifts%uBt(i,j)+&
!!$                zone%species(n)%drifts%uEt(i,j-1)+zone%species(n)%drifts%uBt(i,j-1))*(2.d0*pi*R0/rs0)    
           Fluxes(i,j,3)=0.5d0*(zone%species(n)%drifts%uEt(i,j)+zone%species(n)%drifts%uBt(i,j)+&
                zone%species(n)%drifts%uEt(i,j+1)+zone%species(n)%drifts%uBt(i,j+1))*(2.d0*pi*R0/rs0)
           Fluxes(i,j,4)=0.5d0*(zone%species(n)%drifts%uEt(i,j)+zone%species(n)%drifts%uBt(i,j)+&
                zone%species(n)%drifts%uEt(i,j-1)+zone%species(n)%drifts%uBt(i,j-1))*(2.d0*pi*R0/rs0)    
           Fluxes(i,j,1)=0.5d0*(zone%species(n)%drifts%uEp(i,j)+zone%species(n)%drifts%uBp(i,j)+&
                zone%species(n)%drifts%uEp(i+1,j)+zone%species(n)%drifts%uBp(i+1,j))*(2.d0*pi*R0/rs0)
           Fluxes(i,j,2)=0.5d0*(zone%species(n)%drifts%uEp(i,j)+zone%species(n)%drifts%uBp(i,j)+&
                zone%species(n)%drifts%uEp(i-1,j)+zone%species(n)%drifts%uBp(i-1,j))*(2.d0*pi*R0/rs0)    
        end do
     end do
!     zone%species(n)%fluxes%fluxE=zone%species(n)%fluxes%fluxE+Fluxes
     Source=divergence(zone,Fluxes,Nx,Nz)*zone%species(n)%var(1)%density(i,j)*zone%species(n)%var(1)%temperature(i,j)
     zone%species(n)%sources%SE=zone%species(n)%sources%SE-Source
  end do
end subroutine compute_pdivv_term
