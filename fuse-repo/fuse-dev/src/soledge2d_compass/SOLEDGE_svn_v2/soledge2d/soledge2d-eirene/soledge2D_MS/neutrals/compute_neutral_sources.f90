subroutine compute_neutral_sources(zone)
  use all_variables, only : global_parameters, reference_parameters
  use MZone
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: ne,nn,Te,Gammai,Ti,mi
  real*8 :: compute_FN_ionization_source
  real*8 :: compute_FN_momentum_losses
  real*8 :: compute_FN_E_electron_losses
  real*8 :: compute_FN_E_ion_losses
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        ne=zone%species(0)%var(1)%density(i,j)*reference_parameters%fields%n0
        nn=zone%neutrals%density(i,j)*reference_parameters%fields%n0
        Te=zone%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV
        Ti=zone%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV
        Gammai=zone%species(1)%var(1)%Gamma(i,j)*reference_parameters%fields%c0*reference_parameters%fields%n0
        mi=zone%species(1)%element%mass*m_u
        zone%species(1)%sources%Sn_n(i,j) = compute_FN_ionization_source(ne,nn,Te)
        zone%species(1)%sources%Sn_G(i,j) = compute_FN_momentum_losses(ne,nn,Te,Gammai)
        zone%species(0)%sources%Sn_E(i,j) = compute_FN_E_electron_losses(ne,nn,Te)
        zone%species(1)%sources%Sn_E(i,j) = compute_FN_E_ion_losses(ne,nn,Te,Gammai,Ti,mi)
     end do
  end do
end subroutine compute_neutral_sources


!#######################################################
!############ subfunctions #############################

function compute_FN_ionization_source(ne,nn,Te) result(Sn)
  use all_variables, only : reference_parameters
  real*8,intent(in) :: ne,nn,Te
  real*8 :: Sn
  real*8 :: sigmav, sigmav_r
  real*8 :: Tec
  Tec=max(Te,0.3d0)
  !from NRL
  sigmav=1d-11*(Tec/13.6d0)**(0.5d0)&
       /(13.6d0**1.5d0*(6.d0+Tec/13.6d0))&
       *exp(-13.6d0/(Tec))
  sigmav_r=5.2d-20*sqrt(13.6d0/(Tec))*&
       (0.43d0+0.5d0*log(13.6d0/(Tec))+&
       0.469d0*(13.6d0/(Tec))**(-1.d0/3.d0))
  sigmav_r=sigmav_r+8.75d-39*(Tec)**(-4.5d0)*ne
  Sn=nn*ne*(sigmav-sigmav_r)
  !DD
  Sn=Sn/(reference_parameters%fields%n0/reference_parameters%fields%tau0)
end function compute_FN_ionization_source

function compute_FN_momentum_losses(ne,nn,Te,Gammai) result(SG)
  use all_variables, only : reference_parameters
  real*8,intent(in) :: ne,nn,Te,Gammai
  real*8 :: SG
  real*8 :: sigmav, sigmav_r, sigmav_cx
  real*8 :: Tec
  Tec=max(Te,0.3d0)
  !from NRL
  sigmav=1d-11*(Tec/13.6d0)**(0.5d0)&
       /(13.6d0**1.5d0*(6.d0+Tec/13.6d0))&
       *exp(-13.6d0/(Tec))
  sigmav_r=5.2d-20*sqrt(13.6d0/(Tec))*&
       (0.43d0+0.5d0*log(13.6d0/(Tec))+&
       0.469d0*(13.6d0/(Tec))**(-1.d0/3.d0))
  sigmav_r=sigmav_r+8.75d-39*(Tec)**(-4.5d0)*ne
  sigmav_cx=exp(-0.5D0/Tec)*2.5d-15/exp(-0.5d0)
  SG=-nn*Gammai*sigmav_cx+nn*Gammai*(sigmav-sigmav_r)
  !DD
  SG=SG/(reference_parameters%fields%n0*reference_parameters%fields%c0&
       /reference_parameters%fields%tau0)
end function compute_FN_momentum_losses

function compute_FN_E_electron_losses(ne,nn,Te) result(SE)
  use all_variables, only : reference_parameters
  real*8,intent(in) :: ne,nn,Te
  real*8 :: SE
  real*8 :: sigmav, sigmav_r
  real*8 :: Tec,Tloss,Tloss_r
  Tec=max(Te,0.3d0)
  !from NRL
  sigmav=1d-11*(Tec/13.6d0)**(0.5d0)&
       /(13.6d0**1.5d0*(6.d0+Tec/13.6d0))&
       *exp(-13.6d0/(Tec))
  sigmav_r=5.2d-20*sqrt(13.6d0/(Tec))*&
       (0.43d0+0.5d0*log(13.6d0/(Tec))+&
       0.469d0*(13.6d0/(Tec))**(-1.d0/3.d0))
  sigmav_r=sigmav_r+8.75d-39*(Tec)**(-4.5d0)*ne
  Tloss=25.D0+170.D0*exp(-Tec/2.D0)
  Tloss_r=min(250.D0,8.D0*exp(Tec/9.D0))
  SE=-(nn*ne*sigmav*Tloss&
       +nn*ne*sigmav_r*Tloss_r)
  !DD
  SE=SE/(reference_parameters%fields%n0*reference_parameters%fields%T0eV/&
       reference_parameters%fields%tau0)
end function compute_FN_E_electron_losses

function compute_FN_E_ion_losses(ne,nn,Te,Gammai,Ti,mi) result(SE)
  use all_variables, only : reference_parameters
  use Mphysics
  real*8,intent(in) :: ne,nn,Te,Gammai,Ti,mi
  real*8 :: SE
  real*8 :: sigmav, sigmav_r, sigmav_cx
  real*8 :: Tec
  Tec=max(Te,0.3d0)
  !from NRL
  sigmav=1d-11*(Tec/13.6d0)**(0.5d0)&
       /(13.6d0**1.5d0*(6.d0+Tec/13.6d0))&
       *exp(-13.6d0/(Tec))
  sigmav_r=5.2d-20*sqrt(13.6d0/(Tec))*&
       (0.43d0+0.5d0*log(13.6d0/(Tec))+&
       0.469d0*(13.6d0/(Tec))**(-1.d0/3.d0))
  sigmav_r=sigmav_r+8.75d-39*(Tec)**(-4.5d0)*ne
  sigmav_cx=exp(-0.5D0/Tec)*2.5d-15/exp(-0.5d0)
  SE=-nn*Gammai**2.D0/ne*sigmav_cx*mi/eV&       
       +nn*ne*(sigmav-sigmav_r)*(0.5d0*(Gammai/ne)**2.d0*mi/eV+1.5d0*Ti)
  !DD
  SE=SE/(reference_parameters%fields%n0*reference_parameters%fields%T0eV/&
       reference_parameters%fields%tau0)
end function compute_FN_E_ion_losses
