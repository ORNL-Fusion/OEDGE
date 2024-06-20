subroutine compute_vorticity_RHS4(mode_BC,mode_BC_para)
#include "compile_opt.inc"
  use all_variables, only : global_parameters, zones, reference_parameters, drift_flags
  use Mvorticity_vars
  use Mlist
  use Mphysics
  use Mpastix_solve
  use test_var
  use Moperator
  use MDiffusion_perp
  implicit none
  integer*4,intent(in) :: mode_BC, mode_BC_para
  integer*4 i,j,k
  integer*4 Nx,Nz
  integer*4 East,West,North,South
  real*8 Lambda,Lambdap,Lambdam
  real*8 :: r1, r2, r3, alpha, beta
  real*8 :: R0,c0,rs0,Gamma
  real*8,parameter :: eps=1.d-6

  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  !if(mode_BC_para.eq.LAMBDA_TE) then
  !   write(*,*) 'Parallel BC for potential (RHS) : Lambda*Te'
  !else
  !   write(*,*) 'Parallel BC for potential (RHS) : Bohm current'
  !end if

  CSC_vort%b=0.D0
  do k=1, global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz

     do i=1,Nx
        do j=1,Nz
           if(zones(k)%masks%chi2(i,j).eq.0) then
              CSC_vort%b(Zones(k)%mesh%index(i,j))=zones(k)%electric_fields(1)%RHS(i,j)
           else
              if(zones(k)%masks%chi1(i,j).eq.1) then
                 Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+max(min(Zones(k)%species(1)%var(1)%temperature(i,j-1)&
                      /Zones(k)%species(0)%var(1)%temperature(i,j-1),10.D0),0.1D0)))
                 if(mode_BC_para.eq.LAMBDA_TE) then
                    CSC_vort%b(Zones(k)%mesh%index(i,j))=Lambda*Zones(k)%species(0)%var(1)%temperature(i,j-1)
                 else
                    r1=-(zones(k)%metric_coefficients%G(i,j)+zones(k)%metric_coefficients%G(i,j-1))*0.5D0&
                         /((zones(k)%species(0)%var(1)%density(i,j-1)+zones(k)%species(0)%var(1)%density(i,j))*0.5d0)&
                         *(zones(k)%species(0)%var(1)%temperature(i,j-1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,j)**(1.5d0))*0.5D0&
                         *(zones(k)%species(0)%var(1)%density(i,j)*zones(k)%species(0)%var(1)%temperature(i,j)&
                         -zones(k)%species(0)%var(1)%density(i,j-1)*zones(k)%species(0)%var(1)%temperature(i,j-1))&
                         /(zones(k)%mesh%z(i,j)-zones(k)%mesh%z(i,j-1))
                    alpha=reference_parameters%fields%c0*reference_parameters%fields%W0*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                         *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
                    beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                         *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
                    if(drift_flags%jadvW) then
                       r2=alpha*(zones(k)%species(0)%var(1)%velocity(i,j-1)+&
                            R0/c0*(zones(k)%species(0)%drifts%uEt(i,j-1)+zones(k)%species(0)%drifts%uBt(i,j-1))/(zones(k)%metric_coefficients%G(i,j-1)))&
                            *zones(k)%electric_fields(1)%vorticity(i,j-1)
                    else
                       r2=0.d0
                    end if
                    
                    !                    r3=beta*zones(k)%species(0)%var(1)%Gamma(i,j-1)*(1.D0-exp(Lambda-zones(k)%electric_fields(1)%phi(i,j-1)/zones(k)%species(0)%var(1)%temperature(i,j-1))&
                    !                         *(1.D0+zones(k)%electric_fields(0)%phi(i,j)/zones(k)%species(0)%var(1)%temperature(i,j-1)))
                    Gamma=zones(k)%species(0)%var(1)%Gamma(i,j-1)+(zones(k)%species(0)%drifts%uEt(i,j-1)+zones(k)%species(0)%drifts%uBt(i,j-1))&
                         *(2.d0*pi*R0/rs0)*sqrt(zones(k)%metric_coefficients%ctt(i,j-1))/zones(k)%metric_coefficients%G(i,j-1)&
                         *zones(k)%species(0)%var(1)%density(i,j-1)
                    r3=-beta*Gamma*Lambda
                    CSC_vort%b(Zones(k)%mesh%index(i,j))=r1+r2+r3
                 end if
              else
                 if(zones(k)%masks%chi3(i,j).eq.1) then
                    Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+max(min(Zones(k)%species(1)%var(1)%temperature(i,j+1)&
                         /Zones(k)%species(0)%var(1)%temperature(i,j+1),10.D0),0.1D0)))
                    if(mode_BC_para.eq.LAMBDA_TE) then
                       CSC_vort%b(Zones(k)%mesh%index(i,j))=Lambda*Zones(k)%species(0)%var(1)%temperature(i,j+1)
                    else
                       r1=-(zones(k)%metric_coefficients%G(i,j)+zones(k)%metric_coefficients%G(i,j+1))*0.5D0&
                            /((zones(k)%species(0)%var(1)%density(i,j+1)+zones(k)%species(0)%var(1)%density(i,j))*0.5d0)&
                            *(zones(k)%species(0)%var(1)%temperature(i,j+1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,j)**(1.5d0))*0.5D0&
                            *(zones(k)%species(0)%var(1)%density(i,j)*zones(k)%species(0)%var(1)%temperature(i,j)&
                            -zones(k)%species(0)%var(1)%density(i,j+1)*zones(k)%species(0)%var(1)%temperature(i,j+1))&
                            /(zones(k)%mesh%z(i,j)-zones(k)%mesh%z(i,j+1))
                       alpha=reference_parameters%fields%c0*reference_parameters%fields%W0*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                            *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
                       beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                            *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
                       if(drift_flags%jadvW) then
                          r2=alpha*(zones(k)%species(0)%var(1)%velocity(i,j+1)+&
                               R0/c0*(zones(k)%species(0)%drifts%uEt(i,j+1)+zones(k)%species(0)%drifts%uBt(i,j+1))/(zones(k)%metric_coefficients%G(i,j+1)))&
                               *zones(k)%electric_fields(1)%vorticity(i,j+1)
                       else
                          r2=0.d0
                       end if
                       
                       !                      r3=beta*zones(k)%species(0)%var(1)%Gamma(i,j+1)*(1.D0-exp(Lambda-zones(k)%electric_fields(1)%phi(i,j+1)/zones(k)%species(0)%var(1)%temperature(i,j+1))&
                       !                            *(1.D0+zones(k)%electric_fields(0)%phi(i,j)/zones(k)%species(0)%var(1)%temperature(i,j+1)))
                       Gamma=zones(k)%species(0)%var(1)%Gamma(i,j+1)+(zones(k)%species(0)%drifts%uEt(i,j+1)+zones(k)%species(0)%drifts%uBt(i,j+1))&
                            *(2.d0*pi*R0/rs0)*sqrt(zones(k)%metric_coefficients%ctt(i,j+1))/zones(k)%metric_coefficients%G(i,j+1)&
                            *zones(k)%species(0)%var(1)%density(i,j+1)    
                       r3=-beta*Gamma*Lambda
                       CSC_vort%b(Zones(k)%mesh%index(i,j))=r1+r2+r3
                    end if
                 else
                    if(zones(k)%masks%chi5(i,j).eq.1) then
                       CSC_vort%b(Zones(k)%mesh%index(i,j))=0.D0
                    else
                       if(zones(k)%masks%chi6(i,j).eq.1) then
                          CSC_vort%b(Zones(k)%mesh%index(i,j))=0.D0
                       else
                          CSC_vort%b(Zones(k)%mesh%index(i,j))=0.D0
                       end if
                    end if
                 end if
              end if
           end if
        end do
     end do

     East=Zones(k)%Neighbors(3)
     West=Zones(k)%Neighbors(4)

     do i=1,Nx
        if(West.gt.0) then
           CSC_vort%b(Zones(k)%mesh%index(i,0))=0.D0
        else
           Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+Zones(k)%species(1)%var(1)%temperature(i,0)&
                /Zones(k)%species(0)%var(1)%temperature(i,0)))
           if(mode_BC_para.eq.LAMBDA_TE) then
              CSC_vort%b(Zones(k)%mesh%index(i,0))=Lambda*Zones(k)%species(0)%var(1)%temperature(i,0)
           else
              r1=-(zones(k)%metric_coefficients%G(i,0)+zones(k)%metric_coefficients%G(i,1))*0.5D0&
                   /((zones(k)%species(0)%var(1)%density(i,1)+zones(k)%species(0)%var(1)%density(i,0))*0.5d0)&
                   *(zones(k)%species(0)%var(1)%temperature(i,1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,0)**(1.5d0))*0.5D0&
                   *(zones(k)%species(0)%var(1)%density(i,0)*zones(k)%species(0)%var(1)%temperature(i,0)&
                   -zones(k)%species(0)%var(1)%density(i,1)*zones(k)%species(0)%var(1)%temperature(i,1))&
                   /(zones(k)%mesh%z(i,0)-zones(k)%mesh%z(i,1))
              alpha=reference_parameters%fields%c0*reference_parameters%fields%W0*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                   *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
              beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                   *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
              if(drift_flags%jadvW) then
                 r2=alpha*(zones(k)%species(0)%var(1)%velocity(i,1)+&
                      R0/c0*(zones(k)%species(0)%drifts%uEt(i,1)+zones(k)%species(0)%drifts%uBt(i,1))/(zones(k)%metric_coefficients%G(i,1)))&
                      *zones(k)%electric_fields(1)%vorticity(i,1)
              else
                 r2=0.d0
              end if
              Gamma=zones(k)%species(0)%var(1)%Gamma(i,1)+(zones(k)%species(0)%drifts%uEt(i,1)+zones(k)%species(0)%drifts%uBt(i,1))&
                   *(2.d0*pi*R0/rs0)*sqrt(zones(k)%metric_coefficients%ctt(i,1))/zones(k)%metric_coefficients%G(i,1)&
                   *zones(k)%species(0)%var(1)%density(i,1)
              r3=-beta*Gamma*Lambda
              CSC_vort%b(Zones(k)%mesh%index(i,0))=r1+r2+r3
           end if
        end if
        if(East.gt.0) then
           CSC_vort%b(Zones(k)%mesh%index(i,Nz+1))=0.D0
        else
           Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+Zones(k)%species(1)%var(1)%temperature(i,Nz+1)&
                /Zones(k)%species(0)%var(1)%temperature(i,Nz+1)))
           if(mode_BC_para.eq.LAMBDA_TE) then
              CSC_vort%b(zones(k)%mesh%index(i,Nz+1))=Lambda*Zones(k)%species(0)%var(1)%temperature(i,Nz+1)
           else
              r1=-(zones(k)%metric_coefficients%G(i,Nz+1)+zones(k)%metric_coefficients%G(i,Nz))*0.5D0&
                   /((zones(k)%species(0)%var(1)%density(i,Nz)+zones(k)%species(0)%var(1)%density(i,Nz+1))*0.5d0)&
                   *(zones(k)%species(0)%var(1)%temperature(i,Nz)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,Nz+1)**(1.5d0))*0.5D0&
                   *(zones(k)%species(0)%var(1)%density(i,Nz+1)*zones(k)%species(0)%var(1)%temperature(i,Nz+1)&
                   -zones(k)%species(0)%var(1)%density(i,Nz)*zones(k)%species(0)%var(1)%temperature(i,Nz))&
                   /(zones(k)%mesh%z(i,Nz+1)-zones(k)%mesh%z(i,Nz))
              alpha=reference_parameters%fields%c0*reference_parameters%fields%W0*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                   *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
              beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                   *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
              if(drift_flags%jadvW) then
                 r2=alpha*(zones(k)%species(0)%var(1)%velocity(i,Nz)+&
                      R0/c0*(zones(k)%species(0)%drifts%uEt(i,Nz)+zones(k)%species(0)%drifts%uBt(i,Nz))/(zones(k)%metric_coefficients%G(i,Nz)))&
                      *zones(k)%electric_fields(1)%vorticity(i,Nz)
              else
                 r2=0.d0
              end if
              Gamma=zones(k)%species(0)%var(1)%Gamma(i,Nz)+(zones(k)%species(0)%drifts%uEt(i,Nz)+zones(k)%species(0)%drifts%uBt(i,Nz))&
                   *(2.d0*pi*R0/rs0)*sqrt(zones(k)%metric_coefficients%ctt(i,Nz))/zones(k)%metric_coefficients%G(i,Nz)&
                   *zones(k)%species(0)%var(1)%density(i,Nz)
              r3=-beta*Gamma*Lambda
              CSC_vort%b(Zones(k)%mesh%index(i,Nz+1))=r1+r2+r3
           end if
        end if
     end do

     !corners
     i=0
     CSC_vort%b(zones(k)%mesh%index(i,0))=0.D0
     CSC_vort%b(zones(k)%mesh%index(i,Nz+1))=0.D0
     i=Nx+1
     CSC_vort%b(zones(k)%mesh%index(i,0))=0.D0
     CSC_vort%b(zones(k)%mesh%index(i,Nz+1))=0.D0

     do j=1,Nz
        CSC_vort%b(zones(k)%mesh%index(0,j))=0.D0
        CSC_vort%b(zones(k)%mesh%index(Nx+1,j))=0.D0
     end do

     !North South BC
     South=zones(k)%Neighbors(2)
     North=zones(k)%Neighbors(1)
     if(mode_BC.eq.ZERO_FLUX) then
        if(North.eq.-1) then
           do j=1,Nz
!              CSC_vort%b(zones(k)%mesh%index(Nx+1,j))=0.d0
              CSC_vort%b(zones(k)%mesh%index(Nx+1,j))=-(zones(k)%species(1)%var(2)%density(Nx+1,j)*zones(k)%species(1)%var(2)%temperature(Nx+1,j)&
                   -zones(k)%species(1)%var(2)%density(Nx,j)*zones(k)%species(1)%var(2)%temperature(Nx,j))&
                   /(0.5d0*(zones(k)%species(1)%var(2)%density(Nx,j)+zones(k)%species(1)%var(2)%density(Nx+1,j)))
           end do
        end if
        if(South.eq.-1) then
           do j=1,Nz
!              CSC_vort%b(zones(k)%mesh%index(0,j))=0.d0
              CSC_vort%b(zones(k)%mesh%index(0,j))=-(zones(k)%species(1)%var(2)%density(0,j)*zones(k)%species(1)%var(2)%temperature(0,j)&
                   -zones(k)%species(1)%var(2)%density(1,j)*zones(k)%species(1)%var(2)%temperature(1,j))&
                   /(0.5d0*(zones(k)%species(1)%var(2)%density(0,j)+zones(k)%species(1)%var(2)%density(1,j)))
!              CSC_vort%b(zones(k)%mesh%index(0,j))=0.D0!Source(1,j)
           end do
        end if
     else
        if(North.eq.-1) then
           do j=1,Nz
              CSC_vort%b(zones(k)%mesh%index(Nx+1,j))=test_vort(k)%phi(Nx+1,j)/reference_parameters%fields%phi0
           end do
        end if
        if(South.eq.-1) then
           do j=1,Nz
              CSC_vort%b(zones(k)%mesh%index(0,j))=test_vort(k)%phi(0,j)/reference_parameters%fields%phi0
           end do
        end if
     end if

  end do
end subroutine compute_vorticity_RHS4
