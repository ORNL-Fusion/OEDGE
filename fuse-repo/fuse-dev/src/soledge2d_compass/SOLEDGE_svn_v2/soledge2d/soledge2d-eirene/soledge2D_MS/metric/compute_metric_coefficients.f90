subroutine compute_metric_coefficients()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  real*8 :: DVOL
  real*8 :: R0,rs0
  integer*4 :: Nx,Nz
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     ! computing Jacobian J = det (Jac)
     do i=0,Nx+1
        do j=0,Nz+1
           if((zones(k)%jacobian%dPdR(i,j)*zones(k)%jacobian%dTdz(i,j)&
                -zones(k)%jacobian%dPdz(i,j)*zones(k)%jacobian%dTdR(i,j)).eq.0.D0) then
              zones(k)%metric_coefficients%Jacobian(i,j)=0.D0
           else
              zones(k)%metric_coefficients%Jacobian(i,j)=abs(zones(k)%mesh%Rgeom(i,j)&
                   /(zones(k)%jacobian%dPdR(i,j)*zones(k)%jacobian%dTdz(i,j)&
                   -zones(k)%jacobian%dPdz(i,j)*zones(k)%jacobian%dTdR(i,j)))
           end if
        end do
     end do
     do i=0,Nx+1
        do j=0,Nz+1
           ! ctt = grad(theta*) x grad(theta*) / (2*pi/rs0)^2b   [x denotes scalar product]
           zones(k)%metric_coefficients%ctt(i,j)=(zones(k)%jacobian%dTdz(i,j)**(2.d0)&
                +zones(k)%jacobian%dTdR(i,j)**(2.d0))&
                *rs0**(2.d0)/(2.D0*pi)**2.D0
           ! cpp = grad(psi*) x grad(psi*) / (Dpsi0/rs0)^2
           zones(k)%metric_coefficients%cpp(i,j)=(zones(k)%jacobian%dPdz(i,j)**(2.d0)&
                +zones(k)%jacobian%dPdR(i,j)**(2.d0))&
                *rs0**(2.d0)
           ! cpt = grad(psi*) x grad(theta*) / (2*pi*Dpsi0/rs0^2)
           zones(k)%metric_coefficients%cpt(i,j)=(zones(k)%jacobian%dPdz(i,j)*zones(k)%jacobian%dTdz(i,j)&
                +zones(k)%jacobian%dPdR(i,j)*zones(k)%jacobian%dTdR(i,j))&
                *rs0**(2.d0)/(2.D0*pi)
           if((zones(k)%metric_coefficients%ctt(i,j)*zones(k)%metric_coefficients%cpp(i,j)&
                -zones(k)%metric_coefficients%cpt(i,j)**(2.d0)).ne.0.D0) then
              ! c_tt = dR/dtheta* x dR/dtheta* * (2*pi/rs0)^2        [R=(R,Z) - vector position]
              zones(k)%metric_coefficients%c_tt(i,j)=1.d0/(zones(k)%metric_coefficients%ctt(i,j)*zones(k)%metric_coefficients%cpp(i,j)&
                   -zones(k)%metric_coefficients%cpt(i,j)**(2.d0))*zones(k)%metric_coefficients%cpp(i,j)
              ! c_pp = dR/dpsi* x dR/dpsi* * (Dpsi0/rs0)^2
              zones(k)%metric_coefficients%c_pt(i,j)=-1.d0/(zones(k)%metric_coefficients%ctt(i,j)*zones(k)%metric_coefficients%cpp(i,j)&
                   -zones(k)%metric_coefficients%cpt(i,j)**(2.d0))*zones(k)%metric_coefficients%cpt(i,j)
              ! c_pt = dR/dpsi* x dR/dtheta* * (2*pi*Dpsi0/rs0^2)
              zones(k)%metric_coefficients%c_pp(i,j)=1.d0/(zones(k)%metric_coefficients%ctt(i,j)*zones(k)%metric_coefficients%cpp(i,j)&
                   -zones(k)%metric_coefficients%cpt(i,j)**(2.d0))*zones(k)%metric_coefficients%ctt(i,j)
           end if
           ! G = B x grad(theta*) / |B|  * R0
           if(zones(k)%mesh%B(i,j).ne.0.d0) then
              zones(k)%metric_coefficients%G(i,j)=(zones(k)%mesh%Br(i,j)*zones(k)%jacobian%dtdR(i,j)&
                   +zones(k)%mesh%Bz(i,j)*zones(k)%jacobian%dtdz(i,j))&
                   /(zones(k)%mesh%B(i,j))*R0
           else
              zones(k)%metric_coefficients%G(i,j)=1.d-15
           end if
        end do
     end do
  end do
  call compute_diffareas()
  call compute_dvol()
  call compute_sinepitch()
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx 
        do j=1,Nz
           DVOL=zones(k)%metric_coefficients%Dvol_dd(i,j)
           zones(k)%metric_coefficients%divb(i,j)=(zones(k)%metric_coefficients%ds_east_DD(i,j)*zones(k)%metric_coefficients%sinepitch_east(i,j)&
                -zones(k)%metric_coefficients%ds_west_DD(i,j)*zones(k)%metric_coefficients%sinepitch_west(i,j)&
                )/Dvol*(2.D0*pi*R0/rs0) 
        end do
     end do
  end do
  call MD_broadcast_metric()
  
  call compute_Score()
  call save_metric()
end subroutine compute_metric_coefficients
