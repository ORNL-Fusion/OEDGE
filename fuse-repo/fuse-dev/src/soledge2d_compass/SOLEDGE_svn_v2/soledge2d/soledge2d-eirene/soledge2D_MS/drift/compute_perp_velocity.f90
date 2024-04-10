subroutine compute_perp_velocity(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables, drift_flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: n,Nx,Nz,i,j
  real*8, allocatable :: Pr(:,:)
  real*8, allocatable :: dPhidT(:,:), dPhidP(:,:)
  real*8, allocatable :: dPrdT(:,:), dPrdP(:,:)
  real*8, allocatable :: dBdT(:,:), dBdP(:,:)
  real*8, allocatable :: Jacobian(:,:)
  real*8,parameter :: eps=1.d-6
  real*8 :: c0,R0,rs0
  c0=reference_parameters%fields%c0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0

  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     allocate(Pr(0:Nx+1,0:Nz+1))
     allocate(dPhidP(0:Nx+1,0:Nz+1),dPhidT(0:Nx+1,0:Nz+1))
     allocate(dPrdP(0:Nx+1,0:Nz+1),dPrdT(0:Nx+1,0:Nz+1))
     allocate(dBdP(0:Nx+1,0:Nz+1),dBdT(0:Nx+1,0:Nz+1))
     allocate(Jacobian(0:Nx+1,0:Nz+1))
     Pr(0:Nx+1,0:Nz+1) = zone%species(n)%var(1)%density(0:Nx+1,0:Nz+1)*zone%species(n)%var(1)%temperature(0:Nx+1,0:Nz+1)

!!$     do i=0,Nx+1
!!$        write(5000+zone%number,100) (zone%electric_fields(1)%phi_smooth(i,j),j=0,Nz+1)
!!$100     format(512es15.7)
!!$     end do
     Jacobian(0:Nx+1,0:Nz+1)=(zone%mesh%Rgeom(0:Nx+1,0:Nz+1)&
          /(zone%jacobian%dPdR(0:Nx+1,0:Nz+1)*zone%jacobian%dTdz(0:Nx+1,0:Nz+1)&
          -zone%jacobian%dPdz(0:Nx+1,0:Nz+1)*zone%jacobian%dTdR(0:Nx+1,0:Nz+1)))


     call gradient_centre(Nx,Nz,zone%electric_fields(1)%phi_smooth,dPhidP,dPhidT,zone%mesh%x,2.d0*pi*zone%mesh%z)
     call gradient_centre(Nx,Nz,Pr,dPrdP,dPrdT,zone%mesh%x,2.d0*pi*zone%mesh%z)
     call gradient_centre(Nx,Nz,zone%mesh%B,dBdP,dBdT,zone%mesh%x,2.d0*pi*zone%mesh%z)
     ! Electric drift
     zone%species(n)%drifts%uEp(0:Nx+1,0:Nz+1)=-reference_parameters%fields%phi0*dPhidT(0:Nx+1,0:Nz+1)*zone%mesh%Rgeom(0:Nx+1,0:Nz+1)*&
          zone%mesh%Bphi(0:Nx+1,0:Nz+1)/((Jacobian(0:Nx+1,0:Nz+1)+eps)*(zone%mesh%B(0:Nx+1,0:Nz+1) **(2.d0)+ eps))&
          /(sqrt(zone%metric_coefficients%cpp(0:Nx+1,0:Nz+1))/rs0+eps)
     zone%species(n)%drifts%uEt(0:Nx+1,0:Nz+1)=reference_parameters%fields%phi0*dPhidP(0:Nx+1,0:Nz+1)*zone%mesh%Rgeom(0:Nx+1,0:Nz+1)*&
          zone%mesh%Bphi(0:Nx+1,0:Nz+1)/((Jacobian(0:Nx+1,0:Nz+1)+eps)*(zone%mesh%B(0:Nx+1,0:Nz+1) **(2.d0)+ eps))&
          /(sqrt(zone%metric_coefficients%ctt(0:Nx+1,0:Nz+1))/(rs0/(2.D0*pi))+eps)
     !adim
     zone%species(n)%drifts%uEp(0:Nx+1,0:Nz+1)=zone%species(n)%drifts%uEp(0:Nx+1,0:Nz+1)/c0
     zone%species(n)%drifts%uEt(0:Nx+1,0:Nz+1)=zone%species(n)%drifts%uEt(0:Nx+1,0:Nz+1)/c0

     ! Diamagnetic drift
     zone%species(n)%drifts%udp(0:Nx+1,0:Nz+1)=-reference_parameters%fields%T0eV/(global_parameters%element_list(n)%Z)*dPrdT(0:Nx+1,0:Nz+1)*&
          zone%mesh%Rgeom(0:Nx+1,0:Nz+1)*zone%mesh%Bphi(0:Nx+1,0:Nz+1)/&
          ((zone%species(n)%var(1)%density(0:Nx+1,0:Nz+1)+eps)*(Jacobian(0:Nx+1,0:Nz+1)+eps)*(zone%mesh%B(0:Nx+1,0:Nz+1) **(2.d0)+ eps))&
          /(sqrt(zone%metric_coefficients%cpp(0:Nx+1,0:Nz+1))/rs0+eps)
     zone%species(n)%drifts%udt(0:Nx+1,0:Nz+1)=reference_parameters%fields%T0eV/global_parameters%element_list(n)%Z*dPrdP(0:Nx+1,0:Nz+1)*&
          zone%mesh%Rgeom(0:Nx+1,0:Nz+1)*zone%mesh%Bphi(0:Nx+1,0:Nz+1)/&
          ((zone%species(n)%var(1)%density(0:Nx+1,0:Nz+1)+eps)*(Jacobian(0:Nx+1,0:Nz+1)+eps)*(zone%mesh%B(0:Nx+1,0:Nz+1) **(2.d0)+ eps))&
          /(sqrt(zone%metric_coefficients%ctt(0:Nx+1,0:Nz+1))/(rs0/(2.D0*pi))+eps)
     zone%species(n)%drifts%udp(0:Nx+1,0:Nz+1)=zone%species(n)%drifts%udp(0:Nx+1,0:Nz+1)/c0
     zone%species(n)%drifts%udt(0:Nx+1,0:Nz+1)=zone%species(n)%drifts%udt(0:Nx+1,0:Nz+1)/c0

     ! Gradient/curvature drift
     zone%species(n)%drifts%uBp(0:Nx+1,0:Nz+1)=-2.d0*reference_parameters%fields%T0eV/global_parameters%element_list(n)%Z*&
          zone%species(n)%var(1)%temperature(0:Nx+1,0:Nz+1)*dBdT(0:Nx+1,0:Nz+1)*zone%mesh%Rgeom(0:Nx+1,0:Nz+1)*zone%mesh%Bphi(0:Nx+1,0:Nz+1)/&
          ((Jacobian(0:Nx+1,0:Nz+1)+eps)*(zone%mesh%B(0:Nx+1,0:Nz+1) **(3.d0)+ eps))&
          /(sqrt(zone%metric_coefficients%cpp(0:Nx+1,0:Nz+1))/rs0+eps)
     zone%species(n)%drifts%uBt(0:Nx+1,0:Nz+1)=2.d0*reference_parameters%fields%T0eV/global_parameters%element_list(n)%Z*&
          zone%species(n)%var(1)%temperature(0:Nx+1,0:Nz+1)*dBdP(0:Nx+1,0:Nz+1)*zone%mesh%Rgeom(0:Nx+1,0:Nz+1)*zone%mesh%Bphi(0:Nx+1,0:Nz+1)/&
          ((Jacobian(0:Nx+1,0:Nz+1)+eps)*(zone%mesh%B(0:Nx+1,0:Nz+1) **(3.d0)+ eps))&
          /(sqrt(zone%metric_coefficients%ctt(0:Nx+1,0:Nz+1))/(rs0/(2.D0*pi))+eps)
     zone%species(n)%drifts%uBp(0:Nx+1,0:Nz+1)=zone%species(n)%drifts%uBp(0:Nx+1,0:Nz+1)/c0
     zone%species(n)%drifts%uBt(0:Nx+1,0:Nz+1)=zone%species(n)%drifts%uBt(0:Nx+1,0:Nz+1)/c0

     if(drift_flags%reverse_B) then
        zone%species(n)%drifts%uBt=-zone%species(n)%drifts%uBt
        zone%species(n)%drifts%uEt=-zone%species(n)%drifts%uEt
        zone%species(n)%drifts%uBp=-zone%species(n)%drifts%uBp
        zone%species(n)%drifts%uEp=-zone%species(n)%drifts%uEp
        zone%species(n)%drifts%udp=-zone%species(n)%drifts%udp
        zone%species(n)%drifts%udp=-zone%species(n)%drifts%udp
     end if

     zone%species(n)%drifts%uBts=zone%species(n)%drifts%uBt
     zone%species(n)%drifts%uBps=zone%species(n)%drifts%uBp
     zone%species(n)%drifts%uEts=zone%species(n)%drifts%uEt
     zone%species(n)%drifts%uEps=zone%species(n)%drifts%uEp


     deallocate(Pr,dPhidP,dPhidT,dPrdP,dPrdT,dBdP,dBdT,Jacobian)
     ! broadcast and computation of derivatives on the edges
     !      call gradient_BC_dpsi(n)
     !      call gradient_BC_dtheta(n)    
     if(.not.drift_flags%use_gradB_radial) then
        zone%species(n)%drifts%uBp=0.d0
     end if
     if(.not.drift_flags%use_gradB_poloidal) then
        zone%species(n)%drifts%uBt=0.d0
     end if
     if(.not.drift_flags%use_ExB_radial) then
        zone%species(n)%drifts%uEp=0.d0
     end if
     if(.not.drift_flags%use_ExB_poloidal) then
        zone%species(n)%drifts%uEt=0.d0
     end if


  end do
end subroutine compute_perp_velocity

!!###########################################################################################################################################!!

subroutine gradient_centre(Nx,Nz,X,dXdP,dXdT,P,T)
  use Mphysics
  implicit none
  integer*4 Nx
  integer*4 Nz
  !   real*8,dimension(0:Nx+1,0:Nz+1),intent(in) :: X
  real*8,dimension(0:Nx+1,0:Nz+1) :: X
  real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdP
  real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdT
  real*8,dimension(0:Nx+1,0:Nz+1) :: P
  real*8,dimension(0:Nx+1,0:Nz+1) :: T
  integer*4 i,j
  !second order derivative in the middle
  do i=1,Nx
     do j=1,Nz
        dXdP(i,j)=(X(i+1,j)*(P(i,j)-P(i-1,j))**2.d0&
             -X(i-1,j)*(P(i+1,j)-P(i,j))**2.d0&
             +X(i,j)*((P(i+1,j)-P(i,j))**2.d0-(P(i,j)-P(i-1,j))**2.d0))&
             /((P(i+1,j)-P(i,j))*(P(i,j)-P(i-1,j))*(P(i+1,j)-P(i-1,j)))
        dXdT(i,j)=(X(i,j+1)*(T(i,j)-T(i,j-1))**2.d0&
             -X(i,j-1)*(T(i,j+1)-T(i,j))**2.d0&
             +X(i,j)*((T(i,j+1)-T(i,j))**2.d0&
             -(T(i,j)-T(i,j-1))**2.d0))&
             /((T(i,j+1)-T(i,j))*(T(i,j)-T(i,j-1))&
             *(T(i,j+1)-T(i,j-1)))     
!!$   dXdT(i,j)=(X(i,j+1)*modulo((T(i,j)-T(i,j-1)),2.d0*pi)**2.d0&
!!$             -X(i,j-1)*modulo((T(i,j+1)-T(i,j)),2.d0*pi)**2.d0&
!!$             +X(i,j)*(modulo((T(i,j+1)-T(i,j)),2.d0*pi)**2.d0&
!!$             -modulo((T(i,j)-T(i,j-1)),2.d0*pi)**2.d0))&
!!$             /(modulo((T(i,j+1)-T(i,j)),2.d0*pi)*modulo((T(i,j)-T(i,j-1)),2.d0*pi)&
!!$             *modulo((T(i,j+1)-T(i,j-1)),2.d0*pi))
     end do
  end do
end subroutine gradient_centre

!!###########################################################################################################################################!!

subroutine gradient_BC_dpsi(Nspec)
  use all_variables, only : zones, reference_parameters, global_parameters
  implicit none
  integer*4,intent(inout) :: Nspec
  integer*4 North,South,East,West
  integer*4 mNorth,mSouth,mEast,mWest
  real*8, allocatable :: Pr(:,:)
  real*8, allocatable :: dPhidP(:,:)
  real*8, allocatable :: dPrdP(:,:)
  real*8, allocatable :: dBdP(:,:)
  integer*4 k,Nx,Nz,i
  integer*4 Vois_info(4)
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(Pr(0:Nx+1,0:Nz+1))
     allocate(dPhidP(0:Nx+1,0:Nz+1))!,dPhidT(0:Nx+1,0:Nz+1))
     allocate(dPrdP(0:Nx+1,0:Nz+1))!,dPrdT(0:Nx+1,0:Nz+1))
     allocate(dBdP(0:Nx+1,0:Nz+1))!,dBdT(0:Nx+1,0:Nz+1))
     Vois_info(1)=zones(k)%Neighbors(1)
     Vois_info(2)=zones(k)%Neighbors(2)
     Vois_info(3)=zones(k)%Neighbors(3)
     Vois_info(4)=zones(k)%Neighbors(4)
     Pr(0:Nx+1,0:Nz+1) = zones(k)%species(Nspec)%var(1)%density(0:Nx+1,0:Nz+1)*zones(k)%species(Nspec)%var(1)%temperature(0:Nx+1,0:Nz+1)
     ! on top
     North=zones(k)%Neighbors(1)
     mNorth=zones(k)%MagNeighbors(1)
     if(mNorth.eq.1) then
        call compute_Neighdrift_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,1)
        !!call compute_Neigh_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,1,Vois_info,zones(k)%mesh%cornerX,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dpsi(Nx,Nz,Pr,zones(k)%mesh%x,dPrdP,1)
	call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%B,zones(k)%mesh%x,dBdP,1,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEt(Nx,1:Nz)=reference_parameters%fields%phi0*dPhidP(Nx,1:Nz)*zones(k)%mesh%Rgeom(Nx,1:Nz)*&
		zones(k)%mesh%Bphi(Nx,1:Nz)/(zones(k)%metric_coefficients%Jacobian(Nx,1:Nz)*zones(k)%mesh%B(Nx,1:Nz)**(2.d0))
        zones(k)%species(Nspec)%drifts%uEt(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uEt(Nx,1:Nz)	
        !!zones(k)%species(Nspec)%drifts%uEt(Nx+1,1:Nz)=reference_parameters%fields%phi0*dPhidP(Nx+1,1:Nz)*zones(k)%mesh%Rgeom(Nx+1,1:Nz)*&
	!!	zones(k)%mesh%Bphi(Nx+1,1:Nz)/(zones(k)%metric_coefficients%Jacobian(Nx+1,1:Nz)*zones(k)%mesh%B(Nx+1,1:Nz)**(2.d0))	
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udt(Nx,1:Nz)=reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdP(Nx,1:Nz)*&
		zones(k)%mesh%Rgeom(Nx,1:Nz)*zones(k)%mesh%Bphi(Nx,1:Nz)/&
		(zones(k)%species(Nspec)%var(1)%density(Nx,1:Nz)*zones(k)%metric_coefficients%Jacobian(Nx,1:Nz)*zones(k)%mesh%B(Nx,1:Nz)**(2.d0))
        zones(k)%species(Nspec)%drifts%udt(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%udt(Nx,1:Nz)
	! Curvature drift
        zones(k)%species(Nspec)%drifts%uBt(Nx+1,1:Nz)=2.d0*reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*&
		zones(k)%species(Nspec)%var(1)%temperature(Nx,1:Nz)*dBdP(Nx+1,1:Nz)*zones(k)%mesh%Rgeom(Nx+1,1:Nz)*zones(k)%mesh%Bphi(Nx+1,1:Nz)/&
		(zones(k)%metric_coefficients%Jacobian(Nx+1,1:Nz)*zones(k)%mesh%B(Nx+1,1:Nz)**(3.d0))
     end if
     ! at the bottom
     South=Zones(k)%Neighbors(2)
     mSouth=Zones(k)%MagNeighbors(2)
     if(mSouth.eq.1) then
        call compute_Neighdrift_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,2)
        !!call compute_Neigh_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,2,Vois_info,zones(k)%mesh%cornerX,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dpsi(Nx,Nz,Pr,zones(k)%mesh%x,dPrdP,2)
	call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%B,zones(k)%mesh%x,dBdP,2,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEt(1,1:Nz)=reference_parameters%fields%phi0*dPhidP(1,1:Nz)*zones(k)%mesh%Rgeom(1,1:Nz)*&
		zones(k)%mesh%Bphi(1,1:Nz)/(zones(k)%metric_coefficients%Jacobian(1,1:Nz)*zones(k)%mesh%B(1,1:Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%uEt(0,1:Nz)=zones(k)%species(Nspec)%drifts%uEt(1,1:Nz)		
        !!zones(k)%species(Nspec)%drifts%uEt(0,1:Nz)=reference_parameters%fields%phi0*dPhidP(0,1:Nz)*zones(k)%mesh%Rgeom(0,1:Nz)*&
	!!	zones(k)%mesh%Bphi(0,1:Nz)/(zones(k)%metric_coefficients%Jacobian(0,1:Nz)*zones(k)%mesh%B(0,1:Nz)**(2.d0))
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udt(1,1:Nz)=reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdP(1,1:Nz)*&
		zones(k)%mesh%Rgeom(1,1:Nz)*zones(k)%mesh%Bphi(1,1:Nz)/&
		(zones(k)%species(Nspec)%var(1)%density(1,1:Nz)*zones(k)%metric_coefficients%Jacobian(1,1:Nz)*zones(k)%mesh%B(1,1:Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%udt(0,1:Nz)=zones(k)%species(Nspec)%drifts%udt(1,1:Nz)
	! Curvature drift
        zones(k)%species(Nspec)%drifts%uBt(0,1:Nz)=2.d0*reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*&
		zones(k)%species(Nspec)%var(1)%temperature(0,1:Nz)*dBdP(0,1:Nz)*zones(k)%mesh%Rgeom(0,1:Nz)*zones(k)%mesh%Bphi(0,1:Nz)/&
		(zones(k)%metric_coefficients%Jacobian(0,1:Nz)*zones(k)%mesh%B(0,1:Nz)**(3.d0))
     end if
     ! on the right
     East=Zones(k)%Neighbors(3)
     mEast=Zones(k)%MagNeighbors(3)
     if(mEast.eq.1) then
        call compute_Neighdrift_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,3)
        !!call compute_Neigh_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,3,Vois_info,zones(k)%mesh%cornerX,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dpsi(Nx,Nz,Pr,zones(k)%mesh%x,dPrdP,3)
	call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%B,zones(k)%mesh%x,dBdP,3,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz)=reference_parameters%fields%phi0*dPhidP(1:Nx,Nz)*zones(k)%mesh%Rgeom(1:Nx,Nz)*&
		zones(k)%mesh%Bphi(1:Nx,Nz)/(zones(k)%metric_coefficients%Jacobian(1:Nx,Nz)*zones(k)%mesh%B(1:Nx,Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz)	
        !!zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz+1)=reference_parameters%fields%phi0*dPhidP(1:Nx,Nz+1)*zones(k)%mesh%Rgeom(1:Nx,Nz+1)*&
	!!	zones(k)%mesh%Bphi(1:Nx,Nz+1)/(zones(k)%metric_coefficients%Jacobian(1:Nx,Nz+1)*zones(k)%mesh%B(1:Nx,Nz+1)**(2.d0))	
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udt(1:Nx,Nz)=reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdP(1:Nx,Nz)*&
		zones(k)%mesh%Rgeom(1:Nx,Nz)*zones(k)%mesh%Bphi(1:Nx,Nz)/&
		(zones(k)%species(Nspec)%var(1)%density(1:Nx,Nz)*zones(k)%metric_coefficients%Jacobian(1:Nx,Nz)*zones(k)%mesh%B(1:Nx,Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%udt(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%udt(1:Nx,Nz)
	! Curvature drift
        zones(k)%species(Nspec)%drifts%uBt(1:Nx,Nz+1)=2.d0*reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*&
		zones(k)%species(Nspec)%var(1)%temperature(1:Nx,Nz+1)*dBdP(1:Nx,Nz+1)*zones(k)%mesh%Rgeom(1:Nx,Nz+1)*zones(k)%mesh%Bphi(1:Nx,Nz+1)/&
		(zones(k)%metric_coefficients%Jacobian(1:Nx,Nz+1)*zones(k)%mesh%B(1:Nx,Nz+1)**(3.d0))

     end if
     ! on the left
     West=Zones(k)%Neighbors(4)
     mWest=Zones(k)%MagNeighbors(4)
     if(mWest.eq.1) then
        call compute_Neighdrift_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,4)
        !!call compute_Neigh_dpsi(Nx,Nz,zones(k)%electric_fields(1)%phi,zones(k)%mesh%x,dPhidP,4,Vois_info,zones(k)%mesh%cornerX,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dpsi(Nx,Nz,Pr,zones(k)%mesh%x,dPrdP,4)
	call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%B,zones(k)%mesh%x,dBdP,4,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEt(1:Nx,1)=reference_parameters%fields%phi0*dPhidP(1:Nx,1)*zones(k)%mesh%Rgeom(1:Nx,1)*&
		zones(k)%mesh%Bphi(1:Nx,1)/(zones(k)%metric_coefficients%Jacobian(1:Nx,1)*zones(k)%mesh%B(1:Nx,1)**(2.d0))	
	zones(k)%species(Nspec)%drifts%uEt(1:Nx,0)=zones(k)%species(Nspec)%drifts%uEt(1:Nx,1)	
        !!zones(k)%species(Nspec)%drifts%uEt(1:Nx,0)=reference_parameters%fields%phi0*dPhidP(1:Nx,0)*zones(k)%mesh%Rgeom(1:Nx,0)*&
	!!	zones(k)%mesh%Bphi(1:Nx,0)/(zones(k)%metric_coefficients%Jacobian(1:Nx,0)*zones(k)%mesh%B(1:Nx,0)**(2.d0))
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udt(1:Nx,1)=reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdP(1:Nx,1)*&
		zones(k)%mesh%Rgeom(1:Nx,1)*zones(k)%mesh%Bphi(1:Nx,1)/&
		(zones(k)%species(Nspec)%var(1)%density(1:Nx,1)*zones(k)%metric_coefficients%Jacobian(1:Nx,1)*zones(k)%mesh%B(1:Nx,1)**(2.d0))
	zones(k)%species(Nspec)%drifts%udt(1:Nx,0)=zones(k)%species(Nspec)%drifts%udt(1:Nx,1)
	! Curvature drift
        zones(k)%species(Nspec)%drifts%uBt(1:Nx,0)=2.d0*reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*&
		zones(k)%species(Nspec)%var(1)%temperature(1:Nx,0)*dBdP(1:Nx,0)*zones(k)%mesh%Rgeom(1:Nx,0)*zones(k)%mesh%Bphi(1:Nx,0)/&
		(zones(k)%metric_coefficients%Jacobian(1:Nx,0)*zones(k)%mesh%B(1:Nx,0)**(3.d0))
     end if
  deallocate(Pr,dPhidP,dPrdP,dBdP)
  end do
end subroutine gradient_BC_dpsi

subroutine compute_Neighdrift_dpsi(Nx,Nz,X,P,dXdP,Neigh)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: X
  real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: P
  real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdP
  integer*4 Neigh !1N, 2S, 3E, 4W
  select case(Neigh)
  case(1)        ! on top
     dXdP(Nx,1:Nz)=(X(Nx,1:Nz)-X(Nx-1,1:Nz))/(P(Nx,1:Nz)-P(Nx-1,1:Nz))
  case(2)       ! at the bottom
     dXdP(1,1:Nz)=(X(2,1:Nz)-X(1,1:Nz))/(P(2,1:Nz)-P(1,1:Nz))
  case(3)       ! on the right
     dXdP(2:Nx-1,Nz)=(X(3:Nx,Nz)*(P(2:Nx-1,Nz)-P(1:Nx-2,Nz))**2.d0&
          -X(1:Nx-2,Nz)*(P(3:Nx,Nz)-P(2:Nx-1,Nz))**2.d0&
          +X(2:Nx-1,Nz)*((P(3:Nx,Nz)-P(2:Nx-1,Nz))**2.d0-(P(2:Nx-1,Nz)-P(1:Nx-2,Nz))**2.d0))&
          /((P(3:Nx,Nz)-P(2:Nx-1,Nz))*(P(2:Nx-1,Nz)-P(1:Nx-2,Nz))*(P(3:Nx,Nz)-P(1:Nx-2,Nz)))
     !dXdP(2:Nx-1,Nz) = (X(3:Nx,Nz)-X(1:Nx-2,Nz))/(P(3:Nx,Nz)-P(1:Nx-2,Nz))
     dXdP(1,Nz)=(X(2,Nz)-X(1,Nz))/(P(2,Nz)-P(1,Nz))
     dXdP(Nx,Nz)=(X(Nx,Nz)-X(Nx-1,Nz))/(P(Nx,Nz)-P(Nx-1,Nz))
  case(4)      ! on the left
     dXdP(2:Nx-1,1)=(X(3:Nx,1)*(P(2:Nx-1,1)-P(1:Nx-2,1))**2.d0&
          -X(1:Nx-2,1)*(P(3:Nx,1)-P(2:Nx-1,1))**2.d0&
          +X(2:Nx-1,1)*((P(3:Nx,1)-P(2:Nx-1,1))**2.d0-(P(2:Nx-1,1)-P(1:Nx-2,1))**2.d0))&
          /((P(3:Nx,1)-P(2:Nx-1,1))*(P(2:Nx-1,1)-P(1:Nx-2,1))*(P(3:Nx,1)-P(1:Nx-2,1)))
     !dXdP(2:Nx-1,1) = (X(3:Nx,1)-X(1:Nx-2,1))/(P(3:Nx,1)-P(1:Nx-2,1))
     dXdP(1,1)=(X(2,1)-X(1,1))/(P(2,1)-P(1,1))
     dXdP(Nx,1)=(X(Nx,1)-X(Nx-1,1))/(P(Nx,1)-P(Nx-1,1))
  case default
  end select
end subroutine compute_Neighdrift_dpsi

!!###########################################################################################################################################!!

subroutine gradient_BC_dtheta(Nspec)
  use all_variables, only : zones, reference_parameters, global_parameters
  use Mphysics
  implicit none
  integer*4,intent(inout) :: Nspec
  integer*4 North,South,East,West
  integer*4 mNorth,mSouth,mEast,mWest
  real*8, allocatable :: Pr(:,:)
  real*8, allocatable :: dPhidT(:,:)
  real*8, allocatable :: dPrdT(:,:)
  real*8, allocatable :: dBdT(:,:)
  integer*4 k,Nx,Nz,j
  integer*4 Vois_info(4)
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(Pr(0:Nx+1,0:Nz+1))
     allocate(dPhidT(0:Nx+1,0:Nz+1))
     allocate(dPrdT(0:Nx+1,0:Nz+1))
     allocate(dBdT(0:Nx+1,0:Nz+1))
     Vois_info(1)=zones(k)%Neighbors(1)
     Vois_info(2)=zones(k)%Neighbors(2)
     Vois_info(3)=zones(k)%Neighbors(3)
     Vois_info(4)=zones(k)%Neighbors(4)
     Pr(0:Nx+1,0:Nz+1) = zones(k)%species(Nspec)%var(1)%density(0:Nx+1,0:Nz+1)*zones(k)%species(Nspec)%var(1)%temperature(0:Nx+1,0:Nz+1)
     ! on top
     North=zones(k)%Neighbors(1)
     mNorth=zones(k)%MagNeighbors(1)
     if(mNorth.eq.1) then
        call compute_Neighdrift_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,1)
	!!call compute_Neigh_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,1,Vois_info,&
	!!	zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dtheta(Nx,Nz,Pr,2.d0*pi*zones(k)%mesh%z,dPrdT,1)
	call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%B,2.d0*pi*zones(k)%mesh%z,dBdT,1,Vois_info,zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEp(Nx,1:Nz)=-reference_parameters%fields%phi0*dPhidT(Nx,1:Nz)*zones(k)%mesh%Rgeom(Nx,1:Nz)*&
		zones(k)%mesh%Bphi(Nx,1:Nz)/(zones(k)%metric_coefficients%Jacobian(Nx,1:Nz)*zones(k)%mesh%B(Nx,1:Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%uEp(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uEp(Nx,1:Nz)		
        !!zones(k)%species(Nspec)%drifts%uEp(Nx+1,1:Nz)=-reference_parameters%fields%phi0*dPhidT(Nx+1,1:Nz)*zones(k)%mesh%Rgeom(Nx+1,1:Nz)*&
	!!	zones(k)%mesh%Bphi(Nx+1,1:Nz)/(zones(k)%metric_coefficients%Jacobian(Nx+1,1:Nz)*zones(k)%mesh%B(Nx+1,1:Nz)**(2.d0))
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udp(Nx,1:Nz)=-reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdT(Nx,1:Nz)*&
		zones(k)%mesh%Rgeom(Nx,1:Nz)*zones(k)%mesh%Bphi(Nx,1:Nz)/&
		(zones(k)%species(Nspec)%var(1)%density(Nx,1:Nz)*zones(k)%metric_coefficients%Jacobian(Nx,1:Nz)*zones(k)%mesh%B(Nx,1:Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%udp(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%udp(Nx,1:Nz)
        ! Gradient/curvature drift
        zones(k)%species(Nspec)%drifts%uBp(Nx+1,1:Nz)=-2.d0*reference_parameters%fields%T0eV/(global_parameters%element_list(Nspec)%Z)*&
		zones(k)%species(Nspec)%var(1)%temperature(Nx+1,1:Nz)*dBdT(Nx+1,1:Nz)*zones(k)%mesh%Rgeom(Nx+1,1:Nz)*zones(k)%mesh%Bphi(Nx+1,1:Nz)/&
		(zones(k)%metric_coefficients%Jacobian(Nx+1,1:Nz)*zones(k)%mesh%B(Nx+1,1:Nz)**(3.d0))
     end if
     ! at the bottom
     South=Zones(k)%Neighbors(2)
     mSouth=Zones(k)%MagNeighbors(2)
     if(mSouth.eq.1) then
        call compute_Neighdrift_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,2)
	!!call compute_Neigh_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,2,Vois_info,&
	!!	zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dtheta(Nx,Nz,Pr,2.d0*pi*zones(k)%mesh%z,dPrdT,2)
	call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%B,2.d0*pi*zones(k)%mesh%z,dBdT,2,Vois_info,zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEp(1,1:Nz)=-reference_parameters%fields%phi0*dPhidT(1,1:Nz)*zones(k)%mesh%Rgeom(1,1:Nz)*&
		zones(k)%mesh%Bphi(1,1:Nz)/(zones(k)%metric_coefficients%Jacobian(1,1:Nz)*zones(k)%mesh%B(1,1:Nz)**(2.d0))
        zones(k)%species(Nspec)%drifts%uEp(0,1:Nz)=zones(k)%species(Nspec)%drifts%uEp(1,1:Nz)	
        !!zones(k)%species(Nspec)%drifts%uEp(0,1:Nz)=-reference_parameters%fields%phi0*dPhidT(0,1:Nz)*zones(k)%mesh%Rgeom(0,1:Nz)*&
	!!	zones(k)%mesh%Bphi(0,1:Nz)/(zones(k)%metric_coefficients%Jacobian(0,1:Nz)*zones(k)%mesh%B(0,1:Nz)**(2.d0))	
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udp(1,1:Nz)=-reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdT(1,1:Nz)*&
		zones(k)%mesh%Rgeom(1,1:Nz)*zones(k)%mesh%Bphi(1,1:Nz)/&
		(zones(k)%species(Nspec)%var(1)%density(1,1:Nz)*zones(k)%metric_coefficients%Jacobian(1,1:Nz)*zones(k)%mesh%B(1,1:Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%udp(0,1:Nz)=zones(k)%species(Nspec)%drifts%udp(1,1:Nz)
        ! Gradient/curvature drift
        zones(k)%species(Nspec)%drifts%uBp(0,1:Nz)=-2.d0*reference_parameters%fields%T0eV/(global_parameters%element_list(Nspec)%Z)*&
		zones(k)%species(Nspec)%var(1)%temperature(1,1:Nz)*dBdT(0,1:Nz)*zones(k)%mesh%Rgeom(0,1:Nz)*zones(k)%mesh%Bphi(0,1:Nz)/&
		(zones(k)%metric_coefficients%Jacobian(0,1:Nz)*zones(k)%mesh%B(0,1:Nz)**(3.d0))
     end if
     ! on the right
     East=Zones(k)%Neighbors(3)
     mEast=Zones(k)%MagNeighbors(3)
     if(mEast.eq.1) then
        call compute_Neighdrift_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,3)
	!!call compute_Neigh_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,3,Vois_info,&
	!!	zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dtheta(Nx,Nz,Pr,2.d0*pi*zones(k)%mesh%z,dPrdT,3)
	call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%B,2.d0*pi*zones(k)%mesh%z,dBdT,3,Vois_info,zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz)=-reference_parameters%fields%phi0*dPhidT(1:Nx,Nz)*zones(k)%mesh%Rgeom(1:Nx,Nz)*&
		zones(k)%mesh%Bphi(1:Nx,Nz)/(zones(k)%metric_coefficients%Jacobian(1:Nx,Nz)*zones(k)%mesh%B(1:Nx,Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz)	
        !!zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz+1)=-reference_parameters%fields%phi0*dPhidT(1:Nx,Nz+1)*zones(k)%mesh%Rgeom(1:Nx,Nz+1)*&
	!!	zones(k)%mesh%Bphi(1:Nx,Nz+1)/(zones(k)%metric_coefficients%Jacobian(1:Nx,Nz+1)*zones(k)%mesh%B(1:Nx,Nz+1)**(2.d0))	
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udp(1:Nx,Nz)=-reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdT(1:Nx,Nz)*&
		zones(k)%mesh%Rgeom(1:Nx,Nz)*zones(k)%mesh%Bphi(1:Nx,Nz)/&
		(zones(k)%species(Nspec)%var(1)%density(1:Nx,Nz)*zones(k)%metric_coefficients%Jacobian(1:Nx,Nz)*zones(k)%mesh%B(1:Nx,Nz)**(2.d0))
	zones(k)%species(Nspec)%drifts%udp(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%udp(1:Nx,Nz)
        ! Gradient/curvature drift
        zones(k)%species(Nspec)%drifts%uBp(1:Nx,Nz+1)=-2.d0*reference_parameters%fields%T0eV/(global_parameters%element_list(Nspec)%Z)*&
		zones(k)%species(Nspec)%var(1)%temperature(1:Nx,Nz)*dBdT(1:Nx,Nz+1)*zones(k)%mesh%Rgeom(1:Nx,Nz+1)*zones(k)%mesh%Bphi(1:Nx,Nz+1)/&
		(zones(k)%metric_coefficients%Jacobian(1:Nx,Nz+1)*zones(k)%mesh%B(1:Nx,Nz+1)**(3.d0))
     end if
     ! on the left
     West=Zones(k)%Neighbors(4)
     mWest=Zones(k)%MagNeighbors(4)
     if(mWest.eq.1) then
        call compute_Neighdrift_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,4)
	!!call compute_Neigh_dtheta(Nx,Nz,zones(k)%electric_fields(1)%phi,2.d0*pi*zones(k)%mesh%z,dPhidT,4,Vois_info,&
	!!	zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%electric_fields(1)%CornersPhi)
        call compute_Neighdrift_dtheta(Nx,Nz,Pr,2.d0*pi*zones(k)%mesh%z,dPrdT,4)
	call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%B,2.d0*pi*zones(k)%mesh%z,dBdT,4,Vois_info,zones(k)%mesh%cornerZ*2.d0*pi,zones(k)%mesh%cornersB)
        ! Electric drift
        zones(k)%species(Nspec)%drifts%uEp(1:Nx,1)=-reference_parameters%fields%phi0*dPhidT(1:Nx,1)*zones(k)%mesh%Rgeom(1:Nx,1)*&
		zones(k)%mesh%Bphi(1:Nx,1)/(zones(k)%metric_coefficients%Jacobian(1:Nx,1)*zones(k)%mesh%B(1:Nx,1)**(2.d0))
	zones(k)%species(Nspec)%drifts%uEp(1:Nx,0)=zones(k)%species(Nspec)%drifts%uEp(1:Nx,1)	
        !!zones(k)%species(Nspec)%drifts%uEp(1:Nx,0)=-reference_parameters%fields%phi0*dPhidT(1:Nx,0)*zones(k)%mesh%Rgeom(1:Nx,0)*&
	!!	zones(k)%mesh%Bphi(1:Nx,0)/(zones(k)%metric_coefficients%Jacobian(1:Nx,0)*zones(k)%mesh%B(1:Nx,0)**(2.d0))	
        ! Diamagnetic drift
        zones(k)%species(Nspec)%drifts%udp(1:Nx,1)=-reference_parameters%fields%T0eV/global_parameters%element_list(Nspec)%Z*dPrdT(1:Nx,1)*&
		zones(k)%mesh%Rgeom(1:Nx,1)*zones(k)%mesh%Bphi(1:Nx,1)/&
		(zones(k)%species(Nspec)%var(1)%density(1:Nx,1)*zones(k)%metric_coefficients%Jacobian(1:Nx,1)*zones(k)%mesh%B(1:Nx,1)**(2.d0))
	zones(k)%species(Nspec)%drifts%udp(1:Nx,0)=zones(k)%species(Nspec)%drifts%udp(1:Nx,1)
        ! Gradient/curvature drift
        zones(k)%species(Nspec)%drifts%uBp(1:Nx,0)=-2.d0*reference_parameters%fields%T0eV/(global_parameters%element_list(Nspec)%Z)*&
		zones(k)%species(Nspec)%var(1)%temperature(1:Nx,0)*dBdT(1:Nx,0)*zones(k)%mesh%Rgeom(1:Nx,0)*zones(k)%mesh%Bphi(1:Nx,0)/&
		(zones(k)%metric_coefficients%Jacobian(1:Nx,0)*zones(k)%mesh%B(1:Nx,0)**(3.d0))
     end if
  deallocate(Pr,dPhidT,dPrdT,dBdT)
  end do
end subroutine gradient_BC_dtheta

  subroutine compute_Neighdrift_dtheta(Nx,Nz,X,T,dXdT,Neigh)
    implicit none
    integer*4 Nx,Nz
    real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: X
    real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: T
    real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdT
    integer*4 Neigh !1N, 2S, 3E, 4W
    real*8 pi
    pi=4._8*atan(1._8)
    select case(Neigh)
    case(1) ! on top
       dXdT(Nx,2:Nz-1)=(X(Nx,3:Nz)*(T(Nx,2:Nz-1)-T(Nx,1:Nz-2))**2.d0&
            -X(Nx,1:Nz-2)*(T(Nx,3:Nz)-T(Nx,2:Nz-1))**2.d0&
            +X(Nx,2:Nz-1)*((T(Nx,3:Nz)-T(Nx,2:Nz-1))**2.d0&
            -(T(Nx,2:Nz-1)-T(Nx,1:Nz-2))**2.d0))&
            /((T(Nx,3:Nz)-T(Nx,2:Nz-1))*(T(Nx,2:Nz-1)-T(Nx,1:Nz-2))&
            *(T(Nx,3:Nz)-T(Nx,1:Nz-2)))
       !dXdT(Nx,2:Nz-1) = (X(Nx,3:Nz)-X(Nx,1:Nz-2))/(T(Nx,3:Nz)-T(Nx,1:Nz-2))
       dXdT(Nx,1)=(X(Nx,2)-X(Nx,1))/((T(Nx,2)-T(Nx,1)))
       dXdT(Nx,Nz)=(X(Nx,Nz)-X(Nx,Nz-1))/((T(Nx,Nz)-T(Nx,Nz-1)))
    case(2) ! at the bottom
       dXdT(1,2:Nz-1)=(X(1,3:Nz)*(T(1,2:Nz-1)-T(1,1:Nz-2))**2.d0&
            -X(1,1:Nz-2)*(T(1,3:Nz)-T(1,2:Nz-1))**2.d0&
            +X(1,2:Nz-1)*((T(1,3:Nz)-T(1,2:Nz-1))**2.d0&
            -(T(1,2:Nz-1)-T(1,1:Nz-2))**2.d0))&
            /((T(1,3:Nz)-T(1,2:Nz-1))*(T(1,2:Nz-1)-T(1,1:Nz-2))&
           *(T(1,3:Nz)-T(1,1:Nz-2)))
       !dXdT(1,2:Nz-1) = (X(1,3:Nz)-X(1,1:Nz-2))/(T(1,3:Nz)-T(1,1:Nz-2))
       dXdT(1,1)=(X(1,2)-X(1,1))/((T(1,2)-T(1,1)))
       dXdT(1,Nz)=(X(1,Nz)-X(1,Nz-1))/((T(1,Nz)-T(1,Nz-1)))
    case(3) ! on the right
       dXdT(1:Nx,Nz)=(X(1:Nx,Nz)-X(1:Nx,Nz-1))/(T(1:Nx,Nz)-T(1:Nx,Nz-1))
    case(4) ! on the left
       dXdT(1:Nx,1)=(X(1:Nx,2)-X(1:Nx,1))/(T(1:Nx,2)-T(1:Nx,1))
    case default
    end select
  end subroutine compute_Neighdrift_dtheta




