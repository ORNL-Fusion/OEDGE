subroutine compute_interchange(zone)
  use MZone
  use all_variables, only : reference_parameters, kepsilon_param, global_variables
  implicit none
  Type(Tzone) :: zone
  integer*4 :: i,j,Nx,Nz
  real*8,allocatable :: p_i(:,:)
  real*8,allocatable :: n(:,:)
  real*8,allocatable :: Ti(:,:)
  real*8,allocatable :: Te(:,:)
  real*8 :: temp, temp2
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(p_i(0:Nx+1,0:Nz+1))
  allocate(Ti(0:Nx+1,0:Nz+1))
  allocate(Te(0:Nx+1,0:Nz+1))
  allocate(n(0:Nx+1,0:Nz+1))
  p_i=zone%species(0)%var(1)%density*zone%species(1)%var(1)%temperature
  n=zone%species(0)%var(1)%density
  Ti=zone%species(1)%var(1)%temperature
  Te=zone%species(0)%var(1)%temperature
  do i=1,Nx
     do j=1,Nz
        temp=(n(i,j+1)-n(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *(zone%mesh%B(i,j+1)-zone%mesh%B(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *zone%metric_coefficients%ctt(i,j)
        temp=temp+&
             (n(i+1,j)-n(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *(zone%mesh%B(i+1,j)-zone%mesh%B(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *zone%metric_coefficients%cpp(i,j)
        temp=temp+&
             (n(i+1,j)-n(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *(zone%mesh%B(i,j+1)-zone%mesh%B(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *zone%metric_coefficients%cpt(i,j)
        temp=temp+&
             (n(i,j+1)-n(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *(zone%mesh%B(i+1,j)-zone%mesh%B(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *zone%metric_coefficients%cpt(i,j)
        temp=(temp)&
             /(n(i,j)+1.d-4)/zone%mesh%B(i,j)

        temp2=(Ti(i,j+1)-Ti(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *(zone%mesh%B(i,j+1)-zone%mesh%B(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *zone%metric_coefficients%ctt(i,j)
        temp2=temp2+&
             (Ti(i+1,j)-Ti(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *(zone%mesh%B(i+1,j)-zone%mesh%B(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *zone%metric_coefficients%cpp(i,j)
        temp2=temp2+&
             (Ti(i+1,j)-Ti(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *(zone%mesh%B(i,j+1)-zone%mesh%B(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *zone%metric_coefficients%cpt(i,j)
        temp2=temp2+&
             (Ti(i,j+1)-Ti(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *(zone%mesh%B(i+1,j)-zone%mesh%B(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *zone%metric_coefficients%cpt(i,j)
        temp2=(temp2)&
             /(Ti(i,j)+1.d-4)/zone%mesh%B(i,j)

        zone%kepsilon(1)%interchange(i,j)=((temp+kepsilon_param%gradT_weight*temp2)&
             *zone%mesh%Rgeom(i,j)**2/reference_parameters%geometry%rs0**2&
             -kepsilon_param%th_interchange*(1.d0+Ti(i,j)/Te(i,j)))/zone%mesh%Rgeom(i,j)**2
        zone%kepsilon(1)%interchange(i,j)=sign(sqrt(abs(zone%kepsilon(1)%interchange(i,j))),zone%kepsilon(1)%interchange(i,j))
        zone%kepsilon(1)%interchange(i,j)=zone%kepsilon(1)%interchange(i,j)*sqrt((zone%species(0)%var(1)%temperature(i,j)+&
             zone%species(1)%var(1)%temperature(i,j))/zone%species(1)%element%mass)
        zone%kepsilon(1)%interchange(i,j)=zone%kepsilon(1)%interchange(i,j)*kepsilon_param%C_interchange
        zone%kepsilon(1)%interchange(i,j)=max(zone%kepsilon(1)%interchange(i,j),-1.D0/global_variables%dt*0.5)
        zone%kepsilon(1)%interchange(i,j)=min(zone%kepsilon(1)%interchange(i,j),1.D0/global_variables%dt*0.5)
        zone%kepsilon(1)%interchange(i,j)=zone%kepsilon(1)%interchange(i,j)*(1.D0-zone%masks%chi2(i,j))
        if((zone%masks%chi2(i,j+1).eq.1).or.(zone%masks%chi2(i,j-1).eq.1).or.(zone%masks%chi2(i+1,j).eq.1).or.(zone%masks%chi2(i-1,j).eq.1)) then
           zone%kepsilon(1)%interchange(i,j)=0.D0
        end if
     end do
  end do
  if(zone%Neighbors(1).lt.0) then
     !     zone%kepsilon(1)%interchange(Nx,:)=0.D0
     zone%kepsilon(1)%interchange(Nx,:)=zone%kepsilon(1)%interchange(Nx-1,:)
  end if
  if(zone%Neighbors(2).lt.0) then
     !     zone%kepsilon(1)%interchange(1,:)=0.D0
     zone%kepsilon(1)%interchange(1,:)=zone%kepsilon(1)%interchange(2,:)
  end if
  if(zone%Neighbors(3).lt.0) then
     !     zone%kepsilon(1)%interchange(:,Nz)=0.D0
     zone%kepsilon(1)%interchange(:,Nz)=zone%kepsilon(1)%interchange(:,Nz-1)
  end if
  if(zone%Neighbors(4).lt.0) then
     !     zone%kepsilon(1)%interchange(:,1)=0.D0
     zone%kepsilon(1)%interchange(:,1)=zone%kepsilon(1)%interchange(:,2)
  end if
  deallocate(p_i,n,Ti,Te)


  zone%kepsilon(1)%interchange(:,1)=max(zone%kepsilon(1)%interchange(:,1),0.d0)

end subroutine compute_interchange
