subroutine compute_kh(zone)
  use MZone
  use all_variables, only : reference_parameters, kepsilon_param, global_variables
  implicit none
  Type(Tzone) :: zone
  integer*4 :: i,j,Nx,Nz
  real*8,allocatable :: v(:,:)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(v(0:Nx+1,0:Nz+1))
  v=zone%species(0)%var(1)%velocity*sqrt(zone%mesh%Br**2+zone%mesh%Bz**2)/zone%mesh%B
  do i=1,Nx
     do j=1,Nz
        zone%kepsilon(1)%kh(i,j)=(v(i,j+1)-v(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *(v(i,j+1)-v(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *zone%metric_coefficients%ctt(i,j)
        zone%kepsilon(1)%kh(i,j)=zone%kepsilon(1)%kh(i,j)+&
             (v(i+1,j)-v(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *(v(i+1,j)-v(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *zone%metric_coefficients%cpp(i,j)
        zone%kepsilon(1)%kh(i,j)=zone%kepsilon(1)%kh(i,j)+&
             (v(i+1,j)-v(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *(v(i,j+1)-v(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *zone%metric_coefficients%cpt(i,j)
        zone%kepsilon(1)%kh(i,j)=zone%kepsilon(1)%kh(i,j)+&
             (v(i,j+1)-v(i,j-1))/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))&
             *(v(i+1,j)-v(i-1,j))/(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
             *zone%metric_coefficients%cpt(i,j)
        zone%kepsilon(1)%kh(i,j)=0.5D0*zone%kepsilon(1)%kh(i,j)*zone%species(0)%var(1)%density(i,j)&
             *zone%kepsilon(1)%mu_t(i,j)
        Zone%kepsilon(1)%kh(i,j)=zone%kepsilon(1)%kh(i,j)*reference_parameters%fields%c0**2/&
             reference_parameters%geometry%rs0**2*reference_parameters%fields%tau0/reference_parameters%fields%k0&
             *(reference_parameters%fields%k0**2/reference_parameters%fields%epsilon0)
        Zone%kepsilon(1)%kh(i,j)=Zone%kepsilon(1)%kh(i,j)*(1.D0-zone%masks%chi2(i,j))
        if((zone%masks%chi2(i,j+1).eq.1).or.(zone%masks%chi2(i,j-1).eq.1).or.(zone%masks%chi2(i+1,j).eq.1).or.(zone%masks%chi2(i-1,j).eq.1)) then
           Zone%kepsilon(1)%kh(i,j)=0.D0
        end if
        Zone%kepsilon(1)%kh(i,j)=Zone%kepsilon(1)%kh(i,j)*kepsilon_param%C_kh
     End do
  end do
  if(zone%Neighbors(1).lt.0) then
     zone%kepsilon(1)%kh(Nx,:)=0.D0
  end if
  if(zone%Neighbors(2).lt.0) then
     zone%kepsilon(1)%kh(1,:)=0.D0
  end if
  if(zone%Neighbors(3).lt.0) then
     zone%kepsilon(1)%kh(:,Nz)=0.D0
  end if
  if(zone%Neighbors(4).lt.0) then
     zone%kepsilon(1)%kh(:,1)=0.D0
  end if
  deallocate(v)
end subroutine compute_kh
