subroutine find_delta_r()
  use all_variables, only :global_parameters, zones
  use Msmoothing_vars
  implicit none
  integer*4 :: k,i,j,Nx,Nz
  real*8 :: S
  integer*4 :: Nx_N, Nz_N
  delta_r=1.d10
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(zones(k)%metric_coefficients%delta_r(0:Nx+1,0:Nz+1))
     do i=1,Nx
        do j=1,Nz
           delta_r=min(sqrt(abs((zones(k)%mesh%Rcorner(i,j+1)-zones(k)%mesh%Rcorner(i,j))&
                *(zones(k)%mesh%Zcorner(i+1,j)-zones(k)%mesh%Zcorner(i,j))&
                -(zones(k)%mesh%Zcorner(i,j+1)-zones(k)%mesh%Zcorner(i,j))&
                *(zones(k)%mesh%Rcorner(i+1,j)-zones(k)%mesh%Rcorner(i,j)))),delta_r)
           zones(k)%metric_coefficients%delta_r(i,j)=sqrt((zones(k)%mesh%Rcorner(i,j+1)-zones(k)%mesh%Rcorner(i,j))**2&
                +(zones(k)%mesh%Zcorner(i,j+1)-zones(k)%mesh%Zcorner(i,j))**2)
           zones(k)%metric_coefficients%delta_r(i,j)=min(zones(k)%metric_coefficients%delta_r(i,j),&
                sqrt((zones(k)%mesh%Rcorner(i+1,j)-zones(k)%mesh%Rcorner(i,j))**2&
                +(zones(k)%mesh%Zcorner(i+1,j)-zones(k)%mesh%Zcorner(i,j))**2))
        end do
     end do
  end do
  ! Broadcast
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     !North
     if(zones(k)%Neighbors(1).gt.0) then
        do j=1,Nz
           zones(k)%metric_coefficients%delta_r(Nx+1,j)=zones(zones(k)%Neighbors(1))%metric_coefficients%delta_r(1,j)
        end do
     else
        do j=1,Nz
           zones(k)%metric_coefficients%delta_r(Nx+1,j)=zones(k)%metric_coefficients%delta_r(Nx,j)
        end do
     end if
     !South
     if(zones(k)%Neighbors(2).gt.0) then
        Nx_N=zones(zones(k)%Neighbors(2))%mesh%Nx
        do j=1,Nz
           zones(k)%metric_coefficients%delta_r(0,j)=zones(zones(k)%Neighbors(2))%metric_coefficients%delta_r(Nx_N,j)
        end do
     else
        do j=1,Nz
           zones(k)%metric_coefficients%delta_r(0,j)=zones(k)%metric_coefficients%delta_r(1,j)
        end do
     end if
     !East
     if(zones(k)%Neighbors(3).gt.0) then
        do i=1,Nx
           zones(k)%metric_coefficients%delta_r(i,Nz+1)=zones(zones(k)%Neighbors(3))%metric_coefficients%delta_r(i,1)
        end do
     else
        do i=1,Nx
           zones(k)%metric_coefficients%delta_r(i,Nz+1)=zones(k)%metric_coefficients%delta_r(i,Nz)
        end do
     end if
     !West
     if(zones(k)%Neighbors(4).gt.0) then
        Nz_N=zones(zones(k)%Neighbors(4))%mesh%Nz
        do i=1,Nx
           zones(k)%metric_coefficients%delta_r(i,0)=zones(zones(k)%Neighbors(4))%metric_coefficients%delta_r(i,Nz_N)
        end do
     else
        do i=1,Nx
           zones(k)%metric_coefficients%delta_r(i,0)=zones(k)%metric_coefficients%delta_r(i,1)
        end do
     end if
     !Corners
     zones(k)%metric_coefficients%delta_r(0,0)=0.5d0*(zones(k)%metric_coefficients%delta_r(1,0)+zones(k)%metric_coefficients%delta_r(0,1))
     zones(k)%metric_coefficients%delta_r(0,Nz+1)=0.5d0*(zones(k)%metric_coefficients%delta_r(1,Nz+1)+zones(k)%metric_coefficients%delta_r(0,Nz))
     zones(k)%metric_coefficients%delta_r(Nx+1,0)=0.5d0*(zones(k)%metric_coefficients%delta_r(Nx,0)+zones(k)%metric_coefficients%delta_r(Nx+1,1))
     zones(k)%metric_coefficients%delta_r(Nx+1,Nz+1)=0.5d0*(zones(k)%metric_coefficients%delta_r(Nx,Nz+1)+zones(k)%metric_coefficients%delta_r(Nx+1,Nz))
  end do
end subroutine find_delta_r
