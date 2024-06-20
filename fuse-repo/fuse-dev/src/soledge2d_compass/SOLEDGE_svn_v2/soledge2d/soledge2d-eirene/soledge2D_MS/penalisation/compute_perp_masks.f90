subroutine compute_perp_masks()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     zones(k)%masks%chi5=0.d0
     zones(k)%masks%chi6=0.d0
     do i=1,Nx
        do j=1,Nz
           if((zones(k)%masks%chi2(i,j).eq.1).and.(zones(k)%masks%chi2(i+1,j).eq.0)) then
              zones(k)%masks%chi5(i,j)=1.d0
           end if
           if((zones(k)%masks%chi2(i,j).eq.1).and.(zones(k)%masks%chi2(i-1,j).eq.0)) then
              zones(k)%masks%chi6(i,j)=1.d0
           end if
           if((zones(k)%masks%chi1(i,j).eq.1.d0).or.(zones(k)%masks%chi3(i,j).eq.1.d0)) then
              zones(k)%masks%chi5(i,j)=0.d0
              zones(k)%masks%chi6(i,j)=0.d0
           end if
        end do
     end do
  end do
end subroutine compute_perp_masks
