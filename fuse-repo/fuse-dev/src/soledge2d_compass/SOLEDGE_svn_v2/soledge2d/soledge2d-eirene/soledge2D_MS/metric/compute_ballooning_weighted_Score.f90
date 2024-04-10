subroutine compute_ballooning_weighted_Score()
  use all_variables, only : zones, global_parameters,&
       reference_parameters, ballooning_parameters
  use Mphysics
  implicit none
  integer*4 :: j,k
  integer*4 :: Nx,Nz
  real*8 :: Score
  Score=0.
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     if(zones(k)%Neighbors(1).eq.-1) then
        do j=1,Nz
           Score=Score+zones(k)%metric_coefficients%dS_North_PU(Nx,j)&
                *zones(k)%species(1)%transport_perp%D_p(Nx,j)
        end do
     end if
     if(zones(k)%Neighbors(2).eq.-1) then
        do j=1,Nz
           Score=Score+zones(k)%metric_coefficients%dS_South_PU(1,j)&
                *zones(k)%species(1)%transport_perp%D_p(1,j)
        end do
     end if
  end do
  ballooning_parameters%core_weighted_integral=Score
end subroutine compute_ballooning_weighted_Score
