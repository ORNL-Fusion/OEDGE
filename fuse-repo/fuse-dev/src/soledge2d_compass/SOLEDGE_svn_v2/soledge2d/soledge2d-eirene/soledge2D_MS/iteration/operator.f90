module MOperator

contains

  function divergence(zone,flux,Nx,Nz)
    use all_variables, only : reference_parameters
    use MZone
    use Mphysics
    implicit none
    type(TZone),intent(in) :: zone
    integer*4,intent(in) :: Nx,Nz
    real*8,intent(in) :: flux(Nx,Nz,4)
    real*8 :: divergence(Nx,Nz)
    integer*4 :: i,j
    real*8 :: rs0,DVOL
    rs0=reference_parameters%geometry%rs0
    do i=1,Nx
       do j=1,Nz
          DVOL=zone%metric_coefficients%dvol_DD(i,j)
          divergence(i,j)=flux(i,j,3)&
               *zone%metric_coefficients%ds_east_DD(i,j)&
               -flux(i,j,4)&
               *zone%metric_coefficients%ds_west_DD(i,j)
          divergence(i,j)=divergence(i,j)&
               +(flux(i,j,1)&
               *zone%metric_coefficients%ds_north_DD(i,j)&
               -flux(i,j,2)&
               *zone%metric_coefficients%ds_south_DD(i,j))     
          divergence(i,j)=divergence(i,j)/DVOL
       end do
    end do
  end function divergence

end module MOperator
