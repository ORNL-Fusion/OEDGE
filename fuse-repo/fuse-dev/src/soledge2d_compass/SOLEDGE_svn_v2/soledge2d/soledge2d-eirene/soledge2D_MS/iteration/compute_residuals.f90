subroutine compute_residuals()
  use all_variables, only : global_parameters, zones, global_variables
  implicit none
  integer*4 :: n,k
  real*8,allocatable :: resn(:),resG(:),resT(:)
  allocate(resn(1:global_parameters%N_zones))
  allocate(resG(1:global_parameters%N_zones))
  allocate(resT(1:global_parameters%N_zones))
  resn=0.d0
  resG=0.d0
  resT=0.d0
  do n=0,global_parameters%N_ions
     do k=1,global_parameters%N_zones
        resn(k)=zones(k)%species(n)%residuals%resn
        resG(k)=zones(k)%species(n)%residuals%resG
        resT(k)=zones(k)%species(n)%residuals%resT
     end do
     global_variables%residuals(n)%resn=maxval(resn)
     global_variables%residuals(n)%resG=maxval(resG)
     global_variables%residuals(n)%resT=maxval(resT)
  end do
  deallocate(resn,resG,resT)
end subroutine compute_residuals
