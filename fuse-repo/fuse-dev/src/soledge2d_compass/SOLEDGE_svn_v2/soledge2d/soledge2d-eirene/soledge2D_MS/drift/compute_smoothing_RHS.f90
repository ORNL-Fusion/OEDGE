subroutine compute_smoothing_rhs()
  use all_variables, only : global_parameters, zones, reference_parameters
  use Msmoothing_vars
  use Mlist
  use Mphysics
  use Mpastix_solve
  implicit none
  integer*4 i,j,k
  integer*4 Nx,Nz

  CSC_smoothing%b=0.D0
  do k=1, global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           if(zones(k)%masks%chi2(i,j).eq.0) then
              CSC_smoothing%b(Zones(k)%mesh%index(i,j))=zones(k)%electric_fields(1)%RHS(i,j)
           end if
        end do
     end do
     !zero for the rest should be fine (no need to do anything else)
  end do

end subroutine compute_smoothing_rhs
