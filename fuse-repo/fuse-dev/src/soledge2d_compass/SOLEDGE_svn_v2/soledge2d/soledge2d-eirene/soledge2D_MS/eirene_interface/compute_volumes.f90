subroutine compute_volumes(voltot,volumes,NZones)
  use all_variables, only : global_parameters, zones
  use Meirene_vars
  implicit none
  integer*4,intent(in) :: Nzones
  integer*4 i,j,k
  integer*4 Nx,Nz
  real*8 pi
  real*8, intent(out) :: voltot
  Type(volume) :: volumes(NZones)
  pi=4.D0*atan(1.D0)

  voltot=0.D0

  do k=1,NZones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(volumes(k)%cell(0:Nx+1,0:Nz+1))
     volumes(k)%cell=0.d0
     do i=1,Nx
        do j=1,Nz
           if(Zones(k)%masks%chi2(i,j).eq.0) then
              volumes(k)%cell(i,j) = zones(k)%metric_coefficients%dvol_PU(i,j)
              voltot=voltot + volumes(k)%cell(i,j)          
           end if
        end do
     end do
  end do

end subroutine compute_volumes
