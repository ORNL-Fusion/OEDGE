subroutine current_balance_custom(zone,i)
  use all_variables, only : global_parameters, global_variables
  use MZone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: i
  integer*4 :: j,Nx,Nz
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  Fluxes=-zone%electric_fields(1)%j_perp&
       +zone%electric_fields(1)%j_parallel&
       +zone%electric_fields(1)%j_para_adv_W*global_variables%dt/global_variables%dt_vort&
       +zone%electric_fields(1)%j_diff_W*global_variables%dt/global_variables%dt_vort
  Source=divergence(zone,Fluxes,Nx,Nz)
  open(unit=101,file='north_j',status='unknown')
  open(unit=102,file='south_j',status='unknown')
  open(unit=103,file='east_j',status='unknown')
  open(unit=104,file='west_j',status='unknown')
  open(unit=105,file='div_j',status='unknown')
  do j=1,Nz
     write(101,100) Fluxes(i,j,1)*(1.d0-zone%masks%chi2(i,j))
     write(102,100) Fluxes(i,j,2)*(1.d0-zone%masks%chi2(i,j))
     write(103,100) Fluxes(i,j,3)*(1.d0-zone%masks%chi2(i,j))
     write(104,100) Fluxes(i,j,4)*(1.d0-zone%masks%chi2(i,j))
     write(105,100) Source(i,j)*(1.d0-zone%masks%chi2(i,j))
  end do
  close(101)
  close(102)
  close(103)
  close(104)
  close(105)
100 format(512es15.7)
end subroutine current_balance_custom
