subroutine current_balance(zone)
  use all_variables, only : global_parameters, global_variables, flags, drift_flags
  use MZone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: i,j,Nx,Nz
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  Fluxes=-zone%electric_fields(1)%j_perp&
       +zone%electric_fields(1)%j_parallel&
       +zone%electric_fields(1)%j_para_adv_W*global_variables%dt/global_variables%dt_vort_old&
       +zone%electric_fields(1)%j_diff_W*global_variables%dt/global_variables%dt_vort_old
  if(drift_flags%solve_drift) then
	Fluxes = Fluxes + zone%species(1)%drifts%jExB*global_variables%dt/global_variables%dt_vort_old +&
		 zone%species(1)%drifts%jdiam*global_variables%dt/global_variables%dt_vort_old +&
		 zone%species(1)%drifts%jBxDB*global_variables%dt/global_variables%dt_vort_old
  end if
  Source=divergence(zone,Fluxes,Nx,Nz)
  do i=1,Nx
     write(100+zone%number,100) (Fluxes(i,j,1),j=1,Nz)
     write(200+zone%number,100) (Fluxes(i,j,2),j=1,Nz)
     write(300+zone%number,100) (Fluxes(i,j,3),j=1,Nz)
     write(400+zone%number,100) (Fluxes(i,j,4),j=1,Nz)
     write(500+zone%number,100) (Source(i,j),j=1,Nz)
  end do
100 format(512es15.7)
  deallocate(Fluxes,Source)
end subroutine current_balance
