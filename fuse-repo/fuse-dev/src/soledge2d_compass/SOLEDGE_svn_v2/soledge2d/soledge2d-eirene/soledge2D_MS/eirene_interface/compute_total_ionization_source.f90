subroutine compute_total_ionization_source()
  use all_variables, only : global_parameters, zones, globals, reference_parameters
  use Mphysics
  implicit none
  real*8,allocatable :: source_ionz_tot(:), source_E_ionz(:)
  integer*4 :: i,j,k,Nx,Nz,nion,ind
  real*8 :: Vol
  allocate(source_ionz_tot(0:global_parameters%N_species))
  allocate(source_E_ionz(0:global_parameters%N_species))
  source_ionz_tot = 0.D0
  source_E_ionz = 0.D0
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,Nx,Nz,i,j,Vol,nion,ind)
  !$OMP DO REDUCTION(+:source_ionz_tot,source_E_ionz)
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           Vol=zones(k)%metric_coefficients%dvol_PU(i,j)
           if(Zones(k)%masks%chi2(i,j).eq.0) then
              do nion=1,global_parameters%N_species
                 ind=global_parameters%ind_ion_0(nion)
                 source_ionz_tot(nion)=source_ionz_tot(nion)+&
                      Zones(k)%species(ind)%sources%Sn_n(i,j)*Vol*reference_parameters%fields%n0/reference_parameters%fields%tau0
                 source_E_ionz(nion)=source_E_ionz(nion)+&
                      Zones(k)%species(ind)%sources%Sn_E(i,j)*Vol*reference_parameters%fields%n0&
                      *reference_parameters%fields%T0*kb/reference_parameters%fields%tau0
              end do
              source_E_ionz(0)=source_E_ionz(0)+&
                   Zones(k)%species(0)%sources%Sn_E(i,j)*Vol*reference_parameters%fields%n0&
                   *reference_parameters%fields%T0*kb/reference_parameters%fields%tau0
           end if
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  do nion=1,global_parameters%N_species
     globals%source_ionz_tot(nion)=source_ionz_tot(nion)
     globals%source_E_ionz(nion)=source_E_ionz(nion)
  end do
  globals%source_E_ionz(0)=source_E_ionz(0)
end subroutine compute_total_ionization_source
