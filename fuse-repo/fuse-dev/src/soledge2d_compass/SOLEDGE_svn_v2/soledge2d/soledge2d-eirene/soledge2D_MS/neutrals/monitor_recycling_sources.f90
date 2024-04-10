subroutine monitor_recycling_sources()
  use all_variables, only : zones, reference_parameters, global_parameters, globals
  Implicit none
  integer*4 :: k,i,j
  real*8 :: Sneutrals, Sionization
  Sneutrals=0.D0
  Sionization=0.D0
  do k=1,global_parameters%N_Zones
     do i=1,zones(k)%mesh%Nx
        do j=1,zones(k)%mesh%Nz
           Sneutrals=Sneutrals+zones(k)%neutrals%Sn_nn(i,j)*zones(k)%metric_coefficients%Dvol_PU(i,j)&
                *reference_parameters%fields%n0/reference_parameters%fields%tau0
           Sionization=Sionization+zones(k)%species(1)%sources%Sn_n(i,j)*zones(k)%metric_coefficients%Dvol_PU(i,j)&
                *reference_parameters%fields%n0/reference_parameters%fields%tau0
        end do
     end do
  end do
  write(*,*) 'Neutral source = ', Sneutrals
  write(*,*) 'Ionization source = ', Sionization
  write(*,*) 'ion outflux =',  globals%flux_tot_out_ac(1)
  write(*,*) 'Recycling coefficient =',  Sneutrals/globals%flux_tot_out_ac(1)
  write(*,*) 'Effective recycling coefficient =',  Sionization/globals%flux_tot_out_ac(1)
end subroutine monitor_recycling_sources
