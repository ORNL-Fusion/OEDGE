subroutine monitor_fluxes_and_sources(n_ite)
  use all_variables, only : zones, global_parameters, globals, flags
  implicit none
  integer*4,intent(in) :: n_ite
  integer*4 :: k,nion
  real*8 :: Flux_sum,S_sum
  real*8 :: Flux_loc,S_loc
  character(50) :: filename
  if((modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0)&
       .or.(global_parameters%N_iterations.lt.1000)) then
     call compute_total_radiation()
     do nion=0,global_parameters%N_ions
        !#########################################################
        !##    Particle flux at core interface  ##################
        Flux_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_particle_flux_core(zones(k),nion,Flux_loc)
           Flux_sum=Flux_sum+Flux_loc
        end do
       globals%flux_tot_in_ac(nion)=Flux_sum
        !#########################################################
        !##    Particle flux at wall interface  ##################
        Flux_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_particle_flux_wall(zones(k),nion,Flux_loc)
           Flux_sum=Flux_sum+Flux_loc
        end do
        globals%flux_tot_out_ac(nion)=Flux_sum
        !#########################################################
        !##    Particle source                  ##################
        S_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_particle_source(zones(k),nion,S_loc)
           S_sum=S_sum+S_loc
        end do
        globals%source_n(nion)=S_sum
        !#########################################################
        !#     Particle variations                ##################
        S_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_particle_variation(zones(k),nion,S_loc)
           S_sum=S_sum+S_loc
        end do
        globals%variation_N(nion)=S_sum
        !#########################################################
        !#     Particle stored                    ##################
        S_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_particle_stored(zones(k),nion,S_loc)
           S_sum=S_sum+S_loc
        end do
        globals%stored_N(nion)=S_sum

        !#########################################################
        !##      Energy flux at core interface  ##################
        Flux_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_energy_flux_core(zones(k),nion,Flux_loc)
           Flux_sum=Flux_sum+Flux_loc
        end do
        globals%flux_totE_in_ac(nion)=Flux_sum
        !#########################################################
        !##      Energy flux at wall interface  ##################
        Flux_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_energy_flux_wall(zones(k),nion,Flux_loc)
           Flux_sum=Flux_sum+Flux_loc
        end do
        globals%flux_totE_out_ac(nion)=Flux_sum
        !#########################################################
        !##    Energy source                    ##################
        S_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_energy_source(zones(k),nion,S_loc)
           S_sum=S_sum+S_loc
        end do
        globals%source_E(nion)=S_sum
        !#########################################################
        !#     Energy variations                ##################
        S_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_energy_variation(zones(k),nion,S_loc)
           S_sum=S_sum+S_loc
        end do
        globals%variation_E(nion)=S_sum
        !#########################################################
        !#     Energy stored                    ##################
        S_sum=0.d0
        do k=1,global_parameters%N_zones
           call compute_energy_stored(zones(k),nion,S_loc)
           S_sum=S_sum+S_loc
        end do
        globals%stored_E(nion)=S_sum


        write(filename,"(A9,I0)") "balances_",nion
        open(unit=100,file=trim(adjustl(filename)),status='unknown',access='append')
        write(100,50) globals%flux_tot_in_ac(nion), globals%flux_tot_out_ac(nion), globals%source_n(nion), &
             globals%flux_totE_in_ac(nion), globals%flux_totE_out_ac(nion), globals%source_E(nion)&
             , globals%total_radiation(nion), globals%variation_E(nion), globals%stored_E(nion)&
             , globals%variation_N(nion), globals%stored_N(nion)
50      format(512es15.7)
        close(100)

        if(flags%neutral_model.eq.2) then
           call monitor_recycling_sources()
        end if

     end do
  end if
end subroutine monitor_fluxes_and_sources
