subroutine styx_feed_wall_parameters_to_eirene
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_clgin
  use eirmod_ccona
  use styx2eirene
  use all_variables, only : global_parameters

  implicit none
  integer :: isp,ispz,ipfc,ipump

  real*8 :: tempR

  ! set recycling fluxs (in the future may evolve from call to call)                                                                            
  if (any(pumps(:)%isSpeedSet)) call styx_compute_pump_species_temperature()

  do ispz=1,nspz
     isp=map_species_index(ispz)
     write(*,*) ispz,isp
     if (isp >0) then
        do ipfc=1,n_pfc_types
           RECYCT(ispz,NLIM+1+n_pumps+ipfc)=pfc_models(ipfc)%R(isp)
        enddo
        do ipump=1,n_pumps
              if(pumps(ipump)%isSpeedSet.and.(isp==1)) then
                 if(pumps(ipump)%species_temperature(ispz)>0) then
                    tempR=pumps(ipump)%pumping_Speed*1000.d0&
                         /(pumps(ipump)%Surface*10000.d0*3.638d0&
                         *sqrt(pumps(ipump)%species_temperature(ispz)/&
                         global_parameters%element_list(isp)%mass))
                    RECYCT(ispz,NLIM+1+ipump)=min(max(1.d0-tempR,0.d0),1.d0)
                 else
                    RECYCT(ispz,NLIM+1+ipump)=0.d0
                 end if
              else
                 RECYCT(ispz,NLIM+1+ipump)=pumps(ipump)%R(isp)
              end if
              write(560,*) ipump, ispz,  RECYCT(ispz,NLIM+1+ipump)
        enddo
     else
        ! sputtered species not in Soledge2D, assume averything sticks                                                                            
        RECYCT(ispz,NLIM+1:NLIM+1+n_pumps+n_pfc_types)=0._dp
        RECYCT(ispz,NLIM+1:NLIM+1+n_pumps+n_pfc_types)=0._dp
     endif
  enddo
  ! set wall temperature                                                                                                                        
  do ipfc=1,n_pfc_types
     EWALL(NLIM+1+n_pumps+ipfc)=-1._dp*pfc_models(ipfc)%T*EVKEL
  enddo

end subroutine styx_feed_wall_parameters_to_eirene

