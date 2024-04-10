subroutine styx_feed_plasma_to_eirene
  use all_variables, only : interp_data2, reference_parameters, global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use styx2eirene
  implicit none

  integer :: ipls

  ! first update plasma parameters (interpolation, then feed into eirene)

  call styx_vertices_to_centre_minT(Interp_data2%knots_temperature(:,0,:),eiv_e%T) 
  eiv_e%T = eiv_e%T*reference_parameters%fields%T0eV
  
!  do ipls=1,nplsi ! nplsi does not work if sputtered species not followed
   do ipls=1,global_parameters%n_ions
    call styx_vertices_to_centre(Interp_data2%knots_density(:,ipls,:),eiv_i(ipls)%dens)
    call styx_vertices_to_centre(Interp_data2%knots_density(:,ipls,:)*Interp_data2%knots_velocity(:,ipls,:),eiv_i(ipls)%parmom)
    call styx_vertices_to_centre_minT(Interp_data2%knots_temperature(:,ipls,:),eiv_i(ipls)%T)
    
    ! renormalize to SI units (then EIRENE units in styx_send_plasma)

    eiv_i(ipls)%dens = eiv_i(ipls)%dens*reference_parameters%fields%n0
    eiv_i(ipls)%parmom = eiv_i(ipls)%parmom*reference_parameters%fields%n0&
          *reference_parameters%fields%c0
    
    where (eiv_i(ipls)%dens>0._dp)
      vpar_tri(:,ipls)=eiv_i(ipls)%parmom/eiv_i(ipls)%dens
    elsewhere
       vpar_tri(:,ipls)=0._dp
    end where

    eiv_i(ipls)%vx=vpar_tri(:,ipls)*BX_tri
    eiv_i(ipls)%vy=vpar_tri(:,ipls)*BY_tri
    eiv_i(ipls)%vz=vpar_tri(:,ipls)*BZ_tri

    eiv_i(ipls)%T = eiv_i(ipls)%T*reference_parameters%fields%T0eV
  enddo
  
  ! conversion to cgs units for EIRENE

  call eirene_if1copls

end subroutine styx_feed_plasma_to_eirene
