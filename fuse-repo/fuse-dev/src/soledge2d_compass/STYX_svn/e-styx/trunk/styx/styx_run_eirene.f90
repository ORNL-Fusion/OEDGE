subroutine styx_run_eirene(ntps)
  use all_variables, only : global_parameters, interp_data2, global_variables, reference_parameters
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use eirmod_comprt, only : iunout
  use eirmod_comusr
  use eirmod_ccona
  use eirmod_cpes
  implicit none
  integer, intent(in) :: ntps
  integer :: ier
  real(dp) :: DUMMY,EIRENE_RESET_SECOND
  integer :: n_ion,iatm,itri,itor,icell
  integer :: atom_cnt

  ! tests
  real(dp) :: volpp(Ntor_cells),source_tot(Ntor_cells)

  include 'mpif.h'


  call mpi_barrier(MPI_COMM_WORLD,ier)

  if (ntps > 1) DUMMY=EIRENE_RESET_SECOND()

  if (my_pe==0) then
     write(*,*) '--------------------------------------------------------'
     write(*,*) '  Performing Monte Carlo calculation ... '
     write(*,*) '--------------------------------------------------------'
  endif

  call run_monte_carlo()

  if (my_pe==0) then
     if (direct_coupling) then
        call styx_process_sources_from_eirene(0)
        if (interface_test_mode) call styx_test_preav_sources(0)
     else
        call styx_process_sources_from_eirene(1)
        if (interface_test_mode) call styx_test_preav_sources(1)
     endif

!!!!!!!!!!!!!!!!!!   units at this point   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! atom density in m-3
     ! Sn in m-3.s-1
     ! Sm in (m.s-1.m-3).s-1
     ! SE in Watt.m-3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Interp_data2%tri_Srad(:,0,:)=0._dp

     do itor=1,Ntor_cells
       do itri=1,Ntri_styx
         icell = itri + (itor-1)*(Ntri_styx+1)
         Interp_data2%tri_SE(itri,0,itor)     = SE_tot(icell,0)
         atom_cnt=0
         do n_ion=1,global_parameters%N_species
           Interp_data2%tri_Sn(itri,n_ion,itor)     = Sn_tot(icell,n_ion)
           Interp_data2%tri_SG(itri,n_ion,itor)     = Sm_tot(icell,n_ion)
           Interp_data2%tri_SE(itri,n_ion,itor)     = SE_tot(icell,n_ion)
           Interp_data2%tri_Nn(itri,n_ion,itor)     = atom_density(n_ion,icell)
           Interp_data2%tri_Nm(itri,n_ion,itor)     = mol_density(1,icell)
           Interp_data2%tri_Tn(itri,n_ion,itor)     = T0_at(n_ion,icell)
           Interp_data2%tri_Tm(itri,n_ion,itor)     = T0_mol(1,icell)
           Interp_data2%tri_Srad(itri,n_ion,itor) = Srad_at(icell,n_ion)
         end do
         !!! total radiated power by neutral particle species here !!!!
         do iatm=1,natmi
           Interp_data2%tri_Srad(itri,0,itor)=Interp_data2%tri_Srad(itri,0,itor)+Srad_at(icell,iatm)
         enddo
       enddo
     enddo
     do n_ion=1,global_parameters%N_species
        Interp_data2%neutral_outflux(n_ion) = eos(0)%potat(n_ion,nlim+1) + 2*eos(0)%potml(1,nlim+1)
        ! Patrick: It is assumed that mollecules contain two ions for each species... will have to be improved for multi-species
     enddo
     

     ! check the integrals over PPs (to be removed later)
     
!     Source_tot=0._dp
!     volpp=0._dp
!     do itor=1,Ntor_cells
!       do itri=1,Ntri_styx
!         icell = itri + (itor-1)*(Ntri_styx+1)
!         source_tot(itor)=source_tot(itor)+Interp_data2%tri_Sn(itri,1,itor)*vol_tri_eirene(itri)
!         volpp(itor)=volpp(itor)+vol_tri_eirene(itri)
!       enddo
!     enddo
!     write(*,*) 'volumes of PPs = ',volpp(:)
!     write(*,*) 'total = ',sum(volpp)
!     write(*,*) 'sources per PPs = ',source_tot(:)
!     write(*,*) 'total = ',sum(source_tot)
!
!     ! atom flux on wall
!     write(*,*) 'Atom flux on core interface = ',eos(0)%potat(1,nlim+1),' part/s'
!     write(*,*) ' Sn_intg(0) + flux on core interface = ',sum(source_tot)+eos(0)%potat(1,nlim+1)
!      write(*,*) ' flux in eirene (part/s) = ',fluxt_eirene/elcha
!     !write(*,*) eos(0)%potat(1,nlim+2),eos(0)%potat(1,nlim+3)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         
     ! back to styx units for densities

     Interp_data2%tri_Nn=Interp_data2%tri_Nn/(reference_parameters%fields%n0)
     Interp_data2%tri_Nm=Interp_data2%tri_Nm/(reference_parameters%fields%n0)
     !    Interp_data2%tri_Nti=Interp_data2%tri_Nti/(reference_parameters%fields%n0)

     ! build total sources (directly from EIRENE output, prior to interpolation) for balance checks
     ! mass removed in momentum sources

     Sn_intg = Sn_at_intg + Sn_mol_intg + Sn_tion_intg + Sn_pls_intg
     Sm_intg = Sm_at_intg + Sm_mol_intg + Sm_tion_intg + Sm_pls_intg

     SE_intg = SE_at_intg + SE_mol_intg + SE_tion_intg + SE_pls_intg


     if (Sn_intg(1) > fluxt_eirene/elcha) then
        write(*,*) '----------------------------------------------------'
        write(*,*) ' There are more electrons created than available ...'
        write(*,*) ' fluxt_eirene (part/s) = ',fluxt_eirene/elcha
        write(*,*) ' Sn_intg      (part/s) = ',Sn_intg
        write(*,*) '----------------------------------------------------'
     endif

#ifdef S2D
    call styx_process_fluxes_new()
#endif

  endif




end subroutine styx_run_eirene
