subroutine styx_allocate_eirene_interface_structures
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use eirmod_cpes
  use eirmod_comusr
  implicit none
  integer*4 :: i,ier,nstratat

  include 'mpif.h'

! ion density, momentum, temperature

  if (my_pe ==0) then
    allocate(eiv_i(nplsi))

    do i=1,nplsi
      allocate(eiv_i(i)%dens(Neir_cells))
      allocate(eiv_i(i)%parmom(Neir_cells))
      allocate(eiv_i(i)%vx(Neir_cells))
      allocate(eiv_i(i)%vy(Neir_cells))
      allocate(eiv_i(i)%vz(Neir_cells))
      allocate(eiv_i(i)%T(Neir_cells))

      eiv_i(i)%dens=0._dp
      eiv_i(i)%parmom=0._dp
      eiv_i(i)%vx=0._dp
      eiv_i(i)%vy=0._dp
      eiv_i(i)%vz=0._dp
      eiv_i(i)%T=0._dp
    enddo

! allocate arrays for neutrals

    allocate(Gamma0X_at(natmi,Neir_cells))
    allocate(Gamma0Y_at(natmi,Neir_cells))
    allocate(Gamma0Z_at(natmi,Neir_cells))
    allocate(V0X_at(natmi,Neir_cells))
    allocate(V0Y_at(natmi,Neir_cells))
    allocate(V0Z_at(natmi,Neir_cells))
    allocate(Gammapar0_at(natmi,Neir_cells))
    allocate(V0par_at(natmi,Neir_cells))
    allocate(Gamma0X_mol(nmoli,Neir_cells))
    allocate(Gamma0Y_mol(nmoli,Neir_cells))
    allocate(Gamma0Z_mol(nmoli,Neir_cells))
    allocate(V0X_mol(nmoli,Neir_cells))
    allocate(V0Y_mol(nmoli,Neir_cells))
    allocate(V0Z_mol(nmoli,Neir_cells))
    allocate(Gammapar0_mol(nmoli,Neir_cells))
    allocate(V0par_mol(nmoli,Neir_cells))
    allocate(Gamma0X_tion(nioni,Neir_cells))
    allocate(Gamma0Y_tion(nioni,Neir_cells))
    allocate(Gamma0Z_tion(nioni,Neir_cells))
    allocate(V0X_tion(nioni,Neir_cells))
    allocate(V0Y_tion(nioni,Neir_cells))
    allocate(V0Z_tion(nioni,Neir_cells))
    allocate(Gammapar0_tion(nioni,Neir_cells))
    allocate(V0par_tion(nioni,Neir_cells))
  
! some of these my not be used in direct coupling, check
    allocate(atom_density(natmi,Neir_cells))
    allocate(mol_density(nmoli,Neir_cells))
    allocate(tion_density(nioni,Neir_cells))
    allocate(E0_at(natmi,Neir_cells))
    allocate(E0_mol(nmoli,Neir_cells))
    allocate(E0_tion(nioni,Neir_cells))
    allocate(T0_at(natmi,Neir_cells))
    allocate(T0_mol(nmoli,Neir_cells))
    allocate(T0_tion(nioni,Neir_cells))
    allocate(vpar_styx(Neir_cells))
    allocate(Sn_at(Neir_cells,0:nplsi),Sm_at(Neir_cells,1:nplsi))
    allocate(Sn_mol(Neir_cells,0:nplsi),Sm_mol(Neir_cells,1:nplsi))
    allocate(Sn_tion(Neir_cells,0:nplsi),Sm_tion(Neir_cells,1:nplsi))
    allocate(Sn_pls(Neir_cells,0:nplsi),Sm_pls(Neir_cells,1:nplsi))
    allocate(SE_at(Neir_cells,0:nplsi))
    allocate(SE_mol(Neir_cells,0:nplsi))
    allocate(SE_tion(Neir_cells,0:nplsi))
    allocate(SE_pls(Neir_cells,0:nplsi))
    allocate(Srad_at(Neir_cells,1:natmi))
  
    allocate(Sn_intg(0:nplsi),Sm_intg(1:nplsi),SE_intg(0:nplsi))
    allocate(Sn_at_intg(0:nplsi),Sm_at_intg(1:nplsi),SE_at_intg(0:nplsi))
    allocate(Sn_mol_intg(0:nplsi),Sm_mol_intg(1:nplsi),SE_mol_intg(0:nplsi))
    allocate(Sn_tion_intg(0:nplsi),Sm_tion_intg(1:nplsi),SE_tion_intg(0:nplsi))
    allocate(Sn_pls_intg(0:nplsi),Sm_pls_intg(1:nplsi),SE_pls_intg(0:nplsi))
    allocate(Srad_at_intg(1:natmi))
  
  !allocate(SE_i_io(Neir_cells),SE_i_cx(Neir_cells))

    allocate(Sn_tot(Neir_cells,0:nplsi),Sm_tot(Neir_cells,1:nplsi),SE_tot(Neir_cells,0:nplsi))
  endif ! my_pe == 0

  call mpi_bcast(sc_level,1,mpi_integer,0,mpi_comm_world,ier)
  call mpi_bcast(n_short_cycles,1,mpi_integer,0,mpi_comm_world,ier)
  call mpi_bcast(ns_refresh,1,mpi_integer,0,mpi_comm_world,ier)  
  call mpi_bcast(direct_coupling,1,mpi_logical,0,mpi_comm_world,ier)
  call mpi_bcast(Neir_cells,1,mpi_integer,0,mpi_comm_world,ier)
  call mpi_bcast(NSTRATA,1,mpi_integer,0,mpi_comm_world,ier)
  call mpi_bcast(timedep,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  if (timedep) then
    nstratat=nstrata+1
  else
    nstratat=nstrata
  endif

  allocate(eov(0:NSTRATAt))
  do i=0,NSTRATAt
    allocate(eov(i)%pdena(natmi,Neir_cells))
    allocate(eov(i)%pdenm(nmoli,Neir_cells))
    allocate(eov(i)%pdeni(nioni,Neir_cells))
    allocate(eov(i)%edena(natmi,Neir_cells))
    allocate(eov(i)%edenm(nmoli,Neir_cells))
    allocate(eov(i)%edeni(nioni,Neir_cells))
    allocate(eov(i)%vxdena(natmi,Neir_cells))
    allocate(eov(i)%vxdenm(nmoli,Neir_cells))
    allocate(eov(i)%vxdeni(nioni,Neir_cells))
    allocate(eov(i)%vydena(natmi,Neir_cells))
    allocate(eov(i)%vydenm(nmoli,Neir_cells))
    allocate(eov(i)%vydeni(nioni,Neir_cells))
    allocate(eov(i)%vzdena(natmi,Neir_cells))
    allocate(eov(i)%vzdenm(nmoli,Neir_cells))
    allocate(eov(i)%vzdeni(nioni,Neir_cells))

    eov(i)%pdena=0._dp
    eov(i)%pdenm=0._dp
    eov(i)%pdeni=0._dp
    eov(i)%edena=0._dp
    eov(i)%edenm=0._dp
    eov(i)%edeni=0._dp
    eov(i)%vxdena=0._dp
    eov(i)%vydena=0._dp
    eov(i)%vzdena=0._dp
    eov(i)%vxdenm=0._dp
    eov(i)%vydenm=0._dp
    eov(i)%vzdenm=0._dp
    eov(i)%vxdeni=0._dp
    eov(i)%vydeni=0._dp
    eov(i)%vzdeni=0._dp      
  enddo

! these are the direct coupling data

  do i=0,NSTRATAt
    allocate(eov(i)%pael(Neir_cells))
    allocate(eov(i)%pmel(Neir_cells))
    allocate(eov(i)%piel(Neir_cells))
    allocate(eov(i)%papl(nplsi,Neir_cells))
    allocate(eov(i)%pmpl(nplsi,Neir_cells))
    allocate(eov(i)%pipl(nplsi,Neir_cells))
    allocate(eov(i)%eael(Neir_cells))
    allocate(eov(i)%emel(Neir_cells))
    allocate(eov(i)%eiel(Neir_cells))
    allocate(eov(i)%eapl(Neir_cells))
    allocate(eov(i)%empl(Neir_cells))
    allocate(eov(i)%eipl(Neir_cells))
    allocate(eov(i)%eaplr(nplsi,Neir_cells))
    allocate(eov(i)%emplr(nplsi,Neir_cells))
    allocate(eov(i)%eiplr(nplsi,Neir_cells))
    allocate(eov(i)%mapl(nplsi,Neir_cells))
    allocate(eov(i)%mmpl(nplsi,Neir_cells))
    allocate(eov(i)%mipl(nplsi,Neir_cells))

    eov(i)%pael=0._dp
    eov(i)%pmel=0._dp
    eov(i)%piel=0._dp
    eov(i)%papl=0._dp
    eov(i)%pmpl=0._dp
    eov(i)%pipl=0._dp
    eov(i)%eael=0._dp
    eov(i)%emel=0._dp
    eov(i)%eiel=0._dp
    eov(i)%eapl=0._dp
    eov(i)%empl=0._dp
    eov(i)%eipl=0._dp
    eov(i)%eaplr=0._dp
    eov(i)%emplr=0._dp
    eov(i)%eiplr=0._dp
    eov(i)%mapl=0._dp
    eov(i)%mmpl=0._dp
    eov(i)%mipl=0._dp

  enddo

  ! check whether it is necessary to allocate this for all strata
  allocate(eos(0:NSTRATAt))
  do i=0,NSTRATAt
    allocate(eos(i)%potat(natmi,nsurf_tal))
    allocate(eos(i)%prfaat(natmi,nsurf_tal))
    allocate(eos(i)%prfmat(natmi,nsurf_tal))
    allocate(eos(i)%prfpat(natmi,nsurf_tal))
    allocate(eos(i)%potml(nmoli,nsurf_tal))
    allocate(eos(i)%prfaml(nmoli,nsurf_tal))
    allocate(eos(i)%prfmml(nmoli,nsurf_tal))
    allocate(eos(i)%prfpml(nmoli,nsurf_tal))
    allocate(eos(i)%potio(nioni,nsurf_tal))
    allocate(eos(i)%potpl(nplsi,nsurf_tal))
    allocate(eos(i)%eotat(natmi,nsurf_tal))
    allocate(eos(i)%erfaat(natmi,nsurf_tal))
    allocate(eos(i)%erfmat(natmi,nsurf_tal))
    allocate(eos(i)%erfiat(natmi,nsurf_tal))
    allocate(eos(i)%erfpat(natmi,nsurf_tal))
    allocate(eos(i)%eotml(nmoli,nsurf_tal))
    allocate(eos(i)%erfaml(nmoli,nsurf_tal))
    allocate(eos(i)%erfmml(nmoli,nsurf_tal))
    allocate(eos(i)%erfiml(nmoli,nsurf_tal))
    allocate(eos(i)%erfpml(nmoli,nsurf_tal))
    allocate(eos(i)%eotio(nioni,nsurf_tal))
    allocate(eos(i)%erfaio(nioni,nsurf_tal))
    allocate(eos(i)%erfmio(nioni,nsurf_tal))
    allocate(eos(i)%erfiio(nioni,nsurf_tal))
    allocate(eos(i)%erfpio(nioni,nsurf_tal))
    allocate(eos(i)%eotpl(nplsi,nsurf_tal))
    ! sputtering data
    allocate(eos(i)%sptat(natmi,nsurf_tal))
    allocate(eos(i)%sptml(nmoli,nsurf_tal))
    allocate(eos(i)%sptio(nioni,nsurf_tal))
    allocate(eos(i)%sptpl(nplsi,nsurf_tal))
    allocate(eos(i)%spttot(nsurf_tal))

    allocate(eos(i)%spump(nspz,nsurf_tal))
  enddo


end subroutine styx_allocate_eirene_interface_structures

