subroutine styx_initialize_pfc_models
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comprt, only : iunout
  use styx2eirene
  implicit none
  integer*4 :: ncount,i,isp,j,ipfc
  character(2) :: symbol


  ! turn sputtering on/off

  do ipfc=1,n_pfc_types
    if (pfc_models(ipfc)%sputer_model == 1 .or. pfc_models(ipfc)%sputer_model == 2) then
      pfc_models(ipfc)%sputer_on=.true.
    else
      pfc_models(ipfc)%sputer_on=.false.
    endif
    ! sputering yields are read in % (sputer_model=1, otherwise irrelevant)
    if (pfc_models(ipfc)%sputer_model == 1) then
      pfc_models(ipfc)%sputer_yield_phys=pfc_models(ipfc)%sputer_yield_phys/100._dp
      pfc_models(ipfc)%sputer_yield_chem=pfc_models(ipfc)%sputer_yield_chem/100._dp
    endif
  enddo

  ! check if neutral species are already declared in soledge

  ncount=1
  do i=1,n_pfc_types
    pfc_models(i)%iatm=0
    do isp=2,global_parameters%n_species
      symbol=global_parameters%element_list(isp)%symbol
      if (adjustr(symbol) == adjustr(pfc_models(i)%material)) pfc_models(i)%iatm = isp
    enddo
    if (pfc_models(i)%iatm == 0 .and. pfc_models(i)%sputer_on) then
      write(iunout,*) 'Wall material atom species added only as neutrals in EIRENE, species  = ',pfc_models(i)%material
      pfc_models(i)%iatm=global_parameters%n_species+ncount
      ncount=ncount+1
    endif
  enddo
  
! count and identify material types

  allocate(materials(n_pfc_types))
  do i=1,n_pfc_types
    materials(i)%symbol='  '
    materials(i)%iatm=0
    materials(i)%ipfc=0
  enddo

  n_mat=0
  do i=1,n_pfc_types       
    ! now look in array if this material has already been identified
    ncount=0    
    do j=1,n_mat
      if (materials(j)%symbol == pfc_models(i)%material) then
        ncount=ncount+1
      endif
    enddo
    if (ncount == 0) then
      ! the material has not been identified yet,add it to the list
      n_mat=n_mat+1  
      materials(n_mat)%symbol=pfc_models(i)%material
      materials(n_mat)%iatm=pfc_models(i)%iatm
      materials(n_mat)%ipfc=i             
    endif     
  enddo

  do i=1,n_pumps
    pumps(i)%material='Fe'
  enddo


end subroutine styx_initialize_pfc_models
