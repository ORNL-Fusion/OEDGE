subroutine styx_load_param_eirene()
  use styx2eirene
  use Mreaders
  implicit none
  integer :: finput,ipfc

  finput=10
  open(unit=finput,file='eirene_coupling.txt',status='unknown')

  call read_eirene_am_data_block(finput)
  call read_wall_parameters_block(finput)
  call read_pump_parameters_block(finput)
  call read_gas_puffs_block(finput)
  call read_coupling_mode_block(finput)
  call read_data_histories_block(finput)
  call read_miscellaneous_block(finput)
 
  do ipfc=1,n_pfc_types
    pfc_models(ipfc)%material=adjustl(pfc_models(ipfc)%material)
  
    if ((pfc_models(ipfc)%material /= 'Be').and.(pfc_models(ipfc)%material /= 'C ').and.(pfc_models(ipfc)%material /= 'Fe') &
       .and.(pfc_models(ipfc)%material /= 'Mo').and.(pfc_models(ipfc)%material /='W ')) then
      write(*,*) ' Wall material selection not correct for pfc type #,',ipfc,' : material = ',pfc_models(ipfc)%material
      write(*,*) ' Exit .... '
      call eirene_exit_own(1)  
    endif

    if (pfc_models(ipfc)%sputer_model /= 0 .and. hardwired) then
      write(*,*) '**************************************************************'
      write(*,*) '* hardwired option currently not supported with sputering on *'
      write(*,*) '* hardwired turned off ....                                  *'
      write(*,*) '**************************************************************'
      hardwired=.false.
    endif
  enddo
  
  
  if (sheath_model /=0 .and. sheath_model /= 1) then
     write(*,*) '**************************************************************'
     write(*,*) '* Invalid sheath model, switch to default (0)           .... *'
     write(*,*) '**************************************************************'
     sheath_model=0
  endif
  
  
  close(finput)
  
  ! orphan parameters/statements
  source_species=4
  seed_gap=0
  if (sc_level > 1 .and. ns_refresh > n_short_cycles) ns_refresh=n_short_cycles

  end subroutine styx_load_param_eirene

  subroutine read_eirene_am_data_block(finput)
  use eirmod_precision 
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: finput
  call skip_line(finput,10)
  read(finput,1) am_database
  call skip_line(finput,2)
  read(finput,2) tweak_chemistry
  call skip_line(finput,2)
  read(finput,3) hardwired
  call skip_line(finput,2) 
1 format(13x,i2)
2 format(10x,L1)
3 format(14x,L1)

  end subroutine read_eirene_am_data_block

  subroutine read_wall_parameters_block(finput)
  use all_variables, only : global_parameters
  use eirmod_precision 
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: finput
  integer*4 :: N_species,n,ipfc
  character(128) :: string
  real*8, allocatable :: bufferF(:)
  N_species=global_parameters%n_species
  allocate(bufferF(1:N_species))
  
  call skip_line(finput,5)

  read(finput,1) n_pfc_types
  allocate(pfc_models(n_pfc_types))
  do ipfc=1,n_pfc_types
    allocate(pfc_models(ipfc)%R(N_species))
  enddo
  do ipfc=1,n_pfc_types
    call skip_line(finput,3)
    read(finput,2) pfc_models(ipfc)%material
    call skip_line(finput,2)
    read(finput,3) String
    call parse_line_float(String,N_species,bufferF)
    do n=1,N_species
       pfc_models(ipfc)%R(n)=bufferF(n)
    end do    
    call skip_line(finput,2)
    read(finput,4) pfc_models(ipfc)%T
    call skip_line(finput,2) 
    read(finput,1) pfc_models(ipfc)%sputer_model
    call skip_line(finput,2)
    read(finput,5) pfc_models(ipfc)%sputer_yield_phys
    call skip_line(finput,2)
    read(finput,5) pfc_models(ipfc)%sputer_yield_chem
  enddo

  call skip_line(finput,3)
  read(finput,1) sheath_model
  call skip_line(finput,2) 
  deallocate(bufferF)

1 format(14x,i2) 
2 format(11x,A2)
3 format(A128)
4 format(3x,ES9.2E2)
5 format(8x,ES9.2E2)

  end subroutine read_wall_parameters_block

  subroutine read_pump_parameters_block(finput)
    use all_variables, only : global_parameters
    use eirmod_precision 
    use styx2eirene
    use Mreaders
    implicit none
    integer, intent(in) :: finput
    real*8, allocatable :: bufferF(:)
    logical, allocatable :: bufferL(:)
    character(128) :: string
    integer*4 :: N_species,n,ipump

    N_species=global_parameters%n_species
    allocate(bufferF(1:N_species))
    allocate(bufferL(1:N_species))

    call skip_line(finput,5)

    read(finput,1) n_pumps
    allocate(pumps(n_pumps))

    do ipump=1,n_pumps
       allocate(pumps(ipump)%R(N_species))
       call skip_line(finput,2)
       read(finput,2) String
       call parse_line_float(String,N_species,bufferF)
       do n=1,N_species
          pumps(ipump)%R(n)=bufferF(n)
       end do
       call skip_line(finput,2)
       read(finput,3) pumps(ipump)%isSpeedSet
       call skip_line(finput,2)
       read(finput,4) pumps(ipump)%pumping_speed
    enddo
    call skip_line(finput,2)

1   format(9x,i2)
2   format(A128)
3   format(12x,L1)
4   format(16x,F6.2)
  end subroutine read_pump_parameters_block

  subroutine read_gas_puffs_block(finput)
  use all_variables, only : global_parameters
  use eirmod_precision 
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: finput
  integer*4 :: n,nstrat0
  character(128) :: string
  real*8, allocatable :: bufferF(:)
  integer*4, allocatable :: bufferI(:)
    
  call skip_line(finput,5)
  read(finput,1) Npuffs
  if (Npuffs > 0) then
    allocate(bufferF(1:NPuffs),bufferI(1:Npuffs))
    Nstrat0=2*global_parameters%n_species
    allocate(puffs(Nstrat0+1:Nstrat0+Npuffs))
    call skip_line(finput,2)
    read(finput,2) String
    call parse_line_elements_puff(String,Npuffs,nstrat0)
    call skip_line(finput,2)
    read(finput,2) string
    call parse_line_float(String,Npuffs,bufferF)
    do n=1,Npuffs
      puffs(nstrat0+n)%rate=bufferF(n)
      write(*,*) 'puff', puffs(nstrat0+n)
    enddo 
    call skip_line(finput,2)
    read(finput,2) string
    call parse_line_float(String,Npuffs,bufferF)
    do n=1,Npuffs
      puffs(nstrat0+n)%T0=bufferF(n)
    enddo
    call skip_line(finput,2)
    read(finput,2) string
    call parse_line_float(String,Npuffs,bufferF)
    do n=1,Npuffs
      puffs(nstrat0+n)%divergence=bufferF(n)
    enddo
    call skip_line(finput,2)
    read(finput,2) string
    call parse_line_integer(String,Npuffs,bufferI)
    do n=1,Npuffs
      puffs(nstrat0+n)%itri=bufferI(n)
    enddo
    read(finput,2) string
    call parse_line_integer(String,Npuffs,bufferI)
    do n=1,Npuffs
      puffs(nstrat0+n)%iside=bufferI(n)
    enddo
    read(finput,2) string
    call parse_line_integer(String,Npuffs,bufferI)
    do n=1,Npuffs
      puffs(nstrat0+n)%itor=bufferI(n)
    enddo
    call skip_line(finput,2)
    deallocate(bufferF,bufferI)
  else
    call skip_line(finput,19)
  endif

1 format(8x,i2)
2 format(A128)  
  end subroutine read_gas_puffs_block

  subroutine read_coupling_mode_block(finput)
  use all_variables, only : global_parameters
  use eirmod_precision 
  use styx2eirene
  use Mreaders
  use Meirene_vars
  implicit none
  integer, intent(in) :: finput
  allocate(Npart_eirene(global_parameters%n_species*2+Npuffs))
  allocate(seed_eirene(global_parameters%n_species*2+Npuffs))
  call skip_line(finput,5)
  read(finput,1) direct_coupling
  call skip_line(finput,2)
  read(finput,2) sc_level
  call skip_line(finput,2)
  read(finput,3) n_short_cycles
  call skip_line(finput,2)
  read(finput,4) ns_refresh
  call skip_line(finput,2)
  read(finput,5) timedep
  call skip_line(finput,2)
  read(finput,7) nprnli_styx
  call skip_line(finput,2)
  read(finput,8) nptst_styx
  call skip_line(finput,2)
  read(finput,6) eirene_vars%feedback
  call skip_line(finput,2)
1 format(18x,L1)
2 format(10x,i2)
3 format(16x,i6)
4 format(12x,i6)
5 format(11x,L1)
6 format(10x,i2)
7 format(14x,i8)
8 format(15x,i8)
  end subroutine read_coupling_mode_block
 

  subroutine read_data_histories_block(finput)
  use eirmod_precision 
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: finput
  call skip_line(finput,5)
  read(finput,1) npart_eirene(1)
  call skip_line(finput,2)
  read(finput,2) seed_eirene(1)
  call skip_line(finput,2) 
1 format(14x,i8)
2 format(13x,i6)
  end subroutine read_data_histories_block

  subroutine read_miscellaneous_block(finput)
  use eirmod_precision 
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: finput
  call skip_line(finput,5)
  read(finput,1) nhist_plot
  call skip_line(finput,2)
  read(finput,1) stratum_plot
  call skip_line(finput,2)
  read(finput,2) rad_mat
  call skip_line(finput,2)
  read(finput,3) interface_test_mode
  call skip_line(finput,2)
  read(finput,4) Temin_eirene
1 format(15x,i3)
2 format(10x,L1)
3 format(12x,L1)
4 format(15x,F5.2)  
  end subroutine read_miscellaneous_block
  


