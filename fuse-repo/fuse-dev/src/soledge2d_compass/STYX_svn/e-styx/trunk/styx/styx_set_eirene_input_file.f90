subroutine set_eirene_input_file()
  use all_variables, only : global_parameters
  use eirmod_precision
  use styx2eirene
  use eirmod_comprt, only : iunout
  use Mreaders
  implicit none
  integer :: unit_source, unit_destin
  integer :: nreac,ireac

  if (global_parameters%n_species > 100) then
    write(*,*) 'warning, number of species > 100, eirene input &
                file may be incorretly written ...' 
  endif
  
  unit_source = 2
  unit_destin = 1

  ! checks on AM database for Hydrogen isotopes

  select case (am_database)

    case(1)
      open(unit=unit_source,file='eirene_input_file_AMJUEL_level_1')
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 'EIRENE uses AMJUEL database, level 1                          '
      write(*,*) '--------------------------------------------------------------'
    case(2)
      open(unit=unit_source,file='eirene_input_file_ADAS_atoms')
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 'EIRENE uses ADAS database for atoms                           '
      write(*,*) '--------------------------------------------------------------'
    case(3)
      open(unit=unit_source,file='eirene_input.template')
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 'EIRENE uses a non standard input file, use at your own risk !'
      write(*,*) '--------------------------------------------------------------'
    case(4)
      open(unit=unit_source,file='eirene_input_file_AMJUEL_level_2')
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 'EIRENE uses AMJUEL database, level 2                          '
      write(*,*) '--------------------------------------------------------------'
    case default
      open(unit=unit_source,file='eirene_input_file_AMJUEL_level_1')
      write(*,*) '----------------------------------------------------------------'
      write(*,*) 'Wrong specification of am_database, default AMJUEL level 1 used '
      write(*,*) '----------------------------------------------------------------'
  end select

  if (tweak_chemistry) then
    write(*,*) '----------------------------------------------------------------'
    write(*,*) ' Warning, chemistry tweaks are allowed (see report below)       '
    write(*,*) '----------------------------------------------------------------'
  endif

  
  call write_time_dependent_mode_data(unit_source,unit_destin,1)
  call write_data_for_spatially_resolved_tallies(unit_source,unit_destin)
  call write_string_to_file(unit_destin,'CFILE ADAS ADAS/adf11')
  call write_grid_parameters(unit_source,unit_destin) 
  call write_reflection_model_data(unit_source,unit_destin)
  read(unit_source,'(i6)') nreac
  call write_number_of_AM_processes(unit_source,unit_destin,nreac)
  ireac=nreac
  call write_reaction_cards_for_impurities(ireac,unit_source,unit_destin)
  call write_right_isotope(unit_source,unit_destin,1)  
  ireac=nreac
  call write_atomic_species_cards_for_impurities(ireac,unit_source,unit_destin)
  call write_right_isotope(unit_source,unit_destin,2)
  call write_right_isotope(unit_source,unit_destin,3)
  call write_right_isotope(unit_source,unit_destin,4)
  ireac=nreac  
  call write_bulk_ions_species_cards_for_impurities(ireac,unit_source,unit_destin)
  call write_trim_data_needed(unit_source,unit_destin) 
  call write_dummy_species_sampling_distributions(unit_source,unit_destin)
  call write_strata_data(unit_source,unit_destin,1)
  call write_data_for_ion_resolved_energy_sources(unit_source,unit_destin)  
  call write_strata_data(unit_source,unit_destin,2)
  call write_tallies_to_be_calculated(unit_source,unit_destin)  
  call write_time_dependent_mode_data(unit_source,unit_destin,2)
  close(unit_source)

end subroutine



subroutine copy_file(unit_read,unit_write,num) 
  implicit none
  integer, intent(in) :: unit_read,unit_write,num
  integer :: nc,i
  character(128) :: string
  character(3) :: ncc
  character(20) :: format
  
  do i=1,num
    read(unit_read,'(A128)') string  
    nc=len(trim(string))
    write(ncc,'(i3)') nc
    format='(A'//trim(adjustl(ncc))//')'
    write(unit_write,format) trim(string)
  enddo
end subroutine copy_file

subroutine write_string_to_file(unit_write,string)
  implicit none
  integer, intent(in) :: unit_write
  integer :: nc
  character(*),intent(in) :: string
  character(3) :: ncc
  character(20) :: format

  nc=len(trim(string))
  write(ncc,'(i3)') nc
  format='(A'//trim(adjustl(ncc))//')'
  write(unit_write,format) trim(string)

end subroutine write_string_to_file

subroutine copy_lines_until_finding(unit_read,string,unit_write)
  implicit none
  integer, intent(in) :: unit_read,unit_write
  character(*), intent(in) :: string  
  character(128) :: A
  character(3) :: ncc
  character(20) :: format  
  integer :: nc,lgth

  A='xxx'
  lgth=len(string)
  do while (A(1:lgth) .ne. string) 
    read(unit_read,'(A128)') A
    nc=len(trim(A))
    write(ncc,'(i3)') nc
    format='(A'//trim(adjustl(ncc))//')'
    write(unit_write,format) trim(A)
  enddo

end subroutine copy_lines_until_finding

subroutine write_right_isotope(unit_source,unit_destin,itype)
  use all_variables, only : global_parameters
  use Mreaders
  use styx2eirene, only : pfc_models, n_pfc_types, n_mat, materials
  implicit none
  integer, intent(in) :: itype,unit_source,unit_destin
  character(128) :: string
  integer :: nreacH,add_spc,ipfc,imat
  select case(itype)
    case(1)
      read(unit_source,'(A128)') string
      ! now set the right isotope and mass
      if (string(4:11) .ne. 'H      ' .and. string(12:14).ne.'  1') then
        write(*,*) 'Hydrogen is not species 1, exit ...'
        call eirene_exit_own(1) 
      endif
      if (global_parameters%element_list(1)%symbol == 'D' ) then
      string(4:11) ='D       '
      string(12:14)='  2'
    elseif (global_parameters%element_list(1)%symbol == 'T') then
      string(4:11) ='T       '
      string(12:14)='  3'
    elseif (global_parameters%element_list(1)%symbol == 'H') then
      string(4:11) ='H       '
      string(12:14)='  1'
    endif
    call write_string_to_file(unit_destin,string)   
    ! get number of reactions for hydrogen
    read(string(34:35),'(i2)') nreacH
    call copy_file(unit_source,unit_destin,2*nreacH)
   case(2)
     read(unit_source,'(A128)') string
     ! now set the right isotope and mass
     if (string(4:11) .ne. 'H2     ' .and. string(12:14).ne.'  2') then
       write(*,*) 'Di-hydrogen molecule is not species 1, exit ...'
       call eirene_exit_own(1) 
     endif
    
     if (global_parameters%element_list(1)%symbol == 'D' ) then
      string(4:11) ='D2      '
      string(12:14)='  4'
     elseif (global_parameters%element_list(1)%symbol == 'T') then
      string(4:11) ='T2      '
      string(12:14)='  6'
     elseif (global_parameters%element_list(1)%symbol == 'H') then
      string(4:11) ='H2      '
      string(12:14)='  2'
    endif
    call write_string_to_file(unit_destin,string)
    call copy_lines_until_finding(unit_source,'* TEST ION SPECIES CARDS',unit_destin)
    call copy_file(unit_source,unit_destin,1)
   case(3)
     read(unit_source,'(A128)') string
     ! now set the right isotope and mass
     if (string(4:11) .ne. 'H2+    ' .and. string(12:14).ne.'  2') then
       write(*,*) 'Di-hydrogen molecular ion is not species 1, exit ...'
       call eirene_exit_own(1) 
     endif
     if (global_parameters%element_list(1)%symbol == 'D' ) then
       string(4:11) ='D2+     '
       string(12:14)='  4'
     elseif (global_parameters%element_list(1)%symbol == 'T') then
       string(4:11) ='T2+     '
       string(12:14)='  6'
     elseif (global_parameters%element_list(1)%symbol == 'H') then
       string(4:11) ='H2+     '
       string(12:14)='  2'
     endif
     call write_string_to_file(unit_destin,string)
     call copy_lines_until_finding(unit_source,'* BULK ION SPECIES CARDS',unit_destin)  
     call skip_line(unit_source,1)
     ! add one singly charged ion per material for which species is not followed in styx
     add_spc=0
     do imat=1,n_mat
       ipfc = materials(imat)%ipfc
       if (pfc_models(ipfc)%sputer_on .and. pfc_models(ipfc)%iatm > global_parameters%n_species)&
                   add_spc=add_spc+1
     enddo
     write(unit_destin,'(i5)') global_parameters%n_ions+add_spc
    case(4)
      read(unit_source,'(A128)') string
      ! now set the right isotope and mass
      if (string(4:11) .ne. 'H+     ' .and. string(12:14).ne.'  1') then
        write(*,*) 'hydrogen ion is not species 1, exit ...'
        call eirene_exit_own(1) 
      endif
      if (global_parameters%element_list(1)%symbol == 'D' ) then
        string(4:11) ='D+      '
        string(12:14)='  2'
      elseif (global_parameters%element_list(1)%symbol == 'T') then
        string(4:11) ='T+      '
        string(12:14)='  3'
      elseif (global_parameters%element_list(1)%symbol == 'H') then
        string(4:11) ='H+      '
        string(12:14)='  1'
      endif
      call write_string_to_file(unit_destin,string)
      ! get number of reactions for hydrogen
      read(string(34:35),'(i2)') nreacH
      call copy_file(unit_source,unit_destin,2*nreacH)
  end select
end subroutine write_right_isotope

subroutine write_reaction_cards_for_impurities(ireac,unit_source,unit_destin)
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comprt, only :iunout
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(inout) :: ireac
  integer, intent(in) :: unit_destin,unit_source
  integer :: isp,icx,length,j,add_spc,ipfc,imat
  character(2) :: yearc
  character(2) :: symbol
  character(3) :: ireacc,massc
  character(12) :: chic
  character(128) :: A

!  real*8, allocatable :: single_poly_fit(:,:)

  do isp =2,global_parameters%n_species
    ireac=ireac+1

    ! later to be defined as input
    yearc="96"
    symbol=global_parameters%element_list(isp)%symbol
    call eirene_lowercase(symbol)

    ! ionization rate coefficient
    write(ireacc,'(i3)') ireac
    write(massc,'(i3)') nint(global_parameters%element_list(isp)%mass)
    A=ireacc//' ADAS   H.4 scd'//yearc//'    EI   0'//adjustr(massc)
    call write_string_to_file(unit_destin,A)
    A='    '//symbol//'   1'
    call write_string_to_file(unit_destin,A)

    ! enforce the use of <sigma v>(T) only for pre-averaged coupling mode
    if (direct_coupling) then
      do icx=1,global_parameters%element_list(isp)%ncx
        ireac=ireac+1
        write(ireacc,'(i3)') ireac
        A=ireacc//' '//global_parameters%element_list(isp)%cxdatas(icx)%cx_database1
        call write_string_to_file(unit_destin,A)
        ! same reaction index (cross section + rate coefficient)
        A=ireacc//' '//global_parameters%element_list(isp)%cxdatas(icx)%cx_database2
        call write_string_to_file(unit_destin,A)
      enddo
    else
      ! temporary, integrated in cx database afterwards
      !allocate(single_poly_fit(global_parameters%element_list(isp)%ncx,9))
      !single_poly_fit=0.d0
      ! just a constant for the moment, <sigma v> = 1e-8 cm3.s-1
      !single_poly_fit(:,1)=log(1d-8)

      do icx=1,global_parameters%element_list(isp)%ncx
        ireac=ireac+1
        write(ireacc,'(i3)') ireac
        A=ireacc//' '//global_parameters%element_list(isp)%cxdatas(icx)%cx_database1
        call write_string_to_file(unit_destin,A)
        ! same reaction index (cross section + rate coefficient)
        length=len(global_parameters%element_list(isp)%cxdatas(icx)%cx_database2)
        A=ireacc//' '//'CONST  H.2'//'  FT 000  '//'CX '//&
              global_parameters%element_list(isp)%cxdatas(icx)%cx_database2(length-6+1:length)
        call write_string_to_file(unit_destin,A)
        ! here comes data for the single polynomial fit (WARNING: only 4 digits here !!!)
        write(unit_destin,'(6e12.4)') global_parameters%element_list(isp)%cxdatas(icx)%fit_coeffs(1:6)
        write(unit_destin,'(3e12.4)') global_parameters%element_list(isp)%cxdatas(icx)%fit_coeffs(7:9)
      enddo
!      deallocate(single_poly_fit)
    endif
    ireac=ireac+1
          
    ! recombination rate coefficient
    write(ireacc,'(i3)') ireac 
    A=ireacc//' ADAS   H.4 acd'//yearc//'    RC   0'//adjustr(massc)
    call write_string_to_file(unit_destin,A)

    A='    '//symbol//'   1'
    call write_string_to_file(unit_destin,A)

    ! now data for radiated power

    ireac=ireac+1
    write(ireacc,'(i3)') ireac
    write(chic,'(e12.4)') -1.d0*global_parameters%element_list(isp)%amdatas(1)%ionization_potential
    A=ireacc//' ADAS   H.10plt'//yearc//'    EI   0'//adjustr(massc)//trim(adjustr(chic))//' 0.00000E+00 0.00000E+00'
    call write_string_to_file(unit_destin,A)
    A='    '//symbol//'   1'
    call write_string_to_file(unit_destin,A)

    ireac=ireac+1
    write(ireacc,'(i3)') ireac
    write(chic,'(e12.4)') global_parameters%element_list(isp)%amdatas(1)%ionization_potential 
    A=ireacc//' ADAS   H.10prb'//yearc//'    RC   0'//adjustr(massc)//trim(adjustr(chic))//' 0.00000E+00 0.00000E+00'
    call write_string_to_file(unit_destin,A)
    A='    '//symbol//'   1'
    call write_string_to_file(unit_destin,A)
    

  enddo    

  add_spc=0
  do imat=1,n_mat
    ipfc=materials(imat)%ipfc
    if (adjustr(pfc_models(ipfc)%material) /= adjustr(materials(imat)%symbol)) then
      write(iunout,*) 'inconsistency between pfc_models and materials symbol (set_eirene_input_file)'
      write(iunout,*) 'imat = ',imat,' ipfc = ',ipfc
      write(iunout,*) 'pfc_models(',ipfc,')%material = ',pfc_models(ipfc)%material
      write(iunout,*) 'materials(',imat,')%symbol = ',materials(imat)%symbol
      call eirene_exit_own(1)
    endif  

    if (pfc_models(ipfc)%sputer_on .and. pfc_models(ipfc)%iatm > global_parameters%n_species) then
      add_spc=add_spc+1
      ireac=ireac+1
      write(ireacc,'(i3)') ireac
      select case (pfc_models(ipfc)%material)
        case("Be")
        A= ireacc//' ADAS   H.4 scd96    EI   0  9'
        call write_string_to_file(unit_destin,A)
        A='    be   1'
        call write_string_to_file(unit_destin,A)

        case("C")
        A=ireacc//' ADAS   H.4 scd96    EI   0 12'
        call write_string_to_file(unit_destin,A)
        A='    c   1'
        call write_string_to_file(unit_destin,A)

        case("W")
        A=ireacc//' ADAS   H.4 scd50    EI   0184'
        call write_string_to_file(unit_destin,A)
        A='    w   1'
        call write_string_to_file(unit_destin,A)
      end select
    endif
  enddo

  call copy_lines_until_finding(unit_source,'* NEUTRAL ATOMS SPECIES CARDS',unit_destin)
  call skip_line(unit_source,1)

  write(unit_destin,'(i5)') global_parameters%n_species+add_spc
 
end subroutine write_reaction_cards_for_impurities

subroutine write_atomic_species_cards_for_impurities(ireac,unit_source,unit_destin)
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  implicit none
  integer, intent(inout) :: ireac
  integer, intent(in) :: unit_destin,unit_source
  integer :: isp,ising,icx,nreacisp,ipfc,imat,nions
  character(128) :: A
  character(2) :: ispc,nreacispc,isingc,symbol,iawallc,iionc
  character(3) :: massc,Zc,ireacc
  character(6) :: ibulk,iscd1,iscd2,iscde
  character(12) :: ireacec

  do isp=2,global_parameters%n_species
    ireac=ireac+1
    write(ispc,'(i2)') isp
    symbol=global_parameters%element_list(isp)%symbol
    write(massc,'(i3)') nint(global_parameters%element_list(isp)%mass)
    write(Zc,'(i3)') global_parameters%element_list(isp)%Z
    nreacisp=global_parameters%element_list(isp)%ncx+1
    write(nreacispc,'(i2)') nreacisp
    ! assumes ncx+1 reactions per species (recombination: ion)
    A=adjustr(ispc)//' '//adjustl(symbol)//'      '//adjustr(massc)//adjustr(Zc)//'  1  0'//' '//adjustr(ispc)//' '//adjustr(ispc)//'  0 '//adjustr(nreacispc)//'  0  0'
    call write_string_to_file(unit_destin,A)
    write(ireacc,'(i3)') ireac
    ! here ispc is the index for the corresponding ions
    ising=singly_charged_ions(isp)
    write(isingc,'(i2)') ising
    A='   '//adjustr(ireacc)//'   115  '//adjustr(isingc)//'14     0 30000'
    call write_string_to_file(unit_destin,A)
    write(ireacec,'(e12.5)') float(ireac)+2._dp
    A=ireacec//' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
    call write_string_to_file(unit_destin,A)

    do icx=1,global_parameters%element_list(isp)%ncx
      ireac=ireac+1
      write(ireacc,'(i3)') ireac
      ibulk=global_parameters%element_list(isp)%cxdatas(icx)%ibulk
      iscd1=global_parameters%element_list(isp)%cxdatas(icx)%iscd1
      iscd2=global_parameters%element_list(isp)%cxdatas(icx)%iscd2
      iscde=global_parameters%element_list(isp)%cxdatas(icx)%iscde
      A='   '//ireacc//adjustr(ibulk)//adjustr(iscd1)//adjustr(iscd2)//adjustr(iscde)
      call write_string_to_file(unit_destin,A)
      A=' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
      call write_string_to_file(unit_destin,A)
    enddo
    ! skip recombination and radiated power
    ireac=ireac+3
  enddo

  ! add sputtered neutral species if not there already
  nions=global_parameters%n_ions+1
  do imat=1,n_mat
    ipfc = materials(imat)%ipfc 
    if (pfc_models(ipfc)%sputer_on .and. pfc_models(ipfc)%iatm > global_parameters%n_species) then
      write(iawallc,'(i2)') pfc_models(ipfc)%iatm     
      write(iionc,'(i2)') nions
      ireac=ireac+1
      write(ireacc,'(i3)') ireac

      select case (pfc_models(ipfc)%material)
      ! to be finished - case where sputtered species is not followed by styx
      ! ion index = nion, add 1 and should be the right order: see below)
      ! ionization potential should come from an external database, not be hardwired

      case("Be")      
        A=trim(adjustr(iawallc))//' Be        9  4  1  0 '//trim(adjustr(iawallc))//' '//trim(adjustr(iawallc))//'  0  1  0  0'
        call write_string_to_file(unit_destin,A)
        A='   '//ireacc//'   115  '//adjustr(iionc)//'14     0 00000'
        call write_string_to_file(unit_destin,A)
        A=' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
        call write_string_to_file(unit_destin,A)
      case("C")
        A=trim(adjustr(iawallc))//' C       12  6  1  0 '//trim(adjustr(iawallc))//' '//trim(adjustr(iawallc))//'  0  1  0  0'
        call write_string_to_file(unit_destin,A)
        A='   '//ireacc//'   115  '//adjustr(iionc)//'14     0 00000'
        call write_string_to_file(unit_destin,A)
        A=' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
        call write_string_to_file(unit_destin,A)    
      case("W")
        A=trim(adjustr(iawallc))//' W       184 74  1  0 '//trim(adjustr(iawallc))//' '//trim(adjustr(iawallc))//'  0  1  0  0'
        call write_string_to_file(unit_destin,A)
        A='   '//ireacc//'   115  '//adjustr(iionc)//'14     0 00000'
        call write_string_to_file(unit_destin,A)
        A=' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
        call write_string_to_file(unit_destin,A)
      end select

      nions=nions+1

    endif
  enddo
  call copy_lines_until_finding(unit_source,'* NEUTRAL MOLECULES SPECIES CARDS',unit_destin) 
  call copy_file(unit_source,unit_destin,1)

end subroutine write_atomic_species_cards_for_impurities

subroutine write_bulk_ions_species_cards_for_impurities(ireac,unit_source,unit_destin)
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  implicit none
  integer, intent(inout) :: ireac
  integer, intent(in) :: unit_destin,unit_source
  integer :: iion,isp,iz,i,nions,ipfc,imat
  character(2) :: ispc,iionc,symbol,iiwallc
  character(3) :: massc,Zc,ireacc,izc,iawallc
  character(8) :: symbol_ion
  character(12) :: ireacec,chic
  character(128) :: A

  iion=2

  do isp=2,global_parameters%n_species
    ! skip ionization
    ireac=ireac+1+global_parameters%element_list(isp)%ncx+1
    write(ispc,'(i2)') isp
    write(iionc,'(i2)') iion
    symbol=global_parameters%element_list(isp)%symbol
    write(massc,'(i3)') nint(global_parameters%element_list(isp)%mass)
    write(Zc,'(i3)') global_parameters%element_list(isp)%Z
 
    ! first charge state (recombination into neutrals)
    
    symbol_ion=adjustl(adjustr(symbol)//'+')
  
    ! assumes 1 reactions per species (recombination: ion)
    A=adjustr(iionc)//' '//adjustl(symbol_ion)//adjustr(massc)//adjustr(Zc)//'  1  1'//' '//adjustr(ispc)//' '//adjustr(ispc)//'  0  1  0  0'
    call write_string_to_file(unit_destin,A)   

    write(ireacc,'(i3)') ireac
    ! here ispc is the index for the corresponding ions
    A='   '//adjustr(ireacc)//'   115  '//adjustr(ispc)//'11     0 30000'
    call write_string_to_file(unit_destin,A)

    write(ireacec,'(e12.5)') float(ireac)+2._dp
    !write(chic,'(e12.5)') global_parameters%element_list(isp)%amdatas(1)%ionization_potential
    !A=chic//' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
    A=ireacec//' 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E 00'
    call write_string_to_file(unit_destin,A)

    iion=iion+1

    ! other charge states, just background densities
    do iz=2,global_parameters%element_list(isp)%Z
      ! assumes 1 reactions per species (recombination: ion)
      write(izc,'(i2)') iz
      write(iionc,'(i2)') iion
      if (len(trim(adjustl(izc)))==1) then
        symbol_ion=adjustl(adjustr(symbol)//trim(adjustl(izc))//'+ ')
      else
        symbol_ion=adjustl(adjustr(symbol)//trim(adjustl(izc))//'+')
      endif
      A=adjustr(iionc)//' '//adjustl(symbol_ion)//adjustr(massc)//adjustr(Zc)//'  1'//adjustr(izc)&
                     //' '//adjustr(ispc)//' '//adjustr(ispc)//'  0  0  0  0'
      call write_string_to_file(unit_destin,A)      
      iion=iion+1
    enddo
  enddo

  nions=global_parameters%n_ions
  do imat=1,n_mat
    ipfc = materials(imat)%ipfc
    if (pfc_models(ipfc)%sputer_on .and. pfc_models(ipfc)%iatm > global_parameters%n_species) then
      nions=nions+1
      write(iawallc,'(i3)') pfc_models(ipfc)%iatm
      write(iiwallc,'(i2)') nions
      
      select case (pfc_models(ipfc)%material)
      case("Be")
        A=trim(adjustr(iiwallc))//' Be+       9  4  1  1'//trim(adjustr(iawallc))//trim(adjustr(iawallc))//'  0  0'
        call write_string_to_file(unit_destin,A)
      case("C ")
        A=trim(adjustr(iiwallc))//' C+       12  6  1  1'//trim(adjustr(iawallc))//trim(adjustr(iawallc))//'  0  0'
        call write_string_to_file(unit_destin,A)
      case("W")
        A=trim(adjustr(iiwallc))//' W+      184 74  1  1'//trim(adjustr(iawallc))//trim(adjustr(iawallc))//'  0  0'
        call write_string_to_file(unit_destin,A)
      end select
    endif
  enddo
  call copy_lines_until_finding(unit_source,'*** 6. DATA',unit_destin)
  call copy_file(unit_source,unit_destin,1)
 
end subroutine write_bulk_ions_species_cards_for_impurities

subroutine write_grid_parameters(unit_source,unit_destin)
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: unit_source,unit_destin
  integer :: NTRIAN
  character(6) :: NTTRA_styxc
  character(6) :: NTRIANc,NTorplanesc 
  character(8) :: NRKNOT_styxc
  character(12) :: LZc,ROA_styxc,ang_maxc
  character(128) :: A
 
  call copy_file(unit_source,unit_destin,4)
  call skip_line(unit_source,1)

  NTRIAN=NTRI_styx+1
  write(NTRIANc,'(i6)') NTRIAN
  write(NRKNOT_styxc,'(i8)') NRKNOT_styx

  A=NTRIANc//'     0     0     0'//NRKNOT_styxc//'     0'

  call write_string_to_file(unit_destin,A)

  A='CASE_'//trim(adjustl(fluid_code)) 
  call write_string_to_file(unit_destin,A)
  call skip_line(unit_source,1)

  if (levgeo_styx == 1) then
! cylindrical geometry in EIRENE, LZ = length of the cyclinder
    call copy_file(unit_source,unit_destin,7)
    call skip_line(unit_source,1)
! cm (LZ already in cm, see init_eirene and make_triangles)
    write(LZc,'(e12.4)') LZ
    A='  0.0000E 00  0.0000E 00'//LZc
    call write_string_to_file(unit_destin,A)
  
  elseif (levgeo_styx == 2) then
  
    if (is_3D) then
      call copy_file(unit_source,unit_destin,4)
      A='T'
      call write_string_to_file(unit_destin,A)
      call skip_line(unit_source,1)
    else
      call copy_file(unit_source,unit_destin,5)     
    endif
    A='FTF'
    call write_string_to_file(unit_destin,A)
    call skip_line(unit_source,1)
    write(ROA_styxc,'(e12.4)') ROA_styx

    if (is_3D) then
      write(Ntorplanesc,'(i6)') Ntor_cells+1
      write(ang_maxc,'(e12.4)') ang_max
      A=adjustr(Ntorplanesc)//'     0'//adjustr(Ntorplanesc)
      call write_string_to_file(unit_destin,A)
      call skip_line(unit_source,1)
      A='  0.0000E 00'//ang_maxc//ang_maxc//ang_maxc//ROA_styxc
      call write_string_to_file(unit_destin,A)
      call skip_line(unit_source,2)
    else
      write(NTTRA_styxc,'(i6)') NTTRA_styx
      A='     1     0'//NTTRA_styxc//'     0'
      call write_string_to_file(unit_destin,A)
      A='  0.0000E 00  0.0000E 00  3.6000E 02  3.6000E 02'//ROA_styxc
      call write_string_to_file(unit_destin,A)
      call skip_line(unit_source,2)
    endif
  endif
  call copy_lines_until_finding(unit_source,'*** 3A. DATA FOR NON-DEFAULT',unit_destin)

end subroutine write_grid_parameters

subroutine write_reflection_model_data(unit_source,unit_destin)
  use Mreaders
  use styx2eirene
  implicit none
  integer, intent(in) :: unit_source,unit_destin
  integer :: ipump,ipfc,icard,n_tor_periodicity
  character(2) :: inumc
  character(3) :: icardc
  character(6) :: iawallc
  character(5) :: Ntorc
  character(11) :: ZNMLc,sypc
  character(12) :: sycc,Rc

  n_tor_periodicity=0
  if (is_3D .and. ang_max < 360._dp) n_tor_periodicity=2
  write(unit_destin,'(i6)') 1+n_pumps+n_pfc_types+n_tor_periodicity
  call skip_line(unit_source,1)
  
! first the CEI (Core Edge Interface)
  call write_string_to_file(unit_destin,'* CORE EDGE INTERFACE (CEI) - REFLECTIVE WITH R=1e-30')
  call write_string_to_file(unit_destin,'     1     1     0')
  call write_string_to_file(unit_destin,'     1     0     0     0     0     4')
  call write_string_to_file(unit_destin,'     1     0     0     0')
  call write_string_to_file(unit_destin,' 1.8474E 04  3.5000E-02  0.0000E 00  0.0000E 00  0.0000E 00')
  call write_string_to_file(unit_destin,' 1.0000E 00  1.0000E-30  0.0000E 00')
  call skip_line(unit_source,6)
  icard=1
! now the pumps
  do ipump=1,n_pumps
    icard=icard+1
    write(inumc,'(i2)') ipump
    call write_string_to_file(unit_destin,'* PUMP #'//adjustr(inumc))
    write(inumc,'(i2)') icard
    call write_string_to_file(unit_destin,'    '//trim(adjustr(inumc))//'     1     0')
    call write_string_to_file(unit_destin,'     1     0     0     0     0     4')
    call write_string_to_file(unit_destin,'     1     0     0     0')
    call get_string_characterizing_material(pumps(ipump)%material,ZNMLc)    
    call write_string_to_file(unit_destin,adjustr(ZNMLc)//' -4.0000E-02  0.0000E 00  0.0000E 00  0.0000E 00')
    ! default value of R to first species index
    write(Rc,'(e12.4)') pumps(ipump)%R(1)
    call write_string_to_file(unit_destin,' 1.0000E 00'//adjustr(Rc)//'  0.0000E 00')
  enddo
  call skip_line(unit_source,6)
! and now the pfcs
  do ipfc=1,n_pfc_types
    icard=icard+1
    write(inumc,'(i2)') ipfc
    call write_string_to_file(unit_destin,'* PFC type #'//adjustr(inumc))
    write(inumc,'(i2)') icard
    call write_string_to_file(unit_destin,'    '//trim(adjustr(inumc))//'     1     0')
    call write_string_to_file(unit_destin,'     1     0     0     0     0     4')
    write(iawallc,'(i6)') pfc_models(ipfc)%iatm
    if (pfc_models(ipfc)%sputer_on) then
      if (pfc_models(ipfc)%material /= 'C ') then
        ! no chemical erosion
        if (pfc_models(ipfc)%sputer_model == 1) then
          call write_string_to_file(unit_destin,'     1     1'//trim(adjustr(iawallc))//'     0')
        else
          call write_string_to_file(unit_destin,'     1     2'//trim(adjustr(iawallc))//'     0')
        endif
        call get_string_characterizing_material(pfc_models(ipfc)%material,ZNMLc)
        call write_string_to_file(unit_destin,adjustr(ZNMLc)//' -4.0000E-02  0.0000E 00  0.0000E 00  0.0000E 00')
        call write_string_to_file(unit_destin,' 1.0000E 00  0.0000E 00  0.0000E 00')
        if (pfc_models(ipfc)%sputer_model == 1) then
          write(sypc,'(e11.4)') pfc_models(ipfc)%sputer_yield_phys
          call write_string_to_file(unit_destin,trim(adjustr(sypc))//'  0.0000E 00  0.0000E 00')
        elseif (pfc_models(ipfc)%sputer_model == 2) then
          call write_string_to_file(unit_destin,' 1.0000E 00  0.0000E 00  0.0000E 00')
        endif
      else
        ! carbon = chemical erosion (option A6 see Eirene Manual)
        if (pfc_models(ipfc)%sputer_model == 1) then
          call write_string_to_file(unit_destin,'     1    11'//trim(adjustr(iawallc))//trim(adjustr(iawallc)))
        else
          call write_string_to_file(unit_destin,'     1    22'//trim(adjustr(iawallc))//trim(adjustr(iawallc)))
        endif
        call get_string_characterizing_material(pfc_models(ipfc)%material,ZNMLc)
        call write_string_to_file(unit_destin,adjustr(ZNMLc)//' -4.0000E-02  0.0000E 00  0.0000E 00  0.0000E 00')
        call write_string_to_file(unit_destin,' 1.0000E 00  0.0000E 00  0.0000E 00')
        if (pfc_models(ipfc)%sputer_model == 1) then
          write(sypc,'(e11.4)') pfc_models(ipfc)%sputer_yield_phys
          write(sycc,'(e12.4)') pfc_models(ipfc)%sputer_yield_chem
          call write_string_to_file(unit_destin,trim(adjustr(sypc))//trim(adjustr(sycc))//' 0.0000E 00')
        elseif (pfc_models(ipfc)%sputer_model == 2) then
          call write_string_to_file(unit_destin,' 1.0000E 00  1.0000E 00  0.0000E 00')
        endif
      endif
    else
    ! no sputering    
      call write_string_to_file(unit_destin,'     1     0     0     0')
      call get_string_characterizing_material(pfc_models(ipfc)%material,ZNMLc)
      call write_string_to_file(unit_destin,adjustr(ZNMLc)//' -4.0000E-02  0.0000E 00  0.0000E 00  0.0000E 00')
      call write_string_to_file(unit_destin,' 1.0000E 00  0.0000E 00  0.0000E 00')
    endif
  enddo
  if (n_tor_periodicity == 2) then
    write(Ntorc,'(i5)') Ntor_cells+1
    icard=icard+1
    write(inumc,'(i2)') icard
    call write_string_to_file(unit_destin,'* TOROIDAL SURFACE NO. 1, periodicity with tor. surf. '//trim(adjustl(Ntorc)))
    call write_string_to_file(unit_destin,'    '//trim(adjustr(inumc))//'     3     1')
    !call write_string_to_file(unit_destin,'    3    -2     0     0     0     2     0     0')
    call write_string_to_file(unit_destin,adjustr(Ntorc)//'4    -2     0     0     0     2     0     0')
    call write_string_to_file(unit_destin,'* TOROIDAL SURFACE NO. '//trim(adjustl(Ntorc))//', periodicity with tor. surf. 1')
    icard=icard+1
    write(inumc,'(i2)') icard
    call write_string_to_file(unit_destin,'    '//trim(adjustr(inumc))//'     3 '//adjustr(Ntorc))
    call write_string_to_file(unit_destin,'    14     2     0     0     0     2     0     0')
    !call write_string_to_file(unit_destin,'     3     2     0     0     0     2     0     0')
  endif
  
  call skip_line(unit_source,6)
  call copy_lines_until_finding(unit_source,'* ATOMIC REACTION',unit_destin) 

end subroutine write_reflection_model_data

subroutine write_number_of_AM_processes(unit_source,unit_destin,nreac)
  use all_variables, only : global_parameters
  use styx2eirene, only : pfc_models,n_pfc_types, n_mat, materials
  implicit none
  integer, intent(in) :: unit_source,unit_destin,nreac
  integer :: ncx,isp,ireac_add,ipfc,imat

  ncx=0
  do isp=2,global_parameters%n_species
    ncx=ncx+global_parameters%element_list(isp)%ncx
  enddo
  ireac_add=0
  do imat=1,n_mat
    ipfc = materials(imat)%ipfc
    if (pfc_models(ipfc)%sputer_on .and. pfc_models(ipfc)%iatm > global_parameters%n_species)&
            ireac_add=ireac_add+1
  enddo
  write(unit_destin,'(i6)') nreac + (global_parameters%n_species-1)*4 +ncx + ireac_add
  call copy_file(unit_source,unit_destin,nreac)

end subroutine write_number_of_AM_processes


subroutine write_time_dependent_mode_data(unit_source,unit_destin,ical)
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: unit_source,unit_destin,ical
  character(128) :: A
  character(6) :: ntime_c,nprnlic,nptstc
  character(12) :: dtimvc

  if (ical==1) then
    call copy_file(unit_source,unit_destin,3)
    write(ntime_c,'(i6)') ntime_styx
    if (.not.timedep) then
      A='     1     1 10000     0     0     0     0     0'
    else
      A='     1     1 10000 30000     0     0     0'//trim(ntime_c)
    endif
    call skip_line(unit_source,1)
    call write_string_to_file(unit_destin,A)
  elseif (ical == 2) then
    if (timedep) then
      call skip_line(unit_source,1)
      write(nprnlic,'(i6)') nprnli_styx
      A=trim(nprnlic)
      call write_string_to_file(unit_destin,A)
      write(nptstc,'(i6)') nptst_styx
      A=trim(nptstc)//'     1'
      call write_string_to_file(unit_destin,A)
      write(dtimvc,'(e12.4)') dtimv_styx
      A=trim(dtimvc)//'  0.0000E 00'
      call write_string_to_file(unit_destin,A)
      call write_string_to_file(unit_destin,'** 13.A DATA FOR SNAPSHOT TALLIES')
      write(unit_destin,'(A6)') '     0'
    endif
    call copy_lines_until_finding(unit_source,'*** 14. DATA',unit_destin)
    call copy_file(unit_source,unit_destin,1)
  endif

end subroutine write_time_dependent_mode_data

subroutine write_tallies_to_be_calculated(unit_source,unit_destin)
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: unit_source,unit_destin
  integer :: NLTVOUT,NLTSOUT
  character(6) :: NLTVOUTc,NLTSOUTc
  character(128) :: A

  all_tal=.true.

  if (all_tal) then
! no tallies are switched off
    NLTVOUT=0
    write(NLTVOUTc,'(i6)') NLTVOUT
    A=NLTVOUTc
    write(unit_destin,'(A6)') A

  else

    if (direct_coupling) then
! for atoms, switch off all tallies except (9, 33, 38, 97),(15, 39, 44, 98),(21, 45, 50, 99) ; 100-12=88
      NLTVOUT=88
      write(NLTVOUTc,'(i6)') NLTVOUT
      A=NLTVOUTc
     
      call write_string_to_file(unit_destin,A)
      call write_string_to_file(unit_destin,'    -1    -2    -3    -4    -5    -6    -7    -8   -10   -11   -12   -13')
      call write_string_to_file(unit_destin,'   -14   -16   -17   -18   -19   -20   -22   -23   -24   -25   -26   -27')
      call write_string_to_file(unit_destin,'   -28   -29   -30   -31   -32   -34   -35   -36   -37   -40   -41   -42')
      call write_string_to_file(unit_destin,'   -43   -46   -47   -48   -49   -51   -52   -53   -54   -55   -56   -57')
      call write_string_to_file(unit_destin,'   -58   -59   -60   -61   -62   -63   -64   -65   -66   -67   -68   -69')
      call write_string_to_file(unit_destin,'   -70   -71   -72   -73   -74   -75   -76   -77   -78   -79   -80   -81')
      call write_string_to_file(unit_destin,'   -82   -83   -84   -85   -86   -87   -88   -89   -90   -91   -92   -93')
      call write_string_to_file(unit_destin,'   -94   -95   -96  -100')

    else
! for atoms, switch off all tallies except (1, 5, 85, 89, 93),(2, 6, 86, 90, 94),(3, 7, 87, 91, 95) 100-15=85
      NLTVOUT=95
      write(NLTVOUTc,'(i6)') NLTVOUT
      A=NLTVOUTc
 
      call write_string_to_file(unit_destin,A)
      call write_string_to_file(unit_destin,'    -4    -8   -9    -10   -11   -12   -13   -14   -15   -16   -17   -18')
      call write_string_to_file(unit_destin,'   -19   -20   -21   -22   -23   -24   -25   -26   -27   -28   -29   -30')
      call write_string_to_file(unit_destin,'   -31   -32   -33   -34   -35   -36   -37   -38   -39   -40   -41   -42')
      call write_string_to_file(unit_destin,'   -43   -44   -45   -46   -47   -48   -49   -50   -51   -52   -53   -54')
      call write_string_to_file(unit_destin,'   -55   -56   -57   -58   -59   -60   -61   -62   -63   -64   -65   -66')
      call write_string_to_file(unit_destin,'   -67   -68   -69   -70   -71   -72   -73   -74   -75   -76   -77   -78')
      call write_string_to_file(unit_destin,'   -79   -80   -81   -82   -83   -84   -88   -92   -96   -97   -98   -99')
      call write_string_to_file(unit_destin,'  -100')
    endif
   ! no surface tallies deactivated
   NLTSOUT=0
   write(NLTSOUTc,'(i6)') NLTSOUT
   A=NLTSOUTc
   write(unit_destin,'(A6)') A
  endif
  call skip_line(unit_source,9)
  call copy_lines_until_finding(unit_source,'*** 13. DATA',unit_destin)
end subroutine write_tallies_to_be_calculated

subroutine write_dummy_species_sampling_distributions(unit_source,unit_destin)
  use all_variables, only : global_parameters
  use Mreaders
  implicit none
  integer, intent(in) :: unit_source,unit_destin
  integer, parameter :: nmax_speciesx12=1200
  character(nmax_speciesx12) :: DATDc,DMLDc,DIODc,DPLDc
  integer :: iatm,imol,iion,ipls  

  DATDc = '  0.0000E 00'
  do iatm=1,min(global_parameters%n_species-1,6)
    DATDc=trim(DATDc)//'  0.0000E 00' 
  enddo
  call write_string_to_file(unit_destin,DATDc)

  DMLDc = '  0.0000E 00'

  ! max number set to 6 (format 6664 in input)
  do imol=1,6
    DMLDc=trim(DMLDc)//'  0.0000E 00' 
  enddo
  call write_string_to_file(unit_destin,DMLDc)

  ! max number set to 10 (arbitrarily)
  DIODc = '  0.0000E 00'
  do iion=1,6
    DIODc=DIODc//'  0.0000E 00' 
  enddo
  call write_string_to_file(unit_destin,DIODc)
  
  DPLDc = '  0.0000E 00'
  do ipls=1,min(global_parameters%n_ions-1,6)
    DPLDc=trim(DPLDc)//'  0.0000E 00' 
  enddo
  call write_string_to_file(unit_destin,DPLDc)
  call skip_line(unit_source,4)
  call copy_lines_until_finding(unit_source,'*** 7. DATA',unit_destin)
  
end subroutine write_dummy_species_sampling_distributions

subroutine write_data_for_ion_resolved_energy_sources(unit_source,unit_destin)
  use all_variables, only : global_parameters
  use styx2eirene
  use Mreaders
  implicit none
  integer,intent(in) :: unit_source,unit_destin
  integer :: nadv,i,imat,ipfc
  character(3) :: nadvc
  character(128) :: A
  call skip_line(unit_source,1)
  ispc_add=0
  do imat=1,n_mat
    ipfc = materials(imat)%ipfc
    if (pfc_models(ipfc)%sputer_on .and. pfc_models(ipfc)%iatm > global_parameters%n_species)&
            ispc_add=ispc_add+1
  enddo

  nadv=3*(global_parameters%n_ions+ispc_add)
  write(nadvc,'(i3)') nadv
  A='   '//nadvc//'     0     0     0     0     0'
  call write_string_to_file(unit_destin,A)
  call copy_file(unit_source,unit_destin,1)

  do i=1,global_parameters%n_ions+ispc_add
    A='     3    -1    -1     1'
    call write_string_to_file(unit_destin,A)    
    A='     fatm'
    call write_string_to_file(unit_destin,A)  
    A='     fatm                 fatm'
    call write_string_to_file(unit_destin,A)  
  enddo
  do i=1,global_parameters%n_ions+ispc_add
    A='     3    -1    -1     2'
    call write_string_to_file(unit_destin,A)
    A='     fmol'
    call write_string_to_file(unit_destin,A)
    A='     fmol                 fmol'
    call write_string_to_file(unit_destin,A)  
  enddo
  do i=1,global_parameters%n_ions+ispc_add
    A='     3    -1    -1     3'
    call write_string_to_file(unit_destin,A)
    A='     fion'
    call write_string_to_file(unit_destin,A)  
    A='     fion                 fion'
    call write_string_to_file(unit_destin,A)  
  enddo
  call skip_line(unit_source,1)
  call copy_lines_until_finding(unit_source,'*** 11. DAT',unit_destin)
end subroutine


subroutine write_strata_data(unit_source,unit_destin,ical)
  use eirmod_precision
  use styx2eirene
  use Mreaders
  implicit none
  integer, intent(in) :: unit_source,unit_destin,ical
  integer :: i
  character(NSTRATA+1) :: TRCSRCc
  character(NSTRATA*6) :: INDSRCc
  character(6) :: NSTRATAc
  character(128) :: A
  if (ical==2) then
    call copy_file(unit_source,unit_destin,1)  
    TRCSRCc='T'
    do i=1,NSTRATA
      TRCSRCc=trim(TRCSRCc)//'T'
    enddo
    call write_string_to_file(unit_destin,TRCSRCc) 
    call skip_line(unit_source,1)
    call copy_file(unit_source,unit_destin,2)
  elseif (ical ==1) then
    call skip_line(unit_source,1)
    write(NSTRATAc,'(i6)') NSTRATA
    A=NSTRATAc
    call write_string_to_file(unit_destin,A)
    INDSRCc='     6'
    do i=1,NSTRATA-1
      INDSRCc=trim(INDSRCc)//'     6'
    enddo
    call write_string_to_file(unit_destin,INDSRCc)
    call skip_line(unit_source,1)
    call copy_lines_until_finding(unit_source,'*** 10. DAT',unit_destin)
  endif
end subroutine write_strata_data


subroutine write_data_for_spatially_resolved_tallies(unit_source,unit_destin)
  use Mreaders
  use styx2eirene
  implicit none
  integer,intent(in) :: unit_source,unit_destin
  integer :: NRTAL,NGSTAL
  character(6) :: NRTALc,NGSTALc
  character(128) :: A
  NRTAL=(NTRI_styx+1)*(Ntor_cells+1)
  if (NRTAL > 999999) then
      write(*,*) 'Warning, formating issue with NRTAL > 999999 if it were used'
      !write(*,*) 'This is not yet implemented, too many cells'
  endif
  call skip_line(unit_source,1)
! spatially resolved surface tallies NGSTAL=1 (otherwise 0)
  NGSTAL=1
  write(NRTALc,'(i6)') NRTAL
  write(NGSTALc,'(i6)') NGSTAL
  ! NRTAL -> mesh condensation, may be usefull in the future : deactivate now for safety (formatting issue)
  !A='     1     1     0     1     1     9'//trim(NGSTALc)//trim(NRTALc)//'     0'     !trim(NREAC_ADDc)
  A='     1     1     0     1     1     9'//trim(NGSTALc)//'     0     0'     !trim(NREAC_ADDc)
  call write_string_to_file(unit_destin,A)   
  call copy_file(unit_source,unit_destin,1)
end subroutine

subroutine write_trim_data_needed(unit_source,unit_destin)
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comprt, only : iunout
  use styx2eirene
  implicit none
  integer, intent(in) :: unit_source,unit_destin
  integer :: isp,imat
  character(128) :: A
  character(8) :: trim_file
  logical :: file_exist

  ! path to TRIM database (must be copied locally)
  A='path = TRIM/'
  call write_string_to_file(unit_destin,A)  

  do imat=1,n_mat
    do isp=1,global_parameters%n_species
      trim_file=trim(adjustl(global_parameters%element_list(isp)%symbol))//'_on_'//&
           trim(adjustl(materials(imat)%symbol))
      inquire(file='TRIM/'//trim_file,exist=file_exist)
      if (file_exist) then     
        call write_string_to_file(unit_destin,trim_file)
      else
        write(iunout,*) 'Warning, trim card '//trim_file//' not found, using scaling instead ...'
      endif
    enddo
  enddo

  ! for pumps on Fe : not all TRIM data available, e.g. Be_on_Fe
  !do isp=1,global_parameters%n_species
  !trim_file=trim(adjustl(global_parameters%element_list(isp)%symbol))//'_on_'//&
  !         trim(adjustl(pumps(1)%material))
  !    call write_string_to_file(unit_destin,trim_file)
  !enddo


end subroutine


subroutine get_string_characterizing_material(symbol,ZNMLc)
  use eirmod_comprt, only : iunout
  implicit none
  character(2), intent(in) :: symbol
  character(11), intent(out) :: ZNMLc

  select case (symbol)
    case("Be")
      ZNMLc=' 9.0400E 02'
    case("C ")
      ZNMLc=' 1.2060E 03'
    case("Mo")
      ZNMLc=' 9.6420E 03'
    case("W ")
      ZNMLc=' 1.8474E 04'
    case("Fe")
      ZNMLc=' 5.6260E 03'
    case default
      write(iunout,*) 'Error from get_string_characterizing_wall_material() check material specification in input file'
      write(iunout,*) 'Material should be either Be, C, Fe, Mo, W'
      call eirene_exit_own(1)
  end select

end subroutine




