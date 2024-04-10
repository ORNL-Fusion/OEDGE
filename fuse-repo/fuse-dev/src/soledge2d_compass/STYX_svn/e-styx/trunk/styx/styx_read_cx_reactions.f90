subroutine styx_read_cx_reactions
  use all_variables, only : global_parameters
  use styx2eirene
  implicit none
  integer :: isp,ncx,iion,chrg,icx,iatm_cx
  integer :: iatm_p1,iion_p2,n_species,iatm
  integer :: ispz,iplus,chrg_cx,iion_cx
  integer :: ispz_cx,zion_cx,i
  logical, allocatable :: active_channel(:)
  real*8, allocatable :: ftweak(:)
  real*8, dimension(9) :: fit_coeffs
  character(4) :: symbol_at
  character(5), allocatable :: symbol_ion(:)
  character(11):: buff11
  character(1) :: chrgc
  character(2) :: ncxc,symb,iatm_p1c,iatm_p2c,iion_cxc
  character(20) :: format
  character(8) :: reac_tag
  character(29) :: cx_database1,cx_database2
  character(128) :: string

  ! initialize to zero (summed later on)
  do isp=1,global_parameters%n_species
    global_parameters%element_list(isp)%ncx=0
  enddo
 
  open(unit=666,file='cx_setup',status='old')
  ! number of species for which cx reactions are specified
  read(666,'(A11,i2)') buff11,n_species

  if (n_species > global_parameters%n_species) then
    write(*,*) ' two many species in cx_setup, check settings '
    stop
  endif

  do isp=1,n_species
    read(666,'(A11)') buff11
    read(666,'(A9,A4)') buff11(1:9),symbol_at
    read(666,'(A5,i2)') buff11(1:5),ncx
    
    allocate(symbol_ion(ncx))
    allocate(active_channel(ncx))
    allocate(ftweak(ncx))
    write(ncxc,'(i2)') ncx

    read(666,*)
    read(666,'(A128)') string
    call parse_line_elements_cx(string,ncx,symbol_ion)

    ! look for atomic species number

    symbol_at=adjustl(symbol_at)
    symb=symbol_at(1:2)

    iatm_cx=0
    iatm=1
    do while (iatm_cx==0 .and. iatm <= global_parameters%n_species)
      if (symb == global_parameters%element_list(iatm)%symbol) then
        iatm_cx=iatm
      endif
      iatm=iatm+1
    enddo

    if (iatm_cx == 0) then
      write(*,*) 'cx reactions, atomic species not found ...'
      stop 
    endif

    global_parameters%element_list(iatm_cx)%ncx=ncx

    format=trim(adjustl('(A8,'//ncxc//'L1)'))
    read(666,format) buff11(1:7),(active_channel(i), i=1,ncx)
    format=trim(adjustl('(A6,'//ncxc//'e12.4))'))
    read(666,format) buff11(1:7),(ftweak(i), i=1,ncx)

    allocate(global_parameters%element_list(iatm_cx)%cxdatas(1:ncx))

    ! look for the corresponding ion for all the cx channels
    do icx=1,ncx
      if (active_channel(icx)) then

        ! look for +, get charge
        iplus = index(symbol_ion(icx),'+')
        if (iplus == 0) then
          write(*,*) ' string + not found reading species symbol ',symbol_ion(icx)
          stop
        endif
        chrgc = symbol_ion(icx)(iplus-1:iplus-1)
        if (chrgc== '1' .or. chrgc == '2' .or. chrgc == '3' .or. &
          chrgc == '4' .or. chrgc == '5' .or. chrgc == '6' .or. &
          chrgc == '7' .or. chrgc == '8' .or. chrgc == '9') then
          read(chrgc,'(i1)') chrg_cx
        else
          chrg_cx = 1
        endif

        symbol_ion(icx)=adjustl(symbol_ion(icx))
        if (symbol_ion(icx)(2:2) /=' ' .and. symbol_ion(icx)(2:2) /= '+') then
          symb=symbol_ion(icx)(1:2)
        else
          symb=symbol_ion(icx)(1:1)
        endif   

        iion_cx=0
        do iion=1,global_parameters%n_ions
          ispz=global_parameters%ions_list(iion,1)
          chrg=global_parameters%ions_list(iion,2)
          if (symb == global_parameters%element_list(ispz)%symbol&
                   .and. chrg_cx == chrg ) then
            iion_cx=iion
            ispz_cx=ispz
            zion_cx=chrg
          endif 
        enddo
        if (iion_cx==0) then
          write(*,*) ' element: ',symb,' with charge ',chrg_cx,'not found in species list'
          stop
        endif
        ! iion_cx is now the ion number looked for, and has charge zion_cx
        ! now look for productsi p1,p2
  
        ! neutral, same species as ion iion_cx
        iatm_p1=ispz_cx

        ! ion, same charge as ion iion_cx, same species as reacting atom
        do iion=1,global_parameters%n_ions
          ispz=global_parameters%ions_list(iion,1)
          chrg=global_parameters%ions_list(iion,2)
          if (ispz == iatm_cx .and. chrg == zion_cx) then
            iion_p2=iion
          endif
        enddo
      
        write(iatm_p1c,'(i2)') iatm_p1
        write(iatm_p2c,'(i2)') iion_p2
        write(iion_cxc,'(i2)') iion_cx
 
        !!! fill amdata structures: reference of reactions in AMJUEL + ion spcies + tweaks !!!

        ! first reac data
        reac_tag=trim(adjustl(symbol_at))//'-'//trim(adjustl(symbol_ion(icx)))
        call styx_get_cx_reaction_data(reac_tag,cx_database1,cx_database2,fit_coeffs)

        global_parameters%element_list(iatm_cx)%cxdatas(icx)%cx_database1=cx_database1
        global_parameters%element_list(iatm_cx)%cxdatas(icx)%cx_database2=cx_database2
        if (.not.direct_coupling) then        
          global_parameters%element_list(iatm_cx)%cxdatas(icx)%fit_coeffs = fit_coeffs
        endif
        global_parameters%element_list(iatm_cx)%cxdatas(icx)%ibulk=trim(adjustl(iion_cxc))//'14'
        global_parameters%element_list(iatm_cx)%cxdatas(icx)%iscd1=trim(adjustl(iatm_p1c))//'11'
        global_parameters%element_list(iatm_cx)%cxdatas(icx)%iscd2=trim(adjustl(iatm_p2c))//'14'
        global_parameters%element_list(iatm_cx)%cxdatas(icx)%iscde=' 01001'
      endif
    enddo
    deallocate(symbol_ion,active_channel,ftweak)    
  enddo
  close(666)

end subroutine styx_read_cx_reactions


subroutine  parse_line_elements_cx(String,ncx,symbol_ion)
    use all_variables, only : global_parameters
    use Mphysics
    use styx2eirene, only : puffs
    implicit none
    character(128),intent(in) :: String
    integer*4,intent(in) :: ncx
    character(5), dimension(ncx), intent(out) :: symbol_ion
    character(128) :: buffer
    integer :: pos
    integer*4 i
    character(5) :: dummy_symbol
    buffer=String
    do i=1,ncx-1
       pos=index(buffer,",")
       if(pos.ge.5) then
          dummy_symbol=adjustl(trim(buffer(1:pos-1)))
          symbol_ion(i)=dummy_symbol
          buffer=adjustl(buffer(pos+1:))
       else
          print*,'Unrecognised species list for cx reactions'
          stop
       end if
    end do
    !last species
    dummy_symbol=adjustl(trim(buffer))
    symbol_ion(ncx)=dummy_symbol
  end subroutine parse_line_elements_cx



