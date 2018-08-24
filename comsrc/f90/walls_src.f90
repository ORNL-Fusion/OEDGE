module walls_src

  ! This module contains code and variables related to the DIVIMP option allowing 
  ! for impurity data from previous DIVIMP runs to be used as the sputter source
  ! the current case

  ! Impurity source and energy by charge state and wall element
  ! Source is scaled using the ABSFAC from the previous DIVIMP run

  real :: ws_totals(6)
  real,allocatable :: ws_impflux(:,:),ws_impeng(:,:),ws_impyield(:,:),ws_impfluxyield(:,:)
  real,allocatable :: ws_impeng_dist
  integer,allocatable :: ws_targid(:),ws_wallid(:)
  integer :: ws_targcnt,ws_maxtargcnt

  real*8  :: ws_absfac_neut,ws_absfac_ion
  real*8  :: ws_crmi
  integer :: ws_nizs,ws_cion,ws_nwall



contains



  subroutine read_resolved_deposition_data
    use error_handling
    implicit none
    !      include 'params'
    !      include 'cgeom'
    !      include 'comtor'
    !      include 'div1'
    !      include 'dynam3'
    !
    !integer :: nizs
    !
    !
    integer :: cnt
    real :: r,z,length,ws_eroded,ws_ionized, ws_iondep, ws_neutdep, ws_ionleak
    real,allocatable :: charge(:),flux(:),energy(:)

    integer :: ounit,ierr,iz,in,targid,data_start,taglen
    character*1024 :: fname,line,tag

    fname = 'charge_resolved_deposition_data_in.dat'
    call find_free_unit_number(ounit)

    open(unit=ounit,file=trim(fname),iostat=ierr)

    if (ierr.ne.0) then 
       call errmsg('READ_RESOLVED_DEPOSITION_DATA:PROBLEM OPENING FILE = '//trim(fname),ierr)
    else   
       ! Write out the charge resolved data
       ! Include absfac as well as atomic number and mass of impurity
       ! Include data showing the impurity flux as a fraction of the 
       ! H flux (?) 

       !
       !       Read tagged input file
       !
       ierr = 0

       do while (ierr.eq.0) 

          read(unit=ounit,fmt='(a1024)',iostat=ierr) line

          if (ierr.eq.0) then 

             ! ignore comments
             if (line(1:1).eq.'#') cycle

             call get_tag(line,tag,data_start)
             taglen=len(trim(tag))

             if (trim(tag).eq.'ABSFAC_ION') then 
                read(line(data_start:),*) ws_absfac_ion
             elseif (trim(tag).eq.'ABSFAC_NEUT') then 
                read(line(data_start:),*) ws_absfac_neut
             elseif (trim(tag).eq.'ATOMIC MASS') then 
                read(line(data_start:),*) ws_crmi
             elseif (trim(tag).eq.'ATOMIC NUMBER') then 
                read(line(data_start:),*) ws_cion
             elseif (trim(tag).eq.'CHARGE STATES') then 
                read(line(data_start:),*) ws_nizs
             elseif (trim(tag).eq.'N WALL') then 
                read(line(data_start:),*) ws_nwall
             elseif (trim(tag).eq.'TOTALS') then 
                read(line(data_start:),*) (ws_totals(in),in=1,6)
             elseif (trim(tag).eq.'DATA') then 
                ! allocate storage for data
                ! 1 - wall elements
                ! 2 - charge states
                ! 3 - flux of charge state (particles /m2/s received onto surface, average energies
                call allocate_imp_data(ws_nwall,ws_nizs)

                ! allocate local storage
                if (allocated(charge)) deallocate(charge)
                allocate(charge(ws_nizs))

                if (allocated(flux)) deallocate(flux)
                allocate(flux(ws_nizs))

                if (allocated(energy)) deallocate(energy)
                allocate(energy(ws_nizs))

                do in = 1,ws_nwall
                   read(ounit,'(1x,i8,3(1x,f15.7),i8,512(1x,g18.8))') &
                        cnt, r, z, length, targid, ws_eroded, ws_ionized, ws_iondep, ws_neutdep, &
                        (charge(iz),flux(iz),energy(iz),iz=1,ws_nizs), ws_ionleak

                   ws_targid(in) = targid
                   do iz = 1,ws_nizs
                      ws_impflux(in,iz) = flux(iz)
                      ws_impeng(in,iz) =  energy(iz)
                   end do

                end do

                deallocate(charge)
                deallocate(flux)
                deallocate(energy)

                ! Set ierr =1 to exit loop 
                ierr = 1

             endif

          endif

       end do

       close(ounit)


    end if




    return

  end subroutine read_resolved_deposition_data



  subroutine calc_external_source(maxnds,matp,matt)
    use error_handling
    implicit none
    integer :: matp,matt,maxnds

    ! analyse the external flux data
    ! calculate the yields and effective yields for each wall element
    ! Note that the data read in IS in terms of wall elements and not
    ! target elements ... this will require a mapping back to target data
    ! when the code is asking for target vs wall particle sources
    ! Some of this geometry data may be needed in the data file

    ! Also need to tie into appropriate sputtering models
    ! MATP, MATT ... MATP is set in the input file - needs to match
    ! the species from the DIVIMP run data used


    real,external :: yield
    integer :: in,iz,id
    ! Calculate the total fluxes, individual yields and effective yields for each wall element
    ! zero initial arrays

    if (maxnds.gt.0.and.(.not.allocated(ws_wallid))) then 
       ! allocate and calculated ws_wallid which contains the 
       ! wallindex of the corresponding target index
       if (allocated(ws_wallid)) deallocate(ws_wallid)
       ws_maxtargcnt= maxnds
       allocate(ws_wallid(ws_maxtargcnt))
       do id = 1,ws_nwall
          if (ws_targid(id).gt.0.and.ws_targid(id).le.ws_maxtargcnt) then
             ! create cross index of target -> wall indices
             ws_wallid(ws_targid(id)) = id
          endif
       end do
    endif

    ! calculate the total fluxes, yields, effective yields and average energies
    do in = 1,ws_nwall
       do iz = 1,ws_nizs
          if (ws_impflux(in,iz).gt.0.0) then
             ! non-zero flux for segment and charge ... calculate yield and flux*yield and totals
             ws_impyield(in,iz) = YIELD (MATP,MATT,ws_impeng(in,iz),0.0,0.0)
             ws_impfluxyield(in,iz) = ws_impflux(in,iz) * ws_impyield(in,iz)
             !
             ! accumulate flux weighted totals to get flux weighted average
             !
             ws_impfluxyield(in,ws_nizs+1) = ws_impfluxyield(in,ws_nizs+1) + ws_impfluxyield(in,iz)
             ws_impflux(in,ws_nizs+1) = ws_impflux(in,ws_nizs+1) + ws_impflux(in,iz)
             ws_impeng(in,ws_nizs+1)   = ws_impeng(in,ws_nizs+1) + ws_impeng(in,iz) * ws_impflux(in,iz)
             ws_impyield(in,ws_nizs+1) = ws_impyield(in,ws_nizs+1) + ws_impyield(in,iz) * ws_impflux(in,iz)
             !
          endif
       end do
       ! normalize accumulators
       if (ws_impflux(in,ws_nizs+1).gt.0.0) then 
          ! average energy for all charge states
          ws_impeng(in,ws_nizs+1)   = ws_impeng(in,ws_nizs+1)/ws_impflux(in,ws_nizs+1)
          ! effective yield over all charge states
          ws_impyield(in,ws_nizs+1) = ws_impyield(in,ws_nizs+1)/ ws_impflux(in,ws_nizs+1)
       else
          ws_impeng(in,ws_nizs+1) = 0.0
          ws_impyield(in,ws_nizs+1) = 0.0
       endif

    end do

    return
  end subroutine calc_external_source



  subroutine assign_deposition_data(fydata,maxpts,maxdat,nt,nw)
    use error_handling
    implicit none
    ! nt and nw are used to determine whether target or wall indexing should be used
    integer :: maxpts,nt,nw,maxdat
    real :: fydata(maxpts,maxdat)
    integer :: in,it

    ! assign external flux source data to the divimp fydata arrays
    ! NOTE: sputter option 7 will need to call the get_imp_energy routine to 
    ! obtain a randomly selected impurity impact energy based on the averages for
    ! each charge state and distributed proportional to flux

    ! sputtering yield multipliers are applied in the calling routine
    fydata = 0.0

    ! Calculate data for target source
    if (nw.gt.0.and.ws_nwall.ne.nw) then 
       call errmsg('ERROR: WALLS_SRC.F90: ASSIGN_DEPOSITION_DATA: MISMATCH IN NUMBER OF WALL ELEMENTS:',ws_nwall)
       stop 'ASSIGN_DEPOSITION_DATA'
    endif

    do in = 1,ws_nwall
       ! select target or wall source
       if (nt.gt.0) then 
          it = ws_targid(in)
       else
          it = in
       endif
       if (it.gt.0.and.it.le.nt) then 
          fydata(it,1) = ws_impflux(in,ws_nizs+1) ! flux
          fydata(it,2) = ws_impeng(in,ws_nizs+1)  ! energy
          fydata(it,3) = ws_impeng(in,ws_nizs+1)  ! heat = energy for now
          fydata(it,4) = ws_impyield(in,ws_nizs+1) ! effective yield
          fydata(it,5) = ws_impfluxyield(in,ws_nizs+1) ! total flux*yield for the element
       endif
    end do

    return
  end subroutine assign_deposition_data




  subroutine get_tag(line,tag,data_start)
    implicit none
    character*(*) :: line,tag
    integer :: tag_start,tag_end,data_start

    ! tag limits are defined by the first { and last } ... these should ONLY be used as tag delimiters ... the text between is the tag
    tag_start = index(line,'{')
    tag_end = index(line,'}',back=.true.)
    data_start = 0
    tag = ''

    ! error conditions ... either the delimiters are not found, they are in the wrong order or the tag string is empty
    if (tag_start.eq.0.or.tag_end.eq.0.or.(tag_start+1).gt.(tag_end-1)) return

    data_start = tag_end+1
    tag = line(tag_start+1:tag_end-1)

    return
  end subroutine get_tag


  subroutine allocate_imp_data(nwall,nizs)
    use error_handling
    implicit none
    integer :: nwall,nizs
    integer :: ierr

    ierr = 0
    if (.not.allocated(ws_impflux)) then 
       allocate(ws_impflux(nwall,nizs+1),stat=ierr)
       ws_impflux = 0.0
    endif
    if (.not.allocated(ws_impeng)) then 
       allocate(ws_impeng(nwall,nizs+1),stat=ierr)
       ws_impeng = 0.0
    endif
    if (.not.allocated(ws_impyield)) then 
       allocate(ws_impyield(nwall,nizs+1),stat=ierr)
       ws_impyield = 0.0
    endif
    if (.not.allocated(ws_impfluxyield)) then 
       allocate(ws_impfluxyield(nwall,nizs+1),stat=ierr)
       ws_impfluxyield = 0.0
    endif
    if (.not.allocated(ws_targid)) then 
       allocate(ws_targid(nwall),stat=ierr)
       ws_targid = 0
    endif

    if (ierr.ne.0) then
       call errmsg('ERROR:ALLOCATE_IMP_DATA: ERROR ALLOCATING SPACE ERR =',ierr)
       stop 'ALLOCATE_IMP_DATA'
    endif

    return
  end subroutine allocate_imp_data

  subroutine deallocate_imp_data
    use error_handling
    implicit none
    integer :: ierr

    if (allocated(ws_impflux)) deallocate(ws_impflux)
    if (allocated(ws_impeng)) deallocate(ws_impeng)
    if (allocated(ws_impyield)) deallocate(ws_impyield)
    if (allocated(ws_impfluxyield)) deallocate(ws_impfluxyield)
    if (allocated(ws_targid)) deallocate(ws_targid)
    if (allocated(ws_wallid)) deallocate(ws_wallid)

    return
  end subroutine deallocate_imp_data



  subroutine get_imp_energy(energy,id)
    use error_handling
    implicit none
    real :: energy
    integer :: id
    real,allocatable :: energy_dist(:,:) 
    integer :: iz
    real :: ran
    real,external :: getranf

    ! particle has been sputtered from this wall element
    ! find the impact energy based on the average energies of the charge states striking this element

    if (ws_impfluxyield(id,ws_nizs+1).le.0.0) then 
       call errmsg('ERROR: WALLS_SRC.F90: GET_IMP_ENERGY: NO YIELD FOR SPECIFIED WALL ELEMENT:',id)
       energy = 0.0
       return
    endif

    allocate(energy_dist(ws_nizs+1,2))

    energy = 0.0
    energy_dist = 0.0

    do iz = 1,ws_nizs
       if (ws_impfluxyield(id,iz).gt.0.0) then 
          energy_dist(iz,1) = ws_impeng(id,iz)
          ! generate a cumulative probability distribution - the last bin should have the same value as the total f*y so when normalized will be 1.0
          if (iz.eq.1) then 
             energy_dist(iz,2) = ws_impfluxyield(id,iz)
          else
             energy_dist(iz,2) = energy_dist(iz-1,2) + ws_impfluxyield(id,iz)
          endif
          energy_dist(ws_nizs+1,2) = energy_dist(ws_nizs+1,2) +  ws_impfluxyield(id,iz)
       endif
    end do


    if (energy_dist(ws_nizs+1,2).le.0.0) then ! This should be a bug if it ever happens
       call errmsg('ERROR: WALLS_SRC.F90: GET_IMP_ENERGY: NO CHARGE STATES WITH NON-ZERO YIELD:',id)
       energy = 0.0
       return
    endif

    do iz = 1, ws_nizs
       energy_dist(iz,2) = energy_dist(iz,2)/energy_dist(ws_nizs+1,2)
    end do

    ran = getranf()

    ! find the energy based on the random number

    do iz = 1,ws_nizs
       if (energy_dist(iz,1).le.0.0) cycle
       if (ran.le.energy_dist(iz,2)) then 
          ! found energy based on random number
          energy = energy_dist(iz,1)
          deallocate(energy_dist)
          return
       endif
    end do

    return
  end subroutine get_imp_energy



end module walls_src
