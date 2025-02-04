module ero_interface
  use error_handling

  implicit none

  !
  ! Modeule: ERO_INTERFACE
  !
  ! This module is intended to include all of the interface code for ERO related options.
  !
  ! 1) Background plasma transfer to ERO
  ! 2) Statistics gathering for impurity fluxes to pass to ERO for more detailed modeling
  ! 3) Launch option supporting DIVIMP particle source based on ERO output
  !
  !
  
  private

  integer,parameter :: line_length = 512
  character*6 ,parameter :: line_form='(a512)'
  character*(line_length) :: line


  integer :: ero_opt
  integer :: n_ero_regs

  integer :: ero_export_bg
  integer :: ero_bg_interp_opt
  integer :: ero_bg_extrap_opt
  integer :: ero_shape_opt

  !integer :: ero_record_particle_opt
  integer :: ero_particle_launch_opt
  integer :: ero_remove_src_opt       ! option to turn off DIVIMP surface sources within the ERO volume - defaul off

  integer :: ero_default_nx = 100
  integer :: ero_default_ny = 100

  integer :: erodiv_tor_part_mult


  real*8,allocatable :: ero_coords(:,:,:)  ! ([R,Z],vert,nblock)
  integer,parameter :: ero_nvert = 4          ! immediate support for only one ERO volume
  real :: ero_rvert(ero_nvert),ero_zvert(ero_nvert)
  real*8,allocatable :: ero_offsets(:,:)   ! ([R_off, Z_off], nblock)
  real*8,allocatable :: ero_resolution(:,:) ! ([dx,dy],nblock)
  integer,allocatable :: ero_ncells(:,:) ! ([nx,ny],nblock)
  real*8,allocatable :: ero_tor(:)     ! +/- Toroidal extent in meters for ERO volume


  character*256 :: ero_bg_output_fn_base
  character*256 :: ero_shape_output_fn_base

  character*256 :: ero_bg_command_fn

  !
  ! Variables related to particle track data collection for passing to ERO
  !

  character*256 :: ero_particle_output_fn,ero_particle_transfer_fn
  integer :: ero_part_output_opt
  integer,parameter :: ero_part_output_unit =59
  logical :: ero_in_vol,ero_last_in_vol

  !
  ! Variables related to particle transfer from ERO and calculation of launch source
  !
  
  character*256 :: ero_particle_input_fn
  integer,parameter :: n_ero_data = 12   ! transfer 13 data values : in, charge, r, z, t, vr, vz, vt, v0, energy, atoms, probabilty, cumulative probability
                                         ! last data value is cumulative probability and it is stored separately  
  integer,parameter :: n_extra_data =1   ! record which cells contribute particles

  integer :: n_ero_part
  real*8 :: ero_totsrc,ero_totsrc_raw,ero_totsrc_mtor,ero_src_tor_extent
  real*8, allocatable :: ero_part_data(:,:)
  real, allocatable :: ero_launch_prob(:)

  public write_ero_bg,run_ero_bg,ero_cleanup,read_ero_unstructured_input,init_ero_unstructured_input,load_ero_launch_data,&
       & select_ero_particle,get_ero_part_data,get_ero_part_data2,get_ero_src,init_ero_part_output,close_ero_part_output,ero_init_part_track, &
       & ero_check_part_track_ion, ero_check_part_track_neut,ero_remove_src, ero_opt, ero_part_output_opt, ero_particle_launch_opt,ero_remove_src_opt


contains


  subroutine write_ero_bg
    implicit none

    ! This routine writes an instruction file that is read in by an external 
    ! program which transcribes the background plasma into the appropriate 
    ! format including geometry transformations and extrapolation/interpolation 
    ! if specified. 

    integer :: fnum,in,iv
    character*512 :: ero_bg_output_fn
    character*512 :: ero_shape_output_fn


    call find_free_unit_number(fnum)

    open(unit=fnum,file=trim(ero_bg_command_fn),form='formatted')

    write(0,*)'WRITE_ERO_BG:',trim(ero_bg_command_fn),':'


    write(fnum,'(a,1x,2(1x,i10))') '#NUMBLOCKS ',n_ero_regs
    write(fnum,'(a,1x,a)') '#DIVPLASMAFILE ','divimp_plasma.out'
    write(fnum,'(a,1x,a)') '#DIVGEOFILE ','divimp_grid.out'


    do in = 1,n_ero_regs

       ! IF there is just one ERO region then use the base names as the only name
       if (n_ero_regs.eq.1) then
          ero_bg_output_fn =trim(ero_bg_output_fn_base)//'.m'
          ero_shape_output_fn = trim(ero_shape_output_fn_base)//'.m'
       else
          write(ero_bg_output_fn,'(a,i0.2,a)') trim(ero_bg_output_fn_base),in,'.m'
          write(ero_shape_output_fn,'(a,i0.2,a)') trim(ero_shape_output_fn_base),in,'.m'
       endif

       write(fnum,'(a)') '#EROBLOCK'
       write(fnum,'(a,1x,2(1x,g18.8))') '#OFFSETS',ero_offsets(1,in),ero_offsets(2,in)
       write(fnum,'(a,1x,4(1x,g18.8))') '#RVERTICES',(ero_coords(1,iv,in),iv = 1,4)
       write(fnum,'(a,1x,4(1x,g18.8))') '#ZVERTICES',(ero_coords(2,iv,in),iv = 1,4)
       write(fnum,'(a,1x,4(1x,g18.8))') '#TOR_EXTENT',ero_tor(in)
       write(fnum,'(a,1x,2(1x,g18.8))') '#DELTAXY',ero_resolution(1,in),ero_resolution(2,in)
       write(fnum,'(a,1x,2(1x,i10))') '#NXNY',ero_ncells(1,in),ero_ncells(2,in)
       write(fnum,'(a,1x,2(1x,i10))') '#INTERPOLATE',ero_bg_interp_opt
       write(fnum,'(a,1x,2(1x,i10))') '#EXTRAPOLATE',ero_bg_extrap_opt
       write(fnum,'(a,1x,2(1x,i10))') '#REMAP_BFIELD',1

       write(fnum,'(a,1x,a)') '#EROBGFILE',trim(ero_bg_output_fn)
       write(fnum,'(a,1x,a)') '#EROSHAPEFILE',trim(ero_shape_output_fn)
       write(fnum,'(a,(1x,i10))') '#EROSHAPEOPT',ero_shape_opt

       write(fnum,'(a,(1x,i10))') '#DIVPARTOPT',ero_part_output_opt
       write(fnum,'(a,(1x,i10))') '#ERODIVPARTMULT',erodiv_tor_part_mult
       write(fnum,'(a,(1x,a))') '#DIVPARTFILE',trim(ero_particle_output_fn)
       write(fnum,'(a,(1x,a))') '#EROPARTFILE',trim(ero_particle_transfer_fn)


    end do

    close(fnum)


  end subroutine write_ero_bg


  subroutine run_ero_bg
    implicit none

    ! this routine runs the code to generate the ERO background plasma files from inside divimp


    character*512 :: ero_bg_command
    character*512 :: oedge_plasma_code = 'ero_plasma '


    ero_bg_command = './'//trim(oedge_plasma_code)//' '// trim(ero_bg_command_fn)
    write(0,'(a)') 'ERO_BG_CMD:',trim(ero_bg_command)


    call system(trim(ero_bg_command))

  end subroutine run_ero_bg


  subroutine ero_cleanup
    implicit none

    ! Clean up ERO related stuff after the DIVIMP run. This includes making sure required data is written and writing the background
    ! plasma files 

    write(0,*) 'ERO_CLEANUP:', ero_export_bg

    if (ero_export_bg.eq.1) then 
       ! write files related to ero background plasma file creation
       call write_ero_bg

       ! run the external code to generate the ERO plasma files. 
       call run_ero_bg

    endif

    ! de-allocate any allocated storage

    if (allocated(ero_coords)) deallocate(ero_coords)
    if (allocated(ero_offsets)) deallocate(ero_offsets) 
    if (allocated(ero_resolution)) deallocate(ero_resolution)
    if (allocated(ero_ncells)) deallocate(ero_ncells)
    if (allocated(ero_tor)) deallocate(ero_tor)

    if (allocated(ero_part_data)) deallocate(ero_part_data)
    if (allocated(ero_launch_prob)) deallocate(ero_launch_prob)


  end subroutine ero_cleanup



  subroutine read_ero_unstructured_input(line,tag,fp)
    use mod_io
    implicit none

    ! use *K for ERO unstructured input identifier? 

    INTEGER   fp
    CHARACTER line*72,tag*3
    integer :: in,ierr

    ierr = 0

    !
    !=======================================================================================
    !
    if (tag(1:3).eq.'K01') then 
       !
       ! K01 - ERO interface code activation option 
       !       0 = off
       !       1 = on
       !

       CALL ReadI(line,ero_opt,0,1,'ERO OPTION ON/OFF')

    elseif (tag(1:3).eq.'K02') then 
       !
       ! ERO background options
       !
       ! K02 - ERO export bg 
       !       0 = off
       !       1 = on
       CALL ReadI(line,ero_export_bg,0,1,'EXPORT BG FOR ERO')

    elseif (tag(1:3).eq.'K03') then 
       ! K03 - ERO bg interpolation option 
       !       0 = off
       !       1+= interpolation option
       CALL ReadI(line,ero_bg_interp_opt,0,1,'BG INTERPOLATION OPT FOR ERO')

    elseif (tag(1:3).eq.'K04') then 
       ! K04 - ERO bg extrapolation option 
       !       0 = off
       !       1+= extrapolation option 
       !        
       !       0= extrapolated values are set to 0.0
       !       1= set extrapolated value to value from nearest surface 
       CALL ReadI(line,ero_bg_extrap_opt,0,1,'BG EXTRAPOLATION OPT FOR ERO')

    elseif (tag(1:3).eq.'K05') then 
       ! K05 - ERO bg instruction filename
       !

       CALL ReadC(line,ero_bg_command_fn,'NAME FOR ERO BG PLAMSA COMMAND FILE')


    elseif (tag(1:3).eq.'K06') then 

       ! K06 - ERO bg output filename base
       !
       ! n will be replaced by region number

       CALL ReadC(line,ero_bg_output_fn_base,'BASE NAME FOR ERO BG PLAMSA FILES')


    elseif (tag(1:3).eq.'K07') then 

       ! K07 - ERO region specifications 
       !     - number of regions
       !     - coordinates 
       !     - offsets
       !     - resolution - either dx,dy or nx,ny
       !
       !
       ! Allocate storage
       !
       ! Read in blocks of ERO data

       CALL ReadI(line,n_ero_regs,0,9999,'NUMBER OF ERO REGIONS')

       if (n_ero_regs.gt.1) then
          call errmsg('ERROR: Input K07: Only one ERO region supported at a time',n_ero_regs)
       endif

       call allocate_ero_region_data

       do in = 1,n_ero_regs

          call read_ero_reg_input(in,ierr)
          if (ierr.ne.0) exit

       end do

       ! assign ero_rvert and ero_zvert from the first ero region read
       do in = 1,4
          ero_rvert(in) = ero_coords(1,in,1)
          ero_zvert(in) = ero_coords(2,in,1)
       end do 



    elseif (tag(1:3).eq.'K08') then 

       ! K08 - ERO shape output filename base
       !
       ! n will be replaced by region number

       CALL ReadC(line,ero_shape_output_fn_base,'BASE NAME FOR ERO SURFACE FILES')


    elseif (tag(1:3).eq.'K09') then 

       ! K06 - ERO shape output filename base
       !
       ! n will be replaced by region number

       CALL ReadI(line,ero_shape_opt,0,2,'REQUEST ERO SHAPE GEOMETRY OUTPUT')


       !
       !=======================================================================================
       !

    elseif (tag(1:3).eq.'K20') then 
       ! Generate a file showing particles entering ERO simulation volume
       ! Needs to check neutrals and ions - LAUNCH and DIV
       CALL ReadI(line,ero_part_output_opt,0,1,'COLLECT PARTICLE DATA FOR ERO LAUNCH')


    elseif (tag(1:3).eq.'K21') then 

       ! K21 - DIVIMP particle track output file for ERO
       !
       ! n will be replaced by region number

       CALL ReadC(line,ero_particle_output_fn,'NAME FOR DIVIMP PARTICLE OUTPUT')

    elseif (tag(1:3).eq.'K22') then 

       ! K22 - DIVIMP particle track output file for direct input to ERO - coordinates and units transformed as required
       !
       ! n will be replaced by region number

       CALL ReadC(line,ero_particle_transfer_fn,'NAME FOR ERO PARTICLE INPUT FN')


    elseif (tag(1:3).eq.'K23') then 
       !
       ! DIVIMP records particles at a single toroidal location which is not a good representation
       ! of a toroidally symmetric particle source. To help fix this each DIVIMP particle will be
       ! split into tor_mult_part pieces in ERO and distributed toroidally. Initially this distribution
       ! will be random but other options are possible. 
       !
       
       CALL ReadI(line,erodiv_tor_part_mult,0,100,'NUMBER of ERO particles for each DIVIMP particle')

       !
       !=======================================================================================
       !

    elseif (tag(1:3).eq.'K40') then 
       ! 
       ! K40 - Option to launch particle source loaded from ERO
       !
       
       CALL ReadI(line,ero_particle_launch_opt,0,1,'LAUNCH ERO PARTICLES')

    elseif (tag(1:3).eq.'K41') then 

       ! K08 - ERO shape output filename base
       !
       ! n will be replaced by region number

       CALL ReadC(line,ero_particle_input_fn,'NAME FOR ERO PARTICLE INPUT FILE')


       !
       !=======================================================================================
       !
    elseif (tag(1:3).eq.'K42') then 
       ! 
       ! K40 - Option to launch particle source loaded from ERO
       !
       
       CALL ReadI(line,ero_remove_src_opt,0,3,'Option to turn off DIVIMP surface sources in ERO volume')

    endif



  end subroutine read_ero_unstructured_input




  subroutine allocate_ero_region_data
    implicit none
    ! when number of regions is known the data storage for the region data can be allocated
    integer :: ierr


    if (n_ero_regs.gt.0) then 

       if (allocated(ero_coords)) deallocate(ero_coords)
       allocate(ero_coords(2,4,n_ero_regs),stat=ierr)
       ero_coords = 0.0


       if (allocated(ero_offsets)) deallocate(ero_offsets)
       allocate(ero_offsets(2,n_ero_regs),stat=ierr)
       ero_offsets = 0.0


       if (allocated(ero_resolution)) deallocate(ero_resolution)
       allocate(ero_resolution(2,n_ero_regs),stat=ierr)
       ero_resolution = 0.0


       if (allocated(ero_ncells)) deallocate(ero_ncells)
       allocate(ero_ncells(2,n_ero_regs),stat=ierr)
       ero_ncells(1,:) = ero_default_nx
       ero_ncells(2,:) = ero_default_ny


       if (allocated(ero_tor)) deallocate(ero_tor)
       allocate(ero_tor(n_ero_regs),stat=ierr)
       ! set default to +/- 10cm - can be changed
       ero_tor= 0.1

    endif


  end subroutine allocate_ero_region_data




  subroutine read_ero_reg_input(in,ierr)
    use mod_io
    implicit none

    ! Read a block of ero_reg input - up to the 'END BL' tag
    ! If 'END BL' is not found or another error encountered then exit/stop

    integer :: in
    logical :: finished
    ! use module input buffer called line
    !character*512 :: line
    integer :: ierr
    real*8 :: xlen,ylen

    finished = .false.

    ierr =0 

    do while (.not.finished)

       call rdbuffer(line,'READ BUFFER',ierr)
       
       if (ierr.eq.0) then 


          if (line(2:7).eq.'P1, P2') then 
             ! read points 1,2
             read(line(10:),*) ero_coords(1,1,in),ero_coords(2,1,in),ero_coords(1,2,in),ero_coords(2,2,in)
          elseif (line(2:7).eq.'P4, P3') then 
             ! read points 4,3
             read(line(10:),*) ero_coords(1,4,in),ero_coords(2,4,in),ero_coords(1,3,in),ero_coords(2,3,in)
          elseif (line(2:7).eq.'OFFSET') then 
             ! read offset for output
             read(line(10:),*) ero_offsets(1,in),ero_offsets(2,in)
          elseif (line(2:7).eq.'DX, DY') then 
             ! read grid cell size
             ! takes precedence over grid resolution
             read(line(10:),*) ero_resolution(1,in),ero_resolution(2,in)
          elseif (line(2:7).eq.'NX, NY') then 
             ! read grid number of cells
             ! issue warning if DX, DY are not zero 
             if (ero_resolution(1,in).ne.0.0.and.ero_resolution(2,in).ne.0.0) then
                call errmsg('READ_ERO_REG_INPUT:','NON-ZERO VALUES OF DX,DY HAVE BEEN SPECIFIED - OVER-RIDING NX,NY INPUT')
             endif
             read(line(10:),*) ero_ncells(1,in),ero_ncells(2,in)
          elseif (line(2:11).eq.'TOR_EXTENT') then   
             ! Toridal extent (+/-) in meters of ERO region
             read(line(14:),*) ero_tor(in)
          elseif (line(2:7).eq.'END BL') then 
             ! End of ERO Block - exit 
             finished = .true.
          endif

       else
          ! Error condition - stop the code and issue an error message
          call errmsg('READ_ERO_REG_INPUT: ERROR READING BUFFER:',ierr)
          finished = .true. 
          stop 'ERROR: READ ERO INPUT: EXIT'
       endif

    end do

    if (ierr.eq.0) then 
       ! calculate the resolution and number of cells
       ! The ERO space is assumed rectangular with the X axis between 1,2 and the Y axis 1,4
       xlen = sqrt( (ero_coords(1,2,in)-ero_coords(1,1,in))**2 + (ero_coords(2,2,in)-ero_coords(2,1,in))**2 )
       ylen = sqrt( (ero_coords(1,4,in)-ero_coords(1,1,in))**2 + (ero_coords(2,4,in)-ero_coords(2,1,in))**2 )

       if (ero_resolution(1,in).ne.0.0.and.ero_resolution(2,in).ne.0.0) then 
          ! grid cell resolution specified in spatial coordinates - rounded off to nearest cell size
          ero_ncells(1,in) = int(xlen/ero_resolution(1,in))+1 
          ero_ncells(2,in) = int(ylen/ero_resolution(2,in))+1 

          ! reset resoution to match actual number of cells
          ero_resolution(1,in) = xlen/ero_ncells(1,in)
          ero_resolution(2,in) = ylen/ero_ncells(2,in)

       else
          ! Either the number of cells has been specified or the default values are used

          ! set resoution to match number of cells
          ero_resolution(1,in) = xlen/ero_ncells(1,in)
          ero_resolution(2,in) = ylen/ero_ncells(2,in)

       endif

    endif



  end subroutine read_ero_reg_input


  subroutine init_ero_unstructured_input
    implicit none



    !
    !=======================================================================================
    !
    ! Set default values for K series ERO-DIVIMP interface options
    !
    ! K01 - ERO interface code activation option 
    !       0 = off
    !       1 = on
    !
    ero_opt = 0


    !------------------------- 
    ! ERO background options
    !
    ! K02 - ERO export bg 
    !       0 = off
    !       1 = on
    !
    ero_export_bg = 0

    ! K03 - ERO bg interpolation option 
    !       0 = off
    !       1+= interpolation option
    !
    ero_bg_interp_opt = 0

    ! K04 - ERO bg interpolation option 
    !       0 = off
    !       1+= extrapolation option 
    !        
    !       0= extrapolated values are set to 0.0
    !       1= set extrapolated value to value from nearest surface 
    !
    ero_bg_extrap_opt = 1

    ! K05 - ERO bg instruction filename
    !
    !
    ero_bg_command_fn = 'ero_bg_cmd.dat'

    ! K06 - ERO bg output filename base
    !
    ! n will be replaced by region number
    !
    ! _n.dat
    !
    ero_bg_output_fn_base = 'erodiv_plasma'

    !
    ! K07 - ERO region specifications
    !
    n_ero_regs = 0

    !
    !

    ! K08 - ERO shape output filename base
    !
    ! n will be replaced by region number
    !
    ! _n.dat
    !
    ero_shape_output_fn_base = 'erodiv_geometry'


    ! K09 - ERO shape geometry output option
    !
    ! options : 0,1,2 ... 0 off, 1 cylindrical, 2 toroidal
    !

    ero_shape_opt = 0



    !
    !=======================================================================================
    !
    !--------------------------------------
    !
    ! ERO particle data recording options
    !
    ! Start series at K20 for now
    !
    !
    ! K20 - record particle data for passing to ERO
    !       0 = off
    !       1 = on
    !
    ero_part_output_opt = 0

    !
    ! K21 - file name to store DIVIMP recorded particle data
    !       n will be region number
    !
    ero_particle_output_fn = 'divimp_particle_data.dat'


    !
    ! K22 - file name to store DIVIMP recorded particle data
    !       n will be region number
    !
    ero_particle_transfer_fn = 'divero_particle_data.dat'


    ! K23 - toroidal particle multiplier
    ! DIVIMP records particles at a single toroidal location which is not a good representation
    ! of a toroidally symmetric particle source. To help fix this each DIVIMP particle will be
    ! split into tor_mult_part pieces in ERO and distributed toroidally. Initially this distribution
    ! will be random but other options are possible. 
    !
    erodiv_tor_part_mult = 10

    !
    !=======================================================================================
    !
    !--------------------------------------
    !
    ! ERO particle data launch options
    !
    ! Start series at K40 for now
    !
    ! 
    ! K40 - launch DIVIMP particles based on ERO data
    !       0 = off
    !       1+= launch option (if different options possible)
    !
    ero_particle_launch_opt = 0

    !
    ! K41 - file name to load ERO recorded particle data
    !       n will be region number (if more than one specified). 
    !
    ero_particle_input_fn = 'ero_particle_data.dat'

    !
    ! K42 - option to turn off DIVIMP sources within the ERO source region
    !       0 = off
    !       1 = ON - if wall segment center is in the region then source is turned off
    !       2 = ON - if both wall segment end points are in region then the source is turned off
    !       3 = ON - if any end point is in region then the source is turned off

    ero_remove_src_opt = 0

  end subroutine init_ero_unstructured_input



  subroutine load_ero_launch_data
    implicit none

    integer :: inunit,ierr

    real*8 :: r,z,t,vr,vz,vt,vtot,temp,charge,atoms

    integer :: data_count, in

    integer :: headers_read 

    write(0,'(a,a)') 'Loading ERO particle data:',trim(ero_particle_input_fn)


    call find_free_unit_number(inunit)

    open(unit=inunit,file=trim(ero_particle_input_fn),status='old',iostat=ierr)
    if (ierr.ne.0) then 
       call errmsg('LOAD_ERO_LAUNCH_DATA:ERROR OPENING :'//trim(ero_particle_input_fn),ierr)
       stop 'LOAD_ERO_LAUNCH_DATA'
    endif

    ! Interface has done most of the processing for DIVIMP particle launch

    headers_read = 0

    do while(headers_read<3) 

       read(inunit,line_form,iostat=ierr) line

       if (ierr.eq.0) then 
          
          write(0,'(a,a,a)') 'LINE:',trim(line)

          if (line(1:11) .eq.'#NPARTICLES') then 
             headers_read = headers_read + 1
             read(line(12:),*) data_count

          elseif (line(1:12).eq.'#TOTAL_ATOMS') then 
             headers_read = headers_read + 1
             read(line(13:),*) ero_totsrc_raw, ero_totsrc_mtor
             ! set the ero total particle source to the integrated particles/s/m-tor value
             ero_totsrc=ero_totsrc_mtor

          elseif (line(1:16 ).eq.'#FULL_TOR_EXTENT') then 
             headers_read = headers_read + 1
             read(line(17:),*) ero_src_tor_extent

          endif
       endif

    end do

    if (headers_read.ne.3) then 
       call errmsg('ERROR: Load ERO particle data:','Problem reading headers!')
       stop 'load_analyse_ero_particle_data'
    endif

    ! start at beginning of input file just for completeness
    rewind(inunit)

    ! allocate storage for the distribution arrays

    allocate(ero_part_data(data_count,n_ero_data+n_extra_data),stat=ierr)
    allocate(ero_launch_prob(data_count),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('ERROR: load_analyse_ero_launch_data : allocating ero_part_data:',ierr)
       stop 'load_analyse_ero_launch_data: allocation'
    endif


    n_ero_part = data_count

    data_count =0

    ! generate a list of locations and relative launch probabilities along with launch data including initial charge state

    do while (ierr.eq.0) 

       read(inunit,line_form,iostat=ierr) line


       if (ierr.eq.0) then 
          ! ignore comments
          if (line(1:1).eq.'#') cycle

          data_count = data_count+1
          read(line,*) (ero_part_data(data_count,in),in=1,n_ero_data),ero_launch_prob(data_count)

       endif

    end do

    !write(6,*) 'ERO prob dist:',n_ero_part
    !do in = 1,n_ero_part
    !   write(6,'(a,i6,4(1x,g18.8))') 'ERO_PROB:',in,ero_part_data(in,3),ero_part_data(in,4),ero_part_data(in,12),ero_launch_prob(in)
    !end do

  end subroutine load_ero_launch_data


  subroutine select_ero_particle(r,z,in,rand)
    implicit none
    real :: rand
    real :: r,z,t,vr,vz,vt,temp,charge

    integer, external :: ipos
    integer :: in


    ! transfer 13 data values : in, charge, r, z, t, vr, vz, vt, v0, energy, atoms, probabilty, cumulative probability

    in = ipos(rand,ero_launch_prob,n_ero_part)

    
    !pnum= int(ero_part_data(in,1))
    !charge  = int(ero_part_data(in,2))

    r = ero_part_data(in,3)
    z = ero_part_data(in,4)

    !write(0,'(a,2g18.8,i7,g18.8)') 'Select particle:',r,z,in,rand


    !t = ero_part_data(in,5)
    !vr = ero_part_data(in,6)
    !vz = ero_part_data(in,7)
    !vt = ero_part_data(in,8)
    !vtotal = ero_part_data(in,9)
    !temp = ero_part_data(in,10)
    !atoms = ero_part_data(in,11)
    !prob = ero_part_data(in,12)

    ! increment bin count when selected
    ero_part_data(in,n_ero_data+1) = ero_part_data(in,n_ero_data+1) + 1.0

    return

  end subroutine select_ero_particle

  subroutine get_ero_part_data(in,vr,vz,vt,vtotal,temp,charge)
    implicit none
    ! load ERO particle data based on index
    real :: vr,vz,vt,temp,vtotal
    !real :: r,z,t,vr,vz,vt,temp,charge
    real :: charge
    integer :: in

    !pnum= int(ero_part_data(in,1))
    charge  = int(ero_part_data(in,2))

    !r = ero_part_data(in,3)
    !z = ero_part_data(in,4)

    !t = ero_part_data(in,5)
    vr = ero_part_data(in,6)
    vz = ero_part_data(in,7)
    vt = ero_part_data(in,8)
    vtotal = ero_part_data(in,9)
    temp = ero_part_data(in,10)
    !atoms = ero_part_data(in,11)
    !prob = ero_part_data(in,12)

    return

  end subroutine get_ero_part_data

  subroutine get_ero_part_data2(in,r,z,temp,charge,prob)
    implicit none
    ! load ERO particle data based on index
    real :: r,z,temp,prob
    !real :: r,z,t,vr,vz,vt,temp,charge,prob
    integer :: in
    real :: charge

    !pnum= int(ero_part_data(in,1))
    charge  = int(ero_part_data(in,2))

    r = ero_part_data(in,3)
    z = ero_part_data(in,4)

    !t = ero_part_data(in,5)
    !vr = ero_part_data(in,6)
    !vz = ero_part_data(in,7)
    !vt = ero_part_data(in,8)
    !vtotal = ero_part_data(in,9)
    temp = ero_part_data(in,10)
    !atoms = ero_part_data(in,11)
    prob = ero_part_data(in,12)

    return

  end subroutine get_ero_part_data2



  subroutine get_ero_src(totsrc,totnum)
    implicit none
    real :: totsrc
    integer :: totnum
    
    ! DIVIMP expects the ERO source strength in particles/m-tor/s
    ! ERO returns values in particles/s since each test particle in 
    ! ERO has a variable weight depending on its source.
    ! Scaling this is a challenge. 
    ! 
    write(0,*) 'Getting total source:',n_ero_part,ero_totsrc
    
    totsrc = ero_totsrc
    totnum = n_ero_part

  end subroutine get_ero_src


  subroutine init_ero_part_output
    implicit none
    ! open the output file 
    integer :: ierr

    !call find_free_unit_number(ero_part_output_unit)

    open(ero_part_output_unit,file=trim(ero_particle_output_fn),iostat=ierr)

    write(ero_part_output_unit,'(a)') '# DIVIMP particles entering ERO simulation volume:'
    write(ero_part_output_unit,'(a)') '# R   Z     CHARGE     VR    VZ     VT    VPARA    TEMP   WEIGHT '

    write(0,'(a,1x,i8,1x,a)') 'Initializing DIVIMP->ERO particle output:',ero_part_output_unit,trim(ero_particle_output_fn)

  end subroutine init_ero_part_output


  subroutine close_ero_part_output(absfac,nparticles)
    implicit none
    real :: absfac
    real :: nparticles
    ! write absolute factor and scaling data including weight or atoms/particle
    ! HOW to figure out atoms/particle??
    ! Scaling should be particles/m-tor/s
    ! Need to multiply by the toroidal extent of the ERO simulation to get 
    ! hpw many particles would be entering along that edge of the ERO simulation
    !
    ! If cylindrical option specified use the tor_extent as is.
    ! Tor_length = 2 * tor_extent
    !
    ! If toroidal geometry 
    ! Rprime = sqrt(Rp**2 + Tor_extent**2)
    ! Dist = 2 PI Rprime
    ! angle = 2 x tan(tor_extent/Rp)  (angle in radians)
    ! Tor_length = angle * Rprime          (angle/2PI  * 2PI Rprime)
    !

    !if (ero_shape_opt.eq.0.or.ero_shape_opt.eq.1) then ! cylindrical

    ! for now just use the cylindrical approximation since it won't be off by much
    ! Note: assumes only one ERO block
    !
    write(ero_part_output_unit,'(a,3(1x,g18.8))')  '#ABSFAC',absfac,nparticles,absfac/nparticles*(2.0*ero_tor(1))

    !elseif (ero_shape_opt.eq.2) then ! toroidal
    !   write(ero_part_output_unit,'(a,3(1x,18.8))')  '#ABSFAC:',absfac,real(nparticles),absfac/real(nparticles)*(2,0*tor_extent)

    write(0,*) 'Closing DIVIMP->ERO particle output:',ero_part_output_unit


    close(ero_part_output_unit)
    


  end subroutine close_ero_part_output


  subroutine ero_init_part_track(r,z,ero_record_data)
    implicit none
    real :: r,z
    logical :: ero_record_data
    logical,external :: inpoly

    ! Check to see if particle starts off inside ERO region

    ero_last_in_vol = inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

    ero_in_vol = ero_last_in_vol

    ero_record_data = .true.

    !write(0,'(a,2l6)') 'ERO_START_PARTICLE:',ero_in_vol,ero_record_data


  end subroutine ero_init_part_track



  subroutine ero_check_part_track_ion(ik,ir,iz,r,z,vel,qtim,temi,sputy,ero_record_data)
    use bfield
    implicit none
    real :: r,z,vel,temi,sputy,qtim
    integer :: ik,ir,iz
    logical :: ero_record_data
    logical,external :: inpoly
    real :: b(3)

    ! Check if particle has entered the ERO sample region - if so record and mark it as not to follow any further
    ! Record means to write the particle information to the file
    !
    ! Record current R,Z for now - this could be made more fancy by recording the actual intersection with the boundary but
    ! actually being on the boundary might cause other issues in ERO. 
    !
    
    ero_in_vol = inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

    if (ero_last_in_vol) then 
       if (.not.ero_in_vol) then 
          ! if particle has exited ERO volume set the last flag to false
          ero_last_in_vol = .false. 
       endif
    else
       if (ero_in_vol) then 
          ! record data and set ero_record_data=.false.so that the particle will not be checked again. 
          ! Particle can only enter recording volume once from the outside
          call get_bvector(ik,ir,b)
          
          !write(0,'(a,2(1x,g18.8),1x,i8,1x,6(1x,g18.8))') 'ERO PART ION:',r,z,iz,vel*b(1),vel*b(2),vel*b(3),vel,temi,sputy

          !write(ero_part_output_unit,'(a)') 'TEST!'

          write(ero_part_output_unit,'(2(1x,g18.8),1x,i8,1x,6(1x,g18.8))') r,z,iz,vel*b(1)/qtim,vel*b(2)/qtim,vel*b(3)/qtim,vel/qtim,temi,sputy

          write(0,'(a,i6,2(1x,g18.8),1x,i8,1x,6(1x,g18.8))') 'ERO PART ION:',ero_part_output_unit,r,z,iz,vel/qtim,b(1),b(2),b(3),sqrt(b(1)**2+b(2)**2+b(3)**2)
          ero_record_data = .false. 

       endif

    endif


  end subroutine ero_check_part_track_ion


  subroutine ero_check_part_track_neut(ik,ir,iz,r,z,temn,sputy,vx,vy,fsrate,crmi,ero_record_data)
    implicit none
    real :: r,z,temn,sputy,vx,vy,fsrate,vel,crmi
    integer :: ik,ir,iz
    logical :: ero_record_data
    logical,external :: inpoly
    real,external :: getranf
    real :: vr,vz,vt

    ! Check if particle has entered the ERO sample region - if so record and mark it as not to follow any further
    ! Record means to write the particle information to the file
    !
    ! Record current R,Z for now - this could be made more fancy by recording the actual intersection with the boundary but
    ! actually being on the boundary might cause other issues in ERO. 
    !
    
    ero_in_vol = inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

    if (ero_last_in_vol) then 
       if (.not.ero_in_vol) then 
          ! if particle has exited ERO volume set the last flag to false
          ero_last_in_vol = .false. 
       endif
    else
       if (ero_in_vol) then 
          ! record data and set ero_record_data=.false.so that the particle will not be checked again. 
          ! Particle can only enter recording volume once from the outside
          
          Vel = 1.38E4 * SQRT (temn/CRMI)

          vr = vx/fsrate
          vz = vy/fsrate

          vt = vel**2 - vr**2 - vz**2

          if (vt.gt.0.0) then 
             vt = (2.0 * getranf() -1.0 ) * sqrt (vt)
          else
             vt = 0.0
             vel = sqrt(vr**2 + vz**2) 
          endif

          !write(0,'(a,2(1x,g18.8),1x,i8,1x,6(1x,g18.8))') 'ERO_NEUT:',r,z,iz,vr,vz,vt,vel,temn,sputy
          write(ero_part_output_unit,'(2(1x,g18.8),1x,i8,1x,6(1x,g18.8))') r,z,iz,vr,vz,vt,vel,temn,sputy

          write(0,'(a,2(1x,g18.8),1x,i8,1x,6(1x,g18.8))') 'ERO PART NEUT:',r,z,iz,vr,vz,vt,vel

          ero_record_data = .false. 

       endif

    endif

  end subroutine ero_check_part_track_neut


  subroutine ero_remove_src(wallpts,wallpt,kmfps,kmfcs,kmfss,kmfpws,kmfcws,maxpts,nwall_data,maxnds)
    integer :: wallpts, maxpts, nwall_data, maxnds

    real :: wallpt(maxpts,nwall_data)
    real :: kmfps(maxnds),kmfcs(maxnds),kmfss(maxnds)
    real :: kmfpws(maxpts),kmfcws(maxpts)

    integer :: in,id
    logical :: remove_source
    logical :: testa,testb

    real :: r,z

    logical,external :: inpoly

    ! do nothing if the option is off
    if (ero_remove_src_opt.eq.0) return


    ! Loop over wall and zero out yields where appropriate depending on options
    ! note the target index (if any) for a given wall segment is in wallpt(in,18) 

    do in = 1,wallpts

       if (ero_remove_src_opt.eq.1) then 
          ! if wall center point is in ERO region then zero yields
          r = wallpt(in,1)
          z = wallpt(in,2)

          remove_source= inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

       elseif (ero_remove_src_opt.eq.2) then 
          ! if both end points of the wall segment are in the ERO region then zero yields
          ! start point
          r=wallpt(in,20)
          z=wallpt(in,21)

          testa= inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

          r=wallpt(in,22)
          z=wallpt(in,23)

          testb= inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

          remove_source = (testa.and.testb)

       elseif (ero_remove_src_opt.eq.3) then 
          ! if either end point of the wall segment are in the ERO region then zero yields

          r=wallpt(in,20)
          z=wallpt(in,21)

          testa= inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

          r=wallpt(in,22)
          z=wallpt(in,23)

          testb= inpoly(r,z,ero_nvert,ero_rvert,ero_zvert)

          remove_source = (testa.or.testb)

       elseif (ero_remove_src_opt.eq.4) then 
          ! ** Not yet implemented
          ! This option could set the yield for a segment based on the fraction
          ! of the segment length inside the ERO region. 
          ! However, this would not prevent particles from still being launched inside the ERO 
          ! region ... just scale the source appropriately for wall segments near the edge
          ! of the ERO region. 
          ! The best solution is probably a well resolved grid/wall in the ERO region so that 
          ! edge effects are minimized

       endif

       ! reset the yields related to this wall segment to zero
       if (remove_source) then 

          id = int(wallpt(in,18))
          
          kmfps(id) = 0.0
          kmfcs(id) = 0.0
          kmfss(id) = 0.0

          kmfpws(in) = 0.0
          kmfcws(in) = 0.0

          write(0,'(a,2(1x,g18.8),2i8)') 'ERO_REMOVE_SRC: Zeroed yields for:',r,z,in,id

       endif
   end do
 end subroutine ero_remove_src

end module ero_interface
  
