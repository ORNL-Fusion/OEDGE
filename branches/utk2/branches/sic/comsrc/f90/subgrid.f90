module subgrid

  use error_handling  ! Error message handling routines
  use subgrid_options ! Option data applying to the subgrid

  !
  ! The purpose of this module is to provide configurable support for local grid
  ! resolution enhancement - primarily for recording local densities on a more
  ! finely defined grid. This could be done by overlaying a secondary grid either
  ! in RZ space or by subdividing cells on the grid. Each approach has some 
  ! advantages. Grid refinement makes it easy to match subgrids to plasma conditions
  ! which are recorded on the master grid. Regular underlying grids make it much easier
  ! to calculate cell volumes and cell indices when incrementing density. By using a 
  ! module and allocatable storage it should be possible to support both approaches if 
  ! that is desired. 
  !
  ! Note: The subgrid is an R,Z (or X,Y) overlay on the standard magnetic modeling grid. As
  !       such the locations of the ions need to be mapped to R,Z when their location is being recorded. 
  !       The GETRZ routine should be able to handle most of this requirement.
  !
  !

  implicit none

  private
  save


  !integer :: subgrid_opt
  !integer :: sg_rdim,sg_zdim         ! Number of cells in the R and Z directions and number of different states
  !real :: sg_rmin,sg_rmax,sg_zmin,sg_zmax  ! Window boundaries of the subgrid region
  real :: dr,dz                ! Width of the R and Z cells    


  integer :: sg_hc_opt

  real,allocatable :: subgrid_data(:,:,:)    ! 3D array - R,Z, state
  real,allocatable :: hc_subgrid_data(:,:,:)    ! 3D array - R,Z, state
  !real,allocatable :: raxis(:),zaxis(:)

  integer :: max_state
  integer :: max_hc_state

  !
  ! Set min state index to 0 for now in init_subgrid - can be changed later 
  !
  integer :: min_state
  integer :: min_hc_state

  public :: init_subgrid,update_subgrid,hc_update_subgrid,save_subgrid,norm_subgrid,reload_subgrid
  public :: load_subgrid_xsection,calc_subgrid_axis

contains

  subroutine init_subgrid(global_hc_opt,nstates,hc_states)
    implicit none

    integer :: global_hc_opt
    integer :: nstates,hc_states  ! Number of cells in the R and Z directions

    !
    ! Initializes the subgrid data - allocates arrays - sets internal variables and sets up any extra information required
    ! for locating the cell of the subgrid occupied by the particle.
    !
    real :: rrange, zrange

    !write(0,*) 'Init_subgrid:',global_hc_opt,nstates,hc_states

    !
    ! Assign initial values of min_state and min_hc_state
    !
    min_state = 0
    min_hc_state = 0

    !
    ! Exit initialization if subgrid is off 
    !
    if (subgrid_opt.eq.0) return

    !
    ! Indicates if the HC code is active - in which case additional storage is added to accomodate the HC data 
    !
    sg_hc_opt = global_hc_opt

    max_state = nstates
    max_hc_state = hc_states

    rrange= sg_rmax-sg_rmin
    zrange= sg_zmax-sg_zmin

    dr = rrange/sg_rdim
    dz = zrange/sg_zdim

    call allocate_subgrid

    if (sg_hc_opt.ne.0) then 
       call allocate_hc_subgrid
    endif

  end subroutine init_subgrid


  subroutine get_subgrid_indices(r,z,r_index,z_index)
    implicit none
    real :: r,z
    integer :: r_index,z_index

    !
    ! Use -1 index coordinates to indicate that the particle is not in the subgrid region
    !
    if (r.lt.sg_rmin) then
       r_index = -1
    elseif (r.ge.sg_rmax) then
       r_index = -1
    else
       r_index = int((r-sg_rmin)/dr) + 1
    endif

    if (z.lt.sg_zmin) then
       z_index = -1
    elseif (z.ge.sg_zmax) then
       z_index = -1
    else
       z_index = int((z-sg_zmin)/dz) + 1
    endif

  end subroutine get_subgrid_indices


  subroutine allocate_subgrid
    implicit none
    integer :: iflag

    allocate(subgrid_data(sg_rdim,sg_zdim,min_state:max_state),stat=iflag)

    if (iflag.ne.0) then 
       call errmsg('ALLOCATE_SUBGRID: ALLOCATION ERROR',iflag)
       subgrid_opt = 0
       return
    endif

    subgrid_data = 0.0

  end subroutine allocate_subgrid

  subroutine allocate_hc_subgrid
    implicit none
    integer :: iflag

    allocate(hc_subgrid_data(sg_rdim,sg_zdim,min_hc_state:max_hc_state),stat=iflag)

    if (iflag.ne.0) then 
       call errmsg('ALLOCATE_HC_SUBGRID: ALLOCATION ERROR',iflag)
       sg_hc_opt = 0
       return
    endif

    hc_subgrid_data = 0.0

  end subroutine allocate_hc_subgrid


  subroutine deallocate_subgrids
    implicit none
! slmod begin - tmp
!... Unable to do a proper sync at the moment, so need to add 
!    some checks here, since SUBGRID data is not being allocated 
!    in the first place:

    if (allocated(subgrid_data)) deallocate(subgrid_data)

    if (sg_hc_opt.ne.0) then
       if (allocated(hc_subgrid_data)) deallocate(hc_subgrid_data)
    endif

! slmod end
  end subroutine deallocate_subgrids



  subroutine update_subgrid(r,z,state,sputy)
    implicit none
    real :: r,z
    real :: sputy
    integer :: state

    !
    ! Local variables
    !
    integer :: r_index,z_index

    !
    ! Only record in range states
    !
    if (state.lt.min_state.or.state.gt.max_state) return


    call get_subgrid_indices(r,z,r_index,z_index)

    !
    ! Check to see if particle is on the subgrid.
    !
    if (r_index.lt.0.or.z_index.lt.0) return


    !
    ! Update subgrid data for particles on the grid
    !

    subgrid_data(r_index,z_index,state) = subgrid_data(r_index,z_index,state) + sputy


  end subroutine update_subgrid

  subroutine hc_update_subgrid(r,z,state,sputy)
    implicit none
    real :: r,z
    real :: sputy
    integer :: state

    !
    ! Local variables
    !
    integer :: r_index,z_index


    !
    ! Only record in range states
    !
    if (state.lt.min_hc_state.or.state.gt.max_hc_state) return

    call get_subgrid_indices(r,z,r_index,z_index)

    !
    ! Check to see if particle is on the subgrid.
    !
    if (r_index.lt.0.or.z_index.lt.0) return

    !
    ! Update subgrid data for particles on the grid
    !

    hc_subgrid_data(r_index,z_index,state) = hc_subgrid_data(r_index,z_index,state) + sputy


  end subroutine hc_update_subgrid


  subroutine save_subgrid(ounit)
    implicit none
    integer,intent(in) :: ounit
    !
    ! This routine writes the subgrid data to the RAW DIVIMP output file 
    !
    integer :: ib1,ib2,ib3
    integer :: ios
    !
    ! Always write the subgrid_opt value
    !
    !write(0,*) 'OUNIT=',ounit

    write(ounit,err=100,iostat=ios) subgrid_opt


    !
    ! If the subgrid option is on - write the relevant data - and deallocate arrays afterwards
    !

    if (subgrid_opt.ne.0) then 


       !write(0,*) '2:',subgrid_opt,allocated(subgrid_data),min_state,max_state

       write(ounit,err=100,iostat=ios) sg_rdim,sg_zdim,min_state,max_state,sg_rmin,sg_rmax,sg_zmin,sg_zmax
       write(ounit,err=100,iostat=ios) (((subgrid_data(ib1,ib2,ib3),ib1=1,sg_rdim),ib2=1,sg_zdim),ib3=min_state,max_state)

       ! HC subgrid data if available

       write(ounit,err=100,iostat=ios) sg_hc_opt


       !write(0,*) '3:',sg_hc_opt,allocated(hc_subgrid_data),min_hc_state,max_hc_state

       if (sg_hc_opt.ne.0) then 

          write(ounit,err=100,iostat=ios) min_hc_state,max_hc_state
! slmod begin
          ! Krieger IPP/07 - SUN compiler insists on 132 column limit
          write(ounit,err=100,iostat=ios) &
			 (((hc_subgrid_data(ib1,ib2,ib3),ib1=1,sg_rdim),ib2=1,sg_zdim),ib3=min_hc_state,max_hc_state)
!
!          write(ounit,err=100,iostat=ios) (((hc_subgrid_data(ib1,ib2,ib3),ib1=1,sg_rdim),ib2=1,sg_zdim),ib3=min_hc_state,max_hc_state)
! slmod end

       endif
    endif

    !
    ! Deallocate the subgrids
    !

    call deallocate_subgrids

    return

100 call errmsg('SAVE_SUBGRID:','ERROR WRITING SUBGRID DATA TO RAW FILE')

    return



  end subroutine save_subgrid



  subroutine reload_subgrid(iunit)
    implicit none
    integer,intent(in) :: iunit
    !
    ! This routine is used from within the OUT program to reload the subgrid data from the RAW data file. Option data is first read and then 
    ! the arrays appropriately allocated as required - then the data arrays are loaded. 
    !
    integer :: ib1,ib2,ib3
    integer :: ios
    !
    ! Always write the subgrid_opt value
    !
    read(iunit,err=100,iostat=ios) subgrid_opt

    !
    ! If the subgrid option is on - restore the relevant data
    !

    if (subgrid_opt.ne.0) then 

       !
       ! Load basic data
       !

       read(iunit,err=100,iostat=ios) sg_rdim,sg_zdim,min_state,max_state,sg_rmin,sg_rmax,sg_zmin,sg_zmax

       !
       ! Allocate the subgrid_data array - the specifics have been read into module variables in the last line
       !

       call allocate_subgrid

       read(iunit,err=100,iostat=ios) (((subgrid_data(ib1,ib2,ib3),ib1=1,sg_rdim),ib2=1,sg_zdim),ib3=min_state,max_state)

       ! HC subgrid data if available

       read(iunit,err=100,iostat=ios) sg_hc_opt

       if (sg_hc_opt.ne.0) then 

          read(iunit,err=100,iostat=ios) min_hc_state,max_hc_state
          !
          ! Allocate HC subgrid
          !
          call allocate_hc_subgrid

          !
          ! Read in HC subgrid data
          !
! slmod begin       
          ! Krieger IPP/07 - SUN compiler insists on 132 column limit
          read(iunit,err=100,iostat=ios) &
			 (((hc_subgrid_data(ib1,ib2,ib3),ib1=1,sg_rdim),ib2=1,sg_zdim),ib3=min_hc_state,max_hc_state)
!
!          read(iunit,err=100,iostat=ios) (((hc_subgrid_data(ib1,ib2,ib3),ib1=1,sg_rdim),ib2=1,sg_zdim),ib3=min_hc_state,max_hc_state)
!
! slmod end

       endif
    endif

    !write(0,'(a,10i6)') 'SUBGRID DATA LOADED:',subgrid_opt,sg_hc_opt,sg_rdim,sg_zdim,min_state,max_state,min_hc_state,max_hc_state

    return

100 call errmsg('RELOAD_SUBGRID:','ERROR LOADING SUBGRID DATA FROM RAW FILE')

    return

  end subroutine reload_subgrid




  subroutine norm_subgrid(nizs,cneuta,tneut,tatiz,fsrate,qtim)
    use hc_init_lib_data ! contains charge information for HC states
    implicit none

    integer ::  nizs,cneuta
    real :: tneut,tatiz,fsrate,qtim

    !
    ! Local variables
    !

    real :: total_particles,timestep

    integer :: iz

    !write(0,'(a,3i12,5(1x,g12.5))') 'NORM SUBGRID:',subgrid_opt,nizs,cneuta,tneut,tatiz,fsrate,qtim


    if (subgrid_opt.ne.0) then 


       !write(0,*) 'NORM SUM:',sum(subgrid_data),sum(hc_subgrid_data)

       !
       ! For neutral launches - set the total number of particles to tneut
       ! For ion injections - set it to tatiz
       ! This is required because ion injections may have some neutrals
       ! launched due to recombination and recycling.
       !
       !

       if (cneuta.eq.0) then 
          total_particles = tneut
       elseif (cneuta.eq.1) then 
          total_particles = tatiz
       endif
       
       !write(0,'(a,i5,3(1x,g12.5))') 'NORM TOT:',cneuta,tneut,tatiz,total_particles
       !write(0,'(a,5i10)') 'NORM IZS:',min_state,max_state,min_hc_state,max_hc_state,sg_hc_opt

       !
       ! If the HC subgrid option is ON then we are following hydrocarbons which means that the
       ! general impurity is carbon. This means that the HC CI density data needs to be added
       ! to any density data recorded during the DIVIMP run - perhaps resulting from recombination.
       ! or self-sputtering. 
       !
       ! These data together can then be normalized.
       !
       if (sg_hc_opt.ne.0) then 
          ! Add HC CI data to whatever is recorded in the DIVIMP arrays. 
          ! Neutral carbon is index 2 in the hc_subgrid_data array
          ! Neutral carbon is index 0 in the subgrid_data array
          subgrid_data(:,:,0) = subgrid_data(:,:,0) + hc_subgrid_data(:,:,2)
       endif

       do iz = min_state,max_state

          if (iz.le.0) then 
             timestep = fsrate
          else
             timestep = qtim
          endif

          !
          ! The last argument indicates whether the normalization
          ! should be applied to the regular or hc subgrid
          ! 0 = regular 1=HC
          !

          call normalize_subgrid(iz,timestep,total_particles,0)

       enddo

       if (sg_hc_opt.ne.0) then 

          do iz = min_hc_state,max_hc_state

             !
             ! Note: HC states start indexing at 1 - not 0 - so using a 0 start point for HC is a bug - however to maintain
             !       compatibility with existing output I am not changing min_hc_state at this time
             !
             if (iz.eq.0) then 
                timestep = fsrate
             elseif (get_hc_charge(iz).eq.0) then 
                timestep = fsrate
             else
                timestep = qtim
             endif

             !
             ! The last argument indicates whether the normalization
             ! should be applied to the regular or hc subgrid
             ! 0 = regular 1=HC
             !

             call normalize_subgrid(iz,timestep,total_particles,1)

          enddo


       endif


    endif

  end subroutine norm_subgrid



  subroutine normalize_subgrid(state,timestep,total_particles,subgrid_select)
    implicit none

    !
    ! Subgrid_select choose whether to normalize the regular subgrid (option 0) or the HC subgrid (option 1)
    !

    integer :: state,subgrid_select
    real  :: timestep,total_particles

    !
    ! Local variables
    !

    real :: cellvol ! cell volume  = cell cross sectional area * 1 meter toroidally. 
    real :: norm_fact ! normalization factor

    !
    ! The normalization factor for each state is:
    !     
    !   norm_fact = timestep * 1/total_particles * 1/cellvol
    !

    !
    ! Only do normalizations if the subgrids were turned on - though this routine should not be 
    ! called in the case where subgrid_opt was 0.
    !

    if (subgrid_opt.ne.0) then 

       cellvol = dr * dz

       norm_fact = timestep/(total_particles*cellvol)

       !write(0,'(a,i4,12g12.5)') 'NORMALIZE SUBGRID:',subgrid_select,norm_fact,timestep,total_particles,cellvol,dr,dz


       if (subgrid_select.eq.0.and.state.ge.min_state.and.state.le.max_state) then 

          subgrid_data(:,:,state) = subgrid_data(:,:,state) * norm_fact

       elseif (subgrid_select.eq.1.and.sg_hc_opt.ne.0.and.state.ge.min_hc_state.and.state.le.max_hc_state) then 

          hc_subgrid_data(:,:,state) = hc_subgrid_data(:,:,state) * norm_fact

       endif

    endif

  end subroutine normalize_subgrid

  subroutine subgrid_getstateranges(start_state,end_state,start_hc_state,end_hc_state)
    implicit none

    integer :: start_state,end_state,start_hc_state,end_hc_state

    start_state = min_state
    end_state= max_state
    start_hc_state= min_hc_state
    end_hc_state = max_hc_state

  end subroutine subgrid_getstateranges

  subroutine load_subgrid_xsection(subgrid_xsec,iselect,istate)
    implicit none
    real :: subgrid_xsec(sg_rdim,sg_zdim)
    integer :: iselect,istate

    integer :: ir,iz,is
    !
    ! Out of range state values return the sum over all states
    !


    if (iselect.eq.32.or.iselect.eq.34) then 
    ! 
    ! Load regular data xsection
    !

       if (istate.ge.min_state.and.istate.le.max_state) then

          subgrid_xsec = subgrid_data(:,:,istate)

       else
          !
          ! Return sum over states
          !
          subgrid_xsec = 0.0
          do ir = 1,sg_rdim
             do iz = 1,sg_zdim
                do is = min_state,max_state
                   subgrid_xsec(ir,iz) = subgrid_xsec(ir,iz) + subgrid_data(ir,iz,is)
                enddo
             enddo
           enddo

       endif

    elseif (iselect.eq.33.or.iselect.eq.35) then 
    !
    ! Load HC data xsection
    !
    ! 
    ! Load regular data xsection
    !

       if (istate.ge.min_hc_state.and.istate.le.max_hc_state) then

          subgrid_xsec = hc_subgrid_data(:,:,istate)

       else
          !
          ! Return sum over states
          !
          subgrid_xsec = 0.0
          do ir = 1,sg_rdim
             do iz = 1,sg_zdim
                do is = min_hc_state,max_hc_state
                   subgrid_xsec(ir,iz) = subgrid_xsec(ir,iz) + hc_subgrid_data(ir,iz,is)
                enddo
             enddo
           enddo

       endif

    endif

  end subroutine load_subgrid_xsection

  subroutine calc_subgrid_axis(axis,min,max,npts)
    implicit none
    integer :: npts
    real :: min,max,axis(npts)

    integer :: in
    real :: step
    !
    ! Calculate the axis coordinates based on min, max and npts 
    !

    step = (max-min)/real(npts)

    do in = 1,npts

       axis(in) = min + (real(in) - 0.5) * step

    end do
    
  end subroutine calc_subgrid_axis



end module subgrid
