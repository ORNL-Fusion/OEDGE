module mod_fp_data

  ! 
  ! fp grid variables
  !
  real, allocatable :: fp_density(:,:,:,:)
  real, allocatable :: fp_grid_area(:,:)
  integer, allocatable :: fp_grid_flag(:,:,:)
  real, allocatable :: fp_grid_dist(:,:)

  real, allocatable :: fp_grid_plasma(:,:,:,:)
  
  !integer, parameter :: fp_n_bins = 20
  real,allocatable :: fp_r_bin_width(:)


contains


  subroutine fp_allocate_storage(ierr)
    use mod_params
    use error_handling
    use mod_fperiph_com
    implicit none
    integer :: ierr
    
    ierr = 0

    ! density
    if (allocated(fp_density)) deallocate(fp_density)
    allocate(fp_density(maxnks,fp_n_bins+1,0:maxizs,num_fp_regions),stat=ierr)
    fp_density = 0.0

    ! "area" for each cell ... this is very approximate and just based on the area of the cell at the grid edge modified for change in cell width
    if (allocated(fp_grid_area)) deallocate(fp_grid_area)
    allocate(fp_grid_area(maxnks,num_fp_regions),stat=ierr)
    fp_grid_area = 0.0

    ! this array flags cells in the fp mesh that are outside the wall
    if (allocated(fp_grid_flag)) deallocate(fp_grid_flag)
    allocate(fp_grid_flag(maxnks,fp_n_bins+1,num_fp_regions),stat=ierr)    
    fp_grid_flag = 0

    ! This specifies the radial distance/width of the cells from the grid edge
    if (allocated(fp_grid_dist)) deallocate(fp_grid_dist)
    allocate(fp_grid_dist(fp_n_bins+1,num_fp_regions),stat=ierr)
    fp_grid_dist = 0

    ! R bin width in each region
    if (allocated(fp_r_bin_width)) deallocate(fp_r_bin_width)
    allocate(fp_r_bin_width(num_fp_regions),stat=ierr)
    fp_r_bin_width = 0.0

    ! Grid for detailed FP plasma 
    if (allocated(fp_grid_plasma)) deallocate(fp_grid_plasma)
    ! contains ne,Te,Ti,vb,E,fig,feg data for each cell of peripheral grid
    allocate(fp_grid_plasma(maxnks,fp_n_bins+1,num_fp_regions,7),stat=ierr)
    !allocate(fp_grid_plasma(maxnks,fp_n_bins+1,max_num_fp_regions,7),stat=ierr)
    fp_grid_plasma = 0.0

    if (ierr.ne.0) then 
       call errmsg('MOD_FP_DATA: FP_ALLOCATE_STORAGE: ERROR IN ALLOCATION:',ierr)
    endif

    !write(0,'(a,i8,10(1x,i8))') 'ALLOCATING FP:',ierr,maxnks,fp_n_bins,num_fp_regions,maxizs

  end subroutine fp_allocate_storage



  subroutine fp_deallocate_storage
    implicit none
    if (allocated(fp_density)) deallocate(fp_density)
    if (allocated(fp_grid_area)) deallocate(fp_grid_area)
    if (allocated(fp_grid_flag)) deallocate(fp_grid_flag)
    if (allocated(fp_grid_dist)) deallocate(fp_grid_dist)
    if (allocated(fp_r_bin_width)) deallocate(fp_r_bin_width)
  end subroutine fp_deallocate_storage


  subroutine fp_write_raw(unit)
    use mod_params
    use error_handling
    use mod_fperiph_com
    implicit none
    integer :: unit
    integer :: in

    write(unit) fpopt, num_fp_regions, fp_n_bins


    !
    ! There are calls to STORE before the fp setup code is run. Steve stores intermediate results
    ! when EIRENE finishes for example. However, this means that the following arrays are not
    ! guaranteed to be allocated at this early point. Rather than stopping with an error message
    ! and ending the run. The code will check if fp_density is allocated before saving. 
    !
    
    if ((fpopt.eq.5.or.fpopt.eq.6).and.allocated(fp_density)) then 

       if (allocated(fp_density)) then 
          call rinout('W FP_DEN',fp_density,maxnks*(fp_n_bins+1)*(maxizs+1)*num_fp_regions)
       else
          call errmsg('FP_WRITE_RAW','FP_DENSITY IS NOT ALLOCATED BUT IT SHOULD BE')
          stop 'FP_WRITE_RAW: fp_density not allocated'
       endif

       if (allocated(fp_grid_area)) then 
          call rinout('W FPAREA',fp_grid_area,maxnks*num_fp_regions)
       endif

       if (allocated(fp_grid_flag)) then 
          call rinout('W FPFLAG',fp_grid_flag,maxnks*(fp_n_bins+1)*num_fp_regions)
       endif

       if (allocated(fp_grid_dist)) then 
          call rinout('W FPDIST',fp_grid_dist,(fp_n_bins+1)*num_fp_regions)
       endif

       if (allocated(fp_r_bin_width)) then 
          call rinout('W FPRWID',fp_r_bin_width,num_fp_regions)
       endif

       ! fp_plasma is in fperiph common block and is not allocated
       call rinout('W FPPLAS',fp_plasma,maxnks*max_num_fp_regions*7)

       ! fp_cells is in fperiph common block and is not allocated
       call iinout('W FPCELLS',fp_cells,max_num_fp_regions)

       !write(6,'(a,i8,20(1x,i8))') 'Writing fp_cells:',max_num_fp_regions,(fp_cells(in),in=1,max_num_fp_regions)

       ! Need to change this to num_fp_regions
       if (allocated(fp_grid_plasma)) then 
          call rinout('W FPGPLAS',fp_grid_plasma,maxnks*(fp_n_bins+1)*num_fp_regions*7)
       endif


    endif

  end subroutine fp_write_raw

  subroutine fp_read_raw(unit,version_code,maxrev)
    use mod_params
    use error_handling
    use debug_options
    use mod_fperiph_com
    implicit none
    integer :: unit,ierr
    integer,intent(in) :: version_code,maxrev
    integer :: in,ik,ic,id
    
    read(unit) fpopt, num_fp_regions, fp_n_bins
    
    if (fpopt.eq.5.or.fpopt.eq.6) then 

       call fp_allocate_storage(ierr)

       if (allocated(fp_density)) then 
          call rinout('R FP_DEN',fp_density,maxnks*(fp_n_bins+1)*(maxizs+1)*num_fp_regions)
       else
          call errmsg('FP_READ_RAW','FP_DENSITY IS NOT ALLOCATED BUT IT SHOULD BE')
          stop 'FP_READ_RAW: fp_density not allocated'
       endif

       if (allocated(fp_grid_area)) then 
          call rinout('R FPAREA',fp_grid_area,maxnks*num_fp_regions)
       endif

       if (allocated(fp_grid_flag)) then 
          call rinout('R FPFLAG',fp_grid_flag,maxnks*(fp_n_bins+1)*num_fp_regions)
       endif

       if (allocated(fp_grid_dist)) then 
          call rinout('R FPDIST',fp_grid_dist,(fp_n_bins+1)*num_fp_regions)
       endif

       if (allocated(fp_r_bin_width)) then 
          call rinout('R FPRWID',fp_r_bin_width,num_fp_regions)
       endif

       ! fp_plasma is in fperiph common block and is not allocated
       call rinout('R FPPLAS',fp_plasma,maxnks*max_num_fp_regions*7)

       call pr_trace('FP_READ_RAW','BEFORE READ FP_GRID_PLASMA')

       if (version_code.ge.6*maxrev+49) then

          ! fp_cells is in fperiph common block and is not allocated
          call iinout('R FPCELLS',fp_cells,max_num_fp_regions)

          !write(6,'(a,i8,20(1x,i8))') 'READING fp_cells:',max_num_fp_regions,(fp_cells(in),in=1,max_num_fp_regions)

          if (allocated(fp_grid_plasma)) then 
             call rinout('R FPGPLAS',fp_grid_plasma,maxnks*(fp_n_bins+1)*num_fp_regions*7)
             !call rinout('R FPGPLAS',fp_grid_plasma,maxnks*(fp_n_bins+1)*max_num_fp_regions*7)
          endif

       endif
       call pr_trace('FP_READ_RAW','AFTER READ FP_GRID_PLASMA')

       ! write out fp_grid_plasma

       write(6,'(a,20(1x,i8))') 'WRITING FP PLASMA:',(fp_cells(in),in=1,num_fp_regions)

         do in = 1,num_fp_regions
             do ik =1, fp_cells(in)
               write(6,'(a,2i6,7g12.5)') 'FP PLASMA:',in,ik, (fp_plasma(ik,in,id),id=1,7)
            end do
         end do

         do in = 1,num_fp_regions
            do ic = 1,fp_n_bins
               do ik =1, fp_cells(in)
                 write(6,'(a,3i6,7g12.5)') 'FP GRID PLASMA:',in,ic,ik,(fp_grid_plasma(ik,ic,in,id),id=1,7)
               end do
            end do
         end do



    endif

  end subroutine fp_read_raw


  subroutine fp_write_netcdf
    use mod_params
    use mod_fperiph_com
    use nc_utils_generic
    implicit none
    
    integer :: ierr
    !
    ! Write out raw data
    !
    ! fpopt
    ! num_fp_regions
    ! fp_n_bins
    ! 
    ! if the data is not allocated it is not written
    !

    ierr = write_nc('FPOPT',fpopt,'Ion periphery transport option')
    ierr = write_nc('NUM_FP_REGIONS',num_fp_regions,'Number of peripheral transport regions around the grid')
    
    if (fpopt.eq.5.or.fpopt.eq.6) then 
       ierr = write_nc('FP_N_BINS',fp_n_bins,'Number of radial bins in the periphery')

       if (allocated(fp_density)) then 
          ierr = write_nc('FP_DENSITY',fp_density,['MAXNKS  ','FPNBINP1','MAXIZSP1','NFPREG  '],[maxnks,fp_n_bins+1,maxizs+1,num_fp_regions],'Density in each cell of the periphery grid')
       endif

       if (allocated(fp_grid_area)) then 
          ierr = write_nc('FP_GRID_AREA',fp_grid_area,['MAXNKS','NFPREG'],[maxnks,num_fp_regions],'Area for each "cell" along the ring in the periphery grid')
       endif

       if (allocated(fp_grid_flag)) then 
          ierr = write_nc('FP_GRID_FLAG',fp_grid_flag,['MAXNKS  ','FPNBINP1','NFPREG  '],[maxnks,fp_n_bins+1,num_fp_regions],'Flag indicating whether each associated cell is inside the wall')
       endif

       if (allocated(fp_grid_dist)) then 
          ierr = write_nc('FP_GRID_DIST',fp_grid_dist,['FPNBINP1','NFPREG  '],[fp_n_bins+1,num_fp_regions],'Radial distance from the edge of the grid to the center of ring of cells')
       endif

       if (allocated(fp_r_bin_width)) then 
          ierr = write_nc('FP_R_BIN_WIDTH',fp_r_bin_width,['NFPREG'],[num_fp_regions],'Base radial bin width in each FP region')
       endif

       ! fp_plasma is in fperiph common block and is not allocated
       !call rinout('W FPPLAS',fp_plasma,maxnks*max_num_fp_regions*7)

       ierr = write_nc('FP_PLASMA',fp_plasma,['MAXNKS','MAXNFP','7     '],[maxnks,max_num_fp_regions,7],'Plasma data associated with each FP region')

       if (allocated(fp_grid_plasma)) then 
          ierr = write_nc('FP_GRID_PLASMA',fp_grid_plasma,['MAXNKS  ','FPNBINP1','NFPREG  ','7       '],[maxnks,fp_n_bins+1,num_fp_regions,7],'Plasma data associated with each FP mesh')
       endif


    endif

  end subroutine fp_write_netcdf




  subroutine fp_get_plasma(ik,in,fp_reg,ne,te,ti,vb,ef,tgrade,tgradi)
    use mod_fperiph_com
    implicit none
    integer ik,fp_reg,in
    real,intent(out) :: ne,te,ti,ef,vb,tgrade,tgradi
    !
    !      include 'params'
    !      include 'fperiph_com'
    !
    !     Return the fp plasma for the specified cell on the fp grid
    !
    ne = fp_grid_plasma(ik,in,fp_reg,1)
    te = fp_grid_plasma(ik,in,fp_reg,2)
    ti = fp_grid_plasma(ik,in,fp_reg,3)
    vb = fp_grid_plasma(ik,in,fp_reg,4)
    ef = fp_grid_plasma(ik,in,fp_reg,5)
    tgrade = fp_grid_plasma(ik,in,fp_reg,6)
    tgradi = fp_grid_plasma(ik,in,fp_reg,7)

    !
    ! Temporary fallback for cases where fp_grid_plasma is not availabe
    !

    if (ne.le.0.0.and.te.le.0.0) then
       ne = fp_plasma(ik,fp_reg,1)
       te = fp_plasma(ik,fp_reg,2)
       ti = fp_plasma(ik,fp_reg,3)
       vb = fp_plasma(ik,fp_reg,4)
       ef = fp_plasma(ik,fp_reg,5)
       tgrade = fp_plasma(ik,fp_reg,6)
       tgradi = fp_plasma(ik,fp_reg,7)
    endif


  end subroutine fp_get_plasma





end module mod_fp_data
