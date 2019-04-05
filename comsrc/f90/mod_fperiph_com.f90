module mod_fperiph_com
  use debug_options
  implicit none

  !
  !     these common blocks contain data relevant to the ion far
  !     periphery options.
  !
  !     -*-fortran-*-
  integer,parameter ,public :: max_num_fp_regions = 4
  integer,parameter ,public :: fp_main = 1
  integer,parameter ,public :: fp_pfz  = 2
  ! below are to support double null grids and fp option 5
  integer,parameter ,public :: fp_main2  = 3
  !
  integer,parameter ,public :: fp_pfz2  = 4
  ! common /fp_com/ fpxmaxo,fptimo,fpxmaxi,fptimi,cdperpfp,fpopt,fpropt
  ! save /fp_com/
  integer,public :: fpopt,fpropt
  !
  !
  !
  real,public :: fpxmaxo,fptimo,fpxmaxi,fptimi,cdperpfp
  ! common /fp_neut/ fp_neut_opt,fp_plasma_opt,fp_flow_opt,fp_te,fp_ne,fp_flow_velocity_input,&
  !     fp_neut_ent,fp_neut_ioniz,fp_neut_wall,fp_neut_tmax,fp_neut_targ,fp_neut_data
  !
  ! save /fp_neut/
  integer,public :: fp_neut_opt,fp_plasma_opt,fp_flow_opt
  !
  !
  !
  real,public :: fp_te,fp_ne,fp_neut_ent,fp_neut_ioniz,fp_neut_wall,fp_neut_tmax,fp_neut_targ,&
       fp_flow_velocity_input
  real,public,allocatable :: fp_neut_data(:,:)
  !     >                ,fp_sdrft_start,fp_sdrft_end
  ! common /fp_data/ num_fp_regions,fp_plasma,fp_virt_rings,fp_region,fp_n_bins,fp_grid_width_opt,&
  !     fp_s,fp_rings,fp_cells,fp_walldist,min_fp_walldist,fp_wallcoords,&
  !     fp_flow_velocity,fp_maxdist,fp_cross_tmp
  !
  ! save /fp_data/
  integer,public :: num_fp_regions,fp_region,fp_n_bins,fp_grid_width_opt
  integer,public,allocatable :: fp_virt_rings(:),fp_rings(:),fp_cells(:)
  !     >    ,fp_sdrft_start(num_fp_regions),
  !     >    ,fp_sdrft_end(num_fp_regions)
  real,public :: fp_cross_tmp
  real,public,allocatable :: fp_plasma(:,:,:),min_fp_walldist(:,:),fp_walldist(:,:),&
       fp_wallcoords(:,:,:),fp_flow_velocity(:),fp_maxdist(:)
  
  !
  !     hard coded flag for debugging the new fp options
  !
  real,public,allocatable :: fp_s(:,:)
  !
  !     contents of fp_plasma array:
  !
  !     fp_plasma(ik,ireg,1) = ne
  !     fp_plasma(ik,ireg,2) = te
  !     fp_plasma(ik,ireg,3) = ti
  !     fp_plasma(ik,ireg,4) = vb
  !     fp_plasma(ik,ireg,5) = efield
  !     fp_plasma(ik,ireg,6) = tgrade
  !     fp_plasma(ik,ireg,7) = tgradi
  !
  logical ,public :: debug_fp = .false.

  public :: allocate_mod_fperiph_com,deallocate_mod_fperiph_com

contains

  subroutine allocate_mod_fperiph_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_fperiph_com','ALLOCATE')

    call allocate_array(fp_neut_data,maxpts+1,5,'fp_neut_data',ierr)
    call allocate_array(fp_virt_rings,max_num_fp_regions,'fp_virt_rings',ierr)
    call allocate_array(fp_rings,max_num_fp_regions,'fp_rings',ierr)
    call allocate_array(fp_cells,max_num_fp_regions,'fp_cells',ierr)
    call allocate_array(fp_plasma,maxnks,max_num_fp_regions,7,'fp_plasma',ierr)
    call allocate_array(min_fp_walldist,maxnks,max_num_fp_regions,'min_fp_walldist',ierr)
    call allocate_array(fp_walldist,0,2*maxnks,1,max_num_fp_regions,'fp_walldist',ierr)
    call allocate_array(fp_wallcoords,0,2*maxnks,1,max_num_fp_regions,1,2,'fp_wallcoords',&
         ierr)
    call allocate_array(fp_flow_velocity,max_num_fp_regions,'fp_flow_velocity',ierr)
    call allocate_array(fp_maxdist,max_num_fp_regions,'fp_maxdist',ierr)
    call allocate_array(fp_s,0,2*maxnks,1,max_num_fp_regions,'fp_s',ierr)

  end subroutine allocate_mod_fperiph_com


  subroutine deallocate_mod_fperiph_com
    implicit none

    call pr_trace('mod_fperiph_com','DEALLOCATE')

    if (allocated(fp_neut_data)) deallocate(fp_neut_data)
    if (allocated(fp_virt_rings)) deallocate(fp_virt_rings)
    if (allocated(fp_rings)) deallocate(fp_rings)
    if (allocated(fp_cells)) deallocate(fp_cells)
    if (allocated(fp_plasma)) deallocate(fp_plasma)
    if (allocated(min_fp_walldist)) deallocate(min_fp_walldist)
    if (allocated(fp_walldist)) deallocate(fp_walldist)
    if (allocated(fp_wallcoords)) deallocate(fp_wallcoords)
    if (allocated(fp_flow_velocity)) deallocate(fp_flow_velocity)
    if (allocated(fp_maxdist)) deallocate(fp_maxdist)
    if (allocated(fp_s)) deallocate(fp_s)

  end subroutine deallocate_mod_fperiph_com

end module mod_fperiph_com