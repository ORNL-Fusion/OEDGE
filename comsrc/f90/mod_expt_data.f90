module mod_expt_data
  use mod_params
  use debug_options
  implicit none

  !
  !     this common block contains arrays which hold sets
  !     of experimental data for plotting by the draw routine or
  !     for general use within out.
  !
  !     usually, the experimental data can be loaded directly by
  !     draw - however, the code does not allow arbitrary changes
  !     to be made to the data like adjustment by scaling factors or
  !     variable axis offsets to account to calibration
  !     uncertainties and the like.
  !
  !     this common block allows the plot data to be loaded, modified,
  !     and then made accessible to the plot code. this data will be
  !     used instead of loading data directly from the experimental data
  !     file.
  !
  !     another limitation of the draw routine is that all the data
  !     must share a common axis - this is often not practicable with
  !     experimental data unless interpolation is performed which can
  !     give misleading results. the draw routine should be modified at
  !     some time to allow for different axis data for each set of plot
  !     data.
  !
  !     -*-fortran-*-
  ! common /expt_data/ expt_data_available,expt_data_num,expt_data_ncols,expt_data_axis_type,&
  !     expt_data_axis,expt_data_values,expt_data_title
  !
  ! save /expt_data/
  
  integer,public :: expt_data_maxcols
  integer,public :: expt_data_maxdatx
  integer,public :: expt_dataunit
  !
  !parameter(expt_data_maxcols=1,expt_data_maxdatx=maxdatx, expt_dataunit=exptunit)
  parameter(expt_data_maxcols=1, expt_dataunit=exptunit)
  !
  integer,public :: expt_data_axis_type,expt_data_num,expt_data_ncols,expt_data_available
  real,public,allocatable :: expt_data_axis(:),expt_data_values(:,:)
  !
  character*100,public :: expt_data_title

  public :: allocate_mod_expt_data,deallocate_mod_expt_data

contains

  subroutine allocate_mod_expt_data
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    expt_data_maxdatx = maxdatx
    
    call pr_trace('mod_expt_data','ALLOCATE')

    call allocate_array(expt_data_axis,expt_data_maxdatx,'expt_data_axis',ierr)
    call allocate_array(expt_data_values,expt_data_maxdatx,expt_data_maxcols,'expt_data_values',&
         ierr)

  end subroutine allocate_mod_expt_data


  subroutine deallocate_mod_expt_data
    implicit none

    call pr_trace('mod_expt_data','DEALLOCATE')

    if (allocated(expt_data_axis)) deallocate(expt_data_axis)
    if (allocated(expt_data_values)) deallocate(expt_data_values)

  end subroutine deallocate_mod_expt_data

end module mod_expt_data
