module mod_expt_data

  use mod_params

!c
!c     This Common block contains arrays which hold sets 
!c     of experimental data for plotting by the DRAW routine or
!c     for general use within OUT. 
!c    
!c     Usually, the experimental data can be loaded directly by
!c     DRAW - however, the code does not allow arbitrary changes
!c     to be made to the data like adjustment by scaling factors or
!c     variable axis offsets to account to calibration 
!c     uncertainties and the like. 
!c
!c     This common block allows the plot data to be loaded, modified,
!c     and then made accessible to the plot code. This data will be 
!c     used instead of loading data directly from the experimental data
!c     file. 
!c
!c     Another limitation of the DRAW routine is that all the data
!c     must share a common axis - this is often not practicable with 
!c     experimental data unless interpolation is performed which can 
!c     give misleading results. The DRAW routine should be modified at
!c     some time to allow for different axis data for each set of plot
!c     data.
!c
!      common /expt_data/ expt_data_available,
!     >         expt_data_num,expt_data_ncols,expt_data_axis_type,
!     >         expt_data_axis,expt_data_values,expt_data_title
!c
!        
!      integer expt_data_maxcols
!      integer expt_data_maxdatx  
!      integer expt_dataunit
!      parameter(expt_data_maxcols=1,expt_data_maxdatx=maxdatx,
!     >          expt_dataunit=exptunit)
!c
!      integer expt_data_axis_type,expt_data_num,
!     >           expt_data_ncols,expt_data_available
!c
!      real expt_data_axis(expt_data_maxdatx),
!     >     expt_data_values(expt_data_maxdatx,expt_data_maxcols)
!      character*100 expt_data_title
!c
!c
!c---------------------------------------------------------------------
!c
!c     Experimental data related values 
!c     - used to hold a list of experimantal dataset indices to be 
!c       included on plots if possible
!c
!c---------------------------------------------------------------------
!c
!      integer max_expt_datasets
!      parameter (max_expt_datasets=10)
!      common /expt_data_list/ expt_nsets,expt_datasets
!      integer expt_nsets, expt_datasets(max_expt_datasets)
!c

  implicit none
  private

        
      integer,parameter,public:: expt_data_maxcols=1
      integer,parameter,public:: expt_data_maxdatx=maxdatx
      integer,parameter,public:: expt_dataunit=exptunit
      !parameter(expt_data_maxcols=1,expt_data_maxdatx=maxdatx,expt_dataunit=exptunit)

      integer,public:: expt_data_axis_type,expt_data_num,expt_data_ncols,expt_data_available

      real,public:: expt_data_axis(expt_data_maxdatx),expt_data_values(expt_data_maxdatx,expt_data_maxcols)
      character*100,public:: expt_data_title

      integer,parameter,public::  max_expt_datasets=10
      !parameter (max_expt_datasets=10)

      integer,public:: expt_nsets, expt_datasets(max_expt_datasets)

  
  public :: allocate_mod_expt_data, deallocate_mod_expt_data


contains

  subroutine allocate_mod_expt_data
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_expt_data


  subroutine deallocate_mod_expt_data
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_expt_data



end module mod_expt_data
