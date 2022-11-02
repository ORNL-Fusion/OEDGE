module mod_out_unstruc
  use debug_options
  implicit none

  private

  !
  !     this common block contains variables used in out and read in using the
  !     unstructured input methodology
  !
  !     -*-fortran-*-
  ! common /out_unstructured_input/ new_absfac,psi1_reg,psi2_reg
  ! save /out_unstructured_input/
  real,public :: new_absfac
  integer,public :: absfac_opt,e2dizs_offset  
  real ,public :: psi1_reg,psi2_reg
  integer,public :: absfac_iz, absfac_ir, absfac_ikstart, absfac_ikend
  integer,public :: scale_1d
  
  public :: allocate_mod_out_unstruc,deallocate_mod_out_unstruc,init_out_unstruc_input
  real,public,allocatable :: transport_area(:,:)
  
contains

  subroutine allocate_mod_out_unstruc
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_out_unstruc','ALLOCATE')
    call allocate_array(transport_area,maxnks,maxnrs,'transport_area',ierr)

  end subroutine allocate_mod_out_unstruc


  subroutine deallocate_mod_out_unstruc
    implicit none

    call pr_trace('mod_out_unstruc','DEALLOCATE')
    if (allocated(transport_area)) deallocate(transport_area)

  end subroutine deallocate_mod_out_unstruc



  subroutine init_out_unstruc_input
    use mod_params
    use mod_comtor
    implicit none
    ! 
    !     This routine assigns the default values to any
    !     OUT unstructured input values.The OUT unstructured
    !     input tags start with the letter "O" (oh)
    !
    !     jdemod - it also assigns default values to any other 
    !              unstructured input values shared by DIV and OUT
    !
    !
    !     Note: The netcdf output option A07 from DIVIMP is also supported
    !           since it would be useful to create a netcdf version of the 
    !           raw output from an OUT run
    !
    !     include 'params'
    !
    !     include 'out_unstruc'
    !     include 'comtor'
    !     
    !------------------------------------------------------
    !
    !     O01 - alternate absolute factor specification 
    !         - this option is deactivated by setting the 
    !           default value to zero.  
    !
    new_absfac=0.0
    !
    !
    !------------------------------------------------------
    !
    !     Core fueling code calculates integrated ionization
    !     profiles in the core ... these parameters allow the 
    !     PSIN inner bound of the integration regions to be set
    !     This is used in the pr_eirene_analysis routine
    !
    !     O02 - PSIN bound for calculating core ionization 
    !           profile 1 (psi1_reg)
    !     O03 - PSIN bound for calculating core ionization 
    !           profile 2 (psi2_reg)
    !

    psi1_reg = 0.9
    psi2_reg = 0.95


    ! -----------------------------------------------------------------------
    !
    !     O04:
    !     absfac_opt - options specifying how to calculate scaling factor for DIVIMP
    !                  results in OUT         
    !
    !     This option allows the absolute scaling factor for the DIVIMP
    !     run results to be specified in the OUT routine. It's default
    !     value is zero.
    !      

    absfac_opt = 0

    ! -----------------------------------------------------------------------
    !
    !     O05:
    !     e2dizs_offset - offset to match fluid code impurity charge state results data to 
    !                     DIVIMP results (needed sometimes since the fluid code results
    !                     sometimes contain multiple fluids). 
    !

    e2dizs_offset = 0


    ! -----------------------------------------------------------------------
    !
    !     O06->O09:
    !
    !     Indices to control the ABSFAC calculation
    !     Default values are initially zero. These will be reset to specific default
    !     values depending on the absfac_opt option when the code is executed
    !

    absfac_iz = 0       ! goes to maxizs
    absfac_ir = 0       ! goes to irsep-1
    absfac_ikstart = 0  ! goes to 1
    absfac_ikend = 0    ! goes to nks(absfac_ir) or nks(absfac_ir)-1 inside separatrix

    ! -----------------------------------------------------------------------
    !
    !  O10: Scale plots for 1D (convert from kareas to distance along the field line
    !       0 = off,
    !       1 = on (scale by distance along the field line)
    !       2 = on (scale by transport cell area instead of mesh cell area)
    scale_1d = 0

    
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG A07: 
    !
    !
    !     Option to write a netcdf version of the raw data file 
    !     0 = off   1 = 0n .. default OFF
    !
    !
    netcdf_opt = 0
    !
    !
    !------------------------------------------------------
    !
    return
  end subroutine init_out_unstruc_input




end module mod_out_unstruc
