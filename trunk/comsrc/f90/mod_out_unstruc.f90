module mod_out_unstruc
  use debug_options
  implicit none

  !
  !     this common block contains variables used in out and read in using the
  !     unstructured input methodology
  !
  !     -*-fortran-*-
  ! common /out_unstructured_input/ new_absfac,psi1_reg,psi2_reg
  ! save /out_unstructured_input/
  real,public :: new_absfac
  real*8 ,public :: psi1_reg,psi2_reg

  public :: allocate_mod_out_unstruc,deallocate_mod_out_unstruc

contains

  subroutine allocate_mod_out_unstruc
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_out_unstruc','ALLOCATE')


  end subroutine allocate_mod_out_unstruc


  subroutine deallocate_mod_out_unstruc
    implicit none

    call pr_trace('mod_out_unstruc','DEALLOCATE')


  end subroutine deallocate_mod_out_unstruc

end module mod_out_unstruc