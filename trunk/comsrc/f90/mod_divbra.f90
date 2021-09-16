module mod_divbra
  use debug_options
  implicit none

  !
  !     this file contains the values that are passed from divimp
  !     to b2/eirene. it can be included directly in those source codes
  !     if the varaible declarations do not overlap.
  !
  !     note: the parametric sizes of these arrays must be larger than
  !           the actual values used in divimp for the specific case or
  !           errors will result due to array bounds overruns. one solution
  !           is to ensure that the values used here are larger than the
  !           corresponding divimp parameters.
  !
  !     david elder, sept 6, 1994
  !
  !     -*-fortran-*-
  integer,public :: dbmaxnrs,dbmaxnks,dbmaxizs
  !
  !     ipp/01 krieger: changed dbmaxizs from 18 to 74 (tungsten)
  !
  !      parameter (dbmaxnks=100, dbmaxnrs= 50,dbmaxizs=74)
  !
  !     the following is the declaration of the common block and its contents
  !     that will be passed to b2-eirene
  !
  parameter (dbmaxnks=100, dbmaxnrs= 50,dbmaxizs=18)
  ! common /divbra/ dbni, dbti, dbv, dbpowls, dblines
  !
  ! save /divbra/
  
  !
  double precision,public,allocatable :: dbni(:,:,:),dbti(:,:,:),dbv(:,:,:),dbpowls(:,:,:),&
       dblines(:,:,:)

  public :: allocate_mod_divbra,deallocate_mod_divbra

contains

  subroutine allocate_mod_divbra
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_divbra','ALLOCATE')

    call allocate_array(dbni,1,dbmaxnks,1,dbmaxnrs,-1,dbmaxizs,'dbni',ierr)
    call allocate_array(dbti,1,dbmaxnks,1,dbmaxnrs,-1,dbmaxizs,'dbti',ierr)
    call allocate_array(dbv,1,dbmaxnks,1,dbmaxnrs,-1,dbmaxizs,'dbv',ierr)
    call allocate_array(dbpowls,1,dbmaxnks,1,dbmaxnrs,-1,dbmaxizs,'dbpowls',ierr)
    call allocate_array(dblines,1,dbmaxnks,1,dbmaxnrs,-1,dbmaxizs,'dblines',ierr)

  end subroutine allocate_mod_divbra


  subroutine deallocate_mod_divbra
    implicit none

    call pr_trace('mod_divbra','DEALLOCATE')

    if (allocated(dbni)) deallocate(dbni)
    if (allocated(dbti)) deallocate(dbti)
    if (allocated(dbv)) deallocate(dbv)
    if (allocated(dbpowls)) deallocate(dbpowls)
    if (allocated(dblines)) deallocate(dblines)

  end subroutine deallocate_mod_divbra

end module mod_divbra