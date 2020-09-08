module mod_inel
  use debug_options
  implicit none

  !
  !    the following common block is used to read in the inel rate data
  !
  !     ntev  - number of temperature values - hard coded
  !     nz    - number of charge states - also hard-coded
  !     tevb  - temperatures for rate arrays
  !     rsi   - ionization rate
  !     rre   - recombination rate
  !     rpwr  - radiative power rate
  !     rrcx  - cx recombination rate
  !
  !
  !     -*-fortran-*-
  ! common /inel/ nz,tevb,rsi,rre,rpwr,rrcx
  !
  ! save /inel/
  integer,public :: ntev,nz
  !
  parameter (ntev=101)
  real,public,allocatable :: tevb(:)
  real,public,allocatable :: rsi(:,:)
  real,public,allocatable :: rre(:,:)
  real,public,allocatable :: rpwr(:,:)
  !
  real,public,allocatable :: rrcx(:,:)

  public :: allocate_mod_inel,deallocate_mod_inel

contains

  subroutine allocate_mod_inel
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_inel','ALLOCATE')

    call allocate_array(tevb,ntev,'tevb',ierr)
    call allocate_array(rsi,1,ntev,0,maxizs-1,'rsi',ierr)
    call allocate_array(rre,ntev,maxizs,'rre',ierr)
    call allocate_array(rpwr,1,ntev,0,maxizs,'rpwr',ierr)
    call allocate_array(rrcx,ntev,maxizs,'rrcx',ierr)

  end subroutine allocate_mod_inel


  subroutine deallocate_mod_inel
    implicit none

    call pr_trace('mod_inel','DEALLOCATE')

    if (allocated(tevb)) deallocate(tevb)
    if (allocated(rsi)) deallocate(rsi)
    if (allocated(rre)) deallocate(rre)
    if (allocated(rpwr)) deallocate(rpwr)
    if (allocated(rrcx)) deallocate(rrcx)

  end subroutine deallocate_mod_inel

end module mod_inel