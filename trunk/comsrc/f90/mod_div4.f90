module mod_div4
  use debug_options
  implicit none

  
  
  !     -*-fortran-*-
  ! common /div4c/rizpos,zizpos,recstruk,recmain,recexit,recatiz,recneut,recwalln,reccent,&
  !     rectmax,recfail,implim,mizs,reflect_ion,refflag,refstruk,refmain,refexit,&
  !     refatiz,refneut,refwalln,refcent,reftmax,reffail,mtcrefstruk,mtcrefwalln,refloss,&
  !     refprob,hasleaked,hasleakedcore,cleakp,isol,ifp,irflct,rstart,zstart,totleakcore,&
  !     totsource,totsleak,griderr
  
  ! save /div4c/
  
  
  real,public :: rizpos,zizpos,recstruk,recmain,recexit,recatiz,recneut,recwalln,reccent,&
       rectmax,recfail
  
  
  
  
  
  
  
  
  !
  !     reflected impurities
  !
  integer,public :: implim,mizs
  logical,public :: reflect_ion
  integer,public :: refflag
  real,public :: refstruk,refmain,refexit,refatiz,refneut,refwalln,refcent,reftmax,&
       reffail
  !
  !     leakage variables
  !
  real,public :: mtcrefstruk,mtcrefwalln,refloss,refprob
  logical,public :: hasleaked,hasleakedcore
  integer,public :: cleakp,isol,ifp,irflct
  !
  real,public :: rstart,zstart,totleakcore
  !
  real,public :: totsource,totsleak
  
  logical,public :: griderr

  public :: allocate_mod_div4,deallocate_mod_div4

contains

  subroutine allocate_mod_div4
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div4','ALLOCATE')


  end subroutine allocate_mod_div4


  subroutine deallocate_mod_div4
    implicit none

    call pr_trace('mod_div4','DEALLOCATE')


  end subroutine deallocate_mod_div4

end module mod_div4