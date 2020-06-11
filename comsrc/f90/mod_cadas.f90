module mod_cadas
  !use debug_options
  implicit none

  !
  !     -*-fortran-*-
  !
  ! common /cadas/  year,yeardf,iclass,cdatopt,ptesa,pnesa,pnbs,pnhs,pnzsa,pcoef,iyearh,&
  !     iyearz,useridh,useridz,tmpuserid,adas_iz_rate_mult,adas_rec_rate_mult
  !
  ! save /cadas/
  integer,public :: iclass,cdatopt,iyearz
  ! jdemod -  set the default value of iyearh to -1 to aid in dynamic parameter detection
  ! in LIM. Since iyearh is a required input in both LIM and DIVIMP setting this
  ! value to -1 at this point is not an issue. It may become a problem if the
  ! the code is converted to fully optional input in the future. 
  integer,public :: iyearh = -1
  
  real,public,allocatable :: ptesa(:),pnesa(:),pnbs(:),pnhs(:)
  real,public,allocatable :: pnzsa(:,:),pcoef(:,:)
  real,public :: adas_iz_rate_mult,adas_rec_rate_mult
  
  character,public :: year*2,yeardf*2,useridh*80,useridz*80,tmpuserid*80

  public :: allocate_mod_cadas,deallocate_mod_cadas

contains

  subroutine allocate_mod_cadas(maxpts,maxizs)
    !use mod_params
    use allocate_arrays
    implicit none
    integer :: maxpts
    integer :: maxizs
    integer :: ierr

    !call pr_trace('mod_cadas','ALLOCATE')

    call allocate_array(ptesa,maxpts,'ptesa',ierr)
    call allocate_array(pnesa,maxpts,'pnesa',ierr)
    call allocate_array(pnbs,maxpts,'pnbs',ierr)
    call allocate_array(pnhs,maxpts,'pnhs',ierr)
    call allocate_array(pnzsa,1,maxpts,0,maxizs,'pnzsa',ierr)
    call allocate_array(pcoef,maxpts,maxizs,'pcoef',ierr)

  end subroutine allocate_mod_cadas


  subroutine deallocate_mod_cadas
    implicit none

    !call pr_trace('mod_cadas','DEALLOCATE')

    if (allocated(ptesa)) deallocate(ptesa)
    if (allocated(pnesa)) deallocate(pnesa)
    if (allocated(pnbs)) deallocate(pnbs)
    if (allocated(pnhs)) deallocate(pnhs)
    if (allocated(pnzsa)) deallocate(pnzsa)
    if (allocated(pcoef)) deallocate(pcoef)

  end subroutine deallocate_mod_cadas

end module mod_cadas
