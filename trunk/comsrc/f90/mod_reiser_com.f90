module mod_reiser_com
  use debug_options
  implicit none

  !
  !     variables required for the reiser transport formulation
  !
  !     -*-fortran-*-
  ! common /reiser/ kvhgs,cioptr,lambda1,alpha,alphai,cc1,cc2,qtim2,pi2sqrt,massfactor,&
  !     exp_coeff,fthi,ffi,fvbg,diff,fcell,velavg,aswitch,sk11,sk12,sk13,sd11,sd12,sd13,&
  !     velover,switchb,switchz,switchv,linearpeak,ltb,ltbmin,lvb,lvbmin,coulomb_log,&
  !     combine
  !
  ! save /reiser/
  !
  integer,public :: cioptr,aswitch,sk11,sk12,sk13,sd11,sd12,sd13,velover,switchb,switchz,&
       switchv,linearpeak,combine
  
  real,public :: coulomb_log
  !
  real,public :: qtim2,pi2sqrt,massfactor
  real,public,allocatable :: kvhgs(:,:),lambda1(:),alpha(:,:),alphai(:,:),cc1(:,:),&
       cc2(:,:),exp_coeff(:),ltb(:,:),ltbmin(:,:),lvb(:,:),lvbmin(:,:)
  
  
  
  
  
  
  
  
  !
  !      double precision fthi(maxnks,maxnrs,maxizs),
  !     >       ffi(maxnks,maxnrs,maxizs),fvbg(maxnks,maxnrs,maxizs),
  !     >       fcell(maxnks,maxnrs,maxizs),diff(maxnks,maxnrs,maxizs),
  !     >       velavg(maxnks,maxnrs,maxizs)
  !
  !      common /forceswitches/ aswitch,sff,sfe,sfeg,sfig,sdp,
  !     >                       sk11,sk12,sk13,sd11,sd12,sd13,
  !     >                       datafile,compare,sfftherm,sbeta,
  !     >                       velover,splinesw
  !      integer aswitch,sff,sfe,sfeg,sfig,sdp,
  !     >        sk11,sk12,sk13,sd11,sd12,sd13,
  !     >        datafile,compare,sfftherm,sbeta,velover,
  !     >        splinesw
  !
  real,public,allocatable :: fthi(:,:,:),ffi(:,:,:),fvbg(:,:,:),fcell(:,:,:),diff(:,:,:),&
       velavg(:,:,:)

  public :: allocate_mod_reiser_com,deallocate_mod_reiser_com

contains

  subroutine allocate_mod_reiser_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_reiser_com','ALLOCATE')

    call allocate_array(kvhgs,maxnks,maxnrs,'kvhgs',ierr)
    call allocate_array(lambda1,maxizs,'lambda1',ierr)
    call allocate_array(alpha,maxnks,maxnrs,'alpha',ierr)
    call allocate_array(alphai,maxnks,maxnrs,'alphai',ierr)
    call allocate_array(cc1,maxnks,maxnrs,'cc1',ierr)
    call allocate_array(cc2,maxnks,maxnrs,'cc2',ierr)
    call allocate_array(exp_coeff,8,'exp_coeff',ierr)
    call allocate_array(ltb,maxnks,maxnrs,'ltb',ierr)
    call allocate_array(ltbmin,maxnks,maxnrs,'ltbmin',ierr)
    call allocate_array(lvb,maxnks,maxnrs,'lvb',ierr)
    call allocate_array(lvbmin,maxnks,maxnrs,'lvbmin',ierr)
    call allocate_array(fthi,maxnks,maxnrs,maxizs,'fthi',ierr)
    call allocate_array(ffi,maxnks,maxnrs,maxizs,'ffi',ierr)
    call allocate_array(fvbg,maxnks,maxnrs,maxizs,'fvbg',ierr)
    call allocate_array(fcell,maxnks,maxnrs,maxizs,'fcell',ierr)
    call allocate_array(diff,maxnks,maxnrs,maxizs,'diff',ierr)
    call allocate_array(velavg,maxnks,maxnrs,maxizs,'velavg',ierr)

  end subroutine allocate_mod_reiser_com


  subroutine deallocate_mod_reiser_com
    implicit none

    call pr_trace('mod_reiser_com','DEALLOCATE')

    if (allocated(kvhgs)) deallocate(kvhgs)
    if (allocated(lambda1)) deallocate(lambda1)
    if (allocated(alpha)) deallocate(alpha)
    if (allocated(alphai)) deallocate(alphai)
    if (allocated(cc1)) deallocate(cc1)
    if (allocated(cc2)) deallocate(cc2)
    if (allocated(exp_coeff)) deallocate(exp_coeff)
    if (allocated(ltb)) deallocate(ltb)
    if (allocated(ltbmin)) deallocate(ltbmin)
    if (allocated(lvb)) deallocate(lvb)
    if (allocated(lvbmin)) deallocate(lvbmin)
    if (allocated(fthi)) deallocate(fthi)
    if (allocated(ffi)) deallocate(ffi)
    if (allocated(fvbg)) deallocate(fvbg)
    if (allocated(fcell)) deallocate(fcell)
    if (allocated(diff)) deallocate(diff)
    if (allocated(velavg)) deallocate(velavg)

  end subroutine deallocate_mod_reiser_com

end module mod_reiser_com