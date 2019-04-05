module mod_cedge2d
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /edge2d/ e2dnbs, e2dvhs,e2dtibs,e2dtebs,e2des,e2dion,e2datom,e2drec,e2dz0,&
  !     e2diz0,cre2d,e2dtarg,e2dmach,e2dbvel,e2dflux,dthetag,drho,e2dnzs,e2dpowls,e2dlines,&
  !     cre2dizs,e2dtargopt,e2dareas,e2dvro,e2dgdown,e2dgpara,e2dcxrec,e2dhrec,e2dpedown,&
  !     e2dpepara,e2dvts,e2dpidown,e2dpipara,e2dnes,e2dvzs,fc_target_calc_option,fc_v_calc_opt,&
  !     fc_te_calc_opt,fc_ti_calc_opt,fc_ne_calc_opt,fc_v_interp_opt
  !
  ! slmod end
  !...  need target attribute in order to assign pointer in plot 987:
  ! save /edge2d/
  !
  !      real  e2dnbs(maxnks,maxnrs),e2dvhs(maxnks,maxnrs),
  ! slmod begin
  !
  real, target ,public,allocatable :: e2dnbs(:,:),e2dvhs(:,:),e2dtibs(:,:),e2dtebs(:,:),&
       e2des(:,:),e2dion(:,:),e2dtarg(:,:,:),e2dmach(:,:),e2dbvel(:,:),e2datom(:,:),&
       e2dflux(:,:),e2drec(:,:),e2dcxrec(:,:),e2dhrec(:,:),dthetag(:,:),drho(:,:),e2dnzs(:,:,:),&
       e2dpowls(:,:,:),e2dlines(:,:,:),e2dareas(:,:),e2dvro(:,:),e2dgdown(:,:),&
       e2dgpara(:,:),e2dpedown(:,:),e2dpepara(:,:),e2dpidown(:,:),e2dpipara(:,:),e2dz0(:,:),&
       e2diz0(:,:),e2dnes(:,:),e2dvzs(:,:,:),e2dvts(:,:,:)
  !
  integer,public :: cre2d,cre2dizs,e2dtargopt
  
  
  
  
  
  !
  integer,public :: fc_target_calc_option,fc_v_calc_opt,fc_te_calc_opt,fc_ti_calc_opt,&
       fc_ne_calc_opt,fc_v_interp_opt

  public :: allocate_mod_cedge2d,deallocate_mod_cedge2d

contains

  subroutine allocate_mod_cedge2d
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cedge2d','ALLOCATE')

    call allocate_array(e2dnbs,maxnks,maxnrs,'e2dnbs',ierr)
    call allocate_array(e2dvhs,maxnks,maxnrs,'e2dvhs',ierr)
    call allocate_array(e2dtibs,maxnks,maxnrs,'e2dtibs',ierr)
    call allocate_array(e2dtebs,maxnks,maxnrs,'e2dtebs',ierr)
    call allocate_array(e2des,maxnks,maxnrs,'e2des',ierr)
    call allocate_array(e2dion,maxnks,maxnrs,'e2dion',ierr)
    call allocate_array(e2dtarg,maxnrs,8,2,'e2dtarg',ierr)
    call allocate_array(e2dmach,0,maxnks+1,1,maxnrs,'e2dmach',ierr)
    call allocate_array(e2dbvel,maxnks+1,maxnrs,'e2dbvel',ierr)
    call allocate_array(e2datom,maxnks,maxnrs,'e2datom',ierr)
    call allocate_array(e2dflux,maxnks+1,maxnrs,'e2dflux',ierr)
    call allocate_array(e2drec,maxnks,maxnrs,'e2drec',ierr)
    call allocate_array(e2dcxrec,maxnks,maxnrs,'e2dcxrec',ierr)
    call allocate_array(e2dhrec,maxnks,maxnrs,'e2dhrec',ierr)
    call allocate_array(dthetag,maxnks,maxnrs,'dthetag',ierr)
    call allocate_array(drho,maxnks,maxnrs,'drho',ierr)
    call allocate_array(e2dnzs,1,maxnks,1,maxnrs,0,maxe2dizs,'e2dnzs',ierr)
    call allocate_array(e2dpowls,1,maxnks,1,maxnrs,0,maxe2dizs,'e2dpowls',ierr)
    call allocate_array(e2dlines,1,maxnks,1,maxnrs,0,maxe2dizs,'e2dlines',ierr)
    call allocate_array(e2dareas,maxnks,maxnrs,'e2dareas',ierr)
    call allocate_array(e2dvro,maxnks,maxnrs,'e2dvro',ierr)
    call allocate_array(e2dgdown,maxnks,maxnrs,'e2dgdown',ierr)
    call allocate_array(e2dgpara,maxnks,maxnrs,'e2dgpara',ierr)
    call allocate_array(e2dpedown,maxnks,maxnrs,'e2dpedown',ierr)
    call allocate_array(e2dpepara,maxnks,maxnrs,'e2dpepara',ierr)
    call allocate_array(e2dpidown,maxnks,maxnrs,'e2dpidown',ierr)
    call allocate_array(e2dpipara,maxnks,maxnrs,'e2dpipara',ierr)
    call allocate_array(e2dz0,maxnks,maxnrs,'e2dz0',ierr)
    call allocate_array(e2diz0,maxnks,maxnrs,'e2diz0',ierr)
    call allocate_array(e2dnes,maxnks,maxnrs,'e2dnes',ierr)
    call allocate_array(e2dvzs,1,maxnks,1,maxnrs,0,maxe2dizs,'e2dvzs',ierr)
    call allocate_array(e2dvts,1,maxnrs,1,2,0,maxe2dizs,'e2dvts',ierr)

  end subroutine allocate_mod_cedge2d


  subroutine deallocate_mod_cedge2d
    implicit none

    call pr_trace('mod_cedge2d','DEALLOCATE')

    if (allocated(e2dnbs)) deallocate(e2dnbs)
    if (allocated(e2dvhs)) deallocate(e2dvhs)
    if (allocated(e2dtibs)) deallocate(e2dtibs)
    if (allocated(e2dtebs)) deallocate(e2dtebs)
    if (allocated(e2des)) deallocate(e2des)
    if (allocated(e2dion)) deallocate(e2dion)
    if (allocated(e2dtarg)) deallocate(e2dtarg)
    if (allocated(e2dmach)) deallocate(e2dmach)
    if (allocated(e2dbvel)) deallocate(e2dbvel)
    if (allocated(e2datom)) deallocate(e2datom)
    if (allocated(e2dflux)) deallocate(e2dflux)
    if (allocated(e2drec)) deallocate(e2drec)
    if (allocated(e2dcxrec)) deallocate(e2dcxrec)
    if (allocated(e2dhrec)) deallocate(e2dhrec)
    if (allocated(dthetag)) deallocate(dthetag)
    if (allocated(drho)) deallocate(drho)
    if (allocated(e2dnzs)) deallocate(e2dnzs)
    if (allocated(e2dpowls)) deallocate(e2dpowls)
    if (allocated(e2dlines)) deallocate(e2dlines)
    if (allocated(e2dareas)) deallocate(e2dareas)
    if (allocated(e2dvro)) deallocate(e2dvro)
    if (allocated(e2dgdown)) deallocate(e2dgdown)
    if (allocated(e2dgpara)) deallocate(e2dgpara)
    if (allocated(e2dpedown)) deallocate(e2dpedown)
    if (allocated(e2dpepara)) deallocate(e2dpepara)
    if (allocated(e2dpidown)) deallocate(e2dpidown)
    if (allocated(e2dpipara)) deallocate(e2dpipara)
    if (allocated(e2dz0)) deallocate(e2dz0)
    if (allocated(e2diz0)) deallocate(e2diz0)
    if (allocated(e2dnes)) deallocate(e2dnes)
    if (allocated(e2dvzs)) deallocate(e2dvzs)
    if (allocated(e2dvts)) deallocate(e2dvts)

  end subroutine deallocate_mod_cedge2d

end module mod_cedge2d