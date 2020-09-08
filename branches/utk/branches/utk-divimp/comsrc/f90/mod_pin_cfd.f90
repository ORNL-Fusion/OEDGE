module mod_pin_cfd
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /pin_cfd/ d2ndr, d2ndr_sm, d2ndr_last,pinmom, pinmom_sm, pinmom_last,pinmom_rec,&
  !      pinmom_cx,pinion_sm, pinion_last,pinatom_sm,  pinatom_last,pinmol_sm,  pinmol_last,&
  !     pinena_sm,  pinena_last,pinqi_sm,  pinqi_last,pinqi_1,  pinqi_2, pinqi_3,&
  !     pinqi_4,  pinqi_5, pinqi_6,pinqe_1,  pinqe_2,pinqe_sm,  pinqe_last,pinrec_sm,&
  !      pinrec_last,fract4_new, fract2_new, fract3_new, fract1_new,relax4,     relax1,&
  !      relax2, relax3,fract4_old, fract2_old, fract3_old, fract1_old,f_recomb, oldknbs,&
  !      oldktibs,oldktebs, oldkvhs, old_ne, old_ti,old_te, old_v, old_mach, zhi_sm,fract1_a,&
  !      fract1_b,fract2_a, fract2_b,fract3_a, fract3_b,fract4_a, fract4_b,a_sm_s,&
  !      a_sm_r, a_relax, sum3a, sum3b, sum4a, sum4b,fract3, delta3, fract2_pin, fract3_pin,&
  !      delta2,sum2a, sum2b, sum2c, targ2a, targ2b, targ2c,targ3a, targ3b, targ3c,&
  !      delta, core4a, core4b, core4c,fract4, delta4, fract4_pin, targ4a, targ4b, targ4c,&
  !     core3a, core3b, core3c, rel_min2, rel_min3, rel_min,sum1a, sum1b, sum1c, delta1,&
  !      targ1a, targ1b, targ1c,fract1_pin, err_fl_a, err_fl_b
  !
  ! save /pin_cfd/
  real,public,allocatable :: d2ndr(:,:),d2ndr_sm(:,:)
  real,public,allocatable :: d2ndr_last(:,:),pinmom(:,:)
  real,public,allocatable :: pinmom_sm(:,:),pinmom_last(:,:),pinmom_rec(:,:),pinmom_cx(:,:)
  real,public,allocatable :: pinion_sm(:,:),pinion_last(:,:)
  real,public,allocatable :: pinatom_sm(:,:),pinatom_last(:,:)
  real,public,allocatable :: pinmol_sm(:,:),pinmol_last(:,:)
  real,public,allocatable :: pinena_sm(:,:),pinena_last(:,:)
  real,public,allocatable :: pinqi_sm(:,:),pinqi_last(:,:)
  real,public,allocatable :: pinqi_1(:,:),pinqi_2(:,:)
  real,public,allocatable :: pinqi_3(:,:),pinqi_4(:,:)
  real,public,allocatable :: pinqi_5(:,:),pinqi_6(:,:)
  real,public,allocatable :: pinqe_1(:,:),pinqe_2(:,:)
  real,public,allocatable :: pinqe_sm(:,:),pinqe_last(:,:)
  real,public,allocatable :: pinrec_sm(:,:),pinrec_last(:,:)
  real,public,allocatable :: fract4_new(:)
  real,public,allocatable :: fract2_new(:),fract3_new(:),fract1_new(:)
  real,public,allocatable :: relax4(:)
  real,public,allocatable :: relax1(:),relax2(:),relax3(:)
  real,public,allocatable :: fract4_old(:)
  real,public,allocatable :: fract2_old(:),fract3_old(:),fract1_old(:)
  real,public,allocatable :: f_recomb(:)
  real,public,allocatable :: oldknbs(:,:)
  real,public,allocatable :: oldktibs(:,:)
  real,public,allocatable :: oldktebs(:,:)
  real,public,allocatable :: oldkvhs(:,:)
  real,public,allocatable :: old_ne(:,:)
  real,public,allocatable :: old_ti(:,:)
  real,public,allocatable :: old_te(:,:)
  real,public,allocatable :: old_v(:,:)
  real,public,allocatable :: old_mach(:,:)
  real,public,allocatable :: zhi_sm(:,:)
  real,public :: fract1_a,fract1_b,fract2_a,fract2_b,fract3_a,fract3_b
  real,public :: fract4_a,fract4_b
  real,public :: a_sm_s,a_sm_r,a_relax,sum3a,sum3b,sum3c,fract3,delta3,fract2_pin,fract3_pin,&
       delta2,sum2a,sum2b,sum2c,targ2a,targ2b,targ2c,targ3a,targ3b,targ3c,delta,&
       core3a,core3b,core3c,rel_min2,rel_min3,rel_min,sum1a,sum1b,sum1c,delta1,targ1a,&
       targ1b,targ1c,fract1_pin,sum4a,sum4b,sum4c,fract4,delta4,fract4_pin,targ4a,targ4b,&
       targ4c,core4a,core4b,core4c
  
  
  
  
  integer,public :: err_fl_a,err_fl_b

  public :: allocate_mod_pin_cfd,deallocate_mod_pin_cfd

contains

  subroutine allocate_mod_pin_cfd
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_pin_cfd','ALLOCATE')

    call allocate_array(d2ndr,maxnks,maxnrs,'d2ndr',ierr)
    call allocate_array(d2ndr_sm,maxnks,maxnrs,'d2ndr_sm',ierr)
    call allocate_array(d2ndr_last,maxnks,maxnrs,'d2ndr_last',ierr)
    call allocate_array(pinmom,maxnks,maxnrs,'pinmom',ierr)
    call allocate_array(pinmom_sm,maxnks,maxnrs,'pinmom_sm',ierr)
    call allocate_array(pinmom_last,maxnks,maxnrs,'pinmom_last',ierr)
    call allocate_array(pinmom_rec,maxnks,maxnrs,'pinmom_rec',ierr)
    call allocate_array(pinmom_cx,maxnks,maxnrs,'pinmom_cx',ierr)
    call allocate_array(pinion_sm,maxnks,maxnrs,'pinion_sm',ierr)
    call allocate_array(pinion_last,maxnks,maxnrs,'pinion_last',ierr)
    call allocate_array(pinatom_sm,maxnks,maxnrs,'pinatom_sm',ierr)
    call allocate_array(pinatom_last,maxnks,maxnrs,'pinatom_last',ierr)
    call allocate_array(pinmol_sm,maxnks,maxnrs,'pinmol_sm',ierr)
    call allocate_array(pinmol_last,maxnks,maxnrs,'pinmol_last',ierr)
    call allocate_array(pinena_sm,maxnks,maxnrs,'pinena_sm',ierr)
    call allocate_array(pinena_last,maxnks,maxnrs,'pinena_last',ierr)
    call allocate_array(pinqi_sm,maxnks,maxnrs,'pinqi_sm',ierr)
    call allocate_array(pinqi_last,maxnks,maxnrs,'pinqi_last',ierr)
    call allocate_array(pinqi_1,maxnks,maxnrs,'pinqi_1',ierr)
    call allocate_array(pinqi_2,maxnks,maxnrs,'pinqi_2',ierr)
    call allocate_array(pinqi_3,maxnks,maxnrs,'pinqi_3',ierr)
    call allocate_array(pinqi_4,maxnks,maxnrs,'pinqi_4',ierr)
    call allocate_array(pinqi_5,maxnks,maxnrs,'pinqi_5',ierr)
    call allocate_array(pinqi_6,maxnks,maxnrs,'pinqi_6',ierr)
    call allocate_array(pinqe_1,maxnks,maxnrs,'pinqe_1',ierr)
    call allocate_array(pinqe_2,maxnks,maxnrs,'pinqe_2',ierr)
    call allocate_array(pinqe_sm,maxnks,maxnrs,'pinqe_sm',ierr)
    call allocate_array(pinqe_last,maxnks,maxnrs,'pinqe_last',ierr)
    call allocate_array(pinrec_sm,maxnks,maxnrs,'pinrec_sm',ierr)
    call allocate_array(pinrec_last,maxnks,maxnrs,'pinrec_last',ierr)
    call allocate_array(fract4_new,maxnrs,'fract4_new',ierr)
    call allocate_array(fract2_new,maxnrs,'fract2_new',ierr)
    call allocate_array(fract3_new,maxnrs,'fract3_new',ierr)
    call allocate_array(fract1_new,maxnrs,'fract1_new',ierr)
    call allocate_array(relax4,maxnrs,'relax4',ierr)
    call allocate_array(relax1,maxnrs,'relax1',ierr)
    call allocate_array(relax2,maxnrs,'relax2',ierr)
    call allocate_array(relax3,maxnrs,'relax3',ierr)
    call allocate_array(fract4_old,maxnrs,'fract4_old',ierr)
    call allocate_array(fract2_old,maxnrs,'fract2_old',ierr)
    call allocate_array(fract3_old,maxnrs,'fract3_old',ierr)
    call allocate_array(fract1_old,maxnrs,'fract1_old',ierr)
    call allocate_array(f_recomb,maxnrs,'f_recomb',ierr)
    call allocate_array(oldknbs,maxnks,maxnrs,'oldknbs',ierr)
    call allocate_array(oldktibs,maxnks,maxnrs,'oldktibs',ierr)
    call allocate_array(oldktebs,maxnks,maxnrs,'oldktebs',ierr)
    call allocate_array(oldkvhs,maxnks,maxnrs,'oldkvhs',ierr)
    call allocate_array(old_ne,maxnrs,2,'old_ne',ierr)
    call allocate_array(old_ti,maxnrs,2,'old_ti',ierr)
    call allocate_array(old_te,maxnrs,2,'old_te',ierr)
    call allocate_array(old_v,maxnrs,2,'old_v',ierr)
    call allocate_array(old_mach,maxnrs,2,'old_mach',ierr)
    call allocate_array(zhi_sm,4,maxnrs,'zhi_sm',ierr)

  end subroutine allocate_mod_pin_cfd


  subroutine deallocate_mod_pin_cfd
    implicit none

    call pr_trace('mod_pin_cfd','DEALLOCATE')

    if (allocated(d2ndr)) deallocate(d2ndr)
    if (allocated(d2ndr_sm)) deallocate(d2ndr_sm)
    if (allocated(d2ndr_last)) deallocate(d2ndr_last)
    if (allocated(pinmom)) deallocate(pinmom)
    if (allocated(pinmom_sm)) deallocate(pinmom_sm)
    if (allocated(pinmom_last)) deallocate(pinmom_last)
    if (allocated(pinmom_rec)) deallocate(pinmom_rec)
    if (allocated(pinmom_cx)) deallocate(pinmom_cx)
    if (allocated(pinion_sm)) deallocate(pinion_sm)
    if (allocated(pinion_last)) deallocate(pinion_last)
    if (allocated(pinatom_sm)) deallocate(pinatom_sm)
    if (allocated(pinatom_last)) deallocate(pinatom_last)
    if (allocated(pinmol_sm)) deallocate(pinmol_sm)
    if (allocated(pinmol_last)) deallocate(pinmol_last)
    if (allocated(pinena_sm)) deallocate(pinena_sm)
    if (allocated(pinena_last)) deallocate(pinena_last)
    if (allocated(pinqi_sm)) deallocate(pinqi_sm)
    if (allocated(pinqi_last)) deallocate(pinqi_last)
    if (allocated(pinqi_1)) deallocate(pinqi_1)
    if (allocated(pinqi_2)) deallocate(pinqi_2)
    if (allocated(pinqi_3)) deallocate(pinqi_3)
    if (allocated(pinqi_4)) deallocate(pinqi_4)
    if (allocated(pinqi_5)) deallocate(pinqi_5)
    if (allocated(pinqi_6)) deallocate(pinqi_6)
    if (allocated(pinqe_1)) deallocate(pinqe_1)
    if (allocated(pinqe_2)) deallocate(pinqe_2)
    if (allocated(pinqe_sm)) deallocate(pinqe_sm)
    if (allocated(pinqe_last)) deallocate(pinqe_last)
    if (allocated(pinrec_sm)) deallocate(pinrec_sm)
    if (allocated(pinrec_last)) deallocate(pinrec_last)
    if (allocated(fract4_new)) deallocate(fract4_new)
    if (allocated(fract2_new)) deallocate(fract2_new)
    if (allocated(fract3_new)) deallocate(fract3_new)
    if (allocated(fract1_new)) deallocate(fract1_new)
    if (allocated(relax4)) deallocate(relax4)
    if (allocated(relax1)) deallocate(relax1)
    if (allocated(relax2)) deallocate(relax2)
    if (allocated(relax3)) deallocate(relax3)
    if (allocated(fract4_old)) deallocate(fract4_old)
    if (allocated(fract2_old)) deallocate(fract2_old)
    if (allocated(fract3_old)) deallocate(fract3_old)
    if (allocated(fract1_old)) deallocate(fract1_old)
    if (allocated(f_recomb)) deallocate(f_recomb)
    if (allocated(oldknbs)) deallocate(oldknbs)
    if (allocated(oldktibs)) deallocate(oldktibs)
    if (allocated(oldktebs)) deallocate(oldktebs)
    if (allocated(oldkvhs)) deallocate(oldkvhs)
    if (allocated(old_ne)) deallocate(old_ne)
    if (allocated(old_ti)) deallocate(old_ti)
    if (allocated(old_te)) deallocate(old_te)
    if (allocated(old_v)) deallocate(old_v)
    if (allocated(old_mach)) deallocate(old_mach)
    if (allocated(zhi_sm)) deallocate(zhi_sm)

  end subroutine deallocate_mod_pin_cfd

end module mod_pin_cfd