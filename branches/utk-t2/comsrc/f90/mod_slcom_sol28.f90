module mod_slcom_sol28
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  
  ! common /sol28com/ik0,ik1,ik2,ik3,ik4,ik4i,ik5,ik6,id0,id6,ik12,ik25,ik34,ik45,machcountl2,&
  !     loadsources,imaginaryl2,scalerec,calcmom1,fitline,detach1,detach2,new,&
  !     s28output,s,s12,s25,s34,s45,deltas,n,te,ti,v,teval,timul,isat,p0,p1,p2,p3,p4,p5,&
  !     p6,p12,pmatch,lastp0,lastp5,lastp6,addmp,holdn2,holdv2,lastknbs,lastv1,machno,mach1,&
  !     mach2,machadjustl2,ion01,ion04,rec01,ion05,ion06,rec03,rec04,rec05,rec35,rec06,&
  !     rec36,flx0,flx1,flx2,flx4,flx5,flx6,cfp01,cfp04,cfp05,cfp06,techeat
  integer,public :: ik0,ik1,ik2,ik3,ik4,ik4i,ik5,ik6,id0,id6,ik12,ik25,ik34,ik45,machcountl2
  logical,public :: loadsources,imaginaryl2,scalerec,calcmom1,fitline,detach1,detach2,&
       new,s28output
  
  
  real,public :: s12,s25,s34,s45,deltas,teval,timul,p0,p1,p2,p3,p4,p5,p6,p12,pmatch,&
       lastp0,lastp5,lastp6,holdn2,holdv2,lastv1,machno,mach1,mach2,machadjustl2,ion01,&
       ion04,rec01,ion05,ion06,rec03,rec04,rec05,rec35,rec06,rec36,flx0,flx1,flx2,flx4,&
       flx5,flx6,cfp01,cfp04,cfp05,cfp06,techeat
  real,public,allocatable :: s(:),n(:),te(:),ti(:),v(:),isat(:),addmp(:),lastknbs(:)

  public :: allocate_mod_slcom_sol28,deallocate_mod_slcom_sol28

contains

  subroutine allocate_mod_slcom_sol28
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_slcom_sol28','ALLOCATE')

    call allocate_array(s,0,'s',6,ierr)
    call allocate_array(n,0,'n',6,ierr)
    call allocate_array(te,0,'te',6,ierr)
    call allocate_array(ti,0,'ti',6,ierr)
    call allocate_array(v,0,'v',6,ierr)
    call allocate_array(isat,0,'isat',6,ierr)
    call allocate_array(addmp,maxnks,'addmp',ierr)
    call allocate_array(lastknbs,maxnks,'lastknbs',ierr)

  end subroutine allocate_mod_slcom_sol28


  subroutine deallocate_mod_slcom_sol28
    implicit none

    call pr_trace('mod_slcom_sol28','DEALLOCATE')

    if (allocated(s)) deallocate(s)
    if (allocated(n)) deallocate(n)
    if (allocated(te)) deallocate(te)
    if (allocated(ti)) deallocate(ti)
    if (allocated(v)) deallocate(v)
    if (allocated(isat)) deallocate(isat)
    if (allocated(addmp)) deallocate(addmp)
    if (allocated(lastknbs)) deallocate(lastknbs)

  end subroutine deallocate_mod_slcom_sol28

end module mod_slcom_sol28