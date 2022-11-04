module mod_solswitch
  use debug_options
  implicit none

  !
  !     this contains the common block with the sol option switches
  !
  !     -*-fortran-*-
  integer,public :: maxopts
  !
  parameter (maxopts=45)    ! maxopts needs to be higher than the maximum switch parameter
  !
  integer,public :: swcond,swconv,swprad,swphelp,swpei,swpcx,swmach,swvisc1,swnmom,&
       swion,swe2d,swpow,swgperp,swmajr,swcore,swsmooth,swioni,swerror,swrecom,swionp,swpowp,&
       swqidatiz,swqidmliz,swqidcx,swqidrec,swdetach,swgperpp,swextra,swppion,swppelec,&
       swppress,swqperpe,swqperpi,swepow,swipow
  !
  parameter (swcond=1,swconv=2,swprad=3,swphelp=4,swpei=5,swpcx=6,swmach=7,swvisc1=8,&
       swnmom=9,swion=10,swe2d=11,swpow=12,swgperp=13,swmajr=14,swcore=15,swsmooth=16,&
       swioni=17,swerror=18,swrecom=19,swionp=20,swpowp=21,swqidatiz=22,swqidmliz=23,swqidcx=24,&
       swqidrec=25,swdetach=26,swgperpp=27,swextra=28,swppion=29,swppelec=30,&
       swppress=31,swqperpe=32,swqperpi=33,swepow=34,swipow=35)
  ! common /solswitch/ switch,deflist,ndef
  !
  ! save /solswitch/
  ! common /activesw/ actswcond,actswconv,actswprad,actswphelp,actswpei,actswpcx,actswmach,&
  !     actswvisc1,actswnmom,actswion,actswe2d,actswpow,actswgperp,actswmajr,actswcore,&
  !     actswsmooth,actswioni,actswerror,actswrecom,actswqidatiz,actswqidmliz,actswqidcx,&
  !     actswqidrec,actswdetach,actswppion,actswppelec,actswppress
  !
  ! save /activesw/
  !
  real,public,allocatable :: switch(:),deflist(:,:)
  !
  real,public :: actswcond,actswconv,actswprad,actswphelp,actswpei,actswpcx,actswmach,&
       actswvisc1,actswnmom,actswion,actswe2d,actswpow,actswgperp,actswmajr,actswcore,&
       actswsmooth,actswioni,actswerror,actswrecom,actswqidatiz,actswqidmliz,actswqidcx,&
       actswqidrec,actswdetach,actswppion,actswppelec,actswppress,actswqperpe,actswqperpi,&
       actswepow,actswipow
  
  integer,public :: ndef

  public :: allocate_mod_solswitch,deallocate_mod_solswitch,allocate_mod_solswitch_input

contains

  subroutine allocate_mod_solswitch
    !use mod_params
    use mod_solparams
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_solswitch','ALLOCATE')

    call allocate_array(switch,maxopts,'switch',ierr)
    call allocate_array(deflist,mxspts,3,'deflist',ierr)

  end subroutine allocate_mod_solswitch


  subroutine deallocate_mod_solswitch
    implicit none

    call pr_trace('mod_solswitch','DEALLOCATE')

    if (allocated(switch)) deallocate(switch)
    if (allocated(deflist)) deallocate(deflist)

  end subroutine deallocate_mod_solswitch


  subroutine allocate_mod_solswitch_input
    !use mod_params
    use mod_solparams
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_solswitch','ALLOCATE INPUT')

    ! moved to primary allocation since all of these are allocated before input anyway
    ! - this makes it easier when the allocation code is called in LIM
    !call allocate_array(deflist,mxspts,3,'deflist',ierr)

  end subroutine allocate_mod_solswitch_input

end module mod_solswitch
