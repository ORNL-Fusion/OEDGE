module mod_div1
  use debug_options
  implicit none

  !
  !     this common block contains all of the divimp local variables
  !     many of these are needed in the transport routines. those that are
  !     not will be moved back to local variable status within the div module.
  !
  !     -*-fortran-*-
  !        jdemod - moved vel to particle_specs
  !     >  ssef,yeff,vel,adjust,dcross,quant,
  !        jdemod - moved temi to particle_specs
  !     >  porm,yldtot,yldmax,temi,riz,avvpos,avtpos,avnpos,
  ! common /div1c/ zioniz_tor,partim,timusd,statim,vfluid,twalln,tdep,rfail,tfail,walltotn,&
  !     tmpsrc,tmpdep,tmpmult,tneut,tatiz,twall,rdifft,smax,slast,avatiz,ssef,yeff,&
  !     adjust,dcross,quant,avxpos,avypos,ttmax,tcent,tbyond,tcut,ratiz,rneut,k,rwalln,&
  !     rcent,rtmax,rdep,rwall,avkpos,avspos,avsmax,porm,yldtot,yldmax,riz,avvpos,avtpos,&
  !     avnpos,sputy,energy,rneut1,ryield,ythtot,rmain,tmain,sheath_fraction,fytot,ftot,&
  !     tbelow,spunew,rexit,texit,neut2d_fytot
  
  ! save /div1c/
  !
  real,public :: zioniz_tor
  real,public :: partim,timusd
  real,public :: statim,vfluid,twalln,tdep,rfail,tfail
  real,public :: walltotn,tmpsrc,tmpdep,tmpmult
  real,public :: tneut,tatiz,twall,rdifft,smax,slast
  real,public,allocatable :: avatiz(:)
  !      real      ssef,yeff,vel,adjust,dcross(4),quant
  real,public :: ssef,yeff,adjust,quant
  real,public,allocatable :: dcross(:)
  real,public :: ttmax,tcent,tbyond,tcut,ratiz,rneut
  real,public,allocatable :: avxpos(:),avypos(:)
  real,public :: k,rwalln,rcent,rtmax,rdep,rwall,avkpos,avspos,avsmax
  !      real      porm,yldtot,yldmax,temi,riz,avvpos,avtpos(2),avnpos(2)
  real,public :: porm,yldtot,yldmax,riz,avvpos
  real,public,allocatable :: avtpos(:),avnpos(:)
  real,public :: sputy,energy,rneut1,ryield,ythtot,rmain,tmain
  real,public :: sheath_fraction
  real,public :: fytot,ftot,tbelow,spunew,rexit,texit
  real,public :: neut2d_fytot
  ! moved to local variables and subroutine arguments
  !
  !common /div1ca/ spara,dspara,vpara,dvpara
  !save /div1ca/
  
  !real spara,dspara,vpara,dvpara

  public :: allocate_mod_div1,deallocate_mod_div1

contains

  subroutine allocate_mod_div1
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div1','ALLOCATE')

    call allocate_array(avatiz,2,'avatiz',ierr)
    call allocate_array(dcross,4,'dcross',ierr)
    call allocate_array(avxpos,2,'avxpos',ierr)
    call allocate_array(avypos,2,'avypos',ierr)
    call allocate_array(avtpos,2,'avtpos',ierr)
    call allocate_array(avnpos,2,'avnpos',ierr)

  end subroutine allocate_mod_div1


  subroutine deallocate_mod_div1
    implicit none

    call pr_trace('mod_div1','DEALLOCATE')

    if (allocated(avatiz)) deallocate(avatiz)
    if (allocated(dcross)) deallocate(dcross)
    if (allocated(avxpos)) deallocate(avxpos)
    if (allocated(avypos)) deallocate(avypos)
    if (allocated(avtpos)) deallocate(avtpos)
    if (allocated(avnpos)) deallocate(avnpos)

  end subroutine deallocate_mod_div1

end module mod_div1