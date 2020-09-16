module mod_slcom
  
  use mod_params
  
  !      common /slcom/
  !c
  !     + slopt,
  !     + optdp,
  !     + n2opt,neopt,nieopt,
  !     + n2frac,neng,nieng,
  !     + n2rate,nnrate,nirate,n2cpchs,nncpchs,
  !
  !     + n2_n,n_n,ni_n,
  !     + n2_lim,n_lim,
  !     + n2_x,n2_y,n2_p,n_x,n_y,n_p,n_x_max,n_y_max,n_p_max,
  !     + n2_e,n_e,ni_e,n2_v,n_v,
  !     + n2_t,n_t,
  !     + n_beta,n_psi,
  !
  !     + titl2,lambda,injbint,tag2,mark,
  !     + tloss,wloss,tgloss,lloss,aloss,izloss,tsloss,totvol,totden,
  !     + bsbin,ysbin,izbin,nbin,nsbin,mizs,mimps,matiz,
  !     + svhins,seyins
  !c
  !c ----------------------------------------------------------------------
  !c
  !c declarations:
  !c
  !c output file unit number - set to general debug unit for now
  !c
  !      integer slout
  !      parameter (slout=dbgunit)
  !
  !c
  !c divimp ion profile option:
  !c
  !      integer optdp
  !
  !      integer slopt,n2opt,neopt,nieopt
  !      real    n2frac,neng,nieng
  !      real    n2rate(maxnxs),nnrate(maxnxs),nirate(maxnxs),
  !     +        n2cpchs(maxnxs),nncpchs(maxnxs)
  !c
  !c statistics for n2 break-up:
  !c
  !      real    n2_n,n_n,ni_n
  !      real    n2_lim,n_lim
  !      real    n2_x,n2_y,n2_p,n_x,n_y,n_p,n_x_max,n_y_max,n_p_max
  !      real    n2_e,n_e,ni_e,n2_v,n_v
  !      real    n2_t,n_t
  !      real    n_beta,n_psi
  !c
  !c lim3 variables (for the most part):
  !c
  !      character  titl2*80
  !      real       lambda
  !      integer    injbint(-maxy3d-1:maxy3d+1)
  !c      integer    injbinp(-maxnps:mannps)
  !      real       tag2(maximp)
  !      real       mark
  !c
  !c statisitcs:
  !c
  !      real       tloss(maximp),wloss,tgloss,lloss,aloss,izloss,tsloss
  !      real       totvol,totden(maxizs)
  !c
  !c variables for storing the divimp source distribution:
  !c
  !      real       bsbin(10),ysbin(10),izbin(10)
  !      integer    nbin, nsbin(10,10),mizs,mimps,matiz
  !c
  !c moved from lim3:
  !c
  !      real      svhins(-maxqxs:maxqxs),seyins(-maxqxs:maxqxs,maxizs)
  
  
  implicit none
  private
  
  ! declarations:
  !
  ! output file unit number - set to general debug unit for now
  ! this is set by a call to the initialization routine at the beginning of runlm3
  ! set a default value of 6 which is the baseline output unit number in lim/divimp
  integer,public:: slout = 6
  !parameter (slout=dbgunit)
  
  !
  ! divimp ion profile option:
  !
  integer,public:: optdp
  
  integer,public:: slopt,n2opt,neopt,nieopt
  real,public:: n2frac,neng,nieng
  real,public,allocatable:: n2rate(:),nnrate(:),nirate(:),n2cpchs(:),nncpchs(:)
  !
  ! statistics for n2 break-up:
  !
  real,public:: n2_n,n_n,ni_n
  real,public:: n2_lim,n_lim
  real,public:: n2_x,n2_y,n2_p,n_x,n_y,n_p,n_x_max,n_y_max,n_p_max
  real,public:: n2_e,n_e,ni_e,n2_v,n_v
  real,public:: n2_t,n_t
  real,public:: n_beta,n_psi
  !
  ! lim3 variables (for the most part):
  !
  character,public::  titl2*80
  real,public:: lambda
  integer,public,allocatable:: injbint(:)
  !      integer,public::    injbinp(-maxnps:mannps)
  real,public,allocatable:: tag2(:)
  real,public:: mark
  !
  ! statisitcs:
  !
  real,public:: wloss,tgloss,lloss,aloss,izloss,tsloss
  real,public,allocatable:: tloss(:)
  real,public:: totvol
  real,public,allocatable:: totden(:)
  !
  ! variables for storing the divimp source distribution:
  !
  !      real,public::       bsbin(10),ysbin(10),izbin(10)
  !
  real,public,allocatable:: izbin(:)
  
  
  real,public::  YSBIN(10) = (/ -1.0E+06, -3.2317, -2.7438, -1.9309, &
                    -1.1181, -0.3259,  0.5338,  1.7560,&
                     3.7696, 1.0E+06 /)
  real,public::  BSBIN(10) = (/  -3.1231, -3.0029, -2.5149, -1.7021,&
                     -0.8893, -0.0971,  0.7625,  1.9847,&
                      3.9985, 1.0E+06 /)
  
  integer,public:: nbin = 10
  
  
  
  integer,public:: mizs,mimps,matiz
  integer,public,allocatable:: nsbin(:,:)
  !
  ! moved from lim3:
  !
  real,public,allocatable:: svhins(:),seyins(:,:)
  
  public :: allocate_mod_slcom, deallocate_mod_slcom,set_sl_outunit
  
  
contains
  
  subroutine set_sl_outunit(sloutunit)
    implicit none
    integer :: sloutunit
    slout = sloutunit
  end subroutine set_sl_outunit
  
  
  subroutine allocate_mod_slcom
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(dtev  ,maxnxs,'dtev',ierr)


    call allocate_array(n2rate,maxnxs,'n2rate',ierr)
    call allocate_array(nnrate,maxnxs,'nnrate',ierr)
    call allocate_array(nirate,maxnxs,'nirate',ierr)
    call allocate_array(n2cpchs,maxnxs,'n2cpchs',ierr)
    call allocate_array(nncpchs,maxnxs,'nncpchs',ierr)
    !
    ! jdemod - this array is used in Steve's code from -nys,nys
    !        - it should not be declared based on maxy3d
    !
    !call allocate_array(injbint,-maxy3d-1,'injbint',maxy3d+1,ierr)
    call allocate_array(injbint,-maxnys-1,'injbint',maxnys+1,ierr)
    call allocate_array(tag2,maximp,'tag2',ierr)
    call allocate_array(tloss,maximp,'tloss',ierr)
    call allocate_array(totden,maxizs,'totden',ierr)
    call allocate_array(izbin,10,'izbin',ierr)
    call allocate_array(nsbin,10,10,'nsbin',ierr)
    call allocate_array(svhins,-maxqxs,'svhins',maxqxs,ierr)
    call allocate_array(seyins,-maxqxs,maxqxs,1,maxizs,'seyins',ierr)

  end subroutine allocate_mod_slcom
  
  
  subroutine deallocate_mod_slcom
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

    if (allocated(n2rate)) deallocate(n2rate)
    if (allocated(nnrate)) deallocate(nnrate)
    if (allocated(nirate)) deallocate(nirate)
    if (allocated(n2cpchs)) deallocate(n2cpchs)
    if (allocated(nncpchs)) deallocate(nncpchs)
    if (allocated(injbint)) deallocate(injbint)
    if (allocated(tag2)) deallocate(tag2)
    if (allocated(tloss)) deallocate(tloss)
    if (allocated(totden)) deallocate(totden)
    if (allocated(izbin)) deallocate(izbin)
    if (allocated(nsbin)) deallocate(nsbin)
    if (allocated(svhins)) deallocate(svhins)
    if (allocated(seyins)) deallocate(seyins)

  end subroutine deallocate_mod_slcom
  
  
  
end module mod_slcom
