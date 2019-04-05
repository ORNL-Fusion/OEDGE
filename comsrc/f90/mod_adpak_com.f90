module mod_adpak_com
  use mod_params
  use debug_options
  implicit none

  !
  !     dimension variables and arrays for multi-charge model rates
  !
  !     labelrt - header information from 'b2frates' data file
  !     rtnt    - number of intervals in 'b2frates' temperature data
  !     rtnn    - number of intervals in 'b2frates' density data
  !     rtns    - number of species in a 'b2frates' file
  !     rtnsd   - total number of species from all 'b2frates' files
  !     rtza    - atomic charge state
  !     rtzn    - nuclear charge state
  !     rtza2   - atomic charge state squared
  !     rtt     - [ev] temperature data in 'b2frates' table
  !     rtn     - [/m**3] density data in 'b2frates' table
  !     rtlt    - ln(rtt) where rtt[ev] is 'b2frates' temperature data
  !     rtln    - ln(rtn) where rtn[/m**3] is 'b2frates' density data
  !     rtlsa   - ln(rtsa) where rtsa[m**3/s] is 'b2frates'
  !               data for ionization
  !     rtlra   - ln(rtra) where rtra[m**3/s] is 'b2frates'
  !               data for recombination
  !     rtlqa   - ln(rtqa) where rtqa[ev*m**3/s] is 'b2frates'
  !               data for electron energy loss
  !     rtlcx   - ln(rtcx) where rtcx[m**3/s] is 'b2frates'
  !               data for c-x on neutral hydrogen
  !
  !     -*-fortran-*-
  integer,public :: maxrtnsd,maxrtnt,maxrtnn
  !
  parameter (maxrtnsd=maxizs+1,maxrtnt=50,maxrtnn=50)
  ! common /adpak/ rtnt,rtnn,rtns,rtnsd,rtza,rtzn,rtza2,rtt,rtn,rtlt,rtln,rtlsa,rtlra,&
  !     rtlqa,rtlcx,mcfile,labelrt
  !
  ! save /adpak/
  character*120,public :: labelrt(1:12)
  !
  character*120,public :: mcfile
  integer,public :: rtnt,rtnn,rtns,rtnsd
  real,public,allocatable :: rtza(:),rtzn(:),rtza2(:)
  real,public,allocatable :: rtt(:),rtn(:),rtlt(:),rtln(:)
  real,public,allocatable :: rtlsa(:,:,:),rtlra(:,:,:),rtlqa(:,:,:),rtlcx(:,:,:)

  public :: allocate_mod_adpak_com,deallocate_mod_adpak_com

contains

  subroutine allocate_mod_adpak_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_adpak_com','ALLOCATE')

    call allocate_array(rtza,0,'rtza',maxrtnsd-1,ierr)
    call allocate_array(rtzn,0,'rtzn',maxrtnsd-1,ierr)
    call allocate_array(rtza2,0,'rtza2',maxrtnsd-1,ierr)
    call allocate_array(rtt,0,'rtt',maxrtnt,ierr)
    call allocate_array(rtn,0,'rtn',maxrtnn,ierr)
    call allocate_array(rtlt,0,'rtlt',maxrtnt,ierr)
    call allocate_array(rtln,0,'rtln',maxrtnn,ierr)
    call allocate_array(rtlsa,0,maxrtnt,0,maxrtnn,0,maxrtnsd-1,'rtlsa',ierr)
    call allocate_array(rtlra,0,maxrtnt,0,maxrtnn,0,maxrtnsd-1,'rtlra',ierr)
    call allocate_array(rtlqa,0,maxrtnt,0,maxrtnn,0,maxrtnsd-1,'rtlqa',ierr)
    call allocate_array(rtlcx,0,maxrtnt,0,maxrtnn,0,maxrtnsd-1,'rtlcx',ierr)

  end subroutine allocate_mod_adpak_com


  subroutine deallocate_mod_adpak_com
    implicit none

    call pr_trace('mod_adpak_com','DEALLOCATE')

    if (allocated(rtza)) deallocate(rtza)
    if (allocated(rtzn)) deallocate(rtzn)
    if (allocated(rtza2)) deallocate(rtza2)
    if (allocated(rtt)) deallocate(rtt)
    if (allocated(rtn)) deallocate(rtn)
    if (allocated(rtlt)) deallocate(rtlt)
    if (allocated(rtln)) deallocate(rtln)
    if (allocated(rtlsa)) deallocate(rtlsa)
    if (allocated(rtlra)) deallocate(rtlra)
    if (allocated(rtlqa)) deallocate(rtlqa)
    if (allocated(rtlcx)) deallocate(rtlcx)

  end subroutine deallocate_mod_adpak_com

end module mod_adpak_com
