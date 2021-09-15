module mod_rtheta

  use mod_params


!C
!C  THIS FILE CONTAINS ALL THE DECLARATIONS FOR THE R,THETA
!C  PLOT IMPLEMENTATION.
!C
!
!      COMMON /RTHETA/ NRS,NAS,NSS,RTMESH,RTXCOORD,RTYCOORD,
!     >                SOUTS,TOUTS,ROUTS,SWIDS,TWIDS,RWIDS,
!     >                RTSXY,PSIMIN,PSIMAX,PSIHALF,RTDONE,PSISTEP,
!     >                RS,TS,RVPOS,TVPOS,NRDIV,SOUTS1,PSIRWID,
!     >                TOUTS1,TWIDS1,RWIDM,TWIDM,RTAXY
!C    
!      INTEGER  MAXNAS,MAXNRS,MAXNSS,MAXNRDIV
!      PARAMETER (MAXNAS=100,MAXNRS=50,MAXNSS=240,MAXNRDIV=500)
!      INTEGER  NRS,NAS,NSS,NRDIV
!      REAL     SOUTS(MAXNSS),SWIDS(MAXNSS)
!      REAL     SOUTS1(MAXNSS)
!      REAL     ROUTS(MAXNRS),RWIDS(MAXNRS)
!      REAL     TOUTS(MAXNAS),TWIDS(MAXNAS)
!      REAL     TOUTS1(MAXNAS),TWIDS1(MAXNAS)
!      INTEGER  RTXCOORD(MAXNRS),RTYCOORD(MAXNAS)
!      REAL     RTMESH(MAXNRS,MAXNAS)
!      INTEGER  RTSXY(MAXNRDIV,MAXNSS,2)
!      REAL     RTAXY(MAXNRDIV,MAXNSS,2)
!      REAL     PSIMIN,PSIMAX,PSISTEP,PSIHALF
!      REAL     PSIRWID
!      REAL     RS(MAXNRS),TS(MAXNAS)
!      REAL     RWIDM,TWIDM
!      REAL     RVPOS,TVPOS 
!      LOGICAL  RTDONE
!

  implicit none
  private

      INTEGER,public::  MAXNAS,MAXNRS,MAXNSS,MAXNRDIV
      PARAMETER (MAXNAS=100,MAXNRS=50,MAXNSS=240,MAXNRDIV=500)
      INTEGER,public::  NRS,NAS,NSS,NRDIV
      REAL,public::     SOUTS(MAXNSS),SWIDS(MAXNSS)
      REAL,public::     SOUTS1(MAXNSS)
      REAL,public::     ROUTS(MAXNRS),RWIDS(MAXNRS)
      REAL,public::     TOUTS(MAXNAS),TWIDS(MAXNAS)
      REAL,public::     TOUTS1(MAXNAS),TWIDS1(MAXNAS)
      INTEGER,public::  RTXCOORD(MAXNRS),RTYCOORD(MAXNAS)
      REAL,public::     RTMESH(MAXNRS,MAXNAS)
      INTEGER,public::  RTSXY(MAXNRDIV,MAXNSS,2)
      REAL,public::     RTAXY(MAXNRDIV,MAXNSS,2)
      REAL,public::     PSIMIN,PSIMAX,PSISTEP,PSIHALF
      REAL,public::     PSIRWID
      REAL,public::     RS(MAXNRS),TS(MAXNAS)
      REAL,public::     RWIDM,TWIDM
      REAL,public::     RVPOS,TVPOS 
      LOGICAL,public::  RTDONE
 
  public :: allocate_mod_rtheta, deallocate_mod_rtheta


contains

  subroutine allocate_mod_rtheta
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_rtheta


  subroutine deallocate_mod_rtheta
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_rtheta



end module mod_rtheta
