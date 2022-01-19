module mod_comtau

  use mod_params


!c -*-Fortran-*-
!C                                                                               
!      COMMON /COMTAU/ CTEMI,  CX,     CVABS,  CTBIQX, C215A,  C215B,            
!     >                CMIZS,  CIZ,    CIST,   CIAB,                             
!     >                CIOPTA, CIOPTB, CIOPTC, CIOPTD, CIOPTE,                   
!c
!c jdemod - moved cioptj to comtor since it is really a transport option   
!c    >                CIOPTF, CIOPTG, CIOPTH, CIOPTI, CIOPTJ,                  
!     >                CIOPTF, CIOPTG, CIOPTH, CIOPTI,
!     >                CIOPTK, CIOPTL, CIOPTM, CIOPTN,
!     >                CNEUTA, CNEUTB, CNEUTC, CNEUTD, CNEUTE, CNEUTF,
!     >                CEXNEUT,NRFOPT,NVAOPT
!      REAL            CTEMI,CX,CVABS,CTBIQX,CIST,C215A,C215B                    
!      INTEGER         CIOPTA,CIOPTB,CIOPTC,CIOPTD,CIOPTE,CIOPTF                 
!      INTEGER         CIOPTG,CIOPTH,CIOPTI,CNEUTA,CNEUTB,CNEUTC                 
!c      INTEGER         CMIZS,CIZ,CIAB,CNEUTD,CNEUTE,CNEUTF,CIOPTJ               
!      INTEGER         CMIZS,CIZ,CIAB,CNEUTD,CNEUTE,CNEUTF
!      INTEGER         CIOPTK,CIOPTL,CIOPTM,CIOPTN
!      INTEGER         CEXNEUT,NRFOPT,NVAOPT

  
  
  implicit none
  private


      REAL,public::            CTEMI,CX,CVABS,CTBIQX,C215A,C215B                    

      ! jdemod - Most of LIM assumes that CIST is the elapsed time for the particle from initial injection with t0=0.0.
      ! However, updating the time dependence code allowing injection at different times breaks this assumption
      ! CIST is the elapsed time/particle while RTIME is the time since t=0.0 for the particle. 
      real, public :: cist, rtime

      INTEGER,public::         CIOPTA,CIOPTB,CIOPTC,CIOPTD,CIOPTE,CIOPTF                 
      INTEGER,public::         CIOPTG,CIOPTH,CIOPTI,CNEUTA,CNEUTB,CNEUTC                 
      INTEGER,public::         CMIZS,CIZ,CIAB,CNEUTD,CNEUTE,CNEUTF
      INTEGER,public::         CIOPTK,CIOPTL,CIOPTM,CIOPTN
      INTEGER,public::         CEXNEUT,NRFOPT,NVAOPT

  
  public :: allocate_mod_comtau, deallocate_mod_comtau


contains

  subroutine allocate_mod_comtau
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

  end subroutine allocate_mod_comtau


  subroutine deallocate_mod_comtau
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_comtau



end module mod_comtau
