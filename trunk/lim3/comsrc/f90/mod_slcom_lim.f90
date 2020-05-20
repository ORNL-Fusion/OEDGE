module mod_slcom

  use mod_params

!      COMMON /SLCOM/ 
!c
!     + SLOPT,
!     + optdp,
!     + N2OPT,NEOPT,NIEOPT,
!     + N2FRAC,NENG,NIENG,
!     + N2RATE,NNRATE,NIRATE,N2CPCHS,NNCPCHS,
!
!     + n2_n,n_n,ni_n,
!     + n2_lim,n_lim,      
!     + n2_x,n2_y,n2_p,n_x,n_y,n_p,n_x_max,n_y_max,n_p_max,
!     + n2_e,n_e,ni_e,n2_v,n_v,
!     + n2_t,n_t,
!     + n_beta,n_psi,
!
!     + TITL2,LAMBDA,INJBINT,TAG2,MARK,
!     + TLOSS,WLOSS,TGLOSS,LLOSS,ALOSS,IZLOSS,TSLOSS,TOTVOL,TOTDEN,
!     + BSBIN,YSBIN,IZBIN,NBIN,NSBIN,MIZS,MIMPS,MATIZ,
!     + SVHINS,SEYINS            
!c 
!c ----------------------------------------------------------------------
!c
!c Declarations:
!c
!c Output file unit number - set to general debug unit for now
!c
!      integer slout
!      parameter (slout=dbgunit) 
!
!c
!c DIVIMP ion profile option:
!c
!      INTEGER optdp
!
!      INTEGER SLOPT,N2OPT,NEOPT,NIEOPT
!      REAL    N2FRAC,NENG,NIENG
!      REAL    N2RATE(MAXNXS),NNRATE(MAXNXS),NIRATE(MAXNXS),
!     +        N2CPCHS(MAXNXS),NNCPCHS(MAXNXS)
!c
!c Statistics for N2 break-up:
!c
!      REAL    n2_n,n_n,ni_n
!      REAL    n2_lim,n_lim      
!      REAL    n2_x,n2_y,n2_p,n_x,n_y,n_p,n_x_max,n_y_max,n_p_max
!      REAL    n2_e,n_e,ni_e,n2_v,n_v
!      REAL    n2_t,n_t
!      REAL    n_beta,n_psi
!c
!c LIM3 variables (for the most part):
!c
!      CHARACTER  TITL2*80
!      REAL       LAMBDA
!      INTEGER    INJBINT(-MAXY3D-1:MAXY3D+1)
!c      INTEGER    INJBINP(-MAXNPS:MANNPS)
!      REAL       TAG2(MAXIMP)
!      REAL       MARK
!c
!c Statisitcs:
!c
!      REAL       TLOSS(MAXIMP),WLOSS,TGLOSS,LLOSS,ALOSS,IZLOSS,TSLOSS
!      REAL       TOTVOL,TOTDEN(MAXIZS)
!c
!c Variables for storing the DIVIMP source distribution:
!c
!      REAL       BSBIN(10),YSBIN(10),IZBIN(10)
!      INTEGER    NBIN, NSBIN(10,10),MIZS,MIMPS,MATIZ
!c
!c Moved from LIM3:
!c
!      REAL      SVHINS(-MAXQXS:MAXQXS),SEYINS(-MAXQXS:MAXQXS,MAXIZS)            


  implicit none
  private

! Declarations:
!
! Output file unit number - set to general debug unit for now
! This is set by a call to the initialization routine at the beginning of runlm3
! Set a default value of 6 which is the baseline output unit number in LIM/DIVIMP  
      integer,public:: slout = 6
      !parameter (slout=dbgunit) 

!
! DIVIMP ion profile option:
!
      INTEGER,public:: optdp

      INTEGER,public:: SLOPT,N2OPT,NEOPT,NIEOPT
      REAL,public::    N2FRAC,NENG,NIENG
      REAL,public::    N2RATE(MAXNXS),NNRATE(MAXNXS),NIRATE(MAXNXS),N2CPCHS(MAXNXS),NNCPCHS(MAXNXS)
!
! Statistics for N2 break-up:
!
      REAL,public::    n2_n,n_n,ni_n
      REAL,public::    n2_lim,n_lim      
      REAL,public::    n2_x,n2_y,n2_p,n_x,n_y,n_p,n_x_max,n_y_max,n_p_max
      REAL,public::    n2_e,n_e,ni_e,n2_v,n_v
      REAL,public::    n2_t,n_t
      REAL,public::    n_beta,n_psi
!
! LIM3 variables (for the most part):
!
      CHARACTER,public::  TITL2*80
      REAL,public::       LAMBDA
      INTEGER,public::    INJBINT(-MAXY3D-1:MAXY3D+1)
!      INTEGER,public::    INJBINP(-MAXNPS:MANNPS)
      REAL,public::       TAG2(MAXIMP)
      REAL,public::       MARK
!
! Statisitcs:
!
      REAL,public::       TLOSS(MAXIMP),WLOSS,TGLOSS,LLOSS,ALOSS,IZLOSS,TSLOSS
      REAL,public::       TOTVOL,TOTDEN(MAXIZS)
!
! Variables for storing the DIVIMP source distribution:
!
!      REAL,public::       BSBIN(10),YSBIN(10),IZBIN(10)
!
      REAL,public::       IZBIN(10)


      real,public::  YSBIN(10) = (/ -1.0E+06, -3.2317, -2.7438, -1.9309, &
                    -1.1181, -0.3259,  0.5338,  1.7560,&
                     3.7696, 1.0E+06 /)
      real,public::  BSBIN(10) = (/  -3.1231, -3.0029, -2.5149, -1.7021,&
                     -0.8893, -0.0971,  0.7625,  1.9847,&
                      3.9985, 1.0E+06 /)

      integer,public::  NBIN = 10



      INTEGER,public::    NSBIN(10,10),MIZS,MIMPS,MATIZ
!
! Moved from LIM3:
!
      REAL,public::      SVHINS(-MAXQXS:MAXQXS),SEYINS(-MAXQXS:MAXQXS,MAXIZS)            
 
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

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_slcom


  subroutine deallocate_mod_slcom
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_slcom



end module mod_slcom
