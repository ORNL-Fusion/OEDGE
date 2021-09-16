!     -*-Fortran-*-
!
!
      MODULE mod_filament_params
      IMPLICIT none
      PUBLIC

      INTEGER, PARAMETER ::   
     .  MAX_FIL_NVTX   = 1200,
     .  MAX_FIL_NCELL  = 300,   ! *** HACK ** This is to allow the 1:1 correspondence between vertices and cells, for now...
     .  MAX_FIL_VCELL  = 4,
     .  MAX_FIL_PARAMS = 20,
     .  MAX_IR         = 20

      REAL, PARAMETER :: 
     .  PI  = 3.141593

      REAL*8, PARAMETER :: 
     .  DPI = 3.14159265358979D0 

      END MODULE mod_filament_params
c
c ======================================================================
c
      MODULE mod_filament
      USE mod_filament_params
      IMPLICIT none
      PRIVATE

      TYPE, PUBLIC :: type_filament
        REAL*4    :: version = 1.0
        INTEGER*4 :: status                               ! Is the filament active?
!       Timing:
        REAL*8    :: t_start                              ! Time when the filament was born (s)
        REAL*8    :: t_end                                ! Time it's finished (s)
        REAL*8    :: t_duration                           ! How long the filament lasts (s)
        REAL*8    :: t_last                               ! Global time stamp when filament was last updated (s)
!       Plasma parameters:
!        REAL*4    :: ne(MAX_FIL_NCELL)                    ! Density scaling (m-3)
!        REAL*4    :: vb(MAX_FIL_NCELL)                    ! Parallel flow (m s-1)
!        REAL*4    :: te(MAX_FIL_NCELL)                    ! Te (eV)
!        REAL*4    :: ti(MAX_FIL_NCELL)                    ! Ti (eV)
!        REAL*4    :: ne0                                  ! Principle density parameter (m-3)
!        REAL*4    :: vb0                                  !     "     parallel flow (m s-1)
!        REAL*4    :: te0                                  !     "     Te (eV)
!        REAL*4    :: ti0                                  !     "     Ti (eV)

!       Plasma density:
        INTEGER*4 :: ne_opt                              ! Option (1=...)
        REAL      :: ne_param(MAX_FIL_PARAMS)            ! 
!       Plasma temperature:
        INTEGER*4 :: te_opt                              ! Option (1=...)
        REAL      :: te_param(MAX_FIL_PARAMS)            ! 

!       Cross section parameters:
        INTEGER*4 :: crs_opt                              ! Option (1=circle)
        REAL*8    :: crs_param(MAX_FIL_PARAMS)           ! 
        INTEGER*4 :: crs_plasma                           ! Cross-field plasma model option
!       Parallel evolution along the field line:
        INTEGER*4 :: par_plasma                           ! Parallel plasma model option
        INTEGER*4 :: par_length                           ! Extent of the filament at creation (t = 0) (1-"full", 2-outer)
        INTEGER*4 :: par_growth                           ! Rate of expansion along the field line
        REAL*8    :: par_growth_param(MAX_FIL_PARAMS)    ! ...parameters
        REAL*4    :: par_growth_delay                     ! Time delay to impose after the trigger is registered (s)
!       Toroidal transport:
        INTEGER*4 :: tor_opt                              ! Option
        REAL*8    :: tor_param(MAX_FIL_PARAMS)           ! Parameters for transport function
        REAL*8    :: tor_delay                            ! Time delay before the transport is initiated
        REAL*8    :: tor_position                         ! Current position for filament "center of mass"
        INTEGER*4 :: tor_coordinate                       ! Coordinate used for tracking filament "center of mass"
        INTEGER*4 :: tor_status                           ! Register current transport activity (0-inactive, 1-active)
!       Radial transport:
        INTEGER*4 :: rad_opt                              ! Option
        REAL*8    :: rad_param(MAX_FIL_PARAMS)           ! Parameters for transport function
        REAL*8    :: rad_delay                            ! Time delay before the transport is initiated
        REAL*8    :: rad_position                         ! Current position for filament "center of mass"
        INTEGER*4 :: rad_coordinate                       ! Coordinate used for tracking filament "center of mass"
        INTEGER*4 :: rad_status                           ! Register current transport activity (0-inactive, 1-active)
!       Geometry data:
        INTEGER*4 :: ncell                                ! Number of cross-sectional cells in the filament model  ! INTEGER*2's?
        INTEGER*4 :: nvcell                               ! Number of vertices per cell
        INTEGER*4 :: vcell(MAX_FIL_VCELL,MAX_FIL_NCELL)   ! Index in VTX array of of the vertices for a particular cell
        REAL*8    :: lcell(2,MAX_FIL_NCELL)               ! Length of the cell/tube along the magnetic field
        REAL*8    :: lcell0(2,MAX_FIL_NCELL)              ! Initial length at t = 0
        INTEGER*4 :: nvtx                                 ! Total number of vertices used to represent the cross-section
        REAL*8    :: vtx(3,MAX_FIL_NVTX)                  ! Vertex data
        INTEGER   :: ir_space(0:MAX_IR,MAX_FIL_NCELL)     ! Range of flux rings to search when resolving tetrahedral mesh

      ENDTYPE type_filament
!
!     ------------------------------------------------------------------
      INTEGER, PUBLIC :: 
     .  fil_tor_opt,                     ! Option for toroidal distribution of filaments
     .  fil_tor_n,                       ! Quasi-toroidal mode number
     .  fil_par_opt,                     ! Parallel extent option
     .  fil_time_opt,                    ! Option for keeping track of when filaments are launched 
     .  fil_trace_opt                    ! Option for field line tracing   
      REAL*4, PUBLIC :: 
     .  fil_tor_param(MAX_FIL_PARAMS),  ! Parameters for toroidal distribution function 
     .  fil_par_param(MAX_FIL_PARAMS),  ! Parameters for parallel extent fuction
     .  fil_time_start,                  ! Start time of filament clock for the first iteration of the simulation
     .  fil_time_delay,                  ! Delay before launching filaments
     .  fil_time_incriment,              ! Time increase between interations
     .  fil_time_duration,               ! How long to keep launching filaments
     .  fil_rad_extent,                  ! How big is the filament in the radial direction
     .  fil_tor_extent,                  ! How big in the toroidal direction
     .  fil_par_extent                   ! How big in the parallel direction


      INTEGER, PUBLIC :: nfilament
      TYPE(type_filament), PUBLIC, ALLOCATABLE :: filament(:)
!
!     ------------------------------------------------------------------
      CONTAINS
!
! ======================================================================
! ======================================================================
!
! subroutine: Dummy
!
      SUBROUTINE Dummy
      IMPLICIT none
      RETURN
99    STOP
      END SUBROUTINE Dummy

! ======================================================================
! ======================================================================
      END MODULE MOD_FILAMENT
