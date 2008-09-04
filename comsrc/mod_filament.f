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
!       Timing:
        REAL*8    :: start_time                           ! Time when the filament was born
        REAL*8    :: t                                    ! Local time for filament (absolute time - start time)
!       Plasma parameters:
        REAL*8    :: ne(MAX_FIL_NCELL)                    ! Density scaling (m-3)
        REAL*8    :: vb(MAX_FIL_NCELL)                    ! Parallel flow (m s-1)
        REAL*8    :: te(MAX_FIL_NCELL)                    ! Te (eV)
        REAL*8    :: ti(MAX_FIL_NCELL)                    ! Ti (eV)
!       Cross section parameters:
        INTEGER*4 :: crs_opt                              ! Option (1=circle)
        REAL*8    :: crs_params(MAX_FIL_PARAMS)           ! 

!       Evolution along the field line:
        INTEGER*4 :: par_drift_opt                        ! Parallel momentum/flow option


!       Toroidal transport:


!       Radial transport:
        INTEGER*4 :: rad_opt                              ! Option
        REAL*8    :: rad_params(MAX_FIL_PARAMS)           ! Parameters for transport function
        INTEGER*4 :: rad_trigger_opt                      ! Option for under what conditions to initiate transport
        REAL*8    :: rad_trigger_params(MAX_FIL_PARAMS)   ! Parameters
        REAL*8    :: rad_delay                            ! Time delay before the transport is initiated
        REAL*8    :: rad_position                         ! Current position for filament "center of mass"
        INTEGER*4 :: rad_coordinate                       ! Coordinate used for tracking filament "center of mass"
        INTEGER*4 :: rad_status                           ! Register current transport activity (0-inactive, 1-active)
!       Geometry data:
        INTEGER*4 :: ncell                                ! Number of cross-sectional cells in the filament model  ! INTEGER*2's?
        INTEGER*4 :: nvcell                               ! Number of vertices per cell
        INTEGER*4 :: vcell(MAX_FIL_VCELL,MAX_FIL_NCELL)   ! Index in VTX array of of the vertices for a particular cell
        REAL*8    :: ccell(3,MAX_FIL_NCELL)               ! Cell centre... but should maybe be an INT pointer...
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
     .  fil_tor_params(MAX_FIL_PARAMS),  ! Parameters for toroidal distribution function 
     .  fil_par_params(MAX_FIL_PARAMS),  ! Parameters for parallel extent fuction
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
