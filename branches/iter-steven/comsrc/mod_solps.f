!     -*-Fortran-*-
c
c ======================================================================
c
      MODULE mod_solps_params
      IMPLICIT none
      PUBLIC

      INTEGER, PARAMETER :: 
     .  MAX_SOLPS_DATA = 20,
     .  SOLPS_DENSITY     = 1,
     .  SOLPS_VELOCITY    = 2,
     .  SOLPS_PRESSURE    = 3,
     .  SOLPS_TEMPERATURE = 4

c      CHARACTER*16, PARAMETER, PUBLIC ::
c     .  type_name(4) = ['density','velocity','pressure','temperaure']

      END MODULE mod_solps_params
c
c ======================================================================
c
      MODULE mod_solps
      IMPLICIT none
      PUBLIC

      TYPE :: type_solps_data
         REAL*4              :: version = 1.0
         CHARACTER*1024      :: fname          ! File name
         INTEGER             :: format         ! Format of the data file (1=AK as of 28.07.09)
         INTEGER             :: column         ! Column where the data is located (not including standard data columns to the left of the data area in the file)
         INTEGER             :: type           ! 1=density,2=velocity,3=pressure,4=temperature
         CHARACTER*64        :: tag            ! Descriptor for the particle species 
         INTEGER             :: z              ! Atomic number
         INTEGER             :: a              ! Atomic mass
         INTEGER             :: charge         ! Particle charge (-1=electrons,0=neutrals,+x=ions)
         REAL*4, ALLOCATABLE :: data(:)        ! Data
      ENDTYPE type_solps_data

      INTEGER :: nsolps_data
      TYPE(type_solps_data), PUBLIC, ALLOCATABLE :: solps_data(:)

      INTEGER :: solps_maxik,solps_maxir

      INTEGER, ALLOCATABLE :: map_divimp(:,:),map_osm(:)
      INTEGER, ALLOCATABLE :: solps_ik(:),solps_ir(:)
      REAL   , ALLOCATABLE :: solps_cen(:,:)

      END MODULE mod_solps
c
c ======================================================================
c

