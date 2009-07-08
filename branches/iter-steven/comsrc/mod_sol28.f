!     -*-Fortran-*-
!
! Prefix:
!
! P_
! M_
! TE_
! TI_
!
! PAR
! MOM
! ENE
! ENI
!
! Suffix:
!
! ION
! REC
! ANO
!
!
!
      MODULE mod_sol28_params
      IMPLICIT none
      PUBLIC

      INTEGER, PARAMETER ::   
     .  s28_MAXNKS   = 200, 
     .  s28_MAXNION  =  10, 
     .  s28_MAXPTS   =   4, 
     .  s28_MAXSIDE  =   6, 
     .  s28_MAXXPT   =   2,   
     .  s28_MAXVER   =   4, 
     .  S28_MAXNSRC  =  10,   
     .  ITY_FLUID    =   1,    
     .  ITY_KINETIC  =   2,   
     .  ST_STANDARD  =   1,   
     .  ST_NORMALIZE =   2,
     .  GRD_SOL      =   1,
     .  GRD_PFZ      =   2,
     .  GRD_CORE     =   3,
     .  GRD_TEST     =   4,
     .  GRD_BOUNDARY =  -1,
     .  EIR_MAXNPUFF =  10

      REAL, PUBLIC, PARAMETER ::   
     .  RHI   = 1.0E+37, 
     .  RLO   = 1.0E-37,
     .  V_PI  = 3.141593,
     .  EPS10 = 1.0E-10


      INTEGER, PUBLIC, PARAMETER ::  LO = 1, HI = 2, FULL = 3, TOTAL = 0

      INTEGER, PARAMETER ::   
     .  S28_MAXNNODE = 10

     


      END MODULE mod_sol28_params
c
c ======================================================================
c
      MODULE mod_sol28
      USE mod_sol28_params
      IMPLICIT none
      PRIVATE


      TYPE, PUBLIC :: type_options_osm
         REAL*4    :: version = 1.0

!...     I/O:
         INTEGER   :: log          ! Log file option
         INTEGER   :: logfp        ! File pointer for log file
         INTEGER   :: debug

         INTEGER   :: osm_load        ! Load status
         CHARACTER :: f_osm_load*512  ! Name of file to be loaded
         CHARACTER :: f_osm_dir*512   ! 

!...     Control:
         INTEGER   :: nflukin      ! Number of fluid-kinetic code iterations
         INTEGER   :: eirene       !   execute Eirene
         INTEGER   :: divimp       !   execute DIVIMP

         LOGICAL   :: pin_data     ! Flags whether or not PIN data is available
         INTEGER   :: cosm         ! OSM flux-tube loop internal iteration number for this run
         INTEGER   :: cflukin      ! OSM-kinetic code (Eirene,DIVIMP) iteration number


!...     Fluid solver:  
         INTEGER   :: bc(2)         ! Boundary conditions: 1=targets, 2=upstream  ... targets can be different?


         INTEGER   :: p_ion(2)         ! Ionisation 
         REAL      :: p_ion_frac(2)    ! Imposed ionisation bound relative to half-ring ion sink (fluxes + vol. rec.)
         INTEGER   :: p_rec(2)         ! Volume recombination 
         INTEGER   :: p_ano(2)         ! Anomalous

         INTEGER   :: m_mom(2)         ! ... needs new name...
         INTEGER   :: m_fit(2)         !
         INTEGER   :: m_ano(2)         !
         INTEGER   :: m_ano_dist(2)    ! Distribution along field line of anomalous momentum
         REAL      :: m_ano_exp(2)     ! Distribution exponent

         INTEGER   :: ti(2)            ! Ti calculation
         INTEGER   :: te_ano(2)        ! Distribution of anomalous (cross-field) electron power into the tube
         INTEGER   :: te_ano_psol(2)   ! Magnitude of  anomalous (cross-field) electron power flux
         INTEGER   :: te_fluxlimit(2)  ! Conduction factor
         REAL      :: te_kappa(2)      ! Conduction factor
         INTEGER   :: te_conv(2)       ! Electron energy convection parallel to the magnetic field
         INTEGER   :: te_ion(2)        ! Energy source from neutral code associated with ionisation (not strictly correct)
         INTEGER   :: te_rec(2)        ! Energy source from neutral code associated with volume recombination
         INTEGER   :: ti_ano(2)        ! Distribution of anomalous (cross-field) ion power into the tube
         INTEGER   :: ti_ano_psol(2)   ! Magnitude of anomalous (cross-field) ion power flux
         REAL      :: ti_kappa(2)      ! Conduction factor
         INTEGER   :: ti_conv(2)       ! Ion energy convection 
         REAL      :: ti_ratio(2)      ! Option for setting simple Te to Ti ratio
         INTEGER   :: ti_equil(2)      ! Te,Ti equilibration (fluid quantities only?)
         INTEGER   :: ti_ion(2)        ! Energy source from neutral code associated with ionisation (not strictly correct)
         INTEGER   :: ti_rec(2)        ! Energy source from neutral code associated with volume recombination

         INTEGER   :: super(2)         ! Allow for supersonic flows near targets
 
         REAL      :: s28mode          ! Node assignment routine
      ENDTYPE type_options_osm

      TYPE, PUBLIC :: type_options_eirene
         REAL*4    :: version = 1.0

!...     i/o:

!...     Options:    
         INTEGER   :: geom            ! geometry  2-triangles 3-tetrahedrons (toroidal)
         INTEGER   :: time            ! Execution time (CPU time in seconds)
         REAL      :: dtimv           ! ?
         INTEGER   :: opacity         ! Lyman alpha photon opacity option
         INTEGER   :: photons         ! Photon transport option
         INTEGER   :: trim            ! Fast ion surface collision database
         INTEGER   :: bgk             ! Neutral viscosity
         INTEGER   :: niter           ! Number of Eirene self-iterations
         INTEGER   :: data            ! Eirene input file  1=internal, 2=external 

!...     Particle sources:
         REAL      :: alloc           ! Flux / npts weighting (0.0 = npts only, 1.0 = flux only)
         REAL      :: puff_type(EIR_MAXNPUFF)
         INTEGER   :: puff_npts(EIR_MAXNPUFF)
         REAL      :: puff_flux(EIR_MAXNPUFF)
         REAL      :: puff_frac(EIR_MAXNPUFF)
         INTEGER   :: puff_species(EIR_MAXNPUFF)
         INTEGER   :: puff_index(EIR_MAXNPUFF)
         INTEGER   :: puff_energy(EIR_MAXNPUFF)
         INTEGER   :: puff_cosine_ind(EIR_MAXNPUFF)
         INTEGER   :: puff_cosine_max(EIR_MAXNPUFF)
         INTEGER   :: puff_pos(3,EIR_MAXNPUFF)
         INTEGER   :: puff_vec(3,EIR_MAXNPUFF)
         CHARACTER :: puff_note(EIR_MAXNPUFF)*1024


!...     Geometry / structure:
         INTEGER   :: ntorseg         ! Number of toroidal segments for regular toroidal descretization
         REAL      :: torfrac         ! Fraction of torus included in the simulation 
         INTEGER   :: ntri            ! Number of triangle mesh information lines
         REAL      :: tri             ! Triangle mesh specification
         INTEGER   :: mat1            ! Target material?
         INTEGER   :: mat2            ! Wall material?            
         REAL      :: ctargt          ! Target surface temperature
         REAL      :: cwallt          ! Wall surface temperature
      ENDTYPE type_options_eirene

!
!
!     Interpolation nodes:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_node
         REAL    :: type
         ! Geometry:
         INTEGER :: icell
         REAL    :: s
         ! Plasma:
         REAL    :: jsat(S28_MAXNION)
         REAL    :: ne
         REAL    :: v
         REAL    :: pe
         REAL    :: te
         REAL    :: te_exp
         REAL    :: ni(0:S28_MAXNION)
         REAL    :: pi(0:S28_MAXNION)
         REAL    :: ti(0:S28_MAXNION)
         REAL    :: ti_exp(0:S28_MAXNION)
         REAL    :: machno
         REAL    :: potential
         REAL    :: efield
         ! Interpolation parameters:
         INTEGER   :: tube_range(2)
         INTEGER   :: rad_mode
         REAL      :: rad_coord
         REAL      :: rad_exp
         REAL      :: rad_x
         REAL      :: rad_y
         INTEGER   :: par_mode
         REAL      :: par_exp
         INTEGER   :: par_set
         ! External data file:
         CHARACTER :: file_name*512
         REAL      :: file_shift
         REAL      :: file_scale_ne
         REAL      :: file_scale_M
         REAL      :: file_scale_pe
         REAL      :: file_scale_Te
         REAL      :: file_scale_Ti
         REAL      :: file_scale_V
         ! Radial fit parameters:
         REAL      :: fit_type
         REAL      :: fit_coord
         REAL      :: fit_psin(2)
         REAL      :: fit_shift
         REAL      :: fit_quantity
         REAL      :: fit_p(10)
      ENDTYPE type_node
!
!     Grid:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_grid
        REAL    :: version_grid
        REAL    :: version_tube
        REAL    :: version_cell
        REAL    :: version_pin
        REAL    :: version_photon
        REAL    :: version_kinetic
        REAL    :: version_drift
        REAL    :: version_fluid
        REAL    :: version_impurity
        INTEGER :: n                 ! = NTUBE
        INTEGER :: isep
        INTEGER :: ipfz
        INTEGER :: ikto
        INTEGER :: ikti
      ENDTYPE type_grid
!
!
!     Tubes:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_tube
        INTEGER :: index
        INTEGER :: region                ! core/SOL/PFZ/...
        INTEGER :: ir                    ! Radial index of tube ('ring number' in DIVIMP fluid grid)
        INTEGER :: it                    ! Toroidal index ('segment')
        INTEGER :: type                  ! 'idring' in DIVIMP
!...    Geometry data:
        INTEGER :: n                     ! Number of knots/cells in the flux-tube
        INTEGER :: ikti                       
        INTEGER :: ikto                  ! etc...              
        REAL    :: smax
        REAL    :: pmax
        REAL    :: rho
        REAL    :: psin
        REAL    :: metric(2)             ! Cross-field othogonality metric at ends of tube
!...    DIVIMP geometry data (temporary):
        REAL    :: bratio(2)             ! Field ratio (not totally correct)
        REAL    :: dds   (2)             ! Length of target segment (m)
        REAL    :: rp    (2)             ! Radial position (m) 
        REAL    :: costet(2)             ! Angle of field line intersection (degrees?)
!...    Species parameters:
        INTEGER :: nion                  ! Number of ion species
        REAL    :: zi(s28_MAXNION)       ! Atomic number
        REAL    :: ai(s28_MAXNION)       ! Atomic mass 
        INTEGER :: iontype(s28_MAXNION)  ! Ion data type: 1=fluid, 2=kinetic
        INTEGER :: ionmode(s28_MAXNION)  !          mode: 1=coupled, 2=free
!...    Target boundary conditions:
        REAL    :: jsat(2,s28_MAXNION)
        REAL    :: ne(2)       
        REAL    :: pe(2)       
        REAL    :: te(2)       
        REAL    :: ni(2,s28_MAXNION)       
        REAL    :: vi(2,s28_MAXNION)       
        REAL    :: machno(2)
        REAL    :: pi(2,s28_MAXNION)       
        REAL    :: ti(2,s28_MAXNION)       
!...    Cells:
        INTEGER :: cell_index(2)           ! Index of flux-tube in the CELL array
!...    Store solution parameters for upstream bc work:
        REAL*8 :: parano(3)                !    
        REAL*8 :: momano(3)                !    
        REAL*8 :: eneano(3)                ! 
        REAL*8 :: eniano(3,S28_MAXNION)    ! 
        REAL*8 :: parion(3,S28_MAXNION)    ! Ionisation
        REAL*8 :: parrec(3,S28_MAXNION)    ! Recombination
        REAL*8 :: parvol(3,S28_MAXNION)    ! Recombination
        REAL*8 :: momvol(3,S28_MAXNION)    ! Momentum - neutrals
        REAL*8 :: enerec(3)                ! Electron energy source - neutral volume recombination
        REAL*8 :: eneion(3)                ! Electron energy source - neutral ionisation
        REAL*8 :: enevol(3)                ! Electron energy source - total volume source
        REAL*8 :: enirec(3,S28_MAXNION)    ! Ion - neutral volume recombination
        REAL*8 :: eniion(3,S28_MAXNION)    ! Ion - neutral ionisation
        REAL*8 :: enivol(3,S28_MAXNION)    ! Ion - total
        REAL*8 :: te_kappa(2)
        REAL*8 :: ti_kappa(2,S28_MAXNION)
      ENDTYPE type_tube
!
!     Cells:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_cell
        INTEGER :: ik                               ! Index of cell in standard 2D grid
!...    Geometry data:
        REAL    :: cencar(3)                        ! Cell caresian center in machine coordinates
        REAL    :: centor(3)                        ! Cell toroidal center in 'toroidal' coordinates
        REAL    :: vol                              ! Cell volume 
        INTEGER :: nside                            ! Number of sides for each cell: nside=6 usually
!         ...need some reference to location relative to x-point...
        REAL    :: s                                ! Distance of cell center along the magnetic field line (m), s=0 at inner target (LO index target)
        REAL    :: p                                ! Poloidal distance of cell center from inner target (LO index target)
        REAL    :: sbnd(2)                          ! 's' at ends of each cell
        REAL    :: pbnd(2)                          ! 'p' at ends of each cell
        REAL    :: ds
        REAL    :: dp
        REAL    :: metric                           ! Cross-field othogonality metric
        REAL    :: cfdist(s28_MAXSIDE)              ! Cross-field distance from cell center to side in machine coordinates
!...    Cell side data:   (FOR DELETION!  DO NOT USE!)
        INTEGER :: sidetype   (s28_MAXSIDE)         ! Side geometry category
        INTEGER :: orientation(s28_MAXSIDE)         ! Orientation (?)
        INTEGER :: iside      (s28_MAXSIDE)         ! Side index in surface array
        INTEGER :: sidemap    (s28_MAXSIDE,3)       ! Map to 'best' neighbour (tube index, cell, side)
        REAL    :: sidearea   (s28_MAXSIDE)         ! Geometric area of side
        REAL    :: aproj      (s28_MAXSIDE)         ! Projection of area relative to the magnetic field
!...    Old plasma solution: (?)
      ENDTYPE type_cell

!...  Neutral particle data (including test ions):
!     --------------------------------------------------
      TYPE, PUBLIC :: type_neutral
        REAL :: n_atm                           ! Atom density (m-3)
        REAL :: n_mol                           ! Molecule density (m-3)
        REAL :: ion                             ! ionisation source (m-3 s-1)
        REAL :: rec                             ! Volume recombination sink (m-3 s-1)
        REAL :: mom                             ! Parallel momentum source (...)
        REAL :: qe                              ! Electron energy source (...)
        REAL :: qi                              ! Ion energy source (...)
        REAL :: dalpha(10)                      ! Dalpha emission, 1 = total, 2-6 = components (photons m-3 s-1)
        REAL :: dgamma(10)                      ! Dgamma emission, 1 = total, 2-5 = components (photons m-3 s-1)
      ENDTYPE type_neutral

!...  Kinetic photon data:
!     --------------------------------------------------
      TYPE, PUBLIC :: type_photon
        REAL :: n
      ENDTYPE type_photon

!...  Drift terms:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_drift
        REAL :: radExB
        REAL :: diaExB
      ENDTYPE type_drift

!...  Vector fields (magnetic and electric):
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_field
!...    Magneitic field data at cell center:
        REAL    :: bratio                           ! Ratio of poloidal to total magnetic fields
        REAL    :: bpot                             ! Magnetic potential
        REAL    :: b                                ! Total strength of the magnetic field (Tesla)
        REAL    :: br                               ! Radial component (Tesla)
        REAL    :: bphi                             ! Toroidal component (Tesla)
        REAL    :: bz                               ! Vertical component (Tesla)
!...    Electric field data at cell center:
        REAL    :: epot                             ! Electric potential
        REAL    :: efield                           ! Electric field 
        REAL    :: envec(3)                         ! Caresian basis vector 
        REAL    :: envbf(3)                         ! Basis vector relative to the magnetic field
      ENDTYPE type_field

!...  Impurities:
!     --------------------------------------------------
      TYPE, PUBLIC :: type_impurity
        REAL :: sdlims
      ENDTYPE type_impurity

!...  External fluid code results (most likely imported from 2D fluid code):
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_fluid
!...    Fluid quantities at cell center:
        REAL    :: te                   ! Electron Maxwellian temperature (eV)
        REAL    :: ne                   ! Electron density (m-3)
        REAL    :: ni                   ! Plasma ion density (m-3)
        REAL    :: vi                   ! Plasma ion velocity parallel to the magnetic field
        REAL    :: ti                   ! Plasma ion temperature (eV)
!...    Interpolated fluid quantities at surfaces of 
!       the fluid grid cell (at the surface center, 
!       approximately):
        REAL    :: vesurf(s28_MAXSIDE)
        REAL    :: tisurf(s28_MAXSIDE)
        REAL    :: nisurf(s28_MAXSIDE)  
        REAL    :: visurf(s28_MAXSIDE)  ! Parallel or perpendicular, depending on cell surface
!...    Gradients at surfaces (at the surface 
!       center, approximately):
        REAL    :: tegrad(s28_MAXSIDE)
        REAL    :: negrad(s28_MAXSIDE) 
        REAL    :: vegrad(s28_MAXSIDE)  ! Parallel velocity gradient
        REAL    :: tigrad(s28_MAXSIDE)
        REAL    :: nigrad(s28_MAXSIDE)
        REAL    :: vigrad(s28_MAXSIDE)  ! Parallel velocity gradient
!...    Diffusive cross-field transport 
!       coefficients: 
        REAL    :: dperp                ! Particle
        REAL    :: mperp                ! ??? Momentum 
        REAL    :: chiperp              ! Energy
!...    Fluxes through surfaces:
        REAL    :: parflx(s28_MAXSIDE)  ! Particle 
        REAL    :: momflx(s28_MAXSIDE)  ! Momentum
        REAL    :: eneflx(s28_MAXSIDE)  ! Energy
        REAL    :: eniflx(s28_MAXSIDE)
!...    Volume sources:  
        REAL    :: parion               ! Ionisation
        REAL    :: parrec               ! Recombination
        REAL    :: parvol               ! Net volume particle source
        REAL    :: momvol               ! Momentum
        REAL    :: enerec               ! Energy - electron volume recombination
        REAL    :: eneion               ! Energy - electron ionisation
        REAL    :: enevol               ! Energy - electron total
        REAL    :: enirec               ! Energy - ion volume recombination
        REAL    :: eniion               ! Energy - ion ionisation
        REAL    :: enivol               !        - ion total
!...    Anomalous sources:      
        REAL    :: parano    
        REAL    :: momano  
        REAL    :: eneano  
        REAL    :: eniano  
!...    User sources:
        REAL    :: parusr               ! Particle
        REAL    :: momusr               ! Momentum
        REAL    :: eneusr               ! Energy
        REAL    :: eniusr
!...    Total sources:
        REAL    :: parsrc               ! Particle
        REAL    :: momsrc               ! Momentum
        REAL    :: enesrc               ! Energy
        REAL    :: enisrc
      ENDTYPE type_fluid

!...  Kinetic ion code results (Eirene/DIVIMP):
!     --------------------------------------------------
      TYPE, PUBLIC :: type_kinetic
        REAL :: te
      ENDTYPE type_kinetic

!...  
      REAL, PUBLIC, ALLOCATABLE, SAVE :: tmpflx(:,:)

      END MODULE mod_sol28
!
!
!
!
!
! ----------------------------------------------------------------------
!
      MODULE mod_sol28_global
      USE mod_sol28      
      USE mod_geometry
      IMPLICIT none

      PUBLIC

      INTEGER, PARAMETER :: IND_IK      = 1,  ! Object index parameters
     .                      IND_IR      = 2,  
     .                      IND_IS      = 3,  
     .                      IND_ZONE    = 4,  
     .                      IND_FLUID   = 5,  
     .                      IND_KINETIC = 6,  
     .                      IND_NEUTRAL = 7,  
     .                      IND_FIELD   = 8,  ! Vaccum zone outside standard grid, from external call to TRIANGLE
     .                      IND_CELL    = 9

      INTEGER, PARAMETER :: IND_STDGRD  = 1,  ! Magnetic fluid grid side index, i.e. 12, 23, 34, 41
     .                      IND_TARGET  = 2,  ! Target (block 7 stratum in Eirene input file)
     .                      IND_SURFACE = 3   ! Surface (block 2A non-default surface in Eirene input file)


!...
      TYPE(type_options_osm   ) :: opt
      TYPE(type_options_eirene) :: opt_eir
      INTEGER log, logfp

      INTEGER osmnnode    
      TYPE(type_node) :: osmnode(100)
!...  
      INTEGER, PARAMETER :: MAXNTUBE = 100
      TYPE(type_grid), SAVE :: grid
      INTEGER, SAVE :: ntube
      TYPE(type_tube), ALLOCATABLE, SAVE :: tube(:) 

      INTEGER, PARAMETER :: MAXNCELL = MAXNTUBE * 100
      INTEGER, SAVE :: nion,ncell,nfield,npin,nphoton,nfluid,nkinetic,
     .                 nimpurity,ndrift
      TYPE(type_cell         ), ALLOCATABLE, SAVE :: cell    (:) 
      TYPE(type_field   ), ALLOCATABLE, SAVE :: field   (:) 
      TYPE(type_neutral ), ALLOCATABLE, SAVE :: pin     (:,:) 
      TYPE(type_photon  ), ALLOCATABLE, SAVE :: photon  (:,:) 
      TYPE(type_kinetic ), ALLOCATABLE, SAVE :: kinetic (:,:) 
      TYPE(type_fluid   ), ALLOCATABLE, SAVE :: fluid   (:,:)
      TYPE(type_impurity), ALLOCATABLE, SAVE :: impurity(:,:) 
      TYPE(type_drift   ), ALLOCATABLE, SAVE :: drift   (:,:) 

!...  Reference plasma:
      INTEGER, SAVE :: ref_ntube,ref_nion,ref_nfluid
      TYPE(type_tube), ALLOCATABLE :: ref_tube(:)
      TYPE(type_fluid), ALLOCATABLE, SAVE :: ref_fluid(:,:) 

      END MODULE mod_sol28_global
!
!
!
!
!
! ----------------------------------------------------------------------
!
      MODULE mod_sol28_solver
      USE mod_sol28_params
      USE mod_sol28      
      IMPLICIT none
      PUBLIC


      LOGICAL, SAVE :: output


      REAL*8, PARAMETER :: ECH = 1.602D-19, AMU = 1.67D-27


      TYPE(type_options_osm), SAVE :: opt
      INTEGER :: log,logfp

!...  Local solver options:
      INTEGER, SAVE ::    
c     .  opt_bc(2),             ! Method of assigning boundary conditions (1=target, 2=upstream)
c     .  opt_p_ion(2),             
     .  opt_p_ion_scale(2)
c     .  opt_p_rec(2),             
c     .  opt_p_ano(2),              
c     .  opt_p_radExB(2),   
c     .  opt_m_mom(2),             
c     .  opt_m_fit(2),
c     .  opt_m_ano(2),   
c     .  opt_m_ano_dist(2),   
c     .  opt_te_ano(2),
c     .  opt_te_fluxlimit(2),
c     .  opt_te_conv(2),
c     .  opt_ti(2),
c     .  opt_ti_ano(2),
c     .  opt_ti_conv(2),
c     .  opt_ti_equil(2),
!...   
c     .  opt_super(2)                  ! Allow sonic transition      

c      REAL*8
c     .  opt_p_ion_frac(2),
c     .  opt_m_ano_exp(2),
c     .  opt_ti_ratio(2),
c     .  opt_te_kappa(2),
c     .  opt_ti_kappa(2)


!...  Process control:
      LOGICAL, SAVE ::   
     .  cnt_options,   
     .  cnt_target,
     .  cnt_boundary_conditions,
     .  cnt_super(2),   
     .  cnt_prescription,   
     .  cnt_particles,   
     .  cnt_momentum,
     .  cnt_energy,
     .  cnt_integrate


!...  Solution status:
      INTEGER, SAVE ::   
     .  anl_ic_super(2)

      LOGICAL, SAVE ::   
     .  anl_detached(2),   
     .  anl_imaginary(2)

!...  Manipulation:

!...
      INTEGER, SAVE :: nnode,mnode
      TYPE(type_node), SAVE :: node(S28_MAXNNODE)

      TYPE(type_tube), SAVE :: tube
      TYPE(type_cell        ), SAVE :: cell     (S28_MAXNKS)
      TYPE(type_neutral), SAVE :: pin      (S28_MAXNKS,S28_MAXNION) 
      TYPE(type_fluid  ), SAVE :: fluid    (S28_MAXNKS,S28_MAXNION) 

      INTEGER ref_nion
      TYPE(type_tube) :: ref_tube
      TYPE(type_fluid), SAVE :: ref_fluid(S28_MAXNKS,S28_MAXNION) 

!...  
      INTEGER, SAVE ::   
     .  icmid,                  
     .  icmax,                  
     .  icbnd1(3),              
     .  icbnd2(3),              
     .  ictarg(2),              ! Target cell indeces
     .  intarg(2),              ! Target interpolation node indeces
     .  nion,                   
     .  iontype(S28_MAXNION),   
     .  nsrc(3),                
     .  srctype(S28_MAXNSRC)

      REAL*8, SAVE ::   
     .  tsign(2),
     .  gamma(2)

!...
      REAL*8, SAVE ::   
     .  zi(0:S28_MAXNION),     ! Ion nuclear charge
     .  ai(0:S28_MAXNION),     ! Ion nuclear mass number
     .  mi(0:S28_MAXNION),     ! Ion nuclear mass in kg
     .  ci(0:S28_MAXNION)      ! Ion electron charge state

!...  Solver numerical variables:
      REAL*8, SAVE ::   
!...    Geometry:
     .  smax              

      REAL*8, TARGET, SAVE ::   
!...    Geometry:
     .  sfor(0:S28_MAXNKS+1),   
     .  sbak(0:S28_MAXNKS+1)             


      REAL*8, SAVE ::   
!...    Geometry:
     .  sbnd  (2,S28_MAXNKS),
     .  sdelta(1:S28_MAXNKS),
     .  area(0:S28_MAXNKS+1),   
     .  vol (  S28_MAXNKS)  ,   
!...    Plasma fluid quantities:
     .  isat(0:S28_MAXNKS+1,0:S28_MAXNION),  
     .  ne  (0:S28_MAXNKS+1),                 
     .  ni  (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  vi  (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  pe  (0:S28_MAXNKS+1),                 
     .  pi  (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  te  (0:S28_MAXNKS+1),                 
     .  ti  (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  qe  (0:S28_MAXNKS+1),                 
     .  qi  (0:S28_MAXNKS+1,0:S28_MAXNION),
     .  qcond(0:S28_MAXNKS+1),
     .  qconv(0:S28_MAXNKS+1),
!...
     .  machno(0:S28_MAXNKS+1,0:S28_MAXNION),  
!...    ...:
     .  par(0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  mom(0:S28_MAXNKS+1,0:S28_MAXNION),   
!...    Anomalous sources:  
     .  parano(S28_MAXNKS,0:S28_MAXNION),   
     .  momano(S28_MAXNKS,0:S28_MAXNION),   
     .  eneano(S28_MAXNKS)              ,
     .  eniano(S28_MAXNKS,0:S28_MAXNION),   
!...    User sources:  
     .  parusr(S28_MAXNKS,0:S28_MAXNION),   
     .  momusr(S28_MAXNKS,0:S28_MAXNION),   
     .  eneusr(S28_MAXNKS)              ,
     .  eniusr(S28_MAXNKS,0:S28_MAXNION),   
!...    Cross-field sources:
     .  parflx(S28_MAXNKS,0:S28_MAXNION),   
     .  momflx(S28_MAXNKS,0:S28_MAXNION),   
     .  engflx(S28_MAXNKS,0:S28_MAXNION),   
!...    Volume sources:
     .  parion(S28_MAXNKS,0:S28_MAXNION),   
     .  parrec(S28_MAXNKS,0:S28_MAXNION),   
     .  momvol(S28_MAXNKS,0:S28_MAXNION),   
     .  enerec(S28_MAXNKS,0:S28_MAXNION),
     .  eneion(S28_MAXNKS,0:S28_MAXNION),
     .  enevol(S28_MAXNKS,0:S28_MAXNION),
     .  enirec(S28_MAXNKS,0:S28_MAXNION),   
     .  eniion(S28_MAXNKS,0:S28_MAXNION),   
     .  enivol(S28_MAXNKS,0:S28_MAXNION),   
!...    Net sources:
     .  parsrc(0:S28_MAXNKS+1,S28_MAXNION),   
     .  momsrc(0:S28_MAXNKS+1,S28_MAXNION),   
     .  enesrc(0:S28_MAXNKS+1,0:1),   
     .  enisrc(0:S28_MAXNKS+1,S28_MAXNION),   
!...    Net source integrals along the flux tube:
     .  parint(0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  momint(0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  eneint(0:S28_MAXNKS+1,0:1),
     .  eniint(0:S28_MAXNKS+1,0:S28_MAXNION)


      END MODULE mod_sol28_solver
!
! ----------------------------------------------------------------------
!
      MODULE mod_osm_input
      IMPLICIT none
      PUBLIC

      CHARACTER*256 :: eircpuff(100)      

      END MODULE mod_osm_input
