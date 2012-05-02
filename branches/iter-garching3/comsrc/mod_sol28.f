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
     .  s28_MAXNKS   = 500, 
     .  s28_MAXNION  =  10, 
     .  s28_MAXPTS   =   4, 
     .  s28_MAXSIDE  =   6, 
     .  s28_MAXXPT   =   2,   
     .  s28_MAXVER   =   4, 
     .  S28_MAXNSRC  =  10,   
     .  S28_MAXNTDEL =  100,   
     .  ITY_FLUID    =   1,    
     .  ITY_KINETIC  =   2,   
     .  ST_STANDARD  =   1,   
     .  ST_NORMALIZE =   2,
     .  GRD_SOL      =   1,  ! Should be consistent with the SOL1 definition in SLCOM in DIVIMP
     .  GRD_PFZ      =   2,  !  same
     .  GRD_CORE     =   3,  !  same
     .  GRD_TEST     =   4,
     .  GRD_BOUNDARY =  -1,
     .  EIR_MAXNSTRATA = 100,
     .  EIR_MAXNVOID   = 100,
     .  EIR_MAXNADD    = 500,
     .  EIR_MAXNSPECTRA= 100,
     .  EIR_MAXNSUR    = 100,
     .  EIR_MAXNTET    = 100,
     .  LSND = 1, USND = 2, UDND = 3, LDND = 4, CDND = 5, LINEAR = 6


      REAL, PUBLIC, PARAMETER ::   
     .  RHI   = 1.0E+37, 
     .  RLO   = 1.0E-37,
     .  V_PI  = 3.141593,
     .  EPS10 = 1.0E-10

      REAL*8, PUBLIC, PARAMETER ::   
     .  DPS10 = 1.0D-15


      INTEGER, PUBLIC, PARAMETER ::  LO = 1, HI = 2, FULL = 3, TOTAL = 0

      INTEGER, PARAMETER ::   
     .  S28_MAXNNODE = 10

     


      END MODULE mod_sol28_params

!
! ======================================================================
!
      MODULE mod_sol28_io
      IMPLICIT none
      PUBLIC

      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2, ALL_LINES = 3

      END MODULE mod_sol28_io
c
c ======================================================================
c
      MODULE mod_sol28
      USE mod_sol28_params
      IMPLICIT none
      PRIVATE

      LOGICAL, PUBLIC :: cell_modified, tube_modified  ! *** make PRIVATE when code update to use Alloc_cell, AddCell, etc and they are resident in this module

!     OSM options:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_options_osm
         REAL*4    :: version = 1.0

!...     Applicability:
         INTEGER   :: iteration(2)
         CHARACTER :: tube*128
!         INTEGER   :: tube(2)

!...     I/O:
         INTEGER   :: log          ! Log file option
         INTEGER   :: logfp        ! File pointer for log file
         INTEGER   :: debug

         INTEGER   :: osm_load         ! Load status
         CHARACTER :: f_osm_dir*512    ! 
         CHARACTER :: f_osm_load*512   ! Name of file to be loaded
!...     Grid:
         INTEGER   :: f_grid_format            ! Format of equilibrium grid to be loaded
         INTEGER   :: f_grid_load_method       ! Geometry load scheme 1-DIVIMP files, 2-OSM geometry setup
         CHARACTER :: f_grid_file*512          ! Name of equilibrium grid to be loaded
         INTEGER   :: f_grid_strip             ! Remove boundary cells (1=first and last cells, first and last rings)
         INTEGER   :: grd_ntdel               ! Number of tubes to delete after loading the grid
         INTEGER   :: grd_tdel(S28_MAXNTDEL)  ! List of tubes to delete
!...     Flow control:
         INTEGER   :: nflukin      ! Number of fluid-kinetic code iterations
         INTEGER   :: eirene       !   execute Eirene
         INTEGER   :: divimp       !   execute DIVIMP

         LOGICAL   :: pin_data     ! Flags whether or not PIN data is available
         INTEGER   :: cosm         ! OSM flux-tube loop internal iteration number for this run
         INTEGER   :: cflukin      ! OSM-kinetic code (Eirene,DIVIMP) iteration number
!...     Fluid solver:  
         INTEGER   :: sol_n
         INTEGER   :: sol_tube  (2,100)
         INTEGER   :: sol_option(  100)
         INTEGER   :: bc(2)         ! Boundary conditions: 1=targets, 2=upstream  ... targets can be different?
         INTEGER   :: p_ion(2)         ! Ionisation 
         REAL      :: p_ion_exp(2)     ! Exponent for exponential decay of the ionisation source for P_ION = 3
         REAL      :: p_ion_frac(2)    ! Imposed ionisation bound relative to half-ring ion sink (fluxes + vol. rec.)
         INTEGER   :: p_rec(2)         ! Volume recombination 
         INTEGER   :: p_ano(2)         ! Anomalous
         INTEGER   :: p_ano_dist(2)    ! Distribution along field line of anomalous particle flux
         REAL      :: p_ano_exp(2)     ! Distribution exponent
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

         INTEGER   :: radvel           ! Option for specifying the radial velocity with MODE=7 in S74 (far-SOL)
         REAL      :: radvel_param(2)  ! Parmeters used to specify the radial velocty and distribution

!...     Materials data to load:
         INTEGER   :: mat_opt
         CHARACTER :: mat_file*512
      ENDTYPE type_options_osm


      TYPE, PUBLIC :: type_options_filament
         REAL*4    :: version = 1.0

         INTEGER   :: opt
         INTEGER   :: clip
         INTEGER   :: target_flux
         REAL      :: start_time
         REAL      :: time_step
         REAL      :: scale(3)
         REAL      :: length1
         REAL      :: length2
! Time of interest for non-iterating solution, i.e. advancement to a particular time where everything starts at t=0.0
! Cross field dimention of filament grid refinement: 1 - initial pass, 2 - second pass, 3 - plasma assignment
      ENDTYPE type_options_filament

!     EIRENE options:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_options_eirene
         REAL*4    :: version = 1.0

!...     i/o:
         INTEGER   :: f_eirene_load     ! Load status of reference files
         CHARACTER :: f_eirene_dir*512  ! Directory where EIRENE data files are located 
         CHARACTER :: f_eirene_13*512   ! eirene.13
         CHARACTER :: f_eirene_15*512   ! eirene.15 for time dependent runs

!...     Options:    
         INTEGER   :: geom            ! geometry  2-triangles 3-tetrahedrons (toroidal)
         INTEGER   :: time            ! Execution time (CPU time in seconds)
         REAL      :: dtimv           ! ?
         REAL      :: time0           ! ?
         INTEGER   :: opacity         ! Lyman alpha photon opacity option
         INTEGER   :: photons         ! Photon transport option
         INTEGER   :: trim            ! Fast ion surface collision database
         INTEGER   :: bgk             ! Neutral viscosity
         INTEGER   :: niter           ! Number of Eirene self-iterations
         INTEGER   :: ntime           ! Number of time steps
         INTEGER   :: data            ! Eirene input file  1=internal, 2=external 
         INTEGER   :: ilspt           ! Sputering option
         INTEGER   :: whipe           ! Reduce the plasma density to very low values (testing mode)
!...     3D:
         INTEGER   :: tet_iliin       ! Reflection property for the toroidal boundary surfaces
!        Tetrahedral mesh generation:
         INTEGER       :: tet_n
         REAL          :: tet_type     (EIR_MAXNTET)
         REAL          :: tet_x1       (EIR_MAXNTET)  ! Crop boundary for the grid, if necessary (tet_type = 1.0)
         REAL          :: tet_y1       (EIR_MAXNTET)  !   lower inner point = (x1,y1)  (and for refinement)
         REAL          :: tet_z1       (EIR_MAXNTET)  !     (for refinement only)
         REAL          :: tet_x2       (EIR_MAXNTET)  !   upper outer point = (x2,y2)
         REAL          :: tet_y2       (EIR_MAXNTET)
         REAL          :: tet_z2       (EIR_MAXNTET)  
         INTEGER       :: tet_index    (EIR_MAXNTET)  ! Index of slice (tet_type=2.0 and 3.0)
                                                      ! Tet_type = 2.0, slices:
         INTEGER       :: tet_mode     (EIR_MAXNTET)  !   method for deciding the angular width of the slice
         REAL*8        :: tet_param1   (EIR_MAXNTET)  !   1st parameter used to specify the angular width
         REAL*8        :: tet_param2   (EIR_MAXNTET)  !   2nd parameter
         REAL*8        :: tet_param3   (EIR_MAXNTET)  !   3rd general purpose parameter
         CHARACTER*128 :: tet_del_hole (EIR_MAXNTET)  !   list of holes to apply to the slice (as listed in the additional surfaces)
         CHARACTER*128 :: tet_del_zone (EIR_MAXNTET)  !   list of zones to delete from slice (as specified in when setting the void grid)
                                                      ! Tet_type = 3.0, sectors:
         CHARACTER*128 :: tet_sec_list (EIR_MAXNTET)  !   list of slices to be included in this sector
                                                      ! Tet_type = 4.0, full grid composite:
         CHARACTER*128 :: tet_composite(EIR_MAXNTET)  ! List of sector indices comprising the full grid (tet_type = 4.0)
         REAL          :: tet_offset   (EIR_MAXNTET)  !   angular start location for the grid (-999.0 = symmetric about 0.0)

!...     Surface properties in EIRENE:
         INTEGER       :: sur_n       
         REAL          :: sur_type    (EIR_MAXNSUR)  ! Type: 1.0-non-default index, 1.1-stratum index, 2.0-standard DIVIMP wall index, 3.0-additional DIVIMP wall index
         CHARACTER*128 :: sur_index   (EIR_MAXNSUR)  ! Poloidal index 
         CHARACTER*128 :: sur_sector  (EIR_MAXNSUR)  ! Toroidal sector 
         INTEGER       :: sur_iliin   (EIR_MAXNSUR)  ! Surface transmission option
         INTEGER       :: sur_ilside  (EIR_MAXNSUR)  ! Surface orientation option
         INTEGER       :: sur_ilswch  (EIR_MAXNSUR)  ! Surface index switching option
         REAL          :: sur_tr1     (EIR_MAXNSUR)  ! Surface transparency 1
         REAL          :: sur_tr2     (EIR_MAXNSUR)  ! Surface transparency 2
         REAL          :: sur_recyct  (EIR_MAXNSUR)  ! Recycling fraction
         INTEGER       :: sur_ilspt   (EIR_MAXNSUR)  ! Sputtering option
         INTEGER       :: sur_temp    (EIR_MAXNSUR)  ! Over-ride of globally applied surface temperature
         CHARACTER*64  :: sur_mat     (EIR_MAXNSUR)  ! Material, i.e. C, Be, W, etc.
         REAL          :: sur_coverage(EIR_MAXNSUR)  ! defunct...
         INTEGER       :: sur_hard    (EIR_MAXNSUR)  ! Make sure the specified surface properties are isolated to this surface only
         INTEGER       :: sur_remap   (EIR_MAXNSUR)  ! For the evil remapping scheme, currently used to deal with perfectly conformal surfaces
         CHARACTER*256 :: sur_tag     (EIR_MAXNSUR)
!...     Particle energy spectra:
         INTEGER      :: nadspc
         INTEGER      :: ispsrf     (EIR_MAXNSPECTRA)  ! Surface index, <0=non-default standard, >0=additional surfaces
         CHARACTER*32 :: ispsrf_ref (EIR_MAXNSPECTRA)  ! Which code does the surface index refer to?
         INTEGER      :: iptyp      (EIR_MAXNSPECTRA)  ! Species type eg 1=atoms, 2=molecules, 3=test ions, 4=?
         INTEGER      :: ipsp       (EIR_MAXNSPECTRA)  ! Species sub-index eg, 1=first atom species, 2=second atom species, etc.
         INTEGER      :: isptyp     (EIR_MAXNSPECTRA)  ! Spectrum type wrt units, 1=1/eV/s, 2=1/s
         INTEGER      :: nsps       (EIR_MAXNSPECTRA)  ! Number of bins
         INTEGER      :: isrfcll    (EIR_MAXNSPECTRA)  ! Kind of spectrum, 0=surface flux
         INTEGER      :: idirec     (EIR_MAXNSPECTRA)  ! If >0 then a projection on a direction is used in the statistics (??)
         REAL         :: spcmn      (EIR_MAXNSPECTRA)  ! Lower bound of energy range for spectrum
         REAL         :: spcmx      (EIR_MAXNSPECTRA)  ! Upper bound
         REAL         :: spc_shift  (EIR_MAXNSPECTRA)  ! ??? for future use perhaps
         REAL         :: spcplt_x   (EIR_MAXNSPECTRA)  ! ???
         REAL         :: spcplt_y   (EIR_MAXNSPECTRA)  ! ???
         REAL         :: spcplt_same(EIR_MAXNSPECTRA)  ! ???
         REAL         :: spcvx      (EIR_MAXNSPECTRA)  ! x-direction for IDIREC >0
         REAL         :: spcvy      (EIR_MAXNSPECTRA)  ! y-direction
         REAL         :: spcvz      (EIR_MAXNSPECTRA)  ! z-direction
         REAL         :: spc_p1     (EIR_MAXNSPECTRA,3)  ! staring point of LOS (for getting spectra for all cells that intersect the LOS)
         REAL         :: spc_p2     (EIR_MAXNSPECTRA,3)  ! ending point
         REAL         :: spc_dist   (EIR_MAXNSPECTRA)    ! location / distance of the cell along the LOS
!...     Particle sources:
         REAL      :: alloc           ! Flux / npts weighting (0.0 = npts only, 1.0 = flux only)
!         REAL      :: puff_type   (EIR_MAXNPUFF)
!         INTEGER   :: puff_npts   (EIR_MAXNPUFF)        ! *** IN USE? ***
!         REAL      :: puff_flux   (EIR_MAXNPUFF)
!         REAL      :: puff_frac   (EIR_MAXNPUFF)
!         INTEGER   :: puff_species(EIR_MAXNPUFF)
!         INTEGER   :: puff_index  (EIR_MAXNPUFF)
!         INTEGER   :: puff_energy (EIR_MAXNPUFF)
!         INTEGER   :: puff_cosine_ind(EIR_MAXNPUFF)
!         INTEGER   :: puff_cosine_max(EIR_MAXNPUFF)
!         INTEGER   :: puff_pos (3,EIR_MAXNPUFF)
!         INTEGER   :: puff_vec (3,EIR_MAXNPUFF)
!         CHARACTER :: puff_note(  EIR_MAXNPUFF)*1024

c...    Strata:
        INTEGER   :: nstrata 
        REAL      :: type           (EIR_MAXNSTRATA)
        INTEGER   :: Z              (EIR_MAXNSTRATA)           ! atomic number
        INTEGER   :: A              (EIR_MAXNSTRATA)           ! atomic mass
        INTEGER   :: npts           (EIR_MAXNSTRATA)
        REAL      :: flux           (EIR_MAXNSTRATA)
        REAL      :: flux_fraction  (EIR_MAXNSTRATA)
        CHARACTER*128 :: species_tag(EIR_MAXNSTRATA)
!         CHARACTER*128 :: species_tag*128(EIR_MAXNSTRATA)  ! gfortranH
        INTEGER   :: species        (EIR_MAXNSTRATA)
        INTEGER   :: species_index  (EIR_MAXNSTRATA)
        REAL      :: energy         (EIR_MAXNSTRATA)
        INTEGER   :: target         (EIR_MAXNSTRATA)
        INTEGER   :: range_cell   (2,EIR_MAXNSTRATA)
        INTEGER   :: range_tube   (2,EIR_MAXNSTRATA)
!        REAL      :: cos         
!        REAL      :: cos_max     
!        CHARACTER :: note*512    
!        INTEGER   :: indsrc  ! ...
        CHARACTER*512 :: txtsou  (EIR_MAXNSTRATA)
!         CHARACTER :: txtsou*512   (EIR_MAXNSTRATA)  ! gfortran
!        INTEGER   :: ninitl       
!        INTEGER   :: nemods       
!        CHARACTER :: species_tag*4
!        INTEGER   :: nspez
!        CHARACTER :: distrib*5
!        INTEGER   :: inum
!        INTEGER   :: indim
!        INTEGER   :: insor
!        REAL      :: sorwgt
        REAL      :: sorlim       (EIR_MAXNSTRATA)
        REAL      :: sorind       (EIR_MAXNSTRATA)
!        INTEGER   :: nrsor    
!        INTEGER   :: nasor
        REAL      :: sorad      (6,EIR_MAXNSTRATA)
        REAL      :: sorene       (EIR_MAXNSTRATA)
!        REAL      :: soreni   
        REAL      :: sorcos       (EIR_MAXNSTRATA)
        REAL      :: sormax       (EIR_MAXNSTRATA)

!...     Geometry / structure:
         INTEGER   :: ntorseg         ! Number of toroidal segments for regular toroidal descretization
         REAL      :: torfrac         ! Fraction of torus included in the simulation 
         INTEGER   :: ntri            ! Number of triangle mesh information lines
         REAL      :: tri             ! Triangle mesh specification
         INTEGER   :: mat1            ! Target material?
         INTEGER   :: mat2            ! Wall material?            
         REAL      :: ctargt          ! Target surface temperature
         REAL      :: cwallt          ! Wall surface temperature

         REAL          :: add_version
         INTEGER       :: nadd
         INTEGER       :: add_type    (EIR_MAXNADD)
         INTEGER       :: add_index   (EIR_MAXNADD)
         CHARACTER*128 :: add_file    (EIR_MAXNADD)
         CHARACTER*128 :: add_file_tag(EIR_MAXNADD)
         CHARACTER*128 :: add_tag     (EIR_MAXNADD)
         REAL          :: add_holex   (EIR_MAXNADD)
         REAL          :: add_holey   (EIR_MAXNADD)

!...     Voids between the fluid grid and the wall:
         REAL          :: void_version
         INTEGER       :: nvoid 
         INTEGER       :: void_zone(  EIR_MAXNVOID)
         INTEGER       :: void_grid(2,EIR_MAXNVOID)
         INTEGER       :: void_wall(2,EIR_MAXNVOID)
         INTEGER       :: void_add (2,EIR_MAXNVOID)
         REAL          :: void_res (  EIR_MAXNVOID)
         REAL          :: void_hole(2,EIR_MAXNVOID)
         INTEGER       :: void_code(  EIR_MAXNVOID)
         REAL          :: void_ne  (  EIR_MAXNVOID)
         REAL          :: void_te  (  EIR_MAXNVOID)
         REAL          :: void_ti  (  EIR_MAXNVOID)
         CHARACTER*128 :: void_tag (EIR_MAXNVOID)         

         CHARACTER*128 :: void2_grid(EIR_MAXNVOID)
         CHARACTER*128 :: void2_wall(EIR_MAXNVOID)
         CHARACTER*128 :: void2_add (EIR_MAXNVOID)
         CHARACTER*128 :: void2_hole(EIR_MAXNVOID)

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
         REAL    :: epot
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
         INTEGER   :: file_format
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
         INTEGER   :: fit_width
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
        REAL    :: version_wall
        REAL    :: version_species
        INTEGER :: n                 ! = NTUBE
        INTEGER :: isep
        INTEGER :: isep2
        INTEGER :: ipfz
        INTEGER :: ikto
        INTEGER :: ikti
        REAL*8  :: r0
        REAL*8  :: z0
        REAL*8  :: core_volume_total
        REAL*8  :: core_volume_grid
        REAL*8  :: core_surface_area
        INTEGER :: nxpt
        REAL*8  :: rxpt(2)
        REAL*8  :: zxpt(2)
        INTEGER :: ixpt(2,2)
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
        INTEGER :: type                  ! 'idring' in DIVIMP -- *** useful? ***
!...    Geometry data:
        INTEGER :: n                     ! Number of knots/cells in the flux-tube
        INTEGER :: ikti                       
        INTEGER :: ikto                  ! etc...              
        REAL    :: smax
        REAL    :: pmax
        REAL    :: rho
        REAL    :: psin
        REAL    :: metric(2)             ! Cross-field othogonality metric at ends of tube
        REAL*8  :: volume
        REAL*8  :: surface_area(2)
        REAL*8  :: target_area(2)
!...    3D:
        REAL    :: dangle                ! Toroidal angle step when calculating 3D flux-tube
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
        REAL    :: gamma(2,s28_MAXNION)       
        REAL    :: qe(2,s28_MAXNION)
        REAL    :: te_upstream(2,s28_MAXNION)
        REAL    :: Psol(2,s28_MAXNION)
        REAL    :: efield(2)               ! Electric field strength
        REAL    :: epot  (2)               ! Electrostatic potential
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

      TYPE, PUBLIC :: type_tube2    ! *** NEW TUBE VARIABLES THAT DON'T NEED TO BE SAVED, SO NOT UPDATING MAIN tube SPECIFICATION AT THE MOMENT ***
        INTEGER*4 :: state         ! General info on the "solver state" of the tube, or anything else for that matter
!                           BIT 0 - 1-default symmetry point applied
!                           BIT 1 - 1-solution for ring has been successfully calculated
!                           BIT 2 - 1-node linked to an invalid ring (solution hadn't been calculated yet)
        INTEGER*2 :: target_pe(2)  ! Dynamically specify the target jsat based on the upstream electron pressure
      ENDTYPE type_tube2
!
!     Cells:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_cell
        INTEGER*2 :: ik                             ! Knot index of cell in standard 2D grid
        INTEGER*2 :: ir                             ! Ring index
!...    Geometry data:
        REAL    :: cencar(3)                        ! Cell caresian center in machine coordinates
        REAL    :: centor(3)                        ! Cell toroidal center in 'toroidal' coordinates
        REAL    :: vol                              ! Cell volume 
        INTEGER :: nside                            ! *** DELETE? *** Number of sides for each cell: nside=6 usually
!...    3D:
        REAL    :: s_3D
        REAL    :: sbnd_3D(2)
        REAL    :: area_3D
        REAL    :: volume_3D
!         ...need some reference to location relative to x-point...
        REAL    :: s                                ! Distance of cell center along the magnetic field line (m), s=0 at inner target (LO index target)
        REAL    :: p                                ! Poloidal distance of cell center from inner target (LO index target)
        REAL    :: sbnd(2)                          ! 's' at ends of each cell
        REAL    :: pbnd(2)                          ! 'p' at ends of each cell
        REAL    :: ds                               ! *** delete? ***
        REAL    :: dp                               ! *** delete? ***
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
!
!...  Kinetic ion code results (Eirene/DIVIMP):
!     --------------------------------------------------
      TYPE, PUBLIC :: type_kinetic
        REAL :: te
      ENDTYPE type_kinetic
!
!...  Wall segments:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_wall
        REAL      :: type                 ! Don't know why I used INDEX here, should replace
        INTEGER   :: class                ! everything with a call to the actual name of 
        INTEGER   :: index(10)            ! the variable, i.e. %index(WAL_TUBE) should be
        CHARACTER :: material_tag*128     ! %tube
        INTEGER   :: material
        REAL      :: temperature
        REAL*8    :: v1(2)
        REAL*8    :: v2(2)
        CHARACTER :: file_name*512
        INTEGER   :: file_format
      ENDTYPE type_wall
!
!...  Plasma species:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_species
        REAL      :: version
        CHARACTER :: tag*128
        INTEGER   :: index
        INTEGER   :: type
        CHARACTER :: code_source_tag*128
        INTEGER   :: code_source
        CHARACTER :: code_transport_tag*128
        INTEGER   :: code_transport
        INTEGER   :: maxcharge
        INTEGER   :: nhistory
        INTEGER   :: ntime
      ENDTYPE type_species
!
!...  Simulation code designation:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_code
        REAL      :: version
        CHARACTER :: name*128
        CHARACTER :: home*128
        CHARACTER :: run_script*128
      ENDTYPE type_code
!
!...  Material data:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_material
        REAL      :: version
        INTEGER   :: type
        CHARACTER :: tag*128
        INTEGER   :: A
        INTEGER   :: Z
        REAL      :: mass
        REAL      :: density
        CHARACTER :: crystal_structure*128
        REAL      :: monolayer_thickness
        REAL      :: melting_point
        REAL      :: thermal_conductivity
        INTEGER   :: adas_divimp
        INTEGER   :: adas_eirene
      ENDTYPE type_material
!
!     Target plates:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_target
        REAL      :: version
        INTEGER   :: type         ! ...
        INTEGER   :: location     ! Location relative to the divertor/limiter layout
        INTEGER   :: position     ! Which end of the flux-tube: low index (1, 2, 3...) or high (...n-2, n-1, n) 
        CHARACTER :: tag*1024     
        INTEGER   :: nlist        ! Number of tubes associated with this target
        INTEGER   :: npart        ! Number of parts that the target is made up of if discontinuous
        INTEGER   :: ipeak        ! Index of nominal distribution peak, i.e. usually ISEP or ISEP2
        REAL      :: sep_psin     ! PSIn at the local target/separatrix intersection             
        REAL      :: sep_rho      ! RHO (midplane mapping) at the local target/sep. intersection
        INTEGER   :: ilist(1024)  ! List of the tubes
        INTEGER   :: ipart(1024)  ! Part index
      ENDTYPE type_target
!
!...  
      REAL, PUBLIC, ALLOCATABLE, SAVE :: tmpflx(:,:)

      END MODULE mod_sol28
!
! ----------------------------------------------------------------------
!
      MODULE mod_sol28_targets
      USE mod_sol28      
      IMPLICIT none

      INTEGER :: ntarget
      TYPE(type_target), ALLOCATABLE :: target(:)

      END MODULE mod_sol28_targets
!
! ----------------------------------------------------------------------
!
      MODULE mod_options
      USE mod_sol28      
      IMPLICIT none

      SAVE

      TYPE(type_options_eirene  ) :: opt_eir   ! EIRENE options
      TYPE(type_options_filament) :: opt_fil   ! Filament options

      END MODULE mod_options
!
! ----------------------------------------------------------------------
!
      MODULE mod_sol28_wall
      USE mod_sol28      
      IMPLICIT none

      INTEGER, PARAMETER :: WAL_GROUP  = 1,  ! Group index
     .                      WAL_INDEX  = 2,  ! Index within group
     .                      WAL_TUBE   = 3,  ! Corresponding tube for target segments
     .                      WAL_TARGET = 4,  ! Target index (LO,HI)
     .                      WAL_RANGE1 = 8,  ! Target index (LO,HI)
     .                      WAL_RANGE2 = 9   ! Target index (LO,HI)

      INTEGER, SAVE :: nopt_wall
      TYPE(type_wall), SAVE :: opt_wall(100)

      INTEGER, SAVE :: nwall
      TYPE(type_wall), ALLOCATABLE :: wall(:)

      END MODULE mod_sol28_wall
!
! ----------------------------------------------------------------------
!
      MODULE mod_sol28_reference
      USE mod_sol28      
      USE mod_geometry
      USE mod_options
      IMPLICIT none

      PUBLIC

!...  Reference plasma:
      INTEGER, SAVE :: ref_ntube,ref_nion,ref_ncell,ref_nfluid
      TYPE(type_tube ), ALLOCATABLE       :: ref_tube (:)
      TYPE(type_cell ), ALLOCATABLE, SAVE :: ref_cell (:) 
      TYPE(type_fluid), ALLOCATABLE, SAVE :: ref_fluid(:,:) 

      END MODULE mod_sol28_reference
!
! ----------------------------------------------------------------------
!
      MODULE mod_sol28_global
      USE mod_sol28      
      USE mod_sol28_reference
      USE mod_geometry
      USE mod_options
      IMPLICIT none

      PUBLIC

      INTEGER, PARAMETER :: IND_IK      = 1,  ! Object index parameters
     .                      IND_IR      = 2,  
     .                      IND_IS      = 3,  
     .                      IND_ZONE    = 4,  
     .                      IND_FLUID   = 5,  ! *** NEED TO ADD A FLAG IN THE OBJECT BEING CREATED
     .                      IND_KINETIC = 6,  !     THAT IDENTIFIES WHERE IT WAS GENERATED, SINCE
     .                      IND_NEUTRAL = 7,  !     THESE MAPPINGS ARE SOURCE DEPENDENT, I.E. THEY ARE DIFFERENCE IN mod_sol28 and mod_eirene06...
     .                      IND_FIELD   = 8,  ! Vaccum zone outside standard grid, from external call to TRIANGLE
     .                      IND_CELL    = 9,
     .                      IND_OBJECT  = OBJ_MAXNINDEX+1  ! Just for cell and tube finding in GetCell and and GetTube

      INTEGER, PARAMETER :: IND_STDGRD  = 1,  ! Magnetic fluid grid side index, i.e. 12, 23, 34, 41
     .                      IND_TARGET  = 2,  ! Target (block 7 stratum in Eirene input file)
     .                      IND_SURFACE = 3,  ! Surface (block 2A non-default surface in Eirene input file)
     .                      IND_WALL    = 4   ! Surface (block 2A non-default surface in Eirene input file)


!...
      TYPE(type_options_osm), SAVE :: opt
      INTEGER, PARAMETER :: MAX_NOPT = 10
      INTEGER nopt
      TYPE(type_options_osm), SAVE :: opt_iteration(1:MAX_NOPT)
      INTEGER logop, logfp

!...  Nasty...
      INTEGER, PARAMETER :: MAX_NITERATION = 1000
      INTEGER               iiteration,niteration
      CHARACTER*1024        iteration_buffer(MAX_NITERATION)

      INTEGER osmnnode    
      TYPE(type_node) :: osmnode(100)

      INTEGER store_sopt(1000),store_mnode(1000),store_nnode(1000)
      TYPE(type_node) store_node(20,1000)
!...  
      TYPE(type_grid), SAVE :: grid
      INTEGER, SAVE :: ntube
      TYPE(type_tube ), ALLOCATABLE, SAVE :: tube (:) 
      TYPE(type_tube2), ALLOCATABLE, SAVE :: tube2(:) 

!      INTEGER*4, SAVE, ALLOCATABLE :: tube_state(:),  ! move into the TUBE array eventually, but only local use for now...
!       _state BIT 0 - 1-default symmetry point applied
!              BIT 1 - 1-solution for ring has been successfully calculated
!              BIT 2 - 1-node linked to an invalid ring (solution hadn't been calculated yet)

      INTEGER, SAVE :: nion,ncell,nfield,npin,nphoton,nfluid,nkinetic,
     .                 nimpurity,ndrift
      TYPE(type_cell    ), ALLOCATABLE, SAVE :: cell    (:) 
      TYPE(type_field   ), ALLOCATABLE, SAVE :: field   (:) 
      TYPE(type_neutral ), ALLOCATABLE, SAVE :: pin     (:,:) 
      TYPE(type_photon  ), ALLOCATABLE, SAVE :: photon  (:,:) 
      TYPE(type_kinetic ), ALLOCATABLE, SAVE :: kinetic (:,:) 
      TYPE(type_fluid   ), ALLOCATABLE, SAVE :: fluid   (:,:)
      TYPE(type_impurity), ALLOCATABLE, SAVE :: impurity(:,:) 
      TYPE(type_drift   ), ALLOCATABLE, SAVE :: drift   (:,:) 


!...  Plasma species that are being tracked in the simulation:
      INTEGER, SAVE :: nspecies
      TYPE(type_species), SAVE :: species(100)

!...  Material data:
      INTEGER, SAVE :: nmaterial
      TYPE(type_material), SAVE :: material(100)


!...  Index mapping in GetObject and GetTube:
      INTEGER, ALLOCATABLE, SAVE :: obj_index_map (:,:),
     .                              tube_index_map(:,:)

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


      REAL*8, PARAMETER :: ECH = 1.6022D-19, AMU = 1.67D-27


      TYPE(type_options_osm), SAVE :: opt
      INTEGER :: logop,logfp

      REAL*8, SAVE :: chisq(0:7),t_chisq(0:7),m_chisq(0:7)


!...  Local solver options:
      INTEGER, SAVE ::    
c     .  opt_bc(2),             ! Method of assigning boundary conditions (1=target, 2=upstream)
c     .  opt_p_ion(2),             
     .  sol_option,
     .  opt_p_ion_scale(2),
     .  node_par_mode  (100)
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

      TYPE(type_tube   ), SAVE :: tube
      TYPE(type_cell   ), SAVE :: cell (S28_MAXNKS)
      TYPE(type_neutral), SAVE :: pin  (S28_MAXNKS,S28_MAXNION) 
      TYPE(type_fluid  ), SAVE :: fluid(S28_MAXNKS,S28_MAXNION) 

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
     .  area  (0:S28_MAXNKS+1),   
     .  vol   (  S28_MAXNKS)  ,   
!...    Plasma fluid quantities:
     .  isat  (0:S28_MAXNKS+1,0:S28_MAXNION),  
     .  ne    (0:S28_MAXNKS+1),                 
     .  ni    (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  vi    (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  pe    (0:S28_MAXNKS+1),                 
     .  pi    (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  te    (0:S28_MAXNKS+1),                 
     .  ti    (0:S28_MAXNKS+1,0:S28_MAXNION),   
     .  qe    (0:S28_MAXNKS+1),                 
     .  qi    (0:S28_MAXNKS+1,0:S28_MAXNION),
     .  qcond (0:S28_MAXNKS+1),
     .  qconv (0:S28_MAXNKS+1),
     .  efield(0:S28_MAXNKS+1),                ! Electric field strength
     .  epot  (0:S28_MAXNKS+1),                ! Electrostatic potential
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
