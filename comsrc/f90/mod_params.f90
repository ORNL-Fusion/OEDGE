module mod_params
  use debug_options
  use mod_io_units
  implicit none

  public
  
  !
  !     params: this file contains the maximum array size parameters
  !             and other constants relevant to the divimp code.
  !
  !
  !     parameters for divimp and out
  !
  
  character*5,public :: verson

  ! Memory usage by OEDGE is controlled by the values set in this file until the
  ! time when we transition to using true dynamic allocation with the array sizes
  ! completely specified at run time. Even if we move to such a system the default
  ! values in this file will likely still be relevant for backward compatibility with
  ! cases where the required sizes can't be easily determined from the inputs.
  !
  ! Key values that strongly affect memory utilization are:
  !  - MAXNRS - maximum number of flux surfaces or rings on the computational mesh
  !  - MAXNKS - maximum number of cells or knots along each flux surface or ring
  !  - MAXIMP - the maximum number of impurity particles to be followed
  !  - MAXIZS - the maximum charge state of the impurity species.
  !
  !  - MAXNRS and MAXNKS have to be set large enough for the computational mesh
  !  - MAXIMP needs to be equal to or larger than the number of particles to be followed
  !    in the simulation
  !  - MAXIZS needs to be equal to or greater than the atomic number of the specified impurity species
  !
  ! Key memory utilization parameters

  
  integer,public :: maxnks,maxnrs,maxnds
  integer,public :: maxngs,maxizs,maxins,maximp,isect,maxnts,maxnoc,&
       maxnws,maxnxs,maxnys,maxvmf,maxthe,maxsn,maxpts,maxplrp,msolpt,&
       maxgxs,maxgys,maxch3,maxixs,maxiys,maxseg,maxplts,maxnfla,maxpiniter,&
       mbufle,mbufx,mves
  integer,public :: maxe2dizs
  integer,public :: maxrtnsd,maxvizs

  
  real,public :: hi,lo,root2,pi,raddeg,emi,degrad,ech,amu,machhi,machlo,kboltz,cspeed,&
       eps0
  real, public :: emi_sqrt, larmor_const
  parameter (verson='6a/55' ,root2 =1.414213562, pi=3.141592654, raddeg=57.29577952 ,&
       degrad=1.745329252e-02   ,&
       hi=1.e37    ,lo=1.e-37   ,&
       machhi=1.0e37 ,machlo=1.0e-37  ,cspeed=2.998e8           ,eps0=  8.85e-12)

  parameter(ech=1.602192e-19,amu=1.672614e-27, emi=ech/amu   ,kboltz=1.38e-23)
  
  !
  !     defining poygon sides between vertices.
  !

  integer,public :: inward41,up34,outward23,down12
  !
  !     parameters related to the out program
  !
  parameter(inward41 = 4, up34     = 3, outward23= 2, down12   = 1)
  integer,public :: maxdatx
  !
  !     fortran unit numbers for output
  !
  !
  !    >        , iplot    = 49)
  ! common /outunits/ datunit
  !
  ! slmod begin - temporary
  ! save /outunits/
  integer,parameter,public :: linear_grid = 6        ! tag for uls grids
  integer,parameter,public :: gen_grid    = 7        ! for dave...
  integer,parameter,public :: osm_grid    = 9        ! import geometry from osm
  
  !...  some global options of convenience:
  !     jdemod - toggle this off unless there is a need for it
  !      logical, parameter :: sloutput  = .true.,   ! selects sl screen output
  integer,parameter,public :: ribbon_grid = 8
  ! slmod end
  logical,parameter,public :: sloutput  = .false., ippchange = .true.    ! selects recent ipp garching updates


  public:: initialize_parameters
  
  contains


    subroutine initialize_parameters
      implicit none

      ! This routine is used to initialize paramter values to default values for parameters that could be changed in the input file. This is a step
      ! towards support of full dynamic allocation of some storage to allow for a more flexible build of OEDGE that only utilizes as much memory as needed
      !
      ! This may also be used to define general grid parameters for small, medium and large grids if it is decided to proceed with an intermediate
      ! dynamic implementation for grids. 

      ! Most likely paramters to be present in the input file - grid limits - this may later be replaced by dynamic analysis
      ! of the grid using get_grid_parameters (works only for Sonnet/Carre grids so far.

      maxnrs = 190   ! Z01 : Max number of rings on the grid
      maxnks = 260   ! Z02 : Max number of knots/ring on the grid
      maxpts = 750   ! Z03 : generic variable for max number of wall data/SOL23 info/ADAS variables/etc
      maxseg = 1000  ! Z04 : Max number of wall segments (used for wall definition and wall fluxes)
      maxnws = 10000 ! Z05 : number of walk steps recorded
      msolpt = 100   ! Z06 : Number of points for detailed SOL backgrounds - soledge

      ! Final values derived from case input for tagged input - but need to default to large values for untagged input files
      
      maxpiniter = 25  ! Max number of SOL/PIN iterations
      maximp = 500000  ! Max number of impurities  NIMPS+NIMPS2 (NIMPS,NIMPS2 default:100,0)
      maxizs =  74     ! Max number of impurity charge states (use carbon as default - others will need to be specified)
      maxnts =   1     ! Max number of time slices in time dependent simulation
      maxnfla=  21     ! Max number of fluids in B2 file

      !
      ! These parameters are not currently dynamically definable since the use cases are limited.
      ! It will be easy to add these as needed.
      ! 
      ! Old features that generally won't need to be changed

      maxvmf = 5       ! Max number of VMF blocks (VMF option) - code not used so not worth updating at the present time
      maxplrp= 12      ! OLD code - maximum number of "Particular Line Radiation Profiles"
      isect  = 128     ! Size of blocks of random numbers - choice of 128 likely related to storage or
                       ! number representation in older computers
      mbufle =10  ! Related to code for adding a Baffle on JET grids - not used in decades
      mbufx  =10  ! Related to code for adding a Baffle on JET grids - not used in decades
     
      ! These parameters are unused in DIVIMP - they were used in OUT but have since been deprecated
      maxixs=1     
      maxiys=1
      maxnys=1
      maxnxs=1

      ! used in OUT 
      maxngs=40
      maxch3=300
      maxgxs=201       !MAXGXS/MAXGYS are mostly related to OUT but are used in DIVIMP to modify RMIN/RMAX/ZMIN/ZMAX
      maxgys=200
      maxthe=5000
      maxplts=36
      maxdatx= 10000 ! Z07 : Max number of experimental data points to be loaded - not used in DIVIMP

      ! not used in DIVIMP or OUT currently
      maxnoc=100
      maxsn=5000
      
      
    end subroutine initialize_parameters

    subroutine set_dependent_parameters
      implicit none
      ! This code sets parameters that are derived from other parametric values - splitting this allows
      ! revised values of the parameters to be read from the input file.
      ! Note that most of this duplication is historical and a code clean up could remove most of them
      
      maxvizs = maxizs
      maxrtnsd= maxizs+1
      maxnds = 2 * maxnrs
      maxins=maxnrs
      mves=maxpts   ! Alias for MAXPTS: Max number of wall segments (used for wall definition and wall fluxes) 
      maxe2dizs = maxnfla ! Max number of possible impurity ionization charge states in fluid data - set to number of fluids
      
      emi_sqrt = sqrt(emi)

      ! Restructure to avoid arithmetic underflow in calculation
      !larmor_const = sqrt(2.0 * ech * amu)/ech

      larmor_const = sqrt(2.0 * ech)/ech  * sqrt(amu)
      
    end subroutine set_dependent_parameters
    

    
  
end module mod_params
