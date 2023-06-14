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

      maximp = 500000

      maxizs = 74

      maxvizs = maxizs
      maxrtnsd=maxizs+1

      emi_sqrt = sqrt(emi)

      ! Restructure to avoid arithmetic underflow in calculation
      !larmor_const = sqrt(2.0 * ech * amu)/ech

      larmor_const = sqrt(2.0 * ech)/ech  * sqrt(amu)

      maxnrs = 190
      maxnks = 260
      maxnds = 2 * maxnrs
      
      maxnts=1
      maxnws=10000
      maxnxs=1

      maxnoc=100
      maxnys=1
      maxngs=40


      maxins=maxnrs
      maxvmf=5
      maxpts=750
      maxplrp=12
      msolpt=100
      maxch3=300
      maxgxs=201
      maxgys=200
      maxixs=1
     
      maxiys=1

      isect =128
      maxthe=5000
      maxsn=5000
      maxseg=1000
     
      maxplts=36
      maxnfla=21
      maxpiniter=500
      mbufle =10  
      mbufx  =10
      mves=maxpts

      maxe2dizs=21
      maxdatx = 5000
      
      
    end subroutine initialize_parameters


  
end module mod_params
