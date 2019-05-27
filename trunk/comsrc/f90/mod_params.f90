module mod_params
  use debug_options
  implicit none

  public
  
  !
  !     params: this file contains the maximum array size parameters
  !             and other constants relevant to the divimp code.
  !
  !     -*-fortran-*-
  integer,public :: maxnks,maxnrs,maxnds,maxngs,maxizs,maxins,maximp,isect,maxnts,maxnoc,&
       maxnws,maxnxs,maxnys,maxvmf,maxthe,maxsn,maxpts,maxplrp,inimout,ipinout,msolpt,&
       maxads,maxgxs,maxgys,maxch3,maxixs,maxiys,maxseg,maxplts,maxnfla,maxpiniter,&
       mbufle,mbufx,mves
  !    >        ,iplot
  !
  !
  !       maxgxs and maxgys replace the usage of maxnxs and maxnys in
  !       the plotting routines. the values for these must always
  !       be greater than or equal to their nxs,nys counterparts.
  !       the values of maxixs and maxiys control the space allocated
  !       for the ifxys, ikxys and irxys arrays - if these are not
  !       being used in divimp or out then the values can be set small.
  !       otherwise, they should be set to exactly equal the values of
  !       maxgxs and maxgys.
  !
  !       maxseg is the maximum number of segments in the wall and pump.
  !       this should be set equal to mvesm in pincoms/p01.(or the passing
  !       of flxhw2 will not work)!
  !
  !
  !     support for tungsten and more fluids ...
  !
  !     ipp/01 krieger: change maxizs from 18 to 74 (tungsten)
  !     ipp/01 krieger: change maxnfla from 7 to 19 (h+he+c+ne)
  !     ipp/01 krieger: change maxe2dizs from 6 to 19 (as above)
  !
  !     parameters for divimp and out
  !
  integer,public :: maxe2dizs,ipindat
  real,public :: hi,lo,root2,pi,raddeg,emi,degrad,ech,amu,machhi,machlo,kboltz,cspeed,&
       eps0
  
  character*5,public :: verson
  !     >  maxnks=300  ,maxnrs=200  ,maxnoc=100 ,maxnds=2*maxnrs,maxnys=1,  ! iter
  !     >  maxnks=200  ,maxnrs=100  ,maxnoc=100 ,maxnds=2*maxnrs,maxnys=1,
  !     >  maxnks=200  ,maxnrs=100  ,maxnoc=100 ,maxnds=110   ,maxnys=1,
  !     >  maxngs=40   ,maxizs=74  ,maxins=maxnrs ,maximp=1000000,maxvmf=5,  ! iter
  !     >  maxngs=50   ,maxizs=74   ,maxins=150 ,maximp=100000,maxvmf=5,
  !     >  maxpts=500  ,maxplrp=12  ,msolpt=100 ,maxads=80    ,maxch3=300,
  !
  !     >  maxnks=200  ,maxnrs=400  ,maxnoc=100 ,maxnds=2*maxnrs,maxnys=1,
  !     >  maxngs=50   ,maxizs=74   ,maxins=150 ,maximp=100000,maxvmf=5,
  !     >  maxpts=1000 ,maxplrp=12  ,msolpt=100 ,maxads=80    ,maxch3=300,
  !     >  maxgxs=401  ,maxgys=400  ,maxixs=1   ,maxiys=1     ,ipinout=16,
  !     >  isect =128  ,maxthe=5000 ,maxsn=5000 ,maxseg=1000  ,inimout=37,
  ! slmod end
  !

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

  !parameter (maxnks=260  ,maxnrs=190  ,maximp=10000000,  maxizs=74)   ! extended high res grid, large particles, large max charge state (tungsten)
  parameter (maxnks=200  ,maxnrs=100  ,maximp=500000,  maxizs=74)      ! normal resolution grid, medium particles, tungsten
  !parameter (maxnks=200  ,maxnrs=100  ,maximp=500000,  maxizs=6)      ! normal resolution grid, medium particles, carbon

  parameter (verson='6a/53'  ,maxnts=1   ,maxnws=10000 ,maxnxs=1,&
       maxnoc=100 ,maxnds=2*maxnrs,maxnys=1, maxngs=40   ,&
       maxins=maxnrs ,maxvmf=5, maxpts=750  ,maxplrp=12  ,&
       msolpt=100 ,maxads=80    ,maxch3=300,  maxgxs=201  ,maxgys=200  ,maxixs=1   ,&
       maxiys=1     ,ipinout=16,isect =128  ,maxthe=5000 ,maxsn=5000 ,maxseg=1000  ,&
       inimout=37,maxplts=36  ,maxnfla=21  ,maxpiniter=500           ,ipindat=15,mbufle =10  ,&
       mbufx  =10  ,mves=maxpts,root2 =1.414213562       ,pi=3.141592654,raddeg=57.29577952       ,&
       emi=1.602192e-19/1.672614e-27  ,degrad=1.745329252e-02   ,&
       ech=1.602192e-19,amu=1.672614e-27         ,kboltz=1.38e-23,hi=1.e37    ,lo=1.e-37   ,&
       machhi=1.0e37 ,machlo=1.0e-37  ,cspeed=2.998e8           ,eps0=  8.85e-12)
  !
  !     defining poygon sides between vertices.
  !
  parameter(maxe2dizs=21)
  integer,public :: inward41,up34,outward23,down12
  !
  !     parameters related to the out program
  !
  parameter(inward41 = 4, up34     = 3, outward23= 2, down12   = 1)
  integer,public :: maxdatx
  !
  !     fortran unit numbers for output
  !
  parameter (maxdatx = 5000)
  !
  integer,public :: datunit,htmlunit,pinunit,dbgunit,tmpunit,auxunit,exptunit,tranunit,&
       diagunit
  !    >        , iplot    = 49)
  !
  !     datunit initialization is done at the beginning of the rundiv
  !     module - datunit is in a common block at the moment becasue its
  !     value needs to change for a kludgy work around on reordering
  !     some printed data.
  !
  !     the following is a list of most of the unit numbers referred to
  !     by divimp and eirene - i am trying to centralize this data to
  !     make it easier to select unit numbers when adding output to the code.
  !     much of this information can also be found by looking through the
  !     rundiv and runeire scripts.
  !
  ! unit#  purpose/file name
  !
  !    4   equilibrium grid file
  !    5   standard input  - usually the case file name (.d6i)
  !    6   standard output - usually mapped to debug output (.lim)
  !    7   case print out file (.dat)
  !    8   divimp raw data file (.raw)
  !    9   divimp/out input echo file (.inp)
  !   11   fluid code plasma solution to be read in (if one is specified) ($3$4,$3$4.g80,&
  !     $3.g80))
  !   12   auxiliary data file accompanying the background plasma solution ($3$4.aux,&
  !     $3.aux)
  !   13   experimental data file associated with this shot/grid ($3$4.experimemt,$3.experiment)
  !   14   nimbus - link to fort.5 supplies nimbus namelist input
  !   15   nimbus output file (.pinout)
  !   16   eirene/nimbus file for data transfer to divimp (.pinraw)
  !   17   pin file from divimp supplying plasma background to eirene/nimbus (.pin)
  !   18   ninbus - reserved for pump input files (set in nimbin namelist)
  !   20   eirene99 - link to equilibrium grid file - unit #4
  !  *21   eirene99 - link to trim.dat    ***    divimp - additional sol22 output (.sol22)
  !   22   divimp - html version of the data file (.html)
  !   24   divimp - print-out of pin data (eirene99 or nimbus) (.pinprn)
  !   26   out - column formatted print out of specified xy type plots (.grp)
  !   27   out - signal formatted file for plot outputs
  !   28   divimp/nimbus - jet hybrid wall data - link to shots/hydrid.dat
  !   30   nimbus - nimbus archive file - this is presently turned off
  !   31   eirene99 - link to divimp plasma input - unit #17
  !   32   eirene99 - eirene output file (.eir)
  !   35   nimbus - nimbus input file passed from linkpg (.pinmc)
  !   36   nimbus - nimbus punch file passed back to linkpg
  !   37   nimbus - nimbus print file - scanned by divimp for some data (.pinnim)
  !   41   tran file - data for the jet post processor routines
  !   46   diag - file containing divimp diagnostic and debug information
  !   49   out - additional plot output information - including result for single point data (.plt)
  !  *50   eirene99 - eirene input file    ***   divimp - debugging output file (contents?) (.debug)
  !   52   divimp/eirene99 - eirene geometry file (.eirgeo) - sometimes linked to upgrade.geom
  !   56   out - text print out of reciprocating probe data for specific plots (.probe)
  !   57   out - text print out of erosion/deposition data
  !   60   nimbus - reserved in nimbus for reading adas data
  !   62   divimp - divimp internal format background plasma if requested (.bgp)
  !   71   divimp - additional output from sol23 (.sol23)
  !   73   divimp - excel spreadsheet formatted data from sol23 (.exl23)
  !   74   divimp - link to cfd solution from previous divimp run (<oldcase>.cfd)
  !   75   divimp - cfd plasma solution from current sol23 run (.cfd)
  !   79   divimp - ??? - (.raw.rel)
  !   80   eirene99 - eirene particle tracks (.eirtrc)
  !   81   divimp - input file to eirene generated by divimp - linked to fort.50 by eirene99
  !   85   divimp - ??? - (.g1)
  !   87   divimp - ??? - (.g3)
  !   88   divimp - ??? - (.src)
  !   89   divimp - ??? - (.raw.src)
  !   94   divimp - ??? - (.raw.pla)
  !   95   divimp - ??? - (.raw.geo)
  !   99   eirene99 - temporary file used to add pressure guage data to fort.32
  !
  !
  !      data datunit /7/
  !
  parameter(dbgunit =  6, exptunit= 13,pinunit = 19,htmlunit= 22, tmpunit = 23,tranunit= 41,&
        diagunit= 46, auxunit = 61)
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

end module mod_params
