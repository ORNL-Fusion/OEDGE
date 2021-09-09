module mod_params
  use mod_io_units
  implicit none

  public


  ! change from parameters to variables that can be read in so that LIM can
  ! be built only using the amount of storage needed for a specific case. 
  INTEGER::   MAXNXS = 100
  integer::   MAXNYS = 500
  !integer::   MAXNPS = 31
  integer::   MAXNPS = 20
  integer::   MAXIZS = 74
  integer::   MAXIMP = 500000000
  integer::   MAXQXS = 500
  integer::   MAXQYS = 5000
  INTEGER::   MAXY3D = 500
  integer::   MAXNTS = 1
  integer::   MAXINS = 100

  integer,parameter::   MAXNLS = 8

  integer::   maxpzone = 2   ! maximum number of different poloidal plasma zones (basically each zone has a different
                             ! plasma profile. Each poloidal plane (IP) in the simulation is associated with a poloidal plasma zone
                             ! This would be set greater than 1 if there is a collector probe simulation being run  
                             ! Pzone = 1 is the default plasma zone  
  integer,parameter::   ISECT  = 128
  integer,parameter::   maxput = 1000
  integer,parameter::   MAXOS = 500
  integer,parameter::   MAXLPD = 20
  integer,parameter::   MAXT = 100
  INTEGER,parameter::   MAXLEN = 100

  Logical:: big = .true.

  integer:: max_nsurf = 1 ! maximum number of poloidal elements to limiter surface
  
  ! log version number changes:
  ! 07 - added sngl(ddvs) and sdtimp to raw file and netcdf
  !
  CHARACTER, parameter:: VERSON*5 ='L3/10'
  
  !INTEGER::   MAXLEN,MAXADS
  ! maxads is now defined in mod_cadas2
  ! set default value to 100


  !LOGICAL::   BIG                                                             
  !PARAMETER (MAXNXS=100,  MAXNYS=500, MAXIZS=74,   MAXQYS=5000,   &            
  !     MAXQXS=500,  MAXNLS=8,   MAXNTS=1,    MAXIMP=100000000,&               
  !     MAXY3D=500,  MAXNPS=31,  MAXOS =500,  VERSON='L3/05',&           
  !     MAXINS=100,   ISECT =128, MAXPUT=1000, BIG=.TRUE.,   &
  !     MAXLPD=20,   MAXT=100,    MAXLEN=100,  &
  !     max_nsurf=1) 


  ! Fundamental constants
   REAL::      HI,LO,ROOT2,PI,RADDEG,EMI,DEGRAD,ECH,AMU,machhi,machlo
   parameter (&
       ROOT2 =1.414213562,       PI=3.141592654,                     &
       RADDEG=57.29577952,       EMI=1.602192E-19/1.672614E-27  ,    &
       DEGRAD=1.745329252E-02,   ECH=1.602192E-19,                   &
       AMU=1.672614E-27         ,                                    &
       HI=1.E37    ,LO=1.E-37   ,MACHHI=1.0E37 ,MACHLO=1.0E-37  )    

  
  !     Parameters related specifically to OUT
  !
  !integer  maxthe, maxdatx, maxngs,  maxpts

  !integer,parameter:: maxthe = 1000,  maxdatx = 1000, maxngs = 50, maxpts =  300
  integer,parameter:: maxthe = 1000,  maxdatx = 1000, maxngs = 50, maxpts =  300

  !
  !  Units are now defined in mod_io_units.f90
  !
  !     Some of the logical units used by the LIM code  
  !
  !integer datunit,dbgunit,outunit,exptunit
  !
  !parameter(dbgunit =  6, datunit = 7, exptunit=13, outunit=49)
  !
  !
  ! Unit#  Purpose/File name
  !	    
  !    5   Standard input  - usually the case file name (.d6i)
  !    6   Standard output - usually mapped to debug output (.lim)
  !    7   Case print out file (.dat)
  !    8   LIM raw data file (.raw)
  !    9   LIM/OUT input echo file (.inp)             
  !   13   Experimental Data file associated with this shot/grid ($3$4.experimemt,$3.experiment)
  !   26   OUT .grp output file
  !   49   OUT .plt output file
  ! 
  !
  !         include 'params'
end module mod_params
