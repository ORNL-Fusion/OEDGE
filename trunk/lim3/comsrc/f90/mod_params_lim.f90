module mod_params
  implicit none

  public


  INTEGER::   MAXNXS,MAXNYS,MAXIZS,MAXNLS,MAXQXS,MAXQYS,MAXNTS,MAXIMP         
  INTEGER::   MAXY3D,MAXNPS,MAXINS,ISECT, MAXPUT,MAXOS ,MAXLPD, MAXT
  INTEGER::   MAXLEN,MAXADS
  REAL::      HI,LO,ROOT2,PI,RADDEG,EMI,DEGRAD,ECH,AMU,machhi,machlo
  LOGICAL::   BIG                                                             
  CHARACTER:: VERSON*5                                                        
  PARAMETER (MAXNXS=100,  MAXNYS=500, MAXIZS=74,   MAXQYS=5000,   &            
       MAXQXS=500,  MAXNLS=8,   MAXNTS=1,    MAXIMP=1000000,&               
       MAXY3D=500,  MAXNPS=31,  MAXOS =500,  VERSON='L3/04',&           
       MAXINS=100,   ISECT =128, MAXPUT=1000, BIG=.TRUE.,   &
       MAXLPD=20,   MAXT=10,    MAXLEN=100,  MAXADS=100,     &
       !
       !       Constants
       !
       ROOT2 =1.414213562,       PI=3.141592654,                     &
       RADDEG=57.29577952,       EMI=1.602192E-19/1.672614E-27  ,    &
       DEGRAD=1.745329252E-02,   ECH=1.602192E-19,                   &
       AMU=1.672614E-27         ,                                    &
       HI=1.E37    ,LO=1.E-37   ,MACHHI=1.0E37 ,MACHLO=1.0E-37  )    
  !
  !     Parameters related specifically to OUT
  !
  !integer  maxthe, maxdatx, maxngs,  maxpts
  integer,parameter:: maxthe = 1000,  maxdatx = 1000, maxngs = 50, maxpts =  300
  !
  !     Some of the logical units used by the LIM code  
  !
  integer datunit,dbgunit,outunit,exptunit
  !
  parameter(dbgunit =  6, datunit = 7, exptunit=13, outunit=49)
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
