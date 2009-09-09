!     -*-Fortran-*-
c
c ======================================================================
c
      MODULE mod_legacy
      IMPLICIT none
      PUBLIC
     
      INTEGER, PARAMETER :: MAXNKS = 200, MAXNRS = 200, BOUNDARY = -1, 
     .                      MAXNAS = 100, MAXPUFF=  20, MAXNAS3  = 1000,
     .                      SLOUT = 50  , PINOUT = 88 , MAXASD   = 30,
     .                      IKLO = 1    , IKHI = 2    , MAXNAS2  = 200,
     .                      MAXPTS = 500

      REAL, PARAMETER :: PI = 3.141592, ECH = 1.6022E-19
 
      INTEGER irsep,irtrap,nrs,nlpdato,nlpdati,idring(MAXNRS),
     .        lpdatsw,relmode,rel_step
      LOGICAL connected
      REAL    psitarg(MAXNRS,2),lpdati(MAXNRS,4),lpdato(MAXNRS,4),
     .        te_mult_i,te_mult_o,ti_mult_i,ti_mult_o,n_mult_i,n_mult_o,
     .        relexpdat(2,3,MAXNRS)

c...  Target data:
      INTEGER tarninter(2)
      REAL    tarinter(MAXNRS,10,2),tarintermode(2)                     ! <- load from input file

c...  Grid data:
      INTEGER grdntreg(2),grdntseg(MAXNRS,2),grdtseg(MAXNRS,MAXNRS,2)   
   
      INTEGER osmns28                                                   ! <- load from input file
      REAL    osms28(MAXNRS,20)

c...  Eirene puff sources:
      INTEGER eirnpuff
      REAL    eirpuff(MAXNAS,MAXPUFF)

c...  Vessel structure data:
      INTEGER wallpts
      REAL    wallpt(MAXPTS,32)

      INTEGER eirnasdat
      REAL    eirasdat(MAXNAS2,MAXASD)

      INTEGER eirnspdat
      REAL    eirspdat(MAXNAS3,MAXASD)

      INTEGER nks(MAXNRS),ikins(MAXNKS,MAXNRS),irins(MAXNKS,MAXNRS)



      CONTAINS
c
c ----------------------------------------------------------------------
c
      SUBROUTINE InitializeLegacyVariables
      USE mod_sol28_global

      lpdatsw = 1
 
      rel_step  = 0
      relmode   = 0
      relexpdat = 0.0
       
      te_mult_i = 1.0
      ti_mult_i = 1.0
      n_mult_i  = 1.0
      te_mult_o = 1.0
      ti_mult_o = 1.0
      n_mult_o  = 1.0

      connected = .FALSE.

      tarninter = 0
      osmns28 = 0

      RETURN
      END SUBROUTINE
c
c ----------------------------------------------------------------------
c
      SUBROUTINE AssignLegacyVariables
      USE mod_sol28_global

      irsep  = grid%isep + 1
      irtrap = grid%ipfz + 2
      nrs    = grid%n    + 3

      psitarg(2       :irsep -1,1) = tube(1        :grid%isep-1)%psin
      psitarg(irsep   :irtrap-2,1) = tube(grid%isep:grid%ipfz-1)%psin
      psitarg(irtrap+1:nrs     ,1) = tube(grid%ipfz:grid%n     )%psin

      psitarg(1:MAXNRS,2) = psitarg(1:MAXNRS,1)

      idring(2       :irsep -1) = tube(1        :grid%isep-1)%type
      idring(irsep   :irtrap-2) = tube(grid%isep:grid%ipfz-1)%type
      idring(irtrap+1:nrs     ) = tube(grid%ipfz:grid%n     )%type

      END SUBROUTINE
c
c ----------------------------------------------------------------------
c
      SUBROUTINE LoadLegacyData(fname)
      IMPLICIT none

      CHARACTER*(*) fname
      INTEGER fp,i1,i2,PARAM1
      REAL    version

c...  Load data passed from DIVIMP:
      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
      READ (fp,ERR=98) version
      READ (fp,ERR=98) PARAM1
      IF (PARAM1.GT.MAXNRS) 
     .  CALL ER('LoadLegacyData','Increase MAXNRS',*97)
      READ (fp,ERR=98) (grdntreg(i1),grdntseg(1:PARAM1,i1),
     .                  (grdtseg(1:PARAM1,i2,i1),i2=1,PARAM1),i1=1,2)
      CLOSE (fp)
      
      RETURN
 97   WRITE(0,*) '  PARAM1,MAXNRS:',PARAM1,MAXNRS
      STOP
 98   CALL ER('LoadLegacy','Problem loading legacy data file',*99)
 99   STOP
      END SUBROUTINE
c
c ======================================================================
c
      END MODULE mod_legacy
