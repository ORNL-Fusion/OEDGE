c     -*-Fortran-*-
c
      SUBROUTINE LIM3 (IMODE,NIZS,NIMPS,IMPADD,IMPCF,
     >                 QTIM,CPULIM,PRINPS,         
     >                 XWIDM,YWIDM,FSRATE,IONTIM,NEUTIM,SEED,IGEOM,             
     >                 NTBS,NTIBS,NNBS,NYMFS,NCVS,FACTA,FACTB,ITER,
c slmod
     >                 DEFACT,NRAND,TITLE)           
c     >                 DEFACT,NRAND)           
c slmod end
!      use iter_bm
      use variable_wall
      use yreflection
      IMPLICIT none                                                    
      INCLUDE  'params'                                                         
C     INCLUDE  (PARAMS)                                                         
      INTEGER  IMODE,NIZS,NIMPS,IGEOM,NTBS,NTIBS,NNBS,NYMFS
      INTEGER  ITER,NRAND,NCVS,IMPADD,IMPCF
      REAL     QTIM,CPULIM,XWIDM,YWIDM,FSRATE,IONTIM,NEUTIM                     
      REAL     FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)                                
      DOUBLE PRECISION SEED,DEFACT                                              
      CHARACTER PRINPS(-MAXNPS-1:MAXNPS)*7                                      

C     DUMMY VARIABLES FOR DEBUGGING, WHEN NECESSARY 
C     INTEGER IND123,IND124

C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *   LIM3: MAIN CONTROLLING ROUTINE                                  *        
C  *   ------------------------------                                  *        
C  *                                                                   *        
C  *      THIS  ROUTINE   FOLLOWS THE DIFFUSION WITH TIME              *        
C  *   OF A SET OF INJECTED IMPURITY IONS AND RETURNS                  *        
C  *   EITHER THE STATE OF THE ION CLOUD AT A SET OF TIME              *        
C  *   POINTS (IMPULSE MODE), OR THE STEADY STATE DISTRIBUTION         *        
C  *   (STEADY STATE MODE) OR BOTH.                                    *        
C  *      EACH ION POSSESSES (X,Y,P) COORDINATES, A Y DIRECTION        *        
C  *   VELOCITY, A TEMPERATURE AND AN IONISATION STATE, ALL OF WHICH   *        
C  *   CHANGE WITH TIME.                                               *        
C  *      THE IONS ARE FOLLOWED UNTIL THEY ARE ABSORBED, IONISE        *        
C  *   BEYOND A GIVEN STATE OR UNTIL A CUTOFF TIME IS REACHED.         *        
C  *      NOTE THAT THE TAU FACTORS, THE BACKGROUND TEMPERATURE        *        
C  *   AND DENSITY AND THE LIMITER EDGE ARE FOUND FOR A SET            *        
C  *   OF X POSITIONS BEFORE THE ITERATIONS BEGIN.  DURING THE         *        
C  *   ITERATIONS, THE VALUES NEAREST THE CURRENT (X,Y,P) POSITION     *        
C  *   ARE USED.                                                       *        
C  *      SIMILARLY THE OUTBOARD ELECTRIC FIELD AND DRIFT VELOCITY     *        
C  *   ARE CALCULATED FOR A SET OF Y POSITIONS ALONG WITH A            *        
C  *   SCALING FACTOR FOR EACH X POSITION OUTBOARD SO THAT             *        
C  *   DURING THE ITERATION A VALUE CAN BE CALCULATED FOR ANY          *        
C  *   POSITION BY TAKING THE PRODUCT OF THE APPLICABLE X AND Y        *        
C  *   POSITION VALUES.                                                *        
C  *      THIS REQUIRES MAKING THE ASSUMPTION THAT THE ELECTRIC FIELD  *        
C  *   VALUES ARE PROPORTIONAL TO TEMB/L AND THAT THE BACKGROUND       *        
C  *   VELOCITIES ARE PROPORTIONAL TO SQRT(TEMB).                      *        
C  *      TIMES ARE SCALED BY 1/QTIM. THIS MEANS THAT TIME VALUES      *        
C  *   CAN BE STORED IN INTEGERS. INPUT TIME VALUES ARE ROUNDED TO     *        
C  *   THE NEAREST INTEGER (IE THE NEAREST TIMESTEP).                  *        
C  *      Y DIRECTION VELOCITY VALUES ARE SCALED BY QTIM.              *        
C  *   THIS AGAIN SAVES AN INNER LOOP MULTIPLICATION AS THE CHANGE     *        
C  *   IN Y AT EACH TIME STEP BECOMES Y = Y + VY.                      *        
C  *      LIM MONITORS VARIOUS ASPECTS OF THE DIFFUSION AND            *        
C  *   PRINTS OUT A SET OF DIAGNOSTICS AT THE END OF THE RUN.          *        
C  *                                                                   *        
C  *  ARGUMENTS :-                                                     *        
C  *  IMODE  : SET TO 1 FOR IMPULSE MODE, 2 FOR STEADY STATE MODE      *        
C  *               AND 0 FOR BOTH                                      *        
C  *  NIZS   : MAXIMUM IONIZATION STATE TO BE FOLLOWED                 *        
C  *  NIMPS  : NUMBER OF IMPURITY IONS TO BE USED IN MONTE CARLO       *        
C  *  IMPADD : NUMBER OF ADDITIONAL NEUTRALS TO BE LAUNCHED TO         *
C  *           SIMULATE CROSS-FIELD FLUXES                             *
C  *  IMPCF  : NUMBER OF CROSS-FIELD SPUTTERED PRIMARY IMPURITIES      *
C  *           THE NUMBER IS CALCULATED IN NEUT IF CFBGFF.GT.0.0       *
C  *           OTHERWISE IT IS ZERO (INITIALIZED AT START OF RUNLM3)   *
C  *  QTIM   : SIZE OF QUANTUM TIMESTEP TO BE USED FOLLOWING IONS      *        
C  *  CPULIM : MAXIMUM AMOUNT OF TIME TO BE USED BY ROUTINE (SECS)     *        
C  *           A VALUE OF 0 INDICATES INFINITE TIME                    *        
C  *  XWIDM  : MINIMUM X BIN WIDTH, CALCULATED FROM BIN BOUNDARIES     *        
C  *  YWIDM  : MINIMUM Y BIN WIDTH, CALCULATED FROM BIN BOUNDARIES     *        
C  *  FSRATE : SIZE OF TIMESTEP TO BE USED FOLLOWING NEUTRALS          *        
C  *  IONTIM : CPU TIME SPENT FOLLOWING IONS ACCUMULATOR               *        
C  *  NEUTIM : CPU TIME SPENT FOLLOWING NEUTRALS ACCUMULATOR           *        
C  *  SEED   : CURRENT RANDOM NUMBER SEED VALUE FOR SURAND (D.P)       *        
C  *  IGEOM  : RADIAL GEOMETRY FLAG:  0 SLAB,  1 CYLINDER              *        
C  *                                                                   *        
C  *  DEFACT : FACTOR FOR CONVERTING DDLIMS TO ACTUAL PHYSICAL DENSITY *        
C  *                                                                   *        
C  *                        CHRIS FARRELL (HUNTERSKIL)  MARCH 1988     *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE   'dynam1'                                                        
C     INCLUDE   (DYNAM1)                                                        
      INCLUDE   'dynam3'                                                        
C     INCLUDE   (DYNAM3)                                                        
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
      INCLUDE   'commv'                                                         
C     INCLUDE   (COMMV)                                                         
      INCLUDE   'zommv'                                                         
C     INCLUDE   (ZOMMV)                                                         
      INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
      INCLUDE   'cneut'                                                         
C     INCLUDE   (CNEUT)                                                         
      INCLUDE   'cnoco'                                                         
C     INCLUDE   (CNOCO)                                                         
      INCLUDE   'save'                                                          
C     INCLUDE   (SAVE)                                                          
      INCLUDE   'comnet'                                                        
C     INCLUDE   (COMNET)                                                        
      INCLUDE   'crand'                                                         
C     INCLUDE   (CRAND)                                                         
c slmod begin
      INCLUDE   'slcom'
c slmod end
c
      include 'cadas'
C                                                                               
c      REAL      RADDEG,PI,VY0,VY02,PARTIM,P                            
c
      REAL      VY0,VY02,PARTIM,P                            
      REAL      STATIM,VFLUID,TWALLN,TDEP,ZA02AS,RIZB,RFAIL,TFAIL               
      REAL      SPARA,TIMMAX,TNEUT,TATIZ,TWALL,AVPPOS,RDIFFT,EDGE2              
      REAL      FRQTIM,EDGE1,SSEF,YEFF,YFACT,RRES,TRES,RRES1                    
      REAL      AVXPOS,AVYPOS,TTMAX,TCENT,TBYOND,TCUT,RATIZ,RNEUT,QUANT         
      REAL      TAVXPOS
      REAL      FVYCOL,Y,SVY,ABSY,RWALLN,RCENT,RTMAX,RDEP(2),RWALL(2)    
      REAL      ABSP
      REAL      PORM,OLDY,TEMP(-MAXNYS:MAXNYS),YLDTOT(2),YLDMAX(2)              
      REAL      SPUTY,RMACH,ENERGY,RNEUT1,RYIELD,OLDALP,YTHTOT(2)               
c slmod begin
c
c Moved to common block SLCOM:
c
c      REAL      SVHINS(-MAXQXS:MAXQXS),SEYINS(-MAXQXS:MAXQXS,MAXIZS)            
c
c slmod end
      REAL      YIELD,GYTOT1,GTOT1,TBELOW,SPUNEW,RANDEP                         
      REAL      RSTRUK,TSTRUK,TEMOLD,FACT,RAN,EMAX                              
c
c     jdemod
c
c      REAL      MAT1,MAT2
c
      integer   mat1,mat2
c
      real      tptrac(maxlen,2)
      INTEGER   KKLIM,KK,ICUT(2),NATIZ,NPROD,IP,IFATE,STATUS                    
      INTEGER   IPOS,IQX,IQY,IX,IY,IZ,MAXCIZ,IC,II,IOY,IOD,IO                   
      INTEGER   IMP,IMPLIM,MATLIM,J,JY,JX,IT,MPUT,IN
      CHARACTER WHAT(51)*10,FATE(9)*16,STRING*21                                
      REAL      POLODS(-MAXQXS:MAXQXS),TIMUSD,XM,YM                             
      REAL      SVPOLS(-MAXQXS:MAXQXS)
      REAL      RIONS(MAXIZS),STOTS(20)                                         
      REAL      CISTOT,CISMAX,RSTMIN,TSTEPL,RCONST,SVYBIT,AVAPOS                
      REAL      SDTZS(MAXIZS),QFACT,YYCON,YY,ALPHA                              
      REAL      TSPLIT(MAXINS),TRULET(MAXINS),SDYZS(MAXIZS)                     
      REAL      FACTDEPS
      REAL      SVG,SVYMIN,SVYMOD
      REAL      DPPROB
      INTEGER   NSPLIT(MAXINS),NRULET(MAXINS),IS,IPUT,IGET(0:MAXPUT)            
      integer   traclen 
      LOGICAL   DIFFUS,RESPUT,RES,BIGTRAC                                               
c
c     Add some local variables related to calculating the scaling of the NERODS3 data
c
      real pbnd1,pbnd2,local_pwid
c
c     Add iqy_tmp to support variable wall location
c
      integer :: iqy_tmp

c
c     ADD LOGICAL to record if splitting and rouletting is active to avoid
c     a bug if ALPHA > CXSPLS(IS) = 2*CA in one diffusive step  
c
      logical   split  
c
      DOUBLE PRECISION DSPUTY,DTOTS(20),DTEMI,DQFACT,DELTAX                     
      DOUBLE PRECISION DACT,DEMP(-MAXNYS:MAXNYS,4),DY1,DWOL,DSUM4               
      DOUBLE PRECISION DY2,DSUM1,DSUM2,DSUM3,DIZ,DOUTS(MAXIZS),DIST             
c slmod begin
      REAL       IONCNT,IONPNT
      REAL       RAN1,RAN2,RGAUSS,VPARA,TPARA,VPARAT
      REAL       AVGTRAC    
      REAL       TARGET
      CHARACTER  TITLE*80
c slmod endC
C                                                                               
      DATA  FATE  /'REACHED X=AW',        'HIT Y=0 FROM Y>0',                   
     >             'REACHED Y=2L',        'HIT Y=0 FROM Y<0',                   
     >             'REACHED Y=-2L',       'REACHED TIME CUT',                   
     >             'SPLITTING ION',       'ROULETTE DISCARD',                   
     >             'CHARGE CHECK'/                                              
C                                                                               
      DATA  WHAT  /'  PRIMARY ',   ' SECONDARY',   ' TERTIARY ',                
     >             'QUATERNARY',   '  QUINARY ',   '  SIXTH   ',                
     >             '  SEVENTH ',   '  EIGHTH  ',   '  NINTH   ',                
     >             '  TENTH   ',   ' ELEVENTH ',   ' TWELFTH  ',                
     >             'THIRTEENTH',   'FOURTEENTH',   ' FIFTEENTH',                
     >             ' SIXTEENTH',   ' SEVENTEEN',   'EIGHTEENTH',
     >             'NINETEENTH',   ' TWENTIETH',   'TWENTY-1ST',
     >             'TWENTY-2ND',   'TWENTY-3RD',   'TWENTY-4TH',
     >             'TWENTY-5TH',   'TWENTY-6TH',   'TWENTY-7TH',
     >             'TWENTY-8TH',   'TWENTY-9TH',   ' THIRTIETH',
     >             'THIRTY-1ST',   'THIRTY-2ND',   'THIRTY-3RD',
     >             'THIRTY-4TH',   'THIRTY-5TH',   'THIRTY-6TH',
     >             'THIRTY-7TH',   'THIRTY-8TH',   'THIRTY-9TH',
     >             ' FORTIETH ',   ' FORTY-1ST',   ' FORTY-2ND',
     >             ' FORTY-3RD',   ' FORTY-4TH',   ' FORTY-5TH',
     >             ' FORTY-6TH',   ' FORTY-7TH',   ' FORTY-8TH',
     >             ' FORTY-9TH',   ' FIFTIETH ',   
     >             '    ALL   '/                                                
C                                                                               
c      DATA      RADDEG /57.29577952/, PI     /3.141592654/                      
C                                                                               
C-----------------------------------------------------------------------        
C                   INITIALISATION                                              
C-----------------------------------------------------------------------        
C                                                                               
c slmod
      WRITE(0,*) 'Begin LIM3'

      IF (optdp.EQ.1) THEN
        WRITE(0,*) 'Warning! Hard code adjustment to bin location',
     +             ' for DIVIMP ion profile.'
      ENDIF 

      DO II = 1, NBIN
        BSBIN(II) = BSBIN(II) + 0.5
        YSBIN(II) = YSBIN(II) + 0.5
      ENDDO

      AVGTRAC = 0.0
      TGLOSS  = 0.0
      WLOSS   = 0.0
      LLOSS   = 0.0
      IZLOSS  = 0.0
      TSLOSS  = 0.0
      ALOSS   = 0.0
      MARK    = 0.005
      TARGET  =-5.0

      DO IY=-NYS-1,NYS+1
        INJBINT(IY) = 0
      ENDDO
 
c      DO IP=-MAXNPS,MAXNPS
c        INJBINP(IP) = 0
c      ENDDO
c slmod end 

      TAVXPOS = 0.0
      IF (CPULIM.LE.0.0) CPULIM = 1.0E7                                         
      DWOL = DBLE (CTWOL)                                                       
      IF (CDPSTP.NE.0.0) THEN
         DPPROB = 2.0*QTIM / CDPSTP /CDPSTP
      ELSE
         DPPROB = 0.0
      ENDIF
C                                                                               
C---- PLASMA ELONGATION - SET-UP DELPS ARRAY AND CONO, CONI CONSTANTS.          
C                                                                               
      CONI = (CKI-1.0) / CHALFL                                                 
      CONO = (CKO-1.0) / CHALFL                                                 
      DO 20 IY = 1, NYS                                                         
        YY = MOD (YOUTS(IY), CL)                                                
        IF (YY.GT.CHALFL) YY = CL - YY                                          
        DO 10 IX = 1, NXS                                                       
          IF (XS(IX).GT.0.0) THEN                                               
            DELPS(IX,IY) = 1.0 / (YY*CONI + 1.0)                                
          ELSE                                                                  
            DELPS(IX,IY) = 1.0 / (YY*CONO + 1.0)                                
          ENDIF                                                                 
   10   CONTINUE                                                                
   20 CONTINUE                                                                  
C                                                                               
C---- SET UP FACTORS IN COMMON COMTAU                                           
C                                                                               
      IF (ITER.EQ.1) CALL TAUIN1 (QTIM,NIZS,ICUT,FSRATE,IGEOM,                  
     >                            NTBS,NTIBS,NNBS)                           
C                                                                               
C---- SET UP VFLUID THE FLUID VELOCITY.                                         
C                                                                               
      IF(ABS(CVHYS(1)).EQ.0.0) THEN                                             
         VFLUID = 1.56E4 * SQRT (CTBIN/CRMB)                                    
      ELSE                                                                      
         VFLUID = ABS (CVHYS(1))                                                
      ENDIF                                                                     
      WRITE (6,*) 'LIM3: VFLUID=', VFLUID                                       
C                                                                               
C---- SET UP FVYCOL: POST COLLISION VELOCITY FACTOR    *** NOT USED             
C----                            SO NO NEED FOR AN ARRAY OF IQX VALS !          
C----        VY0   : INITIAL VELOCITY FOR NON-NEUT CASES                        
C----        VY02  : SECOND VELOCITY FOR INJECTION 4
C----        POLODS: POLOIDAL DIFFUSION FACTOR                                  
C----        SVPOLS: SCALED POLOIDAL DRIFT VELOCITY 
C----        SVHINS: SCALED INBOARD PLASMA FLOW VELOCITY                        
C----        SEYINS: SCALED INBOARD ELECTRIC FIELD                              
C----                (DELTAT.DELTAT.ZI.E/MP/MI IS SCALE FACTOR)                 
C                                                                               
c
c     Note: 1.56e4 is associated with a velocity calculated from
c           v = sqrt ( 8 k T / (PI m) ) 
c
c           1.38e4 is associated with a velocity calcualted from 
c           v = sqrt ( 2 kT / m ) 
c
c

      FVYCOL = QTIM * 1.56E+04 / SQRT(CRMI)                                     
      IF (CIOPTE.EQ.1.OR.CIOPTE.EQ.3.OR.CIOPTE.EQ.6.OR.CIOPTE.EQ.8
     >     .OR.CIOPTE.EQ.9.or.ciopte.eq.13) THEN        
        VY0 = 1.56E4 * SQRT (CTEMSC/CRMI)
        
      ELSEIF  (CIOPTE.EQ.12) THEN
        VY0 = -1.56E4 * SQRT (CTEMSC/CRMI)
      ELSEIF  (CIOPTE.EQ.4) THEN
C
C     CTEMSC - USUALLY USED FOR ION TEMPERATURES
C     CENGSC - USUALLY USED FOR NEUTRAL ENERGIES
C     FOR INJECTION OPTION 4 THE ENERGY OF THE IONS IS SPECIFIED
C     BY THE VALUE ENTERED IN CENGSC
C     DAVID ELDER , 1990 , FEB 15 
C
        VY0 =   1.38E4* SQRT(CENGSC/CRMI) 
        VY02 =  1.38E4* SQRT(CEIN2/CRMI)
      ELSE                                                          
        VY0 = 0.0                                                               
      ENDIF                                                                     
C                                                                               
      DO 16 IQX = 1-NQXSO, NQXSI                                                
        POLODS(IQX) = SQRT (2.0 * CDPOL * QTIM * QS(IQX))                       
        SVPOLS(IQX) = QTIM * QS(IQX) * CVPOL
        SVHINS(IQX) = QTIM * QS(IQX) * CVHYIN                                   
        DO 15 IZ = 1, NIZS                                                      
          SEYINS(IQX,IZ) = REAL (IZ) * QTIM * QS(IQX) * QTIM * QS(IQX) *        
     >                     (1.602192E-19/1.672614E-27) * CEYIN / CRMI           
   15   CONTINUE                                                                
   16 CONTINUE                                                                  
c slmod begin - now defunkt I think - April 28, 97
      IF (SLOPT.EQ.1) THEN
        DO ALPHA = 0.0, CA , CA / REAL(NQXSI) 
          IQX = MIN (INT(ALPHA*XSCALI)+1, NQXSI)                          
          WRITE(0,*) 'Before:',IQX,ALPHA,CA,SVHINS(IQX)
          IF (ALPHA.LT.CATIN.OR.ALPHA.GT.0.02) SVHINS(IQX) = 0.0
          WRITE(0,*) 'After :',IQX,ALPHA,CA,SVHINS(IQX)
          WRITE(0,*) ' '
        ENDDO
      ENDIF
c slmod end
C                                                                               
C---- CHECK CPLSMA, POINT OUTSIDE WHICH WE HAVE SPECIAL CHARACTERISTICS         
C---- IF NOT REQUIRED, ARTIFICIALLY SET SO THAT IT NEVER ARISES.                
C                                                                               
      IF (CIOPTB.EQ.0.AND.CIOPTC.EQ.0.AND.CIOPTD.EQ.0) CPLSMA = 2.0*CAW         
      JX = IPOS (CPLSMA*0.999, XS, NXS-1)                                       
C                                                                               
C                     SET UP MONITORING VARIABLES                               
C                                                                               
      RSTMIN = CTIMSC / QTIM                                                    
      IF (NIZS.GT.0) CALL MONINI (RSTMIN, CSTMAX, NIZS, CTBIN)                  
C                                                                               
C---- ZERO ARRAYS                                                               
C                                                                               
c slmod begin
      CALL RZERO (TAG2,    MAXIMP)
c slmod end
      CALL DZERO (DDLIMS, MAXNXS*(2*MAXNYS+1)*(MAXIZS+2))                       
      CALL DZERO (DDTS,   MAXNXS*(2*MAXNYS+1)*MAXIZS)                           
      CALL DZERO (DDYS,   MAXNXS*(2*MAXNYS+1)*MAXIZS)                           
      CALL RZERO (DEPS,   MAXNXS*(MAXIZS+1)*3)
      CALL RZERO (NEROXS, MAXNXS*5*3)                                           
      CALL RZERO (NEROYS, MAXOS*6)                                              
      CALL RZERO (NERODS, MAXOS*5)                                              
      CALL RZERO (NERODS3, MAXOS*5*(2*MAXNPS+1))                                              
      CALL RZERO (WALLS,  (2*MAXNYS+1)*(MAXIZS+4))                              
      CALL RZERO (TIZS,   MAXNXS*(2*MAXNYS+1)*(MAXIZS+2))                       
      CALL RZERO (ZEFFS,  MAXNXS*(2*MAXNYS+1)*6)                                
      CALL DZERO (DDLIM3, MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
      CALL RZERO (TIZ3,   MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
      CALL RZERO (SDTXS,  MAXNXS*MAXIZS)                                        
      CALL RZERO (SDTYS,  (2*MAXNYS+1)*MAXIZS)                                  
      CALL RZERO (SDTZS,  MAXIZS)                                               
      CALL RZERO (SDYXS,  MAXNXS*MAXIZS)                                        
      CALL RZERO (SDYYS,  (2*MAXNYS+1)*MAXIZS)                                  
      CALL RZERO (SDYZS,  MAXIZS)                                               
      CALL DZERO (DOUTS,  MAXIZS)                                               
      CALL RZERO (RIONS,  MAXIZS)                                               
      IF (IMODE.NE.2)                                                           
     >CALL RZERO (LIM5,   MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1)*          
     >                    MAXNTS)                                               
C
      CALL RZERO (SVYBAR,2*MAXQXS+1)
      CALL RZERO (SVYACC,2*MAXQXS+1)          
C                                                                               
C-----------------------------------------------------------------------        
C     PRINT SELECTED PARAMETERS SET UP ON FIRST ITERATION ONLY                  
C-----------------------------------------------------------------------        
C                                                                               
      IF (ITER.EQ.1) THEN                                                       
        IF     (IMODE.EQ.0) THEN                                                
         CALL PRC ('OPERATION MODE 0  COMBINED IMPULSE & STEADY STATE')         
        ELSEIF (IMODE.EQ.1) THEN                                                
         CALL PRC ('OPERATION MODE 1  IMPULSE')                                 
        ELSEIF (IMODE.EQ.2) THEN                                                
         CALL PRC ('OPERATION MODE 2  STEADY STATE')                            
        ENDIF                                                                   
C                                                                               
        CALL PRI ('  NO OF X BINS                     ', NXS)                   
        CALL PRI ('  NO OF Y BINS                     ', 2*NYS)                 
        CALL PRI ('  NO OF Y BINS FOR 3D RESULTS      ', 2*NY3D)                
        CALL PRI ('  NO OF P BINS                     ', 2*MAXNPS+1)            
        CALL PRI ('  NO OF X PTS FOR OUTBOARD FACTORS ', NQXSO)                 
        CALL PRI ('  NO OF X PTS FOR INBOARD FACTORS  ', NQXSI)                 
        CALL PRI ('  NO OF Y POINTS FOR FACTORS       ', NQYS)                  
        CALL PRI ('  NO OF T POINTS FOR TIME RESULTS  ', NTS)                   
        CALL PRI ('  MAXIMUM IONIZATION STATE         ', NIZS)                  
        CALL PRI ('  NO OF IMPURITY IONS TO FOLLOW    ', NIMPS)                 
        IF (IMPADD.GT.0) THEN 
           CALL PRI ('  NO OF EXTRA IMPURITY NEUTALS     ',IMPADD)
           CALL PRR ('     LAUNCHED ON Y =           +/- ',YCFADD)
           IF (CEXNEUT.EQ.0) THEN
              CALL PRC ('     WITH A UNIFORM DISTRIBUTION') 
           ELSEIF (CEXNEUT.EQ.1) THEN 
              CALL PRC ('     WITH A NORMAL DISTRIBUTION') 
           ENDIF
        ENDIF  
        CALL PRR ('  SIZE OF NEUT QUANTUM TIMESTEP (S)', FSRATE)                
        CALL PRR ('  BASE VALUE FOR LIM TIMESTEP   (S)', QTIM)                  
        WRITE (7,'(1X,''  RANDOM NUMBER SEED'',10X,I15)') CISEED                
C                                                                               
        IF     (IGEOM.EQ.0) THEN                                                
          CALL PRC ('  RADIAL GEOMETRY OPTION              SLAB')               
        ELSEIF (IGEOM.EQ.1) THEN                                                
          CALL PRC ('  RADIAL GEOMETRY OPTION              CYLINDER')           
        ENDIF                                                                   
C                                                                               
        IF (CFTCUT.GT.CAW)                                                      
     >    CALL PRR ('  STOP CROSS FIELD TRANSPORT AT X =', CFTCUT)              
C                                                                               
        IF     (CPRINT.EQ.0) THEN                                               
          CALL PRC ('  PRINT OPTION                        REDUCED')            
        ELSEIF (CPRINT.EQ.1) THEN                                               
          CALL PRC ('  PRINT OPTION                        FULL')               
        ENDIF                                                                   
C                                                                               
        IF (ABS(CTHETB-90.0).LT.1.0E-3) THEN                                    
          CALL PRR ('  LIMITER GEOMETRY POLOIDAL: THETAB', CTHETB)              
        ELSE                                                                    
          CALL PRR ('  LIMITER GEOMETRY TOROIDAL: THETAB', CTHETB)              
        ENDIF                                                                   
C                                                                               
        IF (CXSPLS(1).LT.CA) then                                                   
          CALL PRC ('  SPLITTING & ROULETTING IN OPERATION')                    
          split = .true. 
        else
          split = .false.
        endif  
C                                                                               
        CALL PRR ('  ELONGATION PARAMETER OUTBOARD  KO', CKO)                   
        CALL PRR ('  ELONGATION PARAMETER INBOARD   KI', CKI)                   
C                                                                               
        IF (CANAL.LT.CA)                                                        
     >    CALL PRR ('  ANALYTIC EXTENSION INBOARD OF X =', CANAL)               
C                                                                               
        IF (CPRINT.EQ.1) THEN                                                   
          CALL PRR ('  SMALLEST X BIN SIZE  (M)         ', XWIDM)               
          CALL PRR ('  DELTA X INBOARD FOR FACTORS  (M) ', 1.0/XSCALI)          
          CALL PRR ('  DELTA X OUTBOARD FOR FACTORS  (M)', 1.0/XSCALO)          
          CALL PRR ('  SMALLEST Y BIN SIZE  (M)         ', YWIDM)               
          CALL PRR ('  DELTA Y ALONG X=0 FOR FACTORS (M)', 1.0/YSCALE)          
          CALL PRR ('  STOP TIME FOR ITERATION  (S)     ', TIMMAX)              
          CALL PRR ('  ALLOCATED CPU TIME  (S)          ', CPULIM)              
          CALL PRR ('  "NEAR LIMITER" MEANS  0.0 <= X <=', CXNEAR)              
          CALL PRR ('                AND YP- BETWEEN +/-',CYNEAR/CLFACT)        
          CALL PRB                                                              
          CALL PRC ('BOUNDARIES OF X BINS')                                     
          WRITE (7,'(8(1X,F7.4))') CAW,(XS(IX),IX=1,NXS)                        
          CALL PRB                                                              
          CALL PRC ('BOUNDARIES OF Y BINS (MIRRORED FOR -Y REGION)')            
          WRITE (7,'(8F8.4)') 0.0,(YS(IY),IY=1,NYS)                             
          CALL PRB                                                              
          CALL PRC ('BOUNDARIES OF P BINS')                                     
          WRITE (7,'(8(1X,A7))') (PRINPS(IP),IP=-MAXNPS-1,MAXNPS)               
        ENDIF                                                                   
C                                                                               
        IF (IMODE.NE.2) THEN                                                    
          CALL PRB                                                              
          CALL PRC ('DWELL TIME FACTORS AND SPECIFIC OUTPUT TIMES (S)')         
          WRITE (7,9020) (DWELFS(IT),IT=1,NTS)                                  
          DO 30 IZ = 0, NIZS                                                    
            IF (IZ.EQ.0 .AND. CNEUTA.NE.0) GOTO 30                              
            WRITE (7,9021) IZ,(DWELTS(IZ)*DWELFS(IT),IT=1,NTS)                  
            WRITE (6,9022) IZ,(CTIMES(IT,IZ),IT=1,NTS)                          
   30     CONTINUE                                                              
        ENDIF                                                                   
C                                                                               
      RIZB   = REAL (CIZB)                                                      
C
C     PLACE CALL TO INITIALIZE THE SPUTTERING YIELD DATA SO THAT 
C     IT IS AVAILABLE FOR ION AS WELL AS NEUT LAUNCH CASES.
C
C
C     LOAD YIELD COMMON BLOCK WITH APPROPRIATE DATA
C
      IF (CSPUTOPT.EQ.1) THEN
c
c        CALL SYIELD(MATLIM,MAT1,MAT2,CNEUTD,CBOMBF,CBOMBZ,CION,
c     >            CIZB,CRMB,CTSUB)   
c
        CALL SYIELD (MATLIM,MAT1,CNEUTD,
     >               CBOMBF,CBOMBZ,CION,CIZB,CRMB,CEBD)
      ELSE IF (CSPUTOPT.EQ.2) THEN
        CALL SYLD93 (MATLIM,MAT1,CNEUTD,
     >               CBOMBF,CBOMBZ,CION,CIZB,CRMB,CEBD)
      ELSE IF (CSPUTOPT.EQ.3.or.csputopt.eq.4.or.csputopt.eq.5)THEN
        CALL SYLD96 (MATLIM,MAT1,CNEUTD,
     >               CBOMBF,CBOMBZ,CION,CIZB,CRMB,CEBD)
      ENDIF
c

c
c     Set up MAT2 if cneutd.eq.2 - this is code from the original 
c     LIM version of SYIELD
c
      call syield_set_mat2(mat2,cneutd,cbombf,cbombz)
c
      write (0,
     >    '((1x,a,1x,i6),(1x,a,1x,f10.2),3(1x,a,1x,i6),1x,a,1x,g12.5)') 
     >            'Materials: Data Opt=',csputopt,
     >            'Incident Angle=',extra_sputter_angle,
     >            'Plasma Mat1=',mat1,
     >            'Plasma Mat2=',mat2,
     >            'Limiter Mat=',matlim,
     >            'Binding En =',cebd
c
      call test_phys_yld(matlim,mat1)

c
        CQPL =  122. * EXP (-9048./CTSUB)                                       
        CQSL = 1014. * EXP (-9048./CTSUB)                                       
        CALL PRDATA (NIZS,XSCALO,XSCALI)                                       
C                                                                               
C------ CONVERT CSNORM FROM DEGREES INTO RADIANS  (PRDATA PRINTS IT OUT)        
C                                                                               
        CSNORM = CSNORM / RADDEG                                                
        IF (CPRINT.EQ.1) THEN                                                   
          CALL PRB                                                              
          CALL PRC ('SIMPLE FACTORS     ')                                      
          CALL PRR ('  FACTOR FOR POST COLLISION VELOCITY  ', FVYCOL)           
          CALL PRR ('  FACTOR FOR POLOIDAL DIFFUSION       ',                   
     >                      POLODS(0)/SQRT(QS(0)))                              
          CALL PRR ('  FACTOR FOR INBOARD PLASMA FLOW VEL. ',                   
     >                      SVHINS(0)/QS(0))                                    
          CALL PRR ('  FACTOR FOR INBOARD E, IZ STATE 1    ',                   
     >                      SEYINS(0,1)/(QS(0)*QS(0)))                          
        ENDIF                                                                   
C                                                                               
        CALL TAUPR1 (QTIM,NIZS)                                                 
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C    LAUNCH PRIMARY NEUTRALS EN MASSE                                           
C-----------------------------------------------------------------------        
C                                                                               
C---- TNEUT : TOTAL NUMBER OF NEUTRALS LAUNCHED                                 
C---- TATIZ : TOTAL NO OF IONS CREATED                                          
C---- TWALL : TOTAL NO OF IONS REACHING WALL                                    
C---- TWALLN: TOTAL NO OF NEUTRALS REACHING WALL                                
C---- TDEP  : TOTAL NO OF IONS DEPOSITED ON LIMITERS                            
C---- TTMAX : TOTAL NO OF NEUTRALS EXISTING AT TMAX                             
C---- TCENT : TOTAL NO OF NEUTRALS REACHING CENTRE                              
C---- TBYOND: TOTAL NO OF IONS IONISED BEYOND LIMIT                             
C---- TBELOW: TOTAL NO OF IONS RECOMBINING TO NEUTRALS                          
C---- TCUT  : TOTAL NO OF IONS EXISTING AT TMAX                                 
C---- TSTRUK: TOTAL NO OF NEUTRALS STRIKING LIMITER                             
C---- TFAIL : TOTAL NO OF FAILED NEUTRAL LAUNCHES                               
C                                                                               
      TNEUT  = 0.0                                                              
      TRES   = 0.0                                                              
      TATIZ  = 0.0                                                              
      TWALL  = 0.0                                                              
      TWALLN = 0.0                                                              
      TDEP   = 0.0                                                              
      TTMAX  = 0.0                                                              
      TCENT  = 0.0                                                              
      TBYOND = 0.0                                                              
      TBELOW = 0.0                                                              
      TCUT   = 0.0                                                              
      TSTRUK = 0.0                                                              
      TFAIL  = 0.0                                                              
C                                                                               
C---- TSPLIT, TRULET, NSPLIT, NRULET:  FOR TABLE OF SPLITTING DATA              
C---- MPUT: GREATEST SIZE OF SPLITTING BANK.                                    
C                                                                               
      IPUT = 0                                                                  
      MPUT = 0                                                                  
      DO 50 IS = 1, MAXINS                                                      
        TSPLIT(IS) = 0.0                                                        
        TRULET(IS) = 0.0                                                        
        NSPLIT(IS) = 0                                                          
        NRULET(IS) = 0                                                          
   50 CONTINUE                                                                  
C                                                                               
C---- RNEUT1: NUMBER OF PRIMARY NEUTRALS LAUNCHED                               
C---- RNEUT : NUMBER OF NEUTRALS LAUNCHED IN CURRENT BATCH                      
C---- RATIZ : NUMBER OF IONS CREATED IN CURRENT BATCH                           
C---- RSTRUK: NUMBER OF NEUTRALS STRIKING LIMITER                               
C---- STATUS: 1,2,..CMAXGENS FOR PRIMARY,2ND,3RD LAUNCHES ETC, 31 FOR TOTAL. 
C                                                                               
C
C                                                                               
C---- FIT INTERPOLATING CURVE TO SET OF YIELD MODIFIER VALUES                   
C---- CALCULATE INTERPOLATED/EXTRAPOLATED YMF AT EACH OUTBOARD X POSN.          
C---- DIFFERENT YIELD MODIFIERS FOR EACH SIDE OF Y = 0                          
C---- FLAG DETERMINES WHETHER TO APPLY TO PRIMARIES, SECONDARIES, BOTH          
C                                                                               
      DO 55 J = 1, 2                                                            
        DO 55 IQX = 1-NQXSO, 0                                                  
          CYMFPS(IQX,J) = 1.0                                                   
          CYMFSS(IQX,J) = 1.0                                                   
   55 CONTINUE                                                                  
C                                                                               
      IF (CYMFLG.NE.-2) THEN                                                    
        CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,2),                               
     >               NQXSO,QXS(1-NQXSO),CYMFPS(1-NQXSO,1),'LINEAR')             
        CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,3),                               
     >               NQXSO,QXS(1-NQXSO),CYMFPS(1-NQXSO,2),'LINEAR')             
      ENDIF                                                                     
C                                                                               
      IF (CYMFLG.NE.-1) THEN                                                    
        CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,2),                               
     >               NQXSO,QXS(1-NQXSO),CYMFSS(1-NQXSO,1),'LINEAR')             
        CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,3),                               
     >               NQXSO,QXS(1-NQXSO),CYMFSS(1-NQXSO,2),'LINEAR')             
      ENDIF                                                                     
C
C
C

      IF (CNEUTA.EQ.0) THEN                                                     
        STATUS = 1                                                              
        WRITE (0,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'          
        WRITE (6,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'          
        WRITE (7,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'          
c
        CALL NEUT (NATIZ,FSRATE,RRES,                                           
     >             ICUT,MATLIM,MAT1,MAT2,NIMPS,IMPADD,IMPCF,
     >             QTIM,GTOT1,GYTOT1,RSTRUK,     
     >             RATIZ,RNEUT,RWALLN,RCENT,RTMAX,SEED,NRAND,                   
     >             NEUTIM,RFAIL,NYMFS,NCVS,STATUS)         


        IF (NATIZ.EQ.0) GOTO 806                                                
      ELSE                                                                      
        STATUS = 1                                                             
        NATIZ  = NIMPS                                                          
        RATIZ  = REAL (NIMPS)                                                   
        RNEUT  = 0.0                                                            
        RWALLN = 0.0                                                            
        RCENT  = 0.0                                                            
        RTMAX  = 0.0                                                            
        RSTRUK = 0.0                                                            
        RFAIL  = 0.0                                                            
        RRES   = 0.0                                                            
      ENDIF                                                                     

      RNEUT1 = RNEUT                                                            
      RRES1  = RRES                                                             
C
C       SET UP MINIMUM VELOCITY BY MULTIPLYING BY QTIM. THIS 
C       MINIMUM IS CONSTANT OVER SPACE (I.E. NOT MULTIPLIED BY 
C       THE QS(IQX) FACTOR)
C
        SVYMIN = CSVYMIN * QTIM  
C                                                                               
C-----------------------------------------------------------------------        
C    FOLLOW IONS TO ABSORPTION / EVENTUAL FATE ...                              
C    THIS CONTINUATION POINT IS TAKEN AFTER SECONDARY NEUTRALS ARE              
C    LAUNCHED FROM THOSE IONS THAT STRUCK THE LIMITERS.                         
C-----------------------------------------------------------------------        
C                                                                               
  200 CONTINUE                                                                  
c
c  jdemod - moved inside the main loop so that self sputtering is fractioned
c           off properly as well. 
c
c slmod 
      IONCNT  = 0
      if (cneuta.eq.0) then 
         IONPNT  = 0.1 * REAL(NATIZ)
      else
         IONPNT  = 0.1 * REAL(NIMPS)
      endif
c slmod end

C                                                                               
C---- SET UP TAU PARALLEL,HEATING,STOPPING: MAY DEPEND ON CTEMSC                
C---- VALUE CALCULATED IN LAUNCH/NEUT.  PRINT SAMPLE VALUES.                    
C                                                                               
      IF (NIZS.GT.0 .AND. ITER.EQ.1) CALL TAUIN2 (QTIM,NIZS)                    
      IF (NIZS.GT.0) CALL TAUPR2 (QTIM,NIZS)                                    
C                                                                               
C---- RDIFFT: TIME TO FIRST DIFFUSION                                           
C                                                                               
      NPROD  = 0                                                                
      AVXPOS = 0.0                                                              
      AVAPOS = 0.0                                                              
      AVYPOS = 0.0                                                              
      AVPPOS = 0.0                                                              
      DO 201 J = 1, 2                                                           
        RWALL (J) = 0.0                                                         
        RDEP  (J) = 0.0                                                         
        YLDTOT(J) = 0.0                                                         
        YTHTOT(J) = 0.0                                                         
        YLDMAX(J) = 0.0                                                         
  201 CONTINUE                                                                  
      RANDEP = 0.0                                                              
      RDIFFT = 0.0                                                              
      IF (STATUS.LE.50) THEN 
         WRITE (0,9012) '***  FOLLOWING ',WHAT(STATUS),'   IONS    ***'      
         WRITE (6,9012) '***  FOLLOWING ',WHAT(STATUS),'   IONS    ***'      
         WRITE (7,9012) '***  FOLLOWING ',WHAT(STATUS),'   IONS    ***'      
      ELSE 
         WRITE (0,9013) '***  FOLLOWING ',STATUS,' GENERATION IONS ***'  
         WRITE (6,9013) '***  FOLLOWING ',STATUS,' GENERATION IONS ***'  
         WRITE (7,9013) '***  FOLLOWING ',STATUS,' GENERATION IONS ***'  
      ENDIF
C                                                                               
C---- SWITCH ON DEBUG IF REQUIRED ...                                           
C                                                                               
      DEBUGL = .FALSE.                                                          
      IF (CSTEPL.GT.0.0) THEN                                                   
        WRITE (6,9004) NINT(CSTEPL),QTIM                                        
        DEBUGL = .TRUE.                                                         
      ENDIF                                                                     
      debugt = .false.       
      if (cstept.gt.0) then
        WRITE (6,9006) CSTEPT                                       
        DEBUGT = .TRUE.                                                         
      endif    
      write (6,*) 'cstept:' , cstept,debugt,chalfl
      CISTOT = 0.0                                                              
      CISMAX = 0.0                                                              
c slmod begin
      IF (DEBUGL) WRITE(78,'(A15)')
     +  'SEYINS'
      IF (DEBUGL) WRITE(78,'(G15.4)')
     +   SEYINS(1,1)

      WRITE(78,*) ' '

      IF (DEBUGL) WRITE(78,'(A4,A7,13A12)')
     +   'IMP','STEP','Y','SVY','SVYMOD',
     +   'RGAUSS','VPARA','VPARAT','DY1','DY2','QS','DTEMI','QUANT',
     +   'CFSS','YFACT'
c slmod end
C                                                                               
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
C                                                                               
C                     FOR EACH ION DO                                           
C                                                                               
      IMPLIM = 4                                                                
      STATIM = ZA02AS (1)                                                       
      PORM   = -1.0                                                             
      KK     = 1000 * ISECT                                                     
      KKLIM  = KK - 10                                                          
      DO 800  IMP = 1, NATIZ                                                    
        IF (100*(IMP/100).EQ.IMP)                                               
     >    WRITE (6,'('' LIM3: ION'',I6,'' FINISHED'')') IMP                     
        TSTEPL = CSTEPL                                                         
        PORM   = -1.0 * PORM                                                    
        OLDY   = 0.0                                                            
        OLDALP = 0.0                                                            
c
c       Initialize tracking variables - if debugt
c
        if (debugt) then
c
c          Turn track debugging off if it has done the required number 
c          of particles. 
c
           if (imp.gt.cstept) then 
              debugt = .false.
           else
              bigtrac = .false.
              traclen = 1
c slmod begin
              AVGTRAC = 0.0
c slmod end
           endif
           write(6,*) 'imp:', imp,bigtrac,traclen,cstept,debugt 
        endif 
C                                                                               
C------ SET INITIAL CHARACTERISTICS OF ION                                      
C          CX     X POSITION                                                    
C          Y      Y POSITION                                                    
C          P      P POSITION                                                    
C          SVYBIT SCALED VELOCITY, TO BE MULTIPLIED BY QTIM*QS(IQX)             
C                 TO FORM SVY A FEW LINES LATER WHEN IQX IS KNOWN ...           
C          CIZ    IONIZATION LEVEL                                              
C          CIST   SCALED TIME - NOW A REAL                                      
C          ITIME  NUMBER OF NEXT MEASUREMENT TIME                               
C          IX     X BIN ION IN                                                  
C          IY     Y BIN ION IN                                                  
C          IP     P BIN ION IN                                                  
C          IQX    ARRAY OFFSET FOR READING FACTORS                              
C          ALPHA  X POSITION TRANSLATED TO X AXIS USING ELONGATION DATA         
C                                                                               
C  USE DATA FROM NEUT TO INJECT PARTICLES (+/- ALREADY ACCOUNTED FOR)           
C                                                                               
        CTEMI = CTEMSC                                                          
        CIZ   = CIZSC                                                           
        MAXCIZ= CIZSC                                                           
        IF (CNEUTA.EQ.0) THEN                                                   
          CX    = XATIZS(IMP)                                                   
c
c         This option was added to simulate particle production which does 
c         not recycle locally - the entire target production is introduced
c         into the SOL at a specified y coordinate
c
c         Init_y_coord is an optional input with a default value of 0.0
c
          if (init_y_coord.ne.0.0) then 
             Y = sign(init_y_coord ,yatizs(imp))
          else
             Y     = YATIZS(IMP)                                                   
          endif
c
          P     = PATIZS(IMP)                                                   
          SVYBIT= VINS(IMP)                                                     
          SPUTY = SPUTYS(IMP)                                                   
c slmod begin
          IONCNT     = IONCNT + 1.0      
          TLOSS(IMP) = 0.0
            
          IF (IONCNT.GT.IONPNT) THEN
            WRITE(0,'(A,I5,A)') ' ',INT(IONPNT/REAL(NATIZ)*100.0),' %'
            IONPNT = IONPNT + 0.1*REAL(NATIZ)
          ENDIF
c slmod end
C      
C  USE INJECTION COORDINATES SPECIFIED IN DATAFILE                              
C  LAUNCH ALTERNATELY ON +Y,-Y REGIONS.  SET INITIAL IONISATION STATE.          
C                                                                               
        ELSEIF (CNEUTA.EQ.1) THEN                                               
          CX  = CXSC                                                            
C
C         SVYBIT IS SPECIFIED FIRST : FOR INJECTION 4 IF A RANDOM NUMBER
C         IS GREATER THAN THE PROBABILITY THAT THE VELOCITY IS VY0
C         THEN SVYBIT IS SET TO VY02. DIRECTION ALTERNATES +/-. 
C
          if (ciopte.eq.13) then 
             CALL SURAND (SEED,1,RAN)
             SVYBIT= VY0 * sign(1.0,ran-0.5)
          else
             SVYBIT= VY0 * PORM                                                    
          endif   
C
          IF (CIOPTE.LE.1.OR.CIOPTE.EQ.5.OR.CIOPTE.EQ.6.OR.
     >        CIOPTE.EQ.12.or.ciopte.eq.13) THEN                   
            Y = CYSC * PORM                                                     
          ELSEIF (CIOPTE.EQ.4) THEN
            CALL SURAND (SEED,1,RAN)
            Y = CYSC*2.0*(RAN-0.5)
            CALL SURAND (SEED,1,RAN) 
            IF (RAN.GT.CPROB)  SVYBIT = VY02 * PORM 
            NRAND = NRAND +2    
          ELSEIF (CIOPTE.EQ.9) THEN 
C
C           OVERWRITE THE ASSIGNED X INJECTION VALUE OF CXSC,
C           ASSIGN THE Y-INJECTION VALUE
C
            Y = CYSC * PORM
            CALL SURAND(SEED,1,RAN)
            CX = (X0L-X0S)*RAN + X0S 
            NRAND = NRAND +1
c slmod begin - Gaussian ion injection option
          ELSEIF (CIOPTE.EQ.10) THEN

250         CALL SURAND(SEED,1,RAN)
            NRAND = NRAND +1
            Y = 6.0E-2 * RAN * PORM
          
            CALL SURAND(SEED,1,RAN)
            NRAND = NRAND + 1
            IF (RAN.GT.EXP(-((Y/0.02)**2))) GOTO 250

            IONCNT = IONCNT + 1.0      
            
            IF (IONCNT.GT.IONPNT) THEN
              WRITE(0,'(A,I5,A)') ' ',INT(IONPNT/REAL(NIMPS)*100.0),' %'
              IONPNT = IONPNT + 0.1*REAL(NIMPS)
            ENDIF

          ELSEIF (CIOPTE.EQ.11) THEN
c
c 3D Gaussian ion launch
c
            WRITE(0,*) 'Error! CIOPTE = 11 not implimented.'
            STOP
c slmod end
          ELSE                                                                  
            CALL SURAND (SEED, 1, RAN)                                          
            Y = CTWOL * RAN * PORM                                              
            NRAND = NRAND + 1                                                   
            YY = MOD (ABS(Y), CL)                                               
            IF (YY.GT.CHALFL) YY = CL - YY                                      
            IF (CX.GE.0.0) THEN                                                 
              CX = CX * (YY*CONI + 1.0)                                         
            ELSE                                                                
              CX = CX * (YY*CONO + 1.0)                                         
            ENDIF                                                               
          ENDIF                                                                 
          P     = CPSC                                                          
          SPUTY = 1.0                                                           
C                                                                               
C  SIMULATE NEUT USING A RECTANGULAR INJECTION REGION                           
C                                                                               
        ELSEIF (CNEUTA.EQ.2) THEN                                               
          CALL SURAND (SEED, 1, RAN)                                            
          CX  = (X0S + RAN * (X0L-X0S))                                         
          CALL SURAND (SEED, 1, RAN)                                            
          Y   = (Y0S + RAN * (Y0L-Y0S)) * PORM                                  
          P   = CPSC                                                            
          NRAND = NRAND + 2                                                     
          SVYBIT= VY0 * PORM                                                    
          SPUTY = 1.0                                                           
        ENDIF                                                                   
C                                                                               
        ABSY = ABS (Y)                                                          
        YY = MOD (ABSY, CL)                                                     
        IF (YY.GT.CHALFL) YY = CL - YY                                          
        IF (CX.GE.0.0) THEN                                                     
          ALPHA = CX / (YY*CONI + 1.0)                                          
        ELSE                                                                    
          ALPHA = CX / (YY*CONO + 1.0)                                          
        ENDIF                                                                   
        DSPUTY = DBLE (SPUTY)                                                   
C                                                                               
C  SPREADING: SPREAD OUT IONS BY A PARTIAL ITERATION                            
C             OF FRQTIM*QTIM (NOTE MAY BE NEGATIVE)                             
C                                                                               
        IF (CNEUTA.EQ.1 .AND.(CIOPTE.EQ.5.OR.
     >      CIOPTE.EQ.6.OR.CIOPTE.EQ.7.OR.CIOPTE.EQ.8)) THEN

          FRQTIM = (IMP - 0.5) / RATIZ - 0.5                                    
          IF (ALPHA.GE.0.0) THEN                                                
            IQX = MIN (INT(ALPHA*XSCALI)+1, NQXSI)                              
          ELSE                                                                  
            IQX = MAX (INT(ALPHA*XSCALO), 1-NQXSO)                              
          ENDIF                                                                 
          CALL SURAND (SEED, 1, RAN)                                            
C
C         DECIDE WHICH Y-REGION THE PARTICLE IS IN AND THEN USE  
C         THE INDEX TO ACCESS THE X-DIFF DATA FOR THE 
C         APPROPRIATE REGION
C          
C         D.ELDER NOV 23 1990
C
          IF (Y.LE.0.0) THEN                                                    
             IF (Y.GT.-CHALFL) THEN 
                J = 1
             ELSEIF (Y.LT.-C3HALFL) THEN 
                J = 2
             ELSE 
                J = 3
             ENDIF
          ELSE      
             IF (Y.GT.C3HALFL) THEN 
                J = 1
             ELSEIF (Y.LT.CHALFL) THEN 
                J = 2
             ELSE 
                J = 3
             ENDIF
          ENDIF
          CX = CX + FRQTIM *                                                  
     >           (SIGN (CXBFS(IQX,J),CXCFS(IQX,J)-RAN) + CXAFS(IQX,J))          
          NRAND = NRAND + 1                                                     
          Y      = Y + FRQTIM * SVYBIT * QTIM * QS(IQX)                         
          ABSY   = ABS (Y)                                                      
          YY = MOD (ABSY, CL)                                                   
          IF (YY.GT.CHALFL) YY = CL - YY                                        
          IF (CX.GE.0.0) THEN                                                   
            ALPHA = CX / (YY*CONI + 1.0)                                        
          ELSE                                                                  
            ALPHA = CX / (YY*CONO + 1.0)                                        
          ENDIF                                                                 
        ENDIF                                                                   
c
        IX    = IPOS (ALPHA, XS, NXS-1)                                         
        IY    = IPOS (ABSY,  YS, NYS-1)                                         

        IF (Y.LT.0.0) IY = -IY                                                  
        JY    = IABS (IY)                                                      
        IP    = IPOS (P,    PS, 2*MAXNPS) - MAXNPS - 1                          
        IT    = IPOS (RSTMIN, CTIMES(1,CIZ), NTS)                               
        IS    = 1                                                               
c slmod begin - DIVIMP ion profile
        IF (optdp.EQ.1) THEN
          INJBINT(IY) = INJBINT(IY) + 1
c          INJBINP(IP) = INJBINP(IP) + 1
c          CIZ = CIZ + 1
c          WRITE(0 ,*) 'Injection ',Y,IY,EXP(-(Y/0.02)**2),PORM
          WRITE(60,*) 'Injection ',Y,IY,EXP(-(Y/0.02)**2),PORM
          WRITE(79,*) CX,Y,P,IX,IY,IP
        ENDIF
c slmod end
C                                                                               
C------ SET INITIAL IQX,IQY VALUES.                                             
C------ CAN NOW FINISH OFF SVY BY MULTIPLYING BY TIMESTEP FACTOR ...            
C                                                                               
        IF (ALPHA.GE.0.0) THEN                                                  
          IQX = MIN (INT(ALPHA*XSCALI)+1, NQXSI)                                
        ELSE                                                                    
          IQX = MAX (INT(ALPHA*XSCALO), 1-NQXSO)                                
        ENDIF                                                                   
        IQY = 0                                                                 
        SVY = SVYBIT * QTIM                                           
c
C                                                                               
C------ RECORD INJECTION POSITION.                                              
C                                                                               
        AVXPOS = AVXPOS + CX * SPUTY                                            
        AVAPOS = AVAPOS + ALPHA * SPUTY                                         
        AVYPOS = AVYPOS + ABSY * SPUTY                                          
        AVPPOS = AVPPOS + ABS(P) * SPUTY                                        
C                                                                               
C  CHECK IF ION HAS STARTED ABOVE MAX IONIZATION STATE                          
C                                                                               
        IF ((CIZ .GT. CION) .OR. (CIZ .GT. NIZS))  THEN                         
          TBYOND = TBYOND + SPUTY                                               
          IFATE = 9                                                             
          GOTO 790                                                              
        ENDIF                                                                   
C                                                                               
C------ IF SET TI=TB FOR STATE CIZ APPLIES, BETTER DO IT                        
C                                                                               
        IF (CIZ.EQ.CIZSET) CTEMI = MAX (CTEMI,CTEMBS(IX,IY))                    
C                                                                               
C------ CALCULATE POINT AT WHICH DIFFUSION WILL BE FIRST APPLIED.               
C------ DEPENDS ON DIFFUSION OPTION AS FOLLOWS :-                               
C------ 0) IMMEDIATE DIFFUSION                                                  
C------ 1) AFTER A TIME BASED RANDOMLY ON INITIAL TEMPERATURE,                  
C------    -TAUPARA.LOG$/2,  WHERE $ IN (0,1)                                   
C------ 2) AFTER TIME TAUPARA, TAKING INTO ACCOUNT CHANGES IN TAUPARA           
C------    AS ION HEATS UP                                                      
C                                                                               
        RCONST = 1.E20                                                          
        IF (CDIFOP.EQ.0) THEN                                                   
          IF (CIOPTB.NE.1) RCONST = 0.0                                         
        ELSEIF (CDIFOP.EQ.1) THEN                                               
          IF (CFPS(IX,IY,CIZ).GT.0.0) THEN                                      
            NRAND = NRAND + 1                                                   
            CALL SURAND (SEED, 1, RAN)                                          
            RCONST = -CTEMI * LOG (RAN) / CFPS(IX,IY,CIZ) * QS(IQX)             
          ENDIF                                                                 
        ENDIF                                                                   
        SPARA  = 0.0                                                            
        DIFFUS = .FALSE.                                                        
C                                                                               
C------ SET INITIAL (X,Y) COORDINATES & TI IN DOUBLE PRECISION.                 
C------ DY1 ACCUMULATES NON-DIFFUSIVE CHANGES, DY2 DIFFUSION CHANGES.           
C                                                                               
        DY1   = 0.0D0                                                           
        DY2   = DBLE (Y)                                                        
        QUANT = 0.0                                                             
        DTEMI = DBLE (CTEMI)                                                    
C                                                                               
        CALL MONUP (8,SPUTY)                                                    
C                                                                               
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
C                                                                               
C         ITERATE UP TO NEXT EVENT IN VARIABLE STEPS OF QS(IQX)*QTIM            
C         DEPENDENT ON CURRENT X POSITION.  BASED ON MULTIPLES OF               
C         "STANDARD" ITERATION TIME QTIM ITSELF.                                
C         DIST RECORDS THE ITERATION NUMBER, CIST SAME BUT SINGLE PREC.         
C                                                                               
          CIST   = RSTMIN                                                       
          DIST   = DBLE (CIST)                                                  
          QFACT  = QS(IQX)                                                      
          DQFACT = DBLE (QFACT)                                                 

          IF (DEBUGL) THEN                                                      
            WRITE (6,9005)                                                      
            WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,                              
     >        CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,
     >        'ION APPEARED'          
          ENDIF                                                                 

C                                                                               
C-------- FIVE RANDOM NUMBERS ARE ALWAYS USED FOR EACH ITERATION                
C-------- IN LOOP 500.  AN OCCASIONAL EXTRA ONE IS REQUIRED FOR                 
C-------- TESTING RECOMBINATION, FOR ROULETTING OR FOR DETERMINING              
C-------- WHETHER TO APPLY A DELTAY DIFFUSION STEP.  HENCE A MAXIMUM OF         
C-------- EIGHT RANDOMS ARE USED IN EACH ITERATION.  THE RANDOM NUMBERS         
C-------- VECTOR THUS HAS TO BE REGENERATED WHENEVER WE ARE WITHIN 8            
C-------- NUMBERS FROM THE END OF THE VECTOR   (IE. KKLIM).                     
C-------- WHEN WE JUMP OUT OF LOOP 500 A PROPORTION OF THIS VECTOR OF           
C-------- RANDOMS WILL BE WASTED.  TO PREVENT THIS, WE ONLY GENERATE            
C-------- ENOUGH HERE TO REPLACE THOSE USED IN THE LAST PASS THROUGH            
C-------- LOOP 500  (IE THE VALUE OF KK).                                       
C                                                                               
  500     CONTINUE                                                              

c slmod
c          IF (IMP.GT.9.AND.IMP.LT.16) THEN
c            WRITE(99,'(4I4,3F14.8,F14.8)')
c     +         (IMP,CIZ,IX,IY,
c     +         CX,Y,P,
c     +         CTEMI
c          ENDIF
c slmod end
          IF (KK.GT.KKLIM) THEN                                                 
            CALL SURAND (SEED, KK, RANV)                                        
            NRAND = NRAND + KK                                                  
            KK = 0                                                              
          ENDIF                                                                 
c
c         Record particle track if track debugging is turned on.
c
          if (debugt) then 
             tptrac(traclen,1) = cx
             tptrac(traclen,2) = y 
             traclen = traclen + 1
             if (traclen.gt.maxlen) then
                bigtrac = .true.
                traclen = 1
             endif   
c slmod begin - not sure what this does
             IF (TRACLEN.GT.MAXLEN) AVGTRAC = 0.0

             IF (TRACLEN.GT.2) THEN
               AVGTRAC = AVGTRAC + TPTRAC(TRACLEN-1,2) - 
     +                             TPTRAC(TRACLEN-2,2)
             ENDIF
c slmod end
             write(6,'(a,i8,10(1x,g12.5))') 
     >              'trac:',traclen-1,tptrac(traclen-1,1),
     >              tptrac(traclen-1,2),
c slmod
     +             (tptrac(traclen-1,2)-tptrac(traclen-2,2)),
     +             AVGTRAC/(TRACLEN-2)
c     >                  tptrac(traclen-1,2),bigtrac
c slmod end
          endif
C                                                                               
C-------- CHECK FOR CHANGES IN CHARACTERISTIC TIMES DATA. WILL OCCUR            
C-------- WHEN WE USE A SPECIAL PLASMA.  EACH TIME AN ION ENTERS                
C-------- THIS REGION NEW COEFFICIENTS ARE CALCULATED BASED ON ITS              
C-------- TEMPERATURE AT THAT TIME, PROVIDING TEMP IS 10% DIFFERENT             
C-------- FROM THE PREVIOUS TIME THIS REGION WAS ENTERED.                       
C                                                                               


          IF (IX.LE.JX) THEN                                                    
            IF (CIOPTB.GE.2 .OR. CIOPTC.EQ.2 .OR. CIOPTD.EQ.3) THEN             
              TEMOLD = CTOLDS(IX,CIZ)                                           
              IF (CTEMI.GT.1.1*TEMOLD .OR. CTEMI.LT.0.9*TEMOLD) THEN            
                CALL TAUFIX (IX,TEMOLD,CTEMI)                                   
                CTOLDS(IX,CIZ) = CTEMI                                          
C               IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,              
C    >            CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,                 
C    >            'FIX',0,0,TEMOLD                                              
              ENDIF                                                             
            ENDIF                                                               
          ENDIF                                                                 


C                                                                               
C------------ IN MOST CASES, CALCULATE PARALLEL DIFFUSION COEFFICIENT           
C------------ AND MOVE ON.  BUT TO START WITH, EACH ION MUST EXIST FOR          
C------------ "RCONST" ITERATIONS BEFORE PARALLEL DIFFUSION IS APPLIED.         
C------------ THE VALUE OF "RCONST" DEPENDS ON WHICH "FIRST DIFFUSION"          
C------------ OPTION WAS CHOSEN 0,1 OR 2.  ONCE THE ION HAS EXISTED LONG        
C------------ ENOUGH, SET "DIFFUS" FLAG TRUE FOR SUBSEQUENT ITERATIONS.         
C------------ NOTE: CCCFPS = SQRT(4.88E8/(CFPS*CRMI))* QTIM*QS * ...            
C------------ FIXED 14/7/88: IF RCONST < 1  (IE TAUPARA < DELTAT), THEN         
C------------ DIFFUSION SHOULD BE SWITCHED ON STRAIGHT AWAY.                    
C                                                                               
              IF (DIFFUS) THEN                                                  
                SPARA = CTEMI * CCCFPS(IX,IY,CIZ)                               
              ELSE                                                              
                IF (CDIFOP.EQ.2) THEN                                           
                  IF (CFPS(IX,IY,CIZ).GT.0.0) THEN                              
                    RCONST = 2.0 * CTEMI / CFPS(IX,IY,CIZ) * QFACT              
                  ELSE                                                          
                    RCONST = 1.E20                                              
                  ENDIF                                                         
                ENDIF                                                           
                IF (CIST.GE.RCONST .OR. RCONST.LT.1.0) THEN                     
                  RDIFFT = RDIFFT + CIST * QTIM * SPUTY                         
                  DIFFUS = .TRUE.                                               
                  SPARA  = CTEMI * CCCFPS(IX,IY,CIZ)                            
                ELSE                                                            
                  SPARA  = 0.0                                                  
                ENDIF                                                           
              ENDIF                                                             
C                                                                               
              IF (IX.LE.JX) THEN                                                
                IF (CIOPTB.EQ.3) THEN                                           
                  KK = KK + 1                                                   
                  IF (RANV(KK).GT.CFPS(IX,IY,CIZ)/(2.0*CTEMI)) SPARA=0.0        
                ELSEIF (CIOPTB.EQ.4) THEN
                  KK = KK +1
                  IF (RANV(KK).GT.(CFPS(IX,IY,CIZ)/(2.0*CTEMI))) THEN
                    SPARA = 0.0
                  ELSE
                    SPARA = SPARA* SQRT(CTEMI/CTEMSC) 
                    DTEMI = DTEMI + (DBLE(CTEMBSI(IX,IY))-DTEMI) 
     >                    *DMIN1( DBLE(DTEMI/CTEMBSI(IX,IY)),0.5D0)
                    CTEMI = SNGL(DTEMI)
                  ENDIF
                ENDIF                                                           
              ENDIF                                                             
c slmod
              IF (CIOPTB.EQ.13) THEN
c
c               Velocity diffusion:
c
                spara  = 0.0
                vparat = cccfps(ix,iy,ciz)

 7702           nrand = nrand + 1

                call surand(seed,1,ran1)
                if (ran1.eq.0.0) goto 7702
                nrand = nrand + 1

                call surand(seed,1,ran2)
                rgauss = sqrt(-2.0* log(ran1))*cos(2.0*PI*ran2)

                vpara = vparat * rgauss

                SVY = SVY + VPARA * QTIM
              
              ENDIF
c slmod end
C                                                                               
C------------ UPDATE Y POSITION OF ION                                          
C                                                                               
              OLDY  = Y                                                         
              IF (ABSY.LE.CYNEAR .OR. ABSY.GE.CYFAR) THEN                       
                YFACT = CSINTB                                                  
              ELSE                                                              
                YFACT = CSINTB * CLFACT                                         
              ENDIF                                                             
              IF (SVYMIN.EQ.0.0) THEN
                 SVYMOD = SVY
              ELSE
                 SVYMOD = SIGN(MAX(ABS(SVY),SVYMIN),SVY)
                 IF (DEBUGL) 
     >         WRITE(6,*) 'SVYMOD,SVY,SVYMIN:',SVYMOD,SVY,SVYMIN
              ENDIF
C
C             Accumulate velocity information.
C
              svybar(iqx) = svybar(iqx) + abs(svymod * QS(IQX)) * sputy
              svyacc(iqx) = svyacc(iqx) + sputy
C
              DY1 = DY1 + DBLE ((SVYMOD + 0.5 * QUANT)*QS(IQX)*YFACT)             

              IF (SPARA.GT.0.0) THEN                                            
                KK  = KK + 1                                                    
                SPARA = SPARA * YFACT                                           
                DY2 = DY2 + DBLE (SIGN (SPARA,RANV(KK)-0.5))                    
              ENDIF                                                             

c
c             Updating Y coordinate - 
c             NOTE: DY2 contains the initial Y coordinate PLUS all spatial diffusive steps
c                   DY1 contains all forces and velocity diffusive steps   
c
c
              Y     = SNGL (DY1+DY2)                                            
c
c
c             Check for Y-axis reflection 
c
c             cmir_refl_lower < 0
c             cmir_refl_upper > 0
c
              if (yreflection_opt.ne.0) then 

                 call test_reflection(y,oldy,svy)

              endif

c
c slmod begin
              IF (DEBUGL) WRITE(78,'(I4,F7.1,13G12.5)') 
     +          IMP,CIST,
     +          Y,SVY,
     +          SVYMOD,RGAUSS,VPARA*QTIM,VPARAT,DY1,DY2,
     +          QS(IQX),DTEMI,QUANT,CFSS(IX,IY,CIZ),YFACT

c              IF (IONCNT.EQ.10.OR.IONCNT.EQ.20) 
c     +          WRITE(50,*) IONCNT,Y,ABS(Y-OLDY)
c slmod end
              ABSY  = ABS (Y)                                                   
              YY    = MOD (ABSY, CL)                                            
              IF (YY.GT.CHALFL) YY = CL - YY                                    

C                                                                               
C------------ UPDATE X POSITION OF ION, ALLOWING FOR ELONGATION                 
C------------ NOTE YYCON AND ALPHA ARE UPDATED IN INBOARD/OUT LOOP BELOW        
C                                                                               
              OLDALP = ALPHA                                                    
              IF (CX.GE.0.0) THEN                                               
                YYCON = YY*CONI + 1.0                                           
              ELSE                                                              
                YYCON = YY*CONO + 1.0                                           
              ENDIF                                                             
              KK = KK + 1                                                       
C
C             DECIDE WHICH Y-REGION THE PARTICLE IS IN AND THEN USE  
C             THE INDEX TO ACCESS THE X-DIFF DATA FOR THE 
C             APPROPRIATE REGION
C          
C             D.ELDER NOV 23 1990
C
              IF (Y.LE.0.0) THEN                                          
                 IF (Y.GT.-CHALFL) THEN 
                    J = 1
                 ELSEIF (Y.LT.-C3HALFL) THEN 
                    J = 2
                 ELSE 
                    J = 3
                 ENDIF
              ELSE      
                 IF (Y.GT.C3HALFL) THEN 
                    J = 1
                 ELSEIF (Y.LT.CHALFL) THEN 
                    J = 2
                 ELSE 
                    J = 3
                 ENDIF
              ENDIF
C
              IF (CIOPTN.EQ.0) THEN 
                CX = ALPHA * YYCON + CXAFS(IQX,J) +                             
     >                 SIGN (CXBFS(IQX,J),CXCFS(IQX,J)-RANV(KK))                  
              ELSEIF (CIOPTN.EQ.1) THEN 
                CX = ALPHA * YYCON + CXAFS(IQX,J)                              
                KK = KK + 1
                IF (RANV(KK).LT.CXDPS(IQX,J)) THEN 
                   KK = KK +1
                   CX = CX + SIGN(CDPSTP,CXCFS(IQX,J)-RANV(KK))
                ENDIF
              ENDIF                
C                                                                               
C------------ DO NOT NEED THE FOLLOWING TWO LINES FOR THE QUICK STANDARD        
C------------ LIM VERSION WITH NO POLOIDAL DIFFUSION.  THEN KK WILL ONLY        
C------------ BE INCREMENTED 4 TIMES PER LOOP COUNT, BUT THIS IS ALLOWED        
C------------ FOR IN SUBSEQUENT CALLS TO SURAND, WHICH WILL ONLY REPLACE        
C------------ THE KK ACTUAL RANDOM NUMBERS USED.                                
C                                                                               
              IF (BIG) THEN                                                     
                KK = KK + 1                                                     
                P  = P  + SIGN (POLODS(IQX),RANV(KK)-0.5) + SVPOLS(IQX)    
c slmod begin
c                IF (DEBUGL) 
c     +            WRITE(79,*) IQX,IY,IP,P,QTIM,SVPOLS(IQX),POLODS(IQX)
c slmod end
                ABSP = ABS(P) 
              ENDIF                                                             
C                                                                               
C------------ ITERATE CTEMI FOR TEMPERATURE CHANGE                              
C                                                                               
              IF ((CIOPTB.NE.4).OR.(CIOPTB.EQ.4.AND.IX.GT.JX)) THEN
                DTEMI = DTEMI + (DBLE(CTEMBSI(IX,IY))-DTEMI) *                
     >                           DBLE(CFTS(IX,IY,CIZ))                          
                CTEMI = SNGL (DTEMI)                                            
              ENDIF
c slmod begin
              IF (ALPHA.LT.CFTCUT) THEN                                         
                CX    = CFTCUT                                                     
                ALPHA = 0.0                                                   
                IQX   = IPOS (CX, XS, NXS-1)
              ENDIF
c
c Check if ion has penetrated into the DIVIMP grid:
c
              IF (optdp.EQ.1) THEN
                IF (TAG2(IMP).EQ.0.0.AND.ALPHA.GT.MARK) THEN
c
c                 The ion has not entered the DIVIMP grid region yet:
c
                  TAG2(IMP) = 1.0
                  DO II = 1, NBIN                 
                    IF (Y.LT.BSBIN(II)) THEN 
                      NSBIN(II,CIZ) = NSBIN(II,CIZ) + 1 
                      IZBIN(II)     = IZBIN(II)     + CIZ
                      EXIT
                    ENDIF
                  ENDDO
                  TLOSS(IMP) = 0.0 
                ELSEIF (TAG2(IMP).EQ.1.0) THEN
                  IF (ALPHA.LT.MARK) THEN
c
c                   The ion has left the DIVIMP grid region.  Record the
c                   time since leaving the grid to get an estimate for tauFP
c                   to input into DIVIMP:
c
                    TLOSS(IMP) = TLOSS(IMP) + QTIM
                  ELSE            
c
c                   The ion has left the DIVIP grid region and re-entered:
c
                    TLOSS(IMP) = 0.0
                  ENDIF
                ENDIF
              ENDIF
c slmod end

c
c     Add code to check whether the ion has crossed Y=+/-L - if the "shear short circuit" option 
c     is active and the current poloidal position of the particle does not coincide with the limiter -
c     i.e. |P| >CPCO then reset P = (-CPSC,CPSC) randomly distributed. 
c
c     Use oldy and y to determine if the particle cross the chalfl boundaries. 
c
              if (shear_short_circuit_opt.eq.1) then 
                 if (absp.gt.cpco.and.(
     >               (y.lt.-chalfl.and.oldy.gt.-chalfl).or.
     >               (y.gt.-chalfl.and.oldy.lt.-chalfl).or.
     >               (y.lt.chalfl.and.oldy.gt.chalfl).or.
     >               (y.gt.chalfl.and.oldy.lt.chalfl))) then
                    kk = kk+1 
                    p = cpco * (2.0*ranv(kk) -1.0)
c                    write(6,'(a,i10,10g12.5)') 'Shear1:',imp,
c     >                               cx,y,chalfl,oldy,
c     >                               p,cpco,absp
                    absp = abs(p)
c
                 endif
              endif

C                                                                               
C-----------------------------------------------------------------------        
C       ION INBOARD                                                             
C-----------------------------------------------------------------------        
C                                                                               
              IF (CX.GE.0.0) THEN                                               
                YYCON = YY*CONI + 1.0                                           
                ALPHA = CX / YYCON                                              
C                                                                               
C-------------- REFLECT OFF X=A IF REQUIRED; SET IQX POINTER                    
C-------------- ENSURE IQX POINTER NEVER EXCEEDS ARRAY BOUNDS                   
C-------------- (IN CASE OF ROUNDING ERRORS ETC)                                
C                                                                               
                IF (ALPHA.GE.CA) THEN                                           
                  CX    = 2.0 * CA * YYCON - CX                                 
                  ALPHA = CX / YYCON                                            
                  IF (CFLRXA) THEN                                              
                    CICRXA = CICRXA + SPUTY                                     
                    CISRXA = CISRXA + CIST * SPUTY                              
                    CITRXA = CITRXA + CTEMI * SPUTY                             
                    IF (CIST.LT.CIFRXA) CIFRXA = CIST                           
                    CFLRXA = .FALSE.                                            
                  ENDIF                                                         
                ENDIF                                                           
                IQX = MIN (INT(ALPHA*XSCALI)+1, NQXSI)                          
C                                                                               
C-------------- BOUNDARY CONDITION Y>=2L OR Y<=-2L                              
C-------------- IF QTIM IS LARGE ITS POSSIBLE THAT ADDING 2L STILL              
C-------------- LEAVES THE PARTICLE OUTSIDE THE REGION OF INTEREST:             
C-------------- CHECK FOR THIS AND ADD ANOTHER 2L IF NECESSARY ...              
C                                                                               
                IF (Y.LE.-CTWOL) THEN                                           
  401             DY2 = DY2 + 2.0D0 * DWOL                                      
                  Y   = SNGL (DY1+DY2)                                          
                  IF (Y.LE.-CTWOL) GOTO 401                                     
                  ABSY = ABS (Y)                                                
                  IF (CFLY2L) THEN                                              
                    CICY2L = CICY2L + SPUTY                                     
                    IF (CIST.LT.CIFY2L) CIFY2L = CIST                           
                    CFLY2L = .FALSE.                                            
                  ENDIF                                                         
                ELSEIF (Y.GE.CTWOL) THEN                                        
  402             DY2 = DY2 - 2.0D0 * DWOL                                      
                  Y   = SNGL (DY1+DY2)                                          
                  IF (Y.GE.CTWOL) GOTO 402                                      
                  ABSY = ABS (Y)                                                
                  IF (CFLY2L) THEN                                              
                    CICY2L = CICY2L + SPUTY                                     
                    IF (CIST.LT.CIFY2L) CIFY2L = CIST                           
                    CFLY2L = .FALSE.                                            
                  ENDIF                                                         
                ENDIF                                                           
C                                                                               
C-------------- UPDATE ION VELOCITY                                             
C-------------- MAY BE SUBJECT TO BACKGROUND FLOW VELOCITY AND ELECTRIC         
C-------------- FIELDS WHICH ARE SET CONSTANT FOR INBOARD REGION.               
C                                                                               
C--  NOTE 251   IF (Y.GE.0.0) THEN                                              
C-- COMMENT OUT   IF (Y.LT.CL) THEN                                             
C--                 QUANT = SEYINS(IQX,CIZ) +                                   
C--  >                      CFSS(IX,IY,CIZ) * (SVHINS(IQX) - SVY)               
C--               ELSE                                                          
C--------------------------------------------------------------------
C
C   SET ELECTRIC FIELD AND DRIFT VELOCITY EFFECTS TO ZERO OUTSIDE
C   THE SIZE OF THE POLOIDAL EXTENT OF THE LIMITER
C
c       IF ((BIG).AND.(CIOPTJ.EQ.1).AND.(ABSP.GT.CPCO)) THEN
c              QUANT = -CFSS(IX,IY,CIZ)*SVY 
c
c jdemod: I don't understand Steve's comment here - this is required code
c         when using limiters with a limited poloidal extent - though I agree
c         the physics may be incorrect for inboard since it is only applying
c         a frictional force to the particle motion and not including any inboard
c         flows which probabaly should be turned on. However, the code should
c         NOT stop here in any case.
c     
c         For now I will just comment out all of this so that specified poloidal 
c         extent limiters do not affect transport in the confined plasma.
c
c
c slmod begin
c 	      WRITE (0,*) 'Error! Polodal extent.'
c              STOP
c slmod end
c       ELSE
             QUANT =-SEYINS(IQX,CIZ) -                                   
     >           CFSS(IX,IY,CIZ) * (SVHINS(IQX) + SVY)             
c       ENDIF 
c
C--               ENDIF                                                         
C--             ELSE                                                            
C--               IF (Y.GT.-CL) THEN                                            
C--                 QUANT =-SEYINS(IQX,CIZ) -                                   
C--  >                      CFSS(IX,IY,CIZ) * (SVHINS(IQX) + SVY)               
C--               ELSE                                                          
C--                 QUANT = SEYINS(IQX,CIZ) +                                   
C--  >                      CFSS(IX,IY,CIZ) * (SVHINS(IQX) - SVY)               
C--               ENDIF                                                         
C--             ENDIF                                                           
                SVY = SVY + QUANT                                               
C                                                                               
C-----------------------------------------------------------------------        
C       ION OUTBOARD                                                            
C-----------------------------------------------------------------------        
C                                                                               
            ELSE                                                                

              YYCON = YY*CONO + 1.0                                             
              ALPHA = CX / YYCON                                                
C                                                                               
C------------ UPDATE IQX VALUE.                                                 
C                                                                               
              IQX  = MAX (INT(ALPHA*XSCALO), 1-NQXSO)                           
c
c             Calculate an IQY_TMP value to access CAW_QYS which gives the wall distance
c             at a specific value of QYS - this IQY_TMP has a different meaning than the 
c             IQY calculated below (which is relative to the limiter faces). 
c
              if (y.lt.0.0) then 
                 IQY_TMP = max(min(int((y+ctwol)*yscale)+1,nqys),1)
              else
                 IQY_TMP = max(min(int(y*yscale)+1,nqys),1)
              endif

C                                                                               
C------------ NOTE 151,270:  STOPPING CROSS FIELD TRANSPORT.                    
C------------ ORIGINALLY, WHEN X REACHED -4 LAMBDA IT WAS BROUGHT BACK          
C------------ TO 0.  FOR FLEXIBILITY, A CUTOFF IS NOW SPECIFIED IN              
C------------ THE INPUT DATA - CAN BE SWITCHED OFF BY SETTING TO -99.0          
C                                                                               

              IF (ALPHA.LT.CFTCUT) THEN                                         
                CX    = 0.0                                                     
                ALPHA = 0.0                                                     
                IQX   = 0                                                       
C                                                                               
C------------ ION HAS REACHED WALL AND IS ABSORBED                              
C------------ SET TEMPERATURE TO SMALLEST POSSIBLE VALUE                        
C------------ SCORE PARTICLE IN "WALLS" ARRAY (REDIM'D TO -NYS:NYS)             
C                                                                               
c              ELSEIF (ALPHA.LE.CAW) THEN                                        
c
              ELSEIF (ALPHA.LE.CAW_QYS(IQY_TMP)) THEN                                        
                CIAB   = 0                                                      
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(1,IY)                                           
                IF (Y.GT.0.0) THEN                                              
                  CALL MONUP (6,SPUTY)                                          
                  RWALL(2) = RWALL(2) + SPUTY                                   
                ELSE                                                            
                  CALL MONUP (10,SPUTY)                                         
                  RWALL(1) = RWALL(1) + SPUTY                                   
                ENDIF                                                           
                WALLS(IY,CIZ) = WALLS(IY,CIZ) + SPUTY                           
                IFATE = 1                                                       
c slmod begin
c
c Some statistics for the DIVIMP grid source stuff:
c
                IF (optdp.EQ.1) THEN
                  IF (TLOSS(IMP).NE.0.0.AND.IFATE.EQ.1) THEN
                    ALOSS = ALOSS + TLOSS(IMP)
                    WLOSS = WLOSS + 1.0
                  ENDIF
                ENDIF
c slmod end
                GOTO 790                                                        
              ENDIF                                                             

C
C     IF POLOIDAL EXTENT LIMITS ARE IN EFFECT ONLY CHECK FOR COLLISION
C     IF THE P COORDINATE IS INSIDE +/- CPCO.  
C     NOTE: A MAJOR ASSUMPTION IN ALMOST ALL OF THIS CODE IS THAT STEP
C           SIZES ARE SMALL. THUS COLLISIONS OCCUR AT THE APPROXIMATE 
C           POSITIONS OF THE PARTICLES WHEN THEY ARE FOUND TO HAVE HIT 
C           THE LIMITER. FIRST ORDER IS THEN TO RELEASE NEWLY SPUTTERED
C           PARTICLES FROM THE LAST X,Y,P VALUES. SOME REFINING OF THE 
C           X,Y VALUES IS DONE IN HIT AND EDGINT. FOR NOW I WILL LEAVE
C           THE P VALUES UNREFINED.

          IF ((.NOT.BIG).OR.(.NOT.((CIOPTJ.EQ.1).AND.
     >          (ABSP.GT.CPCO)))) THEN

C                                                                               
C------------ ION HAS HIT LIMITER AT Y=0 FROM Y>0 REGION                        
C                                                                               

              EDGE1 = QEDGES(IQX,1)                                             
              EDGE2 = QEDGES(IQX,2)                                             
              IF (OLDY.GT.0.0 .AND. Y.LE.EDGE2) THEN                            
                CIAB   = 1                                                      
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (6,SPUTY)                                            
                IFATE = 2                                                       
                GOTO 780                                                        
C                                                                               
C------------ ION HAS HIT LIMITER AT Y=2L                                       
C                                                                               
              ELSEIF (Y .GE. CTWOL-EDGE1) THEN                                  
                CIAB   = 2                                                      
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (10,SPUTY)                                           
                IFATE = 3                                                       
                GOTO 780                                                        
C                                                                               
C------------ ION HAS HIT LIMITER AT Y=0 FROM Y<0 REGION                        
C                                                                               
              ELSEIF (OLDY.LT.0.0 .AND. Y.GE.-EDGE1) THEN                       
                CIAB   = -1                                                     
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (10,SPUTY)                                           
                IFATE = 4                                                       
                GOTO 780                                                        
C                                                                               
C------------ ION HAS HIT LIMITER AT Y=-2L                                      
C                                                                               
              ELSEIF (Y .LE. EDGE2-CTWOL) THEN                                  
                CIAB   = -2                                                     
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (6,SPUTY)                                            
                IFATE = 5                                                       
                GOTO 780                                                        
              ENDIF                                                             

C           MATCHING ENDIF FOR POLOIDAL EXTENT TEST
            ENDIF  

c slmod tmp
c              WRITE(0,*) 'Error! Ion outboard: ',IMP,CX,Y,ABSP
c              STOP
c slmod end

C                                                                               
C-------------- SET IQY VALUE AND                                               
C-------------- UPDATE ION VELOCITY, AFFECTED BY                                
C-------------- ELECTRIC FIELD AND DRIFT VELOCITY                               
C-------------- ION HAS COMPLETED OUTBOARD STEP   (MONUP1)                      
C                                                                               
c
c           jdemod - WARNING - there is an inconsistency in the 
c                    definition of IQY in the code. IQY is the
c                    index into the underlying QYS and related
c                    arrays. However, it appears that in the code
c                    below IQY was at some point redefined to 
c                    be a number of points between the limiter
c                    surfaces along the field lines - thus the
c                    use of EDGE1, EDGE2 and CYSCLS in the
c                    calculation of IQY here. However, the related
c                    variables that are indexed by IQY were NOT
c                    changed to reflect this usage. 
c                    e.g. CEYS(IQY), CVHYS(IQY) and QYS(IQY) are
c                         all calculated without taking the 
c                         location of the limiter edges into effect.
c                    This means that the IQY value calculated here
c                    does not match properly with these arrays - in
c                    particular QYS. 
c
c                    On the other hand, the 
c                    procedure here will map the IQY range between
c                    limiter surfaces 1:1 onto the CEYS and CVHQYS
c                    arrays thus allowing for a variable limiter
c                    shape still mapping the first data point at IQY=1
c                    to the limiter surface.
c                    The dependencies in CEYS and CVHYS are typically
c                    linear in QYS to CL and so would not change 
c                    much if they were calculated properly for the 
c                    actual limiter surface location. 
c
c                    The biggest concern is QYS - it can not be
c                    properly indexed by this IQY.                     
c
c
            SVG = CALPHE(CIZ) * CTEGS(IX,IY) +
     >            CBETAI(CIZ) * CTIGS(IX,IY) 
c
            IF (Y.GT.0.0) THEN                                              
              IQY   = INT ((Y-EDGE2)  * CYSCLS(IQX)) + 1                    
              IF ((BIG).AND.(CIOPTJ.EQ.1).AND.(ABSP.GT.CPCO)) THEN
                QUANT = -CFSS(IX,IY,CIZ)*SVY
              ELSE 
                QUANT = (CFEXZS(IX,IY,CIZ) * CEYS(IQY)) + SVG +               
     >           (CFSS(IX,IY,CIZ)*(CFVHXS(IX,IY)*CVHYS(IQY)-SVY))  
              ENDIF 
            ELSE                                                            
              IQY   = INT((-Y-EDGE1) * CYSCLS(IQX)) + 1                     
              IF ((BIG).AND.(CIOPTJ.EQ.1).AND.(ABSP.GT.CPCO)) THEN
                QUANT = -CFSS(IX,IY,CIZ)*SVY
              ELSE
                QUANT =-(CFEXZS(IX,IY,CIZ) * CEYS(IQY)) + SVG -              
     >           (CFSS(IX,IY,CIZ)*(CFVHXS(IX,IY)*CVHYS(IQY)+SVY))      
              ENDIF
            ENDIF                                                           
            SVY = SVY + QUANT                                               
            DOUTS(CIZ) = DOUTS(CIZ) + DSPUTY * DQFACT                       
            IF (ALPHA.LT.CRXMIN) CRXMIN = ALPHA                             
          ENDIF                                                             
C                                                                               
C-----------------------------------------------------------------------        
C    BOTH ROUTES CONTINUE HERE.     CHECK FOR COLLISION                         
C    THIS SECTION ACCOUNTS FOR 6% OF CPU TIME AND COULD BE COMMENTED OUT        
C-----------------------------------------------------------------------        
C                                                                               
              KK = KK + 1                                                       
              IF ((CTEMI*RANV(KK)) .LE. CFPS(IX,IY,CIZ)) THEN                   
                CICCOL = CICCOL + SPUTY * QFACT                                 
              ENDIF                                                             
C                                                                               
C-----------------------------------------------------------------------        
C    FIND USER X,Y,P BINS ION LIES WITHIN                                       
C    DO NOT NEED LINES  "430 CONTINUE" THROUGH TO "GOTO 440" FOR                
C    THE QUICKER STANDARD LIM VERSION WITH NO POLOIDAL DIFFUSION.               
C-----------------------------------------------------------------------        
C                                                                               
c
c             jdemod - change structure of these statements to address
c                      intel fortran compiler issue
c
              IF (.NOT.BIG) GOTO 450                                            
  430         CONTINUE                                                          
                IF ((IP.LE.-MAXNPS))GOTO 440                 
                IF ((PS(IP-1).LT.P))GOTO 440                 
                IP = IP - 1                                                     
                GOTO 430                                                        
  440         CONTINUE                                                          
                IF ((IP.GE.MAXNPS)) GOTO 450                 
                IF ((PS(IP).GE.P)) GOTO 450                 
                IP = IP + 1                                                     
                GOTO 440                                                        
  450         CONTINUE                                                          
                IF ((JY.LE.1)) GOTO 460                 
                IF ((YS(JY-1).LT.ABSY)) GOTO 460                 
                JY = JY - 1                                                     
                GOTO 450                                                        
  460         CONTINUE                                                          
                IF ((JY.GE.NYS)) GOTO 470                 
                IF ((YS(JY).GE.ABSY)) GOTO 470                 
                JY = JY + 1                                                     
                GOTO 460                                                        
  470         CONTINUE                                                          
                IF ((IX.LE.1)) GOTO 480                 
                IF ((XS(IX-1).LT.ALPHA)) GOTO 480                 
                IX = IX - 1                                                     
                GOTO 470                                                        
  480         CONTINUE                                                          
                IF ((IX.GE.NXS)) GOTO 490                 
                IF ((XS(IX).GE.ALPHA)) GOTO 490                 
                IX = IX + 1                                                     
                GOTO 480                                                        
  490         CONTINUE                                                          
              IY = JY                                                           
              IF (Y.LT.0.0) IY = -IY                                            
C                                                                               
C-----------------------------------------------------------------------        
C    SCORE PARTICLE IN DDLIMS "NUMBER DENSITY" ARRAY AND IN THE                 
C    CLOUD TEMPERATURES ARRAY DDTS (IN DOUBLE PRECISION)                        
C    ALSO STORE LIM5 COMPONENT, IE TIME DEPENDENT DISTRIBUTIONS, IF             
C    WE HAVE REACHED A GIVEN TIMEPOINT.                                         
C    SOME SECTIONS BELOW ARE ONLY NEEDED FOR 3D AND TIME DEPENDENT              
C    CASES - THEY ARE NOT REQUIRED FOR THE STANDARD VERSION.                    
C    SCORE IN Y STEPSIZES ARRAY THE PARALLEL DIFFUSION COEFFICIENT.             
C-----------------------------------------------------------------------        
C                                                                               
              DDLIMS(IX,IY,CIZ) = DDLIMS(IX,IY,CIZ) + DSPUTY * DQFACT           
c slmod begin
              IF (JY.LE.NY3D.AND.ABS(IP).LT.MAXNPS) then
                 DDLIM3(IX,IY,CIZ,IP) = DDLIM3(IX,IY,CIZ,IP) + 
     +                                 DSPUTY * DQFACT
c              else
c                 write(0,*) 'Warning: IP > MAXNPS or JY > NY3D :',
c     >                      ip,maxnps,jy,ny3d
              endif
c
c              IF (JY.LE.NY3D) DDLIM3(IX,IY,CIZ,IP) =   
c     >           DDLIM3(IX,IY,CIZ,IP) + DSPUTY * DQFACT
c slmod end
c
c slmod begin
c                IF (DEBUGL) WRITE(79,*)
c     +            CX,ALPHA,Y,P,IX,IY,IP,CIZ,DDLIM3(IX,IY,CIZ,IP)
c slmod end
              IF (CIST.GE.CTIMES(IT,CIZ)) THEN                                  
                IF (JY.LE.NY3D)                                                 
     >            LIM5(IX,IY,CIZ,IP,IT) = LIM5(IX,IY,CIZ,IP,IT) + SPUTY         
                IT = IT + 1                                                     
              ENDIF                                                             
              DDTS(IX,IY,CIZ) =DDTS(IX,IY,CIZ)+DSPUTY*   DTEMI   *DQFACT        
              DDYS(IX,IY,CIZ) =DDYS(IX,IY,CIZ)+DSPUTY*DBLE(SPARA)*DQFACT        
C                                                                               
C-----------------------------------------------------------------------        
C    TEST IF IONISATION OR RECOMBINATION OCCURING                               
C    THE SECOND RANDOM NUMBER IS NOT USED IN EVERY ITERATION                    
C    THROUGH LOOP 500, HENCE EXTRA CALL USED FOR SURAND.                        
C-----------------------------------------------------------------------        
C                                                                               
              KK = KK + 1                                                       
              IF (RANV(KK).LE.CPCHS(IX,IY,CIZ)) THEN                            
                KK = KK + 1                                                     
                IF (RANV(KK).LE.CPRCS(IX,IY,CIZ)) THEN                          
C                                                                               
C---------------- EVENT IS A RECOMBINATION.  UPDATE MONITOR VARS                
C---------------- POINTERS ETC.  CHECK FOR C+ --> C EVENT WHICH                 
C---------------- MEANS WE GO ON TO THE NEXT ION.                               
C                                                                               
                  CICRCS(CIZ) = CICRCS(CIZ) + SPUTY                             
                  IF (CIST .LT. CIFRCS(CIZ)) CIFRCS(CIZ) = CIST                 
                  IF (CIST .GT. CILRCS(CIZ)) CILRCS(CIZ) = CIST                 
                  CISRCS(CIZ) = CISRCS(CIZ) + CIST * SPUTY                      
                  CIZ  = CIZ - 1                                                
                  IF (BIG) IT = IPOS (CIST, CTIMES(1,CIZ), NTS)                 
c slmark
                  IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,            
     >              CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,              
     >              'RECOMBINED:',CIZ                                           
                  IF (CIZ.LT.1) THEN                                            
                    TBELOW = TBELOW + SPUTY                                     
                    IFATE = 9                                                   
                    GOTO 790                                                    
                  ENDIF                                                         
                ELSE                                                            
C                                                                               
C---------------- EVENT IS AN IONISATION.    UPDATE MONITOR VARS                
C---------------- POINTERS ETC.  CHECK FOR GOING ABOVE MAXIMUM                  
C---------------- IONISATION STATE - GO ON TO THE NEXT ION.                     
C---------------- RECORD IONISATION POSITION IN "TIZS" ARRAY                    
C---------------- IF SET TI>=TB APPLIES, DO THAT AS WELL.                       
C                                                                               
                  CICIZS(CIZ) = CICIZS(CIZ) + SPUTY                             
                  IF (CIST .LT. CIFIZS(CIZ)) CIFIZS(CIZ) = CIST                 
                  IF (CIST .GT. CILIZS(CIZ)) CILIZS(CIZ) = CIST                 
                  CISIZS(CIZ) = CISIZS(CIZ) + CIST * SPUTY                      
                  TIZS(IX,IY,CIZ) = TIZS(IX,IY,CIZ) + SPUTY                     
                  IF (JY.LE.NY3D)                                               
     >              TIZ3(IX,IY,CIZ,IP) = TIZ3(IX,IY,CIZ,IP) + SPUTY             
                  CIZ  = CIZ + 1                                                
                  IF (BIG) IT = IPOS (CIST, CTIMES(1,CIZ), NTS)                 
                  IF (CIZ.EQ.CIZSET) CTEMI = MAX (CTEMI,CTEMBS(IX,IY))          
                  IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,            
     >              CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,             
     >              'IONISED TO:',CIZ                                           
                  IF (CIZ .GT. NIZS) THEN                                       
                    TBYOND = TBYOND + SPUTY                                     
                    IFATE = 9                                                   
c slmod begin
c
c Check to see if the ion ionised beyond the ionisation state limit before
c it entered the DIVIMP grid.  If this happens a lot, then the validity of the
c ion source for input into DIVIMP would be questionable:
c
                    IF (optdp.EQ.1) THEN
                      IF (TAG2(IMP).NE.1.0) TSLOSS = TSLOSS + 1
                      IZLOSS = IZLOSS + 1.0
                    ENDIF
c slmod end
                    GOTO 790                                                    
                  ENDIF                                                         
                  MAXCIZ = MAX (MAXCIZ, CIZ)                                    
                ENDIF                                                           
              ENDIF                                                             
C                                                                               
C-----------------------------------------------------------------------        
C             SPLITTING PLANE CROSSED - SAVE ION DETAILS, LEAP TO 790           
C-----------------------------------------------------------------------        
C                                                                               
              IF (split.and.ALPHA.GT.CXSPLS(IS)) THEN                                     
                IF (IPUT.EQ.MAXPUT) THEN                                        
                  WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,CX,ALPHA,Y,             
     >              P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'SPLIT FAILED'            
                ELSE                                                            
                  TSPLIT(IS) = TSPLIT(IS) + SPUTY                               
                  NSPLIT(IS) = NSPLIT(IS) + 1                                   
                  IS = IS + 1                                                   
                  IPUT = IPUT + 1                                               
                  MPUT = MAX (MPUT, IPUT)                                       
                  IGET(IPUT) = 1                                                
                  R(1,IPUT) = CX                                                
                  R(2,IPUT) = Y                                                 
                  R(3,IPUT) = P                                                 
                  R(4,IPUT) = QUANT                                             
                  R(5,IPUT) = SVY                                               
                  R(6,IPUT) = CTEMI                                             
                  R(7,IPUT) = CIST                                              
                  R(8,IPUT) = SPUTY                                             
                  R(9,IPUT) = TSTEPL                                            
                  I(1,IPUT) = IQX                                               
                  I(2,IPUT) = IQY                                               
                  I(3,IPUT) = IX                                                
                  I(4,IPUT) = IY                                                
                  I(5,IPUT) = IP                                                
                  I(6,IPUT) = IT                                                
                  I(7,IPUT) = CIZ                                               
                  I(8,IPUT) = MAXCIZ                                            
                  I(9,IPUT) = IS                                                
                  L(1,IPUT) = DIFFUS                                            
                  L(2,IPUT) = CFLRXA                                            
                  L(3,IPUT) = CFLY2L                                            
                  IFATE = 7                                                     
                  GOTO 790                                                      
                ENDIF                                                           
C                                                                               
C-----------------------------------------------------------------------        
C             ROULETTING PLANE CROSSED - DISCARD OR KEEP WITH NEW WEIGHT        
C-----------------------------------------------------------------------        
C                                                                               
              ELSEIF (ALPHA.LT.CXSPLS(IS-1)) THEN                               
                IS = IS - 1                                                     
                KK = KK + 1                                                     
                IF (RANV(KK).GE.CPRUL) THEN                                     
                  IFATE = 8                                                     
                  GOTO 790                                                      
                ELSE                                                            
                  SPUTY  = SPUTY / CPRUL                                        
                  DSPUTY = DBLE (SPUTY)                                         
                  TRULET(IS) = TRULET(IS) + SPUTY                               
                  NRULET(IS) = NRULET(IS) + 1                                   
                  IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,            
     >              CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,               
     >              'ROULETTE KEPT'                                             
                ENDIF                                                           
              ENDIF                                                             
c slmod begin
              IF (optdp.EQ.1.AND.Y.LT.TARGET) THEN
                IF (TLOSS(IMP).NE.0.0) THEN
                  ALOSS  = ALOSS  + TLOSS(IMP)
                  TGLOSS = TGLOSS + 1.0
                ENDIF
c
c Cheating here - not really the wall, may throw off some statistics:
c
                IFATE = 1 
                GOTO 790
              ENDIF
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C       LOOP BACK TO 500 IF ION HAS NOT YET REACHED CUTOFF TIME.                
C       OTHERWISE RECORD THIS FACT WITH MONUP AND STOP ITERATING.               
C       DEBUG PRINTOUT EVERY CSTEPL'TH ITERATION  (EG EVERY 100)                
C       CSTEPL ENTERED AS 0 FOR NO DEBUG OPTION, >0 FOR DEBUG ON.               
C-----------------------------------------------------------------------        
C                                                                               
        IF (DEBUGL) THEN                                                        
          IF (CIST.GE.TSTEPL) THEN                                              
  495       TSTEPL = TSTEPL + CSTEPL                                            
            IF (TSTEPL.LE.CIST) GOTO 495                                        
            WRITE (STRING,'(1X,F10.6,F10.5)') OLDALP,OLDY                       
c slmod
c            WRITE (6,'(I4,F8.3,3I4,5G14.5  )') 
c     +        IMP,CIST,IX,IY,IP,Y,CTEMI,SPARA,CFTS(IX,IY,CIZ),QTIM
            WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,                              
     >        CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,STRING               
c slmod end
          ENDIF                                                                 
        ENDIF                                                                   
C                                                                               
        DIST   = DIST + DQFACT                                                  
        CIST   = SNGL (DIST)                                                    
        QFACT  = QS(IQX)                                                        
        DQFACT = DBLE (QFACT)                                                   
        IF (CIST.LT.CSTMAX) GOTO 500                                            
C                                                                               
        CALL MONUP (7,SPUTY)                                                    
        TCUT = TCUT + SPUTY                                                     
        IFATE = 6                                                               
        GOTO 790                                                                
C                                                                               
C-----------------------------------------------------------------------        
C   ION DEPOSITED ON LIMITER:  STORE DETAILS READY FOR RE-LAUNCH OF             
C   NEUTRAL FRAGMENTS WHEN ALL IONS HAVE BEEN EXHAUSTED.                        
C   EXTRACT EXACT X BIN PARTICLE LIES WITHIN AT THE END OF THIS CURRENT         
C   TIMESTEP  (PRESENT VALUE OF IX INDICATES BIN AT THE START OF THE            
C   CURRENT TIMESTEP, WHICH MIGHT BE FOR X > 0)                                 
C                                                                               
C   BEWARE !!!!! SELF-SPUTTERING COULDN'T BE USED WITH SPLITTING, SINCE         
C   WE WERE OVERWRITING THE ARRAY SPUTYS.  THIS IS ONLY OK SO LONG AS           
C   NPROD <= IMP AT ALL TIMES.  WITH SPLITTING, WE MIGHT HAVE TO RECORD         
C   SAY 3 ABSORPTIONS OF SUB-IONS FOR SOME SPLIT IONS, WHICH WOULD MESS         
C   EVERYTHING UP.  HENCE USE OF SNEWS TEMPORARY ARRAY BELOW.                   
C                                                                               
C   NOTE:  USE OF CX HERE RATHER THAN ALPHA IS DELIBERATE !                     
C   JAN89: USE YMF FACTOR WITH SELF-SPUTTERED NEUTRALS IF FLAG = 1              
C-----------------------------------------------------------------------        
C                                                                               
  780   CONTINUE                                                                
c slmod begin
c
c A little bit of statistics:
c
        IF (optdp.EQ.1) THEN
          IF (TLOSS(IMP).NE.0.0.AND.IFATE.GE.2.AND.IFATE.LE.5) THEN
            ALOSS = ALOSS + TLOSS(IMP)
            LLOSS = LLOSS + 1.0
          ENDIF
        ENDIF
c slmod end
        RMACH = ABS (CVABS/QTIM)                                                
c
c       Impact energy option
c
        if (impact_energy_opt.eq.0) then 
           ENERGY = 3.0 * REAL(CIZ) * CTBIQX +                                     
     >       5.22E-9 * CRMI * RMACH * RMACH + 2.0 * CTEMI                          
        elseif (impact_energy_opt.eq.1) then 
           ENERGY = 3.0 * REAL(CIZ) * CTBIQX + 2.0 * CTEMI                          
        endif
c

        NEROYS(IOY,1) = NEROYS(IOY,1) + SPUTY                                   
        NERODS(IOD,1) = NERODS(IOD,1) + SPUTY                                   
        NERODS3(IOD,IP,1) = NERODS3(IOD,IP,1) + SPUTY                                   
C                                                                               
        KK = KK + 1                                                             
        IF (CIAB.EQ.-1.OR.CIAB.EQ.2) THEN                                       
          DEPS(IX,CIZ,1) = DEPS(IX,CIZ,1) + SPUTY                               
          NEROXS(IX,1,1) = NEROXS(IX,1,1) + SPUTY                               
          RDEP(1)   = RDEP(1) + SPUTY                                     
          IF (CNEUTD.EQ.8)
     >       ENERGY = REAL(CIZ) * CVS(IQX,1)
     >              + 5.22E-9 * CRMI * RMACH * RMACH
     >              + 2.0 * CTEMI
          RYIELD    = YIELD (6, MATLIM, ENERGY,
     >                       ctembs(ix,iy),ctembsi(ix,iy))
     >                      *QMULTS*CYMFSS(IQX,1)           
          RESPUT    = RES (6,MATLIM,RYIELD,.FALSE.,CNEUTD,RANV(KK),
     >                     QMULTS)
          SPUNEW    = SPUTY * RYIELD                                            
C
C         LIMIT THE MAXIMUM VALUE OF A SPUTTERED FRAGMENT TO THE 
C         INPUT VALUE CSPUMAX. THIS CAN BE USED TO CONTROL 
C         SITUATIONS THAT MIGHT DEVELOP INTO A LOCAL RUNAWAY.
C         THE EFFECT CAN BE NEGATED BY MAKING THE VALUE OF 
C         CSPUMAX APPROPRIATELY LARGE. IN MOST CASES A VALUE OF 
C         2.0 WOULD BE ADEQUATE, HOWEVER, A DEFAULT VALUE OF 100.0
C         IS USUALLY ASSIGNED.
C
          SPUNEW = MIN(SPUNEW,CSPUMAX)
C
          YLDTOT(1) = YLDTOT(1) + SPUNEW                                        
          YLDMAX(1) = MAX (YLDMAX(1), SPUNEW)  
        ELSE                                                                    
          DEPS(IX,CIZ,2) = DEPS(IX,CIZ,2) + SPUTY                               
          NEROXS(IX,1,2) = NEROXS(IX,1,2) + SPUTY                               
          RDEP(2)   = RDEP(2) + SPUTY                                           
          IF (CNEUTD.EQ.8)
     >       ENERGY = REAL(CIZ) * CVS(IQX,2)
     >              + 5.22E-9 * CRMI * RMACH * RMACH
     >              + 2.0 * CTEMI
          RYIELD    = YIELD (6, MATLIM, ENERGY,
     >                       ctembs(ix,iy),ctembsi(ix,iy))
     >                      *QMULTS*CYMFSS(IQX,2)           
          RESPUT    = RES (6,MATLIM,RYIELD,.FALSE.,CNEUTD,RANV(KK),
     >                     QMULTS)             
          SPUNEW    = SPUTY * RYIELD                                            
C
C         CSPUMAX - AS ABOVE 
C
          SPUNEW = MIN(SPUNEW,CSPUMAX)
C
          YLDTOT(2) = YLDTOT(2) + SPUNEW                                        
          YLDMAX(2) = MAX (YLDMAX(2), SPUNEW)                                   
        ENDIF                                                                   
        RANDEP = RANDEP + RANV(KK) * SPUTY                                      
C                                                                               
        IF (SPUNEW.GT.CTRESH) THEN                                              
          IF (NPROD.GE.MAXIMP) THEN                                             
            WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,CX,ALPHA,Y,P,
     >        SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'SELFSPUT FAILED'                  
          ELSE                                                                  
            IF (CIAB.EQ.-1.OR.CIAB.EQ.2) THEN                                   
              YTHTOT(1) = YTHTOT(1) + SPUNEW                                    
            ELSE                                                                
              YTHTOT(2) = YTHTOT(2) + SPUNEW                                    
            ENDIF                                                               
            NPROD = NPROD + 1                                                   
            SNEWS(NPROD) = SPUNEW                                               
            IF     (RESPUT) THEN                                                
              RMAXS(NPROD) =-1.0                                                
            ELSEIF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5.OR.               
     >              CNEUTD.EQ.4) THEN                                           
              EMAX = CEMAXF * ENERGY                                            
              RMAXS(NPROD) = 1.0 / ((1.0+CEBD/EMAX) * (1.0+CEBD/EMAX))          
            ELSE                                                                
              RMAXS(NPROD) = 1.0                                                
            ENDIF                                                               
            XPRODS(NPROD) = XM                                                  
            YPRODS(NPROD) = YM                                                  
            PPRODS(NPROD) = P                                                   
          ENDIF                                                                 
        ENDIF                                                                   
C      (GOTO 790)  NOT ACTUALLY NEEDED!                                         
C                                                                               
C-----------------------------------------------------------------------        
C       CURRENT ION OR SUB-ION FINISHED WITH.  UPDATE NUMBER OF IONS            
C       REACHING STATE, ETC.  CARE NEEDED WITH ROULETTED SPUTY VALUES...        
C       SEE IF ANY FURTHER SPLIT IONS EXIST - IF SO RETRIEVE THE DETAILS        
C       OF THE SPLITTING POSITION AND LAUNCH ANOTHER SUB-ION.                   
C       THE RESETTING OF ION DETAILS AT THE POINT OF SPLITTING IS               
C       SOMEWHAT OVER-ENGINEERED BELOW (EG ABSY DOESN'T REALLY NEED             
C       RESETTING), BUT THERE IS NO HARM IN THIS.                               
C-----------------------------------------------------------------------        
C                                                                               
  790   CONTINUE                                                                
        IF (IFATE.NE.7 .AND. IFATE.NE.8) THEN                                   
          CISTOT = CISTOT + CIST * SPUTY                                        
          CISMAX = MAX (CISMAX, CIST)                                           
          DO 792 IZ = CIZSC, MAXCIZ                                             
            RIONS(IZ) = RIONS(IZ) + SPUTY                                       
  792     CONTINUE                                                              
        ENDIF                                                                   
C                                                                               
        IF (DEBUGL) WRITE (6,9003) IMP,CIST+QFACT,IQX,IQY,IX,IY,                
     >    CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,FATE(IFATE)            
c
c       Record particle track to permanent array.)  
c
        if (debugt) then 
c
c          Deal with the last position of the particle
c
           tptrac(traclen,1) = cx
           tptrac(traclen,2) = y 
           traclen = traclen + 1
           if (traclen.gt.maxlen) then
              bigtrac = .true.
              traclen = 1
           endif   
c 
c          Move the data to the permanent array
c 
c          Map the positions to the middle limiter - if necessary
c
           if (bigtrac) then 
c
              if (traclen.eq.1) then
                 if (tptrac(maxlen,2).gt.cl) then
                    do 2040 in = 1,maxlen
                       tptrac(in,2) = tptrac(in,2) - ctwol     
 2040               continue
                 elseif (tptrac(maxlen,2).lt.-cl) then
                    do 2050 in = 1,maxlen
                       tptrac(in,2) = tptrac(in,2) + ctwol     
 2050               continue
                 endif
              elseif (tptrac(traclen-1,2).gt.cl) then
                 do 2060 in = 1,maxlen
                    tptrac(in,2) = tptrac(in,2) - ctwol     
 2060            continue
              elseif (tptrac(traclen-1,2).lt.-cl) then
                 do 2070 in = 1,maxlen
                    tptrac(in,2) = tptrac(in,2) + ctwol     
 2070            continue
              endif
c
              ptracl(imp) = maxlen
c
              do 2000 in = traclen, maxlen
                 ptracs(in-traclen+1,imp,1) = tptrac(in,1)
                 ptracs(in-traclen+1,imp,2) = tptrac(in,2)
 2000         continue
              do 2010 in = 1,traclen-1          
                 ptracs(maxlen-traclen+in+1,imp,1) = tptrac(in,1)
                 ptracs(maxlen-traclen+in+1,imp,2) = tptrac(in,2)
 2010         continue
c
           else
c
              ptracl(imp) = traclen-1
c
              if (tptrac(traclen-1,2).gt.cl) then
                 do 2080 in = 1,traclen -1
                    tptrac(in,2) = tptrac(in,2) - ctwol     
 2080            continue
              elseif (tptrac(traclen-1,2).lt.-cl) then
                 do 2090 in = 1,traclen -1
                    tptrac(in,2) = tptrac(in,2) + ctwol     
 2090            continue
              endif
c
              do 2020 in = 1,traclen-1          
                 ptracs(in,imp,1) = tptrac(in,1)
                 ptracs(in,imp,2) = tptrac(in,2)
 2020         continue
           endif
           write(6,*) 'imp2:',imp,ptracl(imp),bigtrac
     >                   ,traclen
        endif   
C                                                                               
  794   CONTINUE                                                                
        IF (IGET(IPUT).GT.0 .AND. IGET(IPUT).LE.CNSPL) THEN                     
          IF (IGET(IPUT).GT.1) THEN                                             
            CX     = R(1,IPUT)                                                  
            Y      = R(2,IPUT)                                                  
            P      = R(3,IPUT)                                                  
            QUANT  = R(4,IPUT)                                                  
            SVY    = R(5,IPUT)                                                  
            CTEMI  = R(6,IPUT)                                                  
            CIST   = R(7,IPUT)                                                  
            SPUTY  = R(8,IPUT)                                                  
            TSTEPL = R(9,IPUT)                                                  
            IQX    = I(1,IPUT)                                                  
            IQY    = I(2,IPUT)                                                  
            IX     = I(3,IPUT)                                                  
            IY     = I(4,IPUT)                                                  
            IP     = I(5,IPUT)                                                  
            IT     = I(6,IPUT)                                                  
            CIZ    = I(7,IPUT)                                                  
            MAXCIZ = I(8,IPUT)                                                  
            IS     = I(9,IPUT)                                                  
            DIFFUS = L(1,IPUT)                                                  
            CFLRXA = L(2,IPUT)                                                  
            CFLY2L = L(3,IPUT)                                                  
            DY1    = 0.0D0                                                      
            DY2    = DBLE (Y)                                                   
            ABSY   = ABS (Y)                                                    
            YY     = MOD (ABSY, CL)                                             
            IF (YY.GT.CHALFL) YY = CL - YY                                      
            IF (CX.GT.0.0) THEN                                                 
              ALPHA  = CX / (YY*CONI + 1.0)                                     
            ELSE                                                                
              ALPHA  = CX / (YY*CONO + 1.0)                                     
            ENDIF                                                               
            DTEMI  = DBLE (CTEMI)                                               
            DIST   = DBLE (CIST)                                                
            JY     = IABS (IY)                                                  
          ENDIF                                                                 
          SPUTY  = SPUTY / REAL(CNSPL)                                          
          DSPUTY = DBLE (SPUTY)                                                 
          IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,CX,ALPHA,Y,         
     >      P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'SUB-ION:',IGET(IPUT),IPUT       
          IGET(IPUT) = IGET(IPUT) + 1                                           
          GOTO 500                                                              
        ENDIF                                                                   
C                                                                               
        IF (IPUT.GT.0) THEN                                                     
          IPUT = IPUT - 1                                                       
          GOTO 794                                                              
        ENDIF                                                                   
C                                                                               
C-----------------------------------------------------------------------        
C       SEE IF TEST OF CPU TIME USED IS DUE                                     
C       TRAP CASE WHERE TIMUSD=0 OCCURS.                                        
C-----------------------------------------------------------------------        
C                                                                               
        IF (IMP.GE.IMPLIM) THEN                                                 
         TIMUSD = ZA02AS(1) - STATIM                                            
         IF (TIMUSD.GT.0.0) THEN                                                
           PARTIM = (CPULIM-NEUTIM-IONTIM) / TIMUSD                             
         ELSE                                                                   
           PARTIM = 10.0                                                        
         ENDIF                                                                  
         IF (PARTIM.GE.1.05) THEN                                               
           IMPLIM = INT (REAL(IMP) * MIN (4.0,0.25+0.75*PARTIM))                
         ELSE                                                                   
C                                                                               
C--------- HAVE RUN OUT OF CPU TIME, STOP ITERATION                             
C--------- THERE ARE SO MANY COUNTERS IT IS VIRTUALLY IMPOSSIBLE TO             
C--------- WIND UP THE ROUTINE CLEANLY.  JUST WORK OUT HOW MANY IONS            
C--------- HAVE BEEN FOLLOWED IN THE TIME ALLOTTED AND CALL IT A DAY.           
C                                                                               
           IMPLIM = IMP                                                         
           DO 791 J = IMP+1, NATIZ                                              
             RATIZ = RATIZ - SPUTYS(J)                                          
  791      CONTINUE                                                             
           NATIZ = IMP                                                          
           WRITE (6,'('' ERROR:  CPU TIME LIMIT REACHED'')')                    
           WRITE (6,'('' NUMBER OF IONS REDUCED TO'',I5)') NINT(RATIZ)          
           CALL PRB                                                             
           CALL PRC ('ERROR:  CPU TIME LIMIT REACHED')                          
           CALL PRI ('NUMBER OF IMPURITY IONS REDUCED TO ',NINT(RATIZ))         
           CALL PRC ('INCREASE CPU TIME OR INCREASE QUANTUM TIMESTEP')          
           CALL PRC ('RESULTS WILL BE PRINTED BUT THEY SHOULD BE TREATED        
     > WITH CAUTION ...')                                                       
           GOTO 805                                                             
         ENDIF                                                                  
        ENDIF                                                                   
C                                                                               
  800 CONTINUE                                                                  
C                                                                               
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
C                                                                               
  805 CONTINUE                                                                  
C                                                                               
C---- ITERATION OF ALL IONS COMPLETE.  LAUNCH SECONDARY, TERTIARY, ETC          
C---- NEUTRALS AND LEAP BACK TO FOLLOW NEXT LOT OF IONS CREATED.                
C---- YTHTOT AND YLDMAX ARE PRETTY MEANINGLESS IF SPLITTING/ROULETTING          
C---- ARE IN USE SO THERE IS NO POINT IN PRINTING THEM.                         
C---- DON'T BOTHER TO RELAUNCH IF "TOTAL YIELD" IS LESS THAN (SAY) 10.0.        
C                                                                               
      CALL PRB                                                                  
      CALL PRR('AVERAGE X POSITION OF ION INJECTION      ',AVXPOS/RATIZ)        
      CALL PRR('AVERAGE ALPHA POSITION OF ION INJECTION  ',AVAPOS/RATIZ)        
      CALL PRR('AVERAGE ABS(Y) POS OF ION INJECTION      ',AVYPOS/RATIZ)        
      CALL PRR('AVERAGE ABS(P) POS OF ION INJECTION      ',AVPPOS/RATIZ)        
      CALL PRR('AVERAGE TIME TO FIRST DIFFUSION          ',RDIFFT/RATIZ)        
      CALL PRI2('NUMBER OF IONS PLATING OUT ON WALLS',                          
     >                                 NINT(RWALL (1)), NINT(RWALL (2)))        
      CALL PRI2('NUMBER OF IONS PLATING ON LIMITERS ',                          
     >                                 NINT(RDEP  (1)), NINT(RDEP  (2)))        
      IF     (CYMFLG.EQ.-1) THEN                                                
        CALL PRI2('TOTAL SELF-SPUTTERING YIELD        ',                        
     >                                 NINT(YLDTOT(1)), NINT(YLDTOT(2)))        
      ELSE                                                                      
        CALL PRI2('TOTAL SELF-SPUT YIELD  (WITH YMF)  ',                        
     >                                 NINT(YLDTOT(1)), NINT(YLDTOT(2)))        
      ENDIF                                                                     
      CALL PRI2('TOTAL YIELDS ABOVE THRESHOLD YIELD ',                          
     >                                 NINT(YTHTOT(1)), NINT(YTHTOT(2)))        
      CALL PRR2('MAXIMUM YIELD ANY NEUTRAL FRAGMENT ',                          
     >                                              YLDMAX(1),YLDMAX(2))        
      WRITE (6,'(1X,A,I9 )') 'AVERAGE NUMBER OF ITERATIONS PER ION   ',         
     >                                               NINT(CISTOT/RATIZ)         
      WRITE (6,'(1X,A,I13)') 'MAXIMUM NUMBER OF ITERATIONS       ',             
     >                                                     NINT(CISMAX)
      IF ((RDEP(1).NE.0).OR.(RDEP(2).NE.0)) THEN   
        WRITE (6,*) ' AV. RANDOM USED FOR RES ',RANDEP/(RDEP(1)+RDEP(2))
      ENDIF
      
      TATIZ  = TATIZ  + RATIZ                                                   
      TAVXPOS = TAVXPOS + AVXPOS 
      TDEP   = TDEP   + RDEP(1)  + RDEP(2)                                      
      TWALL  = TWALL  + RWALL(1) + RWALL(2)                                     
      TNEUT  = TNEUT  + RNEUT                                                   
      TRES   = TRES   + RRES                                                    
      TWALLN = TWALLN + RWALLN                                                  
      TCENT  = TCENT  + RCENT                                                   
      TTMAX  = TTMAX  + RTMAX                                                   
      TSTRUK = TSTRUK + RSTRUK                                                  
      TFAIL  = TFAIL  + RFAIL                                                   
      IONTIM = IONTIM + ZA02AS (1) - STATIM                                     
C                                                                               
      IF ((CNEUTD.EQ.3.OR.CNEUTD.EQ.4.OR.CNEUTD.GE.5).AND.(NPROD.GT.0)         
     >    .AND.((YTHTOT(1)+YTHTOT(2)).GT.2.0).AND.(STATUS.LT.CMAXGENS)) 
     >  THEN          
        STATUS = STATUS + 1                                                     
        CNEUTA = 0                                                              
        CNEUTB = 2                                                              
        IF (CNEUTD.EQ.4) CNEUTC = 1                                             
        IF (CNEUTE.EQ.1) CNEUTE = 0                                             
        IF (STATUS.LE.50) THEN 
           WRITE(6,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'      
           WRITE(7,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'      
        ELSE
           WRITE (6,9013) '***  LAUNCHING ',STATUS,
     >                    ' GENERATION NEUTRALS  ***'      
           WRITE (7,9013) '***  LAUNCHING ',STATUS,
     >                    ' GENERATION NEUTRALS  ***'      
        ENDIF
        DO 8700 IMP = 1, NPROD                                                  
          SPUTYS(IMP) = SNEWS(IMP)                                              
 8700   CONTINUE                                                                
        CALL LAUNCH (FSRATE,1,NPROD,1,NATIZ,RSTRUK,RRES,                        
     >               RATIZ,RNEUT,RWALLN,RCENT,RTMAX,SEED,NRAND,                 
     >               NEUTIM,RFAIL,STATUS,6,MATLIM,QMULTS)                    
        IF (NATIZ.GT.0) GOTO 200                                                
      ENDIF                                                                     
C                                                                               
      SSEF = 0.0                                                                
      YEFF = 0.0                                                                
      IF (RNEUT1.GT.0.0) SSEF = (TNEUT-RNEUT1) / RNEUT1                         
      IF (GTOT1.GT.0.0)  YEFF = (1.0 + SSEF) * GYTOT1 / GTOT1                   
C                                                                               
  806 CONTINUE                                                                  
      CALL PRB                                                                  
      CALL PRC ('***  S U M M A R Y   D E T A I L S  ***')                      
      CALL PRB                                                                  
      IF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN            
      CALL PRI('NO OF (HIGH) PRIMARY NEUTRALS LAUNCHED   ',                     
     >                                  NINT(RNEUT1-RRES1))                     
      CALL PRI('NO OF (LOW)  PRIMARY NEUTRALS LAUNCHED   ',NINT(RRES1))         
      CALL PRI('NO OF (HIGH) SELF-SPUTTERED NEUTRALS     ',                     
     >                                  NINT(TNEUT-RNEUT1-(TRES-RRES1)))        
      CALL PRI('NO OF (LOW)  SELF-SPUTTERED NEUTRALS     ',                     
     >                                  NINT(TRES-RRES1))                       
      ELSE                                                                      
      CALL PRI('NUMBER OF PRIMARY NEUTRALS LAUNCHED      ',NINT(RNEUT1))        
      CALL PRI('NUMBER OF SELF-SPUTTERED NEUTRALS        ',                     
     >                                  NINT(TNEUT-RNEUT1))                     
      ENDIF                                                                     
      CALL PRI('TOTAL NUMBER OF NEUTRALS LAUNCHED        ',NINT(TNEUT))         
      CALL PRI('TOTAL NO OF NEUTRALS PLATING ON WALLS    ',NINT(TWALLN))        
      CALL PRI('TOTAL NO OF NEUTRALS REACHING CENTRE     ',NINT(TCENT))         
      CALL PRI('TOTAL NO OF NEUTRALS EXISTING AT TMAX    ',NINT(TTMAX))         
      CALL PRI('TOTAL NO OF NEUTRALS STRIKING LIMITER    ',NINT(TSTRUK))        
      CALL PRI('TOTAL NO OF FAILED NEUTRAL LAUNCHES      ',NINT(TFAIL))         
      CALL PRB                                                                  
      CALL PRI('TOTAL NUMBER OF IONS CREATED             ',NINT(TATIZ))         
      CALL PRI('TOTAL NO OF IONS PLATING ON WALLS        ',NINT(TWALL))         
      CALL PRI('TOTAL NO OF IONS PLATING ON LIMITERS     ',NINT(TDEP))          
      CALL PRI('TOTAL NO OF IONS IONISED BEYOND LIMIT    ',NINT(TBYOND))        
      CALL PRI('TOTAL IONS RECOMBINED TO FORM NEUTRALS   ',NINT(TBELOW))        
      CALL PRI('TOTAL NO OF IONS EXISTING AT TMAX        ',NINT(TCUT))          
      CALL PRR('TOTAL AVERAGE X ION INJECTION POSITION   ',
     >          TAVXPOS/TATIZ)
      CALL PRR('SELF-SPUTTERING ENHANCEMENT FACTOR       ',SSEF)                
      CALL PRR('EFFECTIVE YIELD                          ',YEFF)                
C                                                                               
      IF (NSPLIT(1).GT.0) THEN                                                  
        CALL PRB                                                                
        CALL PRB                                                                
        CALL PRC ('*** SPLITTING & ROULETTING ACTIVITY ***')                    
        CALL PRB                                                                
        CALL PRI ('NO OF SUB-IONS TO BE CREATED BY SPLITS   ', CNSPL)           
        CALL PRR ('PROBABILITY OF RETENTION ON ROULETTING   ', CPRUL)           
        CALL PRI ('GREATEST SIZE OF SPLITTING BANK          ', MPUT)            
        WRITE (7,9010)                                                          
        DO 8800 IS = 1, MAXINS                                                  
          IF (NSPLIT(IS).GT.0)                                                  
     >      WRITE (7,9011) IS,CXSPLS(IS),                                       
     >        NINT(TSPLIT(IS)),NSPLIT(IS),NINT(TRULET(IS)),NRULET(IS)           
 8800   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     SCALE ARRAYS DDLIMS, TIZS, ETC   ...                                      
C-----------------------------------------------------------------------        
C                                                                               
      do 3990 iqx = 1-nqxso,nqxsi
        if (svyacc(iqx).gt.0.0) svybar(iqx) = svybar(iqx)/svyacc(iqx)
3990  continue 
C
C---- SET FACTORS ARRAY FOR NORMALISING GRAPHS TO 1 PARTICLE                    
C                                                                               
C
C     SET THE DEFAULT VALUE FOR THESE SCALING CONSTANTS TO 0.0
C     THEN IF THE NUMBER OF NEUTRALS LAUNCHED IS ZERO - AS OCCURS
C     FOR SOME ION INJECTION CASES - ALLOW THE SCALING FACTOR TO 
C     BE SET TO THE INVERSE OF THE NUMBER OF IONS INJECTED.
C
      FACTA(-1) = 0.0                                                           
      IF (RNEUT1.GT.0.0) THEN 
        FACTA(-1) = 1.0 / RNEUT1                               
      ELSEIF (TATIZ.GT.0.0) THEN
        FACTA(-1) = 1.0 / TATIZ
      ENDIF
      FACTB(-1) = FACTA(-1) * FSRATE                                            
      FACTA(0) = 0.0                                                            
      IF (TNEUT.GT.0.0) THEN 
         FACTA(0) = 1.0 / TNEUT                                  
      ELSEIF (TATIZ.GT.0.0) THEN 
         FACTA(0) = 1.0 / TATIZ
      ENDIF  
      FACTB(0) = FACTA(0) * FSRATE                                              
C                                                                               
      IF (NIZS.GT.0) THEN                                                       
        DO 4000 IZ = 1, NIZS                                                    
          IF (TATIZ.GT.0.0) THEN                                                
            FACTA(IZ) = 1.0 / TATIZ                                             
          ELSE                                                                  
            FACTA(IZ) = 0.0                                                     
          ENDIF                                                                 
          FACTB(IZ) = FACTA(IZ) * QTIM                                          
 4000   CONTINUE                                                                
      ENDIF                                                                     
      WRITE (6,*) '1:'
C                                                                               
C-----------------------------------------------------------------------        
C     SECTION FOR SYMMETRIC CONTRIBUTIONS.  NOTES 75,219                        
C-----------------------------------------------------------------------        
C                                                                               
      DO 4130 IZ = -1, NIZS                                                     
       DO 4120 IX = 1, NXS                                                      
        DO 4100 IY = 1, NYS                                                     
          IF (IZ.GT.0) THEN                                                     
            DEMP( IY,1) = DDTS(IX, IY,IZ) + DDTS(IX,-NYS-1+IY,IZ)               
            DEMP(-IY,1) = DDTS(IX,-IY,IZ) + DDTS(IX, NYS+1-IY,IZ)               
            DEMP( IY,2) = DDYS(IX, IY,IZ) + DDYS(IX,-NYS-1+IY,IZ)               
            DEMP(-IY,2) = DDYS(IX,-IY,IZ) + DDYS(IX, NYS+1-IY,IZ)               
          ENDIF                                                                 
          DEMP( IY,3) = DBLE (TIZS(IX, IY,IZ) + TIZS(IX,-NYS-1+IY,IZ))          
          DEMP(-IY,3) = DBLE (TIZS(IX,-IY,IZ) + TIZS(IX, NYS+1-IY,IZ))          
          DEMP( IY,4) = DDLIMS(IX, IY,IZ) + DDLIMS(IX,-NYS-1+IY,IZ)             
          DEMP(-IY,4) = DDLIMS(IX,-IY,IZ) + DDLIMS(IX, NYS+1-IY,IZ)             
 4100   CONTINUE                                                                
        DO 4110 IY = -NYS, NYS                                                  
          IF (IZ.GT.0) THEN                                                     
            DDTS(IX,IY,IZ) = DEMP(IY,1)                                         
            DDYS(IX,IY,IZ) = DEMP(IY,2)                                         
          ENDIF                                                                 
          TIZS(IX,IY,IZ) = SNGL(DEMP(IY,3))                                     
          DDLIMS(IX,IY,IZ) = DEMP(IY,4)                                         
 4110   CONTINUE                                                                
 4120  CONTINUE                                                                 
 4130 CONTINUE                                                                  
      WRITE(6,*) '2:'
C                                                                               
C====================== DDTS AND DDYS ARRAYS ===========================        
C                                                                               
C---- CALCULATE STEADY STATE CLOUD TEMPERATURES AND AVERAGE DELTAY STEPS        
C---- THE VALUES MUST BE STORED IN D.P. TO GET ENOUGH ACCURACY,                 
C---- AND THE SUMMATION IS ALSO DONE IN D.P..                                   
C                                                                               
      DO 4290 IZ = 1, NIZS                                                      
C                                                                               
        DO 4210 IX = 1, NXS                                                     
          DSUM1 = 0.0D0                                                         
          DSUM2 = 0.0D0                                                         
          DSUM3 = 0.0D0                                                         
          DO 4200 IY = -NYS, NYS                                                
            DSUM1 = DSUM1 + DDTS  (IX,IY,IZ)                                    
            DSUM2 = DSUM2 + DDYS  (IX,IY,IZ)                                    
            DSUM3 = DSUM3 + DDLIMS(IX,IY,IZ)                                    
 4200     CONTINUE                                                              
          IF (DSUM3.GT.0.0D0) THEN                                              
            SDTXS(IX,IZ) = DSUM1 / DSUM3                                        
            SDYXS(IX,IZ) = DSUM2 / DSUM3                                        
          ENDIF                                                                 
 4210   CONTINUE                                                                
C                                                                               
        DO 4230 IY = -NYS, NYS                                                  
          DSUM1 = 0.0D0                                                         
          DSUM2 = 0.0D0                                                         
          DSUM3 = 0.0D0                                                         
          DO 4220 IX = 1, NXS                                                   
            DSUM1 = DSUM1 + DDTS  (IX,IY,IZ)                                    
            DSUM2 = DSUM2 + DDYS  (IX,IY,IZ)                                    
            DSUM3 = DSUM3 + DDLIMS(IX,IY,IZ)                                    
 4220     CONTINUE                                                              
          IF (DSUM3.GT.0.0D0) THEN                                              
            SDTYS(IY,IZ) = DSUM1 / DSUM3                                        
            SDYYS(IY,IZ) = DSUM2 / DSUM3                                        
          ENDIF                                                                 
 4230   CONTINUE                                                                
C
C      THE CENTRAL VALUE OF THE AVERAGE TEMPERATURE Y ARRAY IS 
C      NOT SIGNIFICANT BECAUSE THE CENTRAL BIN IS NOT USED. BUT 
C      IT DOES SHOW UP IN PLOTTING. THEREFORE, SET IT TO THE AVERAGE
C      OF THE +1 AND -1  VALUES.
C
C      D. ELDER APRIL 6, 1990
C
        SDTYS(0,IZ) = (SDTYS(1,IZ)+SDTYS(-1,IZ))/2.0
C                                                                               
        DSUM1 = 0.0D0                                                           
        DSUM2 = 0.0D0                                                           
        DSUM3 = 0.0D0                                                           
        DSUM4 = 0.0D0                                                           
        DO 4250 IY = -NYS, NYS                                                  
          DO 4240 IX = 1, NXS                                                   
            DSUM1 = DSUM1 + DDTS  (IX,IY,IZ)                                    
            DSUM3 = DSUM3 + DDLIMS(IX,IY,IZ)                                    
            IF (XS(IX).LE.0.0) THEN                                             
              DSUM2 = DSUM2 + DDYS  (IX,IY,IZ)                                  
              DSUM4 = DSUM4 + DDLIMS(IX,IY,IZ)                                  
            ENDIF                                                               
 4240     CONTINUE                                                              
 4250   CONTINUE                                                                
        IF (DSUM3.GT.0.0D0) SDTZS(IZ) = DSUM1 / DSUM3                           
        IF (DSUM4.GT.0.0D0) SDYZS(IZ) = DSUM2 / DSUM4                           
C                                                                               
        DO 4270 IY = -NYS, NYS                                                  
          DO 4260 IX = 1, NXS                                                   
            IF (DDLIMS(IX,IY,IZ).GT.0.0D0) THEN                                 
              DDTS(IX,IY,IZ) = DDTS(IX,IY,IZ) / DDLIMS(IX,IY,IZ)                
              DDYS(IX,IY,IZ) = DDYS(IX,IY,IZ) / DDLIMS(IX,IY,IZ)                
            ENDIF                                                               
 4260     CONTINUE                                                              
 4270   CONTINUE                                                                
 4290 CONTINUE                                                                  
      WRITE(6,*) '3:'
C                                                                               
C
C
C
C                                                                               
C======================== TIZS AND TIZ3 ARRAYS =========================        
C                                                                               
      DO 4430 IZ = -1, NIZS                                                     
       DO 4420 IX = 1, NXS                                                      
        DO 4410 IY = 1, NYS                                                     
          FACT = FACTA(IZ) /                                                    
     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY))             
          TIZS(IX, IY,IZ) = FACT * TIZS(IX, IY,IZ)                              
          TIZS(IX,-IY,IZ) = FACT * TIZS(IX,-IY,IZ)                              
          IF (IY.LE.NY3D) THEN                                                  
            DO 4400 IP = -MAXNPS, MAXNPS                                        
              TIZ3(IX, IY,IZ,IP) = FACT * TIZ3(IX, IY,IZ,IP)                    
              TIZ3(IX,-IY,IZ,IP) = FACT * TIZ3(IX,-IY,IZ,IP)                    
 4400       CONTINUE                                                            
          ENDIF                                                                 
 4410   CONTINUE                                                                
 4420  CONTINUE                                                                 
 4430 CONTINUE                                                                  
      WRITE(6,*) '4:'
C                                                                               
C======================== LIM5 ARRAY ===================================        
C                                                                               
      IF (IMODE.NE.2) THEN                                                      
       DO 4540 IZ = -1, NIZS                                                    
        DO 4530 IX = 1, NXS                                                     
         DO 4520 IY = 1, NY3D                                                   
          FACT = FACTA(IZ) /                                                    
     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY))             
          DO 4510 IT = 1, NTS                                                   
           DO 4500 IP = -MAXNPS, MAXNPS                                         
            LIM5(IX, IY,IZ,IP,IT) = FACT * LIM5(IX, IY,IZ,IP,IT)                
            LIM5(IX,-IY,IZ,IP,IT) = FACT * LIM5(IX,-IY,IZ,IP,IT)                
 4500      CONTINUE                                                             
 4510     CONTINUE                                                              
 4520    CONTINUE                                                               
 4530   CONTINUE                                                                
 4540  CONTINUE                                                                 
      ENDIF                                                                     
      WRITE(6,*) '5:'
C                                                                               
c slmod begin - total volume
      DO IY = -NYS, NYS
        DO IX = 1, NXS
          IF (IY.NE.0) TOTVOL = TOTVOL + YWIDS(ABS(IY)) * XWIDS(IX)
        ENDDO
      ENDDO
c slmod end
c
C===================== DDLIMS AND DDLIM3 ARRAYS ========================        
C                                                                               
C---- DEAL WITH NEUTRALS FIRST:  SCALE BY BIN WIDTHS, NEUT TIMESTEP             
C---- AND NO. OF NEUTRALS LAUNCHED;                                             
C---- THEN DEAL WITH IONS:  SCALE BY BIN WIDTHS,                                
C---- LIM3 TIMESTEP AND NO. OF IONS FOLLOWED.                                   
C                                                                               
      DO 4630 IZ = -1, NIZS                                                     
       DO 4620 IX = 1, NXS                                                      
        DO 4610 IY = 1, NYS                                                     
          DACT = DBLE (FACTB(IZ) /                                              
     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY)))            
          DDLIMS(IX, IY,IZ) = DACT * DDLIMS(IX, IY,IZ)                          
          DDLIMS(IX,-IY,IZ) = DACT * DDLIMS(IX,-IY,IZ)                          
          IF (IY.LE.NY3D) THEN                                                  
            DO 4600 IP = -MAXNPS, MAXNPS                                        
              DDLIM3(IX, IY,IZ,IP) = DACT * DDLIM3(IX, IY,IZ,IP)                
              DDLIM3(IX,-IY,IZ,IP) = DACT * DDLIM3(IX,-IY,IZ,IP)                
c slmod tmp
              IF (DEBUGL) WRITE(79,*) 
     +          IX,IY,IZ,IP,DDLIM3(IX,IY,IZ,IP),DACT,XWIDS(IX),
     +          YWIDS(IY),XCYLS(IX),DELPS(IX,IY)
c slmod end
 4600       CONTINUE                                                            
          ENDIF                                                                 
 4610   CONTINUE                                                                
 4620  CONTINUE                                                                 
 4630 CONTINUE                                                                  
      WRITE(6,*) '6:'
c slmod begin - density integration over all space
      IF (BIG)THEN
        DO IZ = 1, NIZS              
        
          TOTDEN(IZ) = 0.0

          DO IX = 1, NXS
            DO IY = 1, NYS
              DO IP = -MAXNPS, MAXNPS
                TOTDEN(IZ) = TOTDEN(IZ) + 
     +            DDLIM3(IX, IY,IZ,IP) * XWIDS(IX) * YWIDS(IY)
                TOTDEN(IZ) = TOTDEN(IZ) + 
     +            DDLIM3(IX,-IY,IZ,IP) * XWIDS(IX) * YWIDS(IY)
              ENDDO
            ENDDO
          ENDDO
    
          TOTDEN(IZ) = TOTDEN(IZ) / TOTVOL

        ENDDO
      ENDIF
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C      ANALYTIC EXTENSION NOTE 282 ... APPLIES TO DDLIMS ARRAY ONLY             
C      (BUT THIS THEN AFFECTS POWLS, LINES AND PLRPS ARRAYS BELOW).             
C-----------------------------------------------------------------------        
C                                                                               
       IC = IPOS (CANAL, XS, NXS-1)                                             
       IF (IC.LT.NXS) THEN                                                      
C                                                                               
         CALL RZERO (SAVES, MAXNXS*(MAXIZS+4))                                  
         DO 1890 IZ = -1, NIZS                                                  
           DO 1880 IX = 1, NXS                                                  
             DSUM1 = 0.0D0                                                      
             DO 1870 IY = 1, NYS/2                                              
               DSUM1 = DSUM1 + DBLE(YWIDS(IY)) *                                
     >                 (DDLIMS(IX,IY,IZ) + DDLIMS(IX,-IY,IZ))                   
 1870        CONTINUE                                                           
             SAVES(IX,IZ) = SNGL(DSUM1)                                         
             IF (IZ.GT.0) SAVES(IX,NIZS+1)=SAVES(IX,NIZS+1)+SNGL(DSUM1)         
 1880      CONTINUE                                                             
 1890    CONTINUE                                                               
C                                                                               
         DSUM1 = 0.0D0                                                          
         DO 1910 IZ = 1, NIZS                                                   
           DO 1900 IY = 1, NYS                                                  
             DSUM1 = DSUM1 + DBLE(YWIDS(IY)) *                                  
     >               (DDLIMS(IC,IY,IZ) + DDLIMS(IC,-IY,IZ))                     
 1900      CONTINUE                                                             
 1910    CONTINUE                                                               
         DSUM1 = DSUM1 / (2.0D0 * DWOL) * DBLE (EXP (CVIN *                     
     >         ((CA-XOUTS(IC)) * (CA-XOUTS(IC)) / (CA*CA)) - 1.0))              
         WRITE (6,'('' LIM3: ANLY EXT IC,DSUM1'',I5,G12.4)') IC,DSUM1           
C                                                                               
         DO 1940 IX = IC+1, NXS                                                 
           DO 1930 IY = -NYS, NYS                                               
             DO 1920 IZ = 1, NIZS-1                                             
               DDLIMS(IX,IY,IZ) = 0.0D0                                         
 1920        CONTINUE                                                           
             DDLIMS(IX,IY,NIZS) = DSUM1 * DBLE (EXP (CVIN *                     
     >         (1.0 - (CA-XOUTS(IX)) * (CA-XOUTS(IX)) / (CA*CA))))              
 1930      CONTINUE                                                             
 1940    CONTINUE                                                               
       ENDIF                                                                    
       WRITE(6,*) '7:'
C                                                                               
C-----------------------------------------------------------------------        
C     SCALE DEPOSITION AND NET EROSION QUANTITIES                               
C-----------------------------------------------------------------------        
C                                                                               
C---- SCALE DEPS ARRAY BY DIVIDING THROUGH BY BIN-WIDTHS AND THE TOTAL          
C---- NUMBER OF ABSORBED PARTICLES AT Y=0,2L,-2L                                
C---- SET DEPS(,,3) TO TOTAL FOR BOTH SIDES OF Y=0                              
C                                                                               
C---- NOTE : THE FOLLOWING IS SCALED BY DIVIDING BY THE TOTAL 
C---- NUMBER OF NEUTRALS, HOWEVER FOR ION INJECTION CASES THE 
C---- TOTAL NUMBER OF NEUTRALS IS ZERO AND THE FACTA(0) 
C---- VALUE IS SET TO ZERO. IN ORDER TO KEEP THE DEPOSITION DATA
C---- AVAILABLE FOR PLOTTING, TEST IF TNEUT=0 AND THEN SET THE 
C---- MULTIPLYING FACTOR TO 1.0 FOR THESE CASES. THIS COULD BE
C---- CHANGED TO THE TOTAL NUMBER OF IONS INJECTED.
C
      FACTDEPS = FACTA(0)
      IF (FACTDEPS.EQ.0.0) FACTDEPS=1.0
      IF (NIZS.GT.0) THEN                                                       
        DO 881 IZ = 1, NIZS                                                     
          DO 880 IX = 1, NXS                                                    
            DEPS(IX,IZ,1) = DEPS(IX,IZ,1) / XWIDS(IX) * FACTDEPS                
            DEPS(IX,IZ,2) = DEPS(IX,IZ,2) / XWIDS(IX) * FACTDEPS                
            DEPS(IX,IZ,3) = DEPS(IX,IZ,1) + DEPS(IX,IZ,2)                       
            DEPS(IX,NIZS+1,1) = DEPS(IX,NIZS+1,1) + DEPS(IX,IZ,1)
            DEPS(IX,NIZS+1,2) = DEPS(IX,NIZS+1,2) + DEPS(IX,IZ,2)
            DEPS(IX,NIZS+1,3) = DEPS(IX,NIZS+1,3) + DEPS(IX,IZ,3)
  880     CONTINUE                                                              
  881   CONTINUE                                                                
      ENDIF                                                                     
      WRITE(6,*) '8:'
C                                                                               
C---- SCALE REMOVAL RATES AS STORED IN NEROXS (,2) AND (,3) LOCATIONS           
C---- BY BIN WIDTHS AND TO "1 ATOM LAUNCHED"                                    
C---- CALCULATE TOTAL DEPOSITION, NET EROSION AND NENNL (NOTE 289)              
C---- CALCULATE SUM OVER BOTH SIDES OF Y=0 AND STORE THIS ALSO                  
C                                                                               
      FACT = 0.0                                                                
      IF (TDEP.GT.0.0) FACT = TNEUT / TDEP                                      
      DO 883 J = 1, 2                                                           
       DO 883 IX = 1, NXS                                                       
        NEROXS(IX,1,J) =-NEROXS(IX,1,J) / XWIDS(IX) * FACTA(0)                  
        NEROXS(IX,2,J) = NEROXS(IX,2,J) / XWIDS(IX) * FACTA(-1)                 
        NEROXS(IX,3,J) = NEROXS(IX,3,J) / XWIDS(IX) * FACTA(0)                  
        NEROXS(IX,4,J) = NEROXS(IX,1,J) + NEROXS(IX,3,J)                        
        NEROXS(IX,5,J) = FACT * NEROXS(IX,1,J) + NEROXS(IX,3,J)                 
  883 CONTINUE                                                                  
      DO 884 IX = 1, NXS                                                        
       DO 884 II = 1, 5                                                         
        NEROXS(IX,II,3) = NEROXS(IX,II,1) + NEROXS(IX,II,2)                     
  884 CONTINUE                                                                  
      WRITE(6,*) '9:'
C                                                                               
C---- SCALE REMOVAL RATES IN NEROYS,NERODS ARRAYS ...                           
C                                                                               
      WRITE(6,*) 'NER:',FACTA(0),FACTA(-1)
      DO 885 IO = 1, MAXOS                                                      
       IF (OYWIDS(IO).GT.0.0) THEN
        NEROYS(IO,1) =-NEROYS(IO,1) / OYWIDS(IO) * FACTA(0)                     
        NEROYS(IO,2) = NEROYS(IO,2) / OYWIDS(IO) * FACTA(-1)                    
        NEROYS(IO,3) = NEROYS(IO,3) / OYWIDS(IO) * FACTA(0)                     
       ENDIF
        NEROYS(IO,4) = NEROYS(IO,1) + NEROYS(IO,3)                              
        NEROYS(IO,5) = FACT * NEROYS(IO,1) + NEROYS(IO,3)                       
       IF (ODWIDS(IO).GT.0.0) THEN
        NERODS(IO,1) =-NERODS(IO,1) / ODWIDS(IO) * FACTA(0)                     
        NERODS(IO,2) = NERODS(IO,2) / ODWIDS(IO) * FACTA(-1)                    
        NERODS(IO,3) = NERODS(IO,3) / ODWIDS(IO) * FACTA(0)                     
c
c       Need to scale by the 3D bin width as well taking into account any 
c       limiter poloidal extent.
c
        do ip = -maxnps,maxnps
c
c          Calculate local_pwid - to convert erosion to proper surface density
c           

           pbnd1=min(abs(ps(ip)-pwids(ip)/2.0),
     >               abs(ps(ip)+pwids(ip)/2.0))
           pbnd2=max(abs(ps(ip)-pwids(ip)/2.0),
     >               abs(ps(ip)+pwids(ip)/2.0))

           if (cioptj.eq.1.and.
     >        (cpco.gt.abs(pbnd1).and.cpco.lt.abs(pbnd2))) then
              local_pwid = cpco-abs(pbnd1)
           else
              local_pwid = pwids(ip)
           endif

           NERODS3(IO,IP,1) =-NERODS3(IO,IP,1) / ODWIDS(IO) 
     >                          / local_pwid * FACTA(0)                     
           NERODS3(IO,IP,2) = NERODS3(IO,IP,2) / ODWIDS(IO) 
     >                          / local_pwid * FACTA(-1)                    
           NERODS3(IO,IP,3) = NERODS3(IO,IP,3) / ODWIDS(IO) 
     >                          / local_pwid * FACTA(0)                     
        end do
       ENDIF
        NERODS(IO,4) = NERODS(IO,1) + NERODS(IO,3)                              
        NERODS(IO,5) = FACT * NERODS(IO,1) + NERODS(IO,3)                       
        do ip = -maxnps,maxnps
           NERODS3(IO,IP,4) = NERODS3(IO,IP,1) + NERODS3(IO,IP,3)                              
           NERODS3(IO,IP,5) = FACT * NERODS3(IO,IP,1) + NERODS3(IO,IP,3)                       
        end do
  885 CONTINUE                                                                  

C                                                                               
C---- SCALE WALL DEPOSITION ARRAY TO "1 ATOM LAUNCHED"                          
C---- SET SECONDARY NEUTRALS LINE AND TOTALS LINE.                              
C                                                                               
      DO 889 IY = -NYS, NYS                                                     
        IF (IY.EQ.0) GOTO 889                                                   
        WALLS(IY,-2) = WALLS(IY,0) - WALLS(IY,-1)                               
        DO 888 IZ = -2, NIZS                                                    
          WALLS(IY,IZ) = WALLS(IY,IZ) / YWIDS(IABS(IY)) * FACTA(0)              
          IF (IZ.GT.0) WALLS(IY,NIZS+1)= WALLS(IY,NIZS+1) + WALLS(IY,IZ)        
  888   CONTINUE                                                                
  889 CONTINUE                                                                  
      WRITE(6,*) '11:'
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE RADIATIVE POWER LOSS AND LINE RADIATION                         
C-----------------------------------------------------------------------        
C                                                                               
C---- ZERO ARRAYS.                                                              
C                                                                               
      CALL RZERO (POWLS,  MAXNXS*(2*MAXNYS+1)*(MAXIZS+2))                       
      CALL RZERO (LINES,  MAXNXS*(2*MAXNYS+1)*(MAXIZS+2))                       
C     CALL RZERO (POWL3,  MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
C     CALL RZERO (LINE3,  MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
C                                                                               
C-----------------------------------------------------------------------        
C   FOR 
C-----------------------------------------------------------------------        
C                                                                               
      if (cdatopt.eq.0) then 

      DO 990 IY = -NYS, NYS                                                     
        IF (IY.EQ.0) GOTO 990                                                   
        JY = IABS (IY)                                                          
C                                                                               
C---- CALCULATE PLASMA TEMPERATURE AND ELECTRON DENSITY AT MID POINTS           
C---- OF EACH X BIN.  NOTE THIS INVOLVES A CONVERSION FROM THE REGULAR          
C---- SPACED QXS MESH TO THE USER SUPPLIED XS MESH; AND A CONVERSION            
C---- FROM M**3 TO CM**3 FOR NOCORONA.                                          
C---- THE IQX --> IX INDICES ARE TAKEN FROM COMMON /COMXYT/                     
C                                                                               

      DO 900 IX = 1, NXS                                                        
        PTES(IX) = CTEMBS(IX,IY)                                                
        PNES(IX) = CRNBS(IX,IY) * 1.E-6 * REAL (CIZB)                           
  900 CONTINUE                                                                  
C                                                                               
C                                                                               
C------ CALCULATE BIN VOLUMES CM**3 (ASSUME 1 METRE IN THIRD DIMENSION)         
C------ TRANSFER IMPURITY ION DENSITIES TO PNZS ARRAY FOR NOCORONA              
C------ CONVERTING TO CM**-3                                                    
C                                                                               
        DO 920 IX = 1, NXS                                                      
          PDVOLS(IX) = 1.0E6 * XWIDS(IX) * YWIDS(JY)                            
          DO 910 IZ = 0, NIZS                                                   
            PNZS(IZ+1,1,IX) = 1.0E-6 * SNGL(DDLIMS(IX,IY,IZ))                   
  910     CONTINUE                                                              
  920   CONTINUE                                                                
C                                                                               
C------ CALL ROUTINE FROM NOCORONA PACKAGE TO CALCULATE RADIATIVE               
C------ POWER LOSS (W CM**3) AND LINE RADIATION LOSS (W CM**3)                  
C------ THESE ARE CONVERTED TO UNITS (W M**3)                                   
C------ VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY              
C------ GOES OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.         
C                                                                               
        CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NXS)                          
C                                                                               
C------ COPY INTO "POWLS,LINES" ARRAY:  POWER LOSS; LINE RADIATION LOSS         
C------ WE KEEP TRACK OF GRAND TOTALS FOR POWER LOSS, POWER LOSS IN SOL,        
C------ LINE RADIATION LOSS & LINE RADIATION IN SOL IN DTOTS ARRAY BELOW        
C                                                                               
        DO 940 IX = 1, NXS                                                      
          DO 930 IZ = 0, NIZS                                                   
            POWLS(IX,IY,IZ) = MAX (0.0, PRADIS(1,IZ+1,1,IX)*1.0E6)              
            LINES(IX,IY,IZ) = MAX (0.0, PRADIS(3,IZ+1,1,IX)*1.0E6)              
c           write(6,'(a,3i6,8(1x,g12.5))') 'POWLS:',ix,iy,iz,
c    >         pradis(1,iz+1,1,ix)*1e6,ddlims(ix,iy,iz),
c    >         ptes(ix),pnes(ix),pnzs(iz+1,1,ix),pdvols(ix)

  930     CONTINUE                                                              
  940   CONTINUE                                                                
C                                                                               
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1) LOCATIONS             
C                                                                               
        DO 960 IX = 1, NXS                                                      
          IF (DDLIMS(IX,IY,0).LE.0.0D0) THEN                                    
            POWLS(IX,IY,-1) = 0.0                                               
            LINES(IX,IY,-1) = 0.0                                               
          ELSE                                                                  
            POWLS(IX,IY,-1) = POWLS(IX,IY,0) / SNGL(DDLIMS(IX,IY,0)) *          
     >                                         SNGL(DDLIMS(IX,IY,-1))           
            LINES(IX,IY,-1) = LINES(IX,IY,0) / SNGL(DDLIMS(IX,IY,0)) *          
     >                                         SNGL(DDLIMS(IX,IY,-1))           
          ENDIF                                                                 
  960   CONTINUE                                                                
C                                                                               
C------ EXTRA SECTION FOR 3D ARRAYS POWL3 AND LINE3 ...                         
C                                                                               
C
C       MEMORY RESTRICTIONS ON THE CRAY MAKE IT IMPOSSIBLE TO RUN 3D
C       CASES AND MAINTAIN ALL OF THE 3D DATA ARRAYS
C       THE ARRAYS LINE3, POWL3 ARE ALL JUST CALCULATED 
C       AND THEN WRITTEN TO DISK FOR LATER ANALYSIS
C       TO SAVE STORAGE AT THE COST OF ADDITIONAL/INEFFICIENT 
C       COMPUTATION THE CALCULATION OF THESE QUANTITIES 
C       HAS BEEN MOVED TO THE DMPOUT ROUTINE WHERE THE TIZ3 AND LIM5 
C       ARRAY STORAGE ARE REUSED TO ALLOW THE CALCULATIONS TO PROCEED  
C
C       DAVID ELDER , JAN 29 , 1990 
C
C
C       IF (JY.LE.NY3D) THEN                                                    
C         DO 985 IP = -MAXNPS, MAXNPS                                           
C           DO 970 IX = 1, NXS                                                  
C             DO 970 IZ = 0, NIZS                                               
C               PNZS(IZ+1,1,IX) = 1.0E-6 * SNGL(DDLIM3(IX,IY,IZ,IP))            
C 970       CONTINUE                                                            
C           CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NXS)                      
C           DO 975 IX = 1, NXS                                                  
C             DO 975 IZ = 0, NIZS                                               
C               POWL3(IX,IY,IZ,IP)=MAX(0.0, PRADIS(1,IZ+1,1,IX)*1.0E6)          
C               LINE3(IX,IY,IZ,IP)=MAX(0.0, PRADIS(3,IZ+1,1,IX)*1.0E6)          
C 975       CONTINUE                                                            
C           DO 980 IX = 1, NXS                                                  
C             IF (DDLIM3(IX,IY,0,IP).LE.0.0D0) THEN                             
C               POWL3(IX,IY,-1,IP) = 0.0                                        
C               LINE3(IX,IY,-1,IP) = 0.0                                        
C             ELSE                                                              
C               POWL3(IX,IY,-1,IP) = POWL3(IX,IY,0,IP) /                        
C    >            SNGL(DDLIM3(IX,IY,0,IP)) * SNGL(DDLIM3(IX,IY,-1,IP))          
C               LINE3(IX,IY,-1,IP) = LINE3(IX,IY,0,IP) /                        
C    >            SNGL(DDLIM3(IX,IY,0,IP)) * SNGL(DDLIM3(IX,IY,-1,IP))          
C             ENDIF                                                             
C 980       CONTINUE                                                            
C 985     CONTINUE                                                              
C       ENDIF                                                                   


  990 CONTINUE                                                                  


c
c     USE ADAS DATA
c

      elseif (cdatopt.eq.1) then 



C
      DO 1190 IY = -NYS,NYS
C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
        IF (IY.EQ.0) GOTO 1190                                                   
        JY = IABS (IY)                                                          
C                                                                               


        DO 1100 IX = 1, NXS   
          PTESA(IX) = CTEMBS(IX,IY)
          PNESA(IX) = CRNBS(IX,IY) * RIZB
          PNBS(IX) =  CRNBS(IX,IY)
c
c         Set hydrogen density to zero for now - not available in LIM
c          PNHS(IX) = pinaton(ik,ir)
          pnhs(ix) = 0.0

 1100   CONTINUE
        DO 1120 IX = 1, NXS
          DO 1110 IZ = 0, NIZS
            PNZSA(IX,IZ) = SNGL(DDLIMS(IX,IY,IZ))
 1110     CONTINUE
 1120   CONTINUE
C
C------ GET POWER LOSS FROM ADAS DATA FILES. LOAD TOTAL LINE RADIATION
C------ INTO LINES AND ADD RECOMBINATION AND BREMSSTRAHLUNG POWER TO
C------ GET TOTAL RADIATIVE LOSSES
C
        write(year,'(i2.2)') iyearz
        call xxuid(useridz)
c        YEAR = '89'
c        YEARDF = '89'
        ICLASS = 5
        MIZS = MIN(CION-1,NIZS)
        DO 1130 IZ = 0, MIZS
          CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NXS,PTESA,PNESA,
     +                PCOEF(1,IZ+1))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ+1))
          DO 1135 IX = 1, NXS
            LINES(IX,IY,IZ) = PCOEF(IX,IZ+1)*PNESA(IX)*PNZSA(IX,IZ)
            POWLS(IX,IY,IZ) = LINES(IX,IY,IZ)
c            write (6,'(a,3i5,3g16.8)') 'Debug DIV:',ir,ik,iz,
c     >              pcoef(ik,iz+1),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      DIV:',
c     >            lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
 1135     CONTINUE
 1130   CONTINUE
        ICLASS = 4
        MIZS = MIN(CION,NIZS)
        DO 1140 IZ = 1, MIZS
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NXS,PTESA,PNESA,
     +                PCOEF(1,IZ))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ))
          DO 1145 IX = 1, NXS
            POWLS(IX,IY,IZ) = POWLS(IX,IY,IZ)
     +                        + PCOEF(IX,IZ)*PNESA(IX)*PNZSA(IX,IZ)
 1145     CONTINUE
 1140   CONTINUE
C
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1)
C
        DO 1160 IX = 1, NXS
          IF (DDLIMS(IX,IY,0).LE.0.0) THEN
            POWLS(IX,IY,-1) = 0.0
            LINES(IX,IY,-1) = 0.0
          ELSE
            POWLS(IX,IY,-1) = POWLS(IX,IY,0) *
     +                        SNGL(DDLIMS(IX,IY,-1) / DDLIMS(IX,IY,0))
            LINES(IX,IY,-1) = LINES(IX,IY,0) *
     +                        SNGL(DDLIMS(IX,IY,-1) / DDLIMS(IX,IY,0))
c
c            write (6,'(a,3i5,3g16.8)') 'Debug POW:',ir,ik,iz,
c     >              pcoef(ik,iz),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      POW:',
c     >           lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
          ENDIF
 1160   CONTINUE


 1190 CONTINUE


      endif


      WRITE(6,*) '12:'
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE Z EFFECTIVE ETC OVER "NEAR" REGION... NOTE 107                        WRITE (datunit,'(1X,A,I7,4X,I7)') NAME,I1,I2

C-----------------------------------------------------------------------        
C                                                                               
C---- THREE QUANTITIES PER BIN:  NIE,   ZB.NBT(X),   ZEFF                       
C---- THE AVERAGE VALUES OVER THE "NEAR" REGION ARE SUMMED IN THE DTOTS         
C---- ARRAY BELOW, JUST OVER THE REGION 0:CXNEAR, -CYNEAR:CYNEAR.               
C---- SEE ALSO NOTES 139, 221, 225, 288, 293                                    
C                                                                               
      DEFACT = 0.0D0                                                            
      IF (TNEUT.GT.0.0) DEFACT = DBLE (2.0*GTOT1*YEFF*CSEF*TATIZ/TNEUT)         
c
c     Assign DEFACT to ABSFAC which is in the common include file comtor 
c     This is for compatibility with some DIVIMP code.
c
      ABSFAC = DEFACT
c
c     Calculation of scaling factor:
c
      call prb
      call prc(' CALCULATION OF "ABSOLUTE" FACTOR:')
      call prc(' FORMULA USED: ABSFAC ='//
     >         ' 2.0*GTOT1*YEFF*CSEF*TATIZ/TNEUT')
      call prr(' ABSFAC = ',real(absfac))
      call prr(' GTOT1  = ',gtot1)
      call prr(' YEFF   = ',yeff)
      call prr(' CSEF   = ',csef)
      call pri(' TATIZ  = ',nint(tatiz))
      call pri(' TNEUT  = ',nint(tneut))
c
C                                                                               
      DO 1020 IX = 1, NXS                                                       
        DO 1010 IY = -NYS, NYS                                                  
          JY = IABS(IY)                                                         
C                                                                               
C-------- SUM CLOUD DENSITIES OVER IONISATION STATES                            
C-------- CALCULATE TOTAL IMPURITY INFLUX, MULTIPLY BY SPUTTERING               
C-------- ENHANCEMENT FACTOR...                                                 
C                                                                               
          DSUM1 = 0.0D0                                                         
          DSUM2 = 0.0D0                                                         
          DO 1000 IZ = 1, NIZS                                                  
            DIZ   = DFLOAT (IZ)                                                 
            DSUM1 = DSUM1 +  DIZ  * DDLIMS(IX,IY,IZ)                            
            DSUM2 = DSUM2 +DIZ*DIZ* DDLIMS(IX,IY,IZ)                            
 1000     CONTINUE                                                              
          DSUM1 = DSUM1 * DEFACT                                                
          DSUM2 = DSUM2 * DEFACT                                                
C                                                                               
C-------- CALCULATE AND STORE THE 3 ZEFFS RELATED QUANTITIES.                   
C-------- PRIOR TO NOTE 139, ZEFF WAS CALCULATED FROM THE FORMULA :-            
C-------- ZEFFS(IX,IY,3) = (RIZB * MAX (0.0,ZEFFS(IX,IY,2)) + SUM2) /           
C--------                  (RIZB * CRNBS(IX,IY))                                
C                                                                               
          ZEFFS(IX,IY,1) = SNGL (DSUM1)                                         
          ZEFFS(IX,IY,2) = RIZB * CRNBS(IX,IY) - ZEFFS(IX,IY,1)                 
          IF (ZEFFS(IX,IY,2).GT.0.0) THEN                                       
            ZEFFS(IX,IY,3) = (RIZB * ZEFFS(IX,IY,2) + SNGL(DSUM2)) /            
     >                       (RIZB * CRNBS(IX,IY))                              
          ELSE                                                                  
            ZEFFS(IX,IY,3) = SNGL (DSUM2/DSUM1)                                 
          ENDIF                                                                 
C                                                                               
C-------- NOTE 297.  DUPLICATE SET OF ZEFFS RESULTS BASED ON LT NOT LP+         
C-------- NOTE 299.  JUST USE SIN(THETAB) FACTOR, REMOVE LP+/LP- FACTOR         
C                                                                               
          ZEFFS(IX,IY,4) = SNGL (DSUM1) * CSINTB                                
          ZEFFS(IX,IY,5) = RIZB * CRNBS(IX,IY) - ZEFFS(IX,IY,4)                 
          IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                       
            ZEFFS(IX,IY,6) = (RIZB * ZEFFS(IX,IY,5) + SNGL(DSUM2) *             
     >        CSINTB) / (RIZB * CRNBS(IX,IY))                                   
          ELSE                                                                  
            ZEFFS(IX,IY,6) = SNGL (DSUM2/DSUM1)                                 
          ENDIF                                                                 
C         WRITE (6,'('' LIM3: IX,IY,NB'',2I5,G11.4,'' ==>'',3G11.4)')           
C    >      IX,IY,CRNBS(IX,IY),(ZEFFS(IX,IY,J), J=1,3)                          
 1010   CONTINUE                                                                
 1020 CONTINUE                                                                  
      WRITE(6,*) '13:'
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE TOTALS OVER ENTIRE PLASMA, ETC ....                             
C-----------------------------------------------------------------------        
C                                                                               
      CALL DZERO (DTOTS, 20)                                                    
      DO 4050 IY = -NYS, NYS                                                    
        IF (IY.EQ.0) GOTO 4050                                                  
        JY = IABS (IY)                                                          
        DO 4040 IX = 1, NXS                                                     
          DACT = XWIDS(IX) * XCYLS(IX) * YWIDS(JY)                              
          DTOTS(1) = DTOTS(1) + DACT                                            
          DTOTS(2) = DTOTS(2) + DACT * DBLE(CRNBS(IX,IY))                       
          DO 4030 IZ = 0, NIZS                                                  
            DTOTS(3) = DTOTS(3) + DACT * DDLIMS(IX,IY,IZ)                       
            DTOTS(4) = DTOTS(4) + DACT * DBLE(POWLS(IX,IY,IZ))                  
            DTOTS(5) = DTOTS(5) + DACT * DBLE(LINES(IX,IY,IZ))                  
            IF (XS(IX).LE.0.0) THEN                                             
              DTOTS(6) = DTOTS(6) + DACT                                        
              DTOTS(7) = DTOTS(7) + DACT * DBLE(POWLS(IX,IY,IZ))                
              DTOTS(8) = DTOTS(8) + DACT * DBLE(LINES(IX,IY,IZ))                
            ENDIF                                                               
 4030     CONTINUE                                                              
          IF (XS(IX).GT.0.0.AND.XS(IX).LE.CXNEAR.AND.YS(JY).LE.CYNEAR)          
     >    THEN                                                                  
            DTOTS(9) = DTOTS(9) + DACT                                          
            DTOTS(10)= DTOTS(10)+ DACT * DBLE(ZEFFS(IX,IY,4))                   
            DTOTS(11)= DTOTS(11)+ DACT * DBLE(ZEFFS(IX,IY,5))                   
            DTOTS(12)= DTOTS(12)+ DACT * DBLE(ZEFFS(IX,IY,6))                   
          ENDIF                                                                 
 4040   CONTINUE                                                                
 4050 CONTINUE                                                                  
      DTOTS(10) = DTOTS(10) / DTOTS(9)                                          
      DTOTS(11) = DTOTS(11) / DTOTS(9)                                          
      DTOTS(12) = DTOTS(12) / DTOTS(9)                                          
      DTOTS(13) = DTOTS(2) / DTOTS(1)                                           
      DTOTS(14) = DTOTS(4) / (DTOTS(13) * DTOTS(3))                             
      DO 4090 J = 1, 20                                                         
        STOTS(J) = SNGL(DTOTS(J))                                               
 4090 CONTINUE                                                                  
      WRITE(6,*) '14:'
C                                                                               
C-----------------------------------------------------------------------        
C                     PRINT CLOSING MESSAGE                                     
C-----------------------------------------------------------------------        
C                                                                               
      IF (NIZS.GT.0)                                                            
     >     CALL MONPRI (QTIM,FACTA(1),VFLUID,NIZS,SDTZS,SDYZS,                  
     >           STOTS,DOUTS,RIONS,CTBIN,CRMI)                                  
      CALL PRB                                                                  
      CALL PRI ('NUMBER OF NEUTRALS FOLLOWED   ',NINT(TNEUT))                   
      CALL PRI ('NUMBER OF IONS FOLLOWED       ',NINT(TATIZ))                   
      WRITE (6,'('' NUMBER OF NEUTRALS FOLLOWED   '',G11.4)') TNEUT             
      WRITE (6,'('' NUMBER OF IONS FOLLOWED       '',G11.4)') TATIZ             
C                                                                               
C---- FORMATS ...                                                               
C                                                                               
 9002 FORMAT(1X,I5,F9.1,12X,I4,4X,F10.6,10X,F10.5)                              
c slmod
 9003 FORMAT(1X,I5,F9.1,I4,I4,I4,I4,2F10.6,2F10.5,1P,G11.3,0P,F9.4,            
     >  1P,G10.3,0P,F7.4,3I3,1X,A,:,I3,I4,F8.2)                                 
c
c 9003 FORMAT(1X,I5,F9.1,I6,I6,I4,I4,2F10.6,2F10.5,1P,G11.3,0P,F8.2,            
c     >  1P,G10.3,0P,F7.4,3I3,1X,A,:,I3,I4,F8.2)                                 
c slmod end
 9004 FORMAT(//1X,'LIM DEBUG: DIAGNOSTICS TO BE PRINTED EVERY',I6,              
     >  ' TIMESTEPS  (DELTA T =',G10.3,' SECONDS).',//)                         
 9005 FORMAT(1X,'--ION-----TIME---IQX---IQY--IX--IY-----X-------ALPHA',         
     >  '--------Y---------P-----DRIFT-VEL----TEMP-PARA-DIFF',
     >  '--FRACT-IP-IT-IS',            
     >  12('-'))                                                                
 9006 FORMAT(//1X,'LIM DEBUG: TRACK DIAGNOSTICS TO BE RECORDED FOR FIRST
     >',I6,' PARTICLES')              
 9010 FORMAT(/5X,'                                        IONS SURVIV',         
     >                                          'ING     ',                     
     >       /5X,'SPLITTING    IONS CROSSING INWARDS  AFTER CROSSING ',         
     >                                          'OUTWARDS',                     
     >       /5X,'  PLANE       WEIGHTED  UNWEIGHTED    WEIGHTED  UNW',         
     >                                          'EIGHTED ')                     
 9011 FORMAT(1X,I3,'  X =',F6.3,'M',I9,I11,I13,I11)                             
 9012 FORMAT(/1X,A,A,A)                                                         
 9013 FORMAT(/1X,A,I7,A)
 9020 FORMAT(' IZ',7F8.3,/,(3X,7F8.3))                                          
 9021 FORMAT(I3,1P,1X,7G8.1,/,(4X,7G8.1))                                       
 9022 FORMAT(I3,7F9.1,/,(3X,7F9.1))                                             
c
c slmod begin
      WRITE(0,*) 'Done  LIM3'

      TITL2 = TITLE
      MIZS  = NIZS
      MIMPS = NIMPS
      MATIZ = NATIZ

      CALL OutputDiag
      CALL GetProfiles
c slmod end
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)            
      use error_handling
      IMPLICIT none                                                    
      REAL    OLDALP,ALPHA,OLDY,Y,XM,YM                                         
      INTEGER CIAB,IQX,IX,IOY,IOD                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  HIT:  WE KNOW ROUGHLY WHERE AN ION HAS HIT THE LIMITER.  THIS    *        
C  *  ROUTINE BACKTRACKS ALONG THE FINAL TRAJECTORY TO DETERMINE MORE  *        
C  *  PRECISELY WHERE THE HIT OCCURED.  IT RETURNS THE INDICES IQX,    *        
C  *  IX, IOY AND IOD INDICATING THE POSITION WHERE THE ION STRUCK     *        
C  *  THE SURFACE.  THE ROUTINE RETURNS (XM,YM) THE IMPACT POINT, BUT  *        
C  *  SHOULD HAVE VERY LITTLE EFFECT EXCEPT WHERE THE CROSS FIELD      *        
C  *  DIFFUSION STEPS ARE LARGE  (IE, BIG TIMESTEPS).  WRITTEN IN      *        
C  *  RESPONSE TO NOTE 292.                                            *        
C  *                                                                   *        
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comnet'                                                          
C     INCLUDE (COMNET)                                                          
      REAL    XT,XB,YT,YB,EDGE1,EDGE2,DIST1,DIST2                               
      INTEGER K,IPOS                                                            
      LOGICAL THERE                                                             
c
      real :: xt_org,xb_org
C                                                                               
      XT = OLDALP                                                               
      XB = ALPHA                                                                
c
      xt_org = xt
      xb_org = xb
c
      YT = OLDY                                                                 
      YB = Y                                                                    
      K  = 1                                                                    
      IF (DEBUGL) THEN                                                          
C       WRITE (6,9001) OLDALP,OLDY                                              
C       WRITE (6,9001) ALPHA,Y,IQX,QEDGES(IQX,1),QEDGES(IQX,2),K,.TRUE.         
      ENDIF                                                                     
C                                                                               
  100 CONTINUE                                                                  

      IF (XT.GT.0.0.AND.XB.GT.0.0) THEN 
         write(error_message_data,
     >               '(a,a,g18.10,a,g18.10,a,g18.10,a,g18.10)')
     >  'HIT CALLED WHEN OLDALP and ALPHA both greater than 0.0:',
     >               ' OLDALP =',oldalp,
     >               ' ALPHA =',alpha,
     >               ' OLDY =',oldy,
     >               ' Y =',y
         CALL errmsg('HIT:',error_message_data)


         xm = 0.0 - 1.0e-10
         xt = xm
         xb = xm
         YM = 0.5 * (YT + YB)                                                      
      else
         XM = 0.5 * (XT + XB)                                                      
         YM = 0.5 * (YT + YB)                                                      
      endif

      IQX = MAX (INT (XM*XSCALO), 1-NQXSO)                                      
C                                                                               
      IF (IQX.GT.0) THEN                                                        
        XT    = XM                                                              
        YT    = YM                                                              
        THERE =.FALSE.                                                          
C                                                                               
      ELSE                                                                      
        CALL EDGINT (XM,IQX,1,EDGE1,DIST1)                                      
        CALL EDGINT (XM,IQX,2,EDGE2,DIST2)                                      
        IF ((CIAB.EQ. 1 .AND. YM.LE.EDGE2)        .OR.                          
     >      (CIAB.EQ. 2 .AND. YM.GE.CTWOL-EDGE1)  .OR.                          
     >      (CIAB.EQ.-1 .AND. YM.GE.-EDGE1)       .OR.                          
     >      (CIAB.EQ.-2 .AND. YM.LE.EDGE2-CTWOL)) THEN                          
          XB    = XM                                                            
          YB    = YM                                                            
          THERE =.TRUE.                                                         
        ELSE                                                                    
          XT    = XM                                                            
          YT    = YM                                                            
          THERE =.FALSE.                                                        
        ENDIF                                                                   
      ENDIF                                                                     
C                                                                               
      K = K + 1                                                                 
      IF (K.LE.100) THEN                                                        
C       IF (DEBUGL) WRITE (6,9001) XM,YM,IQX,EDGE1,EDGE2,K,THERE                
        IF ((.NOT.THERE) .OR. (K.LE.5)) GOTO 100                                
      ENDIF                                                                     
C                                                                               
      if (xm.gt.0.0) then 
c
c        correct error by causing impact at limiter tip
c
         write(error_message_data,'(a,5(a,g18.10),a)')
     >        'Calculated XM greater than 0.0:',
     >        ' XM = ',xm,' XB =',xb,' XT = ',xt,
     >        ' XB_ORG =',xb_org,' XT_ORG = ',xt_org,
     >        ' XM=MIN(XT_ORG,XB_ORG) Assigned'

         CALL errmsg('HIT:',error_message_data)
         !xm = 0.0-1.0e-10
         xm = min(xb_org,xt_org)
         IQX = MAX (INT (XM*XSCALO), 1-NQXSO)                                      
         CALL EDGINT (XM,IQX,1,EDGE1,DIST1)                                      
         CALL EDGINT (XM,IQX,2,EDGE2,DIST2)                                      
      endif
c
      IX = IPOS (XM, XS, NXS-1)                                                 
      IF (CIAB.EQ.-1 .OR. CIAB.EQ.2) THEN                                       
        YM  =-EDGE1 - 1.E-10                                                    
        IOY = IPOS (-EDGE1, OYS, MAXOS-1)                                       
        IOD = IPOS (-DIST1, ODS, MAXOS-1)                                       
      ELSE                                                                      
        YM  = EDGE2 + 1.E-10                                                    
        IOY = IPOS ( EDGE2, OYS, MAXOS-1)                                       
        IOD = IPOS ( DIST2, ODS, MAXOS-1)                                       
      ENDIF                                                                     





      RETURN                                                                    
 9001 FORMAT(1X,'HIT: XM',F10.6,' YM',F10.5,:,                                  
     >  ' IQX',I6,' EDGES',2F10.6,I5,L2)                                        
      END                                                                       
