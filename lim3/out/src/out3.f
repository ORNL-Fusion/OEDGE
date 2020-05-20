      PROGRAM OUT3                                                              
      use mod_params
      use mod_comtor
      use error_handling
      use mod_dynam2
      use mod_dynam3
      use mod_comt2
      use mod_comnet
      use mod_comxyt
      use mod_coords
      use mod_rtheta
      use mod_pindata
      use mod_colours
      IMPLICIT  none
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *      THIS PROGRAM READS IN FROM AN EXTERNAL FILE                  *        
C  *    THE RESULTS OF A RUN OF LIM AND DISPLAYS A SELECTION OF PLOTS. *        
C  *                                                                   *        
C  *    C.M.FARRELL   JANUARY 1988                                     *        
C  *                  MASSIVELY UPDATED JULY 88  (LINE-OF-SIGHT ETC).  *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER   MAXQTS,MAXIB                                                    
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      PARAMETER (MAXQTS=(MAXIZS+1)*MAXNTS*2, MAXIB=MAXNXS*2*MAXNYS)             
c      INCLUDE   'dynam2'                                                        
C     INCLUDE   (DYNAM2)                                                        
c      INCLUDE   'dynam3'                                                        
C     INCLUDE   (DYNAM3)                                                        
c      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
c      INCLUDE   'comnet'                                                        
C     INCLUDE   (COMNET)                                                        
C                                                                               
c      INCLUDE   'rtheta'
c      include   'pindata' 
C
      COMMON /NSMOOTH/ NUMSMOOTH
      INTEGER NUMSMOOTH
C
      INTEGER   IPLANE,IX,IY,IZ,NIZS,NLS,IL,PIZS(MAXNLS),NRIGS,ILINE            
      INTEGER   IERR,MAXIZ,ISMOTH,JSMOTH,IALL,IPLOT,IVU,ISTATE,MAXIT            
      INTEGER   MLS,KSMOTH,J,IFOLD,IP,JZ,IRIG,IBAS3D,IVEW3D,NPTS,IQT            
      INTEGER   IPOS,IRX,IRY,IQX,LIMEDG,II,JY,KY,IT,IMODE,IFL,IEXP              
      INTEGER   IPRINT,IPAGE,NCONT,IXMIN,IXMAX,IYMIN,IYMAX,IPMIN,IPMAX       
      INTEGER   IMISC,JXMIN,JXMAX,JYMIN,JYMAX,JX,NYSLIM,JL,ITER,NITERS          
      INTEGER   IGZS(-2:MAXIZS+1),IGLS(MAXNLS),NIN,NQTS,NLOOPS,JMISC            
      INTEGER   IGSS(MAXNLS)
      INTEGER   ITEC,NAVS,IB,IGTS(MAXIB),IXCON,IDUM,NOS,IO,MULT,JO              
      INTEGER   IDUM2,CTIZ,TMPIZ
      INTEGER   MPTS,INLS,GRIND
      INTEGER   istart,istop
      integer   tftrg
C     
C     INDEX VARIABLES FOR R,THETA
C 
      INTEGER   IR
C
      REAL      PMIN,PMAX
      REAL      ystart,ystop,TOTDEPR
      REAL      GRIMIN,GRIMAX 
      REAL      XCON,OMAX,RULT                                                  
      REAL      XINTS(-MAXNYS:MAXNYS,-2:MAXIZS+1),PROJ3D,SUREDG(192,192)        
      REAL      YINTS(MAXNXS,-2:MAXIZS+1),SURFAS(192,192),YMAX,COSGC
      REAL      ZA02AS,TIME1,TIME,XMIN,XMAX,YMIN,ZSCALE,ZMAX,RIZB,SINGC         
      REAL      RV,XV,YV,BIGG,VMAX,VMIN,RIGS(MAXIZS+2),YLIM,YLIMIT,YV2          
      REAL      X,RATIO,POUTS(-MAXNPS:MAXNPS),CLEVLS(20),ZEFMAX,ZEFMIN          
      REAL      XYINTS(-MAXNPS:MAXNPS,-2:MAXIZS+1),TOTALS(4,-2:MAXIZS+1)        
      REAL      XJNTS(-MAXNYS:MAXNYS,MAXNTS+1),TMAX,RV1,XV1,YV1,RV2,XV2         
      REAL      YJNTS(MAXNXS,MAXNTS+1),XYJNTS(-MAXNPS:MAXNPS,MAXNTS+1)          
      REAL      TVALS(0:MAXQTS,-2:MAXIZS+1),QTS(0:MAXQTS),NEXT           
      REAL      QXFUNS(-MAXQXS:MAXQXS,2),AUX1(MAXNXS,1-MAXNYS:MAXNYS)           
      REAL      AUX2(MAXNXS,1-MAXNYS:MAXNYS),YSS(-MAXNYS:MAXNYS)                
      REAL      YOUTS1(-MAXNYS:MAXNYS),FRAC,YWIDSS(-MAXNYS:MAXNYS)              
      REAL      QXWIDS(-MAXQXS:MAXQXS),QTWIDS(0:MAXQTS),FACT,QTNEXT             
      REAL      FACTA(-1:MAXIZS),FACTB(-1:MAXIZS),FP,FS,FT,NP,NS,NT             
      REAL      SUM1(-2:MAXIZS),SUM2(-2:MAXIZS),XFUNS(MAXNXS,2),YY,WMIN         
      REAL      PLAMS(MAXNLS),AVS(0:250),BINA(MAXIB),BINB(MAXIB)                
      REAL      YFUNS(-MAXNYS:MAXNYS),V1,V2,TOTAV                               
      REAL      COORD1(192),COORD2(192)
      REAL      OYVOUT(MAXOS),OYVWID(MAXOS)
      real      tptracx(maxlen),tptracy(maxlen)
c slmod begin
      REAL      TMIN,TMPVOL,NUMSUM,VOLSUM,CPMIN,CPMAX,TSMAX
      INTEGER   I,L
      CHARACTER DUM*72
c slmod end

C
c       jdemod - deposition output variables
c
      character*13 :: value
      character*10000 :: plineout,lineout,pzoneout



c      
C     VARIABLES FOR R,THETA GRAPHS 
C
      REAL SINTS(MAXNSS,-2:MAXIZS+1)
      REAL TINTS(MAXNRS,-2:MAXIZS+1)
      REAL RINTS(MAXNAS,-2:MAXIZS+1)
      REAL SMIN,SMAX     
      CHARACTER*36 RTVIEW,SVIEW 
      CHARACTER*24 RLAB,THLAB,SLAB  
C
      CHARACTER*80   TITLE                                                      
      CHARACTER*24   XLAB,YLAB,PLAB,FLAB,TLAB,FFLAB,YLAB1,YLAB2,YLAB3           
      CHARACTER*24   YLAB4,YLAB5                                                      
      character*24   erolab
      CHARACTER*36   ZLABS(-2:MAXIZS+1),REF,ELABS(17),XVIEW,NVIEW,PLANE         
      CHARACTER*36   INTREF,OLABS(3),LABEL1
      CHARACTER*36   PLABS(MAXNLS),LAB3D1,TLABS(MAXNTS+1),ANLY 
      CHARACTER*72   LAB3D2
      CHARACTER*36   TABLE,NAMES(20),YVIEW,ZLABS1(-2:MAXIZS+1)                  
      CHARACTER*36   LAB3D3                                                     
      CHARACTER*72   JOB,GRAPH,DUMMY,SMOOTH                                     
      CHARACTER*12   INTEGP                                                     
      CHARACTER*7    PRINPS(-MAXNPS-1:MAXNPS)                                   
      CHARACTER*3    BREF                                                       
      CHARACTER*2    COMMAP                                                     
      CHARACTER*1    S(7)                                                       
      DOUBLE PRECISION DSUM(12),DDEPS,DSUM2(6)                               
      LOGICAL        SSS                                                        
      EXTERNAL       IPOS
c     jdemod - 1.0e50 is too large for R4 - use HI
c      DATA    BIGG / 1.0e50 / ,  NIN / 8 /               
      DATA    BIGG / HI / ,  NIN / 8 /               
c
c     Variables for LOS plots  
c
      character*80 graph2
      integer navg,iselect,iexpt,iaxis,iavg,ifact
      real optval  
      integer nplots
c
c     Erosion scaling variable
c     
      real :: ero_scale = 1.0
c
C
C     Initialization
C
      WRITE(0,*) 'Begin OUT3'
c
c     Initialize dynamically allocated storage
c
      call allocate_dynamic_storage
c
c     Initialize plot colours 
c
      call setup_col(16,3)      
c
c     Initialize string variables
c
      anly = ' '
      ref = ' '
      nview = ' '
      plane = ' '
c
      XLAB=' '
      YLAB=' '
      PLAB=' '
      FLAB=' '
      TLAB=' '
      FFLAB=' '
      YLAB1=' '
      YLAB2=' '
      YLAB3=' '
      YLAB4=' '
      YLAB5=' '
c
      JSMOTH   = 99
      NUMSMOOTH=3 
c
c     Initialize relevant PINDATA arrays to zero
c     - although EIRENE can not be invoked by LIM at 
c       the present time - some of the DIVIMP code that
c       has been imported has the ability to use PIN
c       data if it ever becomes available. As a result, 
c       a "pindata" common block is included which 
c       contains references to the variables used by
c       these DIVIMP routines.
c  
      call rzero(pinatom,maxnxs*maxnys)
      call rzero(pinion,maxnxs*maxnys)
      call rzero(pinalpha,maxnxs*maxnys)
C                                                                               
C-----------------------------------------------------------------------        
C     OPEN GHOST GRAPHICS                                                       
C     INITIALISATION: COLLECT VALUES AND SWITCH OFF UNDERFLOW INTERRUPTS        
C     CALL COLECT TO RETRIEVE LIM RESULTS                                       
C     SET FP = NO OF PRIMARIES LAUNCHED / NO. OF SECONDARIES LAUNCHED           
C         FT = TOTAL NEUTRALS LAUNCHED  / NO. OF SECONDARIES LAUNCHED           
C     LINE 10 IS CONTINUATION POINT FOR SUBSEQUENT ITERATIONS ...               
C-----------------------------------------------------------------------        
C                                                                               
      TIME1 = ZA02AS (1)                                                        
      CALL GPSTOP (100)                                                         
      CALL PAPER  (1)                                                           
c
c     jdemod - uncomment printer initialization
c
      CALL HRDLIN(1)
      CALL HRDCHR(1)
c
      CALL XUFLOW (0)                                                           
      IF (MAXY3D.LE.0) THEN                                                     
        WRITE (6,'('' OUT3 ERROR!!!  MAXY3D IS ONLY'',I5)') MAXY3D              
        STOP                                                                    
      ENDIF                                                                     
      REWIND (NIN)                                                              
C                                                                               
   10 CONTINUE                                                                  
      REWIND (5)                                                                
      IERR = 0                                                                  
      CALL COLECT (TITLE,NIZS,NIN,IERR,JOB,IMODE,PLAMS,PIZS,NLS,                
     >             FACTA,FACTB,ITER,NITERS)                                     
      NP = 0.0                                                                  
      NT = 0.0                                                                  
      IF (FACTA(-1).GT.0.0) NP = 1.0 / FACTA(-1)                                
      IF (FACTA(0) .GT.0.0) NT = 1.0 / FACTA(0)                                 
      NS = NT - NP                                                              
      IF (NS.GT.0.001) THEN                                                     
        FP = NP / NS                                                            
        FT = NT / NS                                                            
      ELSE                                                                      
        FP = 0.0                                                                
        FT = 0.0                                                                
      ENDIF                                                                     
      WRITE (6,'('' OUT3: NP,NS,NT,FP,FT='',5F10.3)') NP,NS,NT,FP,FT            
C                                                                               
C-----------------------------------------------------------------------        
C     1/8/88       ***** EXTRA SECTION FOR PLRPS *****                          
C-----------------------------------------------------------------------        
C                                                                               
C---- NOTE 210:  CREATE EXTRA DETAILS FOR SECONDARY NEUTRALS, BY                
C---- SUBTRACTING RESULTS FOR PRIMARIES FROM TOTALS LINES.                      
C---- TO INSERT THE DETAILS, ALL THE OTHERS HAVE TO BE SHIFTED UP               
C---- ONE PLACE - ARRAYS PLRPS, PLRP3, PIZS AND PLAMS ARE AFFECTED.             
C                                                                               
      IL = 1                                                                    
  410 CONTINUE                                                                  
      IF (PIZS(IL).EQ.-1) THEN                                                  
        IF (NLS.EQ.MAXNLS) THEN                                                 
          WRITE (6,'('' ERROR! MAXNLS TOO SMALL CREATING 2ND NEUTS'')')         
          STOP                                                                  
        ENDIF                                                                   
C                                                                              
        DO 430 IX = 1, NXS                                                      
C                                                                               
          DO 430 IY = -NYS, NYS                                                 
            DO 420 JL = NLS, IL, -1                                             
              PLRPS(IX,IY,JL+1) = PLRPS(IX,IY,JL)                               
  420       CONTINUE                                                            
            PLRPS(IX,IY,IL) = FT * PLRPS(IX,IY,IL+2) -                          
     >                        FP * PLRPS(IX,IY,IL+1)                            
  430     CONTINUE                                                              
C                                                                               
C       MOVE UP SPECTROSCOPIC TEMPERATURES
C       TO MATCH EVEN THOUGH THERE IS NO VALUE TO REPLACE THEM
C
        DO 436 JL = NLS,IL,-1
           DO 432 IX = 1,NXS
              SCTXS(IX,JL+1) = SCTXS(IX,JL)
  432      CONTINUE 
           DO 434 IY = -NYS,NYS
              SCTYS(IY,JL+1) = SCTYS(IY,JL)
  434      CONTINUE
  436   CONTINUE  

        DO 460 IX = 1,NXS 
          DO 450 IY = -NY3D, NY3D                                               
            DO 450 IP = -MAXNPS, MAXNPS                                         
              DO 440 JL = NLS, IL, -1                                           
                PLRP3(IX,IY,JL+1,IP) = PLRP3(IX,IY,JL,IP)                       
  440         CONTINUE                                                          
  450     CONTINUE                                                              
  460   CONTINUE                                                                
C
        DO 464 IX = 1,NXS 
          DO 464 IY = -NY3D, NY3D                                               
            DO 464 IP = -MAXNPS, MAXNPS                                         
              PLRP3(IX,IY,IL,IP) = FT * PLRP3(IX,IY,IL+2,IP) -                  
     >                             FP * PLRP3(IX,IY,IL+1,IP)                    
  464     CONTINUE                                                              
C                                                                               
C                                                                               
        DO 470 JL = NLS, IL, -1                                                 
          PIZS(JL+1) = PIZS(JL)                                                 
          PLAMS(JL+1) = PLAMS(JL)                                               
  470   CONTINUE                                                                
C                                                                               
        PIZS(IL) = -2                                                           
        NLS = NLS + 1                                                           
        IL = IL + 1                                                             
      ENDIF                                                                     
C                                                                               
      IL = IL + 1                                                               
      IF (IL.LT.NLS) GOTO 410                                                   
C                                                                               
      DO 480 IX = 1, NXS                                                        
        SAVES(IX,-2) = FT * SAVES(IX,0) - FP * SAVES(IX,-1)                     
  480 CONTINUE                                                                  
C
C-----------------------------------------------------------------------        
C           PRINT HEADING - DIFFERENT FOR SUBSEQUENT ITERATIONS                 
C-----------------------------------------------------------------------        
C                                                                               
      IF (ITER.EQ.1) THEN                                                       
        CALL PRB                                                                
        WRITE (DUMMY,'(''*'',60('' ''),''*'')')                                 
        CALL PRC                                                                
     >('**************************************************************')        
        CALL PRC (DUMMY)                                                        
      WRITE(7,'('' *'',18X,''RUN OF OUT VERSION '',A5,18X,''*'')')VERSON        
      WRITE(7,'('' *'',18X,24(''-''),18X,''*'')')                               
        CALL PRC (DUMMY)                                                        
        WRITE (7,'(1X,''* '',A58,'' *'')') TITLE(1:58)                          
        CALL PRC (DUMMY)                                                        
        WRITE (7,'(1X,''* '',A58,'' *'')') JOB(1:58)                            
        CALL PRC (DUMMY)                                                        
        CALL PRC                                                                
     >('**************************************************************')        
        CALL PRB                                                                
        IF     (IMODE.EQ.0) THEN                                                
         CALL PRC ('OPERATION MODE 0  COMBINED IMPULSE & STEADY STATE')         
        ELSEIF (IMODE.EQ.1) THEN                                                
         CALL PRC ('OPERATION MODE 1  IMPULSE')                                 
        ELSEIF (IMODE.EQ.2) THEN                                                
         CALL PRC ('OPERATION MODE 2  STEADY STATE')                            
        ENDIF                                                                   
C                                                                               
      ELSE                                                                      
        CALL PRB                                                                
        WRITE (DUMMY,'(''*'',60('' ''),''*'')')                                 
        CALL PRC                                                                
     >('**************************************************************')        
        CALL PRC (DUMMY)                                                        
        WRITE (7,'(1X,''* '',17X,A20,21X,'' *'')') TITLE(61:80)                 
        CALL PRC (DUMMY)                                                        
        CALL PRC                                                                
     >('**************************************************************')        
      ENDIF                                                                     
      CALL PRB                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C         READ FIRST LINES FROM DATAFILE, INCLUDING OBSERVATION POSN            
C  LINES BEGINNING '$' ARE TREATED AS COMMENTS AND ECHOED IN RDR ETC            
C-----------------------------------------------------------------------        
C                                                                               
c     Before reading in the files - initialize any OUT unstructured inputs
c
      call InitializeOUTUnstructuredInput
c
c     READ IN the plot input file
c
      CALL RDC (DUMMY,'TITLE LINE',IERR)                                        
C                                                                               
      CALL RDR (RV1  , .TRUE. ,0.0, .FALSE., 0.0 , 'VIEW: RV   ', IERR)         
      CALL RDR (XV1  , .FALSE.,0.0, .FALSE., 0.0 , 'VIEW: XV   ', IERR)         
      CALL RDR (YV1  , .FALSE.,0.0, .FALSE., 0.0 , 'VIEW: YV   ', IERR)         
      CALL RDR (YLIM , .FALSE.,0.0, .FALSE., 0.0 , 'VIEW: YLIMIT',IERR)         
      CALL PRB                                                                  
      NVIEW = ' '                                                               
      RV2   = 100.0                                                             
      COSGC = COS (GC)                                                          
      SINGC = SIN (GC)                                                          
      XV2   = (SC-SV) * COSGC - (TC-TV) * SINGC + RV2                           
      YV2   = (TC-TV) * COSGC + (SC-SV) * SINGC                                 
C                                                                               
C---- IGNORE SELECTION - MAKES GRAPHS CLEARER TO CUT OUT SOME                   
C---- INTERMEDIATE STATES.                                                      
C                                                                               
      DO 121 IZ = -2, MAXIZS+1                                                  
        IGZS(IZ) = 1                                                            
  121 CONTINUE                                                                  
      DO 122 IL = 1, MAXNLS                                                     
        IGLS(IL) = 1                                                            
  122 CONTINUE                                                                  
      DO 123 IT = 1, 100                                                        
        IGTS(IT) = 1                                                            
  123 CONTINUE                                                                  
C                                                                               
      CALL RDRAR (RIGS,NRIGS,MAXIZS+3,-2.0,REAL(NIZS),.TRUE.,                   
     >                                     'STATES TO BE IGNORED',IERR)         
      IF (NRIGS.GT.0) THEN                                                      
        CALL PRB                                                                
        CALL PRC ('IGNORE OPTION REQUESTED :-')                                 
        DO 60 IRIG = 1, NRIGS                                                   
          JZ = NINT(RIGS(IRIG))                                                 
          IF (JZ.LT.-2.OR.JZ.GT.NIZS) THEN                                      
            CALL PRI ('  WARNING!   CANNOT IGNORE IONISATION STATE', JZ)        
          ELSE                                                                  
            CALL PRI ('  PLOTS WILL BE SUPPRESSED FOR CHARGE STATE', JZ)        
            IGZS(JZ) = 0                                                        
            DO 59 IL = 1, NLS                                                   
              IF (PIZS(IL).EQ.JZ) IGLS(IL) = 0                                  
   59       CONTINUE                                                            
          ENDIF                                                                 
   60   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C---- 3D PLOTTING PARAMETERS ...                                                
C                                                                               
      CALL RDI (IVEW3D,.TRUE.,  0 ,.TRUE., 3 ,'3D VIEW FLAG' ,IERR)             
      CALL RDR (PROJ3D,.TRUE.,-90.,.TRUE.,90.,'3D PROJECTION',IERR)             
      CALL RDI (IBAS3D,.TRUE.,  0 ,.TRUE., 2 ,'3D BASE FLAG' ,IERR)             
      CALL RDI (LIMEDG,.TRUE.,  0 ,.TRUE., 2 ,'3D LIM EDGE FLAG' ,IERR)         
      CALL PRB                                                                  
      CALL PRR ('3D VIEWPOINT:  VIEWING ANGLE   ',90.0*REAL(IVEW3D))            
      CALL PRR ('               PROJECTION ANGLE',PROJ3D)                       
      IF (IBAS3D.EQ.0) THEN                                                     
        CALL PRC ('               BASE WILL NOT BE DRAWN')                      
      ELSE                                                                      
        CALL PRC ('               BASE WILL BE DRAWN')                          
      ENDIF                                                                     
      IF (LIMEDG.EQ.0) THEN                                                     
        CALL PRC ('               LIMITER EDGE WILL NOT BE DRAWN')              
      ELSE                                                                      
        CALL PRC ('               LIMITER EDGE WILL BE SUPERIMPOSED')           
      ENDIF                                                                     
C                                                                               
C---- PRINT OPTION FOR TABLES ETC                                               
C                                                                               
      CALL RDI (IPRINT,.TRUE., -1 ,.FALSE.,0 ,'PRINT OPTION  ' ,IERR)           
      IF     (IPRINT.EQ.-1) THEN                                                
       CALL PRC ('PRINT OPTION:  Z EFFECTIVE TABLE (COMPLETE)')                 
      ELSEIF (IPRINT.EQ.0) THEN                                                 
       CALL PRC ('PRINT OPTION:  OFF')                                          
      ELSE                                                                      
       CALL PRC ('PRINT OPTION:  NIE, ZB.NBTRUE, ZEFFS, DENSITIES, ETC')        
       CALL PRI ('               NUMBER OF PAGES OF EACH TO PRINT',             
     >   IPRINT)                                                                
      ENDIF                                                                     
C                                                                               
C---- MISCELLANEOUS PLOT OPTIONS, INCLUDING SMOOTHING ETC                       
C                                                                               
      CALL RDI (NOS  ,.TRUE., 0,.TRUE.,MAXOS,'NPTS FOR NET EROS ',IERR)         
      CALL RDR (OMAX ,.TRUE.,0.0,.TRUE.,CL  ,'Y RANGE FOR EDG,NE',IERR)         
      CALL RDI (IMISC,.TRUE.,0,.FALSE., 0   ,'MISC PLOTS (1)    ',IERR)         
      CALL RDI (JMISC,.TRUE.,0,.FALSE., 0   ,'MISC PLOTS (2)    ',IERR)         
C
C     Specify the range for integration of Total Deposition 
C
      CALL RDR(ystart,.false.,0.0,.FALSE.,0.0,'Min Y for Total Dep',
     >         IERR)         
      CALL RDR(ystop,.false.,0.0,.FALSE.,0.0,'Max Y for Total Dep',
     >         IERR)         
C
      CALL RDR (XCON ,.TRUE.,0.0,.FALSE.,0.0,'XCON FOR AVERAGING',IERR)         
      CALL RDI (ITEC ,.TRUE.,0,.TRUE.,  2   ,'SMOOTH TECHNIQUE  ',IERR)         
      IF (ITEC.EQ.2) THEN
         CALL RDI(NUMSMOOTH,.TRUE.,3,.FALSE.,0,
     >              'NUMBER OF PTS FOR SMOOTH2',IERR) 
      ENDIF 
      CALL RDRAR (AVS,NAVS,251,0.0,1.E9,.FALSE.,'AVERAGE WEIGHTS',IERR)         
      NAVS = NAVS - 1                                                           
      CALL PRR ('CONTOUR PLOT:  AVERAGE OVER ALL Y FOR X >', XCON)              
      IF (ITEC.EQ.0) THEN                                                       
        CALL PRC ('SMOOTHING TECHNIQUE: CUBIC SPLINE')                          
      ELSEIF (ITEC.EQ.1) THEN                                                   
        CALL PRC ('SMOOTHING TECHNIQUE: WEIGHTED AVERAGE WITH FUNCTION')        
        WRITE (7,'((15X,7(F7.3,:,'' :'')))') (AVS(IABS(J)),J=-NAVS,NAVS)        
      ENDIF                                                                     
C
C     READ IN VIEWING POSITION FOR R,THETA PSI PLOTS AS WELL AS
C     THE R,THETA BIN SPECIFICATIONS
C
      CALL RDR(RVPOS,.TRUE.,0.0,.FALSE.,0.0,   'R VIEW FOR PSI',IERR)
      CALL RDR(TVPOS,.TRUE.,0.0,.TRUE. ,360.0, 'THVIEW FOR PSI',IERR)
      TVPOS = DEGRAD* TVPOS
      CALL RDI(NSS,  .TRUE., 0 ,.TRUE.,MAXNSS,  '# PSI POINTS ',IERR)
      CALL RDI(NRDIV,.TRUE., 0 ,.TRUE.,MAXNRDIV,'# PSI R BINS ',IERR)
      CALL RDRAR(RS,NRS,MAXNRS-1,0,CA-CAW,.TRUE.,'R BIN UPBOUNDS',IERR)
      CALL RDRAR(TS,NAS,INT(MAXNAS/2)-1,0,180.0,.TRUE.,
     >             'TH BIN UPBOUNDS',IERR)
c
C-----------------------------------------------------------------------        
c
c     Finished reading input file - process inputs
c
C-----------------------------------------------------------------------        
c
c     Check for new specified value of absfac and assign it if specified
c     

      if (new_absfac.gt.0.0) then 
         absfac = new_absfac
         call errmsg('WARNING:VALUE OF ABSFAC OVER-RIDDEN BY INPUT:',
     >               absfac)
      endif

c      write(0,*) 'ABSFAC:',absfac,new_absfac
C
c
c     If a modified scaling for the net erosion plots has been speficied
c     then calculate it now. 
c
      if (erosion_scaling_opt.eq.0) then 
c
c        Normal scale
c
         ero_scale = 1.0
         erolab = 'M**-1'
c     
      elseif (erosion_scaling_opt.eq.1) then 
c
c        Scaling to mm/hr for Beryllium - 1848 kg/m3
c               1.22e26 Be in 1m2 by 1mm thick layer
c
         if (cion.eq.4) then 
            ero_scale = absfac * 3600.0 / 1.22e26
            erolab = 'mm Be/hr'
c
c        Scaling to mm/hr for Carbon - 2267 kg/m3
c               1.13e26 C in 1m2 by 1mm thick layer
c
         elseif (cion.eq.6) then 
            ero_scale = absfac * 3600.0 / 1.13e26
            erolab = 'mm C/hr'
c
c        Scaling for other impurities not implemented
c
         else
            write(0,*) 'EROSION SCALING NOT'//
     >                 ' IMPLEMENTED FOR THIS SPECIES'
            ero_scale = 1.0
            erolab = 'M**-1'
         endif

c
      endif 
c
c     jdemod - only modify the contents of nerods since that is the
c              only array for which the mm/hr scalings are appropriate.
c
c     Modify contents of NEROXS, NEROYS, NERODS
c
c      do ix = 1,maxnxs
c         do ii = 1,5
c            do j = 1,2
c              neroxs(ix,ii,j) = neroxs(ix,ii,j) * ero_scale
c            end do
c         end do
c      end do
c
c      do io = 1,maxos
c         do ii = 1,6
c            neroys(io,ii) = neroys(io,ii) * ero_scale
c         end do
c      end do
c
      do io = 1,maxos
         do ii = 1,5
            nerods(io,ii) = nerods(io,ii) * ero_scale
            do ip = -maxnps,maxnps
               nerods3(io,ip,ii) = nerods3(io,ip,ii) * ero_scale
            end do
         end do
      end do


c
C-----------------------------------------------------------------------        
C
C     CONVERT INPUT ANGULAR BINS FROM DEGREES TO RADIANS
C
      DO 30 IT = 1,NAS
        TS(IT) = DEGRAD * TS(IT)
 30   CONTINUE
C                                                                               
C-----------------------------------------------------------------------        
C        CALCULATE BINS AND INDEX ARRAYS                                        
C-----------------------------------------------------------------------        
C                                                                               
C---- R VALUES HAVE BEEN READ FOR 0 < R < RW  (EXCLUSIVE).  SET NEXT &         
C---- LAST R VALUE TO THE WALL POSITION.                                                        
C                                                                               
      NRS     = NRS + 1                                                         
      RS(NRS) = CA-CAW                                                              
C                                                                               
C---- THETA(TH)  VALUES HAVE BEEN READ FOR 0 < TH < 2PI (EXCLUSIVE).SET NEXT             
C---- TH VALUE TO 2PI.             
C                                                                               
      NAS     = NAS + 1                                                         
      TS(NAS) = PI
      DO 31 IT = NAS+1,2*NAS-1
         TS(IT) = 2*PI - TS(2*NAS-IT)
   31  CONTINUE
      NAS= 2*NAS
      TS(NAS) = 2*PI
C                                                                               
C---- SET MIDPOINTS ROUTS, TOUTS AND WIDTHS RWIDS,TWIDS                
C                                                                               
C
C---- ALSO CALCULATE CORRESPONDING X,Y COORDINATES FOR THE BIN CENTRES IN THE
C---- R,THETA SYSTEM
C
      RWIDM = CA-CAW
      DO 40 IR = 1, NRS                                                        
        IF (IR.EQ.1) THEN                                                       
          ROUTS(1) = 0.5 * RS(1)                                               
          RWIDS(1) = RS(1)
        ELSE                                                                    
          ROUTS(IR) = 0.5 * (RS(IR) + RS(IR-1))                                 
          RWIDS(IR) = RS(IR)- RS(IR-1)
        ENDIF                                                                   
        RWIDM = MIN(RWIDM,RWIDS(IR))                                          
        RTXCOORD(IR) = IPOS(CA - ROUTS(IR),XS,MAXNXS-1)
  40  CONTINUE                                                                  
C                                                                               
      TWIDM = 2*PI
      DO 50 IT = 1, NAS                                                        
        IF (IT.EQ.1) THEN                                                       
          TOUTS(1) = 0.5 * TS(1)                                             
          TWIDS(1) = TS(1)                                                      
        ELSE                                                                    
          TOUTS(IT) = 0.5 * (TS(IT) + TS(IT-1))                             
          TWIDS(IT) = TS(IT) - TS(IT-1)                                         
        ENDIF                                                                   
        TOUTS1(IT) = TOUTS(IT)
        TWIDS1(IT) = TWIDS(IT)
        TWIDM = MIN (TWIDM, TWIDS(IT))                                          
        RTYCOORD(IT) = IPOS(CL * TOUTS(IT) / PI,YS,2*MAXNYS-1)
  50  CONTINUE                                                                  
C
C     SET RTDONE = .FALSE. TO INDICATE THAT THE VIEWPOINT COORDINATE  
C     TRANSFORMATION HAS NOT BEEN DONE. THIS TRANSFORMATION IN THE RINTS
C     SUBROUTINE CONTAINS A NUMBER OF TRIGONOMETRIC AND ALGEBRAIC
C     FUNCTION CALLS. SINCE THE VIEWPOINT IS FIXED INITIALLY THIS CODE 
C     SHOULD NOT BE EXECUTED FOR EVERY SET OF PLOTS.
C
      RTDONE = .FALSE.
C
C                                                                               
C-----------------------------------------------------------------------        
C                SETUP POSITIONS FOR OUTPUT,  CHARACTER STRINGS ETC             
C-----------------------------------------------------------------------        
C                                                                               
C---- SET X POSITIONS OUTSIDE OF 1:NXS SO THAT THEY WILL NOT CAUSE              
C---- PLOTTING ERRORS ...                                                       
C                                                                               
      IF (NXS.LT.MAXNXS) THEN                                                   
        DO 70 IX = NXS+1, MAXNXS                                                
          XOUTS(IX) = 2.0 * XOUTS(NXS)                                          
   70   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C---- SET QX WIDTHS                                                             
C                                                                               
      DO 72 IQX = -MAXQXS, 0                                                    
        QXWIDS(IQX) = 1.0 / XSCALO                                              
   72 CONTINUE                                                                  
      DO 73 IQX = 1, MAXQXS                                                     
        QXWIDS(IQX) = 1.0 / XSCALI                                              
   73 CONTINUE                                                                  
C                                                                               
      XLAB   = '   X (M) '                                                      
      IF (CTHETB.EQ.90.0) THEN 
         YLAB1  = '   Y-TOROIDAL (M)'                                        
      ELSE 
         WRITE(YLAB1,'(A18,F6.2)') '   Y-POL (M)   TB=',CTHETB
      ENDIF
      YLAB2  = '   T (M) '                                                      
      YLAB3  = '   DIST ALONG LIM SURF'                                         
      YLAB4  = '   Y-CONT (M)'                                                  
      YLAB5  = '   Y-VERTICAL (M)' 
      PLAB   = '   P (M) '                                                      
      TLAB   = '   T (S) '                                                      
      FLAB   = 'M**-1'                                                          
      FFLAB  = ' '                                                              
      TABLE  = 'SYMBOL TABLE'                                                   
      SMOOTH = ' '                                                              
      ANLY   = ' '                                                              
C                                                                               
      YOUTS1(0) = 0.0                                                           
      DO 80 IY = 1, MAXNYS                                                      
        IF (IY.LE.NYS) THEN                                                     
          YOUTS1(IY) = YOUTS(IY)                                                
        ELSE                                                                    
          YOUTS1(IY) = 2.0 * YOUTS(NYS)                                         
        ENDIF                                                                   
        YOUTS1(-IY)=-YOUTS1(IY)                                                 
   80 CONTINUE                                                                  
C                                                                               
      DO 85 IQX = -MAXQXS, MAXQXS                                               
        IF (IQX.LT.1-NQXSO) QXS(IQX) = 2.0 * QXS(1-NQXSO)                       
        IF (IQX.GT.NQXSI-1) QXS(IQX) = 2.0 * QXS(NQXSI-1)                       
   85 CONTINUE                                                                  
C                                                                               
C---- FURTHEST POUT POSITIONS SHOULD REALLY BE +/- INFINITY, BUT PWIDS          
C---- HAVE BEEN SET TO NOMINAL VALUES FOR THESE OUTER BINS.                     
C                                                                               
      DO 90 IP = -MAXNPS, MAXNPS                                                
        POUTS(IP) = PS(IP) - 0.5 * PWIDS(IP)                                    
   90 CONTINUE                                                                  
      WRITE(6,*) 'PS:',(PS(IP),IP=-MAXNPS,MAXNPS)
      WRITE(6,*) 'PWIDS:',(PWIDS(IP),IP=-MAXNPS,MAXNPS)
      WRITE(6,*) 'POUTS:',(POUTS(IP),IP=-MAXNPS,MAXNPS) 
C                                                                               
C---- SET QTS TO SORTED TIME POINTS ARRAY.  INSERT EXTRA POINTS MIDWAY          
C---- BETWEEN TRUE TIMEPOINTS.  THESE CAN THEN BE USED IN CONJUNCTION           
C---- WITH FIDDLED QTWIDS ARRAY FOR CALCULATING AREAS IN ROUTINE DRAW.          
C                                                                               
      NQTS   = 0                                                                
      NLOOPS = 0                                                                
      QTS(0) = 0.0                                                              
      QTWIDS(0) = 0.0                                                           
   95 NLOOPS = NLOOPS + 1                                                       
      QTNEXT = 100.0                                                            
      DO 97 IZ = 0, NIZS                                                        
        DO 97 IT = 1, NTS                                                       
          NEXT = DWELFS(IT)*DWELTS(IZ)                                          
          IF (NEXT.GT.QTS(NQTS)) QTNEXT = MIN (QTNEXT, NEXT)                    
   97 CONTINUE                                                                  
      WRITE (6,*) ' NLOOPS,QTNEXT,NQTS',NLOOPS,QTNEXT,NQTS                      
      IF (QTNEXT.LT.100.0) THEN                                                 
        QTS(NQTS+1) = 0.5 * (QTNEXT + QTS(NQTS))                                
        QTS(NQTS+2) = QTNEXT                                                    
        QTWIDS(NQTS+1) = QTNEXT - QTS(NQTS)                                     
        QTWIDS(NQTS+2) = 0.0                                                    
        NQTS = NQTS + 2                                                         
      ENDIF                                                                     
      IF (NLOOPS.LT.(NIZS+1)*NTS) GOTO 95                                       
      WRITE (6,'('' QTS   '',/1P,(2X,7G11.4))') (QTS   (IQT),IQT=0,NQTS)        
      WRITE (6,'('' QTWIDS'',/1P,(2X,7G11.4))') (QTWIDS(IQT),IQT=0,NQTS)        
C                                                                               
      ELABS(1) = '    TOTAL DEPOSITION'                                         
      ELABS(2) = '    PRIMARY REMOVAL'                                          
      ELABS(3) = '    TOTAL REMOVAL'                                            
      ELABS(4) = '    NET EROSION'                                              
      ELABS(5) = '    NENNL      '                                              
      ELABS(17)= '    PRIM. CF REMOVAL'
      ELABS(6) = '    NIE    '                                                  
      ELABS(7) = '    ZB*NBT '                                                  
      ELABS(8) = '    ZEFF   '                                                  
      ELABS(9) = '    POW TOT'                                                  
      ELABS(10)= '    TB(X) <'                                                  
      ELABS(11)= '    TB(X) >'                                                  
      ELABS(12)= '    NB(X) <'                                                  
      ELABS(13)= '    NB(X) >'                                                  
      ELABS(14)= '    QS(X)  '                                                  
      ELABS(15)= '    EDGE   '                                                  
      ELABS(16)= '    BINS   '                                                  
C
      OLABS(1) = '    PRIMARY BG FLUX'
      OLABS(2) = '    CROSS FIELD BG FLUX'
      OLABS(3) = '    TOTAL BG FLUX'
C                                                                               
      PRINPS(-MAXNPS-1) = '-INFNTY'                                             
      DO 100 IP = -MAXNPS,MAXNPS-1                                              
        WRITE (PRINPS(IP),'(F7.4)') PS(IP)                                      
  100 CONTINUE                                                                  
      PRINPS(MAXNPS) = ' INFNTY'                                                
C                                                                               
C---- SETUP LABELS FOR SYMBOL TABLE ...                                         
C---- FIRST 4 CHARACTERS ARE USED TO IDENTIFY LINES IN GRAPHS...                
C---- FOR IONISATION STATES LIKE 1,2,3 ETC THEY ARE ALTERNATELY PUT             
C---- IN THE FIRST FEW CHARACTERS SO THAT THE POSSIBILITY OF                    
C---- DRAWING THEM AT THE SAME PLACE ON THE PLOTS IS REDUCED - THEY             
C---- MAY STILL BE COINCIDENT OF COURSE, BUT LESS LIKELY.                       
C                                                                               
      ZLABS1(-2)= '  S SECNEUT'                                                 
      ZLABS1(-1)= 'P   PRINEUT'                                                 
      ZLABS1(0) = ' T  TOTNEUT'                                                 
      IF (NIZS.GT.0) THEN                                                       
        DO 110 IZ = 1, NIZS                                                     
          IF     (3*(IZ/3).EQ.IZ) THEN                                          
            WRITE(ZLABS1(IZ), '(I4   ,A5,I2)') IZ,'IONIZ',IZ                    
          ELSEIF (2*(IZ/2).EQ.IZ) THEN                                          
            WRITE(ZLABS1(IZ), '(I3,1X,A5,I2)') IZ,'IONIZ',IZ                    
          ELSE                                                                  
            WRITE(ZLABS1(IZ), '(I2,2X,A5,I2)') IZ,'IONIZ',IZ                    
          ENDIF                                                                 
  110   CONTINUE                                                                
      ENDIF                                                                     
      ZLABS1(NIZS+1) = 'A   ALL IZS'                                            
C                                                                               
C---- SETUP LABELS FOR POSSIBLE PLRPS CASES                                     
C                                                                               
      DO 120 IL = 1, NLS                                                        
        IF     (PIZS(IL).EQ.-2) THEN                                            
          WRITE (PLABS(IL),'(A7,I4)') '  S SEC',INT(PLAMS(IL))                  
        ELSEIF (PIZS(IL).EQ.-1) THEN                                            
          WRITE (PLABS(IL),'(A7,I4)') 'P   PRI',INT(PLAMS(IL))                  
        ELSEIF (PIZS(IL).EQ.0) THEN                                             
          WRITE (PLABS(IL),'(A7,I4)') ' T  TOT',INT(PLAMS(IL))                  
        ELSE                                                                    
          IF     (3*(PIZS(IL)/3).EQ.PIZS(IL)) THEN                              
            WRITE (PLABS(IL),'(I4   ,I2,I5)') PIZS(IL),                         
     >                                  PIZS(IL),INT(PLAMS(IL))                 
          ELSEIF (2*(PIZS(IL)/2).EQ.PIZS(IL)) THEN                              
            WRITE (PLABS(IL),'(I3,1X,I2,I5)') PIZS(IL),                         
     >                                  PIZS(IL),INT(PLAMS(IL))                 
          ELSE                                                                  
            WRITE (PLABS(IL),'(I2,2X,I2,I5)') PIZS(IL),                         
     >                                  PIZS(IL),INT(PLAMS(IL))                 
          ENDIF                                                                 
        ENDIF                                                                   
  120 CONTINUE                                                                  
C                                                                               
C---- DEAL WITH NET EROSION PLOTS AGAINST Y : CHECK NO. OF PTS REQUIRED         
C                                                                               
      MULT = MAXOS / NOS                                                        
      IF (MULT*NOS.NE.MAXOS) THEN                                               
         call errmsg('OUT3: ERROR - MAXOS MUST BE DIVISIBLE BY NOS',nos)
c        WRITE (7,'('' OUT3: ERROR - MAXOS MUST BE DIVISIBLE BY NOS'')')         
        IERR = 1                                                                
      ELSEIF (MULT.GT.1) THEN                                                   
        RULT = REAL (MULT)                                                      
        DO 1124 IO = 1, NOS                                                     
          OYS(IO)    = OYS(MULT*IO)                                             
          OYWIDS(IO) = RULT * OYWIDS(IO)                                        
          OYOUTS(IO) = OYS(IO) - 0.5 * OYWIDS(IO)                               
          ODS(IO)    = ODS(MULT*IO)                                             
          ODWIDS(IO) = RULT * ODWIDS(IO)                                        
          ODOUTS(IO) = ODS(IO) - 0.5 * ODWIDS(IO)                               
          DO 1123 II = 1, 5                                                     
            NEROYS(IO,II) = NEROYS(MULT*(IO-1)+1,II) / RULT                     
            NERODS(IO,II) = NERODS(MULT*(IO-1)+1,II) / RULT                     
            do ip = -maxnps,maxnps
               nerods3(io,ip,ii) = nerods3(MULT*(IO-1)+1,ip,ii)/rult
            end do

            DO 1122 JO = MULT*(IO-1)+2, MULT*IO                                 
              NEROYS(IO,II) = NEROYS(IO,II) + NEROYS(JO,II) / RULT              
              NERODS(IO,II) = NERODS(IO,II) + NERODS(JO,II) / RULT              
              do ip = -maxnps,maxnps
                 nerods3(io,ip,ii) = nerods3(IO,ip,ii)
     >                              +nerods3(jo,ip,ii)/rult
              end do

 1122       CONTINUE                                                            
 1123     CONTINUE                                                              
 1124   CONTINUE                                                                
C       WRITE (6,9004) MAXOS,NOS,(IO,OYWIDS(IO),OYS(IO),OYOUTS(IO),             
C    >    ODWIDS(IO),ODS(IO),ODOUTS(IO),IO=1,MAXOS)                             
      ENDIF                                                                     

C                                                                               
      IF (IERR.NE.0) THEN                                                       
        CALL PRC ('OUT3: INPUT ERRORS FOUND - PROGRAM ABORTING')                
        STOP                                                                    
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C   PRINT Z EFFECTIVES DETAILS, ETC                                             
C-----------------------------------------------------------------------        
C                                                                               
      DO 150 II = 1, 3                                                          
        IF (IPRINT.EQ.0.OR.(IPRINT.EQ.-1.AND.II.NE.3)) GOTO 150                 
        IPAGE = 0                                                               
        ZMAX = 0.0                                                              
        DO 125 IY = -NYS, NYS                                                   
          DO 125 IX = 1, NXS                                                    
            ZMAX = MAX (ZMAX, ZEFFS(IX,IY,II))                                  
  125   CONTINUE                                                                
        ZSCALE = 1.0                                                            
        IF (ZMAX.GT.0.0) ZSCALE = 10.0 ** IEXP (0.0, ZMAX)                      
        DO 140 IY = -NYS, NYS, 24                                               
          IPAGE = IPAGE + 1                                                     
          IF (IPAGE.GT.IPRINT.AND.IPRINT.GT.0) GOTO 140                         
          KY = MIN (NYS, IY+23)                                                 
          WRITE (7,9100) ELABS(II+5)(5:11),ZSCALE                               
          WRITE (7,9101) (YOUTS1(JY),JY=IY+1,KY,2)                              
          WRITE (7,9102) ('-----',JY=IY,KY)                                     
          DO 130 IX = NXS, 1, -1                                                
            WRITE (7,9103) XOUTS(IX),                                           
     >        (ZEFFS(IX,JY,II)/ZSCALE,JY=IY,KY)                                 
  130     CONTINUE                                                              
  140   CONTINUE                                                                
  150 CONTINUE                                                                  
C                                                                               
C---- THIS LOOP GOES FOR IZ = -1,NIZS, SINCE THE -2 LINE DOESN'T EXIST.         
C                                                                               
      DO 250 IZ = -1, NIZS                                                      
        IPAGE = 0                                                               
        ZMAX = 0.0                                                              
        DO 220 IY = -NYS, NYS                                                   
          DO 220 IX = 1, NXS                                                    
            ZMAX = MAX (ZMAX, SDLIMS(IX,IY,IZ))                                 
  220   CONTINUE                                                                
        ZSCALE = 1.0                                                            
        IF (ZMAX.GT.0.0) ZSCALE = 10.0 ** IEXP (0.0, ZMAX)                      
        DO 240 IY = -NYS, NYS, 24                                               
          IPAGE = IPAGE + 1                                                     
          IF (IPAGE.GT.IPRINT) GOTO 240                                         
          KY = MIN (NYS, IY+23)                                                 
          WRITE (7,9100) ZLABS1(IZ)(5:11),ZSCALE                                
          WRITE (7,9101) (YOUTS1(JY),JY=IY+1,KY,2)                              
          WRITE (7,9102) ('-----',JY=IY,KY)                                     
          DO 230 IX = NXS, 1, -1                                                
            WRITE (7,9103) XOUTS(IX),                                           
     >        (SDLIMS(IX,JY,IZ)/ZSCALE,JY=IY,KY)                                
  230     CONTINUE                                                              
  240   CONTINUE                                                                
  250 CONTINUE                                                                  
C                                                                               
 9100 FORMAT('1',///1X,A16,2X,'SCALE FACTOR =',1P,G11.2)                        
 9101 FORMAT(/1X,'  X       Y',F8.3,11F10.3)                                    
 9102 FORMAT(1X,9('-'),24A5)                                                    
 9103 FORMAT(1X,F7.4,2X,24F5.2)                                                 
C                                                                               
C-----------------------------------------------------------------------        
C     PLOT MISCELLANEOUS X FUNCTIONS GRAPHS WITH TB, NB, EDGE, QS, ETC          
C-----------------------------------------------------------------------        
C                                                                               
      CALL PRB                                                                  
      IF (MOD(IMISC/100000000,10).EQ.1) THEN                                    
        CALL PRC('GRAPH OF ELECTRON TEMPERATURE TBE(X) WILL BE PLOTTED')                
        REF = 'PLASMA ELECTRON TEMPERATURE'                                              
        DO IX = 1, NXS                                                      
          XFUNS(IX,1) = CTEMBS(IX,-1)                                           
          XFUNS(IX,2) = CTEMBS(IX,1)                                            
        end do
        CALL LIM_DRAW (XOUTS,XWIDS,XFUNS,MAXNXS,NXS,ANLY,                           
     >    2,99,CAW,CA,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                              
     >    JOB,TITLE,XLAB,FFLAB,ELABS(10),REF,NVIEW,PLANE,TABLE,1,2)             
c
        CALL PRC('GRAPH OF ION TEMPERATURE TBI(X) WILL BE PLOTTED')                
        REF = 'PLASMA ION TEMPERATURE'                                              
        DO IX = 1, NXS                                                      
          XFUNS(IX,1) = CTEMBSI(IX,-1)                                           
          XFUNS(IX,2) = CTEMBSI(IX,1)                                            
        end do
        CALL LIM_DRAW (XOUTS,XWIDS,XFUNS,MAXNXS,NXS,ANLY,                           
     >    2,99,CAW,CA,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                              
     >    JOB,TITLE,XLAB,FFLAB,ELABS(10),REF,NVIEW,PLANE,TABLE,1,2)             
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC/10000000,10).EQ.1) THEN                                     
        CALL PRC ('GRAPH OF DENSITIES NB(X) WILL BE PLOTTED')                   
        REF = 'PLASMA DENSITY'                                                  
        DO 350 IX = 1, NXS                                                      
          XFUNS(IX,1) = CRNBS(IX,-1)                                            
          XFUNS(IX,2) = CRNBS(IX,1)                                             
  350   CONTINUE                                                                
        CALL LIM_DRAW (XOUTS,XWIDS,XFUNS,MAXNXS,NXS,ANLY,                           
     >    2,99,CAW,CA,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                              
     >    JOB,TITLE,XLAB,FFLAB,ELABS(12),REF,NVIEW,PLANE,TABLE,1,2)             
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC/1000000,10).EQ.1) THEN                                      
        CALL PRC ('GRAPH OF TIMESTEP FACTORS Q(X) WILL BE PLOTTED')             
        REF = 'TIMESTEP MULTIPLIERS'                                            
        CALL LIM_DRAW (QXS,QXWIDS,QS,2*MAXQXS+1,2*MAXQXS+1,ANLY,                    
     >    1,99,CAW,CA,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                              
     >    JOB,TITLE,XLAB,FFLAB,ELABS(14),REF,NVIEW,PLANE,TABLE,1,2)             
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC/100000,10).EQ.1.OR.MOD(IMISC/10000,10).EQ.1) THEN           
        DO 400 IQX = 1-NQXSO, 0                                                 
          QXFUNS(IQX-1,1) =-QEDGES(IQX,1)                                       
          QXFUNS(1-IQX,1) = QEDGES(IQX,2)                                       
          QXFUNS(IQX-1,2) = QXS(IQX)                                            
          QXFUNS(1-IQX,2) = QXS(IQX)                                            
  400   CONTINUE                                                                
        IF (MOD(IMISC/100000,10).EQ.1) THEN                                     
          CALL PRC ('GRAPH OF LIMITER EDGE Y(X) WILL BE PLOTTED')               
          REF = 'LIMITER SHAPE'                                                 
          QXFUNS(0,1) = 0.0                                                     
          QXFUNS(0,2) = 0.0                                                     
c          CALL LIM_GRTSET (TITLE,REF,NVIEW,PLANE,JOB,-OMAX,OMAX,                    
c     >     -OMAX/2.5,0.0,TABLE,YLAB1,XLAB,2,SMOOTH,1,ANLY,1)                    
          CALL LIM_GRTSET (TITLE,REF,NVIEW,PLANE,JOB,-OMAX,OMAX,                    
     >      caw,0.0,TABLE,YLAB1,XLAB,2,SMOOTH,1,ANLY,1)                    
          CALL LIM_GRTRAC (QXFUNS(-NQXSO,1),QXFUNS(-NQXSO,2),                       
     >      2*NQXSO+1,ELABS(15),'LINE')                                         
          CALL FRAME                                                            
        ENDIF                                                                   
        IF (MOD(IMISC/10000,10).EQ.1) THEN                                      
          CALL PRC ('GRAPH OF LIMITER TIP ONLY WILL BE PLOTTED')                
          REF = 'LIMITER TIP ONLY'                                              
          QXFUNS(0,1) = OMAX                                                    
          QXFUNS(0,2) = OMAX                                                    
c          CALL LIM_GRTSET (TITLE,REF,NVIEW,PLANE,JOB,-OMAX/2.5,OMAX/2.5,            
c     >     -OMAX/50.0,0.0,TABLE,YLAB1,XLAB,2,SMOOTH,1,ANLY,1)                   
          CALL LIM_GRTSET (TITLE,REF,NVIEW,PLANE,JOB,-OMAX/2.5,OMAX/2.5,            
     >      caw,0.0,TABLE,YLAB1,XLAB,2,SMOOTH,1,ANLY,1)                   
          CALL LIM_GRTRAC (QXFUNS(-NQXSO,1),QXFUNS(-NQXSO,2),                       
     >      2*NQXSO+1,ELABS(15),'POINT')                                        
          CALL FRAME                                                            
        ENDIF                                                                   
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC/1000,10).EQ.1) THEN                                         
        CALL PRC ('GRAPH OF NET EROSION AGAINST Y WILL BE PLOTTED')             
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        CALL DZERO (DSUM, 10)                                                   
        DO 530 IO = 1, NOS                                                      
         DO 530 II = 1, 5                                                       
          IF (NEROYS(IO,II).GT.0.0) THEN                                        
           DSUM(II)   = DSUM(II)   + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ENDIF                                                                 
  530   CONTINUE                                                                
C
C       ONLY CALCULATE INTEGRATION ACROSS FACE IF ystop > ystart
C
        IF (ystop.GT.ystart) THEN
          istart = IPOS(ystart,OYS,NOS)
          istop = IPOS(ystop,OYS,NOS)          
          TOTDEPR = 0.0         
          DO 535 II = istart,istop
             TOTDEPR = TOTDEPR + NEROYS(II,1)*OYWIDS(II)  
 535      CONTINUE
          TOTDEPR = TOTDEPR / (DSUM (1) + DSUM(6))
          WRITE(PLANE,536) ystart,ystop,TOTDEPR
 536      FORMAT('FRACTION:',F7.3,',',F7.3,'=',F8.4)  
        ENDIF
C
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)             
c
c        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c
c
        CALL LIM_DRAW (OYOUTS,OYWIDS,NEROYS,MAXOS,NOS,ANLY,                         
     >    5,99,-OMAX,OMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                        
     >    JOB,TITLE,YLAB1,FLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                 
        plane = ' ' 
        nview  = ' '

c
c     Write out the erosion data
c        
        write(6,'(a)') 'EROSION AS A FUNCTION OF Y:'
        write(6,'(4(2x,a20))') 'DISTANCE','TOTAL DEPOSITION',
     >            'TOTAL REMOVAL','NET EROSION' 
        do io = 1,nos
           write(6,'(4(4x,g18.8))') oyouts(io),neroys(io,1),
     >            neroys(io,3),neroys(io,4)
        end do



      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC/100,10).EQ.1) THEN                                            
        CALL PRC ('NET EROSION AGAINST DISTANCES ALONG SURFACE WILL BE P        
     >LOTTED')                                                                  
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        CALL DZERO (DSUM, 10)                                                   
        DO 630 IO = 1, NOS                                                      
         DO 630 II = 1, 5                                                       
          IF (NERODS(IO,II).GT.0.0) THEN                                        
           DSUM(II)   = DSUM(II)   + DBLE (NERODS(IO,II)*ODWIDS(IO))            
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NERODS(IO,II)*ODWIDS(IO))            
          ENDIF                                                                 
  630   CONTINUE                                                                
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)                
c
        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c
c       Change FLAB to EROLAB
c
        CALL LIM_DRAW (ODOUTS,ODWIDS,NERODS,MAXOS,NOS,ANLY,                         
     >    5,99,-1.5*OMAX,1.5*OMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                
     >    JOB,TITLE,YLAB3,EROLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                 
        nview = ' '
c
c     Write out the erosion data
c        
        write(6,'(a)') 'DEPOSITION ALONG LIMITER SURFACE:'
        write(6,'(6(2x,a20))') 'DISTANCE','WIDTH','TOTAL DEPOSITION',
     >        'PRIMARY REMOVAL','TOTAL REMOVAL','NET EROSION' 
        do io = 1,nos
           write(6,'(6(4x,g18.8))') odouts(io),odwids(io),-nerods(io,1),
     >            -nerods(io,2),-nerods(io,3),-nerods(io,4)

        end do


        write(48,'(a)') ' DEPOSITION ALONG LIMITER SURFACE:'
        write(48,'(12(2x,a20))') 'ABS_DIST_S1','DIST_S1',
     >        'TOTAL_DEPOSITION',
     >        'PRIMARY_REMOVAL','TOTAL_REMOVAL','NET_EROSION',
     >        'ABS_DIST_S2','DIST_S2','TOTAL_DEPOSITION',
     >       'PRIMARY_REMOVAL','TOTAL_REMOVAL','NET_EROSION'
        
        do io = 1,nos/2
           write(48,'(12(4x,g18.8))') abs(odouts(io)),odouts(io),
     >        -nerods(io,1),-nerods(io,2),-nerods(io,3),-nerods(io,4),
     >                           abs(odouts(nos-io+1)),odouts(nos-io+1),
     >          -nerods(nos-io+1,1),-nerods(nos-io+1,2),
     >          -nerods(nos-io+1,3),-nerods(nos-io+1,4)
        end do
        
c        write(6,'(a)') 'EROSION ALONG LIMITER SURFACE:'
c        write(6,'(4(2x,a20))') 'DISTANCE','TOTAL DEPOSITION',
c     >            'TOTAL REMOVAL','NET EROSION' 
c        do io = 1,nos
c           write(6,'(8(4x,g18.8))') odouts(io),
c     >            nerods(io,1),sum(nerods3(io,:,1)),
c     >            nerods(io,3),sum(nerods3(io,:,3)),
c     >            nerods(io,4),sum(nerods3(io,:,4))
c
c        end do

c
c       Write out the poloidally resolved data
c
        write(48,'(a)')
        write(48,'(a)') 'POLOIDALLY RESOLVED DEPOSITION ALONG'//
     >                 ' LIMITER SURFACE:'
        write(48,'(a)')
        write(48,'(a)') 'TOTAL DEPOSITION:'

        plineout = '    OUTS   ABS_OUTS     WIDS  '
        pzoneout = '      0       0          0    '
        do ip = -maxnps,maxnps
c           write(0,*) 'IP:',ip,pzone(ip)
           if (pzone(ip).ne.0) then 
              write(value,'(1x,g12.5)') pouts(ip)
              plineout = trim(plineout)//value
              write(value,'(3x,i7,3x)') pzone(ip)
              pzoneout = trim(pzoneout)//value
           endif
        end do 
              
c        write(0,'(a)') trim(pzoneout)
c        write(0,'(a)') trim(plineout)

        write(48,'(a)') trim(pzoneout)
        write(48,'(a)') trim(plineout)

c        write(0,*) 'NOS:',nos
        do io = 1,nos
           write(lineout,'(3(1x,g12.5))') odouts(io),abs(odouts(io)),
     >                                    odwids(io)
           do ip = -maxnps,maxnps
              if (pzone(ip).ne.0.0) then 
                 write(value,'(1x,g12.5)') -nerods3(io,ip,1)
                 lineout=trim(lineout)//value
              endif
           end do
           write(48,'(a)') trim(lineout)
        end do


        write(48,'(a)')
        write(48,'(a)') 'TOTAL REMOVAL:'
        write(48,'(a)') trim(pzoneout)
        write(48,'(a)') trim(plineout)

        do io = 1,nos
           write(lineout,'(3(1x,g12.5))') odouts(io),abs(odouts(io)),
     >                                    odwids(io)
           do ip = -maxnps,maxnps
              if (pzone(ip).ne.0.0) then 
                 write(value,'(1x,g12.5)') nerods3(io,ip,3)
                 lineout=trim(lineout)//value
              endif
           end do
           write(48,'(a)') trim(lineout)
        end do


        write(48,'(a)')
        write(48,'(a)') 'NET EROSION:'
        write(48,'(a)') trim(pzoneout)
        write(48,'(a)') trim(plineout)

        do io = 1,nos
           write(lineout,'(3(1x,g12.5))') odouts(io),abs(odouts(io)),
     >                                    odwids(io)
           do ip = -maxnps,maxnps
              if (pzone(ip).ne.0.0) then 
                 write(value,'(1x,g12.5)') nerods3(io,ip,4)
                 lineout=trim(lineout)//value
              endif
           end do
           write(48,'(a)') trim(lineout)
        end do

        write(48,'(a)')


      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC/10,10).EQ.1) THEN                                           
        CALL PRC ('GRAPH OF BIN BOUNDARIES WILL BE PLOTTED')                    
        REF = 'BIN BOUNDARIES'                                                  
        IB = 0                                                                  
        DO 720 IX = 1, NXS                                                      
          DO 710 IY = 1, NYS / 2                                                
            YY = MOD (ABS(YOUTS(IY)), CL)                                       
            IF (YY.GT.CHALFL) YY = CL - YY                                      
            BINA(IB+1) = YS(IY)                                                 
            BINB(IB+1) = 0.0                                                    
            BINA(IB+2) = YS(IY)                                                 
            IF (XS(IX).GT.0.0) THEN                                             
              BINB(IB+2) = XS(IX) * (YY*CONI + 1.0)                             
            ELSE                                                                
              BINB(IB+2) = XS(IX) * (YY*CONO + 1.0)                             
            ENDIF                                                               
            BINA(IB+3) = YS(IY) - YWIDS(IY)                                     
            BINB(IB+3) = BINB(IB+2)                                             
            BINA(IB+4) = YS(IY) - YWIDS(IY)                                     
            BINB(IB+4) = 0.0                                                    
            IB = IB + 4                                                         
  710     CONTINUE                                                              
  720   CONTINUE                                                                
        CALL LIM_DRAW (BINA,BINA,BINB,MAXIB,IB,ANLY,                                
     >    1,99,0.0,CL*1.001,CKO*CAW,CKI*CA,IGTS,ITEC,AVS,NAVS,                  
     >    JOB,TITLE,YLAB1,XLAB,ELABS(16),REF,NVIEW,PLANE,TABLE,1,2)             
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(IMISC,10).EQ.1) THEN                                              
        CALL PRC ('DENSITY ALONG INNER ZIG-ZAG WILL BE PLOTTED')                
        REF = 'DENSITY ALONG INNER ZIG-ZAG'                                     
        IF (CANAL.LT.CA)                                                        
     >    WRITE (ANLY,'(''ANALYTIC EXTENSION FOR X >'',F6.3)') CANAL            
        DO 800 IY = 1, NYS                                                      
          YWIDSS(IY)  = YWIDS(IY)                                               
          YWIDSS(-IY) = YWIDS(IY)                                               
          YFUNS(IY)   = SDLIMS(NXS,IY,NIZS)                                     
          YFUNS(-IY)  = SDLIMS(NXS,-IY,NIZS)                                    
  800   CONTINUE                                                                
        YWIDSS(0) = 0.0                                                         
        YFUNS(0)  = 0.5 * (YFUNS(-1) + YFUNS(1))                                
        CALL LIM_DRAW (YOUTS1,YWIDSS,YFUNS,2*MAXNYS+1,2*MAXNYS+1,ANLY,              
     >    1,0,-CL,CL,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                               
     >    JOB,TITLE,YLAB1,FLAB,ZLABS1(NIZS),REF,NVIEW,PLANE,TABLE,1,2)          
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MOD(JMISC,10).EQ.1) THEN                                              
        CALL PRC ('DENSITY PRE-ANALYTIC EXTENSION IS TO BE PLOTTED')            
        REF = 'DENSITY BEFORE ANALYTIC EXT'                                     
        ANLY = ' '                                                              
        DO 900 IZ = -2, NIZS+1                                                  
          ZLABS(IZ) = ZLABS1(IZ)                                                
  900   CONTINUE                                                                
        CALL LIM_DRAW (XOUTS,XWIDS,SAVES,MAXNXS,NXS,ANLY,                           
     >    NIZS+4,99,CAW,CA,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,                     
     >    JOB,TITLE,XLAB,FLAB,ZLABS(-2),REF,NVIEW,PLANE,TABLE,2,2)              
      ENDIF                                                                     
C-----------------------------------------------------------------------
C
C     PLOT CROSS-FIELD PRIMARY REMOVAL. THIS WOULD HAVE BEEN 
C     INCLUDED WITH THE OTHER EROSION PLOTS BUT IT DOES NOT
C     OCCUR VERY OFTEN AND THE AMOUNT OF WORK REQUIRED IN ARRANGING
C     LABELS AND MAKING IT MORE ELEGANT WAS FAR IN EXCESS OF THE 
C     THE TIME AVAILABLE.
C
      IF (MOD(JMISC/10,10).EQ.1) THEN
C
        MULT = MAXOS / NOS       
        RULT = REAL (MULT)                                                      
        DO 1126 IO = 1, NOS                                                     
            NEROYS(IO,6) = NEROYS(MULT*(IO-1)+1,6) / RULT                     
            DO 1125 JO = MULT*(IO-1)+2, MULT*IO                                 
              NEROYS(IO,6) = NEROYS(IO,6) + NEROYS(JO,6) / RULT              
 1125     CONTINUE                                                              
 1126   CONTINUE                        
C
        CALL PRC ('GRAPH OF CF PRIM REMOVAL AGAINST Y WILL BE PLOTTED')        
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        DSUM(11) = 0.0D0
        DSUM(12) = 0.0D0
        DO 1127 IO = 1, NOS                                                 
          IF (NEROYS(IO,6).GT.0.0) THEN                                        
           DSUM(11)   = DSUM(11)   + DBLE (NEROYS(IO,6))            
          ELSE                                                                  
           DSUM(12) = DSUM(12) + DBLE (NEROYS(IO,6))            
          ENDIF                                                                 
1127    CONTINUE                                                                
        WRITE (ELABS(17)(22:30),'(''='',F8.4)')  DSUM(11)+DSUM(12)              
c
c        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c

        CALL LIM_DRAW (OYOUTS,OYWIDS,NEROYS(1,6),MAXOS,NOS,ANLY,                   
     >    1,99,-OMAX,OMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                        
     >    JOB,TITLE,YLAB1,FLAB,ELABS(17),REF,NVIEW,PLANE,TABLE,1,6)            
        nview = ' '

      ENDIF
C-----------------------------------------------------------------------
C
C     PLOT THE DFLUX PARALLEL, CROSS FIELD AND TOTAL 
C     VERSUS Y.
C
      IF (MOD(JMISC/100,100).EQ.1) THEN
C
        MULT = MAXOS / NOS       
        RULT = REAL (MULT)                                                      
        DO 1135 II = 1,3
          DO 1135 IO = 1, NOS                                               
            CDFLUX(IO,II) = CDFLUX(MULT*(IO-1)+1,II) / RULT                     
            DO 1135 JO = MULT*(IO-1)+2, MULT*IO                                 
              CDFLUX(IO,II) = CDFLUX(IO,II) + CDFLUX(JO,II) / RULT              
 1135   CONTINUE                                                              
C
        CALL PRC ('GRAPHS OF D FLUX AGAINST Y WILL BE PLOTTED')        
        WRITE (REF,'(A,I3,A)') 'NET D FLUX  (',NOS,' PTS)'                     
        DO 1137 II = 1,3 
          DSUM2(II) = 0.0D0
          DSUM2(II+1) = 0.0D0
          DO 1138 IO = 1, NOS                                                 
            IF (CDFLUX(IO,II).GT.0.0) THEN                                   
              DSUM2(II)= DSUM2(II)   + DBLE (CDFLUX(IO,II)*OYWIDS(IO))      
            ELSE                                                              
              DSUM2(II+1)=DSUM2(II+1) +DBLE (CDFLUX(IO,II)*OYWIDS(IO))      
            ENDIF                                                             
1138      CONTINUE
          WRITE (OLABS(II)(26:36),'(''='',G9.4)') DSUM2(II)+DSUM2(II+1)        
1137    CONTINUE                                                                
        CALL LIM_DRAW (OYOUTS,OYWIDS,CDFLUX,MAXOS,NOS,ANLY,                   
     >    3,99,-OMAX,OMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                        
     >    JOB,TITLE,YLAB1,FLAB,OLABS,REF,NVIEW,PLANE,TABLE,1,6)            

      ENDIF
C-----------------------------------------------------------------------
C
C     Plot the Net Erosion versus Y as a function of Vertical Y - valid 
C     only for TFTR Limiter - edge option 7.
C
      IF (MOD(JMISC/1000,10).EQ.1) THEN                                         
        CALL PRC ('GRAPH OF NET EROSION AGAINST Y VERTICAL WILL BE'//
     >            ' PLOTTED')             
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        CALL DZERO (DSUM, 10)                                                   
        DO 1200 IO = 1, NOS                                                      
         DO 1200 II = 1, 5                                                       
          IF (NEROYS(IO,II).GT.0.0) THEN                                        
           DSUM(II)   = DSUM(II)   + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ENDIF                                                                 
 1200   CONTINUE                                                                
C
C       ONLY CALCULATE INTEGRATION ACROSS FACE IF ystop > ystart
C
        IF (ystop.GT.ystart) THEN
          istart = IPOS(ystart,OYS,NOS)
          istop = IPOS(ystop,OYS,NOS)          
          TOTDEPR = 0.0         
          DO 1210 II = istart,istop
             TOTDEPR = TOTDEPR + NEROYS(II,1)*OYWIDS(II)  
 1210      CONTINUE
          TOTDEPR = TOTDEPR / (DSUM (1) + DSUM(6))
          WRITE(PLANE,536) CA*(SIN(ystart/CA)),CA*(SIN(ystop/CA)),
     >                     TOTDEPR
        ENDIF
C
C       Convert the plotting coordinates to Y - vertical 
C       Include the effective widths as well - though they should not
C       be used for anything. Note that the widths do not scale linearly-
C       In general OYWIDS contains constant values - OYVWID - will not.
C
        do 1220 ii = 1,nos
          oyvout(ii) = ca*(sin(oyouts(ii)/ca))
1220    continue
        do 1230 ii = 2,nos-1
          oyvwid(ii) = abs ( oyvout(ii+1)-oyvout(ii-1) )/2.0  
1230    continue
        oyvwid(1) = abs(oyvout(2)-oyvout(1))
        oyvwid(nos) = abs(oyvout(nos)-oyvout(nos-1))        
C
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)             
c
c        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c
c       Change FLAB to EROLAB
c
        CALL LIM_DRAW (OYVOUT,OYVWID,NEROYS,MAXOS,NOS,ANLY,                         
     >    5,99,-OMAX,OMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                        
     >    JOB,TITLE,YLAB5,FLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                 
        plane = ' '
        nview = ' '
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF (MOD(JMISC/10000,10).EQ.1) THEN                                      
        CALL PRC ('GRAPH OF VELOCITY SVYBAR(X) WILL BE PLOTTED')             
        REF = 'AVERAGE Y VELOCITY'                                            
        write(14,*) nqxso,nqxsi
        do 1231 ix = -maxqxs,maxqxs
           write(14,*) ix,':',svybar(ix),svyacc(ix)
1231    continue         
        CALL LIM_DRAW (QXS,QXWIDS,SVYBAR,2*MAXQXS+1,2*MAXQXS+1,ANLY,                    
     >    1,99,CAW,CA,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                              
     >    JOB,TITLE,XLAB,' VELOCITY (M/S) ','  SVYBAR(X) ',REF,NVIEW,
     >    PLANE,TABLE,1,2)             
      ENDIF                                                                     
C-----------------------------------------------------------------------
C
C     Plot the Net Erosion versus Y as a function of Vertical Y - valid 
C     only for TFTR Limiter - edge option 7. - Smoothed
C
      IF (MOD(JMISC/100000,10).EQ.1) THEN                                         
        CALL PRC ('GRAPH OF NET EROSION AGAINST Y VERTICAL WILL BE'//
     >            ' PLOTTED')             
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        CALL DZERO (DSUM, 10)                                                   
        DO 1240 IO = 1, NOS                                                      
         DO 1240 II = 1, 5                                                       
          IF (NEROYS(IO,II).GT.0.0) THEN                                        
           DSUM(II)   = DSUM(II)   + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ENDIF                                                                 
 1240   CONTINUE                                                                
C
C       ONLY CALCULATE INTEGRATION ACROSS FACE IF ystop > ystart
C
        IF (ystop.GT.ystart) THEN
          istart = IPOS(ystart,OYS,NOS)
          istop = IPOS(ystop,OYS,NOS)          
          TOTDEPR = 0.0         
          DO 1250 II = istart,istop
             TOTDEPR = TOTDEPR + NEROYS(II,1)*OYWIDS(II)  
 1250      CONTINUE
          TOTDEPR = TOTDEPR / (DSUM (1) + DSUM(6))
          WRITE(PLANE,536) CA*(SIN(ystart/CA)),CA*(SIN(ystop/CA)),
     >                     TOTDEPR
        ENDIF
C
C       Convert the plotting coordinates to Y - vertical 
C       Include the effective widths as well - though they should not
C       be used for anything. Note that the widths do not scale linearly-
C       In general OYWIDS contains constant values - OYVWID - will not.
C
        do 1260 ii = 1,nos
          oyvout(ii) = ca*(sin(oyouts(ii)/ca))
1260    continue
        do 1270 ii = 2,nos-1
          oyvwid(ii) = abs ( oyvout(ii+1)-oyvout(ii-1) )/2.0  
1270    continue
        oyvwid(1) = abs(oyvout(2)-oyvout(1))
        oyvwid(nos) = abs(oyvout(nos)-oyvout(nos-1))        
C
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)             
c
c        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c

        CALL LIM_DRAW (OYVOUT,OYVWID,NEROYS,MAXOS,NOS,ANLY,                         
     >    5,0,-OMAX,OMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                        
     >    JOB,TITLE,YLAB5,FLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                 
        plane = ' '
        nview = ' '
      ENDIF                                                                     
C-----------------------------------------------------------------------
C
C     END OF MISC PLOTS
C
c-----------------------------------------------------------------------
c
c     Initialization 
c
      nplots = 0 

c
c    
c   This is a terrible plotting system which I hate to perpetuate ...
c   ... but here goes ... 
C                                                                               
C#######################################################################        
C    LOS plots
C    ===================                                                        
C    LOS 
C
C#######################################################################        
C                                                                               
      WRITE (7,9198)                                                            
 9198 FORMAT(' ',//2X,                                                          
     >  'REF  TITLE   PLOTTING  INTEGRATION PLOT SMOOTH MAXIZ',     
     >  ' PLANE YFLD ALLZ VU',/1X,71('-'))                                      
 950  CONTINUE                                                                  
c
c      CALL RDG (GRAPH,VMIN,VMAX,
c     >          GRIMIN,GRIMAX,IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,IALL,           
c     >          IVU,'LINE OF GRAPH DETAILS',IERR)                               
c
      CALL RDG_plot(GRAPH,IPLOT,'LOS PLOT SWITCH',IERR)                               
      BREF = GRAPH(1:3)                                                         
c
c      write(0,*) 'RDG_PLOT:',iplot,bref 
c
      IF (IERR.NE.0)          GOTO 9999                                           
      IF (BREF(1:3).NE.'LOS'.and.bref(1:3).ne.'000') GOTO 1000                                           
      IF (IPLOT.EQ.0)         GOTO  950                                           
c
c     Load secondary plot details
c
c
      call rdg_los(graph2,npts,navg,iselect,istate,
     >                   iexpt,iaxis,iavg,ifact,optval,ierr) 
c
c     call plotting routine
c
      nplots = nplots + 1
c
      call plot_los(iselect,istate,npts,navg,
     >                  iexpt,iaxis,iavg,ifact,optval,graph,
     >                  iplot,job,title,table,avs,navs,nplots,
     >                  outunito,nizs,ierr) 
c
      GOTO 950                                                                 
C     *********                                                                 
C                                                                               
C#######################################################################        
C    2 DIMENSIONAL PLOTS                                                        
C    ===================                                                        
C  SOMETIMES A PARAMETER WILL BE CHANGED - MARK IT WITH A '*' IN OUTPUT         
C  LIST FOR POSTERITY.                                                          
C#######################################################################        
C                                                                               
      WRITE (7,9201)                                                            
 9201 FORMAT(' ',//2X,                                                          
     >  'REF  TITLE   PLOTTING  INTEGRATION PLOT SMOOTH MAXIZ',     
     >  ' PLANE YFLD ALLZ VU',/1X,71('-'))                                      
 1000 CONTINUE                                                                  
      CALL RDG (GRAPH,VMIN,VMAX,
     >          GRIMIN,GRIMAX,IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,IALL,           
     >          IVU,'LINE OF GRAPH DETAILS',IERR)                               

      BREF = GRAPH(1:3)                                                         
c      write(0,*) '2D:',iplot,trim(bref)
      IF (IERR.NE.0)        GOTO 9999                                           
      IF (BREF(1:1).NE.'2'.and.bref(1:3).ne.'000') GOTO 1958                                           
      IF (IPLOT.EQ.0)       GOTO 1000                                           
c
      DO 1010 J = 1, 7                                                          
        S(J) = ' '                                                              
 1010 CONTINUE                                                                  
C                                                                               
C---- P PLOTS NOT AVAILABLE IF DPOL=0 - NO POINT TRYING TO PLOT THEM            
C---- PLOTS OF PLRPS SOMETIMES NOT AVAILABLE                                    
C                                                                               
      IF     (IPLOT.NE.0.AND.BREF(3:3).EQ.'T'.AND.IMODE.EQ.2) THEN              
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ELSEIF (IPLOT.NE.0.AND.BREF(3:3).EQ.'P'.AND.CDPOL.LE.0.0) THEN            
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ELSEIF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') THEN                               
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- SMOOTHING NOT ALLOWED SOMETIMES...                                        
C---- NOT ENOUGH POINTS TO SMOOTH P PLOTS  (ONLY 9)                             
C                                                                               
      IF (JSMOTH.NE.99.AND.(BREF(2:2).EQ.'Z'.OR.                                
     >    BREF(3:3).EQ.'P')) THEN                                               
        JSMOTH = 99                                                             
        S(2) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- CAN'T PLOT ABOVE MAX RECORDED IONISATION STATE!                           
C                                                                               
c
c     For 2AD plot the maxiz quantity represents the number of particle 
c     tracks to plot ... so perform a different test. 
c
       IF (BREF.EQ.'2AD') then
          write(6,*) '2ad:',maxiz,cstept
          if (MAXIZ.GT.CSTEPT) THEN
             maxiz = cstept
             s(3) = '*'
          endif
      ELSEIF (MAXIZ.GT.NIZS) THEN                                                   
        MAXIZ = NIZS                                                            
        S(3) = '*'                                                              
      ELSEIF (MAXIZ.LT.-2) THEN                                                 
        MAXIZ = -2                                                              
        S(3) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- 3D RESULTS NOT CALCULATED FOR DEPOSITION, ETC                             
C---- ALSO NOT IF DPOL=0.                                                       
C                                                                               
      IF (IPLANE.NE.99.AND.(BREF(2:2).EQ.'Z'.OR.BREF.EQ.'2DX'.OR.               
     >    BREF.EQ.'2NX'.OR.BREF.EQ.'2WY'.OR.BREF(3:3).EQ.'P'.OR.                
     >    BREF(2:2).EQ.'Y'.OR.BREF(2:2).EQ.'T'.OR.CDPOL.LE.0.0)) THEN           
        IPLANE = 99                                                             
        S(4) = '*'                                                              
      ENDIF                                                                     
      IF (IPLANE.EQ.99) THEN                                                    
        PLANE  = ' '                                                            
        COMMAP = ',P'                                                           
      ELSEIF (IPLANE.GE.-MAXNPS.AND.IPLANE.LE.MAXNPS) THEN                      
        WRITE (PLANE,'(4A)') 'PLOT FOR  ',PRINPS(IPLANE-1),                     
     >    ' < P <',PRINPS(IPLANE)                                               
        COMMAP = ' '                                                            
      ELSE                                                                      
        IPLANE = 99                                                             
        PLANE = ' '                                                             
        COMMAP = ',P'                                                           
        S(4) = '*'                                                              
      ENDIF                                                                     
      ANLY   = ' '                                                              
      IF (CANAL.LT.CA .AND. IPLANE.EQ.99 .AND. (BREF(2:2).EQ.'C'.OR.            
     >  BREF(2:2).EQ.'R'.OR.BREF(2:2).EQ.'L'.OR.BREF(2:2).EQ.'P'.OR.            
     >  BREF(2:2).EQ.'Z') .AND. (BREF(3:3).EQ.'X'.OR.BREF(3:3).EQ.'Y'))         
     >  WRITE (ANLY,'(''ANALYTIC EXTENSION FOR X >'',F6.3)') CANAL              
C                                                                               
C---- SEPARATE RESULTS FOR -Y, +Y NOT ALWAYS AVAILABLE                          
C---- FOLDING N/A FOR X GRAPHS (EXCEPT DEPS/NETE), P GRAPHS                     
C                                                                               
      IF (IFOLD.NE.0.AND.BREF(3:3).NE.'Y'.AND.BREF.NE.'2DX'.AND.                
     >  BREF.NE.'2NX') THEN                                                     
        IFOLD = 0                                                               
        S(5) = '*'                                                              
      ENDIF                                                                     
      IF (IFOLD.NE.0.AND.IFOLD.NE.1) THEN                                       
        IFOLD = 1                                                               
        S(5) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- "ALL IONISATIONS LINE" NOT ALWAYS AVAILABLE - CAN ONLY BE INCLUDED        
C---- IF WE ARE INTEGRATING OVER ALL THE AVAILABLE IONISATION STATES,           
C---- IE. NIZS=MAXIZ.  NOTE THAT THE IGNORE OPTION CAN OF COURSE BE USED        
C---- TO SUPPRESS THE ACTUAL PLOTTING OF INTERMEDIATE LINES IF REQUIRED.        
C                                                                               
      IF (IALL.EQ.1.AND. (BREF.EQ.'2NX'.OR.                    
     >  BREF(2:2).EQ.'Z'.OR.BREF(2:2).EQ.'P'.OR.BREF(2:2).EQ.'T'.OR.            
     >  BREF(2:2).EQ.'Y'.OR.MAXIZ.LT.NIZS)) THEN                                
        IALL = 0                                                                
        S(6) = '*'                                                              
      ENDIF                                                                     
      IF (IALL.NE.0.AND.IALL.NE.1) THEN                                         
        IALL = 0                                                                
        S(6) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- "VIEW POINT FLAG" - SET UP VIEWING POSITION ...                           
C---- SETUP YWIDSS ARRAY TO SAME AS YWIDS ONLY EXTENDING TO -Y REGION.          
C---- (IF VIEWPOINTS NOT USED AND WE'RE USING DDLIMS NOT DDLIM3, THEN           
C---- SET YWIDSS TO YWIDS/2 SO THAT "TOTAL AREA"                                
C---- CALCULATION IN DRAW INTEGRATED OVER -2L:2L BECOMES EFFECTIVELY            
C---- A VALUE FOR -L:L (SIMPLEST WAY TO DEAL WITH THIS I THINK!))               
C---- CHANGED THIS DEC 88 : FIX IFL=7 TO DO THIS                                
C                                                                               
c
c     jdemod - IVU has been given another meaning for the 2NX and 2DX plots
c              it now selects the side of the limiter to plot if non-zero
c
c      IF (IVU.NE.0 .AND. (BREF.EQ.'2DX'.OR.BREF.EQ.'2NX'.OR.                    
      IF (IVU.NE.0 .AND. (                    
     >  BREF(2:2).EQ.'Z'.OR.BREF(3:3).EQ.'P'.OR.BREF(3:3).EQ.'T'.OR.            
     >  BREF.EQ.'2WY'.OR.IFOLD.NE.0.OR.BREF(2:2).EQ.'Y'.OR.                     
     >  BREF(2:2).EQ.'T')) THEN                                                 
        IVU = 0                                                                 
        S(7) = '*'                                                              
      ENDIF                                                                     
C                                                                               
      XVIEW = ' '                                                               
      YVIEW = ' '                                                               
      YLAB = YLAB1                                                              
      DO 1011 IY = -MAXNYS, MAXNYS                                              
        YOUTS(IY) = YOUTS1(IY)                                                  
        IF (IY.EQ.0) THEN                                                       
          YWIDSS(IY) = 0.0                                                      
        ELSE                                                                    
          YWIDSS(IY) = YWIDS(IABS(IY))                                          
        ENDIF                                                                   
 1011 CONTINUE                                                                  
      IFL = 2                                                                   
c
      IF (IVU.EQ.0.AND.IPLANE.EQ.99) IFL = 7                                    
C                                                                               
      IF (IVU.EQ.1) THEN                                                        
        RV = RV1                                                                
        XV = XV1                                                                
        YV = YV1                                                                
        WRITE (XVIEW,9001) RV,XV,YV                                             
        WRITE (YVIEW,9001) RV,XV,YV                                             
      ELSEIF (IVU.GE.2) THEN                                                    
        RV = RV2                                                                
        XV = XV2                                                                
        YV = YV2                                                                
        WRITE (XVIEW,9001) RV,XV,YV                                             
        WRITE (YVIEW,9001) RV,XV,YV                                             
        YLAB = YLAB2                                                            
        IF (GRAPH(20:20).EQ.'Y') GRAPH(20:20) = 'T'                             
      ENDIF                                                                     
C                                                                               
      YLIMIT = MIN (RV/3.0, YLIM)                                               
      NYSLIM = MIN (IPOS(YLIMIT,YS,NYS-1), NYS/2)                               
      WRITE (6,9006) RV,XV,YV,YLIMIT                                            
C                                                                               
C---- ADJUST VMIN,VMAX IF REQUIRED                                              
C                                                                               
C---- FIND INDICES FOR INTEGRATION RANGE
C
      IF     (BREF(3:3).EQ.'X') THEN                                            
        VMIN = MAX (VMIN, CAW)                                                  
        VMAX = MIN (VMAX, CA)                                                   

C
C       INTEGRATE OVER Y      
C
        IF ( ((GRIMAX-GRIMIN).GT.(2.0*CL)) .OR.
     >      ((GRIMAX.EQ.0.0).AND.(GRIMIN.EQ.0.0)) ) THEN 
           GRIMAX = CL
           GRIMIN = -CL
           IYMIN = -NYS/2
           IYMAX = NYS/2
           INTREF = 'INTEGRATED OVER Y'
        ELSE          
           IF (GRIMIN.LT.0.0) THEN 
              GRIND = -1
           ELSE 
              GRIND = 1
           ENDIF  
           IYMIN = SIGN(IPOS(ABS(GRIMIN),YS,NYS-1),GRIND)
           IF (GRIMAX.LT.0.0) THEN 
              GRIND = -1
           ELSE 
              GRIND = 1
           ENDIF             
           IYMAX = SIGN(IPOS(ABS(GRIMAX),YS,NYS-1),GRIND)
           IF (IYMIN.GT.IYMAX) IYMIN=IYMAX
           WRITE(INTREF,'(''INT:Y='',
     >          F6.3,'' TO'',F6.3)') GRIMIN,GRIMAX
        ENDIF      
      ELSEIF (BREF(3:3).EQ.'Y') THEN                                            
        IF (IPLANE.EQ.99) THEN                                                  
          IF (VMIN.LT.0.0) VMIN = MAX (VMIN,-CL)                                
          IF (VMAX.GT.0.0) VMAX = MIN (VMAX, CL)                                
c
c         jdemod - aug/09 - remove limitation for "-"y plots for wall quantities - since +/-y may be assymmetric (?)
c     
c          IF (BREF(2:2).EQ.'W') VMIN = MAX (VMIN, 0.0)                          
c
        ELSE                                                                    
          IF (VMIN.LT.0.0) VMIN = MAX (VMIN,-YS(NY3D))                          
          IF (VMAX.GT.0.0) VMAX = MIN (VMAX, YS(NY3D))                          
        ENDIF                                                                   
C
C       INTEGRATE OVER X
C
        IF (GRIMIN.EQ.0.0.AND.GRIMAX.EQ.0.0) THEN 
           IXMIN = 1
           IXMAX = NXS
           INTREF = 'INTEGRATED OVER X'
        ELSE
           IXMIN = IPOS(GRIMIN,XS,NXS-1)
           IXMAX = IPOS(GRIMAX,XS,NXS-1)
           WRITE(INTREF,'(''INT:X='',
     >          F6.3,'' TO'',F6.3)') GRIMIN,GRIMAX
        ENDIF      
      ELSEIF (BREF(3:3).EQ.'P') THEN                                            
        VMIN = MAX (VMIN, POUTS(1-MAXNPS))                                      
        VMAX = MIN (VMAX, POUTS(MAXNPS-1))                                      
      ELSEIF (BREF(3:3).EQ.'T') THEN                                            
        VMIN = MAX (VMIN, 0.0)                                                  
        TMAX = DWELTS(0) * DWELFS(NTS)                                          
        DO 1015 IZ = 0, MAXIZ                                                   
          TMAX = MAX (TMAX, DWELTS(IZ)*DWELFS(NTS))                             
 1015   CONTINUE                                                                
        VMAX = MIN (VMAX, TMAX)                                                 
      ENDIF                                                                     
C                                                                               
C---- WRITE OUT (POSSIBLY ADJUSTED) LINE OF GRAPH DETAILS                       
C                                                                               
      WRITE (7,'(2X,A12,4F6.2,I4,A1,I6,5(A1,I4),A1)')                           
     >  GRAPH,VMIN,VMAX,GRIMIN,GRIMAX,IPLOT,                             
     >  S(1),JSMOTH,S(2),MAXIZ,S(3),IPLANE,S(4),IFOLD,S(5),IALL,S(6),           
     >  IVU,S(7)                                                                
c      WRITE (0,'(2X,A12,4F6.2,I4,A1,I6,5(A1,I4),A1)')                           
c     >  GRAPH,VMIN,VMAX,GRIMIN,GRIMAX,IPLOT,                             
c     >  S(1),JSMOTH,S(2),MAXIZ,S(3),IPLANE,S(4),IFOLD,S(5),IALL,S(6),           
c     >  IVU,S(7)                                                                
      ISMOTH = JSMOTH + 3                                                       
      IF (ISMOTH.EQ.3) ISMOTH = 1                                               
      IF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') GOTO 1000                              
      IF (VMIN.GE.VMAX)                  GOTO 1000                              
      IF (IPLOT.EQ.0)                    GOTO 1000                              
C                                                                               
C---- SET NUMBER OF LINES TO PLOT FOR PLRPS CASES                               
C---- AND SET SMOOTHING PARAMETER BASED ON THIS.                                
C---- SET SSS FLAG FOR NOTE 172 DISTANCE SQUARED FACTOR FOR XVIEW PLOTS         
C                                                                               
      SSS = .FALSE.                                                             
      IF (BREF(2:2).EQ.'P') THEN                                                
        SSS = .TRUE.                                                            
        MLS = NLS                                                               
        DO 1020 IL = 1, NLS                                                     
          PLABS(IL)(22:30) = '         '                                        
          IF (PIZS(IL).GT.MAXIZ) MLS = MLS - 1                                  
 1020   CONTINUE                                                                
        KSMOTH = MLS + 1                                                        
        DO 1030 IL = MLS, 1, -1                                                 
          IF (PIZS(IL).GE.JSMOTH) KSMOTH = IL                                   
          IF (PIZS(IL).LT.0 .AND. JSMOTH.EQ.0) KSMOTH = IL                      
 1030   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      DO 1040 IZ = -2, MAXIZ                                                    
        ZLABS(IZ) = ZLABS1(IZ)                                                  
 1040 CONTINUE                                                                  
      ZLABS(MAXIZ+1) = ZLABS1(NIZS+1)                                           
C                                                                               
C-----------------------------------------------------------------------        
C  CLOUDS                                                                       
C-----------------------------------------------------------------------        
C                                                                               
      IF     (BREF.EQ.'2CY') THEN                                               
C     ===========================                                               
        CALL RINTX (SDLIMS,SDLIM3,IPLANE,MAXIZ,XINTS,IFOLD,RV,XV,YV,IVU,        
     >              SSS,NYSLIM,FP,FT,1.0,IXMIN,IXMAX)                    
        CALL PROJEC (MAXIZ,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,CSINTB)        
        
        REF = 'IMPURITY ' // INTREF // COMMAP                                     
        CALL LIM_DRAW (YOUTS,YWIDSS,XINTS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,YLAB,FLAB,ZLABS(-2),REF,YVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2CX') THEN                                               
C     ===========================                                               
        CALL RINTY (SDLIMS,SDLIM3,IPLANE,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,        
     >              SSS,NYSLIM,FP,FT,1.0,IYMIN,IYMAX)                  
        REF = 'IMPURITY ' // INTREF // COMMAP                            
        CALL LIM_DRAW (XOUTS,XWIDS,YINTS,MAXNXS,NXS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,XLAB,FLAB,ZLABS(-2),REF,XVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
      ELSEIF (BREF.EQ.'2CP') THEN                                               
C     ===========================                                               
        CALL RINTXY (SDLIMS,SDLIM3,IPLANE,MAXIZ,XYINTS,IFOLD,RV,XV,YV,0,        
     >               SSS,NYSLIM,FP,FT)                                          
        REF = 'IMPURITY INTEGRATED OVER X,Y'                                    
        CALL LIM_DRAW (POUTS,PWIDS,XYINTS,2*MAXNPS+1,2*MAXNPS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,PLAB,FLAB,ZLABS(-2),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
      ELSEIF (BREF.EQ.'2CT') THEN                                               
C     ===========================                                               
        CALL TINTT (LIM5,IPLANE,MAXIZ,QTS,TVALS,IFOLD,RV,XV,YV,0,               
     >              SSS,NYSLIM,FP,FT,MAXQTS,NQTS)                               
        REF = 'IMPURITY INT''D OVER X,Y' // COMMAP                              
        CALL LIM_DRAW (QTS,QTWIDS,TVALS,MAXQTS+1,NQTS+1,ANLY,                       
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,TLAB,FFLAB,ZLABS(-2),REF,NVIEW,PLANE,TABLE,IPLOT,2)         
C                                                                               
C     ELSEIF (BREF.EQ.'2CQ') THEN                                               
C     ===========================                                               
C      EXTRA 5 PLOTS GIVING X DISTRIBUTION AT GIVEN SET OF Y VALUES             
C                                   JAN 89  NOTE 280                            
C      DO 1080 J = 0, 4                                                         
C       IY = IPOS (REAL(J), YS, NYS-1)                                          
C       DO 1070 IX = 1, NXS                                                     
C        DO 1060 IZ = -1, MAXIZ                                                 
C         YINTS(IX,IZ) = SDLIMS(IX,IY,IZ)                                       
C         YINTS(IX,MAXIZ+1) = YINTS(IX,MAXIZ+1) + YINTS(IX,IZ)                  
C1060    CONTINUE                                                               
C        YINTS(IX,-2) = FT * YINTS(IX,0) - FP * YINTS(IX,-1)                    
C1070   CONTINUE                                                                
C       WRITE (REF,'(''IMPURITY AT Y='',F6.2)') YS(IY)                          
C       CALL LIM_DRAW (XOUTS,XWIDS,YINTS,MAXNXS,NXS,ANLY,                           
C    >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
C    >    JOB,TITLE,XLAB,FLAB,ZLABS(-2),REF,XVIEW,PLANE,TABLE,IPLOT,2)          
C1080  CONTINUE                                                                 
C                                                                               
C-----------------------------------------------------------------------        
C  DEPOSITION.  SMOOTHING PARAM NEEDS CAUTION HERE - USE JSMOTH FOR             
C               PLOT WITHOUT ANY NEUTRALS,  NO SMOOTHING (SET ISMOTH=5)         
C               FOR NET EROSION GRAPH.                                          
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2DX') THEN                                               
C     ===========================                                               
         J = 1                                                                   
c
        IF (IFOLD.EQ.1) J = 3                                                   

 1090   IF (J.EQ.1) REF = 'DEPOSITION  Y < 0'                                   
        IF (J.EQ.2) REF = 'DEPOSITION  Y > 0'                                   
        IF (J.EQ.3) REF = 'DEPOSITION  Y<0 + Y>0'                               
        CALL LIM_DRAW (XOUTS,XWIDS,DEPS(1,1,J),MAXNXS,NXS,ANLY,                     
     >    MAXIZ+IALL,JSMOTH,VMIN,VMAX,0.0,BIGG,IGZS(1),ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,FLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,2)           
        IF (J.EQ.1) THEN                                                        
          J = 2                                                                 
          GOTO 1090                                                             
        ENDIF                                                                   
C                                                                               
C---- NET EROSION: SET TOTAL DEPOSITION IN NEROXS(,1)LOCATIONS (ALL <0)         
C---- AND NET EROSION IN NEROXS(,4) LOCATIONS                                   
C---- IPLOT USED SLIGHTLY DIFFERENTLY FOR NET EROSION GRAPHS:                   
C---- IPLOT=1 MEANS PLOT USUAL UNNORMALISED GRAPH                               
C---- IPLOT=2 MEANS PLOT NET EROSION CURVE ONLY (TO FILL WHOLE SCREEN)          
C---- ROUTINE DRAW IS CALLED WITH IFL=6 TO GET A STRAIGHT LINE                  
C---- ALONG F(X) = 0 TO HELP SPOT WHERE NET EROSION IS OCCURING...              
C                                                                               
      ELSEIF (BREF.EQ.'2NX') THEN                                               
C     ===========================                                               
        J = 1                                                                   
c
        IF (IFOLD.EQ.1) J = 3                                                   
 1100   IF (J.EQ.1) REF = 'NET EROSION  Y < 0'                                  
        IF (J.EQ.2) REF = 'NET EROSION  Y > 0'                                  
        IF (J.EQ.3) REF = 'NET EROSION  Y<0 + Y>0'                              
        CALL DZERO (DSUM, 10)                                                   
        DO 1120 IX = 1, NXS                                                     
         DO 1120 II = 1, 5                                                      
          IF (NEROXS(IX,II,J).GT.0.0) THEN                                      
           DSUM(II)   = DSUM(II)   + DBLE (NEROXS(IX,II,J)*XWIDS(IX))           
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NEROXS(IX,II,J)*XWIDS(IX))           
          ENDIF                                                                 
 1120   CONTINUE                                                                
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)                

c
c        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c
c       Change FLAB to EROLAB
c
        IF (IPLOT.EQ.1) THEN                                                    

          CALL LIM_DRAW (XOUTS,XWIDS,NEROXS(1,1,J),MAXNXS,NXS,ANLY,                 
     >      5,JSMOTH,VMIN,VMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                   
     >      JOB,TITLE,XLAB,FLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                
        ELSEIF (IPLOT.EQ.2) THEN                                                
          CALL LIM_DRAW (XOUTS,XWIDS,NEROXS(1,4,J),MAXNXS,NXS,ANLY,                 
     >      2,JSMOTH,VMIN,VMAX,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                   
     >      JOB,TITLE,XLAB,FLAB,ELABS(4),REF,NVIEW,PLANE,TABLE,1,6)             
        ENDIF                                                                   
        IF (J.EQ.1) THEN                                                        
          J = 2                                                                 
          GOTO 1100                                                             
        ENDIF                                                                   
c
        nview = ' '
c
c     Write out the erosion data
c        
        write(6,'(a)') 'EROSION AS A FUNCTION OF X:    Y<0    Y>0'
        write(6,'(7(2x,a20))') 'X','TOTAL DEPOSITION',
     >            'TOTAL REMOVAL','NET EROSION',
     >            'TOTAL DEPOSITION',
     >            'TOTAL REMOVAL','NET EROSION' 
        do io = 1,nxs
           write(6,'(8(4x,g18.8))') xs(io),neroxs(io,1,1),
     >            neroxs(io,3,1),neroxs(io,4,1),neroxs(io,1,2),
     >            neroxs(io,3,2),neroxs(io,4,2)
        end do

C                                                                               
C---- EROSION AND DEPOSITION ALONG THE Y-AXIS: NEROYS
cC                                                                               
C---- IPLOT USED SLIGHTLY DIFFERENTLY FOR NET EROSION GRAPHS:                   
C---- IPLOT=1 MEANS PLOT USUAL UNNORMALISED GRAPH                               
C---- IPLOT=2 MEANS PLOT NET EROSION CURVE ONLY (TO FILL WHOLE SCREEN)          
C---- ROUTINE DRAW IS CALLED WITH IFL=6 TO GET A STRAIGHT LINE                  
C---- ALONG F(X) = 0 TO HELP SPOT WHERE NET EROSION IS OCCURING...              
C                                                                               
      ELSEIF (BREF.EQ.'2NY') THEN                                               
C     ===========================                                               
c
        CALL PRC ('GRAPH OF NET EROSION AGAINST Y WILL BE PLOTTED')             
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        CALL DZERO (DSUM, 10)                                                   
c
        DO IO = 1, NOS                                                      
         DO  II = 1, 5                                                       
          IF (NEROYS(IO,II).GT.0.0) THEN                                        
           DSUM(II)   = DSUM(II)   + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NEROYS(IO,II)*OYWIDS(IO))            
          ENDIF                                                                 
         end do
        end do
C
C       ONLY CALCULATE INTEGRATION ACROSS FACE IF ystop > ystart
C
        IF (vmax.GT.vmin) THEN
          istart = IPOS(vmin,OYS,NOS)
          istop = IPOS(vmax,OYS,NOS)          
          TOTDEPR = 0.0         
          DO II = istart,istop
             TOTDEPR = TOTDEPR + NEROYS(II,1)*OYWIDS(II)  
          end do
          TOTDEPR = TOTDEPR / (DSUM (1) + DSUM(6))
          WRITE(PLANE,537) ystart,ystop,TOTDEPR
 537      FORMAT('FRACTION:',F7.3,',',F7.3,'=',F8.4)  
        ENDIF
C
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)             
c
c        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c
c
        CALL LIM_DRAW (OYOUTS,OYWIDS,NEROYS,MAXOS,NOS,ANLY,                         
     >    5,99,vmin,vmax,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                        
     >    JOB,TITLE,YLAB1,FLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                 
        plane = ' ' 
        nview  = ' '

c
c     Write out the erosion data
c        
        write(6,'(a)') 'EROSION AS A FUNCTION OF Y:'
        write(6,'(4(2x,a20))') 'DISTANCE','TOTAL DEPOSITION',
     >            'TOTAL REMOVAL','NET EROSION' 
        do io = 1,nos
           write(6,'(4(4x,g18.8))') oyouts(io),neroys(io,1),
     >            neroys(io,3),neroys(io,4)

        end do

C                                                                               
C---- EROSION AND DEPOSITION ALONG THE LIMITER SURFACE: NERODS
c
C---- IPLOT USED SLIGHTLY DIFFERENTLY FOR NET EROSION GRAPHS:                   
C---- IPLOT=1 MEANS PLOT USUAL UNNORMALISED GRAPH                               
C---- IPLOT=2 MEANS PLOT NET EROSION CURVE ONLY (TO FILL WHOLE SCREEN)          
C---- ROUTINE DRAW IS CALLED WITH IFL=6 TO GET A STRAIGHT LINE                  
C---- ALONG F(X) = 0 TO HELP SPOT WHERE NET EROSION IS OCCURING...              
C                                                                               
      ELSEIF (BREF.EQ.'2ND') THEN                                               
C     ===========================                                               

        CALL PRC ('NET EROSION AGAINST DISTANCES ALONG SURFACE WILL BE P        
     >LOTTED')                                                                  
        WRITE (REF,'(A,I3,A)') 'NET EROSION  (',NOS,' PTS)'                     
        CALL DZERO (DSUM, 10)                                                   

        DO IO = 1, NOS                                                      
         DO II = 1, 5                                                       
          IF (NERODS(IO,II).GT.0.0) THEN                                        
           DSUM(II)   = DSUM(II)   + DBLE (NERODS(IO,II)*ODWIDS(IO))            
          ELSE                                                                  
           DSUM(II+5) = DSUM(II+5) + DBLE (NERODS(IO,II)*ODWIDS(IO))            
          ENDIF                                                                 
         end do
        end do
c
        WRITE (ELABS(1)(22:30),'(''='',F8.4)')  DSUM(1)+DSUM(6)                 
        WRITE (ELABS(2)(22:30),'(''='',F8.4)')  DSUM(2)+DSUM(7)                 
        WRITE (ELABS(3)(22:30),'(''='',F8.4)')  DSUM(3)+DSUM(8)                 
        WRITE (ELABS(4)(16:30),'(''='',2F7.4)') DSUM(4),DSUM(9)                 
        WRITE (ELABS(5)(16:30),'(''='',2F7.4)') DSUM(5),DSUM(10)                
c
        write(nview,'(a,1x,g12.5)') 'EROSION SCALING=',ero_scale
c
c       Change FLAB to EROLAB
c
        CALL LIM_DRAW (ODOUTS,ODWIDS,NERODS,MAXOS,NOS,ANLY,                         
     >    5,99,vmin,vmax,-BIGG,BIGG,IGTS,ITEC,AVS,NAVS,                
     >    JOB,TITLE,YLAB3,EROLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6)                 
        nview = ' '
c
c     Write out the erosion data
c        
        write(6,'(a)') 'EROSION ALONG LIMITER SURFACE:'
        write(6,'(4(2x,a20))') 'DISTANCE','TOTAL DEPOSITION',
     >            'TOTAL REMOVAL','NET EROSION' 
c
        do io = 1,nos
           write(6,'(4(4x,g18.8))') odouts(io),nerods(io,1),
     >            nerods(io,3),nerods(io,4)

        end do

c        write(6,'(a)') 'EROSION ALONG LIMITER SURFACE:'
c        write(6,'(4(2x,a20))') 'DISTANCE','TOTAL DEPOSITION',
c     >            'TOTAL REMOVAL','NET EROSION' 
c        do io = 1,nos
c           write(6,'(8(4x,g18.8))') odouts(io),
c     >            nerods(io,1),sum(nerods3(io,:,1)),
c     >            nerods(io,3),sum(nerods3(io,:,3)),
c     >            nerods(io,4),sum(nerods3(io,:,4))
c
c        end do

c
c       Write out the poloidally resolved data
c
c        write(6,'(a)') 'POLOIDALLY RESOLVED EROSION ALONG'//
c     >                 ' LIMITER SURFACE:'
c        write(6,'(a)')
c        write(6,'(a)') 'TOTAL DEPOSITION:'
c        write(6,'(13x,100(1x,g12.5))') 
c     >              ((ps(ip)-pwids(ip)/2.0),ip=-maxnps,maxnps)
c        do io = 1,nos
c           write(6,'(101(1x,g12.5))') 
c     >             odouts(io),(nerods3(io,ip,1),ip=-maxnps,maxnps)
c        end do   
c        write(6,'(a)')
c        write(6,'(a)') 'TOTAL REMOVAL:'
c        do io = 1,nos
c           write(6,'(101(1x,g12.5))') 
c     >             odouts(io),(nerods3(io,ip,3),ip=-maxnps,maxnps)
c        end do   
c        write(6,'(a)')
c        write(6,'(a)') 'NET EROSION:'
c        do io = 1,nos
c           write(6,'(101(1x,g12.5))') 
c     >             odouts(io),(nerods3(io,ip,4),ip=-maxnps,maxnps)
c        end do   
c        write(6,'(a)')
c
C                                                                               
      ELSEIF (BREF.EQ.'2WY') THEN                                               
C     ===========================                                               
        REF = 'WALL DEPOSITION'                                                 
        IF (IFOLD.EQ.1) THEN                                                    
          DO 1150 IY = 1, NYS                                                   
            DO 1150 IZ = -2, MAXIZ+1                                            
              WALLS(IY,IZ) = WALLS(IY,IZ) + WALLS(-IY,IZ)                       
 1150     CONTINUE                                                              
        ENDIF                                                                   

        write(6,'(a)') '2WY:'
        write(6,'(a,3i5)') 'WALLS:',nys,maxiz
        do iy = -nys,nys
           write(6,'(a,10(1x,g12.5))') 'WALLS:',youts(iy),
     >          (walls(iy,iz),iz=-2,maxiz+1)
        end do

c
        CALL LIM_DRAW (YOUTS,YWIDS,WALLS,2*MAXNYS+1,2*MAXNYS+1,ANLY,                
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,YLAB,FLAB,ZLABS(-2),REF,YVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  POWER LOSS                                                                   
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2RY') THEN                                               
C     ===========================                                               
        CALL RINTX (POWLS,POWL3,IPLANE,MAXIZ,XINTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM,FP,FT,1.0,IXMIN,IXMAX)                       
        CALL PROJEC (MAXIZ,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,CSINTB)        
        REF = 'POWER LOSS '//INTREF // COMMAP                          
        CALL LIM_DRAW (YOUTS,YWIDSS,XINTS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,YLAB,FLAB,ZLABS(-2),REF,YVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2RX') THEN                                               
C     ===========================                                               
        CALL RINTY (POWLS,POWL3,IPLANE,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM,FP,FT,1.0,IYMIN,IYMAX)                  
        REF = 'POWER LOSS '//INTREF // COMMAP                          
        CALL LIM_DRAW (XOUTS,XWIDS,YINTS,MAXNXS,NXS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,XLAB,FLAB,ZLABS(-2),REF,XVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
      ELSEIF (BREF.EQ.'2RP') THEN                                               
C     ===========================                                               
        CALL RINTXY (POWLS,POWL3,IPLANE,MAXIZ,XYINTS,IFOLD,RV,XV,YV,0,          
     >               SSS,NYSLIM,FP,FT)                                          
        REF = 'POWER LOSS INT''D OVER X,Y'                                      
        CALL LIM_DRAW (POUTS,PWIDS,XYINTS,2*MAXNPS+1,2*MAXNPS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,PLAB,FLAB,ZLABS(-2),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  LINE RADIATION LOSS                                                          
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2LY') THEN                                               
C     ===========================                                               
        CALL RINTX (LINES,LINE3,IPLANE,MAXIZ,XINTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM,FP,FT,1.0,IXMIN,IXMAX)                   
        CALL PROJEC (MAXIZ,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,CSINTB)        
        REF = 'LINE RADIATION '// INTREF // COMMAP                          
        CALL LIM_DRAW (YOUTS,YWIDSS,XINTS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,YLAB,FLAB,ZLABS(-2),REF,YVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2LX') THEN                                               
C     ===========================                                               
        CALL RINTY (LINES,LINE3,IPLANE,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM,FP,FT,1.0,IYMIN,IYMAX)                      
        REF = 'LINE RADIATION '//INTREF // COMMAP                          
        CALL LIM_DRAW (XOUTS,XWIDS,YINTS,MAXNXS,NXS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,XLAB,FLAB,ZLABS(-2),REF,XVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
      ELSEIF (BREF.EQ.'2LP') THEN                                               
C     ===========================                                               
        CALL RINTXY (LINES,LINE3,IPLANE,MAXIZ,XYINTS,IFOLD,RV,XV,YV,0,          
     >               SSS,NYSLIM,FP,FT)                                          
        REF = 'LINE RADIATION INT''D OVER X,Y'                                  
        CALL LIM_DRAW (POUTS,PWIDS,XYINTS,2*MAXNPS+1,2*MAXNPS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,PLAB,FLAB,ZLABS(-2),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  PLRPS:  DON'T PLOT LINES REFERRING TO IZ STATES > REQUESTED                  
C  INTEGRATION ROUTINES CALLED WITH -2 FOR "MINIZ", MEANING THAT THE            
C  SECONDARY NEUTRALS DETAILS HAVE ALREADY BEEN DEALT WITH.                     
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2PY') THEN                                               
C     ===========================                                               
        CALL RINTX (PLRPS,PLRP3,IPLANE,MLS-2,XINTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM,FP,FT,1.0,IXMIN,IXMAX)                          
        CALL PROJEC (MLS-2,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,CSINTB)        
        REF = 'PLRPS '//INTREF // COMMAP                               
        CALL LIM_DRAW (YOUTS,YWIDSS,XINTS(-MAXNYS,-1),2*MAXNYS+1,
     >    2*MAXNYS+1,        
     >    ANLY,MLS,KSMOTH,VMIN,VMAX,0.0,BIGG,IGLS,ITEC,AVS,NAVS,                
     >    JOB,TITLE,YLAB,FLAB,PLABS(1),REF,YVIEW,PLANE,TABLE,IPLOT,IFL)         
C                                                                               
      ELSEIF (BREF.EQ.'2PX') THEN                                               
C     ===========================                                               
        CALL RINTY (PLRPS,PLRP3,IPLANE,MLS-2,YINTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM,FP,FT,1.0,IYMIN,IYMAX)                  
        REF = 'PLRPS '//INTREF // COMMAP                               
        CALL LIM_DRAW (XOUTS,XWIDS,YINTS(1,-1),MAXNXS,NXS,ANLY,                     
     >    MLS,KSMOTH,VMIN,VMAX,0.0,BIGG,IGLS,ITEC,AVS,NAVS,                     
     >    JOB,TITLE,XLAB,FLAB,PLABS(1),REF,XVIEW,PLANE,TABLE,IPLOT,2)           
C                                                                               
      ELSEIF (BREF.EQ.'2PP') THEN                                               
C     ===========================                                               
        CALL RINTXY (PLRPS,PLRP3,IPLANE,MLS-2,XYINTS,IFOLD,RV,XV,YV,0,          
     >               SSS,NYSLIM,FP,FT)                                          
        REF = 'PLRPS INTEGRATED OVER X,Y'                                       
        CALL LIM_DRAW (POUTS,PWIDS,XYINTS(-MAXNPS,-1),2*MAXNPS+1,
     >    2*MAXNPS+1,        
     >    ANLY,MLS,KSMOTH,VMIN,VMAX,0.0,BIGG,IGLS,ITEC,AVS,NAVS,                
     >    JOB,TITLE,PLAB,FLAB,PLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,2)           
C                                                                               
C-----------------------------------------------------------------------        
C  TIZS: POSITIONS AT IONISATION                                                
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2IY') THEN                                               
C     ===========================                                               
        CALL RINTX (TIZS,TIZ3,IPLANE,MAXIZ,XINTS,IFOLD,RV,XV,YV,IVU,            
     >              SSS,NYSLIM,FP,FT,1.0,IXMIN,IXMAX)                    
        CALL PROJEC (MAXIZ,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,CSINTB)        
        REF = 'IONISATION '//INTREF // COMMAP                              
        CALL LIM_DRAW (YOUTS,YWIDSS,XINTS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,YLAB,FLAB,ZLABS(-2),REF,YVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2IX') THEN                                               
C     ===========================                                               
        CALL RINTY (TIZS,TIZ3,IPLANE,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,            
     >              SSS,NYSLIM,FP,FT,1.0,IYMIN,IYMAX)                      
        REF = 'IONISATION '//INTREF // COMMAP                              
        CALL LIM_DRAW (XOUTS,XWIDS,YINTS,MAXNXS,NXS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,XLAB,FLAB,ZLABS(-2),REF,XVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
        DO 1910 IZ = -2, MAXIZ                                                  
          SUM1(IZ) = 0.0                                                        
          SUM2(IZ) = 0.0                                                        
          DO 1900 IX = 1, NXS                                                   
            SUM1(IZ) = SUM1(IZ) + XWIDS(IX) * YINTS(IX,IZ) * XOUTS(IX)          
            SUM2(IZ) = SUM2(IZ) + XWIDS(IX) * YINTS(IX,IZ)                      
 1900     CONTINUE                                                              
          IF (SUM2(IZ).LE.0.0) SUM2(IZ) = 1.E-20                                
 1910   CONTINUE                                                                
        WRITE (6,9010) (ZLABS1(IZ)(5:11),SUM1(IZ)/SUM2(IZ),IZ=-2,MAXIZ)         
C                                                                               
      ELSEIF (BREF.EQ.'2IP') THEN                                               
C     ===========================                                               
        CALL RINTXY (TIZS,TIZ3,IPLANE,MAXIZ,XYINTS,IFOLD,RV,XV,YV,0,            
     >               SSS,NYSLIM,FP,FT)                                          
        REF = 'IONISATION INT''D OVER X,Y'                                      
        CALL LIM_DRAW (POUTS,PWIDS,XYINTS,2*MAXNPS+1,2*MAXNPS+1,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,PLAB,FLAB,ZLABS(-2),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  Z EFFECTIVE ETC.                                                             
C-----------------------------------------------------------------------        
C                                                                               
C---- VIEWPOINT NOT REQUIRED ... USE "NVIEW" STRING WHEN PLOTTING               
C---- INTEGRAL OVER Y INSTEAD OF "XVIEW/YVIEW" IN CALL TO DRAW,                 
C---- AND CALL INTEG WITH "IVU=0" FOR OFF.                                      
C---- THERE IS NO 3D EQUIVALENT OF ZEFFS: PUT DUMMY ARRAY TIZ3 IN               
C---- ARGUMENT LIST FOR RINTX/Y FOLLOWED BY PLANE (WHICH WILL BE 99).           
C---- NOTE 291: ADD DILUTION TO X PLOTS                                         
C                                                                               
      ELSEIF (BREF.EQ.'2ZY') THEN                                               
C     ===========================                                               
        CALL RINTX (ZEFFS,TIZ3,99,3-2,XINTS,IFOLD,RV,XV,YV,0,                   
     >              SSS,NYSLIM,FP,FT,CA-CAW,IXMIN,IXMAX)                    
        REF = 'ZEFFS AVERAGED OVER X'                                           
        CALL LIM_DRAW (YOUTS,YWIDSS,XINTS(-MAXNYS,-1),2*MAXNYS+1,
     >    2*MAXNYS+1,        
     >    ANLY,3,99,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                      
     >    JOB,TITLE,YLAB,FLAB,ELABS(6),REF,NVIEW,PLANE,TABLE,IPLOT,IFL)         
C                                                                               
      ELSEIF (BREF.EQ.'2ZX') THEN                                               
C     ===========================                                               
        CALL RINTY (ZEFFS,TIZ3,99,3-2,YINTS,IFOLD,RV,XV,YV,0,                   
     >              SSS,NYSLIM,FP,FT,CTWOL,IYMIN,IYMAX)                  
        REF = 'ZEFFS AVERAGED OVER Y'                                           
        DO 1950 IX = 1, NXS                                                     
          YINTS(IX,2) = YINTS(IX,0) / CRNBS(IX,1)                               
 1950   CONTINUE                                                                
        ELABS(9) = '    DILUTE'                                                 
        CALL LIM_DRAW (XOUTS,XWIDS,YINTS(1,-1),MAXNXS,NXS,ANLY,                     
     >    4,99,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                           
     >    JOB,TITLE,XLAB,FLAB,ELABS(6),REF,NVIEW,PLANE,TABLE,IPLOT,2)           
        ELABS(9) = '    POW TOT'                                                
C                                                                               
C-----------------------------------------------------------------------        
C  AVERAGE CLOUD TEMPERATURES                                                   
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2TY') THEN                                               
C     ===========================                                               
        REF = 'TEMPERATURE AVERAGED OVER X'                                     
        CALL LIM_DRAW (YOUTS,YWIDSS,SDTYS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXIZ,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(1),ITEC,AVS,NAVS,                
     >    JOB,TITLE,YLAB,FFLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2TX') THEN                                               
C     ===========================                                               
        REF = 'TEMPERATURE AVERAGED OVER Y'                                     
        CALL LIM_DRAW (XOUTS,XWIDS,SDTXS,MAXNXS,NXS,ANLY,                           
     >    MAXIZ,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(1),ITEC,AVS,NAVS,                
     >    JOB,TITLE,XLAB,FFLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  AVERAGE SPECTROSCOPICALLY-WEIGHTED CLOUD TEMPERATURES   
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2SY') THEN                                               
C     ===========================                                               
        DO 1947 IX = 1,NLS
          IF ((PIZS(IX).EQ.0).OR.(PIZS(IX).EQ.-1).OR.
     >        (PIZS(IX).EQ.-2))  THEN
             IGSS(IX) = 0
          ELSE
             WRITE(PLABS(IX),'(A1,I3,F7.2)') 'S',
     >          PIZS(IX),PLAMS(IX) 
             IGSS(IX) = IGLS(IX)
          ENDIF     
 1947   CONTINUE
C
        REF = 'SPECT.TEMP. AVERAGED OVER X'                                 
        CALL LIM_DRAW (YOUTS,YWIDSS,SCTYS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXNLS,ISMOTH,VMIN,VMAX,0.0,BIGG,IGSS(1),ITEC,AVS,NAVS,              
     >    JOB,TITLE,YLAB,FFLAB,PLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2SX') THEN                                               
C     ===========================                                               
        DO 1948 IX = 1,NLS
          IF ((PIZS(IX).EQ.0).OR.(PIZS(IX).EQ.-1).OR.
     >        (PIZS(IX).EQ.-2)) THEN
             IGSS(IX) = 0
          ELSE
             WRITE(PLABS(IX),'(A1,I3,F7.2)') 'S',
     >          PIZS(IX),PLAMS(IX) 
             IGSS(IX) = IGLS(IX)
          ENDIF     
 1948   CONTINUE
C
        REF = 'SPECT.TEMP. AVERAGED OVER Y'                                     
        CALL LIM_DRAW (XOUTS,XWIDS,SCTXS,MAXNXS,NXS,ANLY,                           
     >    MAXNLS,ISMOTH,VMIN,VMAX,0.0,BIGG,IGSS(1),ITEC,AVS,NAVS,             
     >    JOB,TITLE,XLAB,FFLAB,PLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  AVERAGE Y STEPS                                                              
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2YY') THEN                                               
C     ===========================                                               
        REF = 'DELTAY STEPS AVERAGED OVER X'                                    
        CALL LIM_DRAW (YOUTS,YWIDSS,SDYYS,2*MAXNYS+1,2*MAXNYS+1,ANLY,               
     >    MAXIZ,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(1),ITEC,AVS,NAVS,                
     >    JOB,TITLE,YLAB,FFLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,IFL)        
C                                                                               
      ELSEIF (BREF.EQ.'2YX') THEN                                               
C     ===========================                                               
        REF = 'DELTAY STEPS AVERAGED OVER Y'                                    
        CALL LIM_DRAW (XOUTS,XWIDS,SDYXS,MAXNXS,NXS,ANLY,                           
     >    MAXIZ,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(1),ITEC,AVS,NAVS,                
     >    JOB,TITLE,XLAB,FFLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C  DEBUG PARTICLE TRAJECTORIES                                                              
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'2AD') THEN                                               
C     ===========================                                               
c
c       Set plot titles, labels and comments  
c
        ref = 'Particle Trajectories'
        CALL LIM_GRTSET (TITLE,REF,NVIEW,PLANE,JOB,vmin,vmax,                    
     >    grimin,grimax,TABLE,YLAB1,XLAB,2,SMOOTH,1,ANLY,maxiz+1)                    
c
c       Plot limiter edge 
c
        DO 1660 IQX = 1-NQXSO, 0                                                 
          QXFUNS(IQX-1,1) =-QEDGES(IQX,1)                                       
          QXFUNS(1-IQX,1) = QEDGES(IQX,2)                                       
          QXFUNS(IQX-1,2) = QXS(IQX)                                            
          QXFUNS(1-IQX,2) = QXS(IQX)                                            
 1660   CONTINUE                                                                
        CALL PRC ('GRAPH OF LIMITER EDGE Y(X) WILL BE PLOTTED')               
        QXFUNS(0,1) = 0.0                                                     
        QXFUNS(0,2) = 0.0                                                     
        CALL LIM_GRTRAC (QXFUNS(-NQXSO,1),QXFUNS(-NQXSO,2),                       
     >      2*NQXSO+1,ELABS(15),'LINE')                                         
c
c       Plot individual particle trajectories - using different lines
c       It is necessary to reverse the indices because the "X-axis" 
c       as defined in the LIM simulation is the vertical axis from 
c       plasma centre to the wall ... which in plotting coordinates 
c       is the y-axis.
c
        write(6,*) 'Trajectories:',maxiz
        do 1665 ix = 1,maxiz
           write(6,*) 'length:',ptracl(ix)
           do 1670 iy = 1,ptracl(ix)
              tptracx(iy) = ptracs(iy,ix,2)  
              tptracy(iy) = ptracs(iy,ix,1)  
              write(6,*) 'point:',iy,tptracx(iy),tptracy(iy)
 1670      continue
           write(label1,1680) ix,ix  
           write(6,*) 'label1:',label1
           CALL LIM_GRTRAC (tptracx,tptracy,                       
     >         ptracl(ix),label1 ,'LINE')                                         
 1665   continue       
 1680   format('P ',I2,'Particle',i3)  
c
c       Finish plot 
c
        call frame 
c
      ENDIF                                                                     
      GOTO 1000                                                                 
C     *********                                                                 
C                                                                               
C#######################################################################        
C    2 DIMENSIONAL R,THETA PLOTS                                             
C    ===========================                                             
C  SOMETIMES A PARAMETER WILL BE CHANGED - MARK IT WITH A '*' IN OUTPUT         
C  LIST FOR POSTERITY.                                                          
C#######################################################################        
C                                                                               
 1958 CONTINUE                                                                  
      WRITE (7,9220)                                                            
 9220 FORMAT(' ',//2X,                                                          
     >            'REF   TITLE        PLOTTING RANGE INTEG.RANGE PLOT',
     >  ' SMOOTH MAXIZ',        
     >  ' PLANE YFLD ALLZ VU',/1X,71('-'))                                      
 1960 CONTINUE                                                                  
      CALL RDGRT (GRAPH,VMIN,VMAX,SMIN,SMAX,IPLOT,JSMOTH,MAXIZ,
     >            IPLANE,IFOLD,IALL,IVU,'LINE OF GRAPH DETAILS',IERR)           

      BREF = GRAPH(1:3)                                                         
c      write(0,*) '2D RT:',iplot,trim(bref)

      IF (IERR.NE.0)        GOTO 9999                                           
      IF (BREF(1:1).NE.'R'.and.bref(1:3).ne.'000') GOTO 2000                                           
      IF (IPLOT.EQ.0)       GOTO 1960                                           
C                                                                               
      DO 1965 J = 1, 7                                                          
        S(J) = ' '                                                              
 1965 CONTINUE                                                                  
C
C     SOMETIMES NO PLRP PLOTS
C
      IF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') THEN                               
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- SMOOTHING NOT ALLOWED SOMETIMES...                                        
C---- NOT ENOUGH POINTS TO SMOOTH P PLOTS  (ONLY 9)                             
C                                                                               
      IF (JSMOTH.NE.99.AND.(BREF(2:2).EQ.'Z')) THEN                          
        JSMOTH = 99                                                             
        S(2) = '*'                                                             
      ENDIF                                                                     
C                                                                               
C---- CAN'T PLOT ABOVE MAX RECORDED IONISATION STATE!                           
C                                                                               
      IF (MAXIZ.GT.NIZS) THEN                                                   
        MAXIZ = NIZS                                                            
        S(3) = '*'                                                              
      ELSEIF (MAXIZ.LT.-2) THEN                                                 
        MAXIZ = -2                                                              
        S(3) = '*'                                                              
      ENDIF                                                                     
C
C     FOR THESE TYPES OF PLOTS THERE IS CURRENTLY NO 3D SUPPORT
C     SO SET IPLANE TO 99 ALWAYS
C
      IF (IPLANE.NE.99) THEN 
         IPLANE = 99
         S(4) = '*'
      ENDIF
C
C
      PLANE = ' '
      COMMAP = ',P'
C
C     IF (IPLANE.EQ.99) THEN                                                    
C       PLANE  = ' '                                                            
C       COMMAP = ',P'                                                           
C     ELSEIF (IPLANE.GE.-MAXNPS.AND.IPLANE.LE.MAXNPS) THEN                      
C       WRITE (PLANE,'(4A)') 'PLOT FOR  ',PRINPS(IPLANE-1),                     
C    >    ' < P <',PRINPS(IPLANE)                                               
C       COMMAP = ' '                                                            
C     ELSE                                                                      
C       IPLANE = 99                                                             
C       PLANE = ' '                                                             
C       COMMAP = ',P'                                                           
C       S(4) = '*'                                                              
C     ENDIF                                                                     
C
      ANLY   = ' '                                                              
      IF (CANAL.LT.CA .AND. IPLANE.EQ.99 .AND. (BREF(2:2).EQ.'C'.OR.            
     >  BREF(2:2).EQ.'R'.OR.BREF(2:2).EQ.'L'.OR.BREF(2:2).EQ.'P'.OR.            
     >  BREF(2:2).EQ.'Z') .AND. (BREF(3:3).EQ.'R'.OR.BREF(3:3).EQ.'A'))         
     >  WRITE (ANLY,'(''ANALYTIC EXTENSION FOR R <'',F6.3)') CA-CANAL         
C                                                                               
C---- SEPARATE RESULTS FOR -Y, +Y NOT ALWAYS AVAILABLE                          
C---- FOLDING N/A FOR X GRAPHS (EXCEPT DEPS/NETE), P GRAPHS                     
C                                                                               
C
C---- THE THETA DEPENDENCE IS ONLY MAPPED TO 0 -> 2L, HOWEVER Y-FOLDING
C---- IS STILL OBTAINED IN THE RTXY MAPPING ROUTINE. 
C
C     ALLOW IFOLD TO BE PASSED TO ALL R,THETA ROUTINES, IT WILL
C     NOT AFFECT RESULTS AND IS SIGNIFICANT FOR THE THETA AND PSI 
C     INTEGRATION ROUTINES.
C
C     IF (IFOLD.NE.0.AND.BREF(3:3).NE.'A') THEN
C       IFOLD = 0                                                               
C       S(5) = '*'                                                              
C     ENDIF                                                                     
C
      IF (IFOLD.NE.0.AND.IFOLD.NE.1) THEN                                       
        IFOLD = 1                                                               
        S(5) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- "ALL IONISATIONS LINE" NOT ALWAYS AVAILABLE - CAN ONLY BE INCLUDED        
C---- IF WE ARE INTEGRATING OVER ALL THE AVAILABLE IONISATION STATES,           
C---- IE. NIZS=MAXIZ.  NOTE THAT THE IGNORE OPTION CAN OF COURSE BE USED        
C---- TO SUPPRESS THE ACTUAL PLOTTING OF INTERMEDIATE LINES IF REQUIRED.        
C                                                                               
      IF (IALL.EQ.1.AND.(                     
     >  BREF(2:2).EQ.'Z'.OR.BREF(2:2).EQ.'P'            
     >  .OR.MAXIZ.LT.NIZS)) THEN                                
        IALL = 0                                                                
        S(6) = '*'                                                              
      ENDIF                                                                     
      IF (IALL.NE.0.AND.IALL.NE.1) THEN                                         
        IALL = 0                                                                
        S(6) = '*'                                                              
      ENDIF                                                                     
C                        
C     VIEWPOINT FLAG:
C     THIS IS NOT USED FOR R,THETA PLOTS- INSTEAD THE PSI COORDINATE 
C     PLOT TYPE IS SELECTED.  
C                                                       
      IF (IVU.NE.0.) THEN
        IVU = 0
        S(7) = '*'
      ENDIF
C                                                                               
C     SET UP APPROPRIATE GRAPH LABELS
C
      RTVIEW = 'PLASMA CENTRE R,THETA'                                         
      WRITE(SVIEW,'(''VIEWED FROM: '',F6.2,'' M, '',F6.2,'' DEGS'')') 
     >      RVPOS,TVPOS*RADDEG                                            
      RLAB = 'RADIAL DISTANCE (M)'
      THLAB = 'ANGULAR POSITION (DEG)'
      SLAB  = 'VIEW ANGLE PSI (DEG)'                                         
C
C
C      FOR NOW IFL IS CODED AS 2 SINCE THE THETA AXIS ONLY 
C      OCCURS ONCE (I.E. 0->2L NOT -2L->+2L AS FOR Y)
C
C      IFL = 2                                                                  
C      IF (IVU.EQ.0.AND.IPLANE.EQ.99) IFL = 7                                   
C                                                                               
C                                                                               
C---- ADJUST VMIN,VMAX AND SMIN,SMAX IF REQUIRED  
C
C---- FOR S PLOTS THE VMIN,VMAX VALUES MAY BE MODIFIED AGAIN TO COINCIDE 
C---- WITH PSIMIN AND PSIMAX
C                                                             
C                                        
      IF (BREF(3:3).EQ.'R') THEN                                            
        VMIN = MAX (VMIN, 0.0)
        VMAX = MIN (VMAX, CA-CAW )                                              
C
C     ANGULAR LIMITS OF INTEGRATION
C
        IF (SMIN.LT.0.0) THEN
          SMIN = MAX (SMIN, -PI)
          SMAX = MIN (SMAX,  PI)   
        ELSE 
          SMIN = MAX (SMIN, 0.0)
          SMAX = MIN (SMAX,2*PI)   
        ENDIF
      ELSEIF (BREF(3:3).EQ.'A') THEN                                            
        IF (VMIN.LT.0.0) THEN 
C
C       ALLOW THETA PLOTTING ON -PI,+PI OR ON 0,2*PI ONLY
C
          VMIN = MAX (VMIN,-PI)                                
          VMAX = MIN (VMAX,PI)
        ELSE
          VMIN = MAX (VMIN,0.0)
          VMAX = MIN (VMAX,2*PI)  
        ENDIF
C
C       RADIAL LIMITS OF INTEGRATION
C
        SMIN = MAX (SMIN, 0.0)
        SMAX = MIN (SMAX, CA-CAW)
C
C     NOTE: NO INTEGRATION LIMITS ALLOWED FOR S-PLOTS
C
      ELSEIF (BREF(3:3).EQ.'S') THEN
C
C       ALLOW PSI PLOTTING ON -PI,+PI FROM VIEWING AXIS ONLY
C       THESE MAY BE ADJUSTED TO CALCULATED PSIMIN AND PSIMAX VALUES  
C
        VMIN = MAX(VMIN,-PI)
        VMAX = MIN(VMAX,PI)
      ENDIF                                                                     
C     
C     MAKE SURE SMIN < SMAX , IF NOT SET THEM EQUAL
C
      IF (SMAX.LT.SMIN) THEN 
         SMAX = SMIN 
         WRITE (6,'(A32)') 'SMIN,SMAX INVALID: ADJUSTED'
      ENDIF
C
C     USE THE PLANE STRING VARIABLE TO CONTAIN THE INTEGRATION RANGE
C     INFORMATION FOR THE PLOT
C
      IF (BREF(3:3).NE.'S') THEN 
         WRITE(PLANE,'(''INTEGRATION RANGE:'',F7.2,'' TO '',F7.2)')
     >      SMIN*RADDEG,SMAX*RADDEG
      ENDIF
C                                                                               
C---- WRITE OUT (POSSIBLY ADJUSTED) LINE OF GRAPH DETAILS                       
C                                                                               
      WRITE (7,'(2X,A21,4F6.2,I4,A1,I6,5(A1,I4),A1)')                           
     >  GRAPH,VMIN,VMAX,SMIN,SMAX,IPLOT,                                       
     >  S(1),JSMOTH,S(2),MAXIZ,S(3),IPLANE,S(4),IFOLD,S(5),IALL,S(6),           
     >  IVU,S(7)                                                                
      ISMOTH = JSMOTH + 3                                                       
      IF (ISMOTH.EQ.3) ISMOTH = 1                                               
      IF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') GOTO 1960                              
      IF (VMIN.GE.VMAX)                  GOTO 1960                              
      IF (IPLOT.EQ.0)                    GOTO 1960                              
C                                                                               
C---- SET NUMBER OF LINES TO PLOT FOR PLRPS CASES                               
C---- AND SET SMOOTHING PARAMETER BASED ON THIS.                                
C---- SET SSS FLAG FOR NOTE 172 DISTANCE SQUARED FACTOR FOR XVIEW PLOTS         
C                                                                               
      SSS = .FALSE.                                                             
      IF (BREF(2:2).EQ.'P') THEN                                                
        SSS = .TRUE.                                                            
        MLS = NLS                                                               
        DO 1970 IL = 1, NLS                                                     
          PLABS(IL)(22:30) = '         '                                        
          IF (PIZS(IL).GT.MAXIZ) MLS = MLS - 1                                  
 1970   CONTINUE                                                                
        KSMOTH = MLS + 1                                                        
        DO 1975 IL = MLS, 1, -1                                                 
          IF (PIZS(IL).GE.JSMOTH) KSMOTH = IL                                   
          IF (PIZS(IL).LT.0 .AND. JSMOTH.EQ.0) KSMOTH = IL                      
 1975   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      DO 1980 IZ = -2, MAXIZ                                                    
        ZLABS(IZ) = ZLABS1(IZ)                                                  
 1980 CONTINUE                                                                  
      ZLABS(MAXIZ+1) = ZLABS1(NIZS+1)                                           
C
C     RESET THE TOUTS ARRAY TO ITS ORIGINAL VALUES IN CASE IT HAS 
C     BEEN MODIFIED BY PREVIOUS CALLS TO THIS ROUTINE
C
      DO 1981 IT = 1,MAXNAS
         TOUTS(IT) = TOUTS1(IT)
         TWIDS(IT) = TWIDS1(IT)
1981  CONTINUE 
C                                                                               
C-----------------------------------------------------------------------        
C  CLOUDS                                                                       
C-----------------------------------------------------------------------        
C                                                                               
      IF     (BREF.EQ.'RCA') THEN                                               
C     ===========================                                               
        CALL RINTR (SDLIMS,IPLANE,MAXIZ,RINTS,IFOLD,        
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                
        REF = 'IMPURITY INTEGRATED OVER R' // COMMAP
C
C       IT IS NECESSARY TO ADJUST THE OUTPUT ARRAY IF THE PLOTTING RANGE 
C       STARTS AT NEGATIVE THETA. THIS IS DONE SO THAT THE LIMITER WHICH
C       IS AT THETA = 0 MAY BE CENTRED ON THE PLOT. BOTH THE INTEGRATED
C       VALUES (RINTS) AND THE COORDINATES (TOUTS) ARE ADJUSTED
C
        IF (VMIN.LT.0.0) THEN 
          CALL RTORD (RINTS,VMIN)
        ENDIF
C
C       CONVERT TOUTS COORDINATE AXIS, AFTER ADJUSTMENT, BACK TO DEGREES
C
        CALL CVRTDEG ( TOUTS,NAS,VMIN,VMAX)
C                             
        CALL LIM_DRAW (TOUTS,TWIDS,RINTS,MAXNAS,NAS,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,THLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)        
C                                                                               
      ELSEIF (BREF.EQ.'RCR') THEN                                               
C     ===========================                                               
        CALL RINTT (SDLIMS,IPLANE,MAXIZ,TINTS,IFOLD,        
     >              SSS,FP,FT,1.0,SMIN,SMAX)
        REF = 'IMPURITY INTEGRATED OVER A' // COMMAP                            
        CALL LIM_DRAW (ROUTS,RWIDS,TINTS,MAXNRS,NRS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,RLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)         
C                                                                               
      ELSEIF (BREF.EQ.'RCS') THEN                                               
C     ===========================                                               
        CALL RINTPSI (SDLIMS,IPLANE,MAXIZ,SINTS,IFOLD,        
     >              SSS,FP,FT,1.0,VMIN,VMAX,SMIN,SMAX)                        
        REF = 'IMPURITY INTEGRATED FROM S'
           
        CALL LIM_DRAW (SOUTS1,SWIDS,SINTS,MAXNSS,NSS,ANLY,                      
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,SLAB,FLAB,ZLABS(-2),REF,SVIEW,PLANE,TABLE,IPLOT,2)          
C
C-----------------------------------------------------------------------        
C  POWER LOSS                                                                   
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'RRA') THEN                                               
C     ===========================                                               
        CALL RINTR (POWLS,IPLANE,MAXIZ,RINTS,IFOLD,          
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                  
        REF = 'POWER LOSS INTEGRATED OVER R' // COMMAP                          
C                                                                               
C
C       IT IS NECESSARY TO ADJUST THE OUTPUT ARRAY IF THE PLOTTING RANGE 
C       STARTS AT NEGATIVE THETA. THIS IS DONE SO THAT THE LIMITER WHICH
C       IS AT THETA = 0 MAY BE CENTRED ON THE PLOT. 
C
        IF (VMIN.LT.0.0) THEN 
          CALL RTORD (RINTS,VMIN)
        ENDIF
C
C       CONVERT TOUTS COORDINATE AXIS, AFTER ADJUSTMENT, BACK TO DEGREES
C
        CALL CVRTDEG ( TOUTS,NAS,VMIN,VMAX)
C                             
        CALL LIM_DRAW (TOUTS,TWIDS,RINTS,MAXNAS,NAS,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,THLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)        
C
      ELSEIF (BREF.EQ.'RRR') THEN                                               
C     ===========================                                               
        CALL RINTT (POWLS,IPLANE,MAXIZ,TINTS,IFOLD,          
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                 
        REF = 'POWER LOSS INTEGRATED OVER A' // COMMAP                          
        CALL LIM_DRAW (ROUTS,RWIDS,TINTS,MAXNRS,NRS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,RLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)      
C
      ELSEIF (BREF.EQ.'RRS') THEN                                               
C     ===========================                                               
        CALL RINTPSI (POWLS,IPLANE,MAXIZ,SINTS,IFOLD,          
     >              SSS,FP,FT,1.0,VMIN,VMAX,SMIN,SMAX)
        REF = 'POWER LOSS INTEGRATED FROM S'                          
        CALL LIM_DRAW (SOUTS1,SWIDS,SINTS,MAXNSS,NSS,ANLY,                         
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,SLAB,FLAB,ZLABS(-2),REF,SVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C  LINE RADIATION LOSS                                                          
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'RLA') THEN                                               
C     ===========================                                               
        CALL RINTR (LINES,IPLANE,MAXIZ,RINTS,IFOLD,          
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                 
        REF = 'LINE RADIATION INT''D OVER R' // COMMAP                          
C                                                                               
C
C       IT IS NECESSARY TO ADJUST THE OUTPUT ARRAY IF THE PLOTTING RANGE 
C       STARTS AT NEGATIVE THETA. THIS IS DONE SO THAT THE LIMITER WHICH
C       IS AT THETA = 0 MAY BE CENTRED ON THE PLOT. 
C
        IF (VMIN.LT.0.0) THEN 
          CALL RTORD (RINTS,VMIN)
        ENDIF
C
C       CONVERT TOUTS COORDINATE AXIS, AFTER ADJUSTMENT, BACK TO DEGREES
C
        CALL CVRTDEG ( TOUTS,NAS,VMIN,VMAX)
C                             
        CALL LIM_DRAW (TOUTS,TWIDS,RINTS,MAXNAS,NAS,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,THLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)        
C                                                                               
      ELSEIF (BREF.EQ.'RLR') THEN                                               
C     ===========================                                               
        CALL RINTT (LINES,IPLANE,MAXIZ,TINTS,IFOLD,          
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                  
        REF = 'LINE RADIATION INT''D OVER T' // COMMAP                          
        CALL LIM_DRAW (ROUTS,RWIDS,TINTS,MAXNRS,NRS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,RLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)       
C                                                                               
      ELSEIF (BREF.EQ.'RLS') THEN                                               
C     ===========================                                               
        CALL RINTPSI (LINES,IPLANE,MAXIZ,SINTS,IFOLD,          
     >              SSS,FP,FT,1.0,VMIN,VMAX,SMIN,SMAX)
        REF = 'LINE RADIATION INT''D FROM S'
        CALL LIM_DRAW (SOUTS1,SWIDS,SINTS,MAXNSS,NSS,ANLY,                        
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,SLAB,FLAB,ZLABS(-2),REF,SVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C  PLRPS:  DON'T PLOT LINES REFERRING TO IZ STATES > REQUESTED                  
C  INTEGRATION ROUTINES CALLED WITH -2 FOR "MINIZ", MEANING THAT THE            
C  SECONDARY NEUTRALS DETAILS HAVE ALREADY BEEN DEALT WITH.                     
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'RPA') THEN                                               
C     ===========================                                               
        CALL RINTR (PLRPS,IPLANE,MLS-2,RINTS,IFOLD,          
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                  
        REF = 'PLRPS INTEGRATED OVER R' // COMMAP                               
C                                                                               
C
C       IT IS NECESSARY TO ADJUST THE OUTPUT ARRAY IF THE PLOTTING RANGE 
C       STARTS AT NEGATIVE THETA. THIS IS DONE SO THAT THE LIMITER WHICH
C       IS AT THETA = 0 MAY BE CENTRED ON THE PLOT. 
C
        IF (VMIN.LT.0.0) THEN 
          CALL RTORD (RINTS,VMIN)
        ENDIF
C
C       CONVERT TOUTS COORDINATE AXIS, AFTER ADJUSTMENT, BACK TO DEGREES
C
        CALL CVRTDEG ( TOUTS,NAS,VMIN,VMAX)
C                             
        CALL LIM_DRAW (TOUTS,TWIDS,RINTS(1,-1),MAXNAS,NAS,        
     >    ANLY,MLS,KSMOTH,VMIN,VMAX,0.0,BIGG,IGLS,ITEC,AVS,NAVS,                
     >    JOB,TITLE,THLAB,FLAB,PLABS(1),REF,RTVIEW,PLANE,TABLE,IPLOT,2)         
C                                                                               
      ELSEIF (BREF.EQ.'RPR') THEN                                               
C     ===========================                                               
        CALL RINTT (PLRPS,IPLANE,MLS-2,TINTS,IFOLD,          
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                  
        REF = 'PLRPS INTEGRATED OVER A' // COMMAP                               
        CALL LIM_DRAW (ROUTS,RWIDS,TINTS(1,-1),MAXNRS,NRS,ANLY,                     
     >    MLS,KSMOTH,VMIN,VMAX,0.0,BIGG,IGLS,ITEC,AVS,NAVS,                     
     >    JOB,TITLE,RLAB,FLAB,PLABS(1),REF,RTVIEW,PLANE,TABLE,IPLOT,2)   
C                                                                               
      ELSEIF (BREF.EQ.'RPS') THEN                                               
C     ===========================                                               
        WRITE(6,*) 'S:MLS,MLS-2:',MLS,MLS-2 
        CALL RINTPSI (PLRPS,IPLANE,MLS-2,SINTS,IFOLD,          
     >              SSS,FP,FT,1.0,VMIN,VMAX,SMIN,SMAX)
        REF = 'PLRPS INTEGRATED FROM S'                               
        CALL LIM_DRAW (SOUTS1,SWIDS,SINTS(1,-1),MAXNSS,NSS,ANLY,                 
     >    MLS,KSMOTH,VMIN,VMAX,0.0,BIGG,IGLS,ITEC,AVS,NAVS,                     
     >    JOB,TITLE,SLAB,FLAB,PLABS(1),REF,SVIEW,PLANE,TABLE,IPLOT,2)           
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C  TIZS: POSITIONS AT IONISATION                                                
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (BREF.EQ.'RIA') THEN                                               
C     ===========================                                               
        CALL RINTR (TIZS,IPLANE,MAXIZ,XINTS,IFOLD,            
     >              SSS,FP,FT,1.0,SMIN,SMAX)                               
        REF = 'IONISATION INT''D OVER R' // COMMAP                              
C                                                                               
C
C       IT IS NECESSARY TO ADJUST THE OUTPUT ARRAY IF THE PLOTTING RANGE 
C       STARTS AT NEGATIVE THETA. THIS IS DONE SO THAT THE LIMITER WHICH
C       IS AT THETA = 0 MAY BE CENTRED ON THE PLOT. 
C
        IF (VMIN.LT.0.0) THEN 
          CALL RTORD (RINTS,VMIN)
        ENDIF
C
C       CONVERT TOUTS COORDINATE AXIS, AFTER ADJUSTMENT, BACK TO DEGREES
C
        CALL CVRTDEG ( TOUTS,NAS,VMIN,VMAX)
C                             
        CALL LIM_DRAW (TOUTS,TWIDS,RINTS,MAXNAS,NAS,ANLY,               
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,THLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)        
C                                                                               
      ELSEIF (BREF.EQ.'RIR') THEN                                               
C     ===========================                                               
        CALL RINTT (TIZS,PLANE,MAXIZ,TINTS,IFOLD,            
     >              SSS,FP,FT,1.0,SMIN,SMAX)                                  
        REF = 'IONISATION INT''D OVER A' // COMMAP                              
        CALL LIM_DRAW (ROUTS,RWIDS,TINTS,MAXNRS,NRS,ANLY,                           
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,RLAB,FLAB,ZLABS(-2),REF,RTVIEW,PLANE,TABLE,IPLOT,2)         
C                                                                               
C        DO 1985 IZ = -2, MAXIZ                                                 
C          SUM1(IZ) = 0.0                                                       
C          SUM2(IZ) = 0.0                                                       
C          DO 1987 IR = 1, NRS                                                  
C            SUM1(IZ) = SUM1(IZ) + RWIDS(IR) * TINTS(IR,IZ) * ROUTS(IR)         
C            SUM2(IZ) = SUM2(IZ) + RWIDS(IR) * TINTS(IR,IZ)                     
C 1987     CONTINUE                                                             
C          IF (SUM2(IZ).LE.0.0) SUM2(IZ) = 1.E-20                               
C 1985   CONTINUE                                                               
C        WRITE (6,9010) (ZLABS1(IZ)(5:11),SUM1(IZ)/SUM2(IZ),IZ=-2,MAXIZ)        
C                                                                               
      ELSEIF (BREF.EQ.'RIS') THEN                                               
C     ===========================                                               
        CALL RINTPSI (TIZS,IPLANE,MAXIZ,SINTS,IFOLD,        
     >              SSS,FP,FT,1.0,VMIN,VMAX,SMIN,SMAX)
        REF = 'IMPURITY INTEGRATED FROM S'
        CALL LIM_DRAW (SOUTS1,SWIDS,SINTS,MAXNSS,NSS,ANLY,                       
     >    MAXIZ+3+IALL,ISMOTH,VMIN,VMAX,0.0,BIGG,IGZS(-2),ITEC,AVS,NAVS,        
     >    JOB,TITLE,SLAB,FLAB,ZLABS(-2),REF,SVIEW,PLANE,TABLE,IPLOT,2)          
C                                                                               
C-----------------------------------------------------------------------        
C  Z EFFECTIVE ETC.                                                             
C-----------------------------------------------------------------------        
C                                                                               
C---- VIEWPOINT NOT REQUIRED ... USE "NVIEW" STRING WHEN PLOTTING               
C---- INTEGRAL OVER Y INSTEAD OF "XVIEW/YVIEW" IN CALL TO DRAW,                 
C---- AND CALL INTEG WITH "IVU=0" FOR OFF.                                      
C---- THERE IS NO 3D EQUIVALENT OF ZEFFS: PUT DUMMY ARRAY TIZ3 IN               
C---- ARGUMENT LIST FOR RINTX/Y FOLLOWED BY PLANE (WHICH WILL BE 99).           
C---- NOTE 291: ADD DILUTION TO X PLOTS                                         
C                                                                               
      ELSEIF (BREF.EQ.'RZA') THEN                                               
C     ===========================                                               
        CALL RINTR (ZEFFS,99,3-2,RINTS,IFOLD,                   
     >              SSS,FP,FT,CA-CAW,SMIN,SMAX)                             
        REF = 'ZEFFS AVERAGED OVER R'                                           
C                                                                               
C
C       IT IS NECESSARY TO ADJUST THE OUTPUT ARRAY IF THE PLOTTING RANGE 
C       STARTS AT NEGATIVE THETA. THIS IS DONE SO THAT THE LIMITER WHICH
C       IS AT THETA = 0 MAY BE CENTRED ON THE PLOT. 
C
        IF (VMIN.LT.0.0) THEN 
          CALL RTORD (RINTS,VMIN)
        ENDIF
C
C       CONVERT TOUTS COORDINATE AXIS, AFTER ADJUSTMENT, BACK TO DEGREES
C
        CALL CVRTDEG ( TOUTS,NAS,VMIN,VMAX)
C                             
        CALL LIM_DRAW (TOUTS,TWIDS,RINTS(1,-1),MAXNAS,NAS,        
     >    ANLY,3,99,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                      
     >    JOB,TITLE,THLAB,FLAB,ELABS(6),REF,NVIEW,PLANE,TABLE,IPLOT,2)         
C                                                                               
      ELSEIF (BREF.EQ.'RZR') THEN                                               
C     ===========================                                               
        CALL RINTT (ZEFFS,99,3-2,YINTS,IFOLD,                   
     >              SSS,FP,FT,CTWOL,SMIN,SMAX)                                
        REF = 'ZEFFS AVERAGED OVER A'                                           
C        DO 1950 IX = 1, NXS                                                    
C          YINTS(IX,2) = YINTS(IX,0) / CRNBS(IX,1)                              
C 1950   CONTINUE                                                               
C        ELABS(9) = '    DILUTE'                                                
        CALL LIM_DRAW (ROUTS,RWIDS,TINTS(1,-1),MAXNRS,NRS,ANLY,                     
     >    4,99,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                           
     >    JOB,TITLE,RLAB,FLAB,ELABS(6),REF,NVIEW,PLANE,TABLE,IPLOT,2)           
        ELABS(9) = '    POW TOT'                                                
C                                                                               
      ELSEIF (BREF.EQ.'RZS') THEN                                               
C     ===========================                                               
        CALL RINTPSI (ZEFFS,99,3-2,SINTS,IFOLD,                   
     >              SSS,FP,FT,1.0,VMIN,VMAX,SMIN,SMAX)
        REF = 'ZEFFS AVERAGED FROM S'                                           
C        DO 1950 IX = 1, NXS                                                    
C          YINTS(IX,2) = YINTS(IX,0) / CRNBS(IX,1)                              
C 1950   CONTINUE                                                               
C        ELABS(9) = '    DILUTE'                                                
        CALL LIM_DRAW (SOUTS1,SWIDS,SINTS(1,-1),MAXNSS,NSS,ANLY,               
     >    4,99,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                           
     >    JOB,TITLE,SLAB,FLAB,ELABS(6),REF,SVIEW,PLANE,TABLE,IPLOT,2)           
        ELABS(9) = '    POW TOT'                                                
C                                                                               
      ENDIF                                                                     
      GOTO 1960                                                                 
C     *********
C                                                                               
C#######################################################################        
C       3 DIMENSIONAL PLOTS                                                     
C       ===================                                                     
C       LEAP HERE IF BREF STARTS WITH SOMETHIN OTHER THAN "2"                   
C       ASSUME ALL FOLLOWING GRAPH DATA WILL REFER TO 3D CASES.                 
C       UNTIL A REFERENCE BEGINNING T IS FOUND ...                              
C#######################################################################        
C                                                                               
 2000 CONTINUE                                                                  
      WRITE (7,9202)                                                            
 9202 FORMAT(//2X,'REF   TITLE            X RANGE      Y RANGE   NPTS ',        
     >  'STATE PLANE YFOLD ',/1X,71('-'))                                       
C                                                                               
C---- LOOP TO READ IN 3D DETAILS                                                
C                                                                               
 2100 CONTINUE                                                                  
      CALL RDG3D (GRAPH,XMIN,XMAX,YMIN,YMAX,NPTS,ISTATE,IPLANE,IFOLD,           
     >            JSMOTH,'LINE OF 3D GRAPH DETAILS',IERR)                              
      IF (IERR.NE.0) GOTO 9999                                                  

      BREF = GRAPH(1:3)                                                         
c      write(0,*) '3D:',iplot,trim(bref)

      IF (BREF(1:1).NE.'3') GOTO 2500                                           
      IF (NPTS.LE.0)        GOTO 2100                                           
C                                                                               
      DO 2110 J = 1, 4                                                          
        S(J) = ' '                                                              
 2110 CONTINUE                                                                  
C                                                                               
C---- MAXIMUM NO. OF POINTS ALLOWED IS 192  (GHOST80 LIMITATION)                
C                                                                               
      IF (NPTS.GT.192) THEN                                                     
        NPTS = 192                                                              
        S(1) = '*'                                                              
      ELSEIF (NPTS.LT.0) THEN                                                   
        NPTS = 0                                                                
        S(1) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- IONISATION STATE MUST EXIST                                               
C---- SECONDARY NEUTRALS DO NOT EXIST EXCEPT FOR PLRPS !!!                      
C                                                                               
      IF (ISTATE.LT.-2.OR.ISTATE.GT.NIZS+1.OR.                                  
     >                     (ISTATE.EQ.-2.AND.BREF.NE.'3P ')) THEN               
        ISTATE = 99                                                             
        S(2) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- WRITE INTO "PLANE" STRING AS REQUIRED                                     
C                                                                               
      IF (IPLANE.EQ.99) THEN                                                    
        PLANE = ' '                                                             
        INTEGP = 'INT''D OVER P'                                                
      ELSEIF (IPLANE.GE.-MAXNPS.AND.IPLANE.LE.MAXNPS) THEN                      
        WRITE (PLANE,'(4A)') 'PLOT FOR  ',PRINPS(IPLANE-1),                     
     >    ' < P <',PRINPS(IPLANE)                                               
        INTEGP = ' '                                                            
      ELSE                                                                      
        IPLANE = 99                                                             
        INTEGP = 'INT''D OVER P'                                                
        PLANE = ' '                                                             
        S(3) = '*'                                                              
      ENDIF                                                                     
C                                                                               
      TFTRG = 0
      IF (IFOLD.EQ.-1) THEN
        TFTRG = 1
        IFOLD = 0
      ELSEIF (IFOLD.EQ.-2) THEN
        TFTRG = 1 
        IFOLD = 1
      ELSEIF (IFOLD.NE.0.AND.IFOLD.NE.1) THEN                                       
        IFOLD = 1                                                               
        S(4) = '*'                                                              
      ENDIF                                                                     
C
      if (jsmoth.eq.1) then
         smooth = 'Smoothed using neighbourhood averaging.'
      else
         jsmoth = 0
         smooth = ' '  
      endif  
C                                                                               
C---- ADJUST XMIN,XMAX,YMIN,YMAX IF REQUIRED                                    
C                                                                               
      XMIN = MAX (XMIN, CAW)                                                    
      XMAX = MIN (XMAX, CA)                                                     
      IF (IPLANE.EQ.99) THEN                                                    
        IF (YMIN.LT.0.0) YMIN = MAX (YMIN,-CL)                                  
        IF (YMAX.GT.0.0) YMAX = MIN (YMAX, CL)                                  
      ELSE                                                                      
        IF (YMIN.LT.0.0) YMIN = MAX (YMIN,-YS(NY3D))                            
        IF (YMAX.GT.0.0) YMAX = MIN (YMAX, YS(NY3D))                            
      ENDIF                                                                     
      WRITE (LAB3D2,'(F5.2,A6,F5.2,A1,F6.2,A6,F5.2)')                           
     >       XMIN,' < X <',XMAX,' ',YMIN,' < Y <',YMAX                          
C                                                                               
C---- PRINT MESSAGE TO CHANNEL 7 AS PERMANENT RECORD OF GRAPH PARAMS.           
C                                                                               
      write (6,*) 'start:' ,tftrg
      WRITE (7,'(2X,A21,4F6.2,I6,3(A1,I4),A1)') GRAPH,XMIN,XMAX,YMIN,           
     >  YMAX,NPTS,S(1),ISTATE,S(2),IPLANE,S(3),IFOLD,S(4)                       
      IF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') GOTO 2100                              
      IF (ISTATE.EQ.99)                  GOTO 2100                              
      IF (XMIN.GE.XMAX.OR.YMIN.GE.YMAX)  GOTO 2100                              
C                                                                               
C---- SET ITEM TO APPEAR IN "SYMBOL TABLE" LIST                                 
C---- NOTE IT IS ILINE = IL - 2 WE REQUIRE FOR RINT3D, NOT IL - 3  !!           
C                                                                               
      IF (BREF(2:2).EQ.'P') THEN                                                
        DO 2120 IL = 1, NLS                                                     
          IF (PIZS(IL).EQ.ISTATE) THEN                                          
            ILINE = IL - 2                                                      
            LAB3D1= PLABS(IL)(5:11)                                             
          ENDIF                                                                 
 2120   CONTINUE                                                                
      ELSE                                                                      
        LAB3D1= ZLABS1(ISTATE)(5:11)                                            
      ENDIF                                                                     
      ANLY = ' '                                                                
      IF (CANAL.LT.CA .AND. IPLANE.EQ.99 .AND. BREF.NE.'3I ')                   
     >  WRITE (ANLY,'(''ANALYTIC EXTENSION FOR X >'',F6.3)') CANAL              
C                                                                               
C---- CREATE ARRAYS SURFAS/SUREDG TO BE PLOTTED ...                             
C                                                                             
      IF     (BREF.EQ.'3C ') THEN                                               
C     ===========================                                               
        REF = 'IMPURITY ' // INTEGP                                             
        CALL RINT3D (SDLIMS,SDLIM3,ISTATE,IPLANE,IFOLD,NPTS,                    
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,NIZS,SUREDG,QEDGES,TFTRG)                        
      ELSEIF (BREF.EQ.'3R ') THEN                                               
C     ===========================                                               
        REF = 'POWER LOSS ' // INTEGP                                           
        CALL RINT3D (POWLS ,POWL3 ,ISTATE,IPLANE,IFOLD,NPTS,                    
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,NIZS,SUREDG,QEDGES,TFTRG)                        
      ELSEIF (BREF.EQ.'3L ') THEN                                               
C     ===========================                                               
        REF = 'LINE RADIATION ' // INTEGP                                       
        CALL RINT3D (LINES ,LINE3 ,ISTATE,IPLANE,IFOLD,NPTS,                    
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,NIZS,SUREDG,QEDGES,TFTRG)                        
      ELSEIF (BREF.EQ.'3P ') THEN                                               
C     ===========================                                               
        REF = 'PLRPS ' // INTEGP                                                
        CALL RINT3D (PLRPS ,PLRP3 ,ILINE, IPLANE,IFOLD,NPTS,                    
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,NIZS,SUREDG,QEDGES,TFTRG)                        
      ELSEIF (BREF.EQ.'3I ') THEN                                               
C     ===========================                                               
        REF = 'IONISATION ' // INTEGP                                           
        CALL RINT3D (TIZS,  TIZ3,  ISTATE,IPLANE,IFOLD,NPTS,                    
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,NIZS,SUREDG,QEDGES,TFTRG)                        
      ENDIF                                                                     
C
C     SMOOTH THE VALUES IN THE ARRAY SURFAS IF CALLED FOR
C     
      IF (JSMOTH.EQ.1) THEN 
         CALL RSMO3D(SURFAS,NPTS)
      ENDIF  
C                                                                               
C---- PLOT SURFACE INSIDE USUAL BOX, INCLUDING TITLE, JOBNAME, ETC              
C---- GRTSET IS CALLED WITH FLAG=3 TO MEAN "DON'T DRAW ANY AXES"                
C---- THIS SAVES HAVING TO DUPLICATE MOST OF GRTSET BELOW.  THEN                
C---- LIM_GR3D IS CALLED TO DRAW THE SURFACE INSIDE THE FRAMEWORK.                  
C                                                                               
      CALL LIM_GRTSET (TITLE,REF,LAB3D2,PLANE,JOB,0.0,0.0,0.0,1.0,                  
     >    TABLE,XLAB,FLAB,3,SMOOTH,1,ANLY,1)                                    
      CALL LIM_GR3D (SURFAS,NPTS,LAB3D1,IVEW3D,PROJ3D,IBAS3D,                       
     >           SUREDG,LIMEDG)                                                 
      CALL FRAME                                                                
      smooth = ' ' 
C                                                                               
      GOTO 2100                                                                 

C     *********                                                                 
C                                                                               
C#######################################################################        
C       MESH/CONTOUR PLOTS                                                     
C       ===================                                                     
C       LEAP HERE IF BREF STARTS WITH SOMETHIN OTHER THAN "3"                   
C       ASSUME ALL FOLLOWING GRAPH DATA WILL REFER TO MESH CASES.              
C       UNTIL A REFERENCE BEGINNING T IS FOUND ...                              
C#######################################################################        
C                                                                               
C       THIS SECTION WAS ADDED TO SUPPORT CONTOUR AND/OR MESH
C       PLOTTING OF THE PLASMA CHARACTERISTICS ON THE XY XP YP
C       SETS OF AXES. INTEGRATION TAKES PLACE OVER RANGES
C       XMIN,XMAX YMIN,YMAX PMIN,PMAX AS SPECIFIED IN THE INPUT.
C       THE AXIS TO BE INTEGRATED IS ALSO SPECIFIED. THE RESULT IS
C       THE 2D ARRAY SURFAS (NxM) CONTAINING HEIGHTS ABOVE THE PLANE.
C       IF NECESSARY THE INTEGRATED RESULT IS THEN EXPANDED ONTO A
C       LARGER (EVENLY SPACED) GRID. (THIS GRID IS NECESSARY FOR 
C       3D MESH PLOTS. CONTOUR PLOTS ON THE OTHER HAND ONLY REQUIRE
C       THE CO-ORDINATES OF THE IRREGULAR GID SPACING.)
C
C       THESE PLOT OPTIONS ONLY WORK FOR 3D CASES.
C
C       THESE PLOT OPTIONS HAVE NOT BEEN MODIFIED TO WORK  
C       CORRECTLY WITH GHOST
C       
C       DAVID ELDER , 1990 MARCH 14
C
C       *******   

 2500 CONTINUE                                                                  
      WRITE (7,9233)                                                            
 9233 FORMAT(//2X,'REF   TITLE        X RANGE     Y RANGE    P RANGE ',        
     >  'NPTS MPTS STATE YFOLD ',/1X,71('-'))                              

C ======================================================================        
C---- LOOP TO READ IN MESH/CONTOUR DETAILS                                    
C ======================================================================        

 2600 CONTINUE                                                                  
      CALL RDGM (GRAPH,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NPTS,MPTS,ISTATE,
     >          IPLANE,IPLOT,IFOLD, 'LINE OF MESH GRAPH DETAILS',IERR)         
      IF (IERR.NE.0) GOTO 9999                                                  

      BREF = GRAPH(1:3)                                                         
c      write(0,*) 'MESH/CONTOUR:',iplot,trim(bref)
      IF (BREF(1:1).NE.'M') GOTO 3000                                           
      IF ((IPLOT.LE.0).OR.(IPLOT.GT.2))   GOTO 2600             
      IF ((IPLANE.LT.0).OR.(IPLANE.GT.2)) GOTO 2600
      IF (NY3D.LE.1)        GOTO 2600                              
C
C                                                                               
      DO 2610 J = 1, 5                                                          
        S(J) = ' '                                                              
 2610 CONTINUE                                                                  

c     WRITE(0,*) '   Graph: ',npts,mpts,IPLANE,GRAPH

C                                                                               
C---- MAXIMUM NO. OF POINTS ALLOWED IS 192  (GHOST80 LIMITATION)                
C
C---- MAXIMUM ELEMENTS IN MATLAB (NEW ROUTINE CURREMTLY USED FOR OUTPUT)
C     IS 818 THEREFORE LIMIT N,M IN THESE ROUTINES ACCORDINGLY.
C     SAY 90
C                                                                               
      IF (NPTS.GT.90) THEN                                                     
        NPTS = 192                                                              
        S(1) = '*'                                                              
      ELSEIF (NPTS.LT.0) THEN                                                   
        NPTS = 0                                                                
        S(1) = '*'                                                              
      ENDIF                                                                     
      IF (MPTS.GT.90) THEN                                                     
        MPTS = 192                                                              
        S(5) = '*'                                                              
      ELSEIF (MPTS.LT.0) THEN                                                   
        MPTS = 0                                                                
        S(5) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- IONISATION STATE MUST EXIST                                               
C---- SECONDARY NEUTRALS DO NOT EXIST EXCEPT FOR PLRPS !!!                      
C                                                                               
      IF (ISTATE.LT.-2.OR.ISTATE.GT.NIZS+1.OR.                                  
     >                     (ISTATE.EQ.-2.AND.BREF.NE.'MP ')) THEN               
        ISTATE = 99                                                             
        S(2) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- ADJUST XMIN,XMAX,YMIN,YMAX,PMIN,PMAX IF REQUIRED                      
C                                                                               
      XMIN = MAX (XMIN, CAW)                                                    
      XMAX = MIN (XMAX, CA)                                                    
      WRITE(6,*) 'NY3D:',NY3D,YMIN,YMAX,-YS(NY3D),YS(NY3D)
      IF (YMIN.LT.0.0) YMIN = MAX (YMIN,-YS(NY3D))                            
      IF (YMAX.GT.0.0) YMAX = MIN (YMAX, YS(NY3D))                            

C
C     PMIN AND PMAX DO NOT NEED ADJUSTING FOR NOW
C     SINCE THE DEFAULT VALUES ARE INCREDIBLY LARGE
C     I MAY WIND UP USING THE GIVEN VALUES TO ESTIMATE 
C     A SIZE FOR THE OUTER BINS ON INTEGRATION      
C
C                                                                               
C---- WRITE INTO "PLANE" STRING AS REQUIRED                                     
C                                                                               
      plane = ' '
      IF (IPLANE.EQ.0) THEN                                                    
        integp = 'INT''D OVER X'
        write(plane,'(a,2(1x,f8.3))')
     >    integp,xmin,xmax                                                
        write(plane,'(a,2(1x,f8.3))') integp,
     >        xmin,xmax
      ELSEIF (IPLANE.EQ.1) THEN                      
        integp = 'INT''D OVER Y'
        write(plane,'(a,2(1x,f8.3))')
     >     integp,ymin,ymax                                                
        write(plane,'(a,2(1x,f8.3))') integp,
     >        ymin,ymax
      ELSEIF (IPLANE.EQ.2) THEN                                            
        integp = 'INT''D OVER P'
        write(plane,'(a,2(1x,f8.3))')
     >    integp,pmin,pmax                                                
        write(plane,'(a,2(1x,f8.3))') integp,
     >        pmin,pmax
      ENDIF                                                                     
C                                                                               
      IF (IFOLD.NE.0.AND.IFOLD.NE.1) THEN                                       
        IFOLD = 1                                                               
        S(4) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- PRINT MESSAGE TO CHANNEL 7 AS PERMANENT RECORD OF GRAPH PARAMS.           
C                                                                               
      WRITE (7,'(2X,A21,6F6.2,2I6,3(A1,I4),A1)') GRAPH,XMIN,XMAX,YMIN,        
     >  YMAX,PMIN,PMAX,NPTS,MPTS,S(1),ISTATE,S(2),IPLANE,S(3),
     >  IFOLD,S(4)                       
      IF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') GOTO 2600                              
      IF (ISTATE.EQ.99)                  GOTO 2600                              
      IF (XMIN.GE.XMAX.OR.YMIN.GE.YMAX.OR.PMIN.GE.PMAX)  GOTO 2600            
C                                                                               
C---- SET ITEM TO APPEAR IN "SYMBOL TABLE" LIST                                 
C---- NOTE IT IS ILINE = IL - 2 WE REQUIRE FOR RINTM, NOT IL - 3  !!           
C---- THE REASON FOR THIS IS THAT THE DUMMY ARRAY IN RINTM IS DECLARED
C---- FROM -1:NIZS (OR NLS) NOT 1:NIZS (OR NLS) THUS THE INDEX
C---- NUMBER FOR A PARTICULAR LINE MUST BE OFFSET BY TWO.
C
C---- SIMILARLY, IN ORDER TO MAKE THE ADDRESS CALCULATIONS WORK
C---- THE UPPER BOUND MUST BE NLS-2, THUS INLS.
C 
      IF (BREF(2:2).EQ.'P') THEN                                                
        DO 2620 IL = 1, NLS                                                     
          IF (PIZS(IL).EQ.ISTATE) THEN                                          
            ILINE = IL -2                                                      
            LAB3D1= PLABS(IL)(5:11)                                             
          ENDIF                                                                 
 2620   CONTINUE                                                                
        INLS = NLS-2 
      ELSE                                                                      
        LAB3D1= ZLABS1(ISTATE)(5:11)                                            
      ENDIF                                                                     
      ANLY = ' '                                                                
C     IF (CANAL.LT.CA .AND. IPLANE.EQ.99 .AND. BREF.NE.'3I ')                   
C    >  WRITE (ANLY,'(''ANALYTIC EXTENSION FOR X >'',F6.3)') CANAL              
C                                                                               
C---- CREATE ARRAYS SURFAS/SUREDG TO BE PLOTTED ...                             
C                                                                               
      IF     (BREF.EQ.'MC ') THEN                                               
C     ===========================                                               
        REF = 'IMPURITY ' // INTEGP                                             
c slmod begin
c        write(0,*) 'MC:',istate

        IF (BIG) THEN

           ! jdemod - don't understand this since it prevents integration over charge states?
           !IF (ISTATE.GT.NIZS) GOTO 2600
           IF (ISTATE.GT.NIZS+1) GOTO 2600
c
c Find the maximum density and the 10% value of the lowest peak value:
c
          CPMIN =  1.0E30
          CPMAX = -1.0E30

! jdemod - comment out for now ... a lot of overhead just to calculate a common min value
!     if (iz.lt.nizs+1) then


          ! jdemod - I think the point of this code is to use a common minimum value and different
          ! maximum across the series of cloud contour plots
          
!          DO IZ = 0, NIZS
!
!            CALL RINTM (SDLIM3,IZ,IPLANE,IFOLD,NPTS,MPTS,
!     >        SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
!     >        COORD1,COORD2,POUTS,MAXIZS)        
!
!            TMAX = -1.0E30
!
!             DO IX=1,NPTS
!               DO IY=1,MPTS
!                 TMAX = MAX(SURFAS(IX,IY),TMAX)
!               ENDDO
!             ENDDO
!
!             CPMIN = MIN(0.05*TMAX,CPMIN)
!             CPMAX = MAX(     TMAX,CPMAX)
!
!c            IF (ISTATE.EQ.1) WRITE(0,'(A,I6,3G12.4)') 
!c     +        'DEN  ',IZ,TMAX,CPMIN,CPMAX
!           
!             IF (IZ.EQ.ISTATE) TSMAX = TMAX
!
!           ENDDO
!
!           CPMAX = TSMAX
!
!          CALL RINTM (SDLIM3,ISTATE,IPLANE,IFOLD,NPTS,MPTS,
!     >      SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
!     >      COORD1,COORD2,POUTS,MAXIZS)        
!
!

!        else
           ! iz = nizs+1
           CALL RINTM (SDLIM3,ISTATE,IPLANE,IFOLD,NPTS,MPTS,
     >      SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
     >      COORD1,COORD2,POUTS,MAXIZS)        

          
            TMIN = 1.0e30
            TMAX = -1.0E30

            DO IX=1,NPTS
              DO IY=1,MPTS
                IF (SURFAS(IX,IY).GT.0.0) TMIN = MIN(TMIN,SURFAS(IX,IY))
                TMAX = MAX(SURFAS(IX,IY),TMAX)
              ENDDO
            ENDDO

            CPMIN = MIN(0.99*TMIN,CPMIN)
            CPMAX = MAX(     TMAX,CPMAX)
c            write(0,*) 'PLT MAX,MIN:',cpmax,cpmin

c         endif





       ELSE
          WRITE(0,*) 'ERROR: Cross plot contours not available in 2D.'
          STOP
        ENDIF
c
c        CALL RINTM (SDLIM3,ISTATE,IPLANE,IFOLD,NPTS,MPTS,
c     >    SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
c     >    COORD1,COORD2,POUTS,MAXIZS)        
c slmod end
      ELSEIF (BREF.EQ.'MR ') THEN                                               
C     ===========================                                               
        REF = 'POWER LOSS ' // INTEGP                                           
        CALL RINTM (POWL3 ,ISTATE,IPLANE,IFOLD,NPTS,MPTS,                
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
     >    COORD1,COORD2,POUTS,MAXIZS)            
      ELSEIF (BREF.EQ.'ML ') THEN                                               
C     ===========================                                               
        REF = 'LINE RADIATION ' // INTEGP                                       
        CALL RINTM (LINE3 ,ISTATE,IPLANE,IFOLD,NPTS,MPTS,              
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
     >    COORD1,COORD2,POUTS,MAXIZS)       
       ELSEIF (BREF.EQ.'MP ') THEN                                               
C     ===========================                                               
        REF = 'PLRPS ' // INTEGP                                                
c slmod begin
        CALL RINTM (PLRP3 ,ILINE, IPLANE,IFOLD,NPTS,MPTS,            
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,INLS,IPLOT,
     >    COORD1,COORD2,POUTS,MAXNLS-2)             
      
        CPMIN =  1.0E30
        CPMAX = -1.0E30
      
        TMAX = -1.0E30
      
        DO IX=1,NPTS
          DO IY=1,MPTS
            TMAX = MAX(SURFAS(IX,IY),TMAX)
          ENDDO
        ENDDO
      
        CPMIN = MIN(0.05*TMAX,CPMIN)
        CPMAX = MAX(     TMAX,CPMAX)
      
c        WRITE(0,'(A,I6,3G12.4)') 
c     +    'PLRP ',IZ,TMAX,CPMIN,CPMAX
c
c        CALL RINTM (PLRP3 ,ILINE, IPLANE,IFOLD,NPTS,MPTS,            
c     >    SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,INLS,IPLOT,
c     >    COORD1,COORD2,POUTS,MAXNLS-2)             
c slmod end
      ELSEIF (BREF.EQ.'MI ') THEN                                               
C     ===========================                                               
        REF = 'IONISATION ' // INTEGP                                           
        CALL RINTM (TIZ3,  ISTATE,IPLANE,IFOLD,NPTS,MPTS,               
     >    SURFAS,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
     >    COORD1,COORD2,POUTS,MAXIZS)             
      ENDIF                                                                     
C                                                                               
C---- PLOT SURFACE INSIDE USUAL BOX, INCLUDING TITLE, JOBNAME, ETC              
C---- GRTSET IS CALLED WITH FLAG=3 TO MEAN "DON'T DRAW ANY AXES"                
C---- THIS SAVES HAVING TO DUPLICATE MOST OF GRTSET BELOW.  THEN                
C---- LIM_GRM IS CALLED TO DRAW THE SURFACE INSIDE THE FRAMEWORK.                  
C      
      IF (IPLANE.EQ.0) THEN 
         FLAB = 'Y (M)'  
         XLAB = 'P (M)'
      ELSEIF (IPLANE.EQ.1) THEN
         FLAB = 'X (M)'
         XLAB = 'P (M)'
      ELSEIF (IPLANE.EQ.2) THEN 
         FLAB = 'X (M)'
         XLAB = 'Y (M)'
      ENDIF 
C
C
      IX = 3
      IF (IPLOT.EQ.1) IX=2
C
C     DUMMY AXIS MIN,MAX ARE SENT TO GRTSET HERE BECAUSE THEY ARE 
C     USUALLY NOT NEEDED. CAN BE ADDED LATER IF NECESSARY.
C

c slmod begin
      IF (BREF.EQ.'MC '.OR.BREF.EQ.'ML '.OR.BREF.EQ.'MP ') THEN

c        TMIN= 1E30
c        TMAX=-1E30

c        DO IX=1,NPTS
c          DO IY=1,MPTS
c            IF (SURFAS(IX,IY).GT.0.0) TMIN = MIN(TMIN,SURFAS(IX,IY))
c            TMAX = MAX(SURFAS(IX,IY),TMAX)
c          ENDDO
c        ENDDO
c
c This is the correct use for TMIN...?
c 
c        TMIN = 0.0
c
c  Set up labels:
c
        IF (IPLANE.EQ.0) THEN
          FLAB = 'P (M)'  
          XLAB = 'Y (M)'

          CALL LIM_GRTSET (TITLE,REF,LAB3D2,PLANE,JOB,YMIN,YMAX,
     >                PMIN,PMAX,              
     >                TABLE,XLAB,FLAB,IX,SMOOTH,1,ANLY,1)   
        ELSEIF (IPLANE.EQ.1) THEN
          FLAB = 'P (M)'  
          XLAB = 'X (M)'

          CALL LIM_GRTSET (TITLE,REF,LAB3D2,PLANE,JOB,XMIN,XMAX,
     >                PMIN,PMAX,              
     >                TABLE,XLAB,FLAB,IX,SMOOTH,1,ANLY,1)   
        ELSEIF (IPLANE.EQ.2) THEN
          FLAB = 'Y (M)'  
          XLAB = 'X (M)'

          CALL LIM_GRTSET (TITLE,REF,LAB3D2,PLANE,JOB,XMIN,XMAX,
     >                 YMIN,YMAX,              
     >                 TABLE,XLAB,FLAB,IX,SMOOTH,1,ANLY,1)   
         ENDIF

        CALL PSPACE(0.0,1.35,0.0,1.0)
        CALL MAP   (0.0,1.35,0.0,1.0)
        CALL CTRMAG(10)  
c
c PLRP integration:
c
        IF (BREF.EQ.'MP ') THEN

c          WRITE(0,*) 'PLRP integration'

          NUMSUM = 0.0
          VOLSUM = 0.0

          DO IY = -NYS, NYS

            IF (IY.NE.0) THEN

              DO IP = -MAXNPS+1, MAXNPS-1
                DO IX = 1, NXS
                  TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                  NUMSUM = NUMSUM + PLRP3(IX,IY,ILINE+2,IP) * TMPVOL
                  VOLSUM = VOLSUM + TMPVOL
                ENDDO 
              ENDDO

            ENDIF
          ENDDO           
          
          WRITE(DUM,'(F5.1,A,G9.3,A)') 
     +      PLAMS(ILINE+2),' nm, PLRP =',NUMSUM,' photons/s'

          CALL PCSEND (1.33 ,0.50,DUM(1:44))                 
        ENDIF 
c
c Output case parameters:
c
        IF (IPLANE.EQ.0) THEN

          WRITE(DUM,'(A5,F7.1,A3)') 
     +      'TEB =',CTBIN,' eV'
c          WRITE(DUM,'(A5,F7.1,A9,G9.4,A4)') 
c     +      'TEB =',CTBIN,' exp( x /',CLTIN1,') eV'
          CALL PCSEND (1.33 ,0.46,DUM(1:44))                 
        
          WRITE(DUM,'(A5,F7.1,A3)') 
     +      'TIB =',CTIBIN,' eV'
c          WRITE(DUM,'(A5,F7.1,A9,G9.4,A4)') 
c     +      'TIB =',CTIBIN,' exp( x /',CLTIIN1,') eV'
          CALL PCSEND (1.33 ,0.44,DUM(1:44))                 
        
          WRITE(DUM,'(A5,G7.2,A5)') 
     +      'NB  =',CNBIN,' m^-3'
c          WRITE(DUM,'(A5,G8.2,A9,G9.4,A4)') 
c     +      'NB  =',CNBIN,' exp( x /',CLNIN1,') m3'
          CALL PCSEND (1.33 ,0.41,DUM(1:44))                 
        
          WRITE(DUM,'(A24,G12.5,A5)') 
     +      'ELECTRIC FIELD     ',CEYIN  ,' V/m ' 
          CALL PCSEND (1.235,0.38,DUM(1:38))        
          WRITE(DUM,'(A24,G12.5,A5)') 
     +      'PARALLEL FLOW VEL  ',CVHYIN ,' m/s ' 
          CALL PCSEND (1.235,0.36,DUM(1:38))        
        
c          WRITE(DUM,'(A24,G12.5,A5)')
c     +      'POLOIDAL FLOW VEL  ',CVPOL  ,' m/s ' 
c          CALL PCSEND (1.235,0.33,DUM(1:38))        
c          WRITE(DUM,'(A24,G12.5,A5)')
c     +      'DPOL               ',CDPOL  ,' m2/s' 
c          CALL PCSEND (1.235,0.31,DUM(1:38))        
        ENDIF

        WRITE(DUM,'(A24,I9  ,A4)') 'IONIZATION STATE   ',ISTATE ,'    ' 
        CALL PCSEND (1.235,0.28,DUM(1:38))        

!        CLEVLS(1)  = CPMIN
!        CLEVLS(2)  = 0.10*(cpMAX-cpmin) + cpmin
!        CLEVLS(3)  = 0.20*(cpMAX-cpmin) + cpmin
!        CLEVLS(4)  = 0.30*(cpMAX-cpmin) + cpmin
!        CLEVLS(5)  = 0.40*(cpMAX-cpmin) + cpmin
!        CLEVLS(6)  = 0.50*(cpMAX-cpmin) + cpmin 
!        CLEVLS(7)  = 0.60*(cpMAX-cpmin) + cpmin
!        CLEVLS(8)  = 0.70*(cpMAX-cpmin) + cpmin
!        CLEVLS(9)  = 0.80*(cpMAX-cpmin) + cpmin
!        CLEVLS(10) = 0.90*(cpMAX-cpmin) + cpmin

        CLEVLS(1)  = CPMIN
        CLEVLS(2)  = 0.01*(cpMAX-cpmin) + cpmin
        CLEVLS(3)  = 0.05*(cpMAX-cpmin) + cpmin
        CLEVLS(4)  = 0.10*(cpMAX-cpmin) + cpmin
        CLEVLS(5)  = 0.20*(cpMAX-cpmin) + cpmin
        CLEVLS(6)  = 0.30*(cpMAX-cpmin) + cpmin 
        CLEVLS(7)  = 0.40*(cpMAX-cpmin) + cpmin
        CLEVLS(8)  = 0.50*(cpMAX-cpmin) + cpmin
        CLEVLS(9)  = 0.70*(cpMAX-cpmin) + cpmin
        CLEVLS(10) = 0.90*(cpMAX-cpmin) + cpmin
        

!        CLEVLS(1)  = CPMIN
!        CLEVLS(2)  = 0.10*CPMAX
!        CLEVLS(3)  = 0.20*CPMAX
!        CLEVLS(4)  = 0.40*CPMAX
!        CLEVLS(5)  = 0.60*CPMAX
!        CLEVLS(6)  = 0.70*CPMAX 
!        CLEVLS(7)  = 0.80*CPMAX
!        CLEVLS(8)  = 0.90*CPMAX
!        CLEVLS(9)  = 0.95*CPMAX
!        CLEVLS(10) = 0.99*CPMAX
        
c        CLEVLS(1)  = 0.05*(TMAX-TMIN) + TMIN
c        CLEVLS(2)  = 0.10*(TMAX-TMIN) + TMIN
c        CLEVLS(3)  = 0.20*(TMAX-TMIN) + TMIN
c        CLEVLS(4)  = 0.30*(TMAX-TMIN) + TMIN
c        CLEVLS(5)  = 0.40*(TMAX-TMIN) + TMIN
c        CLEVLS(6)  = 0.50*(TMAX-TMIN) + TMIN 
c        CLEVLS(7)  = 0.60*(TMAX-TMIN) + TMIN
c        CLEVLS(8)  = 0.70*(TMAX-TMIN) + TMIN
c        CLEVLS(9)  = 0.80*(TMAX-TMIN) + TMIN
c        CLEVLS(10) = 0.90*(TMAX-TMIN) + TMIN

c        WRITE(0,*) 'Calling LIM_GRCONT95'
 
c        DO IX = 1, NPTS
c          DO IY = 1, MPTS 
c            WRITE(6,'(A,2I3,3G10.3)') 'density contour: ',IX,IY,
c     +        SURFAS(IX,IY),SDLIM3(IX,IY,1,0),SDLIM3(IX,-IY,1,0)
c          ENDDO
c        ENDDO
 
        DO I=1,10
          WRITE (NAMES,'(4X,1P,E8.1)') CLEVLS(I)
          CALL LIM_GRCONT95(SURFAS,1,NPTS,192,1,MPTS,192,CLEVLS(I),
     +                  COORD1,COORD2,NAMES)
        ENDDO
      ELSE
c slmod end
        CALL LIM_GRTSET (TITLE,REF,LAB3D2,PLANE,JOB,0.0,1.0,0.0,1.0,              
     >       TABLE,XLAB,FLAB,IX,SMOOTH,1,ANLY,1)                             
        CALL LIM_GRM (SURFAS,NPTS,MPTS,IPLOT,IPLANE,LAB3D1,IVEW3D,
     >       PROJ3D,IBAS3D,COORD1,COORD2)
c slmod
      ENDIF
c slmod end
      CALL FRAME                                                                
C                                                                               
      GOTO 2600                                                                 
C     *********                                                                 
C                                                                               
C#######################################################################        
C    TIME DEPENDENT PLOTS                                                       
C    ====================                                                       
C  SOMETIMES A PARAMETER WILL BE CHANGED - MARK IT WITH A '*' IN OUTPUT         
C  LIST FOR POSTERITY.                                                          
C-----------------------------------------------------------------------        
C                                                                               
 3000 CONTINUE                                                                  
      WRITE (7,9203)                                                            
 9203 FORMAT(//2X,'REF   TITLE        PLOTTING RANGE PLOT SMOOTH MAXT ',        
     >  'PLANE YFLD STATE VU',/1X,71('-'))                                      
C                                                                               
 3100 CONTINUE                                                                  
      CALL RDGT(GRAPH,VMIN,VMAX,IPLOT,JSMOTH,MAXIT,IPLANE,IFOLD,ISTATE,         
     >          IVU,'LINE OF GRAPH DETAILS',IERR)                               
      IF (IERR.NE.0) GOTO 9999                                                  

      BREF = GRAPH(1:3)                                                         
c      write(0,*) 'TIME DEP:',iplot,trim(bref)

      IF (BREF(1:1).NE.'T'.and.bref(1:3).ne.'000') GOTO 4000                                           
      IF (IPLOT.EQ.0)       GOTO 3100                                           
C                                                                               
      DO 3110 J = 1, 7                                                          
        S(J) = ' '                                                              
 3110 CONTINUE                                                                  
C                                                                               
C---- P PLOTS NOT AVAILABLE IF DPOL=0 - NO POINT TRYING TO PLOT THEM            
C---- PLOTS OF PLRPS SOMETIMES NOT AVAILABLE                                    
C                                                                               
      IF     (IPLOT.NE.0.AND.BREF(3:3).EQ.'P'.AND.CDPOL.LE.0.0) THEN            
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ELSEIF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') THEN                               
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ELSEIF (IPLOT.NE.0.AND.IMODE.EQ.2) THEN                                   
        IPLOT = 0                                                               
        S(1) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- NOT ALLOWED TO SMOOTH P GRAPHS ...                                        
C                                                                               
      IF (JSMOTH.NE.99.AND.BREF(3:3).EQ.'P') THEN                               
        JSMOTH = 99                                                             
        S(2) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- CAN'T PLOT ABOVE MAX RECORDED TIME POINT                                  
C                                                                               
      IF (MAXIT.GT.NTS) THEN                                                    
        MAXIT = NTS                                                             
        S(3) = '*'                                                              
      ELSEIF (MAXIT.LT.1) THEN                                                  
        MAXIT = 1                                                               
        S(3) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- 3D RESULTS NOT CALCULATED WHEN DPOL=0                                     
C                                                                               
      IF (IPLANE.NE.99.AND.(BREF(3:3).EQ.'P'.OR.CDPOL.LE.0.0)) THEN             
        IPLANE = 99                                                             
        S(4) = '*'                                                              
      ENDIF                                                                     
      IF (IPLANE.EQ.99) THEN                                                    
        PLANE = ' '                                                             
        COMMAP = ',P'                                                           
      ELSEIF (IPLANE.GE.-MAXNPS.AND.IPLANE.LE.MAXNPS) THEN                      
        WRITE (PLANE,'(4A)') 'PLOT FOR  ',PRINPS(IPLANE-1),                     
     >    ' < P <',PRINPS(IPLANE)                                               
        COMMAP = ' '                                                            
      ELSE                                                                      
        IPLANE = 99                                                             
        PLANE = ' '                                                             
        S(4) = '*'                                                              
        COMMAP = ',P'                                                           
      ENDIF                                                                     
      ANLY = ' '                                                                
C                                                                               
C---- SEPARATE RESULTS FOR -Y, +Y NOT ALWAYS AVAILABLE                          
C---- FOLDING N/A FOR X GRAPHS, P GRAPHS                                        
C                                                                               
      IF (IFOLD.NE.0.AND.(BREF(3:3).EQ.'X'.OR.BREF(3:3).EQ.'P')) THEN           
        IFOLD = 0                                                               
        S(5) = '*'                                                              
      ENDIF                                                                     
      IF (IFOLD.NE.0.AND.IFOLD.NE.1) THEN                                       
        IFOLD = 1                                                               
        S(5) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- IONISATION STATE MUST EXIST                                               
C                                                                               
      IF (ISTATE.LT.-2.OR.ISTATE.GT.NIZS.OR.                                    
     >                   (ISTATE.EQ.-2.AND.BREF(2:2).NE.'P')) THEN              
        ISTATE = 99                                                             
        S(6) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- "VIEW POINT FLAG" - SET UP VIEWING POSITION ...                           
C                                                                               
      IF (IVU.NE.0 .AND.(BREF(3:3).EQ.'P'.OR.IFOLD.NE.0)) THEN                  
        IVU = 0                                                                 
        S(7) = '*'                                                              
      ENDIF                                                                     
C                                                                               
      XVIEW = ' '                                                               
      YVIEW = ' '                                                               
      YLAB = YLAB1                                                              
      DO 3111 IY = -MAXNYS, MAXNYS                                              
        YOUTS(IY) = YOUTS1(IY)                                                  
        IF (IY.EQ.0) THEN                                                       
          YWIDSS(IY) = 0.0                                                      
        ELSE                                                                    
          YWIDSS(IY)= YWIDS(IABS(IY))                                           
        ENDIF                                                                   
 3111 CONTINUE                                                                  
C                                                                               
      IF (IVU.EQ.1) THEN                                                        
        RV = RV1                                                                
        XV = XV1                                                                
        YV = YV1                                                                
        WRITE (XVIEW,9001) RV,XV,YV                                             
        WRITE (YVIEW,9001) RV,XV,YV                                             
      ELSEIF (IVU.GE.2) THEN                                                    
        RV = RV2                                                                
        XV = XV2                                                                
        YV = YV2                                                                
        WRITE (XVIEW,9001) RV,XV,YV                                             
        WRITE (YVIEW,9001) RV,XV,YV                                             
        YLAB = YLAB2                                                            
        IF (GRAPH(20:20).EQ.'Y') GRAPH(20:20) = 'T'                             
      ENDIF                                                                     
C                                                                               
      YLIMIT = MIN (RV/3.0, YLIM)                                               
      NYSLIM = MIN (IPOS(YLIMIT,YS,NYS-1), NYS/2)                               
      WRITE (6,9006) RV,XV,YV,YLIMIT                                            
C                                                                               
C---- ADJUST VMIN,VMAX IF REQUIRED                                              
C                                                                               
      IF     (BREF(3:3).EQ.'X') THEN                                            
        VMIN = MAX (VMIN,CAW)                                                   
        VMAX = MIN (VMAX, CA)                                                   
      ELSEIF (BREF(3:3).EQ.'Y') THEN                                            
        IF (VMIN.LT.0.0) VMIN = MAX (VMIN,-YS(NY3D))                            
        IF (VMAX.GT.0.0) VMAX = MIN (VMAX, YS(NY3D))                            
      ELSEIF (BREF(3:3).EQ.'P') THEN                                            
        VMIN = MAX (VMIN, POUTS(1-MAXNPS))                                      
        VMAX = MIN (VMAX, POUTS(MAXNPS-1))                                      
      ENDIF                                                                     
C                                                                               
C---- WRITE OUT (POSSIBLY ADJUSTED) LINE OF GRAPH DETAILS                       
C                                                                               
      WRITE (7,'(2X,A21,2F6.2,I4,A1,I6,5(A1,I4),A1)')                           
     >  GRAPH,VMIN,VMAX,IPLOT,                                                  
     >  S(1),JSMOTH,S(2),MAXIT,S(3),IPLANE,S(4),IFOLD,S(5),ISTATE,S(6),         
     >  IVU,S(7)                                                                
      IF (NLS.EQ.0.AND.BREF(2:2).EQ.'P') GOTO 3100                              
      IF (VMIN.GE.VMAX)                  GOTO 3100                              
C                                                                               
C---- SET ITEM TO APPEAR IN "SYMBOL TABLE" LIST                                 
C                                                                               
      IF (BREF(2:2).EQ.'P') THEN                                                
        SSS = .TRUE.                                                            
        DO 3120 IL = 1, NLS                                                     
          IF (PIZS(IL).EQ.ISTATE) THEN                                          
            ILINE = IL - 2                                                      
            TABLE = PLABS(IL)(5:11)                                             
          ENDIF                                                                 
 3120   CONTINUE                                                                
      ELSE                                                                      
        SSS = .FALSE.                                                           
        TABLE = ZLABS1(ISTATE)(5:11)                                            
      ENDIF                                                                     
C                                                                               
C---- SET TIMES ARRAY TO CORRECT TIMES FOR THIS IONISATION STATE                
C                                                                               
      DO 3130 IT = 1, MAXIT                                                     
       IF (ISTATE.LE.0) THEN                                                    
        WRITE (TLABS(IT),'(3X,1P,E8.1,A)') DWELTS(0)*DWELFS(IT),'S'             
       ELSE                                                                     
        WRITE (TLABS(IT),'(3X,1P,E8.1,A)')DWELTS(ISTATE)*DWELFS(IT),'S'         
       ENDIF                                                                    
 3130 CONTINUE                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C  TIME DEPENDENT CLOUDS                                                        
C-----------------------------------------------------------------------        
C                                                                               
      IF     (BREF.EQ.'TCY') THEN                                               
C     ==========================                                                
        CALL TINTX (LIM5,ISTATE,IPLANE,MAXIT,XJNTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM)                                                 
        CALL PROJEC (MAXIZ,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,CSINTB)        
        REF = 'IMPURITY INTEGRATED OVER X' // COMMAP                            
        CALL LIM_DRAW  (YOUTS,YWIDSS,XJNTS,2*MAXNYS+1,2*MAXNYS+1,ANLY,              
     >    MAXIT,JSMOTH,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                   
     >    JOB,TITLE,YLAB,FLAB,TLABS,REF,YVIEW,PLANE,TABLE,IPLOT,2)              
C                                                                               
      ELSEIF (BREF.EQ.'TCX') THEN                                               
C     ==========================                                                
        CALL TINTY (LIM5,ISTATE,IPLANE,MAXIT,YJNTS,IFOLD,RV,XV,YV,IVU,          
     >              SSS,NYSLIM)                                                 
        REF = 'IMPURITY INTEGRATED OVER Y' // COMMAP                            
        CALL LIM_DRAW  (XOUTS,XWIDS,YJNTS,MAXNXS,NXS,ANLY,                          
     >    MAXIT,JSMOTH,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                   
     >    JOB,TITLE,XLAB,FLAB,TLABS,REF,XVIEW,PLANE,TABLE,IPLOT,2)              
C                                                                               
      ELSEIF (BREF.EQ.'TCP') THEN                                               
C     ==========================                                                
        CALL TINTXY (LIM5,ISTATE,IPLANE,MAXIT,XYJNTS,IFOLD,RV,XV,YV,0,          
     >               SSS,NYSLIM)                                                
        REF = 'IMPURITY INTEGRATED OVER X,Y'                                    
        CALL LIM_DRAW  (POUTS,PWIDS,XYJNTS,2*MAXNPS+1,2*MAXNPS+1,ANLY,              
     >    MAXIT,JSMOTH,VMIN,VMAX,0.0,BIGG,IGTS,ITEC,AVS,NAVS,                   
     >    JOB,TITLE,PLAB,FLAB,TLABS,REF,NVIEW,PLANE,TABLE,IPLOT,2)              
C                                                                               
      ENDIF                                                                     
      GOTO 3100                                                                 
C     *********                                                                 
C                                                                               
C#######################################################################        
C       CONTOUR PLOTS                                                           
C       =============                                                           
C       LEAP HERE IF BREF STARTS WITH SOMETHING OTHER THAN "T"                  
C       ASSUME ALL FOLLOWING GRAPH DATA WILL REFER TO CONTOUR PLOTS             
C       UNTIL A REFERENCE BEGINNING OTHER THAN C IS FOUND ...                   
C#######################################################################        
C                                                                               
 4000 CONTINUE                                                                  
      WRITE (7,9204)                                                            
 9204 FORMAT(//2X,'REF   TITLE            X RANGE      Y RANGE   PLOT ',        
     >  'SMOOTH NC DUMMY',/1X,71('-'))                                          
C                                                                               
C---- LOOP TO READ IN 3D DETAILS                                                
C                                                                               
 4100 CONTINUE                                                                  
      CALL RDG3D (GRAPH,XMIN,XMAX,YMIN,YMAX,IPLOT,ISMOTH,NCONT,CTIZ,            
     >            IDUM2,'LINE OF CONTOUR PLOTS DETAILS',IERR)                         
      IF (IERR.NE.0) GOTO 9999                                                  
      BREF = GRAPH(1:3)                                                         
c      write(0,*) '3D CONTOUR:',iplot,trim(bref)
      IF (BREF(1:1).NE.'C'.and.bref(1:3).ne.'000') GOTO 9999                                           
      IF (IPLOT.EQ.0)       GOTO 4100                                           
C                                                                               
      DO 4110 J = 1, 4                                                          
        S(J) = ' '                                                              
 4110 CONTINUE                                                                  
C                                                                               
C---- ONLY ONE TYPE OF PLOT AVAILABLE ...                                       
C                                                                               
      IF (IPLOT.NE.1) THEN                                                      
        IPLOT = 1                                                               
        S(1) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- SMOOTHING PARAMETER MUST BE 0 OR 1 OR 3,5,7,..                            
C---- 0:OFF  1:WEIGHTED AVERAGE ALONG Y ONLY  3,5,7..: AVERAGE IN BOX           
C                                                                               
      IF (ISMOTH.LT.0.OR.2*(ISMOTH/2).EQ.ISMOTH) THEN                           
        IF (ISMOTH.NE.0) S(2) = '*'                                             
        ISMOTH = 0                                                              
      ENDIF                                                                     
C                                                                               
C---- NUMBER OF CONTOURS TO PLOT                                                
C                                                                               
      IF (NCONT.LT.1 .OR. NCONT.GT.20) THEN                                     
        NCONT = 20                                                              
        S(3) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- DUMMY PARAMETER                                                           
C                                                                               
      IF (IDUM2.NE.0) THEN                                                       
        IDUM2 = 0                                                                
        S(5) = '*'                                                              
      ENDIF                                                                     
C                                                                               
C---- SELECT THE QUANTITY REQUIRED:  NIE, ZB.NBTRUE, ZEFF, TOTAL POWER          
C                                                                               
      IF     (BREF.EQ.'CZ1') THEN                                               
        II = 1                                                                  
      ELSEIF (BREF.EQ.'CZ2') THEN                                               
        II = 2                                                                  
      ELSEIF (BREF.EQ.'CZ3') THEN                                               
        II = 3                                                                  
      ELSEIF (BREF.EQ.'CR ') THEN                                               
        II = 4                                                                  
      ELSEIF (BREF.EQ.'CC ') THEN
        II = 5
      ELSEIF (BREF.EQ.'CP ') THEN
        II = 6
      ELSE                                                                      
        II = 0                                                                  
        GRAPH(1:3) = ' * '                                                      
      ENDIF                                                                     
C  
C---  Ionization state for Clouds and Plrps
C
      IF (II.LE.4.AND.CTIZ.NE.0) THEN
        CTIZ = 0
        S(4) = '*'
      ELSEIF (II.EQ.5.AND.(CTIZ.LT.-1.OR.CTIZ.GT.NIZS)) THEN
        CTIZ = 0
        WRITE(6,*) 'INVALID IONIZATION STATE ',CTIZ,
     >             ' FOR CLOUD CONTOUR PLOT' 
        GOTO 4100
      ELSEIF (II.EQ.6.AND.(CTIZ.LT.-2.OR.CTIZ.GT.NIZS)) THEN
        CTIZ = 0
        WRITE(6,*) 'INVALID IONIZATION STATE ',CTIZ,
     >             ' FOR PLRPS CONTOUR PLOT' 
        GOTO 4100
      ENDIF       
C
C     Find a line matching the ionization state specified for 
C     PLRP plots. This will need to be modified if support for
C     multiple lines / ionization state is to be added.
C
      if (ii.eq.6) then
        tmpiz = 0
        do 4115 il = 1,nls
          if (pizs(il).eq.ctiz) then
            tmpiz = il
            goto 4116
          endif 
4115    continue
4116    continue
        if (tmpiz.eq.0) then
          write(6,*) 'No plrp is available for ionization state:',ctiz
          goto 4100
        endif                 
C
C       reset CTIZ to the valid index into the PLRP array
C
        ctiz = tmpiz 
      endif 
C                                                                               
C---- SET AUX1 ARRAY DIMENSIONED WITHOUT THE "HOLE" AT Y=0 TO EITHER            
C---- ONE OF THE ZEFFS QUANTITIES, OR THE TOTAL POWER LOSS VALUE.               
C---- SETUP CORRESPONDING "YOUTS", "YWIDSS" AND "YSS" ARRAYS FOR THIS.          
C---- NOTE 297:  USE SIN(THETAB), LP+/LP- CORRECTED ZEFF VALUES HERE!           
C                                                                               
      CALL RZERO (AUX1, MAXNXS*(2*MAXNYS))                                      
      CALL RZERO (AUX2, MAXNXS*(2*MAXNYS))                                      
      DO 4160 IY = 1, NYS                                                       
        DO 4150 IX = 1, NXS                                                     
          IF (II.LE.3) THEN                                                     
            AUX1(IX,IY)   = ZEFFS(IX,IY,II+3)                                   
            AUX1(IX,1-IY) = ZEFFS(IX,-IY,II+3)                                  
          ELSEIF (II.EQ.4) THEN                                                 
            DO 4140 IZ = 1, NIZS                                                
              AUX1(IX,IY)   = AUX1(IX,IY)   + POWLS(IX,IY,IZ)                   
              AUX1(IX,1-IY) = AUX1(IX,1-IY) + POWLS(IX,-IY,IZ)                  
 4140       CONTINUE                                                            
          ELSEIF (II.EQ.5) THEN
            AUX1(IX,IY)   = SDLIMS(IX,IY,CTIZ)                                   
            AUX1(IX,1-IY) = SDLIMS(IX,-IY,CTIZ)                                  
c slmod - tmp
            WRITE(7,'(4X,A,3I3,2G10.3)') 
     +        'CC Cloud: ',IX,IY,CTIZ,AUX1(IX,IY),AUX1(IX,1-IY)
c slmod end
          ELSEIF (II.EQ.6) THEN
            AUX1(IX,IY)   = PLRPS(IX,IY,CTIZ)                                   
            AUX1(IX,1-IY) = PLRPS(IX,-IY,CTIZ)                                  
          ENDIF                                                                 
 4150   CONTINUE                                                                
        YOUTS(IY)    = YOUTS1(IY)                                               
        YOUTS(1-IY)  =-YOUTS1(IY)                                               
        YWIDSS(IY)   = YWIDS(IY)                                                
        YWIDSS(1-IY) = YWIDS(IY)                                                
        YSS(IY)      = YS(IY)                                                   
        YSS(-IY)     =-YS(IY)                                                   
 4160 CONTINUE                                                                  
      YSS(0) = 0.0                                                              
      YLAB = YLAB1                                                              
C                                                                               
C---- ADJUST XMIN,XMAX,YMIN,YMAX IF REQUIRED; LOCATE XCON IN ARRAY              
C                                                                               
      XMIN = MAX (XMIN, CAW)                                                    
      XMAX = MIN (XMAX, CA)                                                     
      IF (YMIN.LT.0.0) YMIN = MAX (YMIN,-CL)                                    
      IF (YMAX.GT.0.0) YMAX = MIN (YMAX, CL)                                    
      IXMIN = MAX (1, IPOS (XMIN, XS, NXS-1) -1)                                
      IXMAX = IPOS (XMAX, XS, NXS-1)                                            
      IYMIN = MAX (1-NYS, IPOS (YMIN, YSS(1-NYS), 2*NYS-1) - NYS -1)            
      IYMAX = IPOS (YMAX, YSS(1-NYS), 2*NYS-1) - NYS                            
      IXCON = IPOS (XCON, XS, NXS-1)                                            
C                                                                               
C---- NOTE 297 - USE YT VALUES INSTEAD OF YP+ VALUES FOR THESE PLOTS            
C---- NOTE 299 - FACTOR IS JUST SIN(THETAB), NO LP+/LP- TERM                    
C                                                                               
      IF (II.LE.3) THEN                                                         
        YLAB = YLAB4                                                            
        YMIN = YMIN / CSINTB                                                    
        YMAX = YMAX / CSINTB                                                    
        DO 4170 IY = 1-NYS, NYS                                                 
          YOUTS(IY)  = YOUTS(IY)  / CSINTB                                      
          YWIDSS(IY) = YWIDSS(IY) / CSINTB                                      
          YSS(IY)    = YSS(IY)    / CSINTB                                      
 4170   CONTINUE                                                                
      ENDIF                                                                     
C     WRITE (6,'(''   IY       YSS        YOUTS       YWIDSS'',                 
C    >  /,(I5,3F13.5))')                                                        
C    >  (IY,YSS(IY),YOUTS(IY),YWIDSS(IY),IY=1-NYS,NYS)                          
C                                                                               
C---- PRINT MESSAGE TO CHANNEL 7 AS PERMANENT RECORD OF GRAPH PARAMS.           
C                                                                               
      WRITE (7,'(2X,A21,4F6.2,I6,3(A1,I4),A1)') GRAPH,                          
     >  XMIN,XMAX,YMIN,YMAX,IPLOT,S(1),ISMOTH,S(2),NCONT,S(3),IDUM,S(4)         
      IF (XMIN.GE.XMAX.OR.YMIN.GE.YMAX.OR.II.LE.0)     GOTO 4100                
C                                                                               
C-----------------------------------------------------------------------        
C     CASE WHERE NO SMOOTHING REQUIRED                                          
C     ================================                                          
C     SET AUX2 ARRAY TO AUX1 ARRAY WITH NO SMOOTHING ATTEMPT.                   
C-----------------------------------------------------------------------        
C                                                                               
      SMOOTH = ' '                                                              
      IF (ISMOTH.EQ.0) THEN                                                     
        DO 4200 IX = 1, NXS                                                     
          DO 4200 IY = 1-NYS, NYS                                               
            AUX2(IX,IY) = AUX1(IX,IY)                                           
 4200   CONTINUE                                                                
C                                                                               
C-----------------------------------------------------------------------        
C     CASE WHERE "BOX-TYPE" SMOOTHING IS REQUIRED                               
C     ===========================================                               
C     CALCULATE AVERAGE ZEFFS VALUES INSIDE AN (X,Y) BOX USING 3 X BINS         
C     AND "N" Y BINS, GIVEN BY INPUT DATA AND TYPICALLY 3 OR 5.  NOTE           
C     THAT THE RESPECTIVE BIN SIZES MUST BE ACCOUNTED FOR IN THE                
C     AVERAGING.                                                                
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (ISMOTH.GT.1) THEN                                                 
        DO 4370 IX = 1, NXS                                                     
          DO 4370 IY = 1-NYS, NYS                                               
            JXMIN = MAX (1,     IX-1)                                           
            JXMAX = MIN (NXS,   IX+1)                                           
            JYMIN = MAX (1-NYS, IY-ISMOTH/2)                                    
            JYMAX = MIN (NYS,   IY+ISMOTH/2)                                    
            TOTAV = 0.0                                                         
            DO 4360 JX = JXMIN, JXMAX                                           
              DO 4360 JY = JYMIN, JYMAX                                         
                AUX2(IX,IY) = AUX2(IX,IY) +                                     
     >            AUX1(JX,JY) * XWIDS(JX) * YWIDSS(JY)                          
                TOTAV = TOTAV + XWIDS(JX) * YWIDSS(JY)                          
 4360       CONTINUE                                                            
            AUX2(IX,IY) = AUX2(IX,IY) / TOTAV                                   
 4370   CONTINUE                                                                
        WRITE (SMOOTH,'(''AVERAGED OVER 3 X AND'',I2,'' Y BINS   *'')')         
     >    ISMOTH                                                                
C                                                                               
C-----------------------------------------------------------------------        
C     CASE WHERE "WEIGHTED AVERAGE OVER Y RANGE" SMOOTHING REQUIRED             
C     =============================================================             
C     NEW FOR JANUARY 1989!                                                     
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (ISMOTH.EQ.1) THEN                                                 
        TOTAV = 0.0                                                             
        DO 4400 J = -NAVS, NAVS                                                 
          TOTAV = TOTAV + AVS(IABS(J))                                          
 4400   CONTINUE                                                                
        WMIN = 1.E10                                                            
        DO 4410 IY = 1-NYS, NYS                                                 
          WMIN = MIN (WMIN, YWIDSS(IY))                                         
 4410   CONTINUE                                                                
        DO 4420 IX = 1, NXS                                                     
          DO 4420 IY = 1-NYS, NYS                                               
            DO 4420 J = -NAVS, NAVS                                             
              JY = IPOS (YOUTS(IY)+REAL(J)*WMIN,YSS(1-NYS),2*NYS-1) -NYS        
              AUX2(IX,IY) = AUX2(IX,IY) + AUX1(IX,JY) *                         
     >                                            AVS(IABS(J)) / TOTAV          
 4420   CONTINUE                                                                
        WRITE (SMOOTH,'(''SMOOTHED,'',I3,'' WEIGHTS, DY'',F7.4,                 
     >    ''  *'')') 2*NAVS+1,WMIN                                              
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     ALL CASES CONTINUE HERE                                                   
C     =======================                                                   
C     NOTE 282:  FOR X > XCON, ASSUME UNIFORM DISTIBUTION ALONG Y               
C     FIND MIN/MAX VALUES OF (POSSIBLY SMOOTHED) AUX2 ARRAY WITHIN THE          
C     PLOTTING REGION SPECIFIED AND SET CLEVLS ARRAY TO                         
C     PLOT CONTOURS AT MAXIMUM VALUE * 1/4, 1/2, AND 3/4 FOR EXAMPLE.           
C     ROUTINE GRTSET DRAWS THE AXES AND GRAPH FRAMEWORK, TITLES, ETC            
C     AND SUBSEQUENT CALLS TO LIM_GRCONT DRAW THE VARIOUS CONTOUR LINES.            
C-----------------------------------------------------------------------        
C                                                                               
      DO 4480 IX = 1, NXS                                                       
        IF (IX.GT.IXCON) THEN                                                   
          TOTAV = 0.0                                                           
          DO 4470 IY = 1-NYS, NYS                                               
            TOTAV = TOTAV + AUX1(IX,IY) * YWIDSS(IY)                            
 4470     CONTINUE                                                              
          DO 4475 IY = 1-NYS, NYS                                               
            AUX2(IX,IY) = TOTAV / (2.0 * CTWOL)                                 
 4475     CONTINUE                                                              
        ENDIF                                                                   
 4480 CONTINUE                                                                  
C                                                                               
      ZEFMAX =-BIGG                                                             
      ZEFMIN = BIGG                                                             
      DO 4505 IX = 1, NXS                                                       
        IF (XOUTS(IX).LT.XMIN .OR. XOUTS(IX).GT.XMAX) GOTO 4505                 
        DO 4500 IY = 1-NYS, NYS                                                 
          IF (YOUTS(IY).LT.YMIN .OR. YOUTS(IY).GT.YMAX) GOTO 4500               
          ZEFMAX = MAX (ZEFMAX, AUX2(IX,IY))                                    
          ZEFMIN = MIN (ZEFMIN, AUX2(IX,IY))                                    
 4500   CONTINUE                                                                
 4505 CONTINUE                                                                  
      IF (II.EQ.3) ZEFMIN = REAL (CIZB)                                         
C                                                                               
      DO 4510 IL = 1, NCONT                                                     
        FRAC = REAL(IL) / REAL(NCONT+1)                                         
        CLEVLS(IL) = ZEFMIN + FRAC * (ZEFMAX - ZEFMIN)                          
        IF (II.EQ.3) THEN                                                       
        WRITE (NAMES(IL),'(''    ZB'',SP,F5.2,''*(ZMAX-ZB)='',SS,F6.2)')        
     >    FRAC,CLEVLS(IL)                                                       
        ELSE                                                                    
          WRITE (NAMES(IL),'(''       '',1P,G9.2)') CLEVLS(IL)                  
        ENDIF                                                                   
        IF (SMOOTH.NE.' ') NAMES(IL)(31:31) = '*'                               
 4510 CONTINUE                                                                  
C                                                                               
      PLANE = ' '                                                               
      ANLY  = ' '                                                               
      IF (CANAL.LT.CA)                                                          
     >  WRITE (ANLY,'(''ANALYTIC EXTENSION FOR X >'',F6.3)') CANAL              
      IF (II.LE.4) THEN 
        REF = ELABS(5+II)(5:11) // ' CONTOURS'                                    
      ELSEIF (II.EQ.5) THEN     
        WRITE(REF,'(''IONIZ:'',I2,'' CONTOURS'')') CTIZ
      ELSEIF (II.EQ.6) THEN
        WRITE(REF,'(''PLRPS:'',G6.2,'' CONTOURS'')') PLAMS(CTIZ)
      ENDIF  
      IF (IXCON.LT.IXMAX) WRITE (SMOOTH(37:72),'(A,F5.2)')                      
     >  'AVERAGED OVER ALL Y FOR X >',XCON                                      
      CALL LIM_GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XMIN,                              
     >             XMAX,YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NCONT)        





         DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
             IF (AUX2(IX,IY).GT.0.0) TMIN=MIN(TMIN,AUX2(IX,IY))
             TMAX=MAX(AUX2(IX,IY),TMAX)
           ENDDO
         ENDDO

        CLEVLS(1)  = 0.05*(TMAX-TMIN) + TMIN
        CLEVLS(2)  = 0.15*(TMAX-TMIN) + TMIN
        CLEVLS(3)  = 0.25*(TMAX-TMIN) + TMIN
        CLEVLS(4)  = 0.35*(TMAX-TMIN) + TMIN
        CLEVLS(5)  = 0.45*(TMAX-TMIN) + TMIN
        CLEVLS(6)  = 0.55*(TMAX-TMIN) + TMIN 
        CLEVLS(7)  = 0.65*(TMAX-TMIN) + TMIN
        CLEVLS(8)  = 0.75*(TMAX-TMIN) + TMIN
        CLEVLS(9)  = 0.85*(TMAX-TMIN) + TMIN
        CLEVLS(10) = 0.95*(TMAX-TMIN) + TMIN




      DO 4520 IL = 1, NCONT                                                     
c slmod begin
c        WRITE(0,*) IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS,
c     +    CLEVLS(IL)
c     +    XOUTS,
c     +    YOUTS(1-MAXNYS),
c     +    NAMES(IL)
c        WRITE(0,*) ' '
c
        WRITE (NAMES,'(4X,1P,E8.1)') CLEVLS(IL)

        CALL LIM_GRCONT95 (AUX2,IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,                       
     >               MAXNYS,CLEVLS(IL),XOUTS,YOUTS,NAMES)         
c        CALL LIM_GRCONT (AUX2,IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,                       
c     >               MAXNYS,CLEVLS(IL),XOUTS,YOUTS(1-MAXNYS),NAMES(IL))         
c slmod end
4520  CONTINUE                                                                  

      CALL CONTIL (AUX2,IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,                       
     >            MAXNYS,CLEVLS,1,NCONT,XOUTS,YOUTS,NAMES)
      CALL FRAME                                                                
C                                                                               
      GOTO 4100                                                                 
C     *********                                                                 
C                                                                               
C-----------------------------------------------------------------------        
C     FINAL SECTION - PRINT SUMMARY TABLE                                       
C     IF DETAILS WERE RECORDED FOR MORE ITERATIONS, LEAP BACK TO 10.            
C-----------------------------------------------------------------------        
C                                                                               
 9999 CONTINUE                                                                  

      write(0,*) 'Printing Summaries'
      CALL PRB                                                                  
      CALL PRC ('* INDICATES PLOT OPTION ADJUSTED BY PROGRAM')                  
      CALL PRB                                                                  
C                                                                               
C---- CALCULATE TOTAL QUANTITIES AND PRINT TABLE OF RESULTS.                    
C                                                                               
      IPLANE = 0                                                                
      IFOLD  = 0                                                                

      CALL RINTXY (SDLIMS,SDLIM3,IPLANE,NIZS,XYINTS,IFOLD,RV,XV,YV,0,           
     >             SSS,NYSLIM,FP,FT)                                            

      DO 10010 IZ = -2, NIZS+1                                                  
        TOTALS(1,IZ) = XYINTS(0,IZ)                                             
10010 CONTINUE                                                                  
C                                                                               
      CALL RINTXY (POWLS,POWL3,IPLANE,NIZS,XYINTS,IFOLD,RV,XV,YV,0,             
     >             SSS,NYSLIM,FP,FT)                                            
      DO 10020 IZ = -2, NIZS+1                                                  
        TOTALS(2,IZ) = XYINTS(0,IZ)                                             
10020 CONTINUE                                                                  
C                                                                               
      CALL RINTXY (LINES,LINE3,IPLANE,NIZS,XYINTS,IFOLD,RV,XV,YV,0,             
     >             SSS,NYSLIM,FP,FT)                                            
      DO 10030 IZ = -2, NIZS+1                                                  
        TOTALS(3,IZ) = XYINTS(0,IZ)                                             
10030 CONTINUE                                                                  
C                                                                               
      CALL RINTXY (TIZS ,TIZ3 ,IPLANE,NIZS,XYINTS,IFOLD,RV,XV,YV,0,             
     >             SSS,NYSLIM,FP,FT)                                            
      DO 10040 IZ = -2, NIZS+1                                                  
        TOTALS(4,IZ) = XYINTS(0,IZ)                                             
10040 CONTINUE                                                                  
C                                                                               
      WRITE (7,9003)                                                            
     >  (ZLABS1(IZ)(5:11),(TOTALS(IT,IZ),IT=1,4),IZ=-2,NIZS+1)                  
C     ___________________________                                               
      IF (ITER.LT.NITERS) GOTO 10                                               
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                               
      CALL PRB                                                                  
      TIME = ZA02AS (1) - TIME1                                                 
      WRITE (6,'('' OUT3: TOTAL TIME USED ='',G11.4,'' SEC'')') TIME            
      CALL PRR ('TOTAL TIME USED (SECONDS)  ',TIME)                             
      CALL PRB                                                                  
      CALL GREND                                                                

c
c     Clean up dynamically allocated storage 
c     
      call deallocate_dynamic_storage
c

      WRITE(0,*) 'Done  OUT3'

      STOP                                                                      
C                                                                               
 9001 FORMAT('RV=',F7.2,' XV''=',F7.2,' YV''=',F7.2)                            
 9003 FORMAT(//2X,'TOTALS            IMPURITY  POWER LOSS LINE RAD''N',         
     >  ' IONISATION',/1X,60('-'),/,(2X,A16,1P,4(1X,G10.3)))                    
 9004 FORMAT(//1X,'OUT3: MAXOS=',I9,'  NOS=',I9,//1X,                           
     >  '  IO     OYWIDS   OYS    OYOUTS    ODWIDS   ODS    ODOUTS',/1X,        
     >  65('-'),/,(1X,I5,6F9.4))                                                
 9006 FORMAT(1X,'OBSERVATION POINT: RV=',1P,G11.4,',  XV=',G11.4,               
     >  ',  YV=',G11.4,',  YLIM=+/-',0P,F7.3)                                   
 9010 FORMAT(//2X,'IONISATION : AVERAGE DEPTHS',/,                              
     >  (2X,A16,'   X =',F6.3,'M'))                                             
      END                                                                       
c
c
c
      subroutine allocate_dynamic_storage
      ! routine to allocate dynamic storage at fixed sizes - eventually update to allow dynamic size definitions
      use mod_adas_data_spec
      use mod_comnet
      use mod_comt2
      use mod_comtor
      use mod_comvu
      use mod_comxyt
      use mod_coords
      use mod_dynam2
      use mod_dynam3
      use mod_expt_data
      use mod_gcom1
      use mod_reader
      use mod_rtheta
      use mod_slcom
      use mod_colours
      use mod_comgra
      use mod_grminfo
      use mod_limpoly
      use mod_pindata
      use mod_slout
      implicit none


      ! OUT
      call allocate_mod_adas_data_spec
      call allocate_mod_comnet
      call allocate_mod_comt2
      call allocate_mod_comtor
      call allocate_mod_comvu
      call allocate_mod_comxyt
      call allocate_mod_coords
      call allocate_mod_dynam2 
      call allocate_mod_dynam3
      call allocate_mod_expt_data
      call allocate_mod_gcom1
      call allocate_mod_reader
      call allocate_mod_rtheta
      call allocate_mod_slcom
      !call allocate_mod_unstructured
      call allocate_mod_colours
      call allocate_mod_comgra
      call allocate_mod_grminfo
      call allocate_mod_limpoly
      call allocate_mod_pindata
      call allocate_mod_slout

      
      return
      end
c
c
c     
      subroutine deallocate_dynamic_storage
      use mod_adas_data_spec
      use mod_comnet
      use mod_comt2
      use mod_comtor
      use mod_comvu
      use mod_comxyt
      use mod_coords
      use mod_dynam2
      use mod_dynam3
      use mod_expt_data
      use mod_gcom1
      use mod_reader
      use mod_rtheta
      use mod_slcom
      use mod_colours
      use mod_comgra
      use mod_grminfo
      use mod_limpoly
      use mod_pindata
      use mod_slout
      implicit none


      ! OUT
      call deallocate_mod_adas_data_spec
      call deallocate_mod_comnet
      call deallocate_mod_comt2
      call deallocate_mod_comtor
      call deallocate_mod_comvu
      call deallocate_mod_comxyt
      call deallocate_mod_coords
      call deallocate_mod_dynam2 
      call deallocate_mod_dynam3
      call deallocate_mod_expt_data
      call deallocate_mod_gcom1
      call deallocate_mod_reader
      call deallocate_mod_rtheta
      call deallocate_mod_slcom
      !call deallocate_mod_unstructured
      call deallocate_mod_colours
      call deallocate_mod_comgra
      call deallocate_mod_grminfo
      call deallocate_mod_limpoly
      call deallocate_mod_pindata
      call deallocate_mod_slout

      return
      end
