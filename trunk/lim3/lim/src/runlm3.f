c     -*-Fortran-*-
c
      PROGRAM RUNLM3                                                            
      use mod_params
      use yreflection
      use mod_dynam3
      use mod_comtor
      use mod_cadas
      use mod_comtau
      use mod_comxyt
      use mod_coords
      use mod_printr
      IMPLICIT  none
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS PROGRAM READS IN A SET OF INPUT VALUES,                          
C       CALLS LIM3 PASSING IT THE VALUES                                        
C       DUMPS THE OUTPUT IN AN EXTERNAL FILE.                                   
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
c      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
c      INCLUDE 'comtau'                                                          
C     INCLUDE (COMTAU)                                                          
c      INCLUDE 'coords'                                                          
C     INCLUDE (COORDS)                                                          
c      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
c      INCLUDE 'printr'                                                          
c slmod begin
c      INCLUDE 'dynam3'
c      include 'cadas'
c slmod end
C     INCLUDE (PRINTR)                                                          
C                                                                               
      INTEGER        IERR,NM,NC,IPOS,IP,ICHAR,IZ,NTBS,NTIBS,NNBS
      INTEGER        NCVS
      INTEGER        NYMFS,NQS            
      INTEGER        IX,IY,IQX,IQY,NIZS,J,IGEOM,IT,KFAIL(1),NITERS,NOUT         
      INTEGER        IMODE,NIMPS,PIZS(MAXNLS),NLS,ITER,MXXNPS                   
      INTEGER        IMPADD,IMPCF
      INTEGER        KNEUTA,KNEUTB,KNEUTC,KNEUTE,NRAND                          
      REAL           QTIM,CPULIM,PLAMS(MAXNLS),X,Y,TOTLPD                      
      REAL           XWIDM,YWIDM,FSRATE                                         
      REAL           IONTIM,NEUTIM,STATIM,TOTTIM,ZA02AS,PMASS(1)         
c      REAL           IONTIM,NEUTIM,STATIM,TOTTIM,ZA02AS,DEGRAD,PMASS(1)         
      REAL           FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)                
      REAL           MU,CBETA1,CBETA2
      CHARACTER      TITLE*80,JOB*72,JFCB*176,COMENT*77,NUMBER(20)*4            
      CHARACTER*8    SYSTIM,SYSDAT,VSN,DSN(3)                                   
      CHARACTER      PRINPS(-MAXNPS-1:MAXNPS)*7                                 
      DOUBLE PRECISION SEED,DEFACT                                              
c slmod begin
      INTEGER IL
      REAL    NUMSUM,VOLSUM,TMPVOL
c slmod end
c      DATA    DEGRAD / 0.017453292 /,      NOUT / 8 /                           
c
      real ran
c
      DATA    NOUT / 8 /                           
      DATA    NUMBER /' 1ST',' 2ND',' 3RD',' 4TH',' 5TH',' 6TH',' 7TH',         
     >  ' 8TH',' 9TH','10TH','11TH','12TH','13TH','14TH','15TH','16TH',         
     >  '17TH','18TH','19TH','20TH'/                                            
c
c     Runtime debug tracing
c
c      call 

C
C     INITIALIZE VARIABLES THAT REQUIRE IT
C
      IMPCF = 0
c
c     Allocate dynamic storage
c
      call allocate_dynamic_storage
C                                                                               
C-----------------------------------------------------------------------        
C     READ IN DATA - MAKE TEMP COPIES OF SOME INPUT FLAGS OVERWRITTEN           
C     IN LIM3 SELF-SPUTTERING CASES, SINCE THEY MAY BE REQUIRED FOR             
C     SUBSEQUENT ITERATIONS.                                                    
C-----------------------------------------------------------------------        
C                                                                               
      CALL XUFLOW (0)                                                           
      STATIM = ZA02AS (1)                                                       
      IONTIM = 0.0                                                              
      NEUTIM = 0.0                                                              
C                                                                               
      IERR = 0                                                                  
      CALL READIN (TITLE,IGEOM,IMODE,NIZS,NIMPS,IMPADD,                     
     >             FSRATE,QTIM,CPULIM,IERR,NTBS,NTIBS,NNBS,NYMFS,
     >             NCVS,NQS,NITERS)          
C                                                                               
      KNEUTA = CNEUTA                                                           
      KNEUTB = CNEUTB                                                           
      KNEUTC = CNEUTC                                                           
      KNEUTE = CNEUTE                                                           
C
C     SET THE INITIAL NEUTRAL VEL/ANG FLAG IF REQ'D
C
      IF (NVAOPT.EQ.-1) NVAOPT = CNEUTC
C
C
C-----------------------------------------------------------------------
C        CALCULATE CUMULATIVE ATOMIC DISTRIBUTION FUNCTIONS FOR LAUNCH
C        OPTIONS 7 AND 8
C-----------------------------------------------------------------------
C
      WRITE(6,*) 'CLPD:',CLPD               
      IF (CLPD.NE.0) THEN
         TOTLPD = LPDION(1,2) 
         LPDCUM(1) = LPDION(1,2)
         DO 10 IX = 2,CLPD
            TOTLPD = TOTLPD + LPDION(IX,2)
            LPDCUM(IX) = TOTLPD
 10      CONTINUE   
         DO 20 IX = 1, CLPD
            LPDCUM(IX) = LPDCUM(IX) / TOTLPD
            WRITE(6,*) 'LPD:',LPDION(IX,1),LPDION(IX,2),LPDCUM(IX)
 20      CONTINUE  
      ENDIF
C
C---------------------------------------------------------------------
C     CALCULATE ALPHAe AND BETAi COEFFICIENTS FOR TEMPERATURE 
C     GRADIENT FORCES
C
C---------------------------------------------------------------------
C
      MU = CRMI / (CRMI+CRMB)
      CBETA1 = 5.0 * SQRT(2.0) *(1.1*MU**2.5-0.35*MU**1.5)
      CBETA2 = 2.6 - 2*MU + 5.4*MU**2 
      DO 25 IZ = 1,MAXIZS
         IF (CIOPTL.EQ.0) THEN 
            CALPHE(IZ) = 0.0
         ELSEIF (CIOPTL.EQ.1) THEN 
            CALPHE(IZ) = 0.71*IZ*IZ
         ENDIF
         IF (CIOPTM.EQ.0) THEN 
            CBETAI(IZ) = 0.0
         ELSEIF (CIOPTM.EQ.1) THEN 
            CBETAI(IZ) = -3.0*(1.0-MU-CBETA1*IZ*IZ)/CBETA2   
         ENDIF
25    CONTINUE    
C                                                                               
C-----------------------------------------------------------------------        
C        CALCULATE BINS AND INDEX ARRAYS                                        
C-----------------------------------------------------------------------        
C                                                                               
C---- X VALUES HAVE BEEN READ FOR -AW < X < A  (EXCLUSIVE).  SET NEXT &         
C---- LAST X VALUE TO A.                                                        
C                                                                               
      NXS     = NXS + 1                                                         
      XS(NXS) = CA                                                              
C                                                                               
C---- Y VALUES HAVE BEEN READ FOR 0 < Y < L  (EXCLUSIVE).  SET NEXT             
C---- Y VALUE TO L, THEN GENERATE AN EQUIVALENT SPACED MESH FOR THE             
C---- REGION L < Y < 2L.  FINALLY SET LAST Y VALUE TO 2L.                       
C---- NOTE THAT WITH THE "NO MIRROR" FEATURE THE WHOLE Y RANGE WILL BE          
C---- DUPLICATED FOR -2L < Y < 0.                                               
C---- SET NY3D VALUE                                                            
C                                                                               
      NYS     = NYS + 1                                                         
      YS(NYS) = CL                                                              
      DO 100 IY = NYS+1, 2*NYS-1                                                
        YS(IY) = CL+CL - YS(2*NYS-IY)                                           
  100 CONTINUE                                                                  
      NYS     = 2*NYS                                                           
      YS(NYS) = CL+CL                                                           
      NY3D    = MIN (NYS, MAXY3D)                                               
C                                                                               
C---- APPLICATION OF NOTE 182: TOROIDAL LIMITER GEOMETRY OPTION.                
C---- CONNECTION LENGTH NOW BECOMES L.SIN(THETAB), AND IF THETAB                
C---- IS GIVEN AS 90 DEGREES THIS REVERTS TO THE POLOIDAL CASE.                 
C---- ADJUST SOME Y DEPENDENT QUANTITIES ACCORDINGLY.                           
C---- NOTE 293 ETC:  ALL CHANGE, LP+, LP-, LT ETC .......                       
C---- NOTE 349: ENSURE ALL THIS IS IRRELEVANT FOR POLOIDAL GEOMETRY.            
C                                                                               
      IF (CTHETB.EQ.90.0) THEN                                                  
        CAP  = CL                                                               
        CWL  = 0.5                                                              
      ENDIF                                                                     
      CSINTB = SIN (DEGRAD*CTHETB)                                              
      CLFACT = CAP / (2.0*CWL) / (CL*CSINTB)                                    
      DO 105 IY = 1, NYS                                                        
        YS(IY) = YS(IY) * CAP / (2.0*CWL) / CL                                  
  105 CONTINUE                                                                  
      CL     = CAP / (2.0*CWL)                                                  
      CTWOL  = CL+CL                                                            
      CHALFL = 0.5 * CL                                                         
      C3HALFL= 1.5 * CL
c
c     jdemod - initialize the reflection option
c
      call init_reflection(ctwol,ca,caw)
C                                                                               
C---- CALCULATE SET OF X AND Y POSITIONS FOR WHICH                              
C---- FACTORS ARE TO BE CALCULATED                                              
C---- THE X RANGE IS SPLIT INTO INBOARD AND OUTBOARD REGIONS, EACH              
C---- USING A SEPARATE SCALE.  FOR EXAMPLE, IF CAW=-0.1M & CA=1.5M              
C---- & NQXSO=100, THEN WE HAVE 100 DIVISIONS OF 0.001M OUTBOARD AND            
C---- IF NQXSI=100 WE HAVE 100 DIVISIONS OF 0.015M INBOARD.                     
C                                                                               
      XSCALO = -REAL (NQXSO) / CAW                                              
      XSCALI =  REAL (NQXSI) / CA                                               
      YSCALE =  REAL (NQYS) / CTWOL                                             
C                                                                               
      DO 110 IQX = 1-NQXSO, 0                                                   
         QXS(IQX) = (REAL(IQX) - 0.5) / XSCALO                                  
  110 CONTINUE                                                                  
      DO 120 IQX = 1, NQXSI                                                     
         QXS(IQX) = (REAL(IQX) - 0.5) / XSCALI                                  
  120 CONTINUE                                                                  
      DO 130 IQY = 1, NQYS                                                      
         QYS(IQY) = (REAL(IQY) - 0.5) / YSCALE                                  
  130 CONTINUE                                                                  
C                                                                               
C---- SET INDEXING ARRAYS IQXS, IQYS AND MIDPOINTS XOUTS, YOUTS                 
C                                                                               
      DO 140 IX = 1, NXS                                                        
        IF (IX.EQ.1) THEN                                                       
          X = 0.5 * (XS(1) + CAW)                                               
        ELSE                                                                    
          X = 0.5 * (XS(IX) + XS(IX-1))                                         
        ENDIF                                                                   
        XOUTS(IX) = X                                                           
        IF (X.GE.0.0) THEN                                                      
          IQXS(IX) = INT (X * XSCALI) + 1                                       
        ELSE                                                                    
          IQXS(IX) = INT (X * XSCALO)                                           
        ENDIF                                                                   
  140 CONTINUE                                                                  
C                                                                               
      DO 150 IY = 1, NYS                                                        
        IF (IY.EQ.1) THEN                                                       
          Y = 0.5 * YS(1)                                                       
        ELSE                                                                    
          Y = 0.5 * (YS(IY) + YS(IY-1))                                         
        ENDIF                                                                   
        YOUTS(IY) = Y                                                           
        YOUTS(-IY) = -Y                                                         
        IQYS(IY) = INT (Y * YSCALE) + 1                                         
  150 CONTINUE                                                                  
C                                                                               
C---- CALCULATE X BIN WIDTHS: STORE MIN VALUE IN XWIDM                          
C---- CALCULATE CYLINDRICAL GEOMETRY CORRECTION FACTORS FOR X POINTS -          
C---- THESE WILL BE USED IN NORMALISING DDLIMS,DDLIM3,TIZS,ETC                  
C                                                                               
      XWIDM = CA - CAW                                                          
      DO 160 IX = 1, NXS                                                        
        IF (IX.EQ.1) THEN                                                       
          XWIDS(1) = XS(1) - CAW                                                
        ELSE                                                                    
          XWIDS(IX) = XS(IX) - XS(IX-1)                                         
        ENDIF                                                                   
        XWIDM = MIN (XWIDM, XWIDS(IX))                                          


c
c       Check for VALID Xwids values:
c
        if (xwids(ix).le.0.0) then 
           write(0,'(a,2i6,2(1x,g12.5))') 
     >          'LIM ERROR: XWIDS IS <= ZERO: ',
     >           ix,nxs,xs(ix),xs(ix-1)
           write(6,'(a,2i6,2(1x,g12.5))') 
     >          'LIM ERROR: XWIDS IS <= ZERO: ',
     >           ix,nxs,xs(ix),xs(ix-1)
           write(0,'(a,2(1x,g12.5))') 'LIM ERROR: '//
     >                'DO NOT SET BIN BOUNDARIES TO EDGE VALUES',
     >                ca,caw
           write(6,'(a,2(1x,g12.5))') 'LIM ERROR: '//
     >                'DO NOT SET BIN BOUNDARIES TO EDGE VALUES',
     >                ca,caw
           write(0,'(a)') 'LIM IS EXITING'
           write(6,'(a)') 'LIM IS EXITING'
           stop
        endif 


        IF     (IGEOM.EQ.0) THEN                                                
          XCYLS(IX) = 1.0                                                       
        ELSEIF (IGEOM.EQ.1) THEN                                                
          XCYLS(IX) = (CA - (XS(IX)-0.5*XWIDS(IX))) / CA                        
        ENDIF                                                                   
  160 CONTINUE                                                                  
C                                                                               
C---- CALCULATE Y BIN WIDTHS                                                    
C                                                                               
      YWIDM = CTWOL                                                             
      DO 170 IY = 1, NYS                                                        
        IF (IY.EQ.1) THEN                                                       
          YWIDS(1) = YS(1)                                                      
        ELSE                                                                    
          YWIDS(IY) = YS(IY) - YS(IY-1)                                         
        ENDIF                                                                   
        YWIDM = MIN (YWIDM, YWIDS(IY))                                          
c
c       Check for VALID Ywids values:
c
        if (ywids(iy).le.0.0) then 
           write(0,'(a,2i6,2(1x,g12.5))') 
     >          'LIM ERROR: YWIDS IS <= ZERO: ',
     >           iy,nys,ys(iy),ys(iy-1)
           write(6,'(a,2i6,2(1x,g12.5))') 
     >          'LIM ERROR: YWIDS IS <= ZERO: ',
     >           iy,nys,ys(iy),ys(iy-1)
           write(0,'(a,2(1x,g12.5))') 'LIM ERROR: '//
     >                'DO NOT SET BIN BOUNDARIES TO EDGE VALUES',
     >                0.0,CL
           write(6,'(a,2(1x,g12.5))') 'LIM ERROR: '//
     >                'DO NOT SET BIN BOUNDARIES TO EDGE VALUES',
     >                0.0,CL
           write(0,'(a)') 'LIM IS EXITING'
           write(6,'(a)') 'LIM IS EXITING'
           stop
        endif 

  170 CONTINUE                                                                  
C                                                                               
C---- SHIFT "NEAR LIMITER" PARAMETERS TO COINCIDE WITH A BIN BOUNDARY           
C                                                                               
      IX = IPOS (CXNEAR*0.99999, XS, NXS-1)                                     
      IY = IPOS (CYNEAR*0.99999, YS, NYS-1)                                     
      CXNEAR = XS(IX)                                                           
      CYNEAR = YS(IY)                                                           
      CYFAR  = CTWOL - CYNEAR                                                   
C                                                                               
C---- SETUP SYSTEM OF P BINS USING +/-FIRST EXTENT & SUBSEQUENT WIDTH.          
C---- EXAMPLE: CPFIR=1MM, CPSUB=5MM, MAXNPS=4, THEN WE HAVE                     
C----     IP        -4   -3   -2   -1   0   1    2    3    4                    
C----     PS(IP)   -16  -11   -6   -1   1   6   11   16   1.E75MM               
C---- WHERE PS ARE THE UPPER BOUNDS ON EACH P BIN.  HENCE CENTRE BIN            
C---- EXTENDS FROM -1:1MM,  THE FIRST BIN FROM -INFINITY:-16MM,                 
C---- AND THE LAST BIN FROM 16:1.E75MM  (IE. TO +INFINITY)                      
C---- USING "MXXNPS" BELOW PREVENTS IBM COMPILER GENERATING A WARNING           
C---- FOR LOOP 180 IF "MAXNPS" HAPPENS TO BE 1.                                 
C                                                                               

c
      PS(0)  = CPFIR                                                            
      PS(-1) = -CPFIR                                                           
      MXXNPS = MAXNPS - 1                                                       
      DO 180 IP = 1, MXXNPS                                                     
        PS(IP)    = PS(IP-1) + CPSUB                                            
        PS(-1-IP) = PS(-IP)  - CPSUB                                            
  180 CONTINUE                                                                  
      ! jdemod Number too large - should use HI here (or MACHHI)
      ! PS(MAXNPS) = 1.0E75                                                        
      PS(MAXNPS) = HI                                                        
C                                                                               
C---- CALCULATE P BIN WIDTHS.  SET A NOMINAL WIDTH FOR OUTER BINS               
c                        

      PWIDS(-MAXNPS) = 10000.0                                                  
      PWIDS(MAXNPS)  = 10000.0                                                  

      DO 182 IP = 1-MAXNPS, MAXNPS-1                                            
        PWIDS(IP) = PS(IP) - PS(IP-1)                                           
c        POUTS(IP) = (PS(IP)+PS(IP-1))/2.0
  182 CONTINUE                                                                  

c
c     If P reflection option is turned on then set/use the preflect_bound
c     value
c
      if (preflect_opt.eq.1) then 
         ! set the reflection boundary based on the P mesh
         if (preflect_bound.eq.0.0) then 
            preflect_bound = abs(ps(-maxnps))
            write(6,'(a,5(1x,g12.5))') 'Calculating P '//
     >              'reflection boundary:', ps(-maxnps), 
     >               abs(ps(-maxnps)),cpsub,preflect_bound
            ! overwrite the pwids for the two end bins in this case
            pwids(-maxnps) = cpsub
            pwids(maxnps) = cpsub
         endif
      endif

C                                                                               
C---- SET UP PRINPS CHARACTER STINGS FOR OUTPUT P BIN SIZES IN LIM3             
C                                                                               
      PRINPS(-MAXNPS-1) = '-INFNTY'                                             
      DO 184 IP = -MAXNPS,MAXNPS-1                                              
        WRITE (PRINPS(IP),'(F7.4)') PS(IP)                                      
  184 CONTINUE                                                                  
      PRINPS(MAXNPS) = ' INFNTY'                                                
C                                                                               
C---- SETUP QYS POSITIONS FOR PRINTOUTS ETC                                     
C                                                                               
      IQY0   = 1                                                                
      IQYL8  = NQYS/16                                                          
      IQYL4  = NQYS/8                                                           
      IQYL2  = NQYS/4                                                           
      IQYL   = NQYS/2                                                           
      IQY3L2 = (3*NQYS)/4                                                       
      IQY7L4 = (7*NQYS)/8                                                       
      IQY158 = (15*NQYS)/16                                                     
      IQY2L  = NQYS                                                             
      IQYS2  = NINT (0.5 * REAL(NQYS) * CYSTAG / CTWOL)                         
      IQY3S2 = NINT (1.5 * REAL(NQYS) * CYSTAG / CTWOL)                         
      IQY5S2 = NINT (2.5 * REAL(NQYS) * CYSTAG / CTWOL)                         
      IQY5S  = NINT (5.0 * REAL(NQYS) * CYSTAG / CTWOL)                         
      IQY15S = NINT(15.0 * REAL(NQYS) * CYSTAG / CTWOL)                         
C                                                                               
C---- SETUP QXS POSITIONS FOR PRINTOUTS ETC ...                                 
C                                                                               
      IQXA   = NQXSI                                                            
      IQXA2  = NQXSI/2                                                          
      IQXA4  = NQXSI/4                                                          
      IQXIN  = 1                                                                
      IQXOUT = 0                                                                
      IQXAW4 = -NQXSO/4                                                         
      IQXAW2 = -NQXSO/2                                                         
      IQXAW  = 1-NQXSO                                                          
      IQXFAC = 0                                                                
      IF (-CLNOUG.GT.CAW) IQXFAC = NINT (REAL(NQXSO)*CLNOUG/CAW)                
C                                                                               
C---- SETUP YS POSITIONS FOR PRINTOUTS                                          
C                                                                               
      IY0    = 1                                                                
      IY0LT  = -1
      IYL8   = IPOS (0.125*CL,YS, NYS-1)                                        
      IYL4   = IPOS (0.25*CL, YS, NYS-1)                                        
      IYL2   = IPOS (0.5*CL,  YS, NYS-1)                                        
      IYL    = IPOS (CL,      YS, NYS-1)                                        
      IY3L2  = IPOS (1.5*CL,  YS, NYS-1)                                        
      IY7L4  = IPOS (1.75*CL, YS, NYS-1)                                        
      IY158  = IPOS (1.875*CL,YS, NYS-1)                                        
      IY2L   = NYS                                                              
      IYS2   = IPOS (CYSTAG/2.0, YS, NYS-1)                                     
      IY3S2  = IPOS (1.5*CYSTAG, YS, NYS-1)                                     
      IY5S2  = IPOS (2.5*CYSTAG, YS, NYS-1)                                     
C                                                                               
C---- SETUP XS POSITIONS FOR PRINTOUT                                           
C                                                                               
      IXA    = NXS                                                              
      IXA2   = IPOS (0.5*CA,   XS, NXS-1)                                       
      IXA4   = IPOS (0.25*CA,  XS, NXS-1)                                       
      IXIN   = IPOS (1.0E-10,  XS, NXS-1)                                       
      IXOUT  = IPOS (-1.E-10,  XS, NXS-1)                                       
      IXAW4  = IPOS (0.25*CAW, XS, NXS-1)                                       
      IXAW2  = IPOS (0.5*CAW,  XS, NXS-1)                                       
      IXAW   = 1                                                                
      IXFAC  = 0                                                                
      IF (-CLNOUG.GT.CAW) IXFAC = IPOS (-CLNOUG*0.999, XS, NXS-1)               
C                                                                               
C---- FIT INTERPOLATING CURVE TO SET OF TIMESTEP MULTIPLIERS.                   
C---- CALCULATE INTERPOLATED/EXTRAPOLATED QS AT EACH OUTBOARD X POSN.           
C                                                                               
      CALL FITTER (NQS,CQS(1,1),CQS(1,2),                                       
     >             NQXSO+NQXSI,QXS(1-NQXSO),QS(1-NQXSO),'LINEAR')               
C                                                                               
C---- SET UP TIME POINTS FOR RECORDING ION DISTRIBUTIONS IF IMPULSE             
C---- MODE REQUIRED.  ARRAY "CTIMES" WILL HOLD THE "ITERATION NOS AT            
C---- WHICH THE POSITION IS TO BE RECORDED" FOR EACH CHARGE STATE.              
C---- ARRAY "TIMES" WILL HOLD THE ACTUAL TIMEPOINTS (IN SECONDS).               
C---- PRIMARY NEUTRALS WILL BE TREATED AS FOR TOTAL NEUTRALS.                   
C---- SET "CSTMAX" TO ITERATION CUTOFF POINT, EITHER 10 SECONDS (FOR            
C---- STEADY STATE CASE), OR MAX OF GIVEN DWELL TIME * FACTORS                  
C---- TOTALLY IGNORE VARIABLE TIMESTEPS AT THIS STAGE ...                       
C                                                                               
C     CTMAX IS NOW READ IN AS A VARIABLE - USUALLY 10 SECONDS
C     THIS IS AN EFFORT TO SAVE CPU TIME IN CASES WHERE WE ARE 
C     ONLY INTERESTED IN THE BEHAVIOUR OF THE IONS NEAR THEIR 
C     INJECTION POINT.
C
C     DAVID ELDER  1990 JAN 5
C
      CSTMAX = 0                                                                
C
      IF (IMODE.NE.2) THEN                                                      
        DO 190 IT = 1, NTS                                                      
          CTIMES(IT,0) = DWELTS(0)*DWELFS(IT)/FSRATE                            
          CSTMAX = MAX (CSTMAX, DWELTS(0)*DWELFS(IT)/QTIM)                      
          IF (NIZS.GT.0) THEN                                                   
            DO 185 IZ = 1, NIZS                                                 
              CTIMES(IT,IZ) = DWELTS(IZ)*DWELFS(IT)/QTIM                        
              CSTMAX = MAX (CSTMAX, CTIMES(IT,IZ))                              
  185       CONTINUE                                                            
          ENDIF                                                                 
  190   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C---- TAKE MAX OF THIS VALUE WITH 10 SECONDS, IF STEADY STATE REQUIRED          
C                                                                               
      IF (IMODE.NE.1) CSTMAX = MAX (CSTMAX, CTMAX/QTIM)                        
C                                                                               
C---- INTRODUCE AN EXTRA TIMEPOINT WHICH WILL NEVER BE REACHED                  
C---- THIS TIMEPOINT WILL APPLY EQUALLY TO NON-IMPULSE MODE CASES               
C                                                                               
c     NOTE: when using ultra short ion cut off times this value can be 
c           too small for the neutrals - so allow for upper edge of time
c           bin to be large and use proper cut off for stopping following 
c           the ions 
c
      CTIMES(NTS+1,0) = 2.0 * max(CSTMAX,10.0/qtim) * QTIM/FSRATE                              
      IF (NIZS.GT.0) THEN                                                       
        DO 195 IZ = 1, NIZS                                                     
          CTIMES(NTS+1,IZ) = 2.0 * CSTMAX                                       
  195   CONTINUE                                                                
      ENDIF                                                                     
c
c     jdemod - note that NTS must be increased to accomodate this value - NTS
c              should not be zero so that there is always at least one
c              time point - accessing and saving time dependent data should
c              be based on the simulation mode and not the value of NTS.
c
      nts = nts+1
C                                                                               
C-----------------------------------------------------------------------        
C     SET JOB TO INCLUDE DATAFILE,TIME AND DATE TO PROVIDE UNIQUE               
C     REFERENCE FOR THIS RUN.                                                   
C     THE REFERENCE WILL BE PRINTED AS TWO CHAR STRINGS OF LENGTH 36.           
C-----------------------------------------------------------------------        
C                                                                               
      CALL ZA08AS (SYSTIM)                                                      
      CALL ZA09AS (SYSDAT)                                                      
      CALL ZV01AD (5, VSN, DSN, JFCB)                                           
      IF (JFCB(45:52).NE.' ') THEN                                              
        DO 220 NC = 1, 40                                                       
          IF (JFCB(NC:NC).EQ.' ') THEN                                          
            JFCB(NC:NC) = '('                                                   
            DO 210 NM = 1, 8                                                    
              IF (JFCB(44+NM:44+NM).NE.' ') THEN                                
                JFCB(NC+NM:NC+NM) = JFCB(44+NM:44+NM)                           
              ELSE                                                              
                JFCB(NC+NM:NC+NM) = ')'                                         
                GOTO 230                                                        
              ENDIF                                                             
  210       CONTINUE                                                            
            JFCB(NC+9:NC+9) = ')'                                               
            GOTO 230                                                            
          ENDIF                                                                 
  220   CONTINUE                                                                
  230   CONTINUE                                                                
      ENDIF                                                                     
      WRITE (JOB,'(A36,A8,3X,A8,''   LIM'',A5)')                                
     >  JFCB(1:36),SYSDAT,SYSTIM,VERSON                                         
C                                                                               
C-----------------------------------------------------------------------        
C                     PRINT HEADING                                             
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (COMENT,'(''*'',60('' ''),''*'')')                                  
      CALL PRC                                                                  
     >('**************************************************************')        
      CALL PRC (COMENT)                                                         
      WRITE(7,'('' *'',18X,''RUN OF LIM VERSION '',A5,18X,''*'')')VERSON        
      WRITE(7,'('' *'',18X,24(''-''),18X,''*'')')                               
      CALL PRC (COMENT)                                                         
      WRITE (7,'(1X,''* '',A58,'' *'')') TITLE(1:58)                            
      CALL PRC (COMENT)                                                         
      WRITE (7,'(1X,''* '',A58,'' *'')') JOB(1:58)                              
      CALL PRC (COMENT)                                                         
      CALL PRC                                                                  
     >('**************************************************************')        
      CALL PRB                                                                  
      WRITE (6,'(/1X,A,//2X,A,/)') TITLE,JOB                                    
C                                                                               
      IF (IERR.NE.0) GOTO 1002                                                  
C                                                                               
C-----------------------------------------------------------------------        
C  CALCULATE RANDOM NUMBERS SEED FROM DATE AND TIME.                            
C  SURAND REQUIRES A NUMBER IN RANGE 1 < SEED < 2**31-1                         
C  USING SYSTIM AND SYSDAT, EG 23:59.59 AND 09/29/99, WE ADD 1 TO EACH          
C  DIGIT, AND MULTIPLY THEM TOGETHER, GIVING MAXIMUM OF 1296000000.             
C  IF CISEED GIVEN >0, USE THIS INSTEAD AS THE RANDOM NUMBER SEED.              
C-----------------------------------------------------------------------        
C                                                                               
      IF (CISEED.LE.0) THEN                                                     
        CISEED = 1                                                              
        DO 300 J = 1, 8                                                         
          IF (J.EQ.3.OR.J.EQ.6) GOTO 300                                        
          CISEED = CISEED * (ICHAR(SYSTIM(J:J)) - ICHAR('0') + 1) *             
     >                      (ICHAR(SYSDAT(J:J)) - ICHAR('0') + 1)               
  300   CONTINUE                                                                
        CISEED = CISEED - 1                                                     
      ENDIF                                                                     
      SEED = DBLE (REAL (CISEED))                                               
C                                                                               
C---- SET RANDOM NUMBER SEED VIA RANINI DUMMY ROUTINE IF REQUIRED               
C---- FOR NON-IBM APPLICATIONS.                                                 
C---- GET FIRST RANDOM NUMBER IN SEQUENCE (AND DISCARD IT).                     
C                                                                               
      CALL RANINI (CISEED)                                                      
      CALL SURAND (SEED, 1, RAN)                                                
      NRAND = 1                                                                 
C                                                                               
C-----------------------------------------------------------------------        
C     INITIALISE NOCORONA PACKAGE WITH 1 IMPURITY & OUTPUT TO CHANNEL 6         
C-----------------------------------------------------------------------        
C                                                                               
      if (cdatopt.eq.0) then 
         PMASS(1) = CRMI                                                           
         CALL AATOM (PMASS, 1, 6, KFAIL)                                           
         IF (KFAIL(1).NE.0) THEN                                                   
           IF (CIOPTA.GT.2) THEN 
             WRITE (6,9100) CRMI                                                   
             STOP                                                                  
           ELSE
             WRITE (6,9100) CRMI
             WRITE (6,9101)   
           ENDIF
         ENDIF                                                                     
 9100    FORMAT(//1X,'ELEMENT WITH MASS ',G11.4,' IS NOT INCLUDED',                
     >    ' IN THE NOCORONA PACKAGE.',/)                                        
 9101    FORMAT(/1X,'EXECUTION CONTINUES: LINE AND POWER RADIATION',
     >    ' ARE INCORRECT.',/)
      endif
C                                                                               
C-----------------------------------------------------------------------        
C  CALL LIM TO CALCULATE IMPURITY LEVELS                                        
C-----------------------------------------------------------------------        
C                                                                               
      ITER = 1                                                                  
      REWIND (NOUT)                                                             
  500 CALL LIM3 (IMODE,NIZS,NIMPS,IMPADD,IMPCF,
     >           QTIM,CPULIM,PRINPS,              
     >           XWIDM,YWIDM,FSRATE,IONTIM,NEUTIM,SEED,IGEOM,                   
     >           NTBS,NTIBS,NNBS,NYMFS,NCVS,
c slmod
     >           FACTA,FACTB,ITER,DEFACT,NRAND,TITLE)     
c     >           FACTA,FACTB,ITER,DEFACT,NRAND)     
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C   CALCULATE PARTICULAR LINE RADIATION                                         
C-----------------------------------------------------------------------        
C                                                                               
C
      CALL PLRP (NIZS,PLAMS,PIZS,NLS,CION)                                      
c slmod begin
c
c This isn't the place for this but I can tidy things up later.
c
c      WRITE(63,*) ' '
c      WRITE(63,'(A)') 'PLRP integration:'
c      WRITE(63,*) ' '
c      WRITE(63,'(2A4,4A14)') 
c     +  'IL','IZ','PL','INTEG PLRP','TOT VOL','TOT NUM'
c
      WRITE(7,*) ' '
      WRITE(7,'(A)') 'PLRP volume integration (3D):'
      WRITE(7,*) ' '
      WRITE(7,'(A18,A7,2A17)') 
     +  'Ionisation State','Line','Integrated PLRP'
      WRITE(7,'(A18,A7,2A17)') 
     +  ' ',' (nm)','  (photons/s)  '
     
      DO IL = 2, NLS

        NUMSUM = 0.0
        VOLSUM = 0.0

        DO IY = -NY3D, NY3D

          IF (IY.NE.0) THEN

            DO IP = -MAXNPS+1, MAXNPS-1 
     
              DO IX = 1, NXS

                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY)) 
                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY)) * PWIDS(IP)
                NUMSUM = NUMSUM + PLRP3(IX,IY,IL,IP) * TMPVOL
                VOLSUM = VOLSUM + TMPVOL

              ENDDO 
            ENDDO

          ENDIF
        ENDDO           
       
        WRITE(63,'(2I4,4G14.5)')
     +    IL,PIZS(IL),PLAMS(IL),NUMSUM/VOLSUM,VOLSUM,NUMSUM

        WRITE(7,'(I18,F7.1,2G17.5)')
     +    PIZS(IL),PLAMS(IL),NUMSUM

      ENDDO

      WRITE(0,*) 'DEBUG: Integrating PLRP'

c
c     PLRPS instead of PLRP3 
c
      WRITE(7,*) ' '
      WRITE(7,'(A)') 'PLRP volume integration (2D):'
      WRITE(7,*) ' '
      WRITE(7,'(A18,A7,2A17)') 
     +  'Ionisation State','Line','Integrated PLRP'
      WRITE(7,'(A18,A7,2A17)') 
     +  ' ',' (nm)','  (photons/s)  '

      DO IL = 2, NLS

        NUMSUM = 0.0
        VOLSUM = 0.0

        DO IY = -NYS/2, NYS/2

          IF (IY.NE.0) THEN

             DO IX = 1, NXS

                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                NUMSUM = NUMSUM + PLRPS(IX,IY,IL) * TMPVOL
                VOLSUM = VOLSUM + TMPVOL

             ENDDO 

          ENDIF
        ENDDO           
       
        WRITE(63,'(2I4,4G14.5)')
     +    IL,PIZS(IL),PLAMS(IL),NUMSUM/VOLSUM,VOLSUM,NUMSUM

        WRITE(7,'(I18,F7.1,2G17.5)')
     +    PIZS(IL),PLAMS(IL),NUMSUM

      ENDDO



      WRITE(0,*) 'DEBUG: PLRP integration complete'

c slmod end
C
C
C----------------------------------------------------------------------
C     CALCULATE SPECTROSCOPIC TEMPERATURES 
C----------------------------------------------------------------------
C
C     NOTE: THE CALL TO THIS ROUTINE COULD BE PLACED INSIDE PLRP
C           THE SOURCE IS IN THE FILE WITH PLRP 
C
      WRITE(0,*) 'DEBUG: Calling SPECTEMP'
      CALL SPECTEMP (PLAMS,PIZS,NLS)                                  
C                                            
C-----------------------------------------------------------------------        
C   DUMP RESULTS IN AN EXTERNAL FILE                                            
C-----------------------------------------------------------------------        
C                                                                               
      WRITE(0,*) 'DEBUG: Dumping results'
      CALL DMPOUT (TITLE,NIZS,NOUT,IERR,JOB,IMODE,PLAMS,PIZS,NLS,              
     >           FACTA,FACTB,ITER,NITERS)                                       
      IF (IERR.NE.0) GOTO 1003                                                  
C                                                                               
C-----------------------------------------------------------------------        
C  CHECK FOR FURTHER ITERATIONS FOR SELF-CONSISTENT PLASMA                      
C-----------------------------------------------------------------------        
C                                                                               
      WRITE(0,*) 'DEBUG: Checking for further iterations'
      IF (ITER.LT.NITERS .AND. DEFACT.GT.0.0D0) THEN                            
        ITER = ITER + 1                                                         
        TITLE(61:80) = NUMBER(ITER) // ' ITERATION      '                       
        WRITE (COMENT,'(''*'',60('' ''),''*'')')                                
        TOTTIM = ZA02AS (1) - STATIM                                            
        CALL PRI ('TIME USED SO FAR ...    (S)   ',NINT(TOTTIM))                
        CALL PRB                                                                
        CALL PRB                                                                
        CALL PRC                                                                
     >('**************************************************************')        
        CALL PRC (COMENT)                                                       
        WRITE (7,'(1X,''* '',17X,A20,21X,'' *'')') TITLE(61:80)                 
        CALL PRC (COMENT)                                                       
        CALL PRC                                                                
     >('**************************************************************')        
        CALL PRB                                                                
        WRITE (6,'(//1X,130(''*''),//10X,A,//)') TITLE(61:80)                   
        CNEUTA = KNEUTA                                                         
        CNEUTB = KNEUTB                                                         
        CNEUTC = KNEUTC                                                         
        CNEUTE = KNEUTE                                                         
C
C       ANY EXTRA IONS FOR CROSS-FIELD FLUXES MUST BE TAKEN FROM 
C       NIMPS SO THAT A CALL TO NEUT CAN ADD THEM IN AGAIN WITHOUT
C       CAUSING A CASCADE.
C
        NIMPS = NIMPS - IMPADD - IMPCF
C
        CALL TAUIN3 (QTIM,NIZS,DEFACT)                                          
        GOTO 500                                                                
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C                                                                               
      WRITE(0,*) 'DEBUG: Cleaning up'
      TOTTIM = ZA02AS (1) - STATIM                                              
      WRITE (7,'('' TOTAL RANDOM NUMBERS USED'',I12)') NRAND                    
      CALL PRI ('TIME FOLLOWING NEUTRALS (S)   ',NINT(NEUTIM))                  
      CALL PRI ('TIME FOLLOWING IONS     (S)   ',NINT(IONTIM))                  
      CALL PRI ('TOTAL CPU TIME USED     (S)   ',NINT(TOTTIM))                  
      CALL PRB
      CALL PRB   
      WRITE(0,*) 'DEBUG: Getting data and time'
      CALL ZA08AS(SYSTIM)
      CALL ZA09AS(SYSDAT)
      WRITE(0,*) 'DEBUG: Last bit'
      WRITE (COMENT,'(''TIME AT END OF RUN : '',A8,3X,A8)') 
     >                SYSTIM,SYSDAT                         
      CALL PRC(COMENT)
      CALL PRB
      WRITE (6,'('' TOTAL RANDOM NUMBERS USED'',I12)') NRAND                    
      WRITE (6,'('' TIME FOLLOWING NEUTRALS (S)   '',G11.4)') NEUTIM            
      WRITE (6,'('' TIME FOLLOWING IONS     (S)   '',G11.4)') IONTIM            
      WRITE (6,'('' TOTAL CPU TIME USED     (S)   '',G11.4)') TOTTIM            
      WRITE(0,*) 'DEBUG: Done'
c
c     Deallocate dynamic storage
c
      call deallocate_dynamic_storage

      STOP                                                                      
C                                                                               
 1002 CALL PRC ('RUNLIM3: ERROR OCCURED DURING DATA INPUT - ABORTED')           
c
c     Deallocate dynamic storage
c
      call deallocate_dynamic_storage
      STOP                                                                      

 1003 CALL PRC ('RUNLIM3: ERROR OCCURED DURING DUMP. RESULTS NOT SAVED')        
c
c     Deallocate dynamic storage
c
      call deallocate_dynamic_storage
      STOP                                                                      
      END                                                                       
c
c
c
      subroutine allocate_dynamic_storage
      ! routine to allocate dynamic storage at fixed sizes - eventually update to allow dynamic size definitions
      use mod_cadas
      use mod_cadas2
      use mod_cneut
      use mod_cnoco
      use mod_commv
      use mod_comnet
      use mod_comt2
      use mod_comtau
      use mod_comtor
      use mod_comxyt
      use mod_coords
      use mod_crand
      use mod_cyield
      use mod_dynam1

      use mod_dynam3

      use mod_global_options
      use mod_printr
      use mod_save
      use mod_slcom
      use mod_zommv

      implicit none

      ! LIM


      call allocate_mod_cadas
      call allocate_mod_cadas2
      call allocate_mod_cneut
      call allocate_mod_cnoco
      call allocate_mod_commv
      call allocate_mod_comnet
      call allocate_mod_comt2
      call allocate_mod_comtau
      call allocate_mod_comtor
      call allocate_mod_comxyt
      call allocate_mod_coords
      call allocate_mod_crand
      call allocate_mod_cyield
      call allocate_mod_dynam1
      call allocate_mod_dynam3
      call allocate_mod_global_options
      call allocate_mod_printr
      call allocate_mod_save
      call allocate_mod_slcom
      !call allocate_mod_unstructured
      call allocate_mod_zommv

      
      return
      end
c
c
c     
      subroutine deallocate_dynamic_storage
      use mod_cadas
      use mod_cadas2
      use mod_cneut
      use mod_cnoco
      use mod_commv
      use mod_comnet
      use mod_comt2
      use mod_comtau
      use mod_comtor
      use mod_comxyt
      use mod_coords
      use mod_crand
      use mod_cyield
      use mod_dynam1
      use mod_dynam3
      use mod_global_options
      use mod_printr
      use mod_save
      use mod_slcom
      use mod_zommv
      implicit none

      ! LIM
      write(0,*) '1'
      call deallocate_mod_cadas
      write(0,*) '2'
      call deallocate_mod_cadas2
      write(0,*) '3'
      call deallocate_mod_cneut
      write(0,*) '4'
      call deallocate_mod_cnoco
      write(0,*) '5'
      call deallocate_mod_commv
      write(0,*) '6'
      call deallocate_mod_comnet
      write(0,*) '7'
      call deallocate_mod_comt2
      write(0,*) '8'
      call deallocate_mod_comtau
      write(0,*) '9'
      call deallocate_mod_comtor
      write(0,*) '10'
      call deallocate_mod_comxyt
      write(0,*) '11'
      call deallocate_mod_coords
      write(0,*) '12'
      call deallocate_mod_crand
      write(0,*) '13'
      call deallocate_mod_cyield
      write(0,*) '14'
      call deallocate_mod_dynam1
      write(0,*) '15'
      call deallocate_mod_dynam3
      write(0,*) '16'
      call deallocate_mod_global_options
      write(0,*) '17'
      call deallocate_mod_printr
      write(0,*) '18'
      call deallocate_mod_save
      write(0,*) '19'
      call deallocate_mod_slcom
      write(0,*) '20'
      !call deallocate_mod_unstructured
      call deallocate_mod_zommv
      write(0,*) '21'


      return
      end
