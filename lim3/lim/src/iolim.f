c     -*-Fortran-*-
c        
      SUBROUTINE READIN (TITLE,IGEOM,IMODE,NIZS,NIMPS,IMPADD,
     >                   FSRATE,QTIM,CPULIM,IERR,NTBS,NTIBS,NNBS,
     >                   NYMFS,NCVS,NQS,NITERS)                                 
      IMPLICIT  none
      INTEGER   IERR,IGEOM,IMODE,NIZS,NIMPS,NTBS,NTIBS,NNBS,NYMFS           
      INTEGER   IMPADD
      INTEGER   NQS,NITERS,NCVS
      REAL      FSRATE,QTIM,CPULIM                                              
      CHARACTER TITLE*80                                                        
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  READIN:   READS IN THE DATAFILE, PERFORMS VALIDITY CHECKS, ETC.  *        
C  *                                                                   *        
C  *                                                                   *        
C  *  CHRIS FARRELL    FEBRUARY 1988                                   *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comtau'                                                          
C     INCLUDE (COMTAU)                                                          
      INCLUDE 'coords'                                                          
C     INCLUDE (COORDS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
c
c slmod begin
      INCLUDE 'slcom'
c slmog end
      include 'cadas'
c
      INTEGER NDS,ISTEP,NPLANE,JERR                                             
      LOGICAL RLLIM
c slmod
      INTEGER II
c slmod end
C                                                                               
c slmod
c      WRITE(0,*) 'Begin READIN'
c slmod end
c
c jdemod - make sure unstructured input is initialized prior to reading in the input file
c
      call InitializeUnstructuredInput
c
      CALL RDC (TITLE, 'TITLE FOR RUN', IERR)                                   
      call rdi (cdatopt,.true.,0,.true.,1, 'Rad/ioniz data source',ierr)
      call rdc (useridh,'ADAS H userid',ierr)
      call rdi (iyearh,.true., 0,.true.,99,'ADAS H year          ',ierr)
      call rdc (useridz,'ADAS Z userid',ierr)
      call rdi (iyearz,.true., 0,.true.,99,'ADAS Z year          ',ierr)
      CALL RDI (CIOPTA,.TRUE., 0,.TRUE., 6,'IONIZATION OPT       ',IERR)        
      CALL RDI (CIOPTB,.TRUE., 0,.TRUE., 13,'COLLISION OPT       ',IERR)        
      CALL RDI (CIOPTC,.TRUE., 0,.TRUE., 3,'FRICTION OPT         ',IERR)        
      CALL RDI (CIOPTD,.TRUE., 0,.TRUE., 4,'HEATING OPT          ',IERR)        
c slmod
      CALL RDI (CIOPTE,.TRUE., 0,.TRUE.,13,'INJECTION OPT        ',IERR)        
c      CALL RDI (CIOPTE,.TRUE., 0,.TRUE., 9,'INJECTION OPT        ',IERR)        
c slmod end
      CALL RDI (CIOPTF,.TRUE., 0,.TRUE., 9,'SOL OPT              ',IERR)        
      CALL RDI (CIOPTG,.TRUE., 0,.TRUE., 6,'PLASMA DECAY OPT     ',IERR)        
      CALL RDI (CIOPTK,.TRUE.,-1,.TRUE., 6,'PLASMA ION TEMP OPT  ',IERR)
      CALL RDI (CIOPTL,.TRUE., 0,.TRUE., 1,'TEB GRAD COEFF OPT   ',IERR)
      CALL RDI (CIOPTM,.TRUE., 0,.TRUE., 1,'TIB GRAD COEFF OPT   ',IERR)
      CALL RDI (CIOPTH,.TRUE., 0,.TRUE.,11,'LIMITER EDGE OPT     ',IERR)        
      CALL RDI (CIOPTI,.TRUE., 0,.TRUE., 2,'CX RECOMB OPT        ',IERR)        
      CALL RDI (CDIFOP,.TRUE., 0,.TRUE., 2,'FIRST DIFFUSE OPT    ',IERR)        
      CALL RDI (CIOPTN,.TRUE., 0,.TRUE., 1,'DIFFUSION TYPE OPTION',IERR)
      CALL RDI (CDPERP,.TRUE., 0,.TRUE., 1,'DPERP OPTION         ',IERR)
      CALL RDI (CVPOPT,.TRUE., 0,.TRUE., 2,'V PINCH OPTION       ',IERR)
      CALL RDI (CNEUTA,.TRUE., 0,.TRUE., 2,'CONTROL SWITCH       ',IERR)        
c slmod
      CALL RDI (CNEUTB,.TRUE., 0,.TRUE.,10,'LAUNCH OPTION        ',IERR)        
c      CALL RDI (CNEUTB,.TRUE., 0,.TRUE., 8,'LAUNCH OPTION        ',IERR)        
c slmod end
      CALL RDI (CNEUTC,.TRUE., 0,.TRUE.,17,'VEL/ANGLE FLAG       ',IERR)        
      CALL RDI (NVAOPT,.TRUE.,-1,.TRUE.,17,'NEUT VEL/ANGLE FLAG  ',IERR)        
      CALL RDI (CNEUTD,.TRUE., 0,.TRUE., 8,'SPUTTER OPTION       ',IERR)        
      CALL RDI (CNEUTE,.TRUE., 0,.TRUE., 2,'NORMAL OPTION        ',IERR)        
      CALL RDI (CNEUTF,.TRUE., 0,.TRUE., 1,'NEUT SPREADING       ',IERR)        
      CALL RDI (NRFOPT,.TRUE., 0,.TRUE., 1,'NEUT REFLECT OPT     ',IERR)        
      CALL RDI (IGEOM ,.TRUE., 0,.TRUE., 2,'RADIAL GEOM OPT      ',IERR)        

C     POLOIDAL EXTENT OPTION (FOR 3D)
      CALL RDI (CIOPTJ,.TRUE., 0,.TRUE., 1,'POLOIDAL SIZE OPT    ',IERR)
C     READ IN LIMITER POLOIDAL SIZE 
      CALL RDR(CPCO,.TRUE.,0.0,.FALSE.,0.0,'LIMITER POL. EXTENT',  IERR)     
C                                                                               
      CALL RDR (CA,   .TRUE., 0.0,.FALSE.,0.0,'TORUS CENTRE A',    IERR)        
      CALL RDR (CAW,  .FALSE.,0.0,.TRUE., 0.0,'TORUS WALL AW',     IERR)        
      CALL RDR (CL,   .TRUE.,0.01,.FALSE.,0.0,'LIMITER SEP L',     IERR)        
      CALL RDR (CRMB, .TRUE., 0.1,.FALSE.,0.0,'PLASMA ION MASS',   IERR)        
      CALL RDI (CIZB, .TRUE.,  1 ,.FALSE., 0 ,'PLASMA ION CHARGE', IERR)        
C                                                                               
      CALL RDR(CTBOUL,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TBOUT<',       IERR)        
      CALL RDR(CLTOUL,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY LTOUT<', IERR)        
      CALL RDR(CTBOUG,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TBOUT>',       IERR)        
      CALL RDR(CLTOUG,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY LTOUT>', IERR)        
      CALL RDR(CTBIN, .TRUE. ,0.0,.FALSE.,0.0,'TEMP TBIN',         IERR)        
      CALL RDR(CLTIN1,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY LTIN1',  IERR)        
      CALL RDR(CGTIN1,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY GTIN1',  IERR)        
      CALL RDR(CATIN, .TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY ATIN',   IERR)        
      CALL RDR(CLTIN2,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY LTIN2',  IERR)        
      CALL RDR(CGTIN2,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY GTIN2',  IERR)        

      CALL RDRARN(CTBINS,NTBS,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                      1,'SET OF X,TB VALUES',JERR)        
      IF (JERR.NE.0) GOTO 1001                                                  
C
C     READ IN ELECTRON TEMPERATURE GRADIENT INFORMATION, IF ANY
C     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
C     OF THE TEMPERATURE
C
      CALL RDRARN(TMEG,NTEG,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                     1,'SET OF Y,TME VALUES',JERR)        
C
C     READ IN BACKGROUND ION TEMPERATURE CHARACTERISTICS. THERE ARE 
C     CURRENTLY THE SAME OPTIONS ALLOWED FOR BOTH ELECTRON AND ION 
C     TEMPERATURES. THE OPTIONS ARE SELECTED INDEPENDENTLY AND THE 
C     PARAMETERS PASSED TO THE PLASMA SUBROUTINE.
C
C     DAVID ELDER , FEB 5 ,1990
C
      CALL RDR(CTIBOUL,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TIBOUT<',     IERR)       
      CALL RDR(CLTIOUL,.TRUE. ,0.0,.FALSE.,0.0,'TEM DECAY LTIOUT<',IERR)       
      CALL RDR(CTIBOUG,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TIBOUT>',     IERR)       
      CALL RDR(CLTIOUG,.TRUE. ,0.0,.FALSE.,0.0,'TEM DECAY LTIOUT>',IERR)       
      CALL RDR(CTIBIN, .TRUE. ,0.0,.FALSE.,0.0,'TEMP TIBIN',       IERR)       
      CALL RDR(CLTIIN1,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY LTIIN1',IERR)       
      CALL RDR(CGTIIN1,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY GTIIN1',IERR)       
      CALL RDR(CATIIN, .TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY ATIIN', IERR)       
      CALL RDR(CLTIIN2,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY LTIIN2',IERR)       
      CALL RDR(CGTIIN2,.TRUE. ,0.0,.FALSE.,0.0,'TEMP DECAY GTIIN2',IERR)      

      CALL RDRARN(CTIBINS,NTIBS,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,         
     >                                      1,'SET OF X,TB VALUES',JERR)        
      IF (JERR.NE.0) GOTO 1001                                                  
C
C     READ IN ION TEMPERATURE GRADIENT INFORMATION, IF ANY
C     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
C     OF THE TEMPERATURE
C
      CALL RDRARN(TMIG,NTIG,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                     1,'SET OF Y,TMI VALUES',JERR)        
C                                                                               
      CALL RDR(CNBOUL,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY NBOUT<',    IERR)        
      CALL RDR(CLNOUL,.TRUE. ,0.0,.FALSE.,0.0,'DENS DECAY LNOUT<', IERR)        
      CALL RDR(CNBOUG,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY NBOUT>',    IERR)        
      CALL RDR(CLNOUG,.TRUE. ,0.0,.FALSE.,0.0,'DENS DECAY LNOUT>', IERR)        
      CALL RDR(CNBIN, .TRUE. ,0.0,.FALSE.,0.0,'DENSITY NBIN',      IERR)        
      CALL RDR(CLNIN1,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY DECAYLNIN1',IERR)        
      CALL RDR(CGNIN1,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY DECAYGNIN1',IERR)        
      CALL RDR(CANIN, .TRUE. ,0.0,.FALSE.,0.0,'DENSITY DECAY ANIN',IERR)        
      CALL RDR(CLNIN2,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY DECAYLNIN2',IERR)        
      CALL RDR(CGNIN2,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY DECAYGNIN2',IERR)        
      CALL RDR(CNBA  ,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY NBA       ',IERR)        
      CALL RDR(CGAMMA,.TRUE. ,0.0,.FALSE.,0.0,'DENSITY GAMMA     ',IERR)        
      CALL RDRARN(CNBINS,NNBS,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                      1,'SET OF X,NB VALUES',JERR)        
      IF (JERR.NE.0) GOTO 1001                                                  
C                                                                               
      CALL RDR(CVIN,  .FALSE.,0.0,.FALSE.,0.0,'PINCH0,2 PARAMETER',IERR)        
      CALL RDR(VPV0,  .FALSE.,0.0,.FALSE.,0.0,'PINCH1 VELOCITY',   IERR)        
      CALL RDR(VPALPH,.FALSE.,0.0,.FALSE.,0.0,'PINCH1 EXP     ',   IERR)        
      CALL RDR(CVPCUT,.TRUE., 0.0,.FALSE.,0.0,'PINCH2 CUTOFF PT',  IERR)        
      CALL RDR(CVPOUT,.FALSE.,0.0,.FALSE.,0.0,'ARBITRARY V PINCH', IERR) 
      CALL RDR(CVXMIN,.FALSE.,0.0,.FALSE.,0.0,'MIN X for ARB V',   IERR) 
      CALL RDR(CVXMAX,.FALSE.,0.0,.FALSE.,0.0,'MAX X for ARB V',   IERR) 
C
      CALL RDR(CRDXO(1),.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERO<    ',IERR)        
      CALL RDR(CRDXO(2),.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERO>    ',IERR)        
      CALL RDR(CRDXI(1),.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERI<    ',IERR)        
      CALL RDR(CRDXI(2),.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERI>    ',IERR)        
      CALL RDR(CRDXO(3),.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPCO(CENT)',IERR)        
      CALL RDR(CRDXI(3),.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPCI(CENT)',IERR)        
      CALL RDR(CRDD,  .FALSE.,0.0,.FALSE.,0.0,'DIFFUS DECAY FACT', IERR)        
      CALL RDR(DPALPH,.FALSE.,0.0,.FALSE.,0.0,'DIFFUS FACT ALPHA', IERR)
      CALL RDR(DPBETA,.FALSE.,0.0,.FALSE.,0.0,'DIFFUS FACT BETA ', IERR)
C
      CALL RDR(CDPOL, .FALSE.,0.0,.FALSE.,0.0,'POLOIDAL DIFFUS. ', IERR)        
      CALL RDR(CVPOL, .FALSE.,0.0,.FALSE.,0.0,'POLOIDAL DRIFT VEL',IERR)
      CALL RDR(CRMI,  .TRUE. ,0.1,.FALSE.,0.0,'IMPURITY ION MASS', IERR)        
      CALL RDI(CION,  .TRUE. , 1 ,.FALSE., 0 ,'IMP ATOMIC NUMBER', IERR)        
c slmod begin - tmp
      IF (CRMI.EQ.28.0.AND.CION.EQ.7.AND.
     >   (CDATOPT.NE.0.OR.NIZS.GT.1)) THEN
        WRITE(0,*) 'Ionisation data options invalid for molecular N'
        STOP
      ENDIF
c slmod end      
      CALL RDR(CEBD,  .TRUE. ,0.0,.FALSE.,0.0,'CHAR ENERGY EBD',   IERR)        
      CALL RDI(CIZEFF,.TRUE. , 1 ,.FALSE., 0 ,'Z EFFECTIVE (SELF)',IERR)        
      CALL RDR(CTGAS, .FALSE.,0.0,.FALSE.,0.0,'GAS TEMPERATURE',   IERR)        
      CALL RDR(CTIMSC,.FALSE.,0.0,.FALSE.,0.0,'INJECTION TIME',    IERR)        
      CALL RDR(CTEMSC,.TRUE. ,0.0,.FALSE.,0.0,'INJECTION ION TEMP',IERR)        
      CALL RDR(CENGSC,.TRUE. ,0.0,.FALSE.,0.0,'INJ NEUTRAL ENERGY',IERR)        
      CALL RDR(CEIN2, .TRUE. ,0.0,.FALSE.,0.0,'ALTERNATE EIN VAL ',IERR)
C     PROBABILITY OF Ein1 (PEin2 = 1-PEin1) FOR LAUNCH6, INJ9
      CALL RDR(CPROB,.TRUE. ,0.0, .TRUE. ,1.0,'PROB. OF CENGSC   ',IERR)
      CALL RDR(CEMAXF,.FALSE.,0.0,.FALSE.,0.0,'EMAX-FACTOR',       IERR)        
      CALL RDR(CCUT,  .FALSE.,0.0,.FALSE.,0.0,'NEUTRAL X CUTOFF',  IERR)        
C                                                                               
      CALL RDR(CXSC,  .TRUE., CAW,.TRUE., CA, 'INITIAL X POSITION',IERR)        
      CALL RDR(CYSC,  .TRUE.,0.0,.TRUE.,CL+CL,'INITIAL Y POSITION',IERR)        
      CALL RDR(CPSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL P POSITION',IERR)        
      CALL RDI(CIZSC, .TRUE.,  1, .TRUE.,CION,'INITIAL IZ STATE',  IERR)        
C                                                                               
      CALL RDR(CNHC,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHC   ',IERR)        
      CALL RDR(CNHO,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHO   ',IERR)        
      CALL RDR(CLAMHX,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHX ',IERR)        
      CALL RDR(CLAMHY,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHY ',IERR)        
      CALL RDR(CVCX  ,.TRUE., 0.0,.FALSE.,0.0,'CXREC CONSTANT VCX',IERR)        

      CALL RDR(CXNEAR,.TRUE.,0.0,.FALSE. ,CA ,'NEAR PARAM  XNEAR ',IERR)        
c
c     Verify CXNEAR  
c
      if (cxnear.gt.ca) then 
         write(0,'(2(a,1x,g12.5,1x))') 
     >           'INPUT ERROR: CXNEAR SPECIFIED AS',cxnear,
     >           'GREATER THAN CA =',ca
         write(0,'(a)') 'CXNEAR SET TO CA'
         write(datunit,'(2(a,1x,g12.5,1x))') 
     >           'INPUT ERROR: CXNEAR SPECIFIED AS',cxnear,
     >           'GREATER THAN CA =',ca
         write(datunit,'(a)') 'CXNEAR SET TO CA'
         cxnear = ca 
      endif
c
      CALL RDR(CYNEAR,.TRUE.,0.0,.FALSE.,CL+CL,'NEAR PARAM YNEAR ',IERR)        
c
c     Verify CYNEAR 
c
      if (cynear.gt.2.0*cl) then 
         write(0,'(2(a,1x,g12.5,1x))') 
     >           'INPUT ERROR: CYNEAR SPECIFIED AS',cynear,
     >           'GREATER THAN 2.0*CL =',2.0*CL
         write(0,'(a)') 'CYNEAR SET TO CL'
         write(datunit,'(2(a,1x,g12.5,1x))') 
     >           'INPUT ERROR: CYNEAR SPECIFIED AS',cynear,
     >           'GREATER THAN 2.0*CL =',2.0*cl
         write(datunit,'(a)') 'CYEAR SET TO CL'
         cynear = cl
      endif 
c
      CALL RDR(CSNORM,.FALSE.,0.0,.FALSE.,0.0,'MEASURE THETA FROM',IERR)        
      CALL RDR(CPLSMA,.FALSE.,0.0,.FALSE.,0.0,'C PLASMA FOR X <= ',IERR)        
      CALL RDR(CVHYIN,.FALSE.,0.0,.FALSE.,0.0,'INB. PLASMA FLOW V',IERR)        
      CALL RDR(CEYIN, .FALSE.,0.0,.FALSE.,0.0,'INBOARD ELEC FIELD',IERR)        
      CALL RDR(CVHOUT,.FALSE.,0.0,.FALSE.,0.0,'OUT. PLASMA FLOW V',IERR)        
      CALL RDR(CEYOUT,.FALSE.,0.0,.FALSE.,0.0,'OUTBRD ELEC FIELD ',IERR)        
      CALL RDR(CZENH, .FALSE.,0.0,.FALSE.,0.0,'COLLIS ENHANC ZENH',IERR)        
      CALL RDI(CIZSET,.FALSE., 0 ,.FALSE., 0 ,'SET TI=TB AT STATE',IERR)        
      CALL RDR(CYSTAG,.TRUE., 0.0,.TRUE., CL ,'SOL6: YSTAG VALUE ',IERR)        
      CALL RDR(CFIMP, .FALSE.,0.0,.FALSE.,0.0,'SELF-SPUTTER FRACT',IERR)        
      CALL RDR(CTRESH,.FALSE.,0.0,.FALSE.,0.0,'SELF-SPU THRESHOLD',IERR)        
      CALL RDR(CSPUMAX,.FALSE.,0.0,.FALSE.,0.0,'SELF-SPU MAX FRAC',IERR)  
      CALL RDI(CBOMBZ,.TRUE.,  0 ,.FALSE., 0 ,'BOMBION CHARGE    ',IERR)        
      CALL RDI(CBOMBF,.TRUE.,  0 ,.TRUE.,  7 ,'BOMBION FLAG 0:7  ',IERR)        
      CALL RDR(CIRF  ,.TRUE., 0.0,.FALSE.,0.0,'IONISE RATE FACTOR',IERR)        
      CALL RDR(CSEF  ,.TRUE., 0.0,.FALSE.,0.0,'SPUT ENHANCE FACT.',IERR)        
      CALL RDR(CSOLEF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE FACT.',IERR)        
C                                                                               
      CALL RDI(IMODE, .TRUE.,  0, .TRUE.,  2, 'LIM MODE',          IERR)        
      CALL RDI(NIZS,  .TRUE.,  0, .TRUE.,MAXIZS,'MAX IZ STATE',    IERR)        
      CALL RDI(NQXSO, .TRUE.,  1, .TRUE.,MAXQXS,'NO OF QXS OUTBRD',IERR)        
      CALL RDI(NQXSI, .TRUE.,  1, .TRUE.,MAXQXS,'NO OF QXS INBRD ',IERR)        
      CALL RDI(NQYS,  .TRUE.,  1, .TRUE.,MAXQYS,'NO OF QYS POSNS ',IERR)        
      CALL RDI(NIMPS, .TRUE.,  1, .TRUE.,MAXIMP,'NO OF IONS',      IERR)        
      CALL RDR(FSRATE,.TRUE.,1.E-15,.TRUE.,1.0,'NEUT ITERATE TIME',IERR)        
      CALL RDR(QTIM,  .TRUE.,1.E-10,.TRUE.,1.0,'LIM ITERATE TIME', IERR)        
C            
C     EXTRA IONS FOR SIMULATING CROSS-FIELD EFFECTS
C
      CALL RDI(IMPADD,.TRUE.,  0, .TRUE.,MAXIMP-NIMPS,'NO OF EXTRA IONS'
     >          ,IERR)        
      CALL RDR(YCFADD,.TRUE.,0.0,.TRUE.,CL,'LIMIT EXTRA ION REL.' ,IERR)
      CALL RDI(CEXNEUT,.TRUE., 0, .TRUE.,1,'EX.NEUT LAUNCH OPT. ' ,IERR)        
C
C     CROSS-FIELD BACKGROUND FLUX FACTOR FOR EXTRA BACKGROUND ION 
C     SPUTTERING
C
      CALL RDR(CFBGFF,.TRUE.,0.0,.FALSE.,0.0,'CROSS-FIELD FLUX FACT.'
     >         ,IERR)   
C
      CALL RDR(WEDGAN(1),.TRUE. , 0.0,.TRUE. ,90.0,'WEDGE ANGLE <',IERR)        
      CALL RDR(WEDGAN(2),.TRUE. , 0.0,.TRUE. ,90.0,'WEDGE ANGLE >',IERR)        
      CALL RDR(XL1(1),.TRUE. , CAW,.TRUE.,0.0 , 'BLUNT NOSE XL1 <',IERR)        
      CALL RDR(YL1(1),.TRUE. , 0.0,.TRUE.,CL+CL,'BLUNT NOSE YL1 <',IERR)        
      CALL RDR(XL2(1),.TRUE. , CAW,.TRUE.,0.0 , 'BLUNT NOSE XL2 <',IERR)        
      CALL RDR(YL2(1),.TRUE. , 0.0,.TRUE.,CL+CL,'BLUNT NOSE YL2 <',IERR)        
      CALL RDR(XL1(2),.TRUE. , CAW,.TRUE.,0.0 , 'BLUNT NOSE XL1 >',IERR)        
      CALL RDR(YL1(2),.TRUE. , 0.0,.TRUE.,CL+CL,'BLUNT NOSE YL1 >',IERR)        
      CALL RDR(XL2(2),.TRUE. , CAW,.TRUE.,0.0 , 'BLUNT NOSE XL2 >',IERR)        
      CALL RDR(YL2(2),.TRUE. , 0.0,.TRUE.,CL+CL,'BLUNT NOSE YL2 >',IERR)        
C                                                                               
      CALL RDR(X0S,   .FALSE., 0.0,.FALSE.,0.0, 'BOX COORD X0S   ',IERR)        
      CALL RDR(X0L,   .FALSE., 0.0,.FALSE.,0.0, 'BOX COORD X0L   ',IERR)        
      CALL RDR(Y0S,   .FALSE., 0.0,.FALSE.,0.0, 'BOX COORD Y0S   ',IERR)        
      CALL RDR(Y0L,   .FALSE., 0.0,.FALSE.,0.0, 'BOX COORD Y0L   ',IERR)        
C                                                                               
      CALL RDR(CPULIM,.TRUE.,  0.0,.FALSE.,0.0, 'CPU LIMIT',       IERR)        
C                                                                               
      CALL RDR(CPFIR, .TRUE. ,0.0,.FALSE.,0.0,'FIRST +/-P EXTENT ',IERR)        
      CALL RDR(CPSUB, .TRUE. ,0.0,.FALSE.,0.0,'SUBSEQUENT P WIDTH',IERR)        
C                                                                               
C---- READ X VALUES FROM AW TO A (EXCLUSIVE)                                    
C---- READ Y VALUES UP TO L (EXCLUSIVE) ONLY: VALUES FROM L TO 2L WILL          
C---- BE SET SO THAT Y BINS ARE SYMMETRIC ABOUT Y=L                             
C                                                                               
      CALL RDRAR(XS,NXS,MAXNXS-1,  CAW, CA,.TRUE.,'X BIN UPBOUNDS',IERR)        
      CALL RDRAR(YS,NYS,(MAXNYS/2)-1,0.,CL,.TRUE.,'Y BIN UPBOUNDS',IERR)        
C                                                                               
C---- READ IN TIME DEPENDENT DATA  (NOTE 128)                                   
C                                                                               
      CALL RDRAR(DWELTS,NDS,MAXIZS+1,0.0,MACHHI,.FALSE.,
     >           'TAU DWELL',IERR)        
      IF (IMODE.NE.2.AND.NDS.LT.NIZS+1) GOTO 1002                               
      CALL RDRAR(DWELFS,NTS,MAXNTS,  0.0,MACHHI,.TRUE., 
     >           'T FACTORS',IERR)        
C                                                                               
C---- READ IN YIELD MODIFIER FUNCTION AND FLAG                                  
C                                                                               
      CALL RDRARN(CYMFS,NYMFS,MAXINS,-MACHHI,MACHLO,.TRUE.,0.0,MACHHI,            
     >                                      2,'SET OF X,M(X) VALS',IERR)        
      CALL RDI(CYMFLG,.TRUE. ,-2 ,.TRUE. , 0,'YIELD MODIFIER FLAG',IERR)        
C
C---- READ IN Q SPUTTERING PARAMETER MULTIPLIER
C
      CALL RDR(QMULTP,.TRUE., 0.0,.FALSE.,0.0, 'PRIM. Q-MULT SPUT',IERR)        
      CALL RDR(QMULTS,.TRUE., 0.0,.FALSE.,0.0, 'SEC.  Q-MULT SPUT',IERR)        
C                                                                               
C---- READ IN DEBUG OPTS 0:OFF >0 PRINT DIAGNOSTICS EVERY I'TH TIMESTEP         
C---- READ IN RANDOM NUMBER GENERATOR SEED                                      
C                                                                               
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG NEUT ***',IERR)        
      CSTEPN = REAL (ISTEP)                                                     
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG LIM  ***',IERR)        
      CSTEPL = REAL (ISTEP)                                                     
      CALL RDI(ISTEP ,.FALSE., 0,.TRUE. ,MAXT,'***DEBUG TRACKS***',IERR)        
      CSTEPT = ISTEP 
      CALL RDI(CISEED,.FALSE., 0 ,.FALSE., 0 ,'RANDOM NUMBER SEED',IERR)        
      CALL RDI(CPRINT,.FALSE., 0 ,.FALSE., 0 ,'PRINT OPTION (0,1)',IERR)        
      CALL RDR(CTHETB,.FALSE.,0.0,.FALSE.,0.0,'LIMTER GEOM THETAB',IERR)        
C                                                                               
C---- READ IN TIMESTEP MULTIPLICATION FACTORS                                   
C                                                                               
      CALL RDRARN(CQS,NQS,MAXINS,-MACHHI,MACHHI,.TRUE.,0.0,MACHHI,                 
     >                             1,'SET OF TIMESTEP MULTIPLIERS',IERR)        
C                                                                               
C---- READ SPLITTING PLANES                                                     
C                                                                               
      CALL RDRAR(CXSPLS(1),NPLANE,MAXINS,0.0,CA,.TRUE.,'SPLIT PLS',IERR)        
      CXSPLS(0) = 2.0 * CAW                                                     
      CXSPLS(NPLANE+1) = 2.0 * CA                                               
      CALL RDI(CNSPL ,.TRUE. , 1 ,.TRUE., 100,'NUMBER OF SUB-IONS',IERR)        
      CALL RDR(CPRUL ,.TRUE. ,0.0,.TRUE., 1.0,'PROB OF RETENTION ',IERR)        
C                                                                               
      CALL RDR(TC    ,.TRUE. ,0.0,.TRUE.,0.28,'BELT LIMITER "TC" ',IERR)        
      CALL RDR(SC    ,.TRUE.,-.06,.TRUE.,0.22,'BELT LIMITER "SC" ',IERR)        
      CALL RDR(TO    ,.TRUE. ,0.0,.FALSE.,0.0,'BELT LIMITER "TO" ',IERR)        
      CALL RDR(SO    ,.FALSE.,0.0,.TRUE., 0.0,'BELT LIMITER "SO" ',IERR)        
      CALL RDR(TV    ,.FALSE.,0.0,.FALSE.,0.0,'BELT LIMITER "TV" ',IERR)        
      CALL RDR(SV    ,.FALSE.,0.0,.FALSE.,0.0,'BELT LIMITER "SV" ',IERR)        
      CALL RDI(CORECT,.TRUE. , 0 ,.TRUE.,  1 ,'CURVATURE CORRECTN',IERR)        
      IF (SO.NE.SC) GC = ATAN ((TC-TO) / (SO-SC))                               
      RP = SQRT ((TO-TC)*(TO-TC) + (SO-SC)*(SO-SC))                             
      CALL RDR(CLARMR,.FALSE.,0.0,.FALSE.,0.0,'LARMOR RADIUS RL  ',IERR)        
      CALL RDR(CKO   ,.FALSE.,0.0,.FALSE.,0.0,'ELONG OUTBOARD    ',IERR)        
      CALL RDR(CKI   ,.FALSE.,0.0,.FALSE.,0.0,'ELONGATION INBOARD',IERR)        
      CALL RDR(CTBI  ,.FALSE.,0.0,.FALSE.,0.0,'BACK ION TEMP TBI ',IERR)        
      CALL RDI(NITERS,.TRUE. , 1 ,.TRUE. ,20 ,'NO. OF ITERATIONS ',IERR)        
      IF (NITERS.GT.1 .AND. CIOPTB.EQ.3) GOTO 1003                              
c slmod begin
      CALL RDR(CFTCUT,.FALSE.,0.0,.TRUE. ,1.0,'STOP X FIELD TRANS',IERR)        
c      CALL RDR(CFTCUT,.FALSE.,0.0,.TRUE. ,0.0,'STOP X FIELD TRANS',IERR)        
c slmod end
      CALL RDR(CANAL ,.TRUE. ,0.0,.FALSE.,0.0,'ANALYTIC EXTENSION',IERR)        
      CALL RDR(CAP   ,.TRUE. ,0.0,.FALSE.,0.0,'ACT PLAS SURF AREA',IERR)        
      CALL RDR(CWL   ,.TRUE. ,0.0,.FALSE.,0.0,'ACT LIM WETTED LEN',IERR)        
      CALL RDR(CTSUB ,.TRUE. ,0.0,.FALSE.,0.0,'SUBSTRATE TEMP.   ',IERR)        
      CALL RDR(CTMAX ,.TRUE. ,0.0,.FALSE.,0.0,'MAX.ITERATION TIME',IERR)  
      CALL RDR(CXDEEP,.TRUE. ,0.0,.FALSE.,0.0,'SCT CROSSOVER PT. ',IERR)  
C
C     READ IN IONIZATION DATA FOR LAUNCH OPTIONS 7 AND 8 
C                                                                               
      CALL RDRARN(LPDION,CLPD,MAXLPD,-MACHHI,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                      1,'X OR Y,Iprob VALS.',JERR)    
C
C     READ IN ( X(M), Vs(eV) ) DATA VALUES FOR SPUT. OPT 8
C
      CALL RDRARN(CVSA,NCVS,MAXINS,-MACHHI,MACHHI,.TRUE.,0.0,MACHHI,
     >            2,'SET OF X,VS(X) VAL',IERR)
C
C     READ IN STEP LIMITER SHAPE COORDS FOR EDGE OPTION 6
C
      CALL RDR(XST1(1), .TRUE. ,CAW, .TRUE. , 0.0, 'STEP XST1 <', IERR) 
      CALL RDR(YST1(1), .TRUE. ,0.0, .TRUE. ,CL+CL,'STEP YST1 <', IERR) 
      CALL RDR(XST2(1), .TRUE. ,CAW, .TRUE. , 0.0, 'STEP XST2 <', IERR) 
      CALL RDR(YST2(1), .TRUE. ,0.0, .TRUE. ,CL+CL,'STEP YST2 <', IERR) 
      CALL RDR(XST3(1), .TRUE. ,CAW, .TRUE. , 0.0, 'STEP XST3 <', IERR) 
      CALL RDR(YST3(1), .TRUE. ,0.0, .TRUE. ,CL+CL,'STEP YST3 <', IERR) 
      CALL RDR(XST1(2), .TRUE. ,CAW, .TRUE. , 0.0, 'STEP XST1 >', IERR) 
      CALL RDR(YST1(2), .TRUE. ,0.0, .TRUE. ,CL+CL,'STEP YST1 >', IERR) 
      CALL RDR(XST2(2), .TRUE. ,CAW, .TRUE. , 0.0, 'STEP XST2 >', IERR) 
      CALL RDR(YST2(2), .TRUE. ,0.0, .TRUE. ,CL+CL,'STEP YST2 >', IERR) 
      CALL RDR(XST3(2), .TRUE. ,CAW, .TRUE. , 0.0, 'STEP XST3 >', IERR) 
      CALL RDR(YST3(2), .TRUE. ,0.0, .TRUE. ,CL+CL,'STEP YST3 >', IERR) 
C
C     READ IN IQXBRK
C
      CALL RDI(CBRK, .TRUE. , -MAXQXS , .TRUE. , MAXQXS,
     >         'Dpero/Dperi IQXBRK', IERR)
C
C     LIMITER RADIUS FOR EDGE OPTION 7
C
      IF (CIOPTH.EQ.7) THEN 
         RLLIM = .TRUE.
      ELSE 
         RLLIM = .FALSE.
      ENDIF  
      CALL RDR(RLEDGE7,RLLIM,CA,.FALSE.,0.0,'LIMITER CURV.OPT 7',IERR) 
C
C     LIMITER EDGE OPTION 8 - RADIUS OF CURVATURE OF TIP
C
      IF (CIOPTH.EQ.8) THEN 
         RLLIM = .TRUE.
      ELSE 
         RLLIM = .FALSE.
      ENDIF          
      CALL RDR(RLC,.TRUE.,0.0,RLLIM,-CAW,'LIMITER CURV.OPT 8',IERR)        
C
C     MINIMUM VALUE OF ION VELOCITY
C
      CALL RDR(CSVYMIN,.TRUE.,0.0,.FALSE.,0.0, 'ION SVY MINIMUM  ',IERR) 
      CALL RDI(CMAXGENS,.TRUE., 0,.FALSE.,  0, 'MAX. GENERATIONS ',IERR)
      CALL RDI(CDCALC, .TRUE. ,0 ,.TRUE. ,1  , 'CALC D FLUX      ',IERR)
      CALL RDR(CDPSTP, .TRUE.,0.0,.FALSE.,0.0, 'DPERP STEP SIZE  ',IERR)
c slmod
      CALL RDI(optdp , .TRUE.,0  ,.FALSE.,1  , 'DIVIMP ion source',IERR)
      CALL RDI(N2OPT , .TRUE.,0  ,.FALSE.,1  , 'Trac N2 break-up ',IERR)
      CALL RDR(N2FRAC, .TRUE.,0.0,.FALSE.,1.0, 'Fast N fraction  ',IERR)
      CALL RDI(NEOPT,  .TRUE.,1  ,.FALSE.,1  , 'Fast N e-option ',IERR)
      CALL RDR(NENG ,  .TRUE.,0.0,.FALSE.,5.0, 'Fast N energy   ',IERR)
      CALL RDI(NIEOPT, .TRUE.,1  ,.FALSE.,1  , 'N+ e-option  ',IERR)
      CALL RDR(NIENG , .TRUE.,0.0,.FALSE.,5.0, 'N+ energy    ',IERR)
c      CALL RDI(SLOPT,  .TRUE.,0.0,.FALSE.,1  , 'STEVE OPTION     ',IERR)
c      IF (IERR.EQ.1) THEN
c       WRITE(7,*) 'Warning! Special Steve option missing but'//
c    +             ' that is okay.  It is set to zero.'
c        IERR  = 0
        SLOPT = 0
c      ELSE
c        WRITE(0,*) 'SLOPT: ',SLOPT
c      ENDIF
c slmod end
C
C
c slmod
      IF (MAXY3D.LT.2.0*NYS+2) THEN
        WRITE(0,*) 'Warning (READIN): MAXY3D is too small.'
      ENDIF 
     
c      WRITE(0,*) 'Done  READIN'
c slmod end


      RETURN                                                                    
 1001 CALL PRC ('READIN: PLASMA DECAY OPT 3 REQUIRES SETS OF TEMPS   ')         
      CALL PRC ('        AND DENSITIES FOR INTERPOLATION.  THESE ARE ')         
      CALL PRC ('        ENTERED WITH THEIR CORRESPONDING X POSITIONS')         
      CALL PRC ('        SO THAT WE HAVE 2 REAL VALUES PER LINE.  IF ')         
      CALL PRC ('        THIS OPTION NOT USED SET "NO. OF VALS" = 0  ')         
      IERR = 1                                                                  
      RETURN                                                                    
 1002 CALL PRI ('READIN: MUST HAVE AVERAGE DWELL TIMES FOR STATES 0 TO',        
     >  NIZS)                                                                   
      IERR = 1                                                                  
      RETURN                                                                    
 1003 CALL PRC ('READIN: SELF-CONSISTENT PLASMA INCOMPATIBLE WITH COLLIS        
     >ION OPTION 3')                                                            
      IERR = 1                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE PRDATA (NIZS,XSCALO,XSCALI)                                
      use eckstein_2002_yield_data
      use variable_wall
      use iter_bm
C     
      implicit none 
c
      REAL      XSCALO,XSCALI
      INTEGER   NIZS 
C
C                                                                               
C***********************************************************************        
C     THIS ROUTINE PRINTS ALL THE OPTION FLAGS, TORUS ATTRIBUTES ETC            
C     STORED IN COMMONS COMTOR AND COMTAU                                       
C                                                                               
C     C.M.FARRELL    NOVEMBER 1987                                              
C                                                                               
C***********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comtau'                                                          
C     INCLUDE (COMTAU)                                                          
      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
      INCLUDE 'coords'                                                          
C     INCLUDE (COORDS)                                                          
c slmod begin - N2 break
      INCLUDE 'slcom'
c slmod end
      CHARACTER*77 COMENT                                                       
c      DATA RADDEG /57.29577952/                                                 
      INTEGER IX 
      integer iz,icxsc 

C-----------------------------------------------------------------------        
      CALL PRB                                                                  
      CALL PRC  ('TORUS DIMENSIONS')                                            
      WRITE (7,9010) '  X POSITION OF CENTRE         A     (M)    ', CA         

      if (lim_wall_opt.eq.0) then 
         call prc('   WALL OPTION 0: WALL LOCATED A'//
     >            ' CONSTANT DISTANCE FROM LCFS')
       WRITE (7,9010)'  X POSITION OF WALL           AW    (M)    ', CAW        
      elseif (lim_wall_opt.eq.1) then 
         call prc('   WALL OPTION 1: WALL LOCATION'//
     >            ' IS A FUNCTION OF Y')
         call prc('                  WALL IS A STRAIGHT LINE SEGMENT')
         call prc('                  STARTING AT AW AT A SPECIFIED'//
     >            ' Yin')
         call prc('                  REACHING AWMIN AT YHALF (THE')
         call prc('                  HALF WAY POINT BETWEEN LIMITERS')
         call prr('                  Yin    = ',ywall_start)
         call prr('                  AW     = ',caw)
         call prr('                  AW-MIN = ',caw_min)
      endif

      IF (CTHETB.EQ.90.0) THEN                                                  
      WRITE (7,9010) '  HALF LIMITER SEPARATION      L     (M)    ', CL         
      ELSE                                                                      
      WRITE (7,9010) '  HALF LIM SEP (TOROIDAL)      LT    (M)    ',            
     >                                       CAP/(2.0*CWL*CLFACT*CSINTB)        
      WRITE (7,9010) '  HALF LIM SEP (POLOIDAL, X>0) LP+   (M)    ', CL         
      WRITE (7,9010) '  HALF LIM SEP (POLOIDAL, X<0) LP-   (M)    ',            
     >                                              CAP/(2.0*CWL*CLFACT)        
      WRITE (7,9010) '  ACTUAL PLASMA SURFACE AREA   AP    (M*M)  ', CAP        
      WRITE (7,9010) '  ACTUAL LIMITER WETTED LENGTH WL    (M)    ', CWL        
      ENDIF                                                                     
      CALL PRB                                                                  
      CALL PRC  ('PLASMA CHARACTERISTICS')                                      
      CALL PRR  ('  PLASMA ION MASS              MB           ', CRMB)          
      CALL PRI  ('  PLASMA ION CHARGE            ZB           ', CIZB)          
C-----------------------------------------------------------------------        
      IF (CIOPTG.EQ.0) THEN                                                     
       CALL PRR ('  TEMPERATURE OUTBOARD Y < 0   TBOUT<(EV)   ', CTBOUL)        
       CALL PRR ('    OUTBOARD CONST Y > 0       TBOUT>(EV)   ', CTBOUG)        
       CALL PRR ('    INBOARD CONSTANT           TBIN  (EV)   ', CTBIN)         
C                                                                               
      ELSE                                                                      
       CALL PRR ('  TEMPERATURE OUTBOARD Y < 0   TBOUT<(EV)   ', CTBOUL)        
       CALL PRR ('    OUTBOARD EXP DECAY Y < 0   LTOUT<(M)    ', CLTOUL)        
       CALL PRR ('    OUTBOARD BASE VALUE Y > 0  TBOUT>(EV)   ', CTBOUG)        
       CALL PRR ('    OUTBOARD EXP DECAY Y > 0   LTOUT>(M)    ', CLTOUG)        
C                                                                               
      IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.6) THEN                                      
       CALL PRR ('    INBOARD BASE VALUE         TBIN  (EV)   ', CTBIN)         
       IF (CATIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD LINEAR DECAY       GTIN2 (EV)   ', CGTIN2)        
       ELSEIF (CATIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD LINEAR DECAY       GTIN1 (EV)   ', CGTIN1)        
       ELSE                                                                     
        CALL PRR('    0 < X < ATIN LINEAR DECAY  GTIN1 (EV/M) ', CGTIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ATIN  (M)    ', CATIN)         
        CALL PRR('    ATIN < X < A LINEAR DECAY  GTIN2 (EV/M) ', CGTIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.2) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         TBIN  (EV)   ', CTBIN)         
       IF (CATIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LTIN2 (M)    ', CLTIN2)        
       ELSEIF (CATIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LTIN1 (M)    ', CLTIN1)        
       ELSE                                                                     
        CALL PRR('    0<X<ATIN EXPONENTIAL DECAY LTIN1 (M)    ', CLTIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ATIN  (M)    ', CATIN)         
        CALL PRR('    ATIN<X<A EXPONENTIAL DECAY LTIN2 (M)    ', CLTIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.3) THEN                                                 
       CALL PRC ('    INBOARD FITTED TO SET OF GIVEN VALUES')                   
C                                                                               
      ELSEIF (CIOPTG.EQ.4) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         TBIN  (EV)   ', CTBIN)         
       IF (CATIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LTIN2 (M)    ', CLTIN2)        
       ELSEIF (CATIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD LINEAR DECAY       GTIN1 (EV/M) ', CGTIN1)        
       ELSE                                                                     
        CALL PRR('    0 < X < ATIN LINEAR DECAY  GTIN1 (EV/M) ', CGTIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ATIN  (M)    ', CATIN)         
        CALL PRR('    ATIN<X<A EXPONENTIAL DECAY LTIN2 (M)    ', CLTIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.5) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         TBIN  (EV)   ', CTBIN)         
       IF (CATIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD LINEAR DECAY       GTIN2 (EV/M) ', CGTIN2)        
       ELSEIF (CATIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LTIN1 (M)    ', CLTIN1)        
       ELSE                                                                     
        CALL PRR('    0<X<ATIN EXPONENTIAL DECAY LTIN1 (M)    ', CLTIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ATIN  (M)    ', CATIN)         
        CALL PRR('    ATIN < X < A LINEAR DECAY  GTIN2 (EV/M) ', CGTIN2)        
       ENDIF                                                                    
      ENDIF                                                                     
C                                                                               
      ENDIF                                                                     
C
C     ELECTRON TEMPERATURE GRADIENT FUNCTIONS  
C      
      IF (NTEG.NE.0) THEN
        CALL PRC('  TEMPERATURE GRADIENT FUNCTION: ')
        CALL PRC('    POSITION (Y-LCFS)     T MULTIPLIER ') 
        DO 50 IX = 1,NTEG
          WRITE(COMENT,'(7X,F8.3,10X,F8.3)') TMEG(IX,1)*CL,TMEG(IX,2)
          CALL PRC(COMENT)
50      CONTINUE
      ENDIF
C-----------------------------------------------------------------------        
      IF (CIOPTK.EQ.0) THEN                                                     
       CALL PRR ('  ION TEMPERATURE OUTB Y < 0   TIBOUT<(EV)  ',CTIBOUL)       
       CALL PRR ('    OUTBOARD CONST Y > 0       TIBOUT>(EV)  ',CTIBOUG)        
       CALL PRR ('    INBOARD CONSTANT           TIBIN  (EV)  ', CTIBIN)        
C                                                                               
      ELSE                                                                      
       CALL PRR ('  ION TEMPERATURE OUTB Y < 0   TIBOUT<(EV)  ',CTIBOUL)        
       CALL PRR ('    OUTBOARD EXP DECAY Y < 0   LTIOUT<(M)   ',CLTIOUL)        
       CALL PRR ('    OUTBOARD BASE VALUE Y > 0  TIBOUT>(EV)  ',CTIBOUG)        
       CALL PRR ('    OUTBOARD EXP DECAY Y > 0   LTIOUT>(M)   ',CLTIOUG)        
C                                                                               
      IF (CIOPTK.EQ.1.OR.CIOPTK.EQ.6) THEN                                      
       CALL PRR ('    ION INBOARD BASE VALUE     TIBIN  (EV)  ',CTIBIN)         
       IF (CATIIN.LE.0.0) THEN                                                  
        CALL PRR('    ION INBOARD LINEAR DECAY   GTIIN2 (EV)  ',CGTIIN2)        
       ELSEIF (CATIIN.GE.CA) THEN                                             
        CALL PRR('    ION INBOARD LINEAR DECAY   GTIIN1 (EV)  ',CGTIIN1)        
       ELSE                                                                     
        CALL PRR('    0 < X < ATIIN LINEAR DECAY GTIIN1(EV/M) ',CGTIIN1)        
        CALL PRR('    ION INBOARD CROSSOVER      ATIIN (M)    ',CATIIN)         
        CALL PRR('    ATIIN< X < A LINEAR DECAY  GTIIN2(EV/M) ',CGTIIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTK.EQ.2) THEN                                                 
       CALL PRR ('    ION INBOARD BASE VALUE      TIBIN  (EV) ',CTIBIN)         
       IF (CATIIN.LE.0.0) THEN                                                
        CALL PRR('    ION INBOARD EXP. DECAY      LTIIN2 (M)  ',CLTIIN2)        
       ELSEIF (CATIIN.GE.CA) THEN                                               
        CALL PRR('    ION INBOARD EXP. DECAY      LTIIN1 (M)  ',CLTIIN1)        
       ELSE                                                                     
        CALL PRR('    0<X<ATIIN EXPONENTIAL DECAY LTIIN1 (M)  ',CLTIIN1)        
        CALL PRR('    ION INBOARD CROSSOVER       ATIIN  (M)  ',CATIIN)         
        CALL PRR('    ATIIN<X<A EXPONENTIAL DECAY LTIIN2 (M)  ',CLTIIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTK.EQ.3) THEN                                                 
       CALL PRC ('    ION INBOARD FITTED TO SET OF GIVEN VALUES')           
C                                                                               
      ELSEIF (CIOPTK.EQ.4) THEN                                                 
       CALL PRR ('    ION INBOARD BASE VALUE     TIBIN (EV)   ',CTIBIN)        
       IF (CATIIN.LE.0.0) THEN                                               
        CALL PRR('    ION INBOARD EXP. DECAY     LTIIN2 (M)   ',CLTIIN2)        
       ELSEIF (CATIIN.GE.CA) THEN                                               
        CALL PRR('    ION INBOARD LINEAR DECAY   GTIIN1 (EV/M)',CGTIIN1)        
       ELSE                                                                     
        CALL PRR('    0 < X < ATIIN LINEAR DECAY GTIIN1 (EV/M)',CGTIIN1)        
        CALL PRR('    ION INBOARD CROSSOVER      ATIIN  (M)   ',CATIIN)         
        CALL PRR('    ATIIN<X<A EXP. DECAY       LTIIN2 (M)   ',CLTIIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTK.EQ.5) THEN                                                 
       CALL PRR ('    ION INBOARD BASE VALUE     TIBIN  (EV)  ',CTIBIN)         
       IF (CATIIN.LE.0.0) THEN                                                  
        CALL PRR('    ION INBOARD LINEAR DECAY   GTIIN2 (EV/M)',CGTIIN2)        
       ELSEIF (CATIIN.GE.CA) THEN                                               
        CALL PRR('    ION INBOARD EXP. DECAY     LTIIIN1 (M)  ',CLTIIN1)        
       ELSE                                                                     
        CALL PRR('    0<X<ATIIN EXP. DECAY       LTIIN1 (M)   ',CLTIIN1)        
        CALL PRR('    ION INBOARD CROSSOVER      ATIIN  (M)   ',CATIIN)         
        CALL PRR('    ATIIN < X < A LINEAR DECAY GTIIN2 (EV/M)',CGTIIN2)        
       ENDIF                                                                    
      ENDIF                                                                     
C                                                                               
      ENDIF                                                                     
C
C     ION TEMPERATURE GRADIENT FUNCTIONS  
C      
      IF (NTIG.NE.0) THEN
        CALL PRC('  TEMPERATURE GRADIENT FUNCTION: ')
        CALL PRC('    POSITION (Y-LCFS)     T MULTIPLIER ')
        DO 55 IX = 1,NTEG
          WRITE(COMENT,'(7X,F8.3,10X,F8.3)') TMIG(IX,1)*CL,TMIG(IX,2)
          CALL PRC(COMENT)
55      CONTINUE
      ENDIF
C-----------------------------------------------------------------------        
      IF (CIOPTG.EQ.0) THEN                                                     
       CALL PRR ('  ION DENSITY OUTBOARD Y < 0   NBOUT<(M**-3)', CNBOUL)        
       CALL PRR ('    OUTBOARD CONST Y > 0       NBOUT>(M**-3)', CNBOUG)        
       CALL PRR ('    INBOARD CONSTANT           NBIN  (M**-3)', CNBIN)         
C                                                                               
      ELSE                                                                      
       CALL PRR ('  ION DENSITY OUTBOARD Y < 0   NBOUT<(M**-3)', CNBOUL)        
       CALL PRR ('    OUTBOARD EXP DECAY Y < 0   LNOUT<(M)    ', CLNOUL)        
       CALL PRR ('    OUTBOARD BASE VALUE Y > 0  NBOUT>(M**-3)', CNBOUG)        
       CALL PRR ('    OUTBOARD EXP DECAY Y > 0   LNOUT>(M)    ', CLNOUG)        
C                                                                               
      IF (CIOPTG.EQ.1) THEN                                                     
       CALL PRR ('    INBOARD BASE VALUE         NBIN  (M**-3)', CNBIN)         
       IF (CANIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD LINEAR DECAY       GNIN2 (M**-4)', CGNIN2)        
       ELSEIF (CANIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD LINEAR DECAY       GNIN1 (M**-4)', CGNIN1)        
       ELSE                                                                     
        CALL PRR('    0 < X < ANIN LINEAR DECAY  GNIN1 (M**-4)', CGNIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ANIN  (M)    ', CANIN)         
        CALL PRR('    ANIN < X < A LINEAR DECAY  GNIN2 (M**-4)', CGNIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.2) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         NBIN  (M**-3)', CNBIN)         
       IF (CANIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LNIN2 (M)    ', CLNIN2)        
       ELSEIF (CANIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LNIN1 (M)    ', CLNIN1)        
       ELSE                                                                     
        CALL PRR('    0<X<ANIN EXPONENTIAL DECAY LNIN1 (M)    ', CLNIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ANIN  (M)    ', CANIN)         
        CALL PRR('    ANIN<X<A EXPONENTIAL DECAY LNIN2 (M)    ', CLNIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.3) THEN                                                 
       CALL PRC ('    INBOARD FITTED TO SET OF GIVEN VALUES')                   
C                                                                               
      ELSEIF (CIOPTG.EQ.4) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         NBIN  (M**-3)', CNBIN)         
       IF (CANIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LNIN2 (M)    ', CLNIN2)        
       ELSEIF (CANIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD LINEAR DECAY       GNIN1 (M**-4)', CGNIN1)        
       ELSE                                                                     
        CALL PRR('    0 < X < ANIN LINEAR DECAY  GNIN1 (M**-4)', CGNIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ANIN  (M)    ', CANIN)         
        CALL PRR('    ANIN<X<A EXPONENTIAL DECAY LNIN2 (M)    ', CLNIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.5) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         NBIN  (M**-3)', CNBIN)         
       IF (CANIN.LE.0.0) THEN                                                   
        CALL PRR('    INBOARD LINEAR DECAY       GNIN2 (M**-4)', CGNIN2)        
       ELSEIF (CANIN.GE.CA) THEN                                                
        CALL PRR('    INBOARD EXPONENTIAL DECAY  LNIN1 (M)    ', CLNIN1)        
       ELSE                                                                     
        CALL PRR('    0<X<ANIN EXPONENTIAL DECAY LNIN1 (M)    ', CLNIN1)        
        CALL PRR('    INBOARD CROSSOVER POINT    ANIN  (M)    ', CANIN)         
        CALL PRR('    ANIN < X < A LINEAR DECAY  GNIN2 (M**-4)', CGNIN2)        
       ENDIF                                                                    
C                                                                               
      ELSEIF (CIOPTG.EQ.6) THEN                                                 
       CALL PRR ('    INBOARD BASE VALUE         NBIN  (M**-3)', CNBIN)         
       CALL PRR ('    INBOARD DENSITY AT X = A   NBA   (M**-3)', CNBA)          
       CALL PRR ('    INBOARD PARAMETER          GAMMA        ', CGAMMA)        
      ENDIF                                                                     
C                                                                               
      ENDIF                                                                     
C-----------------------------------------------------------------------        
C     CALL PRR  ('  PINCH PARAMETER              S            ', CVIN)          
C
      CALL PRR2 ('  X DIFF RATE -L/2 < Y < 0 OUT/IN (M*M/S)',                
     >                                                CRDXO(1),CRDXI(1))        
      CALL PRR2 ('  X DIFF RATE  0 < Y < L/2 OUT/IN (M*M/S)',                
     >                                                CRDXO(2),CRDXI(2))        
      CALL PRR2 ('  X DIFF RATE CENTRAL REGIONS O/I (M*M/S)',
     >                                                CRDXO(3),CRDXI(3))
      CALL PRI  ('  DPERPO/DPERPI SPLIT          IQXBRK       ', CBRK)
      IF (CBRK.LE.0) THEN 
         CALL PRR  ('                               X-COORD      ',  
     >               FLOAT(CBRK)/XSCALO)
      ELSEIF (CBRK.GT.0) THEN 
         CALL PRR  ('                               X-COORD      ',  
     >               FLOAT(CBRK)/XSCALI)
      ENDIF
      CALL PRR  ('  DIFFUSION DECAY INBOARD      DD    (M*M/S)', CRDD)          
      CALL PRR  ('  POLOIDAL DIFFUSION RATE      DPOL  (M*M/S)', CDPOL)         
      CALL PRR  ('  POLOIDAL DRIFT VELOCITY      VPOL  (M/S)  ', CVPOL)
      CALL PRR  ('  INBOARD PLASMA FLOW VELOCITY VHYIN (M/S)  ', CVHYIN)        
      CALL PRR  ('  INBOARD ELECTRIC FIELD       EYIN  (V/M)  ', CEYIN)         
      CALL PRR  ('  COLLISION ENHANCEMENT FACTOR ZENH         ', CZENH)         
      CALL PRB                                                                  
      CALL PRC  ('IMPURITY ION CHARACTERISTICS')                                
      CALL PRR  ('  IMPURITY ION MASS            MI           ', CRMI)          
      CALL PRI  ('  IMPURITY ATOMIC NUMBER       ZI           ', CION)          
      CALL PRR  ('  CHARACTERISTIC ENERGY        EBD   (EV)   ', CEBD)          
      IF (CSVYMIN.NE.0.0) 
     >CALL PRR  ('  IMPOSED MINIMUM ION VELOCITY SVYMIN (M/S) ',CSVYMIN)
      IF (CIZSET.GE.1.AND.CIZSET.LE.NIZS)                                       
     >CALL PRI  ('  SET TI=TB FOR IONS REACHING STATE         ', CIZSET)        
C-----------------------------------------------------------------------        
      CALL PRB                                                                  
      CALL PRC  ('SOURCE CHARACTERISTICS')                                      
      CALL PRR  ('  INJECTION TIME           (S)     ', CTIMSC)                 
C                                                                               
      ICXSC = INT (CXSC * XSCALO)                                               
      IF     (CNEUTA.EQ.0) THEN                                                 
       CIZSC = 1                                                                
       CALL PRC ('  INITIAL ION TEMPERATURE  (EV)     PASSED BACK FROM N        
     >EUT')                                                                     
       IF (CNEUTB.EQ.0) THEN                                                    
        CALL PRC ('  SOURCE X POSN FOR NEUT   (M)      DISTRIBUTED ACROS        
     >S LIMITER')                                                               
        CALL PRC ('  SOURCE Y POSN FOR NEUT   (M)      DISTRIBUTED ACROS        
     >S LIMITER')                                                               
       ELSEIF (CNEUTB.EQ.1) THEN                                                
        CALL PRR ('  SOURCE X POSN FOR NEUT   (M)     ', CXSC)                  
        CALL PRR ('  SOURCE Y POSN SYMMETRIC  (M)  +/-',QEDGES(ICXSC,1))        
       ELSEIF (CNEUTB.EQ.2) THEN                                                
        CALL PRR ('  SOURCE X POSN FOR NEUT   (M)     ', CXSC)                  
        CALL PRR ('  SOURCE Y POSN ASYMMETRIC (M)    +',QEDGES(ICXSC,2))        
c slmod
       ELSEIF (CNEUTB.EQ.3.OR.CNEUTB.EQ.9.OR.CNEUTB.EQ.10) THEN
c       ELSEIF (CNEUTB.EQ.3) THEN                                                
c slmod end
C
C       MODIFY LAUNCH OPTION 3 TO SUPPORT A LAUNCH FROM ANY POSITION X,Y
C
C       CXSC = MAX (0.0, CXSC)                                                  
C
        CALL PRR ('  SOURCE X POSN FOR NEUT   (M)     ', CXSC)                  
        CALL PRR ('  SOURCE Y POSN FOR NEUT   (M)     ', CYSC)               
       ELSEIF (CNEUTB.EQ.4) THEN                                                
        CALL PRC ('  SOURCE X POSN FOR NEUT   (M)        WALL')                 
        CALL PRC ('  SOURCE Y POSN FOR NEUT   (M)      HOMOGENEOUS')            
       ELSEIF (CNEUTB.EQ.5) THEN
        CALL PRC ('  SOURCE X POSN FOR NEUT (M)  DISTRIBUTED ACROSS LIMI        
     >TER')                                                               
        CALL PRC ('  SOURCE Y POSN FOR NEUT (M)  DISTRIBUTED ACROSS LIMI       
     >TER')                                                               
        CALL PRR ('  SOURCE P POSN FOR NEUT (M)  DISTRIBUTED ACROSS +/-'     
     >,CPSC)                                                            
       ELSEIF (CNEUTB.EQ.6) THEN
        CALL PRR ('  SOURCE X POSN FOR NEUT (M)                    '
     >,CXSC)       
        CALL PRR ('  SOURCE Y POSN FOR NEUT (M)  DISTRIBUTED ON +/-'      
     >,CYSC)                                                               
       ELSEIF (CNEUTB.EQ.7) THEN
        CALL PRC ('  SOURCE X POSN FOR NEUT BASED ON INPUT PROBABILITY D
     >ISTRIBUTION')
        CALL PRR2('  SOURCE Y POSN FOR NEUT DISTRIBUTED ON: '      
     >,Y0S,Y0L)                                                               
       ELSEIF (CNEUTB.EQ.8) THEN
        CALL PRR2('  SOURCE X POSN FOR NEUT DISTRIBUTED ON: '
     >,X0S,X0L)       
        CALL PRC ('  SOURCE Y POSN FOR NEUT BASED ON INPUT PROBABILITY D     
     >ISTRIBUTION')                                                      
       ENDIF                                                                    
c
c     Initial Y coordinate over ride option
c     
       if (init_y_coord.ne.0.0) then 
          call prb
          call prc('  *D98 - SPECIFIED INITIAL Y COORDINATE')
          call prc('       - THIS TAKES'//
     >                      ' PRECEDENCE OVER ANY OTHER OPTIONS')
          call prr('  SOURCE Y POSN IS SEPARATELY SPECIFIED =',
     >      init_y_coord)  
          call prb
       endif


c slmod begin - N2 break
      IF (optdp.EQ.1) THEN
        CALL PRI ('DIVIMP ION SOURCE OPTION          ', optdp)
        CALL PRC ('    Record where ions enter the overlap region ')
        CALL PRC ('    between a LIM and DIVIMP grid to be used as')
        CALL PRC ('    ion source distribution for DIVIMP.  The')
        CALL PRC ('    extent of the DIVIMP grid is currently hard')
        CALL PRC ('    coded in the LIM3 module.')
      ENDIF

      IF (N2OPT.EQ.1) THEN
	WRITE(7,*) ' '
        CALL PRI ('MOLECULAR NITROGEN OPTION         ', N2OPT)
        CALL PRC ('    Nitrogen neutrals are launched as molecules')
        CALL PRC ('    and follow several break-up pathways.      ')
        CALL PRR ('  Fraction of N2 events to give fast neutrals: ',
     +            N2FRAC)
        CALL PRI ('  Fast neutral energy option: ',neopt)
        IF (neopt.EQ.1) THEN
          CALL PRR ('    Fast neutral energy: ',neng)
        ENDIF
        CALL PRI ('  Ion energy option: ',nieopt)
        IF (nieopt.EQ.1) THEN
          CALL PRR ('    Ion energy: ',nieng)
        ENDIF
	WRITE(7,*) ' '
      ENDIF
c slmod end
      ELSEIF (CNEUTA.EQ.1) THEN                                                 
       IF (CIOPTE.EQ.4) THEN
         CALL PRR ('  INITIAL ION ENERGY       (EV)    ', CENGSC)
         call prr ('  INITIAL ION VELOCITY     (M/S)   ',
     >                      1.38E4* SQRT(CENGSC/CRMI))
       ELSE
         CALL PRR ('  INITIAL ION TEMPERATURE  (EV)    ', CTEMSC)           
         call prr ('  INITIAL ION VELOCITY     (M/S)   ',
     >                      1.56E4* SQRT(CTEMSC/CRMI))
       ENDIF
       IF (CIOPTE.EQ.9) THEN 
        CALL PRR2 ('  SOURCE X POSITION UNIFORM ON (M)  ',X0S,X0L)                   
       ELSE
        CALL PRR ('  SOURCE X POSITION        (M)     ', CXSC)                   
       ENDIF
       IF (CIOPTE.EQ.2.OR.CIOPTE.EQ.3.OR.CIOPTE.EQ.7.OR.CIOPTE.EQ.8)
     >    THEN                      
        CALL PRC ('  SOURCE Y POSITION        (M)      HOMOGENEOUS')            
       ELSEIF (CIOPTE.EQ.4) THEN                                              
        CALL PRR ('  SOURCE Y POSN RANDOM     (M)  +/-', CYSC)                  
       ELSE                                                                     
        CALL PRR ('  SOURCE Y POSN SYMMETRIC  (M)  +/-', CYSC)                  
       ENDIF                                                                    
      ELSEIF (CNEUTA.EQ.2) THEN                                                 
       CALL PRR ('  INITIAL ION TEMPERATURE  (EV)     ', CTEMSC)                
       WRITE (COMENT,'(''  SOURCE X POSITION  (M)     UNIFORMLY'',              
     >  '' FROM ('',F7.4,'','',F7.4,'')'')') X0S,X0L                            
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(''  SOURCE Y POSITION  (M)     UNIFORMLY'',              
     >  '' FROM ('',F7.4,'','',F7.4,'')'')') Y0S,Y0L                            
       CALL PRC (COMENT)                                                        
      ENDIF                                                                     
      IF (CNEUTB.NE.5) THEN
        IF (CNEUTC.EQ.11) THEN                                                
          CALL PRR ('  SOURCE P POSITION BETWEEN +/-    ', CPSC)              
        ELSE                                                                 
          CALL PRR ('  SOURCE P POSITION        (M)     ', CPSC)              
        ENDIF                                                                 
      ENDIF 
      CALL PRI ('  INITIAL IONIZATION STATE         ', CIZSC)                   
C-----------------------------------------------------------------------        
      CALL PRB                                                                  
      CALL PRC ('OPTIONS')                                                      
      IF     (CIOPTA.EQ.0) THEN                                                 
       CALL PRC ('  IONISATION OPT   0 : RATES FROM S(Z,TE) DATA')              
      ELSEIF (CIOPTA.EQ.1) THEN                                                 
       CALL PRC ('  IONISATION OPT   1 : RATES FROM S(Z,TE) DATA')              
      ELSEIF (CIOPTA.EQ.2) THEN                                                 
       CALL PRC ('  IONISATION OPT   2 : RATES TAKEN AS MAX (S(Z,TE))')         
      ELSEIF (CIOPTA.EQ.3) THEN                                                 
       CALL PRC ('  IONISATION OPT   3 : RATES FROM ABELS VAN MAANEN WIT        
     >H E-I REC.')                                                              
      ELSEIF (CIOPTA.EQ.4) THEN                                                 
       CALL PRC ('  IONISATION OPT   4 : RATES FROM ABELS VAN MAANEN')          
      ELSEIF (CIOPTA.EQ.5) THEN                                                 
       CALL PRC ('  IONISATION OPT   5 : RATES FROM ABELS VAN MAANEN WIT        
     >H E-I REC.')                                                              
      ELSEIF (CIOPTA.EQ.6) THEN                                                 
       CALL PRC ('  IONISATION OPT   6 : RATES FROM ABELS VAN MAANEN')          
      ENDIF                                                                     
C                                                                               
      IF (NIZS.LT.CION) THEN                                                    
       IF (CIOPTA.EQ.0.OR.CIOPTA.EQ.3.OR.CIOPTA.EQ.4) THEN                      
        WRITE (COMENT,'(23X,''IONS NOT FOLLOWED AFTER STATE'',I3)') NIZS        
        CALL PRC (COMENT)                                                       
       ELSE                                                                     
        WRITE (COMENT,'(23X,''IONISATION DISABLED AFTER STATE'',I3)')           
     >   NIZS                                                                   
        CALL PRC (COMENT)                                                       
       ENDIF                                                                    
      ENDIF                                                                     
      IF (CNEUTA.EQ.0)                                                          
     > CALL PRR ('                       NEUTRAL IONISATION RATE FACTOR'        
     >, CIRF)                                                                   
C-----------------------------------------------------------------------        
      WRITE (COMENT,'(23X,A,1P,G11.4,A)') 'OUTSIDE OF X =',CPLSMA,              
     >  ' ONLY. ELSEWHERE'                                                      
C-----------------------------------------------------------------------        
      IF     (CIOPTB.EQ.0) THEN                                                 
       CALL PRC ('  COLLISION OPTION 0 : TAU PARA = MI.SQRT(TB/MB).TI/')        
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z        
     >ENH.LAM)')                                                                
      ELSEIF (CIOPTB.EQ.1) THEN                                                 
       CALL PRC ('  COLLISION OPTION 1 : TAU PARA = INFINITY, NO DIFFUSI        
     >ON')                                                                      
      ELSEIF (CIOPTB.EQ.2) THEN                                                 
       CALL PRC ('  COLLISION OPTION 2 : TAU PARA = SQRT(MI.TI).TI/')           
       CALL PRC ('                               (6.8E4.NB.ZB.ZEFF.ZI.ZI        
     >.LAM)')                                                                   
       CALL PRI ('                       WHERE ZEFF(SELF) =',CIZEFF)            
      ELSEIF (CIOPTB.EQ.3) THEN                                                
       CALL PRC ('  COLLISION OPTION 3 : TAU PARA = MI.SQRT(TB/MB).TI/')        
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z        
     >ENH.LAM)')                                                                
       CALL PRC ('                       TIME BETWEEN Y DIFF STEPS = TAU        
     > PARA')                                                                   
       IF (CPLSMA.LT.CA) THEN                                                   
       CALL PRC (COMENT)                                                        
       CALL PRC ('                       TIME BETWEEN Y DIFF STEPS = DEL        
     >TA T')                                                                    
       ENDIF
      ELSEIF (CIOPTB.EQ.4) THEN                                                
       CALL PRC ('  COLLISION OPTION 4 : TAU PARA = 0.5.MI.SQRT(TB/MB).T
     >I/')        
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z        
     >ENH.LAM)')                                                                
       CALL PRC ('                       TIME BETWEEN Y DIFF STEPS = TAU        
     > PARA')                                                                   
c slmod
      ELSEIF (CIOPTB.EQ.13) THEN
       CALL PRC ('  COLLISION OPTION 13: PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = Rg * SQRT(2kTI/MI)')
       call prc ('                           * sqrt( dt / (tau para))')
       call prc ('                       At every time step       ')
       call prc ('                       Where Rg = SQRT(-2*ln(X1))* cos
     >(2PI*X2) ')
       call prc ('                       X1,X2 are uniform on [0,1] ')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                               / (1.0+MB/MI) ')
      ENDIF                                                                     
C                                                                               
c slmod end
      IF ((CIOPTB.EQ.1.OR.CIOPTB.EQ.2) .AND. CPLSMA.LT.CA) THEN                 
       CALL PRC (COMENT)                                                        
       CALL PRC ('                       TAU PARA = MI.SQRT(TB/MB).TI/')        
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z        
     >ENH.LAM)')                                                                
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF     (CIOPTC.EQ.0) THEN                                                 
       CALL PRC ('  FRICTION OPTION  0 : TAU STOP = MI.TB.SQRT(TB/MB)/')        
       CALL PRC ('                       (6.8E4.(1+MB/MI).NB.ZB.ZB.ZI.ZI        
     >.ZENH.LAM)')                                                              
      ELSEIF (CIOPTC.EQ.1) THEN                                                 
       CALL PRC ('  FRICTION OPTION  1 : TAU STOP = INFINITY')                  
      ELSEIF (CIOPTC.EQ.2) THEN                                                 
       CALL PRC ('  FRICTION OPTION  2 : TAU STOP = TAU PARALLEL')              
      ELSEIF (CIOPTC.EQ.3) THEN                                                 
       CALL PRC ('  FRICTION OPTION  3 : TAU STOP = 0')              
      ENDIF                                                                     
C                                                                               
      IF (CIOPTC.NE.0.AND.CPLSMA.LT.CA) THEN                                    
       CALL PRC (COMENT)                                                        
       CALL PRC ('                       TAU STOP = MI.TB.SQRT(TB/MB)/')        
       CALL PRC ('                       (6.8E4.(1+MB/MI).NB.ZB.ZB.ZI.ZI        
     >.ZENH.LAM)')                                                              
       ENDIF                                                                    
C-----------------------------------------------------------------------        
      IF     (CIOPTD.EQ.0) THEN                                                 
       CALL PRC ('  HEATING OPTION   0 : TAU HEAT = MI.TB.SQRT(TB/MB)/')        
       CALL PRC ('                               (1.4E5.NB.ZB.ZB.ZI.ZI.Z        
     >ENH.LAM)')                                                                
      ELSEIF (CIOPTD.EQ.1) THEN                                                 
       CALL PRC ('  HEATING OPTION   1 : TAU HEAT = INFINITY')                  
      ELSEIF (CIOPTD.EQ.2) THEN                                                 
       CALL PRC ('  HEATING OPTION   2 : TAU HEAT = ZERO')                      
      ELSEIF (CIOPTD.EQ.3) THEN                                                 
       CALL PRC ('  HEATING OPTION   3 : TAU HEAT =(MI.TB+MB.TI)**3/2 /'        
     >)                                                                         
       CALL PRC ('                       (1.4E5SQRT(MI.MB)NB.ZB.ZB.ZI.ZI        
     >.ZENH.LAM)')                                                              
      ELSEIF (CIOPTD.EQ.4) THEN                                                 
       CALL PRC ('  HEATING OPTION   4 : RETURN NO VALUES')                  
      ENDIF                                                                     
C                                                                               
      IF (CIOPTD.NE.0.AND.CPLSMA.LT.CA) THEN                                    
       CALL PRC (COMENT)                                                        
       CALL PRC ('                       TAU HEAT = MI.TB.SQRT(TB/MB)/')        
       CALL PRC ('                               (1.4E5.NB.ZB.ZB.ZI.ZI.Z        
     >ENH.LAM)')                                                                
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF     (CNEUTA.EQ.0) THEN                                                 
       CALL PRC ('  INJECTION OPT    * : FROM NEUT')                            
      ELSEIF (CIOPTE.EQ.0) THEN                                                 
       CALL PRC ('  INJECTION OPT    0 : INJECT IONS AT (X0,Y0,P0) WITH         
     >V0 = 0')                                                                  
      ELSEIF (CIOPTE.EQ.1) THEN                                                 
       CALL PRC ('  INJECTION OPT    1 : INJECT IONS AT (X0,Y0,P0) ALONG
     > +Y DIRECTION WITH GIVEN V0')                                                                
      ELSEIF (CIOPTE.EQ.2) THEN                                                 
       CALL PRC ('  INJECTION OPT    2 : HOMOGENEOUS INJECTION AT (X0,P0        
     >), V0 = 0')                                                               
      ELSEIF (CIOPTE.EQ.3) THEN                                                 
       CALL PRC ('  INJECTION OPT    3 : HOMOGENEOUS INJECTION AT (X0,P0        
     >) GIVEN V0')                                                              
      ELSEIF (CIOPTE.EQ.4) THEN                                                 
       CALL PRC ('  INJECTION OPT    4 : HOMOGENEOUS INJECTION ON (X0,+/
     >-Y0,P0)  ')          
       CALL PRC ('                       WITH GIVEN V01 (P1) OR V02 (P2)
     >         ')
      ELSEIF (CIOPTE.EQ.5) THEN                                                 
       CALL PRC ('  INJECTION OPT    5 : INJECT IONS AT (X0,Y0,P0) WITH         
     >V0 = 0')                                                                  
      ELSEIF (CIOPTE.EQ.6) THEN                                                 
       CALL PRC ('  INJECTION OPT    6 : INJECT IONS AT (X0,Y0,P0) WITH         
     >GIVEN V0')                                                                
      ELSEIF (CIOPTE.EQ.7) THEN                                                 
       CALL PRC ('  INJECTION OPT    7 : HOMOGENEOUS INJECTION AT (X0,P0        
     >), V0 = 0')                                                               
      ELSEIF (CIOPTE.EQ.8) THEN                                                 
       CALL PRC ('  INJECTION OPT    8 : HOMOGENEOUS INJECTION AT (X0,P0        
     >) GIVEN V0')                                                              
      ELSEIF (CIOPTE.EQ.9) THEN 
       CALL PRC ('  INJECTION OPT    9 : INJECT IONS ON ([X1,X2],Y0,P0)  
     >GIVEN V0')          
c slmod
      ELSEIF (CIOPTE.EQ.10) THEN
       CALL PRC ('  INJECTION OPT   10 : GAUSSIAN LAUNCH CENTERED AT (X0
     >,Y0,P0) ')
      ELSEIF (CIOPTE.EQ.12) THEN
       CALL PRC ('  INJECTION OPT   12 : INJECT IONS AT (X0,Y0,P0) ALONG
     > -Y DIRECTION WITH GIVEN V0')          
c slmod end
      ELSEIF (CIOPTE.EQ.13) THEN                                                 
       CALL PRC ('  INJECTION OPT   13 : INJECT IONS AT (X0,Y0,P0) WITH         
     >V0 = +/-VY0')                                                                  
      ENDIF                                                                     
C                                                                               
      IF (CIOPTE.EQ.5.OR.CIOPTE.EQ.6.OR.CIOPTE.EQ.7.OR.CIOPTE.EQ.8) 
     > CALL PRC ('                       X0 AND Y0 SUBJECT TO SMALL VARI        
     >ATIONS')                                                                  
C-----------------------------------------------------------------------        
      IF     (CIOPTF.EQ.0) THEN                                                 
       CALL PRC ('  SOL OPTION       0 : SOL0,  VHYOUT=0,  EYOUT=0')            
      ELSEIF (CIOPTF.EQ.1) THEN                                                 
       CALL PRC ('  SOL OPTION       1 : SOL1')                                 
      ELSEIF (CIOPTF.EQ.2) THEN                                                 
       CALL PRC ('  SOL OPTION       2 : SOL2')                                 
      ELSEIF (CIOPTF.EQ.3) THEN                                                 
       CALL PRC ('  SOL OPTION       3 : SOL3')                                 
      ELSEIF (CIOPTF.EQ.4) THEN                                                 
       CALL PRC ('  SOL OPTION       4 : SOL4')                                 
      ELSEIF (CIOPTF.EQ.5) THEN                                                 
       WRITE (COMENT,'(''  SOL OPTION       5 : SOL5,  VHYOUT='',               
     >   1P,G11.4,'',  EYOUT='',G11.4)') CVHOUT,CEYOUT                          
       CALL PRC (COMENT)                                                        
      ELSEIF (CIOPTF.EQ.6) THEN                                                 
       WRITE (COMENT,'(''  SOL OPTION       6 : SOL6, STAGNATION POINT Y        
     >STAG='',G11.4)') CYSTAG                                                   
       CALL PRC (COMENT)                                                        
      ELSEIF (CIOPTF.EQ.7) THEN                                                 
       WRITE (COMENT,'(''  SOL OPTION       7 : SOL7, STAGNATION POINT Y        
     >STAG='',G11.4)') CYSTAG                                                   
       CALL PRC (COMENT)                                                        
      ELSEIF (CIOPTF.EQ.8) THEN                                                 
       CALL PRR ('  SOL OPTION       8 : SOL8, SOL ENHANCEMENT FACTOR=',        
     >   CSOLEF)                                                                
      ELSEIF (CIOPTF.EQ.9) THEN                                                 
       CALL PRC ('  SOL OPTION       9 : SOL9')                                 
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF     (CIOPTG.EQ.0) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 0 : CONSTANT OUTBOARD')                    
       CALL PRC ('                       CONSTANT INBOARD')                     
      ELSEIF (CIOPTG.EQ.1) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 1 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       LINEAR / LINEAR INBOARD')              
      ELSEIF (CIOPTG.EQ.2) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 2 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       EXPONENTIAL / EXPONENTIAL INBOA        
     >RD')                                                                      
      ELSEIF (CIOPTG.EQ.3) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 3 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       FITTED TO GIVEN VALUES INBOARD'        
     >)                                                                         
      ELSEIF (CIOPTG.EQ.4) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 4 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       LINEAR / EXPONENTIAL INBOARD')         
      ELSEIF (CIOPTG.EQ.5) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 5 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       EXPONENTIAL / LINEAR INBOARD')         
      ELSEIF (CIOPTG.EQ.6) THEN                                                 
       CALL PRC ('  PLASMA DECAY OPT 6 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       LINEAR TEMPERATURE INBOARD')           
       CALL PRC ('                       "STANDARD JET" DENSITY INBOARD'        
     >)                                                                         
      ENDIF                                                                    
      IF (NTEG.GT.0) 
     > CALL PRC ('                       TEMPERATURE GRADIENT OUTBOARD')     
C-----------------------------------------------------------------------        
      IF  (CTICHG) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT-1 : TI=TE ')
       CALL PRC ('                       ONLY TE DATA USED')         
       CALL PRC ('                       ALL TI INPUT IGNORED')        
      ELSEIF (CIOPTK.EQ.0) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 0 : CONSTANT OUTBOARD')                    
       CALL PRC ('                       CONSTANT INBOARD')                     
      ELSEIF (CIOPTK.EQ.1) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 1 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       LINEAR / LINEAR INBOARD')              
      ELSEIF (CIOPTK.EQ.2) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 2 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       EXPONENTIAL / EXPONENTIAL INBOA        
     >RD')                                                                      
      ELSEIF (CIOPTK.EQ.3) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 3 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       FITTED TO GIVEN VALUES INBOARD'        
     >)                                                                         
      ELSEIF (CIOPTK.EQ.4) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 4 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       LINEAR / EXPONENTIAL INBOARD')         
      ELSEIF (CIOPTK.EQ.5) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 5 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       EXPONENTIAL / LINEAR INBOARD')         
      ELSEIF (CIOPTK.EQ.6) THEN                                                 
       CALL PRC ('  ION TEMP.    OPT 6 : EXPONENTIAL OUTBOARD')                 
       CALL PRC ('                       LINEAR TEMPERATURE INBOARD')           
      ENDIF                                                                     
      IF (NTIG.GT.0) 
     > CALL PRC ('                       TEMPERATURE GRADIENT OUTBOARD')     
C------------------------------------------------------------------------
C     ELECTRON TEMPERATURE GRADIENT COEFFICIENT OPTION
C
      IF (CIOPTL.EQ.0) THEN 
       CALL PRC ('  Te GRAD COEF OPT 0 : ALPHA = 0, ALL Z')
      ELSEIF (CIOPTL.EQ.1) THEN 
       CALL PRC ('  Te GRAD COEF OPT 1 : ALPHA VALUES = .71 Zi^2')    
       WRITE (COMENT,'(22X,8F6.2)')
     >        (CALPHE(IZ),IZ=1,NIZS)
       CALL PRC(COMENT)
      ENDIF          
C-----------------------------------------------------------------------
C     ION TEMPERATURE GRADIENT COEFFICIENT OPTION
C
      IF (CIOPTM.EQ.0) THEN 
       CALL PRC ('  Ti GRAD COEF OPT 0 : BETA = 0 , ALL Z')
      ELSEIF (CIOPTM.EQ.1) THEN 
       CALL PRC ('  Ti GRAD COEF OPT 1 : BETA VALUES')    
       CALL PRC ('                      '//
     >           ' -3(1-u-5*Zi^2*sqrt(2u)u(1.1u-.35))')
       CALL PRC ('                       /(2.6-2u+5.4u^2)'//
     >           '   u=Mi/(Mi+Mb)')
       WRITE (COMENT,'(22X,8F6.2)') 
     >       (CBETAI(IZ),IZ=1,NIZS)
       CALL PRC(COMENT)
      ENDIF          
C-----------------------------------------------------------------------        
      IF     (CIOPTH.EQ.0) THEN                                                 
       CALL PRC ('  LIMITER EDGE OPT 0 : SLAB : EDGE 0 EVERYWHERE')             
      ELSEIF (CIOPTH.EQ.1) THEN                                                 
       WRITE (COMENT,'(A,F7.2,A)') '  LIMITER EDGE OPT 1 : WEDGE, ANGLE         
     >',WEDGAN(1),' DEGS Y<0'                                                   
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(A,F7.2,A)') '                                            
     >',WEDGAN(2),' DEGS Y>0'                                                   
       CALL PRC (COMENT)                                                        
      ELSEIF (CIOPTH.EQ.2) THEN                                                 
       CALL PRC ('  LIMITER EDGE OPT 2 : CIRCLE CENTRE (-0.065,0), RADIU        
     >S 0.065')                                                                 
       CALL PRC ('                       LAUNCH CUTOFF AT X=-0.1')              
      ELSEIF (CIOPTH.EQ.3) THEN                                                 
       CALL PRC ('  LIMITER EDGE OPT 3 : BLUNT NOSE, DEFINED BY THE POIN        
     >TS')                                                                      
       WRITE (COMENT,'(''                       ('',                            
     >  F7.4,'','',F7.4,'') & ('',F7.4,'','',F7.4,'') Y<0'')')                  
     >  XL1(1),-YL1(1),XL2(1),-YL2(1)                                           
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(''                       ('',                            
     >  F7.4,'','',F7.4,'') & ('',F7.4,'','',F7.4,'') Y>0'')')                  
     >  XL1(2),YL1(2),XL2(2),YL2(2)                                             
       CALL PRC (COMENT)                                                        
      ELSEIF (CIOPTH.EQ.4.OR.CIOPTH.EQ.5) THEN                                  
       WRITE (COMENT,'(''  LIMITER EDGE OPT'',I2,'' : BELT, (TC,SC)='',         
     >  2F10.6)') CIOPTH,TC,SC                                                  
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(''                       (TO,SO)='',2F10.6,              
     >  '' GC='',F8.4)') TO,SO,GC*RADDEG                                        
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(''                       CAMERA (TV,SV)='',              
     >   2F10.6)') TV,SV                                                        
       CALL PRC (COMENT)                                                        
       IF (CIOPTH.EQ.5) CALL PRC ('                       SYMMETRICAL VE        
     >RSION OF Y>0 SHAPE')                                                      
      ELSEIF (CIOPTH.EQ.6) THEN 
       CALL PRC ('  LIMITER EDGE OPT 6 : GENERAL STEP ,'//
     >   ' DEFINED BY THE POINTS')
       WRITE (COMENT,'(''                       FOR Y<0  ('',
     >  F7.4,'','',F7.4,'')'')') XST1(1),-YST1(1)
       CALL PRC(COMENT)
       WRITE(COMENT,'(''                       -------  ('',    
     >  F7.4,'','',F7.4,'')'')') XST2(1),-YST2(1)
       CALL PRC(COMENT)
       WRITE(COMENT,'(''                                ('',    
     >  F7.4,'','',F7.4,'')'')') XST3(1),-YST3(1)
       CALL PRC(COMENT)
       WRITE (COMENT,'(''                       FOR Y>0  ('',
     >  F7.4,'','',F7.4,'')'')') XST1(2),YST1(2)
       CALL PRC(COMENT)
       WRITE(COMENT,'(''                       -------  ('',    
     >  F7.4,'','',F7.4,'')'')') XST2(2),YST2(2)
       CALL PRC(COMENT)
       WRITE(COMENT,'(''                                ('',    
     >  F7.4,'','',F7.4,'')'')') XST3(2),YST3(2)
       CALL PRC(COMENT)
      ELSEIF (CIOPTH.EQ.7) THEN  
       CALL PRC ('  LIMITER EDGE OPT 7 : TFTR BELT TYPE, DEFINED'
     >//' BY RADII')
       WRITE(COMENT,'(23X,''LIMITER RADIUS OF CURVATRURE:'',
     >  F7.4)') RLEDGE7
       CALL PRC(COMENT)
       WRITE(COMENT,'(23X,''PLASMA RADIUS:               '',    
     >  F7.4)') CA
       CALL PRC(COMENT)
       WRITE(COMENT,'(23X,''WALL RADIUS:                 '',    
     >  F7.4)') CA-CAW
       CALL PRC(COMENT)
      ELSEIF (CIOPTH.EQ.8) THEN 
       CALL PRC ('  LIMITER EDGE OPT 8 : ROUND TIP,'//
     >   ' DEFINED BY THE RADIUS OF CURVATURE')
       WRITE (COMENT,'(''                       RLC = '',
     >   F7.4,'' M '')') RLC
       CALL PRC(COMENT)
      ELSEIF (CIOPTH.EQ.9) THEN 
       CALL PRC('  LIMITER EDGE OPT 9 : ALT-II PUMPED LIMITER SHAPE')
       CALL PRC('                       DEFINED BY 18 SPECIFIED POINTS')
       CALL PRC('                       USING LINEAR INTERPOLATION TO')
       CALL PRC('                       CALCULATE INTERMEDIATE POINTS')
       CALL PRC('                       AND NORMALS. SYMMETRIC.') 
       CALL PRC('                       EXTENDS STRAIGHT TO WALL')
       CALL PRC('                       BEYOND LAST POINT.')
      elseif (ciopth.eq.10) then
       call prc('  LIMITER EDGE OPT 10: ITER LIMITER SHAPE')
       call prc('                       DEFINED BY:')
       call prc('                       X=-(0.04*Y**2 + 0.09*Y**4)')
      elseif (ciopth.eq.11) then
       call prc('  LIMITER EDGE OPT 11: ITER HALF LIMITER SHAPE')
       call prc('                       DEFINED BY:')
       call prc('                       '//
     >               'X=- lam * ln(1 -/+ C*(Y - Y_re-entrant)/lam ))')
       call prc('                       REMAPPED SO THAT LIMITER TIP')
       call prc('                       IS AT Y=0 INSTEAD OF Y=Y_re')
       call prc('                       ONLY A SINGLE TANGENCY POINT')
       call prc('                       IS SUPPORTED')
       call prc('                       INPUT PARAMETERS:')
       call prr('                       Toroidal setback  = ',
     >                                  rtor_setback)
       call prr('                       Slot     setback  = ',
     >                                  rslot_setback)
       call prr('                       BM Toroidal half width    = ',
     >                                  bm_tor_wid)
       call prr('                       Slot Toroidal half width  = ',
     >                                  slot_tor_wid)
       call prr('                       Lambda_design  = ',
     >                                  lambda_design)
       call prc('                       DERIVED PARAMETERS:')
       call prr('                       T_Re-entrant = Y_re = ', y_re)
       call prr('                       C                   = ', c_lim)
      ENDIF                                                                     
C                                                                               
      call prb
c
      if (lim_wall_opt.eq.0) then 
         call prr('                       LIM WALL OPTION 0: WALL'//
     >            ' IS AT A CONSTANT DISTANCE FROM THE LCFS: AW =',caw)
      elseif (lim_wall_opt.eq.1) then 
         call prc('                       LIM WALL OPTION 1: WALL'//
     >            ' DISTANCE TO LCFS VARIES')
         call prr('                                          WALL'//
     >            ' STARTS AT A DISTANCE OF AW = ',caw)
         call prr('                                          WALL'//
     >            ' DISTANCE THEN CHANGES LINEARLY'//
     >            ' STARTING AT Y = ',ywall_start)
         call prr('                                          WALL'//
     >            ' DISTANCE TO LCFS AT THE MID POINT BETWEEN '//
     >            ' LIMITERS = ',caw_min)
      endif
c
      call prb
c
      if (yreflection.eq.0) then 
         call prc('                        Y REFLECTION OPT 0:'//
     >            ' REFLECTION IN Y-AXIS MIRRORS IS OFF')
      elseif (yreflection.eq.1) then 
         call prc('                        Y REFLECTION OPT 1:'//
     >            ' Y-AXIS REFLECTION IS ON AT TWO MIRROR'//
     >            ' LOCATIONS: ONE EACH FOR Y>0 AND Y<0')
         call prr('                                           '//
     >            ' Y < 0 MIRROR LOCATION = ',cmir_refl_lower)
         call prr('                                           '//
     >            ' Y > 0 MIRROR LOCATION = ',cmir_refl_upper)
      endif


      IF (CORECT.EQ.1)                                                          
     > WRITE (7,'(24X,''CURVATURE CORRECTED, RP='',F10.6)') RP                  
C-----------------------------------------------------------------------        
      IF (CIOPTJ.EQ.1) THEN
         CALL PRR ('                       LIMITER'
     >     //' POLOIDAL EXTENT +/-',CPCO)
c
c        Add print out for shear short circuit option
c
         if (shear_short_circuit_opt.eq.1) then 
          call prc('                       SHEAR SHORT CIRCUIT OPT 1:'//
     >    ' P RESET TO +/-POLOIDAL EXTENT')
          call prc('                                                 '//
     >    ' IF |Y| CROSSES L_CON/2')
         endif
      ENDIF


C-----------------------------------------------------------------------
      IF     (CIOPTI.EQ.0) THEN                                                 
       CALL PRC ('  CX RECOMB OPTION 0 : NO CHARGE EXCHANGE RECOMBINATIO        
     >N')                                                                       
      ELSEIF (CIOPTI.EQ.1) THEN                                                 
       CALL PRC ('  CX RECOMB OPTION 1 : NH=NHC+NHO.EXP(-X/LAMHX)EXP(-Y/        
     >LAMHY) X>0')                                                              
       CALL PRC ('                          NHO.EXP(-Y/LAMHY)                   
     >       X<0')                                                              
       CALL PRC ('                       AND VCX=SQRT(2TB/MB)')                 
       WRITE (COMENT,'(23X,''WHERE  NHC='',1P,G11.4,'',   NHO='',               
     >   G11.4)') CNHC,CNHO                                                     
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(23X,''AND  LAMHX='',1P,G11.4,'', LAMHY='',               
     >   G11.4)') CLAMHX,CLAMHY                                                 
       CALL PRC (COMENT)                                                        
      ELSEIF (CIOPTI.EQ.2) THEN                                                 
       CALL PRC ('  CX RECOMB OPTION 2 : NH=NHC+NHO.EXP(-X/LAMHX)EXP(-Y/        
     >LAMHY) X>0')                                                              
       CALL PRC ('                          NHO.EXP(-Y/LAMHY)                   
     >       X<0')                                                              
       CALL PRR ('                       WITH CONSTANT VCX =', CVCX)            
       WRITE (COMENT,'(23X,''WHERE  NHC='',1P,G11.4,'',   NHO='',               
     >   G11.4)') CNHC,CNHO                                                     
       CALL PRC (COMENT)                                                        
       WRITE (COMENT,'(23X,''AND  LAMHX='',1P,G11.4,'', LAMHY='',               
     >   G11.4)') CLAMHX,CLAMHY                                                 
       CALL PRC (COMENT)                                                        
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (CIOPTB.EQ.1) THEN                                                     
       CALL PRC ('  FIRST DIFFUSION  * : IGNORED')                              
      ELSE                                                                      
      IF     (CDIFOP.EQ.0) THEN                                                 
       CALL PRC ('  FIRST DIFFUSION  0 : DIFFUSION STARTS IMMEDIATELY')         
      ELSEIF (CDIFOP.EQ.1) THEN                                                 
       CALL PRC ('  FIRST DIFFUSION  1 : AFTER T=-TAUPARA.LOG$/2, BASED         
     >ON TI(0)')                                                                
      ELSEIF (CDIFOP.EQ.2) THEN                                                 
       CALL PRC ('  FIRST DIFFUSION  2 : AFTER T=TAUPARA, ITSELF CHANGIN        
     >G WITH T')                                                                
      ENDIF                                                                     
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF     (CIOPTN.EQ.0) THEN                                                 
       CALL PRC ('  DIFFUSION TYPE   0 : DIFFUSION AT EACH TIME STEP')         
      ELSEIF (CIOPTN.EQ.1) THEN                                                
       CALL PRC ('  DIFFUSION TYPE   1 : DIFFUSION OCCURS IN DISCRETE ST
     >EPS')
       CALL PRR ('                       OF SIZE : ',CDPSTP)
       CALL PRC ('                       ON RANDOM TIME STEPS')
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF     (CDPERP.EQ.0) THEN                                                 
       CALL PRC ('  DPERP OPTION     0 : VALUES AS GIVEN - LINEAR INBOAR
     >D') 
       CALL PRR ('                       WITH SLOPE : ', CRDD)       
      ELSEIF (CDPERP.EQ.1) THEN                                                
       CALL PRC ('  DPERP OPTION     1 : VALUES AS GIVEN - DEVELOPMENT I
     >NBOARD')
       CALL PRC ('                       DESCRIBED BY')
       CALL PRC ('                       DP(X) = DP0*(1+a[(CA-x)/CA]^b)'
     >)
       CALL PRR ('                       a = ', DPALPH)
       CALL PRR ('                       b = ', DPBETA)
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF     (CVPOPT.EQ.0) THEN                                                 
       CALL PRC ('  PINCH OPTION     0 : INWARD PINCH VELOCITY DESCRIBED
     >BY:')
       CALL PRC ('                       VIN = 2*S*DP*(CA-X)/CA^2')
       CALL PRR ('                       PINCH PARAMETER S : ',CVIN)         
      ELSEIF (CVPOPT.EQ.1) THEN                                                
       CALL PRC ('  PINCH OPTION     1 : INWARD PINCH VELOCITY DESCRIBED
     >BY:')
       CALL PRC ('                       VIN = V0*[(CA-X)/CA]^a')
       CALL PRR ('                       V0  = ', VPV0)
       CALL PRR ('                       a   = ', VPALPH)
      ELSEIF (CVPOPT.EQ.2) THEN 
       CALL PRC ('  PINCH OPTION     2 : INWARD PINCH VELOCITY DESCRIBED    
     >BY:')
       CALL PRC ('                       VIN = CV*DP*[(1/Neb)*(dNeb/dx)]
     >')
       CALL PRR ('                       CV  = ', CVIN)
       IF (CVPCUT.NE.99.0) 
     >    CALL PRR ('                       VIN =0 FOR X > ',CVPCUT)
      ENDIF
C
C     Arbitrary Vpinch - can be either in or out
C                                                                     
      IF (CVPOUT.NE.0.0) THEN 
       CALL PRR ('  ARBITRARY PINCH    : AN ADDITIONAL PINCH VELOCITY OF    
     >',CVPOUT)
       CALL PRC ('                       IN M/S HAS BEEN SPECIFIED')
       IF (CVPOUT.LT.0.0) THEN 
         CALL PRC ('                       TOWARDS THE WALL. (OUTWARD -
     >V)')
       ELSE
         CALL PRC ('                       TOWARDS THE CORE. (INWARD +V)
     >')
       ENDIF  
       CALL PRR2('                       IN THE X-RANGE: ',CVXMIN,
     >           CVXMAX)    
      ENDIF
C-----------------------------------------------------------------------        
      IF     (CNEUTA.EQ.0) THEN                                                 
       CALL PRC ('  CONTROL SWITCH   0 : NEUT ON: FOLLOW ATOMS TO IONISA        
     >TION POSNS')                                                              
      ELSEIF (CNEUTA.EQ.1) THEN                                                 
       CALL PRC ('  CONTROL SWITCH   1 : NEUT OFF: INJECT IONS AS "INITI        
     >AL STATE"')                                                               
       IF (CNEUTD.LT.3) GOTO 999
       CALL PRC ('                       SPUTTER,VEL/ANG FLAGS APPLY TO 
     >SELF-SPUTTER ONLY')     
       GOTO 990                                                                 
      ELSEIF (CNEUTA.EQ.2) THEN                                                 
       CALL PRC ('  CONTROL SWITCH   2 : NEUT SIMULATION: INJECT IONS WI        
     >THIN AREA')                                                               
       WRITE (COMENT,'(23X,F7.4,''<X<'',F7.4,'', '',                            
     >   F7.4,''<Y<'',F7.4)') X0S,X0L,Y0S,Y0L                                   
       CALL PRC (COMENT)                                                        
       IF (CNEUTD.LT.3) GOTO 999 
       CALL PRC ('                       SPUTTER,VEL/ANG FLAGS APPLY TO 
     >SELF-SPUTTER ONLY')     
       GOTO 990                                              
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF     (CNEUTB.EQ.0) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    0 : DISTRIBUTED LAUNCH ACROSS LIMIT        
     >ER AT P0')                                                                
       CALL PRR ('                       LAUNCH CUTOFF AT X =', CCUT)           
       IF (CLARMR.GT.0.0)                                                       
     > CALL PRR ('                       LARMOR RADIUS RL (M)', CLARMR)         
      ELSEIF (CNEUTB.EQ.1) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    1 : AT (X0,+/-EDGE(X0),P0) SYMMETRI        
     >CALLY X0<0')                                                              
      ELSEIF (CNEUTB.EQ.2) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    2 : AT (X0,+EDGE(X0),P0) ASYMMETRIC        
     >ALLY, X0<0')                                                              
      ELSEIF (CNEUTB.EQ.3) THEN                                                 
       WRITE (COMENT,'(''  LAUNCH OPTION    3 : AT ('',F6.3,                    
     >  '','',F6.3,'',P0), NORMAL=90, P0 GIVEN'')') CXSC,CYSC                
       CALL PRC (COMENT)                                                        
      ELSEIF (CNEUTB.EQ.4) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    4 : HOMOGENEOUSLY AT WALL,P0, NORMA        
     >L=90 DEGS')                                                               
      ELSEIF (CNEUTB.EQ.5) THEN
       CALL PRC ('  LAUNCH OPTION    5 : DISTRIBUTED LAUNCH ACROSS LIMIT        
     >ER TO +/-P0')                                                           
       CALL PRR ('                       LAUNCH CUTOFF AT X=    ',CCUT)        
       CALL PRR ('                       LAUNCH CUTOFF AT P=+/- ',CPSC)
       IF (CLARMR.GT.0.0)                                                       
     > CALL PRR ('                       LARMOR RADIUS RL (M)', CLARMR)         
      ELSEIF (CNEUTB.EQ.6) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    6 : AT (X0,+/-Y0,P0) RANDOMLY DISTR
     >IBUTED ON Y')       
      ELSEIF (CNEUTB.EQ.7) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    7 : AT (X,Y1-Y2,P0)  RANDOMLY DISTR
     >IBUTED ON Y')       
       CALL PRC ('                       INITIAL X-POS BASED ON IONIZATI
     >ON DATA')
       CALL PRC ('                       INITIAL Y-POS DISTRIBUTED ON (M
     >) : ')
       WRITE(7,'(''                       Y1:Y2 '',G12.5,'':'',G12.5)')
     >           Y0S,Y0L
C      CALL PRR2('                       Y1:Y2 ',Y0S,Y0L)
       CALL PRC ('                       INITIAL ENERGY BY INPUT VEL/ANG 
     >FLAG')
       CALL PRC ('                         X-POS      CUMULATIVE PROBABI
     >LITY')
       DO 120 IX = 1,CLPD
         CALL PRR2('                      ',LPDION(IX,1),LPDCUM(IX))
 120   CONTINUE    
      ELSEIF (CNEUTB.EQ.8) THEN                                                 
       CALL PRC ('  LAUNCH OPTION    8 : AT (X1-X2,Y,P0)  RANDOMLY DISTR
     >IBUTED ON Y')       
       CALL PRC ('                       INITIAL Y-POS BASED ON IONIZATI
     >ON DATA')
       CALL PRC ('                       INITIAL X-POS DISTRIBUTED ON (M
     >) : ')
       WRITE(7,'(''                       X1:X2 '',G12.5,'':'',G12.5)')
     >           X0S,X0L
C      CALL PRR2('                       X1:X2 ',X0S,X0L)
       CALL PRC ('                       INITIAL ENERGY BY INPUT VEL/ANG 
     >FLAG')
       CALL PRC ('                         Y-POS      CUMULATIVE PROBABI
     >LITY')
       DO 130 IX = 1,CLPD
         CALL PRR2('                      ',LPDION(IX,1),LPDCUM(IX))
 130   CONTINUE    
c slmod
      ELSEIF (CNEUTB.EQ.9) THEN
       CALL PRC ('  LAUNCH OPTION    9 : 3D Gaussian neutral lanuch.')
      ELSEIF (CNEUTB.EQ.10) THEN
       CALL PRC ('  LAUNCH OPTION   10 : Real 3D Gaussian launch.')
c slmod end
      ENDIF                                                                     
      IF (CFBGFF.GT.0.0) THEN 
       CALL PRR ('                       CROSS-FIELD FACTOR = ',CFBGFF)
      ENDIF  
C-----------------------------------------------------------------------
C
C  NOTE: IF SPUTTER OPT 3 OR MORE HAS BEEN SELECTED LAUNCH 6,7,8 WILL
C        ALLOW PROPER SPUTTERING WITH THE SELECTED VEL/ANG
C        FLAG TO TAKE PLACE
C-----------------------------------------------------------------------
C------- INITIAL NEUTRAL V/A FLAG --------------------------------------
C
C  THIS OPTION IS SET EQUAL TO CNEUTC, IF REQ'D, IN RUNLM3
C  IMMEDIATELY AFTER RETURN FROM THE READIN ROUTINE 
C
      IF (NVAOPT.NE.CNEUTC)  THEN
      IF     (NVAOPT.EQ.0) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 0 : THETA =+/-ASIN($), $ IN (0,1)')        
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-        
     >1))),')                                                                   
       CALL PRC ('                             $ IN (0,1)')                     
      ELSEIF (NVAOPT.EQ.1) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 1 : THETA = ATAN(TAN(BETA)COS(PHI))        
     >')                                                                        
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(        
     >PHI)**2|)')                                                               
      ELSEIF (NVAOPT.EQ.2) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 2 : THETA = ATAN(TAN(BETA)COS(PHI))        
     >')                                                                        
       CALL PRC ('                       VIN = SQRT(2TG/MI).SQRT(ABS(LN(        
     >1-$))),')                                                                 
       CALL PRR ('                             $ IN (0,1),  TG=',CTGAS)         
      ELSEIF (NVAOPT.EQ.3) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 3 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (NVAOPT.EQ.4) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 4 : THETA = ATAN(TAN(BETA)COS(PHI))        
     >')                                                                        
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(        
     >PHI)**2|)')                                                               
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)                 
      ELSEIF (NVAOPT.EQ.5) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 5 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)                 
      ELSEIF (NVAOPT.EQ.6) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 6 : THETA = 0')                            
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (NVAOPT.EQ.7) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 7 : THETA =+/-ACOS((1-$)**1/3), $ I        
     >N (0,1)   "FREE JET"')                                                    
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (NVAOPT.EQ.8) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 8 : THETA = 2PI$,  $ IN (0,1)  "ISO        
     >TROPIC"')                                                                 
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (NVAOPT.EQ.9) THEN                                                 
       CALL PRC ('  INITIAL V/A FLAG 9 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRC ('                       VIN = SQRT(2EIN/MI), WHERE TWO         
     >VALUES')                                                                  
       CALL PRC ('                       ARE USED ALTERNATELY FOR EIN :-        
     >')                                                                        
       WRITE (COMENT,'(23X,''EIN1 ='',1P,G11.4,''  AND EIN2 ='',G11.4)')        
     >   CENGSC,CEIN2                                                           
       CALL PRC (COMENT)                                                        
      ELSEIF (NVAOPT.EQ.10) THEN                                                
       CALL PRC ('  INITIAL V/A FLAG 10: BETA = ACOS((1-$)**1/3)  "3D FR        
     >EE JET"')                                                                 
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')              
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (NVAOPT.EQ.11) THEN                                                
       CALL PRC ('  INITIAL V/A FLAG 11: THETA =+/-ACOS((1-$)**1/3), 0<$        
     ><1   "2.5D FREE JET"')                                                    
       CALL PRR ('                       P0 RANDOMLY TAKEN IN RANGE +/-'        
     >   ,CPSC)                                                                 
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (NVAOPT.EQ.12) THEN
       CALL PRC ('  INITIAL V/A FLAG 12: BETA = ASIN(($)**1/2)  "TRUE 3D        
     >"')                                                                 
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')              
       CALL PRC ('                       + SURFACE NORMAL       ')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
      ELSEIF (NVAOPT.EQ.13) THEN
       CALL PRC ('  INITIAL V/A FLAG 13: THETA=+/-90.0,'//
     >                                  ' Y-DIRECTION ONLY')
       CALL PRR ('                       EIN1 (EV) = ',CENGSC)
       CALL PRR ('                       EIN2 (EV) = ',CEIN2) 
       CALL PRR ('                       PROBABILITY EIN1 = ',CPROB)
       CALL PRR ('                       PROBABILITY EIN2 = ',1.0-CPROB)  
      ELSEIF (NVAOPT.EQ.14) THEN 
       CALL PRC ('  INITIAL V/A FLAG 14: THETA = 2PI$,  $ IN (0,1)  "ISO        
     >TROPIC"')                                                                 
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-        
     >1))),')                                                                   
       CALL PRC ('                             $ IN (0,$MAX)')             
       CALL PRC ('                       MAXIMUM VELOCITY LIMITATION IN 
     >EFFECT')
      ELSEIF (NVAOPT.EQ.15) THEN 
       CALL PRC ('  INITIAL V/A FLAG 15: THETA =+/-90.0, Y-PLANE ONLY')      
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-        
     >1))),')                                                                   
       CALL PRC ('                             $ IN (0,1)')                     
      ELSEIF (NVAOPT.EQ.16) THEN
       CALL PRC ('  INITIAL V/A FLAG 16: THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRR ('                       EIN1 (EV) = ',CENGSC)
       CALL PRR ('                       EIN2 (EV) = ',CEIN2) 
       CALL PRR ('                       PROBABILITY EIN1 = ',CPROB)
       CALL PRR ('                       PROBABILITY EIN2 = ',1.0-CPROB)
      ELSEIF (NVAOPT.EQ.17) THEN                                         
       CALL PRC ('  INITIAL V/A FLAG 17: THETA = ',CIANGN*raddeg)               
       CALL PRC ('                       XY-PLANE ONLY')                  
       CALL PRR ('                       EIN1 (EV) = ',CENGSC)           
       CALL PRR ('                       EIN2 (EV) = ',CEIN2)            
       CALL PRR ('                       PROBABILITY EIN1 = ',CPROB)     
       CALL PRR ('                       PROBABILITY EIN2 = ',1.0-CPROB) 
      ENDIF
      ENDIF
C-----------------------------------------------------------------------        
990   CONTINUE
C-----------------------------------------------------------------------
      IF     (CNEUTC.EQ.0) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   0 : THETA =+/-ASIN($), $ IN (0,1)')        
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-        
     >1))),')                                                                   
       CALL PRC ('                             $ IN (0,1)')                     
      ELSEIF (CNEUTC.EQ.1) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   1 : THETA = ATAN(TAN(BETA)COS(PHI))        
     >')                                                                        
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(        
     >PHI)**2|)')                                                               
      ELSEIF (CNEUTC.EQ.2) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   2 : THETA = ATAN(TAN(BETA)COS(PHI))        
     >')                                                                        
       CALL PRC ('                       VIN = SQRT(2TG/MI).SQRT(ABS(LN(        
     >1-$))),')                                                                 
       CALL PRR ('                             $ IN (0,1),  TG=',CTGAS)         
      ELSEIF (CNEUTC.EQ.3) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   3 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (CNEUTC.EQ.4) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   4 : THETA = ATAN(TAN(BETA)COS(PHI))        
     >')                                                                        
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(        
     >PHI)**2|)')                                                               
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)                 
      ELSEIF (CNEUTC.EQ.5) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   5 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)                 
      ELSEIF (CNEUTC.EQ.6) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   6 : THETA = 0')                            
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
c slmod begin
       CALL PRC ('      *** NOTE ***     dvy and dvp are set to 0.0.')
c slmod end
      ELSEIF (CNEUTC.EQ.7) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   7 : THETA =+/-ACOS((1-$)**1/3), $ I        
     >N (0,1)   "FREE JET"')                                                    
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (CNEUTC.EQ.8) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   8 : THETA = 2PI$,  $ IN (0,1)  "ISO        
     >TROPIC"')                                                                 
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (CNEUTC.EQ.9) THEN                                                 
       CALL PRC ('  VEL/ANGLE FLAG   9 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRC ('                       VIN = SQRT(2EIN/MI), WHERE TWO         
     >VALUES')                                                                  
       CALL PRC ('                       ARE USED ALTERNATELY FOR EIN :-        
     >')                                                                        
       WRITE (COMENT,'(23X,''EIN1 ='',1P,G11.4,''  AND EIN2 ='',G11.4)')        
     >   CENGSC,CEIN2                                                           
       CALL PRC (COMENT)                                                        
      ELSEIF (CNEUTC.EQ.10) THEN                                                
       CALL PRC ('  VEL/ANGLE FLAG  10 : BETA = ACOS((1-$)**1/3)  "3D FR        
     >EE JET"')                                                                 
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')              
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (CNEUTC.EQ.11) THEN                                                
       CALL PRC ('  VEL/ANGLE FLAG  11 : THETA =+/-ACOS((1-$)**1/3), 0<$        
     ><1   "2.5D FREE JET"')                                                    
       CALL PRR ('                       P0 RANDOMLY TAKEN IN RANGE +/-'        
     >   ,CPSC)                                                                 
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',            
     >    CENGSC)                                                               
      ELSEIF (CNEUTC.EQ.12) THEN
       CALL PRC ('  VEL/ANGLE FLAG  12 : BETA = ASIN(($)**1/2)  "TRUE 3D        
     >"')                                                                 
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')              
       CALL PRC ('                       + SURFACE NORMAL       ')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)        
     >) $<$MAX')                                                                
      ELSEIF (CNEUTC.EQ.13) THEN
       CALL PRC ('  VEL/ANGLE FLAG  13 : THETA = +/-90.0, Y-PLANE ONLY')
       CALL PRR ('                       EIN1 (EV) = ',CENGSC)
       CALL PRR ('                       EIN2 (EV) = ',CEIN2) 
       CALL PRR ('                       PROBABILITY EIN1 = ',CPROB)
       CALL PRR ('                       PROBABILITY EIN2 = ',1.0-CPROB)  
      ELSEIF (CNEUTC.EQ.14) THEN 
       CALL PRC ('  VEL/ANGLE FLAG  14 : THETA = 2PI$,  $ IN (0,1)  "ISO        
     >TROPIC"')                                                                 
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-        
     >1))),')                                                                   
       CALL PRC ('                             $ IN (0,$MAX)')              
       CALL PRC ('                       MAXIMUM VELOCITY LIMITATION IN 
     >EFFECT')
      ELSEIF (CNEUTC.EQ.15) THEN 
       CALL PRC ('  VEL/ANGLE FLAG  15 : THETA =+/-90.0, Y-PLANE ONLY')      
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-        
     >1))),')                                                                   
       CALL PRC ('                             $ IN (0,1)')                     
      ELSEIF (CNEUTC.EQ.16) THEN
       CALL PRC ('  VEL/ANGLE FLAG  16 : THETA =+/-ASIN(SQRT($)), $ IN (        
     >0,1)')                                                                    
       CALL PRR ('                       EIN1 (EV) = ',CENGSC)
       CALL PRR ('                       EIN2 (EV) = ',CEIN2) 
       CALL PRR ('                       PROBABILITY EIN1 = ',CPROB)
       CALL PRR ('                       PROBABILITY EIN2 = ',1.0-CPROB)
      ELSEIF (CNEUTC.EQ.17) THEN                                         
       CALL PRC ('  VEL/ANGLE FLAG  17 : THETA = ',CIANGN*raddeg)
       CALL PRC ('                       Y-PLANE ONLY') 
       CALL PRR ('                       EIN1 (EV) = ',CENGSC)           
       CALL PRR ('                       EIN2 (EV) = ',CEIN2)            
       CALL PRR ('                       PROBABILITY EIN1 = ',CPROB)     
       CALL PRR ('                       PROBABILITY EIN2 = ',1.0-CPROB) 
       ENDIF                                                                     
C-----------------------------------------------------------------------
C
C     NEUTRAL REFLECTION OPTION- REFLECT FROM WALLS IN LAUNCH ROUTINE
C     

      WRITE(6,*) 'NRFOPT: ',NRFOPT

      IF (NRFOPT.EQ.1) THEN 
       CALL PRB
       CALL PRC ('  NEUTRAL REFLECTION OPT : ON - NEUTRALS REFLECTED'      
     >//' AT WALL IMPACT')
       CALL PRB
      ENDIF
C-----------------------------------------------------------------------
C
C BRANCH AROUND THE SPUTTER OPTION FOR LAUNCH OPTIONS 6,7,8 IF THE 
C OPTION SPECIFIED DOES NOT INCLUDE SELF-SPUTTERING 
C
      IF (CNEUTB.EQ.3.OR.CNEUTB.EQ.6.OR.CNEUTB.EQ.7.OR.
c slmod
     >    CNEUTB.EQ.8.OR.CNEUTB.EQ.9.OR.CNEUTB.EQ.10) THEN 
c     >    CNEUTB.EQ.8) THEN 
c slmod end
         IF (CNEUTD.LT.3) THEN 
            GOTO 998
         ELSE
           CALL PRB
           CALL PRC ('  SPUTTER OPTION APPLIES TO SELF-SPUTTER ONLY')     
           CALL PRB  
         ENDIF
      ENDIF  
c
C-----------------------------------------------------------------------

      call prc(' PHYSICAL SPUTTERING DATA SOURCE OPTION:')
      IF     (CSPUTOPT.EQ.1) THEN
       CALL PRC ('  SPUTTER SOURCE   1 : Formulation due to Bohdansky')
       CALL PRC ('                       (modified coefficients)')
      ELSEIF (CSPUTOPT.EQ.2) THEN
       CALL PRC ('  SPUTTER SOURCE   2 : Eckstein IPP9/82 (1993)')
      ELSEIF (CSPUTOPT.EQ.3) THEN
       CALL PRC ('  SPUTTER SOURCE   3 : Based on Eckstein IPP9/82 (1993
     >)')
       call prc ('                       Slight changes to H,D,T on C co
     >efficients')
       call prc ('                       From Garcia/Rosales-Roth 1996')
       call prc ('                       Slight changes to D on W co'//
     >'efficients')
      elseif (csputopt.eq.4) then
       CALL PRC ('  SPUTTER SOURCE   4 : Specified CONSTANT Yield Value'
     >)
       call prr ('                       Yield = ',const_yield)
      ELSEIF (CSPUTOPT.EQ.5) THEN
       CALL PRC ('  SPUTTER SOURCE   5 : Defaults to SPUTTER SOURCE'//
     >       ' 3 except for cases with custom data specified instead')
       call prc ('                       Default:'//
     >                            'Based on Eckstein IPP9/82 (1993)')
       call prc ('                       Slight changes to H,D,T on C co
     >efficients')
       call prc ('          List of custom data sets:')
       call prc ('          W Custom Yield Data : K. Krieger')
       call prc ('          D-> Be and Be-> Be  : Ecsktein 2002') 
       call prc ('          D->C and C->C       : Eckstein 2002')
       if ((cion.eq.4.or.cion.eq.6).and.crmb.eq.2.0) then 
       ! Be or C selected
          if (cion.eq.4) then 
             call prc('     BERYLLIUM SPUTTERING YIELD DATA SELECTED:')
          elseif (cion.eq.6) then 
             call prc('     CARBON SPUTTERING YIELD DATA SELECTED:')
          endif
          if (extra_sputter_angle.lt.0.0) then
             call prc('     - ANGLE AVERAGED DATA SELECTED')
          else 
             call prr('     - DATA SELECTED FOR INCIDENT ANGLE =',
     >                  extra_sputter_angle)
          endif
          call print_eck2002_yields(7)
       elseif (cion.eq.74) then 
          ! W selected 
             CALL PRC ('    TUNGSTEN SPUTTERING DATA SELECTED:')
       endif

      ENDIF
C-----------------------------------------------------------------------        
      IF     (CNEUTD.EQ.0) THEN                                                 
       CALL PRC ('  SPUTTER OPTION   0 : SPUTTERING BY BACKGROUND IONS O        
     >NLY')                                                                     
       CALL PRC ('                       EIMP=TB(2+3ZB)')                       
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')         
      ELSEIF (CNEUTD.EQ.1) THEN                                                 
       CALL PRC ('  SPUTTER OPTION   1 : SPUTTERING BY SPECIFIED ION TYP        
     >E ONLY')                                                                  
       CALL PRI ('                       EIMP=TB(2+3ZIMP),  ZIMP=',             
     >   CBOMBZ)                                                                
       CALL PRC ('                       EMAX=EIMP')                            
      ELSEIF (CNEUTD.EQ.2) THEN                                                 
       CALL PRC ('  SPUTTER OPTION   2 : SPUTTERING BY MIXTURE OF BACKGR        
     >OUND IONS')                                                               
       CALL PRR ('                       AND SPECIFIED ION TYPE, FIMP=',        
     >    CFIMP)                                                                
       CALL PRC ('                       AND FLUX2 TAKEN AS FLUX1.FIMP')        
       CALL PRC ('                       EIMP1=TB(2+3ZB)')                      
       CALL PRC ('                       EMAX1=EIMP1.GAMMA(1-GAMMA)-EBD'        
     >)                                                                         
       CALL PRI ('                       EIMP2=TB(2+3ZIMP),  ZIMP=',            
     >   CBOMBZ)                                                                
       CALL PRC ('                       EMAX2=EIMP2')                          
      ELSEIF (CNEUTD.EQ.3) THEN                                                 
       CALL PRC ('  SPUTTER OPTION   3 : INITIAL SPUTTERING BY BACKGND I        
     >ONS ONLY')                                                                
       CALL PRC ('                       EIMP=TB(2+3ZB)')                       
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')         
       CALL PRC ('                       PROPER SELF-SPUTTERING USING AC        
     >TUAL ZIMP')                                                               
       CALL PRC ('                       VALUES ON EXIT AND SAME VEL/ANG        
     > FLAG.')                                                                  
       CALL PRC ('                       EIMP=3TB.ZIMP+5.22E-9.MI.VEXIT.        
     >VEXIT+2TI')                                                               
       CALL PRR ('                       EMAX=EIMP.  THRESHOLD YIELD=',         
     >   CTRESH)                                                                
       CALL PRR ('                       MAX ALLOWED YIELD ANY FRAG.=', 
     >   CSPUMAX)
      ELSEIF (CNEUTD.EQ.4) THEN                                                 
       CALL PRC ('  SPUTTER OPTION   4 : INITIAL SPUTTERING BY BACKGND I        
     >ONS ONLY')                                                                
       CALL PRC ('                       EIMP=TB(2+3ZB)')                       
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')         
       CALL PRC ('                       PROPER SELF-SPUTTERING USING AC        
     >TUAL ZIMP')                                                               
       CALL PRC ('                       VALUES ON EXIT, AND WITH VEL/AN        
     >G FLAG 1')                                                                
       CALL PRC ('                       EIMP=3TB.ZIMP+5.22E-9.MI.VEXIT.        
     >VEXIT+2TI')                                                               
       CALL PRR ('                       EMAX=EIMP.  THRESHOLD YIELD=',         
     >   CTRESH)                                                                
       CALL PRR ('                       MAX ALLOWED YIELD ANY FRAG.=', 
     >   CSPUMAX)
      ELSEIF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN               
       WRITE (7,'(1X,A,I1,A)') '  SPUTTER OPTION   ',CNEUTD,                    
     >   ' : RADIATION ENHANCED SUBLIMATION'                                    
       CALL PRC ('                       V/A FLAG 3 WITH EIN=0.15 FOR RE        
     >S')                                                                       
       CALL PRR ('                       EIMP=TB(2+3ZB),  TSUB=',CTSUB)         
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')         
       CALL PRR ('                       QPL=122EXP(-9048/TSUB)=',CQPL)         
       CALL PRC ('                       RES FOR $ < YPL/(YPH+YPL), $ IN        
     > (0,1)')                                                                  
       CALL PRC ('                       SELF-SPUTTERING USING ACTUAL ZI        
     >MP VALS')                                                                 
       CALL PRC ('                       EIMP=3TB.ZIMP+5.22E-9.MI.VEXIT.        
     >VEXIT+2TI')                                                               
       CALL PRR ('                       EMAX=EIMP.  THRESHOLD YIELD=',         
     >   CTRESH)                                                                
       CALL PRR ('                       MAX ALLOWED YIELD ANY FRAG.=', 
     >   CSPUMAX)
       CALL PRR ('                       QSL=1014EXP(-9048/TSUB)=',CQSL)        
       CALL PRC ('                       RES FOR $ < YSL/(YSH+YSL), $ IN        
     > (0,1)')                                                                  
       IF (CNEUTD.EQ.6) WRITE (7,'(23X,A)') ' ALL HIGH VELOCITIES'              
       IF (CNEUTD.EQ.7) WRITE (7,'(23X,A)') ' ALL LOW VELOCITIES'               
      ELSEIF (CNEUTD.EQ.8) THEN 
       CALL PRC ('  SPUTTER OPTION   8 : INITIAL SPUTTERING BY BACKGND I        
     >ONS ONLY')                                                                
       CALL PRC ('                       EIMP=2TB+ZB.VS)')                    
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')         
       CALL PRC ('                       PROPER SELF-SPUTTERING USING AC        
     >TUAL ZIMP')                                                               
       CALL PRC ('                       VALUES ON EXIT AND SAME VEL/ANG        
     > FLAG.')                                                                  
       CALL PRC ('                       EIMP=ZIMP.VS +5.22E-9.MI.VEXIT.        
     >VEXIT+2TI')                                                               
       CALL PRR ('                       EMAX=EIMP.  THRESHOLD YIELD=',         
     >   CTRESH)                                                                
       CALL PRR ('                       MAX ALLOWED YIELD ANY FRAG.=', 
     >   CSPUMAX)
      ENDIF                                                                     
      CALL PRR ('                       SPUTTERING ENHANCEMENT FACTOR',         
     >   CSEF)                                                                  
      IF (QMULTP.NE.1.0.OR.QMULTS.NE.1.0) THEN 
       CALL PRR('                       Q - MULT    (PRIMARIES)    = ', 
     >     QMULTP)
       CALL PRR('                       Q - MULT    (SECONDARIES)  = ', 
     >     QMULTS)
      ENDIF
c
C-----------------------------------------------------------------------        
c
c     Sputtering particle impact energy option
c
C-----------------------------------------------------------------------        
c
      call prb
c
      if (impact_energy_opt.eq.0) then 
       CALL PRC ('  IMPACT ENERGY OPT 0: IMPURITY PARTICLE IMPACT ENERGY
     > FOR')
       call prc ('                       SELF-SPUTTERING IS CALCULATED A
     >S')
       CALL PRC ('                       STATED IN THE SPUTTER OPTION')
      elseif (impact_energy_opt.eq.1) then
       CALL PRC ('  IMPACT ENERGY OPT 1: IMPURITY PARTICLE IMPACT ENERGY
     > FOR')
       call prc ('                       SELF-SPUTTERING IS CALCULATED A
     >S:')
       CALL PRC ('                       EIMP = ZIMP*3.0*TIB + 2.0*TI')
       call prc ('                       THIS OVER-RIDES ALL SPUTTER OPT
     >IONS - EXCEPT 8')
      endif

C
C-----------------------------------------------------------------------        
C 
C NORMAL OPTION :
C
C-----------------------------------------------------------------------
 998  CONTINUE
      IF     (CNEUTE.EQ.0) THEN                                                 
       CALL PRC ('  NORMAL OPTION    0 : MEASURE THETA FROM SURFACE NORM        
     >AL')                                                                      
      ELSEIF (CNEUTE.EQ.1.OR.CNEUTE.EQ.2) THEN                                  
       WRITE (COMENT,'(A,I1,A,F7.2,A)') '  NORMAL OPTION    ',CNEUTE,          
     >   ' : MEASURE THETA FROM ',CSNORM,' DEGS TO X=0'                         
       CALL PRC (COMENT)                                                        
       IF (CNEUTE.EQ.1) WRITE (7,'(23X,A)') 'FOR PRIMARIES ONLY'                
       IF (CNEUTE.EQ.2) WRITE (7,'(23X,A)') 'FOR PRIMARIES & SELF-SPUT.'        
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF     (CNEUTF.EQ.0) THEN                                                 
       CALL PRC ('  NEUT SPREADING   0 : OFF (LAUNCH AT MESH PTS ONLY)')        
      ELSEIF (CNEUTF.EQ.1) THEN                                                 
       CALL PRC ('  NEUT SPREADING   1 : ON  (LAUNCH BETWEEN MESH PTS)')        
      ENDIF                                                                     
C-----------------------------------------------------------------------        
C                                                                               
C---- CHECK FOR DUBIOUS / NON-IMPLEMENTED  COMBINATIONS ...                     
C                                                                               
  999 CONTINUE                                                                  
      IF (CIOPTI.NE.0.AND.CIZB.NE.1) THEN                                       
       CALL PRB                                                                 
       CALL PRC ('*** WARNING *** C-X RECOM SPECIFIED FOR NON-HYDROGENIC        
     > PLASMA...')                                                              
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
 9010 FORMAT(1X,A,F9.2)                                                         
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE DMPOUT (TITLE,NIZS,NOUT,IERR,JOB,IMODE,PLAMS,PIZS,NLS,        
     >                 FACTA,FACTB,ITER,NITERS)                                 
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  DUMP:  STORE RESULTS OF LIM RUN IN UNFORMATTED FILE "NOUT".      *        
C  *         THIS FILE SHOULD BE REWOUND BEFORE DUMP IS CALLED FOR THE *        
C  *         FIRST TIME.                                               *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      IMPLICIT  none
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'dynam1'                                                        
C     INCLUDE   (DYNAM1)                                                        
      INCLUDE   'dynam3'                                                        
C     INCLUDE   (DYNAM3)                                                        
      
      INCLUDE   'cnoco' 

      CHARACTER TITLE*80,JOB*72                                                 
      INTEGER   NIZS,IMODE,NLS                                                  
      REAL      PLAMS(MAXNLS),FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)                 
      INTEGER   NOUT,IERR,PIZS(MAXNLS),ITER,NITERS,JY                      
C                                                                               
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
      INCLUDE   'comnet'                                                        
C     INCLUDE   (COMNET)                                                        
c
c     jdemod - include ADAS to calculate the 3D POWL and LINE arrays
      include   'cadas'
c
c slmod tmp
      INCLUDE 'comtau'
c slmod end
C                                                                               
C----- DUMP DATASET SHOULD HAVE FORMAT U, RECORD LENGTH 0, BLOCKSIZE            
C----- 6160, PREFERABLY ORGANISED AS A PARTITIONED DATASET.  HENCE              
C----- 1540 4-BYTE WORDS WILL FIT IN EACH BLOCK.                                
C                                                                               
      INTEGER   IOS,JBLOCK,IL,II,IP,J,KBLOCK,IT,IO                              
      INTEGER   IBLOCK,IQX,IX,IY,IZ,IYB,IYE,IZS,IZE,IQS,IQE                     
      DATA      IBLOCK / 1540 /                                                 
      !
      ! jdemod - initialization resulting from complaints by INTEL compiler
      ! 
      ios = 0
      !
c slmod tmp
      IF (CIOPTE.EQ.10) THEN
        WRITE (NOUT,IOSTAT=IOS) 'EMPTY'
        RETURN
      ENDIF
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C     WRITE COMMONS, LINE DATA, SHORT ARRAYS, PARAMETERS, ETC                   
C-----------------------------------------------------------------------        
C                                                                               
C---- FIRST ITEM IS VERSION NUMBER (EG 3I/01) TO BE READ BACK                   
C---- IN COLECT - IF AGREEMENT IS NOT FOUND THEN AN ERROR MESSAGE IS            
C---- ISSUED.                                                                   
C                                                                               
c      write(0,*) 'VERSION:',':',ios,':',trim(verson),':'
c
      WRITE (NOUT,IOSTAT=IOS) VERSON,NY3D,ITER,NITERS,MAXOS                     
      WRITE (NOUT,IOSTAT=IOS)                                                   
     >        NXS,NYS,NQXSO,NQXSI,NQYS,NTS,NIZS,NLS,TITLE,JOB,IMODE             
      WRITE (6,9002)                                                            
     >        NXS,NYS,NQXSO,NQXSI,NQYS,NTS,NIZS,NLS,TITLE,JOB,IMODE,ITER        
c
c     Save the scaling factor for the case - the scaling factor is
c     kept in DEFACT - however, absfac was added for compatibility 
c     with DIVIMP code. ABSFAC was added to the comtor common include
c     file. DEFACT is assigned to ABSFAC in LIM3.L3I
c
      write(NOUT,IOSTAT=IOS)  absfac
C                                                                               
C---- SOME OF COMTOR COMMON, SOME OF COMXYT COMMON, MISCELLANEOUS               
C                                                                                      
      WRITE (NOUT,IOSTAT=IOS)                                                   
     >       CA    ,CAW   ,CL    ,CRMB  ,CTBIN ,CGTIN1,CLTIN1,         
     >       CTBOUL,CLTOUL,CNBIN ,CGNIN1,CLNIN1,CNBOUL,CLNOUL,         
     >       CTIBIN ,CGTIIN1,CLTIIN1,CTIBOUL,CLTIOUL,CSTEPT,         
     >       CVIN  ,CYFAR ,CRDD  ,CRMI  ,CXSC  ,CFIMP ,CONO  ,         
     >       CYSC  ,CTEMSC,CTIMSC,CEBD  ,CTGAS ,CEMAXF,CENGSC,         
     >       CIZB  ,CION  ,CIZSC ,CATIN ,CANIN ,CTRESH,CIZEFF,         
     >       CNHC  ,CNHO  ,CLAMHX,CLAMHY,CXNEAR,CYNEAR,CKO   ,         
     >       CHALFL,CKI   ,CCUT  ,CPFIR ,CPSUB ,CPSC  ,CSNORM,         
     >       CDPOL ,CPLSMA,CEYIN ,CVHYIN,CSTEPN,CIZSET,CZENH ,         
     >       CEYOUT,CVHOUT,CISEED,CDIFOP,CTWOL ,CYSTAG,CONI  ,         
     >       XSCALO,XSCALI,YSCALE,CANAL ,CTHETB,CSINTB,CLFACT,CFBGFF, 
     >       (XS(IX),IX=1,NXS),(YS(IY),IY=1,NYS),                      
     >       (DWELTS(IZ),IZ=0,NIZS),(DWELFS(IT),IT=1,NTS),             
     >       (XWIDS(IX),IX=1,NXS),(YWIDS(IY),IY=1,NYS),                
     >       (PS(IP),IP=-MAXNPS,MAXNPS),(XOUTS(IX),IX=1,NXS),          
     >       (PWIDS(IP),IP=-MAXNPS,MAXNPS),                            
     >       (YOUTS(IY),IY=-NYS,NYS),                                  
     >       (PIZS(IL),IL=1,NLS),(PLAMS(IL),IL=1,NLS)                  
c slmod
     >       ,CLNIN2,CLTIIN2,CLTIN2,CVPOL
c slmod end

c
c     Write out some 3D option information
c
c     CIOPTJ= 3D limiter extent option 
c     CPCO  = 3D extent of limiter
c
c
      write(nout) cioptj,cpco

C                                                                               
C---- SHORT ARRAYS ... BLOCKED I/O USED                                         
C                                                                               
      DO 45 IZS = -2, NIZS+1, 15                                                
        IZE = MIN (IZS+15-1, NIZS+1)                                            
        WRITE (NOUT,IOSTAT=IOS) ((SAVES(IX,IZ),IX=1,NXS),IZ=IZS,IZE)            
   45 CONTINUE                                                                  
C                                                                               
      DO 50 J = 1, 3                                                            
       DO 50 IZS = 1, NIZS+1, 15                                                  
        IZE = MIN (IZS+15-1, NIZS+1)                                              
        WRITE (NOUT,IOSTAT=IOS) ((DEPS(IX,IZ,J),IX=1,NXS),IZ=IZS,IZE)           
   50 CONTINUE                                                                  
C                                                                               
      WRITE (NOUT,IOSTAT=IOS)                                                   
     >  (((NEROXS(IX,II,J),IX=1,NXS),II=1,5),J=1,3)                             
C                                                                               
      DO 55 II = 1, 6                                                           
        WRITE (NOUT,IOSTAT=IOS) (NEROYS(IO,II),IO=1,MAXOS)                      
   55 CONTINUE                                                                  
c
      DO 57 II = 1, 5                                                           
        WRITE (NOUT,IOSTAT=IOS) (NERODS(IO,II),IO=1,MAXOS)                      
   57 CONTINUE                                                                  
c
      DO II = 1, 6                                                           
         do ip = -maxnps,maxnps
            WRITE (NOUT,IOSTAT=IOS) (NERODS3(IO,IP,II),IO=1,MAXOS)                      
         end do
      end do                                                                               
C                                                                               
      DO 60 IZS = -2, NIZS+1, 7                                                 
        IZE = MIN (IZS+7-1, NIZS+1)                                             
        WRITE (NOUT,IOSTAT=IOS) ((WALLS(IY,IZ),IY=-NYS,NYS),IZ=IZS,IZE)         
   60 CONTINUE                                                                  
C                                                                               
      WRITE (NOUT,IOSTAT=IOS) (OYS(IO),ODS(IO),IO=1,MAXOS)                      
      WRITE (NOUT,IOSTAT=IOS) (OYOUTS(IO),ODOUTS(IO),IO=1,MAXOS)                
      WRITE (NOUT,IOSTAT=IOS) (OYWIDS(IO),ODWIDS(IO),IO=1,MAXOS)                
      DO 62 II = 1,3
         WRITE(NOUT,IOSTAT=IOS) (CDFLUX(IO,II),IO=1,MAXOS)
62    CONTINUE
C                                                                               
C---- CHECK FOR ERRORS                                                          
C                                                                               

      IF (IOS.NE.0) STOP                                                        
      JBLOCK = IBLOCK / (NXS+1)                                                 
      KBLOCK = JBLOCK / 2                                                       
C                                                                               
C-----------------------------------------------------------------------        
C     THE ARRAYS MUST BE WRITTEN IN BLOCKS TO PREVENT END-OF                    
C     RECORD ERRORS.  THE NUMBER OF WORDS WRITTEN IN EACH BLOCK WILL            
C     BE TAKEN AS NXS * (IBLOCK/(NXS+1))   (IF THIS IS EXACTLY 1540             
C     END-OF-RECORD ERRORS RESULT, HENCE USE OF NXS+1 IN DENOMINATOR).          
C     FOR DOUBLE PRECISION ARRAYS,                                              
C     A CONVERSION TO SINGLE PRECISION IS MADE FIRST.                           
C-----------------------------------------------------------------------        
C                                                                               
C                                                                               
C================= WRITE DDLIMS ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 100 IZ = -1, NIZS                                                     
        DO 100 IYB = -NYS, NYS, KBLOCK                                          
         IYE = MIN (IYB+KBLOCK-1, NYS)                                          
         WRITE (NOUT) ((SNGL(DDLIMS(IX,IY,IZ)), IX=1,NXS), IY=IYB,IYE)          
  100  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE DDTS ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO IZ =  1, NIZS                                                     
        DO IYB = -NYS, NYS, KBLOCK                                          
         IYE = MIN (IYB+KBLOCK-1, NYS)                                          
         WRITE (NOUT) ((SNGL(DDTS(IX,IY,IZ)), IX=1,NXS), IY=IYB,IYE)          
        end do 
       end do
      ENDIF                                                                     
C                                                                               
C================= WRITE POWLS  ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 200 IZ = -1, NIZS                                                     
        DO 200 IYB = -NYS, NYS, JBLOCK                                          
         IYE = MIN (IYB+JBLOCK-1, NYS)                                          
         WRITE (NOUT) ((POWLS(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                 
  200  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE LINES  ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 300 IZ = -1, NIZS                                                     
        DO 300 IYB = -NYS, NYS, JBLOCK                                          
         IYE = MIN (IYB+JBLOCK-1, NYS)                                          
         WRITE (NOUT) ((LINES(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                 
  300  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE TIZS   ARRAY TO DISC ==========================        
C                                                                               
      DO 500 IZ = -1, NIZS                                                      
       DO 500 IYB = -NYS, NYS, JBLOCK                                           
        IYE = MIN (IYB+JBLOCK-1, NYS)                                           
        WRITE (NOUT) ((TIZS(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                   
  500 CONTINUE                                                                  
C                                                                               
C================= WRITE ZEFFS  ARRAY TO DISC ==========================        
C                                                                               
      DO 600 II = 1, 6                                                          
       DO 600 IYB = -NYS, NYS, JBLOCK                                           
        IYE = MIN (IYB+JBLOCK-1, NYS)                                           
        WRITE (NOUT) ((ZEFFS(IX,IY,II), IX=1,NXS), IY=IYB,IYE)                  
  600 CONTINUE                                                                  
C                                                                               
C================= WRITE DDLIM3 ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 700 IP = -MAXNPS, MAXNPS                                              
        DO 700 IZ = -1, NIZS                                                    
         DO 700 IYB = -NY3D, NY3D, KBLOCK                                       
          IYE = MIN (IYB+KBLOCK-1, NY3D)                                        
          WRITE (NOUT) ((SNGL(DDLIM3(IX,IY,IZ,IP)),IX=1,NXS),IY=IYB,IYE)        
  700  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE TIZ3   ARRAY TO DISC ==========================        
C                                                                               
      IF (CDPOL.GT.0.0) THEN                                                    
       DO 1100 IP = -MAXNPS, MAXNPS                                             
        DO 1100 IZ = -1, NIZS                                                   
         DO 1100 IYB = -NY3D, NY3D, JBLOCK                                      
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          WRITE (NOUT) ((TIZ3(IX,IY,IZ,IP), IX=1,NXS), IY=IYB,IYE)              
 1100  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE LIM5   ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.2) THEN                                                      
       DO 1200 IT = 1, NTS                                                      
        DO 1200 IP = -MAXNPS, MAXNPS                                            
         DO 1200 IZ = -1, NIZS                                                  
          DO 1200 IYB = -NY3D, NY3D, JBLOCK                                     
           IYE = MIN (IYB+JBLOCK-1, NY3D)                                       
           WRITE (NOUT) ((LIM5(IX,IY,IZ,IP,IT), IX=1,NXS), IY=IYB,IYE)          
 1200  CONTINUE                                                                 
      ENDIF                                                                     


C
C                                                                               
C================= WRITE PLRPS  ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 400 IL = 1, NLS                                                       
        DO 400 IYB = -NYS, NYS, JBLOCK                                          
         IYE = MIN (IYB+JBLOCK-1, NYS)                                          
         WRITE (NOUT) ((PLRPS(IX,IY,IL), IX=1,NXS), IY=IYB,IYE)                 
  400  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE PLRP3  ARRAY TO DISC ==========================        
C                                                                               
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 1000 IP = -MAXNPS, MAXNPS                                             
        DO 1000 IL = 1, NLS                                                     
         DO 1000 IYB = -NY3D, NY3D, JBLOCK                                      
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          WRITE (NOUT) ((PLRP3(IX,IY,IL,IP), IX=1,NXS), IY=IYB,IYE)             
 1000  CONTINUE                                                                 
      ENDIF                  
C
C     CALCULATE THE POWL3 AND LINE3 ARRAYS - USE THE TIZ3 AND LIM5
C     ARRAYS FOR TEMPORARY STORAGE.
C     THIS PROCEDURE IS OK BECAUSE THE RESULTS IN THESE ARRAYS ARE 
C     NOT CURRENTLY USED IN ITERATIVE CALCULATIONS. IF THIS CHANGES
C     THIS CODE WILL HAVE TO BE MODIFIED.
C
C     THE ENTIRE REASON BEHIND THESE EFFORTS IS THE CRAY MEMORY 
C     RESTRICTIONS, WHICH CURRENTLY PROHIBIT MOST 3D CASES FROM  
C     RUNNING. THE 3.5 MW (APPROXIMATE) LIMIT IS EXCEEDED FOR 
C     REASONABLY SIZED 3D CASES
C 
C     D. ELDER JAN 29 / 1990
C
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE RADIATIVE POWER LOSS AND LINE RADIATION                         
C-----------------------------------------------------------------------        
C                                                                               
C---- ZERO ARRAYS.                                                              
C                                                                               
      CALL RZERO (TIZ3,  MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
      CALL RZERO (LIM5,  MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
C                                                                               
C-----------------------------------------------------------------------        
C   FOR EACH Y POINT IN TURN ....                                               
C-----------------------------------------------------------------------        
C                                                                               

      if (calc_3d_power.eq.1) then 

!
!     Nocorona
!
      if (cdatopt.eq.0) then 

      DO 990 IY = -NY3D, NY3D                                                
        IF (IY.EQ.0) GOTO 990                                                   
        JY = IABS (IY)                                                          
C                                                                               
C---- CALCULATE PLASMA TEMPERATURE AND ELECTRON DENSITY AT MID POINTS           
C---- OF EACH X BIN.  NOTE THIS INVOLVES A CONVERSION FROM THE REGULAR          
C---- SPACED QXS MESH TO THE USER SUPPLIED XS MESH; AND A CONVERSION            
C---- FROM M**3 TO CM**3 FOR NOCORONA.                                          
C---- THE IQX --> IX INDICES ARE TAKEN FROM COMMON /COMXYT/                     
C                                                                               
      DO 905 IX = 1, NXS                                                        
        PTES(IX) = CTEMBS(IX,IY)                                                
        PNES(IX) = CRNBS(IX,IY) * 1.E-6 * REAL (CIZB)                           
  905 CONTINUE                                                                  
C                                                                               
C                                                                               
C------ CALCULATE BIN VOLUMES CM**3 (ASSUME 1 METRE IN THIRD DIMENSION)         
C------ TRANSFER IMPURITY ION DENSITIES TO PNZS ARRAY FOR NOCORONA              
C------ CONVERTING TO CM**-3                                                    
C                                                                               
        DO 920 IX = 1, NXS                                                      
          PDVOLS(IX) = 1.0E6 * XWIDS(IX) * YWIDS(JY)                            
  920   CONTINUE                                                                
C                                                                               
C------ CALL ROUTINE FROM NOCORONA PACKAGE TO CALCULATE RADIATIVE               
C------ POWER LOSS (W CM**3) AND LINE RADIATION LOSS (W CM**3)                  
C------ THESE ARE CONVERTED TO UNITS (W M**3)                                   
C------ VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY              
C------ GOES OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.         
C                                                                               
C                                                                               
C------ EXTRA SECTION FOR 3D ARRAYS POWL3 AND LINE3 ...                         
C                                                                              
          DO 985 IP = -MAXNPS, MAXNPS                                           
            DO 970 IX = 1, NXS                                                  
              DO 970 IZ = 0, NIZS                                               
                PNZS(IZ+1,1,IX) = 1.0E-6 * SNGL(DDLIM3(IX,IY,IZ,IP))            
  970       CONTINUE                                                            
            CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NXS)                      
            DO 975 IX = 1, NXS                                                  
             DO 975 IZ = 0, NIZS                                               
C              POWL3
               TIZ3(IX,IY,IZ,IP)=MAX(0.0, PRADIS(1,IZ+1,1,IX)*1.0E6)          
C              LINE3
               LIM5(IX,IY,IZ,IP,1)=MAX(0.0, PRADIS(3,IZ+1,1,IX)*1.0E6)          
  975       CONTINUE                                                            
            DO 980 IX = 1, NXS                                                  
              IF (DDLIM3(IX,IY,0,IP).LE.0.0D0) THEN                             
C               POWL3
                TIZ3(IX,IY,-1,IP) = 0.0                                        
C               LINE3
                LIM5(IX,IY,-1,IP,1) = 0.0                                      
              ELSE                                                              
C               POWL3
                TIZ3(IX,IY,-1,IP) = TIZ3(IX,IY,0,IP) /                        
     >            SNGL(DDLIM3(IX,IY,0,IP)) * SNGL(DDLIM3(IX,IY,-1,IP))          
C               LINE3
                LIM5(IX,IY,-1,IP,1) = LIM5(IX,IY,0,IP,1) /                     
     >            SNGL(DDLIM3(IX,IY,0,IP)) * SNGL(DDLIM3(IX,IY,-1,IP))          
              ENDIF                                                             
  980       CONTINUE                                                            
  985     CONTINUE                                                              
  990 CONTINUE                                                                  


!
!     ADAS
!
      elseif (cdatopt.eq.1) then

C
C------ GET POWER LOSS FROM ADAS DATA FILES. LOAD TOTAL LINE RADIATION
C------ INTO LINES AND ADD RECOMBINATION AND BREMSSTRAHLUNG POWER TO
C------ GET TOTAL RADIATIVE LOSSES
C


        write(year,'(i2.2)') iyearz
        call xxuid(useridz)
c
      DO 1190 IY = -NYS,NYS
C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
        IF (IY.EQ.0) GOTO 1190                                                   
        JY = IABS (IY)                                                          
C                                                                               

        DO 1101 IX = 1, NXS   
          PTESA(IX) = CTEMBS(IX,IY)
          PNESA(IX) = CRNBS(IX,IY) * real(cizb)
          PNBS(IX) =  CRNBS(IX,IY)
c
c         Set hydrogen density to zero for now - not available in LIM
c          PNHS(IX) = pinaton(ik,ir)
          pnhs(ix) = 0.0

 1101  CONTINUE


        do ip = -maxnps,maxnps
C

        DO 1120 IX = 1, NXS
          DO 1110 IZ = 0, NIZS
            PNZSA(IX,IZ) = SNGL(DDLIM3(IX,IY,IZ,IP))
 1110     CONTINUE
 1120   CONTINUE



        ICLASS = 5
        !MIZS = MIN(CION-1,NIZS)
        DO 1130 IZ = 0, NIZS
          CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NXS,PTESA,PNESA,
     +                PCOEF(1,IZ+1))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ+1))

          DO 1135 IX = 1, NXS
            ! LINE3 -> TIZ3
            TIZ3(IX,IY,IZ,IP) = PCOEF(IX,IZ+1)*PNESA(IX)*PNZSA(IX,IZ)

            ! POWL3 -> LIM5
            LIM5(IX,IY,IZ,IP,1) = TIZ3(IX,IY,IZ,IP)
c
c            write (6,'(a,3i5,3g16.8)') 'Debug DIV:',ir,ik,iz,
c     >              pcoef(ik,iz+1),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      DIV:',
c     >            lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
 1135     CONTINUE
 1130   CONTINUE


        ICLASS = 4
        !MIZS = MIN(CION,NIZS)
        DO 1140 IZ = 1, NIZS
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NXS,PTESA,PNESA,
     +                PCOEF(1,IZ))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ))
          DO 1145 IX = 1, NXS
            LIM5(IX,IY,IZ,IP,1) = LIM5(IX,IY,IZ,IP,1)
     +                        + PCOEF(IX,IZ)*PNESA(IX)*PNZSA(IX,IZ)
 1145     CONTINUE
 1140   CONTINUE
C
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1)
C
        DO 1160 IX = 1, NXS
          IF (DDLIMS(IX,IY,0).LE.0.0) THEN
            TIZ3(IX,IY,-1,IP) = 0.0
            LIM5(IX,IY,-1,IP,1) = 0.0
          ELSE
            LIM5(IX,IY,-1,IP,1) = LIM5(IX,IY,0,IP,1) *
     +                    SNGL(DDLIM3(IX,IY,-1,IP) / DDLIM3(IX,IY,0,IP))
            TIZ3(IX,IY,-1,IP) = TIZ3(IX,IY,0,IP) *
     +                    SNGL(DDLIM3(IX,IY,-1,IP) / DDLIM3(IX,IY,0,IP))
c
c            write (6,'(a,3i5,3g16.8)') 'Debug POW:',ir,ik,iz,
c     >              pcoef(ik,iz),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      POW:',
c     >           lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
          ENDIF
 1160   CONTINUE

      end do ! end of IP loop

 1190 CONTINUE



      endif

      endif  ! endif for calc_3d_power

                                                   
C                                                                               
C================= WRITE POWL3  ARRAY TO DISC ==========================        
C                  CONTENTS IN TIZ3                                         
C
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 800 IP = -MAXNPS, MAXNPS                                              
        DO 800 IZ = -1, NIZS                                                    
         DO 800 IYB = -NY3D, NY3D, JBLOCK                                       
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          WRITE (NOUT) ((TIZ3(IX,IY,IZ,IP), IX=1,NXS), IY=IYB,IYE)             
  800  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE LINE3  ARRAY TO DISC ==========================        
C                  CONTENTS IN LIM5                                            
C
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 900 IP = -MAXNPS, MAXNPS                                              
        DO 900 IZ = -1, NIZS                                                    
         DO 900 IYB = -NY3D, NY3D, JBLOCK                                       
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          WRITE (NOUT) ((LIM5(IX,IY,IZ,IP,1), IX=1,NXS), IY=IYB,IYE)         
  900  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= WRITE SDTXS,SDTYS,SDYXS,SDYYS =======================        
C                                                                               
      DO 1350 IZ = 1, NIZS                                                      
        WRITE (NOUT) (SDTXS(IX,IZ),IX=1,NXS),(SDTYS(IY,IZ),IY=-NYS,NYS),        
     >               (SDYXS(IX,IZ),IX=1,NXS),(SDYYS(IY,IZ),IY=-NYS,NYS)         
 1350 CONTINUE                                                                  
C                                                                               
C================= WRITE SCTXS,SCTYS ===================================        
C                                                                               
      DO 1360 IZ = 1, NLS                                                     
        WRITE (NOUT) (SCTXS(IX,IZ),IX=1,NXS),(SCTYS(IY,IZ),IY=-NYS,NYS)
 1360 CONTINUE                                                                  
C                                                                               
C============ WRITE ADDITIONAL QUANTITIES - UPDATES, ETC ===============        
C                                                                               


      DO 1500 IQS = -NQXSO, NQXSI, IBLOCK-100                                   
        IQE = MIN (IQS+IBLOCK-100-1, NQXSI)                                     
        WRITE (NOUT,IOSTAT=IOS) (QXS(IQX),IQX=IQS,IQE)                          
        WRITE (NOUT,IOSTAT=IOS) (QS (IQX),IQX=IQS,IQE)                          
        write (nout,iostat=ios) (svybar(iqx),iqx=iqs,iqe)
        write (nout,iostat=ios) (svyacc(iqx),iqx=iqs,iqe)
 1500 CONTINUE                                                                  
C                                                                               
      DO 1600 IQS = -NQXSO, 1, (IBLOCK-50)/2                                    
        IQE = MIN (IQS+(IBLOCK-50)/2-1, 0)                                      
        WRITE (NOUT,IOSTAT=IOS) ((QEDGES(IQX,J),IQX=IQS,IQE),J=1,2)             
        WRITE (NOUT,IOSTAT=IOS) ((QTANS (IQX,J),IQX=IQS,IQE),J=1,2)             
        WRITE (NOUT,IOSTAT=IOS) ((QDISTS(IQX,J),IQX=IQS,IQE),J=1,2)             
        WRITE (NOUT,IOSTAT=IOS) ((QTEMBS(IQX,J),IQX=IQS,IQE),J=1,2)             
        WRITE(NOUT,IOSTAT=IOS) ((QTEMBSI(IQX,J),IQX=IQS,IQE),J=1,2)             
        WRITE (NOUT,IOSTAT=IOS) ((QRNBS (IQX,J),IQX=IQS,IQE),J=1,2)             
 1600 CONTINUE                                                                 
C                                                                               
      WRITE (NOUT,IOSTAT=IOS) (FACTA(IZ),FACTB(IZ),IZ=-1,NIZS)                  
      WRITE (NOUT,IOSTAT=IOS) TC,SC,TO,SO,TV,SV,GC,RP                           
C                                                                               
      DO 2000 IYB = -NYS, NYS, JBLOCK                                           
        IYE = MIN (IYB+JBLOCK-1, NYS)                                           
        WRITE (NOUT) ((CTEMBS(IX,IY), IX=1,NXS), IY=IYB,IYE)                    
        WRITE (NOUT) ((CTEMBSI(IX,IY), IX=1,NXS), IY=IYB,IYE)                  
        WRITE (NOUT) ((CRNBS (IX,IY), IX=1,NXS), IY=IYB,IYE)                    
 2000 CONTINUE                                                                  


c
c     Write out particle tracks for debugging - if cstept is greater 
c     than zero - then particle tracks were accumulated.   
c
      write (6,*) 'cstept2:',cstept 
      if (cstept.gt.0) then  
         write(nout) (ptracl(ix),ix = 1,cstept)         
         write(6,*) 'ptracl:',(ptracl(ix),ix=1,cstept)
         write(nout) (((ptracs(ix,iy,iz),ix=1,ptracl(iy))
     >                  ,iy=1,cstept),iz=1,2)
         do 2100 iy = 1,cstept
            write(6,*) 'trajectory:',iy,ptracl(iy)
            do 2100 ix = 1,ptracl(iy)
               write(6,*) ix,':',(ptracs(ix,iy,iz),iz=1,2)
 2100    continue     
      endif  

C                                                                               
 9999 RETURN                                                                    
 9002 FORMAT(1X,'DUMP:     NXS   NYS  NQXSO NQXSI NQYS   NTS  NIZS  NLS'        
     >     ,/8X,8I6,/8X,A,/8X,A,/1X,'IMODE',I3,/1X,'ITER',I4,/)                 
      END                                                                       
