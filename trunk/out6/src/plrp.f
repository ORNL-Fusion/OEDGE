c     -*-Fortran-*-
c
      SUBROUTINE PLRP (NIZS,PLAMS,PIND,PIZS,PLRPCNT,CION,RIZB,
     >                 plrpopt)               
c
c     ATTENTION: added rizb to parameter list; Krieger, IPP 95
c
      use mod_params
      use mod_dynam2
      use mod_dynam3
      use mod_cadas
      use mod_pindata
      use mod_cgeom
      use mod_cioniz
      IMPLICIT    NONE
c     include 'params'                                              
c     include 'dynam2'                                              
c     include 'dynam3'                                              
c
c     added cadas to get ADAS usage flag;  Krieger, IPP 95
c
c     include 'cadas'
c
c     added pindata to access neutral density; Krieger, IPP 95
c
c     include 'pindata'
      integer maplist(maxizs)
      real    lambda(maxizs)
      real    rizb
c
c     neutral temperature added (may become a poloidal distribution
c     in future; Krieger, IPP/95)
c
      real    tbgr
c
c     functions to calculate boron IV CX; Krieger, IPP/95
c
      real    b4atripl, b4branch, b4fnete

      INTEGER     NIZS,CION,plrpopt
      REAL        PLAMS(-1:MAXPLRP)                                     
      INTEGER     PIND(-1:MAXIZS+1),PLRPCNT                             
      INTEGER     PIZS(-1:MAXPLRP)                                      
C                                                                       
C  *********************************************************************
C  *                                                                   *
C  *  PLRP:  THIS ROUTINE DERIVED FROM NOTES 81,84. EXTRACTS PARTICULAR*
C  *  LINE RADIATION DATA FROM TABLES BELOW AND APPLIES TO THE         *
C  *  NUMBER DENSITY ION DISTRIBUTION PASSED IN AS SDLIMS ARRAY.       *
C  *                                                                   *
C  *  INPUT:  NIZS     NUMBER OF IONISATION STATES USED                *
C  *          SDLIMS   NUMBER DENSITY ION DISTRIBUTION FROM /DYNAM2/   *
C  *       OR TIZS     IONISATION DENSITY DISTRIBUTION FROM /DYNAM3/   *
C  *          SDLIM3   NUMBER DENSITY ION 3D DISTB'N   FROM /DYNAM2/   *
C  *       OR TIZ3     IONISATION DENSITY 3D DISTB'N   FROM /DYNAM3/   *
C  *                                                                   *
C  *  OUTPUT: PLAMS    WAVELENGTH OF EACH LINE                         *
C  *          PLRPS    PARTICULAR LINE RADIATION PROFILES              *
C  *          PIND     Indices of plrps for each ionization state      *
C  *          PIZS     Ionization state attached to each PLRP          *
C  *                                                                   *
C  *                                      C.M.FARRELL   NOVEMBER 1987  *
C  *                                                                   *
C  *  Modified Feb 18, 1992 to support multiple PLRPS for each         *
C  *  ioniztion state.    David Elder                                  *
C  *                                                                   *
C  *********************************************************************
C                                                                       
      INTEGER        L,IZ,JZ,IT,IK,IR                                   
      REAL           KETAS(MAXNKS),VAL                                  
C     INCLUDE        "CGEOM"                                            
c     include 'cgeom'                                            
C     INCLUDE        "CIONIZ"                                           
c     include 'cioniz'                                           
C                                                                       
C---- PARTICULAR LINE RADIATION DATA                                    
C---- PARAMETERS NL= NUMBER OF LINES OF DATA IN FOLLOWING TABLE         
C----            NT= NUMBER OF TEMPERATURE POINTS DATA GIVEN FOR        
C---- NOTE IBM ONLY ALLOWS 19 CONTINUATION LINES IN ANY ONE STATEMENT.  
C---- EACH LINE HAS AT NO., IZ STATE, LAMBDA, 5 VALUES, FLAG            
C---- THE FLAG INDICATES 0: IGNORE THIS LINE  1: USE THIS LINE          
C                                                                       
C---- *** ONLY ONE PLRP IS ALLOWED PER IONISATION STATE ***             
C         Changed to allow more than 1 line for all but                 
C         neutrals.                                                     
C                                                                       
      INTEGER     NL,NT                                                 
      PARAMETER   (NL=41,NT=5)                                          
      REAL        RD(9*NL),LLAMS(NL),LETAS(NT,NL),LTBS(NT)              
c                 added another set of temperatures for tungsten
c                 Krieger IPP/97
      real        ltbs1(nt)
      INTEGER     LIONS(NL),LIZS(NL),LFLAGS(NL)                         

      DATA        LTBS   / 10.,20.,50.,100.,200. /                      
      DATA        LTBS1  / 1., 2.,  4., 12., 40. /                      
C                                                                       
C         *************** HELIUM **********************                 
C         NOTE 150: THERE IS NO HE0 DATA, USE HE+ DATA                  
C         NOTE 186: NOW THERE IS, LAMBDA FROM NOTE 187.                 
C         NOTE 303: NEW HE0 DATA                                        
C                                                                       
      DATA        (RD(L), L=1,19*9) /                                   
     >     2, 0, 587.5, 0.01, .054, 0.22, 0.46, 0.78, 0                 
     >,    2, 0, 667.8, 24. , 52. , 127., 211., 330., 1                 
     >,    2, 1, 468.6, 7.9 , 23. , 50. , 74. , 90.,  1                 
C                                                                       
C         *************** BERYLLIUM *******************                 
C                                                                       
     >,    4, 0, 440.7,3.55e+2,3.77e+2,3.99e+2,3.99e+2,3.99e+2,1        
C    >,    4, 0, 440.7,359.0,393.0,423.0,437.0,455.0, 0                 
     >,    4, 0, 825.4, 12.8, 15.1, 17.4, 19.3, 21.5, 1                 
     >,    4, 1, 313.1,.0104,.0267,.0507,.0670,.0774, 0                 
     >,    4, 1, 527.1, 16.6, 24.9, 33.0, 37.7, 40.0, 1                 
C    >,    4, 1, 527.1, 10.2, 15.9, 23.3, 28.8, 34.0, 1                 
     >,    4, 1, 436.1,9.16e+0,1.43e+1,2.11e+1,2.59e+1,2.95e+1,1        
C    >,    4, 1, 436.0,  6.2, 11.4, 20.3, 28.5, 37.0, 0                 
     >,    4, 2, 372.0,.0056,0.066, 0.58, 2.04, 6.07, 1                 
     >,    4, 3, 465.9, 6.92, 17.4, 45.7, 83.6, 140.0,1                 
C                                                                       
C         *************** CARBON **********************                 
C                                                                       
     >,    6, 0, 909.5, 7.4 , 16. , 35. , 56. , 85.,  1                 
     >,    6, 1, 657.8, 1.4 , 5.5 , 19. , 37. , 63.,  1                 
     >,    6, 1, 514.5, 2.47, 6.8 , 17.8, 29.1, 41.6, 1                 
     >,    6, 1,90.409, 0.14, 0.31, 0.55, 0.69, 0.77, 0                 
     >,    6, 1,133.58,0.099,0.319,0.743,1.056,1.301, 0                 
     >,    6, 1, 68.72,0.826,1.600,2.541,2.988,3.203, 0                 
     >,    6, 1, 85.84,1.203,3.658,9.350,14.42,19.79, 0                 
     >,    6, 1, 687.2, 0.85, 1.61, 2.54, 3.01, 3.22, 1                 
     >,    6, 2, 569.6, 31. , 90. , 189., 287., 386., 0 /               
       data      (RD(L),L=19*9+1,23*9) /                                
     >     6, 2, 464.7, 0.19, 1.0 , 3.9 , 7.4 , 12.,  1                 
     >,    6, 2, 977.0,.0014, .013, .059, 0.11, 0.16, 0                 
     >,    6, 2, 459.6, .028, 0.14, 0.44, 0.72, 0.97, 1                 
     >,    6, 3, 31.24, .046, 0.26, 0.83, 1.3 , 1.64, 1 /               
C                                                                       
C         *************** OXYGEN **********************                 
C         NOTE 287: THERE IS NO O0 DATA, USE C0 DATA                    
C                                                                       
      DATA        (RD(L), L=23*9+1, 9*NL) /                             
     >     8, 0, 909.5, 7.4 , 16. , 35. , 56. , 85.,  1                 
     >,    8, 1, 374.9, 1.4 , 4.8 , 14. , 25. , 38.,  0                 
     >,    8, 1, 435.1, 1.3 , 4.6 , 15. , 27. , 42.,  1                 
     >,    8, 1, 441.5, 0.88, 3.4 , 11. , 21. , 32.,  0                 
     >,    8, 1, 391.2, 2.1 , 8.2 , 27. , 50. , 78.,  0                 
     >,    8, 1, 397.3, 0.88, 3.7 , 13. , 24. , 37.,  0                 
     >,    8, 2, 376.0, 0.7 , 4.4 , 19. , 39. , 65.,  1                 
     >,    8, 2, 559.2, 0.81, 4.9 , 21. , 44. , 71.,  0                 
     >,    8, 2, 370.3, 3.0 , 15. , 57. , 172., 188., 0                 
     >,    8, 2, 703.3,.0074, .093, 0.56, 1.2 , 1.9,  0                 
     >,    8, 3, 789.4,1.9E-4,.0084,0.12, 0.32, 0.58, 1                 
     >,    8, 3, 233.6, .012, .087, 0.40, 0.79, 1.1,  0                 
     >,    8, 4, 629.7,1.5E-6,3.E-4,.0089,.033, 0.07, 1                 
     >,    8, 4, 760.4,2.1E-6,5.E-4,.017, .066, 0.14, 0                 
     >,    8, 5,1032.0,3.4E-8,3.3E-5,.0028,.015,0.04, 1                 
C                                                                       
C         ************** CHLORINE *********************                 
C                                                                       
     >,   17, 1, 489.7, 39. , 92. , 205., 312., 430., 1                 
C                                                                       
C         ************** CHROMIUM *********************                 
C                                                                       
     >,   24, 0, 425.4, 0.34, 0.62, 1.  , 1.3 , 1.6,  1
C                                                                       
C         ************** TUNGSTEN *********************                 
C                                                                       
     >,   74, 0, 400.9, 1.6, 4.3, 8.6, 19.1, 35.4  ,  1  /              
*         data below used for EPS97 results (ugh) (Te=2,4,10,20,40 eV)
*    >,   74, 0, 400.9, 5.0, 10.0, 10.8, 11.1, 11.9,  1  /              

c
c
c     added to include ADAS emissivity data; Krieger, IPP 95
c
      if (plrpopt.eq.0) then

C                                                                       
C-----------------------------------------------------------------------
C     INITIALISE LOCAL ARRAYS.  DIAGNOSTICS TO CHANNEL 6                
C     SET PHOTON EFFICIENCY TO CONSTANT VALUES WHEN EXTRAPOLATING BEYOND
C     LAST SPECIFIED POINTS.                                            
C-----------------------------------------------------------------------
C                                                                       
      WRITE (6,9000)                                                    
      DO 20 L = 1, NL                                                   
        LIONS(L) = NINT (RD(9*(L-1)+1))                                 
        LIZS(L)  = NINT (RD(9*(L-1)+2))                                 
        LLAMS(L) = RD(9*(L-1)+3)                                        
        DO 10 IT = 1, 5                                                 
          LETAS(IT,L) = RD(9*(L-1)+3+IT)                                
   10   CONTINUE                                                        
        LFLAGS(L) = NINT(RD(9*(L-1)+9))                                 
        IF (CION.EQ.LIONS(L).AND.LFLAGS(L).NE.0)                        
     >    WRITE (6,9001) LIONS(L),LIZS(L),LLAMS(L),(LETAS(IT,L),IT=1,5) 
   20 CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C     FOR EACH IONISATION STATE CHECK EACH ROW OF DATA TABLE UNTIL      
C     WE FIND A MATCH.  CALCULATE FITTING DETAILS ALONG REFERENCE LINE  
C     USING FITTER ROUTINE.                                             
C     FOR EACH POINT (X,Y), FIND APPROPRIATE ETA AND CALC PLRPS         
C     FOR PRIMARY NEUTRALS (-1), ENSURE THE (0) LINE DETAILS ARE USED.  
C-----------------------------------------------------------------------
C                                                                       
      PLRPCNT = -1                                                      
      DO 500 IZ = -1, NIZS                                              
       PIND(IZ) = PLRPCNT                                               
       JZ = MAX (IZ, 0)                                                 
       DO 400 L = 1, NL                                                 
        IF (CION.EQ.LIONS(L).AND.JZ.EQ.LIZS(L).AND.LFLAGS(L).NE.0)THEN  
         PLAMS(PLRPCNT) = LLAMS(L)                                      
         PIZS(PLRPCNT) = IZ                                             
         DO 310 IR = 1, NRS                                             
          if (cion.eq.74) then
*           CALL FITTER (5,LTBS1,LETAS(1,L),NKS(IR),                       
*    >                   KTEBS(1,IR),KETAS,'LINEAR')                      
*           IPP Krieger/98  Fit to best available experimental data
            do ik=1,nks(ir)
              ketas(ik)=-0.756+4.3069*alog(1.0+ktebs(1,ir))**1.6224
            end do
          else
            CALL FITTER (5,LTBS,LETAS(1,L),NKS(IR),                       
     >                   KTEBS(1,IR),KETAS,'LINEAR')                      
          endif
          DO 300 IK = 1, NKS(IR)                                        
           VAL = SDLIMS(IK,IR,IZ) / KFIZS(IK,IR,JZ) / KETAS(IK)   
c           VAL = SNGL(DDLIMS(IK,IR,IZ)) / KFIZS(IK,IR,JZ) / KETAS(IK)   
C===>      VAL = TIZS(IK,IR,IZ) / ETA                                   
           PLRPS(IK,IR,PLRPCNT) = VAL                                   
  300     CONTINUE                                                      
  310    CONTINUE                                                       
         PLRPCNT = PLRPCNT + 1                                          
         IF (PLRPCNT.GT.MAXPLRP) THEN                                   
           WRITE(6,*) 'NUMBER OF PLRPS EXCEEDS MAXIMUM: PLRPS INVALID'  
           CALL PRC('****     WARNING     *************************')   
           CALL PRC('NUMBER OF PLRPS EXCEEDS MAXIMUM: PLRPS INVALID')   
           CALL PRC('****     WARNING     *************************')   
           GOTO 600                                                     
         ENDIF                                                          
        ENDIF                                                           
  400  CONTINUE                                                         
  500 CONTINUE                                                          
C                                                                       
  600 CONTINUE                                                          
      PIND(NIZS+1) = PLRPCNT                                            
      PLRPCNT = PLRPCNT -1                                              
                                                                        
      WRITE(6,*) 'PLRPCNT:',PLRPCNT                                     
      WRITE(6,*) 'PIZS:',(PIZS(IZ),PLAMS(IZ),IZ=-1,PLRPCNT)             
      WRITE(6,*) 'PIND:',(PIND(IZ),IZ=-1,NIZS+1)                        
c
c     ADAS PLS based PLRP's
c                                                                        
      elseif (plrpopt.eq.1) then 
c
c     added to include ADAS emissivity data; Krieger, IPP 95
c
c     first, we go through the ADAS pls file and check the entries
c     for the different ionisation stages. Normally, ADAS allows for
c     only one specific line per ionisation state. To extend that,
c     the header line of each data block needs to be reorganized
c     a bit. It needs to be in the format given below:
c-------------------/z1 =  1  908.9nm CI... (other comments)
c     This enables us to map the Z data blocks of a pls file
c     to more than one spectral line of a given ionization state

c     IPP/09 Krieger - fix to use specified year for impurities
c     year='95'
      write(year,'(i2.2)') iyearz
      yeardf='89'
      call ADASCHK(YEAR,YEARDF,cion,maplist,lambda)     
      write(6,*) 'PLRP LINES IN ADAS PLS FILE:'
      write(6,*) 'Z: ',(maplist(iz)-1,iz=1,cion)
      write(6,*) 'LAMBDA (nm): ',(lambda(iz),iz=1,cion)

C-----------------------------------------------------------------------
C     FOR EACH IONISATION STATE CHECK IF THERE ARE MATCHING LINES
C     FOR EACH POINT (X,Y), CALC PLRPS USING ADAS
C     FOR PRIMARY NEUTRALS (-1), ENSURE THE (0) LINE DETAILS ARE USED.  
C-----------------------------------------------------------------------
C                                                                       
      PLRPCNT = -1                                                      
      DO 501 IZ = -1, NIZS                                              
       PIND(IZ) = PLRPCNT                                               
       JZ = MAX (IZ, 0)                                                 
       DO 401 L = 1, CION
        IF (JZ.EQ.(maplist(l)-1))THEN  
         PLAMS(PLRPCNT) = lambda(L)                                      
         PIZS(PLRPCNT) = IZ                                             
*	 write(6,*) 'ik ir  ne        n0        nz        '//
*    >              'nz+1      vale      val0'
         DO 311 IR = 1, NRS                                             
          DO 710 IK = 1, NKS(IR)                                            
            PNESA(IK) = KNBS(IK,IR) * RIZB                                  
            PTESA(IK) = KTEBS(IK,IR)                                        
c           added for CX radiation; Krieger, IPP/95
            pnhs(ik)  = pinatom(ik,ir)
            tbgr = 5.0
  710     CONTINUE                                                          
          ICLASS = 7                                                        
          CALL ADASRD(YEAR,CION,L,ICLASS,NKS(IR),PTESA,PNESA,   
     >              PCOEF(1,L))                                      
c
c          CALL ADASRD(YEAR,YEARDF,CION,L,ICLASS,NKS(IR),PTESA,PNESA,   
c     >                PCOEF(1,L))                                      
c
          DO 301 IK = 1, NKS(IR)                                        
           VAL = SDLIMS(IK,IR,IZ) * PNESA(IK) * PCOEF(IK,L)
c           VAL = SNGL(DDLIMS(IK,IR,IZ)) * PNESA(IK) * PCOEF(IK,L)
           PLRPS(IK,IR,PLRPCNT) = VAL
c          boron IV CX radiation; Krieger, IPP/95
           if (cion.eq.5.and.jz.eq.3) then
             VAL = SDLIMS(IK,IR,IZ+1) * pnhs(ik) * 
     >             b4atripl(tbgr,sdts(ik,ir,iz+1)) *
     >            (b4branch(tbgr,sdts(ik,ir,iz+1)) + 
     >             b4fnete(pnesa(ik),ptesa(ik)))
c             VAL = SNGL(DDLIMS(IK,IR,IZ+1)) * pnhs(ik) * 
c     >             b4atripl(tbgr,sngl(ddts(ik,ir,iz+1))) *
c     >            (b4branch(tbgr,sngl(ddts(ik,ir,iz+1))) + 
c     >             b4fnete(pnesa(ik),ptesa(ik)))
*            if (mod(ik,8).eq.2.and.mod(ir,4).eq.0) then
*              write(6,'(2(1x,i2),1p,6(2x,e8.2))') 
*    >              ik,ir,
*    >              pnesa(ik),pnhs(ik),SDLIMS(IK,IR,IZ),
*    >              SDLIMS(IK,IR,IZ+1),plrps(ik,ir,plrpcnt),
*    >              val
*            endif
             PLRPS(IK,IR,PLRPCNT) = PLRPS(IK,IR,PLRPCNT) + val
           endif
  301     CONTINUE                                                      
  311    CONTINUE                                                       
         PLRPCNT = PLRPCNT + 1                                          
         IF (PLRPCNT.GT.MAXPLRP) THEN                                   
           WRITE(6,*) 'NUMBER OF PLRPS EXCEEDS MAXIMUM: PLRPS INVALID'  
           CALL PRC('****     WARNING     *************************')   
           CALL PRC('NUMBER OF PLRPS EXCEEDS MAXIMUM: PLRPS INVALID')   
           CALL PRC('****     WARNING     *************************')   
           GOTO 601                                                     
         ENDIF                                                          
        ENDIF                                                           
  401  CONTINUE                                                         
  501 CONTINUE                                                          
C                                                                       
  601 CONTINUE                                                          
      PIND(NIZS+1) = PLRPCNT                                            
      PLRPCNT = PLRPCNT -1                                              
                                                                        
      WRITE(6,*) 'PLRPCNT:',PLRPCNT                                     
      WRITE(6,*) 'PIZS:',(PIZS(IZ),PLAMS(IZ),IZ=-1,PLRPCNT)             
      WRITE(6,*) 'PIND:',(PIND(IZ),IZ=-1,NIZS+1)                        
      endif
                                                                        
      RETURN                                                            
C                                                                       
 9000 FORMAT(/1X,'PLRP:  PHOTON EFFICIENCY DATA',//1X,                  
     >  '      SPECIES          LINE                  ELECTRON TEMPER', 
     >  'ATURE  (EV)',/1X,                                              
     >  '    ATOM   CHARGE    WAVELENGTH      10        20        50 ', 
     >  '      100       200',/1X,81('-'))                              
 9001 FORMAT(1X,2I8,3X,F9.1,6X,1P,5G10.3)                               
      END                                                               



      SUBROUTINE ADASCHK(YEAR,YEARDF,IZ0,maplist,lambda)     
C                                                                       
C  READ THE REQUESTED RATE COEFFICIENT FROM THE ADAS MASTER ELEMENT     
C  FILES:                                                               
C        ICLASS = 7: PHOTON EMISSIVITY FOR DIAGNOSTIC LINES 
C                
C  THIS ROUTINE USES THE STANDARD ADAS EXTRACTION ROUTINE D2DATA AND    
C  REALLY ONLY PROVIDES A 'CLEAN' INTERFACE, TAKING CARE OF CHANGES     
C  IN UNITS AND IN PRECISION OF VARIABLES.  IF THE REQUESTED DATA       
C  DOESN'T EXIST (IFAIL=1 RETURNED FROM D2DATA) THE PROGRAM IS STOPPED. 
C                                                                       
      use mod_params
      use mod_cadas2
      IMPLICIT NONE                                                     
C     INCLUDE   "PARAMS"                                                
c     include    'params'                                                
C     INCLUDE   "CADAS2"                                                
c     include    'cadas2'                                                
C                                                                       
      CHARACTER*2 YEAR, YEARDF                                          
      INTEGER IZ0, IZ1, ICLASS, NPTS                                    
      parameter (npts=2)
      REAL TE(NPTS), NE(NPTS)
C                                                                       
      INTEGER I, J                                                      
      integer maplist(maxizs)
      real lambda(maxizs)
      logical isuccess
C                                                                       
      IEVCUT = 0                                                        
      ICLASS=7
C                                                                       
      data ne /1.e19, 1.e19/
      data te /10.,100./
      DO J = 1, NPTS                                                    
        DTEV(J) = DBLE(ALOG10(TE(J)))                                   
        DDENS(J) = DBLE(ALOG10(NE(J)*1.0E-6))                           
      ENDDO                                                             
      do i=0,iz0-1 
        iz1=i+1
        CALL D2INQU(YEAR, YEARDF, TITLF, IFAIL,                           
     >              IZ0, IZ1, ICLASS, NPTS, IEVCUT,                       
     >              MAXADS, ITMAXD, IDMAXD, IZMAXD,                       
     >              DTEV, DDENS,                                          
     >              DTEVD, DDENSD, DRCOFD, ZDATA,                         
     >              DRCOFI)                                               
        isuccess=.false.
        read(titlf,'(26x,i2,1x,f6.1)',err=10) maplist(iz1), lambda(iz1)
        isuccess=.true.
   10   continue
        if (.not.isuccess) then
          maplist(iz1)=-1
          lambda(iz1)=0.
        endif
        IF (IFAIL.EQ.1) THEN                                              
          WRITE(6,1000) IZ0, IZ1,YEAR, YEARDF                             
          return                                                           
        ENDIF                                                             
      enddo
C                                                                       
 1000 FORMAT(' ERROR READING REQUESTED ATOMIC DATA!',/,                 
     >       ' MASTER ELEMENT FILE FOR NUCLEAR CHARGE ',I2,             
     >       ' AND ION CHARGE ',I2,                                     
     >       ' WAS NOT FOUND IN YEAR ',A2,' OR ',A2)                    
C                                                                       
      RETURN                                                            
      END                                                               
c
c
c                                                                        
       SUBROUTINE D2INQU( YEAR   , YEARDF , TITLF  , IFAIL              
     &                  , IZ0    , IZ1    , ICLASS , ITMAX  , IEVCUT    
     &                  , ITDIMD , ITMAXD , IDMAXD , IZMAXD             
     &                  , DTEV   , DDENS                                
     &                  , DTEVD  , DDENSD , DRCOFD , ZDATA              
     &                  , DRCOFI                                        
     &                  )                                               
       IMPLICIT REAL*8(A-H,O-Z)                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C PURPOSE : TO EXTRACT 'SANC0' COLLISIONAL DIELECTRONIC DATA            
C                                                                       
C NOTE    : THE SOURCE DATA IS CONTAINED AS SEQUENTIAL DATASETS         
C           AS FOLLOWS:                                                 
C                                                                       
C                   (1) JETUID.ACD<YR>#<IEL).DATA                       
C                   (2) JETUID.SCD<YR>#<IEL>.DATA                       
C                   (3) JETUID.CCD<YR>#<IEL>.DATA                       
C                   (4) JETUID.PRB<YR>#<IEL>.EV<CUT>.DATA               
C                   (5) JETUID.PLT<YR>#<IEL>.EV<CUT>.DATA               
C                   (6) JETUID.PRC<YR>#<IEL>.EV<CUT>.DATA               
C                   (7) JETUID.PLS<YR>#<IEL>.DATA                       
C                                                                       
C       WHERE, <YR>  = TWO INTEGERS  FOR THE YEAR SELECTED              
C              <IEL> = ELEMENT NAME                                     
C              <CUT> = ENERGY CUT-OFF (EV)                              
C                                                                       
C              IF <CUT> = 0 THEN .EV<CUT> IS DELETED FROM ABOVE FILES.  
C                                                                       
C INPUT  : (C*2)  YEAR      = YEAR OF DATA                              
C          (C*2)  YEARDF    = DEFAULT YEAR OF DATA IF REQUESTED YEAR    
C                             DOES NOT EXIST.                           
C          (I*4)  IZ0       = NUCLEAR CHARGE                            
C          (I*4)  IZ1       = MINIMUM ION CHARGE + 1                    
C          (I*4)  ICLASS    = CLASS OF DATA (1 - 6)                     
C          (I*4)  ITMAX     = NUMBER OF ( DTEV() , DDENS() ) PAIRS      
C          (I*4)  IEVCUT    = ENERGY CUT-OFF (EV)                       
C          (R*8)  DTEV()    = DLOG10(ELECTRON TEMPERATURES (EV))        
C          (R*8)  DDENS()   = DLOG10(ELECTRON DENSITIES (CM-3))         
C                                                                       
C OUTPUT : (C*80) TITLF     = INFORMATION STRING                        
C          (I*4)  ITDIMD    = MAXIMUM NUMBER OF DATA TEMP & DENS        
C          (I*4)  ITMAXD    = NUMBER OF DATA DTEVD()                    
C          (I*4)  IDMAXD    = NUMBER OF DATA DDENS()                    
C          (I*4)  IZMAXD    = NUMBER OF DATA ZDATA()                    
C          (I*4)  ITDIMD    = MAXIMUM NUMBER OF DATA TEMP & DENS        
C          (I*4)  ZDATA()   = Z1 CHARGES IN DATASET                     
C          (I*4)  IFAIL     = -1   IF ROUTINE SUCCESSFUL BUT THE DEFAULT
C                                  YEAR FOR THE DATA WAS USED.          
C                           = 0    IF ROUTINE SUCCESSFUL - DATA FOT THE 
C                                  REQUESTED YEAR USED.                 
C                           = 1    IF ROUTINE OPEN STATEMENT FAILED     
C          (R*8)  DTEVD()   = DLOG10(DATA ELECTRON TEMPERATURES (EV))   
C          (R*8)  DDENSD()  = DLOG10(DATA ELECTRON DENSITIES (CM-3))    
C          (R*8)  DRCOFD()  = DLOG10(DATA RATE COEFFICIENTS (CM-3/S))   
C          (R*8)  DRCOFI()  = INTERPOLATION OF DRCOFD(,,) FOR           
C                             DTEV() & DDENS()                          
C                                                                       
C PROGRAM: (C*2)  SEQUA()   = ION NAMES FOR A PARTICULAR IZ0            
C          (C*2)  ESYM      = ELEMENT SYMBOL FOR NUCLEAR CHARGE IZ0
C          (C*60) USERID    = USER ID UNDER WHICH ADAS DATA IS STORED   
C          (C*60) DSNAME    = FILE NAME ( SEE ABOVE TYPES )             
C          (C*80) STRING    = GENERAL VARIABLE                          
C          (C*80) BLANK     = BLANK STRING                              
C          (C*2)  YEARSV    = LAST YEAR USED IN THIS ROUTINE            
C          (I*4)  IREAD     = INPUT STREAM FOR OPEN STATEMENT           
C          (I*4)  IZ0SV     = LAST IZ0 USED IN THIS ROUTINE             
C          (I*4)  ICLSV     = LAST ICLASS USED IN THIS ROUTINE          
C          (I*4)  INDXZ1    = LOCATION OF IZ1 IN ZDATA()                
C          (I*4)  LCK       = MUST BE GREATER THAN 'ITMAXD' & 'IDMAXD'  
C                             & 'ITMAX' - ARRAY SIZE FOR SPLINE CALCS.  
C          (R*8)  A()       = GENERAL ARRAY                             
C          (R*8)  DRCOF0(,) = INTERPOLATION OF DRCOFD(,,) W.R.T DTEV()  
C          (L*8)  LEXIST    = TRUE --- FILE TO OPEN EXISTS ELSE NOT     
C                                                                       
C PE BRIDEN = ADDED VARIABLES (14/01/91)                                
C                                                                       
C          (I*4)  L1      = PARAMETER = 1                               
C          (I*4)  IOPT    = DEFINES THE BOUNDARY DERIVATIVES FOR THE    
C                             SPLINE ROUTINE 'XXSPLN', SEE 'XXSPLN'.    
C                                                                       
C          (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
C                                      TO X-AXIS.                       
C                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS  
C                                      RELATING TO X-AXIS.              
C                                      (I.E. THEY WERE SET IN A PREVIOUS
C                                            CALL )                     
C                           (VALUE SET TO .FALSE. BY 'XXSPLN')          
C                                                                       
C                                                                       
C          (R*8)  DY()    = SPLINE INTERPOLATED DERIVATIVES             
C                                                                       
C          (R*8 ADAS FUNCTION - 'R8FUN1' ( X -> X) )                    
C                                                                       
C PE BRIDEN = ADDED VARIABLES (23/04/93)                                
C                                                                       
C          (I*4 ADAS FUNCTION - 'I4UNIT' (OUTPUT STREAM))               
C                                                                       
C AUTHOR : JAMES SPENCE (TESSELLA SUPPORT SERVICES PLC)                 
C          K1/0/80                                                      
C          JET  EXT. 4866                                               
C                                                                       
C DATE   : 22/02/90                                                     
C                                                                       
C DATE   : 21/08/90 PE BRIDEN - REVISION: SEQUA(43) CHANGED ('TE'->'TC')
C                                                                       
C DATE   : 08/10/90 PE BRIDEN - REVISION: RENAMED SUBROUTINE            
C                                                                       
C DATE   : 12/11/90 PE BRIDEN - CORRECTION: MOVE THE SETTING OF 'INDXZ1'
C                                           TO AFTER THE  '20 CONTINUE' 
C                                           STATEMENT.   ALSO SAVE  THE 
C                                           VALUE OF 'IZ1MIN'.          
C                                                                       
C DATE   : 14/01/91 PE BRIDEN - ADAS91:     CALLS TO NAG SPLINE ROUTINES
C                                           'E01BAF' & 'E02BBF' REPLACED
C                                           BY  CALLS   TO  ADAS  SPLINE
C                                           ROUTINE 'XXSPLN'.           
C                                                                       
C DATE   : 25/06/91 PE BRIDEN - CORRECTION: CHANGED FOLLOWING DIMENSION:
C                                            'DIMENSION DRCOFI(ITDIMD)' 
C                                           TO                          
C                                            'DIMENSION DRCOFI(ITMAX)'  
C                                                                       
C DATE   : 07/08/91 PE BRIDEN - ADDED ERROR HANDLING IF THE OPEN STATE- 
C                               MENT FAILS. (IFAIL=1 RETURNED)          
C                                                                       
C DATE   : 27/04/92 PE BRIDEN - ADDED DEFAULT YEAR FOR DATA IF REQUESTED
C                               YEAR DOES NOT EXIST. (ADDED 'YEARDF')   
C                               INTRODUCED IFAIL = -1 IF DEFAULT YEAR   
C                               WAS USED AND NOT THE REQUESTED YEAR.    
C                                                                       
C DATE   : 10/03/93 PE BRIDEN - ALLOWED INPUT DATA SETS TO BE ACCESSED  
C                               FROM ANY USERID (DEFAULT = JETSHP)      
C                               - INTRODUCED USERID VARIABLE AND CALL   
C                                 TO XXUID.                             
C                                                                       
C DATE   : 23/04/93 PE BRIDEN - ADDED I4UNIT FUNCTION TO WRITE          
C                               STATEMENTS FOR SCREEN MESSAGES          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
       integer itmax,itdimd,iread,lck,iz0sv,iclsv,ievsv,
     >         ifail,iz0,iz1,iclass,ievcut,itmaxd,idmaxd,izmaxd,
     >         iz1min,izmax,iz1max,id,it,iz,indxz1


       INTEGER   L1, MCLASS                                                     
C                                                                       
       PARAMETER ( L1    =  1 )                                         
       PARAMETER ( IREAD = 12 , LCK = 100 )
       PARAMETER ( MCLASS = 7 )                             
C                                                                       
       INTEGER   LENF1, LENF2, LENF3, LENF4, LENF5, LENF6
c
       INTEGER   I4UNIT                                                 
       INTEGER   IOPT                                                   
       integer   lenu,lenstr                                            
       external  lenstr                                                 
C                                                                       
       real*8 DTEV(ITMAX)   , DDENS(ITMAX)                           
       REAL*8 DTEVD(ITDIMD) , DDENSD(ITDIMD) , ZDATA(ITDIMD)         
       REAL*8 DRCOFD(ITDIMD,ITDIMD,ITDIMD)                           
       REAL*8 DRCOFI(ITMAX)                                          
C                                                                       
       dimension SEQUA(50)                                              
C                                                                       
       REAL*8 A(LCK)                                                 
       REAL*8    DY(LCK)                                                
       REAL*8 DRCOF0(LCK,LCK)                                        
C                                                                       
       CHARACTER DEFADF*5, ESYM*2, XFESYM*2, CLASS(MCLASS)*3
c
       CHARACTER YEAR*2  , YEARDF*2  , YEARSV*2  , TITLF*80  , EVCUT*6  
c
       CHARACTER USERID*60,SEQUA*2 , DSNAME*60 , STRING*80 , BLANKS*80  
       character string1*80
       dimension string1(50)
C                                                                       
       LOGICAL   LEXIST  , LSETX                                        
C                                                                       
       EXTERNAL  R8FUN1                                                 
C                                                                       
C------SET DEFAULT DIRECTORY-----------------
C 
       PARAMETER (DEFADF='adf11')            
c
       SAVE      IZ1MIN                                                 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C       DATA  SEQUA/'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',  
C     &             'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',  
C     &             'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',  
C     &             'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',  
C     &             'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN'/  
       DATA  SEQUA/'h ','he','li','be','b ','c ','n ','o ','f ','ne',   
     &             'na','mg','al','si','p ','s ','cl','ar','k ','ca',   
     &             'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',   
     &             'ga','ge','as','se','br','kr','rb','sr','y ','zr',   
     &             'nb','mo','tc','ru','rh','pd','ag','cd','in','sn'/   
       DATA BLANKS/'                                                    
     &                           '/                                     
C                                                                       
       DATA YEARSV/'  '/                                                
       DATA IZ0SV /0   /                                                
       DATA ICLSV /0   /                                                
       DATA IEVSV /0   /                                                
       DATA CLASS /'acd', 'scd', 'ccd', 'prb', 'plt', 'prc', 'pls'/
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
       IF(YEAR.EQ.YEARSV.AND.IZ0.EQ.IZ0SV.AND.ICLASS.EQ.ICLSV.AND.      
     &    IEVCUT.EQ.IEVSV) GOTO 20                                      
       YEARSV = YEAR                                                    
       IZ0SV  = IZ0                                                     
       ICLSV  = ICLASS                                                  
       IEVSV  = IEVCUT                                                  
       IFAIL  = 0                                                       
C                                                                       
C------ENERGY CUTOFF----------------------------------------------------
C                                                                       
c      IF( IEVCUT.GT.0 ) THEN                                            
c          WRITE(EVCUT,1050) IEVCUT                                      
c          DO 50 I = 1 , 6                                               
c             IF( EVCUT(I:I).EQ.' ' ) LEVMIN = I                         
c   50     CONTINUE                                                      
c          LEVMIN = LEVMIN + 1                                           
c      END IF                                                            
c
C------ENERGY CUTOFF----------------------------------------------------
C                                                                       
      IF( IEVCUT.GT.0 ) THEN                                            
          WRITE(EVCUT,1050) IEVCUT                                      
          CALL XXSLEN(EVCUT,LENF5,LENF6)
      END IF                                                            

C                                                                       
C------GET ADAS DATA SOURCE USERID--------------------------------------
C                                                                       
c
c      Set userid to default userid set in the environment 
c
       USERID = '*'                                                     
       CALL XXUID(USERID)                                               
       CALL XXSLEN(USERID,LENF1,LENF2) 
C                                                                       
C------ELEMENT NAME-----------------------------------------------------
C                                                                       
       ESYM = XFESYM(IZ0)                                 
       CALL XXSLEN(ESYM,LENF3,LENF4)
C                                                                       
C------FILE NAME--------------------------------------------------------
C                                                                       
   10  DSNAME=USERID(LENF1:LENF2)//'/'//DEFADF//'/'//CLASS(ICLASS)//   
     &        YEARSV//'/'//CLASS(ICLASS)//YEARSV//'_'//
     &        ESYM(LENF3:LENF4)//'.dat'                                           
       IF ((ICLASS.GE.4.AND.ICLASS.LE.6) .AND. IEVCUT.NE.0) THEN                                        
         DSNAME=USERID(LENF1:LENF2)//'/'//DEFADF//'/'//CLASS(ICLASS)//   
     &          YEARSV//'/'//CLASS(ICLASS)//YEARSV//'_'//
     &          ESYM(LENF3:LENF4)//'_ev'//EVCUT(LENF5:LENF6)//'.dat'
       ENDIF                                      
c
       write(6,*) 'ADAS PLRP:',dsname 
c
C                                                                       
C------GET ADAS DATA SOURCE USERID--------------------------------------
C                                                                       
c       USERID = '?'                                                     
c       CALL XXUID(USERID)                                               
c       lenu = lenstr(userid)                                            
C                                                                       
C------FILE NAME--------------------------------------------------------
C                                                                       
c   10  IF (SEQUA(IZ0)(2:2).EQ.' ') THEN                                 
c            LEN=1                                                       
c       ELSE                                                             
c            LEN=2                                                       
c       ENDIF                                                            
C                                                                       
c       if     (iclass.eq.1) then                                        
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'acd'                                            
c       elseif (iclass.eq.2) then                                        
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'scd'                                            
c       elseif (iclass.eq.3) then                                        
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'ccd'                                            
c       elseif (iclass.eq.4) then                                        
c         if( ievcut.eq.0 ) then                                         
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'prb'                                            
c         else                                                           
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'prb'//'ev'//evcut(levmin:6)                     
c         end if                                                         
c       elseif (iclass.eq.5) then                                        
c         if( ievcut.eq.0 ) then                                         
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'plt'                                            
c         else                                                           
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'plt'//'ev'//evcut(levmin:6)                     
c         end if                                                         
c       elseif (iclass.eq.6) then                                        
c         if( ievcut.eq.0 ) then                                         
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'prc'                                            
c         else                                                           
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'prc'//'ev'//evcut(levmin:6)                     
c         end if                                                         
c       elseif (iclass.eq.7) then                                        
c           dsname=userid(1:lenu)//'/'//year//'/'//sequa(iz0)(1:len)//   
c     &            '/'//'pls'                                            
c       endif                                                            
c
c       write (6,*) 'dsname:',dsname                                     
C                                                                       
C------PE BRIDEN - MODIFICATION 27/04/92 - INCLUSION OF DEFAULT YEAR -  



C                                                                       
C------DOES FILE TO BE OPEN EXIST OR NOT--------------------------------
C                                                                       
       INQUIRE(FILE=DSNAME,EXIST=LEXIST)                                
C                                                                       
          IF ( (.NOT.LEXIST) .AND. (YEAR.NE.YEARDF) ) THEN              
             WRITE(I4UNIT(0),1060) DSNAME , YEARDF                      
             DSNAME(33:34) = YEARDF                                     
             INQUIRE(FILE=DSNAME,EXIST=LEXIST)                          
             IFAIL  = -1                                                
             YEARSV = YEARDF                                            
          ENDIF                                                         
C                                                                       
       IF( .NOT.LEXIST ) GOTO 9999                                      
C                                                                       
C------PE BRIDEN - END OF MODIFICATION 27/04/92                         
C                                                                       
C------READ FILE # IREAD------------------------------------------------
C                                                                       
C       OPEN(UNIT=IREAD,FILE=DSNAME,ACTION='READ',ERR=9999)             
       OPEN(UNIT=IREAD,FILE=DSNAME,ERR=9999)                            
C                                                                       
       READ(IREAD,1010) IZMAX , IDMAXD , ITMAXD , IZ1MIN , IZ1MAX       
C                                                                       
       READ(IREAD,1020) STRING                                          
       READ(IREAD,1040) ( DDENSD(ID) , ID = 1 , IDMAXD )                
       READ(IREAD,1040) ( DTEVD(IT)  , IT = 1 , ITMAXD )                
C                                                                       
       IZMAXD = 0                                                       
       DO 150 IZ = IZ1MIN , IZ1MAX                                      
          IZMAXD = IZMAXD + 1                                           
          ZDATA(IZMAXD) = IZ                                            
C         IF( IZ .EQ. IZ1 ) INDXZ1 = IZ                                 
          READ(IREAD,1020)STRING                                        
          string1(izmaxd)=string
          DO 100 IT = 1 , ITMAXD                                        
             READ(IREAD,1040) ( DRCOFD(IZMAXD,IT,ID) , ID = 1 , IDMAXD )
  100     CONTINUE                                                      
  150  CONTINUE                                                         
C                                                                       
       CLOSE(12)                                                        
C                                                                       
C------INTERPOLATE USING SPLINES (NAG ALGORITHM)------------------------
C                                                                       
   20  CONTINUE                                                         
C                                                                       
C------PE BRIDEN - CORRECTION 12/11/90 - SET INDXZ1 AFTER '20 CONTINUE'-
C                                                                       
       INDXZ1 = IZ1 - IZ1MIN + 1                                        
       titlf=string1(indxz1)
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C                                                                       
       RETURN
C                                                                       
C-----------------------------------------------------------------------
C DATA SET OPENING/EXISTENCE ERROR HANDLING                             
C-----------------------------------------------------------------------
C                                                                       
 9999  IFAIL  = 1                                                       
       YEARSV = '  '                                                    
       IZ0SV  = 0                                                       
       ICLSV  = 0                                                       
       RETURN                                                           
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
 1000  FORMAT('FILE = ',1A60)                                           
 1010  FORMAT(5I5)                                                      
 1020  FORMAT(1A80)                                                     
 1040  FORMAT(8F10.5)                                                   
 1050  FORMAT(I6)                                                       
 1060  FORMAT(1X,'NOTE: REQUESTED DATASET - ',A60/
     1        7X,'DOES NOT EXIST. USING DEFAULT YEAR (',
     1        A2,') DATASET INSTEAD'/)  
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
       END                                                              



c     added routines for treatment of BoronV charge exchange radiation
c     Krieger, IPP/95

      real function b4fnete(nearg, tearg)
      implicit none
c     ne = electron density (m^-3)
c     te = electron temperature (eV)
c     b4fnete = photons per electron (see below)

      integer dimne, dimte
      parameter (dimne=5)
      parameter (dimte=21)

      real ne_t(dimne), te_t(dimte), fnete_t(dimne,dimte)
      save ne_t, te_t, fnete_t

      real nearg, tearg, ne, te, fnetebuf

      data ne_t /1.0e+17, 1.0e+18, 1.0e+19, 1.0e+20, 1.0e+21/
      data te_t /3.0e+00, 4.0e+00, 5.0e+00, 6.0e+00, 7.0e+00, 8.0e+00,
     >           9.0e+00, 1.0e+01, 1.2e+01, 1.5e+01, 1.8e+01, 2.0e+01,
     >           2.5e+01, 3.0e+01, 3.5e+01, 4.0e+01, 4.5e+01, 5.0e+01,
     >           7.5e+01, 1.0e+02, 2.0e+02/

c     fnete_t = number of emitted photons per electron captured in
c               triplett state of BIV at 282 nm,
c               without primary electron from initial cascade transition
      data fnete_t
     >     /1.4462e+01, 1.2106e+01, 4.6039e+00, 6.3970e-01, 6.6562e-02,
     >      1.5449e+01, 1.2955e+01, 4.9555e+00, 6.9069e-01, 7.1900e-02,
     >      1.6006e+01, 1.3482e+01, 5.2323e+00, 7.3497e-01, 7.6597e-02,
     >      1.6408e+01, 1.3891e+01, 5.4824e+00, 7.7728e-01, 8.1116e-02,
     >      1.6743e+01, 1.4247e+01, 5.7205e+00, 8.1900e-01, 8.5595e-02,
     >      1.7043e+01, 1.4574e+01, 5.9518e+00, 8.6058e-01, 9.0074e-02,
     >      1.7322e+01, 1.4882e+01, 6.1784e+00, 9.0218e-01, 9.4570e-02,
     >      1.7586e+01, 1.5175e+01, 6.4010e+00, 9.4383e-01, 9.9085e-02,
     >      1.8082e+01, 1.5730e+01, 6.8361e+00, 1.0274e+00, 1.0818e-01,
     >      1.8763e+01, 1.6494e+01, 7.4653e+00, 1.1531e+00, 1.2196e-01,
     >      1.9382e+01, 1.7191e+01, 8.0684e+00, 1.2794e+00, 1.3590e-01,
     >      1.9764e+01, 1.7622e+01, 8.4568e+00, 1.3638e+00, 1.4528e-01,
     >      2.0627e+01, 1.8601e+01, 9.3835e+00, 1.5756e+00, 1.6904e-01,
     >      2.1380e+01, 1.9460e+01, 1.0252e+01, 1.7886e+00, 1.9325e-01,
     >      2.2041e+01, 2.0219e+01, 1.1068e+01, 2.0028e+00, 2.1793e-01,
     >      2.2628e+01, 2.0896e+01, 1.1836e+01, 2.2182e+00, 2.4307e-01,
     >      2.3151e+01, 2.1503e+01, 1.2560e+01, 2.4347e+00, 2.6868e-01,
     >      2.3621e+01, 2.2050e+01, 1.3244e+01, 2.6521e+00, 2.9476e-01,
     >      2.5398e+01, 2.4142e+01, 1.6154e+01, 3.7489e+00, 4.3194e-01,
     >      2.6572e+01, 2.5541e+01, 1.8401e+01, 4.8479e+00, 5.7952e-01,
     >      2.8851e+01, 2.8285e+01, 2.3647e+01, 8.9586e+00, 1.2422e+00/

      ne = nearg * 1.e-6				! m^-3->cm^-3
      te = tearg
      ne = min(max(ne, ne_t(1)), ne_t(dimne))
      te = min(max(te, te_t(1)), te_t(dimte))

      call linin2(ne_t,te_t,fnete_t,dimne,dimte,ne,te,fnetebuf)
      b4fnete=fnetebuf

      return

      end

c     ------------------------------------------------------------------

      real function b4branch(t0arg, tzarg)
      implicit none
c     t0 = neutrals temperature (eV)
c     tz = impurity temperature (eV)
c     b4branch = branching ratio (see below)

      integer dimt0, dimtz
      parameter (dimt0=7)
      parameter (dimtz=22)

      real t0_t(dimt0), tz_t(dimtz), branch_t(dimt0,dimtz)
      save t0_t, tz_t, branch_t

      real t0arg, tzarg, t0, tz, branchbuf

      data t0_t /1.0e+00, 2.0e+00, 5.0e+00, 1.0e+01, 2.0e+01, 
     >           5.0e+01, 1.0e+02/
      data tz_t /1.0e+00, 2.0e+00, 3.0e+00, 4.0e+00, 5.0e+00, 
     >           6.0e+00, 7.0e+00, 8.0e+00, 9.0e+00, 1.0e+01, 
     >           2.0e+01, 3.0e+01, 4.0e+01, 5.0e+01, 6.0e+01, 
     >           7.0e+01, 8.0e+01, 9.0e+01, 1.0e+02, 2.0e+02, 
     >           5.0e+02, 1.0e+03/

c     branch_t = percentage of electrons captured into the 3S or 3D state
c                where they can make a transition to 2P
      data branch_t
     >     /9.49069e-01, 9.42006e-01, 9.02994e-01, 7.88099e-01,
     >      5.75389e-01, 2.97882e-01, 1.92816e-01,
     >      9.48021e-01, 9.40408e-01, 8.99590e-01, 7.83552e-01,
     >      5.72411e-01, 2.96504e-01, 1.92832e-01,
     >      9.46819e-01, 9.38704e-01, 8.96042e-01, 7.78878e-01,
     >      5.69229e-01, 2.96076e-01, 1.93337e-01,
     >      9.45511e-01, 9.36963e-01, 8.92454e-01, 7.74334e-01,
     >      5.66256e-01, 2.95263e-01, 1.93154e-01,
     >      9.44092e-01, 9.35079e-01, 8.88745e-01, 7.69835e-01,
     >      5.63175e-01, 2.94536e-01, 1.93019e-01,
     >      9.42616e-01, 9.33104e-01, 8.84950e-01, 7.65290e-01,
     >      5.59968e-01, 2.93697e-01, 1.92881e-01,
     >      9.41078e-01, 9.31020e-01, 8.81044e-01, 7.60680e-01,
     >      5.57045e-01, 2.92822e-01, 1.92742e-01,
     >      9.39422e-01, 9.28825e-01, 8.77042e-01, 7.56284e-01,
     >      5.54150e-01, 2.92021e-01, 1.92594e-01,
     >      9.37680e-01, 9.26512e-01, 8.72989e-01, 7.51702e-01,
     >      5.51222e-01, 2.91220e-01, 1.92447e-01,
     >      9.35845e-01, 9.24044e-01, 8.68887e-01, 7.47256e-01,
     >      5.48261e-01, 2.90442e-01, 1.92330e-01,
     >      9.11327e-01, 8.93185e-01, 8.24866e-01, 7.03916e-01,
     >      5.20335e-01, 2.82788e-01, 1.91032e-01,
     >      8.75469e-01, 8.52683e-01, 7.78973e-01, 6.63484e-01,
     >      4.94678e-01, 2.75737e-01, 1.89860e-01,
     >      8.32168e-01, 8.07463e-01, 7.33979e-01, 6.26036e-01,
     >      4.71192e-01, 2.69030e-01, 1.88730e-01,
     >      7.86251e-01, 7.61606e-01, 6.91406e-01, 5.91616e-01,
     >      4.49855e-01, 2.62867e-01, 1.87711e-01,
     >      7.41003e-01, 7.17460e-01, 6.51859e-01, 5.60016e-01,
     >      4.30090e-01, 2.57054e-01, 1.86745e-01,
     >      6.98039e-01, 6.76031e-01, 6.15371e-01, 5.31122e-01,
     >      4.11976e-01, 2.51669e-01, 1.85799e-01,
     >      6.57990e-01, 6.37637e-01, 5.81844e-01, 5.04595e-01,
     >      3.95288e-01, 2.46655e-01, 1.84972e-01,
     >      6.21041e-01, 6.02323e-01, 5.51144e-01, 4.80330e-01,
     >      3.79929e-01, 2.41900e-01, 1.84179e-01,
     >      5.87018e-01, 5.69851e-01, 5.22967e-01, 4.58184e-01,
     >      3.65788e-01, 2.37502e-01, 1.83448e-01,
     >      3.69337e-01, 3.61991e-01, 3.41731e-01, 3.13084e-01,
     >      2.70320e-01, 2.06189e-01, 1.78635e-01,
     >      1.99158e-01, 1.98220e-01, 1.95595e-01, 1.91797e-01,
     >      1.85988e-01, 1.77898e-01, 1.79075e-01,
     >      1.78397e-01, 1.78496e-01, 1.78807e-01, 1.79376e-01,
     >      1.80673e-01, 1.85499e-01, 1.95238e-01/

      t0 = min(max(t0arg, t0_t(1)), t0_t(dimt0))
      tz = min(max(tzarg, tz_t(1)), tz_t(dimtz))

      call linin2(t0_t,tz_t,branch_t,dimt0,dimtz,t0,tz,branchbuf)
      b4branch=branchbuf

      return

      end


c     ------------------------------------------------------------------

      real function b4atripl(t0arg, tzarg)
      implicit none

c     t0 = neutrals temperature (eV)
c     tz = impurity temperature (eV)
c     b4atripl = triplett charge exchange rate coefficient (m^3/s)

      integer dimt0, dimtz
      parameter (dimt0=7)
      parameter (dimtz=22)

      real t0_t(dimt0), tz_t(dimtz), atripl_t(dimt0,dimtz)
      save t0_t, tz_t, atripl_t

      real t0arg, tzarg, t0, tz, atriplbuf

      data t0_t /1.0e+00, 2.0e+00, 5.0e+00, 1.0e+01, 2.0e+01, 
     >           5.0e+01, 1.0e+02/
      data tz_t /1.0e+00, 2.0e+00, 3.0e+00, 4.0e+00, 5.0e+00, 
     >           6.0e+00, 7.0e+00, 8.0e+00, 9.0e+00, 1.0e+01, 
     >           2.0e+01, 3.0e+01, 4.0e+01, 5.0e+01, 6.0e+01, 
     >           7.0e+01, 8.0e+01, 9.0e+01, 1.0e+02, 2.0e+02, 
     >           5.0e+02, 1.0e+03/

c     atripl_t = charge exchange rate coefficient for process with
c                captured electron in P3 triplett state
      data atripl_t
     >     /2.36899e-15, 2.83932e-15, 3.26673e-15, 3.86757e-15,
     >      5.48030e-15, 1.10863e-14, 1.95567e-14,
     >      2.49591e-15, 2.88760e-15, 3.28575e-15, 3.89383e-15,
     >      5.51401e-15, 1.10982e-14, 1.96449e-14,
     >      2.60007e-15, 2.92856e-15, 3.30368e-15, 3.91942e-15,
     >      5.54627e-15, 1.11375e-14, 1.96712e-14,
     >      2.68604e-15, 2.96608e-15, 3.32404e-15, 3.94617e-15,
     >      5.57956e-15, 1.11721e-14, 1.96985e-14,
     >      2.75680e-15, 2.99884e-15, 3.34350e-15, 3.97306e-15,
     >      5.61244e-15, 1.12073e-14, 1.97292e-14,
     >      2.81749e-15, 3.02835e-15, 3.36326e-15, 3.99987e-15,
     >      5.64421e-15, 1.12413e-14, 1.97581e-14,
     >      2.86899e-15, 3.05522e-15, 3.38311e-15, 4.02611e-15,
     >      5.67737e-15, 1.12741e-14, 1.97871e-14,
     >      2.91329e-15, 3.08007e-15, 3.40204e-15, 4.05400e-15,
     >      5.71094e-15, 1.13090e-14, 1.98162e-14,
     >      2.95196e-15, 3.10337e-15, 3.42322e-15, 4.08124e-15,
     >      5.74410e-15, 1.13432e-14, 1.98446e-14,
     >      2.98617e-15, 3.12383e-15, 3.44422e-15, 4.10895e-15,
     >      5.77715e-15, 1.13772e-14, 1.98738e-14,
     >      3.21737e-15, 3.32016e-15, 3.66836e-15, 4.39313e-15,
     >      6.11076e-15, 1.17179e-14, 2.01612e-14,
     >      3.41149e-15, 3.52607e-15, 3.92005e-15, 4.69219e-15,
     >      6.44883e-15, 1.20576e-14, 2.04458e-14,
     >      3.63062e-15, 3.76108e-15, 4.19229e-15, 5.00113e-15,
     >      6.79045e-15, 1.23930e-14, 2.07322e-14,
     >      3.87794e-15, 4.02124e-15, 4.48152e-15, 5.31950e-15,
     >      7.13645e-15, 1.27284e-14, 2.10146e-14,
     >      4.14782e-15, 4.30085e-15, 4.78350e-15, 5.64434e-15,
     >      7.48317e-15, 1.30608e-14, 2.12958e-14,
     >      4.43452e-15, 4.59506e-15, 5.09553e-15, 5.97595e-15,
     >      7.83193e-15, 1.33921e-14, 2.15726e-14,
     >      4.73445e-15, 4.90096e-15, 5.41597e-15, 6.31211e-15,
     >      8.18152e-15, 1.37207e-14, 2.18522e-14,
     >      5.04517e-15, 5.21658e-15, 5.74365e-15, 6.65341e-15,
     >      8.53216e-15, 1.40472e-14, 2.21301e-14,
     >      5.36435e-15, 5.53983e-15, 6.07690e-15, 6.99755e-15,
     >      8.88222e-15, 1.43708e-14, 2.24052e-14,
     >      8.79131e-15, 8.98033e-15, 9.54641e-15, 1.04842e-14,
     >      1.23257e-14, 1.75059e-14, 2.50717e-14,
     >      1.85686e-14, 1.87285e-14, 1.92053e-14, 1.99894e-14,
     >      2.15204e-14, 2.58474e-14, 3.23590e-14,
     >      3.15656e-14, 3.16902e-14, 3.20624e-14, 3.26778e-14,
     >      3.38901e-14, 3.73948e-14, 4.28731e-14/

      t0 = min(max(t0arg, t0_t(1)), t0_t(dimt0))
      tz = min(max(tzarg, tz_t(1)), tz_t(dimtz))

      call linin2(t0_t,tz_t,atripl_t,dimt0,dimtz,t0,tz,atriplbuf)
      b4atripl=atriplbuf

      return

      end
c
c
c
c      SUBROUTINE ADASRD(YEAR,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)     
c      SUBROUTINE ADASRD(YEAR,YEARDF,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)
C
C  READ THE REQUESTED RATE COEFFICIENT FROM THE ADAS MASTER ELEMENT
C  FILES:
C        ICLASS = 1: RECOMBINATION RATE COEFFICIENT
C                 2: IONISATION RATE COEFFICIENT
C                 3: CHARGE EXCHANGE RECOMBINATION COEFFICIENT
C                 4: POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLU
C                 5: POWER COEFFICIENT FOR LINE RADIATION
C                 6: POWER COEFFICIENT FOR CHARGE EXCHANGE
C                 added by Krieger, IPP 5/95
C                 7: PHOTON EMISSIVITY FOR DIAGNOSTIC LINES
C  THIS ROUTINE USES THE STANDARD ADAS EXTRACTION ROUTINE D2DATA AND
C  REALLY ONLY PROVIDES A 'CLEAN' INTERFACE, TAKING CARE OF CHANGES
C  IN UNITS AND IN PRECISION OF VARIABLES.  IF THE REQUESTED DATA
C  DOESN'T EXIST (IFAIL=1 RETURNED FROM D2DATA) THE PROGRAM IS STOPPED.
C
c      IMPLICIT NONE
C     INCLUDE   "PARAMS"
c      include    'params'
C     INCLUDE   "CADAS2"
c      include    'cadas2'
C
c      CHARACTER*2 YEAR
c      INTEGER IZ0, IZ1, ICLASS, NPTS
c      REAL TE(NPTS), NE(NPTS), COEF(NPTS)
c
c      logical lintrp(maxpts)
C
c      INTEGER I, J
C
c      IEVCUT = 0
C 
c      DO J = 1, NPTS
c
c        DTEV(J) = DBLE(ALOG10(TE(J)))
c        DDENS(J) = DBLE(ALOG10(NE(J)*1.0E-6))
c
c        write(6,*) 'adasrd:',te(j),ne(j),iz0,iz1,year
c        write(6,*) 'adasrd:',npts,iclass,maxads,titlf
c
c      ENDDO
c      CALL D2DATA(YEAR, YEAR, TITLF, IFAIL,
c     >            IZ0, IZ1, ICLASS, NPTS, IEVCUT,
c     >            MAXADS, ITMAXD, IDMAXD, IZMAXD,
c     >            DTEV, DDENS,
c     >            DTEVD, DDENSD, DRCOFD, ZDATA,
c     >            DRCOFI
c     >            , LINTRP
c     >            )
c      IF (IFAIL.EQ.1) THEN
c        WRITE(6,1000) IZ0, IZ1,YEAR
c        STOP
c      ENDIF
C
C  EXTRAPOLATED VALUES ARE RETURNED AS ZERO!
C
c      DO J = 1, NPTS
c        IF (DRCOFI(J).NE.0.0) THEN
c          COEF(J) = 10.**SNGL(DRCOFI(J)) * 1.0E-6
c        ELSE
c
c          write (6,*) 'Extrapolated coefficient:',j,te(j),ne(j)
c
c          COEF(J) = 0.0
c        ENDIF
c      ENDDO
C
c 1000 FORMAT(' ERROR READING REQUESTED ATOMIC DATA!',/,
c     >       ' MASTER ELEMENT FILE FOR NUCLEAR CHARGE ',I2,
c     >       ' AND ION CHARGE ',I2,
c     >       ' WAS NOT FOUND IN YEAR ',A2)
C
c      RETURN
c      END
c 
