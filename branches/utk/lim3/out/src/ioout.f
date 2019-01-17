      SUBROUTINE RDG (GRAPH,VMIN,VMAX,GRIMIN,GRIMAX,
     >                IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,     
     >                IALLIZ,IVU,NAME,IERR)                                     
      use mod_reader
      IMPLICIT  none
      REAL      VMIN,VMAX                                                       
      REAL      GRIMIN,GRIMAX
      INTEGER   IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,IALLIZ,IVU,IERR                 
      CHARACTER NAME*(*),GRAPH*(*)                                              
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS ROUTINE READS IN A LINE OF GRAPH DETAILS                         
C                                                                               
C     PARAMETERS :-                                                             
C     GRAPH  : TITLE OF GRAPH                                                   
C     VMIN,VMAX : PLOTTING LIMITS                   
C     GRIMIN,GRIMAX : INTEGRATION LIMITS                             
C     IPLOT ... IVU : PLOTTING OPTIONS   (7 INTEGERS)                           
C     NAME   : FOR PRINTING IN ANY ERROR MESSAGES                               
C     IERR   : SET TO 1 IF AN ERROR FOUND                                       
C                                                                               
C       DON'T WANT EOF TO BE TREATED AS AN ERROR IN THIS CASE.  STILL           
C     SET IERR TO 1 BUT DON'T PRINT ANY MESSAGES (HENCE GOTO 9998)              
C                                                                               
C        CHRIS FARRELL    JAN 1988                                              
C
C       ADDED INTEGRATION RANGE DATA. VALUES OF 0.0 TO 0.0 ARE 
C     GENERALLY TREATED AS THE FULL RANGE.
C
C       DAVID ELDER,  OCT 29, 1990
C
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'reader'                                                        
C     INCLUDE   (READER)                                                        
      CHARACTER COMENT*72,MESAGE*72                                             
C                                                                               
      MESAGE = 'END OF FILE ON UNIT 5'                                          
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER                  
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG'                                   
      IF (BUFFER(1:1).EQ.'$') GOTO 100                                          
C                                                                               
      MESAGE = 'EXPECTING CHARACTER STRING, 4 REALS AND 7 INTEGERS'             
C
C     REPLACE LIST-DIRECTED I/O ON INTERNAL FILES WITH A BACKSPACE
C     FOLLOWED BY A READ ON THE EXTERNAL FILE THAT IS THE SOURCE FOR 
C     THE BUFFER. THIS IS DONE BECAUSE AT THIS POINT IN TIME THE CFT 
C     AND CFT77 COMPILERS ON THE CRAY DO NOT SUPPORT THIS FEATURE.
C     
C                      DAVID ELDER NOV.15,1989
C
C     READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
C
      BACKSPACE(5)
      READ (5,*,ERR=9999,END=9999) GRAPH                                   
      IF (GRAPH(1:1).NE.'2'.and.graph(1:3).ne.'LOS'.and.
     >    graph(1:3).ne.'000') THEN                       
        IBUF = 1                                                                
      ELSE                                                                      
        IBUF = 0                                                                
C
C       READ (BUFFER,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,IPLOT,JSMOTH,         
C    >    MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                                         
C
        BACKSPACE(5)
        READ (5,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,
     >    GRIMIN,GRIMAX,IPLOT,JSMOTH,         
     >    MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                                         
      ENDIF                                                                     
      RETURN                                                                    
C                                                                               
 9998 IERR = 1                                                                  
      RETURN                                                                    
C                                                                               
 9999 IERR = 1                                                                  
      WRITE (7,'(1X,2A,3(/1X,A))')                                              
     >  'RDG: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER            
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE RDGT(GRAPH,VMIN,VMAX,IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,     
     >                IALLIZ,IVU,NAME,IERR)                                     
      use mod_reader
      IMPLICIT  none
      REAL      VMIN,VMAX                                                       
      INTEGER   IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,IALLIZ,IVU,IERR                 
      CHARACTER NAME*(*),GRAPH*(*)                                              
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS ROUTINE READS IN A LINE OF GRAPH DETAILS                         
C                                                                               
C     PARAMETERS :-                                                             
C     GRAPH  : TITLE OF GRAPH                                                   
C     VMIN,VMAX : PLOTTING LIMITS                                               
C     IPLOT ... IVU : PLOTTING OPTIONS   (7 INTEGERS)                           
C     NAME   : FOR PRINTING IN ANY ERROR MESSAGES                               
C     IERR   : SET TO 1 IF AN ERROR FOUND                                       
C                                                                               
C       DON'T WANT EOF TO BE TREATED AS AN ERROR IN THIS CASE.  STILL           
C     SET IERR TO 1 BUT DON'T PRINT ANY MESSAGES (HENCE GOTO 9998)              
C                                                                               
C        CHRIS FARRELL    JAN 1988                                              
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'reader'                                                        
C     INCLUDE   (READER)                                                        
      CHARACTER COMENT*72,MESAGE*72                                             
C                                                                               
      MESAGE = 'END OF FILE ON UNIT 5'                                          
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER                  
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGT'                                  
      IF (BUFFER(1:1).EQ.'$') GOTO 100                                          
C                                                                               
      MESAGE = 'EXPECTING CHARACTER STRING, 2 REALS AND 7 INTEGERS'             
C
C     REPLACE LIST-DIRECTED I/O ON INTERNAL FILES WITH A BACKSPACE
C     FOLLOWED BY A READ ON THE EXTERNAL FILE THAT IS THE SOURCE FOR 
C     THE BUFFER. THIS IS DONE BECAUSE AT THIS POINT IN TIME THE CFT 
C     AND CFT77 COMPILERS ON THE CRAY DO NOT SUPPORT THIS FEATURE.
C     
C                      DAVID ELDER NOV.15,1989
C
C     READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
C
      BACKSPACE(5)
      READ (5,*,ERR=9999,END=9999) GRAPH                                   
      IF (GRAPH(1:1).NE.'T'.and.graph(1:3).ne.'000') THEN                       
        IBUF = 1                                                                
      ELSE                                                                      
        IBUF = 0                                                                
C
C       READ (BUFFER,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,IPLOT,JSMOTH,         
C    >    MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                                         
C
        BACKSPACE(5)
        READ (5,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,IPLOT,JSMOTH,         
     >    MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                                         
      ENDIF                                                                     
      RETURN                                                                    
C                                                                               
 9998 IERR = 1                                                                  
      RETURN                                                                    
C                                                                               
 9999 IERR = 1                                                                  
      WRITE (7,'(1X,2A,3(/1X,A))')                                              
     > 'RDGT: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER           
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE RDGRT (GRAPH,VMIN,VMAX,SMIN,SMAX,IPLOT,JSMOTH,MAXIZ
     >                  ,IPLANE,IFOLD,IALLIZ,IVU,NAME,IERR)                     
      use mod_reader
      IMPLICIT  none
      REAL      VMIN,VMAX ,SMIN,SMAX                                         
      INTEGER   IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,IALLIZ,IVU,IERR                 
      CHARACTER NAME*(*),GRAPH*(*)                                              
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS ROUTINE READS IN A LINE OF GRAPH DETAILS                         
C                                                                               
C     PARAMETERS :-                                                             
C     GRAPH  : TITLE OF GRAPH                                                   
C     VMIN,VMAX : PLOTTING LIMITS                                               
C     SMIN,SMAX : INTEGRATING RANGE
C     IPLOT ... IVU : PLOTTING OPTIONS   (7 INTEGERS)                           
C     NAME   : FOR PRINTING IN ANY ERROR MESSAGES                               
C     IERR   : SET TO 1 IF AN ERROR FOUND                                       
C                                                                               
C       DON'T WANT EOF TO BE TREATED AS AN ERROR IN THIS CASE.  STILL           
C     SET IERR TO 1 BUT DON'T PRINT ANY MESSAGES (HENCE GOTO 9998)              
C                                                                               
C        DAVID ELDER      JULY 1990                                         
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'reader'                                                        
C     INCLUDE   (READER)                                                        
      CHARACTER COMENT*72,MESAGE*72                                             
      REAL DEGRAD
      DATA DEGRAD /.017453292/
C                                                                               
      MESAGE = 'END OF FILE ON UNIT 5'                                          
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER                  
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGRT'                                
      IF (BUFFER(1:1).EQ.'$') GOTO 100                                          
C                                                                               
      MESAGE = 'EXPECTING CHARACTER STRING, 4 REALS AND 7 INTEGERS'             
C
C     REPLACE LIST-DIRECTED I/O ON INTERNAL FILES WITH A BACKSPACE
C     FOLLOWED BY A READ ON THE EXTERNAL FILE THAT IS THE SOURCE FOR 
C     THE BUFFER. THIS IS DONE BECAUSE AT THIS POINT IN TIME THE CFT 
C     AND CFT77 COMPILERS ON THE CRAY DO NOT SUPPORT THIS FEATURE.
C     
C                      DAVID ELDER NOV.15,1989
C
C     READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
C
      BACKSPACE(5)
      READ (5,*,ERR=9999,END=9999) GRAPH                                   
      IF (GRAPH(1:1).NE.'R'.and.graph(1:3).ne.'000')  THEN                         
        IBUF = 1                                                                
      ELSE                                                                      
        IBUF = 0                                                                
C
C       READ (BUFFER,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,IPLOT,JSMOTH,         
C    >    MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                                         
C
        BACKSPACE(5)
        READ (5,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,SMIN,SMAX,
     >        IPLOT,JSMOTH,MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                      
      ENDIF                                                                     
C
C     ALL THE ANGULAR VALUES READ IN BY THIS ROUTINE ARE EXPECTED
C     IN DEGREES. HOWEVER,  I FELT IT WOULD BE MORE CONSISTENT
C     AS WELL AS SAFER TO MAINTAIN THESE VALUES IN RADIANS 
C     INTERNALLY UNTIL BEFORE THE PLOTTING ROUTINES WHERE THEY 
C     CAN BE RE-CONVERTED BACK TO DEGREES. THIS MEANS THAT ALL
C     ANGULAR VALUES REQUIRED FOR CALCULATIONS CAN BE USED DIRECTLY.
C     THE RECONVERSIONS ARE PERFORMED EITHER IN THE ROUTINE CVRTDEG
C     OR , FOR PSI PLOTS, DIRECTLY IN THE ROUTINE RINTPSI. AT THIS 
C     TIME, THE COORDINATE ARRAY SOUTS FOR PSI PLOTS IS ONLY 
C     MAINTAINED IN DEGREES.
C
C     DAVID ELDER  , AUG 3 1990
C
      IF (GRAPH(3:3).EQ.'A'.OR.GRAPH(3:3).EQ.'S') THEN
         VMIN = VMIN * DEGRAD
         VMAX = VMAX * DEGRAD
      ELSEIF (GRAPH(3:3).EQ.'R') THEN
         SMIN = SMIN * DEGRAD
         SMAX = SMAX * DEGRAD
      ENDIF 
      RETURN                                                                    
C                                                                               
 9998 IERR = 1                                                                  
      RETURN                                                                    
C                                                                               
 9999 IERR = 1                                                                  
      WRITE (7,'(1X,2A,3(/1X,A))')                                              
     >'RDGRT: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER         
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE RDG3D (GRAPH,XMIN,XMAX,YMIN,YMAX,NPTS,ISTATE,IPLANE,           
     >                  IFOLD,JSMOTH,NAME,IERR)                                        
      use mod_reader
      IMPLICIT  none
      REAL      XMIN,XMAX,YMIN,YMAX                                             
      INTEGER   NPTS,ISTATE,IPLANE,IFOLD,JSMOTH,IERR                                   
      CHARACTER NAME*(*),GRAPH*(*)                                              
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS ROUTINE READS IN A LINE OF 3D GRAPH DETAILS                      
C                                                                               
C     PARAMETERS :-                                                             
C     GRAPH  : TITLE OF GRAPH                                                   
C     XMIN,XMAX,YMIN,YMAX : PLOTTING LIMITS                                     
C     NPTS ... IFOLD : PLOTTING OPTIONS   (4 INTEGERS)                          
C     NAME   : FOR PRINTING IN ANY ERROR MESSAGES                               
C     IERR   : SET TO 1 IF AN ERROR FOUND                                       
C                                                                               
C       DON'T WANT EOF TO BE TREATED AS AN ERROR IN THIS CASE.  STILL           
C     SET IERR TO 1 BUT DON'T PRINT ANY MESSAGES (HENCE GOTO 9998)              
C                                                                               
C        CHRIS FARRELL    FEB 1988                                              
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'reader'                                                        
C     INCLUDE   (READER)                                                        
      CHARACTER COMENT*72,MESAGE*72                                             
C                                                                               
      MESAGE = 'END OF FILE ON UNIT 5'                                          
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER                  
c      write(0,*) '3D:',trim(buffer)

      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG3D'                                 
C
      IF (BUFFER(1:1).EQ.'$') GOTO 100                                          
C                                                                               
      MESAGE = 'EXPECTING CHARACTER STRING, 4 REALS AND 5 INTEGERS'             
C
C     REPLACE LIST-DIRECTED I/O ON INTERNAL FILES WITH A BACKSPACE
C     FOLLOWED BY A READ ON THE EXTERNAL FILE THAT IS THE SOURCE FOR 
C     THE BUFFER. THIS IS DONE BECAUSE AT THIS POINT IN TIME THE CFT 
C     AND CFT77 COMPILERS ON THE CRAY DO NOT SUPPORT THIS FEATURE.
C     
C                      DAVID ELDER NOV.15,1989
C
C     Returned to original form for workstation execution. The rest 
C     will be changed when convenient.
C
C     David Elder  Mar 2, 1992
C
C
C      READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
C
c      BACKSPACE(5)
c slmod - all kind of trouble here

c      READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
c
c     jdemod - wtf - backspace was commented out resulting in skipped lines
c                    and breaking plot input
      BACKSPACE(5)
      READ (5,*,ERR=9999,END=9999) GRAPH                                   
C
      IF (GRAPH(1:1).NE.'3'.AND.GRAPH(1:1).NE.'C'.and.
     >    graph(1:3).ne.'000') THEN                         
        IBUF = 1                                                                
      ELSE                                                                      
        IBUF = 0                                                                
C
c        READ (BUFFER,*,ERR=9999,END=9999) GRAPH,XMIN,XMAX,YMIN,YMAX,            
c     >    NPTS,ISTATE,IPLANE,IFOLD,JSMOTH                                              
C
        BACKSPACE(5) 
C
        READ (5,*,ERR=9999,END=9999) GRAPH,XMIN,XMAX,YMIN,YMAX,            
     >    NPTS,ISTATE,IPLANE,IFOLD,JSMOTH                                              
C
      ENDIF                                                                     
      RETURN                                                                    
C                                                                               
 9998 IERR = 1                                                                  
      RETURN                                                                    
C                                                                               
 9999 IERR = 1                                                                  
c      WRITE (7,'(1X,2A,3(/1X,A))')                                              
c     >  'RDG3D: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER          
      WRITE (7,'(4A,A,4F6.2,5I3)')                                              
     >  'RDG3D: ERROR READING ',NAME,MESAGE,' LAST LINE READ :-',
     +  GRAPH,XMIN,XMAX,YMIN,YMAX,NPTS,ISTATE,IPLANE,IFOLD,JSMOTH 
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RDGM (GRAPH,XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NPTS,MPTS,
     >            ISTATE,IPLANE,IPLOT,IFOLD,NAME,IERR)                       
      use mod_reader
      IMPLICIT  none
      REAL      XMIN,XMAX,YMIN,YMAX,PMIN,PMAX                                   
      INTEGER   NPTS,MPTS,ISTATE,IPLANE,IPLOT,IFOLD,IERR                   
      CHARACTER NAME*(*),GRAPH*(*)                                              
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS ROUTINE READS IN A LINE OF MESH/CONTOUR GRAPH DETAILS            
C                                                                               
C     PARAMETERS :-                                                             
C     GRAPH  : TITLE OF GRAPH                                                   
C     XMIN,XMAX,YMIN,YMAX,PMIN,PMAX : PLOTTING LIMITS                      
C     NPTS ... IFOLD : PLOTTING OPTIONS   (4 INTEGERS)                          
C     NAME   : FOR PRINTING IN ANY ERROR MESSAGES                               
C     IERR   : SET TO 1 IF AN ERROR FOUND                                       
C                                                                               
C       DON'T WANT EOF TO BE TREATED AS AN ERROR IN THIS CASE.  STILL           
C     SET IERR TO 1 BUT DON'T PRINT ANY MESSAGES (HENCE GOTO 9998)              
C                                                                               
C        CHRIS FARRELL    FEB 1988, NEW ROUTINE MARCH 15/90 D. ELDER          
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'reader'                                                        
C     INCLUDE   (READER)                                                        
      CHARACTER COMENT*72,MESAGE*72                                             
C                                                                               
      MESAGE = 'END OF FILE ON UNIT 5'                                          
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER                  
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGM'                                 
c      write(0,*) 'M :',trim(buffer)
      IF (BUFFER(1:1).EQ.'$') GOTO 100                                          
C                                                                               
      MESAGE = 'EXPECTING CHARACTER STRING, 2 INTS, 6 REALS, 4 INTS'             
C
C     REPLACE LIST-DIRECTED I/O ON INTERNAL FILES WITH A BACKSPACE
C     FOLLOWED BY A READ ON THE EXTERNAL FILE THAT IS THE SOURCE FOR 
C     THE BUFFER. THIS IS DONE BECAUSE AT THIS POINT IN TIME THE CFT 
C     AND CFT77 COMPILERS ON THE CRAY DO NOT SUPPORT THIS FEATURE.
C     
C                      DAVID ELDER NOV.15,1989
C
C     READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
      BACKSPACE(5)
      READ (5,*,ERR=9999,END=9999) GRAPH                                   
      IF (GRAPH(1:1).NE.'M'.and.graph(1:3).ne.'000') THEN                         
        IBUF = 1                                                                
      ELSE                                                                      
        IBUF = 0                                                                
C
C       READ (BUFFER,*,ERR=9999,END=9999) GRAPH,XMIN,XMAX,YMIN,YMAX,            
C    >    NPTS,ISTATE,IPLANE,IFOLD                                              
C
        BACKSPACE(5) 
        READ (5,*,ERR=9999,END=9999) GRAPH,IPLOT,IPLANE,XMIN,
     >    XMAX,YMIN,YMAX,PMIN,PMAX,NPTS,MPTS,ISTATE,IFOLD
C
C       IT MAY BE NECESSARY TO HAVE 2 LINES OF GRAPH DATA FOR THESE
C       CASES. THE FOLLOWING STATEMENT COULD READ THAT SECOND LINE.
C
C       READ (5,*,ERR=9999,END=9999) NPTS,MPTS,ISTATE,IPLANE,IFOLD
C
      ENDIF                                                                     
      RETURN                                                                    
C                                                                               
 9998 IERR = 1                                                                  
      RETURN                                                                    
C                                                                               
 9999 IERR = 1                                                                  
      WRITE (7,'(1X,2A,3(/1X,A))')                                              
     >  'RDGM: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER          
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE COLECT (TITLE,NIZS,NIN,IERR,JOB,IMODE,PLAMS,PIZS,NLS,          
     >                   FACTA,FACTB,ITER,NITERS)                               
      use mod_params
      use mod_comtor
      use mod_dynam2
      use mod_dynam3
      use mod_comt2
      use mod_comnet
      use mod_comxyt
      use mod_coords
C     
C  *********************************************************************        
C  *                                                                   *        
C  *  COLECT:  FETCH RESULTS OF LIM RUN FROM UNFORMATTED FILE "NIN".   *        
C  *           THIS FILE SHOULD BE REWOUND BEFORE COLECT IS CALLED     *        
C  *           FOR THE FIRST TIME.                                     *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      IMPLICIT  none
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'dynam2'                                                        
C     INCLUDE   (DYNAM2)                                                        
c      INCLUDE   'dynam3'                                                        
C     INCLUDE   (DYNAM3)                                                        
      CHARACTER TITLE*80,JOB*72                                                 
      INTEGER   NIZS,IMODE,NLS,NIN,IERR,ITER,NITERS                             
      REAL      PLAMS(MAXNLS),FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)                 
      INTEGER   PIZS(MAXNLS)                                                    
C                                                                               
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
      INTEGER   IX,IY,IP,IT,IOS,IERR2,JBLOCK,IYB,IYE,IZ,J                       
      INTEGER   IBLOCK,II,IL,IZS,IZE,IQX,KBLOCK,NOS,IO,IQS,IQE                  
      CHARACTER VERSE*5                                                         
      DATA      IBLOCK / 1540 /                                                 
c
c     Version control - calculate a unique and always increasing
c     version code.
c     Maximum revison number for a given version number is maxrev-1 
c
      integer   maxrev,version_code
      integer   vernum, revnum
      parameter (maxrev=100)
c
c  

      write(0,*) 'READING INPUT FILE:'
c
c
C                                                                               
C-----------------------------------------------------------------------        
C     READ  COMMONS, LINE DATA, SHORT ARRAYS, PARAMETERS, ETC                   
C-----------------------------------------------------------------------        
C                                                                               
      READ  (NIN ,IOSTAT=IOS) VERSE,NY3D,ITER,NITERS,NOS                        
c
c     LIM has the main version number in a different location  
c
      write (0,*) 'VERSION:',':',ios,':',trim(verse),':'

      read(verse,'(1x,i1,1x,i2)') vernum,revnum
c
      version_code = vernum * maxrev + revnum
c
      IF (VERSON.NE.VERSE .AND. ITER.EQ.1) THEN                                 
        CALL PRC ('WARNING! PROGRAM HAS BEEN RE-COMPILED SINCE LIM RUN')        
        WRITE (7,'(10X,''ORIGINAL JOB PRODUCED BY LIM'',A5)') VERSE             
        WRITE (7,'(10X,''GRAPHICS PROGRAM IN USE  OUT'',A5)') VERSON            
      ENDIF                                                                     
c
      IF (NY3D.LE.1 .AND. ITER.EQ.1) THEN                                       
        CALL PRC ('WARNING! DUMPFILE CREATED FROM SMALL LIM VERSION')           
        CALL PRC ('         CONTAINS NO 3D OR TIME-DEPENDENT RESULTS.')         
      ENDIF                                                                     
C                                                                               
C---- PARAMETERS NEXT                                                           
C                                                                               
      READ  (NIN ,IOSTAT=IOS)                                                   
     >       NXS,NYS,NQXSO,NQXSI,NQYS,NTS,NIZS,NLS,TITLE,JOB,IMODE              
      WRITE (6,9002)                                                            
     >       NXS,NYS,NQXSO,NQXSI,NQYS,NTS,NIZS,NLS,TITLE,JOB,IMODE,ITER         
      IF (NXS   .GT.MAXNXS) IERR = IERR2 ('NXS',   NXS,   MAXNXS)               
      IF (NYS   .GT.MAXNYS) IERR = IERR2 ('NYS',   NYS,   MAXNYS)               
      IF (NQXSO .GT.MAXQXS) IERR = IERR2 ('NQXSO', NQXSO, MAXQXS)               
      IF (NQXSI .GT.MAXQXS) IERR = IERR2 ('NQXSI', NQXSI, MAXQXS)               
      IF (NQYS  .GT.MAXQYS) IERR = IERR2 ('NQYS',  NQYS,  MAXQYS)               
      IF (NTS   .GT.MAXNTS) IERR = IERR2 ('NTS',   NTS,   MAXNTS)               
      IF (NIZS  .GT.MAXIZS) IERR = IERR2 ('NIZS',  NIZS,  MAXIZS)               
      IF (NLS   .GT.MAXNLS) IERR = IERR2 ('NLS',   NLS,   MAXNLS)               
      IF (NY3D  .GT.MAXY3D) IERR = IERR2 ('NY3D',  NY3D,  MAXY3D)               
      IF (NOS   .GT.MAXOS ) IERR = IERR2 ('NOS',   NOS,   MAXOS )               
c
c     Read in absolute scaling factor
c
      if (version_code.ge.3*maxrev+1) then 
         read(nin,iostat=ios) absfac
      endif
C                                                                               
C---- SOME OF COMTOR COMMON, SOME OF COMXYT COMMON, MISCELLANEOUS               
C                                                                               
      READ  (NIN ,IOSTAT=IOS)                                                   
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

      do iy = 1,nys
         write(6,'(a,2i5,10(1x,g12.5))') 'YS:',iy,nys,ys(iy)
      end do
      do ix = 1,nxs
         write(6,'(a,2i5,10(1x,g12.5))') 'XS:',ix,nxs,xs(ix)
      end do
      do ip = -maxnps,maxnps
         write(6,'(a,2i5,10(1x,g12.5))') 'PS:',ip,maxnps,ps(ip)
      end do



      if (version_code.ge.3*maxrev+5) then  
         read(nin,iostat=ios) (pzone(ip),ip=-maxnps,maxnps)
      else
         ! if pzone isn't available then set every poloidal slice
         ! to be considered associated with a surface
         pzone = 1
      endif
c     
c     Read in some 3D option information
c
c     CIOPTJ= 3D limiter extent option 
c     CPCO  = 3D extent of limiter
c
c
      if (version_code.ge.3*maxrev+3) then  
         read(nin) cioptj,cpco
      endif

C                                                                               
C---- SHORT ARRAYS ... BLOCKED I/O USED                                         
C                                                                               
      DO 45 IZS = -2, NIZS+1, 15                                                
        IZE = MIN (IZS+15-1, NIZS+1)                                            
        READ  (NIN ,IOSTAT=IOS) ((SAVES(IX,IZ),IX=1,NXS),IZ=IZS,IZE)            
   45 CONTINUE                                                                  
C                                                                               
      DO 50 J = 1, 3                                                            
       DO 50 IZS = 1, NIZS +1, 15                                                  
        IZE = MIN (IZS+15-1, NIZS)                                              
        READ  (NIN ,IOSTAT=IOS) ((DEPS(IX,IZ,J),IX=1,NXS),IZ=IZS,IZE)           
   50 CONTINUE                                                                  
C                                                                               
      READ  (NIN ,IOSTAT=IOS)                                                   
     >  (((NEROXS(IX,II,J),IX=1,NXS),II=1,5),J=1,3)                             
C                                                                               
      DO 55 II = 1, 6                                                           
        READ  (NIN ,IOSTAT=IOS) (NEROYS(IO,II),IO=1,MAXOS)                      
   55 CONTINUE                                                                  
C                                                                               
      DO 57 II = 1, 5                                                           
        READ  (NIN ,IOSTAT=IOS) (NERODS(IO,II),IO=1,MAXOS)                      
   57 CONTINUE                                                                  
c
      if (version_code.ge.3*maxrev+2) then  
         DO II = 1, 6                                                           
            do ip = -maxnps,maxnps
               READ (NIN,IOSTAT=IOS) (NERODS3(IO,IP,II),IO=1,MAXOS)                      
            end do
         end do
      endif
C                                                                               
      DO 60 IZS = -2, NIZS+1, 7                                                 
        IZE = MIN (IZS+7-1, NIZS+1)                                             
        READ  (NIN ,IOSTAT=IOS) ((WALLS(IY,IZ),IY=-NYS,NYS),IZ=IZS,IZE)         
   60 CONTINUE                                                                  
C                                                                               
      READ  (NIN ,IOSTAT=IOS) (OYS(IO),ODS(IO),IO=1,MAXOS)                      
      READ  (NIN ,IOSTAT=IOS) (OYOUTS(IO),ODOUTS(IO),IO=1,MAXOS)                
      READ  (NIN ,IOSTAT=IOS) (OYWIDS(IO),ODWIDS(IO),IO=1,MAXOS)                
      DO 62 II= 1,3
         READ (NIN,IOSTAT=IOS) (CDFLUX(IO,II),IO=1,MAXOS)
  62  CONTINUE
C                                                                               
C---- CHECK FOR ERRORS                                                          
C                                                                               
      IF (IOS.NE.0) STOP                                                        
      JBLOCK = IBLOCK / (NXS+1)                                                 
      KBLOCK = JBLOCK / 2                                                       
C                                                                               
C================= READ  SDLIMS ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 100 IZ = -1, NIZS                                                     
        DO 100 IYB = -NYS, NYS, KBLOCK                                          
         IYE = MIN (IYB+KBLOCK-1, NYS)                                          
c
c         WRITE(6,'(a,3I5)') 'SDLIMS:',iz,iyb,iye
c
         READ (NIN) ((SDLIMS(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                  
c         WRITE(6,'(10G13.5)') 
c     >        ((SDLIMS(IX,IY,IZ),IX=1,NXS),IY=IYB,IYE)
c
  100  CONTINUE                                                                 
      ENDIF                                                                     


c
c      do iz = 0,nizs
c         do ix = 1,nxs
c            do iy = 1,nys
c              if (sdlims(ix,iy,iz).gt.0.0) 
c     >          write(6,'(a,3i5,5(1x,g12.5))') 'SDLIMS:',ix,iy,iz,
c     >          xs(ix),xs(ix+1),ys(iy),ys(iy+1),sdlims(ix,iy,iz)
c            end do
c         end do
c      end do
c
C                                                                               
C================= READ  SDTS ARRAY FROM DISC ========================        
C                                                                               
      if (version_code.ge.3*maxrev+1) then  
c
         IF (IMODE.NE.1) THEN                                                      
          DO IZ =  1, NIZS                                                     
           DO IYB = -NYS, NYS, KBLOCK                                          
            IYE = MIN (IYB+KBLOCK-1, NYS)                                          
            READ (NIN) ((SDTS(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                  
           end do
          end do
         ENDIF                                                                     
c
      endif 
C                                                                               
C================= READ  POWLS  ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 200 IZ = -1, NIZS                                                     
        DO 200 IYB = -NYS, NYS, JBLOCK                                          
         IYE = MIN (IYB+JBLOCK-1, NYS)                                          
         READ (NIN) ((POWLS(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                   
  200  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= READ  LINES  ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 300 IZ = -1, NIZS                                                     
        DO 300 IYB = -NYS, NYS, JBLOCK                                          
         IYE = MIN (IYB+JBLOCK-1, NYS)                                          
         READ (NIN) ((LINES(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                   
  300  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= READ  TIZS   ARRAY FROM DISC ========================        
C                                                                               
      DO 500 IZ = -1, NIZS                                                      
       DO 500 IYB = -NYS, NYS, JBLOCK                                           
        IYE = MIN (IYB+JBLOCK-1, NYS)                                           
        READ (NIN) ((TIZS(IX,IY,IZ), IX=1,NXS), IY=IYB,IYE)                     
  500 CONTINUE                                                                  
C                                                                               
C================= READ  ZEFFS  ARRAY FROM DISC ========================        
C                                                                               
      DO 600 II = 1, 6                                                          
       DO 600 IYB = -NYS, NYS, JBLOCK                                           
        IYE = MIN (IYB+JBLOCK-1, NYS)                                           
        READ (NIN) ((ZEFFS(IX,IY,II), IX=1,NXS), IY=IYB,IYE)                    
  600 CONTINUE                                                                  
C                                                                               
C================= READ  SDLIM3 ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 700 IP = -MAXNPS, MAXNPS                                              
        DO 700 IZ = -1, NIZS                                                    
         DO 700 IYB = -NY3D, NY3D, KBLOCK                                       
          IYE = MIN (IYB+KBLOCK-1, NY3D)                                        
          READ (NIN) ((SDLIM3(IX,IY,IZ,IP), IX=1,NXS), IY=IYB,IYE)              
  700  CONTINUE                                                                 
      ENDIF                                                                     
c
c      WRITE(6,'(10G13.5)') ((SDLIM3(1,IY,1,IP),IY=-NY3D,NY3D),
c    >           IP=-MAXNPS,MAXNPS) 
C                                                                               
C================= READ  TIZ3   ARRAY FROM DISC ========================        
C                                                                               
      IF (CDPOL.GT.0.0) THEN                                                    
       DO 1100 IP = -MAXNPS, MAXNPS                                             
        DO 1100 IZ = -1, NIZS                                                   
         DO 1100 IYB = -NY3D, NY3D, JBLOCK                                      
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          READ (NIN) ((TIZ3(IX,IY,IZ,IP), IX=1,NXS), IY=IYB,IYE)                
 1100  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= READ  LIM5   ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.2) THEN                                                      
       DO 1200 IT = 1, NTS                                                      
        DO 1200 IP = -MAXNPS, MAXNPS                                            
         DO 1200 IZ = -1, NIZS                                                  
          DO 1200 IYB = -NY3D, NY3D, JBLOCK                                     
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          READ (NIN) ((LIM5(IX,IY,IZ,IP,IT), IX=1,NXS), IY=IYB,IYE)             
 1200  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= READ  PLRPS  ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1) THEN                                                      
       DO 400 IL = 1, NLS                                                       
        DO 400 IYB = -NYS, NYS, JBLOCK                                          
         IYE = MIN (IYB+JBLOCK-1, NYS)                                          
         READ (NIN) ((PLRPS(IX,IY,IL), IX=1,NXS), IY=IYB,IYE)                   
  400  CONTINUE                                                                 
      ENDIF                                                                     

C                                                                               
C================= READ  PLRP3  ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 1000 IP = -MAXNPS, MAXNPS                                             
        DO 1000 IL = 1, NLS                                                     
         DO 1000 IYB = -NY3D, NY3D, JBLOCK                                      
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          READ (NIN) ((PLRP3(IX,IY,IL,IP), IX=1,NXS), IY=IYB,IYE)               
 1000  CONTINUE                                                                 
      ENDIF                                                                     
C     IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
C      DO 1002 IL = 1,NLS 
C        WRITE(6,*) 'PLAMS:',PLAMS(IL),'PIZS:',PIZS(IL)
C        DO 1002 IX = 15,30                                                     
C          DO 1002 IY = NY3D, -NY3D,-1                                      
C            WRITE(6,'(I7,A,9G12.3)') IY,':',
C    >              (PLRP3(IX,IY,IL,IP), IP=-MAXNPS,MAXNPS)       
C1002  CONTINUE                                                                 
C     ENDIF                                                                     
C                                                                               
C================= READ  POWL3  ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 800 IP = -MAXNPS, MAXNPS                                              
        DO 800 IZ = -1, NIZS                                                    
         DO 800 IYB = -NY3D, NY3D, JBLOCK                                       
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          READ (NIN) ((POWL3(IX,IY,IZ,IP), IX=1,NXS), IY=IYB,IYE)               
  800  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= READ  LINE3  ARRAY FROM DISC ========================        
C                                                                               
      IF (IMODE.NE.1.AND.CDPOL.GT.0.0) THEN                                     
       DO 900 IP = -MAXNPS, MAXNPS                                              
        DO 900 IZ = -1, NIZS                                                    
         DO 900 IYB = -NY3D, NY3D, JBLOCK                                       
          IYE = MIN (IYB+JBLOCK-1, NY3D)                                        
          READ (NIN) ((LINE3(IX,IY,IZ,IP), IX=1,NXS), IY=IYB,IYE)               
  900  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
C================= READ  SDTXS,SDTYS,SDYXS,SDYYS =======================        
C                                                                               
      DO 1350 IZ = 1, NIZS                                                      
        READ (NIN) (SDTXS(IX,IZ),IX=1,NXS),(SDTYS(IY,IZ),IY=-NYS,NYS),          
     >             (SDYXS(IX,IZ),IX=1,NXS),(SDYYS(IY,IZ),IY=-NYS,NYS)           
 1350 CONTINUE                                                                  
C                                                                               
C================= READ  SCTXS,SCTYS ===================================        
C                                                                               
      DO 1360 IZ = 1, NLS                                                      
        READ (NIN) (SCTXS(IX,IZ),IX=1,NXS),(SCTYS(IY,IZ),IY=-NYS,NYS)          
 1360 CONTINUE                                                                  
C                                                                               
C============ READ  ADDITIONAL QUANTITIES - UPDATES, ETC ===============        
C                                                                               
      DO 1500 IQS = -NQXSO, NQXSI, IBLOCK-100                                   
        IQE = MIN (IQS+IBLOCK-100-1, NQXSI)                                     
        READ  (NIN ,IOSTAT=IOS) (QXS(IQX),IQX=IQS,IQE)                          
        READ  (NIN ,IOSTAT=IOS) (QS (IQX),IQX=IQS,IQE)                          
        read  (nin ,iostat=ios) (svybar(iqx),iqx=iqs,iqe)
        read  (nin ,iostat=ios) (svyacc(iqx),iqx=iqs,iqe)
 1500 CONTINUE                                                                  
C                                                                               
      DO 1600 IQS = -NQXSO, 1, (IBLOCK-50)/2                                    
        IQE = MIN (IQS+(IBLOCK-50)/2-1, 0)                                      
        READ  (NIN ,IOSTAT=IOS) ((QEDGES(IQX,J),IQX=IQS,IQE),J=1,2)             
        READ  (NIN ,IOSTAT=IOS) ((QTANS (IQX,J),IQX=IQS,IQE),J=1,2)             
        READ  (NIN ,IOSTAT=IOS) ((QDISTS(IQX,J),IQX=IQS,IQE),J=1,2)             
        READ  (NIN ,IOSTAT=IOS) ((QTEMBS(IQX,J),IQX=IQS,IQE),J=1,2)             
        READ (NIN ,IOSTAT=IOS) ((QTEMBSI(IQX,J),IQX=IQS,IQE),J=1,2)             
        READ  (NIN ,IOSTAT=IOS) ((QRNBS (IQX,J),IQX=IQS,IQE),J=1,2)             
 1600 CONTINUE                                                                  
C                                                                               
      READ  (NIN ,IOSTAT=IOS) (FACTA(IZ),FACTB(IZ),IZ=-1,NIZS)                  
      READ  (NIN ,IOSTAT=IOS) TC,SC,TO,SO,TV,SV,GC,RP                           
C                                                                               
      DO 2000 IYB = -NYS, NYS, JBLOCK                                           
        IYE = MIN (IYB+JBLOCK-1, NYS)                                           
        READ (NIN) ((CTEMBS(IX,IY), IX=1,NXS), IY=IYB,IYE)                      
        READ (NIN) ((CTEMBSI(IX,IY), IX=1,NXS), IY=IYB,IYE)                   
        READ (NIN) ((CRNBS (IX,IY), IX=1,NXS), IY=IYB,IYE)                      
 2000 CONTINUE                                                                  
c
c     Read in particle tracks for debugging - if cstept is greater 
c     than zero - then particle tracks were accumulated.   
c
      write (6,*) 'cstept1:' , cstept
      if (cstept.gt.0) then  
         read(nin) (ptracl(ix),ix = 1,cstept)         
         write(6,*) 'ptracl:',(ptracl(ix),ix=1,cstept)
         read(nin) (((ptracs(ix,iy,iz),ix=1,ptracl(iy))
     >                  ,iy=1,cstept),iz=1,2)
         do 2100 iy = 1,cstept
            write(6,*) 'trajectory:',iy,ptracl(iy)
            do 2100 ix = 1,ptracl(iy)
               write(6,*) ix,':',(ptracs(ix,iy,iz),iz=1,2)
 2100    continue     
      endif  
C                                                                               
C                                                                               
      RETURN                                                                    
 9000 FORMAT(1X,'COLECT: ',A,' VALUES',/,(1X,8F9.4))                            
 9001 FORMAT(1X,'COLECT: ',A,' VALUES',/,(1X,8I9))                              
 9002 FORMAT(1X,'COLECT:   NXS   NYS  NQXSO NQXSI NQYS    NTS  NIZS NLS'        
     >     ,/8X,8I6,/8X,A,/8X,A,/1X,'IMODE',I3,/1X,'ITER',I4,/)                 
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      INTEGER FUNCTION IERR2(STRING, N, MAX)                                    
      CHARACTER STRING*(*)                                                      
      INTEGER   N, MAX                                                          
      WRITE (6,*) ' ERROR: VALUE OF ',STRING,' =',N,' IS GREATER THAN ',        
     > 'PARAMETER VALUE',MAX                                                    
      WRITE (7,*) ' ERROR: VALUE OF ',STRING,' =',N,' IS GREATER THAN ',        
     > 'PARAMETER VALUE',MAX                                                    
      IERR2 = 1                                                                 
      STOP                                                                      
      END                                                                       
