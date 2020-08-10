      PROGRAM FAX                                                       
C                                                                       
C  *********************************************************************
C  *                                                                   *
C  *  FAX: This program takes in the output file from a LIM run, and   *
C  *  rearranges it suitable for photocopying and transmission on      *
C  *  the facsimile machine.  The results take up half as much paper   *
C  *  when processed in this way.                                      *
C  *    The number of lines per page, etc need tuning to the specific  *
C  *  printer being used.  Various minor adjustments may be called for *
C  *  when a new printer becomes available.  Presently, the program is *
C  *  tailored to the IBM3820 High Quality Printer in the DMG at JET.  *
C  *                                                                   *
C  *                                      C.M.Farrell   May 1988       *
C  *                                                                   *
C  *********************************************************************
C                                                                       
      CHARACTER*80 BUFFER(128)                                          
      CHARACTER*34 TITLE                                                
C                                                                       
C---- Set number of lines for first page.  This is the continuation     
C---- point for printing subsequent pages.                              
C                                                                       
      WRITE (6,'('' FAX program started ...'')')                        
      IPAGE = 1                                                         
      NPS = 61                                                          
   50 CONTINUE                                                          
C                                                                       
C---- Read in 2 pages worth of results.  If we are nearing the end of   
C---- file there won't be enough results left - store in mps the actual 
C---- number of lines left.  Although the ERR= and END= parameters are  
C---- coded below the IBM still gives an error message - ignore it !!   
C                                                                       
      MPS = 0                                                           
      DO 100 IP = 1, 2*NPS                                              
        READ (5,'(A80)',ERR=110,END=110) BUFFER(IP)                     
        MPS = MPS + 1                                                   
  100 CONTINUE                                                          
  110 CONTINUE                                                          
      IF (MPS.EQ.0) GOTO 999                                            
C                                                                       
C---- Extract reference for run from 8th line of LIM output.            
C---- Insert page throw and print page no. & title for subsequent pages 
C                                                                       
      IF (IPAGE.EQ.1) THEN                                              
        TITLE = BUFFER(8)(3:36)                                         
        WRITE (7,'(1X)')                                                
      ELSE                                                              
        WRITE (7,'(''1'',/,A,''Page'',I2,/)') TITLE,IPAGE               
      ENDIF                                                             
C                                                                       
C---- Print out 1 page worth of doubled up results, truncating so that  
C---- a total width of 132 characters is obtained.  The final page will 
C---- have some empty space at the end.                                 
C                                                                       
      DO 200 IP = 1, NPS                                                
        IF (IP+NPS.LE.MPS) THEN                                         
          WRITE (7,'(A,A,A)') BUFFER(IP)(1:65),'³',                     
     >                            BUFFER(IP+NPS)(1:65)                  
        ELSEIF (IP.LE.MPS) THEN                                         
          WRITE (7,'( A,A )') BUFFER(IP)(1:65),'³'                      
        ENDIF                                                           
  200 CONTINUE                                                          
      IF (MPS.LT.NPS) GOTO 999                                          
C                                                                       
C---- Leap back to process subsequent pages, with new page length.      
C---- First page generally of different length, whether due to extra    
C---- title being printed, or for a Versatec the initial message, etc.  
C                                                                       
      NPS = 60                                                          
      IPAGE = IPAGE + 1                                                 
      GOTO 50                                                           
  999 WRITE (6,'('' FAX program completed.'')')                         
      STOP                                                              
      END                                                               
