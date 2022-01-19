      SUBROUTINE PRINT (AS,BS,MAXNAS,MAXNBS,NAS,NBS,                            
     >                  ALAB,BLAB,BLABS,REF,IPR)                                
      IMPLICIT  none
      INTEGER   MAXNAS,MAXNBS,NAS,NBS,IPR                                       
      REAL      AS(*),BS(MAXNAS,*)                                              
      CHARACTER ALAB*24,BLAB*24,BLABS(*)*36,REF*36                              
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  PRINT:   THIS ROUTINE PRINTS VALUES OF 2D ARRAY BS AGAINST AS    *        
C  *                                                                   *        
C  *  C.M.FARRELL   OCT 1987                                           *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER IA,IB,JB                                                          
C                                                                               
      IF (IPR.NE.2) RETURN                                                      
      CALL PRC (REF)                                                            
      CALL PRB                                                                  
C                                                                               
      DO 110 IB = 1, NBS, 4                                                     
        WRITE (7,'(2X,A16,'' : '',13X,A16)') ALAB,BLAB                          
        WRITE (7,'(2X,16X,'' :   '',4(A12,2X))')                                
     >    (BLABS(JB), JB = IB, MIN (NBS,IB+3))                                  
        WRITE (7,'(2X,16(''-''),''-:-'',56(''-''))')                            
        DO 100 IA = 1, NAS                                                      
          WRITE (7,'(2X,F12.4,4X,'' :   '',1P,4(G12.5,2X))')                    
     >      AS(IA),(BS(IA,JB), JB = IB, MIN (NBS,IB+3))                         
  100   CONTINUE                                                                
        CALL PRB                                                                
  110 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
