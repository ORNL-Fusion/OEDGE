C                                                                       
C  *********************************************************************
C  *                                                                   *
C  *  SYSIBM: THIS FILE CONTAINS DUMMIES FOR ALL THE EXTERNAL ROUTINES *
C  *  USED THROUGHOUT LIM3.  THE PROGRAM WAS DESIGNED AND OPTIMISED    *
C  *  ON THE IBM3090E AT JET AND MAKES USE OF SOME IBM UTILITIES,      *
C  *  THE HARWELL SUBROUTINE LIBRARY ("HSL") AND THE NOCORONA PACKAGE  *
C  *  ALSO DEVELOPED AT JET.                                           *
C  *    THE CODE COLLECTED HERE WILL REQUIRE MODIFICATIONS WHEN        *
C  *  TRANSPORTING BETWEEN MACHINES, WHILST THE REMAINDER OF THE LIM3  *
C  *  CODE SHOULD REMAIN INTACT.                                       *
C  *                                                                   *
C  *  ROUTINES IN THIS FILE :-                                         *
C  *    XUFLOW: PREVENT UNDERFLOW INTERRUPTS                           *
C  *    SURAND: GET VECTOR OF RANDOM NUMBERS IN (0,1)                  *
C  *    ZA08AS: GET TIME OF DAY AS 8 CHARACTER STRING                  *
C  *    ZA09AS: GET DATE AS 8 CHARACTER STRING                         *
C  *    ZV01AD: GET NAME OF DATAFILE CONNECTED TO CHANNEL 5            *
c  *            - replaced by ioname in rundiv, it is called by        *
c  *              only in the ibm3090 implementation.                  *
C  *    ZA02AS: GET TIME USED SO FAR IN SECONDS                        *
C  *    RANINI: DUMMY INTERFACE TO RANSET RANDOM NUMBER INITIALISER    *
C  *    INVPIN: SPAWN OR CALL PIN AS A ROUTINE                         *
C  *    SURAND2: RETURNS A SINGLE RANDOM NUMBER                        *
c  *    Printerinit: Initialize printer calls for GHOST                *
c  *    killdiv: This subroutine calls an entry point that will allow  *
c  *             cleaner exits on a trapped USR1 signal to the process *
c  *             on a Unix Workstation                                 *
c  *    initkill: This subroutine performs the initialization required *
c  *              for the above subroutine.                            *
c  *    initnc:  Initialize the NOCORONA package.                      *
c  *    ncrrates:Reaction rates fron Nocorona                          *    
c  *    ncrdlong:Radiative data from Nocorona.                         *  
c  *    ioname  : This routine returns the name of a file connected    *
c  *              to a specified i/o stream.                           *
c  *    DMGUID  :Routine to get userid for person executing the code   *
C  *                                                                   *
C  *                                      C.M.FARRELL   JUNE 1988      *
C  *                                      L.D.HORTON    APRIL 1993     *
C  *                                                                   *
C  *********************************************************************
C                                                                       
C     XUFLOW                                                            
C     ======                                                            
C     IBM  : SYSTEM ROUTINE TO PREVENT UNDERFLOW INTERRUPTS OCCURING    
C     CRAY : REPLACE WITH DUMMY ROUTINE HERE.                           
C                                                                       
C     SUBROUTINE XUFLOW (IFLAG)                                         
C     INTEGER IFLAG                                                     
C     WRITE (6,'('' XUFLOW: DUMMIED OUT FOR THIS APPLICATION.'')')      
C     RETURN                                                            
C     END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     SURAND                                                            
C     ======                                                            
C     IBM  : ESSL LIBRARY ROUTINE TO GENERATE VECTOR OF RANDOM NUMBERS  
C     CRAY : REPLACE WITH CALLS TO RANF GENERATOR WITHIN A VECTORISABLE 
C            LOOP.                                                      
C                                                                       
C     SUBROUTINE SURAND (SEED,NRANDS,RANDS)                             
C     DOUBLE PRECISION SEED                                             
C     REAL RANDS(NRANDS)                                                
C     DO 100 J = 1, NRANDS                                              
C       RANDS(J) = RANF ()                                              
C 100 CONTINUE                                                          
C     RETURN                                                            
C     END                                                               
C                                                                       
C                                                                       
C***********************************************************************
C                                                                       
C     SURAND2                                                           
C     ======                                                            
C     IBM  : ESSL LIBRARY ROUTINE TO GENERATE VECTOR OF RANDOM NUMBERS  
C     CRAY : REPLACE WITH CALLS TO RANF GENERATOR WITHIN A VECTORISABLE 
C            LOOP.                                                      
C                                                                       
      SUBROUTINE SURAND2 (SEED,NRANDS,RANDS)                            
      DOUBLE PRECISION SEED                                             
      REAL RANDS(NRANDS)                                                
c                                                                       
      call surand(seed,nrands,rands)                                    
c                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     ZA08AS                                                            
C     ======                                                            
C     IBM  : HARWELL LIBRARY ROUTINE TO EXTRACT TIME IN 8 CHARACTERS    
C     CRAY : REPLACE WITH CALL TO CLOCK SYSTEM ROUTINE.                 
C                                                                       
C     SUBROUTINE ZA08AS (SYSTIM)                                        
C     CHARACTER*8 SYSTIM                                                
C     CALL CLOCK (SYSTIM)                                               
C     RETURN                                                            
C     END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     ZA09AS                                                            
C     ======                                                            
C     IBM  : HARWELL LIBRARY ROUTINE TO EXTRACT DATE IN 8 CHARACTERS    
C     CRAY : REPLACE WITH CALL TO DATE SYSTEM ROUTINE.                  
C                                                                       
C     SUBROUTINE ZA09AS (SYSDAT)                                        
C     CHARACTER*8 SYSDAT                                                
C     CALL DATE (SYSDAT)                                                
C     RETURN                                                            
C     END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     ZV01AD                                                            
C     ======                                                            
C     IBM  : HARWELL LIBRARY ROUTINE TO EXTRACT DATASET NAME CONNECTED  
C            TO CHANNEL 5, AND TO PUT THE NAME IN JFCB(1:44) AND THE    
C            MEMBER IN JFCB(45:52).                                     
C     CRAY : NO EQUIVALENT SYSTEM ROUTINE.                              
C                                                                       
C     SUBROUTINE ZV01AD (IUNIT, VSN, DSN, JFCB)                         
C     CHARACTER VSN*8,DSN(3)*8,JFCB*176                                 
C     WRITE (6,'('' ZV01AD: DUMMIED OUT FOR THIS APPLICATION.'')')      
C     RETURN                                                            
C     END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     ZA02AS                                                            
C     ======                                                            
C     IBM  : HARWELL LIBRARY FUNCTION TO EXTRACT CPU TIME USED SO FAR.  
C     HOT  : FOR HOTSPOT ANALYSIS, DUMMY OUT BY SETTING ZA02AS = 0.0    
C     CRAY : REPLACE WITH SYSTEM FUNCTION SECOND.                       
C                                                                       
C     REAL FUNCTION ZA02AS (IFLAG)                                      
CHOT  ZA02AS = 0.0                                                      
CRAY  ZA02AS = SECOND ()                                                
C     RETURN                                                            
C     END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     RANINI                                                            
C     ======                                                            
C     IBM  : DUMMY ROUTINE - NO NEED TO CALL RANSET                     
C     CRAY : INTERFACE TO RANDOM NO. INITIALISER SYSTEM ROUTINE RANSET  
C                                                                       
      SUBROUTINE RANINI (ISEED)                                         
      INTEGER ISEED                                                     
CRAY  CALL RANSET (ISEED)                                               
      RETURN                                                            
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
C     INVPIN                                                            
C     ======                                                            
C     IBM  : CALL PIN AS A SUBROUTINE                                   
C     CRAY : CALL PIN AS A SUBROUTINE                                   
C     UNIX : START PIN AS AN INDEPENDENT PROCESS                        
C                                                                       
      SUBROUTINE INVOKEPIN(ACTPIN)                                      
      CHARACTER*(*) ACTPIN                                              
C                                                                       
C     DAVID ELDER, FEB 4 1993                                           
C                                                                       
      REAL NIMTIM,ZA02AS                                                
      EXTERNAL ZA02AS                                                   
C                                                                       
CUNIX CALL SYSTEM(ACTPIN)                                               
C                                                                       
      NIMTIM = ZA02AS(1)                                                
      CALL PINPGX                                                       
      NIMTIM = ZA02AS(1) - NIMTIM                                       
      WRITE(6,*) 'TIME USED IN NIMBUS:',NIMTIM, ' (S) '                 
C                                                                       
      RETURN                                                            
      END                                                               
c                                                                       
c                                                                       
c                                                                       
      subroutine printerinit                                            
c                                                                       
c     Sends site dependent GHOST commands to the printer                
c                                                                       
c     These commands are not sent on the JET 3090 but are sent          
c     on the RS6000 in Toronto                                          
c                                                                       
c      call hrdlin(1)                                                   
c      call hrdchr(1)                                                   
c                                                                       
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c                                                                       
      subroutine killdiv                                                
c                                                                       
c     This routine calls the entry point in the div.d3a                 
c     to wrap up processing and then exit.                              
c                                                                       
c     This is not supported on the IBM mainframe - it is only           
c     supported under UNIX.                                             
c                                                                       
c     These are subroutine stubs on the IBM mainframe                   
c                                                                       
c                                                                       
c      call divkill                                                     
c                                                                       
      stop                                                              
      end                                                               
c                                                                       
c                                                                       
c                                                                       
      subroutine initkill                                               
c                                                                       
c     DEFINE the SIGUSR1 kill signal so that the                        
c     signal call can trap it - if it is sent                           
c     to the DIVIMP process. (Workstation only)                         
c                                                                       
c     These are subroutine stubs on the IBM mainframe                   
c                                                                       
c      integer SIGUSR1                                                  
c      parameter (SIGUSR1=30)                                           
c      external killdiv                                                 
c                                                                       
C     SIGUSR1 - comment out for mainframe - or move to routine in the   
c               system module.                                          
C                                                                       
c      call signal(SIGUSR1,killdiv)                                     
c                                                                       
      return                                                            
      end                                                               
c
c
c
      subroutine initnc(crmi)
      implicit none
      real crmi
c
c     This subroutine calls the Nocorona initialization subroutine.
c     It has been separated into the system module in order to allow
c     NOCORONA to be excluded easily .. these can be converted into
c     null subroutines by commenting out the NOCORONA calls.
c                                                                        
      real pmass(1) 
      integer kfail(1)       
c                                                                        
      PMASS(1) = CRMI                                                 
      CALL AATOM (PMASS, 1, 9, KFAIL)                                 
      IF (KFAIL(1).NE.0) THEN                                         
         WRITE (6,9100) CRMI                                           
         STOP                                                          
      ENDIF                                                           
c
 9100 FORMAT(//1X,'ELEMENT WITH MASS ',G11.4,' IS NOT INCLUDED',        
     >    ' IN THE NOCORONA PACKAGE.',/)                                
      return
      end                                                                    
c
c
c
      subroutine ncrrates (nksir)
      implicit none
      integer nksir
      include 'params' 
      include 'cnoco' 
c
c     This subroutine calls the RRATES subroutine in the Nocorona 
c     package. It has been placed in the system module so that
c     dependence on Nocorona can be easily removed by commenting out the
c     call here.
c   
      CALL RRATES (PTES, PNES, PRATES, NKSIR)
c
      return
      end
c
c
c
      subroutine ncrdlong(nksir)
      implicit none
      integer nksir
      include 'params'
      include 'cnoco'         
c
c     This subroutine calls the RDLONG subroutine in the Nocorona 
c     package. It has been placed in the system module so that
c     dependence on Nocorona can be easily removed by commenting out the
c     call here.
c
      CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,nksir)              
c
      return  
      end
c
c
c
      subroutine ioname(iunit,filename,itemp,ierr)
      implicit none
      integer iunit, itemp, ierr
      character filename*(*)
c
      character vsn*8, dsn(3)*8, jfcb*176
c
      ierr = 0
c
      call zv01ad (iunit, vsn, dsn, jfcb)                       
      if (jfcb(45:52).ne.' ') then                                      
        do 220 nc = 1, 40                                               
          if (jfcb(nc:nc).EQ.' ') then                                  
            jfcb(nc:nc) = '('                                           
            do 210 nm = 1, 8                                            
              if (jfcb(44+nm:44+nm).NE.' ') then                        
                jfcb(nc+nm:nc+nm) = jfcb(44+nm:44+nm)                   
              else                                                      
                jfcb(nc+nm:nc+nm) = ')'                                 
                goto 230                                                
              endif                                                     
  210       continue                                                    
            jfcb(nc+9:nc+9) = ')'                                       
            goto 230                                                    
          endif                                                         
  220   continue                                                   
  230   continue                                                        
      endif                                                             
c
      filename = jfcb
      return
      end
c
C================================================================
c
      SUBROUTINE DMGUID(SYSUID,PREFIX)
C
C RETURNS USERID
C
      CHARACTER*(*) SYSUID, PREFIX
C
      CALL GETENV('LOGNAME',SYSUID)
C
      PREFIX = ' '
C
      RETURN
      END
c
C================================================================
c
                                                                        
