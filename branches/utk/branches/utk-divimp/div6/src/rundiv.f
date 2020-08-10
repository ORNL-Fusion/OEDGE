c     -*Fortran*-
c
      PROGRAM RUNDIV                                                    
      use debug_options
      use rand_data
c slmod begin
      use mod_divimp
      use ero_interface
c slmod end
      use mod_fp_data
      IMPLICIT NONE
C                                                                       
C  *********************************************************************
C  *                                                                   *
C  *  RUNDIV:  MONTE CARLO IMPURITY TRANSPORT CODE APPLIED TO THE      *
C  *           JET DIVERTOR.                                           *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *                                                                   *
C  *********************************************************************
C                                                                       
      include 'params'                                                  
      include 'cadas'                                                   
      include 'comtor'                                                  
      include 'cgeom'                                                   
      include 'dynam4'                                                  
      include 'grbound'                                                 
C                                                                       
      INTEGER        IERR,NM,NC,ICHAR,IZ,NYMFS,J,ierr2,in                        
      INTEGER        NIZS,KFAIL(1),NITERS                               
      INTEGER        NIMPS,NIMPS2,ITER,IT                               
      INTEGER        KNEUTA,KNEUTB,KNEUTC,KNEUTE,NRAND                  
      REAL           IONTIM,NEUTIM,STATIM,TOTTIM,ZA02AS,PMASS(1)        
      REAL           CPULIM,RAN                                         
      REAL           FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)                  
c
c     Some local variables for printing
c
      logical chars_left
      integer start_char, end_char, last_end

c
C                                                                       
c     The random number seed generator
c
      integer        calc_random_seed,len,lenstr
      external       calc_random_seed,lenstr
c
c     NOTE: These string variables are written in the RAW file - as
c           a result - any change in their size MUST be accompanied
c           by matching code in the GET routine in IOOUT.D6A
c      
      CHARACTER      TITLE*174,desc*1024,desc_line*58
      character      JOB*72,JFCB*176,COMENT*77,NUMBER(20)*4    
      CHARACTER      EQUIL*60                                           
      CHARACTER*8    SYSTIM,SYSDAT,VSN,DSN(3)                           
      character*5    verse
      DOUBLE PRECISION SEED                                             
      DATA    NUMBER /' 1ST',' 2ND',' 3RD',' 4TH',' 5TH',' 6TH',' 7TH', 
     >  ' 8TH',' 9TH','10TH','11TH','12TH','13TH','14TH','15TH','16TH', 
     >  '17TH','18TH','19TH','20TH'/                                    
c
c     Set hard-coded global trace debugging options
c      call init_trace(0,.true.)
      call init_trace(0,.false.)
      call pr_trace('RUNDIV','BEGIN EXECUTION')
c
c     Iniialize the main .dat file output unit number
c
      datunit = 7
c
c     Initialize some string variables
c
      write (job,'(72x)')
      write (jfcb,'(176x)')
c
c     Initialize the global random number count
c     Found in module rand_data
c
      ran_used = 0.0
c
c     Extract version and revison numbers from the verson string set
c     in the params common block.
c
      verse = verson  
c
      read(verse,'(i1,2x,i2)') vernum,revnum
c
      write (6,*) 'DIVIMP VERSION = ', vernum, ' REVISION = ',revnum
c
      call pr_trace('RUNDIV','VERSION ASSIGNED')
c
C                                                                       
C-----------------------------------------------------------------------
C     READ IN DATA - MAKE TEMP COPIES OF SOME INPUT FLAGS OVERWRITTEN   
C     IN SELF-SPUTTERING CASES, SINCE THEY MAY BE REQUIRED FOR          
C     SUBSEQUENT ITERATIONS.                                            
C-----------------------------------------------------------------------
C                                                                       
      CALL XUFLOW (0)                                                   
      STATIM = ZA02AS (1)                                               
      IONTIM = 0.0                                                      
      NEUTIM = 0.0                                                      
      outgrid = .false.                                                 
c                                                                       
c     This call initializes the trapping of the signal USR1 - on        
c     a UNIX workstation.                                               
c                                                                       
      call initkill                                                     
c
      call pr_trace('RUDNIV','INITKILL RUN')

C                                                                       
c     REPLACED BY NEUTRAL WALL OPTION
c
C     The following variable is used to differentiate between           
c     Wall options 2 and 4, and 3 and 5. If wall option 4 is            
c     selected ... WALLSWCH is set to true and the wall option          
c     is changed to 2. Similarly 5 is set to 3. Wall options 4          
c     and 5 are ONLY used to specifiy where the wall coordinates        
c     are obtained and not how they are treated after loading.          
c                                                                       
c      wallswch = .false.                                                
C                                                                       
      IERR = 0                                                          
      CALL READIN (TITLE,desc,equil,NIZS,NIMPS,NIMPS2,CPULIM,IERR,
     >             NYMFS,NITERS)    
c
      call pr_trace('RUNDIV','AFTER READIN')
c
c      IF (IERR.NE.0) GOTO 1002                                          
c
c      CALL READIN (TITLE,NIZS,NIMPS,NIMPS2,CPULIM,IERR,NYMFS,NITERS)    
C                                                                       
C     PUT pre-loaded geometry data (neutral wall and target coordinates)
c     into appropriate variables if a geometry (one of several specific 
c     shots) has been specified.                                        
c                                                                       
      IF (CGEOOPT.NE.-1) THEN                                           
         CALL LOADGEO                                                   
      ENDIF                                                             
C                                                                       
      KNEUTA = CNEUTA                                                   
      KNEUTB = CNEUTB                                                   
      KNEUTC = CNEUTC                                                   
      KNEUTE = CNEUTE                                                   
C                                                                       
C-----------------------------------------------------------------------
C     SET UP TIME POINTS FOR RECORDING ION DISTRIBUTIONS IF IMPULSE     
C     MODE REQUIRED.  ARRAY "CTIMES" WILL HOLD THE "ITERATION NOS AT    
C     WHICH THE POSITION IS TO BE RECORDED" FOR EACH CHARGE STATE.      
C     ARRAY "TIMES" WILL HOLD THE ACTUAL TIMEPOINTS (IN SECONDS).       
C     PRIMARY NEUTRALS WILL BE TREATED AS FOR TOTAL NEUTRALS.           
C     SET "CSTMAX" TO ITERATION CUTOFF POINT, EITHER THE VALUE OF       
C     CTIMMAX (USUALLY 10 SECONDS) FOR THE STEADY STATE CASE            
C     OR MAX OF GIVEN DWELL TIME * FACTORS                              
C-----------------------------------------------------------------------
C                                                                       
      call rzero(ctimes,(maxnts+2)*(maxizs+2)) 
      CSTMAX = 0                                                        
      IF (IMODE.NE.2) THEN                                              
        DO 190 IT = 1, NTS                                              
          CTIMES(IT,-1) = DWELTS(-1)*DWELFS(IT)/FSRATE                  
          CTIMES(IT, 0) = DWELTS( 0)*DWELFS(IT)/FSRATE                  
          CSTMAX = MAX (CTIMES(IT,-1),CTIMES(IT,0)) * FSRATE / QTIM     
          IF (NIZS.GT.0) THEN                                           
            DO 185 IZ = 1, NIZS                                         
              CTIMES(IT,IZ) = DWELTS(IZ)*DWELFS(IT)/QTIM                
              CSTMAX = MAX (CSTMAX, CTIMES(IT,IZ))                      
  185       CONTINUE                                                    
          ENDIF                                                         
  190   CONTINUE                                                        
c
c       Set a zero'th time point for calculating time-bin widths easily 
c
        do iz = -1,nizs
c
           ctimes(0,iz) = 0.0
c
        end do 
c
      ENDIF                                                             
C                                                                       
C---- TAKE MAX OF THIS VALUE WITH CTIMMAX SECONDS, IF STEADY STATE REQUI
C                                                                       
      IF (IMODE.NE.1) CSTMAX = MAX (CSTMAX, CTIMMAX/QTIM)               
C                                                                       
C---- INTRODUCE AN EXTRA TIMEPOINT WHICH WILL NEVER BE REACHED          
C---- THIS TIMEPOINT WILL APPLY EQUALLY TO NON-IMPULSE MODE CASES       
C                                                                       
      CTIMES(NTS+1,-1) = 2.0 * CSTMAX * QTIM / FSRATE                   
      CTIMES(NTS+1, 0) = 2.0 * CSTMAX * QTIM / FSRATE                   
      IF (NIZS.GT.0) THEN                                               
        DO 195 IZ = 1, NIZS                                             
          CTIMES(NTS+1,IZ) = 2.0 * CSTMAX                               
  195   CONTINUE                                                        
      ENDIF                                                             
C                                                                       
C-----------------------------------------------------------------------
C     SET JOB TO INCLUDE DATAFILE,TIME AND DATE TO PROVIDE UNIQUE       
C     REFERENCE FOR THIS RUN.                                           
C     THE REFERENCE WILL BE PRINTED AS TWO CHAR STRINGS OF LENGTH 36.   
C-----------------------------------------------------------------------
C                                                                       
      CALL ZA08AS (SYSTIM)                                              
      CALL ZA09AS (SYSDAT)                                              
c
c  introduce a system-specific file to get the name of file
c  connected to unit 5
c
c      call ioname(5,jfcb,99,ierr2)
c
      call ioname(jfcb,ierr2)
      if (ierr2.ne.0) then
        write(6,*) ' No file connected to unit 5?!'
        stop
      endif
c
c      CALL ZV01AD (5, VSN, DSN, JFCB)                                   
c      IF (JFCB(45:52).NE.' ') THEN                                      
c        DO 220 NC = 1, 40                                               
c          IF (JFCB(NC:NC).EQ.' ') THEN                                  
c            JFCB(NC:NC) = '('                                           
c            DO 210 NM = 1, 8                                            
c              IF (JFCB(44+NM:44+NM).NE.' ') THEN                        
c                JFCB(NC+NM:NC+NM) = JFCB(44+NM:44+NM)                   
c              ELSE                                                      
c                JFCB(NC+NM:NC+NM) = ')'                                 
c                GOTO 230                                                
c              ENDIF                                                     
c  210       CONTINUE                                                    
c            JFCB(NC+9:NC+9) = ')'                                       
c            GOTO 230                                                    
c          ENDIF                                                         
c  220   CONTINUE                                                        
c  230   CONTINUE                                                        
c      ENDIF                                                             
c
      len = lenstr(jfcb)
c
c     WRITE (JOB,'(A36,1x,A8,1X,A8,'' DIVIMP'',A5)')                       
c     >  JFCB(max(1,len-35):len),SYSDAT,SYSTIM,VERSON
c
      WRITE (JOB,'(36X,1x,A8,1X,A8,'' DIVIMP'',A5)')  
     >  SYSDAT,SYSTIM,VERSON

c
c      WRITE (JOB,'(''VERSION: '',A8,1X,A8,'' DIVIMP'',A5)')
c     >   SYSDAT,SYSTIM,VERSON                                 
C                                                                       
C-----------------------------------------------------------------------
C                     PRINT HEADING                                     
C-----------------------------------------------------------------------
C                                                                       
      WRITE (COMENT,'(''*'',60('' ''),''*'')')                          
      CALL PRC                                                          
     >('**************************************************************')
      CALL PRC (COMENT)                                                 
      WRITE(7,'('' *'',17X,''RUN OF DIVIMP VERSION '',A5,16X,''*'')')   
     >  VERSON                                                          
      WRITE(7,'('' *'',17X,27(''-''),16X,''*'')')                       
      CALL PRC (COMENT)                                                 
c
c     Write out title
c
c
      len = lenstr(title)
      do in = 1,len,58
         WRITE (7,'(1X,''* '',A58,'' *'')') 
     >                          title(in:in+57)                    
c
c     >                          title(in:min(in+57,len))                    
c
      end do  
c
c      WRITE (7,'(1X,''* '',A58,'' *'')') TITLE(1:58)                    
c
c      len = lenstr(title)
c
c      if (len.gt.52) then 
c          WRITE (7,'(1X,''* '',A58,'' *'')') TITLE(59:116)
c      endif
c      if (len.gt.104) then 
c          WRITE (7,'(1X,''* '',A58,'' *'')') TITLE(117:174)
c      endif
c
c     Write out DESCRIPTION
c
      len = lenstr(desc)
c
      if (len.gt.0) then 
         chars_left = .true.
      else
         chars_left = .false.
      endif
c
      start_char= 1
      end_char = 1 
      last_end = 1
c
      do while (chars_left) 
c
         if (end_char.eq.len) then 
c
            desc_line = desc(start_char:end_char)
c
            WRITE (7,'(1X,''* '',A58,'' *'')') 
     >             desc_line                    
            chars_left = .false.        
c
         elseif ((end_char-start_char).gt.57)  then 
c
            if (last_end.eq.start_char) then  
c
               desc_line = desc(start_char:end_char)
c
               WRITE (7,'(1X,''* '',A58,'' *'')') 
     >              desc_line
c
               start_char = end_char+1
               last_end   = start_char 
            else
c
               desc_line = desc(start_char:last_end)
c
               WRITE (7,'(1X,''* '',A58,'' *'')') 
     >              desc_line              
c
               start_char = last_end + 1
               last_end   = start_char
            endif
c
         elseif (desc(end_char:end_char).eq.' ') then 
            last_end = end_char
            end_char = end_char+1 
         else
            end_char = end_char + 1
         endif             
c
      end do

c
c      do in = 1,len,58
c         WRITE (7,'(1X,''* '',A58,'' *'')') 
c     >                          desc(in:in+57)                    
c
c     >                          desc(in:min(in+57,len))                    
c
c      end do  
c
c     Other reference strings 
c
      CALL PRC (COMENT)                                                 
      WRITE (7,'(1X,''* '',A58,'' *'')') JFCB(1:58)                   
      CALL PRC (COMENT)                                                 
      WRITE (7,'(1X,''* '',A53,5x,'' *'')') JOB(37:72)                   
      CALL PRC (COMENT)                                                 
      WRITE (7,'(1X,''* '',A58,'' *'')') EQUIL(1:58)                 
      CALL PRC (COMENT)                                                 
      CALL PRC                                                          
     >('**************************************************************')
      CALL PRB                                                          
      WRITE (6,'(/1X,A,//2X,A,//2x,a,/)') TITLE,JFCB(1:53),JOB                            
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
         ciseed = calc_random_seed(-1)
      ENDIF                                                             

      SEED = DBLE (REAL (CISEED))                                       



C                                                                       
C---- SET RANDOM NUMBER SEED VIA RANINI DUMMY ROUTINE IF REQUIRED       
C---- FOR NON-IBM APPLICATIONS.                                         
C---- GET FIRST RANDOM NUMBER IN SEQUENCE (AND DISCARD IT).             
C                                                                       
      CALL RANINI (CISEED)                                              
      CALL SURAND2 (SEED, 1, RAN)                                       
      NRAND = 1                                                         
C                                                                       
C-----------------------------------------------------------------------
C     INITIALISE NOCORONA PACKAGE WITH 1 IMPURITY & OUTPUT TO CHANNEL 9 
C-----------------------------------------------------------------------
C                                                                       
      if (cdatopt.eq.0) then                                            
c
c       Initialize the No corona package
c
        call initnc(crmi)
c
      elseif (cdatopt.eq.2) then 
c 
c       Initialize the ADPAK database if required.
c
        call readmc 
c
      elseif (cdatopt.eq.3) then 
c 
c       Initialize the INEL database if required.
c
        call inelinput(cion)
c
      endif                                                             
c
c      PMASS(1) = CRMI                                                 
c      CALL AATOM (PMASS, 1, 9, KFAIL)                                 
c      IF (KFAIL(1).NE.0) THEN                                         
c         WRITE (6,9100) CRMI                                           
c         STOP                                                          
c     ENDIF                                                           
c
c 9100 FORMAT(//1X,'ELEMENT WITH MASS ',G11.4,' IS NOT INCLUDED',        
c     >    ' IN THE NOCORONA PACKAGE.',/)                                
c
      call pr_trace('RUNDIV','BEFORE DIV')
C                                                                       
C-----------------------------------------------------------------------
C  CALL DIV TO CALCULATE IMPURITY LEVELS                                
C-----------------------------------------------------------------------
C                                                                       
      ITER = 1                                                          
      REWIND (8)                                                        
  500 CALL DIV (title,equil,
     >          NIZS,NIMPS,NIMPS2,CPULIM,IONTIM,NEUTIM,
     >          SEED,NYMFS,FACTA,FACTB,ITER,NRAND)                           


      call pr_trace('RUNDIV','AFTER DIV')
c
C-----------------------------------------------------------------------
C   DUMP RESULTS IN THE JET TRAN FILE                                    
C-----------------------------------------------------------------------
C                                                                       
c     Print out the TRAN file for JET post-processors if the 
c     option is set. 
c
 
      if (write_tran.eq.1) then 
         call divtrn(nizs,iter,niters,facta,factb,title,job,equil,desc,
     >               jfcb)
      endif

      call pr_trace('RUNDIV','AFTER DIVTRN')
c
C-----------------------------------------------------------------------
C   DUMP RESULTS IN AN EXTERNAL FILE                                    
C-----------------------------------------------------------------------
C                                                                       
      CALL STORE (TITLE,desc,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)

      call pr_trace('RUNDIV','AFTER STORE')
C                                                                       
C-----------------------------------------------------------------------
C  CHECK FOR FURTHER ITERATIONS FOR SELF-CONSISTENT PLASMA              
C-----------------------------------------------------------------------
C                                                                       
      IF (ITER.LT.NITERS .AND. ABSFAC.GT.0.0) THEN                      
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
        CALL TAUIN3 (NIZS,ABSFAC)                                       
        GOTO 500                                                        
      ENDIF                                                             
C                                                                       
C-----------------------------------------------------------------------
c
c     The following subroutine is a "dummy" for now ... it copies
c     the contents of certain arrays into a common block that will
c     be shareable between B2, Eirene and DIVIMP
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
c      call wrtdivbra(nizs)
c
C-----------------------------------------------------------------------
c
      call prb
      call prchtml('CASE EPILOGUE','pr_end','0','B')
      call prb   
      TOTTIM = ZA02AS (1) - STATIM                                      
      WRITE (datunit,'('' TOTAL RANDOM NUMBERS USED'',I12)') NRAND            
      call pri ('TIME SPENT IN PIN/NIMBUS(S)   ',NINT(TOTPINTIM))
      CALL PRI ('TIME FOLLOWING NEUTRALS (S)   ',NINT(NEUTIM))          
      CALL PRI ('TIME FOLLOWING IONS     (S)   ',NINT(IONTIM))          
      CALL PRI ('TOTAL CPU TIME USED     (S)   ',NINT(TOTTIM))          
      CALL PRB                                                          
      WRITE (6,'('' TOTAL RANDOM NUMBERS USED'',I12)') NRAND            
      write (6,'('' TIME SPENT IN PIN/NIMBUS(S)   '',g11.4)') TOTPINTIM
      WRITE (6,'('' TIME FOLLOWING NEUTRALS (S)   '',G11.4)') NEUTIM    
      WRITE (6,'('' TIME FOLLOWING IONS     (S)   '',G11.4)') IONTIM    
      WRITE (6,'('' TOTAL CPU TIME USED     (S)   '',G11.4)') TOTTIM    

      call pr_trace('RUNDIV','AFTER EPILOGUE')

      call create_html(jfcb)

      call pr_trace('RUNDIV','AFTER HTML')

      CALL OutputData  (87,'END OF DIVIMP')
      CALL OutputEIRENE(65,'END OF DIVIMP')

      CALL divClean

      ! Clean up ERO related data and make sure ERO related output is written
      if (ero_opt.ne.0) then 
         CALL ero_cleanup
      endif
 
      ! Put storage de-allocation calls here - just cleanup since end 
      ! of execution should get rid of them anyway
      call fp_deallocate_storage


      STOP 'END OF DIVIMP: NORMAL EXECUTION COMPLETE'

C                                                                      
 1002 CALL PRC ('MAIN: ERROR OCCURED DURING DATA INPUT - RUN ABORTED')      
      STOP 'ERROR OCCURED DURING DATA INPUT - RUN ABORTED'                                                            
      END                                                               
c
c
c
      subroutine wrtdivbra(nizs)
      implicit none
c
      integer nizs
c
c     This subroutine fills the common block that will be used by
c     B2-Eirene - when the codes are called together. 
c
c     It should not be included in the rundiv module but is placed here 
c     for now until a suitable source code module for B2-Eirene 
c     dependent code is implemented.
c
c     David Elder, Sept 6, 1994
c
      include 'params'
      include 'dynam1'
      include 'dynam3' 
      include 'cgeom'
c
      include 'divbra'
c
c     Local variables
c
      integer ik,ir,iz
c
c     Test validity of DIVBRA parameters.
c    
      if (dbmaxnrs.lt.nrs.or.dbmaxnks.lt.nks(irsep).or.
     >    dbmaxizs.lt.nizs) then 
         write(6,*) 'ERROR transcribing to DIVBRA common block'
         write(6,*) 'The dimensions of the B2-Eirene variables'
         write(6,*) 'are too small to conatin the DIVIMP data.'
         write(6,*) 'DB array sizes: ',dbmaxnks,dbmaxnrs,dbmaxizs
         return
      endif
c
c     Transcribe values into the B2-Eirene transfer common block. 
c
      do iz = -1,nizs
         do ir = 1,nrs
            do ik = 1,nks(ir)
               dbni(ik,ir,iz) = ddlims(ik,ir,iz)
               dbti(ik,ir,iz) = ddts(ik,ir,iz)
               dbv(ik,ir,iz)  = dble(sdvs(ik,ir,iz))
               dbpowls(ik,ir,iz) = dble(powls(ik,ir,iz))
               dblines(ik,ir,iz) = dble(lines(ik,ir,iz))
            end do
         end do
      end do
c
c     End of routine
c      
      return
      end
c
c
c
      subroutine create_html(casename)
      implicit none
      include 'params'
      character*(*) casename  
c
c     CREATE_HTML: The purpose of this routine is 
c                  to merge the actual .dat file with the
c                  hypertext references, links and formatting
c                  that has been saved for some lines in the 
c                  code as it went along. The end result will
c                  be an integrated version of the .dat file
c                  containing a selection of hypertext links 
c                  as a table of contents.
c        
      integer len,lenstr,newlen,in
      external lenstr 
      character*512 htmlline
      character*512 datline
      character*530 newline
c
c
c
      call openhtml(casename)
c
c     The .dat file has been tagged so that every html line starts
c     with ()  - this will identify the html lines that must be 
c     substituted in order - the rest of the .dat file is simply
c     copied with the addition of some html.
c
c         
      rewind(tmpunit)
      rewind(datunit)
c
c     Start reading and echoing the data file until the text
c     matches that from the htmlline. 
c
    
 10   read(datunit,'(A512)',end=100,err=200) datline
    
c
      len = lenstr(datline)
c
c     Found HTML reference line - add next from HTML file
c
c     Initialize newline
c
      newline = ' '
c
      if (datline(1:3).eq.'(*)') then 
c
         read(tmpunit,'(A512)',end=200,err=200) htmlline
c 
         len = lenstr(htmlline) 
         write(htmlunit,*) htmlline(1:len)
c
      else
c
         len = lenstr(datline)
c   
c         newline = '<PRE>'//datline(1:len)//'</PRE>'
c
         newlen = 0
c
         do in = 1,len
c
            if (datline(in:in).eq.'>') then 
c
               newline = newline(1:newlen)//'&gt'
               newlen = lenstr(newline)           
c
            elseif (datline(in:in).eq.'<') then 
c
               newline = newline(1:newlen)//'&lt'
               newlen = lenstr(newline)           
c
            else
c
               newline = newline(1:newlen)//datline(in:in)
               newlen = newlen + 1
c
            endif
c
         end do
c
c        newline = datline(1:len)
c
         len = lenstr(newline)
         write(htmlunit,*) newline(1:len)
c
      endif 

      goto 10
c
 100  call closehtml

      return

 200  write(0,*) 'ERROR Creating HTML data file'
      write(6,*) 'ERROR Creating HTML data file'
      call prc('ERROR Creating HTML data file')

      return
      end

