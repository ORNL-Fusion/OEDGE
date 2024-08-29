c     -*-Fortran-*-
c
c     @PROCESS NOOPT
      PROGRAM OUT
      use mod_params
      use error_handling
      use debug_options
      use mod_fp_data
      use mod_params
      use mod_outcom
      use mod_comtor
c slmod begin
      use mod_cgeom
c slmod end
      use allocate_storage_out
      implicit none
C
C  *********************************************************************
C  *                                                                   *
C  *  OUT   :  DRAWS GRAPHS.                                           *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
c     include 'outcom'
c     include 'comtor'
c
      integer iref,iopt,ierr,i1
      character*80 graph,label
      character*1024 desc
c 
      real time,time1,za02as       
c
c     Initialize debug tracing
c
      call init_trace(0,.false.)
c      call init_trace(0,.true.)
      call pr_trace('OUTMAIN','BEGIN EXECUTION')
c
c     Initialize parameter values
c
      call initialize_parameters
c      
c     Allocate dynamic storage
c      
c     jdemod
c     Move dynamic storage allocation into outinit AFTER the parameters
c     have been read from the RAW file in GET      
c
c      call allocate_dynamic_storage
c
      call fp_allocate_storage(ierr)     
C
C-----------------------------------------------------------------------
C     INITIALISATION
C-----------------------------------------------------------------------
C
      TIME1 = ZA02AS (1)

 10   call outinit

      call pr_trace('OUTMAIN','AFTER OUTINIT')

c
c------------------------------------------------------------------
c
c     BEGIN SETUP PROCESSING FOR PLOTTING
c
c------------------------------------------------------------------
c
      call global_plotsetup 

      call pr_trace('OUTMAIN','AFTER GLOBAL_PLOTSETUP')

C
C-----------------------------------------------------------------------
C     PRINT OUTS
C-----------------------------------------------------------------------
C
      call writedata

      call pr_trace('OUTMAIN','AFTER GLOBAL_WRITEDATA')
c
C-----------------------------------------------------------------------
C     HC PLOT INITIALIZATION
C-----------------------------------------------------------------------
c
      call global_hc_plot_init(crmb,crmi)

      call pr_trace('OUTMAIN','AFTER HC_PLOT_INIT')
c
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     M A I N   L O O P   F O R   G R A P H   O P T I O N S
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

  100 CONTINUE

      call plotloopinit(iopt,ierr)

      CALL RDG (IREF,GRAPH,IOPT,IERR)
c
      IF (IERR.NE.0) GOTO 9998
c
      if (iref.lt.10) goto 100 
c

c slmod begin
c...I WILL MOVE THIS TO A SEPARATE ROUTINE.
      IF     (iref.EQ.955) THEN
        IF (stepopt.EQ.0) THEN          
c...      Load the data specifying which steps to plot:
          stepopt = 1
          READ(5,'(A72)') graph
c          WRITE(0,*) 'GRAPH:',graph
          IF (graph(8:11).EQ.'Loop'.OR.graph(8:11).EQ.'LOOP'.OR.
     .        graph(8:11).EQ.'loop') THEN
            READ(graph,*) label,nsteplist,(steplist(i1),i1=1,nsteplist)

c            WRITE(0,*) label,nsteplist,(steplist(i1),i1=1,nsteplist)

c..         Print plots for all steps:
            IF (nsteplist.EQ.-1) nsteplist = 99
          ELSE
            CALL ER('Out','Expecting 955 data line',*9999)
          ENDIF
        ELSE
          IF (stepopt.LT.nsteplist) THEN
c...        Move back up in the input file in order to produce the plots for 
c           the next step in the step list:
            stepopt = stepopt + 1
            graph = '    '
            DO WHILE (graph(2:4).NE.'955')
              BACKSPACE 5
              BACKSPACE 5
              READ(5,'(A72)') graph
c              WRITE(0,*) 'iref up:',graph(1:10),stepopt,nsteplist
            ENDDO
          ELSE
            stepopt   = 0
            nsteplist = 0
          ENDIF
        ENDIF
      ELSEIF (iref.EQ.954.AND.iopt.GT.0) THEN
        CALL SetupSourcePlot(iref,graph,mode)
c...    Now multiplying the velocity by qtim when it is read in:
c        qt = 1.0
      ELSEIF (stepopt.GE.1) THEN
c...    Ignore all non secondary raw file plots while in loop:
        GOTO 100
      ELSEIF (restoresolution.OR.(mode.NE.0.AND.iopt.NE.0)) THEN
c        WRITE(0,*) 'RESTORING BASE SOLUTION DATA'
c        CALL GET (TITLE,desc,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)
c        CALL GET (TITLE,desc,NIZS,JOB,EQUIL,ITER,NITERS)
        CALL GET (desc)

        IF (REFCT.EQ.1) THEN
          CALL REFLECT
          WRITE(6,'('' EQUILIBRIUM GEOMETRY HAS BEEN REFLECTED'')')
        ENDIF

        mode = 0
        restoresolution = .FALSE.
c
c       Save the data read into the JOB string 
c        
        job_saved=job
c
      ENDIF

c
c      if (iref.ne.955.and.(.not.(iref.eq.954.and.iopt.gt.0))) then 
c          IF (stepopt.GE.1) THEN
c...       Ignore all non secondary raw file plots while in loop:
c           GOTO 100
c         elseIF (mode.NE.0.AND.iopt.NE.0) THEN
c           WRITE(0,*) 'RESTORING BASE SOLUTION DATA'
c           CALL GET (TITLE,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)
c           mode = 0
c         ENDIF
c      endif
c slmod end

      call load_additionalplotdata(iref,graph,iopt,ierr)

      if (ierr.ne.0) goto 100

c      write (0,'(a,i7,a,i7)') 'Processing input for plot number:',iref,
c     >                        ' with plot option:',iopt

      write (6,'(a,i7,a,i7)') 'Processing input for plot number:',iref,
     >                        ' with plot option:',iopt
      write (6,'(a,3i5,a,a)') 'IREF:',iref,iopt,ierr,' rest:',graph


      if (iref.gt.10.and.iref.lt.100) then 
         call out000(iref,graph,iopt,ierr)
      elseif (iref.lt.200) then 
         call out100(iref,graph,iopt,ierr)
      elseif (iref.lt.300) then 
         call out200(iref,graph,iopt,ierr)
      elseif (iref.lt.400) then 
         call out300(iref,graph,iopt,ierr)
      elseif (iref.lt.500) then 
         call out400(iref,graph,iopt,ierr)
      elseif (iref.lt.600) then 
         call out500(iref,graph,iopt,ierr)
      elseif (iref.lt.700) then 
         call out600(iref,graph,iopt,ierr)
      elseif (iref.lt.800) then 
         call out700(iref,graph,iopt,ierr)
      elseif (iref.lt.900) then 
         call out800(iref,graph,iopt,ierr)
      elseif (iref.lt.1000) then 
         call out900(iref,graph,iopt,ierr)
      endif

C
C-----------------------------------------------------------------------
C     END OF PLOT JOB
C-----------------------------------------------------------------------
C
      GOTO 100
 9998 continue

c
c     IERR =1 is for end of file condition
c      
      if (ierr.eq.1) then 

         write(0,'(a)') 'OUTMAIN: END OF PLOT FILE REACHED'

      else

         call errmsg('OUTMAIN:PROGRAM OUT: ERROR READING PLOT INPUT'//
     >               ': IERR ',ierr)

      endif

 9999 CONTINUE
C     ___________________________
      IF (ITER.LT.NITERS) GOTO 10
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      TIME = ZA02AS (1) - TIME1
c
c     Place storage deallocation calls just before end of code
c
c     far periphery
c
      call fp_deallocate_storage
c
c     Deallocate dynamic storage
c      
      call deallocate_dynamic_storage
c
      WRITE (6,'(/,'' OUT: TOTAL TIME USED ='',G11.4,'' SEC'',/)') TIME
      CALL GREND
      STOP
C
 9001 FORMAT(/1X,'OUT:    ',A55,'  OPT',I3,' REF ',A3,/1X,79('-'))
 9010 FORMAT(/1X,'   IT    J   IK   JK   IS  ISX    SL        SU    ',
     >  '   KVALS   ',/1X,65('-'))
 9011 FORMAT(1X,6I5,2F10.4,1P,G11.4)
 9012 FORMAT(1X,'PLOT',I3,4X,A)
 9020 FORMAT(/1X,62('*'),/1X,'*',60X,'*',
     >  /1X,'*',18X,'RUN OF OUT VERSION ',A5,18X,'*',
     >  /1X,'*',18X,24('-'),18X,'*',/1X,'*',60X,'*',
     >  /1X,'* ',A58,' *',/1X,'*',60X,'*',/1X,'* ',A53,5X,' *',
     >  /1X,'*',60X,'*',/1X,'* ',A54,4X,' *',
     >  /1X,'*',60X,'*',/1X,62('*'),/)
 9021 FORMAT(/1X,62('*'),/1X,'*',60X,'*',
     >  /1X,'*',18X,A20,22X,'*',/1X,'*',60X,'*',/1X,62('*'),/)
 9031 FORMAT(/1X,' IK IR    R      Z     AREA',10(2X,A7))
 9032 FORMAT(1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,1P,E8.2,10E9.2)
 9034 FORMAT(29X , 1P , 10E9.2 )
 9040 FORMAT(1X,'IZ, WAVELENGTH: ',A)
      END
