c
C  *********************************************************************
C  *                                                                   *
C  *  SYSPGI :THIS FILE CONTAINS DUMMIES FOR ALL THE EXTERNAL ROUTINES *
C  *  USED THROUGHOUT DIV.   THE PROGRAM WAS DESIGNED AND OPTIMISED    *
C  *  ON THE IBM3090E AT JET AND MAKES USE OF SOME IBM UTILITIES,      *
C  *  THE HARWELL SUBROUTINE LIBRARY ("HSL") AND THE NOCORONA PACKAGE  *
C  *  ALSO DEVELOPED AT JET.                                           *
C  *    THE CODE COLLECTED HERE WILL REQUIRE MODIFICATIONS WHEN        *
C  *  TRANSPORTING BETWEEN MACHINES, WHILST THE REMAINDER OF THE DIV   *
C  *  CODE SHOULD REMAIN INTACT.                                       *
C  *                                                                   *
C  *  ROUTINES IN THIS FILE :-                                         *
C  *    XUFLOW: PREVENT UNDERFLOW INTERRUPTS                           *
C  *    SURAND: GET VECTOR OF RANDOM NUMBERS IN (0,1)                  *
C  *    ZA08AS: GET TIME OF DAY AS 8 CHARACTER STRING                  *
C  *    ZA09AS: GET DATE AS 8 CHARACTER STRING                         *
C  *    ZV01AD: GET NAME OF DATAFILE CONNECTED TO CHANNEL 5            *
c  *            - replaced by ioname                                   *
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
c  *    ioname:  This routine finds the case name from an environment  *
c  *             variable                                              *
c  *    CISSUE  :Issue OS command                                      *
c  *    CLOCK2  :Timing function                                       *
c  *    DMGUID  :Routine to get userid for person executing the code   *
C  *                                                                   *
C  *                                      C.M.FARRELL   JUNE 1988      *
C  *                                      L.D.HORTON    APRIL 1993     *
C  *                                      J.D.ELDER                    *
C  *                                                                   *
C  *********************************************************************
C
C     XUFLOW
C     ======
C     IBM  : SYSTEM ROUTINE TO PREVENT UNDERFLOW INTERRUPTS OCCURING
C     CRAY : REPLACE WITH DUMMY ROUTINE HERE.
C
      SUBROUTINE XUFLOW (IFLAG)
      INTEGER IFLAG
      WRITE (6,'('' XUFLOW: DUMMIED OUT FOR THIS APPLICATION.'')')
      RETURN
      END
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
c      SUBROUTINE ZV01AD (IUNIT, VSN, DSN, JFCB)
c      CHARACTER VSN*8,DSN(3)*8,JFCB*176
c      WRITE (6,'('' ZV01AD: DUMMIED OUT FOR THIS APPLICATION.'')')
c      JFCB = 'TORONTO - RS6000 '
c      RETURN
c      END
c
C
C***********************************************************************
C
C     ZA02AS
C     ======
C     IBM  : HARWELL LIBRARY FUNCTION TO EXTRACT CPU TIME USED SO FAR.
C     HOT  : FOR HOTSPOT ANALYSIS, DUMMY OUT BY SETTING ZA02AS = 0.0
C     CRAY : REPLACE WITH SYSTEM FUNCTION SECOND.
C
c slmod begin 
      REAL FUNCTION ZA02AS (IFLAG)
      INTEGER I,MCLOCK,IFLAG

      INTEGER clock_rate,t1
      REAL    t2

      IF     (iflag.EQ.1) THEN
        CALL CPU_TIME(t2)
        ZA02AS = t2
c        za02as = REAL(mclock()) / 100.0
c        rdum1  = etime(val)
c        za02as = val(1) 
      ELSEIF (iflag.EQ.2) THEN
        CALL SYSTEM_CLOCK(t1,clock_rate)
        za02as = REAL(t1) / REAL(clock_rate)
      ELSE
        CALL ER('ZA02AS','Unknown IFLAG value',*99)
      ENDIF

c make this an option with a default set to MCLOCK()
c      write(0,*) 'TIME TEST:',za02as
c
c      REAL FUNCTION ZA02AS (IFLAG)
c      INTEGER I,MCLOCK,IFLAG
c
c      real test
c      I = MCLOCK()
c      ZA02AS = I/100.0
c slmod end

CHOT  ZA02AS = 0.0
C      ZA02AS = SECOND ()
      RETURN
c slmod begin
99    STOP
c slmod end
      END
c
c
c
C
C***********************************************************************
C
C     INVPIN
C     ======
C     IBM  : CALL PIN AS A SUBROUTINE
C     CRAY : CALL PIN AS A SUBROUTINE
C     UNIX : START PIN AS AN INDEPENDENT PROCESS
C
      SUBROUTINE INVOKEPIN(ACTPIN,NIMTIM,retcode)
      CHARACTER*(*) ACTPIN
      real nimtim
      integer retcode
      integer system
C
C     THIS SUBROUTINE PROVIDES THE INTERFACE TO PIN
C     AT THIS TIME IT IS IMPLEMENTED EITHER AS A CALL
C     TO A SUBROUTINE THAT HAS BEEN BOUND TO DIVIMP OR
C     AS A CALL TO A "SYSTEM" SUBROUTINE THAT STARTS UP
C     PIN AS AN INDEPENDENT PROCESS. IT IS IMPORTANT TO
C     NOTE THAT EITHER METHOD CAN BE USED AT VARIOUS SITES.
C     THE CHOICE OF METHOD RELIES ON SYSTEM MEMORY RESTRICTIONS
C     EXECUTION OVERHEAD, AND EASE OF IMPLENTATION AND
C     MAINTENANCE.
C
C     DAVID ELDER, FEB 4 1993
C
C
c     NOTE!: Actpin MUST be converted to a null terminated string
c
      character*255 exeline
c
      integer lenstr,len
      external lenstr
c
      REAL ZA02AS
      EXTERNAL ZA02AS
c
c     Assign the return code to zero for now
c
      retcode = 0
c
c slmod begin
c...  Added optional to ZA02AS that allows a measurement
c     of real time rather than CPU time:
      NIMTIM = ZA02AS(2)
c
c      NIMTIM = ZA02AS(1)
c slmod end
C
C     FOR USE ON A UNIX SYSTEM OR A PROPERLY SET UP MVS SYSTEM
C
c      CALL SYSTEM(ACTPIN)
c
      len = lenstr(actpin)
c
      exeline = actpin(1:len)
c
c     Add NULL termination required for 'C' strings
c
c      exeline(len+1:len+2) = '\0'
c

c      WRITE(0,*) 'INVOKEPIN: "'//exeline(1:len+1)//'"'
c
c      retcode = SYSTEM(exeline(1:len+1))
c
      retcode = SYSTEM(exeline(1:len))
c
c
C     FOR USE IF PIN IS INCLUDED AS A SUBROUTINE.
C
C      CALL PINPGX
C
C
c slmod begin
      NIMTIM = ZA02AS(2) - NIMTIM
c slmod end
      WRITE(6,*) 'TIME USED IN HYDROGEN NEUTRAL CODE:',NIMTIM,' (S)'
      write(6,*) '- PIN RETURN CODE = ',retcode
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
      call hrdlin(1)
      call hrdchr(1)
c
      return
      end
c
c
c
c      subroutine killdiv
c
c     This routine calls the entry point in the div.d3a
c     to wrap up processing and then exit.
c
c     This is not supported on the IBM mainframe - it is only
c     supported under UNIX.
c
c      call divkill
c      stop
c      end
c
c
c
c      subroutine initkill
c
c     DEFINE the SIGUSR1 kill signal so that the
c     signal call can trap it - if it is sent
c     to the DIVIMP process. (Workstation only)
c
c      integer SIGUSR1
c      parameter (SIGUSR1=30)
c      external killdiv
c
C     SIGUSR1 - comment out for mainframe - or move to routine in the
c               system module.
C
c      call signal(SIGUSR1,killdiv)
c      return
c      end
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
      subroutine ioname(filename,ierr)
      implicit none
      integer  ierr
      character filename*(*)
c
c     Get name from environment variable
c     
       integer len1,lenstr
       external lenstr 
c
c
      ierr = 0
c
      filename = ' '

      CALL GetEnv('CASENAME',filename)
c
      len1 = lenstr(filename)
      filename = filename(1:len1)//'.d6i'  
c
      return
c
      end
c
c
c
      subroutine casename(filename,ierr)
      implicit none
      integer  ierr
      character filename*(*)
c
c     Get name from environment variable
c     
      integer len1,lenstr
      external lenstr 
c
      ierr = 0
c
      filename = ' '

      CALL GetEnv('CASENAME',filename)
c
      len1 = lenstr(filename)
      filename = filename(1:len1)
c
      return
c
      end
c
c
c
      subroutine get_div_exec_dir(dirname,ierr)
      implicit none
      integer  ierr
      character dirname*(*)
c
c     Get name from environment variable
c     
      integer len1,lenstr
      external lenstr 
      character*256 :: divhome, divmaindir

c
      ierr = 0
c
      dirname = ' '

c
c     This fix is to work around Steve's scripts where
c     DIVMAINDIR is not yet defined
c
      CALL GetEnv('DIVHOME',divhome)
      CALL GetEnv('DIVMAINDIR',divmaindir)
c
c     If divmaindir is not defined then use divhome
c
      if (len_trim(divmaindir).le.1) then 
         dirname = trim(divhome)
      else
         dirname = trim(divmaindir)
      endif
c
      return
c
      end
c
c
c
      subroutine get_div_data_dir(dirname,ierr)
      implicit none
      integer  ierr
      character dirname*(*)
c
c     Get name from environment variable
c     
      integer len1,lenstr
      external lenstr 
      character*256 :: divdata

c
c     Initial state is invalid data .. if non-zero length string is 
c     returned then the return code is set to 1
c
      ierr = 1
c
      dirname = ' '
c
c     Get data directory from environment
c
      CALL GetEnv('DIVDATDIR',divdata)

      if (len_trim(divdata).gt.0) ierr = 0

      dirname = trim(divdata)

      return
c
      end
c
c
c
      SUBROUTINE run_system_command(cmd,retcode)
      CHARACTER*(*) CMD
      integer retcode
      integer system
C
c     This routine calls the "system" command to issue
c     an OS level command from inside the code. 
c
c     This is mostly useful for file manangement
c
c
c     Initialize the return code to zero for now
c
      retcode = 0
c
c
      retcode = SYSTEM(trim(cmd))
c
c
      RETURN
      END
c
c
c
      subroutine killdiv
c
c     This is SYSTEM specific code that is applicable ONLY to DIVIMP
c
c     However - so that this code may be included in the common modules
c     a routine called divkill - the entry point for this code - has
c     been added to the ioout code module.
c
c
c     This routine calls the entry point in the div and out modules
c     to wrap up processing and then exit.
c
c     This is not supported on the IBM mainframe - it is only
c     supported under UNIX.
c
      call divkill
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
      integer SIGUSR1
      parameter (SIGUSR1=30)
      external killdiv
c
C     SIGUSR1 - comment out for mainframe - or move to routine in the
c               system module.
C
c      call signal(SIGUSR1,killdiv)
      return
      end

      SUBROUTINE CISSUE(command,code)
      IMPLICIT none

      CHARACTER command*(*)
      INTEGER   code

      INTEGER System

      INTEGER   len
      CHARACTER exeline*256

      code = 0

      WRITE(exeline,'(256X)')

      len = LEN_TRIM(command)

      exeline = command(1:len)

c...  Add NULL termination required for 'C' strings:
c      exeline(len+1:len+2) = '\0'

      WRITE(6,*) 'CISSUE: "'//exeline(1:len+1)//'"'

      code = SYSTEM(exeline(1:len+1))

c slmod begin
      IF (code.NE.0) THEN
        WRITE(0,*) 
        WRITE(0,*) '---------------------------------------------------'
        WRITE(0,*) ' NO LONGER STOPPING THE CODE ON FAILED SYSTEM CALL'
        WRITE(0,*) '---------------------------------------------------'
        WRITE(0,*) 
      ENDIF
c
c      IF (code.NE.0) CALL ER('CIssue','Invalid command',*99)
c slmod end

      RETURN
99    STOP
      END


      REAL FUNCTION Clock2()
      IMPLICIT none

      REAL etime,test

      REAL val(2)

c      STOP 'call to CLOCK2, check the code'

      Clock2 = etime(val)

      CALL cpu_time(test)

c      write(0,*) 'TIME TEST:',Clock2,test

      Clock2 = test

c      Clock = val(1) + val(2)

      RETURN
      END
c
c
c
c      REAL FUNCTION ACOSD(val)
c      IMPLICIT none
c
c      REAL val
c
c
c      acosd = ACOS(val) * 180.0 / 3.1415
c
c      RETURN
c      END
c
c
c
c      DOUBLE PRECISION FUNCTION DFLOAT(i)
c
c      IMPLICIT none
c
c      INTEGER i
c
c      dfloat = DBLE(i)
c
c      RETURN
c      END
c
c
c
c      SUBROUTINE GetEnv(param,value)
c      IMPLICIT none
c
c      CHARACTER *(*) param,value
c
c      CALL DOSPARAM@(param,value)
c
c      WRITE(0,*) 'GetEnv: ',param,value
c      RETURN
c      END
c
c
c
      subroutine create_file_name(fname,fbase,image_count,flag)
      implicit none
      character*(*) fname,fbase
      integer image_count,flag
c
c     This routine writes system dependent file names.
c
      integer lenstr,len
      external lenstr
c
      if (flag.eq.1) then 
c
         len = lenstr(fbase)

         if (image_count.lt.10) then 
            write(fname,'(a,i1,a,Z1)') fbase(1:len)//'_image0',
     >                               image_count,'.jpg\0'
         else
            write(fname,'(a,i2,a,Z1)') fbase(1:len)//'_image',
     >                               image_count,'.jpg\0'
         endif
c
      endif

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

c
c     jdemod
c
c     Functions that were previously linked to C code
c
c     DATE and TIME
c

C
C***********************************************************************
C
C     ZA08AS
C     ======
C     IBM  : HARWELL LIBRARY ROUTINE TO EXTRACT TIME IN 8 CHARACTERS
C     CRAY : REPLACE WITH CALL TO CLOCK SYSTEM ROUTINE.
C
      SUBROUTINE ZA08AS (SYSTIM)
      CHARACTER*8 SYSTIM
c
c     jdemod
c
c     Use fortran 90+ intrinsic function
c
      character*10 tim_f

      call date_and_time(time=tim_f)

      systim= tim_f(1:2)//':'//tim_f(3:4)//':'//tim_f(5:6)

c      CALL CLOCKTEST (SYSTIM)

      RETURN
      END
C
C***********************************************************************
C
C     ZA09AS
C     ======
C     IBM  : HARWELL LIBRARY ROUTINE TO EXTRACT DATE IN 8 CHARACTERS
C     CRAY : REPLACE WITH CALL TO DATE SYSTEM ROUTINE.
C
      SUBROUTINE ZA09AS (SYSDAT)
      CHARACTER*8 SYSDAT
c     
c     jdemod
c
c     use fortran 90+ intrinsic function 
c

      character*8 date_f

      call date_and_time(date=date_f)

      sysdat=date_f(5:6)//'/'//date_f(7:8)//'/'//date_f(3:4)

c
c      CALL DATETEST (SYSDAT)
C       SYSDAT = '04/11/88'
      RETURN
      END


C
C     RANDOM NUMBERS
C


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
c
c     Initialization of the intrinsic generator is more complicated
c     than previously implemented - 34 integers are typically 
c     required - these are generated from one input seed value
c     in this routine.
c

c
      integer k,i
      integer,allocatable :: temp_seed(:)

      call random_seed(size=k)

      write(6,*) 'Random seed:',iseed
      write(6,*) 'Random seed size = ',k

      allocate(temp_seed(k))

      do i = 1,k
         temp_seed(i) = int(iseed/i)
      end do

      write(6,*) 'Random Seed Used:',(temp_seed(i),i=1,k)

      call random_seed(put=temp_seed)

      deallocate(temp_seed)

c
c     write(6,*) 'Random seed:',iseed
c     CALL NEWSRAND (ISEED)
C
C     CALL SRAND (ISEED ) FOR CRAY OR IBM DEFAULT GENERATOR
C
      RETURN
      END



C
C***********************************************************************
C
C     SURAND
C     ======
C     IBM  : ESSL LIBRARY ROUTINE TO GENERATE VECTOR OF RANDOM NUMBERS
C     CRAY : REPLACE WITH CALLS TO RANF GENERATOR WITHIN A VECTORISABLE
C            LOOP.
C
      SUBROUTINE SURAND (SEED,NRANDS,CRANDS)
      use rand_data
      implicit none
      DOUBLE PRECISION SEED
      INTEGER NRANDS
      integer j
      REAL CRANDS(NRANDS)
c
c      real newrand
c      EXTERNAL NEWRAND
c

      DO 100 J = 1, NRANDS
          call random_number(crands(j))
c
c          CRANDS(J) = NEWRAND()
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
  100 CONTINUE
      
c
c     Update random number count
c     
      ran_used = ran_used + nrands

      RETURN
      END
c
c
c
      SUBROUTINE SURAND2 (SEED,NRANDS,RANDS)
      use rand_data
      implicit none
c
      DOUBLE PRECISION SEED
      INTEGER NRANDS
      REAL RANDS
c
c      real NEWRAND
c      EXTERNAL NEWRAND
c
      call random_number(rands)
c
c      RANDS = NEWRAND()
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
c
c     Update random number count
c     
      ran_used = ran_used + 1.0
c
      RETURN
      END
c
c
c
      real function getranf ()
      use rand_data
      implicit none
c
c     This routine interfaces with a C language routine which
c     accesses the C-library random number generator. It is intended
c     to ultimately replace surand2 for single random numbers - there 
c     is no reason to pass the seed to the routine generating the 
c     random numbers - only to the initialization routine. SURAND2 was
c     created for historical reasons - since it mimics the calling
c     convention of an IBM mainframe built-in function to generate 
c     random numbers. Neither SURAND nor SURAND2 have any use for
c     a seed value in their current incarnation. 
c
c     The value returned here is between 0.0 and 1.0
c
      REAL NEWRAND
c
c      EXTERNAL NEWRAND
c
      call random_number(newrand)
      GETRANF = NEWRAND
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
c
c     Update random number count
c     
      ran_used = ran_used + 1.0
c
      RETURN
      END
c
c
c
      subroutine getran (rand)
      use rand_data
      implicit none
c
c     This routine interfaces with a C language routine which
c     accesses the C-library random number generator. It is intended
c     to ultimately replace surand2 for single random numbers - there 
c     is no reason to pass the seed to the routine generating the 
c     random numbers - only to the initialization routine. SURAND2 was
c     created for historical reasons - since it mimics the calling
c     convention of an IBM mainframe built-in function to generate 
c     random numbers. Neither SURAND nor SURAND2 have any use for
c     a seed value in their current incarnation. 
c
c     The value returned here is between 0.0 and 1.0
c
c     Subroutine version of getran() function
c
      REAL rand
c
c      real newrand
c      EXTERNAL NEWRAND
c
c
      call random_number(rand)

c
c      RAND = NEWRAND()
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C

c
c     Update random number count
c     
      ran_used = ran_used + 1.0
c
      RETURN
      END
