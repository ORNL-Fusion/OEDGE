C  *********************************************************************
C  *                                                                   *
C  *  SYSWIN32:THIS FILE CONTAINS DUMMIES FOR ALL THE EXTERNAL ROUTINES*
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
C     IBM SYSTEM ROUTINE TO PREVENT UNDERFLOW INTERRUPTS OCCURING,
C     REPLACE BY DUMMY ON OTHER ARCHITECTURES.
C
      SUBROUTINE XUFLOW (IFLAG)
      INTEGER IFLAG
      WRITE (6,'('' XUFLOW: DUMMIED OUT FOR UNIX OPERATING SYSTEM.'')')
      RETURN
      END
C
C***********************************************************************
C
C     ZV01AD
C     ======
C     HARWELL LIBRARY ROUTINE TO EXTRACT DATASET NAME CONNECTED
C     TO CHANNEL 5, AND TO PUT THE NAME IN JFCB(1:44) AND THE
C     MEMBER IN JFCB(45:52).
C     SINCE THERE IS NO EQUIVALENT SYSTEM ROUTINE IN UNIX, THIS
C     CALL SIMPLY RETURNS A FIXED DESCRIPTION STRING.
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
C     HARWELL LIBRARY FUNCTION TO EXTRACT CPU TIME USED SO FAR.
C     FOR HOTSPOT ANALYSIS, DUMMY OUT BY SETTING ZA02AS = 0.0
C     OTHERWISE, REPLACE WITH EQUIVALENT SYSTEM FUNCTION.
C
      REAL FUNCTION ZA02AS (IFLAG)
      USE IFPORT
      REAL(4) TIME(2)
      integer iflag
c
      time(1) = 0.0
      time(2) = 0.0
c
      ZA02AS=ETIME(TIME)
c
      RETURN
      END
C
C***********************************************************************
C
C     RANINI
C     ======
C     INTERFACE TO RANDOM NUMBER INITIALISER SYSTEM ROUTINE
C
c
c     first define modules necessary for MKL and for state keeping
c
      include 'mkl_vsl.fi'

      module randstate
      use MKL_VSL_TYPE
      TYPE (VSL_STREAM_STATE) :: stream
      integer nindex
      end module randstate
     
      SUBROUTINE RANINI (ISEED)
      INTEGER ISEED
c
      USE MKL_VSL_TYPE
      USE MKL_VSL
      USE randstate

      integer brng, status
      write(6,*) 'Random seed:',iseed
c     this is the init for the INTEL MKL random number generator
c
      brng=VSL_BRNG_R250
      status=vslnewstream(stream,brng,iseed)

c     this is a running index for generation of single random numbers

      nindex=0

      RETURN
      END
C
C***********************************************************************
C
C     SURAND
C     ======
C     ESSL LIBRARY ROUTINE TO GENERATE VECTOR OF RANDOM NUMBERS,
C     REPLACED BY INTERFACE TO RANDOM GENERATOR WITHIN A DO LOOP.
C
      SUBROUTINE SURAND (SEED,NRANDS,CRANDS)
      use rand_data
      implicit none
      DOUBLE PRECISION SEED
      INTEGER NRANDS
      REAL CRANDS(NRANDS)

      USE MKL_VSL_TYPE
      USE MKL_VSL
      USE randstate

      INTEGER method, status
      REAL a,b
      
      SAVE a,b,method

      DATA a,b,method /0.0,1.0,0/

      status=vdrnguniform(method,stream,nrands,crands,a,b)
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
      USE MKL_VSL_TYPE
      USE MKL_VSL
      USE randstate
      INTEGER method, status
      REAL randbuf(1000), a, b
      SAVE a,b,method,randbuf

      DATA a,b,method /0.0,1.0,0/

      if (nindex.eq.0) then
        nindex=1000
        status=vdrnguniform(method,stream,nindex,randbuf,a,b)        
      endif
      rands = randbuf(nindex)
      nindex = nindex - 1
 
c     Update random number count
c     
      ran_used = ran_used + 1.0
c
      RETURN
      END


      subroutine killrand
      USE MKL_VSL_TYPE
      USE MKL_VSL
      use randstate
      integer status

      status=vsldeletestream(stream)
      
      return
      end
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
      USE MKL_VSL_TYPE
      USE MKL_VSL
      USE randstate
      INTEGER method, status
      REAL randbuf(1000), a, b
      SAVE a,b,method,randbuf

      DATA a,b,method /0.0,1.0,0/

      if (nindex.eq.0) then
        nindex=1000
        status=vdrnguniform(method,stream,nindex,randbuf,a,b)        
      endif
      getranf = randbuf(nindex)
      nindex = nindex - 1
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
      USE MKL_VSL_TYPE
      USE MKL_VSL
      USE randstate
      INTEGER method, status
      REAL randbuf(1000), a, b
      SAVE a,b,method,randbuf

      DATA a,b,method /0.0,1.0,0/

      if (nindex.eq.0) then
        nindex=1000
        status=vdrnguniform(method,stream,nindex,randbuf,a,b)        
      endif
      rand = randbuf(nindex)
      nindex = nindex - 1
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
      integer time
      external time
c
      REAL ZA02AS
      EXTERNAL ZA02AS
c
c     Assign the return code to zero for now
c
      retcode = 0
c
      NIMTIM = time()
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
      NIMTIM = time() - NIMTIM
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
      USE IFPORT
      implicit none
      integer ierr
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
      USE IFPORT
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
      call killrand
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
      integer ierr
      integer SIGUSR1
      parameter (SIGUSR1=30)
      external killdiv
c
C     SIGUSR1 - comment out for mainframe - or move to routine in the
c               system module.
C
      ierr = signal(SIGUSR1,killdiv, -1)
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

      IF (code.NE.0) CALL ER('CIssue','Invalid command',*99)

      RETURN
99    STOP
      END


      REAL FUNCTION Clock2()
      USE IFPORT
      IMPLICIT none	

      REAL(4) val(2)

      Clock2 = etime(val)

      RETURN
      END
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
c     Create file names for jpeg image output
c
c     Note - null termination is added in the C routine
c
      if (flag.eq.1) then
c
         len = lenstr(fbase)

         if (image_count.lt.10) then
            write(fname,'(a,i1,a,Z1)') fbase(1:len)//'_image0',
     >                               image_count,'.jpg'
         else
            write(fname,'(a,i2,a,Z1)') fbase(1:len)//'_image',
     >                               image_count,'.jpg'
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
C RETURNS USERID (not sure if that works in Win32 environment)
C
      CHARACTER*(*) SYSUID, PREFIX
C
c     CALL GETENV('LOGNAME',SYSUID)
C
      PREFIX = ' '
      SYSUID = ' '
C
      RETURN
      END
c
C================================================================
c
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

