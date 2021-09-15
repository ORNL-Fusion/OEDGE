C  *********************************************************************
C  *                                                                   *
C  *  SYSRS6K:THIS FILE CONTAINS DUMMIES FOR ALL THE EXTERNAL ROUTINES *
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
c  *    ioname:  This routine finds the name of file attached to a     *
c  *             specified i/o stream.                                 *
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
      implicit none
      INTEGER IFLAG
      WRITE (6,'('' XUFLOW: DUMMIED OUT FOR THIS APPLICATION.'')')
      RETURN
      END
C
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
      REAL FUNCTION ZA02AS (IFLAG)
      implicit none
      INTEGER I
      I = MCLOCK()
      ZA02AS = I/100.0
CHOT  ZA02AS = 0.0
C      ZA02AS = SECOND ()
      RETURN
      END
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
      implicit none
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
      NIMTIM = ZA02AS(1)
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
      exeline(len+1:len+2) = x'0'
c
      retcode = SYSTEM(exeline(1:len+1))
c
c
C     FOR USE IF PIN IS INCLUDED AS A SUBROUTINE.
C
C      CALL PINPGX
C
C
      NIMTIM = ZA02AS(1) - NIMTIM
      WRITE(6,*) 'TIME USED IN HYDROGEN NEUTRAL CODE:',NIMTIM,' (S)'
      write(6,*) '- PIN RETURN CODE = ',retcode
C
      RETURN
      END
c
c
c
      subroutine printerinit
      implicit none
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
      use mod_params
      use mod_cnoco
      implicit none
      integer nksir
c     include 'params'
c     include 'cnoco'
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
      use mod_params
      use mod_cnoco
      implicit none
      integer nksir
c     include 'params'
c     include 'cnoco'
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
c     This is somewhat of a kludge that returns the NAME of the 
c     file associated with UNIT 5 by calling UNIX OS level 
c     functions.
c
      integer   l1, l2
      character cmd*80, cmda*40, cmdb*40, cmdc*2, cmdd*2
      data      cmda/'ls -l fort.'/,
     >          cmdb/' 2>/dev/null | cut -f2 -d''->'' > ./fort.'/
c
      ierr = 0
c
      filename = ' '
c
c  only allow i/o streams up to 99 for now
c
      if (iunit.lt.1 .or. iunit.gt.99 .or.
     >    itemp.lt.1 .or. itemp.gt.99) then
        ierr = 1
        return
      endif
c
c  write filename to a temporary file
c
      write(cmdc,'(i2)') iunit
      l1 = 2 - min0(1,iunit/10)
      write(cmdd,'(i2)') itemp
      l2 = 2 - min0(1,itemp/10)
      cmd = cmda(1:11)//cmdc(l1:2)//cmdb(1:39)//cmdd(l2:2)
c
      call system(cmd)
c
c  read filename from temporary file
c
      rewind(itemp)
      read(itemp,'(2x,a80)',end=20) filename
      goto 10
   20 ierr = 2
c
c  delete temporary file
c
   10 close(itemp)
      cmd = 'rm fort.'//cmdd(l2:2)
      call system(cmd)
c
      return
c
      end
c
c     
c
      subroutine killdiv
      implicit none
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
      implicit none
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
      call signal(SIGUSR1,killdiv)
      return
      end

c
C================================================================
c
      SUBROUTINE DMGUID(SYSUID,PREFIX)
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
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
      
