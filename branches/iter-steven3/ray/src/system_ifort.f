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
      INTEGER I,MCLOCK,IFLAG
      I = MCLOCK()
      ZA02AS = I/100.0
CHOT  ZA02AS = 0.0
C      ZA02AS = SECOND ()
      RETURN
      END
c
c
c
c ======================================================================
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
c ======================================================================
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
      filename = TRIM(filename)
c
      return
c
      end
c
c ======================================================================
c
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
c
c ======================================================================
c
      REAL FUNCTION Clock2()
      IMPLICIT none

      REAL etime

      REAL val(2)

      Clock2 = etime(val)

c      Clock = val(1) + val(2)

      RETURN
      END
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
      SUBROUTINE SURAND2 (SEED,NRANDS,RANDS)
      use rand_data
      implicit none
c
      DOUBLE PRECISION SEED
      INTEGER NRANDS
      REAL RANDS

      call random_number(rands)

      RETURN
      END
c slnote - moved here from rundiv.d6a
c
      integer function calc_random_seed(n)
      implicit none
      integer n
c
c     Calculates a random number seed based on the current date and 
c     time. The seed number will be limited to n digits if n is specified
c     greater than zero.
c 
c     Local varaiables
c
      integer j,ciseed
      character*8 systim,sysdat
c
      CALL ZA08AS (SYSTIM)                                              
      CALL ZA09AS (SYSDAT)                                              
c
      CISEED = 1                                                      
      DO 300 J = 1, 8                                                 
         IF (J.EQ.3.OR.J.EQ.6) GOTO 300                                
         CISEED = CISEED * (ICHAR(SYSTIM(J:J)) - ICHAR('0') + 1) *     
     >                      (ICHAR(SYSDAT(J:J)) - ICHAR('0') + 1)       
 300  CONTINUE                                                        
      CISEED = CISEED - 1                                             
c
      if (n.gt.0) then 
c
         do while (ciseed.ge.10**n)
            ciseed = ciseed /10
         end do  
c
      endif 
c
      calc_random_seed = ciseed
c
      return
      end 
c
c
c
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

