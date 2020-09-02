C  *********************************************************************
C  *                                                                   *
C  *  DUMMY: This file contains dummies for all the external routines  *
C  *  used throughout LIM3.  The program was designed and optimised    *
C  *  on the IBM3090E at JET and makes use of some IBM utilities,      *
C  *  the Harwell subroutine library ("HSL") and the NOCORONA package  *
C  *  also developed at JET.                                           *
C  *    The code collected here will require modifications when        *
C  *  transporting between machines, whilst the remainder of the LIM3  *
C  *  code should remain intact.                                       *
C  *                                                                   *
C  *  Routines in this file :-                                         *
C  *    XUFLOW: Prevent underflow interrupts                           *
C  *    SURAND: Get vector of random numbers in (0,1)                  *
C  *    ZA08AS: Get time of day as 8 character string                  *
C  *    ZA09AS: Get date as 8 character string                         *
C  *    ZV01AD: Get name of datafile connected to channel 5            *
C  *    TB04A : Create cubic spline                                    *
C  *    TG01B : Extract value along newly created cubic spline         *
C  *    VC03A : Produce best fit curve through data                    *
C  *    ZA02AS: Get time used so far in seconds                        *
C  *    RANINI: Dummy interface to RANSET random number initialiser    *
C  *                                                                   *
C  *                                      C.M.Farrell   June 1988      *
C  *                                                                   *
C  *********************************************************************
C
C     XUFLOW
C     ======
C     IBM  : System routine to prevent underflow interrupts occuring
C     CRAY : Replace with dummy routine here.
C
      SUBROUTINE XUFLOW (IFLAG)
      implicit none
      INTEGER IFLAG
      WRITE (6,'('' XUFLOW: dummied out for this application.'')')
      RETURN
      END
C
C***********************************************************************
C
C     ZA08AS
C     ======
C     IBM  : Harwell library routine to extract Time in 8 characters
C     CRAY : Replace with call to CLOCK system routine.
C
c      SUBROUTINE ZA08AS (SYSTIM)
c      CHARACTER*8 SYSTIM
c       CALL CLOCKTEST (SYSTIM)
c      RETURN
c      END
      SUBROUTINE ZA08AS (SYSTIM)
      implicit none
      CHARACTER*8 SYSTIM

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
C     IBM  : Harwell library routine to extract Date in 8 characters
C     CRAY : Replace with call to DATE system routine.
C
c      SUBROUTINE ZA09AS (SYSDAT)
c      CHARACTER*8 SYSDAT
c      CALL DATETEST (SYSDAT)
C      SYSDAT = '04/11/88'
c      RETURN
c      END
c
      SUBROUTINE ZA09AS (SYSDAT)
      implicit none
      CHARACTER*8 SYSDAT

      character*8 date_f

      call date_and_time(date=date_f)

      sysdat=date_f(5:6)//'/'//date_f(7:8)//'/'//date_f(3:4)


c      CALL DATETEST (SYSDAT)
C       SYSDAT = '04/11/88'

      RETURN
      END
C
C***********************************************************************
C
C     ZV01AD
C     ======
C     IBM  : Harwell library routine to extract Dataset name connected
C            to channel 5, and to put the name in JFCB(1:44) and the
C            member in JFCB(45:52).
C     CRAY : No equivalent system routine.
C
      SUBROUTINE ZV01AD (IUNIT, VSN, DSN, JFCB)
      implicit none
      integer :: iunit
      CHARACTER VSN*8,DSN(3)*8,JFCB*176
      WRITE (6,'('' ZV01AD: dummied out for this application.'')')
      JFCB = 'SYSTEM.TEST.FOR.IBM(TESTONE)      '
      RETURN
      END
C
C
C***********************************************************************
C
C     ZA02AS
C     ======
C     IBM  : Harwell library function to extract CPU time used so far.
C     HOT  : For hotspot analysis, dummy out by setting ZA02AS = 0.0
C     CRAY : Replace with system function SECOND.
C
      REAL FUNCTION ZA02AS (IFLAG)
      implicit none
      integer :: iflag  ! does nothing ... compatibility
C      integer,external :: mclock
c     mclock is not a universal fortran90 intrinsic
c     change to using etime for now
c      integer,intrinsic :: mclock ! change to intrinsic for gfortran
      real :: etime 
      real vals(2)
      za02as = etime(vals)
c     INTEGER I
c      I = MCLOCK()
C      ZA02AS = I/100.0
      ZA02AS = 0.0
C     ZA02AS = SECOND ()
      RETURN
      END
C
C***********************************************************************
C
C     RANINI
C     ======
C     IBM  : Dummy routine - no need to call RANSET
C     CRAY : Interface to random no. initialiser system routine RANSET
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
