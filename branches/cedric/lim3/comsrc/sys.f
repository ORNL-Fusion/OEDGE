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
      INTEGER I
      I = MCLOCK()
      ZA02AS = I/100.0
CHOT  ZA02AS = 0.0
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
      INTEGER ISEED
      CALL NEWSRAND (ISEED)
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
C     IBM  : ESSL library routine to generate vector of random numbers
C     CRAY : Replace with calls to RANF generator within a vectorisable
C            loop.
C
      SUBROUTINE SURAND (SEED,NRANDS,RANDS)      
      DOUBLE PRECISION SEED
      INTEGER NRANDS  
      REAL RANDS(NRANDS),NEWRAND
      EXTERNAL NEWRAND
      DO 100 J = 1, NRANDS
          RANDS(J) = NEWRAND()
C         
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
  100 CONTINUE
      RETURN
      END
