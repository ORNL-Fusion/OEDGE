      program testran
      implicit none
c
c     Program for testing the random number generator
c
      
      real minran,maxran
c      real*8 minran8,maxran8

      integer,parameter:: ranbins=10000
      
      real randist1(ranbins),randist2(ranbins)
      

      real ranval,lastranval
      real ranres
      real,parameter:: rancut = 0.1
      real getranf
      real*8 seed
      real ranmin,ranmax,ransum1,ransum2
      
      integer in, rancnt,calc_random_seed

      integer,parameter:: ranmaxcnt=100000000
      integer ciseed


      ranres = 1.0/real(ranbins)

      ranmin = 1.0
      ranmax = 0.0

      
      ciseed = calc_random_seed(-1)
      seed = ciseed

      randist1 = 0.0
      randist2 = 0.0

      lastranval = rancut

      
      call ranini(ciseed)


      do rancnt = 1,ranmaxcnt
         

         ranval = getranf()
         
         ranmin = min(ranmin,ranval)
         ranmax = max(ranmax,ranval)

         in = min(int(ranval/ranres) + 1,ranbins) 

         randist1(in) = randist1(in) + 1.0

         if (lastranval.lt.rancut) then 
            randist2(in) = randist2(in) + 1.0
         endif

         lastranval = ranval


      end do

      ransum1 = sum(randist1)
      ransum2 = sum(randist2)

      write(0,'(a,g20.14)') 'RANMIN=',ranmin
      write(0,'(a,g20.14)') 'RANMAX=',ranmax
      write(0,'(a,g20.14)') 'RANSUM1=',ransum1
      write(0,'(a,g20.14)') 'RANSUM2=',ransum2



      do in = 1,ranbins
         write(6,'(i8,3(1x,g18.10))') in,(in-0.5)*ranres,
     >                        randist1(in)/ransum1,
     >                        randist2(in)/ransum2
      end do


c      do in = 1,ranbins
c         write(7,'(i8,3(1x,g18.10))') in,(in-0.5)*ranres,randist1(in),
c     >                              randist2(in)
c      end do


      return
      end

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
      integer k,i
      integer,allocatable :: temp_seed(:)

      call random_seed(size=k)

      write(0,*) 'Random seed:',iseed
      write(0,*) 'Random seed size = ',k

      allocate(temp_seed(k))

      do i = 1,k
         temp_seed(i) = int(iseed/i)
      end do

      write(0,*) 'Random Seed Used:',(temp_seed(i),i=1,k)

      call random_seed(put=temp_seed)

      deallocate(temp_seed)

c
c      CALL NEWSRAND (ISEED)
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
      implicit none
      DOUBLE PRECISION SEED
      INTEGER NRANDS
      integer j
      REAL CRANDS(NRANDS)
c      real newrand
c      EXTERNAL NEWRAND
      DO 100 J = 1, NRANDS
c
          call random_number(crands(j))
c        
c          CRANDS(J) = NEWRAND()
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
  100 CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE SURAND2 (SEED,NRANDS,RANDS)
      implicit none
      DOUBLE PRECISION SEED
      INTEGER NRANDS
      REAL RANDS,NEWRAND
      EXTERNAL NEWRAND
c
      call random_number(rands)
c
c      RANDS = NEWRAND()
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
      RETURN
      END
c
c
c
      real function getranf ()
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
c      EXTERNAL NEWRAND
c
      call random_number(newrand)
c      GETRANF = NEWRAND()
      GETRANF = NEWRAND
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
      RETURN
      END
c
c
c
      subroutine getran (rand)
      implicit none
      real rand
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
c      REAL ranD,newrand
c      EXTERNAL NEWRAND
c
      call random_number(rand)

c      RAND = NEWRAND()
C
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
      RETURN
      END
c
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
      write(0,'(a,a,a)') 'TIME:',systim,':'
      write(0,'(a,a,a)') 'DATE:',sysdat,':'

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

      character*8 date_f

      call date_and_time(date=date_f)

      sysdat=date_f(5:6)//'/'//date_f(7:8)//'/'//date_f(3:4)


c      CALL DATETEST (SYSDAT)
C       SYSDAT = '04/11/88'

      RETURN
      END
