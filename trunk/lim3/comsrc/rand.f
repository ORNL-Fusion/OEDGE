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
C
C***********************************************************************
C
C     SURAND2
C     ======
C     IBM  : ESSL library routine to generate vector of random numbers
C     CRAY : Replace with calls to RANF generator within a vectorisable
C            loop.
C
       SUBROUTINE SURANDD1 (SEED,NRANDS,RANDS)      
       DOUBLE PRECISION SEED
       INTEGER NRANDS 
       REAL RANDS,NEWRAND
       EXTERNAL NEWRAND
       RANDS = NEWRAND()
c         
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
   100 CONTINUE
       RETURN
       END
C
C***********************************************************************
C
C     SURAN2
C     ======
C     IBM  : ESSL library routine to generate vector of random numbers
C     CRAY : Replace with calls to RANF generator within a vectorisable
C            loop.
C
      SUBROUTINE SURAN2 (SEED,NRANDS,RANDS)      
      DOUBLE PRECISION SEED
      INTEGER NRANDS 
      REAL RANDS,NEWRAND
      EXTERNAL NEWRAND
      RANDS = NEWRAND()
C         
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
  100 CONTINUE
      RETURN
      END
C
C
C***********************************************************************
C
C     SURANDTEST
C     ======
C     IBM  : ESSL library routine to generate vector of random numbers
C     CRAY : Replace with calls to RANF generator within a vectorisable
C            loop.
C
      SUBROUTINE SURANDTEST (SEED,NRANDS,RANDS)      
      DOUBLE PRECISION SEED
      INTEGER NRANDS 
      REAL RANDS,NEWRAND
      EXTERNAL NEWRAND
      RANDS = NEWRAND()
C         
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
  100 CONTINUE
      RETURN
      END
