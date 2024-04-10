      SUBROUTINE GZSEED(ISEED)
      INTEGER ISEED
C +---------------------------------------------------------------------
C : Written at: Forschungszentrum KFA Juelich Gmbh
C : Author    : G.Groten, ZAM
C : Update : Busch 16.4.96 REAL*8 -> DOUBLE PRECISION
C : Create the Seed Vector for the KIRKPATRICK-STOLL-GREENWOOD Algorithm
C : If ISEED >= 0, the vector has 532 elements.
C : If ISEED  < 0, the vector has 250 elements.
C : If ISEED is negativ, -ISEED-1 is taken as ISEED.
C : If after this ISEED is 0, 16807 is taken instead of it.
C +---------------------------------------------------------------------
C                          For CMS, MVS, UNICOS, AIX, SunOS, Intel
C The version RASEED KLIBCFT is only for CRAY; but this version is good
C for all machines ( inclusive CRAY).
*
*RASEED/RANGET is a generalized feedback shift register random number
*generator (GFSR-generator) working on 47 bits. Uses vector-code on IBM.
*
*  RASEED(ISEED)  starting subroutine
*                 INTEGER ISEED : Starting Number.
*
*  RANGEN(R,NR)   subroutine to get NR random numbers R between 0 and 1.
*                 INTEGER NR
*                 DOUBLE PRECISION R(NMAX)  with NR <= NMAX.
*
*Each bit has period of 2**532-1 for ISEED >= 0, of 2**250-1 for
*ISEED < 0. If the user gives ISEED=0 or ISEED=-1, an other very good
*starting value is choosen. You can use RANGEN again and again.
*
*You can get the source program with "GIME VMCRAYPR".
*C-------------------------- short example -----------------------------
*      PROGRAM TESRA
*      INTEGER I,ISEED,N,M
*      PARAMETER(M=10,N=100*M+1)
*      DOUBLE PRECISION X(N)
*      ISEED =  0
*      CALL GZSEED(ISEED)
*      CALL GZNGEN(X,N)
*      DO 1 I=1,N,100
*         WRITE(*,'(1X,I4,1X,E22.14,1X,Z16)') I,X(I),X(I)
*    1 CONTINUE
*      END
C-----------------------------------------------------------------------
      INTEGER IEXP2,IPRIME,MULTOR,ONES,DIAG,I,LAST0,IDIV,IBIN,IN1,IN2
      PARAMETER(IEXP2=19, IPRIME=2**IEXP2-1, MULTOR=1387)
      PARAMETER ( IBIN= 24, IDIV= 2**IBIN, LAST0= IDIV-2)
      INTEGER MMM,IDIST,MM1,MM,MM0,M,JSEED,KK
      COMMON /GZCOMM/ KK,MMM,MM1,MM,MM0,M(1064)
CDEC$ PSECT /GZCOMM/ NOSHR
       SAVE /GZCOMM/

      KK = 1
      IF ( ISEED.GE.0 ) THEN
         MMM = 532
         MM = 37
         IDIST = 10
         IF (ISEED .EQ. 0) THEN
            JSEED = 16807
         ELSE
            JSEED = ISEED
         ENDIF
      ELSE
         MMM = 250
         MM = 103
         IDIST = 5
         JSEED = ISEED+1
         IF (JSEED .EQ. 0) THEN
            JSEED = 16807
         ELSE
            JSEED = -JSEED
         ENDIF
      ENDIF
      JSEED = MOD( JSEED , IPRIME )
      MM0 = MMM-MM
      MM1 = MMM+1

      DO 1 I=1,MMM
C------- ersten 24 Bits von 48: Teilen durch QDIV
         JSEED=MOD(JSEED*MULTOR,IPRIME)
         IN1 = JSEED/128
         JSEED=MOD(JSEED*MULTOR,IPRIME)
         IN2 = JSEED/128
         M(I) = IN1*4096 + IN2
C------- die letzten 24 Bits von 48 zu bekommen: Modulo
         JSEED=MOD(JSEED*MULTOR,IPRIME)
         IN1 = JSEED/128
         JSEED=MOD(JSEED*MULTOR,IPRIME)
         IN2 = JSEED/128
         M(I+MMM) = IAND(IN1*4096+IN2,LAST0)
    1 CONTINUE

      ONES = IDIV-1
      DIAG = IDIV/2

      DO 2 I=IDIST,MIN(IBIN,47)*IDIST,IDIST
         M(I) = IOR (M(I),DIAG)
         M(I) = IAND (M(I),ONES)
         ONES = ISHFT(ONES,-1)
         DIAG = ISHFT(DIAG,-1)
    2 CONTINUE

      ONES = IDIV-1
      DIAG = IDIV/2

      DO 3 I=(IBIN+1)*IDIST,47*IDIST,IDIST
         M(I) = 0
         M(I+MMM) = IOR (M(I+MMM),DIAG)
         M(I+MMM) = IAND (M(I+MMM),ONES)
         ONES = ISHFT(ONES,-1)
         DIAG = ISHFT(DIAG,-1)
    3 CONTINUE

      END
CPROCESS DIR('DIR:')
      SUBROUTINE GZNGEN(X,N)
      INTEGER N
      DOUBLE PRECISION X(N)

C  +---------- IBM / 370 -- und RISC 6000 --- SUN ----( CRAY )-- Intel
C  : Random Number Generator (KIRKPATRICK-STOLL)
C  : Es entstehen genau dieselben Zufallszahlen wie bei CRAY.
C  : X muss als DOUBLE PRECISION vereinbart sein - stimmt auch bei CRAY.
C  : Enthalten sind:
C  : Compiler Direktiven zum Vektorisieren fuer CRAY und IBM/370.
C  :
C  : Copyright : Forschungszentrum KFA Juelich
C  : Author: G.Groten
C  +--------------------------------------------------------------------
      DOUBLE PRECISION FA1,FA2
      INTEGER M,K,L,I,KK,IBIN,MMM,MM1,MM,MM0
      PARAMETER ( IBIN = 24, FA1 = 2**IBIN, FA2 = 1D0/FA1**2)
      COMMON /GZCOMM/ KK,MMM,MM1,MM,MM0,M(1064)
CDEC$ PSECT /GZCOMM/ NOSHR
      SAVE /GZCOMM/
C     +-----------------------------------------------------------------
C     : NEXT STATEMENT ONLY TO LOAD INTO REGISTER.
C     +-----------------------------------------------------------------
      K = KK
C     ------------------------------------------------------------------
      I = 0
      L = N+MMM
C     +-----------------------------------------------------------------
C     : BEGIN ROTATIONAL LOOP
C     +-----------------------------------------------------------------
    5 L = L-(MM1-K)
CDIR: IGNORE RECRDEPS
CDIR$ IVDEP
      DO K = K, MIN(L,MM0)
C        +--------------------------------------------------------------
C        : THE M'S TOGETHER ARE EVEN INTEGER RANDOM NUMBERS 0<=M<2**48
C        +--------------------------------------------------------------
         M(K) = IEOR(M(K),M(K+MM))
         M(K+MMM) = IEOR(M(K+MMM),M(K+MMM+MM))
         I = I+1
         X(I)= (M(K)*FA1+(M(K+MMM)+1))*FA2
      end do
C     ------------------------------------------------------------------
CDIR: IGNORE RECRDEPS
CDIR$ IVDEP
      DO K=K,MIN(L,MMM)
         M(K) = IEOR(M(K),M(K-MM0))
         M(K+MMM) = IEOR(M(K+MMM),M(K+MMM-MM0))
         I = I+1
         X(I) = (M(K)*FA1+(M(K+MMM)+1))*FA2
      end do
C     ------------------------------------------------------------------
      IF (K.EQ.MM1) THEN
         K = 1
         GOTO 5
         ENDIF
C     +-----------------------------------------------------------------
C     : END ROTATIONAL LOOP
C     +-----------------------------------------------------------------
      KK = K
      END
