C
C
      SUBROUTINE RDTRIM
C
C  THIS SUBROUTINE READS SELECTIVELY SOME
C  REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C
      USE PRECISION
      USE PARMMOD
      USE CREFMOD
      USE COMPRT, ONLY: IUNOUT
      USE CREF
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, IUN, IFILE, I, IWWW
C
C
      INE=12
      INW=7
      INR=5
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (iunout,*) 
     .    'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EXIT_OWN(1)
      ENDIF
C
      DO 7 IFILE=1,NFLR
!pb        IUN=20
        IUN=28
        OPEN (UNIT=IUN,FILE=REFFIL(IFILE),ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')
C
        READ (IUN,*)
        READ (IUN,*)
        READ (IUN,*)
        READ (IUN,*)
        DO 2 I1=1,INE
          DO 3 I2=1,INW
            READ (IUN,*)
            READ (IUN,*)
            READ (IUN,*) TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
     .                   enar(i1),wiar(i2),HFTR0(I1,I2,IFILE)
C  FIND NEAREST INTEGER FOR NUCLEAR MASS NUMBER
            IWWW=NINT(WM(IFILE))
            WM(IFILE)=IWWW
            READ (IUN,*)
            READ (IUN,*) (HFTR1(I1,I2,I3,IFILE),I3=1,INR)
            READ (IUN,*)
            DO 5 I3=1,INR
              READ (IUN,*) (HFTR2(I1,I2,I3,I4,IFILE),I4=1,INR)
5           CONTINUE
            READ (IUN,*)
            DO 6 I3=1,INR
            DO 6 I4=1,INR
              READ (IUN,*) (HFTR3(I1,I2,I3,I4,I5,IFILE),I5=1,INR)
6           CONTINUE
3         CONTINUE
2       CONTINUE
        CLOSE (UNIT=IUN)
7     CONTINUE
C
      INEM=INE-1
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      PID180=ATAN(1.)/45.
      DO 12 I=1,INW
12      WIAR(I)=COS(WIAR(I)*PID180)
      INWM=INW-1
      DO 13 I=1,INWM
13      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
      INRM=INR-1
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      RETURN
      END
