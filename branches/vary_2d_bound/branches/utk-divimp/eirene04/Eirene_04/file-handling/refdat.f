C
C
      SUBROUTINE REFDAT(TMM,TCC,WMM,WCC)
C
C  THIS SUBROUTINE READS REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C    IFILE=1  H ON FE
C    IFILE=2  D ON FE
C    IFILE=3  H ON C
C    IFILE=4  D ON C
C    IFILE=5  HE ON FE
C    IFILE=6  HE ON C
C    IFILE=7  T ON FE
C    IFILE=8  T ON C
C    IFILE=9  D ON W
C    IFILE=10 HE ON W
C    IFILE=11 H ON W
C    IFILE=12 T ON W
C
      USE PRECISION
      USE PARMMOD
      USE CREFMOD
      USE CREF
      USE CSPEI

      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: TMM(*), TCC(*), WMM(*), WCC(*)
      REAL(DP) :: TML(12), TCL(12), WML(12), WCL(12),
     .          FELD(1092)
      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, NRECL, IUN, IFILE, I, J
C
      NFLR=NHD6
      IF (NHD6.GT.12) THEN
        WRITE (6,*) 'STORAGE ERROR. NHD6 MUST BE LESS OR EQUAL 12 '
        WRITE (6,*) 'NHD6= ',NHD6
        CALL EXIT_OWN(1)
      ENDIF
C
      NRECL=1092
      INE=12
      INEM=INE-1
      TML(1)=1.
      TCL(1)=1.
      WML(1)=56.
      WCL(1)=26.
      TML(2)=2.
      TCL(2)=1.
      WML(2)=56.
      WCL(2)=26.
      TML(3)=1.
      TCL(3)=1.
      WML(3)=12.
      WCL(3)=6.
      TML(4)=2.
      TCL(4)=1.
      WML(4)=12.
      WCL(4)=6.
      TML(5)=4.
      TCL(5)=2.
      WML(5)=56.
      WCL(5)=26.
      TML(6)=4.
      TCL(6)=2.
      WML(6)=12.
      WCL(6)=6.
      TML(7)=3.
      TCL(7)=1.
      WML(7)=56.
      WCL(7)=26.
      TML(8)=3.
      TCL(8)=1.
      WML(8)=12.
      WCL(8)=6.
      TML(9)=2.
      TCL(9)=1.
      WML(9)=184.
      WCL(9)=74.
      TML(10)=4.
      TCL(10)=2.
      WML(10)=184.
      WCL(10)=74.
      TML(11)=1.
      TCL(11)=1.
      WML(11)=184.
      WCL(11)=74.
      TML(12)=3.
      TCL(12)=1.
      WML(12)=184.
      WCL(12)=74.
      DO 10 J=1,NHD6
        TMM(J)=TML(J)
        TCC(J)=TCL(J)
        WMM(J)=WML(J)
        WCC(J)=WCL(J)
10    CONTINUE
      ENAR(1)=1.
      ENAR(2)=2.
      ENAR(3)=5.
      ENAR(4)=10.
      ENAR(5)=20.
      ENAR(6)=50.
      ENAR(7)=100.
      ENAR(8)=200.
      ENAR(9)=500.
      ENAR(10)=1000.
      ENAR(11)=2000.
      ENAR(12)=5000.
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      INW=7
      INWM=INW-1
      PID180=ATAN(1.)/45.
      WIAR(1)=COS(0.D0)
      WIAR(2)=COS(30.*PID180)
      WIAR(3)=COS(45.*PID180)
      WIAR(4)=COS(60.*PID180)
      WIAR(5)=COS(70.*PID180)
      WIAR(6)=COS(80.*PID180)
      WIAR(7)=COS(85.*PID180)
      DO 12 I=1,INWM
12      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
C
      INR=5
      INRM=INR-1
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (6,*) 'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EXIT_OWN(1)
      ENDIF
C
      IUN=21
      REWIND IUN
C
660   FORMAT (1X,1P,10E12.4)
661   FORMAT (4E20.12)
      DO 1 IFILE=1,NFLR
        DO 2 I1=1,INE
          I=0
          READ (IUN,661) (FELD(J),J=1,NRECL)
          DO 3 I2=1,INW
            I=I+1
            HFTR0(I1,I2,IFILE)=FELD(I)
3         CONTINUE
          DO 4 I2=1,INW
            DO 4 I3=1,INR
              I=I+1
              HFTR1(I1,I2,I3,IFILE)=FELD(I)
4         CONTINUE
          DO 5 I2=1,INW
            DO 5 I3=1,INR
              DO 5 I4=1,INR
                I=I+1
                HFTR2(I1,I2,I3,I4,IFILE)=FELD(I)
5         CONTINUE
          DO 6 I2=1,INW
            DO 6 I3=1,INR
              DO 6 I4=1,INR
                DO 6 I5=1,INR
                  I=I+1
                  HFTR3(I1,I2,I3,I4,I5,IFILE)=FELD(I)
6         CONTINUE
2       CONTINUE
1     CONTINUE
C
      RETURN
      END
