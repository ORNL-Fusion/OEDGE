      SUBROUTINE READ_TRIANG (CASENAME)

      USE PRECISION
      USE PARMMOD
      USE CTRIG
      USE CLGIN
      USE CGRID

      IMPLICIT NONE

      CHARACTER*(*), INTENT(IN) :: CASENAME
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: LL, I, IND, IT, IS, IS1, IER, NRK, ISTS, 
     .           ISTMIN, ISTMAX, IC
      
      LL=LEN_TRIM(CASENAME)

      FILENAME=CASENAME(1:LL) // '.npco_char'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '      
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NRK

      IF (NRK /= NRKNOT) THEN
        WRITE (6,*) ' NRKNOT IS WRONG IN EIRENE INPUT FILE'
        WRITE (6,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT_OWN(1)
      END IF
         
      DO I=1,NRKNOT
        READ(30,*) IND, XTRIAN(I), YTRIAN(I)
      END DO
      
!pb      XTRIAN(1:NRKNOT) = XTRIAN(1:NRKNOT) * 100._DP
!pb      YTRIAN(1:NRKNOT) = YTRIAN(1:NRKNOT) * 100._DP

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.elemente'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '      
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NTRII

      IF (NTRII > NR1ST) THEN
        WRITE (6,*) ' NR1ST IS WRONG IN EIRENE INPUT FILE'
        WRITE (6,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT_OWN(1)
      END IF

      DO I=1,NTRII
        READ (30,*) IND, NECKE(1,I), NECKE(2,I), NECKE(3,I)
      END DO

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.neighbors'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '      
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      DO I=1,NTRII
        READ (30,*) IND, NCHBAR(1,I), NSEITE(1,I), INMTI(1,I),
     .                   NCHBAR(2,I), NSEITE(2,I), INMTI(2,I),
     .                   NCHBAR(3,I), NSEITE(3,I), INMTI(3,I),
     .                   IXTRI(I),    IYTRI(I)
        IF (INMTI(1,I) /= 0) INMTI(1,I) = INMTI(1,I) + NLIM
        IF (INMTI(2,I) /= 0) INMTI(2,I) = INMTI(2,I) + NLIM
        IF (INMTI(3,I) /= 0) INMTI(3,I) = INMTI(3,I) + NLIM
      END DO

      CLOSE (UNIT=30)

      IER = 0
      IF ((MAXVAL(NECKE(1:3,1:NTRII)) > NRKNOT) .OR.
     .    (MINVAL(NECKE(1:3,1:NTRII)) <= 0 )) THEN
        WRITE (6,*) ' WRONG COORDINATE NUMBER IS DEFINITION OF',
     .              ' TRIANGLES FOUND '
        IER = 2
      END IF

      DO IT=1,NTRII
        DO IS=1,3
          IF (NCHBAR(IS,IT).EQ.0.AND.INMTI(IS,IT).EQ.0) THEN
            WRITE (6,*) ' ERROR IN READ_TRIANG '
            WRITE (6,*) ' OPEN SIDE OF TRIANGLE ',IT,' SIDE ',IS
            IS1=IS+1
            IF (IS.EQ.3) IS1=1
            WRITE (6,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS,IT)),
     .                                    YTRIAN(NECKE(IS,IT))
            WRITE (6,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS1,IT)),
     .                                    YTRIAN(NECKE(IS1,IT))
            IER = 3
          ENDIF
        ENDDO
      ENDDO
      
      IF (IER /= 0) CALL EXIT_OWN(1)

      ISTMIN=MINVAL(INMTI(1:3,1:NTRII),MASK=(INMTI(1:3,1:NTRII)/=0))
      ISTMAX=MAXVAL(INMTI(1:3,1:NTRII),MASK=(INMTI(1:3,1:NTRII)/=0))

      IC=0
      DO ISTS=ISTMIN,ISTMAX
        DO IS=1,3
          DO IT=1,NTRII
            IF (INMTI(IS,IT) == ISTS) THEN
              IC=IC+1
              INSPAT(IS,IT)=IC
            END IF
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE READ_TRIANG
