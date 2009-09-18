      SUBROUTINE READ_TETRA (CASENAME)

      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE CLGIN
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      CHARACTER*(*), INTENT(IN) :: CASENAME
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: LL, I, IND, IT, IS, IS1, IER, NRK, ISTS,
     .           ISTMIN, ISTMAX, IC, J
      INTEGER :: ITSIDE(3,4)
C
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/

      LL=LEN_TRIM(CASENAME)

      FILENAME=CASENAME(1:LL) // '.npco_char'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NRK

      IF (NRK /= NCOORD) THEN
        WRITE (iunout,*) ' NCOORD IS WRONG IN EIRENE INPUT FILE'
        WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT(1)
      END IF

      DO I=1,NCOORD
        READ(30,*) IND, XTETRA(I), YTETRA(I), ZTETRA(I)
      END DO

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.elemente'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NTET

      IF (NTET > NR1ST) THEN
        WRITE (iunout,*) ' NR1ST IS WRONG IN EIRENE INPUT FILE'
        WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT(1)
      END IF

      DO I=1,NTET
        READ (30,*) IND, NTECK(1,I), NTECK(2,I), NTECK(3,I), NTECK(4,I)
      END DO

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.neighbors'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      DO I=1,NTET
        READ (30,*) IND, NTBAR(1,I), NTSEITE(1,I), INMTIT(1,I),
     .                   NTBAR(2,I), NTSEITE(2,I), INMTIT(2,I),
     .                   NTBAR(3,I), NTSEITE(3,I), INMTIT(3,I),
     .                   NTBAR(4,I), NTSEITE(4,I), INMTIT(4,I)
        IF (INMTIT(1,I) /= 0) INMTIT(1,I) = INMTIT(1,I) + NLIM
        IF (INMTIT(2,I) /= 0) INMTIT(2,I) = INMTIT(2,I) + NLIM
        IF (INMTIT(3,I) /= 0) INMTIT(3,I) = INMTIT(3,I) + NLIM
        IF (INMTIT(4,I) /= 0) INMTIT(4,I) = INMTIT(4,I) + NLIM
      END DO

      CLOSE (UNIT=30)

      IER = 0
      IF ((MAXVAL(NTECK(1:4,1:NTET)) > NCOORD) .OR.
     .    (MINVAL(NTECK(1:4,1:NTET)) <= 0 )) THEN
        WRITE (iunout,*) ' WRONG COORDINATE NUMBER IS DEFINITION OF',
     .              ' TRIANGLES FOUND '
        IER = 2
      END IF

      DO IT=1,NTET
        DO IS=1,4
          IF (NTBAR(IS,IT).EQ.0.AND.INMTIT(IS,IT).EQ.0) THEN
            WRITE (iunout,*) ' ERROR IN READ_TETRA '
            WRITE (iunout,*) ' OPEN SIDE OF TETRAHEDRON ',IT,' SIDE ',IS
            WRITE (iunout,*) ' NTECK ',(NTECK(ITSIDE(J,IS),IT),J=1,3)
            IER = 3
          ENDIF
        ENDDO
      ENDDO

      IF (IER /= 0) CALL EXIT(1)

      RETURN
      END SUBROUTINE READ_TETRA
