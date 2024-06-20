      SUBROUTINE EIRENE_READ_TETRA (CASENAME)
 
!pb 05.12.06: structur coortet is build up from tetrahedra
!pb 07.12.06: set itethand to default value 1
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CTETRA
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      CHARACTER*(*), INTENT(IN) :: CASENAME
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: LL, I, IND, IT, IS, IS1, IER, NRK, ISTS,
     .           ISTMIN, ISTMAX, IC, J, JS, JT, IP1, i1, i2, i3, i4
      INTEGER :: ITSIDE(3,4), IP(3), JP(3)
      TYPE(TET_ELEM), POINTER :: CUR
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
        CALL EIRENE_EXIT_OWN(1)
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
        CALL EIRENE_EXIT_OWN(1)
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
 
      IF (.NOT.ALLOCATED(COORTET)) THEN
        ALLOCATE (COORTET(NCOORD))
        DO I=1,NCOORD
          NULLIFY(COORTET(I)%PTET)
        END DO
      END IF
 
      DO IT=1,NTET
        DO IS=1,4
          IC = NTECK(IS,IT)
          ALLOCATE (CUR)
          CUR%NOTET = IT
          CUR%NEXT_TET => COORTET(IC)%PTET
          COORTET(IC)%PTET => CUR
        ENDDO
      ENDDO
 
!for testing
 
      ntbar = 0
      ntseite = 0
 
      call EIRENE_suche_nachbarn
 
      FILENAME=CASENAME(1:LL) // '.neighbors.out'
      OPEN (UNIT=39,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')
 
      write (39,'(i10)') ntet
 
      DO I=1,NTET
        i1 = 0
        i2 = 0
        i3 = 0
        i4 = 0
        IF (INMTIT(1,I) /= 0) i1 = INMTIT(1,I) - NLIM
        IF (INMTIT(2,I) /= 0) i2 = INMTIT(2,I) - NLIM
        IF (INMTIT(3,I) /= 0) i3 = INMTIT(3,I) - NLIM
        IF (INMTIT(4,I) /= 0) i4 = INMTIT(4,I) - NLIM
        write (39,'(13i10)')
     .        I, NTBAR(1,I), NTSEITE(1,I), i1,
     .           NTBAR(2,I), NTSEITE(2,I), i2,
     .           NTBAR(3,I), NTSEITE(3,I), i3,
     .           NTBAR(4,I), NTSEITE(4,I), i4
      END DO
      close (unit=39)
 
      DO IT=1,NTET
        DO IS=1,4
          IF (NTBAR(IS,IT).EQ.0.AND.INMTIT(IS,IT).EQ.0) THEN
            WRITE (iunout,*)
            WRITE (iunout,*) ' ERROR IN READ_TETRA '
            WRITE (iunout,*) ' OPEN SIDE OF TETRAHEDRON ',IT,' SIDE ',IS
            WRITE (iunout,*) ' NTECK ',(NTECK(ITSIDE(J,IS),IT),J=1,3)
            IER = 3
          ELSE IF (NTBAR(IS,IT) /= 0) THEN
            JT = NTBAR(IS,IT)
            JS =  NTSEITE(IS,IT)
            IF ((NTBAR(JS,JT) /= IT) .OR. (NTSEITE(JS,JT) /= IS)) THEN
              WRITE (iunout,*)
              WRITE (iunout,*) ' ERROR IN READ_TETRA '
              WRITE (iunout,*) ' INCONSISTENCY IN CONNECTION MAP'
              WRITE (iunout,*) ' TETRAHEDRON ',IT,' SIDE ',IS,
     .                         ' HAS NEIGHBOR TET ',JT,' SIDE ',JS
              WRITE (iunout,*) ' BUT '
              WRITE (iunout,*) ' TETRAHEDRON ',JT,' SIDE ',JS,
     .                         ' HAS NEIGHBOR TET',NTBAR(JS,JT),
     .                         ' SIDE ',NTSEITE(JS,JT)
              IER = 4
            END IF
            IP(1)=NTECK(ITSIDE(1,IS),IT)
            IP(2)=NTECK(ITSIDE(2,IS),IT)
            IP(3)=NTECK(ITSIDE(3,IS),IT)
            JP(1)=NTECK(ITSIDE(1,JS),JT)
            JP(2)=NTECK(ITSIDE(2,JS),JT)
            JP(3)=NTECK(ITSIDE(3,JS),JT)
            iloop:do i=1,3
              jloop:do j=i,3
                if (ip(i) == jp(j)) then
                  ip1=jp(j)
                  jp(j)=jp(i)
                  jp(i)=ip1
                  cycle iloop
                endif
              end do jloop
              WRITE (iunout,*)
              write (iunout,*) ' ERROR IN READ_TETRA '
              WRITE (iunout,*) ' INCONSISTENCY IN CONNECTION MAP AND',
     .                         ' ELEMENT DEFINITION '
              WRITE (iunout,*) ' TETRAHEDRON ',IT,' SIDE ',IS,
     .                         ' CONSISTS OF POINTS ', IP
              WRITE (iunout,*) ' NEIGHBOR TET ',JT,' SIDE ',JS,
     .                         ' CONSISTS OF POINTS ', JP
              IER = 5
              exit iloop
            end do iloop
          ENDIF
        ENDDO
      ENDDO
 
      IF (IER /= 0) CALL EIRENE_EXIT_OWN(1)
 
      itethand = 1
 
      RETURN
      END SUBROUTINE EIRENE_READ_TETRA
