      MODULE CTETRA

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CTETRA, DEALLOC_CTETRA, INIT_CTETRA,
     P          TET_ELEM, TET_LISTE

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        RCTET1(:,:), RCTET2(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R  XTETRA(:),  YTETRA(:),  ZTETRA(:),
     R  VTETX(:,:), VTETY(:,:), VTETZ(:,:),
     R  PTETX(:,:), PTETY(:,:), PTETZ(:,:),
     R  XTCEN(:),   YTCEN(:),   ZTCEN(:),
     R  RINCRC(:,:),
     R  RTCEN(:), STCEN(:), TTCEN(:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         ICTET1(:,:), ICTET2(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I  NTECK(:,:),  NTBAR(:,:), NTSEITE(:,:),
     I  INMTIT(:,:), ncltet(:),
     I  NCOOR, NTET, nr1ori, np2ori, nt3ori, ITETHAND, ntet_collaps

      INTEGER, PUBLIC, SAVE :: NCTET1, NCTET2, MCTET1, MCTET2, MCLSTR=0
c slmod begin - tet res
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I INSPATT(:,:)
c slmod end
      TYPE :: TET_ELEM
        INTEGER :: NOTET
        TYPE(TET_ELEM), POINTER :: NEXT_TET
      END TYPE TET_ELEM

      TYPE :: TET_LISTE
        TYPE(TET_ELEM), POINTER :: PTET
      END TYPE TET_LISTE

      TYPE(TET_LISTE), ALLOCATABLE, SAVE, PUBLIC :: COORTET(:)

      CONTAINS


      SUBROUTINE ALLOC_CTETRA

      IF (ALLOCATED(RCTET1)) RETURN

      NCTET1 = (3*6+4*4+3+3)*NTETRA
      NCTET2 = 3*NCOORD
      MCTET1 = (4*4+1)*NTETRA
      MCTET2 = 7

      ALLOCATE (RCTET1(40,NTETRA))
      ALLOCATE (RCTET2(NCTET2))
      ALLOCATE (ICTET1(17,NTETRA))
      ALLOCATE (ICTET2(MCTET2))
c slmod begin - tet res
      ALLOCATE (INSPATT(4,NTETRA))
c slmod end
      WRITE (55,'(A,T25,I15)')
     .       ' CTETRA ',(40*NTETRA+NCTET2)*8 +
     .                  (17*NTETRA+MCTET2)*4

      VTETX  => RCTET1( 1 :  6,:)
      VTETY  => RCTET1( 7 : 12,:)
      VTETZ  => RCTET1(13 : 18,:)
      PTETX  => RCTET1(19 : 22,:)
      PTETY  => RCTET1(23 : 26,:)
      PTETZ  => RCTET1(27 : 30,:)
      XTCEN  => RCTET1(31,:)
      YTCEN  => RCTET1(32,:)
      ZTCEN  => RCTET1(33,:)
      RINCRC => RCTET1(34 : 37,:)
      RTCEN  => RCTET1(38,:)
      STCEN  => RCTET1(39,:)
      TTCEN  => RCTET1(40,:)

      XTETRA => RCTET2(1+0*NCOORD : 1*NCOORD)
      YTETRA => RCTET2(1+1*NCOORD : 2*NCOORD)
      ZTETRA => RCTET2(1+2*NCOORD : 3*NCOORD)

      NTECK   => ICTET1( 1 : 4,:)
      NTBAR   => ICTET1( 5 : 8,:)
      NTSEITE => ICTET1( 9 : 12,:)
      INMTIT  => ICTET1(13 : 16,:)
      ncltet  => ICTET1(17,:)

      NCOOR    => ICTET2(1)
      NTET     => ICTET2(2)
      nr1ori   => ICTET2(3)
      np2ori   => ICTET2(4)
      nt3ori   => ICTET2(5)
      ITETHAND => ICTET2(6)
      ntet_collaps => ICTET2(7)

      CALL INIT_CTETRA
c slmod begin - tet res
      INSPATT = 0
c slmod end
      RETURN
      END SUBROUTINE ALLOC_CTETRA


      SUBROUTINE DEALLOC_CTETRA

      IF (.NOT.ALLOCATED(RCTET1)) RETURN

      DEALLOCATE (RCTET1)
      DEALLOCATE (RCTET2)
      DEALLOCATE (ICTET1)
      DEALLOCATE (ICTET2)
c slmod begin - tet res
      DEALLOCATE (INSPATT)
c slmod end
      RETURN
      END SUBROUTINE DEALLOC_CTETRA


      SUBROUTINE INIT_CTETRA

      RCTET1 = 0._DP
      RCTET2 = 0._DP
      ICTET1 = 0
      ICTET2 = 0

      RETURN
      END SUBROUTINE INIT_CTETRA

      END MODULE CTETRA
