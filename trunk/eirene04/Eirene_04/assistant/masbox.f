C
C
C*DK MASBOX
      SUBROUTINE MASBOX (STRING)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: STRING
      CHARACTER(80) :: AST(5)
      INTEGER :: I, J, ILEN

      DO I=1,5
        DO J=1,80
          AST(I)(J:J)=' '
        ENDDO
      ENDDO

      ILEN=LEN(STRING)
      DO I=1,ILEN+6
        AST(1)(I:I)='*'
        AST(5)(I:I)='*'
      ENDDO
      AST(2)(1:1)='*'
      AST(2)(ILEN+6:ILEN+6)='*'
      AST(3)(1:1)='*'
      AST(3)(ILEN+6:ILEN+6)='*'
      AST(4)(1:1)='*'
      AST(4)(ILEN+6:ILEN+6)='*'
C
      DO I=1,ILEN
        AST(3)(I+3:I+3)=STRING(I:I)
      ENDDO
C
      WRITE (6,'(//1X)')
      DO I=1,5
        WRITE (6,'(A80)') AST(I)
      ENDDO
      WRITE (6,'(//1X)')
C
      RETURN
      END
