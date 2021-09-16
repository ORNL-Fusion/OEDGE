

      SUBROUTINE SUCHE_NACHBARN

      USE CTETRA

      IMPLICIT NONE

      TYPE(TET_ELEM), POINTER :: CUR, CUR2
      INTEGER :: ITET,IS,JTET,JS,MINIS,MAXIS,MITIS,MINJS,MAXJS,MITJS,
     .           IP1,IP2,IP3,JP1,JP2,JP3, IC, i, j
      INTEGER :: JP(3),ip(3)
      INTEGER :: ITSIDE(3,4)
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/

      DO ITET=1,NTET      ! FOR ALL TETRAHEDRONS
        DO IS=1,4         ! AND FOR ALL SIDES OF EACH TETRAHEDRON
          IF (NTBAR(IS,ITET) == 0) THEN   ! IF IT HAS NO NEIGHBOR JET
            IP(1)=NTECK(ITSIDE(1,IS),ITET)
            IP(2)=NTECK(ITSIDE(2,IS),ITET)
            IP(3)=NTECK(ITSIDE(3,IS),ITET)

            CUR => COORTET(IP(1))%PTET
            WHLOOP:DO WHILE (ASSOCIATED(CUR))
              JTET = CUR%NOTET
              IF (JTET /= ITET) THEN  ! OMIT TETRAHEDRON ITET
                JSLOOP:DO JS=1,4                      ! CHECK ALL SIDES
                  IF (NTBAR(JS,JTET) == 0) THEN
                    JP(1)=NTECK(ITSIDE(1,JS),JTET)
                    JP(2)=NTECK(ITSIDE(2,JS),JTET)
                    JP(3)=NTECK(ITSIDE(3,JS),JTET)
                     iloop:do i=1,3
                      jloop:do j=i,3
                        if (ip(i) == jp(j)) then
                          ip1=jp(j)
                          jp(j)=jp(i)
                          jp(i)=ip1
                          cycle iloop
                        endif
                      end do jloop
                      cycle jsloop
                      end do iloop
                          NTBAR(IS,ITET) = JTET ! NEIGHBOR FOUND
                          NTSEITE(IS,ITET) = JS
                          NTBAR(JS,JTET) = ITET
                          NTSEITE(JS,JTET) = IS
                          EXIT WHLOOP
                  END IF
                END DO JSLOOP ! JS
              END IF
              CUR => CUR%NEXT_TET
            END DO WHLOOP ! WHILE
          END IF
        END DO
      END DO

!      DO IC=1,NCOOR
!        CUR => COORTET(IC)%PTET
!        DO WHILE (ASSOCIATED(CUR))
!           CUR2 => CUR
!           CUR => CUR%NEXT_TET
!           DEALLOCATE (CUR2)
!        END DO
!        NULLIFY(COORTET(IC)%PTET)
!      END DO
      WRITE (55,'(A,T25,I15)') ' Nachbar-Liste ',MCLSTR*8

      RETURN
      END







