C
C
      SUBROUTINE EIRENE_SNEIGH
C
C  DETERMINE THE INDICES OF THE FOUR NEIGHBORING CELLS OF
C  AN EIRENE CELL    : NGHPOL
C  AN EIRENE SURFCASE: NGHPLS
C  SIDE NUMBERING:
C
C       (IR+1,IP+1)          (IR,IP+1)
C                       2
C                 +-----------+
C                 |           |
C                 |           |
C               3 |           | 1
C                 |           |
C                 |           |
C                 +-----------+
C                       4
C         (IR+1,IP)          (IR,IP)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_module_avltree
 
      IMPLICIT NONE
 
      INTEGER :: IR, IP, IPART, JP, K, IC, IN
      TYPE(CELL_ELEM), POINTER :: CUR
      type(TAVLTree), pointer :: baum
      logical :: inserted

      IF (LEVGEO == 1) THEN

        IC = 0
        DO IP=1,NP2ND
          DO IR=1,NR1ST
            IC = IC + 1
            INDPOINT(IR,IP) = IC
          ENDDO
        ENDDO
         
      
      ELSE IF ((LEVGEO == 2) .OR. (LEVGEO == 3)) THEN
 
        DO IR=1,NR1ST
          DO IP=1,NP2ND
            NGHPOL(1,IR,IP)=0
            NGHPOL(2,IR,IP)=0
            NGHPOL(3,IR,IP)=0
            NGHPOL(4,IR,IP)=0
            NGHPLS(1,IR,IP)=0
            NGHPLS(2,IR,IP)=0
            NGHPLS(3,IR,IP)=0
            NGHPLS(4,IR,IP)=0
          ENDDO
        ENDDO
 
        DO IR=1,NR1STM
          DO K=1,NPPLG
            DO IP=NPOINT(1,K),NPOINT(2,K)
C  RADIAL NEIGHBORS
              NGHPOL(1,IR,IP)=IR-1
              NGHPOL(3,IR,IP)=MOD(IR+1,NR1ST)
              NGHPLS(1,IR,IP)=IR-1
              NGHPLS(3,IR,IP)=IR
C  POLOIDAL NEIGHBORS
              IF (IP.GT.NPOINT(1,K).AND.IP.LT.NPOINT(2,K)-1) THEN
C  INNER CELL
                NGHPOL(4,IR,IP)=IP-1
                NGHPOL(2,IR,IP)=IP+1
                NGHPLS(4,IR,IP)=IP-1
                NGHPLS(2,IR,IP)=IP
              ELSE IF (IP.EQ.NPOINT(1,K)) THEN
C  FIRST CELL OR SURFACE IN A PART OF A POLYGON
                NGHPOL(2,IR,IP)=IP+1
                NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
                DO IPART=1,NPPLG
                  JP=NPOINT(2,IPART)
                  IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .                ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                    NGHPOL(4,IR,IP)=JP-1
                    NGHPLS(4,IR,IP)=JP-1
                    GOTO 1
                  ENDIF
                ENDDO
              ELSEIF (IP.EQ.NPOINT(2,K)-1) THEN
C  LAST CELL IN A PART OF A POLYGON
                NGHPOL(4,IR,IP)=IP-1
                NGHPLS(4,IR,IP)=IP-1
                NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
                DO IPART=1,NPPLG
                  JP=NPOINT(1,IPART)
                  IF (((XPOL(IR,IP+1)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP+1)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .                ((XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                    NGHPOL(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ENDDO
              ELSEIF (IP.EQ.NPOINT(2,K)) THEN
C  LAST SURFACE IN A PART OF A POLYGON
                NGHPLS(4,IR,IP)=IP-1
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
                DO IPART=1,NPPLG
                  JP=NPOINT(1,IPART)
                  IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .                ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                    NGHPLS(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ENDDO
              ENDIF
1           ENDDO
          ENDDO
        ENDDO
 
        baum => EIRENE_NewTree()
        NNODES=0
        DO IR=1,NR1STM
          DO K=1,NPPLG
            DO IP=NPOINT(1,K),NPOINT(2,K)-1
              IN=IR+(IP-1)*NR1ST
! IR, IP
              IC=NNODES+1
              inserted=.false.
              call EIRENE_insert (baum, xpol(ir,ip), ypol(ir,ip), 
     .                            0._DP, 1._DP, ic, inserted)
              IF (INSERTED) THEN
                NNODES=NNODES+1
                XPOINT(NNODES)=XPOL(IR,IP)
                YPOINT(NNODES)=YPOL(IR,IP)
              END IF
              ALLOCATE (CUR)
              CUR%NOCELL = IN
              CUR%NEXT_CELL => COORCELL(IC)%PCELL
              COORCELL(IC)%PCELL => CUR
              INDPOINT(IR,IP) = IC
! IR+1, IP
              IC=NNODES+1
              inserted=.false.
              call EIRENE_insert (baum, xpol(ir+1,ip), ypol(ir+1,ip),
     .                            0._DP, 1._DP, ic, inserted)
              IF (INSERTED) THEN
                NNODES=NNODES+1
                XPOINT(NNODES)=XPOL(IR+1,IP)
                YPOINT(NNODES)=YPOL(IR+1,IP)
              END IF
              ALLOCATE (CUR)
              CUR%NOCELL = IN
              CUR%NEXT_CELL => COORCELL(IC)%PCELL
              COORCELL(IC)%PCELL => CUR
              INDPOINT(IR+1,IP) = IC
! IR+1, IP+1
              IC=NNODES+1
              inserted=.false.
              call EIRENE_insert (baum, xpol(ir+1,ip+1),ypol(ir+1,ip+1),
     .                            0._DP, 1._DP, ic, inserted)
              IF (INSERTED) THEN
                NNODES=NNODES+1
                XPOINT(NNODES)=XPOL(IR+1,IP+1)
                YPOINT(NNODES)=YPOL(IR+1,IP+1)
              END IF
              ALLOCATE (CUR)
              CUR%NOCELL = IN
              CUR%NEXT_CELL => COORCELL(IC)%PCELL
              COORCELL(IC)%PCELL => CUR
              INDPOINT(IR+1,IP+1) = IC
! IR, IP+1
              IC=NNODES+1
              inserted=.false.
              call EIRENE_insert (baum, xpol(ir,ip+1), ypol(ir,ip+1),
     .                            0._DP, 1._DP, ic, inserted)
              IF (INSERTED) THEN
                NNODES=NNODES+1
                XPOINT(NNODES)=XPOL(IR,IP+1)
                YPOINT(NNODES)=YPOL(IR,IP+1)
              END IF
              ALLOCATE (CUR)
              CUR%NOCELL = IN
              CUR%NEXT_CELL => COORCELL(IC)%PCELL
              COORCELL(IC)%PCELL => CUR
              INDPOINT(IR,IP+1) = IC
            ENDDO
          ENDDO
        ENDDO
        call EIRENE_DestroyTree(baum)
 
C       WRITE (iunout,*) '  IR    IP    S1    S2    S3    S4'
C       DO IR=1,NR1STM
C         DO IP=1,NP2NDM
C           WRITE (iunout,'(6I6)') IR,IP,(NGHPOL(K,IR,IP),K=1,4)
C           WRITE (iunout,'(6I6)') IR,IP,(NGHPLS(K,IR,IP),K=1,4)
C         ENDDO
C       ENDDO
 
      END IF
  
      END
 
 
 
 
