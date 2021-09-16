      SUBROUTINE RPSCUT (AORIG,VALUES)

      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE CTRIG
      USE CCONA
      USE CPLOT
      USE CGEOM
      USE COMUSR
      USE MODULE_AVLTREE

      IMPLICIT NONE

      type(TAVLTree), pointer, save :: baum

      TYPE :: TRIANODE
        INTEGER :: NODES(3)
        INTEGER :: NOTET
        TYPE(TRIANODE), POINTER :: NEXT
      END TYPE TRIANODE

      TYPE :: TRELEM
        INTEGER :: NOTRIAN
        TYPE(TRELEM), POINTER :: NEXT_TRIAN
      END TYPE TRELEM

      TYPE :: TRELEMP
        TYPE(TRELEM), POINTER :: PTREL
      END TYPE TRELEMP
C
      INTERFACE 
        SUBROUTINE INSERT_POINT(baum,X,Y,Z,IC1)
          USE PRECISION
          USE module_avltree
          type(TAVLTree), pointer :: baum
          REAL(DP), INTENT(IN)      :: X, Y, Z
          INTEGER, INTENT(INOUT)  :: IC1
        END SUBROUTINE INSERT_POINT
      END INTERFACE
C
      TYPE(TRIANODE), POINTER :: TRIALIST, CUR
      TYPE(TET_ELEM), POINTER :: CURE
      TYPE(TRELEMP), ALLOCATABLE :: NODELIST(:)
      TYPE(TRELEM), POINTER :: PTRI, PTRI2
      
      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(OUT) :: VALUES(*)
      REAL(DP) :: TET(4,3), CTPNTS(6,3), THIRD, DIST12, DIST23, DIST13,
     .            DIST34, DIST14, DIST, SPP(3), DMAT(6,6), DMIN, DMAX

      INTEGER :: ITET, NCTPNT, IP1, IP2, IP3, IP4, ITRI, NBACK, NKANT,
     .           IC, NCO, IK1, IK2, IKANTE, IP, JC, MINPOS(2), I, J,
     .           IMIN, IMAX, IMID, JMIN, JMAX, JMID
      INTEGER, ALLOCATABLE :: KANTEN(:,:), KANTET(:,:), NUMCO(:)
      INTEGER, ALLOCATABLE, SAVE :: NUMTET(:)
      INTEGER :: KMAT(4,4) = RESHAPE((/ 0, 1, 2, 3,
     .                                  1, 0, 4, 5,
     .                                  2, 4, 0, 6,
     .                                  3, 5, 6, 0 /), (/4, 4/))
      INTEGER, SAVE :: IC1, IC2, IC3
      INTEGER, SAVE :: IFIRST=0
      LOGICAL :: INSERTED, EXI_POINT(6), LINSERT
      LOGICAL, ALLOCATABLE :: VISITED(:),LKANTE(:)

      IF (IFIRST == 0) THEN
        IFIRST = 1
! DETERMINE EDGES OF TETRAHEDRONS
        ALLOCATE (KANTET(6,NTET))
        ALLOCATE (KANTEN(2,6*NTET))
        ALLOCATE (NUMCO(NCOORD))
        ALLOCATE (VISITED(NCOORD))
        KANTEN=0
        KANTET=0
        VISITED=.FALSE.
        NKANT=0
        DO IC=1,NCOORD ! FOR EACH COORDINATE IC THE TET'S BELONGING TO IT ARE STORED
          NCO = 0
          CURE => COORTET(IC)%PTET
          DO WHILE (ASSOCIATED(CURE))
            ITET=CURE%NOTET
            DO IP=1,4 ! FIND ALL COORDINATES CONNECTED TO IC
              IF ((NTECK(IP,ITET) .NE. IC) .AND. 
     .            .NOT.VISITED(NTECK(IP,ITET))) THEN
                NCO=NCO+1
                NUMCO(NCO)=NTECK(IP,ITET)
                VISITED(NTECK(IP,ITET))=.TRUE.
              END IF 
            END DO
            CURE => CURE%NEXT_TET 
          END DO ! WHILE
! CHECK ALL COORDINATES CONNECTED TO IC AND NUMBER TO EDGES
          DO JC=1,NCO
            IKANTE=0
            CURE => COORTET(IC)%PTET
            DO WHILE (ASSOCIATED(CURE))
              ITET=CURE%NOTET
              IK1=0
              IK2=0
              DO IP=1,4
                IF (NTECK(IP,ITET) == IC) IK1=IP
                IF (NTECK(IP,ITET) == NUMCO(JC)) IK2=IP
              END DO
              IF (IK2 .NE. 0) THEN ! COORDINATE JC BELONGS TO TETRAHEDRON ITET
                IF (KANTET(KMAT(IK1,IK2),ITET) == 0) THEN 
! THE EDGE HASN'T BEEN VISITED YET
                  IF (IKANTE == 0) THEN
! EDGE (IK1,IK2) HAS NO NUMBER, ADD IT TO THE LIST OF EDGES                     
                    NKANT=NKANT+1
                    KANTEN(1:2,NKANT)=(/IC,NUMCO(JC)/)
                    KANTET(KMAT(IK1,IK2),ITET)=NKANT
                    IKANTE=NKANT
                  ELSE
! EDGE (IK1,IK2) HAS NUMBER IKANTE, MARK EDGE OF TETRAHEDRON
                    KANTET(KMAT(IK1,IK2),ITET)=IKANTE
                  ENDIF
                ELSE
! THE EDGE HAS BEEN VISITED BEFORE, TAKE THE NUMBER
                  IKANTE=KANTET(KMAT(IK1,IK2),ITET)
                END IF
              END IF
              CURE => CURE%NEXT_TET 
              VISITED(NUMCO(JC))=.FALSE.
            END DO ! WHILE    
          END DO ! JC
        END DO ! IC

!        IF (MINVAL(KANTET) <= 0) THEN
!          WRITE (6,*) ' ERROR IN RPSCUT '
!          WRITE (6,*) ' EDGE OF TETRAHEDRON HAS NO NUMBER '
!        END IF
        
        DEALLOCATE (VISITED)

! DETERMINE WHICH COORDINATES SHOULD BE USED AS TRIANGLE VERTICES
        IF (ABS(CUTPLANE(2)) > EPS10) THEN
          IC1=2
          IC2=3
          IC3=1
        ELSEIF (ABS(CUTPLANE(3)) > EPS10) THEN
          IC1=1
          IC2=3
          IC3=2
        ELSEIF (ABS(CUTPLANE(4)) > EPS10) THEN
          IC1=1
          IC2=2
          IC3=3
        END IF

        baum => NewTree()
        NULLIFY(TRIALIST)
        NRKNOT = 0
        NTRII = 0

        ALLOCATE (LKANTE(NKANT))
        LKANTE=.TRUE.

        DO ITET=1,NTET

          IF (VOL(ITET) < EPS10) CYCLE

          TET(1,1:3) = (/ XTETRA(NTECK(1,ITET)),
     .                    YTETRA(NTECK(1,ITET)),
     .                    ZTETRA(NTECK(1,ITET)) /)
          TET(2,1:3) = (/ XTETRA(NTECK(2,ITET)),
     .                    YTETRA(NTECK(2,ITET)),
     .                    ZTETRA(NTECK(2,ITET)) /)
          TET(3,1:3) = (/ XTETRA(NTECK(3,ITET)),
     .                    YTETRA(NTECK(3,ITET)),
     .                    ZTETRA(NTECK(3,ITET)) /)
          TET(4,1:3) = (/ XTETRA(NTECK(4,ITET)),
     .                    YTETRA(NTECK(4,ITET)),
     .                    ZTETRA(NTECK(4,ITET)) /)
!  CUT TETRAHEDRON TET WITH PLANE CUTPLANE
          NCTPNT=0
          EXI_POINT=.FALSE.
          if (itet == 100988) then
            write (6,*) ' ogottogott '
          end if
          DO I=1,3
            DO J=I+1,4
              IKANTE=KMAT(I,J)
!pb              IF (LKANTE(KANTET(IKANTE,ITET))) THEN
                CALL SCHNITT_GER_EB(SPP, EXI_POINT(IKANTE),
     .                          CUTPLANE, TET(I,:), TET(J,:))
                LKANTE(KANTET(IKANTE,ITET))=EXI_POINT(IKANTE)
                IF (EXI_POINT(IKANTE)) THEN
                  NCTPNT=NCTPNT+1
                  CTPNTS(IKANTE,:)=SPP
                END IF
!pb              END IF
            END DO
          END DO
          IF (NCTPNT < 3) CYCLE
          DO 
            DMAT=HUGE(1._DP)
            DO I=1,6
              DO J=1,6
                IF (EXI_POINT(I) .AND. EXI_POINT(J) .AND. (I.NE.J)) 
     .            DMAT(I,J)=SQRT(SUM((CTPNTS(I,:)-CTPNTS(J,:))**2))
              END DO
            END DO
            DMIN=MINVAL(DMAT)
            DMAX=MAXVAL(DMAT,MASK=DMAT<HUGE(1._DP))
            MINPOS=MINLOC(DMAT)
            IF (NCTPNT > 4) THEN
              EXI_POINT(MINPOS(2))=.FALSE.
              LKANTE(KANTET(MINPOS(2),ITET))=.FALSE.
              NCTPNT=NCTPNT-1
            ELSE IF (NCTPNT == 4) THEN
              IF (DMIN/DMAX < 1.E-2) THEN
                EXI_POINT(MINPOS(2))=.FALSE.
                LKANTE(KANTET(MINPOS(2),ITET))=.FALSE.
                NCTPNT=NCTPNT-1
              ELSE
                EXIT
              END IF
            ELSE ! nctpnt < 4
              IF (DMIN < EPS5) THEN
                EXI_POINT(MINPOS(2))=.FALSE.
                LKANTE(KANTET(MINPOS(2),ITET))=.FALSE.
                NCTPNT=NCTPNT-1
              END IF
!pb              IF (NCTPNT < 3) EXIT
              EXIT
            END IF
          END DO 
          IF (NCTPNT < 3) CYCLE
          IF (DMAX < EPS5) CYCLE
          NCO=1
          DO I=1,6
            IF (EXI_POINT(I)) THEN
              IF (NCO < I) THEN
                CTPNTS(NCO,:)=CTPNTS(I,:)
              END IF
              NCO=NCO+1
            END IF
          END DO
          IF (NCTPNT == 4) call sort_ueberpruef(CTPNTS(1:4,:))
!          CALL TETRA_SCHNITT (CUTPLANE, TET, NCTPNT, CTPNTS)
!pb          IF (NCTPNT > 2) THEN
            DIST12=SQRT((CTPNTS(1,IC1)-CTPNTS(2,IC1))**2 +
     .                  (CTPNTS(1,IC2)-CTPNTS(2,IC2))**2 +
     .                  (CTPNTS(1,IC3)-CTPNTS(2,IC3))**2)
            DIST13=SQRT((CTPNTS(1,IC1)-CTPNTS(3,IC1))**2 +
     .                  (CTPNTS(1,IC2)-CTPNTS(3,IC2))**2 +
     .                  (CTPNTS(1,IC3)-CTPNTS(3,IC3))**2)
            DIST23=SQRT((CTPNTS(2,IC1)-CTPNTS(3,IC1))**2 +
     .                  (CTPNTS(2,IC2)-CTPNTS(3,IC2))**2 +
     .                  (CTPNTS(2,IC3)-CTPNTS(3,IC3))**2)
            DIST=MIN(DIST12,DIST13,DIST23)
!  INTERSECTION WITH TETRAHEDRON IS A TRIANGLE AT LEAST              
            IP1=NRKNOT+1
            inserted=.false.
            CALL INSERT(baum,CTPNTS(1,IC1),CTPNTS(1,IC2),
     .                       CTPNTS(1,IC3),DIST,IP1,INSERTED)
            IF (INSERTED) NRKNOT=NRKNOT+1
            IP2=NRKNOT+1
            inserted=.false.
            CALL INSERT(baum,CTPNTS(2,IC1),CTPNTS(2,IC2),
     .                       CTPNTS(2,IC3),DIST,IP2,INSERTED)
            IF (INSERTED) NRKNOT=NRKNOT+1
            IP3=NRKNOT+1
            inserted=.false.
            CALL INSERT(baum,CTPNTS(3,IC1),CTPNTS(3,IC2),
     .                       CTPNTS(3,IC3),DIST,IP3,INSERTED)
            IF (INSERTED) NRKNOT=NRKNOT+1
            IP4=0

!  STORE TRIANGLE FOR FURTHER USE
            ALLOCATE(CUR)
            CUR%NODES = (/ IP1, IP2, IP3 /)
            CUR%NOTET = ITET
            CUR%NEXT => TRIALIST
            TRIALIST => CUR
            NTRII = NTRII + 1

            IF (NCTPNT == 4) THEN
              DIST34=SQRT((CTPNTS(3,IC1)-CTPNTS(4,IC1))**2 +
     .                    (CTPNTS(3,IC2)-CTPNTS(4,IC2))**2 +
     .                    (CTPNTS(3,IC3)-CTPNTS(4,IC3))**2)
              DIST14=SQRT((CTPNTS(1,IC1)-CTPNTS(4,IC1))**2 +
     .                    (CTPNTS(1,IC2)-CTPNTS(4,IC2))**2 +
     .                    (CTPNTS(1,IC3)-CTPNTS(4,IC3))**2)
              DIST=MIN(DIST34,DIST14,DIST13)
!  INTERSECTION WITH TETRAHEDRON IS A QUADRANGLE, FIRST PART OF 
!  QUADRANGLE IS ALEADY STORED, NOW STORE SECOND PART               
              IP4=NRKNOT+1
              inserted=.false.
              CALL INSERT(baum,CTPNTS(4,IC1),CTPNTS(4,IC2),
     .                         CTPNTS(4,IC3),DIST,IP4,INSERTED)
              IF (INSERTED) NRKNOT=NRKNOT+1
              ALLOCATE(CUR)
              CUR%NODES = (/ IP3, IP4, IP1 /)
              CUR%NOTET = ITET
              CUR%NEXT => TRIALIST
              TRIALIST => CUR
              NTRII = NTRII + 1
            END IF

!pb          END IF
        END DO
        
        NKNOT = NRKNOT

!  CLEAR STORAGE FOR TRIANGLES   
        IF (ALLOCATED(XTRIAN)) THEN
          DEALLOCATE (XTRIAN)
          DEALLOCATE (YTRIAN)
          DEALLOCATE (NECKE)
        END IF
        
        ALLOCATE (XTRIAN(NRKNOT))
        ALLOCATE (YTRIAN(NRKNOT))
        ALLOCATE (NECKE(3,NTRII))
        ALLOCATE (NUMTET(NTRII))
        ALLOCATE (NODELIST(NRKNOT))
        
!  WALK THROUGH COORDINATE TREE AND STORE THE COORDINATES IN XTRIAN, YTRIAN
        CALL TREE_TRAVERSE (BAUM%ROOT)
        
        call DestroyTree(baum)

        DO I=1,NRKNOT
          NULLIFY(NODELIST(I)%PTREL)
        END DO
        
        THIRD = 1._DP / 3._DP
!pb        NBACK=NTRII
        NBACK=0
        DO WHILE (ASSOCIATED(TRIALIST))
          CUR => TRIALIST
          TRIALIST => CUR%NEXT
          IMIN = MINVAL(CUR%NODES)
          IMAX = MAXVAL(CUR%NODES)
          IMID = SUM(CUR%NODES)-IMIN-IMAX
          PTRI => NODELIST(IMIN)%PTREL
          LINSERT=.TRUE.
          DO WHILE (ASSOCIATED(PTRI)) 
            JMIN = MINVAL(NECKE(1:3,PTRI%NOTRIAN)) 
            JMAX = MAXVAL(NECKE(1:3,PTRI%NOTRIAN)) 
            JMID = SUM(NECKE(1:3,PTRI%NOTRIAN))-JMIN-JMAX
            IF ((IMIN == JMIN) .AND. (IMAX == JMAX) .AND. 
     .          (IMID == JMID)) THEN
              LINSERT=.FALSE.
              EXIT
            END IF
            PTRI => PTRI%NEXT_TRIAN 
          END DO
          IF (LINSERT) THEN
            NBACK = NBACK + 1
            NECKE(1:3,NBACK) = CUR%NODES
            NUMTET(NBACK) = CUR%NOTET
            XCOM(NBACK) = THIRD * ( XTRIAN(NECKE(1,NBACK)) + 
     .                              XTRIAN(NECKE(2,NBACK)) +   
     .                              XTRIAN(NECKE(3,NBACK)) )   
            YCOM(NBACK) = THIRD * ( YTRIAN(NECKE(1,NBACK)) + 
     .                              YTRIAN(NECKE(2,NBACK)) +   
     .                              YTRIAN(NECKE(3,NBACK)) )  
            ALLOCATE(PTRI)
            PTRI%NOTRIAN = NBACK
            PTRI%NEXT_TRIAN => NODELIST(CUR%NODES(1))%PTREL
            NODELIST(CUR%NODES(1))%PTREL => PTRI
            ALLOCATE(PTRI)
            PTRI%NOTRIAN = NBACK
            PTRI%NEXT_TRIAN => NODELIST(CUR%NODES(2))%PTREL
            NODELIST(CUR%NODES(2))%PTREL => PTRI
            ALLOCATE(PTRI)
            PTRI%NOTRIAN = NBACK
            PTRI%NEXT_TRIAN => NODELIST(CUR%NODES(3))%PTREL
            NODELIST(CUR%NODES(3))%PTREL => PTRI
          END IF
          DEALLOCATE (CUR)
!pb          NBACK = NBACK - 1
        END DO

        NTRII = NBACK

!PB        IF (NBACK > 0) THEN
!PB          WRITE (6,*) ' NBACK > 0 IN RPSCUT '
!PB          WRITE (6,*) ' SOMETHINGS WRONG HERE '
!PB        END IF

        DO I=1, NRKNOT
          PTRI => NODELIST(I)%PTREL
          DO WHILE (ASSOCIATED(PTRI)) 
            PTRI2 => PTRI
            PTRI => PTRI%NEXT_TRIAN
            DEALLOCATE(PTRI2)
          END DO
        END DO
      END IF  ! IFIRST BLOCK FINISHED

!  FIND VALUES ON TRIANGLES
      
      DO ITRI=1,NTRII
        VALUES(ITRI) = AORIG(NUMTET(ITRI))
      END DO

      RETURN

      CONTAINS

      RECURSIVE SUBROUTINE TREE_Traverse(node)

      type (TAVLNode), pointer :: node, retnode

      IF (ASSOCIATED(node)) THEN
         XTRIAN(NODE%IND) = NODE%XCO
         YTRIAN(NODE%IND) = NODE%YCO
         IF (ASSOCIATED(node%left) .OR. ASSOCIATED(node%right)) THEN
            Call TREE_Traverse(node%left)
            Call TREE_Traverse(node%right)
         END IF
      END IF

      END SUBROUTINE TREE_Traverse


      END SUBROUTINE RPSCUT








