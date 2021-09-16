!****************************************
!* Description:
!*   Implementiert einen AVL-Baum
!*   in Fortran
!* nach einer Vorlage von
!* Author:
!*   Daniel Hottinger <hodaniel@iiic.ethz.ch>
!*   Department of Computer Science, ETH Zurich
!*   SS 2001
!* Licence: GNU GPL v2 or later
!* Created: Mon May 14 20:16:37 CEST 2001
!* Last update: Tue May 29 20:32:47 CEST 2001
!* Changes:
!*   2001-05-29 Added Search
!****************************************

      MODULE module_avltree

      USE PRECISION
      USE CCONA

      implicit none

      private

      public :: NewTree, DestroyTree, Insert, Search, Remove

      integer, parameter :: less = -1,
     .                      equal = 0,
     .                      more = 1

      type, public :: TAVLNode
        integer :: balance
        type(TAVLNode), pointer :: left, right
        real(dp) :: xco, yco, zco, dist
        integer :: ind
      end type TAVLNode

      type, public :: TAVLTree
        type(TAVLNode), pointer :: root
      end type TAVLTree


      contains

!****************************************
!* Speicherverwaltung
!****************************************

      FUNCTION NewNode(X, Y, Z, DST, IND) RESULT(NODE)

      REAL(DP), INTENT(IN) :: X, Y, Z, DST
      INTEGER, INTENT(INOUT) :: IND
      TYPE(TAVLNode), POINTER :: NODE

      allocate (node)
      node%balance = 0
      NULLIFY(node%left)
      NULLIFY(node%right)
      node%xco = x
      node%yco = y
      node%zco = z
      node%dist = dst
      node%ind = ind

      END FUNCTION NewNode


      RECURSIVE SUBROUTINE DestroyNode(node)

      type(tavlnode), pointer :: node
      type(tavlnode), pointer :: cur, l, r

      IF (ASSOCIATED(node)) THEN
         if (associated(node%left)) then
           l => node%left
           call DestroyNode (l)
           nullify(node%left)
         end if
         if (associated(node%right)) then
           r => node%right
           call DestroyNode (r)
           nullify(node%right)
         end if
         cur => node
         if (associated(cur)) DEALLOCATE(cur)
         NULLIFY(node)
      END IF

      END SUBROUTINE DestroyNode


      FUNCTION NewTree () result(tree)

      type (tavltree), pointer :: tree

      ALLOCATE (tree)
      NULLIFY (tree%root)

      END FUNCTION NewTree


      SUBROUTINE DestroyTree (tree)

      type (tavltree), pointer :: tree

      IF (ASSOCIATED(tree)) THEN
         call DestroyNode(tree%root)
         DEALLOCATE(tree)
      END IF
      END SUBROUTINE DestroyTree


!****************************************
!* Rotationen
!****************************************

!****************************************
!* n(l,r(rl,rr)) =>
!* r(n(l,rl),rr)
!*
!*   n                r
!*  / `_           _.' \
!* l    r    =>   n     rr
!*     / \       / \
!*    rl  rr    l   rl
!****************************************

      FUNCTION RotLeft(node) result(right)

      type(TAVLNode), POINTER :: node
      type(TAVLNode), POINTER :: left, right
      integer :: abal, bbal

      left => node%left
      right => node%right

      node%right => right%left
      right%left => node

      abal = node%balance
      bbal = right%balance

      IF (bbal <= EQUAL) THEN
         IF (abal >= MORE) THEN
            right%balance = bbal - 1
         ELSE
            right%balance = abal + bbal - 2
         END IF
         node%balance = abal - 1
      ELSE
         IF (abal <= bbal) THEN
            right%balance = abal - 2
         ELSE
            right%balance = bbal - 1
         END IF
         node%balance = abal - bbal - 1
      END IF

      END FUNCTION RotLeft


!****************************************
!* n(l(ll,lr),r) =>
!* l(ll,n(lr,r))
!*
!*       n        l
!*    _.' \      / `_
!*   l     r => ll   n
!*  / \             / \
!* ll  lr          lr  r
!****************************************

      FUNCTION RotRight(node) RESULT(left)

      type(TAVLNode), pointer :: node
      type(TAVLNode), pointer :: left, right
      integer :: abal, bbal

      left => node%left
      right => node%right

      node%left => left%right
      left%right => node

      abal = node%balance
      bbal = left%balance

      IF (bbal <= EQUAL) THEN
         IF (abal < bbal) THEN
            left%balance = bbal + 1
         ELSE
            left%balance = abal + 2
         END IF
         node%balance = abal - bbal + 1
      ELSE
         IF (abal <= LESS) THEN
            left%balance = bbal + 1
         ELSE
            left%balance = abal + bbal + 2
         END IF
         node%balance = abal + 1
      END IF
      END FUNCTION RotRight


!****************************************
!* Ausbalancieren
!****************************************

!****************************************
!* Beispiel: (Doppelrotation)
!* 2(1,6(5(4),7)) =>
!* 2(1,5(4,6(,7))) =>
!* 5(2(1,4),6(,7))
!*
!*   2             2
!*  / `-._        / `_            5
!* 1      6      1    5          / \
!*       / \  =>     / \    =>  2   6
!*      5   7       4   6       ^    \
!*     /                 \     1 4    7
!*    4                   7
!****************************************

      FUNCTION BalanceNode(node) result(retnode)

      type(TAVLNode), pointer :: node, retnode

      IF (node%balance < LESS) THEN
         IF (node%left%balance > EQUAL) THEN
            node%left => RotLeft(node%left)
         END IF
         retnode => RotRight(node)
      ELSEIF (node%balance > MORE) THEN
         IF (node%right%balance < EQUAL) THEN
            node%right => RotRight(node%right) ! siehe Beispiel
         END IF
         retnode => RotLeft(node)
      END IF
      END FUNCTION BalanceNode


! Called after a node of node^.left was removed
      FUNCTION RestoreLeftBalance(node, oldbalance) result(retnode)

      type(TAVLNode), pointer :: node
      type(TAVLNode), pointer :: retnode
      integer :: oldbalance

      IF (.NOT.ASSOCIATED(node%left)) THEN
         node%balance = node%balance + 1
      ELSEIF ((node%left%balance /= oldbalance) .AND.
     .        (node%left%balance == 0)) THEN
    ! left tree shrunk
         node%balance = node%balance + 1
      END IF
      IF (node%balance > MORE) THEN
         retnode => BalanceNode(node)
      ELSE
         retnode => node
      END IF

      END FUNCTION RestoreLeftBalance


! Called after a node of node^.right was removed
      FUNCTION RestoreRightBalance(node, oldbalance) result(retnode)
      type(TAVLNode), pointer :: node
      integer :: oldbalance
      type(TAVLNode), pointer :: retnode

      IF (.NOT.ASSOCIATED(node%right)) THEN
         node%balance = node%balance - 1
      ELSEIF ((node%right%balance /= oldbalance) .AND.
     .        (node%right%balance == 0)) THEN
         ! right tree shrunk
         node%balance = node%balance - 1
      END IF
      IF (node%balance > LESS) THEN
         retnode => BalanceNode(node)
      ELSE
         retnode => node
      END IF
      END FUNCTION RestoreRightBalance


      RECURSIVE FUNCTION RemoveNodeMostLeft(node, leftmost)
     .          result(retnode)

      type(TAVLNode), pointer :: node, leftmost
      integer :: oldbalance
      type(TAVLNode), pointer :: retnode

      IF (.NOT.ASSOCIATED(node%left)) THEN
         leftmost => node
         retnode => node%right
         RETURN
      END IF

      oldbalance = node%left%balance
      node%left => RemoveNodeMostLeft(node%left, leftmost)
      retnode => RestoreLeftBalance(node, oldbalance)

      END FUNCTION RemoveNodeMostLeft


!****************************************
!* grundlegende Operationen
!****************************************


      RECURSIVE FUNCTION InsertNode(node, x, y, z, dst, ind, inserted)
     .          result(retnode)

      type (TAVLNode), pointer :: node, retnode
      real(dp), intent(in) :: x, y, z, dst
      integer, intent(inout) :: ind
      logical, intent(INOUT) :: INSERTED

      integer :: relation, oldbalance

      IF (.NOT.ASSOCIATED(node)) THEN
         inserted = .TRUE.
         retnode => NewNode(x, y, z, dst, ind)
         RETURN
      END IF

      relation = cmp(x, y, z, node%xco, node%yco, node%zco, node%dist)
      IF (relation == EQUAL) THEN
         ! Don't insert dublicate key/value
         inserted = .FALSE.
         ind = node%ind
         retnode => node
         RETURN
      ELSEIF (relation == LESS) THEN
         IF (ASSOCIATED(node%left)) THEN
            oldbalance = node%left%balance
            node%left => InsertNode(node%left,x,y,z,dst,ind,inserted)
            IF ((oldbalance /= node%left%balance) .AND.
     .          (node%left%balance /= 0)) THEN
            ! Tree has grown
            node%balance = node%balance - 1
            END IF
         ELSE
            inserted = .TRUE.
            node%left => NewNode(x, y, z, dst, ind)
            node%balance = node%balance - 1
         END IF
      ELSEIF (relation == MORE) THEN
         IF (ASSOCIATED(node%right)) THEN
            oldbalance = node%right%balance
            node%right => InsertNode(node%right,x,y,z,dst,ind,inserted)
            IF ((oldbalance /= node%right%balance) .AND.
     .          (node%right%balance /= 0)) THEN
               ! Tree has grown
               node%balance = node%balance + 1
            END IF
         ELSE
            inserted = .TRUE.
            node%right => NewNode(x, y, z, dst, ind)
            node%balance = node%balance + 1
         END IF
      ENDIF

      IF (inserted) THEN
         IF (ABS(node%balance) > 1) THEN
            node => BalanceNode(node)
         END IF
      END IF

      retnode => node

      END Function InsertNode


      SUBROUTINE Insert(tree, x, y, z, dst, ind, inserted)

      type(TAVLTree), pointer :: tree
      real(dp) :: x, y, z, dst
      integer :: ind
      logical :: inserted

      IF (ASSOCIATED(tree)) THEN
         inserted = .FALSE.
         tree%root => InsertNode(tree%root, x, y, z, dst, ind, inserted)
!         call traverse (tree%root,0)
      END IF
      END SUBROUTINE Insert


!pb      RECURSIVE INTEGER Function SearchNode(node, x, y, z)
      RECURSIVE INTEGER Function SearchNode(node, x, y, z) RESULT (nr)

      type (TAVLNode), pointer :: node
      real(dp), intent(in) :: x, y, z

      integer :: relation, ind

      IF (.NOT.ASSOCIATED(node)) THEN
!pb         SearchNode = 0
         nr = 0
         RETURN
      END IF

      relation = cmp(x, y, z, node%xco, node%yco, node%zco, node%dist)
      IF (relation == EQUAL) THEN
         ind = node%ind
      ELSEIF (relation == LESS) THEN
         ind = SearchNode(node%left, x, y, z)
      ELSE
         ind = SearchNode(node%right, x, y, z)
      END IF

!pb      SearchNode = ind
      nr = ind

      END Function SearchNode


      Integer Function Search(tree, x, y, z)

      type (TAVLTree), pointer :: tree
      real(dp), intent(in) :: x, y, z

      Search = 0
      IF (ASSOCIATED(tree)) THEN
         Search = SearchNode(tree%root, x, y, z)
      END IF
      END FUNCTION Search

!****************************************
!* Beispiel:
!* n(l(:,:),r(rl(rll,rlr),rr)) =>
!* n(l(:,:),r(rl(,rlr),rr)) =>
!* rll(l(:,:),r(rl(,rlr),rr))
!*    n                    n                 rll
!*   / `--...___          / `-..__          /   `-..__
!*  l           r        l        r        l          r
!*  ^      __.-' \   =>  ^   __.-' \   =>  ^     __.-' \
!* : :    rl      rr    : : rl      rr    : :   rl      rr
!*      _'  \                 \                   \
!*     rll   rlr               rlr                 rlr
!****************************************

      RECURSIVE Function RemoveNode(node, x, y, z) result(retnode)

      type (TAVLNode), pointer :: node, retnode
      real(dp), intent(in) :: x, y, z

      integer :: relation, oldbalance
      type (TAVLNode), pointer :: garbage, newroot

      IF (.NOT.ASSOCIATED(node)) THEN
         NULLIFY(retnode)
         RETURN
      END IF

      relation = cmp(x, y, z, node%xco, node%yco, node%zco, node%dist)
      IF (relation == EQUAL) THEN
         garbage => node
         IF (.NOT.ASSOCIATED(node%right)) THEN
            retnode => node%left
         ELSE
            oldbalance = node%right%balance
            ! new right node is the leftmost of the right tree
            ! Beispiel
            node%right => RemoveNodeMostLeft(node%right, newroot)
            newroot%left => node%left
            newroot%right => node%right
            newroot%balance = node%balance
            retnode => RestoreRightBalance(newroot, oldbalance)
         END IF
         ! free *only* the removed node
         Nullify(garbage%right)
         Nullify(garbage%left)
         CALL DestroyNode(garbage)
      ELSEIF (relation == LESS) THEN
         IF (ASSOCIATED(node%left)) THEN
            oldbalance = node%left%balance
            node%left => RemoveNode(node%left, x, y, z)
            retnode => RestoreLeftBalance(node, oldbalance)
         END IF
      ELSEIF (relation == MORE) THEN
         IF (ASSOCIATED(node%right)) THEN
            oldbalance = node%right%balance
            node%right => RemoveNode(node%right, x, y, z)
            retnode => RestoreRightBalance(node, oldbalance)
         END IF
      END IF

      END FUNCTION RemoveNode


      Subroutine Remove(tree, x, y, z)

      type (TAVLTree), pointer :: tree
      real(dp), intent(in) :: x, y, z

      IF (ASSOCIATED(tree)) THEN
         tree%root => RemoveNode(tree%root, x, y, z)
      END IF
      END SUBROUTINE Remove


!****************************************
!* Debug
!****************************************

      RECURSIVE SUBROUTINE Traverse(node,iunout)

      type (TAVLNode), pointer :: node, retnode
      integer, intent(in) :: iunout

      IF (ASSOCIATED(node)) THEN
         if (abs(node%balance) > 2) then
            write (iunout,*) node%balance,node%ind
            write (iunout,*) node%xco, node%yco, node%zco
         end if
!         write (iunout,*) node%balance,node%ind, node%xco, node%yco, 
!     .                    node%zco
         IF (ASSOCIATED(node%left) .OR. ASSOCIATED(node%right)) THEN
!            Write (iunout,*) "("
            Call Traverse(node%left,iunout)
!            IF (ASSOCIATED(node%right))  THEN
!               Write (iunout,*) ","
!            END IF
            Call Traverse(node%right,iunout)
!            Write (iunout,*) ")"
         END IF
      END IF

      END SUBROUTINE Traverse


      Subroutine Dump(tree)
      type(TAVLTree), pointer :: tree

      call Traverse(tree%root,0)
      END SUBROUTINE Dump


      integer function cmp (x1, y1, z1, x2, y2, z2, dist)

      real(dp), intent(in) :: x1, y1, z1, x2, y2, z2, dist
      real(dp) :: dx, dy, dz

!      IF ( (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2 < 1.D-10) THEN
!      IF ( SQRT((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)/dist < EPS5) THEN
      IF ( ((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)/dist**2 < EPS10) THEN
         cmp = equal
      else
!         if (abs(x1-x2) > 1.E-5) then
         if (abs(x1-x2)/max(x1,x2,eps10) > eps5) then
            cmp=merge(less,more,x1 < x2)
!         else if (abs(y1-y2) > 1.e-5) then
         else if (abs(y1-y2)/max(y1,y2,eps10) > eps5) then
            cmp=merge(less,more,y1 < y2)
         else
            cmp=merge(less,more,z1 < z2)
         end if
      end if

      return
      end function cmp


      integer function cmp_old (x1, y1, z1, x2, y2, z2)

      real(dp), intent(in) :: x1, y1, z1, x2, y2, z2
      real(dp) :: dx, dy, dz

      IF ( (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2 < 1.D-10) THEN
         cmp_old = equal
      else
         dx=x1-x2
         if (abs(dx) < 1.E-5) then
!  x1 = x2
            dy=y1-y2
            if (abs(dy) < 1.E-5) then
!  x1 = x2 & y1 = y2
               dz=z1-z2
               if (abs(dz) < 1.E-5) then
!  x1 = x2 & y1 = y2 & z1 = z2
                  cmp_old = equal
               else if (dz < 0.d0) then
!  x1 = x2 & y1 = y2 & z1 < z2
                  cmp_old = less
               else
!  x1 = x2 & y1 = y2 & z1 > z2
                  cmp_old = more
               end if
            else if (dy < 0.d0) then
!  x1 = x2 & y1 < y2
               cmp_old = less
            else
!  x1 = x2 & y1 > y2
               cmp_old = more
            end if
         else if (dx < 0.d0) then
!  x1 < x2
            cmp_old = less
         else
!  x1 > x2
            cmp_old = more
         end if
      end if

      end function cmp_old


      END MODULE module_avltree


