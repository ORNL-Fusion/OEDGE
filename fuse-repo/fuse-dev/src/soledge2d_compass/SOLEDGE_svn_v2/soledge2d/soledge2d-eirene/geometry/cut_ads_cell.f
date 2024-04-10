      subroutine EIRENE_cut_ads_cell
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_cgrid
      use EIRMOD_cadgeo
      use EIRMOD_clgin
      use EIRMOD_ctrig
      use EIRMOD_ccona
      use EIRMOD_comusr
      use EIRMOD_comprt
 
      implicit none
      integer :: EIRENE_learc1, ilim, ic1, ic2, is1, is2, izell, izello,
     .           is, j, no, igo, ihit1, ihit2, i
      real(dp) :: xx, yy, richtx, richty, ax, ay, v, vx, vy, t
      logical :: lc1(3), lc2(3), lhit1, lhit2
      TYPE(TRI_ELEM), POINTER :: CUR1, cur2
 
      if (levgeo /= 4) return
      IF (NOPTIM < NR1ST) return
 
 
! by default switch off view to additional surfaces
      IF (NLIMPB >= NLIMPS) THEN
         IGJUM3 = 1
      ELSE
        NO=NOT(0)
! NUMBER OF INTEGERS USED TO STORE SURFACE INFORMATION
        IGO=NLIMI/NBITS
        IF (MOD(NLIMI,NBITS) > 0) IGO = IGO + 1
        IGJUM3(1:NR1ST,1:IGO) = NO
      END IF
 
      do ilim = 1, nlimi
 
! surface is explicitly taken out --> seen from nowhere
        if (igjum0(ilim) == 1) cycle
 
        if ((rlb(ilim) < 2.) .or. (rlb(ilim) >= 3)) then
! algebraic equation, probably curved surface -->
! to be seen from everywhere
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM3(1:NR1ST,ILIM) = 0
          ELSE
            do izell = 1, nr1st
              CALL EIRENE_BITSET (IGJUM3,0,NOPTIM,IZELL,ILIM,0,NBITS)
            end do
          END IF
 
        else
 
! find cell index of start point of additional surface
          ic1 = EIRENE_LEARC1(P1(1,ilim),P1(2,ilim),P1(3,ilim),IS1,1,
     .                 NR1STM,.FALSE.,.FALSE.,0,'CUT_ADS_CELL')
 
! find cell index of end point of additional surface
          ic2 = EIRENE_LEARC1(P2(1,ilim),P2(2,ilim),P2(3,ilim),IS2,1,
     .               NR1STM,.FALSE.,.FALSE.,0,'CUT_ADS_CELL')
 
! check if additional surface starts or end at a cornerpoint of a triangle
          ihit1 = 0
          ihit2 = 0
 
          do i = 1, 3
            lc1(i) = sqrt((xtrian(necke(i,ic1))-p1(1,ilim))**2 +
     .                    (ytrian(necke(i,ic1))-p1(2,ilim))**2) < eps5
            lc2(i) = sqrt((xtrian(necke(i,ic2))-p2(1,ilim))**2 +
     .                    (ytrian(necke(i,ic2))-p2(2,ilim))**2) < eps5
            if (lc1(i)) ihit1 = necke(i,ic1)
            if (lc2(i)) ihit2 = necke(i,ic2)
          end do
 
          lhit1 = ihit1 > 0
          lhit2 = ihit2 > 0
 
          write (0,*) ' ilim, ic1, ic2 ', ilim, ic1, ic2
 
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM3(IC1,ILIM) = 0
            IGJUM3(IC2,ILIM) = 0
          ELSE
            CALL EIRENE_BITSET (IGJUM3,0,NOPTIM,IC1,ILIM,0,NBITS)
            CALL EIRENE_BITSET (IGJUM3,0,NOPTIM,IC2,ILIM,0,NBITS)
          END IF
 
! find cells between ic1 and ic2 along additional surface
          is = 0
          xx = P1(1,ilim)
          yy = P1(2,ilim)
          vx = p2(1,ilim) - p1(1,ilim)
          vy = p2(2,ilim) - p1(2,ilim)
 
          if (lhit1) then
! startpoint on cornerpoint
            cur1 => coortri(ihit1)%ptri
            izell = cur1%notri
          else
! startpoint inside triangle
            izell = ic1
          end if
 
          lalong: do
 
            write (0,*) ' izell ',izell
 
! last cell reached
            if (lhit2) then
              cur2 => coortri(ihit2)%ptri
              do
                if (cur2%notri == izell) exit lalong
                if (.not.associated(cur2%next_tri)) exit
                cur2 => cur2%next_tri
              end do
            else
              if (izell == ic2) exit lalong
            end if
 
            if (izell == 0) exit lalong
 
            do j = 1, 3
 
! don't check entrance side of triangle
              if (j == is) cycle
 
              RICHTX = VTRIX(J,IZELL)
              RICHTY = VTRIY(J,IZELL)
              AX = XTRIAN(NECKE(J,IZELL))-XX
              AY = YTRIAN(NECKE(J,IZELL))-YY
              V  = (AX*VY-AY*VX)/(RICHTY*VX-RICHTX*VY+EPS60)
 
              IF (V.GE.0._DP.AND.V.LE.1._DP) THEN
                IF (ABS(VX).GT.ABS(VY)) THEN
                  T=(AX+V*RICHTX)/VX
                ELSE
                  T=(AY+V*RICHTY)/VY
                ENDIF
                IF (T .LE. 0.) CYCLE
C     ZELLENNUMMER DES NEUEN DREIECKS
                IZELLO=IZELL
                IZELL = NCHBAR(J,IZELLO)
C     SEITENNUMMER DES NEUEN DREIECKS
                IS = NSEITE(J,IZELLO)
                IF (NLIMPB >= NLIMPS) THEN
                  IGJUM3(IZELL,ILIM) = 0
                ELSE
                  CALL EIRENE_BITSET
     .  (IGJUM3,0,NOPTIM,IZELL,ILIM,0,NBITS)
                END IF
                cycle lalong
              end if
 
            end do ! j
 
! wrong triangle check, try next
            if (lhit1) then
              if (.not.associated(cur1%next_tri)) then
                write (iunout,*) ' error in cut_ads_cell '
                write (iunout,*) ' no triangles found along surface ',
     .                           ilim
                exit lalong
              else
                cur1 => cur1%next_tri
                izell = cur1%notri
              end if
            end if
          end do lalong ! loop along additional surface
        end if
      end do ! ilim
 
      return
      end subroutine EIRENE_cut_ads_cell
