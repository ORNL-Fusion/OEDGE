c  modusr routine for EIRENE-OSM-interface (triangles), for iteration
C                 between radiation and non-LTE gas fields.
C  1st attempt: 2 state atoms, one single line, Lyman alpha only,
C               use rates from h2vibr.
c
c  2nd attempt: multi-lines, use col-rad routine here in each cell.
c               this gives identical results to 1st attempt in case of
c               one single line only, but is still more transparent coding

c  a general complication may arise in this iteration scheme, if
c  background tallies are given on a finer grid ("ir-grid") than the
c  test particle tallies ("irt-grid"). One "irt-cell" may comprise
c  several "ir cells". IRT=NCLTAL(IR), if "ir cell" IR belongs to
c  to the coarser "irt cell" IRT.

      subroutine modusr
      use precision
      use parmmod
      use comusr
      use ccona
      use cestim
      use ctrcei
      use comxs
      use cgeom
      use cinit
      use ctrig
      use csdvi
      use clogau
      use comsou
      implicit none
      integer :: IPLS, IR, I, J, IRT, ITRI, MTRI
      REAL(DP) :: DDI4_2, DDI4_3, DDI4_4, DDI4_5, DDI4_6,
     .            X, Y1, Y2, Y3, A, B, C
      REAL(DP) :: RESD_1, RAUSCHD_1, X1, X2, X3, X1Q,
     .            X2Q, X3Q, TOTD_1
      REAL(DP) :: VOLOUT(NRAD), DPLUS(NRAD), D_1(NRAD), TDNE(NRAD),
     .            TETRI(NRAD), DETRI(NRAD),
     .            TABDS_el_imp(NRAD), PDA(NRAD), PDPH(NRAD), 
     .            NETEMIS(NRAD), TABDS_rad(NRAD), 
     .            TDNE2(NRAD), CLST(NRAD)
      real(dp) :: dold_ir,dnew_irt, dnew_ir,
     .            pphpht2_ir,pphpht2_irt,
     .            pphpht3_ir,pphpht3_irt,
     .            pphpht4_ir,pphpht4_irt,
     .            pphpht5_ir,pphpht5_irt,
     .            pphpht6_ir,pphpht6_irt,
     .            sigma_ir,sigma_irt
      real(dp) :: pop0(40), pop1(40), pop2(40), q2(40)
      real(dp) :: alpcr, scr, scrrad
      real(dp) :: e_alpcr, e_scr, e_scrrad
      real(dp) :: e_alpcr_t, e_scr_t, e_scrrad_t
      INTEGER :: NOTRI(NRAD), ICELLST(NRAD), IRINGST(NRAD), 
     .           ICELLRD(NRAD), IRINGRD(NRAD)
c slmod begin
c...  Steve's EIRENE results output routine:
      CALL OUTUSR
      IF (IITER.EQ.NITER) RETURN
c slmod end
c
c  this routine requires multiple ion temperatures.
c  check, if these are indeed available
      if (nplsti.le.1) then
        write (6,*) 'multiple ion temperatures are needed in modusr'
        write (6,*) 'modify input flag indpro(2) (block 5)         '
        write (6,*) 'exit called'
        call exit_own(1)
      endif
c
c  prepare some tallies for output for OSM code, on OSM grid
      NOTRI=0
      VOLOUT=0._DP
      DPLUS=0._DP
      D_1=0._DP
      TDNE=0._DP
      TDNE2=0._DP
      TETRI=0._DP
      DETRI=0._DP
      PDA = 0._DP
      PDPH = 0._DP
      NETEMIS = 0._DP
      CLST = 0._DP
      ICELLST = 0
      IRINGST = 0

      call prousr (clst,1+6*npls,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,
     .             0._dp,ntrii)
      ICELLRD=CLST
      call prousr (clst,2+6*npls,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,
     .             0._dp,ntrii)
      IRINGRD=CLST

      DO IR=1,NTRII
         IF (IYTRI(IR) == -1) THEN
            ITRI = IXTRI(IR)
            NOTRI(ITRI) = NOTRI(ITRI) + 1
            ICELLST(ITRI) = ICELLRD(IR)
            IRINGST(ITRI) = IRINGRD(IR)
! SUMME(VOL)
            VOLOUT(ITRI) = VOLOUT(ITRI) + VOL(IR)
! SUMME(D+ * VOL)
            DPLUS(ITRI) = DPLUS(ITRI) + DIIN(1,IR)*VOL(IR)
! SUMME(D_1 * VOL)
            D_1(ITRI) = D_1(ITRI) + DIIN(2,IR)*VOL(IR)
! SUMME(NE * < >)
            TDNE(ITRI) = TDNE(ITRI) + TABDS1(1,IR)
            TDNE2(ITRI) = TDNE2(ITRI) + TABDS1(2,IR)
! SUMME(TE * NE * VOL)
            TETRI(ITRI) = TETRI(ITRI) + TEIN(IR)*DEIN(IR)*VOL(IR)
! SUMME(NE * VOL)
            DETRI(ITRI) = DETRI(ITRI) + DEIN(IR)*VOL(IR)
! SUMME(PDENA * VOL)
            PDA(ITRI) = PDA(ITRI) + PDENA(1,IR)*VOL(IR)
! SUMME(PDENPH * VOL)
            PDPH(ITRI) = PDPH(ITRI) + PDENPH(2,IR)*VOL(IR)
! SUMME((EPPHT+EPHPHT) * VOL) (ephpht<0)
            NETEMIS(ITRI) = NETEMIS(ITRI) + 
     .                      (EPPHT(IR)+EPHPHT(IR))*VOL(IR)
         END IF
      END DO

      MTRI=COUNT(NOTRI>0)
      DO ITRI=1,MTRI
! SUMME(D+ * VOL) / SUMME(VOL)
         DPLUS(ITRI) = DPLUS(ITRI) / (VOLOUT(ITRI)+EPS60)
! SUMME(D_1 * VOL) / SUMME(VOL)
         D_1(ITRI) = D_1(ITRI) / (VOLOUT(ITRI)+EPS60)
! SUMME(TE * NE * VOL) / SUMME(NE * VOL)
         TETRI(ITRI) = TETRI(ITRI) / (DETRI(ITRI)+EPS60)
! SUMME(NE * VOL) / SUMME(VOL)
         DETRI(ITRI) = DETRI(ITRI) / (VOLOUT(ITRI)+EPS60)
! SUMME(NE * < >) / NE = SUMME(NE*< >) / (SUMME(NE*VOL)/SUMME(VOL))
         TDNE(ITRI) = TDNE(ITRI) / (DETRI(ITRI)+EPS60)
         TDNE2(ITRI) = TDNE2(ITRI) / (DETRI(ITRI)+EPS60)
! SUMME(PDENA * VOL) / SUMME(VOL)
         PDA(ITRI) = PDA(ITRI) / (VOLOUT(ITRI)+EPS60)
! SUMME(PDENPH * VOL) / SUMME(VOL)
         PDPH(ITRI) = PDPH(ITRI) / (VOLOUT(ITRI)+EPS60)
! SUMME((EPPHT+EPHPHT) * VOL) / SUMME(VOL) (ephpht<0)
         NETEMIS(ITRI) = NETEMIS(ITRI) / (VOLOUT(ITRI)+EPS60)
      END DO

c  write output file for OSM code, from last iteration

      WRITE (96,'(I6)') MTRI
      WRITE (96,'(A6,9A12)') 'NO.','VOL ',
     .     'TABDS1','TABDS2', 'CELL INDEX','RING INDEX' 
      DO IR=1,MTRI
         WRITE (96,'(I6,3ES12.4,2I12)') IR, VOLOUT(IR),
     .         TDNE(IR), TDNE2(IR), ICELLST(IR), IRINGST(IR)
      END DO

c  output for OSM from last iteration: done

c  start to work on next iteration

      RESD_1 = 0._DP
      RAUSCHD_1 = 0._DP
      TOTD_1 = 0._DP


C  ITERATE ON BACKGROUND SPECIES
c  be careful: background grid may differ from test particle grid
c              ir, irt.  "ir-grid" may be finer than "irt-grid"
c
      DO IR=1,NRAD
         IRT=NCLTAL(IR)

!  D_1 (ground state atoms. Needed for photon absorption rate)
         DOLD_ir=DIIN(2,IR)
         DNEW_irt=PDENA(1,IRT)
         TIIN(2,IR)=EDENA(1,IRT)/(DNEW_irt+EPS60)/1.5
         DNEW_ir=DNEW_irt
         DIIN(2,IR)=DNEW_ir
c  total ground state particle number
         TOTD_1 = TOTD_1 + DIIN(2,IR)*VOL(IR)
c  residual: rate of change in Lyman-alpha absorption, 1/sec
c            assume: absorption rate per cm**3 on "ir-grid"
c                    is the same as the absorption rate in "irt-grid".
         PPHPHT2_irt=PPHPHT(2,irt)
         PPHPHT2_ir=PPHPHT2_irt
         PPHPHT3_irt=PPHPHT(3,irt)
         PPHPHT3_ir=PPHPHT3_irt
         PPHPHT4_irt=PPHPHT(4,irt)
         PPHPHT4_ir=PPHPHT4_irt
         PPHPHT5_irt=PPHPHT(5,irt)
         PPHPHT5_ir=PPHPHT5_irt
         PPHPHT6_irt=PPHPHT(6,irt)
         PPHPHT6_ir=PPHPHT6_irt

         SIGMA_irt=SIGMA(3,irt)
         SIGMA_ir=SIGMA_irt ! noise is in %, so this is questionable
         RESD_1 = RESD_1 + VOL(IR) *
     .            ABS(DOLD_ir-DNEW_ir)/(DOLD_ir+EPS60)*PPHPHT2_ir/ELCHA
c  statistical noise level in Lyman-alpha absorption rate, 1/sec
         RAUSCHD_1=RAUSCHD_1-SIGMA_ir/100._DP*PPHPHT2_ir*
     .             VOL(IR)/ELCHA
c
         TIIN(3,IR)=TIIN(2,IR)
         TIIN(5,IR)=TIIN(2,IR)
         TIIN(7,IR)=TIIN(2,IR)
         TIIN(9,IR)=TIIN(2,IR)
         TIIN(11,IR)=TIIN(2,IR)

c   save some lines on background tallies, for later plotting and printing
!  Lyman_a
         TIIN(13,IR)=EDENPH(2,IRT)/(PDENPH(2,IRT)+EPS60)/1.5
         DIIN(13,IR)=PDENPH(2,IRT)
!  Lyman_b
         TIIN(14,IR)=EDENPH(3,IRT)/(PDENPH(3,IRT)+EPS60)/1.5
         DIIN(14,IR)=PDENPH(3,IRT)



!  NEXT: NEW EFFECTIVE IONIZATION RATE
c  needed for atom loss due to rad. trapping, from the atoms point of view.

c  this is now done by a call to the full collisional radiation model
c  for atoms. Strictly, only the single right hand side due to rad. trapping
c  needs to be done, because the other contributions (ionising, recomb.,
c  dissoc ex and dissoc.rc.) need not be recalculated.
c
c  alternative (faster): store the coll rad. solutions pop2 and scrrad
c                        for each q2-unit basis vector once, 
c                        and then only form the linear superposition
c                        at each iteration
c
c  done on the finer "ir-grid". Assume here again: pphpht_ir=pphpht_irt
         IF (DIIN(2,IR) > 0._DP) THEN
            DDI4_2=-PPHPHT2_ir/(DOLD_ir+EPS60)/ELCHA
            DDI4_3=-PPHPHT3_ir/(DOLD_ir+EPS60)/ELCHA
            DDI4_4=-PPHPHT4_ir/(DOLD_ir+EPS60)/ELCHA
            DDI4_5=-PPHPHT5_ir/(DOLD_ir+EPS60)/ELCHA
            DDI4_6=-PPHPHT6_ir/(DOLD_ir+EPS60)/ELCHA
         ELSE
            DDI4_2=0._DP
            DDI4_3=0._DP
            DDI4_4=0._DP
            DDI4_5=0._DP
            DDI4_6=0._DP
         END IF

         IF (LGVAC(IR,NPLS+1)) THEN
            TABDS1(1,IR)=0._DP
            TABDS1(2,IR)=0._DP
            TABDS_el_imp(IR)=0._DP
            TABDS_rad(IR)=0._DP
         ELSE
c  q2:  rate for population of n from ground state, through photoexciation
            q2 = 0._dp
            q2(2) = ddi4_2
            q2(3) = ddi4_3
            q2(4) = ddi4_4
            q2(5) = ddi4_5
            q2(6) = ddi4_6
c  pop0: coupling to continuum
c  pop1: coupling to ground state, electron impact
c  pop2: coupling to ground state, photon trapping

            CALL H_COLRAD(TEIN(IR),DEIN(IR),Q2,POP0,POP1,POP2,
     .                      ALPCR,    SCR,    SCRRAD,
     .                    E_ALPCR,  E_SCR,  E_SCRRAD,
     .                    E_ALPCR_T,E_SCR_T,E_SCRRAD_T)

c  provide new ionisation rates for next iteration

c  TABDS_total = TABDS_el_imp + TABDS_rad
c
c  first:
c  conventional opt. thin col rad ionisation rate first
c  nothing to be done here, this is already on tabds1(1,...)
c           TABDS1(1,IR) = TABDS_el_imp(IR)
            TABDS1(1,IR) = SCR*DEIN(IR)

c  next: radiation trapping contribution to effective ionisation rate
c  then: introduce TABDS1(2,IR) = TABDS_rad
            TABDS_rad(IR) = SCRRAD
            TABDS1(2,IR)  = TABDS_rad(IR)
c
c  upper state population, 
c  add contributions from electron impact excitation and from rad. trapping
c  both are coupled to ground state

!  D_n_g=D_n_g^elimp+D_n_g^rad, 
!  D_n_g^elimp: coupling to electron impact excitation from groundstate
!  D_n_g^rad:   coupling to radiation absorption from groundstate
c  needed for emission of radiation from the 5 transitions n-->1 considered
c
            DIIN(3,IR) =(POP1(2)+POP2(2))*DIIN(2,IR)
            DIIN(5,IR) =(POP1(3)+POP2(3))*DIIN(2,IR)
            DIIN(7,IR) =(POP1(4)+POP2(4))*DIIN(2,IR)
            DIIN(9,IR) =(POP1(5)+POP2(5))*DIIN(2,IR)
            DIIN(11,IR)=(POP1(6)+POP2(6))*DIIN(2,IR)
         END IF

      END DO  ! IR loop done

!pb avoid overwriting the new calcutaled density of D_1
      CDENMODEL(2)=REPEAT(' ',10)
      CDENMODEL(3)=REPEAT(' ',10)

c  save cpu-time, do not evaluate densities for coupling
c                 to continuum again in each iteration

      CDENMODEL(4)=REPEAT(' ',10)
      CDENMODEL(5)=REPEAT(' ',10)
      CDENMODEL(6)=REPEAT(' ',10)
      CDENMODEL(7)=REPEAT(' ',10)
      CDENMODEL(8)=REPEAT(' ',10)
      CDENMODEL(9)=REPEAT(' ',10)
      CDENMODEL(10)=REPEAT(' ',10)
      CDENMODEL(11)=REPEAT(' ',10)
      CDENMODEL(12)=REPEAT(' ',10)

      CALL PLASMA_DERIV

!pb switch on strata 4, 6, 8, 10 and 12

      NPTS(4) = HUGE(1) - 1
      NPTS(6) = HUGE(1) - 1
      NPTS(8) = HUGE(1) - 1
      NPTS(10) = HUGE(1) - 1
      NPTS(12) = HUGE(1) - 1
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C
      IF ((NFILEL >=1) .AND. (NFILEL <=5)) THEN
         NFILEL=3
         CALL WRPLAM(TRCFLE,0)
      ELSE
         NFILEL=9
         CALL WRPLAM_XDR(TRCFLE,0)
      END IF


      WRITE (6,*) ' D_1 RESIDUUM OF ITERATION     ',RESD_1
      WRITE (6,*) ' D_1 STATISTICAL NOISE         ',RAUSCHD_1
      WRITE (6,*) ' D_1 TOTAL NUMBER OF PARTICLES ',TOTD_1
      WRITE (6,*) ' D_1 REL. RESIDUUM OF ITERATION',RESD_1/TOTD_1
      WRITE (6,*) ' D_1 REL. STATISTICAL NOISE    ',RAUSCHD_1/TOTD_1

      return
      end
