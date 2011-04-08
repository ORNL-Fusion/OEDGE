C     EIRENE04 COMPILATION
C ===== SOURCE: broad_usr.f


      SUBROUTINE BROAD_USR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: diagno.f


c      SUBROUTINE DIAGNO
c      IMPLICIT NONE
c      RETURN
c      END
C ===== SOURCE: geomd.f
C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J

      ndxa=0
      ndya=0

      read (30,*)
      read (30,*)
      read (30,*)
      read (30,*)

1     continue
      read (30,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a110)') zeile

      read (30,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=br(ix,iy,1)
          Y1(IX)=bz(ix,iy,1)
          X2(IX)=br(ix,iy,2)
          Y2(IX)=bz(ix,iy,2)
          X3(IX)=br(ix,iy,4)
          Y3(IX)=bz(ix,iy,4)
          X4(IX)=br(ix,iy,3)
          Y4(IX)=bz(ix,iy,3)
1014    CONTINUE
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IPART=1
      NPOINT(1,IPART)=1
      IY=1
      DO IX=1,NDXA
        DX=BR(IX+1,IY,1)-BR(IX,IY,2)
        DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
        IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
          NPOINT(2,IPART)=IX+1+IPART-1
          IPART=IPART+1
          NPOINT(1,IPART)=IX+1+IPART-1
        ENDIF
      ENDDO
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1
c
C     do j=1,ndya+1
C       write (6,*)
C       write (6,*) 'in geomd polygon ',j
C       write (6,'(1p,6e12.4)') (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD//
      END
C ===== SOURCE: geousr.f
C
C
C
      SUBROUTINE GEOUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE COMPRT
      USE CTRCEI
      USE CCONA
      USE CGEOM
      USE CGRID
      USE CLGIN
      USE CINIT
      USE CPOLYG
      USE CTRIG
      IMPLICIT NONE
      CHARACTER(80) :: ZEILE
      REAL(DP) :: XCOOR, YCOOR, ZCOOR, xan, xen, yan, yen, zan, zen 
      INTEGER :: NADMOD, NASMOD, NORMOD, NRS, IPUNKT, I,NSSIR, NSSIP,
     .           IDIR, IR, IP, IT, IC, IN, NAS
      logical :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4
C
C MODIFY GEOMETRY
C
C
      READ (IUNIN,'(A80)') ZEILE
      READ (IUNIN,'(3I6)') NADMOD,NASMOD,NORMOD

      DO I=1,NADMOD
        READ (IUNIN,'(2I6,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR

        SELECT CASE(IPUNKT)

        CASE DEFAULT
           WRITE (6,*) 'WRONG POINTNUMBER IN ADDUSR '
           WRITE (6,*) 'INPUT LINE READING'
           WRITE (6,'(2I6,1P,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
           WRITE (6,*) ' IS IGNORED '

        CASE (1)
           P1(1,NRS)=XCOOR
           P1(2,NRS)=YCOOR
           P1(3,NRS)=ZCOOR

        CASE (2)
           P2(1,NRS)=XCOOR
           P2(2,NRS)=YCOOR
           P2(3,NRS)=ZCOOR

        CASE (3)
           P3(1,NRS)=XCOOR
           P3(2,NRS)=YCOOR
           P3(3,NRS)=ZCOOR

        CASE (4)
           P4(1,NRS)=XCOOR
           P4(2,NRS)=YCOOR
           P4(3,NRS)=ZCOOR

        CASE (5)
           P5(1,NRS)=XCOOR
           P5(2,NRS)=YCOOR
           P5(3,NRS)=ZCOOR

        CASE (6)
           P6(1,NRS)=XCOOR
           P6(2,NRS)=YCOOR
           P6(3,NRS)=ZCOOR

        END SELECT
      ENDDO

      DO I=1,NASMOD
        READ (IUNIN,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
        IF (IPUNKT.EQ.1) THEN
          P1(1,NAS)=XPOL(NSSIR,NSSIP)
          P1(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSEIF (IPUNKT.EQ.2) THEN
          P2(1,NAS)=XPOL(NSSIR,NSSIP)
          P2(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSE
          WRITE (6,*) 'WRONG POINTNUMBER IN ADDUSR '
          WRITE (6,*) 'INPUT LINE READING'
          WRITE (6,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
          WRITE (6,*) ' IS IGNORED '
        ENDIF
      ENDDO

      DO I = 1, NORMOD
        READ (IUNIN,'(5I6)') IDIR,IR,IP
        IF (IDIR == 1) THEN
          PLNX(IR,IP) = -PLNX(IR,IP)
          PLNY(IR,IP) = -PLNY(IR,IP)
        ELSE IF (IDIR == 2) THEN
          PPLNX(IR,IP) = -PPLNX(IR,IP)
          PPLNY(IR,IP) = -PPLNY(IR,IP)
        ELSE
          WRITE (6,*) ' IDIR =',IDIR,' NOT FORESEEN IN GEOUSR '
          WRITE (6,*) IDIR,IR,IP
          WRITE (6,*) ' IS IGNORED '
        END IF
      END DO
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2 SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C
C
C
C  SET SOME VOLUMES EXPLIZIT
C
C
C  MODIFY REFLECTION MODEL AT TARGET PLATES
C
C     do i=1,nlimps
C       do isp=1,natmi+nmoli+nioni
C         recyct(isp,i)=1.
C       enddo
C     enddo

      if (nlimi >= 50) then
      p1(1,1) = xtrian(necke(2,3254))
      p1(2,1) = ytrian(necke(2,3254))

      p2(1,6) = xtrian(necke(3,3395))
      p2(2,6) = ytrian(necke(3,3395))
      

      p2(1,50) = xtrian(necke(3,3083))
      p2(2,50) = ytrian(necke(3,3083))
      
      p1(1,24) = xtrian(necke(2,3251))
      p1(2,24) = ytrian(necke(2,3251))

      p2(1,23) = xtrian(necke(3,3856))
      p2(2,23) = ytrian(necke(3,3856))
      
      p1(1,7) = xtrian(necke(2,3884))
      p1(2,7) = ytrian(necke(2,3884))
      end if

      RETURN
      END
C ===== SOURCE: iniusr.f


      SUBROUTINE iniUSR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: leausr.f


      FUNCTION LEAUSR(A,B,C)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A, B, C
      INTEGER :: LEAUSR
      LEAUSR=1
      RETURN
      END
C ===== SOURCE: modusr.f
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
c...  Steve's EIRENE results output routine:ain
      WRITE(0,*) 'HERE IN MODUSR',iiter,niter,nphoti

      IF (NPLSI.EQ.5.OR.NPLSI.EQ.18) THEN
c      IF (NPHOTI.EQ.0) THEN                       ! Okay to call MODBGK on the last iteration, any point? 
c      IF (NPHOTI.EQ.0.AND.NITER.GE.1) THEN                   
c      IF (NPHOTI.EQ.0.AND.NITER.GE.1.AND.IITER.LT.NITER) THEN
        WRITE(0,*) 'CALLING MODBGK'
        CALL MODBGK
c        RETURN
      ENDIF

      IF (IITER.GE.NITER) THEN
c        WRITE(0,*) 'CALLING OUTUSR'
c        CALL OUTUSR
      ENDIF

      IF (NPHOTI.EQ.0) RETURN
      WRITE(0,*) 'WHOA! CALLING PHOTON CODE!'
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

      STOP 'DEBUG: COMPILER PROBLEM' ! SEE ALL IYTRY REFERENCES BELOW

      DO IR=1,NTRII
         IF (.TRUE.) THEN
c         IF (IYTRI(IR) == -1) THEN   ! COMPILER PROBLEM IT IYTRY
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
C ===== SOURCE: mshadj.f
C
C
      SUBROUTINE MSHADJ (X1,Y1,X2,Y2,XPLG,YPLG,NPLP,NPOINT,M1,NDX,IR)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*)
      INTEGER, INTENT(IN) ::M1, NDX, IR
      REAL(DP), INTENT(INOUT) :: XPLG(M1,*),YPLG(M1,*)
      INTEGER, INTENT(INOUT) :: NPOINT(2,*), NPLP
      REAL(DP) :: EPS, D, D1
      INTEGER :: NPRT, NP, I
      LOGICAL :: LDUMCL,LPRT
C  LPRT=.TRUE.  : VALID PART
C  LDUMCL=.TRUE.: DUMMY CELL
C
      EPS=1.E-20
      NPRT=0
      NP=0
      LDUMCL=.FALSE.
      LPRT=.FALSE.
C
      DO 1 I=1,NDX
        D=ABS((X1(I)-X2(I))**2+(Y1(I)-Y2(I))**2)
        IF (D.LE.EPS) THEN
          IF (LPRT) THEN
C  ENDE DER VALID ZELLEN
            NPOINT(2,NPRT)=NP
            LPRT=.FALSE.
            LDUMCL=.TRUE.
          ENDIF
        ELSEIF (D.GT.EPS) THEN
C  STARTPUNKT DIESES POLYGONS = ERSTER VALID PUNKT?
          IF (.NOT.LPRT.AND..NOT.LDUMCL) THEN
            LPRT=.TRUE.
            NPRT=NPRT+1
            NP=NP+1
            XPLG(IR,NP)=X1(I)
            YPLG(IR,NP)=Y1(I)
            NPOINT(1,NPRT)=NP
          ELSEIF (.NOT.LPRT.AND.LDUMCL) THEN
            D1=ABS((X1(I+1)-X2(I+1))**2+(Y1(I+1)-Y2(I+1))**2)
            IF (D1.GT.EPS) THEN
C  ENDE DER DUMMYZELLEN: SUCHE NAECHSTE ZELLE MIT D1.GT.0.
              LDUMCL=.FALSE.
              LPRT=.TRUE.
              NPRT=NPRT+1
              NPOINT(1,NPRT)=NP
            ENDIF
          ENDIF
        ENDIF
        NP=NP+1
        XPLG(IR,NP)=X2(I)
        YPLG(IR,NP)=Y2(I)
1     CONTINUE
      NPOINT(2,NPRT)=NP
      NPLP=NPRT
      RETURN
      END
C ===== SOURCE: outusr.f

      SUBROUTINE OUTUSR
      USE precision
      USE parmmod
      USE comusr
      USE ccona
      USE cestim
      USE ctrcei
      USE comxs
      USE cgeom
      USE cinit
      USE ctrig
      USE csdvi
      USE clogau
      USE comsou
      USE coutau
      USE ctetra  ! TET diff
      IMPLICIT none

      REAL*8    :: FTABRC1,FEELRC1

      INTEGER   :: FP,IR,ITRI,MTRI,IIRC,IRRC,IPLS,I1,IADD,IPLSTI,IPLSV,
     .             INC,IFL,NOTRI(NRAD),ICOUNT,ITALLY,IATM,IMOL,IION,I,J,
     .             IADV,IS,NSUR,JJ,IPHOT,NIPLS,IAEI,IRDS,NSIDE,
     .             ISTRAA,ISTRAE,NEW_ITER, JFEXMN, JFEXMX,MSURFG
      INTEGER   :: PROB1,PROB2,IPANU,SAVE_IPANU(100)
      LOGICAL   :: OUTPUT1,BULK_SOURCES(NPLSI)
      REAL*8    :: DDUM(20),RECPAR,RECMOM,RECENG,RECELE,RECADD,EEADD,
     .             SIGNUM,RH2PH2(0:8,0:8),DEF,TEF,RATIO,DEJ,TEI,  
     .             FPRM(6),RCMIN, RCMAX, CONV
      CHARACTER :: FILNAM*8,H123*4,REAC*9,CRC*3

c...  Dynamicallocationize:
      REAL*8  :: CLST(NRAD)
c      REAL*8  :: CLST(NRAD),ADDV2(7,NRAD)
      INTEGER :: ICELLRD(NRAD), IRINGRD(NRAD)

      REAL*8 :: SUMION

c...  Iteration data:
      INTEGER :: ITERNO 
      INTEGER IT,IG,NGAUGE
      REAL*8 XGAUGE(3),YGAUGE(3),RGAUGE(3),DIST,X1,X2,X3,Y1,Y2,Y3,FACT
      REAL*8, ALLOCATABLE :: XCEN(:),YCEN(:),ZCEN(:),T_VOL(:),
     .                       T_PDENA(:,:,:),T_EDENA(:,:,:),
     .                       T_PDENM(:,:,:),T_EDENM(:,:,:)
c     .                       T_PDENM(:,:,:),T_EDENM(:,:,:)

      INTEGER :: logfp = 6

      DATA ITERNO / 1 /

      SAVE



      ! TETRAHEDRON
      ntrii = ntet  ! TET diff



      OUTPUT1=.TRUE.
      PROB1=10000000
      PROB2=0

      CLST=0.0D0
      call prousr (clst,1+6*npls,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,
     .             0._dp,ntrii)
      ICELLRD=CLST
      CLST=0.0D0
      call prousr (clst,2+6*npls,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,
     .             0._dp,ntrii)
      IRINGRD=CLST

      NOTRI = 0
c      DO IR=1,NTRII
c         IF (IYTRI(IR) == -1) THEN
c            ITRI = IXTRI(IR)
c            NOTRI(ITRI) = NOTRI(ITRI) + 1
cc            ICELLST(ITRI) = ICELLRD(IR)
cc            IRINGST(ITRI) = IRINGRD(IR)
c         ENDIF
c      ENDDO
c      MTRI=COUNT(NOTRI>0)

      WRITE(0,*) 'HERE IN OUTUSR'

      FP = 99
      OPEN(UNIT=FP,FILE='eirene.transfer',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=98)



      WRITE(FP,80) '* NTRII  :',NTRII,NRAD,MTRI
      WRITE(FP,80) '* NSTORDR:',NSTORDR
      WRITE(FP,80) '* NATMI  :',NATMI,NMOLI,NIONI,NPHOTI
      WRITE(FP,80) '* NPLSI  :',NPLSI
      WRITE(FP,80) '* MPLSTI :',MPLSTI
      WRITE(FP,80) '* MPLSV  :',MPLSV

c...  Misc. electron data, VOL, VOLTAL? (are they the same?), NCLTAL(IR)=IR?

c...  Determine the number of surfaces:
      NSIDE = 0 ! TET diff

      NSUR=0
      DO IR=1,NTRII  
        DO IS=1,NSIDE
          IF (INSPAT(IS,IR).NE.0) NSUR=NSUR+1
        ENDDO
      ENDDO
c
c     ----------------------------------------------------------------------
c     Bulk species:
      BULK_SOURCES = .FALSE.

      NIPLS=1
      BULK_SOURCES(1) = .TRUE.
      IF (NPHOTI.EQ.0.AND.NITER.GE.1) NIPLS=3      ! Include D and D2 BGK bulk "ion" species...
      IF (NPHOTI.EQ.0.AND.NITER.EQ.0) THEN
        NIPLS=NPLSI  
        BULK_SOURCES = .TRUE.
      ENDIF 

      DO IPLS=1,NIPLS          ! Photons causing crash at COPV at the moment... this is caused by 
c      DO IPLS=1,MIN(NPLSI,1)  ! specifying NITER > 0 when not BGK or photons...

        IPLSTI=MPLSTI(IPLS)
        IPLSV =MPLSV (IPLS)
        ICOUNT=39
        ITALLY=16

        WRITE(logfp,*) 'COPV:',ipls,lcopv 
        WRITE(logfp,*) 'LCOPV:',lcopv

        WRITE(FP,80) '* BULK PARTICLES - VOLUME TALLIES',IPLS
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NTRII
        WRITE(FP,83) (0,I1=1,ITALLY)
        SUMION=0.0D0
        DO IR=1,NTRII
          DDUM=0.0D0
          IF (BULK_SOURCES(IPLS)) THEN                                   ! Only calculated for legit plasma species, IPLS=1 for now...
c...        Ionisation source from test particles [s-1]:
            DDUM(1)=PAPL (IPLS,IR)*VOL(IR)/ELCHA     ! Atoms
            DDUM(2)=PMPL (IPLS,IR)*VOL(IR)/ELCHA     ! Molecules
            DDUM(3)=PIPL (IPLS,IR)*VOL(IR)/ELCHA     ! Ions
            IF (NPHOTI.GT.0) DDUM(4)=PPHPL(IPLS,IR)*VOL(IR)/ELCHA     ! Photons
            DDUM(5)=DDUM(1)+DDUM(2)+DDUM(3)+DDUM(4)  ! Total
            SUMION=SUMION+DDUM(5)                    ! *TEMP*
c...        Volume recombination sink [???]:
            RECPAR=0.0D0
            RECMOM=0.0D0
            RECENG=0.0D0
            RECELE=0.0D0
            DO IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              IF (NSTORDR >= NRAD) THEN
                RECADD=TABRC1(IRRC,IR)*DIIN(IPLS,IR)  ! TURNING THIS *OFF* BECAUSE I DON'T KNOW WHERE IT CAME FROM...
                EEADD= EELRC1(IRRC,IR)*DIIN(IPLS,IR)
              ELSE
                RECADD=FTABRC1(IRRC,IR)*DIIN(IPLS,IR)
                EEADD= FEELRC1(IRRC,IR)*DIIN(IPLS,IR)
              ENDIF
              RECPAR=RECPAR+RECADD
              RECMOM=RECMOM+RECADD*PARMOM(IPLS,IR)
              RECENG=RECENG+RECADD*
     .                      (1.5D0*TIIN(IPLSTI,IR)+EDRIFT(IPLS,IR))
              RECELE=RECELE+EEADD
            ENDDO
            DDUM(6)=RECPAR*VOL(IR)
            RECENG = 0.0D0
            RECELE = 0.0D0
c...        Momentum source [???]:
            IADD=0                   ! Having IADD.NE.0 in couple_B2.5.f gives improved statistics? Ask...
            INC=NCLTAL(IR)
            SIGNUM=SIGN(1.0D0,BVIN(IPLSV,IR))
c           ***HACK***
c           Need to keep IPLS=1 or this dies...
            DDUM(7)=COPV(IADD+MIN(1,IPLS),INC)*VOL(IR)*1.0D-05*
     .              SIGNUM/ELCHA ! Was VOLTAL(IR)
c            DDUM(7)=COPV(IADD+IPLS,INC)*VOL(IR)*1.0D-05*SIGNUM/ELCHA ! Was VOLTAL(IR)
c            WRITE(6,'(A,2I6,1P,2E10.2)') 
c     .        'COPV:',inc,ir,COPV(IADD+IPLS,INC),vol(ir)
c            DDUM(7)=(COPV(IADD+IPLS,INC)+RECMOM)*VOL(IR)*1.D-5*SIGNUM/ELCHA 
            IF (DABS(DDUM(7)).LT. 1.0D-37) DDUM(7)= 0.0D+00 ! Make pretty
            IF (     DDUM(7) .LT.-1.0D+37) DDUM(7)=-1.0D+37
            IF (     DDUM(7) .GT. 1.0D+37) DDUM(7)= 1.0D+37     
c...        Energy sources [???]:
            DDUM(8)=(EAPL(IR)+EMPL(IR)+EIPL(IR)+RECENG)*VOL(IR)!*ELCHA  ! Ions       RECENG=0.0 ABOVE!
            DDUM(9)=(EAEL(IR)+EMEL(IR)+EIEL(IR)+RECELE)*VOL(IR)!*ELCHA  ! Elections  RECELE=0.0 ABOVE!
          ENDIF
c...      Plasma quantities:
          DDUM(10)=DIIN(IPLS,IR)*1.0D+06                             ! Density [m-3]
          DDUM(11)=VXIN(IPLSV,IR)                                    ! Velocity [...]
          DDUM(12)=VYIN(IPLSV,IR)
          DDUM(13)=VZIN(IPLSV,IR)
          DDUM(14)=BVIN(IPLSV,IR)                                    ! Velocity parallel to B[...]
          DDUM(15)=TIIN(IPLSTI,IR)                                   ! Temperature (yes?) [...]
          DDUM(16)=EDRIFT(IPLS,IR)                                   ! Drift energy [...]
c...      Output:
          ICOUNT=ICOUNT+1
          IF (ICOUNT.EQ.40) THEN
            WRITE(FP,'(A,7X,20(I14))') '*',(I1,I1=1,ITALLY)
            ICOUNT=0
          ENDIF
          WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
        ENDDO
c...    Surfaces fluxes:
        ICOUNT=39
        ITALLY=3
        WRITE(FP,80) '* BULK PARTICLES - SURFACE FLUXES',IPLS
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NSUR
        WRITE(FP,85) (I1,I1=1,ITALLY)
        CONV = 1.602176D-19
        DO IR=1,NTRII 
          DO IS=1,NSIDE
            IF (INSPAT(IS,IR).EQ.0) CYCLE
            DDUM=0.0D0
c...        Output:
            MSURFG=NLIM+NSTS+INSPAT(IS,IR)
            DDUM(1)=DBLE(IS)                  ! Side index of the triangle
            DDUM(2)=POTPL(IPLS,MSURFG)/CONV   ! D+ particle flux (s-1)
            DDUM(3)=EOTPL(IPLS,MSURFG)/CONV   ! D+ energy flux (eV s-1)
            ICOUNT=ICOUNT+14
            IF (ICOUNT.EQ.40) THEN
              WRITE(FP,'(A,11X,20(I12))') '*',(I1,I1=1,ITALLY)
              ICOUNT=0
            ENDIF
            WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
          ENDDO
        ENDDO
      ENDDO
 80   FORMAT(A,20(I8))
 81   FORMAT(20(I8))
 82   FORMAT(I8,1P,20(E14.6))
 83   FORMAT(6X,1P,20(I14))
 84   FORMAT(2I8,1P,20(E14.6))
 85   FORMAT(12X,1P,20(I14))

      WRITE(logfp,*) 'SUMION:',sumion,sumion*ELCHA  ! *TEMP*
c
c     ----------------------------------------------------------------------
c...  Test atoms:
      WRITE(logfp,80) '* TEST ATOMS',IATM  ! *TEMP*
      DO IATM=1,NATMI
c...    Volume averaged tallies:
        ICOUNT=39
        ITALLY=7
        WRITE(FP,80) '* TEST ATOMS - VOLUME TALLIES',IATM
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NTRII
        WRITE(FP,83) (I1,I1=1,ITALLY)
        DO IR=1,NTRII 
          DDUM=0.0D0
          DDUM(1)=PDENA(IATM,IR)*VOL(IR)                 ! Density [particles]
          DDUM(2)=0.0D0                                  ! Velocity [...]
          DDUM(3)=0.0D0
          DDUM(4)=0.0D0
          DDUM(5)=0.0D0                                  ! Velocity parallel to B [...]
          DDUM(6)=EDENA(IATM,IR)*VOL(IR)                 ! Energy [eV]
          DDUM(7)=EDENA(IATM,IR)/(PDENA(IATM,IR)+EPS10)  ! Energy per particle [eV]
c...      Output:
          ICOUNT=ICOUNT+1
          IF (ICOUNT.EQ.40) THEN
            WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
            ICOUNT=0
          ENDIF
          WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
        ENDDO
c...    Surfaces fluxes:
        ICOUNT=39
        ITALLY=8
        WRITE(FP,80) '* TEST ATOMS - SURFACE FLUXES',IATM
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NSUR
        WRITE(FP,85) (I1,I1=1,ITALLY)
        CONV = 1.602176D-19
        DO IR=1,NTRII 
          DO IS=1,NSIDE
            IF (INSPAT(IS,IR).EQ.0) CYCLE
            DDUM=0.0D0
c...        Output:
            MSURFG=NLIM+NSTS+INSPAT(IS,IR)
            DDUM(1)=DBLE(IS)                   ! Side index of the triangle
            DDUM(2)=POTAT  (IATM,MSURFG)/CONV  ! Incident atom particle flux (s-1)
            DDUM(3)=EOTAT  (IATM,MSURFG)/CONV  ! Incident atom energy   flux (eV s-1)
            DDUM(4)=PRFAAT (IATM,MSURFG)/CONV  ! Emitted atom flux from incident atoms     (s-1)
            DDUM(5)=PRFMAT (IATM,MSURFG)/CONV  ! Emitted atom flux from incident mols.     (s-1)
            DDUM(6)=PRFIAT (IATM,MSURFG)/CONV  ! Emitted atom flux from incident test ions (s-1)
c            DDUM(7)=PRFPHAT(IATM,MSURFG)/CONV  ! Emitted atom flux from incident photons   (s-1)
            DDUM(7)=PRFPAT (IATM,MSURFG)/CONV  ! Emitted atom flux from incident bulk ions (s-1)
            DDUM(8)=ERFAAT (IATM,MSURFG)/CONV  ! ???
            ICOUNT=ICOUNT+14
            IF (ICOUNT.EQ.40) THEN
              WRITE(FP,'(A,11X,20(I12))') '*',(I1,I1=1,ITALLY)
              ICOUNT=0
            ENDIF
            WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
          ENDDO
        ENDDO
      ENDDO
c
c     ----------------------------------------------------------------------=
c...  Test molecules:
      ICOUNT=39
      ITALLY=7
      WRITE(logfp,80) '* TEST MOLECULES',IMOL
      DO IMOL=1,NMOLI
        WRITE(FP,80) '* TEST MOLECULES - VOLUME TALLIES',IMOL
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NTRII
        WRITE(FP,83) (I1,I1=1,ITALLY)
        DO IR=1,NTRII 
          DDUM=0.0D0
          DDUM(1)=PDENM(IMOL,IR)*VOL(IR)                 ! Density [particles]
          DDUM(2)=0.0D0                                  ! Velocity [...]
          DDUM(3)=0.0D0
          DDUM(4)=0.0D0
          DDUM(5)=0.0D0                                  ! Velocity parallel to B [...]
          DDUM(6)=EDENM(IMOL,IR)*VOL(IR)                 ! Energy [eV]
          DDUM(7)=EDENM(IMOL,IR)/(PDENM(IMOL,IR)+EPS10)  ! Energy per particle [eV]
c...      Output:
          ICOUNT=ICOUNT+1
          IF (ICOUNT.EQ.40) THEN
            WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
            ICOUNT=0
          ENDIF
          WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
        ENDDO
c...    Surfaces fluxes:
        ICOUNT=39
        ITALLY=8
        WRITE(FP,80) '* TEST MOLECULES - SURFACE FLUXES',IMOL
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NSUR
        WRITE(FP,85) (I1,I1=1,ITALLY)
        CONV = 1.602176D-19
        DO IR=1,NTRII 
          DO IS=1,NSIDE
            IF (INSPAT(IS,IR).EQ.0) CYCLE
            DDUM=0.0D0
c...        Output:
            MSURFG=NLIM+NSTS+INSPAT(IS,IR)
            DDUM(1)=DBLE(IS)                   ! Side index of the triangle
            DDUM(2)=POTML  (IMOL,MSURFG)/CONV  ! Incident molecule particle flux (s-1)
            DDUM(3)=EOTML  (IMOL,MSURFG)/CONV  ! Incident molecule energy   flux (eV s-1)
            DDUM(4)=PRFAML (IMOL,MSURFG)/CONV  ! Emitted molecule flux from incident atoms     (s-1)
            DDUM(5)=PRFMML (IMOL,MSURFG)/CONV  ! Emitted molecule flux from incident mols.     (s-1)
            DDUM(6)=PRFIML (IMOL,MSURFG)/CONV  ! Emitted molecule flux from incident test ions (s-1)
c            DDUM(7)=PRFPHML(IMOL,MSURFG)/CONV  ! Emitted molecule flux from incident photons   (s-1)
            DDUM(7)=PRFPML (IMOL,MSURFG)/CONV  ! Emitted molecule flux from incident bulk ions (s-1)
            DDUM(8)=ERFAML (IMOL,MSURFG)/CONV  ! ???
            ICOUNT=ICOUNT+14
            IF (ICOUNT.EQ.40) THEN
              WRITE(FP,'(A,11X,20(I12))') '*',(I1,I1=1,ITALLY)
              ICOUNT=0
            ENDIF
            WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
          ENDDO
        ENDDO
      ENDDO
c
c     ----------------------------------------------------------------------
c...  Test ions:
      ICOUNT=39
      ITALLY=9
c...  Load rate coefficient data if calculating H2+ from CR-model:
      IF (NIONI.NE.1.OR.NMOLI.NE.1) THEN
        WRITE(logfp,*) 'ERROR: POSSIBLE PROBLEM IN OUTUSR.  IT IS '
        WRITE(logfp,*) '       ASSUMED THAT IMOL=1 IS D2 and IION=1'
        WRITE(logfp,*) '       IS D2+ FOR THE CALCULATION OF THE '
        WRITE(logfp,*) '       D2+ DENSITY FROM THE CR-MODEL.  THIS'
        WRITE(logfp,*) '       ASSUMPTION NEEDS TO BE REMOVED.'
        STOP
      ENDIF
      FILNAM='AMJUEL  '
      H123  ='H.12'
      REAC  ='2.0c     '
      CRC   ='OT '
      FPRM = 0._DP
      RCMIN = -HUGE(1._DP)
      RCMAX =  HUGE(1._DP)
      JFEXMN = 0
      JFEXMX = 0
      REACDAT(NREACI+1)%LOTH = .FALSE.
      CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .            RCMIN, RCMAX, FPRM, JFEXMN, JFEXMX,'  ',0)
      DO I=1,9
        DO J=1,9
          RH2PH2(I-1,J-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(I,J)
        ENDDO
      ENDDO
c      CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC)
c      DO I=1,9
c        DO J=1,9
c          RH2PH2(I-1,J-1)=CREAC(I,J,NREACI+1)
c        ENDDO
c      ENDDO
c...  Loop over grid elements:
      WRITE(logfp,80) '* TEST IONS - VOLUME TALLIES',IION
      DO IION=1,NIONI
        WRITE(FP,80) '* TEST IONS - VOLUME TALLIES',IION
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NTRII
        WRITE(FP,83) (I1,I1=1,ITALLY)
        DO IR=1,NTRII 
          DDUM=0.0D0
c...      Some work to get n_H2+ from CR-model rates:
          IF (.NOT.LGVAC(IR,NPLS+1)) THEN
c          IF (.NOT.LGVAC(IR,0)) THEN
c          IF (.NOT.LGVAC(IR,0).AND.NITER.LE.1) THEN
            DEF=LOG(DEIN(IR)*1.0D-08+EPS10)
            TEF=LOG(TEIN(IR)        +EPS10)
            RATIO=0.0D0
            DO J=0,8
              DEJ=DEF**J
              DO I=0,8
                TEI=TEF**I
                RATIO=RATIO+RH2PH2(I,J)*TEI*DEJ
              ENDDO
            ENDDO
            IF (RATIO.LT.80.0D0) THEN
              RATIO=DEXP(RATIO)
            ELSE
              IF (OUTPUT1) THEN
                WRITE(logfp,*) 'WARNING: EXCESSIVE RATIO IN OUTUSR WHEN'
                WRITE(logfp,*) '         CALCULATING D2+ DENSITY'
                WRITE(6,*) 'WARNING: EXCESSIVE RATIO IN OUTUSR WHEN'
                WRITE(6,*) '         CALCULATING D2+ DENSITY'
                OUTPUT1=.FALSE.
              ENDIF
              PROB1=MIN(IR,PROB1)
              PROB2=MAX(IR,PROB2)
              RATIO=0.0D0
            ENDIF
          ELSE
            RATIO=0.0D0
          ENDIF
c          WRITE(6,*) 'RATIO:',RATIO,LGVAC(IR,0)
c          RATIO=1.0D0
          IMOL=1
          DDUM(1)=PDENM(IMOL,IR)*VOL(IR)*RATIO           ! Density from CR-model [particles]
          DDUM(2)=RATIO                                  ! ...temp...
c...      Done determining n_H2+, load remaining quantities:
          DDUM(3)=PDENI(IION,IR)*VOL(IR)                 ! Density (statistical)
          DDUM(4)=0.0D0                                  ! Velocity [...]
          DDUM(5)=0.0D0
          DDUM(6)=0.0D0
          DDUM(7)=0.0D0                                  ! Velocity parallel to B [...]
          DDUM(8)=EDENI(IION,IR)*VOL(IR)                 ! Energy [eV]
          DDUM(9)=EDENI(IION,IR)/(PDENI(IION,IR)+EPS10)  ! Energy per particle [eV]
c...      Output:
          ICOUNT=ICOUNT+1
          IF (ICOUNT.EQ.40) THEN
            WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
            ICOUNT=0
          ENDIF
          WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
        ENDDO
c...    Surfaces fluxes:
        ICOUNT=39
        ITALLY=3
        WRITE(FP,80) '* TEST IONS - SURFACE FLUXES',IION
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NSUR
        WRITE(FP,85) (I1,I1=1,ITALLY)
        CONV = 1.602176D-19
        DO IR=1,NTRII 
          DO IS=1,NSIDE
            IF (INSPAT(IS,IR).EQ.0) CYCLE
            DDUM=0.0D0
c...        Output:
            MSURFG=NLIM+NSTS+INSPAT(IS,IR)
            DDUM(1)=DBLE(IS)                  ! Side index of the triangle
            DDUM(2)=POTIO(IION,MSURFG)/CONV   ! Ion particle flux (s-1)
            DDUM(3)=EOTIO(IION,MSURFG)/CONV   ! Ion energy flux (eV s-1)
            ICOUNT=ICOUNT+14
            IF (ICOUNT.EQ.40) THEN
              WRITE(FP,'(A,11X,20(I12))') '*',(I1,I1=1,ITALLY)
              ICOUNT=0
            ENDIF
            WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
          ENDDO
        ENDDO
      ENDDO
c
c     ----------------------------------------------------------------------
c...  Photons:
      WRITE(logfp,80) '* PHOTONS',NPHOTI
      DO IPHOT=1,NPHOTI
c...    Volume averaged tallies:
        ICOUNT=39
        ITALLY=3
        WRITE(FP,80) '* TEST PHOTONS - VOLUME TALLIES',IPHOT
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NTRII
        WRITE(FP,83) (I1,I1=1,ITALLY)
        DO IR=1,NTRII 
          DDUM=0.0D0
          DDUM(1)=PDENPH(IPHOT,IR)*VOL(IR)                   ! Density [particles]
          DDUM(2)=EDENPH(IPHOT,IR)*VOL(IR)                   ! Energy [eV]
          DDUM(3)=EDENPH(IPHOT,IR)/(PDENPH(IPHOT,IR)+EPS10)  ! Energy per particle [eV]
c          DDUM(2)=0.0D0                                      ! Velocity [...]
c          DDUM(3)=0.0D0
c          DDUM(4)=0.0D0  
c          DDUM(5)=0.0D0                                      ! Velocity parallel to B [...]
c          DDUM(6)=EDENPH(IPHOT,IR)*VOL(IR)                   ! Energy [eV]
c          DDUM(7)=EDENPH(IPHOT,IR)/(PDENPH(IPHOT,IR)+EPS10)  ! Energy per particle [eV]
c...      Output:
          ICOUNT=ICOUNT+1
          IF (ICOUNT.EQ.40) THEN
            WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
            ICOUNT=0
          ENDIF
          WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
        ENDDO
      ENDDO
c
c     ----------------------------------------------------------------------
c...  Balmer lines:
      ICOUNT=39
      ITALLY=7
c...  Halpha:
      ADDV=0.0D0
      CALL HALPHA(0,1,2,3,4,5,7,6)
c      DO IR=1,NTRII 
c        DO IADV=1,7
c          ADDV2(IADV,IR)=ADDV(IADV,IR)  ! Shorthand for this..?
c        ENDDO
c      ENDDO
      WRITE(FP,80) '* LINE EMISSION: HALPHA'
      WRITE(FP,81) ITALLY
      WRITE(FP,81) NTRII
      WRITE(FP,83) (I1,I1=1,ITALLY)
      DO IR=1,NTRII 
        DDUM=0.0D0
c        DO IADV=1,7
c          DDUM(IADV)  =ADDV2(IADV,IR)*VOL(IR)  ! Dalpha [photons s-1]
c        ENDDO
        DDUM(1:7)=ADDV(1:7,IR)*VOL(IR)  ! Dalpha [photons s-1]
c...    Output:
        ICOUNT=ICOUNT+1
        IF (ICOUNT.EQ.40) THEN
          WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
          ICOUNT=0
        ENDIF
        DO I1=1,ITALLY
          IF (DDUM(I1).LT.1.0D-30) DDUM(I1)=1.0D-30
          IF (DDUM(I1).GT.1.0D+30) DDUM(I1)=1.0D+30
        ENDDO
        WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
      ENDDO
c...  Hgamma:
      ICOUNT=39
      ADDV=0.0D0
      CALL HGAMMA(0,1,2,3,4,5,6)
c      ADDV2(IADV,IR)=ADDV(IADV,IR)
c      DO IR=1,NTRII 
c        DO IADV=1,6
c
c        ENDDO
c      ENDDO
      WRITE(FP,80) '* LINE EMISSION: HGAMMA'
      WRITE(FP,81) ITALLY
      WRITE(FP,81) NTRII
      WRITE(FP,83) (I1,I1=1,ITALLY)
      DO IR=1,NTRII 
        DDUM=0.0D0
c        DO IADV=1,6
c          DDUM(IADV+7)=ADDV2(IADV,IR)*VOL(IR)  ! Dgamma [photons s-1]
c        ENDDO
        DDUM(1:6)=ADDV(1:6,IR)*VOL(IR)  ! Dalpha [photons s-1]
c...    Output:
        ICOUNT=ICOUNT+1
        IF (ICOUNT.EQ.40) THEN
          WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
          ICOUNT=0
        ENDIF
        DO I1=1,ITALLY
          IF (DDUM(I1).LT.1.0D-30) DDUM(I1)=1.0D-30
          IF (DDUM(I1).GT.1.0D+30) DDUM(I1)=1.0D+30
        ENDDO
        WRITE(FP,82) IR,(DDUM(I1),I1=1,ITALLY)
      ENDDO
c
c     ----------------------------------------------------------------------
c...  Misc. electron data, VOL, VOLTAL? (are they the same?), NCLTAL(IR)=IR?
      WRITE(logfp,80) '* MISC'
      IF (.TRUE.) THEN
        ITALLY=13
        ICOUNT=39
        WRITE(FP,80) '* MISC'
        WRITE(FP,81) ITALLY
        WRITE(FP,81) NTRII
        WRITE(FP,83) (I1,I1=1,ITALLY)
        DO IR=1,NTRII 
          DDUM=0.0D0
          DDUM(1) =DBLE(NCLTAL(IR))   ! NCLTAL
          DDUM(2) =DBLE(ICELLRD(IR))  ! Index 1?
          DDUM(3) =DBLE(IRINGRD(IR))  ! Index 2?
          DDUM(4) =VOL(IR)         
          DDUM(5) =VOLTAL(IR)
          DDUM(6) =DEIN(IR)*1.0D+06   ! Electrons [m-3]
          DDUM(7) =TEIN(IR)
          DDUM(8) =BXIN(IR)           ! B-field
          DDUM(9) =BYIN(IR)
          DDUM(10)=BZIN(IR)
          DDUM(11)=BFIN(IR)
          IATM=1                      ! *Total* ionisation rate output for IATM=1, assumed to be D...
          DO IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (ir.EQ.1) WRITE(logfp,*) 'TABDS1:',IAEI,IRDS  ! *TEMP*
            DDUM(11+IRDS)=DDUM(11+IRDS)+TABDS1(IRDS,IR)
          ENDDO 
c...      Output:
          ICOUNT=ICOUNT+1
          IF (ICOUNT.EQ.40) THEN
            WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
            ICOUNT=0
          ENDIF
          WRITE(FP,90) IR,(DDUM(I1),I1=1,ITALLY)
        ENDDO
      ENDIF
 90   FORMAT(I8,3F12.0,1P,20(E12.4))
c
c     ----------------------------------------------------------------------
c...  Particle source:
      WRITE(logfp,80) '* PARTICLE SOURCES'
      ICOUNT=39
      WRITE(FP,80) '* PARTICLE SOURCES'
      WRITE (FP,'(I6)') NSTRAI
      DO IS=1,NSTRAI
        WRITE (FP,'(I6,I8,1P,3E14.6,0P)') 
     .    IS,SAVE_IPANU(IS),FLUXT(IS),PTRASH(IS),ETRASH(IS)
      ENDDO
c
c     ----------------------------------------------------------------------
c...  Pumped fluxes - a bit of a hack, but the easiest way to do things 
c     at the moment, pulled from output.f:
      WRITE(logfp,80) '* PUMPED FLUX'
      ICOUNT=39
      WRITE(FP,80) '* PUMPED FLUX'
      IF (LSPUMP) THEN
        DO J=1,NLIMPS
          JJ=J
          IF (J.GT.NLIM) JJ=-(J-NLIM)
          DO IS=1,NSPTOT
            IF (SPUMP(IS,J).GT.0.0D0) THEN
              WRITE (FP,91) JJ,''''//TEXTS(IS)//'''',SPUMP(IS,J)
            ENDIF
          ENDDO
        ENDDO
      ELSE
      ENDIF
c
c ----------------------------------------------------------------------
c...  End of volume/surface data for sum over strata marker:
      WRITE(FP,80) '* DONE (SUM OVER STRATA)'
c
c ----------------------------------------------------------------------
c...  Insert iteration data, if any:
      ITALLY=7
      FACT=ELCHA*0.66667D0*7.502D0*1.0D+06
      WRITE(logfp,80) '* ITERATION DATA'
      WRITE(logfp,*) 'ITERNO?',iterno
      DO IT=1,ITERNO-1
        WRITE(FP,80) '* ITERATION DATA',IT
        WRITE(FP,81) NGAUGE
        DO IG=1,NGAUGE
          ICOUNT=39
          WRITE(FP,80) '*   PRESSURE GAUGE',IG
          WRITE(FP,81) ITALLY
          WRITE(FP,81) NSTRAI+1
          WRITE(FP,83) (I1,I1=1,ITALLY)
          DO IS=1,NSTRAI+1
            DDUM=0.0D0
            IF (IS.EQ.NSTRAI+1) THEN
              DDUM(1)=T_VOL(IG)
c             Atoms:
              DDUM(2)=SUM(T_EDENA(1:NSTRAI,IG,IT))*FACT     ! mTorr
              DDUM(3)=SUM(T_PDENA(1:NSTRAI,IG,IT))*1.0E+06  ! particles m-3
              DDUM(4)=SUM(T_EDENA(1:NSTRAI,IG,IT))*1.0E+06  ! eV m-3
c             Molecules:
              DDUM(5)=SUM(T_EDENM(1:NSTRAI,IG,IT))*FACT     
              DDUM(6)=SUM(T_PDENM(1:NSTRAI,IG,IT))*1.0E+06  
              DDUM(7)=SUM(T_EDENM(1:NSTRAI,IG,IT))*1.0E+06  
            ELSE
              DDUM(1)=T_VOL(IG)
c             Atoms:
              DDUM(2)=T_EDENA(IS,IG,IT)*FACT
              DDUM(3)=T_PDENA(IS,IG,IT)*1.0E+06
              DDUM(4)=T_EDENA(IS,IG,IT)*1.0E+06
c             Molecules:
              DDUM(5)=T_EDENM(IS,IG,IT)*FACT
              DDUM(6)=T_PDENM(IS,IG,IT)*1.0E+06
              DDUM(7)=T_EDENM(IS,IG,IT)*1.0E+06
            ENDIF
c...        Output:
            ICOUNT=ICOUNT+1
            IF (ICOUNT.EQ.40) THEN
              WRITE(FP,'(A,5X,20(I14))') '*',(I1,I1=1,ITALLY)
              ICOUNT=0
            ENDIF
            WRITE(FP,82) IS,(DDUM(I1),I1=1,ITALLY)
          ENDDO
        ENDDO
      ENDDO
cc ----------------------------------------------------------------------
cc...  Insert iteration data, if any:
c      ICOUNT=39
c      ITALLY=4
c      FACT=ELCHA*0.667D0*7.502D0*1.0D+06
c      WRITE(logfp,80) '* ITERATION DATA'
c      WRITE(logfp,*) 'ITERNO?',iterno
c      DO IT=1,ITERNO-1
c        WRITE(FP,80) '* ITERATION DATA',IT
c        DO IG=1,NGAUGE
c          WRITE(FP,80) '*   PRESSURE GAUGE',IG
c          WRITE(FP,81) ITALLY
c          WRITE(FP,81) NSTRAI+1
c          WRITE(FP,83) (I1,I1=1,ITALLY)
c          DO IS=1,NSTRAI+1
c            DDUM=0.0D0
c            IF (IS.EQ.NSTRAI+1) THEN
c              DDUM(1)=T_VOL(IG)
c              DDUM(2)=SUM(T_EDENM(1:NSTRAI,IG,IT))*FACT       ! mTorr
c              DDUM(3)=SUM(T_PDENM(1:NSTRAI,IG,IT))*T_VOL(IG)  ! Total particles 
c              DDUM(4)=SUM(T_EDENM(1:NSTRAI,IG,IT))*T_VOL(IG)  ! Total energy
c            ELSE
c              DDUM(1)=T_VOL(IG)
c              DDUM(2)=T_EDENM(IS,IG,IT)*FACT
c              DDUM(3)=T_PDENM(IS,IG,IT)*T_VOL(IG)
c              DDUM(4)=T_EDENM(IS,IG,IT)*T_VOL(IG)
c            ENDIF
cc...        Output:
c            ICOUNT=ICOUNT+1
c            IF (ICOUNT.EQ.40) THEN
c              WRITE(FP,'(A,7X,20(I12))') '*',(I1,I1=1,ITALLY)
c              ICOUNT=0
c            ENDIF
c            WRITE(FP,82) IS,(DDUM(I1),I1=1,ITALLY)
c          ENDDO
c        ENDDO
c      ENDDO
c
c ----------------------------------------------------------------------
c...  EOF marker:
      WRITE(FP,80) '* DONE (FOR GOOD THIS TIME)'

c...  Clear arrays:
      IF (ALLOCATED(T_VOL  )) DEALLOCATE(T_VOL  )
      IF (ALLOCATED(T_PDENA)) DEALLOCATE(T_PDENA)
      IF (ALLOCATED(T_EDENA)) DEALLOCATE(T_EDENA)
      IF (ALLOCATED(T_PDENM)) DEALLOCATE(T_PDENM)
      IF (ALLOCATED(T_EDENM)) DEALLOCATE(T_EDENM)


 91   FORMAT(I6,2X,A,2X,1P,E12.4,0P)


      ! TETRAHEDRON
      ntrii = 0

      RETURN
c...
c     ------------------------------------------------------------------
c
c      ENTRY OUTUS1(ISTRAA,ISTRAE,NEW_ITER)

      ENTRY OUTUS1(ISTRAA,ISTRAE,NEW_ITER,IPANU)

      WRITE(logfp,*) 'OUTUS1:',ISTRAA,ISTRAE,NEW_ITER,NTRII
      WRITE(logfp,*) 'OUTUS1:',fluxt(istraa)

      SAVE_IPANU(ISTRAA) = IPANU
           
c...  Calculate the center (or so) of each triangle:
      ALLOCATE(XCEN(NTRII))
      ALLOCATE(YCEN(NTRII))
      ALLOCATE(ZCEN(NTRII))
      XCEN = 0.0D0
      YCEN = 0.0D0
      DO IR=1,NTRII 
        X1=XTRIAN(NECKE(1,IR))
        Y1=YTRIAN(NECKE(1,IR))
        X2=XTRIAN(NECKE(2,IR))
        Y2=YTRIAN(NECKE(2,IR))
        X3=XTRIAN(NECKE(3,IR))
        Y3=YTRIAN(NECKE(3,IR))
C       CENTER OF GRAVITY IN TRIANGLE 1 AND TRIANGLE 2
        XCEN(IR)=(X1+X2+X3)/3._DP
        YCEN(IR)=(Y1+Y2+Y3)/3._DP         
        ZCEN(IR)=0.0D0
      ENDDO

c...  Check for triangle centers that are inside the "pressure
c     gauge" region:
      NGAUGE = 3
      XGAUGE(1) =  845.0D0  ! ITER midplane
      YGAUGE(1) =   62.0D0
      RGAUGE(1) =    8.0D0
      XGAUGE(2) =  475.0D0  ! ITER divertor, below dome
      YGAUGE(2) = -420.0D0
      RGAUGE(2) =   20.0D0
      XGAUGE(3) =    4.0D0
      YGAUGE(3) =   17.5D0
      RGAUGE(3) =    2.0D0

      IS=ISTRAA
      IT=ITERNO

      WRITE(logfp,*) '  CHECKING PRESSURE GAUGES',IS,IT

      IF (.NOT.ALLOCATED(T_VOL)) THEN
        ALLOCATE(T_VOL(NGAUGE))
        ALLOCATE(T_PDENA(NSTRAI,NGAUGE,NITER+1))
        ALLOCATE(T_EDENA(NSTRAI,NGAUGE,NITER+1))
        ALLOCATE(T_PDENM(NSTRAI,NGAUGE,NITER+1))
        ALLOCATE(T_EDENM(NSTRAI,NGAUGE,NITER+1))
      ENDIF

      T_VOL=EPS60       
      T_PDENA(IS,1:NGAUGE,IT)=0.0D0
      T_EDENA(IS,1:NGAUGE,IT)=0.0D0
      T_PDENM(IS,1:NGAUGE,IT)=0.0D0
      T_EDENM(IS,1:NGAUGE,IT)=0.0D0
      IMOL=1  ! Not reliable...
      DO IG=1,NGAUGE
        DO IR=1,NTRII
          DIST = DSQRT( (XCEN(IR)-XGAUGE(IG))**2 + 
     .                  (YCEN(IR)-YGAUGE(IG))**2)
          IF (DIST.LE.RGAUGE(IG)) THEN
            T_VOL  (   IG   )=T_VOL  (   IG   )+               VOL(IR)
            T_PDENA(IS,IG,IT)=T_PDENA(IS,IG,IT)+PDENA(IMOL,IR)*VOL(IR)
            T_EDENA(IS,IG,IT)=T_EDENA(IS,IG,IT)+EDENA(IMOL,IR)*VOL(IR)
            T_PDENM(IS,IG,IT)=T_PDENM(IS,IG,IT)+PDENM(IMOL,IR)*VOL(IR)
            T_EDENM(IS,IG,IT)=T_EDENM(IS,IG,IT)+EDENM(IMOL,IR)*VOL(IR)
          ENDIF
        ENDDO
        T_PDENA(IS,IG,IT)=T_PDENA(IS,IG,IT)/T_VOL(IG)
        T_EDENA(IS,IG,IT)=T_EDENA(IS,IG,IT)/T_VOL(IG)
        T_PDENM(IS,IG,IT)=T_PDENM(IS,IG,IT)/T_VOL(IG)
        T_EDENM(IS,IG,IT)=T_EDENM(IS,IG,IT)/T_VOL(IG)
        WRITE(logfp,*) 'GAUGE:',IS,IT,IG,
     .                 T_EDENA(IS,IG,IT),T_EDENM(IS,IG,IT)
      ENDDO

c...  Clear arrays:
      IF (ALLOCATED(XCEN)) DEALLOCATE(XCEN)
      IF (ALLOCATED(YCEN)) DEALLOCATE(YCEN)
      IF (ALLOCATED(ZCEN)) DEALLOCATE(ZCEN)

      IF (IS.EQ.NSTRAI) ITERNO=ITERNO+1

c      ENTRY OUTUS1(ISTRAA,ISTRAE,NEW_ITER,IPANU)
c
c      WRITE(logfp,*) 'OUTUS1:',ISTRAA,ISTRAE,NEW_ITER,NTRII
c      WRITE(logfp,*) 'OUTUS1:',fluxt(istraa)
c
c      SAVE_IPANU(ISTRAA) = IPANU
c           
cc...  Calculate the center (or so) of each triangle:
c      ALLOCATE(XCEN(NTRII))
c      ALLOCATE(YCEN(NTRII))
c      ALLOCATE(ZCEN(NTRII))
c      XCEN = 0.0D0
c      YCEN = 0.0D0
c      DO IR=1,NTRII 
c        X1=XTRIAN(NECKE(1,IR))
c        Y1=YTRIAN(NECKE(1,IR))
c        X2=XTRIAN(NECKE(2,IR))
c        Y2=YTRIAN(NECKE(2,IR))
c        X3=XTRIAN(NECKE(3,IR))
c        Y3=YTRIAN(NECKE(3,IR))
cC       CENTER OF GRAVITY IN TRIANGLE 1 AND TRIANGLE 2
c        XCEN(IR)=(X1+X2+X3)/3._DP
c        YCEN(IR)=(Y1+Y2+Y3)/3._DP         
c        ZCEN(IR)=0.0D0
c      ENDDO
c
cc...  Check for triangle centers that are inside the "pressure
cc     gauge" region:
c      NGAUGE = 3
c      XGAUGE(1) =   65.0D0
c      YGAUGE(1) =  -66.0D0
c      RGAUGE(1) =    2.0D0
c      XGAUGE(2) =  190.0D0
c      YGAUGE(2) =    0.0D0
c      RGAUGE(2) =   10.0D0
c      XGAUGE(3) =  475.0D0  ! ITER divertor, below dome
c      YGAUGE(3) = -420.0D0
c      RGAUGE(3) =   30.0D0
c
c      IS=ISTRAA
c      IT=ITERNO
c
c      WRITE(logfp,*) '  CHECKING PRESSURE GAUGES',IS,IT
c
c      IF (.NOT.ALLOCATED(T_VOL)) THEN
c        ALLOCATE(T_VOL(NGAUGE))
c        ALLOCATE(T_PDENM(NSTRAI,NGAUGE,NITER+1))
c        ALLOCATE(T_EDENM(NSTRAI,NGAUGE,NITER+1))
c      ENDIF
c
c      T_VOL=0.0D0       
c      T_EDENM(IS,1:NGAUGE,IT)=0.0D0
c      IMOL=1  ! Not reliable...
c      DO IG=1,NGAUGE
c        DO IR=1,NTRII
c          DIST = DSQRT( (XCEN(IR)-XGAUGE(IG))**2 + 
c     .                  (YCEN(IR)-YGAUGE(IG))**2)
c          IF (DIST.LE.RGAUGE(IG)) THEN
c            T_VOL  (   IG   )=T_VOL  (   IG   )+               VOL(IR)
c            T_PDENM(IS,IG,IT)=T_PDENM(IS,IG,IT)+PDENM(IMOL,IR)*VOL(IR)
c            T_EDENM(IS,IG,IT)=T_EDENM(IS,IG,IT)+EDENM(IMOL,IR)*VOL(IR)
c          ENDIF
c        ENDDO
c        T_PDENM(IS,IG,IT)=T_PDENM(IS,IG,IT)/T_VOL(IG)
c        T_EDENM(IS,IG,IT)=T_EDENM(IS,IG,IT)/T_VOL(IG)
c        WRITE(logfp,*) 'GAUGE:',IS,IT,IG,T_EDENM(IS,IG,IT)
c      ENDDO
c
cc...  Clear arrays:
c      IF (ALLOCATED(XCEN)) DEALLOCATE(XCEN)
c      IF (ALLOCATED(YCEN)) DEALLOCATE(YCEN)
c      IF (ALLOCATED(ZCEN)) DEALLOCATE(ZCEN)
c
c      IF (IS.EQ.NSTRAI) ITERNO=ITERNO+1
c
      WRITE(logfp,*) '  DONE'


      RETURN
 98   WRITE(0,*) 'ERROR: UNABLE TO OPEN DATA TRANSFER FILE'
      WRITE(6,*) 'ERROR: UNABLE TO OPEN DATA TRANSFER FILE'
 99   STOP
      END
C ===== SOURCE: plausr.f


      SUBROUTINE PLAUSR
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CSTEP
      USE CGRID
      USE CINIT
      USE COMSOU
      USE CTETRA  ! TET diff
      IMPLICIT NONE
      REAL(DP) :: FACTOR, STEP
      INTEGER :: NLINES, ITET, ISIDE, I, NBIN, ISTRA, ISRFS, ISOR, 
     .           ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS
      INTEGER :: IDEZ
      INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
      REAL(DP) :: FLX, TE, TI, DE, DELR, FL
      CHARACTER(100) :: ZEILE
c slmod begin
      INTEGER :: I1,I2,I3
      REAL(DP) :: ARTRI3
      INTEGER ITSIDE(3,4)
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
c slmod end

      WRITE(0,*) 'DEBUG: HERE IN PLAUSR'

      ZEILE='*   '      
      DO WHILE (ZEILE(1:1) == '*')
         READ (31,'(A100)',END=99) ZEILE
      END DO

      READ (ZEILE,*) NLINES

      CALL ALLOC_CSTEP
      ALLOCATE (KSTEP(NSTEP))
      ALLOCATE (INOSRC(NSTEP))
      ALLOCATE (IPLAN(NSTEP))
      ALLOCATE (IPLEN(NSTEP))
      KSTEP = 0
      INOSRC = 0
      IPLAN = 0
      IPLEN = 0

      DO I=1, NLINES
        READ (31,*) ITET, ISIDE, FLX, TE, TI, DE
        STRATA:DO ISTRA=1,NSTRAI
          IF (.NOT.NLSRF(ISTRA)) CYCLE
          DO ISRFS=1,NSRFSI(ISTRA)
            ISOR=SORLIM(ISRFS,ISTRA)
            ITEC1=IDEZ(ISOR,1,4)
            ITEC2=IDEZ(ISOR,2,4)
            ITEC3=IDEZ(ISOR,3,4)
c slmod begin - debug
c            IF (istra.EQ.2) THEN
c              WRITE(0,*) 'ITEC1-3:',itec1,itec2,itec3
c            ENDIF
c slmod end            
            IF ((ITEC1 /= 4).AND.(ITEC2 /= 4).AND.(ITEC3 /= 4)) CYCLE
            ISTEP=SORIND(ISRFS,ISTRA)
            IF (ISTEP.EQ.0) THEN
              WRITE (6,*) 'ERROR IN PRIMARY SOURCE DATA '
              WRITE (6,*) 'STEPFUNCTION REQUESTED FOR SOURCE SURFACE '
              WRITE (6,*) 'NO. ', INSOR(ISRFS,ISTRA),' BUT SORIND.EQ.0.'
              CALL EXIT(1)
            ELSEIF (ISTEP.GT.NSTEP) THEN
              CALL MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
              CALL EXIT(1)
            ENDIF
            INDSRF=INSOR(ISRFS,ISTRA)
            IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)
            IF (INMTIT(ISIDE,ITET) == INDSRF) THEN
              IF (KSTEP(ISTEP) == 0) RRSTEP(ISTEP,1) = 0._DP
              KSTEP(ISTEP) = KSTEP(ISTEP) + 1
              INOSRC(ISTEP) = ISTRA
              IS1 = ISIDE + 1
              IF (IS1.GT.3) IS1=1
              IF (NSPEZ(ISTRA) <= 0) THEN
                IPLAN(ISTEP)=1
                IPLEN(ISTEP)=NPLSI
              ELSE
                IPLAN(ISTEP)=NSPEZ(ISTRA)
                IPLEN(ISTEP)=NSPEZ(ISTRA)
              END IF
              IRSTEP(ISTEP,KSTEP(ISTEP))=ITET
              IPSTEP(ISTEP,KSTEP(ISTEP))=ISIDE
              ITSTEP(ISTEP,KSTEP(ISTEP))=1
              IASTEP(ISTEP,KSTEP(ISTEP))=0
              IBSTEP(ISTEP,KSTEP(ISTEP))=1
c slmod begin - bug?
c... Not sure what to say here, but if only X,Y points are used then there are
c    cases where deltaX,Y is zero, i.e. the points are displaced from each other
c    in the Z direction only.  For the current case (m-tet-0003d) this appears to
c    happen only for the high index target, presumably because the points are 
c    ordered differently at the low index target, avoiding this problem (by chance).
c    What I'm not sure about is whether DELR has any physical meaning here and
c    whether or not it's important, i.e. should it be the area of the surface?
c
              I1=NTECK(ITSIDE(1,ISIDE),ITET)
              I2=NTECK(ITSIDE(2,ISIDE),ITET)
              I3=NTECK(ITSIDE(3,ISIDE),ITET)
              DELR=ARTRI3(XTETRA(I1),YTETRA(I1),ZTETRA(I1),
     .                    XTETRA(I2),YTETRA(I2),ZTETRA(I2),
     .                    XTETRA(I3),YTETRA(I3),ZTETRA(I3))
c              DELR =  SQRT(
c     .           (XTETRA(NTECK(ISIDE,ITET))-XTETRA(NTECK(IS1,ITET)))**2+
c     .           (YTETRA(NTECK(ISIDE,ITET))-YTETRA(NTECK(IS1,ITET)))**2+
c     .           (ZTETRA(NTECK(ISIDE,ITET))-ZTETRA(NTECK(IS1,ITET)))**2)
c
c              DELR =  SQRT(
c     .           (XTETRA(NTECK(ISIDE,ITET))-XTETRA(NTECK(IS1,ITET)))**2+
c     .           (YTETRA(NTECK(ISIDE,ITET))-YTETRA(NTECK(IS1,ITET)))**2)
c slmod end
              RRSTEP(ISTEP,KSTEP(ISTEP)+1)=
     .           RRSTEP(ISTEP,KSTEP(ISTEP)) + DELR         
              TESTEP(ISTEP,KSTEP(ISTEP)) = TE    
c slmod begin - debug
c              IF (istra.EQ.2) THEN
c                WRITE(0,'(A,3I9,2F10.2)') 
c     .            '       :',ipls,istep,kstep(istep),delr,ABS(flx)
c              ENDIF
c slmod end            
              DO IPLS=IPLAN(ISTEP), IPLEN(ISTEP)
                TISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = TI
                DISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = DE
                VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VXIN(IPLS,ITET) ! Approximation
                VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VYIN(IPLS,ITET)
                VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZIN(IPLS,ITET)
                FLSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ABS(FLX)/DELR
              END DO
              EXIT STRATA
            END IF
          END DO 
        END DO STRATA
      END DO

      DO ISTEP = 1, NSTEP
        IF (KSTEP(ISTEP) > 0) THEN
          NBIN=KSTEP(ISTEP)+1
          FL=STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP)
          FLUX(INOSRC(ISTEP))=FL
c slmod begin - debug
c          WRITE(0,*) 'FLUX:',istep,inosrc(istep),fl
c          WRITE(0,*) '    :',iplan(istep),iplen(istep),nbin
c slmod end
        END IF
      END DO

      DEALLOCATE (KSTEP)
      DEALLOCATE (INOSRC)

 99   RETURN
      END
C ===== SOURCE: pltusr.f
C
C
      SUBROUTINE PLTUSR(PLABLE,J)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE CCONA
      USE CLGIN
      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: PLABLE
      INTEGER, INTENT(IN) :: J
      INTEGER :: IB, MERK, IER
      REAL(DP) :: XS, YS, ZX0, ZY0, ZZ0, CX, CY, CZ, RZYLB, B0B, B1B,
     .          B2B, B3B, F0B, F1B, F2B, F3B, X0, Y0, B0, B1, B2, Z1, Z2
      INTEGER :: IN
      RETURN
      END
C ===== SOURCE: prousr.f
CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  TEXTOR  *  (HELIUM INCLUDED)
C           ************
C
      SUBROUTINE PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C     P0 : CENTRAL VALUE
C     P1 : STARTING RADIUS FOR POLYNOMIAL
C     P2 : SWITCH FOR PHASE 1: OH-PHASE
C                           2: NI-PHASE
C     P3 : FACTOR FOR TI: TI=K*TE
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CGRID
      USE CGEOM
      USE CINIT
      USE CCONA
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: NTR, I, J, LL
      REAL(DP), ALLOCATABLE, SAVE :: PLAS(:,:)
      REAL(DP) :: bnorm
      INTEGER, SAVE :: INDAR(12)=(/ (0, i=1,12) /) 


c from user_tri
      IF (.NOT.ALLOCATED(PLAS)) THEN

c        WRITE(0,*) 'DEBUG: LOADING DATA'
c slmod begin
        INDAR = 0
c slmod end

        LL=LEN_TRIM(CASENAME)

        FILENAME=CASENAME(1:LL) // '.plasma'
        OPEN (UNIT=31,FILE=FILENAME,ACCESS='SEQUENTIAL',
c slmod begin
c     .        FORM='FORMATTED',POSITION='REWIND')  
c                                                  
     .        FORM='FORMATTED')
c slmod end

        ZEILE='*   '      
        DO WHILE (ZEILE(1:1) == '*')
           READ (31,'(A100)') ZEILE
        END DO

        READ (ZEILE,*) NTR

        IF (NTR /= NR1ST-1) THEN
          WRITE (6,*) ' WRONG NUMBER OF TRIANGLES IN PLASMA FILE'
          WRITE (6,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
          CALL EXIT(1)
        END IF

!        ALLOCATE (PLAS(2+4*NPLS,NRAD))
!        ALLOCATE (PLAS(9,NRAD))
        ALLOCATE (PLAS(10,NRAD))  ! sl
        PLAS=0._DP
        DO I=1, NTR
          READ (31,*) J, PLAS(:,I)
          bnorm=sqrt(plas(7,i)**2 + plas(8,i)**2 + plas(9,I)**2)
          if (bnorm < eps12) then
             plas(7,i) = 0._dp
             plas(8,i) = 0._dp
             plas(9,i) = 1._dp
          else
             plas(7:9,i) = plas(7:9,i)/bnorm
          end if
        END DO
        plas(9,ntr+1:nrad) = 1._dp
      END IF
      
c      WRITE(0,*) 'DEBUG: PROUSR ',indx,npls
c      WRITE(0,*) 'DEBUG: INDAR  ',indar(1:9)
c slmod begin
      IF (INDX == -999 ) THEN
c        IF (ALLOCATED(PLAS)) DEALLOCATE(PLAS)
        INDAR=0 
      ELSEIF (INDX == -100 ) THEN
c...    Load ionisation data from previous run that included photon transport: 
        WRITE(88,*) 'PROUSR: LOADING TABDS1(2,*)'
        PRO(1:N) = PLAS(10,1:N)
      ELSEIF (INDX <= 1) THEN
c
c      IF (INDX <= 1) THEN
c slmod end
        PRO(1:N) = PLAS(INDX+1,1:N)
        INDAR(INDX+1) = INDAR(INDX+1)+1
      ELSEIF (INDX == 1+1*NPLS) THEN
        IF (INDAR(3) == 0) THEN
           PRO(1:N) = PLAS(3,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(3) = INDAR(3) + 1
      ELSEIF (INDX == 1+2*NPLS) THEN
        IF (INDAR(4) == 0) THEN
           PRO(1:N) = PLAS(4,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(4) = INDAR(4) + 1
      ELSEIF (INDX == 1+3*NPLS) THEN
        IF (INDAR(5) == 0) THEN
           PRO(1:N) = PLAS(5,1:N)     
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(5) = INDAR(5) + 1
      ELSEIF (INDX == 1+4*NPLS) THEN
        IF (INDAR(6) == 0) THEN
           PRO(1:N) = PLAS(6,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(6) = INDAR(6) + 1
      ELSEIF (INDX == 1+5*NPLS) THEN
        IF (INDAR(7) == 0) THEN
           PRO(1:N) = PLAS(7,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(7) = INDAR(7) + 1
      ELSEIF (INDX == 2+5*NPLS) THEN
        IF (INDAR(8) == 0) THEN
           PRO(1:N) = PLAS(8,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(8) = INDAR(8) + 1
      ELSEIF (INDX == 3+5*NPLS) THEN
        IF (INDAR(9) == 0) THEN
           PRO(1:N) = PLAS(9,1:N)
        ELSE
           PRO(1:N) = 1._DP
        END IF
        INDAR(9) = INDAR(9) + 1
      ELSE
        WRITE (6,*) ' PROUSR: NO DATA PROVIDED FOR INDEX ',INDX
        PRO(1:N) = 0._DP
      END IF

      RETURN
c end of code from user_tri
      
      IF (.NOT.ALLOCATED(PLAS)) THEN

        LL=LEN_TRIM(CASENAME)

        FILENAME=CASENAME(1:LL) // '.plasma'
        OPEN (UNIT=31,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')

        ZEILE='*   '      
        DO WHILE (ZEILE(1:1) == '*')
           READ (31,'(A100)') ZEILE
        END DO

        READ (ZEILE,*) NTR

        IF (NTR /= NR1ST-1) THEN
          WRITE (6,*) ' WRONG NUMBER OF TRIANGLES IN PLASMA FILE'
          WRITE (6,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
          CALL EXIT(1)
        END IF

!        ALLOCATE (PLAS(2+4*NPLS,NRAD))
!        ALLOCATE (PLAS(9,NRAD))
!        ALLOCATE (PLAS(10,NRAD))  ! sl
        ALLOCATE (PLAS(12,NRAD))  ! pb
        PLAS=0._DP
        DO I=1, NTR
          READ (31,*) J, PLAS(:,I)
          bnorm=sqrt(plas(7,i)**2 + plas(8,i)**2 + plas(9,I)**2)
          if (bnorm < eps12) then
             plas(7,i) = 0._dp
             plas(8,i) = 0._dp
             plas(9,i) = 1._dp
          else
             plas(7:9,i) = plas(7:9,i)/bnorm
          end if
        END DO
        plas(9,ntr+1:nrad) = 1._dp
      END IF
      
c      WRITE(0,*) 'PROUSR:',indx,npls
c slmod begin
      IF (INDX == -100 ) THEN
c...    Load ionisation data from previous run that included photon transport: 
        WRITE(0,*) 'PROUSR: LOADING TABDS1(2,*)'
        PRO(1:N) = PLAS(10,1:N)
      ELSEIF (INDX <= 1) THEN
c
c      IF (INDX <= 1) THEN
c slmod end
        PRO(1:N) = PLAS(INDX+1,1:N)
        INDAR(INDX+1) = INDAR(INDX+1)+1
      ELSEIF (INDX == 1+1*NPLS) THEN
        IF (INDAR(3) == 0) THEN
           PRO(1:N) = PLAS(3,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(3) = INDAR(3) + 1
      ELSEIF (INDX == 1+2*NPLS) THEN
        IF (INDAR(4) == 0) THEN
           PRO(1:N) = PLAS(4,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(4) = INDAR(4) + 1
      ELSEIF (INDX == 1+3*NPLS) THEN
        IF (INDAR(5) == 0) THEN
           PRO(1:N) = PLAS(5,1:N)     
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(5) = INDAR(5) + 1
      ELSEIF (INDX == 1+4*NPLS) THEN
        IF (INDAR(6) == 0) THEN
           PRO(1:N) = PLAS(6,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(6) = INDAR(6) + 1
      ELSEIF (INDX == 1+5*NPLS) THEN
        IF (INDAR(7) == 0) THEN
           PRO(1:N) = PLAS(7,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(7) = INDAR(7) + 1
      ELSEIF (INDX == 2+5*NPLS) THEN
        IF (INDAR(8) == 0) THEN
           PRO(1:N) = PLAS(8,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(8) = INDAR(8) + 1
      ELSEIF (INDX == 3+5*NPLS) THEN
        IF (INDAR(9) == 0) THEN
           PRO(1:N) = PLAS(9,1:N)
        ELSE
           PRO(1:N) = 1._DP
        END IF
        INDAR(9) = INDAR(9) + 1
      ELSEIF (INDX == 1+6*NPLS) THEN
        PRO(1:N) = PLAS(11,1:N)
        INDAR(11) = INDAR(11) + 1
      ELSEIF (INDX == 2+6*NPLS) THEN
        PRO(1:N) = PLAS(12,1:N)
        INDAR(12) = INDAR(12) + 1
      ELSE
        WRITE (6,*) ' PROUSR: NO DATA PROVIDED FOR INDEX ',INDX
        PRO(1:N) = 0._DP
      END IF

      RETURN
      END
C ===== SOURCE: refusr.f


      SUBROUTINE REFUSR
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW,XCW,XMP,XCP,ZCOS,ZSIN,EXPI,RPROB,
     .                        E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST
      ENTRY RF0USR
      ENTRY SPTUSR
      ENTRY SP0USR
      ENTRY SP1USR
      ENTRY RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
      RETURN
      END
C ===== SOURCE: retusr.f
c
c
      subroutine retusr(sig)
      USE PRECISION
      implicit none
      real(dp), intent(in) :: sig
      return
      end
C ===== SOURCE: samusr.f
C
C
      SUBROUTINE SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      IMPLICIT NONE
      REAL(DP) :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6
      REAL(DP) :: X0,Y0,Z0,TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,
     .                       WEISPZ
      INTEGER :: NLSF,is1, is2
      INTEGER :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IER

      entry sm0usr (is1,is2,sorad1,sorad2,sorad3,sorad4,sorad5,sorad6)
      return

      entry SM1USR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)

      RETURN
      END
C ===== SOURCE: sigusr.f


      SUBROUTINE SIGUSR(IFIRST,JJJ,ZDS,DUMMY1,PSIG,DUMMY2,ARGST,
     .                  XD0,YD0,ZD0,XD1,YD1,ZD1)
C
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H ALPHA FLUX PSIG(I),I=0,4 CONTRIBUTIONS
C          FROM ATOMS, MOLECULES, TEST IONS AND BULK IONS
C          THE INTEGRANT ARGST SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CSPEI
      USE COMSOU
      USE CLOGAU
      USE COMXS
      USE COMSIG
      USE CUPD
      USE CZT1
      USE CTRCEI
      USE CGRID
      USE CCONA
      USE CADGEO
      USE CLGIN
      USE CSDVI
      USE COUTAU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JJJ
      INTEGER, INTENT(INOUT) :: IFIRST
      REAL(DP), INTENT(INOUT) :: PSIG(0:NSPZ),ARGST(0:NSPZ,NRAD)
      REAL(DP), INTENT(IN) :: ZDS,DUMMY1,DUMMY2,XD0,YD0,ZD0,XD1,YD1,ZD1
      RETURN
      END
C ===== SOURCE: talusr.f
c
c
      subroutine talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE PRECISION
      USE PARMMOD
      USE CGRID
      USE CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: i

      ilast=1
      return 1
      end
C ===== SOURCE: timusr.f


      SUBROUTINE TIMUSR(N,X,Y,Z,VX,VY,VZ,N1,N2,T,IC,IE,NP,NL)
      USE PRECISION
      USE PARMMOD
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: X,Y,Z,VX,VY,VZ,T,cx,cy,cz,sc
      INTEGER, INTENT(IN) :: N, N1, N2, IC, IE, NP, IS, NRCELL
      LOGICAL :: NL

      ENTRY NORUSR(is,x,y,z,cx,cy,cz,sc,VX,VY,VZ,NRCELL)

      RETURN
      END
C ===== SOURCE: tmsusr.f


      SUBROUTINE TMSUSR (T0)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: T0
      RETURN
      END
C ===== SOURCE: upcusr.f
C
C
      SUBROUTINE UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE COMXS
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
C  FOR PARTICLE DENSITY IN CELL NO. NCELL
C     COLV(1,NCELL)=COLV(1,NCELL)+WS
C
      RETURN
      END
C ===== SOURCE: upnusr.f
c
c
      subroutine upnusr
      return
      end
C ===== SOURCE: upsusr.f
C
C
      SUBROUTINE UPSUSR(WT,IND)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WT
      INTEGER, INTENT(IN) :: IND
      RETURN
      END
C ===== SOURCE: uptusr.f
C
C
      SUBROUTINE UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CUPD
      USE COMXS
      USE CSPEZ
      USE CGRID
      USE CLOGAU
      USE CCONA
      USE CPOLYG
      USE CZT1
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: CNDYNA(NATM),CNDYNP(NPLS)
CDR
      REAL(DP) :: VPX(NRAD),VPY(NRAD),VRX(NRAD),VRY(NRAD)
CDR
      INTEGER :: IFIRST, IAT, IPL, I, IR, IP, IRD, IA1, IA2, IA3, NA4,
     .           INDEXM, INDEXF
      DATA IFIRST/0/
      SAVE

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        DO IAT=1,NATMI
          CNDYNA(IAT)=1.D3*AMUA*RMASSA(IAT)
        END DO
        DO IPL=1,NPLSI
          CNDYNP(IPL)=1.D3*AMUA*RMASSP(IPL)
        END DO
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        DO I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
        END DO
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
          END DO
        END DO
        IA1=NATMI+NMOLI
        IA2=2*IA1
        IA3=3*IA1
        NA4=4*IA1
        INDEXM=NPLSI
        INDEXF=2*NPLSI
      ENDIF
      RETURN
      END



C ===== SOURCE: usrtrc.f
      SUBROUTINE USRTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CSTEP
      USE CGRID
      USE CINIT
      USE COMSOU
      USE CTRIG
      USE COMPRT
      USE CUPD
      IMPLICIT NONE

      INTEGER  IFLAG,ISYM,count
      REAL(DP) XPLO,YPLO,ZPLO

      DATA count /0/
     
      IF (IFLAG.EQ.0) COUNT = COUNT + 1

      WRITE(90,'(I8,3F9.3,3I4,I4,I7,1P,1E10.2,0P,F8.4,I6,7I5,2X,5I5,
     .  1P,E10.2,0P)') 
     .  count,XPLO,YPLO,ZPLO,IFLAG,ISYM,istra,0,
c     .  count,XPLO,YPLO,ZPLO,IFLAG,ISYM,istra,ntrseg,
     .  npanu,e0,weight,nacell,ifpath,iupdte,
     .  ityp,nblock,masurf,msurf,mtsurf,
     .  nrcell,npcell,nblock,nntcll,ntcell,
     .  time

      

      RETURN
      END
C ===== SOURCE: vdion.f


      FUNCTION VDION (I)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP) :: VDION
      VDION=0.
      RETURN
      END
C ===== SOURCE: vecusr.f


      SUBROUTINE VECUSR (I,VX,VY,VZ,IPLS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      RETURN
      END
C ===== SOURCE: volusr.f


      SUBROUTINE VOLUSR(N,A)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: A(*)
      INTEGER, INTENT(IN) :: N
      RETURN
      END
C ===== SOURCE: uptcop.f
C
C
      SUBROUTINE UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CSPEZ
      USE CZT1
      USE CGEOM
      USE COMPRT
      USE CSDVI
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                      XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: P, WTRSIG, EION, V0_PARB, PARMOM_0, DIST, WTR
      INTEGER :: IAEL, IREL, IPL2, IAEI, IRDS, IBGK, IICX, IIEI, IIEL,
     .           IMEL, IPL1, I, IPL, IIO, IRD, IP, IR, IML, IAT, IFIRST,
     .           IRCX, IADD, ICOU, IACX, IRDD, IMCX, IMEI, IPLSTI,
     .           IPLSV, IPLV
      INTEGER, SAVE :: NMTSP
      REAL(DP), ALLOCATABLE, SAVE ::
     . CNDYNA(:), CNDYNM(:), CNDYNI(:)
CDR
      REAL(DP), ALLOCATABLE, SAVE ::
     . VPX(:),    VPY(:),    VRX(:),    VRY(:)
CDR
      DATA IFIRST/0/
!pb      SAVE
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNA(NATM))
        ALLOCATE (CNDYNM(NMOL))
        ALLOCATE (CNDYNI(NION))
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        DO 11 IAT=1,NATMI
11        CNDYNA(IAT)=AMUA*RMASSA(IAT)
        DO 12 IML=1,NMOLI
12        CNDYNM(IML)=AMUA*RMASSM(IML)
        DO 13 IIO=1,NIONI
13        CNDYNI(IIO)=AMUA*RMASSI(IIO)
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
C
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        DO 1 I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
1       CONTINUE
        DO 2 IR=1,NR1STM
          DO 2 IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
2       CONTINUE
C
        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI
C
      ENDIF
C
C  WV=WEIGHT/VEL
C
C  ATOMS
      IF (ITYP.EQ.1) THEN
        DO 20 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 20
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C  1,NPLSI:
C              PARTICLE CHARGE EXCHANGE RATE DUE TO IPLS: #/S
C              WITH ATOM SPECIES IATM=1,NATMI, PER ION
C  EACH RATE IS WEIGHTED WITH THE FACTOR (E0/EI-1), E0 BEING
C  THE NEUTRAL PARTICLE ENERGY, EI THE MEAN PLASMA ION ENERGY
C  THESE RATES ARE SCALED IN THE SHORT CYCLE WITH EI*NI
C
C
          IF (NCPVI.LT.NPLSI) GOTO 20
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 51
          DO 52 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 52
            IPLSTI = MPLSTI(IPLS)
            EION=1.5*TIIN(IPLSTI,IRD)+EDRIFT(IPLS,IRD)
            WTRSIG=WTR*SIGVCX(IRCX)/DIIN(IPLS,IRD)
            COPV(IPLS,IRDD)=COPV(IPLS,IRDD)+WTRSIG*(E0/EION-1.)
            LMETSP(NMTSP+IPLS)=.TRUE.
52        CONTINUE
51        CONTINUE
C
C.........................................
C
C   MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C.........................................
C
C
C  CONTRIBUTIONS FROM ATOMS
C  NPLSI+1, 2*NPLSI:
C
          IF (NCPVI.LT.2*NPLSI) GOTO 20
C
          IADD=NPLSI
          V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD))
          PARMOM_0=V0_PARB*CNDYNA(IATM)
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 59
          DO 56 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 56
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 56
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHA
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM(IPL,IRD)
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                        SIGN(1._DP,BVIN(IPLV,IRD))
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
56        CONTINUE
59        CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 61 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 62 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IRDS)*P
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  IPLV=MPLSV(IPL)
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                            SIGN(1._DP,BVIN(IPLV,IRD))
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
62            CONTINUE
            ENDIF
61        CONTINUE
C
C  ION IMPACT IONIZATION CONTRIBUTION: NOT INCLUDED
C
C
C  ELASTIC CONTRIBUTION FROM ATOMS
C
C
          IF (LGAEL(IATM,0,0).EQ.0) GOTO 80
C  DEFAULT TRACKLENGTH ESTIMATOR (BGK APPROXIMATION)
          DO 81  IAEL=1,NAELI(IATM)
            IREL=LGAEL(IATM,IAEL,0)
            IPLS=LGAEL(IATM,IAEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 81
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 81
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB*
     .                      SIGN(1._DP,BVIN(IPLSV,IRD))
            LMETSP(NMTSP+IPL2)=.TRUE.
81        CONTINUE
80      CONTINUE
C
20      CONTINUE
C
C  MOLECULES
      ELSEIF (ITYP.EQ.2) THEN
C
        DO 200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 200
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C
C  CONTRIBUTIONS FROM MOLECULES
C  2*NPLSI+1, 3*NPLSI:
C
          IF (NCPVI.LT.3*NPLSI) GOTO 200
C
          IADD=2*NPLSI
          V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD))
          PARMOM_0=V0_PARB*CNDYNM(IMOL)
C
          IF (LGMCX(IMOL,0,0).EQ.0) GOTO 590
          DO 560 IMCX=1,NMCXI(IMOL)
            IRCX=LGMCX(IMOL,IMCX,0)
            IPLS=LGMCX(IMOL,IMCX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 560
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 560
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHM
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM(IPL,IRD)
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                        SIGN(1._DP,BVIN(IPLV,IRD))
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
560       CONTINUE
590       CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 610 IMEI=1,NMDSI(IMOL)
            IRDS=LGMEI(IMOL,IMEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 620 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IRDS)*P
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  IPLV=MPLSV(IPL)
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                            SIGN(1._DP,BVIN(IPLV,IRD))
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
620           CONTINUE
            ENDIF
610       CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM MOLECULES
C
C
          IF (LGMEL(IMOL,0,0).EQ.0) GOTO 800
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 810 IMEL=1,NMELI(IMOL)
            IREL=LGMEL(IMOL,IMEL,0)
            IPLS=LGMEL(IMOL,IMEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 810
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 810
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLSV,IRD))
            LMETSP(NMTSP+IPL2)=.TRUE.
810       CONTINUE
800     CONTINUE
C
C
200     CONTINUE
C
C  TEST IONS
C
      ELSEIF (ITYP.EQ.3) THEN
C
        DO 2000 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 2000
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C
C  CONTRIBUTIONS FROM TEST IONS
C  3*NPLSI+1, 4*NPLSI:
C
          IF (NCPVI.LT.4*NPLSI) GOTO 2000
C
          IADD=3*NPLSI
          V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD))
          PARMOM_0=V0_PARB*CNDYNI(IION)
C
          IF (LGICX(IION,0,0).EQ.0) GOTO 5900
          DO 5600 IICX=1,NICXI(IION)
            IRCX=LGICX(IION,IICX,0)
            IPLS=LGICX(IION,IICX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 5600
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 5600
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHI
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
C
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM(IPL,IRD)
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              IPL2=IADD+IPL
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                        SIGN(1._DP,BVIN(IPLV,IRD))
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
5600      CONTINUE
5900      CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 6100 IIEI=1,NIDSI(IION)
            IRDS=LGIEI(IION,IIEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 6200 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IRDS)*P
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  IPLV=MPLSV(IPL)
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                            SIGN(1._DP,BVIN(IPLV,IRD))
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
6200          CONTINUE
            ENDIF
6100      CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM TEST IONS
C
          IF (LGIEL(IION,0,0).EQ.0) GOTO 8000
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 8100 IIEL=1,NIELI(IION)
            IREL=LGIEL(IION,IIEL,0)
            IPLS=LGIEL(IION,IIEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 8100
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 8100
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*PARMOM(IPLS,IRD)
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLSV,IRD))
            LMETSP(NMTSP+IPL2)=.TRUE.
8100      CONTINUE
8000    CONTINUE
C
C
2000    CONTINUE
C
C
      ENDIF
C
      RETURN
      END
