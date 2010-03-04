c     -*-Fortran-*-
c ======================================================================
c
      SUBROUTINE CalcMetricQuickandDirty
      IMPLICIT none

      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'slcom'

      REAL    CALCTHETA
      REAL    FACTOR,F1,F2

      INTEGER IKR,IRR,IKTMP,errcnt,ir,ik,id,i,ik1
      LOGICAL cont
      REAL thesft,themin,extthe,shift,delta,len1,len2,theta2(MAXNKS)



      WRITE(0,*) 'Q&D: CALCULATING METRIC'

c...  Determine the orthogonality of all cells that are complete polygons:
      CALL CalcOrth

      CALL SetupGrid

      CALL BuildMap

c...  Some out of place initialization:
      DO IR=1,NRS
         DO IK=1,NKS(IR)
            TAGDV(IK,IR) = 1
         ENDDO
      ENDDO

      thetag = -99.0

      DO ik = 1, nks(irsep)
        thetag(ik,irsep) = kps(ik,irsep) * 10.0
      ENDDO

      dthetg = thetag(ikti-1,irsep) - thetag(ikto+1,irsep) 

c...  SOL:
      DO ir = irsep + 1, irwall-1
        IF (ringtype(ir).EQ.SOL1) THEN
          DO ik = 1, nks(ir)
            thetag(ik,ir) = thetag(ikins(ik,ir),irins(ik,ir))
          ENDDO
c...      Finish off rings:
          DO ik = 2, nks(ir)
            IF (thetag(ik,ir).EQ.-99.0.AND.thetag(ik-1,ir).NE.-99.0) 
     .        thetag(ik,ir) = thetag(ik-1,ir) + 
     .                        (kps(ik,ir) - kps(ik-1,ir)) * 10.0
          ENDDO
          DO ik = nks(ir)-1, 1, -1
            IF (thetag(ik,ir).EQ.-99.0.AND.thetag(ik+1,ir).NE.-99.0) 
     .        thetag(ik,ir) = thetag(ik+1,ir) -
     .                        (kps(ik+1,ir) - kps(ik,ir)) * 10.0
          ENDDO
        ENDIF
      ENDDO


c...  Core:
      DO ir = irsep-1, 2, -1
        DO ik = 1, nks(ir)
          thetag(ik,ir) = thetag(ikouts(ik,ir),irouts(ik,ir))
        ENDDO
      ENDDO

c...  Primary PFZ:
      DO ir = nrs, irtrap+1, -1
        DO ik = 1, nks(ir)
          thetag(ik,ir) = thetag(ikouts(ik,ir),irouts(ik,ir))

          IF (ir.EQ.nrs.AND.ik.GE.ikti2(ir)) 
     .      thetag(ik,ir) = thetag(ik,ir) - dthetg
        ENDDO
      ENDDO

c...  Secondary PFZ, if there is one:
      cont = .TRUE.
      DO WHILE (cont) 
        cont = .FALSE.
        DO ir = irsep, irwall-1
c          WRITE(0,*) 'SEARCHING FOR 2nd PFZ'
          DO ik = 1, nks(ir)
            IF (thetag(ik           ,ir           ).EQ.-99.0.AND.
     .          thetag(ikouts(ik,ir),irouts(ik,ir)).NE.-99.0) THEN
c              WRITE(0,*) '  WORKING HARD FOR ',ik,ir
              thetag(ik,ir) = thetag(ikouts(ik,ir),irouts(ik,ir))
              cont = .TRUE.
            ENDIF
          ENDDO
        ENDDO  
      ENDDO

c...  Clean up:
      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        theta2(1:nks(ir)) = thetag(1:nks(ir),ir)
        DO ik = 2, nks(ir)-1
          IF (theta2(ik).EQ.theta2(ik+1)) THEN
            thetag(ik  ,ir) = -99.0 
            thetag(ik+1,ir) = -99.0 
          ENDIF
          IF (theta2(ik).EQ.theta2(ik-1)) THEN
            thetag(ik  ,ir) = -99.0 
            thetag(ik-1,ir) = -99.0 
          ENDIF
        ENDDO 
c...    Finish off rings:
        DO ik = 2, nks(ir)
          IF (thetag(ik,ir).EQ.-99.0.AND.thetag(ik-1,ir).NE.-99.0) 
     .      thetag(ik,ir) = thetag(ik-1,ir) + 
     .                      (kps(ik,ir) - kps(ik-1,ir)) * 10.0
        ENDDO
        DO ik = nks(ir)-1, 1, -1
          IF (thetag(ik,ir).EQ.-99.0.AND.thetag(ik+1,ir).NE.-99.0) 
     .      thetag(ik,ir) = thetag(ik+1,ir) -
     .                      (kps(ik+1,ir) - kps(ik,ir)) * 10.0
        ENDDO
      ENDDO



      CALL RepairMetric




c
c     Boundary rings:
c
      DO ik = 1, nks(1)
        thetag(ik,1) = thetag(ikouts(ik,1),irouts(ik,1))
      ENDDO

      DO ik = 1, nks(irwall)
        thetag(ik,irwall) = thetag(ikins(ik,irwall),irins(ik,irwall))
      ENDDO

      DO ik = 1, nks(irtrap)
        thetag(ik,irtrap) = thetag(ikouts(ik,irtrap),irouts(ik,irtrap))
      ENDDO


      IF (CTARGOPT.EQ.6) THEN
        DO id = 1, nds
          ik = ikds(id)
          ir = irds(id)

          IF (ik.EQ.1) THEN
            ik1 = ik + 1
          ELSE
            ik1 = ik - 1
          ENDIF

c          WRITE(0,*) id,ik,nks(ir),ir

          delta = SQRT((rs(ik,ir) - rs(ik1,ir))**2 +
     .                 (zs(ik,ir) - zs(ik1,ir))**2)

          shift = SQRT((rp(id) - rs(ik,ir))**2 +
     .                 (zp(id) - zs(ik,ir))**2) + delta


          thetat(id) = thetag(ik1,ir) +
     .      (thetag(ik,ir) - thetag(ik1,ir)) *
     .      (shift / delta)

c          WRITE(0,*) id,ik,ik1,ir,thetag(ik,ir),thetag(ik1,ir),
c     .               shift/delta,thetat(id)
        ENDDO
      ELSE
        STOP 'OBSOLETE TARGET OPTION'
      ENDIF



c...  Intialize the metric on boundary rings:
      CALL ResetRing(1     ,2)
      CALL ResetRing(irtrap,irtrap+1)
c...  BUG: Referencing IRWALL-1 is no good if IRWALL-1 is broken:
      CALL ResetRing(irwall,-1)


c
c     The THETA values are shifted so that the minimum THETA value on the
c     grid is 1.  This is important for some of the code in the DIV module
c     which assumes that THETA less than zero implies particle impact at the
c     target:
c
      themin =  1.0

      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          IF (thetag(ik,ir).LT.themin) themin = thetag(ik,ir)
        ENDDO
      ENDDO

      IF (ctargopt.EQ.1.OR.ctargopt.EQ.6) THEN
        DO I = 1, nds
          IF (thetat(i).LT.themin) themin = thetat(i)
        ENDDO
      ENDIF

      IF (themin.LT.1.0) THEN
        thesft = -themin + 1.0
      ELSE
        thesft = 0.0
      ENDIF

      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          thetag(ik,ir) = thetag(ik,ir) + thesft
        ENDDO
      ENDDO

      DO i = 1, nds
        thetat(i) = thetat(i) + thesft
      ENDDO


c      WRITE(0,*) 'THESHIFT=',thesft



c      CALL OutputData(85,'THETA GORE')
c      STOP 'sdfsdfsdl'

      RETURN
c
c     Error output:
99    CONTINUE
      STOP
      END
c
c ======================================================================
c
c subroutine: RepairMetric
c
      SUBROUTINE RepairMetric
      IMPLICIT none

      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'slcom'

      INTEGER ik,ir,ik1,ik3,count
      REAL    theta1,theta2,theta3,frac

      count = 0

      CALL OutputGrid(86,'Before metric repair')

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        ik = 1
        DO WHILE (ik.LT.nks(ir)-1)
          ik = ik + 1

          theta1 = thetag(ik-1,ir)
          theta2 = thetag(ik  ,ir)

          IF (ik.LT.nks(ir)-3) THEN         
            theta3 = 0.0
            DO ik1 = ik+1, nks(ir)-2
              IF (thetag(ik1  ,ir).LT.thetag(ik1+1,ir).AND.
     .            thetag(ik1+1,ir).LT.thetag(ik1+2,ir).AND.
     .            thetag(ik1  ,ir).GT.theta1.AND.
     .            theta3.EQ.0.0) THEN
                ik3 = ik1
                theta3 = thetag(ik3,ir)
              ENDIF
            ENDDO
          ELSE
            ik3 = ik + 1
            theta3 = thetag(ik3,ir)
          ENDIF


c          IF (ir.EQ.20) 
c            WRITE(0,'(A,3I6,3F12.4)') 'METRIC:',count,ik,ir,
c     .                 theta1,theta2,theta3

          IF    ((theta1.LT.theta2.AND.theta2.GT.theta3.AND.
     .            theta1.LT.theta3).OR.
     .           (theta1.GT.theta2.AND.theta2.LT.theta3.AND.
     .            theta1.LT.theta3)) THEN

            frac = (kps(ik ,ir) - kps(ik-1,ir)) /
     .             (kps(ik3,ir) - kps(ik-1,ir))

            thetag(ik,ir) = theta1 + frac * (theta3 - theta1)

            ik = 1
            count = count + 1

          ELSEIF (theta2.GE.theta3.AND.ik.EQ.nks(ir)-1) THEN
            frac = (kps(ik3,ir) - kps(ik-1,ir)) /
     .             (kps(ik ,ir) - kps(ik-1,ir))

            thetag(ik+1,ir) = thetag(ik-1,ir) +
     .                        frac * (thetag(ik,ir) - thetag(ik-1,ir))
            ik = 1
            count = count + 1

 

          ELSEIF (theta1.GT.theta2.AND.theta2.GE.theta3) THEN
            IF (ik.EQ.2) THEN
              CALL WN('RepairMetric','Starting point '//
     .                'for metric is poorly defined')

              thetag(ik,ir) = theta1 + 0.5 * ABS(theta1 - theta2)

            ELSE

              frac = (kps(ik  ,ir) - kps(ik-2,ir)) /
     .               (kps(ik-1,ir) - kps(ik-2,ir))
 
              thetag(ik,ir) = thetag(ik-2,ir) +
     .                        frac * (thetag(ik-1,ir) - thetag(ik-2,ir))
            ENDIF

            ik = 1
            count = count + 1

          ELSEIF (theta1.LT.theta2.AND.theta2.LT.theta3) THEN
          ELSE
            CALL ER('RepairMetric','Unrecognized metric problem '//
     .              'sequence',*98)
c     .              'sequence',*99)

          ENDIF

c          IF (count.EQ.10) STOP 'COUTING>>FSDFD'

        ENDDO

      ENDDO

c...  Check 1st and last points:
      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        IF (thetag(1,ir).GT.thetag(2,ir)) THEN
          count = count + 1

          thetag(1,ir) = 2.0 * thetag(2,ir) - 
     .                         thetag(3,ir)
        ENDIF

        IF (thetag(nks(ir),ir).LT.thetag(nks(ir)-1,ir)) THEN
          count = count + 1

          thetag(nks(ir),ir) = 2.0 * thetag(nks(ir)-1,ir) - 
     .                               thetag(nks(ir)-2,ir)
        ENDIF
       
      ENDDO

      IF (count.GT.0) THEN
        CALL WN('RepairMetric','Metric repairs attempted')
        WRITE(SLOUT,*) '  COUNT = ',count
        WRITE(SLOUT,*)
      ENDIF

      RETURN
98    IF (grdnmod.NE.0) THEN
        IF (sloutput) THEN
          WRITE(0,*)
          WRITE(0,*) '-------------------------------------------------'
          WRITE(0,*) '(NOT STOPPING THE CODE -- ASSUME 2 X-POINT ISSUE)'
          WRITE(0,*) '-------------------------------------------------'
          WRITE(0,*)
        ENDIF 
        RETURN
      ENDIF
99    CONTINUE
      CALL OutputData(85,'Theta error reported')
      WRITE(SLOUT,*) 'CURRENT INDEX:',ik,ik3
      WRITE(SLOUT,*) 'THETA1-3:',theta1,theta2,theta3
      WRITE(SLOUT,*) 'LOGIC:',theta1.LT.theta2,theta2.LT.theta3
      DO ik = 1, nks(ir)
        WRITE(SLOUT,'(A,2I6,F12.4)') 'THETAG:',ik,ir,thetag(ik,ir)
      ENDDO
      CALL DumpGrid('PROBLEMS CALCULATING THETAG')
      END
c
c ======================================================================
c
c subroutine: CalcMetric
c
c THETA Generation Code
c
c Variables:
c
c THECRE - line parameter for reference ring
c THECCU - line parameter for current ring
c THEA   - orthogonal reference r co-ordinate
c THEB   - orthogonal reference z co-ordinate
c THESFT - the shift in THEAT as a result of non-orthogonality
c THEX   - displacement between centerpoints of IK and IK-1 on IR
c IKR - reference cell for THETA
c IRR - reference ring for THETA
c THETOL - non-orthogonality tolerance for choosing a separatrix
c          reference knot
c THESCL - THETA scale for extrapolating THETA into triangular cells
c
c ======================================================================
c
      SUBROUTINE CalcMetric

      IMPLICIT none

      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'slcom'

      REAL    CALCTHETA
      REAL    FACTOR,F1,F2

      INTEGER IKR,IRR,IKTMP,errcnt,ir,ik,id,i,ik1,ir1,ir2
      REAL thesft,themin,extthe,shift,delta,len1,len2



c      WRITE(0,*) 'CAlculating metric'


      IF (stopopt.EQ.14) THEN
        WRITE(0,*) 
        WRITE(0,*) '*** NOT CALCULATING METRIC ***'
        WRITE(0,*) 
        RETURN
      ENDIF
c
c Determine the orthogonality of all cells that are complete polygons:
c
      CALL CalcOrth
      CALL SetupGrid
c
c Some out of place initialization:
c
         DO IR=1,NRS
            DO IK=1,NKS(IR)
               TAGDV(IK,IR) = 1
            ENDDO
         ENDDO
c
c     Separatix:
c
c      DO ik = 1, nks(irsep)
c        thetag(ik,irsep) = REAL(ik)
c      ENDDO

      CALL MS('CalcMetric','Using KPS to initialize separatrix metric')

      DO ik = 1, nks(irsep)
        thetag(ik,irsep) = kps(ik,irsep) * 100.0
      ENDDO


      dthetg = thetag(ikti - 1,irsep) -
     +         thetag(ikto + 1,irsep) + 1.0

      CALL OutputGrid(85,'Before calculating metric')
c
c     Scrape off layer:
c     ------------------------------------------------------------------
      DO ir = irsep+1, irwall-1
        DO ik = 1, nks(ir)
          thetag(ik,ir) = CalcTheta(ik,ir)
        ENDDO
      ENDDO
c
c     Secondary PFZ for double-null grids:
c     ------------------------------------------------------------------
      DO ir = irsep+1, irwall-1
        IF (ringtype(ir).EQ.PFZ) EXIT
      ENDDO
      IF (ir.NE.irwall) THEN
c...    Identify the rings just outside the ring that's just inside the 
c       secondary separatrix
        ir1 = irouts(1          ,irsep2)
        ir2 = irouts(nks(irsep2),irsep2)  
c...    Extend THETAG on these rings as per the definition on the
c       primary separatrix ring:
        ir = ir1
        DO ik = 1, nks(ir)
          IF (irins(ik,ir).NE.irsep2) THEN
            thetag(ik,ir) =
     .        thetag(ik-1,ir) + (kps(ik  ,ir) - kps(ik-1,ir)) * 100.0
          ENDIF
        ENDDO
        ir = ir2
        DO ik = nks(ir), 1, -1
          IF (irins(ik,ir).NE.irsep2) THEN
            thetag(ik,ir) = 
     .        thetag(ik+1,ir) - (kps(ik+1,ir) - kps(ik  ,ir)) * 100.0
          ENDIF
        ENDDO
        write(0,*) 'IR1,2=',ir1,ir2

        DO ir = ir1+1, irwall-1
          IF (ir.EQ.ir2.OR.ringtype(ir).EQ.PFZ) CYCLE
          DO ik = 1, nks(ir)
            thetag(ik,ir) = CalcTheta(ik,ir)
          ENDDO        
        ENDDO

        DO ir = irwall-1, ir1+1, -1
          IF (ir.EQ.ir2.OR.ringtype(ir).NE.PFZ) CYCLE
          DO ik = 1, nks(ir)
            thetag(ik,ir) = CalcTheta(ik,ir)
          ENDDO        
        ENDDO

      ENDIF

c
c     Private plasma:
c     ------------------------------------------------------------------
      DO ir = nrs, irtrap+1, -1
        DO ik = 1, nks(ir)
          thetag(ik,ir) = CalcTheta(ik,ir)
        ENDDO

c...    Seach for lost cells:
        DO ik = ikmids(ir), 1, -1
c          IF (ir.EQ.38) WRITE(0,*) '--->',ik,ir,thetag(ik,ir)
          IF (thetag(ik,ir).EQ.1.0) THEN
c...        Assume that the theta value can be extrapolated:
            thetag(ik,ir) = thetag(ik+2,ir) + 
     .                      (thetag(ik+1,ir) - thetag(ik+2,ir)) * 
     .                      (kps(ik  ,ir) - kps(ik+2,ir)) / 
     .                      (kps(ik+1,ir) - kps(ik+2,ir))
          ENDIF
        ENDDO
        DO ik = ikmids(ir)+1, nks(ir)
c          IF (ir.EQ.38) WRITE(0,*) '--->',ik,ir,thetag(ik,ir)
          IF (thetag(ik,ir).EQ.1.0) THEN
c...        Assume that the theta value can be extrapolated:
            thetag(ik,ir) = thetag(ik-2,ir) + 
     .                      (thetag(ik-1,ir) - thetag(ik-2,ir)) * 
     .                      (kps(ik  ,ir) - kps(ik-2,ir)) / 
     .                      (kps(ik-1,ir) - kps(ik-2,ir))
          ENDIF
        ENDDO
      ENDDO

c...  There is sometimes a problem with the near-x-point cells of the
c     PFZ do to the approximate nature of DTHETG -- just smooth the metric:
      IF (.NOT.nopriv) THEN
        IF (thetag(ikto2(nrs),nrs).GT.thetag(ikti2(nrs),nrs)) THEN
          ir = nrs
          len1 = kps(ikti2(ir)  ,ir) - kps(ikto2(ir),ir)
          len2 = kps(ikti2(ir)+1,ir) - kps(ikto2(ir),ir)
          thetag(ikti2(ir),ir) = thetag(ikto2(ir),ir) +
     .      len1 / len2 * (thetag(ikti2(ir)+1,ir) - 
     .                     thetag(ikto2(ir)  ,ir))
        ENDIF
      ENDIF
c
c     Core:
c
      DO ir = irsep-1, 2, -1
        DO ik = 1, nks(ir)-1
          thetag(ik,ir) = CalcTheta(ik,ir)
        ENDDO

        thetag(nks(ir),ir) = thetag(nks(ir)-1,ir) + 1.0
      ENDDO


      CALL RepairMetric

c
c     Boundary rings:
c
      DO ik = 1, nks(1)
        thetag(ik,1) = thetag(ikouts(ik,1),irouts(ik,1))
      ENDDO

      DO ik = 1, nks(irwall)
        thetag(ik,irwall) = thetag(ikins(ik,irwall),irins(ik,irwall))
      ENDDO

      DO ik = 1, nks(irtrap)
        thetag(ik,irtrap) = thetag(ikouts(ik,irtrap),irouts(ik,irtrap))
      ENDDO
c
c ----------------------------------------------------------------------
c Initialize THETAT for CTARGOPT 1&6 - other options are in WALLS.D4A:
c ----------------------------------------------------------------------
c
c slnote 06/07/95
c
c As far as I can tell, THETAT is only defined for CTARGOPT.EQ.1,
c but it is used in DIV depending on NORTHOPT, not CTARGOPT.  This
c was why I was running into some trouble with cross field transport near
c the target. Maybe a check should be added to identify illegal
c CTARGOPT and NORTHOPT combinations.
c
c Is THETAT defined in the core?
c
      IF (CTARGOPT.EQ.1) THEN
        CALL WN('CalcMetric','THETAT questionable for target option 1')

         id = 0
         do ir = irwall, irsep, -1
            id = id + 1
            THETAT(ID)=THETAG(NKS(IR),IR)
         enddo
         do ir = nrs, irtrap, -1
            id = id + 1
            THETAT(ID)=THETAG(NKS(IR),IR)
         enddo
         do ir = irtrap, nrs
            id = id + 1
            THETAT(ID)=THETAG(1,IR)
         enddo
         do ir = irsep, irwall
            id = id + 1
            THETAT(ID)=THETAG(1,IR)
         enddo

      ELSEIF (CTARGOPT.EQ.6) THEN
c
        DO id = 1, nds
          ik = ikds(id)
          ir = irds(id)

          IF (ik.EQ.1) THEN
            ik1 = ik + 1
          ELSE
            ik1 = ik - 1
          ENDIF

          delta = SQRT((rs(ik,ir) - rs(ik1,ir))**2 +
     .                 (zs(ik,ir) - zs(ik1,ir))**2)

          shift = SQRT((rp(id) - rs(ik,ir))**2 +
     .                 (zp(id) - zs(ik,ir))**2) + delta

c          WRITE(0,*) id,ik,ik1,ir,thetag(ik,ir),thetag(ik1,ir),
c     .               shift/delta

          thetat(id) = thetag(ik1,ir) +
     .      (thetag(ik,ir) - thetag(ik1,ir)) *
     .      (shift / delta)
        ENDDO

      ENDIF
c
c     Intialize the metric on boundary rings:
c
      CALL ResetRing(1     ,2)
      CALL ResetRing(irtrap,irtrap+1)
c...  BUG: Referencing IRWALL-1 is no good if IRWALL-1 is broken:
      CALL ResetRing(irwall,-1)
c      CALL ResetRing(irwall,irwall-1)


c
c     The THETA values are shifted so that the minimum THETA value on the
c     grid is 1.  This is important for some of the code in the DIV module
c     which assumes that THETA less than zero implies particle impact at the
c     target:
c
      themin =  1.0

      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          IF (thetag(ik,ir).LT.themin) themin = thetag(ik,ir)
        ENDDO
      ENDDO

      IF (ctargopt.EQ.1.OR.ctargopt.EQ.6) THEN
        DO I = 1, nds
          IF (thetat(i).LT.themin) themin = thetat(i)
        ENDDO
      ENDIF

      IF (themin.LT.1.0) THEN
        thesft = -themin + 1.0
      ELSE
        thesft = 0.0
      ENDIF

      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          thetag(ik,ir) = thetag(ik,ir) + thesft
        ENDDO
      ENDDO

      DO i = 1, nds
        thetat(i) = thetat(i) + thesft
      ENDDO

      CALL OutputGrid(87,'After calculating metric')

c      CALL OutputData(85,'THETA PATCHING')
c      CALL DumpGrid('THETA PATCHING')
      RETURN
c
c     Error output:
999   CONTINUE
      STOP
      END
c
c ======================================================================
c FUNCTION: CALCTHETA
c ======================================================================
c
      REAL FUNCTION CalcTheta(ikcell,ircell)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

c     Input:
      INTEGER ikcell,ircell

      INTEGER CalcPoint

      INTEGER   RANGE
c      PARAMETER (RANGE = 20)


      INTEGER  IK,IR,IKRE,IRRE,IKR,IRR,IDR,ISIDE,nV1,nV2,IK1,IK2,
     .         ika,ikb,id,ii,ii1,ikstart,ikend,ii2,ik3,ii3,irref,ii4,
     .         dum,ikack,ikmin,ikmax,id1,id3
      REAL     THECRE,THECCU,THEA,THEB,THESFT,THEX,THESCL,THETA
      REAL     DIST,DISTB,DISTF,ra,za,rb,zb,distmin,tmin,
     .  shift,shiftmin,r1,r2,z1,z2,t1,t2,r,z,t3,z3,r3,delta

      DOUBLE PRECISION ar,az,br,bz,cr,cz,t

      range = 10

      ir  = ircell
      ik1 = NULL

      IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) THEN
        CalcTheta = 0
        RETURN
      ENDIF

      IF (ir.EQ.irsep)
     .  CALL ER('CalcTheta','Separatrix must be pre-initialised',*99)
c
c
c
      IF     (ir.LT.irsep) THEN
        ikack = ikouts(ikcell,ir)
        irref = irouts(ikcell,ir)
        nv1   = 2
        nv2   = 3
        ikmin = 1
        ikmax = nks(ir) - 1
      ELSEIF (ir.GT.irtrap.OR.ringtype(ir).EQ.PFZ) THEN
        ikack = ikouts(ikcell,ir)
        irref = irouts(ikcell,ir)
        nv1   = 2
        nv2   = 3
        IF (ikcell.LE.ikto2(ir)) THEN
          ikmin = 1
          ikmax = ikto2(ir)
        ELSE
          ikmin = ikti2(ir)
          ikmax = nks(ir)
        ENDIF
c...BROKEN PFZ:
        IF (irref.EQ.irwall) THEN
          CalcTheta = 1.0
          RETURN
        ENDIF
      ELSE
        ikack = ikins(ikcell,ir)
        irref = irins(ikcell,ir)
        nv1   = 1
        nv2   = 4
        ikmin = 1
        ikmax = nks(ir)
      ENDIF

10    CONTINUE
      ikstart = MAX(ikmin,ikcell-RANGE)
      ikend   = MIN(ikmax,ikcell+RANGE)
c
c     Find the minimum perpendicular distance from the center of
c     cell IKCELL,IR to the polygon side of a cell in the
c     relevant cross-field direction.  As default, take the
c     perpendicular distance to the polygon side of IKCELL,IR:
c
      distmin = SQRT((rs(ikcell,ir) - rs(ikack,irref))**2.0 +
     .               (zs(ikcell,ir) - zs(ikack,irref))**2.0)
c Find a better scale distance...
c      distmin = 0.25 * r0
      id      = korpg(ikcell,ir)
      cr      = rs   (ikcell,ir)
      cz      = zs   (ikcell,ir)

      DO ik = ikstart, ikend
        id = korpg(ik,ir)
        ar = rvertp(nv1,id)
        az = zvertp(nv1,id)
        br = rvertp(nv2,id)
        bz = zvertp(nv2,id)

        IF (CalcPoint(ar,az,br,bz,cr,cz,t).EQ.2) THEN
          r    = ar + t * (br - ar)
          z    = az + t * (bz - az)
          dist = SQRT((r - cr)**2 + (z - cz)**2)

          IF (dist.LT.distmin) THEN
            distmin = dist
            ik1     = ik
            t1      = t
            r1      = r
            z1      = z
          ENDIF
        ENDIF
      ENDDO

      IF (ik1.EQ.NULL) THEN
        IF (range.LT.nks(ir)/2) THEN
c May not need to do this anymore since scale distance is improved...
          range = range + 5
          GOTO 10
        ELSE
c
c This sucks... another day...
c
          tmin    = HI
          distmin = SQRT((rs(ikcell,ir) - rs(ikack,irref))**2.0 +
     .                   (zs(ikcell,ir) - zs(ikack,irref))**2.0)

          cr = rs(ikcell,ir)
          cz = zs(ikcell,ir)

          DO ik = MAX(ikmin,ikcell-nks(ir)/2),
     .            MIN(ikmax,ikcell+nks(ir)/2)
            id = korpg(ik,ir)
            ar = rvertp(nv1,id)
            az = zvertp(nv1,id)
            br = rvertp(nv2,id)
            bz = zvertp(nv2,id)

            r  = 0.5 * (ar + br)
            z  = 0.5 * (az + bz)

            dist = SQRT((r - cr)**2.0 + (z - cz)**2.0)

c            WRITE(50,'(5X,I4,2G18.6,I6,2G18.6)')
c     .        ik,t,ABS(t-0.5),ik1,dist,distmin

            IF (dist.LT.distmin) THEN
              distmin = dist

              dum = CalcPoint(ar,az,br,bz,cr,cz,t)

              r1  = ar + t * (br - ar)
              z1  = az + t * (bz - az)
              ik1 = ik
              t1  = t
            ENDIF
          ENDDO

          IF (rel_count.LT.0) WRITE(SLOUT,'(A)') 'Despirate focus...'

          IF (ik1.EQ.NULL) THEN
            IF (ir.GE.irsep) THEN
              CALL MS('CalcTheta','Lost cell')
              CalcTheta = 1.0
              RETURN
c              CALL ER('CalcTheta','Lost cell',*99)
            ELSE
              CALL MS('CalcTheta','Core THETA error')
              CalcTheta = 1.0
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
c
c     For each cell on the reference ring (within +-RANGE) find the
c     minimum perpendicular distance from that cell to a polygon
c     side between the reference ring and the ring in question (IRCELL).
c     Polygon data from the ring IRCELL is used, although data
c     from either ring is valid:
c
20    CONTINUE
      ikstart = MAX(ikmin,ikcell-RANGE)
      ikend   = MIN(ikmax,ikcell+RANGE)

      shiftmin = HI
      ik3      = NULL

      DO ii = MAX(1,ikack-RANGE), MIN(nks(irref),ikack+RANGE)

        distmin = SQRT((rs(ikcell,ir) - rs(ikack,irref))**2.0 +
     .                 (zs(ikcell,ir) - zs(ikack,irref))**2.0)
c Find a better scale distance...
c        distmin = 0.25 * r0
        ik2     = NULL

        DO ik = ikstart, ikend
c
c         For each polygon side between IRREF and IRCELL in the
c         vicinity of IKCELL, find the minimum distance from a
c         valid perpendicular intersection (for cell centers on
c         IRREF) to R1,Z1:
c
          id = korpg (ik,ir)
          ar = rvertp(nv1,id)
          az = zvertp(nv1,id)
          br = rvertp(nv2,id)
          bz = zvertp(nv2,id)
          cr = rs(ii,irref)
          cz = zs(ii,irref)

          IF (CalcPoint(ar,az,br,bz,cr,cz,t).EQ.2) THEN
            r    = ar + t * (br - ar)
            z    = az + t * (bz - az)
            dist = SQRT((r - cr)**2 + (z - cz)**2)

            IF (dist.LT.distmin) THEN
              distmin = dist
              ik2     = ik
              t2      = t
              r2      = r
              z2      = z
            ENDIF
          ENDIF


c            WRITE(50,'(5X,2I4,G18.6,I6,2G18.6)')
c     .        ik,ii,t,ik2,dist,distmin

        ENDDO

        IF (ik2.GT.NULL) THEN
          shift = SQRT((r2 - r1)**2 + (z2 - z1)**2)

          IF (ir.EQ.2) THEN

c            WRITE(50,'(5X,2G18.6,2I6)')
c     .        shift,shiftmin,ii,ik2

          ENDIF

          IF (shift.LT.shiftmin) THEN
            shiftmin = shift
            ii3      = ii
            ik3      = ik2
            t3       = t2
            r3       = r2
            z3       = z2
          ENDIF
        ENDIF
      ENDDO
c
c
c
      IF (ik3.EQ.NULL) THEN
        IF (range.LT.nks(ir)/2) THEN
          range = range + 5
          GOTO 20
        ELSE
          id = korpg(ikcell,ir)
          ar = rvertp(nv1,id)
          az = zvertp(nv1,id)
          br = rvertp(nv2,id)
          bz = zvertp(nv2,id)
          cr = rs(ikack,irref)
          cz = zs(ikack,irref)

          dum = CalcPoint(ar,az,br,bz,cr,cz,t)

          ii3 = ikack
          r3  = ar + t * (br - ar)
          z3  = az + t * (bz - az)
          ik3 = ikcell
          t3  = t

          IF (rel_count.LT.0)
     .      WRITE(SLOUT,'(A)') 'Despirate reference...'
        ENDIF
      ENDIF

c
c     Determine where the target cell is with respect to the reference
c     cell:
c
      id1 = korpg(ik1,ir)
      id3 = korpg(ik3,ir)

      IF     (ik1.LT.ik3) THEN
        shift = SQRT((r1 - rvertp(nv2,id1))**2.0 +
     .               (z1 - zvertp(nv2,id1))**2.0) +
     .          SQRT((r3 - rvertp(nv1,id3))**2.0 +
     .               (z3 - zvertp(nv1,id3))**2.0)

c        WRITE(50,'(A,2I5,F10.6)') '1: ',ik1,ik3,shift
      ELSEIF (ik3.LT.ik1) THEN
        shift = SQRT((r3 - rvertp(nv2,id3))**2.0 +
     .               (z3 - zvertp(nv2,id3))**2.0) +
     .          SQRT((r1 - rvertp(nv1,id1))**2.0 +
     .               (z1 - zvertp(nv1,id1))**2.0)

c        WRITE(50,'(A,2I5,F10.6)') '2: ',ik1,ik3,shift
      ELSE
        shift = SQRT((r3 - r1)**2 + (z3 - z1)**2)

c        WRITE(50,'(A,2I5,F10.6)') '3: ',ik1,ik3,shift
      ENDIF

      DO ik = MIN(ik1+1,ik3+1), MAX(ik1-1,ik3-1)
        id = korpg(ik,ir)
        shift = shift + SQRT((rvertp(nv2,id) - rvertp(nv1,id))**2.0 +
     .                       (zvertp(nv2,id) - zvertp(nv1,id))**2.0)

c        WRITE(50,'(A,2I5,F10.6,I5)') '4: ',ik1,ik3,shift,ik
      ENDDO





c This is really messy.. clean it up...
      IF (ii3.EQ.1.AND.
     .    (ik1.LT.ik3.OR.(ik1.EQ.ik3.AND.t1.LT.t3))) THEN
c
c       Extrapolate reference metric beyond end of ring:
c
        delta = SQRT((rs(ii3,irref) - rs(ii3+1,irref))**2 +
     .               (zs(ii3,irref) - zs(ii3+1,irref))**2)

        CalcTheta = thetag(ii3+1,irref) + (shift + delta) / delta *
     .              (thetag(ii3,irref) - thetag(ii3+1,irref))

      ELSEIF (ii3.EQ.nks(irref).AND.
     .    (ik1.GT.ik3.OR.(ik1.EQ.ik3.AND.t1.GT.t3))) THEN
c
c       Extrapolate reference metric beyond end of ring:
c
        delta = SQRT((rs(ii3,irref) - rs(ii3-1,irref))**2 +
     .               (zs(ii3,irref) - zs(ii3-1,irref))**2)

        CalcTheta = thetag(ii3-1,irref) + (shift + delta) / delta *
     .              (thetag(ii3,irref) - thetag(ii3-1,irref))

      ELSE
c
c       Interpolate reference metric:
c
        IF (ik1.LT.ik3.OR.(ik1.EQ.ik3.AND.t1.LT.t3)) THEN
          ii4 = ii3 - 1
        ELSE
          ii4 = ii3 + 1
        ENDIF

        delta = SQRT((rs(ii3,irref) - rs(ii4,irref))**2 +
     .               (zs(ii3,irref) - zs(ii4,irref))**2)

        CalcTheta = thetag(ii3,irref) + shift / delta *
     .              (thetag(ii4,irref) - thetag(ii3,irref))

c        IF (ir.EQ.37.AND.ii3.EQ.64) THEN
c          WRITE(SLOUT,*) 'shift=',shift,delta,shift/delta
c          WRITE(SLOUT,*) 'shift ',ii3,ii4,irref
c          WRITE(SLOUT,*) 'shift ',thetag(ii3,irref),thetag(ii4,irref)
c        ENDIF

      ENDIF

      IF (ir.EQ.nrs.AND.ikcell.GE.ikti2(ir))
     .  CalcTheta = CalcTheta - dthetg

c     Debug:
      IF (rel_count.LT.0) THEN
        WRITE(SLOUT,'(2X,A,3I4,F8.3,2X,3I4,F8.3,2F10.5)')
     .    'm: IR,IKc,IK1,T1 IRe,II3,IK3,T3,SHIFT the = ',
     .    ir,ikcell,ik1,t1,irref,ii3,ik3,t3,shift,CalcTheta
      ENDIF

      RETURN
c
c     Error output:
99    CONTINUE
      WRITE(EROUT,'(5X,A,8I4,F10.4)')
     .  'IKCELL,IR,IRREF,IKACK,NV1,NV2,IKSTART,IKEND,SHIFT = ',
     .  ikcell,ir,irref,ikack,nv1,nv2,ikstart,ikend,shift
      WRITE(EROUT,'(5X,A,I4,1P,3E13.6)') 'IK1,T1,R1,Z1 = ',ik1,t1,r1,z1
      WRITE(EROUT,'(5X,A,I4,1P,3E13.6)') 'IK2,T2,R2,Z2 = ',ik2,t2,r2,z2
      WRITE(EROUT,'(5X,A,I4,1P,3E13.6)') 'II3,T3,R3,Z3 = ',ii3,t3,r3,z3
      STOP
      END
c
c FUNCTION END: CALCTHETA ===============================================
c
c ======================================================================
c FUNCTION: EXTTHE(IK1,IR1,IDIR)
c ======================================================================
c
      REAL FUNCTION EXTTHE(IK,IR,IDIR)

      INTEGER IK,IR,IDIR

C     INCLUDE "PARAMS"
      INCLUDE 'params'
C     INCLUDE "CGEOM"
      INCLUDE 'cgeom'
C     INCLUDE "COMTOR"
      INCLUDE 'comtor'
      INCLUDE 'slcom'



      REAL    RT,ZT,L1,L2
      REAL    CENLEN
      INTEGER IN

      IF (CTARGOPT.NE.6) THEN
         CALL WARN('ExtThe','Target option not set to 6 as required',
     +             0,0,0)
         STOP
      ENDIF
c
c Extapolate backward:
c
      IF (IDIR.EQ.0) THEN

         IF (THETAG(IK+1,IR).EQ.-1) THEN
            WRITE( 0,*) 'Error   (EXTTHE): Bad cell reference:',IK+1,IR
c            WRITE(SLOUT,*)
c     .        'Error   (EXTTHE): Bad cell reference:',IK+1,IR
            STOP
         ENDIF

         IN = KORPG(IK,IR)

         RT = 0.5 * (RVERTP(1,IN) + RVERTP(2,IN))
         ZT = 0.5 * (ZVERTP(1,IN) + ZVERTP(2,IN))

         L1 = SQRT((RS(IK,IR) - RT)**2 + (ZS(IK,IR) - ZT)**2)
         L2 = CENLEN(IK,IR,IK+1,IR)

         EXTTHE = THETAG(IK,IR)-(THETAG(IK+1,IR)-THETAG(IK,IR))*L1/L2
c
c Extrapolate forward:
c
      ELSEIF (IDIR.EQ.1) THEN

         IF (THETAG(IK-1,IR).EQ.-1) THEN
            WRITE( 0,*) 'Error   (EXTTHE): Bad cell reference:',IK-1,IR
c            WRITE(SLOUT,*)
c     .        'Error   (EXTTHE): Bad cell reference:',IK-1,IR
            STOP
         ENDIF

         IN = KORPG(IK,IR)

         RT = 0.5 * (RVERTP(3,IN) + RVERTP(4,IN))
         ZT = 0.5 * (ZVERTP(3,IN) + ZVERTP(4,IN))

         L1 = SQRT((RS(IK,IR)-RT)**2 + (ZS(IK,IR)-ZT)**2)
         L2 = CENLEN(IK,IR,IK-1,IR)

         EXTTHE = THETAG(IK,IR)+(THETAG(IK,IR)-THETAG(IK-1,IR))*L1/L2

      ELSE
         WRITE( 0,*) 'Error   (EXTTHE): Illegal IDIR value.'
c         WRITE(SLOUT,*) 'Error   (EXTTHE): Illegal IDIR value.'
         STOP
      ENDIF

c      WRITE(SLOUT,'(A,2I4,4F10.5)') 'Diag: ',IK,IR,L1,L2,EXTTHE,
c     +                            THETAG(IK,IR)

      RETURN
      END
c
c FUNCTION END: EXTTHE =================================================
c

c
c ======================================================================
c SUBROUTINE: WARN
c ======================================================================
c
      SUBROUTINE WARN(MODULE,MESSAGE,PAR1,PAR2,PAR3)

      CHARACTER MODULE*(*),MESSAGE*(*)
      INTEGER   PAR1,PAR2,PAR3

      WRITE(51,'(4A,::)') ' Warning (',MODULE,'): ',MESSAGE
      IF (par3.NE.-1) WRITE( 0,'(4A,::)')
     +  ' Warning (',MODULE,'): ',MESSAGE

      IF (PAR1.NE.0) THEN
        WRITE(51,'(I4,1X,::)') PAR1
        IF (par3.NE.-1) WRITE( 0,'(I4,1X,::)') PAR1
      ENDIF

      IF (PAR2.NE.0) THEN
        WRITE(51,'(I4,1X,::)') PAR2
        IF (par3.NE.-1) WRITE( 0,'(I4,1X,::)') PAR2
      ENDIF

      IF (PAR3.NE.0) THEN
        WRITE(51,'(I4,1X,::)') PAR3
        IF (par3.NE.-1) WRITE( 0,'(I4,1X,::)') PAR3
      ENDIF

      WRITE(51,*) ' '
      IF (par3.NE.-1) WRITE( 0,*) ' '

      RETURN
      END
c
c SUBROUTINE END: WARN ================================================
c


c
c ======================================================================
c
c subroutine GetPoint
c
c
c Don't think I'll need this afterall.
c
c ======================================================================
c
      SUBROUTINE GetPoint(ik,ir,p,r,z)
c
c Input:
c
      REAL    p
      INTEGER ik,ir
c
c Output:
c
      REAL    r,z

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER id
      REAL    r1,r2,z1,z2,p1,p2,t
c
c Check that p is not out of bounds:
c



      IF (p.EQ.kps(ik,ir)) THEN
c
c Return value for cell center:
c
        r = rs(ik,ir)
        z = zs(ik,ir)
      ELSE
c
c Otherwise:
c
        id = korpg(ik,ir)

        IF (p.LE.kpb(ik-1,ir)) THEN
          p1 = kpb(ik,ir)
          p2 = kps(ik,ir)

          r1 = 0.5 * (rvertp(3,id) + rvertp(4,id))
          z1 = 0.5 * (zvertp(3,id) + zvertp(4,id))

          r2 = rs(ik,ir)
          z2 = zs(ik,ir)
         ELSE
          p1 = kps(ik  ,ir)
          p2 = kpb(ik-1,ir)

          r1 = rs(ik,ir)
          z1 = zs(ik,ir)

          r2 = 0.5 * (rvertp(1,id) + rvertp(2,id))
          z2 = 0.5 * (zvertp(1,id) + zvertp(2,id))
        ENDIF

        t = (p - p2) / (p1 - p2)
c
c Check to confirm that 0 < t < 1:
c


        r = r2 + t * (r1 - r2)
        z = z2 + t * (z1 - z2)
      ENDIF

c      WRITE(SLOUT,'(A,2I3,10E11.5)')
c     .  'GetPoint  : ',
c     .  ik,ir,
c     .  r,z,p,t,p1,p2,r1,r2,z1,z2

      RETURN
      END
c
c ======================================================================
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c ======================================================================
c
c ======================================================================
c
c ======================================================================
c
c ======================================================================
c
c
c Old THETAG calculation code.  To be deleted.
c
c
c ======================================================================
c
c SUBROUTINE: CALCTHETA
c
c THETA Generation Code
c
c Variables:
c
c THECRE - line parameter for reference ring
c THECCU - line parameter for current ring
c THEA   - orthogonal reference r co-ordinate
c THEB   - orthogonal reference z co-ordinate
c THESFT - the shift in THEAT as a result of non-orthogonality
c THEX   - displacement between centerpoints of IK and IK-1 on IR
c IKR - reference cell for THETA
c IRR - reference ring for THETA
c THETOL - non-orthogonality tolerance for choosing a separatrix
c          reference knot
c THESCL - THETA scale for extrapolating THETA into triangular cells
c
c ======================================================================
c
      SUBROUTINE CALCTHETA_Old
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cioniz'
      include 'reader'
      include 'dynam5'

      REAL    CALTHETA_Old
      REAL    FACTOR,F1,F2

      INTEGER IKR,IRR,IKTMP
c
      integer ik,ir,id,i
      real extthe_Old,themin,thesft
c
c Determine the orthogonality of all cells that are complete polygons:
c
      CALL CALCORTH
c
c Some out of place initialization:
c
         DO IR=1,NRS
            DO IK=1,NKS(IR)
               TAGDV(IK,IR) = 1
            ENDDO
         ENDDO
c
c Initialize THETAG along the separatix:
c
      DO IK = 1, NKS(IRSEP)
        THETAG(IK,IRSEP) = REAL(IK)
      ENDDO
c
c ----------------------------------------------------------------------
c For cells in the Scrape Off Layer:
c ----------------------------------------------------------------------
c
      DO IR = IRSEP+1, IRWALL

        FACTOR = 360.0

        IRR = IRINS(INT(0.5*NKS(IR)),IR)

        DO IK = 1, NKS(IR)

          F1 = ABS(ALPH(IK          ,IR )-90.0)
          F2 = ABS(ALPH(IKINS(IK,IR),IRR)-90.0)

c          IF (F1.LT.1.0.AND.F2.LT.1.0) THEN
c            IKR = IKOUTS(IK,IR)
c            GOTO 10
c          ELSE
            IF ((F1+F2).LT.FACTOR) THEN
              IKR    = IK
              FACTOR = F1 + F2
            ENDIF
c          ENDIF

        ENDDO

        IF (FACTOR.GT.2.0)
     +    CALL WARN('Calctheta','Bad SOL THETA reference for ring',
     +              IR,0,0)

10      CONTINUE

        WRITE(50,'(A,4I4,2F8.2)') 'SOL  THETA reference:',
     +        IKR,NKS(IR),IR,IRR,ALPH(IKR,IR),ALPH(IKINS(IKR,IR),IRR)

        THETAG(IKR,IR) = THETAG(IKINS(IKR,IR),IRR)

        DO IK = IKR-1, 1, -1
          THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKINS(IK,IR),IRR,0)
        ENDDO

        DO IK = IKR+1, NKS(IR)
          THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKINS(IK,IR),IRR,0)
        ENDDO
      ENDDO
c
c ----------------------------------------------------------------------
c For cells in the Private Flux Zone:
c ----------------------------------------------------------------------
c
c slnote 2/6/95
c
c There is the same assumption made here as that in
c the TRUE NEAREST NEIGHBOURS code for the calculation of DTHETG,
c namely that the grid is orthogonal in the x-point region.
c
c slnote 31/05/95
c
c This assumption of orthogonality if ALPH-90 +- 2 degrees
c is probably excessive and may be leading to the not-quite
c orthogonal motion in the least orthogonal regions of the
c private plasma. However, the large 'tolerance' is needed
c in order for the THETA generation routine to behave well. - FIXED,
c now set to 1 degree.
c
c slnote
c
c There should also be a check here for the case where there
c is no reference knot but where THEA and THEB are not zero.
c I am not sure that such a case could exist with the current
c C-Mod grids (May '95) but it is something to keep in mind.
c
c It is assumed that there are no triangular cells that are not
c near the ends of the rings.  There is also the assuption that
c the cells along the 'symmetry' point in the private plasma
c are orthogonal.
c
c ASSUMPTION - the grid is orthogonal in the x-point region when
c              calculating DTHETG
c
c
c NOTE (Aug 8, 96) - Check that IKTO is okay for JET CMOD grids
c

c
c NOTE: Check that there is a PFZ in this grid - do not run code if there isn't one.
c       Also - need to avoid assumptions about grid geometry - NRS may not
c       be next to IRSEP for all grids - double nulls are another example. Code
c       needs to be rewritten to use IRINS/IROUTS exclusively.
c      
c
c     Small bug fix: IKTO is less than IKTI so the original test was always
c     true which meant that the following code was still executed for grids
c     without a PFZ. The desginations of IKTO and IKTI came from X-point up JET
c     grids where IK=1 corresponds to the outer target. 
c
      if (ikto.ge.1.and.ikti.le.nks(irsep)) then 
c
c      if (ikti.ge.1.and.ikto.le.nks(irsep)) then 
c

      DTHETG = THETAG(IKTI  - 1,IRSEP) -
     +         THETAG(IKTO  + 1,IRSEP) + 1.0

      IR  = NRS
      IRR = IRSEP
      IKR = IKTO

      WRITE(50,'(A,4I4,2F8.2)') 'PFZ  THETA reference:',
     +      IKR,NKS(IR),IR,IRR,ALPH(IKR,IR),ALPH(IKINS(IKR,IR),IRR)

      THETAG(IKR,NRS) = THETAG(IKOUTS(IKR,NRS),IRSEP)

      DO IK = IKR - 1, 1, -1
        THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKOUTS(IK,IR),IRR,1)
      ENDDO

      THETAG(IKR+1,NRS) = THETAG(IKOUTS(IKR+1,NRS),IRSEP) - DTHETG

      DO IK = IKR + 2, NKS(IR)
        THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKOUTS(IK,IR),IRR,1) - DTHETG
      ENDDO

      DO IR = NRS-1, IRWALL+1, -1
         IRR = IROUTS(1,IR)
         IKR = IKTO

         WRITE(50,'(A,4I4,2F8.2)') 'PFZ  THETA reference:',
     +     IKR,NKS(IR),IR,IRR,ALPH(IKR,IR),ALPH(IKINS(IKR,IR),IRR)

         DO IK = IKR-1, 1, -1
            THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKOUTS(IK,IR),IRR,1)
         ENDDO

         DO IK = IKR, NKS(IR)
            THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKOUTS(IK,IR),IRR,1)
         ENDDO
      ENDDO
c
c     End of test for existence of PFZ region on the grid
c
      endif


c
c ----------------------------------------------------------------------
c Core Region
c ----------------------------------------------------------------------
c

      DO IR = IRSEP-1, 1, -1

        FACTOR = 360.0

        IRR = IROUTS(INT(0.5*NKS(IR)),IR)

        DO IK = 2, NKS(IR) - 1

          F1 = ABS(ALPH(IK           ,IR ) - 90.0)
          F2 = ABS(ALPH(IKOUTS(IK,IR),IRR) - 90.0)

          IF ((F1+F2).LT.FACTOR) THEN
            IKR    = IK
            FACTOR = F1 + F2
          ENDIF

        ENDDO

        IF (FACTOR.GT.2.0)
     +    CALL WARN('Calctheta','Bad core THETA reference for ring',
     +              IR,0,0)

30      CONTINUE

        WRITE(50,'(A,4I4,2F8.2)') 'Core THETA reference:',
     +    IKR,NKS(IR),IR,IRR,ALPH(IKR,IR),ALPH(IKOUTS(IKR,IR),IRR)

        THETAG(IKR,IR) = THETAG(IKOUTS(IKR,IR),IRR)

        DO IK = IKR-1, 1, -1
          THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKOUTS(IK,IR),IRR,1)
        ENDDO

        DO IK = IKR+1, NKS(IR)-1
          THETAG(IK,IR) = CALTHETA_Old(IK,IR,IKOUTS(IK,IR),IRR,1)
        ENDDO
c
c Estimate THETA value for the degenerate cells at the end of the
c core rings - nothing fancy here:
c
        THETAG(NKS(IR),IR) = THETAG(NKS(IR)-1,IR) + 1
      ENDDO
c
c ----------------------------------------------------------------------
c Initialize THETAT for CTARGOPT 1&6 - other options are in WALLS.D6A:
c ----------------------------------------------------------------------
c
c slnote 06/07/95
c
c As far as I can tell, THETAT is only defined for CTARGOPT.EQ.1,
c but it is used in DIV depending on NORTHOPT, not CTARGOPT.  This
c was why I was running into some trouble with cross field transport near
c the target. Maybe a check should be added to identify illegal
c CTARGOPT and NORTHOPT combinations.
c
c Is THETAT defined in the core?
c No
c
      IF (CTARGOPT.EQ.1) THEN
        CALL WARN('Calctheta','THETAT questionable for target option 1',
     +            0,0,0)

         id = 0
         do ir = irwall, irsep, -1
            id = id + 1
            THETAT(ID)=THETAG(NKS(IR),IR)
         enddo
         do ir = nrs, irtrap, -1
            id = id + 1
            THETAT(ID)=THETAG(NKS(IR),IR)
         enddo
         do ir = irtrap, nrs
            id = id + 1
            THETAT(ID)=THETAG(1,IR)
         enddo
         do ir = irsep, irwall
            id = id + 1
            THETAT(ID)=THETAG(1,IR)
         enddo

      ELSEIF (CTARGOPT.EQ.6) THEN

        ID = 0

        DO IR = IRWALL, IRSEP, -1
          ID = ID + 1
          THETAT(ID) = EXTTHE_OLD(NKS(IR),IR,1)
        ENDDO

         do ir = nrs, irtrap, -1
            id = id + 1
            THETAT(ID)=EXTTHE_OLD(NKS(IR),IR,1)
         enddo

         do ir = irtrap, nrs
            id = id + 1
            THETAT(ID)=EXTTHE_OLD(1,IR,0)
         enddo

         do ir = irsep, irwall
            id = id + 1
            THETAT(ID)=EXTTHE_OLD(1,IR,0)
         enddo

      ENDIF
c
c ----------------------------------------------------------------------
c Output THETAG and THETAT:
c ----------------------------------------------------------------------
c
c The THETA values are shifted so that the minimum THETA value on the
c grid is 1.  This is important for some of the code in the DIV module
c which assumes that THETA less than 0 implies particle impact at the
c target.
c
      THEMIN =  1.0

      DO IR = 1, NRS
         DO IK = 1, NKS(IR)
            IF (THETAG(IK,IR).LT.THEMIN) THEMIN = THETAG(IK,IR)
         ENDDO
      ENDDO

      IF (CTARGOPT.EQ.1.OR.CTARGOPT.EQ.6) THEN
         DO I = 1, ID
            IF (THETAT(I).LT.THEMIN) THEMIN = THETAT(I)
          ENDDO
      ENDIF

      IF (THEMIN.LT.1.0) THEN
         THESFT = -THEMIN + 1.0
      ELSE
         THESFT = 0.0
      ENDIF

      DO IR = 1, NRS
         DO IK = 1, NKS(IR)
            THETAG(IK,IR) = THETAG(IK,IR) + THESFT
         ENDDO
      ENDDO

c
c     Write data to main debug file if print option is set
c
      if (cprint.eq.3.or.cprint.eq.9) then 
c
         WRITE (6,*) ' '
         WRITE (6,*) 'THETAG values generated by DIVIMP:'

         DO IR=1,NRS
            WRITE (6,*) ' '
            WRITE (6,'(2A6,4A15)') 'IK','IR','THETAG IN','THETAG',
     +                         'THETAG OUT','ANGLE'
            DO IK=1,NKS(IR)
               WRITE (6,'(2I6,4G15.7)') IK,IR,
     +                               THETAG(IKINS(IK,IR),IRINS(IK,IR)),
     +                               THETAG(IK,IR),
     +                               THETAG(IKOUTS(IK,IR),IROUTS(IK,IR))
     +                               ,ALPH(IK,IR)

            ENDDO
         ENDDO
c
      endif

      WRITE (50,*) ' '
      WRITE (50,*) 'THETAG values generated by DIVIMP:'

      DO IR=1,NRS
         WRITE (50,*) ' '
         WRITE (50,'(2A6,4A15)') 'IK','IR','THETAG IN','THETAG',
     +                         'THETAG OUT','ANGLE'
         DO IK=1,NKS(IR)
            WRITE (50,'(2I6,4G15.7)') IK,IR,
     +                               THETAG(IKINS(IK,IR),IRINS(IK,IR)),
     +                               THETAG(IK,IR),
     +                               THETAG(IKOUTS(IK,IR),IROUTS(IK,IR))
     +                               ,ALPH(IK,IR)

c
c Check for decreasing THETAG for increaing IK - core not checked:
c
            IF (ik.NE.1) THEN
              IF (THETAG(IK,IR).LT.THETAG(IK-1,IR).AND.
     +            IR.GE.IRSEP.AND.THETAG(IK,IR).NE.-1)
     +          CALL WARN('Calctheta','THETA sequence error',IK,IR,0)
            ENDIF
         ENDDO
      ENDDO

      IF (CTARGOPT.EQ.1.OR.CTARGOPT.EQ.6) THEN

        WRITE (50,*) ' '
        WRITE (50,*) 'THETAT information: '
        WRITE (50,*) ' '

        WRITE (50,'(A4,A15)') 'ID','THETAT'
        DO I = 1, ID
c
c THETAT may have values less than 1.0:
c
          THETAT(I) = THETAT(I) + THESFT

          WRITE (50,'(I4,F15.8)') I,THETAT(I)
        ENDDO

       ENDIF

      END
c
c = SUBROUTINE END: CALCTHETA ==========================================
c
c
c ======================================================================
c FUNCTION: CALTHETA_Old
c ======================================================================
c
      REAL FUNCTION CALTHETA_Old(IK,IR,IKRE,IRRE,ISIDE)
      implicit none
c     include "params"
      include 'params'
c     include "cgeom"
      include 'cgeom'
c     include "comtor"
      include 'comtor'

      INTEGER  IK,IR,IKRE,IRRE,IKR,IRR,IDR,ISIDE,IV1,IV2,IK1,IK2

      REAL     THECRE,THECCU,THEA,THEB,THESFT,THEX,THESCL,THETA
      REAL     DIST,DISTB,DISTF
c
      real     cenlen_Old
c
      IKR = IKRE
      IRR = IRRE
c
c Check that the reference cell is on the grid:
c
      IF (IKR.LT.1.OR.IKR.GT.NKS(IRR)) THEN
         WRITE(50,'(A,4I4)') 'Error (CalTheta): Bad reference cell:',
     +                        IK,IR,IKR,IRR
         WRITE( 0,'(A,4I4)') 'Error (CalTheta): Bad reference cell:',
     +                        IK,IR,IKR,IRR
         STOP
      ENDIF
c
c Trivial case - both cells are orthogonal:
c
      IF (ABS(COSALPH(IK  ,IR)).LT.1.8E-02.AND.
     +    ABS(COSALPH(IKR,IRR)).LT.1.8E-02) THEN
         CALTHETA_Old = THETAG(IKR,IRR)
         RETURN
      ENDIF
c
c slnote - Replace the '1', '2', etc. with parameters
c
      IF (ISIDE.EQ.1) THEN
         IV2 = 4
         IV1 = 1
      ELSE
         IV2 = 3
         IV1 = 2
      ENDIF

      IDR = KORPG(IKR,IRR)

10    THEA = RVERTP(IV2,IDR) - RVERTP(IV1,IDR)
      THEB = ZVERTP(IV2,IDR) - ZVERTP(IV1,IDR)
c
c For non-triangular cells:
c
        THECRE = ((RS(IKR,IRR) - RVERTP(IV1,IDR)) * THEA +
     +            (ZS(IKR,IRR) - ZVERTP(IV1,IDR)) * THEB) /
     +           (THEA**2 + THEB**2)

        THECCU = ((RS(IK,IR) - RVERTP(IV1,IDR)) * THEA +
     +            (ZS(IK,IR) - ZVERTP(IV1,IDR)) * THEB) /
     +           (THEA**2 + THEB**2)

c
c Check to see if IKR needs to be shifted (this is necessary if
c another cell on the reference ring would be a better reference
c cell than the current one):
c
        IF (((THECRE.LT.0.0).OR.(THECRE.GT.1.0)).OR.
     +      ((THECCU.LT.0.0).OR.(THECCU.GT.1.0))) THEN

          DIST  = CENLEN_OLD(IK,IR,IKR,IRR)
          DISTF = CENLEN_OLD(IK,IR,MIN(IKR+1,NKS(IRR)),IRR)
          DISTB = CENLEN_OLD(IK,IR,MAX(IKR-1,1)       ,IRR)

          IF (DISTB.LT.DIST.AND.DISTB.LT.DISTF) THEN

            IF (IK.GT.1) THEN
c              IF (IR.GE.27.AND.IK.GT.NKS(IR)-10) THEN
c                CALL WARN('CalTheta','Backward adjustment at',IK,IR,0)
c              ENDIF
              IKR = IKR-1
              GOTO 10
            ENDIF

          ELSEIF (DISTF.LT.DIST.AND.DISTF.LT.DISTB) THEN

            IF (IKR.LT.NKS(IRR)) THEN
c              IF (IR.GE.27.AND.IK.GT.NKS(IR)-10) THEN
c                CALL WARN('CalTheta','Forward adjustment at',IK,IR,0)
c              ENDIF
              IKR = IKR + 1
              GOTO 10
            ENDIF

          ENDIF

        ENDIF
c
c Everything is all set - calculate THETA:
c
        THESFT = SQRT(THEA**2 * (THECRE - THECCU)**2 +
     +                THEB**2 * (THECRE - THECCU)**2)

         IF (THECCU.LT.THECRE) THEN
           IF (IKR.GT.1) THEN
             IK1    = IKR - 1
           ELSE
             IK1    = IKR + 1
             THESFT = -THESFT
           ENDIF
         ELSE
           IF (IKR.LT.NKS(IRR)) THEN
             IK1    = IKR + 1
           ELSE
             IK1    = IKR - 1
             THESFT = -THESFT
           ENDIF
         ENDIF

         CALTHETA_Old = THETAG(IKR,IRR) + THESFT *
     +              (THETAG(IK1,IRR) - THETAG(IKR,IRR)) /
     +              CENLEN_OLD(IK1,IRR,IKR,IRR)

c         IF (IR.GE.27.AND.IK.GT.NKS(IR)-10) THEN
c           WRITE(50,'(A,6I6,2F7.2,3F15.7)')
c     +        'Theta: ',IK,IR,IKR,IK1,IRR,IKRE,
c     +        THECCU,THECRE,THEA,THEB,CALTHETA_Old
c         ENDIF

      RETURN
      END
c
c FUNCTION END: CALTHETA_Old ===============================================
c
c
c ======================================================================
c FUNCTION: CENLEN_OLD
c ======================================================================
c
      REAL FUNCTION CENLEN_OLD(IK1,IR1,IK2,IR2)
      implicit none
c     include "params"
      include 'params'
c     include "cgeom"
      include 'cgeom'

      INTEGER IK1,IK2,IR1,IR2
      REAL    RCEN,ZCEN,DELTAL
c
c If the cells are not adjacent and on the same ring then the distance
c returned is simply the R,Z displacement between cell centers:
c
      IF (ABS(IK1-IK2).GT.1.OR.IR1.NE.IR2) THEN
         CENLEN_OLD = SQRT((RS(IK1,IR1) - RS(IK2,IR2))**2 +
     +                 (ZS(IK1,IR1) - ZS(IK2,IR2))**2)
      ELSE
c
c slnote - Add an option reference
c
        IF (IK1.LT.IK2) THEN
          RCEN = 0.5 * (RVERTP(3,KORPG(IK1,IR1)) +
     +                  RVERTP(4,KORPG(IK1,IR1)))
          ZCEN = 0.5 * (ZVERTP(3,KORPG(IK1,IR1)) +
     +                  ZVERTP(4,KORPG(IK1,IR1)))
        ELSE
          RCEN = 0.5 * (RVERTP(3,KORPG(IK2,IR2)) +
     +                  RVERTP(4,KORPG(IK2,IR2)))
          ZCEN = 0.5 * (ZVERTP(3,KORPG(IK2,IR2)) +
     +                  ZVERTP(4,KORPG(IK2,IR2)))
        ENDIF

        CENLEN_OLD=SQRT((RS(IK1,IR1)-RCEN)**2 + (ZS(IK1,IR1)-ZCEN)**2) +
     +         SQRT((RS(IK2,IR2)-RCEN)**2 + (ZS(IK2,IR2)-ZCEN)**2)
      ENDIF

      RETURN
      END
c
c FUNCTION END: CENLEN_OLD ===============================================
c
c
c ======================================================================
c FUNCTION: EXTTHE_OLD(IK1,IR1,IDIR)
c ======================================================================
c
      REAL FUNCTION EXTTHE_OLD(IK,IR,IDIR)
      implicit none
      INTEGER IK,IR,IDIR

c     include "params"
      include 'params'
c     include "cgeom"
      include 'cgeom'
c     include "comtor"
      include 'comtor'

c      COMMON /SLCOM/ IKBRKS(4,MAXNRS),IKBRKE(4,MAXNRS),NBR,
c     +               NBREAK(MAXNRS),NBS,IRORG(MAXNRS)
c     +
c     +              ,NCELL,IRWAL(MAXNKS),IKWAL(MAXNKS)
c      INTEGER NCELL,IRWAL,IKWAL

c      INTEGER IKBRKS,IKBRKE,NBR,NBREAK
c      INTEGER NBS,IRORG

      REAL    RT,ZT,L1,L2
      REAL    CENLEN_OLD
      INTEGER IN

      IF (CTARGOPT.NE.6) THEN
         CALL WARN('Extthe_Old','Target option not set to 6 as '//
     .             'required',
     +             0,0,0)
         STOP
      ENDIF
c
      IN = KORPG(IK,IR)
c
c     IF there is NO polygon reference for this cell - then exit
c
      if (in.eq.0) then
c
         write (50,*) 'NO polygon information for cell:',ik,ir
c
c        Return base theta value
c
         extthe_Old = thetag(ik,ir)
c
         return
c
      endif

c
c Extapolate backward:
c
      IF (IDIR.EQ.0) THEN

         IF (THETAG(IK+1,IR).EQ.-1) THEN
            WRITE( 0,*) 'Error   (EXTTHE_OLD): Bad cell reference:',
     .                   IK+1,IR
            WRITE(50,*) 'Error   (EXTTHE_OLD): Bad cell reference:',
     .                   IK+1,IR
            STOP
         ENDIF

         IN = KORPG(IK,IR)

         RT = 0.5 * (RVERTP(1,IN) + RVERTP(2,IN))
         ZT = 0.5 * (ZVERTP(1,IN) + ZVERTP(2,IN))

         L1 = SQRT((RS(IK,IR) - RT)**2 + (ZS(IK,IR) - ZT)**2)
         L2 = CENLEN_OLD(IK,IR,IK+1,IR)

         EXTTHE_OLD = THETAG(IK,IR)-(THETAG(IK+1,IR)-THETAG(IK,IR))*
     .                L1/L2
c
c Extrapolate forward:
c
      ELSEIF (IDIR.EQ.1) THEN

         IF (THETAG(IK-1,IR).EQ.-1) THEN
            WRITE( 0,*) 'Error   (EXTTHE_OLD): Bad cell reference:',
     .                  IK-1,IR
            WRITE(50,*) 'Error   (EXTTHE_OLD): Bad cell reference:',
     .                  IK-1,IR
            STOP
         ENDIF

         IN = KORPG(IK,IR)

         RT = 0.5 * (RVERTP(3,IN) + RVERTP(4,IN))
         ZT = 0.5 * (ZVERTP(3,IN) + ZVERTP(4,IN))

         L1 = SQRT((RS(IK,IR)-RT)**2 + (ZS(IK,IR)-ZT)**2)
         L2 = CENLEN_OLD(IK,IR,IK-1,IR)

         EXTTHE_OLD = THETAG(IK,IR)+(THETAG(IK,IR)-THETAG(IK-1,IR))*
     .                L1/L2

      ELSE
         WRITE( 0,*) 'Error   (EXTTHE_OLD): Illegal IDIR value.'
         WRITE(50,*) 'Error   (EXTTHE_OLD): Illegal IDIR value.'
         STOP
      ENDIF

c      WRITE(50,'(A,2I4,4F10.5)') 'Diag: ',IK,IR,L1,L2,EXTTHE_OLD,
c     +                            THETAG(IK,IR)

      RETURN
      END
c
c FUNCTION END: EXTTHE_OLD =================================================
c
c
c ======================================================================
c
c SUBROUTINE: GETRZ
c
c Given IK,IR,S and CROSS, GETRZ returns an approximate R,Z co-ordinate
c for the impurity... currently the R,Z position is returned assuming
c that CROSS is zero.
c
c ======================================================================
c
      SUBROUTINE GETRZ_Old(IK,IR,S,CROSS,R,Z)
      implicit none
c     include "params"
      include 'params'
c     include "cgeom"
      include 'cgeom'
c     include "comtor"
      include 'comtor'
c
c     Local Variables:
c
      INTEGER IK,IR,ID
      REAL    S,CROSS,R,Z,HOLDR,HOLDZ
      REAL    THETA1,THETA2,DS,SR,SZ,DCROSS,CROSSR,CROSSZ,DVR,DVZ
      REAL    RCEN,ZCEN
      REAL    DELTAD,DELTAS
c
      real    atan2c,cenlen_Old
      external atan2c,cenlen_Old
c
      HOLDR=R
      HOLDZ=Z

      ID = KORPG(IK,IR)
c
c Temporary fix? - this was not working for the JET CMOD grid!
c
      IF (KAREAS(IK,IR).EQ.0.0.OR.
     +    IR.EQ.IRWALL.OR.IR.EQ.IRTRAP.OR.IR.EQ.1) THEN
         R = RS(IK,IR)
         Z = ZS(IK,IR)
         RETURN
      ENDIF
c
c slnote - Need an option setting here?
c
      IF (S.LT.KSS(IK,IR)) THEN

        IF (IK.EQ.1) THEN
c
c Should be able to use target S values here:
c
          IF (IR.GE.IRSEP) THEN
            R = RS(IK,IR)
            Z = ZS(IK,IR)
            RETURN
          ELSE
c
c slnote - Not sure that this is the best place for the correction
c   for the first knot of the core rings.
c
c Needs some fixing...
c
             CALL WARN('GetRZ','CORE! Last cell DS calculation',IK,IR,0)

             DELTAD = CENLEN_OLD(NKS(IR),IR,NKS(IR)-1,IR)
             DELTAS = KBACDS(NKS(IR),IR)
           ENDIF

         ELSE

           IF (PDOPT.EQ.1) THEN

             RCEN = 0.5 * (RVERTP(1,ID) + RVERTP(2,ID))
             ZCEN = 0.5 * (ZVERTP(1,ID) + ZVERTP(2,ID))

             DELTAD = SQRT((RCEN - RS(IK,IR))**2 +
     +                     (ZCEN - ZS(IK,IR))**2)
             DELTAS = KSS(IK,IR) - KSB(IK-1,IR)

           ELSE
             DELTAD = CENLEN_OLD(IK,IR,IK-1,IR)
             DELTAS = KBACDS(IK,IR)
           ENDIF

         ENDIF

         DS = (S - KSS(IK,IR)) * DELTAD / DELTAS

      ELSE

        IF (IK.EQ.NKS(IR)) THEN

          IF (IR.GE.IRSEP) THEN
            R = RS(IK,IR)
            Z = ZS(IK,IR)
            RETURN
          ELSE
c
c Needs some fixin...
c
            CALL WARN('GetRZ','CORE! Last cell DS calculation',IK,IR,0)

            DS = (S - KSS(IK,IR)) * CENLEN_OLD(1,IR,2,IR) / KFORDS(1,IR)
          ENDIF

        ELSE

          IF (PDOPT.EQ.1) THEN

            RCEN = 0.5 * (RVERTP(3,ID) + RVERTP(4,ID))
            ZCEN = 0.5 * (ZVERTP(3,ID) + ZVERTP(4,ID))

            DS = (S - KSS(IK,IR)) *
     +         SQRT((RCEN - RS(IK,IR))**2 + (ZCEN - ZS(IK,IR))**2) /
     +         (KSB(IK,IR) - KSS(IK,IR))

           ELSE

             DS = (S - KSS(IK,IR)) *
     +            CENLEN_OLD(IK,IR,IK+1,IR) / KFORDS(IK,IR)

           ENDIF

         ENDIF

      ENDIF

      DCROSS = 0.0

      DVZ = 0.5 * (ZVERTP(3,ID) + ZVERTP(4,ID))-
     +      0.5 * (ZVERTP(1,ID) + ZVERTP(2,ID))

      DVR = 0.5 * (RVERTP(3,ID) + RVERTP(4,ID))-
     +      0.5 * (RVERTP(1,ID) + RVERTP(2,ID))
c
      THETA1 = ATAN2C(ABS(DVZ),ABS(DVR))
      THETA2 = PI / 2.0 - THETA1

      SR = DS * COS(THETA1)
      SZ = DS * SIN(THETA1)

      CROSSR = DCROSS * COS(THETA2)
      CROSSZ = DCROSS * SIN(THETA2)

      IF     (DVR.GE.0.0.AND.DVZ.GE.0.0) THEN
         CROSSZ = -CROSSZ
      ELSEIF (DVR.GE.0.0.AND.DVZ.LT.0.0) THEN
         SZ     = -SZ
         CROSSR = -CROSSR
         CROSSZ = -CROSSZ
      ELSEIF (DVR.LT.0.0.AND.DVZ.GE.0.0) THEN
         SR     = -SR
      ELSEIF (DVR.LT.0.0.AND.DVZ.LT.0.0) THEN
         SR     = -SR
         SZ     = -SZ
         CROSSR = -CROSSR
      ENDIF

      R = RS(IK,IR) + SR + CROSSR
      Z = ZS(IK,IR) + SZ + CROSSZ

      RETURN
      END
c
c
c END SUBROUTINE: GETRZ ================================================
c

c
c ======================================================================
c SUBROUTINE: WARN
c ======================================================================
c
      SUBROUTINE WARN_OLD(MODULE,MESSAGE,PAR1,PAR2,PAR3)
      implicit none
      CHARACTER MODULE*(*),MESSAGE*(*)
      INTEGER   PAR1,PAR2,PAR3

      WRITE(50,'(4A,$)') ' Warning (',MODULE,'): ',MESSAGE
c      WRITE( 0,'(4A,$)') ' Warning (',MODULE,'): ',MESSAGE

      IF (PAR1.NE.0) THEN
        WRITE(50,'(I4,1X,$)') PAR1
c        WRITE( 0,'(I4,1X,$)') PAR1
      ENDIF

      IF (PAR2.NE.0) THEN
        WRITE(50,'(I4,1X,$)') PAR2
c        WRITE( 0,'(I4,1X,$)') PAR2
      ENDIF

      IF (PAR3.NE.0) THEN
        WRITE(50,'(I4,1X,$)') PAR3
c        WRITE( 0,'(I4,1X,$)') PAR3
      ENDIF

      WRITE(50,*) ' '
c      WRITE( 0,*) ' '

      RETURN
      END
c
c SUBROUTINE END: WARN ================================================
c


c
c ======================================================================
c SUBROUTINE: GETWID_OLD
c ======================================================================
c
      REAL FUNCTION GETWID_Old(IK,IR)
      implicit none
      include 'params'
      include 'cgeom'
      integer ik,ir
      REAL    WIDTOP,WIDBOT
      INTEGER ID

      IF (IR.EQ.1.OR.IR.EQ.IRWALL.OR.IR.EQ.IRTRAP) THEN

        GETWID_OLD = 0.0

      ELSE

        ID = KORPG(IK,IR)
        WIDBOT = SQRT((ZVERTP(1,ID) - ZVERTP(2,ID))**2.0 +
     +                (RVERTP(1,ID) - RVERTP(2,ID))**2.0)

        WIDTOP = SQRT((ZVERTP(3,ID) - ZVERTP(4,ID))**2.0 +
     +                (RVERTP(3,ID) - RVERTP(4,ID))**2.0)

        GETWID_OLD = 0.5 * (WIDTOP + WIDBOT)
      ENDIF

      RETURN
      END
c
c SUBROUTINE END: GETWID_OLD  ==============================================
c

